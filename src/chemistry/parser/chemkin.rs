/// Parser for CHEMKIN-II format:
///   - chem.inp: element/species/reaction declarations
///   - therm.dat: NASA polynomial thermodynamic data
///   - tran.dat:  transport parameters
///
/// CHEMKIN-II file formats:
///
/// **mech.inp** structure:
///   ELEMENTS
///     H O N AR
///   END
///   SPECIES
///     H2 H O O2 OH H2O HO2 H2O2 N2 AR
///   END
///   REACTIONS   [UNITS]
///     <species equation>  A  b  Ea
///     LOW/  A  b  Ea  /
///     TROE/ a  T3  T1  [T2] /
///     REV/  A  b  Ea  /
///     <species name>/  efficiency  /   (third-body efficiencies)
///     DUPLICATE
///   END
///
/// **therm.dat** — fixed-column, 4 lines per species:
///   Line 1: name(1-18), date(19-24), elements+counts(25-44), phase(45),
///            T_low(46-55), T_high(56-65), T_common(66-75), "1" at col 80
///   Line 2: a1..a5 (high range), each 15 chars wide, "2" at col 80
///   Line 3: a6 a7 (high range) + a1 a2 a3 (low range), "3" at col 80
///   Line 4: a4..a7 (low range), "4" at col 80
///
/// **tran.dat** — one line per species:
///   name  geometry  eps/k  sigma  dipole  alpha  Zrot
///   (geometry: 0=atom, 1=linear, 2=nonlinear)

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;
use crate::chemistry::mechanism::{Arrhenius, Mechanism, Reaction, RateType, ThirdBodySpec, TroeParams};
use crate::chemistry::species::{GeometryType, NasaPoly, Species, TransportParams};

// cal/mol → J/mol
const CAL_TO_J: f64 = 4.184;

/// Parse a CHEMKIN-II mechanism from text content.
///
/// # Arguments
/// * `mech_src`  — contents of the mechanism file (ELEMENTS/SPECIES/REACTIONS blocks)
/// * `therm_src` — contents of therm.dat; if None the THERMO block inside `mech_src` is used
/// * `tran_src`  — contents of tran.dat; if None transport data are left at defaults
pub fn parse_mechanism(
    mech_src: &str,
    therm_src: Option<&str>,
    tran_src: Option<&str>,
) -> Result<Mechanism> {
    let _elements = parse_elements(mech_src)?;
    let species_names = parse_species_list(mech_src)?;

    // Prefer separate therm.dat; fall back to embedded THERMO block in mech_src
    let thermo_text = therm_src.unwrap_or(mech_src);
    let mut species = parse_thermo(thermo_text, &species_names)?;

    if let Some(tran) = tran_src {
        attach_transport(&mut species, tran)?;
    }

    let reactions = parse_reactions(mech_src, &species)?;
    Ok(Mechanism { species, reactions })
}

/// Convenience wrapper: read files from disk and call `parse_mechanism`.
pub fn parse_file(
    mech_path: &str,
    therm_path: Option<&str>,
    tran_path: Option<&str>,
) -> Result<Mechanism> {
    let mech_src = std::fs::read_to_string(mech_path)
        .map_err(|e| anyhow!("Cannot read {mech_path}: {e}"))?;
    let therm_src = therm_path
        .map(|p| std::fs::read_to_string(p).map_err(|e| anyhow!("Cannot read {p}: {e}")))
        .transpose()?;
    let tran_src = tran_path
        .map(|p| std::fs::read_to_string(p).map_err(|e| anyhow!("Cannot read {p}: {e}")))
        .transpose()?;
    parse_mechanism(&mech_src, therm_src.as_deref(), tran_src.as_deref())
}

// ---------------------------------------------------------------------------
// Block extraction helpers
// ---------------------------------------------------------------------------

/// Strip inline comments (! to end of line) and return the cleaned string.
fn strip_comments(src: &str) -> String {
    src.lines()
        .map(|line| {
            if let Some(pos) = line.find('!') {
                &line[..pos]
            } else {
                line
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Extract the text content between `keyword … END` (case-insensitive).
/// Returns None if the keyword is not found.
fn extract_block<'a>(src: &'a str, keyword: &str) -> Option<&'a str> {
    let kw_upper = keyword.to_uppercase();
    // Find keyword at start of a word (may have trailing space or newline)
    let src_upper = src.to_uppercase();
    let kw_pos = src_upper.find(kw_upper.as_str())?;
    // Skip past the keyword line
    let after_kw = &src[kw_pos + keyword.len()..];
    // Find the next END token (whole word)
    let end_pos = {
        let upper_after = after_kw.to_uppercase();
        find_end_token(&upper_after)?
    };
    Some(after_kw[..end_pos].trim())
}

/// Find the position of a standalone `END` token (not part of a longer word).
fn find_end_token(upper: &str) -> Option<usize> {
    let mut pos = 0;
    while pos < upper.len() {
        if let Some(idx) = upper[pos..].find("END") {
            let abs = pos + idx;
            // Check that it's not part of a longer identifier
            let before_ok = abs == 0
                || !upper.as_bytes()[abs - 1].is_ascii_alphanumeric();
            let after_ok = abs + 3 >= upper.len()
                || !upper.as_bytes()[abs + 3].is_ascii_alphanumeric();
            if before_ok && after_ok {
                return Some(abs);
            }
            pos = abs + 3;
        } else {
            break;
        }
    }
    None
}

// ---------------------------------------------------------------------------
// ELEMENTS block
// ---------------------------------------------------------------------------

fn parse_elements(src: &str) -> Result<Vec<String>> {
    let clean = strip_comments(src);
    let block = match extract_block(&clean, "ELEMENTS") {
        Some(b) => b.to_string(),
        None => return Ok(vec![]),
    };
    let elements: Vec<String> = block
        .split_whitespace()
        .map(|s| s.to_uppercase())
        .collect();
    Ok(elements)
}

// ---------------------------------------------------------------------------
// SPECIES block
// ---------------------------------------------------------------------------

fn parse_species_list(src: &str) -> Result<Vec<String>> {
    let clean = strip_comments(src);
    let block = match extract_block(&clean, "SPECIES") {
        Some(b) => b.to_string(),
        None => bail!("No SPECIES block found in mechanism"),
    };
    let names: Vec<String> = block.split_whitespace().map(str::to_string).collect();
    if names.is_empty() {
        bail!("SPECIES block is empty");
    }
    Ok(names)
}

// ---------------------------------------------------------------------------
// THERMO / therm.dat — NASA 7-coefficient fixed-column format
// ---------------------------------------------------------------------------

/// Atom weights in kg/mol (upper-case keys).
fn atom_weights() -> HashMap<&'static str, f64> {
    [
        ("H",  1.008e-3),
        ("O",  15.999e-3),
        ("N",  14.007e-3),
        ("C",  12.011e-3),
        ("AR", 39.948e-3),
        ("HE", 4.003e-3),
        ("S",  32.06e-3),
        ("CL", 35.45e-3),
        ("F",  18.998e-3),
        ("BR", 79.904e-3),
    ]
    .into()
}

/// Parse a f64 from a fixed-width column, trimming whitespace. Returns 0.0 if blank.
fn parse_col(line: &str, start: usize, end: usize) -> f64 {
    let end = end.min(line.len());
    if start >= end {
        return 0.0;
    }
    line[start..end].trim().parse::<f64>().unwrap_or(0.0)
}

/// Parse NASA thermo data from therm.dat (or the THERMO block of mech.inp).
/// Only species listed in `species_names` are returned, in that order.
fn parse_thermo(src: &str, species_names: &[String]) -> Result<Vec<Species>> {
    let aw = atom_weights();

    // Collect all NASA entries keyed by upper-case name.
    let mut nasa_map: HashMap<String, Species> = HashMap::new();

    // If the source has a THERMO block, extract it; otherwise treat whole src as therm.dat.
    let clean = strip_comments(src);
    let body = extract_block(&clean, "THERMO").unwrap_or(clean.as_str());

    // Iterate over lines, looking for groups of 4 that form a NASA entry.
    let lines: Vec<&str> = body.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let l1 = lines[i];
        // The record indicator is in column 79 (0-based).  Line-1 has "1".
        if l1.len() < 79 {
            i += 1;
            continue;
        }
        let indicator = l1.chars().nth(79).unwrap_or(' ');
        if indicator != '1' {
            i += 1;
            continue;
        }
        // Need at least 3 more lines
        if i + 3 >= lines.len() {
            i += 1;
            continue;
        }
        let l2 = lines[i + 1];
        let l3 = lines[i + 2];
        let l4 = lines[i + 3];

        // Species name: columns 0-17 (trim whitespace)
        let name_raw = if l1.len() >= 18 { l1[..18].trim() } else { l1.trim() };
        let name = name_raw.to_string();

        // Temperature bounds: T_low(45-54), T_high(55-64), T_common(65-74)
        let t_low    = parse_col(l1, 45, 55);
        let t_high   = parse_col(l1, 55, 65);
        let t_common = parse_col(l1, 65, 75);

        // Parse element composition from columns 24-43: 4 pairs of (symbol, count)
        // Each pair occupies 5 chars: 2 for symbol, 3 for count
        let mut composition: HashMap<String, f64> = HashMap::new();
        for k in 0..4 {
            let col_start = 24 + k * 5;
            let col_end   = col_start + 5;
            if l1.len() < col_end {
                break;
            }
            let sym_raw = l1[col_start..col_start + 2].trim().to_uppercase();
            let cnt_raw = l1[col_start + 2..col_end].trim();
            if sym_raw.is_empty() {
                continue;
            }
            let cnt: f64 = cnt_raw.parse().unwrap_or(0.0);
            if cnt > 0.0 {
                composition.insert(sym_raw, cnt);
            }
        }

        // Molecular weight
        let mut mw = 0.0_f64;
        for (el, &cnt) in &composition {
            if let Some(&w) = aw.get(el.as_str()) {
                mw += cnt * w;
            }
        }

        // High-T coefficients (7): lines 2 and 3
        // Line 2: a1..a5 (cols 0..74, each 15 chars wide)
        // Line 3: a6, a7 (first 30 chars) then a1_low, a2_low, a3_low
        let high_coeffs = [
            parse_col(l2,  0, 15),
            parse_col(l2, 15, 30),
            parse_col(l2, 30, 45),
            parse_col(l2, 45, 60),
            parse_col(l2, 60, 75),
            parse_col(l3,  0, 15),
            parse_col(l3, 15, 30),
        ];

        // Low-T coefficients: last 3 on line 3, then 4 on line 4
        let low_coeffs = [
            parse_col(l3, 30, 45),
            parse_col(l3, 45, 60),
            parse_col(l3, 60, 75),
            parse_col(l4,  0, 15),
            parse_col(l4, 15, 30),
            parse_col(l4, 30, 45),
            parse_col(l4, 45, 60),
        ];

        let nasa_low = NasaPoly {
            t_low,
            t_high: t_common,
            coeffs: low_coeffs,
        };
        let nasa_high = NasaPoly {
            t_low: t_common,
            t_high,
            coeffs: high_coeffs,
        };

        // Default transport (will be overridden by tran.dat if provided)
        let transport = TransportParams {
            geometry: GeometryType::Nonlinear,
            well_depth: 0.0,
            diameter: 3.0,
            dipole_moment: 0.0,
            polarizability: 0.0,
            rot_relax: 1.0,
        };

        let sp = Species {
            name: name.clone(),
            molecular_weight: mw,
            composition,
            nasa_high,
            nasa_low,
            transport,
        };

        nasa_map.insert(name.to_uppercase(), sp);
        i += 4;
    }

    // Return species in the order declared in the SPECIES block.
    let mut out = Vec::new();
    for sp_name in species_names {
        let key = sp_name.to_uppercase();
        match nasa_map.remove(&key) {
            Some(sp) => out.push(sp),
            None => bail!("No thermo data found for species '{sp_name}'"),
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// tran.dat — transport parameters
// ---------------------------------------------------------------------------

fn attach_transport(species: &mut Vec<Species>, tran: &str) -> Result<()> {
    let clean = strip_comments(tran);

    // Build index map for O(1) lookup
    let idx_map: HashMap<String, usize> = species
        .iter()
        .enumerate()
        .map(|(i, sp)| (sp.name.to_uppercase(), i))
        .collect();

    for line in clean.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut tokens = line.split_whitespace();
        let name = match tokens.next() {
            Some(n) => n.to_uppercase(),
            None => continue,
        };
        if let Some(&idx) = idx_map.get(&name) {
            // geometry  eps/k  sigma  dipole  alpha  Zrot
            let geom_code: i32 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(2);
            let eps_k:  f64 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(0.0);
            let sigma:  f64 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(3.0);
            let dipole: f64 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(0.0);
            let alpha:  f64 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(0.0);
            let zrot:   f64 = tokens.next().and_then(|t| t.parse().ok()).unwrap_or(1.0);

            let geometry = match geom_code {
                0 => GeometryType::Atom,
                1 => GeometryType::Linear,
                _ => GeometryType::Nonlinear,
            };

            species[idx].transport = TransportParams {
                geometry,
                well_depth: eps_k,
                diameter: sigma,
                dipole_moment: dipole,
                polarizability: alpha,
                rot_relax: zrot,
            };
        }
        // Species not in list → silently skip
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// REACTIONS block
// ---------------------------------------------------------------------------

/// Tokenise the REACTIONS block into logical "reaction records".
/// A record is the stoichiometry+Arrhenius line plus any modifier lines
/// (LOW/, TROE/, REV/, efficiency entries, DUPLICATE).
struct ReactionRecord {
    /// The main reaction equation line, e.g. "H + O2 <=> OH + O  3.52e+15  -0.7  17069"
    equation_line: String,
    low:       Option<[f64; 3]>,
    troe:      Option<[f64; 4]>,      // [a, T3, T1, T2] (T2 may be 0 if absent)
    rev:       Option<[f64; 3]>,
    efficiencies: Vec<(String, f64)>, // (species_name, efficiency)
    duplicate: bool,
}

fn parse_reactions(src: &str, species: &[Species]) -> Result<Vec<Reaction>> {
    let clean = strip_comments(src);

    // Find the REACTIONS block.  It may have optional unit keywords after it.
    // We need the raw block body (between REACTIONS line and END).
    let block = match extract_block(&clean, "REACTIONS") {
        Some(b) => b.to_string(),
        None => return Ok(vec![]),
    };

    let sp_index = |name: &str| -> Option<usize> {
        species.iter().position(|s| s.name.eq_ignore_ascii_case(name))
    };

    let records = collect_reaction_records(&block)?;

    let mut out = Vec::new();
    for rec in records {
        let rxn = build_reaction(&rec, species, &sp_index)?;
        out.push(rxn);
    }
    Ok(out)
}

/// Walk through lines of the REACTIONS block and group them into records.
fn collect_reaction_records(block: &str) -> Result<Vec<ReactionRecord>> {
    let lines: Vec<&str> = block.lines().collect();
    let mut records: Vec<ReactionRecord> = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();
        i += 1;
        if line.is_empty() {
            continue;
        }

        // Skip unit declarations that may appear after REACTIONS keyword
        let upper = line.to_uppercase();
        if upper.starts_with("CAL") || upper.starts_with("KCAL")
            || upper.starts_with("JOUL") || upper.starts_with("KELVINS")
            || upper.starts_with("EVOLTS") || upper.starts_with("MOLES")
            || upper.starts_with("MOLECULES")
        {
            continue;
        }

        // Must contain either <=> or => to be an equation line
        if !line.contains("=>") {
            continue;
        }

        let mut rec = ReactionRecord {
            equation_line: line.to_string(),
            low: None,
            troe: None,
            rev: None,
            efficiencies: Vec::new(),
            duplicate: false,
        };

        // Consume modifier lines
        while i < lines.len() {
            let mline = lines[i].trim();
            let mupper = mline.to_uppercase();

            if mupper.starts_with("LOW") && mline.contains('/') {
                rec.low = Some(parse_slash_params(mline, 3)?);
                i += 1;
            } else if mupper.starts_with("TROE") && mline.contains('/') {
                let p = parse_slash_params_up_to_4(mline);
                rec.troe = Some(p);
                i += 1;
            } else if mupper.starts_with("REV") && mline.contains('/') {
                rec.rev = Some(parse_slash_params(mline, 3)?);
                i += 1;
            } else if mupper == "DUPLICATE" || mupper == "DUP" {
                rec.duplicate = true;
                i += 1;
            } else if mline.contains('/') && !mline.contains("=>") {
                // Third-body efficiency line: "H2O/ 6.0/ H2/ 2.0/ ..."
                let effs = parse_efficiency_line(mline);
                rec.efficiencies.extend(effs);
                i += 1;
            } else {
                break;
            }
        }

        records.push(rec);
    }
    Ok(records)
}

/// Parse `/  a  b  c  /` content from a keyword line.
fn parse_slash_params(line: &str, expected: usize) -> Result<[f64; 3]> {
    let inner = slash_inner(line);
    let vals: Vec<f64> = inner
        .split_whitespace()
        .filter_map(|t| t.parse().ok())
        .collect();
    if vals.len() < expected {
        bail!("Expected {expected} parameters in '{line}', got {}", vals.len());
    }
    Ok([vals[0], vals[1], vals[2]])
}

/// Parse TROE/ a T3 T1 [T2] / — up to 4 values, T2 defaults to 0.
fn parse_slash_params_up_to_4(line: &str) -> [f64; 4] {
    let inner = slash_inner(line);
    let vals: Vec<f64> = inner
        .split_whitespace()
        .filter_map(|t| t.parse().ok())
        .collect();
    [
        vals.first().copied().unwrap_or(0.0),
        vals.get(1).copied().unwrap_or(0.0),
        vals.get(2).copied().unwrap_or(0.0),
        vals.get(3).copied().unwrap_or(0.0),
    ]
}

/// Extract the text between the first and last `/` on a line.
fn slash_inner(line: &str) -> &str {
    let first = line.find('/').map(|p| p + 1).unwrap_or(0);
    let last  = line.rfind('/').unwrap_or(line.len());
    if first < last { &line[first..last] } else { "" }
}

/// Parse a third-body efficiency line: `H2/ 2.0/ H2O/ 6.0/ AR/ 0.7/`
fn parse_efficiency_line(line: &str) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    let parts: Vec<&str> = line.split('/').collect();
    // pairs: (name, value, name, value, …)
    let mut k = 0;
    while k + 1 < parts.len() {
        let name = parts[k].trim().to_string();
        let val: f64 = parts[k + 1].trim().parse().unwrap_or(1.0);
        if !name.is_empty() {
            out.push((name, val));
        }
        k += 2;
    }
    out
}

/// Build a `Reaction` from a parsed record.
fn build_reaction(
    rec: &ReactionRecord,
    _species: &[Species],
    sp_index: &impl Fn(&str) -> Option<usize>,
) -> Result<Reaction> {
    // Split equation and Arrhenius params.
    // The Arrhenius parameters are the last 3 whitespace-separated tokens on the equation line;
    // everything before them is the stoichiometric equation.
    let (eq_str, arrh) = split_equation_arrhenius(&rec.equation_line)?;

    let reversible = eq_str.contains("<=>");
    let separator = if reversible { "<=>" } else { "=>" };
    let parts: Vec<&str> = eq_str.splitn(2, separator).collect();
    if parts.len() != 2 {
        bail!("Cannot split equation: '{eq_str}'");
    }
    let reactants = parse_stoich(parts[0], sp_index)?;
    let products  = parse_stoich(parts[1], sp_index)?;

    // Determine if there is a third-body marker
    let has_m = eq_str.contains("(+M)") || eq_str.contains("(+m)")
        || {
            // bare "+ M" (not part of a species name)
            let eq_upper = eq_str.to_uppercase();
            eq_upper.contains("+ M ") || eq_upper.ends_with("+ M")
                || eq_upper.contains("+M ") || eq_upper.ends_with("+M")
        };

    let high = Arrhenius { a: arrh[0], b: arrh[1], ea: arrh[2] * CAL_TO_J };

    let rate = if let Some(low_params) = rec.low {
        let low = Arrhenius {
            a: low_params[0],
            b: low_params[1],
            ea: low_params[2] * CAL_TO_J,
        };
        let troe = rec.troe.map(|p| TroeParams {
            a:  p[0],
            t3: p[1],
            t1: p[2],
            t2: if p[3] != 0.0 { Some(p[3]) } else { None },
        });
        RateType::Falloff { high, low, troe }
    } else {
        RateType::Arrhenius { a: high.a, b: high.b, ea: high.ea }
    };

    // Third-body spec
    let third_body = if has_m || !rec.efficiencies.is_empty() {
        let efficiencies = rec
            .efficiencies
            .iter()
            .filter_map(|(name, eff)| {
                sp_index(name).map(|idx| (idx, *eff))
            })
            .collect();
        Some(ThirdBodySpec { efficiencies })
    } else {
        None
    };

    Ok(Reaction {
        reactants,
        products,
        rate,
        reversible,
        third_body,
        duplicate: rec.duplicate,
    })
}

/// Split the equation line into the stoichiometric string and [A, b, Ea].
/// The last 3 whitespace tokens are the Arrhenius params; the rest is the equation.
fn split_equation_arrhenius(line: &str) -> Result<(String, [f64; 3])> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 4 {
        bail!("Reaction line too short: '{line}'");
    }
    let n = tokens.len();
    let ea: f64 = tokens[n - 1].parse().map_err(|_| anyhow!("Cannot parse Ea in '{line}'"))?;
    let b:  f64 = tokens[n - 2].parse().map_err(|_| anyhow!("Cannot parse b in '{line}'"))?;
    let a:  f64 = tokens[n - 3].parse().map_err(|_| anyhow!("Cannot parse A in '{line}'"))?;
    // Re-join equation portion (everything before the last 3 tokens)
    let eq = tokens[..n - 3].join(" ");
    Ok((eq, [a, b, ea]))
}

/// Parse a stoichiometric expression (one side of the equation).
fn parse_stoich(
    side: &str,
    sp_index: &impl Fn(&str) -> Option<usize>,
) -> Result<Vec<(usize, f64)>> {
    let mut out = Vec::new();
    // Remove third-body markers
    let side = side
        .replace("(+M)", "")
        .replace("(+m)", "");
    for token in side.split('+') {
        let token = token.trim();
        if token.is_empty() || token.eq_ignore_ascii_case("M") {
            continue;
        }
        let (coeff, name) = parse_stoich_token(token);
        match sp_index(name) {
            Some(idx) => out.push((idx, coeff)),
            None => bail!("Unknown species in equation: '{name}'"),
        }
    }
    Ok(out)
}

fn parse_stoich_token(token: &str) -> (f64, &str) {
    let trimmed = token.trim();
    let end = trimmed.find(|c: char| c.is_alphabetic()).unwrap_or(0);
    if end == 0 {
        (1.0, trimmed)
    } else {
        let coeff_str = trimmed[..end].trim();
        let name = trimmed[end..].trim();
        let coeff = coeff_str.parse::<f64>().unwrap_or(1.0);
        (coeff, name)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn test_data_dir() -> String {
        let manifest = env!("CARGO_MANIFEST_DIR");
        format!("{manifest}/data/chemkin_test")
    }

    // ------------------------------------------------------------------
    // Helper: build a tiny in-memory H2/O2 mechanism for unit tests.
    // ------------------------------------------------------------------

    const MECH_SRC: &str = r#"
ELEMENTS
  H O N AR
END

SPECIES
  H2 O2 OH H2O H O HO2 AR N2
END

REACTIONS  CALORIES/MOLE
! H + O2 chain-branching step
  H + O2 <=> O + OH     3.547E+15  -0.406  1.6599E+04
! O + H2
  O + H2 <=> H + OH     5.080E+04   2.670  6.292E+03
! OH + H2
  OH + H2 <=> H + H2O   2.160E+08   1.510  3.430E+03
! H + OH + M falloff (illustrative)
  H + O2 (+M) <=> HO2 (+M)   1.475E+12   0.60   0.00
    LOW/ 3.482E+16  -0.411  -1.115E+03 /
    TROE/ 0.5  1.0E-30  1.0E+30 /
    H2O/ 10.6/ H2/ 1.5/ AR/ 0.67/
  2H + M <=> H2 + M     7.000E+17  -1.000  0.000
    H2/ 2.5/ H2O/ 12.0/ AR/ 0.0/
END
"#;

    const THERM_SRC: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/chemkin_test/therm.dat"
    ));

    const TRAN_SRC: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/chemkin_test/tran.dat"
    ));

    // ------------------------------------------------------------------
    // Block parser tests
    // ------------------------------------------------------------------

    #[test]
    fn test_parse_elements() {
        let els = parse_elements(MECH_SRC).unwrap();
        assert!(els.contains(&"H".to_string()));
        assert!(els.contains(&"O".to_string()));
        assert!(els.contains(&"AR".to_string()));
    }

    #[test]
    fn test_parse_species_list() {
        let names = parse_species_list(MECH_SRC).unwrap();
        assert_eq!(names[0], "H2");
        assert_eq!(names[1], "O2");
        assert!(names.contains(&"AR".to_string()));
        assert_eq!(names.len(), 9);
    }

    // ------------------------------------------------------------------
    // Thermo tests
    // ------------------------------------------------------------------

    #[test]
    fn test_thermo_h2_coefficients() {
        let names: Vec<String> = vec!["H2".to_string()];
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let sp = &sp_list[0];
        assert_eq!(sp.name, "H2");
        // Low-T first coeff: 2.34433112
        assert!((sp.nasa_low.coeffs[0] - 2.34433112).abs() < 1e-6,
            "low a1 = {}", sp.nasa_low.coeffs[0]);
        // High-T first coeff: 3.33727920
        assert!((sp.nasa_high.coeffs[0] - 3.33727920).abs() < 1e-6,
            "high a1 = {}", sp.nasa_high.coeffs[0]);
        // T ranges
        assert!((sp.nasa_low.t_low  - 200.0).abs() < 0.1);
        assert!((sp.nasa_low.t_high - 1000.0).abs() < 0.1);
        assert!((sp.nasa_high.t_high - 3500.0).abs() < 0.1);
    }

    #[test]
    fn test_thermo_molecular_weights() {
        let names: Vec<String> = ["H2", "O2", "H2O", "OH", "H", "O", "AR"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let mw = |n: &str| sp_list.iter().find(|s| s.name == n).unwrap().molecular_weight;

        assert!((mw("H2")  - 2.016e-3).abs()  < 1e-6, "H2 mw");
        assert!((mw("O2")  - 31.998e-3).abs() < 1e-6, "O2 mw");
        assert!((mw("H2O") - 18.015e-3).abs() < 1e-6, "H2O mw");
        assert!((mw("OH")  - 17.007e-3).abs() < 1e-6, "OH mw");
        assert!((mw("H")   - 1.008e-3).abs()  < 1e-6, "H mw");
        assert!((mw("O")   - 15.999e-3).abs() < 1e-6, "O mw");
        assert!((mw("AR")  - 39.948e-3).abs() < 1e-6, "AR mw");
    }

    #[test]
    fn test_thermo_missing_species_error() {
        let names: Vec<String> = vec!["NONEXISTENT".to_string()];
        assert!(parse_thermo(THERM_SRC, &names).is_err());
    }

    // ------------------------------------------------------------------
    // Transport tests
    // ------------------------------------------------------------------

    #[test]
    fn test_transport_geometry() {
        let names: Vec<String> = ["H2", "O2", "H2O", "H"]
            .iter().map(|s| s.to_string()).collect();
        let mut sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        attach_transport(&mut sp_list, TRAN_SRC).unwrap();

        let sp = |n: &str| sp_list.iter().find(|s| s.name == n).unwrap();
        assert_eq!(sp("H").transport.geometry,   GeometryType::Atom,     "H geometry");
        assert_eq!(sp("H2").transport.geometry,  GeometryType::Linear,   "H2 geometry");
        assert_eq!(sp("H2O").transport.geometry, GeometryType::Nonlinear,"H2O geometry");
    }

    #[test]
    fn test_transport_values() {
        let names: Vec<String> = vec!["N2".to_string()];
        let mut sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        attach_transport(&mut sp_list, TRAN_SRC).unwrap();
        let n2 = &sp_list[0].transport;
        assert!((n2.well_depth  - 97.53).abs()  < 0.01, "N2 eps/k");
        assert!((n2.diameter    - 3.621).abs()  < 0.001,"N2 sigma");
        assert!((n2.rot_relax   - 4.0).abs()    < 0.01, "N2 Zrot");
    }

    // ------------------------------------------------------------------
    // Reaction parsing tests
    // ------------------------------------------------------------------

    #[test]
    fn test_reaction_count() {
        let names: Vec<String> = ["H2","O2","OH","H2O","H","O","HO2","AR","N2"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let rxns = parse_reactions(MECH_SRC, &sp_list).unwrap();
        assert_eq!(rxns.len(), 5, "Expected 5 reactions, got {}", rxns.len());
    }

    #[test]
    fn test_arrhenius_reaction() {
        let names: Vec<String> = ["H2","O2","OH","H2O","H","O","HO2","AR","N2"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let rxns = parse_reactions(MECH_SRC, &sp_list).unwrap();
        // First reaction: H + O2 <=> O + OH, A=3.547e15, b=-0.406, Ea=16599 cal
        let rxn = &rxns[0];
        assert!(rxn.reversible);
        match rxn.rate {
            RateType::Arrhenius { a, b, ea } => {
                assert!((a - 3.547e15).abs() / 3.547e15 < 1e-4, "A = {a}");
                assert!((b - (-0.406)).abs() < 1e-4, "b = {b}");
                let ea_expected = 1.6599e4 * CAL_TO_J;
                assert!((ea - ea_expected).abs() / ea_expected < 1e-4, "Ea = {ea}");
            }
            _ => panic!("Expected Arrhenius for reaction 1"),
        }
    }

    #[test]
    fn test_falloff_reaction() {
        let names: Vec<String> = ["H2","O2","OH","H2O","H","O","HO2","AR","N2"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let rxns = parse_reactions(MECH_SRC, &sp_list).unwrap();
        // Reaction 4: H + O2 (+M) <=> HO2 (+M) with LOW and TROE
        let rxn = &rxns[3];
        assert!(rxn.third_body.is_some(), "Should have third body");
        match &rxn.rate {
            RateType::Falloff { high, low, troe } => {
                assert!((high.a - 1.475e12).abs() / 1.475e12 < 1e-4, "high.A");
                assert!((low.a  - 3.482e16).abs() / 3.482e16 < 1e-4, "low.A");
                let tr = troe.as_ref().expect("TROE params");
                assert!((tr.a - 0.5).abs() < 1e-6, "troe.a");
            }
            _ => panic!("Expected Falloff"),
        }
    }

    #[test]
    fn test_third_body_efficiencies() {
        let names: Vec<String> = ["H2","O2","OH","H2O","H","O","HO2","AR","N2"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let rxns = parse_reactions(MECH_SRC, &sp_list).unwrap();
        // Reaction 4: H + O2 (+M) — efficiency for H2O should be 10.6
        let tb = rxns[3].third_body.as_ref().expect("third body");
        let h2o_idx = sp_list.iter().position(|s| s.name == "H2O").unwrap();
        let eff = tb.efficiencies.iter().find(|&&(i, _)| i == h2o_idx);
        assert!(eff.is_some(), "H2O efficiency should be present");
        assert!((eff.unwrap().1 - 10.6).abs() < 0.01, "H2O eff = {}", eff.unwrap().1);
    }

    #[test]
    fn test_three_body_reaction() {
        let names: Vec<String> = ["H2","O2","OH","H2O","H","O","HO2","AR","N2"]
            .iter().map(|s| s.to_string()).collect();
        let sp_list = parse_thermo(THERM_SRC, &names).unwrap();
        let rxns = parse_reactions(MECH_SRC, &sp_list).unwrap();
        // Reaction 5: 2H + M <=> H2 + M — bare third body, Arrhenius
        let rxn = &rxns[4];
        assert!(rxn.third_body.is_some(), "Should have third body");
        match rxn.rate {
            RateType::Arrhenius { a, .. } => {
                assert!((a - 7.0e17).abs() / 7.0e17 < 1e-4, "A = {a}");
            }
            _ => panic!("Expected Arrhenius for three-body reaction"),
        }
    }

    // ------------------------------------------------------------------
    // Full round-trip from test data files
    // ------------------------------------------------------------------

    #[test]
    fn test_parse_file_roundtrip() {
        let dir = test_data_dir();
        let mech = parse_file(
            &format!("{dir}/mech.inp"),
            Some(&format!("{dir}/therm.dat")),
            Some(&format!("{dir}/tran.dat")),
        )
        .expect("parse test mechanism files");
        assert_eq!(mech.n_species(), 9, "species count");
        assert_eq!(mech.n_reactions(), 5, "reaction count");
        // Spot check: H2O molecular weight
        let h2o = mech.species_index("H2O").expect("H2O");
        let mw = mech.species[h2o].molecular_weight;
        assert!((mw - 18.015e-3).abs() < 1e-5, "H2O mw = {mw}");
        // Transport: H2O is nonlinear
        assert_eq!(
            mech.species[h2o].transport.geometry,
            GeometryType::Nonlinear,
            "H2O geometry"
        );
    }
}
