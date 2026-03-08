/// Parser for Cantera YAML mechanism format (.yaml).
///
/// Cantera YAML structure (relevant sections):
///   phases:
///     - name: gas
///       species: [...]
///       kinetics: gas
///       transport: mixture-averaged
///   species:
///     - name: H2
///       composition: {H: 2}
///       thermo:
///         model: NASA7
///         temperature-ranges: [200, 1000, 3500]
///         data:
///           - [2.34, ...]   # low range
///           - [3.34, ...]   # high range
///       transport:
///         model: gas
///         geometry: linear
///         well-depth: 38.0
///         diameter: 2.92
///         ...
///   reactions:
///     - equation: H + O2 <=> O + OH
///       rate-constant: {A: 2.644e+16, b: -0.6707, Ea: 17041 cal/mol}

use anyhow::Result;
use crate::chemistry::mechanism::Mechanism;

/// Parse a Cantera YAML mechanism from a file path.
pub fn parse_file(path: &str) -> Result<Mechanism> {
    let content = std::fs::read_to_string(path)?;
    parse_str(&content)
}

/// Parse a Cantera YAML mechanism from a string.
pub fn parse_str(content: &str) -> Result<Mechanism> {
    let doc: serde_yaml::Value = serde_yaml::from_str(content)?;
    let species = parse_species(&doc)?;
    let reactions = parse_reactions(&doc, &species)?;
    Ok(Mechanism { species, reactions })
}

fn parse_species(doc: &serde_yaml::Value) -> Result<Vec<crate::chemistry::species::Species>> {
    use crate::chemistry::species::{GeometryType, NasaPoly, Species, TransportParams};

    let species_list = doc["species"]
        .as_sequence()
        .ok_or_else(|| anyhow::anyhow!("Missing 'species' section in YAML"))?;

    let mut out = Vec::new();
    for sp_node in species_list {
        let name = sp_node["name"].as_str().unwrap_or("").to_string();

        // Molecular weight from composition
        let mw = molecular_weight_from_composition(&sp_node["composition"])?;

        // NASA7 coefficients
        let thermo = &sp_node["thermo"];
        let t_ranges = thermo["temperature-ranges"]
            .as_sequence()
            .ok_or_else(|| anyhow::anyhow!("Missing temperature-ranges for {name}"))?;
        let t_low = t_ranges[0].as_f64().unwrap_or(200.0);
        let t_mid = t_ranges[1].as_f64().unwrap_or(1000.0);
        let t_high = t_ranges[2].as_f64().unwrap_or(3500.0);

        let data = thermo["data"]
            .as_sequence()
            .ok_or_else(|| anyhow::anyhow!("Missing thermo data for {name}"))?;
        let low_coeffs = parse_7_coeffs(&data[0])?;
        let high_coeffs = parse_7_coeffs(&data[1])?;

        // Transport
        let tr = &sp_node["transport"];
        let geometry = match tr["geometry"].as_str().unwrap_or("nonlinear") {
            "atom"     => GeometryType::Atom,
            "linear"   => GeometryType::Linear,
            _          => GeometryType::Nonlinear,
        };
        let transport = TransportParams {
            geometry,
            well_depth:      tr["well-depth"].as_f64().unwrap_or(0.0),
            diameter:        tr["diameter"].as_f64().unwrap_or(3.0),
            dipole_moment:   tr["dipole"].as_f64().unwrap_or(0.0),
            polarizability:  tr["polarizability"].as_f64().unwrap_or(0.0),
            rot_relax:       tr["rotational-relaxation"].as_f64().unwrap_or(1.0),
        };

        out.push(Species {
            name,
            molecular_weight: mw,
            nasa_low:  NasaPoly { t_low, t_high: t_mid, coeffs: low_coeffs },
            nasa_high: NasaPoly { t_low: t_mid, t_high, coeffs: high_coeffs },
            transport,
        });
    }
    Ok(out)
}

fn parse_7_coeffs(node: &serde_yaml::Value) -> Result<[f64; 7]> {
    let seq = node.as_sequence()
        .ok_or_else(|| anyhow::anyhow!("NASA coeff block is not a sequence"))?;
    if seq.len() < 7 {
        anyhow::bail!("Expected 7 NASA coefficients, got {}", seq.len());
    }
    Ok([
        seq[0].as_f64().unwrap_or(0.0),
        seq[1].as_f64().unwrap_or(0.0),
        seq[2].as_f64().unwrap_or(0.0),
        seq[3].as_f64().unwrap_or(0.0),
        seq[4].as_f64().unwrap_or(0.0),
        seq[5].as_f64().unwrap_or(0.0),
        seq[6].as_f64().unwrap_or(0.0),
    ])
}

/// Compute molecular weight [kg/mol] from Cantera composition map.
/// Element atomic weights (most common elements only).
fn molecular_weight_from_composition(comp: &serde_yaml::Value) -> Result<f64> {
    let atom_weights: std::collections::HashMap<&str, f64> = [
        ("H", 1.008e-3), ("O", 15.999e-3), ("N", 14.007e-3),
        ("C", 12.011e-3), ("AR", 39.948e-3), ("Ar", 39.948e-3),
        ("HE", 4.003e-3), ("He", 4.003e-3),
        ("S", 32.06e-3),  ("CL", 35.45e-3), ("Cl", 35.45e-3),
    ].into();

    let map = comp.as_mapping()
        .ok_or_else(|| anyhow::anyhow!("Composition is not a mapping"))?;
    let mut mw = 0.0;
    for (el, count) in map {
        let el_str = el.as_str().unwrap_or("");
        let w = atom_weights.get(el_str)
            .ok_or_else(|| anyhow::anyhow!("Unknown element: {el_str}"))?;
        let n = count.as_f64().unwrap_or(0.0);
        mw += n * w;
    }
    Ok(mw)
}

fn parse_reactions(
    doc: &serde_yaml::Value,
    species: &[crate::chemistry::species::Species],
) -> Result<Vec<crate::chemistry::mechanism::Reaction>> {
    use crate::chemistry::mechanism::*;

    let rxn_list = match doc["reactions"].as_sequence() {
        Some(s) => s,
        None => return Ok(vec![]),
    };

    let sp_index = |name: &str| -> Option<usize> {
        species.iter().position(|s| s.name == name)
    };

    let mut out = Vec::new();
    for rxn_node in rxn_list {
        let eq = rxn_node["equation"].as_str().unwrap_or("");
        let reversible = eq.contains("<=>");
        let (lhs, rhs) = if reversible {
            let p: Vec<_> = eq.splitn(2, "<=>").collect();
            (p[0].trim(), p[1].trim())
        } else {
            let p: Vec<_> = eq.splitn(2, "=>").collect();
            (p[0].trim(), p[1].trim())
        };

        let reactants = parse_stoich(lhs, &sp_index)?;
        let products  = parse_stoich(rhs, &sp_index)?;

        // Rate constant
        let rate_node = &rxn_node["rate-constant"];
        let rate = if !rate_node.is_null() {
            let a  = rate_node["A"].as_f64().unwrap_or(0.0);
            let b  = rate_node["b"].as_f64().unwrap_or(0.0);
            let ea = parse_energy(rate_node)?;
            RateType::Arrhenius { a, b, ea }
        } else {
            // Falloff / pressure-dependent
            let high_node = &rxn_node["high-P-rate-constant"];
            let low_node  = &rxn_node["low-P-rate-constant"];
            let high = Arrhenius {
                a:  high_node["A"].as_f64().unwrap_or(0.0),
                b:  high_node["b"].as_f64().unwrap_or(0.0),
                ea: parse_energy(high_node)?,
            };
            let low = Arrhenius {
                a:  low_node["A"].as_f64().unwrap_or(0.0),
                b:  low_node["b"].as_f64().unwrap_or(0.0),
                ea: parse_energy(low_node)?,
            };
            let troe = if !rxn_node["Troe"].is_null() {
                let tr = &rxn_node["Troe"];
                Some(TroeParams {
                    a:  tr["A"].as_f64().unwrap_or(0.0),
                    t3: tr["T3"].as_f64().unwrap_or(1e30),
                    t1: tr["T1"].as_f64().unwrap_or(1e30),
                    t2: tr["T2"].as_f64(),
                })
            } else {
                None
            };
            RateType::Falloff { high, low, troe }
        };

        // Third body
        let third_body = if eq.contains("(+M)") || eq.contains("+ M") {
            let effs = if let Some(m) = rxn_node["efficiencies"].as_mapping() {
                m.iter().filter_map(|(k, v)| {
                    let name = k.as_str()?;
                    let eff = v.as_f64()?;
                    let idx = sp_index(name)?;
                    Some((idx, eff))
                }).collect()
            } else {
                vec![]
            };
            Some(ThirdBodySpec { efficiencies: effs })
        } else {
            None
        };

        out.push(Reaction { reactants, products, rate, reversible, third_body });
    }
    Ok(out)
}

fn parse_stoich(
    side: &str,
    sp_index: &impl Fn(&str) -> Option<usize>,
) -> Result<Vec<(usize, f64)>> {
    let mut out = Vec::new();
    // Remove third-body markers like "(+M)" or "+ M"
    let side = side.replace("(+M)", "").replace("(+m)", "");
    for token in side.split('+') {
        let token = token.trim();
        if token.is_empty() || token == "M" || token == "m" {
            continue;
        }
        // Extract optional leading coefficient, e.g. "2 OH" or "2OH"
        let (coeff, name) = parse_stoich_token(token);
        if let Some(idx) = sp_index(name) {
            out.push((idx, coeff));
        } else {
            anyhow::bail!("Unknown species in equation: '{name}'");
        }
    }
    Ok(out)
}

fn parse_stoich_token(token: &str) -> (f64, &str) {
    // Try to split leading float from species name
    let trimmed = token.trim();
    let end = trimmed
        .find(|c: char| c.is_alphabetic())
        .unwrap_or(0);
    if end == 0 {
        (1.0, trimmed)
    } else {
        let coeff_str = trimmed[..end].trim();
        let name = trimmed[end..].trim();
        let coeff = coeff_str.parse::<f64>().unwrap_or(1.0);
        (coeff, name)
    }
}

/// Parse activation energy from a rate node; convert cal/mol → J/mol if needed.
fn parse_energy(node: &serde_yaml::Value) -> Result<f64> {
    let ea_raw = node["Ea"].as_f64().unwrap_or(0.0);
    // Cantera YAML uses J/kmol by default in newer versions, but many files use cal/mol.
    // Heuristic: if Ea > 1e6, assume J/kmol; if < 1e5, assume cal/mol.
    // Proper handling requires reading the units field.
    let units = node["Ea-units"].as_str().unwrap_or("cal/mol");
    let ea_j_mol = match units {
        "cal/mol" | "cal" => ea_raw * 4.184,      // cal/mol → J/mol
        "kcal/mol"        => ea_raw * 4184.0,
        "J/mol"           => ea_raw,
        "kJ/mol"          => ea_raw * 1000.0,
        "K"               => ea_raw * 8.314462618, // K * R
        _                 => ea_raw * 4.184,       // assume cal/mol
    };
    Ok(ea_j_mol)
}
