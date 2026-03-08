/// Parser for Cantera YAML mechanism format (.yaml).
///
/// Cantera YAML structure (relevant sections):
///   units: {length: cm, time: s, quantity: mol, activation-energy: cal/mol}
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
    let ea_unit = parse_ea_unit(&doc);
    let species = parse_species(&doc)?;
    let reactions = parse_reactions(&doc, &species, &ea_unit)?;
    Ok(Mechanism { species, reactions })
}

/// Extract the activation-energy unit from the top-level `units:` block.
fn parse_ea_unit(doc: &serde_yaml::Value) -> String {
    doc["units"]["activation-energy"]
        .as_str()
        .unwrap_or("cal/mol")
        .to_string()
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
        let t_low  = t_ranges[0].as_f64().unwrap_or(200.0);
        let t_mid  = t_ranges[1].as_f64().unwrap_or(1000.0);
        let t_high = t_ranges[2].as_f64().unwrap_or(3500.0);

        let data = thermo["data"]
            .as_sequence()
            .ok_or_else(|| anyhow::anyhow!("Missing thermo data for {name}"))?;
        let low_coeffs  = parse_7_coeffs(&data[0])?;
        let high_coeffs = parse_7_coeffs(&data[1])?;

        // Transport
        let tr = &sp_node["transport"];
        let geometry = match tr["geometry"].as_str().unwrap_or("nonlinear") {
            "atom"   => GeometryType::Atom,
            "linear" => GeometryType::Linear,
            _        => GeometryType::Nonlinear,
        };
        let transport = TransportParams {
            geometry,
            well_depth:     tr["well-depth"].as_f64().unwrap_or(0.0),
            diameter:       tr["diameter"].as_f64().unwrap_or(3.0),
            dipole_moment:  tr["dipole"].as_f64().unwrap_or(0.0),
            polarizability: tr["polarizability"].as_f64().unwrap_or(0.0),
            rot_relax:      tr["rotational-relaxation"].as_f64().unwrap_or(1.0),
        };

        out.push(Species {
            name,
            molecular_weight: mw,
            nasa_low:  NasaPoly { t_low, t_high: t_mid,  coeffs: low_coeffs },
            nasa_high: NasaPoly { t_low: t_mid, t_high,  coeffs: high_coeffs },
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

/// Compute molecular weight [kg/mol] from a Cantera composition map.
/// Element symbol lookup is case-insensitive (normalised to Title Case).
fn molecular_weight_from_composition(comp: &serde_yaml::Value) -> Result<f64> {
    // Canonical map: Title-case symbol → atomic weight [kg/mol]
    let atom_weights: std::collections::HashMap<&str, f64> = [
        ("H",  1.008e-3),
        ("O",  15.999e-3),
        ("N",  14.007e-3),
        ("C",  12.011e-3),
        ("Ar", 39.948e-3),
        ("He", 4.003e-3),
        ("S",  32.06e-3),
        ("Cl", 35.45e-3),
        ("F",  18.998e-3),
        ("Br", 79.904e-3),
        ("I",  126.904e-3),
        ("Si", 28.085e-3),
        ("P",  30.974e-3),
    ].into();

    let map = comp.as_mapping()
        .ok_or_else(|| anyhow::anyhow!("Composition is not a mapping"))?;
    let mut mw = 0.0;
    for (el, count) in map {
        let el_str = el.as_str().unwrap_or("");
        // Normalise to Title Case: first char upper, rest lower.
        let canonical = normalise_element(el_str);
        let w = atom_weights.get(canonical.as_str())
            .ok_or_else(|| anyhow::anyhow!("Unknown element: {el_str}"))?;
        let n = count.as_f64().unwrap_or(0.0);
        mw += n * w;
    }
    Ok(mw)
}

/// Normalise an element symbol to Title Case (e.g. "AR" → "Ar", "ar" → "Ar", "H" → "H").
fn normalise_element(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => {
            let upper: String = first.to_uppercase().collect();
            let rest: String = chars.flat_map(|c| c.to_lowercase()).collect();
            upper + &rest
        }
    }
}

fn parse_reactions(
    doc: &serde_yaml::Value,
    species: &[crate::chemistry::species::Species],
    ea_unit: &str,
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

        let duplicate = rxn_node["duplicate"].as_bool().unwrap_or(false);

        // Determine rate type: Arrhenius, Falloff, or PLOG
        let rate_node  = &rxn_node["rate-constant"];
        let plog_node  = &rxn_node["pressure-dependent-Arrhenius"];
        let high_node  = &rxn_node["high-P-rate-constant"];

        let rate = if !rate_node.is_null() {
            let a  = rate_node["A"].as_f64().unwrap_or(0.0);
            let b  = rate_node["b"].as_f64().unwrap_or(0.0);
            let ea = parse_energy(rate_node, ea_unit)?;
            RateType::Arrhenius { a, b, ea }
        } else if !plog_node.is_null() {
            let entries = plog_node.as_sequence()
                .ok_or_else(|| anyhow::anyhow!("PLOG: pressure-dependent-Arrhenius is not a sequence in '{eq}'"))?;
            let mut rates = Vec::new();
            for entry in entries {
                let p_pa = parse_pressure(&entry["P"])?;
                let a  = entry["A"].as_f64().unwrap_or(0.0);
                let b  = entry["b"].as_f64().unwrap_or(0.0);
                let ea = parse_energy(entry, ea_unit)?;
                rates.push((p_pa, Arrhenius { a, b, ea }));
            }
            RateType::Plog { rates }
        } else {
            // Falloff (high-P / low-P)
            let low_node = &rxn_node["low-P-rate-constant"];
            let high = Arrhenius {
                a:  high_node["A"].as_f64().unwrap_or(0.0),
                b:  high_node["b"].as_f64().unwrap_or(0.0),
                ea: parse_energy(high_node, ea_unit)?,
            };
            let low = Arrhenius {
                a:  low_node["A"].as_f64().unwrap_or(0.0),
                b:  low_node["b"].as_f64().unwrap_or(0.0),
                ea: parse_energy(low_node, ea_unit)?,
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
                    let eff  = v.as_f64()?;
                    let idx  = sp_index(name)?;
                    Some((idx, eff))
                }).collect()
            } else {
                vec![]
            };
            Some(ThirdBodySpec { efficiencies: effs })
        } else {
            None
        };

        out.push(Reaction { reactants, products, rate, reversible, third_body, duplicate });
    }
    Ok(out)
}

fn parse_stoich(
    side: &str,
    sp_index: &impl Fn(&str) -> Option<usize>,
) -> Result<Vec<(usize, f64)>> {
    let mut out = Vec::new();
    // Remove third-body markers like "(+M)" or "(+m)"
    let side = side.replace("(+M)", "").replace("(+m)", "");
    for token in side.split('+') {
        let token = token.trim();
        if token.is_empty() || token == "M" || token == "m" {
            continue;
        }
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
    let trimmed = token.trim();
    let end = trimmed
        .find(|c: char| c.is_alphabetic())
        .unwrap_or(0);
    if end == 0 {
        (1.0, trimmed)
    } else {
        let coeff_str = trimmed[..end].trim();
        let name      = trimmed[end..].trim();
        let coeff     = coeff_str.parse::<f64>().unwrap_or(1.0);
        (coeff, name)
    }
}

/// Parse a pressure value that may be a bare float (Pa) or a string with units,
/// e.g. `1.0 atm`, `101325`, `1 bar`.
fn parse_pressure(node: &serde_yaml::Value) -> Result<f64> {
    if let Some(v) = node.as_f64() {
        return Ok(v); // bare number → Pa
    }
    let s = node.as_str()
        .ok_or_else(|| anyhow::anyhow!("PLOG pressure is not a number or string"))?;
    let s = s.trim();
    // Try to split into value + unit
    let (val_str, unit) = if let Some(pos) = s.find(|c: char| c.is_alphabetic()) {
        (s[..pos].trim(), s[pos..].trim())
    } else {
        (s, "Pa")
    };
    let val: f64 = val_str.parse()
        .map_err(|_| anyhow::anyhow!("Cannot parse pressure value '{val_str}'"))?;
    let pa = match unit.to_lowercase().as_str() {
        "pa"   => val,
        "atm"  => val * 101325.0,
        "bar"  => val * 1.0e5,
        "torr" => val * 133.322,
        other  => anyhow::bail!("Unknown pressure unit '{other}'"),
    };
    Ok(pa)
}

/// Parse activation energy from a rate node and convert to J/mol.
/// `ea_unit` is taken from the top-level `units: {activation-energy: ...}` block.
fn parse_energy(node: &serde_yaml::Value, ea_unit: &str) -> Result<f64> {
    let ea_raw = node["Ea"].as_f64().unwrap_or(0.0);
    let ea_j_mol = match ea_unit {
        "cal/mol" | "cal"  => ea_raw * 4.184,       // cal/mol → J/mol
        "kcal/mol"         => ea_raw * 4184.0,
        "J/mol"            => ea_raw,
        "kJ/mol"           => ea_raw * 1000.0,
        "K"                => ea_raw * 8.314462618,  // K * R
        "J/kmol"           => ea_raw * 1e-3,
        _                  => ea_raw * 4.184,        // assume cal/mol
    };
    Ok(ea_j_mol)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Path to the bundled Cantera H2/O2 mechanism.
    fn h2o2_path() -> String {
        // Relative to the crate root (where Cargo.toml lives)
        let manifest = env!("CARGO_MANIFEST_DIR");
        format!("{manifest}/data/h2o2.yaml")
    }

    #[test]
    fn test_h2o2_species_count() {
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        // h2o2.yaml contains: H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2
        assert_eq!(mech.n_species(), 10, "expected 10 species");
    }

    #[test]
    fn test_h2o2_reaction_count() {
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        // h2o2.yaml has 29 reactions
        assert_eq!(mech.n_reactions(), 29, "expected 29 reactions");
    }

    #[test]
    fn test_h2o2_molecular_weights() {
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        let sp = |name: &str| mech.species.iter().find(|s| s.name == name).unwrap();

        // H2: 2 * 1.008e-3 = 2.016e-3 kg/mol
        approx::assert_abs_diff_eq!(sp("H2").molecular_weight, 2.016e-3, epsilon = 1e-6);
        // O2: 2 * 15.999e-3 = 31.998e-3
        approx::assert_abs_diff_eq!(sp("O2").molecular_weight, 31.998e-3, epsilon = 1e-6);
        // H2O: 2*1.008e-3 + 15.999e-3 = 18.015e-3
        approx::assert_abs_diff_eq!(sp("H2O").molecular_weight, 18.015e-3, epsilon = 1e-6);
        // AR: 39.948e-3 (case-insensitive: YAML uses {Ar: 1})
        approx::assert_abs_diff_eq!(sp("AR").molecular_weight, 39.948e-3, epsilon = 1e-6);
        // N2: 2 * 14.007e-3 = 28.014e-3
        approx::assert_abs_diff_eq!(sp("N2").molecular_weight, 28.014e-3, epsilon = 1e-6);
    }

    #[test]
    fn test_h2o2_activation_energy_units() {
        // Units block says cal/mol; Ea should be converted to J/mol.
        // Reaction 3: O + H2 <=> H + OH, Ea = 6260 cal/mol
        // Expected: 6260 * 4.184 = 26191.84 J/mol
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        let rxn = &mech.reactions[2]; // 0-indexed
        if let crate::chemistry::mechanism::RateType::Arrhenius { ea, .. } = rxn.rate {
            approx::assert_abs_diff_eq!(ea, 6260.0 * 4.184, epsilon = 1e-3);
        } else {
            panic!("Reaction 3 should be Arrhenius");
        }
    }

    #[test]
    fn test_h2o2_duplicate_flag() {
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        // Reactions 24–29 (0-indexed 23–28) are marked duplicate
        for idx in 23..29 {
            assert!(
                mech.reactions[idx].duplicate,
                "reaction {idx} should have duplicate=true"
            );
        }
        // Reaction 1 (idx 0) is not duplicate
        assert!(!mech.reactions[0].duplicate, "reaction 0 should not be duplicate");
    }

    #[test]
    fn test_h2o2_falloff_reaction() {
        // Reaction 22 (idx 21): 2 OH (+M) <=> H2O2 (+M), type falloff with Troe
        let mech = parse_file(&h2o2_path()).expect("failed to parse h2o2.yaml");
        let rxn = &mech.reactions[21];
        match &rxn.rate {
            crate::chemistry::mechanism::RateType::Falloff { troe: Some(t), high, low } => {
                approx::assert_abs_diff_eq!(t.a, 0.7346, epsilon = 1e-6);
                approx::assert_abs_diff_eq!(t.t3, 94.0, epsilon = 1e-6);
                approx::assert_abs_diff_eq!(t.t1, 1756.0, epsilon = 1e-3);
                // high-P: A=7.4e13, b=-0.37, Ea=0
                approx::assert_abs_diff_eq!(high.a, 7.4e13, epsilon = 1e6);
                // low-P: A=2.3e18, b=-0.9, Ea=-1700 cal/mol → -7112.8 J/mol
                approx::assert_abs_diff_eq!(low.ea, -1700.0 * 4.184, epsilon = 1e-3);
            }
            _ => panic!("Reaction 22 should be Falloff with Troe"),
        }
        assert!(rxn.third_body.is_some(), "Reaction 22 should have third body");
    }

    #[test]
    fn test_normalise_element() {
        assert_eq!(normalise_element("AR"), "Ar");
        assert_eq!(normalise_element("ar"), "Ar");
        assert_eq!(normalise_element("Ar"), "Ar");
        assert_eq!(normalise_element("H"),  "H");
        assert_eq!(normalise_element("he"), "He");
        assert_eq!(normalise_element("CL"), "Cl");
    }

    #[test]
    fn test_plog_parse() {
        // Minimal synthetic PLOG mechanism
        let yaml = r#"
units: {activation-energy: cal/mol}
species:
- name: H
  composition: {H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, 25473.0, -0.447]
    - [2.5, 0.0, 0.0, 0.0, 0.0, 25473.0, -0.447]
  transport:
    model: gas
    geometry: atom
    well-depth: 145.0
    diameter: 2.05
- name: O
  composition: {O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, 29230.0, 4.92]
    - [2.5, 0.0, 0.0, 0.0, 0.0, 29230.0, 4.92]
  transport:
    model: gas
    geometry: atom
    well-depth: 80.0
    diameter: 2.75
- name: OH
  composition: {O: 1, H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [4.125, -3.225e-03, 6.527e-06, -5.798e-09, 2.062e-12, 3346.0, -0.69]
    - [2.869, 1.013e-03, -2.276e-07, 2.174e-11, -5.126e-16, 3886.0, 5.70]
  transport:
    model: gas
    geometry: linear
    well-depth: 80.0
    diameter: 2.75
reactions:
- equation: H + O <=> OH
  pressure-dependent-Arrhenius:
  - {P: 1.0 atm, A: 1.0e+13, b: 0.0, Ea: 1000.0}
  - {P: 10.0 atm, A: 2.0e+13, b: 0.0, Ea: 1200.0}
"#;
        let mech = parse_str(yaml).expect("failed to parse PLOG yaml");
        assert_eq!(mech.n_reactions(), 1);
        match &mech.reactions[0].rate {
            crate::chemistry::mechanism::RateType::Plog { rates } => {
                assert_eq!(rates.len(), 2);
                approx::assert_abs_diff_eq!(rates[0].0, 101325.0, epsilon = 1.0);
                approx::assert_abs_diff_eq!(rates[1].0, 1013250.0, epsilon = 10.0);
                approx::assert_abs_diff_eq!(rates[0].1.ea, 1000.0 * 4.184, epsilon = 1e-3);
            }
            _ => panic!("Expected Plog rate type"),
        }
    }
}
