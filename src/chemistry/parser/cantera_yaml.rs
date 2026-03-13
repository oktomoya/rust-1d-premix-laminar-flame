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

    // Read activation-energy unit from the top-level units: block (default cal/mol).
    let ea_units = doc["units"]["activation-energy"]
        .as_str()
        .unwrap_or("cal/mol")
        .to_string();

    let species = parse_species(&doc)?;
    let reactions = parse_reactions(&doc, &species, &ea_units)?;
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

        // Molecular weight from composition; also store the element map
        let (mw, composition) = molecular_weight_from_composition(&sp_node["composition"])?;

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
            composition,
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
/// Also returns the normalised element map (upper-case keys, atom counts).
/// Element lookup is case-insensitive (normalised to upper-case).
fn molecular_weight_from_composition(
    comp: &serde_yaml::Value,
) -> Result<(f64, std::collections::HashMap<String, f64>)> {
    let atom_weights: std::collections::HashMap<&str, f64> = [
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
    ].into();

    let map = comp.as_mapping()
        .ok_or_else(|| anyhow::anyhow!("Composition is not a mapping"))?;
    let mut mw = 0.0;
    let mut composition = std::collections::HashMap::new();
    for (el, count) in map {
        let el_raw = el.as_str().unwrap_or("");
        let el_upper = el_raw.to_uppercase();
        let w = atom_weights.get(el_upper.as_str())
            .ok_or_else(|| anyhow::anyhow!("Unknown element: {el_raw}"))?;
        let n = count.as_f64().unwrap_or(0.0);
        mw += n * w;
        composition.insert(el_upper, n);
    }
    Ok((mw, composition))
}

fn parse_reactions(
    doc: &serde_yaml::Value,
    species: &[crate::chemistry::species::Species],
    ea_units: &str,
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

        // duplicate flag
        let duplicate = rxn_node["duplicate"].as_bool().unwrap_or(false);

        // Determine reaction type
        let rxn_type = rxn_node["type"].as_str().unwrap_or("");

        let rate = if rxn_type == "pressure-dependent-Arrhenius" {
            // PLOG reaction
            let rate_list = rxn_node["rate-constants"]
                .as_sequence()
                .ok_or_else(|| anyhow::anyhow!("PLOG reaction missing rate-constants: {eq}"))?;
            let mut rates = Vec::new();
            for entry in rate_list {
                let p_pa = parse_pressure(entry)?;
                let a  = entry["A"].as_f64().unwrap_or(0.0);
                let b  = entry["b"].as_f64().unwrap_or(0.0);
                let ea = parse_energy(entry, ea_units)?;
                rates.push((p_pa, Arrhenius { a, b, ea }));
            }
            RateType::Plog { rates }
        } else {
            let rate_node = &rxn_node["rate-constant"];
            if rate_node.is_null() {
                // Falloff
                let high_node = &rxn_node["high-P-rate-constant"];
                let low_node  = &rxn_node["low-P-rate-constant"];
                let high = Arrhenius {
                    a:  high_node["A"].as_f64().unwrap_or(0.0),
                    b:  high_node["b"].as_f64().unwrap_or(0.0),
                    ea: parse_energy(high_node, ea_units)?,
                };
                let low = Arrhenius {
                    a:  low_node["A"].as_f64().unwrap_or(0.0),
                    b:  low_node["b"].as_f64().unwrap_or(0.0),
                    ea: parse_energy(low_node, ea_units)?,
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
            } else {
                let a  = rate_node["A"].as_f64().unwrap_or(0.0);
                let b  = rate_node["b"].as_f64().unwrap_or(0.0);
                let ea = parse_energy(rate_node, ea_units)?;
                RateType::Arrhenius { a, b, ea }
            }
        };

        // Third body: explicit type=three-body, or "(+M)"/"+ M" in equation
        let third_body = if rxn_type == "three-body"
            || eq.contains("(+M)")
            || eq.contains("+ M")
        {
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

        out.push(Reaction { reactants, products, rate, reversible, third_body, duplicate });
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

/// Parse activation energy from a rate-constant node and convert to J/mol.
/// `ea_units` comes from the top-level `units: activation-energy:` field.
fn parse_energy(node: &serde_yaml::Value, ea_units: &str) -> Result<f64> {
    let ea_raw = node["Ea"].as_f64().unwrap_or(0.0);
    let ea_j_mol = match ea_units {
        "cal/mol" | "cal"  => ea_raw * 4.184,
        "kcal/mol"         => ea_raw * 4184.0,
        "J/mol"            => ea_raw,
        "kJ/mol"           => ea_raw * 1000.0,
        "K"                => ea_raw * 8.314462618,
        "J/kmol"           => ea_raw * 1e-3,
        _                  => ea_raw * 4.184,  // default: cal/mol
    };
    Ok(ea_j_mol)
}

/// Parse pressure from a PLOG rate-constants entry.
/// Accepts a plain float (assumed Pa) or a string "value unit".
fn parse_pressure(node: &serde_yaml::Value) -> Result<f64> {
    let p_val = &node["P"];
    if let Some(p) = p_val.as_f64() {
        return Ok(p);
    }
    if let Some(s) = p_val.as_str() {
        let parts: Vec<&str> = s.split_whitespace().collect();
        if parts.len() == 2 {
            let val: f64 = parts[0].parse()
                .map_err(|_| anyhow::anyhow!("Cannot parse pressure value: {s}"))?;
            let pa = match parts[1] {
                "atm"  => val * 101325.0,
                "bar"  => val * 1e5,
                "kPa"  => val * 1e3,
                "MPa"  => val * 1e6,
                "Pa"   => val,
                other  => anyhow::bail!("Unknown pressure unit: {other}"),
            };
            return Ok(pa);
        }
    }
    anyhow::bail!("Cannot parse pressure field in PLOG entry")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h2o2_path() -> String {
        let manifest = env!("CARGO_MANIFEST_DIR");
        format!("{manifest}/data/h2o2.yaml")
    }

    #[test]
    fn test_parse_h2o2_species_count() {
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        // h2o2.yaml has 10 species: H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2
        assert_eq!(mech.n_species(), 10);
    }

    #[test]
    fn test_parse_h2o2_reaction_count() {
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        // h2o2.yaml has 29 reactions
        assert_eq!(mech.n_reactions(), 29);
    }

    #[test]
    fn test_molecular_weights() {
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        let sp_idx = |name: &str| mech.species_index(name).expect(name);

        // H2: 2 * 1.008e-3
        let mw_h2 = mech.species[sp_idx("H2")].molecular_weight;
        assert!((mw_h2 - 2.016e-3).abs() < 1e-6, "H2 mw = {mw_h2}");

        // O2: 2 * 15.999e-3
        let mw_o2 = mech.species[sp_idx("O2")].molecular_weight;
        assert!((mw_o2 - 31.998e-3).abs() < 1e-6, "O2 mw = {mw_o2}");

        // AR (Ar in YAML): 39.948e-3
        let mw_ar = mech.species[sp_idx("AR")].molecular_weight;
        assert!((mw_ar - 39.948e-3).abs() < 1e-6, "AR mw = {mw_ar}");

        // N2: 2 * 14.007e-3
        let mw_n2 = mech.species[sp_idx("N2")].molecular_weight;
        assert!((mw_n2 - 28.014e-3).abs() < 1e-6, "N2 mw = {mw_n2}");
    }

    #[test]
    fn test_activation_energy_units_cal_mol() {
        // Reaction 3: O + H2 <=> H + OH, Ea = 6260 cal/mol → 26191.84 J/mol
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        let rxn = &mech.reactions[2];
        if let crate::chemistry::mechanism::RateType::Arrhenius { ea, .. } = rxn.rate {
            let expected = 6260.0 * 4.184;
            assert!((ea - expected).abs() < 1.0, "Ea = {ea}, expected {expected}");
        } else {
            panic!("Reaction 3 should be Arrhenius");
        }
    }

    #[test]
    fn test_duplicate_flag() {
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        // Reactions 24-29 (indices 23-28) are marked duplicate
        for i in 23..29 {
            assert!(mech.reactions[i].duplicate,
                "Reaction {} should have duplicate=true", i + 1);
        }
        // Reaction 1 (index 0) is not duplicate
        assert!(!mech.reactions[0].duplicate);
    }

    #[test]
    fn test_falloff_reaction_22() {
        // Reaction 22: 2 OH (+M) <=> H2O2 (+M), type falloff with Troe
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        let rxn = &mech.reactions[21];
        assert!(rxn.reversible);
        match &rxn.rate {
            crate::chemistry::mechanism::RateType::Falloff { high, low, troe } => {
                assert!((high.a - 7.4e13).abs() / 7.4e13 < 1e-6,
                    "high.A = {}", high.a);
                assert!((low.a - 2.3e18).abs() / 2.3e18 < 1e-6,
                    "low.A = {}", low.a);
                let tr = troe.as_ref().expect("Troe params");
                assert!((tr.a - 0.7346).abs() < 1e-6, "Troe.A = {}", tr.a);
            }
            other => panic!("Expected Falloff, got {other:?}"),
        }
    }

    #[test]
    fn test_three_body_reaction_1() {
        // Reaction 1: 2 O + M <=> O2 + M, type three-body
        let mech = parse_file(&h2o2_path()).expect("parse h2o2.yaml");
        let rxn = &mech.reactions[0];
        assert!(rxn.third_body.is_some(), "Reaction 1 should have third body");
        let tb = rxn.third_body.as_ref().unwrap();
        // efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
        assert!(!tb.efficiencies.is_empty());
    }

    #[test]
    fn test_case_insensitive_elements() {
        // Parse a minimal YAML with lowercase element names
        let yaml = r#"
units: {activation-energy: cal/mol}
species:
- name: TEST
  composition: {h: 2, o: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.0, 0.0, 0.0, 0.0, 0.0, -1000.0, 1.0]
    - [3.0, 0.0, 0.0, 0.0, 0.0, -1000.0, 1.0]
  transport:
    model: gas
    geometry: nonlinear
    well-depth: 100.0
    diameter: 3.0
reactions: []
"#;
        let mech = parse_str(yaml).expect("parse lowercase elements");
        let mw = mech.species[0].molecular_weight;
        // H2O: 2*1.008e-3 + 15.999e-3 = 18.015e-3
        assert!((mw - 18.015e-3).abs() < 1e-6, "mw = {mw}");
    }

    #[test]
    fn test_plog_parsing() {
        let yaml = r#"
units: {activation-energy: cal/mol}
species:
- name: H
  composition: {H: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, 25000.0, -0.45]
    - [2.5, 0.0, 0.0, 0.0, 0.0, 25000.0, -0.45]
  transport:
    model: gas
    geometry: atom
    well-depth: 145.0
    diameter: 2.05
- name: O2
  composition: {O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.78, -3.0e-03, 9.8e-06, -9.7e-09, 3.3e-12, -1063.9, 3.66]
    - [3.28, 1.48e-03, -7.58e-07, 2.09e-10, -2.17e-14, -1088.5, 5.45]
  transport:
    model: gas
    geometry: linear
    well-depth: 107.4
    diameter: 3.46
- name: HO2
  composition: {H: 1, O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [4.3, -4.7e-03, 2.1e-05, -2.4e-08, 9.3e-12, 264.0, 3.71]
    - [4.02, 2.24e-03, -6.34e-07, 1.14e-10, -1.08e-14, 111.9, 3.79]
  transport:
    model: gas
    geometry: nonlinear
    well-depth: 107.4
    diameter: 3.46
reactions:
- equation: H + O2 (+M) <=> HO2 (+M)
  type: pressure-dependent-Arrhenius
  rate-constants:
  - {P: 101325.0, A: 1.0e+13, b: 0.0, Ea: 0.0}
  - {P: 1013250.0, A: 2.0e+13, b: 0.0, Ea: 0.0}
"#;
        let mech = parse_str(yaml).expect("parse PLOG");
        assert_eq!(mech.n_reactions(), 1);
        match &mech.reactions[0].rate {
            crate::chemistry::mechanism::RateType::Plog { rates } => {
                assert_eq!(rates.len(), 2);
                assert!((rates[0].0 - 101325.0).abs() < 1.0);
                assert!((rates[1].0 - 1013250.0).abs() < 1.0);
            }
            other => panic!("Expected Plog, got {other:?}"),
        }
    }
}
