use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{enthalpy_molar, entropy_species, R_UNIVERSAL};

/// Arrhenius rate coefficient: k = A * T^b * exp(-Ea / (R*T))
/// Units depend on reaction order; Ea in J/mol.
pub fn arrhenius(a: f64, b: f64, ea: f64, t: f64) -> f64 {
    a * t.powf(b) * (-ea / (R_UNIVERSAL * t)).exp()
}

/// Equilibrium constant Kc (concentration-based) for reaction i.
/// Kc = Kp * (P_atm / (R*T))^(Σνk_products - Σνk_reactants)
fn equilibrium_constant(mech: &Mechanism, rxn_idx: usize, t: f64) -> f64 {
    let rxn = &mech.reactions[rxn_idx];
    let p_atm = 101325.0_f64; // Pa

    // ΔG°/RT = Σ νk * (g_k / RT) = Σ νk * (h_k/RT - s_k/R)
    let mut delta_g_over_rt = 0.0_f64;
    let mut delta_nu = 0.0_f64;

    for &(k, nu) in &rxn.products {
        let sp = &mech.species[k];
        let h = enthalpy_molar(sp, t);
        let s = entropy_species(sp, t) * sp.molecular_weight; // back to J/(mol·K)
        delta_g_over_rt += nu * (h / (R_UNIVERSAL * t) - s / R_UNIVERSAL);
        delta_nu += nu;
    }
    for &(k, nu) in &rxn.reactants {
        let sp = &mech.species[k];
        let h = enthalpy_molar(sp, t);
        let s = entropy_species(sp, t) * sp.molecular_weight;
        delta_g_over_rt -= nu * (h / (R_UNIVERSAL * t) - s / R_UNIVERSAL);
        delta_nu -= nu;
    }

    // Kp = exp(-ΔG°/RT)
    let kp = (-delta_g_over_rt).exp();
    // Kc = Kp * (c° )^delta_nu  where c° = P_atm/(R*T) in mol/cm³ (CGS, matching stored A values)
    // P_atm/(R*T) gives mol/m³; multiply by 1e-6 to convert to mol/cm³.
    kp * (p_atm / (R_UNIVERSAL * t) * 1e-6).powf(delta_nu)
}

/// Evaluate molar production rates ωk [mol/(m³·s)] for all species.
/// `concentrations` [mol/m³]: c[k] = ρ * Yk / Wk
/// `pressure` [Pa]: required for PLOG pressure-dependent reactions.
pub fn production_rates(mech: &Mechanism, t: f64, concentrations: &[f64], pressure: f64) -> Vec<f64> {
    let nk = mech.species.len();
    let mut wdot = vec![0.0_f64; nk];

    for (i, rxn) in mech.reactions.iter().enumerate() {
        let kf = match &rxn.rate {
            crate::chemistry::mechanism::RateType::Arrhenius { a, b, ea } => {
                arrhenius(*a, *b, *ea, t)
            }
            crate::chemistry::mechanism::RateType::Falloff { high, low, troe } => {
                falloff_rate(high, low, troe.as_ref(), t, concentrations, mech, i)
            }
            crate::chemistry::mechanism::RateType::Plog { rates } => {
                plog_rate(rates, t, pressure)
            }
        };

        // Third-body efficiency.
        // Falloff reactions already incorporate [M] inside falloff_rate() via Pr,
        // so they must NOT be multiplied by m_conc again.
        // Only plain three-body Arrhenius reactions need the extra [M] factor.
        let kf_eff = match &rxn.rate {
            crate::chemistry::mechanism::RateType::Falloff { .. } => kf,
            _ => {
                if mech.reactions[i].third_body.is_some() {
                    let m_conc = third_body_concentration(mech, i, concentrations);
                    kf * m_conc
                } else {
                    kf
                }
            }
        };

        let kc = equilibrium_constant(mech, i, t);
        // Use kf_eff (includes [M] for three-body) so forward and reverse are symmetric.
        let kr = kf_eff / kc.max(1e-300);

        // Forward rate of progress
        let qf: f64 = rxn.reactants.iter()
            .map(|&(k, nu)| concentrations[k].max(0.0).powf(nu))
            .product();
        // Reverse rate of progress
        let qr: f64 = rxn.products.iter()
            .map(|&(k, nu)| concentrations[k].max(0.0).powf(nu))
            .product();

        let q = kf_eff * qf - kr * qr;

        for &(k, nu) in &rxn.products {
            wdot[k] += nu * q;
        }
        for &(k, nu) in &rxn.reactants {
            wdot[k] -= nu * q;
        }
    }

    wdot
}

/// PLOG rate: log-linear interpolation between pressure-bracketing entries.
fn plog_rate(rates: &[(f64, crate::chemistry::mechanism::Arrhenius)], t: f64, pressure: f64) -> f64 {
    if rates.is_empty() {
        return 0.0;
    }
    // Clamp to lowest/highest pressure entry
    if pressure <= rates[0].0 {
        let r = &rates[0].1;
        return arrhenius(r.a, r.b, r.ea, t);
    }
    if pressure >= rates[rates.len() - 1].0 {
        let r = &rates[rates.len() - 1].1;
        return arrhenius(r.a, r.b, r.ea, t);
    }
    // Find bracketing pair
    let idx = rates.partition_point(|&(p, _)| p < pressure);
    let (p1, r1) = &rates[idx - 1];
    let (p2, r2) = &rates[idx];
    let k1 = arrhenius(r1.a, r1.b, r1.ea, t).max(1e-300);
    let k2 = arrhenius(r2.a, r2.b, r2.ea, t).max(1e-300);
    let log_k = k1.ln() + (pressure.ln() - p1.ln()) * (k2.ln() - k1.ln()) / (p2.ln() - p1.ln());
    log_k.exp()
}

fn falloff_rate(
    high: &crate::chemistry::mechanism::Arrhenius,
    low: &crate::chemistry::mechanism::Arrhenius,
    troe: Option<&crate::chemistry::mechanism::TroeParams>,
    t: f64,
    concentrations: &[f64],
    mech: &Mechanism,
    rxn_idx: usize,
) -> f64 {
    let k_inf = arrhenius(high.a, high.b, high.ea, t);
    let k0 = arrhenius(low.a, low.b, low.ea, t);
    let m = third_body_concentration(mech, rxn_idx, concentrations);
    let pr = k0 * m / k_inf.max(1e-300);

    let f = if let Some(tr) = troe {
        troe_broadening(tr, t, pr)
    } else {
        1.0 // Lindemann
    };

    k_inf * (pr / (1.0 + pr)) * f
}

fn troe_broadening(troe: &crate::chemistry::mechanism::TroeParams, t: f64, pr: f64) -> f64 {
    let f_cent = (1.0 - troe.a) * (-t / troe.t3).exp()
        + troe.a * (-t / troe.t1).exp()
        + troe.t2.map(|t2| (-t2 / t).exp()).unwrap_or(0.0);
    let log_f_cent = f_cent.max(1e-300).log10();
    let c = -0.4 - 0.67 * log_f_cent;
    let n = 0.75 - 1.27 * log_f_cent;
    let d = 0.14;
    let log_pr = pr.max(1e-300).log10();
    let f1 = (log_pr + c) / (n - d * (log_pr + c));
    10.0_f64.powf(log_f_cent / (1.0 + f1 * f1))
}

fn third_body_concentration(mech: &Mechanism, rxn_idx: usize, concentrations: &[f64]) -> f64 {
    let m_default: f64 = concentrations.iter().sum();
    match &mech.reactions[rxn_idx].third_body {
        None => 0.0,
        Some(tb) => {
            let mut m = m_default;
            for &(k, eff) in &tb.efficiencies {
                m += (eff - 1.0) * concentrations[k];
            }
            m
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;

    fn h2o2_mech() -> Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    fn check(label: &str, got: f64, expected: f64, rtol: f64) {
        let rel = (got - expected).abs() / expected.abs().max(1e-300);
        assert!(
            rel < rtol,
            "{label}: got {got:.6e}, expected {expected:.6e}, rel err {rel:.2e}"
        );
    }

    // -----------------------------------------------------------------------
    // 1. Equilibrium constant Kc
    //
    // Reaction 11 (index 10): H + O2 <=> O + OH  Δν = 0 (dimensionless Kc)
    // Reference: Cantera 3.1.0 equilibrium_constants[10]
    // -----------------------------------------------------------------------
    #[test]
    fn test_kc_rxn11_1000k() {
        let mech = h2o2_mech();
        let kc = equilibrium_constant(&mech, 10, 1000.0);
        check("Kc rxn11 1000K", kc, 3.979362891605e-3, 1e-6);
    }

    #[test]
    fn test_kc_rxn11_2000k() {
        let mech = h2o2_mech();
        let kc = equilibrium_constant(&mech, 10, 2000.0);
        check("Kc rxn11 2000K", kc, 2.310976953973e-1, 1e-6);
    }

    // -----------------------------------------------------------------------
    // 2. Troe broadening factor F
    //
    // Reaction 22 (index 21): 2 OH (+M) <=> H2O2 (+M)
    // Troe: A=0.7346, T3=94, T1=1756, T2=5182
    // Reference: Python implementation of Cantera's TroeRate::updateTemp formula
    // -----------------------------------------------------------------------
    #[test]
    fn test_troe_broadening() {
        let troe = crate::chemistry::mechanism::TroeParams {
            a: 0.7346, t3: 94.0, t1: 1756.0, t2: Some(5182.0),
        };

        // (T, Pr, expected_F)
        let cases = [
            (500.0,  0.1,  0.74076387999260_f64),
            (500.0,  1.0,  0.56737659656155),
            (500.0,  10.0, 0.69693521179172),
            (1000.0, 0.1,  0.59868084443837),
            (1000.0, 1.0,  0.42639156438466),
            (1000.0, 10.0, 0.58081881002988),
            (1500.0, 0.1,  0.49907753234010),
            (1500.0, 1.0,  0.34587257114316),
            (1500.0, 10.0, 0.50597956826065),
        ];
        for (t, pr, expected) in cases {
            let f = troe_broadening(&troe, t, pr);
            check(&format!("Troe F T={t} Pr={pr}"), f, expected, 1e-10);
        }
    }

    // -----------------------------------------------------------------------
    // 3. Molar production rates ωk
    //
    // State: T=1500 K, P=101325 Pa, X = H2:2, O2:1, N2:3.76 (no radicals)
    //
    // Unit convention: A values in the YAML are in CGS (cm³/mol/s).
    // Concentrations must therefore be supplied in mol/cm³ so that
    // k × c₁ × c₂ yields mol/(cm³·s).
    //
    // Reference concentrations and wdot from Cantera 3.1.0 (converted to CGS):
    //   conc [mol/cm³] = Cantera [kmol/m³] × 1e-3
    //   wdot [mol/(cm³·s)] = Cantera [kmol/(m³·s)] × 1e-3
    // -----------------------------------------------------------------------
    #[test]
    fn test_production_rates_h2o2_state() {
        let mech = h2o2_mech();
        let nk = mech.n_species();

        // Concentrations in mol/cm³ (CGS, consistent with stored A values)
        let mut conc = vec![0.0_f64; nk];
        let idx = |name: &str| mech.species_index(name).unwrap();
        conc[idx("H2")] = 2.403667924005e-6;
        conc[idx("O2")] = 1.201833962002e-6;
        conc[idx("N2")] = 4.518895697129e-6;

        let wdot = production_rates(&mech, 1500.0, &conc, 101325.0);

        // Dominant reaction at this state: reverse of rxn 17 (H + HO2 <=> H2 + O2)
        // → H2 + O2 → H + HO2
        check("ω_H2",  wdot[idx("H2")],  -2.297824021187e-6, 1e-4);
        check("ω_H",   wdot[idx("H")],    2.297863386999e-6, 1e-4);
        check("ω_O2",  wdot[idx("O2")],  -2.297784802860e-6, 1e-4);
        check("ω_HO2", wdot[idx("HO2")],  2.297784655376e-6, 1e-4);

        // Species absent in initial mixture should have near-zero production
        // (H2O, H2O2 are not produced at this state)
        assert!(wdot[idx("H2O")].abs()  < 1e-20, "ω_H2O  = {}", wdot[idx("H2O")]);
        assert!(wdot[idx("H2O2")].abs() < 1e-20, "ω_H2O2 = {}", wdot[idx("H2O2")]);
    }
}
