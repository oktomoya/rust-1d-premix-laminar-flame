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
    // Kc = Kp * (P_atm / (R*T))^delta_nu
    kp * (p_atm / (R_UNIVERSAL * t)).powf(delta_nu)
}

/// Evaluate molar production rates ωk [mol/(m³·s)] for all species.
/// `concentrations` [mol/m³]: c[k] = ρ * Yk / Wk
pub fn production_rates(mech: &Mechanism, t: f64, concentrations: &[f64]) -> Vec<f64> {
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
        };

        // Third-body efficiency
        let m_conc = third_body_concentration(mech, i, concentrations);
        let kf_eff = if mech.reactions[i].third_body.is_some() {
            kf * m_conc
        } else {
            kf
        };

        let kc = equilibrium_constant(mech, i, t);
        let kr = kf / kc.max(1e-300);

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
