/// Mixture-averaged transport properties.

use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{cp_species, mean_molecular_weight};
use crate::transport::species_props::{binary_diffusion, thermal_conductivity, viscosity};

/// Mixture viscosity [Pa·s] via Wilke mixing rule.
pub fn mixture_viscosity(mech: &Mechanism, mole_fractions: &[f64], t: f64) -> f64 {
    let nk = mech.n_species();
    let mu: Vec<f64> = mech.species.iter().map(|s| viscosity(s, t)).collect();

    let mut mu_mix = 0.0_f64;
    for k in 0..nk {
        if mole_fractions[k] < 1e-20 {
            continue;
        }
        let mut phi_k = 0.0_f64;
        for j in 0..nk {
            if mole_fractions[j] < 1e-20 {
                continue;
            }
            let w_ratio = (mech.species[j].molecular_weight / mech.species[k].molecular_weight).powf(0.25);
            let phi_kj = (1.0 + (mu[k] / mu[j]).sqrt() * w_ratio).powi(2)
                / (8.0 * (1.0 + mech.species[k].molecular_weight / mech.species[j].molecular_weight)).sqrt();
            phi_k += mole_fractions[j] * phi_kj;
        }
        mu_mix += mole_fractions[k] * mu[k] / phi_k;
    }
    mu_mix
}

/// Mixture thermal conductivity [W/(m·K)] using the arithmetic-harmonic mean.
pub fn mixture_thermal_conductivity(
    mech: &Mechanism,
    mole_fractions: &[f64],
    mass_fractions: &[f64],
    t: f64,
    pressure: f64,
) -> f64 {
    let _ = mass_fractions;
    let nk = mech.n_species();
    let mut sum1 = 0.0_f64;
    let mut sum2 = 0.0_f64;
    for k in 0..nk {
        let mu_k = viscosity(&mech.species[k], t);
        let cp_k = cp_species(&mech.species[k], t);
        let lam_k = thermal_conductivity(&mech.species[k], mu_k, cp_k, t, pressure);
        sum1 += mole_fractions[k] * lam_k;
        sum2 += mole_fractions[k] / lam_k.max(1e-30);
    }
    0.5 * (sum1 + 1.0 / sum2)
}

/// Mixture-averaged diffusion coefficients Dkm [m²/s] for each species k.
///
/// Matches Cantera's `getMixDiffCoeffs` (used with molar flux basis):
///   Dkm = (1 - Yk) / Σ_{j≠k} (Xj / Dkj)
///
/// where Yk is the mass fraction. Cantera's formula is
/// `(W̄ - Xk·Wk) / (W̄ · Σ Xj/Dkj) = (1 - Yk) / Σ Xj/Dkj`.
pub fn mixture_diffusion_coefficients(
    mech: &Mechanism,
    mole_fractions: &[f64],
    mass_fractions: &[f64],
    t: f64,
    pressure: f64,
) -> Vec<f64> {
    let nk = mech.n_species();
    let mut dkm = vec![0.0_f64; nk];

    for k in 0..nk {
        let yk = mass_fractions[k];
        let sum: f64 = (0..nk)
            .filter(|&j| j != k)
            .map(|j| {
                let dkj = binary_diffusion(&mech.species[k], &mech.species[j], t, pressure);
                mole_fractions[j] / dkj.max(1e-30)
            })
            .sum();
        dkm[k] = if sum > 1e-30 {
            (1.0 - yk).max(0.0) / sum
        } else {
            0.0
        };
    }
    dkm
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
        let rel = (got - expected).abs() / expected.abs();
        assert!(rel < rtol, "{label}: got {got:.6e}, expected {expected:.6e}, rel err {rel:.2e}");
    }

    fn air_mole_fractions(mech: &Mechanism) -> Vec<f64> {
        let nk = mech.n_species();
        let mut x = vec![0.0_f64; nk];
        x[mech.species_index("O2").unwrap()] = 0.21;
        x[mech.species_index("N2").unwrap()] = 0.79;
        x
    }

    // -----------------------------------------------------------------------
    // Wilke mixture viscosity for air (O2:0.21, N2:0.79 mole fractions)
    // Reference: Cantera 3.1.0 (same h2o2.yaml transport parameters).
    // Residual ~0.06% vs Cantera; 0.5% tolerance is ample.
    // -----------------------------------------------------------------------
    #[test]
    fn test_wilke_air_300k() {
        let mech = h2o2_mech();
        let x = air_mole_fractions(&mech);
        check("mu_air 300K", mixture_viscosity(&mech, &x, 300.0), 1.863048267765e-5, 5e-3);
    }

    #[test]
    fn test_wilke_air_1000k() {
        let mech = h2o2_mech();
        let x = air_mole_fractions(&mech);
        check("mu_air 1000K", mixture_viscosity(&mech, &x, 1000.0), 4.285066600992e-5, 5e-3);
    }

    // -----------------------------------------------------------------------
    // Wilke mixture viscosity for H2/N2 (X_H2=0.5, X_N2=0.5)
    // Large MW contrast (W_H2=2 vs W_N2=28) exposes the wrong-exponent bug
    // that was invisible in the air test (W_O2≈W_N2).
    // Reference: Cantera 3.1.0 (same h2o2.yaml transport parameters).
    // -----------------------------------------------------------------------
    #[test]
    fn test_wilke_h2n2_300k() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let mut x = vec![0.0_f64; nk];
        x[mech.species_index("H2").unwrap()] = 0.5;
        x[mech.species_index("N2").unwrap()] = 0.5;
        // Reference: analytical Wilke formula using μ_H2=9.78e-6, μ_N2=1.809e-5 (both
        // validated against Cantera 3.1.0 at 0.5% rtol). The buggy formula
        // ((Wk/Wj)^0.125 instead of (Wj/Wk)^0.25) gives ~1.50e-5 here, a 13%
        // error that this test catches. 2% tolerance is tighter than the bug.
        check("mu_H2N2 300K", mixture_viscosity(&mech, &x, 300.0), 1.727e-5, 2e-2);
    }
}

/// Convert mass fractions to mole fractions.
pub fn mass_to_mole_fractions(mech: &Mechanism, y: &[f64]) -> Vec<f64> {
    let w_mean = mean_molecular_weight(&mech.species, y);
    mech.species.iter().zip(y.iter())
        .map(|(s, &yk)| yk * w_mean / s.molecular_weight)
        .collect()
}
