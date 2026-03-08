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
            let w_ratio = (mech.species[k].molecular_weight / mech.species[j].molecular_weight).sqrt();
            let phi_kj = (1.0 + (mu[k] / mu[j]).sqrt() * w_ratio.powf(0.25)).powi(2)
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
) -> f64 {
    let _ = mass_fractions;
    let nk = mech.n_species();
    let mut sum1 = 0.0_f64;
    let mut sum2 = 0.0_f64;
    for k in 0..nk {
        let mu_k = viscosity(&mech.species[k], t);
        let cp_k = cp_species(&mech.species[k], t);
        let lam_k = thermal_conductivity(&mech.species[k], mu_k, cp_k, t);
        sum1 += mole_fractions[k] * lam_k;
        sum2 += mole_fractions[k] / lam_k.max(1e-30);
    }
    0.5 * (sum1 + 1.0 / sum2)
}

/// Mixture-averaged diffusion coefficients Dkm [m²/s] for each species k.
/// Dkm = (1 - Xk) / Σ_{j≠k} (Xj / Dkj)
pub fn mixture_diffusion_coefficients(
    mech: &Mechanism,
    mole_fractions: &[f64],
    t: f64,
    pressure: f64,
) -> Vec<f64> {
    let nk = mech.n_species();
    let mut dkm = vec![0.0_f64; nk];

    for k in 0..nk {
        let xk = mole_fractions[k];
        let sum: f64 = (0..nk)
            .filter(|&j| j != k)
            .map(|j| {
                let dkj = binary_diffusion(&mech.species[k], &mech.species[j], t, pressure);
                mole_fractions[j] / dkj.max(1e-30)
            })
            .sum();
        dkm[k] = if sum > 1e-30 {
            (1.0 - xk).max(0.0) / sum
        } else {
            0.0
        };
    }
    dkm
}

/// Convert mass fractions to mole fractions.
pub fn mass_to_mole_fractions(mech: &Mechanism, y: &[f64]) -> Vec<f64> {
    let w_mean = mean_molecular_weight(&mech.species, y);
    mech.species.iter().zip(y.iter())
        .map(|(s, &yk)| yk * w_mean / s.molecular_weight)
        .collect()
}
