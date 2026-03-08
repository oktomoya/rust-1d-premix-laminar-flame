/// Species-level transport properties via Chapman-Enskog kinetic theory.

use crate::chemistry::species::Species;
use crate::transport::collision_integrals::{omega11, omega22};

/// Species viscosity [Pa·s] via Chapman-Enskog theory.
/// μk = (5/16) * sqrt(π·Wk·kb·T) / (π·σk²·Ω*(2,2))
/// Simplified (pre-factors absorbed into the constant 2.6693e-6):
///   μk [Pa·s] = 2.6693e-6 * sqrt(Wk[g/mol] * T) / (σk[Å]² * Ω*(2,2)(T*))
pub fn viscosity(species: &Species, t: f64) -> f64 {
    let t_star = t / species.transport.well_depth;
    let om22 = omega22(t_star.max(0.1));
    let w_g_mol = species.molecular_weight * 1000.0; // kg/mol → g/mol
    2.6693e-6 * (w_g_mol * t).sqrt() / (species.transport.diameter.powi(2) * om22)
    // [Pa·s] — the 2.6693e-6 factor gives SI directly when σ is in Angstrom
}

/// Species thermal conductivity [W/(m·K)] via modified Eucken approximation.
///   λk = μk * cp_trans / Wk    (translational + internal contributions)
/// Using Mason-Monchick: λk = μk/Wk * (R * (f_trans * 5/2 + f_int * cv_int))
/// Simple Eucken: λk = μk * (cp_k + 1.25 * R/Wk)
pub fn thermal_conductivity(species: &Species, mu_k: f64, cp_k: f64, t: f64) -> f64 {
    use crate::chemistry::thermo::R_UNIVERSAL;
    let _ = t;
    // Modified Eucken formula
    mu_k * (cp_k + 1.25 * R_UNIVERSAL / species.molecular_weight)
}

/// Binary diffusion coefficient Dij [m²/s] between species i and j.
/// Dij = 3/(16·n) * sqrt(2·π·kbT / (μij)) / (π·σij²·Ω*(1,1))
/// Simplified: Dij = 2.6280e-5 * T^1.5 / (P * σij² * Ω*(1,1) * sqrt(Wij))
/// where σij = (σi + σj)/2, Wij = 2·Wi·Wj/(Wi+Wj), P in Pa.
pub fn binary_diffusion(sp_i: &Species, sp_j: &Species, t: f64, pressure: f64) -> f64 {
    let sigma_ij = 0.5 * (sp_i.transport.diameter + sp_j.transport.diameter); // Angstrom
    let eps_ij   = (sp_i.transport.well_depth * sp_j.transport.well_depth).sqrt(); // K
    let wi_g = sp_i.molecular_weight * 1000.0; // g/mol
    let wj_g = sp_j.molecular_weight * 1000.0;
    let w_ij = 2.0 * wi_g * wj_g / (wi_g + wj_g); // reduced molecular weight [g/mol]

    let t_star = t / eps_ij.max(1.0);
    let om11 = omega11(t_star.max(0.1));

    // Pressure in atm for the standard formula
    let p_atm = pressure / 101325.0;
    // Dij [cm²/s] = 1.858e-3 * T^1.5 / (P_atm * σij² * Ω*(1,1) * sqrt(Wij))
    let d_cm2_s = 1.858e-3 * t.powf(1.5) / (p_atm * sigma_ij.powi(2) * om11 * w_ij.sqrt());
    d_cm2_s * 1e-4 // cm²/s → m²/s
}
