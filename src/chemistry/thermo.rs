/// Universal gas constant [J/(mol·K)]
pub const R_UNIVERSAL: f64 = 8.314462618;

use crate::chemistry::species::Species;

/// Evaluate cp/R for a species at temperature t using NASA polynomials.
pub fn cp_over_r(species: &Species, t: f64) -> f64 {
    let a = &species.nasa_poly(t).coeffs;
    a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])))
}

/// Specific heat at constant pressure [J/(kg·K)] for species k.
pub fn cp_species(species: &Species, t: f64) -> f64 {
    cp_over_r(species, t) * R_UNIVERSAL / species.molecular_weight
}

/// Specific enthalpy [J/kg] for species k at temperature t.
pub fn enthalpy_species(species: &Species, t: f64) -> f64 {
    let a = &species.nasa_poly(t).coeffs;
    let h_over_rt = a[0] + t * (a[1] / 2.0
        + t * (a[2] / 3.0
        + t * (a[3] / 4.0
        + t * a[4] / 5.0)))
        + a[5] / t;
    h_over_rt * R_UNIVERSAL * t / species.molecular_weight
}

/// Molar enthalpy [J/mol] for species k at temperature t.
pub fn enthalpy_molar(species: &Species, t: f64) -> f64 {
    let a = &species.nasa_poly(t).coeffs;
    let h_over_rt = a[0] + t * (a[1] / 2.0
        + t * (a[2] / 3.0
        + t * (a[3] / 4.0
        + t * a[4] / 5.0)))
        + a[5] / t;
    h_over_rt * R_UNIVERSAL * t
}

/// Standard entropy [J/(kg·K)] for species k at temperature t (1 atm reference).
pub fn entropy_species(species: &Species, t: f64) -> f64 {
    let a = &species.nasa_poly(t).coeffs;
    let s_over_r = a[0] * t.ln() + t * (a[1] + t * (a[2] / 2.0 + t * (a[3] / 3.0 + t * a[4] / 4.0))) + a[6];
    s_over_r * R_UNIVERSAL / species.molecular_weight
}

/// Mixture-averaged cp [J/(kg·K)] given mass fractions y[k].
pub fn cp_mixture(species: &[Species], y: &[f64], t: f64) -> f64 {
    species.iter().zip(y.iter()).map(|(s, &yk)| yk * cp_species(s, t)).sum()
}

/// Mixture mean molecular weight [kg/mol] given mass fractions y[k].
pub fn mean_molecular_weight(species: &[Species], y: &[f64]) -> f64 {
    let sum_y_over_w: f64 = species.iter().zip(y.iter())
        .map(|(s, &yk)| yk / s.molecular_weight)
        .sum();
    1.0 / sum_y_over_w
}

/// Mixture density [kg/m³] via ideal gas law.
pub fn density(pressure: f64, t: f64, mean_mw: f64) -> f64 {
    pressure * mean_mw / (R_UNIVERSAL * t)
}
