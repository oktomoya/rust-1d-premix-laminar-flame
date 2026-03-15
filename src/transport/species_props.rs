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

/// Species thermal conductivity [W/(m·K)] via the full Mason-Monchick (1962) formula.
///
/// Matches Cantera's `GasTransport::fitProperties` exactly, including the
/// rotational relaxation correction using `Zrot` (rotational collision number).
///
/// cv_rot (in units of R):  0.0 (atom) | 1.0 (linear) | 1.5 (nonlinear)
/// cv_int (vibrational, in units of R):  cp/R - 2.5 - cv_rot
///
/// Rotational relaxation factor:
///   fz(T*) = 1 + π^1.5/√T* · (0.5 + 1/T*) + (π²/4 + 2)/T*
///   A  = 2.5 - f_int
///   B  = Zrot · fz(298/ε) / fz(T/ε)  +  (2/π) · (5/3 · cv_rot + f_int)
///   c1 = (2/π) · A / B
///   f_rot   = f_int · (1 + c1)
///   f_trans = 2.5 · (1 − c1 · cv_rot / 1.5)
///   λk = (μk/Wk) · R · (f_trans·1.5 + f_rot·cv_rot + f_int·cv_int)
pub fn thermal_conductivity(species: &Species, mu_k: f64, cp_k: f64, t: f64, pressure: f64) -> f64 {
    use crate::chemistry::species::GeometryType;
    use crate::chemistry::thermo::R_UNIVERSAL;
    use std::f64::consts::PI;

    let r_over_w = R_UNIVERSAL / species.molecular_weight;

    // cv_rot in units of R: 0 (atom), 1 (linear), 3/2 (nonlinear)
    let cv_rot: f64 = match species.transport.geometry {
        GeometryType::Atom     => 0.0,
        GeometryType::Linear   => 1.0,
        GeometryType::Nonlinear => 1.5,
    };

    // cp/R (dimensionless); vibrational internal modes (also in units of R)
    let cp_r   = cp_k / r_over_w;
    let cv_int = (cp_r - 2.5 - cv_rot).max(0.0);

    // f_int = ρk * D_kk / μk,  where ρk = P·W/(R·T) for pure species k
    let d_kk  = binary_diffusion(species, species, t, pressure);
    let rho_k = pressure * species.molecular_weight / (R_UNIVERSAL * t);
    let f_int = rho_k * d_kk / mu_k;

    // fz(T*) — temperature-dependent rotational relaxation factor
    let fz = |ts: f64| -> f64 {
        let ts = ts.max(0.01);
        1.0 + PI.powf(1.5) / ts.sqrt() * (0.5 + 1.0 / ts)
            + (0.25 * PI * PI + 2.0) / ts
    };
    let eps_k     = species.transport.well_depth.max(1.0); // ε/kb [K]
    let t_star    = t / eps_k;
    let t_star_298 = 298.0 / eps_k;
    let zrot      = species.transport.rot_relax;

    let a_factor = 2.5 - f_int;
    let b_factor = zrot * fz(t_star_298) / fz(t_star)
        + (2.0 / PI) * (5.0 / 3.0 * cv_rot + f_int);
    let c1 = if b_factor.abs() > 1e-30 {
        (2.0 / PI) * a_factor / b_factor
    } else {
        0.0
    };

    let f_rot   = f_int * (1.0 + c1);
    let f_trans = 2.5 * (1.0 - c1 * cv_rot / 1.5);

    // λk = (μk/Wk) · R · (f_trans·1.5 + f_rot·cv_rot + f_int·cv_int)
    mu_k * r_over_w * (f_trans * 1.5 + f_rot * cv_rot + f_int * cv_int)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;

    fn h2o2_mech() -> crate::chemistry::mechanism::Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    fn check(label: &str, got: f64, expected: f64, rtol: f64) {
        let rel = (got - expected).abs() / expected.abs();
        assert!(rel < rtol, "{label}: got {got:.6e}, expected {expected:.6e}, rel err {rel:.2e}");
    }

    // -----------------------------------------------------------------------
    // Viscosity — Chapman-Enskog, Neufeld Ω*(2,2)
    // Reference: Cantera 3.1.0 (same h2o2.yaml transport parameters).
    // Residual ~0.06% from Cantera's use of a bilinear-interpolation table
    // vs our Neufeld polynomial; 0.5% tolerance is ample.
    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // Thermal conductivity — Mason-Monchick formula
    //
    // References for sanity-check comparisons:
    //   NIST dilute-gas values (ideal-gas limit):
    //     N2 300 K: 25.93 mW/(m·K)  H2 300 K: 186.7 mW/(m·K)  Ar 300 K: 17.72 mW/(m·K)
    //   Chapman-Enskog + Mason-Monchick typically overestimates λ by 2–8% vs NIST for
    //   simple non-polar species; tolerance set to 10% to accommodate this.
    //
    //   Note: reference values below should be cross-validated against Cantera 3.1.0
    //   once available (see scripts/verify_transport.py).
    // -----------------------------------------------------------------------
    #[test]
    fn test_thermal_conductivity_n2_300k() {
        let mech = h2o2_mech();
        use crate::chemistry::thermo::cp_species;
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        let mu = viscosity(sp, 300.0);
        let cp = cp_species(sp, 300.0);
        // Reference: Cantera 3.1.0 pure-N2, mixture-averaged, h2o2.yaml at 300 K.
        check("lambda_N2 300K", thermal_conductivity(sp, mu, cp, 300.0, 101325.0), 2.646311e-2, 1e-2);
    }

    #[test]
    fn test_thermal_conductivity_h2_300k() {
        let mech = h2o2_mech();
        use crate::chemistry::thermo::cp_species;
        let sp = &mech.species[mech.species_index("H2").unwrap()];
        let mu = viscosity(sp, 300.0);
        let cp = cp_species(sp, 300.0);
        // Reference: Cantera 3.1.0 pure-H2, mixture-averaged, h2o2.yaml at 300 K.
        check("lambda_H2 300K", thermal_conductivity(sp, mu, cp, 300.0, 101325.0), 1.869231e-1, 1e-2);
    }

    #[test]
    fn test_thermal_conductivity_ar_300k() {
        let mech = h2o2_mech();
        use crate::chemistry::thermo::cp_species;
        let sp = &mech.species[mech.species_index("AR").unwrap()];
        let mu = viscosity(sp, 300.0);
        let cp = cp_species(sp, 300.0);
        // Reference: Cantera 3.1.0 pure-Ar, mixture-averaged, h2o2.yaml at 300 K.
        check("lambda_AR 300K", thermal_conductivity(sp, mu, cp, 300.0, 101325.0), 1.805973e-2, 1e-2);
    }

    #[test]
    fn test_viscosity_n2_300k() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        check("mu_N2 300K", viscosity(sp, 300.0), 1.808570419229e-5, 5e-3);
    }

    #[test]
    fn test_viscosity_n2_1000k() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        check("mu_N2 1000K", viscosity(sp, 1000.0), 4.149871864818e-5, 5e-3);
    }

    // -----------------------------------------------------------------------
    // Binary diffusion — Chapman-Enskog BSL formula, Neufeld Ω*(1,1)
    // D_H2N2 in the dilute limit (X_H2 → 0) compared against
    // Cantera mix_diff_coeffs with X_H2=0.001, X_N2=0.999.
    // -----------------------------------------------------------------------
    #[test]
    fn test_binary_diffusion_h2n2_300k() {
        let mech = h2o2_mech();
        let h2 = &mech.species[mech.species_index("H2").unwrap()];
        let n2 = &mech.species[mech.species_index("N2").unwrap()];
        check("D_H2N2 300K", binary_diffusion(h2, n2, 300.0, 101325.0), 7.796992738673e-5, 5e-3);
    }

    #[test]
    fn test_binary_diffusion_h2n2_1000k() {
        let mech = h2o2_mech();
        let h2 = &mech.species[mech.species_index("H2").unwrap()];
        let n2 = &mech.species[mech.species_index("N2").unwrap()];
        check("D_H2N2 1000K", binary_diffusion(h2, n2, 1000.0, 101325.0), 5.856224715919e-4, 5e-3);
    }
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
    // BSL formula: Dij [cm²/s] = 1.858e-3 * T^1.5 * sqrt(1/Wi + 1/Wj) / (P_atm * σij² * Ω*(1,1))
    // Equivalently: 1.858e-3 * sqrt(2/Wij) = 2.6280e-3 / sqrt(Wij)
    let d_cm2_s = 2.6280e-3 * t.powf(1.5) / (p_atm * sigma_ij.powi(2) * om11 * w_ij.sqrt());
    d_cm2_s * 1e-4 // cm²/s → m²/s
}
