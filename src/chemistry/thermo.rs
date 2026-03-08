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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;

    fn h2o2_mech() -> crate::chemistry::mechanism::Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    /// Relative tolerance for comparison against Cantera reference values.
    /// Both sides evaluate the same NASA7 polynomial with the same coefficients
    /// and the same T_mid boundary condition (T <= T_mid → low polynomial).
    /// Residual discrepancy is pure f64 rounding (~1e-14); 1e-8 is ample.
    const RTOL: f64 = 1e-8;

    fn check(label: &str, got: f64, expected: f64) {
        let rel = (got - expected).abs() / expected.abs().max(1e-10);
        assert!(
            rel < RTOL,
            "{label}: got {got:.6e}, expected {expected:.6e}, rel error {rel:.2e}"
        );
    }

    // -----------------------------------------------------------------------
    // Reference values from Cantera 3.1.0 using the same h2o2.yaml NASA7
    // polynomials.  All per-mass quantities (cp, h, s) are in J/kg or
    // J/(kg·K).  Enthalpy reference is the standard-state formation enthalpy
    // embedded in the NASA7 a[5] coefficient (298.15 K, 1 atm).
    // -----------------------------------------------------------------------

    #[test]
    fn test_cp_h2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2").unwrap()];
        check("cp H2 300K",  cp_species(sp,  300.0), 14310.9052553698);
        check("cp H2 1000K", cp_species(sp, 1000.0), 14961.8782033635);
        check("cp H2 2000K", cp_species(sp, 2000.0), 16992.5956055308);
    }

    #[test]
    fn test_cp_o2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("O2").unwrap()];
        check("cp O2 300K",  cp_species(sp,  300.0),  918.4346250542);
        check("cp O2 1000K", cp_species(sp, 1000.0), 1090.1610871509);
        check("cp O2 2000K", cp_species(sp, 2000.0), 1181.2113707811);
    }

    #[test]
    fn test_cp_h2o() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2O").unwrap()];
        check("cp H2O 300K",  cp_species(sp,  300.0), 1864.9154285171);
        check("cp H2O 1000K", cp_species(sp, 1000.0), 2292.2422463756);
        check("cp H2O 2000K", cp_species(sp, 2000.0), 2872.7122283858);
    }

    #[test]
    fn test_cp_n2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        check("cp N2 300K",  cp_species(sp,  300.0), 1037.8911357957);
        check("cp N2 1000K", cp_species(sp, 1000.0), 1169.4847572643);
        check("cp N2 2000K", cp_species(sp, 2000.0), 1284.6545256479);
    }

    #[test]
    fn test_enthalpy_h2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2").unwrap()];
        check("h H2 300K",  enthalpy_species(sp,  300.0),    26468.504563);
        check("h H2 1000K", enthalpy_species(sp, 1000.0), 10261177.530487);
        check("h H2 2000K", enthalpy_species(sp, 2000.0), 26260574.993899);
    }

    #[test]
    fn test_enthalpy_o2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("O2").unwrap()];
        check("h O2 300K",  enthalpy_species(sp,  300.0),    1698.818008);
        check("h O2 1000K", enthalpy_species(sp, 1000.0),  709632.193256);
        check("h O2 2000K", enthalpy_species(sp, 2000.0), 1850273.617578);
    }

    #[test]
    fn test_enthalpy_h2o() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2O").unwrap()];
        check("h H2O 300K",  enthalpy_species(sp,  300.0), -13420065.305308);
        check("h H2O 1000K", enthalpy_species(sp, 1000.0), -11980133.500956);
        check("h H2O 2000K", enthalpy_species(sp, 2000.0),  -9369299.575998);
    }

    #[test]
    fn test_enthalpy_n2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        check("h N2 300K",  enthalpy_species(sp,  300.0),    1970.993858);
        check("h N2 1000K", enthalpy_species(sp, 1000.0),  766397.701133);
        check("h N2 2000K", enthalpy_species(sp, 2000.0), 2003721.968647);
    }

    #[test]
    fn test_entropy_h2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2").unwrap()];
        check("s H2 300K",  entropy_species(sp,  300.0), 64910.0638546894);
        check("s H2 1000K", entropy_species(sp, 1000.0), 82458.1945217282);
        check("s H2 2000K", entropy_species(sp, 2000.0), 93466.2619901014);
    }

    #[test]
    fn test_entropy_o2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("O2").unwrap()];
        check("s O2 300K",  entropy_species(sp,  300.0), 6416.9652759636);
        check("s O2 1000K", entropy_species(sp, 1000.0), 7612.5505786534);
        check("s O2 2000K", entropy_species(sp, 2000.0), 8399.5934864246);
    }

    #[test]
    fn test_entropy_h2o() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("H2O").unwrap()];
        check("s H2O 300K",  entropy_species(sp,  300.0), 10493.2462571728);
        check("s H2O 1000K", entropy_species(sp, 1000.0), 12918.9567443583);
        check("s H2O 2000K", entropy_species(sp, 2000.0), 14705.2885288374);
    }

    #[test]
    fn test_entropy_n2() {
        let mech = h2o2_mech();
        let sp = &mech.species[mech.species_index("N2").unwrap()];
        check("s N2 300K",  entropy_species(sp,  300.0), 6842.7243797601);
        check("s N2 1000K", entropy_species(sp, 1000.0), 8141.9484588735);
        check("s N2 2000K", entropy_species(sp, 2000.0), 8994.9805533253);
    }
}
