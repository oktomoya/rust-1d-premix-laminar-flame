/// End-to-end H2/air flame validation (issue #14).
///
/// Runs the full solver (PT-skipped, multi-pass Newton + adaptive grid
/// refinement) starting from the Cantera reference profile.  Verifies:
///   - Newton converges in each pass (no divergence)
///   - The laminar flame speed Su is within 0.5% of the Cantera reference
///     value of 2.3354 m/s (mixture-averaged transport, mole-fraction basis,
///     h2o2.yaml, φ=1.0, T_u=300 K, P=1 atm).
///   - T_max is between 2000 K and 2600 K.
///
/// Remaining error (~0.3%): grid resolution and upwind convection scheme
/// differences vs Cantera's adaptive solver.

use std::io::Write;
use tempfile::NamedTempFile;

use premix1d::{
    io::input::FlameConfig,
    flame::solver_driver::run_flame,
};

const CANTERA_SU_REFERENCE: f64 = 2.3354; // m/s, mixture-averaged, h2o2.yaml, φ=1 H2/air

fn manifest() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

#[test]
fn test_h2air_flame_speed_within_5pct_of_cantera() {
    let manifest = manifest();

    // Write a temporary output CSV (avoid cluttering the workspace).
    let out_file = NamedTempFile::new().expect("tmp output file");
    let out_path = out_file.path().to_str().unwrap().to_string();

    // Build TOML config that matches hydrogen_air.toml but uses the Cantera
    // initial profile and writes output to a temp file.
    let toml_content = format!(
        r#"
[mechanism]
file = "{manifest}/data/h2o2.yaml"
format = "cantera_yaml"

[flame]
pressure = 101325.0
fuel = {{ H2 = 1.0 }}
oxidizer = {{ O2 = 0.21, N2 = 0.79 }}
equivalence_ratio = 1.0
t_unburned = 300.0
domain_length = 0.02

[grid]
initial_points = 100
max_points = 1000
grad = 0.05
curv = 0.10

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
time_steps = 500
dt_initial = 1.0e-7
su_initial_guess = 2.3354
initial_profile = "{manifest}/data/cantera_h2air_initial.csv"

[output]
file = "{out_path}"
"#
    );

    let config: FlameConfig = toml::from_str(&toml_content)
        .expect("parse config");

    run_flame(&config).expect("run_flame must not fail");

    // Read back the CSV output to extract Su.
    let csv_content = std::fs::read_to_string(&out_path)
        .expect("read output CSV");
    let mut lines = csv_content.lines();
    let _header = lines.next().expect("header");
    // z [m], T [K], u [m/s], rho [kg/m3], hrr [W/m3], X_*, Y_*
    // We need u [m/s] at z=0 (left boundary) which equals Su.
    let first_row = lines.next().expect("first data row");
    let cols: Vec<f64> = first_row.split(',')
        .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
        .collect();
    let su = cols[2]; // u_m_s column

    let rel_err = (su - CANTERA_SU_REFERENCE).abs() / CANTERA_SU_REFERENCE;

    // Also extract T_max from all rows.
    let t_max = csv_content.lines().skip(1)
        .filter_map(|l| l.split(',').nth(1)?.trim().parse::<f64>().ok())
        .fold(0.0_f64, f64::max);

    println!(
        "H2/air validation: Su = {su:.4} m/s  (Cantera ref = {CANTERA_SU_REFERENCE:.4} m/s, \
         err = {:.1}%,  T_max = {t_max:.0} K)",
        rel_err * 100.0
    );

    assert!(
        (2000.0..2600.0).contains(&t_max),
        "T_max = {t_max:.1} K, expected 2000–2600 K"
    );

    assert!(
        rel_err < 0.005,
        "Su = {su:.4} m/s  (Cantera reference = {CANTERA_SU_REFERENCE:.4} m/s, \
         relative error = {:.2}%, limit = 0.5%)",
        rel_err * 100.0
    );
}
