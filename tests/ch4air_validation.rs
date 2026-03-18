/// End-to-end CH4/air flame validation (issue #15).
///
/// Runs the full solver starting from the Cantera reference profile.
/// Verifies:
///   - The laminar flame speed Su is within 5% of the Cantera reference
///     value of 0.3814 m/s (mixture-averaged, gri30.yaml, φ=1.0,
///     T_u=300 K, P=1 atm).
///   - T_max is between 2000 K and 2400 K.

use tempfile::NamedTempFile;

use premix1d::{
    io::input::FlameConfig,
    flame::solver_driver::run_flame,
};

const CANTERA_SU_REFERENCE: f64 = 0.3814; // m/s, mixture-averaged, gri30.yaml, φ=1 CH4/air

fn manifest() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

#[test]
fn test_ch4air_flame_speed_within_5pct_of_cantera() {
    let manifest = manifest();

    let out_file = NamedTempFile::new().expect("tmp output file");
    let out_path = out_file.path().to_str().unwrap().to_string();

    let toml_content = format!(
        r#"
[mechanism]
file = "{manifest}/data/gri30.yaml"
format = "cantera_yaml"

[flame]
pressure = 101325.0
fuel = {{ CH4 = 1.0 }}
oxidizer = {{ O2 = 1.0, N2 = 3.76 }}
equivalence_ratio = 1.0
t_unburned = 300.0
domain_length = 0.03

[grid]
initial_points = 100
max_points = 1000
grad = 0.07
curv = 0.14

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
time_steps = 500
dt_initial = 1.0e-7
su_initial_guess = 0.3814
initial_profile = "{manifest}/data/cantera_ch4air_initial.csv"

[output]
file = "{out_path}"
"#
    );

    let config: FlameConfig = toml::from_str(&toml_content)
        .expect("parse config");

    let t_start = std::time::Instant::now();
    run_flame(&config).expect("run_flame must not fail");
    let elapsed = t_start.elapsed();
    println!("CH4/air solver wall-clock time: {:.1}s", elapsed.as_secs_f64());

    let csv_content = std::fs::read_to_string(&out_path)
        .expect("read output CSV");
    let mut lines = csv_content.lines();
    let header = lines.next().expect("header");
    let headers: Vec<&str> = header.split(',').collect();
    let u_col = headers.iter().position(|h| h.trim() == "u [m/s]").expect("u [m/s] column");
    let first_row = lines.next().expect("first data row");
    let cols: Vec<f64> = first_row.split(',')
        .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
        .collect();
    let su = cols[u_col];

    let rel_err = (su - CANTERA_SU_REFERENCE).abs() / CANTERA_SU_REFERENCE;

    let t_max = csv_content.lines().skip(1)
        .filter_map(|l| l.split(',').nth(1)?.trim().parse::<f64>().ok())
        .fold(0.0_f64, f64::max);

    println!(
        "CH4/air validation: Su = {su:.4} m/s  (Cantera ref = {CANTERA_SU_REFERENCE:.4} m/s, \
         err = {:.1}%,  T_max = {t_max:.0} K)",
        rel_err * 100.0
    );

    assert!(
        (2000.0..2400.0).contains(&t_max),
        "T_max = {t_max:.1} K, expected 2000–2400 K"
    );

    assert!(
        rel_err < 0.05,
        "Su = {su:.4} m/s  (Cantera reference = {CANTERA_SU_REFERENCE:.4} m/s, \
         relative error = {:.2}%, limit = 5%)",
        rel_err * 100.0
    );
}
