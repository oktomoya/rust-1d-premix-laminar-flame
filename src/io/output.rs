use anyhow::Result;
use crate::chemistry::kinetics::production_rates;
use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{density, enthalpy_molar, mean_molecular_weight};
use crate::flame::domain::Grid;
use crate::flame::state::FlameState;

/// Write flame solution to a CSV file.
/// Columns: z [m], T [K], u [m/s], rho [kg/m³], hrr [W/m³],
///          X_species1, …, Y_species1, …
pub fn write_csv(
    path: &str,
    x: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    pressure: f64,
) -> Result<()> {
    let mut wtr = csv::Writer::from_path(path)?;
    let nk = mech.n_species();
    let nj = grid.n_points();

    // Header
    let mut header = vec![
        "z_m".to_string(), "T_K".to_string(), "u_m_s".to_string(),
        "rho_kg_m3".to_string(), "hrr_W_m3".to_string(),
    ];
    for sp in &mech.species { header.push(format!("X_{}", sp.name)); }
    for sp in &mech.species { header.push(format!("Y_{}", sp.name)); }
    wtr.write_record(&header)?;

    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();

    for j in 0..nj {
        let t = state.temperature(j);
        let y: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
        let w = mean_molecular_weight(&mech.species, &y);
        let rho = density(pressure, t, w);
        let u = m / rho.max(1e-30);

        // Mole fractions: Xk = (Yk/Wk) / Σj(Yj/Wj)
        let sum_y_over_w: f64 = (0..nk)
            .map(|k| y[k] / mech.species[k].molecular_weight)
            .sum();
        let xk: Vec<f64> = (0..nk)
            .map(|k| y[k] / mech.species[k].molecular_weight / sum_y_over_w.max(1e-300))
            .collect();

        // Heat release rate: HRR = -Σk ωk·hk  [W/m³]
        // concentrations in mol/cm³ (CGS, matching stored A values)
        let conc: Vec<f64> = (0..nk)
            .map(|k| rho * y[k] / mech.species[k].molecular_weight * 1e-6)
            .collect();
        let wdot = production_rates(mech, t, &conc, pressure); // mol/(cm³·s)
        let hrr = -(0..nk)
            .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t))
            .sum::<f64>() * 1e6; // mol/(cm³·s)×J/mol×1e6 → W/m³

        let mut row = vec![
            format!("{:.6e}", grid.z[j]),
            format!("{:.4}", t),
            format!("{:.6e}", u),
            format!("{:.6e}", rho),
            format!("{:.6e}", hrr),
        ];
        for k in 0..nk { row.push(format!("{:.6e}", xk[k])); }
        for k in 0..nk { row.push(format!("{:.6e}", y[k])); }
        wtr.write_record(&row)?;
    }

    wtr.flush()?;
    println!("Solution written to {path}");
    Ok(())
}

/// Print a summary of the flame solution.
pub fn print_summary(x: &[f64], mech: &Mechanism, grid: &Grid, pressure: f64) {
    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();
    let nk = mech.n_species();

    // Flame speed = M / ρ_unburned
    let j = 0;
    let t0 = state.temperature(j);
    let y0: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
    let w0 = mean_molecular_weight(&mech.species, &y0);
    let rho0 = density(pressure, t0, w0);
    let su = m / rho0.max(1e-30);

    // Maximum temperature
    let t_max = (0..grid.n_points())
        .map(|j| state.temperature(j))
        .fold(f64::NEG_INFINITY, f64::max);

    println!("---------------------------------------------------");
    println!("  Laminar flame speed Su = {su:.4} m/s");
    println!("  Max temperature       = {t_max:.1} K");
    println!("  Grid points           = {}", grid.n_points());
    println!("  Mass flux M           = {m:.4e} kg/(m²·s)");
    println!("---------------------------------------------------");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;
    use crate::flame::domain::Grid;
    use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};

    fn h2o2_mech() -> Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    // Build a uniform N2-only solution and verify mole fractions == mass fractions
    // (since X = Y when only one species is present), and HRR ≈ 0 (no reactions).
    #[test]
    fn test_mole_fractions_pure_n2() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nv = natj(&mech);
        let nj = 4;
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();

        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            x[idx_t(nv, j)] = 1000.0;
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = 0.1;

        let state = FlameState::new(&x, &mech, &grid);
        let pressure = 101325.0;

        for j in 0..nj {
            let y: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
            let sum_y_over_w: f64 = (0..nk)
                .map(|k| y[k] / mech.species[k].molecular_weight)
                .sum();
            let xk: Vec<f64> = (0..nk)
                .map(|k| y[k] / mech.species[k].molecular_weight / sum_y_over_w.max(1e-300))
                .collect();

            // Pure N2: X_N2 should be 1.0
            assert!((xk[n2_idx] - 1.0).abs() < 1e-10,
                "X_N2[{j}] = {:.6}, expected 1.0", xk[n2_idx]);

            // Mole fractions sum to 1
            let sum_x: f64 = xk.iter().sum();
            assert!((sum_x - 1.0).abs() < 1e-10,
                "Σ Xk[{j}] = {sum_x:.6}, expected 1.0");

            // HRR ≈ 0 for pure N2 (no reactions involve only N2)
            let t = state.temperature(j);
            let w = mean_molecular_weight(&mech.species, &y);
            let rho = density(pressure, t, w);
            let conc: Vec<f64> = (0..nk)
                .map(|k| rho * y[k] / mech.species[k].molecular_weight * 1e-6)
                .collect();
            let wdot = production_rates(&mech, t, &conc, pressure);
            let hrr = -(0..nk)
                .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t))
                .sum::<f64>() * 1e6;
            assert!(hrr.abs() < 1.0, "HRR[{j}] = {hrr:.3e} W/m³, expected ≈ 0");
        }
    }

    // Verify write_csv produces a file with the expected number of columns.
    #[test]
    fn test_write_csv_column_count() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nv = natj(&mech);
        let nj = 3;
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();

        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            x[idx_t(nv, j)] = 1000.0;
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = 0.1;

        let path = format!("/tmp/test_output_{}.csv", std::process::id());
        write_csv(&path, &x, &mech, &grid, 101325.0).expect("write_csv must succeed");

        let mut rdr = csv::Reader::from_path(&path).expect("read csv");
        let headers = rdr.headers().expect("headers").clone();
        // z, T, u, rho, hrr, X_k×nk, Y_k×nk
        let expected_cols = 5 + 2 * nk;
        assert_eq!(headers.len(), expected_cols,
            "expected {expected_cols} columns, got {}: {:?}", headers.len(), headers);

        let _ = std::fs::remove_file(&path);
    }
}
