use anyhow::Result;
use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{density, mean_molecular_weight};
use crate::flame::domain::Grid;
use crate::flame::state::{idx_t, idx_y, natj, FlameState};

/// Write flame solution to a CSV file.
/// Columns: z [m], T [K], u [m/s], rho [kg/m³], Y_species1, Y_species2, ...
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
    let nv = natj(mech);

    // Header
    let mut header = vec!["z_m".to_string(), "T_K".to_string(), "u_m_s".to_string(), "rho_kg_m3".to_string()];
    for sp in &mech.species {
        header.push(format!("Y_{}", sp.name));
    }
    wtr.write_record(&header)?;

    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();

    for j in 0..nj {
        let t = state.temperature(j);
        let y: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
        let w = mean_molecular_weight(&mech.species, &y);
        let rho = density(pressure, t, w);
        let u = m / rho.max(1e-30);

        let mut row = vec![
            format!("{:.6e}", grid.z[j]),
            format!("{:.4}", t),
            format!("{:.6e}", u),
            format!("{:.6e}", rho),
        ];
        for k in 0..nk {
            row.push(format!("{:.6e}", y[k]));
        }
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
