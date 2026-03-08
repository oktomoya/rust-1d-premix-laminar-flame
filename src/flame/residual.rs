/// Residual evaluation F(x) for the 1D freely propagating premixed flame.
///
/// Governing equations (SI units throughout):
///
/// Interior points j = 1 … J-2:
///   Species k (k = 0 … K-2):
///     F_yk = M * (Yk_j - Yk_{j-1})/dz_m  +  (jk_{j-1/2} - jk_{j-3/2}) / dz_av  -  ωk * Wk
///   Species K-1: sum closure  F_yK = ΣYk - 1
///   Energy:
///     F_T = -M * cp * (T_j - T_{j-1})/dz_m
///           + (λ_{j+1/2}*(T_{j+1}-T_j)/dz_p - λ_{j-1/2}*(T_j-T_{j-1})/dz_m) / dz_av
///           - Σk (jk_{j-1/2} + jk_{j+1/2})/2 * cpk * (T_j - T_{j-1})/dz_m  (enthalpy transport)
///           - Σk ωk * hk                                       (heat release [W/m³])
///
/// Left boundary (j = 0):
///   F_yk = M * (Yk_0 - Yku_k) + jk_{1/2}    (inlet flux condition)
///   F_T  = T_0 - T_unburned
///
/// Right boundary (j = J-1):
///   F_yk = Yk_{J-1} - Yk_{J-2}   (zero gradient)
///   F_T  = T_{J-1} - T_{J-2}
///
/// Eigenvalue closure: at the fixed point j = j_fix:
///   F_M = T_{j_fix} - T_fix   (replaces the trivial M equation at that point)
///
/// Pseudo-transient embedding (optional):
///   When x_old and rdt = 1/dt are provided, adds rdt * (x - x_old) to all
///   interior equations, implementing backward-Euler time-stepping.

use crate::chemistry::kinetics::production_rates;
use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{
    cp_mixture, cp_species, density, enthalpy_molar, mean_molecular_weight,
};
use crate::flame::domain::Grid;
use crate::flame::state::{idx_m, idx_t, idx_y, natj, FlameState};
use crate::transport::mixture::{
    mass_to_mole_fractions, mixture_diffusion_coefficients, mixture_thermal_conductivity,
};

pub struct FlameConfig {
    pub pressure: f64,        // Pa
    pub t_unburned: f64,      // K
    pub y_unburned: Vec<f64>, // mass fractions at inlet
    pub z_fix: f64,           // position where T is pinned [m]
    pub t_fix: f64,           // pinned temperature [K]
}

/// Evaluate residual F(x) → rhs (length = solution_length).
///
/// `x_old` and `rdt` enable pseudo-transient continuation:
///   rhs[i] += rdt * (x[i] - x_old[i])   for all interior equations.
/// Pass `x_old = None` (and `rdt = 0.0`) for steady-state evaluation.
pub fn eval_residual(
    x: &[f64],
    rhs: &mut [f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    x_old: Option<&[f64]>,
    rdt: f64,
) {
    let nk = mech.n_species();
    let nj = grid.n_points();
    let nv = natj(mech);
    let p = config.pressure;

    rhs.iter_mut().for_each(|r| *r = 0.0);

    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();

    let dz = grid.dz();

    // Pre-compute midpoint transport: λ[j] and Dk[k][j] at interval j (between j and j+1)
    let mut lambda_mid = vec![0.0_f64; nj - 1];
    let mut dk_mid = vec![vec![0.0_f64; nj - 1]; nk];
    let mut rho_mid = vec![0.0_f64; nj - 1];

    for j in 0..nj - 1 {
        let t_av = 0.5 * (state.temperature(j) + state.temperature(j + 1));
        let y_av: Vec<f64> = (0..nk)
            .map(|k| 0.5 * (state.species(k, j) + state.species(k, j + 1)))
            .collect();
        let x_av = mass_to_mole_fractions(mech, &y_av);
        let w_mean = mean_molecular_weight(&mech.species, &y_av);
        rho_mid[j] = density(p, t_av, w_mean);
        lambda_mid[j] = mixture_thermal_conductivity(mech, &x_av, &y_av, t_av, p);
        let dk = mixture_diffusion_coefficients(mech, &x_av, t_av, p);
        for k in 0..nk {
            dk_mid[k][j] = dk[k];
        }
    }

    // Diffusion fluxes jk [kg/(m²·s)] at each midpoint j (between j and j+1).
    // Mixture-averaged with correction velocity to enforce Σjk = 0.
    let mut jk_mid = vec![vec![0.0_f64; nj - 1]; nk];

    for j in 0..nj - 1 {
        let mut jk_raw = vec![0.0_f64; nk];
        let mut sum_jk = 0.0_f64;
        for k in 0..nk {
            let dy_dz = (state.species(k, j + 1) - state.species(k, j)) / dz[j];
            jk_raw[k] = -rho_mid[j] * dk_mid[k][j] * dy_dz;
            sum_jk += jk_raw[k];
        }
        for k in 0..nk {
            let yk_av = 0.5 * (state.species(k, j) + state.species(k, j + 1));
            jk_mid[k][j] = jk_raw[k] - yk_av * sum_jk;
        }
    }

    // --- Left boundary (j = 0) ---
    {
        let j = 0;
        rhs[idx_t(nv, j)] = state.temperature(j) - config.t_unburned;

        for k in 0..nk {
            let yk = state.species(k, j);
            // Inlet flux BC: M*(Yk - Yku) + jk_{1/2} = 0
            rhs[idx_y(nv, j, k)] = m * (yk - config.y_unburned[k]) + jk_mid[k][0];
        }
    }

    // --- Interior points j = 1 … J-2 ---
    let j_fix = find_fixed_point(grid, config.z_fix);

    for j in 1..nj - 1 {
        let t_j   = state.temperature(j);
        let t_jm1 = state.temperature(j - 1);
        let t_jp1 = state.temperature(j + 1);
        let dz_m = dz[j - 1]; // dz between j-1 and j
        let dz_p = dz[j];     // dz between j and j+1
        let dz_av = 0.5 * (dz_m + dz_p);

        let y_j: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
        let w_mean = mean_molecular_weight(&mech.species, &y_j);
        let rho_j = density(p, t_j, w_mean);

        // Chemical production rates: concentrations in mol/cm³ (CGS, consistent with A values)
        let conc: Vec<f64> = (0..nk)
            .map(|k| rho_j * y_j[k] / mech.species[k].molecular_weight * 1e-6)
            .collect();
        let wdot = production_rates(mech, t_j, &conc, p); // mol/(cm³·s)

        // Species equations
        let mut sum_yk = 0.0_f64;
        for k in 0..nk - 1 {
            let yk_j   = state.species(k, j);
            let yk_jm1 = state.species(k, j - 1);
            let convec = m * (yk_j - yk_jm1) / dz_m;
            let diffus = (jk_mid[k][j] - jk_mid[k][j - 1]) / dz_av;
            // wdot mol/(cm³·s) × 1e6 → mol/(m³·s) × Wk → kg/(m³·s)
            let source = wdot[k] * mech.species[k].molecular_weight * 1e6;
            rhs[idx_y(nv, j, k)] = convec + diffus - source;
            sum_yk += yk_j;
        }
        // Last species: mass fraction sum closure
        rhs[idx_y(nv, j, nk - 1)] = sum_yk + state.species(nk - 1, j) - 1.0;

        // Energy equation
        let cp_j       = cp_mixture(&mech.species, &y_j, t_j);
        let lambda_j   = lambda_mid[j];      // midpoint between j and j+1
        let lambda_jm1 = lambda_mid[j - 1];  // midpoint between j-1 and j

        let conduction = (lambda_j * (t_jp1 - t_j) / dz_p
            - lambda_jm1 * (t_j - t_jm1) / dz_m) / dz_av;

        // Enthalpy transport: Σk jk * cpk * dT/dz  evaluated at point j.
        // Average the diffusion fluxes at the two adjacent midpoints (j-1/2 and j+1/2)
        // for second-order accuracy, consistent with Cantera's FreeFlame formulation.
        let dt_dz_m = (t_j - t_jm1) / dz_m;
        let mut enthalpy_transport = 0.0_f64;
        for k in 0..nk {
            let cp_k = cp_species(&mech.species[k], t_j);
            let jk_avg = 0.5 * (jk_mid[k][j - 1] + jk_mid[k][j]);
            enthalpy_transport += jk_avg * cp_k * dt_dz_m;
        }

        // Heat release: Σk ωk * hk [W/m³]
        // wdot mol/(cm³·s) × 1e6 → mol/(m³·s); hk J/mol → W/m³
        let heat_release: f64 = (0..nk)
            .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t_j))
            .sum::<f64>() * 1e6;

        let convec_t = m * cp_j * (t_j - t_jm1) / dz_m;
        rhs[idx_t(nv, j)] = -convec_t + conduction - enthalpy_transport - heat_release;

        // Eigenvalue closure: fix temperature at j_fix to locate the flame
        if j == j_fix {
            rhs[idx_m(nv, nj)] = t_j - config.t_fix;
        }
    }

    // --- Right boundary (j = J-1) ---
    {
        let j = nj - 1;
        rhs[idx_t(nv, j)] = state.temperature(j) - state.temperature(j - 1);
        for k in 0..nk {
            rhs[idx_y(nv, j, k)] = state.species(k, j) - state.species(k, j - 1);
        }
    }

    // --- Pseudo-transient embedding ---
    // Adds rdt*(x - x_old) to all interior variable equations.
    if let Some(x_old) = x_old {
        if rdt > 0.0 {
            for j in 1..nj - 1 {
                let it = idx_t(nv, j);
                rhs[it] += rdt * (x[it] - x_old[it]);
                for k in 0..nk {
                    let iy = idx_y(nv, j, k);
                    rhs[iy] += rdt * (x[iy] - x_old[iy]);
                }
            }
        }
    }
}

/// Find the grid index closest to z_fix.
fn find_fixed_point(grid: &Grid, z_fix: f64) -> usize {
    grid.z.iter().enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - z_fix).abs()).partial_cmp(&((*b - z_fix).abs())).unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(grid.n_points() / 2)
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

    // -----------------------------------------------------------------------
    // Residual ≈ 0 for a uniform steady-state profile with no reactions.
    //
    // State: T = 1000 K (uniform), Y_N2 = 1 (all other species = 0),
    //        M = 0.2 kg/(m²·s).
    //
    // Why this should give F ≈ 0:
    //   - Uniform T and Y → all convection and diffusion terms = 0
    //   - N2 is not a reactant/product in h2o2.yaml → ωk = 0 for all k
    //   - Left BC: M*(Y_N2 - 1) + jk_0 = 0 + 0 = 0
    //   - Right BC: zero gradient trivially satisfied
    //   - Eigenvalue: T(j_fix) - T_fix = 0
    // -----------------------------------------------------------------------
    #[test]
    fn test_residual_zero_for_uniform_n2() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 6;
        let nv = natj(&mech);
        let grid = Grid::uniform(0.02, nj);

        let n2_idx = mech.species_index("N2").unwrap();

        // Build solution vector: T = 1000 K, Y_N2 = 1, M = 0.2
        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            x[idx_t(nv, j)] = 1000.0;
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = 0.2;

        // y_unburned = pure N2 at 1000 K; z_fix and t_fix consistent with profile
        let mut y_unburned = vec![0.0_f64; nk];
        y_unburned[n2_idx] = 1.0;
        let config = FlameConfig {
            pressure: 101325.0,
            t_unburned: 1000.0,
            y_unburned,
            z_fix: grid.z[nj / 2],
            t_fix: 1000.0,
        };

        let mut rhs = vec![0.0_f64; solution_length(&mech, nj)];
        eval_residual(&x, &mut rhs, &mech, &grid, &config, None, 0.0);

        let max_abs = rhs.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_abs < 1e-8,
            "Max |F| = {max_abs:.3e} (expected < 1e-8 for uniform steady N2 profile)"
        );
    }

    // -----------------------------------------------------------------------
    // Pseudo-transient term: for x_old = x, rdt * (x - x_old) = 0.
    // So residual with pseudo-transient should equal steady-state residual.
    // -----------------------------------------------------------------------
    #[test]
    fn test_pseudo_transient_zero_when_x_eq_x_old() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 6;
        let nv = natj(&mech);
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();

        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            x[idx_t(nv, j)] = 1000.0;
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = 0.2;

        let mut y_unburned = vec![0.0_f64; nk];
        y_unburned[n2_idx] = 1.0;
        let config = FlameConfig {
            pressure: 101325.0,
            t_unburned: 1000.0,
            y_unburned,
            z_fix: grid.z[nj / 2],
            t_fix: 1000.0,
        };

        let mut rhs_steady = vec![0.0_f64; solution_length(&mech, nj)];
        let mut rhs_pt     = vec![0.0_f64; solution_length(&mech, nj)];
        eval_residual(&x, &mut rhs_steady, &mech, &grid, &config, None, 0.0);
        eval_residual(&x, &mut rhs_pt,     &mech, &grid, &config, Some(&x.clone()), 1e4);

        for (a, b) in rhs_steady.iter().zip(rhs_pt.iter()) {
            assert!((a - b).abs() < 1e-15, "PT residual differs from steady: {a:.3e} vs {b:.3e}");
        }
    }
}
