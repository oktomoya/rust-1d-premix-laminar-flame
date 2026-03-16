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
    cp_mixture, density, enthalpy_molar, mean_molecular_weight,
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

    // Pre-compute midpoint transport: λ[j], Dk[k][j], ρ[j], W̄[j]
    // at interval j (between grid points j and j+1).
    let mut lambda_mid  = vec![0.0_f64; nj - 1];
    let mut dk_mid      = vec![vec![0.0_f64; nj - 1]; nk];
    let mut rho_mid     = vec![0.0_f64; nj - 1];
    let mut wmean_mid   = vec![0.0_f64; nj - 1];

    for j in 0..nj - 1 {
        let t_av = 0.5 * (state.temperature(j) + state.temperature(j + 1));
        let y_av: Vec<f64> = (0..nk)
            .map(|k| 0.5 * (state.species(k, j) + state.species(k, j + 1)))
            .collect();
        let x_av = mass_to_mole_fractions(mech, &y_av);
        let w_mean = mean_molecular_weight(&mech.species, &y_av);
        rho_mid[j]   = density(p, t_av, w_mean);
        wmean_mid[j] = w_mean;
        lambda_mid[j] = mixture_thermal_conductivity(mech, &x_av, &y_av, t_av, p);
        let dk = mixture_diffusion_coefficients(mech, &x_av, &y_av, t_av, p);
        for k in 0..nk {
            dk_mid[k][j] = dk[k];
        }
    }

    // Mole fractions at each grid point (needed for mole-basis flux gradient).
    let xk_grid: Vec<Vec<f64>> = (0..nj)
        .map(|j| {
            let y_j: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
            mass_to_mole_fractions(mech, &y_j)
        })
        .collect();

    // Diffusion fluxes jk [kg/(m²·s)] at each midpoint j (between j and j+1).
    //
    // Mole-fraction basis (matches PREMIX and Cantera):
    //   jk_raw = -ρ (Wk/W̄) Dkm dXk/dz
    //
    // Correction velocity enforces Σjk = 0:
    //   jk = jk_raw - Yk * Σj jk_raw_j
    let mut jk_mid = vec![vec![0.0_f64; nj - 1]; nk];

    for j in 0..nj - 1 {
        let mut jk_raw = vec![0.0_f64; nk];
        let mut sum_jk = 0.0_f64;
        for k in 0..nk {
            let dx_dz = (xk_grid[j + 1][k] - xk_grid[j][k]) / dz[j];
            let wk_over_wmean = mech.species[k].molecular_weight / wmean_mid[j];
            jk_raw[k] = -rho_mid[j] * wk_over_wmean * dk_mid[k][j] * dx_dz;
            sum_jk += jk_raw[k];
        }
        for k in 0..nk {
            // Use left-cell Yk for the correction, consistent with upwind convection.
            jk_mid[k][j] = jk_raw[k] - state.species(k, j) * sum_jk;
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

        // Enthalpy transport: Σk jk * dhk/dz  (matches Cantera Flow1D exactly)
        //
        // dhk/dz uses upwind differencing (u > 0 for free flame):
        //   dhk/dz ≈ (hk(j) - hk(j-1)) / dz_m   [J/(kg·m)]
        //
        // jk averaged over adjacent midpoints (j-1/2) and (j+1/2).
        let mut enthalpy_transport = 0.0_f64;
        for k in 0..nk {
            let hk_j   = enthalpy_molar(&mech.species[k], t_j)   / mech.species[k].molecular_weight;
            let hk_jm1 = enthalpy_molar(&mech.species[k], t_jm1) / mech.species[k].molecular_weight;
            let dhk_dz = (hk_j - hk_jm1) / dz_m;
            let jk_avg = 0.5 * (jk_mid[k][j - 1] + jk_mid[k][j]);
            enthalpy_transport += jk_avg * dhk_dz;
        }

        // Heat release: Σk ωk * hk [W/m³]
        // wdot mol/(cm³·s) × 1e6 → mol/(m³·s); hk J/mol → W/m³
        let heat_release: f64 = (0..nk)
            .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t_j))
            .sum::<f64>() * 1e6;

        let convec_t = m * cp_j * (t_j - t_jm1) / dz_m;
        rhs[idx_t(nv, j)] = -convec_t + conduction - enthalpy_transport - heat_release;

        // Eigenvalue closure: fix temperature at j_fix to locate the flame.
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
    // Backward-Euler unsteady term for interior variables:
    //   T:   ρ·cp · (T - T_old) / dt   [W/m³]   → adds  rdt·ρ·cp·ΔT
    //   Yk:  ρ    · (Yk - Yk_old) / dt [kg/(m³·s)] → adds  rdt·ρ·ΔYk
    // The density and cp factors are essential for dimensional consistency.
    if let Some(x_old) = x_old {
        if rdt > 0.0 {
            for j in 1..nj - 1 {
                let t_j = x[idx_t(nv, j)];
                let y_j: Vec<f64> = (0..nk).map(|k| x[idx_y(nv, j, k)]).collect();
                let w_j = mean_molecular_weight(&mech.species, &y_j);
                let rho_j = density(p, t_j, w_j);
                let cp_j = cp_mixture(&mech.species, &y_j, t_j);

                let it = idx_t(nv, j);
                rhs[it] += rdt * rho_j * cp_j * (x[it] - x_old[it]);
                for k in 0..nk {
                    let iy = idx_y(nv, j, k);
                    rhs[iy] += rdt * rho_j * (x[iy] - x_old[iy]);
                }
            }
        }
    }
}

/// Evaluate residuals for only a contiguous range of grid points.
///
/// Only grid points `j_lo..j_hi` (exclusive) are computed and written to `rhs`.
/// Transport midpoints adjacent to the range are computed on-demand; all other
/// midpoints are skipped.  Entries of `rhs` outside the range are **not** touched.
///
/// Used by the Jacobian builder to limit work to the local stencil: when column
/// `col` corresponds to grid point `gp = col/nv`, only points `gp−1..gp+2` can
/// have nonzero partial derivatives.
///
/// The caller must ensure `rhs[row_min..row_max]` is pre-filled with the baseline
/// residual `f0` so that entries for rows outside the range correctly evaluate to 0
/// when differenced against `f0`.
pub fn eval_residual_range(
    x: &[f64],
    rhs: &mut [f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    j_lo: usize,
    j_hi: usize,
) {
    let nk = mech.n_species();
    let nj = grid.n_points();
    let nv = natj(mech);
    let p = config.pressure;
    let dz = grid.dz();

    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();

    // Midpoint range: for interior point j we need midpoints j-1 and j.
    // Covering [j_lo, j_hi) requires midpoints [m_lo, m_hi).
    let m_lo: usize = j_lo.saturating_sub(1);
    let m_hi: usize = j_hi.min(nj - 1); // exclusive; nj-1 midpoints total
    let nm = m_hi.saturating_sub(m_lo);

    // Compact arrays indexed by (midpoint_index - m_lo).
    let mut lambda_local = vec![0.0_f64; nm];
    let mut dk_local     = vec![vec![0.0_f64; nm]; nk];
    let mut rho_local    = vec![0.0_f64; nm];
    let mut wmean_local  = vec![0.0_f64; nm];

    for jm in m_lo..m_hi {
        let mi = jm - m_lo;
        let t_av = 0.5 * (state.temperature(jm) + state.temperature(jm + 1));
        let y_av: Vec<f64> = (0..nk)
            .map(|k| 0.5 * (state.species(k, jm) + state.species(k, jm + 1)))
            .collect();
        let x_av = mass_to_mole_fractions(mech, &y_av);
        let w_mean = mean_molecular_weight(&mech.species, &y_av);
        rho_local[mi]    = density(p, t_av, w_mean);
        wmean_local[mi]  = w_mean;
        lambda_local[mi] = mixture_thermal_conductivity(mech, &x_av, &y_av, t_av, p);
        let dk = mixture_diffusion_coefficients(mech, &x_av, &y_av, t_av, p);
        for k in 0..nk { dk_local[k][mi] = dk[k]; }
    }

    // Mole fractions for grid points [m_lo, m_hi+1) — both endpoints of each midpoint.
    let g_lo = m_lo;
    let g_hi = (m_hi + 1).min(nj);
    let ng = g_hi - g_lo;
    let mut xk_local: Vec<Vec<f64>> = Vec::with_capacity(ng);
    for j in g_lo..g_hi {
        let y_j: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
        xk_local.push(mass_to_mole_fractions(mech, &y_j));
    }
    let xkg = |j: usize| &xk_local[j - g_lo];

    // Diffusion fluxes jk_mid[k][mi] at midpoint m_lo + mi.
    let mut jk_local = vec![vec![0.0_f64; nm]; nk];
    for jm in m_lo..m_hi {
        let mi = jm - m_lo;
        let mut jk_raw = vec![0.0_f64; nk];
        let mut sum_jk = 0.0_f64;
        for k in 0..nk {
            let dx_dz = (xkg(jm + 1)[k] - xkg(jm)[k]) / dz[jm];
            let wk_over_wmean = mech.species[k].molecular_weight / wmean_local[mi];
            jk_raw[k] = -rho_local[mi] * wk_over_wmean * dk_local[k][mi] * dx_dz;
            sum_jk += jk_raw[k];
        }
        for k in 0..nk {
            jk_local[k][mi] = jk_raw[k] - state.species(k, jm) * sum_jk;
        }
    }
    // Accessor: diffusion flux for species k at midpoint jm (global index).
    let jkm = |k: usize, jm: usize| jk_local[k][jm - m_lo];

    let j_fix = find_fixed_point(grid, config.z_fix);

    for j in j_lo..j_hi {
        // Zero only the entries for this grid point.
        rhs[j * nv..(j + 1) * nv].fill(0.0);

        if j == 0 {
            // Left boundary
            rhs[idx_t(nv, 0)] = state.temperature(0) - config.t_unburned;
            for k in 0..nk {
                rhs[idx_y(nv, 0, k)] =
                    m * (state.species(k, 0) - config.y_unburned[k]) + jkm(k, 0);
            }
        } else if j == nj - 1 {
            // Right boundary (zero gradient)
            rhs[idx_t(nv, j)] = state.temperature(j) - state.temperature(j - 1);
            for k in 0..nk {
                rhs[idx_y(nv, j, k)] = state.species(k, j) - state.species(k, j - 1);
            }
        } else {
            // Interior
            let t_j   = state.temperature(j);
            let t_jm1 = state.temperature(j - 1);
            let t_jp1 = state.temperature(j + 1);
            let dz_m = dz[j - 1];
            let dz_p = dz[j];
            let dz_av = 0.5 * (dz_m + dz_p);

            let y_j: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
            let w_mean = mean_molecular_weight(&mech.species, &y_j);
            let rho_j  = density(p, t_j, w_mean);

            let conc: Vec<f64> = (0..nk)
                .map(|k| rho_j * y_j[k] / mech.species[k].molecular_weight * 1e-6)
                .collect();
            let wdot = production_rates(mech, t_j, &conc, p);

            let mut sum_yk = 0.0_f64;
            for k in 0..nk - 1 {
                let convec = m * (state.species(k, j) - state.species(k, j - 1)) / dz_m;
                let diffus = (jkm(k, j) - jkm(k, j - 1)) / dz_av;
                let source = wdot[k] * mech.species[k].molecular_weight * 1e6;
                rhs[idx_y(nv, j, k)] = convec + diffus - source;
                sum_yk += state.species(k, j);
            }
            rhs[idx_y(nv, j, nk - 1)] = sum_yk + state.species(nk - 1, j) - 1.0;

            let cp_j       = cp_mixture(&mech.species, &y_j, t_j);
            let lambda_j   = lambda_local[j - m_lo];
            let lambda_jm1 = lambda_local[j - 1 - m_lo];
            let conduction = (lambda_j * (t_jp1 - t_j) / dz_p
                - lambda_jm1 * (t_j - t_jm1) / dz_m) / dz_av;

            let mut enthalpy_transport = 0.0_f64;
            for k in 0..nk {
                let hk_j   = enthalpy_molar(&mech.species[k], t_j)   / mech.species[k].molecular_weight;
                let hk_jm1 = enthalpy_molar(&mech.species[k], t_jm1) / mech.species[k].molecular_weight;
                let dhk_dz = (hk_j - hk_jm1) / dz_m;
                let jk_avg = 0.5 * (jkm(k, j - 1) + jkm(k, j));
                enthalpy_transport += jk_avg * dhk_dz;
            }

            let heat_release: f64 = (0..nk)
                .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t_j))
                .sum::<f64>() * 1e6;

            let convec_t = m * cp_j * (t_j - t_jm1) / dz_m;
            rhs[idx_t(nv, j)] = -convec_t + conduction - enthalpy_transport - heat_release;

            if j == j_fix {
                rhs[idx_m(nv, nj)] = t_j - config.t_fix;
            }
        }
    }
}

/// Return the solution-vector column index of T at the fixed-point grid node.
///
/// Used by the Jacobian rank-2 Woodbury correction to handle the out-of-band
/// entry ∂F[n-1]/∂T[j_fix] = 1 in the eigenvalue row.
pub fn j_fix_t_col(mech: &Mechanism, grid: &Grid, config: &FlameConfig) -> usize {
    let nv = natj(mech);
    let j_fix = find_fixed_point(grid, config.z_fix);
    idx_t(nv, j_fix)
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
