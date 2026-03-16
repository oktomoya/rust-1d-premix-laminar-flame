/// Finite-difference Jacobian for the bordered-banded flame solver.
///
/// The flame system has n = nj·nv + 1 unknowns (state variables + mass-flux
/// eigenvalue M).  The equations are:
///   F[0..n-2]: interior/boundary residuals (local stencil, bandwidth kl=ku=2·nv)
///   F[n-1]:    eigenvalue closure  T[j_fix] = t_fix
///
/// This module builds the (n-1)×(n-1) **interior Jacobian** J_int = ∂F[0..n-2]/∂x[0..n-2]
/// and the separate M-column vector m_col = ∂F[0..n-2]/∂M.  The eigenvalue equation
/// and M are handled by `bordered_solve_m` using a Schur-complement bordered solve,
/// which avoids the Sherman-Morrison singularity that arises when the M-row diagonal
/// is exactly zero.

use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::residual::{eval_residual, eval_residual_range, j_fix_t_col, FlameConfig};
use crate::flame::state::solution_length;
use crate::solver::banded::BandedMatrix;

/// Build the full (n-1)×(n-1) interior Jacobian via full `eval_residual` per column.
/// O(n²) cost — only used for cross-checking.
#[cfg(test)]
fn numerical_jacobian_full(
    x: &[f64],
    f0: &[f64],
    mech: &crate::chemistry::mechanism::Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) -> Vec<Vec<f64>> {
    use crate::flame::state::solution_length;
    let n = solution_length(mech, grid.n_points());
    let n_int = n - 1;
    let eps_rel = (2.0_f64 * f64::EPSILON).sqrt();
    let mut jac = vec![vec![0.0_f64; n_int]; n_int];
    let mut x_pert = x.to_vec();
    let mut f_pert = vec![0.0_f64; n];
    for col in 0..n_int {
        let x_col = x[col];
        let h = eps_rel * x_col.abs().max(1.0);
        x_pert[col] = x_col + h;
        eval_residual(&x_pert, &mut f_pert, mech, grid, config, None, 0.0);
        for i in 0..n_int {
            jac[i][col] = (f_pert[i] - f0[i]) / h;
        }
        x_pert[col] = x_col;
    }
    jac
}

/// Build the (n-1)×(n-1) interior Jacobian J_int and the M-column vector.
///
/// Returns `(jac, m_col, jfc)` where:
/// - `jac`:   banded (n-1)×(n-1) matrix  ∂F[0..n-2]/∂x[0..n-2].
/// - `m_col`: vector of length n-1,       ∂F[0..n-2]/∂M.
/// - `jfc`:   column index of T[j_fix] within x[0..n-2] (needed for the Schur step).
pub fn numerical_jacobian(
    x: &[f64],
    f0: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) -> (BandedMatrix, Vec<f64>, usize) {
    let n = solution_length(mech, grid.n_points());
    let n_int = n - 1; // interior unknowns (all except M)
    let nj = grid.n_points();
    let nv = mech.n_species() + 1;
    let kl = 2 * nv;
    let ku = 2 * nv;

    let mut jac = BandedMatrix::new(n_int, kl, ku);
    let mut x_pert = x.to_vec();
    let mut f_pert = f0.to_vec();
    let mut m_col = vec![0.0_f64; n_int];

    let eps_rel = (2.0_f64 * f64::EPSILON).sqrt();
    let m_col_idx = n - 1; // index of M in x

    for col in 0..n_int {
        let x_col = x[col];
        let h = eps_rel * x_col.abs().max(1.0);

        let row_min = col.saturating_sub(ku);
        let row_max = (col + kl + 1).min(n_int);

        // Reset within-band interior rows to f0.
        f_pert[row_min..row_max].copy_from_slice(&f0[row_min..row_max]);
        // Also reset the eigenvalue row so it doesn't bleed between columns.
        f_pert[m_col_idx] = f0[m_col_idx];

        x_pert[col] = x_col + h;

        let gp = col / nv;
        let j_lo = gp.saturating_sub(1);
        let j_hi = (gp + 2).min(nj);
        eval_residual_range(&x_pert, &mut f_pert, mech, grid, config, j_lo, j_hi);

        for i in row_min..row_max {
            jac.set(i, col, (f_pert[i] - f0[i]) / h);
        }

        x_pert[col] = x_col;
    }

    // M column: one full residual eval with M perturbed.
    {
        let x_m = x[m_col_idx];
        let h = eps_rel * x_m.abs().max(1.0);
        x_pert[m_col_idx] = x_m + h;
        eval_residual(&x_pert, &mut f_pert, mech, grid, config, None, 0.0);
        for i in 0..n_int {
            m_col[i] = (f_pert[i] - f0[i]) / h;
        }
        x_pert[m_col_idx] = x_m;
    }

    let jfc = j_fix_t_col(mech, grid, config);
    (jac, m_col, jfc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;
    use crate::flame::domain::Grid;
    use crate::flame::residual::FlameConfig;
    use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};

    fn h2o2_mech() -> crate::chemistry::mechanism::Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    /// Verify that `BandedMatrix::solve_factored` gives the correct solution for the
    /// realistic flame Jacobian (n≈90, kl=ku=2*nv=22).  Builds J, computes b=J*x_exact
    /// for a known x_exact, then solves J*x=b and checks x==x_exact.
    #[test]
    fn test_banded_solve_correct_for_flame_jacobian() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 8; // small but realistic bandwidth
        let nv = natj(&mech);
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();
        let h2_idx = mech.species_index("H2").unwrap();

        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            let s = j as f64 / (nj - 1) as f64;
            x[idx_t(nv, j)] = 300.0 + s * 1500.0;
            x[idx_y(nv, j, n2_idx)] = 0.7 - 0.2 * s;
            x[idx_y(nv, j, h2_idx)] = 0.3 + 0.2 * s;
        }
        x[idx_m(nv, nj)] = 0.3;

        let mut y_unburned = vec![0.0_f64; nk];
        y_unburned[n2_idx] = 0.7;
        y_unburned[h2_idx] = 0.3;
        let config = FlameConfig {
            pressure: 101325.0,
            t_unburned: 300.0,
            y_unburned,
            z_fix: grid.z[nj / 2],
            t_fix: 900.0,
        };

        let n = solution_length(&mech, nj);
        let n_int = n - 1;
        let mut f0 = vec![0.0_f64; n];
        eval_residual(&x, &mut f0, &mech, &grid, &config, None, 0.0);

        let (mut jac, _, _) = numerical_jacobian(&x, &f0, &mech, &grid, &config);

        // Build a known rhs: b = J * x_exact where x_exact = [1, 2, 3, ...]
        let x_exact: Vec<f64> = (0..n_int).map(|i| (i as f64 + 1.0) * 0.001).collect();

        // Compute b = J * x_exact using the dense reference.
        let jac_ref = numerical_jacobian_full(&x, &f0, &mech, &grid, &config);
        let mut b = vec![0.0_f64; n_int];
        for i in 0..n_int {
            for j in 0..n_int {
                b[i] += jac_ref[i][j] * x_exact[j];
            }
        }

        // Solve J * x_sol = b using the banded factorization.
        jac.factor_in_place().expect("factorization must succeed");
        let mut x_sol = b.clone();
        jac.solve_factored(&mut x_sol);

        // Check x_sol ≈ x_exact.
        let max_err = x_sol.iter().zip(x_exact.iter())
            .map(|(s, e)| ((s - e) / e.abs().max(1e-10)).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_err < 1e-6,
            "Banded solve relative error = {max_err:.3e} (should be < 1e-6). \
             Banded LU factorization may be incorrect for realistic flame Jacobian."
        );
    }

    /// Verify that `numerical_jacobian` (using `eval_residual_range`) produces
    /// the same entries as the reference dense Jacobian built from full `eval_residual`.
    #[test]
    fn test_jacobian_range_matches_full() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 8;
        let nv = natj(&mech);
        // Use a non-trivial profile with a temperature gradient so Jacobian entries
        // are non-zero across the bandwidth.
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();
        let h2_idx = mech.species_index("H2").unwrap();

        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            let s = j as f64 / (nj - 1) as f64;
            x[idx_t(nv, j)] = 300.0 + s * 1500.0;  // 300→1800 K ramp
            x[idx_y(nv, j, n2_idx)] = 0.7 - 0.2 * s;
            x[idx_y(nv, j, h2_idx)] = 0.3 + 0.2 * s;
        }
        x[idx_m(nv, nj)] = 0.3;

        let mut y_unburned = vec![0.0_f64; nk];
        y_unburned[n2_idx] = 0.7;
        y_unburned[h2_idx] = 0.3;
        let config = FlameConfig {
            pressure: 101325.0,
            t_unburned: 300.0,
            y_unburned,
            z_fix: grid.z[nj / 2],
            t_fix: 900.0,
        };

        let n = solution_length(&mech, nj);
        let n_int = n - 1;
        let mut f0 = vec![0.0_f64; n];
        eval_residual(&x, &mut f0, &mech, &grid, &config, None, 0.0);

        // Build banded Jacobian via range (current approach).
        let (jac_band, _, _) = numerical_jacobian(&x, &f0, &mech, &grid, &config);

        // Build dense reference via full eval_residual.
        let jac_ref = numerical_jacobian_full(&x, &f0, &mech, &grid, &config);

        // Compare in-band entries.
        let kl = 2 * nv;
        let ku = 2 * nv;
        let mut max_rel_err = 0.0_f64;
        let mut worst = (0usize, 0usize);
        for col in 0..n_int {
            let row_min = col.saturating_sub(ku);
            let row_max = (col + kl + 1).min(n_int);
            for row in row_min..row_max {
                let banded = jac_band.get(row, col);
                let reference = jac_ref[row][col];
                let scale = reference.abs().max(1.0);
                let rel_err = (banded - reference).abs() / scale;
                if rel_err > max_rel_err {
                    max_rel_err = rel_err;
                    worst = (row, col);
                }
            }
        }
        assert!(
            max_rel_err < 1e-6,
            "Max relative error in banded Jacobian = {max_rel_err:.3e} at ({}, {}): \
             banded={:.6e}, ref={:.6e}",
            worst.0, worst.1,
            jac_band.get(worst.0, worst.1),
            jac_ref[worst.0][worst.1],
        );
    }
}

/// Bordered-banded solve for the mass-flux eigenvalue M.
///
/// The full n×n flame system is split as a bordered banded system:
///
///   [ J_int  | m_col ] [ step_int ]   [ -F_int ]
///   [ e_jfc^T|  d    ] [ δM       ] = [ -F_ev  ]
///
/// where:
///   J_int   = (n-1)×(n-1) interior Jacobian (already LU-factored)
///   m_col   = ∂F_int/∂M (M column for interior rows)
///   e_jfc^T = eigenvalue row (only entry at column jfc = T[j_fix] index)
///   d       = ∂F_ev/∂M  (0 for steady Newton; rdt for pseudo-transient)
///   -F_ev   = -(T[j_fix] − t_fix)  [+ rdt·(M − M_old) for PT]
///
/// On entry  `step_int` = J_int^{-1} · (−F_int)  (already solved by caller).
/// On exit   `step_int` is updated with the M correction, and δM is returned.
///
/// **Algorithm** (Schur complement):
///   z = J_int^{-1} · m_col
///   S = d − e_jfc^T · z  =  d − z[jfc]
///   δM = (−F_ev − step_int[jfc]) / S
///   step_int -= z · δM
pub fn bordered_solve_m(
    step_int: &mut [f64],
    m_col: &[f64],
    jac: &BandedMatrix,
    jfc: usize,
    neg_f_ev: f64,
    d: f64,
) -> f64 {
    let mut z = m_col.to_vec();
    jac.solve_factored(&mut z);

    let s = d - z[jfc];
    if s.abs() < 1e-300 { return 0.0; }

    let dm = (neg_f_ev - step_int[jfc]) / s;
    for i in 0..step_int.len() {
        step_int[i] -= z[i] * dm;
    }
    dm
}
