/// Damped Newton iteration with backtracking line search.

use anyhow::{bail, Result};
use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::jacobian::numerical_jacobian_dense;
use crate::flame::residual::{eval_residual, FlameConfig};
use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};
use crate::solver::dense::DenseMatrix;

pub struct NewtonConfig {
    pub atol: f64,
    pub rtol: f64,
    pub max_iter: usize,
    /// Max number of iterations before the Jacobian is recomputed.
    pub max_jac_age: usize,
    /// Accept the current solution if the line search fails to improve the
    /// residual for this many consecutive iterations (coarse-grid stall).
    pub max_stall_iter: usize,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        NewtonConfig {
            atol: 1e-9,
            rtol: 1e-6,
            max_iter: 50,
            max_jac_age: 3,
            max_stall_iter: 5,
        }
    }
}

/// Solve F(x) = 0 starting from initial guess x (modified in-place).
pub fn solve(
    x: &mut Vec<f64>,
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    newton: &NewtonConfig,
) -> Result<()> {
    let n = solution_length(mech, grid.n_points());
    let mut f = vec![0.0_f64; n];
    // Cached Jacobian; reused for up to max_jac_age iterations.
    let mut cached_jac: Option<DenseMatrix> = None;
    let mut jac_age = newton.max_jac_age + 1; // force evaluation on first iteration
    let mut stall_count = 0usize;
    let mut best_norm_seen = f64::INFINITY;
    let mut prev_norm_f = f64::INFINITY;

    for iter in 0..newton.max_iter {
        eval_residual(x, &mut f, mech, grid, config, None, 0.0);

        let norm_f = norm2(&f);
        eprintln!("Newton iter {iter:3}: ‖F‖ = {norm_f:.3e}");

        // Scaled convergence: ‖F‖ < atol·√n + rtol·‖x‖
        if norm_f < newton.atol * (n as f64).sqrt() + newton.rtol * norm2(x) {
            eprintln!("Newton converged in {iter} iterations");
            return Ok(());
        }

        // Force Jacobian recompute if the residual dropped by >10x since the
        // last evaluation: the old Jacobian is from a very different linearisation
        // point and would give a poor (possibly diverging) step.
        if norm_f < prev_norm_f / 10.0 {
            jac_age = newton.max_jac_age + 1;
        }
        prev_norm_f = norm_f;

        // (Re)compute Jacobian when stale or absent.
        if jac_age >= newton.max_jac_age || cached_jac.is_none() {
            cached_jac = Some(numerical_jacobian_dense(x, &f, mech, grid, config));
            jac_age = 0;
        }

        // Solve J·step = -F.
        let mut step: Vec<f64> = f.iter().map(|v| -v).collect();
        cached_jac.as_ref().unwrap().solve(&mut step)?;

        // Clip the step to physical scale limits (T ≤ 500 K, Y ≤ 0.5, etc.)
        // before the line search.  This prevents the first trial point from
        // being wildly unphysical when the Jacobian is far from the solution.
        let clip = step_clip_factor(x, &step, mech, grid);
        if clip < 1.0 {
            for s in step.iter_mut() { *s *= clip; }
        }

        // Backtracking line search.
        let (alpha, improved) = line_search(x, &step, mech, grid, config, norm_f)?;
        for i in 0..n {
            x[i] += alpha * step[i];
        }

        // Prevent NaN propagation from large/unphysical steps.
        clamp_physical(x, mech, grid);

        // Invalidate cached Jacobian if line search cut the step or no improvement.
        if !improved || alpha < 0.5 {
            jac_age = newton.max_jac_age + 1;
        } else {
            jac_age += 1;
        }

        // Stall detection: if the residual hasn't improved by >1% relative to
        // the best seen in the last max_stall_iter iterations, the solver has
        // reached the truncation-error floor for this grid.  Accept the current
        // solution and let the caller refine the grid.
        if norm_f < best_norm_seen * (1.0 - 0.01) {
            best_norm_seen = norm_f;
            stall_count = 0;
        } else {
            stall_count += 1;
            if stall_count >= newton.max_stall_iter {
                eprintln!("Newton: stalled for {} iterations, accepting ‖F‖={norm_f:.3e}",
                          newton.max_stall_iter);
                return Ok(());
            }
        }
    }

    bail!("Newton solver did not converge in {} iterations", newton.max_iter)
}

/// Returns (alpha, improved) where improved=true iff ‖F(x+α·step)‖ < ‖F(x)‖.
fn line_search(
    x: &[f64],
    step: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    norm_f0: f64,
) -> Result<(f64, bool)> {
    let n = x.len();
    let mut alpha = 1.0_f64;
    let mut f_new = vec![0.0_f64; n];
    let mut x_new = x.to_vec();
    let mut best_alpha = 0.0_f64;  // 0 = no step is safer than an unphysical step
    let mut best_norm = norm_f0;   // initialise to current norm, not INFINITY

    for _ in 0..12 {
        for i in 0..n {
            x_new[i] = x[i] + alpha * step[i];
        }
        // Clamp before evaluation to prevent NaN from unphysical temperatures.
        clamp_physical(&mut x_new, mech, grid);
        eval_residual(&x_new, &mut f_new, mech, grid, config, None, 0.0);
        let norm_new = norm2(&f_new);
        if norm_new.is_finite() && norm_new < norm_f0 {
            return Ok((alpha, true)); // first improvement found
        }
        if norm_new.is_finite() && norm_new < best_norm {
            best_norm = norm_new;
            best_alpha = alpha;
        }
        alpha *= 0.5;
    }
    // If best_alpha == 0 here, no step improved the residual, so we stay put.
    Ok((best_alpha, false))
}

/// Returns the largest factor ≤ 1 such that no component of the clipped step
/// violates physical scale limits:
///   T: max change 500 K
///   Y: max change 0.5
///   M: max change 2× |M| (but at least 1 kg/(m²·s))
fn step_clip_factor(x: &[f64], step: &[f64], mech: &Mechanism, grid: &Grid) -> f64 {
    let nk = mech.n_species();
    let nv = natj(mech);
    let nj = grid.n_points();
    let mut factor = 1.0_f64;

    for j in 0..nj {
        let dt = step[idx_t(nv, j)].abs();
        if dt > 500.0 { factor = factor.min(500.0 / dt); }
        for k in 0..nk {
            let dy = step[idx_y(nv, j, k)].abs();
            if dy > 0.5 { factor = factor.min(0.5 / dy); }
        }
    }
    let im = idx_m(nv, nj);
    let dm = step[im].abs();
    let m_scale = x[im].abs().max(1.0) * 2.0;
    if dm > m_scale { factor = factor.min(m_scale / dm); }

    factor
}

/// Clip variables to physically valid ranges to prevent NaN propagation.
fn clamp_physical(x: &mut Vec<f64>, mech: &Mechanism, grid: &Grid) {
    let nk = mech.n_species();
    let nv = natj(mech);
    let nj = grid.n_points();

    for j in 0..nj {
        let it = idx_t(nv, j);
        if !x[it].is_finite() || x[it] < 1.0 { x[it] = 1.0; }
        if x[it] > 10000.0 { x[it] = 10000.0; }
        for k in 0..nk {
            let iy = idx_y(nv, j, k);
            // Clip to physical mass-fraction range [0, 1].  Negative or > 1
            // values cause W_mean → 0 or density overflow → NaN.
            if !x[iy].is_finite() || x[iy] < 0.0 { x[iy] = 0.0; }
            if x[iy] > 1.0 { x[iy] = 1.0; }
        }
    }
    let im = idx_m(nv, nj);
    if !x[im].is_finite() || x[im] <= 0.0 { x[im] = 1e-6; }
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;
    use crate::flame::domain::Grid;
    use crate::flame::residual::FlameConfig;
    use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};

    fn h2o2_mech() -> Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    // The uniform N2 profile is an exact steady-state solution (‖F‖ ≈ 0).
    // The scaled convergence criterion should fire at iteration 0, before
    // any Jacobian is ever built (avoiding the M-column degeneracy).
    #[test]
    fn test_newton_converges_immediately_at_zero_residual() {
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

        // Very loose tolerances — ‖F‖ ≈ 0 should satisfy them at iter 0.
        let newton_cfg = NewtonConfig {
            atol: 1e-6,
            rtol: 1e-4,
            max_iter: 1,
            max_jac_age: 5,
            max_stall_iter: 5,
        };
        solve(&mut x, &mech, &grid, &config, &newton_cfg)
            .expect("Newton should satisfy convergence at iter 0 for zero-residual solution");
    }

    // Verify the scaled convergence criterion: ‖F‖ < atol·√n + rtol·‖x‖.
    // With ‖F‖ = 0 and any positive tolerances this must hold.
    #[test]
    fn test_newton_scaled_convergence_criterion() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 4;
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

        // Evaluate residual directly and check against the criterion.
        let n = solution_length(&mech, nj);
        let mut f = vec![0.0_f64; n];
        eval_residual(&x, &mut f, &mech, &grid, &config, None, 0.0);
        let norm_f = norm2(&f);
        let threshold = 1e-6 * (n as f64).sqrt() + 1e-4 * norm2(&x);

        assert!(
            norm_f < threshold,
            "‖F‖ = {norm_f:.3e} should be < threshold {threshold:.3e} for exact N2 profile"
        );
    }
}
