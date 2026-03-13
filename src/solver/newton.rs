/// Damped Newton iteration with backtracking line search.

use anyhow::{bail, Result};
use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::jacobian::numerical_jacobian;
use crate::flame::residual::{eval_residual, FlameConfig};
use crate::flame::state::solution_length;

pub struct NewtonConfig {
    pub atol: f64,
    pub rtol: f64,
    pub max_iter: usize,
    /// Max number of iterations before the Jacobian is recomputed.
    pub max_jac_age: usize,
}

impl Default for NewtonConfig {
    fn default() -> Self {
        NewtonConfig {
            atol: 1e-9,
            rtol: 1e-6,
            max_iter: 50,
            max_jac_age: 5,
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
    // Cached factored Jacobian; reused for up to max_jac_age iterations.
    let mut cached_jac: Option<crate::solver::banded::BandedMatrix> = None;
    let mut jac_age = newton.max_jac_age + 1; // force evaluation on first iteration

    for iter in 0..newton.max_iter {
        eval_residual(x, &mut f, mech, grid, config, None, 0.0);

        let norm_f = norm2(&f);
        eprintln!("Newton iter {iter:3}: ‖F‖ = {norm_f:.3e}");

        // Scaled convergence: ‖F‖ < atol·√n + rtol·‖x‖
        if norm_f < newton.atol * (n as f64).sqrt() + newton.rtol * norm2(x) {
            eprintln!("Newton converged in {iter} iterations");
            return Ok(());
        }

        // (Re)compute and factor Jacobian when stale or absent.
        if jac_age >= newton.max_jac_age || cached_jac.is_none() {
            let mut jac = numerical_jacobian(x, &f, mech, grid, config);
            jac.factor_in_place()?;
            cached_jac = Some(jac);
            jac_age = 0;
        }

        // Solve cached J·step = -F using the stored factorization.
        let mut step: Vec<f64> = f.iter().map(|v| -v).collect();
        cached_jac.as_ref().unwrap().solve_factored(&mut step);

        // Backtracking line search.
        let alpha = line_search(x, &f, &step, mech, grid, config, norm_f)?;
        for i in 0..n {
            x[i] += alpha * step[i];
        }

        // Invalidate cached Jacobian if line search had to cut the step heavily.
        if alpha < 0.1 {
            jac_age = newton.max_jac_age + 1;
        } else {
            jac_age += 1;
        }
    }

    bail!("Newton solver did not converge in {} iterations", newton.max_iter)
}

fn line_search(
    x: &[f64],
    _f0: &[f64],
    step: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    norm_f0: f64,
) -> Result<f64> {
    let n = x.len();
    let mut alpha = 1.0_f64;
    let mut f_new = vec![0.0_f64; n];
    let mut x_new = x.to_vec();

    for _ in 0..10 {
        for i in 0..n {
            x_new[i] = x[i] + alpha * step[i];
        }
        eval_residual(&x_new, &mut f_new, mech, grid, config, None, 0.0);
        if norm2(&f_new) < norm_f0 {
            return Ok(alpha);
        }
        alpha *= 0.5;
    }
    Ok(alpha) // accept small step rather than failing
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
