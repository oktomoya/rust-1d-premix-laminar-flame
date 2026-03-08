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
    /// Max number of Jacobian reuses before re-evaluation
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
    let mut jac_age = newton.max_jac_age + 1; // force Jacobian evaluation on first iteration

    for iter in 0..newton.max_iter {
        eval_residual(x, &mut f, mech, grid, config);

        let norm_f = norm2(&f);
        eprintln!("Newton iter {iter:3}: ‖F‖ = {norm_f:.3e}");

        if norm_f < newton.atol + newton.rtol * norm2(x) {
            eprintln!("Newton converged in {iter} iterations");
            return Ok(());
        }

        // Build or reuse Jacobian
        if jac_age >= newton.max_jac_age {
            // (Re)compute Jacobian
            let mut jac = numerical_jacobian(x, &f, mech, grid, config);
            let mut step = f.clone(); // solve J·Δx = -F
            step.iter_mut().for_each(|v| *v = -*v);
            jac.solve(&mut step)?;

            // Backtracking line search
            let alpha = line_search(x, &f, &step, mech, grid, config, norm_f)?;
            for i in 0..n {
                x[i] += alpha * step[i];
            }
            jac_age = 0;
        } else {
            // Would reuse old Jacobian here; for simplicity always recompute
            jac_age = newton.max_jac_age + 1;
        }

        jac_age += 1;
    }

    bail!("Newton solver did not converge in {} iterations", newton.max_iter)
}

fn line_search(
    x: &[f64],
    f0: &[f64],
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
        eval_residual(&x_new, &mut f_new, mech, grid, config);
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
