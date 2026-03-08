/// Pseudo-transient continuation (time-stepping) to get Newton in the basin of attraction.
///
/// Adds a transient term (x - x_old)/dt to the residual:
///   F_transient(x) = F_steady(x) + (x - x_old) / dt
/// Solving this converges x toward the steady-state solution as dt → ∞.

use anyhow::Result;
use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::jacobian::numerical_jacobian;
use crate::flame::residual::{eval_residual, FlameConfig};
use crate::flame::state::solution_length;

pub struct PseudoTransientConfig {
    pub n_steps: usize,
    pub dt_initial: f64,
    pub dt_max: f64,
    /// Grow dt each step by this factor if the step succeeded
    pub dt_grow: f64,
}

impl Default for PseudoTransientConfig {
    fn default() -> Self {
        PseudoTransientConfig {
            n_steps: 100,
            dt_initial: 1e-7,
            dt_max: 1e-3,
            dt_grow: 1.5,
        }
    }
}

/// Apply pseudo-transient time stepping to x (modified in-place).
pub fn step(
    x: &mut Vec<f64>,
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    pt: &PseudoTransientConfig,
) -> Result<()> {
    let n = solution_length(mech, grid.n_points());
    let mut f = vec![0.0_f64; n];
    let mut dt = pt.dt_initial;

    for step_idx in 0..pt.n_steps {
        eval_residual(x, &mut f, mech, grid, config);
        let norm_f = norm2(&f);

        eprintln!("PT step {step_idx:4}: ‖F‖ = {norm_f:.3e},  dt = {dt:.2e}");

        if norm_f < 1e-3 {
            eprintln!("Pseudo-transient: switching to Newton (‖F‖ < 1e-3)");
            break;
        }

        // Build Jacobian of transient system: J_pt = J_steady + I/dt
        let mut jac = numerical_jacobian(x, &f, mech, grid, config);
        for i in 0..n {
            let val = jac.get(i, i) + 1.0 / dt;
            jac.set(i, i, val);
        }

        // RHS = -F_steady  (transient term: (x - x_old)/dt = 0 initially)
        let mut rhs: Vec<f64> = f.iter().map(|v| -v).collect();
        jac.solve(&mut rhs)?;

        for i in 0..n {
            x[i] += rhs[i];
        }

        dt = (dt * pt.dt_grow).min(pt.dt_max);
    }

    Ok(())
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
