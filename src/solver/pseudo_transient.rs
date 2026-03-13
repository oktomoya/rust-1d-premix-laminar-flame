/// Pseudo-transient continuation (time-stepping) to get Newton in the basin of attraction.
///
/// At each step, solves the implicit system:
///   (J_steady + I/dt) · Δx = -F_steady(x_old)
/// where x_old is the state at the beginning of the step.  This is equivalent to
/// backward-Euler integration of dx/dt = -F_steady(x) starting from x_old.

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
        // Save x_old explicitly — the transient term is (x - x_old)/dt.
        // Evaluating F_pt at x = x_old gives F_steady(x_old) (PT term = 0),
        // which is the correct RHS for the backward-Euler linear system.
        let x_old = x.clone();
        let rdt = 1.0 / dt;

        eval_residual(x, &mut f, mech, grid, config, Some(&x_old), rdt);
        let norm_f = norm2(&f);

        eprintln!("PT step {step_idx:4}: ‖F‖ = {norm_f:.3e},  dt = {dt:.2e}");

        if norm_f < 1e-3 {
            eprintln!("Pseudo-transient: switching to Newton (‖F‖ < 1e-3)");
            break;
        }

        // Jacobian of PT system: J_pt = J_steady + I/dt.
        // numerical_jacobian uses F_steady, so we add I/dt manually.
        let mut jac = numerical_jacobian(x, &f, mech, grid, config);
        for i in 0..n {
            let val = jac.get(i, i) + rdt;
            jac.set(i, i, val);
        }

        // RHS = -F_pt(x_old) = -F_steady(x_old)  (PT term is zero at x = x_old).
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

    // Starting from a perturbed N2 profile, PT steps should reduce the
    // steady-state residual norm toward zero.
    // (The +I/dt diagonal prevents zero-pivot even when J_steady is singular.)
    #[test]
    fn test_pseudo_transient_reduces_residual() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nj = 6;
        let nv = natj(&mech);
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();

        // Uniform N2 profile — exact steady state.
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

        // Perturb one interior temperature node.
        x[idx_t(nv, nj / 2)] += 50.0;

        // Measure initial steady-state residual.
        let n = solution_length(&mech, nj);
        let mut f_init = vec![0.0_f64; n];
        eval_residual(&x, &mut f_init, &mech, &grid, &config, None, 0.0);
        let norm_init = norm2(&f_init);

        // Run a single PT step (large dt to get near-Newton behaviour).
        // One step is sufficient to verify the residual decreases; multi-step
        // behaviour depends on the problem conditioning.
        let pt_cfg = PseudoTransientConfig {
            n_steps: 1,
            dt_initial: 1e-4,
            ..Default::default()
        };
        step(&mut x, &mech, &grid, &config, &pt_cfg).expect("PT must not fail");

        // Measure final steady-state residual.
        let mut f_final = vec![0.0_f64; n];
        eval_residual(&x, &mut f_final, &mech, &grid, &config, None, 0.0);
        let norm_final = norm2(&f_final);

        assert!(
            norm_final < norm_init,
            "PT should reduce residual: initial {norm_init:.3e}, final {norm_final:.3e}"
        );
    }

    // When x_old = x, F_pt(x) = F_steady(x) (transient term is zero).
    // This is the same property verified in residual.rs but tested here
    // through the PT solver's explicit x_old tracking.
    #[test]
    fn test_pseudo_transient_preserves_steady_solution() {
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

        // PT on the exact steady state: ‖F‖ < 1e-3 at step 0 → exits immediately.
        // x should not change (step = 0 because RHS ≈ 0).
        let x_before = x.clone();
        let pt_cfg = PseudoTransientConfig {
            n_steps: 10,
            dt_initial: 1e-6,
            ..Default::default()
        };
        step(&mut x, &mech, &grid, &config, &pt_cfg).expect("PT must not fail");

        // x should be essentially unchanged.
        let max_delta = x.iter().zip(x_before.iter()).map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_delta < 1e-10,
            "PT on exact steady state should not change x (max Δ = {max_delta:.3e})"
        );
    }
}
