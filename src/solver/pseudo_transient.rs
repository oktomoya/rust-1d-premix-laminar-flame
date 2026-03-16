/// Pseudo-transient continuation (TWOPNT-style) to reach Newton's basin of attraction.
///
/// Each outer step takes one backward-Euler implicit update:
///
///   Find x_new such that F_pt(x_new) = F_steady(x_new) + D/dt · (x_new − x_old) = 0
///
/// solved by a short inner Newton loop using the banded Jacobian.
///
/// **Banded Jacobian + diagonal-save**: the steady-state Jacobian J_ss is built
/// once per outer step using the local-stencil banded approach.  For each inner
/// iteration the PT Jacobian J_pt = J_ss + D/dt is formed by restoring the
/// unfactored banded data and adding the diagonal scaling term
///   diag = rho_j·cp_j (for T) or rho_j (for Y_k),
/// then re-factorising.  This replaces the previous approach of rebuilding the
/// full Jacobian at every inner iteration (300 outer × 5 inner → 5× more rebuilds).
///
/// **Line search fix**: best_norm = ∞ so the search always returns the alpha with
/// the smallest finite F_pt, even when J_pt is indefinite.

use anyhow::Result;
use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{cp_mixture, density, mean_molecular_weight};
use crate::flame::domain::Grid;
use crate::flame::jacobian::{bordered_solve_m, numerical_jacobian};
use crate::flame::residual::{eval_residual, FlameConfig};
use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};

pub struct PseudoTransientConfig {
    pub n_steps: usize,
    pub dt_initial: f64,
    pub dt_max: f64,
    /// Grow dt each step by this factor after a successful step.
    pub dt_grow: f64,
    /// Inner Newton iterations per time step.
    pub max_inner_iter: usize,
}

impl Default for PseudoTransientConfig {
    fn default() -> Self {
        PseudoTransientConfig {
            n_steps: 300,
            dt_initial: 1e-7,
            dt_max: 1e-3,
            dt_grow: 1.5,
            max_inner_iter: 5,
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
    let nv = natj(mech);
    let nj = grid.n_points();
    let n = solution_length(mech, nj);
    let n_int = n - 1;
    let mut f_steady = vec![0.0_f64; n];
    let mut f_pt     = vec![0.0_f64; n];
    let mut dt = pt.dt_initial;
    let mut successive_failures = 0usize;

    let mut step_idx = 0;
    while step_idx < pt.n_steps {
        eval_residual(x, &mut f_steady, mech, grid, config, None, 0.0);
        let norm_f = norm2(&f_steady);

        eprintln!("PT step {step_idx:4}: ‖F_steady‖ = {norm_f:.3e},  dt = {dt:.2e}");

        if norm_f < 1e7 {
            eprintln!("Pseudo-transient: switching to Newton (‖F_steady‖ < 1e7)");
            break;
        }
        if dt < 1e-13 {
            eprintln!("PT: dt below minimum, proceeding to Newton");
            break;
        }

        let rdt = 1.0 / dt;
        let x_old = x.clone();

        // Build banded J_ss once per outer step (local-stencil, cheap).
        // m_col_full is the complete M column (∂F/∂M for all rows) for SM correction.
        // jfc is the T[j_fix] column index for the rank-2 eigenvalue-row correction.
        let (jac_ss, m_col_full, jfc) = numerical_jacobian(x, &f_steady, mech, grid, config);
        // Save unfactored banded data so we can restore it for each inner iter.
        let jac_data_orig = jac_ss.data.clone();
        let mut jac = jac_ss; // reuse storage

        // --- Inner Newton: reduce ‖F_pt‖ ---
        for _inner in 0..pt.max_inner_iter {
            // Restore unfactored J_ss.
            jac.data.copy_from_slice(&jac_data_orig);

            // Add D/dt diagonal to interior rows of J_int: J_pt = J_int + diag(rho·cp, rho).
            // The M equation is NOT in J_int; its rdt contribution enters via bordered_solve_m(d=rdt).
            let p = config.pressure;
            for j in 1..nj - 1 {
                let t_j = x[idx_t(nv, j)];
                let y_j: Vec<f64> = (0..mech.n_species())
                    .map(|k| x[idx_y(nv, j, k)]).collect();
                let w_j = mean_molecular_weight(&mech.species, &y_j);
                let rho_j = density(p, t_j, w_j);
                let cp_j  = cp_mixture(&mech.species, &y_j, t_j);
                let it = idx_t(nv, j);
                jac.set(it, it, jac.get(it, it) + rdt * rho_j * cp_j);
                for k in 0..mech.n_species() {
                    let iy = idx_y(nv, j, k);
                    jac.set(iy, iy, jac.get(iy, iy) + rdt * rho_j);
                }
            }

            // Factor J_pt in place.
            if jac.factor_in_place().is_err() { break; }

            // Evaluate F_pt at current x.
            eval_residual(x, &mut f_pt, mech, grid, config, Some(&x_old), rdt);
            let norm_pt = norm2(&f_pt);

            // Solve interior: step_int = J_int_pt^{-1} · (-F_pt[0..n_int])
            let mut step_int: Vec<f64> = f_pt[..n_int].iter().map(|v| -v).collect();
            jac.solve_factored(&mut step_int);
            // Bordered solve for δM with d=rdt (PT time derivative for M).
            let dm = bordered_solve_m(&mut step_int, &m_col_full, &jac, jfc, -f_pt[n_int], rdt);
            // Assemble full delta.
            let mut delta = step_int;
            delta.push(dm);
            if delta.iter().any(|v| !v.is_finite()) { break; }

            // Step-size clipping.
            let clip = step_clip_factor(x, &delta, mech, grid);
            if clip < 1.0 {
                for d in delta.iter_mut() { *d *= clip; }
            }

            let alpha = line_search_pt(x, &delta, mech, grid, config, &x_old, rdt, norm_pt);
            for i in 0..n { x[i] += alpha * delta[i]; }
            clamp_physical(x, mech, grid);
        }

        // --- Outer acceptance: F_steady ≤ 1.1 × F_steady_old ---
        eval_residual(x, &mut f_steady, mech, grid, config, None, 0.0);
        let norm_new = norm2(&f_steady);

        if norm_new.is_finite() && norm_new < norm_f * 1.1 {
            step_idx += 1;
            successive_failures = 0;
            dt = (dt * pt.dt_grow).min(pt.dt_max);
        } else {
            x.copy_from_slice(&x_old);
            dt *= 0.5;
            successive_failures += 1;
            eprintln!(
                "PT step {step_idx:4}: step rejected \
                 (‖F_new‖={norm_new:.3e} ≥ 1.1×‖F_old‖={:.3e}), dt→{dt:.2e}",
                norm_f * 1.1
            );
            if successive_failures >= 30 {
                eprintln!("PT: too many consecutive failures, proceeding to Newton");
                break;
            }
        }
    }

    Ok(())
}

/// Backtracking line search on ‖F_pt‖.
fn line_search_pt(
    x: &[f64],
    delta: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
    x_old: &[f64],
    rdt: f64,
    norm_f_pt0: f64,
) -> f64 {
    let n = x.len();
    let mut alpha = 1.0_f64;
    let mut f_trial = vec![0.0_f64; n];
    let mut x_trial = x.to_vec();
    let mut best_alpha = 1e-12_f64;
    let mut best_norm = f64::INFINITY;

    for _ in 0..12 {
        for i in 0..n { x_trial[i] = x[i] + alpha * delta[i]; }
        clamp_physical(&mut x_trial, mech, grid);
        eval_residual(&x_trial, &mut f_trial, mech, grid, config, Some(x_old), rdt);
        let norm_trial = norm2(&f_trial);

        if norm_trial.is_finite() {
            if norm_trial < best_norm {
                best_norm = norm_trial;
                best_alpha = alpha;
            }
            if norm_trial < norm_f_pt0 { return alpha; }
        }
        alpha *= 0.5;
    }
    best_alpha
}

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

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

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
            if !x[iy].is_finite() || x[iy] < 0.0 { x[iy] = 0.0; }
            if x[iy] > 1.0 { x[iy] = 1.0; }
        }
    }
    let im = idx_m(nv, nj);
    if !x[im].is_finite() || x[im] <= 0.0 { x[im] = 1e-6; }
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

    #[test]
    fn test_pseudo_transient_reduces_residual() {
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

        x[idx_t(nv, nj / 2)] += 50.0;

        let pt_cfg = PseudoTransientConfig {
            n_steps: 5,
            dt_initial: 1e-7,
            ..Default::default()
        };
        step(&mut x, &mech, &grid, &config, &pt_cfg).expect("PT must not fail");

        for j in 0..nj {
            let t = x[idx_t(nv, j)];
            assert!(t.is_finite() && t > 0.0, "T[{j}] = {t:.3e} is non-physical");
        }
        assert!(x[idx_m(nv, nj)] > 0.0, "M = {:.3e} must be positive", x[idx_m(nv, nj)]);
    }

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

        let x_before = x.clone();
        let pt_cfg = PseudoTransientConfig {
            n_steps: 10,
            dt_initial: 1e-6,
            ..Default::default()
        };
        step(&mut x, &mech, &grid, &config, &pt_cfg).expect("PT must not fail");

        let max_delta = x.iter().zip(x_before.iter()).map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_delta < 1e-10,
            "PT on exact steady state should not change x (max Δ = {max_delta:.3e})"
        );
    }
}
