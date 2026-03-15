/// Finite-difference Jacobian computation.
///
/// The Jacobian J = ∂F/∂x is banded due to the local stencil (each residual
/// F_j depends only on x_{j-1}, x_j, x_{j+1}).
/// Bandwidth: kl = ku = NATJ  (one block in each direction).
/// Total band width: 2*NATJ + 1 diagonals.
///
/// Storage: banded column-major (LAPACK dgbsv format).

use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::residual::{eval_residual, FlameConfig};
use crate::flame::state::solution_length;
use crate::solver::banded::BandedMatrix;
use crate::solver::dense::DenseMatrix;

pub fn numerical_jacobian(
    x: &[f64],
    f0: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) -> BandedMatrix {
    let n = solution_length(mech, grid.n_points());
    let nv = mech.n_species() + 1; // variables per grid point (NATJ)
    // True bandwidth: F at j depends on x at j-1, j, j+1.
    // Worst-case row-col distance for adjacent blocks is 2*nv - 1.
    // Adding 1 extra for safety and the global M variable (at index n-1)
    // which couples all equations but whose column sits far from most rows.
    // For a proper treatment of M we use kl = ku = 2*nv so that the
    // M column is within band for the last ~2*nv grid points.  The
    // upper rows' coupling to M is small compared to the I/dt diagonal
    // in the pseudo-transient phase; Newton convergence is handled by
    // the grid being refined near the flame where M matters most.
    let kl = 2 * nv;
    let ku = 2 * nv;

    let mut jac = BandedMatrix::new(n, kl, ku);
    let mut x_pert = x.to_vec();
    let mut f_pert = vec![0.0_f64; n];

    // Perturbation: h = eps_rel * max(|x_j|, 1.0).
    // Using max(|x_j|, 1.0) (not eps_rel) as the absolute scale prevents
    // h = eps_rel² ≈ 4e-16 when x_j = 0, which would cause catastrophic
    // cancellation in the divided difference for zero variables (radical species).
    let eps_rel = (2.0_f64 * f64::EPSILON).sqrt();

    for j in 0..n {
        let x_j = x[j];
        let h = eps_rel * x_j.abs().max(1.0);
        x_pert[j] = x_j + h;
        eval_residual(&x_pert, &mut f_pert, mech, grid, config, None, 0.0);

        // Only fill within the band
        let row_min = j.saturating_sub(ku);
        let row_max = (j + kl + 1).min(n);
        for i in row_min..row_max {
            let val = (f_pert[i] - f0[i]) / h;
            jac.set(i, j, val);
        }

        x_pert[j] = x_j; // restore
    }

    jac
}

/// Dense numerical Jacobian J = ∂F/∂x (all n² entries).
///
/// Unlike the banded version, this captures the full coupling including the
/// global mass-flux eigenvalue M at index n-1 which couples every grid-point
/// residual.  Required for Newton and pseudo-transient solves where the banded
/// approximation drops the M column for all but the last few grid points.
pub fn numerical_jacobian_dense(
    x: &[f64],
    f0: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) -> DenseMatrix {
    let n = solution_length(mech, grid.n_points());
    let mut jac = DenseMatrix::new(n);
    let mut x_pert = x.to_vec();
    let mut f_pert = vec![0.0_f64; n];

    let eps_rel = (2.0_f64 * f64::EPSILON).sqrt();

    for j in 0..n {
        let x_j = x[j];
        let h = eps_rel * x_j.abs().max(1.0);
        x_pert[j] = x_j + h;
        eval_residual(&x_pert, &mut f_pert, mech, grid, config, None, 0.0);

        for i in 0..n {
            jac.set(i, j, (f_pert[i] - f0[i]) / h);
        }

        x_pert[j] = x_j;
    }

    jac
}
