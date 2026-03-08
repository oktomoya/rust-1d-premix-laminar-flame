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

pub fn numerical_jacobian(
    x: &[f64],
    f0: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) -> BandedMatrix {
    let n = solution_length(mech, grid.n_points());
    let nv = mech.n_species() + 1; // variables per grid point (NATJ)
    let kl = nv;     // lower bandwidth
    let ku = nv;     // upper bandwidth

    let mut jac = BandedMatrix::new(n, kl, ku);
    let mut x_pert = x.to_vec();
    let mut f_pert = vec![0.0_f64; n];

    // Relative + absolute perturbation (SQRT of machine epsilon)
    let eps_rel = (2.0_f64 * f64::EPSILON).sqrt();
    let eps_abs = eps_rel;

    for j in 0..n {
        let x_j = x[j];
        let h = eps_rel * x_j.abs().max(eps_abs);
        x_pert[j] = x_j + h;
        eval_residual(&x_pert, &mut f_pert, mech, grid, config);

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
