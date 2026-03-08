/// Adaptive grid refinement based on GRAD and CURV criteria (as in PREMIX).
///
/// A new point is inserted between j and j+1 if:
///   GRAD: |φ_{j+1} - φ_j| > grad * (max|φ| - min|φ|) + eps
///   CURV: |κ_j|             > curv * max|dφ/dz| + eps
/// where κ_j is the second derivative at the midpoint.

use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;
use crate::flame::state::{idx_t, idx_y, natj, solution_length};

pub struct RefineCriteria {
    /// Gradient criterion (0–1; smaller = finer grid)
    pub grad: f64,
    /// Curvature criterion (0–1; smaller = finer grid)
    pub curv: f64,
    /// Maximum allowed grid points
    pub max_points: usize,
}

impl Default for RefineCriteria {
    fn default() -> Self {
        RefineCriteria { grad: 0.05, curv: 0.10, max_points: 500 }
    }
}

/// Returns a list of interval indices (j) where a new midpoint should be inserted.
/// Caller is responsible for inserting in reverse order to preserve indices.
pub fn find_refinement_points(
    x: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    criteria: &RefineCriteria,
) -> Vec<usize> {
    let nk = mech.n_species();
    let nv = natj(mech);
    let nj = grid.n_points();
    let dz = grid.dz();

    // Collect the field variables we check (T and all species)
    // phi[var][j]
    let mut phi: Vec<Vec<f64>> = Vec::new();
    // Temperature
    phi.push((0..nj).map(|j| x[idx_t(nv, j)]).collect());
    // Species
    for k in 0..nk {
        phi.push((0..nj).map(|j| x[idx_y(nv, j, k)]).collect());
    }

    let mut insert = vec![false; nj - 1];

    for var in &phi {
        let phi_min = var.iter().cloned().fold(f64::INFINITY, f64::min);
        let phi_max = var.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let phi_range = (phi_max - phi_min).max(1e-30);

        // First derivatives at midpoints
        let dphi_dz: Vec<f64> = (0..nj - 1)
            .map(|j| (var[j + 1] - var[j]) / dz[j])
            .collect();
        let max_dphi = dphi_dz.iter().cloned().map(f64::abs).fold(0.0_f64, f64::max).max(1e-30);

        for j in 0..nj - 1 {
            // GRAD criterion
            if (var[j + 1] - var[j]).abs() > criteria.grad * phi_range {
                insert[j] = true;
            }
            // CURV criterion (second derivative estimate at j+1, using adjacent midpoints)
            if j + 1 < nj - 1 {
                let curv = (dphi_dz[j + 1] - dphi_dz[j]).abs()
                    / (0.5 * (dz[j] + dz[j + 1]));
                if curv > criteria.curv * max_dphi {
                    insert[j] = true;
                    insert[j + 1] = true;
                }
            }
        }
    }

    insert.iter().enumerate()
        .filter(|(_, &b)| b)
        .map(|(j, _)| j)
        .collect()
}

/// Refine the grid, inserting midpoints and linearly interpolating the solution.
/// Returns the new grid and the interpolated solution vector.
pub fn refine(
    x: &[f64],
    mech: &Mechanism,
    grid: &Grid,
    criteria: &RefineCriteria,
) -> Option<(Grid, Vec<f64>)> {
    let nv = natj(mech);
    let nj = grid.n_points();

    if nj >= criteria.max_points {
        return None;
    }

    let mut insert_at = find_refinement_points(x, mech, grid, criteria);
    if insert_at.is_empty() {
        return None;
    }

    // Remove duplicate indices and sort in reverse for safe insertion
    insert_at.sort_unstable();
    insert_at.dedup();

    // Limit total number of new points
    let max_new = criteria.max_points.saturating_sub(nj);
    if insert_at.len() > max_new {
        insert_at.truncate(max_new);
    }

    // Build new grid and solution
    let mut new_z = grid.z.clone();
    let mut new_x = x.to_vec();
    let m = x[idx_m_raw(nv, nj)]; // preserve mass flux

    // Insert in reverse order so earlier indices remain valid
    for &j in insert_at.iter().rev() {
        let z_new = 0.5 * (new_z[j] + new_z[j + 1]);
        new_z.insert(j + 1, z_new);

        // Linearly interpolate all variables at the new point
        let start = j * nv;
        let end   = (j + 1) * nv;
        let left: Vec<f64>  = new_x[start..start + nv].to_vec();
        let right: Vec<f64> = new_x[end..end + nv].to_vec();
        let mid: Vec<f64>   = left.iter().zip(right.iter())
            .map(|(l, r)| 0.5 * (l + r))
            .collect();
        new_x.splice((j + 1) * nv..(j + 1) * nv, mid);
    }

    // Fix mass flux location (appended at end)
    let new_nj = new_z.len();
    if new_x.len() == nv * new_nj {
        // Mass flux is already embedded; nothing to do
    } else {
        // Append mass flux if missing
        new_x.push(m);
    }

    let new_grid = Grid { z: new_z };
    Some((new_grid, new_x))
}

fn idx_m_raw(nv: usize, nj: usize) -> usize {
    nv * nj
}
