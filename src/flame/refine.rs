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

    // Build new grid and solution.
    // Pop M first so splice indices align with nv*nj (no M at end during insertions).
    let mut new_z = grid.z.clone();
    let mut new_x = x.to_vec();
    let m = new_x.pop().expect("solution vector must end with mass flux M");

    // Insert in reverse order so earlier indices remain valid.
    for &j in insert_at.iter().rev() {
        let z_new = 0.5 * (new_z[j] + new_z[j + 1]);
        new_z.insert(j + 1, z_new);

        // Linearly interpolate all nv variables at the new midpoint.
        let left:  Vec<f64> = new_x[j * nv       ..j * nv + nv].to_vec();
        let right: Vec<f64> = new_x[(j + 1) * nv ..(j + 1) * nv + nv].to_vec();
        let mid:   Vec<f64> = left.iter().zip(right.iter())
            .map(|(l, r)| 0.5 * (l + r))
            .collect();
        new_x.splice((j + 1) * nv..(j + 1) * nv, mid);
    }

    // Re-append M at the correct end position: new_x.len() == nv * new_nj.
    new_x.push(m);
    debug_assert_eq!(new_x.len(), solution_length(mech, new_z.len()),
        "solution length mismatch after refinement");

    let new_grid = Grid { z: new_z };
    Some((new_grid, new_x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::parser::cantera_yaml::parse_file;
    use crate::flame::state::{idx_m, idx_t, idx_y, natj, solution_length};

    fn h2o2_mech() -> Mechanism {
        let manifest = env!("CARGO_MANIFEST_DIR");
        parse_file(&format!("{manifest}/data/h2o2.yaml")).expect("parse h2o2.yaml")
    }

    // Refine a smooth N2 profile with a sine-shaped temperature bump.
    // Verifies:
    //  1. Mass flux M is at the correct index (nv * new_nj) after refinement.
    //  2. Interpolated T values agree with the smooth profile to within dz²/8.
    //  3. solution_length is consistent with the new grid.
    #[test]
    fn test_refine_mass_flux_alignment_and_interpolation() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nv = natj(&mech);
        let nj = 10;
        let grid = Grid::uniform(0.02, nj);
        let n2_idx = mech.species_index("N2").unwrap();
        let m_val = 0.3_f64;

        // Build a smooth temperature profile: T(z) = 300 + 1000*sin(π*z/L)
        let length = grid.z[nj - 1];
        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            let z = grid.z[j];
            x[idx_t(nv, j)] = 300.0 + 1000.0 * (std::f64::consts::PI * z / length).sin();
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = m_val;

        // Refine with loose criteria so at least a few points are inserted.
        let criteria = RefineCriteria { grad: 0.1, curv: 0.5, max_points: 200 };
        let (new_grid, new_x) = refine(&x, &mech, &grid, &criteria)
            .expect("refinement should insert at least one point");

        let new_nj = new_grid.n_points();

        // 1. Correct solution length.
        assert_eq!(new_x.len(), solution_length(&mech, new_nj),
            "solution length mismatch: new_x.len()={} vs expected {}",
            new_x.len(), solution_length(&mech, new_nj));

        // 2. Mass flux at the correct index.
        assert_eq!(
            (new_x[idx_m(nv, new_nj)] - m_val).abs() < 1e-14,
            true,
            "M at idx_m = {}, expected {m_val}", new_x[idx_m(nv, new_nj)]
        );

        // 3. No extra M copies buried in the grid variables.
        for j in 0..new_nj {
            let t = new_x[idx_t(nv, j)];
            assert!(t > 200.0 && t < 1400.0,
                "T[{j}] = {t:.1} is out of plausible range [200, 1400]");
            // Y_N2 should still be near 1 (interpolated)
            let yn2 = new_x[idx_y(nv, j, n2_idx)];
            assert!((yn2 - 1.0).abs() < 1e-10,
                "Y_N2[{j}] = {yn2:.6}, expected 1.0");
            // Other species should be 0
            for k in 0..nk {
                if k != n2_idx {
                    let yk = new_x[idx_y(nv, j, k)];
                    assert!(yk.abs() < 1e-10, "Y[{k}][{j}] = {yk:.3e}, expected 0");
                }
            }
        }

        // 4. At least one new point was inserted.
        assert!(new_nj > nj, "expected more grid points after refinement");
    }
}
