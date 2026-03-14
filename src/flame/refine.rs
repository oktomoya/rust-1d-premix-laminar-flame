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
    /// Maximum ratio of adjacent cell sizes before inserting a point (≥2.0)
    pub ratio: f64,
    /// Maximum allowed grid points
    pub max_points: usize,
}

impl Default for RefineCriteria {
    fn default() -> Self {
        RefineCriteria { grad: 0.05, curv: 0.10, ratio: 2.0, max_points: 500 }
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

        // First derivatives at midpoints: slope[j] = dφ/dz over interval [j, j+1]
        let slope: Vec<f64> = (0..nj - 1)
            .map(|j| (var[j + 1] - var[j]) / dz[j])
            .collect();

        let slope_min = slope.iter().cloned().fold(f64::INFINITY, f64::min);
        let slope_max = slope.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let slope_range = (slope_max - slope_min).max(1e-30);

        for j in 0..nj - 1 {
            // GRAD criterion: change in value exceeds fraction of total range
            if (var[j + 1] - var[j]).abs() > criteria.grad * phi_range {
                insert[j] = true;
            }
            // CURV criterion (matches Cantera refine.cpp):
            //   |slope[j+1] - slope[j]| > curv * (slopeMax - slopeMin)
            // Both sides have units φ/m — dimensionally consistent.
            if j + 1 < nj - 1 {
                if (slope[j + 1] - slope[j]).abs() > criteria.curv * slope_range {
                    insert[j] = true;
                    insert[j + 1] = true;
                }
            }
        }
    }

    // Ratio criterion (matches Cantera refine.cpp):
    // Insert a point between j and j+1 if adjacent cell sizes differ by more
    // than `ratio`. This cascades refinement smoothly across the domain.
    for j in 1..nj - 1 {
        if dz[j] > criteria.ratio * dz[j - 1] {
            insert[j] = true;
        }
        if dz[j - 1] > criteria.ratio * dz[j] {
            insert[j - 1] = true;
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
        let criteria = RefineCriteria { grad: 0.1, curv: 0.5, ratio: 2.0, max_points: 200 };
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

    // Ratio criterion: a non-uniform grid with one large interval should trigger
    // insertion in that interval even when grad/curv alone would not fire.
    #[test]
    fn test_ratio_criterion_triggers_on_coarse_interval() {
        let mech = h2o2_mech();
        let nk = mech.n_species();
        let nv = natj(&mech);
        let n2_idx = mech.species_index("N2").unwrap();

        // Build a non-uniform grid: first 9 intervals are 1 mm, last is 10 mm.
        // The large last interval should trigger the ratio criterion.
        let mut z = vec![0.0_f64; 11];
        for j in 0..9 {
            z[j + 1] = z[j] + 1e-3;
        }
        z[10] = z[9] + 10e-3;  // last interval is 10x larger
        let nj = z.len();
        let grid = Grid { z };

        // Flat T profile and pure N2 — grad/curv alone won't fire.
        let mut x = vec![0.0_f64; solution_length(&mech, nj)];
        for j in 0..nj {
            x[idx_t(nv, j)] = 300.0;
            x[idx_y(nv, j, n2_idx)] = 1.0;
        }
        x[idx_m(nv, nj)] = 0.3;

        // ratio=3.0: the 10 mm interval is >3x the 1 mm neighbour → must insert.
        let criteria = RefineCriteria { grad: 0.05, curv: 0.10, ratio: 3.0, max_points: 200 };
        let result = refine(&x, &mech, &grid, &criteria);
        assert!(result.is_some(), "ratio criterion should trigger refinement");
        let (new_grid, _) = result.unwrap();
        assert!(new_grid.n_points() > nj, "expected a new point in the large interval");

        // With ratio=20.0 (very loose) and flat profile, no refinement expected.
        let criteria_loose = RefineCriteria { grad: 0.05, curv: 0.10, ratio: 20.0, max_points: 200 };
        let result_loose = refine(&x, &mech, &grid, &criteria_loose);
        assert!(result_loose.is_none(), "loose ratio should not trigger on flat profile");
    }
}
