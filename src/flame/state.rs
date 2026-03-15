/// Solution state and flat vector layout for the Newton solver.
///
/// Solution vector layout (column-major by grid point):
///   x = [T_0, Y_{0,0}, ..., Y_{K-1,0},   T_1, Y_{0,1}, ...,   ...,   T_{J-1}, ..., M]
///   Length = (K+1)*J + 1
///   - Variables per point: NATJ = K + 1  (temperature + K species)
///   - Index of T at point j:          j * NATJ
///   - Index of Yk at point j:         j * NATJ + 1 + k
///   - Index of mass flux M:           NATJ * J   (last element)

use anyhow::Result;
use crate::chemistry::mechanism::Mechanism;
use crate::flame::domain::Grid;

/// Number of unknowns per grid point: T + K species = K+1.
pub fn natj(mech: &Mechanism) -> usize {
    mech.n_species() + 1
}

/// Total length of the solution vector.
pub fn solution_length(mech: &Mechanism, n_points: usize) -> usize {
    natj(mech) * n_points + 1
}

/// Index of temperature at grid point j.
#[inline]
pub fn idx_t(natj: usize, j: usize) -> usize {
    j * natj
}

/// Index of species k mass fraction at grid point j.
#[inline]
pub fn idx_y(natj: usize, j: usize, k: usize) -> usize {
    j * natj + 1 + k
}

/// Index of the mass flux M (eigenvalue).
#[inline]
pub fn idx_m(natj: usize, n_points: usize) -> usize {
    natj * n_points
}

/// High-level view of the solution (borrows from a flat Vec<f64>).
pub struct FlameState<'a> {
    pub x: &'a [f64],
    pub mech: &'a Mechanism,
    pub grid: &'a Grid,
    pub natj: usize,
}

impl<'a> FlameState<'a> {
    pub fn new(x: &'a [f64], mech: &'a Mechanism, grid: &'a Grid) -> Self {
        FlameState { x, mech, grid, natj: natj(mech) }
    }

    pub fn temperature(&self, j: usize) -> f64 {
        self.x[idx_t(self.natj, j)]
    }

    pub fn species(&self, k: usize, j: usize) -> f64 {
        self.x[idx_y(self.natj, j, k)]
    }

    pub fn mass_flux(&self) -> f64 {
        self.x[idx_m(self.natj, self.grid.n_points())]
    }

    /// Return mass fraction slice Y[0..K] at point j.
    pub fn y_slice(&self, j: usize) -> &[f64] {
        let start = idx_y(self.natj, j, 0);
        &self.x[start..start + self.mech.n_species()]
    }
}

/// Mutable solution access.
pub struct FlameStateMut<'a> {
    pub x: &'a mut Vec<f64>,
    pub mech: &'a Mechanism,
    pub grid: &'a Grid,
    pub natj: usize,
}

impl<'a> FlameStateMut<'a> {
    pub fn new(x: &'a mut Vec<f64>, mech: &'a Mechanism, grid: &'a Grid) -> Self {
        let n = natj(mech);
        FlameStateMut { x, mech, grid, natj: n }
    }

    pub fn set_temperature(&mut self, j: usize, t: f64) {
        let idx = idx_t(self.natj, j);
        self.x[idx] = t;
    }

    pub fn set_species(&mut self, k: usize, j: usize, yk: f64) {
        let idx = idx_y(self.natj, j, k);
        self.x[idx] = yk;
    }

    pub fn set_mass_flux(&mut self, m: f64) {
        let idx = idx_m(self.natj, self.grid.n_points());
        self.x[idx] = m;
    }
}

/// Build an initial guess solution vector.
/// - Linear temperature profile from T_unburned to T_burned
/// - Linear species profile from reactants to products
/// - Initial mass flux estimate based on flame speed guess
pub fn initial_guess(
    mech: &Mechanism,
    grid: &Grid,
    t_unburned: f64,
    t_burned: f64,
    y_reactants: &[f64],
    y_products: &[f64],
    mass_flux_guess: f64,
    z_center: f64,   // position of flame center (steep profile midpoint)
    z_width: f64,    // width of transition region
) -> Vec<f64> {
    let nk = mech.n_species();
    let nj = grid.n_points();
    let nv = natj(mech);
    let mut x = vec![0.0_f64; solution_length(mech, nj)];

    for j in 0..nj {
        let zj = grid.z[j];
        // Sigmoid-like weighting to create a smooth initial profile
        let s = sigmoid((zj - z_center) / z_width.max(1e-10));
        let t = t_unburned + s * (t_burned - t_unburned);
        x[idx_t(nv, j)] = t;
        for k in 0..nk {
            x[idx_y(nv, j, k)] = y_reactants[k] + s * (y_products[k] - y_reactants[k]);
        }
    }
    x[idx_m(nv, nj)] = mass_flux_guess;
    x
}

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Load an initial guess from a CSV file (e.g. produced by a Cantera solve).
///
/// Expected columns: `z [m]`, `T [K]`, `M [kg/m2/s]`, then `Y_<name>` for each
/// species in the mechanism (additional columns are ignored).
/// The profile is linearly interpolated onto `grid`.
pub fn initial_guess_from_csv(
    path: &str,
    mech: &Mechanism,
    grid: &Grid,
) -> Result<Vec<f64>> {
    let nk = mech.n_species();
    let nj = grid.n_points();
    let nv = natj(mech);

    // --- Parse CSV ---
    let mut rdr = csv::Reader::from_path(path)?;
    let headers: Vec<String> = rdr.headers()?.iter().map(|s| s.to_string()).collect();

    let col_z = headers.iter().position(|h| h == "z [m]")
        .ok_or_else(|| anyhow::anyhow!("CSV missing column 'z [m]'"))?;
    let col_t = headers.iter().position(|h| h == "T [K]")
        .ok_or_else(|| anyhow::anyhow!("CSV missing column 'T [K]'"))?;
    let col_m = headers.iter().position(|h| h == "M [kg/m2/s]")
        .ok_or_else(|| anyhow::anyhow!("CSV missing column 'M [kg/m2/s]'"))?;

    // Find Y column for each species (None → species absent in CSV, stays 0)
    let col_y: Vec<Option<usize>> = (0..nk)
        .map(|k| {
            let want = format!("Y_{}", mech.species[k].name);
            headers.iter().position(|h| h == &want)
        })
        .collect();

    let mut csv_z: Vec<f64> = Vec::new();
    let mut csv_t: Vec<f64> = Vec::new();
    let mut csv_m: Vec<f64> = Vec::new();
    let mut csv_y: Vec<Vec<f64>> = vec![Vec::new(); nk];

    for result in rdr.records() {
        let rec = result?;
        let parse = |c: usize| rec.get(c).unwrap_or("0").parse::<f64>().unwrap_or(0.0);
        csv_z.push(parse(col_z));
        csv_t.push(parse(col_t));
        csv_m.push(parse(col_m));
        for k in 0..nk {
            csv_y[k].push(col_y[k].map(|c| parse(c)).unwrap_or(0.0));
        }
    }
    anyhow::ensure!(!csv_z.is_empty(), "CSV file has no data rows");

    let m_val = csv_m.iter().copied().sum::<f64>() / csv_m.len() as f64;

    // --- Interpolate onto target grid ---
    let mut x = vec![0.0_f64; solution_length(mech, nj)];
    for j in 0..nj {
        let zj = grid.z[j];
        let (t_j, y_j) = interp_profile(zj, &csv_z, &csv_t, &csv_y, nk);
        x[idx_t(nv, j)] = t_j.max(1.0);
        // Renormalise species
        let y_sum: f64 = y_j.iter().sum();
        for k in 0..nk {
            x[idx_y(nv, j, k)] = if y_sum > 0.0 { (y_j[k] / y_sum).max(0.0) } else { 0.0 };
        }
    }
    x[idx_m(nv, nj)] = m_val.max(1e-6);
    Ok(x)
}

/// Linear interpolation at position `z` from sorted profile vectors.
fn interp_profile(
    z: f64,
    csv_z: &[f64],
    csv_t: &[f64],
    csv_y: &[Vec<f64>],
    nk: usize,
) -> (f64, Vec<f64>) {
    let n = csv_z.len();
    // Clamp to endpoints
    if z <= csv_z[0] {
        let y: Vec<f64> = (0..nk).map(|k| csv_y[k][0]).collect();
        return (csv_t[0], y);
    }
    if z >= csv_z[n - 1] {
        let y: Vec<f64> = (0..nk).map(|k| csv_y[k][n - 1]).collect();
        return (csv_t[n - 1], y);
    }
    // Binary search for bracket
    let i = csv_z.partition_point(|&zv| zv < z) - 1;
    let dz = csv_z[i + 1] - csv_z[i];
    let alpha = if dz > 0.0 { (z - csv_z[i]) / dz } else { 0.0 };
    let t = csv_t[i] + alpha * (csv_t[i + 1] - csv_t[i]);
    let y: Vec<f64> = (0..nk)
        .map(|k| csv_y[k][i] + alpha * (csv_y[k][i + 1] - csv_y[k][i]))
        .collect();
    (t, y)
}
