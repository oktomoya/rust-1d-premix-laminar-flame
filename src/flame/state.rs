/// Solution state and flat vector layout for the Newton solver.
///
/// Solution vector layout (column-major by grid point):
///   x = [T_0, Y_{0,0}, ..., Y_{K-1,0},   T_1, Y_{0,1}, ...,   ...,   T_{J-1}, ..., M]
///   Length = (K+1)*J + 1
///   - Variables per point: NATJ = K + 1  (temperature + K species)
///   - Index of T at point j:          j * NATJ
///   - Index of Yk at point j:         j * NATJ + 1 + k
///   - Index of mass flux M:           NATJ * J   (last element)

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
        self.x[idx_y(self.natj, k, j)]
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
