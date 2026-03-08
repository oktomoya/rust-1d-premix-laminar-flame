/// 1D spatial domain (non-uniform grid).

#[derive(Debug, Clone)]
pub struct Grid {
    /// Grid point positions [m]
    pub z: Vec<f64>,
}

impl Grid {
    /// Create a uniform grid from z=0 to z=length with n_points points.
    pub fn uniform(length: f64, n_points: usize) -> Self {
        let z = (0..n_points)
            .map(|i| length * i as f64 / (n_points - 1) as f64)
            .collect();
        Grid { z }
    }

    pub fn n_points(&self) -> usize {
        self.z.len()
    }

    /// Cell widths: dz[j] = z[j+1] - z[j], length = n_points - 1
    pub fn dz(&self) -> Vec<f64> {
        self.z.windows(2).map(|w| w[1] - w[0]).collect()
    }

    /// Midpoint positions: z_mid[j] = (z[j] + z[j+1]) / 2, length = n_points - 1
    pub fn z_mid(&self) -> Vec<f64> {
        self.z.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect()
    }

    /// Insert a new point between j and j+1 (midpoint).
    pub fn insert_midpoint(&mut self, j: usize) {
        let z_new = 0.5 * (self.z[j] + self.z[j + 1]);
        self.z.insert(j + 1, z_new);
    }
}
