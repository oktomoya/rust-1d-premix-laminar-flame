/// Dense matrix storage and LU factorization for Newton solve.
///
/// Uses faer's partial-pivot LU to correctly handle dense coupling,
/// including the global mass-flux eigenvalue M which couples all
/// grid-point residuals and cannot be captured by a banded layout.

use anyhow::Result;
use faer::Mat;
use faer::prelude::Solve;

pub struct DenseMatrix {
    pub n: usize,
    data: Mat<f64>,
}

impl DenseMatrix {
    pub fn new(n: usize) -> Self {
        DenseMatrix {
            n,
            data: Mat::zeros(n, n),
        }
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        self.data[(row, col)] = val;
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[(row, col)]
    }

    /// Solve A·x = b in-place using faer's partial-pivot LU.
    /// b is overwritten with the solution.
    pub fn solve(&self, b: &mut Vec<f64>) -> Result<()> {
        let n = self.n;
        // Build a column vector from b.
        let rhs = Mat::from_fn(n, 1, |i, _| b[i]);
        // Factor and solve.
        let lu = self.data.partial_piv_lu();
        let sol = lu.solve(&rhs);
        for i in 0..n {
            b[i] = sol[(i, 0)];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dense_tridiagonal_solve() {
        // Same 5×5 tridiagonal system as banded tests.
        // Exact solution: u_j = j*(6-j), i.e. [5, 8, 9, 8, 5].
        let n = 5;
        let mut mat = DenseMatrix::new(n);
        for j in 0..n {
            mat.set(j, j, 2.0);
            if j > 0     { mat.set(j, j - 1, -1.0); }
            if j < n - 1 { mat.set(j, j + 1, -1.0); }
        }
        let mut b = vec![2.0_f64; n];
        mat.solve(&mut b).unwrap();
        let expected = [5.0_f64, 8.0, 9.0, 8.0, 5.0];
        for (j, (&u, &ex)) in b.iter().zip(expected.iter()).enumerate() {
            assert!((u - ex).abs() < 1e-10, "u[{j}] = {u:.6}, expected {ex:.6}");
        }
    }

    #[test]
    fn test_dense_with_global_coupling() {
        // A system where the last column couples all rows — like the M variable.
        // 4×4: A[i][i]=2 for all i, A[i][3]=1 for all i (global last-column coupling).
        // For i=3: the loop sets A[3][3]=2 then A[3][3]=1; restore to 2 explicitly.
        // Final matrix rows: [2,0,0,1], [0,2,0,1], [0,0,2,1], [0,0,0,2].
        // b = [3, 3, 3, 3].
        // Row 3: 2*x3 = 3 → x3 = 1.5.
        // Row i<3: 2*xi + 1.5 = 3 → xi = 0.75.
        let n = 4;
        let mut mat = DenseMatrix::new(n);
        for i in 0..n {
            mat.set(i, i, 2.0);
            mat.set(i, n - 1, 1.0); // global coupling column (overwrites diagonal for i=n-1)
        }
        mat.set(n - 1, n - 1, 2.0); // restore diagonal for last row
        let mut b = vec![3.0_f64; n];
        mat.solve(&mut b).unwrap();
        let expected = [0.75_f64, 0.75, 0.75, 1.5];
        for (j, (&u, &ex)) in b.iter().zip(expected.iter()).enumerate() {
            assert!((u - ex).abs() < 1e-10, "x[{j}] = {u:.6}, expected {ex:.6}");
        }
    }
}
