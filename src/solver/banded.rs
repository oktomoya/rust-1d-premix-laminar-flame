/// Banded matrix storage and LU factorization for Newton solve.
///
/// Layout: LAPACK-compatible general banded format (GB).
/// Storage: data[col * ldab + (kl + ku + row - col)] = A[row][col]
///   for max(0, col-ku) <= row <= min(n-1, col+kl).
///
/// `ldab = 2*kl + ku + 1` allocates an extra `kl` rows per column for the
/// fill-in created by partial pivoting (same convention as LAPACK `dgbtrf`).

pub struct BandedMatrix {
    pub n: usize,
    pub kl: usize,    // lower bandwidth
    pub ku: usize,    // upper bandwidth
    pub data: Vec<f64>,
    pub ldab: usize,
    ipiv: Vec<usize>, // pivot row indices; filled by factor_in_place
}

impl BandedMatrix {
    pub fn new(n: usize, kl: usize, ku: usize) -> Self {
        let ldab = 2 * kl + ku + 1;
        BandedMatrix {
            n,
            kl,
            ku,
            data: vec![0.0_f64; ldab * n],
            ldab,
            ipiv: vec![0usize; n],
        }
    }

    #[inline]
    fn storage_idx(&self, row: usize, col: usize) -> usize {
        col * self.ldab + (self.kl + self.ku + row - col)
    }

    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        let idx = self.storage_idx(row, col);
        self.data[idx] = val;
    }

    pub fn get(&self, row: usize, col: usize) -> f64 {
        let idx = self.storage_idx(row, col);
        self.data[idx]
    }

    /// Factor the banded matrix in-place as P·A = L·U using Gaussian elimination
    /// with partial pivoting (max-element within the lower band of each column).
    /// Pivot row indices are stored in `self.ipiv` for use by `solve_factored`.
    ///
    /// Cost: O(n · kl · ku).
    pub fn factor_in_place(&mut self) -> anyhow::Result<()> {
        let n  = self.n;
        let kl = self.kl;
        let ku = self.ku;
        self.ipiv.resize(n, 0);

        for k in 0..n {
            // Find pivot: row with max |A[i][k]| in i = k..min(k+kl, n-1).
            let row_end = (k + kl + 1).min(n);
            let pivot_row = (k..row_end)
                .max_by(|&i, &j| {
                    self.get(i, k).abs().partial_cmp(&self.get(j, k).abs())
                        .unwrap_or(std::cmp::Ordering::Less)
                })
                .unwrap_or(k);
            self.ipiv[k] = pivot_row;

            // Swap rows k and pivot_row in band storage.
            // Range: columns k..min(pivot_row + ku + 1, n).
            // All storage_idx computations stay in [0, ldab-1] because:
            //   pivot_row - k <= kl, so indices shift by at most kl within ldab.
            if pivot_row != k {
                let swap_end = (pivot_row + ku + 1).min(n);
                for j in k..swap_end {
                    let idx_k = self.storage_idx(k, j);
                    let idx_p = self.storage_idx(pivot_row, j);
                    self.data.swap(idx_k, idx_p);
                }
            }

            let akk = self.get(k, k);
            if akk.abs() < 1e-300 {
                anyhow::bail!("Zero pivot in banded LU at column {k}");
            }

            // Eliminate rows k+1..row_end using the pivot row.
            // After pivoting, row k may have entries up to column k+ku+kl
            // (fill-in from the pivot swap), so col_end extends to ku+kl.
            let col_end = (k + ku + kl + 1).min(n);
            for i in k + 1..row_end {
                let factor = self.get(i, k) / akk;
                self.set(i, k, factor); // store L multiplier in place
                for j in k + 1..col_end {
                    let update = self.get(i, j) - factor * self.get(k, j);
                    self.set(i, j, update);
                }
            }
        }
        Ok(())
    }

    /// Solve P·L·U·x = b using the factored matrix from `factor_in_place`.
    /// `b` is overwritten with the solution.
    pub fn solve_factored(&self, b: &mut [f64]) {
        let n  = self.n;
        let kl = self.kl;
        let ku = self.ku;

        // Apply row permutation P to b.
        for k in 0..n {
            b.swap(k, self.ipiv[k]);
        }
        // Forward substitution: L·y = P·b  (L is unit lower triangular).
        for k in 0..n {
            let row_end = (k + kl + 1).min(n);
            for i in k + 1..row_end {
                b[i] -= self.get(i, k) * b[k];
            }
        }
        // Backward substitution: U·x = y.
        // U may have fill-in extending to column k+ku+kl after pivoting.
        for k in (0..n).rev() {
            let col_end = (k + ku + kl + 1).min(n);
            for j in k + 1..col_end {
                b[k] -= self.get(k, j) * b[j];
            }
            b[k] /= self.get(k, k);
        }
    }

    /// Factor and solve A·x = b in-place.
    /// Equivalent to calling `factor_in_place` then `solve_factored`.
    pub fn solve(&mut self, b: &mut Vec<f64>) -> anyhow::Result<()> {
        self.factor_in_place()?;
        self.solve_factored(b);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_banded_tridiagonal_solve() {
        // Solve -u_{j-1} + 2*u_j - u_{j+1} = 2 for j = 0..5,
        // with u_{-1} = u_5 = 0 (Dirichlet, absorbed into RHS).
        // Exact solution: u_j = j*(6-j), i.e. [5, 8, 9, 8, 5].
        let n = 5;
        let mut mat = BandedMatrix::new(n, 1, 1);
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
    fn test_banded_wider_band() {
        // 4×4 system with kl=2, ku=1; exact solution x = [1, 1, 1, 1].
        // A = [[4,1,0,0],[2,5,2,0],[1,3,6,1],[0,2,4,7]], b = [5, 9, 11, 13].
        let n = 4;
        let mut mat = BandedMatrix::new(n, 2, 1);
        mat.set(0, 0, 4.0); mat.set(0, 1, 1.0);
        mat.set(1, 0, 2.0); mat.set(1, 1, 5.0); mat.set(1, 2, 2.0);
        mat.set(2, 0, 1.0); mat.set(2, 1, 3.0); mat.set(2, 2, 6.0); mat.set(2, 3, 1.0);
                             mat.set(3, 1, 2.0); mat.set(3, 2, 4.0); mat.set(3, 3, 7.0);
        let mut b = vec![5.0_f64, 9.0, 11.0, 13.0];
        mat.solve(&mut b).unwrap();
        for (j, &u) in b.iter().enumerate() {
            assert!((u - 1.0).abs() < 1e-10, "x[{j}] = {u:.6}, expected 1.0");
        }
    }

    #[test]
    fn test_banded_pivoting_required() {
        // A[0][0] = 0 → without pivoting this would immediately fail.
        // kl = ku = 1, n = 3.
        // A = [[0, 2, 0], [1, 0, 2], [0, 1, 3]], b = [4, 5, 8].
        // Exact solution: x = [1, 2, 2].
        let n = 3;
        let mut mat = BandedMatrix::new(n, 1, 1);
        mat.set(0, 1, 2.0);
        mat.set(1, 0, 1.0);
        mat.set(1, 2, 2.0);
        mat.set(2, 1, 1.0);
        mat.set(2, 2, 3.0);
        let mut b = vec![4.0_f64, 5.0, 8.0];
        mat.solve(&mut b).unwrap();
        let expected = [1.0_f64, 2.0, 2.0];
        for (j, (&u, &ex)) in b.iter().zip(expected.iter()).enumerate() {
            assert!((u - ex).abs() < 1e-10, "x[{j}] = {u:.6}, expected {ex:.6}");
        }
    }
}
