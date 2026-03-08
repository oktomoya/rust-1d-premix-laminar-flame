/// Banded matrix storage and LU factorization for Newton solve.
///
/// Layout: LAPACK-compatible general banded format (GB).
/// Storage array ab[kl+ku+1][n], where:
///   ab[ku + i - j][j] = A[i][j]   for  max(0, j-ku) <= i <= min(n-1, j+kl)

pub struct BandedMatrix {
    pub n: usize,
    pub kl: usize,    // lower bandwidth
    pub ku: usize,    // upper bandwidth
    /// Column-major storage: data[col * (2*kl+ku+1) + (ku + row - col)]
    pub data: Vec<f64>,
    pub ldab: usize,
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
        }
    }

    #[inline]
    fn storage_idx(&self, row: usize, col: usize) -> usize {
        // LAPACK band storage: ab[ku + kl + row - col][col]
        // Extra kl rows for LAPACK pivot space
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

    /// Solve A·x = b in-place (b is overwritten with solution).
    /// Expands the band to a dense matrix and applies Gaussian elimination
    /// with partial pivoting. Suitable for correctness; replace with a true
    /// banded LU (LAPACK dgbsv) for performance on large systems.
    pub fn solve(&mut self, b: &mut Vec<f64>) -> anyhow::Result<()> {
        let n = self.n;

        // Expand band storage → dense matrix (row-major)
        let mut a = vec![0.0_f64; n * n];
        for col in 0..n {
            let row_min = col.saturating_sub(self.ku);
            let row_max = (col + self.kl + 1).min(n);
            for row in row_min..row_max {
                a[row * n + col] = self.get(row, col);
            }
        }

        // Gaussian elimination with partial pivoting
        let mut piv = vec![0usize; n];
        for k in 0..n {
            // Find pivot
            let (max_row, _) = (k..n)
                .map(|i| (i, a[i * n + k].abs()))
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .unwrap_or((k, 0.0));
            piv[k] = max_row;

            if max_row != k {
                for j in 0..n { a.swap(k * n + j, max_row * n + j); }
                b.swap(k, max_row);
            }

            let akk = a[k * n + k];
            if akk.abs() < 1e-300 {
                anyhow::bail!("Singular matrix in banded solve at column {k}");
            }
            for i in k + 1..n {
                let factor = a[i * n + k] / akk;
                a[i * n + k] = 0.0;
                for j in k + 1..n {
                    let v = a[k * n + j];
                    a[i * n + j] -= factor * v;
                }
                b[i] -= factor * b[k];
            }
        }

        // Back substitution
        for k in (0..n).rev() {
            for j in k + 1..n {
                b[k] -= a[k * n + j] * b[j];
            }
            b[k] /= a[k * n + k];
        }

        Ok(())
    }
}
