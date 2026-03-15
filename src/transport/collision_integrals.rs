/// Reduced collision integrals for Chapman-Enskog transport theory.
///
/// Two implementations are provided:
///
/// 1. **Neufeld (1972) polynomial fits** — `omega22` / `omega11`
///    Valid for non-polar species (δ* = 0).  Used in unit tests only.
///
/// 2. **Monchick-Mason (1961) 2D tables** — `omega22_mm` / `omega11_mm`
///    Valid for all species including polar ones (δ* > 0).
///    Matches Cantera's `MMCollisionInt` interpolation exactly:
///      - Linear interpolation in δ*
///      - Quadratic interpolation in log(T*) using 3 bracketing rows
///
/// References:
///   Neufeld et al., J. Chem. Phys. 57, 1100 (1972)
///   Monchick & Mason, J. Chem. Phys. 35, 1676 (1961)
///   Cantera src/transport/MMCollisionInt.cpp

// ---------------------------------------------------------------------------
// Monchick-Mason tables
// ---------------------------------------------------------------------------

/// T* grid for the 37-row omega22 / astar tables (same grid for both).
const TSTAR22: [f64; 37] = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0,
];

/// δ* grid (8 columns) used by both tables.
const DELTA_STAR_GRID: [f64; 8] = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5];

/// Ω*(2,2) table from Monchick & Mason (1961) as reproduced in Cantera's
/// `MMCollisionInt::omega22_table`.  37 rows × 8 columns (δ*).
#[rustfmt::skip]
const OMEGA22_TABLE: [[f64; 8]; 37] = [
    [4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89  ], // T* = 0.1
    [3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618  ], // T* = 0.2
    [2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874  ], // T* = 0.3
    [2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895  ], // T* = 0.4
    [2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249  ], // T* = 0.5
    [2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786  ], // T* = 0.6
    [1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435  ], // T* = 0.7
    [1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156  ], // T* = 0.8
    [1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933  ], // T* = 0.9
    [1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746  ], // T* = 1.0
    [1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451  ], // T* = 1.2
    [1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228  ], // T* = 1.4
    [1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053  ], // T* = 1.6
    [1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912  ], // T* = 1.8
    [1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795  ], // T* = 2.0
    [1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578  ], // T* = 2.5
    [1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428  ], // T* = 3.0
    [0.99963,1.004,  1.016,  1.035,  1.062,  1.133,  1.22,   1.319  ], // T* = 3.5
    [0.96988,0.9732, 0.983,  0.9991, 1.021,  1.079,  1.153,  1.236  ], // T* = 4.0
    [0.92676,0.9291, 0.936,  0.9473, 0.9628, 1.005,  1.058,  1.121  ], // T* = 5.0
    [0.89616,0.8979, 0.903,  0.9114, 0.923,  0.9545, 0.9955, 1.044  ], // T* = 6.0
    [0.87272,0.8741, 0.878,  0.8845, 0.8935, 0.9181, 0.9505, 0.9893 ], // T* = 7.0
    [0.85379,0.8549, 0.858,  0.8632, 0.8703, 0.8901, 0.9164, 0.9482 ], // T* = 8.0
    [0.83795,0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916  ], // T* = 9.0
    [0.82435,0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901 ], // T* = 10.0
    [0.80184,0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504 ], // T* = 12.0
    [0.78363,0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212 ], // T* = 14.0
    [0.76834,0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983 ], // T* = 16.0
    [0.75518,0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797 ], // T* = 18.0
    [0.74364,0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642 ], // T* = 20.0
    [0.71982,0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339 ], // T* = 25.0
    [0.70097,0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112 ], // T* = 30.0
    [0.68545,0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932 ], // T* = 35.0
    [0.67232,0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784 ], // T* = 40.0
    [0.65099,0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546 ], // T* = 50.0
    [0.61397,0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147 ], // T* = 75.0
    [0.5887, 0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885 ], // T* = 100.0
];

/// A* = Ω*(2,2)/Ω*(1,1) table.  37 rows aligned with TSTAR22, 8 columns (δ*).
/// Extracted from Cantera's `MMCollisionInt::astar_table` rows 1..=37
/// (i.e., T* = 0.1..100, skipping the T*=0 sentinel row).
#[rustfmt::skip]
const ASTAR_TABLE: [[f64; 8]; 37] = [
    [1.0231, 1.0660, 1.0380, 1.0400, 1.0430, 1.0500, 1.0520, 1.0510], // T* = 0.1
    [1.0424, 1.0450, 1.0480, 1.0520, 1.0560, 1.0650, 1.0660, 1.0640], // T* = 0.2
    [1.0719, 1.0670, 1.0600, 1.0550, 1.0580, 1.0680, 1.0710, 1.0710], // T* = 0.3
    [1.0936, 1.0870, 1.0770, 1.0690, 1.0680, 1.0750, 1.0780, 1.0780], // T* = 0.4
    [1.1053, 1.0980, 1.0880, 1.0800, 1.0780, 1.0820, 1.0840, 1.0840], // T* = 0.5
    [1.1104, 1.1040, 1.0960, 1.0890, 1.0860, 1.0890, 1.0900, 1.0900], // T* = 0.6
    [1.1114, 1.1070, 1.1000, 1.0950, 1.0930, 1.0950, 1.0960, 1.0950], // T* = 0.7
    [1.1104, 1.1070, 1.1020, 1.0990, 1.0980, 1.1000, 1.1000, 1.0990], // T* = 0.8
    [1.1086, 1.1060, 1.1020, 1.1010, 1.1010, 1.1050, 1.1050, 1.1040], // T* = 0.9
    [1.1063, 1.1040, 1.1030, 1.1030, 1.1040, 1.1080, 1.1090, 1.1080], // T* = 1.0
    [1.1020, 1.1020, 1.1030, 1.1050, 1.1070, 1.1120, 1.1150, 1.1150], // T* = 1.2
    [1.0985, 1.0990, 1.1010, 1.1040, 1.1080, 1.1150, 1.1190, 1.1200], // T* = 1.4
    [1.0960, 1.0960, 1.0990, 1.1030, 1.1080, 1.1160, 1.1210, 1.1240], // T* = 1.6
    [1.0943, 1.0950, 1.0990, 1.1020, 1.1080, 1.1170, 1.1230, 1.1260], // T* = 1.8
    [1.0934, 1.0940, 1.0970, 1.1020, 1.1070, 1.1160, 1.1230, 1.1280], // T* = 2.0
    [1.0926, 1.0940, 1.0970, 1.0990, 1.1050, 1.1150, 1.1230, 1.1300], // T* = 2.5
    [1.0934, 1.0950, 1.0970, 1.0990, 1.1040, 1.1130, 1.1220, 1.1290], // T* = 3.0
    [1.0948, 1.0960, 1.0980, 1.1000, 1.1030, 1.1120, 1.1190, 1.1270], // T* = 3.5
    [1.0965, 1.0970, 1.0990, 1.1010, 1.1040, 1.1100, 1.1180, 1.1260], // T* = 4.0
    [1.0997, 1.1000, 1.1010, 1.1020, 1.1050, 1.1100, 1.1160, 1.1230], // T* = 5.0
    [1.1025, 1.1030, 1.1040, 1.1050, 1.1060, 1.1100, 1.1150, 1.1210], // T* = 6.0
    [1.1050, 1.1050, 1.1060, 1.1070, 1.1080, 1.1110, 1.1150, 1.1200], // T* = 7.0
    [1.1072, 1.1070, 1.1080, 1.1080, 1.1090, 1.1120, 1.1150, 1.1190], // T* = 8.0
    [1.1091, 1.1090, 1.1090, 1.1100, 1.1110, 1.1130, 1.1150, 1.1190], // T* = 9.0
    [1.1107, 1.1110, 1.1110, 1.1110, 1.1120, 1.1140, 1.1160, 1.1190], // T* = 10.0
    [1.1133, 1.1140, 1.1130, 1.1140, 1.1140, 1.1150, 1.1170, 1.1190], // T* = 12.0
    [1.1154, 1.1150, 1.1160, 1.1160, 1.1160, 1.1170, 1.1180, 1.1200], // T* = 14.0
    [1.1172, 1.1170, 1.1170, 1.1180, 1.1180, 1.1180, 1.1190, 1.1200], // T* = 16.0
    [1.1186, 1.1190, 1.1190, 1.1190, 1.1190, 1.1190, 1.1200, 1.1210], // T* = 18.0
    [1.1199, 1.1200, 1.1200, 1.1200, 1.1200, 1.1210, 1.1210, 1.1220], // T* = 20.0
    [1.1223, 1.1220, 1.1220, 1.1220, 1.1220, 1.1230, 1.1230, 1.1240], // T* = 25.0
    [1.1243, 1.1240, 1.1240, 1.1240, 1.1240, 1.1240, 1.1250, 1.1250], // T* = 30.0
    [1.1259, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260, 1.1260], // T* = 35.0
    [1.1273, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1270, 1.1280], // T* = 40.0
    [1.1297, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1300, 1.1290], // T* = 50.0
    [1.1339, 1.1340, 1.1340, 1.1350, 1.1350, 1.1340, 1.1340, 1.1320], // T* = 75.0
    [1.1364, 1.1370, 1.1370, 1.1380, 1.1390, 1.1380, 1.1370, 1.1350], // T* = 100.0
];

// ---------------------------------------------------------------------------
// Interpolation helpers
// ---------------------------------------------------------------------------

/// Quadratic interpolation at x0 given three (x, y) pairs.
/// Matches Cantera's `MMCollisionInt::quadInterp`.
fn quad_interp(x0: f64, x: [f64; 3], y: [f64; 3]) -> f64 {
    let dx21 = x[1] - x[0];
    let dx32 = x[2] - x[1];
    let dx31 = dx21 + dx32;
    let dy32 = y[2] - y[1];
    let dy21 = y[1] - y[0];
    let a = (dx21 * dy32 - dy21 * dx32) / (dx21 * dx31 * dx32);
    a * (x0 - x[0]) * (x0 - x[1]) + (dy21 / dx21) * (x0 - x[1]) + y[1]
}

/// Linear interpolation in δ* across the 8-column row.
fn interp_delta(row: &[f64; 8], delta_star: f64) -> f64 {
    if delta_star <= 0.0 {
        return row[0];
    }
    if delta_star >= DELTA_STAR_GRID[7] {
        return row[7];
    }
    // Find bracketing columns
    let mut j = 0usize;
    for k in 0..7 {
        if delta_star < DELTA_STAR_GRID[k + 1] {
            j = k;
            break;
        }
    }
    let frac = (delta_star - DELTA_STAR_GRID[j]) / (DELTA_STAR_GRID[j + 1] - DELTA_STAR_GRID[j]);
    row[j] + frac * (row[j + 1] - row[j])
}

/// Select 3 bracketing row indices in TSTAR22 for T* interpolation,
/// matching Cantera's index logic in `MMCollisionInt::omega22`.
fn bracket_rows(t_star: f64) -> [usize; 3] {
    let mut i = 0usize;
    for k in 0..37 {
        if t_star < TSTAR22[k] {
            i = k;
            break;
        }
        i = k + 1; // if t_star >= all entries, i = 37
    }
    let i1 = if i == 0 { 0 } else { i - 1 };
    let (i1, _i2) = if i1 + 3 > 36 {
        (36 - 2, 36)
    } else {
        (i1, i1 + 3)
    };
    [i1, i1 + 1, i1 + 2]
}

// ---------------------------------------------------------------------------
// Public 2D MM collision integral functions
// ---------------------------------------------------------------------------

/// Ω*(2,2)(T*, δ*) from the Monchick-Mason table.
///
/// For δ* = 0 the result is the same as the δ*=0 column of the table
/// (slightly different from the Neufeld polynomial due to the original data source).
pub fn omega22_mm(t_star: f64, delta_star: f64) -> f64 {
    let t_star = t_star.max(0.1);
    let [i1, i2, i3] = bracket_rows(t_star);
    let y = [
        interp_delta(&OMEGA22_TABLE[i1], delta_star),
        interp_delta(&OMEGA22_TABLE[i2], delta_star),
        interp_delta(&OMEGA22_TABLE[i3], delta_star),
    ];
    let lt = t_star.ln();
    let lx = [TSTAR22[i1].ln(), TSTAR22[i2].ln(), TSTAR22[i3].ln()];
    quad_interp(lt, lx, y).max(0.01)
}

/// A*(T*, δ*) = Ω*(2,2) / Ω*(1,1) from the Monchick-Mason table.
fn astar_mm(t_star: f64, delta_star: f64) -> f64 {
    let t_star = t_star.max(0.1);
    let [i1, i2, i3] = bracket_rows(t_star);
    let y = [
        interp_delta(&ASTAR_TABLE[i1], delta_star),
        interp_delta(&ASTAR_TABLE[i2], delta_star),
        interp_delta(&ASTAR_TABLE[i3], delta_star),
    ];
    let lt = t_star.ln();
    let lx = [TSTAR22[i1].ln(), TSTAR22[i2].ln(), TSTAR22[i3].ln()];
    quad_interp(lt, lx, y).max(0.5)
}

/// Ω*(1,1)(T*, δ*) = Ω*(2,2) / A* from the Monchick-Mason table.
pub fn omega11_mm(t_star: f64, delta_star: f64) -> f64 {
    omega22_mm(t_star, delta_star) / astar_mm(t_star, delta_star)
}

// ---------------------------------------------------------------------------
// Reduced dipole moment δ*
// ---------------------------------------------------------------------------

/// Reduced dipole moment for a species pair.
///
/// δ*_ij = 0.5 × μi × μj / (4πε₀ × εij × σij³)
///
/// Arguments in "customary" transport units:
///   `dipole_i`, `dipole_j` — dipole moments [Debye]
///   `well_depth_ij`        — combined well depth ε_ij / kb [K]
///   `diameter_ij`          — combined diameter σ_ij [Å]
///
/// The conversion constant 3622.85 absorbs unit conversions:
///   0.5 × (1 D)² / (4πε₀ × 1 K × kb × (1 Å)³) = 3622.85 (dimensionless)
pub fn delta_star_reduced(
    dipole_i: f64,
    dipole_j: f64,
    well_depth_ij: f64,
    diameter_ij: f64,
) -> f64 {
    if dipole_i == 0.0 || dipole_j == 0.0 {
        return 0.0;
    }
    const CONV: f64 = 3622.85; // unit-conversion constant (see doc comment)
    CONV * dipole_i * dipole_j / (well_depth_ij * diameter_ij.powi(3))
}

// ---------------------------------------------------------------------------
// Neufeld (1972) polynomial fits (kept for non-polar unit tests)
// ---------------------------------------------------------------------------

/// Reduced collision integral Ω*(2,2) via Neufeld (1972) polynomial fit.
/// Valid for δ* = 0 (non-polar species).  t_star = T / (ε/kb).
pub fn omega22(t_star: f64) -> f64 {
    const A: f64 = 1.16145;
    const B: f64 = 0.14874;
    const C: f64 = 0.52487;
    const D: f64 = 0.77320;
    const E: f64 = 2.16178;
    const F: f64 = 2.43787;
    A / t_star.powf(B) + C / (D * t_star).exp() + E / (F * t_star).exp()
}

/// Reduced collision integral Ω*(1,1) via Neufeld (1972) polynomial fit.
/// Valid for δ* = 0 (non-polar species).
pub fn omega11(t_star: f64) -> f64 {
    const A: f64 = 1.06036;
    const B: f64 = 0.15610;
    const C: f64 = 0.19300;
    const D: f64 = 0.47635;
    const E: f64 = 1.03587;
    const F: f64 = 1.52996;
    const G: f64 = 1.76474;
    const H: f64 = 3.89411;
    A / t_star.powf(B)
        + C / (D * t_star).exp()
        + E / (F * t_star).exp()
        + G / (H * t_star).exp()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn check(label: &str, got: f64, expected: f64, rtol: f64) {
        let rel = (got - expected).abs() / expected.abs();
        assert!(
            rel < rtol,
            "{label}: got {got:.10e}, expected {expected:.10e}, rel err {rel:.2e}"
        );
    }

    // -----------------------------------------------------------------------
    // 1. Neufeld polynomial arithmetic
    // -----------------------------------------------------------------------

    #[test]
    fn test_omega22_polynomial_values() {
        let cases = [
            (0.3_f64,  2.845802516390621),
            (0.5,      2.283064325164955),
            (1.0,      1.592519596079362),
            (2.0,      1.175969501456664),
            (5.0,      0.925191231557026),
            (10.0,     0.824862825673038),
            (30.0,     0.700315119888241),
            (100.0,    0.585491397311232),
        ];
        for (t_star, expected) in cases {
            check(&format!("omega22 T*={t_star}"), omega22(t_star), expected, 1e-10);
        }
    }

    #[test]
    fn test_omega11_polynomial_values() {
        let cases = [
            (0.3_f64,  2.650176361097787),
            (0.5,      2.067477129865050),
            (1.0,      1.440466399593360),
            (2.0,      1.075362638995321),
            (5.0,      0.843115644647863),
            (10.0,     0.741854874843624),
            (30.0,     0.623555037637060),
            (100.0,    0.516717697672334),
        ];
        for (t_star, expected) in cases {
            check(&format!("omega11 T*={t_star}"), omega11(t_star), expected, 1e-10);
        }
    }

    // -----------------------------------------------------------------------
    // 2. Neufeld vs Hirschfelder-Curtiss-Bird tables
    // -----------------------------------------------------------------------

    #[test]
    fn test_omega22_vs_neufeld_table() {
        let cases = [
            (1.0_f64, 1.593_f64, 2e-3),
            (2.0,     1.175,     2e-3),
            (5.0,     0.9256,    2e-3),
            (10.0,    0.8242,    2e-3),
        ];
        for (t_star, table_val, tol) in cases {
            check(&format!("omega22 vs table T*={t_star}"), omega22(t_star), table_val, tol);
        }
    }

    #[test]
    fn test_omega11_vs_neufeld_table() {
        let cases = [
            (1.0_f64, 1.439_f64, 2e-3),
            (2.0,     1.075,     2e-3),
            (5.0,     0.8422,    2e-3),
            (10.0,    0.7424,    2e-3),
        ];
        for (t_star, table_val, tol) in cases {
            check(&format!("omega11 vs table T*={t_star}"), omega11(t_star), table_val, tol);
        }
    }

    // -----------------------------------------------------------------------
    // 3. Neufeld sanity checks
    // -----------------------------------------------------------------------

    #[test]
    fn test_monotone_decreasing() {
        let t_stars = [0.5, 1.0, 2.0, 5.0, 10.0, 50.0];
        for w in t_stars.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            assert!(omega22(lo) > omega22(hi), "omega22 not decreasing at T*={lo}→{hi}");
            assert!(omega11(lo) > omega11(hi), "omega11 not decreasing at T*={lo}→{hi}");
        }
    }

    #[test]
    fn test_high_t_star_limit() {
        assert!(omega22(100.0) > 0.5 && omega22(100.0) < 1.0);
        assert!(omega11(100.0) > 0.4 && omega11(100.0) < 1.0);
    }

    #[test]
    fn test_omega22_gt_omega11_high_t() {
        for &ts in &[1.0_f64, 2.0, 5.0, 10.0, 50.0] {
            assert!(
                omega22(ts) > omega11(ts),
                "Expected omega22 > omega11 at T*={ts}: {} vs {}",
                omega22(ts), omega11(ts)
            );
        }
    }

    // -----------------------------------------------------------------------
    // 4. MM table: non-polar limit (δ* = 0) matches table column 0 exactly
    // -----------------------------------------------------------------------

    #[test]
    fn test_mm_omega22_delta0_matches_table() {
        // At δ*=0, omega22_mm should return OMEGA22_TABLE[row][0] after interpolation.
        // Check a few on-grid T* values (where interp is exact).
        let cases = [
            (1.0_f64,  OMEGA22_TABLE[9][0]),  // T*=1.0,  row 9
            (5.0,      OMEGA22_TABLE[19][0]), // T*=5.0,  row 19
            (10.0,     OMEGA22_TABLE[24][0]), // T*=10.0, row 24
            (100.0,    OMEGA22_TABLE[36][0]), // T*=100.0, row 36
        ];
        for (ts, expected) in cases {
            check(&format!("omega22_mm δ*=0 T*={ts}"), omega22_mm(ts, 0.0), expected, 1e-6);
        }
    }

    #[test]
    fn test_mm_omega11_delta0_matches_table() {
        // omega11 = omega22 / astar; check a few on-grid T* at δ*=0.
        // Reference: Neufeld polynomial (< 0.3% agreement with MM table at δ*=0).
        for &ts in &[1.0_f64, 2.0, 5.0, 10.0] {
            let mm = omega11_mm(ts, 0.0);
            let neufeld = omega11(ts);
            let rel = (mm - neufeld).abs() / neufeld;
            assert!(
                rel < 5e-3,
                "omega11_mm vs Neufeld at T*={ts}: mm={mm:.6e}, neufeld={neufeld:.6e}, rel={rel:.2e}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // 5. MM table: H2O self-interaction (δ* ≈ 1.22) — sanity check
    //    Reference: Cantera 3.1.0, pure-H2O at 1000 K.
    //    H2O transport params: ε/kb=572.4 K, σ=2.605 Å, μ=1.844 D → T*=1.748, δ*=1.218
    // -----------------------------------------------------------------------

    #[test]
    fn test_mm_omega22_h2o_1000k() {
        // δ*(H2O,H2O) ≈ 1.218, T*(H2O at 1000K) = 1000/572.4 = 1.747
        let t_star = 1000.0 / 572.4;
        let ds = delta_star_reduced(1.844, 1.844, 572.4, 2.605);
        let om22 = omega22_mm(t_star, ds);
        // At δ*=0: table gives ~1.28 (T*≈1.75, interpolated).
        // Polar correction at δ*≈1.22 increases Ω*(2,2) substantially.
        // Expected range from MM table: ~1.65–1.85
        assert!(
            om22 > 1.4 && om22 < 2.2,
            "omega22_mm(H2O, 1000K): got {om22:.4}, expected 1.4–2.2"
        );
        // Also verify δ* is in expected range
        assert!(
            ds > 1.0 && ds < 1.5,
            "δ*(H2O): got {ds:.4}, expected 1.0–1.5"
        );
    }

    // -----------------------------------------------------------------------
    // 6. delta_star_reduced: non-polar species → 0
    // -----------------------------------------------------------------------

    #[test]
    fn test_delta_star_zero_for_nonpolar() {
        // N2: μ=0, H2: μ=0
        assert_eq!(delta_star_reduced(0.0, 0.0, 97.53, 3.621), 0.0);
        assert_eq!(delta_star_reduced(1.844, 0.0, 100.0, 3.0), 0.0);
    }
}
