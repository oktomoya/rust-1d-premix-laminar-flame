/// Neufeld polynomial fits for reduced collision integrals.
/// Source: Neufeld et al., J. Chem. Phys. 57, 1100 (1972).

/// Reduced collision integral Ω*(2,2) for viscosity/thermal conductivity.
/// t_star = T / (ε/kb)  (reduced temperature)
pub fn omega22(t_star: f64) -> f64 {
    const A: f64 = 1.16145;
    const B: f64 = 0.14874;
    const C: f64 = 0.52487;
    const D: f64 = 0.77320;
    const E: f64 = 2.16178;
    const F: f64 = 2.43787;
    A / t_star.powf(B) + C / (D * t_star).exp() + E / (F * t_star).exp()
}

/// Reduced collision integral Ω*(1,1) for diffusion coefficients.
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
    // 1. Polynomial arithmetic — compare against Python reference values
    //    computed from the same Neufeld (1972) formula.  rtol = 1e-10 ensures
    //    the Rust implementation evaluates the polynomial exactly.
    // -----------------------------------------------------------------------

    #[test]
    fn test_omega22_polynomial_values() {
        // Reference: Python math library, identical formula
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
    // 2. Neufeld (1972) Table I cross-validation
    //    Reference values from Hirschfelder, Curtiss & Bird (1954) Table A.3-II
    //    as reproduced in the Neufeld polynomial fit paper.
    //    Neufeld claims accuracy within 0.2% for 0.3 ≤ T* ≤ 100.
    //    We test over the range where the fit is well-validated (T* ≥ 1).
    // -----------------------------------------------------------------------

    #[test]
    fn test_omega22_vs_neufeld_table() {
        // (T*, Ω*(2,2) table value, tolerance)
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
    // 3. Sanity checks
    // -----------------------------------------------------------------------

    /// Both integrals must decrease monotonically with T*.
    #[test]
    fn test_monotone_decreasing() {
        let t_stars = [0.5, 1.0, 2.0, 5.0, 10.0, 50.0];
        for w in t_stars.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            assert!(omega22(lo) > omega22(hi), "omega22 not decreasing at T*={lo}→{hi}");
            assert!(omega11(lo) > omega11(hi), "omega11 not decreasing at T*={lo}→{hi}");
        }
    }

    /// At high T*, both integrals approach ~1 from above (hard-sphere limit).
    #[test]
    fn test_high_t_star_limit() {
        assert!(omega22(100.0) > 0.5 && omega22(100.0) < 1.0);
        assert!(omega11(100.0) > 0.4 && omega11(100.0) < 1.0);
    }

    /// omega22 > omega11 for T* ≥ 1 (standard kinetic-theory result).
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
}
