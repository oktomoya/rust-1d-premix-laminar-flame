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
    use approx::assert_relative_eq;

    #[test]
    fn omega22_known_value() {
        // At T* = 1.0, Ω*(2,2) ≈ 1.593 (from literature tables)
        assert_relative_eq!(omega22(1.0), 1.593, epsilon = 0.01);
    }
}
