/// Geometry type for transport property calculation.
#[derive(Debug, Clone, PartialEq)]
pub enum GeometryType {
    Atom,
    Linear,
    Nonlinear,
}

/// Transport parameters sourced from a CHEMKIN-style transport database.
#[derive(Debug, Clone)]
pub struct TransportParams {
    pub geometry: GeometryType,
    /// Lennard-Jones well depth ε/kb [K]
    pub well_depth: f64,
    /// Lennard-Jones collision diameter σ [Angstrom]
    pub diameter: f64,
    /// Dipole moment [Debye]
    pub dipole_moment: f64,
    /// Polarizability [Angstrom³]
    pub polarizability: f64,
    /// Rotational relaxation collision number at 298 K
    pub rot_relax: f64,
}

/// NASA 7-coefficient polynomial for one temperature range.
/// cp/R = a[0] + a[1]*T + a[2]*T² + a[3]*T³ + a[4]*T⁴
/// h/RT = a[0] + a[1]/2*T + a[2]/3*T² + a[3]/4*T³ + a[4]/5*T⁴ + a[5]/T
/// s/R  = a[0]*ln(T) + a[1]*T + a[2]/2*T² + a[3]/3*T³ + a[4]/4*T⁴ + a[6]
#[derive(Debug, Clone)]
pub struct NasaPoly {
    pub t_low: f64,
    pub t_high: f64,
    pub coeffs: [f64; 7],
}

/// A chemical species with thermodynamic and transport data.
#[derive(Debug, Clone)]
pub struct Species {
    pub name: String,
    /// Molecular weight [kg/mol]
    pub molecular_weight: f64,
    /// High-temperature polynomial (T_mid … T_high)
    pub nasa_high: NasaPoly,
    /// Low-temperature polynomial (T_low … T_mid)
    pub nasa_low: NasaPoly,
    /// Transport parameters
    pub transport: TransportParams,
}

impl Species {
    /// Select the NASA polynomial valid at temperature `t`.
    /// At exactly T_mid, the low polynomial is used (matching Cantera: `T <= T_mid → low`).
    pub fn nasa_poly(&self, t: f64) -> &NasaPoly {
        if t > self.nasa_low.t_high {
            &self.nasa_high
        } else {
            &self.nasa_low
        }
    }
}
