/// Residual evaluation F(x) for the 1D freely propagating premixed flame.
///
/// Governing equations (SI units throughout):
///
/// Interior points j = 1 … J-2:
///   Species k (k = 0 … K-2):
///     F_yk = M * (Yk_j - Yk_{j-1})/dz_m  -  (jk_j - jk_{j-1}) / dz_av  -  ωk * Wk
///   Species K-1: closure  F_yK = 1 - ΣYk - Y_{K,j}   or same as others with correction vel.
///   Energy:
///     F_T = -M * cp * (T_j - T_{j-1})/dz_m
///           + (λ_j * (T_{j+1}-T_j)/dz_p - λ_{j-1} * (T_j-T_{j-1})/dz_m) / dz_av
///           - Σk jk_{j-1/2} * cpk * dT/dz (enthalpy transport term)
///           - Σk ωk * hk * Wk               (heat release)
///   Continuity (eigenvalue):
///     F_M = M_j - M_{j-1}   (trivially = 0 enforces dM/dz = 0)
///
/// Left boundary (j = 0):
///   F_yk = M * Yk_0 + jk_0 - M * Yku_k    (inlet flux condition)
///   F_T  = T_0 - T_unburned
///
/// Right boundary (j = J-1):
///   F_yk = Yk_{J-1} - Yk_{J-2}   (zero gradient)
///   F_T  = T_{J-1} - T_{J-2}
///
/// Eigenvalue closure: at the fixed point j = j_fix:
///   F_M = T_{j_fix} - T_fix   (replaces the trivial M equation at that point)

use crate::chemistry::kinetics::production_rates;
use crate::chemistry::mechanism::Mechanism;
use crate::chemistry::thermo::{
    cp_mixture, cp_species, density, enthalpy_molar, mean_molecular_weight,
};
use crate::flame::domain::Grid;
use crate::flame::state::{idx_m, idx_t, idx_y, natj, FlameState};
use crate::transport::mixture::{
    mass_to_mole_fractions, mixture_diffusion_coefficients, mixture_thermal_conductivity,
};

pub struct FlameConfig {
    pub pressure: f64,       // Pa
    pub t_unburned: f64,     // K
    pub y_unburned: Vec<f64>, // mass fractions at inlet
    pub z_fix: f64,          // position where T is pinned [m]
    pub t_fix: f64,          // pinned temperature [K]
}

/// Evaluate residual F(x) → rhs (length = solution_length).
pub fn eval_residual(
    x: &[f64],
    rhs: &mut [f64],
    mech: &Mechanism,
    grid: &Grid,
    config: &FlameConfig,
) {
    let nk = mech.n_species();
    let nj = grid.n_points();
    let nv = natj(mech);
    let p = config.pressure;

    // Zero out residual
    rhs.iter_mut().for_each(|r| *r = 0.0);

    let state = FlameState::new(x, mech, grid);
    let m = state.mass_flux();

    // Pre-compute midpoint transport properties: λ[j] and Dk[k][j] at interval midpoints
    // (between j and j+1), indexed 0..nj-1.
    let dz = grid.dz();

    let mut lambda_mid = vec![0.0_f64; nj - 1]; // thermal conductivity at midpoints
    let mut dk_mid = vec![vec![0.0_f64; nj - 1]; nk]; // diffusion coeff at midpoints
    let mut rho_mid = vec![0.0_f64; nj - 1];

    for j in 0..nj - 1 {
        // Average T and Y between j and j+1
        let t_av = 0.5 * (state.temperature(j) + state.temperature(j + 1));
        let y_av: Vec<f64> = (0..nk)
            .map(|k| 0.5 * (state.species(k, j) + state.species(k, j + 1)))
            .collect();
        let x_av = mass_to_mole_fractions(mech, &y_av);
        let w_mean = mean_molecular_weight(&mech.species, &y_av);
        rho_mid[j] = density(p, t_av, w_mean);
        lambda_mid[j] = mixture_thermal_conductivity(mech, &x_av, &y_av, t_av);
        let dk = mixture_diffusion_coefficients(mech, &x_av, t_av, p);
        for k in 0..nk {
            dk_mid[k][j] = dk[k];
        }
    }

    // Compute diffusion fluxes jk = ρ*Yk*Vk [kg/(m²·s)] at each midpoint j (between j and j+1).
    // Mixture-averaged with correction velocity:
    //   jk_raw = -ρ * Dk * dYk/dz
    //   Vc = Σk jk_raw  (correction to ensure Σjk = 0)
    //   jk = jk_raw - Yk * ρ * Vc / ... (correction velocity term)
    let mut jk_mid = vec![vec![0.0_f64; nj - 1]; nk]; // jk[k][j]

    for j in 0..nj - 1 {
        let mut jk_raw = vec![0.0_f64; nk];
        let mut sum_jk = 0.0_f64;
        for k in 0..nk {
            let dy_dz = (state.species(k, j + 1) - state.species(k, j)) / dz[j];
            jk_raw[k] = -rho_mid[j] * dk_mid[k][j] * dy_dz;
            sum_jk += jk_raw[k];
        }
        // Correction velocity: Vc such that Σ(jk + Yk * ρ * Vc) = 0
        // Applied as: jk_corrected = jk_raw - Yk_av * sum_jk
        for k in 0..nk {
            let yk_av = 0.5 * (state.species(k, j) + state.species(k, j + 1));
            jk_mid[k][j] = jk_raw[k] - yk_av * sum_jk;
        }
    }

    // --- Left boundary (j = 0) ---
    {
        let j = 0;
        let t_j = state.temperature(j);
        rhs[idx_t(nv, j)] = t_j - config.t_unburned;

        let mut sum_y = 0.0_f64;
        for k in 0..nk - 1 {
            let yk = state.species(k, j);
            // Inlet flux BC: M * Yk - jk_{0} = M * Yku_k → F = M*(Yk - Yku) + jk
            // For a no-diffusion-at-inlet approximation: F = Yk - Yku
            rhs[idx_y(nv, j, k)] = yk - config.y_unburned[k];
            sum_y += yk;
        }
        // Last species: sum constraint
        rhs[idx_y(nv, j, nk - 1)] = state.species(nk - 1, j) - config.y_unburned[nk - 1];
    }

    // --- Interior points j = 1 … J-2 ---
    let j_fix = find_fixed_point(grid, config.z_fix);

    for j in 1..nj - 1 {
        let t_j = state.temperature(j);
        let t_jm1 = state.temperature(j - 1);
        let t_jp1 = state.temperature(j + 1);
        let dz_m = dz[j - 1];
        let dz_p = dz[j];
        let dz_av = 0.5 * (dz_m + dz_p);

        let y_j: Vec<f64> = (0..nk).map(|k| state.species(k, j)).collect();
        let w_mean = mean_molecular_weight(&mech.species, &y_j);
        let rho_j = density(p, t_j, w_mean);
        let x_j = mass_to_mole_fractions(mech, &y_j);

        // Chemical production rates ωk [mol/(m³·s)]
        let conc: Vec<f64> = (0..nk)
            .map(|k| rho_j * y_j[k] / mech.species[k].molecular_weight)
            .collect();
        let wdot = production_rates(mech, t_j, &conc);

        // Species equations
        let mut sum_y = 0.0_f64;
        for k in 0..nk - 1 {
            let yk_j  = state.species(k, j);
            let yk_jm1 = state.species(k, j - 1);
            let convec = m * (yk_j - yk_jm1) / dz_m;
            let diffus = (jk_mid[k][j] - jk_mid[k][j - 1]) / dz_av;
            let source = wdot[k] * mech.species[k].molecular_weight;
            rhs[idx_y(nv, j, k)] = convec + diffus - source;
            sum_y += yk_j;
        }
        // Last species: correction velocity closure
        rhs[idx_y(nv, j, nk - 1)] = sum_y + state.species(nk - 1, j) - 1.0;

        // Energy equation
        let cp_j = cp_mixture(&mech.species, &y_j, t_j);
        let lambda_j   = lambda_mid[j];      // midpoint j (between j and j+1)
        let lambda_jm1 = lambda_mid[j - 1];  // midpoint j-1 (between j-1 and j)

        let conduction = (lambda_j * (t_jp1 - t_j) / dz_p
            - lambda_jm1 * (t_j - t_jm1) / dz_m) / dz_av;

        // Enthalpy transport: Σk jk_{av} * cpk * dT/dz
        let dt_dz = (t_jp1 - t_jm1) / (dz_m + dz_p);
        let mut enthalpy_transport = 0.0_f64;
        for k in 0..nk {
            let jk_av = 0.5 * (jk_mid[k][j].min_or(j, &jk_mid) + jk_mid[k][j - 1]);
            let cp_k = cp_species(&mech.species[k], t_j);
            enthalpy_transport += jk_av * cp_k * dt_dz;
        }

        // Heat release: Σk ωk * hk * Wk  (note hk in J/mol, Wk in kg/mol)
        let heat_release: f64 = (0..nk)
            .map(|k| wdot[k] * enthalpy_molar(&mech.species[k], t_j))
            .sum();

        // Residual: ρ*cp*∂T/∂t = ...  but steady state → F_T = 0
        // Written as: 0 = -M*cp*dT/dz_upwind + conduction - enthalpy_transport - heat_release
        let convec_t = m * cp_j * (t_j - t_jm1) / dz_m;
        rhs[idx_t(nv, j)] = -convec_t + conduction - enthalpy_transport - heat_release;

        // Continuity / eigenvalue
        if j == j_fix {
            // Fix temperature to locate flame
            rhs[idx_m(nv, nj)] = t_j - config.t_fix;
        }
    }

    // --- Right boundary (j = J-1) ---
    {
        let j = nj - 1;
        rhs[idx_t(nv, j)] = state.temperature(j) - state.temperature(j - 1);
        for k in 0..nk {
            rhs[idx_y(nv, j, k)] = state.species(k, j) - state.species(k, j - 1);
        }
    }

    // Mass flux equations: dM/dz = 0 → M_j = M_{j-1} (except at j_fix)
    // The global M is the last unknown; individual continuity eqs are trivially satisfied.
    // (This simplified layout uses a single M; see DESIGN.md section 3.1 for discussion.)
}

/// Find the grid index closest to z_fix.
fn find_fixed_point(grid: &Grid, z_fix: f64) -> usize {
    grid.z.iter().enumerate()
        .min_by(|(_, a), (_, b)| {
            ((*a - z_fix).abs()).partial_cmp(&((*b - z_fix).abs())).unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(grid.n_points() / 2)
}

// Helper trait to avoid borrow issues in inner loop
trait MinOr {
    fn min_or(&self, j: usize, full: &Vec<Vec<f64>>) -> f64;
}
impl MinOr for f64 {
    fn min_or(&self, _j: usize, _full: &Vec<Vec<f64>>) -> f64 { *self }
}
