/// Top-level solver driver: loads mechanism, builds initial guess, iterates
/// Newton + grid refinement until converged.

use anyhow::Result;
use crate::chemistry::parser::cantera_yaml;
use crate::chemistry::thermo::{density, enthalpy_molar, mean_molecular_weight};
use crate::flame::domain::Grid;
use crate::flame::refine::{refine, RefineCriteria};
use crate::flame::residual::FlameConfig as ResidualConfig;
use crate::flame::state::{initial_guess, solution_length};
use crate::io::input::FlameConfig;
use crate::io::output::{print_summary, write_csv};
use crate::solver::newton::{solve as newton_solve, NewtonConfig};
use crate::solver::pseudo_transient::{step as pt_step, PseudoTransientConfig};

pub fn run_flame(config: &FlameConfig) -> Result<()> {
    // --- Load mechanism ---
    let mech = match config.mechanism.format.as_str() {
        "cantera_yaml" | "yaml" => cantera_yaml::parse_file(&config.mechanism.file)?,
        other => anyhow::bail!("Unsupported mechanism format: {other}"),
    };
    println!("Loaded mechanism: {} species, {} reactions",
             mech.n_species(), mech.n_reactions());

    // --- Compute unburned mixture composition ---
    let (y_unburned, y_burned) = compute_compositions(&mech, config)?;
    let t_unburned = config.flame.t_unburned;
    let p = config.flame.pressure;

    // Estimate adiabatic flame temperature
    let t_ad = estimate_t_adiabatic(&mech, &y_unburned, t_unburned, p);
    println!("Estimated adiabatic flame temperature: {t_ad:.1} K");

    // --- Build initial grid and solution ---
    let mut grid = Grid::uniform(config.flame.domain_length, config.grid.initial_points);
    let z_center = 0.5 * config.flame.domain_length;  // put flame in the middle
    let z_width  = 0.1 * config.flame.domain_length;

    // Initial mass flux guess: Su ≈ 0.4 m/s (methane-air), ρ_unburned
    let w_u = mean_molecular_weight(&mech.species, &y_unburned);
    let rho_u = density(p, t_unburned, w_u);
    let su_guess = 0.4_f64; // m/s
    let m_guess = rho_u * su_guess;

    let mut x = initial_guess(
        &mech, &grid, t_unburned, t_ad,
        &y_unburned, &y_burned,
        m_guess, z_center, z_width,
    );

    // --- Residual configuration ---
    let t_fix = t_unburned + 200.0;
    let res_config = ResidualConfig {
        pressure: p,
        t_unburned,
        y_unburned: y_unburned.clone(),
        z_fix: z_center,
        t_fix,
    };

    // --- Phase 1: Pseudo-transient continuation ---
    let pt_config = PseudoTransientConfig {
        n_steps: config.solver.time_steps,
        dt_initial: config.solver.dt_initial,
        ..Default::default()
    };
    println!("--- Phase 1: Pseudo-transient continuation ---");
    pt_step(&mut x, &mech, &grid, &res_config, &pt_config)?;

    // --- Phase 2: Newton solve + adaptive refinement loop ---
    let newton_config = NewtonConfig {
        atol: config.solver.atol,
        rtol: config.solver.rtol,
        max_iter: config.solver.max_newton_iter,
        ..Default::default()
    };
    let refine_criteria = RefineCriteria {
        grad: config.grid.grad,
        curv: config.grid.curv,
        max_points: config.grid.max_points,
    };
    let mut res_config = res_config;

    for refine_pass in 0..20 {
        println!("--- Newton solve, pass {refine_pass} ({} grid points) ---", grid.n_points());
        newton_solve(&mut x, &mech, &grid, &res_config, &newton_config)?;

        // Update z_fix: find grid point where T first crosses t_fix
        res_config.z_fix = find_z_fix(&x, &mech, &grid, res_config.t_fix);

        // Adaptive refinement
        match refine(&x, &mech, &grid, &refine_criteria) {
            Some((new_grid, new_x)) => {
                println!("  Refined: {} → {} grid points", grid.n_points(), new_grid.n_points());
                grid = new_grid;
                x = new_x;
            }
            None => {
                println!("  Grid converged at {} points", grid.n_points());
                break;
            }
        }
    }

    // --- Output ---
    print_summary(&x, &mech, &grid, p);
    write_csv(&config.output.file, &x, &mech, &grid, p)?;

    Ok(())
}

/// Build reactant and product mass fractions from equivalence ratio +
/// fuel/oxidizer specifications.
///
/// Algorithm:
/// 1. Normalise fuel mole fractions to sum = 1.
/// 2. Compute stoichiometric O2 needed per mole of fuel mixture from element
///    counts (C→CO2 uses 1 O2/C, H→H2O uses 0.25 O2/H, O in fuel reduces need).
/// 3. Scale oxidizer so that O2_actual = O2_stoich / φ.
/// 4. Mix and convert to mass fractions.
fn compute_compositions(
    mech: &crate::chemistry::mechanism::Mechanism,
    config: &FlameConfig,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let nk = mech.n_species();
    let phi = config.flame.equivalence_ratio;

    // --- Normalise fuel (sum = 1) ---
    let fuel_total: f64 = config.flame.fuel.values().sum();
    anyhow::ensure!(fuel_total > 0.0, "Fuel mole fractions sum to zero");
    let mut x_fuel = vec![0.0_f64; nk];
    for (name, &frac) in &config.flame.fuel {
        if let Some(k) = mech.species_index(name) {
            x_fuel[k] = frac / fuel_total;
        }
    }

    // --- Normalise oxidizer (sum = 1) ---
    let ox_total: f64 = config.flame.oxidizer.values().sum();
    anyhow::ensure!(ox_total > 0.0, "Oxidizer mole fractions sum to zero");
    let mut x_ox = vec![0.0_f64; nk];
    for (name, &frac) in &config.flame.oxidizer {
        if let Some(k) = mech.species_index(name) {
            x_ox[k] = frac / ox_total;
        }
    }

    // --- Stoichiometric O2 per mole of fuel mixture ---
    // For species with composition {C:a, H:b, O:c, S:d}:
    //   O2_needed = a + b/4 - c/2 + d   (complete combustion to CO2 + H2O + SO2)
    let stoich_o2_per_fuel: f64 = x_fuel.iter().enumerate()
        .map(|(k, &xk)| {
            if xk == 0.0 { return 0.0; }
            let sp = &mech.species[k];
            let c = sp.composition.get("C").copied().unwrap_or(0.0);
            let h = sp.composition.get("H").copied().unwrap_or(0.0);
            let o = sp.composition.get("O").copied().unwrap_or(0.0);
            let s = sp.composition.get("S").copied().unwrap_or(0.0);
            xk * (c + h / 4.0 - o / 2.0 + s)
        })
        .sum();

    // --- O2 mole fraction in the oxidizer stream ---
    let o2_idx = mech.species_index("O2");
    let x_o2_in_ox = o2_idx.map(|k| x_ox[k]).unwrap_or(0.0);
    anyhow::ensure!(x_o2_in_ox > 0.0, "Oxidizer contains no O2");

    // --- Mixing ratio: n_fuel = 1, n_ox = stoich_o2 / (x_o2_in_ox * phi) ---
    let n_ox = if stoich_o2_per_fuel > 0.0 {
        stoich_o2_per_fuel / (x_o2_in_ox * phi)
    } else {
        // No combustible content — just mix 1:1 as fallback
        1.0
    };
    let total = 1.0 + n_ox;

    // --- Mixed mole fractions ---
    let mut x_mol = vec![0.0_f64; nk];
    for k in 0..nk {
        x_mol[k] = (x_fuel[k] + n_ox * x_ox[k]) / total;
    }

    // --- Mole → mass fractions ---
    let w_mean: f64 = mech.species.iter().zip(x_mol.iter())
        .map(|(s, &x)| x * s.molecular_weight)
        .sum();
    let y_unburned: Vec<f64> = mech.species.iter().zip(x_mol.iter())
        .map(|(s, &x)| x * s.molecular_weight / w_mean)
        .collect();

    let y_burned = estimate_burned_composition(mech, &x_mol);

    Ok((y_unburned, y_burned))
}

/// Estimate burned gas composition from element balance (complete combustion).
///
/// Given the unburned mole fractions, count atoms and distribute to products:
///   all C → CO2,  all H → H2O,  all N → N2,  remaining O2 stays,
///   inert species (Ar, He) pass through unchanged.
/// For rich flames (φ > 1) there may be insufficient O2; this estimate
/// still distributes available O2 to CO2 then H2O and leaves the rest as
/// CO / H2 if oxygen-limited — sufficient for an initial guess.
fn estimate_burned_composition(
    mech: &crate::chemistry::mechanism::Mechanism,
    x_unburned: &[f64],
) -> Vec<f64> {
    let nk = mech.n_species();

    // Count atoms in unburned mixture [mol per mol of mixture]
    let (mut n_c, mut n_h, mut n_o, mut n_n, mut n_s, mut n_ar, mut n_he) =
        (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for k in 0..nk {
        let xk = x_unburned[k];
        if xk == 0.0 { continue; }
        let sp = &mech.species[k];
        n_c  += xk * sp.composition.get("C" ).copied().unwrap_or(0.0);
        n_h  += xk * sp.composition.get("H" ).copied().unwrap_or(0.0);
        n_o  += xk * sp.composition.get("O" ).copied().unwrap_or(0.0);  // already ×2 for O2
        n_n  += xk * sp.composition.get("N" ).copied().unwrap_or(0.0);
        n_s  += xk * sp.composition.get("S" ).copied().unwrap_or(0.0);
        n_ar += xk * sp.composition.get("AR").copied().unwrap_or(0.0);
        n_he += xk * sp.composition.get("HE").copied().unwrap_or(0.0);
    }
    // O from fuel goes into the O pool (available O atoms)
    let mut o_avail = n_o; // O atoms

    // Build product mole amounts (mol per mol of original mixture)
    let mut x_prod = vec![0.0_f64; nk];
    let idx = |name: &str| mech.species_index(name);

    // CO2: each C needs 2 O atoms
    let co2_moles = if let Some(k) = idx("CO2") {
        let moles = n_c.min(o_avail / 2.0);
        x_prod[k] = moles;
        o_avail -= 2.0 * moles;
        n_c -= moles;
        moles
    } else { 0.0 };
    let _ = co2_moles;

    // H2O: each 2 H needs 1 O atom
    if let Some(k) = idx("H2O") {
        let h2o_moles = (n_h / 2.0).min(o_avail);
        x_prod[k] = h2o_moles;
        o_avail -= h2o_moles;
        n_h -= 2.0 * h2o_moles;
    }

    // Remaining C → CO (if O available) or soot (just leave as CO here)
    if n_c > 0.0 {
        if let Some(k) = idx("CO") {
            let co_moles = n_c.min(o_avail);
            x_prod[k] = co_moles;
            o_avail -= co_moles;
            n_c -= co_moles;
        }
    }

    // Remaining H → H2
    if n_h > 0.0 {
        if let Some(k) = idx("H2") {
            x_prod[k] = n_h / 2.0;
        }
    }

    // Excess O2
    if o_avail > 0.0 {
        if let Some(k) = idx("O2") {
            x_prod[k] = o_avail / 2.0;
        }
    }

    // N → N2
    if let Some(k) = idx("N2") {
        x_prod[k] = n_n / 2.0;
    }

    // S → SO2
    if n_s > 0.0 {
        if let Some(k) = idx("SO2") {
            x_prod[k] = n_s;
        }
    }

    // Inerts pass through
    if let Some(k) = idx("AR")  { x_prod[k] = n_ar; }
    if let Some(k) = idx("HE")  { x_prod[k] = n_he; }

    // Normalize to sum = 1, then convert to mass fractions
    let x_sum: f64 = x_prod.iter().sum();
    if x_sum > 0.0 {
        x_prod.iter_mut().for_each(|v| *v /= x_sum);
    }
    let w_mean: f64 = mech.species.iter().zip(x_prod.iter())
        .map(|(s, &x)| x * s.molecular_weight)
        .sum();
    if w_mean > 0.0 {
        mech.species.iter().zip(x_prod.iter_mut())
            .for_each(|(s, x)| *x *= s.molecular_weight / w_mean);
    }
    x_prod
}

/// Find the position z where T first crosses t_fix (left → right).
/// Used to update z_fix dynamically after each Newton pass.
fn find_z_fix(
    x: &[f64],
    mech: &crate::chemistry::mechanism::Mechanism,
    grid: &Grid,
    t_fix: f64,
) -> f64 {
    use crate::flame::state::{idx_t, natj};
    let nv = natj(mech);
    let nj = grid.n_points();
    for j in 0..nj - 1 {
        let t_j   = x[idx_t(nv, j)];
        let t_jp1 = x[idx_t(nv, j + 1)];
        if t_j <= t_fix && t_jp1 > t_fix {
            // Linear interpolation
            let alpha = (t_fix - t_j) / (t_jp1 - t_j);
            return grid.z[j] + alpha * (grid.z[j + 1] - grid.z[j]);
        }
    }
    // Fallback: midpoint
    grid.z[nj / 2]
}

/// Estimate adiabatic flame temperature using enthalpy balance.
/// H_reactants(T_u) = H_products(T_ad) — solved iteratively.
fn estimate_t_adiabatic(
    mech: &crate::chemistry::mechanism::Mechanism,
    y_unburned: &[f64],
    t_u: f64,
    _p: f64,
) -> f64 {
    // Rough estimate: T_ad ≈ T_u + ΔH_combustion / cp_products
    // Use NASA cp integrated from T_u upward
    let h_u: f64 = mech.species.iter().zip(y_unburned.iter())
        .map(|(s, &yk)| yk * crate::chemistry::thermo::enthalpy_molar(s, t_u) / s.molecular_weight)
        .sum();

    // Iterate T_ad
    let mut t_ad = t_u + 1500.0; // initial guess
    for _ in 0..100 {
        let h_ad: f64 = mech.species.iter().zip(y_unburned.iter())
            .map(|(s, &yk)| yk * crate::chemistry::thermo::enthalpy_molar(s, t_ad) / s.molecular_weight)
            .sum();
        let cp: f64 = crate::chemistry::thermo::cp_mixture(&mech.species, y_unburned, t_ad);
        let dt = (h_u - h_ad) / cp.max(1.0);
        t_ad += dt;
        if dt.abs() < 1.0 { break; }
    }
    t_ad.clamp(t_u + 100.0, 6000.0)
}
