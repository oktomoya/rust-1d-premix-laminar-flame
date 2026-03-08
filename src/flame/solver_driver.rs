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

    for refine_pass in 0..20 {
        println!("--- Newton solve, pass {refine_pass} ({} grid points) ---", grid.n_points());
        newton_solve(&mut x, &mech, &grid, &res_config, &newton_config)?;

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

/// Build reactant and (equilibrium-estimated) product mass fractions
/// from equivalence ratio + fuel/oxidizer specifications.
fn compute_compositions(
    mech: &crate::chemistry::mechanism::Mechanism,
    config: &FlameConfig,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let nk = mech.n_species();
    let phi = config.flame.equivalence_ratio;
    let p = config.flame.pressure;

    // Build mole fractions from fuel + oxidizer
    let mut x_mol = vec![0.0_f64; nk];
    let mut total = 0.0_f64;

    // Stoichiometric O2 for fuel (simplified: assume CxHyOz fuel)
    // For a general implementation, this requires knowing fuel composition.
    // Here we use the user-supplied mole fractions directly.
    for (name, &frac) in &config.flame.fuel {
        if let Some(k) = mech.species_index(name) {
            x_mol[k] += frac;
            total += frac;
        }
    }
    for (name, &frac) in &config.flame.oxidizer {
        if let Some(k) = mech.species_index(name) {
            x_mol[k] += frac;
            total += frac;
        }
    }
    // Normalize
    if total > 0.0 {
        x_mol.iter_mut().for_each(|x| *x /= total);
    }

    // Convert to mass fractions
    let w_mean: f64 = mech.species.iter().zip(x_mol.iter())
        .map(|(s, &x)| x * s.molecular_weight)
        .sum();
    let y_unburned: Vec<f64> = mech.species.iter().zip(x_mol.iter())
        .map(|(s, &x)| x * s.molecular_weight / w_mean)
        .collect();

    // For products, use a simple complete-combustion estimate
    // (a proper equilibrium solver would be needed for accuracy)
    let y_burned = estimate_burned_composition(mech, &y_unburned);

    Ok((y_unburned, y_burned))
}

/// Very rough burned gas estimate: convert fuel → CO2 + H2O, rest unchanged.
/// A proper implementation uses chemical equilibrium (CEA or element balance).
fn estimate_burned_composition(
    mech: &crate::chemistry::mechanism::Mechanism,
    y_unburned: &[f64],
) -> Vec<f64> {
    // Just return the unburned composition as placeholder
    // TODO: replace with proper equilibrium or element balance
    y_unburned.to_vec()
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
