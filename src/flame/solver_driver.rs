/// Top-level solver driver: loads mechanism, builds initial guess, iterates
/// Newton + grid refinement until converged.

use anyhow::Result;
use crate::chemistry::parser::cantera_yaml;
use crate::chemistry::thermo::{density, mean_molecular_weight};
use crate::flame::domain::Grid;
use crate::flame::refine::{refine, RefineCriteria};
use crate::flame::residual::FlameConfig as ResidualConfig;
use crate::flame::residual::eval_residual;
use crate::flame::state::{initial_guess, initial_guess_from_csv, solution_length};
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
    // compute_compositions returns complete-combustion products (H₂O + N₂, no radicals).
    let (y_unburned, y_burned_complete) = compute_compositions(&mech, config)?;
    let t_unburned = config.flame.t_unburned;
    let p = config.flame.pressure;

    // Estimate T_ad from complete-combustion products.  The radical-seeded
    // composition (see below) has higher stored enthalpy and would give a
    // spuriously low T_ad if used here.
    let t_ad = estimate_t_adiabatic(&mech, &y_unburned, &y_burned_complete, t_unburned, p);
    println!("Estimated adiabatic flame temperature: {t_ad:.1} K");

    // Seed equilibrium radicals into the burned composition for the initial
    // profile.  Without H/OH, chain-branching in the PT phase cannot start
    // (initiation is ~10¹² × slower than branching).
    let y_burned = seed_radicals(&mech, y_burned_complete);

    // --- Build initial grid and solution ---
    let mut grid = Grid::uniform(config.flame.domain_length, config.grid.initial_points);

    // Place the initial profile center at 40% of the domain.
    let z_center = 0.40 * config.flame.domain_length;

    // Initial flame width: smooth sigmoid spanning ~5 grid cells.
    // δ_T(H2/air) ≈ 0.5 mm; 5 × dz spans approximately one flame thickness.
    let dz = config.flame.domain_length / config.grid.initial_points as f64;
    let z_width = 5.0 * dz;

    // Eigenvalue reference point: temperature equals t_fix at z_fix.
    //
    // Strategy: place z_fix at the initial flame centre (z_center) and set
    // t_fix to the midpoint temperature between T_u and T_ad.  Advantages:
    //   1. The constraint T(z_fix) = t_fix is EXACTLY satisfied by the initial
    //      sigmoid guess (sigmoid at z_center gives T = midpoint).
    //   2. z_fix is in the reaction zone (where T changes fastest) →
    //      eigenvalue equation is well-conditioned throughout the solve.
    //   3. Placing z_fix at 80% domain with t_fix = T_ad can force M to drift
    //      as T(z=80%) evolves during PT, creating tension that stalls convergence.
    let t_fix = t_unburned + 0.5 * (t_ad - t_unburned);
    let z_fix = z_center;

    // Initial mass flux guess: Su ≈ 0.4 m/s (methane-air), ρ_unburned
    let w_u = mean_molecular_weight(&mech.species, &y_unburned);
    let rho_u = density(p, t_unburned, w_u);
    let (su_guess, su_source) = match config.solver.su_initial_guess {
        Some(v) => (v, "user-specified"),
        None    => (0.4_f64, "default (methane-air)"),
    };
    println!("Initial flame speed guess: {su_guess:.3} m/s  ({su_source})");
    let m_guess = rho_u * su_guess;

    // Build initial guess: either from a Cantera CSV profile or from a sigmoid.
    let mut x = if let Some(ref csv_path) = config.solver.initial_profile {
        println!("Loading initial profile from: {csv_path}");
        initial_guess_from_csv(csv_path, &mech, &grid)?
    } else {
        initial_guess(
            &mech, &grid, t_unburned, t_ad,
            &y_unburned, &y_burned,
            m_guess, z_center, z_width,
        )
    };

    // --- Residual configuration ---
    let res_config = ResidualConfig {
        pressure: p,
        t_unburned,
        y_unburned: y_unburned.clone(),
        z_fix,
        t_fix,
    };

    // --- Phase 1: Pseudo-transient continuation ---
    // Skip PT when starting from an external initial profile (e.g. Cantera CSV):
    // the profile is already physically reasonable and close to the Newton basin,
    // so PT is counterproductive (it uses a coarse grid approximation and the
    // J_pt eigenvalue analysis doesn't apply to an already-converged profile).
    let pt_config = PseudoTransientConfig {
        n_steps: config.solver.time_steps,
        dt_initial: config.solver.dt_initial,
        dt_max: 1e-4,  // cap at 100 µs: keeps implicit steps ≤ ~10 K per cell
        ..Default::default()
    };
    if config.solver.initial_profile.is_none() {
        println!("--- Phase 1: Pseudo-transient continuation ---");
        pt_step(&mut x, &mech, &grid, &res_config, &pt_config)?;
    } else {
        println!("--- Phase 1: Skipped (starting from external initial profile) ---");
    }

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
        ratio: config.grid.ratio,
        max_points: config.grid.max_points,
    };
    let mut res_config = res_config;

    for refine_pass in 0..20 {
        println!("--- Newton solve, pass {refine_pass} ({} grid points) ---", grid.n_points());

        // Measure F before Newton so we can check whether it actually helped.
        let n = solution_length(&mech, grid.n_points());
        let mut f_buf = vec![0.0_f64; n];
        eval_residual(&x, &mut f_buf, &mech, &grid, &res_config, None, 0.0);
        let norm_before = f_buf.iter().map(|v| v * v).sum::<f64>().sqrt();

        newton_solve(&mut x, &mech, &grid, &res_config, &newton_config)?;

        eval_residual(&x, &mut f_buf, &mech, &grid, &res_config, None, 0.0);
        let norm_after = f_buf.iter().map(|v| v * v).sum::<f64>().sqrt();
        println!("  ‖F‖ after Newton: {norm_after:.3e} (was {norm_before:.3e})");

        // Update z_fix: find grid point where T first crosses t_fix.
        res_config.z_fix = find_z_fix(&x, &mech, &grid, res_config.t_fix);

        // Refinement gating:
        //   - Always refine if F < 1e6 (cleanly converged on this grid).
        //   - Also refine if Newton improved F by ≥10% and F is finite: Newton
        //     has reached the truncation-error floor for this grid and a finer
        //     grid will let it converge further.  This is the normal path when
        //     starting from an external (e.g. Cantera) initial profile.
        //   - Fall back to PT only when neither condition holds (Newton diverged
        //     or made no progress from a bad initial state).
        let converged = norm_after
            < newton_config.atol * (n as f64).sqrt()
                + newton_config.rtol * norm_after;
        let at_truncation_floor = norm_after.is_finite()
            && norm_after < norm_before * 0.9;
        let good_progress = norm_after < 1e6 || at_truncation_floor;

        if !good_progress {
            println!("  Newton made insufficient progress; running more PT before refining.");
            pt_step(&mut x, &mech, &grid, &res_config, &pt_config)?;
            // Don't refine this pass; try Newton again next iteration.
            continue;
        }
        let _ = converged;

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

/// Build reactant and product mass fractions from either:
///   - Equivalence-ratio mode: `fuel` + `oxidizer` + `equivalence_ratio` in config, or
///   - Direct composition mode: `composition` mole fractions in config.
fn compute_compositions(
    mech: &crate::chemistry::mechanism::Mechanism,
    config: &FlameConfig,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let nk = mech.n_species();

    // Select mode based on which config fields are present.
    let x_mol: Vec<f64> = if let Some(comp) = &config.flame.composition {
        // --- Direct composition mode ---
        let total: f64 = comp.values().sum();
        anyhow::ensure!(total > 0.0, "composition mole fractions sum to zero");
        let mut x = vec![0.0_f64; nk];
        for (name, &frac) in comp {
            if let Some(k) = mech.species_index(name) {
                x[k] = frac / total;
            }
        }
        x
    } else {
        // --- Equivalence-ratio mode ---
        let fuel = config.flame.fuel.as_ref()
            .ok_or_else(|| anyhow::anyhow!("Must specify either `composition` or `fuel`+`oxidizer`+`equivalence_ratio`"))?;
        let oxidizer = config.flame.oxidizer.as_ref()
            .ok_or_else(|| anyhow::anyhow!("Must specify `oxidizer` when using equivalence-ratio mode"))?;
        let phi = config.flame.equivalence_ratio
            .ok_or_else(|| anyhow::anyhow!("Must specify `equivalence_ratio` when using equivalence-ratio mode"))?;

        // Normalise fuel (sum = 1)
        let fuel_total: f64 = fuel.values().sum();
        anyhow::ensure!(fuel_total > 0.0, "Fuel mole fractions sum to zero");
        let mut x_fuel = vec![0.0_f64; nk];
        for (name, &frac) in fuel {
            if let Some(k) = mech.species_index(name) {
                x_fuel[k] = frac / fuel_total;
            }
        }

        // Normalise oxidizer (sum = 1)
        let ox_total: f64 = oxidizer.values().sum();
        anyhow::ensure!(ox_total > 0.0, "Oxidizer mole fractions sum to zero");
        let mut x_ox = vec![0.0_f64; nk];
        for (name, &frac) in oxidizer {
            if let Some(k) = mech.species_index(name) {
                x_ox[k] = frac / ox_total;
            }
        }

        // Stoichiometric O2 per mole of fuel mixture:
        //   O2_needed = C + H/4 - O/2 + S  (complete combustion to CO2 + H2O + SO2)
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

        let o2_idx = mech.species_index("O2");
        let x_o2_in_ox = o2_idx.map(|k| x_ox[k]).unwrap_or(0.0);
        anyhow::ensure!(x_o2_in_ox > 0.0, "Oxidizer contains no O2");

        let n_ox = if stoich_o2_per_fuel > 0.0 {
            stoich_o2_per_fuel / (x_o2_in_ox * phi)
        } else {
            1.0  // no combustible content — mix 1:1 as fallback
        };
        let total = 1.0 + n_ox;

        let mut x = vec![0.0_f64; nk];
        for k in 0..nk {
            x[k] = (x_fuel[k] + n_ox * x_ox[k]) / total;
        }
        x
    };

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

/// Seed the burned composition with equilibrium radicals.
///
/// Complete-combustion products (H₂O + N₂) contain no radicals.
/// Without radicals, chain branching cannot start: initiation reactions are
/// ~10¹² × slower than chain-branching, so pseudo-transient cannot ignite
/// the flame.  Cantera uses `equilibrate('HP')` to obtain the correct burned
/// state including minor species; we approximate this by partially
/// dissociating H₂O via element-conserving reactions:
///
///   H₂O → H  + OH            (Kp determines dissociation fraction)
///   2H₂O → 2H₂ + O₂         (reverse combustion, Kp-based)
///
/// The fractions are tuned to match the approximate equilibrium mole
/// fractions at T_ad ≈ 2500 K for H₂/air (φ = 1):
///   x_H ≈ 0.6%, x_OH ≈ 3.8%, x_H₂ ≈ 0.8%, x_O₂ ≈ 0.4%.
/// These values seed the chain-branching reactions and are orders of
/// magnitude more effective than initiation reactions at T ≈ 1000–2000 K.
fn seed_radicals(
    mech: &crate::chemistry::mechanism::Mechanism,
    mut y: Vec<f64>,
) -> Vec<f64> {
    // Find the relevant species indices; bail out if H₂O is absent.
    let h2o = match mech.species_index("H2O") { Some(k) => k, None => return y };
    let h   = mech.species_index("H");
    let oh  = mech.species_index("OH");
    let h2  = mech.species_index("H2");
    let o2  = mech.species_index("O2");

    let y_h2o = y[h2o];
    if y_h2o < 1e-6 { return y; }

    // Reaction 1: H₂O → H + OH   (dissociation fraction ~6%)
    //   W(H₂O)=18, W(H)=1, W(OH)=17; mass is conserved exactly.
    let f1 = 0.06_f64;
    let d1 = f1 * y_h2o;
    y[h2o] -= d1;
    if let Some(k) = h  { y[k] += d1 * (1.0 / 18.0); }
    if let Some(k) = oh { y[k] += d1 * (17.0 / 18.0); }

    // Reaction 2: 2H₂O → 2H₂ + O₂  (reverse combustion, fraction ~2%)
    //   Per 2 mol H₂O (36 g) → 2 mol H₂ (4 g) + 1 mol O₂ (32 g): mass conserved.
    let f2 = 0.02_f64;
    let d2 = f2 * y[h2o]; // use updated y[h2o]
    y[h2o] -= d2;
    if let Some(k) = h2 { y[k] += d2 * (4.0 / 36.0); }
    if let Some(k) = o2 { y[k] += d2 * (32.0 / 36.0); }

    // Re-normalise to ensure Σ Yk = 1 (rounding from discrete transfers).
    let sum: f64 = y.iter().sum();
    if sum > 0.0 { y.iter_mut().for_each(|v| *v /= sum); }
    y
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
/// Solves H_reactants(T_u) = H_products(T_ad) iteratively with Newton steps.
fn estimate_t_adiabatic(
    mech: &crate::chemistry::mechanism::Mechanism,
    y_unburned: &[f64],
    y_burned: &[f64],
    t_u: f64,
    _p: f64,
) -> f64 {
    use crate::chemistry::thermo::{cp_mixture, enthalpy_molar};

    // Specific enthalpy of reactants at T_u [J/kg]
    let h_u: f64 = mech.species.iter().zip(y_unburned.iter())
        .map(|(s, &yk)| yk * enthalpy_molar(s, t_u) / s.molecular_weight)
        .sum();

    // Newton iteration: find T_ad such that h_products(T_ad) = h_u
    let mut t_ad = t_u + 1500.0;
    for _ in 0..200 {
        let h_ad: f64 = mech.species.iter().zip(y_burned.iter())
            .map(|(s, &yk)| yk * enthalpy_molar(s, t_ad) / s.molecular_weight)
            .sum();
        let cp: f64 = cp_mixture(&mech.species, y_burned, t_ad);
        let dt = (h_u - h_ad) / cp.max(1.0);
        t_ad += dt;
        if dt.abs() < 0.1 { break; }
    }
    t_ad.clamp(t_u + 100.0, 6000.0)
}
