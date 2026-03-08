# 1D Freely Propagating Premixed Laminar Flame Solver — Task Description & TODOs

## Project Description

Implement a 1D steady-state freely propagating premixed laminar flame solver in Rust.
Given a fuel-oxidizer mixture at unburned conditions (temperature, pressure, equivalence ratio),
the code computes the full spatial profiles of temperature and species mass fractions across the
flame, as well as the laminar flame speed (eigenvalue of the problem).

The solver follows the classical PREMIX/Cantera approach:
- **Governing equations**: species conservation + energy equation on a non-uniform 1D grid,
  coupled to detailed chemical kinetics and mixture-averaged transport.
- **Discretization**: finite differences (upwind convection, central diffusion).
- **Solver**: damped Newton iteration with pseudo-transient continuation for robustness,
  plus adaptive grid refinement (GRAD/CURV criteria).
- **Input**: Cantera YAML or CHEMKIN-II mechanism files; TOML config.
- **Output**: CSV profile of z, T, u, ρ, Yk for all species.

Reference implementations: PREMIX (Kee et al., Sandia 1985) and Cantera `FreeFlame`.

---

## TODOs

### Phase 1 — Chemistry (prerequisite for everything)

- [ ] **`chemistry/parser/cantera_yaml.rs`** — complete Cantera YAML parser
  - [ ] Fix `molecular_weight_from_composition`: handle case-insensitive element names (e.g. `N` vs `n`)
  - [ ] Parse activation energy units from the `units:` block at the top of the YAML file (currently hardcoded to cal/mol)
  - [ ] Handle `duplicate` reaction flag
  - [ ] Handle `PLOG` pressure-dependent reactions
  - [ ] Test against `data/h2o2.yaml` (copy from cantera reference repo)

- [ ] **`chemistry/parser/chemkin.rs`** — implement CHEMKIN-II parser
  - [ ] Parse `ELEMENTS … END` block
  - [ ] Parse `SPECIES … END` block
  - [ ] Parse `therm.dat`: NASA 7-coeff blocks (fixed-column format)
  - [ ] Parse `tran.dat`: one line per species (geometry, ε/k, σ, μ, α, Zrot)
  - [ ] Parse `REACTIONS … END`: Arrhenius A/b/Ea, `LOW/`, `TROE/`, `REV/`, `DUPLICATE`

- [ ] **`chemistry/thermo.rs`** — unit tests
  - [ ] Verify cp(T), h(T), s(T) against tabulated JANAF values for H2, O2, H2O, N2

- [ ] **`chemistry/kinetics.rs`** — unit tests
  - [ ] Verify Kc and ωk for a simple reaction against known values
  - [ ] Test Troe falloff broadening against reference tables

### Phase 2 — Transport

- [ ] **`transport/collision_integrals.rs`**
  - [ ] Unit test `omega11` and `omega22` against Neufeld (1972) table values

- [ ] **`transport/species_props.rs`**
  - [ ] Verify μk for N2 at 300 K and 1000 K against JANAF/literature
  - [ ] Verify Dij(H2, N2) against known values

- [ ] **`transport/mixture.rs`**
  - [ ] Verify Wilke mixture viscosity for air (O2/N2) against experiment

### Phase 3 — Flame Residual (core correctness)

- [ ] **`flame/residual.rs`** — fix and complete
  - [ ] Fix left-boundary species BC: use proper inlet flux form
    `F_yk = M·(Yk - Yku) + jk_{1/2}` (not just `Yk - Yku`)
  - [ ] Fix enthalpy transport term: remove erroneous `MinOr` trait hack;
    use `jk_mid[k][j-1]` directly for the average
  - [ ] Remove unused `sum_y` variable in left BC
  - [ ] Remove unused `x_j` variable
  - [ ] Add time-derivative term `rdt * (x - x_old)` support for pseudo-transient embedding
  - [ ] Unit test: verify residual is ≈ 0 for a known analytic steady profile

- [ ] **`flame/state.rs`**
  - [ ] Implement a proper equilibrium estimate for `y_burned` in `initial_guess`
    (currently a placeholder returning reactant composition)

- [ ] **`flame/solver_driver.rs`**
  - [ ] Implement `compute_compositions`: proper mole-fraction-to-mass-fraction conversion
    respecting equivalence ratio (fuel/oxidizer stoichiometry)
  - [ ] Implement `estimate_burned_composition`: element-balance complete combustion estimate
    for CO2, H2O, N2 products
  - [ ] Pass `ResidualConfig` `z_fix` dynamically (track where T crosses `t_fix` after each Newton pass)

### Phase 4 — Solver

- [ ] **`solver/banded.rs`**
  - [ ] Replace the dense LU fallback with a true banded LU factorization
    (LAPACK `dgbsv` via `lapack` crate, or a hand-rolled band LU)
    to avoid O(n³) cost for large grids
  - [ ] Add a test: solve a small banded system and verify against known solution

- [ ] **`solver/newton.rs`**
  - [ ] Implement Jacobian reuse (age-based): only recompute after `max_jac_age` steps
  - [ ] Add convergence criterion based on scaled residual norm (not just absolute)
  - [ ] Test on a small nonlinear system (e.g., 2×2 Newton problem)

- [ ] **`solver/pseudo_transient.rs`**
  - [ ] Pass `x_old` correctly for time derivative: `(x - x_old)/dt` residual contribution
  - [ ] Test standalone on a simple ODE

### Phase 5 — Grid Refinement & IO

- [ ] **`flame/refine.rs`**
  - [ ] Fix `new_x` length management: ensure mass flux `M` is always the last element
    after splice insertions (currently may be misaligned)
  - [ ] Add a test: refine a simple profile and verify interpolation accuracy

- [ ] **`io/input.rs`**
  - [ ] Validate config on load (e.g., equivalence_ratio > 0, domain_length > 0)
  - [ ] Support reading unburned composition directly as mass/mole fractions (bypass φ)

- [ ] **`io/output.rs`**
  - [ ] Add heat release rate column: `HRR = -Σk ωk·hk·Wk` [W/m³]
  - [ ] Add mole fraction columns in addition to mass fractions

### Phase 6 — Integration Tests & Validation

- [ ] Copy `data/h2o2.yaml` from the Cantera reference repo into the `data/` directory
- [ ] Run H2/air flame (φ=1, 300 K, 1 atm) end-to-end:
  - [ ] Target: Su ≈ 2.1 m/s (compare against Cantera `FreeFlame`)
  - [ ] Compare T profile and species profiles (H2, O2, H2O, OH) against Cantera output
- [ ] Run CH4/air flame (φ=1, 300 K, 1 atm) with GRI-Mech 3.0:
  - [ ] Target: Su ≈ 0.37 m/s
- [ ] Profile performance on GRI (53 species, ~320 reactions): identify bottlenecks

### Phase 7 — Optional Enhancements

- [ ] Multicomponent transport (full diffusion matrix Dkj, Soret effect)
- [ ] CHEMKIN-II `.inp` + `.dat` parser completion
- [ ] Sensitivity analysis: ∂Su/∂Ai for reaction rate parameters
- [ ] Radiation heat loss model (optically thin, CO2/H2O)
- [ ] Parallel Jacobian evaluation (rayon)
