# Design: 1D Freely Propagating Premixed Laminar Flame Solver in Rust

## 1. Problem Statement

Compute the temperature profile, species mass fraction profiles, and laminar flame speed of a
1D steady-state freely propagating premixed flame. The code solves the full species transport
and energy equations coupled to detailed chemical kinetics.

Reference implementations:
- Cantera `FreeFlame` — `cantera/src/onD/Flow1D.cpp`

---

## 2. Governing Equations

The steady-state 1D flame is described in lab coordinates with spatial variable z [m].
For a **freely propagating** adiabatic flame the mass flux M = ρu is constant (the eigenvalue).

### 2.1 Notation

| Symbol | Meaning | Units (SI) |e
|--------|---------|------------|
| z      | spatial coordinate | m |
| T      | temperature | K |
| Yk     | mass fraction of species k | — |
| M      | mass flux = ρu (eigenvalue) | kg/(m²·s) |
| ρ      | mixture density | kg/m³ |
| u      | velocity | m/s |
| λ      | thermal conductivity | W/(m·K) |
| cp     | mixture specific heat | J/(kg·K) |
| cpk    | species specific heat | J/(kg·K) |
| hk     | species specific enthalpy | J/kg |
| Wk     | species molecular weight | kg/mol |
| W      | mixture mean molecular weight | kg/mol |
| jk     | diffusive mass flux of species k | kg/(m²·s) |
| ωk     | molar production rate of species k | mol/(m³·s) |
| Dk     | mixture-averaged diffusion coeff | m²/s |
| P      | pressure (constant) | Pa |

### 2.2 Species Conservation (k = 1 … K)

```
M · dYk/dz  =  -d(jk)/dz  +  ωk · Wk
```

Diffusive flux using mixture-averaged model with correction velocity:

```
jk = -ρ · Dk · dYk/dz  +  Yk · Vc · ρ
Vc = Σk ρ·Dk·dYk/dz / ρ           (correction velocity, ensures Σ jk = 0)
```

Alternatively, in PREMIX's formulation (which we follow):
`jk = ρ · Yk · Vk`  where Vk is the diffusion velocity.

### 2.3 Energy Equation

```
M · cp · dT/dz  =  d(λ · dT/dz)/dz  -  Σk jk · cpk · dT/dz  -  Σk ωk · hk · Wk
```

(The last term `Σk ωk · hk · Wk` is the heat release rate; note hk here is J/mol.)

### 2.4 Continuity

```
dM/dz = 0   →   M = const   (unknown flame eigenvalue)
```

### 2.5 Equation of State

```
P = ρ · (R / W) · T        (ideal gas)
ρ = P · W / (R · T)
```

### 2.6 Boundary Conditions

**Left boundary (z = 0, unburned reactants):**
- `T(0) = Tu`  (fresh gas temperature)
- `Yk(0) = Yku`  for k = 1 … K  (reactant composition)
- Species BC enforced via inlet flux: `M·Yk - jk = M·Yku`

**Right boundary (z = L, burned products):**
- Zero-gradient: `dT/dz = 0`, `dYk/dz = 0`

**Flame location (eigenvalue closure):**
- Fix temperature at one interior point: `T(z_fix) = T_fix`
  (typically T_fix ≈ Tu + 200 K; this prevents translation of the solution)

---

## 3. Numerical Method

### 3.1 Spatial Discretization

Non-uniform 1D grid: `z_1 < z_2 < … < z_J`

Variables at each grid point j: `(T_j, Y_{1,j}, … Y_{K,j})`  — total `(K+1)·J` unknowns.
The mass flux M is a global scalar unknown → total unknowns: `(K+1)·J + 1`.

**Finite difference stencils** (following PREMIX/Cantera):

Convection term at j (upwind):
```
M · dYk/dz |_j  ≈  M · (Yk_j - Yk_{j-1}) / (z_j - z_{j-1})
```

Diffusion term at j (central, midpoint conductances):
```
d(jk)/dz |_j  ≈  (jk_{j+1/2} - jk_{j-1/2}) / ((z_{j+1} - z_{j-1})/2)
```
where midpoint fluxes use arithmetic averages of adjacent point properties.

Conduction term at j:
```
d(λ dT/dz)/dz |_j  ≈  (λ_{j+1/2}·(T_{j+1}-T_j)/dz_p  -  λ_{j-1/2}·(T_j-T_{j-1})/dz_m)
                        / ((z_{j+1} - z_{j-1})/2)
```

### 3.2 Newton Solver (Damped Newton + Pseudo-Transient Continuation)

The discretized system is `F(x) = 0` where x is the full solution vector.

**Algorithm:**

1. Start with initial guess from a simple linear profile.
2. **Pseudo-transient continuation** (time-stepping): add `(x - x_old)/Δt` to residual, solve for several time steps until solution is close to steady state.
3. Switch to **damped Newton** iteration:
   - Assemble Jacobian J = ∂F/∂x (numerically via finite differences, banded structure)
   - Solve J·Δx = -F  using banded LU factorization
   - Line search: find step size λ ∈ (0,1] such that ‖F(x + λΔx)‖ < ‖F(x)‖
   - Update x ← x + λΔx
4. Convergence: ‖F(x)‖ < tol

**Jacobian structure:** The system has a block-banded structure with bandwidth proportional to
(K+1) (number of variables per point). Using LAPACK's `dgbsv` or a custom banded solver.

### 3.3 Adaptive Grid Refinement (GRAD/CURV)

After each Newton solve, refine the grid by inserting points where:
- **GRAD**: `|φ_{j+1} - φ_j| / max|φ| > GRAD` (large gradient)
- **CURV**: `|(φ_{j+1} - φ_j)/dz_{j+1/2} - (φ_j - φ_{j-1})/dz_{j-1/2}| / max|dφ/dz| > CURV`

(where φ ranges over T and significant species Yk)

Refinement adds a midpoint between j and j+1; solution is interpolated linearly.
Repeat solve → refine until no more points are added (or max grid points reached).

### 3.4 Solution Strategy

```
1. Solve with fixed T profile (energy equation off) → get species profiles
2. Turn on energy equation → full solve
3. Refine grid, resolve
4. Repeat until converged
```

---

## 4. Thermodynamic Properties (NASA 7-Coefficient Polynomials)

For each species k, two temperature ranges [T_low, T_mid] and [T_mid, T_high]:

```
cp_k(T)/R = a1 + a2·T + a3·T² + a4·T³ + a5·T⁴
h_k(T)/RT = a1 + a2/2·T + a3/3·T² + a4/4·T³ + a5/5·T⁴ + a6/T
s_k(T)/R  = a1·ln(T) + a2·T + a3/2·T² + a4/3·T³ + a5/4·T⁴ + a7
```

Mixture properties:
```
cp_mix = Σk Yk · cpk
h_mix  = Σk Yk · hk
```

---

## 5. Chemical Kinetics

Elementary reactions with Arrhenius rate coefficients:

```
kf = A · T^b · exp(-Ea / (R·T))
```

For each reaction i, forward rate `kf_i` and reverse rate `kr_i = kf_i / Kc_i`.
Species production rates:
```
ωk = Σi νki · (kf_i · Πreactants [Cj]^νji  -  kr_i · Πproducts [Cj]^νji)
```

Support: third-body reactions, pressure-dependent (Lindemann/Troe) falloff.

---

## 6. Transport Properties (Mixture-Averaged)

**Species viscosity** (Chapman-Enskog kinetic theory):
```
μk = 2.6693e-6 · sqrt(Wk · T) / (σk² · Ω^(2,2)*(T*))
```

**Species thermal conductivity:**
```
λk = μk / Wk · (R · (f_trans·cv_trans + f_rot·cv_rot + f_vib·cv_vib))
```
(Eucken modified formula)

**Mixture viscosity** (Wilke mixing rule):
```
μ_mix = Σk Xk·μk / Σj Xj·Φkj
```

**Mixture thermal conductivity** (arithmetic mean approximation):
```
λ_mix = 0.5 · (Σk Xk·λk  +  1 / Σk Xk/λk)
```

**Species diffusion coefficients** (binary diffusion + mixture-averaged):
```
Dkm = (1 - Xk) / Σ_{j≠k} Xj / Dkj
```

Transport parameters (ε/kb, σ, μ_dipole, α_polarizability, Z_rot) from CHEMKIN/Cantera
transport database (`.dat` or `.yaml` format).

---

## 7. Software Architecture

```
rust-1d-premix-laminar-flame/
├── Cargo.toml
├── DESIGN.md                      ← this file
├── src/
│   ├── main.rs                    ← CLI entry point: parse args, run flame
│   ├── lib.rs                     ← public API
│   │
│   ├── chemistry/
│   │   ├── mod.rs
│   │   ├── species.rs             ← Species struct (name, Wk, NASA coeffs, transport params)
│   │   ├── thermo.rs              ← NASA polynomial eval: cp, h, s, g
│   │   ├── kinetics.rs            ← reaction rates: Arrhenius, falloff, equilibrium constants
│   │   ├── mechanism.rs           ← Mechanism struct holding all species + reactions
│   │   └── parser/
│   │       ├── mod.rs
│   │       ├── chemkin.rs         ← Parse CHEMKIN-II .inp + .dat format
│   │       └── cantera_yaml.rs    ← Parse Cantera .yaml mechanism format
│   │
│   ├── transport/
│   │   ├── mod.rs
│   │   ├── collision_integrals.rs ← Omega*(T*) polynomial fits
│   │   ├── species_props.rs       ← μk, λk, Dkj for species pairs
│   │   └── mixture.rs             ← μ_mix, λ_mix, Dkm_mix
│   │
│   ├── flame/
│   │   ├── mod.rs
│   │   ├── domain.rs              ← Grid: z[], dz[], midpoint arrays; insert/refine
│   │   ├── state.rs               ← Solution vector x: T[], Yk[][], M; indexing helpers
│   │   ├── residual.rs            ← F(x): species + energy + continuity + BCs
│   │   ├── jacobian.rs            ← Finite-difference Jacobian (banded)
│   │   └── refine.rs              ← GRAD/CURV refinement criterion; interpolation
│   │
│   ├── solver/
│   │   ├── mod.rs
│   │   ├── newton.rs              ← Damped Newton with line search
│   │   ├── pseudo_transient.rs    ← Time-stepping (pseudo-transient continuation)
│   │   └── banded.rs              ← Banded LU factorization and solve
│   │
│   └── io/
│       ├── mod.rs
│       ├── input.rs               ← Read flame.toml config file
│       └── output.rs              ← Write CSV: z, T, u, Yk, HRR
│
├── examples/
│   ├── methane_air.toml           ← Example input: CH4/air, φ=1, 1 atm
│   └── hydrogen_air.toml          ← Example input: H2/air, φ=1, 1 atm
│
└── data/
    └── h2o2.yaml                  ← Cantera mechanism (copy from reference)
```

---

## 8. Data Structures

### 8.1 Mechanism

```rust
pub struct Species {
    pub name: String,
    pub molecular_weight: f64,          // kg/mol
    pub nasa_low: [f64; 7],             // coefficients for T < T_mid
    pub nasa_high: [f64; 7],            // coefficients for T >= T_mid
    pub nasa_t_mid: f64,
    pub transport: TransportParams,
}

pub struct TransportParams {
    pub geometry: GeometryType,         // Atom, Linear, Nonlinear
    pub well_depth: f64,                // ε/kb [K]
    pub diameter: f64,                  // σ [Angstrom]
    pub dipole_moment: f64,             // μ [Debye]
    pub polarizability: f64,            // α [Angstrom^3]
    pub rot_relax: f64,                 // Z_rot
}

pub struct Reaction {
    pub reactants: Vec<(usize, f64)>,   // (species_index, stoich_coeff)
    pub products: Vec<(usize, f64)>,
    pub rate: RateType,
    pub third_body: Option<ThirdBodySpec>,
}

pub enum RateType {
    Arrhenius { a: f64, b: f64, ea: f64 },
    Falloff { high: Arrhenius, low: Arrhenius, troe: Option<TroeParams> },
}

pub struct Mechanism {
    pub species: Vec<Species>,
    pub reactions: Vec<Reaction>,
}
```

### 8.2 Flame State

```rust
pub struct Grid {
    pub z: Vec<f64>,          // grid point positions [m]
}

pub struct FlameState {
    pub grid: Grid,
    pub temperature: Vec<f64>,          // T[j], length = npoints
    pub mass_fractions: Vec<Vec<f64>>,  // Y[k][j]
    pub mass_flux: f64,                 // M = ρu [kg/(m²·s)], the eigenvalue
    pub pressure: f64,                  // [Pa]
}
```

Solution vector x (flat, column-major by point):
```
x = [T_1, Y_{1,1}, …, Y_{K,1},   T_2, Y_{1,2}, …, Y_{K,2},   …,   T_J, …, Y_{K,J},   M]
```
Length: `(K+1)·J + 1`

### 8.3 Residual

```rust
pub struct FlameResidual<'a> {
    pub mech: &'a Mechanism,
    pub config: &'a FlameConfig,
    pub z_fixed: f64,       // position where T is pinned
    pub t_fixed: f64,       // pinned temperature
}

impl FlameResidual<'_> {
    pub fn eval(&self, state: &FlameState, rhs: &mut Vec<f64>);
    pub fn jacobian_banded(&self, state: &FlameState) -> BandedMatrix;
}
```

---

## 9. Input File Format (TOML)

```toml
# examples/methane_air.toml

[mechanism]
file = "data/gri30.yaml"   # Cantera YAML or CHEMKIN chem.inp+therm.dat+tran.dat
format = "cantera_yaml"    # or "chemkin"

[flame]
pressure = 101325.0        # Pa
fuel = { CH4 = 1.0 }
oxidizer = { O2 = 0.21, N2 = 0.79 }
equivalence_ratio = 1.0
T_unburned = 300.0         # K
domain_length = 0.02       # m

[grid]
initial_points = 20
max_points = 500
grad = 0.05
curv = 0.10

[solver]
atol = 1.0e-9
rtol = 1.0e-6
max_newton_iter = 50
time_steps = 100
dt_initial = 1.0e-7        # s

[output]
file = "flame_solution.csv"
```

---

## 10. Output

CSV with columns:
```
z [m], T [K], u [m/s], rho [kg/m³], HRR [W/m³], Y_CH4, Y_O2, Y_CO2, Y_H2O, ...
```

Scalar summary printed to stdout:
```
Laminar flame speed Su = X.XXX m/s
Adiabatic flame temperature Tad = XXXX K
Grid points: XXX
```

---

## 11. Implementation Plan

### Phase 1 — Chemistry (no transport, no solver)
1. `chemistry/species.rs` — Species struct, NASA poly evaluation
2. `chemistry/kinetics.rs` — Arrhenius rate, equilibrium constant, production rates
3. `chemistry/parser/chemkin.rs` — Parse `chem.inp` (mechanism) and `therm.dat`
4. Test: reproduce SENKIN equilibrium calculation

### Phase 2 — Transport
5. `transport/collision_integrals.rs` — Omega fits (Neufeld polynomials)
6. `transport/species_props.rs` — μk, λk, Dkj
7. `transport/mixture.rs` — mixture-averaged μ, λ, Dk
8. Parse `tran.dat` (CHEMKIN transport)

### Phase 3 — Flame Domain and Residual
9. `flame/domain.rs` — Grid struct, midpoint helpers
10. `flame/state.rs` — Solution vector layout, indexing
11. `flame/residual.rs` — F(x): species + energy BVP + BCs
12. Unit test: verify residual is zero for a known analytic profile

### Phase 4 — Solver
13. `solver/banded.rs` — Banded matrix storage (LU, solve)
14. `flame/jacobian.rs` — Finite-difference Jacobian
15. `solver/pseudo_transient.rs` — Time-stepping loop
16. `solver/newton.rs` — Damped Newton with line search
17. Integration test: solve H2/O2 flame (small mechanism, fast)

### Phase 5 — Grid Refinement and IO
18. `flame/refine.rs` — GRAD/CURV refinement, interpolation
19. `io/input.rs` — TOML config parser
20. `io/output.rs` — CSV writer
21. `main.rs` — CLI glue

### Phase 6 — Validation
22. Compare against Cantera `FreeFlame` for:
    - H2/air, φ=1.0, 1 atm → Su ≈ 2.1 m/s
    - CH4/air, φ=1.0, 1 atm → Su ≈ 0.37 m/s
23. Profile and optimize hot spots

---

## 12. Key External Crates

| Crate | Purpose |
|-------|---------|
| `serde` + `toml` | Input file parsing |
| `serde_yaml` | Cantera YAML mechanism parsing |
| `ndarray` | Multi-dimensional arrays (species matrix Y[k][j]) |
| `faer` | Dense/banded linear algebra (LU solve for Newton step) |
| `anyhow` | Error handling |
| `clap` | CLI argument parsing |
| `csv` | Output writing |

---

## 13. Validation References

- Smooke (1991), "Reduced Kinetic Mechanisms and Asymptotic Approximations for Methane-Air Flames"
- Kee, Grcar, Smooke, Miller (1985), "PREMIX: A Fortran Program for Modeling Steady Laminar One-Dimensional Premixed Flames", Sandia Report SAND85-8240
- Kee, Rupley, Miller (1989), "Chemkin-II: A Fortran Chemical Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics", Sandia Report SAND89-8009
- Cantera documentation: https://cantera.org/stable/userguide/flames.html
