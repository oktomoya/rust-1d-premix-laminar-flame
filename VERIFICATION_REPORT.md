# Cantera YAML Parser Verification Report

**Date:** 2026-03-08
**Checked by:** Cross-referencing against Cantera source in `/Users/okuno/Projects/combustion/cantera`
**Mechanism file:** `data/h2o2.yaml` (H2/O2 sub-mechanism, 10 species, 29 reactions)

---

## 1. What Was Checked

| Item | Source consulted |
|------|-----------------|
| PLOG interpolation formula | `cantera/src/kinetics/PlogRate.cpp` |
| Troe falloff formula (Fcent, F broadening) | `cantera/src/kinetics/Falloff.cpp` – `TroeRate` class |
| Energy unit conversion table | `cantera/src/base/Units.cpp` |
| Reaction/falloff data parsing | `cantera/src/kinetics/Reaction.cpp`, `Falloff.cpp` |
| Molecular weight calculation (g/mol vs kg/mol) | Our `molecular_weight_from_composition` vs h2o2.yaml |
| Species/reaction count in h2o2.yaml | Manual count of YAML |
| Duplicate flag parsing | h2o2.yaml reactions 24–29 |
| Third-body efficiency parsing | h2o2.yaml reactions 1, 2, 6, 12, 15, 22 |
| Stoichiometric coefficient parsing | All 29 reactions |

### Cantera Python

`python3 -c 'import cantera as ct'` → **not installed**.
The script `scripts/verify_cantera.py` was written and will run when cantera is available.

---

## 2. Bugs Found and Fixed

### Bug 1 — `is_badvalue()` does not exist in serde_yaml (compilation error)

**File:** `src/chemistry/parser/cantera_yaml.rs`, lines 210 and 224

**Problem:** The parser called `rate_node.is_badvalue()` and `rxn_node["Troe"].is_badvalue()`.
`serde_yaml::Value` has no such method; it is a Cantera C++ AnyMap concept.
In serde_yaml, indexing a mapping with a missing key returns `Value::Null`, so
`is_null()` is sufficient and the code did not compile.

**Before (line 210):**
```rust
if rate_node.is_null() || rate_node.is_badvalue() {
```

**Before (line 224):**
```rust
let troe = if !rxn_node["Troe"].is_null() && !rxn_node["Troe"].is_badvalue() {
```

**After (line 210):**
```rust
if rate_node.is_null() {
```

**After (line 224):**
```rust
let troe = if !rxn_node["Troe"].is_null() {
```

---

### Bug 2 — Falloff reactions double-counted third-body concentration M (logic error)

**File:** `src/chemistry/kinetics.rs`, `production_rates` function

**Problem:** For falloff reactions, `falloff_rate()` already incorporates `[M]`
through the reduced pressure `Pr = k₀ [M] / k∞`. The code then unconditionally
multiplied the returned rate by `m_conc` a second time for any reaction that has
`third_body.is_some()`, inflating the falloff rate by a factor of `[M]`.

Reactions 22 in h2o2.yaml (`2 OH (+M) <=> H2O2 (+M)`) was affected.

**Before:**
```rust
let m_conc = third_body_concentration(mech, i, concentrations);
let kf_eff = if mech.reactions[i].third_body.is_some() {
    kf * m_conc
} else {
    kf
};
```

**After:**
```rust
// Falloff reactions already incorporate [M] inside falloff_rate() via Pr,
// so they must NOT be multiplied by m_conc again.
// Only plain three-body Arrhenius reactions need the extra [M] factor.
let kf_eff = match &rxn.rate {
    crate::chemistry::mechanism::RateType::Falloff { .. } => kf,
    _ => {
        if mech.reactions[i].third_body.is_some() {
            let m_conc = third_body_concentration(mech, i, concentrations);
            kf * m_conc
        } else {
            kf
        }
    }
};
```

---

## 3. Verification Results

### a. Molecular weight calculation

- Units are **kg/mol** throughout (elements in the table are stored as `e-3` fractions).
- Element lookup is case-insensitive: keys are uppercased before hash-map lookup.
- `{Ar: 1}` composition → `"Ar".to_uppercase()` = `"AR"` → 39.948e-3 kg/mol ✓
- Test `test_molecular_weights` checks H2, O2, AR, N2. **PASS**

Reference values from h2o2.yaml:

| Species | Composition | Expected MW (kg/mol) | Our result |
|---------|-------------|----------------------|-----------|
| H2  | {H:2}       | 2.016e-3  | ✓ |
| O2  | {O:2}       | 31.998e-3 | ✓ |
| AR  | {Ar:1}      | 39.948e-3 | ✓ |
| N2  | {N:2}       | 28.014e-3 | ✓ |

### b. Activation energy unit conversion

Global units block in h2o2.yaml: `activation-energy: cal/mol`

Conversion factor verified against `cantera/src/base/Units.cpp`:
```
{"cal", Units(4.184, 1, 2, -2)}   →  1 cal = 4.184 J
```

Our conversion: `ea_raw * 4.184` (cal/mol → J/mol). ✓
Full table implemented:

| Unit    | Factor to J/mol |
|---------|-----------------|
| cal/mol | × 4.184 |
| kcal/mol| × 4184.0 |
| J/mol   | × 1.0 |
| kJ/mol  | × 1000.0 |
| K       | × 8.314462618 (= R) |
| J/kmol  | × 1e-3 |

Test `test_activation_energy_units_cal_mol` verifies reaction 3 (Ea=6260 cal/mol → 26191.84 J/mol). **PASS**

### c. PLOG interpolation formula

Cantera (`PlogRate.cpp`) uses log-linear interpolation:
```
log k = log k₁ + (log P - log P₁) * (log k₂ - log k₁) / (log P₂ - log P₁)
```

Our implementation in `kinetics.rs::plog_rate`:
```rust
let log_k = k1.ln() + (pressure.ln() - p1.ln()) * (k2.ln() - k1.ln()) / (p2.ln() - p1.ln());
log_k.exp()
```
Matches exactly. ✓
Pressure field in PLOG entries: parsed with unit conversion (atm, bar, kPa, MPa, Pa). ✓
Test `test_plog_parsing` verifies two-entry PLOG reaction. **PASS**

### d. Troe falloff parameters

Cantera YAML format: `Troe: {A: ..., T3: ..., T1: ..., T2: ...}`

Cantera's `TroeRate::setParameters` reads them in exactly this order (A, T3, T1, T2).

Our `TroeParams` struct: `a`, `t3`, `t1`, `t2` — matching order. ✓

Cantera's `Fcent` formula (`Falloff.cpp`, `TroeRate::updateTemp`):
```cpp
Fcent = (1 - A)*exp(-T/T3) + A*exp(-T/T1) + exp(-T2/T)
```

Our Rust:
```rust
let f_cent = (1.0 - troe.a) * (-t / troe.t3).exp()
    + troe.a * (-t / troe.t1).exp()
    + troe.t2.map(|t2| (-t2 / t).exp()).unwrap_or(0.0);
```
Matches. ✓

Cantera's `TroeRate::F` broadening:
```cpp
double cc = -0.4 - 0.67 * log10(Fcent);
double nn = 0.75 - 1.27 * log10(Fcent);
double f1 = (lpr + cc) / (nn - 0.14*(lpr + cc));
return pow(10, log10(Fcent) / (1 + f1*f1));
```

Our `troe_broadening` uses `c`, `n`, `d=0.14` with identical expressions. ✓

Test `test_falloff_reaction_22` verifies high.A=7.4e13, low.A=2.3e18, Troe.A=0.7346. **PASS**

### e. Third-body efficiency parsing

Format in h2o2.yaml: `efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}`

Our code stores efficiencies as `Vec<(usize, f64)>` (species index, enhancement factor).
The default efficiency of 1.0 is implicit; only non-default species are stored.
When computing `[M]_eff`:
```rust
let mut m = total_concentration;
for &(k, eff) in &tb.efficiencies {
    m += (eff - 1.0) * concentrations[k];
}
```
This correctly applies per-species enhancement above/below unity. ✓

Test `test_three_body_reaction_1` checks reaction 1 has non-empty efficiencies. **PASS**

### f. Stoichiometric coefficient parsing

Handles:
- `2 OH` (space-separated integer coefficient)
- `2OH` (coefficient directly adjacent to name)
- Implicit coefficient of 1 when no leading digits

Verified for h2o2.yaml reactions with `2 O`, `2 H`, `2 OH`, `H + 2 O2`. ✓

Third-body marker `(+M)` stripped before parsing; standalone `M` tokens skipped. ✓

### g. Duplicate flag handling

h2o2.yaml has 6 duplicate reactions: reactions 24–29 (0-indexed: 23–28).
Parser reads `rxn_node["duplicate"].as_bool().unwrap_or(false)`.

Test `test_duplicate_flag` asserts reactions[23..28].duplicate == true and reactions[0].duplicate == false. **PASS**

### h. Reaction count for h2o2.yaml

Manual count of h2o2.yaml reactions section: **29 reactions** (as labelled # Reaction 1 … # Reaction 29).

Test `test_parse_h2o2_reaction_count` asserts `mech.n_reactions() == 29`. **PASS**

---

## 4. cargo test Output

```
running 10 tests
test chemistry::parser::cantera_yaml::tests::test_case_insensitive_elements ... ok
test transport::collision_integrals::tests::omega22_known_value ... ok
test chemistry::parser::cantera_yaml::tests::test_plog_parsing ... ok
test chemistry::parser::cantera_yaml::tests::test_parse_h2o2_species_count ... ok
test chemistry::parser::cantera_yaml::tests::test_three_body_reaction_1 ... ok
test chemistry::parser::cantera_yaml::tests::test_falloff_reaction_22 ... ok
test chemistry::parser::cantera_yaml::tests::test_molecular_weights ... ok
test chemistry::parser::cantera_yaml::tests::test_activation_energy_units_cal_mol ... ok
test chemistry::parser::cantera_yaml::tests::test_duplicate_flag ... ok
test chemistry::parser::cantera_yaml::tests::test_parse_h2o2_reaction_count ... ok

test result: ok. 10 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.02s
```

All 10 tests pass.

---

## 5. Remaining Concerns

1. **No cantera Python installed** — `scripts/verify_cantera.py` was written but could not be run. Numerical cross-validation of NASA7 coefficients, Arrhenius A/b/Ea values, and Troe parameters against cantera's Python API is deferred until `pip install cantera` is available.

2. **`entropy_species` naming** (`thermo.rs`) — The function is named `entropy_species` and documented as J/(kg·K), but its caller in `kinetics.rs::equilibrium_constant` multiplies by `molecular_weight` to recover J/(mol·K). This is correct but fragile: if a future caller does not perform the conversion the units will be wrong. Adding a `entropy_molar` counterpart would reduce the risk.

3. **Specific-third-body reactions** (reactions 7–10, 13–14) — These have explicit collision partners (O2, H2O, N2, AR) that appear on both sides of the equation. Our parser correctly treats them as plain `Arrhenius` reactions with no `ThirdBodySpec`. This matches Cantera's treatment of "explicit third body" partners, which are handled as ordinary stoichiometric participants.

4. **`+ M` detection** — The check `eq.contains("+ M")` is redundant for h2o2.yaml since all `M`-style three-body reactions also carry `type: three-body`. The check would, however, catch mechanisms that omit the `type` field but include ` + M` in the equation string. It is harmless for h2o2.yaml.

5. **Chemically-activated reactions** — Cantera's `FalloffRate` supports a `chemically-activated` variant where the low-P limit is used differently. Our parser does not distinguish this case; it would be parsed as a plain falloff. No such reactions exist in h2o2.yaml.
