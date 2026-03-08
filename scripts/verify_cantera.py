"""
Cantera YAML parser verification script.

Loads data/h2o2.yaml via the cantera Python package and prints reference values
that can be cross-checked against the Rust parser output.

Usage:
    python3 scripts/verify_cantera.py
"""

import sys
try:
    import cantera as ct
except ImportError:
    print("cantera not installed. Install with: pip install cantera")
    sys.exit(1)

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
YAML = os.path.join(ROOT, "data", "h2o2.yaml")

gas = ct.Solution(YAML)

print("=" * 60)
print("SPECIES")
print("=" * 60)
for sp in gas.species():
    print(f"  {sp.name:8s}  MW = {sp.molecular_weight:.6f} g/mol")

print()
print("=" * 60)
print("NASA7 COEFFICIENTS")
print("=" * 60)
for name in ["H2", "O2", "H2O"]:
    sp = gas.species(name)
    lo = sp.thermo.coeffs[1:8]
    hi = sp.thermo.coeffs[8:15]
    t_ranges = sp.thermo.min_temp, sp.thermo.coeffs[0], sp.thermo.max_temp
    print(f"  {name}: T = {t_ranges}")
    print(f"    low:  {lo}")
    print(f"    high: {hi}")

print()
print("=" * 60)
print(f"REACTIONS (total = {gas.n_reactions})")
print("=" * 60)
for i, rxn in enumerate(gas.reactions()):
    rtype = rxn.reaction_type  # Cantera 3.x: "Arrhenius", "three-body-Arrhenius", "falloff-Troe", etc.
    line = f"  [{i+1:2d}] {rxn.equation:<45s} {rtype}"
    if rxn.duplicate:
        line += "  [duplicate]"
    print(line)
    if rtype in ("Arrhenius", "three-body-Arrhenius"):
        r = rxn.rate
        print(f"        A={r.pre_exponential_factor:.4e}  b={r.temperature_exponent:.4f}  "
              f"Ea={r.activation_energy:.4f} J/mol")
        if rxn.third_body and rxn.third_body.efficiencies:
            print(f"        efficiencies: {rxn.third_body.efficiencies}")
    elif rtype.startswith("falloff-"):
        rate = rxn.rate
        rhi = rate.high_rate
        rlo = rate.low_rate
        print(f"        high: A={rhi.pre_exponential_factor:.4e}  b={rhi.temperature_exponent:.4f}  "
              f"Ea={rhi.activation_energy:.4f} J/mol")
        print(f"        low:  A={rlo.pre_exponential_factor:.4e}  b={rlo.temperature_exponent:.4f}  "
              f"Ea={rlo.activation_energy:.4f} J/mol")
        fc = rate.falloff_coeffs
        print(f"        Troe: A={fc[0]}  T3={fc[1]}  T1={fc[2]}  T2={fc[3] if len(fc)>3 else 'N/A'}")
        if rxn.third_body and rxn.third_body.efficiencies:
            print(f"        efficiencies: {rxn.third_body.efficiencies}")
    elif rtype == "pressure-dependent-Arrhenius":
        for p, r in rxn.rate.rates:
            print(f"        P={p:.4e} Pa  A={r.pre_exponential_factor:.4e}  "
                  f"b={r.temperature_exponent:.4f}  Ea={r.activation_energy:.4f} J/mol")

print()
print("=" * 60)
print("DUPLICATE REACTIONS")
print("=" * 60)
for i, rxn in enumerate(gas.reactions()):
    if rxn.duplicate:
        print(f"  [{i+1}] {rxn.equation}")


print()
print("=" * 60)
print("TRANSPORT — thermal conductivity [W/(m·K)]")
print("=" * 60)
gas_mix = ct.Solution(YAML, transport_model="mixture-averaged")
for name in ["H2", "O2", "H2O", "N2", "AR"]:
    for T in [300.0, 1000.0]:
        gas_mix.TPX = T, 101325.0, f"{name}:1"
        print(f"  lambda_{name}_{T:.0f}K = {gas_mix.thermal_conductivity:.12e}")

print()
print("=" * 60)
print("TRANSPORT — mixture thermal conductivity, air 300K/1000K")
print("=" * 60)
for T in [300.0, 1000.0]:
    gas_mix.TPX = T, 101325.0, "O2:0.21, N2:0.79"
    print(f"  lambda_air_{T:.0f}K = {gas_mix.thermal_conductivity:.12e}")
