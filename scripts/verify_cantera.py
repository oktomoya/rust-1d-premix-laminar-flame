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
    rtype = type(rxn).__name__
    line = f"  [{i+1:2d}] {rxn.equation:<45s} {rtype}"
    if hasattr(rxn, 'duplicate') and rxn.duplicate:
        line += "  [duplicate]"
    print(line)
    if rtype in ("ElementaryReaction", "ThreeBodyReaction"):
        r = rxn.rate
        print(f"        A={r.pre_exponential_factor:.4e}  b={r.temperature_exponent:.4f}  "
              f"Ea={r.activation_energy:.4f} J/mol")
    elif rtype == "FalloffReaction":
        rhi = rxn.high_rate
        rlo = rxn.low_rate
        print(f"        high: A={rhi.pre_exponential_factor:.4e}  b={rhi.temperature_exponent:.4f}  "
              f"Ea={rhi.activation_energy:.4f} J/mol")
        print(f"        low:  A={rlo.pre_exponential_factor:.4e}  b={rlo.temperature_exponent:.4f}  "
              f"Ea={rlo.activation_energy:.4f} J/mol")
        fc = rxn.falloff.parameters
        print(f"        Troe: A={fc[0]}  T3={fc[1]}  T1={fc[2]}  T2={fc[3] if len(fc)>3 else 'N/A'}")
        if rxn.efficiencies:
            print(f"        efficiencies: {rxn.efficiencies}")
    elif rtype == "PlogReaction":
        for p, r in rxn.rates:
            print(f"        P={p:.4e} Pa  A={r.pre_exponential_factor:.4e}  "
                  f"b={r.temperature_exponent:.4f}  Ea={r.activation_energy:.4f} J/mol")

print()
print("=" * 60)
print("DUPLICATE REACTIONS")
print("=" * 60)
for i, rxn in enumerate(gas.reactions()):
    if rxn.duplicate:
        print(f"  [{i+1}] {rxn.equation}")
