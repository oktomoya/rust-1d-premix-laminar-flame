"""
Solve H2/air stoichiometric freely-propagating flame with Cantera,
then export a remapped profile as initial condition for the Rust solver.

The Cantera solution (50 mm domain) is clipped to a window centered on
the flame and remapped to the Rust domain [0, domain_length_mm].

Output CSV columns: z_m, T_K, M_kg_m2s, Y_H2, Y_H, Y_O, Y_O2, Y_OH, Y_H2O, Y_HO2, Y_H2O2, Y_AR, Y_N2

Usage:
  .venv/bin/python3 scripts/cantera_h2air_solution.py
"""

import sys, os
import numpy as np
from scipy.interpolate import interp1d
import cantera as ct

# ---------------------------------------------------------------------------
# Flame configuration  (must match hydrogen_air.toml)
# ---------------------------------------------------------------------------
P_PA        = 101325.0
T_IN        = 300.0
PHI         = 1.0
RUST_DOMAIN = 0.020   # m  — must match [flame] domain_length in TOML
FLAME_POS   = 0.40    # flame center at 40% of domain (matches solver_driver.rs)
CT_WIDTH    = 0.050   # m  — wider Cantera domain so both edges are equilibrated

SPECIES_ORDER = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'AR', 'N2']

mech_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'h2o2.yaml')
gas = ct.Solution(mech_file)
gas.set_equivalence_ratio(PHI, fuel={'H2': 1.0}, oxidizer={'O2': 0.21, 'N2': 0.79})
gas.TP = T_IN, P_PA

print(f"Unburned mixture: {gas.mole_fraction_dict()}")
print(f"T_in = {T_IN} K,  P = {P_PA} Pa,  phi = {PHI}")

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
f = ct.FreeFlame(gas, width=CT_WIDTH)
f.set_refine_criteria(ratio=3, slope=0.05, curve=0.10)
f.transport_model = 'mixture-averaged'
print("\nSolving with Cantera...")
f.solve(loglevel=1, auto=True)

Su   = f.velocity[0]
Tmax = f.T.max()
M_ct = (f.density * f.velocity).mean()
print(f"\nFlame speed Su = {Su:.4f} m/s")
print(f"T_max = {Tmax:.1f} K")
print(f"Mass flux M = {M_ct:.6f} kg/(m²·s)")

# ---------------------------------------------------------------------------
# Find flame center in the Cantera domain
# ---------------------------------------------------------------------------
T_arr = f.T
T_mid = T_arr[0] + 0.5 * (T_arr[-1] - T_arr[0])
i_mid = np.argmin(np.abs(T_arr - T_mid))
z_flame_ct = f.grid[i_mid]
print(f"Flame center (T={T_mid:.0f} K) at z = {z_flame_ct*1000:.2f} mm in Cantera domain")

# ---------------------------------------------------------------------------
# Define a window in Cantera coordinates that maps onto the Rust domain
#   z_rust = z_flame_pos * RUST_DOMAIN  (flame center position in Rust domain)
#   z_win_left  = z_flame_ct - z_flame_pos * RUST_DOMAIN
#   z_win_right = z_flame_ct + (1 - z_flame_pos) * RUST_DOMAIN
# ---------------------------------------------------------------------------
z_flame_rust = FLAME_POS * RUST_DOMAIN
z_win_left   = z_flame_ct - z_flame_rust
z_win_right  = z_win_left + RUST_DOMAIN
print(f"Cantera window: [{z_win_left*1000:.2f}, {z_win_right*1000:.2f}] mm → mapped to [0, {RUST_DOMAIN*1000:.0f}] mm")

# ---------------------------------------------------------------------------
# Build interpolators on the Cantera grid
# ---------------------------------------------------------------------------
z_ct = f.grid
interp_T   = interp1d(z_ct, f.T,        kind='linear', fill_value='extrapolate')
interp_M   = interp1d(z_ct, f.density * f.velocity, kind='linear', fill_value='extrapolate')

Y_ct = f.Y  # shape: (n_species, n_points)
species_ct = [gas.species_name(k) for k in range(gas.n_species)]
interp_Y = {}
for s in SPECIES_ORDER:
    k_ct = gas.species_index(s)
    interp_Y[s] = interp1d(z_ct, Y_ct[k_ct, :], kind='linear', fill_value='extrapolate')

# ---------------------------------------------------------------------------
# Sample onto the Rust grid (100 uniform points by default, but we just
# output the full Cantera grid remapped — the Rust code interpolates onto
# whatever grid it's using)
# Use 200 points for a smooth initial condition.
# ---------------------------------------------------------------------------
N_OUT = 200
z_rust = np.linspace(0.0, RUST_DOMAIN, N_OUT)
z_cantera_mapped = z_win_left + z_rust  # linear mapping

T_out = interp_T(z_cantera_mapped)
M_out = interp_M(z_cantera_mapped)
Y_out = {s: np.clip(interp_Y[s](z_cantera_mapped), 0.0, 1.0) for s in SPECIES_ORDER}

# Renormalize Y to sum=1 at each point
Y_arr = np.array([Y_out[s] for s in SPECIES_ORDER])  # (10, N_OUT)
Y_sum = Y_arr.sum(axis=0)
Y_arr /= np.where(Y_sum > 0, Y_sum, 1.0)

# Clamp temperatures to physical range
T_out = np.clip(T_out, T_IN, Tmax * 1.05)

# ---------------------------------------------------------------------------
# Write output CSV
# ---------------------------------------------------------------------------
out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'cantera_h2air_initial.csv')
with open(out_path, 'w') as fh:
    y_cols = ','.join(f'Y_{s}' for s in SPECIES_ORDER)
    fh.write(f'z_m,T_K,M_kg_m2s,{y_cols}\n')
    for j in range(N_OUT):
        y_vals = ','.join(f'{Y_arr[k,j]:.10e}' for k in range(len(SPECIES_ORDER)))
        fh.write(f'{z_rust[j]:.10e},{T_out[j]:.6f},{M_out[j]:.10e},{y_vals}\n')

print(f"\nSaved to: {out_path}")
print(f"  {N_OUT} grid points, z = [0, {RUST_DOMAIN*1000:.0f}] mm")
print(f"  T range: [{T_out[0]:.1f}, {T_out[-1]:.1f}] K")

# ---------------------------------------------------------------------------
# Print TOML values
# ---------------------------------------------------------------------------
print(f"\n=== Recommended hydrogen_air.toml values ===")
print(f"su_initial_guess = {Su:.4f}   # m/s  (Cantera mixture-averaged)")
print(f"# T_ad = {Tmax:.1f} K,  M = {M_ct:.6f} kg/(m²·s)")
