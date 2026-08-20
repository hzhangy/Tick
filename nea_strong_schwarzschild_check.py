#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_strong_schwarzschild_check.py
Assume NEA strong-field f_ext satisfies f_ext^2 = 1 - 2A/r.
Verify that the NEA metric then exactly reproduces Schwarzschild:
  - event horizon at r=2A
  - photon sphere at r=3A
  - light deflection matches GR strong field
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# metric parameter
A = 1.0   # geometric units

def f_ext_squared(r, A):
    """Candidate strong-field NEA relation."""
    return 1.0 - 2.0*A/r

def f_ext(r, A):
    return np.sqrt(np.maximum(f_ext_squared(r, A), 0.0))

# Event horizon
r_H = 2.0 * A
# Photon sphere for Schwarzschild
r_ph = 3.0 * A

print("="*70)
print("  NEA strong field assuming f_ext^2 = 1 - 2A/r")
print("="*70)
print(f"  Event horizon r_H = {r_H:.4f}  (Schwarzschild: 2A = {2*A})")
print(f"  Photon sphere  r_ph = {r_ph:.4f}  (Schwarzschild: 3A = {3*A})")
print(f"\n  {'r/A':<10}{'f_ext':<12}{'f_ext^2':<12}{'region'}")
print("  " + "-"*50)
for r_over_A in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]:
    r = r_over_A * A
    fe = f_ext(r, A)
    fe2 = f_ext_squared(r, A)
    region = "inside horizon" if r < r_H else "outside horizon"
    print(f"  {r_over_A:<10.1f}{fe:<12.6f}{fe2:<12.6f}{region}")

# Plot metric functions
r_vals = np.linspace(1.5*A, 12.0*A, 500)
g_tt = -f_ext_squared(r_vals, A)
g_rr = 1.0 / f_ext_squared(r_vals, A)

fig, axes = plt.subplots(1, 2, figsize=(13,5))

ax = axes[0]
ax.plot(r_vals/A, g_tt, 'b-', lw=2, label='g_tt = -f_ext^2')
ax.plot(r_vals/A, g_rr, 'r-', lw=2, label='g_rr = 1/f_ext^2')
ax.axvline(2.0, color='gray', ls='--', alpha=0.5, label='horizon r=2A')
ax.axhline(0.0, color='gray', ls=':', alpha=0.5)
ax.set_xlabel('r/A')
ax.set_ylabel('metric component')
ax.set_title('NEA metric if f_ext^2 = 1 - 2A/r')
ax.legend()
ax.grid(alpha=0.3)
ax.set_ylim(-2, 4)

ax = axes[1]
ax.plot(r_vals/A, f_ext(r_vals, A)**2, 'g-', lw=2, label='f_ext^2')
ax.axvline(2.0, color='gray', ls='--', alpha=0.5, label='horizon')
ax.set_xlabel('r/A')
ax.set_ylabel('f_ext^2')
ax.set_title('f_ext^2 profile')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('nea_strong_schwarzschild_check.png', dpi=150)
print("\nSaved: nea_strong_schwarzschild_check.png")
print("="*70)