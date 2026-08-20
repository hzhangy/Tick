#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_strong_geometry.py
Strong force from exact K4 tetrahedron geometry.
No free parameters.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

a = 1.0                       # exchange range / lattice spacing
DELTA = 1.0 - np.sqrt(3)/2    # bandwidth deficit per locked node

# K4 vertices (alternate vertices of a cube of side a)
K4_VERTICES = np.array([
    [1, 1, 1],
    [1,-1,-1],
    [-1,1,-1],
    [-1,-1,1]
]) * (a/2)

def distance_matrix(d):
    """Return distances between vertices of K4 at origin and K4 shifted by (d,0,0)."""
    V2 = K4_VERTICES + np.array([d, 0, 0])
    dists = []
    for v1 in K4_VERTICES:
        for v2 in V2:
            dists.append(np.linalg.norm(v1 - v2))
    return np.array(dists)

def strong_potential(d):
    """Potential V(d) based on exact vertex distances."""
    dists = distance_matrix(d)
    # hard core: any vertex pair very close
    if np.any(dists < 0.5 * a):
        return np.inf
    # cross edges: vertex pairs within exchange range a
    k = np.sum(dists <= a)
    return -DELTA * k

# scan distances from 0.4 to 3.0 fm
d_vals = np.linspace(0.4, 3.0, 300)
V_pot = []
k_vals = []

for d in d_vals:
    Vd = strong_potential(d)
    if np.isinf(Vd):
        V_pot.append(np.nan)   # hard core region
        k_vals.append(0)
    else:
        V_pot.append(Vd)
        dists = distance_matrix(d)
        k_vals.append(np.sum(dists <= a))

V_pot = np.array(V_pot)
valid = ~np.isnan(V_pot)
d_valid = d_vals[valid]
V_valid = V_pot[valid]

print("="*70)
print("  Exact K4 tetrahedron geometry: strong potential")
print("="*70)
print(f"  Delta = {DELTA:.6f}")

if len(V_valid) > 0:
    d_min = d_valid[np.argmin(V_valid)]
    V_min = V_valid.min()
    print(f"  Minimum potential at d_eq = {d_min:.3f} fm, V_min = {V_min:.6f}")
    if (~valid).any():
        print(f"  Hard core region starts below d = {d_vals[~valid][-1]:.3f} fm")
    print(f"\n  {'d/fm':<10}{'cross edges k':<15}{'V (units Delta)':<15}")
    for d in [0.6, 0.8, d_min, 1.2, 1.5, 2.0]:
        Vd = strong_potential(d)
        if np.isinf(Vd):
            print(f"  {d:<10.3f}{'hard core':<15}{'inf':<15}")
        else:
            dists = distance_matrix(d)
            k = np.sum(dists <= a)
            print(f"  {d:<10.3f}{k:<15}{Vd:<15.6f}")

    # plot
    fig, axes = plt.subplots(1, 2, figsize=(13,5))
    ax = axes[0]
    ax.plot(d_vals, k_vals, 'b-', lw=2, label='cross edges k')
    ax.set_xlabel('d (fm)')
    ax.set_ylabel('k')
    ax.set_title('Cross-edge count from K4 geometry')
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(d_vals, V_pot, 'g-', lw=2, label='V(d)')
    ax.axvline(d_min, color='red', ls='--', alpha=0.5, label=f'd_eq={d_min:.3f} fm')
    ax.set_xlabel('d (fm)')
    ax.set_ylabel('V (units of Delta)')
    ax.set_title('Strong force potential (exact K4 geometry)')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('nea_strong_geometry.png', dpi=150)
    print("\nSaved: nea_strong_geometry.png")
else:
    print("All distances in hard core region; check geometry.")
print("="*70)