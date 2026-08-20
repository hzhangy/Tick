#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_em_charge_propagation.py
Electromagnetic force as K4 unpaired direction phase on C8 network.

Build a 3D C8 grid, place a K4 defect with one unpaired phase,
solve discrete phase propagation, check Coulomb 1/r emerges.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ================ 3D C8 grid ================
L = 40
N = L**3

def idx(x, y, z):
    return (x * L + y) * L + z

def build_adjacency(L):
    """Build 6-neighbor C8 grid."""
    adj = [[] for _ in range(L**3)]
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = idx(x, y, z)
                for dx, dy, dz in [(1,0,0),(0,1,0),(0,0,1)]:
                    nx, ny, nz = x+dx, y+dy, z+dz
                    if nx < L and ny < L and nz < L:
                        j = idx(nx, ny, nz)
                        adj[i].append(j)
                        adj[j].append(i)
    return adj

def solve_phase_field(adj, source_idx, max_iter=5000):
    """
    Solve discrete Poisson equation on C8 grid.
    Source = unpaired direction phase at K4 defect.
    """
    n = len(adj)
    phi = np.zeros(n)
    rho = np.zeros(n)
    rho[source_idx] = 1.0   # unit charge from unpaired phase

    for _ in range(max_iter):
        phi_new = phi.copy()
        for i in range(n):
            # Jacobi update: phi_i = (sum_j phi_j + rho_i) / degree
            deg = len(adj[i])
            if deg > 0:
                neighbor_sum = sum(phi[j] for j in adj[i])
                phi_new[i] = (neighbor_sum + rho[i]) / deg
        # boundaries fixed at zero (implicit via not updating boundaries)
        # interior only
        for x in range(1, L-1):
            for y in range(1, L-1):
                for z in range(1, L-1):
                    i = idx(x, y, z)
                    deg = len(adj[i])
                    neighbor_sum = sum(phi[j] for j in adj[i])
                    phi_new[i] = (neighbor_sum + rho[i]) / deg
        phi = phi_new
    return phi

print("="*70)
print("  EM: K4 unpaired phase on C8 -> Coulomb potential")
print("="*70)

adj = build_adjacency(L)
source = idx(L//2, L//2, L//2)
phi = solve_phase_field(adj, source, max_iter=3000)

# Check 1/r decay along x-axis
print(f"  Source at center, charge=1")
print(f"\n  {'r':<8}{'phi(r)':<14}{'phi*r':<14}")
print("  " + "-"*40)
r_vals = []
phi_r_vals = []
for r in range(1, 12):
    i = idx(L//2 + r, L//2, L//2)
    val = phi[i]
    r_vals.append(r)
    phi_r_vals.append(val)
    print(f"  {r:<8}{val:<14.6f}{val*r:<14.6f}")

# Fit phi ~ A/r + C
A_coeff = np.mean(np.array(phi_r_vals) * np.array(r_vals))
print(f"\n  Average phi*r = {A_coeff:.6f}")
print(f"  (Constant phi*r indicates clean 1/r Coulomb potential)")
print("="*70)