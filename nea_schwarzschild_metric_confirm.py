#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_schwarzschild_metric_confirm.py
Correct NEA strong-field metric:
  ds^2 = -f^2 dt^2 + dr^2/f^2 + r^2 dOmega^2
  with f^2 = 1 - 2A/r

Confirm it reproduces Schwarzschild light deflection.
"""
import numpy as np
from scipy.integrate import quad

A = 1.0  # geometric mass parameter

def f2(r, A):
    return 1.0 - 2.0*A/r

def integrand_schw(u, b, A):
    """du/dphi = sqrt(1/b^2 - u^2(1 - 2Au)) for Schwarzschild."""
    return 1.0/np.sqrt(np.maximum(1.0/b**2 - u**2*(1.0-2.0*A*u), 1e-20))

def find_u0(b, A):
    """Find u0 where du/dphi = 0."""
    func = lambda u: 1.0/b**2 - u**2*(1.0-2.0*A*u)
    # search root
    u_arr = np.linspace(1e-12, 0.99/A, 10000)
    vals = np.array([func(u) for u in u_arr])
    for i in range(len(u_arr)-1):
        if vals[i]*vals[i+1] <= 0:
            return (u_arr[i]*vals[i+1] - u_arr[i+1]*vals[i])/(vals[i+1]-vals[i])
    return None

def deflection_schw(b, A):
    """Light deflection angle for Schwarzschild."""
    u0 = find_u0(b, A)
    if u0 is None:
        return None
    integral, _ = quad(integrand_schw, 1e-12, u0, args=(b,A), limit=200)
    return 2.0*integral - np.pi

print("="*70)
print("  Correct NEA strong-field metric = Schwarzschild")
print("="*70)
print(f"  f_ext^2 = 1 - 2A/r, with A = {A}")
print(f"  {'b/A':<8}{'alpha (rad)':<16}{'alpha (deg)':<16}")
print("  " + "-"*50)

for br in [2.6, 2.7, 2.8, 3.0, 3.5, 4.0, 5.0, 10.0]:
    b = br*A
    alpha = deflection_schw(b, A)
    if alpha is None:
        print(f"  {br:<8}{'captured':<16}{'---':<16}")
    else:
        print(f"  {br:<8}{alpha:<16.6f}{np.degrees(alpha):<16.4f}")
print("="*70)