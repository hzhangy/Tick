#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_em_alpha_zk.py
Check if electromagnetic coupling can be expressed via NEA constants:
alpha^{-1} = 25*sqrt(3)*pi + 1, and ZY energy unit.
"""
import numpy as np

# NEA constants
alpha_inv_nea = 25.0 * np.sqrt(3) * np.pi + 1
alpha_nea = 1.0 / alpha_inv_nea
ZY_MeV = 0.4066
hbar_c_MeV_fm = 197.32698

# Standard electromagnetic fine structure constant
alpha_exp = 1.0 / 137.035999

print("="*70)
print("  Electromagnetic coupling from NEA alpha and ZY")
print("="*70)
print(f"  alpha^{-1} NEA = {alpha_inv_nea:.6f}")
print(f"  alpha^{-1} exp = {137.035999:.6f}")
print(f"  alpha NEA = {alpha_nea:.8f}")
print(f"  alpha exp = {alpha_exp:.8f}")
print(f"  relative deviation = {(alpha_nea-alpha_exp)/alpha_exp*100:.4f}%")

# Standard Coulomb force constant: F = e^2 / r^2
# e^2 = alpha * hbar * c
e2_exp = alpha_exp * hbar_c_MeV_fm
e2_nea = alpha_nea * hbar_c_MeV_fm

print(f"\n  Standard e^2 (MeV*fm) = {e2_exp:.6f}")
print(f"  NEA e^2 (MeV*fm)      = {e2_nea:.6f}")
print(f"  ratio NEA/exp         = {e2_nea/e2_exp:.6f}")

# Check if kappa_e can be alpha * ZY or e^2 / something
kappa_candidate_1 = alpha_nea * ZY_MeV
kappa_candidate_2 = e2_nea   # already MeV*fm
print(f"\n  Candidate kappa_e = alpha*ZY = {kappa_candidate_1:.6f} MeV")
print(f"  Candidate kappa_e = e^2 = {e2_nea:.6f} MeV*fm")
print(f"  These are the only natural NEA electromagnetic coupling candidates.")
print("="*70)