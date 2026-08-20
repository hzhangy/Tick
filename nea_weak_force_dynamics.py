#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_weak_force_dynamics.py
弱力动力学：混合角跑动、费米耦合、W宽度、中子寿命
"""
import numpy as np

# NEA 基础常数
alpha_low_inv = 25*np.sqrt(3)*np.pi + 1.0
alpha_MZ = 1/127.9
Delta = 1 - np.sqrt(3)/2
pi = np.pi

# 1) 混合角跑动
sin2_nea = 1/(1+pi) - Delta/(4*pi)
sin2_exp = 0.2312
print(f"sin^2θ_W(M_Z): 预言={sin2_nea:.5f}  实验={sin2_exp:.5f}  "
      f"偏差={abs(sin2_nea-sin2_exp)/sin2_exp*100:.2f}%")

# 2) 用修正后的混合角重算 W/Z
m_p = 0.938272  # GeV
v_h = m_p * alpha_low_inv * 6.0 / pi
m_W = v_h * np.sqrt(pi * alpha_MZ / sin2_nea)
m_Z = m_W / np.sqrt(1 - sin2_nea)
print(f"\nv_h = {v_h:.3f} GeV  (实验 246.22)")
print(f"m_W = {m_W:.3f} GeV  (实验 80.379)  偏差={abs(m_W-80.379)/80.379*100:.2f}%")
print(f"m_Z = {m_Z:.3f} GeV  (实验 91.188)  偏差={abs(m_Z-91.188)/91.188*100:.2f}%")

# 3) 费米耦合常数
G_F_nea = pi * alpha_MZ / (np.sqrt(2) * sin2_nea * m_W**2)
G_F_exp = 1.166e-5
print(f"\nG_F = {G_F_nea:.3e} GeV^-2  (实验 {G_F_exp:.3e})  偏差={abs(G_F_nea-G_F_exp)/G_F_exp*100:.2f}%")

# 4) W 玻色子宽度
# Standard W width: Gamma = G_F * m_W^3 / (6*sqrt(2)*pi) * N_channels
# N_channels = 3 lepton + 2 quark flavors × 3 colors × (1 + alpha_s/pi)
alpha_s_MZ = 0.118
N_ch = 3 + 2*3*(1 + alpha_s_MZ/np.pi)   # ≈ 9.22
Gamma_W_nea = G_F_nea * m_W**3 / (6*np.sqrt(2)*np.pi) * N_ch
Gamma_W_exp = 2.09
print(f"\nΓ_W = {Gamma_W_nea:.2f} GeV  (实验 {Gamma_W_exp})  偏差={abs(Gamma_W_nea-Gamma_W_exp)/Gamma_W_exp*100:.1f}%")
tau_W = 6.58e-25 / Gamma_W_nea
print(f"τ_W = {tau_W:.2e} s  (实验约 3e-25 s)")

# 5) 中子 beta 衰变
Vud = 0.9736
dm = 1.293e-3  # GeV
Gamma_n = G_F_nea**2 * Vud**2 * dm**5 / (60*pi**3)
tau_n = 6.58e-25 / Gamma_n
print(f"\nΓ_n = {Gamma_n:.2e} GeV")
print(f"τ_n = {tau_n:.0f} s  (实验 880 s)  量级={'正确' if 300<tau_n<3000 else '需修正'}")