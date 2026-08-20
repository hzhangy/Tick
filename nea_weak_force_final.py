#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_weak_force_final.py
修正 W/Z 质量计算：使用弱尺度跑动耦合 α(M_Z) 和实验混合角。
检验 NEA 公式链的最终精度。
"""
import numpy as np

# 基本常数
m_p = 0.938272  # GeV
alpha_low_inv = 25*np.sqrt(3)*np.pi + 1.0  # 低能精细结构常数倒数
alpha_low = 1.0 / alpha_low_inv

# v_h 用低能 α 计算（NEA 公式）
v_h = m_p * alpha_low_inv * 6.0 / np.pi  # GeV

# 弱尺度输入：α(M_Z) 和 sin^2θ_W(M_Z)
alpha_MZ = 1.0 / 127.9       # 实验弱尺度耦合
sin2_thetaW_MZ = 0.2312      # 实验弱混合角

# NEA 公式：m_W = v_h * sqrt(π α(M_Z) / sin^2θ_W(M_Z))
m_W = v_h * np.sqrt(np.pi * alpha_MZ / sin2_thetaW_MZ)
m_Z = m_W / np.sqrt(1.0 - sin2_thetaW_MZ)

# 实验值
m_W_exp = 80.379
m_Z_exp = 91.188

print("="*70)
print("  最终弱力 W/Z 质量计算（用弱尺度耦合）")
print("="*70)
print(f"  v_h = {v_h:.3f} GeV  (实验 246.22 GeV, 偏差 {abs(v_h-246.22)/246.22*100:.2f}%)")
print(f"  α(M_Z) = 1/127.9 = {alpha_MZ:.6f}")
print(f"  sin^2θ_W(M_Z) = {sin2_thetaW_MZ}")
print()
print(f"  m_W = {m_W:.3f} GeV   (实验 {m_W_exp} GeV, 偏差 {abs(m_W-m_W_exp)/m_W_exp*100:.2f}%)")
print(f"  m_Z = {m_Z:.3f} GeV   (实验 {m_Z_exp} GeV, 偏差 {abs(m_Z-m_Z_exp)/m_Z_exp*100:.2f}%)")
print(f"  m_W/m_Z = {m_W/m_Z:.4f}  (实验 {m_W_exp/m_Z_exp:.4f})")
print("="*70)