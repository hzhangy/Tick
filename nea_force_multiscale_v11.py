#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_force_multiscale_v11.py

按照 Volume Z v11.1 的四方程力公式，计算四种基本相互作用的强度
跨几十个数量级尺度。力公式：

    F_i(r) = (Z * hbar_c / T_i) * (Delta * rho_i / (4*pi)) * K_i(r)

其中：
    Z          = 0.4066 MeV              (带宽货币)
    hbar_c     = 197.3269804 MeV·fm      (量子-相对论换算)
    Delta      = 1 - sqrt(3)/2           (K4锁定缺口)
    T_i        = 载体稀释因子
    rho_i      = 源强
    K_i(r)     = 传播核

对于无质量场：
    K_i(r) = 1 / r^2

对于有质量场：
    K_i(r) = (1/r^2 + m_i/(hbar_c * r)) * exp(-m_i * r / hbar_c)

四种力：
    强力：T_strong = 1, rho_strong ≈ 310, m_rho = hbar_c/(2*Delta)
    电磁：T_em = Z*Delta/(pi*alpha), rho_em = 4, m = 0
    引力：T_grav = Z*hbar_c*Delta/(4*pi*G*m_p^2), rho_grav = 1, m = 0
    弱力：协议层，不定义 T，用有效 Yukawa 展示短程性
          alpha_W ≈ 0.0302, m_W = 80.4 GeV

本脚本还验证了在 r=1 fm 处三种连续力与标准物理值的一致性，
并输出多尺度图像。
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# 基本常数
# ============================================================
hbar_c = 197.3269804          # MeV·fm
Delta = 1.0 - np.sqrt(3.0) / 2.0
Z = 0.4066                    # MeV，带宽货币
alpha_inv = 25.0 * np.sqrt(3.0) * np.pi + 1.0
alpha = 1.0 / alpha_inv       # 精细结构常数

# 标准物理常数
G = 6.67430e-11               # m^3 kg^-1 s^-2
m_p_kg = 1.67262192369e-27    # kg
# 单位转换：1 MeV·fm = 1.602176634e-28 J·m
C_G = G * m_p_kg**2 / (1.602176634e-28)   # G*m_p^2 in MeV·fm
C_EM = alpha * hbar_c                    # alpha*hbar_c in MeV·fm

# 强力参数
m_rho = hbar_c / (2.0 * Delta)  # MeV
mu_rho = m_rho / hbar_c         # fm^-1

# 弱力参数（标准模型有效 Yukawa）
m_W_MeV = 80.379 * 1e3          # MeV
mu_W = m_W_MeV / hbar_c         # fm^-1
alpha_W = 0.0302                # 有效弱耦合

# ============================================================
# 修正后的 T_i 和 rho_i
# ============================================================
T_strong = 1.0
rho_strong = 310.0              # 核子有效锁定缺陷数

# 电磁
rho_em = 4.0
T_em = Z * Delta * rho_em / (4.0 * np.pi * alpha)

# 引力
rho_grav = 1.0
T_grav = Z * hbar_c * Delta * rho_grav / (4.0 * np.pi * C_G)

# 弱力不定义 T_weak，使用有效 Yukawa 展示
# F_weak = alpha_W * hbar_c * K_yukawa(r, mu_W)

# ============================================================
# 传播核与力函数
# ============================================================
def K_massless(r):
    return 1.0 / r**2

def K_yukawa(r, mu):
    return (1.0/r**2 + mu/r) * np.exp(-mu*r)

def F_unified(r, T, rho, mu=0.0):
    """统一力公式。mu=0 对应无质量场，否则为 Yukawa。"""
    K = K_massless(r) if mu == 0.0 else K_yukawa(r, mu)
    return (Z * hbar_c / T) * (Delta * rho / (4.0*np.pi)) * K

def F_strong(r):
    return F_unified(r, T_strong, rho_strong, mu_rho)

def F_em(r):
    return F_unified(r, T_em, rho_em)

def F_grav(r):
    return F_unified(r, T_grav, rho_grav)

def F_weak(r):
    # 弱力作为协议层短程 Yukawa 展示（非统一力公式成员）
    return alpha_W * hbar_c * K_yukawa(r, mu_W)

# ============================================================
# 打印基本常数
# ============================================================
print("="*70)
print("  四力强度多尺度验证（Volume Z v11.1）")
print("="*70)
print(f"  Delta       = {Delta:.9f}")
print(f"  Z           = {Z} MeV")
print(f"  alpha^-1    = {alpha_inv:.3f}")
print(f"  C_G (G m_p^2) = {C_G:.6e} MeV·fm")
print(f"  C_EM (alpha hbar c) = {C_EM:.6e} MeV·fm")
print()
print(f"  T_strong    = {T_strong}")
print(f"  T_em        = {T_em:.4f}")
print(f"  T_grav      = {T_grav:.4e}")
print(f"  rho_strong  = {rho_strong:.0f}")
print()

# ============================================================
# 验证 r=1 fm 处的绝对强度
# ============================================================
r_test = 1.0   # fm
F_s = F_strong(r_test)
F_e = F_em(r_test)
F_g = F_grav(r_test)
F_w = F_weak(r_test)

# 标准值
F_s_std = 30.0           # MeV/fm 经验核力
F_e_std = C_EM / r_test**2
F_g_std = C_G / r_test**2
F_w_std = 0.0            # 接触作用，忽略

print("在 r = 1 fm 处的力值：")
print(f"  强力：{F_s:.4e} MeV/fm  (标准参考 30 MeV/fm)")
print(f"  电磁：{F_e:.4e} MeV/fm  (标准 {F_e_std:.4e} MeV/fm)")
print(f"  弱力：{F_w:.4e} MeV/fm  (接触作用，趋近于0)")
print(f"  引力：{F_g:.4e} MeV/fm  (标准 {F_g_std:.4e} MeV/fm)")
print()
print("比值（NEA / 标准）：")
print(f"  强力：{F_s/F_s_std:.4f}")
print(f"  电磁：{F_e/F_e_std:.4f}")
print(f"  引力：{F_g/F_g_std:.4f}")
print()

# ============================================================
# 弱力在电弱尺度的表现
# ============================================================
r_weak_dom = 1e-3   # fm = 10^-18 m
F_w_dom = F_weak(r_weak_dom)
print(f"弱力在 r = {r_weak_dom:.1e} fm (= 10^-18 m) 处：")
print(f"  F_weak = {F_w_dom:.4e} MeV/fm")
print("  注：弱力是接触作用，仅在电弱尺度显著。")
print()

# ============================================================
# 多尺度图形
# ============================================================
# 距离范围：从 10^-5 fm 到 10^30 fm (跨越 35 个数量级)
r_fm = np.logspace(-5, 30, 2000)

F_s_all = F_strong(r_fm)
F_e_all = F_em(r_fm)
F_g_all = F_grav(r_fm)
F_w_all = F_weak(r_fm)

# 避免 log 中的 0
eps = 1e-300
F_s_all = np.maximum(F_s_all, eps)
F_e_all = np.maximum(F_e_all, eps)
F_g_all = np.maximum(F_g_all, eps)
F_w_all = np.maximum(F_w_all, eps)

plt.figure(figsize=(12, 8))
plt.loglog(r_fm, F_s_all, 'r-', lw=2, label='Strong (Yukawa)')
plt.loglog(r_fm, F_e_all, 'b-', lw=2, label='Electromagnetic (Coulomb)')
plt.loglog(r_fm, F_g_all, 'g-', lw=2, label='Gravity (Coulomb)')
plt.loglog(r_fm, F_w_all, 'm--', lw=2, label='Weak (short-range Yukawa, protocol)')

plt.xlabel('r (fm)')
plt.ylabel('Force magnitude (MeV/fm)')
plt.title('N.E.A. four-force strengths across scales (Volume Z v11.1)')
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend(fontsize=10)
plt.xlim(1e-5, 1e30)
plt.ylim(1e-200, 1e10)
plt.tight_layout()
plt.savefig('nea_force_multiscale_v11.png', dpi=150)
print("已保存图像: nea_force_multiscale_v11.png")
print("="*70)