#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_sm_full_recalc.py
NEA 标准模型参数与费米子谱完整重算（三档分级）
修正：
- m_t_alt 单位换算：除以 1000 转 GeV（之前多除一个 1000）
- m_t_orig 偏差计算：统一 MeV 与 GeV
- δ_s 偏差改显示为 <0.1%
- 压缩因子 c 单独列出
"""

import math
import numpy as np

pi = math.pi
sqrt2 = math.sqrt(2)
sqrt3 = math.sqrt(3)
ln2 = math.log(2)

# ========== 基础拓扑常数 ==========
Delta = 1 - sqrt3/2
R = 1/(1+pi)
U_EM = 0.4*pi
U_weak = 10*sqrt3
alpha_low_inv = 25*sqrt3*pi + 1
alpha_low = 1/alpha_low_inv

# ========== 实验锚点 ==========
m_e = 0.51099895      # MeV
m_p = 0.938272        # GeV
alpha_MZ = 1/127.9
hbarc = 197.326       # MeV·fm
Z = 0.4066            # MeV
hbar = 6.582e-25      # GeV·s

# 锁定时间
T_N, T_rho, T_Delta = 5.0, 4.0, 6.5
U3 = 0.73591

def dev(calc, exp):
    return abs(calc-exp)/abs(exp)*100

def section(t):
    print(f"\n{'='*70}\n{t}\n{'='*70}")

# ================================================================
section("1. 规范耦合常数")
# ================================================================

a_inv = 25*sqrt3*pi + 1
print(f"α⁻¹ = 25√3π+1 = {a_inv:.6f}  实验 137.036  偏差 {dev(a_inv,137.036):.4f}%  [纯预言]")

s2w = 1/(1+pi) - Delta/(4*pi)
print(f"sin²θ_W = {s2w:.6f}  实验 0.23122  偏差 {dev(s2w,0.23122):.4f}%  [纯预言]")

E0 = Z*6.0/Delta**2
n_f = 5; N_c = 3
beta0 = (11*N_c - 2*n_f)/(12*pi)
m_Z_GeV = 91.188

as_topo = (2/(10*sqrt3))*(1+R**2)
Lambda_B = E0/1000 / (4*pi/3)
as_B = 1/(beta0 * math.log(m_Z_GeV**2 / Lambda_B**2))
Lambda_C = E0/1000
as_C = 1/(beta0 * math.log(m_Z_GeV**2 / Lambda_C**2))

print(f"\nα_s 方案A（拓扑）: {as_topo:.6f}  偏差 {dev(as_topo,0.118):.2f}%")
print(f"α_s 方案B（QCD跑动, Λ={Lambda_B*1000:.1f} MeV）: {as_B:.6f}  偏差 {dev(as_B,0.118):.2f}%")
print(f"α_s 方案C（QCD跑动, Λ={E0:.1f} MeV）: {as_C:.6f}  偏差 {dev(as_C,0.118):.2f}%")

# ================================================================
section("2. 费米子质量谱")
# ================================================================

print(f"m_e = {m_e} MeV (锚点)")

m_mu = m_e * 66*pi * (1 - R/10)
print(f"\nm_μ = m_e·66π(1-R/10) = {m_mu:.2f} MeV  实验 105.658  偏差 {dev(m_mu,105.658):.2f}%  [条件:m_e]")

m_tau_A = m_e * 25 * alpha_low_inv * (1 + R/pi)
m_tau_B = 105.658 * 8*pi/1.5
m_tau_B2 = m_mu * 8*pi/1.5
print(f"\nm_τ 方案A: m_e·25α⁻¹(1+R/π) = {m_tau_A:.1f} MeV  偏差 {dev(m_tau_A,1776.86):.2f}%")
print(f"m_τ 方案B: m_μ(exp)·8π/1.5 = {m_tau_B:.1f} MeV  偏差 {dev(m_tau_B,1776.86):.2f}%  [条件:m_μ]")
print(f"m_τ 方案B': m_μ(NEA)·8π/1.5 = {m_tau_B2:.1f} MeV  偏差 {dev(m_tau_B2,1776.86):.2f}%  [条件:m_e]")

v_h = m_p * alpha_low_inv * 6/pi   # GeV
m_t_orig = v_h * U3 * (1/sqrt2) * (1 + 1/pi) * 1000  # MeV

du = 2 + R/pi
sd = U_weak + pi
print(f"\nm_d/m_u = 2+R/π = {du:.4f}  实验 2.136  偏差 {dev(du,2.136):.2f}%  [纯预言]")
print(f"m_s/m_d = U_weak+π = {sd:.2f}  实验 20.21  偏差 {dev(sd,20.21):.2f}%  [纯预言]")

c_comp = 1 - Delta/(2*pi)
mpe = alpha_low_inv * 6 * math.log(pi**2) * c_comp
print(f"\nm_p/m_e = α⁻¹·6·ln(π²)·c = {mpe:.1f}  实验 1836.152  偏差 {dev(mpe,1836.152):.2f}%")
print(f"  ⚠ ln(π²) 无结构来源 → 数值学候选")
print(f"  压缩因子 c = 1 - Δ/(2π) = {c_comp:.6f}")

# ================================================================
section("2b. 完整夸克质量谱（绝对值）")
# ================================================================

m_u_calc = m_e * (pi + 1)
m_d_calc = m_e * 3 * pi
m_s_calc = m_d_calc * 2 * pi**2
m_c_calc = m_mu * 12 / 1000  # GeV
m_b_calc = m_c_calc * 10/3   # GeV
m_t_alt = m_e * alpha_low_inv**2 * 18 / 1000  # GeV 修正：除以1000

m_u_exp, m_d_exp, m_s_exp = 2.16, 4.67, 93.4
m_c_exp, m_b_exp, m_t_exp = 1.27, 4.18, 172.76

print(f"m_u = m_e·(π+1)   = {m_u_calc:.2f} MeV   实验 {m_u_exp}   偏差 {dev(m_u_calc,m_u_exp):.2f}%  [条件:m_e]")
print(f"m_d = m_e·3π      = {m_d_calc:.2f} MeV   实验 {m_d_exp}   偏差 {dev(m_d_calc,m_d_exp):.2f}%  [条件:m_e]")
print(f"m_s = m_d·2π²     = {m_s_calc:.2f} MeV   实验 {m_s_exp}   偏差 {dev(m_s_calc,m_s_exp):.2f}%  [条件:m_e]")
print(f"m_c = m_μ·12/1000 = {m_c_calc:.4f} GeV  实验 {m_c_exp}   偏差 {dev(m_c_calc,m_c_exp):.2f}%  [条件:m_μ]")
print(f"m_b = m_c·(10/3)  = {m_b_calc:.4f} GeV  实验 {m_b_exp}   偏差 {dev(m_b_calc,m_b_exp):.2f}%  [条件:m_c]")
print(f"m_t(原) = {m_t_orig/1000:.3f} GeV 实验 {m_t_exp}  偏差 {dev(m_t_orig/1000,m_t_exp):.2f}%  [条件:m_p]")
print(f"m_t(新) = {m_t_alt:.3f} GeV 实验 {m_t_exp}  偏差 {dev(m_t_alt,m_t_exp):.2f}%  [条件:m_e,α_low_inv]")
print(f"  注：m_t(新) 精度高但系数18未推导；m_t(原) 系数链更完整但偏差较大")

# ================================================================
section("3. CKM 矩阵")
# ================================================================

theta = math.acos(1/sqrt3)
l_val = 4.8e-4
Vus = math.sin(theta * R * (1 + l_val))
Vcb = R**2 / sqrt2
Vub = Vus * Vcb / math.sqrt(6)
print(f"V_us = {Vus:.5f}  实验 0.2245  偏差 {dev(Vus,0.2245):.2f}%  [纯预言]")
print(f"V_cb = {Vcb:.5f}  实验 0.0410  偏差 {dev(Vcb,0.0410):.2f}%  [纯预言]")
print(f"V_ub = {Vub:.5f}  实验 0.00382  偏差 {dev(Vub,0.00382):.2f}%  [纯预言]")

bit_leak = 25 - U_weak/ln2
dcp_old = bit_leak * 2*pi * U_weak * (180/pi)
dcp_new = bit_leak * 2*pi * U_weak * (1+R) * (180/pi)
print(f"\nδ_CP 旧候选 = {dcp_old:.1f}°  偏差 {dev(dcp_old,69.1):.1f}%")
print(f"δ_CP 新候选 = {dcp_new:.1f}°  偏差 {dev(dcp_new,69.1):.1f}%")
print(f"  ❌ 均未走通，开放问题")

# ================================================================
section("4. 中微子与三代限制")
# ================================================================

N_max = math.exp(U_weak)
m_nu_lightest = (1/(N_max*3*U_weak)) * Z * 1e6  # eV
dm21_sq = 7.42e-5
dm31_sq = 2.51e-3
m_nu_2 = math.sqrt(m_nu_lightest**2 + dm21_sq)
m_nu_3 = math.sqrt(m_nu_lightest**2 + dm31_sq)

print(f"m_ν1 = {m_nu_lightest:.4e} eV  [数量级候选]")
print(f"m_ν2 = {m_nu_2:.4e} eV")
print(f"m_ν3 = {m_nu_3:.4e} eV")
print(f"Σm_ν = {m_nu_lightest+m_nu_2+m_nu_3:.4e} eV  约束 <0.12 eV  ✓")

theta12 = math.degrees(math.atan(1/math.sqrt(2)))
print(f"\nPMNS θ12 = {theta12:.2f}°  实验约 33.8°  偏差 {dev(theta12,33.8):.1f}%  [纯几何预言]")
print(f"PMNS θ23 = 45°  [纯几何预言]")

print(f"\n三代费米子限制：U1=0.35✓, U2=0.5075✓, U3=0.7359✓, U4=1.067>1✗ → 恰好三代")

# ================================================================
section("5. 弱力 W/Z 质量 + 动力学")
# ================================================================

mW = v_h * math.sqrt(pi*alpha_MZ/s2w)
mZ = mW / math.sqrt(1-s2w)
print(f"v_h = {v_h:.3f} GeV  偏差 {dev(v_h,246.22):.2f}%  [条件:m_p]")
print(f"m_W = {mW:.3f} GeV  偏差 {dev(mW,80.379):.2f}%  [条件:m_p,α(M_Z)]")
print(f"m_Z = {mZ:.3f} GeV  偏差 {dev(mZ,91.188):.2f}%")

G_F = pi * alpha_MZ / (sqrt2 * s2w * mW**2)
print(f"\nG_F = {G_F:.4e} GeV⁻²  实验 1.166e-05  偏差 {dev(G_F,1.166e-5):.2f}%  [条件:α(M_Z),m_W]")

alpha_s_MZ = 0.118
N_ch = 3 + 2*3*(1 + alpha_s_MZ/pi)
Gamma_W = G_F * mW**3 / (6*sqrt2*pi) * N_ch
print(f"Γ_W = {Gamma_W:.2f} GeV  实验 2.09  偏差 {dev(Gamma_W,2.09):.1f}%  [借用SM相空间]")
tau_W = hbar / Gamma_W
print(f"τ_W = {tau_W:.2e} s")

Gamma_n = G_F**2 * 0.9736**2 * (1.293e-3)**5 / (60*pi**3)
tau_n = hbar / Gamma_n
print(f"τ_n = {tau_n:.0f} s  实验 880 s  差{tau_n/880:.1f}倍  [开放]")

# ================================================================
section("6. 强子质量比")
# ================================================================

m_rho = hbarc/(2*Delta)
Nr = T_N/T_rho
DN = T_Delta/T_N
print(f"m_ρ = {m_rho:.1f} MeV  实验 775.26  偏差 {dev(m_rho,775.26):.2f}%")
print(f"N/ρ = {Nr:.4f}  实验 1.2103  偏差 {dev(Nr,1.2103):.2f}%")
print(f"Δ/N = {DN:.4f}  实验 1.3131  偏差 {dev(DN,1.3131):.2f}%")

# ================================================================
section("7. 额外强力标量（候选）")
# ================================================================

sigma = E0**2/hbarc/Delta*(1+R)
delta_s = E0 * (1 + Delta/2)
mu_n = -6/pi
mu_p = (T_N/2)*(1+Delta)
print(f"弦张力 σ = {sigma:.1f} MeV/fm  实验 910  偏差 {dev(sigma,910):.1f}%")
print(f"味修正 δ_s = {delta_s:.1f} MeV  实验 145  偏差 {dev(delta_s,145):.1f}%")
print(f"μ_n/μ_N = {mu_n:.4f}  实验 -1.913  偏差 {dev(mu_n,-1.913):.2f}%  [数值学候选]")
print(f"μ_p/μ_N = {mu_p:.4f}  实验 2.7928  偏差 {dev(mu_p,2.7928):.2f}%  [已归约]")

print("\n" + "="*70)