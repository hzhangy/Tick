#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_quark_mass_coefficients.py
尝试为绝对夸克质量中的数值系数寻找 NEA 拓扑来源。

候选表达式：
  m_u = m_e / R                （R 的倒数，R 是几何投影残差）
  m_d = m_e · 3π              （3 个空间方向 × π）
  m_s = m_d · 2π²             （2 个环路 × π²）
  m_c = m_μ · 12              （12 = 3D 亲吻数 / C8 边数）
  m_b = m_c · (10/3)          （10 = Stride-10，3 = K4 圈数/颜色数）
  m_t = m_e · α⁻² · 18        （18 = 3 个 K4 的总边数 3×6）

所有系数都尝试从 K4/C8、带宽、寻址等已有结构解读。
"""
import math

pi = math.pi
sqrt3 = math.sqrt(3)
Delta = 1 - sqrt3/2
R = 1/(1+pi)
alpha_inv = 25*sqrt3*pi + 1
U_weak = 10*sqrt3

# 实验锚点
m_e = 0.51099895      # MeV
m_mu_exp = 105.658    # MeV

# 实验夸克质量（MS-bar，MeV）
m_u_exp = 2.16
m_d_exp = 4.67
m_s_exp = 93.4
m_c_exp = 1270        # MeV
m_b_exp = 4180        # MeV
m_t_exp = 172760      # MeV

def dev(calc, exp):
    return abs(calc-exp)/abs(exp)*100

print("="*72)
print("  绝对夸克质量系数：拓扑来源探索")
print("="*72)

# ---------- m_u：候选 m_e / R ----------
m_u = m_e / R
print(f"\n[m_u]  m_e / R")
print(f"  R = 1/(1+π) = {R:.6f}")
print(f"  1/R = 1+π = {1/R:.6f}")
print(f"  m_u = {m_u:.3f} MeV  实验 {m_u_exp}  偏差 {dev(m_u,m_u_exp):.2f}%")
print(f"  拓扑含义：上夸克质量 = 电子质量 × 投影残差的倒数")
print(f"  R 是 1D→2D 投影阻力，1/R 是 2D→1D 增强因子")

# ---------- m_d：候选 m_e · 3π ----------
m_d = m_e * 3 * pi
print(f"\n[m_d]  m_e · 3π")
print(f"  m_d = {m_d:.3f} MeV  实验 {m_d_exp}  偏差 {dev(m_d,m_d_exp):.2f}%")
print(f"  拓扑含义：3 个空间方向 × 相位环 π")

# ---------- m_s：候选 m_d · 2π² ----------
m_s = m_d * 2 * pi**2
print(f"\n[m_s]  m_d · 2π²")
print(f"  m_s = {m_s:.2f} MeV  实验 {m_s_exp}  偏差 {dev(m_s,m_s_exp):.2f}%")
print(f"  拓扑含义：两个环路 × π²（2D 编织面的两个独立相位环）")

# ---------- m_c：候选 m_μ · 12 ----------
# 使用 NEA 自身的 m_μ 以保持零实验输入
m_mu_nea = m_e * 66 * pi * (1 - R/10)
m_c = m_mu_nea * 12 / 1000  # MeV -> GeV
m_c_MeV = m_c * 1000
print(f"\n[m_c]  m_μ · 12")
print(f"  m_μ(NEA) = {m_mu_nea:.2f} MeV")
print(f"  m_c = {m_c_MeV:.1f} MeV  实验 {m_c_exp}  偏差 {dev(m_c_MeV,m_c_exp):.2f}%")
print(f"  拓扑含义：12 = 3D 亲吻数 = C8 边数 = K12 顶点数")
print(f"  粲夸克是第二代：亲吻数使第二代锁定通道饱和")

# ---------- m_b：候选 m_c · (10/3) ----------
m_b = m_c * (10/3) * 1000  # MeV
print(f"\n[m_b]  m_c · (10/3)")
print(f"  m_b = {m_b:.1f} MeV  实验 {m_b_exp}  偏差 {dev(m_b,m_b_exp):.2f}%")
print(f"  拓扑含义：10 = Stride-10 寻址宽度，3 = K4 圈数/颜色数")

# ---------- m_t：候选 m_e · α⁻² · 18 ----------
m_t = m_e * alpha_inv**2 * 18 / 1000 * 1000  # MeV
print(f"\n[m_t]  m_e · α⁻² · 18")
print(f"  m_t = {m_t:.0f} MeV  实验 {m_t_exp}  偏差 {dev(m_t,m_t_exp):.2f}%")
print(f"  拓扑含义：α⁻² 是电磁交换率平方，18 = 3 个 K4 的总边数 3×6")

print("\n" + "="*72)
print("  总结")
print("="*72)
print(f"""
系数拓扑来源候选：
  1/R          ← 投影残差倒数（唯一出现在 m_u）
  3π           ← 3 个空间方向 × π
  2π²          ← 两个相位环
  12           ← 亲吻数 / C8 边数 / K12 顶点数
  10/3         ← Stride-10 ÷ K4 圈数
  18           ← 3 个 K4 总边数

状态：
  ✓ 每个系数都能指认一个 NEA 拓扑量
  ✓ 偏差全部 <4%
  ⚠ 但“指认”还不是“推导”：尚需证明这些拓扑量以该组合出现
    而不是事后拼凑。下一步需要从 K4/C8 动力学推出这些组合。
""")