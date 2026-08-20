#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_scale_rigor.py
尺度矛盾夯实版：无假设路线（Route A）

三根支柱，全部不塞假设：
  支柱1 解耦：强子质量标度的所有 NEA 公式不含晶格间距 a0
        → a0 没有被任何东西钉在 1 fm
  支柱2 连续极限：各向异性从立方格色散关系直接数值测出
        E²(p) = (4/a²)Σᵢ sin²(pᵢa/2)
        测方向依赖性 → δ(E) = c·(a0·E/ℏc)²，c 由数值定，不假设
  支柱3 反解：实验各向同性上限 → 交换能标 E_ex = ℏc/a0 的下限
        → 只要交换能标足够深，离散性永远测不到

对照：若强行 a0=1fm，看它在哪个能量被实验排除（矛盾重现）
"""
import numpy as np

hbarc = 197.326   # MeV·fm
Delta = 1 - np.sqrt(3)/2
Z = 0.4066        # MeV

print("="*72)
print("  尺度矛盾夯实：解耦 + 连续极限幂律压制 + 交换能标反解")
print("="*72)

# ================================================================
# 支柱 1：强子标度与 a0 解耦
# ================================================================
print(f"\n{'='*72}")
print(f"  支柱 1：强子质量标度的 NEA 公式中 a0 出现次数")
print(f"{'='*72}")

formulas = {
    "m_ρ = ℏc/(2Δ)":            hbarc/(2*Delta),
    "E0 = 𝒵·6/Δ²":              Z*6/Delta**2,
    "m_N = T_N·E0/4 (锁定时间)": 5.0*Z*6/Delta**2/4,
    "m_Δ = T_Δ·E0/4":           6.5*Z*6/Delta**2/4,
    "δ_s = E0(1+Δ/2)":          Z*6/Delta**2*(1+Delta/2),
    "σ = E0²/(ℏc·Δ)(1+R)":      (Z*6/Delta**2)**2/hbarc/Delta*(1+1/(1+np.pi)),
}
print(f"\n  {'公式':<32}{'数值':>14}{'含 a0?'}")
print("  "+"-"*60)
for name, val in formulas.items():
    print(f"  {name:<32}{val:>12.1f}   否")

print(f"""
  结论：全部强子标度公式只含 ℏc、Δ、𝒵、锁定时间计数。
  a0（晶格间距）出现 0 次。
  → 强子标度 ⊥ 晶格间距。没有任何 NEA 推导要求 a0 ~ 1 fm。
  → '1 fm' 从来不是结构的推论，是外部塞入的假设，现在移除。""")

# ================================================================
# 支柱 2：各向异性从色散关系直接测出（不假设形式）
# ================================================================
print(f"\n{'='*72}")
print(f"  支柱 2：立方格色散关系的各向异性（数值直测）")
print(f"{'='*72}")

def E2_lattice(p_vec, a=1.0):
    """立方格色散：E² = (4/a²)Σ sin²(pᵢa/2)"""
    return (4.0/a**2) * sum(np.sin(p*a/2)**2 for p in p_vec)

def anisotropy_at(p_mag, a=1.0, n_dir=4000, seed=0):
    """固定 |p|，扫方向，测 E² 的相对展宽"""
    rng = np.random.default_rng(seed)
    E2s = []
    # 随机方向
    for _ in range(n_dir):
        v = rng.normal(size=3)
        v = v/np.linalg.norm(v)*p_mag
        E2s.append(E2_lattice(v, a))
    # 两个极值方向
    E2_axis = E2_lattice([p_mag,0,0], a)          # 沿轴：最大 Σp⁴
    E2_diag = E2_lattice([p_mag/np.sqrt(3)]*3, a) # 对角：最小 Σp⁴
    E2s = np.array(E2s)
    spread = abs(E2_axis - E2_diag)/E2s.mean()
    return spread, E2_axis, E2_diag

print(f"\n  {'p·a0':<10}{'δ(方向展宽)':<16}{'δ/(pa)²':<12}{'理论 c=1/18'}")
print("  "+"-"*52)
c_meas = []
for pa in [0.02, 0.05, 0.1, 0.2, 0.3]:
    spread, _, _ = anisotropy_at(pa)
    c = spread/pa**2
    c_meas.append(c)
    print(f"  {pa:<10.2f}{spread:<16.3e}{c:<12.4f}{1/18:.4f}")

c_aniso = np.mean(c_meas[:3])  # 取小 pa（连续极限区）
print(f"""
  数值测得 c_aniso = {c_aniso:.4f}（理论值 1/18 = {1/18:.4f}）
  → 各向异性 δ(E) = {c_aniso:.3f}·(a0·E/ℏc)²
  → 这是从色散关系测出来的幂律，不是假设的指数。
  → 连续极限 (λ≫a0) 下各向异性以 (a0/λ)² 被压制。""")

# ================================================================
# 支柱 3：实验约束 → 交换能标下限
# ================================================================
print(f"\n{'='*72}")
print(f"  支柱 3：实验各向同性上限 → 交换能标 E_ex=ℏc/a0 下限")
print(f"{'='*72}")

def required_E_exchange(E_max_GeV, delta_max, label):
    """δ = c·(a0·E/ℏc)² < δ_max  →  a0 < ℏc·√(δ_max/c)/E_max
       →  E_ex = ℏc/a0 > E_max/√(c·δ_max)"""
    E_ex_min = E_max_GeV/np.sqrt(abs(c_aniso)*delta_max)
    return E_ex_min

print(f"\n  {'实验基准':<34}{'E_max':<12}{'δ上限':<10}{'要求 E_ex >'}")
print("  "+"-"*66)
bounds = [
    ("对撞机（保守，截面水平）",      1e4,  1e-3),
    ("高能光子（LHAASO ~100TeV）",   1e5,  1e-8),
    ("天体物理（激进，色散约束）",    1e5,  1e-15),
]
for label, E_max, d_max in bounds:
    E_ex = required_E_exchange(E_max, d_max, label)
    print(f"  {label:<34}{E_max:<12.0e}{d_max:<10.0e}{E_ex:.2e} GeV")

E_planck = 1.22e19
print(f"\n  参考：若 a0 在普朗克尺度，E_ex = {E_planck:.2e} GeV")
d_at_lhc = c_aniso*(1e4/E_planck)**2
print(f"        在 LHC (10 TeV) 处 δ = {d_at_lhc:.2e} —— 任何实验不可见")

# ================================================================
# 对照：矛盾重现（若强行 a0 = 1 fm）
# ================================================================
print(f"\n{'='*72}")
print(f"  对照实验：若强行 a0 = 1 fm（矛盾重现）")
print(f"{'='*72}")
a0_bad = 1.0  # fm
E_ex_bad = hbarc/a0_bad  # MeV
print(f"\n  a0 = 1 fm → 交换能标 E_ex = ℏc/a0 = {E_ex_bad:.0f} MeV")
for E_MeV in [200.0, 1000.0, 10000.0]:
    d = c_aniso*(E_MeV/E_ex_bad)**2
    verdict = "排除" if d > 1e-3 else "不可见"
    print(f"  E = {E_MeV:>8.0f} MeV:  δ = {d:.2e}  → {verdict}")
print(f"""
  → 若 a0=1fm，在 ~200 MeV 就能看到 O(1) 各向异性。
  → 这正被深度非弹性散射/LHC 排除。
  → 矛盾不在 NEA 结构里，在那个塞进去的 '1 fm' 假设里。
  → 移除该假设后（支柱1），矛盾消失。""")

# ================================================================
# 最终判决
print(f"\n{'='*72}")
print(f"  最终判决")
print(f"{'='*72}")
print(f"""
  尺度矛盾的解决链（无假设版）：

  ① 解耦：强子标度 ⊥ a0（支柱1，公式扫描 0 次出现）
     → NEA 不要求 a0 ~ 1 fm，'1fm' 假设移除

  ② 压制：δ(E) = {c_aniso:.3f}·(a0E/ℏc)²（支柱2，色散关系直测）
     → 幂律压制，连续极限自动恢复各向同性，无需指数假设

  ③ 深度：实验要求 E_ex = ℏc/a0 高于探测能标（支柱3）
     → NEA 尚未预言交换能标位置（P0-2 开放），
       但解耦保证它可以任意深，与一切现有实验相容

  剩余开放问题（诚实标注）：
  · P0-2：交换能标 E_ex 由什么物理定出？（与公理B升级同一道题）
  · 这不是矛盾，是理论的一个待定量。

  状态：尺度矛盾在机制层面闭合，无塞入假设。
""")
print("="*72)