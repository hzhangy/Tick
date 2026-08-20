#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_kappa_lock_v2.py
Close the last gap: kappa has correct dimension.

微观规则：
  K4 闭环锁定 → f_ext = sqrt(3)/2
  带宽缺失量 Delta = 1 - f_ext = 1 - sqrt(3)/2

宏观引力：
  锁定节点密度 rho_lock
  场方程：nabla^2(1-f_ext) = -kappa rho_lock
  其中 kappa = Delta * a（a 是外部尺度锚，fm）

因此：
  总源强度 A_gr = N_lock * kappa / (4*pi)
  这是唯一确定的值，不是自由参数。

弱场度规：
  ds^2 = -f_ext^2 dt^2 + dl^2/f_ext^2
  f_ext = 1 - A_gr/r
  → n_eff = 1/f_ext^2 → alpha = 4 A_gr / b

输出：
  1. 理论 A_gr
  2. 数值光线偏折 alpha*b 是否等于 4 A_gr
  3. 引力时间膨胀是否与弱场 GR 一致
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ========================== 物理尺度锚 ==========================
a = 1.0           # fm，外部标定尺度

# ========================== 微观常数 ==========================
DELTA = 1.0 - np.sqrt(3)/2      # 单个锁定节点的带宽缺失量
F_EXT_K4 = np.sqrt(3)/2         # K4 锁定后的外部带宽

print("="*72)
print("  NEA gravity with kappa dimensionally correct")
print("="*72)
print(f"  a (fm)                = {a}")
print(f"  DELTA                 = {DELTA:.6f}")
print(f"  kappa = DELTA * a     = {DELTA*a:.6f} fm")

# ========================== 宏观质量 ==========================
# 一个 K4 锁定核包含 N_lock 个锁定节点
N_lock = 27

# 总源强度 A_gr 是唯一确定的，不是手动参数
A_gr = N_lock * DELTA * a / (4.0 * np.pi)
print(f"\n  Locked nodes in core  = {N_lock}")
print(f"  A_gr = N_lock*kappa/(4pi) = {A_gr:.6f} fm")

# ========================== 模拟参数 ==========================
L = 2000.0
cx = cy = L / 2.0
ds = 0.05

offsets = [(dx, dy, dz) for dx in (-1,0,1) for dy in (-1,0,1) for dz in (-1,0,1)]

def f_ext_at(x, y):
    """f_ext 由 27 个锁定节点的 1/r 场叠加决定，零自由参数。"""
    s = 0.0
    for dx, dy, dz in offsets:
        r = np.sqrt((x-(cx+dx))**2 + (y-(cy+dy))**2 + dz**2) + 1e-9
        s += DELTA * a / (4.0 * np.pi * r)
    return 1.0 - s

def n_full(x, y):
    f = f_ext_at(x, y)
    return 1.0 / (f * f)

def grad_n(x, y, h=0.02):
    nx = (n_full(x+h, y) - n_full(x-h, y)) / (2*h)
    ny = (n_full(x, y+h) - n_full(x, y-h)) / (2*h)
    return nx, ny

def shoot_ray(y_start):
    x = 0.0
    y = float(y_start)
    vx, vy = 1.0, 0.0
    steps = 0
    while x < L-1.0 and 0.0 < y < L:
        gnx, gny = grad_n(x, y)
        n = n_full(x, y)
        px = n * vx
        py = n * vy
        px += gnx * ds
        py += gny * ds
        pnorm = np.sqrt(px*px + py*py)
        vx, vy = px/pnorm, py/pnorm
        x += vx * ds
        y += vy * ds
        steps += 1
        if steps > 200000:
            break
    return np.arctan2(vy, vx)

print(f"\n{'='*72}")
print("  Light bending: numerical alpha*b vs 4*A_gr")
print("="*72)
print(f"  {'b':<8}{'alpha (deg)':<16}{'alpha*b':<16}{'4A/b (deg)':<16}")
b_vals = [20.0, 40.0, 80.0, 120.0, 160.0, 200.0]
alpha_b_vals = []
for b in b_vals:
    alpha = shoot_ray(cy + b)
    alpha_b_vals.append(alpha * b)
    gr_pred = -4.0 * A_gr / b
    print(f"  {b:<8.0f}{np.degrees(alpha):<16.8f}{alpha*b:<16.8f}"
          f"{np.degrees(gr_pred):<16.8f}")

mean_ab = np.mean(alpha_b_vals)
disp = np.std(alpha_b_vals) / abs(mean_ab)
print(f"\n  Numerical alpha*b mean  = {mean_ab:.8f}")
print(f"  Theoretical -4*A_gr    = {-4.0*A_gr:.8f}")
print(f"  Ratio num/theory       = {mean_ab/(-4.0*A_gr):.6f}")
print(f"  Dispersion             = {disp:.6f}")

# ========================== 时间膨胀 ==========================
print(f"\n{'='*72}")
print("  Gravitational time dilation: dtau/dt = f_ext")
print("="*72)
print(f"  {'r':<8}{'f_ext':<14}{'GR weak field':<16}{'ratio'}")
for r in [50, 100, 200, 400]:
    fe = f_ext_at(cx, cy+r)
    gr_val = 1.0 - A_gr/r
    ratio = fe/gr_val
    print(f"  {r:<8}{fe:<14.8f}{gr_val:<16.8f}{ratio:<8.6f}")

print("\nDone.")
print("="*72)