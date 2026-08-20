#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_kappa_lock.py

Close the last gap: kappa is the locked node's bandwidth deficit.

微观规则：
  K4 闭环锁定 → f_ext = sqrt(3)/2
  带宽缺失量 Delta = 1 - f_ext = 1 - sqrt(3)/2

宏观引力：
  锁定节点密度 rho_lock
  场方程：nabla^2(1-f_ext) = -kappa rho_lock
  其中 kappa = Delta (每个锁定节点的源强度)

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

# ========================== 微观常数 ==========================
DELTA = 1.0 - np.sqrt(3)/2      # 单个锁定节点的带宽缺失量
print("="*72)
print("  NEA gravity with kappa locked by K4 bandwidth deficit")
print("="*72)
print(f"  Locked node f_ext     = sqrt(3)/2 = {np.sqrt(3)/2:.6f}")
print(f"  Bandwidth deficit     = 1 - f_ext = {DELTA:.6f}")

# ========================== 宏观质量 ==========================
# 一个 K4 锁定核包含 N_lock 个锁定节点
N_lock = 27

# 总源强度 A_gr 是唯一确定的，不是手动参数
A_gr = N_lock * DELTA / (4.0 * np.pi)
print(f"\n  Locked nodes in core  = {N_lock}")
print(f"  Total source strength A_gr = N_lock*Delta/(4pi)")
print(f"                            = {A_gr:.6f}")

# ========================== 模拟参数 ==========================
L = 2000.0
cx = cy = L / 2.0
ds = 0.05

# 锁定核中心 3x3x3 的节点位置（相对中心）
offsets = [(dx, dy, dz) for dx in (-1,0,1) for dy in (-1,0,1) for dz in (-1,0,1)]

def f_ext_at(x, y):
    """f_ext 由 27 个锁定节点的 1/r 场叠加决定，零自由参数。"""
    # 每个锁定节点贡献 Δ/(4π r)
    s = 0.0
    for dx, dy, dz in offsets:
        r = np.sqrt((x-(cx+dx))**2 + (y-(cy+dy))**2 + dz**2) + 1e-9
        s += DELTA / (4.0 * np.pi * r)
    # 弱场下 f_ext = 1 - s
    return 1.0 - s

def n_full(x, y):
    """完整 NEA 度规的有效折射率：n = 1/f_ext^2"""
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
print(f"  Expected dispersion < 0.05 for clean alpha ~ 1/b")

# ========================== 时间膨胀 ==========================
print(f"\n{'='*72}")
print("  Gravitational time dilation: dtau/dt = f_ext")
print("="*72)
print(f"  {'r':<8}{'f_ext':<14}{'GR weak field':<16}{'ratio'}")
for r in [50, 100, 200, 400]:
    fe = f_ext_at(cx, cy+r)  # f_ext at radius r along y-axis
    gr_val = 1.0 - A_gr/r   # GR weak field dtau/dt
    ratio = fe/gr_val
    print(f"  {r:<8}{fe:<14.8f}{gr_val:<16.8f}{ratio:<8.6f}")

# ========================== 绘图 ==========================
fig, axes = plt.subplots(2,2, figsize=(14, 12))

# (a) f_ext radial profile
ax = axes[0,0]
r_arr = np.linspace(5, 500, 300)
fe_arr = [f_ext_at(cx, cy+r) for r in r_arr]
ax.plot(r_arr, fe_arr, 'r-', lw=2, label='f_ext (numerical)')
ax.plot(r_arr, 1.0 - A_gr/r_arr, 'b--', lw=1.5, label='1 - A_gr/r')
ax.set_xlabel('r')
ax.set_ylabel('f_ext')
ax.set_title('Radial f_ext profile')
ax.legend()
ax.grid(alpha=0.3)

# (b) alpha*b
ax = axes[0,1]
ax.plot(b_vals, alpha_b_vals, 'ro-', ms=6, lw=2, label='numerical alpha*b')
ax.axhline(-4*A_gr, color='blue', ls='--', lw=1.5, label='GR: -4A_gr')
ax.set_xlabel('b')
ax.set_ylabel('alpha*b')
ax.set_title(f'Light bending (A_gr={A_gr:.5f})')
ax.legend()
ax.grid(alpha=0.3)

# (c) bending angles
ax = axes[1,0]
b_fine = np.linspace(15, 250, 40)
alpha_fine = [abs(shoot_ray(cy+b)) for b in b_fine]
ax.loglog(b_fine, alpha_fine, 'ro-', ms=4, label='numerical')
ax.loglog(b_fine, 4*A_gr/b_fine, 'b--', lw=2, label='GR: 4A/b')
ax.set_xlabel('b')
ax.set_ylabel('|alpha| (rad)')
ax.set_title('Bending angle scaling')
ax.legend()
ax.grid(alpha=0.3)

# (d) time dilation
ax = axes[1,1]
ax.plot(r_arr, 1.0 - A_gr/r_arr, 'b-', lw=2, label='GR: 1-A/r')
ax.plot(r_arr, [f_ext_at(cx, cy+r) for r in r_arr], 'r--', lw=1.5, label='NEA: f_ext')
ax.set_xlabel('r')
ax.set_ylabel('d tau / dt')
ax.set_title('Gravitational time dilation')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('nea_kappa_lock.png', dpi=150)
print("\nSaved: nea_kappa_lock.png")

print(f"\n{'='*72}")
print("  Summary")
print("="*72)
print(f"""
  kappa = Delta = 1 - sqrt(3)/2 = {DELTA:.6f}
  A_gr  = N_lock * kappa / (4*pi)
        = {N_lock} * {DELTA:.6f} / (4*pi)
        = {A_gr:.6f}

  弱场引力完全由微观 K4 锁定深度决定，无自由参数。

  Numerical alpha*b ratio to GR = {mean_ab/(-4*A_gr):.4f}
  Dispersion                     = {disp:.6f}

  NEA degree of freedom count: 0 free parameters for weak-field gravity.
""")
print("="*72)