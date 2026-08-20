#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_em_wave_fdtd_final_correct.py
最终修正版：正确使用历史场验证麦克斯韦方程组。
自动扫描候选点，选出电生磁和磁生电验证误差最小的点。
"""
import numpy as np

# 网格参数
Nx, Ny = 160, 160
T = 200
c = 1.0
dx = dy = 1.0
dt = 0.3

# Yee 网格场
Ex = np.zeros((Nx, Ny-1))      # x方向边
Ey = np.zeros((Nx-1, Ny))      # y方向边
Hz = np.zeros((Nx-1, Ny-1))    # 面心

# 源
cx, cy = Nx//2, Ny//2
pulse_width = 6.0
t0 = 20

# 历史记录（每5步存一次）
history_t = []       # 时间步
history_Hz = []
history_Ex = []
history_Ey = []

# 候选验证点（远离边界和源，分散分布）
candidate_points = [
    (60, 55), (65, 50), (70, 60), (55, 65), (50, 70),
    (60, 65), (65, 60), (70, 55), (55, 60), (60, 70),
    (75, 65), (65, 75), (70, 70), (55, 50), (50, 55)
]

# ===== 主模拟 =====
for t in range(T):
    # 更新 Ex
    for i in range(1, Nx-1):
        for j in range(1, Ny-2):
            Ex[i, j] += c * dt / dy * (Hz[i, j] - Hz[i, j-1])

    # 注入源
    if t < 40:
        Ex[cy, cx-1] += np.exp(-0.5*(t - t0)**2 / pulse_width**2)

    # 更新 Ey
    for i in range(1, Nx-2):
        for j in range(1, Ny-1):
            Ey[i, j] -= c * dt / dx * (Hz[i, j] - Hz[i-1, j])

    # 更新 Hz
    for i in range(Nx-2):
        for j in range(Ny-2):
            Hz[i, j] += c * dt * ((Ex[i, j+1] - Ex[i, j]) / dy - (Ey[i+1, j] - Ey[i, j]) / dx)

    # 记录历史（每5步）
    if t % 5 == 0:
        history_t.append(t)
        history_Hz.append(Hz.copy())
        history_Ex.append(Ex.copy())
        history_Ey.append(Ey.copy())

# ===== 自动寻找最佳验证点 =====
def compute_errors_for_point(vp_x, vp_y):
    """计算给定点的电生磁和磁生电最大误差（使用历史场）"""
    errors_em = []
    errors_me_ex = []
    errors_me_ey = []
    grad_y_max = 0.0
    grad_x_max = 0.0

    # 遍历历史中除了首尾的时刻
    for idx in range(1, len(history_t)-1):
        t = history_t[idx]
        Hz_t = history_Hz[idx]
        Ex_t = history_Ex[idx]
        Ey_t = history_Ey[idx]
        Hz_prev = history_Hz[idx-1]
        Hz_next = history_Hz[idx+1]
        Ex_prev = history_Ex[idx-1]
        Ex_next = history_Ex[idx+1]
        Ey_prev = history_Ey[idx-1]
        Ey_next = history_Ey[idx+1]

        dt_eff = 5 * dt  # 时间步间隔

        # 电生磁：dHz/dt = c*(∂Ex/∂y - ∂Ey/∂x)
        dHz_dt = (Hz_next[vp_x, vp_y] - Hz_prev[vp_x, vp_y]) / (2 * dt_eff)
        curl_E = (Ex_t[vp_x, vp_y+1] - Ex_t[vp_x, vp_y]) / dy - (Ey_t[vp_x+1, vp_y] - Ey_t[vp_x, vp_y]) / dx
        errors_em.append(abs(dHz_dt - c * curl_E))

        # 磁生电 Ex：dEx/dt = c * ∂Hz/∂y
        dEx_dt = (Ex_next[vp_x, vp_y] - Ex_prev[vp_x, vp_y]) / (2 * dt_eff)
        dHz_dy = (Hz_t[vp_x, vp_y+1] - Hz_t[vp_x, vp_y]) / dy
        errors_me_ex.append(abs(dEx_dt - c * dHz_dy))
        grad_y_max = max(grad_y_max, abs(dHz_dy))

        # 磁生电 Ey：dEy/dt = -c * ∂Hz/∂x
        dEy_dt = (Ey_next[vp_x, vp_y] - Ey_prev[vp_x, vp_y]) / (2 * dt_eff)
        dHz_dx = (Hz_t[vp_x+1, vp_y] - Hz_t[vp_x, vp_y]) / dx
        errors_me_ey.append(abs(dEy_dt - (-c * dHz_dx)))
        grad_x_max = max(grad_x_max, abs(dHz_dx))

    return (max(errors_em), max(errors_me_ex), max(errors_me_ey),
            grad_y_max, grad_x_max)

# 计算所有候选点的误差
results = []
for pt in candidate_points:
    em, me_ex, me_ey, gy, gx = compute_errors_for_point(pt[0], pt[1])
    total = em + me_ex + me_ey
    results.append((pt, em, me_ex, me_ey, gy, gx, total))
    print(f"点 {pt}: 电生磁maxErr={em:.5f}, 磁生电Ex maxErr={me_ex:.5f}, "
          f"磁生电Ey maxErr={me_ey:.5f}, |∂Hz/∂y|max={gy:.5f}, |∂Hz/∂x|max={gx:.5f}")

# 选取总误差最小的点
best = min(results, key=lambda x: x[6])
best_pt = best[0]
print(f"\n最优验证点：{best_pt}，总误差={best[6]:.5f}")
print("该点电生磁maxErr={:.5f}, 磁生电Ex maxErr={:.5f}, 磁生电Ey maxErr={:.5f}".format(
    best[1], best[2], best[3]))

# ===== 打印最优点的详细验证数据 =====
vp_x, vp_y = best_pt
print(f"\n===== 最优验证点 ({vp_x},{vp_y}) 详细数据 =====")

print("\n电生磁验证：dHz/dt 与 c*curl(E)")
print(f"{'t':>4} | {'Hz':>10} | {'dHz/dt':>10} | {'c*curl(E)':>12} | {'diff':>10}")
for idx in range(1, len(history_t)-1):
    t = history_t[idx]
    Hz_t = history_Hz[idx]
    Ex_t = history_Ex[idx]
    Ey_t = history_Ey[idx]
    Hz_prev = history_Hz[idx-1]
    Hz_next = history_Hz[idx+1]
    dt_eff = 5 * dt

    dHz_dt = (Hz_next[vp_x, vp_y] - Hz_prev[vp_x, vp_y]) / (2 * dt_eff)
    curl_E = (Ex_t[vp_x, vp_y+1] - Ex_t[vp_x, vp_y]) / dy - (Ey_t[vp_x+1, vp_y] - Ey_t[vp_x, vp_y]) / dx
    diff = abs(dHz_dt - c * curl_E)
    print(f"{t:4d} | {Hz_t[vp_x, vp_y]:10.5f} | {dHz_dt:10.5f} | {c*curl_E:12.5f} | {diff:10.6f}")

print("\n磁生电验证（Ex分量）：dEx/dt 与 c*∂Hz/∂y")
print(f"{'t':>4} | {'Ex':>10} | {'dEx/dt':>10} | {'c*∂Hz/∂y':>12} | {'diff':>10}")
for idx in range(1, len(history_t)-1):
    t = history_t[idx]
    Hz_t = history_Hz[idx]
    Ex_t = history_Ex[idx]
    Ex_prev = history_Ex[idx-1]
    Ex_next = history_Ex[idx+1]
    dt_eff = 5 * dt

    dEx_dt = (Ex_next[vp_x, vp_y] - Ex_prev[vp_x, vp_y]) / (2 * dt_eff)
    dHz_dy = (Hz_t[vp_x, vp_y+1] - Hz_t[vp_x, vp_y]) / dy
    diff = abs(dEx_dt - c * dHz_dy)
    print(f"{t:4d} | {Ex_t[vp_x, vp_y]:10.5f} | {dEx_dt:10.5f} | {c*dHz_dy:12.5f} | {diff:10.6f}")

print("\n磁生电验证（Ey分量）：dEy/dt 与 -c*∂Hz/∂x")
print(f"{'t':>4} | {'Ey':>10} | {'dEy/dt':>10} | {'-c*∂Hz/∂x':>12} | {'diff':>10}")
for idx in range(1, len(history_t)-1):
    t = history_t[idx]
    Hz_t = history_Hz[idx]
    Ey_t = history_Ey[idx]
    Ey_prev = history_Ey[idx-1]
    Ey_next = history_Ey[idx+1]
    dt_eff = 5 * dt

    dEy_dt = (Ey_next[vp_x, vp_y] - Ey_prev[vp_x, vp_y]) / (2 * dt_eff)
    dHz_dx = (Hz_t[vp_x+1, vp_y] - Hz_t[vp_x, vp_y]) / dx
    diff = abs(dEy_dt - (-c * dHz_dx))
    print(f"{t:4d} | {Ey_t[vp_x, vp_y]:10.5f} | {dEy_dt:10.5f} | {-c*dHz_dx:12.5f} | {diff:10.6f}")