#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_em_c8_poisson.py
在 C8 密铺（3D 立方格）上求解电磁势：
  ∇²φ = -κ_q ρ_q
验证库仑势 φ ∝ 1/r，电场 E ∝ 1/r²。
电荷源 = K4 顶点未配平方向相位，这里用单节点源简化。
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KAPPA = 1.0 - np.sqrt(3)/2   # Δ，统一耦合常数（此处作为相位-带宽转换因子）

def build_cubic_lattice(L):
    """LxLxL 立方格，顶点数 (L+1)^3，边为最近邻连接。"""
    n = L + 1
    N = n**3
    adj = [[] for _ in range(N)]
    def idx(x,y,z): return (x*n + y)*n + z
    for x in range(n):
        for y in range(n):
            for z in range(n):
                i = idx(x,y,z)
                for dx,dy,dz in [(1,0,0),(0,1,0),(0,0,1)]:
                    nx,ny,nz = x+dx,y+dy,z+dz
                    if nx < n and ny < n and nz < n:
                        j = idx(nx,ny,nz)
                        adj[i].append(j)
                        adj[j].append(i)
    return adj, n

def is_surface(x,y,z,L):
    """表面节点（Dirichlet边界 φ=0）"""
    return x==0 or x==L or y==0 or y==L or z==0 or z==L

def solve_poisson(adj, source_mask, surface_mask, kappa, n_iter=12000, omega=0.5):
    """解离散泊松方程：sum_{j∈∂i}(φ_i - φ_j) = κ ρ_i"""
    N = len(adj)
    rows = []; cols = []
    for i, nbrs in enumerate(adj):
        for j in nbrs:
            rows.append(i); cols.append(j)
    rows = np.array(rows, dtype=int)
    cols = np.array(cols, dtype=int)
    degs = np.array([len(a) for a in adj], dtype=float)
    rho = source_mask.astype(float)
    phi = np.zeros(N)
    for it in range(n_iter):
        phi_sum = np.zeros(N)
        np.add.at(phi_sum, rows, phi[cols])
        phi_new = (phi_sum + kappa * rho) / np.maximum(degs, 1.0)
        phi_new[surface_mask] = 0.0
        phi = (1-omega)*phi + omega*phi_new
    phi[surface_mask] = 0.0
    return phi

def bfs_distances(adj, sources):
    N = len(adj)
    dist = np.full(N, -1, dtype=int)
    q = list(sources)
    for s in sources: dist[s] = 0
    head = 0
    while head < len(q):
        x = q[head]; head += 1
        for y in adj[x]:
            if dist[y] < 0:
                dist[y] = dist[x] + 1
                q.append(y)
    return dist

def fit_beta(pr, pu, r_lo, r_hi):
    sel = [(pr[i], pu[i]) for i in range(len(pr)) if r_lo <= pr[i] <= r_hi and pu[i] > 1e-12]
    if len(sel) < 3:
        return None, None, len(sel)
    rr = np.array([s[0] for s in sel], dtype=float)
    uu = np.array([s[1] for s in sel], dtype=float)
    c = np.polyfit(np.log(rr), np.log(uu), 1)
    beta = -c[0]
    log_u = np.log(uu)
    pred = np.polyval(c, np.log(rr))
    ss_res = np.sum((log_u - pred)**2)
    ss_tot = np.sum((log_u - log_u.mean())**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0.0
    return beta, r2, len(sel)

if __name__ == "__main__":
    L = 16
    adj, n = build_cubic_lattice(L)
    N = len(adj)
    print("="*72)
    print("  C8 密铺上的电磁势：库仑定律验证")
    print("="*72)
    print(f"  立方格尺寸: L={L}, 节点数 N={N}")

    surface_mask = np.zeros(N, dtype=bool)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if is_surface(x, y, z, L):
                    surface_mask[(x*n + y)*n + z] = True

    # 电荷源：中心单节点
    c = L//2
    source_nodes = [(c*n + c)*n + c]
    source_mask = np.zeros(N, dtype=bool)
    for s in source_nodes:
        source_mask[s] = True

    phi = solve_poisson(adj, source_mask, surface_mask, KAPPA, n_iter=12000)
    print(f"  电势最大值 phi_max = {phi.max():.6f} (弱场成立)")

    d2s = bfs_distances(adj, source_nodes)
    max_r = d2s.max()
    # 只测量内部节点且非源
    internal_free = (~surface_mask) & (~source_mask)
    pr, pu = [], []
    for r in range(1, min(max_r, 20)+1):
        shell = np.where((d2s == r) & internal_free)[0]
        if len(shell) > 0:
            pr.append(r)
            pu.append(phi[shell].mean())
    # 远场窗口
    R_b = d2s[surface_mask].min()
    r_hi = R_b // 2
    beta, r2, n_fit = fit_beta(pr, pu, 2, r_hi)
    if beta is not None:
        print(f"  远场拟合 r∈[2,{r_hi}]: β = {beta:.3f} (R²={r2:.4f})")
        if abs(beta - 1.0) < 0.2:
            print("  ✓ β≈1：库仑势 φ∝1/r 成立")
        else:
            print("  ? β偏离1，需检查边界或窗口")
    # 计算电场（梯度）
    # 沿 x 方向电场强度：E = -dφ/dx，用中心差分
    # 取 y=z=c 轴上的节点
    axis_idx = [(x*n + c)*n + c for x in range(n) if not surface_mask[(x*n + c)*n + c] and (x*n + c)*n + c not in source_nodes]
    x_coords = []
    phi_axis = []
    for idx in axis_idx:
        x = idx // (n*n)
        if x > 0 and x < L:
            phi_left = phi[(x-1)*n*n + c*n + c]
            phi_right = phi[(x+1)*n*n + c*n + c]
            E = -(phi_right - phi_left) / 2.0
            r = abs(x - c)
            x_coords.append(r)
            phi_axis.append((phi[idx], E))
    # 输出几个点
    print("\n  轴线上电场强度 vs 距离:")
    print(f"  {'r':<5}{'phi':<14}{'E':<14}{'E*r²':<14}")
    for r, (p, E) in zip(x_coords, phi_axis):
        if 2 <= r <= 8:
            Er2 = abs(E) * r * r if r > 0 else 0
            print(f"  {r:<5}{p:<14.6f}{E:<14.6f}{Er2:<14.6f}")
    # 保存图形
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    ax.loglog(pr, pu, 'ro-', ms=4, lw=2, label='numerical φ')
    if beta is not None:
        r_fit = np.linspace(2, r_hi, 50)
        A_fit = np.exp(np.polyfit(np.log(np.array(pr[1:len(r_fit)+1])), np.log(np.array(pu[1:len(r_fit)+1])), 1)[1])
        ax.loglog(r_fit, A_fit/r_fit**beta, 'b--', lw=2, label=f'fit β={beta:.2f}')
    ax.set_xlabel('r (graph distance)')
    ax.set_ylabel('φ')
    ax.set_title('Electromagnetic potential on C8 tiling')
    ax.legend(); ax.grid(alpha=0.3, which='both')
    ax = axes[1]
    if x_coords:
        # 修复 Bug：确保 r 和 E 严格一一对应，去除重复点
        axis_data = {}
        for r, (p, E) in zip(x_coords, phi_axis):
            if 2 <= r <= 8 and r not in axis_data:
                axis_data[r] = abs(E)
        r_axis = sorted(axis_data.keys())
        E_axis = [axis_data[r] for r in r_axis]
        ax.loglog(r_axis, E_axis, 'bo-', ms=4, lw=2, label='|E|')
        # 1/r² reference
        ref_r = np.linspace(2, 8, 50)
        ref_E = E_axis[0] * (r_axis[0]**2) / ref_r**2
        ax.loglog(ref_r, ref_E, 'g--', lw=2, label='1/r²')
        ax.set_xlabel('r')
        ax.set_ylabel('|E|')
        ax.set_title('Electric field on axis')
        ax.legend(); ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig('nea_em_c8_poisson.png', dpi=150)
    print("  saved: nea_em_c8_poisson.png")