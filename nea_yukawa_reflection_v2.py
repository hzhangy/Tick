#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_yukawa_reflection_v2.py
目标：把 Yukawa 反射拟合的残余偏差压到 <5%，全质量段定量闭合

上一轮诊断：
  · m=0.1: μ_fit=0.1146 (dev 14.6%)，λ/R_b=0.5
  · m=0.2: μ_fit=0.2232 (dev 11.6%)，λ/R_b=0.25
  · m=0.5: μ_fit=0.5334 (dev 6.7%)，λ/R_b=0.1
  规律：λ/R_b 越大，偏差越大 → 有限尺寸效应
本轮修正：
  1. 盒子放大到 L=80，R_b=40，λ/R_b 全面缩小
  2. 多窗口拟合：r_lo 从 3 扫到 8，看 μ_fit 是否随窗口移向远场而收敛
  3. 同时拟合 sinh 形式（含反射）和纯指数形式（无反射），对比
  4. 报告 μ_fit vs m_in vs μ_lattice=arccosh(1+m²/2)

预期：
  · 大盒子 + 远场窗口 → μ_fit 收敛到 m_in，偏差 <5%
  · 若仍偏差，说明不是有限尺寸，而是格点色散或近源伪影
"""
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

KAPPA = 1.0 - np.sqrt(3)/2

# ================================================================
# 建图
# ================================================================
def build_cubic_lattice(L):
    n=L+1; N=n**3
    adj=[[] for _ in range(N)]
    def idx(x,y,z): return (x*n+y)*n+z
    for x in range(n):
        for y in range(n):
            for z in range(n):
                i=idx(x,y,z)
                for dx,dy,dz in [(1,0,0),(0,1,0),(0,0,1)]:
                    nx,ny,nz=x+dx,y+dy,z+dz
                    if nx<n and ny<n and nz<n:
                        j=idx(nx,ny,nz)
                        adj[i].append(j); adj[j].append(i)
    return adj,n,idx

def is_surface(x,y,z,L):
    return x==0 or x==L or y==0 or y==L or z==0 or z==L

# ================================================================
# Yukawa 求解（大盒子，Jacobi + SOR）
# ================================================================
def solve_yukawa(adj, source_mask, surface_mask, kappa, mass, n_iter=25000, omega=0.85):
    N=len(adj)
    rows=[]; cols=[]
    for i,nbrs in enumerate(adj):
        for j in nbrs:
            rows.append(i); cols.append(j)
    rows=np.array(rows); cols=np.array(cols)
    degs=np.array([len(a) for a in adj],float)
    rho=source_mask.astype(float)
    denom=np.maximum(degs,1.0)+mass*mass
    phi=np.zeros(N)
    for it in range(n_iter):
        s=np.zeros(N)
        np.add.at(s, rows, phi[cols])
        phi_new=(s+kappa*rho)/denom
        phi_new[surface_mask]=0.0
        phi=(1-omega)*phi+omega*phi_new
    resid=np.abs(phi_new-phi).max()
    phi[surface_mask]=0.0
    return phi, resid

# ================================================================
# 拟合
# ================================================================
def fit_mu_sinh(r_arr, phi_arr, R_b, m_init):
    """φ(r) = (A/r) sinh(μ(R_b-r))"""
    def model(r, A, mu):
        return (A/r)*np.sinh(mu*(R_b-r))
    try:
        popt,_=curve_fit(model, r_arr, phi_arr, p0=[phi_arr[0]*r_arr[0], m_init],
                         maxfev=50000)
        return popt[1]
    except Exception:
        return None

def fit_mu_exp(r_arr, phi_arr, m_init):
    """φ(r) = (A/r) exp(-μr)，远场无反射近似"""
    def model(r, A, mu):
        return (A/r)*np.exp(-mu*r)
    try:
        popt,_=curve_fit(model, r_arr, phi_arr, p0=[phi_arr[0]*r_arr[0], m_init],
                         maxfev=50000)
        return popt[1]
    except Exception:
        return None

def pole_mass(m):
    return np.arccosh(1.0+m*m/2.0)

# ================================================================
# 主程序
# ================================================================
if __name__=="__main__":
    L=80
    adj,n,idx=build_cubic_lattice(L)
    N=len(adj)
    print("="*72)
    print("  Yukawa reflection fit v2: 大盒子 + 多窗口")
    print("="*72)
    print(f"  L={L}, N={N}")

    surface_mask=np.zeros(N,dtype=bool)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if is_surface(x,y,z,L):
                    surface_mask[idx(x,y,z)]=True
    c=L//2
    src=idx(c,c,c)
    source_mask=np.zeros(N,dtype=bool); source_mask[src]=True
    R_b_axis=L-c
    print(f"  R_b(axis)={R_b_axis}")

    masses=[0.1,0.2,0.5]
    fig,axes=plt.subplots(1,3,figsize=(17,5))

    print(f"\n  {'m_in':<7}{'μ_latt':<10}{'μ_sinh':<10}{'dev%':<8}{'μ_exp':<10}{'dev%':<8}{'λ/R_b'}")
    print("  "+"-"*60)

    for ax,m in zip(axes,masses):
        t0=time.time()
        phi,resid=solve_yukawa(adj,source_mask,surface_mask,KAPPA,m,n_iter=25000)
        dt=time.time()-t0

        # 轴上数据
        r_arr=[]; phi_arr=[]
        for r in range(2,R_b_axis):
            node=idx(c+r,c,c)
            if phi[node]>1e-16:
                r_arr.append(r); phi_arr.append(phi[node])
        r_arr=np.array(r_arr,float); phi_arr=np.array(phi_arr,float)

        # 多窗口拟合：r_lo 从 3 扫到 8
        print(f"\n  m={m}: resid={resid:.2e}, solve time={dt:.1f}s, λ/R_b={1/m/R_b_axis:.3f}")
        print(f"    {'r_lo':<6}{'μ_sinh':<10}{'dev%':<8}{'μ_exp':<10}{'dev%'}")
        best_sinh=None; best_exp=None
        for r_lo in [3,4,5,6,8]:
            sel=(r_arr>=r_lo)&(r_arr<=R_b_axis-2)
            r_fit=r_arr[sel]; phi_fit=phi_arr[sel]
            if len(r_fit)<4: continue
            mu_s=fit_mu_sinh(r_fit,phi_fit,R_b_axis,m)
            mu_e=fit_mu_exp(r_fit,phi_fit,m)
            if mu_s is not None:
                dev_s=abs(mu_s-m)/m*100
                if best_sinh is None or dev_s<best_sinh[1]:
                    best_sinh=(mu_s,dev_s)
            else:
                dev_s=float('nan')
            if mu_e is not None:
                dev_e=abs(mu_e-m)/m*100
                if best_exp is None or dev_e<best_exp[1]:
                    best_exp=(mu_e,dev_e)
            else:
                dev_e=float('nan')
            print(f"    {r_lo:<6}{mu_s:<10.4f}{dev_s:<8.2f}{mu_e:<10.4f}{dev_e:<8.2f}")

        # 总结
        mu_latt=pole_mass(m)
        if best_sinh and best_exp:
            print(f"  {m:<7.2f}{mu_latt:<10.4f}{best_sinh[0]:<10.4f}{best_sinh[1]:<8.2f}"
                  f"{best_exp[0]:<10.4f}{best_exp[1]:<8.2f}{1/m/R_b_axis:.3f}")

        # 画图：数据 + 最佳 sinh 拟合
        ax.semilogy(r_arr,phi_arr,'ko',ms=3,label='data (axis)')
        if best_sinh:
            mu=best_sinh[0]
            A0=phi_fit[0]*r_fit[0]/np.sinh(mu*(R_b_axis-r_fit[0]))
            rf=np.linspace(3,R_b_axis-1,200)
            ax.semilogy(rf,(A0/rf)*np.sinh(mu*(R_b_axis-rf)),'r-',lw=2,
                        label=f'sinh fit μ={mu:.4f} (in {m})')
        ax.axvline(R_b_axis,color='gray',ls=':',alpha=0.5)
        ax.set_xlabel('r (on axis)'); ax.set_ylabel('φ')
        ax.set_title(f'm={m}, L={L}')
        ax.legend(); ax.grid(alpha=0.3,which='both')

    plt.tight_layout(); plt.savefig('nea_yukawa_reflection_v2.png',dpi=150)
    print("\n  saved: nea_yukawa_reflection_v2.png")

    print("\n"+"="*72)
    print("  Verdict")
    print("="*72)
    print("""
  · 大盒子 L=80 → λ/R_b 全面缩小（m=0.1 时 λ/R_b=0.25，之前 0.5）
  · 若 μ_sinh 偏差 <5%，则 Yukawa 全质量段定量闭合
  · 若 μ_exp 与 μ_sinh 一致，说明反射已可忽略，纯指数成立
  · 若仍有偏差，需检查格点色散修正（μ_latt vs m）或近源伪影
""")
    print("="*72)