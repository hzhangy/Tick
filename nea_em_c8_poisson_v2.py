#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_em_c8_poisson_v2.py
v1 修正：画图 x/y 维度不匹配的崩溃 + 测量升级
  1. L: 16→32，源到边界 R_b: 8→16，远场窗口 r∈[2,8]（v1 只有 r∈[2,4] 三个点）
  2. E×r² 在更宽范围验证 1/r² 平台
  3. 画图 r_axis/E_axis 改为同一循环构建
注：泊松用 Jacobi 迭代求解，不做特征分解，L=32 (N=35937) 可行。
"""
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KAPPA = 1.0 - np.sqrt(3)/2

def build_cubic_lattice(L):
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
                        adj[i].append(j); adj[j].append(i)
    return adj, n

def is_surface(x,y,z,L):
    return x==0 or x==L or y==0 or y==L or z==0 or z==L

def solve_poisson(adj, source_mask, surface_mask, kappa, n_iter=10000, omega=0.6):
    N = len(adj)
    rows=[]; cols=[]
    for i,nbrs in enumerate(adj):
        for j in nbrs:
            rows.append(i); cols.append(j)
    rows=np.array(rows); cols=np.array(cols)
    degs=np.array([len(a) for a in adj],float)
    rho=source_mask.astype(float)
    phi=np.zeros(N)
    for it in range(n_iter):
        s=np.zeros(N)
        np.add.at(s, rows, phi[cols])
        phi_new=(s + kappa*rho)/np.maximum(degs,1.0)
        phi_new[surface_mask]=0.0
        phi=(1-omega)*phi + omega*phi_new
    resid=np.abs(phi_new-phi).max()
    phi[surface_mask]=0.0
    return phi, resid

def bfs_distances(adj, sources):
    N=len(adj)
    dist=np.full(N,-1,dtype=int)
    q=list(sources)
    for s in sources: dist[s]=0
    head=0
    while head<len(q):
        x=q[head]; head+=1
        for y in adj[x]:
            if dist[y]<0:
                dist[y]=dist[x]+1; q.append(y)
    return dist

def fit_beta(pr,pu,r_lo,r_hi):
    sel=[(pr[i],pu[i]) for i in range(len(pr)) if r_lo<=pr[i]<=r_hi and pu[i]>1e-12]
    if len(sel)<3: return None,None,len(sel)
    rr=np.array([s[0] for s in sel],float)
    uu=np.array([s[1] for s in sel],float)
    c=np.polyfit(np.log(rr),np.log(uu),1)
    beta=-c[0]
    pred=np.polyval(c,np.log(rr))
    r2=1-np.sum((np.log(uu)-pred)**2)/np.sum((np.log(uu)-np.log(uu).mean())**2)
    return beta,r2,len(sel)

if __name__=="__main__":
    L=32
    adj,n=build_cubic_lattice(L)
    N=len(adj)
    print("="*72)
    print("  C8 密铺上的电磁势 v2 (L=32, 宽远场窗口)")
    print("="*72)
    print(f"  L={L}, N={N}")

    surface_mask=np.zeros(N,dtype=bool)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if is_surface(x,y,z,L):
                    surface_mask[(x*n+y)*n+z]=True
    c=L//2
    src=(c*n+c)*n+c
    source_mask=np.zeros(N,dtype=bool); source_mask[src]=True

    t0=time.time()
    phi,resid=solve_poisson(adj,source_mask,surface_mask,KAPPA,n_iter=10000)
    print(f"  Jacobi: {time.time()-t0:.1f}s, 残差={resid:.2e}, phi_max={phi.max():.5f} (弱场)")

    d2s=bfs_distances(adj,[src])
    R_b=int(d2s[surface_mask].min())
    r_hi=R_b//2
    print(f"  R_b={R_b}, 远场窗口 r∈[2,{r_hi}]")

    internal_free=(~surface_mask)&(~source_mask)
    pr,pu=[],[]
    for r in range(1,min(int(d2s.max()),25)+1):
        shell=np.where((d2s==r)&internal_free)[0]
        if len(shell)>0:
            pr.append(r); pu.append(phi[shell].mean())

    print(f"\n  {'r':<5}{'phi':<14}{'phi×r':<12}")
    for i in range(len(pr)):
        print(f"  {pr[i]:<5}{pu[i]:<14.7f}{pu[i]*pr[i]:<12.7f}")

    beta,r2,nfit=fit_beta(pr,pu,2,r_hi)
    print(f"\n  远场拟合 r∈[2,{r_hi}]: β={beta:.3f} (R²={r2:.4f}, n={nfit})")
    if beta is not None and abs(beta-1)<0.2:
        print("  ✓ β≈1：库仑 1/r 确认")
    else:
        print("  ? β 偏离 1，需检查")

    # 轴线 E×r²
    axis_pts=[]
    for x in range(1,L):
        if x==c: continue
        i=(x*n+c)*n+c
        r=abs(x-c)
        El=phi[((x-1)*n+c)*n+c]; Er=phi[((x+1)*n+c)*n+c]
        E=-(Er-El)/2.0
        axis_pts.append((r,phi[i],E))

    print(f"\n  轴线 E×r²（1/r² 平台检验, r≤{r_hi+2}）:")
    ur=[]
    for r,p,E in axis_pts:
        if 2<=r<=r_hi+2:
            print(f"    r={r:<3} E×r²={abs(E)*r*r:.6f}")
            ur.append(abs(E)*r*r)
    if ur:
        ur=np.array(ur)
        print(f"    均值={ur.mean():.6f}, 离散度={ur.std()/ur.mean():.3f}")

    # 画图（修正：同一循环构建 r_axis/E_axis）
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    ax=axes[0]
    ax.loglog(pr,pu,'ro-',ms=4,lw=2,label='numerical φ')
    if beta is not None:
        rr=np.array([pr[i] for i in range(len(pr)) if 2<=pr[i]<=r_hi],float)
        uu=np.array([pu[i] for i in range(len(pr)) if 2<=pr[i]<=r_hi],float)
        cc=np.polyfit(np.log(rr),np.log(uu),1)
        rf=np.linspace(2,r_hi,50)
        ax.loglog(rf,np.exp(cc[1])/rf**beta,'b--',lw=2,label=f'fit β={beta:.2f}')
    ax.axvline(r_hi,color='green',ls=':',alpha=0.6,label=f'window r={r_hi}')
    ax.set_xlabel('r'); ax.set_ylabel('φ')
    ax.set_title('Coulomb potential on C8 tiling')
    ax.legend(); ax.grid(alpha=0.3,which='both')

    ax=axes[1]
    pts=[(r,abs(E)) for r,p,E in axis_pts if 2<=r<=r_hi+2]
    if pts:
        ax.loglog([q[0] for q in pts],[q[1] for q in pts],'bo-',ms=4,label='|E|')
        r0,E0=pts[0]
        rf=np.linspace(2,r_hi+2,50)
        ax.loglog(rf,E0*r0*r0/rf**2,'g--',lw=2,label='1/r²')
    ax.set_xlabel('r'); ax.set_ylabel('|E|')
    ax.set_title('Electric field on axis')
    ax.legend(); ax.grid(alpha=0.3,which='both')

    plt.tight_layout(); plt.savefig('nea_em_c8_poisson_v2.png',dpi=150)
    print("  saved: nea_em_c8_poisson_v2.png")