#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_3d_gravity_v5.py
决定性 β 测量：把源推离边界，只在干净远场窗口拟合

v4诊断：β=3.13 是边界污染。源离接地边界只有8格距，
        r>6 后场被 u=0 边界拽下，拟合窗口含尾部 → 假 β=3.13。
        u×r² 在 r=1..6 线性上涨(=1/r特征)，r>6 掉头(=边界拽)。
v5修正：
  · L=32（源到边界 R_b=16），远场窗口更宽
  · 只在 r ∈ [r_lo, R_b/2] 拟合（避开近场格子伪影 + 边界污染）
  · 同时打印全范围拟合做对照，透明展示范围依赖性
  · 用 u×r 平坦性做独立判据（1/r ⟺ u×r=常数）
"""
import numpy as np
import time
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

KAPPA = 1.0 - np.sqrt(3)/2   # Δ

# ============================================================
# 建纯立方格子（无动力学、无加边）→ 谱维度=3
# ============================================================
def build_cubic_graph(L):
    N = L**3
    G = [[] for _ in range(N)]
    def idx(x,y,z): return (x*L + y)*L + z
    for x in range(L):
        for y in range(L):
            for z in range(L):
                i = idx(x,y,z)
                for dx,dy,dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
                    nx,ny,nz = x+dx,y+dy,z+dz
                    if 0<=nx<L and 0<=ny<L and 0<=nz<L:
                        G[i].append(idx(nx,ny,nz))
    return G, idx

def is_surface(x,y,z,L):
    return x==0 or x==L-1 or y==0 or y==L-1 or z==0 or z==L-1

# ============================================================
# 向量化 Jacobi 求解 ∇²u=-κρ，表面 Dirichlet u=0
# ============================================================
def build_edge_arrays(G):
    rows=[]; cols=[]
    for i,nbrs in enumerate(G):
        for j in nbrs:
            rows.append(i); cols.append(j)
    return np.array(rows,dtype=int), np.array(cols,dtype=int)

def solve_field(G, source_mask, surface_mask, kappa, n_iter=9000, omega=0.6):
    N=len(G)
    rows,cols = build_edge_arrays(G)
    degs = np.array([len(G[i]) for i in range(N)],dtype=float)
    rho = source_mask.astype(float)
    u = np.zeros(N)
    for it in range(n_iter):
        u_sum = np.zeros(N)
        np.add.at(u_sum, rows, u[cols])
        u_new = (u_sum + kappa*rho)/np.maximum(degs,1.0)
        u_new[surface_mask]=0.0
        u = (1-omega)*u + omega*u_new
        if it%3000==0:
            resid = np.abs(u_new-u).max()
    u[surface_mask]=0.0
    resid = np.abs(u_new-u).max()
    return u, resid

# ============================================================
# BFS 距离
# ============================================================
def multi_source_bfs(G, sources):
    N=len(G)
    dist=np.full(N,-1,dtype=int)
    q=list(sources)
    for s in sources: dist[s]=0
    head=0
    while head<len(q):
        x=q[head]; head+=1
        for y in G[x]:
            if dist[y]<0:
                dist[y]=dist[x]+1; q.append(y)
    return dist

# ============================================================
# 幂律拟合（可指定窗口）
# ============================================================
def fit_beta(pr, pu, r_lo, r_hi):
    sel = [(pr[i],pu[i]) for i in range(len(pr))
           if r_lo<=pr[i]<=r_hi and pu[i]>1e-12]
    if len(sel)<3: return None,None,None,len(sel)
    rr=np.array([s[0] for s in sel],float)
    uu=np.array([s[1] for s in sel],float)
    c=np.polyfit(np.log(rr),np.log(uu),1)
    beta=-c[0]
    log_u=np.log(uu); pred=np.polyval(c,np.log(rr))
    r2=1-np.sum((log_u-pred)**2)/np.sum((log_u-log_u.mean())**2)
    return beta, np.exp(c[1]), r2, len(sel)

# ============================================================
# 主程序
# ============================================================
print("="*72)
print("  v5: 决定性 β 测量（源远离边界 + 窗口拟合）")
print("="*72)

L = 32            # ← 想要更宽窗口可改成 40 或 48（更慢）
N = L**3
c = L//2
print(f"  L={L}, N={N}, 源在中心 ({c},{c},{c})")

G, idx = build_cubic_graph(L)
surface_mask = np.zeros(N,dtype=bool)
for x in range(L):
    for y in range(L):
        for z in range(L):
            if is_surface(x,y,z,L):
                surface_mask[idx(x,y,z)]=True
n_surf=int(surface_mask.sum())

# 源：中心1个K4（4节点）
source_nodes=[idx(c,c,c),idx(c+1,c,c),idx(c,c+1,c),idx(c,c,c+1)]
source_mask=np.zeros(N,dtype=bool)
for s in source_nodes: source_mask[s]=True

# 源到边界的图距离 R_b
d_src_bfs = multi_source_bfs(G, source_nodes)
surf_nodes = np.where(surface_mask)[0]
R_b = int(d_src_bfs[surf_nodes].min())
print(f"  表面节点={n_surf}")
print(f"  源到边界距离 R_b = {R_b}")
print(f"  干净远场窗口取 r ∈ [3, {R_b//2}]")

# 求解
print(f"\n  场求解 (κ=Δ={KAPPA:.5f})")
t0=time.time()
u,resid = solve_field(G, source_mask, surface_mask, KAPPA, n_iter=9000, omega=0.6)
print(f"    耗时={time.time()-t0:.1f}s, 收敛残差={resid:.2e}")
print(f"    u_max={u.max():.6f} ({'✓弱场' if u.max()<0.5 else '✗强场'})")

# 径向剖面
max_r=int(d_src_bfs.max())
pr=[];pu=[];pn=[]
for r in range(0, min(max_r,R_b+6)+1):
    shell=np.where((d_src_bfs==r)&~source_mask&~surface_mask)[0]
    if len(shell)>0:
        pr.append(r); pu.append(u[shell].mean()); pn.append(len(shell))

print(f"\n  径向剖面（R_b={R_b}，边界在 r≈{R_b}）")
print(f"    {'r':<4}{'u':<14}{'u×r':<14}{'u×r²':<12}{'n'}")
for i in range(len(pr)):
    r=pr[i]
    ur = pu[i]*r if r>0 else 0
    ur2= pu[i]*r*r if r>0 else 0
    tag=""
    if r==R_b//2: tag=" ← 窗口上界"
    if r==R_b:    tag=" ← 边界"
    print(f"    {r:<4}{pu[i]:<14.7f}{ur:<14.7f}{ur2:<12.6f}{pn[i]}{tag}")

# 三种拟合，透明对比
print(f"\n  幂律拟合 u∝r^(-β) —— 展示范围依赖性")
b_full,_,r2_full,n_full = fit_beta(pr,pu,2,max_r)
b_win, A_win, r2_win, n_win = fit_beta(pr,pu,3,R_b//2)
b_win2,_,r2_win2,n_win2 = fit_beta(pr,pu,2,R_b//2)
print(f"    全范围 r∈[2,{max_r}]:        β={b_full:.3f}  (n={n_full}, R²={r2_full:.3f})  ← 含边界尾部，被污染")
print(f"    远场窗 r∈[3,{R_b//2}]:        β={b_win:.3f}  (n={n_win}, R²={r2_win:.3f})  ← 干净窗口")
print(f"    远场窗 r∈[2,{R_b//2}]:        β={b_win2:.3f}  (n={n_win2}, R²={r2_win2:.3f})")

# u×r 平坦性（1/r 的独立判据）
ur_win=[pu[i]*pr[i] for i in range(len(pr)) if 3<=pr[i]<=R_b//2 and pu[i]>1e-12]
if len(ur_win)>3:
    ur_arr=np.array(ur_win)
    disp=ur_arr.std()/ur_arr.mean()
    print(f"\n  u×r 平坦性（窗口内）: 均值={ur_arr.mean():.6f}, 离散度={disp:.3f}")
    print(f"    {'✓ u×r≈常数 → u∝1/r 确认' if disp<0.25 else '✗ u×r 不平坦'}")

print(f"\n  ★ 判据: 远场窗 β≈1 → 3D引力从交换格子涌现")
print(f"          全范围 β≫1 → 边界污染（预期，非物理）")

# ============================================================
# 绘图
# ============================================================
fig,axes=plt.subplots(2,2,figsize=(14,12))

# (a) log-log 剖面 + 窗口拟合线
ax=axes[0,0]
ax.loglog(pr,pu,'ko-',ms=4,lw=1.5,label='numerical')
ax.axvline(R_b,color='red',ls=':',alpha=0.6,label=f'boundary r={R_b}')
ax.axvline(R_b//2,color='green',ls=':',alpha=0.6,label=f'window edge r={R_b//2}')
if b_win is not None:
    rf=np.linspace(2,R_b//2,50)
    ax.loglog(rf, A_win/rf**b_win,'g--',lw=2,label=f'window fit β={b_win:.2f}')
ax.set_xlabel('graph distance r'); ax.set_ylabel('u')
ax.set_title('Field decay (log-log), green=clean window')
ax.legend(); ax.grid(alpha=0.3,which='both')

# (b) u×r 平坦性
ax=axes[0,1]
ur_plot=[(pr[i],pu[i]*pr[i]) for i in range(len(pr)) if pr[i]>0 and pu[i]>1e-12]
ax.plot([x[0] for x in ur_plot],[x[1] for x in ur_plot],'ro-',ms=4)
ax.axvline(R_b//2,color='green',ls=':',label=f'window edge')
ax.axvline(R_b,color='red',ls=':',label='boundary')
ax.set_xlabel('r'); ax.set_ylabel('u × r')
ax.set_title('1/r test: u×r flat ⟹ u∝1/r')
ax.legend(); ax.grid(alpha=0.3)

# (c) u×r² 线性性
ax=axes[1,0]
ur2_plot=[(pr[i],pu[i]*pr[i]**2) for i in range(len(pr)) if pr[i]>0 and pu[i]>1e-12]
ax.plot([x[0] for x in ur2_plot],[x[1] for x in ur2_plot],'bo-',ms=4)
ax.axvline(R_b//2,color='green',ls=':',label='window edge')
ax.axvline(R_b,color='red',ls=':',label='boundary')
ax.set_xlabel('r'); ax.set_ylabel('u × r²')
ax.set_title('u×r² linear before boundary, falls after')
ax.legend(); ax.grid(alpha=0.3)

# (d) 场切片
ax=axes[1,1]
slice_u=np.zeros((L,L))
for x in range(L):
    for y in range(L):
        slice_u[x,y]=u[idx(x,y,c)]
im=ax.imshow(slice_u.T,origin='lower',cmap='hot',extent=[0,L,0,L])
ax.plot(c+0.5,c+0.5,'c*',ms=15,label='source')
ax.set_title(f'u field (z={c} slice), L={L}')
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.legend()
plt.colorbar(im,ax=ax,shrink=0.8)

plt.tight_layout(); plt.savefig('nea_3d_gravity_v5.png',dpi=150)
print(f"\n  saved: nea_3d_gravity_v5.png")
print("="*72)