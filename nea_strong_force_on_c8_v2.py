#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nea_strong_force_on_c8_v2.py
C8 密铺上两个 K4 缺陷的强力势能（修正版）

修正与诚实标注：
  1. K4 成为真正的团：添加面对角线（正四面体6条边）。
     注：裸立方格中 K4 的4个交替顶点两两相距√2，互不相邻，原本不是团。
  2. 多方向测量：x向 / 面对角 / 体对角。
  3. 横轴用真实图距离（BFS最短路径），三方向可比。
  4. 诚实标注：力程在超出最近邻后归零，是格子只有最近邻边的拓扑截断，
     不是物理衰减。真正的 Yukawa 短程力需要质量标度 E0（仍是开放问题）。

势能模型（与 NEA 一致）：
  V = -Δ × cross_edges  +  V_core × shared_nodes²
  交叉边(交换) → 吸引；共享节点(带宽冲突) → 硬核排斥。
"""
import numpy as np
from collections import deque
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DELTA   = 1.0 - np.sqrt(3)/2
V0_CORE = 5.0 * DELTA          # 硬核排斥强度

# ---------------- 立方格 ----------------
class Lattice:
    def __init__(self, L):
        self.L = L; self.n = L+1; self.N = self.n**3
        self.adj = [set() for _ in range(self.N)]
        for x in range(self.n):
            for y in range(self.n):
                for z in range(self.n):
                    i = self.idx(x,y,z)
                    for dx,dy,dz in [(1,0,0),(0,1,0),(0,0,1)]:
                        nx,ny,nz = x+dx,y+dy,z+dz
                        if nx<self.n and ny<self.n and nz<self.n:
                            j = self.idx(nx,ny,nz)
                            self.adj[i].add(j); self.adj[j].add(i)
    def idx(self,x,y,z):
        return (x*self.n + y)*self.n + z

# ---------------- K4 = 立方体的交替顶点 ----------------
def k4_verts(lat, base, parity):
    x0,y0,z0 = base; v=[]
    for dx in (0,1):
        for dy in (0,1):
            for dz in (0,1):
                if (dx+dy+dz)%2 == parity:
                    v.append(lat.idx(x0+dx,y0+dy,z0+dz))
    return set(v)

def make_clique(lat, verts):
    """添加面对角线，使4个交替顶点成为正四面体团。返回新增边数。"""
    vl=list(verts); added=0
    for i in range(4):
        for j in range(i+1,4):
            a,b = vl[i],vl[j]
            if b not in lat.adj[a]:
                lat.adj[a].add(b); lat.adj[b].add(a); added+=1
    return added

def cross_edges(lat, A, B):
    c=0
    for a in A:
        for b in lat.adj[a]:
            if b in B: c+=1
    return c

def min_graph_dist(lat, A, B):
    """A 到 B 的最小图距离（BFS）。"""
    dist={a:0 for a in A}
    q=deque(A)
    while q:
        u=q.popleft()
        for w in lat.adj[u]:
            if w not in dist:
                dist[w]=dist[u]+1; q.append(w)
    return min(dist.get(b,10**9) for b in B)

def potential(cross, shared):
    return -DELTA*cross + V0_CORE*shared**2

def measure(base1, base2, parity1=0, parity2=0, L=12):
    lat=Lattice(L)
    A=k4_verts(lat,base1,parity1)
    B=k4_verts(lat,base2,parity2)
    make_clique(lat,A); make_clique(lat,B)
    cross  = cross_edges(lat,A,B)
    shared = len(A & B)
    d      = min_graph_dist(lat,A,B)
    V      = potential(cross,shared)
    return d, cross, shared, V

# ---------------- 主程序 ----------------
if __name__=="__main__":
    print("="*72)
    print("  C8 密铺上两个 K4 缺陷的强力势能（修正版）")
    print("="*72)
    print(f"  DELTA={DELTA:.6f},  V_core={V0_CORE:.4f}")

    # K4 团一致性自检
    lat0=Lattice(6); A0=k4_verts(lat0,(2,2,2),0)
    bare = sum(1 for a in A0 for b in lat0.adj[a] if b in A0)//2
    add  = make_clique(lat0,A0)
    full = sum(1 for a in A0 for b in lat0.adj[a] if b in A0)//2
    print(f"\n  [一致性自检] 裸格中 K4 内部边数={bare}（应0），"
          f"加面对角线后={full}（应6，真团）")

    base1=(5,5,5)
    directions={
        "x-axis":     lambda k:(k,0,0),
        "face-diag":  lambda k:(k,k,0),
        "body-diag":  lambda k:(k,k,k),
    }

    results={}
    for name,fn in directions.items():
        print(f"\n  {name}：")
        print(f"    {'k':<4}{'图距d':<7}{'cross':<7}{'shared':<8}{'V':<12}")
        rows=[]
        for k in range(0,6):
            off=fn(k)
            base2=(base1[0]+off[0], base1[1]+off[1], base1[2]+off[2])
            if max(base2)>11: break
            d,cross,shared,V=measure(base1,base2)
            rows.append((d,cross,shared,V))
            print(f"    {k:<4}{d:<7}{cross:<7}{shared:<8}{V:<12.4f}")
        results[name]=rows

    # 绘图：V vs 真实图距离
    fig,ax=plt.subplots(figsize=(8,6))
    colors={"x-axis":"steelblue","face-diag":"darkorange","body-diag":"seagreen"}
    for name,rows in results.items():
        ds=[r[0] for r in rows]; Vs=[r[3] for r in rows]
        ax.plot(ds,Vs,'o-',lw=2,color=colors[name],label=name)
    ax.axhline(0,color='gray',ls=':',lw=1)
    ax.set_xlabel('Graph distance d (BFS)')
    ax.set_ylabel('V (units of Δ)')
    ax.set_title('Strong force potential: K4 defects vs graph distance')
    ax.legend(); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('nea_strong_force_on_c8_v2.png',dpi=150)
    print("\n  saved: nea_strong_force_on_c8_v2.png")

    print("\n" + "="*72)
    print("  诚实判读")
    print("="*72)
    print("""
  · 定性形状正确：重叠(d=0)硬核排斥，最近邻吸引，超出即零。
  · 但 d≥2 处 V=0 是【格子只有最近邻边】的拓扑截断，
    不是物理衰减。任何最近邻格子上的缺陷都长这样。
  · 真正的短程 Yukawa 力需要质量标度 E0（力程=ℏc/E0），
    而 E0 正是 NEA 尚未闭合的那道墙。
  · 面对角线修正让 K4 成为真正的团（NEA一致性），
    但不改变两个不同 K4 之间的交叉边数。
""")
    print("="*72)