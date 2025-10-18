#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np, math
from numpy.linalg import eig, norm, det
from pymatgen.symmetry.groups import PointGroup
from pymatgen.core import Lattice
import random
from .eqvpg2label import get_std_pg
from .small_func import get_rotation,orthonormal_basis_from_vector
TOL = 1e-4

def plane_basis(u):
    u = np.asarray(u, float); u /= norm(u)
    a = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(a, u)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    e1 = a - np.dot(a, u) * u
    e1 /= norm(e1)
    e2 = np.cross(u, e1)
    return e1, e2, u

def restrict_to_plane(M, B):

    Mp = B.T @ M @ B
    return Mp[:2, :2]

def solve_Q2_ordered(R2_list, S2_list, atol=TOL):

    I2 = np.eye(2)
    blocks = [np.kron(R.T, I2) - np.kron(I2, S) for R, S in zip(R2_list, S2_list)]
    M = np.vstack(blocks)
    try:
        _, _, Vt = np.linalg.svd(M)
    except np.linalg.LinAlgError:
        return None
    q = Vt[-1, :]
    Qhat = q.reshape(2, 2, order='F')
    U, _, Vt = np.linalg.svd(Qhat)
    Q2 = U @ Vt            
    
    for Qcand in (Q2, -Q2):
        ok = all(np.allclose(Qcand @ R @ Qcand.T, S, atol=atol)
                 for R, S in zip(R2_list, S2_list))
        if ok:
            return Qcand
    return None



def _match_ok(P, inp_R, std_R, ordered: bool):
    def close(A,B): return np.allclose(A, B, atol=TOL)

    if ordered:
        for i in range(len(inp_R)):
            PRP = P @ inp_R[i] @ P.T
            S = std_R[i]
            if close(PRP, S):                
                continue
            return False
        return True
    else:
        for R in inp_R:
            PRP = P @ R @ P.T
            ok = any(close(PRP, S)
                     for S in std_R)
            if not ok:
                return False
        return True

N_BY_SYM = {
    "1":1, "m":1, "-1":1,
    "2":2, "2/m":2, "mm2":2, "222":2, "mmm":2,
    "4":4, "-4":4, "4/m":4, "4mm":4, "-42m":4, "422":4, "4/mmm":4,
    "3":3, "-3":3, "3m":3, "-3m":3, "32":3,
    "6":6, "-6":6, "6/m":6, "6mm":6, "-6m2":6, "622":6, "6/mmm":6,
    "23":3, "m-3":3, "432":4, "-43m":3, "m-3m":4,
}
PG_LIST = list(N_BY_SYM.keys())
def axis_from_rot(R: np.ndarray) -> np.ndarray:
    v = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]], float)
    if norm(v) < 1e-8:
        vals, vecs = eig(R)
        idx = next((i for i,val in enumerate(vals) if abs(val-1) < 1e-6), None)
        if idx is None: return np.zeros(3)
        v = np.real(vecs[:, idx])
    n = norm(v)
    if n < 1e-12: return np.zeros(3)
    v = v/n
    idx_nz = np.nonzero(np.abs(v) > 1e-6)[0]
    if idx_nz.size>0 and v[idx_nz[0]] < 0: v = -v
    return v

def _solve_P_ordered(inp_R, std_R, atol=TOL):
    pairs = [(R, S) for (R, S) in zip(inp_R, std_R)
             if not np.allclose(R, np.eye(3), atol=atol)]
    if len(pairs) < 2:
        pairs = list(zip(inp_R, std_R))

    blocks = []
    I = np.eye(3)
    for R, S in pairs:
        blocks.append(np.kron(R.T, I) - np.kron(I, S))
    M = np.vstack(blocks)

    try:
        _, _, Vt = np.linalg.svd(M)
    except np.linalg.LinAlgError:
        return None
    v = Vt[-1, :]
    Phat = v.reshape(3, 3, order='F')

    U, _, Vt = np.linalg.svd(Phat)
    P0 = U @ Vt
    
    cands = [P0, -P0]
    for P in cands:
        if all(np.allclose(P @ inp_R[i] @ P.T, std_R[i], atol=atol) for i in range(len(inp_R))):
            return P

    return None

def collect_axes(rots):
    for sgn in (+1, -1):
        axes=[]
        for R in rots:
            if (sgn>0 and det(R)<0) or (sgn<0 and det(R)>0): continue
            if np.allclose(R, np.eye(3), atol=TOL): continue
            v = axis_from_rot(R)
            if norm(v) < 1e-8: continue
            if not any(np.allclose(v,u,atol=1e-6) or np.allclose(v,-u,atol=1e-6) for u in axes):
                axes.append(v)
        if axes: return axes
    return []

def rot_about_axis(axis, angle):
    axis = np.array(axis,float); axis /= norm(axis)
    x,y,z = axis; c,s = math.cos(angle), math.sin(angle); t = 1-c
    return np.array([[t*x*x+c, t*x*y-s*z, t*x*z+s*y],
                     [t*x*y+s*z, t*y*y+c, t*y*z-s*x],
                     [t*x*z-s*y, t*y*z+s*x, t*z*z+c]], float)

def rot_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a/norm(a); b = b/norm(b)
    v = np.cross(a,b); s = norm(v); c = float(np.dot(a,b))
    if s < 1e-10:
        if c > 0: return np.eye(3)  
        tmp = np.array([1.0,0.0,0.0])
        if abs(np.dot(tmp,a)) > 0.9: tmp = np.array([0.0,1.0,0.0])
        axis = np.cross(a,tmp); 
        if norm(axis) < 1e-12: axis = np.cross(a, np.array([0,0,1.0]))
        axis /= norm(axis)
        return rot_about_axis(axis, math.pi)
    K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]], float)
    return np.eye(3) + K + K@K*((1-c)/(s*s))

def mirror_through_axis_plane(axis):
    axis = np.asarray(axis, float); axis /= norm(axis)
    v = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(v, axis)) > 0.9:
        v = np.array([0.0, 1.0, 0.0])
    n = np.cross(axis, v); n /= norm(n)

    return np.eye(3) - 2.0 * np.outer(n, n)

def _axis_order(axis, R_list):
    for R in R_list:
        v = axis_from_rot(R)
        if norm(v)<1e-6 or det(R)<=0: continue
        if np.allclose(abs(np.dot(v,axis)), 1.0, atol=1e-3):
            th = math.acos(max(min((np.trace(R)-1)/2,1),-1))
            if th < 1e-6: continue
            n = int(round(2*math.pi/th))
            if n>=1: return n
    return None

def _find_similarity_matrix(sym, inp_R, std_R, ordered=False):
    if ordered:
        P_lin = _solve_P_ordered(inp_R, std_R, atol=TOL)
        if P_lin is not None:
            return P_lin
    axes_in  = collect_axes(inp_R)
    axes_std = collect_axes(std_R)
    if not axes_std: return None

    if len(axes_std) >= 2:
        in_ord  = [ _axis_order(a, inp_R)  for a in axes_in ]
        std_ord = [ _axis_order(a, std_R)  for a in axes_std ]
        for i, a_in in enumerate(axes_in):
            oi = in_ord[i]
            if oi is None: continue
            for j, a_std in enumerate(axes_std):
                oj = std_ord[j]
                if oj is None or oi!=oj: continue
                R0 = rot_between(a_in, a_std)
                for k, b_in in enumerate(axes_in):
                    if k==i or abs(np.dot(b_in,a_in))>1-1e-3: continue
                    ok = in_ord[k]
                    if ok is None: continue
                    for l, b_std in enumerate(axes_std):
                        if l==j or abs(np.dot(b_std,a_std))>1-1e-3: continue
                        ol = std_ord[l]
                        if ol is None or ol!=ok: continue
                        w1 = R0@b_in
                        w1p = w1 - np.dot(w1,a_std)*a_std
                        v2p = b_std - np.dot(b_std,a_std)*a_std
                        if norm(w1p)<1e-8 or norm(v2p)<1e-8: continue
                        w1u, v2u = w1p/norm(w1p), v2p/norm(v2p)
                        cosang = max(min(np.dot(w1u,v2u),1),-1)
                        ang = math.acos(cosang)
                        if np.dot(np.cross(w1u,v2u), a_std) < 0: ang = -ang
                        P = rot_about_axis(a_std, ang) @ R0
                        # if det(P) < 0: continue
                        if _match_ok(P, inp_R, std_R, ordered):
                            return P
                        M = mirror_through_axis_plane(a_std)
                        P2 = M @ P
                        if _match_ok(P2, inp_R, std_R, ordered):
                            return P2
        return None

    axis_std = axes_std[0]
    axis_in  = axes_in[0]
    R_align_axis = rot_between(axis_in, axis_std)

    def R_axis(angle):
        return rot_about_axis(axis_std, angle)

    # Pre-align using the mirror normal (key step that previously worked)
    if any(det(R) < 0 for R in std_R):
        Rin_m  = next((R for R in inp_R if det(R) < 0), None)
        Rstd_m = next((R for R in std_R if det(R) < 0), None)
        if Rin_m is not None and Rstd_m is not None:
            Rin_m_al = R_align_axis @ Rin_m @ R_align_axis.T

            def mirror_normal(M):
                vals, vecs = eig(M + np.eye(3))
                for i, val in enumerate(vals):
                    if abs(val) < 1e-6:
                        n = np.real(vecs[:, i])
                        if norm(n) > 1e-8:
                            return n / norm(n)
                return None

            n_in  = mirror_normal(Rin_m_al)
            n_std = mirror_normal(Rstd_m)
            if n_in is not None and n_std is not None:
                n_in_p  = n_in  - np.dot(n_in , axis_std) * axis_std
                n_std_p = n_std - np.dot(n_std, axis_std) * axis_std
                if norm(n_in_p) > 1e-8 and norm(n_std_p) > 1e-8:
                    a = n_in_p / norm(n_in_p); b = n_std_p / norm(n_std_p)
                    cosang = max(min(np.dot(a, b), 1), -1)
                    ang = math.acos(cosang)
                    if np.dot(np.cross(a, b), axis_std) < 0:
                        ang = -ang
                    Pm = R_axis(ang) @ R_align_axis
                    if _match_ok(Pm, inp_R, std_R, ordered):
                        return Pm

    M0 = mirror_through_axis_plane(axis_std)
    n  = max(1, N_BY_SYM.get(sym, 1))

    for k in range(n):
        P = R_axis(2 * math.pi * k / n) @ R_align_axis
        if _match_ok(P, inp_R, std_R, ordered):
            return P

    for k in range(2 * n):
        P = R_axis(math.pi * k / n) @ R_align_axis
        if _match_ok(P, inp_R, std_R, ordered):
            return P

    for k in range(n):
        Prot = R_axis(2 * math.pi * k / n) @ R_align_axis
        P = M0 @ Prot
        if _match_ok(P, inp_R, std_R, ordered):
            return P

    for k in range(2 * n):
        Prot = R_axis(math.pi * k / n) @ R_align_axis
        P = M0 @ Prot
        if _match_ok(P, inp_R, std_R, ordered):
            return P

    for k_theta in range(2 * n):               
        Prot = R_axis(math.pi * k_theta / n) @ R_align_axis
        for k_phi in range(8 * n):              
            phi  = math.pi * k_phi / (8.0 * n)
            Mphi = R_axis(phi) @ M0 @ R_axis(-phi)
            P    = Mphi @ Prot
            if _match_ok(P, inp_R, std_R, ordered):
                return P

    
    inp_aligned = [R_align_axis @ R @ R_align_axis.T for R in inp_R]
    e1, e2, u = plane_basis(axis_std)
    B = np.column_stack([e1, e2, u])


    R2_list = [restrict_to_plane(M, B) for M in inp_aligned]
    S2_list = [restrict_to_plane(M, B) for M in std_R]

    Q2 = solve_Q2_ordered(R2_list, S2_list, atol=TOL)
    if Q2 is not None:

        Q3 = np.eye(3)
        Q3[:2, :2] = Q2
        P = B @ Q3 @ B.T @ R_align_axis
        if _match_ok(P, inp_R, std_R, ordered):
            return P
    return None

def identify_pointgroup(rots, ordered=False):
    for sym in PG_LIST:
        std_R = pg_ops_cartesian(sym)
        if len(rots) != len(std_R):
            continue
        P = _find_similarity_matrix(sym, rots, std_R, ordered=ordered)
        if P is not None:
            return sym, P
    return None, None

def rand_SO3():
    q=np.random.randn(4); q/=norm(q); w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]], float)

def self_test(sym, rounds=10,order=False):
    random.seed()
    std_R = pg_ops_cartesian(sym)
    for _ in range(rounds):
        Q = rand_SO3()
        rots = [Q@R@Q.T for R in std_R]
        if not order:
            random.shuffle(rots)
        
        s_out,P = identify_pointgroup(rots, ordered=order)
        if s_out != sym or P is None:
            print(f"**Test failed:** point group {sym} not identified, got {s_out}")
            return
        if order:
            err = max(
                np.linalg.norm(P @ rots[i] @ P.T - std_R[i], ord='fro')
                for i in range(len(std_R))
            )
        else:
            err = max(min(norm(P@R@P.T - S) for S in std_R) for R in rots)
        if err > 1e-5:
            print(f"**Test failed:** {sym} residual {err}")
            return
    print(f"Point group {sym} passed the test in {rounds} rounds.")

def pg_ops_cartesian(sym: str, *, a=1.0, b=None, c=None, bravais=None):
    # b = random.random()
    # c = random.random()
    if bravais is None:
        if sym in {"432","m-3m","-43m","m-3"}:
            bravais = "cubic"
        elif sym in {"4","-4","4/m","4mm","-42m","422","4/mmm"}:
            bravais = "tetragonal"
        elif sym in {"3","-3","32","3m","-3m"}:
            bravais = "hexagonal"  
        elif sym in {"6","-6","6/m","622","6mm","-6m2","6/mmm"}:
            bravais = "hexagonal"
        elif sym in {"222","mm2","mmm"}:
            bravais = "orthorhombic"
        elif sym in {"2","m","2/m"}:
            bravais = "monoclinic"
        else:
            bravais = "cubic"

    if bravais == "cubic":
        lat = Lattice.cubic(a)
    elif bravais == "tetragonal":
        lat = Lattice.tetragonal(a, c or 1.5*a)
    elif bravais == "hexagonal":
        lat = Lattice.hexagonal(a, c or 1.6*a)
    elif bravais == "orthorhombic":
        lat = Lattice.orthorhombic(a, b or 1.2*a, c or 1.4*a)
    elif bravais == "monoclinic":
        lat = Lattice.monoclinic(a, b or 1.2*a, c or 1.4*a, 110)
    else:
        lat = Lattice.from_parameters(a, b or 1.2*a, c or 1.4*a, 95, 102, 107)

    A = lat.matrix.T
    Ainv = np.linalg.inv(A)

    ops = PointGroup(sym).symmetry_ops
    R_frac = [np.rint(op.rotation_matrix).astype(int) for op in ops]
    R_cart = [A @ R @ Ainv for R in R_frac]

    return R_cart


def spin_axis(spin_rot_list,symbol):
    symbol_sf = get_std_pg(symbol)[1]
    if symbol_sf in ['C1','Ci']:
        axis = [np.array([0,0,1]),np.array([1,0,0])]
        return axis
    element = []
    for i in range(len(spin_rot_list)):
        element.append(get_rotation(spin_rot_list[i]))
        
    axis = []
    
    if symbol_sf[0] in ['C','S']:   
          
        if 'v' in symbol_sf or 'h' in symbol_sf:
            order = eval(symbol_sf[1:-1])
        elif 's' in symbol_sf:
            order = 2
        elif 'S' in symbol_sf:
            order = eval(symbol_sf[1:])/2
        else:
            order = eval(symbol_sf[1:])
        
        for i in element:
            if i[2] > 0:
                if abs(360/order - i[0]) < 1e-2:
                    axis.append(i[1])
                    break

        if 's' in symbol_sf:
            for i in element:
                if i[2] < 0:
                    if abs(180 - i[0]) < 1e-2:
                        axis.append(i[1])
                        break
                    
        if 's' in symbol_sf or 'v' in symbol_sf:
            for i in element:
                if i[2] < 0:
                    if abs(180 - i[0]) < 1e-2:
                        axis.append(i[1])
                        break
        else:
            axis.append(orthonormal_basis_from_vector(axis[0]))
        
    elif symbol_sf[0] == 'D':
        

        if 'd' in symbol_sf or 'h' in symbol_sf:
            order = eval(symbol_sf[1:-1])
        else:
            order = eval(symbol_sf[1:])
            
        for i in element:
            if symbol_sf == 'D2d':
                if i[2] < 0:
                    if abs(90 - i[0]) < 1e-2:
                        axis.append(i[1])
                        break
            else:
                if i[2] > 0:
                    if abs(360/order - i[0]) < 1e-2:
                        axis.append(i[1])
                        break
        
        for i in element:
            if abs(180 - i[0]) < 1e-2 and abs(np.dot(i[1],axis[0])) < 1e-2:
                axis.append(i[1])
                break
            
    elif symbol_sf[0] in ['T','I']:
        for i in element:
            if i[2] > 0:
                if abs(180 - i[0]) < 1e-2:
                    axis.append(i[1])
                    break
            
        for i in element:
            if i[2] > 0:
                if abs(180 - i[0]) < 1e-2 and abs(np.dot(i[1],axis[0])) < 1e-2:
                    axis.append(i[1])
                    break
    
    elif symbol_sf[0] == 'O':
        for i in element:
            if i[2] > 0:
                if abs(90 - i[0]) < 1e-2:
                    axis.append(i[1])
                    break
            
        for i in element:
            if i[2] > 0:
                if abs(90 - i[0]) < 1e-2 and abs(np.dot(i[1],axis[0])) < 1e-2:
                    axis.append(i[1])
                    break
    else:
        print('Unknown point group:',symbol_sf)
    axis[0] = axis[0] / norm(axis[0])
    axis[1] = axis[1] / norm(axis[1])
    return axis




if __name__ == "__main__":

    self_test('4mm',rounds=10,order=True)
