import numpy as np
from numpy.linalg import norm, inv
import os
from .small_func import *
from .SG_utils import identify_SG_lattice


script_dir = os.path.dirname(os.path.abspath(__file__))

data1_path = os.path.join(script_dir, 'ssg_data/type1sg_ops.npy')

SG_data = np.load(data1_path, allow_pickle=True)


def get_sg_iso_tau(rot_list, tau_list, iso_R, tol=1e-6):
    '''
    For symorphic SG and W=(A, 0), check if A @ G @ A^-1 = G, where the ordering of ops may change, 
    but the matrix form of R cannot change.
    For non-symorphic SG, check first if the point group part satisfies A @ P_G @ A^-1 = P_G,
    and then solve if there exist t, s.t W=(A, t) satisfies W @ G @ W^-1 = G
    
    return: is_iso: True/False, iso_tau, and 
            iso_mapping=[j1, j2, ..., jn]: f(G) = W G W^-1, map nth op to jth op (for sorted G ops)
            (rot_trans_list[i] = rot_orig_list[iso_mapping[i]])

    1. compare G and G': first sort ops using (angle, axis, det), then compare sorted op directly
    2. input rot_list and tau_list are sorted
    3. input ops are in G prim latt, iso_R also in G prim latt
    '''
    is_iso, iso_tau, iso_mapping = False, None, []
    nop = len(rot_list)
    rot_trans_list = [round_mat(iso_R @ R @ inv(iso_R)) for R in rot_list]
    tau_trans_list = [iso_R @ t for t in tau_list]

    sorted_rot_trans_list, sorted_tau_trans_list, sort_order = sort_rot(rot_trans_list, tau_trans_list)
    if all([np.allclose(R1, R2) for R1, R2 in zip(rot_list, sorted_rot_trans_list)]):
        iso_mapping = [find_mat(R_trans, rot_list) for R_trans in rot_trans_list]

        # if transformed tau are identical, then iso_tau is zero
        if all([norm(t1 - t2) < tol for t1, t2 in zip(tau_list, sorted_tau_trans_list)]):
            is_iso, iso_tau = True, [np.zeros(3)]
        else:
            # compute iso_tau using sorted ops, (iso_R, t) @ (R1, t1) = (R2, t2) @ (iso_R, t)
            # ==> (R2 - 1) t = iso_R @ t1 - t2 % 1, transform into M @ t = w % [1]
            M = np.zeros((3 * nop, 3))
            w = np.zeros(3 * nop)
            for iop, R1, t1, iop2, R2 in zip(range(nop), rot_list, tau_list, iso_mapping, rot_trans_list):
                t2 = tau_list[iop2]
                M[iop * 3: (iop + 1) * 3, :] = R2 - np.eye(3)
                w[iop * 3: (iop + 1) * 3] = iso_R @ t1 - t2

            # solve M @ t = w % [1]
            trial_tlist = solve_linear_inhomogeneous(M, w, mod_num=[1,1,1], allow_infinite_sol=True)
            if len(trial_tlist) > 0:
                is_iso = True
                iso_tau = trial_tlist
                
                # check result
                # tau_trans_list = [iso_tau + iso_R @ t - iso_R @ R @ inv(iso_R) @ iso_tau
                #                     for R, t in zip(rot_list, tau_list)]
                # assert all([np.allclose(rot_trans_list[i], rot_list[iso_mapping[i]]) for i in range(nop)]), (
                #             rot_list, rot_trans_list, iso_mapping)
                # assert all([norm(np.mod(round_vec(tau_trans_list[i] - tau_list[iso_mapping[i]]),
                #                         [1,1,1])) < tol for i in range(nop)]), (tau_trans_list, tau_list)
                
    return is_iso, iso_tau, iso_mapping


def find_sg_iso_transform(gid, tol=1e-6, test=0):
    '''
    find (A, t) G (A, t)^-1 = G (the ordering of ops may change, need to find group isomorphism)
    Method:
        1. For a given crystal system, the isomorphism of rotation part is easy to obtain (given manually) 
        2. For each sg, first check the rotation isomorphism, then solve translation part.
    Output:
        a dictionary {'iso_rot': [R, ...], 'iso_tau': [t, ...]}
    '''
    gid_data = SG_data[gid - 1]
    nop = len(gid_data['rotP']) // 2
    sorted_rot_list, sorted_tau_list, sort_order = sort_rot(gid_data['rotP'][0: nop], gid_data['tauP'][0: nop])
    prim_basis = identify_SG_lattice(gid)[1] # each col is a prim basis vector

    if gid in [1, 2]:
        # any SL(3, Z) is an isomorphism
        iso_sg = 207 
    elif gid in range(3, 16):
        # any SL(2, Z) along x, z axis is an isomorphism, and C2x, C2z
        # here use 207 to include C2x, C2z (extra ops will be abandoned)
        iso_sg = 207 
    elif gid in range(16, 75):
        # PG: 432 (sg 207)
        iso_sg = 207
    elif gid in range(75, 143):
        # PG: 422 (sg 89)
        iso_sg = 89
    elif gid in range(143, 168):
        # PG: 32 (sg 149)
        iso_sg = 149
    elif gid in range(168, 195):
        # PG: 622 (sg 177)
        iso_sg = 177
    elif gid in range(195, 231):
        # PG: 432 (sg 207)
        iso_sg = 207

    iso_trans = {'iso_rot': [], 'iso_tau': [], 'iso_mapping': []}
    iso_rot_candidates = SG_data[iso_sg - 1]['rotC']
    nop_iso = len(iso_rot_candidates) // 2
    iso_rot_candidates = iso_rot_candidates[0: nop_iso]

    if gid in [1, 2]:
        # triclinic system
        iso_rot_candidates.extend(find_triclinic_trial_iso())

    if gid in range(3, 16):
       #if not monoclinic_simplify: 
            # When finding eqv 3D reps, does not need to consider trial iso, and monoclinic_simplify=True
            # otherwise monoclinic_simplify=False and consider trial isomorphism with bound [-Max, Max+1)
        iso_rot_candidates.extend(find_monoclinic_trial_iso(Max=2)) # 50
       #iso_rot_candidates.extend(find_monoclinic_trial_iso(Max=12)) # 50
       #iso_rot_candidates.extend(find_monoclinic_trial_iso(Max=36)) # 50
       #iso_rot_candidates.extend(find_monoclinic_trial_iso(Max=50)) # 50

    print('conv candidates:', iso_rot_candidates) if test else None
    if norm(prim_basis - np.eye(3)) > tol:
        # transform iso_R to G prim basis, and keep only integer matrices
        iso_rot_candidates = [round_mat(inv(prim_basis) @ R @ prim_basis) for R in iso_rot_candidates]
        iso_rot_candidates = [R for R in iso_rot_candidates if vec_is_int(R)]
    print('num of iso_rot_candidates:', len(iso_rot_candidates)) if test else None
    print('candidates', iso_rot_candidates) if test else None

    for ith, iso_R in enumerate(iso_rot_candidates):        
        print('\n', ith, iso_R.tolist()) if test else None
        is_iso, iso_tau, iso_mapping = get_sg_iso_tau(sorted_rot_list, sorted_tau_list, iso_R)
        if is_iso:
            iso_R = prim_basis @ iso_R @ inv(prim_basis)
            conv_iso_tau = []
            for t in iso_tau:
                conv_iso_tau.append(prim_basis @ t)
            iso_trans['iso_rot'].append(iso_R)  # iso_rot in G prim latt
            iso_trans['iso_tau'].append(conv_iso_tau)  # iso_tau in G prim latt
            iso_trans['iso_mapping'].append(iso_mapping)
            print('find isomorphism:\n', iso_R, iso_tau) if test else None
        else:
            print('Not isomorphism!') if test else None

    print('total sol:', len(iso_trans['iso_rot'])) if test else None
    return iso_trans

