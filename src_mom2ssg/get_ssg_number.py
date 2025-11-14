import os
import sys
import numpy as np
from numpy.linalg import norm, det
from spglib import get_symmetry_dataset

from .small_func import *
from .SG_utils import identify_SG_lattice
from .SG_isomorphism import find_sg_iso_transform
from .find_ssg_operation import findAllOp, findAllOp_v2
from .load_ssgdata import load_ssg_list
from .eqvpg2label import hm_to_schoenflies, schoenflies_to_hm
from .poscar_io import read_poscar_no_elements

def search4ssg(cell, ssg_list, tol = 1e-4, tolm=1e-4):
    ssg_dict = findAllOp(cell, tol, tolm)
    # print(ssg_dict)
    search = str(ssg_dict['Gnum']) + '.' + str(ssg_dict['Ik']) + '.' + str(ssg_dict['It'])
    # Normalize equivalent HM forms; use Schoenflies for non-crystal groups.
    eqvPg = ssg_dict['QLabel']
    CRYSTAL_SCH = {
        'C1','Ci','C2','Cs','C2h','D2','C2v','D2h',
        'C4','S4','C4h','D4','C4v','D2d','D4h',
        'C3','C3i','D3','C3v','D3d',
        'C6','C3h','C6h','D6','C6v','D3h','D6h',
        'T','Th','O','Td','Oh'
    }
    try:
        sch = hm_to_schoenflies(str(eqvPg))
    except Exception:
        sch = None
    # Only convert to Schoenflies for non-crystal point groups; keep crystal HM unchanged
    if sch is not None and sch not in CRYSTAL_SCH:
        eqvPg = sch
    dim = findDimension(cell,tolm=tolm)
    H = ssg_dict['Hnum']
    ssg_maybe = []
    # print(search, H, dim, eqvPg)
    flag_maybe = True
    for s in ssg_list:
        if s['dim'] == dim and s['search'] == search:
            if s['Hid'] == int(H) and s['eqvPg'] == str(eqvPg):
                ssg_maybe.append(s)
                flag_maybe = False
    if flag_maybe:
        ssg_dict = findAllOp_v2(cell, tol,tolm)
        # print(ssg_dict)
        search = str(ssg_dict['Gnum']) + '.' + str(ssg_dict['Ik']) + '.' + str(ssg_dict['It'])
        eqvPg = ssg_dict['QLabel']
        try:
            sch = hm_to_schoenflies(str(eqvPg))
        except Exception:
            sch = None
        if sch is not None and sch not in CRYSTAL_SCH:
            eqvPg = sch
        dim = findDimension(cell,tolm=tolm)
        H = ssg_dict['Hnum']
        ssg_maybe = []
        # print(search, H, dim, eqvPg)
        for s in ssg_list:
            if s['dim'] == dim and s['search'] == search:
                if s['Hid'] == int(H) and s['eqvPg'] == str(eqvPg):
                    ssg_maybe.append(s)
    
    # print(search, H, dim, eqvPg)
    # print(len(ssg_maybe))
    if len(ssg_maybe) == 1:
        return ssg_maybe[0]['ssgNum']
    ssg_q = []
    out = []
    transform = ssg_dict['transformation_matrix']
    shift = ssg_dict['original_shift']
    Gnum = ssg_dict['Gnum']

    shift_list = [shift]

    trans_dict = find_sg_iso_transform(Gnum)
    w_iso_list = trans_dict['iso_rot']
    # print(w_iso_list)
    shift_iso_list = trans_dict['iso_tau']
    # print(shift_iso_list)
    transformation = ssg_dict['transformation_matrix']
    shift = ssg_dict['original_shift']
    # print(shift)
    # print(len(ssg_maybe))
    for ss in ssg_maybe:
        flag_s1 = True
        if flag_s1:
            for cnt, w_iso in enumerate(w_iso_list):
                transform = w_iso @ transformation
                tau_iso_list = shift_iso_list[cnt]
                for sh in shift_list:
                    if flag_s1:
                    # print(transform,sh)
                        checkq = checkQ(ss,ssg_dict, dim, sh, transform)
                        checkh = checkH(ss,ssg_dict, sh, transform)
                        # print(ss['ssgNum'],'Q and H',checkq,checkh)
                        if checkq and checkh:
                            # print(transform,sh)
                            flag_s1 = False
                            ssg_q.append(ss['ssgNum'])
                            out.append(ss)
    ssg_q = unique_list(ssg_q)
    # print(len(ssg_q))
    if len(ssg_q) == 1:
        return ssg_q[0]


    return 'need more loop'



def checkQ(ss, ssg_dict, dim, shift, transform):
    def check_rot(rot1, rot2, dim):
        if dim == 1:
            if abs(det(rot1) - det(rot2)) < 1e-2:
                return True
            else:
                return False

        diff_trace = abs(np.trace(rot1) - np.trace(rot2))
        if diff_trace < 1e-2:
            return True
        else:
            return False
    QRotC_std = ss['QRotC']
    QTauC_std = ss['QTauC']
    URot_std_list = ss['URot']
    QRotC = []
    QTauC = []
    # print('check q of ', ss['ssgNum'])
    # trnasform = ssg_dict['transformation_matrix']
    transform_inv = np.linalg.inv(transform)
    # shift = ssg_dict['original_shift']
    Rot_list = ssg_dict['RotC']
    Tau_list = ssg_dict['TauC']
    spin_list = ssg_dict['spin']

    for cnt, rot in enumerate(QRotC_std):
        tau = QTauC_std[cnt]
        QRotC.append(transform_inv@rot@transform)
        QTauC.append(transform_inv@(rot@shift+tau-shift))


    rep_T = 0
    for URot_std in URot_std_list:
        flag2 = True
        for cnt1, rot1 in enumerate(QRotC):
            flag = False
            # flag2 = True
            for cnt2, rot2 in enumerate(Rot_list):
                t1 = QTauC[cnt1]
                t2 = Tau_list[cnt2]
                mod_t1 = change_tau(t1)
                mod_t2 = change_tau(t2)
                if identity_vec(mod_t1, mod_t2) and norm(rot1-rot2) < 1e-3:
                    flag = True
                    spin_rot1 = spin_list[cnt2]
                    spin_rot2 = URot_std[cnt1]
                    spin_rot1 = np.array(spin_rot1)
                    spin_rot2 = np.array(spin_rot2)
                    if not check_rot(spin_rot1, spin_rot2, dim):
                        flag2 = False
            if not flag:
                return False
        if flag2:
            rep_T = rep_T + 1
    if rep_T == 1:
        return True
    if rep_T > 1:
        return True
    return False

def checkH(ss,ssg_dict, shift, transform):
    gid = ssg_dict['Gnum']
    HRotC_std = ss['HRotC']
    HTauC_std = ss['HTauC']
    tau_col = ss['superCell']
    t1 = np.array([tau_col[0][0], tau_col[1][0], tau_col[2][0]]) 
    t2 = np.array([tau_col[0][1], tau_col[1][1], tau_col[2][1]]) 
    t3 = np.array([tau_col[0][2], tau_col[1][2], tau_col[2][2]])
    prim_vec = identify_SG_lattice(gid)[1] # each col is a prim basis vector               
    for t_append in [t1, t2, t3]:
        HRotC_std.append(np.eye(3))
        HTauC_std.append(prim_vec @ t_append)
    transform_inv = np.linalg.inv(transform)
    Rot_list = ssg_dict['HRotC']
    Tau_list = ssg_dict['HTauC']
    HRotC = []
    HTauC = []
    for cnt, rot in enumerate(HRotC_std):
        tau = HTauC_std[cnt]
        HRotC.append(transform_inv@rot@transform)
        HTauC.append(transform_inv@(rot@shift+tau-shift))
    for cnt1, rot1 in enumerate(HRotC):
        flag = False
        for cnt2, rot2 in enumerate(Rot_list):
            t1 = HTauC[cnt1]
            t2 = Tau_list[cnt2]
            mod_t1 = change_tau(t1)
            mod_t2 = change_tau(t2)
            if identity_vec(mod_t1, mod_t2) and norm(rot1-rot2) < 1e-3:
                flag = True
        if not flag:
            return False
    return True
