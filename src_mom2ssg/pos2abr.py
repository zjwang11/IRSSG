import os
import numpy as np
from scipy.io import FortranFile 
from pymatgen.io.vasp import Poscar
from pymatgen.core import Lattice, Structure
from spglib import get_symmetry_dataset
import struct
import sys
import copy
def nint_fortran(x: np.ndarray) -> np.ndarray:
    """
    Mimic Fortran's NINT:  half-away-from-zero rounding.
    """
    return np.where(x >= 0.0,
                    np.floor(x + 0.5),
                    np.ceil (x - 0.5)).astype(int)
    
def analyze_symmetry(cell):
    if len(cell) == 5:
        lattice, positions, numbers, elements, mag = cell
        cell_spglib = (lattice, positions, numbers, mag)
    else:
        cell_spglib = copy.deepcopy(cell)
    
    dataset     = get_symmetry_dataset(cell_spglib,symprec=1e-2)
    spg_number  = int(dataset["number"])
    spg_symbol  = dataset["international"] 

    rotations   = dataset["rotations"]
    translations = dataset["translations"]


    return spg_number, spg_symbol, rotations, translations

TOL = 1.0e-4
INV3  = 1.0 / 3.0
DINV3 = 2.0 / 3.0
Pabc = np.array([[ 1.0 ,  0.0 ,  0.0 ],
                 [ 0.0 ,  1.0 ,  0.0 ],
                 [ 0.0 ,  0.0 ,  1.0 ]]).T

Cabc = np.array([[ 0.5 , -0.5 ,  0.0 ],
                 [ 0.5 ,  0.5 ,  0.0 ],
                 [ 0.0 ,  0.0 ,  1.0 ]]).T

Cabc68 = np.array([[ 0.5 ,  0.5 ,  0.0 ],
                   [-0.5 ,  0.5 ,  0.0 ],
                   [ 0.0 ,  0.0 ,  1.0 ]]).T

Babc = np.array([[ 0.5 ,  0.5 ,  0.0 ],
                 [-0.5 ,  0.5 ,  0.0 ],
                 [ 0.0 ,  0.0 ,  1.0 ]]).T

Aabc = np.array([[ 1.0 ,  0.0 ,  0.0 ],
                 [ 0.0 ,  0.5 ,  0.5 ],
                 [ 0.0 , -0.5 ,  0.5 ]]).T

Rabc = np.array([[  DINV3 ,  INV3 ,  INV3 ],
                 [ -INV3 ,  INV3 ,  INV3 ],
                 [ -INV3 , -DINV3,  INV3 ]]).T

Fabc = np.array([[ 0.0 , 0.5 , 0.5 ],
                 [ 0.5 , 0.0 , 0.5 ],
                 [ 0.5 , 0.5 , 0.0 ]]).T

Iabc = np.array([[-0.5 , 0.5 , 0.5 ],
                 [ 0.5 , -0.5, 0.5 ],
                 [ 0.5 , 0.5 , -0.5]]).T


def read_fortran_unformatted(filename, num_sym):
    rot_table = np.zeros((3,3,num_sym), dtype=np.int32)
    trans_table = np.zeros((3,num_sym), dtype=np.float64)

    with open(filename, 'rb') as f:
        for ir in range(num_sym):
            reclen_bytes = f.read(4)
            if not reclen_bytes:
                break
            reclen = struct.unpack('i', reclen_bytes)[0]

            rot_data = np.fromfile(f, dtype=np.int32, count=9)
            trans_data = np.fromfile(f, dtype=np.float64, count=3)
            rot_table[:, :, ir] = rot_data.reshape(3, 3).T
            trans_table[:, ir] = trans_data

            f.read(4)

    return rot_table.T, trans_table.T


def crys_conv(space_group      : int,
              international    : str,
              positions        : np.ndarray,   # shape (3, num_atom)
              rotations        : np.ndarray,   # shape (3,3,num_sym)  (int)
              translations     : np.ndarray,   # shape (3,num_sym)
              tables_dir       : str = "ssg_data/TablesSymElemSG"):

    positions   = np.asarray(positions,   dtype=np.float64)
    rotations   = np.asarray(rotations,   dtype=np.int64)
    translations = np.asarray(translations, dtype=np.float64)

    num_atom = positions.shape[1]
    num_sym  = rotations.shape[0]

    csgn = f"{space_group:d}".rjust(3).strip()  
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir,os.path.join(tables_dir, f"SymElemSG_{csgn}.data"))

    rot_table   = np.empty_like(rotations)
    trans_table = np.empty_like(translations)

    
    rot_table, trans_table = read_fortran_unformatted(path, num_sym)


    fc = international[0].upper()

    Kc2p = {'P': Pabc, 'C': Cabc, 'B': Babc, 'A': Aabc,
            'R': Rabc, 'F': Fabc, 'I': Iabc}.get(fc)


    if Kc2p is None:
        raise ValueError(f"Un-recognised Bravais symbol '{fc}'")

    if space_group == 68 and fc == 'C':
        Kc2p = Cabc68.copy()

    latt_center = np.zeros((4, 3))
    if fc == 'C':
        latt_center[3] = (0.5, 0.5, 0.0)
    elif fc == 'B':
        latt_center[3] = (0.5, 0.0, 0.5)
    elif fc == 'A':
        latt_center[3] = (0.0, 0.5, 0.5)
    elif fc == 'R':
        latt_center[2] = (2/3, 1/3, 1/3)
        latt_center[3] = (1/3, 2/3, 2/3)
    elif fc == 'F':
        latt_center[1] = (0.0, 0.5, 0.5)
        latt_center[2] = (0.5, 0.0, 0.5)
        latt_center[3] = (0.5, 0.5, 0.0)
    elif fc == 'I':
        latt_center[3] = (0.5, 0.5, 0.5)

    p2cR = np.linalg.inv(Kc2p)


    translations = (translations + 1e-7) % 1.0

    trans_conv = np.matmul(Kc2p, translations.T).T

    rot_conv = np.zeros_like(rotations)
    for i in range(num_sym):
        rot_conv[i]   = np.rint(Kc2p@rotations[i]@p2cR).astype(int)
        
    spg2tab = np.zeros(num_sym, dtype=int)
    for ir in range(num_sym):
        for jr in range(num_sym):
            if np.array_equal(rot_table[jr,...], rot_conv[ir,...]):
                spg2tab[ir] = jr
                break

    spginv_loc = next((ir for ir in range(num_sym)
                       if abs(rot_conv[ir,0,0] + rot_conv[ir,1,1] + rot_conv[ir,2,2] + 3) < 1e-3), None)
    tabinv_loc = spg2tab[spginv_loc] if spginv_loc is not None else None
    

    shift = np.zeros(3)
    if spginv_loc is not None:
        shift = 0.5 * (trans_table[tabinv_loc,:] - trans_conv[spginv_loc,:])

    choices = np.array([[0,0,0], [0.5,0,0], [0,0.5,0], [0.5,0.5,0],
                        [0,0,0.5], [0.5,0,0.5], [0,0.5,0.5], [0.5,0.5,0.5]])
    possible_shift = shift[None,:] + choices   # (3,8)

    TOL = 1e-5
    match_found = False
    
    for sh in possible_shift:                # 8 个候选 shift
        Rs  = (rot_conv @ sh)                # (N,3)  每行 R·s
        T_p = ((trans_conv + sh - Rs + 1e-8) % 1.0 )  # (N,3)

        for ir in range(num_sym):            # 每个对称操作
            tgt = trans_table[spg2tab[ir]]   # 目标 (3,)

            ok = False
            for iv in range(4):
                # diff = latt_center[iv] + T_p[ir] - tgt   # (3,)
                # err  = np.abs(np.round(diff) - diff).sum()
                diffvec = latt_center[iv] + T_p[ir] - tgt   # (3,)
                
                nint    = nint_fortran(diffvec)
                
                err     = np.abs(nint - diffvec).sum()         
                if err <= TOL:       
                    ok = True
                    break             
            if not ok:
                break                  
        else:
            shift = sh.copy()
            match_found = True
            break

    if not match_found:
        raise RuntimeError("Error!")


    pos_conv  = Kc2p @ positions.T
    pos_conv += shift[:,None]
    pos_shift = p2cR @ pos_conv

    return pos_shift,p2cR @ shift

def pos2abr(cell):    
    spg_no, spg_sym, R, T = analyze_symmetry(cell)

    frac_coords = cell[1] 


    pos_shift,shift = crys_conv(spg_no,spg_sym,frac_coords,R,T)
    if len(cell) == 5:
        return (cell[0],pos_shift.T,cell[2],cell[3],cell[4]),shift
    else:
        return (cell[0],pos_shift.T,cell[2],cell[3]),shift
