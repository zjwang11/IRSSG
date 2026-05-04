from spglib import get_symmetry_dataset,find_primitive
import numpy as np
from numpy.linalg import inv, det

# ************************ function std_cell_pos2ssg ****************************
# 
# > Input a cell
# get the standard conventional cell of SG, without supercell 
# 
# ***********************************************************************
def std_cell_pos2ssg(cell,tol=1e-2): 
    lattice, positions, numbers, mag = cell
    cell_spglib = (lattice, positions, numbers, mag)
    dataset = get_symmetry_dataset(cell_spglib,symprec=tol)
    # print(dataset['number'])
    # mag here not matter
    mag = np.zeros(np.shape(dataset['std_positions']))
    cell_std = (dataset['std_lattice'], dataset['std_positions'], dataset['std_types'], mag)
    return cell_std


# ************************ function prim_cell ****************************
# 
# > Input a cell
# get the standard primitive cell of SG, without supercell 
# 
# ***********************************************************************
def prim_cell(cell, tol=1e-3): # FIX 
    lattice, positions, numbers, mag = cell
    cell_spglib = (lattice, positions, numbers, mag)
    cell_p = find_primitive(cell_spglib, symprec=tol)
    lattice_p, positions_p, numbers_p = cell_p
    mag = np.zeros(np.shape(positions_p))
    cell_prim = (lattice_p, positions_p, numbers_p, mag)
    return cell_prim


# ************************ function extend_cell ****************************
# 
# > from a primitive Cell
# get the standard POSCAR after considering the supercell
# 
# ***********************************************************************
# =========================================================
def extend_cell(cell_prim, superCell):
    tau_col = superCell
    t1 = np.array([tau_col[0][0], tau_col[1][0], tau_col[2][0]]) 
    t2 = np.array([tau_col[0][1], tau_col[1][1], tau_col[2][1]]) 
    t3 = np.array([tau_col[0][2], tau_col[1][2], tau_col[2][2]])
    lattice_p, positions_p, numbers_p, mag = cell_prim
    # superCell lattice
    a_p, b_p, c_p = lattice_p[0], lattice_p[1], lattice_p[2]
    a_s = tau_col[0][0] * a_p + tau_col[1][0] * b_p + tau_col[2][0] * c_p
    b_s = tau_col[0][1] * a_p + tau_col[1][1] * b_p + tau_col[2][1] * c_p
    c_s = tau_col[0][2] * a_p + tau_col[1][2] * b_p + tau_col[2][2] * c_p
    lattice_s = [a_s, b_s, c_s]
    # Convert fractional coordinates in primitive cell to Cartesian coordinates
    # pos_cart = [pos[0] * a_p + pos[1] * b_p + pos[2] * c_p for pos in positions_p]

    # Calculate determinant of the supercell matrix (volume ratio)
    detSuperCell = round(det(superCell))

    # Precompute supercell lattice matrix for fractional conversion
    lattice_matrix = np.array(lattice_s).T
    eps = 1e-5

    # Enumerate all possible positions within supercell bounds
    positions_s = []
    numbers_s = []
    count = 0
    for pos, num in zip(positions_p, numbers_p):
        # print('==========')
        for ix in range(-(detSuperCell + 2), detSuperCell + 3):
            for iy in range(-(detSuperCell + 2), detSuperCell + 3):
                for iz in range(-(detSuperCell + 2), detSuperCell + 3):
                    pp = [pos[0] +ix, pos[1] + iy, pos[2] + iz]
                    # print(pp)
                    # print(lattice_matrix)
                    fractional = inv(superCell) @ pp
                    # print(fractional)

                    # Clean near-zero components for numerical stability
                    fractional[np.abs(fractional) < eps] = 0.0
                    # print(fractional)

                    # Check if atom lies within the [0,1) range in all directions
                    if np.all((fractional >= -eps) & (fractional < 1.0 - eps)):
                        positions_s.append(fractional)
                        # print('========================')
                        # print(fractional)
                        numbers_s.append(num)
                        count += 1

    # Validate expected number of atoms in supercell
    # print(lattice_s, positions_s, numbers_s)
    if count != detSuperCell * len(numbers_p):
        print(count)
        print(detSuperCell, len(numbers_p))
        raise ValueError('Number of atoms in the supercell is incorrect! Check function extend_cell.')

    # Initialize magnetization (placeholder: all zero)
    mag_s = np.zeros(np.shape(positions_s))

    cell_s = (lattice_s, positions_s, numbers_s, mag_s)
    return cell_s

