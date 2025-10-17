from spglib import standardize_cell,get_symmetry_dataset
import numpy as np
# ************************ function prim_cell ****************************
# 
# > Input a cell
# get the standard primitive cell of SG, without supercell 
# 
# ***********************************************************************
def prim_cell(cell,tol=1e-2): # FIX 
    lattice, positions, numbers, elements, mag = cell
    cell_spglib = (lattice, positions, numbers)
    cell_p = standardize_cell(cell_spglib, to_primitive=True, no_idealize=False,symprec=tol)
    # cell_p = find_primitive(cell_spglib,symprec=1e-2)
    lattice_p, positions_p, numbers_p = cell_p
    mag = np.zeros(np.shape(positions_p))
    cell_prim = (lattice_p, positions_p, numbers_p, elements, mag)

    return cell_prim



# ************************ function std_cell ****************************
# 
# > Input a cell
# get the standard conventional cell of SG, without supercell 
# 
# ***********************************************************************
def std_cell(cell,tol=1e-2): 
    lattice, positions, numbers, elements, mag = cell
    cell_spglib = (lattice, positions, numbers, mag)
    dataset = get_symmetry_dataset(cell_spglib,symprec=tol)
    # mag here not matter
    mag = np.zeros(np.shape(dataset['std_positions']))
    cell_std = (dataset['std_lattice'], dataset['std_positions'], dataset['std_types'], elements, mag)
    return dataset['number'],cell_std
