from .small_func import *
from .eqvpg2label import get_std_pg
from pymatgen.core import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
import numpy as np
from numpy.linalg import norm, inv, det
from spglib import *

def generate_coord(rot_list, tau_list):
    """
    Generate coordinate list using symmetry operations {R|t}.
    
    Applies rotation matrices and translation vectors to generic test points
    to generate a complete set of equivalent positions.
    
    Args:
        rot_list: List of 3x3 rotation matrices
        tau_list: List of 3-element translation vectors
        
    Returns:
        list: List of generated 3D coordinates
    """
    # Generic test points to avoid special positions
    gen_pts = [
        np.array([0.1722, 0.6933, 0.9344]),
        np.array([0.8399, 0.5677, 0.0655]),
        np.array([0.1234, 0.4567, 0.0876]),
        np.array([0.2468, 0.5721, 0.7834])
    ]
    
    pos_list = []
    for R, t in zip(rot_list, tau_list):
        for p in gen_pts:
            pos_list.append(R @ p + t)
    
    return pos_list

def findQlabel(rot_list):
    """
    Find point group label from rotation matrices using molecular symmetry analysis.
    
    Args:
        rot_list: List of 3x3 rotation matrices
        
    Returns:
        str: Schoenflies point group symbol
    """
    rot_list = unique_matrix(rot_list)
    
    # Create molecular structure for point group analysis
    coord = [[0, 0, 0]]  # Origin point
    elements = ['C']     # Carbon at origin
    
    # Generic test points to avoid special positions
    gen_pts = [
        np.array([0.1722, 0.6933, 0.9344]),
        np.array([0.8399, 0.5677, 0.0655]),
        np.array([0.1234, 0.4567, 0.0876]),
        np.array([0.2468, 0.5721, 0.7834])
    ]
    
    # Apply rotations to generate equivalent positions
    for R in rot_list:
        for p in gen_pts:
            coord.append((R @ p.T).T)
            elements.append('H')
    
    # Analyze point group using pymatgen
    mol = Molecule(elements, coord)
    try:
        point_group = PointGroupAnalyzer(mol, tolerance=0.0001, 
                                    eigen_tolerance=0.001, 
                                    matrix_tolerance=0.001)
    except:
        point_group = PointGroupAnalyzer(mol, tolerance=0.0001, 
                                    eigen_tolerance=0.001)
    
    return point_group.sch_symbol
    



def findQ(cell,tolm=1e-4): # find the largest possilbel PG of the quotient group , get all the operations
    dim = findDimension(cell,tolm=tolm)

    mag = cell[3].copy()
    numbers = cell[2].copy()
    elements = []
    coord = []
    num = np.size(numbers, 0)
    for at in range(num):
        if not identity_vec(mag[at], [0, 0, 0]):
            coord.append(mag[at])
    # append [0, 0, 0] to the mass center

    coord.append(np.array([0., 0., 0.]))
    # coord = np.unique(coord, axis = 0)
    tol_mag = 1e-5 # tol for magnetism
    coord = np.array(coord)
    rounded_coord = np.round(coord / tol_mag) * tol_mag
    unique_coord, indices = np.unique(rounded_coord, axis=0, return_index=True)
    coord = coord[sorted(indices)]

    num_uni = np.size(coord, axis = 0)
    for i in range(num_uni):
        elements.append('H')

    # num2ele = {1: 'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6: 'C', 7:'N', 8:'O', 9:'F'}
    # for num in new_numbers:
    #     elements.append(num2ele[num])
    mol = Molecule(elements, coord)
    # print(elements, coord)
    try:
        point_group = PointGroupAnalyzer(mol, tolerance=tolm, eigen_tolerance=0.0001, matrix_tolerance=0.00001)
    except:
        point_group = PointGroupAnalyzer(mol, tolerance=tolm, eigen_tolerance=0.0001)
    # print('posible Q')
    # print(point_group.sch_symbol)
    op = point_group.get_symmetry_operations()
    op_list = []
    for i in op:
        op_list.append(i.rotation_matrix) # rotation operation list
    label = point_group.sch_symbol
    # print(op_list)
    E3 = np.eye(3)
    if dim == 2 or 3: #FIXME
        return op_list
    else:
        return [E3, -E3]



def findG(cell, tol):
    """
    Find non-magnetic space group information.
    
    Args:
        cell: Tuple of (lattice, positions, numbers, magmoms)
        tol: Tolerance for symmetry analysis
        
    Returns:
        list: [Gnum, rots, trans, supercell_multiple]
            - Gnum: Space group number
            - rots: Rotation matrices
            - trans: Translation vectors  
            - supercell_multiple: Volume ratio relative to primitive cell
    """
    # Create non-magnetic cell (ignore magnetic moments)
    cellG = (cell[0].copy(), cell[1].copy(), cell[2].copy())
    
    # Get space group number and symmetry operations
    Gnum = get_spacegroup(cellG)
    Gop = get_symmetry(cellG, symprec=tol)
    dataset = get_symmetry_dataset(cellG, symprec=tol)
    
    # Calculate supercell multiple from transformation matrix determinant
    supercell_multiple = abs(det(dataset['transformation_matrix']))
    
    return [Gnum, Gop['rotations'], Gop['translations'], supercell_multiple]



def change_shift(shift_orignal, sg): #change the std of spglib to the std of bilbao
    shift = np.zeros(3)
    if sg in [48, 86, 126, 210]:
        shift = np.array([-0.25, -0.25, -0.25])
    elif sg in [70, 201]:
        shift = np.array([-0.375, -0.375, -0.375])
    elif sg in [85, 129, 130]:
        shift = np.array([0.25, -0.25, 0])
    elif sg in [50, 59, 125]:
        shift = np.array([-0.25, -0.25, 0])
    elif sg in [133, 134, 137]:
        shift = np.array([0.25, -0.25, 0.25])
    elif sg in [141, 142]:
        shift = np.array([0.5, 0.25, 0.125])
    elif sg == 68:
        shift = np.array([0, -0.25, 0.25])
    elif sg == 88:
        shift = np.array([0, 0.25, 0.125])
    elif sg == 138:
        shift = np.array([0.25, -0.25, -0.25])
    elif sg in [222, 224]:
        shift = np.array([0.25, 0.25, 0.25])
    elif sg == 227:
        shift = np.array([0.125, 0.125, 0.125])
    else:
        shift = 0
    out = shift_orignal - shift
    return out

def findH(cell, tol):
    """
    Find the pure lattice space group number by grouping atoms with identical magnetic moments.
    
    Args:
        cell: Tuple of (lattice, positions, numbers, magmoms)
        tol: Tolerance for magnetic moment comparison
        
    Returns:
        list: [H_number, Hrot, Htran]
            - H_number: Pure lattice space group number
            - Hrot: Rotation matrices
            - Htran: Translation vectors
    """
    lattice, position, numbers, mag = cell
    lattice = lattice.copy()
    position = position.copy()
    numbers = numbers.copy()
    mag = mag.copy()
    
    num_atom = len(numbers)
    new_numbers = []
    numbers_basket = []
    mag_basket = []
    num_list = []
    num = 1
    
    # Group atoms by element type and magnetic moment
    for at in range(num_atom):
        if at == 0:
            # First atom
            numbers_basket.append(numbers[at])
            mag_basket.append(mag[at])
            new_numbers.append(1)
            num_list.append(1)
            num += 1
        else:
            # Check if atom matches existing group
            found_match = False
            for i in range(len(numbers_basket)):
                if (numbers[at] == numbers_basket[i] and 
                    identity_vec(mag[at], mag_basket[i], tol)):
                    new_numbers.append(num_list[i])
                    found_match = True
                    break
            
            if not found_match:
                # Create new group
                numbers_basket.append(numbers[at])
                mag_basket.append(mag[at])
                new_numbers.append(num)
                num_list.append(num)
                num += 1
    
    # Analyze symmetry of grouped structure
    cellH = (lattice, position, new_numbers)
    H_dataset = get_symmetry_dataset(cellH, symprec=tol)
    H_number = H_dataset['number']
    Hrot = H_dataset['rotations']
    Htran = H_dataset['translations']
    
    # Generate coordinates and verify
    pos_gen = generate_coord(Hrot, Htran)
    numbers = np.arange(4).tolist() * len(Hrot)
    cellH_new = (lattice, pos_gen, numbers)
    H_dataset = get_symmetry_dataset(cellH_new)
    H_number = H_dataset['number']
    
    return [H_number, Hrot, Htran]

def findH_v2(out_spin, rot, tau, cell):
    lattice = cell[0].copy()
    rot_list = []
    tau_list = []
    for cnt, spin in enumerate(out_spin):
        if norm(spin - np.eye(3)) < 1e-3:
            rot_list.append(rot[cnt])
            tau_list.append(tau[cnt])
    position = generate_coord(rot_list, tau_list)
    numbers = np.arange(4).tolist() * np.size(rot_list, axis=0)
    cellH= (lattice, position, numbers)
    dataset = get_symmetry_dataset(cellH,symprec=1e-5)
    # print(dataset)
    # print(cellG)
    # print(rot_list)
    # print(trans_list)
    Gnum = int(dataset['number'])
    return [Gnum, rot_list, tau_list]


def std_cell(cell): # standard the cell by sort the position on each atom
    def position_and_magmom(positions, magmoms):
        assert np.size(positions) == np.size(magmoms)
        num_atom = np.size(positions, 0)
        pos_mag = np.zeros((num_atom, 6))
        pos_mag[:, 0:3] = positions
        pos_mag[:,3:] = magmoms
        return pos_mag 
    def mod_pos(cell):
        lattice = cell[0].copy()
        position = cell[1].copy()
        numbers = cell[2].copy()
        mag = cell[3].copy()
        pos_mag = position_and_magmom(position, mag)
        num_atom = np.size(position, 0)
        poss = []
        for at in range(num_atom):
            pos = position[at]
            mod_position = [(pos[0]+1e-6)%1, (pos[1]+1e-6)%1, (pos[2]+1e-6)%1]
            poss.append(mod_position)

        poss = np.array(poss)
        cell_new = (lattice, poss, numbers, mag)
        return cell_new
    def get_deg(numbers):
        unique_numbers = sorted(set(numbers))
        num_deg = []
        number_ar = np.array(numbers)
        for n in unique_numbers:
            num_deg.append(len(np.nonzero(number_ar == n)[0]))
        return num_deg
    cell = mod_pos(cell)
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[3].copy()
    pos_mag = position_and_magmom(position, mag)
    

    unique_numbers = sorted(set(numbers))
    sorted_positions = []
    sorted_magmoms = []
    sorted_numbers = []
    
    for atomic_num in unique_numbers:

        indices = [i for i, num in enumerate(numbers) if num == atomic_num]

        atom_positions = position[indices]

        if isinstance(mag, list):
            atom_magmoms = np.array([mag[i] for i in indices])
        else:
            atom_magmoms = mag[indices]

        atom_pos_mag = position_and_magmom(atom_positions, atom_magmoms)
        sorted_atom_pos_mag = sort_array(atom_pos_mag)

        sorted_positions.extend(sorted_atom_pos_mag[:, 0:3])
        sorted_magmoms.extend(sorted_atom_pos_mag[:, 3:])
        sorted_numbers.extend([atomic_num] * len(indices))
    

    sorted_positions = np.array(sorted_positions)
    sorted_magmoms = np.array(sorted_magmoms)
    sorted_numbers = np.array(sorted_numbers)

    cell_standard = (lattice, sorted_positions, sorted_numbers, sorted_magmoms)

    return cell_standard



def findGnum(rot_list, trans_list, lattice):
    """
    Find space group number and transformation information from symmetry operations.
    
    Args:
        rot_list: List of rotation matrices
        trans_list: List of translation vectors
        lattice: Lattice vectors (3x3 matrix)
        
    Returns:
        list: [Gnum, transform_matrix, shift]
            - Gnum: Space group number
            - transform_matrix: Transformation matrix to standard setting
            - shift: Origin shift vector
    """
    # Generate coordinates from symmetry operations
    position = generate_coord(rot_list, trans_list)
    numbers = np.arange(4).tolist() * np.size(rot_list, axis=0)
    cellG = (lattice, position, numbers)
    
    # Get symmetry dataset from spglib
    dataset = get_symmetry_dataset(cellG, symprec=1e-4)
    
    Gnum = int(dataset['number'])
    transform_matrix = dataset['transformation_matrix']
    shift = change_shift(dataset['origin_shift'], Gnum)
    
    return [Gnum, transform_matrix, shift]



def identity_cell(cell1, cell2, tolm=1e-4): 
    cell1_std = std_cell(cell1)
    cell2_std = std_cell(cell2)

    mag1 = cell1_std[3].copy()
    mag2 = cell2_std[3].copy()
    for cnt, m1 in enumerate(mag1):
        m2 = mag2[cnt]
        if not identity_vec(m1, m2, tolm):
            return False
    return True


def twoD2threeD(cell):
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[3].copy()
    nonzero = []
    for m in mag:
        if not identity_vec(m, [0, 0, 0]):
            nonzero.append(m)
    nonzero = unique_matrix(nonzero)
    assert np.size(nonzero, axis=0) > 1
    for i in nonzero[1:]: 
        p_vec = np.cross(nonzero[0], i)
        if not identity_vec(p_vec, [0, 0, 0]):
            p_vec = p_vec * norm(i) / norm(p_vec)
            break
    mag_new = []
    for m in mag:
        mag_new.append(m + p_vec)
    cell_new = (lattice, position, numbers, mag_new)
    # print(p_vec)
    return cell_new



def findAllOp(cell, tol=1e-3,tolm=1e-4):
    dim = findDimension(cell,tolm=tolm)
    # id dim = 2, transfer it into a 3d which will maintain all the operation, but still .P
    if dim == 2:
        cell = twoD2threeD(cell)
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[3].copy()
    out_spin = []
    out_rot = []
    out_tran = []
    spinOp = findQ(cell,tolm=1e-4)  # find the spin-point group operations
    latticeOp = findG(cell, tol)
    non_mag_G = latticeOp[0]
    rot = latticeOp[1]
    tran = latticeOp[2]
    supercell_multiple = latticeOp[3]
    assert np.size(rot, 0) == np.size(tran, 0)
    sg_op_num = np.size(rot, 0)
    for i in range(sg_op_num):
        r = rot[i]
        t = tran[i]
        for spin_rot in spinOp:
            pos_new = []
            mag_new = []
            for pos in position:
                pos_new.append((r @ pos.T + t).T)
            for m in mag:
                mag_new.append((spin_rot @ m.T).T)
            # print(t)
            # print(spin_rot)
            cell_new = (lattice, pos_new, numbers, mag_new)
            if identity_cell(cell, cell_new, tolm):
                out_spin.append(spin_rot)
                out_rot.append(r)
                out_tran.append(t)
    # using the dimension to identify the ssg
    fH = findH(cell, tol)
    Hnum = fH[0]
    Hrot = fH[1]
    Htran = fH[2]
    H_op_num = np.size(Hrot, axis=0)
    G_op_num = np.size(out_rot, axis=0)
    H_point_num = np.size(unique_matrix(Hrot), axis=0)
    G_point_num = np.size(unique_matrix(out_rot), axis=0)
    # print(H_point_num, G_point_num)
    # get G number
    fG = findGnum(out_rot, out_tran, lattice)
    Gnum = fG[0]
    transform = fG[1]
    shift = fG[2]
    # get ik and it
    It = int(round(G_point_num/H_point_num))
    IKIT = G_op_num/H_op_num
    assert abs(IKIT - round(IKIT)) < 1e-3
    IKIT = round(IKIT)
    Ik = IKIT/It
    assert abs(Ik - round(Ik)) < 1e-3
    Ik = round(Ik)
    # find Q
    # print(out_spin)
    spin_uniqe = unique_matrix(out_spin)
    # print(spin_uniqe)
    # when dim == 1, it can only be fermmomagnetic or anti-ferromagnetic
    if dim == 1:
        if np.size(spin_uniqe, axis=0) == 1:
            Gnum = Hnum
            assert It == 1
            Ik = 1
            spin_label = 'C1'
        else:
            assert np.size(spin_uniqe, axis=0) == 2
            assert It == 1 or 2
            Ik = int(2/It)
            spin_label = 'Cs'

    else:
        spin_label = findQlabel(out_spin)
        # print(spin_label)

    return {'spin': out_spin, 'It': It, 'Hnum': Hnum, 'Ik': Ik, 'Gnum': Gnum, 'QLabel': get_std_pg(spin_label)[0], 'RotC': out_rot,
            'TauC': out_tran, 'transformation_matrix': transform, 'original_shift': shift, 'HRotC': Hrot, 'HTauC': Htran}

def findAllOp_v2(cell, tol,tolm=1e-4):
    dim = findDimension(cell,tolm=tolm)
    # id dim = 2, transfer it into a 3d which will maintain all the operation, but still .P
    if dim == 2:
        cell = twoD2threeD(cell)
    lattice = cell[0].copy()
    position = cell[1].copy()
    numbers = cell[2].copy()
    mag = cell[3].copy()
    out_spin = []
    out_rot = []
    out_tran = []
    spinOp = findQ(cell,tolm)
    latticeOp = findG(cell, tol)
    non_mag_G = latticeOp[0]
    rot = latticeOp[1]
    tran = latticeOp[2]
    supercell_multiple = latticeOp[3]
    assert np.size(rot, 0) == np.size(tran, 0)
    sg_op_num = np.size(rot, 0)
    for i in range(sg_op_num):
        r = rot[i]
        t = tran[i]
        for spin_rot in spinOp:
            pos_new = []
            mag_new = []
            for pos in position:
                pos_new.append((r @ pos.T + t).T)
            for m in mag:
                mag_new.append((spin_rot @ m.T).T)
            # print(t)
            # print(spin_rot)
            cell_new = (lattice, pos_new, numbers, mag_new)
            if identity_cell(cell, cell_new,tolm):
                out_spin.append(spin_rot)
                out_rot.append(r)
                out_tran.append(t)

    fH = findH_v2(out_spin, out_rot, out_tran, cell)
    Hnum = fH[0]
    Hrot = fH[1]
    Htran = fH[2]
    H_op_num = np.size(Hrot, axis=0)
    G_op_num = np.size(out_rot, axis=0)
    H_point_num = np.size(unique_matrix(Hrot), axis=0)
    G_point_num = np.size(unique_matrix(out_rot), axis=0)
    # print(H_point_num, G_point_num)
    # get G number
    fG = findGnum(out_rot, out_tran, lattice)
    Gnum = fG[0]
    transform = fG[1]
    shift = fG[2]
    # get ik and it
    It = int(round(G_point_num/H_point_num))
    IKIT = G_op_num/H_op_num
    assert abs(IKIT - round(IKIT)) < 1e-3
    IKIT = round(IKIT)
    Ik = IKIT/It
    assert abs(Ik - round(Ik)) < 1e-3
    Ik = round(Ik)
    # find Q
    # print(out_spin)
    spin_uniqe = unique_matrix(out_spin)
    # print(spin_uniqe)
    # when dim == 1, it can only be fermmomagnetic or anti-ferromagnetic
    if dim == 1:
        if np.size(spin_uniqe, axis=0) == 1:
            Gnum = Hnum
            assert It == 1
            Ik = 1
            spin_label = 'C1'
        else:
            assert np.size(spin_uniqe, axis=0) == 2
            assert It == 1 or 2
            Ik = int(2/It)
            spin_label = 'Cs'

    else:
        spin_label = findQlabel(out_spin)
        # print(spin_label)

    return {'spin': out_spin, 'It': It, 'Hnum': Hnum, 'Ik': Ik, 'Gnum': Gnum, 'QLabel': get_std_pg(spin_label)[0], 'RotC': out_rot,
            'TauC': out_tran, 'transformation_matrix': transform, 'original_shift': shift, 'HRotC': Hrot, 'HTauC': Htran}

