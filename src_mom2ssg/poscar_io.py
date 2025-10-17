import numpy as np
from numpy.linalg import norm
from pymatgen.io.cif import CifParser
from collections import OrderedDict
# ************************ function read_poscar ****************************
# 
# > Input POSCAR, can with MAGMOM   x  y  z  mz  my  mz
# > (lattice, positions, numbers, mag) 
# > (lattice, positions, numbers, mag) = read_poscar() 
# 
# ***********************************************************************
def read_poscar(file_name = 'POSCAR'): 
    lattice  = np.zeros((3, 3))
    deg_of_points = []
    numbers = []
    mag = []
    # file_name = 'POSCAR_test'
    with open(file_name, 'r') as f:
        f_lines = f.readlines()
        for cnt, line in enumerate(f_lines):
            # cnt = 0,1 are meaningless
            # cnt = 2,3,4 show the basis vectors
            if cnt in [2, 3, 4]:
                line = line.strip().split()
                lattice[cnt-2, :] = [float(c) for c in line]
            if cnt  == 5:
                elements = line.strip().split()
            if cnt  == 6:
                line = line.strip().split()
                deg_of_points = [int(c) for c in line]
                num_of_points = sum(deg_of_points)
                num_of_types = np.size(deg_of_points)
                # get the numbers of cell
                for m in range(num_of_types):
                    for n in range(deg_of_points[m]):
                        numbers.append(m+1)
                positions = np.zeros((num_of_points, 3))
            if cnt > 7 and cnt < (7 + num_of_points + 1):
                line = line.strip().split()
                if len(line) == 3:
                    positions[cnt-8, :] = [float(c) for c in line]
                    mag.append([0, 0, 0])
                elif len(line) == 6:
                    positions[cnt-8, :] = [float(c) for c in line[:3]]
                    mag.append([float(c) for c in line[3:]])
                else:
                    # print('line', cnt+1 ,'not valid POSCAR input !!!')
                    raise ValueError('Invalid input, please check POSCAR line ' + str(cnt+1) )
        cell = (lattice, positions, numbers, elements, np.array(mag)) 
    return cell


def read_poscar_no_elements(file_name='POSCAR'):
    """
    Read POSCAR file to extract crystal structure and magnetic moments.
    
    Args:
        file_name: Path to POSCAR file (default: 'POSCAR')
        
    Returns:
        tuple: (cell, elements)
            - cell: Tuple of (lattice, positions, numbers, magmoms)
            - elements: List of element symbols
    """
    lattice = np.zeros((3, 3))
    deg_of_points = []
    numbers = []
    mag = []
    
    with open(file_name, 'r') as f:
        f_lines = f.readlines()
        
        for cnt, line in enumerate(f_lines):
            # Lines 2-4: Lattice vectors
            if cnt in [2, 3, 4]:
                line = line.strip().split()
                lattice[cnt-2, :] = [float(c) for c in line]
            
            # Line 5: Element symbols
            elif cnt == 5:
                elements = line.strip().split()
            
            # Line 6: Number of atoms per element
            elif cnt == 6:
                line = line.strip().split()
                deg_of_points = [int(c) for c in line]
                num_of_points = sum(deg_of_points)
                num_of_types = len(deg_of_points)
                
                # Generate atom type numbers
                for m in range(num_of_types):
                    for n in range(deg_of_points[m]):
                        numbers.append(m + 1)
                
                positions = np.zeros((num_of_points, 3))
            
            # Lines 8+: Atomic positions and magnetic moments
            elif cnt > 7 and cnt < (7 + num_of_points + 1):
                line = line.strip().split()
                atom_idx = cnt - 8
                
                if len(line) == 3:
                    # Position only (non-magnetic)
                    positions[atom_idx, :] = [float(c) for c in line]
                    mag.append([0, 0, 0])
                elif len(line) == 6:
                    # Position + magnetic moment
                    positions[atom_idx, :] = [float(c) for c in line[:3]]
                    mag.append([float(c) for c in line[3:]])
                else:
                    raise ValueError(f'Invalid input format at line {cnt+1}')
    
    cell = (lattice, positions, numbers, np.array(mag))
    return cell, elements


# ************************ function write_poscar ****************************
# 
# > input cell = (lattice, positions, numbers, mag) 
# > write a POSCAR file using the cell 
# 
# ***********************************************************************
def write_poscar(cell, file_name='POSCAR.pos2ssg'):

    lattice, positions, numbers,element, mag = cell
    
    positions = np.array(positions)
    mag = np.array(mag)

    # get numbers of the cell
    unique_numbers = []
    for n in numbers:
        if np.int32(n) not in unique_numbers:
            unique_numbers.append(np.int32(n))
    
    numbers = [[np.int32(nums)] for nums in numbers]

    # numbers of elements
    element_counts = [numbers.count(t) for t in unique_numbers]
    # write POSCAR
    with open(file_name, 'w') as f:
        f.write("Generated by write_poscar\n")
        f.write("1.0\n")
        for vec in lattice:
            f.write("  {:.16f}  {:>.16f}  {:>.16f}\n".format(*vec))
        f.write("  " + "  ".join(t for t in element) + "\n")
        f.write("  " + "  ".join(str(c) for c in element_counts) + "\n")
        f.write("Direc\n")

        for i in range(len(positions)):
            pos_str = "  {: .16f}  {: .16f}  {: .16f}".format(*positions[i])
            if mag is not None and len(mag) > 0:
                mag_str = "  {: .6f}  {: .6f}  {: .6f}".format(*mag[i])
                if norm(mag[i]) > 1e-5:
                    f.write(pos_str + mag_str + "\n")
                else:
                    f.write(pos_str + "\n")
            else:
                f.write(pos_str + "\n")


def mcif2cell(file_name):
    cc = CifParser(file_name)
    structure = cc.get_structures(primitive=False)[0]
    atom_types = [site.species_string for site in structure]
    atom_map = list(OrderedDict.fromkeys(atom_types))
    element_map = {element: i+1 for i, element in enumerate(atom_map)}
    numbers = [element_map[site.species_string] for site in structure]
    elements = list(structure.symbol_set)

    lattice = structure.lattice.matrix
    position = structure.frac_coords
    m = structure.site_properties['magmom']
    mag = []
    for i in m:
        mag.append(i.moment)
    
    # Get unique element symbols (no duplicates)
    element_symbols = list(OrderedDict.fromkeys([site.species_string for site in structure]))
    
    cell = (lattice, position, numbers, mag)
    return cell, element_symbols