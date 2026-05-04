
# Element periodic table mapping (atomic number to element symbol)
ELEMENT_TABLE = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
    19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
    49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
    55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr'
}

# Reverse mapping (element symbol to atomic number)
ELEMENT_TO_NUMBER = {v: k for k, v in ELEMENT_TABLE.items()}

magmoms_3d = {'Ti':3.0, 'V':3.0, 'Cr':3.0, 'Mn':3.0, 'Fe':3.0, 'Co':3.0, 'Ni':3.0, 'Cu':3.0}
magmoms_4d = {'Y':3.0, 'Ru':3.0, 'Rh':3.0, 'Os':2.0, 'Ir':0.5}
magmoms_4f = {'Ce':3.0, 'Pr':3.0, 'Nd':3.0, 'Pm':3.0, 'Sm':3.0, 'Eu':7.0, 'Gd':8.0,
                'Tb':8.0, 'Dy':8.0, 'Ho':8.0, 'Er':6.0, 'Tm':6.0, 'Yb':6.0}
magmoms_5f = {'U':2.0, 'Np':2.0, 'Pu':5.0}
magmom_dict_default = {**magmoms_3d, **magmoms_4d, **magmoms_4f, **magmoms_5f}


def construct_magmom_dict(atom_list,mag_list):
    if atom_list == []:
        return magmom_dict_default
    if len(atom_list) != len(mag_list):
        raise ValueError("atom_list and mag_list must have the same length")
    magmom_dict = {}
    for atom,mag in zip(atom_list,mag_list):
        magmom_dict[atom] = mag
    return magmom_dict

# ************** function analyze_magnetic_atoms ************************
# 
# from the cell get the magnetic atoms
# get the magnetitute of the magnetic moments
#
# ***********************************************************************
def analyze_magnetic_atoms(cell,magmom_dict = magmom_dict_default):
    lattice, positions, numbers, mag = cell

    # mapping atomic number to element symbol
    type_id_list = sorted(set(numbers))
    typeid_to_element = {atomic_num: ELEMENT_TABLE[atomic_num] for atomic_num in type_id_list}

    # find the magnetic moments type
    magnetic_species = {}
    magnetic_index = []
    for type_id, elem in typeid_to_element.items():
        if elem in magmom_dict:
            magnetic_species[type_id] = magmom_dict[elem]
            magnetic_index.append(type_id)

    # find which magnetic moments
    magnetic_atom_indices = []
    magnetic_atom_moments = []

    for i, type_id in enumerate(numbers):
        if type_id in magnetic_species:
            magnetic_atom_indices.append(i)
            magnetic_atom_moments.append(magnetic_species[type_id])

    return {
        'magnetic_species_count': len(magnetic_species),
        'magnetic_index': magnetic_index,
        'magnetic_species_moments': magnetic_species,  # {type_id: mag}
        'magnetic_atom_indices': magnetic_atom_indices,  # atom index
        'magnetic_atom_moments': magnetic_atom_moments  # same length as indices
    }
