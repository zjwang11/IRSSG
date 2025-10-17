#!/usr/bin/env python3
"""
SSG Label Module

This module provides the ssgGroup class for handling Spin Space Group (SSG) operations
and labeling. It manages SSG data structures, symmetry operations, and provides
functionality for generating SSG labels and operations.

Main functionality:
- ssgGroup class for SSG data management
- SSG operation loading and processing
- SU2 matrix calculations for spin operations
- SSG labeling and output generation
- Anti-unitary operation handling

Author: MOM2SSG Team
"""

import numpy as np
from numpy.linalg import norm, det
from math import acos, sin, pi

from .SG_utils import *
from .eqvpg2label import get_std_pg
from .small_func import SU2,round_vec

def identify_lattice_type(gid):
    # identify Bravais lattice for a given sg, return: Brav_latt
    SGTricP = [1, 2]
    SGMonoP = [3, 4, 6, 7, 10, 11, 13, 14]
    SGMonoB = [5, 8, 9, 12, 15]
    SGOrthP = [16, 17, 18, 19, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 47,
               48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]
    SGOrthB1 = [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68]
    SGOrthB2 = [38, 39, 40, 41]
    SGOrthF = [22, 42, 43, 69, 70]
    SGOrthI = [23, 24, 44, 45, 46, 71, 72, 73, 74]
    SGTetrP = [75, 76, 77, 78, 81, 83, 84, 85, 86, 89, 90, 91, 92, 93, 94,
               95, 96, 99, 100, 101, 102, 103, 104, 105, 106, 111, 112,
               113, 114, 115, 116, 117, 118, 123, 124, 125, 126, 127,
               128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138]
    SGTetrI = [79, 80, 82, 87, 88, 97, 98, 107, 108, 109, 110, 119, 120,
               121, 122, 139, 140, 141, 142]
    SGTrigP = [146, 148, 155, 160, 161, 166, 167]
    SGHexaP = [143, 144, 145, 147, 149, 150, 151, 152, 153, 154, 156,
               157, 158, 159, 162, 163, 164, 165, 168, 169, 170, 171,
               172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
               183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194]
    SGCubcP = [195, 198, 200, 201, 205, 207, 208, 212, 213, 215, 218,
               221, 222, 223, 224]
    SGCubcF = [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228]
    SGCubcI = [197, 199, 204, 206, 211, 214, 217, 220, 229, 230]

    # each row of prim_vec is a primitive base vector, and each column of prim_vec^-1 is a primitive base vec of BZ.
    if gid in SGTricP + SGMonoP + SGOrthP + SGTetrP + SGHexaP + SGCubcP:
        lattice_type  = 'P'
    elif gid in SGMonoB + SGOrthB1:
        lattice_type  = 'C'
    elif gid in SGOrthB2:
        lattice_type  = 'A'
    elif gid in SGOrthI + SGTetrI + SGCubcI:
        lattice_type  = 'I'
    elif gid in SGOrthF + SGCubcF:
        lattice_type  = 'F'
    elif gid in SGTrigP:
        lattice_type  = 'R'
    else:
        raise ValueError('Wrong gid!', gid)
    # return transposed prim_vec, i.e., each col a prim vector
    return lattice_type

class ssgGroup:
    """
    Spin Space Group (SSG) class for managing SSG operations and data.
    
    This class handles the complete SSG data structure including symmetry operations,
    spin rotations, SU2 matrices, and anti-unitary operations. It provides methods
    for loading SSG data, processing operations, and generating labels.
    
    Attributes:
        group_type (str): Type of group ('double' or 'single')
        ssgNum (str): SSG number identifier
        kvec (list): Translation factors
        b1, b2, b3 (list): Reciprocal lattice vectors
        pure_T (list): Pure translation operations
        superCell (list): Supercell information
        Gid (str): Group identifier
        rotC (list): Rotation matrices for all operations
        tauC (list): Translation vectors for all operations
        spin (list): Spin rotation matrices
        su2s (list): SU2 matrices for spin operations
        time_reversal (list): Time reversal factors (1 or -1)
        anti_spin (list): Anti-unitary spin operations
        anti_rotC (list): Anti-unitary rotation matrices
        anti_tau (list): Anti-unitary translation vectors
        mul_table (list): Multiplication table for operations
        factor_su2 (list): SU2 factors
        factor_trans (list): Translation factors
    """
    
    def __init__(self, ssgNum, group_type):
        """
        Initialize the SSG group.
        
        Args:
            ssgNum (str): SSG number identifier
            group_type (str): Type of group ('double' or 'single')
        """
        assert group_type == 'double' or group_type == 'single'
        self.group_type = group_type
        self.ssgNum = ssgNum  # SSG number identifier
        
        # Translation and reciprocal lattice vectors
        self.kvec = []  # Translation factors
        self.b1 = []    # First reciprocal lattice vector
        self.b2 = []    # Second reciprocal lattice vector
        self.b3 = []    # Third reciprocal lattice vector
        self.pure_T = []  # Pure translation operations
        
        # SSG structure information
        self.superCell = []  # Supercell information
        self.Gid = ''        # Group identifier
        
        # Symmetry operations
        self.rotC = []  # Rotation matrices for all operations
        self.tauC = []  # Translation vectors for all operations
        self.spin = []  # Spin rotation matrices
        self.su2s = []  # SU2 matrices for spin operations
        self.time_reversal = []  # Time reversal factors (1 or -1)
        
        # Anti-unitary operations (can be chosen arbitrarily)
        self.anti_spin = []   # Anti-unitary spin operations
        self.anti_rotC = []   # Anti-unitary rotation matrices
        self.anti_tau = []    # Anti-unitary translation vectors
        
        # Group structure
        self.mul_table = []    # Multiplication table for operations
        self.factor_su2 = []   # SU2 factors
        self.factor_trans = [] # Translation factors

    def load_ssg(self, ssg_dic):
        """
        Load SSG data from dictionary and process operations.
        
        This method loads SSG data from a dictionary and processes the symmetry
        operations, including unitary and anti-unitary operations. It handles
        different group types and constructs the appropriate operation sets.
        
        Args:
            ssg_dic (dict): Dictionary containing SSG data with keys:
                - 'URot': Spin rotation matrices
                - 'QRotC': Q-space rotation matrices
                - 'QTauC': Q-space translation vectors
                - 'superCell': Supercell information
                - 'Gid': Group identifier
                - Other SSG-specific data
        """
        def uni_or_anti(ssg_dic):
            """
            Determine if the group contains anti-unitary operations.
            
            Returns:
                list: [uni, spin, QRot, QTau] where uni=1 for unitary, -1 for anti-unitary
            """
            URot = ssg_dic['URot'][-1]
            QRot = ssg_dic['QRotC']
            QTau = ssg_dic['QTauC']
            uni = 1  # -1 if anti-unitary
            
            # Check for anti-unitary operations (determinant < 0)
            for i, spin in enumerate(URot):
                if det(spin) < 0:
                    uni = -1
                    return [uni, spin, QRot[i], QTau[i]]
            return [uni]
        self.superCell = ssg_dic['superCell']
        self.Gid = ssg_dic['Gid']
        HRot = ssg_dic['HRotC']
        HTau = ssg_dic['HTauC']
        QRot = ssg_dic['QRotC']
        QTau = ssg_dic['QTauC']
        URot = ssg_dic['URot'][-1]
        # print_matlist(URot)
        # print_matlist(QTau)
        uni = uni_or_anti(ssg_dic)
        if uni[0] == -1:
            A_spin = uni[1]
            A_rot = uni[2]
            A_tau = uni[3]
            self.anti_spin = A_spin
            self.anti_rotC = A_rot
            self.anti_tau = A_tau
        else:
            self.anti_spin = 0
            self.anti_rotC = 0
            self.anti_tau = 0
        # generate a new element list
        space_rot_uni = []
        space_tau_uni = []
        spin_uni = []
        su2_uni = []
        space_rot_tot = []
        space_tau_tot = []
        spin_tot = []
        su2_tot = []
        time_reversal = []
        for i, hrot in enumerate(HRot):
            htau = HTau[i]
            for j, spin in enumerate(URot):
                if det(spin) > 0:
                    qrot = QRot[j]
                    qtau = QTau[j]
                    # r1 t1 * r2 t2 = r1r2| r1t2 + t1
                    # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
                    rot_new = qrot @ hrot
                    tau_new = qrot @ htau + qtau
                    su2_new = SU2(spin)
                    space_rot_uni.append(rot_new)
                    space_tau_uni.append(tau_new)
                    spin_uni.append(spin)
                    su2_uni.append(su2_new)
                    time_reversal.append(1)
                    space_rot_tot.append(rot_new)
                    space_tau_tot.append(tau_new)
                    spin_tot.append(spin)
                    su2_tot.append(su2_new)

        # here the space_rot, space_tau, spin, su2 are the unitary subgroup of an antiunitary group
        # use the antiunitary element A to generate all the elements G=H+AH
        if uni[0] == -1:
            Time_R = -1j * np.array([[0, -1j], [1j, 0]])
            A_su2 = Time_R @ np.conj(SU2(-A_spin))
            for i,rot in enumerate(space_rot_uni):
                tau = space_tau_uni[i]
                spin = spin_uni[i]
                su2 =  su2_uni[i]
                space_rot_tot.append(A_rot @ rot)
                space_tau_tot.append(A_rot @ tau + A_tau)
                spin_tot.append(A_spin @ spin)
                su2_tot.append(A_su2 @ np.conj(su2))
                time_reversal.append(-1)
        self.rotC = space_rot_tot
        self.tauC = space_tau_tot
        self.spin = spin_tot
        self.su2s = su2_tot
        self.time_reversal = time_reversal

    def load_ssg_2d(self, ssg_dic):
        self.anti_spin = np.array([[1,0,0],[0,1,0],[0,0,-1]])
        self.anti_rotC = np.array([[1,0,0],[0,1,0],[0,0,1]])
        self.anti_tau = np.array([0,0,0])
        self.superCell = ssg_dic['superCell']
        self.Gid = ssg_dic['Gid']
        HRot = ssg_dic['HRotC']
        HTau = ssg_dic['HTauC']
        QRot = ssg_dic['QRotC']
        QTau = ssg_dic['QTauC']
        URot = ssg_dic['URot'][-1]

        space_rot_uni = []
        space_tau_uni = []
        spin_uni = []
        su2_uni = []
        space_rot_tot = []
        space_tau_tot = []
        spin_tot = []
        su2_tot = []
        time_reversal = []
        for i, hrot in enumerate(HRot):
            htau = HTau[i]
            for j, spin in enumerate(URot):
                if det(spin) > 0:
                    qrot = QRot[j]
                    qtau = QTau[j]
                    # r1 t1 * r2 t2 = r1r2| r1t2 + t1
                    # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
                    rot_new = qrot @ hrot
                    tau_new = qrot @ htau + qtau
                    su2_new = SU2(spin)
                else:
                    spin = np.array([[1,0,0],[0,1,0],[0,0,-1]]) @ spin
                    qrot = QRot[j]
                    qtau = QTau[j]
                    # r1 t1 * r2 t2 = r1r2| r1t2 + t1
                    # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
                    rot_new = qrot @ hrot
                    tau_new = qrot @ htau + qtau
                    su2_new = SU2(spin)
                space_rot_uni.append(rot_new)
                space_tau_uni.append(tau_new)
                spin_uni.append(spin)
                su2_uni.append(su2_new)
                time_reversal.append(1)
                space_rot_tot.append(rot_new)
                space_tau_tot.append(tau_new)
                spin_tot.append(spin)
                su2_tot.append(su2_new)

        # here the space_rot, space_tau, spin, su2 are the unitary subgroup of an antiunitary group
        # use the antiunitary element A to generate all the elements G=H+AH
        A_spin = np.array([[1,0,0],[0,1,0],[0,0,-1]])
        A_rot = np.eye(3)
        A_tau = np.array([0,0,0])
        Time_R = -1j * np.array([[0, -1j], [1j, 0]])
        A_su2 = Time_R @ np.conj(SU2(-A_spin))
        # A_su2 = SU2(A_spin)
        for i,rot in enumerate(space_rot_uni):
            tau = space_tau_uni[i]
            spin = spin_uni[i]
            su2 =  su2_uni[i]
            space_rot_tot.append(A_rot @ rot)
            space_tau_tot.append(A_rot @ tau + A_tau)
            spin_tot.append(A_spin @ spin)
            su2_tot.append(A_su2 @ np.conj(su2))
            time_reversal.append(-1)
        self.rotC = space_rot_tot
        self.tauC = space_tau_tot
        self.spin = spin_tot
        self.su2s = su2_tot
        self.time_reversal = time_reversal

    def load_ssg_1d(self, ssg_dic):
        self.anti_spin = np.array([[-1,0,0],[0,1,0],[0,0,1]])
        self.anti_rotC = np.array([[1,0,0],[0,1,0],[0,0,1]])
        self.anti_tau = np.array([0,0,0])
        self.superCell = ssg_dic['superCell']
        self.Gid = ssg_dic['Gid']
        HRot = ssg_dic['HRotC']
        HTau = ssg_dic['HTauC']
        QRot = ssg_dic['QRotC']
        QTau = ssg_dic['QTauC']
        URot = ssg_dic['URot'][-1]
        # uni = uni_or_anti(ssg_dic)
        # assert uni[0] == -1
        # A_spin = uni[1]
        # A_rot = uni[2]
        # A_tau = uni[3]
        # generate a new element list
        space_rot_uni = []
        space_tau_uni = []
        spin_uni = []
        su2_uni = []
        space_rot_tot = []
        space_tau_tot = []
        spin_tot = []
        su2_tot = []
        time_reversal = []
        for i, hrot in enumerate(HRot):
            htau = HTau[i]
            for j, spin in enumerate(URot):
                if det(spin) > 0:
                    qrot = QRot[j]
                    qtau = QTau[j]
                    # r1 t1 * r2 t2 = r1r2| r1t2 + t1
                    # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
                    rot_new = qrot @ hrot
                    tau_new = qrot @ htau + qtau
                    su2_new = SU2(spin)
                else:
                    spin = np.array([[-1,0,0],[0,1,0],[0,0,1]]) @ spin
                    qrot = QRot[j]
                    qtau = QTau[j]
                    # r1 t1 * r2 t2 = r1r2| r1t2 + t1
                    # qrot|qtau * hrot|htau = qrot*hrot| qrot*htau + qtau
                    rot_new = qrot @ hrot
                    tau_new = qrot @ htau + qtau
                    su2_new = SU2(spin)
                spin_c2z = np.array([[-1,0,0],[0,-1,0],[0,0,1]]) @ spin
                su2_c2z = SU2(spin_c2z)
                space_rot_uni.append(rot_new)
                space_tau_uni.append(tau_new)
                spin_uni.append(spin)
                su2_uni.append(su2_new)
                time_reversal.append(1)
                space_rot_tot.append(rot_new)
                space_tau_tot.append(tau_new)
                spin_tot.append(spin)
                su2_tot.append(su2_new)
                # double the unitary part by C2z
                space_rot_uni.append(rot_new)
                space_tau_uni.append(tau_new)
                spin_uni.append(spin_c2z)
                su2_uni.append(su2_c2z)
                time_reversal.append(1)
                space_rot_tot.append(rot_new)
                space_tau_tot.append(tau_new)
                spin_tot.append(spin_c2z)
                su2_tot.append(su2_c2z)

        # here the space_rot, space_tau, spin, su2 are the unitary subgroup of an antiunitary group
        # use the antiunitary element A to generate all the elements G=H+AH
        A_spin = np.array([[-1,0,0],[0,1,0],[0,0,1]])
        A_rot = np.eye(3)
        A_tau = np.array([0,0,0])
        Time_R = -1j * np.array([[0, -1j], [1j, 0]])
        A_su2 = Time_R @ np.conj(SU2(-A_spin))
        # A_su2 = SU2(A_spin)
        for i,rot in enumerate(space_rot_uni):
            tau = space_tau_uni[i]
            spin = spin_uni[i]
            su2 =  su2_uni[i]
            space_rot_tot.append(A_rot @ rot)
            space_tau_tot.append(A_rot @ tau + A_tau)
            spin_tot.append(A_spin @ spin)
            su2_tot.append(A_su2 @ np.conj(su2))
            time_reversal.append(-1)
        self.rotC = space_rot_tot
        self.tauC = space_tau_tot
        self.spin = spin_tot
        self.su2s = su2_tot
        self.time_reversal = time_reversal

    def ssg_bz(self):
        def calculate_reciprocal_lattice(pure_t):
            # cal the reciprocal basis from pure_t
            a1, a2, a3 = pure_t
            V = np.dot(a1, np.cross(a2, a3))
            b1 = (2 * np.pi * np.cross(a2, a3)) / V
            b2 = (2 * np.pi * np.cross(a3, a1)) / V
            b3 = (2 * np.pi * np.cross(a1, a2)) / V
            
            return b1, b2, b3
        tau_col = self.superCell
        gid = self.Gid
        prim_vec = identify_SG_lattice(gid)[1] # each col is a prim basis vector   
        t1 = np.array([tau_col[0][0], tau_col[1][0], tau_col[2][0]]) 
        t2 = np.array([tau_col[0][1], tau_col[1][1], tau_col[2][1]]) 
        t3 = np.array([tau_col[0][2], tau_col[1][2], tau_col[2][2]])
        pure_t = []
        for t_append in [t1, t2, t3]:
            pure_t.append(prim_vec @ t_append)
        self.pure_T = pure_t
        # pure t form the a_{i,j,k} in the basis of conventional basis of G
        b1, b2, b3 = calculate_reciprocal_lattice(pure_t)
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3


def loadSsgGroup(ssgNum, kvec, group_type, ssg_dic):
    """
    Load and initialize a complete SSG group with all operations.
    
    This function creates an SSG group and loads the appropriate operations
    based on the SSG number type (P for coplanar, L for linear, or general).
    
    Args:
        ssgNum (str): SSG number identifier
        kvec (list): Translation vector
        group_type (str): Type of group ('double' or 'single')
        ssg_dic (dict): SSG data dictionary
    
    Returns:
        ssgGroup: Initialized SSG group with all operations loaded
    """
    ssg = ssgGroup(ssgNum, group_type)
    
    # Load operations based on SSG type
    if 'P' in ssgNum:
        # Coplanar condition - 2D operations
        ssg.load_ssg_2d(ssg_dic)
    elif 'L' in ssgNum:
        # Linear condition - 1D operations
        ssg.load_ssg_1d(ssg_dic)
    else:
        # General 3D operations
        ssg.load_ssg(ssg_dic)
    
    # Initialize Brillouin zone and group structure
    ssg.ssg_bz()
    ssg.get_table()
    ssg.get_trans_factor(kvec)
    return ssg


def loadSsgOps(ssgNum, ssg_dic):
    """
    Load SSG operations for single group type.
    
    This function creates a simplified SSG group for loading operations only,
    without the full group structure.
    
    Args:
        ssgNum (str): SSG number identifier
        ssg_dic (dict): SSG data dictionary
    
    Returns:
        ssgGroup: SSG group with operations loaded
    """
    ssg = ssgGroup(ssgNum, 'single')
    ssg.load_ssg(ssg_dic)
    ssg.ssg_bz()
    return ssg


def identity_ops(gens, angle, axis, det, tau, Gid, supercell):
    """
    Check if operations are identical based on generators.
    
    This function compares symmetry operations to determine if they are
    identical based on their generators, angles, axes, determinants,
    and translation vectors.
    
    Args:
        gens (list): Generator parameters
        angle (float): Rotation angle
        axis (list): Rotation axis
        det (float): Determinant
        tau (list): Translation vector
        Gid (str): Group identifier
        supercell (list): Supercell information
    
    Returns:
        bool: True if operations are identical, False otherwise
    """
    diff1 = gens[0] - det
    axis = [axis[0], axis[1], axis[2]]
    diff2 = norm(np.cross(np.array([gens[1], gens[2], gens[3]]), np.array(axis)))
    diff3 = np.mod(angle, 360) - gens[4]
    tau2 = np.array([gens[5], gens[6], gens[7]])
    diff4 = identity_tau(tau, tau2, supercell, Gid)
    
    if norm(diff1) < 1e-4 and norm(diff2) < 1e-4 and norm(diff3) < 1e-4 and diff4[0]:
        return True
    return False

def is_int(num):
    """
    Check if a number is close enough to an integer.
    
    Args:
        num (float): Number to check
    
    Returns:
        int: 1 if close to integer, 0 otherwise
    """
    return 1 if abs(num - round(num)) < 1e-3 else 0



def Fraction(dec, tol=1e-3, max_down=1000):
    """
    Convert decimal to fraction representation.
    
    This function finds the best fractional representation of a decimal
    number within the given tolerance.
    
    Args:
        dec (float): Decimal number to convert (must be between 0 and 1)
        tol (float): Tolerance for matching (default: 1e-3)
        max_down (int): Maximum denominator to search (default: 1000)
    
    Returns:
        tuple: (numerator, denominator) of the best fraction
    
    Raises:
        ValueError: If dec is not between 0 and 1, or if no fraction is found
    """
    if dec < 0 or dec > 1:
        raise ValueError('Wrong rotation!')
    
    for down in range(max_down):
        dd = down + 1
        for up in range(dd):
            test = up / dd
            if abs(test - dec) < tol:
                return (up, dd)
    
    raise ValueError('can not find fractional rotation')


def get_op(rot_list):
    """
    Generate symmetry operation label from rotation parameters.
    
    This function converts rotation parameters (determinant, angle, axis)
    into a standard symmetry operation label (e.g., '2', '3', '4', '6', '-1', etc.).
    
    Args:
        rot_list (list): [determinant, angle, axis0, axis1, axis2]
            - determinant: 1 for proper rotation, -1 for improper
            - angle: Rotation angle in degrees
            - axis0, axis1, axis2: Rotation axis components
    
    Returns:
        str: Symmetry operation label
    """
    [det, angle, axis0, axis1, axis2] = rot_list
    
    def axis2int(axis0, axis1, axis2):
        """
        Convert axis components to integers if possible.
        
        Returns:
            tuple: Integer axis components if possible, otherwise original values
        """
        if is_int(axis0) * is_int(axis1) * is_int(axis2) == 0:
            for milti in range(6):
                mul = milti + 1
                if is_int(mul * axis0) * is_int(mul * axis1) * is_int(mul * axis2) == 1:
                    return round_vec((mul * axis0, mul * axis1, mul * axis2))
            return (axis0, axis1, axis2)
        else:
            return (axis0, axis1, axis2)
    # print(axis0, axis1, axis2)
    (axis0, axis1, axis2) = axis2int(axis0, axis1, axis2)
    if det == 1:
        detstr = ''
    if det == -1:
        detstr = '-'
    if angle < 0.1:
        if det == 1:
            return '1'
        if det == -1:
            return '-1'
    if angle < 181:
        frac = Fraction(angle/360)
        # print(angle)
        # print(frac)
        if frac == (1, 2) and det == 1:
            rtn = '2'
            return rtn
        if frac == (1, 2) and det == -1:
            rtn = 'm'
            return rtn
        if frac[0] == 1:
            rtn = detstr + str(frac[1]) 
            return rtn
        else:
            rtn = detstr + str(frac[1]) + '^{'+ str(frac[0])+'}' 
            return rtn
    else:
        frac = Fraction(angle/360)
        # [axis0, axis1, axis2] = -[axis0, axis1, axis2]
        a0 = frac[0]
        a1 = frac[1]
        a2 = a1 - a0
        if a2 == a1 + 1:
            rtn = detstr + str(a1) + '^-' 
            return rtn
        else:
            rtn = detstr + str(a1) + '^{-'+ str(a2) + '}'
            return rtn

def get_rotation1(R, return_list=True):
    if norm(R- np.eye(3)) < 1e-4:
        return [1, 0, 0, 0, 1]
    if norm(R+ np.eye(3)) < 1e-4:
        return [-1, 0, 0, 0, 1]
    det = np.linalg.det(R)
    assert np.allclose(abs(det), 1), det
    det = round(det)
    tmpR = det * R
    arg = (np.trace(tmpR) - 1) / 2
    if arg > 1:
        arg = 1
    elif arg < -1:
        arg = -1
    angle = acos(arg)
    axis = np.zeros(3)
    if abs(abs(angle) - pi) < 1e-4:
        for i in range(3):
            axis[i] = 1
            axis = axis + np.dot(tmpR, axis)
            if max(abs(axis)) > 1e-1:
                break
        assert max(abs(axis)) > 1e-1, 'can\'t find axis'
    elif abs(angle) > 1e-3:
        # standard case, see Altmann's book
        axis[0] = tmpR[2, 1] - tmpR[1, 2]
        axis[1] = tmpR[0, 2] - tmpR[2, 0]
        axis[2] = tmpR[1, 0] - tmpR[0, 1]
        axis = axis / sin(angle) / 2
    elif abs(angle) < 1e-4:
        axis[0] = 1
    # for non-orthogonal coordinates, axis may have norm>1, need to normalize
    axis = axis / np.linalg.norm(axis)
    axis = round_vec(axis)
    if axis[2] != 0:
        if axis[2] < 0:
            angle = 2*pi-angle
        axis = [axis[0]/axis[2], axis[1]/axis[2], 1]
        axis = round_vec(axis)
    elif axis[1] != 0:
        if axis[1] < 0:
            angle = 2*pi-angle
        axis = [axis[0]/axis[1], 1, 0]
        axis = round_vec(axis)
    else:
        if axis[0] < 0:
            angle = 2*pi-angle
        axis = [1, 0, 0]
    angle = angle / pi * 180
    if return_list:
        return [det, angle, axis[0], axis[1], axis[2]]
    else:
        return angle, axis, det



def get_SSG_label_list(ssgnum,ssg_list):
    # ssgnum = '229.2.1.1'
    # # read the SSG
    slash_num = [*range(10,16), *range(83,89), *range(123,143), *range(175,177), *range(191, 195)]

    ssgdic_list = []
    find_ssg = 0
    for ssg in ssg_list:
        num = ssg['ssgNum']
        if ssg['ssgNum'] == ssgnum:
            a = ssg
            find_ssg = 1
            appendix = [get_std_pg(a['eqvPg'])[1]]
            break
    if not find_ssg:
        raise ValueError('not a valid SSG number in the database, please check the input SSG number.')
    # print(a)
    Gid = int(ssgnum.split('.')[0])
    # print(gen_generators(Gid))
    label, SG_gens = gen_generators(Gid)
    Ssg = loadSsgOps(ssgnum, a)
    rotC = Ssg.rotC
    tauC = Ssg.tauC
    spinC =  Ssg.spin
    pure_T = Ssg.pure_T
    # print(Ssg.pure_T)
    
    output = [identify_lattice_type(Gid)]
    for lab, gens in zip(label, SG_gens):
        find = 0
        for rot, spin, tau in zip(rotC, spinC, tauC):
            det, angle, axis0, axis1, axis2 = get_rotation1(rot)
            # angle = angle/np.pi * 180
            axis = [axis0, axis1, axis2]
            ar = [det, axis0, axis1, axis2, angle, tau[0], tau[1], tau[2]]
            # print(ar)
            if identity_ops(gens, angle, axis, det, tau, Gid, Ssg.superCell):
                # print(get_op(get_rotation1(spin)), '||', lab)
                output.append(get_op(get_rotation1(spin))+'||'+lab)
                find = 1
        if find == 0:
            print('something wrong')

    # then the translation part
    lattice, tau_col = identify_SG_lattice(Gid)
    t1 = [1, 0, 0, 0, 0, tau_col[0][0], tau_col[1][0], tau_col[2][0]] 
    t2 = [1, 0, 0, 0, 0, tau_col[0][1], tau_col[1][1], tau_col[2][1]]
    t3 = [1, 0, 0, 0, 0, tau_col[0][2], tau_col[1][2], tau_col[2][2]]
    label = ['tau1', 'tau2', 'tau3']
    SG_gens = [t1, t2, t3]
    
    for lab, gens in zip(label, SG_gens):
        find = 0
        for rot, spin, tau in zip(rotC, spinC, tauC):
            det, angle, axis0, axis1, axis2 = get_rotation1(rot)
            # angle = angle/np.pi * 180
            axis = [axis0, axis1, axis2]
            ar = [det, axis0, axis1, axis2, angle, tau[0], tau[1], tau[2]]
            # print(ar)
            if identity_ops(gens, angle, axis, det, tau, Gid, Ssg.superCell):
                # print(get_op(get_rotation1(spin)), '||', lab)
                output.append(get_op(get_rotation1(spin))+'||'+lab)

                find = 1
        if find == 0:
            print('something wrong')

    if Gid in slash_num:
        output.insert(2,'/')
    
    if ssgnum[-1] == 'L':
        appendix.append('I')
    elif ssgnum[-1] == 'P':
        appendix.append('II')
    else:
        appendix.append('III')
    return [output,appendix]



def get_SSG_label(ssgnum,ssg_list):
    label_list,appendix = get_SSG_label_list(ssgnum,ssg_list)
    # print(label_list)
    label_str = label_list[0]+'^{'+label_list[-3].split('||')[0]\
        +','+label_list[-2].split('||')[0]\
            +','+label_list[-1].split('||')[0]+'} '
    for i in range(1, len(label_list)-3):
        if '||' in label_list[i]:
            j = label_list[i].split('||')
            if len(j[0]) > 1:
                label_str = label_str + j[1] + '^{' + j[0]+ '}' + ' '
            else:
                label_str = label_str + j[1] + '^' + j[0]+ ' '
        else:
            label_str = label_str + label_list[i] + ' '
    return label_str[:-1]+' ('+appendix[0]+'^'+appendix[1]+')'
