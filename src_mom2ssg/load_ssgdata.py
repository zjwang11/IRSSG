#!/usr/bin/env python3
"""
SSG Data Loading Module

This module provides functions for loading and processing Spin Space Group (SSG) data
from pickle files. It handles the loading of SSG lists and individual SSG data,
as well as the construction of standardized symmetry operations.

Main functionality:
- Load complete SSG data lists from pickle files
- Load individual SSG data by SSG number
- Construct standardized symmetry operations from SSG data
- Handle spin-dependent and spin-independent operations

Author: MOM2SSG Team
"""

import os
import pickle
import numpy as np
def load_ssg_list():
    """
    Load the complete list of SSG data from the pickle file.
    
    This function loads all available Spin Space Group data from the identify.pkl file
    located in the ssg_data directory. The pickle file contains a list of dictionaries,
    each representing a different SSG with its symmetry operations and properties.
    
    Returns:
        list: List of dictionaries containing SSG data. Each dictionary contains:
            - 'ssgNum': SSG number (e.g., '1.1', '1.2L', '1.3P')
            - 'HRotC': H-space rotation matrices
            - 'HTauC': H-space translation vectors
            - 'QRotC': Q-space rotation matrices
            - 'QTauC': Q-space translation vectors
            - 'URot': Spin rotation matrices
            - Other SSG-specific properties
    
    Raises:
        FileNotFoundError: If the identify.pkl file is not found
        pickle.PickleError: If the pickle file is corrupted or invalid
    """
    # Get the directory containing this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct path to the SSG data file
    pkl_file_path = os.path.join(script_dir, 'ssg_data', 'identify.pkl')
    
    # Load the SSG data from pickle file
    with open(pkl_file_path, 'rb') as f:
        ssg_list = pickle.load(f)
    
    return ssg_list

def load_one_ssg(ssgnum):
    """
    Load a specific SSG data by its SSG number.
    
    This function searches through the complete SSG data list to find and return
    the data for a specific Spin Space Group identified by its SSG number.
    
    Args:
        ssgnum (str): SSG number to load (e.g., '1.1', '1.2L', '1.3P')
    
    Returns:
        dict: Dictionary containing the SSG data for the specified SSG number.
              The dictionary contains the same structure as described in load_ssg_list().
    
    Raises:
        FileNotFoundError: If the identify.pkl file is not found
        ValueError: If the specified SSG number is not found in the data
        pickle.PickleError: If the pickle file is corrupted or invalid
    
    Example:
        >>> ssg_data = load_one_ssg('1.1')
        >>> print(ssg_data['ssgNum'])
        '1.1'
    """
    # Get the directory containing this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct path to the SSG data file
    file = os.path.join(script_dir, 'ssg_data', 'identify.pkl')
    
    # Check if the data file exists
    if not os.path.exists(file):
        raise FileNotFoundError(f"SSG data file not found: {file}")
    
    # Load the complete SSG data list
    with open(file, 'rb') as f:
        ssg_list = pickle.load(f)
    
    # Search for the specific SSG number
    for ssg in ssg_list:
        if ssg['ssgNum'] == ssgnum:
            return ssg
    
    # Raise error if SSG number not found
    raise ValueError(f'Cannot load SSG number {ssgnum} from the data')


def construct_std_operations(ssg_dic, spin_id=None):
    """
    Construct standardized symmetry operations from SSG data.
    
    This function processes the raw SSG data to construct standardized symmetry
    operations that can be used for crystal structure analysis. It combines
    H-space and Q-space operations to generate the complete set of symmetry
    operations for the given SSG.
    
    The function handles two modes:
    1. spin_id specified: Returns operations for a specific spin configuration
    2. spin_id=None: Returns operations for all possible spin configurations
    
    Args:
        ssg_dic (dict): SSG data dictionary containing:
            - 'HRotC': H-space rotation matrices (list of 3x3 arrays)
            - 'HTauC': H-space translation vectors (list of 3x1 arrays)
            - 'QRotC': Q-space rotation matrices (list of 3x3 arrays)
            - 'QTauC': Q-space translation vectors (list of 3x1 arrays)
            - 'URot': Spin rotation matrices (list of lists of 3x3 arrays)
        spin_id (int, optional): Index of specific spin configuration to use.
                                If None, processes all spin configurations.
    
    Returns:
        tuple: (rots, taus, spins) containing:
            - rots (list): List of 3x3 rotation matrices for all operations
            - taus (list): List of 3x1 translation vectors for all operations
            - spins (list): List of spin rotation matrices:
                - If spin_id is specified: List of 3x3 spin rotation matrices
                - If spin_id is None: List of lists, each containing spin rotations
                  for one spin configuration
    
    Note:
        The function applies the transformation: (R, t) = (Q*H, Q*τ_H + τ_Q)
        where H and Q are H-space and Q-space operations respectively.
        The operations are then standardized using the inverse transformation matrix S.
    """
    # Extract H-space operations (crystal structure symmetries)
    HRot = ssg_dic['HRotC']  # H-space rotation matrices
    HTau = ssg_dic['HTauC']  # H-space translation vectors
    
    # Extract Q-space operations (magnetic symmetries)
    QRot = ssg_dic['QRotC']  # Q-space rotation matrices
    QTau = ssg_dic['QTauC']  # Q-space translation vectors
    
    # Handle spin operations based on spin_id parameter
    if spin_id is not None:
        # Use specific spin configuration
        URot = ssg_dic['URot'][spin_id]
    else:
        # Use all spin configurations
        URots = ssg_dic['URot']
    
    # Standardization transformation matrix (identity for now)
    S = np.eye(3)
    invS = np.linalg.inv(S)
    
    # Initialize output lists
    rots = []  # Rotation matrices
    taus = []  # Translation vectors
    spins = []  # Spin rotation matrices
    
    # Initialize spin lists for all configurations if spin_id is None
    if spin_id is None:
        for URot in URots:
            spins.append([])
    
    # Construct operations by combining H-space and Q-space operations
    for i, hrot in enumerate(HRot):
        htau = HTau[i]
        
        if spin_id is not None:
            # Process specific spin configuration
            for j, spin in enumerate(URot):
                qrot = QRot[j]
                qtau = QTau[j]
                
                # Combine operations: (R, t) = (Q*H, Q*τ_H + τ_Q)
                rot_new = qrot @ hrot
                tau_new = qrot @ htau + qtau
                
                # Store spin rotation and standardized operations
                spins.append(spin)
                rots.append(invS @ rot_new @ S)
                taus.append(invS @ tau_new)
        else:
            # Process all spin configurations
            for j, _ in enumerate(URots[0]):
                qrot = QRot[j]
                qtau = QTau[j]
                
                # Combine operations: (R, t) = (Q*H, Q*τ_H + τ_Q)
                rot_new = qrot @ hrot
                tau_new = qrot @ htau + qtau
                
                # Store standardized operations
                rots.append(invS @ rot_new @ S)
                taus.append(invS @ tau_new)
                
                # Store spin rotations for all configurations
                for k in range(len(URots)):
                    spins[k].append(URots[k][j])
    
    return rots, taus, spins