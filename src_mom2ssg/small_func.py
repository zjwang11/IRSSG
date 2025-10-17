import numpy as np
from numpy.linalg import norm, inv, det
from itertools import product
# from smith_form.gauss_elim import gauss_elim_np
# from find_minimal_latt import standardize_prim_basis

from math import cos, sin, acos, pi

# Try to import the C implementation, fall back to Python implementation if not available
try:
    from .smith_form.smith_form_C import smith_form
except (OSError, ImportError):
    # Fall back to pure Python implementation
    from .smith_form.smith_tool import smith_form_np
    
    def smith_form(M):
        """
        Wrapper for the pure Python Smith form implementation
        Returns: rank, M1, L, LI, R, RI
        """
        n1, n2 = M.shape
        if M.size == 0:
            R = np.eye(n2, dtype='int')
            RI = np.eye(n2, dtype='int')
            L = np.eye(n1, dtype='int')
            LI = np.eye(n1, dtype='int')
            Lmd = np.zeros([n1, n2], dtype='int')
            rank = 0
        else:
            SM, L, R, rank = smith_form_np(M)
            Lmd = SM
            LI = np.linalg.inv(L).astype(int)
            RI = np.linalg.inv(R).astype(int)
        
        return rank, Lmd, L, LI, R, RI
    
def round_vec(vec, tol=1e-6):
    """
    Round vector elements to integers if they are close to integers.
    
    Args:
        vec: Input vector
        tol: Tolerance for rounding (default: 1e-6)
        
    Returns:
        np.array: Vector with rounded elements
    """
    return np.array([np.round(v) if abs(v - np.round(v)) < tol else v for v in vec]) 

def vec_is_int(vec, tol=1e-6):
    return True if norm(vec - np.round(vec)) < tol else False

def round_num(num, tol=1e-6):
    return int(round(num)) if abs(num - round(num)) < tol else num
    
def round_mat(m, tol=1e-6):
    return np.array(np.round(m), dtype=int) if norm(m - np.round(m)) < tol else m

def find_mat(m, m_list, tol=1e-6):
    tmp = [norm(m - R) < tol for R in m_list]
    assert sum(tmp) == 1, (m, m_list)
    return tmp.index(1)

def print_mat(mat, wfile=None, print1=True):
    dim = np.shape(mat)[0]
    strg = ''
    for row in range(dim):
        strg += '[ '
        for col in range(dim):
            num = mat[row, col]
            if abs(num) < 1e-8:
                strg += '     0      '
            elif abs(num.real) < 1e-8:
                strg += ' %4.3f *I ' % (num.imag)
            elif abs(num.imag) < 1e-8:
                strg += ' %4.3f ' % (num.real)
            else:
                strg += ' (%4.3f + %4.3f *I) ' % (num.real, num.imag)
        strg += ' ]\n'
    if print1:
        print(strg, file=wfile)
    else:
        return strg


def find_monoclinic_trial_iso(Max=50):
    trial_iso = [np.array([[comb[0], 0, comb[1]],[0, 1, 0],[comb[2], 0, comb[3]]], dtype=int)
                      for comb in product(*np.tile(np.arange(-Max, Max + 1), 4).reshape((4, 2*Max+1)))
                      if comb[0] * comb[3] - comb[1] * comb[2] == 1]
    return trial_iso

def find_triclinic_trial_iso():
    # for triclinic SG 1-2, give possible isomorphisms. 
    # For SL(3, Z) matrices, Max=2 ==> 67704. Max=3 ==> N=640824
    Max = 2
   #trial_iso1 = [np.array(comb, dtype=int).reshape((3,3))
   #                  for comb in product(*np.tile(np.arange(0, Max + 1), 9).reshape((9, Max + 1)))
   #                  if det3(np.array(comb).reshape((3,3))) == 1]
    trial_iso1 = [np.array(comb, dtype=int).reshape((3,3))
                      for comb in product(*np.tile(np.arange(-Max, Max + 1), 9).reshape((9, 2*Max+1)))
                      if det3(np.array(comb).reshape((3,3))) == 1]

    # consider SL(2, Z) matrices along y,z axis, as supercells in triclinic systems are along y,z direction
    Max = 24 # 12 
    trial_iso2 = [np.array([[1, 0, 0], [0, comb[0], comb[1]], [0, comb[2], comb[3]]], dtype=int)
                      for comb in product(*np.tile(np.arange(-Max, Max + 1), 4).reshape((4, 2*Max+1)))
                      if comb[0] * comb[3] - comb[1] * comb[2] == 1]
    
    return trial_iso1 + trial_iso2

def solve_linear_inhomogeneous(M, w, mod_num=[1,1,1], tol=1e-6, allow_infinite_sol=False, test=False):
    # find v s.t. M @ v = w mod n
    # 1. compute smith form, i.e., M = L^-1 @ SM @ R^-1
    # 2. compute all possible v, s.t. SM @ v1 = L @ w, and SM @ R^-1 @ v = L @ w (v=R@v1)

    assert np.allclose(M, np.array(M).astype(int)), M
    M = np.array(M).astype(int)

    rank, SM, L, LI, R, RI = smith_form(M)
    # rank, SM, L, LI, R, RI = smith_form_sage(M)

    diag = SM.diagonal()
    w_prime = L @ w
    nop = len(diag) // 3
    mod_num_all = mod_num * (np.shape(w)[0] // 3)  # mod_num for each vi
    print('M, w:\n', M, w, '\ndiag', SM.diagonal()) if test else None 
    print('w_prime', w_prime) if test else None
    print('R', R) if test else None
    assert np.allclose(M, LI @ SM @ RI), (M, LI, SM, RI)

    if not all([round_num(wi) % mod_i == 0 for irow, wi, mod_i in zip(range(len(w_prime)), w_prime, mod_num_all) 
                if irow >= len(diag) or diag[irow] == 0]):
        # no soluation for Mv=w
       #print('mod 1 of w_prime')
       #print([round_num(wi) % mod_i  for irow, wi, mod_i in zip(range(len(w_prime)), w_prime, mod_num_all) 
       #        if irow >= len(diag) or diag[irow] == 0])
        return []

    if allow_infinite_sol:
        # if allow infinite solutions, then the diagonal elements can be zero, 
        # and are replaced by 1 to obtain finite solutions.
        diag = [1 if d == 0 else d for d in diag]
    else:
        assert all([d > 0 for d in diag]), (diag, M, w)

    vplist = [] # SM @ vp = w_prime
    for combination in product(*[range(d) for d in diag]):
        v = np.array([float(wi) / di + mod_i * ci / di for ci, di, wi, mod_i in zip(combination, diag, w_prime, mod_num_all)])
        vplist.append(v)

    vlist = [round_vec(R @ v) for v in vplist]
    assert all([norm(np.mod(round_vec(M @ v - np.array(w).astype(float).reshape(np.shape(w)[0])), mod_num_all)) < tol for v in vlist]),\
                (vlist, [round_vec(M @ v - w) for v in vlist], mod_num,
                 [round_vec(SM @ vp - w_prime) for vp in vplist])
    return vlist

def identity_vec(v1, v2, tol=1e-4):
    """
    Check if two vectors are identical within tolerance.
    
    Args:
        v1: First vector
        v2: Second vector
        tol: Tolerance for comparison (default: 1e-4)
        
    Returns:
        bool: True if vectors are identical within tolerance
    """
    diff = np.linalg.norm(np.array(v1) - np.array(v2))
    return diff < tol


def unique_list(list1):
    """
    Remove duplicate elements from a list while preserving order.
    
    Args:
        list1: Input list (can contain any hashable elements)
        
    Returns:
        list: List with unique elements in original order
    """
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list


def unique_matrix(matrix_list, tol=1e-3):
    unique_list = [matrix_list[0]]  # add the first
    for mat in matrix_list[1:]:
        is_unique = True
        for uni in unique_list:
            if np.linalg.norm(mat - uni) < tol:
                is_unique = False
                break
        if is_unique:
            unique_list.append(mat)
    return unique_list


def sort_array(target):
    """
    Sort a 2D array by rows with numerical precision control.
    
    Args:
        target: 2D numpy array to be sorted
        
    Returns:
        np.array: Sorted 2D array
    """
    size_1, size_2 = target.shape
    
    # Round to 5 decimal places for numerical stability
    for m in range(size_1):
        for n in range(size_2):
            target[m, n] = round(target[m, n]*1e+5)/ 1e+5
    
    # Convert to list, sort, and convert back to array
    sort_list = target.tolist()
    sort_list.sort()
    return np.array(sort_list)

def count_identity_pairs(int_matrices, float_matrices, tol=1e-3):
    """
    Count the number of pairs where both matrices are identity matrices.
    
    Args:
        int_matrices: List of 3x3 integer matrices
        float_matrices: List of 3x3 float matrices  
        tol: Tolerance for float matrix comparison
        
    Returns:
        int: Number of pairs where both matrices are identity matrices
    """
    identity_count = 0
    identity_3x3 = np.eye(3)
    
    for i in range(len(int_matrices)):
        int_mat = np.array(int_matrices[i])
        float_mat = np.array(float_matrices[i])
        
        # Check if integer matrix is identity
        int_is_identity = identity_vec(int_mat.flatten(), identity_3x3.flatten(), tol)
        
        # Check if float matrix is identity  
        float_is_identity = identity_vec(float_mat.flatten(), identity_3x3.flatten(), tol)
        
        if int_is_identity and float_is_identity:
            identity_count += 1
            
    return identity_count

def change_tau(tau):
    """
    Normalize translation vector to fractional coordinates [0,1).
    
    Args:
        tau: Translation vector (3-element array)
        
    Returns:
        np.array: Normalized translation vector with values in [0,1)
    """
    tau_new = [(round(tau[0]*100)/100)%1, (round(tau[1]*100)/100)%1, (round(tau[2]*100)/100)%1]
    return np.array(tau_new) 

def findDimension(cell,tolm=1e-4): # collinear or coplanar or no-coplanar, return 3, 2, 1
    mag = cell[3].copy()
    rank = np.linalg.matrix_rank(mag,tol=tolm)
    return rank

def is_integer_matrix(A, tol=1e-4):

    A = np.asarray(A, dtype=float)
    return np.all(np.abs(A - np.rint(A)) <= tol)


def l1_order_in_l2(L1,
                   L2,
                   atol: float = 1e-4,
                   rtol: float = 1e-4):

    if len(L1) != len(L2):
        return False

    def equal_mod_int(a, b, tol=1e-8) -> bool:
        a = np.asarray(a); b = np.asarray(b)
        if a.shape != b.shape:
            return False
        r = (a - b) - np.rint(a - b)
        return np.all(np.abs(r) <= tol)

    def triple_match(t1, t2) -> bool:
        for i, (a, b) in enumerate(zip(t1, t2)):
            a = np.asarray(a); b = np.asarray(b)
            if a.shape != b.shape:
                return False
            if i == 1:
                if not equal_mod_int(a, b, atol):
                    return False
            else:
                if not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
                    return False
        return True

    used = [False] * len(L2)
    order = [-1] * len(L1)

    for i, t1 in enumerate(L1):
        found = False
        for j, t2 in enumerate(L2):
            if not used[j] and triple_match(t1, t2):
                used[j] = True
                order[i] = j
                found = True
                break
        if not found:
            miss_L1 = [k for k, idx in enumerate(order) if idx == -1] + [i]
            miss_L2 = [k for k, u in enumerate(used) if not u]
            return False

    return order

        
def _rot_about_axis(axis, angle):
    axis = np.asarray(axis, float)
    n = np.linalg.norm(axis)
    if n < 1e-5:  
        return np.eye(3)
    x,y,z = axis / n
    c,s = np.cos(angle), np.sin(angle)
    t = 1 - c
    return np.array([[t*x*x + c,     t*x*y - s*z,  t*x*z + s*y],
                     [t*x*y + s*z,   t*y*y + c,    t*y*z - s*x],
                     [t*x*z - s*y,   t*y*z + s*x,  t*z*z + c]], float)


def dedup_real_3x3(L, tol=1e-4):
    A = np.array([np.asarray(x, float).reshape(3,3) for x in L])   
    N = len(A)
    Q = np.round(A / tol).astype(np.int64).reshape(N, -1)           
    key_to_idx, uniq_indices = {}, []
    inv_map = np.empty(N, dtype=int)
    for i in range(N):
        key = tuple(Q[i])                                         
        if key in key_to_idx:
            inv_map[i] = key_to_idx[key]
        else:
            key_to_idx[key] = len(uniq_indices)
            inv_map[i] = key_to_idx[key]
            uniq_indices.append(i)
    uniq = A[uniq_indices]
    return uniq

def wrap01(f):
    f = np.asarray(f, float)
    return f - np.floor(f)

def dedup_vectors(V, tol=None):
    A = np.asarray(V, dtype=float)
    if A.ndim != 2:  
        A = np.vstack([np.asarray(x, float).ravel() for x in V])

    K = A if tol is None else np.round(A / tol).astype(np.int64)
    Kc = np.ascontiguousarray(K)

    seen = set()
    keep_idx = []
    for i in range(Kc.shape[0]):
        key = Kc[i].tobytes()  
        if key not in seen:
            seen.add(key)
            keep_idx.append(i)
    return A[keep_idx]




def get_rotation(R, return_list=False):
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



def SU2(so3):
    """
    Convert SO(3) rotation matrix to SU(2) matrix.
    
    This function converts a 3D rotation matrix to its corresponding 2x2 SU(2) 
    representation using the standard mapping from SO(3) to SU(2).
    
    Args:
        so3: 3x3 SO(3) rotation matrix
        
    Returns:
        np.array: 2x2 SU(2) matrix (complex)
    """
    # Pauli matrices
    sigma0 = np.array([[1, 0], [0, 1]], dtype=complex)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
    
    # Extract rotation angle and axis
    angle, axis, det = get_rotation(so3)
    
    # Construct SU(2) matrix using Rodrigues' formula
    su2 = (cos(angle / 2) * sigma0 - 
           1j * sin(angle / 2) * (axis[0] * sigma1 + axis[1] * sigma2 + axis[2] * sigma3))
    
    return su2


def sort_rot(rot_list, tau_list):
    # sort rot according to det, angle, axis.
    nop = len(rot_list)
    det_angle_axis_list = []
    for rot in rot_list:
        angle, axis, det = get_rotation(rot)
        det_angle_axis_list.append([det, angle, axis[0], axis[1], axis[2]])

    det_angle_axis_list_orig = det_angle_axis_list.copy()
    det_angle_axis_list.sort(key=lambda x:(-x[0],x[1],x[2],x[3],x[4]))
    #print('sorted:',det_angle_axis_list)
    sorted_rot, sorted_tau, sort_order, sort_order_inv = [], [], np.zeros(nop, dtype=int), []
    for ith, item in enumerate(det_angle_axis_list):
        index = det_angle_axis_list_orig.index(item)
        sorted_rot.append(rot_list[index])
        sorted_tau.append(tau_list[index])
        sort_order[index] = ith  # i-th op in orig ops is sorted to sort_order[i]-th op
        sort_order_inv.append(index)  # i-th op in sorted ops is sort_order_inv[i]-th orig op

    return sorted_rot, sorted_tau, sort_order 


def _plane_normal_from_vectors(V, tol=1e-6):
    V = np.asarray(V, float)
    if V.ndim == 1: V = V[None, :]
    nz = [v for v in V if np.linalg.norm(v) > tol]
    v1 = nz[0]
    for v2 in nz[1:]:
        c = np.cross(v1, v2)
        
        if np.linalg.norm(c) > tol:
            return c / np.linalg.norm(c)
    raise ValueError("Error: collinear!!!")

def line_normal_vector(moments, tolm=1e-4):
    v = np.asarray(moments, dtype=float)
    v = v[np.linalg.norm(v, axis=1) > tolm]
    if v.size == 0:
        raise ValueError("Zero magnetic moment")
    _, s, vh = np.linalg.svd(v)
    if len(s)>1:
        if (s[1] / s[0]) > 1e-2:      
            raise ValueError("Not collinear")
    d = vh[0]                     
    d /= np.linalg.norm(d)      

    return d


def generate_normal_vector(magnetic_moments, tolm=1e-4):
    nonzero_moments = [m for m in magnetic_moments if not np.linalg.norm(m - [0,0,0]) < tolm]

    moments_matrix = np.array(nonzero_moments)
    
    _, _, vh = np.linalg.svd(moments_matrix)
    
    normal_vector = vh[-1]
    
    normal_vector = normal_vector / np.linalg.norm(normal_vector)
    
    return normal_vector

def align_axis_to_z(a):
    a = np.asarray(a, float)
    a = a / np.linalg.norm(a)
    h = np.array([1.0, 0.0, 0.0]) if abs(np.dot(a, [1,0,0])) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = h - np.dot(h, a) * a
    x /= np.linalg.norm(x)
    y = np.cross(a, x)
    M = np.column_stack([x, y, a])  # columns = x', y', z'
    return M  # use R_new = M.T @ R @ M


def orthonormal_basis_from_vector(n):
    n = np.asarray(n, dtype=float).ravel()
    n_norm = np.linalg.norm(n)
    n = n / n_norm                       

    idx = np.argmin(np.abs(n))
    e = np.eye(3)[idx]                   

    # 2. u = n × e
    u = np.cross(n, e)
    u /= np.linalg.norm(u)

    # 3. v = n × u
    v = np.cross(n, u)
    v /= np.linalg.norm(v)

    return u, v


def _make_basis_lock_z(xp, zp, tol=1e-3):
    """
    Build an orthonormal right-handed basis with a locked z':
      - z' is kept as given (only normalized).
      - x' is orthogonalized to z' (Gram–Schmidt).
      - y' = z' × x'.
    Returns:
        C : (3,3) whose columns are unit vectors [x', y', z'] expressed in the OLD basis.
             det(C) = +1 (right-handed).
    """
    xp = np.asarray(xp, float).reshape(3)
    zp = np.asarray(zp, float).reshape(3)

    # Normalize locked z'
    nz = np.linalg.norm(zp)
    if nz < tol:
        raise ValueError("z' is too small in magnitude.")
    z = zp / nz

    # Orthogonalize x' to z'
    x_ortho = xp - z * np.dot(xp, z)
    nx = np.linalg.norm(x_ortho)
    if nx < tol:
        # Fallback if xp is nearly colinear with z': pick a seed least aligned with z
        e = np.eye(3)
        seed = e[np.argmin(np.abs(e @ z))]
        x_ortho = seed - z * np.dot(seed, z)
        nx = np.linalg.norm(x_ortho)
        if nx < tol:
            raise ValueError("Degenerate input: cannot construct x' orthogonal to z'.")

    x = x_ortho / nx

    # Right-handed basis: y' = z' × x'
    y = np.cross(z, x)
    ny = np.linalg.norm(y)
    if ny < tol:
        raise ValueError("Degenerate basis: y' has near-zero norm (numerical issue).")
    y /= ny

    C = np.column_stack([x, y, z])

    # Ensure right-handedness; if det < 0, flip x' (do not touch locked z')
    if np.linalg.det(C) < 0:
        x = -x
        y = np.cross(z, x)
        y /= np.linalg.norm(y)
        C = np.column_stack([x, y, z])

    return C


def change_basis_O3(Rs, xprime, zprime, tol=1e-3):
    """
    Convert one or many O(3) matrices from OLD-basis representation to the NEW basis
    defined by x', z' (given in OLD coordinates, not necessarily orthonormal).
    Formula: R_new = C^T R_old C, where columns of C are [x', y', z'] in OLD coords.
    Args:
        Rs     : array of shape (3,3) or (...,3,3)
        xprime : array-like shape (3,), x' in OLD coords
        zprime : array-like shape (3,), z' in OLD coords (locked direction)
    Returns:
        Array with same shape as Rs, now expressed in the NEW basis.
    """
    C = _make_basis_lock_z(xprime, zprime, tol=tol)
    A = np.asarray(Rs, float)
    if A.shape == (3, 3):
        return C.T @ A @ C
    if A.ndim >= 2 and A.shape[-2:] == (3, 3):
        return np.einsum('ij,...jk,kl->...il', C.T, A, C)
    raise ValueError("Rs must have shape (3,3) or (...,3,3).")


def transform_vectors_to_new_coords(vecs, xprime, zprime, tol=1e-3):
    """
    Transform vectors from OLD to NEW coordinates: v_new = C^T v_old.
    Args:
        vecs   : array shape (3,) or (...,3)
        xprime : x' in OLD coords
        zprime : z' in OLD coords (locked direction)
    Returns:
        Array of same shape as vecs, expressed in the NEW basis.
    """
    C = _make_basis_lock_z(xprime, zprime, tol=tol)
    V = np.asarray(vecs, float)
    if V.shape == (3,):
        return C.T @ V
    if V.ndim >= 1 and V.shape[-1] == 3:
        return np.einsum('ij,...j->...i', C.T, V)
    raise ValueError("vecs must have shape (3,) or (...,3).")
