import numpy as np
from scipy.spatial.transform import Rotation as Rot
from math import cos, sin, acos, pi
from numpy.linalg import norm, inv, det
from .small_func import findDimension,round_vec,line_normal_vector,generate_normal_vector,orthonormal_basis_from_vector

def get_rotation(R):
    det = np.linalg.det(R)
    tmpR = det * R
    arg = (np.trace(tmpR) - 1) / 2
    if arg > 1:
        arg = 1
    elif arg < -1:
        arg = -1
    angle = acos(arg)
    axis = np.zeros((3, 1))
    if abs(abs(angle) - pi) < 1e-4:
        for i in range(3):
            axis[i] = 1
            axis = axis + np.dot(tmpR, axis)
            if max(abs(axis)) > 1e-1:
                break
        assert max(abs(axis)) > 1e-1, 'can\'t find axis'
        axis = axis / np.linalg.norm(axis)
    elif abs(angle) > 1e-3:
        # standard case, see Altmann's book
        axis[0] = tmpR[2, 1] - tmpR[1, 2]
        axis[1] = tmpR[0, 2] - tmpR[2, 0]
        axis[2] = tmpR[1, 0] - tmpR[0, 1]
        axis = axis / sin(angle) / 2
    elif abs(angle) < 1e-4:
        axis[0] = 1

    return angle, axis, det

def SU2(so3):  
    # SU2 matrix calculated from SO3 matrix, may be different from su2s read from Bilbao
    sigma0 = np.array([[1, 0], [0, 1]], dtype=complex)
    sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)

    angle, axis, det = get_rotation(so3)
    # print(angle, axis)
    # print(so3)
    su2 = cos(angle / 2) * sigma0 - 1j * sin(angle / 2) * (axis[0] * sigma1 + axis[1] * sigma2 + axis[2] * sigma3)
    # selsu2c.append(su2)
    return su2

def axis_angle_to_so3_scipy(axis, angle):
    axis = np.asarray(axis, dtype=float)
    rotvec = axis / np.linalg.norm(axis) * angle 
    R = Rot.from_rotvec(rotvec).as_matrix()      
    return R


            
def get_SU2_list(matlist):
    su2_list = []
    for i, mat in enumerate(matlist):
        mat_ = SU2(np.array(mat))
        su2_list.append(mat_)
    return su2_list


def judge_time_reversal(matlist):
    T_list = []
    for mat in matlist:
        det = np.linalg.det(mat)
        if det<0:
            T_list.append(-1)
        else:
            T_list.append(1)
    return T_list



def generate_irssg_in(spg, ssgnum, msgnum, cell, mag, operations, msg_operations, tolm=1e-4):
    dim_mag = findDimension(cell,tolm=tolm)
    
    wbfile = open('ssg.data','wb')
    # print(type(spg))
    np.array([spg],dtype=np.int32).tofile(wbfile)
    b   = ssgnum.encode('utf-8')
    np.array([len(b)], dtype=np.int32).tofile(wbfile)
    wbfile.write(b)
    
    spin_only_T_list = []
    spin_only_O3_list = []
    spin_only_SU2_list = []
    
    if dim_mag == 3:
        np.array([1],dtype=np.int32).tofile(wbfile)
        T_list = [1]
        np.array([np.eye(3)], dtype=np.float64).tofile(wbfile)
        np.array([np.eye(2,dtype=np.complex128)], dtype=np.complex128).tofile(wbfile)
        np.array(T_list,dtype=np.int32).tofile(wbfile)
        spin_only_O3_list = [np.eye(3)]
        spin_only_SU2_list = [np.eye(2,dtype=np.complex128)]
    elif dim_mag == 2:
        np.array([2],dtype=np.int32).tofile(wbfile)
        T_list = [1,-1]
        normal_vector = generate_normal_vector(mag,tolm=tolm)
        spin_only_Mz_O3 = axis_angle_to_so3_scipy(normal_vector,np.pi) * (-1)
        spin_only_Mz_SU2 = SU2(axis_angle_to_so3_scipy(normal_vector,np.pi))
        np.array([np.eye(3,dtype=np.float64),spin_only_Mz_O3], dtype=np.float64).tofile(wbfile)
        np.array([np.eye(2,dtype=np.complex128),spin_only_Mz_SU2], dtype=np.complex128).tofile(wbfile)
        np.array(T_list,dtype=np.int32).tofile(wbfile)
        
        spin_only_O3_list = [np.eye(3,dtype=np.float64),spin_only_Mz_O3]
        spin_only_SU2_list = [np.eye(2,dtype=np.complex128),spin_only_Mz_SU2]
        
    else:
        np.array([4],dtype=np.int32).tofile(wbfile)
        T_list = [1,1,-1,-1]

        normal_vector = line_normal_vector(mag,tolm=tolm)

        spin_only_C2z_O3 = axis_angle_to_so3_scipy(normal_vector,np.pi)
        spin_only_C2z_SU2 = SU2(axis_angle_to_so3_scipy(normal_vector,np.pi))

        v1, v2 = orthonormal_basis_from_vector(normal_vector)
        spin_only_Mv1_O3 = axis_angle_to_so3_scipy(v1,np.pi)*(-1)
        spin_only_Mv1_SU2 =  SU2(axis_angle_to_so3_scipy(v1,np.pi))
        spin_only_Mv2_O3 = axis_angle_to_so3_scipy(v2,np.pi)*(-1)
        spin_only_Mv2_SU2 =  SU2(axis_angle_to_so3_scipy(v2,np.pi))

        np.array([np.eye(3,dtype=np.float64),spin_only_C2z_O3,spin_only_Mv1_O3,spin_only_Mv2_O3], dtype=np.float64).tofile(wbfile)
        np.array([np.eye(2,dtype=np.complex128),spin_only_C2z_SU2,spin_only_Mv1_SU2,spin_only_Mv2_SU2], dtype=np.complex128).tofile(wbfile)
        np.array(T_list,dtype=np.int32).tofile(wbfile)
        
        spin_only_O3_list = [np.eye(3,dtype=np.float64),spin_only_C2z_O3,spin_only_Mv1_O3,spin_only_Mv2_O3]
        spin_only_SU2_list = [np.eye(2,dtype=np.complex128),spin_only_C2z_SU2,spin_only_Mv1_SU2,spin_only_Mv2_SU2]
        
    spin_only_T_list = T_list.copy()

    
    
    
    np.array([len(operations['spin'])], dtype=np.int32).tofile(wbfile)
    np.array(operations['spin']).tofile(wbfile)        
    np.array(get_SU2_list(operations['spin']),dtype=np.complex128).tofile(wbfile)
    T_list = judge_time_reversal(operations['spin'])
    np.array(T_list,dtype=np.int32).tofile(wbfile)
    np.rint(np.array(operations['RotC'])).astype(np.int32).tofile(wbfile)
    np.array(operations['TauC']).tofile(wbfile)
    wbfile.close()
    
    T_list_all = []
    O3_list_all = []
    SU2_list_all = []
    for i in range(len(spin_only_T_list)):
        for j in range(len(operations['spin'])):
            T_list_all.append(spin_only_T_list[i]*T_list[j] < 0)
            O3_list_all.append(operations['RotC'][j])
            SU2_list_all.append(spin_only_SU2_list[i]@get_SU2_list(operations['spin'])[j])
    
    # for wannier symmetrization
    np.save('ssgop_wansym.npy',{'T_list':T_list_all,'O3_list':O3_list_all,'SU2_list':SU2_list_all},allow_pickle=True)
    
    # for msg
    wbfile = open('msg.data','wb')
    np.array([spg],dtype=np.int32).tofile(wbfile)
    b   = msgnum.encode('utf-8')
    np.array([len(b)], dtype=np.int32).tofile(wbfile)
    wbfile.write(b)
    
    np.array([1],dtype=np.int32).tofile(wbfile)
    np.array([np.eye(3)], dtype=np.float64).tofile(wbfile)
    np.array([np.eye(2,dtype=np.complex128)], dtype=np.complex128).tofile(wbfile)
    np.array([1],dtype=np.int32).tofile(wbfile)
    
    np.array([len(msg_operations['spin'])], dtype=np.int32).tofile(wbfile)
    np.array(msg_operations['spin']).tofile(wbfile)        
    np.array(get_SU2_list(msg_operations['spin']),dtype=np.complex128).tofile(wbfile)
    T_list = judge_time_reversal(msg_operations['spin'])
    np.array(T_list,dtype=np.int32).tofile(wbfile)
    np.rint(np.array(msg_operations['RotC'])).astype(np.int32).tofile(wbfile)
    np.array(msg_operations['TauC']).tofile(wbfile)
    
    wbfile.close()