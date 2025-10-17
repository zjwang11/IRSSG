import numpy as np
import spglib



def identify_SG_lattice(gid):
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
        latt = 'P'
        prim_vec = np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    elif gid in SGMonoB + SGOrthB1:
        latt = 'B'
        prim_vec = np.array([[1. / 2, -1. / 2, 0.], [1. / 2, 1. / 2, 0.], [0., 0., 1.]])
    elif gid in SGOrthB2:
        latt = 'B2'
        prim_vec = np.array([[1., 0., 0.], [0., 1. / 2, 1. / 2], [0., -1. / 2, 1. / 2]])
    elif gid in SGOrthI + SGTetrI + SGCubcI:
        latt = 'I'
        prim_vec = np.array([[-1. / 2, 1. / 2, 1. / 2], [1. / 2, -1. / 2, 1. / 2], [1. / 2, 1. / 2, -1. / 2]])
    elif gid in SGOrthF + SGCubcF:
        latt = 'F'
        prim_vec = np.array([[0., 1. / 2, 1. / 2], [1. / 2, 0., 1. / 2], [1. / 2, 1. / 2, 0.]])
    elif gid in SGTrigP:
        latt = 'R'
        prim_vec = np.array([[2. / 3, 1. / 3, 1. / 3], [-1. / 3, 1. / 3, 1. / 3], [-1. / 3, -2. / 3, 1. / 3]])
    else:
        raise ValueError('Wrong gid!', gid)
    # return transposed prim_vec, i.e., each col a prim vector
    return latt, prim_vec.T


def generate_coord(rot_list, tau_list):
    # generate coord list using {rot|tau}, starting from a given generic position
    gen_pts = [np.array([0.1722, 0.6933, 0.9344]),
                np.array([0.8399, 0.5677, 0.0655]),
                np.array([0.1234, 0.4567, 0.0876]),
                np.array([0.2468, 0.5721, 0.7834])]
    pos_list = []
    for R, t in zip(rot_list, tau_list):
        for p in gen_pts:
           #pos_list.append(latt_home(R @ p + t))
            pos_list.append(R @ p + t)

    return pos_list, len(gen_pts)



def gen_generators(Gid):
    def comb(*ops):
        return [sum(values) for values in zip(*ops)]
    # ops = [det, axis_x, axis_y, axis_z, angel, taux, tauy, tauz]
    E = [1, 0, 0, 1, 0, 0, 0, 0]
    P = [-1, 0, 0, 1, 0, 0, 0, 0]
    C2z = [1, 0, 0, 1, 180, 0, 0, 0]
    C2x = [1, 1, 0, 0, 180, 0, 0, 0]
    C2y = [1, 0, 1, 0, 180, 0, 0, 0]
    C2xy = [1, 1, 1, 0, 180, 0, 0, 0]
    C2Tri = [1, 2, 1, 0, 180, 0, 0, 0]
    MTri = [-1, 2, 1, 0, 180, 0, 0, 0]
    C3z = [1, 0, 0, 1, 120, 0, 0, 0]
    PC3z = [-1, 0, 0, 1, 120, 0, 0, 0]
    C3h = [1, 1, 1, 1, 120, 0, 0, 0]
    PC3h = [-1, 1, 1, 1, 120, 0, 0, 0]
    C4z = [1, 0, 0, 1, 90, 0, 0, 0]
    C4x = [1, 1, 0, 0, 90, 0, 0, 0]
    PC4z = [-1, 0, 0, 1, 90, 0, 0, 0]
    PC4x = [-1, 1, 0, 0, 90, 0, 0, 0]
    C6z = [1, 0, 0, 1, 60, 0, 0, 0]
    PC6z = [-1, 0, 0, 1, 60, 0, 0, 0]
    Mz = [-1, 0, 0, 1, 180, 0, 0, 0]
    Mx = [-1, 1, 0, 0, 180, 0, 0, 0]
    My = [-1, 0, 1, 0, 180, 0, 0, 0]
    Mxy = [-1, 1, 1, 0, 180, 0, 0, 0]
    taux = [0,0,0,0,0,0.5,0,0]
    tauy = [0,0,0,0,0,0,0.5,0]
    tauz = [0,0,0,0,0,0,0,0.5]
    Fz = [0,0,0,0,0,0,0,0.25]
    Fxy = [0,0,0,0,0,0.25,0.25,0]
    Fxz = [0,0,0,0,0,0.25,0,0.25]
    Fyz = [0,0,0,0,0,0,0.25,0.25]
    Fxyz = [0,0,0,0,0,0.25,0.25,0.25]
    tau6 = [0,0,0,0,0,0,0,1/6]
    tau3 = [0,0,0,0,0,0,0,1/3]
    if Gid == 1:
        label = ['1']
        ops = [E]
    if Gid == 2:
        label = ['-1']
        ops = [P]
    if Gid == 3:
        label = ['2']
        ops = [C2y]
    if Gid == 4:
        label = ['2_1']
        ops = [[1, 0, 1, 0, 180, 0, 0.5, 0]]
    if Gid == 5:
        label = ['2']
        ops = [C2y]
    if Gid == 6:
        label = ['m']
        ops = [My]
    if Gid == 7:
        label = ['c']
        ops = [[-1, 0, 1, 0, 180, 0, 0, 0.5]]
    if Gid == 8:
        label = ['m']
        ops = [My]
    if Gid == 9:
        label = ['c']
        ops = [[-1, 0, 1, 0, 180, 0, 0, 0.5]]
    if Gid == 10:
        label = ['2', 'm']
        ops = [C2y, My]
    if Gid == 11:
        label = ['2_1', 'm']
        ops = [[1, 0, 1, 0, 180, 0, 0.5, 0],[-1, 0, 1, 0, 180, 0, 0.5, 0]]
    if Gid == 12:
        label = ['2', 'm']
        ops = [C2y, My]
    if Gid == 13:
        label = ['2', 'c']
        ops = [[1, 0, 1, 0, 180, 0, 0, 0.5],[-1, 0, 1, 0, 180, 0, 0, 0.5]]
    if Gid == 14:
        label = ['2_1', 'c']
        ops = [[1, 0, 1, 0, 180, 0, 0.5, 0.5],[-1, 0, 1, 0, 180, 0, 0.5, 0.5]]
    if Gid == 15:
        label = ['2', 'c']
        ops = [[1, 0, 1, 0, 180, 0, 0, 0.5],[-1, 0, 1, 0, 180, 0, 0, 0.5]]
    if Gid == 16:
        label = ['2', '2', '2']
        ops = [C2x, C2y, C2z]
    if Gid == 17:
        label = ['2', '2', '2_1']
        ops = [C2x, [1, 0, 1, 0, 180, 0, 0, 0.5], [1, 0, 0, 1, 180, 0, 0, 0.5]]
    if Gid == 18:
        label = ['2_1', '2_1', '2']
        ops = [[1, 1, 0, 0, 180, 0.5, 0.5, 0], [1, 0, 1, 0, 180, 0.5, 0.5, 0], C2z]
    if Gid == 19:
        label = ['2_1', '2_1', '2_1']
        ops = [[1, 1, 0, 0, 180, 0.5, 0.5, 0], [1, 0, 1, 0, 180, 0, 0.5, 0.5], [1, 0, 0, 1, 180, 0.5, 0, 0.5]]
    if Gid == 20:
        label = ['2', '2', '2_1']
        ops = [C2x, [1, 0, 1, 0, 180, 0, 0, 0.5], [1, 0, 0, 1, 180, 0, 0, 0.5]]
    if Gid == 21 or Gid == 22 or Gid == 23:
        label = ['2', '2', '2']
        ops = [C2x, C2y, C2z]
    if Gid == 24:
        label = ['2_1', '2_1', '2_1']
        ops = [[1, 1, 0, 0, 180, 0.5, 0.5, 0], [1, 0, 1, 0, 180, 0, 0.5, 0.5], [1, 0, 0, 1, 180, 0.5, 0, 0.5]]
    if Gid == 25 or Gid == 35 or Gid == 38 or Gid == 42 or Gid == 44:
        label = ['m', 'm', '2']
        ops = [Mx, My, C2z]
    if Gid == 26 or Gid == 36:
        label = ['m', 'c', '2_1']
        ops = [Mx, comb(My, tauz), comb(C2z, tauz)]
    if Gid == 27 or Gid == 37:
        label = ['c', 'c', '2']
        ops = [comb(Mx, tauz), comb(My, tauz), C2z]
    if Gid == 28 or Gid == 40 or Gid == 46:
        label = ['m', 'a', '2']
        ops = [comb(Mx, taux), comb(My, taux), C2z]
    if Gid == 29:
        label = ['c', 'a', '2_1']
        ops = [comb(Mx, taux, tauz), comb(My, taux), comb(C2z, tauz)]
    if Gid == 30:
        label = ['n', 'c', '2']
        ops = [comb(Mx, tauy, tauz), comb(My, tauy, tauz), C2z]
    if Gid == 31:
        label = ['m', 'n', '2_1']
        ops = [Mx, comb(My, taux, tauz), comb(C2z, taux, tauz)]
    if Gid == 32 or Gid == 45:
        label = ['b', 'a', '2']
        ops = [comb(Mx, taux, tauy), comb(My, taux, tauy), C2z]
    if Gid == 33:
        label = ['n', 'a', '2_1']
        ops = [comb(Mx, taux, tauy, tauz), comb(My, taux, tauy), comb(C2z, tauz)]
    if Gid == 34:
        label = ['n', 'n', '2']
        ops = [comb(Mx, taux, tauy, tauz), comb(My, taux, tauy, tauz), C2z]
    if Gid == 30:
        label = ['n', 'c', '2']
        ops = [comb(Mx, tauy, tauz), comb(My, tauy, tauz), C2z]
    if Gid == 39:
        label = ['e', 'm', '2']
        ops = [comb(Mx, tauy), comb(My, tauy), C2z]
    if Gid == 41:
        label = ['e', 'a', '2']
        ops = [comb(Mx, taux, tauy), comb(My, taux, tauy), C2z]
    if Gid == 43:
        label = ['d', 'd', '2']
        ops = [comb(Mx, Fxyz), comb(My, Fxyz), C2z]
    if Gid == 47 or Gid == 65 or Gid == 69 or Gid == 71:
        label = ['m', 'm', 'm']
        ops = [Mx, My, Mz]
    if Gid == 48:
        label = ['n', 'n', 'n']
        ops = [comb(Mx, tauy, tauz), comb(My, taux, tauz), comb(Mz, taux, tauy)]
    if Gid == 49 or Gid == 66:
        label = ['c', 'c', 'm']
        ops = [comb(Mx, tauz), comb(My, tauz), Mz]
    if Gid == 50:
        label = ['b', 'a', 'n']
        ops = [comb(Mx, tauy), comb(My, taux), comb(Mz, taux, tauy)]
    if Gid == 51:
        label = ['m', 'm', 'a']
        ops = [comb(Mx, taux), My, comb(Mz, taux)]
    if Gid == 52:
        label = ['n', 'n', 'a']
        ops = [comb(Mx, tauy, tauz), comb(My, taux, tauy, tauz), comb(Mz, taux)]
    if Gid == 53:
        label = ['m', 'n', 'a']
        ops = [Mx, comb(My, taux, tauz), comb(Mz, taux, tauz)]
    if Gid == 54:
        label = ['c', 'c', 'a']
        ops = [comb(Mx, taux, tauz), comb(My, tauz), comb(Mz, taux)]
    if Gid == 55 or Gid == 72:
        label = ['b', 'a', 'm']
        ops = [comb(Mx, taux, tauy), comb(My, taux, tauy), Mz]
    if Gid == 56:
        label = ['c', 'c', 'n']
        ops = [comb(Mx, taux, tauz), comb(My, tauy, tauz), comb(Mz, taux, tauy)]
    if Gid == 57:
        label = ['b', 'c', 'm']
        ops = [comb(Mx, tauy), comb(My, tauy, tauz), comb(Mz, tauz)]
    if Gid == 58:
        label = ['n', 'n', 'm']
        ops = [comb(Mx, taux, tauy, tauz), comb(My, taux, tauy, tauz), Mz]
    if Gid == 59:
        label = ['m', 'm', 'n']
        ops = [comb(Mx, taux), comb(My, tauy), comb(Mz, taux, tauy)]
    if Gid == 60:
        label = ['b', 'c', 'n']
        ops = [comb(Mx, taux, tauy), comb(My, tauz), comb(Mz, taux, tauy, tauz)]
    if Gid == 61 or Gid == 73:
        label = ['b', 'c', 'a']
        ops = [comb(Mx, taux, tauy), comb(My, tauy, tauz), comb(Mz, taux, tauz)]
    if Gid == 62:
        label = ['n', 'm', 'a']
        ops = [comb(Mx, taux, tauy, tauz), comb(My, tauy), comb(Mz, taux, tauz)]
    if Gid == 63:
        label = ['m', 'c', 'm']
        ops = [Mx, comb(My, tauz), comb(Mz, tauz)]
    if Gid == 64:
        label = ['m', 'c', 'e']
        ops = [Mx, comb(My, tauy, tauz), comb(Mz, tauy, tauz)]
    if Gid == 67:
        label = ['m', 'm', 'e']
        ops = [Mx, comb(My, tauy), comb(Mz, tauy)]
    if Gid == 68:
        label = ['c', 'c', 'e']
        ops = [comb(Mx, taux, tauz), comb(My, tauz), comb(Mz, taux)]
    if Gid == 70:
        label = ['d', 'd', 'd']
        ops = [comb(Mx, Fyz), comb(My, Fxz), comb(Mz, Fxy)]
    if Gid == 74:
        label = ['m', 'm', 'a']
        ops = [comb(Mx), comb(My, tauy), comb(Mz, tauy)]
    # tetragonal part
    if Gid == 75 or Gid == 79:
        label = ['4']
        ops = [C4z]
    if Gid == 76:
        label = ['4_1']
        ops = [comb(C4z, Fz)]
    if Gid == 77:
        label = ['4_2']
        ops = [comb(C4z, tauz)]
    if Gid == 78:
        label = ['4_3']
        ops = [comb(C4z, Fz, tauz)]
    if Gid == 80:
        label = ['4_1']
        ops = [comb(C4z, Fz, tauy)]
    if Gid == 81 or Gid == 82:
        label = ['-4']
        ops = [PC4z]
    if Gid == 83 or Gid == 87:
        label = ['4', 'm']
        ops = [C4z, Mz]
    if Gid == 84:
        label = ['4_2', 'm']
        ops = [comb(C4z, tauz), Mz]
    if Gid == 85:
        label = ['4', 'n']
        ops = [comb(C4z, taux), comb(Mz, taux, tauy)]
    if Gid == 86:
        label = ['4_2', 'n']
        ops = [comb(C4z, tauy, tauz), comb(Mz, taux, tauy)]
    if Gid == 88:
        label = ['4_1', 'a']
        ops = [comb(C4z, taux, Fxyz), comb(Mz, taux, tauz)]
    if Gid == 89 or Gid == 97:
        label = ['4', '2', '2']
        ops = [C4z, C2x, C2xy]
    if Gid == 90:
        label = ['4', '2_1', '2']
        ops = [comb(C4z, taux, tauy), comb(C2x, taux, tauy), C2xy]
    if Gid == 91:
        label = ['4_1', '2', '2']
        ops = [comb(C4z, Fz), comb(C2x, tauz), comb(C2xy, Fz, tauz)]
    if Gid == 92:
        label = ['4_1', '2_1', '2']
        ops = [comb(C4z, Fz, taux, tauy), comb(C2x, Fz, taux, tauy, tauz), C2xy]
    if Gid == 93:
        label = ['4_2', '2', '2']
        ops = [comb(C4z, tauz), C2x, comb(C2xy, tauz)]
    if Gid == 94:
        label = ['4_2', '2_1', '2']
        ops = [comb(C4z, taux, tauy, tauz), comb(C2x, taux, tauy, tauz), C2xy]
    if Gid == 95:
        label = ['4_3', '2', '2']
        ops = [comb(C4z, Fz, tauz), comb(C2x, tauz), comb(C2xy, Fz)]
    if Gid == 96:
        label = ['4_3', '2_1', '2']
        ops = [comb(C4z, Fz, tauz, taux, tauy), comb(C2x, taux, tauy, Fz), C2xy]
    if Gid == 98:
        label = ['4_1', '2', '2']
        ops = [comb(C4z, Fz, tauy), comb(C2x, Fz, tauy), comb(C2xy, taux, tauy, tauz)]
    if Gid == 99 or Gid == 107:
        label = ['4', 'm', 'm']
        ops = [C4z, Mx, Mxy]
    if Gid == 100:
        label = ['4', 'b', 'm']
        ops = [C4z, comb(Mx, taux, tauy), comb(Mxy, taux, tauy)]
    if Gid == 101:
        label = ['4_2', 'c', 'm']
        ops = [comb(C4z, tauz), comb(Mx, tauz), Mxy]
    if Gid == 102:
        label = ['4_2', 'n', 'm']
        ops = [comb(C4z, taux, tauy, tauz), comb(Mx, taux, tauy, tauz), Mxy]
    if Gid == 103:
        label = ['4', 'c', 'c']
        ops = [C4z, comb(Mx, tauz), comb(Mxy, tauz)]
    if Gid == 104:
        label = ['4', 'n', 'c']
        ops = [C4z, comb(Mx, taux, tauy, tauz), comb(Mxy, taux, tauy, tauz)]
    if Gid == 105:
        label = ['4_2', 'm', 'c']
        ops = [comb(C4z, tauz), Mx, comb(Mxy, tauz)]
    if Gid == 106:
        label = ['4_2', 'b', 'c']
        ops = [comb(C4z, tauz), comb(Mx, taux, tauy), comb(Mxy, taux, tauy, tauz)]
    if Gid == 108:
        label = ['4', 'c', 'm']
        ops = [C4z, comb(Mx, tauz), comb(Mxy, tauz)]
    if Gid == 109:
        label = ['4_1', 'm', 'd']
        ops = [comb(C4z, tauy, Fz), comb(Mx, taux, tauy, tauz), comb(Mxy, tauy, Fz)]
    if Gid == 110:
        label = ['4_1', 'c', 'd']
        ops = [comb(C4z, tauy, Fz), comb(Mx, taux, tauy), comb(Mxy, tauy, Fz, tauz)]
    if Gid == 111 or Gid == 121:
        label = ['-4', '2', 'm']
        ops = [PC4z, C2x, Mxy]
    if Gid == 112:
        label = ['-4', '2', 'c']
        ops = [PC4z, comb(C2x, tauz), comb(Mxy, tauz)]
    if Gid == 113:
        label = ['-4', '2_1', 'm']
        ops = [PC4z, comb(C2x, taux, tauy), comb(Mxy, taux, tauy)]
    if Gid == 114:
        label = ['-4', '2_1', 'c']
        ops = [PC4z, comb(C2x, taux, tauy, tauz), comb(Mxy, taux, tauy, tauz)]
    if Gid == 115 or Gid == 119:
        label = ['-4', 'm', '2']
        ops = [PC4z, Mx, C2xy]
    if Gid == 116 or Gid == 120:
        label = ['-4', 'c', '2']
        ops = [PC4z, comb(Mx, tauz), comb(C2xy, tauz)]
    if Gid == 117:
        label = ['-4', 'b', '2']
        ops = [PC4z, comb(Mx, taux, tauy), comb(C2xy, taux, tauy)]
    if Gid == 118:
        label = ['-4', 'n', '2']
        ops = [PC4z, comb(Mx, taux, tauy, tauz), comb(C2xy, taux, tauy, tauz)]
    if Gid == 122:
        label = ['-4', '2', 'd']
        ops = [PC4z, comb(C2x, tauz, taux, Fz), comb(Mxy, tauz, taux, Fz)]
    if Gid == 123 or Gid == 139:
        label = ['4', 'm', 'm', 'm']
        ops = [C4z, Mz, Mx, Mxy]
    if Gid == 124:
        label = ['4', 'm', 'c', 'c']
        ops = [C4z, Mz, comb(Mx, tauz), comb(Mxy, tauz)]
    if Gid == 125:
        label = ['4', 'n', 'b', 'm']
        ops = [comb(C4z, taux), comb(Mz, taux, tauy), comb(Mx, tauy), Mxy]
    if Gid == 126:
        label = ['4', 'n', 'n', 'c']
        ops = [comb(C4z, taux), comb(Mz, taux, tauy), comb(Mx, tauy, tauz), comb(Mxy, tauz)]
    if Gid == 127:
        label = ['4', 'm', 'b', 'm']
        ops = [C4z, Mz, comb(Mx, taux, tauy), comb(Mxy, taux, tauy)]
    if Gid == 128:
        label = ['4', 'm', 'n', 'c']
        ops = [C4z, Mz, comb(Mx, taux, tauy, tauz), comb(Mxy, taux, tauy, tauz)]
    if Gid == 129:
        label = ['4', 'n', 'm', 'm']
        ops = [comb(C4z, taux), comb(Mz, taux, tauy), comb(Mx, taux), comb(Mxy, taux, tauy)]
    if Gid == 130:
        label = ['4', 'n', 'c', 'c']
        ops = [comb(C4z, taux), comb(Mz, taux, tauy), comb(Mx, taux, tauz), comb(Mxy, taux, tauy, tauz)]
    if Gid == 131:
        label = ['4_2', 'm', 'm', 'c']
        ops = [comb(C4z, tauz), Mz, Mx, comb(Mxy, tauz)]
    if Gid == 132:
        label = ['4_2', 'm', 'c', 'm']
        ops = [comb(C4z, tauz), Mz, comb(Mx, tauz), Mxy]
    if Gid == 133:
        label = ['4_2', 'n', 'b', 'c']
        ops = [comb(C4z, taux, tauz), comb(Mz, taux, tauy), comb(Mx, tauy), comb(Mxy, tauz)]
    if Gid == 134:
        label = ['4_2', 'n', 'n', 'm']
        ops = [comb(C4z, taux, tauz), comb(Mz, taux, tauy), comb(Mx, tauy, tauz), Mxy]
    if Gid == 135:
        label = ['4_2', 'm', 'b', 'c']
        ops = [comb(C4z, tauz), Mz, comb(Mx, taux, tauy), comb(Mxy, taux, tauy, tauz)]
    if Gid == 136:
        label = ['4_2', 'm', 'n', 'm']
        ops = [comb(C4z, taux, tauy, tauz), Mz, comb(Mx, taux, tauy, tauz), Mxy]
    if Gid == 137:
        label = ['4_2', 'n', 'm', 'c']
        ops = [comb(C4z, taux, tauz), comb(Mz, taux, tauy), comb(Mx, taux), comb(Mxy, taux, tauy, tauz)]
    if Gid == 138:
        label = ['4_2', 'n', 'c', 'm']
        ops = [comb(C4z, taux, tauz), comb(Mz, taux, tauy), comb(Mx, taux, tauz), comb(Mxy, taux, tauy)]
    if Gid == 140:
        label = ['4', 'm', 'c', 'm']
        ops = [C4z, Mz, comb(Mx, tauz), comb(Mxy, tauz)]
    if Gid == 141:
        label = ['4_1', 'a', 'm', 'd']
        ops = [comb(C4z, Fxyz, tauy), comb(Mz, taux, tauz), Mx, comb(Mxy, taux, tauz, Fxyz)]
    if Gid == 142:
        label = ['4_1', 'a', 'c', 'd']
        ops = [comb(C4z, Fxyz, tauy), comb(Mz, taux, tauz), comb(Mx, tauz), comb(Mxy, taux, Fxyz)]
    # trigonal part
    if Gid == 143 or Gid == 146:
        label = ['3']
        ops = [C3z]
    if Gid == 144:
        label = ['3_1']
        ops = [comb(C3z, tau3)]
    if Gid == 145:
        label = ['3_2']
        ops = [comb(C3z, tau3, tau3)]
    if Gid == 147 or Gid == 148:
        label = ['-3']
        ops = [PC3z]
    if Gid == 149:
        label = ['3', '1', '2']
        ops = [C3z, E, C2Tri]
    if Gid == 150:
        label = ['3', '2', '1']
        ops = [C3z, C2x, E]
    if Gid == 151:
        label = ['3_1', '1', '2']
        ops = [comb(C3z, tau3), E, C2Tri]
    if Gid == 152:
        label = ['3_1', '2', '1']
        ops = [comb(C3z, tau3), comb(C2x, tau3, tau3), E]
    if Gid == 153:
        label = ['3_2', '1', '2']
        ops = [comb(C3z, tau3, tau3), E, C2Tri]
    if Gid == 154:
        label = ['3_2', '2', '1']
        ops = [comb(C3z, tau3, tau3), comb(C2x, tau3), E]
    if Gid == 155:
        label = ['3', '2']
        ops = [C3z, C2x]
    if Gid == 156:
        label = ['3', 'm', '1']
        ops = [C3z, Mx, E]
    if Gid == 157:
        label = ['3', '1', 'm']
        ops = [C3z, E, MTri]
    if Gid == 158:
        label = ['3', 'c', '1']
        ops = [C3z, comb(Mx, tauz), E]
    if Gid == 159:
        label = ['3', '1', 'c']
        ops = [C3z, E, comb(MTri, tauz)]
    if Gid == 160:
        label = ['3', 'm']
        ops = [C3z, Mx]
    if Gid == 161:
        label = ['3', 'c']
        ops = [C3z, comb(Mx, tauz)]
    if Gid == 162:
        label = ['-3', '1', 'm']
        ops = [PC3z, E, MTri]
    if Gid == 163:
        label = ['-3', '1', 'c']
        ops = [PC3z, E, comb(MTri, tauz)]
    if Gid == 164:
        label = ['-3', 'm', '1']
        ops = [PC3z, Mx, E]
    if Gid == 165:
        label = ['-3', 'c', '1']
        ops = [PC3z, comb(Mx, tauz), E]
    if Gid == 166:
        label = ['-3', 'm']
        ops = [PC3z, Mx]
    if Gid == 167:
        label = ['-3', 'c']
        ops = [PC3z, comb(Mx, tauz)]
    # hexagonal part
    if Gid == 168:
        label = ['6']
        ops = [C6z]
    if Gid == 169:
        label = ['6_1']
        ops = [comb(C6z, tau6)]
    if Gid == 170:
        label = ['6_5']
        ops = [comb(C6z, tau3, tauz)]
    if Gid == 171:
        label = ['6_2']
        ops = [comb(C6z, tau3)]
    if Gid == 172:
        label = ['6_4']
        ops = [comb(C6z, tau6, tauz)]
    if Gid == 173:
        label = ['6_3']
        ops = [comb(C6z, tauz)]
    if Gid == 174:
        label = ['-6']
        ops = [PC6z]
    if Gid == 175:
        label = ['6', 'm']
        ops = [C6z, Mz]
    if Gid == 176:
        label = ['6_3', 'm']
        ops = [comb(C6z, tauz), comb(Mz, tauz)]
    if Gid == 177:
        label = ['6', '2', '2']
        ops = [C6z, C2x, C2Tri]
    if Gid == 178:
        label = ['6_1', '2', '2']
        ops = [comb(C6z, tau6), C2x, comb(C2Tri, tau6)]
    if Gid == 179:
        label = ['6_5', '2', '2']
        ops = [comb(C6z, tau3, tauz), C2x, comb(C2Tri, tau3, tauz)]
    if Gid == 180:
        label = ['6_2', '2', '2']
        ops = [comb(C6z, tau3), C2x, comb(C2Tri, tau3)]
    if Gid == 181:
        label = ['6_4', '2', '2']
        ops = [comb(C6z, tau6, tauz), C2x, comb(C2Tri, tau6, tauz)]
    if Gid == 182:
        label = ['6_3', '2', '2']
        ops = [comb(C6z, tauz), C2x, comb(C2Tri, tauz)]
    if Gid == 183:
        label = ['6', 'm', 'm']
        ops = [C6z, Mx, MTri]
    if Gid == 184:
        label = ['6', 'c', 'c']
        ops = [C6z, comb(Mx, tauz), comb(MTri, tauz)]
    if Gid == 185:
        label = ['6_3', 'c', 'm']
        ops = [comb(C6z, tauz), comb(Mx, tauz), MTri]
    if Gid == 186:
        label = ['6_3', 'm', 'c']
        ops = [comb(C6z, tauz), Mx, comb(MTri, tauz)]
    if Gid == 187:
        label = ['-6', 'm', '2']
        ops = [PC6z, Mx, C2Tri]
    if Gid == 188:
        label = ['-6', 'c', '2']
        ops = [comb(PC6z, tauz), comb(Mx, tauz), C2Tri]
    if Gid == 189:
        label = ['-6', '2', 'm']
        ops = [PC6z, C2x, MTri]
    if Gid == 190:
        label = ['-6', '2', 'c']
        ops = [comb(PC6z, tauz), C2x, comb(MTri, tauz)]
    if Gid == 191:
        label = ['6', 'm', 'm', 'm']
        ops = [C6z, Mz, Mx, MTri]
    if Gid == 192:
        label = ['6', 'm', 'c', 'c']
        ops = [C6z, Mz, comb(Mx, tauz), comb(MTri, tauz)]
    if Gid == 193:
        label = ['6_3', 'm', 'c', 'm']
        ops = [comb(PC6z, tauz), comb(Mz, tauz), comb(Mx, tauz), MTri]
    if Gid == 194:
        label = ['6_3', 'm', 'm', 'c']
        ops = [comb(PC6z, tauz), comb(Mz, tauz), Mx, comb(MTri, tauz)]
    # hexagonal part
    if Gid == 195 or Gid == 196 or Gid == 197:
        label = ['2', '3']
        ops = [C2x, C3h]
    if Gid == 198 or Gid == 199:
        label = ['2_1', '3']
        ops = [comb(C2x, taux, tauy), C3h]
    if Gid == 200 or Gid == 202 or Gid == 204:
        label = ['m', '-3']
        ops = [Mx, PC3h]
    if Gid == 201:
        label = ['n', '-3']
        ops = [comb(Mx, tauy, tauz), PC3h]
    if Gid == 203:
        label = ['d', '-3']
        ops = [comb(Mx, Fyz), PC3h]
    if Gid == 205 or Gid == 206:
        label = ['a', '-3']
        ops = [comb(Mx, taux, tauy), PC3h]
    if Gid == 207 or Gid == 209 or Gid == 211:
        label = ['4', '3', '2']
        ops = [C4x, C3h, C2xy]
    if Gid == 208:
        label = ['4_2', '3', '2']
        ops = [comb(C4x, taux, tauy, tauz), C3h, comb(C2xy, taux, tauy, tauz)]
    if Gid == 210:
        label = ['4_1', '3', '2']
        ops = [comb(C4x, Fxyz, tauy, tauz), C3h, comb(C2xy, Fxyz, taux, tauz)]
    if Gid == 212:
        label = ['4_3', '3', '2']
        ops = [comb(C4x, Fxyz, taux, tauy), C3h, comb(C2xy, Fxyz, tauy, tauz)]
    if Gid == 213 or Gid == 214:
        label = ['4_1', '3', '2']
        ops = [comb(C4x, Fxyz, tauz), C3h, comb(C2xy, Fxyz, taux)]
    if Gid == 215 or Gid == 216 or Gid == 217:
        label = ['-4', '3', 'm']
        ops = [PC4x, C3h, Mxy]
    if Gid == 218:
        label = ['-4', '3', 'n']
        ops = [comb(PC4x, taux, tauy, tauz), C3h, comb(Mxy, taux, tauy, tauz)]
    if Gid == 219:
        label = ['-4', '3', 'c']
        ops = [comb(PC4x, taux, tauy, tauz), C3h, comb(Mxy, taux, tauy, tauz)]
    if Gid == 220:
        label = ['-4', '3', 'd']
        ops = [comb(PC4x, taux, tauy, Fxyz), C3h, comb(Mxy, Fxyz, tauy, tauz)]
    if Gid == 221 or Gid == 225 or Gid == 229:
        label = ['m', '-3', 'm']
        ops = [Mx, PC3h, Mxy]
    if Gid == 222:
        label = ['n', '-3', 'n']
        ops = [comb(Mx, tauy, tauz), PC3h, comb(Mxy, tauz)]
    if Gid == 223:
        label = ['m', '-3', 'n']
        ops = [Mx, PC3h, comb(Mxy, taux, tauy, tauz)]
    if Gid == 224:
        label = ['n', '-3', 'm']
        ops = [comb(Mx, tauy, tauz), PC3h, comb(Mxy, taux, tauy)]
    if Gid == 226:
        label = ['m', '-3', 'c']
        ops = [Mx, PC3h, comb(Mxy, taux, tauy, tauz)]
    if Gid == 227:
        label = ['d', '-3', 'm']
        ops = [comb(Mx, taux, Fyz, tauz), PC3h, comb(Mxy, Fxy, tauy, tauz)]
    if Gid == 228:
        label = ['d', '-3', 'c']
        ops = [comb(Mx, taux, Fyz, tauy), PC3h, comb(Mxy, Fxy, tauy)]
    if Gid == 230:
        label = ['a', '-3', 'd']
        ops = [comb(Mx, taux, tauy), PC3h, comb(Mxy, Fxyz, tauy, tauz)]
    return label, ops



def sg_symbol_from_number(sg_number: int) -> str:
    """
    Convert space group number to international symbol.
    
    Args:
        sg_number: Space group number (1-230)
        
    Returns:
        str: International space group symbol (short form)
    """
    if not isinstance(sg_number, int) or sg_number < 1 or sg_number > 230:
        raise ValueError(f"Invalid space group number: {sg_number}. Must be integer between 1-230.")
    
    for hall in range(1, 531):
        t = spglib.get_spacegroup_type(hall)
        if t is None:
            continue
        number = getattr(t, "number", t["number"])
        if number == sg_number:
            return getattr(t, "international_short", t["international_short"])
    
    raise ValueError(f"Space group number {sg_number} not found in spglib database.")

def identity_tau(tau1, tau2, supercell, gid):
    """
    Check if two tau vectors are identical considering supercell structure.
    
    Two tau vectors are considered identical if their difference can be expressed
    as an integer linear combination of supercell basis vectors.
    
    Args:
        tau1: First tau vector
        tau2: Second tau vector  
        supercell: Supercell matrix
        gid: Space group ID
        
    Returns:
        list: [is_identical, l1, l2, l3]
            - is_identical: Boolean indicating if vectors are identical
            - l1, l2, l3: Integer coefficients if identical, 0 otherwise
    """
    # Get primitive basis vectors for the space group
    prim_vec = identify_SG_lattice(gid)[1]  # each col is a prim basis vector   
    
    # Extract supercell basis vectors
    t1 = np.array([supercell[0][0], supercell[1][0], supercell[2][0]]) 
    t2 = np.array([supercell[0][1], supercell[1][1], supercell[2][1]]) 
    t3 = np.array([supercell[0][2], supercell[1][2], supercell[2][2]])
    
    # Transform to primitive coordinates
    s1 = prim_vec @ t1
    s2 = prim_vec @ t2
    s3 = prim_vec @ t3
    
    # Set up linear system: A * [l1, l2, l3] = tau1 - tau2
    A = np.array([[s1[0], s2[0], s3[0]], 
                  [s1[1], s2[1], s3[1]], 
                  [s1[2], s2[2], s3[2]]])
    B = tau1 - tau2
    
    try:
        solve = np.linalg.solve(A, B)
        l1, l2, l3 = solve[0], solve[1], solve[2]
        
        # Check if solutions are integers
        for l in [l1, l2, l3]:
            if np.abs(l - np.round(l)) > 1e-4:
                return [False, 0, 0, 0]
        
        return [True, int(round(l1)), int(round(l2)), int(round(l3))]
    except np.linalg.LinAlgError:
        return [False, 0, 0, 0]
    
    
if __name__ == '__main__':
    gens = gen_generators(200)
    print(gens)
