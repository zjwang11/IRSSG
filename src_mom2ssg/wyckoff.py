import numpy as np


def wrap01(f):
    f = np.asarray(f, float)
    return f - np.floor(f)

def get_swyckoff(cell, operations, tol=1e-4):

    lattice, position, numbers, elements, mag = cell
    N = len(position)

    def wrap_m05(x): return (np.asarray(x, float) + 0.5) % 1.0 - 0.5
    def same_mod1(a,b,atol=tol):
        d = np.asarray(a) - np.asarray(b); d -= np.rint(d)
        return np.max(np.abs(d)) < atol
    def l2_mod1(a,b):
        d = np.asarray(a) - np.asarray(b); d -= np.rint(d)
        return np.linalg.norm(d)

    def build_orbits(pos_i, ops):
        m = len(pos_i); used = np.zeros(m, bool); orbits = []
        for s in range(m):
            if used[s]: continue
            o = {s}; stack = [s]
            while stack:
                p = stack.pop()
                x = pos_i[p]
                for iop in range(len(ops['RotC'])):
                    R = np.asarray(ops['RotC'][iop]); t = np.asarray(ops['TauC'][iop])
                    y = wrap_m05(x @ R.T + t)
                    diffs = pos_i - y; diffs -= np.rint(diffs)
                    hit = np.where(np.max(np.abs(diffs),axis=1) < tol)[0]
                    for q in hit:
                        if q not in o:
                            o.add(q); stack.append(q)
            idx = np.array(sorted(o), int); used[idx]=True; orbits.append(idx)
        return orbits

    def stabilizer(x, ops):
        idx = []
        for iop in range(len(ops['RotC'])):
            R = np.asarray(ops['RotC'][iop]); t = np.asarray(ops['TauC'][iop])
            if same_mod1(wrap_m05(x @ R.T + t), x):
                idx.append(iop)
        return idx

    def symmetrize_rep(xrep, ops):
        idx = stabilizer(xrep, ops)
        if len(idx) <= 1:
            return wrap01(xrep)
        base = np.asarray(xrep, float)
        acc = np.zeros(3, float); cnt = 0
        for iop in idx:
            R = np.asarray(ops['RotC'][iop]); t = np.asarray(ops['TauC'][iop])
            y = wrap_m05(base @ R.T + t)
            d = y - base; d -= np.rint(d)
            acc += d; cnt += 1
        return wrap01(base + acc / max(cnt,1))

    def generate_member_from_rep(xrep_sym, x_target, ops):
        best = None; best_d = 1e9; best_op = None
        for iop in range(len(ops['RotC'])):
            R = np.asarray(ops['RotC'][iop]); t = np.asarray(ops['TauC'][iop])
            cand = wrap01(wrap_m05(xrep_sym) @ R.T + t)
            d = l2_mod1(cand, x_target)
            if d < best_d - 1e-12:
                best_d = d; best = cand; best_op = iop
        return best, best_op

    numbers = np.asarray(numbers)  
    pos_out = wrap01(position.copy())
    mag_out = np.array(mag, float, copy=True)

    wyckoff_out = {}

    idx_by_num = {}
    for i, z in enumerate(numbers):
        idx_by_num.setdefault(int(z), []).append(i)

    for z, idxs in idx_by_num.items():
        pos_i = wrap01(position[idxs])
        mag_i = np.asarray(mag[idxs], float)

        orbits = build_orbits(pos_i, operations)
        wy_dict = {}

        for orb in orbits:
            rep_loc = int(orb[0])
            xrep0   = pos_i[rep_loc]
            xrep_sym = symmetrize_rep(xrep0, operations)

            block_pos = []
            block_mag = []
            for k in orb:
                x_new, op_k = generate_member_from_rep(xrep_sym, pos_i[k], operations)
                block_pos.append(x_new)

                block_mag.append(mag_i[k])

                pos_out[idxs[k]] = x_new
                mag_out[idxs[k]] = block_mag[-1]

            wy_dict[idxs[rep_loc]] = {
                "position": np.array(block_pos, float),
                "magmom":   np.array(block_mag, float),
                "symmetry": [  
                    [np.asarray(operations['RotC'][iop]),
                     np.asarray(operations['TauC'][iop]),
                     np.asarray(operations['spin'][iop])]
                    for iop in stabilizer(xrep_sym, operations)
                ],
            }
        wyckoff_out[z] = wy_dict

    elements_out = elements

    cell_new = [np.array(lattice, float),
                pos_out,                 
                numbers.copy(),
                elements_out,
                mag_out]
    return cell_new, wyckoff_out


def output_wyckoff(cell, wyckoff, filename='swyckoff.out'):
    elements = cell[3]
    with open(filename, 'w') as f:
        for element_id, wyckoff_dict in wyckoff.items():
            element = elements[element_id-1]
            
            
            for data in wyckoff_dict.values():
                for i in range(len(data['position'])):
                    if i == 0:
                        f.write(f"{len(data['position']):5d} {element:>5}   ")
                    else:
                        f.write("              ")
                    f.write(f" {data['position'][i][0]:>12.6f}   {data['position'][i][1]:>12.6f}   {data['position'][i][2]:>12.6f}   {data['magmom'][i][0]:>12.6f}   {data['magmom'][i][1]:>12.6f}   {data['magmom'][i][2]:>12.6f}\n")

                # f.write("\n")
