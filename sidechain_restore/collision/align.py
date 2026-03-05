import numpy as np

def get_ca_coords(atoms, chain, ranges):
    coords = []
    reskeys = []
    for a in atoms:
        if a["chain"] != chain:
            continue
        if a["atom"].strip() != "CA":
            continue
        r = a["resi"]
        if any(lo <= r <= hi for lo, hi in ranges):
            coords.append([a["x"], a["y"], a["z"]])
            reskeys.append(r)
    return np.array(coords, dtype=float), reskeys

def kabsch(P, Q):
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    C = P0.T @ Q0
    V, S, Wt = np.linalg.svd(C)
    d = np.linalg.det(V @ Wt)
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    t = Qc - R @ Pc
    return R, t

def apply_rt(coords, R, t):
    return (R @ coords.T).T + t

def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))
