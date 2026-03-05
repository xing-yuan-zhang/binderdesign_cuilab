import numpy as np

def compute_Rg(coords):
    if coords.shape[0] == 0:
        return 0.0
    center = coords.mean(axis=0)
    diffs = coords - center
    return float(np.sqrt((diffs * diffs).sum(axis=1).mean()))

def dihedral(p0, p1, p2, p3):
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b1 /= np.linalg.norm(b1)

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))
