import numpy as np
from . import config
from .geom import dihedral

def build_backbone_dict(chain_atoms):
    bb = {}
    for a in chain_atoms:
        name = a["atom"].strip()
        if name not in ("N", "CA", "C"):
            continue
        resi = a["resi"]
        bb.setdefault(resi, {})
        bb[resi][name] = np.array([a["x"], a["y"], a["z"]], dtype=float)
    return bb

def compute_phi_psi_for_chain(chain_atoms):
    bb = build_backbone_dict(chain_atoms)
    res_indices = sorted(bb.keys())
    phi = {}
    psi = {}

    for i, resi in enumerate(res_indices):
        if i > 0:
            prev_resi = res_indices[i - 1]
            if prev_resi == resi - 1:
                prev = bb.get(prev_resi, {})
                cur  = bb.get(resi, {})
                if "C" in prev and all(k in cur for k in ("N","CA","C")):
                    phi[resi] = dihedral(prev["C"], cur["N"], cur["CA"], cur["C"])
        if i < len(res_indices) - 1:
            next_resi = res_indices[i + 1]
            if next_resi == resi + 1:
                cur = bb.get(resi, {})
                nxt = bb.get(next_resi, {})
                if all(k in cur for k in ("N","CA","C")) and "N" in nxt:
                    psi[resi] = dihedral(cur["N"], cur["CA"], cur["C"], nxt["N"])

    return res_indices, phi, psi

def helix_segments_for_binder(binder_atoms):
    res_indices, phi, psi = compute_phi_psi_for_chain(binder_atoms)
    if not res_indices:
        return [], 0, 0

    helix_flags = {}
    for resi in res_indices:
        if resi not in phi or resi not in psi:
            helix_flags[resi] = False
            continue
        ph = phi[resi]
        ps = psi[resi]
        helix_flags[resi] = (
            config.phi_helix_min <= ph <= config.phi_helix_max and
            config.psi_helix_min <= ps <= config.psi_helix_max
        )

    segments = []
    current = []
    for resi in res_indices:
        if helix_flags.get(resi, False):
            if not current:
                current = [resi, resi]
            else:
                if resi == current[1] + 1:
                    current[1] = resi
                else:
                    segments.append(tuple(current))
                    current = [resi, resi]
        else:
            if current:
                segments.append(tuple(current))
                current = []
    if current:
        segments.append(tuple(current))

    filtered = []
    total_helix_res = 0
    for start, end in segments:
        length = end - start + 1
        if length >= config.min_helix_len:
            filtered.append((start, end, length))
            total_helix_res += length

    return filtered, len(filtered), total_helix_res
