from dataclasses import dataclass

import numpy as np

from .structure import atoms_by_chain, atoms_to_array, parse_structure_atoms


def compute_rg(coords):
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


@dataclass
class GeometryConfig:
    receptor_chain: str = "A"
    binder_chain: str = "B"
    iface_cutoff: float = 5.0
    min_iface_contacts: int = 40
    hotspots: tuple[str, ...] = ()
    hotspot_cutoff: float = 8.0
    min_hotspots_contacted: int = 3
    min_allowed_dist: float = 2.0
    binder_len_min: int = 60
    binder_len_max: int = 100
    rg_min: float = 7.0
    rg_max: float = 20.0
    bb_contact_cutoff: float = 8.0
    min_internal_ca_contacts: int = 80
    min_neighbors_per_residue: int = 2
    use_helix_filter: bool = True
    min_helices: int = 3
    min_helix_len: int = 5
    phi_helix_min: float = -120.0
    phi_helix_max: float = -30.0
    psi_helix_min: float = -80.0
    psi_helix_max: float = 45.0


def parse_hotspots(spec):
    if isinstance(spec, (list, tuple)):
        return tuple(str(item).strip() for item in spec if str(item).strip())
    return tuple(token.strip() for token in str(spec).split(",") if token.strip())


def compute_score(metrics):
    return 100.0 * float(metrics["hotspots_contacted"]) + float(metrics["iface_contacts"])


def _chain_ca_atoms(atoms, chain_id):
    return [atom for atom in atoms if atom["chain"] == chain_id and atom["atom"].strip() == "CA"]


def _get_hotspot_atoms(atoms, hotspot_res):
    chain = hotspot_res[0]
    resi = int(hotspot_res[1:])
    return [atom for atom in atoms if atom["chain"] == chain and atom["resi"] == resi]


def _build_backbone_dict(chain_atoms):
    bb = {}
    for atom in chain_atoms:
        name = atom["atom"].strip()
        if name not in {"N", "CA", "C"}:
            continue
        bb.setdefault(atom["resi"], {})
        bb[atom["resi"]][name] = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
    return bb


def _compute_phi_psi_for_chain(chain_atoms):
    bb = _build_backbone_dict(chain_atoms)
    residues = sorted(bb.keys())
    phi = {}
    psi = {}

    for index, resi in enumerate(residues):
        if index > 0:
            prev_resi = residues[index - 1]
            if prev_resi == resi - 1:
                prev = bb.get(prev_resi, {})
                cur = bb.get(resi, {})
                if "C" in prev and all(key in cur for key in ("N", "CA", "C")):
                    phi[resi] = dihedral(prev["C"], cur["N"], cur["CA"], cur["C"])
        if index < len(residues) - 1:
            next_resi = residues[index + 1]
            if next_resi == resi + 1:
                cur = bb.get(resi, {})
                nxt = bb.get(next_resi, {})
                if all(key in cur for key in ("N", "CA", "C")) and "N" in nxt:
                    psi[resi] = dihedral(cur["N"], cur["CA"], cur["C"], nxt["N"])

    return residues, phi, psi


def helix_segments_for_binder(binder_atoms, config):
    residues, phi, psi = _compute_phi_psi_for_chain(binder_atoms)
    if not residues:
        return [], 0, 0

    helix_flags = {}
    for resi in residues:
        if resi not in phi or resi not in psi:
            helix_flags[resi] = False
            continue
        helix_flags[resi] = (
            config.phi_helix_min <= phi[resi] <= config.phi_helix_max
            and config.psi_helix_min <= psi[resi] <= config.psi_helix_max
        )

    segments = []
    current = []
    for resi in residues:
        if helix_flags.get(resi, False):
            if not current:
                current = [resi, resi]
            elif resi == current[1] + 1:
                current[1] = resi
            else:
                segments.append(tuple(current))
                current = [resi, resi]
        elif current:
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


def triage_design(path, config):
    atoms = parse_structure_atoms(path)

    rec_atoms = atoms_by_chain(atoms, config.receptor_chain)
    bind_atoms = atoms_by_chain(atoms, config.binder_chain)

    if not bind_atoms:
        return False, f"no binder chain {config.binder_chain}", {}
    if not rec_atoms:
        return False, f"no receptor chain {config.receptor_chain}", {}

    rec_coords = atoms_to_array(rec_atoms)
    bind_coords = atoms_to_array(bind_atoms)

    diff = rec_coords[:, None, :] - bind_coords[None, :, :]
    dist2 = np.sum(diff * diff, axis=2)
    dists = np.sqrt(dist2)

    min_ab_dist = float(dists.min())
    iface_contacts = int((dists < config.iface_cutoff).sum())
    clash = min_ab_dist < config.min_allowed_dist

    hotspots_contacted = 0
    for hotspot in config.hotspots:
        hotspot_atoms = _get_hotspot_atoms(rec_atoms, hotspot)
        if not hotspot_atoms:
            continue
        hotspot_coords = atoms_to_array(hotspot_atoms)
        diff_hb = hotspot_coords[:, None, :] - bind_coords[None, :, :]
        d2_hb = np.sum(diff_hb * diff_hb, axis=2)
        min_hb = float(np.sqrt(d2_hb.min()))
        if min_hb < config.hotspot_cutoff:
            hotspots_contacted += 1

    bind_ca_atoms = _chain_ca_atoms(atoms, config.binder_chain)
    bind_ca_coords = atoms_to_array(bind_ca_atoms)
    binder_len = bind_ca_coords.shape[0]
    if binder_len == 0:
        return False, f"no binder CA atoms for chain {config.binder_chain}", {}

    rg = compute_rg(bind_ca_coords)

    diff_bb = bind_ca_coords[:, None, :] - bind_ca_coords[None, :, :]
    dist2_bb = np.sum(diff_bb * diff_bb, axis=2)
    np.fill_diagonal(dist2_bb, np.inf)
    d_bb = np.sqrt(dist2_bb)

    internal_contacts = int((d_bb < config.bb_contact_cutoff).sum() / 2)
    neighbors_per_res = (d_bb < config.bb_contact_cutoff).sum(axis=1)
    min_neighbors = int(neighbors_per_res.min())

    _, helix_count, total_helix_res = helix_segments_for_binder(bind_atoms, config)

    reasons = []
    if not (config.binder_len_min <= binder_len <= config.binder_len_max):
        reasons.append(f"binder_len={binder_len} outside [{config.binder_len_min},{config.binder_len_max}]")
    if not (config.rg_min <= rg <= config.rg_max):
        reasons.append(f"Rg={rg:.2f} outside [{config.rg_min},{config.rg_max}]")
    if internal_contacts < config.min_internal_ca_contacts:
        reasons.append(f"internal_CA_contacts={internal_contacts} < {config.min_internal_ca_contacts}")
    if min_neighbors < config.min_neighbors_per_residue:
        reasons.append(f"min_neighbors_per_residue={min_neighbors} < {config.min_neighbors_per_residue}")
    if iface_contacts < config.min_iface_contacts:
        reasons.append(f"iface_contacts={iface_contacts} < {config.min_iface_contacts}")
    if hotspots_contacted < config.min_hotspots_contacted:
        reasons.append(f"hotspots_contacted={hotspots_contacted} < {config.min_hotspots_contacted}")
    if clash:
        reasons.append(f"clash: min_AB_dist={min_ab_dist:.2f} < {config.min_allowed_dist}")
    if config.use_helix_filter and helix_count < config.min_helices:
        reasons.append(f"helix_count={helix_count} < {config.min_helices}")

    metrics = {
        "binder_len": binder_len,
        "Rg": rg,
        "internal_CA_contacts": internal_contacts,
        "min_neighbors_per_residue": min_neighbors,
        "iface_contacts": iface_contacts,
        "hotspots_contacted": hotspots_contacted,
        "min_AB_dist": min_ab_dist,
        "helix_count": helix_count,
        "total_helix_res": total_helix_res,
    }
    return len(reasons) == 0, ("PASS" if not reasons else "; ".join(reasons)), metrics
