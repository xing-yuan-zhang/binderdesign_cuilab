import json
from collections import OrderedDict
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
}


def _np():
    import numpy as np

    return np


def _parse_pdb_atoms(path):
    atoms = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            altloc = line[16].strip()
            if altloc not in {"", "A"}:
                continue
            try:
                atoms.append(
                    {
                        "line": line.rstrip("\n"),
                        "record": line[0:6].strip(),
                        "atom": line[12:16].strip(),
                        "resn": line[17:20].strip(),
                        "chain": line[21].strip() or "_",
                        "resi": int(line[22:26]),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "bfactor": float(line[60:66] or 0.0),
                    }
                )
            except ValueError:
                continue
    return atoms


def _parse_cif_atoms(path):
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("gemmi is required to parse mmCIF files") from exc

    structure = gemmi.read_structure(str(path))
    atoms = []
    if len(structure) == 0:
        return atoms
    model = structure[0]
    for chain in model:
        chain_name = (chain.name or "_").strip() or "_"
        for residue in chain:
            try:
                seqid = residue.seqid.num
            except Exception:
                continue
            resn = residue.name.strip()
            for atom in residue:
                altloc = getattr(atom, "altloc", "") or ""
                if str(altloc).strip() not in {"", "A"}:
                    continue
                pos = atom.pos
                atoms.append(
                    {
                        "line": "",
                        "record": "ATOM",
                        "atom": atom.name.strip(),
                        "resn": resn,
                        "chain": chain_name,
                        "resi": int(seqid),
                        "x": float(pos.x),
                        "y": float(pos.y),
                        "z": float(pos.z),
                        "bfactor": float(getattr(atom, "b_iso", 0.0) or 0.0),
                    }
                )
    return atoms


def parse_structure_atoms(path):
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return _parse_cif_atoms(path)
    return _parse_pdb_atoms(path)


def parse_pdb_atoms(path):
    return parse_structure_atoms(path)


def atoms_by_chain(atoms, chain):
    return [atom for atom in atoms if atom["chain"] == chain]


def atoms_to_array(atoms):
    np = _np()
    if not atoms:
        return np.zeros((0, 3), dtype=float)
    return np.array([[atom["x"], atom["y"], atom["z"]] for atom in atoms], dtype=float)


def residues_for_chain(atoms, chain):
    residues = OrderedDict()
    for atom in atoms:
        if atom["chain"] != chain:
            continue
        key = atom["resi"]
        residues.setdefault(
            key,
            {
                "resi": atom["resi"],
                "resn": atom["resn"],
                "atoms": [],
            },
        )
        residues[key]["atoms"].append(atom)
    return list(residues.values())


def chain_sequence(atoms, chain):
    residues = residues_for_chain(atoms, chain)
    sequence = "".join(AA3_TO_1.get(residue["resn"], "X") for residue in residues)
    return sequence, residues


def chain_ca_coords_by_order(atoms, chain, positions=None):
    np = _np()
    residues = residues_for_chain(atoms, chain)
    coords = []
    residue_positions = []
    wanted = set(positions) if positions else None
    for index, residue in enumerate(residues, start=1):
        if wanted is not None and index not in wanted:
            continue
        ca = None
        for atom in residue["atoms"]:
            if atom["atom"] == "CA":
                ca = atom
                break
        if ca is None:
            continue
        coords.append([ca["x"], ca["y"], ca["z"]])
        residue_positions.append(index)
    return np.array(coords, dtype=float), residue_positions


def chain_break_metrics(atoms, chain, max_ca_gap=4.5):
    np = _np()
    coords, _ = chain_ca_coords_by_order(atoms, chain)
    if len(coords) < 2:
        return 0, 0.0
    diffs = coords[1:] - coords[:-1]
    gaps = np.sqrt((diffs * diffs).sum(axis=1))
    chain_breaks = int((gaps > max_ca_gap).sum())
    max_gap = float(gaps.max()) if len(gaps) else 0.0
    return chain_breaks, max_gap


def min_interchain_distance(atoms, chain_a, chain_b):
    np = _np()
    coords_a = atoms_to_array(atoms_by_chain(atoms, chain_a))
    coords_b = atoms_to_array(atoms_by_chain(atoms, chain_b))
    if coords_a.shape[0] == 0 or coords_b.shape[0] == 0:
        return None
    delta = coords_a[:, None, :] - coords_b[None, :, :]
    dist2 = np.sum(delta * delta, axis=2)
    return float(np.sqrt(dist2.min()))


def atom_contact_count(atoms, chain_a, chain_b, cutoff=5.0):
    np = _np()
    coords_a = atoms_to_array(atoms_by_chain(atoms, chain_a))
    coords_b = atoms_to_array(atoms_by_chain(atoms, chain_b))
    if coords_a.shape[0] == 0 or coords_b.shape[0] == 0:
        return 0
    delta = coords_a[:, None, :] - coords_b[None, :, :]
    dist2 = np.sum(delta * delta, axis=2)
    return int((dist2 < cutoff * cutoff).sum())


def kabsch(P, Q):
    np = _np()
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    covariance = P0.T @ Q0
    V, _, Wt = np.linalg.svd(covariance)
    d = np.linalg.det(V @ Wt)
    D = np.diag([1.0, 1.0, d])
    rotation = V @ D @ Wt
    translation = Qc - rotation @ Pc
    return rotation, translation


def apply_rt(coords, rotation, translation):
    return (rotation @ coords.T).T + translation


def rmsd(P, Q):
    np = _np()
    if len(P) == 0 or len(P) != len(Q):
        return None
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))


def tm_score(P, Q):
    np = _np()
    if len(P) == 0 or len(P) != len(Q):
        return None
    length = len(P)
    if length < 3:
        return 1.0 if length else None
    d0 = max(0.5, 1.24 * np.cbrt(length - 15.0) - 1.8)
    distances = np.sqrt(((P - Q) ** 2).sum(axis=1))
    return float(np.mean(1.0 / (1.0 + (distances / d0) ** 2)))


def _trim_pair(left, right):
    np = _np()
    size = min(len(left), len(right))
    if size == 0:
        return np.zeros((0, 3), dtype=float), np.zeros((0, 3), dtype=float)
    return left[:size], right[:size]


def parse_motif_positions(spec):
    if not spec:
        return []
    positions = []
    for chunk in str(spec).split(","):
        token = chunk.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            positions.extend(range(int(start), int(end) + 1))
        else:
            positions.append(int(token))
    return sorted(set(positions))


def binder_metrics_after_receptor_alignment(
    reference_pdb,
    model_pdb,
    receptor_chain="A",
    binder_chain="B",
    motif_positions=None,
):
    ref_atoms = parse_structure_atoms(reference_pdb)
    model_atoms = parse_structure_atoms(model_pdb)

    ref_rec, _ = chain_ca_coords_by_order(ref_atoms, receptor_chain)
    mod_rec, _ = chain_ca_coords_by_order(model_atoms, receptor_chain)
    ref_bind, _ = chain_ca_coords_by_order(ref_atoms, binder_chain)
    mod_bind, _ = chain_ca_coords_by_order(model_atoms, binder_chain)

    ref_rec, mod_rec = _trim_pair(ref_rec, mod_rec)
    ref_bind, mod_bind = _trim_pair(ref_bind, mod_bind)

    if len(ref_rec) >= 3:
        rotation, translation = kabsch(mod_rec, ref_rec)
        mod_bind_aligned = apply_rt(mod_bind, rotation, translation)
    elif len(ref_bind) >= 3:
        rotation, translation = kabsch(mod_bind, ref_bind)
        mod_bind_aligned = apply_rt(mod_bind, rotation, translation)
    else:
        return {
            "af_rmsd_binder": None,
            "af_tm_binder": None,
            "af_motif_rmsd": None,
        }

    metrics = {
        "af_rmsd_binder": rmsd(mod_bind_aligned, ref_bind),
        "af_tm_binder": tm_score(mod_bind_aligned, ref_bind),
        "af_motif_rmsd": None,
    }

    motif_positions = motif_positions or []
    if motif_positions:
        ref_motif, _ = chain_ca_coords_by_order(ref_atoms, binder_chain, positions=motif_positions)
        mod_motif, _ = chain_ca_coords_by_order(model_atoms, binder_chain, positions=motif_positions)
        ref_motif, mod_motif = _trim_pair(ref_motif, mod_motif)
        if len(ref_motif) and len(mod_motif):
            mod_motif_aligned = apply_rt(mod_motif, rotation, translation)
            metrics["af_motif_rmsd"] = rmsd(mod_motif_aligned, ref_motif)

    return metrics


def mean_plddt_from_pdb(path):
    np = _np()
    atoms = parse_structure_atoms(path)
    values = [atom["bfactor"] for atom in atoms if atom["atom"] == "CA"]
    if not values:
        values = [atom["bfactor"] for atom in atoms]
    if not values:
        return None
    return float(np.mean(values))


def load_json_candidates(folder):
    payloads = []
    for path in sorted(Path(folder).glob("*.json")):
        try:
            with open(path, encoding="utf-8") as handle:
                payloads.append((path, json.load(handle)))
        except (json.JSONDecodeError, OSError):
            continue
    return payloads


def _walk_json_values(payload):
    if isinstance(payload, dict):
        for key, value in payload.items():
            yield key, value
            yield from _walk_json_values(value)
    elif isinstance(payload, list):
        for value in payload:
            yield from _walk_json_values(value)


def _first_numeric(payload, keys):
    wanted = {item.lower() for item in keys}
    for key, value in _walk_json_values(payload):
        if str(key).lower() not in wanted:
            continue
        if isinstance(value, (float, int)):
            return float(value)
    return None


def _first_array(payload, keys):
    wanted = {item.lower() for item in keys}
    for key, value in _walk_json_values(payload):
        if str(key).lower() not in wanted:
            continue
        if isinstance(value, list) and value:
            return value
    return None


def extract_af_json_metrics(folder, receptor_len=None, binder_len=None):
    np = _np()
    metrics = {
        "af_plddt": None,
        "af_iptm": None,
        "af_ptm": None,
        "af_pae_interface": None,
    }
    for _, payload in load_json_candidates(folder):
        if not isinstance(payload, dict):
            continue

        if metrics["af_plddt"] is None:
            plddt = _first_array(payload, ["plddt"])
            if isinstance(plddt, list) and plddt:
                try:
                    metrics["af_plddt"] = float(np.mean(np.array(plddt, dtype=float)))
                except Exception:
                    pass

        if metrics["af_iptm"] is None:
            metrics["af_iptm"] = _first_numeric(payload, ["iptm", "ipTM"])

        if metrics["af_ptm"] is None:
            metrics["af_ptm"] = _first_numeric(payload, ["ptm", "pTM"])

        if receptor_len and binder_len and metrics["af_pae_interface"] is None:
            pae = _first_array(payload, ["pae", "predicted_aligned_error", "pae_matrix"])
            if isinstance(pae, list) and pae:
                try:
                    pae_array = np.array(pae, dtype=float)
                    left = pae_array[:receptor_len, receptor_len : receptor_len + binder_len]
                    right = pae_array[receptor_len : receptor_len + binder_len, :receptor_len]
                    values = []
                    if left.size:
                        values.append(left.mean())
                    if right.size:
                        values.append(right.mean())
                    if values:
                        metrics["af_pae_interface"] = float(np.mean(values))
                except Exception:
                    pass
    return metrics


def global_sequence_identity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    rows = len(seq_a) + 1
    cols = len(seq_b) + 1
    score = [[0] * cols for _ in range(rows)]
    trace = [[None] * cols for _ in range(rows)]

    for i in range(1, rows):
        score[i][0] = -i
        trace[i][0] = "U"
    for j in range(1, cols):
        score[0][j] = -j
        trace[0][j] = "L"

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = 2 if seq_a[i - 1] == seq_b[j - 1] else -1
            diag = score[i - 1][j - 1] + match_score
            up = score[i - 1][j] - 1
            left = score[i][j - 1] - 1
            best = max(diag, up, left)
            score[i][j] = best
            trace[i][j] = "D" if best == diag else ("U" if best == up else "L")

    i = len(seq_a)
    j = len(seq_b)
    matches = 0
    aligned = 0
    while i > 0 or j > 0:
        step = trace[i][j]
        if step == "D":
            aligned += 1
            if seq_a[i - 1] == seq_b[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif step == "U":
            aligned += 1
            i -= 1
        else:
            aligned += 1
            j -= 1

    return matches / aligned if aligned else 0.0
