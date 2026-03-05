import os, glob, shutil
import numpy as np

from . import config
from .pdbio import (
    parse_pdb_atoms, atoms_to_array, get_chain_atoms,
    get_chain_CA_atoms, get_hotspot_atoms
)
from .geom import compute_Rg
from .helix import helix_segments_for_binder

def compute_score(metrics):
    return 100.0 * float(metrics["hotspots_contacted"]) + float(metrics["iface_contacts"])

def triage_one(path):
    atoms = parse_pdb_atoms(path)

    rec_atoms  = get_chain_atoms(atoms, config.receptor_chain)
    bind_atoms = get_chain_atoms(atoms, config.binder_chain)

    if not bind_atoms:
        return (False, f"no binder chain {config.binder_chain}", {})

    rec_coords  = atoms_to_array(rec_atoms)
    bind_coords = atoms_to_array(bind_atoms)

    if rec_coords.shape[0] == 0:
        return (False, "no receptor atoms", {})

    diff = rec_coords[:, None, :] - bind_coords[None, :, :]
    dist2 = np.sum(diff * diff, axis=2)
    dists = np.sqrt(dist2)

    min_AB_dist     = float(dists.min())
    iface_contacts  = int((dists < config.iface_cutoff).sum())
    clash           = min_AB_dist < config.min_allowed_dist

    hotspots_contacted = 0
    for h in config.hotspots:
        h_atoms = get_hotspot_atoms(rec_atoms, h)
        if not h_atoms:
            continue
        h_coords = atoms_to_array(h_atoms)
        diff_hb  = h_coords[:, None, :] - bind_coords[None, :, :]
        d2_hb    = np.sum(diff_hb * diff_hb, axis=2)
        min_hb   = float(np.sqrt(d2_hb.min()))
        if min_hb < config.hotspot_cutoff:
            hotspots_contacted += 1

    bind_CA_atoms   = get_chain_CA_atoms(atoms, config.binder_chain)
    bind_CA_coords  = atoms_to_array(bind_CA_atoms)
    binder_len      = bind_CA_coords.shape[0]

    if binder_len == 0:
        return (False, "no binder CA atoms", {})

    Rg = compute_Rg(bind_CA_coords)

    diff_bb = bind_CA_coords[:, None, :] - bind_CA_coords[None, :, :]
    dist2_bb = np.sum(diff_bb * diff_bb, axis=2)
    np.fill_diagonal(dist2_bb, np.inf)
    d_bb = np.sqrt(dist2_bb)

    internal_contacts = int((d_bb < config.bb_contact_cutoff).sum() / 2)
    neighbors_per_res = (d_bb < config.bb_contact_cutoff).sum(axis=1)
    min_neighbors     = int(neighbors_per_res.min())

    helix_segments, helix_count, total_helix_res = helix_segments_for_binder(bind_atoms)

    reasons = []

    if not (config.binder_len_min <= binder_len <= config.binder_len_max):
        reasons.append(f"binder_len={binder_len} outside [{config.binder_len_min},{config.binder_len_max}]")

    if not (config.Rg_min <= Rg <= config.Rg_max):
        reasons.append(f"Rg={Rg:.2f} outside [{config.Rg_min},{config.Rg_max}]")

    if internal_contacts < config.min_internal_CA_contacts:
        reasons.append(f"internal_CA_contacts={internal_contacts} < {config.min_internal_CA_contacts}")

    if min_neighbors < config.min_neighbors_per_residue:
        reasons.append(f"min_neighbors_per_residue={min_neighbors} < {config.min_neighbors_per_residue}")

    if iface_contacts < config.min_iface_contacts:
        reasons.append(f"iface_contacts={iface_contacts} < {config.min_iface_contacts}")

    if hotspots_contacted < config.min_hotspots_contacted:
        reasons.append(f"hotspots_contacted={hotspots_contacted} < {config.min_hotspots_contacted}")

    if clash:
        reasons.append(f"clash: min_AB_dist={min_AB_dist:.2f} < {config.min_allowed_dist}")

    if config.use_helix_filter and helix_count < config.min_helices:
        reasons.append(f"helix_count={helix_count} < {config.min_helices}")

    ok = len(reasons) == 0

    metrics = {
        "binder_len": binder_len,
        "Rg": Rg,
        "internal_CA_contacts": internal_contacts,
        "min_neighbors_per_residue": min_neighbors,
        "iface_contacts": iface_contacts,
        "hotspots_contacted": hotspots_contacted,
        "min_AB_dist": min_AB_dist,
        "helix_count": helix_count,
        "total_helix_res": total_helix_res,
    }

    reason_str = "PASS" if ok else "; ".join(reasons)
    return ok, reason_str, metrics

def main(folder):
    folder = os.path.abspath(folder)
    keepers_dir = os.path.join(folder, "keepers")
    os.makedirs(keepers_dir, exist_ok=True)

    pdbs = sorted(glob.glob(os.path.join(folder, "*.pdb")))
    if not pdbs:
        print(f"No PDBs found in {folder}")
        return

    header = (
        "design\tstatus\tscore\tbinder_len\tRg\tiface_contacts\t"
        "hotspots_contacted\tmin_AB_dist\tinternal_CA_contacts\t"
        "min_neighbors\thelix_count\ttotal_helix_res\treason"
    )
    print(header)
    rows = [header]

    for pdb_path in pdbs:
        base = os.path.basename(pdb_path)
        ok, reason, m = triage_one(pdb_path)

        if not m:
            row = f"{base}\tFAIL\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t{reason}"
            print(row)
            rows.append(row)
            continue

        score = compute_score(m)

        row = (
            f"{base}\t"
            f"{'PASS' if ok else 'FAIL'}\t"
            f"{score:.3f}\t"
            f"{m['binder_len']}\t"
            f"{m['Rg']:.2f}\t"
            f"{m['iface_contacts']}\t"
            f"{m['hotspots_contacted']}\t"
            f"{m['min_AB_dist']:.2f}\t"
            f"{m['internal_CA_contacts']}\t"
            f"{m['min_neighbors_per_residue']}\t"
            f"{m['helix_count']}\t"
            f"{m['total_helix_res']}\t"
            f"{reason}"
        )

        print(row)
        rows.append(row)

        if ok:
            out_path = os.path.join(keepers_dir, base)
            shutil.copy2(pdb_path, out_path)

    tsv_path = os.path.join(folder, "triage_scores.tsv")
    with open(tsv_path, "w") as f:
        for r in rows:
            f.write(r + "\n")

    print(f"\n[info] Wrote triage scores to: {tsv_path}")
    print(f"[info] PASS designs copied to: {keepers_dir}")
