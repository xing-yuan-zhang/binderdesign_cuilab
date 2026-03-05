import os
import numpy as np

from . import config
from .pdbio import parse_pdb_atoms, write_pdb_lines
from .align import get_ca_coords, kabsch, apply_rt, rmsd

def restore_one(template_atoms, template_pdb_path, design_pdb, out_pdb):
    des_atoms = parse_pdb_atoms(design_pdb)

    if not any(a["chain"] == config.binder_chain for a in des_atoms):
        print(f"[skip] {design_pdb}: no chain {config.binder_chain} found")
        return

    tpl_CA, tpl_keys = get_ca_coords(template_atoms, config.template_chain, config.contig_ranges)
    des_CA, des_keys = get_ca_coords(des_atoms, config.template_chain, config.contig_ranges)

    if len(tpl_CA) == 0 or len(tpl_CA) != len(des_CA):
        print(f"[warn] {design_pdb}: CA count mismatch or zero (tpl={len(tpl_CA)}, des={len(des_CA)}); skipping.")
        return

    if tpl_keys != des_keys:
        print(f"[warn] {design_pdb}: residue numbering mismatch between template and design; skipping.")
        return

    R, t = kabsch(des_CA, tpl_CA)

    out_lines = []
    serial = 1

    with open(template_pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                break
            if line.startswith("MODEL"):
                out_lines.append(line.rstrip("\n"))

    for a in template_atoms:
        if a["chain"] != config.template_chain:
            continue
        line = list(a["line"])
        line[6:11] = list(f"{serial:5d}")
        line[30:54] = list(f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}")
        out_lines.append("".join(line))
        serial += 1

    for a in des_atoms:
        if a["chain"] != config.binder_chain:
            continue
        v = np.array([a["x"], a["y"], a["z"]], dtype=float)
        v2 = R @ v + t
        x2, y2, z2 = v2.tolist()

        line = list(a["line"])
        line[6:11] = list(f"{serial:5d}")
        line[21] = config.binder_chain
        line[30:54] = list(f"{x2:8.3f}{y2:8.3f}{z2:8.3f}")
        out_lines.append("".join(line))
        serial += 1

    out_lines.append("END")
    write_pdb_lines(out_pdb, out_lines)

    des_CA_aligned = apply_rt(des_CA, R, t)
    r = rmsd(des_CA_aligned, tpl_CA)
    print(f"[ok] {design_pdb} -> {out_pdb} (CA RMSD={r:.3f} Å)")

def main(template_pdb, output_root, seed_dirs):
    template_atoms = parse_pdb_atoms(template_pdb)
    os.makedirs(output_root, exist_ok=True)

    for seed_dir in seed_dirs:
        seed_dir = os.path.abspath(seed_dir)
        if not os.path.isdir(seed_dir):
            print(f"[warn] {seed_dir} is not a directory, skipping.")
            continue

        seed_tag = os.path.basename(seed_dir)

        for fname in os.listdir(seed_dir):
            if not fname.lower().endswith(".pdb"):
                continue
            design_pdb = os.path.join(seed_dir, fname)

            if os.path.abspath(design_pdb) == os.path.abspath(template_pdb):
                continue

            base, ext = os.path.splitext(fname)
            out_name = f"{base}_restored_{seed_tag}{ext}"
            out_pdb = os.path.join(output_root, out_name)

            restore_one(template_atoms, template_pdb, design_pdb, out_pdb)
