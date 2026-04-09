import argparse
from pathlib import Path

from .common import (
    build_rank_key,
    coerce_float,
    infer_config_from_parts,
    infer_seed_from_parts,
    list_json_files,
    merge_rows,
    read_json,
    read_tsv,
    row_matches,
    sanitize_token,
    truthy,
    write_tsv,
)


def _structure():
    from . import structure

    return structure


def _iter_design_files(root, patterns):
    seen = set()
    for pattern in patterns:
        for path in sorted(root.rglob(pattern)):
            token = str(path.resolve())
            if token in seen or not path.is_file():
                continue
            seen.add(token)
            yield path


def cmd_scan_designs(args):
    structure = _structure()
    root = Path(args.design_root).resolve()
    exclude_dirs = set(args.exclude_dir or [])
    rows = []
    seen_ids = set()

    for path in _iter_design_files(root, args.pattern):
        if any(part in exclude_dirs for part in path.parts):
            continue
        atoms = structure.parse_structure_atoms(path)
        binder_seq, _ = structure.chain_sequence(atoms, args.binder_chain)
        seed = infer_seed_from_parts(path.parts)
        config_name = infer_config_from_parts(path.parts)
        design_name = path.stem
        design_id = sanitize_token(
            args.id_format.format(
                config_name=config_name,
                seed=seed or "noseed",
                stem=design_name,
            )
        )
        if design_id in seen_ids:
            suffix = 2
            while f"{design_id}_{suffix}" in seen_ids:
                suffix += 1
            design_id = f"{design_id}_{suffix}"
        seen_ids.add(design_id)
        rows.append(
            {
                "design_id": design_id,
                "design_name": design_name,
                "config_name": config_name,
                "seed": seed,
                "source_pdb": str(path),
                "binder_chain": args.binder_chain,
                "receptor_chain": args.receptor_chain,
                "binder_sequence": binder_seq,
                "binder_len_scan": str(len(binder_seq)),
            }
        )

    write_tsv(args.output, rows)


def _parse_columns_arg(columns):
    if not columns:
        return []
    return [token.strip() for token in columns.split(",") if token.strip()]


def cmd_export_manifest(args):
    rows = read_tsv(args.master)
    selected = [row for row in rows if row_matches(row, args.require or [], args.require_empty or [])]
    if args.sort_by:
        selected.sort(key=lambda row: row.get(args.sort_by, ""))
    if args.limit:
        selected = selected[: args.limit]
    columns = _parse_columns_arg(args.columns)
    if columns:
        selected = [{column: row.get(column, "") for column in columns} for row in selected]
    write_tsv(args.output, selected)


def cmd_merge_json_dir(args):
    base_rows = read_tsv(args.master)
    updates = [read_json(path) for path in list_json_files(args.result_dir)]
    write_tsv(args.output, merge_rows(base_rows, updates, key=args.key))


def cmd_merge_tsv(args):
    base_rows = read_tsv(args.master)
    updates = read_tsv(args.update)
    write_tsv(args.output, merge_rows(base_rows, updates, key=args.key))


def cmd_count_rows(args):
    print(len(read_tsv(args.table)))


def cmd_print_row(args):
    row = read_tsv(args.table)[args.index]
    print("\t".join(str(row.get(column, "")) for column in args.columns))


def _load_sequence_from_fasta(path):
    seq = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def cmd_write_af_fastas(args):
    structure = _structure()
    rows = read_tsv(args.master)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.receptor_sequence:
        receptor_seq = args.receptor_sequence.strip()
    elif args.receptor_fasta:
        receptor_seq = _load_sequence_from_fasta(args.receptor_fasta)
    elif args.receptor_pdb:
        receptor_seq, _ = structure.chain_sequence(structure.parse_structure_atoms(args.receptor_pdb), args.receptor_chain)
    else:
        raise SystemExit("Provide --receptor-sequence, --receptor-fasta, or --receptor-pdb")

    results = []
    for row in rows:
        binder_seq = row.get(args.binder_sequence_column, "")
        if not binder_seq and row.get("source_pdb"):
            binder_seq, _ = structure.chain_sequence(structure.parse_structure_atoms(row["source_pdb"]), row.get("binder_chain", "B"))
        if not binder_seq:
            continue
        fasta_path = output_dir / f"{row['design_id']}.fa"
        sequence = f"{receptor_seq}:{binder_seq}" if args.multimer else binder_seq
        with open(fasta_path, "w", encoding="utf-8") as handle:
            handle.write(f">{row['design_id']}\n{sequence}\n")
        results.append(
            {
                "design_id": row["design_id"],
                "af_fasta": str(fasta_path),
                "af_mode": "multimer" if args.multimer else "monomer",
                "af_backend": args.backend,
            }
        )
    write_tsv(args.output, results)


def _parse_fasta_records(path):
    records = []
    header = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def _parse_header_metric(header, key):
    key_eq = f"{key}="
    for token in header.replace(",", " ").split():
        if token.startswith(key_eq):
            return token.split("=", 1)[1]
    return ""


def _find_mpnn_fasta(folder):
    for pattern in ("*.fa", "*.fasta", "seqs/*.fa", "seqs/*.fasta"):
        candidates = sorted(folder.glob(pattern))
        if candidates:
            return candidates[0]
    return None


def _extract_mpnn_binder_sequence(sequence):
    if "/" in sequence:
        parts = [part.strip() for part in sequence.split("/") if part.strip()]
        if parts:
            return parts[-1]
    return sequence.strip()


def cmd_parse_proteinmpnn(args):
    rows = read_tsv(args.master)
    results = []
    for row in rows:
        design_id = row["design_id"]
        out_dir = Path(args.mpnn_root) / design_id
        payload = {
            "design_id": design_id,
            "binder_sequence": "",
            "mpnn_fasta": "",
            "mpnn_score": "",
            "mpnn_global_score": "",
            "mpnn_temperature": "",
        }
        fasta_path = _find_mpnn_fasta(out_dir)
        if fasta_path is None:
            results.append(payload)
            continue

        records = _parse_fasta_records(fasta_path)
        chosen_header = ""
        chosen_sequence = ""
        for header, sequence in records:
            if "native" in header.lower():
                continue
            chosen_header = header
            chosen_sequence = sequence
            break
        if not chosen_sequence and records:
            chosen_header, chosen_sequence = records[0]

        payload.update(
            {
                "binder_sequence": _extract_mpnn_binder_sequence(chosen_sequence),
                "mpnn_fasta": str(fasta_path),
                "mpnn_score": _parse_header_metric(chosen_header, "score"),
                "mpnn_global_score": _parse_header_metric(chosen_header, "global_score"),
                "mpnn_temperature": _parse_header_metric(chosen_header, "T"),
            }
        )
        results.append(payload)
    write_tsv(args.output, results)


def _find_best_af_model(folder):
    patterns = [
        "ranked_0*.pdb",
        "ranked_0*.cif",
        "ranked_0*.mmcif",
        "model_0*.pdb",
        "model_0*.cif",
        "*.pdb",
        "*.cif",
        "*.mmcif",
    ]
    for pattern in patterns:
        candidates = sorted(folder.glob(pattern))
        if candidates:
            return candidates[0]
    return None


def cmd_parse_af(args):
    structure = _structure()
    rows = read_tsv(args.master)
    motif_positions = structure.parse_motif_positions(args.motif_positions)
    results = []

    for row in rows:
        design_id = row["design_id"]
        af_dir = Path(args.af_root) / design_id
        payload = {
            "design_id": design_id,
            "af_pass": 0,
            "af_reason": "missing_af_dir",
        }
        if not af_dir.exists():
            results.append(payload)
            continue

        model_path = _find_best_af_model(af_dir)
        if model_path is None:
            payload["af_reason"] = "missing_ranked_model"
            results.append(payload)
            continue

        reference_pdb = row.get(args.reference_column) or row.get(args.fallback_reference_column)
        if not reference_pdb:
            payload["af_reason"] = "missing_reference_pdb"
            results.append(payload)
            continue

        pred_atoms = structure.parse_structure_atoms(model_path)
        receptor_chain = row.get("receptor_chain", args.receptor_chain) or args.receptor_chain
        binder_chain = row.get("binder_chain", args.binder_chain) or args.binder_chain
        _, receptor_residues = structure.chain_sequence(pred_atoms, receptor_chain)
        _, binder_residues = structure.chain_sequence(pred_atoms, binder_chain)

        payload.update(
            {
                "af_backend": row.get("af_backend", ""),
                "af_model_path": str(model_path),
                "af_model_pdb": str(model_path),
                "af_iface_contacts": structure.atom_contact_count(
                    pred_atoms,
                    receptor_chain,
                    binder_chain,
                    cutoff=args.iface_cutoff,
                ),
            }
        )
        payload.update(
            structure.binder_metrics_after_receptor_alignment(
                reference_pdb,
                str(model_path),
                receptor_chain=receptor_chain,
                binder_chain=binder_chain,
                motif_positions=motif_positions,
            )
        )
        payload.update(
            structure.extract_af_json_metrics(
                af_dir,
                receptor_len=len(receptor_residues),
                binder_len=len(binder_residues),
            )
        )
        if payload["af_plddt"] is None:
            payload["af_plddt"] = structure.mean_plddt_from_pdb(str(model_path))

        reasons = []
        af_rmsd = coerce_float(payload.get("af_rmsd_binder"))
        af_plddt = coerce_float(payload.get("af_plddt"))
        af_iptm = coerce_float(payload.get("af_iptm"))
        af_pae = coerce_float(payload.get("af_pae_interface"))
        af_motif = coerce_float(payload.get("af_motif_rmsd"))

        if af_rmsd is None or af_rmsd > args.max_rmsd:
            reasons.append(f"rmsd={af_rmsd if af_rmsd is not None else 'NA'}")
        if af_plddt is None or af_plddt < args.min_plddt:
            reasons.append(f"plddt={af_plddt if af_plddt is not None else 'NA'}")
        if af_iptm is not None and af_iptm < args.min_iptm:
            reasons.append(f"iptm={af_iptm:.3f}")
        if af_pae is not None and af_pae > args.max_interface_pae:
            reasons.append(f"iface_pae={af_pae:.3f}")
        if af_motif is not None and af_motif > args.max_motif_rmsd:
            reasons.append(f"motif_rmsd={af_motif:.3f}")
        if int(payload["af_iface_contacts"]) < args.min_iface_contacts:
            reasons.append(f"iface_contacts={payload['af_iface_contacts']}")

        payload["af_pass"] = 0 if reasons else 1
        payload["af_reason"] = "PASS" if not reasons else "; ".join(reasons)
        results.append(payload)

    write_tsv(args.output, results)


def cmd_cluster_sequences(args):
    structure = _structure()
    rows = [row for row in read_tsv(args.master) if row.get(args.sequence_column, "")]
    rank_specs = args.rank_by or [
        "af_pass",
        "geometry_pass",
        "robustness_pass",
        "energy_pass",
        "ddg_pass",
        "af_plddt",
        "-af_rmsd_binder",
        "geometry_score",
    ]
    ranked = sorted(rows, key=lambda row: build_rank_key(row, rank_specs))
    assigned = set()
    clusters = []

    for row in ranked:
        if row["design_id"] in assigned:
            continue
        cluster = [row]
        assigned.add(row["design_id"])
        ref_seq = row[args.sequence_column]
        for other in ranked:
            if other["design_id"] in assigned:
                continue
            if structure.global_sequence_identity(ref_seq, other[args.sequence_column]) >= args.identity_threshold:
                cluster.append(other)
                assigned.add(other["design_id"])
        clusters.append(cluster)

    output_rows = []
    for idx, cluster in enumerate(clusters, start=1):
        cluster_id = f"cluster_{idx:04d}"
        cluster_sorted = sorted(cluster, key=lambda row: build_rank_key(row, rank_specs))
        for member_idx, row in enumerate(cluster_sorted):
            output_rows.append(
                {
                    "design_id": row["design_id"],
                    "cluster_id": cluster_id,
                    "cluster_size": len(cluster_sorted),
                    "representative": 1 if member_idx < args.representatives_per_cluster else 0,
                }
            )
    write_tsv(args.output, output_rows)


def cmd_pick_final(args):
    rows = read_tsv(args.master)
    selected = [row for row in rows if row_matches(row, args.require or [], [])]
    rank_specs = args.rank_by or [
        "representative",
        "af_pass",
        "geometry_pass",
        "robustness_pass",
        "energy_pass",
        "ddg_pass",
        "af_plddt",
        "-af_rmsd_binder",
        "geometry_score",
    ]
    selected = sorted(selected, key=lambda row: build_rank_key(row, rank_specs))

    counts = {}
    final_rows = []
    for row in selected:
        cluster_id = row.get("cluster_id") or row["design_id"]
        counts.setdefault(cluster_id, 0)
        if counts[cluster_id] >= args.per_cluster:
            continue
        counts[cluster_id] += 1
        final_rows.append({"design_id": row["design_id"], "cluster_id": cluster_id, "final_pick": 1})
        if len(final_rows) >= args.max_picks:
            break
    write_tsv(args.output, final_rows)


def cmd_summarize_master(args):
    rows = read_tsv(args.master)
    print(f"total\t{len(rows)}")
    for column in args.columns:
        passed = sum(1 for row in rows if truthy(row.get(column)))
        print(f"{column}\t{passed}")


def build_parser():
    parser = argparse.ArgumentParser(description="BinderDesign-CUILab master-table operations")
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser("scan-designs")
    scan.add_argument("--design-root", required=True)
    scan.add_argument("--output", required=True)
    scan.add_argument("--pattern", action="append", default=["*.pdb", "*.cif", "*.mmcif"])
    scan.add_argument("--binder-chain", default="B")
    scan.add_argument("--receptor-chain", default="A")
    scan.add_argument("--exclude-dir", action="append", default=["keepers", "logs", "results"])
    scan.add_argument("--id-format", default="{config_name}__seed{seed}__{stem}")
    scan.set_defaults(func=cmd_scan_designs)

    export_manifest = subparsers.add_parser("export-manifest")
    export_manifest.add_argument("--master", required=True)
    export_manifest.add_argument("--output", required=True)
    export_manifest.add_argument("--require", action="append")
    export_manifest.add_argument("--require-empty", action="append")
    export_manifest.add_argument("--columns")
    export_manifest.add_argument("--sort-by")
    export_manifest.add_argument("--limit", type=int)
    export_manifest.set_defaults(func=cmd_export_manifest)

    merge_json = subparsers.add_parser("merge-json-dir")
    merge_json.add_argument("--master", required=True)
    merge_json.add_argument("--result-dir", required=True)
    merge_json.add_argument("--output", required=True)
    merge_json.add_argument("--key", default="design_id")
    merge_json.set_defaults(func=cmd_merge_json_dir)

    merge_tsv = subparsers.add_parser("merge-tsv")
    merge_tsv.add_argument("--master", required=True)
    merge_tsv.add_argument("--update", required=True)
    merge_tsv.add_argument("--output", required=True)
    merge_tsv.add_argument("--key", default="design_id")
    merge_tsv.set_defaults(func=cmd_merge_tsv)

    count_rows = subparsers.add_parser("count-rows")
    count_rows.add_argument("--table", required=True)
    count_rows.set_defaults(func=cmd_count_rows)

    print_row = subparsers.add_parser("print-row")
    print_row.add_argument("--table", required=True)
    print_row.add_argument("--index", type=int, required=True)
    print_row.add_argument("--columns", nargs="+", required=True)
    print_row.set_defaults(func=cmd_print_row)

    af_fasta = subparsers.add_parser("write-af-fastas")
    af_fasta.add_argument("--master", required=True)
    af_fasta.add_argument("--output-dir", required=True)
    af_fasta.add_argument("--output", required=True)
    af_fasta.add_argument("--receptor-sequence")
    af_fasta.add_argument("--receptor-fasta")
    af_fasta.add_argument("--receptor-pdb")
    af_fasta.add_argument("--receptor-chain", default="A")
    af_fasta.add_argument("--binder-sequence-column", default="binder_sequence")
    af_fasta.add_argument("--multimer", action="store_true")
    af_fasta.add_argument("--backend", choices=["af2", "af3"], default="af2")
    af_fasta.set_defaults(func=cmd_write_af_fastas)

    mpnn = subparsers.add_parser("parse-proteinmpnn")
    mpnn.add_argument("--master", required=True)
    mpnn.add_argument("--mpnn-root", required=True)
    mpnn.add_argument("--output", required=True)
    mpnn.set_defaults(func=cmd_parse_proteinmpnn)

    parse_af = subparsers.add_parser("parse-af")
    parse_af.add_argument("--master", required=True)
    parse_af.add_argument("--af-root", required=True)
    parse_af.add_argument("--output", required=True)
    parse_af.add_argument("--reference-column", default="restored_pdb")
    parse_af.add_argument("--fallback-reference-column", default="source_pdb")
    parse_af.add_argument("--receptor-chain", default="A")
    parse_af.add_argument("--binder-chain", default="B")
    parse_af.add_argument("--motif-positions")
    parse_af.add_argument("--max-rmsd", type=float, default=3.0)
    parse_af.add_argument("--min-plddt", type=float, default=80.0)
    parse_af.add_argument("--min-iptm", type=float, default=0.0)
    parse_af.add_argument("--max-interface-pae", type=float, default=15.0)
    parse_af.add_argument("--max-motif-rmsd", type=float, default=1.5)
    parse_af.add_argument("--min-iface-contacts", type=int, default=20)
    parse_af.add_argument("--iface-cutoff", type=float, default=5.0)
    parse_af.set_defaults(func=cmd_parse_af)

    cluster = subparsers.add_parser("cluster-sequences")
    cluster.add_argument("--master", required=True)
    cluster.add_argument("--output", required=True)
    cluster.add_argument("--sequence-column", default="binder_sequence")
    cluster.add_argument("--identity-threshold", type=float, default=0.7)
    cluster.add_argument("--representatives-per-cluster", type=int, default=1)
    cluster.add_argument("--rank-by", nargs="*")
    cluster.set_defaults(func=cmd_cluster_sequences)

    pick = subparsers.add_parser("pick-final")
    pick.add_argument("--master", required=True)
    pick.add_argument("--output", required=True)
    pick.add_argument("--require", action="append")
    pick.add_argument("--per-cluster", type=int, default=1)
    pick.add_argument("--max-picks", type=int, default=20)
    pick.add_argument("--rank-by", nargs="*")
    pick.set_defaults(func=cmd_pick_final)

    summary = subparsers.add_parser("summarize-master")
    summary.add_argument("--master", required=True)
    summary.add_argument(
        "--columns",
        nargs="+",
        default=[
            "qc_pass",
            "af_pass",
            "geometry_pass",
            "energy_pass",
            "robustness_pass",
            "ddg_pass",
            "final_pick",
        ],
    )
    summary.set_defaults(func=cmd_summarize_master)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
