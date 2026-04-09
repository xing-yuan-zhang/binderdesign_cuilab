import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from .common import sanitize_token


def _parse_chain_residues(path, chain):
    from .structure import parse_structure_atoms

    residues = []
    seen = set()
    for atom in parse_structure_atoms(path):
        if atom["chain"] != chain:
            continue
        resi = atom["resi"]
        if resi not in seen:
            seen.add(resi)
            residues.append(resi)
    return residues


def _parse_binding_region(spec, default_chain):
    hotspots = []
    for chunk in spec.split(","):
        token = chunk.strip()
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            left_chain = left[0] if left and left[0].isalpha() else default_chain
            right_chain = right[0] if right and right[0].isalpha() else left_chain
            if left_chain != right_chain:
                raise ValueError(f"Cross-chain hotspot range is not supported: {token}")
            left_resi = int(left[1:] if left[0].isalpha() else left)
            right_resi = int(right[1:] if right[0].isalpha() else right)
            for resi in range(left_resi, right_resi + 1):
                hotspots.append(f"{left_chain}{resi}")
            continue
        chain = token[0] if token[0].isalpha() else default_chain
        resi = int(token[1:] if token[0].isalpha() else token)
        hotspots.append(f"{chain}{resi}")

    ordered = []
    seen = set()
    for hotspot in hotspots:
        if hotspot not in seen:
            seen.add(hotspot)
            ordered.append(hotspot)
    return ordered


def _compress_ranges(numbers):
    if not numbers:
        return []
    ranges = []
    start = numbers[0]
    end = numbers[0]
    for value in numbers[1:]:
        if value == end + 1:
            end = value
        else:
            ranges.append((start, end))
            start = value
            end = value
    ranges.append((start, end))
    return ranges


def _build_contig(residues, hotspots, target_chain, binder_min, binder_max, crop_padding, full_chain):
    if not residues:
        raise ValueError(f"No residues found for chain {target_chain}; target file must contain ATOM/HETATM coordinates")

    hotspot_numbers = [int(item[1:]) for item in hotspots]
    if full_chain:
        selected = residues
    else:
        index_by_resi = {resi: idx for idx, resi in enumerate(residues)}
        missing = [resi for resi in hotspot_numbers if resi not in index_by_resi]
        if missing:
            raise ValueError(f"Hotspots not found in target chain {target_chain}: {missing}")
        left = max(0, min(index_by_resi[resi] for resi in hotspot_numbers) - crop_padding)
        right = min(len(residues) - 1, max(index_by_resi[resi] for resi in hotspot_numbers) + crop_padding)
        selected = residues[left : right + 1]

    segments = [f"{target_chain}{start}-{end}" for start, end in _compress_ranges(selected)]
    return "[" + " ".join(segments) + f"/0 {binder_min}-{binder_max}]"


def _relative_paths():
    return {
        "inputs_dir": "inputs",
        "configs_dir": "configs",
        "tables_dir": "tables",
        "stage_tables_dir": "stage_tables",
        "manifests_dir": "manifests",
        "stage_results_dir": "stage_results",
        "logs_dir": "logs",
        "af_fastas_dir": "af_fastas",
        "rfdiffusion_out": "outputs/rfdiffusion",
        "restored_out": "outputs/restored",
        "mpnn_out": "outputs/mpnn",
        "af_out": "outputs/alphafold",
        "rosetta_interface_work": "outputs/rosetta/interface",
        "rosetta_fastrelax_work": "outputs/rosetta/fastrelax",
        "rosetta_ddg_work": "outputs/rosetta/ddg",
        "robustness_raw_dir": "stage_tables/robustness_raw",
    }


def _thresholds():
    return {
        "af_max_rmsd": 3.0,
        "af_min_plddt": 80.0,
        "af_min_iptm": 0.0,
        "af_max_interface_pae": 15.0,
        "af_max_motif_rmsd": 1.5,
        "af_min_iface_contacts": 20,
        "ddg_max_value": 0.0,
        "cluster_identity_threshold": 0.7,
        "final_max_picks": 20,
    }


def _resolve(project_root, relative_path):
    return project_root / relative_path


def _ensure_dirs(project_root, paths):
    for relative_path in paths.values():
        _resolve(project_root, relative_path).mkdir(parents=True, exist_ok=True)


def _write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _load_project(path):
    path = Path(path).resolve()
    with open(path, encoding="utf-8") as handle:
        project = json.load(handle)
    return path, project


def _make_project(args):
    workspace_root = Path(args.workspace_root).resolve()
    project_name = sanitize_token(args.name)
    project_root = workspace_root / project_name
    if project_root.exists() and not args.overwrite:
        raise SystemExit(f"Project already exists: {project_root}")

    paths = _relative_paths()
    _ensure_dirs(project_root, paths)

    target_copy = project_root / paths["inputs_dir"] / Path(args.target_pdb).name
    shutil.copy2(args.target_pdb, target_copy)

    hotspots = _parse_binding_region(args.binding_region, args.target_chain)
    residues = _parse_chain_residues(target_copy, args.target_chain)
    contig = _build_contig(
        residues=residues,
        hotspots=hotspots,
        target_chain=args.target_chain,
        binder_min=args.binder_length_min,
        binder_max=args.binder_length_max,
        crop_padding=args.crop_padding,
        full_chain=args.full_chain,
    )

    rfd_config = project_root / paths["configs_dir"] / "rfdiffusion.txt"
    with open(rfd_config, "w", encoding="utf-8") as handle:
        handle.write(f"CONTIG={contig}\n")
        handle.write("HOTSPOTS=[" + ",".join(f'"{item}"' for item in hotspots) + "]\n")

    project = {
        "project_name": project_name,
        "target_chain": args.target_chain,
        "target_pdb": (Path(paths["inputs_dir"]) / target_copy.name).as_posix(),
        "binding_region_input": args.binding_region,
        "hotspots": hotspots,
        "binder_length_min": args.binder_length_min,
        "binder_length_max": args.binder_length_max,
        "crop_padding": args.crop_padding,
        "full_chain": bool(args.full_chain),
        "contig": contig,
        "rfdiffusion_config": (Path(paths["configs_dir"]) / "rfdiffusion.txt").as_posix(),
        "num_seeds": args.num_seeds,
        "num_designs_per_seed": args.num_designs_per_seed,
        "mpnn_num_seq_per_target": args.mpnn_num_seq_per_target,
        "alphafold_backend": args.af_backend,
        "paths": paths,
        "thresholds": _thresholds(),
    }

    project_json = project_root / "project.json"
    _write_json(project_json, project)
    return project_json


def _repo_root():
    return Path(__file__).resolve().parents[1]


def _stage_paths(project_root):
    return {
        "master0": project_root / "tables" / "master.tsv",
        "master2": project_root / "tables" / "master.stage2.tsv",
        "master2b": project_root / "tables" / "master.stage2b.tsv",
        "master2c": project_root / "tables" / "master.stage2c.tsv",
        "master3": project_root / "tables" / "master.stage3.tsv",
        "master4": project_root / "tables" / "master.stage4.tsv",
        "master5": project_root / "tables" / "master.stage5.tsv",
        "master6": project_root / "tables" / "master.stage6.tsv",
        "master7": project_root / "tables" / "master.stage7.tsv",
        "master8": project_root / "tables" / "master.stage8.tsv",
        "master_final": project_root / "tables" / "master.final.tsv",
        "light_qc_manifest": project_root / "manifests" / "light_qc.tsv",
        "mpnn_manifest": project_root / "manifests" / "mpnn.tsv",
        "af_manifest": project_root / "manifests" / "alphafold.tsv",
        "geometry_manifest": project_root / "manifests" / "geometry.tsv",
        "energy_manifest": project_root / "manifests" / "energy.tsv",
        "robustness_manifest": project_root / "manifests" / "robustness.tsv",
        "ddg_manifest": project_root / "manifests" / "ddg.tsv",
        "mpnn_tsv": project_root / "stage_tables" / "mpnn.tsv",
        "af_fastas_tsv": project_root / "stage_tables" / "af_fastas.tsv",
        "af_metrics_tsv": project_root / "stage_tables" / "af_metrics.tsv",
        "robustness_tsv": project_root / "stage_tables" / "robustness.tsv",
        "clusters_tsv": project_root / "stage_tables" / "clusters.tsv",
        "final_picks_tsv": project_root / "stage_tables" / "final_picks.tsv",
    }


def _run_local(cmd, dry_run=False, label="local"):
    printable = " ".join(str(part) for part in cmd)
    if dry_run:
        print(f"[{label}] {printable}")
        return
    subprocess.run(cmd, check=True)


def _sbatch(account, partition, script, script_args, array=None, dry_run=False, label="job"):
    cmd = ["sbatch", "--parsable", "--wait", "-A", account, "-p", partition]
    if array:
        cmd.extend(["--array", array])
    cmd.append(str(script))
    cmd.extend(str(arg) for arg in script_args)
    if dry_run:
        print(f"[{label}] {' '.join(cmd)}")
        return "DRYRUN"
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = result.stdout.strip().split(";")[0].strip()
    print(f"[{label}] job_id={job_id}")
    return job_id


def _table_count(path):
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]
    return max(0, len(lines) - 1)


def cmd_init_project(args):
    print(_make_project(args))


def cmd_submit_marlowe(args):
    project_json, project = _load_project(args.project)
    project_root = project_json.parent
    paths = _stage_paths(project_root)

    account = args.account or os.environ.get("MARLOWE_ACCOUNT")
    partition = args.partition or os.environ.get("MARLOWE_PARTITION", "preempt")
    if not account:
        raise SystemExit("Set --account or MARLOWE_ACCOUNT before submitting")

    repo_root = _repo_root()
    python_exe = sys.executable or "python"
    rel_paths = project["paths"]
    thresholds = project["thresholds"]

    target_pdb = _resolve(project_root, project["target_pdb"])
    rfd_config = _resolve(project_root, project["rfdiffusion_config"])
    rfd_out = _resolve(project_root, rel_paths["rfdiffusion_out"])
    restored_out = _resolve(project_root, rel_paths["restored_out"])
    mpnn_out = _resolve(project_root, rel_paths["mpnn_out"])
    af_fastas_dir = _resolve(project_root, rel_paths["af_fastas_dir"])
    af_out = _resolve(project_root, rel_paths["af_out"])
    stage_results = _resolve(project_root, rel_paths["stage_results_dir"])
    light_qc_results = stage_results / "light_qc"
    geometry_results = stage_results / "geometry"
    energy_results = stage_results / "energy"
    ddg_results = stage_results / "ddg"
    rosetta_interface_work = _resolve(project_root, rel_paths["rosetta_interface_work"])
    rosetta_fastrelax_work = _resolve(project_root, rel_paths["rosetta_fastrelax_work"])
    rosetta_ddg_work = _resolve(project_root, rel_paths["rosetta_ddg_work"])
    robustness_raw_dir = _resolve(project_root, rel_paths["robustness_raw_dir"])

    def cli(*parts):
        return [python_exe, "-m", "binderdesign_cuilab.cli", "table", *[str(part) for part in parts]]

    def run_cli(label, *parts):
        _run_local(cli(*parts), dry_run=args.dry_run, label=label)

    def submit_array(label, script_name, manifest_path, extra_args):
        count = 1 if args.dry_run and not Path(manifest_path).exists() else _table_count(manifest_path)
        if count <= 0:
            raise SystemExit(f"No rows available for {label}: {manifest_path}")
        return _sbatch(
            account=account,
            partition=partition,
            script=repo_root / "scripts" / "marlowe" / script_name,
            script_args=extra_args,
            array=f"0-{count - 1}",
            dry_run=args.dry_run,
            label=label,
        )

    _sbatch(
        account=account,
        partition=partition,
        script=repo_root / "scripts" / "marlowe" / "run_rfdiffusion_array.sh",
        script_args=[target_pdb, rfd_config, rfd_out, int(project["num_designs_per_seed"])],
        array=f"0-{int(project['num_seeds']) - 1}",
        dry_run=args.dry_run,
        label="rfdiffusion",
    )

    run_cli("scan_designs", "scan-designs", "--design-root", rfd_out, "--output", paths["master0"])
    run_cli("export_light_qc", "export-manifest", "--master", paths["master0"], "--output", paths["light_qc_manifest"])

    submit_array(
        "light_qc",
        "run_light_qc_array.sh",
        paths["light_qc_manifest"],
        [paths["light_qc_manifest"], target_pdb, restored_out, light_qc_results],
    )

    run_cli("merge_qc", "merge-json-dir", "--master", paths["master0"], "--result-dir", light_qc_results, "--output", paths["master2"])
    run_cli("export_mpnn", "export-manifest", "--master", paths["master2"], "--output", paths["mpnn_manifest"], "--require", "qc_pass=1")

    submit_array("proteinmpnn", "run_proteinmpnn_array.sh", paths["mpnn_manifest"], [paths["mpnn_manifest"], mpnn_out])

    run_cli("parse_mpnn", "parse-proteinmpnn", "--master", paths["master2"], "--mpnn-root", mpnn_out, "--output", paths["mpnn_tsv"])
    run_cli("merge_mpnn", "merge-tsv", "--master", paths["master2"], "--update", paths["mpnn_tsv"], "--output", paths["master2b"])
    run_cli(
        "write_af_fastas",
        "write-af-fastas",
        "--master",
        paths["master2b"],
        "--output-dir",
        af_fastas_dir,
        "--output",
        paths["af_fastas_tsv"],
        "--receptor-pdb",
        target_pdb,
        "--receptor-chain",
        project["target_chain"],
        "--multimer",
        "--backend",
        project["alphafold_backend"],
    )
    run_cli("merge_af_fastas", "merge-tsv", "--master", paths["master2b"], "--update", paths["af_fastas_tsv"], "--output", paths["master2c"])
    run_cli(
        "export_af",
        "export-manifest",
        "--master",
        paths["master2c"],
        "--output",
        paths["af_manifest"],
        "--require",
        "qc_pass=1",
        "--require",
        "binder_sequence",
    )

    submit_array(
        "alphafold",
        "run_alphafold_array.sh",
        paths["af_manifest"],
        [paths["af_manifest"], project["alphafold_backend"], af_out],
    )

    run_cli(
        "parse_af",
        "parse-af",
        "--master",
        paths["master2c"],
        "--af-root",
        af_out,
        "--output",
        paths["af_metrics_tsv"],
        "--max-rmsd",
        thresholds["af_max_rmsd"],
        "--min-plddt",
        thresholds["af_min_plddt"],
        "--min-iptm",
        thresholds["af_min_iptm"],
        "--max-interface-pae",
        thresholds["af_max_interface_pae"],
        "--max-motif-rmsd",
        thresholds["af_max_motif_rmsd"],
        "--min-iface-contacts",
        thresholds["af_min_iface_contacts"],
    )
    run_cli("merge_af", "merge-tsv", "--master", paths["master2c"], "--update", paths["af_metrics_tsv"], "--output", paths["master3"])
    run_cli(
        "export_geometry",
        "export-manifest",
        "--master",
        paths["master3"],
        "--output",
        paths["geometry_manifest"],
        "--require",
        "qc_pass=1",
        "--require",
        "af_pass=1",
    )

    submit_array(
        "geometry",
        "run_geometry_array.sh",
        paths["geometry_manifest"],
        [
            paths["geometry_manifest"],
            geometry_results,
            ",".join(project["hotspots"]),
            project["binder_length_min"],
            project["binder_length_max"],
        ],
    )

    run_cli("merge_geometry", "merge-json-dir", "--master", paths["master3"], "--result-dir", geometry_results, "--output", paths["master4"])
    run_cli(
        "export_energy",
        "export-manifest",
        "--master",
        paths["master4"],
        "--output",
        paths["energy_manifest"],
        "--require",
        "qc_pass=1",
        "--require",
        "af_pass=1",
        "--require",
        "geometry_pass=1",
    )

    submit_array(
        "rosetta_interface",
        "run_rosetta_interface_array.sh",
        paths["energy_manifest"],
        [paths["energy_manifest"], rosetta_interface_work, energy_results],
    )

    run_cli("merge_energy", "merge-json-dir", "--master", paths["master4"], "--result-dir", energy_results, "--output", paths["master5"])
    run_cli(
        "export_robustness",
        "export-manifest",
        "--master",
        paths["master5"],
        "--output",
        paths["robustness_manifest"],
        "--require",
        "af_pass=1",
        "--require",
        "geometry_pass=1",
        "--require",
        "energy_pass=1",
    )

    submit_array(
        "rosetta_fastrelax",
        "run_rosetta_fastrelax_array.sh",
        paths["robustness_manifest"],
        [paths["robustness_manifest"], rosetta_fastrelax_work, robustness_raw_dir],
    )
    _sbatch(
        account=account,
        partition=partition,
        script=repo_root / "scripts" / "marlowe" / "aggregate_robustness.sh",
        script_args=[robustness_raw_dir, paths["robustness_tsv"]],
        dry_run=args.dry_run,
        label="aggregate_robustness",
    )

    run_cli("merge_robustness", "merge-tsv", "--master", paths["master5"], "--update", paths["robustness_tsv"], "--output", paths["master6"])
    run_cli(
        "export_ddg",
        "export-manifest",
        "--master",
        paths["master6"],
        "--output",
        paths["ddg_manifest"],
        "--require",
        "af_pass=1",
        "--require",
        "geometry_pass=1",
        "--require",
        "energy_pass=1",
        "--require",
        "robustness_pass=1",
    )

    submit_array("ddg", "run_rosetta_ddg_array.sh", paths["ddg_manifest"], [paths["ddg_manifest"], rosetta_ddg_work, ddg_results])
    run_cli("merge_ddg", "merge-json-dir", "--master", paths["master6"], "--result-dir", ddg_results, "--output", paths["master7"])
    run_cli(
        "cluster",
        "cluster-sequences",
        "--master",
        paths["master7"],
        "--output",
        paths["clusters_tsv"],
        "--identity-threshold",
        thresholds["cluster_identity_threshold"],
    )
    run_cli("merge_clusters", "merge-tsv", "--master", paths["master7"], "--update", paths["clusters_tsv"], "--output", paths["master8"])
    run_cli(
        "pick_final",
        "pick-final",
        "--master",
        paths["master8"],
        "--output",
        paths["final_picks_tsv"],
        "--require",
        "af_pass=1",
        "--require",
        "geometry_pass=1",
        "--require",
        "energy_pass=1",
        "--require",
        "robustness_pass=1",
        "--require",
        "ddg_pass=1",
        "--max-picks",
        thresholds["final_max_picks"],
    )
    run_cli("merge_final", "merge-tsv", "--master", paths["master8"], "--update", paths["final_picks_tsv"], "--output", paths["master_final"])
    run_cli("summary", "summarize-master", "--master", paths["master_final"])

    print(f"final_project_table={paths['master_final']}")
    print(f"final_picks={paths['final_picks_tsv']}")


def cmd_run_project(args):
    project_json = _make_project(args)
    print(f"project={project_json}")
    if args.submit_marlowe:
        submit_args = argparse.Namespace(
            project=str(project_json),
            account=args.account,
            partition=args.partition,
            dry_run=args.dry_run,
        )
        cmd_submit_marlowe(submit_args)


def build_parser():
    parser = argparse.ArgumentParser(description="BinderDesign-CUILab end-to-end project orchestration")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--name", required=True)
    common.add_argument("--target-pdb", required=True)
    common.add_argument("--target-chain", default="A")
    common.add_argument("--binding-region", required=True)
    common.add_argument("--binder-length-min", type=int, default=60)
    common.add_argument("--binder-length-max", type=int, default=80)
    common.add_argument("--crop-padding", type=int, default=20)
    common.add_argument("--full-chain", action="store_true")
    common.add_argument("--num-seeds", type=int, default=4)
    common.add_argument("--num-designs-per-seed", type=int, default=100)
    common.add_argument("--mpnn-num-seq-per-target", type=int, default=1)
    common.add_argument("--af-backend", choices=["af2", "af3"], default="af2")
    common.add_argument("--workspace-root", default="projects")
    common.add_argument("--overwrite", action="store_true")

    init = subparsers.add_parser("init-project", parents=[common])
    init.set_defaults(func=cmd_init_project)

    submit = subparsers.add_parser("submit-marlowe")
    submit.add_argument("--project", required=True)
    submit.add_argument("--account")
    submit.add_argument("--partition")
    submit.add_argument("--dry-run", action="store_true")
    submit.set_defaults(func=cmd_submit_marlowe)

    run = subparsers.add_parser("run-project", parents=[common])
    run.add_argument("--submit-marlowe", action="store_true")
    run.add_argument("--account")
    run.add_argument("--partition")
    run.add_argument("--dry-run", action="store_true")
    run.set_defaults(func=cmd_run_project)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
