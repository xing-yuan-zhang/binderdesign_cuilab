import os
from pathlib import Path
import statistics
import subprocess

import numpy as np

from .common import coerce_float, read_tsv, write_json, write_tsv
from .geometry import GeometryConfig, compute_score, parse_hotspots, triage_design
from .structure import (
    apply_rt,
    chain_break_metrics,
    kabsch,
    min_interchain_distance,
    parse_structure_atoms,
    residues_for_chain,
    rmsd,
)


def _repo_root():
    return Path(__file__).resolve().parents[1]


def _row(manifest, index):
    rows = read_tsv(manifest)
    return rows[index]


def _require_env(name):
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"Missing environment variable: {name}")
    return value


def _parse_config_value(path, key):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(f"{key}="):
                return line.strip().split("=", 1)[1]
    return ""


def _apptainer_binds():
    binds = []
    for token in [
        os.environ.get("ROOT"),
        os.environ.get("MARLOWE_SCRATCH"),
        str(_repo_root()),
        str(Path.cwd()),
    ]:
        if token and token not in binds and Path(token).exists():
            binds.append(token)
    extra = os.environ.get("BINDERDESIGN_CUILAB_APPTAINER_BIND", "")
    for token in extra.split(","):
        token = token.strip()
        if token and token not in binds:
            binds.append(token)
    return binds


def _containerize_if_needed(command, sif_envs, gpu=False):
    sif = ""
    for env_name in sif_envs:
        sif = os.environ.get(env_name, "").strip()
        if sif:
            break
    if not sif:
        return [str(part) for part in command]

    wrapped = ["apptainer", "exec"]
    if gpu:
        wrapped.append("--nv")
    for bind in _apptainer_binds():
        wrapped.extend(["--bind", bind])
    wrapped.append(sif)
    wrapped.extend(str(part) for part in command)
    return wrapped


def _run(cmd):
    subprocess.run([str(part) for part in cmd], check=True)


def _write_pdb_lines(path, lines):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def _ca_by_resi(atoms, chain):
    mapping = {}
    for atom in atoms:
        if atom["chain"] == chain and atom["atom"] == "CA":
            mapping[atom["resi"]] = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
    return mapping


def _restore_complex(template_pdb, design_pdb, out_pdb, receptor_chain="A", binder_chain="B"):
    if Path(template_pdb).suffix.lower() not in {".pdb", ".ent"}:
        return False, "template target must be a PDB file for restore"
    if Path(design_pdb).suffix.lower() not in {".pdb", ".ent"}:
        return False, "design structure must be a PDB file for restore"

    template_atoms = parse_structure_atoms(template_pdb)
    design_atoms = parse_structure_atoms(design_pdb)
    if not any(atom["chain"] == binder_chain for atom in design_atoms):
        return False, f"missing binder chain {binder_chain}"

    template_ca = _ca_by_resi(template_atoms, receptor_chain)
    design_ca = _ca_by_resi(design_atoms, receptor_chain)
    common = sorted(set(template_ca).intersection(design_ca))
    if len(common) < 3:
        return False, f"insufficient common CA residues: {len(common)}"

    design_coords = np.array([design_ca[resi] for resi in common], dtype=float)
    template_coords = np.array([template_ca[resi] for resi in common], dtype=float)
    rotation, translation = kabsch(design_coords, template_coords)

    out_lines = []
    serial = 1
    with open(template_pdb, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                break
            if line.startswith("MODEL"):
                out_lines.append(line.rstrip("\n"))

    for atom in template_atoms:
        if atom["chain"] != receptor_chain or not atom["line"]:
            continue
        line = list(atom["line"])
        line[6:11] = list(f"{serial:5d}")
        line[30:54] = list(f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}")
        out_lines.append("".join(line))
        serial += 1

    for atom in design_atoms:
        if atom["chain"] != binder_chain or not atom["line"]:
            continue
        transformed = rotation @ np.array([atom["x"], atom["y"], atom["z"]], dtype=float) + translation
        line = list(atom["line"])
        line[6:11] = list(f"{serial:5d}")
        line[21] = binder_chain
        line[30:54] = list(f"{transformed[0]:8.3f}{transformed[1]:8.3f}{transformed[2]:8.3f}")
        out_lines.append("".join(line))
        serial += 1

    out_lines.append("END")
    _write_pdb_lines(out_pdb, out_lines)
    aligned = apply_rt(design_coords, rotation, translation)
    return True, f"restored_rmsd={rmsd(aligned, template_coords):.3f}"


def _qc_metrics_for_path(path, receptor_chain, binder_chain, min_dist=2.0, max_gap=4.5):
    atoms = parse_structure_atoms(path)
    binder_residues = residues_for_chain(atoms, binder_chain)
    binder_len = len(binder_residues)
    min_ab = min_interchain_distance(atoms, receptor_chain, binder_chain)
    chain_breaks, max_ca_gap = chain_break_metrics(atoms, binder_chain, max_ca_gap=max_gap)
    reasons = []

    if binder_len == 0:
        reasons.append(f"no binder chain {binder_chain}")
    if min_ab is None:
        reasons.append("missing receptor/binder atoms")
    elif min_ab < min_dist:
        reasons.append(f"min_AB_dist={min_ab:.2f} < {min_dist:.2f}")
    if chain_breaks > 0:
        reasons.append(f"chain_breaks={chain_breaks}")

    return {
        "qc_pass": 0 if reasons else 1,
        "qc_reason": "PASS" if not reasons else "; ".join(reasons),
        "qc_min_ab_dist": None if min_ab is None else round(min_ab, 3),
        "qc_chain_breaks": chain_breaks,
        "qc_max_ca_gap": round(max_ca_gap, 3),
        "qc_binder_len": binder_len,
    }


def _rosetta_binary(rosetta_bin, *candidates):
    for candidate in candidates:
        path = Path(rosetta_bin) / candidate
        if path.exists():
            return path
    return Path(rosetta_bin) / candidates[0]


def _run_rosetta(command):
    _run(_containerize_if_needed(command, ("ROSETTA_SIF",), gpu=False))


def run_rfdiffusion_task(target_pdb, config_file, output_root, num_designs, index):
    python_exe = _require_env("RFDIFFUSION_PYTHON")
    script = _require_env("RFDIFFUSION_SCRIPT")
    contig = _parse_config_value(config_file, "CONTIG")
    hotspots = _parse_config_value(config_file, "HOTSPOTS")

    outdir = Path(output_root) / f"seed_{index}"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        python_exe,
        script,
        "--config-name",
        "base",
        f"inference.input_pdb={target_pdb}",
        f"inference.output_prefix={outdir / 'design'}",
        f"inference.num_designs={num_designs}",
        f"contigmap.contigs={contig}",
        f"+ppi.hotspot_res={hotspots}",
        f"+seed={index}",
    ]
    if os.environ.get("RFDIFFUSION_MODEL_DIR"):
        cmd.append(f"inference.model_directory_path={os.environ['RFDIFFUSION_MODEL_DIR']}")
    if os.environ.get("RFDIFFUSION_CKPT_PATH"):
        cmd.append(f"+inference.ckpt_path={os.environ['RFDIFFUSION_CKPT_PATH']}")
    if os.environ.get("RFDIFFUSION_EXTRA_ARGS"):
        cmd.extend(os.environ["RFDIFFUSION_EXTRA_ARGS"].split())
    _run(_containerize_if_needed(cmd, ("RFDIFFUSION_SIF",), gpu=True))


def run_light_qc_task(manifest, index, template_pdb, output_root, result_dir):
    row = _row(manifest, index)
    design_id = row["design_id"]
    source_pdb = row["source_pdb"]
    output_root = Path(output_root)
    result_dir = Path(result_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    restored_pdb = output_root / f"{design_id}_restored.pdb"

    ok, reason = _restore_complex(
        template_pdb,
        source_pdb,
        restored_pdb,
        row.get("receptor_chain", "A"),
        row.get("binder_chain", "B"),
    )
    if not ok:
        payload = {"design_id": design_id, "restored_pdb": str(restored_pdb), "qc_pass": 0, "qc_reason": reason}
    else:
        payload = {"design_id": design_id, "restored_pdb": str(restored_pdb)}
        payload.update(
            _qc_metrics_for_path(
                str(restored_pdb),
                receptor_chain=row.get("receptor_chain", "A"),
                binder_chain=row.get("binder_chain", "B"),
                min_dist=float(os.environ.get("QC_MIN_AB_DIST", "2.0")),
                max_gap=float(os.environ.get("QC_MAX_CA_GAP", "4.5")),
            )
        )
    write_json(result_dir / f"{design_id}.json", payload)


def run_geometry_task(manifest, index, result_dir, hotspots, binder_length_min, binder_length_max, pdb_column="restored_pdb"):
    row = _row(manifest, index)
    design_id = row["design_id"]
    hotspot_list = parse_hotspots(hotspots)
    config = GeometryConfig(
        receptor_chain=row.get("receptor_chain", "A") or "A",
        binder_chain=row.get("binder_chain", "B") or "B",
        hotspots=hotspot_list,
        min_hotspots_contacted=min(int(os.environ.get("GEOMETRY_MIN_HOTSPOTS_CONTACTED", "3")), len(hotspot_list))
        if hotspot_list
        else 0,
        binder_len_min=int(binder_length_min),
        binder_len_max=int(binder_length_max),
        min_iface_contacts=int(os.environ.get("GEOMETRY_MIN_IFACE_CONTACTS", "40")),
        hotspot_cutoff=float(os.environ.get("GEOMETRY_HOTSPOT_CUTOFF", "8.0")),
        iface_cutoff=float(os.environ.get("GEOMETRY_IFACE_CUTOFF", "5.0")),
        min_allowed_dist=float(os.environ.get("GEOMETRY_MIN_ALLOWED_DIST", "2.0")),
        rg_min=float(os.environ.get("GEOMETRY_RG_MIN", "7.0")),
        rg_max=float(os.environ.get("GEOMETRY_RG_MAX", "20.0")),
        min_internal_ca_contacts=int(os.environ.get("GEOMETRY_MIN_INTERNAL_CA_CONTACTS", "80")),
        min_neighbors_per_residue=int(os.environ.get("GEOMETRY_MIN_NEIGHBORS", "2")),
        use_helix_filter=os.environ.get("GEOMETRY_USE_HELIX_FILTER", "1") not in {"0", "false", "False"},
        min_helices=int(os.environ.get("GEOMETRY_MIN_HELICES", "3")),
        min_helix_len=int(os.environ.get("GEOMETRY_MIN_HELIX_LEN", "5")),
    )
    ok, reason, metrics = triage_design(row[pdb_column], config)
    payload = {"design_id": design_id, "geometry_pass": 1 if ok else 0, "geometry_reason": reason}
    if metrics:
        payload.update(
            {
                "geometry_score": round(compute_score(metrics), 3),
                "binder_len": metrics["binder_len"],
                "binder_rg": round(metrics["Rg"], 3),
                "iface_contacts": metrics["iface_contacts"],
                "hotspots_contacted": metrics["hotspots_contacted"],
                "min_ab_dist": round(metrics["min_AB_dist"], 3),
                "internal_ca_contacts": metrics["internal_CA_contacts"],
                "min_neighbors_per_residue": metrics["min_neighbors_per_residue"],
                "helix_count": metrics["helix_count"],
                "total_helix_res": metrics["total_helix_res"],
            }
        )
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    write_json(result_dir / f"{design_id}.json", payload)


def run_proteinmpnn_task(manifest, index, output_root, pdb_column="restored_pdb"):
    python_exe = _require_env("PROTEINMPNN_PYTHON")
    script = _require_env("PROTEINMPNN_SCRIPT")
    row = _row(manifest, index)
    design_id = row["design_id"]
    pdb_path = row[pdb_column]
    outdir = Path(output_root) / design_id
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        python_exe,
        script,
        "--pdb_path",
        pdb_path,
        "--out_folder",
        str(outdir),
        "--num_seq_per_target",
        os.environ.get("MPNN_NUM_SEQ_PER_TARGET", "1"),
        "--sampling_temp",
        os.environ.get("MPNN_SAMPLING_TEMP", "0.1"),
        "--batch_size",
        os.environ.get("MPNN_BATCH_SIZE", "1"),
        "--pdb_path_chains",
        os.environ.get("MPNN_DESIGN_CHAINS", "B"),
    ]
    if os.environ.get("MPNN_EXTRA_ARGS"):
        cmd.extend(os.environ["MPNN_EXTRA_ARGS"].split())
    _run(_containerize_if_needed(cmd, ("PROTEINMPNN_SIF",), gpu=False))


def run_alphafold_task(manifest, index, output_root, backend=None):
    row = _row(manifest, index)
    design_id = row["design_id"]
    fasta = row["af_fasta"]
    outdir = Path(output_root) / design_id
    outdir.mkdir(parents=True, exist_ok=True)

    backend = (backend or row.get("af_backend") or os.environ.get("AF_BACKEND") or "af2").strip().lower()
    if backend == "af2":
        batch_cmd = os.environ.get("AF2_BATCH_CMD") or os.environ.get("COLABFOLD_CMD") or "colabfold_batch"
        cmd = [
            batch_cmd,
            "--model-type",
            os.environ.get("AF2_MODEL_TYPE", os.environ.get("COLABFOLD_MODEL_TYPE", "alphafold2_multimer_v3")),
        ]
        extra = os.environ.get("AF2_EXTRA_ARGS") or os.environ.get("COLABFOLD_EXTRA_FLAGS")
        if extra:
            cmd.extend(extra.split())
        cmd.extend([fasta, str(outdir)])
        _run(_containerize_if_needed(cmd, ("AF2_SIF", "COLABFOLD_SIF"), gpu=True))
        return

    if backend == "af3":
        af3_cmd = os.environ.get("AF3_CMD", "").strip()
        if af3_cmd:
            cmd = [af3_cmd]
        else:
            cmd = [_require_env("AF3_PYTHON"), _require_env("AF3_SCRIPT")]
        if os.environ.get("AF3_MODEL_TYPE"):
            cmd.extend(["--model-type", os.environ["AF3_MODEL_TYPE"]])
        if os.environ.get("AF3_EXTRA_ARGS"):
            cmd.extend(os.environ["AF3_EXTRA_ARGS"].split())
        cmd.extend([fasta, str(outdir)])
        _run(_containerize_if_needed(cmd, ("AF3_SIF",), gpu=True))
        return

    raise SystemExit(f"Unsupported AlphaFold backend: {backend}")


def run_af2_task(manifest, index, output_root):
    run_alphafold_task(manifest, index, output_root, backend="af2")


def _parse_score_sc(path):
    header = None
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line.startswith("SCORE:"):
                continue
            parts = line.split()
            if header is None and "description" in parts:
                header = parts[1:]
                continue
            if header is None:
                continue
            values = parts[1:]
            if len(values) != len(header):
                continue
            rows.append(dict(zip(header, values)))
    return rows


def _score_value(row, names):
    for name in names:
        value = row.get(name, "")
        if value not in {"", "NA", "nan"}:
            return value
    return ""


def _interface_payload(design_id, row):
    return {
        "design_id": design_id,
        "fa_rep": _score_value(row, ["fa_rep"]),
        "interface_dg": _score_value(row, ["dG_separated", "interface_dG", "interface_delta"]),
        "buried_unsat": _score_value(row, ["buried_unsat", "buried_unsat_hbonds", "buns_all_heavy"]),
        "packstat": _score_value(row, ["packstat"]),
        "total_score": _score_value(row, ["total_score", "score"]),
    }


def _apply_interface_thresholds(payload):
    reasons = []
    fa_rep = coerce_float(payload.get("fa_rep"))
    interface_dg = coerce_float(payload.get("interface_dg"))
    buried_unsat = coerce_float(payload.get("buried_unsat"))
    packstat = coerce_float(payload.get("packstat"))

    if fa_rep is not None and fa_rep > float(os.environ.get("ENERGY_MAX_FA_REP", "100")):
        reasons.append(f"fa_rep={fa_rep:.3f}")
    if buried_unsat is not None and buried_unsat > float(os.environ.get("ENERGY_MAX_BURIED_UNSAT", "10")):
        reasons.append(f"buried_unsat={buried_unsat:.3f}")
    if packstat is not None and packstat < float(os.environ.get("ENERGY_MIN_PACKSTAT", "0.5")):
        reasons.append(f"packstat={packstat:.3f}")
    if interface_dg is not None and interface_dg > float(os.environ.get("ENERGY_MAX_INTERFACE_DG", "0.0")):
        reasons.append(f"interface_dg={interface_dg:.3f}")

    payload["energy_pass"] = 0 if reasons else 1
    payload["energy_reason"] = "PASS" if not reasons else "; ".join(reasons)
    return payload


def run_rosetta_interface_task(manifest, index, work_root, result_dir, pdb_column="restored_pdb"):
    rosetta_bin = Path(_require_env("ROSETTA_BIN"))
    rosetta_db = Path(_require_env("ROSETTA_DB"))
    row = _row(manifest, index)
    design_id = row["design_id"]
    pdb_path = row[pdb_column]
    workdir = Path(work_root) / design_id
    workdir.mkdir(parents=True, exist_ok=True)

    score_bin = _rosetta_binary(rosetta_bin, "score_jd2.default.linuxgccrelease", "score_jd2.static.linuxgccrelease")
    ia_bin = _rosetta_binary(rosetta_bin, "InterfaceAnalyzer.default.linuxgccrelease", "InterfaceAnalyzer.static.linuxgccrelease")
    score_sc = workdir / "score.sc"
    ia_sc = workdir / "interface.sc"

    _run_rosetta([score_bin, "-database", rosetta_db, "-s", pdb_path, "-out:path:all", workdir, "-out:file:scorefile", score_sc])
    _run_rosetta(
        [
            ia_bin,
            "-database",
            rosetta_db,
            "-s",
            pdb_path,
            "-interface",
            os.environ.get("ROSETTA_INTERFACE", "A_B"),
            "-pack_input",
            "true",
            "-pack_separated",
            "true",
            "-compute_packstat",
            "true",
            "-out:path:all",
            workdir,
            "-out:file:scorefile",
            ia_sc,
        ]
    )

    merged = {}
    for scorefile in (score_sc, ia_sc):
        rows = _parse_score_sc(scorefile)
        if rows:
            merged.update(rows[-1])
    payload = _apply_interface_thresholds(_interface_payload(design_id, merged))
    write_json(Path(result_dir) / f"{design_id}.json", payload)


def _rep_seed(rep):
    return 1001 + rep * 1001


def run_rosetta_fastrelax_task(manifest, index, work_root, raw_dir, pdb_column="restored_pdb"):
    rosetta_bin = Path(_require_env("ROSETTA_BIN"))
    rosetta_db = Path(_require_env("ROSETTA_DB"))
    row = _row(manifest, index)
    design_id = row["design_id"]
    pdb_path = row[pdb_column]
    nstruct = int(os.environ.get("ROBUSTNESS_NSTRUCT", "8"))

    workdir = Path(work_root) / design_id / "fastrelax"
    workdir.mkdir(parents=True, exist_ok=True)
    raw_tsv = Path(raw_dir) / f"{design_id}.tsv"
    raw_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(raw_tsv, "w", encoding="utf-8") as handle:
        handle.write("design_id\treplicate_id\ttotal_score\tinterface_dg\n")

    rosetta_scripts = Path(os.environ.get("ROSETTA_ROSETTASCRIPTS_BIN", _rosetta_binary(rosetta_bin, "rosetta_scripts.default.linuxgccrelease", "rosetta_scripts.static.linuxgccrelease")))
    relax_xml = Path(os.environ.get("ROSETTA_RELAX_XML", _repo_root() / "rosetta" / "relax_interface.xml"))
    relax_bin = Path(os.environ.get("ROSETTA_FASTRELAX_BIN", _rosetta_binary(rosetta_bin, "relax.default.linuxgccrelease", "relax.static.linuxgccrelease")))
    ia_bin = _rosetta_binary(rosetta_bin, "InterfaceAnalyzer.default.linuxgccrelease", "InterfaceAnalyzer.static.linuxgccrelease")

    for rep in range(nstruct):
        repdir = workdir / f"rep{rep}"
        repdir.mkdir(parents=True, exist_ok=True)
        score_sc = repdir / "score.sc"
        if rosetta_scripts.exists() and relax_xml.exists():
            _run_rosetta(
                [
                    rosetta_scripts,
                    "-database",
                    rosetta_db,
                    "-parser:protocol",
                    relax_xml,
                    "-in:file:s",
                    pdb_path,
                    "-ignore_unrecognized_res",
                    "true",
                    "-run:constant_seed",
                    "-run:jran",
                    _rep_seed(rep),
                    "-nstruct",
                    "1",
                    "-out:path:all",
                    repdir,
                    "-out:file:scorefile",
                    score_sc,
                ]
            )
        else:
            _run_rosetta(
                [
                    relax_bin,
                    "-database",
                    rosetta_db,
                    "-s",
                    pdb_path,
                    "-nstruct",
                    "1",
                    "-run:constant_seed",
                    "-run:jran",
                    _rep_seed(rep),
                    "-relax:constrain_relax_to_start_coords",
                    "-relax:coord_constrain_sidechains",
                    "-out:path:all",
                    repdir,
                    "-out:file:scorefile",
                    score_sc,
                ]
            )

        score_rows = _parse_score_sc(score_sc)
        total_score = ""
        if score_rows:
            total_score = _score_value(score_rows[-1], ["total_score", "score"])

        model_pdbs = sorted(repdir.glob("*.pdb"))
        interface_dg = ""
        if model_pdbs:
            ia_sc = repdir / "interface.sc"
            _run_rosetta(
                [
                    ia_bin,
                    "-database",
                    rosetta_db,
                    "-s",
                    model_pdbs[0],
                    "-interface",
                    os.environ.get("ROSETTA_INTERFACE", "A_B"),
                    "-pack_input",
                    "true",
                    "-pack_separated",
                    "true",
                    "-compute_packstat",
                    "true",
                    "-out:path:all",
                    repdir,
                    "-out:file:scorefile",
                    ia_sc,
                ]
            )
            ia_rows = _parse_score_sc(ia_sc)
            if ia_rows:
                interface_dg = _score_value(ia_rows[-1], ["dG_separated", "interface_dG", "interface_delta"])

        with open(raw_tsv, "a", encoding="utf-8") as handle:
            handle.write(f"{design_id}\trep{rep}\t{total_score}\t{interface_dg}\n")


def run_rosetta_ddg_task(manifest, index, work_root, result_dir, pdb_column="restored_pdb"):
    rosetta_bin = Path(_require_env("ROSETTA_BIN"))
    rosetta_db = Path(_require_env("ROSETTA_DB"))
    row = _row(manifest, index)
    design_id = row["design_id"]
    pdb_path = row[pdb_column]
    workdir = Path(work_root) / design_id / "ddg"
    workdir.mkdir(parents=True, exist_ok=True)

    ddg_bin = os.environ.get("ROSETTA_DDG_BIN", "").strip()
    ddg_sc = workdir / "ddg.sc"
    ddg_value = None

    if ddg_bin:
        command = [ddg_bin]
        if os.environ.get("ROSETTA_DDG_EXTRA_ARGS"):
            command.extend(os.environ["ROSETTA_DDG_EXTRA_ARGS"].split())
        command.extend(["-database", rosetta_db, "-s", pdb_path, "-out:path:all", workdir, "-out:file:scorefile", ddg_sc])
        _run_rosetta(command)
        rows = _parse_score_sc(ddg_sc)
        if rows:
            ddg_value = coerce_float(_score_value(rows[-1], ["ddg", "interface_ddg", "total_score", "score"]))
    else:
        ia_bin = _rosetta_binary(rosetta_bin, "InterfaceAnalyzer.default.linuxgccrelease", "InterfaceAnalyzer.static.linuxgccrelease")
        _run_rosetta(
            [
                ia_bin,
                "-database",
                rosetta_db,
                "-s",
                pdb_path,
                "-interface",
                os.environ.get("ROSETTA_INTERFACE", "A_B"),
                "-pack_input",
                "true",
                "-pack_separated",
                "true",
                "-compute_packstat",
                "true",
                "-out:path:all",
                workdir,
                "-out:file:scorefile",
                ddg_sc,
            ]
        )
        rows = _parse_score_sc(ddg_sc)
        if rows:
            ddg_value = coerce_float(_score_value(rows[-1], ["dG_separated", "interface_dG", "interface_delta"]))

    max_ddg = float(os.environ.get("DDG_MAX_VALUE", "0.0"))
    reason = "PASS" if ddg_value is not None and ddg_value <= max_ddg else f"ddg={ddg_value}"
    if not ddg_bin and ddg_value is not None:
        reason = f"{reason}; fallback=InterfaceAnalyzer" if reason != "PASS" else "PASS; fallback=InterfaceAnalyzer"
    payload = {
        "design_id": design_id,
        "ddg_value": "" if ddg_value is None else round(ddg_value, 4),
        "ddg_pass": 1 if ddg_value is not None and ddg_value <= max_ddg else 0,
        "ddg_reason": reason,
    }
    write_json(Path(result_dir) / f"{design_id}.json", payload)


def aggregate_robustness(raw_dir, output_tsv):
    rows = []
    for path in sorted(Path(raw_dir).glob("*.tsv")):
        rows.extend(read_tsv(path))

    grouped = {}
    for row in rows:
        grouped.setdefault(row["design_id"], []).append(row)

    summary = []
    for design_id, members in grouped.items():
        dg_values = [coerce_float(member.get("interface_dg")) for member in members]
        dg_values = [value for value in dg_values if value is not None]
        score_values = [coerce_float(member.get("total_score")) for member in members]
        score_values = [value for value in score_values if value is not None]

        dg_std = statistics.pstdev(dg_values) if len(dg_values) >= 2 else None
        score_std = statistics.pstdev(score_values) if len(score_values) >= 2 else None
        reasons = []
        if len(members) < int(os.environ.get("ROBUSTNESS_MIN_REPLICATES", "5")):
            reasons.append(f"replicates={len(members)}")
        if dg_std is not None and dg_std > float(os.environ.get("ROBUSTNESS_MAX_DG_STD", "2.0")):
            reasons.append(f"dg_std={dg_std:.3f}")
        if score_std is not None and score_std > float(os.environ.get("ROBUSTNESS_MAX_SCORE_STD", "5.0")):
            reasons.append(f"score_std={score_std:.3f}")
        summary.append(
            {
                "design_id": design_id,
                "robustness_n": len(members),
                "dg_std": "" if dg_std is None else round(dg_std, 4),
                "score_std": "" if score_std is None else round(score_std, 4),
                "robustness_pass": 0 if reasons else 1,
                "robustness_reason": "PASS" if not reasons else "; ".join(reasons),
            }
        )

    write_tsv(output_tsv, summary)
