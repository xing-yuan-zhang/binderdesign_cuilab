import csv
import json
import math
import re
from pathlib import Path


DEFAULT_MASTER_COLUMNS = [
    "design_id",
    "design_name",
    "config_name",
    "seed",
    "source_pdb",
    "binder_chain",
    "receptor_chain",
    "binder_sequence",
    "binder_len_scan",
    "restored_pdb",
    "qc_pass",
    "qc_reason",
    "qc_min_ab_dist",
    "qc_chain_breaks",
    "qc_max_ca_gap",
    "qc_binder_len",
    "geometry_pass",
    "geometry_reason",
    "geometry_score",
    "binder_len",
    "binder_rg",
    "iface_contacts",
    "hotspots_contacted",
    "min_ab_dist",
    "internal_ca_contacts",
    "min_neighbors_per_residue",
    "helix_count",
    "total_helix_res",
    "mpnn_fasta",
    "mpnn_score",
    "mpnn_global_score",
    "mpnn_temperature",
    "af_fasta",
    "af_mode",
    "af_backend",
    "af_pass",
    "af_reason",
    "af_model_path",
    "af_model_pdb",
    "af_plddt",
    "af_iptm",
    "af_ptm",
    "af_pae_interface",
    "af_iface_contacts",
    "af_rmsd_binder",
    "af_tm_binder",
    "af_motif_rmsd",
    "energy_pass",
    "energy_reason",
    "fa_rep",
    "interface_dg",
    "buried_unsat",
    "packstat",
    "total_score",
    "robustness_pass",
    "robustness_reason",
    "robustness_n",
    "dg_std",
    "score_std",
    "hotspot_consistency",
    "ddg_pass",
    "ddg_reason",
    "ddg_value",
    "cluster_id",
    "cluster_size",
    "representative",
    "md_pass",
    "md_reason",
    "final_pick",
]


TRUE_VALUES = {"1", "true", "t", "yes", "y", "pass", "passed", "ok"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "fail", "failed"}
SEED_RE = re.compile(r"^seed[_-]?(\d+)$", re.IGNORECASE)


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_tsv(path, rows, fieldnames=None):
    rows = list(rows)
    ensure_parent(path)
    if fieldnames is None:
        fieldnames = []
        seen = set()
        for name in DEFAULT_MASTER_COLUMNS:
            fieldnames.append(name)
            seen.add(name)
        for row in rows:
            for key in row.keys():
                if key not in seen:
                    fieldnames.append(key)
                    seen.add(key)

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path, payload):
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def sanitize_token(value):
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("._") or "item"


def infer_seed_from_parts(parts):
    for part in parts:
        match = SEED_RE.match(part)
        if match:
            return match.group(1)
    return ""


def infer_config_from_parts(parts):
    for idx, part in enumerate(parts):
        if SEED_RE.match(part) and idx > 0:
            return parts[idx - 1]
    return parts[-2] if len(parts) >= 2 else parts[-1]


def truthy(value):
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if not text:
        return False
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    try:
        return float(text) != 0.0
    except ValueError:
        return True


def coerce_float(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text or text.lower() in {"none", "nan", "na", "-"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_condition(expr):
    if "!=" in expr:
        key, value = expr.split("!=", 1)
        return key.strip(), "!=", value.strip()
    if "=" in expr:
        key, value = expr.split("=", 1)
        return key.strip(), "=", value.strip()
    return expr.strip(), "truthy", ""


def match_condition(row, expr):
    key, op, expected = parse_condition(expr)
    actual = row.get(key, "")
    if op == "truthy":
        return truthy(actual)
    if op == "=":
        return str(actual) == expected
    if op == "!=":
        return str(actual) != expected
    raise ValueError(f"Unsupported condition: {expr}")


def row_matches(row, conditions=None, require_empty=None):
    conditions = conditions or []
    require_empty = require_empty or []
    for expr in conditions:
        if not match_condition(row, expr):
            return False
    for key in require_empty:
        if str(row.get(key, "")).strip():
            return False
    return True


def merge_rows(base_rows, update_rows, key="design_id"):
    base_rows = [dict(row) for row in base_rows]
    index = {row.get(key, ""): row for row in base_rows}
    for update in update_rows:
        token = update.get(key, "")
        if token in index:
            index[token].update(update)
    return base_rows


def list_json_files(folder):
    folder = Path(folder)
    if not folder.exists():
        return []
    return sorted(str(path) for path in folder.rglob("*.json"))


def numeric_sort_value(value, descending):
    number = coerce_float(value)
    if number is None:
        return math.inf if descending else -math.inf
    return -number if descending else number


def text_sort_value(value, descending):
    text = str(value or "").lower()
    if descending:
        return "".join(chr(255 - ord(char)) for char in text)
    return text


def build_rank_key(row, specs):
    key = []
    for spec in specs:
        descending = True
        column = spec
        if spec.startswith("-"):
            descending = False
            column = spec[1:]
        value = row.get(column, "")
        number = coerce_float(value)
        if number is not None:
            key.append(numeric_sort_value(number, descending))
            continue
        if str(value).strip().lower() in TRUE_VALUES.union(FALSE_VALUES):
            bool_number = 1.0 if truthy(value) else 0.0
            key.append(numeric_sort_value(bool_number, descending))
            continue
        key.append(text_sort_value(value, descending))
    return tuple(key)
