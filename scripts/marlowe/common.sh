#!/usr/bin/env bash

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export ROOT

: "${MARLOWE_PROJECT_ID:?Set MARLOWE_PROJECT_ID, for example m223813}"

export MARLOWE_PARTITION="${MARLOWE_PARTITION:-preempt}"
export MARLOWE_SCRATCH="${MARLOWE_SCRATCH:-/scratch/${MARLOWE_PROJECT_ID}/${USER}/mm-pbsa}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${MARLOWE_SCRATCH}/apptainercache}"
export HF_HOME="${HF_HOME:-${MARLOWE_SCRATCH}/cache/hf}"
export TORCH_HOME="${TORCH_HOME:-${MARLOWE_SCRATCH}/cache/torch}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${MARLOWE_SCRATCH}/cache/xdg}"
export TMPDIR="${TMPDIR:-${MARLOWE_SCRATCH}/tmp}"

mkdir -p \
  "$APPTAINER_CACHEDIR" \
  "$HF_HOME" \
  "$TORCH_HOME" \
  "$XDG_CACHE_HOME" \
  "$TMPDIR" \
  "$MARLOWE_SCRATCH/logs"

cd "$ROOT"

module load slurm 2>/dev/null || true

ensure_apptainer() {
  if command -v apptainer >/dev/null 2>&1; then
    return
  fi
  module load apptainer 2>/dev/null || true
}

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "[ERROR] missing file: $1" >&2
    exit 2
  fi
}

require_dir() {
  if [[ ! -d "$1" ]]; then
    echo "[ERROR] missing directory: $1" >&2
    exit 2
  fi
}
