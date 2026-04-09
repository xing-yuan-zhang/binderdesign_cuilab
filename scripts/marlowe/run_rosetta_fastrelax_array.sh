#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for Rosetta FastRelax robustness sampling.

#SBATCH --job-name=ros_relax
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
MANIFEST=$1
WORK_ROOT=$2
RAW_DIR=$3

: "${ROSETTA_BIN:?Set ROSETTA_BIN}"
: "${ROSETTA_DB:?Set ROSETTA_DB}"

require_file "$MANIFEST"
require_dir "$ROSETTA_BIN"
require_dir "$ROSETTA_DB"
mkdir -p "$WORK_ROOT" "$RAW_DIR"

if [[ -n "${ROSETTA_SIF:-}" ]]; then
  ensure_apptainer
fi

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task rosetta-fastrelax \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --work-root "$WORK_ROOT" \
  --raw-dir "$RAW_DIR"
