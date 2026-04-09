#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for Rosetta DDG / interface dG negative filtering.

#SBATCH --job-name=ros_ddg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
MANIFEST=$1
WORK_ROOT=$2
RESULT_DIR=$3

: "${ROSETTA_BIN:?Set ROSETTA_BIN}"
: "${ROSETTA_DB:?Set ROSETTA_DB}"

require_file "$MANIFEST"
require_dir "$ROSETTA_BIN"
require_dir "$ROSETTA_DB"
mkdir -p "$WORK_ROOT" "$RESULT_DIR"

if [[ -n "${ROSETTA_SIF:-}" ]]; then
  ensure_apptainer
fi

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task rosetta-ddg \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --work-root "$WORK_ROOT" \
  --result-dir "$RESULT_DIR"
