#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for task geometry filtering.

#SBATCH --job-name=geometry
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
MANIFEST=$1
RESULT_DIR=$2
HOTSPOTS=$3
BINDER_LENGTH_MIN=$4
BINDER_LENGTH_MAX=$5

require_file "$MANIFEST"
mkdir -p "$RESULT_DIR"

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task geometry \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --result-dir "$RESULT_DIR" \
  --hotspots "$HOTSPOTS" \
  --binder-length-min "$BINDER_LENGTH_MIN" \
  --binder-length-max "$BINDER_LENGTH_MAX"
