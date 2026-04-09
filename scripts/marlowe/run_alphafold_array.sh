#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for AlphaFold2/3.
# Usage:
#   sbatch ... run_alphafold_array.sh MANIFEST.tsv af2|af3 OUTPUT_ROOT

#SBATCH --job-name=alphafold
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
MANIFEST=$1
BACKEND=$2
OUTPUT_ROOT=$3

require_file "$MANIFEST"
mkdir -p "$OUTPUT_ROOT"

if [[ -n "${AF2_SIF:-}" || -n "${AF3_SIF:-}" || -n "${COLABFOLD_SIF:-}" ]]; then
  ensure_apptainer
fi

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task alphafold \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --backend "$BACKEND" \
  --output-root "$OUTPUT_ROOT"
