#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for ProteinMPNN.

#SBATCH --job-name=proteinmpnn
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

: "${PROTEINMPNN_PYTHON:?Set PROTEINMPNN_PYTHON}"
: "${PROTEINMPNN_SCRIPT:?Set PROTEINMPNN_SCRIPT}"
BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"

MANIFEST=$1
OUTPUT_ROOT=$2

require_file "$MANIFEST"
require_file "$PROTEINMPNN_SCRIPT"
mkdir -p "$OUTPUT_ROOT"

if [[ -n "${PROTEINMPNN_SIF:-}" ]]; then
  ensure_apptainer
fi

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task proteinmpnn \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --output-root "$OUTPUT_ROOT" \
  --pdb-column "${MPNN_PDB_COLUMN:-restored_pdb}"
