#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for light QC / restore.

#SBATCH --job-name=light_qc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
MANIFEST=$1
TEMPLATE_PDB=$2
OUTPUT_ROOT=$3
RESULT_DIR=$4

require_file "$MANIFEST"
require_file "$TEMPLATE_PDB"
mkdir -p "$OUTPUT_ROOT" "$RESULT_DIR"

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task light-qc \
  --manifest "$MANIFEST" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}" \
  --template-pdb "$TEMPLATE_PDB" \
  --output-root "$OUTPUT_ROOT" \
  --result-dir "$RESULT_DIR"
