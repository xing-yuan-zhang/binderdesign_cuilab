#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for RFdiffusion.

#SBATCH --job-name=rfdiffusion
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

: "${RFDIFFUSION_PYTHON:?Set RFDIFFUSION_PYTHON}"
: "${RFDIFFUSION_SCRIPT:?Set RFDIFFUSION_SCRIPT}"
BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"

TARGET_PDB=$1
CONFIG_FILE=$2
OUTPUT_ROOT=$3
NUM_DESIGNS=$4

require_file "$TARGET_PDB"
require_file "$CONFIG_FILE"
require_file "$RFDIFFUSION_SCRIPT"
mkdir -p "$OUTPUT_ROOT"

if [[ -n "${RFDIFFUSION_SIF:-}" ]]; then
  ensure_apptainer
fi

"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task rfdiffusion \
  --target-pdb "$TARGET_PDB" \
  --config-file "$CONFIG_FILE" \
  --output-root "$OUTPUT_ROOT" \
  --num-designs "$NUM_DESIGNS" \
  --index "${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"
