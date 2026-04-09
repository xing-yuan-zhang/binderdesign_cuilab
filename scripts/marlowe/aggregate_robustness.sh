#!/usr/bin/env bash
# BinderDesign-CUILab stage wrapper for robustness aggregation.

#SBATCH --job-name=agg_robust
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

set -euo pipefail

source "$(dirname "$0")/common.sh"

BINDERDESIGN_CUILAB_PYTHON="${BINDERDESIGN_CUILAB_PYTHON:-python}"
RAW_DIR=$1
OUTPUT_TSV=$2

require_dir "$RAW_DIR"
"$BINDERDESIGN_CUILAB_PYTHON" -m binderdesign_cuilab.cli task aggregate-robustness \
  --raw-dir "$RAW_DIR" \
  --output "$OUTPUT_TSV"
