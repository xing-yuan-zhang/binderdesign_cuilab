#!/usr/bin/env bash
#SBATCH --job-name=rfdiffusion
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)

SIF="$ROOT/images/rfd.sif"
IN="$ROOT/inputs/FBP17_2EFL.pdb"
CKPT="$ROOT/models/Complex_base_ckpt.pt"

CONFIG=$1
NAME=$(basename $CONFIG .txt)

SEED="${SLURM_ARRAY_TASK_ID}"

OUTDIR="$ROOT/outputs/${NAME}/seed_${SEED}"
OUTPREFIX="$OUTDIR/${NAME}"

mkdir -p "$ROOT/logs" "$ROOT/cache" "$ROOT/envs" "$ROOT/schedules" "$OUTDIR"

export DGLBACKEND=pytorch
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export HF_HOME="$ROOT/cache/hf"
export TORCH_HOME="$ROOT/cache/torch"
export XDG_CACHE_HOME="$ROOT/cache/xdg"

echo "[INFO] job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID"
echo "[INFO] config=$CONFIG"

module load apptainer 2>/dev/null || true
nvidia-smi || true

for f in "$SIF" "$IN" "$CKPT" "$CONFIG"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] missing file: $f" >&2
    exit 2
  fi
done

CONTIG=$(grep CONTIG $CONFIG | cut -d= -f2)
HOTSPOTS=$(grep HOTSPOTS $CONFIG | cut -d= -f2)

apptainer exec --nv \
  --bind "$ROOT:/work" \
  --bind "$ROOT/schedules:/opt/RFdiffusion/schedules" \
  "$SIF" \
  bash -lc "
set -euo pipefail

export PYTHONNOUSERSITE=1
export DGLBACKEND=pytorch
export PYTHONPATH=/opt/RFdiffusion:\${PYTHONPATH:-}

VENV_PY=/work/envs/rfdpy/bin/python

echo '[INFO] python:' \$($VENV_PY -V)

$VENV_PY -c \"import torch;print('torch',torch.__version__, 'cuda',torch.cuda.is_available())\"
$VENV_PY -c \"import rfdiffusion;print('rfdiffusion ok')\"

echo '[INFO] starting inference'

$VENV_PY /opt/RFdiffusion/scripts/run_inference.py \
  --config-name base \
  inference.model_directory_path=/work/models \
  inference.input_pdb=/work/inputs/FBP17_2EFL.pdb \
  +inference.ckpt_path=/work/models/Complex_base_ckpt.pt \
  inference.output_prefix=/work/outputs/${NAME}/seed_${SEED}/${NAME} \
  inference.num_designs=200 \
  contigmap.contigs=$CONTIG \
  +ppi.hotspot_res=$HOTSPOTS \
  +seed=$SEED

echo '[INFO] finished inference'
"