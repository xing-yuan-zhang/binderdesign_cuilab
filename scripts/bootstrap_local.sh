#!/usr/bin/env bash

set -euo pipefail

PYTHON_EXE="${1:-python3}"
VENV_DIR="${VENV_DIR:-.venv}"

"$PYTHON_EXE" -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install --upgrade pip
"$VENV_DIR/bin/pip" install -e .

mkdir -p projects logs

echo
echo "Bootstrap complete."
echo "Activate with:"
echo "  source $VENV_DIR/bin/activate"
echo
echo "Then you can run:"
echo "  binderdesign-cuilab --help"
echo "  binderdesign-cuilab run-project --help"
echo "  binderdesign-cuilab table --help"
