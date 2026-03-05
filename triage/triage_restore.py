"""
Triage binder designs from restored complexes (v5, ≥3 helices).

- Input: folder containing restored PDBs
  (chain A = receptor, chain B = binder backbone).

- Output:
  - Prints per-design metrics + score to stdout.
  - Copies designs that pass all filters into <folder>/keepers.
  - Writes a tab-separated table with all designs + scores to:
        <folder>/triage_scores.tsv

Usage:
    python triage_restore.py /path/to/restored_folder
"""

import sys
from helix.triage import main

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)
    main(sys.argv[1])
