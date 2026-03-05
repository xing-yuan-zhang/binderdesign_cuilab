"""
Batch restore backbone for RFdiffusion designs in specific seed folders.

Usage:
    python sidechain_restore.py \
        TEMPLATE.pdb \
        OUTPUT_ROOT \
        seed_0 [seed_1 ...]

Args:
- template_pdb:
    template PDB (receptor in chain A, same contigs as RFdiffusion)

- output_root:
    Folder where restored PDBs will be written.
    All designs from all seeds go into this one folder.

- seed_dirs (1+):
    One or more folders like seed_0, seed_1, seed_2,
    each containing RFdiffusion PDBs:
      - chain A: cropped receptor
      - chain B: designed binder backbone

For each design PDB in each seed folder:
    output file name = <original_base>_restored_<seed_basename>.pdb
"""

import sys, os
from collision.restore import main

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(1)

    template_pdb = os.path.abspath(sys.argv[1])
    output_root  = os.path.abspath(sys.argv[2])
    seed_dirs    = sys.argv[3:]
    main(template_pdb, output_root, seed_dirs)
