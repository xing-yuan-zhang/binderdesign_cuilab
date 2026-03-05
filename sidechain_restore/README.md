# Description

Batch restore TG2 backbone for RFdiffusion designs in specific seed folders.

## Run (same interface)

```bash
python sidechain_restore.py template.pdb output_root seed_0 [seed_1 ...]
```

## Library

```python
from collision.restore import restore_one, main
from collision import config
```

## Customize contigs / chains

Edit `tg2_restore/config.py`.
