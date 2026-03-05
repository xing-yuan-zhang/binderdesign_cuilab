# Description

Triage binder designs from restored complexes (v5, ≥3 helices).

## Run

```bash
python triage_restored_helix.py /path/to/restored_folder
```

## Library

```python
from helix.triage import triage_one, main
from helix import config

ok, reason, metrics = triage_one("design.pdb")
score = 100.0 * metrics["hotspots_contacted"] + metrics["iface_contacts"]
```

## Customize thresholds

Edit `helix/config.py`.