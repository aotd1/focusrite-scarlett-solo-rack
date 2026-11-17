# Build

The model files use millimetres. Each models/ subfolder contains:

- <variant>.3mf — print-oriented geometry, front face on bed, rotated 45°.
- <variant>.stl — the same print orientation.
- service.stl — the rack-oriented geometry, with the front face at Y=0.

Use Python 3.14 and a virtual environment:

```sh
python -m pip install -r source/requirements.txt
python scripts/build.py
python scripts/check.py
```

The build includes section, mounting-profile, slot, cup and export checks.
It writes temporary measurement reports and an optional fit coupon under build/.
Original-source attribution is in ATTRIBUTION.md.

To prepare a release ZIP from the tracked tree:

```sh
python scripts/package.py
```

This does not print or publish to MakerWorld.
