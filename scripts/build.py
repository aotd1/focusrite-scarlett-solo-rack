"""Rebuild and check the three rack variants. No printer or network access."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parents[1]
for script in [
    'source/knob-cup/build.py',
    'source/knob-cup/verify_panel_profile.py',
    'source/knob-cup/verify_knob_cup.py',
    'source/knob-cup/verify_geometry.py',
    'source/centered/build.py',
    'source/centered/verify.py',
]:
    result=subprocess.run([sys.executable,script],cwd=root,capture_output=True,text=True)
    if result.returncode:raise SystemExit(script+'\n'+result.stderr)
    print('PASS:',script)
