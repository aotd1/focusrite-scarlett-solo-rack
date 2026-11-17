"""Package the committed project and prepare named 3MF/media release assets."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

root=Path(__file__).resolve().parents[1]
config=json.loads((root/'publication.json').read_text())
version=config['version']
output=root/'build/release';output.mkdir(parents=True,exist_ok=True)
archive=output/f"focusrite-scarlett-solo-rack-v{version}.zip"
subprocess.run(['git','archive','--format=zip',f'--prefix=focusrite-scarlett-solo-rack-v{version}/',
                '-o',str(archive),'HEAD'],cwd=root,check=True)
files=[archive]
for variant in config['variants']:
    target=output/(variant+'.3mf')
    shutil.copyfile(root/'models'/variant/(variant+'.3mf'),target);files.append(target)
for entry in config['print_profiles']:
    source=root/entry['file']; profile=output/source.name
    shutil.copyfile(source,profile); files.append(profile)
for name in ['knob-demo.gif','knob-demo.mp4','installed-photo.png']:
    target=output/name;shutil.copyfile(root/'assets/knob'/name,target);files.append(target)
(output/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
print('Release files:',output.relative_to(root))
