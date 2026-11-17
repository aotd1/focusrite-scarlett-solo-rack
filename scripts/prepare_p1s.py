"""Prepare a published variant with a portable P1S PETG setup.

Run Bambu Studio on build/p1s/<variant>/<variant>-project.3mf to slice and save it.
No printer connection or print command is used.
"""
from pathlib import Path
import argparse
import json
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from PIL import Image
import io

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('variant', nargs='?', default='knob-cup',
                    choices=['knob-cup', 'centered-solid', 'centered-front-vents'])
variant = parser.parse_args().variant
label = {'knob-cup':'Knob cup', 'centered-solid':'Centered solid',
         'centered-front-vents':'Centered front vents'}[variant]
OUT = ROOT / 'build/p1s' / variant
OUT.mkdir(parents=True, exist_ok=True)
mesh = trimesh.load_mesh(ROOT / 'models' / variant / (variant+'.stl'))
assert mesh.is_watertight and mesh.is_winding_consistent
assert len(mesh.split()) == 1 and abs(mesh.extents[2] - 84.2) < .002
center = mesh.bounds.mean(axis=0)
ns = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
ET.register_namespace('', ns)
model = ET.Element(f'{{{ns}}}model', {'unit': 'millimeter'})
for name, value in {
    'Application': 'BambuStudio-02.09.01.70',
    'BambuStudio:3mfVersion': '1',
    'Title': f'Scarlett Solo 4th Gen - {label} - P1S PETG',
    'Description': '84.2 mm depth; front face down; 45 degree bed rotation.',
    'License': 'CC BY-NC 4.0',
    'Copyright': 'Remix of Jonathan Lemoine, Printables model 1565705.',
    'Thumbnail_Middle': '/Metadata/plate_1.png',
    'Thumbnail_Small': '/Metadata/plate_1_small.png',
}.items():
    ET.SubElement(model, f'{{{ns}}}metadata', {'name': name}).text = value
resources = ET.SubElement(model, f'{{{ns}}}resources')
obj = ET.SubElement(resources, f'{{{ns}}}object', {'id': '1', 'type': 'model'})
geometry = ET.SubElement(obj, f'{{{ns}}}mesh')
vertices = ET.SubElement(geometry, f'{{{ns}}}vertices')
for v in mesh.vertices - center:
    ET.SubElement(vertices, f'{{{ns}}}vertex', dict(zip('xyz', map(lambda x: f'{x:.9f}', v))))
triangles = ET.SubElement(geometry, f'{{{ns}}}triangles')
for face in mesh.faces:
    ET.SubElement(triangles, f'{{{ns}}}triangle', dict(zip(('v1','v2','v3'), map(str, face))))
build = ET.SubElement(model, f'{{{ns}}}build')
ET.SubElement(build, f'{{{ns}}}item', {'objectid': '1', 'printable': '1',
    'transform': '1 0 0 0 1 0 0 0 1 ' + ' '.join(f'{x:.9f}' for x in center)})
config = ET.Element('config')
obj = ET.SubElement(config, 'object', {'id': '1'})
def meta(parent, key, value):
    ET.SubElement(parent, 'metadata', {'key':key, 'value':str(value)})
meta(obj, 'name', label+' - 84.2 mm')
meta(obj, 'extruder', 1)
ET.SubElement(obj, 'metadata', {'face_count':str(len(mesh.faces))})
part = ET.SubElement(obj, 'part', {'id':'1', 'subtype':'normal_part'})
meta(part, 'name', label+' - 84.2 mm')
meta(part, 'source_file', variant+'.stl')
meta(part, 'matrix', '1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1')
for axis, value in zip('xyz', center):
    meta(part, f'source_offset_{axis}', f'{value:.9f}')
plate = ET.SubElement(config, 'plate')
for key, value in {'plater_id':1, 'plater_name':label, 'locked':'true',
                   'filament_map_mode':'Auto For Flush', 'filament_maps':'1',
                   'filament_volume_maps':'0', 'gcode_file':''}.items():
    meta(plate, key, value)
instance = ET.SubElement(plate, 'model_instance')
for key, value in {'object_id':1, 'instance_id':0, 'identify_id':1}.items():
    meta(instance, key, value)
image = Image.open(ROOT / 'assets/gallery' / (variant+'.png')).convert('RGBA')
image.thumbnail((512,512))
preview = Image.new('RGBA', (512,512), '#f4f5f7')
preview.alpha_composite(image, ((512-image.width)//2,(512-image.height)//2))
png = io.BytesIO(); preview.save(png, format='PNG')
small = io.BytesIO(); preview.resize((128,128)).save(small, format='PNG')
files = {
    '[Content_Types].xml': b'''<?xml version="1.0"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/><Default Extension="png" ContentType="image/png"/><Default Extension="gcode" ContentType="text/x.gcode"/></Types>''',
    '_rels/.rels': b'''<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/><Relationship Target="/Metadata/plate_1.png" Id="rel-2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/thumbnail"/></Relationships>''',
    '3D/3dmodel.model': ET.tostring(model, encoding='utf-8', xml_declaration=True),
    'Metadata/model_settings.config': ET.tostring(config, encoding='utf-8', xml_declaration=True),
    'Metadata/project_settings.config': (ROOT/'profiles/p1s/settings.json').read_bytes(),
    'Metadata/plate_1.png': png.getvalue(),
    'Metadata/plate_1_small.png': small.getvalue(),
}
with zipfile.ZipFile(OUT/(variant+'-project.3mf'), 'w', zipfile.ZIP_DEFLATED) as archive:
    for name, data in files.items():
        archive.writestr(name, data)
print('Prepared', (OUT/(variant+'-project.3mf')).relative_to(ROOT))
