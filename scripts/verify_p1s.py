"""Check the exported P1S project and its actual Bambu Studio toolpaths."""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import json
import math
import re
import zipfile
import numpy as np
import trimesh
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('variant', nargs='?', default='knob-cup',
                    choices=['knob-cup','centered-solid','centered-front-vents'])
variant = parser.parse_args().variant
project = ROOT/'profiles/p1s'/(variant+'-p1s-petg.3mf')
gcode_path = ROOT/'build/p1s'/variant/'verified/plate_1.gcode'
gcode = gcode_path.read_text()
source = trimesh.load_mesh(ROOT/'models'/variant/(variant+'.stl'))
scene = trimesh.load(project, force='scene')
mesh = scene.to_mesh()
assert mesh.is_watertight and mesh.is_winding_consistent and len(mesh.split()) == 1
error = max(cKDTree(source.vertices).query(mesh.vertices)[0].max(),
            cKDTree(mesh.vertices).query(source.vertices)[0].max())
assert error < .002, error
assert np.allclose(mesh.extents, source.extents, atol=.002)
assert abs(mesh.extents[2]-84.2) < .002
with zipfile.ZipFile(project) as archive:
    settings = json.loads(archive.read('Metadata/project_settings.config'))
    assert settings['filament_ids'] == ['GFG99']
    assert settings['filament_settings_id'] == ['Generic PETG']
    assert settings['printer_settings_id'] == 'Bambu Lab P1S 0.4 nozzle'
    process_changes = set(settings['different_settings_to_system'][0].split(';'))
    filament_changes = set(settings['different_settings_to_system'][1].split(';'))
    assert {'enable_support','support_type','support_on_build_plate_only',
            'sparse_infill_density','sparse_infill_pattern'} <= process_changes
    assert {'nozzle_temperature','textured_plate_temp'} <= filament_changes
    for name in archive.namelist():
        if name.endswith(('.model','.config','.json')):
            data = archive.read(name)
            assert not any(s in data for s in (b'/Users/',b'file://',b'access_token')), name
cfg = dict(re.findall(r'^; ([a-zA-Z_0-9]+) = (.*)$', gcode, re.M))
expected = {'printer_model':'Bambu Lab P1S','nozzle_diameter':'0.4',
    'layer_height':'0.2','initial_layer_print_height':'0.2','wall_loops':'2',
    'top_shell_layers':'5','bottom_shell_layers':'5','sparse_infill_density':'20%',
    'sparse_infill_pattern':'gyroid','enable_support':'1','support_type':'normal(auto)',
    'support_on_build_plate_only':'1','support_top_z_distance':'0.2',
    'support_interface_top_layers':'2','enable_prime_tower':'0','brim_type':'no_brim',
    'filament_type':'PETG','nozzle_temperature':'250','nozzle_temperature_initial_layer':'250',
    'textured_plate_temp':'75','textured_plate_temp_initial_layer':'75',
    'curr_bed_type':'Textured PEI Plate','filament_max_volumetric_speed':'12'}
for key, value in expected.items():
    assert cfg.get(key) == value, (key,cfg.get(key),value)
    saved = settings.get(key)
    if saved is not None:
        if isinstance(saved,list): saved = saved[0]
        assert str(saved) == value, (key,saved,value)
heights = np.array([float(x) for x in re.findall(r'^; Z_HEIGHT: ([\d.]+)',gcode,re.M)])
assert len(heights) == 421 and np.allclose(heights,np.arange(1,422)*.2)
assert re.search(r'^M104 S0\b',gcode,re.M) and re.search(r'^M140 S0\b',gcode,re.M)
# Transform the CAD cup datum to the fixed, print-oriented coordinate system.
service = trimesh.load_mesh(ROOT/'models/knob-cup/service.stl')
rotation = (trimesh.transformations.rotation_matrix(np.pi/4,[0,0,1]) @
            trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
turned = service.copy(); turned.apply_transform(rotation)
shift = source.bounds[0]-turned.bounds[0]
center = (rotation @ [-80.725,0,20.45,1])[:3]+shift
center[:2] -= np.array([float(v) for v in settings['extruder_offset'][0].split('x')])
xy = np.zeros(2); height=0.; feature='Custom'; relative=True; epos=0.; feed=0.
nozzle=bed=None; temperatures=set(); paths=[]; flows=[]; roles=Counter()
membrane=Counter(); pocket=Counter(); support_interface=set()
def distance(a,b):
    d=b-a; den=np.dot(d,d)
    t=np.clip(np.dot(center[:2]-a,d)/den,0,1) if den else 0
    return np.linalg.norm(a+t*d-center[:2])
for raw in gcode.splitlines():
    if raw.startswith('; Z_HEIGHT:'): height=round(float(raw.split(':')[1]),3)
    if raw.startswith('; FEATURE:'): feature=raw.split(': ',1)[1]
    words=raw.split(';',1)[0].split()
    if not words: continue
    cmd=words[0]
    v={w[0]:float(w[1:]) for w in words[1:] if re.fullmatch(r'[XYZEFIJSP][-+]?\d*\.?\d+',w)}
    if cmd in ('M104','M109') and 'S' in v: nozzle=v['S']
    if cmd in ('M140','M190') and 'S' in v: bed=v['S']
    if cmd=='M83':relative=True
    if cmd=='M82':relative=False
    if cmd=='G92' and 'E' in v:epos=v['E']
    if cmd not in ('G0','G1','G2','G3'):continue
    end=np.array([v.get('X',xy[0]),v.get('Y',xy[1])])
    extrusion=v.get('E',0) if relative else v.get('E',epos)-epos
    if 'E' in v:epos=v['E']
    feed=v.get('F',feed); points=np.array([xy,end]); length=np.linalg.norm(end-xy)
    if cmd in ('G2','G3') and ('I' in v or 'J' in v):
        arc=xy+[v.get('I',0),v.get('J',0)];radius=np.linalg.norm(xy-arc)
        start=math.atan2(*(xy-arc)[::-1]); stop=math.atan2(*(end-arc)[::-1])
        angle=(stop-start)%(2*np.pi) if cmd=='G3' else -((start-stop)%(2*np.pi))
        if abs(angle)<1e-10:angle=2*np.pi if cmd=='G3' else -2*np.pi
        length=abs(angle)*radius
        theta=np.linspace(start,start+angle,max(2,int(np.ceil(length/.2))+1))
        points=arc+radius*np.column_stack((np.cos(theta),np.sin(theta)))
        points[0]=xy;points[-1]=end
    if extrusion>0 and length>.01 and feature!='Custom' and height>0:
        paths.extend(points);temperatures.add((nozzle,bed));roles[feature]+=1
        if length>=.1:flows.append(extrusion*np.pi*(1.75/2)**2/length*feed/60)
        near=min(distance(a,b) for a,b in zip(points[:-1],points[1:]))
        support='Support' in feature
        if feature=='Support interface' and 12.4<=height<=13.05 and near<19:
            support_interface.add(height)
        if not support and 13.05<=height<=13.45 and near<2:membrane[height]+=1
        if not support and 13.45<=height<=15.4 and near<2:pocket[height]+=1
    xy=end
paths=np.array(paths)
assert paths.min()>=18 and paths.max()<=238, (paths.min(0),paths.max(0))
assert temperatures=={(250.,75.)},temperatures
assert max(flows)<12.1,max(flows)
assert roles['Sparse infill'] and roles['Support'] and roles['Support interface']
if variant == 'knob-cup':
    assert support_interface
    assert membrane[13.2] and membrane[13.4],membrane
    assert not any(pocket.values()),pocket
estimate = re.search(r'total estimated time: ([^\n;]+)',gcode).group(1)
weight = float(re.search(r'total filament weight \[g\] : ([\d.]+)',gcode).group(1))
report = {'printer':'Bambu Lab P1S','variant':variant,'verified_settings':expected,
    'mesh_vertex_max_deviation_mm':float(error),'size_mm':mesh.extents.tolist(),
    'watertight':True,'components':1,'layers':len(heights),
    'last_layer_mm':float(heights[-1]),
    'toolpath_xy_bounds_mm':[paths.min(0).tolist(),paths.max(0).tolist()],
    'max_flow_mm3_s_on_moves_ge_0p1mm':max(flows), 'estimated_time':estimate.strip(),
    'estimated_filament_g':weight,'sha256':hashlib.sha256(project.read_bytes()).hexdigest()}
if variant == 'knob-cup':
    report.update(cup_support_interface_layers_mm=sorted(support_interface),
                  cup_membrane_extrusion_layers_mm=sorted(membrane), magnet_pocket_center_clear=True)
report_name = 'validation.json' if variant == 'knob-cup' else variant+'-validation.json'
(ROOT/'profiles/p1s'/report_name).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
