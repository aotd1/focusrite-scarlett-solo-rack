# CC BY-NC 4.0. Original model: Jonathan Lemoine.
# Modified mesh source/verification/rendering; see ATTRIBUTION.md and LICENSE.md.
"""Reproducible vented Scarlett remix with mirrored layout and blind knob cup.
Preserves rack holes, upper edge and Scarlett inner width; sleeve depth is configurable.
All inputs and geometry units are millimetres. No printer calls.
"""
from pathlib import Path
import json,hashlib,math,io,zipfile
from xml.etree import ElementTree as ET
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUTPUT=PROJECT/'models/knob-cup'
STATE=PROJECT/'build/knob-cup'
OUTPUT.mkdir(parents=True,exist_ok=True)
STATE.mkdir(parents=True,exist_ok=True)
cfg=json.loads((ROOT/'parameters.json').read_text())
source=PROJECT/'source/original.stl'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='2d9ca6b355928fdd7d8014ae271b9dde45f5652f91894ac32cf06bf1c6dd52a4'
original=trimesh.load_mesh(source,process=True)
original.apply_translation([0,102.5999984741211,-181.125])
assert original.is_watertight and original.is_winding_consistent
# Measured planes from source mesh; front face Y=0, original rack bottom Z=0.
xl,xli,xri,xr=-94.125,-91.125,52.875,55.875
floor,roof=1.3000030517578125,42.69999694824219
target=cfg['measured_body_height_mm']+cfg['vertical_clearance_total_mm']
delta=target-(roof-floor)
assert 0<delta<8 and cfg['remove_rubber_feet']
source_depth=float(original.bounds[1,1]);depth=source_depth;split=cfg['stretch_split_height_mm']
def box(lo,hi):
 lo=np.array(lo,dtype=float);hi=np.array(hi,dtype=float)
 b=trimesh.creation.box(hi-lo);b.apply_translation((lo+hi)/2);return b
def boolean(operation,meshes):
 result=getattr(trimesh.boolean,operation)(meshes,engine='manifold')
 assert result.is_watertight and result.is_winding_consistent
 return result
roi=box([xl,-1,-20],[xr,depth+1,split])
lower=boolean('intersection',[original,roi]);upper=boolean('difference',[original,roi]);lower.apply_translation([0,0,-delta])
bridge=[box([xl,0,split-delta],[xli,depth,split]),box([xri,0,split-delta],[xr,depth,split])]
remix=boolean('union',[upper,lower,*bridge])
# Thin only the rear sleeve; leave 24 mm near rack ears and fan-control carrier unchanged.
t=cfg['rear_side_wall_mm'];front=cfg['retain_front_reinforced_depth_mm']
cuts=[box([xl-.01,front,-20],[xli-t,depth+1,60]),box([xri+t,front,-20],[xr+.01,depth+1,60])]
remix=boolean('difference',[remix,*cuts])
# Rebuild only the 2.6mm faceplate from its section: extend its rounded
# outer silhouette downward everywhere, and reuse every opening verbatim.
# This avoids filling any of the newly enlarged Scarlett aperture.
if cfg.get('full_width_lower_panel',False) or cfg.get('lower_extension_width_mm'):
 import manifold3d
 section_y=1.0;d=cfg['front_panel_depth_mm']
 source_loops=original.section(plane_normal=[0,1,0],plane_origin=[0,section_y,20]).discrete
 outer=max(source_loops,key=lambda v:np.ptp(v[:,0]))[:,[0,2]].copy()
 extension_width=cfg.get('lower_extension_width_mm')
 if not extension_width:outer[outer[:,1]<22,1]-=delta
 existing_loops=remix.section(plane_normal=[0,1,0],plane_origin=[0,section_y,20]).discrete
 holes=[v[:,[0,2]] for v in existing_loops if np.ptp(v[:,0])<200]
 if extension_width:
  # Restore the original ear silhouette and add only the centered lower lip.
  half=extension_width/2
  extension=np.array([[-half,-delta],[half,-delta],[half,4],[-half,4]])
  cross=(manifold3d.CrossSection([outer],manifold3d.FillRule.EvenOdd)+manifold3d.CrossSection([extension],manifold3d.FillRule.EvenOdd))-manifold3d.CrossSection(holes,manifold3d.FillRule.EvenOdd)
 else:cross=manifold3d.CrossSection([outer,*holes],manifold3d.FillRule.EvenOdd)
 raw=cross.extrude(d).to_mesh()
 panel=trimesh.Trimesh(vertices=raw.vert_properties[:,:3],faces=raw.tri_verts,process=True)
 panel.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]));panel.apply_translation([0,d,0])
 rear=boolean('intersection',[remix,box([-140,d,-20],[140,depth+1,60])])
 remix=boolean('union',[panel,rear])
# Shorten the sleeve from the rear, preserving the full front profile and cup.
target_depth=cfg.get('focusrite_depth_mm',source_depth)
# Decimal 84.2 and the source STL plane differ by ~1.1e-6 mm.
# Snap only within numerical tolerance; do not translate or scale the sleeve.
assert cfg['retain_front_reinforced_depth_mm']<target_depth<=source_depth+1e-5
if abs(target_depth-source_depth)<=1e-5:
 target_depth=source_depth
if target_depth<source_depth:
 # Retain the source's shaped rear edge, moving its solid collar forward.
 # A simple cut would lose the small end features of the original sleeve.
 collar=target_depth-max(y[1] for y in cfg['ventilation']['depth_ranges_y_mm'])
 assert collar>=5
 rear_collar=boolean('intersection',[remix,box([-140,source_depth-collar,-20],[140,source_depth+.01,60])])
 rear_collar.apply_translation([0,target_depth-source_depth,0])
 shortened=boolean('intersection',[remix,box([-140,-1,-20],[140,target_depth-collar+.01,60])])
 remix=boolean('union',[shortened,rear_collar])
depth=target_depth
cup=cfg.get('knob_cup',{})
cup_report=None
if cup.get('enabled',False):
 # Remove the old rear fan-controller carrier; keep the full-width fascia.
 d=cfg['front_panel_depth_mm'];cx=cup['center_x_before_mirror_mm'];cz=(44-delta)/2
 remix=boolean('difference',[remix,box([xr,d,-10],[110,26,60])])
 inside_radius=(cup['device_diameter_mm']+cup['diameter_clearance_mm'])/2
 outer_radius=inside_radius+cup['wall_mm'];inside_depth=cup['inside_depth_mm'];outside_depth=inside_depth+cup['bottom_mm']
 def cylinder_y(radius,y0,y1):
  c=trimesh.creation.cylinder(radius=radius,height=y1-y0,sections=256);c.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]));c.apply_translation([cx,(y0+y1)/2,cz]);return c
 remix=boolean('union',[remix,cylinder_y(outer_radius,0,outside_depth)])
 remix=boolean('difference',[remix,cylinder_y(inside_radius,-1,inside_depth)])
 cup_report={'inside_diameter_mm':inside_radius*2,'inside_depth_mm':inside_depth,'wall_mm':cup['wall_mm'],'bottom_mm':cup['bottom_mm'],'outside_diameter_mm':outer_radius*2,'outside_depth_mm':outside_depth,'center_x_before_mirror_mm':cx,'center_z_mm':cz,'closed_bottom':True,'mouth_plane_y_mm':0,'inside_floor_plane_y_mm':inside_depth,'outer_back_plane_y_mm':outside_depth}
 magnet=cup.get('magnet_pocket')
 if magnet:
  pocket_floor=outside_depth-magnet['depth_mm']
  assert pocket_floor>inside_depth and magnet['diameter_mm']<inside_radius*2
  remix=boolean('difference',[remix,cylinder_y(magnet['diameter_mm']/2,pocket_floor,outside_depth+.1)])
  cup_report['magnet_pocket']={**magnet,'floor_plane_y_mm':pocket_floor,'rear_plane_y_mm':outside_depth,'remaining_membrane_mm':pocket_floor-inside_depth}
vent=cfg.get('ventilation',{})
vent_report=None
if vent.get('enabled',False):
 # Slots run along Y, the vertical print direction; rounded ends span 6mm.
 before_vent_volume=float(remix.volume);before_vent_euler=int(remix.euler_number)
 radius=vent['slot_width_mm']/2;tools=[];locations=[]
 def capsule(center,y0,y1,axis):
  if axis=='z':
   lo,hi=-8,50
   pieces=[box([center-radius,y0+radius,lo],[center+radius,y1-radius,hi])]
   for y in [y0+radius,y1-radius]:
    c=trimesh.creation.cylinder(radius=radius,height=hi-lo,sections=32);c.apply_translation([center,y,(lo+hi)/2]);pieces.append(c)
  else:
   lo,hi=xl-5,xr+5
   pieces=[box([lo,y0+radius,center-radius],[hi,y1-radius,center+radius])]
   for y in [y0+radius,y1-radius]:
    c=trimesh.creation.cylinder(radius=radius,height=hi-lo,sections=32);c.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[0,1,0]));c.apply_translation([(lo+hi)/2,y,center]);pieces.append(c)
  return boolean('union',pieces)
 xs=np.linspace(xli+vent['top_bottom_center_inset_mm'],xri-vent['top_bottom_center_inset_mm'],vent['top_bottom_columns'])
 zs=np.linspace(floor-delta+vent['side_center_inset_mm'],roof-vent['side_center_inset_mm'],vent['side_columns'])
 for y0,y1 in vent['depth_ranges_y_mm']:
  assert np.isclose(y1-y0,vent['slot_length_mm']) and y0>front and y1<depth-5
  for x in xs:
   tools.append(capsule(float(x),y0,y1,'z'))
   for face,z in [('top',roof+.45),('bottom',floor-delta-.45)]:locations.append({'face':face,'center_mm':[float(x),(y0+y1)/2,float(z)]})
  for z in zs:
   tools.append(capsule(float(z),y0,y1,'x'))
   for face,x in [('left',xli-t/2),('right',xri+t/2)]:locations.append({'face':face,'center_mm':[float(x),(y0+y1)/2,float(z)]})
 remix=boolean('difference',[remix,*tools])
 assert remix.euler_number==before_vent_euler-2*len(locations)
 vent_report={'slot_width_mm':vent['slot_width_mm'],'slot_length_mm':vent['slot_length_mm'],'hole_count':len(locations),'holes':locations,'before_volume_mm3':before_vent_volume,'removed_volume_mm3':before_vent_volume-float(remix.volume),'before_euler':before_vent_euler,'after_euler':int(remix.euler_number),'front_solid_band_mm':min(y[0] for y in vent['depth_ranges_y_mm']),'rear_solid_band_mm':depth-max(y[1] for y in vent['depth_ranges_y_mm'])}
mirrored=cfg.get('mirror_x_for_focusrite_on_right',False)
if mirrored:
 remix.apply_scale([-1,1,1])
 if vent_report:
  for hole in vent_report['holes']:
   hole['center_mm'][0]*=-1
   hole['face']={'left':'right','right':'left'}.get(hole['face'],hole['face'])
if cup_report:cup_report['center_x_mm']=-cup_report['center_x_before_mirror_mm'] if mirrored else cup_report['center_x_before_mirror_mm']
assert len(remix.split())==1
remix.export(OUTPUT/'service.stl')
# Place the front panel on the print bed; 45deg in XY fits a 256mm plate.
print_mesh=remix.copy();print_mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]));print_mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/4,[0,0,1]));print_mesh.apply_translation(-print_mesh.bounds[0]);print_mesh.apply_translation([(256-print_mesh.extents[0])/2,(256-print_mesh.extents[1])/2,0]);assert (print_mesh.extents[:2]<250).all()
# STL float32 coordinates collapse a few submicron boolean slivers after rotation.
# Remove those degenerate/duplicate facets before exporting either print format.
print_mesh.vertices=print_mesh.vertices.astype(np.float32).astype(np.float64)
print_mesh.merge_vertices(digits_vertex=5)
print_mesh.update_faces(print_mesh.nondegenerate_faces())
print_mesh.update_faces(print_mesh.unique_faces())
print_mesh.remove_unreferenced_vertices()
assert print_mesh.is_watertight and print_mesh.is_winding_consistent
print_mesh.export(OUTPUT/'knob-cup.stl')
assert trimesh.load_mesh(OUTPUT/'knob-cup.stl').is_watertight
scene=trimesh.Scene();scene.add_geometry(print_mesh,geom_name='Scarlett Solo 4Gen D84.2 K40.8 mm - vented remix',node_name='Scarlett Solo + Zigbee knob')
raw_3mf=trimesh.exchange.threemf.export_3MF(scene)
ns='http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
metadata={
 'Title':'Scarlett Solo 4th Gen + Zigbee knob - 10-inch rack remix, D84.2 K40.8',
 'Designer':'Remix of Jonathan Lemoine (@JonathanLemoi_601472)',
 'Description':'Modified Gen3 mount: feet-off Gen4 opening, mirrored, vented, blind knob cup and rear magnet recess; physical fit untested.',
 'Copyright':'Original model: Jonathan Lemoine; modified version: remix contributors',
 'LicenseTerms':'CC BY-NC 4.0 - https://creativecommons.org/licenses/by-nc/4.0/ ; source license locally recorded, live recheck pending',
 'Source':'https://www.printables.com/model/1565705-focusrite-scarlettsolo-3rd-gen-10-inch-rack-mount',
}
buffer=io.BytesIO()
with zipfile.ZipFile(io.BytesIO(raw_3mf)) as zin,zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_DEFLATED) as zout:
 for name in zin.namelist():
  data=zin.read(name)
  if name=='3D/3dmodel.model':
   xml=ET.fromstring(data)
   for i,(key,value) in enumerate(metadata.items()):
    tag=ET.Element('{'+ns+'}metadata',{'name':key});tag.text=value;xml.insert(i,tag)
   data=ET.tostring(xml,encoding='utf-8',xml_declaration=True)
  zout.writestr(name,data)
(OUTPUT/'knob-cup.3mf').write_bytes(buffer.getvalue())
# Short ring for a cheap fit check, with exactly the final opening profile.
coupon_xlo,coupon_xhi=(-xr-.001,-xl+.001) if mirrored else (xl-.001,xr+.001)
coupon=boolean('intersection',[remix,box([coupon_xlo,-.01,-20],[coupon_xhi,cfg['fit_coupon_depth_mm'],60])]);coupon.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]));coupon.apply_translation(-coupon.bounds[0]);coupon.export(STATE/'fit-coupon.stl')
# Validate opening at front, reinforced sleeve and thinned rear sleeve.
sections=[]
for y in [1,12,27,depth-2]:
 sec=remix.section(plane_normal=[0,1,0],plane_origin=[0,y,20]);loops=sec.discrete;inner=[v for v in loops if np.isclose(np.ptp(v[:,0]),144,atol=.002) and np.isclose(np.ptp(v[:,2]),target,atol=.002)];assert len(inner)==1,(y,[(np.ptp(v[:,0]),np.ptp(v[:,2])) for v in loops]);v=inner[0];sections.append({'y_mm':y,'opening_width_mm':float(np.ptp(v[:,0])),'opening_height_mm':float(np.ptp(v[:,2])),'floor_z_mm':float(v[:,2].min()),'roof_z_mm':float(v[:,2].max())})
# Verify no source ear material was removed after layout mirroring.
# Existing faceplate rebuild adds material at some bevels; section tests verify
# hole/outline profiles, and review compares entire ears to the latest model.
checks=[]
protected=[('left_rack_ear_above_corner',box([-130,-1,4],[-110,8,60])),('right_rack_ear_above_corner',box([110,-1,4],[130,8,60]))]
if not cup.get('enabled',False):protected.append(('fan_controller_above_original_bottom',box([57.875,-1,0],[110,24,60])))
for name,region in protected:
 ear_reference=original.copy()
 if mirrored:ear_reference.apply_scale([-1,1,1])
 a=boolean('intersection',[ear_reference,region]);b=boolean('intersection',[remix,region]);removed=float(trimesh.boolean.difference([a,b],engine='manifold').volume);added=float(trimesh.boolean.difference([b,a],engine='manifold').volume);assert abs(removed)<.01,(name,removed);checks.append({'region':name,'original_material_removed_mm3':removed,'material_added_mm3':added})
report={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_author':'Jonathan Lemoine','license':'CC BY-NC 4.0','parameters':cfg,'target_opening_mm':[144,target],'source_opening_mm':[144,roof-floor],'opening_enlarged_down_mm':delta,'overall_size_mm':remix.extents.tolist(),'service_bounds_mm':remix.bounds.tolist(),'upper_edge_unchanged':bool(np.isclose(remix.bounds[1,2],44)),'source_ear_material_retained':checks,'sections':sections,'watertight':bool(remix.is_watertight),'winding_consistent':bool(remix.is_winding_consistent),'connected_components':len(remix.split()),'volume_mm3':float(remix.volume),'original_volume_mm3':float(original.volume),'minimum_top_bottom_wall_reference_mm':0.9,'rear_side_wall_mm':t,'front_reinforced_length_mm':front,'print_bounds_mm':print_mesh.bounds.tolist(),'print_size_mm':print_mesh.extents.tolist(),'coupon_size_mm':coupon.extents.tolist(),'fit_not_physically_tested':True,'printing_started':False,'modified_stl_sha256':hashlib.sha256((OUTPUT/'knob-cup.stl').read_bytes()).hexdigest()}
report['ventilation']=vent_report
report['knob_cup']=cup_report
report['mirrored_x']=mirrored
report['lower_extension_width_mm']=cfg.get('lower_extension_width_mm',float(remix.extents[0]))
report['focusrite_depth_mm']=depth
report['source_depth_mm']=source_depth
report['source_metadata_status']='locally recorded; live source page recheck required before upload'
report['source_model_url']='https://www.printables.com/model/1565705-focusrite-scarlettsolo-3rd-gen-10-inch-rack-mount'
report['depth_datum']='front face Y=0 to rear sleeve edge'
report['front_layout']={'focusrite':'right','knob_cup':'left'} if mirrored else {'focusrite':'left','knob_cup':'right'}
(STATE/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps(report,ensure_ascii=False,indent=2))
