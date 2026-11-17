# CC BY-NC 4.0. Original model: Jonathan Lemoine.
# Modified mesh source/verification/rendering; see ATTRIBUTION.md and LICENSE.md.
from pathlib import Path
import json,numpy as np,trimesh
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUTPUT=PROJECT/'models/knob-cup'
STATE=PROJECT/'build/knob-cup'
OUTPUT.mkdir(parents=True,exist_ok=True)
STATE.mkdir(parents=True,exist_ok=True)
m=trimesh.load_mesh(OUTPUT/'service.stl');r=json.loads((STATE/'validation.json').read_text());c=r['knob_cup'];cx=c['center_x_mm'];cz=c['center_z_mm'];p=c['magnet_pocket'];inner_d=c['inside_diameter_mm'];outer_d=c['outside_diameter_mm'];assert np.isclose(inner_d,40.8) and np.isclose(outer_d,44.8)
assert cx<0 and r['front_layout']=={'focusrite':'right','knob_cup':'left'}
assert m.is_watertight and m.is_winding_consistent and len(m.split())==1
sections=[]
for y in [1,6,12.95,13.05,13.35,13.45,15.35,15.45]:
 circles=[]
 for v in m.section(plane_normal=[0,1,0],plane_origin=[0,y,20]).discrete:
  if np.allclose((v[:,[0,2]].min(0)+v[:,[0,2]].max(0))/2,[cx,cz],atol=.002):circles.append(np.ptp(v[:,[0,2]],axis=0))
 assert any(np.allclose(d,[inner_d,inner_d],atol=.002) for d in circles)==(y<13)
 pocket_air=p['floor_plane_y_mm']<y<c['outside_depth_mm']
 assert any(np.allclose(d,[6,6],atol=.002) for d in circles)==pocket_air,(y,circles)
 if 2.6<y<c['outside_depth_mm']:assert any(np.allclose(d,[outer_d,outer_d],atol=.002) for d in circles)
 sections.append({'depth_y_mm':y,'cup_open':y<13,'magnet_pocket_open':pocket_air,'circle_dimensions_mm':[d.tolist() for d in circles]})
def cy(radius,y0,y1):
 b=trimesh.creation.cylinder(radius=radius,height=y1-y0,sections=128);b.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]));b.apply_translation([cx,(y0+y1)/2,cz]);return b
def volume_intersection(probe):
 x=trimesh.boolean.intersection([m,probe],engine='manifold');return float(x.volume) if len(x.faces) else 0.
proxy=cy(20,.01,12.99);interference=volume_intersection(proxy);assert interference<.01
membrane=cy(2.9,13.01,13.39);assert np.isclose(volume_intersection(membrane),membrane.volume,rtol=1e-5)
pocket=cy(2.98,13.41,15.39);assert volume_intersection(pocket)<.001
# Independently intersect actual triangles with rays through and beside the recess.
def ray_y(x,z):
 tri=m.triangles;origin=np.array([x,-1.,z]);direction=np.array([0.,1.,0.]);e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h);valid=np.abs(det)>1e-10;inv=np.zeros_like(det);inv[valid]=1/det[valid];s=origin-tri[:,0];u=np.einsum('ij,ij->i',s,h)*inv;q=np.cross(s,e1);v=np.einsum('j,ij->i',direction,q)*inv;t=np.einsum('ij,ij->i',e2,q)*inv;valid&=(u>=-1e-8)&(v>=-1e-8)&(u+v<=1+1e-8)&(t>=0);return np.unique(np.round(t[valid]-1,5))
central=ray_y(cx+.137,cz+.211);outside=ray_y(cx+5.137,cz+.211)
assert np.allclose(central,[13,13.4],atol=.002),central
assert np.allclose(outside,[13,15.4],atol=.002),outside
assert np.isclose(outside[-1]-central[-1],2,atol=.002)
probe=trimesh.creation.box([54.124,9,60]);probe.apply_translation([-82.938,20.5,20]);assert volume_intersection(probe)<.01
loops=m.section(plane_normal=[0,1,0],plane_origin=[0,1,20]).discrete
opening=[v for v in loops if np.allclose(np.ptp(v[:,[0,2]],axis=0),[144,44.5],atol=.002)];assert len(opening)==1;center=(opening[0][:,[0,2]].min(0)+opening[0][:,[0,2]].max(0))/2;assert center[0]>0
c.update(actual_sections=sections,nominal_40mm_knob_proxy_interference_mm3=interference,magnet_pocket_air_probe_mm3=volume_intersection(pocket),central_ray_solid_interval_y_mm=central.tolist(),off_center_ray_solid_interval_y_mm=outside.tolist(),actual_membrane_mm=float(np.diff(central)[0]),actual_full_bottom_mm=float(np.diff(outside)[0]),focusrite_opening_center_xz_mm=center.tolist(),fit_not_physically_tested=True)
(STATE/'validation.json').write_text(json.dumps(r,ensure_ascii=False,indent=2));print(json.dumps({'inside_diameter_mm':inner_d,'inside_depth_mm':13,'bottom_mm':c['actual_full_bottom_mm'],'magnet_diameter_mm':6,'magnet_depth_mm':float(outside[-1]-central[-1]),'remaining_membrane_mm':c['actual_membrane_mm'],'watertight':True},indent=2))
