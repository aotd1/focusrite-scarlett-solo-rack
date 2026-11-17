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
m=trimesh.load_mesh(OUTPUT/'service.stl')
old=trimesh.load_mesh(PROJECT/'source/original.stl');old.apply_translation([0,102.5999984741211,-181.125])
cfg=json.loads((ROOT/'parameters.json').read_text());report=json.loads((STATE/'validation.json').read_text())
if cfg.get('mirror_x_for_focusrite_on_right'):old.apply_scale([-1,1,1])
def loops(mesh):return [v[:,[0,2]] for v in mesh.section(plane_normal=[0,1,0],plane_origin=[0,1,20]).discrete]
def dist_to_poly(points,poly):
 a=poly[:-1];ab=poly[1:]-a;valid=(ab*ab).sum(1)>1e-14;a=a[valid];ab=ab[valid]
 ap=points[:,None,:]-a[None,:,:];t=np.clip((ap*ab).sum(2)/(ab*ab).sum(1),0,1)
 return float(np.sqrt(((ap-t[:,:,None]*ab)**2).sum(2)).min(1).max())
newloops=loops(m);oldholes=[v for v in loops(old) if np.ptp(v[:,0])<12 and np.ptp(v[:,1])<12];newholes=[v for v in newloops if np.ptp(v[:,0])<12 and np.ptp(v[:,1])<12]
if cfg.get('knob_cup',{}).get('enabled'):
 oldholes=[v for v in oldholes if abs((v[:,0].min()+v[:,0].max())/2)>110]
 newholes=[v for v in newholes if abs((v[:,0].min()+v[:,0].max())/2)>110]
assert len(oldholes)==len(newholes)==(4 if cfg.get('knob_cup',{}).get('enabled') else 5)
checks=[]
for p in oldholes:
 c=(p.min(0)+p.max(0))/2;q=min(newholes,key=lambda q:np.linalg.norm((q.min(0)+q.max(0))/2-c));err=max(dist_to_poly(p,q),dist_to_poly(q,p));assert err<.002,err
 checks.append({'center_xz_mm':c.tolist(),'profile_max_deviation_mm':err})
outer=max(newloops,key=lambda v:np.ptp(v[:,0]));oldouter=max(loops(old),key=lambda v:np.ptp(v[:,0]))
def bottom_at(poly,x):
 zs=[]
 for a,b in zip(poly[:-1],poly[1:]):
  if abs(b[0]-a[0])>1e-8 and min(a[0],b[0])-1e-8<=x<=max(a[0],b[0])+1e-8:zs.append(a[1]+(x-a[0])*(b[1]-a[1])/(b[0]-a[0]))
 return min(zs)
width=cfg.get('lower_extension_width_mm',240);half=width/2
bottom=[bottom_at(outer,x) for x in np.linspace(-half+.01,half-.01,49)]
assert np.ptp(bottom)<.002 and np.isclose(bottom[0],-report['opening_enlarged_down_mm'],atol=.002)
# At both side ears, including rounded corners, the lower outline is stock.
ear_samples=[-127,-126,-125,-120,-110,110,120,125,126,127]
ear_error=max(abs(bottom_at(outer,x)-bottom_at(oldouter,x)) for x in ear_samples)
assert ear_error<.002,ear_error
bottom_vertices=outer[np.isclose(outer[:,1],-report['opening_enlarged_down_mm'],atol=.002)]
assert np.allclose([bottom_vertices[:,0].min(),bottom_vertices[:,0].max()],[-half,half],atol=.002)
assert m.is_watertight and m.is_winding_consistent and len(m.split())==1
assert np.allclose(m.extents,[255.050003,cfg['focusrite_depth_mm'],47.100006],atol=.002)
report.update(full_width_lower_panel=False,lower_extension_width_mm=width,stock_side_outline_max_deviation_mm=ear_error,front_panel_depth_mm=cfg['front_panel_depth_mm'],existing_mounting_hole_profiles=checks,lower_edge_verified_x_range_mm=[-half,half],lower_edge_z_mm=float(bottom[0]),lower_edge_max_deviation_mm=float(np.ptp(bottom)))
(STATE/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps({'holes_unchanged':len(checks),'lower_extension_width_mm':width,'stock_side_outline_max_deviation_mm':ear_error,'watertight':True},indent=2))
