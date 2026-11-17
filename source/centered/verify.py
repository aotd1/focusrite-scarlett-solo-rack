# CC BY-NC 4.0. Original model: Jonathan Lemoine; see ATTRIBUTION.md.
"""Measure final mesh sections and compare variants and original mounting geometry."""
from pathlib import Path
import json
from collections import Counter
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from build import STEMS, box

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUTPUT=PROJECT/'models'
STATE=PROJECT/'build/centered'
STATE.mkdir(parents=True,exist_ok=True)

def clip(mesh,region):return trimesh.boolean.intersection([mesh,region],engine='manifold')
def volume(mesh):return abs(float(mesh.volume)) if len(mesh.faces) else 0.
def difference(a,b):return volume(trimesh.boolean.difference([a,b],engine='manifold'))
def symmetric_difference(a,b):return difference(a,b)+difference(b,a)
def loops(mesh,y):return mesh.section(plane_normal=[0,1,0],plane_origin=[0,y,20]).discrete
def profile_error(p,q):
    a=q[:-1];v=q[1:]-a;valid=np.linalg.norm(v,axis=1)>1e-9;a=a[valid];v=v[valid]
    w=p[:,None]-a[None]
    t=np.clip((w*v).sum(2)/(v*v).sum(1),0,1)
    return float(np.linalg.norm(w-t[:,:,None]*v,axis=2).min(1).max())

def verify():
    report=json.loads((STATE/'validation.json').read_text())
    meshes={key:trimesh.load_mesh(OUTPUT/key/'service.stl') for key,stem in STEMS.items()}
    source=trimesh.load_mesh(PROJECT/'source/original.stl')
    source.apply_translation([0,102.5999984741211,-181.125]);source.apply_scale([-1,1,1])
    source_mounts=[v[:,[0,2]] for v in loops(source,1) if np.ptp(v[:,0])<12 and abs((v[:,0].min()+v[:,0].max())/2)>110]
    assert len(source_mounts)==4
    for key,m in meshes.items():
        assert m.is_watertight and m.is_winding_consistent and len(m.split())==1 and m.volume>0
        assert np.allclose(m.extents,[255.05,84.2,47.1],atol=.002)
        assert np.allclose(m.bounds[:,1],[0,84.2],atol=.00002)
        centers=[]
        for y in [.1,1,2.5,6,12,27,82.2,84.19]:
            openings=[v for v in loops(m,y) if np.allclose(np.ptp(v[:,[0,2]],axis=0),[144,44.5],atol=.002)]
            assert len(openings)==1,(key,y)
            opening=openings[0]
            center=(opening[:,0].min()+opening[:,0].max())/2
            assert abs(center)<.002
            centers.append(dict(y_mm=y,opening_mm=[144,44.5],center_x_mm=float(center)))
        contours=loops(m,1)
        mounts=[v[:,[0,2]] for v in contours if np.ptp(v[:,0])<12 and abs((v[:,0].min()+v[:,0].max())/2)>110]
        assert len(mounts)==4
        errors=[]
        for p in source_mounts:
            q=min(mounts,key=lambda q:np.linalg.norm(q.mean(0)-p.mean(0)))
            e=max(profile_error(p,q),profile_error(q,p));assert e<.002;errors.append(e)
        # Neither knob cup nor fan carrier remains behind either front wing.
        no_pod=0.
        for xlo,xhi in [(-106.49,-75.01),(75.01,106.49)]:
            no_pod+=volume(clip(m,box([xlo,2.6001,-10],[xhi,84.21,60])))
        assert no_pod<.001,(key,no_pod)
        # Actual wall sections verify all 56 rounded through-slots at new positions.
        counts={}
        for face,axis,coord,axes,expected in [('top',2,43.15,[0,1],22),('bottom',2,-2.25,[0,1],22),
                                              ('left',0,-73,[1,2],6),('right',0,73,[1,2],6)]:
            origin=np.zeros(3);origin[axis]=coord
            found=[v for v in m.section(np.eye(3)[axis],origin).discrete
                   if np.allclose(sorted(np.ptp(v[:,axes],axis=0)),[6,20],atol=.002)]
            assert len(found)==expected,(key,face,len(found))
            for v in found:
                assert tuple(np.round([v[:,1].min(),v[:,1].max()],3)) in [(30.,50.),(56.,76.)]
                p=(v.min(0)+v.max(0))/2;half=np.full(3,.1);half[axis]=2.5
                assert volume(clip(m,box(p-half,p+half)))<1e-5
            counts[face]=len(found)
        # Check 2 mm rear walls from full closed contours and original end plane.
        rear=loops(m,82.2);inner=min(rear,key=lambda v:np.ptp(v[:,0]));outer=max(rear,key=lambda v:np.ptp(v[:,0]))
        walls=[inner[:,0].min()-outer[:,0].min(),outer[:,0].max()-inner[:,0].max()]
        assert np.allclose(walls,[2,2],atol=.002)
        actual_front=[]
        for y in [.1,1,2.5]:
            slots=[v for v in loops(m,y) if np.allclose(np.ptp(v[:,[0,2]],axis=0),[6,26],atol=.002)]
            assert len(slots)==(6 if key.endswith('front-vents') else 0)
            for v in slots:
                center=(v[:,[0,2]].min(0)+v[:,[0,2]].max(0))/2
                assert min(abs(center[0]-np.array([-102,-92,-82,82,92,102])))<.002
                assert abs(center[1]-20.449996948242188)<.002
            actual_front.append(dict(y_mm=y,slots=len(slots)))
        stem=STEMS[key]
        printed=trimesh.load_mesh(OUTPUT/key/(key+'.stl'))
        mf=trimesh.load(OUTPUT/key/(key+'.3mf')).to_mesh()
        for exported in [printed,mf]:
            assert exported.is_watertight and exported.is_winding_consistent and len(exported.split())==1
            assert np.isclose(exported.volume,m.volume,rtol=1e-5)
            assert np.allclose(exported.extents,printed.extents,atol=.002)
        assert np.isclose(printed.bounds[0,2],0,atol=.00002)
        assert np.isclose(printed.extents[2],84.2,atol=.00002)
        assert np.all(printed.bounds[0]>=0) and np.all(printed.bounds[1,:2]<=256)
        assert cKDTree(printed.vertices).query(mf.vertices)[0].max()<.00002
        report['variants'][key]['independent_checks']=dict(passed=True,opening_sections=centers,
            mounting_profile_errors_mm=errors,no_pod_material_mm3=no_pod,ventilation_counts=counts,
            rear_side_walls_mm=[float(x) for x in walls],front_slot_sections=actual_front,
            stl_3mf_readback_passed=True)
    plain,vented=meshes.values()
    assert difference(vented,plain)<.001
    rear_region=box([-140,2.6001,-10],[140,85,60])
    rear_error=symmetric_difference(clip(plain,rear_region),clip(vented,rear_region))
    assert rear_error<.001
    ear_errors=[]
    for a,b in [(-140,-106.501),(106.501,140)]:
        region=box([a,-.01,-10],[b,85,60])
        error=symmetric_difference(clip(plain,region),clip(vented,region));assert error<.001;ear_errors.append(error)
    # Full-thickness air at each front slot and retained 4 mm solid ribs.
    for x in [-102,-92,-82,82,92,102]:
        assert volume(clip(vented,box([x-.1,-.1,20.3],[x+.1,2.7,20.6])))<1e-5
    for x in [-97,-87,87,97]:
        p=box([x-.2,.1,20],[x+.2,2.5,21])
        assert np.isclose(volume(clip(vented,p)),p.volume,atol=1e-5)
    assert vented.euler_number==plain.euler_number-12
    report['variant_comparison']=dict(passed=True,only_front_material_removed=True,
        rear_geometry_difference_mm3=rear_error,ear_geometry_differences_mm3=ear_errors,
        added_front_holes=6,front_rib_nominal_width_mm=4,
        added_front_open_area_mm2=float((plain.volume-vented.volume)/2.6),
        cooling_performance_measured=False)
    (STATE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['variant_comparison'],indent=2))

if __name__=='__main__':verify()
