"""Check the distributed meshes without rebuilding or changing model files."""
from pathlib import Path
import numpy as np
import trimesh

root=Path(__file__).resolve().parents[1]
for variant,cx in [('knob-cup',19.125),('centered-solid',0),('centered-front-vents',0)]:
    directory=root/'models'/variant
    service=trimesh.load_mesh(directory/'service.stl')
    printed=trimesh.load_mesh(directory/(variant+'.stl'))
    packed=trimesh.load(directory/(variant+'.3mf')).to_mesh()
    for mesh in [service,printed,packed]:
        assert mesh.is_watertight and mesh.is_winding_consistent
        assert len(mesh.split())==1 and mesh.volume>0
        assert np.isclose(mesh.volume,service.volume,rtol=1e-5)
    assert np.allclose(service.extents,[255.05,84.2,47.1],atol=.002)
    assert np.allclose(service.bounds[:,1],[0,84.2],atol=.00002)
    assert np.allclose(printed.extents,packed.extents,atol=.002)
    assert np.isclose(printed.bounds[0,2],0,atol=.00002)
    for y in [1,12,27,82.2]:
        loops=service.section([0,1,0],[0,y,20]).discrete
        openings=[v for v in loops if np.allclose(np.ptp(v[:,[0,2]],axis=0),[144,44.5],atol=.002)]
        assert len(openings)==1
        assert np.isclose((openings[0][:,0].min()+openings[0][:,0].max())/2,cx,atol=.002)
    print('PASS:',variant)
