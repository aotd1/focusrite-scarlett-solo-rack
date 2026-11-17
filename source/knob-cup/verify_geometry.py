# CC BY-NC 4.0. Original model: Jonathan Lemoine.
# Modified mesh source/verification/rendering; see ATTRIBUTION.md and LICENSE.md.
"""Independent section, solid-probe, export and source-profile checks (mm).

Run after prepare_remix.py, verify_panel_profile.py and verify_knob_cup.py.
No previous model or printer connection is needed.
"""
from pathlib import Path
from collections import Counter
import json
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import trimesh
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUTPUT=PROJECT/'models/knob-cup'
STATE=PROJECT/'build/knob-cup'
OUTPUT.mkdir(parents=True,exist_ok=True)
STATE.mkdir(parents=True,exist_ok=True)
m = trimesh.load_mesh(OUTPUT / 'service.stl')
r = json.loads((STATE / 'validation.json').read_text())
source = trimesh.load_mesh(PROJECT / 'source/original.stl')
source.apply_translation([0, 102.5999984741211, -181.125])
source.apply_scale([-1, 1, 1])
depth = float(source.bounds[1, 1])
assert np.allclose(m.bounds[:, 1], [0, depth], atol=2e-5)
assert abs(depth - 84.2) < 2e-5
assert np.allclose(m.extents, [255.05, 84.2, 47.1], atol=.002)
assert m.is_watertight and m.is_winding_consistent and m.volume > 0
assert len(m.split()) == 1

def box(lo, hi):
    lo, hi = np.array(lo), np.array(hi)
    b = trimesh.creation.box(hi - lo)
    b.apply_translation((hi + lo) / 2)
    return b

def material(probe):
    a = trimesh.boolean.intersection([m, probe], engine='manifold')
    return abs(float(a.volume)) if len(a.faces) else 0.

def distance(points, poly):
    a, ab = poly[:-1], poly[1:] - poly[:-1]
    valid = np.linalg.norm(ab, axis=1) > 1e-9
    a, ab = a[valid], ab[valid]
    ap = points[:, None] - a[None]
    t = np.clip((ap * ab).sum(2) / (ab * ab).sum(1), 0, 1)
    return float(np.linalg.norm(ap - t[:, :, None] * ab, axis=2).min(1).max())

# Compare the actual original shaped end to the restored end. Only the outer
# side walls are deliberately thinned. The horizontal skins retain the source
# planes with the lower half translated by the known extension. The previous
# remix rebuilt the inner side-wall transitions; these are checked against the
# side-wall transitions are modified; they are not claimed stock.
rear_checks = []
for y in [depth - 1.5, depth - .5, depth - .01]:
    loops = m.section([0, 1, 0], [0, y, 20]).discrete
    inner = min(loops, key=lambda v: np.ptp(v[:, 0]))
    outer = max(loops, key=lambda v: np.ptp(v[:, 0]))
    assert np.allclose(np.ptp(inner[:, [0, 2]], axis=0), [144, 44.5], atol=.002)
    walls = [inner[:, 0].min() - outer[:, 0].min(), outer[:, 0].max() - inner[:, 0].max()]
    assert np.allclose(walls, [2, 2], atol=.002)
    source_loops = source.section([0, 1, 0], [0, y, 20]).discrete
    src_inner = min(source_loops, key=lambda v: np.ptp(v[:, 0]))[:, [0, 2]].copy()
    src_inner[src_inner[:, 1] < 20, 1] -= r['opening_enlarged_down_mm']
    actual = inner[:, [0, 2]]
    actual_horizontal = actual[(actual[:, 0] > -45) & (actual[:, 0] < 83)]
    source_horizontal = src_inner[(src_inner[:, 0] > -45) & (src_inner[:, 0] < 83)]
    # Add interpolated middle points: triangulation may leave no interior
    # vertices on long straight edges.
    for poly, points in [(actual, 'actual'), (src_inner, 'source')]:
        midpoints = (poly[1:] + poly[:-1]) / 2
        selected = midpoints[(midpoints[:, 0] > -45) & (midpoints[:, 0] < 83)]
        if points == 'actual': actual_horizontal = np.vstack([actual_horizontal, selected])
        else: source_horizontal = np.vstack([source_horizontal, selected])
    assert len(actual_horizontal) >= 2 and len(source_horizontal) >= 2
    deviation = max(distance(actual_horizontal, src_inner), distance(source_horizontal, actual))
    assert deviation < .002, (y, deviation)
    rear_checks.append(dict(y_mm=y, opening_mm=np.ptp(actual, axis=0).tolist(),
                            side_walls_mm=walls, source_horizontal_profile_error_mm=deviation))
end_faces = m.triangles[np.all(np.isclose(m.triangles[:, :, 1], depth, atol=.00002), axis=1)]
assert len(end_faces) > 0

# Inspect actual closed slot contours on all four skins. Verify positions and
# dimensions from geometry, independently of the builder's location report.
sections = []
holes = []
for face, axis, coord, axes, count in [
    ('top', 2, 43.15, [0, 1], 22),
    ('bottom', 2, -2.25, [0, 1], 22),
    ('left', 0, -53.875, [1, 2], 6),
    ('right', 0, 92.125, [1, 2], 6),
]:
    origin = np.zeros(3)
    origin[axis] = coord
    loops = m.section(np.eye(3)[axis], origin).discrete
    vents = [v for v in loops if np.allclose(sorted(np.ptp(v[:, axes], axis=0)), [6, 20], atol=.002)]
    assert len(vents) == count, (face, len(vents))
    found_ranges = set()
    for v in vents:
        lo, hi = v.min(0), v.max(0)
        yrange = tuple(np.round([lo[1], hi[1]], 3))
        assert yrange in [(30., 50.), (56., 76.)]
        found_ranges.add(yrange)
        center = (lo + hi) / 2
        # Probe the full wall thickness at the slot center and near both ends;
        # these must all be air, not merely shallow decorative recesses.
        for y in [lo[1] + 1, center[1], hi[1] - 1]:
            p = center.copy()
            p[1] = y
            half = np.full(3, .12)
            half[axis] = 2.5
            assert material(box(p - half, p + half)) < 1e-5
        holes.append(dict(face=face, center_mm=center.tolist(), y_range_mm=list(yrange)))
    assert found_ranges == {(30., 50.), (56., 76.)}
    sections.append(dict(face=face, actual_closed_slot_profiles=len(vents), slot_mm=[6, 20]))
assert Counter(h['face'] for h in holes) == {'top':22, 'bottom':22, 'left':6, 'right':6}
assert m.euler_number == r['ventilation']['before_euler'] - 112

# Reinforced front and unvented rear band must retain solid walls. Probe all
# four skins at multiple depths, using source-derived surface coordinates.
for y in [3, 12, 23.9, 27, 78, 82, depth - .1]:
    for p in [[10,y,43.15], [10,y,-2.25], [-53.875,y,20], [92.125,y,20]]:
        p=np.array(p)
        probe=box(p-.1,p+.1)
        assert np.isclose(material(probe), probe.volume, atol=1e-5), (y,p)

# Read exported files back; compare surface samples after reversing placement.
printed=trimesh.load_mesh(OUTPUT/'knob-cup.stl')
scene=trimesh.load(OUTPUT/'knob-cup.3mf')
mf=scene.to_mesh()
expected=m.copy()
expected.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
expected.apply_transform(trimesh.transformations.rotation_matrix(np.pi/4,[0,0,1]))
expected.apply_translation(-expected.bounds[0])
expected.apply_translation([(256-expected.extents[0])/2,(256-expected.extents[1])/2,0])
def vertex_error(a,b):
    return float(max(cKDTree(a.vertices).query(b.vertices)[0].max(),
                     cKDTree(b.vertices).query(a.vertices)[0].max()))
placement_error=vertex_error(expected,printed)
format_error=vertex_error(printed,mf)
assert placement_error<.0001 and format_error<.0001
for exported in [printed,mf]:
    assert exported.is_watertight and exported.is_winding_consistent and exported.volume>0
    assert len(exported.split())==1
    assert np.allclose(exported.extents, printed.extents, atol=.002)
    assert np.isclose(exported.volume,m.volume,rtol=1e-5)
assert np.all(printed.bounds[0]>=0) and np.all(printed.bounds[1,:2]<256)
assert np.isclose(printed.bounds[0,2],0,atol=.00002)
assert np.isclose(printed.extents[2],depth,atol=.00002)
with zipfile.ZipFile(OUTPUT/'knob-cup.3mf') as z:
    model=ET.fromstring(z.read('3D/3dmodel.model'))
    assert model.attrib['unit']=='millimeter'
    metadata={x.attrib['name']:x.text for x in model if x.tag.endswith('metadata')}
    assert 'Jonathan Lemoine' in metadata['Designer']
    assert 'CC BY-NC 4.0' in metadata['LicenseTerms']
    assert not any('gcode' in n.lower() for n in z.namelist())
r['independent_geometry_checks'] = dict(
    passed=True, rear_sections=rear_checks, rear_end_face_triangles=len(end_faces),
    ventilation_sections=sections, holes_measured_from_mesh=holes,
    through_wall_air_probes=168, solid_band_probe_count=28,
    print_and_3mf_readback_passed=True, source_end_plane_mm=depth,
    stl_units='millimetres by convention; 3MF declares millimeter',
    print_extents_mm=printed.extents.tolist(),
    print_placement_vertex_error_mm=placement_error, stl_3mf_vertex_error_mm=format_error,
)
(STATE/'validation.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r['independent_geometry_checks'],indent=2))
