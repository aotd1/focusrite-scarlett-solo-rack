# CC BY-NC 4.0. Original model: Jonathan Lemoine; see ATTRIBUTION.md.
"""Rebuild both centered variants directly from the original STL. Units: mm."""
from pathlib import Path
import io
import hashlib
import json
import zipfile
from xml.etree import ElementTree as ET
import numpy as np
import trimesh
import manifold3d

ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUTPUT=PROJECT/'models'
STATE=PROJECT/'build/centered'
STATE.mkdir(parents=True,exist_ok=True)
SOURCE=PROJECT/'source/original.stl'
SOURCE_SHA='2d9ca6b355928fdd7d8014ae271b9dde45f5652f91894ac32cf06bf1c6dd52a4'
STEMS={'centered-solid':'Scarlett_Solo_4Gen_Centered_D84p2',
       'centered-front-vents':'Scarlett_Solo_4Gen_Centered_FrontVents_D84p2'}

def box(lo,hi):
    lo,hi=np.array(lo,dtype=float),np.array(hi,dtype=float)
    mesh=trimesh.creation.box(hi-lo);mesh.apply_translation((hi+lo)/2)
    return mesh

def boolean(operation,meshes):
    result=getattr(trimesh.boolean,operation)(meshes,engine='manifold')
    assert result.is_watertight and result.is_winding_consistent
    return result

def extrude_section(cross,depth):
    raw=cross.extrude(depth).to_mesh()
    mesh=trimesh.Trimesh(vertices=raw.vert_properties[:,:3],faces=raw.tri_verts,process=True)
    mesh.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
    mesh.apply_translation([0,depth,0])
    return mesh

def capsule(center,start,end,axis,span):
    """Rounded slot along Y through Z or X (rear sleeve)."""
    radius=3.;lo,hi=span
    if axis=='z':
        pieces=[box([center-radius,start+radius,lo],[center+radius,end-radius,hi])]
        for y in [start+radius,end-radius]:
            c=trimesh.creation.cylinder(radius=radius,height=hi-lo,sections=32)
            c.apply_translation([center,y,(hi+lo)/2]);pieces.append(c)
    else:
        pieces=[box([lo,start+radius,center-radius],[hi,end-radius,center+radius])]
        for y in [start+radius,end-radius]:
            c=trimesh.creation.cylinder(radius=radius,height=hi-lo,sections=32)
            c.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[0,1,0]))
            c.apply_translation([(lo+hi)/2,y,center]);pieces.append(c)
    return boolean('union',pieces)

def front_slot(x,z,width,length,depth):
    radius=width/2
    pieces=[box([x-radius,-1,z-length/2+radius],[x+radius,depth+1,z+length/2-radius])]
    for zz in [z-length/2+radius,z+length/2-radius]:
        c=trimesh.creation.cylinder(radius=radius,height=depth+2,sections=64)
        c.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
        c.apply_translation([x,depth/2,zz]);pieces.append(c)
    return boolean('union',pieces)

def export_model(mesh,directory,stem):
    directory.mkdir(exist_ok=True)
    mesh.export(directory/'service.stl')
    printed=mesh.copy()
    printed.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
    printed.apply_transform(trimesh.transformations.rotation_matrix(np.pi/4,[0,0,1]))
    printed.apply_translation(-printed.bounds[0])
    printed.apply_translation([(256-printed.extents[0])/2,(256-printed.extents[1])/2,0])
    printed.vertices=printed.vertices.astype(np.float32).astype(np.float64)
    printed.merge_vertices(digits_vertex=5)
    printed.update_faces(printed.nondegenerate_faces())
    printed.update_faces(printed.unique_faces())
    printed.remove_unreferenced_vertices()
    assert printed.is_watertight and printed.is_winding_consistent
    assert len(printed.split())==1 and (printed.extents[:2]<250).all()
    printed.export(directory/(directory.name+'.stl'))
    raw=trimesh.exchange.threemf.export_3MF(trimesh.Scene(printed))
    metadata={
        'Title':stem.replace('_',' '),
        'Designer':'Remix of Jonathan Lemoine (@JonathanLemoi_601472)',
        'Description':'Centered Scarlett Solo 4th Gen; feet removed; no knob cup; physical fit untested.',
        'LicenseTerms':'CC BY-NC 4.0 - https://creativecommons.org/licenses/by-nc/4.0/ ; recorded source license; live recheck pending',
        'Source':'https://www.printables.com/model/1565705-focusrite-scarlettsolo-3rd-gen-10-inch-rack-mount',
    }
    ns='http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
    buffer=io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(raw)) as zin,zipfile.ZipFile(buffer,'w',compression=zipfile.ZIP_DEFLATED) as zout:
        for name in zin.namelist():
            data=zin.read(name)
            if name=='3D/3dmodel.model':
                xml=ET.fromstring(data)
                for i,(key,value) in enumerate(metadata.items()):
                    tag=ET.Element('{'+ns+'}metadata',{'name':key});tag.text=value;xml.insert(i,tag)
                data=ET.tostring(xml,encoding='utf-8',xml_declaration=True)
            zout.writestr(name,data)
    (directory/(directory.name+'.3mf')).write_bytes(buffer.getvalue())
    return dict(watertight=bool(mesh.is_watertight),winding_consistent=bool(mesh.is_winding_consistent),
                connected_components=len(mesh.split()),bounds_mm=mesh.bounds.tolist(),size_mm=mesh.extents.tolist(),
                volume_mm3=float(mesh.volume),print_size_mm=printed.extents.tolist(),
                service_sha256=hashlib.sha256((directory/'service.stl').read_bytes()).hexdigest())

def build():
    cfg=json.loads((ROOT/'parameters.json').read_text())
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
    src=trimesh.load_mesh(SOURCE)
    src.apply_translation([0,102.5999984741211,-181.125])
    xl,xli,xri,xr=-94.125,-91.125,52.875,55.875
    floor,roof=1.3000030517578125,42.69999694824219
    delta=cfg['body_height_mm']+cfg['vertical_clearance_mm']-(roof-floor)
    depth=float(src.bounds[1,1]);d=cfg['front_panel_depth_mm']
    assert abs(depth-cfg['depth_mm'])<1e-5 and cfg['remove_rubber_feet']
    # Same downward opening enlargement and rear thinning as the cup variant.
    roi=box([xl,-1,-20],[xr,depth+1,20])
    lower=boolean('intersection',[src,roi]);lower.apply_translation([0,0,-delta])
    upper=boolean('difference',[src,roi])
    bridges=[box([xl,0,20-delta],[xli,depth,20]),box([xri,0,20-delta],[xr,depth,20])]
    modified=boolean('union',[lower,upper,*bridges])
    wall=cfg['rear_side_wall_mm'];band=cfg['front_reinforced_depth_mm']
    modified=boolean('difference',[modified,
        box([xl-.01,band,-20],[xli-wall,depth+1,60]),
        box([xri+wall,band,-20],[xr+.01,depth+1,60])])
    # Detach the actual sleeve from its asymmetric fascia, discard all fan pod
    # geometry, then translate the complete sleeve rather than only its opening.
    # Clip on the exact side planes: extending this crop even 0.001 mm would
    # retain a thin wafer of the adjacent original fan-controller carrier.
    sleeve=boolean('intersection',[modified,box([xl,d,-20],[xr,depth+1,60])])
    shift=-(xli+xri)/2
    sleeve.apply_translation([shift,0,0])
    # Keep the original rear ear backing, in its original rack coordinates.
    ear_source=boolean('difference',[src,box([xr,d,-10],[110,26,60])])
    ears=[boolean('intersection',[ear_source,box([a,d,-10],[b,depth+1,60])])
          for a,b in [(-140,-106.5),(106.5,140)]]
    loops=src.section(plane_normal=[0,1,0],plane_origin=[0,1,20]).discrete
    outer=max(loops,key=lambda v:np.ptp(v[:,0]))[:,[0,2]]
    mounts=[v[:,[0,2]] for v in loops if np.ptp(v[:,0])<12 and abs((v[:,0].min()+v[:,0].max())/2)>110]
    assert len(mounts)==4
    loops=modified.section(plane_normal=[0,1,0],plane_origin=[0,1,20]).discrete
    openings=[v[:,[0,2]].copy() for v in loops if np.allclose(np.ptp(v[:,[0,2]],axis=0),[144,44.5],atol=.002)]
    assert len(openings)==1
    opening=openings[0];opening[:,0]+=shift
    half=cfg['lower_extension_width_mm']/2
    extension=np.array([[-half,-delta],[half,-delta],[half,4],[-half,4]])
    cs=lambda polys:manifold3d.CrossSection(polys,manifold3d.FillRule.EvenOdd)
    panel=extrude_section(cs([outer])+cs([extension])-cs([*mounts,opening]),d)
    # Cut the same 56 sleeve slots. Coordinates shifted with the centered body.
    vent=cfg['sleeve_slots'];assert vent['width_mm']==6 and vent['length_mm']==20
    tools=[]
    for y0,y1 in vent['depth_ranges_mm']:
        for x in np.linspace(xli+12,xri-12,vent['top_bottom_columns'])+shift:
            tools.append(capsule(x,y0,y1,'z',(-8,50)))
        for z in np.linspace(floor-delta+10,roof-10,vent['side_columns']):
            tools.append(capsule(z,y0,y1,'x',(-80,80)))
    sleeve=boolean('difference',[sleeve,*tools])
    plain=boolean('union',[panel,sleeve,*ears])
    # Retain the previous layout's complete X reflection, including ear backing.
    plain.apply_scale([-1,1,1])
    front=cfg['front_slots']
    front_tools=[front_slot(x,front['center_z_mm'],front['width_mm'],front['length_mm'],d)
                 for x in front['centers_x_mm']]
    slotted=boolean('difference',[plain,*front_tools])
    report=dict(source_sha256=SOURCE_SHA,units='mm',opening_mm=[144,44.5],opening_center_x_mm=0,
                depth_datum='Y=0 front face to original rear edge',source_depth_mm=depth,
                lower_extension_mm=delta,sleeve_vent_count=56,front_vent_count=6,
                parameters=cfg,physical_fit_tested=False,print_tested=False,variants={})
    for key,mesh in [('centered-solid',plain),('centered-front-vents',slotted)]:
        assert mesh.is_watertight and mesh.is_winding_consistent and len(mesh.split())==1
        report['variants'][key]=export_model(mesh,OUTPUT/key,STEMS[key])
    (STATE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':build()
