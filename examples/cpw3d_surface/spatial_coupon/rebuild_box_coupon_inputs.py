#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Rebuild complete box-boundary inputs from a retained campaign model.

All retained nodal degrees of freedom survive, with explicit index mappings.
Missing corner degrees of freedom are added, never approximated or truncated.
The original configs, sources, material geometry, solver orders and tolerances
are retained unchanged; corrected configs are written into a new directory.
"""
import argparse
import copy
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from box_trace import complete_box_trace
from export_legacy_trace_mesh import trace_values
from generate_spatial_response import (
    conductor_at_points, conductor_trace_lifts, open_paths,
    write_basis, write_surface_trace,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_geometry(root):
    edges = []
    with (root/'mesh-signature.csv').open() as stream:
        for row in csv.DictReader(stream):
            r = {k:float(v) for k,v in row.items()}
            edges.append({"Conductor":int(r['Conductor']), "Point":[r['Px'],r['Py'],r['Pz']],
                          "Tangent":[r['Tx'],r['Ty'],r['Tz']], "GapDirection":[r['Gx'],r['Gy'],r['Gz']],
                          "ProcessNormal":[0.,0.,r['Nz']], "Interval":[r['S0'],r['S1']],
                          "VertexArm":bool(r['VertexArm'])})
    facets = {}
    with (root/'plan-view-mask.csv').open() as stream:
        for r in csv.DictReader(stream):
            key=int(r['Facet'])
            item=facets.setdefault(key,{"Conductor":int(r['Conductor']), "Plane":float(r['Plane']),"Points":[]})
            item['Points'].append([float(r['X']),float(r['Y'])])
    return edges,list(facets.values())


def rebuild(reference_manifest, model_name, output):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError(f"Refuse to reuse input directory {output}")
    reference = json.loads(Path(reference_manifest).read_text())
    cases = {c['Kind']:c for c in reference['Cases'] if c['Model']==model_name}
    if set(cases) != {'thin','fabricated'}:
        raise ValueError("Expected both retained model cases")
    configs = {kind:json.loads(Path(c['SourceConfig']).read_text()) for kind,c in cases.items()}
    source_root = Path(cases['thin']['SourceConfig']).parent
    source_library = Path(cases['thin']['SourceLibrary'])
    library = json.loads(source_library.read_text())
    model = copy.deepcopy(next(m for m in library['Models'] if m['Name']==model_name))
    basis_path = source_library.parent/model['BasisPoints']
    canonical = np.loadtxt(basis_path, delimiter=',', skiprows=1, ndmin=2)
    nold = len(canonical)
    old_sources = configs['thin']['Boundaries']['PrescribedPotential']
    if [s['Index'] for s in old_sources] != list(range(1,len(old_sources)+1)):
        raise ValueError("Expected contiguous retained source indices")
    for source, other in zip(old_sources, configs['fabricated']['Boundaries']['PrescribedPotential'], strict=True):
        if any(source.get(key) != other.get(key) for key in ('Index','Attributes','DataFile')):
            raise ValueError("Thin/fabricated source trace contracts differ")
    values = [trace_values(Path(s['DataFile'])) for s in old_sources]
    raw_basis = []
    for table in values[:nold]:
        ones = [p for p,v in table.items() if abs(v-1.)<1e-12]
        if len(ones)!=1:
            raise ValueError("Retained nodal source has no unique unit vertex")
        raw_basis.append(ones[0])
    raw_basis=np.asarray(raw_basis)
    design=np.column_stack((raw_basis,np.ones(nold)))
    transform=np.linalg.lstsq(design,canonical,rcond=None)[0]
    residual=float(np.max(np.abs(design@transform-canonical)))
    if residual>1e-8 or np.max(np.abs(transform[:3]@transform[:3].T-np.eye(3)))>1e-8:
        raise ValueError("Retained coordinate transform is not a rigid frame")
    retained_points=sorted(values[0])
    if any(set(v)!=set(retained_points) for v in values):
        raise ValueError("Retained sources use different geometric vertices")
    points,triangles,groups,geometry_report=complete_box_trace(retained_points)
    edges,facets=read_geometry(source_root)
    radius=float(library['MatchingRadius'])
    fabrication=library['Fabrication']
    process={"Units":fabrication['LengthUnit'],"Radius":radius,
             "MetalThickness":float(fabrication['MetalThickness']),
             "Overetch":float(fabrication['OveretchDepth']),
             "SidewallAngle":float(fabrication['SidewallAngleDegrees']),
             "TopRounding":float(fabrication['TopRoundingRadius']),
             "TrenchRounding":float(fabrication['BottomRoundingRadius'])}
    if (process['Units']!='um' or process['SidewallAngle']!=90. or
            process['TopRounding']!=0. or process['TrenchRounding']!=0. or
            not all(np.isfinite(value) for name,value in process.items() if name!='Units') or
            not 0<=process['Overetch']<radius or process['MetalThickness']<=0):
        raise ValueError('This compatibility adapter requires an explicit sharp vertical um process')
    labels=conductor_at_points(points,edges,radius,process['MetalThickness'],process['SidewallAngle'],facets)
    conductors=max(edge['Conductor'] for edge in edges)
    if len(old_sources)-nold != conductors-1:
        raise ValueError("Unexpected retained conductor-state count")
    active=np.arange(len(points)) if conductors==1 else np.flatnonzero(labels==0)
    active_lookup={int(v):i+1 for i,v in enumerate(active)}
    # Match old nodes through the explicit tolerance-bounded corner snap map.
    old_geometry={point:i for i,point in enumerate(retained_points)}
    old_to_new={}
    for index,point in enumerate(raw_basis,start=1):
        vertex=geometry_report['RetainedVertexMap'][old_geometry[tuple(point)]]
        if vertex not in active_lookup:
            raise ValueError("A retained free basis node became conductor constrained")
        old_to_new[index]=active_lookup[vertex]
    if len(set(old_to_new.values()))!=nold:
        raise ValueError("A retained degree of freedom was dropped")
    zero=(np.flatnonzero(labels)+1).tolist() if conductors==1 else []
    if conductors==1:
        if {old_to_new[i] for i in model.get('ZeroTraceIndices',[])} != set(zero)&set(old_to_new.values()):
            raise ValueError("Retained zero-trace constraints changed")
    paths=open_paths(groups,labels,active.tolist(),conductors) if conductors>1 else []
    output.mkdir(parents=True)
    for name in ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv'):
        (output/name).write_bytes((source_root/name).read_bytes())
    (output/'process.toml').write_text('\n'.join(f'{k} = {json.dumps(v)}' for k,v in process.items())+'\n')
    traces,zero_trace=write_basis(output,points,triangles,active.tolist())
    transformed=np.column_stack((points,np.ones(len(points))))@transform
    np.savetxt(output/'basis-points.csv',transformed[active],delimiter=',',header='x,y,z',comments='',fmt='%.16e')
    with (output/'trace-vertices.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(['vertex','x','y','z','basis','conductor'])
        for i,p in enumerate(transformed):writer.writerow([i+1,*[f'{x:.16e}' for x in p],active_lookup.get(i,0),int(labels[i])])
    with (output/'trace-triangles.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(['triangle','vertex_i','vertex_j','vertex_k'])
        for i,t in enumerate(triangles,start=1):writer.writerow([i,*[int(v)+1 for v in t]])
    lift_paths={}
    for conductor,lift in conductor_trace_lifts(points,labels,conductors).items():
        path=output/f'conductor-{conductor}.csv';write_surface_trace(path,points,triangles,lift);lift_paths[conductor]=path
    for key in ('ContourGroups','OpenContourPaths','ZeroTraceIndices'):
        model.pop(key,None)
    if paths:model['OpenContourPaths']=paths
    else:
        model['ContourGroups']=groups
        if zero:model['ZeroTraceIndices']=zero
    model['BasisPoints']='basis-points.csv'
    model['TraceMesh']={'Vertices':'trace-vertices.csv','Triangles':'trace-triangles.csv'}
    for kind,cfg in configs.items():
        state_templates=copy.deepcopy(cfg['Boundaries']['PrescribedPotential'][nold:])
        cfg['Boundaries']['PrescribedPotential']=[{'Index':i,'Attributes':old_sources[0]['Attributes'],'DataFile':str(p)} for i,p in enumerate(traces,start=1)]
        for offset,conductor in enumerate(range(2,conductors+1)):
            source=state_templates[offset];source['Index']=len(active)+offset+1;source['DataFile']=str(lift_paths[conductor]);cfg['Boundaries']['PrescribedPotential'].append(source)
            old_to_new[nold+offset+1]=source['Index']
        cfg['Problem']['Output']=str(output/f'reference-{kind}')
        cfg['Problem']['OutputFormats']={'Paraview':False,'GridFunction':False}
        cfg['Solver']['Electrostatic'].update(Save=0,ResponseMatrix=True,AggregateResponseMatrix=True)
        (output/f'spatial_{kind}.json').write_text(json.dumps(cfg,indent=2)+'\n')
        field='Thin' if kind=='thin' else 'Fabricated'
        model[field+'Matrix']=f'generated/{kind}/domain-response-matrix.csv'
        model[field+'SurfaceMatrix']=f'generated/{kind}/surface-response-matrix.csv'
    updated=copy.deepcopy(library);updated['Models']=[model];updated['TraceCapVersion']=2;updated['MatchingBoxVersion']=2
    (output/'process-library.json').write_text(json.dumps(updated,indent=2)+'\n')
    report={'Version':1,'Model':model_name,'SourceDefinitionChanged':True,'AllRetainedDegreesOfFreedomPreserved':True,
            'OriginalSources':len(old_sources),'Sources':len(configs['thin']['Boundaries']['PrescribedPotential']),
            'OriginalContourDOFs':nold,'ContourDOFs':len(active),'ConductorStates':conductors-1,
            'OldToNewSourceIndices':old_to_new,'FrameFitResidual':residual,'Geometry':geometry_report,
            'OriginalConfigSHA256':{k:sha(c['SourceConfig']) for k,c in cases.items()},
            'OriginalLibrarySHA256':sha(source_library),'Fabrication':fabrication,
            'GeometryInputsSHA256':{name:sha(source_root/name) for name in
                                    ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv')},
            'OriginalSourcesSHA256':{s['DataFile']:sha(s['DataFile']) for s in old_sources},
            'OutputConfigSHA256':{k:sha(output/f'spatial_{k}.json') for k in cases},
            'OutputSourceSHA256':{s['DataFile']:sha(s['DataFile']) for s in configs['thin']['Boundaries']['PrescribedPotential']},
            'ZeroTraceIndices':zero,'LibraryQualified':False}
    (output/'basis-contract.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k in ('Model','OriginalSources','Sources','ContourDOFs','ConductorStates','FrameFitResidual')},indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference_manifest',type=Path)
    parser.add_argument('model')
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    rebuild(args.reference_manifest,args.model,args.output)


if __name__=='__main__':main()
