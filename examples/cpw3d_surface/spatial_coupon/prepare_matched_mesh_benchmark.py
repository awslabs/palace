#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare a portable, matched p5 pilot from prism controls and frozen-surface tet variants.

The complete pilot basis consists of every conductor state and two affine matching
traces. This is not a truncated production library: it is an explicitly small test
basis used to compare meshing cost and accuracy under identical solve settings.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tomllib

HERE=Path(__file__).resolve().parent


def family(attribute):
    return attribute//1000 if attribute!=1 else 1


def grouped_areas(areas):
    result={}
    for key,value in areas.items():
        a=int(key);f=family(a)
        key=(f,a%100 if f in (4,5,6) else 0)
        result[str(key)]=result.get(str(key),0.)+value
    return result


def prepare(args):
    args.output.mkdir(parents=True,exist_ok=False)
    for d in ('meshes','configs','traces','input','audit','tools'):(args.output/d).mkdir()
    for name in ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv','process.toml'):
        shutil.copy2(args.input/name,args.output/'input'/name)
    cases=[];bounds=None;targets={}
    for kind in ('thin','fabricated'):
        study_path=args.mesh_root/f'{args.prefix}-{kind}.msh.volume-study.toml'
        study=tomllib.load(study_path.open('rb'))
        for suffix in ('.volume-study.toml','.sizing.txt','.process.toml','.geometry-gate.txt'):
            source=args.mesh_root/f'{args.prefix}-{kind}.msh{suffix}'
            if source.exists():shutil.copy2(source,args.output/'audit'/source.name)
        if bounds is None:bounds=(study['Lower'],study['Upper'])
        if bounds!=(study['Lower'],study['Upper']):raise ValueError('Thin/fabricated outer domains differ')
        sources=[('prism',args.mesh_root/f'{args.prefix}-prism-{kind}{"-full" if kind=="fabricated" else ""}.msh',None)]
        sources += [(v['Name'],Path(v['Mesh']),dict(v,SharedSurfaceSeconds=study['SurfaceSeconds'],
                                                   IndependentMeshingSecondsEstimate=study['SurfaceSeconds']+v['VolumeTotalSeconds']))
                    for v in study['Variants']]
        for strategy,mesh,metadata in sources:
            name=f'{kind}-{strategy}';dest=args.output/'meshes'/f'{name}.msh'
            shutil.copy2(mesh,dest)
            measure=args.output/'audit'/f'{name}-measures.json'
            quality=args.output/'audit'/f'{name}-quality.json'
            subprocess.run([str(args.measures_bin),str(dest),str(measure)],check=True,timeout=40)
            subprocess.run([str(args.quality_bin),str(dest),str(args.output/'input/plan-view-boundary.csv'),str(int(kind=='fabricated')),str(quality)],check=True,timeout=40)
            measures=json.loads(measure.read_text());qualities=json.loads(quality.read_text())
            attributes=sorted(int(a) for a in measures['BoundaryAreas'])
            geometry={'Volumes':measures['MaterialVolumes'],'Areas':grouped_areas(measures['BoundaryAreas'])}
            if kind not in targets:targets[kind]=geometry
            error=0.
            for group in ('Volumes','Areas'):
                if geometry[group].keys()!=targets[kind][group].keys():raise ValueError(f'{name}: physical families differ')
                error=max(error,max(abs(v/targets[kind][group][k]-1) for k,v in geometry[group].items()))
            if error>1e-6:raise ValueError(f'{name}: prism/tet geometry differs by {error}')
            if metadata is None:
                prism_metadata=Path(str(mesh.resolve())+'.metadata.json')
                if prism_metadata.exists():metadata=json.loads(prism_metadata.read_text())
            cases.append(dict(Name=name,Kind=kind,Strategy=strategy,Mesh=str(dest.relative_to(args.output)),
                              MeshSHA256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                              Elements=qualities['Elements'],Kappa=qualities['KappaPercentiles'],
                              GeometryRelativeDifference=error,Attributes=attributes,
                              FrozenSurfaceSHA256=metadata.get('BoundarySHA256') if metadata else None,
                              MeshGeneration=metadata))
    lower,upper=bounds
    vertices=[(x,y,z) for z in (lower[2],upper[2]) for y in (lower[1],upper[1]) for x in (lower[0],upper[0])]
    faces=[(0,1,3,2),(4,5,7,6),(0,1,5,4),(2,3,7,6),(0,2,6,4),(1,3,7,5)]
    for label in ('zero','affine-x','affine-y'):
        with (args.output/'traces'/f'{label}.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['x','y','z','V','triangle'])
            index=0
            for a,b,c,d in faces:
                for tri in ((a,b,c),(a,c,d)):
                    index+=1
                    for vertex in tri:
                        point=vertices[vertex]
                        axis=0 if label=='affine-x' else 1
                        value=0. if label=='zero' else (point[axis]-(lower[axis]+upper[axis])/2)/((upper[axis]-lower[axis])/2)
                        writer.writerow([*point,value,index])
    for case in cases:
        attributes=case['Attributes'];kind=case['Kind']
        metals=[a for a in attributes if family(a) in (4,5,6)]
        conductors=sorted({a%100 for a in metals})
        sources=[dict(Index=i,Attributes=[1],DataFile='traces/zero.csv',TerminalAttributes=[a for a in metals if a%100==conductor])
                 for i,conductor in enumerate(conductors,1)]
        labels=[f'conductor-{c}' for c in conductors]
        for label in ('affine-x','affine-y'):
            sources.append(dict(Index=len(sources)+1,Attributes=[1],DataFile=f'traces/{label}.csv'));labels.append(label)
        dielectric=[]
        for index,(name,eps) in enumerate((('MA',10.),('MS',11.47),('SA',4.)),1):
            f={'MA':6 if kind=='fabricated' else 4,'MS':5 if kind=='fabricated' else 4,'SA':3}[name]
            dielectric.append(dict(Index=index,Attributes=[a for a in attributes if family(a)==f],Type=name,
                                   Thickness=.002,Permittivity=eps,LossTan=.001,LocalizeEdgeEnergy=True,
                                   SaveLocalEdgeEnergy=False,EdgeAttributes=[a for a in attributes if family(a)==3],
                                   EdgeDistances=[.5],EdgeFrameNormal=[0.,0.,1.],EdgeExcludeAttributes=[1]))
        config=dict(Problem=dict(Type='Electrostatic',Verbose=1,OutputFormats=dict(Paraview=False,GridFunction=False)),
                    Model=dict(Mesh=case['Mesh'],L0=1e-6,Refinement=dict(MaxIts=0)),
                    Domains=dict(Materials=[dict(Attributes=[1],Permittivity=11.47),dict(Attributes=[2],Permittivity=1.)]),
                    Boundaries=dict(Ground=dict(Attributes=metals),PrescribedPotential=sources,Postprocessing=dict(Dielectric=dielectric)),
                    Solver=dict(Order=5,Electrostatic=dict(Save=0,ResponseMatrix=True,AggregateResponseMatrix=True),
                                Linear=dict(Type='BoomerAMG',KSPType='CG',Tol=1e-8,MaxIts=1000,
                                            EstimatorTol=.5,EstimatorMaxIts=5,EstimatorMG=True)))
        case['Sources']=labels
        for stage in ('worker','reducer'):
            config['Problem']['Output']=f'postpro/{case["Name"]}-{stage}'
            path=args.output/'configs'/f'{case["Name"]}-{stage}.json';path.write_text(json.dumps(config,indent=2)+'\n')
    manifest=dict(Version=1,Scope='Matched pilot: complete conductor states plus two affine traces, not full production trace basis',
                  Order=5,SolverTolerance=1e-8,GeometryOrder=1,InterfaceTypes={1:'MA',2:'MS',3:'SA'},Bounds=dict(Lower=lower,Upper=upper),
                  Notes=['All conductor masks must lie strictly inside the matching box for the zero-trace conductor states.',
                         'Tet variants share an identical surface triangulation per kind. Prism surface mesh is independent.',
                         'Thin raw SPR is not an acceptance metric. Interface validation uses fixed whole-family integrals.'],
                  Cases=cases)
    # No conductor may intersect the outer matching boundary in this pilot: that
    # would make the zero outer trace incompatible with a nonzero conductor state.
    with (args.input/'plan-view-boundary.csv').open() as stream:
        for row in csv.DictReader(stream):
            x,y=float(row['X']),float(row['Y'])
            if not lower[0]<x<upper[0] or not lower[1]<y<upper[1]:
                raise ValueError('Conductor touches matching boundary; this pilot needs a compatible trace lift')
    manifest['InputSHA256']={str(p.relative_to(args.output)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for folder in ('input','traces','configs','audit') for p in (args.output/folder).iterdir()}
    for name in ('run_bounded_mesher.py','run_matched_mesh_benchmark.py'):
        shutil.copy2(HERE/name,args.output/'tools'/name)
    manifest['ToolSHA256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.output/'tools').iterdir()}
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'Cases':len(cases),'SourceCount':len(cases[0]['Sources']),'Output':str(args.output)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--mesh-root',type=Path,required=True)
    p.add_argument('--prefix',default='strip')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--quality-bin',type=Path,required=True)
    p.add_argument('--measures-bin',type=Path,required=True)
    prepare(p.parse_args())
