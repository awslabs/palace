#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Generate, label and audit one isolated graded coupon before any response solve."""
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def physical_key(dim,attribute):
    if dim==3 or attribute==1:return dim,attribute
    family=attribute//1000
    if family==3:return dim,3100 if attribute>=3100 else 3000
    if family not in (4,5,6):raise ValueError('Unexpected physical family')
    return dim,1000*family+attribute%100


def grouped(entries):
    result={}
    for dim,attr,value in entries:
        key=physical_key(dim,attr);result[key]=result.get(key,0.)+value
    return result


def reference_entries(data):
    return [(3,int(k),v) for k,v in data['MaterialVolumes'].items()]+[(2,int(k),v) for k,v in data['BoundaryAreas'].items()]


def compare_geometry(reference,actual):
    a,b=grouped(reference),grouped(actual)
    if a.keys()!=b.keys():raise ValueError(f'Physical families differ: missing={a.keys()-b.keys()}, extra={b.keys()-a.keys()}')
    worst=max(abs(a[k]-b[k])/max(abs(a[k]),1e-300) for k in a)
    if worst>1e-5:raise ValueError(f'Physical geometry changed by {worst}')
    return worst


def mesh_recipe(root,case):
    inputs=Path(case['InputDirectory']);tools=Path(root)/'tools'
    files=[inputs/name for name in ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv','process.toml')]
    files.append(inputs/'traces/basis-0001.csv')
    if case.get('EtchBoundary'):files.append(Path(case['EtchBoundary']))
    if case.get('MeshStatisticsBinary'):files.append(Path(case['MeshStatisticsBinary']))
    files.extend(tools/name for name in ('run_graded_library_mesh.py','mesh_spatial_coupon.jl','mesh_graded_tet_experiment.jl',
        'frozen_volume_study.jl','graded_curve_distance.jl','graded_size_points.jl','interface_ownership.jl',
        'ownership_bernstein.jl','label_interface_patches.jl','relabel_frozen_interface_mesh.jl'))
    return {'Files':{str(path):sha(path) for path in files},
            'Settings':{name:case.get(name) for name in ('Kind','SurfaceSize','TraceSize','SurfaceAlgorithm','TraceConstraintMode')},
            'Volume':{'MinimumSize':.002,'NearGrowth':.5,'FarGrowth':2.,'TransitionDistance':.03,'MaximumSize':.5,'OptimizeVolume':False},
            'MeshFormat':'2.2 binary'}


def main():
    root=Path(sys.argv[1]).resolve();key=sys.argv[2]
    manifest=json.loads((root/'campaign.json').read_text());case=next(c for c in manifest['Cases'] if c['Key']==key)
    if not case['MeshRequired']:raise ValueError('Not a graded mesh case')
    directory=Path(case['Mesh']).parent;directory.mkdir(parents=True,exist_ok=True)
    recipe=mesh_recipe(root,case)
    status_path=directory/'mesh-state.json'
    if status_path.exists():
        old=json.loads(status_path.read_text())
        if old['State']=='Completed' and old.get('Recipe')==recipe and sha(case['Mesh'])==old['MeshSHA256']:
            print('Previously validated mesh:',case['Mesh']);return
        raise RuntimeError('Partial mesh stages retained; use a separate recovery attempt')
    status={'Case':key,'State':'Running','StartUTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'Stages':[], 'Recipe':recipe}
    save(status_path,status);started=time.monotonic()
    tools=root/'tools';inputs=Path(case['InputDirectory']);env=dict(os.environ)
    env.update(TET_GEOMETRY_ORDER='1',TET_SURFACE_ALGORITHM=str(case['SurfaceAlgorithm']),TET_ALGORITHM3D='10',
               TET_HXT_QUALITY='.1',TET_TRACE_CONSTRAINT_MODE=case.get('TraceConstraintMode','all'),TET_TRACE_SIZE=str(case['TraceSize']),
               TET_CAP_SIZE='0',TET_VERBOSITY='3',JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    julia=['julia','--startup-file=no','--project='+manifest['JuliaProject']]
    base=[*julia,str(tools/'mesh_graded_tet_experiment.jl'),str(inputs),case['Kind']]
    arguments=['1',str(case['SurfaceSize']),'.5','--process',str(inputs/'process.toml'),
               '--matching-trace',str(inputs/'traces/basis-0001.csv')]
    if case.get('EtchBoundary'):
        if sha(case['EtchBoundary'])!=case['EtchBoundarySHA256']:raise ValueError('Etch footprint changed')
        arguments+=['--etch-boundary',case['EtchBoundary']]
    def run(name,command,seconds,memory):
        log=directory/(name+'.log')
        wrapper=[manifest['Python'],str(tools/'run_bounded_mesher.py'),'--seconds',str(seconds),
                 '--memory-gib',str(memory),'--log',str(log),'--',*map(str,command)]
        subprocess.run(wrapper,env=env,check=True)
        record=json.loads(Path(str(log)+'.json').read_text());record['Stage']=name
        status['Stages'].append(record);save(status_path,status)
    try:
        for path,digest in case['SourceHashes'].items():
            if sha(path)!=digest:raise ValueError('A source trace changed')
        reference=json.loads(Path(case['ReferenceMeasures']).read_text())
        run('cad',[*base,str(directory/'cad.msh'),*arguments,'--geometry-only'],300,16)
        cad=directory/'cad.msh.cad-measures.csv'
        entries=[(int(r['dimension']),int(r['attribute']),float(r['measure'])) for r in csv.DictReader(cad.open())]
        status['CADMaximumRelativeGeometryDifference']=compare_geometry(reference_entries(reference),entries)
        save(status_path,status)
        settings=directory/'volume-campaign.toml'
        settings.write_text('MaxElements = 16000000\nMaxNodes = 5000000\nOptimizeVolume = false\n'
                            '[[Variants]]\nName = "selected"\nMinimumSize = 0.002\n'
                            'NearGrowth = 0.5\nFarGrowth = 2.0\nTransitionDistance = 0.03\nMaximumSize = 0.5\n')
        frozen=directory/'family-selected.msh'
        if frozen.exists():
            raise ValueError('Unbound prebuilt pilot mesh: preserve it and use a fresh campaign mesh directory')
        run('volume',[*base,str(directory/'family.msh'),*arguments,'--volume-study',str(settings),
                      '--reference-measures',str(cad)],18000,300 if os.environ.get('PBS_JOBID') else 110)
        log=(directory/'volume.log').read_text()
        if case.get('TraceConstraintMode','all')=='all':
            skipped=re.findall(r'off-box segments=(\d+)',log)
            if not skipped or any(int(x) for x in skipped):raise ValueError('Matching trace does not cover the box correctly')
        # Source geometry is independently complete even when only level curves
        # are constrained in the mesh. PDE accuracy remains a separate gate.
        from box_trace import validate_box_trace
        from audit_trace_continuity import read_trace
        triangles=read_trace(inputs/'traces/basis-0001.csv')
        points=sorted({p for tri in triangles.values() for p,_ in tri});lookup={p:i for i,p in enumerate(points)}
        status['TraceGeometry']=validate_box_trace(points,[[lookup[p] for p,_ in tri] for tri in triangles.values()])
        if mesh_recipe(root,case)!=recipe:raise ValueError('Meshing inputs changed during generation')
        expected=directory/'expected-attributes.csv'
        expected.write_text('dimension,attribute,measure\n'+''.join(f'2,{a},{v}\n' for a,v in reference['BoundaryAreas'].items()))
        run('relabel',[*julia,str(tools/'relabel_frozen_interface_mesh.jl'),str(inputs),case['Kind'],str(frozen),case['Mesh'],
                       '--expected-measures',str(expected)],2400,280 if os.environ.get('PBS_JOBID') else 110)
        with open(case['Mesh'],'rb') as f:
            if f.readline().strip()!=b'$MeshFormat' or f.readline().strip()!=b'2.2 1 8':raise ValueError('Required MSH 2.2 binary output missing')
        run('measures',[root/'bin/audit_mesh_measures',case['Mesh'],directory/'measures.json'],2400,32)
        actual=json.loads((directory/'measures.json').read_text())
        status['MeshMaximumRelativeGeometryDifference']=compare_geometry(reference_entries(reference),reference_entries(actual))
        if set(actual['BoundaryAreas'])!=set(reference['BoundaryAreas']):raise ValueError('Missing/extra slot or conductor attribute')
        run('quality',[root/'bin/audit_coupon_mesh',case['Mesh'],inputs/'plan-view-boundary.csv',
                       '1' if case['Kind']=='fabricated' else '0',directory/'quality.json'],1200,32)
        quality=json.loads((directory/'quality.json').read_text())
        if quality['NonpositiveCenterJacobians'] or quality['GeometryCounts']!={'Tetrahedron':quality['Elements']}:
            raise ValueError('Invalid or non-tetrahedral mesh')
        stats_binary=Path(case['MeshStatisticsBinary'])
        if sha(stats_binary)!=case['MeshStatisticsBinarySHA256']:raise ValueError('Mesh-statistics executable changed')
        # Count H1 after the real Palace loading, cracking and initial-refinement
        # path. Raw serialized topology is NOT a bound on the processed mesh.
        probe=json.loads((Path(case['CaseDirectory'])/'worker.json').read_text())
        probe['Problem']['Output']=str(directory/'processed-statistics')
        probe_config=directory/'mesh-statistics-input.json';save(probe_config,probe)
        run('processed-statistics',['mpirun','-n','1',stats_binary,'--mesh-statistics',probe_config],2400,100)
        stats_path=directory/'processed-statistics/mesh-statistics.json'
        processed=json.loads(stats_path.read_text())
        if processed['Order']!=case['Order'] or processed['H1TrueDOFs']<=0:
            raise ValueError('Invalid processed H1 statistics')
        for name in ('CrackInternalBoundaryElements','RefineCrackElements'):
            if processed[name]!=probe['Model'].get(name,True):raise ValueError('Mesh preprocessing policy differs')
        status.update(ProcessedH1TrueDOFs=processed['H1TrueDOFs'],
                      ProcessedStatisticsSHA256=sha(stats_path),MeshStatisticsBinarySHA256=sha(stats_binary))
        status.update(State='Completed',MeshSHA256=sha(case['Mesh']),Elements=quality['Elements'],
                      MaximumKappa=quality['KappaPercentiles']['100.000000'],
                      GeometryAndCoverageVerified=True,SlotPartitionQualified=False,
                      ElapsedSeconds=time.monotonic()-started)
    except Exception as error:
        status.update(State='Failed',Error=str(error),ElapsedSeconds=time.monotonic()-started)
        save(status_path,status);raise
    save(status_path,status)
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
