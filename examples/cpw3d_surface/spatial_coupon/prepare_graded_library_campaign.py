#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare an isolated, mesh-gated complete library campaign from retained provenance."""
import argparse
import copy
import hashlib
import json
import shutil
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--reference-manifest',type=Path,required=True)
    parser.add_argument('--geometry-audit-root',type=Path,required=True)
    parser.add_argument('--python',required=True)
    parser.add_argument('--julia-project',required=True)
    parser.add_argument('--etch-footprints',type=Path,help='JSON map from case key to retained footprint CSV')
    parser.add_argument('--trace-mode',choices=('all','levels'),required=True)
    parser.add_argument('--mesh-statistics-binary',type=Path,required=True)
    args=parser.parse_args()
    footprints=json.loads(args.etch_footprints.read_text()) if args.etch_footprints else {}
    root=args.root.resolve();ref=json.loads(args.reference_manifest.read_text())
    if (root/'campaign.json').exists():parser.error('Refuse to overwrite campaign')
    original_library=root/'reference-library/process-library.json'
    library=json.loads(original_library.read_text());candidate=copy.deepcopy(library)
    candidate['Name']='single-transmon-graded-complete-box'
    candidate['BenchmarkStatus']='PendingGeneration'
    candidate['TraceCapVersion']=2;candidate['MatchingBoxVersion']=2
    new=root/'new-library';new.mkdir()
    records=[]
    for old in ref['Cases']:
        key=old['Key'];kind=old['Kind'];number=key[:2]
        spatial=old['Model'].startswith('spatialedgecluster_')
        input_path=root/'inputs'/number/f'spatial_{kind}.json' if spatial else Path(old['SourceConfig'])
        config=json.loads(input_path.read_text());d=root/'cases'/key;d.mkdir(parents=True)
        raw=json.loads(Path(old['SourceConfig']).read_text())
        if config['Solver']['Order']!=raw['Solver']['Order'] or config['Solver']['Linear']!=raw['Solver']['Linear']:
            raise ValueError('Order/tolerances/linear settings changed')
        if config['Domains']!=raw['Domains'] or config['Boundaries'].get('Postprocessing')!=raw['Boundaries'].get('Postprocessing'):
            raise ValueError('Material or interface settings changed')
        config['Problem']['OutputFormats']={'Paraview':False,'GridFunction':False}
        config['Solver']['Electrostatic'].update(Save=0,ResponseMatrix=True,AggregateResponseMatrix=True)
        mesh=root/'production-meshes'/key/'coupon.msh' if spatial else Path(old['Mesh'])
        config['Model']['Mesh']=str(mesh)
        for phase in ('worker','reducer'):
            p=copy.deepcopy(config);p['Problem']['Output']=str(d/phase);save(d/(phase+'.json'),p)
        source_hashes={s['DataFile']:sha(s['DataFile']) for s in config['Boundaries']['PrescribedPotential']}
        old_model=next(m for m in library['Models'] if m['Name']==old['Model'])
        field='Thin' if kind=='thin' else 'Fabricated'
        record=copy.deepcopy(old)
        record.update(CaseDirectory=str(d),Mesh=str(mesh),RetainedMesh=old['Mesh'],
                      RetainedMeshSHA256=ref['Inputs'][old['Mesh']]['SHA256'],
                      OriginalSources=old['Sources'],Sources=len(config['Boundaries']['PrescribedPotential']),
                      SourceIndices=[s['Index'] for s in config['Boundaries']['PrescribedPotential']],
                      SourceHashes=source_hashes,InputConfig=str(input_path),InputConfigSHA256=sha(input_path),
                      MeshRequired=spatial,InputDirectory=str(root/'inputs'/number),
                      ReferenceDomain=str(root/'reference-library'/old_model[field+'Matrix']),
                      ReferenceSurface=str(root/'reference-library'/old_model[field+'SurfaceMatrix']),
                      DestinationDomain=str(new/old_model[field+'Matrix']),
                      DestinationSurface=str(new/old_model[field+'SurfaceMatrix']),
                      SurfaceSize=.002 if kind=='thin' else .0005,TraceSize=.01,
                      SurfaceAlgorithm=5,TraceConstraintMode=args.trace_mode,
                      MeshStatisticsBinary=str(args.mesh_statistics_binary.resolve()),
                      MeshStatisticsBinarySHA256=sha(args.mesh_statistics_binary),
                      WorkerSHA256=sha(d/'worker.json'),ReducerSHA256=sha(d/'reducer.json'))
        if spatial:
            metrics=args.geometry_audit_root/'audit'/f'{key}-reference-measures.json'
            dest=root/'audit'/f'{key}-reference-measures.json';shutil.copy2(metrics,dest)
            record['ReferenceMeasures']=str(dest)
            if key in footprints:
                footprint=Path(footprints[key])
                destination=root/'inputs'/number/'retained-etch.csv';shutil.copy2(footprint,destination)
                record.update(EtchBoundary=str(destination),EtchBoundarySHA256=sha(destination))
            else:record['EtchBoundary']=None
        records.append(record)
        if kind!='thin':continue
        dest_model=next(m for m in candidate['Models'] if m['Name']==old['Model'])
        folder=Path(old_model['BasisPoints']).parent
        if spatial:
            corrected=json.loads((root/'inputs'/number/'process-library.json').read_text())['Models'][0]
            for name in ('ContourGroups','OpenContourPaths','ZeroTraceIndices'):
                dest_model.pop(name,None)
                if name in corrected:dest_model[name]=corrected[name]
            dest_model['BasisPoints']=str(folder/'basis-points.csv')
            dest_model['TraceMesh']={'Vertices':str(folder/'trace-vertices.csv'),'Triangles':str(folder/'trace-triangles.csv')}
            for name in ('basis-points.csv','trace-vertices.csv','trace-triangles.csv'):
                path=new/folder/name;path.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/'inputs'/number/name,path)
        else:
            for name in [old_model['BasisPoints'],*old_model.get('TraceMesh',{}).values()]:
                dest=new/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/'reference-library'/name,dest)
    binary=root/'bin/palace-reference'
    shutil.copy2(ref['BuildBinary'],binary)
    if sha(binary)!=ref['BinarySHA256']:raise ValueError('Reference executable changed')
    save(root/'candidate-library.json',candidate)
    save(root/'campaign.json',{'Version':1,'Root':str(root),'ReferenceManifest':str(args.reference_manifest.resolve()),
                              'ReferenceManifestSHA256':sha(args.reference_manifest),'BuildBinary':str(binary),
                              'BinarySHA256':sha(binary),'Python':args.python,'JuliaProject':args.julia_project,
                              'Cases':records,'CandidateLibrary':str(root/'candidate-library.json'),
                              'OriginalSourceCount':sum(c['Sources'] for c in ref['Cases']),
                              'SourceCount':sum(c['Sources'] for c in records),'LibraryQualified':False,
                              'Scope':'Corrected complete box basis and graded spatial meshes; historical timing is not a pure mesh-only baseline.'})
    jobs=root/'jobs';jobs.mkdir(exist_ok=True)
    for c in records:
        for mesh_job in ([True,False] if c['MeshRequired'] else [False]):
            tag=('mesh-' if mesh_job else 'response-')+c['Key']
            instance='c8g.48xlarge' if mesh_job else 'r8g.48xlarge'
            nodes=1 if mesh_job else c['Nodes']
            script='run_graded_library_mesh.py' if mesh_job else 'run_graded_library_case.py'
            log=root/'cases'/c['Key']/(('mesh-' if mesh_job else '')+'pbs.log')
            text=f'''#!/bin/bash
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#PBS -N graded-{tag}
#PBS -q normal-g
#PBS -P DS-EM-FEM
#PBS -l select={nodes}:ncpus=192:mpiprocs=192
#PBS -l place=scatter:excl
#PBS -l instance_type={instance}
#PBS -l efa_support=True,subnet_id=subnet-0c98d793bbcebb39a
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o {log}
set -euo pipefail
source /etc/profile.d/modules.sh
module load gcc openmpi arm/armpl julia/1.10.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
{args.python} {root}/tools/{script} {root} {c['Key']}
'''
            (jobs/(tag+'.pbs')).write_text(text)
    print(json.dumps({'Cases':len(records),'Sources':sum(c['Sources'] for c in records),'Root':str(root)},indent=2))


if __name__=='__main__':main()
