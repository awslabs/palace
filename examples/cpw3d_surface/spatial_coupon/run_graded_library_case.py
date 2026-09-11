#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run a complete retained-dimensional response case; never promote accuracy implicitly."""
import csv
import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from archive_storage import reserve, seal


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(path,value):
    temp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def read_rows(path):
    with open(path,newline='') as stream:
        rows=[{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(stream,skipinitialspace=True)]
    if not rows or not all(math.isfinite(v) for row in rows for v in row.values()):
        raise ValueError(f'Empty/nonfinite response table {path}')
    return rows


def matrix_contract(path,kind):
    rows=read_rows(path);columns=set(rows[0])
    required={'basis_i','basis_j','Q_ij (J)'} if kind=='domain' else {'interface','edge','basis_i','basis_j','Q_total_ij (J)'}
    if not required<=columns:raise ValueError(f'Missing mandatory {kind} matrix columns: {required-columns}')
    group_names=[key for key in ('interface','edge','R (m)') if key in columns]
    return {'Columns':sorted(columns),'GroupNames':group_names,
            'Groups':sorted({tuple(row[k] for k in group_names) for row in rows})}


def validate_matrix(path,indices,kind='domain',expected_contract=None):
    contract=matrix_contract(path,kind)
    if expected_contract is not None and contract!=expected_contract:
        raise ValueError('Matrix schema or interface/edge/radius coverage differs from the explicit contract')
    rows=read_rows(path);lookup={i:k for k,i in enumerate(indices)};groups=defaultdict(list)
    key_names=contract['GroupNames']
    for row in rows:groups[tuple(row[key] for key in key_names)].append(row)
    expected={(a,b) for pos,a in enumerate(indices) for b in indices[pos:]}
    output={}
    for key,entries in groups.items():
        seen=set();matrices={name:np.zeros((len(indices),len(indices))) for name in rows[0] if name not in (*key_names,'basis_i','basis_j')}
        for row in entries:
            a,b=int(row['basis_i']),int(row['basis_j'])
            if a!=row['basis_i'] or b!=row['basis_j'] or (a,b) not in expected or (a,b) in seen:
                raise ValueError('Invalid/duplicate response-matrix basis indices')
            seen.add((a,b));i,j=lookup[a],lookup[b]
            for name,matrix in matrices.items():matrix[i,j]=matrix[j,i]=row[name]
        if seen!=expected:raise ValueError('Incomplete full response matrix')
        summary={}
        for name,matrix in matrices.items():
            eigen=np.linalg.eigvalsh(matrix);scale=max(np.max(np.abs(eigen)),1e-300)
            if name in ('Q_ij (J)','Q_total_ij (J)') and eigen[0]<-1e-8*scale:
                raise ValueError(f'Nonpositive response energy matrix {key} {name}: {eigen[0]/scale}')
            summary[name]={'Frobenius':float(np.linalg.norm(matrix)), 'MinimumEigenvalueRelative':float(eigen[0]/scale)}
        output[str(key)]=summary
    return {'Rows':len(rows),'BasisSize':len(indices),'Groups':output}


def compare_unchanged_control(actual,reference):
    a,b=read_rows(actual),read_rows(reference)
    if set(a[0])!=set(b[0]):raise ValueError('Unchanged control matrix schema differs')
    keys=[key for key in ('interface','edge','R (m)','basis_i','basis_j') if key in a[0]]
    left={tuple(row[k] for k in keys):row for row in a}
    right={tuple(row[k] for k in keys):row for row in b}
    if len(left)!=len(a) or len(right)!=len(b) or left.keys()!=right.keys():
        raise ValueError('Unchanged control matrix keys differ')
    reports={}
    for column in a[0]:
        if column in keys:continue
        squared=sum((row[column]-right[key][column])**2 for key,row in left.items())
        norm=sum(row[column]**2 for row in right.values())
        error=math.sqrt(squared/max(norm,1e-300))
        if error>1e-5:raise ValueError(f'Unchanged control differs in {column}: {error}')
        reports[column]=error
    return reports


def compact_corner(source,reference,destination):
    rows=read_rows(source);old=read_rows(reference)
    columns=['interface','edge','basis_i','basis_j','Q_total_ij (J)']
    if set(old[0])!=set(columns):return source
    if {r['edge'] for r in rows}!={1.} or len({r['R (m)'] for r in rows})!=1:
        raise ValueError('Unexpected corner localization convention')
    with open(destination,'w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,lineterminator='\n');writer.writeheader()
        for row in rows:
            record={name:row[name] for name in columns};record['Q_total_ij (J)']=row['Q_ij (J)'];writer.writerow(record)
    return destination


def main():
    root=Path(sys.argv[1]).resolve();key=sys.argv[2]
    manifest=json.loads((root/'campaign.json').read_text());case=next(c for c in manifest['Cases'] if c['Key']==key)
    directory=Path(case['CaseDirectory']);status_path=directory/'status.json'
    if status_path.exists():raise RuntimeError('Case already attempted; preserve partial artifacts and recover explicitly')
    status={'Case':key,'State':'Running','StartUTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            'PBSJob':os.environ.get('PBS_JOBID'),'Stages':[],'LibraryQualified':False}
    save(status_path,status);started=time.monotonic()
    def stage(name,environment):
        env=dict(os.environ)
        for name_ in ('PALACE_RESPONSE_ARCHIVE_DIR','PALACE_RESPONSE_ARCHIVE_ONLY','PALACE_RESPONSE_REDUCE_ONLY',
                      'PALACE_RESPONSE_RECYCLE_INITIAL_GUESS','PALACE_RESPONSE_BLOCK_SIZE'):env.pop(name_,None)
        env.update(environment)
        export=[part for name_ in environment for part in ('-x',name_)]
        command=['/usr/bin/time','-v','-o',str(directory/(name+'.time')),'mpirun',*export,
                 '-n',str(case['Ranks']),'--hostfile',os.environ['PBS_NODEFILE'],manifest['BuildBinary'],str(directory/(name+'.json'))]
        before=time.monotonic()
        with (directory/(name+'.log')).open('w') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,env=env)
        text=(directory/(name+'.log')).read_text(errors='replace')
        record={'Stage':name,'WallSeconds':time.monotonic()-before,'ReturnCode':result.returncode,'Command':command,
                'Iterations':[int(v) for v in re.findall(r'PCG solver converged in (\d+) iterations',text)],
                'DOFs':re.findall(r'H1 \(p = \d+\): (\d+), ND \(p = \d+\): (\d+), RT \(p = \d+\): (\d+)',text)}
        status['Stages'].append(record);save(status_path,status)
        if result.returncode or 'did NOT converge' in text:raise RuntimeError(f'{name} failed or contains unconverged sources')
        if name=='worker':
            actual=[int(v) for v in re.findall(r'It \d+/\d+: Index = (\d+)',text)]
            if actual!=case['SourceIndices']:raise RuntimeError('Not every requested source was processed')
    try:
        if not os.environ.get('PBS_JOBID'):raise ValueError('Response stages require a bounded PBS allocation')
        for path,digest in [(manifest['BuildBinary'],manifest['BinarySHA256']),
                            (str(directory/'worker.json'),case['WorkerSHA256']),
                            (str(directory/'reducer.json'),case['ReducerSHA256']),*case['SourceHashes'].items()]:
            if sha(path)!=digest:raise RuntimeError(f'Input changed: {path}')
        configuration=json.loads((directory/'worker.json').read_text())
        if any(item.get('FluxRecovery',False) for item in configuration['Boundaries']['Postprocessing']['Dielectric']):
            raise ValueError('This campaign archive budget supports potential-only response archives')
        mesh_sha=sha(case['Mesh'])
        if case['MeshRequired']:
            mesh=json.loads((Path(case['Mesh']).parent/'mesh-state.json').read_text())
            from run_graded_library_mesh import mesh_recipe
            if (mesh['State']!='Completed' or mesh['MeshSHA256']!=mesh_sha or
                    mesh.get('Recipe')!=mesh_recipe(root,case)):
                raise RuntimeError('Mesh has not passed its current recipe gates')
            stats_path=Path(case['Mesh']).parent/'processed-statistics/mesh-statistics.json'
            if (sha(stats_path)!=mesh['ProcessedStatisticsSHA256'] or
                    sha(case['MeshStatisticsBinary'])!=case['MeshStatisticsBinarySHA256']):
                raise ValueError('Processed-mesh statistics provenance changed')
            processed=json.loads(stats_path.read_text())
            if processed['Order']!=case['Order'] or processed['Mesh']!=case['Mesh']:
                raise ValueError('Processed-mesh statistics belong to a different case')
            ndofs=processed['H1TrueDOFs']
            if ndofs!=mesh['ProcessedH1TrueDOFs'] or ndofs<=0:raise ValueError('Invalid processed H1 count')
            estimate=int(1.1*(8*ndofs*case['Sources']+48*case['Ranks']*case['Sources']))+2**30
            status['ExpectedH1DOFs']=ndofs
        else:
            if mesh_sha!=case['RetainedMeshSHA256']:raise ValueError('Retained control mesh changed')
            estimate=case['EstimatedArchiveBytes']
        reserve(root,key,estimate)
        status['ArchiveReservationBytes']=estimate
        status['MeshSHA256']=mesh_sha
        archive=directory/'archive';archive.mkdir()
        stage('worker',{'PALACE_RESPONSE_ARCHIVE_DIR':str(archive),'PALACE_RESPONSE_ARCHIVE_ONLY':'1'})
        if case['MeshRequired']:
            dofs=status['Stages'][-1]['DOFs']
            if not dofs or int(dofs[-1][0])!=status['ExpectedH1DOFs']:
                raise ValueError('Solver H1 size differs from its post-preprocessing statistics probe')
        for index in case['SourceIndices']:
            for rank in range(case['Ranks']):
                if not (archive/f'source-{index:06d}-rank-{rank:06d}-V.bin').is_file():
                    raise ValueError('Missing archived potential source/rank')
        status['ArchiveBytes']=sum(p.stat().st_size for p in archive.iterdir())
        seal(root,key,status['ArchiveBytes']);status['ArchiveSealed']=True;save(status_path,status)
        stage('reducer',{'PALACE_RESPONSE_ARCHIVE_DIR':str(archive),'PALACE_RESPONSE_REDUCE_ONLY':'1',
                         'PALACE_RESPONSE_BLOCK_SIZE':str(case['BlockSize'])})
        domain=directory/'reducer/domain-response-matrix.csv'
        surface=directory/'reducer/surface-response-matrix.csv'
        surface=compact_corner(surface,case['ReferenceSurface'],directory/'reducer/surface-response-matrix-compact.csv')
        domain_contract=matrix_contract(case['ReferenceDomain'],'domain')
        surface_contract=matrix_contract(case['ReferenceSurface'],'surface')
        interface_position=surface_contract['GroupNames'].index('interface')
        expected_interfaces={item['Index'] for item in configuration['Boundaries']['Postprocessing']['Dielectric']}
        if {g[interface_position] for g in surface_contract['Groups']}!=expected_interfaces:
            raise ValueError('Reference matrix interface contract is incomplete')
        checks={'Domain':validate_matrix(domain,case['SourceIndices'],'domain',domain_contract),
                'Surface':validate_matrix(surface,case['SourceIndices'],'surface',surface_contract)}
        if not case['MeshRequired']:
            checks['UnchangedControl']={
                'Domain':compare_unchanged_control(domain,case['ReferenceDomain']),
                'Surface':compare_unchanged_control(surface,case['ReferenceSurface'])}
        save(directory/'matrix-validation.json',checks)
        for src,dst in ((domain,case['DestinationDomain']),(surface,case['DestinationSurface'])):
            dest=Path(dst);dest.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists():raise ValueError('Refuse to replace an existing candidate matrix')
            shutil.copy2(src,dest)
        status.update(State='Completed',ElapsedSeconds=time.monotonic()-started,FullBasisComplete=True,
                      DomainSHA256=sha(case['DestinationDomain']),SurfaceSHA256=sha(case['DestinationSurface']),
                      AccuracyQualification='Pending full-library and device comparisons')
        save(status_path,status)
        # Successful source archives are retained until the full library is validated.
        with (root/'publish.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            complete=all((Path(c['CaseDirectory'])/'status.json').exists() and
                         json.loads((Path(c['CaseDirectory'])/'status.json').read_text())['State']=='Completed'
                         for c in manifest['Cases'])
            if complete:
                candidate=json.loads((root/'candidate-library.json').read_text())
                candidate['BenchmarkStatus']='GeneratedCandidateNotAccuracyQualified'
                save(root/'new-library/process-library.json',candidate)
    except Exception as error:
        status.update(State='Failed',Error=str(error),ElapsedSeconds=time.monotonic()-started)
        save(status_path,status);raise
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
