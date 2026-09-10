#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Execute every mesh pilot with one executable, host, MPI count and p5/tolerance.
Run from a dedicated host/allocation. No scheduling or resource auto-selection.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda:stream.read(8*1024*1024),b''):h.update(data)
    return h.hexdigest()


def timing(log):
    text=log.read_text()
    match=re.search(r'H1 \(p = 5\):\s*(\d+)',text)
    sources=[]
    for m in re.finditer(r'Response source timing: index=(\d+), iterations=(\d+), solve_seconds=([^,]+), total_seconds=([^\s]+)',text):
        sources.append(dict(Index=int(m[1]),Iterations=int(m[2]),SolveSeconds=float(m[3]),TotalSeconds=float(m[4])))
    table={};active=False;category=''
    for line in text.splitlines():
        if line.startswith('Elapsed Time Report'):active=True;continue
        if active and line.startswith('Peak Memory'):break
        if active:
            m=re.match(r'^(\s*)([A-Za-z][A-Za-z0-9 /_-]*?)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*$',line)
            if m:
                name=m[2].strip()
                if not m[1]:category=name;key=name
                else:key=category+'/'+name
                table[key]={'Min':float(m[3]),'Max':float(m[4]),'Mean':float(m[5])}
    return dict(H1Dofs=int(match[1]) if match else None,Sources=sources,ElapsedTimerSeconds=table)


def matrix(path,interfaces=False):
    result={}
    with path.open() as stream:
        for raw in csv.DictReader(stream):
            row={k.strip():float(v) for k,v in raw.items()}
            group=int(row['interface']) if interfaces else 0
            key=(int(row['basis_i']),int(row['basis_j']))
            value=row['Q_total_ij (J)' if interfaces else 'Q_ij (J)']
            if not math.isfinite(value):raise ValueError('Non-finite matrix entry')
            entries=result.setdefault(group,{})
            if key in entries and entries[key]!=value:raise ValueError('Conflicting duplicate matrix entry')
            entries[key]=value
    # Palace writes compact symmetric matrices (upper triangle). Reconstruct the
    # lower triangle before Frobenius norms, defects, or spectral comparisons.
    for entries in result.values():
        for (i,j),value in list(entries.items()):
            reverse=(j,i)
            if reverse in entries and not math.isclose(entries[reverse],value,rel_tol=1e-10,abs_tol=0.):
                raise ValueError('Nonsymmetric response matrix')
            entries.setdefault(reverse,value)
    return result


def norm(matrix):return math.sqrt(sum(v*v for v in matrix.values()))

def difference(a,b):
    if a.keys()!=b.keys():raise ValueError('Matrix basis indices differ')
    delta={k:b[k]-v for k,v in a.items()}
    return norm(delta)/norm(a) if norm(a) else None


def diagonal_differences(a,b):
    if a.keys()!=b.keys():raise ValueError('Matrix basis indices differ')
    return {str(i):abs(b[(i,i)]-v)/abs(v) if v else (0. if b[(i,i)]==0 else None)
            for (i,j),v in a.items() if i==j}


def compare(root,manifest,records):
    result={}
    for case in manifest['Cases']:
        name=case['Name'];kind=case['Kind']
        if name not in records or not records[name].get('Completed'):continue
        baseline=f'{kind}-prism'
        if baseline not in records or not records[baseline].get('Completed'):continue
        actual=matrix(root/'postpro'/f'{name}-reducer/domain-response-matrix.csv')[0]
        reference=matrix(root/'postpro'/f'{baseline}-reducer/domain-response-matrix.csv')[0]
        entry={'DomainMatrixRelativeDifference':difference(reference,actual),
               'DomainDiagonalRelativeDifferences':diagonal_differences(reference,actual),
               'WorkerPlusReducerSpeedup':records[baseline]['TotalSeconds']/records[name]['TotalSeconds'],
               'H1DofReduction':records[baseline]['Worker']['H1Dofs']/records[name]['Worker']['H1Dofs']}
        if kind=='fabricated':
            a=matrix(root/'postpro'/f'{baseline}-reducer/surface-response-matrix.csv',True)
            b=matrix(root/'postpro'/f'{name}-reducer/surface-response-matrix.csv',True)
            entry['InterfaceMatrixRelativeDifferences']={str(i):difference(a[i],b[i]) for i in a}
            entry['InterfaceDiagonalRelativeDifferences']={str(i):diagonal_differences(a[i],b[i]) for i in a}
        result[name]=entry
    defects={}
    strategies=sorted({c['Strategy'] for c in manifest['Cases']})
    matrices={}
    for strategy in strategies:
        names=[f'{kind}-{strategy}' for kind in ('thin','fabricated')]
        if all(n in records and records[n].get('Completed') for n in names):
            t,f=[matrix(root/'postpro'/f'{n}-reducer/domain-response-matrix.csv')[0] for n in names]
            matrices[strategy]=({k:f[k]-v for k,v in t.items()},t)
    if 'prism' in matrices:
        reference,thin=matrices['prism']
        for strategy,(value,_) in matrices.items():
            delta={k:value[k]-v for k,v in reference.items()}
            defects[strategy]={'RelativeToDefect':difference(reference,value),
                               'RelativeToThinDomain':norm(delta)/norm(thin)}
    return {'PerMesh':result,'DomainDefect':defects,'ThinRawSPRGated':False,'LibraryQualified':False}


def run(args):
    root=args.root.resolve();binary=args.palace.resolve();os.chdir(root)
    manifest=json.loads(Path('manifest.json').read_text())
    for name,expected in manifest['InputSHA256'].items():
        if digest(Path(name))!=expected:raise ValueError(f'Changed benchmark input {name}')
    binary_hash=digest(binary)
    summary={'Version':1,'Host':platform.node(),'Machine':platform.machine(),'Python':platform.python_version(),
             'Executable':str(binary),'ExecutableSHA256':binary_hash,'Ranks':args.ranks,'Binding':args.binding,
             'Order':manifest['Order'],'SolverTolerance':manifest['SolverTolerance'],
             'Limits':{'WorkerSeconds':args.worker_seconds,'ReducerSeconds':args.reducer_seconds,'MemoryGiB':args.memory_gib},
             'Cases':{},'Scope':manifest['Scope'],'LibraryQualified':False}
    try:summary['CPUInformation']=subprocess.run(['lscpu'],capture_output=True,text=True,timeout=10).stdout
    except FileNotFoundError:summary['CPUInformation']='unavailable'
    output=Path(args.results)
    if output.exists():raise ValueError('Refuse to overwrite benchmark results')
    def save():
        temporary=output.with_suffix('.tmp');temporary.write_text(json.dumps(summary,indent=2)+'\n');temporary.replace(output)
    base_env={k:v for k,v in os.environ.items() if not k.startswith('PALACE_RESPONSE_')}
    for case in manifest['Cases']:
        name=case['Name']
        if args.case and name not in args.case:continue
        record={'Completed':False};summary['Cases'][name]=record;save()
        try:
            if digest(binary)!=binary_hash or digest(Path(case['Mesh']))!=case['MeshSHA256']:
                raise ValueError('Executable or mesh changed during benchmark')
            archive=root/'archives'/name;archive.mkdir(parents=True,exist_ok=False)
            for stage in ('worker','reducer'):
                log=root/'logs'/f'{name}-{stage}.log'
                env=dict(base_env,PALACE_RESPONSE_ARCHIVE_DIR=str(archive))
                if stage=='worker':env.update(PALACE_RESPONSE_ARCHIVE_ONLY='1',PALACE_RESPONSE_SOURCE_TIMING='1')
                else:env.update(PALACE_RESPONSE_REDUCE_ONLY='1',PALACE_RESPONSE_BLOCK_SIZE=str(len(case['Sources'])))
                exports=[arg for key in env if key.startswith('PALACE_RESPONSE_') for arg in ('-x',key)]
                seconds=args.worker_seconds if stage=='worker' else args.reducer_seconds
                command=['mpirun','--bind-to',args.binding,'-n',str(args.ranks),*exports,str(binary),f'configs/{name}-{stage}.json']
                subprocess.run([sys.executable,'tools/run_bounded_mesher.py','--seconds',str(seconds),
                                '--memory-gib',str(args.memory_gib),'--log',str(log),'--',*command],
                               env=env,check=True,timeout=seconds+20)
                resource=json.loads(Path(str(log)+'.json').read_text())
                record[stage.title()]={**resource,**timing(log)};save()
            if len(record['Worker']['Sources'])!=len(case['Sources']):
                raise ValueError('Missing per-source timing evidence')
            record['Completed']=True
            record['TotalSeconds']=record['Worker']['Seconds']+record['Reducer']['Seconds']
            record['PeakProcessTreeRSSBytes']=max(record['Worker']['PeakProcessTreeRSSBytes'],record['Reducer']['PeakProcessTreeRSSBytes'])
            record['NonSourceWorkerWallSeconds']=record['Worker']['Seconds']-sum(s['TotalSeconds'] for s in record['Worker']['Sources'])
        except (ValueError,subprocess.SubprocessError) as error:
            record['Error']=str(error)
            save()
            # A failed solve is not a numerical comparison. Preserve archives and
            # stop the sequential benchmark rather than burning the next case.
            return False
        summary['Comparison']=compare(root,manifest,summary['Cases']);save()
        print(json.dumps({'Case':name,'Dofs':record['Worker']['H1Dofs'],'Seconds':record['TotalSeconds'],
                          'Iterations':[s['Iterations'] for s in record['Worker']['Sources']]}),flush=True)
    return True


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--palace',type=Path,required=True)
    p.add_argument('--ranks',type=int,required=True)
    p.add_argument('--binding',choices=('core','none'),default='core')
    p.add_argument('--worker-seconds',type=float,default=480)
    p.add_argument('--reducer-seconds',type=float,default=240)
    p.add_argument('--memory-gib',type=float,default=100)
    p.add_argument('--case',action='append')
    p.add_argument('--results',default='timings.json',help='Separate result ledger for selected follow-up cases')
    args=p.parse_args()
    if args.ranks<1 or not all(math.isfinite(v) and v>0 for v in (args.worker_seconds,args.reducer_seconds,args.memory_gib)):
        p.error('positive resource limits are required')
    raise SystemExit(0 if run(args) else 1)
