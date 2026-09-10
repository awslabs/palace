#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Partition quadrature sweep on frozen archived p5 fields; no new PDE solves.

The unmasked interface matrix and domain matrix are controls. Every slot of each
physical family uses the same quadrature, so sample-wise ownership is exhaustive.
Agreement between orders is an empirical diagnostic, not a rigorous error bound.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

HERE=Path(__file__).resolve().parent


def rows(path):
    with path.open() as stream:
        return [{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def matrices(path,interface=False):
    data=rows(path)
    ids=sorted({int(r['basis_i']) for r in data}|{int(r['basis_j']) for r in data})
    index={i:k for k,i in enumerate(ids)}
    result={}
    for r in data:
        key=int(r['interface']) if interface else 0
        result.setdefault(key,np.zeros((len(ids),len(ids))))[index[int(r['basis_i'])],index[int(r['basis_j'])]]=r['Q_total_ij (J)' if interface else 'Q_ij (J)']
    return result


def relative(a,b):
    norm=np.linalg.norm(a)
    return float(np.linalg.norm(a-b)/norm) if norm else (0. if np.linalg.norm(b)==0 else None)


def run(args):
    args.output.mkdir(parents=True,exist_ok=False)
    level=args.case/f'level{args.level}'
    original=json.loads((level/'probe.json').read_text())
    mapping=json.loads((level/'interface-map.json').read_text())
    archive=level/'archives'
    if not archive.exists():raise ValueError('Missing archived fields')
    groups={}
    for info in mapping.values():
        original_attribute=info['CanonicalAttribute']
        family=original_attribute//1000
        slot=original_attribute%100 if family==3 else (original_attribute%1000)//100
        group=groups.setdefault(info['Group'],{'Attributes':set(),'Slots':set(),
                                               'OwnershipGroup':0 if family==3 else original_attribute%100,
                                               'Type':{3:'SA',5:'MS',6:'MA'}[family],
                                               'Permittivity':{3:4.,5:11.47,6:10.}[family]})
        group['Attributes'].add(info['Attribute']);group['Slots'].add(slot)
    tools=args.output/'tools';tools.mkdir()
    names=['mesh_spatial_coupon.jl','export_quadrature_ownership.jl','run_bounded_mesher.py','run_ownership_quadrature.py']
    for name in names:shutil.copy2(HERE/name,tools/name)
    env={k:v for k,v in os.environ.items() if not k.startswith('PALACE_RESPONSE_')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',JULIA_NUM_THREADS='1')
    def bounded(command,log,seconds=None):
        seconds=args.seconds if seconds is None else seconds
        subprocess.run([sys.executable,str(tools/'run_bounded_mesher.py'),'--seconds',str(seconds),
                        '--memory-gib','10','--log',str(log),'--',*command],env=env,check=True,timeout=seconds+15,
                       stdout=subprocess.DEVNULL)
    ownership=args.output/'ownership.csv'
    bounded([args.julia,'--startup-file=no',f'--project={args.julia_project}',
             str(tools/'export_quadrature_ownership.jl'),str(args.case/'input'),str(ownership)],
            args.output/'export.log',25)
    worker=json.loads((level/'worker.log.json').read_text())['Command']
    ranks=int(worker[worker.index('-n')+1])
    env.update(PALACE_RESPONSE_ARCHIVE_DIR=str(archive),PALACE_RESPONSE_REDUCE_ONLY='1',
               PALACE_RESPONSE_BLOCK_SIZE='2')
    reference_domain=matrices(level/'postpro/domain-response-matrix.csv')[0]
    summary={'Scope':'Fixed-FE-field quadrature convergence, not library/PDE qualification',
             'LibraryQualified':False,'LcOverride':args.lc,'Source':str(level),'Ranks':ranks,'Steps':[],
             'PalaceSHA256':hashlib.sha256(args.palace.read_bytes()).hexdigest(),
             'OwnershipSHA256':hashlib.sha256(ownership.read_bytes()).hexdigest()}
    previous=None
    for order in args.orders:
        output=args.output/f'q{order}'
        config=json.loads(json.dumps(original));config['Problem']['Output']=str(output)
        if args.lc is not None:
            config['Model']['Lc']=args.lc
        dielectric=[];indices={};index=0
        for key,group in sorted(groups.items()):
            base=dict(Attributes=sorted(group['Attributes']),Type=group['Type'],Thickness=.002,
                      Permittivity=group['Permittivity'],LossTan=.001,
                      LocalizeEdgeEnergy=True,SaveLocalEdgeEnergy=False,
                      EdgeAttributes=sorted(a for g in groups.values() if g['Type']=='SA' for a in g['Attributes']),
                      EdgeExcludeAttributes=[1],EdgeFrameNormal=[0.,0.,1.],EdgeDistances=[.5])
            index+=1;dielectric.append(dict(base,Index=index));indices[(key,None)]=index
            for slot in sorted(group['Slots']):
                index+=1;indices[(key,slot)]=index
                dielectric.append(dict(base,Index=index,OwnershipDataFile=str(ownership),
                                       OwnershipGroup=group['OwnershipGroup'],OwnershipSlot=slot,
                                       OwnershipQuadratureOrder=order))
        config['Boundaries']['Postprocessing']['Dielectric']=dielectric
        path=args.output/f'q{order}.json';path.write_text(json.dumps(config,indent=2)+'\n')
        bounded(['mpirun','-n',str(ranks),'-x','PALACE_RESPONSE_ARCHIVE_DIR','-x','PALACE_RESPONSE_REDUCE_ONLY',
                 '-x','PALACE_RESPONSE_BLOCK_SIZE',str(args.palace),str(path)],args.output/f'q{order}.log')
        data=matrices(output/'surface-response-matrix.csv',True)
        domain=matrices(output/'domain-response-matrix.csv')[0]
        conservation={};slots={};minimum=float('inf')
        for key,group in groups.items():
            total=data[indices[(key,None)]]
            summed=np.zeros_like(total)
            for slot in group['Slots']:
                m=data[indices[(key,slot)]];summed+=m
                slots[f'{key}/slot-{slot}']=m
                scale=np.linalg.norm(m)
                if scale:
                    minimum=min(minimum,float(np.linalg.eigvalsh((m+m.T)/2).min()/scale))
                else:
                    minimum=min(minimum,0.)
            conservation[key]=relative(total,summed)
        change={key:relative(previous[key],m) for key,m in slots.items()} if previous is not None else None
        step={'QuadratureOrder':order,'DomainMatrixRelativeChange':relative(reference_domain,domain),
              'FamilyConservationRelativeError':conservation,'MinimumNormalizedEigenvalue':minimum,
              'SlotMatrixChange':change,'SlotMatrices':{k:v.tolist() for k,v in slots.items()}}
        summary['Steps'].append(step)
        (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps({k:v for k,v in step.items() if k!='SlotMatrices'}),flush=True)
        if step['DomainMatrixRelativeChange']>1e-10 or max(conservation.values())>1e-9 or minimum< -1e-10:
            raise ValueError('Domain/conservation/positive-energy control failed')
        previous=slots


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True,help='Completed partition-study case directory')
    parser.add_argument('--level',type=int,default=0)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--palace',type=Path,required=True)
    parser.add_argument('--julia',default='julia')
    parser.add_argument('--julia-project',type=Path,required=True)
    parser.add_argument('--orders',type=int,nargs='+',default=[12,20,32,48])
    parser.add_argument('--seconds',type=float,default=90,help='Per-reducer time limit; archived fields are never re-solved')
    parser.add_argument('--lc',type=float,help='Explicit characteristic length for a unit-scaling control on unchanged topology')
    args=parser.parse_args()
    if args.orders!=sorted(set(args.orders)) or not all(1<=q<=100 for q in args.orders):
        parser.error('orders must be unique, increasing, and in [1,100]')
    if not np.isfinite(args.seconds) or args.seconds<=0:
        parser.error('seconds must be finite and positive')
    if args.lc is not None and (not np.isfinite(args.lc) or args.lc<=0):
        parser.error('Lc must be finite and positive')
    run(args)
