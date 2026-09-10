#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bounded p5 study of interface ownership, separate from mesh/PDE accuracy.

Unresolved triangles are tagged in a diagnostic copy, never removed. Their
positive interface energy gives an upper bound on fixed-FE-field redistribution
between slots of the same physical family/conductor. It is not a bound on PDE
error, the domain defect, or final device observables.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
OFFSET=10000


def numeric_rows(path):
    with path.open() as stream:
        return [{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def physical_group(attribute):
    family=attribute//1000
    if family in (5,6):
        return f"{'MS' if family==5 else 'MA'}:conductor-{attribute%100}"
    if family==3:
        return f"SA:{'etched' if attribute>=3100 else 'unetched'}"
    raise ValueError(f"Unsupported fabricated interface {attribute}")


def make_probe(mesh,attributes,output,prescribed=False):
    canonical=lambda a:a%OFFSET
    conductors=sorted({canonical(a)%100 for a in attributes if canonical(a)//1000 in (5,6)})
    if not conductors:
        raise ValueError("No conductors")
    dielectric=[]
    mapping={}
    for index,attribute in enumerate(sorted(attributes),1):
        original=canonical(attribute)
        family=original//1000
        name,permittivity={3:("SA",4.),5:("MS",11.47),6:("MA",10.)}[family]
        dielectric.append(dict(Index=index,Attributes=[attribute],Type=name,Thickness=.002,
                               Permittivity=permittivity,LossTan=.001))
        mapping[index]=dict(Attribute=attribute,CanonicalAttribute=original,
                            Group=physical_group(original),Unresolved=attribute>=OFFSET)
    config={
        "Problem":{"Type":"Electrostatic","Verbose":1,"Output":str(output),
                   "OutputFormats":{"Paraview":False,"GridFunction":False}},
        "Model":{"Mesh":str(mesh),"L0":1e-6,"Refinement":{"MaxIts":0}},
        "Domains":{"Materials":[{"Attributes":[1],"Permittivity":11.47},
                                   {"Attributes":[2],"Permittivity":1.}]},
        "Boundaries":{"Ground":{"Attributes":[1]},
                      "Terminal":[{"Index":i,"Attributes":[a for a in attributes
                                    if canonical(a)//1000 in (5,6) and canonical(a)%100==c]}
                                  for i,c in enumerate(conductors,1)],
                      "Postprocessing":{"Dielectric":dielectric}},
        "Solver":{"Order":5,"Electrostatic":{"Save":0},
                  "Linear":{"Type":"BoomerAMG","KSPType":"CG","Tol":1e-8,"MaxIts":500,
                            "EstimatorTol":.5,"EstimatorMaxIts":5,"EstimatorMG":True}},
    }
    if prescribed:
        zero=output.parent/'zero-trace.csv'
        zero.write_text('x,y,z,V,triangle\n0,0,0,0,1\n1,0,0,0,1\n0,1,0,0,1\n')
        terminals=config['Boundaries'].pop('Terminal')
        config['Boundaries']['Ground']['Attributes']=sorted(a for a in attributes if canonical(a)//1000 in (5,6))
        config['Boundaries']['PrescribedPotential']=[dict(Index=t['Index'],Attributes=[1],
            TerminalAttributes=t['Attributes'],DataFile=str(zero)) for t in terminals]
        config['Solver']['Electrostatic'].update(ResponseMatrix=True,AggregateResponseMatrix=True)
        # The response-matrix writer requires localization metadata. The study
        # reads Q_total_ij (whole-interface energy), not this optional core term.
        edge_attributes=[a for a in attributes if canonical(a)//1000==3]
        for interface in config['Boundaries']['Postprocessing']['Dielectric']:
            interface.update(LocalizeEdgeEnergy=True,SaveLocalEdgeEnergy=False,
                             EdgeAttributes=edge_attributes,EdgeDistances=[0.5],
                             EdgeFrameNormal=[0.,0.,1.],EdgeExcludeAttributes=[1])
    return config,mapping


def energies(output,mapping):
    if (output/'domain-E.csv').exists():
        domain={int(r['i']):r['E_elec (J)'] for r in numeric_rows(output/'domain-E.csv')}
        surfaces={int(r['i']):{index:r[f'p_surf[{index}]']*domain[int(r['i'])] for index in mapping}
                  for r in numeric_rows(output/'surface-Q.csv')}
    else:
        domain={int(r['basis_i']):r['Q_ij (J)'] for r in numeric_rows(output/'domain-response-matrix.csv')
                if r['basis_i']==r['basis_j']}
        surfaces={source:{} for source in domain}
        for r in numeric_rows(output/'surface-response-matrix.csv'):
            if r['basis_i']==r['basis_j']:
                surfaces[int(r['basis_i'])][int(r['interface'])]=r['Q_total_ij (J)']
    if not domain or domain.keys()!=surfaces.keys() or not all(math.isfinite(q) and q>0 for q in domain.values()):
        raise ValueError('Missing, non-finite, or zero-energy conductor state')
    result={}
    for source,row in surfaces.items():
        if set(row)!=set(mapping):
            raise ValueError('Surface interface index set differs from the probe configuration')
        slots={};groups={};uncertain={}
        for index,info in mapping.items():
            q=row[index]
            if not math.isfinite(q) or q<0:
                raise ValueError("Invalid interface energy cannot bound redistribution")
            a=str(info['CanonicalAttribute']);g=info['Group']
            slots[a]=slots.get(a,0.)+q
            groups[g]=groups.get(g,0.)+q
            uncertain[g]=uncertain.get(g,0.)+(q if info['Unresolved'] else 0.)
        result[str(source)]={"DomainEnergy":domain[source],"Slots":slots,"Groups":groups,
                             "UnresolvedEnergy":uncertain,
                             "GroupNormalizedPartitionBound":{g:uncertain[g]/q if q else 0. for g,q in groups.items()},
                             "SlotRelativePartitionBounds":{a:uncertain[physical_group(int(a))]/q if q else None for a,q in slots.items()}}
    return result


def relative(x,y):
    return abs(y-x)/abs(x) if x else (0. if y==0 else None)


def compare_steps(first,second):
    if first.keys()!=second.keys():
        raise ValueError("Source sets changed")
    report={"Domain":0.,"Groups":0.,"Slots":0.}
    for source,a in first.items():
        b=second[source]
        report['Domain']=max(report['Domain'],relative(a['DomainEnergy'],b['DomainEnergy']))
        for group in ('Groups','Slots'):
            if a[group].keys()!=b[group].keys():
                raise ValueError("Physical group/slot set changed")
            changes=[relative(x,b[group][key]) for key,x in a[group].items()]
            if any(v is None for v in changes):
                raise ValueError("A zero-energy slot became nonzero")
            report[group]=max(report[group],max(changes,default=0.))
    return report


def run(args):
    args.root.mkdir(parents=True,exist_ok=False)
    tools=args.root/'tools';tools.mkdir()
    names=['mesh_spatial_coupon.jl','mesh_graded_tet_experiment.jl','frozen_volume_study.jl','label_interface_patches.jl',
           'interface_ownership.jl','ownership_bernstein.jl','graded_curve_distance.jl','graded_size_points.jl',
           'tag_partition_uncertainty.jl','run_bounded_mesher.py','run_partition_refinement.py']
    for name in names:shutil.copy2(HERE/name,tools/name)
    summary={"Version":1,"Scope":"Fixed-FE-field partition bounds and independent convergence observations",
             "LibraryQualified":False,"Order":5,"SolveTolerance":1e-8,"Cases":[],
             "Ranks":args.ranks,"ArchiveStates":args.archive_states,
             "ResourceLimits":{"SolveSeconds":args.solve_seconds,"SolveMemoryGiB":args.solve_memory_gib},
             "Criteria":{"GroupNormalizedPartitionBound":args.partition_tol,
                         "ObservedDomainChange":5e-4,"ObservedGroupAndSlotChange":5e-3},
             "PartitionBoundDefinition":"Each slot absolute ownership error <= unresolved group energy, for the fixed FE field and interface functional; not a PDE-error bound",
             "SHA256":{name:hashlib.sha256((tools/name).read_bytes()).hexdigest() for name in names}}
    summary['SHA256']['PalaceBinary']=hashlib.sha256(args.palace.read_bytes()).hexdigest()
    base_env={k:v for k,v in os.environ.items() if not k.startswith(('TET_','PALACE_RESPONSE_'))}
    base_env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',JULIA_NUM_THREADS='1',
                    TET_GEOMETRY_ORDER='1',TET_EDGE_TANGENT_SIZE='.05',TET_SURFACE_ALGORITHM='5',
                    TET_HXT_QUALITY='.1',TET_SLOT_MINIMUM_SIZE=str(args.slot_minimum))
    def persist():
        (args.root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    def bounded(command,log,env,seconds,memory):
        subprocess.run([sys.executable,str(tools/'run_bounded_mesher.py'),'--seconds',str(seconds),
                        '--memory-gib',str(memory),'--log',str(log),'--',*command],
                       env=env,check=True,timeout=seconds+15,stdout=subprocess.DEVNULL)
    baseline=None
    if args.baseline:
        baseline=json.loads((args.baseline/'summary.json').read_text())
        if baseline['Order']!=5 or baseline['SolveTolerance']!=1e-8 or baseline.get('ArchiveStates',False)!=args.archive_states:
            raise ValueError('Incompatible baseline solve settings')
        summary['Baseline']=str(args.baseline)
        summary['SHA256']['BaselineSummary']=hashlib.sha256((args.baseline/'summary.json').read_bytes()).hexdigest()
    for name in args.case:
        source=args.inputs/name;case=args.root/name;case.mkdir()
        inputs=case/'input';inputs.mkdir()
        for f in ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv','process.toml',
                  'fabricated-expected-interfaces.csv','fabricated-geometry.msh.cad-measures.csv'):
            shutil.copy2(source/f,inputs/f)
        lane={'Case':name,'Steps':[]};summary['Cases'].append(lane)
        hints={};previous=None;start=0
        if baseline:
            old=next(c for c in baseline['Cases'] if c['Case']==name)['Steps'][0]
            if old['Iteration']!=0:
                raise ValueError('Expected a level-zero baseline')
            for f in inputs.iterdir():
                if f.read_bytes()!=(args.baseline/name/'input'/f.name).read_bytes():
                    raise ValueError('Baseline geometry/process inputs differ')
            shutil.copytree(args.baseline/name/'level0',case/'level0')
            lane['Steps'].append(dict(old,ReusedFrom=str(args.baseline/name/'level0')))
            previous=old['Energies'];start=1
            with (case/'level0/fabricated.msh.interface-partition.csv.refine.csv').open() as stream:
                for r in csv.DictReader(stream):
                    p=tuple(float(r[d]) for d in ('x','y','z'));h=float(r['size'])
                    hints[p]=min(hints.get(p,h),h)
        try:
            for iteration in range(start,args.iterations):
                step=case/f'level{iteration}';step.mkdir();env=dict(base_env)
                if hints:
                    path=step/'refinement-points.csv'
                    with path.open('w',newline='') as stream:
                        writer=csv.writer(stream);writer.writerow(['x','y','z','size'])
                        writer.writerows((*p,h) for p,h in sorted(hints.items()))
                    env['TET_SLOT_REFINEMENT_POINTS']=str(path)
                mesh=step/'fabricated.msh'
                command=[args.julia,'--startup-file=no',f'--project={args.julia_project}',
                         str(tools/'mesh_graded_tet_experiment.jl'),str(inputs),'fabricated',str(mesh),
                         '1.0','.02','.25','--process',str(inputs/'process.toml'),
                         '--element-interface-slots','--reference-measures',str(inputs/'fabricated-geometry.msh.cad-measures.csv'),
                         '--expected-interfaces',str(inputs/'fabricated-expected-interfaces.csv')]
                bounded(command,step/'mesher.log',env,60,6)
                diagnostic=step/'partition-audit.msh'
                bounded([args.julia,'--startup-file=no',f'--project={args.julia_project}',
                         str(tools/'tag_partition_uncertainty.jl'),str(mesh),str(diagnostic)],
                        step/'tagger.log',env,20,4)
                with Path(str(diagnostic)+'.attributes.csv').open() as stream:
                    attributes=[int(r['attribute']) for r in csv.DictReader(stream)]
                config,mapping=make_probe(diagnostic,attributes,step/'postpro',args.archive_states)
                config_path=step/'probe.json';config_path.write_text(json.dumps(config,indent=2)+'\n')
                (step/'interface-map.json').write_text(json.dumps(mapping,indent=2)+'\n')
                if args.archive_states:
                    archive=step/'archives';archive.mkdir()
                    worker_config=json.loads(json.dumps(config))
                    worker_config['Problem']['Output']=str(step/'worker')
                    worker_path=step/'worker.json';worker_path.write_text(json.dumps(worker_config,indent=2)+'\n')
                    worker_env=dict(env,PALACE_RESPONSE_ARCHIVE_DIR=str(archive),PALACE_RESPONSE_ARCHIVE_ONLY='1')
                    bounded(['mpirun','-n',str(args.ranks),'-x','PALACE_RESPONSE_ARCHIVE_DIR','-x','PALACE_RESPONSE_ARCHIVE_ONLY',
                             str(args.palace),str(worker_path)],step/'worker.log',worker_env,args.solve_seconds,args.solve_memory_gib)
                    reducer_env=dict(env,PALACE_RESPONSE_ARCHIVE_DIR=str(archive),PALACE_RESPONSE_REDUCE_ONLY='1',PALACE_RESPONSE_BLOCK_SIZE='2')
                    bounded(['mpirun','-n',str(args.ranks),'-x','PALACE_RESPONSE_ARCHIVE_DIR','-x','PALACE_RESPONSE_REDUCE_ONLY',
                             '-x','PALACE_RESPONSE_BLOCK_SIZE',str(args.palace),str(config_path)],step/'reducer.log',reducer_env,args.solve_seconds,args.solve_memory_gib)
                else:
                    bounded(['mpirun','-n',str(args.ranks),str(args.palace),str(config_path)],
                            step/'solver.log',env,args.solve_seconds,args.solve_memory_gib)
                result=energies(step/'postpro',mapping)
                metadata=json.loads(Path(str(mesh)+'.metadata.json').read_text())
                with Path(str(mesh)+'.interface-partition.csv').open() as stream:
                    partition=list(csv.DictReader(stream))
                record={'Iteration':iteration,'Elements':metadata['VolumeElementCount'],
                        'RefinementPoints':len(hints),'Energies':result,
                        'MaxUnresolvedAreaFraction':max(float(r['unresolved_fraction']) for r in partition),
                        'MaxSampledAmbiguousAreaFraction':max(float(r['ambiguous_fraction']) for r in partition),
                        'MaxGroupNormalizedPartitionBound':max(v for s in result.values() for v in s['GroupNormalizedPartitionBound'].values())}
                if previous is not None:record['Change']=compare_steps(previous,result)
                record['PartitionBoundPassed']=record['MaxGroupNormalizedPartitionBound']<=args.partition_tol
                change=record.get('Change')
                record['ObservedProbeConvergence']=bool(change and change['Domain']<=5e-4 and
                                                        change['Groups']<=5e-3 and change['Slots']<=5e-3)
                # Once per case, independently check that diagnostic splitting did
                # not change the PDE or the sum of per-slot interface energies.
                if iteration==0:
                    original_attributes=sorted({a%OFFSET for a in attributes})
                    original_config,original_map=make_probe(mesh,original_attributes,step/'control',args.archive_states)
                    path=step/'control.json';path.write_text(json.dumps(original_config,indent=2)+'\n')
                    bounded(['mpirun','-n',str(args.ranks),str(args.palace),str(path)],step/'control.log',env,args.solve_seconds,args.solve_memory_gib)
                    changes=compare_steps(energies(step/'control',original_map),result)
                    record['DiagnosticControlChange']=changes
                    if max(changes.values())>1e-6:raise ValueError('Diagnostic tags changed the solve/energies')
                lane['Steps'].append(record);persist();print(json.dumps({k:v for k,v in record.items() if k!='Energies'}),flush=True)
                if record['PartitionBoundPassed'] and record['ObservedProbeConvergence']:
                    lane['StudyCriteriaReached']=True
                    break
                with Path(str(mesh)+'.interface-partition.csv.refine.csv').open() as stream:
                    for r in csv.DictReader(stream):
                        p=tuple(float(r[d]) for d in ('x','y','z'));h=float(r['size'])
                        hints[p]=min(hints.get(p,h),h)
                previous=result
            lane['Completed']=True
        except (ValueError,subprocess.SubprocessError,KeyError) as error:
            lane['Error']=str(error);lane['Completed']=False
        persist()
    return all(lane.get('Completed') for lane in summary['Cases'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,help='Reuse an identical, completed level zero without another solve')
    parser.add_argument('--archive-states',action='store_true',help='Exact streaming/reduction of prescribed one-volt conductor states')
    parser.add_argument('--palace',type=Path,required=True)
    parser.add_argument('--julia',default='julia')
    parser.add_argument('--julia-project',type=Path,required=True)
    parser.add_argument('--case',action='append',required=True)
    parser.add_argument('--iterations',type=int,default=3,choices=range(1,5))
    parser.add_argument('--ranks',type=int,default=2,choices=range(1,7))
    parser.add_argument('--slot-minimum',type=float,default=.002)
    parser.add_argument('--partition-tol',type=float,default=.005)
    parser.add_argument('--solve-seconds',type=float,default=120)
    parser.add_argument('--solve-memory-gib',type=float,default=24)
    args=parser.parse_args()
    for value,name in ((args.partition_tol,'partition tolerance'),(args.solve_seconds,'solve seconds'),
                       (args.solve_memory_gib,'solve memory')):
        if not math.isfinite(value) or value<=0:
            parser.error(f'{name} must be finite and positive')
    raise SystemExit(0 if run(args) else 1)
