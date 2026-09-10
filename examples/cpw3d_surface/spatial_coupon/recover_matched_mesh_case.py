#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Recover interrupted pilot sources without re-solving certified completed sources.

Only the controlled pilot's common matching boundary and grounded-conductor state
layout is supported. Recovered matrices are valid, but interrupted timing is NOT
reported as an uninterrupted generation time or an exact speedup.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from run_matched_mesh_benchmark import digest,timing


def run(args):
    root=args.root.resolve();os.chdir(root)
    results=json.loads(Path('timings.json').read_text())
    manifest=json.loads(Path('manifest.json').read_text())
    case=next(c for c in manifest['Cases'] if c['Name']==args.case)
    if results['Cases'][args.case].get('Completed'):raise ValueError('Case already completed')
    binary=Path(results['Executable'])
    if digest(binary)!=results['ExecutableSHA256'] or digest(Path(case['Mesh']))!=case['MeshSHA256']:
        raise ValueError('Executable or mesh changed')
    config=json.loads(Path(f'configs/{args.case}-worker.json').read_text())
    sources=config['Boundaries']['PrescribedPotential']
    ground=set(config['Boundaries']['Ground']['Attributes'])
    if not all(s['Attributes']==[1] and set(s.get('TerminalAttributes',[]))<=ground for s in sources):
        raise ValueError('Source subsetting would change the essential-boundary contract')
    old_log=Path(f'logs/{args.case}-worker.log')
    resource=json.loads(Path(str(old_log)+'.json').read_text())
    if resource['StopReason']!='timeout':raise ValueError('Only an explicit timed-out worker is supported')
    completed={s['Index']:dict(s,AfterRestart=False) for s in timing(old_log)['Sources']}
    archive=root/'archives'/args.case
    ranks=results['Ranks']
    for index in completed:
        if not all((archive/f'source-{index:06d}-rank-{rank:06d}-V.bin').exists() for rank in range(ranks)):
            raise ValueError('Missing archive for a timing-certified source')
    output=root/'recovery'/args.case;output.mkdir(parents=True,exist_ok=False)
    report={'Case':args.case,'Completed':False,'InterruptedWorkerSeconds':resource['Seconds'],
            'InterruptedWorker':resource,
            'WorkerSecondsLowerBound':resource['Limits']['Seconds'],
            'Sources':list(completed.values()),'Attempts':[],'ExactTotalSecondsAvailable':False}
    (output/'worker-original.json').write_text(json.dumps(config,indent=2)+'\n')
    def save(): (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    env={k:v for k,v in os.environ.items() if not k.startswith('PALACE_RESPONSE_')}
    def execute(path,log,stage):
        child=dict(env,PALACE_RESPONSE_ARCHIVE_DIR=str(archive))
        if stage=='worker':child.update(PALACE_RESPONSE_ARCHIVE_ONLY='1',PALACE_RESPONSE_SOURCE_TIMING='1')
        else:child.update(PALACE_RESPONSE_REDUCE_ONLY='1',PALACE_RESPONSE_BLOCK_SIZE=str(len(sources)))
        exports=[x for k in child if k.startswith('PALACE_RESPONSE_') for x in ('-x',k)]
        command=['mpirun','--bind-to','core','-n',str(ranks),*exports,str(binary),str(path)]
        subprocess.run([sys.executable,'tools/run_bounded_mesher.py','--seconds',str(args.seconds),
                        '--memory-gib',str(results['Limits']['MemoryGiB']),'--log',str(log),'--',*command],
                       env=child,check=True,timeout=args.seconds+20)
        return {**json.loads(Path(str(log)+'.json').read_text()),**timing(log)}
    save()
    for source in sources:
        index=source['Index']
        if index in completed:continue
        partial=output/f'old-source-{index}';partial.mkdir()
        for path in archive.glob(f'source-{index:06d}-*'):shutil.move(str(path),partial/path.name)
        c=json.loads(json.dumps(config));c['Boundaries']['PrescribedPotential']=[source]
        c['Problem']['Output']=str(output/f'worker-{index}')
        path=output/f'worker-{index}.json';path.write_text(json.dumps(c,indent=2)+'\n')
        measured=execute(path,output/f'worker-{index}.log','worker')
        if len(measured['Sources'])!=1 or measured['Sources'][0]['Index']!=index:
            raise ValueError('Missing recovered source timing')
        report['Attempts'].append(measured)
        completed[index]=dict(measured['Sources'][0],AfterRestart=True)
        report['Sources']=[completed[i] for i in sorted(completed)];save()
    reducer=execute(root/f'configs/{args.case}-reducer.json',output/'reducer.log','reducer')
    report['Reducer']=reducer;report['Completed']=True
    report['TotalSecondsLowerBound']=report['WorkerSecondsLowerBound']+reducer['Seconds']
    report['ActualSecondsIncludingInterruptedWork']=resource['Seconds']+sum(a['Seconds'] for a in report['Attempts'])+reducer['Seconds']
    save()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--case',required=True)
    p.add_argument('--seconds',type=float,default=600)
    run(p.parse_args())
