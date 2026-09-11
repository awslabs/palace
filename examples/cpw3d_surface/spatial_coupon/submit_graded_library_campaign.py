#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Submit unique mesh/response jobs, keeping the user's total PBS job count below 40."""
import argparse
import fcntl
import getpass
import json
import re
import subprocess
import time
from pathlib import Path
from run_graded_library_case import save


def active_jobs():
    data=json.loads(subprocess.check_output(['qstat','-f','-F','json'],text=True))
    owner=getpass.getuser()+'@'
    return {key:job for key,job in data.get('Jobs',{}).items()
            if job.get('Job_Owner','').startswith(owner) and job.get('job_state')!='F'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--phase',choices=('controls','meshes','responses','all'),required=True)
    parser.add_argument('--keys',help='Optional comma-separated case keys')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args();root=args.root.resolve()
    manifest=json.loads((root/'campaign.json').read_text())
    keys=set(args.keys.split(',')) if args.keys else None
    selected=[c for c in manifest['Cases'] if keys is None or c['Key'] in keys]
    if keys and keys!={c['Key'] for c in selected}:raise ValueError('Unknown case key')
    metadata=root/'audit/device-preflight/surface-response-requirements.json'
    if not json.loads(metadata.read_text())['Complete']:raise ValueError('Strict metadata preflight incomplete')
    with (root/'submit.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=root/'submissions.json'
        submitted=json.loads(path.read_text()) if path.exists() else {}
        actions=[]
        if args.phase in ('meshes','all'):
            for case in selected:
                if not case['MeshRequired']:continue
                state=Path(case['Mesh']).parent/'mesh-state.json'
                if state.exists() and json.loads(state.read_text())['State']=='Completed':continue
                tag='mesh-'+case['Key']
                if tag not in submitted:actions.append((tag,None))
        if args.phase in ('controls','responses','all'):
            for case in selected:
                if args.phase=='controls' and case['MeshRequired']:continue
                tag='response-'+case['Key']
                if tag in submitted:continue
                dependency=None
                if case['MeshRequired']:
                    state=Path(case['Mesh']).parent/'mesh-state.json'
                    if not state.exists() or json.loads(state.read_text())['State']!='Completed':
                        dependency='mesh-'+case['Key']
                        if dependency not in submitted and dependency not in {tag for tag,_ in actions}:
                            raise ValueError(f'No completed or scheduled mesh for {case["Key"]}')
                actions.append((tag,dependency))
        current=active_jobs()
        if len(current)+len(actions)>40:raise RuntimeError(f'Would exceed 40 user jobs: {len(current)} active/queued + {len(actions)} new')
        if args.dry_run:
            print(json.dumps({'UserJobs':len(current),'Actions':actions},indent=2));return
        for tag,dependency in actions:
            if len(active_jobs())>=40:raise RuntimeError('User job limit reached; remaining jobs not submitted')
            command=['qsub']
            if dependency:command+=['-W','depend=afterok:'+submitted[dependency]['Job']]
            command.append(str(root/'jobs'/(tag+'.pbs')))
            job=subprocess.check_output(command,text=True).strip()
            if not re.fullmatch(r'\d+\.[A-Za-z0-9._-]+',job):raise RuntimeError('Ambiguous qsub receipt; reconcile before retry')
            submitted[tag]={'Job':job,'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'Dependency':dependency,'Command':command}
            save(path,submitted);print(tag,job,flush=True)


if __name__=='__main__':main()
