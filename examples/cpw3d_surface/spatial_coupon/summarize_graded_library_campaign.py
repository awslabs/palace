#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Summarize generation state without converting numerical completion into qualification."""
import json
import sys
from pathlib import Path
from run_graded_library_case import save


def load(path,default):
    return json.loads(path.read_text()) if path.exists() else default


def main():
    root=Path(sys.argv[1]).resolve();manifest=load(root/'campaign.json',{})
    reference=load(Path(manifest['ReferenceManifest']),{})
    old={c['Key']:c for c in reference['Cases']}
    cases=[];new_node_hours=0.;old_node_hours=0.
    for c in manifest['Cases']:
        state=load(Path(c['CaseDirectory'])/'status.json',{'State':'NotStarted'})
        mesh=load(Path(c['Mesh']).parent/'mesh-state.json',{'State':'NotStarted'}) if c['MeshRequired'] else {'State':'RetainedControl'}
        previous=load(Path(old[c['Key']]['CaseDirectory'])/'status.json',{})
        stages={s['Stage']:s for s in state.get('Stages',[])}
        original={s['Stage']:s for s in previous.get('Stages',[])}
        elapsed=sum(s.get('WallSeconds',0.) for k,s in stages.items() if k in ('worker','reducer'))
        baseline=sum(s.get('WallSeconds',0.) for k,s in original.items() if k in ('worker','reducer'))
        if state['State']=='Completed':new_node_hours+=elapsed*c['Nodes']/3600
        old_node_hours+=baseline*c['Nodes']/3600
        row={'Case':c['Key'],'Order':c['Order'],'Sources':c['Sources'],'OriginalSources':c['OriginalSources'],
             'Nodes':c['Nodes'],'Ranks':c['Ranks'],'MeshState':mesh['State'],'ResponseState':state['State'],
             'MeshElements':mesh.get('Elements'),'MaximumKappa':mesh.get('MaximumKappa'),
             'WorkerSeconds':stages.get('worker',{}).get('WallSeconds'),
             'ReducerSeconds':stages.get('reducer',{}).get('WallSeconds'),
             'OriginalWorkerReducerSeconds':baseline,'Error':state.get('Error',mesh.get('Error'))}
        if state['State']=='Completed':
            row['WorkerReducerSeconds']=elapsed
            row['HistoricalCaseSpeedRatio']=baseline/elapsed
        iterations=stages.get('worker',{}).get('Iterations',[])
        if iterations:row['Iterations']={'Min':min(iterations),'Max':max(iterations),'Mean':sum(iterations)/len(iterations)}
        cases.append(row)
    report={'Version':1,'Complete':all(c['ResponseState']=='Completed' for c in cases),
            'CompletedCases':sum(c['ResponseState']=='Completed' for c in cases),
            'CaseCount':len(cases),'SourceCount':manifest['SourceCount'],
            'OriginalSourceCount':manifest['OriginalSourceCount'],'Cases':cases,
            'CompletedResponseStageNodeHours':new_node_hours,'OriginalResponseStageNodeHours':old_node_hours,
            'LibraryManifestPublished':(root/'new-library/process-library.json').exists(),
            'LibraryQualified':False,
            'ComparisonScope':'Complete box/cap source correction changes the workload; historical ratios are not pure mesh-only speedups. Meshing and control-pilot costs are separate.'}
    save(root/'summary.json',report)
    for c in cases:print(c['Case'],c['MeshState'],c['ResponseState'],c.get('WorkerReducerSeconds'),c['Error'] or '')
    print('Completed',report['CompletedCases'],'/',report['CaseCount'],'cases; sources',report['SourceCount'])


if __name__=='__main__':main()
