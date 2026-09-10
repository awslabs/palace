#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Combine same-host benchmark ledgers without treating interrupted times as exact."""
import argparse
import json
from pathlib import Path
import numpy as np
from run_matched_mesh_benchmark import matrix,norm,difference,diagonal_differences


def worst_energy_difference(reference,candidate):
    ids=sorted({i for i,j in reference})
    a=np.array([[reference[(i,j)] for j in ids] for i in ids])
    b=np.array([[candidate[(i,j)] for j in ids] for i in ids])
    values,vectors=np.linalg.eigh(a)
    # Do not silently truncate a gauge or weak mode to claim a worst-case bound.
    if not len(values) or values[-1]<=0 or values[0]<=1e-12*values[-1]:return None
    inverse_root=vectors/np.sqrt(values)
    return float(np.max(np.abs(np.linalg.eigvalsh(inverse_root.T@(b-a)@inverse_root))))


def summarize(roots):
    cases={};metadata={};contract=None;input_contract=None
    for root in roots:
        manifest=json.loads((root/'manifest.json').read_text())
        immutable={k:v for k,v in manifest['InputSHA256'].items() if k.startswith(('input/','traces/'))}
        if input_contract is not None and immutable!=input_contract:raise ValueError('Geometry/process/excitations differ')
        input_contract=immutable
        for case in manifest['Cases']:
            if case['Name'] in metadata and metadata[case['Name']]['MeshSHA256']!=case['MeshSHA256']:
                raise ValueError('Conflicting mesh for the same benchmark case')
            metadata[case['Name']]=case
        for filename in ('timings.json','timings-followup.json'):
            path=root/filename
            if not path.exists():continue
            ledger=json.loads(path.read_text())
            identity=tuple(ledger[k] for k in ('Host','ExecutableSHA256','Ranks','Order','SolverTolerance'))+(ledger.get('Binding','core'),)
            if contract is not None and identity!=contract:raise ValueError('Hardware/executable/solver contract differs')
            contract=identity
            for name,data in ledger['Cases'].items():
                if not data.get('Completed'):continue
                cases[name]={'Root':root,'Worker':data['Worker'],'Reducer':data['Reducer'],
                             'TotalSeconds':data['TotalSeconds'],'TimingIsLowerBound':False,
                             'PeakProcessTreeRSSBytes':data['PeakProcessTreeRSSBytes']}
        for path in (root/'recovery').glob('*/summary.json') if (root/'recovery').exists() else []:
            data=json.loads(path.read_text())
            if not data.get('Completed'):continue
            name=data['Case']
            if name in cases:raise ValueError('Recovery would overwrite completed exact timing')
            original_resource=root/'logs'/f'{name}-worker.log.json'
            original=json.loads(original_resource.read_text()) if original_resource.exists() else data.get('InterruptedWorker',{})
            cases[name]={'Root':root,'Worker':{'H1Dofs':data['Reducer']['H1Dofs'],'Sources':data['Sources']},
                         'Reducer':data['Reducer'],'TotalSeconds':data['TotalSecondsLowerBound'],
                         'TimingIsLowerBound':True,'Recovery':data,
                         'PeakProcessTreeRSSBytes':max([original.get('PeakProcessTreeRSSBytes',0),data['Reducer']['PeakProcessTreeRSSBytes']]+[a['PeakProcessTreeRSSBytes'] for a in data['Attempts']])}
    report={'Scope':'Same-host p5 pilot, not complete production library generation',
            'Contract':dict(zip(('Host','ExecutableSHA256','Ranks','Order','SolverTolerance','Binding'),contract or ())),
            'Cases':{},'ThinRawSPRGated':False,'LibraryQualified':False}
    domains={}
    for name,data in cases.items():
        case=metadata[name];root=data['Root']
        domains[name]=matrix(root/'postpro'/f'{name}-reducer/domain-response-matrix.csv')[0]
        entry={'Kind':case['Kind'],'Strategy':case['Strategy'],'Elements':case['Elements'],
               'KappaMax':case['Kappa']['100.000000'],'H1Dofs':data['Worker']['H1Dofs'],
               'SourceTimings':data['Worker']['Sources'],'ReducerSeconds':data['Reducer']['Seconds'],
               'TotalSeconds':data['TotalSeconds'],'TimingIsLowerBound':data['TimingIsLowerBound'],
               'PeakProcessTreeRSSBytes':data['PeakProcessTreeRSSBytes']}
        if 'Seconds' in data['Worker']:entry['WorkerSeconds']=data['Worker']['Seconds']
        report['Cases'][name]=entry
    for name,data in cases.items():
        case=metadata[name];kind=case['Kind'];baseline=f'{kind}-prism'
        if baseline not in cases:continue
        entry=report['Cases'][name];reference=cases[baseline]
        entry['DomainMatrixRelativeDifference']=difference(domains[baseline],domains[name])
        entry['DomainWorstEnergyRelativeDifference']=worst_energy_difference(domains[baseline],domains[name])
        entry['DomainDiagonalRelativeDifferences']=diagonal_differences(domains[baseline],domains[name])
        entry['Speedup']=reference['TotalSeconds']/data['TotalSeconds'] if not data['TimingIsLowerBound'] or name==baseline else None
        entry['SpeedupIsLowerBound']=reference['TimingIsLowerBound'] and name!=baseline and entry['Speedup'] is not None
        entry['H1DofReduction']=reference['Worker']['H1Dofs']/data['Worker']['H1Dofs']
        if kind=='fabricated':
            a=matrix(reference['Root']/'postpro'/f'{baseline}-reducer/surface-response-matrix.csv',True)
            b=matrix(data['Root']/'postpro'/f'{name}-reducer/surface-response-matrix.csv',True)
            entry['InterfaceMatrixRelativeDifferences']={str(i):difference(a[i],b[i]) for i in a}
            entry['InterfaceWorstEnergyRelativeDifferences']={str(i):worst_energy_difference(a[i],b[i]) for i in a}
            entry['InterfaceDiagonalRelativeDifferences']={str(i):diagonal_differences(a[i],b[i]) for i in a}
    if 'thin-prism' in domains and 'fabricated-prism' in domains:
        reference={k:domains['fabricated-prism'][k]-v for k,v in domains['thin-prism'].items()}
        report['DomainDefect']={}
        for strategy in sorted({c['Strategy'] for c in metadata.values()}):
            thin,fab=f'thin-{strategy}',f'fabricated-{strategy}'
            if thin not in domains or fab not in domains:continue
            value={k:domains[fab][k]-v for k,v in domains[thin].items()}
            report['DomainDefect'][strategy]={'RelativeToDefect':difference(reference,value),
                'RelativeToThinDomain':norm({k:value[k]-v for k,v in reference.items()})/norm(domains['thin-prism'])}
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite analysis')
    result=summarize(args.root)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for name,r in result['Cases'].items():
        print(name,r['H1Dofs'],('>' if r['TimingIsLowerBound'] else '')+f"{r['TotalSeconds']:.2f}s",
              (('>' if r.get('SpeedupIsLowerBound') else '')+f"{r['Speedup']:.2f}x") if r.get('Speedup') is not None else 'unknown speedup',r.get('DomainMatrixRelativeDifference'),r.get('InterfaceMatrixRelativeDifferences'))
