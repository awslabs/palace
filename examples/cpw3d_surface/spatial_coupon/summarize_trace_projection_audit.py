#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Group all 135 audited columns using preexisting corrected-model constraints.

Column-normalized nodal Gram spectra are numerical diagnostics, not rigorous rank
certificates. No column is removed or relabeled based on its numerical values.
"""
import argparse
import csv
import json
from pathlib import Path
import shutil

import numpy as np

from prepare_trace_projection_audit import sha


def fixed_groups(contract, library):
    models = [m for m in library['Models'] if m['Name'] == contract['Model']]
    if len(models) != 1 or models[0]['ZeroTraceIndices'] != contract['ZeroTraceIndices']:
        raise ValueError('Corrected model / source-contract constrained mapping mismatch')
    constrained = contract['ZeroTraceIndices']
    if len(constrained) != 40 or len(set(constrained)) != 40 or not set(constrained) <= set(range(1, 136)):
        raise ValueError('Expected fixed 40-constrained/95-free partition')
    return {'all135': list(range(1, 136)),
            'free95': [i for i in range(1, 136) if i not in constrained],
            'constrained40': constrained}


def rank_diagnostic(gram):
    # Nonzero columns are normalized, zero columns retained. This is not null filtering.
    diagonal = np.diag(gram)
    if np.any(diagonal < 0):
        raise ValueError('Negative nodal Gram diagonal')
    scale = np.sqrt(np.where(diagonal > 0, diagonal, 1.))
    normalized = gram/scale[:, None]/scale[None, :]
    spectrum = np.linalg.eigvalsh((normalized+normalized.T)/2)
    largest = max(float(spectrum[-1]), 0.)
    return {'NormalizedGramEigenvalues': spectrum.tolist(),
            'RanksByRelativeGramEigenvalueThreshold': {
                str(tol): int(np.count_nonzero(spectrum > tol*largest))
                for tol in (1e-10, 1e-12, 1e-14)},
            'ZeroDiagonalColumns': int(np.count_nonzero(diagonal == 0)),
            'NodalGramDiagonalRange': [float(diagonal.min()), float(diagonal.max())]}


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def summarize(projection, inputs, library_path):
    contract = json.loads((inputs/'basis-contract.json').read_text())
    library = json.loads(library_path.read_text())
    groups = fixed_groups(contract, library)
    report = {'CorrectedModel': contract['Model'], 'CorrectedModelSHA256': sha(library_path),
              'FixedGroups': groups, 'RankMethod': 'Unweighted nodal Gram, columns normalized to unit Euclidean norm; zero columns retained; Gram thresholds square singular-value thresholds; no rigorous rank certificate',
              'Meshes': {}, 'OmittedSources': []}
    for mesh in ('prism', 'tet14'):
        norm_rows = rows(projection/f'{mesh}-norms.csv')
        final = [r for r in norm_rows if int(r['quadrature_order']) == 20]
        nodal = rows(projection/f'{mesh}-nodal-stats.csv')
        if [int(r['source']) for r in final] != list(range(1, 136)):
            raise ValueError('Incomplete final source results')
        grams = {kind: np.loadtxt(projection/f'{mesh}-{kind}-gram.csv', delimiter=',')
                 for kind in ('raw', 'grounded')}
        result = {'Groups': {}, 'QuadratureRefinement': {}}
        for name, indices in groups.items():
            selected = [final[i-1] for i in indices]
            nd = [nodal[i-1] for i in indices]
            raw = np.array([float(r['raw_interpolation_error_l2'])/float(r['target_l2']) for r in selected])
            grounded = np.array([float(r['grounded_relative_error']) for r in selected])
            def stats(a):
                return {'Min': float(a.min()), 'Median': float(np.median(a)),
                        'Max': float(a.max()), 'WorstSource': indices[int(a.argmax())]}
            select = np.array(indices)-1
            result['Groups'][name] = {
                'RawRelativeError': stats(raw), 'GroundedRelativeError': stats(grounded),
                'RawContactSplitMax': max(float(r['raw_contact_split_max']) for r in nd),
                'RawContactSplitL2Max': max(float(r['raw_contact_l2_split']) for r in nd),
                'RawContactNonzeroAbove1e-10': [int(r['source']) for r in nd if float(r['raw_contact_split_max']) > 1e-10],
                'ZeroGroundedColumns': [int(r['source']) for r in nd if int(r['grounded_nonzero_dofs']) == 0],
                'RankDiagnostics': {kind: rank_diagnostic(g[np.ix_(select, select)]) for kind, g in grams.items()}}
        for a, b in ((12, 16), (16, 20)):
            refinement = {}
            for key in ('target_l2', 'raw_interpolation_error_l2', 'grounded_error_l2'):
                x = np.array([float(r[key]) for r in norm_rows if int(r['quadrature_order']) == a])
                y = np.array([float(r[key]) for r in norm_rows if int(r['quadrature_order']) == b])
                refinement[key] = {'MaxAbsoluteChange': float(np.max(abs(x-y))),
                                   'MaxRelativeChange': float(np.max(abs(x-y)/y))}
            result['QuadratureRefinement'][f'{a}->{b}'] = refinement
        report['Meshes'][mesh] = result
    provenance = projection/'provenance'
    provenance.mkdir(exist_ok=True)
    shutil.copyfile(library_path, provenance/'process-library.json')
    (projection/'grouped-summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({mesh: {name: {'GroundedRelativeError': group['GroundedRelativeError'],
                                   'RawContactSplitMax': group['RawContactSplitMax'],
                                   'RawRank': group['RankDiagnostics']['raw']['RanksByRelativeGramEigenvalueThreshold'],
                                   'GroundedRank': group['RankDiagnostics']['grounded']['RanksByRelativeGramEigenvalueThreshold']}
                            for name, group in data['Groups'].items()}
                      for mesh, data in report['Meshes'].items()}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('projection', type=Path)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('corrected_library', type=Path)
    args = parser.parse_args()
    summarize(args.projection, args.inputs, args.corrected_library)
