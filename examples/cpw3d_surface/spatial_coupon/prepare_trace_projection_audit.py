#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validate the complete trace bank and expose the frozen production coefficient to an audit.

Writes only to a caller-specified scratch directory. Does not alter production sources.
The generated header is the actual class, with private changed to public for audit access.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_trace(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows or set(rows[0]) != {'x', 'y', 'z', 'V', 'triangle'}:
        raise ValueError('Expected x,y,z,V,triangle data')
    groups = {}
    for row in rows:
        data = np.array([float(row[k]) for k in ('x', 'y', 'z', 'V', 'triangle')])
        if not np.isfinite(data).all() or data[4] <= 0 or data[4] != int(data[4]):
            raise ValueError('Nonfinite data or invalid triangle ID')
        groups.setdefault(int(data[4]), []).append(data[:4])
    if any(len(v) != 3 for v in groups.values()):
        raise ValueError('Every triangle must have exactly three rows')
    ids = sorted(groups)
    return ids, np.array([groups[k] for k in ids])


def clip_polygon(poly, triangle):
    """Convex 2D clipping; used only for source-data overlap validation."""
    def cross(a, b):
        return a[0]*b[1]-a[1]*b[0]
    orientation = np.sign(cross(triangle[1]-triangle[0], triangle[2]-triangle[0]))
    for a, b in zip(triangle, np.roll(triangle, -1, axis=0)):
        result = []
        if not len(poly):
            break
        p = poly[-1]
        dp = orientation*cross(b-a, p-a)
        for q in poly:
            dq = orientation*cross(b-a, q-a)
            if (dp >= 0) != (dq >= 0):
                result.append(p + dp/(dp-dq)*(q-p))
            if dq >= 0:
                result.append(q)
            p, dp = q, dq
        poly = result
    return np.asarray(poly)


def area(poly):
    if len(poly) < 3:
        return 0.
    p = np.asarray(poly)-poly[0]
    return abs(np.sum(p[:, 0]*np.roll(p[:, 1], -1)-p[:, 1]*np.roll(p[:, 0], -1)))/2


def validate_bank(geometry, values):
    """values shape (triangle, vertex, source); reject ambiguous P1 source banks."""
    if not np.isfinite(geometry).all() or not np.isfinite(values).all():
        raise ValueError('Nonfinite data')
    normals = np.cross(geometry[:, 1]-geometry[:, 0], geometry[:, 2]-geometry[:, 0])
    if np.any(np.linalg.norm(normals, axis=1) <= 1e-14):
        raise ValueError('Degenerate triangle')
    axes = np.argmax(abs(normals), axis=1)
    for i, tri in enumerate(geometry):
        axis = axes[i]
        if np.ptp(tri[:, axis]) > 1e-11:
            raise ValueError('This bounded audit requires axis-aligned planar source faces')
        xy = np.delete(tri, axis, axis=1)
        for j in range(i):
            other = geometry[j]
            if axes[j] != axis or abs(other[0, axis]-tri[0, axis]) > 1e-11:
                continue
            uv = np.delete(other, axis, axis=1)
            if np.any(xy.max(axis=0) <= uv.min(axis=0)) or np.any(uv.max(axis=0) <= xy.min(axis=0)):
                continue
            if area(clip_polygon(list(xy), uv)) > 1e-12:
                raise ValueError(f'Overlapping source triangle interiors: {j+1}, {i+1}')
    points, point_values = [], []
    seen = {}
    for p, v in zip(geometry.reshape(-1, 3), values.reshape(-1, values.shape[-1])):
        key = tuple(np.round(p, 11))
        if key in seen:
            if np.max(abs(point_values[seen[key]]-v)) > 1e-10:
                raise ValueError('Inconsistent shared vertex values')
        else:
            seen[key] = len(points)
            points.append(p)
            point_values.append(v)
    points, point_values = np.array(points), np.array(point_values)
    # Include T-junction vertices, not just matching triangle vertex IDs.
    max_jump = 0.
    for tri, val, normal in zip(geometry, values, normals):
        ab = np.column_stack((tri[1]-tri[0], tri[2]-tri[0]))
        uv = np.linalg.lstsq(ab, (points-tri[0]).T, rcond=None)[0].T
        bary = np.column_stack((1-uv.sum(axis=1), uv))
        plane = abs((points-tri[0])@normal)/np.linalg.norm(normal)
        inside = (plane <= 1e-11) & (bary.min(axis=1) >= -1e-11)
        max_jump = max(max_jump, float(np.max(abs(bary[inside]@val-point_values[inside]))))
    if max_jump > 1e-9:
        raise ValueError(f'Inconsistent source values at T-junctions: {max_jump}')
    return {'Triangles': len(geometry), 'UniqueVertices': len(points),
            'MaximumVertexOrTJunctionJump': max_jump, 'PositiveAreaOverlaps': 0}


def prepare(inputs, palace, output):
    # Historical model-specific diagnostic.  It is intentionally outside the
    # geometry-independence gate; new generic producers use trace_audit_contract.py.
    config = json.loads((inputs/'template.json').read_text())
    sources = config['Boundaries']['PrescribedPotential']
    if len(sources) != 135 or [s['Index'] for s in sources] != list(range(1, 136)):
        raise ValueError('Requires the complete ordered 135-source bank')
    if config['Boundaries']['Ground']['Attributes'] != [5001, 6001]:
        raise ValueError('Unexpected ground attributes')
    if config['Solver']['Order'] != 4 or config['Boundaries'].get('Terminal'):
        raise ValueError('Unexpected solver order or terminal states')
    all_values, mapping, geometry, triangle_ids = [], [], None, None
    contract = json.loads((inputs/'basis-contract.json').read_text())
    reverse = {v: int(k) for k, v in contract['OldToNewSourceIndices'].items()}
    for s in sources:
        if s['Attributes'] != [1] or s.get('TerminalAttributes'):
            raise ValueError('Unexpected source attributes or conductor activation')
        path = inputs/'traces'/Path(s['DataFile']).name
        if path.name != f"basis-{s['Index']:04}.csv":
            raise ValueError('Unexpected source filename mapping')
        ids, data = read_trace(path)
        if geometry is None:
            triangle_ids, geometry = ids, data[:, :, :3]
        if ids != triangle_ids or not np.array_equal(geometry, data[:, :, :3]):
            raise ValueError('Source geometry/order differs between CSVs')
        all_values.append(data[:, :, 3])
        mapping.append({'Index': s['Index'], 'CSV': str(path), 'SHA256': sha(path),
                        'Original130Index': reverse.get(s['Index']),
                        'ContractZeroTraceIndex': s['Index'] in contract['ZeroTraceIndices']})
    values = np.stack(all_values, axis=2)
    validation = validate_bank(geometry, values)
    with (inputs/'trace-vertices.csv').open() as stream:
        vertices = list(csv.DictReader(stream))
    with (inputs/'trace-triangles.csv').open() as stream:
        triangles = list(csv.DictReader(stream))
    vertex_map = {int(v['vertex']): v for v in vertices}
    max_frame_error = 0.
    for t, actual, val in zip(triangles, geometry, values):
        if int(t['triangle']) != triangle_ids[int(t['triangle'])-1]:
            raise ValueError('Metadata triangle mapping mismatch')
        for j, k in enumerate(('vertex_i', 'vertex_j', 'vertex_k')):
            vertex = vertex_map[int(t[k])]
            expected = np.array([-float(vertex['x']), float(vertex['z']), float(vertex['y'])])
            max_frame_error = max(max_frame_error, float(np.linalg.norm(expected-actual[j])))
            basis = int(vertex['basis'])
            expected_values = np.zeros(135)
            if basis:
                expected_values[basis-1] = 1.
            if not np.array_equal(expected_values, val[j]):
                raise ValueError('Metadata basis index does not match CSV values')
    if len(triangles) != len(geometry) or max_frame_error > 1e-10:
        raise ValueError('Metadata/CSV coordinate frame mismatch')
    validation['MetadataToCSVTransform'] = '(-x,z,y)'
    validation['MetadataToCSVMaxResidual'] = max_frame_error
    production = palace/'palace/models/laplaceoperator.cpp'
    text = production.read_text()
    start = text.index('class TracePotentialCoefficient :')
    end = text.index('\n};', start)+3
    header = text[start:end].replace('private:', 'public:', 1)
    output.mkdir(parents=True, exist_ok=True)
    (output/'production_trace_coefficient.hpp').write_text(
        '// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.\n'
        '// SPDX-License-Identifier: Apache-2.0\n'
        '// Audit extraction from frozen production; only private -> public changed.\n'+header+'\n')
    report = {'InputHashes': {str(p): sha(p) for p in sorted(inputs.glob('*')) if p.is_file()},
              'ProductionCoefficientSHA256': sha(production), 'Sources': mapping,
              'Validation': validation, 'OmittedSourceIndices': [],
              'AuditCoordinates': 'Original mesh coordinate units (micrometers); volt CSV values',
              'CoefficientHeaderSHA256': sha(output/'production_trace_coefficient.hpp')}
    (output/'input-audit.json').write_text(json.dumps(report, indent=2)+'\n')
    (output/'source-files.txt').write_text(''.join(m['CSV']+'\n' for m in mapping))
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path)
    parser.add_argument('palace', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    prepare(args.inputs, args.palace, args.output)
