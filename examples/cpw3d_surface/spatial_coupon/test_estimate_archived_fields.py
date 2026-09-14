#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Bounded native archive regression; outputs only under ARCHIVE_DIAGNOSTIC_SCRATCH.

Set ARCHIVE_DIAGNOSTIC_EXE to a separately built Palace executable. The synthetic
8-hex p4 problem has 729 H1 DOFs and four prescribed fields (affine, kink, zero,
and an explicit nonzero combination). No production data is touched.
"""
import copy
import csv
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import estimate_archived_fields as diagnostic


def box_mesh(path):
    points = [(i / 2, j / 2, k / 2) for k in range(3) for j in range(3) for i in range(3)]
    def idx(i, j, k):
        return i + 3 * j + 9 * k
    cells, faces = [], []
    for k in range(2):
        for j in range(2):
            for i in range(2):
                v = [idx(i, j, k), idx(i+1, j, k), idx(i+1, j+1, k), idx(i, j+1, k),
                     idx(i, j, k+1), idx(i+1, j, k+1), idx(i+1, j+1, k+1), idx(i, j+1, k+1)]
                cells.append((1 if k == 0 else 2, v))
                for active, attr, local in [(k == 0, 3, [0, 3, 2, 1]),
                                             (k == 1, 1, [4, 5, 6, 7]),
                                             (j == 0, 2, [0, 1, 5, 4]),
                                             (i == 1, 2, [1, 2, 6, 5]),
                                             (j == 1, 2, [2, 3, 7, 6]),
                                             (i == 0, 2, [3, 0, 4, 7]),
                                             (k == 0, 4, [4, 5, 6, 7])]:
                    if active:
                        faces.append((attr, [v[n] for n in local]))
    text = 'MFEM mesh v1.0\n\ndimension\n3\n\nelements\n8\n'
    text += ''.join(f'{a} 5 ' + ' '.join(map(str, c)) + '\n' for a, c in cells)
    text += f'\nboundary\n{len(faces)}\n'
    text += ''.join(f'{a} 3 ' + ' '.join(map(str, f)) + '\n' for a, f in faces)
    text += '\nvertices\n27\n3\n' + ''.join(' '.join(map(str, p))+'\n' for p in points)
    path.write_text(text)


def source_triangles():
    triangles = []
    # Boundary rectangles on a source lattice with an off-mesh kink x=.37.
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        u = [0., .37, 1.] if other[0] == 0 else [0., 1.]
        v = [0., .37, 1.] if other[1] == 0 else [0., 1.]
        for side in (0., 1.):
            for a, b in zip(u, u[1:]):
                for c, d in zip(v, v[1:]):
                    corners = []
                    for s, t in [(a, c), (b, c), (b, d), (a, d)]:
                        p = [0., 0., 0.]
                        p[axis], p[other[0]], p[other[1]] = side, s, t
                        corners.append(p)
                    triangles.extend([[corners[i] for i in face] for face in [(0, 1, 2), (0, 2, 3)]])
    return triangles


def fixture(root, polarized=False):
    box_mesh(root / 'box.mesh')
    for index in range(1, 5):
        rows = ['x,y,z,V,triangle']
        for t, triangle in enumerate(source_triangles(), 1):
            for x, y, z in triangle:
                affine = x + (.23*z if polarized else 0.)
                kink = abs(x-.37) + (.4*z if polarized else 0.)
                value = {1: affine, 2: kink, 3: 0., 4: 2*affine-.4*kink}[index]
                rows.append(f'{x},{y},{z},{value:.17g},{t}')
        (root / f'basis-{index:04d}.csv').write_text('\n'.join(rows)+'\n')
    surfaces = [dict(Index=i, Attributes=[attr], Type=kind, Thickness=.002,
                     Permittivity=4., LossTan=.001, EdgeAttributes=[2], EdgeDistances=[.2],
                     LocalizeEdgeEnergy=True, SaveLocalEdgeEnergy=False,
                     EdgeFrameNormal=[0., 0., 1.]) for i, (kind, attr) in enumerate([('MA', 1), ('MS', 3), ('SA', 4)], 1)]
    return {'Problem': {'Type': 'Electrostatic', 'Verbose': 1, 'Output': str(root / 'ordinary'),
                        'OutputFormats': {'Paraview': False, 'GridFunction': False}},
            'Model': {'Mesh': str(root / 'box.mesh'), 'L0': 1.e-6, 'Refinement': {'MaxIts': 0}},
            'Domains': {'Materials': [{'Attributes': [1], 'Permittivity': 2.},
                                      {'Attributes': [2], 'Permittivity': 1.}]},
            'Boundaries': {'PrescribedPotential': [dict(Index=i, Attributes=[1, 2, 3],
                                                        DataFile=str(root / f'basis-{i:04d}.csv')) for i in range(1, 5)],
                           'Postprocessing': {'Dielectric': surfaces}},
            'Solver': {'Order': 4, 'Electrostatic': {'Save': 0, 'ResponseMatrix': True,
                                                    'AggregateResponseMatrix': True},
                       'Linear': {'Type': 'BoomerAMG', 'KSPType': 'CG', 'Tol': 1.e-13,
                                  'MaxIts': 300, 'EstimatorTol': 1.e-12,
                                  'EstimatorMaxIts': 300, 'EstimatorMG': False}}}


def request_for(config):
    sources = config['Boundaries']['PrescribedPotential']
    contract = Path(sources[0]['DataFile']).parent / 'basis-contract.json'
    diagnostic.write(contract, {'Sources': len(sources), 'ZeroTraceIndices': [],
                                'OutputSourceSHA256': {s['DataFile']: diagnostic.digest(s['DataFile']) for s in sources}})
    return {'Version': 1, 'SourceIds': [1, 2, 3, 4],
            'BasisContract': {'Path': str(contract), 'SHA256': diagnostic.digest(contract)},
            'SourceSHA256': [diagnostic.digest(s['DataFile']) for s in config['Boundaries']['PrescribedPotential']],
            'ZeroTraceIndices': [], 'WriteElementIndicators': True,
            'Validation': {'RelResidualTol': 1.e-10, 'AbsResidualTol': 1.e-14,
                           'BCAbsTolV': 1.e-12, 'MinEnergyJ': 1.e-30, 'MinCancellationRatio': 1.e-8},
            'Excitations': [{'Name': f'source{i+1}', 'Coefficients': [int(j == i) for j in range(4)]}
                            for i in range(4)] + [{'Name': 'combination', 'Coefficients': [2., -.4, 0., 0.]},
                                                  {'Name': 'zero', 'Coefficients': [0., 0., 0., 0.]},
                                                  {'Name': 'tiny', 'Coefficients': [1.e-10, 0., 0., 0.]},
                                                  {'Name': 'cancel', 'Coefficients': [2., -.4, 0., -1.]}]}


@unittest.skipUnless(os.getenv('ARCHIVE_DIAGNOSTIC_EXE'), 'requires separately built diagnostic binary')
class ArchiveDiagnosticsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.exe = Path(os.environ['ARCHIVE_DIAGNOSTIC_EXE']).resolve()
        cls.root = Path(tempfile.mkdtemp(prefix='archive-regression-', dir=os.environ['ARCHIVE_DIAGNOSTIC_SCRATCH']))
        cls.config = fixture(cls.root)
        cls.request = request_for(cls.config)
        diagnostic.write(cls.root / 'ordinary.json', cls.config)
        diagnostic.write(cls.root / 'request.json', cls.request)
        cls.env = {k: v for k, v in os.environ.items() if not k.startswith('PALACE_RESPONSE_')}
        cls.env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        cls.archive = cls.root / 'archive'
        cls.run_native(cls.root / 'ordinary.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(cls.archive)}, 'ordinary.log')
        cls.before = {str(p): diagnostic.digest(p) for p in cls.archive.iterdir()}
        cls.result = cls.estimate('estimate', launcher=True)
        cls.rows = {r['Name']: r for r in cls.result['Excitations']}

    @classmethod
    def run_native(cls, config, extra, log, ok=True, ranks=1):
        env = dict(cls.env, **extra)
        with (cls.root / log).open('w') as stream:
            result = subprocess.run(['mpirun', '-n', str(ranks), str(cls.exe), str(config)], env=env,
                                    stdout=stream, stderr=subprocess.STDOUT, timeout=90)
        text = (cls.root / log).read_text()
        if (result.returncode == 0) != ok:
            raise AssertionError(f'unexpected exit {result.returncode}: {text[-5000:]}')
        return text

    @classmethod
    def estimate(cls, name, request=None, config=None, archive=None, ok=True, flags=None,
                 launcher=False, ranks=1):
        config = copy.deepcopy(config or cls.config)
        config['Problem']['Output'] = str(cls.root / name / 'postpro')
        diagnostic.write(cls.root / f'{name}.json', config)
        diagnostic.write(cls.root / f'{name}-request.json', request or cls.request)
        if launcher:
            command = [sys.executable, str(Path(diagnostic.__file__)), str(cls.root / f'{name}.json'),
                       '--request', str(cls.root / f'{name}-request.json'), '--archive', str(archive or cls.archive),
                       '--binary', str(cls.exe), '--binary-sha256', diagnostic.digest(cls.exe),
                       '--output', str(cls.root / name), '--ranks', str(ranks), '--seconds', '90', '--memory-gib', '4']
            with (cls.root / f'{name}-launcher.log').open('w') as stream:
                result = subprocess.run(command, env=cls.env, stdout=stream, stderr=subprocess.STDOUT, timeout=100)
            if (result.returncode == 0) != ok:
                raise AssertionError((cls.root / f'{name}-launcher.log').read_text() + '\n' +
                                     (cls.root / name / 'solver.log').read_text()[-5000:])
            text = ''
        else:
            text = cls.run_native(cls.root / f'{name}.json',
                                  dict(PALACE_RESPONSE_ARCHIVE_DIR=str(archive or cls.archive),
                                       PALACE_RESPONSE_ESTIMATE_ONLY='1',
                                       PALACE_RESPONSE_ESTIMATE_REQUEST=str(cls.root / f'{name}-request.json'),
                                       **(flags or {})), f'{name}.log', ok=ok, ranks=ranks)
        return diagnostic.load(cls.root / name / 'postpro/archive-estimates.json') if ok else text

    def test_analytic_zero_combination_and_ordinary_agreement(self):
        self.assertEqual(self.result['SampleCount'], 8)
        self.assertEqual(self.result['PDESolves'], 0)
        self.assertEqual(self.result['GlobalH1TrueDofs'], 729)
        # U = .5 eps0 <eps_r> * volume * (1 volt / 1 micron)^2; <eps_r>=1.5.
        self.assertAlmostEqual(self.rows['source1']['EnergyJ'] / (.75 * 8.8541878128e-12 * 1.e-6), 1., delta=2.e-9)
        self.assertLess(self.rows['source1']['Eta'], 2.e-10)
        self.assertGreater(self.rows['source2']['Eta'], .001)
        for name in ['source3', 'zero']:
            row = self.rows[name]
            self.assertEqual(row['EnergyJ'], 0.)
            self.assertEqual(row['EtaRawSqrtJ'], 0.)
            self.assertIsNone(row['Eta'])
            self.assertTrue(row['RecoveryConverged'])
            self.assertEqual(row['RecoveryIterations'], 0)
        self.assertIsNone(self.rows['tiny']['Eta'])
        self.assertEqual(self.rows['cancel']['CancellationStatus'], 'strong_cancellation')
        self.assertFalse(self.rows['cancel']['NormalizationUsable'])
        self.assertGreater(self.rows['tiny']['EnergyJ'], 0.)
        for key in ['EnergyJ', 'Eta', 'EtaRawSqrtJ']:
            self.assertAlmostEqual(self.rows['combination'][key] / self.rows['source4'][key], 1., delta=1.e-8)
        for a, b in zip(self.rows['combination']['Interfaces'], self.rows['source4']['Interfaces']):
            self.assertAlmostEqual(a['EnergyJ'] / b['EnergyJ'], 1., delta=1.e-8)
            self.assertAlmostEqual(a['Windows'][0]['InsideJ'] / b['Windows'][0]['InsideJ'], 1., delta=1.e-8)
        with (self.root / 'ordinary/error-indicators.csv').open() as stream:
            ordinary = list(csv.reader(stream))
        # Historical diagnostic averages normalized per-source squared norms, including zero.
        rms = math.sqrt(sum((self.rows[f'source{i}']['Eta'] or 0.)**2 for i in range(1, 5)) / 4)
        self.assertAlmostEqual(rms / float(ordinary[1][0]), 1., delta=1.e-9)
        self.assertTrue(all(c['Passed'] for c in self.result['Checks']))
        self.assertTrue(all(r['RecoveryConverged'] for r in self.rows.values()))
        # Quadratic form against independently assembled ordinary domain response matrix.
        with (self.root / 'ordinary/domain-response-matrix.csv').open() as stream:
            entries = list(csv.reader(stream))[1:]
        c = [2., -.4, 0., 0.]
        q = sum(float(row[2]) * c[int(float(row[0]))-1] * c[int(float(row[1]))-1] * (1 if row[0] == row[1] else 2)
                for row in entries)
        self.assertAlmostEqual(q / self.rows['combination']['EnergyJ'], 1., delta=1.e-9)
        # Same quadratic forms for all three whole and window interface energies.
        with (self.root / 'ordinary/surface-response-matrix.csv').open() as stream:
            entries = list(csv.reader(stream))[1:]
        for interface in self.rows['combination']['Interfaces']:
            sums = [0., 0.]
            for row in entries:
                if int(float(row[0])) != interface['Index']:
                    continue
                i, j = int(float(row[3]))-1, int(float(row[4]))-1
                factor = c[i] * c[j] * (1 if i == j else 2)
                sums[0] += factor * float(row[5])
                sums[1] += factor * float(row[8])
            self.assertAlmostEqual(sums[0] / interface['Windows'][0]['InsideJ'], 1., delta=1.e-9)
            self.assertAlmostEqual(sums[1] / interface['EnergyJ'], 1., delta=1.e-9)
        self.assert_all_quadratic_responses(self.result, self.request, self.root / 'ordinary')
        for n, row in enumerate(self.result['Excitations']):
            with (self.root / 'estimate/postpro' / f'archive-elements-case-{n:04d}-rank-000000.csv').open() as stream:
                local = list(csv.DictReader(stream))
            self.assertEqual(len(local), 8)
            self.assertEqual(len({tuple(float(r[f'center_{axis}_m']) for axis in 'xyz') for r in local}), 8)
            self.assertAlmostEqual(sum(float(r['eta_raw_squared_J']) for r in local), row['EtaRawSqrtJ']**2,
                                   delta=max(1.e-60, row['EtaRawSqrtJ']**2 * 1.e-12))

    def assert_all_quadratic_responses(self, result, request, ordinary):
        with (ordinary / 'domain-response-matrix.csv').open() as stream:
            domain = list(csv.reader(stream))[1:]
        with (ordinary / 'surface-response-matrix.csv').open() as stream:
            surface = list(csv.reader(stream))[1:]
        components = ['InsideJ', 'InsideNormalJ', 'InsideTangentialJ',
                      'TotalJ', 'TotalNormalJ', 'TotalTangentialJ']
        coefficients = {r['Name']: r['Coefficients'] for r in request['Excitations']}
        for excitation in result['Excitations']:
            c = coefficients[excitation['Name']]
            def quadratic(rows, ii, jj, column):
                terms = []
                for r in rows:
                    i, j = int(float(r[ii]))-1, int(float(r[jj]))-1
                    terms.append(c[i]*c[j]*(1 if i == j else 2)*float(r[column]))
                return sum(terms), max(1.e-60, sum(abs(t) for t in terms)*1.e-9)
            q, tol = quadratic(domain, 0, 1, 2)
            self.assertAlmostEqual(excitation['EnergyJ'], q, delta=tol)
            self.assertEqual(len(excitation['InterfaceResponses']), 3)
            for response in excitation['InterfaceResponses']:
                rows = [r for r in surface if int(float(r[0])) == response['Index']]
                self.assertTrue(rows)
                for column, name in enumerate(components, 5):
                    q, tol = quadratic(rows, 3, 4, column)
                    self.assertAlmostEqual(response[name], q, delta=tol,
                                           msg=f"{excitation['Name']}:{response['Index']}:{name}")

    def test_nontrivial_polarizations(self):
        for recovered in [False, True]:
            root = self.root / f'polarized-input-{recovered}'
            root.mkdir()
            config = fixture(root, polarized=True)
            for interface in config['Boundaries']['Postprocessing']['Dielectric']:
                interface['FluxRecovery'] = recovered
            request = request_for(config)
            diagnostic.write(root / 'ordinary.json', config)
            archive = root / 'archive'
            self.run_native(root / 'ordinary.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive)},
                            f'polarized-ordinary-{recovered}.log')
            result = self.estimate(f'polarized-{recovered}', request=request, config=config, archive=archive)
            self.assert_all_quadratic_responses(result, request, root / 'ordinary')
            first = result['Excitations'][0]['InterfaceResponses']
            self.assertGreater(first[0]['TotalNormalJ'], 1.e-25)
            self.assertGreater(first[1]['TotalNormalJ'], 1.e-25)
            self.assertGreater(first[2]['TotalNormalJ'], 1.e-25)
            self.assertGreater(first[2]['TotalTangentialJ'], 1.e-25)
            self.assertEqual(result['Excitations'][0]['SurfaceUsesRecoveredFlux'], recovered)

    def test_read_only_export_and_optional_inputs_rejected(self):
        sentinel = self.root / 'box_preprocessed.mesh'
        sentinel.write_bytes(b'preserved neighboring mesh\n')
        original = sentinel.read_bytes()
        for option in ['ExportPrerefinedMesh', 'Partitioning', 'OwnershipDataFile']:
            for launcher in [False, True]:
                config = copy.deepcopy(self.config)
                if option == 'OwnershipDataFile':
                    config['Boundaries']['Postprocessing']['Dielectric'][0][option] = str(sentinel)
                else:
                    config['Model'][option] = True if option == 'ExportPrerefinedMesh' else str(sentinel)
                name = f'reject-{option}-{launcher}'
                text = self.estimate(name, config=config, launcher=launcher, ok=False)
                if launcher:
                    text = (self.root / f'{name}-launcher.log').read_text()
                self.assertIn(option, text)
                self.assertEqual(sentinel.read_bytes(), original)
                self.assertFalse((self.root / name).exists())

    def test_ground_precedence_and_ground_only_corruption(self):
        root = self.root / 'grounded-input'
        root.mkdir()
        config = fixture(root)
        config['Boundaries']['Ground'] = {'Attributes': [3]}
        for source in config['Boundaries']['PrescribedPotential']:
            source['Attributes'] = [1, 2]
        request = request_for(config)
        request['Excitations'] = request['Excitations'][:1]
        diagnostic.write(root / 'ordinary.json', config)
        for ranks in [1, 2]:
            archive = root / f'archive{ranks}'
            self.run_native(root / 'ordinary.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive)},
                            f'grounded-ordinary{ranks}.log', ranks=ranks)
            result = self.estimate(f'grounded{ranks}', request=request, config=config,
                                   archive=archive, ranks=ranks)
            categories = result['BoundaryCategories']
            self.assertEqual(categories['MatchingOnlyTrueDofs'], 305)
            self.assertEqual(categories['PhysicalGroundTrueDofs'], 81)
            self.assertEqual(categories['IntersectionTrueDofs'], 32)
            for check in result['Checks']:
                self.assertEqual(check['MatchingOnlyBCMaxV'], 0.)
                self.assertEqual(check['PhysicalGroundBCMaxV'], 0.)
                self.assertEqual(check['IntersectionBCMaxV'], 0.)
        bad = root / 'bad-archive'
        shutil.copytree(root / 'archive1', bad)
        path = bad / 'source-000001-rank-000000-V.bin'
        data = bytearray(path.read_bytes())
        # Vertex 1 is (.5,0,0): intended x=.5 V, imposed ground=0 at intersection.
        self.assertEqual(struct.unpack_from('=d', data, 48+8)[0], 0.)
        struct.pack_into('=d', data, 48+8, .125/result['VoltageScaleV'])
        path.write_bytes(data)
        self.estimate('ground-corrupt', request=request, config=config, archive=bad, ok=False)
        failed = diagnostic.load(self.root / 'ground-corrupt/postpro/archive-estimates.json')['Checks'][0]
        self.assertEqual(failed['MatchingOnlyBCMaxV'], 0.)
        self.assertAlmostEqual(failed['PhysicalGroundBCMaxV'], .125)
        self.assertAlmostEqual(failed['IntersectionBCMaxV'], .125)
        self.assertFalse(failed['Passed'])

    def test_unused_source_assignment_swap_rejected_before_launch(self):
        config = copy.deepcopy(self.config)
        request = copy.deepcopy(self.request)
        request['Excitations'] = request['Excitations'][:1]
        sources = config['Boundaries']['PrescribedPotential']
        sources[1]['DataFile'], sources[2]['DataFile'] = sources[2]['DataFile'], sources[1]['DataFile']
        request['SourceSHA256'] = [diagnostic.digest(s['DataFile']) for s in sources]
        self.assertEqual(diagnostic.validate_request(request, sources), {1})
        self.estimate('unused-swap', request=request, config=config, launcher=True, ok=False)
        self.assertIn('ordered source index/name/hash mapping',
                      (self.root / 'unused-swap-launcher.log').read_text())
        self.assertFalse((self.root / 'unused-swap').exists())

    def test_surface_sweep_preserves_field_residual_recovery_and_primary_observables(self):
        for number, extras in enumerate([[0, 1.5], [0, 4, 4], [0, 13]]):
            bad = copy.deepcopy(self.request)
            bad['Excitations'] = bad['Excitations'][:1]
            bad['SurfaceQuadratureExtras'] = extras
            self.assertIn('SurfaceQuadratureExtras', self.estimate(f'bad-surface-order{number}', request=bad, ok=False))
        request = copy.deepcopy(self.request)
        request['SurfaceQuadratureExtras'] = [0, 4, 8, 12]
        result = self.estimate('surface-sweep', request=request)
        self.assertEqual(result['Checks'], self.result['Checks'])
        self.assertEqual(result['Recovery'], self.result['Recovery'])
        for before, after in zip(self.result['Excitations'], result['Excitations']):
            primary = {k: v for k, v in after.items() if k not in ['SurfaceQuadratureSweeps', 'SurfaceQuadratureControls']}
            self.assertEqual(primary, before)
            self.assertEqual(after['SurfaceQuadratureSweeps'][0]['InterfaceResponses'], after['InterfaceResponses'])
            for sweep in after['SurfaceQuadratureSweeps']:
                for old, new in zip(before['InterfaceResponses'], sweep['InterfaceResponses']):
                    for key in ['TotalJ', 'TotalNormalJ', 'TotalTangentialJ']:
                        self.assertAlmostEqual(old[key], new[key], delta=1.e-35+1.e-9*abs(old[key]))
        rules = diagnostic.load(self.root / 'surface-sweep/postpro/archive-quadrature-case-0000-rank-000000.json')
        self.assertEqual({r['ExtraOrder'] for r in rules}, {0, 4, 8, 12})
        for rule in rules:
            self.assertGreater(rule['MinimumReferenceWeight'], 0.)
            self.assertGreater(rule['PointCount'], 0)
            self.assertEqual(rule['RequestedOrder'], 8+rule['ExtraOrder'])
            self.assertGreaterEqual(rule['ActualRuleOrder'], rule['RequestedOrder'])

    @unittest.skipUnless(os.getenv('ARCHIVE_DIAGNOSTIC_BASELINE_EXE'), 'requires frozen baseline binary')
    def test_surface_baseline_against_frozen_binary(self):
        new_exe = type(self).exe
        try:
            type(self).exe = Path(os.environ['ARCHIVE_DIAGNOSTIC_BASELINE_EXE'])
            baseline = self.estimate('frozen-binary-baseline')
        finally:
            type(self).exe = new_exe
        self.assertEqual(baseline['Checks'], self.result['Checks'])
        for before, after in zip(baseline['Excitations'], self.result['Excitations']):
            for key in ['EnergyJ', 'EtaRawSqrtJ', 'NormalizationSqrtJ', 'RecoveryIterations', 'RecoveryConverged']:
                self.assertEqual(before[key], after[key])
            for a, b in zip(before['InterfaceResponses'], after['InterfaceResponses']):
                for key in ['InsideJ', 'InsideNormalJ', 'InsideTangentialJ', 'TotalJ', 'TotalNormalJ', 'TotalTangentialJ']:
                    self.assertAlmostEqual(a[key], b[key], delta=1.e-35+1.e-10*abs(a[key]))

    def test_surface_polynomial_and_independent_split_window(self):
        root = self.root / 'surface-polynomial-input'
        root.mkdir()
        config = fixture(root, polarized=True)
        for material in config['Domains']['Materials']:
            material['Permittivity'] = 1.
        request = request_for(config)
        request['Excitations'] = request['Excitations'][:1]
        request['SurfaceQuadratureExtras'] = [0, 4, 8, 12]
        diagnostic.write(root / 'ordinary.json', config)
        archive = root / 'archive'
        self.run_native(root / 'ordinary.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive)}, 'surface-polynomial-ordinary.log')
        result = self.estimate('surface-polynomial', request=request, config=config, archive=archive)
        source = result['Excitations'][0]
        expected_total = .5 * 8.8541878128e-12 * .002e-6 * .23**2 / 4.
        evidence = []
        for sweep in source['SurfaceQuadratureSweeps']:
            ma = next(r for r in sweep['InterfaceResponses'] if r['Index'] == 1)
            self.assertAlmostEqual(ma['TotalNormalJ']/expected_total, 1., delta=2.e-9)
            # Independent GL implementation, and a separate exact geometric split:
            # a .2-wide perimeter window has area1-(1-2*.2)^2=.64 on the unit square.
            n = (8+sweep['ExtraOrder'])//2+1
            outside = 0.
            for i in range(n):
                x = math.cos(math.pi*(i+.75)/(n+.5))
                for iteration in range(100):
                    a, b = 1., x
                    for k in range(2, n+1):
                        a, b = b, ((2*k-1)*x*b-(k-1)*a)/k
                    derivative = n*(x*b-a)/(x*x-1)
                    step = b/derivative
                    x -= step
                    if abs(step) < 1.e-15:
                        break
                else:
                    self.fail('independent GL root iteration failed')
                weight = 2/((1-x*x)*derivative*derivative)
                for origin in [0., .5]:
                    coordinate = origin+.25*(x+1)
                    if .2 < coordinate < .8:
                        outside += .25*weight
            expected_rule_fraction = 1-outside**2
            measured = ma['InsideNormalJ']/ma['TotalNormalJ']
            self.assertAlmostEqual(measured, expected_rule_fraction, delta=1.e-9)
            evidence.append({'ExtraOrder': sweep['ExtraOrder'], 'MeasuredWindowFraction': measured,
                             'IndependentRuleFraction': expected_rule_fraction, 'ExactSplitFraction': .64,
                             'WindowQuadratureErrorFraction': measured-.64})
        self.assertGreater(max(abs(r['WindowQuadratureErrorFraction']) for r in evidence), .001)
        diagnostic.write(self.root / 'independent-split-window.json', evidence)

    def test_triangle_rule_metadata_and_negative_rule_rejection(self):
        for order in [4, 8]:
            root = self.root / f'triangle-quadrature-input-p{order}'
            root.mkdir()
            config = fixture(root, polarized=True)
            config['Model']['MakeSimplex'] = True
            config['Solver']['Order'] = order
            for material in config['Domains']['Materials']:
                material['Permittivity'] = 1.
            request = request_for(config)
            request['Excitations'] = request['Excitations'][:1]
            request['SurfaceQuadratureExtras'] = [0, 4, 8, 12]
            diagnostic.write(root / 'ordinary.json', config)
            archive = root / 'archive'
            self.run_native(root / 'ordinary.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive)}, f'triangle-ordinary-p{order}.log')
            name = f'triangle-quadrature-p{order}'
            result = self.estimate(name, request=request, config=config, archive=archive, ok=order==4)
            if order == 8:
                self.assertIn('Nonpositive surface matrix quadrature weight', result)
                self.assertIn('no rule substitution', result)
            else:
                rules = diagnostic.load(self.root / name / 'postpro/archive-quadrature-case-0000-rank-000000.json')
                self.assertTrue(rules)
                self.assertTrue(all(r['Geometry'] == 2 and r['MinimumReferenceWeight'] > 0 for r in rules))

    def test_recovery_tolerance_convergence(self):
        values = []
        for name, tol, its in [('loose', .5, 5), ('medium', 1.e-7, 300), ('tight', 1.e-11, 300), ('unconverged', 1.e-14, 1)]:
            config = copy.deepcopy(self.config)
            config['Solver']['Linear'].update(EstimatorTol=tol, EstimatorMaxIts=its)
            request = copy.deepcopy(self.request)
            request['Excitations'] = [request['Excitations'][1]]
            row = self.estimate(name, request, config)['Excitations'][0]
            values.append(row)
        self.assertLess(abs(values[1]['Eta'] / values[2]['Eta'] - 1), 1.e-6)
        self.assertGreater(abs(values[0]['Eta'] / values[2]['Eta'] - 1), .001)
        self.assertFalse(values[3]['RecoveryConverged'])
        self.assertEqual(values[3]['Status'], 'recovery_unconverged')

    def test_malformed_and_permuted_archives_rejected(self):
        original = (self.archive / 'source-000001-rank-000000-V.bin').read_bytes()
        cases = {'truncated': original[:-1], 'trailing': original + b'0',
                 'header': b'badmagic' + original[8:]}
        for key, offset, value in [('source', 16, 2), ('rank', 24, 1), ('ranks', 32, 2), ('dofs', 40, 728),
                                    ('huge', 40, 2**40)]:
            data = bytearray(original)
            struct.pack_into('=q', data, offset, value)
            cases[key] = bytes(data)
        data = bytearray(original)
        struct.pack_into('=d', data, 48, float('nan'))
        cases['nonfinite'] = bytes(data)
        struct.pack_into('=d', data, 48, float('inf'))
        cases['infinite'] = bytes(data)
        values = list(struct.unpack(f'={int((len(original)-48)/8)}d', original[48:]))
        cases['permuted'] = original[:48] + struct.pack(f'={len(values)}d', *values[::-1])
        # Swap two interior DOFs, leaving every prescribed boundary value untouched.
        data = bytearray(original)
        # Hex H1 element interior true DOFs are last; two interior values differ in x.
        a, b = len(values)-27, len(values)-25
        struct.pack_into('=d', data, 48+8*a, values[b])
        struct.pack_into('=d', data, 48+8*b, values[a])
        self.assertNotEqual(values[a], values[b])
        cases['interior-permuted'] = bytes(data)
        request = copy.deepcopy(self.request)
        request['Excitations'] = [request['Excitations'][0]]
        for name, content in cases.items():
            with self.subTest(name=name):
                archive = self.root / f'bad-{name}'
                archive.mkdir()
                (archive / 'source-000001-rank-000000-V.bin').write_bytes(content)
                text = self.estimate(f'reject-{name}', request, archive=archive, ok=False)
                self.assertIn('Verification failed', text)
                if name == 'interior-permuted':
                    report = diagnostic.load(self.root / f'reject-{name}/postpro/archive-estimates.json')
                    self.assertEqual(report['Checks'][0]['BCMaxV'], 0.)
                    self.assertGreater(report['Checks'][0]['ResidualNorm'], report['Checks'][0]['ResidualThreshold'])

    def test_incompatible_flags_output_reuse_constraints_and_nonmutation(self):
        for n, flag in enumerate(diagnostic.CONFLICTS):
            text = self.estimate(f'flag{n}', flags={flag: '1'}, ok=False)
            self.assertIn('conflicts', text)
            self.assertFalse((self.root / f'flag{n}/postpro').exists())
        config = copy.deepcopy(self.config)
        config['Model']['Refinement']['MaxIts'] = 1
        self.assertIn('no AMR', self.estimate('amr', config=config, ok=False))
        config = copy.deepcopy(self.config)
        config['Solver']['Order'] = 3
        self.assertIn('archive header', self.estimate('wrong-order', config=config, ok=False))
        request = copy.deepcopy(self.request)
        request['ZeroTraceIndices'] = [1]
        self.assertIn('constrained', self.estimate('constrained', request, ok=False))
        request = copy.deepcopy(self.request)
        request['SourceIds'] = [2, 1, 3, 4]
        self.assertIn('SourceIds', self.estimate('source-order', request, ok=False))
        # Direct executable reuse is rejected before touching any ordinary output.
        env = dict(PALACE_RESPONSE_ARCHIVE_DIR=str(self.archive), PALACE_RESPONSE_ESTIMATE_ONLY='1',
                   PALACE_RESPONSE_ESTIMATE_REQUEST=str(self.root / 'request.json'))
        text = self.run_native(self.root / 'ordinary.json', env, 'reuse.log', ok=False)
        self.assertIn('must not already exist', text)
        config = copy.deepcopy(self.config)
        config['Problem']['Output'] = str(self.archive)
        diagnostic.write(self.root / 'overlap.json', config)
        self.assertIn('overlaps input', self.run_native(self.root / 'overlap.json', env, 'overlap.log', ok=False))
        self.assertEqual(self.before, {str(p): diagnostic.digest(p) for p in self.archive.iterdir()})
        self.assertEqual(diagnostic.load(self.root / 'estimate/provenance-after.json')['ChangedInputs'], [])

    def test_selected_source_with_multigrid_and_no_element_output(self):
        config = copy.deepcopy(self.config)
        config['Solver']['Linear']['EstimatorMG'] = True
        request = copy.deepcopy(self.request)
        request['Excitations'] = request['Excitations'][1:2]
        request['WriteElementIndicators'] = False
        archive = self.root / 'selected-archive'
        archive.mkdir()
        source = self.archive / 'source-000002-rank-000000-V.bin'
        shutil.copy2(source, archive / source.name)
        result = self.estimate('selected-mg', request=request, config=config, archive=archive)
        row = result['Excitations'][0]
        self.assertEqual(result['SampleCount'], 1)
        self.assertTrue(row['RecoveryConverged'])
        self.assertAlmostEqual(row['Eta'] / self.rows['source2']['Eta'], 1., delta=1.e-9)
        self.assertFalse(list((self.root / 'selected-mg/postpro').glob('archive-elements-*')))

    def test_streaming_and_reducer_paths_unchanged(self):
        config = copy.deepcopy(self.config)
        config['Problem']['Output'] = str(self.root / 'streaming')
        diagnostic.write(self.root / 'streaming.json', config)
        archive = self.root / 'streaming-archive'
        self.run_native(self.root / 'streaming.json',
                        {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive), 'PALACE_RESPONSE_ARCHIVE_ONLY': '1'},
                        'streaming.log')
        for source in self.archive.iterdir():
            self.assertEqual(diagnostic.digest(source), diagnostic.digest(archive / source.name))
        config['Problem']['Output'] = str(self.root / 'reducer')
        diagnostic.write(self.root / 'reducer.json', config)
        self.run_native(self.root / 'reducer.json',
                        {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive), 'PALACE_RESPONSE_REDUCE_ONLY': '1',
                         'PALACE_RESPONSE_BLOCK_SIZE': '2'}, 'reducer.log')
        for filename in ['domain-response-matrix.csv', 'surface-response-matrix.csv']:
            with (self.root / 'ordinary' / filename).open() as stream:
                a = sorted([list(map(float, r)) for r in list(csv.reader(stream))[1:]])
            with (self.root / 'reducer' / filename).open() as stream:
                b = sorted([list(map(float, r)) for r in list(csv.reader(stream))[1:]])
            self.assertEqual(len(a), len(b))
            for x, y in zip(a, b):
                for xx, yy in zip(x, y):
                    self.assertAlmostEqual(xx, yy, delta=max(1.e-40, abs(xx)*1.e-9))

    def test_two_rank_layout_and_archive_agreement(self):
        config = copy.deepcopy(self.config)
        config['Problem']['Output'] = str(self.root / 'ordinary2')
        diagnostic.write(self.root / 'ordinary2.json', config)
        archive = self.root / 'archive2'
        self.run_native(self.root / 'ordinary2.json', {'PALACE_RESPONSE_ARCHIVE_DIR': str(archive)}, 'ordinary2.log', ranks=2)
        result = self.estimate('estimate2', archive=archive, ranks=2)
        self.assertEqual(result['Ranks'], 2)
        # A bad shard on just one rank must fail promptly, not hang in a collective.
        bad = self.root / 'archive2-bad'
        shutil.copytree(archive, bad)
        path = bad / 'source-000001-rank-000001-V.bin'
        path.write_bytes(path.read_bytes()[:-1])
        request = copy.deepcopy(self.request)
        request['Excitations'] = request['Excitations'][:1]
        self.estimate('reject-one-rank', request=request, archive=bad, ok=False, ranks=2)
        for row in result['Excitations']:
            reference = self.rows[row['Name']]
            if reference['EnergyJ'] and reference['NormalizationUsable']:
                self.assertAlmostEqual(row['EnergyJ'] / reference['EnergyJ'], 1., delta=1.e-9)
            if reference['Eta'] and reference['Eta'] > 1.e-8:
                self.assertAlmostEqual(row['Eta'] / reference['Eta'], 1., delta=1.e-8)


class RequestValidationTest(unittest.TestCase):
    def test_hash_constraints_and_nonfinite(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('ARCHIVE_DIAGNOSTIC_SCRATCH')) as temp:
            config = fixture(Path(temp))
            request = request_for(config)
            sources = config['Boundaries']['PrescribedPotential']
            self.assertEqual(diagnostic.validate_request(request, sources), {1, 2, 3, 4})
            diagnostic.validate_basis_contract(request, sources)
            bad = copy.deepcopy(request)
            bad['ZeroTraceIndices'] = [3]
            with self.assertRaisesRegex(ValueError, 'fixed ZeroTraceIndices'):
                diagnostic.validate_basis_contract(bad, sources)
            for extras in [[], [4], [0, -1], [0, 13], [0, 4, 4], [0, 1.5]]:
                bad = copy.deepcopy(request)
                bad['SurfaceQuadratureExtras'] = extras
                with self.assertRaisesRegex(ValueError, 'SurfaceQuadratureExtras'):
                    diagnostic.validate_request(bad, sources)
            for mutate in [lambda r: r['SourceSHA256'].__setitem__(0, '0'*64),
                           lambda r: r.update(ZeroTraceIndices=[1]),
                           lambda r: r['Excitations'][0]['Coefficients'].__setitem__(0, float('inf'))]:
                bad = copy.deepcopy(request)
                mutate(bad)
                with self.assertRaises(ValueError):
                    diagnostic.validate_request(bad, sources)


if __name__ == '__main__':
    unittest.main()
