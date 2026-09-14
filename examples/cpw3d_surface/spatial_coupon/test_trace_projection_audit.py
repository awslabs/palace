#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bounded synthetic tests. Set TRACE_PROJECTION_EXE and TRACE_PROJECTION_SCRATCH.

No build/install implicit in tests; all files go under the supplied audit scratch.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from numpy.polynomial import Polynomial

from prepare_trace_projection_audit import read_trace, validate_bank
from summarize_trace_projection_audit import fixed_groups, rank_diagnostic


TRIANGLES = np.array([[[0., 0, 1], [1, 0, 1], [1, 1, 1]],
                      [[0., 0, 1], [1, 1, 1], [0, 1, 1]]])


class TraceProjectionDataTests(unittest.TestCase):
    def test_affine_and_continuous_kink_data(self):
        for value in (TRIANGLES[:, :, 0]+2*TRIANGLES[:, :, 1],
                      abs(TRIANGLES[:, :, 0]-TRIANGLES[:, :, 1])):
            self.assertEqual(validate_bank(TRIANGLES, value[:, :, None])['PositiveAreaOverlaps'], 0)

    def test_overlapping_interiors(self):
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            validate_bank(np.array([TRIANGLES[0], TRIANGLES[0]]), np.zeros((2, 3, 1)))

    def test_inconsistent_shared_vertices(self):
        values = np.zeros((2, 3, 1))
        values[1, 0, 0] = 1
        with self.assertRaisesRegex(ValueError, 'Inconsistent'):
            validate_bank(TRIANGLES, values)

    def test_degenerate_triangle(self):
        with self.assertRaisesRegex(ValueError, 'Degenerate'):
            validate_bank(np.zeros((1, 3, 3)), np.zeros((1, 3, 1)))

    def test_fixed_model_mapping_not_numerical_filtering(self):
        constrained = list(range(1, 41))
        contract = {'Model': 'test', 'ZeroTraceIndices': constrained}
        library = {'Models': [{'Name': 'test', 'ZeroTraceIndices': constrained}]}
        groups = fixed_groups(contract, library)
        self.assertEqual(groups['all135'], list(range(1, 136)))
        self.assertEqual(groups['free95'], list(range(41, 136)))
        library['Models'][0]['ZeroTraceIndices'] = list(range(2, 42))
        with self.assertRaisesRegex(ValueError, 'mapping mismatch'):
            fixed_groups(contract, library)

    def test_gram_rank_retains_zero_columns(self):
        result = rank_diagnostic(np.diag([4., 0., 9.]))
        self.assertEqual(result['ZeroDiagonalColumns'], 1)
        self.assertEqual(len(result['NormalizedGramEigenvalues']), 3)
        self.assertEqual(set(result['RanksByRelativeGramEigenvalueThreshold'].values()), {2})

    def test_malformed_csv(self):
        root = os.environ.get('TRACE_PROJECTION_SCRATCH')
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            path = Path(tmp)/'malformed.csv'
            for text in ('x,y,z,V,triangle\n0,0,0,1,1\n',
                         'x,y,z,V,triangle\n0,0,nan,1,1\n',
                         'x,y,z,V,triangle\n0,0,0,1,0.5\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_trace(path)


@unittest.skipUnless(os.environ.get('TRACE_PROJECTION_EXE') and
                     os.environ.get('TRACE_PROJECTION_SCRATCH'),
                     'TRACE_PROJECTION_EXE and TRACE_PROJECTION_SCRATCH are required')
class TraceProjectionMFEMTests(unittest.TestCase):
    def setUp(self):
        self.exe = os.environ['TRACE_PROJECTION_EXE']
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ['TRACE_PROJECTION_SCRATCH'])
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def run_case(self, function, aligned=False, offset=0):
        trace = self.root/'source.csv'
        with trace.open('w') as stream:
            writer = csv.writer(stream)
            writer.writerow(['x', 'y', 'z', 'V', 'triangle'])
            for i, triangle in enumerate(TRIANGLES):
                for p in triangle:
                    writer.writerow([p[0], p[1], p[2]+offset, function(*p), i+1])
        manifest = self.root/'manifest.txt'
        manifest.write_text(str(trace)+'\n')
        mesh = self.root/'square.msh'
        # Every open edge in this synthetic square is explicitly a physical ground.
        elements = ['1 2 0 1 2', '1 2 0 2 3'] if aligned else ['1 3 0 1 2 3']
        # MFEM native format keeps embedding dimension 3 for a single flat face;
        # the Gmsh reader reduces constant-z surfaces to dimension 2.
        mesh.write_text('MFEM mesh v1.0\n\ndimension\n2\n\nelements\n'
                        f'{len(elements)}\n'+'\n'.join(elements)+'\n\nboundary\n4\n'
                        '5001 1 0 1\n6001 1 1 2\n6001 1 2 3\n6001 1 3 0\n'
                        '\nvertices\n4\n3\n0 0 1\n1 0 1\n1 1 1\n0 1 1\n')
        prefix = self.root/'result'
        result = subprocess.run([self.exe, str(mesh), str(manifest), str(prefix)],
                                capture_output=True, text=True, timeout=15)
        if offset:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Source coverage/shared-node', result.stderr)
            return
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        with (self.root/'result-norms.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        return {int(r['quadrature_order']): {k: float(v) for k, v in r.items()} for r in rows}

    def test_globally_affine_reproduced_roundoff(self):
        for aligned in (False, True):
            data = self.run_case(lambda x, y, z: 1+2*x-3*y, aligned)
            self.assertLess(data[20]['raw_interpolation_error_l2'], 1e-13)

    def test_aligned_p1_kink(self):
        data = self.run_case(lambda x, y, z: abs(x-y), True)
        self.assertLess(data[20]['raw_interpolation_error_l2'], 1e-13)
        self.assertEqual(data[20]['kink_crossed_faces'], 0)

    def test_nonaligned_p1_kink(self):
        data = self.run_case(lambda x, y, z: abs(x-y))
        self.assertGreater(data[20]['raw_interpolation_error_l2'], 1e-3)
        self.assertEqual(data[20]['kink_crossed_faces'], 1)
        self.assertAlmostEqual(data[16]['raw_interpolation_error_l2'],
                               data[20]['raw_interpolation_error_l2'], places=12)

    def test_physical_ground_contact_override(self):
        data = self.run_case(lambda x, y, z: 1.)
        self.assertLess(data[20]['raw_interpolation_error_l2'], 1e-13)
        self.assertEqual(data[20]['maximum_ground_override'], 1.)
        self.assertEqual(data[20]['maximum_retained_nodal_value'], 1.)
        # Independent tensor-product p4 GLL polynomial: zero at endpoints, one at
        # all interior nodes. The 2D projection is g(x)g(y), with corners zero.
        nodes = np.array([0., (1-np.sqrt(3/7))/2, .5, (1+np.sqrt(3/7))/2, 1.])
        g = Polynomial.fit(nodes, [0, 1, 1, 1, 0], 4).convert()
        integral = g.integ()(1)-g.integ()(0)
        squared = (g*g).integ()(1)-(g*g).integ()(0)
        expected = np.sqrt(1-2*integral**2+squared**2)
        self.assertAlmostEqual(data[20]['grounded_error_l2'], expected, places=12)
        self.assertAlmostEqual(data[20]['projected_l2'], squared, places=12)
        with (self.root/'result-nodal-stats.csv').open() as stream:
            nodal = next(csv.DictReader(stream))
        self.assertAlmostEqual(float(nodal['raw_contact_l2_split']), 2., places=12)
        self.assertEqual(float(nodal['raw_contact_split_max']), 1.)
        self.assertEqual(int(nodal['raw_nonzero_dofs']), 25)
        self.assertEqual(int(nodal['grounded_nonzero_dofs']), 9)

    def test_genuinely_off_surface_fails_closed(self):
        self.run_case(lambda x, y, z: 1., offset=1e-5)


if __name__ == '__main__':
    unittest.main()
