# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('matched_benchmark',Path(__file__).with_name('run_matched_mesh_benchmark.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class MatchedMeshBenchmarkTest(unittest.TestCase):
    def test_source_and_hierarchical_setup_timing(self):
        with tempfile.TemporaryDirectory() as tmp:
            log=Path(tmp)/'worker.log'
            log.write_text(' H1 (p = 5): 12345, ND (p = 5): 45678\n'
                'Response source timing: index=2, iterations=7, solve_seconds=1.200000000e+00, total_seconds=1.300000000e+00\n'
                'Elapsed Time Report (s) Min. Max. Avg.\n'
                'Operator Construction 1.0 2.0 1.5\n'
                'Linear Solve 4.0 6.0 5.0\n'
                '  Setup 2.0 3.0 2.5\n'
                'Total 10.0 11.0 10.5\n'
                'Peak Memory Per-Node Total\n')
            result=module.timing(log)
            self.assertEqual(result['H1Dofs'],12345)
            self.assertEqual(result['Sources'][0]['Iterations'],7)
            self.assertEqual(result['Sources'][0]['TotalSeconds'],1.3)
            self.assertEqual(result['ElapsedTimerSeconds']['Linear Solve/Setup']['Max'],3.0)

    def test_upper_triangle_is_reconstructed_for_matrix_norms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'matrix.csv'
            path.write_text('basis_i,basis_j,Q_ij (J)\n1,1,2\n1,2,1\n2,2,3\n')
            values=module.matrix(path)[0]
            self.assertEqual(values[(2,1)],1.)
            self.assertAlmostEqual(module.norm(values)**2,15.)

    def test_matrix_comparison_retains_defect_and_excludes_thin_spr_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);manifest={'Cases':[]};records={}
            for kind in ('thin','fabricated'):
                for strategy in ('prism','tet'):
                    name=f'{kind}-{strategy}';path=root/'postpro'/f'{name}-reducer';path.mkdir(parents=True)
                    value=(1.0 if kind=='thin' else .9)+(.01 if strategy=='tet' else 0)
                    (path/'domain-response-matrix.csv').write_text(f'basis_i,basis_j,Q_ij (J)\n1,1,{value}\n')
                    (path/'surface-response-matrix.csv').write_text('interface,basis_i,basis_j,Q_total_ij (J)\n1,1,1,0.1\n')
                    manifest['Cases'].append(dict(Name=name,Kind=kind,Strategy=strategy))
                    records[name]=dict(Completed=True,TotalSeconds=10 if strategy=='prism' else 5,
                                       Worker=dict(H1Dofs=100 if strategy=='prism' else 20))
            result=module.compare(root,manifest,records)
            self.assertFalse(result['ThinRawSPRGated'])
            self.assertFalse(result['LibraryQualified'])
            self.assertEqual(result['PerMesh']['thin-tet']['WorkerPlusReducerSpeedup'],2)
            self.assertEqual(result['PerMesh']['fabricated-tet']['H1DofReduction'],5)
            self.assertNotIn('InterfaceMatrixRelativeDifferences',result['PerMesh']['thin-tet'])
            self.assertLess(result['DomainDefect']['tet']['RelativeToThinDomain'],1e-14)


if __name__=='__main__':unittest.main()
