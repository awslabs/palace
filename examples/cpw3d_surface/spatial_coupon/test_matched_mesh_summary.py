# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import json
from pathlib import Path
import tempfile
import unittest
from summarize_matched_mesh_benchmark import summarize,worst_energy_difference


class MatchedSummaryTest(unittest.TestCase):
    def test_worst_energy_keeps_weak_modes(self):
        reference={(1,1):2.,(1,2):0.,(2,1):0.,(2,2):1.}
        candidate=dict(reference);candidate[(1,1)]=2.2
        self.assertAlmostEqual(worst_energy_difference(reference,candidate),.1)
        reference[(2,2)]=0.
        self.assertIsNone(worst_energy_difference(reference,candidate))

    def test_recovered_reference_is_only_a_speedup_lower_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            cases=[dict(Name=f'fabricated-{s}',Kind='fabricated',Strategy=s,MeshSHA256=s,
                        Elements=20,Kappa={'100.000000':10}) for s in ('prism','tet')]
            (root/'manifest.json').write_text(json.dumps(dict(Cases=cases,InputSHA256={})))
            reducer=dict(H1Dofs=100,Seconds=10,PeakProcessTreeRSSBytes=1000)
            ledger=dict(Host='test',ExecutableSHA256='binary',Ranks=2,Order=5,SolverTolerance=1e-8,
                        Cases={'fabricated-tet':dict(Completed=True,Worker=dict(H1Dofs=10,Seconds=6,Sources=[]),
                                Reducer=dict(reducer,Seconds=4),TotalSeconds=10,PeakProcessTreeRSSBytes=1000)})
            (root/'timings.json').write_text(json.dumps(ledger))
            recovery=root/'recovery/fabricated-prism';recovery.mkdir(parents=True)
            (recovery/'summary.json').write_text(json.dumps(dict(Case='fabricated-prism',Completed=True,
                Reducer=reducer,Sources=[],Attempts=[],TotalSecondsLowerBound=110,
                ActualSecondsIncludingInterruptedWork=200)))
            for case in cases:
                out=root/'postpro'/f"{case['Name']}-reducer";out.mkdir(parents=True)
                (out/'domain-response-matrix.csv').write_text('basis_i,basis_j,Q_ij (J)\n1,1,1\n')
                (out/'surface-response-matrix.csv').write_text('interface,basis_i,basis_j,Q_total_ij (J)\n1,1,1,.1\n')
            result=summarize([root]);tet=result['Cases']['fabricated-tet']
            self.assertTrue(tet['SpeedupIsLowerBound'])
            self.assertEqual(tet['Speedup'],11)
            self.assertTrue(result['Cases']['fabricated-prism']['TimingIsLowerBound'])
            self.assertFalse(result['LibraryQualified'])


if __name__=='__main__':unittest.main()
