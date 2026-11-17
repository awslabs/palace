# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('partition_study',Path(__file__).with_name('run_partition_refinement.py'))
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PartitionStudyTest(unittest.TestCase):
    def test_all_conductor_parts_remain_in_prescribed_boundary_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            attributes=[3100,13100,5001,15001,6001,16001,5002,6002]
            config,mapping=module.make_probe(root/'mesh.msh',attributes,root/'postpro',True)
            boundary=config['Boundaries']
            self.assertEqual(set(boundary['Ground']['Attributes']),set(attributes)-{3100,13100})
            states=boundary['PrescribedPotential']
            self.assertEqual(set(states[0]['TerminalAttributes']),{5001,15001,6001,16001})
            self.assertEqual(set(states[1]['TerminalAttributes']),{5002,6002})
            self.assertTrue(all(p['Attributes']==[1] for p in states))
            self.assertEqual(len(mapping),len(attributes))
            self.assertTrue(config['Solver']['Electrostatic']['ResponseMatrix'])

    def test_disjoint_energy_parts_and_conservative_group_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'domain-E.csv').write_text('i,E_elec (J)\n1,2\n')
            (root/'surface-Q.csv').write_text('i,p_surf[1],p_surf[2],p_surf[3]\n1,0.1,0.2,0.3\n')
            mapping={i:dict(CanonicalAttribute=a,Group='SA:etched',Unresolved=u)
                     for i,a,u in [(1,3100,False),(2,3100,True),(3,3101,False)]}
            result=module.energies(root,mapping)['1']
            self.assertAlmostEqual(result['Groups']['SA:etched'],1.2)
            self.assertAlmostEqual(result['UnresolvedEnergy']['SA:etched'],.4)
            self.assertAlmostEqual(result['Slots']['3100'],.6)
            self.assertAlmostEqual(result['GroupNormalizedPartitionBound']['SA:etched'],1/3)
            self.assertAlmostEqual(result['SlotRelativePartitionBounds']['3100'],2/3)
            (root/'domain-E.csv').write_text('i,E_elec (J)\n1,0\n')
            with self.assertRaises(ValueError):module.energies(root,mapping)

    def test_tool_snapshot_is_self_contained(self):
        # SL9 review m-4 (upstream integration P2): the archived run_partition_refinement.py imports
        # crack_policy; the snapshot carries it, and the archived copy runs from the snapshot alone.
        with tempfile.TemporaryDirectory() as tmp:
            tools=Path(tmp)/'archive'/'tools';tools.mkdir(parents=True)
            names=module.snapshot_tools(tools)
            self.assertIn('crack_policy.py',names)
            self.assertEqual(sorted(p.name for p in tools.iterdir()),sorted(names))
            env={k:v for k,v in os.environ.items() if k!='PYTHONPATH'}
            env['PYTHONDONTWRITEBYTECODE']='1'
            result=subprocess.run([sys.executable,str(tools/'run_partition_refinement.py'),'--help'],cwd=tmp,env=env,
                                  capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('--crack-policy',result.stdout)


if __name__=='__main__':
    unittest.main()
