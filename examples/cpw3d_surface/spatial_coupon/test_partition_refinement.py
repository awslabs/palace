# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import importlib.util
from pathlib import Path
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


if __name__=='__main__':
    unittest.main()
