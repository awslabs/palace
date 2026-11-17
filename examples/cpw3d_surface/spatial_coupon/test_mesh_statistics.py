# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


@unittest.skipUnless(os.environ.get('PALACE_MESH_STATS_BINARY'),'Set the native Palace mesh-statistics executable')
class MeshStatisticsTest(unittest.TestCase):
    def test_internal_pec_crack_policy_matches_solver_h1_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);mesh=root/'two-tetrahedra.mesh'
            mesh.write_text('''MFEM mesh v1.0

dimension
3

elements
2
1 4 0 1 2 3
2 4 0 2 1 4

boundary
7
4 2 0 1 2
1 2 0 3 1
1 2 1 3 2
1 2 2 3 0
1 2 0 1 4
1 2 1 2 4
1 2 2 0 4

vertices
5
3
0 0 0
1 0 0
0 1 0
0 0 1
0 0 -1
''')
            cfg={'Problem':{'Type':'Electrostatic','Output':str(root/'statistics')},
                 'Model':{'Mesh':str(mesh),'L0':1e-6,'Refinement':{'MaxIts':0}},
                 'Domains':{'Materials':[{'Attributes':[1,2],'Permittivity':1.0}]},
                 'Boundaries':{'Ground':{'Attributes':[4]},'Terminal':[{'Index':1,'Attributes':[1]}]},
                 'Solver':{'Order':2,'Electrostatic':{'Save':0},
                           'Linear':{'Type':'BoomerAMG','KSPType':'CG','Tol':1e-8,'MaxIts':100}}}
            executable=os.environ['PALACE_MESH_STATS_BINARY']
            # The raw conforming mesh has 5 vertices + 9 edges = 14 p2 DOFs. The legacy key cracks the
            # full internal Ground sheet, which separates two p2 tetrahedra (20); with the key unset
            # (upstream's default since PR 994, recorded as null: SD-001) a Ground sheet is not cracked.
            for tag,crack,record,dofs in (('legacy',True,True,20),('default',None,None,14)):
                cfg['Model'].pop('CrackInternalBoundaryElements',None)
                if crack is not None:cfg['Model']['CrackInternalBoundaryElements']=crack
                cfg['Problem']['Output']=str(root/f'statistics-{tag}')
                config=root/f'input-{tag}.json';config.write_text(json.dumps(cfg))
                subprocess.run([executable,'--mesh-statistics',str(config)],check=True,capture_output=True,text=True)
                result=json.loads((root/f'statistics-{tag}/mesh-statistics.json').read_text())
                self.assertIs(result['CrackInternalBoundaryElements'],record)
                self.assertTrue(result['RefineCrackElements'])
                self.assertEqual(result['H1TrueDOFs'],dofs)
                cfg['Problem']['Output']=str(root/f'solve-{tag}');config.write_text(json.dumps(cfg))
                run=subprocess.run([executable,str(config)],check=True,capture_output=True,text=True)
                actual=int(re.search(r'H1 \(p = 2\): (\d+)',run.stdout)[1])
                self.assertEqual(actual,result['H1TrueDOFs'])


if __name__=='__main__':unittest.main()
