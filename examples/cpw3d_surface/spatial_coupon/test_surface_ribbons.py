# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from run_general_mesh_suite import fixtures, geometry_gate, measures

ROOT=Path(__file__).resolve().parent


class SurfaceRibbonTest(unittest.TestCase):
    def test_rows_preserve_two_conductor_geometry(self):
        julia=shutil.which('julia')
        if not julia:self.skipTest('Julia unavailable')
        project=Path(os.environ.get('PALACE_JULIA_PROJECT',ROOT.parents[2]/'test/examples'))
        probe=subprocess.run([julia,f'--project={project}','-e','using Gmsh'],capture_output=True)
        if probe.returncode:self.skipTest('No instantiated Gmsh project')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);case=next(p for p in fixtures(root) if p.name=='two-conductors')
            env=dict(os.environ,JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                     TET_GEOMETRY_ORDER='1',TET_SURFACE_RIBBON_ROWS='2',TET_SURFACE_RIBBON_ASPECT='2')
            base=[julia,'--startup-file=no',f'--project={project}',str(ROOT/'mesh_graded_tet_experiment.jl'),str(case),'fabricated']
            options=['1','.01','.4','--process',str(case/'process.toml')]
            subprocess.run([*base,str(root/'cad.msh'),*options,'--geometry-only'],env=env,
                           check=True,capture_output=True,text=True,timeout=90)
            import json
            expected=json.loads((case/'expected.json').read_text())
            geometry_gate(expected,measures(root/'cad.msh.cad-measures.csv'),'fabricated')
            study=root/'volume.toml'
            study.write_text('MaxElements=200000\nMaxNodes=80000\nOptimizeVolume=false\nRepairInvalidVolume=true\n'
                             '[[Variants]]\nName="test"\nMinimumSize=0.01\nNearGrowth=0.5\nFarGrowth=2.0\n'
                             'TransitionDistance=0.03\nMaximumSize=0.4\n')
            subprocess.run([*base,str(root/'ribbon.msh'),*options,'--reference-measures',
                            str(root/'cad.msh.cad-measures.csv'),'--volume-study',str(study)],env=env,
                           check=True,capture_output=True,text=True,timeout=90)
            with (root/'ribbon-test.msh').open('rb') as stream:
                self.assertEqual(stream.readline().strip(),b'$MeshFormat')
                self.assertEqual(stream.readline().strip(),b'2.2 1 8')


if __name__=='__main__':unittest.main()
