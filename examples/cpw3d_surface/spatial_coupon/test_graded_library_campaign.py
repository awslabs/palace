# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import tempfile
import unittest
from pathlib import Path
from run_graded_library_case import validate_matrix, compact_corner, compare_unchanged_control, matrix_contract
from run_graded_library_mesh import compare_geometry, mesh_recipe
from archive_storage import reserve, seal, raw_tetrahedral_h1_dofs


class GradedCampaignTest(unittest.TestCase):
    def test_complete_symmetric_psd_matrix_and_missing_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            p=Path(temporary)/'matrix.csv'
            p.write_text('basis_i,basis_j,Q_ij (J)\n1,1,2\n1,3,-1\n3,3,2\n')
            self.assertEqual(validate_matrix(p,[1,3])['BasisSize'],2)
            self.assertEqual(compare_unchanged_control(p,p),{'Q_ij (J)':0.})
            with self.assertRaises(ValueError):validate_matrix(p,[1,2,3])
            p.write_text('basis_i,basis_j,Q_ij (J)\n1,1,1\n1,3,2\n3,3,1\n')
            with self.assertRaisesRegex(ValueError,'Nonpositive'):validate_matrix(p,[1,3])

    def test_corner_core_convention_and_lf(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);a,b,c=[root/name for name in ('raw.csv','reference.csv','compact.csv')]
            a.write_text('interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)\n1,1,2e-6,1,1,3,8\n')
            b.write_text('interface,edge,basis_i,basis_j,Q_total_ij (J)\n1,1,1,1,3\n')
            self.assertEqual(compact_corner(a,b,c),c)
            self.assertNotIn(b'\r',c.read_bytes())
            self.assertEqual(compare_unchanged_control(c,b),{'Q_total_ij (J)':0.})

    def test_missing_energy_column_or_whole_interface_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);reference=root/'ref.csv';actual=root/'actual.csv'
            reference.write_text('interface,edge,basis_i,basis_j,Q_total_ij (J)\n1,1,1,1,2\n2,1,1,1,3\n')
            contract=matrix_contract(reference,'surface')
            actual.write_text('interface,edge,basis_i,basis_j,Q_total_ij (J)\n1,1,1,1,2\n')
            with self.assertRaises(ValueError):validate_matrix(actual,[1],'surface',contract)
            actual.write_text('basis_i,basis_j\n1,1\n')
            with self.assertRaises(ValueError):validate_matrix(actual,[1])
            with self.assertRaises(ValueError):compare_unchanged_control(actual,reference)
            actual.write_text('basis_i,basis_j,Q_ij (J)\n1,1,2\n')
            with self.assertRaises(ValueError):validate_matrix(actual,[1],'surface',contract)

    def test_storage_admissions_account_for_outstanding_workers(self):
        with tempfile.TemporaryDirectory() as root:
            reserve(root,'a',200,margin=0,free_bytes=350)
            with self.assertRaises(RuntimeError):reserve(root,'b',200,margin=0,free_bytes=350)
            seal(root,'a',200)
            with self.assertRaises(RuntimeError):reserve(root,'b',200,margin=0,free_bytes=150)
            reserve(root,'b',100,margin=0,free_bytes=150)
            with self.assertRaises(ValueError):reserve(root,'b',100,margin=0,free_bytes=150)

    def test_recipe_changes_when_trace_bytes_or_size_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);inputs=root/'inputs';tools=root/'tools';inputs.mkdir();tools.mkdir();(inputs/'traces').mkdir()
            for name in ('mesh-signature.csv','plan-view-mask.csv','plan-view-boundary.csv','process.toml','traces/basis-0001.csv'):
                (inputs/name).write_text('original')
            for name in ('run_graded_library_mesh.py','mesh_spatial_coupon.jl','mesh_graded_tet_experiment.jl',
                         'frozen_volume_study.jl','graded_curve_distance.jl','graded_size_points.jl','interface_ownership.jl',
                         'ownership_bernstein.jl','label_interface_patches.jl','relabel_frozen_interface_mesh.jl'):
                (tools/name).write_text('tool')
            case={'InputDirectory':str(inputs),'Kind':'thin','SurfaceSize':.002,'TraceSize':.01,'SurfaceAlgorithm':5,'TraceConstraintMode':'all'}
            before=mesh_recipe(root,case)
            (inputs/'traces/basis-0001.csv').write_text('changed connectivity')
            self.assertNotEqual(before,mesh_recipe(root,case))
            (inputs/'traces/basis-0001.csv').write_text('original');case['TraceSize']=.02
            self.assertNotEqual(before,mesh_recipe(root,case))

    def test_raw_tetrahedral_h1_topology_diagnostic(self):
        counts={'Vertices':4,'Edges':6,'Faces':4,'Elements':1}
        self.assertEqual(raw_tetrahedral_h1_dofs(counts,4),35)
        self.assertEqual(raw_tetrahedral_h1_dofs(counts,5),56)

    def test_geometry_groups_slots_but_not_process_planes_or_conductors(self):
        ref=[(3,1,10.),(2,3000,1.),(2,3001,2.),(2,3100,4.),(2,5001,6.),(2,5101,3.)]
        cad=[(3,1,10.),(2,3000,3.),(2,3101,4.),(2,5001,9.)]
        self.assertEqual(compare_geometry(ref,cad),0.)
        with self.assertRaises(ValueError):compare_geometry(ref,[(d,a+(1 if a==5001 else 0),v) for d,a,v in cad])
        with self.assertRaises(ValueError):compare_geometry(ref,[(d,3000 if a==3101 else a,v) for d,a,v in cad])


if __name__=='__main__':unittest.main()
