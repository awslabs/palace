# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import tempfile
import json
import os
import subprocess
import sys
from unittest.mock import patch
import unittest
from pathlib import Path
from run_graded_library_case import validate_matrix, compact_corner, compare_unchanged_control, matrix_contract
from run_graded_library_mesh import compare_geometry, mesh_recipe, trace_policy_environment, mesh_environment, main as mesh_main
from archive_storage import reserve, seal, raw_tetrahedral_h1_dofs
from submit_graded_library_campaign import main as submit_main


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
                         'frozen_volume_study.jl','graded_curve_distance.jl','graded_size_points.jl','graded_trace_size.jl',
                         'surface_ribbon_constraints.jl','interface_ownership.jl',
                         'ownership_bernstein.jl','label_interface_patches.jl','relabel_frozen_interface_mesh.jl'):
                (tools/name).write_text('tool')
            case={'InputDirectory':str(inputs),'Kind':'thin','SurfaceSize':.002,'TraceSize':.01,'TraceSizeScope':'matching','SurfaceAlgorithm':5,'TraceConstraintMode':'all'}
            before=mesh_recipe(root,case)
            (inputs/'traces/basis-0001.csv').write_text('changed connectivity')
            self.assertNotEqual(before,mesh_recipe(root,case))
            (inputs/'traces/basis-0001.csv').write_text('original')
            for field,value in (('TraceSize',.02),('TraceSizeScope','matching-and-volume'),
                                ('TraceConstraintMode','levels'),('TraceRelativeSize',.2),('TraceSurfaceGrowth',2.)):
                with self.subTest(field=field):
                    self.assertNotEqual(before,mesh_recipe(root,dict(case,**{field:value})))
            for name in ('graded_trace_size.jl','mesh_spatial_coupon.jl'):
                (tools/name).write_text('changed policy implementation')
                self.assertNotEqual(before,mesh_recipe(root,case))
                (tools/name).write_text('tool')
            del case['TraceSizeScope']
            with self.assertRaisesRegex(ValueError,'explicit TraceSizeScope'):mesh_recipe(root,case)

    def test_trace_policy_defaults_and_rejections(self):
        self.assertEqual(trace_policy_environment({})['TET_TRACE_SIZE_SCOPE'],'off')
        for scope in ('matching','matching-and-volume','legacy-global'):
            valid={'TraceSize':.01,'TraceSizeScope':scope,'TraceConstraintMode':'levels'}
            env=trace_policy_environment(valid)
            self.assertEqual(env['TET_TRACE_SIZE_SCOPE'],scope)
            self.assertEqual(float(env['TET_TRACE_SURFACE_GROWTH']),.5 if scope=='legacy-global' else 4.)
            with self.assertRaises(ValueError):trace_policy_environment(dict(valid,TraceSize=0.))
            with self.assertRaises(ValueError):trace_policy_environment(dict(valid,TraceConstraintMode='none'))
            del valid['TraceConstraintMode']
            with self.assertRaisesRegex(ValueError,'explicit TraceConstraintMode'):trace_policy_environment(valid)
        for bad in ({'TraceSize':.01},{'TraceSize':float('nan')},{'TraceSize':-1},
                    {'TraceRelativeSize':.2},{'TraceSizeScope':'global'},{'TraceConstraintMode':'bad'},
                    {'TraceSize':True},{'TraceRelativeSize':float('inf')},{'TraceSurfaceGrowth':0},
                    {'TraceSize':.01,'TraceSizeScope':'legacy-global','TraceConstraintMode':'all','TraceSurfaceGrowth':4}):
            with self.subTest(bad=bad),self.assertRaises(ValueError):trace_policy_environment(bad)
        self.assertEqual(trace_policy_environment({'TraceConstraintMode':'none'})['TET_TRACE_SIZE_SCOPE'],'off')

    def test_campaign_environment_does_not_inherit_trace_controls(self):
        case={'SurfaceAlgorithm':5,'TraceSize':0,'TraceSizeScope':'off','TraceConstraintMode':'none'}
        with patch.dict(os.environ,{'TET_TRACE_SIZE':'10','TET_TRACE_SIZE_SCOPE':'legacy-global',
                                   'TET_TRACE_SURFACE_GROWTH':'99','TET_TRACE_RELATIVE_SIZE':'.3',
                                   'TET_EDGE_TANGENT_SIZE':'20'}):
            env=mesh_environment(case)
        self.assertEqual(env['TET_TRACE_SIZE'],'0')
        self.assertEqual(env['TET_TRACE_SIZE_SCOPE'],'off')
        self.assertEqual(env['TET_TRACE_CONSTRAINT_MODE'],'none')
        self.assertEqual(env['TET_TRACE_SURFACE_GROWTH'],'4.0')
        self.assertEqual(env['TET_TRACE_RELATIVE_SIZE'],'0.0')
        self.assertNotIn('TET_EDGE_TANGENT_SIZE',env)

    def test_unsafe_historical_manifest_fails_before_cache_or_job(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);directory=root/'mesh';directory.mkdir()
            state=directory/'mesh-state.json'
            state.write_text(json.dumps({'State':'Completed','Recipe':{},'MeshSHA256':'old'}))
            original=state.read_bytes()
            (root/'campaign.json').write_text(json.dumps({'Cases':[{'Key':'old','MeshRequired':True,
                'Mesh':str(directory/'coupon.msh'),'TraceSize':.01,'TraceConstraintMode':'all'}]}))
            with patch.object(sys,'argv',['run',str(root),'old']),patch('subprocess.run') as run:
                with self.assertRaisesRegex(ValueError,'explicit TraceSizeScope'):mesh_main()
                run.assert_not_called()
            self.assertEqual(state.read_bytes(),original)

    def test_unsafe_manifest_rejected_before_scheduler_access(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            (root/'campaign.json').write_text(json.dumps({'Cases':[{'Key':'old',
                'MeshRequired':True,'TraceSize':.01,'TraceConstraintMode':'levels'}]}))
            for phase in ('meshes','responses','all'):
                with self.subTest(phase=phase), patch.object(sys,'argv',
                        ['submit',str(root),'--phase',phase]), \
                        patch('subprocess.check_output') as scheduler:
                    with self.assertRaisesRegex(ValueError,'explicit TraceSizeScope'):
                        submit_main()
                    scheduler.assert_not_called()
                    self.assertFalse((root/'submit.lock').exists())
                    self.assertFalse((root/'submissions.json').exists())

    def test_new_campaign_requires_explicit_policy_before_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'new'
            command=[sys.executable,str(Path(__file__).with_name('prepare_graded_library_campaign.py')),
                     str(root),'--reference-manifest','missing','--geometry-audit-root','missing',
                     '--python',sys.executable,'--julia-project','missing','--trace-mode','levels',
                     '--mesh-statistics-binary','missing']
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('--trace-size-scope',result.stderr)
            self.assertFalse(root.exists())
            result=subprocess.run(command+['--trace-size-scope','matching'],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('positive TraceSize',result.stderr)
            self.assertFalse(root.exists())

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
