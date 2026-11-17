# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""run_graded_library_mesh.check_mesh_preprocessing_policy (R4 T3 #4; SD-001; decision 1570): an unset
crack key is upstream's default (null in mesh-statistics.json), not true."""
import unittest

from run_graded_library_mesh import check_mesh_preprocessing_policy


class MeshPreprocessingPolicyTest(unittest.TestCase):
    def test_unset_crack_key_matches_null_only(self):
        check_mesh_preprocessing_policy({'CrackInternalBoundaryElements':None,'RefineCrackElements':True},{})
        with self.assertRaises(ValueError):
            check_mesh_preprocessing_policy({'CrackInternalBoundaryElements':True,'RefineCrackElements':True},{})

    def test_explicit_values_must_match(self):
        for value in (True,False):
            check_mesh_preprocessing_policy({'CrackInternalBoundaryElements':value,'RefineCrackElements':False},
                                            {'CrackInternalBoundaryElements':value,'RefineCrackElements':False})
        with self.assertRaises(ValueError):
            check_mesh_preprocessing_policy({'CrackInternalBoundaryElements':None,'RefineCrackElements':True},
                                            {'CrackInternalBoundaryElements':True})
        with self.assertRaises(ValueError):
            check_mesh_preprocessing_policy({'CrackInternalBoundaryElements':True,'RefineCrackElements':True},
                                            {'CrackInternalBoundaryElements':True,'RefineCrackElements':False})

    def test_a_missing_record_field_fails_closed(self):
        with self.assertRaises(KeyError):
            check_mesh_preprocessing_policy({'RefineCrackElements':True},{})


if __name__=='__main__':unittest.main()
