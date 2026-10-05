#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The frozen tree's vendored MFEM patches are never silently dropped by build_from_cmake_tree.py.

    python3 -m unittest test_build_from_cmake_tree
"""
from pathlib import Path
import tempfile
import unittest

from build_from_cmake_tree import MFEM_PATCH_DIRECTORY, frozen_mfem_patches, resolve_mfem_patches

PR5494 = f"{MFEM_PATCH_DIRECTORY}/mfem_pr5494.diff"
OTHER = f"{MFEM_PATCH_DIRECTORY}/mfem_pr9999.diff"


class ResolveMFEMPatches(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.source = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def freeze(self, *patches):
        """A frozen tree carrying the given diffs; returns its manifest file set."""
        (self.source / MFEM_PATCH_DIRECTORY).mkdir(parents=True, exist_ok=True)
        for relative in patches:
            (self.source / relative).write_text("diff --git a/x b/x\n")
        return {p: "sha" for p in patches} | {"palace/CMakeLists.txt": "sha"}

    def test_tree_without_patches_builds_against_the_prefix_mfem(self):
        files = self.freeze()
        self.assertEqual(frozen_mfem_patches(self.source), [])
        self.assertEqual(resolve_mfem_patches(self.source, files, []), [])

    def test_old_qsub_line_without_mfem_patch_applies_the_frozen_patches(self):
        files = self.freeze(PR5494)
        self.assertEqual(frozen_mfem_patches(self.source), [PR5494])
        self.assertEqual(resolve_mfem_patches(self.source, files, []), [PR5494])

    def test_explicit_full_list_is_kept_in_the_given_order(self):
        files = self.freeze(PR5494, OTHER)
        self.assertEqual(resolve_mfem_patches(self.source, files, [OTHER, PR5494]), [OTHER, PR5494])

    def test_explicit_list_omitting_a_frozen_patch_is_refused(self):
        files = self.freeze(PR5494, OTHER)
        with self.assertRaises(SystemExit) as refused:
            resolve_mfem_patches(self.source, files, [PR5494])
        self.assertIn(OTHER, str(refused.exception))

    def test_patch_outside_the_manifest_is_refused(self):
        files = self.freeze(PR5494)
        stray = f"{MFEM_PATCH_DIRECTORY}/stray.diff"
        (self.source / stray).write_text("diff --git a/y b/y\n")
        with self.assertRaises(SystemExit) as refused:
            resolve_mfem_patches(self.source, files, [])
        self.assertIn(stray, str(refused.exception))

    def test_only_diff_files_of_the_mfem_patch_directory_count(self):
        files = self.freeze(PR5494)
        (self.source / MFEM_PATCH_DIRECTORY / "README.md").write_text("not a patch\n")
        (self.source / "extern/patch/mumps").mkdir(parents=True)
        (self.source / "extern/patch/mumps/patch_build.diff").write_text("diff --git a/z b/z\n")
        self.assertEqual(resolve_mfem_patches(self.source, files, []), [PR5494])


if __name__ == "__main__":
    unittest.main()
