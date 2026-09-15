# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Canonical cache, bound native-hmax, and exact rigid publication tests."""
import copy
import json
import math
from pathlib import Path
import shutil
import stat
import tempfile
import unittest

import meshio
import numpy as np

from canonical_mesh_build import (build_record, same_canonical_build,
                                  validate_build_record)
from publish_rigid_coupon_mesh import (_exact_mesh_structure, sha256,
                                       transform_gmsh22)
from run_native_mmg_adaptation import effective_far_size, run
from transform_coupon_source_contract import validate_rigid_transform


class CanonicalBuildTest(unittest.TestCase):
    def inputs(self):
        return {name: str(index) * 64 for index, name in enumerate(
            ("Signature", "Boundary", "Mask", "Process", "SemanticContract", "MeshRecipe"), 1)}

    def test_cache_reuse_requires_exact_source_gate_tool_and_artifact_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); mesh = root / "canonical.msh"; mesh.write_text("mesh")
            artifacts = {"candidate-mesh": {"Path": str(mesh), "SHA256": sha256(mesh)}}
            record = build_record(self.inputs(), {"MaximumElements": 4_000_000},
                                  {"adapter-mmg": "a" * 64}, artifacts)
            validate_build_record(record, self.inputs(), {"MaximumElements": 4_000_000},
                                  {"adapter-mmg": "a" * 64})
            self.assertTrue(same_canonical_build(record, copy.deepcopy(record)))
            for field, replacement in (("SourceSHA256", "b" * 64),
                                       ("CanonicalToolSHA256", "c" * 64)):
                changed = copy.deepcopy(record)
                key = next(iter(changed["CanonicalCacheKey"][field]))
                changed["CanonicalCacheKey"][field][key] = replacement
                self.assertFalse(same_canonical_build(record, changed))
                with self.assertRaises(ValueError):
                    validate_build_record(changed, self.inputs(),
                                          {"MaximumElements": 4_000_000},
                                          {"adapter-mmg": "a" * 64}, check_files=False)
            changed = copy.deepcopy(record)
            changed["CanonicalArtifacts"]["candidate-mesh"]["SHA256"] = "d" * 64
            self.assertFalse(same_canonical_build(record, changed))


class NativeHmaxTest(unittest.TestCase):
    def recipe(self, root):
        path = root / "recipe.json"
        path.write_text(json.dumps({
            "FarSize": 0.32,
            "FarFieldBudgetPolicy": {"Name": "seed-fraction-far-field-v1",
                "RequestedFarSize": 0.16, "Pressure": 2.0,
                "EffectiveFarSize": 0.32}}))
        return path

    def test_effective_far_policy_is_actual_adapter_hmax(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("seed.mesh", "metric.f64", "pins.txt", "fixed.txt"):
                (root / name).write_text(name)
            adapter = root / "adapter"
            adapter.write_text("#!/usr/bin/env python3\nimport pathlib,shutil,sys\n"
                               "shutil.copyfile(sys.argv[1],sys.argv[4])\n")
            adapter.chmod(adapter.stat().st_mode | stat.S_IXUSR)
            receipt = run(adapter, root / "seed.mesh", root / "metric.f64",
                          root / "pins.txt", self.recipe(root), root / "adapted.meshb",
                          root / "receipt.json", hmin=.025, hgrad=1.3,
                          mode="freeze-selected", fixed_triangles=root / "fixed.txt")
            self.assertEqual(float(receipt["Command"][6]), 0.32)
            self.assertEqual(receipt["EffectiveFarSize"], 0.32)

    def test_recorded_policy_not_consumed_and_tampered_hmax_are_rejected(self):
        recipe = {"FarSize": .16, "FarFieldBudgetPolicy": {
            "Name": "seed-fraction-far-field-v1", "RequestedFarSize": .16,
            "Pressure": 2., "EffectiveFarSize": .32}}
        with self.assertRaisesRegex(ValueError, "not consumed"):
            effective_far_size(recipe)
        recipe["FarSize"] = .31
        recipe["FarFieldBudgetPolicy"]["EffectiveFarSize"] = .31
        with self.assertRaisesRegex(ValueError, "internally inconsistent"):
            effective_far_size(recipe)


class ExactRigidGmshTest(unittest.TestCase):
    def mesh(self, path):
        points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        cells = [("triangle", np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])),
                 ("tetra", np.array([[0, 1, 2, 3]]))]
        data = {"gmsh:physical": [np.full(4, 3), np.array([7])],
                "gmsh:geometrical": [np.full(4, 13), np.array([17])]}
        field = {"interface": np.array([3, 2]), "vacuum": np.array([7, 3])}
        meshio.write(path, meshio.Mesh(points, cells, cell_data=data, field_data=field),
                     file_format="gmsh22", binary=True)

    def test_identity_is_fresh_and_rotation_preserves_exact_structure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root / "source.msh"; self.mesh(source)
            identity = np.eye(4)
            identity_output = root / "identity.msh"
            before, after, error = transform_gmsh22(source, identity_output, identity)
            self.assertEqual(source.read_bytes(), identity_output.read_bytes())
            self.assertEqual(before["tags"], after["tags"]); self.assertEqual(error, 0.)
            angle = .63
            matrix = validate_rigid_transform([
                math.cos(angle), -math.sin(angle), 0, 1.2,
                math.sin(angle), math.cos(angle), 0, -.7,
                0, 0, 1, .9, 0, 0, 0, 1])
            rotated = root / "rotated.msh"
            before, after, error = transform_gmsh22(source, rotated, matrix)
            self.assertNotEqual(sha256(source), sha256(rotated))
            self.assertEqual(before["tags"], after["tags"])
            self.assertLessEqual(error, 1e-15)
            left, right = meshio.read(source), meshio.read(rotated)
            self.assertTrue(_exact_mesh_structure(left, right))
            expected = left.points @ np.asarray(matrix)[:3, :3].T + np.asarray(matrix)[:3, 3]
            np.testing.assert_allclose(right.points, expected, rtol=0, atol=1e-15)

    def test_reflection_and_preexisting_output_are_rejected(self):
        with self.assertRaises(ValueError):
            validate_rigid_transform([-1., 0, 0, 0, 0, 1., 0, 0,
                                      0, 0, 1., 0, 0, 0, 0, 1.])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root / "source.msh"; self.mesh(source)
            output = root / "output.msh"; output.write_text("occupied")
            with self.assertRaises(ValueError):
                transform_gmsh22(source, output, np.eye(4))


if __name__ == "__main__":
    unittest.main()
