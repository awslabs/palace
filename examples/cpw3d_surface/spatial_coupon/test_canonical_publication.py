# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Canonical cache, bound native-hmax, and exact rigid publication tests."""
import atexit
import copy
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest

import meshio
import numpy as np

from canonical_mesh_build import (CANONICAL_ARTIFACT_ROLES, CANONICAL_TOOL_ROLES,
                                  build_record, canonical_sha256, same_canonical_build,
                                  validate_build_record)
from mesh_stage_contract import CANONICAL_STAGE_ORDER
from publish_rigid_coupon_mesh import (_exact_mesh_structure, sha256,
                                       transform_gmsh22)
from run_native_mmg_adaptation import effective_far_size, resolve_mmg_library, run
from testdata.build_tiny_native_adapter import build as build_native_fixture, compiler
from transform_coupon_source_contract import validate_rigid_transform


HERE = Path(__file__).resolve().parent
_NATIVE_FIXTURE_DIRECTORY = tempfile.TemporaryDirectory(prefix="tiny-native-adapter-")
atexit.register(_NATIVE_FIXTURE_DIRECTORY.cleanup)
if compiler() is not None:
    ADAPTER, MMG_LIBRARY = build_native_fixture(_NATIVE_FIXTURE_DIRECTORY.name)
else:
    ADAPTER, MMG_LIBRARY = None, None


class CanonicalBuildTest(unittest.TestCase):
    GATES = {"MaximumElements": 4_000_000}

    def inputs(self):
        return {name: str(index) * 64 for index, name in enumerate(
            ("Signature", "Boundary", "Mask", "Process", "SemanticContract", "MeshRecipe"), 1)}

    def tools(self):
        return {role: canonical_sha256(role) for role in sorted(CANONICAL_TOOL_ROLES)}

    def artifacts(self, root):
        result = {}
        for role in sorted(CANONICAL_ARTIFACT_ROLES):
            path = root / f"{role}.bin"; path.write_text(role)
            result[role] = {"Path": str(path), "SHA256": sha256(path)}
        return result

    def reports(self):
        return {stage: canonical_sha256(stage) for stage in CANONICAL_STAGE_ORDER}

    def test_cache_reuse_requires_exact_source_gate_tool_and_artifact_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = build_record(self.inputs(), self.GATES, self.tools(),
                                  self.artifacts(root), self.reports())
            validate_build_record(record, self.inputs(), self.GATES, self.tools())
            self.assertTrue(same_canonical_build(record, copy.deepcopy(record)))
            self.assertIn("native-adaptation-mmg/mmg-library",
                          record["CanonicalCacheKey"]["CanonicalToolSHA256"])
            for field, key, replacement in (
                    ("SourceSHA256", "Signature", "b" * 64),
                    ("CanonicalToolSHA256", "native-adaptation-mmg/mmg-library", "c" * 64),
                    ("CanonicalToolSHA256", "native-adaptation-mmg/adapter-mmg", "e" * 64)):
                changed = copy.deepcopy(record)
                changed["CanonicalCacheKey"][field][key] = replacement
                self.assertFalse(same_canonical_build(record, changed))
                with self.assertRaisesRegex(ValueError, "cache key differs"):
                    validate_build_record(changed, self.inputs(), self.GATES, self.tools(),
                                          check_files=False)
            changed = copy.deepcopy(record)
            changed["CanonicalArtifacts"]["canonical-candidate-mesh"]["SHA256"] = "d" * 64
            self.assertFalse(same_canonical_build(record, changed))
            with self.assertRaisesRegex(ValueError, "immutable payload"):
                validate_build_record(changed, self.inputs(), self.GATES, self.tools(),
                                      check_files=False)
            changed = copy.deepcopy(record)
            changed["CanonicalStageReportSHA256"]["native-adaptation-mmg"] = "f" * 64
            self.assertFalse(same_canonical_build(record, changed))

    def test_changed_mmg_library_hash_is_a_different_canonical_build(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); artifacts = self.artifacts(root)
            record = build_record(self.inputs(), self.GATES, self.tools(), artifacts,
                                  self.reports())
            replaced = dict(self.tools())
            replaced["native-adaptation-mmg/mmg-library"] = "a" * 64
            other = build_record(self.inputs(), self.GATES, replaced, artifacts, self.reports())
            self.assertNotEqual(record["CanonicalBuildId"], other["CanonicalBuildId"])
            self.assertFalse(same_canonical_build(record, other))
            with self.assertRaisesRegex(ValueError, "cache key differs"):
                validate_build_record(record, self.inputs(), self.GATES, replaced)
            incomplete = {role: digest for role, digest in self.tools().items()
                          if role != "native-adaptation-mmg/mmg-library"}
            with self.assertRaisesRegex(ValueError, "exact canonical schema"):
                build_record(self.inputs(), self.GATES, incomplete, artifacts, self.reports())

    def test_every_canonical_artifact_role_is_required_and_no_extra_is_allowed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); artifacts = self.artifacts(root)
            self.assertIn("adaptation-receipt", CANONICAL_ARTIFACT_ROLES)
            record = build_record(self.inputs(), self.GATES, self.tools(), artifacts,
                                  self.reports())
            for role in sorted(CANONICAL_ARTIFACT_ROLES):
                with self.subTest(role=role):
                    missing = {name: item for name, item in artifacts.items() if name != role}
                    with self.assertRaisesRegex(ValueError, "exact six-stage output schema"):
                        build_record(self.inputs(), self.GATES, self.tools(), missing,
                                     self.reports())
                    deleted = copy.deepcopy(record)
                    del deleted["CanonicalArtifacts"][role]
                    payload = {key: value for key, value in deleted.items()
                               if key != "CanonicalBuildSHA256"}
                    deleted["CanonicalBuildSHA256"] = canonical_sha256(payload)
                    with self.assertRaisesRegex(ValueError, "exact six-stage output schema"):
                        validate_build_record(deleted, self.inputs(), self.GATES, self.tools(),
                                              check_files=False)
            extra = dict(artifacts)
            extra["unreviewed-extra"] = artifacts["adaptation-receipt"]
            with self.assertRaisesRegex(ValueError, "exact six-stage output schema"):
                build_record(self.inputs(), self.GATES, self.tools(), extra, self.reports())
            for stage in CANONICAL_STAGE_ORDER:
                with self.subTest(stage=stage):
                    reports = {name: digest for name, digest in self.reports().items()
                               if name != stage}
                    with self.assertRaisesRegex(ValueError, "exact canonical schema"):
                        build_record(self.inputs(), self.GATES, self.tools(), artifacts, reports)


@unittest.skipUnless(ADAPTER is not None, "native fixture adapter requires a C compiler (cc)")
class NativeHmaxTest(unittest.TestCase):
    def recipe(self, root, required_count=3):
        path = root / "recipe.json"
        path.write_text(json.dumps({
            "FarSize": 0.32, "Tetrahedra": 10,
            "RequiredTetrahedra": {"Count": required_count},
            "FarFieldBudgetPolicy": {"Name": "seed-fraction-far-field-v1",
                "RequestedFarSize": 0.16, "Pressure": 2.0,
                "EffectiveFarSize": 0.32}}))
        return path

    def inputs(self, root, required="1\n4\n9\n"):
        for name in ("seed.mesh", "metric.f64", "pins.txt", "fixed.txt"):
            (root / name).write_text(name)
        (root / "required.txt").write_text(required)
        return dict(hmin=.025, hgrad=1.3, mode="freeze-selected",
                    fixed_triangles=root / "fixed.txt",
                    required_tetrahedra=root / "required.txt")

    def test_effective_far_policy_is_actual_adapter_hmax_and_library_is_bound(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); controls = self.inputs(root)
            receipt = run(ADAPTER, root / "seed.mesh", root / "metric.f64",
                          root / "pins.txt", self.recipe(root), root / "adapted.meshb",
                          root / "receipt.json", mmg_library=MMG_LIBRARY, **controls)
            self.assertEqual(float(receipt["Command"][6]), 0.32)
            self.assertEqual(receipt["EffectiveFarSize"], 0.32)
            self.assertEqual(receipt["Command"][-2:],
                             ["--required-tetrahedra", str((root / "required.txt").resolve())])
            self.assertEqual(receipt["RequiredTetrahedra"], 3)
            self.assertEqual(receipt["RequiredTetrahedraSHA256"], sha256(root / "required.txt"))
            self.assertEqual(receipt["AdapterSHA256"], sha256(ADAPTER))
            self.assertEqual(receipt["MMGLibrarySHA256"], sha256(MMG_LIBRARY))
            self.assertEqual(Path(receipt["MMGLibraryPath"]), Path(MMG_LIBRARY).resolve())
            self.assertEqual(resolve_mmg_library(ADAPTER)["SHA256"], sha256(MMG_LIBRARY))

    def test_library_not_resolved_from_adapter_link_rpath_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); controls = self.inputs(root)
            copied = root / Path(MMG_LIBRARY).name; shutil.copyfile(MMG_LIBRARY, copied)
            with self.assertRaisesRegex(ValueError, "different MMG library"):
                run(ADAPTER, root / "seed.mesh", root / "metric.f64", root / "pins.txt",
                    self.recipe(root), root / "adapted.meshb", root / "receipt.json",
                    mmg_library=copied, **controls)
            script = HERE / "testdata" / "tiny_native_adapter.py"
            with self.assertRaisesRegex(ValueError, "exactly one MMG3D"):
                run(script, root / "seed.mesh", root / "metric.f64", root / "pins.txt",
                    self.recipe(root), root / "adapted.meshb", root / "receipt.json",
                    mmg_library=MMG_LIBRARY, **controls)
            self.assertFalse((root / "adapted.meshb").exists())

    def test_required_tetrahedron_list_is_mandatory_and_checked_against_the_recipe(self):
        for required, count, message in (("", 3, "non-empty"), ("1\n4\n4\n", 3, "duplicated"),
                                         ("0\n4\n", 2, "out of range"), ("1\n11\n", 2, "out of range"),
                                         ("1\n4\n", 3, "differs from the recipe"),
                                         ("1\nx\n", 2, "positive integers")):
            with self.subTest(required=required), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary); controls = self.inputs(root, required)
                with self.assertRaisesRegex(ValueError, message):
                    run(ADAPTER, root / "seed.mesh", root / "metric.f64", root / "pins.txt",
                        self.recipe(root, count), root / "adapted.meshb", root / "receipt.json",
                        mmg_library=MMG_LIBRARY, **controls)
                self.assertFalse((root / "adapted.meshb").exists())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); controls = self.inputs(root)
            (root / "required.txt").unlink()
            with self.assertRaisesRegex(ValueError, "must exist"):
                run(ADAPTER, root / "seed.mesh", root / "metric.f64", root / "pins.txt",
                    self.recipe(root), root / "adapted.meshb", root / "receipt.json",
                    mmg_library=MMG_LIBRARY, **controls)
            recipe = root / "recipe.json"
            recipe.write_text(json.dumps({"FarSize": 0.32, "Tetrahedra": 10,
                "FarFieldBudgetPolicy": {"Name": "seed-fraction-far-field-v1",
                    "RequestedFarSize": 0.16, "Pressure": 2.0, "EffectiveFarSize": 0.32}}))
            (root / "required.txt").write_text("1\n")
            with self.assertRaisesRegex(ValueError, "required-tetrahedra record"):
                run(ADAPTER, root / "seed.mesh", root / "metric.f64", root / "pins.txt",
                    recipe, root / "adapted.meshb", root / "receipt.json",
                    mmg_library=MMG_LIBRARY, **controls)

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
