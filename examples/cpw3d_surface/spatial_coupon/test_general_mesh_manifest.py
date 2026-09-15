# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import copy
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import tempfile
import unittest

import meshio
import numpy as np

from audit_edge_metric_mesh import analyze
from general_mesh_audit_producer import (KINDS, _footprint_boundary_comparison,
                                         _footprint_boundary_distance,
                                         _global_diagonal_bands,
                                         _normalized_footprint_boundary,
                                         _protected_surface_report,
                                         produce as produce_audit)
from general_mesh_manifest import run_manifest, sha256, validate_manifest
from mesh_array_io import read_mesh
from mesh_stage_contract import validate_tool_invocation
from normalize_general_mesh_evidence import normalize
from semantic_mesh_contract import (derive_feature_topology, validate_semantic_contract)


HERE = Path(__file__).resolve().parent
MESHER = HERE / "testdata" / "tiny_mesh_audit_producer.py"
AUDITOR = HERE / "general_mesh_audit_producer.py"
BOUNDED = HERE / "run_bounded_mesher.py"
STAGER = HERE / "testdata" / "tiny_mesh_stage.py"
IDENTITY = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
ANGLE = 0.63
ROTATION = [math.cos(ANGLE), -math.sin(ANGLE), 0, 0,
            math.sin(ANGLE), math.cos(ANGLE), 0, 0,
            0, 0, 1, 0, 0, 0, 0, 1]


class GeneralMeshManifestTest(unittest.TestCase):
    def write_case(self, root, name, edges, *, scale):
        directory = root / name; directory.mkdir()
        signature = directory / "signature.csv"
        with signature.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Index", "Slot", "Conductor", "Px", "Py", "Pz",
                             "Tx", "Ty", "Tz", "S0", "S1"])
            for i in range(edges):
                if name == "subdivided":
                    slot, conductor, y = 0, 1, 0
                else:
                    slot, conductor, y = i % 2, i // 2 + 1, i
                writer.writerow([i + 1, slot, conductor, (i + .5) * scale,
                                 y * scale, 0, 1, 0, 0, -.5 * scale, .5 * scale])
        (directory / "boundary.csv").write_text(
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
            f"1,1,1,0,0,Physical,{10 * scale},0\n"
            f"1,2,1,0,0,Continuation,{2 * scale},0\n")
        (directory / "mask.csv").write_text("mask\n")
        (directory / "process.toml").write_text(
            f'Units = "um"\nMetalThickness = {0.004 * scale}\n')
        recipe = {"Version": 1, "GeometryOrder": 1,
                  "NormalSizeOverThickness": .25,
                  "TangentialSizeOverNormalSize": 4.0}
        (directory / "recipe.json").write_text(json.dumps(recipe) + "\n")
        labels = [{"Attribute": 1, "Role": "air-outer", "AdjacentMaterials": [1],
                   "Protected": True},
                  {"Attribute": 2, "Role": "substrate-outer", "AdjacentMaterials": [7],
                   "Protected": True},
                  {"Attribute": 3, "Role": "interface", "AdjacentMaterials": [1, 7],
                   "Protected": True}]
        corners = [[10 * scale, 0, 0]]
        contract = {"Version": 1,
                    "VolumeMaterials": [{"Attribute": 1, "Material": "substrate"},
                                        {"Attribute": 7, "Material": "vacuum"}],
                    "BoundaryLabels": labels, "SemanticCorners": corners,
                    "ProtectedSupports": [item["Role"] for item in labels],
                    "MetricSurfaceRoles": ["air-outer", "substrate-outer"],
                    "CutSurfaceRoles": ["air-outer"], "UnmatchedPolicy": "Error",
                    "FeatureTopology": derive_feature_topology(
                        signature, directory / "boundary.csv", corners)}
        (directory / "semantic.json").write_text(json.dumps(contract) + "\n")
        names = {"Signature": "signature.csv", "Boundary": "boundary.csv", "Mask": "mask.csv",
                 "Process": "process.toml", "SemanticContract": "semantic.json",
                 "MeshRecipe": "recipe.json"}
        files = {role: {"Name": filename, "SHA256": sha256(directory / filename)}
                 for role, filename in names.items()}
        return {"Id": name, "Variants": [{"Id": "identity", "Transform": list(IDENTITY)},
                                           {"Id": "rotate-z-0.63", "Transform": list(ROTATION)}],
                "TransformComparison": {"Reference": "identity",
                                        "Transformed": "rotate-z-0.63",
                                        "MaximumRelativeInvariantError": 1e-8},
                "Source": {"Directory": name, "SignatureRole": "Signature",
                           "SignatureColumns": ["Index", "Slot", "Conductor"],
                           "Files": files}, "TestScale": scale}

    def make_suite(self, root):
        cases = [self.write_case(root, "base", 1, scale=1),
                 self.write_case(root, "subdivided", 2, scale=2),
                 self.write_case(root, "six-edge-supplemental", 6, scale=3)]
        tools = [(MESHER.name, MESHER), (STAGER.name, STAGER),
                 (AUDITOR.name, AUDITOR), (BOUNDED.name, BOUNDED)]
        manifest = {"Version": 2, "RepositoryRoot": ".",
            "Gates": {"CornerTolerance": 1e-8, "MaximumNormalFactor": 2,
                      "MinimumAchievedAspect": 1.5, "MaximumCornerAspect": 3,
                      "MinimumNoncornerAspect": 10,
                      "MaximumProtectedMeasureError": 1e-12,
                      "MinimumScaledJacobian": .01,
                      "MaximumJacobianCondition": 100, "MaximumSeconds": 10,
                      "MaximumRSSGiB": 1, "MaximumElements": 1000},
            "Tools": [{"Name": name, "Path": str(path), "SHA256": sha256(path)}
                      for name, path in tools],
            "StageToolSHA256": {
                "seed-generation": {"runtime": sha256(sys.executable),
                                    "mesher": sha256(MESHER)},
                "metric-preparation": {"runtime": sha256(sys.executable),
                                       "metric-preparer": sha256(STAGER)},
                "native-adaptation-mmg": {"runtime": sha256(STAGER),
                                          "adapter-mmg": sha256(STAGER)},
                "label-restoration": {"runtime": sha256(sys.executable),
                                      "label-restorer": sha256(STAGER)},
                "final-gmsh-publication": {"runtime": sha256(sys.executable),
                                           "publisher": sha256(STAGER)}},
            "ScalingComparisons": [
                {"Id": "subdivision", "Kind": "cad-subdivision-sensitivity",
                 "Reference": ["base", "identity"], "Compared": ["subdivided", "identity"],
                 "MaximumNormalizedDOFRatio": 1.1},
                {"Id": "features", "Kind": "feature-scaling",
                 "Reference": ["base", "identity"],
                 "Compared": ["six-edge-supplemental", "identity"],
                 "MaximumNormalizedDOFRatio": 10}], "Cases": cases}
        path = root / "suite.json"; path.write_text(json.dumps(manifest, indent=2) + "\n")
        return path, manifest

    def args(self, manifest, output, *, audit=None, preflight=False):
        return SimpleNamespace(manifest=manifest, input=[], root=output,
                               preflight_only=preflight, audit_root=audit)

    def produce_matrix(self, root, manifest_path, manifest):
        audits = root / "audits"; audits.mkdir()
        for case in manifest["Cases"]:
            directory = root / case["Id"]
            inputs = {role: item["SHA256"] for role, item in case["Source"]["Files"].items()}
            inputs_path = directory / "input-hashes.json"
            inputs_path.write_text(json.dumps(inputs))
            identity_mesh = None
            for variant in case["Variants"]:
                variant_id = variant["Id"]; stem = f"{case['Id']}--{variant_id}"
                transform = directory / f"{variant_id}-transform.json"
                transform.write_text(json.dumps(variant["Transform"]))
                seed = root / f"{stem}-seed.msh"
                metric = root / f"{stem}-metric.json"
                mmg_seed = root / f"{stem}-mmg-seed.msh"
                pins = root / f"{stem}-pins.txt"
                fixed_triangles = root / f"{stem}-fixed-triangles.txt"
                restoration_recipe = root / f"{stem}-restoration-recipe.json"
                adapted = root / f"{stem}-adapted.msh"
                restored = root / f"{stem}-restored.msh"
                mesh = root / f"{stem}.msh"
                ownership = root / f"{stem}-ownership.csv"
                stage_reports = {}
                def launch(stage, inputs, artifacts, tools, command):
                    log = root / f"{stem}-{stage}.log"
                    invocation = [sys.executable, str(BOUNDED), "--seconds", "10",
                                  "--memory-gib", "1", "--log", str(log),
                                  "--stage", stage]
                    for name, path in inputs.items():
                        invocation += ["--input", f"{name}={path}"]
                    for name, path in artifacts.items():
                        invocation += ["--artifact", f"{name}={path}"]
                    for role, path in tools.items():
                        invocation += ["--tool", f"{role}={path}"]
                    subprocess.run([*invocation, "--", *command], check=True,
                                   stdout=subprocess.DEVNULL)
                    stage_reports[stage] = log.with_suffix(".log.json")
                launch("seed-generation", {}, {"seed-mesh": seed},
                       {"runtime": sys.executable, "mesher": MESHER},
                       [sys.executable, str(MESHER), str(seed), str(transform),
                        "--scale", str(case["TestScale"])])
                launch("metric-preparation", {"seed-mesh": seed},
                       {"metric": metric, "mmg-seed": mmg_seed, "pins": pins,
                        "fixed-triangles": fixed_triangles,
                        "restoration-recipe": restoration_recipe},
                       {"runtime": sys.executable, "metric-preparer": STAGER},
                       [sys.executable, str(STAGER), "metric", str(seed), str(metric),
                        "--mmg-seed", str(mmg_seed), "--pins", str(pins),
                        "--fixed-triangles", str(fixed_triangles),
                        "--recipe", str(restoration_recipe)])
                adapter = root / "tiny-native-adapter"
                if not adapter.exists():
                    shutil.copyfile(STAGER, adapter); adapter.chmod(0o755)
                launch("native-adaptation-mmg",
                       {"mmg-seed": mmg_seed, "metric": metric, "pins": pins,
                        "fixed-triangles": fixed_triangles},
                       {"adapted-mesh": adapted},
                       {"runtime": adapter, "adapter-mmg": adapter},
                       [str(adapter), "adapt", str(mmg_seed), str(adapted),
                        "--metric", str(metric), "--pins", str(pins),
                        "--fixed-triangles", str(fixed_triangles)])
                launch("label-restoration", {"adapted-mesh": adapted,
                                              "restoration-recipe": restoration_recipe},
                       {"restored-mesh": restored},
                       {"runtime": sys.executable, "label-restorer": STAGER},
                       [sys.executable, str(STAGER), "restore", str(adapted), str(restored),
                        "--recipe", str(restoration_recipe)])
                launch("final-gmsh-publication", {"restored-mesh": restored},
                       {"candidate-mesh": mesh, "ownership-partition": ownership},
                       {"runtime": sys.executable, "publisher": STAGER},
                       [sys.executable, str(STAGER), "publish", str(restored), str(mesh),
                        "--ownership", str(ownership)])
                if identity_mesh is None: identity_mesh = mesh
                records = {}
                for kind in KINDS:
                    record = root / f"{stem}-{kind}.json"
                    produce_audit(kind, case["Id"], variant_id, mesh, inputs_path,
                                  transform, record, contract=directory / "semantic.json",
                                  recipe=directory / "recipe.json", process=directory / "process.toml",
                                  signature=directory / "signature.csv",
                                  identity_mesh=identity_mesh, stage_reports=stage_reports,
                                  command=[str(AUDITOR), kind, stem])
                    records[kind] = record
                normalize(manifest_path, case["Id"], variant_id, mesh, records,
                          audits / f"{stem}.json")
        return audits

    def test_staged_gmsh_evidence_accepts_complete_fixture_matrix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            self.assertTrue(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
            summary = json.loads((root / "out" / "summary.json").read_text())
            self.assertTrue(summary["Passed"])
            self.assertEqual(len(next(c for c in summary["Cases"]
                                      if c["Id"].startswith("six"))["VariantResults"]), 2)

    def test_mismatched_invalid_and_duplicate_mesh_content_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            evidence_path = audits / "base--identity.json"
            evidence = json.loads(evidence_path.read_text())
            evidence["Mesh"] = json.loads((audits / "subdivided--identity.json").read_text())["Mesh"]
            evidence_path.write_text(json.dumps(evidence))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            evidence_path = audits / "base--identity.json"
            evidence = json.loads(evidence_path.read_text())
            invalid = root / "invalid.msh"; invalid.write_text("not a Gmsh mesh\n")
            evidence["Mesh"] = {"Path": str(invalid), "SHA256": sha256(invalid)}
            for item in evidence["AuditRecords"]:
                record_path = Path(item["Path"])
                if not record_path.is_absolute(): record_path = audits / record_path
                record = json.loads(record_path.read_text())
                record["MeshSHA256"] = sha256(invalid)
                record_path.write_text(json.dumps(record))
                item["SHA256"] = sha256(record_path)
            evidence_path.write_text(json.dumps(evidence))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            identity = json.loads((audits / "base--identity.json").read_text())
            rotated_path = audits / "base--rotate-z-0.63.json"
            rotated = json.loads(rotated_path.read_text())
            copied = root / "copied-identity.msh"
            identity_mesh = Path(identity["Mesh"]["Path"])
            if not identity_mesh.is_absolute(): identity_mesh = audits / identity_mesh
            shutil.copyfile(identity_mesh, copied)
            rotated["Mesh"] = {"Path": str(copied), "SHA256": sha256(copied)}
            for item in rotated["AuditRecords"]:
                record_path = Path(item["Path"])
                if not record_path.is_absolute(): record_path = audits / record_path
                record = json.loads(record_path.read_text())
                record["MeshSHA256"] = sha256(copied)
                record_path.write_text(json.dumps(record))
                item["SHA256"] = sha256(record_path)
            rotated_path.write_text(json.dumps(rotated))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

    def test_duplicate_audit_content_and_identity_as_rotation_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            evidence_path = audits / "base--identity.json"
            evidence = json.loads(evidence_path.read_text())
            source = Path(evidence["AuditRecords"][0]["Path"])
            if not source.is_absolute(): source = audits / source
            copied = root / "copied-record.json"; shutil.copyfile(source, copied)
            evidence["AuditRecords"][1] = {"Path": str(copied),
                                            "SHA256": sha256(copied),
                                            "Kind": evidence["AuditRecords"][1]["Kind"]}
            evidence_path.write_text(json.dumps(evidence))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            case = manifest["Cases"][0]; directory = root / "base"
            with self.assertRaises(ValueError):
                produce_audit("variant-transform", "base", "rotate-z-0.63",
                    root / "base--identity.msh", directory / "input-hashes.json",
                    directory / "rotate-z-0.63-transform.json", root / "bad.json",
                    identity_mesh=root / "base--identity.msh")

    def test_all_audit_bindings_and_raw_as_own_audit_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            stem = "base--identity"
            records = {kind: root / f"{stem}-{kind}.json" for kind in KINDS}
            source = json.loads(records["mesh-invariants"].read_text())
            mutations = {
                "mesh": lambda record: record.__setitem__("MeshSHA256", "0" * 64),
                "input": lambda record: record["InputSHA256"].__setitem__(
                    "Signature", "0" * 64),
                "transform": lambda record: record.__setitem__("TransformSHA256", "0" * 64),
                "command": lambda record: record.__setitem__("Command", []),
                "environment": lambda record: record.__setitem__("Environment", None),
                "producer": lambda record: record["Producer"].__setitem__("SHA256", "0" * 64),
            }
            for name, mutate in mutations.items():
                with self.subTest(name=name):
                    record = copy.deepcopy(source); mutate(record)
                    path = root / f"mutated-{name}.json"
                    path.write_text(json.dumps(record))
                    selected = dict(records); selected["mesh-invariants"] = path
                    with self.assertRaises(ValueError):
                        normalize(manifest_path, "base", "identity",
                                  root / f"{stem}.msh", selected,
                                  root / f"normalized-{name}.json")
            bounded_source = json.loads(records["bounded-run"].read_text())
            for name, mutate in {
                    "stage-command": lambda value: value["BoundedStages"][
                        "native-adaptation-mmg"].__setitem__("Command", []),
                    "uninvolved-adapter": lambda value: value["BoundedStages"][
                        "native-adaptation-mmg"]["Tools"]["adapter-mmg"].__setitem__(
                            "Path", str(MESHER))}.items():
                with self.subTest(name=name):
                    record = copy.deepcopy(bounded_source); mutate(record)
                    path = root / f"mutated-{name}.json"
                    path.write_text(json.dumps(record))
                    selected = dict(records); selected["bounded-run"] = path
                    with self.assertRaises(ValueError):
                        normalize(manifest_path, "base", "identity",
                                  root / f"{stem}.msh", selected,
                                  root / f"normalized-{name}.json")
            same_record = {kind: records["bounded-run"] for kind in KINDS}
            with self.assertRaises(ValueError):
                normalize(manifest_path, "base", "identity", root / f"{stem}.msh",
                          same_record, root / "raw-as-own-audit.json")

    def test_non_native_primary_tool_must_be_first_executed_script(self):
        tools = {"runtime": sys.executable, "metric-preparer": STAGER}
        validate_tool_invocation(
            "metric-preparation", [sys.executable, "-B", str(STAGER)], tools)
        with self.assertRaisesRegex(ValueError, "executed script position"):
            validate_tool_invocation(
                "metric-preparation",
                [sys.executable, "-B", str(MESHER), str(STAGER)], tools)

    def test_uninvoked_adapter_mmg_stage_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed = root / "seed"; seed.write_text("seed")
            metric = root / "metric"; metric.write_text("metric")
            pins = root / "pins"; pins.write_text("pins")
            fixed = root / "fixed"; fixed.write_text("fixed")
            result = subprocess.run(
                [sys.executable, str(BOUNDED), "--seconds", "1", "--memory-gib", "1",
                 "--log", str(root / "bad.log"), "--stage", "native-adaptation-mmg",
                 "--input", f"mmg-seed={seed}", "--input", f"metric={metric}",
                 "--input", f"pins={pins}", "--input", f"fixed-triangles={fixed}",
                 "--artifact", f"adapted-mesh={root / 'adapted.msh'}",
                 "--tool", f"runtime={sys.executable}",
                 "--tool", f"adapter-mmg={STAGER}", "--",
                 sys.executable, str(MESHER), str(root / "unused")],
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Native adapter must be both runtime and adapter-mmg", result.stderr)

    def test_preexisting_stage_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "seed.msh"
            output.write_text("preexisting")
            result = subprocess.run(
                [sys.executable, str(BOUNDED), "--seconds", "1", "--memory-gib", "1",
                 "--log", str(root / "seed.log"), "--stage", "seed-generation",
                 "--artifact", f"seed-mesh={output}",
                 "--tool", f"runtime={sys.executable}", "--tool", f"mesher={MESHER}",
                 "--", sys.executable, str(MESHER), str(output), str(root / "missing")],
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stage output must be absent", result.stderr)

    def test_same_area_displaced_protected_support_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            reference = read_mesh(root / "base--identity-seed.msh")
            moved = read_mesh(root / "base--identity.msh")
            moved.points += np.array([0.01, 0.0, 0.0])
            contract = json.loads((root / "base/semantic.json").read_text())
            result = _protected_surface_report(reference, moved, contract)
            self.assertGreater(result["MaximumSupportVertexDistance"], 0)
            self.assertFalse(result["PlaneSupportsMatch"])

    def test_protected_footprint_accepts_refinement_and_rejects_same_area_reshape(self):
        square = np.array([[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]],
                           [[0., 0., 0.], [1., 1., 0.], [0., 1., 0.]]])
        boundary = np.array([[0., 0., 0.], [.5, 0., 0.], [1., 0., 0.],
                             [1., .5, 0.], [1., 1., 0.], [.5, 1., 0.],
                             [0., 1., 0.], [0., .5, 0.]])
        center = np.array([.5, .5, 0.])
        refined = np.array([[boundary[i], boundary[(i + 1) % len(boundary)], center]
                            for i in range(len(boundary))])
        distance, left_topology, right_topology = _footprint_boundary_comparison(
            square, refined)
        self.assertLessEqual(distance, 1e-15)
        self.assertEqual(left_topology, right_topology)
        self.assertEqual(left_topology,
                         {"Components": 1, "BoundaryLoops": 1, "Holes": 0,
                          "LoopsPerComponent": [1]})
        reshaped_points = np.array([[0., 0., 0.], [1., 0., 0.],
                                    [1.2, 1., 0.], [.2, 1., 0.]])
        reshaped = reshaped_points[[[0, 1, 2], [0, 2, 3]]]
        square_area = np.linalg.norm(np.cross(square[:, 1] - square[:, 0],
                                               square[:, 2] - square[:, 0]), axis=1).sum() / 2
        reshaped_area = np.linalg.norm(np.cross(reshaped[:, 1] - reshaped[:, 0],
                                                 reshaped[:, 2] - reshaped[:, 0]), axis=1).sum() / 2
        self.assertAlmostEqual(square_area, reshaped_area)
        self.assertGreater(_footprint_boundary_distance(square, reshaped), .1)

    def test_protected_footprint_rejects_same_vertices_equal_area_rewiring(self):
        points = np.array([[0., 0., 0.], [0., 1., 0.],
                           [1., 1., 0.], [3., 2., 0.]])
        first = points[[[0, 1, 2], [0, 2, 3]]]
        rewired = points[[[0, 1, 2], [1, 3, 2]]]
        area = lambda xyz: np.linalg.norm(
            np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1).sum() / 2
        self.assertAlmostEqual(area(first), area(rewired))
        self.assertEqual(set(map(tuple, first.reshape(-1, 3))),
                         set(map(tuple, rewired.reshape(-1, 3))))
        distance, first_topology, rewired_topology = _footprint_boundary_comparison(
            first, rewired)
        self.assertEqual(first_topology, rewired_topology)
        self.assertGreater(distance, .2)

    def test_protected_footprint_records_component_and_hole_topology(self):
        ring_points = np.array([[0., 0., 0.], [3., 0., 0.], [3., 3., 0.], [0., 3., 0.],
                                [1., 1., 0.], [2., 1., 0.], [2., 2., 0.], [1., 2., 0.]])
        ring = ring_points[[[0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                            [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]]]
        _, topology = _normalized_footprint_boundary(ring)
        self.assertEqual(topology,
                         {"Components": 1, "BoundaryLoops": 2, "Holes": 1,
                          "LoopsPerComponent": [2]})

    def test_diagonal_detector_includes_declared_maximum_and_rejects_semantic_lines(self):
        def strip(angle, length=1.0, width=.4):
            columns = int(round(length / .05)) + 1
            rotation = np.array([[math.cos(angle), -math.sin(angle)],
                                 [math.sin(angle), math.cos(angle)]])
            points = []
            for transverse in (-width / 2, 0., width / 2):
                for index in range(columns):
                    xy = rotation @ np.array([index * .05, transverse])
                    points.append([xy[0], xy[1], 0.])
            triangles = []
            for row in range(2):
                for index in range(columns - 1):
                    first = row * columns + index
                    last = (row + 1) * columns + index
                    triangles.extend([[first, first + 1, last + 1],
                                      [first, last + 1, last]])
            return np.asarray(points), np.asarray(triangles)

        points, triangles = strip(math.pi / 6)
        mesh = meshio.Mesh(points, [("triangle", triangles)],
                           cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
        arbitrary = _global_diagonal_bands(
            mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025)
        self.assertEqual(arbitrary["ShortEdgeThreshold"], .05)
        self.assertGreater(arbitrary["ShortInternalEdges"], 0)
        self.assertEqual(arbitrary["GlobalDiagonalBands"], 1)
        direction = [math.cos(math.pi / 6), math.sin(math.pi / 6), 0.]
        physical = _global_diagonal_bands(mesh, [[[0., 0., 0.], direction]], .025)
        self.assertEqual(physical["GlobalDiagonalBands"], 0)
        self.assertTrue(physical["LongShortEdgeComponents"][0][
            "AlignedWithPhysicalSegment"])

    def test_diagonal_detector_does_not_aggregate_orthogonal_supports(self):
        # Each support carries only a sub-half-diameter chain.  They meet along
        # one ridge, but combining chains from orthogonal planes would invent a
        # diagonal direction with no geometric meaning.
        points = np.array([[x, y, 0.] for y in (-.2, 0., .2)
                           for x in np.linspace(0., .2, 5)] +
                          [[.2, y, z] for z in (.2, .4)
                           for y in (-.2, 0., .2)])
        triangles = []
        for row in range(2):
            for column in range(4):
                first = row * 5 + column; last = (row + 1) * 5 + column
                triangles.extend([[first, first + 1, last + 1],
                                  [first, last + 1, last]])
        # The second plane reuses the x=.2 ridge (indices 4, 9, 14).
        rows = [[4, 9, 14], [15, 16, 17], [18, 19, 20]]
        for row in range(2):
            for column in range(2):
                first, right = rows[row][column], rows[row][column + 1]
                last, diagonal = rows[row + 1][column], rows[row + 1][column + 1]
                triangles.extend([[first, right, diagonal], [first, diagonal, last]])
        triangles = np.asarray(triangles)
        mesh = meshio.Mesh(points, [("triangle", triangles)],
                           cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
        report = _global_diagonal_bands(
            mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025)
        self.assertEqual(report["GlobalDiagonalBands"], 0)

    def test_measured_semantic_gates_reject_bound_failures(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            evidence_path = audits / "subdivided--identity.json"
            original_evidence = json.loads(evidence_path.read_text())
            item = next(value for value in original_evidence["AuditRecords"]
                        if value["Kind"] == "mesh-topology-quality")
            record_path = Path(item["Path"])
            if not record_path.is_absolute(): record_path = audits / record_path
            original_record = json.loads(record_path.read_text())
            mutations = {
                "corner-under-resolution": lambda value: value["CornerNeighborhoods"][0].__setitem__(
                    "MaximumAspect", 1000),
                "subdivision-isotropy": lambda value: value["SubdivisionNeighborhoods"][0].__setitem__(
                    "MaximumAspect", 1),
                "altered-protected-support": lambda value: value["ProtectedSurfaces"].__setitem__(
                    "MaximumRelativeMeasureError", .5),
                "unresolved-ownership": lambda value: value["OwnershipClosure"].__setitem__(
                    "Unmatched", 1),
                "short-edge-diagonal-band": lambda value: value["TraceDiagonal"].__setitem__(
                    "GlobalDiagonalBands", 1),
            }
            for index, (name, mutate) in enumerate(mutations.items()):
                with self.subTest(name=name):
                    evidence = copy.deepcopy(original_evidence)
                    record = copy.deepcopy(original_record)
                    mutate(evidence); mutate(record["Measurements"])
                    record_path.write_text(json.dumps(record))
                    record_item = next(value for value in evidence["AuditRecords"]
                                       if value["Kind"] == "mesh-topology-quality")
                    record_item["SHA256"] = sha256(record_path)
                    evidence_path.write_text(json.dumps(evidence))
                    self.assertFalse(run_manifest(self.args(
                        manifest_path, root / f"semantic-failure-{index}", audit=audits)))
            record_path.write_text(json.dumps(original_record))
            evidence_path.write_text(json.dumps(original_evidence))

    def test_normalized_measurement_tampering_and_negative_resources_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            path = audits / "base--identity.json"; evidence = json.loads(path.read_text())
            evidence["MeshQuality"]["MinimumScaledJacobian"] = .9
            path.write_text(json.dumps(evidence))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            path = audits / "base--identity.json"; evidence = json.loads(path.read_text())
            record_item = next(item for item in evidence["AuditRecords"]
                               if item["Kind"] == "bounded-run")
            record_path = Path(record_item["Path"])
            if not record_path.is_absolute(): record_path = audits / record_path
            record = json.loads(record_path.read_text())
            record["Measurements"]["Resources"]["Seconds"] = -1
            record_path.write_text(json.dumps(record))
            record_item["SHA256"] = sha256(record_path)
            evidence["Resources"]["Seconds"] = -1
            path.write_text(json.dumps(evidence))
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

    def test_manifest_rejects_identity_as_rotation_and_nonrigid_transform(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            manifest["Cases"][0]["Variants"][1]["Transform"] = IDENTITY
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                validate_manifest(manifest, manifest_path)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            manifest["Cases"][0]["Variants"][1]["Transform"][0] = 2
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                validate_manifest(manifest, manifest_path)

    def test_actual_spatial_coupon_output_uses_frozen_material_identity(self):
        mesh = os.environ.get("ACTUAL_SPATIAL_COUPON_MESH")
        stage_report = os.environ.get("ACTUAL_SPATIAL_COUPON_STAGE_REPORT")
        required = os.environ.get("REQUIRE_ACTUAL_SPATIAL_COUPON", "0") == "1"
        if not mesh or not stage_report:
            if required:
                self.fail("qualification requires a fresh mesh and seed-stage report")
            self.skipTest("actual four-edge producer qualification was not requested")
        mesh, stage_report = Path(mesh), Path(stage_report)
        source = HERE / "testdata/four-edge-9d2cb9bbb3fe"
        from mesh_stage_contract import validate_stage_report
        bounded = validate_stage_report(json.loads(stage_report.read_text()), "seed-generation")
        self.assertEqual(bounded["Artifacts"]["seed-mesh"]["SHA256"], sha256(mesh))
        self.assertEqual(Path(bounded["Tools"]["mesher"]["Path"]),
                         (HERE / "mesh_spatial_coupon.jl").resolve())
        command = bounded["Command"]
        for name in ("mesh-signature.csv", "plan-view-mask.csv", "plan-view-boundary.csv"):
            path = (source / name).resolve()
            self.assertTrue(any(str(path) == str(Path(item).resolve()) or
                                (not Path(item).is_absolute() and
                                 str(path).endswith("/" + item)) for item in command))
        contract = source / "semantic-contract.json"
        report, _ = analyze(read_mesh(mesh), json.loads(contract.read_text()),
                            require_material_names=True)
        self.assertEqual(report["PhysicalVolumeNames"], {1: "substrate", 2: "vacuum"})
        self.assertEqual(report["BoundaryAdjacency"][5001], [1])
        self.assertEqual(report["BoundaryAdjacency"][6001], [2])
        wrong = json.loads(contract.read_text())
        wrong["VolumeMaterials"][0]["Material"] = "vacuum"
        wrong["VolumeMaterials"][1]["Material"] = "substrate"
        with self.assertRaises(ValueError):
            analyze(read_mesh(mesh), wrong, require_material_names=True)

    def test_remote_verified_contracts_match_source_model_and_physical_corners(self):
        manifest = json.loads((HERE / "geometry-independence-suite.json").read_text())
        for case_id, expected_model, expected_edges in (
                ("four-edge-9d2cb9bbb3fe",
                 "spatialedgecluster_edgecount-4_9d2cb9bbb3fe", 4),
                ("ten-edge-6791f1c84123",
                 "spatialedgecluster_edgecount-10_6791f1c84123", 10)):
            with self.subTest(case=case_id):
                case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
                directory = Path(manifest["RepositoryRoot"])
                directory = (HERE / directory / case["Source"]["Directory"]).resolve()
                process = json.loads((directory / "process-library.json").read_text())
                contract = json.loads((directory / "semantic-contract.json").read_text())
                with (directory / "mesh-signature.csv").open(newline="") as stream:
                    signature = list(csv.DictReader(stream))
                with (directory / "plan-view-boundary.csv").open(newline="") as stream:
                    boundary = list(csv.DictReader(stream))
                self.assertEqual(case["InventoryStatus"], "RemoteVerified")
                self.assertEqual(process["Name"], expected_model)
                self.assertEqual(process["Models"][0]["Name"], expected_model)
                self.assertEqual(len(signature), expected_edges)
                self.assertEqual(len(process["Models"][0]["Edges"]), expected_edges)
                expected_corners = [[float(row["X"]), float(row["Y"]),
                                     float(row["Plane"])]
                                    for row in boundary if row["Class"] == "Physical"]
                self.assertEqual(contract["SemanticCorners"], expected_corners)
                materials = {item["Attribute"]: item["Material"]
                             for item in contract["VolumeMaterials"]}
                self.assertEqual(materials, {1: "substrate", 2: "vacuum"})
                pairs = {(int(row["Slot"]), int(row["Conductor"]))
                         for row in signature}
                roles = {item["Role"] for item in contract["BoundaryLabels"]}
                for slot, conductor in pairs:
                    self.assertIn(f"conductor-{conductor}-slot-{slot}-ms", roles)
                    self.assertIn(f"conductor-{conductor}-slot-{slot}-ma", roles)
                for item in contract["BoundaryLabels"]:
                    if item["Role"].endswith("-ms"):
                        self.assertEqual(item["AdjacentMaterials"], [1])
                    if item["Role"].endswith("-ma"):
                        self.assertEqual(item["AdjacentMaterials"], [2])
        producer = (HERE / "mesh_spatial_coupon.jl").read_text()
        self.assertIn('addPhysicalGroup(3, substrate_tags, 1, "substrate")', producer)
        self.assertIn('addPhysicalGroup(3, vacuum_tags, 2, "vacuum")', producer)

    def test_omitted_variant_and_checked_in_preflight_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            (audits / "base--rotate-z-0.63.json").unlink()
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        manifest_path = HERE / "geometry-independence-suite.json"
        manifest = json.loads(manifest_path.read_text())
        _, _, matrix = validate_manifest(manifest, manifest_path)
        self.assertEqual(len(matrix), 2 * len(manifest["Cases"]))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "preflight"
            self.assertTrue(run_manifest(self.args(manifest_path, output, preflight=True)))
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual({case["Id"] for case in summary["Cases"] if not case["Passed"]},
                             set())
            remote = {case["Id"]: case for case in summary["Cases"]
                      if case["Id"] in {"four-edge-9d2cb9bbb3fe",
                                        "ten-edge-6791f1c84123"}}
            self.assertEqual(remote["four-edge-9d2cb9bbb3fe"]["DiscoveredEdgeCount"], 4)
            self.assertEqual(remote["ten-edge-6791f1c84123"]["DiscoveredEdgeCount"], 10)


class SemanticContractTest(unittest.TestCase):
    def test_finite_segments_distinguish_subdivision_from_separated_collinear_features(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            signature = root / "signature.csv"
            boundary = root / "boundary.csv"
            boundary.write_text(
                "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
                "1,1,1,0,0,Physical,0,0\n"
                "1,2,1,0,0,Continuation,3,0\n")
            header = "Index,Slot,Conductor,Px,Py,Pz,Tx,Ty,Tz,S0,S1\n"
            signature.write_text(header +
                "1,0,1,.5,0,0,1,0,0,-.5,.5\n"
                "2,0,1,2.5,0,0,1,0,0,-.5,.5\n")
            separated = derive_feature_topology(signature, boundary, [[10, 0, 0]])
            self.assertEqual(separated["PhysicalFeatureCount"], 2)
            self.assertEqual(separated["CADSubdivisionCount"], 0)
            self.assertEqual(separated["BoundaryContinuationVertexCount"], 1)
            signature.write_text(header +
                "1,0,1,.5,0,0,1,0,0,-.5,.5\n"
                "2,0,1,1.5,0,0,1,0,0,-.5,.5\n")
            subdivided = derive_feature_topology(signature, boundary, [[10, 0, 0]])
            self.assertEqual(subdivided["PhysicalFeatureCount"], 1)
            self.assertEqual(subdivided["CADSubdivisionCount"], 1)
            self.assertEqual(subdivided["CADSubdivisionEndpoints"], [[1.0, 0.0, 0.0]])

    def test_vacuous_contract_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_semantic_contract({"Version": 1, "VolumeMaterials": [],
                "BoundaryLabels": [], "SemanticCorners": [], "ProtectedSupports": [],
                "MetricSurfaceRoles": [], "CutSurfaceRoles": [], "UnmatchedPolicy": "Error"})


if __name__ == "__main__":
    unittest.main()
