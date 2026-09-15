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

from general_mesh_audit_producer import KINDS, produce as produce_audit
from general_mesh_manifest import run_manifest, sha256, validate_manifest
from normalize_general_mesh_evidence import normalize
from semantic_mesh_contract import validate_semantic_contract


HERE = Path(__file__).resolve().parent
MESHER = HERE / "testdata" / "tiny_mesh_audit_producer.py"
AUDITOR = HERE / "general_mesh_audit_producer.py"
BOUNDED = HERE / "run_bounded_mesher.py"
ADAPTOR = HERE / "audit_edge_metric_mesh.py"
MMG_FIXTURE = HERE / "mesh_array_io.py"
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
                writer.writerow([i + 1, slot, conductor, i * scale, y * scale, 0,
                                 1, 0, 0, -i * scale, (2 - i) * scale])
        (directory / "boundary.csv").write_text("boundary\n")
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
        contract = {"Version": 1,
                    "VolumeMaterials": [{"Attribute": 1, "Material": "air"},
                                        {"Attribute": 7, "Material": "substrate"}],
                    "BoundaryLabels": labels, "SemanticCorners": [[0, 0, 0]],
                    "ProtectedSupports": [item["Role"] for item in labels],
                    "MetricSurfaceRoles": ["air-outer", "substrate-outer"],
                    "CutSurfaceRoles": ["air-outer"], "UnmatchedPolicy": "Error"}
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
        tools = [(MESHER.name, MESHER), (AUDITOR.name, AUDITOR), (BOUNDED.name, BOUNDED)]
        manifest = {"Version": 2, "RepositoryRoot": ".",
            "Gates": {"CornerTolerance": 1e-8, "MaximumNormalFactor": 2,
                      "MinimumAchievedAspect": 1.5, "MinimumScaledJacobian": .01,
                      "MaximumJacobianCondition": 100, "MaximumSeconds": 10,
                      "MaximumRSSGiB": 1, "MaximumElements": 1000},
            "Tools": [{"Name": name, "Path": str(path), "SHA256": sha256(path)}
                      for name, path in tools],
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
                mesh, log = root / f"{stem}.msh", root / f"{stem}.log"
                command = [sys.executable, str(MESHER), str(mesh), str(transform),
                           "--scale", str(case["TestScale"])]
                subprocess.run([sys.executable, str(BOUNDED), "--seconds", "10",
                                "--memory-gib", "1", "--log", str(log),
                                "--artifact", str(mesh), "--require-complete-toolchain",
                                "--tool", f"runtime={sys.executable}",
                                "--tool", f"mesher={MESHER}",
                                "--tool", f"adaptor={ADAPTOR}",
                                "--tool", f"mmg={MMG_FIXTURE}",
                                "--", *command], check=True,
                               stdout=subprocess.DEVNULL)
                if identity_mesh is None: identity_mesh = mesh
                records = {}
                for kind in KINDS:
                    record = root / f"{stem}-{kind}.json"
                    produce_audit(kind, case["Id"], variant_id, mesh, inputs_path,
                                  transform, record, contract=directory / "semantic.json",
                                  recipe=directory / "recipe.json", process=directory / "process.toml",
                                  signature=directory / "signature.csv",
                                  identity_mesh=identity_mesh,
                                  bounded_report=log.with_suffix(".log.json"),
                                  command=[str(AUDITOR), kind, stem])
                    records[kind] = record
                normalize(manifest_path, case["Id"], variant_id, mesh, records,
                          audits / f"{stem}.json")
        return audits

    def test_real_gmsh_producer_chain_accepts_complete_matrix(self):
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
                    "launcher-command": lambda value: value["BoundedLauncher"].__setitem__(
                        "Command", []),
                    "toolchain-digest": lambda value: value["BoundedLauncher"]["Toolchain"][
                        "mmg"].__setitem__("SHA256", "0" * 64)}.items():
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
                pairs = {(int(row["Slot"]), int(row["Conductor"]))
                         for row in signature}
                roles = {item["Role"] for item in contract["BoundaryLabels"]}
                for slot, conductor in pairs:
                    self.assertIn(f"conductor-{conductor}-slot-{slot}-ms", roles)
                    self.assertIn(f"conductor-{conductor}-slot-{slot}-ma", roles)

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
    def test_vacuous_contract_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_semantic_contract({"Version": 1, "VolumeMaterials": [],
                "BoundaryLabels": [], "SemanticCorners": [], "ProtectedSupports": [],
                "MetricSurfaceRoles": [], "CutSurfaceRoles": [], "UnmatchedPolicy": "Error"})


if __name__ == "__main__":
    unittest.main()
