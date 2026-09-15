# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

from general_mesh_manifest import run_manifest, sha256, validate_manifest
from normalize_general_mesh_evidence import normalize
from semantic_mesh_contract import validate_semantic_contract
from testdata.tiny_mesh_audit_producer import produce


HERE = Path(__file__).resolve().parent
PRODUCER = HERE / "testdata" / "tiny_mesh_audit_producer.py"
IDENTITY = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
ANGLE = 0.63
ROTATION = [math.cos(ANGLE), -math.sin(ANGLE), 0, 0,
            math.sin(ANGLE), math.cos(ANGLE), 0, 0,
            0, 0, 1, 0, 0, 0, 0, 1]


class GeneralMeshManifestTest(unittest.TestCase):
    def write_case(self, root, name, edges, labels, *, features, subdivisions):
        directory = root / name
        directory.mkdir()
        signature = directory / "signature.csv"
        with signature.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Index", "Slot", "Conductor"])
            writer.writerows((i + 1, i % 2, i % 3 + 1) for i in range(edges))
        (directory / "boundary.csv").write_text("boundary\n")
        (directory / "mask.csv").write_text("mask\n")
        (directory / "process.toml").write_text('Units = "um"\n')
        (directory / "recipe.json").write_text('{"Version": 1}\n')
        contract = {
            "Version": 1,
            "VolumeMaterials": [{"Attribute": 1, "Material": "air"},
                                {"Attribute": 7, "Material": "substrate"}],
            "BoundaryLabels": labels,
            "SemanticCorners": [[0, 0, 0]],
            "ProtectedSupports": [item["Role"] for item in labels],
            "MetricSurfaceRoles": [labels[-1]["Role"]],
            "CutSurfaceRoles": [labels[0]["Role"]],
            "UnmatchedPolicy": "Error",
        }
        (directory / "semantic.json").write_text(json.dumps(contract) + "\n")
        names = {"Signature": "signature.csv", "Boundary": "boundary.csv", "Mask": "mask.csv",
                 "Process": "process.toml", "SemanticContract": "semantic.json",
                 "MeshRecipe": "recipe.json"}
        files = {role: {"Name": filename, "SHA256": sha256(directory / filename)}
                 for role, filename in names.items()}
        variants = [{"Id": "identity", "Transform": IDENTITY},
                    {"Id": "rotate-z-0.63", "Transform": ROTATION}]
        case = {"Id": name, "Variants": variants,
                "TransformComparison": {"Reference": "identity",
                                        "Transformed": "rotate-z-0.63",
                                        "MaximumRelativeInvariantError": 1e-8},
                "Source": {"Directory": name, "SignatureRole": "Signature",
                           "SignatureColumns": ["Index", "Slot", "Conductor"],
                           "Files": files},
                "TestComplexity": {"H1DOFs": features * 60 + subdivisions,
                                   "FeatureCount": features,
                                   "CADSubdivisionCount": subdivisions}}
        return case

    def make_suite(self, root):
        simple = [
            {"Attribute": 1, "Role": "outer", "AdjacentMaterials": [1], "Protected": True},
            {"Attribute": 5001, "Role": "conductor-1-slot-0", "AdjacentMaterials": [7],
             "Protected": True},
        ]
        six_labels = [
            {"Attribute": 1, "Role": "outer", "AdjacentMaterials": [1], "Protected": True},
            {"Attribute": 3100, "Role": "matching-0", "AdjacentMaterials": [1, 7],
             "Protected": True},
            {"Attribute": 3101, "Role": "matching-1", "AdjacentMaterials": [1, 7],
             "Protected": True},
        ]
        for family, material in ((5000, 1), (6000, 7)):
            for slot in (0, 1):
                for conductor in (1, 2):
                    six_labels.append({"Attribute": family + 100 * slot + conductor,
                                       "Role": f"surface-{family}-{slot}-{conductor}",
                                       "AdjacentMaterials": [material], "Protected": True})
        cases = [self.write_case(root, "base", 1, simple, features=2, subdivisions=2),
                 self.write_case(root, "subdivided", 2, simple, features=2, subdivisions=5),
                 self.write_case(root, "six-edge-supplemental", 6, six_labels,
                                 features=6, subdivisions=6)]
        manifest = {
            "Version": 2, "RepositoryRoot": ".",
            "Gates": {"CornerTolerance": 1e-8, "MaximumNormalFactor": 2,
                      "MinimumAchievedAspect": 1.5, "MinimumScaledJacobian": 0.01,
                      "MaximumJacobianCondition": 100, "MaximumSeconds": 10,
                      "MaximumRSSGiB": 1, "MaximumElements": 1000},
            "Tools": [{"Name": "tiny-producer", "Path": str(PRODUCER),
                       "SHA256": sha256(PRODUCER)}],
            "ScalingComparisons": [
                {"Id": "subdivision", "Kind": "cad-subdivision-sensitivity",
                 "Reference": ["base", "identity"], "Compared": ["subdivided", "identity"],
                 "MaximumNormalizedDOFRatio": 1.1},
                {"Id": "features", "Kind": "feature-scaling",
                 "Reference": ["base", "identity"],
                 "Compared": ["six-edge-supplemental", "identity"],
                 "MaximumNormalizedDOFRatio": 1.1},
            ], "Cases": cases}
        path = root / "suite.json"
        path.write_text(json.dumps(manifest, indent=2) + "\n")
        return path, manifest

    def args(self, manifest, output, *, audit=None, preflight=False):
        return SimpleNamespace(manifest=manifest, input=[], root=output,
                               preflight_only=preflight, audit_root=audit)

    def produce_matrix(self, root, manifest_path, manifest):
        audits = root / "audits"
        audits.mkdir()
        for case in manifest["Cases"]:
            contract = root / case["Id"] / "semantic.json"
            complexity = case["TestComplexity"]
            for variant in ("identity", "rotate-z-0.63"):
                stem = f"{case['Id']}--{variant}"
                raw, mesh = root / f"{stem}-raw.json", root / f"{stem}.msh"
                produce(contract, case["Id"], variant, raw, mesh,
                        dofs=complexity["H1DOFs"], features=complexity["FeatureCount"],
                        subdivisions=complexity["CADSubdivisionCount"])
                normalize(manifest_path, case["Id"], variant, raw, mesh, [raw],
                          audits / f"{stem}.json")
        return audits

    def test_producer_backed_matrix_accepts_six_edge_multiconductor_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            self.assertTrue(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
            summary = json.loads((root / "out" / "summary.json").read_text())
            self.assertTrue(summary["Passed"])
            six = next(case for case in summary["Cases"] if case["Id"].startswith("six"))
            self.assertEqual(six["DiscoveredEdgeCount"], 6)
            self.assertEqual(six["DiscoveredConductors"], [1, 2, 3])
            self.assertEqual(len(six["VariantResults"]), 2)

    def test_wrong_case_variant_input_process_recipe_tool_mesh_and_stale_audit_fail(self):
        mutations = {
            "case": lambda e: e.__setitem__("CaseId", "wrong"),
            "variant": lambda e: e.__setitem__("Variant", "wrong"),
            "transform": lambda e: e.__setitem__("TransformSHA256", "0" * 64),
            "input": lambda e: e["InputSHA256"].__setitem__("Signature", "0" * 64),
            "process": lambda e: e.__setitem__("ProcessSHA256", "0" * 64),
            "recipe": lambda e: e.__setitem__("RecipeSHA256", "0" * 64),
            "tool": lambda e: e.__setitem__("ToolSHA256", {"wrong": "0" * 64}),
            "mesh": lambda e: e["Mesh"].__setitem__("SHA256", "0" * 64),
            "empty-materials": lambda e: e.__setitem__("ActualVolumeMaterials", []),
            "empty-labels": lambda e: e.__setitem__("ActualBoundaryAttributes", []),
            "empty-adjacency": lambda e: e.__setitem__("ActualAdjacency", {}),
            "empty-corners": lambda e: e.__setitem__("ActualSemanticCorners", []),
            "empty-protected": lambda e: e["ProtectedSurfaces"].__setitem__("Actual", []),
        }
        for name, mutation in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                manifest_path, manifest = self.make_suite(root)
                audits = self.produce_matrix(root, manifest_path, manifest)
                path = audits / "base--identity.json"
                evidence = json.loads(path.read_text())
                mutation(evidence)
                path.write_text(json.dumps(evidence) + "\n")
                self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            raw = root / "base--identity-raw.json"
            raw.write_text(raw.read_text() + " ")
            self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

    def test_empty_contracts_and_omitted_variant_or_supplemental_fail_closed(self):
        invalid_fields = {"VolumeMaterials": [], "BoundaryLabels": [],
                          "SemanticCorners": [], "ProtectedSupports": []}
        for field, value in invalid_fields.items():
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                manifest_path, manifest = self.make_suite(root)
                path = root / "base" / "semantic.json"
                contract = json.loads(path.read_text())
                contract[field] = value
                path.write_text(json.dumps(contract))
                manifest["Cases"][0]["Source"]["Files"]["SemanticContract"]["SHA256"] = sha256(path)
                manifest_path.write_text(json.dumps(manifest))
                self.assertFalse(run_manifest(
                    self.args(manifest_path, root / "out", preflight=True)))
                summary = json.loads((root / "out" / "summary.json").read_text())
                failed = next(case for case in summary["Cases"] if case["Id"] == "base")
                self.assertIn(field, failed["Error"])
        for missing in ("base--rotate-z-0.63.json", "six-edge-supplemental--identity.json"):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                manifest_path, manifest = self.make_suite(root)
                audits = self.produce_matrix(root, manifest_path, manifest)
                (audits / missing).unlink()
                self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

    def test_negative_resources_bad_quality_and_scaling_fail(self):
        mutations = {
            "negative-seconds": ("Resources", "Seconds", -1),
            "negative-memory": ("Resources", "PeakRSSGiB", -1),
            "negative-elements": ("Resources", "Elements", -1),
            "quality": ("MeshQuality", "MinimumScaledJacobian", 0),
            "jacobian": ("MeshQuality", "MaximumJacobianCondition", 1000),
            "negative-subdivisions": ("Complexity", "CADSubdivisionCount", -1),
            "scaling": ("Complexity", "H1DOFs", 10000),
        }
        for name, (section, key, value) in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                manifest_path, manifest = self.make_suite(root)
                audits = self.produce_matrix(root, manifest_path, manifest)
                target = "six-edge-supplemental--identity.json" if name == "scaling" else "base--identity.json"
                path = audits / target
                evidence = json.loads(path.read_text())
                evidence[section][key] = value
                path.write_text(json.dumps(evidence) + "\n")
                self.assertFalse(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

    def test_checked_in_manifest_schema_and_real_preflight_fail_only_unavailable_inputs(self):
        manifest_path = HERE / "geometry-independence-suite.json"
        manifest = json.loads(manifest_path.read_text())
        _, _, matrix = validate_manifest(manifest, manifest_path)
        self.assertEqual(len(matrix), 2 * len(manifest["Cases"]))
        ordinary = {case["Id"] for case in manifest["Cases"]}
        self.assertTrue({"hole", "rounded-strip", "opposed-layers", "concave-multislot"} <= ordinary)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "preflight"
            self.assertFalse(run_manifest(self.args(manifest_path, output, preflight=True)))
            summary = json.loads((output / "summary.json").read_text())
            failed = {case["Id"] for case in summary["Cases"] if not case["Passed"]}
            self.assertEqual(failed, {"four-edge-9d2cb9bbb3fe", "ten-edge-6791f1c84123"})


class SemanticContractTest(unittest.TestCase):
    def test_vacuous_contract_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_semantic_contract({"Version": 1, "VolumeMaterials": [],
                                        "BoundaryLabels": [], "SemanticCorners": [],
                                        "ProtectedSupports": [], "MetricSurfaceRoles": [],
                                        "CutSurfaceRoles": [], "UnmatchedPolicy": "Error"})


if __name__ == "__main__":
    unittest.main()
