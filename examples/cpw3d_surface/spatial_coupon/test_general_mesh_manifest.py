# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import atexit
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
from general_mesh_audit_producer import (KINDS, VARIANT_AUDITS_KIND,
                                         _footprint_boundary_comparison,
                                         _footprint_boundary_distance,
                                         _global_diagonal_bands,
                                         _normalized_footprint_boundary,
                                         _ownership_report, _planar_diameter,
                                         _protected_surface_report,
                                         produce as produce_audit, produce_variant_audits)
from canonical_mesh_build import (CANONICAL_ARTIFACT_ROLES, build_record,
                                  build_record_from_stage_reports)
from audit_edge_metric_mesh import (ANISOTROPY_GATE_APPLIED, ANISOTROPY_GATE_NOT_APPLICABLE,
                                    LAYER_ADJACENT_BAND_RULE)
from general_mesh_manifest import (_physical_comparison_failures,
                                   _validate_source_transformation, audit_manifest_evidence, case_gates,
                                   layer_covered_band, run_manifest, sha256, validate_manifest,
                                   validate_production_recipe, validate_production_recipe_commands)
from mesh_array_io import read_mesh
from mesh_stage_contract import (CANONICAL_STAGE_ORDER, validate_stage_report,
                                 validate_tool_invocation)
from normalize_general_mesh_evidence import normalize
from semantic_mesh_contract import (derive_feature_topology, load_semantic_contract,
                                    validate_semantic_contract)
from testdata.build_tiny_native_adapter import build as build_native_fixture, compiler


HERE = Path(__file__).resolve().parent
MESHER = HERE / "testdata" / "tiny_mesh_audit_producer.py"
AUDITOR = HERE / "general_mesh_audit_producer.py"
BOUNDED = HERE / "run_bounded_mesher.py"
STAGER = HERE / "testdata" / "tiny_mesh_stage.py"
SCRIPT_ADAPTER = HERE / "testdata" / "tiny_native_adapter.py"
WRAPPER = HERE / "run_native_mmg_adaptation.py"
# Seed-side required-region gates (decision 30): the seed and label-restoration
# commands carry them with equal values (the fixture manifest's corner gate is 3).
SEED_QUALITY_GATE_OPTIONS = ("--maximum-corner-aspect", "3", "--minimum-scaled-jacobian", ".01",
                             "--maximum-jacobian-condition", "1000",
                             "--maximum-quality-displacement-over-normal", ".75")
# The native fixture adapter links a fixture libmmg3d through rpath so the
# wrapper's runtime-library resolution is exercised; it needs a C compiler.
_NATIVE_FIXTURE_DIRECTORY = tempfile.TemporaryDirectory(prefix="tiny-native-adapter-")
atexit.register(_NATIVE_FIXTURE_DIRECTORY.cleanup)
if compiler() is not None:
    ADAPTER, MMG_LIBRARY = build_native_fixture(_NATIVE_FIXTURE_DIRECTORY.name)
else:
    ADAPTER, MMG_LIBRARY = None, None


def require_native_fixture():
    if ADAPTER is None:
        raise unittest.SkipTest("native fixture adapter requires a C compiler (cc)")
TRANSFORMER = HERE / "transform_coupon_source_contract.py"
IDENTITY = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
# Fixture corner-isotropy prescription shared by the seed and metric stages.
NORMAL_SIZE = 0.1
CORNER_ISOTROPY_RADIUS = 0.5
# Fixture trace basis: immutable roles, their stage input names/options and the
# dimensionless TraceBasisSizeRatio both stages pass.
TRACE_BASIS_FILES = {"BasisContract": "basis-contract.json", "TraceVertices": "trace-vertices.csv",
                     "TraceTriangles": "trace-triangles.csv", "ProcessLibrary": "process-library.json"}
TRACE_BASIS_BINDINGS = {"BasisContract": ("source-basis-contract", "--trace-basis-contract"),
                        "TraceVertices": ("source-trace-vertices", "--trace-vertices"),
                        "TraceTriangles": ("source-trace-triangles", "--trace-triangles"),
                        "ProcessLibrary": ("source-process-library", "--process-library")}
TRACE_BASIS_SIZE_RATIO = 1.0
ANGLE = 0.63
ROTATION = [math.cos(ANGLE), -math.sin(ANGLE), 0, 0,
            math.sin(ANGLE), math.cos(ANGLE), 0, 0,
            0, 0, 1, 0, 0, 0, 0, 1]


class GeneralMeshManifestTest(unittest.TestCase):
    def write_case(self, root, name, edges, *, scale, retained_etch=False, trace_basis=False):
        directory = root / name; directory.mkdir()
        signature = directory / "signature.csv"
        with signature.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Index", "Slot", "Conductor", "Px", "Py", "Pz",
                             "Gx", "Gy", "Gz", "Tx", "Ty", "Tz", "Nz", "S0", "S1"])
            for i in range(edges):
                if name == "subdivided":
                    slot, conductor, y = 0, 1, 0
                else:
                    slot, conductor, y = i % 2, i // 2 + 1, i
                writer.writerow([i + 1, slot, conductor, (i + .5) * scale,
                                 y * scale, 0, 0, 1, 0, 1, 0, 0, 1,
                                 -.5 * scale, .5 * scale])
        (directory / "boundary.csv").write_text(
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
            f"1,1,1,0,0,Physical,{10 * scale},0\n"
            f"1,2,1,0,0,Continuation,{2 * scale},0\n")
        (directory / "mask.csv").write_text(
            "Facet,Vertex,Conductor,Plane,X,Y\n"
            f"1,1,1,0,0,0\n1,2,1,0,{10 * scale},0\n1,3,1,0,0,{10 * scale}\n")
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
        source_extra = {"EtchFootprint": "producer-default"}
        if retained_etch:
            # A device etch footprint bound as an immutable source (fixture content).
            (directory / "retained-etch.csv").write_text("Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n")
            names["RetainedEtch"] = "retained-etch.csv"
            source_extra = {}
        if trace_basis:
            # A trace basis bound as immutable sources (fixture content, hashed only).
            for role, filename in TRACE_BASIS_FILES.items():
                (directory / filename).write_text(f"{role} fixture of {name}\n")
                names[role] = filename
        files = {role: {"Name": filename, "SHA256": sha256(directory / filename)}
                 for role, filename in names.items()}
        return {"Id": name, "Variants": [{"Id": "identity", "Transform": list(IDENTITY)},
                                           {"Id": "rotate-z-0.63", "Transform": list(ROTATION)}],
                "TransformComparison": {"Reference": "identity",
                    "Transformed": "rotate-z-0.63", "MaximumRelativeVolumeError": 1e-8,
                    "MaximumRelativeSurfaceMeasureError": 1e-8,
                    "MaximumProtectedSupportHausdorff": 1e-8,
                    "MaximumProtectedMeasureError": 1e-8,
                    "MaximumQualityDistributionRelativeError": 1e-8,
                    "MaximumAnisotropyRelativeError": 1e-8,
                    "MaximumComplexityRatio": 1.01},
                "Source": {"Directory": name, "SignatureRole": "Signature",
                           "SignatureColumns": ["Index", "Slot", "Conductor"],
                           "Files": files, **source_extra}, "TestScale": scale}

    def make_suite(self, root):
        require_native_fixture()
        cases = [self.write_case(root, "base", 1, scale=1),
                 self.write_case(root, "subdivided", 2, scale=2, retained_etch=True),
                 self.write_case(root, "six-edge-supplemental", 6, scale=3, trace_basis=True)]
        tools = [(MESHER.name, MESHER), (STAGER.name, STAGER),
                 (TRANSFORMER.name, TRANSFORMER),
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
                "canonical-source-validation": {"runtime": sha256(sys.executable),
                    "source-validator": sha256(TRANSFORMER)},
                "seed-generation": {"runtime": sha256(sys.executable),
                                    "mesher": sha256(MESHER)},
                "metric-preparation": {"runtime": sha256(sys.executable),
                                       "metric-preparer": sha256(STAGER)},
                "native-adaptation-mmg": {"runtime": sha256(sys.executable),
                    "adaptation-wrapper": sha256(WRAPPER), "adapter-mmg": sha256(ADAPTER),
                    "mmg-library": sha256(MMG_LIBRARY)},
                "label-restoration": {"runtime": sha256(sys.executable),
                                      "label-restorer": sha256(STAGER)},
                "canonical-gmsh-publication": {"runtime": sha256(sys.executable),
                                               "publisher": sha256(STAGER)},
                "proper-rigid-publication": {"runtime": sha256(sys.executable),
                    "rigid-publisher": sha256(STAGER),
                    "ownership-runtime": sha256(sys.executable),
                    "ownership-auditor": sha256(STAGER)}},
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

    def produce_matrix(self, root, manifest_path, manifest, compare_consolidated=False):
        audits = root / "audits"; audits.mkdir()
        for case in manifest["Cases"]:
            directory = root / case["Id"]
            inputs = {role: item["SHA256"] for role, item in case["Source"]["Files"].items()}
            inputs_path = directory / "input-hashes.json"; inputs_path.write_text(json.dumps(inputs))
            canonical_stem = f"{case['Id']}--canonical"
            identity_transform = directory / "canonical-identity-transform.json"
            identity_transform.write_text(json.dumps(IDENTITY))
            seed = root / f"{canonical_stem}-seed.msh"
            metric = root / f"{canonical_stem}-metric.json"
            mmg_seed = root / f"{canonical_stem}-mmg-seed.msh"
            pins = root / f"{canonical_stem}-pins.txt"
            fixed = root / f"{canonical_stem}-fixed.txt"
            required = root / f"{canonical_stem}-required.txt"
            recipe = root / f"{canonical_stem}-recipe.json"
            canonical_semantic = root / f"{canonical_stem}-semantic.json"
            canonical_supports = root / f"{canonical_stem}-supports.json"
            adapted = root / f"{canonical_stem}-adapted.msh"
            adaptation_receipt = root / f"{canonical_stem}-adaptation.json"
            local_restored = root / f"{canonical_stem}-local-restored.msh"
            restored = root / f"{canonical_stem}-restored.msh"
            canonical_mesh = root / f"{canonical_stem}.msh"
            canonical_ownership = root / f"{canonical_stem}-ownership.csv"
            canonical_quadrature = root / f"{canonical_stem}-ownership.quadrature.csv"
            canonical_reports = {}
            def launch(stem, reports, stage, stage_inputs, artifacts, tools, command):
                log = root / f"{stem}-{stage}.log"
                invocation = [sys.executable, str(BOUNDED), "--seconds", "10",
                              "--memory-gib", "1", "--log", str(log), "--stage", stage]
                for name, path in stage_inputs.items(): invocation += ["--input", f"{name}={path}"]
                for name, path in artifacts.items(): invocation += ["--artifact", f"{name}={path}"]
                for role, path in tools.items(): invocation += ["--tool", f"{role}={path}"]
                subprocess.run([*invocation, "--", *command], check=True,
                               stdout=subprocess.DEVNULL)
                reports[stage] = log.with_suffix(".log.json")
            launch(canonical_stem, canonical_reports, "canonical-source-validation", {
                "source-semantic-contract": directory / "semantic.json",
                "source-signature": directory / "signature.csv",
                "source-boundary": directory / "boundary.csv",
                "source-mask": directory / "mask.csv", "canonical-transform": identity_transform},
                {"canonical-semantic-contract": canonical_semantic,
                 "canonical-supports": canonical_supports},
                {"runtime": sys.executable, "source-validator": TRANSFORMER},
                [sys.executable, str(TRANSFORMER), str(directory), str(identity_transform),
                 str(canonical_semantic), str(canonical_supports), "--semantic-input",
                 str(directory / "semantic.json"), "--signature", str(directory / "signature.csv"),
                 "--boundary", str(directory / "boundary.csv"), "--mask", str(directory / "mask.csv")])
            census = root / f"{canonical_stem}-corner-census.json"
            etch = (directory / case["Source"]["Files"]["RetainedEtch"]["Name"]
                    if "RetainedEtch" in case["Source"]["Files"] else None)
            # The trace basis is bound to the seed and the metric stage when declared.
            basis_inputs, basis_options = {}, []
            if "BasisContract" in case["Source"]["Files"]:
                for role, (name, option) in TRACE_BASIS_BINDINGS.items():
                    path = directory / case["Source"]["Files"][role]["Name"]
                    basis_inputs[name] = path; basis_options += [option, str(path)]
                basis_options += ["--trace-basis-size-ratio", str(TRACE_BASIS_SIZE_RATIO)]
            launch(canonical_stem, canonical_reports, "seed-generation",
                {"source-signature": directory / "signature.csv",
                 "source-boundary": directory / "boundary.csv",
                 "source-mask": directory / "mask.csv",
                 "canonical-semantic-contract": canonical_semantic,
                 **({"source-retained-etch": etch} if etch is not None else {}),
                 **basis_inputs},
                {"seed-mesh": seed, "seed-corner-census": census},
                {"runtime": sys.executable, "mesher": MESHER},
                [sys.executable, str(MESHER), str(seed), str(identity_transform),
                 "--scale", str(case["TestScale"]), "--signature", str(directory / "signature.csv"),
                 "--mask", str(directory / "mask.csv"), "--boundary", str(directory / "boundary.csv"),
                 "--semantic-contract", str(canonical_semantic),
                 "--corner-isotropy-radius", str(CORNER_ISOTROPY_RADIUS),
                 "--lc-fine", str(NORMAL_SIZE), "--corner-census", str(census),
                 *(["--etch-boundary", str(etch)] if etch is not None else []),
                 *basis_options, *SEED_QUALITY_GATE_OPTIONS])
            launch(canonical_stem, canonical_reports, "metric-preparation",
                {"seed-mesh": seed, "seed-corner-census": census,
                 "canonical-semantic-contract": canonical_semantic,
                 "canonical-supports": canonical_supports, **basis_inputs},
                {"metric": metric, "mmg-seed": mmg_seed, "pins": pins,
                 "fixed-triangles": fixed, "required-tetrahedra": required,
                 "restoration-recipe": recipe},
                {"runtime": sys.executable, "metric-preparer": STAGER},
                [sys.executable, str(STAGER), "metric", str(seed), str(metric),
                 "--mmg-seed", str(mmg_seed), "--pins", str(pins), "--fixed-triangles", str(fixed),
                 "--required-tetrahedra", str(required), "--recipe", str(recipe), "--semantic-contract", str(canonical_semantic),
                 "--transformed-supports", str(canonical_supports),
                 "--seed-census", str(census),
                 "--normal", str(NORMAL_SIZE), "--tangent", str(CORNER_ISOTROPY_RADIUS),
                 *basis_options])
            launch(canonical_stem, canonical_reports, "native-adaptation-mmg",
                {"mmg-seed": mmg_seed, "metric": metric, "pins": pins,
                 "fixed-triangles": fixed, "required-tetrahedra": required,
                 "restoration-recipe": recipe},
                {"adapted-mesh": adapted, "adaptation-receipt": adaptation_receipt},
                {"runtime": sys.executable, "adaptation-wrapper": WRAPPER, "adapter-mmg": ADAPTER,
                 "mmg-library": MMG_LIBRARY},
                [sys.executable, str(WRAPPER), str(mmg_seed), str(metric), str(pins), str(recipe),
                 str(adapted), str(adaptation_receipt), "--adapter", str(ADAPTER),
                 "--mmg-library", str(MMG_LIBRARY), "--hmin", ".1",
                 "--hgrad", "1.3", "--fixed-triangles", str(fixed),
                 "--required-tetrahedra", str(required)])
            launch(canonical_stem, canonical_reports, "label-restoration",
                {"adapted-mesh": adapted, "restoration-recipe": recipe},
                {"source-local-restored-mesh": local_restored, "restored-mesh": restored},
                {"runtime": sys.executable, "label-restorer": STAGER},
                [sys.executable, str(STAGER), "restore", str(adapted), str(restored),
                 "--recipe", str(recipe), "--source-local-output", str(local_restored),
                 *SEED_QUALITY_GATE_OPTIONS])
            launch(canonical_stem, canonical_reports, "canonical-gmsh-publication",
                {"restored-mesh": restored, "source-process": directory / "process.toml",
                 "source-signature": directory / "signature.csv",
                 "source-boundary": directory / "boundary.csv"},
                {"canonical-candidate-mesh": canonical_mesh,
                 "canonical-ownership-partition": canonical_ownership,
                 "canonical-ownership-quadrature-partition": canonical_quadrature},
                {"runtime": sys.executable, "publisher": STAGER},
                [sys.executable, str(STAGER), "publish", str(restored), str(canonical_mesh),
                 "--ownership", str(canonical_ownership),
                 "--ownership-quadrature", str(canonical_quadrature),
                 "--process", str(directory / "process.toml"),
                 "--signature", str(directory / "signature.csv"),
                 "--boundary", str(directory / "boundary.csv")])
            canonical_tools = {f"{stage}/{role}": digest for stage in CANONICAL_STAGE_ORDER
                               for role, digest in manifest["StageToolSHA256"][stage].items()}
            canonical_record = root / f"{canonical_stem}-build.json"
            canonical_record.write_text(json.dumps(build_record_from_stage_reports(
                inputs, manifest["Gates"], canonical_tools, canonical_reports)) + "\n")
            for variant in case["Variants"]:
                variant_id = variant["Id"]; stem = f"{case['Id']}--{variant_id}"
                transform = directory / f"{variant_id}-transform.json"
                transform.write_text(json.dumps(variant["Transform"]))
                mesh = root / f"{stem}.msh"; transformed_semantic = root / f"{stem}-semantic.json"
                transformed_supports = root / f"{stem}-supports.json"
                ownership = root / f"{stem}-ownership.csv"
                quadrature = root / f"{stem}-ownership.quadrature.csv"
                receipt = root / f"{stem}-transform-receipt.json"
                stage_reports = dict(canonical_reports)
                launch(stem, stage_reports, "proper-rigid-publication", {
                    "canonical-candidate-mesh": canonical_mesh,
                    "canonical-build-record": canonical_record, "placement-transform": transform,
                    "source-semantic-contract": directory / "semantic.json",
                    "source-signature": directory / "signature.csv",
                    "source-boundary": directory / "boundary.csv",
                    "source-mask": directory / "mask.csv", "source-process": directory / "process.toml"},
                    {"candidate-mesh": mesh, "transformed-semantic-contract": transformed_semantic,
                     "transformed-supports": transformed_supports, "ownership-partition": ownership,
                     "ownership-quadrature-partition": quadrature, "transform-receipt": receipt},
                    {"runtime": sys.executable, "rigid-publisher": STAGER,
                     "ownership-runtime": sys.executable, "ownership-auditor": STAGER},
                    [sys.executable, str(STAGER), "rigid", str(canonical_mesh), str(mesh),
                     "--transform", str(transform),
                     "--semantic-input", str(directory / "semantic.json"), "--signature",
                     str(directory / "signature.csv"), "--boundary", str(directory / "boundary.csv"),
                     "--mask", str(directory / "mask.csv"), "--process", str(directory / "process.toml"),
                     "--canonical-build-record", str(canonical_record),
                     "--transformed-semantic", str(transformed_semantic),
                     "--transformed-supports", str(transformed_supports), "--receipt", str(receipt),
                     "--ownership", str(ownership), "--ownership-quadrature", str(quadrature),
                     "--ownership-runtime", sys.executable, "--ownership-auditor", str(STAGER)])
                records = {}
                for kind in KINDS:
                    record = root / f"{stem}-{kind}.json"
                    produce_audit(kind, case["Id"], variant_id, mesh, inputs_path,
                        transform, record, contract=directory / "semantic.json",
                        recipe=directory / "recipe.json", process=directory / "process.toml",
                        signature=directory / "signature.csv", identity_mesh=canonical_mesh,
                        identity_seed_mesh=seed, stage_reports=stage_reports,
                        command=[str(AUDITOR), kind, stem])
                    records[kind] = record
                if compare_consolidated:
                    # One process, the mesh read once, the same producer functions: every
                    # record equals its standalone twin except the recorded Command.
                    consolidated_command = [str(AUDITOR), VARIANT_AUDITS_KIND, stem]
                    consolidated = {kind: root / f"{stem}-consolidated-{kind}.json" for kind in KINDS}
                    produce_variant_audits(case["Id"], variant_id, mesh, inputs_path, transform,
                        consolidated, contract=directory / "semantic.json",
                        recipe=directory / "recipe.json", process=directory / "process.toml",
                        signature=directory / "signature.csv", identity_mesh=canonical_mesh,
                        identity_seed_mesh=seed, stage_reports=stage_reports,
                        command=consolidated_command)
                    for kind in KINDS:
                        standalone = json.loads(records[kind].read_text())
                        merged = json.loads(consolidated[kind].read_text())
                        self.assertEqual(merged.pop("Command"), consolidated_command)
                        self.assertEqual(standalone.pop("Command"), [str(AUDITOR), kind, stem])
                        self.assertEqual(merged, standalone, kind)
                    with self.assertRaisesRegex(ValueError, "must be fresh"):
                        produce_variant_audits(case["Id"], variant_id, mesh, inputs_path,
                                               transform, consolidated, stage_reports=stage_reports)
                    with self.assertRaisesRegex(ValueError, "per audit kind"):
                        produce_variant_audits(case["Id"], variant_id, mesh, inputs_path,
                                               transform, {"bounded-run": root / "one.json"})
                normalize(manifest_path, case["Id"], variant_id, mesh, records,
                          audits / f"{stem}.json")
        return audits

    def test_variant_audits_match_standalone_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest, compare_consolidated=True)
            self.assertTrue(run_manifest(self.args(manifest_path, root / "out", audit=audits)))

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
                stage_contract = __import__("mesh_stage_contract")
                reports = {stage: root / f"base--canonical-{stage}.log.json"
                           for stage in stage_contract.CANONICAL_STAGE_ORDER}
                reports["proper-rigid-publication"] = (
                    root / "base--identity-proper-rigid-publication.log.json")
                produce_audit("variant-transform", "base", "rotate-z-0.63",
                    root / "base--identity.msh", directory / "input-hashes.json",
                    directory / "rotate-z-0.63-transform.json", root / "bad.json",
                    contract=directory / "semantic.json",
                    identity_mesh=root / "base--canonical.msh",
                    identity_seed_mesh=root / "base--canonical-seed.msh",
                    stage_reports=reports)

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

    def test_metric_recipe_must_consume_exact_reconstructed_supports(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            case = manifest["Cases"][0]
            source_paths = {role: root / case["Id"] / item["Name"]
                            for role, item in case["Source"]["Files"].items()}
            stage_contract = __import__("mesh_stage_contract")
            reports = {stage: json.loads((root / f"base--canonical-{stage}.log.json").read_text())
                       for stage in stage_contract.CANONICAL_STAGE_ORDER}
            reports["proper-rigid-publication"] = json.loads(
                (root / "base--identity-proper-rigid-publication.log.json").read_text())
            binding = {"Transform": IDENTITY,
                       "InputSHA256": {role: item["SHA256"]
                                       for role, item in case["Source"]["Files"].items()}}
            _validate_source_transformation(reports, binding, source_paths)
            recipe_path = Path(reports["metric-preparation"]["Artifacts"]
                               ["restoration-recipe"]["Path"])
            original = json.loads(recipe_path.read_text())
            mutations = {
                "omitted": lambda value: value.pop("TransformedSupports"),
                "alternate-path": lambda value: value.__setitem__(
                    "TransformedSupportsArtifact", str(root / "alternate.json")),
                "tampered": lambda value: value["TransformedSupports"]["Edges"][0]
                    ["Point"].__setitem__(0, 123.0),
            }
            for name, mutate in mutations.items():
                with self.subTest(name=name):
                    candidate = copy.deepcopy(original); mutate(candidate)
                    recipe_path.write_text(json.dumps(candidate))
                    with self.assertRaisesRegex(ValueError, "did not consume"):
                        _validate_source_transformation(reports, binding, source_paths)
            recipe_path.write_text(json.dumps(original))

    def test_etch_footprint_is_a_bound_source_or_an_explicit_producer_default(self):
        from mesh_stage_contract import (validate_canonical_dag, validate_command_bindings,
                                         validate_seed_corner_isotropy)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            default_case, etch_case = manifest["Cases"][0], manifest["Cases"][1]
            self.assertEqual(default_case["Source"]["EtchFootprint"], "producer-default")
            self.assertIn("RetainedEtch", etch_case["Source"]["Files"])
            # A case declaring neither, both, or another footprint word fails closed.
            for description, mutate in (
                    ("neither", lambda source: source.pop("EtchFootprint")),
                    ("both", lambda source: source["Files"].__setitem__(
                        "RetainedEtch", etch_case["Source"]["Files"]["RetainedEtch"])),
                    ("other", lambda source: source.__setitem__("EtchFootprint", "device"))):
                candidate = copy.deepcopy(manifest); mutate(candidate["Cases"][0]["Source"])
                with self.assertRaisesRegex(ValueError, "etch footprint", msg=description):
                    validate_manifest(candidate, manifest_path)
            validate_manifest(manifest, manifest_path)
            stage_contract = __import__("mesh_stage_contract")

            def reports_for(case_id):
                reports = {stage: json.loads(
                    (root / f"{case_id}--canonical-{stage}.log.json").read_text())
                           for stage in stage_contract.CANONICAL_STAGE_ORDER}
                reports["proper-rigid-publication"] = json.loads(
                    (root / f"{case_id}--identity-proper-rigid-publication.log.json").read_text())
                case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
                paths = {role: root / case_id / item["Name"]
                         for role, item in case["Source"]["Files"].items()}
                binding = {"Transform": IDENTITY, "InputSHA256": {
                    role: item["SHA256"] for role, item in case["Source"]["Files"].items()}}
                return reports, binding, paths

            # The bound footprint is consumed: recorded in the seed census by hash.
            reports, binding, paths = reports_for("subdivided")
            _validate_source_transformation(reports, binding, paths)
            seed = reports["seed-generation"]
            census = validate_seed_corner_isotropy(
                seed, reports["metric-preparation"]["Artifacts"]["restoration-recipe"]["Path"])
            self.assertEqual(census["EtchBoundarySHA256"],
                             etch_case["Source"]["Files"]["RetainedEtch"]["SHA256"])
            # Substituted footprint: the immutable hash differs from what the seed bound.
            substituted = copy.deepcopy(binding)
            substituted["InputSHA256"]["RetainedEtch"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "retained etch"):
                _validate_source_transformation(reports, substituted, paths)
            # Omitted: the case declares a footprint the seed stage did not bind.
            omitted = copy.deepcopy(reports)
            del omitted["seed-generation"]["Inputs"]["source-retained-etch"]
            with self.assertRaisesRegex(ValueError, "retained etch"):
                _validate_source_transformation(omitted, binding, paths)
            # Tampered census: the recorded footprint hash differs from the bound input.
            census_item = seed["Artifacts"]["seed-corner-census"]
            tampered = json.loads(Path(census_item["Path"]).read_text())
            tampered["EtchBoundarySHA256"] = "1" * 64
            tampered_path = root / "census-etch-tampered.json"
            tampered_path.write_text(json.dumps(tampered))
            tampered_report = copy.deepcopy(seed)
            tampered_report["Artifacts"]["seed-corner-census"] = {
                "Path": str(tampered_path), "SHA256": sha256(tampered_path)}
            with self.assertRaisesRegex(ValueError, "etch footprint"):
                validate_seed_corner_isotropy(
                    tampered_report,
                    reports["metric-preparation"]["Artifacts"]["restoration-recipe"]["Path"])
            # The census must measure exactly the contract's boundary labels (the
            # labels the seed is written with): a missing, extra or duplicated label
            # fails closed even when every recorded area is positive.
            recipe_path = reports["metric-preparation"]["Artifacts"]["restoration-recipe"]["Path"]
            contract_labels = {item["Attribute"] for item in json.loads(
                Path(recipe_path).read_text())["SemanticContract"]["BoundaryLabels"]}
            self.assertEqual({row["Attribute"] for row in census["InterfaceAreas"]},
                             contract_labels)
            for description, mutate in (
                    ("missing", lambda rows: rows.pop()),
                    ("extra", lambda rows: rows.append(
                        {"Attribute": 9999, "Name": "surface_9999", "Triangles": 1,
                         "Area": 1.0})),
                    ("duplicate", lambda rows: rows.append(dict(rows[0])))):
                mutated = json.loads(Path(census_item["Path"]).read_text())
                mutate(mutated["InterfaceAreas"])
                mutated_path = root / f"census-labels-{description}.json"
                mutated_path.write_text(json.dumps(mutated))
                mutated_report = copy.deepcopy(seed)
                mutated_report["Artifacts"]["seed-corner-census"] = {
                    "Path": str(mutated_path), "SHA256": sha256(mutated_path)}
                with self.assertRaisesRegex(ValueError, "labels differ", msg=description):
                    validate_seed_corner_isotropy(mutated_report, recipe_path)
            # The simplified footprint polygons are bound to the shared tolerance and
            # to their own record: another tolerance, a deviation beyond tolerance x
            # local scale, an inconsistent record or a stale summary fails closed.
            def polygon_mutation(**changes):
                def mutate(data):
                    for key, value in changes.items():
                        if key == "FootprintCollinearTolerance":
                            data[key] = value
                        elif key == "FootprintPolygons":
                            data[key] = value
                        elif key == "FootprintSimplification":
                            data[key].update(value)
                        else:
                            data["FootprintPolygons"][0]["Simplification"][key] = value
                return mutate
            for description, mutate, message in (
                    ("tolerance", polygon_mutation(FootprintCollinearTolerance=1e-5),
                     "COPLANAR_TOLERANCE"),
                    ("deviation", polygon_mutation(MaximumDeviation=1e-3,
                                                  MaximumRelativeDeviation=1e-4),
                     "exceeds the collinearity tolerance"),
                    ("count", polygon_mutation(RemovedVertexCount=2), "inconsistent"),
                    ("index", polygon_mutation(RemovedVertexIndices=[9]), "inconsistent"),
                    ("empty", polygon_mutation(FootprintPolygons=[]), "lacks the simplified"),
                    ("summary", polygon_mutation(FootprintSimplification={"RemovedVertices": 7}),
                     "summary differs")):
                mutated = json.loads(Path(census_item["Path"]).read_text())
                mutate(mutated)
                mutated_path = root / f"census-footprint-{description}.json"
                mutated_path.write_text(json.dumps(mutated))
                mutated_report = copy.deepcopy(seed)
                mutated_report["Artifacts"]["seed-corner-census"] = {
                    "Path": str(mutated_path), "SHA256": sha256(mutated_path)}
                with self.assertRaisesRegex(ValueError, message, msg=description):
                    validate_seed_corner_isotropy(mutated_report, recipe_path)
            # The metric recipe records the census's simplified footprint edges with
            # their provenance; a missing, re-sourced or edited record fails closed,
            # and the metric stage must bind exactly the seed's census.
            from mesh_stage_contract import footprint_segments, validate_footprint_segments
            recipe_data = json.loads(Path(recipe_path).read_text())
            self.assertEqual(recipe_data["FootprintSegments"]["Provenance"],
                             etch_case["Source"]["Files"]["RetainedEtch"]["SHA256"])
            self.assertEqual(recipe_data["FootprintSegments"]["Segments"],
                             footprint_segments(census))
            self.assertEqual(len(recipe_data["FootprintSegments"]["Segments"]), 4)
            validate_footprint_segments(recipe_data, census)
            for description, mutate in (
                    ("missing", lambda data: data.pop("FootprintSegments")),
                    ("provenance", lambda data: data["FootprintSegments"].__setitem__(
                        "Provenance", "producer-default")),
                    ("segment", lambda data: data["FootprintSegments"]["Segments"][0].__setitem__(
                        0, 0.5)),
                    ("dropped", lambda data: data["FootprintSegments"]["Segments"].pop()),
                    ("tolerance", lambda data: data["FootprintSegments"].__setitem__(
                        "Tolerance", 1e-5))):
                mutated = copy.deepcopy(recipe_data); mutate(mutated)
                with self.assertRaisesRegex(ValueError, "FootprintSegments differ",
                                            msg=description):
                    validate_footprint_segments(mutated, census)
                mutated_path = root / f"recipe-footprint-{description}.json"
                mutated_path.write_text(json.dumps(mutated))
                with self.assertRaisesRegex(ValueError, "FootprintSegments differ",
                                            msg=description):
                    validate_seed_corner_isotropy(seed, mutated_path)
            metric = reports["metric-preparation"]
            self.assertEqual(metric["Inputs"]["seed-corner-census"]["SHA256"],
                             census_item["SHA256"])
            other_census = copy.deepcopy(reports)
            other_census["metric-preparation"]["Inputs"]["seed-corner-census"] = {
                "Path": str(tampered_path), "SHA256": sha256(tampered_path)}
            paths_by_stage = {stage: root / f"subdivided--canonical-{stage}.log.json"
                              for stage in stage_contract.CANONICAL_STAGE_ORDER}
            rebound = root / "metric-other-census.log.json"
            rebound.write_text(json.dumps(other_census["metric-preparation"]))
            # (rejected by the command binding before the digest-chain link is reached)
            with self.assertRaisesRegex(ValueError, "seed-census"):
                validate_canonical_dag({**paths_by_stage, "metric-preparation": rebound},
                                       root / "subdivided--canonical.msh")
            metric_inputs = {name: item["Path"] for name, item in metric["Inputs"].items()}
            metric_artifacts = {name: item["Path"] for name, item in metric["Artifacts"].items()}
            validate_command_bindings("metric-preparation", metric["Command"], metric_inputs,
                                      metric_artifacts, metric["WorkingDirectory"])
            with self.assertRaises(ValueError):
                validate_command_bindings(
                    "metric-preparation", metric["Command"],
                    dict(metric_inputs, **{"seed-corner-census": str(tampered_path)}),
                    metric_artifacts, metric["WorkingDirectory"])
            # The seed command must pass exactly the bound path, and never an
            # undeclared footprint.
            input_paths = {name: item["Path"] for name, item in seed["Inputs"].items()}
            artifact_paths = {name: item["Path"] for name, item in seed["Artifacts"].items()}
            validate_command_bindings("seed-generation", seed["Command"], input_paths,
                                      artifact_paths, seed["WorkingDirectory"])
            other = dict(input_paths, **{"source-retained-etch": str(root / "base" / "boundary.csv")})
            with self.assertRaises(ValueError):
                validate_command_bindings("seed-generation", seed["Command"], other,
                                          artifact_paths, seed["WorkingDirectory"])
            undeclared = {name: path for name, path in input_paths.items()
                          if name != "source-retained-etch"}
            with self.assertRaisesRegex(ValueError, "without a bound"):
                validate_command_bindings("seed-generation", seed["Command"], undeclared,
                                          artifact_paths, seed["WorkingDirectory"])
            # A producer-default case: the seed binds no footprint and says so.
            reports, binding, paths = reports_for("base")
            _validate_source_transformation(reports, binding, paths)
            default_census = validate_seed_corner_isotropy(
                reports["seed-generation"],
                reports["metric-preparation"]["Artifacts"]["restoration-recipe"]["Path"])
            self.assertEqual(default_census["EtchBoundary"], "producer-default")
            undeclared_bound = copy.deepcopy(reports)
            undeclared_bound["seed-generation"]["Inputs"]["source-retained-etch"] = \
                copy.deepcopy(seed["Inputs"]["source-retained-etch"])
            with self.assertRaisesRegex(ValueError, "does not declare"):
                _validate_source_transformation(undeclared_bound, binding, paths)

    def _stage_report(self, root, name):
        path = root / f"{name}.log.json"
        return path, json.loads(path.read_text())

    def _substitute(self, command, option, replacement):
        command = list(command); command[command.index(option) + 1] = str(replacement)
        return command

    def test_trace_basis_is_bound_to_seed_and_metric_together_or_not_at_all(self):
        from mesh_stage_contract import (bound_trace_basis, validate_canonical_dag,
                                         validate_trace_basis_sizing)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            plain_case, basis_case = manifest["Cases"][0], manifest["Cases"][2]
            self.assertTrue(set(TRACE_BASIS_FILES) <= set(basis_case["Source"]["Files"]))
            self.assertFalse(set(TRACE_BASIS_FILES) & set(plain_case["Source"]["Files"]))
            # A case freezing some but not all trace basis roles fails preflight.
            for role in TRACE_BASIS_FILES:
                candidate = copy.deepcopy(manifest)
                del candidate["Cases"][2]["Source"]["Files"][role]
                with self.assertRaisesRegex(ValueError, "trace basis roles", msg=role):
                    validate_manifest(candidate, manifest_path)
            validate_manifest(manifest, manifest_path)

            def reports_for(case_id):
                reports = {stage: json.loads(
                    (root / f"{case_id}--canonical-{stage}.log.json").read_text())
                           for stage in CANONICAL_STAGE_ORDER}
                reports["proper-rigid-publication"] = json.loads(
                    (root / f"{case_id}--identity-proper-rigid-publication.log.json").read_text())
                case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
                paths = {role: root / case_id / item["Name"]
                         for role, item in case["Source"]["Files"].items()}
                binding = {"Transform": IDENTITY, "InputSHA256": {
                    role: item["SHA256"] for role, item in case["Source"]["Files"].items()}}
                return reports, binding, paths

            reports, binding, paths = reports_for(basis_case["Id"])
            _validate_source_transformation(reports, binding, paths)
            seed, metric = reports["seed-generation"], reports["metric-preparation"]
            recipe_path = Path(metric["Artifacts"]["restoration-recipe"]["Path"])
            census_path = Path(seed["Artifacts"]["seed-corner-census"]["Path"])
            recipe, census = json.loads(recipe_path.read_text()), json.loads(census_path.read_text())
            expected = {role: basis_case["Source"]["Files"][role]["SHA256"] for role in TRACE_BASIS_FILES}
            self.assertEqual(bound_trace_basis(seed), expected)
            self.assertEqual(bound_trace_basis(metric), expected)
            record = validate_trace_basis_sizing(seed, metric, recipe, census)
            self.assertEqual(record["Ratio"], TRACE_BASIS_SIZE_RATIO)
            self.assertTrue(record["RatioIsDimensionless"])
            self.assertEqual(record["InputSHA256"], expected)
            self.assertEqual(census["TraceBasisSizing"]["InputSHA256"], expected)
            stage_paths = {stage: root / f"{basis_case['Id']}--canonical-{stage}.log.json"
                           for stage in CANONICAL_STAGE_ORDER}
            validate_canonical_dag(stage_paths, root / f"{basis_case['Id']}--canonical.msh")

            def rewrite(report, changes, name):
                path = root / f"trace-basis-{name}.json"; path.write_text(json.dumps(changes))
                return path

            def rejected(seed_report=seed, metric_report=metric, recipe_data=recipe,
                         census_data=census, message="", description=""):
                with self.assertRaisesRegex(ValueError, message, msg=description):
                    validate_trace_basis_sizing(seed_report, metric_report, recipe_data, census_data)

            def without_option(report, option, count=2):
                report = copy.deepcopy(report); index = report["Command"].index(option)
                del report["Command"][index:index + count]; return report

            def with_option_value(report, option, value):
                report = copy.deepcopy(report)
                report["Command"][report["Command"].index(option) + 1] = value; return report

            # Seed without the basis while the metric binds it (and vice versa).
            seed_unbound = copy.deepcopy(seed)
            for name, _ in TRACE_BASIS_BINDINGS.values():
                del seed_unbound["Inputs"][name]
            rejected(seed_report=seed_unbound, message="different trace bases", description="seed unbound")
            metric_unbound = copy.deepcopy(metric)
            for name, _ in TRACE_BASIS_BINDINGS.values():
                del metric_unbound["Inputs"][name]
            rejected(metric_report=metric_unbound, message="different trace bases", description="metric unbound")
            # Partial binding, substituted digest, different ratios, missing ratio option.
            partial = copy.deepcopy(seed); del partial["Inputs"]["source-trace-vertices"]
            rejected(seed_report=partial, message="incomplete trace basis")
            substituted = copy.deepcopy(metric)
            substituted["Inputs"]["source-basis-contract"]["SHA256"] = "0" * 64
            rejected(metric_report=substituted, message="different trace bases")
            rejected(metric_report=with_option_value(metric, "--trace-basis-size-ratio", "2.0"),
                     message="ratio")
            rejected(seed_report=without_option(seed, "--trace-basis-size-ratio"), message="ratio")
            rejected(seed_report=with_option_value(seed, "--trace-basis-size-ratio", "-1"),
                     message="positive finite")
            # Records: ratio, digests, triangles and statistics must match the binding.
            for description, changes in (
                    ("ratio", {"Ratio": 2 * TRACE_BASIS_SIZE_RATIO}),
                    ("dimensionless", {"RatioIsDimensionless": False}),
                    ("digest", {"InputSHA256": {**expected, "TraceVertices": "1" * 64}}),
                    ("triangles", {"MeshFrameTriangles": [[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]]}),
                    ("count", {"Triangles": 2}),
                    ("statistics", {"CutSurfaceSize": None}),
                    ("edges", {"BasisEdgesBelowFarSize": -1})):
                rejected(recipe_data={**recipe, "TraceBasisSizing": {**recipe["TraceBasisSizing"], **changes}},
                         description=f"recipe {description}")
            rejected(recipe_data={**recipe, "TraceBasisSizing": None}, description="recipe omitted")
            # The audit's source-driven band lines: exactly the basis edges, bound by digest.
            edges = recipe["TraceBasisEdges"]
            self.assertEqual(edges["Count"], len(edges["Segments"]))
            self.assertEqual(edges["InputSHA256"], expected)
            for description, changes in (
                    ("digest", {"InputSHA256": {**expected, "BasisContract": "4" * 64}}),
                    ("dropped", {"Segments": edges["Segments"][1:], "Count": edges["Count"] - 1}),
                    ("count", {"Count": edges["Count"] + 1}),
                    # One-to-one: a duplicated segment must not cover a missing one.
                    ("duplicated", {"Segments": [edges["Segments"][0]] + edges["Segments"][:-1]}),
                    ("moved", {"Segments": [[v + 1.0 for v in segment] for segment in edges["Segments"]]})):
                rejected(recipe_data={**recipe, "TraceBasisEdges": {**edges, **changes}},
                         message="trace basis edge", description=f"edges {description}")
            rejected(recipe_data={k: v for k, v in recipe.items() if k != "TraceBasisEdges"},
                     message="trace basis edge", description="edges omitted")
            for description, changes in (
                    ("ratio", {"Ratio": 2 * TRACE_BASIS_SIZE_RATIO}),
                    ("digest", {"InputSHA256": {**expected, "ProcessLibrary": "2" * 64}}),
                    ("triangles", {"MeshFrameTriangles": [[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]]})):
                rejected(census_data={**census, "TraceBasisSizing": {**census["TraceBasisSizing"], **changes}},
                         description=f"census {description}")
            rejected(census_data={**census, "TraceBasisSizing": None}, description="census omitted")
            # The DAG validator applies the same rule through the bound reports.
            tampered_recipe = rewrite(metric, {**recipe, "TraceBasisSizing": None}, "recipe-omitted")
            tampered_metric = copy.deepcopy(metric)
            tampered_metric["Artifacts"]["restoration-recipe"] = {
                "Path": str(tampered_recipe), "SHA256": sha256(tampered_recipe)}
            tampered_metric_path = rewrite(metric, tampered_metric, "metric-report")
            with self.assertRaises(ValueError):
                validate_canonical_dag({**stage_paths, "metric-preparation": tampered_metric_path},
                                       root / f"{basis_case['Id']}--canonical.msh")
            # Manifest binding: the seed/metric must consume exactly the declared basis.
            undeclared = copy.deepcopy(binding)
            for role in TRACE_BASIS_FILES:
                del undeclared["InputSHA256"][role]
            with self.assertRaisesRegex(ValueError, "does not declare"):
                _validate_source_transformation(reports, undeclared, paths)
            swapped = copy.deepcopy(binding); swapped["InputSHA256"]["TraceTriangles"] = "3" * 64
            with self.assertRaisesRegex(ValueError, "immutable trace basis"):
                _validate_source_transformation(reports, swapped, paths)
            # A case without a basis: nothing bound, nothing recorded, no ratio option.
            plain_reports, plain_binding, plain_paths = reports_for(plain_case["Id"])
            plain_seed, plain_metric = plain_reports["seed-generation"], plain_reports["metric-preparation"]
            plain_recipe = json.loads(Path(plain_metric["Artifacts"]["restoration-recipe"]["Path"]).read_text())
            plain_census = json.loads(Path(plain_seed["Artifacts"]["seed-corner-census"]["Path"]).read_text())
            self.assertIsNone(bound_trace_basis(plain_seed))
            self.assertIsNone(validate_trace_basis_sizing(plain_seed, plain_metric, plain_recipe, plain_census))
            self.assertIsNone(plain_census["TraceBasisSizing"])
            self.assertNotIn("TraceBasisSizing", plain_recipe)
            self.assertNotIn("TraceBasisEdges", plain_recipe)
            with_ratio = copy.deepcopy(plain_seed); with_ratio["Command"] += ["--trace-basis-size-ratio", "1.0"]
            rejected(seed_report=with_ratio, metric_report=plain_metric, recipe_data=plain_recipe,
                     census_data=plain_census, message="without a bound trace basis")
            rejected(seed_report=plain_seed, metric_report=plain_metric,
                     recipe_data={**plain_recipe, "TraceBasisSizing": recipe["TraceBasisSizing"]},
                     census_data=plain_census, message="without a bound trace basis")
            bound_plain = copy.deepcopy(plain_reports)
            bound_plain["seed-generation"]["Inputs"].update(
                {name: {"Path": str(paths[role]), "SHA256": expected[role]}
                 for role, (name, _) in TRACE_BASIS_BINDINGS.items()})
            with self.assertRaisesRegex(ValueError, "does not declare"):
                _validate_source_transformation(bound_plain, plain_binding, plain_paths)
            # Stage command bindings: the options must name the bound inputs exactly.
            from mesh_stage_contract import validate_command_bindings
            inputs = {name: item["Path"] for name, item in seed["Inputs"].items()}
            artifacts = {name: item["Path"] for name, item in seed["Artifacts"].items()}
            validate_command_bindings("seed-generation", seed["Command"], inputs, artifacts,
                                      seed["WorkingDirectory"])
            with self.assertRaisesRegex(ValueError, "differs from its bound input"):
                validate_command_bindings(
                    "seed-generation", with_option_value(seed, "--trace-vertices", "/other.csv")["Command"],
                    inputs, artifacts, seed["WorkingDirectory"])
            with self.assertRaisesRegex(ValueError, "without a bound"):
                validate_command_bindings("seed-generation", seed["Command"],
                                          {k: v for k, v in inputs.items() if k != "source-trace-vertices"},
                                          artifacts, seed["WorkingDirectory"])

    def test_stage_commands_must_consume_exactly_the_bound_source_and_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            alternate = root / "subdivided"
            substitutions = {
                "seed-generation": ("base--canonical-seed-generation",
                                    ["--mask", "--boundary", "--signature"]),
                "canonical-gmsh-publication": ("base--canonical-canonical-gmsh-publication",
                                               ["--process", "--signature", "--boundary"]),
                "proper-rigid-publication": ("base--identity-proper-rigid-publication",
                                             ["--semantic-input", "--signature", "--boundary",
                                              "--mask", "--process"]),
            }
            for stage, (name, options) in substitutions.items():
                _, report = self._stage_report(root, name)
                validate_stage_report(report, stage)
                for option in options:
                    with self.subTest(stage=stage, option=option):
                        original = Path(report["Command"][report["Command"].index(option) + 1])
                        changed = copy.deepcopy(report)
                        changed["Command"] = self._substitute(
                            report["Command"], option, alternate / original.name)
                        # Named options must equal the binding; the seeder's signature
                        # is consumed positionally and must occur in argv.
                        with self.assertRaisesRegex(
                                ValueError, "differs from its bound|did not consume its bound"):
                            validate_stage_report(changed, stage)
                        dropped = copy.deepcopy(report)
                        position = dropped["Command"].index(option)
                        del dropped["Command"][position:position + 2]
                        with self.assertRaisesRegex(
                                ValueError, "exactly one|did not consume its bound"):
                            validate_stage_report(dropped, stage)
            _, placement = self._stage_report(root, "base--identity-proper-rigid-publication")
            other = root / "base--rotate-z-0.63-ownership.csv"
            for option in ("--ownership", "--ownership-quadrature", "--canonical-build-record",
                           "--transformed-semantic", "--transformed-supports"):
                with self.subTest(option=option):
                    changed = copy.deepcopy(placement)
                    changed["Command"] = self._substitute(placement["Command"], option, other)
                    with self.assertRaisesRegex(ValueError, "differs from its bound"):
                        validate_stage_report(changed, "proper-rigid-publication")
            with self.assertRaisesRegex(ValueError, "did not consume its bound"):
                changed = copy.deepcopy(placement)
                changed["Command"] = [value if value != str(root / "base--canonical.msh") else
                                      str(root / "subdivided--canonical.msh")
                                      for value in placement["Command"]]
                validate_stage_report(changed, "proper-rigid-publication")
            self.assertNotIn("--source", placement["Command"])

    def test_ownership_audit_must_consume_bound_source_and_write_bound_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            _, report = self._stage_report(root, "base--identity-proper-rigid-publication")
            receipt_path = Path(report["Artifacts"]["transform-receipt"]["Path"])
            original = json.loads(receipt_path.read_text())
            def rebind(receipt):
                receipt_path.write_text(json.dumps(receipt))
                changed = copy.deepcopy(report)
                changed["Artifacts"]["transform-receipt"]["SHA256"] = sha256(receipt_path)
                return changed
            try:
                for option in ("--process", "--signature", "--boundary"):
                    with self.subTest(option=option):
                        receipt = copy.deepcopy(original)
                        receipt["OwnershipCommand"] = self._substitute(
                            original["OwnershipCommand"], option,
                            root / "subdivided" / Path(original["OwnershipCommand"][
                                original["OwnershipCommand"].index(option) + 1]).name)
                        with self.assertRaisesRegex(ValueError, "Ownership audit command"):
                            validate_stage_report(rebind(receipt), "proper-rigid-publication")
                receipt = copy.deepcopy(original)
                receipt["OwnershipCommand"] = [
                    str(root / "base--rotate-z-0.63-ownership.csv")
                    if value == report["Artifacts"]["ownership-partition"]["Path"] else value
                    for value in original["OwnershipCommand"]]
                with self.assertRaisesRegex(ValueError, "Ownership audit command"):
                    validate_stage_report(rebind(receipt), "proper-rigid-publication")
                receipt = copy.deepcopy(original)
                receipt["OwnershipCommand"][1] = str(MESHER)
                with self.assertRaisesRegex(ValueError, "Ownership audit command"):
                    validate_stage_report(rebind(receipt), "proper-rigid-publication")
                for key in ("OwnershipSHA256", "OwnershipAuditorSHA256"):
                    with self.subTest(key=key):
                        receipt = copy.deepcopy(original); receipt[key] = "0" * 64
                        with self.assertRaisesRegex(ValueError, "differ from bindings"):
                            validate_stage_report(rebind(receipt), "proper-rigid-publication")
                receipt = copy.deepcopy(original)
                receipt["SourceInputSHA256"]["source-process"] = "0" * 64
                with self.assertRaisesRegex(ValueError, "differ from bindings"):
                    validate_stage_report(rebind(receipt), "proper-rigid-publication")
            finally:
                receipt_path.write_text(json.dumps(original) + "\n")

    def test_changed_mmg_library_hash_is_rejected_everywhere(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            _, report = self._stage_report(root, "base--canonical-native-adaptation-mmg")
            validate_stage_report(report, "native-adaptation-mmg")
            self.assertEqual(report["Tools"]["mmg-library"]["SHA256"], sha256(MMG_LIBRARY))
            receipt_path = Path(report["Artifacts"]["adaptation-receipt"]["Path"])
            original = json.loads(receipt_path.read_text())
            self.assertEqual(original["MMGLibrarySHA256"], sha256(MMG_LIBRARY))
            try:
                receipt = copy.deepcopy(original); receipt["MMGLibrarySHA256"] = "0" * 64
                receipt_path.write_text(json.dumps(receipt))
                changed = copy.deepcopy(report)
                changed["Artifacts"]["adaptation-receipt"]["SHA256"] = sha256(receipt_path)
                with self.assertRaisesRegex(ValueError, "adapter/library differ"):
                    validate_stage_report(changed, "native-adaptation-mmg")
            finally:
                receipt_path.write_text(json.dumps(original, indent=2) + "\n")
            replaced = copy.deepcopy(manifest)
            replaced["StageToolSHA256"]["native-adaptation-mmg"]["mmg-library"] = "0" * 64
            replaced_path = root / "replaced-library-suite.json"
            replaced_path.write_text(json.dumps(replaced))
            stem = "base--identity"
            records = {kind: root / f"{stem}-{kind}.json" for kind in KINDS}
            with self.assertRaisesRegex(ValueError, "mmg-library"):
                normalize(replaced_path, "base", "identity", root / f"{stem}.msh", records,
                          root / "replaced-library-normalized.json")
            changed = copy.deepcopy(report)
            changed["Tools"]["mmg-library"]["Path"] = str(MESHER)
            with self.assertRaisesRegex(ValueError, "mmg-library"):
                validate_stage_report(changed, "native-adaptation-mmg")
            self.assertEqual(
                report["Command"][report["Command"].index("--mmg-library") + 1],
                str(MMG_LIBRARY))

    def test_quadrature_partition_rejects_changed_missing_and_duplicate_owner(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summary = root / "ownership.csv"
            header = ("attribute,elements,ambiguous_fraction,unresolved_elements,"
                      "unresolved_fraction,quadrature_rule,quadrature_order,"
                      "quadrature_points,quadrature_whole_measure,quadrature_owned_measure,"
                      "quadrature_relative_closure,quadrature_closure_tolerance,"
                      "quadrature_unmatched,quadrature_overlaps,quadrature_positive_weights\n")
            suffix = ",1,0,0,0,Gauss4,4,8,3,3,0,1e-12,0,0,1\n"
            summary.write_text(header + "2" + suffix + "3" + suffix)
            quadrature = root / "quadrature.csv"
            contract = {"CutSurfaceRoles": ["matching"], "BoundaryLabels": [
                {"Attribute": 1, "Role": "matching"},
                {"Attribute": 2, "Role": "slot-0"},
                {"Attribute": 3, "Role": "slot-1"}]}
            quadrature.write_text("attribute,measure\n2,2\n3,1\n")
            report = _ownership_report(summary, quadrature, contract)
            self.assertEqual(report["ResponseOwnership"]["OwnerAttributes"], [2, 3])
            for name, value in {
                    "changed": "attribute,measure\n2,2\n4,1\n",
                    "missing": "attribute,measure\n2,3\n",
                    "duplicate": "attribute,measure\n2,2\n2,1\n"}.items():
                with self.subTest(name=name):
                    quadrature.write_text(value)
                    with self.assertRaises(ValueError):
                        _ownership_report(summary, quadrature, contract)

    def test_non_native_primary_tool_must_be_first_executed_script(self):
        tools = {"runtime": sys.executable, "metric-preparer": STAGER}
        validate_tool_invocation(
            "metric-preparation", [sys.executable, "-B", str(STAGER)], tools)
        with self.assertRaisesRegex(ValueError, "executed script position"):
            validate_tool_invocation(
                "metric-preparation",
                [sys.executable, "-B", str(MESHER), str(STAGER)], tools)

    def test_relative_tool_argv_must_resolve_to_the_bound_tool_in_its_working_directory(self):
        tools = {"runtime": sys.executable, "metric-preparer": STAGER}
        relative = str(STAGER.relative_to(HERE.parent.parent.parent))
        validate_tool_invocation(
            "metric-preparation", [sys.executable, relative], tools,
            working_directory=HERE.parent.parent.parent)
        with tempfile.TemporaryDirectory() as temporary:
            # Same repository-relative spelling, but launched from a directory
            # holding a different file at that path: the executed script is
            # not the bound tool even though the suffix matches.
            copied = Path(temporary) / relative
            copied.parent.mkdir(parents=True)
            copied.write_text("print('impostor')\n")
            with self.assertRaisesRegex(ValueError, "executed script position"):
                validate_tool_invocation(
                    "metric-preparation", [sys.executable, relative], tools,
                    working_directory=temporary)
            with self.assertRaisesRegex(ValueError, "runtime is not argv"):
                validate_tool_invocation(
                    "metric-preparation",
                    [Path(sys.executable).name, str(STAGER)], tools,
                    working_directory=temporary)

    def test_uninvoked_adapter_mmg_stage_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed = root / "seed"; seed.write_text("seed")
            metric = root / "metric"; metric.write_text("metric")
            pins = root / "pins"; pins.write_text("pins")
            fixed = root / "fixed"; fixed.write_text("fixed")
            required = root / "required"; required.write_text("1\n")
            recipe = root / "recipe"; recipe.write_text("{}")
            result = subprocess.run(
                [sys.executable, str(BOUNDED), "--seconds", "1", "--memory-gib", "1",
                 "--log", str(root / "bad.log"), "--stage", "native-adaptation-mmg",
                 "--input", f"mmg-seed={seed}", "--input", f"metric={metric}",
                 "--input", f"pins={pins}", "--input", f"fixed-triangles={fixed}",
                 "--input", f"required-tetrahedra={required}",
                 "--input", f"restoration-recipe={recipe}",
                 "--artifact", f"adapted-mesh={root / 'adapted.msh'}",
                 "--artifact", f"adaptation-receipt={root / 'receipt.json'}",
                 "--tool", f"runtime={sys.executable}",
                 "--tool", f"adaptation-wrapper={WRAPPER}",
                 "--tool", f"adapter-mmg={SCRIPT_ADAPTER}",
                 "--tool", f"mmg-library={SCRIPT_ADAPTER}", "--",
                 sys.executable, str(WRAPPER), str(seed), str(metric), str(pins),
                 str(recipe), str(root / "adapted.msh"), str(root / "receipt.json"),
                 "--mmg-library", str(SCRIPT_ADAPTER), "--fixed-triangles", str(fixed),
                 "--required-tetrahedra", str(required)],
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--adapter", result.stderr)
            # The required-tetrahedra list is a bound adaptation input: a stage that
            # neither binds nor passes it fails closed before launching.
            for omit_input, omit_option in ((True, False), (False, True)):
                command = [sys.executable, str(BOUNDED), "--seconds", "1", "--memory-gib", "1",
                           "--log", str(root / f"bad-{omit_input}-{omit_option}.log"),
                           "--stage", "native-adaptation-mmg",
                           "--input", f"mmg-seed={seed}", "--input", f"metric={metric}",
                           "--input", f"pins={pins}", "--input", f"fixed-triangles={fixed}",
                           *([] if omit_input else ["--input", f"required-tetrahedra={required}"]),
                           "--input", f"restoration-recipe={recipe}",
                           "--artifact", f"adapted-mesh={root / 'adapted.msh'}",
                           "--artifact", f"adaptation-receipt={root / 'receipt.json'}",
                           "--tool", f"runtime={sys.executable}",
                           "--tool", f"adaptation-wrapper={WRAPPER}",
                           "--tool", f"adapter-mmg={SCRIPT_ADAPTER}",
                           "--tool", f"mmg-library={SCRIPT_ADAPTER}", "--",
                           sys.executable, str(WRAPPER), str(seed), str(metric), str(pins),
                           str(recipe), str(root / "adapted.msh"), str(root / "receipt.json"),
                           "--adapter", str(SCRIPT_ADAPTER), "--mmg-library", str(SCRIPT_ADAPTER),
                           "--fixed-triangles", str(fixed),
                           *([] if omit_option else ["--required-tetrahedra", str(required)])]
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue("required-tetrahedra" in result.stderr or
                                "frozen stage contract" in result.stderr)
                self.assertFalse((root / "adapted.msh").exists())

    def _seed_stage_launch(self, root, output, *, semantic_option=True):
        sources = {name: root / f"{name}.csv" for name in ("signature", "boundary", "mask")}
        for path in sources.values(): path.write_text(path.name)
        semantic = root / "semantic.json"; semantic.write_text("{}")
        census = root / "census.json"
        return [sys.executable, str(BOUNDED), "--seconds", "1", "--memory-gib", "1",
                "--log", str(root / "seed.log"), "--stage", "seed-generation",
                *[value for name, path in sources.items()
                  for value in ("--input", f"source-{name}={path}")],
                "--input", f"canonical-semantic-contract={semantic}",
                "--artifact", f"seed-mesh={output}", "--artifact", f"seed-corner-census={census}",
                "--tool", f"runtime={sys.executable}", "--tool", f"mesher={MESHER}",
                "--", sys.executable, str(MESHER), str(output), str(root / "missing"),
                *[value for name, path in sources.items()
                  for value in (f"--{name}", str(path))],
                *(["--semantic-contract", str(semantic)] if semantic_option else []),
                "--corner-isotropy-radius", "0.5", "--lc-fine", "0.1",
                "--corner-census", str(census)]

    def test_preexisting_stage_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "seed.msh"
            output.write_text("preexisting")
            result = subprocess.run(self._seed_stage_launch(root, output),
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stage output must be absent", result.stderr)

    def test_seed_stage_without_the_bound_semantic_contract_option_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "seed.msh"
            result = subprocess.run(
                self._seed_stage_launch(root, output, semantic_option=False),
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exactly one --semantic-contract", result.stderr)
            self.assertFalse(output.exists())

    def test_seed_corner_isotropy_must_equal_the_recipe_prescription(self):
        from mesh_stage_contract import validate_canonical_dag, validate_seed_corner_isotropy
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            reports = {stage: root / f"base--canonical-{stage}.log.json"
                       for stage in CANONICAL_STAGE_ORDER}
            recipe = root / "base--canonical-recipe.json"
            seed_report = json.loads(reports["seed-generation"].read_text())
            validate_canonical_dag(reports, root / "base--canonical.msh")
            census = validate_seed_corner_isotropy(seed_report, recipe)
            self.assertEqual(census["CornerIsotropyRadius"], CORNER_ISOTROPY_RADIUS)
            self.assertEqual(census["IsotropicSize"], NORMAL_SIZE)
            self.assertEqual(census["SemanticCorners"],
                             json.loads(recipe.read_text())["TruePhysicalCorners"])

            def tampered(mutate, description):
                report = copy.deepcopy(seed_report); mutate(report)
                path = root / f"tampered-{description}.json"
                path.write_text(json.dumps(report))
                with self.assertRaises(ValueError):
                    validate_seed_corner_isotropy(report, recipe)
                with self.assertRaises(ValueError):
                    validate_canonical_dag({**reports, "seed-generation": path},
                                           root / "base--canonical.msh")

            def replace_option(option, value):
                def mutate(report):
                    index = report["Command"].index(option) + 1
                    report["Command"][index] = value
                return mutate

            def remove_option(option):
                def mutate(report):
                    index = report["Command"].index(option)
                    del report["Command"][index:index + 2]
                return mutate

            def rewrite_census(**changes):
                def mutate(report):
                    item = report["Artifacts"]["seed-corner-census"]
                    data = json.loads(Path(item["Path"]).read_text()); data.update(changes)
                    path = root / f"census-{'-'.join(changes)}.json"
                    path.write_text(json.dumps(data))
                    item["Path"] = str(path); item["SHA256"] = sha256(path)
                return mutate

            tampered(replace_option("--corner-isotropy-radius", str(2 * CORNER_ISOTROPY_RADIUS)),
                     "radius-mismatch")
            tampered(replace_option("--lc-fine", str(2 * NORMAL_SIZE)), "size-mismatch")
            tampered(remove_option("--corner-isotropy-radius"), "radius-omitted")
            tampered(remove_option("--lc-fine"), "size-omitted")
            tampered(rewrite_census(CornerIsotropyRadius=2 * CORNER_ISOTROPY_RADIUS),
                     "census-radius")
            tampered(rewrite_census(IsotropicSize=2 * NORMAL_SIZE), "census-size")
            tampered(rewrite_census(SemanticCorners=[[9., 9., 9.]]), "census-corners")
            tampered(rewrite_census(LongitudinalFaces=None), "census-faces-omitted")
            tampered(rewrite_census(LongitudinalFaces=[{"Surface": 1}]), "census-faces-incomplete")
            # The recipe's junction segments must be the census's CAD junction curves.
            census_data = json.loads(
                Path(seed_report["Artifacts"]["seed-corner-census"]["Path"]).read_text())
            curves = census_data["JunctionCurves"]
            tampered(rewrite_census(JunctionCurves=None), "census-junction-omitted")
            tampered(rewrite_census(JunctionCurves={**curves, "TotalLength": 2 * curves["TotalLength"]}),
                     "census-junction-length")
            tampered(rewrite_census(JunctionCurves={**curves, "Count": 0}), "census-junction-count")

            # The metric stage must have frozen the seed corner balls it received.
            from mesh_stage_contract import validate_protected_corner_balls
            recipe_data = json.loads(recipe.read_text())
            balls = validate_protected_corner_balls(recipe_data)
            self.assertEqual(balls["Radius"], recipe_data["CornerIsotropyRadius"])
            self.assertEqual(len(balls["PerCorner"]), len(recipe_data["TruePhysicalCorners"]))

            def recipe_without(**changes):
                data = copy.deepcopy(recipe_data)
                for key, value in changes.items():
                    if value is None: data.pop(key, None)
                    else: data[key] = value
                return data

            from mesh_stage_contract import validate_junction_segments
            junction = recipe_data["JunctionSegments"]
            self.assertEqual(validate_junction_segments(recipe_data, census_data), junction)
            self.assertEqual(junction["CutSurfaceAttributes"], [1])
            self.assertEqual(junction["MaterialInterfaceAttributes"], [3])
            for description, data in (
                    ("junction-omitted", recipe_without(JunctionSegments=None)),
                    ("junction-dropped", recipe_without(JunctionSegments={
                        **junction, "Segments": junction["Segments"][1:],
                        "Count": junction["Count"] - 1})),
                    ("junction-count", recipe_without(JunctionSegments={
                        **junction, "Count": junction["Count"] + 1})),
                    ("junction-length", recipe_without(JunctionSegments={
                        **junction, "TotalLength": 2 * junction["TotalLength"]})),
                    ("junction-degenerate", recipe_without(JunctionSegments={
                        **junction, "Segments": junction["Segments"][:-1] + [[0., 0., 0., 0., 0., 0.]]})),
                    ("junction-cut-labels", recipe_without(JunctionSegments={
                        **junction, "CutSurfaceAttributes": [3]})),
                    ("junction-interface-labels", recipe_without(JunctionSegments={
                        **junction, "MaterialInterfaceAttributes": [1, 3]}))):
                with self.assertRaises(ValueError, msg=description):
                    validate_junction_segments(data, census_data)
                path = root / f"recipe-{description}.json"
                path.write_text(json.dumps(data))
                with self.assertRaises(ValueError, msg=description):
                    validate_seed_corner_isotropy(seed_report, path)

            for description, data in (
                    ("omitted", recipe_without(ProtectedCornerBalls=None)),
                    ("radius", recipe_without(ProtectedCornerBalls={
                        **recipe_data["ProtectedCornerBalls"],
                        "Radius": 2 * recipe_data["CornerIsotropyRadius"]})),
                    ("corners", recipe_without(ProtectedCornerBalls={
                        **recipe_data["ProtectedCornerBalls"],
                        "PerCorner": [{"Point": [9., 9., 9.], "FrozenTriangles": 1}]})),
                    ("uncovered", recipe_without(ProtectedCornerBalls={
                        **recipe_data["ProtectedCornerBalls"], "FrozenTriangles": 10 ** 6})),
                    ("count", recipe_without(ProtectedCornerBalls={
                        **recipe_data["ProtectedCornerBalls"], "FrozenTriangles": -1}))):
                with self.assertRaises(ValueError, msg=description):
                    validate_protected_corner_balls(data)
                path = root / f"recipe-{description}.json"
                path.write_text(json.dumps(data))
                with self.assertRaises(ValueError, msg=description):
                    validate_seed_corner_isotropy(seed_report, path)

    @staticmethod
    def _star_mesh(points, faces, materials):
        """Closed labeled surface -> tetrahedral mesh (one star per material body).

        `faces` are (triangle, label, material-or-materials): every triangle is
        joined to the centroid of each listed body (two bodies for an interface),
        so the surface is exactly the mesh boundary/interface set."""
        points = [list(map(float, point)) for point in points]
        bodies = [(t, m) for t, _, ms in faces
                  for m in (ms if isinstance(ms, tuple) else (ms,))]
        tetrahedra, tetrahedron_materials, centroids = [], [], {}
        for triangle, material in bodies:
            if material not in centroids:
                body = np.unique([t for t, m in bodies if m == material])
                centroids[material] = len(points)
                points.append(np.mean(np.asarray(points)[body], axis=0).tolist())
            tetrahedra.append([centroids[material], *triangle])
            tetrahedron_materials.append(material)
        xyz = np.asarray(points)[np.asarray(tetrahedra)]
        signed = np.einsum("ij,ij->i", np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]),
                           xyz[:, 3] - xyz[:, 0])
        tetrahedra = np.asarray(tetrahedra)
        tetrahedra[signed < 0] = tetrahedra[signed < 0][:, [0, 2, 1, 3]]
        triangles = np.asarray([triangle for triangle, _, _ in faces])
        labels = np.asarray([label for _, label, _ in faces])
        names = {1: "substrate", 2: "vacuum"}
        return meshio.Mesh(np.asarray(points), [("triangle", triangles), ("tetra", tetrahedra)],
                           cell_data={"gmsh:physical": [labels, np.asarray(tetrahedron_materials)]},
                           field_data={names[m]: np.array([m, 3]) for m in centroids})

    @staticmethod
    def _slot_union_contract(roles, materials=(1,)):
        labels = [{"Attribute": 1, "Role": "matching-surface", "AdjacentMaterials": [1],
                   "Protected": True}]
        labels += [{"Attribute": attribute, "Role": role, "AdjacentMaterials": list(adjacent),
                    "Protected": True} for attribute, (role, adjacent) in roles.items()]
        return {"VolumeMaterials": [{"Attribute": 1, "Material": "substrate"}] + (
                    [{"Attribute": 2, "Material": "vacuum"}] if 2 in materials else []),
                "BoundaryLabels": labels, "MetricSurfaceRoles": [], "CutSurfaceRoles": []}

    def _slot_split_cubes(self, bottom_labels, split_candidate=True, transform=None):
        """Unit cube whose bottom face carries two labels (reference: two triangles
        along the diagonal; candidate: four triangles around the face center, one of
        the first label) and five matching-surface faces."""
        corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
        sides = [([0, 1, 5], 1), ([0, 5, 4], 1), ([1, 2, 6], 1), ([1, 6, 5], 1),
                 ([2, 3, 7], 1), ([2, 7, 6], 1), ([3, 0, 4], 1), ([3, 4, 7], 1),
                 ([4, 5, 6], 1), ([4, 6, 7], 1)]
        first, second = bottom_labels
        reference_faces = [([0, 2, 1], first), ([0, 3, 2], second)] + sides
        center = np.array([[.5, .5, 0.]])
        candidate_faces = [([0, 8, 1], first), ([1, 8, 2], second), ([2, 8, 3], second),
                           ([3, 8, 0], second)] + sides if split_candidate else reference_faces
        def mesh(points, faces):
            if transform is not None:
                matrix = np.asarray(transform, dtype=float).reshape(4, 4)
                points = points @ matrix[:3, :3].T + matrix[:3, 3]
            return self._star_mesh(points, [(t, label, 1) for t, label in faces], (1,))
        return (mesh(corners, reference_faces),
                mesh(np.vstack((corners, center)), candidate_faces))

    def test_coplanar_slot_seam_labels_are_audited_as_one_union_patch(self):
        contract = self._slot_union_contract({
            3100: ("etched-substrate-vacuum-slot-0", [1]),
            3101: ("etched-substrate-vacuum-slot-1", [1])})
        reference, candidate = self._slot_split_cubes((3100, 3101))
        result = _protected_surface_report(reference, candidate, contract)
        self.assertTrue(result["PlaneSupportsMatch"]); self.assertTrue(result["TopologyMatches"])
        # Six audited patches: five matching faces and ONE union of the slot labels.
        self.assertEqual(result["PatchCount"], 6); self.assertEqual(result["LabelPatchCount"], 7)
        self.assertEqual(result["CoplanarSlotUnions"], 1)
        self.assertLessEqual(result["MaximumRelativeMeasureError"], 1e-12)
        self.assertLessEqual(result["MaximumSupportVertexDistance"], 1e-12)
        key = next(key for key in result["ByPlaneSupport"] if "+" in key)
        self.assertTrue(key.startswith(
            "etched-substrate-vacuum-slot-0+etched-substrate-vacuum-slot-1 "))
        union = result["ByPlaneSupport"][key]["CoplanarSlotUnion"]
        self.assertEqual(union["Labels"], [3100, 3101])
        self.assertEqual(union["ReferenceLabelAreas"], {"3100": .5, "3101": .5})
        self.assertAlmostEqual(union["CandidateLabelAreas"]["3100"], .25)
        self.assertAlmostEqual(union["CandidateLabelAreas"]["3101"], .75)
        self.assertAlmostEqual(union["SeamLength"], math.sqrt(2.0))
        # Rotation covariance: the same audit in a rotated frame.
        rotated_reference, rotated_candidate = self._slot_split_cubes((3100, 3101),
                                                                       transform=ROTATION)
        rotated = _protected_surface_report(rotated_reference, rotated_candidate, contract)
        for name in ("PatchCount", "LabelPatchCount", "CoplanarSlotUnions", "TopologyMatches",
                     "PlaneSupportsMatch"):
            self.assertEqual(rotated[name], result[name], name)
        self.assertLessEqual(rotated["MaximumRelativeMeasureError"], 1e-12)
        self.assertLessEqual(rotated["MaximumSupportVertexDistance"], 1e-12)
        rotated_union = next(value["CoplanarSlotUnion"] for value in
                             rotated["ByPlaneSupport"].values() if "CoplanarSlotUnion" in value)
        self.assertEqual(rotated_union["Labels"], union["Labels"])
        self.assertAlmostEqual(rotated_union["SeamLength"], union["SeamLength"])

    def test_coplanar_adjacent_labels_of_different_roles_are_never_unioned(self):
        # Same geometry, but the roles differ by more than the slot index: each
        # label keeps its own 1e-8 measure/boundary check, which the re-partitioned
        # candidate fails.
        for roles in ((("substrate-vacuum-slot-0", [1]), ("conductor-1-slot-0-ms", [1])),
                      (("etched-substrate-vacuum-slot-0", [1]), ("substrate-outer", [1]))):
            contract = self._slot_union_contract({3000: roles[0], 5001: roles[1]})
            reference, candidate = self._slot_split_cubes((3000, 5001))
            result = _protected_surface_report(reference, candidate, contract)
            self.assertTrue(result["PlaneSupportsMatch"], roles)
            self.assertEqual(result["PatchCount"], 7); self.assertEqual(result["CoplanarSlotUnions"], 0)
            self.assertFalse(any("+" in key for key in result["ByPlaneSupport"]))
            self.assertAlmostEqual(result["MaximumRelativeMeasureError"], .5)
            self.assertGreater(result["MaximumSupportVertexDistance"], .1)
            unchanged, _ = self._slot_split_cubes((3000, 5001), split_candidate=False)
            same = _protected_surface_report(reference, unchanged, contract)
            self.assertLessEqual(same["MaximumRelativeMeasureError"], 1e-12)

    def test_slot_labels_meeting_along_a_feature_line_are_not_unioned(self):
        # Two substrate bodies separated by a vertical interface wall: the bottom
        # faces 3100 (left) and 3101 (right) are coplanar and edge-adjacent, but the
        # shared edge carries a third, non-coplanar surface triangle (the wall), so it
        # is a feature line and the slot labels stay separate patches.
        points = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [0, 1, 0],
                           [0, 0, 1], [1, 0, 1], [2, 0, 1], [2, 1, 1], [1, 1, 1], [0, 1, 1]],
                          dtype=float)
        left = [([0, 4, 1], 3100, 1), ([0, 5, 4], 3100, 1),     # bottom
                ([0, 1, 7], 1, 1), ([0, 7, 6], 1, 1), ([4, 5, 11], 1, 1), ([4, 11, 10], 1, 1),
                ([5, 0, 6], 1, 1), ([5, 6, 11], 1, 1), ([6, 7, 10], 1, 1), ([6, 10, 11], 1, 1),
                ([1, 4, 10], 5001, (1, 2)), ([1, 10, 7], 5001, (1, 2))]  # interface wall
        right = [([1, 3, 2], 3101, 2), ([1, 4, 3], 3101, 2),
                 ([1, 2, 8], 1, 2), ([1, 8, 7], 1, 2), ([2, 3, 9], 1, 2), ([2, 9, 8], 1, 2),
                 ([3, 4, 10], 1, 2), ([3, 10, 9], 1, 2), ([7, 8, 9], 1, 2), ([7, 9, 10], 1, 2)]
        contract = self._slot_union_contract({
            3100: ("etched-substrate-vacuum-slot-0", [1]),
            3101: ("etched-substrate-vacuum-slot-1", [2]),
            5001: ("conductor-1-slot-0-ms", [1, 2])}, materials=(1, 2))
        contract["BoundaryLabels"][0]["AdjacentMaterials"] = [1, 2]
        contract["BoundaryLabels"][0]["AdjacentMaterialSets"] = [[1], [2]]
        mesh = self._star_mesh(points, left + right, (1, 2))
        result = _protected_surface_report(mesh, mesh, contract)
        self.assertTrue(result["PlaneSupportsMatch"])
        self.assertEqual(result["CoplanarSlotUnions"], 0)
        self.assertEqual(result["PatchCount"], result["LabelPatchCount"])
        self.assertIn("etched-substrate-vacuum-slot-0 0.0 0.0 1.0 0.0", result["ByPlaneSupport"])
        self.assertIn("etched-substrate-vacuum-slot-1 0.0 0.0 1.0 0.0", result["ByPlaneSupport"])

    def test_same_area_displaced_protected_support_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            self.produce_matrix(root, manifest_path, manifest)
            reference = read_mesh(root / "base--canonical-seed.msh")
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

    def test_protected_footprint_resolves_pinch_vertices_per_wedge(self):
        """Two sub-regions of one label touching at a vertex are a legitimate planar
        topology: each wedge gets its own vertex copy, loops close, geometry is kept."""
        from general_mesh_audit_producer import _split_pinch_vertices
        points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [-1., 0., 0.],
                           [0., -1., 0.], [1., 1., 0.], [-1., 1., 0.], [1., -1., 0.]])
        bow_tie = points[[[0, 1, 2], [0, 3, 4]]]
        pinches = []
        segments, topology = _normalized_footprint_boundary(bow_tie, pinches)
        self.assertEqual(topology, {"Components": 2, "BoundaryLoops": 2, "Holes": 0,
                                    "LoopsPerComponent": [1, 1]})
        self.assertEqual(pinches, [[0., 0., 0.]])
        self.assertEqual(len(segments), 6)
        three_fans = points[[[0, 1, 5], [0, 6, 3], [0, 4, 7]]]
        pinches = []
        _, topology = _normalized_footprint_boundary(three_fans, pinches)
        self.assertEqual(topology, {"Components": 3, "BoundaryLoops": 3, "Holes": 0,
                                    "LoopsPerComponent": [1, 1, 1]})
        self.assertEqual(pinches, [[0., 0., 0.]])
        # The pinch vertex stays a segment endpoint even where the boundary runs
        # straight through it, so the geometry comparison sees it.
        straight = points[[[0, 1, 5], [0, 6, 3]]]
        segments, _ = _normalized_footprint_boundary(straight)
        self.assertTrue(any(np.array_equal(end, [0., 0., 0.])
                            for segment in segments for end in segment))
        diagnostics = {}
        distance, left, right = _footprint_boundary_comparison(bow_tie, bow_tie, diagnostics)
        self.assertEqual(distance, 0.)
        self.assertEqual(left, right)
        self.assertEqual(diagnostics["PinchVertices"],
                         {"Reference": 1, "Candidate": 1,
                          "ReferenceCoordinates": [[0., 0., 0.]],
                          "CandidateCoordinates": [[0., 0., 0.]]})
        # Non-pinched patches are untouched by the split.
        square = np.array([[[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]],
                           [[0., 0., 0.], [1., 1., 0.], [0., 1., 0.]]])
        unique, inverse = np.unique(square.reshape(-1, 3), axis=0, return_inverse=True)
        split_points, split_triangles, coordinates = _split_pinch_vertices(
            unique, inverse.reshape(-1, 3))
        np.testing.assert_array_equal(split_points, unique)
        np.testing.assert_array_equal(split_triangles, inverse.reshape(-1, 3))
        self.assertEqual(coordinates, [])
        segments, topology = _normalized_footprint_boundary(square)
        self.assertEqual(len(segments), 4)
        self.assertEqual(topology, {"Components": 1, "BoundaryLoops": 1, "Holes": 0,
                                    "LoopsPerComponent": [1]})
        # A nonmanifold patch (an edge shared by three triangles) still fails closed.
        nonmanifold = np.array([[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
                                [[0., 0., 0.], [1., 0., 0.], [0., -1., 0.]],
                                [[0., 0., 0.], [1., 0., 0.], [0., 0., 1.]]])
        with self.assertRaisesRegex(ValueError, "Nonmanifold"):
            _normalized_footprint_boundary(nonmanifold)

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
        # A band along a simplified etch footprint edge is a feature band too (a
        # physical dielectric step edge), not a diagonal: aligned with the footprint
        # only -> not counted; aligned with neither -> counted.
        footprint = _global_diagonal_bands(mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025,
                                           footprint_segments=[[[0., 0., 0.], direction]])
        self.assertEqual(footprint["GlobalDiagonalBands"], 0)
        component = footprint["LongShortEdgeComponents"][0]
        self.assertFalse(component["AlignedWithPhysicalSegment"])
        self.assertTrue(component["AlignedWithFootprintSegment"])
        self.assertTrue(component["AlignedWithFeature"])
        self.assertEqual(footprint["FeatureSegments"]["Signature"], 1)
        self.assertEqual(footprint["FeatureSegments"]["Footprint"], 1)
        neither = _global_diagonal_bands(mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025,
                                         footprint_segments=[[[0., 0., 0.], [0., 1., 0.]]])
        self.assertEqual(neither["GlobalDiagonalBands"], 1)
        self.assertFalse(neither["LongShortEdgeComponents"][0]["AlignedWithFeature"])
        # A band along a cut-surface/material-interface junction line carries the same
        # band legitimately: aligned with the junction only -> not counted.
        junction = _global_diagonal_bands(mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025,
                                          footprint_segments=[[[0., 0., 0.], [0., 1., 0.]]],
                                          junction_segments=[[[0., 0., 0.], direction]])
        self.assertEqual(junction["GlobalDiagonalBands"], 0)
        component = junction["LongShortEdgeComponents"][0]
        self.assertFalse(component["AlignedWithPhysicalSegment"])
        self.assertFalse(component["AlignedWithFootprintSegment"])
        self.assertTrue(component["AlignedWithJunctionSegment"])
        self.assertTrue(component["AlignedWithFeature"])
        self.assertEqual(junction["FeatureSegments"]["Junction"], 1)
        self.assertEqual(neither["FeatureSegments"]["Junction"], 0)
        # Decision 21: a band lying ON a bound trace-basis edge (direction aligned and
        # both endpoints within 2 x threshold of the edge segment) is source-driven and
        # reported separately; the same band beside a parallel basis edge, or with no
        # basis bound, is still a diagonal over-refinement.
        unaligned = [[[0., 0., 0.], [1., 0., 0.]]]
        far_edge = [[[-2., 0., 0.], [-2., 0., 0.] + 3 * np.asarray(direction)]]
        on_edge = [[[-0.05, -0.02, 0.], (np.asarray([-0.05, -0.02, 0.]) + 1.1 * np.asarray(direction)).tolist()]]
        accepted = _global_diagonal_bands(mesh, unaligned, .025, trace_basis_edges=on_edge)
        self.assertEqual(accepted["GlobalDiagonalBands"], 0)
        component = accepted["LongShortEdgeComponents"][0]
        self.assertTrue(component["OnTraceBasisEdge"]); self.assertTrue(component["AlignedWithFeature"])
        self.assertFalse(component["AlignedWithPhysicalSegment"])
        self.assertEqual(accepted["TraceBasisEdgeBands"]["Count"], 1)
        self.assertAlmostEqual(accepted["TraceBasisEdgeBands"]["TotalSpan"], component["Span"])
        self.assertEqual(accepted["LineLikeBandsAlignedWith"],
                         {"Signature": 0, "Footprint": 0, "Junction": 0, "TraceBasis": 1})
        self.assertEqual(accepted["FeatureSegments"]["TraceBasis"], 1)
        # Parallel but 2 um away: direction matches, position does not.
        shifted = np.asarray(on_edge) + np.array([0., 2., 0.])
        rejected = _global_diagonal_bands(mesh, unaligned, .025, trace_basis_edges=shifted.tolist())
        self.assertEqual(rejected["GlobalDiagonalBands"], 1)
        self.assertFalse(rejected["LongShortEdgeComponents"][0]["OnTraceBasisEdge"])
        self.assertEqual(rejected["TraceBasisEdgeBands"]["Count"], 0)
        # An edge on the band's line but ending before the band: the segment, not the line.
        short_edge = [[[-0.05, -0.02, 0.], (np.asarray([-0.05, -0.02, 0.]) + 0.4 * np.asarray(direction)).tolist()]]
        self.assertEqual(_global_diagonal_bands(mesh, unaligned, .025,
                                                trace_basis_edges=short_edge)["GlobalDiagonalBands"], 1)
        self.assertEqual(_global_diagonal_bands(mesh, unaligned, .025,
                                                trace_basis_edges=far_edge)["GlobalDiagonalBands"], 1)
        self.assertEqual(neither["FeatureSegments"]["TraceBasis"], 0)
        self.assertEqual(neither["TraceBasisEdgeBands"]["Count"], 0)
        with self.assertRaisesRegex(ValueError, "Degenerate trace basis edge"):
            _global_diagonal_bands(mesh, unaligned, .025, trace_basis_edges=[[[0., 0., 0.], [0., 0., 0.]]])
        # Rotation covariance of the on-edge classification.
        angle = 0.9; q = np.array([[math.cos(angle), -math.sin(angle), 0.],
                                   [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]])
        shift = np.array([3., -1., 2.])
        rotated_mesh = meshio.Mesh(points @ q.T + shift, [("triangle", triangles)],
                                   cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
        rotated = _global_diagonal_bands(
            rotated_mesh, (np.asarray(unaligned) @ q.T + shift).tolist(), .025,
            trace_basis_edges=(np.asarray(on_edge) @ q.T + shift).tolist())
        self.assertEqual(rotated["GlobalDiagonalBands"], 0)
        self.assertTrue(rotated["LongShortEdgeComponents"][0]["OnTraceBasisEdge"])
        rotated_shifted = _global_diagonal_bands(
            rotated_mesh, (np.asarray(unaligned) @ q.T + shift).tolist(), .025,
            trace_basis_edges=(shifted @ q.T + shift).tolist())
        self.assertEqual(rotated_shifted["GlobalDiagonalBands"], 1)
        with self.assertRaisesRegex(ValueError, "Degenerate feature segment"):
            _global_diagonal_bands(mesh, [[[0., 0., 0.], [1., 0., 0.]]], .025,
                                   footprint_segments=[[[0., 0., 0.], [0., 0., 0.]]])

    def test_diagonal_detector_alignment_within_the_band_resolvability(self):
        # Supervisor decision 36: a band's direction is resolvable only to RMSWidth /
        # Span, so alignment with a feature segment is judged within that ratio (the
        # 1e-6 cosine floor kept for degenerate widths); the gate is unchanged.
        def grid(points, triangles, origin, angle, columns, rows, spacing, transverse):
            rotation = np.array([[math.cos(angle), -math.sin(angle)],
                                 [math.sin(angle), math.cos(angle)]])
            base = len(points)
            points.extend([*(rotation @ np.array([index * spacing, offset]) + origin), 0.]
                          for offset in transverse for index in range(columns))
            for row in range(rows - 1):
                for index in range(columns - 1):
                    first = base + row * columns + index; last = first + columns
                    triangles.extend([[first, first + 1, last + 1], [first, last + 1, last]])
        def strip(angle, length):
            # A 10/12.5 nm layer-like band (four rows, 30 nm wide) on a plane whose
            # ordinary triangulation is the 50 nm grid (so the short-edge threshold is
            # the production 0.6 x median = 30 nm), as on a metal face.
            points, triangles = [], []
            columns = int(round(length / .0125)) + 1
            grid(points, triangles, np.zeros(2), angle, columns, 4, .0125, (0., .01, .02, .03))
            coarse = int(round(length / .05)) + 1
            grid(points, triangles, np.array([0., .1]), 0., coarse, 23, .05, [.05 * k for k in range(23)])
            triangles = np.asarray(triangles)
            return meshio.Mesh(np.asarray(points), [("triangle", triangles)],
                               cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
        segment = [[[0., 0., 0.], [1., 0., 0.]]]
        def bands(angle, length):
            report = _global_diagonal_bands(strip(angle, length), segment, .025)
            components = [item for item in report["LongShortEdgeComponents"] if item["LineLike"]]
            self.assertEqual(len(components), 1)
            return report["GlobalDiagonalBands"], components[0]
        # The 1 um edge-layer band tilted by 5 mrad (the ten-edge x = -1 bands: 2.2-2.9
        # mrad at width/span 1.5e-2, rejected by the fixed 1.4 mrad tolerance).
        count, band = bands(.005, 1.)
        self.assertEqual(count, 0); self.assertTrue(band["AlignedWithPhysicalSegment"])
        self.assertAlmostEqual(band["AlignmentAngles"]["Signature"], .005, delta=1e-3)
        self.assertGreater(band["DirectionResolvability"], .005)
        self.assertLess(band["PhysicalSegmentAlignment"], 1 - 1e-6)     # the old rule rejected it
        # A short band whose angle exceeds its resolvability is still rejected ...
        count, band = bands(.02, 1.)
        self.assertEqual(count, 1); self.assertFalse(band["AlignedWithFeature"])
        self.assertLess(band["DirectionResolvability"], math.sin(.02))
        # ... as is a 45-degree diagonal.
        count, band = bands(math.pi / 4, 1.)
        self.assertEqual(count, 1); self.assertAlmostEqual(band["AlignmentAngles"]["Signature"], math.pi / 4, delta=1e-3)
        # Long bands are unchanged: aligned at zero angle; the same 5 mrad tilt exceeds
        # a 10 um band's resolvability (1.1e-3) and is rejected.
        count, band = bands(0., 10.)
        self.assertEqual(count, 0); self.assertTrue(band["AlignedWithPhysicalSegment"])
        self.assertLess(band["AlignmentAngles"]["Signature"], 1e-4)          # within the 1e-6 cosine floor
        count, band = bands(.005, 10.)
        self.assertEqual(count, 1); self.assertFalse(band["AlignedWithPhysicalSegment"])
        self.assertIn("decision 36", _global_diagonal_bands(strip(0., 1.), segment, .025)
                      ["FeatureSegments"]["Alignment"])

    def test_diagonal_detector_threshold_tolerates_construction_roundoff(self):
        # Seed grid edges sit at exactly 2 x NormalSize up to construction roundoff
        # (~1e-12 relative on the coupon).  A band of such edges must be classified
        # the same way whether the roundoff falls below or above the threshold, so
        # a rigid placement cannot split it (the recorded P2 fragility); a
        # geometrically longer edge is still above the threshold.
        columns = 21
        points = np.array([[index * .05, transverse, 0.] for transverse in (-.2, 0., .2)
                           for index in range(columns)])
        triangles = []
        for row in range(2):
            for index in range(columns - 1):
                first = row * columns + index; last = (row + 1) * columns + index
                triangles.extend([[first, first + 1, last + 1], [first, last + 1, last]])
        triangles = np.asarray(triangles)

        def bands(stretch):
            stretched = points * np.array([stretch, 1., 1.])
            mesh = meshio.Mesh(stretched, [("triangle", triangles)],
                               cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
            report = _global_diagonal_bands(mesh, [[[0., 0., 0.], [0., 1., 0.]]], .025)
            self.assertEqual(report["ShortEdgeThreshold"], .05)
            return report["GlobalDiagonalBands"], report["ShortInternalEdges"]

        # The 20 interior row edges (the boundary rows have one owner) form one band.
        self.assertEqual(bands(1.0), (1, 20))
        self.assertEqual(bands(1.0 - 1e-9), (1, 20))
        self.assertEqual(bands(1.0 + 1e-9), (1, 20))
        self.assertEqual(bands(1.0 + 1e-3), (0, 0))

    def test_diagonal_detector_half_diameter_rule_is_rotation_covariant(self):
        # The four-edge rotate-z placement lost its two basis bands because the
        # patch "diameter" was an axis-aligned bounding-box diagonal (22.6 -> 31.6
        # under the 0.63 rad rotation) while the band span (11.44) is invariant.
        # The in-plane convex-hull diameter is a rigid-motion invariant, so the
        # same patch rotated by a non-axis angle gives the same diameter, span and
        # band classification.
        from scipy.spatial import Delaunay
        coarse = [[x, y] for x in np.linspace(0., 1., 11) for y in np.linspace(0., 1., 11)
                  if abs(y - .5) > 1e-9]
        band = [[x, .5] for x in np.arange(0., .8 + 1e-9, .05)]
        planar = np.asarray(coarse + band)
        triangles = Delaunay(planar).simplices
        points = np.column_stack((planar, np.zeros(len(planar))))
        angle = .63; q = np.array([[math.cos(angle), -math.sin(angle), 0.],
                                   [math.sin(angle), math.cos(angle), 0.], [0., 0., 1.]])
        shift = np.array([2., -3., 1.])
        self.assertAlmostEqual(_planar_diameter(points, [0., 0., 1.]), math.sqrt(2.), places=12)
        self.assertAlmostEqual(_planar_diameter(points @ q.T + shift, q @ [0., 0., 1.]),
                               math.sqrt(2.), places=12)
        # The bounding-box diagonal of the rotated square would be 1.98, whose half
        # (0.99) exceeds the 0.8 band span: the band would vanish under rotation.
        self.assertGreater(.5 * float(np.linalg.norm(np.ptp(points @ q.T, axis=0))), .8)
        transverse = [[[0., 0., 0.], [0., 1., 0.]]]
        reports = []
        for rotation, offset in ((np.eye(3), np.zeros(3)), (q, shift)):
            mesh = meshio.Mesh(points @ rotation.T + offset, [("triangle", triangles)],
                               cell_data={"gmsh:physical": [np.ones(len(triangles), int)]})
            reports.append(_global_diagonal_bands(
                mesh, (np.asarray(transverse) @ rotation.T + offset).tolist(), .025))
        for report in reports:
            self.assertEqual(report["ShortEdgeThreshold"], .05)
            self.assertEqual(report["GlobalDiagonalBands"], 1)
            self.assertEqual(len(report["LongShortEdgeComponents"]), 1)
        spans = [report["LongShortEdgeComponents"][0]["Span"] for report in reports]
        self.assertAlmostEqual(spans[0], .8, places=9)
        self.assertAlmostEqual(spans[0], spans[1], places=9)
        self.assertEqual(reports[0]["ShortInternalEdges"], reports[1]["ShortInternalEdges"])

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

    def test_physical_covariance_rejects_transformed_geometry_material_quality_anisotropy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            audits = self.produce_matrix(root, manifest_path, manifest)
            reference = json.loads((audits / "base--identity.json").read_text())
            transformed = json.loads((audits / "base--rotate-z-0.63.json").read_text())
            comparison = manifest["Cases"][0]["TransformComparison"]
            self.assertEqual(_physical_comparison_failures(
                reference, transformed, comparison), [])
            mutations = {
                "protected support geometry/topology": lambda value: value[
                    "PhysicalCovariance"]["ProtectedSurfaces"].__setitem__(
                        "MaximumSupportVertexDistance", 1.0),
                "physical labels/material adjacency": lambda value: value[
                    "PhysicalCovariance"].__setitem__("LabelsMaterialsAdjacencyMatch", False),
                "orientation/quality distribution": lambda value: value[
                    "PhysicalCovariance"]["TransformedQuality"].__setitem__(
                        "PositiveOrientation", False),
                "local-frame anisotropy": lambda value: value[
                    "AchievedAnisotropy"].__setitem__("TangentialP50", 100.0),
            }
            for expected, mutate in mutations.items():
                with self.subTest(expected=expected):
                    candidate = copy.deepcopy(transformed); mutate(candidate)
                    self.assertIn(expected, _physical_comparison_failures(
                        reference, candidate, comparison))

    def test_manifest_rejects_identity_rotation_nonrigid_and_reflection(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest = self.make_suite(root)
            manifest["Cases"][0]["Variants"][1]["Transform"] = IDENTITY
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                validate_manifest(manifest, manifest_path)
        for first_axis in (2, -1):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary); manifest_path, manifest = self.make_suite(root)
                manifest["Cases"][0]["Variants"][1]["Transform"] = list(IDENTITY)
                manifest["Cases"][0]["Variants"][1]["Transform"][0] = first_axis
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

    def test_calibration_manifest_gate_relaxation_never_reaches_production(self):
        # Supervisor decision 22: the MA/MS calibration manifest, and only it, carries
        # MinimumAchievedAspect 0.9 (anisotropy-design gate); the production suite keeps
        # 1.5, has no Calibration block anywhere, and every calibration case is labeled.
        production_path = HERE / "geometry-independence-suite.json"
        calibration_path = HERE / "geometry-independence-calibration-ma.json"
        production = json.loads(production_path.read_text())
        calibration = json.loads(calibration_path.read_text())
        self.assertEqual(production["Gates"]["MinimumAchievedAspect"], 1.5)
        self.assertNotIn("Calibration", production)
        for case in production["Cases"]:
            self.assertNotIn("Calibration", case)
            self.assertNotIn("calib-ma", case["Id"])
            self.assertNotEqual(case["InventoryStatus"], "Calibration")
        deviations = calibration["Calibration"]["GateDeviations"]
        self.assertEqual(set(deviations), {"MinimumAchievedAspect", "EdgeLayerQualityRule",
                                           "MaximumElements"})
        self.assertEqual(deviations["MinimumAchievedAspect"]["Production"], 1.5)
        self.assertEqual(deviations["MinimumAchievedAspect"]["Calibration"], 0.9)
        self.assertIn("FORBIDDEN", deviations["MinimumAchievedAspect"]["ProductionUse"])
        self.assertEqual(calibration["Gates"]["MinimumAchievedAspect"], 0.9)
        # Decision 32: the layer-local quality rule is a calibration-only gate; the
        # production suite has no such gate and no case with a layer.
        self.assertNotIn("EdgeLayerQualityRule", production["Gates"])
        rule = calibration["Gates"]["EdgeLayerQualityRule"]
        self.assertEqual(rule["MaximumEdgeAspect"], 100.0)   # 2 x lc_tangent / EdgeSize
        self.assertEqual(rule["ScaledJacobianRoundoffFloor"], 1e-12)
        self.assertIsNone(deviations["EdgeLayerQualityRule"]["Production"])
        self.assertEqual(deviations["EdgeLayerQualityRule"]["Calibration"], 100.0)
        self.assertIn("FORBIDDEN", deviations["EdgeLayerQualityRule"]["ProductionUse"])
        deviating = {"MinimumAchievedAspect", "EdgeLayerQualityRule"}
        self.assertEqual({key: value for key, value in calibration["Gates"].items()
                          if key not in deviating},
                         {key: value for key, value in production["Gates"].items()
                          if key not in deviating})
        # Decision 33: the element cap 5,000,000 is a per-case calibration deviation
        # of the 1 nm aspect-4 case only; both manifests' Gates keep 4,000,000 (the
        # production value and the build cache key) and production cases never
        # carry a cap.
        self.assertEqual(production["Gates"]["MaximumElements"], 4000000)
        self.assertEqual(calibration["Gates"]["MaximumElements"], 4000000)
        cap = deviations["MaximumElements"]
        self.assertEqual((cap["Production"], cap["Calibration"]), (4000000, 5000000))
        self.assertEqual(cap["Cases"], ["four-edge-calib-ma-el1c"])
        self.assertIn("FORBIDDEN", cap["ProductionUse"])
        for case in production["Cases"]:
            self.assertEqual(case_gates(production, case)["MaximumElements"], 4000000)
        for case in calibration["Cases"]:
            expected = 5000000 if case["Id"] in cap["Cases"] else 4000000
            self.assertEqual(case_gates(calibration, case)["MaximumElements"], expected)
            self.assertEqual(case["Calibration"].get("MaximumElements"),
                             5000000 if case["Id"] in cap["Cases"] else None)
        capped = next(case for case in calibration["Cases"] if case["Id"] == cap["Cases"][0])
        over = {"Resources": {"ExitCode": 0, "Seconds": 1., "PeakRSSGiB": 1., "Elements": 4500000,
                              "CanonicalBuild": {"Seconds": 1., "PeakRSSGiB": 1.},
                              "PlacementPublication": {"Seconds": 1., "PeakRSSGiB": 1.}}}
        contract = load_semantic_contract(
            HERE / capped["Source"]["Directory"].split("spatial_coupon/", 1)[1] / "semantic-contract.json")
        binding = {"CaseId": capped["Id"], "Variant": "identity",
                   "Transform": capped["Variants"][0]["Transform"], "TransformSHA256": "x",
                   "InputSHA256": {"Process": "p", "SemanticContract": "s", "MeshRecipe": "r"},
                   "ToolSHA256": {}, "StageToolSHA256": {}}
        self.assertNotIn("bounded-resources", audit_manifest_evidence(
            over, case_gates(calibration, capped), contract, binding))
        self.assertIn("bounded-resources", audit_manifest_evidence(
            over, case_gates(calibration, calibration["Cases"][0]), contract, binding))
        self.assertIn("bounded-resources", audit_manifest_evidence(
            over, production["Gates"], contract, binding))
        with tempfile.TemporaryDirectory() as temporary:
            forged = Path(temporary) / "geometry-independence-suite.json"
            def rejected(mutate, message):
                manifest = copy.deepcopy(production); mutate(manifest)
                forged.write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, message):
                    validate_manifest(manifest, forged, check_available_files=False)
            # A production manifest cannot carry the rule ...
            rejected(lambda m: m["Gates"].__setitem__("EdgeLayerQualityRule", rule),
                     "production manifest cannot carry")
            # ... nor a case with an EdgeLayer or Calibration block (preflight fails).
            rejected(lambda m: m["Cases"][0].__setitem__("EdgeLayer", {"EdgeSize": .001}),
                     "calibration or edge-layer block in a production manifest")
            rejected(lambda m: m["Cases"][0].__setitem__("Calibration",
                                                         calibration["Cases"][0]["Calibration"]),
                     "calibration or edge-layer block in a production manifest")
            # In the calibration manifest the rule must be labeled and consistent.
            def rejected_calibration(mutate, message):
                manifest = copy.deepcopy(calibration); mutate(manifest)
                forged.write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, message):
                    validate_manifest(manifest, forged, check_available_files=False)
            rejected_calibration(lambda m: m["Calibration"]["GateDeviations"].pop("EdgeLayerQualityRule"),
                                 "not labeled")
            rejected_calibration(lambda m: m["Gates"]["EdgeLayerQualityRule"].__setitem__(
                "MaximumEdgeAspect", 90.), "not labeled")
            rejected_calibration(lambda m: m["Gates"]["EdgeLayerQualityRule"].__setitem__(
                "ScaledJacobianRoundoffFloor", .01), "not labeled")
            rejected_calibration(lambda m: m["Cases"][-1].__setitem__("EdgeLayer", {}),
                                 "must declare its edge layer under Calibration")
            # The element cap must be labeled: deviation present, values equal, case
            # named, above the manifest gate; the label names exactly the declaring cases.
            def capped_case(m):
                return next(case for case in m["Cases"] if case["Id"] == cap["Cases"][0])
            rejected_calibration(lambda m: m["Calibration"]["GateDeviations"].pop("MaximumElements"),
                                 "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: capped_case(m)["Calibration"].__setitem__("MaximumElements", 6000000),
                                 "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: capped_case(m)["Calibration"].__setitem__("MaximumElements", 3000000),
                                 "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: m["Calibration"]["GateDeviations"]["MaximumElements"].__setitem__(
                "Production", 5000000), "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: m["Calibration"]["GateDeviations"]["MaximumElements"].__setitem__(
                "Cases", []), "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: m["Cases"][0]["Calibration"].__setitem__("MaximumElements", 5000000),
                                 "not labeled as a calibration-only deviation")
            rejected_calibration(lambda m: capped_case(m)["Calibration"].pop("MaximumElements"),
                                 "name exactly the cases declaring the cap")
            rejected(lambda m: m["Cases"][0].__setitem__("Calibration", {"MaximumElements": 5000000}),
                     "calibration or edge-layer block in a production manifest")
        self.assertEqual(calibration["Tools"], production["Tools"])
        self.assertEqual(calibration["StageToolSHA256"], production["StageToolSHA256"])
        base = next(item for item in production["Cases"] if item["Id"] == "four-edge-9d2cb9bbb3fe")
        for case in calibration["Cases"]:
            self.assertIn("calib-ma", case["Id"])
            self.assertEqual(case["InventoryStatus"], "Calibration")
            self.assertEqual(case["Calibration"]["BaseCase"], base["Id"])
            self.assertEqual(case["Source"], base["Source"])
            self.assertEqual(case["Variants"], base["Variants"])
            self.assertEqual(case["TransformComparison"], base["TransformComparison"])
        validate_manifest(calibration, calibration_path)

    def test_refreeze_manifest_tools_keeps_both_manifests_current_together(self):
        # The in-repo refreeze recomputes every repository-tool digest of the production
        # manifest and mirrors Tools / StageToolSHA256 into the calibration manifest; the
        # committed manifests must be current (fails closed when a refreeze was forgotten).
        import refreeze_manifest_tools as refreezer
        production_path = HERE / "geometry-independence-suite.json"
        calibration_path = HERE / "geometry-independence-calibration-ma.json"
        self.assertEqual(refreezer.refreeze(production_path, calibration_path, check_only=True),
                         ([], False))
        production = json.loads(production_path.read_text())
        for (stage, role), name in refreezer.STAGE_REPOSITORY_TOOLS.items():
            self.assertEqual(production["StageToolSHA256"][stage][role], sha256(HERE / name),
                             f"{stage}/{role}")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "examples" / "cpw3d_surface" / "spatial_coupon"
            root.mkdir(parents=True)
            for name in [tool["Path"].rsplit("/", 1)[1] for tool in production["Tools"]] + list(
                    refreezer.STAGE_REPOSITORY_TOOLS.values()):
                shutil.copy(HERE / name, root / name)
            stale = copy.deepcopy(production)
            stale["Tools"][0]["SHA256"] = "0" * 64
            stale["StageToolSHA256"]["seed-generation"]["mesher"] = "1" * 64
            stale_production = root / "geometry-independence-suite.json"
            stale_calibration = root / "geometry-independence-calibration-ma.json"
            stale_production.write_text(json.dumps(stale))
            # The calibration copy still mirrors the (current) production digests.
            stale_calibration.write_text(json.dumps(json.loads(calibration_path.read_text())))
            changes, mirror_stale = refreezer.refreeze(stale_production, stale_calibration,
                                                      check_only=True)
            self.assertEqual([(name, old) for name, old, _ in changes],
                             [(production["Tools"][0]["Name"], "0" * 64),
                              ("seed-generation/mesher", "1" * 64)])
            self.assertFalse(mirror_stale)
            self.assertEqual(json.loads(stale_production.read_text()), stale)  # check-only
            changes, _ = refreezer.refreeze(stale_production, stale_calibration,
                                            check_only=False)
            self.assertEqual(len(changes), 2)
            refrozen = json.loads(stale_production.read_text())
            self.assertEqual(refrozen, production)
            mirrored = json.loads(stale_calibration.read_text())
            self.assertEqual(mirrored["Tools"], production["Tools"])
            self.assertEqual(mirrored["StageToolSHA256"], production["StageToolSHA256"])
            self.assertEqual(mirrored["Gates"]["MinimumAchievedAspect"], 0.9)
            self.assertEqual(refreezer.refreeze(stale_production, stale_calibration,
                                                check_only=True), ([], False))
            # A calibration mirror that drifted from an otherwise current production
            # manifest is stale on its own and is restored by the refreeze.
            drifted = json.loads(stale_calibration.read_text())
            drifted["StageToolSHA256"]["native-adaptation-mmg"]["adapter-mmg"] = "2" * 64
            stale_calibration.write_text(json.dumps(drifted))
            self.assertEqual(refreezer.refreeze(stale_production, stale_calibration,
                                                check_only=True), ([], True))
            refreezer.refreeze(stale_production, stale_calibration, check_only=False)
            self.assertEqual(json.loads(stale_calibration.read_text())["StageToolSHA256"],
                             production["StageToolSHA256"])
            # Machine-bound identities are never recomputed.
            for stage, roles in production["StageToolSHA256"].items():
                for role in roles:
                    if (stage, role) not in refreezer.STAGE_REPOSITORY_TOOLS:
                        self.assertEqual(json.loads(stale_production.read_text())
                                         ["StageToolSHA256"][stage][role], roles[role])
            # The reviewed adapter is machine-bound: its digest changes only through an
            # explicit --adapter-mmg naming the executable of the recorded build whose
            # source digest is the repository's adapt_edge_metric.cpp.
            (root / "testdata").mkdir()
            shutil.copy(HERE / "adapt_edge_metric.cpp", root / "adapt_edge_metric.cpp")
            adapter = root / "adapt_edge_metric"; adapter.write_bytes(b"adapter build")
            record = {"ExecutableSHA256": sha256(adapter),
                      "SourceSHA256": sha256(root / "adapt_edge_metric.cpp")}
            (root / "testdata" / "adapter-build.json").write_text(json.dumps(record))
            changes, _ = refreezer.refreeze(stale_production, stale_calibration, check_only=False,
                                            adapter=adapter)
            self.assertEqual(changes, [("native-adaptation-mmg/adapter-mmg",
                                        production["StageToolSHA256"]["native-adaptation-mmg"]
                                        ["adapter-mmg"], sha256(adapter))])
            for path in (stale_production, stale_calibration):
                self.assertEqual(json.loads(path.read_text())["StageToolSHA256"]
                                 ["native-adaptation-mmg"]["adapter-mmg"], sha256(adapter))
            other = root / "other_adapter"; other.write_bytes(b"unrecorded build")
            with self.assertRaisesRegex(ValueError, "differ from the recorded build"):
                refreezer.refreeze(stale_production, stale_calibration, check_only=True,
                                   adapter=other)
            (root / "adapt_edge_metric.cpp").write_text("// edited source\n")
            with self.assertRaisesRegex(ValueError, "differ from the recorded build"):
                refreezer.refreeze(stale_production, stale_calibration, check_only=True,
                                   adapter=adapter)
            # The Julia launcher is machine-bound likewise: --julia-runtime refreezes the
            # three Julia runtime roles together.
            launcher = root / "julialauncher"; launcher.write_bytes(b"julia launcher")
            changes, _ = refreezer.refreeze(stale_production, stale_calibration, check_only=False,
                                            julia_runtime=launcher)
            self.assertEqual(sorted(name for name, _, _ in changes),
                             sorted(f"{stage}/{role}" for stage, role in refreezer.JULIA_RUNTIME_ROLES))
            for stage, role in refreezer.JULIA_RUNTIME_ROLES:
                self.assertEqual(json.loads(stale_calibration.read_text())["StageToolSHA256"]
                                 [stage][role], sha256(launcher))
            with self.assertRaisesRegex(ValueError, "does not exist"):
                refreezer.refreeze(stale_production, stale_calibration, check_only=True,
                                   julia_runtime=root / "missing")
            # The committed record names the committed source and the frozen adapter digest.
            committed = json.loads((HERE / "testdata" / "adapter-build.json").read_text())
            self.assertEqual(committed["SourceSHA256"], sha256(HERE / "adapt_edge_metric.cpp"))
            self.assertEqual(committed["ExecutableSHA256"],
                             production["StageToolSHA256"]["native-adaptation-mmg"]["adapter-mmg"])
            self.assertEqual(committed["MMGLibrarySHA256"],
                             production["StageToolSHA256"]["native-adaptation-mmg"]["mmg-library"])

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


class AchievedAnisotropyDesignGateTest(unittest.TestCase):
    """Decisions 34B/35: the production manifest records the EL4c recipe; the
    achieved-anisotropy design gate judges the one-NormalSize band outside a recorded
    edge layer and is not applicable by construction on a layer-covered band, whose
    layer-adjacent band is recorded and never gated."""

    def setUp(self):
        self.production = json.loads((HERE / "geometry-independence-suite.json").read_text())
        self.gates = self.production["Gates"]
        case = self.production["Cases"][0]
        self.contract = load_semantic_contract(HERE / "testdata" / "one-edge-semantic.json")
        self.binding = {"CaseId": case["Id"], "Variant": "identity",
                        "Transform": case["Variants"][0]["Transform"],
                        "TransformSHA256": "x", "InputSHA256": {"Process": "p", "SemanticContract": "s",
                                                                "MeshRecipe": "r"},
                        "ToolSHA256": {}, "StageToolSHA256": {}}

    @staticmethod
    def covered(**overrides):
        # The EL4c identity measurement (decision 35 record).
        adjacent = {"Cells": 5955, "DistanceCutoff": .025 * 3.,
                    "NearestSpanVertexDistance": {"Minimum": .0334, "Maximum": .083},
                    "TangentialP50": .05, "Transverse1P90": .0405, "Transverse2P90": .0665,
                    "TransverseP90OverNormalSize": .0665 / .025, "Rule": LAYER_ADJACENT_BAND_RULE}
        layer = {"Cells": 26870, "Reach": .0334, "EdgeSize": .004, "Aspect": 4.,
                 "TangentialP50": .0187, "Transverse1P90": .0246, "Transverse2P90": .0469, "Rule": "r"}
        record = {"Samples": 0, "DistanceCutoff": .025, "NormalTarget": .025, "TangentialP50": None,
                  "Transverse1P90": None, "Transverse2P90": None, "ExcludedEdgeLayerCells": 26870,
                  "Gate": ANISOTROPY_GATE_NOT_APPLICABLE, "LayerAdjacentBand": adjacent,
                  "EdgeLayer": layer}
        record.update(overrides)
        return record

    def failures(self, widths, gates=None):
        names = audit_manifest_evidence({"AchievedAnisotropy": widths}, gates or self.gates,
                                        self.contract, self.binding)
        return "achieved-anisotropy" in names

    def test_production_manifest_records_the_el4c_recipe(self):
        recipe = validate_production_recipe(self.production)
        self.assertEqual(recipe["SeedCommandOptions"],
                         {"--lc-tangent": .05, "--edge-size": .004, "--edge-growth-ratio": 2.,
                          "--edge-layer-aspect": 4., "--corner-size": .004})
        self.assertEqual(recipe["MetricCommandOptions"],
                         {"--far-growth": .5, "--edge-size": .004, "--edge-growth-ratio": 2.,
                          "--edge-layer-aspect": 4., "--corner-size": .004})
        self.assertEqual(recipe["AdaptationCommandOptions"], {"--hmin": .004})
        # The physical gates and the design gate values are unchanged.
        self.assertEqual(self.gates["MinimumAchievedAspect"], 1.5)
        self.assertEqual(self.gates["MaximumNormalFactor"], 2.0)
        self.assertEqual(self.gates["MinimumScaledJacobian"], .01)
        self.assertEqual(self.gates["MaximumJacobianCondition"], 1000.)
        self.assertEqual(self.gates["MaximumCornerAspect"], 4.)
        self.assertEqual(self.gates["MaximumProtectedMeasureError"], 1e-8)
        self.assertEqual(self.gates["MaximumElements"], 4000000)
        self.assertEqual((self.gates["MaximumSeconds"], self.gates["MaximumRSSGiB"]), (1800, 8.))
        # The calibration manifest never carries the block and its cases keep their
        # pre-34B production values; EL4c is labeled as the adopted recipe.
        calibration = json.loads((HERE / "geometry-independence-calibration-ma.json").read_text())
        self.assertNotIn("ProductionRecipe", calibration)
        with self.assertRaisesRegex(ValueError, "cannot carry a production recipe"):
            validate_production_recipe({**calibration, "ProductionRecipe": recipe})
        el4c = next(case for case in calibration["Cases"] if case["Id"] == "four-edge-calib-ma-el4c")
        self.assertIn("AdoptedAsProductionRecipe", el4c["Calibration"])
        for case in calibration["Cases"]:
            self.assertNotIn("ProductionValues", case["Calibration"])
            self.assertEqual(case["Calibration"]["ProductionValuesBefore34B"]["--lc-tangent"], .1)
        self.assertEqual(el4c["Calibration"]["ProductionValuesBefore34B"]["--edge-size"], 0.)
        self.assertEqual(el4c["Calibration"]["SeedCommandOptions"]["--lc-tangent"],
                         recipe["SeedCommandOptions"]["--lc-tangent"])
        self.assertEqual(el4c["Calibration"]["MetricCommandOptions"]["--far-growth"],
                         recipe["MetricCommandOptions"]["--far-growth"])
        self.assertEqual(el4c["Calibration"]["AdaptationCommandOptions"], recipe["AdaptationCommandOptions"])
        # Malformed blocks are rejected.
        for broken in ({**recipe, "AdaptationCommandOptions": {}},
                       {**recipe, "SeedCommandOptions": {"--lc-tangent": "0.05"}},
                       {**recipe, "MetricCommandOptions": {"far-growth": .5}}, []):
            with self.assertRaisesRegex(ValueError, "Production recipe must bind"):
                validate_production_recipe({**self.production, "ProductionRecipe": broken})

    def test_production_recipe_commands_are_bound_exactly_once(self):
        stages = {"seed-generation": {"Command": ["julia", "seed.jl", "--lc-tangent", ".05",
                                                  "--edge-size", "0.004", "--edge-growth-ratio", "2",
                                                  "--edge-layer-aspect", "4", "--corner-size", ".004"]},
                  "metric-preparation": {"Command": ["python3", "metric.py", "--far-growth", "0.5",
                                                     "--edge-size", ".004", "--edge-growth-ratio", "2.0",
                                                     "--edge-layer-aspect", "4.0", "--corner-size", "0.004"]},
                  "native-adaptation-mmg": {"Command": ["python3", "adapt.py", "--hmin", ".004"]}}
        case = self.production["Cases"][0]
        validate_production_recipe_commands(self.production, case, stages)
        validate_production_recipe_commands({}, case, stages)                 # no recipe block
        validate_production_recipe_commands(self.production, {**case, "Calibration": {}}, stages)
        for stage, command in (("seed-generation", ["julia", "seed.jl", "--lc-tangent", ".1",
                                                    "--edge-size", "0.004", "--edge-growth-ratio", "2",
                                                    "--edge-layer-aspect", "4", "--corner-size", ".004"]),
                               ("metric-preparation", ["python3", "metric.py", "--edge-size", ".004",
                                                       "--edge-growth-ratio", "2", "--edge-layer-aspect",
                                                       "4", "--corner-size", "0.004"]),
                               ("native-adaptation-mmg", ["python3", "adapt.py", "--hmin", ".025"]),
                               ("native-adaptation-mmg", ["python3", "adapt.py", "--hmin", ".004",
                                                          "--hmin", ".004"])):
            with self.assertRaisesRegex(ValueError, f"{stage} command does not execute the production recipe"):
                validate_production_recipe_commands(self.production, case, {**stages, stage: {"Command": command}})

    def test_layer_covered_band_is_not_applicable_and_the_adjacent_band_is_informational(self):
        self.assertTrue(layer_covered_band(self.covered()))
        self.assertFalse(self.failures(self.covered()))
        # The layer-adjacent band's values do not gate: aspect 0.75 and transverse P90
        # 2.66 x NormalSize are recorded, not judged.
        loose = self.covered(); loose["LayerAdjacentBand"] = dict(
            loose["LayerAdjacentBand"], TangentialP50=.01, Transverse1P90=.2, Transverse2P90=.3,
            TransverseP90OverNormalSize=.3 / .025)
        self.assertFalse(self.failures(loose))
        # Negatives: anything less than the complete not-applicable record is judged as
        # an ordinary band sample and fails on Samples 0 or the gate values.
        negatives = {
            "band cells outside the layer": {"Samples": 5, "TangentialP50": .05, "Transverse1P90": .04,
                                             "Transverse2P90": .066},
            "no layer cells": {"ExcludedEdgeLayerCells": 0},
            "gate not declared": {"Gate": ANISOTROPY_GATE_APPLIED},
            "no layer record": {"EdgeLayer": None},
            "no adjacent band": {"LayerAdjacentBand": None},
            "adjacent band without cells": {"LayerAdjacentBand": dict(self.covered()["LayerAdjacentBand"],
                                                                      Cells=0)},
            "adjacent band statistic missing": {"LayerAdjacentBand": dict(
                self.covered()["LayerAdjacentBand"], Transverse2P90=None)},
            "normal factor inconsistent": {"LayerAdjacentBand": dict(
                self.covered()["LayerAdjacentBand"], TransverseP90OverNormalSize=2.)},
            "adjacent band rule": {"LayerAdjacentBand": dict(self.covered()["LayerAdjacentBand"], Rule="x")},
            "adjacent band cutoff": {"LayerAdjacentBand": dict(self.covered()["LayerAdjacentBand"],
                                                               DistanceCutoff=.05)},
            "gated sample cutoff": {"DistanceCutoff": .075},
            "layer without cells": {"EdgeLayer": dict(self.covered()["EdgeLayer"], Cells=0)},
            "no span distance": {"LayerAdjacentBand": dict(self.covered()["LayerAdjacentBand"],
                                                           NearestSpanVertexDistance=None)},
        }
        for name, overrides in negatives.items():
            with self.subTest(name=name):
                record = self.covered(**overrides)
                self.assertFalse(layer_covered_band(record))
                self.assertTrue(self.failures(record))
        # A band with cells outside the layer within one NormalSize is judged by 1.5 and
        # the normal factor exactly as before (production four-edge: 0.0946 / 0.0500).
        judged = {"Samples": 2103, "DistanceCutoff": .025, "NormalTarget": .025, "TangentialP50": .0946,
                  "Transverse1P90": .0417, "Transverse2P90": .05, "ExcludedEdgeLayerCells": 0,
                  "Gate": ANISOTROPY_GATE_APPLIED, "LayerAdjacentBand": None, "EdgeLayer": None}
        self.assertFalse(self.failures(judged))
        self.assertTrue(self.failures({**judged, "TangentialP50": .057}))        # V2: 1.14 < 1.5
        self.assertTrue(self.failures({**judged, "Transverse2P90": .0665}))       # > 2 x NormalSize
        self.assertTrue(self.failures({**judged, "Gate": None}))
        self.assertTrue(self.failures({**judged, "Gate": ANISOTROPY_GATE_NOT_APPLICABLE}))
        self.assertTrue(self.failures({**judged, "Samples": 0}))

    def test_covariance_compares_the_adjacent_band_when_the_gate_is_not_applicable(self):
        comparison = {"MaximumRelativeVolumeError": 1e-8, "MaximumRelativeSurfaceMeasureError": 1e-8,
                      "MaximumProtectedSupportHausdorff": 1e-8, "MaximumProtectedMeasureError": 1e-8,
                      "MaximumQualityDistributionRelativeError": 1e-8,
                      "MaximumAnisotropyRelativeError": .35, "MaximumComplexityRatio": 1.1}
        physical = {"LabelsMaterialsAdjacencyMatch": True,
                    "ReferenceInvariants": {"Volume:1": 1., "Area:1": 1.},
                    "TransformedInvariants": {"Volume:1": 1., "Area:1": 1.},
                    "ProtectedSurfaces": {"PlaneSupportsMatch": True, "TopologyMatches": True,
                                          "MaximumSupportVertexDistance": 0., "MaximumRelativeMeasureError": 0.},
                    "ReferenceQuality": {"PositiveOrientation": True, "ScaledJacobianQuantiles": [.1],
                                         "JacobianConditionQuantiles": [2.]},
                    "TransformedQuality": {"PositiveOrientation": True, "ScaledJacobianQuantiles": [.1],
                                           "JacobianConditionQuantiles": [2.]}}
        def evidence(widths):
            return {"PhysicalCovariance": physical, "AchievedAnisotropy": widths,
                    "Complexity": {"H1DOFs": 100}, "Resources": {"Elements": 50}}
        reference, transformed = evidence(self.covered()), evidence(self.covered())
        self.assertEqual(_physical_comparison_failures(reference, transformed, comparison), [])
        # The adjacent band is what covariance compares ...
        rotated = self.covered(); rotated["LayerAdjacentBand"] = dict(
            rotated["LayerAdjacentBand"], TangentialP50=.1)
        self.assertEqual(_physical_comparison_failures(reference, evidence(rotated), comparison),
                         ["local-frame anisotropy"])
        # ... and the two variants must agree on whether the gate applied.
        applied = self.covered(Samples=5, TangentialP50=.05, Transverse1P90=.04, Transverse2P90=.066,
                               Gate=ANISOTROPY_GATE_APPLIED)
        self.assertEqual(_physical_comparison_failures(reference, evidence(applied), comparison),
                         ["local-frame anisotropy"])


class EdgeLayerQualityGateTest(unittest.TestCase):
    """Decision 32: under Gates.EdgeLayerQualityRule the recorded layer's cells are
    judged by orientation and edge aspect, every other cell by the production
    quality gates; without the rule the whole mesh is judged as before."""

    def setUp(self):
        manifest = json.loads((HERE / "geometry-independence-calibration-ma.json").read_text())
        self.gates = manifest["Gates"]
        self.production_gates = json.loads(
            (HERE / "geometry-independence-suite.json").read_text())["Gates"]
        case = manifest["Cases"][0]
        self.contract = load_semantic_contract(
            HERE / case["Source"]["Directory"].split("spatial_coupon/", 1)[1] / "semantic-contract.json")
        self.binding = {"CaseId": case["Id"], "Variant": "identity",
                        "Transform": case["Variants"][0]["Transform"],
                        "TransformSHA256": "x", "InputSHA256": {"Process": "p", "SemanticContract": "s",
                                                                "MeshRecipe": "r"},
                        "ToolSHA256": {}, "StageToolSHA256": {}}

    def quality(self, layer_scaled=4e-4, layer_aspect=70.7, layer_oriented=True,
                outside_scaled=.02, outside_condition=200.):
        from edge_volume_metric import EDGE_LAYER_QUALITY_RULE
        return {"Samples": 1000, "PositiveOrientation": layer_oriented,
                "MinimumScaledJacobian": min(layer_scaled, outside_scaled),
                "MaximumJacobianCondition": max(outside_condition, 900.),
                "OutsideEdgeLayer": {"Samples": 900, "PositiveOrientation": True,
                                     "MinimumScaledJacobian": outside_scaled,
                                     "MaximumJacobianCondition": outside_condition},
                "EdgeLayer": {"Cells": 100, "Reach": .03355, "PositiveOrientation": layer_oriented,
                              "MinimumScaledJacobian": layer_scaled, "MinimumDeterminant": 1e-12,
                              "MaximumEdgeAspect": layer_aspect, "Rule": EDGE_LAYER_QUALITY_RULE}}

    def failures(self, quality, gates=None):
        names = audit_manifest_evidence({"MeshQuality": quality}, gates or self.gates,
                                        self.contract, self.binding)
        return {name for name in names if name in ("edge-layer-quality", "mesh-quality-jacobian")}

    def test_layer_rule_judges_layer_cells_and_production_gates_judge_the_rest(self):
        # A 1 nm x 50 nm layer (scaled Jacobian 4e-4, aspect 70.7) passes under the rule ...
        self.assertEqual(self.failures(self.quality()), set())
        # ... and fails the whole-mesh scaled-Jacobian gate under the production gates.
        self.assertEqual(self.failures(self.quality(), self.production_gates),
                         {"mesh-quality-jacobian"})
        # Negatives: a layer cell with negative orientation; a layer cell flat to
        # roundoff; an aspect above the bound; a non-layer cell below 0.01; a
        # non-layer condition above the gate; inconsistent sample counts.
        self.assertEqual(self.failures(self.quality(layer_oriented=False)), {"edge-layer-quality"})
        self.assertEqual(self.failures(self.quality(layer_scaled=1e-13)), {"edge-layer-quality"})
        self.assertEqual(self.failures(self.quality(layer_aspect=100.5)), {"edge-layer-quality"})
        self.assertEqual(self.failures(self.quality(layer_aspect=100.)), set())
        self.assertEqual(self.failures(self.quality(outside_scaled=.009)), {"mesh-quality-jacobian"})
        self.assertEqual(self.failures(self.quality(outside_condition=1001.)), {"mesh-quality-jacobian"})
        inconsistent = self.quality(); inconsistent["OutsideEdgeLayer"]["Samples"] = 899
        self.assertEqual(self.failures(inconsistent), {"edge-layer-quality"})
        # A mesh without a recorded layer is judged whole under either gate set.
        plain = self.quality(); plain["EdgeLayer"] = None; plain["MinimumScaledJacobian"] = .02
        plain["MaximumJacobianCondition"] = 200.
        self.assertEqual(self.failures(plain), set())
        self.assertEqual(self.failures(plain, self.production_gates), set())

    def test_rule_is_case_keyed_on_the_calibration_declaration(self):
        """case_gates: the rule judges only a case declaring
        Calibration.EdgeLayerQualityRule (its stages executed the bound); the 4 nm
        layer case, which declares none, is judged by MinimumScaledJacobian on its
        whole mesh (its record), and a production case never sees the rule."""
        manifest = json.loads((HERE / "geometry-independence-calibration-ma.json").read_text())
        declared = next(case for case in manifest["Cases"]
                        if case["Calibration"].get("EdgeLayerQualityRule") is not None)
        undeclared = next(case for case in manifest["Cases"]
                          if case["Id"] == "four-edge-calib-ma-edge-layer-4nm")
        self.assertIn("EdgeLayer", undeclared["Calibration"])
        self.assertNotIn("EdgeLayerQualityRule", undeclared["Calibration"])
        self.assertEqual(case_gates(manifest, declared), manifest["Gates"])
        without = case_gates(manifest, undeclared)
        self.assertNotIn("EdgeLayerQualityRule", without)
        self.assertEqual({key: value for key, value in manifest["Gates"].items()
                          if key != "EdgeLayerQualityRule"}, without)
        # The same 4e-4 layer passes for the declared case and fails the whole-mesh
        # scaled-Jacobian gate for the undeclared one (never "edge-layer-quality").
        self.assertEqual(self.failures(self.quality(), case_gates(manifest, declared)), set())
        self.assertEqual(self.failures(self.quality(), without), {"mesh-quality-jacobian"})
        production = json.loads((HERE / "geometry-independence-suite.json").read_text())
        for case in production["Cases"]:
            self.assertEqual(case_gates(production, case), production["Gates"])
        # The manifest Gates (the build cache key) are not mutated.
        self.assertIn("EdgeLayerQualityRule", manifest["Gates"])


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


class RequiredRegionContractTest(unittest.TestCase):
    """The metric stage's required tetrahedra are the recipe's corner balls and edge
    layer, and the seed stage gated the same region with the restorer's values."""

    GATES = ["--maximum-corner-aspect", "4", "--minimum-scaled-jacobian", ".01",
             "--maximum-jacobian-condition", "1000",
             "--maximum-quality-displacement-over-normal", ".75"]

    def fixture(self, root, layer=True):
        corners = [[0., 0., 0.], [10., 0., 0.]]
        spans = [[1., 0., .1, 9., 0., .1]]
        recipe = {"Tetrahedra": 100, "CornerIsotropyRadius": .1, "TruePhysicalCorners": corners,
                  "RequiredTetrahedra": {"Count": 5 if layer else 3, "CornerRadius": .1, "IndexBase": 1,
                      "PerCorner": [{"Point": corners[0], "Tetrahedra": 2},
                                    {"Point": corners[1], "Tetrahedra": 1}],
                      "LayerRequiredReach": .028 * 1.05 + .004 if layer else None,
                      "PerSpan": [{"Span": spans[0], "Tetrahedra": 3}] if layer else []}}
        if layer:
            recipe["EdgeLayer"] = {"EdgeSize": .004, "LayerThickness": .028, "RowZigzag": .05,
                                   "RequiredReach": .028 * 1.05 + .004,
                                   "Spans": spans, "SpanCount": 1}
        # The census counts the set recomputed on the final seed positions: exactly
        # the recipe count (the pre-move set may differ and is reported).
        census = {"SeedQualityOptimization": {
            "MaximumCornerAspect": 4., "MinimumScaledJacobian": .01,
            "MaximumJacobianCondition": 1000.,
            "DisplacementBoundOverNormal": .75, "RequiredTetrahedra": 5 if layer else 3,
            "RequiredTetrahedraBeforeMoves": 6 if layer else 3,
            "RequiredCellsBelowGateAfter": 0, "CornerAspectsAfter": [3.4, 3.74],
            "RequiredMinimumScaledJacobianAfter": .02,
            "RequiredCellsAboveConditionAfter": 0,
            "RequiredMaximumJacobianConditionAfter": 412.}}
        seed = {"Command": ["julia", "mesh_spatial_coupon.jl", "sig.csv", "fabricated", "seed.msh",
                            *self.GATES]}
        restored = root / f"restored-{layer}.msh"
        restored.with_suffix(".projection.json").write_text(json.dumps(
            {"RequiredTetrahedra": 5 if layer else 3, "RequiredVertices": 12,
             "RequiredMaximumJacobianCondition": 412.}))
        restoration = {"Command": ["python3", "restore_planar_metric_mesh.py", "adapted.meshb",
                                   "recipe.json", str(restored), *self.GATES],
                       "Artifacts": {"restored-mesh": {"Path": str(restored), "SHA256": "0" * 64}}}
        required = root / f"required-tetrahedra-{layer}.txt"
        required.write_text("3\n7\n40\n41\n99\n" if layer else "3\n7\n40\n")
        return seed, restoration, recipe, census, required

    def test_required_region_is_bound_and_seed_gated(self):
        from mesh_stage_contract import validate_required_region
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed, restoration, recipe, census, required = self.fixture(root)
            self.assertEqual(validate_required_region(seed, restoration, recipe, census, required),
                             recipe["RequiredTetrahedra"])
            seed0, restoration0, recipe0, census0, required0 = self.fixture(root, layer=False)
            validate_required_region(seed0, restoration0, recipe0, census0, required0)

            def rejected(message, seed_report=seed, restoration_report=restoration,
                         recipe_data=recipe, census_data=census, path=required):
                with self.assertRaisesRegex(ValueError, message):
                    validate_required_region(seed_report, restoration_report, recipe_data,
                                             census_data, path)

            def with_record(**changes):
                data = copy.deepcopy(recipe); data["RequiredTetrahedra"].update(changes); return data

            def with_quality(**changes):
                data = copy.deepcopy(census); data["SeedQualityOptimization"].update(changes); return data

            def with_value(report, option, value):
                report = copy.deepcopy(report)
                report["Command"][report["Command"].index(option) + 1] = value; return report

            def listing(name, text):
                path = root / name; path.write_text(text); return path

            rejected("lacks the required-tetrahedra record",
                     recipe_data={k: v for k, v in recipe.items() if k != "RequiredTetrahedra"})
            rejected("differs from the recipe record", recipe_data=with_record(Count=4))
            rejected("sorted, unique and in range", path=listing("dup.txt", "3\n3\n7\n40\n99\n"))
            rejected("sorted, unique and in range", path=listing("range.txt", "0\n3\n7\n40\n99\n"))
            rejected("sorted, unique and in range", path=listing("over.txt", "3\n7\n40\n99\n101\n"))
            rejected("non-empty positive integers", path=listing("empty.txt", ""))
            rejected("non-empty positive integers", path=listing("text.txt", "3\nx\n"))
            rejected("radius differs", recipe_data=with_record(CornerRadius=.2))
            rejected("differ from the recipe semantic corners",
                     recipe_data=with_record(PerCorner=[{"Point": [0., 0., 0.], "Tetrahedra": 2}]))
            rejected("differ from the recipe semantic corners",
                     recipe_data=with_record(PerCorner=[{"Point": [0., 0., 0.], "Tetrahedra": 2},
                                                        {"Point": [10., 0., 0.], "Tetrahedra": 0}]))
            rejected("differs from the recipe edge layer", recipe_data=with_record(LayerRequiredReach=.032))
            recorded_reach = copy.deepcopy(recipe); recorded_reach["EdgeLayer"]["RequiredReach"] = .032
            rejected("differs from the recipe edge layer", recipe_data=recorded_reach)
            rejected("differs from the recipe edge layer", recipe_data=with_record(PerSpan=[]))
            rejected("differs from the recipe edge layer",
                     recipe_data=with_record(PerSpan=[{"Span": [0., 0., 0., 1., 0., 0.], "Tetrahedra": 3}]))
            rejected("do not cover the list", recipe_data=with_record(Count=7),
                     path=listing("seven.txt", "1\n2\n3\n4\n5\n6\n7\n"))
            without_layer = {k: v for k, v in recipe.items() if k != "EdgeLayer"}
            rejected("without a recipe edge layer", recipe_data=without_layer)
            rejected("differs from the label-restoration command",
                     seed_report=with_value(seed, "--maximum-corner-aspect", "5"))
            rejected("differs from the label-restoration command",
                     restoration_report=with_value(restoration, "--minimum-scaled-jacobian", ".02"))
            rejected("must provide --maximum-quality-displacement-over-normal",
                     seed_report={"Command": seed["Command"][:-2]})
            rejected("does not record a gated", census_data={})
            rejected("does not record a gated", census_data=with_quality(MaximumCornerAspect=3.8))
            rejected("does not record a gated", census_data=with_quality(RequiredCellsBelowGateAfter=1))
            rejected("does not record a gated", census_data=with_quality(CornerAspectsAfter=[3.4, 4.1]))
            rejected("does not record a gated", census_data=with_quality(CornerAspectsAfter=[3.4]))
            rejected("below the scaled-Jacobian gate",
                     census_data=with_quality(RequiredMinimumScaledJacobianAfter=.009))
            # The Jacobian condition gate (decision 34): the seed carries the restorer's
            # bound, records it, gates every scaled-Jacobian-gated required cell by it,
            # and the restorer's required cells are within it.
            rejected("differs from the label-restoration command",
                     seed_report=with_value(seed, "--maximum-jacobian-condition", "500"))
            rejected("must provide --maximum-jacobian-condition",
                     seed_report={"Command": [token for token in seed["Command"]
                                              if token not in ("--maximum-jacobian-condition", "1000")]},
                     restoration_report={"Command": [token for token in restoration["Command"]
                                                     if token not in ("--maximum-jacobian-condition", "1000")],
                                         "Artifacts": restoration["Artifacts"]})
            rejected("does not record a gated", census_data=with_quality(MaximumJacobianCondition=900.))
            rejected("does not record a gated", census_data=with_quality(RequiredCellsAboveConditionAfter=1))
            rejected("above the Jacobian condition gate",
                     census_data=with_quality(RequiredMaximumJacobianConditionAfter=1000.5))
            (root / "restored-above.projection.json").write_text(json.dumps(
                {"RequiredTetrahedra": 5, "RequiredVertices": 12,
                 "RequiredMaximumJacobianCondition": 1200.}))
            above = copy.deepcopy(restoration)
            above["Artifacts"]["restored-mesh"]["Path"] = str(root / "restored-above.msh")
            rejected("required tetrahedra above the Jacobian condition gate", restoration_report=above)
            # The seed must have gated exactly the listed set (the pre-move count is
            # only reported) and the restorer must have found exactly that many
            # required cells in the adapted mesh.
            rejected("does not record a gated", census_data=with_quality(RequiredTetrahedra=6))
            rejected("does not record a gated", census_data=with_quality(RequiredTetrahedra=4))
            stale = copy.deepcopy(restoration)
            stale["Artifacts"]["restored-mesh"]["Path"] = str(root / "stale.msh")
            (root / "stale.projection.json").write_text(json.dumps({"RequiredTetrahedra": 4}))
            rejected("different number of required tetrahedra", restoration_report=stale)
            (root / "stale.projection.json").write_text(json.dumps({"RequiredVertices": 4}))
            rejected("Restored required tetrahedra", restoration_report=stale)

    def test_edge_layer_quality_rule_is_bound_between_seed_and_restorer(self):
        from edge_volume_metric import EDGE_LAYER_QUALITY_RULE
        from mesh_stage_contract import validate_required_region
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seed, restoration, recipe, census, required = self.fixture(root)
            option = ["--edge-layer-maximum-aspect", "100"]
            seed_rule = {"MaximumEdgeAspect": 100., "ScaledJacobianRoundoffFloor": 1e-12,
                         "LayerCells": 3, "CellsAboveBoundAfter": 0, "CellsBelowRoundoffFloorAfter": 0,
                         "MaximumEdgeAspectAfter": 71.2, "MinimumScaledJacobian": 4e-4}
            restorer_rule = {"Passes": True, "MaximumEdgeAspectBound": 100.,
                             "ScaledJacobianRoundoffFloor": 1e-12, "MaximumEdgeAspect": 71.2,
                             "PositiveOrientation": True, "Cells": 3, "Rule": EDGE_LAYER_QUALITY_RULE}
            def ruled(seed_option=option, restorer_option=option, seed_record=seed_rule,
                      restorer_record=restorer_rule, recipe_data=recipe):
                seed_report = copy.deepcopy(seed); seed_report["Command"] += seed_option
                restoration_report = copy.deepcopy(restoration)
                restoration_report["Command"] += restorer_option
                restored = root / "ruled.msh"
                restoration_report["Artifacts"]["restored-mesh"]["Path"] = str(restored)
                report = {"RequiredTetrahedra": 5, "EdgeLayerQuality": restorer_record}
                restored.with_suffix(".projection.json").write_text(json.dumps(report))
                census_data = copy.deepcopy(census)
                census_data["SeedQualityOptimization"]["EdgeLayerQualityRule"] = seed_record
                return validate_required_region(seed_report, restoration_report, recipe_data,
                                                census_data, required)
            self.assertEqual(ruled(), recipe["RequiredTetrahedra"])
            def rejected(message, **changes):
                with self.assertRaisesRegex(ValueError, message):
                    ruled(**changes)
            rejected("differs from the label-restoration command", restorer_option=[])
            rejected("differs from the label-restoration command",
                     restorer_option=["--edge-layer-maximum-aspect", "90"])
            rejected("recorded without the bound option", seed_option=[], restorer_option=[])
            rejected("does not record a gated edge-layer", seed_record=None)
            rejected("does not record a gated edge-layer", seed_record={**seed_rule, "MaximumEdgeAspect": 90.})
            rejected("does not record a gated edge-layer", seed_record={**seed_rule, "CellsAboveBoundAfter": 1})
            rejected("does not record a gated edge-layer",
                     seed_record={**seed_rule, "CellsBelowRoundoffFloorAfter": 1})
            rejected("does not record a gated edge-layer", seed_record={**seed_rule, "LayerCells": 0})
            rejected("does not record a passing edge-layer", restorer_record=None)
            rejected("does not record a passing edge-layer", restorer_record={**restorer_rule, "Passes": False})
            rejected("does not record a passing edge-layer",
                     restorer_record={**restorer_rule, "MaximumEdgeAspect": 100.5})
            rejected("does not record a passing edge-layer",
                     restorer_record={**restorer_rule, "MaximumEdgeAspectBound": 90.})
            rejected("finite bound above 1", seed_option=["--edge-layer-maximum-aspect", "1"],
                     restorer_option=["--edge-layer-maximum-aspect", "1"])
            without_layer = {k: v for k, v in recipe.items() if k != "EdgeLayer"}
            without_layer["RequiredTetrahedra"] = {**recipe["RequiredTetrahedra"],
                                                   "LayerRequiredReach": None, "PerSpan": [],
                                                   "PerCorner": [{"Point": [0., 0., 0.], "Tetrahedra": 3},
                                                                 {"Point": [10., 0., 0.], "Tetrahedra": 2}]}
            rejected("without a recipe edge layer", recipe_data=without_layer)


class EdgeLayerContractTest(unittest.TestCase):
    """The seed, metric and adapter commands honor one edge layer or none."""

    @staticmethod
    def layer_fixture():
        offsets = [.001, .003, .007, .015, .031]
        layer = {"EdgeSize": .001, "GrowthRatio": 2.0, "Aspect": 4.0, "Reach": .024, "Layers": 5,
                 "RowOffsets": offsets, "LayerThickness": .031, "Spans": [[1, 0, .1, 9, 0, .1]],
                 "SpanCount": 1, "TotalSpanLength": 8.0, "SeedRows": 10, "SeedRowNodes": 200}
        recipe = {"NormalSize": .025, "SurfaceProtectionRadius": .05, "EdgeLayer": layer}
        census = {"EdgeLayer": {"EdgeSize": .001, "GrowthRatio": 2.0, "Aspect": 4.0,
                                "NormalSize": .025, "Layers": 5, "RowOffsets": offsets,
                                "LayerThickness": .031, "TotalSpanLength": 8.0, "Rows": 10,
                                "RowNodes": 200, "Curves": [{"Curve": 7}]}}
        seed = {"Command": ["julia", "mesh_spatial_coupon.jl", "sig.csv", "fabricated", "seed.msh",
                            "--lc-fine", ".025", "--edge-size", ".001", "--edge-growth-ratio", "2",
                            "--edge-layer-aspect", "4"]}
        metric = {"Command": ["python3", "prepare_edge_metric_scout.py", "seed.msh", "metric",
                              "--normal", ".025", "--edge-size", ".001", "--edge-growth-ratio", "2.0",
                              "--edge-layer-aspect", "4.0"]}
        adaptation = {"Command": ["python3", "run_native_mmg_adaptation.py", "--hmin", ".001",
                                  "--hgrad", "1.15"]}
        return seed, metric, adaptation, recipe, census

    def test_edge_layer_is_bound_across_seed_metric_and_adapter_or_absent(self):
        from mesh_stage_contract import validate_edge_layer
        seed, metric, adaptation, recipe, census = self.layer_fixture()
        self.assertEqual(validate_edge_layer(seed, metric, adaptation, recipe, census),
                         recipe["EdgeLayer"])
        # The ratio and aspect may be left at their documented defaults.
        defaults = copy.deepcopy(seed); del defaults["Command"][-4:]
        default_metric = copy.deepcopy(metric); del default_metric["Command"][-4:]
        validate_edge_layer(defaults, default_metric, adaptation, recipe, census)

        def rejected(message, seed_report=seed, metric_report=metric, adaptation_report=adaptation,
                     recipe_data=recipe, census_data=census):
            with self.assertRaisesRegex(ValueError, message):
                validate_edge_layer(seed_report, metric_report, adaptation_report, recipe_data,
                                    census_data)

        def with_value(report, option, value):
            report = copy.deepcopy(report)
            report["Command"][report["Command"].index(option) + 1] = value; return report

        def with_layer(record, **changes):
            record = copy.deepcopy(record); record["EdgeLayer"].update(changes); return record

        rejected("Seed command --edge-size", seed_report=with_value(seed, "--edge-size", ".002"))
        rejected("Seed command --edge-growth-ratio", seed_report=with_value(seed, "--edge-growth-ratio", "1.5"))
        rejected("Seed command --edge-layer-aspect", seed_report=with_value(seed, "--edge-layer-aspect", "8"))
        rejected("Metric command --edge-size", metric_report=with_value(metric, "--edge-size", ".0005"))
        rejected("Metric command --edge-layer-aspect", metric_report=with_value(metric, "--edge-layer-aspect", "5"))
        rejected("--hmin differs", adaptation_report=with_value(adaptation, "--hmin", ".025"))
        rejected("census edge layer EdgeSize", census_data=with_layer(census, EdgeSize=.002))
        rejected("census edge layer Rows", census_data=with_layer(census, Rows=9))
        rejected("census edge layer rows differ", census_data=with_layer(census, RowNodes=199))
        rejected("do not follow", recipe_data=with_layer(recipe, RowOffsets=[.001, .003, .007, .015, .03],
                                                         LayerThickness=.03),
                 census_data=with_layer(census, RowOffsets=[.001, .003, .007, .015, .03], LayerThickness=.03))
        rejected("do not follow", recipe_data=with_layer(recipe, Reach=.02))
        rejected("thicker than the frozen band", recipe_data=dict(recipe, SurfaceProtectionRadius=.03))
        rejected("spans differ", recipe_data=with_layer(recipe, SpanCount=2))
        rejected("sizes are invalid", recipe_data=with_layer(recipe, EdgeSize=.03),
                 census_data=with_layer(census, EdgeSize=.03))
        rejected("record is missing", census_data={})
        # Without a recipe layer: no seed layer, no metric layer options, hmin = NormalSize.
        plain_recipe = {k: v for k, v in recipe.items() if k != "EdgeLayer"}
        plain_seed = copy.deepcopy(seed); del plain_seed["Command"][-6:]
        plain_metric = copy.deepcopy(metric); del plain_metric["Command"][-6:]
        plain_adaptation = with_value(adaptation, "--hmin", ".025")
        self.assertIsNone(validate_edge_layer(plain_seed, plain_metric, plain_adaptation, plain_recipe, {}))
        validate_edge_layer(plain_seed, plain_metric, plain_adaptation, plain_recipe,
                            {"EdgeLayer": {"EdgeSize": 0.0}})
        rejected("Seed command seeds an edge layer", seed_report=seed, metric_report=plain_metric,
                 adaptation_report=plain_adaptation, recipe_data=plain_recipe, census_data={})
        rejected("Metric command binds an edge layer", seed_report=plain_seed, metric_report=metric,
                 adaptation_report=plain_adaptation, recipe_data=plain_recipe, census_data={})
        rejected("census records an edge layer", seed_report=plain_seed, metric_report=plain_metric,
                 adaptation_report=plain_adaptation, recipe_data=plain_recipe)
        rejected("--hmin differs", seed_report=plain_seed, metric_report=plain_metric,
                 adaptation_report=adaptation, recipe_data=plain_recipe, census_data={})


class CornerGradingContractTest(unittest.TestCase):
    """Decision 33: the seed and the metric stage prescribe one corner grading
    (CornerSize growing by the edge-layer ratio to NormalSize inside the corner
    balls) or none; with a layer the rows reach the ball boundary."""

    @staticmethod
    def fixture():
        grading = {"CornerSize": .004, "GrowthRatio": 2.0, "NormalSize": .025, "Radius": .1,
                   "Reach": .021, "ShellRadii": [.004, .012, .028, .1], "ShellSizes": [.004, .008, .016, .025]}
        recipe = {"NormalSize": .025, "CornerIsotropyRadius": .1, "CornerGrading": dict(grading),
                  "TruePhysicalCorners": [[0., 0., 0.], [10., 0., 0.]]}
        census = {"CornerGrading": dict(grading),
                  "EdgeLayer": {"EdgeSize": .004, "LayerReachesCornerBall": True, "TaperSubdivisions": [],
                                "UnlayeredEdgeLengthPerCorner": [
                                    {"Corner": 0, "UnlayeredLengths": [.1, .1], "Maximum": .1},
                                    {"Corner": 1, "UnlayeredLengths": [.1], "Maximum": .1}]}}
        seed = {"Command": ["julia", "mesh_spatial_coupon.jl", "sig.csv", "fabricated", "seed.msh",
                            "--lc-fine", ".025", "--edge-size", ".004", "--edge-growth-ratio", "2",
                            "--corner-size", ".004"]}
        metric = {"Command": ["python3", "prepare_edge_metric_scout.py", "seed.msh", "metric",
                              "--normal", ".025", "--edge-size", ".004", "--corner-size", "0.004"]}
        return seed, metric, recipe, census

    def test_corner_grading_is_bound_across_seed_metric_census_and_recipe_or_absent(self):
        from mesh_stage_contract import validate_corner_grading
        seed, metric, recipe, census = self.fixture()
        self.assertEqual(validate_corner_grading(seed, metric, recipe, census), recipe["CornerGrading"])

        def rejected(message, seed_report=seed, metric_report=metric, recipe_data=recipe,
                     census_data=census):
            with self.assertRaisesRegex(ValueError, message):
                validate_corner_grading(seed_report, metric_report, recipe_data, census_data)

        def with_value(report, option, value):
            report = copy.deepcopy(report)
            report["Command"][report["Command"].index(option) + 1] = value; return report

        def with_grading(record, **changes):
            record = copy.deepcopy(record); record["CornerGrading"].update(changes); return record

        rejected("differ from the recipe CornerSize", seed_report=with_value(seed, "--corner-size", ".002"))
        rejected("differ from the recipe CornerSize", metric_report=with_value(metric, "--corner-size", ".008"))
        rejected("differ from the recipe CornerSize", census_data=with_grading(census, CornerSize=.002))
        rejected("--edge-growth-ratio differs", seed_report=with_value(seed, "--edge-growth-ratio", "1.5"))
        rejected("do not follow", recipe_data=with_grading(recipe, ShellRadii=[.004, .012, .028]))
        rejected("do not follow", recipe_data=with_grading(recipe, Reach=.02))
        # The grading must reach NormalSize inside the ball (radius 0.02 < reach 0.021).
        rejected("do not follow", recipe_data={**with_grading(recipe, Radius=.02), "CornerIsotropyRadius": .02},
                 census_data=with_grading(census, Radius=.02))
        rejected("census corner grading differs", census_data=with_grading(census, ShellRadii=[.004, .1]))
        rejected("do not follow", recipe_data=with_grading(recipe, ShellSizes=[.004, .008, .016, .02]))
        rejected("sizes are invalid", recipe_data=with_grading(recipe, CornerSize=.03),
                 census_data=with_grading(census, CornerSize=.03),
                 seed_report=with_value(seed, "--corner-size", ".03"),
                 metric_report=with_value(metric, "--corner-size", ".03"))
        # With a layer the rows must reach the ball (no taper, un-layered lengths recorded).
        tapered = copy.deepcopy(census); tapered["EdgeLayer"]["LayerReachesCornerBall"] = False
        rejected("reaching the corner balls", census_data=tapered)
        tapered = copy.deepcopy(census); tapered["EdgeLayer"]["TaperSubdivisions"] = [2]
        rejected("reaching the corner balls", census_data=tapered)
        tapered = copy.deepcopy(census); del tapered["EdgeLayer"]["UnlayeredEdgeLengthPerCorner"][1]
        rejected("reaching the corner balls", census_data=tapered)
        # Without a layer the grading stands alone.
        self.assertIsNotNone(validate_corner_grading(seed, metric, recipe,
                                                     {"CornerGrading": census["CornerGrading"]}))
        # Without a recipe record: no option anywhere, census none (0 or absent).
        plain_recipe = {k: v for k, v in recipe.items() if k != "CornerGrading"}
        plain_seed = copy.deepcopy(seed); del plain_seed["Command"][-2:]
        plain_metric = copy.deepcopy(metric); del plain_metric["Command"][-2:]
        self.assertIsNone(validate_corner_grading(plain_seed, plain_metric, plain_recipe, {}))
        self.assertIsNone(validate_corner_grading(plain_seed, plain_metric, plain_recipe,
                                                  {"CornerGrading": {"CornerSize": 0.0}}))
        rejected("without a recipe CornerGrading", seed_report=seed, metric_report=plain_metric,
                 recipe_data=plain_recipe, census_data={})
        rejected("without a recipe CornerGrading", seed_report=plain_seed, metric_report=metric,
                 recipe_data=plain_recipe, census_data={})
        rejected("without a recipe CornerGrading", seed_report=plain_seed, metric_report=plain_metric,
                 recipe_data=plain_recipe)


if __name__ == "__main__":
    unittest.main()
