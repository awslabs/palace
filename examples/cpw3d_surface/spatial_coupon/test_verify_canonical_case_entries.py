# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from general_mesh_manifest import sha256
# Imported as a module so unittest does not re-collect the fixture suite's tests here.
import test_general_mesh_manifest as fixture_suite
from verify_canonical_case_entries import (covariance_failures, validate_calibration_commands,
                                           verify_case)

HERE = Path(__file__).resolve().parent


class VerifyCanonicalCaseEntriesTest(unittest.TestCase):
    """The per-entry verifier judges every variant and the covariance comparison
    even when one variant fails, and reports every failure together."""

    def fixture(self, root):
        harness = fixture_suite.GeneralMeshManifestTest(
            "test_staged_gmsh_evidence_accepts_complete_fixture_matrix")
        manifest_path, manifest = harness.make_suite(root)
        audits = harness.produce_matrix(root, manifest_path, manifest)
        return manifest_path, manifest, audits

    @staticmethod
    def _break_evidence(audits, case_id, variant_id):
        """Edit a recorded measurement of one variant's normalized evidence so the
        bound producer records no longer agree with it (a binding failure)."""
        path = audits / f"{case_id}--{variant_id}.json"
        evidence = json.loads(path.read_text())
        evidence["MeshQuality"]["MinimumScaledJacobian"] = 0.0
        path.write_text(json.dumps(evidence))

    def test_passing_case_reports_every_variant_and_the_covariance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest, audits = self.fixture(root)
            report = verify_case(manifest_path, audits, "base")
            self.assertTrue(report["Passed"]); self.assertEqual(report["Failures"], [])
            self.assertEqual(set(report["Entries"]), {"identity", "rotate-z-0.63"})
            for entry in report["Entries"].values():
                self.assertTrue(entry["Passed"]); self.assertEqual(entry["GateFailures"], [])
                self.assertIsNone(entry["Error"])
                self.assertEqual(entry["CanonicalBuildId"], report["SharedCanonicalBuildId"])
                self.assertIn("MaximumRelativeMeasureError", entry["ProtectedSurfaces"])
                self.assertNotIn("LongShortEdgeComponents", entry["TraceDiagonal"])
            self.assertEqual(report["TransformComparisonFailures"], [])
            self.assertEqual(report["CanonicalReuseFailures"], [])
            self.assertEqual(report["CanonicalBuildIds"], [report["SharedCanonicalBuildId"]])
            self.assertEqual(report["Manifest"]["SHA256"], sha256(manifest_path))
            # The CLI writes the same report and exits 0.
            output = root / "verification.json"
            result = subprocess.run(
                [sys.executable, str(HERE / "verify_canonical_case_entries.py"),
                 str(manifest_path), str(audits), str(output), "base"],
                capture_output=True, text=True, check=False, cwd=HERE)
            self.assertEqual(result.returncode, 0, result.stderr)
            written = json.loads(output.read_text())
            self.assertTrue(written["Passed"])
            self.assertEqual(written["Entries"].keys(), report["Entries"].keys())
            with self.assertRaises(subprocess.CalledProcessError):
                subprocess.run(
                    [sys.executable, str(HERE / "verify_canonical_case_entries.py"),
                     str(manifest_path), str(audits), str(output), "base"],
                    capture_output=True, text=True, check=True, cwd=HERE)

    def test_failing_identity_does_not_stop_the_other_variant_or_the_covariance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest, audits = self.fixture(root)
            self._break_evidence(audits, "base", "identity")
            report = verify_case(manifest_path, audits, "base")
            self.assertFalse(report["Passed"])
            identity, rotated = report["Entries"]["identity"], report["Entries"]["rotate-z-0.63"]
            self.assertFalse(identity["Passed"]); self.assertIsNotNone(identity["Error"])
            self.assertIn("bound producer records", identity["Error"])
            # The second variant was still judged completely ...
            self.assertTrue(rotated["Passed"]); self.assertIsNone(rotated["Error"])
            self.assertEqual(rotated["GateFailures"], [])
            self.assertIn("CanonicalBuildId", rotated)
            # ... and the covariance comparison was still evaluated on both evidences.
            self.assertEqual(report["TransformComparisonFailures"], [])
            self.assertIsNotNone(report["TransformMaximumCoordinateError"])
            self.assertEqual([f for f in report["Failures"] if f.startswith("identity:")],
                             [f"identity: {identity['Error']}"])
            # The CLI still writes the full report and exits nonzero.
            output = root / "verification.json"
            result = subprocess.run(
                [sys.executable, str(HERE / "verify_canonical_case_entries.py"),
                 str(manifest_path), str(audits), str(output), "base"],
                capture_output=True, text=True, check=False, cwd=HERE)
            self.assertEqual(result.returncode, 1, result.stderr)
            written = json.loads(output.read_text())
            self.assertFalse(written["Passed"])
            self.assertEqual(set(written["Entries"]), {"identity", "rotate-z-0.63"})

    def test_gate_failures_of_every_variant_and_covariance_failures_are_all_reported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest, audits = self.fixture(root)
            # A gate both variants miss (bindings intact): each entry lists it and the
            # covariance is still compared.
            for variant_id in ("identity", "rotate-z-0.63"):
                path = audits / f"base--{variant_id}.json"
                evidence = json.loads(path.read_text()); evidence["Version"] = 2
                path.write_text(json.dumps(evidence))
            report = verify_case(manifest_path, audits, "base")
            self.assertFalse(report["Passed"])
            for variant_id, entry in report["Entries"].items():
                self.assertEqual(entry["GateFailures"], ["provenance-binding"], variant_id)
                self.assertIsNone(entry["Error"], variant_id)
                self.assertIsNotNone(entry.get("CanonicalBuildId"), variant_id)
            self.assertEqual(
                sorted(report["Failures"]),
                sorted(f"{variant}: gate failure: provenance-binding"
                       for variant in ("identity", "rotate-z-0.63")))
            self.assertEqual(report["TransformComparisonFailures"], [])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest, audits = self.fixture(root)
            # A transformed variant bound to another seed: its entry errors AND the
            # covariance comparison reports the broken source-seed covariance.
            path = audits / "base--rotate-z-0.63.json"
            evidence = json.loads(path.read_text())
            evidence["IdentitySeedMeshSHA256"] = "0" * 64
            path.write_text(json.dumps(evidence))
            report = verify_case(manifest_path, audits, "base")
            self.assertFalse(report["Passed"])
            self.assertTrue(report["Entries"]["identity"]["Passed"])
            self.assertIsNotNone(report["Entries"]["rotate-z-0.63"]["Error"])
            self.assertEqual(report["TransformComparisonFailures"],
                             ["exact source-seed covariance"])
            self.assertIn("covariance: exact source-seed covariance", report["Failures"])
            # Missing evidence for the transformed variant: judged as missing, the
            # identity entry is still complete.
            path.unlink()
            report = verify_case(manifest_path, audits, "base")
            self.assertFalse(report["Passed"])
            self.assertTrue(report["Entries"]["identity"]["Passed"])
            self.assertIsNotNone(report["Entries"]["rotate-z-0.63"]["Error"])
            self.assertEqual(report["TransformComparisonFailures"],
                             ["missing transform evidence"])
            self.assertIn("covariance: missing transform evidence", report["Failures"])
            with self.assertRaisesRegex(ValueError, "not in the manifest"):
                verify_case(manifest_path, audits, "no-such-case")

    def test_covariance_comparison_reports_seed_and_physical_failures(self):
        manifest = {"Gates": {"CornerTolerance": 1e-8}}
        case = {"TransformComparison": {"Reference": "identity", "Transformed": "rotated",
                                        "MaximumRelativeVolumeError": 1e-8,
                                        "MaximumRelativeSurfaceMeasureError": 1e-8,
                                        "MaximumProtectedSupportHausdorff": 1e-8,
                                        "MaximumProtectedMeasureError": 1e-8,
                                        "MaximumQualityDistributionRelativeError": 1e-8,
                                        "MaximumAnisotropyRelativeError": 1e-8,
                                        "MaximumComplexityRatio": 1.01}}
        physical = {"LabelsMaterialsAdjacencyMatch": True,
                    "ReferenceInvariants": {"Volume:1": 1.0, "Area:1": 2.0},
                    "TransformedInvariants": {"Volume:1": 1.0, "Area:1": 2.0},
                    "ProtectedSurfaces": {"PlaneSupportsMatch": True, "TopologyMatches": True,
                                          "MaximumSupportVertexDistance": 0.0,
                                          "MaximumRelativeMeasureError": 0.0},
                    "ReferenceQuality": {"PositiveOrientation": True,
                                         "ScaledJacobianQuantiles": [.5],
                                         "JacobianConditionQuantiles": [2.]},
                    "TransformedQuality": {"PositiveOrientation": True,
                                           "ScaledJacobianQuantiles": [.5],
                                           "JacobianConditionQuantiles": [2.]}}
        anisotropy = {"TangentialP50": .1, "Transverse1P90": .02, "Transverse2P90": .02}
        identity = {"Mesh": {"SHA256": "a" * 64}, "IdentityMeshSHA256": "a" * 64,
                    "IdentitySeedMeshSHA256": "b" * 64, "AchievedAnisotropy": anisotropy,
                    "Complexity": {"H1DOFs": 100}, "Resources": {"Elements": 50}}
        rotated = {"Mesh": {"SHA256": "c" * 64}, "IdentityMeshSHA256": "a" * 64,
                   "IdentitySeedMeshSHA256": "b" * 64, "TransformMaximumCoordinateError": 0.0,
                   "PhysicalCovariance": physical, "AchievedAnisotropy": anisotropy,
                   "Complexity": {"H1DOFs": 100}, "Resources": {"Elements": 50}}
        failures, error = covariance_failures(manifest, case,
                                              {"identity": identity, "rotated": rotated})
        self.assertEqual((failures, error), ([], 0.0))
        other_seed = dict(rotated, IdentitySeedMeshSHA256="d" * 64)
        failures, _ = covariance_failures(manifest, case,
                                          {"identity": identity, "rotated": other_seed})
        self.assertEqual(failures, ["exact source-seed covariance"])
        heavier = dict(rotated, Complexity={"H1DOFs": 200})
        failures, _ = covariance_failures(manifest, case,
                                          {"identity": identity, "rotated": heavier})
        self.assertEqual(failures, ["complexity ratio"])
        failures, error = covariance_failures(manifest, case, {"identity": identity})
        self.assertEqual((failures, error), (["missing transform evidence"], None))

    @staticmethod
    def _declare_calibration(manifest_path, case_id, seed, metric, production, adaptation=None):
        """Label `case_id` as a calibration case with the given declared option values."""
        manifest = json.loads(manifest_path.read_text())
        # A case may carry a Calibration block only in a labeled calibration manifest.
        manifest.setdefault("Calibration", {"Purpose": "fixture calibration manifest"})
        case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
        case["Calibration"] = {"Label": "fixture calibration", "BaseCase": case_id,
                               "SeedCommandOptions": seed, "MetricCommandOptions": metric,
                               "ProductionValues": production}
        if adaptation is not None:
            case["Calibration"]["AdaptationCommandOptions"] = adaptation
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    def test_calibration_option_declarations_are_bound_to_the_recorded_commands(self):
        # The fixture seed executes --lc-fine NORMAL_SIZE, the metric --normal
        # NORMAL_SIZE / --tangent CORNER_ISOTROPY_RADIUS and the adaptation --hmin .1; a
        # calibration label is accepted only when the recorded commands executed exactly
        # what it declares.
        normal, radius = fixture_suite.NORMAL_SIZE, fixture_suite.CORNER_ISOTROPY_RADIUS
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest_path, manifest, audits = self.fixture(root)
            self._declare_calibration(manifest_path, "base", {"--lc-fine": normal},
                                      {"--normal": normal},
                                      {"--lc-fine": 2 * normal, "--normal": 2 * normal,
                                       "--hmin": 0.2}, adaptation={"--hmin": 0.1})
            report = verify_case(manifest_path, audits, "base")
            self.assertTrue(report["Passed"], report["Failures"])
            rejected = (
                # Labeled as a finer variant than the one that was built (a "V2" label
                # on a root built with the "V1" options): the declared pair is absent
                # and the recorded value is the production value.
                ({"--lc-fine": normal / 2}, {}, {"--lc-fine": normal}, None,
                 "does not execute calibration option --lc-fine"),
                # A declared metric option the recorded command never executed.
                ({}, {"--far-growth": 0.5}, {"--far-growth": 1.0}, None,
                 "does not execute calibration option --far-growth"),
                # The adaptation hmin labeled halved while the root kept it.
                ({}, {}, {"--hmin": 0.1}, {"--hmin": 0.05},
                 "does not execute calibration option --hmin"),
                # Declared at the production value: not a calibration variant.
                ({"--lc-fine": normal}, {}, {"--lc-fine": normal}, None,
                 "declared at its production value"),
                # An undeclared production option executed away from production.
                ({"--lc-fine": normal}, {}, {"--lc-fine": 2 * normal, "--tangent": 2 * radius},
                 None, "undeclared calibration option --tangent"),
                # A declared option without a production value.
                ({"--lc-fine": normal}, {}, {}, None, "has no production value"))
            for seed, metric, production, adaptation, message in rejected:
                self._declare_calibration(manifest_path, "base", seed, metric, production,
                                          adaptation=adaptation)
                report = verify_case(manifest_path, audits, "base")
                self.assertFalse(report["Passed"], (seed, metric, production, adaptation))
                for entry in report["Entries"].values():
                    self.assertIn(message, entry["Error"] or "",
                                  (seed, metric, production, adaptation))
            # Without a Calibration block nothing is required.
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
            self.assertTrue(verify_case(manifest_path, audits, "base")["Passed"])

    def test_validate_calibration_commands_reads_the_recorded_argv_numerically(self):
        # Recorded argv tokens are strings (".05" and "0.05" are the same value); a
        # repeated option is not "exactly once"; option=value tokens are not executed.
        stages = {"seed-generation": {"Command": ["julia", "mesher.jl", "--lc-tangent", ".05",
                                                  "--lc-fine", ".025"]},
                  "metric-preparation": {"Command": ["python3", "metric.py", "--normal",
                                                     ".025", "--far-growth", "0.5"]},
                  "native-adaptation-mmg": {"Command": ["python3", "adapt.py", "--hmin",
                                                        ".025"]},
                  "label-restoration": {"Command": ["python3", "restore.py",
                                                    "--minimum-scaled-jacobian", ".01"]}}
        case = {"Calibration": {"SeedCommandOptions": {"--lc-tangent": 0.05},
                                "MetricCommandOptions": {"--far-growth": 0.5},
                                "ProductionValues": {"--lc-tangent": 0.1, "--far-growth": 1.0}}}
        validate_calibration_commands(case, stages)
        validate_calibration_commands({}, stages)
        # An adaptation block binds --hmin; without the block hmin is not declared.
        halved = {"Calibration": {**case["Calibration"],
                                  "AdaptationCommandOptions": {"--hmin": 0.0125},
                                  "ProductionValues": {**case["Calibration"]["ProductionValues"],
                                                       "--hmin": 0.025}}}
        with self.assertRaisesRegex(ValueError, "--hmin=0.0125 exactly once"):
            validate_calibration_commands(halved, stages)
        validate_calibration_commands(halved, {
            **stages, "native-adaptation-mmg": {"Command": ["python3", "adapt.py", "--hmin",
                                                            "0.0125"]}})
        with self.assertRaisesRegex(ValueError, "exactly once"):
            validate_calibration_commands(case, {
                **stages, "seed-generation": {"Command": stages["seed-generation"]["Command"]
                                              + ["--lc-tangent", ".1"]}})
        with self.assertRaisesRegex(ValueError, "exactly once"):
            validate_calibration_commands(case, {
                **stages, "metric-preparation": {"Command": ["python3", "metric.py",
                                                             "--far-growth=0.5"]}})
        # An option shared by two stage commands is declared for each with one value;
        # every declaring stage executes it and a differing value is rejected.
        shared = {"Calibration": {"SeedCommandOptions": {"--edge-size": 0.004},
                                  "MetricCommandOptions": {"--edge-size": 0.004},
                                  "ProductionValues": {"--edge-size": 0.0}}}
        layered = {**stages,
                   "seed-generation": {"Command": ["julia", "seed.jl", "--edge-size", ".004"]},
                   "metric-preparation": {"Command": ["python3", "metric.py", "--edge-size", "0.004"]}}
        validate_calibration_commands(shared, layered)
        with self.assertRaisesRegex(ValueError, "declared with two values"):
            validate_calibration_commands({"Calibration": {
                **shared["Calibration"], "MetricCommandOptions": {"--edge-size": 0.002}}}, layered)
        with self.assertRaisesRegex(ValueError, "exactly once"):
            validate_calibration_commands(shared, {
                **layered, "metric-preparation": {"Command": ["python3", "metric.py"]}})
        with self.assertRaisesRegex(ValueError, "away from its production value"):
            validate_calibration_commands({"Calibration": {
                **shared["Calibration"], "MetricCommandOptions": {}}}, layered)
        with self.assertRaisesRegex(ValueError, "ends with option"):
            validate_calibration_commands(case, {
                **stages, "metric-preparation": {"Command": ["python3", "--far-growth"]}})
        # The label-restoration command is bound the same way (RestorationCommandOptions).
        restored = {"Calibration": {"RestorationCommandOptions": {"--edge-layer-maximum-aspect": 100.},
                                    "SeedCommandOptions": {"--edge-layer-maximum-aspect": 100.},
                                    "ProductionValues": {"--edge-layer-maximum-aspect": 0.}}}
        ruled = {**stages,
                 "seed-generation": {"Command": ["julia", "seed.jl", "--edge-layer-maximum-aspect", "100"]},
                 "label-restoration": {"Command": ["python3", "restore.py", "--edge-layer-maximum-aspect", "100"]}}
        validate_calibration_commands(restored, ruled)
        with self.assertRaisesRegex(ValueError, "label-restoration command does not execute"):
            validate_calibration_commands(restored, {**ruled, "label-restoration": stages["label-restoration"]})

    def test_edge_layer_quality_rule_binds_seed_and_restorer_to_the_manifest_gate(self):
        from verify_canonical_case_entries import validate_edge_layer_quality_rule_binding
        rule = {"MaximumEdgeAspect": 100., "ScaledJacobianRoundoffFloor": 1e-12}
        calibration = {"Gates": {"EdgeLayerQualityRule": rule}}
        production = {"Gates": {}}
        def stages(seed=(), restorer=()):
            return {"seed-generation": {"Command": ["julia", "seed.jl", *seed]},
                    "label-restoration": {"Command": ["python3", "restore.py", *restorer]}}
        option = ("--edge-layer-maximum-aspect", "100")
        declaring = {"Calibration": {"EdgeLayerQualityRule": {"MaximumEdgeAspect": 100.}}}
        plain = {"Calibration": {}}
        # Declaring case: both stages execute the manifest bound exactly once.
        self.assertEqual(validate_edge_layer_quality_rule_binding(
            calibration, declaring, stages(option, option)), 100.)
        # Non-declaring cases (4 nm layer, production cases) execute it in neither stage.
        self.assertIsNone(validate_edge_layer_quality_rule_binding(calibration, plain, stages()))
        self.assertIsNone(validate_edge_layer_quality_rule_binding(production, {}, stages()))
        def rejected(message, manifest, case, bounded):
            with self.assertRaisesRegex(ValueError, message):
                validate_edge_layer_quality_rule_binding(manifest, case, bounded)
        rejected("exactly once", calibration, declaring, stages(option, ()))
        rejected("exactly once", calibration, declaring, stages((), option))
        rejected("exactly once", calibration, declaring,
                 stages(option, ("--edge-layer-maximum-aspect", "90")))
        rejected("exactly once", calibration, declaring, stages(option + option, option))
        rejected("exactly once", calibration,
                 {"Calibration": {"EdgeLayerQualityRule": {"MaximumEdgeAspect": 90.}}},
                 stages(option, option))
        rejected("without a declared edge-layer quality rule", calibration, plain, stages(option, option))
        rejected("without a declared edge-layer quality rule", production, {}, stages(option, option))
        rejected("manifest gates lack", production, declaring, stages(option, option))


if __name__ == "__main__":
    unittest.main()
