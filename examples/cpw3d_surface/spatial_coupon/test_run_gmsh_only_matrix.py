# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""run_gmsh_only_matrix: the library-build record is read from the per-case root's
records alone (schema, stop attribution, headroom flags, totals); with Julia, the
`coupon-library build` command end to end on two small fixtures - a registered
temporary copy of a gallery case (real two-pass probe) and one-edge-straight."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

from run_gmsh_only_matrix import (HEADROOM_FRACTION, LIBRARY_BUILD_RECORD, STATUS_BUILT, STATUS_FAILED,
                                  STATUS_UNSUPPORTED, case_record, library_totals, run_matrix)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PRODUCTION_MANIFEST = HERE / "geometry-independence-suite.json"
GIB = 2**30
# The manifest's bounded-stage limits (decision 61c: 3600 s / 12 GiB; the synthetic stage
# reports below carry them as run_bounded_mesher.py records them).
PRODUCTION_GATES = json.loads(PRODUCTION_MANIFEST.read_text())["Gates"]
STAGE_LIMITS = (float(PRODUCTION_GATES["MaximumSeconds"]), float(PRODUCTION_GATES["MaximumRSSGiB"]))


def stage_report(root, name, seconds, peak_gib, *, code=0, reason=None, limits=STAGE_LIMITS, audit=False):
    suffix = ".audit.log.json" if audit else ".log.json"
    (root / f"{name}{suffix}").write_text(json.dumps({
        "Version": 3, "Command": ["x"], "Seconds": seconds, "PeakProcessTreeRSSBytes": int(peak_gib * GIB),
        "ReturnCode": code, "StopReason": reason, "Limits": {"Seconds": limits[0], "MemoryGiB": limits[1]},
        "Artifacts": {}}))


class LibraryBuildRecordTest(unittest.TestCase):
    """The per-case record from synthetic root records (no mesh is built)."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.manifest_path = PRODUCTION_MANIFEST
        self.manifest = json.loads(self.manifest_path.read_text())
        self.case = next(case for case in self.manifest["Cases"] if "BasisContract" not in case["Source"]["Files"])
        self.gates = self.manifest["Gates"]

    def record(self, root, **overrides):
        options = dict(h1_order=4, driver_return_code=0, wall_seconds=12.5)
        options.update(overrides)
        return case_record(self.manifest, self.manifest_path, self.case, root, **options)

    def write_built_root(self, root, *, elements=(400, 40, 4), verification_passed=True, estimate=500.0,
                         gmsh_seconds=100.0, gmsh_gib=2.0):
        root.mkdir(parents=True)
        (root / "build-summary.json").write_text(json.dumps(
            {"Case": self.case["Id"], "Status": STATUS_BUILT if verification_passed else STATUS_FAILED,
             "Stage": "per-entry-verification", "ReturnCode": 0 if verification_passed else 1,
             "ScopeGuard": None, "Message": None if verification_passed else "per-entry verification failed"}))
        cap = self.gates["MaximumElements"]
        (root / "build-cost-estimate.json").write_text(json.dumps(
            {"EstimatedElements": estimate, "EstimateOverCap": estimate / cap, "MaximumElements": cap,
             "Passed": estimate <= cap}))
        (root / "build-census.json").write_text(json.dumps({"PrismTubes": {"Quality": {
            kind: {"Count": count} for kind, count in zip(("Tetrahedron", "Prism", "Pyramid"), elements)}}}))
        (root / "canonical-build.json").write_text(json.dumps({"CanonicalBuildId": "c" * 64}))
        (root / "audits").mkdir()
        for variant in self.case["Variants"]:
            (root / "audits" / f"{self.case['Id']}--{variant['Id']}.json").write_text(json.dumps(
                {"Mesh": {"Path": f"{root}/{variant['Id']}.msh", "SHA256": variant["Id"][0] * 64}}))
        (root / "per-entry-verification.json").write_text(json.dumps(
            {"Passed": verification_passed, "Failures": [] if verification_passed else ["identity: gate failure: x"],
             "SharedCanonicalBuildId": "c" * 64 if verification_passed else None}))
        stage_report(root, "canonical-source", 0.1, 0.0)
        stage_report(root, "gmsh-build", gmsh_seconds, gmsh_gib)
        stage_report(root, "canonical-publish", 10.0, 1.0)
        stage_report(root, "identity-publication", 20.0, 1.5)
        stage_report(root, "identity-variant-audits", 30.0, 1.0, audit=True)
        stage_report(root, "per-entry-verification", 200.0, 0.9, audit=True)

    def test_built_and_passed_record_schema(self):
        root = self.tmp / "built"
        self.write_built_root(root)
        record = self.record(root)
        self.assertEqual(set(record), {"Case", "InventoryStatus", "FixtureVersion", "Calibration", "Scope", "Status", "Passed",
                                       "StoppedBy", "CanonicalBuildId", "Variants", "Elements", "H1", "Estimate",
                                       "Stages", "Verification", "HeadroomFlags", "Root", "DriverReturnCode",
                                       "WallSeconds"})
        self.assertEqual((record["Status"], record["Passed"], record["StoppedBy"]), (STATUS_BUILT, True, None))
        self.assertEqual(record["Scope"]["UnsupportedClasses"], [])
        self.assertTrue(record["Scope"]["ExhibitedClasses"])
        self.assertEqual(record["CanonicalBuildId"], "c" * 64)
        self.assertEqual(record["Elements"], {"Tetrahedron": 400, "Prism": 40, "Pyramid": 4, "Total": 444})
        self.assertEqual(record["Variants"]["identity"]["SHA256"], "i" * 64)
        self.assertEqual(set(record["Variants"]), {variant["Id"] for variant in self.case["Variants"]})
        self.assertAlmostEqual(record["Estimate"]["EstimateOverActual"], 500.0 / 444)
        self.assertEqual(record["Estimate"]["ActualElements"], 444)
        self.assertEqual(record["Estimate"]["MaximumElements"], self.gates["MaximumElements"])
        self.assertEqual(set(record["Stages"]), {"canonical-source", "gmsh-build", "canonical-publish",
                                                 "identity-publication", "identity-variant-audits",
                                                 "per-entry-verification"})
        self.assertEqual(record["Stages"]["gmsh-build"]["PeakGiB"], 2.0)
        self.assertEqual(record["Stages"]["gmsh-build"]["Limits"], {"Seconds": 3600.0, "MemoryGiB": 12.0})
        self.assertTrue(record["Verification"]["Passed"])
        self.assertEqual(record["HeadroomFlags"], [])
        self.assertIsNone(record["H1"])   # no identity mesh in the synthetic root
        self.assertEqual(record["WallSeconds"], 12.5)

    def test_verification_failure_gate_guard_and_stage_stops_are_attributed_exactly(self):
        root = self.tmp / "verification"
        self.write_built_root(root, verification_passed=False)
        record = self.record(root, driver_return_code=1)
        self.assertEqual((record["Status"], record["Passed"]), (STATUS_FAILED, False))
        self.assertEqual(record["StoppedBy"]["Kind"], "Verification")
        self.assertEqual(record["StoppedBy"]["Failures"], ["identity: gate failure: x"])
        cap = self.gates["MaximumElements"]
        gate = self.tmp / "gate"; gate.mkdir()
        (gate / "build-summary.json").write_text(json.dumps(
            {"Status": STATUS_FAILED, "Stage": "headroom-gate", "ReturnCode": 1, "ScopeGuard": None,
             "Message": "pre-build element estimate exceeds the cap"}))
        (gate / "build-cost-estimate.json").write_text(json.dumps(
            {"EstimatedElements": 1.5 * cap, "EstimateOverCap": 1.5, "MaximumElements": cap, "Passed": False}))
        record = self.record(gate, driver_return_code=1)
        self.assertEqual(record["StoppedBy"]["Kind"], "HeadroomGate")
        self.assertEqual(record["StoppedBy"]["Id"], "MaximumElements")
        self.assertIsNone(record["Elements"])
        self.assertFalse(record["Estimate"]["Passed"])
        self.assertIn(f"estimate {1.5 * cap:.0f} >= {HEADROOM_FRACTION} x MaximumElements {cap}", record["HeadroomFlags"])
        guard = self.tmp / "guard"; guard.mkdir()
        (guard / "build-summary.json").write_text(json.dumps(
            {"Status": STATUS_UNSUPPORTED, "Stage": "gmsh-build", "ReturnCode": 1, "ScopeGuard": "TopRounding",
             "Message": "unsupported class TopRounding"}))
        stage_report(guard, "gmsh-build", 3.0, 0.5, code=1)
        record = self.record(guard, driver_return_code=1)
        self.assertEqual(record["Status"], STATUS_UNSUPPORTED)
        self.assertEqual(record["StoppedBy"], {"Kind": "ScopeGuard", "Id": "TopRounding", "Stage": "gmsh-build",
                                               "Message": "unsupported class TopRounding"})
        timeout = self.tmp / "timeout"; timeout.mkdir()
        (timeout / "build-summary.json").write_text(json.dumps(
            {"Status": STATUS_FAILED, "Stage": "gmsh-build", "ReturnCode": 124, "ScopeGuard": None,
             "Message": "stage gmsh-build failed rc=124"}))
        stage_report(timeout, "gmsh-build", STAGE_LIMITS[0] + 0.2, 3.0, code=-15, reason="timeout")
        record = self.record(timeout, driver_return_code=1)
        self.assertEqual((record["StoppedBy"]["Kind"], record["StoppedBy"]["Id"], record["StoppedBy"]["StopReason"]),
                         ("Stage", "gmsh-build", "timeout"))
        self.assertEqual(len(record["HeadroomFlags"]), 1)
        empty = self.tmp / "empty"; empty.mkdir()
        record = self.record(empty, driver_return_code=2)
        self.assertEqual((record["Status"], record["StoppedBy"]["Kind"]), (STATUS_FAILED, "Driver"))

    def test_headroom_flags_and_library_totals(self):
        cap, bound = self.gates["MaximumElements"], self.gates["MaximumRSSGiB"]
        root = self.tmp / "tight"
        count = int(HEADROOM_FRACTION * cap)
        self.write_built_root(root, elements=(count, 0, 0), estimate=0.5 * cap, gmsh_gib=HEADROOM_FRACTION * bound + 0.01)
        record = self.record(root)
        self.assertEqual(len(record["HeadroomFlags"]), 2, record["HeadroomFlags"])
        self.assertTrue(record["HeadroomFlags"][0].startswith(f"elements {count} >="))
        self.assertIn("stage gmsh-build", record["HeadroomFlags"][1])
        self.assertIn(f"{HEADROOM_FRACTION} x {float(bound)} GiB", record["HeadroomFlags"][1])
        loose = self.tmp / "loose"
        self.write_built_root(loose, elements=(count - 1, 0, 0), gmsh_gib=HEADROOM_FRACTION * bound - 0.01)
        self.assertEqual(self.record(loose)["HeadroomFlags"], [])
        guard = self.tmp / "guard"; guard.mkdir()
        (guard / "build-summary.json").write_text(json.dumps(
            {"Status": STATUS_UNSUPPORTED, "Stage": "gmsh-build", "ReturnCode": 1, "ScopeGuard": "TrenchRounding",
             "Message": "unsupported class TrenchRounding"}))
        records = [record, self.record(loose), self.record(guard, driver_return_code=1),
                   self.record(self.tmp / "missing", driver_return_code=1)]
        totals = library_totals(records, jobs=2, wall_seconds=99.0, manifest=self.manifest,
                                manifest_path=self.manifest_path, commit="abc")
        self.assertEqual({key: totals[key] for key in ("CasesAttempted", "CasesBuilt", "CasesPassed",
                                                       "CasesUnsupported", "CasesFailed")},
                         {"CasesAttempted": 4, "CasesBuilt": 2, "CasesPassed": 2, "CasesUnsupported": 1,
                          "CasesFailed": 1})
        self.assertEqual(totals["FlaggedCases"], [self.case["Id"]])
        self.assertEqual(totals["Bounds"], {"MaximumElements": cap, "MaximumSeconds": self.gates["MaximumSeconds"],
                                            "MaximumRSSGiB": bound})
        self.assertEqual((totals["Jobs"], totals["WallClockSeconds"], totals["Commit"]), (2, 99.0, "abc"))
        self.assertEqual(totals["Manifest"]["Path"], str(self.manifest_path))

    def test_matrix_rejects_unknown_cases_legacy_manifests_and_bad_limits(self):
        with self.assertRaises(ValueError):
            run_matrix(self.manifest_path, ["no-such-case"], root=self.tmp / "root")
        # A labeled Gmsh-only calibration manifest is accepted (decision 53 builds a labeled
        # ring-set variant through the matrix); the legacy-pipeline calibration manifest is not.
        with self.assertRaises(ValueError):
            run_matrix(HERE / "geometry-independence-calibration-ma.json", root=self.tmp / "root")
        with self.assertRaises(ValueError):
            run_matrix(self.manifest_path, [self.case["Id"]], root=self.tmp / "root", jobs=0)
        with self.assertRaises(ValueError):
            run_matrix(self.manifest_path, [self.case["Id"]], root=self.tmp / "root", build_limit=-1)


class LibraryBuildEndToEndTest(unittest.TestCase):
    """`coupon-library build` end to end: a temporary copy of the smallest gallery case
    (no repository contract copied) is registered through the real two-pass probe and
    built next to one-edge-straight as a pool of two; the record is kept under /tmp
    (`coupon-library-e2e-*`, the root printed) for the evidence chain."""

    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        probe = subprocess.run([cls.julia, "--startup-file=no", f"--project={REPO / 'test' / 'examples'}", "-e",
                                "import Gmsh: gmsh; gmsh.initialize(); gmsh.finalize()"],
                               capture_output=True, text=True)
        if probe.returncode != 0:
            raise unittest.SkipTest("the test/examples Julia project cannot load Gmsh")

    def test_register_and_build_two_small_fixtures(self):
        import coupon_library
        from register_case import REGISTRATION_RECORD
        manifest = json.loads(PRODUCTION_MANIFEST.read_text())
        # The smallest production cases: the gallery case with the fewest signature rows
        # among those binding a trace basis (registered as a copy) and the repository
        # fixture the scaling comparisons reference (built as itself).
        with_basis = [case for case in manifest["Cases"] if "BasisContract" in case["Source"]["Files"]]
        smallest = min(with_basis, key=lambda case: len(
            (REPO / case["Source"]["Directory"] / case["Source"]["Files"]["Signature"]["Name"]).read_text().splitlines()))
        fixture = manifest["ScalingComparisons"][0]["Reference"][0]
        root = Path(tempfile.mkdtemp(prefix=f"coupon-library-e2e-{time.strftime('%Y%m%d-%H%M%S')}-", dir="/tmp"))
        manifest["RepositoryRoot"] = str(REPO)
        manifest_path = root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        source = root / "source"
        shutil.copytree(REPO / smallest["Source"]["Directory"], source)
        (source / "semantic-contract.json").unlink()
        copy_id = f"{smallest['Id']}-copy"
        recipe = smallest["Source"]["Files"]["MeshRecipe"]
        argv = ["build", "--manifest", str(manifest_path), "--root", str(root / "library"), "--jobs", "2",
                "--register", f"{copy_id}={source}", "--footprint", "producer-default",
                "--inventory-status", "RepositoryAssessmentFixture", "--work", str(root / "register"),
                "--case", fixture, "--julia", self.julia]
        if "RepositoryPath" in recipe:
            argv += ["--mesh-recipe", recipe["RepositoryPath"]]
        print(f"\ncoupon-library e2e root {root}", flush=True)
        code = coupon_library.main(argv)
        record = json.loads((root / "library" / LIBRARY_BUILD_RECORD).read_text())
        self.assertEqual(code, 0, json.dumps(record["Library"], indent=1))
        # Registration: the real two-pass derivation reproduces the repository contract.
        registration = json.loads((root / "register" / copy_id / REGISTRATION_RECORD).read_text())
        self.assertEqual((registration["Status"], registration["FixtureVersion"]), ("registered", 1))
        frozen = json.loads((REPO / smallest["Source"]["Directory"] / "semantic-contract.json").read_text())
        derived = json.loads((source / "semantic-contract.json").read_text())
        for key in frozen:
            if key != "Derivation":
                self.assertEqual(derived[key], frozen[key], key)
        registered = json.loads(manifest_path.read_text())["Cases"][-1]
        self.assertEqual(registered["Id"], copy_id)
        # Build: both cases built and verified as one library of two jobs.
        totals = record["Library"]
        self.assertEqual({key: totals[key] for key in ("CasesAttempted", "CasesBuilt", "CasesPassed",
                                                       "CasesUnsupported", "CasesFailed")},
                         {"CasesAttempted": 2, "CasesBuilt": 2, "CasesPassed": 2, "CasesUnsupported": 0,
                          "CasesFailed": 0})
        self.assertEqual(totals["Jobs"], 2)
        by_case = {case["Case"]: case for case in record["Cases"]}
        self.assertEqual(set(by_case), {fixture, copy_id})
        for case in record["Cases"]:
            self.assertTrue(case["Passed"], case["StoppedBy"])
            self.assertIsNone(case["StoppedBy"])
            self.assertEqual(len(case["CanonicalBuildId"]), 64)
            for variant in case["Variants"].values():
                self.assertEqual(len(variant["SHA256"]), 64)
                self.assertTrue(Path(variant["Path"]).is_file())
            self.assertGreater(case["Elements"]["Prism"], 0)
            self.assertGreater(case["Elements"]["Pyramid"], 0)
            self.assertEqual(case["Elements"]["Total"], sum(case["Elements"][k] for k in ("Tetrahedron", "Prism", "Pyramid")))
            self.assertEqual(case["H1"]["Order"], 4)
            self.assertGreater(case["H1"]["DOFs"], case["Elements"]["Total"])
            counts = case["H1"]["EntityCounts"]
            self.assertEqual((counts["Tetrahedra"], counts["Prisms"], counts["Pyramids"]),
                             (case["Elements"]["Tetrahedron"], case["Elements"]["Prism"], case["Elements"]["Pyramid"]))
            from mixed_mesh import h1_dofs_from_counts
            self.assertEqual(h1_dofs_from_counts(counts, 4), case["H1"]["DOFs"])
            self.assertLess(case["Estimate"]["EstimateOverCap"], 1.0)
            self.assertGreater(case["Estimate"]["EstimateOverActual"], 0.0)
            self.assertTrue(case["Verification"]["Passed"])
            for stage in ("canonical-source", "gmsh-build", "canonical-publish", "identity-publication",
                          "identity-variant-audits", "per-entry-verification"):
                self.assertIn(stage, case["Stages"])
                self.assertEqual(case["Stages"][stage]["Limits"],
                                 {"Seconds": float(manifest["Gates"]["MaximumSeconds"]),
                                  "MemoryGiB": float(manifest["Gates"]["MaximumRSSGiB"])})
                self.assertEqual(case["Stages"][stage]["ReturnCode"], 0)
        # The copy binds exactly the gallery case's sources (the contract is the derived one).
        self.assertEqual(json.loads((root / "library" / copy_id / "input-hashes.json").read_text())["Signature"],
                         smallest["Source"]["Files"]["Signature"]["SHA256"])


if __name__ == "__main__":
    unittest.main()
