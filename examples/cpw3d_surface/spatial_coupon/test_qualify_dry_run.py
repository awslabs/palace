# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""coupon-library qualify --dry-run end to end on the four-edge case with the recorded
physics-11 configuration (controls, stage prefix) and on the gallery-06 case: the
generated configs equal the recorded worker / reducer / local-edge configs apart from
paths, the plan equals the recorded plan in stages and pins; the analysis of the
recorded results through the same records gives the recorded verdicts; the per-coupon
fail-closed stops.  Needs a local identity mesh of each case and the assessment tree."""
import argparse
import glob
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import estimate_stages  # noqa: E402
import gates  # noqa: E402
import qualify_library  # noqa: E402
from mixed_mesh import h1_dofs_from_counts  # noqa: E402
from run_gmsh_only_matrix import sha256  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
MANIFEST = HERE / "geometry-independence-suite.json"
BINARY_SHA256 = "b28f089ae12c25863493566b2b8ca11af2c8ffb0e273e7aa67a2b42046eacf27"
CASES = {
    "four-edge-9d2cb9bbb3fe": {"Campaign": "four-edge-physics-11", "Prefix": "va", "Controls": [1, 7, 23, 26, 34, 35, 48, 80],
                               "Sources": 80, "Stages": ["va-p4", "va-p5-control", "va-p3-control", "va-p4-local-edge"]},
    "three-edge-419576fdab24": {"Campaign": "gallery-physics-06b", "Prefix": "g06b", "Controls": [14, 25, 26, 33, 52, 99, 133, 134],
                                "Sources": 135, "Stages": ["g06b-p4", "g06b-p5-control", "g06b-p3-control", "g06b-p4-local-edge"]},
}


def local_identity_mesh(case_id):
    """The newest local Gmsh-only root of the case with an identity mesh (None when absent)."""
    roots = sorted(glob.glob(f"/tmp/coupon-gmsh-only-{case_id}-*/identity.msh"), key=os.path.getmtime)
    return Path(roots[-1]) if roots else None


def available():
    return all(local_identity_mesh(case_id) is not None and (ASSESSMENT / spec["Campaign"] / "results" / "main" / "status.json").is_file()
               for case_id, spec in CASES.items())


def strip_paths(value):
    if isinstance(value, dict):
        return {key: strip_paths(item) for key, item in value.items()}
    if isinstance(value, list):
        return [strip_paths(item) for item in value]
    if isinstance(value, str) and "/" in value:
        return "<path>"
    return value


@unittest.skipUnless(available(), "local identity meshes and the assessment campaigns are needed")
class QualifyDryRunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coupon-qualify-test-"))
        cases = []
        for case_id in CASES:
            mesh = local_identity_mesh(case_id)
            counts = estimate_stages.entity_counts_of_mesh(mesh)
            cases.append({"Case": case_id, "Status": "built", "Passed": True, "CanonicalBuildId": None,
                          "Variants": {"identity": {"Path": str(mesh), "SHA256": sha256(mesh)}},
                          "Elements": {"Tetrahedron": counts["Tetrahedra"], "Prism": counts["Prisms"],
                                       "Pyramid": counts["Pyramids"],
                                       "Total": counts["Tetrahedra"] + counts["Prisms"] + counts["Pyramids"]},
                          "H1": {"Order": 4, "DOFs": h1_dofs_from_counts(counts, 4), "EntityCounts": counts},
                          "StoppedBy": None, "HeadroomFlags": [], "Root": str(mesh.parent)})
        commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True).strip()
        cls.build_record = cls.tmp / "library-build.json"
        cls.build_record.write_text(json.dumps(
            {"Version": 1, "Command": "coupon-library build", "Root": str(cls.tmp), "Cases": cases,
             "Library": {"Commit": commit, "Manifest": {"Path": str(MANIFEST), "SHA256": sha256(MANIFEST)}}}, indent=2))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, True)

    def dry_run(self, case_id, root, extra=()):
        spec = CASES[case_id]
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", str(ASSESSMENT / spec["Campaign"] / "reference"),
                   "--remote", "soca-green-job:/data/home/simlap/coupon_accuracy_assessment_20260913",
                   "--orders", "p4", "--controls", "p3,p5", "--max-jobs", "2", "--frozen-binary-sha256", BINARY_SHA256,
                   "--case", case_id, "--stage-prefix", spec["Prefix"], "--root", str(root), "--dry-run", *extra]
        for control in spec["Controls"]:
            command += ["--control-source", str(control)]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads((root / "library-qualification.json").read_text())

    def check_against_campaign(self, case_id):
        spec = CASES[case_id]
        campaign = ASSESSMENT / spec["Campaign"]
        root = self.tmp / f"dry-{spec['Prefix']}"
        record = self.dry_run(case_id, root)
        case = record["Cases"][0]
        self.assertEqual(case["Status"], "planned", case.get("StoppedBy"))
        self.assertIsNone(case["StoppedBy"])
        self.assertEqual(case["Sources"]["Count"], spec["Sources"])
        self.assertEqual(case["Controls"]["Indices"], spec["Controls"])
        self.assertTrue(case["Estimate"]["FitsOneJob"])
        self.assertTrue(case["Mesh"]["Verified"])
        # Configs equal the recorded ones apart from paths (mesh, output, trace directory).
        for stage in spec["Stages"]:
            names = ("config.json",) if stage.endswith("local-edge") else ("worker.json", "reducer.json")
            for name in names:
                recorded = strip_paths(json.loads((campaign / "main" / stage / name).read_text()))
                generated = strip_paths(json.loads((root / case_id / "main" / stage / name).read_text()))
                self.assertEqual(generated, recorded, f"{stage}/{name}")
        # The plan equals the recorded one in stages (names, config file, environment,
        # dependencies, order) and pins (the same files; every trace digest equal; the mesh pin
        # is the build record's).
        recorded_plan = json.loads((campaign / "main" / "plan.json").read_text())
        plan = json.loads((root / case_id / "main" / "plan.json").read_text())

        def stage_view(stage):
            environment = {key: (Path(value).name if key == "PALACE_RESPONSE_ARCHIVE_DIR" else value)
                           for key, value in stage["Environment"].items()}
            return (stage["Name"], Path(stage["Config"]).name, environment, stage["Requires"])
        self.assertEqual([stage_view(s) for s in plan["Stages"]], [stage_view(s) for s in recorded_plan["Stages"]])
        recorded_pins = {Path(key).name: value for key, value in recorded_plan["PinnedSHA256"].items()}
        pins = {Path(key).name: value for key, value in plan["PinnedSHA256"].items()}
        self.assertEqual(len(pins), len(recorded_pins))
        traces = {name: digest for name, digest in recorded_pins.items() if name.startswith("basis-")}
        self.assertEqual(len(traces), spec["Sources"])
        self.assertEqual({name: pins[name] for name in traces}, traces)
        self.assertEqual(set(pins) - set(recorded_pins), {Path(plan["MeshRemote"]).name})
        self.assertEqual(plan["MeshSHA256"], sha256(local_identity_mesh(case_id)))
        self.assertEqual(plan["Ranks"], recorded_plan["Ranks"])
        self.assertEqual(plan["DeadlineSeconds"], recorded_plan["DeadlineSeconds"])
        self.assertEqual(plan["MinimumMemAvailableBytes"], recorded_plan["MinimumMemAvailableBytes"])
        self.assertEqual(plan["BinarySHA256"], BINARY_SHA256)
        self.assertTrue(plan["Binary"].endswith(f"palace-archive-estimate-{BINARY_SHA256}.bin"))
        for stage in plan["Stages"]:
            self.assertGreaterEqual(stage["CapSeconds"], stage["MinimumSeconds"])
            self.assertLessEqual(stage["CapSeconds"], plan["DeadlineSeconds"])
        self.assertTrue((root / case_id / "main" / "job.pbs").is_file())
        self.assertTrue((root / "qualification-gates.json").is_file())
        self.assertEqual(record["Gates"]["SHA256"], sha256(HERE / "qualify" / "qualification-gates.json"))
        self.assertTrue(record["Library"]["DryRun"])
        self.assertEqual(record["Library"]["JobsSubmitted"], 0)
        self.assertEqual(record["Library"]["CouponsPlanned"], 1)
        return root, record

    def test_four_edge_reproduces_physics_11(self):
        self.check_against_campaign("four-edge-9d2cb9bbb3fe")

    def test_gallery_06_reproduces_physics_06b(self):
        self.check_against_campaign("three-edge-419576fdab24")

    def test_controls_by_class_when_not_named(self):
        root = self.tmp / "dry-by-class"
        spec = CASES["three-edge-419576fdab24"]
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", str(ASSESSMENT / spec["Campaign"] / "reference"), "--frozen-binary-sha256", BINARY_SHA256,
                   "--case", "three-edge-419576fdab24", "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        case = json.loads((root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(len(case["Controls"]["Indices"]), 8)
        self.assertIn("by class", case["Controls"]["Rule"])
        self.assertEqual(len(set(case["Controls"]["Classes"].values())), 7)   # every non-zero-trace class of the case
        self.assertEqual(case["StagePrefix"], "three-edge-419576fdab24")
        self.assertEqual(case["Plan"]["StageNames"][0], "three-edge-419576fdab24-p4-worker")

    def analysis_context(self, case_id, reference, controls=None):
        spec = CASES[case_id]
        build = json.loads(self.build_record.read_text())
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        args = argparse.Namespace(reference=reference, control_source=controls or spec["Controls"], control_count=8,
                                  stage_prefix=spec["Prefix"], orders=[4], controls=[3, 5], frozen_binary_sha256=BINARY_SHA256)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        model = estimate_stages.load_cost_model()
        table, digest = gates.load_gates()
        root = self.tmp / f"analysis-{spec['Prefix']}-{Path(str(reference)).name}"
        root.mkdir(exist_ok=True)
        case = next(item for item in build["Cases"] if item["Case"] == case_id)
        record, context = qualify_library.prepare_case(case, manifest_path=MANIFEST, manifest=manifest, args=args, root=root,
                                                       remote={"Host": "h", "Root": "/r"}, profile=profile, cost_model=model,
                                                       gates=table, gates_digest=digest)
        return record, context, root, table, digest, profile, manifest

    def assert_campaign_untouched(self, campaign, before):
        after = {path: path.stat().st_mtime_ns for path in (campaign / "results").rglob("*") if path.is_file()}
        self.assertEqual(after, before, "the recorded campaign's results tree must stay read only")

    def test_analysis_of_the_recorded_results_reproduces_the_verdicts(self):
        for case_id, spec in CASES.items():
            campaign = ASSESSMENT / spec["Campaign"]
            before = {path: path.stat().st_mtime_ns for path in (campaign / "results").rglob("*") if path.is_file()}
            record, context, root, table, digest, profile, manifest = self.analysis_context(case_id, campaign / "reference")
            gate_record = qualify_library.analyze_case(record, context, campaign / "results", gates=table, gates_digest=digest,
                                                       profile=profile)
            self.assert_campaign_untouched(campaign, before)
            self.assertEqual(Path(record["Cost"]["Path"]), root / case_id / "cost-summary.json")
            self.assertTrue((root / case_id / "cost-summary.json").is_file())
            self.assertEqual(gate_record["Verdict"], gates.VERDICT_PASSED, gate_record["Reason"])
            self.assertEqual(record["Status"], "qualified")
            self.assertEqual(record["Qualification"]["ReferenceAnchor"], "vs p4 anchor")
            self.assertAlmostEqual(record["Cost"]["MainStage"]["NodeHours"],
                                   0.508 if case_id.startswith("four") else 0.974, places=3)
            self.assertLess(record["Cost"]["MainStageOverReference"], 0.3)
            recorded_plan = json.loads((campaign / "main" / "plan.json").read_text())
            if recorded_plan["MeshSHA256"] == record["Mesh"]["SHA256"]:
                # Same mesh as the campaign: the Palace-printed H1 equals the closed-form estimate.
                self.assertEqual(record["Cost"]["MainStage"]["H1"], record["Estimate"]["H1ByOrder"]["p4"])
            comparison = root / case_id / "comparison"
            for name in ("source-classes.csv", "class-statistics.md", "ma-ms-offsets.json", "p-sequence-controls.json",
                         "key-sources.md", f"{spec['Prefix']}-p4-vs-reference.json"):
                self.assertTrue((comparison / name).is_file(), name)
            library = qualify_library.process_library_entries([record], {case_id: context}, manifest_path=MANIFEST,
                                                              manifest=manifest, root=root)
            self.assertEqual(len(library["Models"]), 1)
            self.assertTrue(library["Models"][0]["LibraryQualified"])
            self.assertEqual(library["Models"][0]["Qualification"]["Verdict"], gates.VERDICT_PASSED)
            self.assertEqual(library["Models"][0]["CouponMesh"]["SHA256"], record["Mesh"]["SHA256"])

    def test_without_reference_matrices_the_verdict_is_pending(self):
        case_id = "four-edge-9d2cb9bbb3fe"
        campaign = ASSESSMENT / CASES[case_id]["Campaign"]
        inputs_only = self.tmp / "inputs-only-campaign"
        if not inputs_only.exists():
            inputs_only.mkdir()
            os.symlink(campaign / "reference" / "inputs-07", inputs_only / "inputs-07")
        record, context, root, table, digest, profile, manifest = self.analysis_context(case_id, inputs_only)
        self.assertIsNone(record["Reference"]["Results"])
        before = {path: path.stat().st_mtime_ns for path in (campaign / "results").rglob("*") if path.is_file()}
        gate_record = qualify_library.analyze_case(record, context, campaign / "results", gates=table, gates_digest=digest,
                                                   profile=profile)
        self.assert_campaign_untouched(campaign, before)
        self.assertEqual(gate_record["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(set(gate_record["Gates"]), {"PSequenceControls"})
        self.assertEqual(record["Status"], "pending-qualification")
        self.assertIsNone(record["Cost"]["ReferenceNodeHours"])
        library = qualify_library.process_library_entries([record], {case_id: context}, manifest_path=MANIFEST,
                                                          manifest=manifest, root=root)
        self.assertFalse(library["Models"][0]["LibraryQualified"])

    def test_fail_closed_stops(self):
        case_id = "four-edge-9d2cb9bbb3fe"
        campaign = ASSESSMENT / CASES[case_id]["Campaign"]
        build = json.loads(self.build_record.read_text())
        # A coupon the build did not pass, a coupon with no reference inputs, a corrupt mesh digest
        # and a coupon that does not fit the walltime are recorded stops, never crashes.
        cases = json.loads(json.dumps(build["Cases"]))
        cases[0].update(Passed=False, Status="failed")
        cases[1]["Variants"]["identity"]["SHA256"] = "0" * 64
        stopped = self.tmp / "library-build-stopped.json"
        stopped.write_text(json.dumps({**build, "Cases": cases}, indent=2))
        root = self.tmp / "dry-stopped"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(stopped),
                   "--reference", str(campaign / "reference"), "--frozen-binary-sha256", BINARY_SHA256, "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        record = json.loads((root / "library-qualification.json").read_text())
        by_case = {case["Case"]: case for case in record["Cases"]}
        self.assertEqual(by_case["four-edge-9d2cb9bbb3fe"]["StoppedBy"]["Kind"], "Build")
        self.assertEqual(by_case["four-edge-9d2cb9bbb3fe"]["Status"], "skipped")
        self.assertEqual(by_case["three-edge-419576fdab24"]["StoppedBy"]["Kind"], "Mesh")
        self.assertEqual(by_case["three-edge-419576fdab24"]["Status"], "failed")
        self.assertEqual(record["Library"]["CouponsSkipped"], 1)
        self.assertEqual(record["Library"]["CouponsFailed"], 1)
        # No reference inputs for the case (another campaign) -> skipped with the reason.
        root = self.tmp / "dry-no-inputs"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", str(ASSESSMENT / "gallery-physics-06b" / "reference"), "--frozen-binary-sha256", BINARY_SHA256,
                   "--case", case_id, "--root", str(root), "--dry-run"]
        subprocess.run(command, text=True, capture_output=True)
        case = json.loads((root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(case["StoppedBy"]["Kind"], "Reference")
        self.assertIn("no reference campaign inputs", case["StoppedBy"]["Message"])
        # A walltime the coupon cannot fit -> the estimate gate stops it before any plan.
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        profile["WalltimeSeconds"] = 600
        short = self.tmp / "short-profile.json"
        short.write_text(json.dumps(profile))
        root = self.tmp / "dry-short"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", str(campaign / "reference"), "--frozen-binary-sha256", BINARY_SHA256, "--case", case_id,
                   "--root", str(root), "--dry-run", "--cluster-profile", str(short)]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 1)
        case = json.loads((root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(case["StoppedBy"]["Kind"], "Estimate")
        self.assertIn("does NOT fit", case["StoppedBy"]["Message"])
        self.assertFalse((root / case_id / "main" / "plan.json").exists())
        # Without --dry-run the remote is mandatory.
        result = subprocess.run([sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                                 "--frozen-binary-sha256", BINARY_SHA256, "--root", str(self.tmp / "no-remote")],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--remote", result.stderr)


if __name__ == "__main__":
    unittest.main()
