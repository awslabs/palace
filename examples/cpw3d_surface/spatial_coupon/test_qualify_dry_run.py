# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""coupon-library qualify --dry-run end to end on the four-edge case with the recorded
physics-11 configuration (controls, stage prefix) and on the gallery-06 case: the
generated configs equal the recorded worker / reducer / local-edge configs apart from
paths, the plan equals the recorded plan in stages and pins; the analysis of the
recorded results through the same records gives the recorded verdicts (gallery-10: the
reference order p5 added as a main stage and gated, p_SA not applicable); the concurrent
job scheduler against a fake remote replaying the recorded trees; the per-coupon
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


# The recorded two-edge campaign: reference at p5 (a main stage the command adds), no SA interface.
GALLERY_10 = {"Case": "two-edge-8dd4bc70f183", "Campaign": "gallery-physics-10", "Prefix": "g10",
              "Controls": [1, 7, 21, 25, 26, 35, 43, 78], "Sources": 78}


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

    @unittest.skipUnless(local_identity_mesh(GALLERY_10["Case"]) is not None
                         and (ASSESSMENT / GALLERY_10["Campaign"] / "results" / "main" / "status.json").is_file(),
                         "the two-edge mesh and the gallery-10 campaign are needed")
    def test_recorded_gallery_10_is_gated_at_the_reference_order_without_sa(self):
        """gallery-physics-10 (reference at p5; MA / MS only): --orders p4 runs p4 AND p5 as main
        stages, the p5 (same-order) comparison is gated and reproduces RESULTS.md (E 78/78,
        p_MA 33/59/78 with the strongest-20 failing at 53 / 58, p_MS 75/78/78), p_SA is
        NotApplicable (not a failure), the p4 comparison is informational, both costs recorded."""
        spec = GALLERY_10
        mesh = local_identity_mesh(spec["Case"])
        counts = estimate_stages.entity_counts_of_mesh(mesh)
        build = json.loads(self.build_record.read_text())
        build["Cases"] = [{"Case": spec["Case"], "Status": "built", "Passed": True, "CanonicalBuildId": None,
                           "Variants": {"identity": {"Path": str(mesh), "SHA256": sha256(mesh)}},
                           "Elements": {"Tetrahedron": counts["Tetrahedra"], "Prism": counts["Prisms"], "Pyramid": counts["Pyramids"],
                                        "Total": counts["Tetrahedra"] + counts["Prisms"] + counts["Pyramids"]},
                           "H1": {"Order": 4, "DOFs": h1_dofs_from_counts(counts, 4), "EntityCounts": counts},
                           "StoppedBy": None, "HeadroomFlags": [], "Root": str(mesh.parent)}]
        campaign = ASSESSMENT / spec["Campaign"]
        args = argparse.Namespace(reference=campaign / "reference", control_source=spec["Controls"], control_count=8,
                                  stage_prefix=spec["Prefix"], orders=[4], controls=[3, 5], frozen_binary_sha256=BINARY_SHA256)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        model = estimate_stages.load_cost_model()
        table, digest = gates.load_gates()
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        root = self.tmp / "analysis-g10"
        root.mkdir(exist_ok=True)
        before = {path: path.stat().st_mtime_ns for path in (campaign / "results").rglob("*") if path.is_file()}
        record, context = qualify_library.prepare_case(build["Cases"][0], manifest_path=MANIFEST, manifest=manifest, args=args,
                                                       root=root, remote={"Host": "h", "Root": "/r"}, profile=profile,
                                                       cost_model=model, gates=table, gates_digest=digest)
        self.assertEqual(record["Orders"]["Main"], ["p4", "p5"])
        self.assertEqual(record["Orders"]["Gated"], "p5")
        self.assertEqual(record["Reference"]["Interfaces"], ["MA", "MS"])
        self.assertEqual(record["Plan"]["StageNames"], ["g10-p4-worker", "g10-p4-reducer", "g10-p5-worker", "g10-p5-reducer",
                                                        "g10-p5-control-worker", "g10-p5-control-reducer",
                                                        "g10-p3-control-worker", "g10-p3-control-reducer", "g10-p4-local-edge"])
        recorded_plan = json.loads((campaign / "main" / "plan.json").read_text())
        self.assertEqual(record["Plan"]["StageNames"], [stage["Name"] for stage in recorded_plan["Stages"]])
        gate_record = qualify_library.analyze_case(record, context, campaign / "results", gates=table, gates_digest=digest,
                                                   profile=profile)
        self.assert_campaign_untouched(campaign, before)
        self.assertEqual(gate_record["Verdict"], gates.VERDICT_FAILED)
        self.assertEqual(gate_record["Reason"], "failing gates ['p_MA'] (p_SA not applicable)")
        self.assertEqual(gate_record["GatedOrder"], 5)
        self.assertEqual(record["Qualification"]["GatedStage"], "g10-p5")
        self.assertEqual(record["Qualification"]["ReferenceAnchor"], "vs p5 anchor")
        self.assertEqual(record["Qualification"]["NotApplicable"], ["p_SA"])
        self.assertTrue(gate_record["GatesPassed"]["p_SA"])
        self.assertTrue(gate_record["GatesPassed"]["PSequenceControls"])
        self.assertEqual(gate_record["Gates"]["PSequenceControls"]["NotApplicableObservables"], ["p_SA"])
        energy, ma, ms = gate_record["Gates"]["E"], gate_record["Gates"]["p_MA"], gate_record["Gates"]["p_MS"]
        self.assertEqual((energy["AllFree"]["n"], energy["AllFree"]["within_1pct"]), (78, 78))
        self.assertEqual((ma["Free"]["within_1pct"], ma["Free"]["within_2pct"], ma["Free"]["within_5pct"]), (33, 59, 78))
        self.assertEqual(ma["StrongestFailing"], [53, 58])
        self.assertAlmostEqual(ma["Free"]["signed_median"], 0.0128, places=4)
        self.assertEqual((ms["Free"]["within_1pct"], ms["Free"]["within_2pct"]), (75, 78))
        self.assertEqual(list(record["Qualification"]["Informational"]), ["g10-p4-vs-reference"])
        self.assertAlmostEqual(record["Cost"]["MainStages"]["g10-p4"]["NodeHours"], 0.118, places=3)
        self.assertAlmostEqual(record["Cost"]["MainStages"]["g10-p5"]["NodeHours"], 0.3835, places=3)
        # The recorded campaign ran the a22b471c1 mesh (7,915,021 H1 at p4); the local production mesh may differ.
        self.assertEqual(record["Cost"]["MainStage"]["H1"], 7915021)
        self.assertEqual(record["Status"], "failed")
        library = qualify_library.process_library_entries([record], {spec["Case"]: context}, manifest_path=MANIFEST,
                                                          manifest=manifest, root=root)
        self.assertFalse(library["Models"][0]["LibraryQualified"])

    def test_jobs_run_concurrently_up_to_max_jobs(self):
        """Two planned coupons, --max-jobs 2, a fake remote that replays the recorded trees: both
        jobs are submitted before either is polled done, every active job is polled each
        round, each coupon is fetched / verified / analyzed when its job leaves the queue,
        and the totals carry the measured critical path and the job count."""
        events = []
        finish_after = {"four-edge-9d2cb9bbb3fe": 2, "three-edge-419576fdab24": 1}
        polls = {}

        def fake_upload(record, context, *, remote, profile):
            events.append(("upload", record["Case"]))
            return {"Commands": [], "UTC": "fake"}

        def fake_submit(host, pbs_bin, script, cwd, *, job_cap, user=None):
            case = Path(cwd).parts[-2]
            events.append(("submit", case))
            return {"Job": f"{len(events)}.fake", "UTC": qualify_library.remote_side.utc(), "UserJobsBefore": 0, "JobCap": job_cap,
                    "Command": "qsub"}

        def fake_poll(host, pbs_bin, job_id, status_path):
            case = Path(status_path).parts[-3]
            polls[case] = polls.get(case, 0) + 1
            events.append(("poll", case))
            state = "F" if polls[case] >= finish_after[case] else "R"
            return {"UTC": "fake", "JobState": state, "QStat": "", "Status": None}

        def fake_fetch(host, remote_directory, local_directory):
            case = Path(remote_directory).parts[-2]
            events.append(("fetch", case))
            shutil.copytree(ASSESSMENT / CASES[case]["Campaign"] / "results" / "main", local_directory, dirs_exist_ok=True)
            return ["rsync", "fake"]

        def fake_remote_sha256(host, paths):
            digests = {}
            for path in paths:
                case = Path(path).parts[-4] if "reducer" in path else Path(path).parts[-5]
                for candidate in ("four-edge-9d2cb9bbb3fe", "three-edge-419576fdab24"):
                    if candidate in path:
                        case = candidate
                local = self.tmp / "concurrent" / case / "results" / "main" / path.split("/main/", 1)[1]
                digests[path] = sha256(local)
            return digests

        def fake_delete(host, archives):
            events.append(("delete", Path(archives[0]).parts[-4]))
            return {"Archives": list(archives), "SizesBeforeDeletion": "0", "DeletedUTC": "fake", "Remaining": ""}

        fakes = {"submit": fake_submit, "poll": fake_poll, "fetch": fake_fetch, "remote_sha256": fake_remote_sha256,
                 "delete_archives": fake_delete, "qstat_history": lambda host, pbs_bin, job: "job_state = F"}
        saved = {name: getattr(qualify_library.remote_side, name) for name in fakes}
        saved_upload = qualify_library.upload_case
        controls = {case_id: spec["Controls"] for case_id, spec in CASES.items()}
        saved_prepare = qualify_library.prepare_case

        def prepare_with_recorded_controls(case_record, **kwargs):
            spec = CASES[case_record["Case"]]
            kwargs["args"].control_source = spec["Controls"]
            kwargs["args"].stage_prefix = spec["Prefix"]
            return saved_prepare(case_record, **kwargs)

        try:
            for name, fake in fakes.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = fake_upload
            qualify_library.prepare_case = prepare_with_recorded_controls
            args = argparse.Namespace(build_record=self.build_record, reference=None, remote="h:/r", orders=[4], controls=[3, 5],
                                      control_count=8, control_source=None, max_jobs=2, frozen_binary_sha256=BINARY_SHA256,
                                      stage_prefix=None, case=None, root=self.tmp / "concurrent", dry_run=False, resume=False,
                                      monitor_interval=0, monitor_polls=10, cluster_profile=HERE / "qualify" / "cluster-profile.json",
                                      cost_model=estimate_stages.COST_MODEL, gates=gates.GATES_FILE)
            # Each coupon binds its own reference campaign: a directory holding both inputs trees.
            reference = self.tmp / "both-references"
            if not reference.exists():
                reference.mkdir()
                for spec in CASES.values():
                    for entry in (ASSESSMENT / spec["Campaign"] / "reference").iterdir():
                        if entry.name.startswith(("inputs-", "case-")):
                            os.symlink(entry, reference / entry.name)
            args.reference = reference
            record = qualify_library.run_qualify(args, log=lambda message: None)
            # --resume on the same root adopts the recorded job ids: no upload, no qsub, the
            # plans are re-derived byte-identical, both coupons are polled / fetched / analyzed again.
            first_events = list(events)
            events.clear()
            polls.clear()
            args.resume = True
            resumed = qualify_library.run_qualify(args, log=lambda message: None)
        finally:
            for name, fake in saved.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = saved_upload
            qualify_library.prepare_case = saved_prepare
        del controls
        events, resumed_events = first_events, events
        kinds = [kind for kind, _ in events]
        self.assertEqual(kinds[:4], ["upload", "submit", "upload", "submit"], events)
        first_fetch = kinds.index("fetch")
        self.assertEqual(kinds.count("submit"), 2)
        self.assertLess(kinds.index("submit", kinds.index("submit") + 1), first_fetch, "both jobs queued before any fetch")
        # Round 1 polls both (three-edge done -> fetched); round 2 polls the four-edge job alone.
        self.assertEqual([case for kind, case in events if kind == "poll"],
                         ["four-edge-9d2cb9bbb3fe", "three-edge-419576fdab24", "four-edge-9d2cb9bbb3fe"])
        self.assertEqual([case for kind, case in events if kind == "fetch"], ["three-edge-419576fdab24", "four-edge-9d2cb9bbb3fe"])
        self.assertEqual([case for kind, case in events if kind == "delete"], ["three-edge-419576fdab24", "four-edge-9d2cb9bbb3fe"])
        totals = record["Library"]
        self.assertEqual(totals["JobsSubmitted"], 2)
        self.assertEqual(totals["MaxJobs"], 2)
        self.assertEqual(totals["CouponsQualified"], 2)
        self.assertIsNotNone(totals["CriticalPathSeconds"])
        self.assertGreaterEqual(totals["CriticalPathSeconds"], 0.0)
        self.assertEqual(set(totals["JobWallSeconds"]), set(CASES))
        self.assertAlmostEqual(totals["NodeHours"], sum(case["Cost"]["JobNodeHours"] for case in record["Cases"]))
        for case in record["Cases"]:
            self.assertEqual(case["Status"], "qualified", case.get("StoppedBy"))
            self.assertEqual(case["Qualification"]["Verdict"], gates.VERDICT_PASSED)
            self.assertEqual(case["Monitor"]["LastJobState"], "F")
            self.assertTrue(all(entry["OK"] for entry in case["ResultDigests"].values()))
            self.assertIn("ArchiveDeletion", case)
        self.assertTrue((self.tmp / "concurrent" / "process-library.json").is_file())
        resumed_kinds = [kind for kind, _ in resumed_events]
        self.assertNotIn("upload", resumed_kinds)
        self.assertNotIn("submit", resumed_kinds)
        self.assertEqual(resumed_kinds.count("fetch"), 2)
        self.assertEqual(resumed["Library"]["JobsSubmitted"], 2)
        self.assertEqual(resumed["Library"]["CouponsQualified"], 2)
        for case in resumed["Cases"]:
            self.assertTrue(case["Monitor"]["Resumed"])
            self.assertTrue(case["Upload"]["Resumed"])
            self.assertEqual(case["Submission"], json.loads((self.tmp / "concurrent" / case["Case"] / "submission.json").read_text()))
        self.assertIsNotNone(resumed["Library"]["CriticalPathSeconds"])

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
