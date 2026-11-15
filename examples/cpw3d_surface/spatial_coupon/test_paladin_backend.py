# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (USER decision 1674 (2); Y-SPEC s2; decision 1799 (2)):
qualify/paladin_backend.py runs a run_stages plan as paladin runs (one per stage, one node).
Checked here against the FAKE paladin CLI fake_paladin_cli.py and fake pin / schema / zero-hole
tools (HANDBOOK s18: these are synthetic; the live probe and the rehearsal are the real rows):
the variable -> key table; the Paladin block and the one-route rule; the E-2 sizing; the lookup
decision; and the drive: order, checks, outputs, status.json, resume by metadata (a driver killed
after a submit, while waiting, after a run completed), the five-condition completion read, the
cap cancel, out-of-memory, transport failures (one re-submit), refusals, the pin / schema /
archive / configuration checks, two runs of one stage; and, when the records of key
274fe1a150c0 are present on this host, the dry run of its fab plan of record."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import paladin_backend  # noqa: E402

ARCHIVE = "/data/case/p4/archive"
WORKER = {"PALACE_RESPONSE_ARCHIVE_DIR": ARCHIVE, "PALACE_RESPONSE_ARCHIVE_ONLY": "1",
          "PALACE_RESPONSE_SOURCE_TIMING": "1"}
REDUCER = {"PALACE_RESPONSE_ARCHIVE_DIR": ARCHIVE, "PALACE_RESPONSE_REDUCE_ONLY": "1",
           "PALACE_RESPONSE_BLOCK_SIZE": "48"}
RECORDS = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/"
               "readiness-20261029/stageF/chainport/records-274f/qualify-274fe1a150c0/"
               "fab-spatial-3-edge-a61f856dce13/spatial-3-edge-a61f856dce13/main")
REMOTE = ("/data/home/simlap/coupon_accuracy_assessment_20260913/e-remote/"
          "fab-spatial-3-edge-a61f856dce13/spatial-3-edge-a61f856dce13/main")
FAKE = HERE / "fake_paladin_cli.py"
FAKE_TOOLS = '''import json, os, sys
with open(os.environ["FAKE_TOOL_LOG"], "a") as f:
    f.write(json.dumps(sys.argv[1:]) + "\\n")
if sys.argv[1] == "check":
    print("FAIL (1)" if os.environ.get("FAKE_PIN_FAIL") == "check" else "PASS fake tip=0 objects=3 check=now")
    sys.exit(1 if os.environ.get("FAKE_PIN_FAIL") == "check" else 0)
if sys.argv[1] == "postrun":
    failed = os.environ.get("FAKE_PIN_FAIL") == "postrun"
    print("POSTRUN-V2 %s run=%s" % ("FAIL" if failed else "PASS", sys.argv[sys.argv.index("--run-id") + 1]))
    sys.exit(1 if failed else 0)
failed = os.environ.get("FAKE_SCHEMA_FAIL")
print(("INVALID %s: x" if failed else "VALID %s") % sys.argv[2])
'''
FAKE_ZERO_HOLE = '''import json, os, sys
json.dump({"Verdict": "CLEAN"}, open(sys.argv[sys.argv.index("--out") + 1], "w"))
sys.exit(int(os.environ.get("FAKE_ZERO_HOLE_RC", "0")))
'''
IMAGE_BINARY_SHA256 = "b" * 64


def fake_block(root):
    root = Path(root)
    (root / "tools.py").write_text(FAKE_TOOLS)
    (root / "zero_hole.py").write_text(FAKE_ZERO_HOLE)
    (root / "pin.json").write_text(json.dumps({"selector": "dev-fake-pin-0123456789"}))
    cli = root / "paladin"
    cli.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{FAKE}" "$@"\n')
    cli.chmod(0o755)
    return {"Selector": "dev-fake-pin-0123456789", "PinRecord": str(root / "pin.json"), "PinCheck": str(root / "tools.py"),
            "SchemaCheck": str(root / "tools.py"), "ToolPython": sys.executable, "ImageSHA256": "a" * 64,
            "ImageBinary": "/opt/view/bin/palace-arm64.bin", "ImageBinarySHA256": IMAGE_BINARY_SHA256, "CLI": str(cli),
            "Host": "local", "Lane": "chainport", "MemTotalGiB": {"m8g.48xlarge": 742.49, "r8g.48xlarge": 1450.3},
            "ZeroHoleScan": str(root / "zero_hole.py")}


def write_cookie(path, seconds_ahead):
    """A Netscape cookie file in ~/.midway/cookie's form: one session entry expiring `seconds_ahead` from now."""
    expiry = int(time.time() + seconds_ahead)
    Path(path).write_text("# Netscape HTTP Cookie File\n"
                          f"#HttpOnly_.midway.example\tTRUE\t/\tTRUE\t{expiry}\tsession\tVALUE-NOT-TO-BE-READ\n"
                          f".midway.example\tTRUE\t/\tTRUE\t{expiry + 10**6}\tother\tx\n")
    return str(path)


class RecordingWatchers:
    """A stand-in for paladin_backend.RunWatchers: records each watcher start (argv, directory), none runs."""

    def __init__(self):
        self.started, self.live = [], set()

    def alive(self, directory):
        return directory in self.live

    def start(self, argv, directory):
        self.started.append((argv, directory))
        self.live.add(directory)
        return 4242

    def runs(self):
        return [int(argv[argv.index("--run-id") + 1]) for argv, _ in self.started]


def three_stage_plan(root, *, nodes=None):
    """A worker, its reducer and a response-matrix stage; configurations under root/remote/main."""
    main = Path(root) / "remote" / "main"
    archive = str(main / "p4" / "archive")
    mesh = Path(root) / "local" / "mesh.msh"
    mesh.parent.mkdir(parents=True, exist_ok=True)
    mesh.write_text("$MeshFormat\n")
    stages = []
    for name, environment, requires, output in (
            ("w", {**WORKER, "PALACE_RESPONSE_ARCHIVE_DIR": archive}, [], "p4/worker"),
            ("r", {**REDUCER, "PALACE_RESPONSE_ARCHIVE_DIR": archive}, ["w"], "p4/reducer"),
            ("local-edge", {}, [], "local-edge/output")):
        path = main / (name + "-in.json")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({
            "Problem": {"Type": "Electrostatic", "Output": str(main / output)}, "Model": {"Mesh": "/data/case/mesh.msh"},
            "Boundaries": {"PrescribedPotential": [{"Index": i, "Attributes": [1]} for i in (1, 2, 3)]},
            "Solver": {"Electrostatic": {"ResponseMatrix": True}}}))
        stages.append({"Name": name, "Config": str(path), "Environment": environment, "CapSeconds": 900,
                       "Requires": requires})
    plan = {"Case": "case", "Ranks": 4, "Instance": {"Type": "m8g.48xlarge", "EstimatedPalacePeakGiB": 112.4},
            "BinarySHA256": IMAGE_BINARY_SHA256, "MeshRemote": "/data/case/mesh.msh", "MeshLocal": str(mesh),
            "MeshSHA256": hashlib.sha256(mesh.read_bytes()).hexdigest(), "Stages": stages}
    if nodes:
        plan["Nodes"] = nodes
    (main / "plan.json").write_text(json.dumps(plan))
    return main / "plan.json", archive


class ResponseArchiveKeys(unittest.TestCase):
    def test_the_stage_forms_of_build_plan(self):
        self.assertEqual(paladin_backend.response_archive_keys(WORKER),
                         {"Mode": "Archive", "Directory": ARCHIVE, "SourceTiming": True})
        self.assertEqual(paladin_backend.response_archive_keys(REDUCER),
                         {"Mode": "Reduce", "Directory": ARCHIVE, "BlockSize": 48})
        self.assertEqual(paladin_backend.response_archive_keys({}), {})
        self.assertEqual(paladin_backend.response_archive_keys(
            {"PALACE_RESPONSE_ESTIMATE_ONLY": "1", "PALACE_RESPONSE_ESTIMATE_REQUEST": "/data/r.json",
             "PALACE_RESPONSE_ARCHIVE_DIR": ARCHIVE, "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS": "1"}),
            {"Mode": "Estimate", "EstimateRequest": "/data/r.json", "Directory": ARCHIVE, "RecycleInitialGuess": True})

    def test_refusals(self):
        for environment in ({"PALACE_RESPONSE_ARCHIVE_DIR": "archive"},
                            {"PALACE_RESPONSE_ARCHIVE_ONLY": "1", "PALACE_RESPONSE_REDUCE_ONLY": "1"},
                            {"PALACE_RESPONSE_ARCHIVE_ONLY": "0"}, {"PALACE_RESPONSE_SOURCE_TIMING": "yes"},
                            {"PALACE_RESPONSE_BLOCK_SIZE": "0"}, {"PALACE_RESPONSE_DERIVED_CACHE": "0"},
                            {"OMP_NUM_THREADS": "1"}):
            with self.subTest(environment=environment), self.assertRaises(ValueError):
                paladin_backend.response_archive_keys(environment)

    def test_run_config(self):
        config = {"Problem": {"Type": "Electrostatic"}, "Solver": {"Electrostatic": {"ResponseMatrix": True}}}
        out = paladin_backend.run_config(config, WORKER)
        self.assertEqual(out["Solver"]["Electrostatic"],
                         {"ResponseMatrix": True,
                          "ResponseArchive": {"Mode": "Archive", "Directory": ARCHIVE, "SourceTiming": True}})
        self.assertNotIn("ResponseArchive", config["Solver"]["Electrostatic"])
        self.assertEqual(paladin_backend.run_config(config, {}), config)
        with self.assertRaises(ValueError):
            paladin_backend.run_config(out, REDUCER)


class BlockSizingLookup(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.block = fake_block(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_the_block_and_the_one_route_rule(self):
        self.assertIsNone(paladin_backend.backend_of("pbs", {"Launcher": {"Kind": "apptainer"}}))
        self.assertEqual(paladin_backend.backend_of("paladin", {"Paladin": self.block})["Selector"], self.block["Selector"])
        refused = [{"Paladin": {**self.block, "Extra": 1}},
                   {"Paladin": {key: value for key, value in self.block.items() if key != "ZeroHoleScan"}},
                   {"Paladin": self.block, "Launcher": {"Kind": "apptainer"}},
                   {"Paladin": {**self.block, "ImageSHA256": "x"}},
                   {"Paladin": {**self.block, "PinCheck": "relative.py"}},
                   {"Paladin": {**self.block, "MemTotalGiB": {}}},
                   {"Paladin": {**self.block, "Selector": "dev-other-pin-9999999999"}}]  # not the pin record's
        for profile in refused:
            with self.subTest(profile=sorted(profile)), self.assertRaises(ValueError):
                paladin_backend.backend_of("paladin", profile)
        with self.assertRaises(ValueError):
            paladin_backend.backend_of("pbs", {"Paladin": self.block})
        with self.assertRaises(ValueError):
            paladin_backend.backend_of("paladin", {})
        line = paladin_backend.route_line(self.block)
        self.assertEqual(len(line.splitlines()), 1)
        self.assertIn("route paladin: selector dev-fake-pin-0123456789, pin record " + self.block["PinRecord"], line)

    def test_e2_sizing(self):
        memory = self.block["MemTotalGiB"]
        chosen = paladin_backend.instance_for(112.4366594, memory)  # PBS 118804's thin plan
        self.assertEqual(chosen["Type"], "m8g.48xlarge")
        self.assertAlmostEqual(chosen["FigureGiB"], 146.17, places=2)
        self.assertAlmostEqual(chosen["BoundGiB"], 519.743, places=3)
        self.assertEqual(paladin_backend.instance_for(380.0, memory)["Type"], "r8g.48xlarge")  # 494 <= 520, margin 5 %
        self.assertEqual(paladin_backend.instance_for(420.0, memory)["Type"], "r8g.48xlarge")  # above m8g's bound
        with self.assertRaises(paladin_backend.PaladinError):
            paladin_backend.instance_for(800.0, memory)  # 1040 > 0.70 x 1450.3

    def test_lookup_decision(self):
        no_run = {"run_id": 1, "status": "submission_failed", "job_id": None}
        live = {"run_id": 2, "status": "queued", "job_id": "7"}
        self.assertEqual(paladin_backend.lookup_decision([]), ("submit", 0))
        self.assertEqual(paladin_backend.lookup_decision([no_run]), ("submit", 1))
        self.assertEqual(paladin_backend.lookup_decision([live, no_run]), ("adopt", live))
        for rows in ([no_run, dict(no_run, run_id=3)], [live, dict(live, run_id=4, status="completed")],
                     [dict(no_run, job_id="9"), live]):
            with self.subTest(rows=rows), self.assertRaises(paladin_backend.PaladinError):
                paladin_backend.lookup_decision(rows)
        self.assertEqual([paladin_backend.run_class(s) for s in ("running", "canceling", "completed", "failed", "canceled",
                                                                 "out_of_memory", "weird")],
                         ["live", "live", "completed", "failed", "cancelled", "out-of-memory", "unknown"])


class FakeCLICase(unittest.TestCase):
    """The fake cluster, tools and records of the drive tests (no test of its own)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.block = fake_block(self.root)
        self.state = self.root / "fake-state.json"
        self.argv_log = self.root / "fake-argv.jsonl"
        self.tool_log = self.root / "fake-tools.jsonl"
        (self.root / "runs").mkdir()
        self.out = self.root / "out"
        self.plan_path, self.archive = three_stage_plan(self.root)
        self.environment = mock.patch.dict(os.environ, {"FAKE_PALADIN_STATE": str(self.state),
                                                        "FAKE_PALADIN_LOG": str(self.argv_log),
                                                        "FAKE_TOOL_LOG": str(self.tool_log)})
        self.environment.start()
        self.slept = []  # the waits between repeated lookups (INTENT_RULE), recorded, not slept
        self.pause = mock.patch.object(paladin_backend, "pause", side_effect=self.slept.append)
        self.pause.start()
        self.watchers = RecordingWatchers()
        self.cookie = write_cookie(self.root / "cookie", 86400)
        self.script()

    def tearDown(self):
        self.pause.stop()
        self.environment.stop()
        self.tmp.cleanup()

    def script(self, plan=None, lookups=None):
        with open(self.state, "w") as f:
            json.dump({"Root": str(self.root / "runs"), "Plan": plan or {}, "Lookups": lookups or []}, f)

    def lookups(self, answers):
        """The next my-jobs answers of the fake (its "Lookups"), without touching its runs."""
        state = json.loads(self.state.read_text())
        state["Lookups"] = list(answers)
        self.state.write_text(json.dumps(state))

    def job(self, now=paladin_backend.utc_now):
        return paladin_backend.PaladinJob(str(self.plan_path), block=self.block, remote_dir=str(self.plan_path.parent),
                                          local_dir=str(self.out), archive_estimates_gb={self.archive: 0.001},
                                          log=lambda message: None, now=now, watchers=self.watchers,
                                          session_cookie=self.cookie)

    def drive(self, now=paladin_backend.utc_now):
        return paladin_backend.drive(self.job(now), poll_seconds=0, sleep=lambda seconds: None)

    def calls(self):
        if not self.argv_log.exists():
            return []
        return [json.loads(line) for line in self.argv_log.read_text().splitlines()]

    def submitted(self):
        return [call[call.index("--label") + 1] for call in self.calls() if call[:2] == ["palace", "submit-detailed"]]

    def records(self):
        return json.loads((self.out / "paladin-status.json").read_text())["Stages"]


class DriveWithFakeCLI(FakeCLICase):
    def test_runs_in_order_with_every_check_and_record(self):
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        submit = next(call for call in self.calls() if call[:2] == ["palace", "submit-detailed"])
        for flag, value in (("--instance-type", "m8g.48xlarge"), ("--nodes", "1"), ("--ranks-per-node", "4"),
                            ("--palace", "dev-fake-pin-0123456789")):
            self.assertEqual(submit[submit.index(flag) + 1], value)
        for token in ("programme=spr", "lane=chainport", "case=case", "stage=w"):
            self.assertIn(token, submit)
        self.assertTrue(any(part.startswith("plan=sha256:") for part in submit))
        records = state["Stages"]
        for name, record in records.items():
            self.assertEqual(record["State"], "complete")
            self.assertTrue(record["PinCheckBefore"].startswith("PASS "))
            self.assertTrue(record["PinCheckAfter"].startswith("POSTRUN-V2 PASS "))
            self.assertTrue(record["SchemaCheck"].startswith("VALID "))
            self.assertEqual(record["Instance"]["Type"], "m8g.48xlarge")
            self.assertAlmostEqual(record["NodeHours"], 600.0 / 3600.0)
            self.assertEqual(sorted(record["LeafDiff"]), [".Model.Mesh", ".Model.Refinement.SaveAdaptMesh",
                                                          ".Problem.Output", ".Solver.Device"])
            self.assertEqual(record["QdashUpload"]["Artifacts"][0], ["log", "FILE_UPLOAD"])
            self.assertEqual((record["Parsed"]["H1"], record["Parsed"]["PCG"]), (1234, [7]))
        self.assertEqual(records["w"]["ArchiveGuard"]["FilesBefore"], 0)
        self.assertEqual((records["w"]["Archive"]["Found"], records["w"]["Archive"]["OK"]), (12, True))
        self.assertNotIn("Archive", records["r"])
        reducer_outputs = self.plan_path.parent / "p4" / "reducer"
        self.assertTrue((reducer_outputs / "surface-response-matrix.csv").is_file())
        with open(self.out / "paladin" / "w.json") as f:
            config = json.load(f)
        self.assertEqual(config["Solver"]["Electrostatic"]["ResponseArchive"],
                         {"Mode": "Archive", "Directory": self.archive, "SourceTiming": True})
        self.assertEqual(config["Model"]["Mesh"], str(self.root / "local" / "mesh.msh"))
        status = json.loads((self.plan_path.parent / "status.json").read_text())
        self.assertEqual((status["Backend"], status["State"]), ("paladin", "complete"))
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["complete"] * 3)
        self.assertEqual(status["Preflight"]["ImageBinarySHA256"], IMAGE_BINARY_SHA256)
        self.assertEqual(status["Preflight"]["Selector"], "dev-fake-pin-0123456789")
        self.assertAlmostEqual(status["TotalSeconds"], 1800.0)
        count = len(self.calls())
        self.drive()
        self.assertEqual(len(self.calls()), count)  # a completed job: no paladin call at all

    def killed_after_the_submit(self):
        real = paladin_backend.paladin_call

        def call(cli, arguments, **options):
            data = real(cli, arguments, **options)
            if arguments[:2] == ["palace", "submit-detailed"]:
                raise KeyboardInterrupt("the driver is killed after the submit returned")
            return data

        with mock.patch.object(paladin_backend, "paladin_call", side_effect=call), self.assertRaises(KeyboardInterrupt):
            self.drive()
        self.assertNotIn("RunID", self.records().get("w", {}))

    def test_a_driver_killed_after_the_submit_adopts_the_run(self):
        self.killed_after_the_submit()
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["w"]["Adopted"]["Status"], "queued")
        self.assertEqual(state["Stages"]["w"]["State"], "complete")

    def test_a_run_that_completed_while_the_driver_was_dead_is_adopted(self):
        self.killed_after_the_submit()
        run_id = json.loads(self.state.read_text())["NextID"] - 1
        for _ in range(2):  # the run starts and ends while no driver watches
            paladin_backend.paladin_call(self.block["CLI"], ["palace", "describe-run", str(run_id), "--json"])
        state = self.drive()
        self.assertEqual(state["Stages"]["w"]["Adopted"]["Status"], "completed")
        self.assertEqual(self.submitted().count("chainport-w"), 1)
        self.assertEqual(state["Stages"]["local-edge"]["State"], "complete")

    def test_a_driver_killed_in_its_wait_confirms_the_recorded_run(self):
        def killed(seconds):
            raise KeyboardInterrupt("the driver is killed in its wait")

        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=killed)
        run_id = self.records()["w"]["RunID"]
        state = self.drive()
        self.assertEqual(state["Stages"]["w"]["RunID"], run_id)
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_a_row_failed_by_the_upload_alone_reads_completed_by_the_five_conditions(self):
        self.script({"chainport-r": "failed-upload"})
        state = self.drive()
        record = state["Stages"]["r"]
        self.assertEqual((record["State"], record["RowStatus"], record["CompletedByNote1"]), ("complete", "failed", True))
        self.assertTrue(all(record["CompletionRead"]["Conditions"].values()))
        self.assertIn("harness upload", record["CompletionRead"]["ErrorLines"][0])
        self.assertIn("ArrowTypeError", record["CompletionRead"]["TracebackLast"])

    def assert_stops_and_never_resubmits(self, plan, state_expected, *, now=paladin_backend.utc_now):
        self.script(plan)
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(now)
        self.assertEqual(self.records()["w"]["State"], state_expected)
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(now)  # the recorded State stops it
        (self.out / "paladin-status.json").unlink()  # a driver without the local record: the lookup stops it
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(now)
        self.assertEqual(self.submitted(), ["chainport-w"])
        self.assertNotIn("r", self.records())

    def test_a_failed_run_stops_and_is_never_resubmitted(self):
        self.assert_stops_and_never_resubmits({"chainport-w": "failed"}, "failed")
        record = self.records()["w"]
        self.assertFalse(record["CompletionRead"]["Conditions"]["1-ElapsedTimeReport"])
        self.assertNotIn("Outputs", record)  # no read of its data

    def test_a_copy_that_differs_from_its_source_stops(self):
        with mock.patch.object(paladin_backend, "manifest", side_effect=[{"palace.json": "1"}, {"palace.json": "2"}]), \
                self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "failed")
        self.assertNotIn("r", self.records())

    def test_a_run_cancelled_by_the_cap_stops_as_timed_out(self):
        later = lambda: paladin_backend.utc_now() + datetime.timedelta(seconds=10**4)
        self.assert_stops_and_never_resubmits({"chainport-w": "hang"}, "timed-out", now=later)
        cancels = [call for call in self.calls() if call[:2] == ["palace", "cancel"]]
        self.assertEqual(len(cancels), 1)

    def test_out_of_memory_stops(self):
        self.assert_stops_and_never_resubmits({"chainport-w": "oom"}, "out-of-memory")

    def test_one_resubmit_after_a_transport_failure(self):
        self.script({"chainport-w": "transport-once"})
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["w"]["NoRunRowsBeforeSubmit"], 1)
        self.assertEqual(len(state["Stages"]["w"]["SubmitErrors"]), 1)

    def test_a_second_transport_failure_stops(self):
        self.script({"chainport-w": "transport"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-w"])

    def test_a_refused_submit_stops_without_a_run(self):
        self.script({"chainport-w": "refuse"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), ["chainport-w"])
        self.assertNotIn("RunID", self.records().get("w", {}))

    def test_the_checks_before_a_submit_fail_closed(self):
        for variable, value in (("FAKE_PIN_FAIL", "check"), ("FAKE_SCHEMA_FAIL", "1")):
            with self.subTest(variable=variable), mock.patch.dict(os.environ, {variable: value}):
                with self.assertRaises(paladin_backend.PaladinError):
                    self.drive()
                self.assertEqual(self.submitted(), [])
        Path(self.archive).mkdir(parents=True, exist_ok=True)
        (Path(self.archive) / "source-000001-rank-000000-V.bin").write_bytes(b"old")
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), [])

    def test_a_postrun_pin_failure_stops_before_any_read(self):
        with mock.patch.dict(os.environ, {"FAKE_PIN_FAIL": "postrun"}), self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertFalse((self.plan_path.parent / "p4" / "worker").exists())  # nothing read or copied
        self.assertNotIn("r", self.records())

    def assert_worker_fails_after_its_run(self, plan, environment=None):
        self.script(plan)
        with mock.patch.dict(os.environ, environment or {}), self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "failed")
        self.assertNotIn("r", self.records())
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_short_archive_stops(self):
        self.assert_worker_fails_after_its_run({"chainport-w": "short-archive"})
        self.assertEqual(self.records()["w"]["Archive"]["MissingCount"], 1)

    def test_an_archive_with_a_zero_hole_stops(self):
        self.assert_worker_fails_after_its_run({}, {"FAKE_ZERO_HOLE_RC": "1"})

    def test_a_configuration_paladin_changed_beyond_its_own_paths_stops(self):
        self.assert_worker_fails_after_its_run({"chainport-w": "drift"})
        self.assertIn(".Solver.Order", self.records()["w"]["LeafDiff"])

    def test_two_runs_of_one_stage_stop(self):
        self.killed_after_the_submit()
        state = json.loads(self.state.read_text())
        twin = dict(state["Runs"][str(state["NextID"] - 1)])
        state["Runs"][str(state["NextID"])] = twin
        self.state.write_text(json.dumps(state))
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), ["chainport-w"])

    # Decision 1820 (1), MAJOR-1 of review r1799: the SubmitIntent, failed lookups, two consistent lookups.
    def test_the_intent_is_on_disk_before_the_submit_and_its_nonce_is_in_the_metadata(self):
        self.killed_after_the_submit()  # the driver ended between the submit and the save of its RunID
        intent = self.records()["w"]["SubmitIntents"][-1]
        self.assertNotIn("Outcome", intent)
        self.assertEqual((intent["Stage"], len(intent["Nonce"])), ("w", 32))
        submit = next(call for call in self.calls() if call[:2] == ["palace", "submit-detailed"])
        self.assertIn(f"nonce={intent['Nonce']}", submit)

    def test_r1_an_open_intent_and_one_empty_answer_adopts_without_a_second_submit(self):
        self.killed_after_the_submit()
        self.lookups(["empty"])
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        intent = state["Stages"]["w"]["SubmitIntents"][0]
        self.assertEqual((intent["Outcome"].split()[0], intent["NonceMatches"]), ("adopted", True))
        self.assertEqual(state["Stages"]["w"]["State"], "complete")

    def test_r2_r3_a_partial_or_keyless_answer_is_a_failed_lookup_and_never_a_submit(self):
        for answer in ("omitted", "nokey", "fail"):
            with self.subTest(answer=answer):
                self.setUp_fresh()
                self.killed_after_the_submit()
                self.lookups([answer])
                with self.assertRaises(paladin_backend.PaladinError) as stop:
                    self.drive()
                if answer != "fail":
                    self.assertIn("FAILED lookup", str(stop.exception))
                self.assertEqual(self.submitted(), ["chainport-w"])
                self.drive()  # the next driver, with a true listing, adopts the one run
                self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def setUp_fresh(self):
        """A new fake cluster and an empty local record in the same test."""
        for path in (self.state, self.argv_log):
            if path.exists():
                path.unlink()
        for directory in (self.out, Path(self.archive), self.root / "runs"):
            shutil.rmtree(directory, ignore_errors=True)
        (self.root / "runs").mkdir()
        self.script()

    def test_an_open_intent_with_no_run_found_stops_after_three_lookups_over_60_s(self):
        self.killed_after_the_submit()
        del self.slept[:]
        self.lookups(["empty"] * 4)
        with self.assertRaises(paladin_backend.PaladinError) as stop:
            self.drive()
        self.assertIn("submit intent without an outcome", str(stop.exception))
        self.assertEqual(self.slept, [30.0, 30.0])
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_lost_record_and_one_empty_answer_adopts_by_the_second_lookup(self):
        self.killed_after_the_submit()
        (self.out / "paladin-status.json").unlink()  # the lost record: no intent, no RunID
        self.lookups(["empty"])
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["w"]["Adopted"]["Status"], "queued")
        self.assertIn(20.0, self.slept)

    def test_a_first_submit_needs_two_consistent_empty_lookups_20_s_apart(self):
        self.drive()
        listings = [call for call in self.calls() if call[:2] == ["palace", "my-jobs"]]
        self.assertEqual(self.slept, [20.0] * 3)  # one gap before each of the three submits
        self.assertEqual(len(listings), 9)  # per stage: two before the submit, one after it

    def test_two_lookups_that_disagree_stop_without_a_submit(self):
        self.script({"chainport-w": "transport-once"})
        real = paladin_backend.find_runs
        answers = iter([[], [{"run_id": 1, "status": "submission_failed", "job_id": None}]])
        with mock.patch.object(paladin_backend, "find_runs", side_effect=lambda cli, identity: next(answers)), \
                self.assertRaises(paladin_backend.PaladinError) as stop:
            self.drive()
        self.assertIn("disagree", str(stop.exception))
        self.assertEqual(self.submitted(), [])
        self.assertIs(paladin_backend.find_runs, real)

    def test_two_runs_of_one_stage_end_the_newer_early_and_mark_the_archive_unusable(self):
        self.killed_after_the_submit()
        state = json.loads(self.state.read_text())
        first = state["NextID"] - 1
        twin = dict(state["Runs"][str(first)], run_dir=str(self.root / "runs" / "twin"))
        state["Runs"][str(first + 1)], state["NextID"] = twin, first + 2
        self.state.write_text(json.dumps(state))
        with self.assertRaises(paladin_backend.DuplicateRuns):
            self.drive()
        record = self.records()["w"]
        duplicate = record["DuplicateRuns"][-1]
        self.assertEqual((duplicate["Kept"], duplicate["Newer"][0]["RunID"], duplicate["Newer"][0]["EarlyEnded"]),
                         (first, first + 1, True))
        cancels = [call for call in self.calls() if call[:2] == ["palace", "cancel"]]
        self.assertEqual(cancels, [["palace", "cancel", str(first + 1), "--json"]])
        runs = json.loads(self.state.read_text())["Runs"]
        self.assertEqual((runs[str(first)]["status"], runs[str(first + 1)]["status"]), ("queued", "canceled"))
        marker = Path(self.archive) / paladin_backend.UNUSABLE_MARKER
        self.assertTrue(marker.is_file())
        self.assertEqual(self.submitted(), ["chainport-w"])
        job = self.job()
        with self.assertRaises(paladin_backend.PaladinError):
            job.refuse_unusable_archive("r", self.archive)  # a reducer refuses it
        with self.assertRaises(paladin_backend.PaladinError):
            job.archive_guard(self.archive)  # a writer refuses it (not empty)

    # Decision 1820 (2), MAJOR-2 of review r1799: one driver per root (LOCK_RULE).
    def test_r8_a_second_driver_of_one_plan_is_refused_while_the_first_lives(self):
        real = paladin_backend.find_runs
        outcomes = {}

        def second_driver():
            job = paladin_backend.PaladinJob(str(self.plan_path), block=self.block, remote_dir=str(self.plan_path.parent),
                                             local_dir=str(self.root / "out-B"), archive_estimates_gb={self.archive: 0.001},
                                             log=lambda message: None)
            try:
                paladin_backend.drive(job, poll_seconds=0, sleep=lambda seconds: None)
                outcomes["B"] = "drove"
            except paladin_backend.PaladinError as error:
                outcomes["B"] = str(error)

        def find(cli, identity):
            if "B" not in outcomes:  # the first driver's first lookup: the second driver starts now, in its own thread
                thread = threading.Thread(target=second_driver)
                thread.start()
                thread.join(60)
            return real(cli, identity)

        with mock.patch.object(paladin_backend, "find_runs", side_effect=find):
            state = self.drive()
        self.assertIn("held by a live driver", outcomes["B"])
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["local-edge"]["State"], "complete")
        self.assertFalse((self.root / "out-B" / "paladin-status.json").exists())

    def test_a_lock_held_by_another_live_process_is_refused(self):
        lock = self.plan_path.parent / "paladin-driver.lock"
        holder = subprocess.Popen([sys.executable, "-c", "import fcntl, os, sys; d = os.open(sys.argv[1], os.O_RDWR | os.O_CREAT); "
                                   "fcntl.flock(d, fcntl.LOCK_EX); print('held', flush=True); sys.stdin.read()", str(lock)],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        try:
            self.assertEqual(holder.stdout.readline().strip(), "held")
            with self.assertRaises(paladin_backend.PaladinError) as stop:
                self.drive()
            self.assertIn("held by a live driver", str(stop.exception))
            self.assertEqual(self.calls(), [])  # no paladin call at all
        finally:
            holder.stdin.close()
            holder.wait(30)
            holder.stdout.close()
        self.drive()  # the holder ended: the lock is free
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_a_stale_lock_is_taken_over_with_a_recorded_line(self):
        ended = subprocess.run([sys.executable, "-c", "import os; print(os.getpid())"], capture_output=True, text=True)
        stale = json.dumps({"PID": int(ended.stdout), "Host": "this host", "StartUTC": "2026-11-15T00:00:00Z", "Argv": []})
        self.out.mkdir(parents=True, exist_ok=True)
        (self.out / "paladin-driver.lock").write_text(stale)  # a record under a free lock: its driver ended
        self.drive()
        lines = (self.out / "paladin-driver.lock.takeovers").read_text().splitlines()
        self.assertEqual(len(lines), 1)
        self.assertIn("took over the stale lock of " + stale, lines[0])
        self.assertEqual((self.out / "paladin-driver.lock").read_text(), "")  # released at the end of the drive
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    # Decision 1820 (3), MAJOR-3 of review r1799: a transport answer during the wait is not a state (WAIT_RULE).
    def describes(self, answers):
        state = json.loads(self.state.read_text())
        state["Describes"] = list(answers)
        self.state.write_text(json.dumps(state))

    def test_r4_one_empty_answer_while_running_is_re_asked_and_the_run_completes(self):
        for answer in ("empty", "nostatus", "otherid", "fail"):
            with self.subTest(answer=answer):
                self.setUp_fresh()
                self.describes(["true", answer])  # the run's second describe-run (running -> completed) answers badly
                del self.slept[:]
                state = self.drive()
                self.assertEqual(state["Stages"]["w"]["State"], "complete")
                self.assertNotIn("WaitTransport", state["Stages"]["w"])
                self.assertIn(10.0, self.slept)  # the back-off
                self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_r5_an_explicit_unlisted_status_is_unknown(self):
        self.describes(["unlisted"])
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "unknown")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_r6_a_completed_answer_without_its_run_directory_is_re_asked(self):
        self.describes(["true", "truncated"])
        state = self.drive()
        self.assertEqual(state["Stages"]["w"]["State"], "complete")
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_ten_failed_polls_over_600_s_stop_with_no_state_and_a_resume_continues(self):
        clock = [paladin_backend.utc_now()]

        def now():
            clock[0] += datetime.timedelta(seconds=20)
            return clock[0]

        self.describes(["true"] + ["empty"] * 60)
        with self.assertRaises(paladin_backend.PaladinError) as stop:
            paladin_backend.drive(self.job(now), poll_seconds=0, sleep=lambda seconds: None)
        self.assertIn("no State recorded", str(stop.exception))
        record = self.records()["w"]
        self.assertNotIn("State", record)
        self.assertGreaterEqual(record["WaitTransport"]["Count"], paladin_backend.WAIT_TRANSPORT_POLLS)
        self.assertEqual([call for call in self.calls() if call[:2] == ["palace", "cancel"]], [])
        self.describes([])
        state = self.drive()  # the resume is not refused: it confirms the run and reads it
        self.assertEqual(state["Stages"]["w"]["State"], "complete")
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_a_running_answer_without_started_counts_the_cap_from_the_submit_time(self):
        clock = [paladin_backend.utc_now()]

        def now():  # 1000 s per reading: the cap (900 s) is past at the first poll after the submit
            clock[0] += datetime.timedelta(seconds=1000)
            return clock[0]

        self.script({"chainport-w": "hang"})
        self.describes(["truncated"])  # the first describe-run (queued -> running) answers without 'started'
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(now)
        record = self.records()["w"]
        self.assertTrue(record["CapCancel"]["From"].startswith("submitted"))
        self.assertEqual(record["State"], "timed-out")
        self.assertEqual(len([call for call in self.calls() if call[:2] == ["palace", "cancel"]]), 1)

    # Decision 1820 (4), MAJOR-4 of review r1799: every run has an unattended run watcher (WATCHER_RULE).
    def live_runs_without_a_watcher(self):
        runs = json.loads(self.state.read_text())["Runs"]
        live = [int(run_id) for run_id, run in runs.items() if paladin_backend.run_class(run["status"]) == "live"]
        return sorted(set(live) - set(self.watchers.runs()))

    def test_every_submitted_run_gets_a_watcher_with_its_stage_cap(self):
        state = self.drive()
        self.assertEqual(self.watchers.runs(), [state["Stages"][name]["RunID"] for name in ("w", "r", "local-edge")])
        argv, directory = self.watchers.started[0]
        for flag, value in (("--cap-seconds", "900.0"), ("--host", "local"), ("--bar", "0.85"), ("--grace", "120.0"),
                            ("--job", state["Stages"]["w"]["JobID"])):
            self.assertEqual(argv[argv.index(flag) + 1], value)
        self.assertEqual(directory, str(self.out / "watchers" / f"w-{state['Stages']['w']['RunID']}"))

    def test_a_resume_re_attaches_to_a_live_watcher_and_restarts_a_missing_one(self):
        def ended(seconds):
            raise KeyboardInterrupt("the driver process ends in its wait")

        self.script({"chainport-w": "hang"})
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        self.assertEqual(len(self.watchers.started), 1)
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        self.assertEqual(len(self.watchers.started), 1)  # the resume re-attached: the live watcher is kept
        self.watchers.live.clear()  # the watcher process is gone (a host reboot)
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        self.assertEqual([directory for _, directory in self.watchers.started], [self.watchers.started[0][1]] * 2)
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_every_stop_path_leaves_each_live_run_with_a_live_watcher(self):
        def ended(seconds):
            raise KeyboardInterrupt("the driver process is ended (an RSS stop of its orchestrator)")

        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        clock = [paladin_backend.utc_now()]

        def now():
            clock[0] += datetime.timedelta(seconds=20)
            return clock[0]

        self.describes(["empty"] * 60)  # the MAJOR-3 STOP (failed polls) with the run live
        with self.assertRaises(paladin_backend.PaladinError):
            paladin_backend.drive(self.job(now), poll_seconds=0, sleep=lambda seconds: None)
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        self.describes([])
        self.setUp_fresh()
        self.watchers = RecordingWatchers()
        self.killed_after_the_submit()  # the window of MAJOR-1: the run is live, its RunID never saved
        self.assertEqual(self.live_runs_without_a_watcher(), [int(json.loads(self.state.read_text())["NextID"]) - 1])
        state = json.loads(self.state.read_text())
        first = state["NextID"] - 1
        state["Runs"][str(first + 1)] = dict(state["Runs"][str(first)], run_dir=str(self.root / "runs" / "twin"),
                                             job_id="twin-job")
        state["NextID"] = first + 2
        self.state.write_text(json.dumps(state))
        with mock.patch.dict(os.environ, {"FAKE_CANCEL_FAIL": "1"}), self.assertRaises(paladin_backend.DuplicateRuns):
            self.drive()  # two runs: the newer one's early end fails
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        caps = {int(argv[argv.index("--run-id") + 1]): argv[argv.index("--cap-seconds") + 1] for argv, _ in self.watchers.started}
        self.assertEqual(caps, {first: "900.0", first + 1: "0.0"})

    def test_a_run_whose_watcher_does_not_start_is_ended_early(self):
        self.watchers.start = lambda argv, directory: None
        with self.assertRaises(paladin_backend.WatcherNotStarted):
            self.drive()
        record = self.records()["w"]
        self.assertTrue(record["Watchers"][-1]["EarlyEnded"])
        runs = json.loads(self.state.read_text())["Runs"]
        self.assertEqual([run["status"] for run in runs.values()], ["canceled"])
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_cap_end_by_the_watcher_is_read_as_timed_out(self):
        def ended(seconds):
            raise KeyboardInterrupt("the driver process ends in its wait")

        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        record = self.records()["w"]
        Path(record["Watchers"][-1]["Dir"]).mkdir(parents=True, exist_ok=True)
        (Path(record["Watchers"][-1]["Dir"]) / "end.json").write_text(json.dumps({"Kind": "cap"}))
        paladin_backend.paladin_call(self.block["CLI"], ["palace", "cancel", str(record["RunID"]), "--json"])
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "timed-out")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_stage_without_a_cap_is_not_submitted(self):
        plan = json.loads(self.plan_path.read_text())
        plan["Stages"][0]["CapSeconds"] = None
        self.plan_path.write_text(json.dumps(plan))
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), [])

    # Decision 1822 (b): the session check before every submit (SESSION_RULE).
    def test_a_session_that_expires_between_stages_stops_before_the_next_submit_and_resumes(self):
        real = paladin_backend.PaladinJob.finish

        def finish(job, stage, record, data):
            real(job, stage, record, data)
            if stage["Name"] == "w":
                write_cookie(self.cookie, 600)  # the session now expires in 10 min: < 900 + 3600 s

        with mock.patch.object(paladin_backend.PaladinJob, "finish", finish), \
                self.assertRaises(paladin_backend.SessionNotValid) as stop:
            self.drive()
        self.assertIn("SHORT", str(stop.exception))
        records = self.records()
        self.assertEqual(records["w"]["State"], "complete")
        self.assertNotIn("SubmitIntents", records["r"])  # stopped BEFORE the submit: resumable
        self.assertNotIn("State", records["r"])
        self.assertEqual(self.submitted(), ["chainport-w"])
        self.assertNotIn("VALUE-NOT-TO-BE-READ", (self.out / "paladin-status.json").read_text())
        write_cookie(self.cookie, 86400)  # mwinit
        state = self.drive()
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["r"]["SessionChecks"][-1]["Verdict"], "OK")
        self.assertEqual(state["Stages"]["r"]["SessionChecks"][-1]["NeedSeconds"], 4500.0)

    def test_no_session_or_an_expired_one_refuses_the_first_submit(self):
        for content in ("", None, -60):
            with self.subTest(content=content):
                if content is None:
                    os.unlink(self.cookie)
                elif content == "":
                    Path(self.cookie).write_text("# Netscape HTTP Cookie File\n")
                else:
                    write_cookie(self.cookie, content)
                with self.assertRaises(paladin_backend.SessionNotValid):
                    self.drive()
                self.assertEqual(self.submitted(), [])
        check = paladin_backend.session_check(write_cookie(self.cookie, 7200), 1800)
        self.assertEqual((check["Verdict"], check["Entries"]), ("OK", 1))  # the 'other' entry is not a session entry

    # Decision 1820 m-1: the completion read = PREREG-ESSURR-NOTE-1 exactly.
    def test_an_error_line_that_is_not_the_upload_harness_s_fails_condition_3(self):
        self.script({"chainport-w": "failed-other-error"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        record = self.records()["w"]
        self.assertFalse(record["CompletionRead"]["Conditions"]["3-OnlyUploadErrors"])
        self.assertEqual(record["State"], "failed")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_missing_table_of_a_full_solve_fails_condition_4(self):
        self.script({"chainport-local-edge": "missing-table"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        record = self.records()["local-edge"]
        self.assertEqual(record["CompletionRead"]["MissingOutputs"], ["error-indicators.csv"])
        self.assertEqual(record["State"], "failed")

    def test_a_run_its_watcher_ended_fails_condition_5(self):
        self.script({"chainport-w": "failed-upload"})
        real = paladin_backend.PaladinJob.watcher_end
        with mock.patch.object(paladin_backend.PaladinJob, "watcher_end",
                               lambda job, record: {"Kind": "memory"} if record["Name"] == "w" else real(job, record)), \
                self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertFalse(self.records()["w"]["CompletionRead"]["Conditions"]["5-NoStop"])

    def test_the_expected_outputs_per_stage_type(self):
        main = self.plan_path.parent
        configs = {name: paladin_backend.run_config(json.loads((main / f"{name}-in.json").read_text()), environment)
                   for name, environment in (("w", WORKER), ("r", REDUCER), ("local-edge", {}))}
        self.assertEqual(paladin_backend.expected_outputs(configs["w"]), ["palace.json", "palace_config_resolved.json"])
        self.assertEqual(paladin_backend.expected_outputs(configs["r"]),
                         ["domain-response-matrix.csv", "palace.json", "palace_config_resolved.json", "surface-response-matrix.csv"])
        self.assertIn("error-indicators.csv", paladin_backend.expected_outputs(configs["local-edge"]))
        estimate = paladin_backend.run_config(json.loads((main / "w-in.json").read_text()),
                                              {"PALACE_RESPONSE_ESTIMATE_ONLY": "1", "PALACE_RESPONSE_ARCHIVE_DIR": ARCHIVE})
        with self.assertRaises(paladin_backend.PaladinError):
            paladin_backend.expected_outputs(estimate)  # no list: not submitted

    REHEARSAL = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/readiness-20261029/stageF/chainport/"
                     "yb-rehearsal-1/driver/thin-spatial-3-edge-a61f856dce13/spatial-3-edge-a61f856dce13-thin")

    @unittest.skipUnless(REHEARSAL.is_dir(), "the Y-b rehearsal's records (2026-11-15) are absent on this host")
    def test_the_expected_outputs_equal_the_rehearsal_s_real_rows(self):
        """The three stage types of the rehearsal (paladin runs 14813302 / 14813397 / 14813418): the names the rule
        expects = the names each real row wrote (its fetched postpro copy)."""
        for stage, output in (("p4-worker", "p4/worker"), ("p4-reducer", "p4/reducer"), ("p4-local-edge", "p4-local-edge/output")):
            with self.subTest(stage=stage):
                prefix = "spatial-3-edge-a61f856dce13-thin-"
                config = json.loads((self.REHEARSAL / "main" / "paladin" / f"{prefix}{stage}.json").read_text())
                directory = self.REHEARSAL / "results" / "main" / (prefix + output.split("/")[0]) / "/".join(output.split("/")[1:])
                written = sorted(path.name for path in directory.iterdir() if path.is_file())
                self.assertEqual(paladin_backend.expected_outputs(config), written)

    def test_a_multi_node_plan_is_refused(self):
        plan_path, _ = three_stage_plan(self.root / "multi", nodes=2)
        with self.assertRaises(paladin_backend.PaladinError):
            paladin_backend.PaladinJob(str(plan_path), block=self.block, remote_dir=str(plan_path.parent),
                                       local_dir=str(self.out))

    def test_main_dry_run_and_submit(self):
        profile = self.root / "profile.json"
        profile.write_text(json.dumps({"Paladin": self.block}))
        common = [str(self.plan_path), "--profile", str(profile), "--out", str(self.out)]
        with mock.patch("sys.stdout"), mock.patch("sys.stderr"), \
                mock.patch.object(paladin_backend, "RunWatchers", RecordingWatchers), \
                mock.patch.object(paladin_backend, "MIDWAY_COOKIE", self.cookie):
            self.assertEqual(paladin_backend.main(common + ["--dry-run"]), 0)
            self.assertEqual(self.calls(), [])
            self.script({"chainport-local-edge": "failed"})
            submit = ["--submit", "--remote-dir", str(self.plan_path.parent), "--archive-gb", f"{self.archive}=0.001",
                      "--poll-seconds", "0"]
            for refused in ([], ["--lane-root", "relative/root"], ["--lane-root", str(self.root / "elsewhere")]):
                with self.subTest(refused=refused), self.assertRaises(SystemExit):
                    paladin_backend.main(common + submit + refused)  # m-3: no lane root / the archive outside it
            self.assertEqual(self.calls(), [])
            self.assertEqual(paladin_backend.main(common + submit + ["--lane-root", str(self.root)]), 3)
        states = {name: record["State"] for name, record in self.records().items()}
        self.assertEqual(states, {"w": "complete", "r": "complete", "local-edge": "failed"})
        self.assertTrue(all(record["DataGuard"]["FreeBytes"] > 0 for record in self.records().values()))

    def test_m2_the_data_guard_runs_before_every_submit(self):
        with mock.patch.object(paladin_backend, "ARCHIVE_FREE_BYTES", 1e18), self.assertRaises(paladin_backend.PaladinError) as stop:
            self.drive()
        self.assertIn("/data guard refuses", str(stop.exception))
        self.assertEqual(self.submitted(), [])
        real = paladin_backend.PaladinJob.data_guard
        guarded = []

        def data_guard(job):
            guarded.append(job.remote_dir)
            return real(job)

        with mock.patch.object(paladin_backend.PaladinJob, "data_guard", data_guard):
            self.drive()
        self.assertEqual(len(guarded), 3)  # one per submit, the non-writers included

    def test_m3_an_archive_outside_the_lane_root_is_refused_before_its_submit(self):
        job = paladin_backend.PaladinJob(str(self.plan_path), block=self.block, remote_dir=str(self.plan_path.parent),
                                         local_dir=str(self.out), archive_estimates_gb={self.archive: 0.001},
                                         log=lambda message: None, watchers=self.watchers, session_cookie=self.cookie,
                                         lane_root=str(self.root / "elsewhere"))
        with self.assertRaises(paladin_backend.PaladinError) as stop:
            paladin_backend.drive(job, poll_seconds=0, sleep=lambda seconds: None)
        self.assertIn("outside the lane root", str(stop.exception))
        self.assertEqual(self.submitted(), [])
        self.assertTrue(paladin_backend.under("/data/a/b", "/data/a") and not paladin_backend.under("/data/ab", "/data/a"))
        self.assertFalse(paladin_backend.under("/data/a/../b", "/data/a"))


class RunWatcherProcesses(FakeCLICase):
    """Decision 1820 (4): REAL paladin_run_watcher.py processes (nohup, their own session) against the fake
    CLI and a fake PBS (qstat / qdel) and exec-host ssh on PATH; the watcher's interval 1 s, grace 5 s,
    the worker's cap 3 s."""

    def setUp(self):
        super().setUp()
        fake_bin = self.root / "fake-bin"
        fake_bin.mkdir()
        for name, mode in (("qstat", "--pbs qstat"), ("qdel", "--pbs qdel"), ("ssh", "--ssh")):
            (fake_bin / name).write_text(f'#!/bin/sh\nexec "{sys.executable}" "{FAKE}" {mode} "$@"\n')
            (fake_bin / name).chmod(0o755)
        self.fake_bin = fake_bin
        self.path = mock.patch.dict(os.environ, {"PATH": f"{fake_bin}:{os.environ['PATH']}"})
        self.path.start()
        self.timing = [mock.patch.object(paladin_backend, "WATCHER_INTERVAL_SECONDS", 1.0),
                       mock.patch.object(paladin_backend, "WATCHER_GRACE_SECONDS", 5.0)]
        for patch in self.timing:
            patch.start()
        plan = json.loads(self.plan_path.read_text())
        plan["Stages"][0]["CapSeconds"] = 3
        self.plan_path.write_text(json.dumps(plan))
        self.watchers = paladin_backend.RunWatchers()

    def tearDown(self):
        for patch in self.timing:
            patch.stop()
        self.path.stop()
        for directory in (self.out / "watchers").glob("*"):
            self.wait_for(lambda: not self.watchers.alive(str(directory)), 30)
        for process in self.watchers.processes:
            process.wait(30)
        super().tearDown()

    def job(self, now=paladin_backend.utc_now):
        return paladin_backend.PaladinJob(str(self.plan_path), block=self.block, remote_dir=str(self.plan_path.parent),
                                          local_dir=str(self.out), archive_estimates_gb={self.archive: 0.001},
                                          log=lambda message: None, now=now, watchers=self.watchers,
                                          pbs_bin=str(self.fake_bin), session_cookie=self.cookie)

    @staticmethod
    def wait_for(condition, seconds):
        deadline = time.time() + seconds
        while time.time() < deadline:
            if condition():
                return True
            time.sleep(0.2)
        return condition()

    def ended_in_its_wait(self):
        """The driver submits the worker (plan 'hang'), polls once (the run starts) and its process ends."""
        def ended(seconds):
            raise KeyboardInterrupt("the driver process ends")

        self.script({"chainport-w": "hang"})
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        record = self.records()["w"]
        directory = Path(record["Watchers"][-1]["Dir"])
        self.assertTrue(self.watchers.alive(str(directory)))  # it outlives the driver
        self.assertTrue(self.wait_for(lambda: (directory / "end.json").is_file(), 30))
        self.assertTrue(self.wait_for(lambda: not self.watchers.alive(str(directory)), 30))
        return record, json.loads((directory / "end.json").read_text()), (directory / "events.log").read_text()

    def test_the_watcher_ends_the_run_at_its_cap_after_the_driver_ended(self):
        record, end, events = self.ended_in_its_wait()
        self.assertEqual((end["Kind"], end["Cancel"]["OK"], end["TerminalAfterCancel"]), ("cap", True, True))
        self.assertIn("from the PBS start", end["Reason"])
        self.assertIn("CAP ", events)
        self.assertIn("END exit 11", events)
        self.assertEqual(json.loads(self.state.read_text())["Runs"][str(record["RunID"])]["status"], "canceled")
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()  # the resumed driver classes it as a cap stop and submits nothing
        self.assertEqual(self.records()["w"]["State"], "timed-out")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_the_watcher_ends_the_run_at_the_memory_bar(self):
        with mock.patch.dict(os.environ, {"FAKE_MEMINFO": "1000 100"}):  # 0.90 used >= 0.85
            plan = json.loads(self.plan_path.read_text())
            plan["Stages"][0]["CapSeconds"] = 3600
            self.plan_path.write_text(json.dumps(plan))
            record, end, events = self.ended_in_its_wait()
        self.assertEqual(end["Kind"], "memory")
        self.assertIn("used 0.9000 of MemTotal >= bar 0.85", end["Reason"])
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "memory-stop")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_a_failed_cancel_falls_back_to_qdel_by_the_job_id(self):
        with mock.patch.dict(os.environ, {"FAKE_CANCEL_FAIL": "1"}):
            record, end, events = self.ended_in_its_wait()
        self.assertFalse(end["Cancel"]["OK"])
        self.assertEqual(end["Fallback"]["Command"], f"{self.fake_bin}/qdel {record['JobID']}")
        self.assertEqual((end["Fallback"]["ReturnCode"], end["Fallback"]["JobEnded"]), (0, True))
        self.assertIn("FALLBACK", events)
        self.assertEqual(json.loads(self.state.read_text())["Runs"][str(record["RunID"])]["status"], "canceled")

    def test_a_watcher_ends_with_its_run(self):
        plan = json.loads(self.plan_path.read_text())
        plan["Stages"][0]["CapSeconds"] = 3600
        self.plan_path.write_text(json.dumps(plan))
        state = self.drive()
        for name in ("w", "r", "local-edge"):
            directory = Path(state["Stages"][name]["Watchers"][-1]["Dir"])
            self.assertTrue(self.wait_for(lambda: not self.watchers.alive(str(directory)), 30))
            self.assertFalse((directory / "end.json").exists())
            self.assertRegex((directory / "events.log").read_text(), r"END exit 0: the (row is terminal|PBS job has ended)")
        self.assertEqual(state["Stages"]["local-edge"]["State"], "complete")


@unittest.skipUnless((RECORDS / "plan.json").exists(), "records of 274fe1a150c0 absent")
class PlanOfRecord(unittest.TestCase):
    def test_fab_plan_of_274fe1a150c0(self):
        with tempfile.TemporaryDirectory() as out:
            block = fake_block(out)
            job = paladin_backend.PaladinJob(str(RECORDS / "plan.json"), block={**block, "ImageBinarySHA256": "c" * 64},
                                             remote_dir=REMOTE, local_dir=os.path.join(out, "job"),
                                             mapping=[(REMOTE, str(RECORDS))])
            runs = job.dry_run()
        self.assertEqual(len(runs), 7)
        modes = [run["ResponseArchive"].get("Mode") for run in runs]
        self.assertEqual(modes, ["Archive", "Reduce"] * 3 + [None])
        for worker, reducer in zip(runs[0:6:2], runs[1:6:2]):
            self.assertEqual(worker["ResponseArchive"]["Directory"], reducer["ResponseArchive"]["Directory"])
            self.assertEqual(reducer["Requires"], [worker["Stage"]])
            self.assertTrue(worker["ResponseArchive"]["SourceTiming"])
            self.assertEqual((worker["Nodes"], worker["RanksPerNode"], worker["Instance"]["Type"]), (1, 192, "m8g.48xlarge"))


if __name__ == "__main__":
    unittest.main()
