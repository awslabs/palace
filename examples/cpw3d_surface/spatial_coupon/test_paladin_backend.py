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
import types
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


START_LIMIT = 5  # decision 1830 (2): a test that requests a 5th watcher start FAILS there (no unbounded starts)
REAL_RUN_WATCHERS = paladin_backend.RunWatchers  # the default starter (real processes), kept for the one opt-in class
REAL_SESSION_DIRECTORY = os.path.realpath(os.path.expanduser("~/.midway"))


def under_the_real_session_directory(path):
    path = os.path.realpath(os.path.expanduser(str(path)))
    return path == REAL_SESSION_DIRECTORY or path.startswith(REAL_SESSION_DIRECTORY + os.sep)


class GuardedTestCase(unittest.TestCase):
    """Decision 1836 m-7 / ruling 1837: the base of every test class here and in test_qualify_paladin.  A test FAILS when
    it reaches the DEFAULT watcher starter (a job built without an explicit starter that starts a watcher), a real
    watcher start (unless its class sets REAL_WATCHER_PROCESSES: RunWatcherProcesses only), or a session file under the
    real home's ~/.midway.  Each violation is recorded and raised; a cleanup fails the test on any recorded one (also
    when the code under test swallowed the raise)."""

    REAL_WATCHER_PROCESSES = False

    def setUp(self):
        super().setUp()
        self.guard_violations = []
        guard = self

        def violation(what):
            guard.guard_violations.append(what)
            raise AssertionError(f"test guard (decision 1836 m-7): {what}")

        class DefaultStarterReached:
            def __init__(self, *args, **options):
                violation("the DEFAULT watcher starter was reached (a job built without an explicit starter)")

        def real_start(starter, argv, directory, **options):
            violation(f"a REAL watcher start in a class without REAL_WATCHER_PROCESSES: {directory}")

        real_session_check = paladin_backend.session_check

        def session_check(cookie, need_seconds, **options):
            if under_the_real_session_directory(cookie):
                violation(f"the real session file was read: {cookie}")
            return real_session_check(cookie, need_seconds, **options)

        patches = [mock.patch.object(paladin_backend, "RunWatchers", DefaultStarterReached),
                   mock.patch.object(paladin_backend, "session_check", session_check)]
        if not self.REAL_WATCHER_PROCESSES:
            patches.append(mock.patch.object(REAL_RUN_WATCHERS, "start", real_start))
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)
        self.addCleanup(lambda: self.assertEqual(self.guard_violations, [], "test guard (decision 1836 m-7)"))


class RecordingWatchers:
    """A stand-in for paladin_backend.RunWatchers: records each watcher start (argv, directory), none runs;
    the START_LIMIT-th start request fails the test."""

    def __init__(self):
        self.started, self.live = [], set()

    def alive(self, directory):
        return directory in self.live

    END_LINE = "END exit 0: the row is terminal (the recording stand-in)"

    def start(self, argv, directory):
        if len(self.started) + 1 >= START_LIMIT:
            raise AssertionError(f"watcher start request {len(self.started) + 1}: the test's bound (decision 1830 (2))")
        self.started.append((argv, directory))
        self.live.add(directory)
        os.makedirs(directory, exist_ok=True)  # its events.log: a watcher that ended with its run (COMPLETION_RULE (5))
        with open(os.path.join(directory, "events.log"), "a") as f:
            f.write(f"2026-11-15T00:00:00Z START (the recording stand-in)\n2026-11-15T00:00:00Z {self.END_LINE}\n")
        return 4242

    def runs(self):
        return [int(argv[argv.index("--run-id") + 1]) for argv, _ in self.started]

    def covered(self, runs):
        """The fake's run ids with a recorded watcher: by --run-id, or by --nonce (paladin_run_watcher's NONCE_RULE:
        the nonce watcher finds its run; WatcherInProcess shows that part)."""
        started = [(argv, directory) for argv, directory in self.started if directory in self.live]  # the live ones
        nonces = {argv[argv.index("--nonce") + 1] for argv, _ in started if "--nonce" in argv}
        return ({int(argv[argv.index("--run-id") + 1]) for argv, _ in started} - {0} |
                {int(run_id) for run_id, run in runs.items() if (run.get("metadata") or {}).get("nonce") in nonces})


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


class ResponseArchiveKeys(GuardedTestCase):
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


class BlockSizingLookup(GuardedTestCase):
    def setUp(self):
        super().setUp()
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


class FakeCLICase(GuardedTestCase):
    """The fake cluster, tools and records of the drive tests (no test of its own)."""

    def setUp(self):
        super().setUp()
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

    def watchers_root(self):
        """WATCHER_RULE: one root per plan (decision 1836 m-8)."""
        return Path(os.path.realpath(self.plan_path.parent)) / "paladin-watchers"


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
        self.watchers = RecordingWatchers()  # a fresh recording stand-in (its start bound counts per fake cluster)

    # Decision 1836 m-2 / m-3 (review r1833 MINOR-2 / MINOR-3): an unknown receipt keeps the intent open; the close tool.
    def submit_without_an_answer(self, *, answer=None):
        """The submit call makes the run and then gives no usable answer (the 600-s timeout form; or `answer`)."""
        real = paladin_backend.paladin_call

        def call(cli, arguments, **options):
            if arguments[:2] == ["palace", "submit-detailed"]:
                real(cli, arguments, **options)
                if answer is not None:
                    return answer
                raise paladin_backend.PaladinUnanswered("palace submit-detailed: no answer in 600 s")
            return real(cli, arguments, **options)

        return mock.patch.object(paladin_backend, "paladin_call", side_effect=call)

    def test_b2_a_submit_without_an_answer_keeps_the_intent_open_and_never_submits_again(self):
        """Review r1833 B2: the submit made the run and then timed out; the listings answer empty (was: Outcome 'error',
        then a resume with two empty answers submitted a SECOND copy)."""
        self.lookups(["true", "true", "empty", "empty", "empty"])  # the two lookups before the submit, three after it
        del self.slept[:]
        with self.submit_without_an_answer(), self.assertRaises(paladin_backend.PaladinError) as stop:
            self.drive()
        self.assertIn("receipt is unknown", str(stop.exception))
        self.assertIn("RECOVERY", str(stop.exception))
        intent = self.records()["w"]["SubmitIntents"][-1]
        self.assertNotIn("Outcome", intent)
        self.assertEqual(self.records()["w"]["SubmitErrors"][-1]["Receipt"], "unknown")
        self.assertEqual(self.slept, [20.0, 30.0, 30.0])
        self.lookups(["empty"] * 4)  # a resume meets empty answers again (its first lookup + 3): STOP, no second submit
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), ["chainport-w"])
        state = self.drive()  # a true listing: the run is adopted
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual((state["Stages"]["w"]["SubmitIntents"][-1]["NonceMatches"], state["Stages"]["w"]["State"]), (True, "complete"))

    def test_a_submit_answer_without_a_run_id_keeps_the_intent_open(self):
        with self.submit_without_an_answer(answer={"jobs": []}):
            state = self.drive()  # the run made by the call is found by the repeated lookup and adopted
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertTrue(state["Stages"]["w"]["SubmitIntents"][-1]["Outcome"].startswith("adopted"))

    def test_b3_an_intent_whose_submit_was_never_sent_is_closed_by_the_tool_after_recorded_lookups(self):
        """Review r1833 B3: every resume STOPs with no recovery named, and only a hand edit could close the intent."""
        real = paladin_backend.paladin_call

        def ended_before_the_call(cli, arguments, **options):
            if arguments[:2] == ["palace", "submit-detailed"]:
                raise KeyboardInterrupt("the driver process ends after the intent, before the call")
            return real(cli, arguments, **options)

        with mock.patch.object(paladin_backend, "paladin_call", side_effect=ended_before_the_call), \
                self.assertRaises(KeyboardInterrupt):
            self.drive()
        for _ in range(2):
            with self.assertRaises(paladin_backend.PaladinError) as stop:
                self.drive()
            self.assertIn("RECOVERY", str(stop.exception))
            self.assertIn("--close-intent w", str(stop.exception))
        self.assertEqual(self.submitted(), [])
        profile = self.root / "profile.json"
        profile.write_text(json.dumps({"Paladin": self.block}))
        del self.slept[:]
        with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
            code = paladin_backend.main([str(self.plan_path), "--profile", str(profile), "--out", str(self.out), "--close-intent", "w"])
        self.assertEqual(code, 0)
        intent = self.records()["w"]["SubmitIntents"][-1]
        self.assertEqual((intent["Outcome"], len(intent["ClosedBy"]["Lookups"])), ("closed: no run after the recorded lookups", 3))
        self.assertEqual(self.slept, [30.0, 30.0])
        state = self.drive()  # the resume submits the stage once
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])
        self.assertEqual(state["Stages"]["w"]["State"], "complete")

    def test_the_close_tool_refuses_when_a_run_of_the_stage_is_listed(self):
        self.killed_after_the_submit()  # an open intent AND a live run
        profile = self.root / "profile.json"
        profile.write_text(json.dumps({"Paladin": self.block}))
        with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
            code = paladin_backend.main([str(self.plan_path), "--profile", str(profile), "--out", str(self.out), "--close-intent", "w"])
        self.assertEqual(code, 3)
        self.assertNotIn("Outcome", self.records()["w"]["SubmitIntents"][-1])

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
                                             log=lambda message: None, watchers=RecordingWatchers(),
                                             session_cookie=self.cookie)  # decision 1836 m-7: never the host's own
            try:
                paladin_backend.drive(job, poll_seconds=0, sleep=lambda seconds: None)
                outcomes["B"] = "drove"
            except paladin_backend.PaladinError as error:
                outcomes["B"] = str(error)

        starts = []

        def find(cli, identity):
            if not outcomes:  # the first driver's first lookup: the second driver starts now, in its own thread
                starts.append(1)
                self.assertLess(len(starts), START_LIMIT)  # decision 1830 (2)
                outcomes["B"] = "started"
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

    def test_one_failed_poll_below_the_bound_is_re_polled_and_the_run_completes(self):
        self.describes(["true", "empty", "empty", "empty"])  # three bad answers in a row = ONE failed poll
        state = self.drive()
        self.assertEqual(state["Stages"]["w"]["State"], "complete")
        self.assertEqual(self.slept.count(10.0), 1)
        self.assertEqual(self.submitted(), ["chainport-w", "chainport-r", "chainport-local-edge"])

    def test_a_stop_right_after_a_submit_leaves_the_run_with_its_watcher(self):
        recording = self.watchers.start

        def the_submit_window_watcher_ends_at_once(argv, directory):
            pid = recording(argv, directory)
            if "/nonce-" in directory:
                self.watchers.live.discard(directory)  # it ended (e.g. its search met a failed listing first)
            return pid

        self.watchers.start = the_submit_window_watcher_ends_at_once
        self.lookups(["true", "true", "fail"])  # the lookup after the submit fails: STOP
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.submitted(), ["chainport-w"])
        self.assertEqual(self.live_runs_without_a_watcher(), [])

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
        return sorted(set(live) - self.watchers.covered(runs))

    def test_every_submitted_run_gets_a_watcher_with_its_stage_cap(self):
        """Decision 1836 m-1: the watcher of each run is started BEFORE its submit call, keyed by the intent's nonce
        (NONCE_RULE); it learns the run id from the hint file; no second watcher is started for the run."""
        state = self.drive()
        self.assertEqual(len(self.watchers.started), 3)
        for (argv, directory), name in zip(self.watchers.started, ("w", "r", "local-edge")):
            intent = state["Stages"][name]["SubmitIntents"][-1]
            for flag, value in (("--run-id", "0"), ("--nonce", intent["Nonce"]), ("--cap-seconds", "900.0"), ("--host", "local"),
                                ("--bar", "0.85"), ("--grace", "120.0"), ("--submitted", intent["UTC"])):
                self.assertEqual(argv[argv.index(flag) + 1], value)
            self.assertEqual(json.loads(argv[argv.index("--identity") + 1])["stage"], name)
            self.assertEqual(directory, str(self.watchers_root() / f"nonce-{intent['Nonce']}"))
            self.assertEqual((Path(directory) / "run-id").read_text().strip(), str(state["Stages"][name]["RunID"]))
            self.assertEqual(intent["RunID"], state["Stages"][name]["RunID"])

    def test_b1_a_driver_ended_in_the_submit_window_leaves_its_run_watched(self):
        """Review r1833 B1: the driver ends between the submit and the save, the resume meets a rows_omitted lookup
        and STOPs: the live run is watched by the nonce watcher started before the submit (was: no watcher)."""
        self.killed_after_the_submit()
        self.lookups(["omitted"])
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        argv, _ = self.watchers.started[0]
        self.assertEqual(argv[argv.index("--nonce") + 1], self.records()["w"]["SubmitIntents"][-1]["Nonce"])

    def test_a_submit_window_watcher_that_does_not_start_means_no_submit_call(self):
        self.watchers.start = lambda argv, directory: None
        with self.assertRaises(paladin_backend.WatcherNotStarted):
            self.drive()
        self.assertEqual(self.submitted(), [])
        intent = self.records()["w"]["SubmitIntents"][-1]
        self.assertEqual((intent["Outcome"], intent["Watcher"]["Started"]), ("not sent: the submit-window watcher did not start", 2))

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
        run_id = self.records()["w"]["RunID"]
        self.assertEqual([Path(directory).name for _, directory in self.watchers.started],
                         [f"nonce-{self.records()['w']['SubmitIntents'][-1]['Nonce']}", f"run-{run_id}"])
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
        self.assertEqual(self.live_runs_without_a_watcher(), [])  # decision 1836 m-1: its nonce watcher (was: unwatched)
        state = json.loads(self.state.read_text())
        first = state["NextID"] - 1
        state["Runs"][str(first + 1)] = dict(state["Runs"][str(first)], run_dir=str(self.root / "runs" / "twin"),
                                             job_id="twin-job")
        state["NextID"] = first + 2
        self.state.write_text(json.dumps(state))
        with mock.patch.dict(os.environ, {"FAKE_CANCEL_FAIL": "1"}), self.assertRaises(paladin_backend.DuplicateRuns):
            self.drive()  # two runs: the newer one's early end fails
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        caps = {int(argv[argv.index("--run-id") + 1]): argv[argv.index("--cap-seconds") + 1] for argv, _ in self.watchers.started
                if argv[argv.index("--run-id") + 1] != "0"}  # the run watchers (the nonce watcher of the window aside)
        self.assertEqual(caps, {first: "900.0", first + 1: "0.0"})

    def test_b4_a_watcher_that_dies_while_the_driver_lives_is_restarted_at_the_next_poll(self):
        """Review r1833 B4: one start over 6 polls with the watcher dead (was); now the next poll restarts it, and a
        watcher that keeps dying is restarted at most 3 times (then STOP)."""
        self.script({"chainport-w": "hang"})
        polls = []

        def dies_each_time(seconds):
            polls.append(seconds)
            self.watchers.live.clear()  # the watcher process is gone after every poll

        with mock.patch.object(paladin_backend, "DRIVE_MAX_ADVANCES", 6), self.assertRaises(paladin_backend.PaladinError) as stop:
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=dies_each_time)
        self.assertIn("3 watcher starts in this driver process", str(stop.exception))
        run_id = self.records()["w"]["RunID"]
        self.assertEqual([Path(directory).name for _, directory in self.watchers.started],
                         [f"nonce-{self.records()['w']['SubmitIntents'][-1]['Nonce']}"] + [f"run-{run_id}"] * 3)
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_no_watcher_restart_after_its_end_json(self):
        self.script({"chainport-w": "hang"})
        polls = []

        def ended_then_dead(seconds):
            polls.append(seconds)
            record = self.records()["w"]
            directory = self.watchers_root() / f"run-{record['RunID']}"
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "end.json").write_text(json.dumps({"Kind": "cap", "Verified": True}))
            self.watchers.live.clear()

        with mock.patch.object(paladin_backend, "DRIVE_MAX_ADVANCES", 4), self.assertRaises(paladin_backend.PaladinError) as stop:
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended_then_dead)
        self.assertIn("4 advances without completion", str(stop.exception))
        self.assertEqual(len(self.watchers.started), 1)  # the nonce watcher only

    def test_two_out_directories_of_one_plan_share_one_watcher_per_run(self):
        """Review r1833 MINOR-8: two --out directories adopting one run gave two watchers; now one root per plan."""
        def ended(seconds):
            raise KeyboardInterrupt("the driver process ends in its wait")

        self.script({"chainport-w": "hang"})
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        run_id = self.records()["w"]["RunID"]
        self.watchers.live.add(str(self.watchers_root() / f"run-{run_id}"))  # the nonce watcher holds run-<id> (NONCE_RULE)
        other = paladin_backend.PaladinJob(str(self.plan_path), block=self.block, remote_dir=str(self.plan_path.parent),
                                           local_dir=str(self.root / "out-B"), archive_estimates_gb={self.archive: 0.001},
                                           log=lambda message: None, watchers=self.watchers, session_cookie=self.cookie)
        self.assertEqual(other.watchers_root(), self.job().watchers_root())
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(other, poll_seconds=0, sleep=ended)  # the second --out directory adopts the run
        self.assertEqual(json.loads((self.root / "out-B" / "paladin-status.json").read_text())["Stages"]["w"]["RunID"], run_id)
        self.assertEqual(len(self.watchers.started), 1)  # no second watcher for the run
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_an_adoption_starts_the_missing_watcher_before_any_poll(self):
        self.killed_after_the_submit()
        self.watchers.live.clear()  # the submit-window watcher is gone

        def ended(job, stage, record):
            raise KeyboardInterrupt("the driver process ends right after the adoption")

        with mock.patch.object(paladin_backend.PaladinJob, "poll", ended), self.assertRaises(KeyboardInterrupt):
            self.drive()
        self.assertEqual(self.live_runs_without_a_watcher(), [])
        self.assertEqual(Path(self.watchers.started[-1][1]).name, f"run-{self.records()['w']['RunID']}")

    def test_a_run_whose_watcher_does_not_start_is_ended_early(self):
        recording = self.watchers.start
        self.watchers.start = lambda argv, directory: None if "/run-" in directory else recording(argv, directory)
        self.killed_after_the_submit()
        self.watchers.live.clear()  # the nonce watcher is gone: the adoption must start a run watcher, which does not start
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
        directory = self.watchers_root() / f"run-{record['RunID']}"  # where the nonce watcher writes (NONCE_RULE)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "end.json").write_text(json.dumps({"Kind": "cap"}))
        paladin_backend.paladin_call(self.block["CLI"], ["palace", "cancel", str(record["RunID"]), "--json"])
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        self.assertEqual(self.records()["w"]["State"], "timed-out")
        self.assertEqual(self.submitted(), ["chainport-w"])

    def test_an_end_its_watcher_could_not_verify_stops_every_advance(self):
        """Decision 1836 (1)(b): END-UNVERIFIED.json of a run watcher: the driver and every resume STOP, naming the run."""
        def ended(seconds):
            raise KeyboardInterrupt("the driver process ends in its wait")

        self.script({"chainport-w": "hang"})
        with self.assertRaises(KeyboardInterrupt):
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=ended)
        record = self.records()["w"]
        directory = self.watchers_root() / f"run-{record['RunID']}"
        directory.mkdir(parents=True, exist_ok=True)  # the recording stand-in writes nothing
        (directory / paladin_backend.UNVERIFIED_MARKER).write_text(json.dumps(
            {"RunID": record["RunID"], "Job": record["JobID"], "Check": "by hand: paladin palace status ..."}))
        for _ in range(2):  # the driver, then a resume (a bound of 3 advances: without the STOP the test fails fast)
            count = len(self.calls())
            with mock.patch.object(paladin_backend, "DRIVE_MAX_ADVANCES", 3), \
                    self.assertRaises(paladin_backend.PaladinError) as stop:
                self.drive()
            self.assertIn("NOT VERIFIED", str(stop.exception))
            self.assertIn(str(record["RunID"]), str(stop.exception))
            self.assertEqual(len(self.calls()), count)  # nothing asked, nothing submitted
        self.assertEqual(json.loads((self.out / "paladin-status.json").read_text())["UnverifiedEnds"][0]["RunID"], record["RunID"])
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

    # Decision 1836 m-6 (review r1833 MINOR-6): the completion read = the five conditions exactly.
    def test_every_error_line_is_judged_not_only_the_first_twenty(self):
        self.script({"chainport-w": "failed-deep-error"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        read = self.records()["w"]["CompletionRead"]
        self.assertFalse(read["Conditions"]["3-OnlyUploadErrors"])
        self.assertEqual(read["ErrorCount"], 26)
        self.assertIn("artifact index", read["Offending"][0])

    def test_the_traceback_s_last_line_is_judged(self):
        for end, offending in (("failed-stray-exception", "ValueError: the result table was truncated"),
                               ("failed-foreign-traceback", "Traceback (most recent call last):")):
            with self.subTest(end=end):
                self.setUp_fresh()
                self.script({"chainport-w": end})
                with self.assertRaises(paladin_backend.PaladinError):
                    self.drive()
                read = self.records()["w"]["CompletionRead"]
                self.assertFalse(read["Conditions"]["3-OnlyUploadErrors"])
                self.assertEqual(read["Offending"][0], offending)
                self.assertEqual(self.records()["w"]["State"], "failed")

    def test_the_harness_s_own_traceback_passes_and_is_quoted(self):
        holds, offending, last = paladin_backend.judge_console(
            [(1, "ERROR qdash_das.palace.runner.process_harness: harness upload final tick failed"),
             (2, "Traceback (most recent call last):"), (9, "KeyError: 'x'"), (10, "Traceback (most recent call last):"),
             (14, "pyarrow.lib.ArrowTypeError: (1)")], 1)
        self.assertEqual((holds, offending, last), (True, [], "pyarrow.lib.ArrowTypeError: (1)"))  # a chained exception
        self.assertFalse(paladin_backend.judge_console([(1, "ERROR x: harness upload")], 2)[0])  # not all ERROR lines read

    def test_every_name_difference_is_reported(self):
        self.script({"chainport-w": "failed-upload-extra", "chainport-local-edge": "completed-extra"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive()
        read = self.records()["w"]["CompletionRead"]
        self.assertEqual((read["ExtraOutputs"], read["Conditions"]["4-ExpectedOutputs"]), (["surprise.csv"], False))  # NOTE-1
        self.setUp_fresh()
        self.script({"chainport-local-edge": "completed-extra"})
        state = self.drive()  # a completed row: the extra name is reported, not a stop
        self.assertEqual(state["Stages"]["local-edge"]["CompletionRead"]["NameDifferences"],
                         {"Missing": [], "Extra": ["surprise.csv"]})
        self.assertEqual(state["Stages"]["local-edge"]["State"], "complete")

    def test_the_watcher_s_end_line_is_read(self):
        for line, holds in ((RecordingWatchers.END_LINE, True), ("END exit 11: verified (the row terminal)", False),
                            ("CAP 1000 s from the PBS start >= cap 900 s -> the end request", False), (None, False)):
            with self.subTest(line=line):
                self.setUp_fresh()
                del self.slept[:]
                self.script({"chainport-w": "failed-upload"})
                RecordingWatchers.END_LINE, saved = (line or "START again"), RecordingWatchers.END_LINE
                try:
                    if holds:
                        self.drive()
                    else:
                        with self.assertRaises(paladin_backend.PaladinError):
                            self.drive()
                finally:
                    RecordingWatchers.END_LINE = saved
                read = self.records()["w"]["CompletionRead"]
                self.assertEqual(read["Conditions"]["5-NoStop"], holds)
                if line is None:
                    self.assertEqual(self.slept.count(paladin_backend.WATCHER_INTERVAL_SECONDS), 3)  # waited for the END line

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

    def test_a_stage_type_without_an_output_list_is_not_submitted(self):
        plan = json.loads(self.plan_path.read_text())
        plan["Stages"][0]["Environment"] = {"PALACE_RESPONSE_ESTIMATE_ONLY": "1", "PALACE_RESPONSE_ARCHIVE_DIR": self.archive}
        self.plan_path.write_text(json.dumps(plan))
        with self.assertRaises(paladin_backend.PaladinError) as stop:
            self.drive()
        self.assertIn("no expected-output list", str(stop.exception))
        self.assertEqual(self.submitted(), [])

    def test_an_archive_marked_unusable_does_not_verify(self):
        self.drive()
        job = self.job()
        config = job.stage_run_config(job.plan["Stages"][0])
        self.assertTrue(job.archive_verify(self.archive, config, {"Name": "w"})["OK"])
        (Path(self.archive) / paladin_backend.UNUSABLE_MARKER).write_text("{}")
        verification = job.archive_verify(self.archive, config, {"Name": "w"})
        self.assertEqual((verification["OK"], verification["Unusable"]), (False, True))

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


class CountingRunWatchers(REAL_RUN_WATCHERS):
    """The real watcher starts, counted: the START_LIMIT-th start request fails the test (decision 1830 (2))."""

    def __init__(self):
        super().__init__()
        self.requests = 0

    def start(self, argv, directory, **options):
        self.requests += 1
        if self.requests >= START_LIMIT:
            raise AssertionError(f"real watcher start request {self.requests}: the test's bound (decision 1830 (2))")
        return super().start(argv, directory, **options)


class RunWatcherProcesses(FakeCLICase):
    """Decision 1820 (4): REAL paladin_run_watcher.py processes (nohup, their own session) against the fake
    CLI and a fake PBS (qstat / qdel) and exec-host ssh on PATH; the watcher's interval 1 s, grace 5 s,
    the worker's cap 3 s.  THE ONE OPT-IN CLASS of real watcher processes (ruling 1837): the counting starter
    (the 5th start request fails the test), only the fake CLI / PBS / session file in a temp directory, and every
    started process ended by its recorded pid in tearDown (a leftover fails the test)."""

    REAL_WATCHER_PROCESSES = True

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
        self.watchers = CountingRunWatchers()
        for path in (self.root, self.cookie, self.out, fake_bin, self.block["CLI"], self.state):  # ruling 1837 (1)
            self.assertFalse(under_the_real_session_directory(path), path)

    def tearDown(self):
        for patch in self.timing:
            patch.stop()
        self.path.stop()
        leftovers = []
        for process in self.watchers.processes:  # ruling 1837 (2): by the recorded pid, never by name or group
            if not self.wait_for(lambda: process.poll() is not None, 30):
                leftovers.append(process.pid)
                process.terminate()
                try:
                    process.wait(10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(10)
        super().tearDown()
        self.assertEqual(leftovers, [], "watcher processes left running at the end of the test (ended by pid)")

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
        directory = self.watchers_root() / f"run-{record['RunID']}"  # the nonce watcher moves there (NONCE_RULE)
        self.assertTrue(self.wait_for(lambda: self.watchers.alive(str(directory)), 30))  # it outlives the driver
        self.assertTrue(self.wait_for(lambda: (directory / "end.json").is_file(), 30))
        self.assertTrue(self.wait_for(lambda: not self.watchers.alive(str(directory)), 30))
        return record, json.loads((directory / "end.json").read_text()), (directory / "events.log").read_text()

    def test_the_watcher_ends_the_run_at_its_cap_after_the_driver_ended(self):
        record, end, events = self.ended_in_its_wait()
        self.assertEqual((end["Kind"], end["Rounds"][0]["Cancel"]["OK"], end["Verified"]), ("cap", True, True))
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
        self.assertFalse(end["Rounds"][0]["Cancel"]["OK"])
        self.assertEqual(end["Rounds"][0]["Qdel"]["Command"], f"{self.fake_bin}/qdel {record['JobID']}")
        self.assertEqual((end["Rounds"][0]["Qdel"]["ReturnCode"], end["Verified"]), (0, True))
        self.assertIn("qdel", events)
        self.assertEqual(json.loads(self.state.read_text())["Runs"][str(record["RunID"])]["status"], "canceled")

    def test_a_watcher_ends_with_its_run(self):
        plan = json.loads(self.plan_path.read_text())
        plan["Stages"][0]["CapSeconds"] = 3600
        self.plan_path.write_text(json.dumps(plan))
        state = self.drive()
        for name in ("w", "r", "local-edge"):
            nonce = Path(state["Stages"][name]["SubmitIntents"][-1]["Watcher"]["Dir"])
            self.assertTrue(self.wait_for(lambda: not self.watchers.alive(str(nonce)), 30))
            directory = self.watchers_root() / f"run-{state['Stages'][name]['RunID']}"
            self.assertTrue(self.wait_for(lambda: not self.watchers.alive(str(directory)), 30))
            self.assertFalse((directory / "end.json").exists())
            self.assertRegex((directory / "events.log").read_text(), r"END exit 0: the (row is terminal|PBS job has ended)")
        self.assertEqual(state["Stages"]["local-edge"]["State"], "complete")


class Bounds(FakeCLICase):
    """Decision 1830 (3): the bounds in code (no real process is started here)."""

    def test_at_most_three_watcher_starts_per_run_id_per_driver_process(self):
        self.watchers.start = lambda argv, directory: self.watchers.started.append(argv) or None  # never starts
        job = self.job()
        with self.assertRaises(paladin_backend.WatcherNotStarted):
            job.ensure_watcher("w", 999, job_id=None, submitted=paladin_backend.utc_now().isoformat(), cap_seconds=60)
        with self.assertRaises(paladin_backend.PaladinError) as stop:
            job.ensure_watcher("w", 999, job_id=None, submitted=paladin_backend.utc_now().isoformat(), cap_seconds=60)
        self.assertIn("3 watcher starts in this driver process", str(stop.exception))
        self.assertEqual(len(self.watchers.started), 3)
        self.assertEqual(job.watcher_starts, {999: 3})

    def test_one_watcher_directory_per_run_id(self):
        job = self.job()
        job.ensure_watcher("w", 7, job_id="j", submitted=paladin_backend.utc_now().isoformat(), cap_seconds=60)
        job.ensure_watcher("other-name", 7, job_id="j", submitted=paladin_backend.utc_now().isoformat(), cap_seconds=60)
        self.assertEqual([directory for _, directory in self.watchers.started], [str(self.watchers_root() / "run-7")])

    def test_the_drive_loop_is_bounded_and_sleeps_at_least_its_minimum(self):
        self.script({"chainport-w": "hang"})
        slept = []
        with mock.patch.object(paladin_backend, "DRIVE_MAX_ADVANCES", 4), self.assertRaises(paladin_backend.PaladinError) as stop:
            paladin_backend.drive(self.job(), poll_seconds=0, sleep=slept.append)
        self.assertIn("4 advances without completion", str(stop.exception))
        self.assertEqual(slept, [paladin_backend.DRIVE_MIN_POLL_SECONDS] * 4)


class FakeTime:
    """A clock for the in-process watcher: sleep() advances time() (no real wait); the sleeps recorded."""

    def __init__(self):
        self.now, self.slept = time.time(), []

    def time(self):
        return self.now

    def sleep(self, seconds):
        self.slept.append(seconds)
        self.now += seconds


class WatcherInProcess(GuardedTestCase):
    """paladin_run_watcher.Watcher run IN THIS PROCESS against scripted reads and a fake clock (no process started;
    decision 1830 (2)): the cap, the memory bar, the qdel fallback, the read bound, one watcher per directory; and
    decision 1836 (review r1833 MAJOR-1): W1 / W2 / W3, the expired session, an unverified end, the blind clock."""

    def setUp(self):
        super().setUp()
        import paladin_run_watcher
        self.module = paladin_run_watcher
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name) / "watchers" / "run-7"
        self.dir.mkdir(parents=True)
        self.clock = FakeTime()
        self.calls, self.status, self.cancel = [], "running", "ok"  # cancel: ok | refused | accepted (the row stays)
        self.status_fails = False  # the row's read missing
        self.qdel = []  # the outcomes of the next qdel calls: ok | rc1 | ssh; none left = ok
        self.pbs = "R"  # R | Q | NOTFOUND | F | ssh (the cluster read fails)
        self.start_ago, self.start_at, self.memory = 100, None, (1000, 900)  # the PBS start: now - start_ago, or start_at
        self.commands = []

        def paladin_call(cli, arguments, **options):
            self.calls.append(arguments[:2])
            if len(self.calls) > 200:  # the fake's own bound: a watcher without its read bound fails here, it never spins
                raise AssertionError("200 reads: the watcher did not end")
            if arguments[1] == "cancel":
                if self.cancel == "refused":
                    raise paladin_backend.PaladinError("palace cancel: 'failed transport: the session has expired; run mwinit'")
                if self.cancel == "ok":
                    self.status = "canceled"
                return {}
            if self.status_fails:
                raise paladin_backend.PaladinError("palace status: no answer")
            return {"run_id": 7, "job_id": "j7", "status": self.status,
                    "terminal": self.status in ("canceled", "completed", "failed")}

        def remote_shell(host, command, **options):
            self.calls.append(command.split()[0].split("/")[-1])
            self.commands.append(command)
            if "qdel" in command:
                outcome = self.qdel.pop(0) if self.qdel else "ok"
                if outcome == "ssh":
                    raise paladin_backend.PaladinError("ssh soca-green-job: no answer in 120 s")
                if outcome == "ok":
                    self.pbs = "F"
                return subprocess.CompletedProcess(["qdel"], 0 if outcome == "ok" else 1, stdout="",
                                                   stderr="" if outcome == "ok" else "qdel: Server connection failed")
            if self.pbs == "ssh":
                raise paladin_backend.PaladinError("ssh soca-green-job: no answer in 120 s")
            kept = command.split()[5:]  # bash -s -- JOB PBS_BIN [KEPT_HOSTS...]
            stime = self.start_at or int(self.clock.now) - self.start_ago
            lines = [f"STATE {self.pbs}", f"STIME {stime if self.pbs == 'R' else '-'}",
                     f"EXIT {0 if self.pbs == 'F' else '-'}"]
            hosts = ["node"] if self.pbs == "R" else kept if self.pbs == "NOTFOUND" else []
            if self.pbs == "R":
                lines.append("HOSTS node")
            lines += [f"MEM {host} {self.memory[0]} {self.memory[1]}" for host in hosts]
            return subprocess.CompletedProcess(["bash"], 0, stdout="\n".join(lines) + "\n", stderr="")

        self.listings = []  # the next my-jobs answers of the nonce phase: (rows, warnings); none left = no rows

        def paladin_envelope(cli, arguments, **options):
            self.calls.append(arguments[:2])
            rows, warnings = self.listings.pop(0) if self.listings else ([], [])
            return {"ok": True, "data": {"jobs": rows}, "warnings": warnings}

        self.patches = [mock.patch.object(paladin_backend, "paladin_call", paladin_call),
                        mock.patch.object(paladin_backend, "paladin_envelope", paladin_envelope),
                        mock.patch.object(paladin_backend, "remote_shell", remote_shell),
                        mock.patch.object(self.module, "time", self.clock)]
        for patch in self.patches:
            patch.start()

    def tearDown(self):
        for patch in self.patches:
            patch.stop()
        self.tmp.cleanup()

    def run_watcher(self, cap=50.0, bar=0.85, submitted_ago=0.0, interval=1.0):
        submitted = paladin_backend.utc_now() - datetime.timedelta(seconds=submitted_ago)
        args = ["--run-id", "7", "--job", "j7", "--cli", "paladin", "--host", "local", "--pbs-bin", "/pbs",
                "--cap-seconds", str(cap), "--submitted", submitted.isoformat(), "--dir", str(self.dir),
                "--bar", str(bar), "--interval", str(interval), "--grace", "1"]
        return self.module.main(args)

    def end(self):
        return json.loads((self.dir / "end.json").read_text())

    def test_the_cap_from_the_pbs_start(self):
        self.assertEqual(self.run_watcher(cap=50), 11)
        self.assertEqual((self.end()["Kind"], self.end()["Rounds"][0]["Cancel"]["OK"], self.end()["Verified"]), ("cap", True, True))
        self.assertIn("from the PBS start", self.end()["Reason"])
        self.assertIn("row terminal", self.end()["VerifiedBy"])

    def test_the_memory_bar(self):
        self.memory = (1000, 100)
        self.assertEqual(self.run_watcher(cap=10**6), 10)
        self.assertEqual(self.end()["Kind"], "memory")

    def test_a_failed_cancel_falls_back_to_qdel(self):
        self.cancel = "refused"
        self.assertEqual(self.run_watcher(cap=50), 11)
        round_one = self.end()["Rounds"][0]
        self.assertEqual((round_one["Qdel"]["Command"], self.end()["Verified"]), ("/pbs/qdel j7", True))
        self.assertIn("PBS job_state F", self.end()["VerifiedBy"])

    def test_w1_an_empty_qstat_answer_throughout_keeps_the_cap_in_force_by_the_row_s_status(self):
        """Review r1833 W1: qstat answers empty (NOTFOUND) throughout while the row runs 1000 s past a 50-s cap."""
        self.pbs = "NOTFOUND"
        self.assertEqual(self.run_watcher(cap=50, submitted_ago=1000), 11)
        self.assertEqual(self.calls.count(["palace", "cancel"]), 1)
        self.assertIn("from the submit time", self.end()["Reason"])
        self.assertIn("running by the row's status (a missing PBS read)", self.end()["Reason"])
        self.assertTrue(self.end()["Verified"])

    def test_w2_an_empty_qstat_answer_after_a_failed_qdel_is_not_an_end(self):
        """Review r1833 W2: the end request refused, qdel rc 1 ('Server connection failed'), then an empty qstat answer:
        NOT ended (was: JobEnded true, exit 11); the second round's qdel is verified by the finished-job answer."""
        self.cancel, self.qdel = "refused", ["rc1", "ok"]
        real = paladin_backend.remote_shell

        def empty_after_the_first_qdel(host, command, **options):
            result = real(host, command, **options)
            if "qdel" in command and result.returncode == 1:
                self.pbs = "NOTFOUND"  # the next reads answer empty, the job still runs
            return result

        with mock.patch.object(paladin_backend, "remote_shell", empty_after_the_first_qdel):
            self.assertEqual(self.run_watcher(cap=50), 11)
        rounds = self.end()["Rounds"]
        self.assertEqual([entry["VerifiedBy"] for entry in rounds][0], None)
        self.assertEqual((len(rounds), rounds[1]["Qdel"]["ReturnCode"]), (2, 0))
        self.assertIn("PBS job_state F", self.end()["VerifiedBy"])
        self.assertEqual(self.calls.count(["palace", "cancel"]), 2)  # the end request again, then the qdel again
        self.assertIn(1.0, self.clock.slept)  # the back-off of round 1 (grace 1 s x 2^0)

    def test_w3_a_qdel_that_meets_an_ssh_error_is_repeated(self):
        """Review r1833 W3: the end request refused, then an ssh error on the qdel (was: exit 11, no retry)."""
        self.cancel, self.qdel = "refused", ["ssh", "ssh", "ok"]
        self.assertEqual(self.run_watcher(cap=50), 11)
        rounds = self.end()["Rounds"]
        self.assertEqual([("Error" in entry["Qdel"], entry["VerifiedBy"] is None) for entry in rounds],
                         [(True, True), (True, True), (False, False)])
        self.assertIn(2.0, self.clock.slept)  # the back-off of round 2 (grace 1 s x 2^1)

    def test_an_expired_session_ends_the_run_by_the_qdel_by_job_id(self):
        """Decision 1836 (1)(c): a run that queued beyond the session's life: the end request is refused at its cap
        (an expired session), the qdel by the PBS job id (no session) ends it, verified."""
        self.cancel = "refused"
        self.assertEqual(self.run_watcher(cap=50), 11)
        rounds = self.end()["Rounds"]
        self.assertEqual(len(rounds), 1)
        self.assertIn("session has expired", rounds[0]["Cancel"]["Error"])
        self.assertEqual(rounds[0]["Qdel"]["ReturnCode"], 0)
        self.assertTrue(self.end()["Verified"])

    def test_an_end_not_verified_at_the_bound_leaves_a_marker_and_exit_12(self):
        self.cancel, self.qdel = "refused", ["rc1"] * self.module.END_MAX_ROUNDS
        self.assertEqual(self.run_watcher(cap=50), 12)
        self.assertEqual((self.end()["Verified"], len(self.end()["Rounds"])), (False, self.module.END_MAX_ROUNDS))
        marker = json.loads((self.dir / self.module.UNVERIFIED_MARKER).read_text())
        self.assertEqual((marker["RunID"], marker["Job"]), (7, "j7"))
        self.assertEqual([seconds for seconds in self.clock.slept if seconds > 1.0][:3], [2.0, 4.0, 8.0])
        self.assertIn("END exit 12", (self.dir / "events.log").read_text())

    def test_the_end_request_accepted_but_not_verified_is_followed_by_the_qdel(self):
        self.cancel = "accepted"  # the request is accepted, the row keeps running
        self.assertEqual(self.run_watcher(cap=50), 11)
        self.assertEqual(self.end()["Rounds"][0]["Cancel"]["OK"], True)
        self.assertIn("PBS job_state F", self.end()["VerifiedBy"])

    def test_no_read_for_the_blind_time_judges_the_cap_by_the_clock(self):
        self.status_fails, self.pbs = True, "ssh"
        real = paladin_backend.paladin_call

        def recovers_at_the_end_request(cli, arguments, **options):
            if arguments[1] == "cancel":
                self.status_fails = False
            return real(cli, arguments, **options)

        with mock.patch.object(paladin_backend, "paladin_call", recovers_at_the_end_request):
            self.assertEqual(self.run_watcher(cap=50, submitted_ago=1000, interval=30), 11)
        self.assertIn("by the clock", self.end()["Reason"])
        self.assertGreaterEqual(sum(self.clock.slept), self.module.BLIND_SECONDS)

    def test_a_missing_pbs_read_keeps_the_memory_bar_on_the_kept_hosts(self):
        real = paladin_backend.remote_shell
        reads = []

        def empty_after_one_good_read(host, command, **options):
            result = real(host, command, **options)
            if "qdel" not in command:
                reads.append(command)
                self.pbs, self.memory = "NOTFOUND", (1000, 100)  # after the first read: empty answers, the memory high
            return result

        with mock.patch.object(paladin_backend, "remote_shell", empty_after_one_good_read):
            self.assertEqual(self.run_watcher(cap=10**6), 10)
        self.assertTrue(reads[1].endswith(" node"))  # the kept exec host is read under the missing PBS read
        self.assertEqual(self.end()["Kind"], "memory")

    def test_a_queued_row_under_an_empty_qstat_answer_is_not_capped(self):
        self.status, self.pbs = "queued", "NOTFOUND"
        with mock.patch.object(self.module, "WATCHER_LIFE_BEYOND_CAP_SECONDS", 300.0), \
                mock.patch.object(self.module, "QUEUE_LIFE_SECONDS", 0.0):
            self.assertEqual(self.run_watcher(cap=5, interval=30), 5)
        self.assertEqual(self.calls.count(["palace", "cancel"]), 0)

    def test_the_read_bound(self):
        self.status, self.pbs = "queued", "Q"
        with mock.patch.object(self.module, "WATCHER_LIFE_BEYOND_CAP_SECONDS", 0.0), \
                mock.patch.object(self.module, "QUEUE_LIFE_SECONDS", 0.0):
            self.assertEqual(self.run_watcher(cap=5), 5)
        # the life bound: the submit time + 0 + 5 s + 0 (the submit time is read from the real clock, the fake clock
        # started a moment before it: 5 or 6 reads)
        self.assertIn(self.calls.count(["palace", "status"]), (5, 6))
        self.assertFalse((self.dir / "end.json").exists())

    def test_the_read_count_bound_holds_without_a_clock(self):
        """Decision 1830 (3): the count bound holds by itself (a clock that never advances: the life deadline is never met)."""
        self.status, self.pbs = "queued", "Q"
        frozen = types.SimpleNamespace(time=lambda: self.clock.now, sleep=lambda seconds: None)
        with mock.patch.object(self.module, "time", frozen), \
                mock.patch.object(self.module, "WATCHER_LIFE_BEYOND_CAP_SECONDS", 0.0), \
                mock.patch.object(self.module, "QUEUE_LIFE_SECONDS", 0.0):
            self.assertEqual(self.run_watcher(cap=5), 5)
        self.assertEqual(self.calls.count(["palace", "status"]), 15)  # int((0 + 5 + 0) / 1) + 10 reads

    def test_a_run_queued_beyond_cap_plus_one_day_keeps_its_watcher_and_is_capped_from_its_pbs_start(self):
        """Review r1833 MINOR-4 (b): the life bound counted from the watcher's start ended it while queued (exit 5)."""
        self.status, self.pbs, self.start_ago = "queued", "Q", 0
        real = paladin_backend.remote_shell
        reads = []

        def starts_after_a_long_queue(host, command, **options):
            reads.append(1)
            if len(reads) == 30:  # 30 reads x 30 s = 900 s queued > cap 50 + the life beyond it 600 s (patched)
                self.status, self.pbs, self.start_at = "running", "R", int(self.clock.now)
            return real(host, command, **options)

        with mock.patch.object(paladin_backend, "remote_shell", starts_after_a_long_queue), \
                mock.patch.object(self.module, "WATCHER_LIFE_BEYOND_CAP_SECONDS", 600.0):
            self.assertEqual(self.run_watcher(cap=50, interval=30), 11)
        self.assertIn("from the PBS start", self.end()["Reason"])

    # Decision 1836 m-1 (review r1833 MINOR-1): the nonce phase of the submit-window watcher (NONCE_RULE).
    def run_nonce(self, search=600.0):
        self.nonce_dir = self.dir.parent / "nonce-abc"
        args = ["--run-id", "0", "--nonce", "abc", "--identity", json.dumps({"stage": "w"}), "--search-seconds", str(search),
                "--cli", "paladin", "--host", "local", "--pbs-bin", "/pbs", "--cap-seconds", "50",
                "--submitted", paladin_backend.utc_now().isoformat(), "--dir", str(self.nonce_dir), "--interval", "30",
                "--grace", "1"]
        return self.module.main(args)

    ROW = {"run_id": 7, "job_id": "j7", "status": "running", "metadata": {"stage": "w", "nonce": "abc"}}

    def test_the_nonce_watcher_finds_its_run_by_the_driver_s_hint_and_watches_it(self):
        (self.dir.parent / "nonce-abc").mkdir()
        (self.dir.parent / "nonce-abc" / "run-id").write_text("7\n")
        self.assertEqual(self.run_nonce(), 11)
        self.assertEqual(self.end()["RunID"], 7)  # end.json in run-7, the run's own directory
        self.assertIn("found by the driver's hint", (self.nonce_dir / "events.log").read_text())
        self.assertNotIn(["palace", "my-jobs"], self.calls)

    def test_the_nonce_watcher_finds_its_run_by_the_lookup_after_a_missing_read(self):
        self.listings = [([], [{"code": "rows_omitted"}]), ([dict(self.ROW, metadata={"nonce": "other"})], []), ([self.ROW], [])]
        self.assertEqual(self.run_nonce(), 11)
        self.assertEqual(self.calls.count(["palace", "my-jobs"]), 3)
        events = (self.nonce_dir / "events.log").read_text()
        self.assertIn("READERR palace my-jobs", events)
        self.assertIn("found by the lookup by the nonce", events)

    def test_the_nonce_watcher_ends_on_a_no_run_row(self):
        self.listings = [([dict(self.ROW, status="submission_failed", job_id=None)], [])]
        self.assertEqual(self.run_nonce(), 0)
        self.assertIn("made no run", (self.nonce_dir / "events.log").read_text())
        self.assertNotIn(["palace", "status"], self.calls)

    def test_the_nonce_watcher_ends_when_no_run_appears_in_its_search_time(self):
        self.assertEqual(self.run_nonce(search=300), 0)
        self.assertIn("no run with this nonce in 300 s", (self.nonce_dir / "events.log").read_text())
        self.assertEqual(self.calls.count(["palace", "my-jobs"]), 11)  # every 30 s over 300 s
        self.assertLessEqual(sum(self.clock.slept), 300)

    def test_the_nonce_watcher_leaves_a_run_another_watcher_holds(self):
        held = paladin_backend.DriverLock(self.dir / "watcher.lock", log=lambda line: None).acquire()
        try:
            self.listings = [([self.ROW], [])]
            self.assertEqual(self.run_nonce(), 0)
        finally:
            held.release()
        self.assertIn("watched by another watcher", (self.nonce_dir / "events.log").read_text())
        self.assertNotIn(["palace", "status"], self.calls)

    def test_one_watcher_per_directory(self):
        held = paladin_backend.DriverLock(self.dir / "watcher.lock", log=lambda line: None).acquire()
        try:
            self.assertEqual(self.run_watcher(), 3)
        finally:
            held.release()
        self.assertEqual(self.calls, [])


class GuardSelfTest(unittest.TestCase):
    """Decision 1836 m-7: the guard of GuardedTestCase fails each forbidden form (nested test runs; no process is
    started and the real session file is not opened: the guard raises before the read)."""

    def run_nested(self, body, *, real_processes=False):
        class Nested(GuardedTestCase):
            REAL_WATCHER_PROCESSES = real_processes

            def test_body(self):
                body(self)

        result = unittest.TestResult()
        Nested("test_body").run(result)
        return result

    def test_each_forbidden_form_fails_the_test(self):
        def default_starter_through_a_job(test):
            with tempfile.TemporaryDirectory() as root:
                plan_path, _ = three_stage_plan(root)
                job = paladin_backend.PaladinJob(str(plan_path), block=fake_block(root), remote_dir=root,
                                                 local_dir=os.path.join(root, "out"))
                job.dry_run()  # a dry run does not reach the default starter ...
                job.watchers  # noqa: B018 (... its first use does)

        def real_session_file(test):
            paladin_backend.session_check(paladin_backend.MIDWAY_COOKIE, 1800)

        def real_start_swallowed(test):
            try:
                REAL_RUN_WATCHERS().start(["true"], "/nonexistent/watchers/run-1")
            except AssertionError:
                pass  # swallowed by the code under test: the cleanup still fails the test

        for body in (default_starter_through_a_job, real_session_file, real_start_swallowed):
            with self.subTest(body=body.__name__):
                result = self.run_nested(body)
                self.assertFalse(result.wasSuccessful())
                self.assertIn("decision 1836 m-7", str(result.failures + result.errors))
        allowed = self.run_nested(lambda test: None, real_processes=True)
        self.assertTrue(allowed.wasSuccessful())
        self.assertIs(paladin_backend.RunWatchers, REAL_RUN_WATCHERS)  # every patch undone


@unittest.skipUnless((RECORDS / "plan.json").exists(), "records of 274fe1a150c0 absent")
class PlanOfRecord(GuardedTestCase):
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
