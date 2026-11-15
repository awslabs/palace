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


class DriveWithFakeCLI(unittest.TestCase):
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
                                          log=lambda message: None, now=now)

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

    def test_a_multi_node_plan_is_refused(self):
        plan_path, _ = three_stage_plan(self.root / "multi", nodes=2)
        with self.assertRaises(paladin_backend.PaladinError):
            paladin_backend.PaladinJob(str(plan_path), block=self.block, remote_dir=str(plan_path.parent),
                                       local_dir=str(self.out))

    def test_main_dry_run_and_submit(self):
        profile = self.root / "profile.json"
        profile.write_text(json.dumps({"Paladin": self.block}))
        common = [str(self.plan_path), "--profile", str(profile), "--out", str(self.out)]
        with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
            self.assertEqual(paladin_backend.main(common + ["--dry-run"]), 0)
            self.assertEqual(self.calls(), [])
            self.script({"chainport-local-edge": "failed"})
            self.assertEqual(paladin_backend.main(common + ["--submit", "--remote-dir", str(self.plan_path.parent),
                                                            "--archive-gb", f"{self.archive}=0.001",
                                                            "--poll-seconds", "0"]), 3)
        states = {name: record["State"] for name, record in self.records().items()}
        self.assertEqual(states, {"w": "complete", "r": "complete", "local-edge": "failed"})


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
