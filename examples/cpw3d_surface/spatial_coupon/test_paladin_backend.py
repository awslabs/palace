# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (USER decision 1674 (2); Y-SPEC s2):
qualify/paladin_backend.py turns a run_stages plan into paladin runs (one per stage) whose
configurations carry the stage variables as Solver.Electrostatic.ResponseArchive keys, and
drives them (submit / wait / fetch).
Checked here: the variable -> key table and its refusals, the run configurations, the submit
command, the Requires order and the partition check; the drive against the FAKE paladin CLI
fake_paladin_cli.py (the order, the run_id recorded at the submit and never resubmitted, a
failed run, the cap -> cancel, a refused submit, the fetch and the log parse); and, when the
records of key 274fe1a150c0 are present on this host, its fab plan of record."""
import datetime
import json
import os
from pathlib import Path
import sys
import tempfile
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


def three_stage_plan(root):
    """A worker, its reducer and a stage without archive, configurations written to root."""
    stages = []
    for name, environment, requires in (("w", WORKER, []), ("r", REDUCER, ["w"]),
                                        ("local-edge", {}, [])):
        path = os.path.join(root, name + "-in.json")
        with open(path, "w") as f:
            json.dump({"Problem": {"Type": "Electrostatic"}, "Model": {"Mesh": "m.msh"},
                       "Solver": {"Electrostatic": {}}}, f)
        stages.append({"Name": name, "Config": path, "Environment": environment,
                       "CapSeconds": 900, "Requires": requires})
    return {"Case": "case", "Ranks": 192, "Instance": {"Type": "m8g.48xlarge"},
            "Stages": stages}


class ResponseArchiveKeys(unittest.TestCase):
    def test_the_stage_forms_of_build_plan(self):
        self.assertEqual(paladin_backend.response_archive_keys(WORKER),
                         {"Mode": "Archive", "Directory": ARCHIVE, "SourceTiming": True})
        self.assertEqual(paladin_backend.response_archive_keys(REDUCER),
                         {"Mode": "Reduce", "Directory": ARCHIVE, "BlockSize": 48})
        self.assertEqual(paladin_backend.response_archive_keys({}), {})
        self.assertEqual(paladin_backend.response_archive_keys(
            {"PALACE_RESPONSE_ESTIMATE_ONLY": "1",
             "PALACE_RESPONSE_ESTIMATE_REQUEST": "/data/r.json",
             "PALACE_RESPONSE_ARCHIVE_DIR": ARCHIVE,
             "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS": "1"}),
            {"Mode": "Estimate", "EstimateRequest": "/data/r.json", "Directory": ARCHIVE,
             "RecycleInitialGuess": True})

    def test_refusals(self):
        for environment in (
                {"PALACE_RESPONSE_ARCHIVE_DIR": "archive"},
                {"PALACE_RESPONSE_ARCHIVE_ONLY": "1", "PALACE_RESPONSE_REDUCE_ONLY": "1"},
                {"PALACE_RESPONSE_ARCHIVE_ONLY": "0"},
                {"PALACE_RESPONSE_SOURCE_TIMING": "yes"},
                {"PALACE_RESPONSE_BLOCK_SIZE": "0"},
                {"PALACE_RESPONSE_DERIVED_CACHE": "0"},
                {"OMP_NUM_THREADS": "1"}):
            with self.subTest(environment=environment):
                with self.assertRaises(ValueError):
                    paladin_backend.response_archive_keys(environment)

    def test_run_config(self):
        config = {"Problem": {"Type": "Electrostatic"},
                  "Solver": {"Electrostatic": {"ResponseMatrix": True}}}
        out = paladin_backend.run_config(config, WORKER)
        self.assertEqual(out["Solver"]["Electrostatic"],
                         {"ResponseMatrix": True,
                          "ResponseArchive": {"Mode": "Archive", "Directory": ARCHIVE,
                                              "SourceTiming": True}})
        self.assertNotIn("ResponseArchive", config["Solver"]["Electrostatic"])
        self.assertEqual(paladin_backend.run_config(config, {}), config)
        with self.assertRaises(ValueError):
            paladin_backend.run_config(out, REDUCER)


class PlanRuns(unittest.TestCase):
    def plan(self, root):
        return three_stage_plan(root)

    def test_runs_and_submit_commands(self):
        with tempfile.TemporaryDirectory() as root:
            runs = paladin_backend.plan_runs(
                self.plan(root), plan_sha256="ab" * 32, selector="dev-x", lane="chainport",
                nodes=1, out=root)
            self.assertEqual([run["Stage"] for run in runs], ["w", "r", "local-edge"])
            self.assertEqual(runs[1]["Requires"], ["w"])
            self.assertEqual(runs[1]["ResponseArchive"]["Mode"], "Reduce")
            self.assertEqual(runs[2]["ResponseArchive"], {})
            with open(runs[0]["RunConfig"]) as f:
                self.assertEqual(json.load(f)["Solver"]["Electrostatic"]["ResponseArchive"],
                                 {"Mode": "Archive", "Directory": ARCHIVE,
                                  "SourceTiming": True})
            submit = runs[0]["Submit"]
            self.assertEqual(submit[:4], ["paladin", "palace", "submit-detailed",
                                          runs[0]["RunConfig"]])
            for flag, value in (("--instance-type", "m8g.48xlarge"), ("--nodes", "1"),
                                ("--ranks-per-node", "192"), ("--palace", "dev-x")):
                self.assertEqual(submit[submit.index(flag) + 1], value)
            self.assertIn("plan=abababababababab", submit)
            self.assertEqual(submit[-1], "--json")
            runs2 = paladin_backend.plan_runs(
                self.plan(root), plan_sha256="ab" * 32, selector="dev-x", lane="chainport",
                nodes=4, out=root)
            self.assertEqual(runs2[0]["RanksPerNode"], 48)

    def test_refusals(self):
        with tempfile.TemporaryDirectory() as root:
            plan = self.plan(root)
            with self.assertRaises(ValueError):
                paladin_backend.plan_runs(plan, plan_sha256="0" * 64, selector="s",
                                          lane="l", nodes=5, out=root)
            plan["Stages"][0]["Requires"] = ["r"]
            with self.assertRaises(ValueError):
                paladin_backend.plan_runs(plan, plan_sha256="0" * 64, selector="s",
                                          lane="l", nodes=1, out=root)


class DriveWithFakeCLI(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.cli = os.path.join(self.root, "paladin")
        with open(self.cli, "w") as f:
            f.write(f'#!/bin/sh\nexec "{sys.executable}" "{FAKE}" "$@"\n')
        os.chmod(self.cli, 0o755)
        self.state = os.path.join(self.root, "fake-state.json")
        self.log = os.path.join(self.root, "fake-argv.jsonl")
        self.out = os.path.join(self.root, "out")
        os.makedirs(os.path.join(self.root, "remote"))
        os.makedirs(self.out)
        self.status_path = os.path.join(self.out, "status.json")
        self.environment = mock.patch.dict(os.environ, {"FAKE_PALADIN_STATE": self.state,
                                                        "FAKE_PALADIN_LOG": self.log})
        self.environment.start()

    def tearDown(self):
        self.environment.stop()
        self.tmp.cleanup()

    def script(self, plan=None):
        with open(self.state, "w") as f:
            json.dump({"Root": os.path.join(self.root, "remote"), "Plan": plan or {}}, f)

    def runs(self):
        return paladin_backend.plan_runs(three_stage_plan(self.root), plan_sha256="cd" * 32,
                                         selector="dev-x", lane="chainport", nodes=1,
                                         out=self.out)

    def drive(self, runs, now=paladin_backend.utc_now):
        return paladin_backend.drive(runs, status_path=self.status_path, case="case",
                                     cli=self.cli, out=self.out, remote_host=None,
                                     poll_seconds=0, now=now, sleep=lambda seconds: None)

    def calls(self):
        if not os.path.isfile(self.log):
            return []
        with open(self.log) as f:
            return [json.loads(line) for line in f]

    def submitted(self):
        return [call[call.index("--label") + 1] for call in self.calls()
                if call[:2] == ["palace", "submit-detailed"]]

    def status(self):
        with open(self.status_path) as f:
            return json.load(f)

    def test_runs_in_order_and_fetches(self):
        self.script()
        status = self.drive(self.runs())
        self.assertEqual(self.submitted(), ["w", "r", "local-edge"])
        # Each run ends before the next submit: submit, describe (running), describe (end).
        self.assertEqual([" ".join(call[:2]) for call in self.calls()],
                         ["palace submit-detailed", "palace describe-run",
                          "palace describe-run"] * 3)
        self.assertEqual(status, self.status())
        self.assertEqual(status["Backend"], "paladin")
        records = status["Stages"]
        self.assertEqual(len({record["RunID"] for record in records.values()}), 3)
        for name, record in records.items():
            self.assertEqual(record["State"], "completed")
            self.assertAlmostEqual(record["NodeHours"], 600.0 / 3600.0)
            self.assertTrue(os.path.isfile(os.path.join(
                record["Fetched"], "postpro", "domain-response-matrix.csv")))
            self.assertEqual((record["Log"]["H1"], record["Log"]["PCG"]), (1234, [7]))
        with open(os.path.join(records["r"]["Fetched"], "postpro", "run-config.json")) as f:
            self.assertEqual(json.load(f)["Solver"]["Electrostatic"]["ResponseArchive"],
                             {"Mode": "Reduce", "Directory": ARCHIVE, "BlockSize": 48})
        count = len(self.calls())
        self.drive(self.runs())
        self.assertEqual(len(self.calls()), count)

    def test_a_failed_run_stops_and_is_never_resubmitted(self):
        self.script({"w": "failed"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(self.runs())
        records = self.status()["Stages"]
        self.assertEqual(records["w"]["State"], "failed")
        self.assertNotIn("r", records)
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(self.runs())
        self.assertEqual(self.submitted(), ["w"])

    def test_a_crash_in_the_wait_keeps_the_run_id(self):
        self.script()

        def killed(seconds):
            raise RuntimeError("the driver was killed in its wait")

        with self.assertRaises(RuntimeError):
            paladin_backend.drive(self.runs(), status_path=self.status_path, case="case",
                                  cli=self.cli, out=self.out, remote_host=None,
                                  poll_seconds=0, sleep=killed)
        first = self.status()["Stages"]["w"]
        self.assertNotIn("State", first)
        status = self.drive(self.runs())
        self.assertEqual(self.submitted(), ["w", "r", "local-edge"])
        self.assertEqual(status["Stages"]["w"]["RunID"], first["RunID"])
        self.assertEqual(status["Stages"]["w"]["State"], "completed")

    def test_the_cap_cancels(self):
        self.script({"w": "hang"})
        later = lambda: paladin_backend.utc_now() + datetime.timedelta(seconds=10**4)
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(self.runs(), now=later)
        self.assertEqual(self.status()["Stages"]["w"]["State"], "timed-out")
        self.assertIn(["palace", "cancel", str(self.status()["Stages"]["w"]["RunID"]),
                       "--json"], self.calls())
        self.assertEqual(self.submitted(), ["w"])

    def test_the_reducer_reads_the_workers_partition(self):
        self.script()
        runs = self.runs()
        runs[1]["RanksPerNode"] = 96
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(runs)
        self.assertEqual(self.submitted(), ["w"])

    def test_a_refused_submit_records_no_run(self):
        self.script({"w": "refuse"})
        with self.assertRaises(paladin_backend.PaladinError):
            self.drive(self.runs())
        self.assertFalse(os.path.exists(self.status_path))

    def test_main_submit(self):
        self.script({"local-edge": "failed"})
        plan = os.path.join(self.root, "plan.json")
        with open(plan, "w") as f:
            json.dump(three_stage_plan(self.root), f)
        arguments = [plan, "--submit", "--selector", "dev-x", "--lane", "chainport",
                     "--out", self.out, "--paladin-cli", self.cli, "--remote-host", "",
                     "--poll-seconds", "0"]
        with mock.patch("sys.stdout"), mock.patch("sys.stderr"):
            self.assertEqual(paladin_backend.main(arguments), 3)
        states = {name: record["State"] for name, record in self.status()["Stages"].items()}
        self.assertEqual(states, {"w": "completed", "r": "completed", "local-edge": "failed"})


@unittest.skipUnless((RECORDS / "plan.json").exists(), "records of 274fe1a150c0 absent")
class PlanOfRecord(unittest.TestCase):
    def test_fab_plan_of_274fe1a150c0(self):
        with tempfile.TemporaryDirectory() as out:
            code = paladin_backend.main([str(RECORDS / "plan.json"), "--dry-run",
                                         "--selector", "dev-simlapointe-sp-pin-03e576ff3e",
                                         "--lane", "chainport", "--out", out,
                                         "--map", f"{REMOTE}={RECORDS}"])
            self.assertEqual(code, 0)
            with open(os.path.join(out, "paladin-runs.json")) as f:
                runs = json.load(f)
        self.assertEqual(len(runs), 7)
        modes = [run["ResponseArchive"].get("Mode") for run in runs]
        self.assertEqual(modes, ["Archive", "Reduce"] * 3 + [None])
        for worker, reducer in zip(runs[0:6:2], runs[1:6:2]):
            self.assertEqual(worker["ResponseArchive"]["Directory"],
                             reducer["ResponseArchive"]["Directory"])
            self.assertEqual(reducer["Requires"], [worker["Stage"]])
            self.assertTrue(worker["ResponseArchive"]["SourceTiming"])
            self.assertEqual((worker["Nodes"], worker["RanksPerNode"]), (1, 192))


if __name__ == "__main__":
    unittest.main()
