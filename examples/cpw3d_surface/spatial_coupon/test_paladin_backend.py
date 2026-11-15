# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (USER decision 1674 (2); Y-SPEC s2), first piece:
qualify/paladin_backend.py turns a run_stages plan into paladin runs (one per stage) whose
configurations carry the stage variables as Solver.Electrostatic.ResponseArchive keys.
Checked here: the variable -> key table and its refusals, the run configurations, the submit
command, the Requires order and the partition check; and, when the records of key
274fe1a150c0 are present on this host, its fab plan of record."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

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
