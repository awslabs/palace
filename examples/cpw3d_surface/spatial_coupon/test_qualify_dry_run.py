# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""coupon-library qualify --dry-run end to end on the four-edge case with the recorded
physics-11 configuration (controls, stage prefix) and on the gallery-06 case: the
generated configs equal the recorded worker / reducer / local-edge configs apart from
paths, the plan equals the recorded plan in stages and pins; the run config of every
gallery case derived from the case's own sources equals the reference's config apart
from Mesh / Output / DataFile directory (and the recipe-bound Order / Tol); the analysis
of the recorded results through the same records gives the recorded verdicts
(gallery-10: the reference order p5 added as a main stage and gated, p_SA not
applicable); --reference none plans the coupon on its own inputs; the concurrent job
scheduler against a fake remote replaying the recorded trees; the per-coupon fail-closed
stops.  Needs a local identity mesh of each case and the assessment tree."""
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
import build_plan  # noqa: E402
import case_inputs  # noqa: E402
import compare_split_matrices  # noqa: E402
import estimate_stages  # noqa: E402
import gates  # noqa: E402
import job_split  # noqa: E402
import qualify_library  # noqa: E402
import summarize_cost  # noqa: E402
from mixed_mesh import h1_dofs_from_counts  # noqa: E402
from run_gmsh_only_matrix import sha256  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
# A read-only mirror of graded_v2 inputs the assessment tree does not hold (the ten-edge input 09).
REFERENCE_MIRROR = Path(os.environ.get("COUPON_REFERENCE_MIRROR", "/tmp/library-device-path-reference"))
MANIFEST = HERE / "geometry-independence-suite.json"
# The five gallery cases and the config their graded_v2 reference ran (worker.json of the
# recorded campaign, or the producer's spatial_fabricated.json of the mirrored inputs).
GALLERY_REFERENCE_CONFIGS = {
    "four-edge-9d2cb9bbb3fe": ASSESSMENT / "four-edge-physics-13" / "reference" / "case-07-fabricated" / "worker.json",
    "three-edge-419576fdab24": ASSESSMENT / "gallery-physics-06b" / "reference" / "case-06-fabricated" / "worker.json",
    "two-edge-8dd4bc70f183": ASSESSMENT / "gallery-physics-10" / "reference" / "case-10-fabricated" / "worker.json",
    "two-edge-3f8992613e95": ASSESSMENT / "gallery-physics-05" / "reference" / "case-05-fabricated" / "worker.json",
    "ten-edge-6791f1c84123": REFERENCE_MIRROR / "inputs-09" / "spatial_fabricated.json",
}
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


class MonitorTransportFailureTest(unittest.TestCase):
    """A lost ssh connection during monitoring is a transport failure, not "the job left the
    queue" (the 2026-09-21 device run: a poll that returned nothing was read as done, the
    fetch's rsync failed and the driver crashed with the second job still running)."""

    def test_unreachable_poll_keeps_the_job_active(self):
        remote_side = qualify_library.remote_side
        saved = remote_side.ssh
        try:
            remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 255, stdout="", stderr="timed out")
            unreachable = remote_side.poll("h", "/pbs", "1.h", "/r/case/main/status.json")
            remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(
                ["ssh"], 1, stdout="---STATUS---\n\n" + remote_side.POLL_MARKER + "\n", stderr="")
            finished = remote_side.poll("h", "/pbs", "1.h", "/r/case/main/status.json")
            remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(
                ["ssh"], 0, stdout="job_state = R\n---STATUS---\n{\"State\": \"running\", \"Stages\": []}\n" + remote_side.POLL_MARKER + "\n", stderr="")
            running = remote_side.poll("h", "/pbs", "1.h", "/r/case/main/status.json")
        finally:
            remote_side.ssh = saved
        self.assertFalse(unreachable["Reachable"])
        self.assertEqual((unreachable["JobState"], unreachable["SSHReturnCode"]), (None, 255))
        self.assertTrue(finished["Reachable"])
        self.assertIsNone(finished["JobState"])          # reachable and not in the queue: done
        self.assertEqual((running["Reachable"], running["JobState"], running["Status"]["State"]), (True, "R", "running"))
        record = {"Case": "c", "Submission": {"Job": "1.h"}, "Remote": {"Case": "/r/case"},
                  "Monitor": {"Polls": 0, "LastJobState": "R", "LastStages": None}}
        logged = []
        saved_poll = remote_side.poll
        try:
            remote_side.poll = lambda *args: dict(unreachable)
            done_unreachable = qualify_library.poll_case(record, remote={"Host": "h"}, profile={"PBSBin": "/pbs"}, log=logged.append)
            state_after_unreachable = record["Monitor"]["LastJobState"]
            remote_side.poll = lambda *args: dict(finished)
            done_finished = qualify_library.poll_case(record, remote={"Host": "h"}, profile={"PBSBin": "/pbs"}, log=logged.append)
        finally:
            remote_side.poll = saved_poll
        self.assertFalse(done_unreachable)
        self.assertEqual(record["Monitor"]["TransportFailures"], 1)
        self.assertEqual(state_after_unreachable, "R")     # the unreachable poll changes no job state
        self.assertIn("unreachable", logged[0])
        self.assertTrue(done_finished)
        self.assertEqual(record["Monitor"]["Polls"], 2)

    def test_held_job_is_still_in_the_queue(self):
        """PBS H (the SOCA dispatcher's capacity retry: 'Job held by root', error_message
        CF:ROLLBACK_COMPLETE:retry=1, released at retry_eligible_after) is not "left the queue"
        (the 2026-09-21 split acceptance: the reducer job 46339 was held 4 min after its
        submission and the driver stopped the coupon with 'no status.json fetched')."""
        remote_side = qualify_library.remote_side
        self.assertEqual([remote_side.in_queue(state) for state in ("Q", "R", "E", "H", "W", "T", "S", "B", "F", "X", None)],
                         [True] * 8 + [False] * 3)
        held = {"UTC": "u", "JobState": "H", "Reachable": True, "SSHReturnCode": 0, "Status": None,
                "QStat": "job_state = H\nHold_Types = u\ncomment = Job held by root on Mon Sep 21 21:17:12 2026\n"
                         "Resource_List.error_message = CF:ROLLBACK_COMPLETE:retry=1"}
        record = {"Case": "c", "Submission": {"Job": "1.h"}, "Remote": {"Case": "/r/case"},
                  "Monitor": {"Polls": 0, "LastJobState": "Q", "LastStages": None}}
        logged = []
        saved_poll = remote_side.poll
        try:
            remote_side.poll = lambda *args: dict(held)
            done = qualify_library.poll_case(record, remote={"Host": "h"}, profile={"PBSBin": "/pbs"}, log=logged.append)
        finally:
            remote_side.poll = saved_poll
        self.assertFalse(done)
        self.assertEqual(record["Monitor"]["LastJobState"], "H")
        self.assertIn("held, still queued", logged[0])
        self.assertIn("Job held by root", logged[0])
        commands = []
        saved_ssh = remote_side.ssh
        try:
            remote_side.ssh = lambda host, command, check=True: (commands.append(command) or subprocess.CompletedProcess(
                ["ssh"], 0, stdout="job_state = H\nHold_Types = u\ncomment = Job held by root\nResource_List.error_message = CF:ROLLBACK_COMPLETE:retry=1\n"
                                   "---STATUS---\n\n" + remote_side.POLL_MARKER + "\n", stderr=""))
            polled = remote_side.poll("h", "/pbs", "1.h", "/r/case/main/status.json")
        finally:
            remote_side.ssh = saved_ssh
        self.assertIn("Hold_Types|error_message", commands[0])          # the poll's qstat grep carries the hold and its reason
        self.assertEqual(polled["JobState"], "H")
        self.assertIn("CF:ROLLBACK_COMPLETE", polled["QStat"])

    def test_fetch_transport_failure_is_a_recorded_stop(self):
        remote_side = qualify_library.remote_side
        saved = remote_side.fetch

        def failing_fetch(host, remote_directory, local_directory):
            raise subprocess.CalledProcessError(255, ["rsync", remote_directory, str(local_directory)])
        tmp = Path(tempfile.mkdtemp(prefix="qualify-fetch-stop-"))
        record = {"Case": "c", "Root": str(tmp), "Remote": {"Case": "/r/case"}, "Submission": {"Job": "1.h"}}
        context = {"jobs": [{"Name": "single", "Kind": "single", "RemoteDirectory": "/r/case/main", "Submission": {"Job": "1.h"}}]}
        try:
            remote_side.fetch = failing_fetch
            with self.assertRaises(qualify_library.CaseStop) as stop:
                qualify_library.finish_case(record, context, remote={"Host": "h"}, profile={"PBSBin": "/pbs"})
        finally:
            remote_side.fetch = saved
            shutil.rmtree(tmp, True)
        self.assertEqual(stop.exception.record["Kind"], "Fetch")
        self.assertIn("--resume", stop.exception.record["Message"])
        self.assertEqual(stop.exception.record["Command"][0], "rsync")
        self.assertNotIn("Fetch", record)


class JobSplitTest(unittest.TestCase):
    """The per-coupon source split policy (decision 61b) on the recorded device coupon
    spatial-3-edge-5d3b5e644745 (225 sources, H1 42.7M at p4: 29,272 s single-job estimate
    against the 21,600 s walltime - the coupon the first device library run failed closed)."""

    @classmethod
    def setUpClass(cls):
        build = json.loads((HERE / "qualify" / "device-library-20260921" / "library-build.json").read_text())
        cls.counts = next(case for case in build["Cases"] if case["Case"] == "spatial-3-edge-5d3b5e644745")["H1"]["EntityCounts"]
        cls.model = estimate_stages.load_cost_model()
        cls.profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        cls.indices = list(range(1, 226))
        cls.layout = qualify_library.stage_layout("c", [4], [3, 5], 225, 8)
        stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in cls.layout if item["Kind"] == "response"]
        local = next(item for item in cls.layout if item["Kind"] == "local-edge")
        # The recorded 7f03 estimate and split (decision 61b) ran at the cost model's
        # measured block size 6 with the pair-scaled reducer formula (before decision
        # 62(1) split the block-pair part into evaluation + Gram): reproduced as recorded.
        cls.estimate = estimate_stages.estimate(cls.counts, stages, model={**cls.model, "ReducerEvaluationFraction": 0.0},
                                                profile=cls.profile, block_size=6,
                                                local_edge=(local["EstimateKey"], local["Order"], local["Sources"]))
        cls.estimate_default = estimate_stages.estimate(cls.counts, stages, model=cls.model, profile=cls.profile,
                                                        local_edge=(local["EstimateKey"], local["Order"], local["Sources"]))

    def split(self, mode, max_jobs, fixed=None):
        policy = job_split.normalize_policy(mode, max_jobs=max_jobs, walltime_seconds=self.profile["WalltimeSeconds"],
                                            fixed_jobs=fixed, user_job_cap=self.profile["UserJobCap"])
        return job_split.plan_split(indices=self.indices, layout=self.layout, estimate=self.estimate, policy=policy,
                                    model=self.model, profile=self.profile)

    def test_single_job_estimate_is_the_recorded_fail_closed_one(self):
        self.assertFalse(self.estimate["FitsOneJob"])
        self.assertAlmostEqual(self.estimate["JobSecondsEstimateWithPreflightAndMargin"]["2.0"], 29271.5, delta=1.0)
        # At the decision-62(1) default (b = 48) the reducer estimate shrinks (its evaluation
        # part 38 -> 5 block rows) but 7f03 still does not fit one job: fail closed unchanged.
        self.assertEqual(self.estimate_default["ReducerBlockSize"], 48)
        self.assertLess(self.estimate_default["Stages"]["p4-225"]["ReducerSecondsEstimate"],
                        self.estimate["Stages"]["p4-225"]["ReducerSecondsEstimate"])
        self.assertFalse(self.estimate_default["FitsOneJob"])
        self.assertGreater(self.estimate_default["JobSecondsEstimateWithPreflightAndMargin"]["2.0"], self.profile["WalltimeSeconds"])
        single = self.split("fixed", 4, fixed=1)
        self.assertFalse(single["Fits"])
        self.assertIsNone(single["N"])
        self.assertIn("even the maximal split N = 1", single["Decision"])
        self.assertAlmostEqual(single["Candidates"][0]["LongestJobSeconds"],
                               self.estimate["JobSecondsEstimateWithPreflightAndMargin"]["2.0"], places=6)

    def test_frugal_is_the_fewest_jobs_that_fit(self):
        record = self.split("frugal", 4)
        self.assertEqual(record["N"], 2)
        self.assertTrue(record["Fits"])
        self.assertEqual(sum(len(block) for block in record["Blocks"]), 225)
        self.assertEqual([index for block in record["Blocks"] for index in block], self.indices)   # contiguous, in order
        self.assertLess(len(record["Blocks"][0]), len(record["Blocks"][1]))   # job 1 carries the controls + local-edge
        self.assertEqual([job["Kind"] for job in record["Jobs"]], ["worker", "worker", "reducer"])
        worst = record["WorstPCGFactor"]
        for job in record["Jobs"]:
            self.assertLess(job["SecondsEstimateWithPreflightAndMargin"][worst], self.profile["WalltimeSeconds"])
        # Balanced: the two worker jobs end within one source's cost of each other.
        per_source = self.estimate["Stages"]["p4-225"]["ByPCGFactor"][worst]["PerSourceSecondsEstimate"] * self.model["PreflightAndMarginFactor"]
        w1, w2 = (job["SecondsEstimateWithPreflightAndMargin"][worst] for job in record["Jobs"][:2])
        self.assertLess(abs(w1 - w2), per_source)
        self.assertEqual(record["ControlsJob"], "worker-1")
        self.assertEqual([item["N"] for item in record["Candidates"]], [1, 2, 3, 4])
        self.assertEqual([item["Fits"] for item in record["Candidates"]], [False, True, True, True])

    def test_speed_minimizes_the_estimated_critical_path_within_max_jobs(self):
        record = self.split("speed", 4)
        self.assertEqual(record["N"], 4)
        self.assertEqual(record["Blocks"], [self.indices[:17], self.indices[17:87], self.indices[87:156], self.indices[156:]])
        worst = record["WorstPCGFactor"]
        paths = [item["CriticalPathSeconds"] for item in record["Candidates"]]
        self.assertEqual(paths, sorted(paths, reverse=True))   # more jobs, shorter path (the reducer job is the floor)
        self.assertEqual(record["CriticalPathEstimateSeconds"][worst], min(paths))
        reducer = record["Jobs"][-1]
        self.assertEqual(reducer["Kind"], "reducer")
        self.assertEqual(reducer["Sources"], self.indices)
        self.assertAlmostEqual(record["CriticalPathEstimateSeconds"][worst],
                               max(job["SecondsEstimateWithPreflightAndMargin"][worst] for job in record["Jobs"][:-1])
                               + reducer["SecondsEstimateWithPreflightAndMargin"][worst])
        # Node time grows with N (one preflight and non-source setup per job): recorded per candidate.
        nodes = [item["NodeSeconds"] for item in record["Candidates"]]
        self.assertEqual(nodes, sorted(nodes))
        # speed never uses more jobs than shorten the path: with a walltime the reducer job
        # dominates, two candidates tie and the smaller N is chosen.
        capped = self.split("speed", 2)
        self.assertEqual(capped["N"], 2)

    def test_fixed_and_the_fail_closed_bounds(self):
        record = self.split("fixed", 4, fixed=3)
        self.assertEqual((record["N"], len(record["Jobs"])), (3, 4))
        with self.assertRaisesRegex(ValueError, "exceeds the largest split"):
            self.split("fixed", 4, fixed=5)
        with self.assertRaisesRegex(ValueError, "FixedJobs"):
            job_split.normalize_policy("fixed", max_jobs=4, walltime_seconds=1.0)
        with self.assertRaisesRegex(ValueError, "not one of"):
            job_split.normalize_policy("fast", max_jobs=4, walltime_seconds=1.0)
        # Even the maximal split does not fit a short walltime: fail closed with the reason.
        policy = job_split.normalize_policy("speed", max_jobs=4, walltime_seconds=3600.0, user_job_cap=40)
        short = job_split.plan_split(indices=self.indices, layout=self.layout, estimate=self.estimate, policy=policy,
                                     model=self.model, profile=self.profile)
        self.assertFalse(short["Fits"])
        self.assertIn("even the maximal split N = 4", short["Decision"])
        self.assertIn("fail closed", short["Decision"])

    def test_controls_take_a_separate_first_job_when_they_leave_no_room(self):
        # 8 sources split in 4: the controls + local-edge alone outweigh a 2-source block, so
        # job 1 carries them alone and the blocks go to jobs 2..4 (recorded).
        indices = list(range(1, 9))
        layout = qualify_library.stage_layout("c", [4], [3, 5], 8, 8)
        stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in layout if item["Kind"] == "response"]
        local = next(item for item in layout if item["Kind"] == "local-edge")
        estimate = estimate_stages.estimate(self.counts, stages, model=self.model, profile=self.profile,
                                            local_edge=(local["EstimateKey"], local["Order"], local["Sources"]))
        policy = job_split.normalize_policy("fixed", max_jobs=4, walltime_seconds=self.profile["WalltimeSeconds"], fixed_jobs=4)
        record = job_split.plan_split(indices=indices, layout=layout, estimate=estimate, policy=policy, model=self.model,
                                      profile=self.profile)
        self.assertEqual(record["ControlsJob"], "separate")
        self.assertEqual(record["Blocks"][0], [])
        self.assertEqual(sum(len(block) for block in record["Blocks"]), 8)
        self.assertEqual([len(block) for block in record["Blocks"][1:]], [3, 3, 2])

    def test_block_plan_pins_and_caps(self):
        record = self.split("speed", 4)
        digests = {"c-p4": {f"worker-block{k}.json": f"{k}" * 64 for k in range(1, 5)} | {"reducer.json": "r" * 64},
                   "c-p5-control": {"worker.json": "a" * 64, "reducer.json": "b" * 64},
                   "c-p3-control": {"worker.json": "c" * 64, "reducer.json": "d" * 64},
                   "c-p4-local-edge": {"config.json": "e" * 64}}
        traces = {f"/r/case/inputs/traces/basis-{i:04d}.csv": "t" * 64 for i in self.indices}
        common = dict(case_id="c", remote_case_root="/r/case", mesh={"Remote": "/r/case/mesh/m.msh", "SHA256": "m" * 64, "Local": "/l/m.msh"},
                      stage_layout=self.layout, estimate=self.estimate, config_digests=digests, trace_pins=traces,
                      profile=self.profile, binary="/r/b.bin", binary_sha256="0" * 64, mpiexec="/r/mpiexec_bound.sh",
                      purpose="test", factors=["1.0", "1.5", "2.0"])
        plans = {job["Name"]: build_plan.build_job_plan(job_name=job["Name"], split_job=job, **common) for job in record["Jobs"]}
        self.assertEqual(plans["worker-1"]["StageNames"], ["c-p4-worker-block1", "c-p5-control-worker", "c-p5-control-reducer",
                                                           "c-p3-control-worker", "c-p3-control-reducer", "c-p4-local-edge"])
        self.assertEqual(plans["worker-2"]["StageNames"], ["c-p4-worker-block2"])
        self.assertEqual(plans["reducer"]["StageNames"], ["c-p4-reducer"])
        self.assertEqual(plans["reducer"]["Stages"][0]["Requires"], [])
        archive = "/r/case/main/c-p4/archive"
        for name in ("worker-1", "worker-2", "worker-3", "worker-4"):
            block = plans[name]["Stages"][0]
            self.assertEqual(block["Environment"]["PALACE_RESPONSE_ARCHIVE_DIR"], archive)
            self.assertEqual(block["Environment"]["PALACE_RESPONSE_ARCHIVE_ONLY"], "1")
            self.assertEqual(block["Config"], f"/r/case/main/c-p4/worker-block{name[-1]}.json")
            self.assertGreaterEqual(block["CapSeconds"], block["MinimumSeconds"])
            self.assertEqual(len(plans[name]["PinnedSHA256"]), 1 + 225 + (6 if name == "worker-1" else 1))
        self.assertEqual(plans["reducer"]["Stages"][0]["Environment"]["PALACE_RESPONSE_ARCHIVE_DIR"], archive)
        self.assertEqual(len(plans["reducer"]["PinnedSHA256"]), 1 + 225 + 1)
        # A block's cap follows its own source count: block 2 (70 sources) above block 1 (17).
        self.assertGreater(plans["worker-2"]["Stages"][0]["CapSeconds"], plans["worker-1"]["Stages"][0]["CapSeconds"])
        script = build_plan.render_job_script(profile=self.profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                              job_name="j", walltime_seconds=21600, job_directory="/r/case/main/jobs/worker-2")
        self.assertIn("D=/r/case/main/jobs/worker-2\n", script)
        self.assertIn("#PBS -o /r/case/main/jobs/worker-2/pbs.log", script)
        single = build_plan.render_job_script(profile=self.profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                              job_name="j", walltime_seconds=21600)
        self.assertIn("D=/r/case/main\n", single)

    def test_manifest_job_policy_default_and_command_line_override(self):
        """The production recipe records the default policy (frugal) under
        PhysicsRun.JobPolicy; case_inputs and general_mesh_manifest validate it; the
        command line overrides it and the origin is recorded."""
        from general_mesh_manifest import validate_physics_run
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        physics_run = case_inputs.physics_run_parameters(manifest)
        self.assertEqual(physics_run["JobPolicy"]["Mode"], "frugal")
        self.assertIsNone(physics_run["JobPolicy"]["FixedJobs"])
        self.assertIn("61b", physics_run["JobPolicy"]["Rule"])
        validate_physics_run(manifest["ProductionRecipe"])
        for broken in ({"Mode": "fast", "Rule": "x"}, {"Mode": "fixed", "Rule": "x"}, {"Mode": "frugal", "FixedJobs": 2, "Rule": "x"},
                       {"Mode": "frugal"}, "frugal"):
            recipe = json.loads(json.dumps(manifest["ProductionRecipe"]))
            recipe["PhysicsRun"]["JobPolicy"] = broken
            with self.assertRaisesRegex(ValueError, "JobPolicy"):
                validate_physics_run(recipe)
            with self.assertRaisesRegex(case_inputs.CaseInputError, "JobPolicy"):
                case_inputs.physics_run_parameters({**manifest, "ProductionRecipe": recipe})
        recipe = json.loads(json.dumps(manifest["ProductionRecipe"]))
        recipe["PhysicsRun"]["JobPolicy"] = {"Mode": "fixed", "FixedJobs": 3, "Rule": "x"}
        self.assertEqual(case_inputs.physics_run_parameters({**manifest, "ProductionRecipe": recipe})["JobPolicy"],
                         {"Mode": "fixed", "FixedJobs": 3, "Rule": "x"})
        args = argparse.Namespace(job_policy=None, fixed_jobs=None, max_jobs=4)
        policy = qualify_library.job_policy_of(args, physics_run, self.profile)
        self.assertEqual((policy["Mode"], policy["MaxJobs"], policy["WalltimeSeconds"], policy["UserJobCap"]),
                         ("frugal", 4, self.profile["WalltimeSeconds"], self.profile["UserJobCap"]))
        self.assertIn("manifest", policy["Origin"])
        policy = qualify_library.job_policy_of(argparse.Namespace(job_policy="speed", fixed_jobs=None, max_jobs=4), physics_run, self.profile)
        self.assertEqual((policy["Mode"], policy["Origin"]), ("speed", "--job-policy"))
        policy = qualify_library.job_policy_of(argparse.Namespace(job_policy=None, fixed_jobs=None, max_jobs=2), {"Order": 4}, self.profile)
        self.assertEqual((policy["Mode"], policy["Origin"]), ("frugal", "built-in default frugal"))
        with self.assertRaises(qualify_library.CaseStop):
            qualify_library.job_policy_of(argparse.Namespace(job_policy="fixed", fixed_jobs=None, max_jobs=2), physics_run, self.profile)

    def test_manifest_reducer_block_size_default_and_command_line_override(self):
        """Decision 62(1): the production recipe records PhysicsRun.ReducerBlockSize 48
        (PreviousValue 6); case_inputs and general_mesh_manifest validate it; the command
        line overrides it and the origin is recorded; the plan's reducer stages carry it."""
        from general_mesh_manifest import validate_physics_run
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        block = manifest["ProductionRecipe"]["PhysicsRun"]["ReducerBlockSize"]
        self.assertEqual((block["Value"], block["PreviousValue"]), (48, 6))
        self.assertIn("62(1)", block["Provenance"])
        physics_run = case_inputs.physics_run_parameters(manifest)
        self.assertEqual(physics_run["ReducerBlockSize"], 48)
        validate_physics_run(manifest["ProductionRecipe"])
        for broken in ({"Value": 0, "Rule": "x"}, {"Value": 6}, {"Value": "6", "Rule": "x"}, {"Value": True, "Rule": "x"}, 6):
            recipe = json.loads(json.dumps(manifest["ProductionRecipe"]))
            recipe["PhysicsRun"]["ReducerBlockSize"] = broken
            with self.assertRaisesRegex(ValueError, "ReducerBlockSize"):
                validate_physics_run(recipe)
            with self.assertRaisesRegex(case_inputs.CaseInputError, "ReducerBlockSize"):
                case_inputs.physics_run_parameters({**manifest, "ProductionRecipe": recipe})
        chosen = qualify_library.reducer_block_size_of(argparse.Namespace(reducer_block_size=None), physics_run)
        self.assertEqual((chosen["Value"], chosen["Origin"]), (48, "manifest ProductionRecipe.PhysicsRun.ReducerBlockSize"))
        chosen = qualify_library.reducer_block_size_of(argparse.Namespace(reducer_block_size=6), physics_run)
        self.assertEqual((chosen["Value"], chosen["Origin"]), (6, "--reducer-block-size"))
        chosen = qualify_library.reducer_block_size_of(argparse.Namespace(reducer_block_size=None), {"Order": 4})
        self.assertEqual((chosen["Value"], chosen["Origin"]), (48, "built-in default build_plan.DEFAULT_REDUCER_BLOCK_SIZE"))
        self.assertIn("14 MB", chosen["Rule"])
        with self.assertRaises(qualify_library.CaseStop):
            qualify_library.reducer_block_size_of(argparse.Namespace(reducer_block_size=0), physics_run)
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        self.assertIsNone(parser.parse_args(["--build-record", "b", "--reference", "none", "--remote", "h:/r", "--frozen-binary-sha256", "f"]).reducer_block_size)
        self.assertEqual(parser.parse_args(["--build-record", "b", "--reference", "none", "--remote", "h:/r", "--frozen-binary-sha256", "f",
                                            "--reducer-block-size", "12"]).reducer_block_size, 12)

    def test_compare_split_matrices_reports_roundoff_and_differences(self):
        tmp = Path(tempfile.mkdtemp(prefix="split-compare-"))
        try:
            for name, rows in (("a", [("1", "1", "+1.000000000000e-17"), ("1", "2", "+2.000000000000e-18"), ("2", "2", "+3.000000000000e-17")]),
                               ("b", [("1", "2", "+2.000000000000e-18"), ("2", "2", "+3.000000000001e-17"), ("1", "1", "+1.000000000000e-17")]),
                               ("c", [("1", "1", "+1.000000000000e-17"), ("1", "2", "+2.000000000000e-18"), ("2", "2", "+3.300000000000e-17")])):
                (tmp / name).mkdir()
                (tmp / name / "domain-response-matrix.csv").write_text(
                    "  basis_i,  basis_j,                   Q_ij (J)\n" + "".join(f" {i}.00e+00, {j}.00e+00, {q}\n" for i, j, q in rows))
                (tmp / name / "surface-response-matrix.csv").write_text(
                    "interface, edge, R (m), basis_i, basis_j, Q_ij (J)\n" + "".join(f" 1.00e+00, 1.00e+00, +2.0e-06, {i}.00e+00, {j}.00e+00, {q}\n" for i, j, q in rows))
            same = compare_split_matrices.compare(tmp / "a", tmp / "b")
            self.assertTrue(same["Equal"])
            self.assertAlmostEqual(same["MaxRelativeToLargestEntry"], 1e-12 / 3, delta=1e-14)   # roundoff in the last digit, row order free
            self.assertEqual(same["Matrices"]["domain"]["Columns"]["Q_ij (J)"]["ExactText"], 2)
            different = compare_split_matrices.compare(tmp / "a", tmp / "c")
            self.assertFalse(different["Equal"])
            self.assertAlmostEqual(different["MaxRelativeDifference"], 0.3e-17 / 3.3e-17, places=6)   # |a - b| / max(|a|, |b|)
            self.assertEqual(different["Matrices"]["surface"]["Columns"]["Q_ij (J)"]["Worst"]["Key"], (1, 1, 2, 2))
            per_source = different["Matrices"]["domain"]["Columns"]["Q_ij (J)"]["PerSourceMaxRelativeDifference"]
            self.assertEqual(list(per_source), ["1", "2"])
            self.assertEqual(per_source["1"], 0.0)
            self.assertAlmostEqual(per_source["2"], 0.3e-17 / 3.3e-17, places=6)
            with self.assertRaisesRegex(ValueError, "different row keys"):
                (tmp / "d").mkdir()
                for name in ("domain-response-matrix.csv", "surface-response-matrix.csv"):
                    text = (tmp / "a" / name).read_text().splitlines()
                    (tmp / "d" / name).write_text("\n".join(text[:-1]) + "\n")
                compare_split_matrices.compare(tmp / "a", tmp / "d")
        finally:
            shutil.rmtree(tmp, True)

    def test_merged_split_statuses_read_as_one_stage(self):
        def timing(indices):
            return [{"Index": i, "Iterations": 20, "SolveSeconds": 10.0, "TotalSeconds": 12.0} for i in indices]
        statuses = {
            "worker-1": {"PBSJobID": "1.h", "Host": "n1", "StartUTC": "2026-09-21T10:00:00Z", "EndUTC": "2026-09-21T10:30:00Z",
                         "TotalSeconds": 1800.0, "State": "complete",
                         "Stages": [{"Name": "c-p4-worker-block1", "State": "complete", "WallSeconds": 100.0, "NodePeakUsedBytesSampled": 5,
                                     "MaxSingleProcessRSSBytes": 7,
                                     "Parsed": {"Order": 4, "H1": 100, "SourceTiming": timing([1, 2]), "PCG": [20, 20], "Nonconvergence": [],
                                                "PalaceTotalSeconds": 60.0, "PalacePeakMemory": {"Total": "10.0G"}}},
                                    {"Name": "c-p3-control-worker", "State": "complete", "WallSeconds": 5.0,
                                     "Parsed": {"SourceTiming": timing([1]), "PCG": [20], "Nonconvergence": [], "PalaceTotalSeconds": 4.0}},
                                    {"Name": "c-p3-control-reducer", "State": "complete", "WallSeconds": 3.0,
                                     "Parsed": {"SourceTiming": [], "PCG": [], "Nonconvergence": [], "PalaceTotalSeconds": 2.0}},
                                    {"Name": "c-p4-local-edge", "State": "complete", "WallSeconds": 9.0, "Parsed": {"PCG": [1], "Nonconvergence": []}}]},
            "worker-2": {"PBSJobID": "2.h", "Host": "n2", "StartUTC": "2026-09-21T10:00:00Z", "EndUTC": "2026-09-21T10:40:00Z",
                         "TotalSeconds": 2400.0, "State": "complete",
                         "Stages": [{"Name": "c-p4-worker-block2", "State": "complete", "WallSeconds": 200.0, "NodePeakUsedBytesSampled": 9,
                                     "MaxSingleProcessRSSBytes": 3,
                                     "Parsed": {"Order": 4, "H1": 100, "SourceTiming": timing([3, 4, 5]), "PCG": [20, 20, 20], "Nonconvergence": [],
                                                "PalaceTotalSeconds": 96.0, "PalacePeakMemory": {"Total": "12.0G"}}}]},
            "reducer": {"PBSJobID": "3.h", "Host": "n3", "StartUTC": "2026-09-21T11:00:00Z", "EndUTC": "2026-09-21T11:20:00Z",
                        "TotalSeconds": 1200.0, "State": "complete",
                        "Stages": [{"Name": "c-p4-reducer", "State": "complete", "WallSeconds": 500.0,
                                    "Parsed": {"SourceTiming": [], "PCG": [], "Nonconvergence": [], "PalaceTotalSeconds": 450.0}}]}}
        jobs = [{"Name": "worker-1", "Kind": "worker"}, {"Name": "worker-2", "Kind": "worker"}, {"Name": "reducer", "Kind": "reducer"}]
        merged = summarize_cost.merge_split_statuses(statuses, jobs)
        self.assertEqual(merged["TotalSeconds"], 5400.0)
        self.assertEqual((merged["StartUTC"], merged["EndUTC"]), ("2026-09-21T10:00:00Z", "2026-09-21T11:20:00Z"))
        names = [stage["Name"] for stage in merged["Stages"]]
        self.assertEqual(sorted(names), sorted(["c-p4-worker", "c-p4-reducer", "c-p3-control-worker", "c-p3-control-reducer", "c-p4-local-edge"]))
        worker = next(stage for stage in merged["Stages"] if stage["Name"] == "c-p4-worker")
        self.assertEqual(worker["WallSeconds"], 300.0)
        self.assertEqual([t["Index"] for t in worker["Parsed"]["SourceTiming"]], [1, 2, 3, 4, 5])
        self.assertEqual(worker["Parsed"]["PalaceTotalSeconds"], 156.0)
        self.assertEqual(worker["Parsed"]["PalacePeakMemory"]["Total"], "12.0G")
        self.assertEqual((worker["NodePeakUsedBytesSampled"], worker["MaxSingleProcessRSSBytes"], worker["Blocks"]), (9, 7, 2))
        summary = summarize_cost.summarize(merged, full_sources=5, nodes=1)
        self.assertAlmostEqual(summary["JobNodeHours"], 1.5)
        self.assertEqual(summary["Stages"]["c-p4"]["Sources"], [1, 2, 3, 4, 5])
        self.assertEqual(summary["Stages"]["c-p4"]["StageWallSeconds"], 800.0)
        self.assertAlmostEqual(summary["Stages"]["c-p4"]["WorkerNonSourceSeconds"], 156.0 - 60.0)
        self.assertEqual({name: job["NodeHours"] for name, job in summary["Jobs"].items()}, {"worker-1": 0.5, "worker-2": 2400 / 3600, "reducer": 1200 / 3600})


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
        # The recorded campaigns reduced at PALACE_RESPONSE_BLOCK_SIZE 6 (before decision
        # 62(1)): the command-line override reproduces their plans; the origin is recorded.
        record = self.dry_run(case_id, root, extra=("--reducer-block-size", "6"))
        case = record["Cases"][0]
        self.assertEqual(case["Status"], "planned", case.get("StoppedBy"))
        self.assertIsNone(case["StoppedBy"])
        self.assertEqual((case["ReducerBlockSize"]["Value"], case["ReducerBlockSize"]["Origin"]), (6, "--reducer-block-size"))
        self.assertEqual(record["Library"]["ReducerBlockSize"]["CommandLine"], 6)
        self.assertEqual(case["Estimate"]["ReducerBlockSize"], 6)
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
        # dependencies, order) and pins the same files: every trace name of the recorded plan
        # (the traces are regenerated from the case's basis - their digests are the run's
        # own, recorded under Inputs.Sources); the mesh pin is the build record's.
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
        regenerated = {source["Name"]: source["SHA256"] for source in case["Inputs"]["Sources"]}
        self.assertEqual({name: pins[name] for name in traces}, {name: regenerated[name] for name in traces})
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
                                  stage_prefix=spec["Prefix"], orders=[], controls=[3, 5], frozen_binary_sha256=BINARY_SHA256,
                                  max_jobs=2, job_policy=None, fixed_jobs=None, reducer_block_size=None)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        model = estimate_stages.load_cost_model()
        table, digest = gates.load_gates()
        root = self.tmp / f"analysis-{spec['Prefix']}-{Path(str(reference)).name if reference is not None else 'none'}"
        root.mkdir(exist_ok=True)
        case = next(item for item in build["Cases"] if item["Case"] == case_id)
        record, context = qualify_library.prepare_case(case, manifest_path=MANIFEST, manifest=manifest, args=args, root=root,
                                                       remote={"Host": "h", "Root": "/r"}, profile=profile, cost_model=model,
                                                       gates=table, gates_digest=digest)
        return record, context, root, table, digest, profile, manifest

    @staticmethod
    def snapshot(campaign):
        """mtime of every file of the recorded campaign (results, reference, main, ...)."""
        return {path: path.stat().st_mtime_ns for path in Path(campaign).rglob("*") if path.is_file()}

    def assert_campaign_untouched(self, campaign, before):
        self.assertEqual(self.snapshot(campaign), before, "the recorded campaign directory must stay read only")

    def test_analysis_of_the_recorded_results_reproduces_the_verdicts(self):
        for case_id, spec in CASES.items():
            campaign = ASSESSMENT / spec["Campaign"]
            before = self.snapshot(campaign)
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
                                  stage_prefix=spec["Prefix"], orders=[], controls=[3, 5], frozen_binary_sha256=BINARY_SHA256,
                                  max_jobs=2, job_policy=None, fixed_jobs=None, reducer_block_size=None)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        model = estimate_stages.load_cost_model()
        table, digest = gates.load_gates()
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        root = self.tmp / "analysis-g10"
        root.mkdir(exist_ok=True)
        before = self.snapshot(campaign)
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
            args = argparse.Namespace(build_record=self.build_record, reference=None, remote="h:/r", orders=[], controls=[3, 5],
                                      control_count=8, control_source=None, max_jobs=2, frozen_binary_sha256=BINARY_SHA256,
                                      stage_prefix=None, case=None, root=self.tmp / "concurrent", dry_run=False, resume=False,
                                      monitor_interval=0, monitor_polls=10, cluster_profile=HERE / "qualify" / "cluster-profile.json",
                                      cost_model=estimate_stages.COST_MODEL, gates=gates.GATES_FILE, job_policy=None, fixed_jobs=None,
                                      merge_into=None, reducer_block_size=None)
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
            # A coupon that fits one job is one job (N = 1, the recorded campaigns' layout).
            self.assertEqual((case["Split"]["N"], case["JobPolicy"]["Mode"], case["Split"]["ControlsJob"]), (1, "frugal", "worker-1"))
            self.assertEqual([job["Kind"] for job in case["Jobs"]], ["single"])
            self.assertEqual(case["Cost"]["Split"]["N"], 1)
            self.assertEqual(case["Cost"]["Jobs"]["single"]["ActualSeconds"], case["Cost"]["JobTotalSeconds"])
            self.assertIsNotNone(case["Cost"]["CriticalPathSeconds"])
        self.assertEqual({case: split["N"] for case, split in totals["Splits"].items()}, {case: 1 for case in CASES})
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

    def test_split_coupon_runs_worker_jobs_then_the_reducer_job(self):
        """--job-policy fixed --fixed-jobs 2 on the four-edge coupon against a fake remote: both
        worker jobs are submitted at once, each is completed from its own status.json when it
        leaves the queue, the archive union is counted before the one reducer job is submitted,
        the coupon is fetched / verified / analyzed after the reducer job, the merged cost reads
        the split as one main stage (per-source timings in block order) and records every job."""
        case_id = "four-edge-9d2cb9bbb3fe"
        spec = CASES[case_id]
        recorded = json.loads((ASSESSMENT / spec["Campaign"] / "results" / "main" / "status.json").read_text())
        events = []
        holder = {}

        def job_status(name):
            blocks = holder["context"]["split"]["Blocks"]
            stages = {stage["Name"]: stage for stage in recorded["Stages"]}
            worker = stages[f"{spec['Prefix']}-p4-worker"]

            def block_stage(k):
                block = set(blocks[k - 1])
                timings = [t for t in worker["Parsed"]["SourceTiming"] if t["Index"] in block]
                positions = [i for i, t in enumerate(worker["Parsed"]["SourceTiming"]) if t["Index"] in block]
                parsed = dict(worker["Parsed"], SourceTiming=timings, PCG=[worker["Parsed"]["PCG"][i] for i in positions],
                              PalaceTotalSeconds=worker["Parsed"]["PalaceTotalSeconds"] * len(block) / 80)
                return dict(worker, Name=f"{spec['Prefix']}-p4-worker-block{k}", Parsed=parsed,
                            WallSeconds=worker["WallSeconds"] * len(block) / 80)
            if name == "worker-1":
                selected = [block_stage(1)] + [stages[n] for n in stages if "control" in n or n.endswith("local-edge")]
            elif name == "worker-2":
                selected = [block_stage(2)]
            else:
                selected = [stages[f"{spec['Prefix']}-p4-reducer"]]
            return dict(recorded, Stages=selected, PBSJobID=f"{name}.fake",
                        TotalSeconds=sum(stage["WallSeconds"] for stage in selected) + 60.0)

        def fake_upload(record, context, *, remote, profile):
            events.append(("upload", record["Case"]))
            return {"Commands": [], "UTC": "fake"}

        def fake_submit(host, pbs_bin, script, cwd, *, job_cap, user=None):
            events.append(("submit", Path(cwd).name))
            return {"Job": f"{Path(cwd).name}.fake", "UTC": qualify_library.remote_side.utc(), "UserJobsBefore": 0, "JobCap": job_cap,
                    "Command": "qsub"}

        def fake_poll(host, pbs_bin, job_id, status_path):
            events.append(("poll", Path(status_path).parts[-2]))
            return {"UTC": "fake", "JobState": "F", "QStat": "", "Status": None, "Reachable": True}

        def fake_read_json(host, path):
            events.append(("status", Path(path).parts[-2]))
            return job_status(Path(path).parts[-2])

        def fake_count(host, directory):
            events.append(("archive-count", Path(directory).parts[-2]))
            return 80 * 192

        def fake_fetch(host, remote_directory, local_directory):
            events.append(("fetch", Path(remote_directory).parts[-2]))
            shutil.copytree(ASSESSMENT / spec["Campaign"] / "results" / "main", local_directory, dirs_exist_ok=True)
            for name in ("worker-1", "worker-2", "reducer"):
                (Path(local_directory) / "jobs" / name).mkdir(parents=True, exist_ok=True)
                (Path(local_directory) / "jobs" / name / "status.json").write_text(json.dumps(job_status(name)))
            return ["rsync", "fake"]

        def fake_remote_sha256(host, paths):
            return {path: sha256(self.tmp / "split" / case_id / "results" / "main" / path.split("/main/", 1)[1]) for path in paths}

        def fake_delete(host, archives):
            events.append(("delete", tuple(Path(a).parts[-2] for a in archives)))
            return {"Archives": list(archives), "SizesBeforeDeletion": "0", "DeletedUTC": "fake", "Remaining": ""}

        fakes = {"submit": fake_submit, "poll": fake_poll, "fetch": fake_fetch, "remote_sha256": fake_remote_sha256,
                 "delete_archives": fake_delete, "qstat_history": lambda host, pbs_bin, job: "job_state = F",
                 "read_json": fake_read_json, "count_archive_potentials": fake_count}
        saved = {name: getattr(qualify_library.remote_side, name) for name in fakes}
        saved_upload, saved_prepare = qualify_library.upload_case, qualify_library.prepare_case

        def prepare_with_recorded_controls(case_record, **kwargs):
            kwargs["args"].control_source = spec["Controls"]
            kwargs["args"].stage_prefix = spec["Prefix"]
            record, context = saved_prepare(case_record, **kwargs)
            holder["record"], holder["context"] = record, context
            return record, context
        try:
            for name, fake in fakes.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = fake_upload
            qualify_library.prepare_case = prepare_with_recorded_controls
            manifest = json.loads(MANIFEST.read_text())
            case_library = qualify_library.source_directory(MANIFEST, manifest, qualify_library.manifest_case(manifest, case_id)) \
                / qualify_library.manifest_case(manifest, case_id)["Source"]["Files"]["ProcessLibrary"]["Name"]
            model_name = json.loads(case_library.read_text())["Models"][0]["Name"]
            previous_library = self.tmp / "previous-process-library.json"
            previous_library.write_text(json.dumps({"Version": 1, "Root": "/tmp/previous", "Models": [
                {"Name": "other-coupon", "LibraryQualified": False},
                {"Name": model_name, "LibraryQualified": False, "Stale": True}]}))
            args = argparse.Namespace(build_record=self.build_record, reference=ASSESSMENT / spec["Campaign"] / "reference", remote="h:/r",
                                      orders=[], controls=[3, 5], control_count=8, control_source=None, max_jobs=2,
                                      frozen_binary_sha256=BINARY_SHA256, stage_prefix=None, case=[case_id], root=self.tmp / "split",
                                      dry_run=False, resume=False, monitor_interval=0, monitor_polls=10,
                                      cluster_profile=HERE / "qualify" / "cluster-profile.json", cost_model=estimate_stages.COST_MODEL,
                                      gates=gates.GATES_FILE, job_policy="fixed", fixed_jobs=2, merge_into=previous_library,
                                      reducer_block_size=None)
            record = qualify_library.run_qualify(args, log=lambda message: None)
        finally:
            for name, fake in saved.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = saved_upload
            qualify_library.prepare_case = saved_prepare
        case = record["Cases"][0]
        self.assertEqual(case["Status"], "qualified", case.get("StoppedBy"))
        self.assertEqual((case["Split"]["N"], case["JobPolicy"]["Mode"], case["JobPolicy"]["FixedJobs"]), (2, "fixed", 2))
        blocks = holder["context"]["split"]["Blocks"]
        self.assertEqual(case["Split"]["Blocks"], [len(blocks[0]), len(blocks[1])])
        self.assertEqual(blocks[0] + blocks[1], list(range(1, 81)))
        kinds = [kind for kind, _ in events]
        self.assertEqual(events[:3], [("upload", case_id), ("submit", "worker-1"), ("submit", "worker-2")])
        self.assertLess(kinds.index("poll"), kinds.index("archive-count"))
        self.assertEqual([name for kind, name in events if kind == "status"], ["worker-1", "worker-2"])
        self.assertEqual(events.index(("archive-count", f"{spec['Prefix']}-p4")) + 1, events.index(("submit", "reducer")))
        self.assertLess(events.index(("submit", "reducer")), events.index(("fetch", case_id)))
        self.assertEqual([name for kind, name in events if kind == "submit"], ["worker-1", "worker-2", "reducer"])
        self.assertEqual(kinds.count("fetch"), 1)
        self.assertEqual(kinds.count("delete"), 1)
        self.assertEqual(case["ArchiveUnion"]["Stages"][f"{spec['Prefix']}-p4"]["OK"], True)
        self.assertEqual(case["ArchiveUnion"]["Stages"][f"{spec['Prefix']}-p4"]["Expected"], 80 * 192)
        self.assertEqual([job["Name"] for job in case["Jobs"]], ["worker-1", "worker-2", "reducer"])
        self.assertEqual(case["Jobs"][2]["Requires"], ["worker-1", "worker-2"])
        for job in case["Jobs"]:
            self.assertEqual(job["Submission"]["Job"], f"{job['Name']}.fake")
            self.assertTrue(Path(job["SubmissionRecord"]).is_file())
            self.assertEqual(job["Monitor"]["LastJobState"], "F")
        self.assertEqual(case["Jobs"][0]["Status"]["State"], "complete")
        self.assertEqual(case["Jobs"][0]["StageNames"][1:], [f"{spec['Prefix']}-p5-control-worker", f"{spec['Prefix']}-p5-control-reducer",
                                                              f"{spec['Prefix']}-p3-control-worker", f"{spec['Prefix']}-p3-control-reducer",
                                                              f"{spec['Prefix']}-p4-local-edge"])
        # Verdict and gate values from the same reducer matrices as the single-job record.
        self.assertEqual(case["Qualification"]["Verdict"], gates.VERDICT_PASSED)
        cost = case["Cost"]
        self.assertEqual(cost["Split"], {**cost["Split"], "N": 2, "Policy": "fixed", "ControlsJob": "worker-1"})
        self.assertEqual(set(cost["Jobs"]), {"worker-1", "worker-2", "reducer"})
        for name, job in cost["Jobs"].items():
            self.assertAlmostEqual(job["ActualSeconds"], job_status(name)["TotalSeconds"])
            self.assertGreater(job["EstimateSecondsWithPreflightAndMargin"], 0.0)
            self.assertEqual(job["PBSJobID"], f"{name}.fake")
        self.assertAlmostEqual(cost["JobTotalSeconds"], sum(job["ActualSeconds"] for job in cost["Jobs"].values()))
        self.assertAlmostEqual(cost["JobNodeHours"], cost["JobTotalSeconds"] / 3600.0)
        self.assertIsNotNone(cost["CriticalPathSeconds"])
        summary = json.loads((self.tmp / "split" / case_id / "cost-summary.json").read_text())
        main = summary["Stages"][f"{spec['Prefix']}-p4"]
        self.assertEqual(main["Sources"], list(range(1, 81)))
        self.assertEqual(len(main["PCGIterations"]), 80)
        recorded_main = {stage["Name"]: stage for stage in recorded["Stages"]}
        self.assertAlmostEqual(main["WorkerWallSeconds"], recorded_main[f"{spec['Prefix']}-p4-worker"]["WallSeconds"])
        self.assertAlmostEqual(main["ReducerWallSeconds"], recorded_main[f"{spec['Prefix']}-p4-reducer"]["WallSeconds"])
        self.assertTrue((self.tmp / "split" / case_id / "results" / "main" / "jobs" / "reducer" / "status.json").is_file())
        self.assertTrue((self.tmp / "split" / case_id / "archive-union.json").is_file())
        totals = record["Library"]
        self.assertEqual(totals["JobsSubmitted"], 3)
        self.assertEqual(totals["Splits"][case_id]["N"], 2)
        self.assertEqual(totals["JobPolicy"]["Mode"], "fixed")
        self.assertAlmostEqual(totals["NodeHours"], cost["JobNodeHours"])
        # --merge-into: the previous library's other model kept first, the same Name replaced by this run's, provenance recorded.
        merged = json.loads((self.tmp / "split" / "process-library.json").read_text())
        self.assertEqual([model["Name"] for model in merged["Models"]], ["other-coupon", model_name])
        self.assertNotIn("Stale", merged["Models"][1])
        self.assertEqual(merged["Models"][1]["CouponMesh"]["SHA256"], case["Mesh"]["SHA256"])
        self.assertEqual((merged["MergedFrom"]["Kept"], merged["MergedFrom"]["Replaced"]), (["other-coupon"], [model_name]))
        self.assertEqual(merged["MergedFrom"]["SHA256"], sha256(previous_library))
        self.assertIn("Stale", json.loads(previous_library.read_text())["Models"][1])   # the previous file is not modified
        # --resume of the split coupon: the three recorded submissions are adopted (the plans
        # are byte-identical; the identity check reads both sides in one order - the 2026-09-21
        # split acceptance resume stopped with 'plans differ' on identical plans: glob order
        # reducer / worker-1 / worker-2 against job order worker-1 / worker-2 / reducer).
        try:
            for name, fake in fakes.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = fake_upload
            qualify_library.prepare_case = prepare_with_recorded_controls
            args.resume = True
            del events[:]
            resumed = qualify_library.run_qualify(args, log=lambda message: None)
        finally:
            for name, fake in saved.items():
                setattr(qualify_library.remote_side, name, fake)
            qualify_library.upload_case = saved_upload
            qualify_library.prepare_case = saved_prepare
        resumed_case = resumed["Cases"][0]
        self.assertIsNone(resumed_case.get("StoppedBy"), resumed_case.get("StoppedBy"))
        self.assertEqual(resumed_case["Status"], "qualified")
        self.assertEqual([kind for kind, _ in events if kind in ("upload", "submit")], [])
        self.assertEqual([job["Monitor"]["Resumed"] for job in resumed_case["Jobs"]], [True, True, True])
        self.assertEqual([job["Submission"]["Job"] for job in resumed_case["Jobs"]], ["worker-1.fake", "worker-2.fake", "reducer.fake"])

    def test_without_reference_matrices_the_verdict_is_pending(self):
        case_id = "four-edge-9d2cb9bbb3fe"
        campaign = ASSESSMENT / CASES[case_id]["Campaign"]
        inputs_only = self.tmp / "inputs-only-campaign"
        if not inputs_only.exists():
            inputs_only.mkdir()
            os.symlink(campaign / "reference" / "inputs-07", inputs_only / "inputs-07")
        record, context, root, table, digest, profile, manifest = self.analysis_context(case_id, inputs_only)
        self.assertIsNone(record["Reference"]["Results"])
        before = self.snapshot(campaign)
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
        # Without --dry-run the remote is mandatory; --reference must be spelled ('none' included).
        result = subprocess.run([sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                                 "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, "--root", str(self.tmp / "no-remote")],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--remote", result.stderr)
        result = subprocess.run([sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                                 "--frozen-binary-sha256", BINARY_SHA256, "--root", str(self.tmp / "no-reference"), "--dry-run"],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--reference", result.stderr)
        # A reference whose config differs from the case's own sources fails closed.
        altered = self.tmp / "altered-reference"
        if not altered.exists():
            altered.mkdir()
            shutil.copytree(campaign / "reference" / "inputs-07", altered / "inputs-07")
            producer = json.loads((altered / "inputs-07" / "spatial_fabricated.json").read_text())
            producer["Domains"]["Materials"][0]["Permittivity"] = 9.0
            (altered / "inputs-07" / "spatial_fabricated.json").write_text(json.dumps(producer, indent=2))
        root = self.tmp / "dry-altered"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", str(altered), "--frozen-binary-sha256", BINARY_SHA256, "--case", case_id, "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        case = json.loads((root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(case["StoppedBy"]["Kind"], "Reference")
        self.assertIn("Permittivity", str(case["StoppedBy"]["Differences"]))

    def test_run_config_derived_from_the_case_equals_every_gallery_reference(self):
        """Decision 52: for all five gallery cases the config derived from the case's own
        sources (process library, trace basis, mesh $PhysicalNames, recipe PhysicsRun)
        equals the config the graded_v2 reference ran apart from Model.Mesh, Problem.Output
        and the DataFile directory; Solver.Order / Linear.Tol are the recipe's (the
        references at p5 / Tol 1e-8 - gallery 10 and the ten-edge - differ there only)."""
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        physics_run = case_inputs.physics_run_parameters(manifest)
        self.assertEqual((physics_run["Order"], physics_run["LinearTol"]), (4, 1e-10))
        checked = []
        for case_id, reference_path in GALLERY_REFERENCE_CONFIGS.items():
            mesh = local_identity_mesh(case_id)
            if mesh is None or not reference_path.is_file():
                continue
            case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
            directory = qualify_library.source_directory(MANIFEST, manifest, case)
            out = self.tmp / "derived" / case_id
            config, record = case_inputs.derive(case, directory, mesh_path=mesh, physics_run=physics_run, out_dir=out)
            reference = json.loads(reference_path.read_text())
            self.assertEqual(case_inputs.config_differences(config, reference, ignore_solver=("Order", "Linear.Tol")), [], case_id)
            self.assertEqual(config["Solver"]["Order"], 4)
            self.assertEqual(config["Solver"]["Linear"]["Tol"], 1e-10)
            self.assertEqual(len(config["Boundaries"]["PrescribedPotential"]), len(reference["Boundaries"]["PrescribedPotential"]))
            self.assertEqual(record["Interfaces"], case_inputs.interface_types(reference))
            self.assertTrue(record["AttributeCheck"]["Passed"])
            # Every regenerated trace has the reference trace's name; the digests differ from the
            # producer's files by the canonical-frame round trip only (FrameFitResidual, recorded).
            self.assertEqual([source["Name"] for source in record["Sources"]],
                             [Path(entry["DataFile"]).name for entry in reference["Boundaries"]["PrescribedPotential"]])
            self.assertIsNotNone(record["Traces"]["FrameFitResidual"])
            reference_traces = reference_path.parent / "traces" if (reference_path.parent / "traces").is_dir() else \
                reference_path.parents[1] / f"inputs-{reference_path.parent.name.split('-')[1]}" / "traces"
            if reference_traces.is_dir():
                for source in record["Sources"][:3] + record["Sources"][-1:]:
                    candidate = reference_traces / source["Name"]
                    if not candidate.is_file():
                        candidate = reference_traces.parent / source["Name"]
                    ours = [line.split(",") for line in Path(source["Path"]).read_text().splitlines()[1:]]
                    theirs = [line.split(",") for line in candidate.read_text().splitlines()[1:]]
                    self.assertEqual(len(ours), len(theirs))
                    for a, b in zip(ours, theirs):
                        self.assertEqual((a[3], a[4]), (b[3], b[4]))     # V and triangle columns identical
                        self.assertTrue(all(abs(float(x) - float(y)) <= 1e-13 for x, y in zip(a[:3], b[:3])), (source["Name"], a, b))
            checked.append(case_id)
        self.assertGreaterEqual(len(checked), 4, checked)
        if len(checked) < 5:
            self.skipTest(f"derived-config equality proven for {checked}; the others lack a local mesh or reference config")

    def test_reference_none_plans_the_coupon_on_its_own_inputs(self):
        case_id = "four-edge-9d2cb9bbb3fe"
        root = self.tmp / "dry-none"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, "--case", case_id, "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = json.loads((root / "library-qualification.json").read_text())
        case = record["Cases"][0]
        self.assertEqual(case["Status"], "planned", case.get("StoppedBy"))
        self.assertIn("--reference none", case["Reference"]["Rule"])
        self.assertEqual(case["Orders"], {**case["Orders"], "Main": ["p4"], "Gated": "p4", "ReferenceOrder": None, "RecipeOrder": "p4"})
        self.assertEqual(case["Sources"]["Count"], 80)
        self.assertEqual(case["Inputs"]["Origin"], "case")
        self.assertEqual(case["Configs"]["ReferenceConfig"], None)
        self.assertEqual(len(case["Inputs"]["Sources"]), 80)
        self.assertTrue((root / case_id / "inputs" / "traces" / "basis-0080.csv").is_file())
        worker = json.loads((root / case_id / "main" / f"{case_id}-p4" / "worker.json").read_text())
        self.assertEqual(worker["Solver"]["Linear"]["Tol"], 1e-10)
        self.assertEqual(worker["Solver"]["Order"], 4)
        self.assertEqual(worker["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeFrameNormal"], [0.0, 0.0, 1.0])
        # The analysis without a reference: the recorded results give PendingQualification.
        campaign = ASSESSMENT / CASES[case_id]["Campaign"]
        record_, context, root_, table, digest, profile, manifest = self.analysis_context(case_id, None)
        before = self.snapshot(campaign)
        gate_record = qualify_library.analyze_case(record_, context, campaign / "results", gates=table, gates_digest=digest,
                                                   profile=profile)
        self.assert_campaign_untouched(campaign, before)
        self.assertEqual(gate_record["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(record_["Status"], "pending-qualification")
        self.assertIsNone(record_["Cost"]["ReferenceNodeHours"])


if __name__ == "__main__":
    unittest.main()
