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
        context = {"jobs": [{"Name": "single", "Kind": "single", "RemoteDirectory": "/r/case/main", "Submission": {"Job": "1.h"},
                             "ExitCheck": {"OK": True, "Reasons": []}}]}
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


class RemoteFailClosedTest(unittest.TestCase):
    """Decision 485 (b), the b-batch1 O1 attempt 1 (PBS 57745 / 57749) reproduced: the shared remote root
    still held fab-/thin-spatial-2-edge-252cfb0e368d from an Oct-5 run; the job script died in 1 s at
    `mkdir "$D/tmp"` (pbs-status.json ExitCode 1), and the driver fetched the Oct-5 status.json (PBSJobID
    56380) as the job's results.  Now a pre-existing remote case directory stops the coupon before any
    upload unless --adopt-remote-case, and every job's exit record + status.json job id are checked before
    anything is fetched."""

    SUBMISSION = {"Job": "57745.ip-192-168-54-24.us-west-2.compute.internal", "UTC": "2026-10-07T03:11:13Z"}
    STALE_STATUS = {"Version": 2, "Case": "spatial-2-edge-252cfb0e368d", "PBSJobID": "56380.ip-192-168-54-24.us-west-2.compute.internal",
                    "StartUTC": "2026-10-05T09:25:05Z", "State": "complete", "TotalSeconds": 3816.0,
                    "Stages": [{"Name": "spatial-2-edge-252cfb0e368d-p4-worker", "State": "complete", "Parsed": {}}]}

    def fake_remote(self, files, reachable=True):
        remote_side = qualify_library.remote_side
        saved = {name: getattr(remote_side, name) for name in ("read_json_reachable", "fetch", "path_exists", "upload", "upload_file",
                                                               "qstat_history")}
        calls = []

        def read_json_reachable(host, path):
            calls.append(("read", path))
            return (True, files.get(path)) if reachable else (False, None)

        def fetch(host, remote_directory, local_directory):
            calls.append(("fetch", remote_directory))
            raise AssertionError("fetched after a failed exit check")
        remote_side.read_json_reachable, remote_side.fetch = read_json_reachable, fetch
        self.addCleanup(lambda: [setattr(remote_side, name, fake) for name, fake in saved.items()])
        return calls

    def test_job_exit_record_and_status_job_id_are_checked_before_the_fetch(self):
        remote_side = qualify_library.remote_side
        directory = "/r/fab-spatial-2-edge-252cfb0e368d/spatial-2-edge-252cfb0e368d/main"
        job_id = self.SUBMISSION["Job"]
        # The reproduced failure: the trap's exit record reads 1 and the status.json is another job's.
        calls = self.fake_remote({f"{directory}/pbs-status.json": {"ExitCode": 1, "JobID": job_id, "UTC": "2026-10-07T03:14:11Z"},
                                  f"{directory}/status.json": self.STALE_STATUS})
        check = remote_side.job_exit("h", directory, job_id)
        self.assertFalse(check["OK"])
        self.assertEqual(len(check["Reasons"]), 2)
        self.assertIn("ExitCode 1", check["Reasons"][0])
        self.assertIn("56380", check["Reasons"][1])
        self.assertEqual(check["Rule"], remote_side.JOB_EXIT_RULE)
        tmp = Path(tempfile.mkdtemp(prefix="qualify-exit-stop-"))
        self.addCleanup(shutil.rmtree, tmp, True)
        record = {"Case": "spatial-2-edge-252cfb0e368d", "Root": str(tmp), "Remote": {"Case": directory.rsplit("/", 1)[0]},
                  "Submission": dict(self.SUBMISSION)}
        job = {"Name": "single", "Kind": "single", "RemoteDirectory": directory, "Submission": dict(self.SUBMISSION)}
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.finish_case(record, {"jobs": [job]}, remote={"Host": "h"}, profile={"PBSBin": "/pbs"})
        self.assertEqual(stop.exception.record["Kind"], "JobExit")
        self.assertIn("nothing fetched", stop.exception.record["Message"])
        self.assertEqual(stop.exception.record["ExitCheck"]["StatusPBSJobID"], self.STALE_STATUS["PBSJobID"])
        self.assertFalse(job["ExitCheck"]["OK"])
        self.assertEqual([kind for kind, _ in calls], ["read", "read", "read", "read"])   # never "fetch"
        self.assertNotIn("Fetch", record)
        # A worker job of a split coupon is checked the same way before its status is read.
        worker = {"Name": "worker-1", "Kind": "worker", "RemoteDirectory": directory, "Submission": dict(self.SUBMISSION)}
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.complete_worker_job(record, {"jobs": [worker]}, worker, remote={"Host": "h"}, profile={})
        self.assertEqual(stop.exception.record["Kind"], "JobExit")
        # No exit record at all (the trap never ran) and a missing status.json: both named.
        self.fake_remote({})
        check = remote_side.job_exit("h", directory, job_id)
        self.assertFalse(check["OK"])
        self.assertTrue(check["Reachable"])
        self.assertEqual(len(check["Reasons"]), 2)
        self.assertIn("no pbs-status.json", check["Reasons"][0])
        # A lost connection is a TRANSPORT failure, not a dead job (decision 495 (3)): not OK, Reachable
        # False, the reason names --resume; the driver stops the coupon as Transport, not JobExit.
        calls = self.fake_remote({f"{directory}/pbs-status.json": {"ExitCode": 0, "JobID": job_id}}, reachable=False)
        check = remote_side.job_exit("h", directory, job_id)
        self.assertEqual((check["OK"], check["Reachable"]), (False, False))
        self.assertEqual(len(check["Reasons"]), 1)
        self.assertIn("transport failure", check["Reasons"][0])
        self.assertIn("--resume", check["Reasons"][0])
        self.assertNotIn("EXIT trap", check["Reasons"][0])
        self.assertEqual([kind for kind, _ in calls], ["read"])   # the status read is not attempted after the first round trip did not answer
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.check_job_exit(record, job, remote={"Host": "h"})
        self.assertEqual(stop.exception.record["Kind"], "Transport")
        self.assertFalse(job["ExitCheck"]["Reachable"])
        # A clean exit of the submitted job with its own status.json passes.
        self.fake_remote({f"{directory}/pbs-status.json": {"ExitCode": 0, "JobID": job_id, "UTC": "2026-10-07T04:19:00Z"},
                          f"{directory}/status.json": {**self.STALE_STATUS, "PBSJobID": job_id}})
        check = remote_side.job_exit("h", directory, job_id)
        self.assertTrue(check["OK"])
        self.assertEqual(check["ExitRecord"]["ExitCode"], 0)
        self.assertEqual(qualify_library.check_job_exit(record, job, remote={"Host": "h"})["OK"], True)
        self.assertTrue(job["ExitCheck"]["OK"])
        # The fetched status.json must also be the submitted job's (the post-fetch side of the check).
        results = tmp / "results" / "main"
        results.mkdir(parents=True)
        (results / "status.json").write_text(json.dumps(self.STALE_STATUS))
        remote_side.fetch = lambda host, remote_directory, local_directory: ["rsync", "fake"]
        remote_side.qstat_history = lambda host, pbs_bin, job_id: "Exit_status = 0"
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.finish_case(record, {"jobs": [job]}, remote={"Host": "h"}, profile={"PBSBin": "/pbs"})
        self.assertEqual(stop.exception.record["Kind"], "JobExit")
        self.assertIn("56380", stop.exception.record["Message"])

    def test_pre_existing_remote_case_directory_is_refused_unless_adopted(self):
        remote_side = qualify_library.remote_side
        tmp = Path(tempfile.mkdtemp(prefix="qualify-remote-case-"))
        self.addCleanup(shutil.rmtree, tmp, True)
        (tmp / "main").mkdir()
        (tmp / "main" / "plan.json").write_text("{}")
        mesh, trace = tmp / "identity.msh", tmp / "basis-0001.csv"
        mesh.write_text("mesh\n")
        trace.write_text("trace\n")
        remote_case = "/r/fab-spatial-2-edge-252cfb0e368d/spatial-2-edge-252cfb0e368d"
        record = {"Case": "spatial-2-edge-252cfb0e368d", "Root": str(tmp),
                  "Remote": {"Case": remote_case, "Run": remote_case.rsplit("/", 1)[0], "Mesh": f"{remote_case}/mesh/identity-0.msh"}}
        context = {"identity": {"Path": str(mesh)}, "sources": [{"Path": str(trace), "Name": "basis-0001.csv"}]}
        uploads = []
        saved = {name: getattr(remote_side, name) for name in ("ssh", "upload", "upload_file")}
        self.addCleanup(lambda: [setattr(remote_side, name, fake) for name, fake in saved.items()])
        remote_side.upload = lambda host, local, remote: uploads.append(("upload", remote)) or ["rsync", remote]
        remote_side.upload_file = lambda host, local, remote: uploads.append(("file", remote)) or ["rsync", remote]
        # The shared root still holds the Oct-5 case directory: refused before any upload.
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 0, stdout=remote_side.EXISTS_MARKER + "\n", stderr="")
        self.assertTrue(remote_side.path_exists("h", remote_case))
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.upload_case(record, context, remote={"Host": "h"}, profile={})
        self.assertEqual(stop.exception.record["Kind"], "RemoteCase")
        self.assertIn("--adopt-remote-case", stop.exception.record["Message"])
        self.assertEqual(stop.exception.record["Rule"], qualify_library.REMOTE_CASE_RULE)
        self.assertEqual(uploads, [])
        self.assertFalse((tmp / "upload").exists())
        # Adopted explicitly: uploaded, the adoption recorded.
        upload = qualify_library.upload_case(record, context, remote={"Host": "h"}, profile={}, adopt_remote_case=True)
        self.assertEqual((upload["RemoteCaseExisted"], upload["AdoptedRemoteCase"]), (True, True))
        self.assertEqual([kind for kind, _ in uploads], ["upload", "file"])
        # A fresh case directory uploads as before, nothing adopted.
        uploads.clear()
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 0, stdout=remote_side.ABSENT_MARKER + "\n", stderr="")
        upload = qualify_library.upload_case(record, context, remote={"Host": "h"}, profile={})
        self.assertEqual((upload["RemoteCaseExisted"], upload["AdoptedRemoteCase"]), (False, False))
        self.assertEqual(len(uploads), 2)
        # A lost connection is no answer: the coupon stops (never "absent" from silence).
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 255, stdout="", stderr="timed out")
        with self.assertRaisesRegex(RuntimeError, "no answer"):
            remote_side.path_exists("h", remote_case)
        # The marker-based JSON reader (decision 495 (3)): silence -> (False, None); an absent file -> (True, None); content -> parsed.
        self.assertEqual(remote_side.read_json_reachable("h", "/r/x.json"), (False, None))
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 0, stdout=remote_side.READ_MARKER + "\n", stderr="")
        self.assertEqual(remote_side.read_json_reachable("h", "/r/x.json"), (True, None))
        self.assertIsNone(remote_side.read_json("h", "/r/x.json"))
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 0, stdout='{"ExitCode": 0}\n' + remote_side.READ_MARKER + "\n", stderr="")
        self.assertEqual(remote_side.read_json_reachable("h", "/r/x.json"), (True, {"ExitCode": 0}))
        self.assertEqual(remote_side.read_json("h", "/r/x.json"), {"ExitCode": 0})
        remote_side.ssh = lambda host, command, check=True: subprocess.CompletedProcess(["ssh"], 255, stdout="", stderr="timed out")
        with self.assertRaises(qualify_library.CaseStop) as stop:
            qualify_library.upload_case(record, context, remote={"Host": "h"}, profile={})
        self.assertEqual(stop.exception.record["Kind"], "RemoteCase")
        # The dry run never asks (parse: the flag exists and defaults off).
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        self.assertFalse(parser.parse_args(["--build-record", "b", "--reference", "none", "--dry-run"]).adopt_remote_case)
        self.assertTrue(parser.parse_args(["--build-record", "b", "--reference", "none", "--adopt-remote-case"]).adopt_remote_case)


class ReusedMainTest(unittest.TestCase):
    """--controls-only --reuse-main ROOT (decisions 474 (A) / 479 / 485 (c)): the stored run's main stage is
    reused only when this run's inputs are identical to it - the identity mesh digest, the run config apart
    from the remote paths, every regenerated trace against the stored plan's pins, and the stored reducer
    CSVs at the digests the stored run verified against its remote; every mismatch fails closed, the stored
    root is never written, and the splice copies the CSVs byte-identically (re-hashed)."""

    CASE = "spatial-8-edge-05c322f6cda8"
    PREFIX = "spatial-8-edge-05c322f6cda8-p4"

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="qualify-reuse-main-"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.root = self.tmp / f"fab-{self.CASE}"
        case_root = self.root / self.CASE
        self.config = {"Problem": {"Type": "Electrostatic", "Output": "/old/root/main"}, "Model": {"Mesh": "/old/root/mesh/identity-cdb48405d063.msh"},
                       "Solver": {"Order": 4, "Linear": {"Tol": 1e-10}},
                       "Boundaries": {"PrescribedPotential": [{"Index": k, "DataFile": f"/old/root/inputs/traces/basis-{k:04d}.csv"} for k in (1, 2, 3)]}}
        (case_root / "inputs").mkdir(parents=True)
        (case_root / "inputs" / "run-config.json").write_text(json.dumps(self.config, indent=2) + "\n")
        self.traces = {f"basis-{k:04d}.csv": f"{k:064x}" for k in (1, 2, 3)}
        (case_root / "main" / "jobs" / "worker-1").mkdir(parents=True)
        (case_root / "main" / "jobs" / "worker-1" / "plan.json").write_text(json.dumps(
            {"PinnedSHA256": {**{f"/old/root/inputs/traces/{name}": digest for name, digest in self.traces.items()},
                              "/old/root/mesh/identity-cdb48405d063.msh": "c" * 64}}))
        reducer = case_root / "results" / "main" / self.PREFIX / "reducer"
        reducer.mkdir(parents=True)
        (reducer / "domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
        (reducer / "surface-response-matrix.csv").write_text("interface,edge,R (m),basis_i,basis_j,Q_ij (J)\n2,0,2e-6,1,1,0.5\n")
        digests = {name: sha256(reducer / name) for name in qualify_library.REUSED_MAIN_CSVS}
        self.stored = {"Case": self.CASE, "Status": "failed", "StoppedBy": None,
                       "Mesh": {"Local": "/old/identity.msh", "SHA256": "c" * 64, "Verified": True},
                       "Inputs": {"ConfigSHA256": sha256(case_root / "inputs" / "run-config.json")},
                       "Stages": [{"Prefix": self.PREFIX, "Role": "main", "Order": 4}, {"Prefix": f"{self.CASE}-p5-control", "Role": "control", "Order": 5}],
                       "Controls": {"Indices": [1, 2, 3, 8, 10, 31, 32, 108]},
                       "Jobs": [{"Name": "worker-1", "Kind": "worker", "Nodes": 1, "Ranks": 192, "Submission": {"Job": "57756.h"}},
                                {"Name": "reducer", "Kind": "reducer", "Nodes": 1, "Ranks": 192, "Submission": {"Job": "57843.h"}}],
                       "Nodes": {"MultiNode": False, "Main": 1, "Fixed": 1},
                       "ResultDigests": {f"main/{self.PREFIX}/reducer/{name}": {"Local": digest, "Remote": digest, "OK": True}
                                         for name, digest in digests.items()},
                       "Qualification": {"Verdict": "Failed", "UnjudgedTypes": []},
                       "Cost": {"MainStage": {"NodeHours": 1.2}, "JobNodeHours": 2.3}}
        (self.root / "library-qualification.json").write_text(json.dumps({"Version": 1, "ToolCommit": "7b95205da0", "Gates": {"SHA256": "g" * 64},
                                                                          "Cases": [self.stored]}))
        self.derived = json.loads(json.dumps(self.config))
        self.derived["Model"]["Mesh"] = "/new/root/mesh/identity-cdb48405d063.msh"
        self.derived["Problem"]["Output"] = "/new/root/main"
        for entry in self.derived["Boundaries"]["PrescribedPotential"]:
            entry["DataFile"] = entry["DataFile"].replace("/old/root", "/new/root")

    def load(self, **overrides):
        kwargs = dict(identity_sha256="c" * 64, run_config=self.derived, trace_digests=dict(self.traces), main_prefix=self.PREFIX)
        kwargs.update(overrides)
        return qualify_library.load_reused_main(self.root, self.CASE, **kwargs)

    def test_identical_inputs_reuse_the_stored_main_stage(self):
        before = {path: path.stat().st_mtime_ns for path in self.root.rglob("*") if path.is_file()}
        reused = self.load()
        self.assertEqual((reused["Case"], reused["MainPrefix"], reused["Root"]), (self.CASE, self.PREFIX, str(self.root)))
        self.assertTrue(reused["Mesh"]["Identical"] and reused["Traces"]["Identical"] and reused["Reducer"]["VerifiedAgainstStoredRecord"])
        self.assertEqual(reused["Traces"], {"Count": 3, "PinnedByStoredPlans": 3, "Identical": True})
        self.assertEqual(reused["RunConfig"]["IdenticalApartFrom"], list(case_inputs.PATH_FIELDS))
        self.assertEqual((reused["Reducer"]["Job"], reused["Reducer"]["PBSJobID"], reused["Reducer"]["Nodes"], reused["Reducer"]["Ranks"]),
                         ("reducer", "57843.h", 1, 192))
        self.assertEqual(sorted(reused["Reducer"]["SHA256"]), sorted(qualify_library.REUSED_MAIN_CSVS))
        self.assertEqual((reused["StoredVerdict"], reused["StoredControls"]), ("Failed", [1, 2, 3, 8, 10, 31, 32, 108]))
        self.assertEqual((reused["StoredMainStageCost"], reused["StoredJobNodeHours"]), ({"NodeHours": 1.2}, 2.3))
        self.assertEqual(reused["Record"]["ToolCommit"], "7b95205da0")
        self.assertEqual(reused["Rule"], job_split.CONTROLS_ONLY_RULE)
        # The splice: byte-identical copies in this run's results, re-hashed; the stored root untouched.
        results = self.tmp / "new" / "results"
        splice = qualify_library.splice_reused_main(results, reused)
        self.assertEqual(splice["SHA256"], reused["Reducer"]["SHA256"])
        for name in qualify_library.REUSED_MAIN_CSVS:
            self.assertEqual((results / "main" / self.PREFIX / "reducer" / name).read_bytes(),
                             (Path(reused["Reducer"]["Directory"]) / name).read_bytes())
        self.assertEqual({path: path.stat().st_mtime_ns for path in self.root.rglob("*") if path.is_file()}, before)
        self.assertEqual(qualify_library.reducer_ranks_of({"ReusedMain": reused, "Jobs": []}), 192)

    def test_every_identity_mismatch_fails_closed(self):
        def stop(message, **overrides):
            with self.assertRaises(qualify_library.CaseStop) as context:
                self.load(**overrides)
            self.assertEqual(context.exception.record["Kind"], "ReuseMain")
            self.assertIn(message, context.exception.record["Message"])
        stop("not the same coupon mesh", identity_sha256="d" * 64)
        other_order = json.loads(json.dumps(self.derived))
        other_order["Solver"]["Order"] = 5
        stop("differs from the stored run's", run_config=other_order)
        other_trace = dict(self.traces)
        other_trace["basis-0002.csv"] = "e" * 64
        stop("1 of 3 regenerated traces differ", trace_digests=other_trace)
        stop("are not this run's", main_prefix="other-p4")
        with self.assertRaises(qualify_library.CaseStop) as context:
            qualify_library.load_reused_main(self.root, "spatial-8-edge-818f8956d075", identity_sha256="c" * 64, run_config=self.derived,
                                             trace_digests=self.traces, main_prefix=self.PREFIX)
        self.assertIn("holds no case", context.exception.record["Message"])
        with self.assertRaises(qualify_library.CaseStop) as context:
            qualify_library.load_reused_main(self.tmp, self.CASE, identity_sha256="c" * 64, run_config=self.derived,
                                             trace_digests=self.traces, main_prefix=self.PREFIX)
        self.assertIn("holds no library-qualification.json", context.exception.record["Message"])
        # A rewritten stored reducer CSV (the decision-475 process rule's failure mode) is not reusable ...
        reducer = self.root / self.CASE / "results" / "main" / self.PREFIX / "reducer"
        original = (reducer / "domain-response-matrix.csv").read_text()
        (reducer / "domain-response-matrix.csv").write_text(original + "1,2,0.1\n")
        stop("is not the one the stored run fetched and verified")
        (reducer / "domain-response-matrix.csv").write_text(original)
        self.load()
        # ... as is a rewritten stored run config, a stored run that never reached its analysis, or one
        # whose plans pin no traces.
        config_path = self.root / self.CASE / "inputs" / "run-config.json"
        config_path.write_text(json.dumps(self.config))
        stop("the stored root was rewritten")
        config_path.write_text(json.dumps(self.config, indent=2) + "\n")
        record_path = self.root / "library-qualification.json"
        record = json.loads(record_path.read_text())
        record["Cases"][0]["Qualification"] = None
        record_path.write_text(json.dumps(record))
        stop("was not analyzed")
        record["Cases"][0]["Qualification"] = {"Verdict": "Failed"}
        record_path.write_text(json.dumps(record))
        self.load()
        (self.root / self.CASE / "main" / "jobs" / "worker-1" / "plan.json").write_text(json.dumps({"PinnedSHA256": {}}))
        stop("0 pinned")

    def test_prepare_case_never_writes_the_shared_args(self):
        """Decision 495 (1), MAJOR-1 of the review: prepare_case once wrote args.control_amplitudes (the one Namespace
        shared by every coupon of the run), so in a multi-coupon controls-only run every later coupon took the
        FIRST coupon's reducer as its prior main stage.  The prior main stage is a local value now; no attribute
        of `args` is assigned anywhere in qualify_library (the AST is the guard; the two-case dry run
        test_controls_only_two_cases_take_their_own_reducers is the behavioural one)."""
        import ast
        tree = ast.parse(Path(qualify_library.__file__).read_text())
        writes = []
        for node in ast.walk(tree):
            targets = []
            if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                for leaf in ast.walk(target):
                    if isinstance(leaf, ast.Attribute) and isinstance(leaf.value, ast.Name) and leaf.value.id == "args":
                        writes.append((leaf.lineno, ast.unparse(leaf)))
        self.assertEqual(writes, [])
        source = Path(qualify_library.__file__).read_text()
        self.assertNotIn("args.control_amplitudes =", source)
        self.assertIn("control_amplitudes = args.control_amplitudes", source)

    def test_controls_only_and_reuse_main_go_together(self):
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        args = parser.parse_args(["--build-record", "b", "--reference", "none", "--controls-only", "--reuse-main", str(self.root)])
        self.assertTrue(args.controls_only)
        self.assertEqual(args.reuse_main, self.root)
        build = self.tmp / "library-build.json"
        build.write_text(json.dumps({"Cases": [], "Library": {"Commit": "x", "Manifest": {"Path": str(MANIFEST), "SHA256": sha256(MANIFEST)}}}))
        for extra in (["--controls-only"], ["--reuse-main", str(self.root)]):
            args = parser.parse_args(["--build-record", str(build), "--reference", "none", "--dry-run", "--root", str(self.tmp / "run"), *extra])
            with self.assertRaisesRegex(ValueError, "go together"):
                qualify_library.run_qualify(args, log=lambda message: None)
        with self.assertRaisesRegex(ValueError, "go together"):
            qualify_library.prepare_case({"Case": "c"}, manifest_path=MANIFEST, manifest={}, args=argparse.Namespace(controls_only=True, reuse_main=None),
                                         root=self.tmp, remote=None, profile={}, cost_model={}, gates={}, gates_digest="")


class JobSplitTest(unittest.TestCase):
    """The per-coupon source split policy (decision 61b) on the recorded device coupon
    spatial-3-edge-5d3b5e644745 (225 sources, H1 42.7M at p4: 29,272 s single-job estimate
    against the 21,600 s walltime - the coupon the first device library run failed closed)."""

    @classmethod
    def setUpClass(cls):
        build = json.loads((HERE / "qualify" / "device-library-20260921" / "library-build.json").read_text())
        cls.counts = next(case for case in build["Cases"] if case["Case"] == "spatial-3-edge-5d3b5e644745")["H1"]["EntityCounts"]
        # The recorded decision-61b split ran on the physics-11 model (PREVIOUS_COST_MODEL since
        # the decision-64a refit): reproduced with it; the refit model's own 7f03 outcome below.
        cls.model = estimate_stages.load_cost_model(estimate_stages.PREVIOUS_COST_MODEL)
        # The decision-64a device refit (the stage-2 plans of record); the decision-457 (2)
        # measured refit is the default since (test_refit_cost_model.MeasuredRefitTest).
        cls.model_refit = estimate_stages.load_cost_model(estimate_stages.DEVICE_COST_MODEL)
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

    def test_refit_model_keeps_7f03_fail_closed_as_one_job_and_reproduces_the_recorded_speed_split(self):
        # The decision-64a refit model (qualify/cost-model-device-20260922.json, from the 2026-09-22 run): 7f03
        # still does not fit one 6 h job at 2.0x PCG with preflight and margin, and the speed
        # policy at --max-jobs 6 gives the split the run recorded (PBS 47214-47300: six
        # worker jobs, the controls + local-edge alone in job 1, blocks of 45 sources).
        stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in self.layout if item["Kind"] == "response"]
        local = next(item for item in self.layout if item["Kind"] == "local-edge")
        estimate = estimate_stages.estimate(self.counts, stages, model=self.model_refit, profile=self.profile,
                                            local_edge=(local["EstimateKey"], local["Order"], local["Sources"]))
        self.assertEqual(estimate["ReducerBlockSize"], 48)
        self.assertFalse(estimate["FitsOneJob"])
        self.assertGreater(estimate["JobSecondsEstimateWithPreflightAndMargin"]["2.0"], self.profile["WalltimeSeconds"])
        # The refit's reducer rate (streaming Gram at b = 48) is far below the b28 / b = 6 one.
        self.assertLess(estimate["Stages"]["p4-225"]["ReducerSecondsEstimate"],
                        0.5 * self.estimate_default["Stages"]["p4-225"]["ReducerSecondsEstimate"])
        policy = job_split.normalize_policy("speed", max_jobs=6, walltime_seconds=self.profile["WalltimeSeconds"],
                                            user_job_cap=self.profile["UserJobCap"])
        record = job_split.plan_split(indices=self.indices, layout=self.layout, estimate=estimate, policy=policy,
                                      model=self.model_refit, profile=self.profile)
        self.assertEqual((record["N"], record["ControlsJob"]), (6, "separate"))
        self.assertEqual([len(block) for block in record["Blocks"]], [0, 45, 45, 45, 45, 45])
        self.assertEqual([item["Fits"] for item in record["Candidates"]], [False, True, True, True, True, True])

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

    def test_controls_only_plans_one_job_of_the_fixed_stages_alone(self):
        """Decisions 474 (A) / 479 / 485 (c): the controls-only re-qualification is ONE job of the control +
        local-edge stages (no source block, no reducer) whose estimate is the separate job 1's of a split, on
        the Fixed node count; its plan (build_job_plan) pins the mesh, every trace and the fixed stages'
        configs only - the b-batch1 lane script's controls-only plan, now the driver's."""
        policy = job_split.normalize_policy("frugal", max_jobs=4, walltime_seconds=self.profile["WalltimeSeconds"], user_job_cap=40)
        record = job_split.plan_split(indices=self.indices, layout=self.layout, estimate=self.estimate_default, policy=policy,
                                      model=self.model_refit, profile=self.profile, controls_only=True)
        self.assertTrue(record["ControlsOnly"])
        self.assertEqual(record["Rule"], job_split.CONTROLS_ONLY_RULE)
        self.assertEqual((record["N"], record["Blocks"], record["ControlsJob"]), (1, [[]], "controls-only"))
        self.assertEqual([job["Kind"] for job in record["Jobs"]], ["controls-only"])
        job = record["Jobs"][0]
        self.assertEqual((job["Name"], job["Sources"], job["Block"]), ("controls-only", [], 1))
        self.assertTrue(job["Fits"] and record["Fits"])
        self.assertIn("controls-only: ONE job", record["Decision"])
        self.assertIn("the main stages are reused", record["Decision"])
        # The same seconds as a separate job 1 of a split that leaves it no block (8 sources in 4 jobs).
        layout = qualify_library.stage_layout("c", [4], [3, 5], 8, 8)
        stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in layout if item["Kind"] == "response"]
        local = next(item for item in layout if item["Kind"] == "local-edge")
        estimate = estimate_stages.estimate(self.counts, stages, model=self.model_refit, profile=self.profile,
                                            local_edge=(local["EstimateKey"], local["Order"], local["Sources"]))
        fixed = job_split.normalize_policy("fixed", max_jobs=4, walltime_seconds=self.profile["WalltimeSeconds"], fixed_jobs=4)
        split = job_split.plan_split(indices=list(range(1, 9)), layout=layout, estimate=estimate, policy=fixed, model=self.model_refit,
                                     profile=self.profile)
        controls = job_split.plan_split(indices=list(range(1, 9)), layout=layout, estimate=estimate, policy=policy, model=self.model_refit,
                                        profile=self.profile, controls_only=True)
        self.assertEqual(split["ControlsJob"], "separate")
        self.assertEqual(controls["Jobs"][0]["SecondsEstimateWithPreflightAndMargin"], split["Jobs"][0]["SecondsEstimateWithPreflightAndMargin"])
        self.assertNotIn("Nodes", controls["Jobs"][0])
        multi = job_split.plan_split(indices=list(range(1, 9)), layout=layout, estimate=estimate, policy=policy, model=self.model_refit,
                                     profile=self.profile, controls_only=True, nodes={"Main": 4, "Fixed": 2})
        self.assertEqual(multi["Jobs"][0]["Nodes"], 2)
        self.assertEqual(multi["NodeSecondsEstimate"]["2.0"], 2 * multi["Jobs"][0]["SecondsEstimateWithPreflightAndMargin"]["2.0"])
        # The plan of the controls-only job: the fixed stages only, their configs + the mesh + every trace pinned.
        digests = {"c-p4": {"worker-block1.json": "a" * 64, "reducer.json": "b" * 64},
                   "c-p5-control": {"worker.json": "c" * 64, "reducer.json": "d" * 64},
                   "c-p3-control": {"worker.json": "e" * 64, "reducer.json": "f" * 64},
                   "c-p4-local-edge": {"config.json": "0" * 64}}
        pins = {f"/r/case/inputs/traces/basis-{k:04d}.csv": f"{k:064x}" for k in range(1, 9)}
        plan = build_plan.build_job_plan(case_id="c", job_name="controls-only", remote_case_root="/r/case",
                                         mesh={"Remote": "/r/case/mesh/m.msh", "SHA256": "9" * 64, "Local": "/l/m.msh"}, stage_layout=layout,
                                         split_job={**controls["Jobs"][0], "Kind": "worker"}, estimate=estimate, config_digests=digests,
                                         trace_pins=pins, profile=self.profile, binary="/r/p.bin", binary_sha256="1" * 64, mpiexec="/r/mpi",
                                         purpose="t", factors=["1.0", "1.5", "2.0"])
        self.assertEqual(plan["StageNames"], ["c-p5-control-worker", "c-p5-control-reducer", "c-p3-control-worker", "c-p3-control-reducer",
                                              "c-p4-local-edge"])
        self.assertEqual(plan["BlockSources"], [])
        self.assertEqual(sorted(plan["PinnedSHA256"]), sorted(list(pins) + ["/r/case/mesh/m.msh", "/r/case/main/c-p5-control/worker.json",
                                                                              "/r/case/main/c-p5-control/reducer.json", "/r/case/main/c-p3-control/worker.json",
                                                                              "/r/case/main/c-p3-control/reducer.json", "/r/case/main/c-p4-local-edge/config.json"]))
        self.assertNotIn("/r/case/main/c-p4/reducer.json", plan["PinnedSHA256"])

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
                                              job_name="j", walltime_seconds=21600, job_directory="/r/case/main/jobs/worker-2",
                                              instance_type=plans["worker-2"]["Instance"]["Type"])
        self.assertIn("D=/r/case/main/jobs/worker-2\n", script)
        self.assertIn("#PBS -o /r/case/main/jobs/worker-2/pbs.log", script)
        self.assertIn("#PBS -l instance_type=m8g.48xlarge\n", script)
        single = build_plan.render_job_script(profile=self.profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                              job_name="j", walltime_seconds=21600, instance_type="r8g.48xlarge")
        self.assertIn("D=/r/case/main\n", single)
        self.assertIn("#PBS -l instance_type=r8g.48xlarge\n", single)
        with self.assertRaisesRegex(ValueError, "c8g.48xlarge.*not an instance"):
            build_plan.render_job_script(profile=self.profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                         job_name="j", walltime_seconds=21600, instance_type="c8g.48xlarge")

    def test_instance_per_job_by_estimated_peak(self):
        """USER decision 2026-09-22: m8g.48xlarge (768 GiB) is the default instance, r8g.48xlarge
        (1,536 GiB) the fallback when a job's estimated Palace peak exceeds MemoryFitFraction
        (0.6) x the m8g memory; the plan records the instance and its admission guard
        MinimumMemAvailableBytes = 0.6 x the chosen memory; every job of the device library
        (largest estimate 338 GB) runs on m8g."""
        self.assertEqual([item["Type"] for item in self.profile["Instances"]], ["m8g.48xlarge", "r8g.48xlarge"])
        self.assertNotIn("c8g.48xlarge", [item["Type"] for item in self.profile["Instances"]])
        self.assertIn("c8g.48xlarge (384 GiB) is deliberately absent", self.profile["InstanceRule"])
        gb_per_gib = self.model["PalaceGBPerGiB"]
        small = build_plan.select_instance(self.profile, 338.0, gb_per_gib)
        self.assertEqual((small["Type"], small["MemoryGiB"], small["Fits"]), ("m8g.48xlarge", 768, True))
        self.assertEqual(small["MinimumMemAvailableBytes"], int(0.6 * 768 * 1024 ** 3))
        self.assertAlmostEqual(small["EstimatedPalacePeakGiB"], 338.0 / gb_per_gib)
        edge = build_plan.select_instance(self.profile, 0.6 * 768 * gb_per_gib, gb_per_gib)
        self.assertEqual(edge["Type"], "m8g.48xlarge")
        large = build_plan.select_instance(self.profile, 0.6 * 768 * gb_per_gib * 1.001, gb_per_gib)
        self.assertEqual((large["Type"], large["Fits"]), ("r8g.48xlarge", True))
        self.assertEqual(large["MinimumMemAvailableBytes"], int(0.6 * 1536 * 1024 ** 3))
        huge = build_plan.select_instance(self.profile, 2000.0 * gb_per_gib, gb_per_gib)
        self.assertEqual((huge["Type"], huge["Fits"]), ("r8g.48xlarge", False))
        self.assertEqual(build_plan.largest_node_gib(self.profile), 1485.13)
        # Every job plan of the split carries its own instance from the stages it runs: the
        # reducer job's peak is the reducer estimate, a block worker's the worker estimate,
        # job 1 also the controls and the local-edge stage.
        record = self.split("speed", 4)
        digests = {"c-p4": {f"worker-block{k}.json": f"{k}" * 64 for k in range(1, 5)} | {"reducer.json": "r" * 64},
                   "c-p5-control": {"worker.json": "a" * 64, "reducer.json": "b" * 64},
                   "c-p3-control": {"worker.json": "c" * 64, "reducer.json": "d" * 64},
                   "c-p4-local-edge": {"config.json": "e" * 64}}
        common = dict(case_id="c", remote_case_root="/r/case", mesh={"Remote": "/r/case/mesh/m.msh", "SHA256": "m" * 64, "Local": "/l/m.msh"},
                      stage_layout=self.layout, estimate=self.estimate, config_digests=digests, trace_pins={},
                      profile=self.profile, binary="/r/b.bin", binary_sha256="0" * 64, mpiexec="/r/mpiexec_bound.sh",
                      purpose="test", factors=["1.0", "1.5", "2.0"])
        plans = {job["Name"]: build_plan.build_job_plan(job_name=job["Name"], split_job=job, **common) for job in record["Jobs"]}
        stages = self.estimate["Stages"]
        main, p5, p3, local = (next(item["EstimateKey"] for item in self.layout if item["Prefix"] == prefix)
                               for prefix in ("c-p4", "c-p5-control", "c-p3-control", "c-p4-local-edge"))
        self.assertEqual(plans["reducer"]["Instance"]["EstimatedPalacePeakGB"], stages[main]["ReducerPalacePeakGBEstimate"])
        self.assertEqual(plans["worker-2"]["Instance"]["EstimatedPalacePeakGB"], stages[main]["WorkerPalacePeakGBEstimate"])
        self.assertEqual(plans["worker-1"]["Instance"]["EstimatedPalacePeakGB"],
                         max(stages[main]["WorkerPalacePeakGBEstimate"], stages[p5]["WorkerPalacePeakGBEstimate"],
                             stages[p5]["ReducerPalacePeakGBEstimate"], stages[p3]["WorkerPalacePeakGBEstimate"],
                             stages[p3]["ReducerPalacePeakGBEstimate"], stages[local]["PalacePeakGBEstimate"]))
        for plan in plans.values():
            self.assertEqual(plan["Instance"]["Type"], "m8g.48xlarge")
            self.assertEqual(plan["MinimumMemAvailableBytes"], plan["Instance"]["MinimumMemAvailableBytes"])
        whole = build_plan.build_plan(**common)
        self.assertEqual(whole["Instance"]["EstimatedPalacePeakGB"], self.estimate["MaxPalacePeakGBEstimate"])
        self.assertEqual(whole["Instance"]["Type"], "m8g.48xlarge")

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

    def test_manifest_frozen_executable_default_and_command_line_override(self):
        """Decision 63: the production recipe records PhysicsRun.FrozenExecutable 170439c4...
        (PreviousSHA256 b28f089a..., the streaming one-pass Gram executable of decision 62(4));
        case_inputs and general_mesh_manifest validate it; --frozen-binary-sha256 overrides
        it and the origin is recorded; the ReducerBlockSize rule records that b = N is
        permitted with the streaming executable."""
        from general_mesh_manifest import validate_physics_run
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        block = manifest["ProductionRecipe"]["PhysicsRun"]["FrozenExecutable"]
        self.assertEqual(block["SHA256"], "9ef5256bc9ca9abc109954fb97b096879b3b472f18483f867b8de0b73a832058")
        self.assertEqual(block["PreviousSHA256"], "170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b")
        self.assertEqual(block["CostModelSHA256"], build_plan.COST_MODEL_FROZEN_BINARY_SHA256)
        self.assertEqual(block["SHA256"], build_plan.DEFAULT_FROZEN_BINARY_SHA256)
        self.assertEqual(block["PreviousSHA256"], build_plan.PREVIOUS_FROZEN_BINARY_SHA256)
        self.assertIn("decision 69", block["Provenance"])
        self.assertIn("62(4)", block["Provenance"])
        self.assertIn("b = N", manifest["ProductionRecipe"]["PhysicsRun"]["ReducerBlockSize"]["Rule"])
        self.assertIn("b = N", build_plan.REDUCER_BLOCK_SIZE_RULE)
        physics_run = case_inputs.physics_run_parameters(manifest)
        self.assertEqual(physics_run["FrozenExecutableSHA256"], block["SHA256"])
        validate_physics_run(manifest["ProductionRecipe"])
        for broken in ({"SHA256": "abc", "Rule": "x", "Provenance": "y"}, {"SHA256": block["SHA256"]},
                       {"SHA256": block["SHA256"], "PreviousSHA256": "b28", "Rule": "x", "Provenance": "y"}, block["SHA256"]):
            recipe = json.loads(json.dumps(manifest["ProductionRecipe"]))
            recipe["PhysicsRun"]["FrozenExecutable"] = broken
            with self.assertRaisesRegex(ValueError, "FrozenExecutable"):
                validate_physics_run(recipe)
            if not (isinstance(broken, dict) and broken.get("Rule") and broken.get("PreviousSHA256") == "b28"):
                with self.assertRaisesRegex(case_inputs.CaseInputError, "FrozenExecutable"):
                    case_inputs.physics_run_parameters({**manifest, "ProductionRecipe": recipe})
        chosen = qualify_library.frozen_binary_of(argparse.Namespace(frozen_binary_sha256=None), physics_run)
        self.assertEqual((chosen["SHA256"], chosen["Origin"]), (block["SHA256"], "manifest ProductionRecipe.PhysicsRun.FrozenExecutable"))
        chosen = qualify_library.frozen_binary_of(argparse.Namespace(frozen_binary_sha256=BINARY_SHA256), physics_run)
        self.assertEqual((chosen["SHA256"], chosen["Origin"]), (BINARY_SHA256, "--frozen-binary-sha256"))
        chosen = qualify_library.frozen_binary_of(argparse.Namespace(frozen_binary_sha256=None), {"Order": 4})
        self.assertEqual((chosen["SHA256"], chosen["Origin"]),
                         (build_plan.DEFAULT_FROZEN_BINARY_SHA256, "built-in default build_plan.DEFAULT_FROZEN_BINARY_SHA256"))
        self.assertIn("8.4e-13", chosen["Rule"])
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        self.assertIsNone(parser.parse_args(["--build-record", "b", "--reference", "none", "--remote", "h:/r"]).frozen_binary_sha256)
        self.assertEqual(parser.parse_args(["--build-record", "b", "--reference", "none", "--remote", "h:/r",
                                            "--frozen-binary-sha256", "f"]).frozen_binary_sha256, "f")

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
        # The recorded campaign ran on r8g.48xlarge with a 1,000 GiB guard; since the instance
        # decision the plan chooses its instance by the estimated peak (m8g here) and its guard
        # is 0.6 x that instance's memory.
        self.assertEqual(recorded_plan["MinimumMemAvailableBytes"], 1073741824000)
        self.assertEqual(plan["Instance"]["Type"], "m8g.48xlarge")
        self.assertEqual(plan["MinimumMemAvailableBytes"], plan["Instance"]["MinimumMemAvailableBytes"])
        self.assertEqual(plan["MinimumMemAvailableBytes"], int(0.6 * 768 * 1024 ** 3))
        self.assertIn("#PBS -l instance_type=m8g.48xlarge\n", (root / case_id / "main" / "job.pbs").read_text())
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

        def fake_upload(record, context, *, remote, profile, adopt_remote_case=False):
            events.append(("upload", record["Case"]))
            return {"Commands": [], "UTC": "fake"}

        def fake_submit(host, pbs_bin, script, cwd, *, job_cap, user=None):
            case = Path(cwd).parts[-2]
            events.append(("submit", case))
            return {"Job": f"{len(events)}.fake", "UTC": qualify_library.remote_side.utc(), "UserJobsBefore": 0, "JobCap": job_cap,
                    "Command": "qsub"}

        def fake_job_exit(host, remote_job_directory, job_id):
            events.append(("exit-check", Path(remote_job_directory).parts[-2]))
            return {"OK": True, "Reasons": [], "ExitRecord": {"ExitCode": 0, "JobID": job_id}, "StatusPBSJobID": job_id}

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
            # The replayed status.json is stamped with this run's job id (the recorded one would read as stale).
            status_path = Path(local_directory) / "status.json"
            submission = json.loads((Path(local_directory).parents[1] / "submission.json").read_text())
            status_path.write_text(json.dumps({**json.loads(status_path.read_text()), "PBSJobID": submission["Job"]}))
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
                 "delete_archives": fake_delete, "qstat_history": lambda host, pbs_bin, job: "job_state = F", "job_exit": fake_job_exit}
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

        def fake_upload(record, context, *, remote, profile, adopt_remote_case=False):
            events.append(("upload", record["Case"]))
            return {"Commands": [], "UTC": "fake"}

        def fake_job_exit(host, remote_job_directory, job_id):
            events.append(("exit-check", Path(remote_job_directory).name))
            return {"OK": True, "Reasons": [], "ExitRecord": {"ExitCode": 0, "JobID": job_id}, "StatusPBSJobID": job_id}

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
                 "read_json": fake_read_json, "count_archive_potentials": fake_count, "job_exit": fake_job_exit}
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
            # A previous run's kept model: the same source model under another name, its
            # fabricated matrices recorded absolute (the campaign's reducer), no thin matrices.
            other = dict(json.loads(case_library.read_text())["Models"][0])
            reducer = ASSESSMENT / spec["Campaign"] / "results" / "main" / f"{spec['Prefix']}-p4" / "reducer"
            other.update({"Name": "other-coupon", "LibraryQualified": False,
                          "FabricatedMatrix": str(reducer / "domain-response-matrix.csv"),
                          "FabricatedSurfaceMatrix": str(reducer / "surface-response-matrix.csv"),
                          "Qualification": {"Verdict": "PendingQualification", "Record": "/tmp/previous/other/qualification.json"},
                          "SourceProcessLibrary": {"Path": str(case_library), "SHA256": sha256(case_library)}})
            previous_library = self.tmp / "previous-process-library.json"
            previous_library.write_text(json.dumps({"Version": 1, "Root": "/tmp/previous", "Models": [
                other, {"Name": model_name, "LibraryQualified": False, "Stale": True}]}))
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
        self.assertEqual(merged["Version"], 3)
        self.assertEqual(merged["Models"][0]["FabricatedMatrix"], "models/other-coupon/fabricated-domain-response-matrix.csv")
        self.assertIsNone(merged["Models"][0]["ThinMatrix"])
        self.assertEqual(merged["Loadable"]["NotLoadable"], ["other-coupon", model_name])
        preflight = json.loads((self.tmp / "split" / "process-library-preflight.json").read_text())
        self.assertTrue(preflight["PreflightOnly"])
        self.assertEqual(preflight["Models"][1]["ThinMatrix"], f"models/{qualify_library.model_slug(model_name)}/thin-domain-response-matrix.csv")
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

    def test_controls_only_two_cases_take_their_own_reducers(self):
        """Decision 495 (1): TWO cases in ONE stored root (a --reference none dry run of the four-edge and the
        three-edge cases, completed with their campaigns' recorded reducer CSVs): ONE controls-only run of both
        takes EACH coupon's own reducer as its prior main stage (ControlAmplitudes.PriorMainStage =
        ReusedMain.Reducer.Directory) - the first coupon's reducer leaks into no other coupon; the shared args
        Namespace is never written (the in-process prepare_case sequence of run_qualify)."""
        cases = {"four-edge-9d2cb9bbb3fe": "four-edge-physics-11", "three-edge-419576fdab24": "gallery-physics-06b"}
        prefixes = {"four-edge-9d2cb9bbb3fe": "va", "three-edge-419576fdab24": "g06b"}   # the campaigns' recorded stage prefixes
        stored = self.tmp / "stored-two"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, *[token for case_id in cases for token in ("--case", case_id)],
                   "--root", str(stored), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record_path = stored / "library-qualification.json"
        record = json.loads(record_path.read_text())
        self.assertEqual([case["Status"] for case in record["Cases"]], ["planned", "planned"])
        own = {}
        for case, campaign in zip(record["Cases"], (ASSESSMENT / cases[c] for c in cases)):
            case_id = case["Case"]
            reducer = stored / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"
            reducer.mkdir(parents=True)
            digests = {}
            for name in qualify_library.REUSED_MAIN_CSVS:
                shutil.copyfile(campaign / "results" / "main" / f"{prefixes[case_id]}-p4" / "reducer" / name, reducer / name)
                digests[name] = sha256(reducer / name)
            case["ResultDigests"] = {f"main/{case_id}-p4/reducer/{name}": {"Local": digest, "Remote": digest, "OK": True}
                                     for name, digest in digests.items()}
            case["Jobs"][0]["Submission"] = {"Job": f"{case_id}.fake"}
            case["Qualification"] = {"Verdict": "Failed", "UnjudgedTypes": []}
            case["Cost"] = {"MainStage": {"NodeHours": 1.0}, "JobNodeHours": 1.5}
            own[case_id] = str(reducer)
        record_path.write_text(json.dumps(record, indent=2))
        self.assertNotEqual(own["four-edge-9d2cb9bbb3fe"], own["three-edge-419576fdab24"])
        # One controls-only run of both coupons.
        root = self.tmp / "controls-only-two"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, *[token for case_id in cases for token in ("--case", case_id)],
                   "--controls", "p3,p5", "--controls-only", "--reuse-main", str(stored), "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        planned = {case["Case"]: case for case in json.loads((root / "library-qualification.json").read_text())["Cases"]}
        for case_id in cases:
            case = planned[case_id]
            self.assertEqual(case["Status"], "planned", case.get("StoppedBy"))
            self.assertEqual(case["ControlAmplitudes"]["PriorMainStage"], own[case_id], case_id)
            self.assertEqual(case["ReusedMain"]["Reducer"]["Directory"], own[case_id])
            self.assertEqual(case["ReusedMain"]["Case"], case_id)
            self.assertIn("decision 474 (A)", case["Controls"]["Rule"])
            self.assertEqual(len(case["Controls"]["Indices"]), 8)
        self.assertNotEqual(planned["four-edge-9d2cb9bbb3fe"]["Controls"]["Indices"], planned["three-edge-419576fdab24"]["Controls"]["Indices"])
        # The same two coupons through ONE prepare_case sequence sharing ONE args Namespace: args untouched.
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        args = parser.parse_args(["--build-record", str(self.build_record), "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256,
                                  "--controls", "p3,p5", "--controls-only", "--reuse-main", str(stored), "--root", str(root / "shared"), "--dry-run"])
        build = json.loads(self.build_record.read_text())
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        cost_model = estimate_stages.load_cost_model(args.cost_model)
        table, digest = gates.load_gates(args.gates)
        (root / "shared").mkdir(parents=True, exist_ok=True)
        before = dict(vars(args))
        priors = {}
        for case_id in cases:
            case_record = next(item for item in build["Cases"] if item["Case"] == case_id)
            planned_record, _ = qualify_library.prepare_case(case_record, manifest_path=MANIFEST, manifest=manifest, args=args, root=root / "shared",
                                                             remote=None, profile=profile, cost_model=cost_model, gates=table, gates_digest=digest)
            priors[case_id] = planned_record["ControlAmplitudes"]["PriorMainStage"]
            self.assertEqual(dict(vars(args)), before)
        self.assertEqual(priors, own)

    def test_controls_only_reuses_the_stored_main_stage_and_plans_the_control_job(self):
        """--controls-only --reuse-main (decisions 474 (A) / 479 / 485 (c)) on a stored run synthesized from a
        --reference none dry run of the four-edge case with the recorded physics-11 reducer CSVs: the
        identity checks pass (the same mesh, config and traces), the controls are the amplitude-informed
        choice from the stored reducer, one controls-only job of the fixed stages is planned, and the
        reused main stage is recorded with the stored reducer's digests and partition."""
        case_id = "four-edge-9d2cb9bbb3fe"
        campaign = ASSESSMENT / CASES[case_id]["Campaign"]
        stored = self.tmp / "stored-run"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, "--case", case_id, "--stage-prefix", "va",
                   "--root", str(stored), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        # The stored run: its planned record completed by hand with the recorded reducer CSVs, their digests
        # and a single job (what an analyzed run records).
        record_path = stored / "library-qualification.json"
        record = json.loads(record_path.read_text())
        case = record["Cases"][0]
        reducer = stored / case_id / "results" / "main" / "va-p4" / "reducer"
        reducer.mkdir(parents=True)
        digests = {}
        for name in qualify_library.REUSED_MAIN_CSVS:
            shutil.copyfile(campaign / "results" / "main" / "va-p4" / "reducer" / name, reducer / name)
            digests[name] = sha256(reducer / name)
        case["ResultDigests"] = {f"main/va-p4/reducer/{name}": {"Local": digest, "Remote": digest, "OK": True} for name, digest in digests.items()}
        case["Jobs"][0]["Submission"] = {"Job": "1.fake"}
        case["Qualification"] = {"Verdict": "Failed", "UnjudgedTypes": []}
        case["Cost"] = {"MainStage": {"NodeHours": 1.0}, "JobNodeHours": 1.5}
        record_path.write_text(json.dumps(record, indent=2))
        root = self.tmp / "controls-only"
        command = [sys.executable, str(HERE / "coupon_library.py"), "qualify", "--build-record", str(self.build_record),
                   "--reference", "none", "--frozen-binary-sha256", BINARY_SHA256, "--case", case_id, "--stage-prefix", "va",
                   "--controls", "p3,p5", "--controls-only", "--reuse-main", str(stored), "--root", str(root), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        planned = json.loads((root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(planned["Status"], "planned", planned.get("StoppedBy"))
        reused = planned["ReusedMain"]
        self.assertEqual((reused["Root"], reused["MainPrefix"], reused["Case"]), (str(stored), "va-p4", case_id))
        self.assertEqual(reused["Reducer"]["SHA256"], digests)
        self.assertEqual((reused["Reducer"]["Job"], reused["Reducer"]["PBSJobID"], reused["Reducer"]["Nodes"]), ("single", "1.fake", 1))
        self.assertEqual(reused["Traces"], {"Count": 80, "PinnedByStoredPlans": 80, "Identical": True})
        self.assertEqual(reused["StoredVerdict"], "Failed")
        self.assertEqual(reused["ReusedStages"], ["va-p4"])
        self.assertEqual(reused["SolvedStages"], ["va-p5-control", "va-p3-control", "va-p4-local-edge"])
        self.assertEqual(planned["ControlAmplitudes"]["PriorMainStage"], str(reducer))
        self.assertIn("decision 474 (A)", planned["Controls"]["Rule"])
        self.assertEqual(len(planned["Controls"]["Indices"]), 8)
        self.assertEqual((planned["Split"]["N"], planned["Split"]["ControlsJob"], planned["Split"]["Blocks"]), (1, "controls-only", [0]))
        self.assertEqual([(job["Name"], job["Kind"], job["Sources"]) for job in planned["Jobs"]], [("controls-only", "controls-only", [])])
        plan = json.loads(Path(planned["Jobs"][0]["Plan"]).read_text())
        self.assertEqual(plan["JobKind"], "controls-only")
        self.assertEqual(plan["StageNames"], ["va-p5-control-worker", "va-p5-control-reducer", "va-p3-control-worker", "va-p3-control-reducer",
                                              "va-p4-local-edge"])
        self.assertEqual(plan["ReusedMain"]["Reducer"]["SHA256"], digests)
        self.assertNotIn(f"{planned['Remote']['Case']}/main/va-p4/worker.json", plan["PinnedSHA256"])
        self.assertEqual(sum(1 for path in plan["PinnedSHA256"] if "/inputs/traces/" in path), 80)
        self.assertTrue((root / case_id / "main" / "jobs" / "controls-only" / "job.pbs").is_file())
        self.assertEqual(planned["Nodes"]["MainOrigin"], "the stored run's reducer job (--reuse-main)")
        # A stored run of another mesh is refused (ReuseMain, the coupon stops before any plan).
        record["Cases"][0]["Mesh"]["SHA256"] = "0" * 64
        record_path.write_text(json.dumps(record, indent=2))
        other = self.tmp / "controls-only-other-mesh"
        command[command.index(str(root))] = str(other)
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        stopped = json.loads((other / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual((stopped["Status"], stopped["StoppedBy"]["Kind"]), ("failed", "ReuseMain"))
        self.assertIn("not the same coupon mesh", stopped["StoppedBy"]["Message"])


if __name__ == "__main__":
    unittest.main()
