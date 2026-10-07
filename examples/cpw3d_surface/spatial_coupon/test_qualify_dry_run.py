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
stops.  Needs a local identity mesh of each case and the assessment tree.

The fixture is shell-aware (decision 500 (4)): every identity mesh the current mesher publishes
carries the decision-61a per-ring MA shells, so the synthetic build record binds the root's
RadialShells census (the identity receipt's, the census file beside the mesh); the byte-identity
tests compare the BASE (pre-shell-expansion) stage configs with the stored pre-shell campaign
configs and the shell expansion separately against the decision-61a rule.  The result /
scheduler / controls-only replays run on shell-era stored roots of (b) batch 1 (Batch1ReplayTest)."""
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
import build_configs  # noqa: E402
import build_plan  # noqa: E402
import case_inputs  # noqa: E402
import compare_split_matrices  # noqa: E402
import estimate_stages  # noqa: E402
import gates  # noqa: E402
import job_split  # noqa: E402
import qualify_library  # noqa: E402
import summarize_cost  # noqa: E402
from mixed_mesh import SHELL_LABEL_STRIDE, h1_dofs_from_counts  # noqa: E402
from run_gmsh_only_matrix import radial_shells_record, sha256  # noqa: E402

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


SHELL_CENSUS_NAME = "identity.msh.radial-shells.json"


def local_identity_mesh(case_id):
    """The newest local Gmsh-only root of the case with an identity mesh (None when absent)."""
    roots = sorted(glob.glob(f"/tmp/coupon-gmsh-only-{case_id}-*/identity.msh"), key=os.path.getmtime)
    return Path(roots[-1]) if roots else None


def shell_census_binding(mesh):
    """The build record's RadialShells binding of a local identity mesh (decision 61a): the root's
    identity receipt (run_gmsh_only_matrix.radial_shells_record) with Shells.Path re-pointed to the
    census file beside the mesh, which must carry the receipt's SHA256 and bind this mesh.  None
    when the root carries no receipt / census (a root rebuilt without its records)."""
    root = Path(mesh).parent
    binding = radial_shells_record(root)
    census_path = root / SHELL_CENSUS_NAME
    if binding is None or not census_path.is_file():
        return None
    if sha256(census_path) != binding["Shells"]["SHA256"]:
        raise AssertionError(f"{census_path} is not the census the identity receipt binds ({binding['Shells']['SHA256']})")
    census = json.loads(census_path.read_text())
    if census["Mesh"]["SHA256"] != sha256(mesh):
        raise AssertionError(f"{census_path} binds mesh {census['Mesh']['SHA256']}, not {mesh}")
    return {**binding, "Shells": {"Path": str(census_path), "SHA256": binding["Shells"]["SHA256"]}}


def available():
    return all(local_identity_mesh(case_id) is not None and shell_census_binding(local_identity_mesh(case_id)) is not None
               and (ASSESSMENT / spec["Campaign"] / "main" / "plan.json").is_file() for case_id, spec in CASES.items())


def snapshot(directory):
    """mtime of every file of a stored directory (a recorded campaign, a stored qualify root)."""
    return {path: path.stat().st_mtime_ns for path in Path(directory).rglob("*") if path.is_file()}


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
        # Decision 513 (1): the decision-500 (3) usage restriction is recorded as LIFTED (after the S3p / S2p e2e),
        # in the record and in the rule; the restriction wording is gone.
        self.assertEqual(reused["Usage"], job_split.CONTROLS_ONLY_USAGE_RECORD)
        self.assertEqual((reused["Usage"]["Restriction"], reused["Usage"]["LiftedBy"]), ("decision 500 (3)", "decision 513 (1)"))
        self.assertIn("LIFTED by decision 513", reused["Rule"])
        self.assertNotIn("USAGE RESTRICTION", reused["Rule"])
        self.assertNotIn("UsageRestriction", reused)
        self.assertFalse(hasattr(job_split, "CONTROLS_ONLY_USAGE_RESTRICTION"))
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


@unittest.skipUnless(available(), "local identity meshes with their radial-shell census and the assessment campaigns are needed")
class QualifyDryRunTest(unittest.TestCase):
    """The gallery cases on their local identity meshes: the one-node plans against the stored pre-shell
    campaigns (base configs byte-identical, the shell expansion by the decision-61a rule), the control
    choice, the fail-closed stops and the derived-config identity."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coupon-qualify-test-"))
        cases = []
        for case_id in CASES:
            mesh = local_identity_mesh(case_id)
            counts = estimate_stages.entity_counts_of_mesh(mesh)
            cases.append({"Case": case_id, "Status": "built", "Passed": True, "CanonicalBuildId": None,
                          "Variants": {"identity": {"Path": str(mesh), "SHA256": sha256(mesh)}},
                          "RadialShells": shell_census_binding(mesh),
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
        # The mesh carries the decision-61a shells: the record binds the census, the run inputs expand
        # the base config (parent labels, the reference comparison view) into one MA entry per shell.
        census = json.loads(Path(case["Mesh"]["RadialShells"]["Shells"]["Path"]).read_text())
        self.assertEqual((case["Mesh"]["RadialShells"]["Binding"], case["Mesh"]["RadialShells"]["ShellCount"]), ("RadialShells", 15))
        self.assertEqual(census["Mesh"]["SHA256"], case["Mesh"]["SHA256"])
        base_config = json.loads((root / case_id / "inputs" / "base-run-config.json").read_text())
        run_config = json.loads((root / case_id / "inputs" / "run-config.json").read_text())
        self.assertEqual(case["Inputs"]["BaseConfigSHA256"], sha256(root / case_id / "inputs" / "base-run-config.json"))
        self.assertEqual(case["Reference"]["ConfigEqualsDerived"]["ComparedConfig"],
                         "the base config (parent labels) of the radial-shell relabel")
        self.assert_shell_expansion_follows_decision_61a(base_config, run_config, census, case["Inputs"]["RadialShells"])
        # The BASE stage configs (derived from the base run config exactly as the driver derives the
        # stage configs from the run config) equal the stored pre-shell campaign's apart from paths;
        # the written (shell-expanded) stage configs are the base stage configs expanded by the rule.
        plan = json.loads((root / case_id / "main" / "plan.json").read_text())
        remote_case, remote_mesh, remote_traces = case["Remote"]["Case"], plan["MeshRemote"], f"{case['Remote']['Case']}/inputs/traces"
        indices = case["Sources"]["Indices"]
        for stage in spec["Stages"]:
            order = int(stage.rsplit("-p", 1)[1].split("-")[0])
            subset = indices if stage == spec["Stages"][0] else spec["Controls"]
            if stage.endswith("local-edge"):
                config, _ = build_configs.derive(base_config, remote_mesh, f"{remote_case}/main/{stage}", subset, remote_traces,
                                                 order=order, save_local_edge_energy=True)
                config["Problem"]["Output"] = f"{remote_case}/main/{stage}/output"
                base_stages = {"config.json": config}
            else:
                worker, reducer = build_configs.derive(base_config, remote_mesh, f"{remote_case}/main/{stage}", subset, remote_traces,
                                                       order=order)
                base_stages = {"worker.json": worker, "reducer.json": reducer}
            for name, base_stage in base_stages.items():
                recorded = strip_paths(json.loads((campaign / "main" / stage / name).read_text()))
                self.assertEqual(strip_paths(base_stage), recorded, f"{stage}/{name}: base stage config vs the stored campaign")
                generated = json.loads((root / case_id / "main" / stage / name).read_text())
                expanded, _ = case_inputs.expand_radial_shells(base_stage, census)
                self.assertEqual(generated, expanded, f"{stage}/{name}: the written stage config vs the expanded base stage config")
                self.assertNotEqual(generated, base_stage)
        # The plan equals the recorded one in stages (names, config file, environment,
        # dependencies, order) and pins the same files: every trace name of the recorded plan
        # (the traces are regenerated from the case's basis - their digests are the run's
        # own, recorded under Inputs.Sources); the mesh pin is the build record's.
        recorded_plan = json.loads((campaign / "main" / "plan.json").read_text())

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

    def assert_shell_expansion_follows_decision_61a(self, base_config, run_config, census, shell_record):
        """Decision 61a: the run config is the base config with every attribute list naming a parent MA
        label naming its shell labels (10000 x ordinal + parent; ordinal 1 far, 1 + k / 1 + K + k the top /
        bottom ring k), and the one MA Dielectric entry replaced by one entry per shell ordinal (fresh
        indices after the base maximum, the same Type / layer / edge settings); every other field is
        untouched; the recorded shell map names every shell entry with its ring radii."""
        rings = len(census["RingRadii"])
        parents = sorted({int(shell["Parent"]) for shell in census["Shells"]})
        self.assertEqual(len(census["Shells"]), (1 + 2 * rings) * len(parents))
        self.assertEqual(census["LabelStride"], SHELL_LABEL_STRIDE)
        for shell in census["Shells"]:
            self.assertEqual(int(shell["Label"]), SHELL_LABEL_STRIDE * int(shell["Ordinal"]) + int(shell["Parent"]))
        base_entries = base_config["Boundaries"]["Postprocessing"]["Dielectric"]
        run_entries = run_config["Boundaries"]["Postprocessing"]["Dielectric"]
        ma_entries = [entry for entry in base_entries if entry["Type"] == "MA"]
        self.assertEqual(len(ma_entries), 1)
        self.assertEqual(sorted(int(a) for a in ma_entries[0]["Attributes"]), parents)
        others = [entry for entry in base_entries if entry["Type"] != "MA"]
        self.assertEqual([entry for entry in run_entries if entry["Type"] != "MA"], others)
        shells = [entry for entry in run_entries if entry["Type"] == "MA"]
        self.assertEqual(len(shells), 1 + 2 * rings)
        next_index = max(int(entry["Index"]) for entry in base_entries) + 1
        self.assertEqual([int(entry["Index"]) for entry in shells], list(range(next_index, next_index + len(shells))))
        for ordinal, entry in enumerate(shells, start=1):
            self.assertEqual(sorted(entry["Attributes"]), sorted(SHELL_LABEL_STRIDE * ordinal + parent for parent in parents))
            self.assertEqual({key: value for key, value in entry.items() if key not in ("Index", "Attributes")},
                             {key: value for key, value in ma_entries[0].items() if key not in ("Index", "Attributes")})
            recorded = shell_record["Interfaces"][str(entry["Index"])]
            self.assertEqual((recorded["Type"], recorded["Ordinal"], recorded["BaseIndex"]), ("MA", ordinal, int(ma_entries[0]["Index"])))
            if ordinal == 1:
                self.assertEqual(recorded["Kind"], "far")
            else:
                ring = (ordinal - 2) % rings + 1
                self.assertEqual((recorded["Kind"], recorded["Ring"]), ("top" if ordinal <= 1 + rings else "bottom", ring))
                self.assertEqual(recorded["OuterRadius"], census["RingRadii"][ring - 1])
        self.assertEqual(shell_record["RingRadii"], census["RingRadii"])
        # Ground names the shells in place of the parents; nothing outside Boundaries changes.
        self.assertEqual([a for a in run_config["Boundaries"]["Ground"]["Attributes"] if int(a) < SHELL_LABEL_STRIDE],
                         [a for a in base_config["Boundaries"]["Ground"]["Attributes"] if int(a) not in parents])
        self.assertEqual(sorted(a for a in run_config["Boundaries"]["Ground"]["Attributes"] if int(a) >= SHELL_LABEL_STRIDE),
                         sorted(int(shell["Label"]) for shell in census["Shells"]))
        self.assertEqual({key: value for key, value in run_config.items() if key != "Boundaries"},
                         {key: value for key, value in base_config.items() if key != "Boundaries"})
        expanded, shell_map = case_inputs.expand_radial_shells(base_config, census)
        self.assertEqual(expanded, run_config)
        self.assertEqual({str(index): value for index, value in shell_map.items()}, shell_record["Interfaces"])

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

    @unittest.skipUnless(local_identity_mesh(GALLERY_10["Case"]) is not None
                         and (ASSESSMENT / GALLERY_10["Campaign"] / "results" / "main" / "status.json").is_file(),
                         "the two-edge mesh and the gallery-10 campaign are needed (no local two-edge root was rebuilt; the stored "
                         "gallery-10 results predate the decision-61a shells, so its analysis half needs a shell-era record too)")
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
                           "Variants": {"identity": {"Path": str(mesh), "SHA256": sha256(mesh)}}, "RadialShells": shell_census_binding(mesh),
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
        before = snapshot(campaign)
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
        self.assertEqual(snapshot(campaign), before, "the recorded campaign directory must stay read only")
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
        """Decision 52: for every gallery case with a local identity mesh the BASE config derived
        from the case's own sources (process library, trace basis, mesh $PhysicalNames, recipe
        PhysicsRun; the parent MA labels of the mesh's decision-61a shell census) equals the
        config the graded_v2 reference ran apart from Model.Mesh, Problem.Output and the
        DataFile directory; Solver.Order / Linear.Tol are the recipe's (the references at p5 /
        Tol 1e-8 - gallery 10 and the ten-edge - differ there only); the run config is the base
        expanded by the census.  The two cases of CASES are required; the others are checked
        when their mesh and reference config are present (the gap is the explicit skip of
        test_every_gallery_case_has_a_local_root_and_reference_config)."""
        manifest = json.loads(MANIFEST.read_text())
        manifest["Path"] = str(MANIFEST)
        physics_run = case_inputs.physics_run_parameters(manifest)
        self.assertEqual((physics_run["Order"], physics_run["LinearTol"]), (4, 1e-10))
        checked, unchecked = [], []
        for case_id, reference_path in GALLERY_REFERENCE_CONFIGS.items():
            mesh = local_identity_mesh(case_id)
            if mesh is None or not reference_path.is_file():
                unchecked.append(case_id)
                continue
            case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
            directory = qualify_library.source_directory(MANIFEST, manifest, case)
            out = self.tmp / "derived" / case_id
            census_binding = shell_census_binding(mesh)
            self.assertIsNotNone(census_binding, f"{mesh} carries no radial-shell census beside it")
            census = json.loads(Path(census_binding["Shells"]["Path"]).read_text())
            run_config, record = case_inputs.derive(case, directory, mesh_path=mesh, physics_run=physics_run, out_dir=out,
                                                    radial_shells=census)
            config = record["BaseConfig"]
            self.assertEqual(case_inputs.expand_radial_shells(config, census)[0], run_config)
            reference = json.loads(reference_path.read_text())
            self.assertEqual(case_inputs.config_differences(config, reference, ignore_solver=("Order", "Linear.Tol")), [], case_id)
            self.assertEqual(config["Solver"]["Order"], 4)
            self.assertEqual(config["Solver"]["Linear"]["Tol"], 1e-10)
            self.assertEqual(len(config["Boundaries"]["PrescribedPotential"]), len(reference["Boundaries"]["PrescribedPotential"]))
            self.assertEqual(record["BaseInterfaces"], case_inputs.interface_types(reference))
            self.assertEqual(len(record["Interfaces"]), len(record["BaseInterfaces"]) - 1 + census_binding["ShellCount"])
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
        self.assertTrue(set(CASES) <= set(checked), checked)
        self.assertEqual(sorted(checked + unchecked), sorted(GALLERY_REFERENCE_CONFIGS))

    def test_every_gallery_case_has_a_local_root_and_reference_config(self):
        """The coverage marker of the derived-config identity (decision 513 (3)): PASSES only when every one of
        the five gallery cases has a local identity root with its census and a reference config; otherwise an
        explicit SKIP naming the unchecked cases (today the two-edge / ten-edge roots: decision 496 rebuilt the
        four-edge and three-edge roots only), so the gap stays visible in the counts."""
        unchecked = {}
        for case_id, reference_path in GALLERY_REFERENCE_CONFIGS.items():
            mesh = local_identity_mesh(case_id)
            missing = [name for name, present in (("local identity root", mesh is not None),
                                                  ("radial-shell census", mesh is not None and shell_census_binding(mesh) is not None),
                                                  ("reference config", reference_path.is_file())) if not present]
            if missing:
                unchecked[case_id] = missing
        if unchecked:
            self.skipTest(f"derived-config identity unchecked for {unchecked}: the other {len(GALLERY_REFERENCE_CONFIGS) - len(unchecked)} "
                          f"gallery cases are checked by test_run_config_derived_from_the_case_equals_every_gallery_reference")

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
        # The run config carries the shells: the base MA entry and 15 shell entries (the analysis of a
        # --reference none run on stored results is Batch1ReplayTest's S4 replay).
        self.assertEqual([entry["Type"] for entry in worker["Boundaries"]["Postprocessing"]["Dielectric"]].count("MA"), 15)
        self.assertEqual(case["Inputs"]["RadialShells"]["RingRadii"], json.loads(Path(case["Mesh"]["RadialShells"]["Shells"]["Path"]).read_text())["RingRadii"])



# The shell-era replay inputs (decisions 500 (4) / 508): the (b) batch-1 stored roots of record (S4 one node, S3p /
# S2p split + re-qualified) with their identity meshes and reducer CSVs copied from the cluster lane directory,
# every file sha256-identical on both ends (never committed; the stored records keep their cluster paths, which the
# fixture relocates in its own copies - nothing under the root is written).
BATCH1 = Path(os.environ.get("COUPON_BATCH1_REPLAY_ROOT", "/tmp/tools-followups2-inputs/b-batch1"))
BATCH1_CLUSTER_LANE = "/data/home/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/curved-clusters-20261005/b-batch1"
BATCH1_BINARY_SHA256 = "cc7c4091fa47bde8739fe57678d1b72a22a4aee51b0de8da18ccb9912a5aec88"
BATCH1_WINDOWS = {
    "sct002-S4": {"Case": "spatial-3-edge-7391dde39d93", "Sources": 148, "Controls": [1, 2, 3, 4, 5, 25, 26, 75]},
    "sct002-S3p": {"Case": "spatial-8-edge-05c322f6cda8", "Sources": 258, "Controls": [1, 2, 3, 8, 10, 31, 32, 108],
                   "RequalifiedControls": [1, 2, 3, 15, 18, 32, 48, 108], "Model": "spatialedgecluster_edgecount-8_e7f44561bbf0"},
    "sct002-S2p": {"Case": "spatial-8-edge-818f8956d075", "Sources": 252, "Controls": [1, 2, 3, 7, 20, 31, 32, 93],
                   "RequalifiedControls": [1, 2, 3, 7, 20, 32, 33, 93], "Model": "spatialedgecluster_edgecount-8_d19e77e8fe71"},
}
# The decision-479 control sets of record of the two 8-edge re-qualifications.
BATCH1_REQUALIFICATION_CONTROLS = BATCH1 / "records" / "eight-edge-requalification-controls.json"


def batch1_stored_root(window, kind="fab"):
    """The stored qualify root of a window (qualify/<W>/fab-<case> or thin-<case>)."""
    case = BATCH1_WINDOWS[window]["Case"]
    return BATCH1 / "qualify" / window / f"{kind}-{case}"


def batch1_available():
    needed = [BATCH1_REQUALIFICATION_CONTROLS]
    for window, spec in BATCH1_WINDOWS.items():
        case = spec["Case"]
        needed += [BATCH1 / "registration" / window / "manifest-main.json",
                   BATCH1 / "registration" / window / "build" / "library-build.json",
                   BATCH1 / "registration" / window / "build" / case / "identity.msh",
                   batch1_stored_root(window) / "library-qualification.json",
                   batch1_stored_root(window) / case / "results" / "main" / f"{case}-p4" / "reducer" / "surface-response-matrix.csv"]
        if "RequalifiedControls" in spec:
            needed += [BATCH1 / "requalify" / window / "library-qualification.json",
                       BATCH1 / "requalify" / window / case / "results" / "main" / f"{case}-p5-control" / "reducer" / "surface-response-matrix.csv",
                       batch1_stored_root(window, "thin") / f"{case}-thin" / "results" / "main" / f"{case}-thin-p4" / "reducer" / "surface-response-matrix.csv",
                       BATCH1 / "f-qualification" / case / "spatial-qualification-requal.json",
                       BATCH1 / "f-qualification" / case / "dense" / "dense-traces.json",
                       BATCH1 / "preflight" / window / "c0-p4" / "postpro" / "palace.json",
                       BATCH1 / "library" / f"b1-requal-interim-{window.split('-')[1]}" / "process-library.json"]
    needed.append(BATCH1 / "registration" / "sct002-S4" / "build" / f"{BATCH1_WINDOWS['sct002-S4']['Case']}-thin" / "identity.msh")
    return all(path.is_file() for path in needed)


def relocated(text):
    return text.replace(BATCH1_CLUSTER_LANE, str(BATCH1))


def assert_numbers_reproduced(test, stored, actual, path="", *, ignore=(), rel=1e-12):
    """Every value of a stored record is reproduced by the actual one: floats to `rel`, everything
    else equal; keys only in the actual record are allowed (the tool's newer fields); the paths in
    `ignore` (dotted) are not compared (host-specific paths, UTC stamps)."""
    if path in ignore:
        return
    if isinstance(stored, dict):
        test.assertIsInstance(actual, dict, path)
        for key, value in stored.items():
            test.assertIn(key, actual, f"{path}.{key}")
            assert_numbers_reproduced(test, value, actual[key], f"{path}.{key}" if path else key, ignore=ignore, rel=rel)
    elif isinstance(stored, list):
        test.assertEqual(len(stored), len(actual), path)
        for index, (left, right) in enumerate(zip(stored, actual)):
            assert_numbers_reproduced(test, left, right, f"{path}[{index}]", ignore=ignore, rel=rel)
    elif isinstance(stored, float) and isinstance(actual, (int, float)) and not isinstance(actual, bool):
        test.assertTrue(abs(stored - actual) <= rel * max(abs(stored), abs(actual), 1e-300) or stored == actual, f"{path}: {stored} != {actual}")
    else:
        test.assertEqual(stored, actual, path)


@unittest.skipUnless(batch1_available(), f"the (b) batch-1 replay inputs under {BATCH1} are needed (decision 508)")
class Batch1ReplayTest(unittest.TestCase):
    """The shell-era replays on the (b) batch-1 stored roots: the analysis of the stored S4 results
    reproduces its record, the concurrent scheduler against a fake remote replaying S4 fab + thin,
    the split-coupon run replaying the S3p 3-job root, and the controls-only re-qualification of
    S3p / S2p (decisions 474 (A) / 479 / 485 (c)) end to end against the records of record.  The
    stored roots are read only (mtimes snapshotted); the driver runs its real paths (the mesh
    hashed, the identity checks, the control choice, the analysis, the (F) evaluate)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coupon-batch1-replay-"))
        cls.build_records = {}
        for window in BATCH1_WINDOWS:
            local = cls.tmp / window
            local.mkdir()
            manifest = local / "manifest-main.json"
            manifest.write_text(relocated((BATCH1 / "registration" / window / "manifest-main.json").read_text()))
            build = json.loads(relocated((BATCH1 / "registration" / window / "build" / "library-build.json").read_text()))
            build["Library"]["Manifest"] = {"Path": str(manifest), "SHA256": sha256(manifest)}
            (local / "library-build.json").write_text(json.dumps(build, indent=2))
            cls.build_records[window] = local / "library-build.json"
        cls.before = snapshot(BATCH1)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, True)

    def tearDown(self):
        self.assertEqual(snapshot(BATCH1), self.before, "the stored batch-1 roots must stay read only")

    def args(self, window, root, *cases, extra=()):
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        tokens = ["--build-record", str(self.build_records[window]), "--reference", "none", "--remote", "h:/r", "--max-jobs", "2",
                  "--job-policy", "frugal", "--frozen-binary-sha256", BATCH1_BINARY_SHA256, "--root", str(root), "--monitor-interval", "0",
                  "--monitor-polls", "10", *[token for case in cases for token in ("--case", case)], *extra]
        return parser.parse_args(tokens)

    def build_case(self, window, case_id):
        return next(item for item in json.loads(self.build_records[window].read_text())["Cases"] if item["Case"] == case_id)

    def prepare(self, window, case_id, root, *, extra=()):
        args = self.args(window, root, case_id, extra=extra)
        build = json.loads(self.build_records[window].read_text())
        manifest = json.loads(Path(build["Library"]["Manifest"]["Path"]).read_text())
        manifest["Path"] = build["Library"]["Manifest"]["Path"]
        profile = json.loads(Path(args.cluster_profile).read_text())
        table, digest = gates.load_gates(args.gates)
        root.mkdir(parents=True, exist_ok=True)
        record, context = qualify_library.prepare_case(self.build_case(window, case_id), manifest_path=Path(manifest["Path"]),
                                                       manifest=manifest, args=args, root=root, remote={"Host": "h", "Root": "/r"},
                                                       profile=profile, cost_model=estimate_stages.load_cost_model(args.cost_model),
                                                       gates=table, gates_digest=digest)
        return record, context, table, digest, profile

    def assert_plan_reproduces_stored(self, record, stored_case, stored_plans, *, stage_configs):
        """The plan of this run equals the stored run's: controls, stage names, every regenerated trace
        digest = the stored plan's pin, the stage configs apart from the remote paths, the mesh."""
        self.assertEqual(record["Controls"]["Indices"], stored_case["Controls"]["Indices"])
        self.assertEqual(record["Plan"]["StageNames"], stored_case["Plan"]["StageNames"])
        self.assertEqual(record["Mesh"]["SHA256"], stored_case["Mesh"]["SHA256"])
        self.assertEqual(record["Sources"]["Count"], stored_case["Sources"]["Count"])
        self.assertEqual(record["Inputs"]["RadialShells"]["Interfaces"], stored_case["Inputs"]["RadialShells"]["Interfaces"])
        stored_pins = {}
        for plan in stored_plans:
            stored_pins.update({Path(key).name: value for key, value in plan["PinnedSHA256"].items() if "/inputs/traces/" in key})
        regenerated = {source["Name"]: source["SHA256"] for source in record["Inputs"]["Sources"]}
        self.assertEqual(len(stored_pins), stored_case["Sources"]["Count"])
        self.assertEqual(regenerated, stored_pins)
        for stage, names in stage_configs.items():
            for name in names:
                generated = strip_paths(json.loads((Path(record["Root"]) / "main" / stage / name).read_text()))
                stored = strip_paths(json.loads((Path(stored_case["Root"].replace(BATCH1_CLUSTER_LANE, str(BATCH1))) / "main" / stage / name)
                                                .read_text()))
                self.assertEqual(generated, stored, f"{stage}/{name}")

    def replay_remote(self, results_of, events, *, archive_count=None):
        """A fake remote replaying the stored results of each case: submit / poll / exit check / fetch
        (the stored results/main copied, every status.json stamped with this run's job id) / digests /
        archive deletion; `results_of` = case id -> stored results/main directory."""
        remote_side = qualify_library.remote_side

        def job_id_of(directory):
            return f"{Path(directory).name}.fake"

        def fake_submit(host, pbs_bin, script, cwd, *, job_cap, user=None):
            events.append(("submit", Path(cwd).name))
            return {"Job": job_id_of(cwd), "UTC": remote_side.utc(), "UserJobsBefore": 0, "JobCap": job_cap, "Command": "qsub"}

        def fake_job_exit(host, remote_job_directory, job_id):
            events.append(("exit-check", Path(remote_job_directory).name))
            return {"OK": True, "Reasons": [], "ExitRecord": {"ExitCode": 0, "JobID": job_id}, "StatusPBSJobID": job_id}

        def fake_poll(host, pbs_bin, job_id, status_path):
            events.append(("poll", Path(status_path).parts[-2]))
            return {"UTC": "fake", "JobState": "F", "QStat": "", "Status": None, "Reachable": True}

        def stamped(path, job_directory):
            status = json.loads(Path(path).read_text())
            return {**status, "PBSJobID": job_id_of(job_directory)}

        def fake_read_json(host, path):
            case = Path(path).parts[-5] if Path(path).parts[-3] == "jobs" else Path(path).parts[-3]
            events.append(("status", Path(path).parts[-2]))
            relative = Path(path).relative_to(Path(path).parents[2]) if Path(path).parts[-3] == "jobs" else Path(path).name
            return stamped(results_of[case] / relative, Path(path).parent)

        def fake_fetch(host, remote_directory, local_directory):
            case = Path(remote_directory).parts[-2]
            events.append(("fetch", case))
            local = Path(local_directory)
            shutil.copytree(results_of[case], local, dirs_exist_ok=True)
            for status_path in [local / "status.json"] + sorted(local.glob("jobs/*/status.json")):
                if status_path.is_file():
                    status_path.write_text(json.dumps(stamped(status_path, status_path.parent)))
            return ["rsync", "fake"]

        def fake_remote_sha256(host, paths):
            digests = {}
            for path in paths:
                case = next(case for case in results_of if f"/{case}/main/" in path)
                digests[path] = sha256(Path(self.run_root) / case / "results" / "main" / path.split("/main/", 1)[1])
            return digests

        def fake_delete(host, archives):
            events.append(("delete", tuple(Path(a).parts[-2] for a in archives)))
            return {"Archives": list(archives), "SizesBeforeDeletion": "0", "DeletedUTC": "fake", "Remaining": ""}

        def fake_count(host, directory):
            events.append(("archive-count", Path(directory).parts[-2]))
            return archive_count

        def fake_upload(record, context, *, remote, profile, adopt_remote_case=False):
            events.append(("upload", record["Case"]))
            return {"Commands": [], "UTC": "fake"}
        fakes = {"submit": fake_submit, "poll": fake_poll, "fetch": fake_fetch, "remote_sha256": fake_remote_sha256,
                 "delete_archives": fake_delete, "qstat_history": lambda host, pbs_bin, job: "job_state = F", "job_exit": fake_job_exit,
                 "read_json": fake_read_json, "count_archive_potentials": fake_count}
        saved = {name: getattr(remote_side, name) for name in fakes}
        saved_upload = qualify_library.upload_case
        for name, fake in fakes.items():
            setattr(remote_side, name, fake)
        qualify_library.upload_case = fake_upload

        def restore():
            for name, fake in saved.items():
                setattr(remote_side, name, fake)
            qualify_library.upload_case = saved_upload
        self.addCleanup(restore)
        return restore

    def test_analysis_of_the_recorded_results_reproduces_the_verdicts(self):
        """The S4 fab coupon (one node, 148 sources, --reference none): prepare_case on the stored build
        record reproduces the stored plan (controls by class, stage names, every trace pin, the stage
        configs, the shell map), and analyze_case on the stored results reproduces the stored record:
        the verdict (PendingQualification: controls passed, no reference), every p-sequence step of every
        control, the MA_sharp tail summary and the main-stage cost, the stored root untouched."""
        window, spec = "sct002-S4", BATCH1_WINDOWS["sct002-S4"]
        case_id = spec["Case"]
        stored_root = batch1_stored_root(window)
        stored_case = json.loads((stored_root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(stored_case["Controls"]["Indices"], spec["Controls"])
        record, context, table, digest, profile = self.prepare(window, case_id, self.tmp / "s4-analysis")
        self.assertEqual(record["Mesh"]["RadialShells"]["Binding"], "RadialShells")
        self.assertEqual(record["Mesh"]["RadialShells"]["ShellCount"], 15)
        self.assertEqual((record["Split"]["N"], record["JobPolicy"]["Mode"]), (1, "frugal"))
        self.assertEqual((record["ReducerBlockSize"]["Value"], record["ReducerBlockSize"]["Origin"]), (48, "manifest ProductionRecipe.PhysicsRun.ReducerBlockSize"))
        self.assert_plan_reproduces_stored(record, stored_case, [json.loads((stored_root / case_id / "main" / "plan.json").read_text())],
                                           stage_configs={f"{case_id}-p4": ("worker.json", "reducer.json"),
                                                          f"{case_id}-p5-control": ("worker.json", "reducer.json"),
                                                          f"{case_id}-p3-control": ("worker.json", "reducer.json"),
                                                          f"{case_id}-p4-local-edge": ("config.json",)})
        gate_record = qualify_library.analyze_case(record, context, stored_root / case_id / "results", gates=table, gates_digest=digest,
                                                   profile=profile)
        stored_gate = json.loads((stored_root / case_id / "qualification.json").read_text())
        self.assertEqual(gate_record["Verdict"], stored_gate["Verdict"])
        self.assertEqual(gate_record["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(gate_record["GatesPassed"], stored_gate["GatesPassed"])
        self.assertEqual(gate_record["Reason"], stored_gate["Reason"])
        self.assertEqual(gate_record["Interfaces"], stored_gate["Interfaces"])
        # Every p-step of every control to the stored digit; the stored run (gate table V5, before the
        # amplitude floor) judged every observable, the gate of record leaves a below-floor Type of a
        # control unjudged (decision 474 (A): Passed None, BelowAmplitudeFloor recorded) - never the other way.
        below_floor = self.assert_p_steps_reproduced(stored_gate, gate_record)
        judged = gate_record["Gates"]["PSequenceControls"]["JudgedControls"]
        self.assertEqual(judged["E"], 8)
        self.assertEqual({name: 8 - sum(1 for _, observable in below_floor if observable == name) for name in judged}, judged)
        self.assertEqual(gate_record["Gates"]["PSequenceControls"]["UnjudgedTypes"], [])
        self.assertEqual(record["Status"], "pending-qualification")
        self.assertIsNone(record["Cost"]["ReferenceNodeHours"])
        assert_numbers_reproduced(self, stored_case["Cost"]["MainStage"], record["Cost"]["MainStage"], "Cost.MainStage")
        self.assertAlmostEqual(record["Cost"]["MainStage"]["NodeHours"], 0.59996, places=5)
        # The MA_sharp tail (ma_tail.py unchanged since the stored run): the fitted estimators agree to the
        # libm level of the two hosts (the stored record was computed on the cluster login node; the Fit2-4
        # deficit quartiles differ at 3e-11 relative), the ring sums and raw / sharp participations exactly.
        stored_tail = json.loads((stored_root / case_id / "comparison" / "ma-tail.json").read_text())
        tail = json.loads((Path(record["Root"]) / "comparison" / "ma-tail.json").read_text())
        assert_numbers_reproduced(self, stored_tail["Orders"]["p4"]["Summary"], tail["Orders"]["p4"]["Summary"], "MATail.p4.Summary", rel=1e-9)
        assert_numbers_reproduced(self, stored_tail["Orders"]["p4"]["PerSource"], tail["Orders"]["p4"]["PerSource"], "MATail.p4.PerSource", rel=1e-9)
        for index, per_source in stored_tail["Orders"]["p4"]["PerSource"].items():
            self.assertEqual(per_source["Q_MA_raw"], tail["Orders"]["p4"]["PerSource"][index]["Q_MA_raw"], index)
        self.assertEqual(stored_tail["RingRadii"], tail["RingRadii"])
        manifest_path = self.build_records[window].parent / "manifest-main.json"
        library = qualify_library.process_library_entries([record], {case_id: context}, manifest_path=manifest_path,
                                                          manifest=json.loads(manifest_path.read_text()) | {"Path": str(manifest_path)},
                                                          root=Path(record["Root"]).parent)
        self.assertEqual(len(library["Models"]), 1)
        self.assertFalse(library["Models"][0]["LibraryQualified"])
        self.assertEqual(library["Models"][0]["Qualification"]["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(library["Models"][0]["CouponMesh"]["SHA256"], stored_case["Mesh"]["SHA256"])

    def assert_p_steps_reproduced(self, stored_gate, gate_record):
        """StepToHigherOrder / StepFromLowerOrder / Bound of every control observable equal the stored
        ones; Passed equal unless the gate of record left it unjudged below the amplitude floor; returns
        the (control, observable) pairs left unjudged."""
        below_floor = set()
        for control, observables in stored_gate["Gates"]["PSequenceControls"]["Controls"].items():
            for name, stored_step in observables.items():
                step = gate_record["Gates"]["PSequenceControls"]["Controls"][control][name]
                assert_numbers_reproduced(self, {key: stored_step[key] for key in ("StepToHigherOrder", "StepFromLowerOrder", "Bound")},
                                          step, f"control {control} {name}")
                if step["Passed"] is None:
                    self.assertFalse(step["Judged"], (control, name))
                    self.assertIn("BelowAmplitudeFloor", step, (control, name))
                    below_floor.add((control, name))
                else:
                    self.assertEqual(step["Passed"], stored_step["Passed"], (control, name))
        return below_floor

    def test_without_reference_matrices_the_verdict_is_pending(self):
        """--reference none (the batch-1 device coupons have no reference campaign): the S4 record carries the
        rule, the reference cost is None and the verdict can only be PendingQualification or Failed - the
        stored S4 fab run (controls passed) reads PendingQualification, the stored S4 thin run (--controls "",
        no p-sequence control) Failed; both reproduced from the stored results."""
        window, spec = "sct002-S4", BATCH1_WINDOWS["sct002-S4"]
        for kind, expected in (("fab", gates.VERDICT_PENDING), ("thin", gates.VERDICT_FAILED)):
            case_id = spec["Case"] if kind == "fab" else f"{spec['Case']}-thin"
            stored_root = batch1_stored_root(window, kind)
            stored_case = json.loads((stored_root / "library-qualification.json").read_text())["Cases"][0]
            record, context, table, digest, profile = self.prepare(window, case_id, self.tmp / f"s4-{kind}-none",
                                                                   extra=() if kind == "fab" else ("--controls", ""))
            self.assertIn("--reference none", record["Reference"]["Rule"])
            self.assertEqual(record["Kind"], "fabricated" if kind == "fab" else "thin")
            self.assertEqual(record["Orders"]["Controls"], ["p3", "p5"] if kind == "fab" else [])
            gate_record = qualify_library.analyze_case(record, context, stored_root / case_id / "results", gates=table, gates_digest=digest,
                                                       profile=profile)
            self.assertEqual((gate_record["Verdict"], stored_case["Qualification"]["Verdict"]), (expected, expected), kind)
            self.assertEqual(set(gate_record["Gates"]), {"PSequenceControls"})
            self.assertEqual(record["Status"], "pending-qualification" if kind == "fab" else "failed")
            self.assertIsNone(record["Cost"]["ReferenceNodeHours"])
            self.assertEqual(record["Qualification"]["Reason"], stored_case["Qualification"]["Reason"])

    def test_jobs_run_concurrently_up_to_max_jobs(self):
        """Two planned coupons (the S4 fab and its thin twin, one job each), --max-jobs 2, a fake remote
        replaying the stored results: both jobs are submitted before either is polled done, each coupon is
        fetched / verified / analyzed when its job leaves the queue, the fab reads PendingQualification and
        the thin Failed as recorded, the thin response is attached to the fab model, and --resume re-derives
        the plans byte-identical and adopts the recorded submissions (no upload, no qsub)."""
        window, spec = "sct002-S4", BATCH1_WINDOWS["sct002-S4"]
        fab, thin = spec["Case"], f"{spec['Case']}-thin"
        events = []
        self.run_root = self.tmp / "s4-concurrent"
        results_of = {fab: batch1_stored_root(window, "fab") / fab / "results" / "main",
                      thin: batch1_stored_root(window, "thin") / thin / "results" / "main"}
        self.replay_remote(results_of, events)
        saved_prepare = qualify_library.prepare_case

        def prepare_with_the_recorded_controls(case_record, **kwargs):
            # The stored thin twin ran --controls "" (decision 290); the one shared Namespace is the test's.
            kwargs["args"].controls = [] if case_record["Case"] == thin else [3, 5]
            return saved_prepare(case_record, **kwargs)
        qualify_library.prepare_case = prepare_with_the_recorded_controls
        self.addCleanup(setattr, qualify_library, "prepare_case", saved_prepare)
        args = self.args(window, self.run_root, fab, thin)
        record = qualify_library.run_qualify(args, log=lambda message: None)
        first_events = list(events)
        events.clear()
        args.resume = True
        resumed = qualify_library.run_qualify(args, log=lambda message: None)
        kinds = [kind for kind, _ in first_events]
        self.assertEqual(kinds[:4], ["upload", "submit", "upload", "submit"], first_events)
        self.assertLess(kinds.index("submit", kinds.index("submit") + 1), kinds.index("fetch"), "both jobs queued before any fetch")
        self.assertEqual([case for kind, case in first_events if kind == "fetch"], [fab, thin])
        totals = record["Library"]
        self.assertEqual((totals["JobsSubmitted"], totals["MaxJobs"], totals["CouponsPending"], totals["CouponsFailed"]), (2, 2, 1, 1))
        self.assertEqual(set(totals["JobWallSeconds"]), {fab, thin})
        by_case = {case["Case"]: case for case in record["Cases"]}
        stored = {case_id: json.loads((batch1_stored_root(window, kind) / "library-qualification.json").read_text())["Cases"][0]
                  for case_id, kind in ((fab, "fab"), (thin, "thin"))}
        for case_id, case in by_case.items():
            self.assertEqual(case["Status"], stored[case_id]["Status"], case.get("StoppedBy"))
            self.assertEqual(case["Qualification"]["Verdict"], stored[case_id]["Qualification"]["Verdict"])
            self.assertEqual(case["Monitor"]["LastJobState"], "F")
            self.assertTrue(all(entry["OK"] for entry in case["ResultDigests"].values()))
            self.assertEqual((case["Split"]["N"], case["JobPolicy"]["Mode"]), (1, "frugal"))
            self.assertEqual([job["Kind"] for job in case["Jobs"]], ["single"])
            self.assertAlmostEqual(case["Cost"]["JobTotalSeconds"], stored[case_id]["Cost"]["JobTotalSeconds"])
            self.assertEqual(case["Mesh"]["SHA256"], stored[case_id]["Mesh"]["SHA256"])
        self.assertEqual(by_case[fab]["Controls"]["Indices"], spec["Controls"])
        library = json.loads((self.run_root / "process-library.json").read_text())
        self.assertEqual(len(library["Models"]), 1)
        model = library["Models"][0]
        self.assertEqual(model["QualificationStatus"], "PendingQualification")
        self.assertTrue(model["ThinMatrix"] and model["FabricatedMatrix"])
        self.assertEqual(model["CouponMesh"]["SHA256"], stored[fab]["Mesh"]["SHA256"])
        resumed_kinds = [kind for kind, _ in events]
        self.assertNotIn("upload", resumed_kinds)
        self.assertNotIn("submit", resumed_kinds)
        self.assertEqual(resumed_kinds.count("fetch"), 2)
        for case in resumed["Cases"]:
            self.assertTrue(case["Monitor"]["Resumed"])
            self.assertEqual(case["Submission"], json.loads((self.run_root / case["Case"] / "submission.json").read_text()))

    def test_split_coupon_runs_worker_jobs_then_the_reducer_job(self):
        """The S3p first pass replayed: --job-policy frugal splits the 258 sources into 2 worker jobs ([94, 164]
        under the cost model of record; the first pass split [101, 157] under cost model V2) + the reducer job; both workers are submitted at once, each completed from its own stored
        status.json, the archive union counted before the reducer job, the coupon fetched / verified / analyzed
        after the reducer job; the first-pass controls [1, 2, 3, 8, 10, 31, 32, 108] and the stored plans' trace
        pins are reproduced; the gate reads the stored matrices under the gate table of record (the decision-474
        amplitude floor: a Type with no judged control is Unjudged, the lift of the first pass's Failed verdict is
        the decision-479 re-qualification's business)."""
        window, spec = "sct002-S3p", BATCH1_WINDOWS["sct002-S3p"]
        case_id = spec["Case"]
        stored_root = batch1_stored_root(window)
        stored_case = json.loads((stored_root / "library-qualification.json").read_text())["Cases"][0]
        events = []
        self.run_root = self.tmp / "s3p-split"
        self.replay_remote({case_id: stored_root / case_id / "results" / "main"}, events, archive_count=spec["Sources"] * 192)
        args = self.args(window, self.run_root, case_id)
        record = qualify_library.run_qualify(args, log=lambda message: None)
        case = record["Cases"][0]
        self.assertIsNone(case.get("StoppedBy"), case.get("StoppedBy"))
        # The split of record under the cost model of record (V3, decision 457 (2)): [94, 164] = the
        # decision-479 re-qualification's own record of the stored plan; the first pass ran cost model V2
        # (f695da1b...) and split [101, 157] - the stored per-block statuses merge to the same stage totals.
        requal_case = json.loads((BATCH1 / "requalify" / window / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual((case["Split"]["N"], case["JobPolicy"]["Mode"], case["Split"]["Blocks"]), (2, "frugal", requal_case["Split"]["Blocks"]))
        self.assertEqual((case["Split"]["Blocks"], stored_case["Split"]["Blocks"]), ([94, 164], [101, 157]))
        self.assertEqual(case["Controls"]["Indices"], spec["Controls"])
        self.assertEqual(case["Plan"]["StageNames"], stored_case["Plan"]["StageNames"])
        stored_plans = [json.loads(path.read_text()) for path in sorted((stored_root / case_id / "main" / "jobs").glob("*/plan.json"))]
        self.assert_plan_reproduces_stored(case, stored_case, stored_plans, stage_configs={
            f"{case_id}-p4": ("reducer.json",), f"{case_id}-p5-control": ("worker.json", "reducer.json"),
            f"{case_id}-p3-control": ("worker.json", "reducer.json"), f"{case_id}-p4-local-edge": ("config.json",)})
        blocks = [json.loads((Path(case["Root"]) / "main" / f"{case_id}-p4" / f"worker-block{k}.json").read_text()) for k in (1, 2)]
        self.assertEqual([entry["Index"] for block in blocks for entry in block["Boundaries"]["PrescribedPotential"]], list(range(1, spec["Sources"] + 1)))
        kinds = [kind for kind, _ in events]
        self.assertEqual(events[:3], [("upload", case_id), ("submit", "worker-1"), ("submit", "worker-2")])
        self.assertLess(kinds.index("poll"), kinds.index("archive-count"))
        self.assertEqual([name for kind, name in events if kind == "status"], ["worker-1", "worker-2"])
        self.assertEqual(events.index(("archive-count", f"{case_id}-p4")) + 1, events.index(("submit", "reducer")))
        self.assertLess(events.index(("submit", "reducer")), events.index(("fetch", case_id)))
        self.assertEqual(kinds.count("fetch"), 1)
        self.assertEqual([job["Name"] for job in case["Jobs"]], ["worker-1", "worker-2", "reducer"])
        self.assertEqual(case["Jobs"][2]["Requires"], ["worker-1", "worker-2"])
        self.assertEqual(case["ArchiveUnion"]["Stages"][f"{case_id}-p4"]["Expected"], spec["Sources"] * 192)
        for job, stored_job in zip(case["Jobs"], stored_case["Jobs"]):
            self.assertEqual(job["Monitor"]["LastJobState"], "F")
            self.assertEqual(job["StageNames"], stored_job["StageNames"])
        assert_numbers_reproduced(self, stored_case["Cost"]["MainStage"], case["Cost"]["MainStage"], "Cost.MainStage")
        self.assertAlmostEqual(case["Cost"]["JobTotalSeconds"], stored_case["Cost"]["JobTotalSeconds"])
        self.assertEqual(record["Library"]["JobsSubmitted"], 3)
        # The stored first pass (gate table V5, before the amplitude floor) read Failed on below-floor controls;
        # the gate of record judges only above-floor Types: the verdict is Failed or PendingQualification with
        # the below-floor Types recorded, never Passed / Qualified (decisions 474 / 477 (1)).
        self.assertIn(case["Qualification"]["Verdict"], (gates.VERDICT_FAILED, gates.VERDICT_PENDING))
        self.assertEqual(stored_case["Qualification"]["Verdict"], gates.VERDICT_FAILED)
        gate_record = json.loads((self.run_root / case_id / "qualification.json").read_text())
        stored_gate = json.loads((stored_root / case_id / "qualification.json").read_text())
        for control, observables in stored_gate["Gates"]["PSequenceControls"]["Controls"].items():
            for name, stored_step in observables.items():
                step = gate_record["Gates"]["PSequenceControls"]["Controls"][control][name]
                assert_numbers_reproduced(self, {key: stored_step[key] for key in ("StepToHigherOrder", "StepFromLowerOrder", "Bound")},
                                          step, f"control {control} {name}")
        library = json.loads((self.run_root / "process-library.json").read_text())
        self.assertFalse(library["Models"][0]["LibraryQualified"])
        self.assertEqual(library["Models"][0]["QualificationStatus"], "Failed" if case["Qualification"]["Verdict"] == gates.VERDICT_FAILED
                         else "PendingQualification")
        # --resume adopts the three recorded submissions (plans byte-identical).
        del events[:]
        args.resume = True
        resumed = qualify_library.run_qualify(args, log=lambda message: None)
        resumed_case = resumed["Cases"][0]
        self.assertIsNone(resumed_case.get("StoppedBy"), resumed_case.get("StoppedBy"))
        self.assertEqual([kind for kind, _ in events if kind in ("upload", "submit")], [])
        self.assertEqual([job["Submission"]["Job"] for job in resumed_case["Jobs"]], ["worker-1.fake", "worker-2.fake", "reducer.fake"])

    def merged_stored_root(self, windows):
        """ONE stored root holding the first-pass runs of several windows (their records concatenated, the
        case directories symlinked read-only): the decision-495 (1) two-cases fixture."""
        merged = self.tmp / ("stored-" + "-".join(windows))
        if merged.exists():
            return merged
        merged.mkdir()
        records = [json.loads((batch1_stored_root(window) / "library-qualification.json").read_text()) for window in windows]
        for window in windows:
            case_id = BATCH1_WINDOWS[window]["Case"]
            os.symlink(batch1_stored_root(window) / case_id, merged / case_id)
        (merged / "library-qualification.json").write_text(json.dumps({**records[0], "Cases": [case for record in records for case in record["Cases"]]}))
        return merged

    def controls_only_args(self, window, root, case_id, reuse_root, build_record=None):
        args = self.args(window, root, case_id, extra=("--controls", "p3,p5", "--controls-only", "--reuse-main", str(reuse_root)))
        if build_record is not None:
            args.build_record = build_record
        return args

    def test_controls_only_two_cases_take_their_own_reducers(self):
        """Decision 495 (1): the S3p and S2p first-pass runs in ONE stored root, ONE controls-only run of both
        (one build record holding both cases): each coupon's PriorMainStage is its OWN stored reducer, the
        control sets are the decision-479 sets of record, and the shared args Namespace is never written."""
        windows = ("sct002-S3p", "sct002-S2p")
        stored = self.merged_stored_root(windows)
        # One build record of both cases (the two windows' records, the same manifest recipe).
        builds = [json.loads(self.build_records[window].read_text()) for window in windows]
        both = self.tmp / "both-build.json"
        both.write_text(json.dumps({**builds[0], "Cases": [case for build in builds for case in build["Cases"] if case["Case"] in
                                                           {BATCH1_WINDOWS[w]["Case"] for w in windows}]}))
        root = self.tmp / "controls-only-two"
        root.mkdir()
        cases = [BATCH1_WINDOWS[window]["Case"] for window in windows]
        parser = argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        args = parser.parse_args(["--build-record", str(both), "--reference", "none", "--frozen-binary-sha256", BATCH1_BINARY_SHA256,
                                  "--controls", "p3,p5", "--controls-only", "--reuse-main", str(stored), "--root", str(root), "--dry-run",
                                  "--job-policy", "frugal", *[token for case in cases for token in ("--case", case)]])
        manifests = {window: json.loads(Path(self.build_records[window].parent / "manifest-main.json").read_text())
                     | {"Path": str(self.build_records[window].parent / "manifest-main.json")} for window in windows}
        profile = json.loads(Path(args.cluster_profile).read_text())
        cost_model = estimate_stages.load_cost_model(args.cost_model)
        table, digest = gates.load_gates(args.gates)
        before = dict(vars(args))
        priors, controls = {}, {}
        for window in windows:
            case_id = BATCH1_WINDOWS[window]["Case"]
            # Each case's sources live under its own window's registration: its own manifest copy.
            record, _ = qualify_library.prepare_case(self.build_case(window, case_id), manifest_path=Path(manifests[window]["Path"]),
                                                     manifest=manifests[window], args=args, root=root, remote=None, profile=profile,
                                                     cost_model=cost_model, gates=table, gates_digest=digest)
            self.assertEqual(dict(vars(args)), before)
            priors[case_id] = record["ControlAmplitudes"]["PriorMainStage"]
            controls[case_id] = record["Controls"]["Indices"]
            self.assertEqual(record["ReusedMain"]["Reducer"]["Directory"], priors[case_id])
            self.assertEqual(record["ReusedMain"]["Case"], case_id)
            self.assertEqual(record["ReusedMain"]["StoredControls"], BATCH1_WINDOWS[window]["Controls"])
            self.assertEqual(record["ReusedMain"]["StoredVerdict"], gates.VERDICT_FAILED)
            self.assertIn("decision 474 (A)", record["Controls"]["Rule"])
        self.assertEqual(priors, {BATCH1_WINDOWS[window]["Case"]: str(stored / BATCH1_WINDOWS[window]["Case"] / "results" / "main"
                                                                        / f"{BATCH1_WINDOWS[window]['Case']}-p4" / "reducer") for window in windows})
        self.assertEqual(controls, {BATCH1_WINDOWS[window]["Case"]: BATCH1_WINDOWS[window]["RequalifiedControls"] for window in windows})

    def test_controls_only_reuses_the_stored_main_stage_and_plans_the_control_job(self):
        """Decisions 474 (A) / 479 / 485 (c) / 500 (3): the controls-only re-qualification of the two 8-edge
        coupons S3p and S2p replayed END TO END from the stored inputs (no solve) - --controls-only --reuse-main
        on the stored first-pass root reproduces the decision-479 control sets [1, 2, 3, 15, 18, 32, 48, 108] /
        [1, 2, 3, 7, 20, 32, 33, 93] and the job plan of record (stages, caps, ranks, the trace pins and the
        configs the job ran); the controls-only job's stored results, fetched through the fake remote and the
        stored reducer spliced in, give the re-qualification verdict of record (PendingQualification: every judged
        control passed) with every p-sequence step to the stored digit; and the (F) evaluate on this run's main
        reducer with the stored dense traces reproduces the -requal (F) record to the printed digit."""
        controls_of_record = json.loads(BATCH1_REQUALIFICATION_CONTROLS.read_text())
        for window in ("sct002-S3p", "sct002-S2p"):
            with self.subTest(window=window):
                self.replay_requalification(window, controls_of_record[f"fab-{BATCH1_WINDOWS[window]['Case']}"])

    def replay_requalification(self, window, controls_record):
        spec = BATCH1_WINDOWS[window]
        case_id = spec["Case"]
        stored_root = batch1_stored_root(window)
        requalify_root = BATCH1 / "requalify" / window
        requal_case = json.loads((requalify_root / "library-qualification.json").read_text())["Cases"][0]
        self.assertEqual(controls_record["AmplitudeInformedControls"], spec["RequalifiedControls"])
        self.assertEqual(controls_record["RecordedControls"], spec["Controls"])
        self.assertEqual(requal_case["Controls"]["Indices"], spec["RequalifiedControls"])
        events = []
        self.run_root = self.tmp / f"{window}-requalify"
        self.replay_remote({case_id: requalify_root / case_id / "results" / "main"}, events)
        args = self.controls_only_args(window, self.run_root, case_id, stored_root)
        record = qualify_library.run_qualify(args, log=lambda message: None)
        case = record["Cases"][0]
        self.assertIsNone(case.get("StoppedBy"), case.get("StoppedBy"))
        # (1) the decision-479 control set and the amplitude-informed rule from the stored reducer.
        self.assertEqual(case["Controls"]["Indices"], spec["RequalifiedControls"])
        self.assertEqual(case["ControlAmplitudes"]["PriorMainStage"], str(stored_root / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"))
        assert_numbers_reproduced(self, {key: requal_case["ControlAmplitudes"][key] for key in ("FloorRatio", "Observables", "BelowFloorClasses")},
                                  case["ControlAmplitudes"], "ControlAmplitudes")
        reused = case["ReusedMain"]
        self.assertEqual((reused["Case"], reused["MainPrefix"], reused["StoredControls"], reused["StoredVerdict"]),
                         (case_id, f"{case_id}-p4", spec["Controls"], gates.VERDICT_FAILED))
        self.assertEqual(reused["Traces"], {"Count": spec["Sources"], "PinnedByStoredPlans": spec["Sources"], "Identical": True})
        self.assertEqual((reused["Reducer"]["Job"], reused["Reducer"]["Nodes"]), ("reducer", 1))
        # (2) the job plan of record: one controls-only job of the fixed stages; stages / caps / ranks identical
        # to the stored requalification plan, every pin in common identical (the stored plan - the lane script's -
        # pinned the p4 reducer.json it never ran instead of the 5 configs the job runs; the instance is the
        # driver's choice by the fixed stages' estimated peak: both recorded in the qualify-tools lane REPORT).
        self.assertEqual([(job["Name"], job["Kind"], job["Sources"]) for job in case["Jobs"]], [("controls-only", "controls-only", [])])
        plan = json.loads(Path(case["Jobs"][0]["Plan"]).read_text())
        stored_plan = json.loads((requalify_root / case_id / "main" / "jobs" / "controls-only" / "plan.json").read_text())
        self.assertEqual(plan["JobKind"], stored_plan["JobKind"])

        def stage_view(stage):
            environment = {key: (Path(value).name if key == "PALACE_RESPONSE_ARCHIVE_DIR" else value) for key, value in stage["Environment"].items()}
            return (stage["Name"], Path(stage["Config"]).name, environment, stage["Requires"], stage["CapSeconds"], stage["MinimumSeconds"])
        self.assertEqual([stage_view(s) for s in plan["Stages"]], [stage_view(s) for s in stored_plan["Stages"]])
        self.assertEqual((plan["Ranks"], plan["BinarySHA256"], plan["MeshSHA256"]), (stored_plan["Ranks"], stored_plan["BinarySHA256"], stored_plan["MeshSHA256"]))
        pins = {key.split(f"/{case_id}/", 1)[1]: value for key, value in plan["PinnedSHA256"].items()}
        stored_pins = {key.split(f"/{case_id}/", 1)[1]: value for key, value in stored_plan["PinnedSHA256"].items()}
        common = set(pins) & set(stored_pins)
        self.assertEqual(len(common), spec["Sources"] + 1)       # every trace + the mesh
        self.assertEqual({key: pins[key] for key in common}, {key: stored_pins[key] for key in common})
        self.assertEqual(set(pins) - set(stored_pins), {f"main/{case_id}-p3-control/worker.json", f"main/{case_id}-p3-control/reducer.json",
                                                         f"main/{case_id}-p5-control/worker.json", f"main/{case_id}-p5-control/reducer.json",
                                                         f"main/{case_id}-p4-local-edge/config.json"})
        self.assertEqual(set(stored_pins) - set(pins), {f"main/{case_id}-p4/reducer.json"})
        for stage, names in {f"{case_id}-p5-control": ("worker.json", "reducer.json"), f"{case_id}-p3-control": ("worker.json", "reducer.json"),
                             f"{case_id}-p4-local-edge": ("config.json",)}.items():
            for name in names:
                generated = strip_paths(json.loads((Path(case["Root"]) / "main" / stage / name).read_text()))
                stored = strip_paths(json.loads((requalify_root / case_id / "main" / stage / name).read_text()))
                self.assertEqual(generated, stored, f"{stage}/{name}")
        # (3) the stored controls-only results through the real fetch / splice / analysis path: the verdict of record.
        self.assertEqual([kind for kind, _ in events if kind in ("upload", "submit", "fetch")], ["upload", "submit", "fetch"])
        self.assertEqual(case["ReusedMain"]["Splice"]["SHA256"], reused["Reducer"]["SHA256"])
        self.assertEqual(case["Status"], requal_case["Status"])
        self.assertEqual(case["Qualification"]["Verdict"], requal_case["Qualification"]["Verdict"])
        self.assertEqual(case["Qualification"]["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(case["Qualification"]["GatesPassed"], {"PSequenceControls": True})
        gate_record = json.loads((self.run_root / case_id / "qualification.json").read_text())
        stored_gate = json.loads((requalify_root / case_id / "qualification.json").read_text())
        assert_numbers_reproduced(self, stored_gate["Gates"]["PSequenceControls"]["Controls"], gate_record["Gates"]["PSequenceControls"]["Controls"],
                                  "PSequenceControls.Controls")
        self.assertEqual(gate_record["Verdict"], stored_gate["Verdict"])
        self.assertEqual(case["Cost"]["ReusedMainStage"]["Cost"], reused["StoredMainStageCost"])
        self.assertLess(case["Cost"]["JobNodeHours"], reused["StoredJobNodeHours"])
        library = json.loads((self.run_root / "process-library.json").read_text())
        model = library["Models"][0]
        self.assertEqual((model["Name"], model["QualificationStatus"]), (spec["Model"], "PendingQualification"))
        self.assertEqual(model["ControlsOnlyRequalification"]["Controls"], spec["RequalifiedControls"])
        # Decision 513 (1): the record of record and the library model carry the LIFTED usage state, not the restriction.
        for usage in (case["ReusedMain"]["Usage"], model["ControlsOnlyRequalification"]["Usage"]):
            self.assertEqual(usage, job_split.CONTROLS_ONLY_USAGE_RECORD)
            self.assertEqual(usage["LiftedBy"], "decision 513 (1)")
            self.assertIn("LIFTED by decision 513", usage["Rule"])
        self.assertNotIn("UsageRestriction", model["ControlsOnlyRequalification"])
        # (4) the (F) evaluate (decision 477 / 485 (a)) on THIS run's main reducer (the spliced stored one) with the
        # stored dense traces, thin matrices, gate solves and the requalification library: the -requal (F) record.
        f_root = BATCH1 / "f-qualification" / case_id
        stored_f = json.loads((f_root / "spatial-qualification-requal.json").read_text())
        output = self.run_root / "spatial-qualification-requal.json"
        # The stored dense-traces manifest and the twin configs it names carry the lane's cluster paths:
        # relocated copies under this run's root (the stored twin outputs are read where they are).
        manifest = json.loads(relocated((f_root / "dense" / "dense-traces.json").read_text()))
        for name, config_path in manifest["Configs"].items():
            local_config = self.run_root / "dense" / name / "config.json"
            local_config.parent.mkdir(parents=True, exist_ok=True)
            local_config.write_text(relocated(Path(config_path).read_text()))
            manifest["Configs"][name] = str(local_config)
        dense_traces = self.run_root / "dense-traces.json"
        dense_traces.write_text(json.dumps(manifest, indent=2))
        command = [sys.executable, str(HERE / "coupon_library.py"), "spatial-qualify", "evaluate", "--dense-traces", str(dense_traces),
                   "--fabricated-matrices", str(self.run_root / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"),
                   "--thin-matrices", str(batch1_stored_root(window, "thin") / f"{case_id}-thin" / "results" / "main" / f"{case_id}-thin-p4" / "reducer"),
                   *[token for order in (3, 4, 5) for token in ("--gate", str(BATCH1 / "preflight" / window / f"c0-p{order}" / "postpro" / "palace.json"))],
                   "--library", str(BATCH1 / "library" / f"b1-requal-interim-{window.split('-')[1]}" / "process-library.json"),
                   "--output", str(output), "--dry-run"]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        evaluated = json.loads(output.read_text())
        self.assertEqual(evaluated["Status"], stored_f["Status"])
        self.assertEqual(evaluated["Status"], "Qualified")
        assert_numbers_reproduced(self, stored_f, evaluated, rel=0.0, ignore=("Gate.Solves[0].Source", "Gate.Solves[1].Source", "Gate.Solves[2].Source",
                                                                              "PreviousStatus", "Rule", "UnjudgedTypesRule"))
        self.assertEqual(evaluated["ControlVerdict"], gates.VERDICT_PENDING)


if __name__ == "__main__":
    unittest.main()
