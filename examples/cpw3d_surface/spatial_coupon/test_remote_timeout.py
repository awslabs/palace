# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""OI-REMOTE-SUBMIT-SSH-NO-TIMEOUT (decision 1246; fix round P2, decision 1260 Q-E): every
remote round trip of qualify/remote.py carries a timeout (overridable by the environment); a
READ-ONLY round trip that times out is retried once; qsub is never retried and fails closed
with the reconciliation named; the default command lines are unchanged."""
import importlib
import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import remote  # noqa: E402


class FakeRun:
    """A subprocess.run stand-in: the first `timeouts` calls raise TimeoutExpired, the next
    returns `result`; every call is recorded with its keyword arguments."""

    def __init__(self, timeouts, result=None):
        self.timeouts = timeouts
        self.result = result
        self.calls = []

    def __call__(self, command, **kwargs):
        self.calls.append((command, kwargs))
        if len(self.calls) <= self.timeouts:
            raise subprocess.TimeoutExpired(command, kwargs.get("timeout"))
        return self.result


def completed(stdout, returncode=0):
    return subprocess.CompletedProcess(["ssh"], returncode, stdout=stdout, stderr="")


class RemoteTimeoutTest(unittest.TestCase):
    def test_default_ssh_call_carries_the_timeout_and_nothing_else_changes(self):
        fake = FakeRun(0, completed("ok\n"))
        with mock.patch.object(remote.subprocess, "run", fake):
            result = remote.ssh("host", "echo ok")
        self.assertEqual(result.stdout, "ok\n")
        command, kwargs = fake.calls[0]
        self.assertEqual(command, ["ssh", "host", "echo ok"])
        self.assertEqual(kwargs, {"text": True, "capture_output": True, "check": True,
                                  "timeout": remote.SSH_TIMEOUT_SECONDS})
        self.assertEqual(remote.SSH_TIMEOUT_SECONDS, 300.0)
        self.assertEqual(remote.RSYNC_TIMEOUT_SECONDS, 1800.0)

    def test_environment_overrides_the_timeouts(self):
        with mock.patch.dict(os.environ, {"COUPON_LIBRARY_SSH_TIMEOUT_SECONDS": "7.5",
                                          "COUPON_LIBRARY_RSYNC_TIMEOUT_SECONDS": "42"}):
            overridden = importlib.reload(remote)
            self.assertEqual(overridden.SSH_TIMEOUT_SECONDS, 7.5)
            self.assertEqual(overridden.RSYNC_TIMEOUT_SECONDS, 42.0)
        restored = importlib.reload(remote)
        self.assertEqual(restored.SSH_TIMEOUT_SECONDS, 300.0)

    def test_read_only_round_trip_is_retried_once_then_fails_closed(self):
        # One timeout, then an answer: two calls, the answer returned.
        fake = FakeRun(1, completed("abc  /x\n"))
        with mock.patch.object(remote.subprocess, "run", fake):
            digests = remote.remote_sha256("host", ["/x"])
        self.assertEqual(len(fake.calls), 2)
        self.assertEqual(digests, {"/x": "abc"})
        # Two timeouts: the retry budget (READ_ONLY_RETRIES = 1) is spent; check=True raises.
        fake = FakeRun(2)
        with mock.patch.object(remote.subprocess, "run", fake):
            with self.assertRaises(remote.RemoteTimeout) as raised:
                remote.remote_sha256("host", ["/x"])
        self.assertEqual(len(fake.calls), 2)
        self.assertIn("did not answer within", str(raised.exception))

    def test_read_only_probes_report_unreachable_on_a_timeout(self):
        fake = FakeRun(2)
        with mock.patch.object(remote.subprocess, "run", fake):
            record = remote.poll("host", "/opt/pbs/bin", "1.host", "/case/status.json")
        self.assertFalse(record["Reachable"])
        self.assertIsNone(record["JobState"])
        self.assertEqual(record["SSHReturnCode"], 124)
        self.assertEqual(len(fake.calls), 2)
        fake = FakeRun(2)
        with mock.patch.object(remote.subprocess, "run", fake):
            reachable, content = remote.read_json_reachable("host", "/case/status.json")
        self.assertFalse(reachable)
        self.assertIsNone(content)
        fake = FakeRun(2)
        with mock.patch.object(remote.subprocess, "run", fake):
            with self.assertRaises(RuntimeError) as raised:
                remote.path_exists("host", "/case")
        self.assertIn("returned no answer", str(raised.exception))
        fake = FakeRun(2)
        with mock.patch.object(remote.subprocess, "run", fake):
            self.assertEqual(remote.count_archive_potentials("host", "/case/archive"), 0)

    def test_qsub_is_never_retried_and_fails_closed_with_the_reconciliation(self):
        qstat = completed('{"Jobs": {}}')
        fake = FakeRun(0, qstat)
        calls = []

        def run(command, **kwargs):
            calls.append(command)
            if "qsub" in command[2]:
                raise subprocess.TimeoutExpired(command, kwargs.get("timeout"))
            return qstat

        with mock.patch.object(remote.subprocess, "run", run):
            with self.assertRaises(remote.RemoteTimeout) as raised:
                remote.submit("host", "/opt/pbs/bin", "job.pbs", "/case", job_cap=40, user="u")
        self.assertEqual(sum("qsub" in c[2] for c in calls), 1)  # one qsub attempt only
        message = str(raised.exception)
        self.assertIn("receipt is UNKNOWN", message)
        self.assertIn("reconcile with qstat before any retry", message)
        self.assertIn("decision 1246", message)
        del fake

    def test_rsync_transfer_times_out_without_a_retry(self):
        fake = FakeRun(1)
        with mock.patch.object(remote.subprocess, "run", fake):
            with self.assertRaises(remote.RemoteTimeout) as raised:
                remote.fetch("host", "/remote/case", "/tmp/local")
        self.assertEqual(len(fake.calls), 1)
        command, kwargs = fake.calls[0]
        self.assertEqual(command, remote.fetch_command("host", "/remote/case", "/tmp/local"))
        self.assertEqual(kwargs, {"check": True, "timeout": remote.RSYNC_TIMEOUT_SECONDS})
        self.assertIn("re-runnable", str(raised.exception))

    def test_state_changing_round_trips_are_not_retried(self):
        # mkdir -p before an upload and rm -rf of the archives: one attempt each.
        fake = FakeRun(1)
        with mock.patch.object(remote.subprocess, "run", fake):
            with self.assertRaises(remote.RemoteTimeout):
                remote.upload("host", "/tmp/local", "/remote/case")
        self.assertEqual(len(fake.calls), 1)
        self.assertEqual(fake.calls[0][1]["timeout"], remote.SSH_TIMEOUT_SECONDS)


if __name__ == "__main__":
    unittest.main()
