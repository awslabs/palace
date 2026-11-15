# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Decision 1799 (l): the --backend switch of `coupon-library qualify` (qualify_library.py).
Checked here: --backend pbs is the default and refuses a profile with a Paladin block;
--backend paladin refuses a profile without one, one with a Launcher block too, and one with an
unknown key, before anything is planned; it prints ONE route line (selector, pin record);
executables_of names the pinned image's executable; and the job wrappers submit_paladin_job /
poll_paladin_job drive a plan through the FAKE paladin CLI (a paladin STOP stops the coupon).
The fetch / evaluate / merge of a paladin coupon is the PBS code path with three branches
(finish_case): it is exercised end to end by the live rehearsal only (HANDBOOK s18)."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import paladin_backend  # noqa: E402
import qualify_library  # noqa: E402
import test_paladin_backend  # noqa: E402  (the fake block and plan)

PROFILE = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())


class BackendSwitch(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.block = test_paladin_backend.fake_block(self.root)
        manifest = self.root / "manifest.json"
        manifest.write_text("{}")
        self.build = self.root / "build.json"
        self.build.write_text(json.dumps({"Library": {"Manifest": {"Path": str(manifest),
                                                                   "SHA256": hashlib.sha256(b"{}").hexdigest()},
                                                      "Commit": "c0ffee"}, "Cases": []}))

    def tearDown(self):
        self.tmp.cleanup()

    def arguments(self, profile, *extra):
        path = self.root / "profile.json"
        path.write_text(json.dumps(profile))
        parser = qualify_library.argparse.ArgumentParser()
        qualify_library.add_arguments(parser)
        return parser.parse_args(["--build-record", str(self.build), "--reference", "none", "--dry-run",
                                  "--cluster-profile", str(path), "--root", str(self.root / "run"), *extra])

    def run_qualify(self, args):
        lines = []
        with contextlib.redirect_stdout(io.StringIO()):
            qualify_library.run_qualify(args, log=lines.append)
        return lines

    def test_the_default_is_pbs_and_prints_no_route(self):
        args = self.arguments(PROFILE)
        self.assertEqual(args.backend, "pbs")
        self.assertFalse([line for line in self.run_qualify(args) if line.startswith("route ")])

    def test_paladin_prints_one_route_line(self):
        lines = self.run_qualify(self.arguments({**PROFILE, "Paladin": self.block}, "--backend", "paladin"))
        routes = [line for line in lines if line.startswith("route ")]
        self.assertEqual(len(routes), 1)
        self.assertIn(f"route paladin: selector {self.block['Selector']}, pin record {self.block['PinRecord']}", routes[0])

    def test_refusals_before_anything_is_planned(self):
        for profile, backend in (({**PROFILE, "Paladin": self.block}, "pbs"),
                                 (PROFILE, "paladin"),
                                 ({**PROFILE, "Paladin": self.block, "Launcher": {"Kind": "apptainer"}}, "paladin"),
                                 ({**PROFILE, "Paladin": {**self.block, "Unknown": 1}}, "paladin")):
            with self.subTest(backend=backend, keys=sorted(set(profile) - set(PROFILE))):
                with self.assertRaises(ValueError):
                    self.run_qualify(self.arguments(profile, "--backend", backend))
                self.assertFalse((self.root / "run").exists())

    def test_a_paladin_run_holds_the_root_lock_and_a_second_driver_is_refused(self):
        # Decision 1820 (2): one driver per root (paladin_backend.LOCK_RULE); a dry run takes no lock.
        def arguments():
            args = self.arguments({**PROFILE, "Paladin": self.block}, "--backend", "paladin", "--remote", "host:/data/r")
            args.dry_run = False
            return args

        lock = self.root / "run" / "chain-driver.lock"
        held = paladin_backend.DriverLock(lock, log=lambda line: None).acquire()
        try:
            with self.assertRaises(paladin_backend.PaladinError) as stop:
                self.run_qualify(arguments())
            self.assertIn("held by a live driver", str(stop.exception))
        finally:
            held.release()
        self.assertFalse((self.root / "run" / qualify_library.GATES_COPY).exists())  # refused before anything is written
        observed = []
        real = paladin_backend.DriverLock.acquire

        def acquire(lock_object):
            observed.append(lock_object.path)
            return real(lock_object)

        with mock.patch.object(paladin_backend.DriverLock, "acquire", acquire):
            try:
                self.run_qualify(arguments())
            except Exception:  # noqa: BLE001  (no cluster here: only the lock's life is checked)
                pass
        self.assertEqual(observed, [str(lock)])
        self.assertEqual(lock.read_text(), "")  # released at the end, whatever the end
        self.run_qualify(self.arguments({**PROFILE, "Paladin": self.block}, "--backend", "paladin"))  # a dry run: no lock
        self.assertEqual(lock.read_text(), "")

    def test_executables_of_names_the_pinned_image(self):
        profile = {**PROFILE, "Paladin": self.block}
        args = types.SimpleNamespace(frozen_binary_remote_path=None, mpiexec_remote_path=None)
        self.assertEqual(qualify_library.executables_of(args, profile, "/r", self.block["ImageBinarySHA256"]),
                         (self.block["ImageBinary"], "paladin"))
        with self.assertRaises(ValueError):
            qualify_library.executables_of(args, profile, "/r", "0" * 64)
        with self.assertRaises(ValueError):
            qualify_library.executables_of(types.SimpleNamespace(frozen_binary_remote_path="/x", mpiexec_remote_path=None),
                                           profile, "/r", self.block["ImageBinarySHA256"])


class PaladinJobWrappers(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.block = test_paladin_backend.fake_block(self.root)
        self.plan_path, archive = test_paladin_backend.three_stage_plan(self.root)
        (self.root / "runs").mkdir()
        self.state = self.root / "fake-state.json"
        self.state.write_text(json.dumps({"Root": str(self.root / "runs"), "Plan": {}}))
        self.environment = mock.patch.dict(os.environ, {"FAKE_PALADIN_STATE": str(self.state),
                                                        "FAKE_PALADIN_LOG": str(self.root / "argv.jsonl"),
                                                        "FAKE_TOOL_LOG": str(self.root / "tools.jsonl")})
        self.environment.start()
        self.pause = mock.patch.object(paladin_backend, "pause")  # the waits between lookups (INTENT_RULE)
        self.pause.start()
        self.watchers = mock.patch.object(paladin_backend, "RunWatchers", test_paladin_backend.RecordingWatchers)
        self.watchers.start()
        self.cookie = mock.patch.object(paladin_backend, "MIDWAY_COOKIE",
                                        test_paladin_backend.write_cookie(self.root / "cookie", 86400))
        self.cookie.start()
        estimate = self.root / "stage-estimate.json"
        estimate.write_text(json.dumps({"Stages": {"k4": {"ArchiveGBEstimate": 0.001}}}))
        main = self.plan_path.parent
        self.local = self.root / "local-case" / "main"
        self.local.mkdir(parents=True)
        self.record = {"Case": "case", "Root": str(self.local.parent), "Remote": {"Case": str(main.parent)},
                       "Estimate": {"Path": str(estimate)}, "Upload": {"Test": "no upload"}}
        self.job = {"Name": "single", "Kind": "single", "Plan": str(self.plan_path), "RemoteDirectory": str(main),
                    "Directory": str(self.local), "SubmissionRecord": str(self.local.parent / "submission.json")}
        self.context = {"paladin": self.block, "layout": [{"Prefix": "p4", "Kind": "response", "EstimateKey": "k4"}],
                        "jobs": [self.job]}
        # The configurations are read from the local case root: here the same files.
        os.rmdir(self.local)
        os.symlink(main, self.local)

    def tearDown(self):
        self.cookie.stop()
        self.watchers.stop()
        self.pause.stop()
        self.environment.stop()
        self.tmp.cleanup()

    def test_submit_and_poll_to_completion(self):
        lines = []
        qualify_library.submit_paladin_job(self.record, self.context, self.job, remote=None, profile=PROFILE,
                                           log=lines.append)
        submission = json.loads(Path(self.job["SubmissionRecord"]).read_text())
        self.assertEqual((submission["Backend"], submission["Job"]), ("paladin", "paladin:case/single"))
        polls = 1
        while not qualify_library.poll_paladin_job(self.record, self.context, self.job, log=lines.append):
            polls += 1
            self.assertLess(polls, 20)
        status = json.loads((self.plan_path.parent / "status.json").read_text())
        self.assertEqual((status["State"], status["PaladinRunIDs"]), ("complete", list(self.job["Paladin"]["Runs"].values())))
        self.assertEqual(self.job["Monitor"]["LastJobState"], "complete")

    def test_a_paladin_stop_stops_the_coupon(self):
        self.state.write_text(json.dumps({"Root": str(self.root / "runs"), "Plan": {"chainport-w": "failed"}}))
        qualify_library.submit_paladin_job(self.record, self.context, self.job, remote=None, profile=PROFILE,
                                           log=lambda line: None)
        with self.assertRaises(qualify_library.CaseStop) as stop:
            for _ in range(10):
                qualify_library.poll_paladin_job(self.record, self.context, self.job, log=lambda line: None)
        self.assertEqual(stop.exception.record["Kind"], "Paladin")


if __name__ == "__main__":
    unittest.main()
