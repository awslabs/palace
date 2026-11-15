# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form X of `coupon-library qualify` (USER decision 1674 (2)): a cluster profile with a
Launcher block runs every stage inside a published Palace image through apptainer on ONE
node.  Checked here: the block's validation, the plan fields and their refusals, the job
script (no module line, every other line as the PBS form), the executables and the remote
image check of the driver, the dense-twin planner, run_stages' launch command and preflight
(a fake apptainer that records its argv and environment), and palace_in_image.py.  Without
the block the PBS form is unchanged (the byte-equality of the plans and job scripts of
record is proved by replay on the records)."""
import contextlib
import io
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import time
import unittest
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import build_plan  # noqa: E402
import dense_twin_plan  # noqa: E402
import palace_in_image  # noqa: E402
import qualify_library  # noqa: E402
import run_stages  # noqa: E402

DIGEST_IMAGE = "1" * 64
DIGEST_BINARY = "2" * 64


def launcher(**changes):
    block = {"Kind": "apptainer", "Apptainer": "/bin/apptainer", "Image": "/data/sif/palace.sif", "ImageSHA256": DIGEST_IMAGE,
             "Binds": ["/data/home/simlap"], "ImageMPIExec": "/opt/view/bin/mpirun",
             "ImageBinary": "/opt/view/bin/palace-arm64.bin", "ImageBinarySHA256": DIGEST_BINARY,
             "MPIExecArguments": ["--prefix", "/opt/view", "--map-by", "numa", "--bind-to", "core"],
             "Environment": {"GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"}}
    block.update(changes)
    return block


def profile_x(**changes):
    profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
    profile["Launcher"] = launcher(**changes)
    return profile


def render(profile, nodes=1):
    return build_plan.render_job_script(profile=profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run_stages.py",
                                        job_name="coupon-qualify-case", walltime_seconds=12600, instance_type="m8g.48xlarge",
                                        nodes=nodes)


class LauncherBlockTest(unittest.TestCase):
    def test_absent_block_is_the_pbs_form(self):
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        self.assertIsNone(build_plan.launcher_of(profile))
        self.assertEqual(build_plan.launcher_plan_fields(profile, binary="b", binary_sha256="d", mpiexec="m", nodes=4), {})

    def test_block_needs_exactly_its_keys_and_kind_apptainer(self):
        self.assertEqual(build_plan.launcher_of(profile_x())["Image"], "/data/sif/palace.sif")
        broken = profile_x()
        del broken["Launcher"]["Binds"]
        with self.assertRaisesRegex(ValueError, "must have the keys"):
            build_plan.launcher_of(broken)
        extra = profile_x()
        extra["Launcher"]["Modules"] = []
        with self.assertRaisesRegex(ValueError, "must have the keys"):
            build_plan.launcher_of(extra)
        optional = build_plan.launcher_of(profile_x(ReadOnlyBinds=["/data/in"], MemoryStopFraction=0.8))
        self.assertEqual(optional["ReadOnlyBinds"], ["/data/in"])
        with self.assertRaisesRegex(ValueError, "MemoryStopFraction"):
            build_plan.launcher_of(profile_x(MemoryStopFraction=1.5))
        with self.assertRaisesRegex(ValueError, "ReadOnlyBinds"):
            build_plan.launcher_of(profile_x(ReadOnlyBinds="/data/in"))
        with self.assertRaisesRegex(ValueError, "not supported"):
            build_plan.launcher_of(profile_x(Kind="docker"))
        with self.assertRaisesRegex(ValueError, "must be lists"):
            build_plan.launcher_of(profile_x(Environment=[]))

    def test_plan_fields_carry_the_block_on_one_node_with_the_image_executables(self):
        profile = profile_x()
        fields = build_plan.launcher_plan_fields(profile, binary="/opt/view/bin/palace-arm64.bin", binary_sha256=DIGEST_BINARY,
                                                 mpiexec="/opt/view/bin/mpirun", nodes=1)
        self.assertEqual(fields, {"Launcher": profile["Launcher"]})
        with self.assertRaisesRegex(ValueError, "ONE node"):
            build_plan.launcher_plan_fields(profile, binary="/opt/view/bin/palace-arm64.bin", binary_sha256=DIGEST_BINARY,
                                            mpiexec="/opt/view/bin/mpirun", nodes=2)
        with self.assertRaisesRegex(ValueError, "must name the image"):
            build_plan.launcher_plan_fields(profile, binary="/r/palace-x.bin", binary_sha256=DIGEST_BINARY,
                                            mpiexec="/opt/view/bin/mpirun", nodes=1)

    def test_job_script_drops_only_the_module_lines(self):
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        pbs = render(profile).splitlines()
        x = render(profile_x()).splitlines()
        modules = ["source /etc/profile.d/modules.sh", "module load " + " ".join(profile["Modules"])]
        self.assertTrue(all(line in pbs for line in modules))
        self.assertFalse(any("form X" in line for line in pbs))
        self.assertFalse(any(line in x for line in modules))
        comment = [line for line in x if line.startswith("# form X:")]
        self.assertEqual(len(comment), 1)
        self.assertIn(DIGEST_IMAGE, comment[0])
        self.assertEqual([line for line in pbs if line not in modules], [line for line in x if line != comment[0]])
        with self.assertRaisesRegex(ValueError, "ONE node"):
            render(profile_x(), nodes=2)


class DriverExecutablesTest(unittest.TestCase):
    def args(self, **values):
        return mock.Mock(spec=["frozen_binary_remote_path", "mpiexec_remote_path"],
                         frozen_binary_remote_path=values.get("binary"), mpiexec_remote_path=values.get("mpiexec"))

    def test_executables_are_the_image_ones(self):
        profile = profile_x()
        self.assertEqual(qualify_library.executables_of(self.args(), profile, "/r", DIGEST_BINARY),
                         ("/opt/view/bin/palace-arm64.bin", "/opt/view/bin/mpirun"))
        with self.assertRaisesRegex(ValueError, "do not apply"):
            qualify_library.executables_of(self.args(binary="/r/p.bin"), profile, "/r", DIGEST_BINARY)
        with self.assertRaisesRegex(ValueError, "not the form-X image executable"):
            qualify_library.executables_of(self.args(), profile, "/r", "3" * 64)

    def test_remote_image_check_fails_closed(self):
        profile = profile_x()
        triple = ("/opt/view/bin/palace-arm64.bin", DIGEST_BINARY, "/opt/view/bin/mpirun")
        calls = []

        def verify(exists=True, image=DIGEST_IMAGE, inner=DIGEST_BINARY, executables=(triple,)):
            return qualify_library.verify_remote_executables(
                "host", list(executables), profile=profile, remote_root="/r",
                path_exists=lambda host, path: exists, remote_sha256=lambda host, paths: {paths[0]: image},
                image_sha256=lambda host, apptainer, sif, path: calls.append((apptainer, sif, path)) or inner)

        record = verify()
        self.assertEqual(record["Image"], {"Path": "/data/sif/palace.sif", "SHA256": DIGEST_IMAGE, "Apptainer": "/bin/apptainer"})
        self.assertEqual(calls[-1], ("/bin/apptainer", "/data/sif/palace.sif", "/opt/view/bin/palace-arm64.bin"))
        with self.assertRaisesRegex(RuntimeError, "does not exist"):
            verify(exists=False)
        with self.assertRaisesRegex(RuntimeError, "has digest"):
            verify(image="0" * 64)
        with self.assertRaisesRegex(RuntimeError, "inside"):
            verify(inner="0" * 64)
        with self.assertRaisesRegex(RuntimeError, "not the image's"):
            verify(executables=[("/r/palace.bin", DIGEST_BINARY, "/r/mpiexec_bound.sh")])

    def test_dense_twin_executables_are_the_image_ones(self):
        profile = profile_x()
        found = dense_twin_plan.executables_of_record(profile, "/r", DIGEST_BINARY, exists=lambda path: True,
                                                      digest=lambda path: DIGEST_IMAGE)
        self.assertEqual(found, ("/opt/view/bin/palace-arm64.bin", "/opt/view/bin/mpirun"))
        with self.assertRaisesRegex(RuntimeError, "does not exist at digest"):
            dense_twin_plan.executables_of_record(profile, "/r", DIGEST_BINARY, exists=lambda path: True,
                                                  digest=lambda path: "0" * 64)
        with self.assertRaisesRegex(RuntimeError, "do not apply"):
            dense_twin_plan.executables_of_record(profile, "/r", DIGEST_BINARY, binary_path="/r/p.bin",
                                                  exists=lambda path: True, digest=lambda path: DIGEST_IMAGE)


FAKE_APPTAINER = '''#!/usr/bin/env python3
import json, os, sys, time
args = sys.argv[1:]
if args == ["--version"]:
    print("apptainer version 0.0-fake")
    sys.exit(0)
assert args[0] == "exec", args
index, binds, flags = 1, [], []
while args[index].startswith("--"):
    if args[index] == "--bind":
        binds.append(args[index + 1])
        index += 2
    else:
        flags.append(args[index])
        index += 1
image, rest = args[index], args[index + 1:]
# --cleanenv: the image sees only the APPTAINERENV_* variables (prefix removed).
inner = ({key[len("APPTAINERENV_"):]: value for key, value in os.environ.items() if key.startswith("APPTAINERENV_")}
         if "--cleanenv" in flags else dict(os.environ))
if rest[:1] == ["sha256sum"]:
    print(os.environ.get("FAKE_INNER_DIGEST", "%s") + "  " + rest[1])
    sys.exit(0)
if rest == ["env"]:
    for key, value in sorted(inner.items()):
        print(key + "=" + value)
    sys.exit(0)
if os.environ.get("FAKE_IGNORE_TERM_STAGE") == rest[-1]:
    # Decision 1721 F-1: a stage that ignores the cap's SIGTERM (timeout's KILL escalation ends it).
    import signal
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    open(os.environ["FAKE_IGNORE_TERM_COUNT"], "a").write("x")
    if os.environ.get("FAKE_IGNORE_TERM_FI_DOMAIN"):
        print("Open MPI failed an OFI Libfabric library call (fi_domain).  This is highly")
        print("  Error: Cannot allocate memory (1)", flush=True)
    time.sleep(60)
if os.environ.get("FAKE_SLEEP_STAGE") == rest[-1]:
    time.sleep(30)
if os.environ.get("FAKE_FI_DOMAIN_STAGE") == rest[-1]:
    counter = os.environ["FAKE_FI_DOMAIN_COUNT"]
    count = int(open(counter).read()) if os.path.exists(counter) else 0
    if count < int(os.environ["FAKE_FI_DOMAIN_TIMES"]):
        open(counter, "w").write(str(count + 1))
        if os.environ.get("FAKE_RESULT_FILE"):
            os.makedirs(os.path.dirname(os.environ["FAKE_RESULT_FILE"]), exist_ok=True)
            open(os.environ["FAKE_RESULT_FILE"], "w").write("partial\\n")
        if os.environ.get("FAKE_PALACE_STARTED"):
            print("Running with 4 MPI processes")
        print("Open MPI failed an OFI Libfabric library call (fi_domain).  This is highly")
        print("  Error: Cannot allocate memory (1)")
        sys.exit(14)
if os.environ.get("FAKE_FLAG"):
    open(os.environ["FAKE_FLAG"], "w").close()
    time.sleep(30)
with open(os.environ["FAKE_RECORD"], "a") as stream:
    stream.write(json.dumps({"argv": args, "env": inner}) + "\\n")
print("H1 (p = 3): 1,000, ND (p = 3): 2,000, RT (p = 3): 3,000")
''' % DIGEST_BINARY


class RunStagesLauncherTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.apptainer = self.root / "apptainer"
        self.apptainer.write_text(FAKE_APPTAINER)
        self.apptainer.chmod(0o755)
        self.image = self.root / "palace.sif"
        self.image.write_bytes(b"not an image")
        self.record = self.root / "launches.jsonl"
        self.saved = dict(os.environ)

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.saved)
        self.tmp.cleanup()

    def plan(self, *, nodes=None, stage_binary=False, image_sha256=None):
        stages = []
        for name, environment, requires in (("p3-worker", {"PALACE_RESPONSE_ARCHIVE_ONLY": "1",
                                                           "PALACE_RESPONSE_ARCHIVE_DIR": str(self.root / "archive")}, []),
                                            ("p3-reducer", {"PALACE_RESPONSE_REDUCE_ONLY": "1", "PALACE_RESPONSE_BLOCK_SIZE": "3",
                                                            "PALACE_RESPONSE_ARCHIVE_DIR": str(self.root / "archive")}, ["p3-worker"])):
            config = self.root / f"{name}.json"
            config.write_text(json.dumps({"Problem": {"Output": str(self.root / f"out-{name}")}}))
            stages.append({"Name": name, "Config": str(config), "Environment": environment, "CapSeconds": 60,
                           "MinimumSeconds": 0, "Requires": requires})
        if stage_binary:
            stages[0].update(Binary="/r/other.bin", BinarySHA256="4" * 64)
        block = launcher(Apptainer=str(self.apptainer), Image=str(self.image),
                         ImageSHA256=image_sha256 or run_stages.sha(self.image))
        plan = {"Version": 1, "Case": "x", "Ranks": 4, "DeadlineSeconds": 600, "DeadlineMarginSeconds": 1,
                "Instance": {"Type": "m8g.48xlarge"}, "MinimumMemAvailableBytes": 0, "Binary": block["ImageBinary"],
                "BinarySHA256": DIGEST_BINARY, "MPIExec": block["ImageMPIExec"], "PinnedSHA256": {}, "Stages": stages,
                "Launcher": block}
        if nodes:
            plan["Nodes"] = nodes
        directory = self.root / "main"
        directory.mkdir()
        (directory / "plan.json").write_text(json.dumps(plan))
        return directory / "plan.json"

    def run_plan(self, plan_path, hosts=None):
        host = socket.gethostname().split(".")[0]
        (self.root / "nodes").write_text("\n".join(hosts or [host]) + "\n")
        os.environ.update(PBS_NODEFILE=str(self.root / "nodes"), PBS_JOBID="1.fake", FAKE_RECORD=str(self.record),
                          PALACE_RESPONSE_DERIVED_CACHE="/inherited")
        real_command = run_stages.launcher_command

        def command(*args, **kwargs):
            # A host without GNU time (/usr/bin/time: the SOCA nodes have it) runs the launch
            # without its "/usr/bin/time -v -o FILE" part; test_capped_command checks that part.
            parts = real_command(*args, **kwargs)
            return parts if Path("/usr/bin/time").exists() else parts[:4] + parts[8:]

        with mock.patch.object(run_stages, "conflicting_processes", return_value=[]), \
                mock.patch.object(run_stages, "launcher_command", side_effect=command), \
                contextlib.redirect_stdout(io.StringIO()):
            return run_stages.main(["run_stages.py", str(plan_path)])

    def test_capped_command(self):
        block = launcher()
        parts = run_stages.launcher_command(block, cap=60.7, timefile=Path("/t/s.time"), exports=["-x", "PALACE_RESPONSE_X"],
                                            ranks=8, config="c.json")
        self.assertEqual(parts, ["timeout", "-k", "30", "60", "/usr/bin/time", "-v", "-o", "/t/s.time",
                                 *run_stages.launcher_launch(block, exports=["-x", "PALACE_RESPONSE_X"], ranks=8,
                                                             arguments=["c.json"])])

    def test_stages_run_inside_the_image_with_their_variables(self):
        plan_path = self.plan()
        self.assertEqual(self.run_plan(plan_path), 0)
        status = json.loads((plan_path.parent / "status.json").read_text())
        self.assertEqual(status["State"], "complete")
        self.assertEqual(status["Preflight"]["BinarySHA256"], DIGEST_BINARY)
        self.assertEqual(status["Preflight"]["Launcher"]["ImageSHA256"], run_stages.sha(self.image))
        launches = [json.loads(line) for line in self.record.read_text().splitlines()]
        self.assertEqual(len(launches), 2)
        worker, reducer = launches
        block = json.loads(plan_path.read_text())["Launcher"]
        self.assertEqual(worker["argv"], ["exec", "--cleanenv", "--no-home", "--bind", "/data/home/simlap", str(self.image),
                                          "/opt/view/bin/mpirun",
                                          *block["MPIExecArguments"], "-x", "GLIBC_TUNABLES",
                                          "-x", "PALACE_RESPONSE_ARCHIVE_ONLY", "-x", "PALACE_RESPONSE_ARCHIVE_DIR",
                                          "-n", "4", "/opt/view/bin/palace-arm64.bin", str(self.root / "p3-worker.json")])
        passed = {key: value for key, value in self.saved.items() if run_stages.LAUNCHER_PASSED_PATTERN.match(key)}
        self.assertEqual(worker["env"], {"PALACE_RESPONSE_ARCHIVE_ONLY": "1", "PALACE_RESPONSE_ARCHIVE_DIR": str(self.root / "archive"),
                                         "GLIBC_TUNABLES": "glibc.malloc.hugetlb=1", **passed})
        self.assertEqual(status["Preflight"]["Launcher"]["ApptainerVersion"], "apptainer version 0.0-fake")
        self.assertEqual(status["Preflight"]["Launcher"]["ImageEnvironment"], {"GLIBC_TUNABLES": "glibc.malloc.hugetlb=1", **passed})
        self.assertEqual(status["Preflight"]["MemoryStopBarBytes"], int(0.85 * status["Preflight"]["MemTotalBytes"]))
        self.assertEqual(reducer["env"]["PALACE_RESPONSE_REDUCE_ONLY"], "1")
        self.assertNotIn("PALACE_RESPONSE_ARCHIVE_ONLY", reducer["env"])
        self.assertEqual(status["Stages"][0]["Parsed"]["H1"], 1000)

    def test_refusals(self):
        host = socket.gethostname().split(".")[0]
        with self.assertRaisesRegex(SystemExit, "one node only"):
            self.run_plan(self.plan(nodes=2), hosts=[host, "other-node"])
        self.tearDown()
        self.setUp()
        with self.assertRaisesRegex(SystemExit, "its own executable"):
            self.run_plan(self.plan(stage_binary=True))
        self.tearDown()
        self.setUp()
        with self.assertRaisesRegex(SystemExit, "Image hash mismatch"):
            self.run_plan(self.plan(image_sha256="0" * 64))
        self.tearDown()
        self.setUp()
        os.environ["FAKE_INNER_DIGEST"] = "0" * 64
        with self.assertRaisesRegex(SystemExit, "Executable hash mismatch"):
            self.run_plan(self.plan())
        self.assertFalse(self.record.exists())


    def test_in_image_digest_is_checked_against_the_launcher_pin(self):
        plan_path = self.plan()
        plan = json.loads(plan_path.read_text())
        plan["Launcher"]["ImageBinarySHA256"] = "5" * 64
        plan_path.write_text(json.dumps(plan))
        with self.assertRaisesRegex(SystemExit, "the Launcher pin"):
            self.run_plan(plan_path)
        self.assertFalse(self.record.exists())

    def test_read_only_binds(self):
        prefix = run_stages.launcher_launch_prefix(launcher(ReadOnlyBinds=["/data/in"]))
        self.assertEqual(prefix, ["/bin/apptainer", "exec", "--cleanenv", "--no-home", "--bind", "/data/home/simlap",
                                  "--bind", "/data/in:/data/in:ro", "/data/sif/palace.sif"])

    def test_memory_stop_rule(self):
        block = launcher()
        self.assertEqual(run_stages.memory_stop_bar(block, 1000), 850)
        self.assertEqual(run_stages.memory_stop_bar(launcher(MemoryStopFraction=0.5), 1000), 500)
        self.assertIsNone(run_stages.memory_stop_index([100, 850, 200], 850))
        self.assertEqual(run_stages.memory_stop_index([100, 851, 900], 850), 1)

    def fake_meminfo(self, flag):
        total = 100 << 30
        return lambda: {"MemTotal": total, "MemAvailable": total - ((90 << 30) if flag.exists() else (1 << 30))}

    def test_memory_stop_ends_the_running_stage_and_the_job(self):
        flag = self.root / "stage-started"
        os.environ["FAKE_FLAG"] = str(flag)
        plan_path = self.plan()
        with mock.patch.object(run_stages, "meminfo", side_effect=self.fake_meminfo(flag)):
            started = time.monotonic()
            self.assertEqual(self.run_plan(plan_path), 3)
            self.assertLess(time.monotonic() - started, 25)
        status = json.loads((plan_path.parent / "status.json").read_text())
        self.assertEqual(status["State"], "memory-stop")
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["memory-stop"])
        self.assertNotEqual(status["Stages"][0]["ReturnCode"], 0)
        stop = json.loads((plan_path.parent / "memory-stop.json").read_text())
        self.assertEqual(stop["Stage"], "p3-worker")
        self.assertEqual(stop["BarBytes"], int(0.85 * (100 << 30)))
        self.assertGreater(stop["UsedBytes"], stop["BarBytes"])
        self.assertTrue(stop["Samples"])
        self.assertFalse(self.record.exists())

    # Decision 1710: OI-AMRFORM-EFA-AFTER-KILL under form X (PBS 118751): the wait after a killed
    # stage, its 300-s stop, and the ONE retry of a stage that failed at MPI start in fi_domain.
    FI_DOMAIN_LOG = ("Open MPI failed an OFI Libfabric library call (fi_domain).  This is highly\n"
                     "  Error: Cannot allocate memory (98011153694732)\n")

    def test_session_processes_reads_the_session_field(self):
        proc = self.root / "proc"
        for pid, stat in (("10", "10 (a b) S 1 10 42 0"), ("11", "11 (palace) R 10 11 42 0"),
                          ("12", "12 (x)) y) S 1 12 7 0"), ("self", "13 (s) S 1 13 42 0")):
            (proc / pid).mkdir(parents=True)
            (proc / pid / "stat").write_text(stat + "\n")
        (proc / "14").mkdir()  # exited while the table was read: no stat file
        self.assertEqual(run_stages.session_processes(42, proc), [10, 11])
        self.assertEqual(run_stages.session_processes(7, proc), [12])

    def test_wait_for_session_exit_and_its_300_s_stop(self):
        now = [0.0]

        def sleep(seconds):
            now[0] += seconds

        answers = [[5, 6], [5], []]
        exited = run_stages.wait_for_session_exit(42, processes=lambda session: answers.pop(0), sleep=sleep,
                                                  clock=lambda: now[0])
        self.assertEqual(exited, {"Session": 42, "ExitWaitSeconds": 2.0, "Remaining": []})
        now[0] = 0.0
        stuck = run_stages.wait_for_session_exit(42, processes=lambda session: [5], sleep=sleep, clock=lambda: now[0])
        self.assertEqual((stuck["Remaining"], stuck["ExitWaitSeconds"]), ([5], 300.0))
        self.assertEqual((run_stages.AFTER_KILL_EXIT_SECONDS, run_stages.AFTER_KILL_PAUSE_SECONDS), (300, 90))

    def test_retry_conditions_one_by_one(self):
        log = "PRTE has detected that the head of the session directory tree\n" + self.FI_DOMAIN_LOG
        self.assertIn("fi_domain", run_stages.retry_reason(14, log, []))
        self.assertIsNone(run_stages.retry_reason(0, log, []))  # rc 0
        self.assertIsNone(run_stages.retry_reason(124, log, []))  # ended by its cap
        self.assertIsNone(run_stages.retry_reason(14, self.FI_DOMAIN_LOG.splitlines()[1] + "\n", []))  # no fi_domain line
        self.assertIsNone(run_stages.retry_reason(14, self.FI_DOMAIN_LOG.splitlines()[0] + "\n", []))  # no allocation line
        self.assertIsNone(run_stages.retry_reason(14, "Running with 192 MPI processes\n" + log, []))  # Palace had started
        self.assertIsNone(run_stages.retry_reason(14, log, ["/out/palace.json"]))  # a result file exists

    def test_stage_result_files(self):
        out, archive, config = self.root / "out", self.root / "archive", self.root / "stage.json"
        config.write_text(json.dumps({"Problem": {"Output": str(out)}}))
        (out / "sub").mkdir(parents=True)
        (out / "sub" / "a.csv").write_text("x")
        archive.mkdir()
        (archive / "source-000001-rank-000000-V.bin").write_bytes(b"1")
        writer = {"Config": str(config), "Environment": {"PALACE_RESPONSE_ARCHIVE_ONLY": "1",
                                                         "PALACE_RESPONSE_ARCHIVE_DIR": str(archive)}}
        reader = {"Config": str(config), "Environment": {"PALACE_RESPONSE_REDUCE_ONLY": "1",
                                                         "PALACE_RESPONSE_ARCHIVE_DIR": str(archive)}}
        self.assertEqual(len(run_stages.stage_result_files(writer)), 2)
        self.assertEqual(run_stages.stage_result_files(reader), [str(out / "sub" / "a.csv")])

    def fi_domain_plan(self, times, **extra):
        plan_path = self.plan()
        worker_config = json.loads(plan_path.read_text())["Stages"][0]["Config"]
        os.environ.update(FAKE_FI_DOMAIN_STAGE=worker_config, FAKE_FI_DOMAIN_COUNT=str(self.root / "fi-count"),
                          FAKE_FI_DOMAIN_TIMES=str(times), **extra)
        return plan_path

    def status_of(self, plan_path):
        return json.loads((plan_path.parent / "status.json").read_text())

    def test_one_retry_after_a_fi_domain_failure_at_mpi_start(self):
        plan_path = self.fi_domain_plan(1)
        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 0):
            self.assertEqual(self.run_plan(plan_path), 0)
        status = self.status_of(plan_path)
        worker = status["Stages"][0]
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["complete", "complete"])
        self.assertEqual((worker["Attempts"], worker["Attempt1"]["ReturnCode"]), (2, 14))
        self.assertIn("fi_domain", worker["RetryReason"])
        self.assertEqual(worker["Attempt1"]["Exit"]["Remaining"], [])
        log, kept = plan_path.parent / "p3-worker.log", plan_path.parent / "p3-worker-attempt1.log"
        self.assertEqual(worker["Attempt1"]["Kept"][str(log)], str(kept))
        self.assertIn("fi_domain", kept.read_text())
        self.assertNotIn("fi_domain", log.read_text())
        self.assertNotIn("Attempts", status["Stages"][1])
        self.assertEqual(len(self.record.read_text().splitlines()), 2)

    def test_a_second_fi_domain_failure_is_final(self):
        plan_path = self.fi_domain_plan(5)
        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 0):
            self.assertEqual(self.run_plan(plan_path), 1)
        status = self.status_of(plan_path)
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["failed", "skipped-dependency"])
        self.assertEqual(status["Stages"][0]["Attempts"], 2)
        self.assertEqual((self.root / "fi-count").read_text(), "2")

    def assert_no_retry(self, plan_path):
        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 0):
            self.assertEqual(self.run_plan(plan_path), 1)
        status = self.status_of(plan_path)
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["failed", "skipped-dependency"])
        self.assertNotIn("Attempts", status["Stages"][0])
        self.assertEqual((self.root / "fi-count").read_text(), "1")

    def test_no_retry_when_palace_had_started(self):
        self.assert_no_retry(self.fi_domain_plan(1, FAKE_PALACE_STARTED="1"))

    def test_no_retry_when_a_result_file_exists(self):
        self.assert_no_retry(self.fi_domain_plan(1, FAKE_RESULT_FILE=str(self.root / "out-p3-worker" / "palace.json")))

    def capped_plan(self):
        plan_path = self.plan()
        plan = json.loads(plan_path.read_text())
        plan["Stages"][0]["CapSeconds"] = 2
        plan["Stages"][1]["Requires"] = []
        plan_path.write_text(json.dumps(plan))
        os.environ["FAKE_SLEEP_STAGE"] = plan["Stages"][0]["Config"]
        return plan_path

    def test_after_a_cap_kill_the_session_is_empty_and_the_next_launch_waits(self):
        plan_path = self.capped_plan()
        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 1):
            self.assertEqual(self.run_plan(plan_path), 1)
        worker, reducer = self.status_of(plan_path)["Stages"]
        self.assertEqual((worker["State"], worker["ReturnCode"], worker["AfterKill"]["Remaining"]), ("timed-out", 124, []))
        self.assertNotIn("Attempts", worker)
        self.assertEqual((reducer["State"], reducer["AfterKillPauseSeconds"]), ("complete", 1))
        self.assertEqual(len(self.record.read_text().splitlines()), 1)

    def test_processes_left_after_a_kill_stop_the_job(self):
        plan_path = self.capped_plan()
        with mock.patch.object(run_stages, "AFTER_KILL_EXIT_SECONDS", 1), \
                mock.patch.object(run_stages, "session_processes", return_value=[99999]):
            self.assertEqual(self.run_plan(plan_path), 4)
        status = self.status_of(plan_path)
        self.assertEqual((status["State"], status["AfterKillStop"]["Remaining"]), ("after-kill-stop", [99999]))
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["timed-out"])
        self.assertFalse(self.record.exists())

    # Decision 1721 F-1: a cap is recognised by the wall, whatever the rc the kill produced.
    def test_stage_capped_rule(self):
        self.assertTrue(run_stages.stage_capped(124, 20.0, 20))
        self.assertTrue(run_stages.stage_capped(-9, 50.1, 20))
        self.assertTrue(run_stages.stage_capped(137, 50.0, 20.7))  # the cap given to timeout = 20
        self.assertTrue(run_stages.stage_capped(1, 20.3, 20))  # a stage at its cap, whatever its rc
        self.assertFalse(run_stages.stage_capped(0, 20.01, 20))  # rc 0: the command ended before the cap
        self.assertFalse(run_stages.stage_capped(14, 3.0, 20))
        self.assertEqual([run_stages.cap_end(code) for code in (124, -9, 137, 1)], ["SIGTERM", "SIGKILL", "SIGKILL", "rc 1"])
        self.assertIsNone(run_stages.retry_reason(-9, self.FI_DOMAIN_LOG, [], capped=True))
        self.assertIsNotNone(run_stages.retry_reason(-9, self.FI_DOMAIN_LOG, [], capped=False))

    def term_ignoring_plan(self, **extra):
        plan_path = self.capped_plan()
        os.environ.pop("FAKE_SLEEP_STAGE")
        os.environ.update(FAKE_IGNORE_TERM_STAGE=json.loads(plan_path.read_text())["Stages"][0]["Config"],
                          FAKE_IGNORE_TERM_COUNT=str(self.root / "term-count"), **extra)
        return plan_path

    def assert_kill_escalation_cap(self, plan_path):
        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 1):
            self.assertEqual(self.run_plan(plan_path), 1)
        worker, reducer = self.status_of(plan_path)["Stages"]
        self.assertIn(worker["ReturnCode"], (-9, 137))
        self.assertGreaterEqual(worker["WallSeconds"], 30)
        self.assertEqual((worker["State"], worker["TimedOut"], worker["CapEnd"]), ("timed-out", True, "SIGKILL"))
        self.assertEqual(worker["AfterKill"]["Remaining"], [])
        self.assertNotIn("Attempts", worker)
        self.assertEqual((self.root / "term-count").read_text(), "x")
        self.assertEqual((reducer["State"], reducer["AfterKillPauseSeconds"]), ("complete", 1))

    def test_a_cap_ended_by_the_kill_escalation_is_a_cap(self):
        self.assert_kill_escalation_cap(self.term_ignoring_plan())

    def test_a_cap_kill_with_the_fi_domain_text_is_never_retried(self):
        self.assert_kill_escalation_cap(self.term_ignoring_plan(FAKE_IGNORE_TERM_FI_DOMAIN="1"))

    # Decision 1721 F-2 / F-5: the memory stop fails closed.
    def test_a_sampler_failure_kills_the_stage_and_ends_the_job(self):
        flag = self.root / "stage-started"
        os.environ["FAKE_FLAG"] = str(flag)
        plan_path = self.plan()
        total = 100 << 30

        def failing_meminfo():
            if flag.exists():
                raise OSError("a sampler read fails")
            return {"MemTotal": total, "MemAvailable": total - (1 << 30)}

        with mock.patch.object(run_stages, "meminfo", side_effect=failing_meminfo):
            started = time.monotonic()
            self.assertEqual(self.run_plan(plan_path), 5)
            self.assertLess(time.monotonic() - started, 25)
        status = self.status_of(plan_path)
        self.assertEqual(status["State"], "sampler-failed")
        self.assertIn("a sampler read fails", status["SamplerFailed"]["Reason"])
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["sampler-failed"])
        self.assertNotEqual(status["Stages"][0]["ReturnCode"], 0)
        self.assertFalse(self.record.exists())

    def test_a_sampler_not_alive_blocks_every_launch(self):
        plan_path = self.plan()
        with mock.patch.object(run_stages.threading.Thread, "is_alive", return_value=False):
            self.assertEqual(self.run_plan(plan_path), 5)
        status = self.status_of(plan_path)
        self.assertEqual(status["State"], "sampler-failed")
        self.assertEqual([(stage["State"], stage.get("NotLaunched")) for stage in status["Stages"]], [("sampler-failed", True)])
        self.assertIn("not alive", status["SamplerFailed"]["Reason"])
        self.assertFalse(self.record.exists())

    def test_a_stop_file_write_failure_still_kills_the_stage(self):
        flag = self.root / "stage-started"
        os.environ["FAKE_FLAG"] = str(flag)
        plan_path = self.plan()
        (plan_path.parent / "memory-stop.json").mkdir()  # the write of the stop file fails
        with mock.patch.object(run_stages, "meminfo", side_effect=self.fake_meminfo(flag)):
            started = time.monotonic()
            self.assertEqual(self.run_plan(plan_path), 3)
            self.assertLess(time.monotonic() - started, 25)
        status = self.status_of(plan_path)
        self.assertEqual([stage["State"] for stage in status["Stages"]], ["memory-stop"])
        self.assertNotEqual(status["Stages"][0]["ReturnCode"], 0)
        self.assertIn("IsADirectoryError", status["MemoryStop"]["StopFileError"])
        self.assertFalse(self.record.exists())

    def test_a_crossing_during_the_after_kill_pause_blocks_the_next_launch(self):
        plan_path = self.capped_plan()
        flag, real_sleep = self.root / "in-pause", time.sleep

        def sleep(seconds):
            if seconds == 6:  # the after-kill pause
                flag.touch()
            real_sleep(seconds)

        with mock.patch.object(run_stages, "AFTER_KILL_PAUSE_SECONDS", 6), \
                mock.patch.object(run_stages, "meminfo", side_effect=self.fake_meminfo(flag)), \
                mock.patch.object(run_stages.time, "sleep", side_effect=sleep):
            self.assertEqual(self.run_plan(plan_path), 3)
        status = self.status_of(plan_path)
        self.assertEqual([(stage["State"], stage.get("NotLaunched")) for stage in status["Stages"]],
                         [("timed-out", None), ("memory-stop", True)])
        self.assertEqual(json.loads((plan_path.parent / "memory-stop.json").read_text())["Stage"], "after-kill-pause")
        self.assertFalse(self.record.exists())

    def test_the_plan_and_the_runner_are_hashed(self):
        plan_path = self.plan()
        self.assertEqual(self.run_plan(plan_path), 0)
        preflight = self.status_of(plan_path)["Preflight"]
        self.assertEqual(preflight["PlanSHA256"], run_stages.sha(plan_path))
        self.assertEqual(preflight["RunnerSHA256"], run_stages.sha(Path(run_stages.__file__).resolve()))


class PalaceInImageTest(unittest.TestCase):
    def test_more_than_one_node_is_refused(self):
        with tempfile.TemporaryDirectory() as root:
            nodefile = Path(root) / "nodes"
            nodefile.write_text("node-a\nnode-a\nnode-b\n")
            with self.assertRaisesRegex(SystemExit, "2 nodes in PBS_NODEFILE"):
                palace_in_image.image_launch(profile_x(), 96, ["c.json"], {"PBS_NODEFILE": str(nodefile)},
                                             digest=lambda path: DIGEST_IMAGE,
                                             run=lambda command: f"{DIGEST_BINARY}  {command[-1]}\n")

    def test_launch_command_and_environment(self):
        profile = profile_x()
        command, environment, check = palace_in_image.image_launch(
            profile, 96, ["--surface-response-preflight", "c0-p4.json"],
            {"PBS_JOBID": "1.fake", "PALACE_RESPONSE_ARCHIVE_DIR": "/a", "HOME": "/h"},
            digest=lambda path: DIGEST_IMAGE, run=lambda command: f"{DIGEST_BINARY}  {command[-1]}\n")
        self.assertEqual(command, ["/bin/apptainer", "exec", "--cleanenv", "--no-home", "--bind", "/data/home/simlap",
                                   "/data/sif/palace.sif",
                                   "/opt/view/bin/mpirun", *profile["Launcher"]["MPIExecArguments"], "-x", "GLIBC_TUNABLES",
                                   "-x", "PALACE_RESPONSE_ARCHIVE_DIR", "-n", "96", "/opt/view/bin/palace-arm64.bin",
                                   "--surface-response-preflight", "c0-p4.json"])
        self.assertEqual(environment, {"PALACE_RESPONSE_ARCHIVE_DIR": "/a", "HOME": "/h", "GLIBC_TUNABLES": "glibc.malloc.hugetlb=1",
                                       "APPTAINERENV_PALACE_RESPONSE_ARCHIVE_DIR": "/a",
                                       "APPTAINERENV_GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"})
        self.assertEqual(check["BinarySHA256"], DIGEST_BINARY)
        with self.assertRaisesRegex(SystemExit, "Executable hash mismatch"):
            palace_in_image.image_launch(profile, 96, ["c.json"], {}, digest=lambda path: DIGEST_IMAGE,
                                         run=lambda command: f"{'0' * 64}  {command[-1]}\n")

    def test_profile_without_launcher_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "profile.json"
            path.write_text((HERE / "qualify" / "cluster-profile.json").read_text())
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(palace_in_image.main(["palace_in_image.py", str(path), "8", "c.json"]), 2)


if __name__ == "__main__":
    unittest.main()
