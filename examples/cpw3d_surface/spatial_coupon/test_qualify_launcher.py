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
        with self.assertRaisesRegex(ValueError, "exactly the keys"):
            build_plan.launcher_of(broken)
        extra = profile_x()
        extra["Launcher"]["Modules"] = []
        with self.assertRaisesRegex(ValueError, "exactly the keys"):
            build_plan.launcher_of(extra)
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
import json, os, sys
args = sys.argv[1:]
assert args[0] == "exec", args
index, binds = 1, []
while args[index] == "--bind":
    binds.append(args[index + 1])
    index += 2
image, rest = args[index], args[index + 1:]
if rest[:1] == ["sha256sum"]:
    print(os.environ.get("FAKE_INNER_DIGEST", "%s") + "  " + rest[1])
    sys.exit(0)
with open(os.environ["FAKE_RECORD"], "a") as stream:
    stream.write(json.dumps({"argv": args, "env": {key: value for key, value in os.environ.items()
                                                   if key.startswith(("PALACE_RESPONSE_", "PBS_", "GLIBC_TUNABLES"))}}) + "\\n")
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
        self.assertEqual(worker["argv"], ["exec", "--bind", "/data/home/simlap", str(self.image), "/opt/view/bin/mpirun",
                                          *block["MPIExecArguments"], "-x", "GLIBC_TUNABLES",
                                          "-x", "PALACE_RESPONSE_ARCHIVE_ONLY", "-x", "PALACE_RESPONSE_ARCHIVE_DIR",
                                          "-n", "4", "/opt/view/bin/palace-arm64.bin", str(self.root / "p3-worker.json")])
        self.assertEqual(worker["env"], {"PALACE_RESPONSE_ARCHIVE_ONLY": "1", "PALACE_RESPONSE_ARCHIVE_DIR": str(self.root / "archive"),
                                         "GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"})
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


class PalaceInImageTest(unittest.TestCase):
    def test_launch_command_and_environment(self):
        profile = profile_x()
        command, environment, check = palace_in_image.image_launch(
            profile, 96, ["--surface-response-preflight", "c0-p4.json"],
            {"PBS_JOBID": "1.fake", "PALACE_RESPONSE_ARCHIVE_DIR": "/a", "HOME": "/h"},
            digest=lambda path: DIGEST_IMAGE, run=lambda command: f"{DIGEST_BINARY}  {command[-1]}\n")
        self.assertEqual(command, ["/bin/apptainer", "exec", "--bind", "/data/home/simlap", "/data/sif/palace.sif",
                                   "/opt/view/bin/mpirun", *profile["Launcher"]["MPIExecArguments"], "-x", "GLIBC_TUNABLES",
                                   "-x", "PALACE_RESPONSE_ARCHIVE_DIR", "-n", "96", "/opt/view/bin/palace-arm64.bin",
                                   "--surface-response-preflight", "c0-p4.json"])
        self.assertEqual(environment, {"PALACE_RESPONSE_ARCHIVE_DIR": "/a", "HOME": "/h", "GLIBC_TUNABLES": "glibc.malloc.hugetlb=1"})
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
