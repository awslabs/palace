# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The Julia runtime resolver: the pin names the binary a launcher runs (decision 1014)."""
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from julia_runtime import JULIA_RUNTIME_QUERY, resolve_julia_runtime


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class ResolveJuliaRuntimeTest(unittest.TestCase):
    def test_fake_launcher_resolves_to_the_binary_it_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "julia-9.9.9" / "bin" / "julia"
            binary.parent.mkdir(parents=True)
            binary.write_bytes(b"julia binary")
            launcher = root / "julialauncher"
            launcher.write_text(f"#!/bin/sh\nprintf '%s\\n9.9.9\\n' '{binary}'\n")
            launcher.chmod(0o755)
            self.assertEqual(resolve_julia_runtime(launcher), (binary.resolve(), "9.9.9"))
            # A symlink to the launcher (the PATH `julia` of a juliaup install) resolves alike.
            link = root / "julia"
            link.symlink_to(launcher)
            self.assertEqual(resolve_julia_runtime(link), (binary.resolve(), "9.9.9"))
            # The launcher's own bytes are not the identity: a replaced launcher naming the
            # same binary resolves to the same binary.
            launcher.write_text(f"#!/bin/sh\n# self-updated\nprintf '%s\\n9.9.9\\n' '{binary}'\n")
            self.assertEqual(resolve_julia_runtime(link), (binary.resolve(), "9.9.9"))
            self.assertNotEqual(sha256(launcher), sha256(binary))

    def test_errors_are_value_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, "does not exist"):
                resolve_julia_runtime(root / "missing")
            launcher = root / "julialauncher"
            launcher.write_text(f"#!/bin/sh\nprintf '%s\\n9.9.9\\n' '{root / 'absent'}'\n")
            launcher.chmod(0o755)
            with self.assertRaisesRegex(ValueError, "missing executable"):
                resolve_julia_runtime(launcher)
            launcher.write_text("#!/bin/sh\nexit 3\n")
            with self.assertRaisesRegex(ValueError, "did not report"):
                resolve_julia_runtime(launcher)
            launcher.write_text("#!/bin/sh\necho only-one-line\n")
            with self.assertRaisesRegex(ValueError, "did not report"):
                resolve_julia_runtime(launcher)
            not_executable = root / "plain"
            not_executable.write_bytes(b"julia launcher")
            with self.assertRaisesRegex(ValueError, "did not report"):
                resolve_julia_runtime(not_executable)

    def test_the_machine_julia_resolves_to_a_fixed_point(self):
        julia = shutil.which(os.environ.get("JULIA", "julia"))
        if julia is None:
            self.skipTest("Julia is not available")
        executable, version = resolve_julia_runtime(julia)
        self.assertTrue(executable.is_file())
        self.assertEqual(executable, executable.resolve())
        # The binary resolves to itself, and reports the version `julia --version` prints.
        self.assertEqual(resolve_julia_runtime(executable), (executable, version))
        printed = subprocess.run([julia, "--startup-file=no", "--version"], capture_output=True,
                                 text=True, check=True).stdout.strip()
        self.assertEqual(printed, f"julia version {version}")
        # The query is what the binary itself reports for its command.
        command = subprocess.run([str(executable), "--startup-file=no", "-e",
                                  "print(Base.julia_cmd()[1])"], capture_output=True, text=True,
                                 check=True).stdout.strip()
        self.assertEqual(Path(command).resolve(), executable)
        self.assertIn("Sys.BINDIR", JULIA_RUNTIME_QUERY)


if __name__ == "__main__":
    unittest.main()
