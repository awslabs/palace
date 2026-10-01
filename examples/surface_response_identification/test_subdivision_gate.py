# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""subdivision_gate exit code (decision 214 (i)): the verdict under EVERY variant. A layout
passing the collinear refinements but failing the rotation variant fails the gate (exit 1)
while the collinear-only line still counts it; every variant passing exits 0. The per-layout
identification is replaced by canned rows (no Palace run)."""
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import subdivision_gate as G  # noqa: E402

LABELS = ["refine-0.5", "refine-0.37", "rotate-37", "rotate-37-refine-0.5"]


def row(name, failing=()):
    variants = {"base": {"Digest": "d0"}}
    for label in LABELS:
        variants[label] = {"Digest": "d0", "Pass": label not in failing}
    return {"Layout": name, "Mesh": name + ".msh2", "Variants": variants, "Pass": not failing, "Seconds": 0.1}


class SubdivisionGateExitCodeTest(unittest.TestCase):
    def run_gate(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            meshes = os.path.join(directory, "meshes")
            os.makedirs(meshes)
            for r in rows:
                with open(os.path.join(meshes, r["Layout"] + ".msh2"), "w") as target:
                    target.write("")
            canned = {r["Layout"]: r for r in rows}
            output = os.path.join(directory, "out")
            printed = io.StringIO()
            with mock.patch.object(G, "run_layout", lambda name, path, lay, args: canned[name]), contextlib.redirect_stdout(printed):
                code = G.main(["--output", output, "--meshes", meshes])
            with open(os.path.join(output, "results.json")) as source:
                results = json.load(source)
        return code, results, printed.getvalue()

    def test_rotation_failure_fails_the_gate_but_not_the_collinear_line(self):
        code, results, printed = self.run_gate([row("stack-curved-k4-rho3", failing=("rotate-37", "rotate-37-refine-0.5")), row("fillet")])
        self.assertEqual(code, 1)
        self.assertEqual((results["SubdivisionPass"], results["Pass"], results["Total"]), (2, 1, 2))
        self.assertEqual(results["PerVariant"], {"refine-0.5": 2, "refine-0.37": 2, "rotate-37": 1, "rotate-37-refine-0.5": 1})
        self.assertIn("SUBDIVISION GATE 2 / 2 layouts PASS under collinear subdivision (refine-0.5, refine-0.37); 1 / 2 under every variant", printed)

    def test_every_variant_passing_exits_zero(self):
        code, results, printed = self.run_gate([row("stack-curved-k4-rho3"), row("fillet")])
        self.assertEqual(code, 0)
        self.assertEqual((results["SubdivisionPass"], results["Pass"], results["Total"]), (2, 2, 2))
        self.assertIn("2 / 2 under every variant", printed)

    def test_collinear_failure_fails_both_counts(self):
        code, results, _ = self.run_gate([row("joint", failing=("refine-0.37",))])
        self.assertEqual(code, 1)
        self.assertEqual((results["SubdivisionPass"], results["Pass"]), (0, 0))


if __name__ == "__main__":
    unittest.main()
