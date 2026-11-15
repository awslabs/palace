#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""check_surface_response_output.py on synthetic Maxwell trees (decision 1847 E-1 (c)): a valid corrected tree passes (exit 0);
the MarkNotApplied form of the same tree (corrected cells nan, their values in the ' partial (uncovered part raw)' twins) reads
NOT-APPLIED (exit 3, never the 'valid' line; CorrectionNotApplied is a RuntimeError, so callers of check() fail closed); a nan
without a twin keeps the 'Invalid ...' error."""

import contextlib
import io
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

CPW2D = Path(__file__).resolve().parent
sys.path.insert(0, str(CPW2D))
import check_surface_response_output as CSRO  # noqa: E402

CORRECTED = ("E_norm corrected (J)", "E_surf corrected (J)", "p_surf corrected", "Q_surf corrected")


def write_tree(root, marked=False, nan_without_twin=False):
    """One frequency, one excitation, interfaces 1 / 2. marked: the MarkNotApplied form (every corrected column nan + its twin)."""
    post = Path(root) / "postpro"
    post.mkdir(parents=True)
    norm_raw, norm_corr = 1.0e-12, 1.1e-12
    e_raw = {1: 1.0e-16, 2: 2.0e-16}
    e_corr = {1: 1.2e-16, 2: 2.3e-16}
    raw_header = ["f (GHz)"] + [f"{q}[{k}]" for k in (1, 2) for q in ("p_surf", "Q_surf")]
    raw_row = [5.0] + [v for k in (1, 2) for v in (e_raw[k] / norm_raw, norm_raw / e_raw[k])]
    (post / "surface-Q.csv").write_text(",".join(raw_header) + "\n" + ",".join(f"{v:+.12e}" for v in raw_row) + "\n")
    header = ["f (GHz)", "exc", "interface", "E_norm raw (J)", "E_norm corrected (J)", "E_surf raw (J)", "p_surf raw", "Q_surf raw",
              "E_surf corrected (J)", "p_surf corrected", "Q_surf corrected"]
    rows = []
    for k in (1, 2):
        rows.append([5.0, 1, k, norm_raw, norm_corr, e_raw[k], e_raw[k] / norm_raw, norm_raw / e_raw[k], e_corr[k],
                     e_corr[k] / norm_corr, norm_corr / e_corr[k]])
    cols = [header.index(h) for h in CORRECTED]
    out_header = list(header)
    if marked:
        out_header += [header[j] + CSRO.PARTIAL_COLUMN_SUFFIX for j in cols]
    text = [",".join(out_header)]
    for r in rows:
        cells = [str(v) if isinstance(v, int) else f"{v:+.12e}" for v in r]
        if marked or nan_without_twin:
            twins = [cells[j] for j in cols]
            cells = ["nan" if j in cols else c for j, c in enumerate(cells)]
            if marked:
                cells += twins
        text.append(",".join(cells))
    (post / "surface-Q-corrected.csv").write_text("\n".join(text) + "\n")
    return post


class NotAppliedTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="csro-"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_main(self, post):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = CSRO.main([str(post), str(post)])
        return rc, buf.getvalue()

    def test_valid_tree_passes(self):
        rc, text = self.run_main(write_tree(self.tmp / "warn"))
        self.assertEqual(rc, 0)
        self.assertEqual(text, "Raw outputs are unchanged and corrected surface output is valid.\n")

    def test_marked_tree_reads_not_applied(self):
        post = write_tree(self.tmp / "mark", marked=True)
        rc, text = self.run_main(post)
        self.assertEqual(rc, CSRO.NOT_APPLIED_EXIT)
        self.assertNotEqual(rc, 0)
        self.assertTrue(text.startswith("NOT-APPLIED: 'E_norm corrected (J)' in corrected row 1, interface 1.0"), text)
        self.assertNotIn("valid", text)
        with self.assertRaises(RuntimeError) as caught:
            CSRO.check(post, post, 1.0e-10)
        self.assertIsInstance(caught.exception, CSRO.CorrectionNotApplied)

    def test_nan_without_twin_is_invalid(self):
        post = write_tree(self.tmp / "nan", nan_without_twin=True)
        with self.assertRaises(RuntimeError) as caught:
            CSRO.check(post, post, 1.0e-10)
        self.assertNotIsInstance(caught.exception, CSRO.CorrectionNotApplied)
        self.assertIn("Invalid corrected normalization energy", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
