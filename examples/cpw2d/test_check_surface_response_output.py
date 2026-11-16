#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""check_surface_response_output.py on synthetic Maxwell trees (decisions 1847 E-1 (c) / 1908 / 1914): a valid corrected tree passes
(exit 0); the MarkNotApplied form of the same tree (corrected cells nan, their values in the ' partial (uncovered part raw)' twins,
the explicit 'coverage mark E_norm' / 'coverage mark E_surf' columns 1) reads NOT-APPLIED (exit 3, never the 'valid' line;
CorrectionNotApplied is a RuntimeError, so callers of check() fail closed); a MarkNotApplied tree with every mark 0 is checked as its
Warn twin (the 1892 false positive: a by-design nan beside its twin is not a mark); an unmarked non-finite / text cell or a missing
mark column on a marked tree is a REFUSAL (exit 5); a text token on a Warn tree is a REFUSAL; a nan without any signal keeps the
'Invalid ...' error."""

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

cm = CSRO.cm
CORRECTED = ("E_norm corrected (J)", "E_surf corrected (J)", "p_surf corrected", "Q_surf corrected")


def write_tree(root, marked=False, nan_without_twin=False, explicit=True, marks=None, token="nan", edits=None):
    """One frequency, one excitation, interfaces 1 / 2. marked: the MarkNotApplied form (every corrected column + its twin; the
    cells of the rows `marks` selects ({interface: (E_norm mark, E_surf mark)}; default every row (1, 1)) print `token`, their
    values in the twins); `explicit` appends the mark columns (decision 1914; False = the form before it); `edits` =
    {(interface, column): token} applied last."""
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
    marks = marks if marks is not None else {1: (1, 1), 2: (1, 1)}
    out_header = list(header)
    if marked:
        out_header += [header[j] + CSRO.PARTIAL_COLUMN_SUFFIX for j in cols]
        if explicit:
            out_header += [cm.MARK_E_NORM, cm.MARK_E_SURF]
    text = [",".join(out_header)]
    for r in rows:
        k = r[2]
        cells = [str(v) if isinstance(v, int) else f"{v:+.12e}" for v in r]
        if nan_without_twin:
            cells = ["nan" if j in cols else c for j, c in enumerate(cells)]
        elif marked:
            norm_mark, surf_mark = marks.get(k, (0, 0))
            moved = {j for j in cols if (norm_mark if header[j].startswith("E_norm") else surf_mark)}
            twins = [cells[j] if j in moved else "nan" for j in cols]
            cells = [token if j in moved else c for j, c in enumerate(cells)] + twins
            if explicit:
                cells += [str(norm_mark), str(surf_mark)]
        for (iface, column), tok in (edits or {}).items():
            if iface == k:
                cells[out_header.index(column)] = tok
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
        for token in ("nan", "", "NA"):
            post = write_tree(self.tmp / f"mark-{token or 'empty'}", marked=True, token=token)
            rc, text = self.run_main(post)
            self.assertEqual(rc, CSRO.NOT_APPLIED_EXIT, token)
            self.assertTrue(text.startswith("NOT-APPLIED: corrected row 1, interface 1.0 is marked not applied"), text)
            self.assertIn("'coverage mark E_norm' = 1", text)
            self.assertNotIn("is valid", text)
            with self.assertRaises(RuntimeError) as caught:
                CSRO.check(post, post, 1.0e-10)
            self.assertIsInstance(caught.exception, CSRO.CorrectionNotApplied)

    def test_unmarked_marknotapplied_tree_reads_as_its_warn_twin(self):
        """Decision 1892's false positive gone: every mark 0 (a 0-mark MarkNotApplied tree whose twins exist on every row): valid."""
        rc, text = self.run_main(write_tree(self.tmp / "zero", marked=True, marks={}))
        self.assertEqual(rc, 0)
        self.assertEqual(text, "Raw outputs are unchanged and corrected surface output is valid.\n")

    def test_zero_mark_tree_with_a_by_design_infinity_beside_its_twin_is_valid(self):
        """Decision 1934 (the round-2 cs4 mutation's discriminating form): a 0-mark MarkNotApplied tree whose interface 2 has E = p = 0
        and Q = +inf BY DESIGN (Q_surf's own rule), its twin column present: valid, as its Warn twin (a twin beside a non-finite
        token is not a mark; only the mark column says so)."""
        edits = {(2, "E_surf corrected (J)"): "+0.000000000000e+00", (2, "p_surf corrected"): "+0.000000000000e+00",
                 (2, "Q_surf corrected"): "inf"}
        rc_w, text_w = self.run_main(write_tree(self.tmp / "inf-warn", edits=edits))
        rc, text = self.run_main(write_tree(self.tmp / "inf-zero", marked=True, marks={}, edits=edits))
        self.assertEqual((rc, text), (rc_w, text_w))
        self.assertEqual(rc, 0, text)

    def test_duplicate_mark_column_and_the_one_mark_token(self):
        """Decision 1934: a second 'coverage mark E_surf' column (R-3), and a mark token '1.0' (R-4): REFUSALS (exit 5)."""
        post = write_tree(self.tmp / "dup", marked=True, marks={1: (0, 1), 2: (0, 0)})
        path = post / "surface-Q-corrected.csv"
        lines = path.read_text().splitlines()
        path.write_text("\n".join([lines[0] + "," + cm.MARK_E_SURF] + [ln + ",0" for ln in lines[1:]]) + "\n")
        rc, text = self.run_main(post)
        self.assertEqual(rc, cm.REFUSAL_EXIT, text)
        self.assertIn("a duplicate mark / twin / exposure column name", text)
        post = write_tree(self.tmp / "tok", marked=True, marks={1: (0, 1), 2: (0, 0)}, edits={(1, cm.MARK_E_SURF): "1.0"})
        rc, text = self.run_main(post)
        self.assertEqual(rc, cm.REFUSAL_EXIT, text)
        self.assertIn("neither the token 0 nor the token 1", text)

    def test_refusals(self):
        cases = {
            "unmarked nan": write_tree(self.tmp / "r1", marked=True, marks={1: (0, 1)}, edits={(2, "E_surf corrected (J)"): "nan"}),
            "unmarked text": write_tree(self.tmp / "r2", marked=True, marks={}, edits={(1, "p_surf corrected"): "NA"}),
            "no mark column": write_tree(self.tmp / "r3", marked=True, explicit=False),
            "Warn text token": write_tree(self.tmp / "r4", edits={(1, "E_norm corrected (J)"): ""}),
        }
        for name, post in cases.items():
            rc, text = self.run_main(post)
            self.assertEqual(rc, cm.REFUSAL_EXIT, (name, text))
            self.assertTrue(text.startswith("REFUSED (decision 1908)"), (name, text))
        # a refusal outranks a NOT-APPLIED row on the same tree (never dropped)
        rc, text = self.run_main(write_tree(self.tmp / "r5", marked=True, marks={1: (1, 1), 2: (1, 0)},
                                            edits={(2, "p_surf corrected"): ""}))
        self.assertEqual(rc, cm.REFUSAL_EXIT, text)

    def test_nan_without_twin_is_invalid(self):
        post = write_tree(self.tmp / "nan", nan_without_twin=True)
        with self.assertRaises(RuntimeError) as caught:
            CSRO.check(post, post, 1.0e-10)
        self.assertNotIsInstance(caught.exception, CSRO.CorrectionNotApplied)
        self.assertIn("Invalid corrected normalization energy", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
