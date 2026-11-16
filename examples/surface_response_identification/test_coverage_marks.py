#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for coverage_marks.py (decisions 1908 / 1914): the explicit mark columns are matched EXACTLY (never by prefix), the
MarkNotApplied signals of a table header / a tree, and read_cell's rule (mark 1 NOT-APPLIED whatever the token; mark 0 and a finite
number read; any other unmarked token or a missing mark column a Refusal; a Warn tree refuses a text token only)."""

import json
import math
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_marks as cm  # noqa: E402


class HeaderTest(unittest.TestCase):
    def test_mark_headers_exact(self):
        for h in ("coverage mark E_norm", "coverage mark E_surf", "coverage mark E_elec", "coverage mark measurement",
                  "coverage mark[3]", "coverage mark[i][12]"):
            self.assertTrue(cm.is_mark_header(h), h)
        for h in ("coverage mark E_surf2", "coverage mark", "coverage mark[3] (J)", "coverage mark[i][2]x", " coverage mark[3]",
                  "coverage marks[3]", "coverage mark[i]", "E_surf corrected (J)"):
            self.assertFalse(cm.is_mark_header(h), h)
        self.assertEqual(cm.MARK_INTERFACE.format(7), "coverage mark[7]")
        self.assertEqual(cm.MARK_MATRIX_ENTRY.format(2), "coverage mark[i][2]")

    def test_header_signals(self):
        self.assertEqual(cm.header_signals(["i", "E_elec raw (J)"]), [])
        self.assertEqual(cm.header_signals(["E_elec corrected (J)", "E_elec corrected (J)" + cm.PARTIAL_COLUMN_SUFFIX]), ["twin column"])
        self.assertEqual(cm.header_signals(["not-applied uncovered raw share of E_norm"]), ["exposure column"])
        self.assertEqual(cm.header_signals(["coverage mark[2]"]), ["mark column"])

    def test_policies_and_tree_signals(self):
        self.assertEqual(cm.config_policy({"Solver": {"SurfaceResponseCorrection": {"Coverage": {"Policy": "Warn"}}}}), "Warn")
        self.assertEqual(cm.config_policy({"Solver": {"Electrostatic": {"ResponseCorrection": {"Coverage": {"Policy": "MarkNotApplied"}}}}}),
                         "MarkNotApplied")
        self.assertIsNone(cm.config_policy({}))
        self.assertEqual(cm.palace_policy({"SurfaceResponse": {"Coverage": {"Policy": "MarkNotApplied"}}}), "MarkNotApplied")
        tmp = tempfile.mkdtemp()
        try:
            os.makedirs(os.path.join(tmp, "iteration1"))
            with open(os.path.join(tmp, "a.csv"), "w") as f:
                f.write("# a comment\ni,E (J)\n1,2\n")
            self.assertEqual(cm.tree_signals(tmp), [])
            with open(os.path.join(tmp, "iteration1", "b.csv"), "w") as f:
                f.write("i,coverage mark E_elec\n1,0\n")
            self.assertEqual(cm.tree_signals(tmp), [f"mark column in {os.path.join('iteration1', 'b.csv')}"])
            self.assertEqual(cm.tree_signals(tmp, palace={"SurfaceResponse": {"Coverage": {"Policy": "MarkNotApplied"}}})[-1],
                             "palace.json SurfaceResponse.Coverage.Policy MarkNotApplied")
            self.assertEqual(cm.tree_signals(tmp, config={"Solver": {"SurfaceResponseCorrection": {"Coverage": {"Policy": "MarkNotApplied"}}}})[-1],
                             "config_resolved.json Coverage.Policy MarkNotApplied")
            with open(os.path.join(tmp, "p.json"), "w") as f:
                json.dump({"x": 1}, f)
            self.assertEqual(cm.read_json_if_present(os.path.join(tmp, "p.json")), {"x": 1})
            self.assertIsNone(cm.read_json_if_present(os.path.join(tmp, "absent.json")))
        finally:
            shutil.rmtree(tmp)


class ReadCellTest(unittest.TestCase):
    def test_marked_tree(self):
        m = cm.MARK_E_SURF
        for token in ("nan", "", "NA", "null", "+1.0e-16", 1.0e-16, float("nan")):
            self.assertIsNone(cm.read_cell({"x": token, m: "1"}, "x", m, True, "w"), token)
            self.assertIsNone(cm.read_cell({"x": token, m: 1.0}, "x", m, True, "w"), token)
        self.assertEqual(cm.read_cell({"x": "+1.5e-16", m: "0"}, "x", m, True, "w"), 1.5e-16)
        self.assertEqual(cm.read_cell({"x": 2.5, m: 0.0}, "x", m, True, "w"), 2.5)
        for token in ("nan", "", "NA", "null", "+inf", "-inf", float("nan"), float("inf"), "word"):
            with self.assertRaises(cm.Refusal, msg=repr(token)) as caught:
                cm.read_cell({"x": token, m: "0"}, "x", m, True, "here")
            self.assertIn("not a finite number and not explicitly marked", str(caught.exception))
            self.assertTrue(str(caught.exception).startswith("REFUSED (decision 1908): here 'x'"))
        with self.assertRaises(cm.Refusal) as caught:
            cm.read_cell({"x": "1.0"}, "x", m, True, "w")
        self.assertIn("no explicit mark column", str(caught.exception))
        for bad in ("2", "", "yes", "nan"):
            with self.assertRaises(cm.Refusal, msg=bad):
                cm.read_cell({"x": "1.0", m: bad}, "x", m, True, "w")
        self.assertTrue(cm.mark_is_set({m: " 1 "}, m, "w"))
        self.assertFalse(cm.mark_is_set({m: "0"}, m, "w"))
        with self.assertRaises(cm.Refusal):
            cm.mark_is_set({}, m, "w")

    def test_warn_tree(self):
        m = cm.MARK_E_SURF
        self.assertTrue(math.isnan(cm.read_cell({"x": "nan"}, "x", m, False, "w")))
        self.assertTrue(math.isinf(cm.read_cell({"x": "+inf"}, "x", m, False, "w")))
        self.assertTrue(math.isnan(cm.read_cell({"x": float("nan")}, "x", m, False, "w")))
        self.assertEqual(cm.read_cell({"x": "+1.0e-16"}, "x", m, False, "w"), 1.0e-16)
        for token in ("", "NA", "null", "word"):
            with self.assertRaises(cm.Refusal, msg=token) as caught:
                cm.read_cell({"x": token}, "x", m, False, "w")
            self.assertIn("a text token in a corrected cell", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
