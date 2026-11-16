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
            self.assertIsNone(cm.read_cell({"x": token, m: "   1"}, "x", m, True, "w"), token)
        self.assertEqual(cm.read_cell({"x": "+1.5e-16", m: "0"}, "x", m, True, "w"), 1.5e-16)
        self.assertEqual(cm.read_cell({"x": 2.5, m: "0"}, "x", m, True, "w"), 2.5)
        for token in ("nan", "", "NA", "null", "+inf", "-inf", float("nan"), float("inf"), "word"):
            with self.assertRaises(cm.Refusal, msg=repr(token)) as caught:
                cm.read_cell({"x": token, m: "0"}, "x", m, True, "here")
            self.assertIn("not a finite number and not explicitly marked", str(caught.exception))
            self.assertTrue(str(caught.exception).startswith("REFUSED (decision 1908): here 'x'"))
        with self.assertRaises(cm.Refusal) as caught:
            cm.read_cell({"x": "1.0"}, "x", m, True, "w")
        self.assertIn("no explicit mark column", str(caught.exception))
        for bad in ("2", "", "yes", "nan", "1.0", "0.0", "true", "+1", 1.0, 0.0, 1, None):
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


class OneTokenOneReadingTest(unittest.TestCase):
    """Decision 1934 (R-4): a mark token reads ONE way: the text '1' / '0' (blank padding stripped); '1.0', a parsed float, 'true'
    are refusals in read_cell and mark_is_set alike (the review's N2: coverage_report read '1.0' as 1, grade_cells refused it)."""

    def test_mark_token(self):
        m = cm.MARK_INTERFACE.format(1)
        self.assertTrue(cm.mark_token({m: "                   1"}, m, "w"))
        self.assertFalse(cm.mark_token({m: "0 "}, m, "w"))
        for bad in ("1.0", "0.0", "1e0", "true", "", "01", 1.0, 1, 0.0, None):
            with self.assertRaises(cm.Refusal, msg=repr(bad)) as caught:
                cm.mark_token({m: bad}, m, "w")
            self.assertIn("neither the token 0 nor the token 1", str(caught.exception))
            with self.assertRaises(cm.Refusal, msg=repr(bad)):
                cm.mark_is_set({m: bad}, m, "w")
            with self.assertRaises(cm.Refusal, msg=repr(bad)):
                cm.read_cell({"x": "nan", m: bad}, "x", m, True, "w")


class DuplicateColumnTest(unittest.TestCase):
    """Decision 1934 (R-3): a header naming a mark, twin or exposure column twice is a Refusal (the review's N3); other duplicates
    and distinct names are not this rule's."""

    def test_duplicates_refused(self):
        twin = "E_surf postprocessed fixed-trace[1] (J)" + cm.PARTIAL_COLUMN_SUFFIX
        for name in ("coverage mark[1]", "coverage mark E_elec", "coverage mark[i][2]", twin, cm.EXPOSURE_HEADERS[1]):
            with self.assertRaises(cm.Refusal, msg=name) as caught:
                cm.refuse_duplicate_columns(["i", name, "x", name], "t.csv")
            self.assertIn("a duplicate mark / twin / exposure column name", str(caught.exception))
            self.assertIn(f"'{name}' = '2 columns'", str(caught.exception))
        cm.refuse_duplicate_columns(["i", "coverage mark[1]", "coverage mark[2]", twin, "x", "x"], "t.csv")


class ProvenanceTest(unittest.TestCase):
    """Decision 1934 (R-1 / R-2): the provenance object, and accept_report's two arms and every refusal."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cm-prov-")
        self.writer = os.path.join(self.tmp, "writer.py")
        with open(self.writer, "w") as f:
            f.write("# a writer\n")
        self.allowed = {cm.sha256_file(self.writer): cm.sha256_file(cm.__file__)}

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def tree(self, name, header="i,E_elec postprocessed fixed-trace (J)\n", palace=None, config=None):
        post = os.path.join(self.tmp, name, "postpro")
        os.makedirs(post)
        with open(os.path.join(post, "surface-Q-corrected.csv"), "w") as f:
            f.write(header + "1,1.0\n")
        if palace is not None:
            with open(os.path.join(post, "palace.json"), "w") as f:
                f.write(palace)
        if config is not None:
            with open(os.path.join(post, "config_resolved.json"), "w") as f:
                json.dump(config, f)
        return post

    def test_provenance_object(self):
        p = cm.provenance(self.writer)
        self.assertEqual(p, {"writer": "writer.py", "writer_sha256": cm.sha256_file(self.writer),
                             "coverage_marks_sha256": cm.sha256_file(cm.__file__), "explicit_marks": True})

    def test_allowed_provenance_accepted_whatever_the_tree(self):
        rep = {"postpro": "/nonexistent/postpro", cm.PROVENANCE_KEY: cm.provenance(self.writer)}
        self.assertIn("(allowed)", cm.accept_report(rep, "r.json", self.allowed))

    def test_provenance_not_allowed_refused(self):
        good = cm.provenance(self.writer)
        unmarked = self.tree("warn")
        for name, prov in (("other writer", dict(good, writer_sha256="0" * 64)), ("not an object", "coverage_report.py"),
                           ("module differs", dict(good, coverage_marks_sha256="1" * 64)),
                           ("explicit_marks false", dict(good, explicit_marks=False)), ("explicit_marks absent",
                           {k: v for k, v in good.items() if k != "explicit_marks"}), ("no sha", {"writer": "coverage_report.py"})):
            with self.assertRaises(cm.ReportRefusal, msg=name) as caught:
                cm.accept_report({"postpro": unmarked, cm.PROVENANCE_KEY: prov}, "r.json", self.allowed)
            self.assertIn(cm.REPORT_REFUSAL_TEXT, str(caught.exception))
            self.assertTrue(str(caught.exception).startswith("REFUSED (decision 1934): r.json: "))

    def test_no_provenance_unmarked_tree_on_host_accepted(self):
        for name, kw in (("plain", {}), ("warn config", {"config": {"Solver": {"Electrostatic": {"ResponseCorrection": {"Coverage": {"Policy": "Warn"}}}}}}),
                         ("warn palace", {"palace": '{"SurfaceResponse": {"Coverage": {"Policy": "Warn"}}}'})):
            post = self.tree(name.replace(" ", "-"), **kw)
            self.assertIn("reads unmarked", cm.accept_report({"postpro": post}, "r.json", self.allowed), name)
            self.assertIn("reads unmarked", cm.accept_report({"postpro": os.path.dirname(post)}, "r.json", self.allowed), name)

    def test_no_provenance_tree_absent_refused(self):
        for rep in ({"postpro": os.path.join(self.tmp, "gone", "postpro")}, {}, {"postpro": None}, []):
            with self.assertRaises(cm.ReportRefusal) as caught:
                cm.accept_report(rep, "r.json", self.allowed)
            self.assertIn(cm.REPORT_REFUSAL_TEXT, str(caught.exception))

    def test_no_provenance_marked_tree_refused(self):
        forms = {
            "mark column": {"header": "i,E_elec postprocessed fixed-trace (J),coverage mark E_elec\n"},
            "twin column": {"header": "i,E_elec postprocessed fixed-trace (J)" + cm.PARTIAL_COLUMN_SUFFIX + "\n"},
            "comment then mark": {"header": "# a comment line\ni,coverage mark[1]\n"},
            "config": {"config": {"Solver": {"SurfaceResponseCorrection": {"Coverage": {"Policy": "MarkNotApplied"}}}}},
            "palace key": {"palace": '{"Problem": {}, "X": {"MarkedCells": 3}}'},
            "palace policy": {"palace": '{"SurfaceResponse": {"Coverage": {"Policy": "MarkNotApplied"}}}'},
        }
        for name, kw in forms.items():
            post = self.tree(name.replace(" ", "-"), **kw)
            with self.assertRaises(cm.ReportRefusal, msg=name) as caught:
                cm.accept_report({"postpro": post}, "r.json", self.allowed)
            self.assertIn("is MARKED", str(caught.exception), name)

    def test_no_provenance_tree_fail_closed_refused(self):
        forms = {"palace truncated": {"palace": '{"SurfaceResponse": {"Cov'},
                 "palace unknown policy": {"palace": '{"SurfaceResponse": {"Coverage": {"Policy": "Abort"}}}'},
                 "config unknown policy": {"config": {"Solver": {"SurfaceResponseCorrection": {"Coverage": {"Policy": "Abort"}}}}}}
        for name, kw in forms.items():
            post = self.tree(name.replace(" ", "-"), **kw)
            with self.assertRaises(cm.ReportRefusal, msg=name) as caught:
                cm.accept_report({"postpro": post}, "r.json", self.allowed)
            self.assertIn("cannot be read", str(caught.exception), name)

    def test_large_palace_scan_and_its_tail_check(self):
        """The guard's large-file path (GUARD_PARSE_LIMIT lowered for the test): a scanned file that decodes but is truncated after
        its Coverage object is refused by the tail check alone (the round-2 g6 mutation's discriminating form)."""
        saved = cm.GUARD_PARSE_LIMIT
        cm.GUARD_PARSE_LIMIT = 16
        try:
            whole = '{\n "SurfaceResponse": {\n  "Coverage": {"Policy": "Warn"}\n },\n "Pad": "' + "x" * 64 + '"\n}\n'
            post = self.tree("large-ok", palace=whole)
            self.assertIn("reads unmarked", cm.accept_report({"postpro": post}, "r.json", self.allowed))
            post = self.tree("large-truncated", palace=whole[: whole.index('"Pad"')])
            with self.assertRaises(cm.ReportRefusal) as caught:
                cm.accept_report({"postpro": post}, "r.json", self.allowed)
            self.assertIn("does not end with a column-0", str(caught.exception))
        finally:
            cm.GUARD_PARSE_LIMIT = saved


if __name__ == "__main__":
    unittest.main()
