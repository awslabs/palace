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
    """Decisions 1934 (R-1) / 1942: the provenance object; accept_report accepts an allowed provenance only; every refusal."""

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
            self.assertTrue(str(caught.exception).startswith("REFUSED (decision 1942): r.json: "))

    def test_no_provenance_report_refused_whatever_its_tree(self):
        """Decision 1942 (R-2 (b) removed): a report WITHOUT provenance is refused whatever its tree: absent, on this host and
        unmarked (a postpro or a run directory; a Warn config / palace.json), marked (a mark / twin column, the config, palace.json),
        unreadable, or a report that borrowed a MARKED ladder (block_D.ladder.source) beside an unmarked 'postpro'."""
        unmarked = self.tree("warn", config={"Solver": {"Electrostatic": {"ResponseCorrection": {"Coverage": {"Policy": "Warn"}}}}},
                             palace='{"SurfaceResponse": {"Coverage": {"Policy": "Warn"}}}')
        marked = self.tree("mark", header="i,E_elec postprocessed fixed-trace (J),coverage mark E_elec\n")
        twin = self.tree("twin", header="i,E_elec postprocessed fixed-trace (J)" + cm.PARTIAL_COLUMN_SUFFIX + "\n")
        config = self.tree("config", config={"Solver": {"SurfaceResponseCorrection": {"Coverage": {"Policy": "MarkNotApplied"}}}})
        truncated = self.tree("truncated", palace='{"SurfaceResponse": {"Cov')
        forms = {"tree absent": {"postpro": os.path.join(self.tmp, "gone", "postpro")}, "no postpro key": {}, "postpro null": {"postpro": None},
                 "unmarked postpro on host": {"postpro": unmarked}, "unmarked run directory on host": {"postpro": os.path.dirname(unmarked)},
                 "mark column": {"postpro": marked}, "twin column": {"postpro": twin}, "config MarkNotApplied": {"postpro": config},
                 "palace.json truncated": {"postpro": truncated},
                 "borrowed marked ladder": {"postpro": unmarked, "block_D": {"ladder": {"source": marked, "borrowed": True}}}}
        for name, rep in forms.items():
            with self.assertRaises(cm.ReportRefusal, msg=name) as caught:
                cm.accept_report(rep, "r.json", self.allowed)
            self.assertEqual(str(caught.exception), f"REFUSED (decision 1942): r.json: {cm.REPORT_REFUSAL_TEXT} (no provenance)", name)
        with self.assertRaises(cm.ReportRefusal) as caught:
            cm.accept_report([], "r.json", self.allowed)
        self.assertEqual(str(caught.exception), f"REFUSED (decision 1942): r.json: {cm.REPORT_REFUSAL_TEXT} (the report is not a JSON object)")

    def test_accept_report_reads_no_file(self):
        """Decision 1944 (3c): accept_report decides from the report object alone: every file / directory access raises here, and
        the accepted and the refused forms give the same answers as without the patch."""
        import builtins
        from unittest import mock

        def forbidden(*args, **kwargs):
            raise AssertionError(f"accept_report touched the file system: {args[:1]}")

        good = {"postpro": self.tree("w"), cm.PROVENANCE_KEY: cm.provenance(self.writer)}
        bad = {"postpro": good["postpro"]}
        with mock.patch.object(builtins, "open", forbidden), mock.patch.object(os, "listdir", forbidden), \
                mock.patch.object(os, "stat", forbidden), mock.patch.object(os, "scandir", forbidden), \
                mock.patch.object(os.path, "isdir", forbidden), mock.patch.object(os.path, "isfile", forbidden), \
                mock.patch.object(os.path, "exists", forbidden):
            self.assertIn("(allowed)", cm.accept_report(good, "r.json", self.allowed))
            with self.assertRaises(cm.ReportRefusal):
                cm.accept_report(bad, "r.json", self.allowed)

    def test_the_tree_rule_of_round_3_is_gone(self):
        """Decision 1944 (3): the guard's tree rule (R-2 (b)'s only use) is no longer in this module."""
        for name in ("guard_rule_signals", "palace_signals", "GuardRuleRefuse", "GUARD_PARSE_LIMIT", "GUARD_MARK_KEYS"):
            self.assertFalse(hasattr(cm, name), name)

if __name__ == "__main__":
    unittest.main()
