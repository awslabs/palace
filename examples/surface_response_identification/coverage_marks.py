#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""coverage_marks.py -- the explicit per-cell marks of Coverage.Policy MarkNotApplied, read by the grading readers.

THE RULE (decisions 1908 / 1914; one module, imported by grade_cells.py, coverage_report.py and
examples/cpw2d/check_surface_response_output.py): a reader never infers a mark from a token. Under MarkNotApplied every
corrected table of Palace appends, last, integer 0 / 1 MARK COLUMNS (surfaceresponsecoverage.hpp; docs/src/reference.md),
each the writer's own mark of the cells it covers; their headers are matched EXACTLY here:

  surface-Q-corrected.csv (Maxwell)   MARK_E_NORM: the row's E_norm corrected / corrected fixed-flux / self-consistent
                                      MARK_E_SURF: the row interface's E / p / Q corrected, fixed-flux, self-consistent
  surface-Q-corrected.csv (ES)        MARK_E_ELEC: the source's E_elec postprocessed fixed-trace / fixed-flux / corrected
                                      MARK_INTERFACE.format(k): interface k's E / p / Q postprocessed fixed-trace / ff / corrected
  surface-response-confidence.csv     MARK_MEASUREMENT: the row's frequency / magnetic cells
  terminal-C / Cm / Cinv-corrected    MARK_MATRIX_ENTRY.format(j) on row i: entry (i, j), fixed-trace and corrected

A TREE is MarkNotApplied when ANY signal says so (tree_signals): a table header with a twin column (`... partial (uncovered
part raw)`), an exposure column (EXPOSURE_HEADERS) or a mark column; config_resolved.json's Coverage.Policy; palace.json's
SurfaceResponse.Coverage.Policy. On such a tree a corrected / normalization cell a reader turns into a reading is:
  mark 1                      -> NOT-APPLIED (no grade, band, number, sum or step; never its twin);
  mark 0 and a finite number  -> read as printed (the Warn twin's value);
  mark 0 and any other token  -> REFUSAL (nan, inf, '', NA, null, a word);
  its mark column absent      -> REFUSAL ('no explicit mark': a policy binary before the mark columns, or a re-writer).
On a Warn / pre-policy tree every reading stays frozen, except that a TEXT token (one float() refuses) in a read corrected
cell is a REFUSAL (decision 1914 (2)); nan / inf there keep the reader's frozen reading.
A REFUSAL raises Refusal; every reader prints its one line and exits REFUSAL_EXIT before it opens any output file.
Pure Python."""

import json
import math
import os
import re

PARTIAL_COLUMN_SUFFIX = " partial (uncovered part raw)"
MARK_E_NORM = "coverage mark E_norm"
MARK_E_SURF = "coverage mark E_surf"
MARK_E_ELEC = "coverage mark E_elec"
MARK_INTERFACE = "coverage mark[{}]"
MARK_MEASUREMENT = "coverage mark measurement"
MARK_MATRIX_ENTRY = "coverage mark[i][{}]"
MARK_HEADERS = (MARK_E_NORM, MARK_E_SURF, MARK_E_ELEC, MARK_MEASUREMENT)
MARK_INDEXED_RE = re.compile(r"coverage mark\[(?:i\]\[)?\d+\]")
EXPOSURE_HEADERS = ("not-applied uncovered raw share of E_norm", "not-applied uncovered raw share of E_elec")
MARK_NOT_APPLIED = "MarkNotApplied"
REFUSAL_EXIT = 5
NOT_APPLIED = "NOT-APPLIED"
# USER 1867 Q3 (keep + flag): the kept value of an unmarked Type on a row whose normalization mark is 1.
FLAG_Q3 = "FLAG Q3: normalization marked"


class Refusal(Exception):
    """A corrected / normalization cell that is not a finite number and is not explicitly marked, or a marked tree without
    the explicit mark column (decision 1908): the reader prints str(self) and exits REFUSAL_EXIT, writing nothing."""

    def __init__(self, where, column, token, why):
        self.where, self.column, self.token, self.why = where, column, token, why
        super().__init__(f"REFUSED (decision 1908): {where} '{column}' = '{token}': {why}")


def is_mark_header(header):
    """An explicit mark column (exact names; the indexed forms by a full match)."""
    return header in MARK_HEADERS or MARK_INDEXED_RE.fullmatch(header) is not None


def header_signals(header):
    """The MarkNotApplied signals of one table header (a list of stripped names): twin, exposure and mark columns."""
    out = []
    if any(h.endswith(PARTIAL_COLUMN_SUFFIX) for h in header):
        out.append("twin column")
    if any(h in EXPOSURE_HEADERS for h in header):
        out.append("exposure column")
    if any(is_mark_header(h) for h in header):
        out.append("mark column")
    return out


def csv_header(path):
    """The stripped header names of a Palace CSV (the first non-comment line), or None when the file is absent."""
    if not os.path.isfile(path):
        return None
    with open(path, "rb") as handle:
        for raw in handle:
            line = raw.decode("utf-8", "replace").rstrip("\r\n")
            if line.strip() and not line.lstrip().startswith("#"):
                return [c.strip() for c in line.split(",")]
    return []


def config_policy(config):
    """Coverage.Policy of a config dict (Solver.SurfaceResponseCorrection or Solver.Electrostatic.ResponseCorrection)."""
    solver = (config or {}).get("Solver") or {}
    for section in (solver.get("SurfaceResponseCorrection"), (solver.get("Electrostatic") or {}).get("ResponseCorrection")):
        policy = ((section or {}).get("Coverage") or {}).get("Policy") if isinstance(section, dict) else None
        if policy is not None:
            return policy
    return None


def palace_policy(palace):
    """SurfaceResponse.Coverage.Policy of a parsed palace.json (None: absent)."""
    coverage = ((palace or {}).get("SurfaceResponse") or {}).get("Coverage") or {}
    return coverage.get("Policy") if isinstance(coverage, dict) else None


def level_dirs(postpro):
    """postpro/ and its iterationNN/ levels."""
    out = [postpro]
    if os.path.isdir(postpro):
        for name in sorted(os.listdir(postpro)):
            if re.fullmatch(r"iteration\d+", name) and os.path.isdir(os.path.join(postpro, name)):
                out.append(os.path.join(postpro, name))
    return out


def tree_signals(postpro, config=None, palace=None):
    """Every MarkNotApplied signal of a postpro tree: the CSV headers of postpro/ and its levels, the config, palace.json."""
    out = []
    for d in level_dirs(postpro):
        for name in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            if name.endswith(".csv"):
                for s in header_signals(csv_header(os.path.join(d, name)) or []):
                    out.append(f"{s} in {os.path.relpath(os.path.join(d, name), postpro)}")
    if config_policy(config) == MARK_NOT_APPLIED:
        out.append("config_resolved.json Coverage.Policy MarkNotApplied")
    if palace_policy(palace) == MARK_NOT_APPLIED:
        out.append("palace.json SurfaceResponse.Coverage.Policy MarkNotApplied")
    return out


def read_json_if_present(path):
    if not path or not os.path.isfile(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def is_number_token(token):
    try:
        float(token)
    except (TypeError, ValueError):
        return False
    return True


def read_cell(row, column, mark_column, marked_tree, where):
    """The reading of ONE corrected / normalization cell of `row` ({stripped header: token}). -> a float, or None for a
    NOT-APPLIED (explicitly marked) cell. Raises Refusal. A Warn / pre-policy tree returns float(token) (nan / inf kept:
    the frozen reading) and refuses a text token only."""
    token = row[column]
    token = token.strip() if isinstance(token, str) else token
    if not marked_tree:
        if isinstance(token, float):
            return token
        if not is_number_token(token):
            raise Refusal(where, column, token, "a text token in a corrected cell")
        return float(token)
    if mark_column not in row:
        raise Refusal(where, column, token, f"no explicit mark column '{mark_column}' on a MarkNotApplied tree")
    mark = row[mark_column]
    mark = mark.strip() if isinstance(mark, str) else mark
    if mark in ("1", 1.0, 1):
        return None
    if mark not in ("0", 0.0, 0):
        raise Refusal(where, mark_column, mark, "a mark column holds neither 0 nor 1")
    value = token if isinstance(token, float) else (float(token) if is_number_token(token) else None)
    if value is None or not math.isfinite(value):
        raise Refusal(where, column, token, "not a finite number and not explicitly marked")
    return value


def mark_is_set(row, mark_column, where):
    """True / False for a mark column of `row` (0 / 1); Refusal when absent or another token."""
    if mark_column not in row:
        raise Refusal(where, mark_column, "", f"no explicit mark column '{mark_column}' on a MarkNotApplied tree")
    mark = row[mark_column]
    mark = mark.strip() if isinstance(mark, str) else mark
    if mark in ("1", 1.0, 1):
        return True
    if mark in ("0", 0.0, 0):
        return False
    raise Refusal(where, mark_column, mark, "a mark column holds neither 0 nor 1")
