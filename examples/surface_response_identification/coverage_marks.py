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

ROUND 3 (decision 1934): (R-3) a table whose header names a mark, twin or exposure column TWICE is a REFUSAL in every reader
(refuse_duplicate_columns: two statements of one cell); (R-4) a mark token reads ONE way in every reader: the text '1' or '0'
(blank-padded as Palace prints it) exactly; anything else (1.0, true, '', a float a reader parsed) is a REFUSAL; (R-1) every
allowed report writer writes provenance(its own file) into its JSON: its name and sha256, the sha256 of this file and
"explicit_marks": true; (R-2) every report consumer calls accept_report: a report is read only when its provenance names a
writer on the consumer's allowed list (by sha256), or when it has NO provenance and its tree (the report's 'postpro') is on this
host and reads UNMARKED by the guard's rule (guard_rule_signals = marked_tree_guard_v2's table / config / palace.json signals);
else a ReportRefusal (exit REFUSAL_EXIT, one line). Provenance guards against the ACCIDENTAL use of a reader that is not on the
list (it writes no provenance); it is not a signature: a provenance object pasted into another JSON cannot be told from the
writer's own by its content.
Pure Python."""

import hashlib
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


def is_signal_column(header):
    """A mark, twin or exposure column (the columns that state or carry a cell's mark)."""
    return is_mark_header(header) or header.endswith(PARTIAL_COLUMN_SUFFIX) or header in EXPOSURE_HEADERS


def refuse_duplicate_columns(header, where):
    """Decision 1934 (R-3): a header (stripped names) naming a mark, twin or exposure column more than once is a Refusal (a
    dict of the row would keep only the last statement)."""
    counts = {}
    for name in header:
        if is_signal_column(name):
            counts[name] = counts.get(name, 0) + 1
    duplicates = sorted(name for name, n in counts.items() if n > 1)
    if duplicates:
        raise Refusal(where, duplicates[0], f"{counts[duplicates[0]]} columns",
                      "a duplicate mark / twin / exposure column name: two statements of one cell (decision 1934)")


def mark_token(row, mark_column, where):
    """Decision 1934 (R-4): the ONE reading of a mark token: the text '1' -> True, '0' -> False (blank padding stripped);
    any other token (a number a reader parsed, '1.0', 'true', '', a short row's None) is a Refusal."""
    mark = row[mark_column]
    if isinstance(mark, str):
        if mark.strip() == "1":
            return True
        if mark.strip() == "0":
            return False
    raise Refusal(where, mark_column, mark, "a mark column holds neither the token 0 nor the token 1")


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
    if mark_token(row, mark_column, where):
        return None
    value = token if isinstance(token, float) else (float(token) if is_number_token(token) else None)
    if value is None or not math.isfinite(value):
        raise Refusal(where, column, token, "not a finite number and not explicitly marked")
    return value


def mark_is_set(row, mark_column, where):
    """True / False for a mark column of `row` (the token 1 / 0); Refusal when absent or another token."""
    if mark_column not in row:
        raise Refusal(where, mark_column, "", f"no explicit mark column '{mark_column}' on a MarkNotApplied tree")
    return mark_token(row, mark_column, where)


# ---- decision 1934: report provenance (R-1) and the report consumers' acceptance rule (R-2) ----
PROVENANCE_KEY = "provenance"
REPORT_REFUSAL_TEXT = "report without an allowed provenance and its tree cannot be shown unmarked"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def provenance(writer_path):
    """The provenance object an allowed report writer puts into its JSON (R-1): the writer's file name and sha256, the sha256 of
    this module and the statement that the writer reads the explicit marks."""
    return {"writer": os.path.basename(writer_path), "writer_sha256": sha256_file(writer_path),
            "coverage_marks_sha256": sha256_file(__file__), "explicit_marks": True}


class ReportRefusal(Refusal):
    """A report a consumer may not read (R-2): printed as one line, exit REFUSAL_EXIT, no output file."""

    def __init__(self, report_path, why):
        self.where, self.column, self.token, self.why = report_path, PROVENANCE_KEY, "", why
        Exception.__init__(self, f"REFUSED (decision 1934): {report_path}: {REPORT_REFUSAL_TEXT} ({why})")


class GuardRuleRefuse(Exception):
    """The guard's fail-closed cases (a palace.json that does not parse or is truncated, an unknown Policy): the tree cannot be
    shown unmarked."""


# The guard's palace.json rule (P), = marked_tree_guard_v2.palace_signals (9f1ca5bb66a6a9fc) with its constants.
GUARD_MARK_KEYS = ("MarkedStatuses", "MarkedCells", "MarkedMeasurements", "MarkedWithoutCell", "Marked")
GUARD_KNOWN_POLICIES = ("Warn", MARK_NOT_APPLIED)
GUARD_PARSE_LIMIT = 64 << 20
GUARD_CHUNK = 64 << 20
GUARD_CARRY = 4096
GUARD_KEY_RE = re.compile(rb'"(MarkedStatuses|MarkedCells|MarkedMeasurements|MarkedWithoutCell|Marked)"\s*:')
GUARD_POLICY_RE = re.compile(rb'"Policy"\s*:\s*"MarkNotApplied"')
GUARD_COVERAGE_RE = re.compile(rb'"Coverage"\s*:\s*\{')


def _walk_marks(obj, found):
    if isinstance(obj, dict):
        if obj.get("Policy") == MARK_NOT_APPLIED:
            found.add("Policy MarkNotApplied")
        found.update(f"key {k}" for k in GUARD_MARK_KEYS if k in obj)
        for v in obj.values():
            _walk_marks(v, found)
    elif isinstance(obj, list):
        for v in obj:
            _walk_marks(v, found)


def _scan_chunks(path):
    with open(path, "rb") as handle:
        carry, offset = b"", 0
        while True:
            block = handle.read(GUARD_CHUNK)
            if not block:
                return
            buf = carry + block
            yield buf, offset
            keep = min(GUARD_CARRY, len(buf))
            carry, offset = buf[-keep:], offset + len(buf) - keep


def palace_signals(path):
    """(P) of the guard: -> (the signals found (sorted), SurfaceResponse.Coverage.Policy or None, how). Raises GuardRuleRefuse."""
    size = os.path.getsize(path)
    if size <= GUARD_PARSE_LIMIT:
        try:
            with open(path) as handle:
                doc = json.load(handle)
        except (ValueError, UnicodeDecodeError) as exc:
            raise GuardRuleRefuse(f"palace.json {path} does not parse ({type(exc).__name__}: {str(exc)[:80]})")
        found = set()
        _walk_marks(doc, found)
        sr = doc.get("SurfaceResponse") if isinstance(doc, dict) else None
        cov = sr.get("Coverage") if isinstance(sr, dict) else None
        if cov is not None and not isinstance(cov, dict):
            raise GuardRuleRefuse(f"palace.json {path}: SurfaceResponse.Coverage is not an object")
        policy = cov.get("Policy") if isinstance(cov, dict) else None
        return sorted(found), policy, f"json.load ({size} B)"
    with open(path, "rb") as handle:
        handle.seek(max(0, size - 64))
        tail = handle.read()
    if not re.search(rb"\n\}\s*\Z", tail):
        raise GuardRuleRefuse(f"palace.json {path} ({size} B) does not end with a column-0 '}}' (truncated)")
    found, policies, seen = set(), [], set()
    decoder = json.JSONDecoder()
    for buf, off in _scan_chunks(path):
        for m in GUARD_KEY_RE.finditer(buf):
            found.add(f"key {m.group(1).decode()}")
        if GUARD_POLICY_RE.search(buf):
            found.add("Policy MarkNotApplied")
        for m in GUARD_COVERAGE_RE.finditer(buf):
            at = off + m.start()
            if at in seen:
                continue
            seen.add(at)
            with open(path, "rb") as g:
                g.seek(at + m.end() - m.start() - 1)
                window = g.read(GUARD_CHUNK).decode("utf-8", "replace")
            try:
                obj, _ = decoder.raw_decode(window)
            except ValueError:
                raise GuardRuleRefuse(f"palace.json {path}: a Coverage object at byte {at} does not decode")
            if isinstance(obj, dict) and "Policy" in obj:
                policies.append(obj["Policy"])
    unknown = [p for p in policies if p not in GUARD_KNOWN_POLICIES]
    if unknown:
        raise GuardRuleRefuse(f"palace.json {path}: unknown Coverage.Policy {unknown[0]!r}")
    return sorted(found), (MARK_NOT_APPLIED if MARK_NOT_APPLIED in policies else (policies[0] if policies else None)), \
        f"full scan ({size} B; {len(seen)} Coverage object(s) decoded)"


def guard_rule_signals(postpro):
    """The MarkNotApplied signals of a tree by the guard's rule (marked_tree_guard_v2.decide): (T) the tables of postpro/ and its
    levels, (C) config_resolved.json's Coverage.Policy, (P) postpro/palace.json. [] = the tree reads UNMARKED (a tree without
    palace.json: T and C decide, the guard's written rule). Raises GuardRuleRefuse (fail-closed)."""
    out = tree_signals(postpro)
    path = os.path.join(postpro, "config_resolved.json")
    if os.path.isfile(path):
        try:
            with open(path) as handle:
                cfg = json.load(handle)
        except (ValueError, UnicodeDecodeError) as exc:
            raise GuardRuleRefuse(f"config_resolved.json {path} does not parse ({type(exc).__name__})")
        solver = (cfg.get("Solver") if isinstance(cfg, dict) else None) or {}
        for section in (solver.get("SurfaceResponseCorrection"), (solver.get("Electrostatic") or {}).get("ResponseCorrection")):
            if isinstance(section, dict) and isinstance(section.get("Coverage"), dict) and "Policy" in section["Coverage"]:
                policy = section["Coverage"]["Policy"]
                if policy == MARK_NOT_APPLIED:
                    out.append("config_resolved.json Coverage.Policy MarkNotApplied")
                elif policy not in GUARD_KNOWN_POLICIES:
                    raise GuardRuleRefuse(f"config_resolved.json: unknown Coverage.Policy {policy!r}")
                break
    pj = os.path.join(postpro, "palace.json")
    if os.path.isfile(pj):
        found, policy, _ = palace_signals(pj)
        if policy is not None and policy not in GUARD_KNOWN_POLICIES:
            raise GuardRuleRefuse(f"palace.json {pj}: unknown Coverage.Policy {policy!r}")
        out += [f"palace.json {s}" for s in found]
    return out


def accept_report(report, report_path, allowed_writers):
    """R-2: may a consumer read this report? `allowed_writers` = {writer sha256: the coverage_marks sha256 that writer requires}.
    -> the acceptance text; raises ReportRefusal. Accepted ONLY (a) with a provenance naming an allowed writer (and its module, and
    explicit_marks true), or (b) WITHOUT a provenance when the report's tree ('postpro') is on this host and reads unmarked by the
    guard's rule (an old report of an old tree stays readable)."""
    if not isinstance(report, dict):
        raise ReportRefusal(report_path, "the report is not a JSON object")
    if PROVENANCE_KEY in report:
        prov = report[PROVENANCE_KEY]
        if not isinstance(prov, dict):
            raise ReportRefusal(report_path, "its provenance is not an object")
        writer = prov.get("writer_sha256")
        if not isinstance(writer, str) or writer not in allowed_writers:
            raise ReportRefusal(report_path, f"its provenance names {prov.get('writer')!r} sha256 {str(writer)[:16]}: not on the allowed list")
        if prov.get("coverage_marks_sha256") != allowed_writers[writer] or prov.get("explicit_marks") is not True:
            raise ReportRefusal(report_path, f"its provenance does not state the explicit marks of the allowed writer {writer[:16]}")
        return f"provenance {prov.get('writer')} {writer[:16]} (allowed)"
    postpro = report.get("postpro")
    if not isinstance(postpro, str) or not os.path.isdir(postpro):
        raise ReportRefusal(report_path, f"no provenance; its tree {postpro} is not on this host")
    if os.path.isdir(os.path.join(postpro, "postpro")):  # the guard's reading of a run directory
        postpro = os.path.join(postpro, "postpro")
    try:
        signals = guard_rule_signals(postpro)
    except (GuardRuleRefuse, OSError) as exc:
        raise ReportRefusal(report_path, f"no provenance; its tree {postpro} cannot be read: {exc}")
    if signals:
        raise ReportRefusal(report_path, f"no provenance; its tree {postpro} is MARKED ({signals[0]})")
    return f"no provenance; its tree {postpro} reads unmarked (the guard's rule)"
