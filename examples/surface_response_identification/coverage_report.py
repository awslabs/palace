#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Coverage / confidence report for a corrected thin-mesh surface-response run (readiness STAGE D).

Reads ONLY what a Palace run already writes (its postpro directory: config_resolved.json,
palace.json, surface-Q-corrected / -edge, the surface-response-* CSVs, the iterationNN/ ladder)
plus the preflight manifest (surface-response-requirements.json) and writes
coverage-report.{json,md} with four blocks (DESIGN.md, readiness-20261029/coverage):

  A  COVERAGE: per interface Type the within-R raw surface energy W_t partitioned by what the
     library did with it (Exact / Interpolated / near-matched models vs the RAW classes: Missing,
     SCP holes, DomainBoundary, ConductorConsistency, GeometryFallbacks per class, the Maxwell
     loop remainder, the UNCLAIMED remainder as its own row: OI-COVERAGE-EXCLUSION-ENERGY),
     in ENERGY (the PatchEnergy raw partition; a PROXY by length when absent) and in LENGTH.
  B  the NAMED RESIDUAL TABLE: the terms the band does not calibrate, with their bounds of record.
  C  the PREDICTED ERROR BAND per Type from a calibration file (coverage_bands_v2.json, written by
     coverage_calibration.py) and ONE label HIGH / MEDIUM / LOW with its reasons.
  D  QUOTABILITY: the AMR ladder (ft per Type and its steps), the DEPTH rule, the closure spread.

Pure Python (no numpy); exit 0 when the report was written, 2 when a required input is missing.
"""

import argparse
import glob
import hashlib
import json
import math
import os
import re
import sys

TYPES = ("SA", "MS", "MA")
MODEL_STATUSES = ("Exact", "Interpolated", "QuantumNearMatched")
# The manifest counts LegacyContract requirements as corrected by exact-contract models.
STATUS_ROW = {"Exact": "Exact", "LegacyContract": "Exact", "Interpolated": "Interpolated",
              "QuantumNearMatched": "QuantumNearMatched"}
RAW_ROWS = ("Missing", "SCP", "DomainBoundary", "ConductorConsistency", "GeometryFallbacks",
            "LoopRemainder", "Unclaimed")
FALLBACK_CLASSES = ("OverlapExcluded", "UnlocatedPatch", "InterfaceMappingConflict", "KeyMismatch",
                    "NoMatchedFeature", "AllPatchesOutside", "NoPerimeter", "ContinuationUnresolved")
# DESIGN section 4: the class of a model = its requirement Topology, SpatialEdgeCluster split by
# EdgeCount band (the PLAN A.1.3 cost classes S / M / L use the same cuts); the jjsa key named.
SEC_BANDS = ((2, 3, "SpatialEdgeCluster-2-3"), (4, 8, "SpatialEdgeCluster-4-8"),
             (9, 20, "SpatialEdgeCluster-9-20"), (21, 10 ** 9, "SpatialEdgeCluster-gt20"))
JUNCTION_LEAD_KEYS = ("a16b8a6772d3",)
GROWN_COUPON_KEYS = ("a5e14f9da16f", "278e279fac2b", "85a2520aa7d3", "a16b8a6772d3")
# Label rules (PLAN section 4 item 3 + AMENDMENT 1 MINOR-7; DESIGN 1.3).
HIGH_EXACT_SHARE = 0.95
CLASS_MATERIAL_SHARE = 0.05
RAW_CLASS_HIGH_BAR = 0.01
RAW_CLASS_LOW_BAR = 0.10
CALIBRATED_MIN_WINDOWS = 3
# The DEPTH rule (DESIGN 1.2): the last ft step > 0.3 pt is not quotable; a step ratio >= 0.7 over the
# last three solves is not shrinking (DEPTH-UNBOUNDED); the Richardson tail is capped at 3 x the step.
QUOTABLE_LAST_STEP_PT = 0.3
DEPTH_RATIO_UNBOUNDED = 0.7
DEPTH_TAIL_CAP_FACTOR = 3.0
INDUCTIVE_EMAG_OVER_EELEC = 1.0
INDUCTIVE_LOOP_SHARE = 0.05
ENERGY_BOOKKEEPING_REL_TOL = 1e-9
LENGTH_BOOKKEEPING_TOL_UM = 1e-6
DEFAULT_BANDS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coverage_bands_v2.json")


class MissingInput(Exception):
    """A required input of the report is absent (exit code 2; the file is named)."""


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    if not os.path.isfile(path):
        raise MissingInput(path)
    with open(path) as f:
        return json.load(f)


def read_csv(path, required=False):
    """Palace CSV -> list of dict rows (stripped headers; '#' comment lines skipped; floats parsed)."""
    if not os.path.isfile(path):
        if required:
            raise MissingInput(path)
        return None
    rows = []
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cells = [c.strip() for c in line.split(",")]
            if header is None:
                header = cells
                continue
            row = {}
            for key, cell in zip(header, cells):
                row[key] = _parse_cell(cell)
            rows.append(row)
    return rows


def _parse_cell(cell):
    try:
        return float(cell)
    except ValueError:
        return cell


def indexed_columns(row, prefix):
    """{interface index: value} for columns like 'E_surf raw[3] (J)' / 'raw energy[3] (J)' / 'share[3]'."""
    out = {}
    pat = re.compile(r"^" + re.escape(prefix) + r"\[(\d+)\]")
    for key, value in row.items():
        m = pat.match(key)
        if m:
            out[int(m.group(1))] = value
    return out


def is_number(x):
    return isinstance(x, float) and not math.isnan(x)


def fmt_pct(x, digits=2):
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{100.0 * x:+.{digits}f} %"


def fmt_pt(x, digits=2):
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:+.{digits}f} pt"


# --------------------------------------------------------------------------------------------------
# The run (one Palace postpro directory, ES or Maxwell)
# --------------------------------------------------------------------------------------------------


def load_interfaces(config):
    """Dielectric[] -> {index: {"Type", "R_um", "LocalizeEdgeEnergy"}} (R = EdgeDistances.back())."""
    out = {}
    for entry in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", []):
        if "Index" not in entry or "Type" not in entry:
            continue
        distances = entry.get("EdgeDistances") or []
        out[int(entry["Index"])] = {
            "Type": entry["Type"],
            "R_um": float(distances[-1]) if distances else None,
            "LocalizeEdgeEnergy": bool(entry.get("LocalizeEdgeEnergy", False)),
        }
    return out


def solver_kind(config):
    problem = config.get("Problem", {})
    kind = problem.get("Type", "")
    if kind == "Electrostatic":
        return "ES"
    if kind == "Driven":
        return "Maxwell"
    return kind or "unknown"


def target_interfaces(config, kind):
    """The ResponseCorrection section: ES Solver.Electrostatic.ResponseCorrection; Maxwell
    Solver.SurfaceResponseCorrection. Returns (TargetInterfaces, the section)."""
    solver = config.get("Solver", {})
    if kind == "ES":
        rc = solver.get("Electrostatic", {}).get("ResponseCorrection", {})
    else:
        rc = solver.get("SurfaceResponseCorrection", {})
    return [int(i) for i in rc.get("TargetInterfaces", [])], rc


def load_surface_energies(postpro, kind):
    """Per interface index: E_surf raw / fixed-trace (ft) / fixed-flux (ff) / self-consistent (sc) in J,
    the closure columns; evaluation = the FIRST row (ES: the one excitation; Maxwell: the first f / exc)."""
    rows = read_csv(os.path.join(postpro, "surface-Q-corrected.csv"), required=True)
    if not rows:
        raise MissingInput(os.path.join(postpro, "surface-Q-corrected.csv") + " (no rows)")
    out = {"interfaces": {}, "closure": {}}
    if kind == "ES":
        row = rows[0]
        raw = indexed_columns(row, "E_surf raw")
        ft = indexed_columns(row, "E_surf postprocessed fixed-trace")
        ff = indexed_columns(row, "E_surf postprocessed fixed-flux")
        sc = indexed_columns(row, "E_surf corrected")
        for i in raw:
            out["interfaces"][i] = {"raw": raw[i], "ft": ft.get(i), "ff": ff.get(i), "sc": sc.get(i)}
        out["E_elec_raw"] = row.get("E_elec raw (J)")
        out["E_elec_ft"] = row.get("E_elec postprocessed fixed-trace (J)")
        out["closure"] = {
            "max trace closure spread": row.get("max trace closure spread"),
            "response-weighted local trace closure spread": row.get("response-weighted local trace closure spread"),
            "confidence pass": row.get("confidence pass"),
        }
    else:
        first = rows[0]
        key = (first.get("f (GHz)"), first.get("exc"))
        for row in rows:
            if (row.get("f (GHz)"), row.get("exc")) != key:
                continue
            i = int(row["interface"])
            out["interfaces"][i] = {"raw": row.get("E_surf raw (J)"), "ft": row.get("E_surf corrected (J)"),
                                    "ff": row.get("E_surf corrected fixed-flux (J)"),
                                    "sc": row.get("E_surf self-consistent (J)")}
        out["E_elec_raw"] = first.get("E_norm raw (J)")
        out["E_elec_ft"] = first.get("E_norm corrected (J)")
        out["evaluation_key"] = {"f_GHz": key[0], "exc": key[1]}
        conf = read_csv(os.path.join(postpro, "surface-response-confidence.csv"))
        if conf:
            c = conf[0]
            out["closure"] = {
                "max trace closure spread": c.get("max trace closure spread"),
                "response-weighted local trace closure spread": c.get("response-weighted local trace closure spread"),
                "confidence pass": c.get("confidence pass"),
                "loop residual": c.get("loop residual"),
                "response-weighted loop residual": c.get("response-weighted loop residual"),
            }
    return out


def load_edge_energies(postpro):
    """surface-Q-edge.csv -> {interface: E_out (J)} at the library radius (the first evaluation)."""
    rows = read_csv(os.path.join(postpro, "surface-Q-edge.csv"), required=True)
    out = {}
    key = None
    for row in rows:
        k = (row.get("i", row.get("f (GHz)")), row.get("exc"))
        if key is None:
            key = k
        if k != key:
            break
        out[int(row["interface"])] = {"E_out": row.get("E_out (J)"), "R_m": row.get("R (m)")}
    return out


def within_r_by_type(interfaces, targets, energies, edges):
    """W_t = sum over the Type's target interfaces of (E_surf raw[i] - E_out[i]); plus the ft / ff / sc sums."""
    out = {}
    for i in targets:
        if i not in interfaces or i not in energies["interfaces"]:
            continue
        t = interfaces[i]["Type"]
        e = energies["interfaces"][i]
        e_out = edges.get(i, {}).get("E_out", 0.0) or 0.0
        acc = out.setdefault(t, {"W": 0.0, "E_raw": 0.0, "E_out": 0.0, "ft": 0.0, "ff": 0.0, "sc": 0.0,
                                 "interfaces": []})
        acc["W"] += (e["raw"] or 0.0) - e_out
        acc["E_raw"] += e["raw"] or 0.0
        acc["E_out"] += e_out
        for k in ("ft", "ff", "sc"):
            acc[k] = (acc[k] + e[k]) if (is_number(e.get(k)) and acc[k] is not None) else None
        acc["interfaces"].append(i)
    return out


def load_catalog(palace):
    """palace.json SurfaceResponse.ModelCatalog -> {Index: {"Name", "Topology", "PatchCount"}}."""
    sr = palace.get("SurfaceResponse", {})
    out = {}
    for m in sr.get("ModelCatalog", []):
        out[int(m["Index"])] = {"Name": m.get("Name"), "Topology": m.get("Topology"),
                                "PatchCount": m.get("PatchCount"), "PatchWeight": m.get("PatchWeight")}
    return out


# --------------------------------------------------------------------------------------------------
# The preflight manifest (surface-response-requirements.json): classes and the LENGTH partition
# --------------------------------------------------------------------------------------------------


def requirement_class(req):
    """DESIGN section 4: Topology, SpatialEdgeCluster split by EdgeCount band, the jjsa key named."""
    topology = req.get("Topology") or "Unknown"
    key = (req.get("Hash") or "")[:12]
    if topology == "SpatialEdgeCluster":
        if key in JUNCTION_LEAD_KEYS or any(k in (s.get("Name") or "") for s in req.get("SelectedModels", [])
                                            for k in JUNCTION_LEAD_KEYS):
            return "JunctionLead"
        count = (req.get("Signature") or {}).get("EdgeCount")
        if isinstance(count, (int, float)):
            for lo, hi, name in SEC_BANDS:
                if lo <= count <= hi:
                    return name
        return "SpatialEdgeCluster-unknown"
    return topology


def requirement_types(req, all_types):
    types = sorted({i.get("Type") for i in req.get("Interfaces", []) if i.get("Type")})
    return types or list(all_types)


def load_manifest(path):
    """-> {"models": {name: {Status, Class, Topology, Grown, Unboxable, Hash}}, "requirements": [...], raw manifest}."""
    manifest = read_json(path)
    models = {}
    for req in manifest.get("Requirements", []):
        cls = requirement_class(req)
        grown = bool(((req.get("SpatialSupport") or {}).get("CornerContainment") or {}).get("Growth")) or \
            (req.get("Hash") or "")[:12] in GROWN_COUPON_KEYS
        for sel in req.get("SelectedModels", []):
            name = sel.get("Name")
            if not name:
                continue
            models.setdefault(name, {"Status": req.get("Status"), "Class": cls, "Topology": req.get("Topology"),
                                     "Grown": grown, "Hash": (req.get("Hash") or "")[:12],
                                     "Unboxable": bool((req.get("Signature") or {}).get("Unboxable"))})
    by_topology = {}
    for req in manifest.get("Requirements", []):
        if req.get("Status") not in STATUS_ROW:
            continue
        topo = req.get("Topology")
        length = float(req.get("TotalEdgeLength") or 0.0)
        if topo not in by_topology or length > by_topology[topo]["Length"]:
            by_topology[topo] = {"Status": req.get("Status"), "Class": requirement_class(req), "Length": length,
                                 "Hash": (req.get("Hash") or "")[:12]}
    return {"models": models, "by_topology": by_topology, "manifest": manifest, "path": path, "sha256": sha256_file(path)}


def length_partition(loaded, types, discovery=None):
    """Block A LENGTH column per Type (um): the model rows by Status x class, the RAW rows by class, the
    bookkeeping vs Identification.Totals.AssignedLength."""
    manifest = loaded["manifest"]
    unbuildable = set()
    if discovery:
        for key in discovery.get("UnbuildableRequirements", []) or []:
            unbuildable.add((key.get("Hash") if isinstance(key, dict) else str(key))[:12])
    out = {}
    for t in types:
        out[t] = {"models": {}, "classes": {}, "raw": {}, "missing_sub": {"Unboxable": 0.0, "Unbuildable": 0.0,
                                                                            "buildable": 0.0}, "grown": 0.0}
    for req in manifest.get("Requirements", []):
        length = float(req.get("TotalEdgeLength") or 0.0)
        status = req.get("Status")
        cls = requirement_class(req)
        key = (req.get("Hash") or "")[:12]
        grown = bool(((req.get("SpatialSupport") or {}).get("CornerContainment") or {}).get("Growth")) or \
            key in GROWN_COUPON_KEYS
        for t in requirement_types(req, types):
            if t not in out:
                continue
            if status in STATUS_ROW:
                row = STATUS_ROW[status]
                out[t]["models"][row] = out[t]["models"].get(row, 0.0) + length
                out[t]["classes"][cls] = out[t]["classes"].get(cls, 0.0) + length
                if grown:
                    out[t]["grown"] += length
            elif status == "Missing":
                out[t]["raw"]["Missing"] = out[t]["raw"].get("Missing", 0.0) + length
                sub = "Unboxable" if (req.get("Signature") or {}).get("Unboxable") else \
                    ("Unbuildable" if key in unbuildable else "buildable")
                out[t]["missing_sub"][sub] += length
            elif status == "DomainBoundary":
                out[t]["raw"]["DomainBoundary"] = out[t]["raw"].get("DomainBoundary", 0.0) + length
    summary = manifest.get("Summary", {})
    scp_outside = float((summary.get("SpatialClusterPriority") or {}).get("OutsideLength") or 0.0)
    db_raw = float((summary.get("DomainBoundary") or {}).get("RawPortionLength") or 0.0)
    fallbacks = (manifest.get("Identification", {}).get("Diagnostics", {}) or {}).get("GeometryFallbacks") or {}
    for t in types:
        if scp_outside:
            out[t]["raw"]["SCP"] = scp_outside
        if db_raw and "DomainBoundary" not in out[t]["raw"]:
            out[t]["raw"]["DomainBoundary"] = db_raw
        for entry in fallbacks.get("Entries", []) or []:
            if entry.get("Interfaces") and t not in entry["Interfaces"]:
                continue
            name = "GeometryFallbacks:" + str(entry.get("Class"))
            out[t]["raw"][name] = out[t]["raw"].get(name, 0.0) + float(entry.get("Length") or 0.0)
        total = sum(out[t]["models"].values()) + sum(out[t]["raw"].values())
        out[t]["total"] = total
        for part in ("models", "classes", "raw"):
            out[t][part + "_share"] = {k: (v / total if total else 0.0) for k, v in out[t][part].items()}
        out[t]["raw_share_total"] = sum(out[t]["raw"].values()) / total if total else 0.0
        out[t]["grown_share"] = out[t]["grown"] / total if total else 0.0
    assigned = (manifest.get("Identification", {}).get("Totals") or {}).get("AssignedLength")
    status_sum = sum(float(r.get("TotalEdgeLength") or 0.0) for r in manifest.get("Requirements", []))
    bookkeeping = None if assigned is None else status_sum - float(assigned)
    exclusions = [{"Class": e.get("Class"), "Count": e.get("Count"), "Length": e.get("Length")}
                  for e in manifest.get("Identification", {}).get("Exclusions", []) or []]
    return {"types": out, "assigned_length": assigned, "requirements_length": status_sum,
            "bookkeeping_miss_um": bookkeeping, "exclusions": exclusions,
            "fallback_locations": [{"Class": e.get("Class"), "What": e.get("What"), "Site": e.get("Site"),
                                    "ExampleOrigin": e.get("ExampleOrigin"), "Length": e.get("Length")}
                                   for e in fallbacks.get("Entries", []) or []][:20],
            "uncovered_by_type": (summary.get("Uncovered") or {}).get("ByType"),
            "counts": summary.get("Counts")}


# --------------------------------------------------------------------------------------------------
# Block A ENERGY: the PatchEnergy raw partition by model -> Status / class, the RAW class CSVs
# --------------------------------------------------------------------------------------------------

RAW_CSV_FILES = {
    "Missing": ("surface-response-uncovered-energy.csv", "uncovered raw energy"),
    "DomainBoundary": ("surface-response-domain-boundary-energy.csv", "domainboundary raw energy"),
    "SCP": ("surface-response-spatial-cluster-priority-energy.csv", "spatialclusterpriority raw energy"),
    "ConductorConsistency": ("surface-response-conductor-consistency-energy.csv", "conductorconsistency raw energy"),
}


def _eval0(rows):
    return [r for r in (rows or []) if r.get("evaluation", 0.0) == 0.0]


def _first_eval_rows(rows):
    """The evaluation-0 rows of the FIRST (f, exc) key (Maxwell has one block per excitation)."""
    rows = _eval0(rows)
    if not rows:
        return rows
    key = (rows[0].get("f (GHz)"), rows[0].get("exc"))
    return [r for r in rows if (r.get("f (GHz)"), r.get("exc")) == key]


def _sum_by_type(index_values, interfaces, targets):
    out = {}
    for i, v in index_values.items():
        if i in targets and i in interfaces and is_number(v):
            t = interfaces[i]["Type"]
            out[t] = out.get(t, 0.0) + v
    return out


def _accumulate(dst, add):
    for t, v in add.items():
        dst[t] = dst.get(t, 0.0) + v


def raw_class_energies(postpro, interfaces, targets, kind):
    """The recorded RAW classes' within-R raw energy per Type: {row: {Type: J}} plus their locations."""
    rows_out = {}
    locations = {}
    for row_name, (fname, prefix) in RAW_CSV_FILES.items():
        rows = read_csv(os.path.join(postpro, fname))
        if rows is None:
            continue
        acc = {}
        for r in _first_eval_rows(rows):
            _accumulate(acc, _sum_by_type(indexed_columns(r, prefix), interfaces, targets))
            if r.get("portions", 0) or r.get("length (m)", 0):
                locations.setdefault(row_name, []).append({"type": r.get("type"), "portions": r.get("portions"),
                                                           "length_m": r.get("length (m)")})
        rows_out[row_name] = acc
    fallbacks = read_csv(os.path.join(postpro, "surface-response-geometry-fallbacks.csv"))
    if fallbacks is not None:
        for r in fallbacks:
            name = "GeometryFallbacks:" + str(r.get("class"))
            t = r.get("type")
            acc = rows_out.setdefault(name, {})
            # columns: raw energy (J) = the fallback portions' own within-R raw energy; within-R energy (J) = the
            # Type's whole within-R energy (the denominator of 'share')
            if t in TYPES and is_number(r.get("raw energy (J)")):
                acc[t] = acc.get(t, 0.0) + r["raw energy (J)"]
                locations.setdefault(name, []).append({"type": t, "count": r.get("count"),
                                                       "length_m": r.get("length (m)")})
    if kind == "Maxwell":
        loops = read_csv(os.path.join(postpro, "surface-response-loop-remainder-energy.csv"))
        if loops is not None:
            acc = {}
            for r in _first_eval_rows(loops):
                _accumulate(acc, _sum_by_type(indexed_columns(r, "remainder energy"), interfaces, targets))
            rows_out["LoopRemainder"] = acc
    return rows_out, locations


def model_energies(postpro, interfaces, targets):
    """surface-response-model-energy.csv evaluation 0 -> {model index: {Type: fabricated energy (J)}}."""
    rows = read_csv(os.path.join(postpro, "surface-response-model-energy.csv"), required=True)
    out = {}
    for r in _first_eval_rows(rows):
        idx = int(r["model"])
        out[idx] = _sum_by_type(indexed_columns(r, "fabricated surface energy"), interfaces, targets)
    return out


def raw_partition(postpro, interfaces, targets):
    """surface-response-model-raw-energy.csv (PatchEnergy) -> ({label: {Type: J}}, {label: [model indices]},
    {label: weight length (m)}) for evaluation 0; None when the run has no raw partition (the PROXY path)."""
    rows = read_csv(os.path.join(postpro, "surface-response-model-raw-energy.csv"))
    if rows is None:
        return None
    labels = {}
    for r in _first_eval_rows(rows):
        labels[str(r["label"])] = _sum_by_type(indexed_columns(r, "raw energy"), interfaces, targets)
    groups = read_csv(os.path.join(postpro, "surface-response-model-raw-groups.csv")) or []
    members = {}
    weights = {}
    for g in groups:
        label = str(g["label"])
        members.setdefault(label, []).append(int(g["model"]))
        weights[(label, int(g["model"]))] = g.get("weight length (m)", 0.0)
    for label in labels:
        if label != "unclaimed" and label not in members:
            members[label] = [int(m[1:]) for m in label.split("+") if m.startswith("m")]
    return {"labels": labels, "members": members, "weights": weights}


# A solve-time catalog name absent from the manifest is the FIRST-ORDER CURVATURE TERM of a straight-like feature (the
# kappa0.1 node models of the curvature family; surfaceresponseoperator.cpp FindFirstOrderNodes: "described by its
# straight model (its anchor)", an additive term on the anchor's patches). Decision 1052 (the ANCHOR rule): such a
# model is credited to the Status / class of the manifest-known model it shares a raw-partition label with (its
# anchor). A model with no anchor in any label of the run falls back to the catalog Topology -> the LONGEST manifest
# requirement of that topology (the first topology present wins); every mapping is RECORDED (block_A.anchor_mapped_models
# with the energy share it carries). The catalog exports no requirement hash (a candidate export: ModelCatalog[].RequirementHash).
# Decision 1053 (the TIE-BREAK, order-independent over EVERY label the model appears in): prefer the STRAIGHT-LIKE known
# member (isolated edge / strip / gap: the families FindFirstOrderNodes anchors); among several, the member sharing the
# most raw energy with the model (summed over their shared labels and Types); then the lowest model index.
TOPOLOGY_FALLBACK = {"curved edge": ("CurvedEdge", "IsolatedEdge"), "isolated edge": ("IsolatedEdge",),
                     "convex corner": ("ConvexCorner",), "concave corner": ("ConcaveCorner",),
                     "parallel-edge cluster": ("ParallelEdgeCluster",), "spatial edge cluster": ("SpatialEdgeCluster",),
                     "same-conductor strip": ("SameConductorStrip",), "same-conductor gap": ("SameConductorGap",),
                     "different-conductor gap": ("DifferentConductorGap",)}
STRAIGHT_LIKE_TOPOLOGIES = frozenset({"isolatededge", "sameconductorstrip", "sameconductorgap", "differentconductorgap"})


def is_straight_like(topology):
    """True for the straight families the C++ anchors a curvature node to (isolated edge / strip / gap), whichever
    spelling the source uses (the manifest's IsolatedEdge, the catalog's 'isolated edge' / 'same-conductor strip')."""
    return re.sub(r"[^a-z]", "", str(topology or "").lower()) in STRAIGHT_LIKE_TOPOLOGIES


def anchor_models(partition, catalog, manifest_models):
    """{manifest-absent model index: the index of the manifest-known model it shares a raw label with} (the ANCHOR
    rule, decision 1052). Candidates = the known members of EVERY label the model appears in, ranked by the decision-1053
    tie-break: straight-like first, then the most raw energy shared with the model, then the lowest index (independent
    of the CSV row order)."""
    anchors = {}
    if partition is None:
        return anchors
    candidates = {}  # absent model -> {known member: shared raw energy (J)}
    for label, members in partition["members"].items():
        known = [m for m in members if catalog.get(m, {}).get("Name") in manifest_models]
        if not known:
            continue
        energy = sum(partition["labels"].get(label, {}).values())
        for m in members:
            if m in known:
                continue
            shared = candidates.setdefault(m, {})
            for k in known:
                shared[k] = shared.get(k, 0.0) + energy
    for m, shared in candidates.items():
        def rank(k):
            topology = manifest_models[catalog[k]["Name"]].get("Topology") or catalog[k].get("Topology")
            return (not is_straight_like(topology), -shared[k], k)
        anchors[m] = min(shared, key=rank)
    return anchors


def classify_model(idx, catalog, manifest_models, by_topology=None, mapped=None, anchors=None):
    """-> (status row, class, requirement info) of a catalog model via its manifest requirement (by name; else via its
    ANCHOR; else by topology fallback); a name with no requirement at all -> ("Unknown", Topology, None)."""
    entry = catalog.get(idx, {})
    info = manifest_models.get(entry.get("Name"))
    if info is None:
        topo = str(entry.get("Topology") or "")
        anchor = (anchors or {}).get(idx)
        if anchor is not None:
            info = manifest_models.get(catalog.get(anchor, {}).get("Name"))
            if mapped is not None:
                mapped.append({"model": idx, "name": entry.get("Name"), "topology": topo, "via": "anchor",
                               "anchor_model": anchor, "anchor_name": catalog.get(anchor, {}).get("Name"),
                               "mapped_to": info["Topology"], "status": info["Status"], "class": info["Class"]})
        elif by_topology:
            for candidate in TOPOLOGY_FALLBACK.get(topo, (topo,)):
                if candidate in by_topology:
                    info = by_topology[candidate]
                    if mapped is not None:
                        mapped.append({"model": idx, "name": entry.get("Name"), "topology": topo, "via": "topology",
                                       "anchor_model": None, "anchor_name": None, "mapped_to": candidate,
                                       "status": info["Status"], "class": info["Class"]})
                    break
    if info is None:
        return "Unknown", str(entry.get("Topology") or "Unknown"), None
    return STATUS_ROW.get(info["Status"], str(info["Status"])), info["Class"], info


def energy_partition(postpro, interfaces, targets, kind, catalog, manifest_models, W, lengths, by_topology=None):
    """Block A ENERGY column per Type. ENERGY mode when the raw partition exists; else the PROXY
    (model rows = W_t x the LENGTH share; the recorded RAW rows' own energies; Unclaimed = the remainder)."""
    partition = raw_partition(postpro, interfaces, targets)
    fab = model_energies(postpro, interfaces, targets)
    raw_rows, locations = raw_class_energies(postpro, interfaces, targets, kind)
    mode = "ENERGY" if partition is not None else "PROXY"
    anchors = anchor_models(partition, catalog, manifest_models)
    out = {"mode": mode, "types": {}, "locations": locations, "anchor_mapped_models": []}
    first_type = next(iter(W), None)
    for t, w in W.items():
        W_t = w["W"]
        acc = {"W": W_t, "E_out": w["E_out"], "models": {}, "classes": {}, "raw": {}, "unclaimed_label": None,
               "fab_total": 0.0, "raw_models_total": 0.0, "g": None, "grown": 0.0, "mapped_models": 0.0}
        fab_t = sum(v.get(t, 0.0) for v in fab.values())
        acc["fab_total"] = fab_t
        if partition is not None:
            for label, by_type in partition["labels"].items():
                e = by_type.get(t, 0.0)
                if label == "unclaimed":
                    acc["unclaimed_label"] = e
                    continue
                members = partition["members"].get(label, [])
                weights = [partition["weights"].get((label, m), 0.0) for m in members]
                wsum = sum(weights) or float(len(members))
                for m, wgt in zip(members, weights):
                    share = (wgt / wsum) if wsum else 1.0 / max(len(members), 1)
                    status, cls, info = classify_model(m, catalog, manifest_models, by_topology,
                                                       out["anchor_mapped_models"] if t == first_type else None, anchors)
                    acc["models"][status] = acc["models"].get(status, 0.0) + e * share
                    acc["classes"][cls] = acc["classes"].get(cls, 0.0) + e * share
                    acc["raw_models_total"] += e * share
                    if (info or {}).get("Grown"):
                        acc["grown"] += e * share
                    if catalog.get(m, {}).get("Name") not in manifest_models:
                        acc["mapped_models"] += e * share
            for row_name, by_type in raw_rows.items():
                acc["raw"][row_name] = by_type.get(t, 0.0)
            recorded_raw = sum(acc["raw"].values())
            acc["raw"]["Unclaimed"] = (acc["unclaimed_label"] or 0.0) - recorded_raw
            if acc["raw_models_total"] > 0.0:
                acc["g"] = (fab_t - acc["raw_models_total"]) / acc["raw_models_total"]
            total = sum(acc["models"].values()) + sum(acc["raw"].values())
            acc["bookkeeping_miss_rel"] = (total - W_t) / W_t if W_t else 0.0
        else:
            L = lengths["types"].get(t, {})
            for status, share in L.get("models_share", {}).items():
                acc["models"][status] = W_t * share
            for cls, share in L.get("classes_share", {}).items():
                acc["classes"][cls] = W_t * share
            for row_name, by_type in raw_rows.items():
                acc["raw"][row_name] = by_type.get(t, 0.0)
            for row_name, share in L.get("raw_share", {}).items():
                if row_name not in acc["raw"]:
                    acc["raw"][row_name] = W_t * share
            acc["raw"]["Unclaimed"] = W_t - sum(acc["models"].values()) - sum(acc["raw"].values())
            acc["raw_models_total"] = sum(acc["models"].values())
            acc["grown"] = W_t * L.get("grown_share", 0.0)
            acc["bookkeeping_miss_rel"] = 0.0
        for part in ("models", "classes", "raw"):
            acc[part + "_share"] = {k: (v / W_t if W_t else 0.0) for k, v in acc[part].items()}
        acc["raw_share_total"] = sum(acc["raw"].values()) / W_t if W_t else 0.0
        acc["exact_share"] = acc["models_share"].get("Exact", 0.0)
        acc["grown_share"] = acc["grown"] / W_t if W_t else 0.0
        acc["unclaimed_share"] = acc["raw_share"].get("Unclaimed", 0.0)
        acc["mapped_models_share"] = acc["mapped_models"] / W_t if W_t else 0.0
        out["types"][t] = acc
    return out


# --------------------------------------------------------------------------------------------------
# Block D: the AMR ladder (ft per Type per solve), the steps, the DEPTH rule
# --------------------------------------------------------------------------------------------------


def iteration_dirs(base):
    """iterationNN/ (zero-padded or not: 'iteration1', 'iteration02') in numeric order."""
    found = []
    for d in glob.glob(os.path.join(base, "iteration*")):
        m = re.match(r"iteration(\d+)$", os.path.basename(d))
        if m and os.path.isdir(d):
            found.append((int(m.group(1)), d))
    return [d for _, d in sorted(found)]


def ladder(postpro, interfaces, targets, kind, ladder_dir=None):
    """Every iterationNN/ + the root (the FINAL solve) -> per Type the ft fraction ft / E_elec_ft and its
    step in pt relative to the final; the DOF ladder from palace.json where present."""
    base = ladder_dir or postpro
    dirs = iteration_dirs(base)
    solves = []
    for d in dirs + [base]:
        if not os.path.isfile(os.path.join(d, "surface-Q-corrected.csv")):
            continue
        try:
            e = load_surface_energies(d, kind)
        except MissingInput:
            continue
        edges = load_edge_energies(d) if os.path.isfile(os.path.join(d, "surface-Q-edge.csv")) else {}
        W = within_r_by_type(interfaces, targets, e, edges)
        denom = e.get("E_elec_ft")
        entry = {"dir": os.path.basename(d) if d != base else "final", "E_elec_ft": denom, "ft": {}, "raw": {}}
        for t, w in W.items():
            entry["ft"][t] = (w["ft"] / denom) if (is_number(denom) and denom and w["ft"] is not None) else None
            entry["raw"][t] = (w["E_raw"] / denom) if (is_number(denom) and denom) else None
        solves.append(entry)
    dofs = None
    pj = os.path.join(base, "palace.json")
    if os.path.isfile(pj):
        amr = read_json(pj).get("AdaptiveMeshRefinement")
        if isinstance(amr, list):
            dofs = [a.get("DegreesOfFreedom") for a in amr if isinstance(a, dict)] or None
        elif isinstance(amr, dict):
            dofs = amr.get("DegreesOfFreedom") if isinstance(amr.get("DegreesOfFreedom"), list) else None
    return {"source": base, "borrowed": ladder_dir is not None, "solves": solves, "dofs": dofs,
            "types": depth_rule(solves)}


def depth_rule(solves):
    """Per Type: the last ft step (pt of the final ft), the step ratio over the last three solves, the
    quotability verdict, and the Richardson tail (DESIGN 1.2 DEPTH)."""
    out = {}
    types = set()
    for s in solves:
        types.update(s["ft"].keys())
    for t in sorted(types):
        values = [s["ft"].get(t) for s in solves if s["ft"].get(t) is not None]
        rec = {"solves": len(values), "last_step_pt": None, "step_ratio": None, "quotable": None,
               "tail_pt": 0.0, "status": "ladder: not in this run" if len(values) < 2 else "ok",
               "steps_pt": []}
        if len(values) >= 2 and values[-1]:
            steps = [100.0 * (values[i] - values[i - 1]) / values[-1] for i in range(1, len(values))]
            rec["steps_pt"] = steps
            rec["last_step_pt"] = steps[-1]
            if len(steps) >= 2 and steps[-2] != 0.0:
                rec["step_ratio"] = abs(steps[-1]) / abs(steps[-2])
            rec["quotable"] = abs(steps[-1]) <= QUOTABLE_LAST_STEP_PT
            ratio = rec["step_ratio"]
            if rec["quotable"]:
                rec["status"] = "quotable"
            elif ratio is not None and ratio >= DEPTH_RATIO_UNBOUNDED:
                rec["status"] = "DEPTH-UNBOUNDED"
                rec["tail_pt"] = math.copysign(DEPTH_TAIL_CAP_FACTOR * abs(steps[-1]), steps[-1])
            else:
                rec["status"] = "DEPTH"
                tail = abs(steps[-1]) / (1.0 - ratio) if ratio is not None else abs(steps[-1])
                tail = min(tail, DEPTH_TAIL_CAP_FACTOR * abs(steps[-1]))
                rec["tail_pt"] = math.copysign(tail, steps[-1])
        out[t] = rec
    return out


# --------------------------------------------------------------------------------------------------
# Block B: the NAMED RESIDUAL TABLE (frozen data; a row prints when its trigger fires; decision 475)
# --------------------------------------------------------------------------------------------------

NAMED_RESIDUALS = (
    {"term": "grown-coupon (T5) excess", "decisions": "988 / 990",
     "bound": "+4.30 % of the grown box = +0.48 pt of the device MA_sharp on the transmon; scaled by the grown coupons' MA within-R share",
     "trigger": "a model of a grown key (a5e14f9da16f / 278e279fac2b / 85a2520aa7d3 / a16b8a6772d3) or a requirement with CornerContainment.Growth",
     "band_use": "ADDED one-sided to the MA upper edge", "key": "grown"},
    {"term": "the cut convention (reference side)", "decisions": "988",
     "bound": "direct - B = -2.71 pt of the box", "trigger": "always listed; fires on no thin run",
     "band_use": "REPORTED, never added", "key": "cut"},
    {"term": "the MA reading convention (raw vs sharp)", "decisions": "PLAN item 2; 994",
     "bound": "MA raw vs MA_sharp: the protocol band [-2, +9] %; the transmon raw +6.38 vs sharp +5.91",
     "trigger": "an MA Type whose config has no MA shell form", "band_use": "the MA band labelled 'MA raw' (LOW on MA)",
     "key": "ma_raw"},
    {"term": "p3-REFERENCE / u_p", "decisions": "973 / 978 / 994",
     "bound": "u_p 2.012 % is an ESTIMATE, not a bound; the transmon MA_sharp +5.91 % [+4.95, +6.57]",
     "trigger": "--reference-order 3 (the MA_sharp reading graded against a p3 reference; DESIGN 5.1)",
     "band_use": "the band WIDENED by u_p on the p3-graded Type (MA) only; labelled p3-REFERENCE",
     "key": "p3"},
    {"term": "the C form terms", "decisions": "921 / 936 / 938 / 949",
     "bound": "R-50 class +0.31 / +0.32 / +0.37 / +0.98 % (C4 / O2 / S5 / S6, Thomson); junction class -0.14..-0.26 % (Faraday EMF)",
     "trigger": "a Maxwell run: LumpedPort R = 50 ohm (R-50 class) or a junction-driven port (L_J / C_J)",
     "band_use": "the C column labelled by class; C never in the label", "key": "c_form"},
    {"term": "Maxwell INDUCTIVE", "decisions": "920 / 923 / 950 / 952",
     "bound": "C1: E_mag / E_elec 11.5, the [2] identity +4.25 % NOT MET (973 / 995)",
     "trigger": "E_mag / E_elec > 1 OR the max loop share > 0.05", "band_use": "labelled INDUCTIVE (OPERATOR); confidence capped at MEDIUM",
     "key": "inductive"},
    {"term": "DEPTH (ES / Maxwell)", "decisions": "1016 / 1024 / 1045 / 1047",
     "bound": "S7-rc SA -6.02 F (stopped at the 11th solve); S1b-rc SA -6.48 / MS -5.19 F (the AMR does not refine at junction footprints)",
     "trigger": "the run's own ladder (ES or Maxwell): the last ft step > 0.3 pt or the step ratio >= 0.7",
     "band_use": "NOT QUOTABLE: the band extended by the Richardson tail (ratio < 0.7) or DEPTH-UNBOUNDED (LOW)", "key": "depth"},
    {"term": "Maxwell FORM", "decisions": "929 D / 1024",
     "bound": "C_eff* +1.42 F [FORM / OPERATOR] on S7-rc (the S6 class); OI-MX-PEMBED-SLIVER-CONDITIONING",
     "trigger": "a P-embedded R2a port (--port-form r2a-pembed, or the config's LumpedPort Length override)",
     "band_use": "the C column labelled FORM; no Type band change",
     "key": "form"},
    {"term": "the junction-footprint AMR", "decisions": "1045",
     "bound": "OI-MX-JUNCTION-FOOTPRINT-AMR: the AMR error ~70 % at the 3 junction footprints",
     "trigger": "a Maxwell run with a junction LumpedPort and an error indicator that does not fall between solves",
     "band_use": "listed; the DEPTH rule carries it", "key": "jj_amr"},
    {"term": "the depth tail (ES / Maxwell)", "decisions": "PLAN item 2; 951 / 953",
     "bound": "+-0.2 pt ft; the transmon's last steps +0.43 / the W-i tail +0.50; S7-rc last ft steps SA +0.46 / MS +0.70 / MA +0.21",
     "trigger": "every run with >= 2 AMR solves", "band_use": "the Richardson tail of the DEPTH row", "key": "tail"},
    {"term": "the isolated-edge MS reading", "decisions": "994",
     "bound": "S1 MS +1.76 after stretchfix; the transmon MS -1.73",
     "trigger": "Type MS with > 50 % of its within-R energy under IsolatedEdge models", "band_use": "listed as the MS bias context",
     "key": "ms_isolated"},
)

GROWN_EXCESS_PT_PER_SHARE = 4.30  # +4.30 % of the grown box's MA energy (988 / 990)
P3_UP_PCT = 2.0122280303430604
# DESIGN 5.1: the transmon of record grades MA_sharp alone against the 973 r2.5 p3 reference (SA / MS / C vs the thick
# ES best), so --reference-order 3 widens the MA band only (DESIGN 1.2: "WIDENED by u_p on that Type").
P3_GRADED_TYPES = ("MA",)
PORT_FORMS = ("r2a-pembed",)


def maxwell_flags(postpro, config, energies, port_form=None):
    """INDUCTIVE (E_mag / E_elec, the max loop share), the C class (R-50 / junction), FORM (a P-embedded R2a port:
    --port-form r2a-pembed or the config's LumpedPort Length override; DESIGN 1.2 / decision 929 item D)."""
    flags = {"inductive": False, "E_mag_over_E_elec": None, "max_loop_share": None, "c_class": None,
             "form": port_form in PORT_FORMS, "port_form": port_form, "junction_port": False}
    ports = config.get("Boundaries", {}).get("LumpedPort", []) or []
    excited = [p for p in ports if p.get("Excitation")] or ports
    for p in excited:
        if (p.get("R") or 0.0) == 50.0:
            flags["c_class"] = "R-50"
        if (p.get("L") or 0.0) != 0.0 or (p.get("C") or 0.0) != 0.0:
            flags["junction_port"] = True
            flags["c_class"] = flags["c_class"] or "junction"
        if p.get("Length") is not None:
            flags["form"] = True
    fit = read_csv(os.path.join(postpro, "surface-response-loop-fit.csv"))
    if fit:
        flags["max_loop_share"] = _eval0(fit)[0].get("max loop share") if _eval0(fit) else None
    rows = read_csv(os.path.join(postpro, "domain-E.csv"))
    if rows:
        r = rows[0]
        e_mag = next((v for k, v in r.items() if k.startswith("E_mag")), None)
        e_elec = next((v for k, v in r.items() if k.startswith("E_elec")), None)
        if is_number(e_mag) and is_number(e_elec) and e_elec:
            flags["E_mag_over_E_elec"] = e_mag / e_elec
    if (flags["E_mag_over_E_elec"] or 0.0) > INDUCTIVE_EMAG_OVER_EELEC or \
            (flags["max_loop_share"] or 0.0) > INDUCTIVE_LOOP_SHARE:
        flags["inductive"] = True
    return flags


def indicator_falls(postpro):
    """error-indicators.csv of every iterationNN/ + the root: does the Norm fall between solves?"""
    norms = []
    for d in iteration_dirs(postpro) + [postpro]:
        rows = read_csv(os.path.join(d, "error-indicators.csv"))
        if rows and is_number(rows[0].get("Norm")):
            norms.append(rows[0]["Norm"])
    if len(norms) < 2:
        return None
    return norms[-1] < norms[0]


def named_residuals(kind, energy, lengths, ladder_types, flags, reference_order, ma_shell, manifest_models):
    """-> list of rows {term, fired, detail}; the fired terms carry their band extensions per Type."""
    rows = []
    # the grown coupons' within-R share: ENERGY where the run has the raw partition, else LENGTH (the proxy)
    grown_share = {t: acc.get("grown_share", 0.0) for t, acc in energy["types"].items()}
    for r in NAMED_RESIDUALS:
        fired = False
        detail = ""
        ext = {}
        k = r["key"]
        if k == "grown":
            fired = any(v > 0.0 for v in grown_share.values()) or any(m.get("Grown") for m in manifest_models.values())
            if fired:
                ma = grown_share.get("MA", 0.0)
                ext = {"MA": {"upper_pt": GROWN_EXCESS_PT_PER_SHARE * ma}}
                detail = f"grown coupons = {fmt_pct(ma)} of the MA within-R energy ({energy['mode']}) -> +{GROWN_EXCESS_PT_PER_SHARE * ma:.2f} pt (one-sided)"
        elif k == "ma_raw":
            fired = ("MA" in energy["types"]) and not ma_shell
            detail = "MA shells not in the config: the MA band is the sharp one, the run reads MA raw" if fired else ""
        elif k == "p3":
            fired = reference_order == 3
            if fired:
                ext = {t: {"half_pt": P3_UP_PCT} for t in energy["types"] if t in P3_GRADED_TYPES}
                detail = f"u_p {P3_UP_PCT:.3f} % added to the half-width of the p3-graded Type ({' / '.join(ext) or 'none in this run'})"
        elif k == "c_form":
            fired = kind == "Maxwell" and flags.get("c_class") is not None
            detail = f"C class {flags.get('c_class')}" if fired else ""
        elif k == "inductive":
            fired = bool(flags.get("inductive"))
            detail = f"E_mag / E_elec {flags.get('E_mag_over_E_elec')}; max loop share {flags.get('max_loop_share')}" if fired else ""
        elif k == "depth":
            bad = {t: v for t, v in ladder_types.items() if v["status"] in ("DEPTH", "DEPTH-UNBOUNDED")}
            fired = bool(bad)
            ext = {t: {"tail_pt": v["tail_pt"], "status": v["status"]} for t, v in bad.items()}
            detail = "; ".join(f"{t} {v['status']} last step {fmt_pt(v['last_step_pt'])} ratio {v['step_ratio']}" for t, v in bad.items())
        elif k == "form":
            fired = kind == "Maxwell" and bool(flags.get("form"))
            detail = (f"P-embedded R2a port ({flags.get('port_form') or 'the config Length override'}): the C column is FORM"
                      if fired else "")
        elif k == "jj_amr":
            fired = kind == "Maxwell" and bool(flags.get("junction_port")) and flags.get("indicator_falls") is False
            detail = "the error indicator does not fall between solves" if fired else ""
        elif k == "tail":
            fired = any(v["solves"] >= 2 for v in ladder_types.values())
            detail = "; ".join(f"{t} last step {fmt_pt(v['last_step_pt'])}" for t, v in ladder_types.items() if v["solves"] >= 2)
        elif k == "ms_isolated":
            ms = energy["types"].get("MS", {})
            share = ms.get("classes_share", {}).get("IsolatedEdge", 0.0)
            fired = share > 0.5
            detail = f"IsolatedEdge {fmt_pct(share)} of the MS within-R energy" if fired else ""
        rows.append({"term": r["term"], "decisions": r["decisions"], "bound": r["bound"], "trigger": r["trigger"],
                     "band_use": r["band_use"], "fired": fired, "detail": detail, "extensions": ext})
    return rows


# --------------------------------------------------------------------------------------------------
# Block C: the predicted error band per Type from the calibration file, and the label
# --------------------------------------------------------------------------------------------------

CORNER_CLASSES = ("ConvexCorner", "ConcaveCorner")


def class_predictors(lengths_t):
    """The M2 levers (LENGTH shares): the SpatialEdgeCluster share and the corner share."""
    shares = lengths_t.get("classes_share", {}) if lengths_t else {}
    sec = sum(v for k, v in shares.items() if k.startswith("SpatialEdgeCluster") or k == "JunctionLead")
    corner = sum(v for k, v in shares.items() if k in CORNER_CLASSES)
    return {"sec_share": sec, "corner_share": corner}


def predict_bias(band_t, predictors):
    """b_t(x) for the selected model: M0 beta0; M1 beta0 + beta1 x last step; M2 beta0 + beta1 sec + beta2 corner."""
    coef = band_t.get("coefficients", {})
    model = band_t.get("model", "M0")
    b = coef.get("beta0", 0.0)
    if model == "M1":
        step = predictors.get("last_step_pt")
        b += coef.get("beta1", 0.0) * (step if step is not None else 0.0)
    elif model == "M2":
        b += coef.get("beta1", 0.0) * predictors.get("sec_share", 0.0) + coef.get("beta2", 0.0) * predictors.get("corner_share", 0.0)
    return b


def load_bands(path):
    bands = read_json(path)
    bands["_path"] = path
    bands["_sha256"] = sha256_file(path)
    return bands


def band_table(bands, family, energy, lengths, ladder_types, residuals, proxy):
    """Per Type: b +- h in pt with the gap term r_t |g_t| and the fired extensions; C beside (M0 only)."""
    fam = (bands.get("families") or {}).get(family) or {}
    types_cal = fam.get("types") or {}
    out = {}
    for t, acc in energy["types"].items():
        cal = types_cal.get(t)
        lad = ladder_types.get(t, {})
        pred = class_predictors(lengths["types"].get(t))
        pred["last_step_pt"] = lad.get("last_step_pt")
        row = {"calibrated": cal is not None, "predictors": pred}
        if cal is None:
            row["note"] = f"no calibration for {family} / {t}"
            out[t] = row
            continue
        b = predict_bias(cal, pred)
        H = cal.get("H_pt", 0.0) * cal.get("widened", 1.0)
        g = acc.get("g")
        g_source = "the run's own correction ratio"
        if g is None or proxy:
            g = cal.get("median_g", 0.0)
            g_source = "the calibration median (no raw partition)"
        r = acc.get("raw_share_total", 0.0)
        gap = 100.0 * r * abs(g)
        lo = b - H - gap
        hi = b + H + gap
        if proxy or acc.get("g") is None:
            pass  # two-sided: g is the calibration median, not the run's own (DESIGN 1.3)
        else:
            # one-sided in the direction of the THIN bias (PLAN section 4 item 3; decision 1052): the uncorrected share r
            # misses the library's correction g, so thin + library is LOW by r g when g > 0 and the discrepancy
            # (thin + library - thick) / thick can be MORE NEGATIVE: the band extends downward; g < 0 the mirror
            if g > 0:
                hi = b + H
            else:
                lo = b - H
        ext_notes = []
        for res in residuals:
            if not res["fired"]:
                continue
            e = res["extensions"].get(t)
            if not e:
                continue
            if "half_pt" in e:
                lo -= e["half_pt"]
                hi += e["half_pt"]
                ext_notes.append(f"{res['term']}: +-{e['half_pt']:.2f} pt")
            if "upper_pt" in e:
                hi += e["upper_pt"]
                ext_notes.append(f"{res['term']}: +{e['upper_pt']:.2f} pt (upper)")
            if "tail_pt" in e:
                # a positive step = the ft still RISING with depth = the thin result is LOW vs the converged
                # one, so the discrepancy (thin - thick) can be MORE NEGATIVE: the band extends downward
                if e["tail_pt"] >= 0:
                    lo -= e["tail_pt"]
                else:
                    hi -= e["tail_pt"]
                ext_notes.append(f"{res['term']}: {e['status']} tail {fmt_pt(e['tail_pt'])} (band extended {'down' if e['tail_pt'] >= 0 else 'up'})")
        row.update({"model": cal.get("model"), "coefficients": cal.get("coefficients"), "bias_pt": b, "H_pt": cal.get("H_pt"),
                    "widened": cal.get("widened", 1.0), "gap_pt": gap, "g": g, "g_source": g_source,
                    "raw_share": r, "low_pt": lo, "high_pt": hi, "extensions": ext_notes,
                    "loo_coverage": cal.get("loo_coverage"), "loo_coverage_unwidened": unwidened_loo(cal), "n": cal.get("n"),
                    "median_residual_pt": cal.get("median_residual_pt")})
        out[t] = row
    c = types_cal.get("C")
    if c:
        # DESIGN 5.5: C gets the M0 band (coverage_calibration.py fits C with M0 only); the band printed IS the band
        # validated: beta0 (the median) +- H x widened. A C entry of another model is a calibration-file error.
        if c.get("model") != "M0":
            raise ValueError(f"the C band must be M0 (DESIGN 5.5); the band file carries {c.get('model')}")
        beta0 = c.get("coefficients", {}).get("beta0", 0.0)
        half = c.get("H_pt", 0.0) * c.get("widened", 1.0)
        out["C"] = {"calibrated": True, "model": c.get("model"), "bias_pt": beta0, "H_pt": c.get("H_pt"),
                    "widened": c.get("widened", 1.0), "loo_coverage": c.get("loo_coverage"), "loo_coverage_unwidened": unwidened_loo(c),
                    "n": c.get("n"), "low_pt": beta0 - half, "high_pt": beta0 + half,
                    "note": "C is printed beside; it never enters the label (DESIGN 5.5)"}
    return out


def unwidened_loo(cal):
    """The selected model's LOO coverage BEFORE the post-hoc widening (the candidates' record), when widened."""
    if cal.get("widened", 1.0) <= 1.0:
        return None
    return ((cal.get("candidates") or {}).get(cal.get("model")) or {}).get("loo_coverage")


def label_run(energy, lengths, ladder_types, residuals, flags, bands_family, proxy, ma_reading):
    """HIGH / MEDIUM / LOW with the list of reasons (DESIGN 1.3; PLAN section 4 item 3 + MINOR-7)."""
    census = (bands_family or {}).get("class_census") or {}
    reasons_low = []
    reasons_not_high = []
    for t, acc in energy["types"].items():
        column = acc["classes_share"] if not proxy else lengths["types"].get(t, {}).get("classes_share", {})
        exact = acc["exact_share"] if not proxy else lengths["types"].get(t, {}).get("models_share", {}).get("Exact", 0.0)
        if exact < HIGH_EXACT_SHARE:
            reasons_not_high.append(f"{t}: Exact share {fmt_pct(exact)} < 95 %")
        for cls, share in column.items():
            if share > CLASS_MATERIAL_SHARE and census.get(cls, 0) < CALIBRATED_MIN_WINDOWS:
                reasons_low.append(f"{t}: UNCALIBRATED class {cls} {fmt_pct(share)} > 5 % (n_c = {census.get(cls, 0)})")
        for row_name, share in acc["raw_share"].items():
            if share > RAW_CLASS_LOW_BAR:
                reasons_low.append(f"{t}: RAW {row_name} {fmt_pct(share)} > 10 %")
            elif share > RAW_CLASS_HIGH_BAR:
                reasons_not_high.append(f"{t}: RAW {row_name} {fmt_pct(share)} > 1 %")
        if acc["raw_share_total"] > RAW_CLASS_LOW_BAR:
            reasons_low.append(f"{t}: RAW total {fmt_pct(acc['raw_share_total'])} > 10 %")
        lad = ladder_types.get(t, {})
        if lad.get("status") == "DEPTH-UNBOUNDED":
            reasons_low.append(f"{t}: DEPTH-UNBOUNDED (step ratio {lad.get('step_ratio')})")
        elif lad.get("status") == "DEPTH":
            reasons_not_high.append(f"{t}: not quotable (last step {fmt_pt(lad.get('last_step_pt'))})")
        elif lad.get("status") == "ladder: not in this run":
            reasons_not_high.append(f"{t}: ladder not in this run")
    for name in lengths.get("fallback_locations", []):
        if name.get("What") == "run":
            reasons_low.append(f"GeometryFallback What = run ({name.get('Class')}): 100 % raw")
    if ma_reading == "raw" and "MA" in energy["types"]:
        reasons_low.append("MA read raw where the band is the sharp one (the MA_sharp transformation not supplied)")
    if flags.get("inductive"):
        reasons_not_high.append("INDUCTIVE (OPERATOR, 950): confidence capped at MEDIUM")
    if proxy:
        reasons_not_high.append("ENERGY = PROXY (by length; PatchEnergy off): capped at MEDIUM")
    label = "LOW" if reasons_low else ("MEDIUM" if reasons_not_high else "HIGH")
    return {"label": label, "reasons_low": reasons_low, "reasons_not_high": reasons_not_high}


# --------------------------------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------------------------------


def build_report(postpro, manifest_path, bands=None, discovery_path=None, ladder_dir=None, reference_order=None,
                 ma_reading="raw", label=None, port_form=None):
    """The full report (dict) for one run. Raises MissingInput when a required input is absent."""
    config = read_json(os.path.join(postpro, "config_resolved.json"))
    palace = read_json(os.path.join(postpro, "palace.json"))
    kind = solver_kind(config)
    family = "ES" if kind == "ES" else "Maxwell"
    interfaces = load_interfaces(config)
    targets, rc = target_interfaces(config, kind)
    targets = targets or sorted(interfaces)
    types = [t for t in TYPES if any(interfaces[i]["Type"] == t for i in targets if i in interfaces)]
    untargeted = sorted(i for i in interfaces if i not in targets)
    energies = load_surface_energies(postpro, kind)
    edges = load_edge_energies(postpro)
    W = within_r_by_type(interfaces, targets, energies, edges)
    catalog = load_catalog(palace)
    loaded = load_manifest(manifest_path)
    discovery = read_json(discovery_path) if discovery_path else None
    lengths = length_partition(loaded, types, discovery)
    energy = energy_partition(postpro, interfaces, targets, kind, catalog, loaded["models"], W, lengths, loaded["by_topology"])
    proxy = energy["mode"] == "PROXY"
    lad = ladder(postpro, interfaces, targets, kind, ladder_dir)
    flags = maxwell_flags(postpro, config, energies, port_form) if kind == "Maxwell" else {}
    if kind == "Maxwell":
        flags["indicator_falls"] = indicator_falls(ladder_dir or postpro)
    residuals = named_residuals(kind, energy, lengths, lad["types"], flags, reference_order, ma_reading == "sharp",
                               loaded["models"])
    localize = any(v["LocalizeEdgeEnergy"] for v in interfaces.values())
    report = {
        "tool": "coverage_report.py", "label_name": label or os.path.basename(os.path.dirname(postpro.rstrip("/"))) or postpro,
        "postpro": os.path.abspath(postpro), "manifest": {"path": os.path.abspath(manifest_path), "sha256": loaded["sha256"]},
        "solver": kind, "family": family, "dofs": (palace.get("Problem") or {}).get("DegreesOfFreedom"),
        "interfaces": {str(i): v for i, v in interfaces.items()}, "targets": targets, "untargeted": untargeted,
        "types": types, "energy_mode": energy["mode"], "localize_edge_energy": localize,
        "W": {t: {"W_J": w["W"], "E_raw_J": w["E_raw"], "E_out_J": w["E_out"], "ft_J": w["ft"], "ff_J": w["ff"],
                  "interfaces": w["interfaces"]} for t, w in W.items()},
        "E_elec_raw_J": energies.get("E_elec_raw"), "E_elec_ft_J": energies.get("E_elec_ft"),
        "block_A": {"energy": energy["types"], "length": lengths["types"], "energy_locations": energy["locations"],
                    "anchor_mapped_models": energy["anchor_mapped_models"],
                    "mapped_models_share": {t: acc.get("mapped_models_share", 0.0) for t, acc in energy["types"].items()},
                    "length_bookkeeping_miss_um": lengths["bookkeeping_miss_um"], "assigned_length_um": lengths["assigned_length"],
                    "exclusions": lengths["exclusions"], "fallback_locations": lengths["fallback_locations"],
                    "uncovered_by_type": lengths["uncovered_by_type"], "counts": lengths["counts"]},
        "block_B": residuals, "block_D": {"ladder": lad, "closure": energies.get("closure"), "flags": flags,
                                          "ff_vs_ft": {t: ((w["ff"] - w["ft"]) / w["ft"] if (w["ft"] and w["ff"] is not None) else None)
                                                       for t, w in W.items()}},
        "ma_reading": ma_reading, "reference_order": reference_order, "port_form": port_form,
    }
    if bands is not None:
        fam = (bands.get("families") or {}).get(family) or {}
        report["block_C"] = band_table(bands, family, energy, lengths, lad["types"], residuals, proxy)
        if "C" in report["block_C"] and kind == "Maxwell":
            # DESIGN 1.2: the C column labelled by class (R-50 / junction) and FORM (a P-embedded R2a port)
            report["block_C"]["C"]["c_class"] = flags.get("c_class")
            report["block_C"]["C"]["form"] = bool(flags.get("form"))
        report["bands"] = {"path": bands["_path"], "sha256": bands["_sha256"], "version": bands.get("version"),
                           "fit_date": bands.get("fit_date"), "provenance": bands.get("provenance"),
                           "dev_set_caveat": bands.get("dev_set_caveat")}
        report["label"] = label_run(energy, lengths, lad["types"], residuals, flags, fam, proxy, ma_reading)
    else:
        report["block_C"] = None
        report["label"] = {"label": "n/a (no band file)", "reasons_low": [], "reasons_not_high": []}
    return report


def build_preflight_report(manifest_path, palace_path=None, discovery_path=None, label=None):
    """The preflight-only form (the full chips): block A LENGTH only; the READY criterion-1 reading by length
    (<= 5 % uncorrected within-R share per Type: Missing + fallback + SCP outside + DomainBoundary; PLAN A.8)."""
    loaded = load_manifest(manifest_path)
    discovery = read_json(discovery_path) if discovery_path else None
    manifest = loaded["manifest"]
    types = sorted({i.get("Type") for r in manifest.get("Requirements", []) for i in r.get("Interfaces", []) if i.get("Type")}) or list(TYPES)
    lengths = length_partition(loaded, types, discovery)
    ready = {}
    for t in types:
        share = lengths["types"][t]["raw_share_total"]
        ready[t] = {"uncorrected_share_by_length": share, "ready": share <= CLASS_MATERIAL_SHARE}
    report = {"tool": "coverage_report.py", "form": "preflight", "label_name": label or os.path.basename(os.path.dirname(manifest_path)),
              "manifest": {"path": os.path.abspath(manifest_path), "sha256": loaded["sha256"]}, "types": types,
              "block_A": {"length": lengths["types"], "length_bookkeeping_miss_um": lengths["bookkeeping_miss_um"],
                          "assigned_length_um": lengths["assigned_length"], "exclusions": lengths["exclusions"],
                          "fallback_locations": lengths["fallback_locations"], "uncovered_by_type": lengths["uncovered_by_type"],
                          "counts": lengths["counts"]},
              "block_B": "not available on a preflight", "block_C": "not available on a preflight (coverage only; no band)",
              "block_D": "not available on a preflight", "ready_criterion_1": ready,
              "label": {"label": "READY" if all(v["ready"] for v in ready.values()) else "NOT READY",
                        "reasons_low": [f"{t}: uncorrected {fmt_pct(v['uncorrected_share_by_length'])} by length > 5 %"
                                        for t, v in ready.items() if not v["ready"]], "reasons_not_high": []}}
    if palace_path and os.path.isfile(palace_path):
        pj = read_json(palace_path)
        report["dofs"] = (pj.get("Problem") or {}).get("DegreesOfFreedom")
        gf = ((pj.get("SurfaceResponse") or {}).get("Diagnostics") or {}).get("GeometryFallbacks")
        if gf:
            report["block_A"]["geometry_fallbacks_record"] = {"Classes": gf.get("Classes"), "Count": gf.get("Count"),
                                                              "RawShareByLength": gf.get("RawShareByLength")}
    return report


# --------------------------------------------------------------------------------------------------
# Output: coverage-report.{json,md} and the one-line summary
# --------------------------------------------------------------------------------------------------


def _row_names(energy_t, length_t):
    names = list(MODEL_STATUSES)
    for k in list(energy_t.get("models", {})) + list(length_t.get("models", {})):
        if k not in names:
            names.append(k)
    raw = []
    for k in list(energy_t.get("raw", {})) + list(length_t.get("raw", {})):
        if k not in raw:
            raw.append(k)
    if "Unclaimed" in raw:
        raw.remove("Unclaimed")
        raw.append("Unclaimed")
    return names, raw


def render_markdown(report):
    lines = [f"# COVERAGE REPORT — {report['label_name']}", ""]
    lab = report["label"]
    lines.append(f"**LABEL: {lab['label']}**")
    for r in lab.get("reasons_low", []):
        lines.append(f"- LOW: {r}")
    for r in lab.get("reasons_not_high", []):
        lines.append(f"- not HIGH: {r}")
    lines.append("")
    if report.get("form") == "preflight":
        lines += ["Preflight form: block A LENGTH only; blocks B / C / D not available on a preflight (coverage only; no band).", ""]
        lines.append("| Type | uncorrected by length | READY (<= 5 %) |")
        lines.append("|---|---|---|")
        for t, v in report["ready_criterion_1"].items():
            lines.append(f"| {t} | {fmt_pct(v['uncorrected_share_by_length'])} | {'yes' if v['ready'] else 'NO'} |")
        lines.append("")
    else:
        lines.append(f"Solver {report['solver']}; DOFs {report.get('dofs')}; ENERGY = {report['energy_mode']}"
                     + (" (by length; PatchEnergy off): the label is capped at MEDIUM" if report["energy_mode"] == "PROXY" else "")
                     + ("; OI-LOCALEDGE: uncovered energy may be low under AMR" if report.get("localize_edge_energy") else ""))
        if report.get("bands"):
            b = report["bands"]
            lines.append(f"Bands: {os.path.basename(b['path'])} version {b.get('version')} fit {b.get('fit_date')} sha256 {b['sha256'][:16]}…"
                         f" — {b.get('dev_set_caveat') or ''}")
        lines.append("")
    lines.append("## A. COVERAGE by Type (ENERGY | LENGTH)")
    A = report["block_A"]
    for t in report["types"]:
        e = (A.get("energy") or {}).get(t, {})
        L = A["length"].get(t, {})
        lines.append("")
        W = e.get("W")
        lines.append(f"### {t}" + (f" — W_t = {W:.4e} J (beyond-R E_out {e.get('E_out', 0.0):.3e} J)" if W is not None else "")
                     + f"; assigned length {L.get('total', 0.0):.2f} um")
        lines.append("| row | ENERGY share | LENGTH share |")
        lines.append("|---|---|---|")
        model_rows, raw_rows = _row_names(e, L)
        for k in model_rows:
            lines.append(f"| corrected by {k} models | {fmt_pct(e.get('models_share', {}).get(k, 0.0)) if e else 'n/a'} | {fmt_pct(L.get('models_share', {}).get(k, 0.0))} |")
        for k in raw_rows:
            lines.append(f"| RAW: {k} | {fmt_pct(e.get('raw_share', {}).get(k, 0.0)) if e else 'n/a'} | {fmt_pct(L.get('raw_share', {}).get(k, 0.0))} |")
        if e:
            lines.append(f"| RAW total | {fmt_pct(e.get('raw_share_total'))} | {fmt_pct(L.get('raw_share_total', 0.0))} |")
            miss = e.get("bookkeeping_miss_rel", 0.0)
            if abs(miss) > ENERGY_BOOKKEEPING_REL_TOL and report["energy_mode"] == "ENERGY":
                lines.append(f"| BOOKKEEPING MISS (energy) | {miss:+.3e} rel | |")
        sub = L.get("missing_sub", {})
        if any(sub.values()):
            lines.append(f"| Missing sub-classes (length) | | Unboxable {sub.get('Unboxable', 0):.2f} / Unbuildable {sub.get('Unbuildable', 0):.2f} / buildable {sub.get('buildable', 0):.2f} um |")
        lines.append("")
        lines.append("classes (ENERGY | LENGTH): " + "; ".join(
            f"{c} {fmt_pct(e.get('classes_share', {}).get(c, 0.0)) if e else 'n/a'} | {fmt_pct(L.get('classes_share', {}).get(c, 0.0))}"
            for c in sorted(set(list(L.get("classes_share", {})) + list((e or {}).get("classes_share", {}))))))
    if A.get("length_bookkeeping_miss_um") is not None and abs(A["length_bookkeeping_miss_um"]) > LENGTH_BOOKKEEPING_TOL_UM:
        lines.append(f"\nBOOKKEEPING MISS (length): requirements - AssignedLength = {A['length_bookkeeping_miss_um']:+.4f} um")
    lines.append("\nExclusions (length by class; no per-class energy export: OI-COVERAGE-EXCLUSION-ENERGY): "
                 + "; ".join(f"{x['Class']} {x['Count']} x {x['Length']:.1f} um" for x in A.get("exclusions", [])))
    if A.get("fallback_locations"):
        lines.append("GeometryFallbacks: " + "; ".join(f"{f['Class']} ({f['What']}) at {f['ExampleOrigin']} {f.get('Length', 0.0):.2f} um"
                                                      for f in A["fallback_locations"]))
    if A.get("anchor_mapped_models"):
        shares = A.get("mapped_models_share") or {}
        lines.append("Manifest-absent catalog models (the first-order curvature node models; decision 1052 ANCHOR rule): "
                     + "; ".join(f"m{m['model']} {m['name']} -> {m['status']} / {m['class']} via {m['via']}"
                                 + (f" (anchor m{m['anchor_model']} {m['anchor_name']})" if m['via'] == 'anchor' else f" ({m['mapped_to']})")
                                 for m in A["anchor_mapped_models"])
                     + "; energy share they carry: " + ", ".join(f"{t} {fmt_pct(s)}" for t, s in shares.items()))
    if report.get("form") == "preflight":
        return "\n".join(lines) + "\n"
    lines += ["", "## B. NAMED RESIDUALS"]
    lines.append("| term | decisions | bound of record | fired | detail |")
    lines.append("|---|---|---|---|---|")
    for r in report["block_B"]:
        lines.append(f"| {r['term']} | {r['decisions']} | {r['bound']} | {'YES' if r['fired'] else 'no (' + r['trigger'] + ')'} | {r['detail']} |")
    lines += ["", "## C. PREDICTED ERROR BAND per Type (predicted (thin + library - thick) / thick, pt)"]
    C = report.get("block_C")
    if C:
        lines.append("| Type | model | bias | H (x widened) | gap r|g| | band | LOO coverage (n) | extensions |")
        lines.append("|---|---|---|---|---|---|---|---|")
        form = (report.get("block_D") or {}).get("flags", {}).get("form")
        for t, row in C.items():
            if not row.get("calibrated"):
                lines.append(f"| {t} | — | | | | {row.get('note')} | | |")
                continue
            loo = row.get("loo_coverage") or {}
            widened = row.get("widened", 1.0)
            widen_txt = f"{row.get('H_pt', 0.0):.2f}" + (f" x {widened:.4f} WIDENED" if widened > 1.0 else " x 1 (not widened)")
            loo_txt = f"{loo.get('coverage', 'n/a')} CI {loo.get('ci90')} (n {row.get('n')})"
            if row.get("loo_coverage_unwidened") is not None:
                loo_txt += f"; unwidened {row['loo_coverage_unwidened']:.3f}: the >= 90 % is reached by the widening"
            ext = list(row.get("extensions", []))
            if t == "C" and form:
                ext.append("Maxwell FORM (P-embedded R2a port): the C column is FORM")
            lines.append(f"| {t} | {row.get('model')} | {fmt_pt(row.get('bias_pt'))} | {widen_txt} | "
                         f"{row.get('gap_pt', 0.0):.2f} | [{row['low_pt']:+.2f}, {row['high_pt']:+.2f}] | {loo_txt} | {'; '.join(ext)} |")
    else:
        lines.append("no band file")
    lines += ["", "## D. QUOTABILITY"]
    D = report["block_D"]
    for t, v in D["ladder"]["types"].items():
        lines.append(f"- {t}: {v['solves']} solves; steps (pt of the final ft) {', '.join(f'{s:+.3f}' for s in v['steps_pt'])}; "
                     f"last {fmt_pt(v['last_step_pt'], 3)}; ratio {v['step_ratio']}; {v['status']}; tail {fmt_pt(v['tail_pt'], 3)}")
    if D["ladder"].get("borrowed"):
        lines.append(f"- ladder borrowed from {D['ladder']['source']}")
    lines.append(f"- closure: {D.get('closure')}")
    lines.append(f"- ff vs ft: " + "; ".join(f"{t} {fmt_pct(v, 3)}" for t, v in D["ff_vs_ft"].items()))
    if D.get("flags"):
        lines.append(f"- Maxwell flags: {D['flags']}")
    return "\n".join(lines) + "\n"


def summary_line(report):
    parts = []
    A = report["block_A"]
    for t in report["types"]:
        if report.get("form") == "preflight":
            parts.append(f"{t} uncorrected {fmt_pct(report['ready_criterion_1'][t]['uncorrected_share_by_length'])} by length")
            continue
        e = A["energy"].get(t, {})
        band = (report.get("block_C") or {}).get(t, {})
        b = f" band [{band['low_pt']:+.2f}, {band['high_pt']:+.2f}]" if band.get("calibrated") else ""
        parts.append(f"{t} exact {fmt_pct(e.get('exact_share'))} / raw {fmt_pct(e.get('raw_share_total'))}{b}")
    fired = [r["term"] for r in report.get("block_B", []) if isinstance(r, dict) and r.get("fired")]
    return f"COVERAGE {report['label']['label']}: " + "; ".join(parts) + (f"; flags {', '.join(fired)}" if fired else "")


def write_report(report, out_dir, stem="coverage-report"):
    os.makedirs(out_dir, exist_ok=True)
    jpath = os.path.join(out_dir, stem + ".json")
    with open(jpath, "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")
    with open(os.path.join(out_dir, stem + ".md"), "w") as f:
        f.write(render_markdown(report))
    return jpath


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("postpro", nargs="?", help="the run's postpro directory (omit with --preflight)")
    parser.add_argument("--manifest", required=True, help="the preflight manifest surface-response-requirements.json")
    parser.add_argument("--bands", default=DEFAULT_BANDS, help="the calibration band file (default: coverage_bands_v2.json)")
    parser.add_argument("--no-bands", action="store_true", help="coverage only (blocks A / B / D; no band, no label)")
    parser.add_argument("--discovery", help="the discovery closure manifest (UnbuildableRequirements keys)")
    parser.add_argument("--ladder-run", help="a follower borrows the AMR ladder of this run directory")
    parser.add_argument("--reference-order", type=int,
                        help="3 = the MA_sharp reading is graded against a p3 reference (p3-REFERENCE: u_p added to the MA band; DESIGN 5.1)")
    parser.add_argument("--port-form", choices=PORT_FORMS,
                        help="r2a-pembed = a P-embedded R2a port (Maxwell FORM: the C column labelled FORM; decision 929 item D)")
    parser.add_argument("--ma-reading", choices=("raw", "sharp"), default="raw",
                        help="sharp = the MA_sharp transformation (thin_ma_sharp) is applied to this run's reading")
    parser.add_argument("--preflight", action="store_true", help="the preflight-only form (manifest + palace.json)")
    parser.add_argument("--palace-json", help="palace.json of the preflight (optional with --preflight)")
    parser.add_argument("--label", help="the run's name in the report")
    parser.add_argument("--output", help="output directory (default: beside the postpro)")
    args = parser.parse_args(argv)
    try:
        if args.preflight:
            report = build_preflight_report(args.manifest, args.palace_json, args.discovery, args.label)
            out_dir = args.output or os.path.dirname(os.path.abspath(args.manifest))
        else:
            if not args.postpro:
                parser.error("postpro is required without --preflight")
            bands = None if args.no_bands else load_bands(args.bands)
            report = build_report(args.postpro, args.manifest, bands, args.discovery, args.ladder_run, args.reference_order,
                                  args.ma_reading, args.label, args.port_form)
            out_dir = args.output or args.postpro
    except MissingInput as exc:
        print(f"coverage_report: MISSING INPUT {exc}", file=sys.stderr)
        return 2
    write_report(report, out_dir)
    print(summary_line(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
