#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Held-out interpolation gate of an angle-interpolated corner family (USER decision 121 (C);
gate form USER decision 149 (5), 2026-09-29): the runtime's kink-aware stencil
(MatchCornerFamily / SelectCornerFamilyStencil, mirrored by corner_family_interpolation.
select_stencil: segments of coupons sharing TraceBasis ConnectivityAngleDegrees, Lagrange on the
segment's nodes nearest to the angle, never across a knot-corner passage) applied to the node
coupons' response matrices and compared with the coupon actually built at every held-out angle.

The gate is the CornerFamilyInterpolation entry of spatial_coupon/qualify/qualification-gates.json:
PARTICIPATION-REFERENCED — for every held-out coupon and interface SA / MS / MA the residual of
the interpolated fabricated surface matrix on the held-out trace, referred to the held-out
coupon's own fabricated energy of that interface, within MaximumParticipationReferencedResidual.
Reported alongside, never gated: the DEFECT-REFERENCED residuals (fabricated - thin, relative to
the defect energy), the corrected participation's error (the defect residual over the
fabricated energy) and the Frobenius residuals. A held-out angle the stencil rule refuses fails
the check.

The held-out TRACE the residuals are evaluated on is recorded (review of decision 154, MAJOR-1:
the verdict depends on it). The GATING trace is the gate file's GatingTrace = option-c (USER
decision 161 (2), 2026-09-30): a run whose trace forms are exactly that one is the family's
verdict (JSON Gating true), any other run (the band trace of a recorded cache) is reported
alongside. --trace recorded (the default) uses the cache's heldout-coefficients.csv
and classifies it from basis-points.csv + coupon-spec.json against the generator's two forms —
TraceForm "band" (the pre-149 (6) metal-band cutoff: exactly zero on both metal rings, blind
to the free knots next to the metal) or "option-c" (decision 149 (6): zero on the PEC part
only) — or "recorded-unclassified" when it matches neither or the cache has no basis points.
--trace option-c recomputes the option-(c) coefficients from basis-points.csv + the coupon
spec with generate_corner_response.heldout_potential (the committed generator), so the verdict
is a property of the family and not of the generator version that wrote the cache. Which trace
GATES a family is a USER decision; the tool records the form it used (TraceForm in the JSON
record, trace_form in the CSV) and never chooses. Fail closed: a held-out coefficient vector
whose length differs from the matrix size, or a zero fabricated held-out energy, is an error.

    corner_family_heldout_check.py OUT.csv COUPON_DIR... [--nodes 75:82.5 90:82.5 ...]
                                   [--trace recorded|option-c]
                                   [--gates qualification-gates.json] [--json OUT.json]

COUPON_DIR = a corner coupon cache directory (coupon-spec.json, process-library.json,
heldout-coefficients.csv, basis-points.csv, postpro/{fabricated,thin}/{domain-response-matrix,
surface-response-matrix-aggregate}.csv). Without --nodes the nodes are the coupons whose model
carries TraceBasis.ConnectivityAngleDegrees; with it, ANGLE:CONNECTIVITY specs select the nodes
(a recorded coupon without the record is stamped with the connectivity, the caller's
responsibility: no basis event between the two angles, which the C++ load check verifies). An
AllRingsFollowMetal family (TraceBasis.RingLayout, corner-basis refinement 2026-09-30) has no
connectivity records and no events: its nodes are given as plain ANGLE specs (--nodes 75 90 ...),
the stencil is the cubic sliding window of the one segment."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import corner_family_interpolation as family  # noqa: E402
import generate_corner_response as generator  # noqa: E402

GATES_FILE = HERE.parent / "spatial_coupon" / "qualify" / "qualification-gates.json"
GATE = "CornerFamilyInterpolation"
SIGNATURE_ANGLE_TOLERANCE_DEGREES = family.SIGNATURE_ANGLE_TOLERANCE_DEGREES
INTERFACES = {1: "SA", 2: "MS", 3: "MA"}
DEFECTS = ("domain", "SA", "MS", "MA")
FABRICATED = ("domain_fab", "SA_fab", "MS_fab", "MA_fab")
GATED = ("SA_fab", "MS_fab", "MA_fab")
TRACE_SOURCES = ("recorded", "option-c")
TRACE_FORM_BAND, TRACE_FORM_OPTION_C, TRACE_FORM_UNCLASSIFIED = "band", "option-c", "recorded-unclassified"
TRACE_CLASSIFICATION_TOLERANCE = 1.0e-9


def read_domain(path):
    rows = list(csv.reader(open(path)))[1:]
    n = max(int(float(r[1])) for r in rows)
    m = np.zeros((n, n))
    for r in rows:
        i, j, q = int(float(r[0])) - 1, int(float(r[1])) - 1, float(r[2])
        m[i, j] = m[j, i] = q
    return m


def read_surface(path):
    rows = list(csv.reader(open(path)))
    header = [h.strip() for h in rows[0]]
    ci, cj, cq = header.index("basis_i"), header.index("basis_j"), header.index("Q_total_ij (J)")
    cint, cedge = header.index("interface"), header.index("edge")
    out, seen = {}, set()
    for r in rows[1:]:
        if not r:
            continue
        key = int(float(r[cint]))
        i, j = int(float(r[ci])) - 1, int(float(r[cj])) - 1
        dedup = (key, int(float(r[cedge])), i, j)
        if dedup in seen:
            continue
        seen.add(dedup)
        q = float(r[cq])
        m = out.setdefault(key, {})
        m[(i, j)] = m.get((i, j), 0.0) + q
        if i != j:
            m[(j, i)] = m.get((j, i), 0.0) + q
    n = max(max(i, j) for m in out.values() for (i, j) in m) + 1
    return {k: np.array([[m.get((i, j), 0.0) for j in range(n)] for i in range(n)]) for k, m in out.items()}


def generator_traces(d, spec):
    """The generator's two held-out traces at the coupon's recorded basis points, or None when
    the cache carries no basis-points.csv (a synthetic coupon)."""
    basis_points = d / "basis-points.csv"
    if not basis_points.is_file():
        return None
    points = np.loadtxt(basis_points, delimiter=",", skiprows=1, ndmin=2)
    radius = float(spec["MatchingRadius"])
    thickness = float(spec["Fabrication"]["metal_thickness"])
    angle, topology = float(spec["AngleDegrees"]), spec["Topology"]
    polynomial = generator.heldout_polynomial(points, radius)
    return {
        TRACE_FORM_BAND: generator.metal_band_cutoff(points, radius, thickness) * polynomial,
        TRACE_FORM_OPTION_C: generator.heldout_potential(points, radius, thickness, angle, topology),
    }


def heldout_trace(d, spec, trace_source):
    """(trace, TraceForm) for one coupon directory under the requested trace source."""
    forms = generator_traces(d, spec)
    if trace_source == "option-c":
        if forms is None:
            raise FileNotFoundError(f"{d}: --trace option-c needs basis-points.csv to recompute the held-out trace")
        return forms[TRACE_FORM_OPTION_C], TRACE_FORM_OPTION_C
    if trace_source != "recorded":
        raise ValueError(f"unknown trace source {trace_source!r}; one of {TRACE_SOURCES}")
    trace = np.loadtxt(d / "heldout-coefficients.csv", delimiter=",", skiprows=1, ndmin=1)
    form = TRACE_FORM_UNCLASSIFIED
    if forms is not None:
        for name, reference in forms.items():
            if len(reference) == len(trace) and np.allclose(trace, reference, rtol=TRACE_CLASSIFICATION_TOLERANCE,
                                                            atol=TRACE_CLASSIFICATION_TOLERANCE * np.max(np.abs(reference))):
                form = name
    return trace, form


def coupon(directory, trace_source="recorded"):
    d = Path(directory)
    spec = json.load(open(d / "coupon-spec.json"))
    fab_dom = read_domain(d / "postpro/fabricated/domain-response-matrix.csv")
    thin_dom = read_domain(d / "postpro/thin/domain-response-matrix.csv")
    fab = read_surface(d / "postpro/fabricated/surface-response-matrix-aggregate.csv")
    thin = read_surface(d / "postpro/thin/surface-response-matrix-aggregate.csv")
    quantities = {"domain": fab_dom - thin_dom, "domain_fab": fab_dom}
    for k, name in INTERFACES.items():
        quantities[name] = fab[k] - thin[k]
        quantities[name + "_fab"] = fab[k]
    trace, trace_form = heldout_trace(d, spec, trace_source)
    size = fab_dom.shape[0]
    if any(q.shape != (size, size) for q in quantities.values()):
        raise ValueError(f"{d}: the domain and surface response matrices differ in size")
    if len(trace) != size:
        raise ValueError(f"{d}: {len(trace)} held-out coefficients for {size} basis functions")
    connectivity = spec.get("ConnectivityAngleDegrees")
    rule = None
    library = d / "process-library.json"
    if library.exists():
        model = json.load(open(library))["Models"][0]
        record = model.get("TraceBasis") or {}
        connectivity = record.get("ConnectivityAngleDegrees", connectivity)
        if "MetalInteriorKnots" in record and "FreeKnots" in record:
            rule = generator.TraceBasisRule.from_record(record)
    return {
        "topology": spec["Topology"],
        "angle": float(spec["AngleDegrees"]),
        "connectivity": None if connectivity is None else float(connectivity),
        "rule": rule,
        "quantities": quantities,
        "trace": trace,
        "trace_form": trace_form,
        "dir": str(d),
    }


def family_rule(fam):
    """The one trace basis rule of a family's coupons (None for a family without TraceBasis
    records); mixed rules are an error (the runtime refuses them too)."""
    rules = [c["rule"] for c in fam if c["rule"] is not None]
    if not rules:
        return None
    if any(r != rules[0] for r in rules) or len(rules) != len(fam):
        raise ValueError("the family's coupons are not built on one TraceBasis rule")
    return rules[0]


def convexity_of(topology):
    return "convex" if topology.lower().startswith("convex") else "concave"


def select_nodes(fam, node_specs):
    """[(angle, connectivity, index, coupon)] and the held-out coupons of one topology."""
    nodes, used = [], set()
    if node_specs is None:
        for c in sorted(fam, key=lambda c: (c["angle"], c["connectivity"] or 0.0)):
            if c["connectivity"] is not None:
                nodes.append((c["angle"], c["connectivity"], len(nodes), c))
                used.add(c["dir"])
    else:
        for angle, connectivity in node_specs:
            candidates = [c for c in fam if abs(c["angle"] - angle) <= SIGNATURE_ANGLE_TOLERANCE_DEGREES
                          and c["dir"] not in used
                          and (c["connectivity"] is None or connectivity is None
                               or abs(c["connectivity"] - connectivity) <= 1e-9)]
            exact = [c for c in candidates if c["connectivity"] is not None]
            pick = exact[0] if exact else (candidates[0] if candidates else None)
            if pick is None:
                print(f"WARNING: no {convexity_of(fam[0]['topology'])} coupon for node {angle:g}:{connectivity}")
                continue
            nodes.append((angle, connectivity, len(nodes), pick))
            used.add(pick["dir"])
    return nodes, [c for c in fam if c["dir"] not in used]


def evaluate(coupons, node_specs, gate):
    """Rows of the residual table and the per-family verdict record."""
    limit = 100.0 * gate["MaximumParticipationReferencedResidual"]
    table = [("topology", "heldout_angle", "rule", "connectivity_angle", "abscissae_angles", "weights", "quantity",
              "frobenius_residual_%", "heldout_energy_residual_%", "heldout_energy_residual_over_fabricated_%", "gated",
              "trace_form")]
    families = {}
    for topology in sorted({c["topology"] for c in coupons}):
        fam = [c for c in coupons if c["topology"] == topology]
        rule = family_rule(fam)
        if rule is not None and rule.all_rings and node_specs is None:
            raise ValueError(f"{topology}: an AllRingsFollowMetal family has no connectivity records to tell "
                             "the nodes from the held-out coupons: pass --nodes ANGLE ...")
        nodes, held_out = select_nodes(fam, node_specs)
        stencil_nodes = [(a, k, i) for a, k, i, _ in nodes]
        record = {"Topology": topology, "Nodes": [{"AngleDegrees": a, "ConnectivityAngleDegrees": k, "Dir": c["dir"]}
                                                   for a, k, _, c in nodes],
                  "HeldOut": [], "Refused": [], "TraceForms": sorted({h["trace_form"] for h in held_out}),
                  "WorstParticipationReferencedPercent": {},
                  "WorstDefectReferencedPercent": {}, "WorstCorrectedParticipationPercent": {}}
        for h in sorted(held_out, key=lambda c: c["angle"]):
            if not stencil_nodes:
                stencil = {"reason": "the family has no nodes (no coupon with a TraceBasis.ConnectivityAngleDegrees "
                                     "record and no --nodes list): a legacy family is never interpolated"}
            else:
                stencil = family.select_stencil(stencil_nodes, h["angle"], convexity_of(topology),
                                                ring_size=(rule or generator.LEGACY_RULE).ring_size, rule=rule)
            if "reason" in stencil:
                record["Refused"].append({"AngleDegrees": h["angle"], "Reason": stencil["reason"], "Dir": h["dir"]})
                table.append((topology, f"{h['angle']:g}", "refused: " + stencil["reason"], "", "", "", "", "", "", "", "",
                              h["trace_form"]))
                continue
            window = [nodes[i][3] for i, _ in stencil["nodes"]]
            weights = [w for _, w in stencil["nodes"]]
            entry = {"AngleDegrees": h["angle"], "Rule": stencil["rule"], "Dir": h["dir"], "TraceForm": h["trace_form"],
                     "ConnectivityAngleDegrees": stencil["connectivity_angle_degrees"],
                     "Stencil": [c["angle"] for c in window], "Weights": weights,
                     "ParticipationReferencedPercent": {}, "DefectReferencedPercent": {},
                     "CorrectedParticipationPercent": {}, "FrobeniusPercent": {}}
            for quantity in DEFECTS + FABRICATED:
                actual = h["quantities"][quantity]
                interp = sum(w * c["quantities"][quantity] for w, c in zip(weights, window))
                frob = 100.0 * np.linalg.norm(interp - actual) / np.linalg.norm(actual)
                v = h["trace"]
                e_actual = v @ actual @ v
                fab_quantity = quantity if quantity.endswith("_fab") else quantity + "_fab"
                e_fabricated = v @ h["quantities"][fab_quantity] @ v
                if e_actual == 0.0 or e_fabricated == 0.0:
                    raise ValueError(f"{h['dir']}: zero held-out energy of {quantity} (trace form {h['trace_form']}): "
                                     "the residual is undefined")
                energy = 100.0 * (v @ interp @ v - e_actual) / e_actual
                over_fab = 100.0 * (v @ interp @ v - e_actual) / e_fabricated
                entry["FrobeniusPercent"][quantity] = frob
                if quantity.endswith("_fab"):
                    entry["ParticipationReferencedPercent"][quantity[:-4]] = over_fab
                else:
                    entry["DefectReferencedPercent"][quantity] = energy
                    entry["CorrectedParticipationPercent"][quantity] = over_fab
                gated = quantity in GATED
                connectivity_text = ("" if stencil["connectivity_angle_degrees"] is None
                                     else f"{stencil['connectivity_angle_degrees']:g}")
                table.append((topology, f"{h['angle']:g}", stencil["rule"], connectivity_text,
                              "/".join(f"{c['angle']:g}" for c in window), "/".join(f"{w:+.4f}" for w in weights),
                              quantity, f"{frob:.4f}", f"{energy:+.4f}", f"{over_fab:+.4f}", "yes" if gated else "no",
                              h["trace_form"]))
            entry["Passed"] = all(abs(entry["ParticipationReferencedPercent"][x]) <= limit for x in ("SA", "MS", "MA"))
            record["HeldOut"].append(entry)
            for key, source in (("WorstParticipationReferencedPercent", "ParticipationReferencedPercent"),
                                ("WorstDefectReferencedPercent", "DefectReferencedPercent"),
                                ("WorstCorrectedParticipationPercent", "CorrectedParticipationPercent")):
                for quantity, value in entry[source].items():
                    if abs(value) >= abs(record[key].get(quantity, (0.0, None))[0]):
                        record[key][quantity] = (value, h["angle"])
        record["Passed"] = bool(record["HeldOut"]) and not record["Refused"] and all(e["Passed"] for e in record["HeldOut"])
        families[topology] = record
    return table, families


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out")
    parser.add_argument("coupons", nargs="+")
    parser.add_argument("--nodes", nargs="*", default=None,
                        help="node coupons as ANGLE:CONNECTIVITY or ANGLE (legacy); default: the coupons with a "
                             "TraceBasis.ConnectivityAngleDegrees record")
    parser.add_argument("--trace", choices=TRACE_SOURCES, default="recorded",
                        help="the held-out trace: the cache's heldout-coefficients.csv (classified band | option-c | "
                             "recorded-unclassified) or the option-(c) coefficients recomputed from basis-points.csv + "
                             "the coupon spec; the family's VERDICT is the run on the gate file's GatingTrace "
                             "(option-c, USER decision 161 (2)), any other run is reported alongside")
    parser.add_argument("--gates", type=Path, default=GATES_FILE)
    parser.add_argument("--json", type=Path, help="the verdict record (the gate file's digest, both forms, the verdict)")
    args = parser.parse_args()
    gates = json.load(open(args.gates))
    gate = gates["Gates"][GATE]
    node_specs = None
    if args.nodes is not None:
        node_specs = []
        for spec in args.nodes:
            angle, _, connectivity = spec.partition(":")
            node_specs.append((float(angle), float(connectivity) if connectivity else None))
    table, families = evaluate([coupon(d, args.trace) for d in args.coupons], node_specs, gate)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    csv.writer(open(args.out, "w")).writerows(table)
    for row in table:
        print(",".join(row))
    limit = 100.0 * gate["MaximumParticipationReferencedResidual"]
    passed = bool(families) and all(f["Passed"] for f in families.values())
    for topology, f in families.items():
        fmt = lambda d: {q: f"{v:+.3f} @ {a:g}" for q, (v, a) in sorted(d.items())}  # noqa: E731
        print(f"{topology}: participation-referenced (GATED, |residual| <= {limit:g} % of the fabricated energy): "
              f"{fmt(f['WorstParticipationReferencedPercent'])}")
        print(f"{topology}: corrected participation error (defect residual / fabricated energy, reported): "
              f"{fmt(f['WorstCorrectedParticipationPercent'])}")
        print(f"{topology}: defect-referenced (reported): {fmt(f['WorstDefectReferencedPercent'])}")
        if f["Refused"]:
            print(f"{topology}: REFUSED held-out angles: {[(r['AngleDegrees'], r['Reason']) for r in f['Refused']]}")
        print(f"{topology}: {len(f['HeldOut'])} held-out coupons, {len(f['Nodes'])} nodes, held-out trace form "
              f"{'/'.join(f['TraceForms'])} -> {'PASS' if f['Passed'] else 'FAIL'}")
    gating_trace = gate.get("GatingTrace", "option-c")
    forms_used = sorted({form for f in families.values() for form in f["TraceForms"]})
    gating = forms_used == [gating_trace]
    print(f"{GATE} (participation-referenced, {limit:g} %; USER decision 149 (5)) on the {args.trace} held-out trace "
          f"(forms {'/'.join(forms_used)}; {'the GATING trace' if gating else 'NOT the gating trace ' + gating_trace + ': reported alongside'}): "
          f"{'PASS' if passed else 'FAIL'} -> {args.out}")
    if args.json:
        record = {"Version": 3, "Gate": GATE, "GatesFile": str(args.gates.resolve()),
                  "TraceSource": args.trace, "GatingTrace": gating_trace, "Gating": gating,
                  "TraceForms": sorted({form for f in families.values() for form in f["TraceForms"]}),
                  "GatesSHA256": hashlib.sha256(args.gates.read_bytes()).hexdigest(),
                  "GatesVersion": gates.get("Version"),
                  "MaximumParticipationReferencedResidual": gate["MaximumParticipationReferencedResidual"],
                  "Families": families, "Passed": passed}
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(record, indent=2, default=str) + "\n")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
