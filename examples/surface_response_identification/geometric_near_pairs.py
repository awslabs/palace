#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The geometric near pairs of ANY surface-response library and their CORRECTION DIFFERENCE (identification fix round
DESIGN section 5 / A.4 / A.5 (b); decisions 1059 Q1 / Q4, 1079 Q1 (ii), 1097 Q4): the library tool that turns the
idfix probes v3j_fb2_response_calibration_v2.py / v3j_pec_nearpair_calibration_v1.py into a reusable reader.

A PAIR = two library models of one class (SpatialEdgeCluster / ParallelEdgeCluster) whose signatures lie inside the
geometric near-match band, so a Missing key within the band of both would be served by the nearer one. Its CORRECTION
DIFFERENCE per interface Type: D = sum over the Type's interfaces of (Q_fab - Q_thin)_ij, the within-R correction the
placement applies, compared by ||D_A - D_B||_F / ||D_A||_F when the two trace bases have the same dimension (Frobenius)
and, basis-size free, by the relative difference of trace(D) and of the top-10 |eigenvalues| (symmetrised). The bar
GEOMETRIC_NEAR_PAIR_CORRECTION_BAR (2 % per Type; A.4 (c), set AFTER the v3-j measurement = a calibration) flags a
pair: FLAGGED pairs are RECORDED and listed by the coverage report (an information row), never excluded.

THE PAIR LIST: with --manifest, the preflight manifest's Library.GeometricNearPairs (the C++'s in-band pairs of
record: the one-to-one context assignment of decision 1079 Q2) is THE authority. Without one the tool enumerates
candidates itself in the PROBE FORM (labelled "ProbeForm": a ParallelEdgeCluster pair by the same EdgeCount / GapSide
/ Conductor pattern with every |dOffset| <= GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS x R; a SpatialEdgeCluster pair by
the flattened same-field-set signature deviation: claims (Portions / Box / Vertices) <= the claims band and Context
<= the context band, SAME-ORDER fields) — NOT the C++'s one-to-one assignment: a report from it says so.

Usage: geometric_near_pairs.py --library process-library.json [--manifest surface-response-requirements.json]
       --out near-pairs.json [--radius R]
Record test: on the v3-j library of record the SpatialEdgeCluster pair a2f53cd8630e ~ 7d88636452b7 reads SA Frobenius
1.028e-3, 7d88 ~ 5bda 8.485e-3 and the ParallelEdgeCluster pair af972e5a4a5d ~ 186d141371db 9.935e-3 (the probes' .out).
"""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import identification_constants as ic  # noqa: E402

PAIR_CLASSES = ("SpatialEdgeCluster", "ParallelEdgeCluster")
TOP_EIGENVALUES = 10


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(path, library_path):
    """A matrix path of the library (absolute, or relative to the library file)."""
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(library_path)), path))


def load_surface_matrix(path, interface_types):
    """surface-response-matrix.csv -> ({Type: {(basis_i, basis_j): sum of Q_ij over the Type's interfaces}}, n):
    columns by NAME (interface, basis_i, basis_j, Q_ij (J)); the probes' form."""
    acc = {}
    n = 0
    with open(path) as f:
        reader = csv.reader(f)
        header = [c.strip() for c in next(reader)]
        col = {name: i for i, name in enumerate(header)}
        for key in ("interface", "basis_i", "basis_j", "Q_ij (J)"):
            if key not in col:
                raise ValueError(f"{path}: column {key!r} missing (columns {header[:6]}...)")
        for row in reader:
            if len(row) <= col["Q_ij (J)"]:
                continue
            interface = int(float(row[col["interface"]]))
            bi = int(float(row[col["basis_i"]]))
            bj = int(float(row[col["basis_j"]]))
            q = float(row[col["Q_ij (J)"]])
            t = interface_types.get(interface, f"I{interface}")
            acc.setdefault(t, {})
            acc[t][(bi, bj)] = acc[t].get((bi, bj), 0.0) + q
            n = max(n, bi, bj)
    return acc, n


def correction_matrices(model, library_path):
    """({Type: D = Q_fab - Q_thin}, basis size) of a model from its surface matrices."""
    types = {int(i["Coupon"]): i["Type"] for i in model.get("Interfaces", [])}
    fab, nf = load_surface_matrix(resolve(model["FabricatedSurfaceMatrix"], library_path), types)
    thin, nt = load_surface_matrix(resolve(model["ThinSurfaceMatrix"], library_path), types)
    if nf != nt:
        raise ValueError(f"{model['Name']}: fabricated basis {nf} != thin basis {nt}")
    D = {}
    for t in sorted(set(fab) | set(thin)):
        a = fab.get(t, {})
        b = thin.get(t, {})
        D[t] = {k: a.get(k, 0.0) - b.get(k, 0.0) for k in set(a) | set(b)}
    return D, nf


def frobenius(a, b=None):
    keys = set(a) | (set(b) if b else set())
    return math.sqrt(sum((a.get(k, 0.0) - (b.get(k, 0.0) if b else 0.0)) ** 2 for k in keys))


def eigenvalues(a, n):
    import numpy as np
    M = np.zeros((n, n))
    for (i, j), v in a.items():
        M[i - 1, j - 1] = v
    M = 0.5 * (M + M.T)
    return np.sort(np.linalg.eigvalsh(M))


def correction_difference(Da, na, Db, nb):
    """Per Type: {Frobenius (same basis only), Trace, Top10, Max}; Max = the largest of the applicable metrics."""
    import numpy as np
    out = {}
    for t in sorted(set(Da) & set(Db)):
        ea, eb = eigenvalues(Da[t], na), eigenvalues(Db[t], nb)
        ta, tb = float(np.sum(ea)), float(np.sum(eb))
        trace = abs(ta - tb) / max(abs(ta), 1e-300)
        ka = ea[np.argsort(-np.abs(ea))][:TOP_EIGENVALUES]
        kb = eb[np.argsort(-np.abs(eb))][:TOP_EIGENVALUES]
        m = min(len(ka), len(kb))
        top = float(np.max(np.abs(ka[:m] - kb[:m])) / max(float(np.max(np.abs(ka))), 1e-300))
        entry = {"Trace": trace, "Top10": top, "SameBasis": na == nb}
        if na == nb:
            entry["Frobenius"] = frobenius(Da[t], Db[t]) / max(frobenius(Da[t]), 1e-300)
        entry["Max"] = max(v for k, v in entry.items() if k in ("Frobenius", "Trace", "Top10"))
        out[t] = entry
    return out


# --------------------------------------------------------------------------------------------------
# The pair list: the manifest's Library.GeometricNearPairs (the authority) or the PROBE-FORM enumeration
# --------------------------------------------------------------------------------------------------


def flatten(obj, prefix=""):
    """(path, value) leaves of a JSON object in sorted key order (the probes' same-field-set reading)."""
    out = []
    if isinstance(obj, dict):
        for k in sorted(obj):
            out += flatten(obj[k], prefix + "/" + k)
    elif isinstance(obj, list):
        for i, x in enumerate(obj):
            out += flatten(x, prefix + f"[{i}]")
    elif isinstance(obj, bool) or not isinstance(obj, (int, float)):
        out.append((prefix, str(obj)))
    else:
        out.append((prefix, float(obj)))
    return out


def probe_form_pairs(library):
    """The in-band pairs by the probes' own reading (NOT the C++'s one-to-one context assignment)."""
    R = float(library["MatchingRadius"])
    pairs = []
    pec = [m for m in library["Models"] if m.get("Topology") == "ParallelEdgeCluster" and m.get("Edges")]
    for a, b in itertools.combinations(pec, 2):
        ea, eb = a["Edges"], b["Edges"]
        if len(ea) != len(eb) or [e.get("GapDirection") for e in ea] != [e.get("GapDirection") for e in eb] \
                or [e.get("Conductor", 1) for e in ea] != [e.get("Conductor", 1) for e in eb]:
            continue
        dev = max(abs(float(x["Offset"]) - float(y["Offset"])) / R for x, y in zip(ea, eb))
        if dev <= ic.GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS:
            pairs.append({"Class": "ParallelEdgeCluster", "ModelA": a["Name"], "ModelB": b["Name"],
                          "MaxOffsetDeviationOverR": dev, "Form": "ProbeForm"})
    sec = [m for m in library["Models"] if m.get("Topology") == "SpatialEdgeCluster" and isinstance(m.get("Signature"), dict)]
    for a, b in itertools.combinations(sec, 2):
        fa, fb = dict(flatten(a["Signature"])), dict(flatten(b["Signature"]))
        if set(fa) != set(fb):
            continue
        claims = context = 0.0
        same_text = True
        for p in fa:
            if isinstance(fa[p], float) and isinstance(fb[p], float):
                d = abs(fa[p] - fb[p])
                if "/Context" in p:
                    context = max(context, d)
                else:
                    claims = max(claims, d)
            elif fa[p] != fb[p]:
                same_text = False
        if same_text and claims <= ic.GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS and context <= ic.GEOMETRIC_NEAR_MATCH_CONTEXT_OVER_RADIUS:
            pairs.append({"Class": "SpatialEdgeCluster", "ModelA": a["Name"], "ModelB": b["Name"],
                          "ClaimsDeviationOverR": claims, "ContextDeviationOverR": context, "Form": "ProbeForm"})
    return pairs


def manifest_pairs(manifest):
    pairs = ((manifest.get("Library") or {}).get("GeometricNearPairs")) or []
    out = []
    for p in pairs:
        entry = dict(p)
        entry["Form"] = "Manifest"
        out.append(entry)
    return out


def near_pairs_record(library_path, manifest_path=None, bar=ic.GEOMETRIC_NEAR_PAIR_CORRECTION_BAR):
    with open(library_path) as f:
        library = json.load(f)
    manifest = None
    if manifest_path:
        with open(manifest_path) as f:
            manifest = json.load(f)
    pairs = manifest_pairs(manifest) if manifest is not None else probe_form_pairs(library)
    models = {m["Name"]: m for m in library["Models"]}
    cache = {}

    def matrices(name):
        if name not in cache:
            cache[name] = correction_matrices(models[name], library_path)
        return cache[name]

    rows = []
    for p in pairs:
        row = dict(p)
        missing = [n for n in (p["ModelA"], p["ModelB"]) if n not in models]
        if missing:
            row["CorrectionDifference"] = None
            row["Note"] = f"model(s) not in the library: {', '.join(missing)}"
            rows.append(row)
            continue
        Da, na = matrices(p["ModelA"])
        Db, nb = matrices(p["ModelB"])
        diff = correction_difference(Da, na, Db, nb)
        row["CorrectionDifference"] = diff
        row["BasisSizes"] = [na, nb]
        row["MaxOverTypes"] = max((v["Max"] for v in diff.values()), default=None)
        row["Flagged"] = bool(diff) and row["MaxOverTypes"] > bar
        rows.append(row)
    return {
        "Tool": "geometric_near_pairs.py",
        "Library": {"Path": os.path.abspath(library_path), "Name": library.get("Name"), "SHA256": sha256_file(library_path),
                    "MatchingRadius": library.get("MatchingRadius"), "Models": len(library["Models"])},
        "Manifest": ({"Path": os.path.abspath(manifest_path), "SHA256": sha256_file(manifest_path)} if manifest_path else None),
        "PairSource": "Manifest Library.GeometricNearPairs (the C++'s in-band pairs of record, decision 1079 Q1 (ii))" if manifest is not None
                      else "ProbeForm: the tool's own enumeration (PEC offsets within the claims band; SEC flattened same-order signature "
                           "deviation within the claims / context bands) - NOT the C++'s one-to-one context assignment",
        "Bands": {"ClaimsOverR": ic.GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS, "ContextOverR": ic.GEOMETRIC_NEAR_MATCH_CONTEXT_OVER_RADIUS},
        "CorrectionBar": bar,
        "Rule": "DESIGN A.4 (c) / A.5 (b): per Type D = sum over the Type's interfaces of (Q_fab - Q_thin); Frobenius where the bases agree, "
                "else the trace and the top-10 |eigenvalue| spectrum; a pair over the bar is RECORDED and listed (information), never excluded",
        "Pairs": rows,
        "Flagged": [[r["ModelA"], r["ModelB"], r["MaxOverTypes"]] for r in rows if r.get("Flagged")],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", required=True, help="process-library.json")
    parser.add_argument("--manifest", help="a preflight manifest: its Library.GeometricNearPairs is the pair list (else the ProbeForm enumeration)")
    parser.add_argument("--out", required=True, help="the record to write (JSON)")
    parser.add_argument("--bar", type=float, default=ic.GEOMETRIC_NEAR_PAIR_CORRECTION_BAR, help="the per-Type correction-difference bar (default 0.02)")
    args = parser.parse_args(argv)
    record = near_pairs_record(args.library, args.manifest, args.bar)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(record, f, indent=1, sort_keys=True)
        f.write("\n")
    for r in record["Pairs"]:
        diff = r.get("CorrectionDifference") or {}
        print(f"{r['Class']} {r['ModelA'][-12:]} ~ {r['ModelB'][-12:]} ({r['Form']}): "
              + " ".join(f"{t} {v['Max']:.3e}" for t, v in sorted(diff.items())) + (" FLAGGED" if r.get("Flagged") else ""))
    print(f"{len(record['Pairs'])} pair(s), {len(record['Flagged'])} over the {record['CorrectionBar']:.0%} bar -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
