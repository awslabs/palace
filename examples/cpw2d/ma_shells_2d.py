#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Per-shell MA and the sharp-edge extrapolation of a 2D StraightEdgeBuilder response run
(decision 66 part C: the 2D models of a process library report their MA with the same
cutoff bookkeeping as the 3D spatial coupons).

A --edge-distances run of generate_edge_response.py / generate_edge_pair_response.py
localizes every interface energy at the ring radii R_1 < ... < R_K (< 0.2 um, the historical
coupon radius, kept as the last entry): surface-response-matrix.csv holds, per (interface,
edge, R, i, j), Q_ij = the energy of the quadrature points whose NEAREST metal edge is
`edge` and whose distance to it is below R (cumulative), and Q_total_ij = the whole energy
attributed to that edge.  The ring k energy is Q(R_k) - Q(R_(k-1)) with R_0 = 0: ring 1
[0, R_1) holds the singularity the elements cannot represent (R_1 = the mesh's edge size =
the recorded MA cutoff), rings 2..K are resolved - the 3D radial shells' definition
(relabel_radial_ma_shells.py: shell = ring interval of the centroid distance to the nearest
metal edge line, kind = that line).

Fabricated coupon: the MA edge points are the sidewall endpoints (ma_edge_attributes), Palace
orders them by (x, y) so the bottom corner (the metal / trench edge at the process plane)
precedes the top corner (the sharp 90-degree metal edge) of each sidewall: odd CSV edges are
`bottom`, even ones `top`.  Per source the same-kind edges' rings are summed (as the 3D shells
sum every top line into one ring label) and the spatial coupons' estimator is applied verbatim
(radial_ma_profile.profile_of + ma_tail's Consistent rule: top Theory@2, bottom Fit2-4):
MA_raw = the sum of the edges' Q_total, MA_tail = the remainders, MA_sharp = MA_raw + MA_tail.

Thin coupon: the sheet edge carries MS below and MA above one surface; its MA is
log-divergent in the cutoff by construction (decision 66: the thin coupon is the device-
consistent thin-metal term, recorded at its cutoff, never extrapolated): the per-shell
profile and MA_raw are reported with the cutoff R_1, no tail.

usage: ma_shells_2d.py --postpro DIR --kind fabricated|thin [--ma-interface 3] --out-json PATH [--out-md PATH]
"""
import argparse
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
QUALIFY = HERE.parent / "cpw3d_surface" / "spatial_coupon" / "qualify"
sys.path.insert(0, str(QUALIFY))
import compare_matrices  # noqa: E402
import ma_tail  # noqa: E402
import radial_ma_profile as profile  # noqa: E402

COUPON_RADIUS = 0.2
KINDS = ("fabricated", "thin")
EDGE_KIND_RULE = ("fabricated MA edge points = the sidewall endpoints in Palace's (x, y) order: CSV edge 2m-1 = the bottom "
                  "(process-plane) corner, edge 2m = the top corner of sidewall m")
THIN_RULE = ("thin coupon: MA_raw at the recorded cutoff (ring 1 = the sheet edge's element size), no extrapolation - the "
             "thin sheet edge's MA is log-divergent in the cutoff by construction (decision 66)")


def edge_kind(edge, kind):
    if kind == "thin":
        return "sheet"
    return "bottom" if edge % 2 == 1 else "top"


def ring_energies(surface, interface, i):
    """{edge: [(inner, outer, Q_ring)] over the radii} of the diagonal source i; the last
    radius (the coupon radius) is the far ring."""
    by_edge = {}
    for (index, edge, radius, a, b), values in surface.items():
        if index != interface or a != i or b != i:
            continue
        by_edge.setdefault(edge, {})[radius] = values
    out = {}
    for edge, per_radius in sorted(by_edge.items()):
        radii = sorted(per_radius)
        rings, previous, cumulative = [], 0.0, 0.0
        for radius in radii:
            within = per_radius[radius]["Q_ij (J)"]
            rings.append((previous * 1.0e6, radius * 1.0e6, within - cumulative))
            previous, cumulative = radius, within
        out[edge] = {"Rings": rings, "Q_total": per_radius[radii[-1]]["Q_total_ij (J)"]}
    return out


def analyze(postpro, *, kind, ma_interface=3):
    """Per source (positive domain energy): the per-kind ring profile, MA_raw, MA_tail
    (fabricated), MA_sharp, the deficits of the spread estimators; plus the summary."""
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    postpro = Path(postpro)
    domain, surface = compare_matrices.load(postpro)
    radii = sorted({key[2] for key in surface if key[0] == ma_interface})
    if len(radii) < 2 or abs(radii[-1] * 1.0e6 - COUPON_RADIUS) > 1.0e-9:
        raise ValueError(f"{postpro}: the MA interface {ma_interface} carries no ring radii ending at the coupon radius "
                         f"{COUPON_RADIUS} um (found {[r * 1.0e6 for r in radii]} um): run with --edge-distances")
    ring_radii = [r * 1.0e6 for r in radii[:-1]]
    cutoff = ring_radii[0]
    per_source = {}
    for (i, j), energy in sorted(domain.items()):
        if i != j or energy <= 0.0:
            continue
        edges = ring_energies(surface, ma_interface, i)
        if not edges:
            raise ValueError(f"{postpro}: no MA rows for source {i}")
        if kind == "fabricated" and len(edges) % 2:
            raise ValueError(f"{postpro}: {len(edges)} MA edges - the fabricated MA edge points come in (bottom, top) pairs")
        raw = math.fsum(entry["Q_total"] for entry in edges.values())
        by_kind = {}
        for edge, entry in edges.items():
            rings = by_kind.setdefault(edge_kind(edge, kind), {})
            for ring, (inner, outer, q) in enumerate(entry["Rings"], start=1):
                inner0, outer0, q0 = rings.get(ring, (inner, outer, 0.0))
                rings[ring] = (inner, outer, q0 + q)
        record = {"E": energy, "Q_MA_raw": raw, "p_MA_raw": raw / energy, "Cutoff": cutoff, "Kinds": {}}
        remainders = {name: 0.0 for name in ma_tail.SPREAD_ESTIMATORS}
        for name, rings in by_kind.items():
            # The far ring (the coupon radius) is outside the tube-ring profile.
            fitted_rings = {ring: value for ring, value in rings.items() if ring <= len(ring_radii)}
            share = sum(v[2] for v in rings.values()) / raw if raw else None
            kind_record = {"Share": share, "Rings": {str(ring): {"InnerRadius": v[0], "OuterRadius": v[1], "Q": v[2]}
                                                      for ring, v in sorted(rings.items())}}
            if kind == "fabricated":
                fitted = profile.profile_of(fitted_rings)
                consistent = fitted["Estimates"].get(profile.kind_estimator(ma_tail.ESTIMATOR, name))
                kind_record.update({"Estimator": profile.kind_estimator(ma_tail.ESTIMATOR, name),
                                    "Remainder": consistent["Remainder"] if consistent else 0.0,
                                    "Estimates": {est: fitted["Estimates"].get(est) for est in ("Fit2-4", "Theory@2")},
                                    "LocalSlopes": fitted["LocalSlopes"], "CleanPowerLaw": fitted["CleanPowerLaw"]})
                for est in ma_tail.SPREAD_ESTIMATORS:
                    estimate = fitted["Estimates"].get(profile.kind_estimator(est, name))
                    if estimate:
                        remainders[est] += estimate["Remainder"]
            else:
                ordered = [fitted_rings[k] for k in sorted(fitted_rings)]
                kind_record["LocalSlopes"] = profile.local_slopes(ordered)
            record["Kinds"][name] = kind_record
        if kind == "fabricated":
            top = record["Kinds"].get("top", {}).get("Estimates", {})
            alpha, ring1 = top.get(ma_tail.ALPHA_ESTIMATOR), top.get(ma_tail.RING1_ESTIMATOR)
            tail = remainders[ma_tail.ESTIMATOR]
            record.update({"Q_MA_tail": tail, "Q_MA_sharp": raw + tail, "p_MA_sharp": (raw + tail) / energy,
                           "Alpha": alpha["Alpha"] if alpha else None, "AlphaSE": alpha["AlphaSE"] if alpha else None,
                           "Ring1Factor": ring1["Ring1ResolvedOverModel"] if ring1 else None,
                           "Deficit": {name: (value / raw if raw else None) for name, value in remainders.items()}})
        per_source[i] = record
    out = {"Version": 1, "Kind": kind, "Postpro": str(postpro), "MAInterface": ma_interface, "RingRadii": ring_radii,
           "Cutoff": cutoff, "CouponRadius": COUPON_RADIUS, "LengthUnit": "um",
           "Rule": (ma_tail.RULE + "; " + EDGE_KIND_RULE) if kind == "fabricated" else THIN_RULE,
           "PerSource": {str(i): record for i, record in per_source.items()}}
    if kind == "fabricated":
        out["Summary"] = ma_tail.summary(per_source, sorted(per_source), at=())
        out["Summary"]["p_MA"] = {"RawMedian": ma_tail._median([r["p_MA_raw"] for r in per_source.values()]),
                                  "SharpMedian": ma_tail._median([r["p_MA_sharp"] for r in per_source.values()])}
    else:
        slopes = [s for r in per_source.values() for s in r["Kinds"]["sheet"]["LocalSlopes"] if s is not None]
        out["Summary"] = {"Sources": len(per_source), "Estimator": None,
                          "p_MA": {"RawMedian": ma_tail._median([r["p_MA_raw"] for r in per_source.values()])},
                          "LocalSlopeMedian": statistics.median(slopes) if slopes else None}
    return out


def markdown(record):
    lines = [f"# 2D MA shells - {record['Kind']} ({record['Postpro']})", "",
             f"Ring radii (um): {record['RingRadii']}; cutoff (ring 1 outer radius) {record['Cutoff'] * 1000:g} nm; "
             f"coupon radius {record['CouponRadius']} um; sources {record['Summary']['Sources']}.", ""]
    summary = record["Summary"]
    if record["Kind"] == "fabricated":
        deficit = summary["Deficit"]
        lines += [f"- p_MA median raw {summary['p_MA']['RawMedian']:.4e} -> sharp {summary['p_MA']['SharpMedian']:.4e}",
                  f"- deficit ({summary['Estimator']}) median {ma_tail.pc(deficit[summary['Estimator']]['Median'])}, "
                  f"quartiles {deficit[summary['Estimator']]['Quartiles']}; spread Fit2-4 {ma_tail.pc(deficit['Fit2-4']['Median'])} "
                  f"/ Theory@2 {ma_tail.pc(deficit['Theory@2']['Median'])}",
                  f"- alpha top (Fit2-4) median {ma_tail.f3(summary['Alpha']['Median'])} (theory {summary['Alpha']['Theoretical']:.3f}), "
                  f"ring-1 factor median {ma_tail.f3(summary['Ring1Factor']['Median'])}",
                  "", "Rule: " + record["Rule"]]
    else:
        lines += [f"- p_MA median raw {summary['p_MA']['RawMedian']:.4e} at the cutoff {record['Cutoff'] * 1000:g} nm",
                  f"- local slope median of the sheet-edge MA density {ma_tail.f3(summary['LocalSlopeMedian'])}",
                  "", "Rule: " + record["Rule"]]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--postpro", type=Path, required=True)
    parser.add_argument("--kind", choices=KINDS, required=True)
    parser.add_argument("--ma-interface", type=int, default=3)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-md", type=Path)
    args = parser.parse_args(argv)
    record = analyze(args.postpro, kind=args.kind, ma_interface=args.ma_interface)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(record, indent=2) + "\n")
    if args.out_md:
        args.out_md.write_text(markdown(record))
    print(markdown(record), end="")


if __name__ == "__main__":
    main()
