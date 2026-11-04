#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Signature-only process library from a version-2 preflight manifest.

One model per COUPON of ``Identification.Features`` (decision 85(1), the library contract):
features of one type whose signature TOPOLOGY (the signature with every continuous parameter
removed) is equal and whose parameters agree within the recorded tolerance
(``SignatureParameterToleranceOverR`` = 1e-3 R for offsets / separations / radii,
``SignatureAngleToleranceDegrees`` = 1e-2 deg) are grouped by single linkage into one coupon
whose ``Signature`` is the group's representative (every parameter the midpoint of its range,
in the orientation nearest to the lexicographically smallest member — a set function; the same
rule as ``RepresentativeSignature`` in surfaceresponseidentification.cpp), with placeholder
matrix paths: enough for the patch dry run (``palace --surface-response-preflight`` builds the
patches without reading any matrix; the matcher accepts a model within the tolerance and takes
the nearest), never for a field solve. A dry run with this library must patch every feature of
the manifest exactly once, i.e. cover the whole perimeter minus the recorded exclusions.

    python3 -m surface_response_identification.signature_library MANIFEST OUTPUT.json
"""

import argparse
import json
import os
import sys

# Round-3 class (10) (decisions 492 / 493 / 510): the arc chording below feeds the content-hashed
# coupon sources, so its trigonometry is the correctly rounded deterministic_math (platform-
# independent by definition), not libm. Imported both as a package module and as a plain
# module on sys.path (the two ways this module is loaded).
try:
    from . import deterministic_math
except ImportError:
    import deterministic_math

# Feature types the library cannot model (an UnclassifiedParallelPair is an exclusion
# described as a feature so that the partition stays exact).
UNMODELLED_TYPES = ("UnclassifiedParallelPair", "CurvedUnclassifiedParallelPair")

PAIR_TYPES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip", "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip")
LONGITUDINAL_TYPES = ("IsolatedEdge", "CurvedEdge", "ParallelEdgeCluster") + PAIR_TYPES
VERTEX_TYPES = ("ConvexCorner", "ConcaveCorner", "Endpoint", "Junction")


def conductor_count(feature):
    """Canonical conductor labels of a feature (1 for single-conductor features). A
    SpatialEdgeCluster signature of contract v3 labels its conductors by first appearance over
    the sorted Portions THEN the sorted Context (a foreign conductor touching no claim gets the
    next label), so the count runs over both lists."""
    signature = feature["Signature"]
    if feature["Type"] == "SpatialEdgeCluster":
        return max(int(p["Conductor"]) for p in signature["Portions"] + signature.get("Context", []))
    if feature["Type"] == "ParallelEdgeCluster" or feature["Type"] in PAIR_TYPES:
        return max(int(e["Conductor"]) for e in signature["Edges"])
    return 1


# Continuous signature entries (surfaceresponseidentification.cpp IsLengthParameter /
# IsAngleParameter) and the recorded tolerances (Identification.Conventions).
LENGTH_KEYS = ("OffsetOverR", "SeparationOverR", "RadiusOverR", "CornerRadiusOverR")
ANGLE_KEYS = ("AngleDegrees", "ArmAnglesDegrees")
PARAMETER_TOLERANCE_OVER_R = 1.0e-3
ANGLE_TOLERANCE_DEGREES = 1.0e-2
LENGTH_QUANTUM_OVER_R = 1.0e-6
ANGLE_QUANTUM_DEGREES = 1.0e-6
# Block (b) DESIGN section 4 (decision 303; surfaceresponseidentification.hpp
# kClusterQuantumNearMatchMaxQuanta): two SpatialEdgeCluster signatures of one topology key
# (every number of Portions / Context [].P / Arc / Gap, Box and Vertices[].P / TurnDegrees
# nulled, entry order preserved) whose numbers agree within this many signature quanta (1e-6 R
# / 1e-6 deg) are ONE geometry at the grid: the matcher resolves the feature to the model
# (Match.QuantumNearMatch recorded), the library groups them into one coupon keyed by the
# lexicographically smallest member (the others listed under NearKeys).
CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA = 4
# Decisions 287 (b) / 288 (2) / 317 MINOR-1 (kClusterQuantumInclusiveMargin): every quantum
# threshold is read HALF-QUANTUM INCLUSIVE - an on-grid difference of k quanta computes to
# k +- 1e-9 in floating point, so a feature matches iff max |delta| <= k + 1/2 quanta and two
# models (or two span-cap allowances) are one geometry iff max |delta| <= 2 k + 1/2.
CLUSTER_QUANTUM_INCLUSIVE_MARGIN = 0.5
CLUSTER_LENGTH_KEYS = ("P", "Arc", "Gap", "Box")
CLUSTER_ANGLE_KEYS = ("TurnDegrees",)


def is_cluster_signature(signature):
    return isinstance(signature, dict) and signature.get("Type") == "SpatialEdgeCluster"


def split_cluster_parameters(signature):
    """(topology key, lengths, angles, length paths, angle paths) of a SpatialEdgeCluster
    signature: the C++ SplitSignatureParameters cluster branch (traversal order = sorted keys
    per object, list order preserved; the paths name every number, e.g. Portions[3].P[1])."""
    lengths, angles, length_paths, angle_paths = [], [], [], []

    def walk(node, path):
        if isinstance(node, dict):
            out = {}
            for key in sorted(node):
                value = node[key]
                entry = f"{path}.{key}" if path else key
                if key in CLUSTER_LENGTH_KEYS or key in CLUSTER_ANGLE_KEYS:
                    target, paths = (lengths, length_paths) if key in CLUSTER_LENGTH_KEYS else (angles, angle_paths)
                    if isinstance(value, list):
                        for k, v in enumerate(value):
                            target.append(float(v))
                            paths.append(f"{entry}[{k}]")
                    else:
                        target.append(float(value))
                        paths.append(entry)
                    out[key] = None
                else:
                    out[key] = walk(value, entry)
            return out
        if isinstance(node, list):
            return [walk(v, f"{path}[{k}]") for k, v in enumerate(node)]
        return node
    topology = json.dumps(walk(signature, ""), sort_keys=True)
    return topology, lengths, angles, length_paths, angle_paths


def cluster_quantum_difference(a, b):
    """(max |delta| in quanta, differing paths) of two SpatialEdgeCluster signatures of one
    topology key (the C++ ClusterSignatureQuantumDifference), or None when the topology keys
    differ (a permuted entry order included)."""
    ta, la, aa, lpa, apa = split_cluster_parameters(a)
    tb, lb, ab, _, _ = split_cluster_parameters(b)
    if ta != tb or len(la) != len(lb) or len(aa) != len(ab):
        return None
    worst, paths = 0.0, []
    for x, y, path in zip(la, lb, lpa):
        quanta = abs(x - y) / LENGTH_QUANTUM_OVER_R
        if quanta > 0.5:
            paths.append(path)
        worst = max(worst, quanta)
    for x, y, path in zip(aa, ab, apa):
        quanta = abs(x - y) / ANGLE_QUANTUM_DEGREES
        if quanta > 0.5:
            paths.append(path)
        worst = max(worst, quanta)
    return worst, paths


def split_parameters(signature):
    """(topology key, lengths, angles) of a signature: the continuous entries replaced by null in
    the serialised topology key (SpatialEdgeCluster: the quantum near-match branch, every number
    of P / Arc / Gap / Box / TurnDegrees a parameter on its quantum)."""
    lengths, angles = [], []
    if is_cluster_signature(signature):
        topology, lengths, angles, _, _ = split_cluster_parameters(signature)
        return topology, lengths, angles

    def walk(node):
        if isinstance(node, dict):
            out = {}
            for key in sorted(node):
                value = node[key]
                if key in LENGTH_KEYS or key in ANGLE_KEYS:
                    target = lengths if key in LENGTH_KEYS else angles
                    target.extend(float(v) for v in (value if isinstance(value, list) else [value]))
                    out[key] = None
                else:
                    out[key] = walk(value)
            return out
        if isinstance(node, list):
            return [walk(v) for v in node]
        return node
    return json.dumps(walk(signature), sort_keys=True), lengths, angles


def substitute_parameters(signature, lengths, angles):
    """The signature with its continuous entries replaced in the same traversal order (a
    SpatialEdgeCluster unchanged: the representative of a near-matching group is its
    lexicographically smallest member, never a midpoint)."""
    if is_cluster_signature(signature):
        return signature
    lengths, angles = list(lengths), list(angles)

    def walk(node):
        if isinstance(node, dict):
            out = {}
            for key in sorted(node):
                value = node[key]
                if key in LENGTH_KEYS or key in ANGLE_KEYS:
                    source, quantum = (lengths, LENGTH_QUANTUM_OVER_R) if key in LENGTH_KEYS else (angles, ANGLE_QUANTUM_DEGREES)
                    if isinstance(value, list):
                        out[key] = [round(source.pop(0) / quantum) * quantum for _ in value]
                    else:
                        out[key] = round(source.pop(0) / quantum) * quantum
                else:
                    out[key] = walk(value)
            return out
        if isinstance(node, list):
            return [walk(v) for v in node]
        return node
    out = walk(signature)
    assert not lengths and not angles
    return out


def mirror_translational(signature):
    """The mirror orientation of an ``Edges`` signature (reversed, gap sides negated, offsets from
    the top, conductors relabelled by first appearance); other signatures unchanged."""
    if not isinstance(signature, dict) or not signature.get("Edges"):
        return signature
    edges = signature["Edges"]
    span = float(edges[-1]["OffsetOverR"])
    labels, out = {}, []
    for edge in reversed(edges):
        e = dict(edge)
        e["OffsetOverR"] = round((span - float(edge["OffsetOverR"])) / LENGTH_QUANTUM_OVER_R) * LENGTH_QUANTUM_OVER_R
        e["GapSide"] = -int(edge["GapSide"])
        e["Conductor"] = labels.setdefault(int(edge["Conductor"]), len(labels) + 1)
        out.append(e)
    mirror = dict(signature)
    mirror["Edges"] = out
    return mirror


def within_cluster_quantum_near_match(max_delta_quanta):
    """The half-quantum inclusive near-match rule: max |delta| <= k + 1/2 quanta."""
    return max_delta_quanta <= CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA + CLUSTER_QUANTUM_INCLUSIVE_MARGIN


def cluster_quantum_duplicate(max_delta_quanta):
    """Two cluster signatures are one geometry (refused as two models / two allowances) iff
    max |delta| <= 2 k + 1/2 quanta (half-quantum inclusive)."""
    return max_delta_quanta <= 2 * CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA + CLUSTER_QUANTUM_INCLUSIVE_MARGIN


def signature_deviation(a, b):
    """max |difference| / tolerance over the parameters (both orientations of b; the smaller;
    a SpatialEdgeCluster: max |delta| / ((CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA + 1/2) quanta),
    the half-quantum inclusive rule), or None when the topologies differ. Within tolerance
    iff <= 1."""
    if is_cluster_signature(a) or is_cluster_signature(b):
        difference = cluster_quantum_difference(a, b)
        if difference is None:
            return None
        return difference[0] / (CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA + CLUSTER_QUANTUM_INCLUSIVE_MARGIN)
    ta, la, aa = split_parameters(a)
    best = None
    for candidate in (b, mirror_translational(b)):
        tb, lb, ab = split_parameters(candidate)
        if ta != tb or len(la) != len(lb) or len(aa) != len(ab):
            continue
        d = max([abs(x - y) / PARAMETER_TOLERANCE_OVER_R for x, y in zip(la, lb)] + [abs(x - y) / ANGLE_TOLERANCE_DEGREES for x, y in zip(aa, ab)] + [0.0])
        best = d if best is None else min(best, d)
    return best


def representative_signature(signatures):
    """The coupon's signature of a group of instances: the lexicographically smallest instance
    with every parameter replaced by the midpoint of its range over the group (each instance in
    the orientation nearest to that lead)."""
    lead = min(signatures, key=lambda s: json.dumps(s, sort_keys=True))
    t_lead, l_lead, a_lead = split_parameters(lead)
    lo_l, hi_l, lo_a, hi_a = list(l_lead), list(l_lead), list(a_lead), list(a_lead)
    for signature in signatures:
        best = None
        for candidate in (signature, mirror_translational(signature)):
            t, lens, angs = split_parameters(candidate)
            if t != t_lead or len(lens) != len(l_lead) or len(angs) != len(a_lead):
                continue
            d = max([abs(x - y) for x, y in zip(lens, l_lead)] + [abs(x - y) for x, y in zip(angs, a_lead)] + [0.0])
            if best is None or d < best[0]:
                best = (d, lens, angs)
        assert best is not None, "representative over different topologies"
        lo_l = [min(x, y) for x, y in zip(lo_l, best[1])]
        hi_l = [max(x, y) for x, y in zip(hi_l, best[1])]
        lo_a = [min(x, y) for x, y in zip(lo_a, best[2])]
        hi_a = [max(x, y) for x, y in zip(hi_a, best[2])]
    return substitute_parameters(lead, [0.5 * (x + y) for x, y in zip(lo_l, hi_l)], [0.5 * (x + y) for x, y in zip(lo_a, hi_a)])


def group_features(features):
    """Coupons of a feature list: per (type, topology key incl. both orientations) the distinct
    signatures grouped by single linkage at the tolerance (order independent); returns a list of
    (representative signature, member features, parameter spread)."""
    bases = {}
    for feature in features:
        signature = feature["Signature"]
        keys = [feature["Type"] + "|" + split_parameters(s)[0] for s in (signature, mirror_translational(signature))]
        bases.setdefault(min(keys), {}).setdefault(json.dumps(signature, sort_keys=True), (signature, []))[1].append(feature)
    groups = []
    for base_key in sorted(bases):
        distinct = [bases[base_key][k] for k in sorted(bases[base_key])]
        parent = list(range(len(distinct)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i in range(len(distinct)):
            for j in range(i + 1, len(distinct)):
                d = signature_deviation(distinct[i][0], distinct[j][0])
                if d is not None and d <= 1.0:
                    parent[find(i)] = find(j)
        members = {}
        for i in range(len(distinct)):
            members.setdefault(find(i), []).append(i)
        for root in sorted(members, key=lambda r: json.dumps(distinct[r][0], sort_keys=True)):
            signatures = [distinct[i][0] for i in members[root]]
            representative = representative_signature(signatures)
            spread = max((signature_deviation(representative, s) or 0.0) for s in signatures)
            groups.append((representative, [f for i in members[root] for f in distinct[i][1]], spread))
    return groups


# Canonical chording of the arc portions of a SpatialEdgeCluster signature (option A,
# Identification.Conventions ClusterArcChordStepDegrees / ClusterArcChordMaxLengthOverR): the
# coupon builder receives the arcs and chords them at this step, so that the coupon geometry
# is a function of the signature and not of the device mesh.
CLUSTER_ARC_CHORD_STEP_DEGREES = 5.0
CLUSTER_ARC_CHORD_MAX_LENGTH_OVER_R = 0.25


def cluster_plan_view_edges(signature, radius, step_degrees=CLUSTER_ARC_CHORD_STEP_DEGREES, max_chord_over_R=CLUSTER_ARC_CHORD_MAX_LENGTH_OVER_R, include_context=False):
    """Straight plan-view edges of a SpatialEdgeCluster signature in its canonical frame
    (mesh units = the signature's units of R times ``radius``): every straight portion as one
    edge ``{"P0", "P1", "Gap", "Conductor", "Interfaces", "Law", "Portion"}``; every arc
    portion (``Arc`` = centre + midpoint, ``GapRadial``) chorded into n equal chords, n =
    max(ceil(sweep / step), ceil(arc length / (max_chord_over_R R)), 1), each chord's gap
    direction the arc's radial direction at the chord's middle (times ``GapRadial``). A closed
    circle (equal ends) is chorded over 2 pi. With ``include_context`` the ``Context`` entries of
    a contract-v3 signature (the device plan clipped to the ``Box``: the continuation chains,
    ``Chain`` true, and the foreign edges) follow the claims in the same encoding, each carrying
    ``"Context": True`` and its ``Chain`` flag (``Portion`` indexes the Context list).

    Every float of the result is produced by CPython scalar arithmetic on the serialised
    coordinates and by deterministic_math's correctly rounded atan2 / cos / sin (round-3 class
    (10)): a straight portion involves no trigonometry and is returned bitwise as before."""
    import math
    edges = []
    entries = [(index, portion, False) for index, portion in enumerate(signature["Portions"])]
    if include_context:
        entries += [(index, portion, True) for index, portion in enumerate(signature.get("Context", []))]
    for index, portion, context in entries:
        p = [float(v) * radius for v in portion["P"]]
        a, b = (p[0], p[1]), (p[2], p[3])
        common = {"Conductor": int(portion["Conductor"]), "Interfaces": portion.get("Interfaces", []), "Law": portion.get("Law"), "Portion": index}
        if context:
            common.update({"Context": True, "Chain": bool(portion.get("Chain", False))})
        if "Arc" not in portion:
            edges.append(dict(common, P0=a, P1=b, Gap=(float(portion["Gap"][0]), float(portion["Gap"][1]))))
            continue
        arc = [float(v) * radius for v in portion["Arc"]]
        c, m = (arc[0], arc[1]), (arc[2], arc[3])
        r = math.hypot(a[0] - c[0], a[1] - c[1])
        angle = lambda q: deterministic_math.atan2(q[1] - c[1], q[0] - c[0])
        ta, tb, tm = angle(a), angle(b), angle(m)
        closed = math.hypot(a[0] - b[0], a[1] - b[1]) <= 1.0e-9 * max(r, 1.0)
        if closed:
            sweep = 2.0 * math.pi
        else:
            # The arc from a to b through m: the counterclockwise sweep from a to b holds m, or
            # the clockwise one does.
            ccw = (tb - ta) % (2.0 * math.pi)
            if ((tm - ta) % (2.0 * math.pi)) <= ccw + 1.0e-12:
                sweep = ccw
            else:
                sweep = ccw - 2.0 * math.pi
        n = max(int(math.ceil(abs(sweep) / math.radians(step_degrees) - 1.0e-9)), int(math.ceil(r * abs(sweep) / (max_chord_over_R * radius) - 1.0e-9)), 1)
        sign = int(portion["GapRadial"])
        for k in range(n):
            t0, t1 = ta + sweep * k / n, ta + sweep * (k + 1) / n
            q0 = (c[0] + r * deterministic_math.cos(t0), c[1] + r * deterministic_math.sin(t0))
            q1 = (c[0] + r * deterministic_math.cos(t1), c[1] + r * deterministic_math.sin(t1))
            tmid = 0.5 * (t0 + t1)
            edges.append(dict(common, P0=q0, P1=q1, Gap=(sign * deterministic_math.cos(tmid), sign * deterministic_math.sin(tmid)),
                              Chord=k, Chords=n))
    return edges


# Spatial-support contract v3 (USER decision 281, supervisor decision 282; surfaceresponse-
# identification.hpp kSupport*): the claims-derived support box of a SpatialEdgeCluster in
# its canonical frame, units of R — the coupon generator's coupon_bounds / edge_rows /
# extended_interval rule (claim-cut ends lengthened to >= R, every row end at or beyond R
# continued by 2R, rows widened by R on both sides, padding R) evaluated on the serialised
# signature, arcs chorded as cluster_plan_view_edges chords them; the same numbers as the C++
# SupportBoxFromSignature (its ClaimsBox record), quantised on the 1e-6 R grid.
SUPPORT_CONTINUATION_OVER_R = 2.0
SUPPORT_PADDING_OVER_R = 1.0
SUPPORT_END_COINCIDENCE_OVER_R = 1.0e-5
SUPPORT_FACE_SNAP_OVER_R = 1.0e-3
SUPPORT_FACE_CLEARANCE_OVER_R = 0.25
SUPPORT_FACE_GROWTH_STEP_OVER_R = 0.25
SUPPORT_FACE_GROWTH_MAX_STEPS = 12
SUPPORT_SPAN_CAP_OVER_R = 16.0
# Decision 287 (b): the C++ face rules (T1 snap, T2 clearance and crossing separation) compare
# on this grid — a value within half a quantum of a threshold is AT it and takes the rule's
# inclusive side. The face rules themselves are not mirrored here (the C++ record is read).
SUPPORT_COMPARISON_QUANTUM_OVER_R = 1.0e-6


def cluster_support_box(signature):
    """[x0, y0, x1, y1] in units of R of the claims-derived support box (rule B2) of a
    SpatialEdgeCluster signature in its own frame (the ``Box`` of a contract-2 signature, the
    ``ClaimsBox`` of the manifest record before any T3 growth)."""
    import math
    edges = cluster_plan_view_edges(signature, 1.0)
    vertices = [(float(v["P"][0]), float(v["P"][1])) for v in signature.get("Vertices", [])]

    def near(a, b, tolerance):
        return math.hypot(a[0] - b[0], a[1] - b[1]) <= tolerance

    def connected(i, point):
        if any(near(point, v, SUPPORT_END_COINCIDENCE_OVER_R) for v in vertices):
            return True
        return any(near(point, q, SUPPORT_END_COINCIDENCE_OVER_R) for j, e in enumerate(edges) if j != i for q in (e["P0"], e["P1"]))

    x0 = y0 = math.inf
    x1 = y1 = -math.inf
    for i, edge in enumerate(edges):
        gx, gy = edge["Gap"]
        norm = math.hypot(gx, gy)
        gap = (gx / norm, gy / norm)
        tangent = (gap[1], -gap[0])
        p0, p1 = edge["P0"], edge["P1"]
        length = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
        midpoint = (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]))
        begin_free, end_free = not connected(i, p0), not connected(i, p1)
        forward_is_p1 = (p1[0] - midpoint[0]) * tangent[0] + (p1[1] - midpoint[1]) * tangent[1] > 0.0
        end_is_free = end_free if forward_is_p1 else begin_free
        begin_is_free = begin_free if forward_is_p1 else end_free
        half = 0.5 * length
        begin = -max(half, 1.0) if begin_is_free else -half
        end = max(half, 1.0) if end_is_free else half
        if begin <= -1.0 + 1.0e-10:
            begin -= SUPPORT_CONTINUATION_OVER_R
        if end >= 1.0 - 1.0e-10:
            end += SUPPORT_CONTINUATION_OVER_R
        for coordinate in (begin, end):
            for side in (-1.0, 1.0):
                px = midpoint[0] + coordinate * tangent[0] + side * gap[0]
                py = midpoint[1] + coordinate * tangent[1] + side * gap[1]
                x0, y0, x1, y1 = min(x0, px), min(y0, py), max(x1, px), max(y1, py)

    def quantize(v):
        q = round(v / LENGTH_QUANTUM_OVER_R) * LENGTH_QUANTUM_OVER_R
        return 0.0 if q == 0.0 else q
    return [quantize(x0 - SUPPORT_PADDING_OVER_R), quantize(y0 - SUPPORT_PADDING_OVER_R), quantize(x1 + SUPPORT_PADDING_OVER_R), quantize(y1 + SUPPORT_PADDING_OVER_R)]


# --------------------------------------------------------------------------------------------------------------------
# The C++ writer's float print (fix round, OI-SIGNATURE-FLOAT-PRINT, decisions 1157 Q4 / 1228 M2 / 1233 (4)): palace hashes
# the KeyText `nlohmann::json::dump()` prints (SignatureKeyAndHash), whose doubles come from nlohmann/json 3.12.0
# detail/conversions/to_chars.hpp = the GRISU2 algorithm (Loitsch 2010) with nlohmann's cached powers and format_buffer. Grisu2
# is NOT the shortest round-trip print (Python's repr): 2.9192109999999998 (nlohmann) vs 2.919211 (repr) on W06's manifest. The
# functions below are a faithful port (64-bit integer arithmetic; the same table, rounding and formatting), so that every
# Python-side serialisation of a Signature / KeyText (cpp_json_dump) is BIT-EQUAL to the C++ writer: 0 device re-keys and 0
# library-signature changes by construction (AMENDMENT 1228). The equivalence fixture (test_signature_library.py) re-prints
# every KeyText of the manifests of record.
# --------------------------------------------------------------------------------------------------------------------

_U64 = (1 << 64) - 1
_GRISU_ALPHA = -60
_GRISU_GAMMA = -32
_CACHED_POWERS_MIN_DEC_EXP = -300
_CACHED_POWERS_DEC_STEP = 8
# nlohmann dtoa_impl::kCachedPowers (79 entries: f, e, k) -- 10^k ~= f x 2^e, f normalized.
_CACHED_POWERS = (
    (0xAB70FE17C79AC6CA, -1060, -300), (0xFF77B1FCBEBCDC4F, -1034, -292), (0xBE5691EF416BD60C, -1007, -284), (0x8DD01FAD907FFC3C, -980, -276),
    (0xD3515C2831559A83, -954, -268), (0x9D71AC8FADA6C9B5, -927, -260), (0xEA9C227723EE8BCB, -901, -252), (0xAECC49914078536D, -874, -244),
    (0x823C12795DB6CE57, -847, -236), (0xC21094364DFB5637, -821, -228), (0x9096EA6F3848984F, -794, -220), (0xD77485CB25823AC7, -768, -212),
    (0xA086CFCD97BF97F4, -741, -204), (0xEF340A98172AACE5, -715, -196), (0xB23867FB2A35B28E, -688, -188), (0x84C8D4DFD2C63F3B, -661, -180),
    (0xC5DD44271AD3CDBA, -635, -172), (0x936B9FCEBB25C996, -608, -164), (0xDBAC6C247D62A584, -582, -156), (0xA3AB66580D5FDAF6, -555, -148),
    (0xF3E2F893DEC3F126, -529, -140), (0xB5B5ADA8AAFF80B8, -502, -132), (0x87625F056C7C4A8B, -475, -124), (0xC9BCFF6034C13053, -449, -116),
    (0x964E858C91BA2655, -422, -108), (0xDFF9772470297EBD, -396, -100), (0xA6DFBD9FB8E5B88F, -369, -92), (0xF8A95FCF88747D94, -343, -84),
    (0xB94470938FA89BCF, -316, -76), (0x8A08F0F8BF0F156B, -289, -68), (0xCDB02555653131B6, -263, -60), (0x993FE2C6D07B7FAC, -236, -52),
    (0xE45C10C42A2B3B06, -210, -44), (0xAA242499697392D3, -183, -36), (0xFD87B5F28300CA0E, -157, -28), (0xBCE5086492111AEB, -130, -20),
    (0x8CBCCC096F5088CC, -103, -12), (0xD1B71758E219652C, -77, -4), (0x9C40000000000000, -50, 4), (0xE8D4A51000000000, -24, 12),
    (0xAD78EBC5AC620000, 3, 20), (0x813F3978F8940984, 30, 28), (0xC097CE7BC90715B3, 56, 36), (0x8F7E32CE7BEA5C70, 83, 44),
    (0xD5D238A4ABE98068, 109, 52), (0x9F4F2726179A2245, 136, 60), (0xED63A231D4C4FB27, 162, 68), (0xB0DE65388CC8ADA8, 189, 76),
    (0x83C7088E1AAB65DB, 216, 84), (0xC45D1DF942711D9A, 242, 92), (0x924D692CA61BE758, 269, 100), (0xDA01EE641A708DEA, 295, 108),
    (0xA26DA3999AEF774A, 322, 116), (0xF209787BB47D6B85, 348, 124), (0xB454E4A179DD1877, 375, 132), (0x865B86925B9BC5C2, 402, 140),
    (0xC83553C5C8965D3D, 428, 148), (0x952AB45CFA97A0B3, 455, 156), (0xDE469FBD99A05FE3, 481, 164), (0xA59BC234DB398C25, 508, 172),
    (0xF6C69A72A3989F5C, 534, 180), (0xB7DCBF5354E9BECE, 561, 188), (0x88FCF317F22241E2, 588, 196), (0xCC20CE9BD35C78A5, 614, 204),
    (0x98165AF37B2153DF, 641, 212), (0xE2A0B5DC971F303A, 667, 220), (0xA8D9D1535CE3B396, 694, 228), (0xFB9B7CD9A4A7443C, 720, 236),
    (0xBB764C4CA7A44410, 747, 244), (0x8BAB8EEFB6409C1A, 774, 252), (0xD01FEF10A657842C, 800, 260), (0x9B10A4E5E9913129, 827, 268),
    (0xE7109BFBA19C0C9D, 853, 276), (0xAC2820D9623BF429, 880, 284), (0x80444B5E7AA7CF85, 907, 292), (0xBF21E44003ACDD2D, 933, 300),
    (0x8E679C2F5E44FF8F, 960, 308), (0xD433179D9C8CB841, 986, 316), (0x9E19DB92B4E31BA9, 1013, 324),
)


def _diyfp_mul(xf, xe, yf, ye):
    """diyfp::mul: the upper 64 bits of the 128-bit product, rounded (ties up) at bit 63 of the low half."""
    u_lo, u_hi = xf & 0xFFFFFFFF, xf >> 32
    v_lo, v_hi = yf & 0xFFFFFFFF, yf >> 32
    p0, p1, p2, p3 = u_lo * v_lo, u_lo * v_hi, u_hi * v_lo, u_hi * v_hi
    q = (p0 >> 32) + (p1 & 0xFFFFFFFF) + (p2 & 0xFFFFFFFF)
    q += 1 << 31
    h = p3 + (p2 >> 32) + (p1 >> 32) + (q >> 32)
    return h & _U64, xe + ye + 64


def _diyfp_normalize(f, e):
    while (f >> 63) == 0:
        f <<= 1
        e -= 1
    return f, e


def _grisu2_boundaries(value):
    """compute_boundaries for a positive finite double -> (v, m_minus, m_plus) as (f, e) pairs, normalized, m+/m- sharing e."""
    import struct
    bits = struct.unpack("<Q", struct.pack("<d", value))[0]
    precision = 53
    bias = 1024 - 1 + (precision - 1)
    hidden = 1 << (precision - 1)
    E = bits >> (precision - 1)
    F = bits & (hidden - 1)
    if E == 0:
        vf, ve = F, 1 - bias
    else:
        vf, ve = F + hidden, E - bias
    lower_closer = F == 0 and E > 1
    mpf, mpe = 2 * vf + 1, ve - 1
    if lower_closer:
        mmf, mme = 4 * vf - 1, ve - 2
    else:
        mmf, mme = 2 * vf - 1, ve - 1
    wpf, wpe = _diyfp_normalize(mpf, mpe)
    delta = mme - wpe
    wmf = mmf << delta
    v = _diyfp_normalize(vf, ve)
    return v, (wmf, wpe), (wpf, wpe)


def _cached_power_for_binary_exponent(e):
    f = _GRISU_ALPHA - e - 1
    k = ((f * 78913) >> 18) + (1 if f > 0 else 0)
    # C++ integer division truncates toward zero; f * 78913 / 2^18 with f < 0 must round toward zero, not floor
    if f * 78913 < 0:
        k = -((-f * 78913) >> 18) + (1 if f > 0 else 0)
    index = (-_CACHED_POWERS_MIN_DEC_EXP + k + (_CACHED_POWERS_DEC_STEP - 1)) // _CACHED_POWERS_DEC_STEP
    return _CACHED_POWERS[index]


def _grisu2_round(buf, dist, delta, rest, ten_k):
    while rest < dist and delta - rest >= ten_k and (rest + ten_k < dist or dist - rest > rest + ten_k - dist):
        buf[-1] -= 1
        rest += ten_k


def _grisu2_digit_gen(buf, decimal_exponent, m_minus, w, m_plus):
    """grisu2_digit_gen: the decimal digits of V = buf x 10^decimal_exponent with M- <= V <= M+ (returns the exponent)."""
    delta = (m_plus[0] - m_minus[0]) & _U64
    dist = (m_plus[0] - w[0]) & _U64
    e = m_plus[1]
    one_f = 1 << -e
    p1 = m_plus[0] >> -e
    p2 = m_plus[0] & (one_f - 1)
    pow10 = 1
    k = 1
    for digits, p in ((10, 1000000000), (9, 100000000), (8, 10000000), (7, 1000000), (6, 100000), (5, 10000), (4, 1000), (3, 100), (2, 10)):
        if p1 >= p:
            pow10, k = p, digits
            break
    n = k
    while n > 0:
        d, r = divmod(p1, pow10)
        buf.append(d)
        p1 = r
        n -= 1
        rest = (p1 << -e) + p2
        if rest <= delta:
            decimal_exponent += n
            _grisu2_round(buf, dist, delta, rest, pow10 << -e)
            return decimal_exponent
        pow10 //= 10
    m = 0
    while True:
        p2 *= 10
        d = p2 >> -e
        r = p2 & (one_f - 1)
        buf.append(d)
        p2 = r
        m += 1
        delta *= 10
        dist *= 10
        if p2 <= delta:
            break
    decimal_exponent -= m
    _grisu2_round(buf, dist, delta, p2, one_f)
    return decimal_exponent


def _grisu2(value):
    """grisu2(double) -> (digits list, decimal_exponent) with value = digits x 10^decimal_exponent."""
    v, m_minus, m_plus = _grisu2_boundaries(value)
    cf, ce, ck = _cached_power_for_binary_exponent(m_plus[1])
    w = _diyfp_mul(v[0], v[1], cf, ce)
    w_minus = _diyfp_mul(m_minus[0], m_minus[1], cf, ce)
    w_plus = _diyfp_mul(m_plus[0], m_plus[1], cf, ce)
    big_minus = (w_minus[0] + 1, w_minus[1])
    big_plus = (w_plus[0] - 1, w_plus[1])
    buf = []
    exponent = _grisu2_digit_gen(buf, -ck, big_minus, w, big_plus)
    return buf, exponent


def _format_buffer(digits, decimal_exponent, min_exp=-4, max_exp=15):
    """nlohmann dtoa_impl::format_buffer (kMinExp -4, kMaxExp digits10 = 15): fixed notation in [10^min_exp, 10^max_exp), else
    the exponent form d.ddde+XX (at least two exponent digits)."""
    text = "".join(chr(48 + d) for d in digits)
    k = len(digits)
    n = k + decimal_exponent
    if k <= n <= max_exp:
        return text + "0" * (n - k) + ".0"
    if 0 < n <= max_exp:
        return text[:n] + "." + text[n:]
    if min_exp < n <= 0:
        return "0." + "0" * (-n) + text
    mantissa = text if k == 1 else text[0] + "." + text[1:]
    e = n - 1
    sign = "-" if e < 0 else "+"
    e = abs(e)
    return mantissa + "e" + sign + (f"{e:02d}" if e < 100 else f"{e:d}")


def cpp_to_chars(value):
    """nlohmann::detail::to_chars(double): the text `json::dump()` writes for a finite double (the C++ writer's print)."""
    import math
    if not isinstance(value, float):
        value = float(value)
    if not math.isfinite(value):
        return "null"
    sign = "-" if math.copysign(1.0, value) < 0 else ""
    value = abs(value)
    if value == 0.0:
        return sign + "0.0"
    digits, exponent = _grisu2(value)
    return sign + _format_buffer(digits, exponent)


def cpp_json_string(text):
    """nlohmann's string escaping with ensure_ascii = false: the two-character escapes, \\u00XX for the other control characters,
    every other code point written as is (UTF-8)."""
    out = ['"']
    for ch in text:
        code = ord(ch)
        if ch == '"':
            out.append('\\"')
        elif ch == "\\":
            out.append("\\\\")
        elif ch == "\b":
            out.append("\\b")
        elif ch == "\f":
            out.append("\\f")
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\r":
            out.append("\\r")
        elif ch == "\t":
            out.append("\\t")
        elif code <= 0x1F:
            out.append(f"\\u{code:04x}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def cpp_json_dump(value):
    """`nlohmann::json::dump()` (no indent) of a JSON value parsed by Python's json: objects with their keys in std::map order
    (byte order of the UTF-8 key; Python's sort_keys order for ASCII keys), ',' / ':' without spaces, integers as integers,
    doubles by cpp_to_chars, true / false / null. A Python float that is integral prints "N.0" as nlohmann does for a double."""
    if value is True:
        return "true"
    if value is False:
        return "false"
    if value is None:
        return "null"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return cpp_to_chars(value)
    if isinstance(value, str):
        return cpp_json_string(value)
    if isinstance(value, dict):
        items = sorted(value.items(), key=lambda kv: kv[0].encode("utf-8"))
        return "{" + ",".join(cpp_json_string(k) + ":" + cpp_json_dump(v) for k, v in items) + "}"
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(cpp_json_dump(v) for v in value) + "]"
    raise TypeError(f"cpp_json_dump: unsupported value {type(value).__name__}")


def signature_key_text(signature, feature_type=None):
    """The exact text palace hashes into a feature's key (SignatureKeyAndHash: the Signature with "Type" set, dumped by the C++
    writer); `feature_type` None keeps the Signature's own Type."""
    body = dict(signature)
    if feature_type is not None:
        body["Type"] = feature_type
    return cpp_json_dump(body)


def signature_hash(signature):
    """sha256 of the C++ writer's print of the signature (cpp_json_dump; decision 1233 (4): the legacy float-hash path now equals
    the C++ hash where Python's repr and nlohmann's Grisu2 differed)."""
    import hashlib
    return hashlib.sha256(cpp_json_dump(signature).encode("utf-8")).hexdigest()


class DuplicateSignatureError(ValueError):
    """Two models of one library carry one canonical Signature (OI-LIBRARY-DUPLICATE-SIGNATURE; fail closed, by name)."""


def duplicate_signature_models(models):
    """{signature hash: [model Names]} for every canonical Signature (json sort_keys) that two or more models of one
    library carry; a model without a Signature is keyed by its Name alone and never collides here."""
    by_signature = {}
    for model in models:
        signature = model.get("Signature")
        if signature is None:
            continue
        by_signature.setdefault(signature_hash(signature), []).append(str(model.get("Name")))
    return {key: names for key, names in by_signature.items() if len(names) > 1}


def refuse_duplicate_signatures(models, *, library="library"):
    """OI-LIBRARY-DUPLICATE-SIGNATURE (decisions 1103 / 1126 (5) / 1128 Q4): a library holds AT MOST ONE model per signature,
    whatever its form. Two models of one key tie on (deviation, rank) in the C++ matcher and the lexicographically smaller
    Name serves the key SILENTLY (surfaceresponseoperator.cpp; proved on S1 by the lane-C confirmation): RG-7's one-reference
    and -pieces3 forms of a670c49a86d8 live in separate libraries. Every library ASSEMBLY tool calls this on the models it is
    about to write and raises DuplicateSignatureError naming the models; returns the models unchanged otherwise. The refusal is
    on EXACT canonical equality of the Signature (json sort_keys): two Signatures that differ in any value, however close (the
    one-quantum PEC pair of DESIGN A.5 (b), OffsetOverR 1e-6 apart), are two keys here; their tie in the matcher's
    tolerance-based exact match is the C++ load-time refusal's (decision 1103 (b), post-freeze) - no tolerance is applied here."""
    duplicates = duplicate_signature_models(models)
    if duplicates:
        raise DuplicateSignatureError(
            f"{library}: {len(duplicates)} signature(s) carried by more than one model (OI-LIBRARY-DUPLICATE-SIGNATURE, decision "
            "1103: at most one model per signature per library): "
            + "; ".join(f"{key[:12]} <- {', '.join(sorted(names))}" for key, names in sorted(duplicates.items())))
    return models


class SignatureKeyError(ValueError):
    """A requirement's key does not agree with its Signature / KeyText (fail closed, by name)."""


KEY_PATH_KEYTEXT = "KeyText"
KEY_PATH_LEGACY = "LegacyFloatHash"
KEY_PREFIX_MIN_HEX = 12


def verify_signature_key(signature, requirement_key, key_text=None, *, name="requirement"):
    """The key of a Signature-bearing requirement, verified against its Signature (decision 605
    (2) F-4 (c); the 438 (3) rule): returns {"Path", "KeyText", "SignatureHash", "RequirementKey",
    "Equal": True} or raises SignatureKeyError.

    - ``key_text`` present (a manifest written by a binary that records KeyText, the exact text
      palace hashed, Type included): sha256(key_text) must equal the requirement key (or carry
      it as a >= 12-hex prefix) AND json.loads(key_text) must equal ``signature`` (parsed-dict
      equality; no float is re-serialised, so nlohmann's Grisu2 lexemes never matter).
    - ``key_text`` None (a record-era manifest of the older binaries): the legacy float hash
      ``signature_hash(signature)`` is compared with the key, unchanged behaviour.
    - ``requirement_key`` None: no digest comparison (there is no key to compare); the parse
      check still runs and the record's SignatureHash is the digest of the text (or the legacy
      hash). Palace always writes Hash beside KeyText, so this is the "no key" caller's path
      (e.g. a generated basis), not a manifest record's.
    """
    import hashlib
    key = None if requirement_key is None else str(requirement_key)

    def matches(digest):
        return key is None or digest == key or (len(key) >= KEY_PREFIX_MIN_HEX and digest.startswith(key))

    if key_text is None:
        digest = signature_hash(signature)
        if not matches(digest):
            raise SignatureKeyError(f"{name}: signature_hash(Signature) {digest[:16]}… != the requirement key {key[:16]}… "
                                    f"(legacy float-hash path: the record carries no KeyText)")
        return {"Path": KEY_PATH_LEGACY, "KeyText": None, "SignatureHash": digest, "RequirementKey": key, "Equal": True}
    if not isinstance(key_text, str) or not key_text:
        raise SignatureKeyError(f"{name}: KeyText must be a non-empty string, not {type(key_text).__name__}")
    digest = hashlib.sha256(key_text.encode()).hexdigest()
    if not matches(digest):
        raise SignatureKeyError(f"{name}: sha256(KeyText) {digest[:16]}… != the requirement key {key[:16]}…")
    try:
        parsed = json.loads(key_text)
    except ValueError as error:
        raise SignatureKeyError(f"{name}: KeyText is not JSON ({error})") from error
    if parsed != signature:
        raise SignatureKeyError(f"{name}: json.loads(KeyText) != Signature (the key text names another signature)")
    return {"Path": KEY_PATH_KEYTEXT, "KeyText": key_text, "SignatureHash": digest, "RequirementKey": key, "Equal": True}


def feature_signature_key(feature):
    """The key of a manifest feature / requirement record: its verified Hash when the record
    carries KeyText (verify_signature_key, fail closed), else the legacy float hash of its
    Signature (a record-era manifest)."""
    if feature.get("KeyText") is None:
        return signature_hash(feature["Signature"])
    return verify_signature_key(feature["Signature"], feature.get("Hash"), feature["KeyText"],
                                name=f"feature {feature.get('Id', feature.get('Hash'))}")["SignatureHash"]


def context_digest(signature):
    """sha256 of the serialised {"Box", "Context"} of a contract-3 SpatialEdgeCluster signature
    (surfaceresponseidentification SpatialSupportContextDigest; the manifest records it as
    SpatialSupport.ContextDigest); '' for a claims-only signature. A Python float re-dump of
    the C++-hashed text (the same class as the pre-KeyText signature_hash): the manifest records
    no ContextText, so legacy_contract_alias cross-checks this digest against the recorded one
    and fails closed by name when a double's lexeme differs."""
    import hashlib
    if "Box" not in signature:
        return ""
    context = {"Box": signature["Box"], "Context": signature.get("Context", [])}
    return hashlib.sha256(json.dumps(context, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def legacy_contract_alias(feature, reason):
    """The library-side legacy-contract alias (USER decision 283) a legacy SpatialEdgeCluster
    model lists under ``LegacyContractAliases`` so that the contract-3 key of this manifest
    feature resolves to it: {"Key", "ContextDigest", "Reason", "Context", "ClaimsKey"}. The
    feature's claims-only key (the manifest's SpatialSupport.ClaimsKey, returned as ClaimsKey)
    must be the legacy model's key (checked by Palace at match time, fail closed); the digest
    is cross-checked against the manifest's record when present."""
    signature = feature["Signature"]
    if "Box" not in signature:
        raise ValueError(f"feature {feature.get('Id')} has a claims-only (contract-2) signature: nothing to alias")
    digest = context_digest(signature)
    recorded = (feature.get("SpatialSupport") or {}).get("ContextDigest")
    if recorded is not None and recorded != digest:
        raise ValueError(f"feature {feature.get('Id')}: the Python context digest {digest[:12]} differs from the manifest's {str(recorded)[:12]}")
    return {"Key": feature["Hash"], "ContextDigest": digest, "Reason": reason,
            "Context": {"Box": signature["Box"], "Context": signature.get("Context", [])},
            "ClaimsKey": (feature.get("SpatialSupport") or {}).get("ClaimsKey")}


def signature_model(feature_type, signature, radius, model_name, matrix_directory="signature-only-matrices"):
    """A geometry-only library model keyed by its canonical ``Signature`` (placeholder matrix
    paths; the version-1 geometry parameters the library reader validates are derived from
    the signature, the matching itself is by Signature). Shared by the signature-only library
    and the discovery placeholders (discover_surface_response_requirements.py)."""
    feature = {"Type": feature_type, "Signature": signature}
    model = {
        "Name": model_name,
        "Topology": feature_type,
        "Signature": signature,
        "FabricatedMatrix": f"{matrix_directory}/{model_name}-fabricated.csv",
        "ThinMatrix": f"{matrix_directory}/{model_name}-thin.csv",
        "BasisPoints": f"{matrix_directory}/{model_name}-basis-points.csv",
    }
    if feature_type in LONGITUDINAL_TYPES:
        # Longitudinal coupon depth: the matching radius (any positive depth; the dry run
        # records it with every patch weight).
        model["CouponDepth"] = radius
    if feature_type in PAIR_TYPES:
        model["Separation"] = float(signature["SeparationOverR"]) * radius
    elif feature_type in ("ConvexCorner", "ConcaveCorner"):
        model["Angle"] = float(signature["AngleDegrees"])
        model["CornerRadius"] = float(signature["CornerRadiusOverR"]) * radius
    elif feature_type == "Junction":
        angles, total = [], 0.0
        for difference in signature["ArmAnglesDegrees"]:
            angles.append(total)
            total += float(difference)
        model["ArmAngles"] = angles
    elif feature_type == "ParallelEdgeCluster":
        model["Edges"] = [{"Offset": float(e["OffsetOverR"]) * radius, "GapDirection": int(e["GapSide"]), "Conductor": int(e["Conductor"])} for e in signature["Edges"]]
    conductors = conductor_count(feature)
    if conductors > 1:
        # One reference per canonical conductor label; positions are placeholders.
        model["ConductorReferences"] = [[0.0, 0.0, -radius * (k + 1)] for k in range(conductors)]
    if feature_type != "SpatialEdgeCluster":
        law = json.loads(signature["Law"]) if "Law" in signature else {"Type": "PEC"}
        if law != {"Type": "PEC"}:
            model["BoundaryCondition"] = law
    return model


def build_signature_library(manifest, name="signature-only", matrix_directory="signature-only-matrices"):
    identification = manifest["Identification"]
    radius = float(identification["MatchingRadius"])
    models = []
    # An UnboxableFeature key (decision 282, "Unboxable": true) is a Missing placeholder no
    # builder makes and no library may serve: skipped, listed under Unboxable.
    unboxable = [f["Id"] for f in identification["Features"] if f.get("Signature", {}).get("Unboxable")]
    modelled = [f for f in identification["Features"]
                if f["Type"] not in UNMODELLED_TYPES and not f.get("Signature", {}).get("Unboxable")]
    for representative, members, spread in group_features(modelled):
        feature_type = members[0]["Type"]
        # The model key: a group whose representative IS its members' signature, on a manifest
        # carrying KeyText, takes the members' recorded key verified against that text (the C++
        # key, never a float re-dump; decision 605 (2) F-4 (c)); otherwise (a record-era manifest,
        # or a representative formed over several distinct signatures) the key is re-derived here
        # as before (legacy float hash).
        if members[0].get("KeyText") is not None and all(f["Signature"] == representative for f in members):
            representative_key = feature_signature_key(members[0])
        else:
            representative_key = signature_hash(representative)
        model_name = f"{feature_type}-{representative_key[:12]}"
        model = signature_model(feature_type, representative, radius, model_name, matrix_directory)
        model["Instances"] = len(members)
        model["DistinctSignatures"] = len({json.dumps(f["Signature"], sort_keys=True) for f in members})
        model["ParameterSpread"] = spread
        if feature_type == "SpatialEdgeCluster":
            near_keys = sorted({feature_signature_key(f) for f in members} - {representative_key})
            if near_keys:
                # The members' keys the matcher resolves to this model by the quantum
                # near-match (block (b) DESIGN section 4): recorded, never compared.
                model["NearKeys"] = near_keys
        models.append(model)
    library = {"Version": 2, "Name": name, "MatchingRadius": radius, "TraceLiftVersion": 2, "Models": models}
    if unboxable:
        library["UnboxableFeatures"] = unboxable
    return library


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("manifest")
    parser.add_argument("output")
    parser.add_argument("--name", default="signature-only")
    args = parser.parse_args(argv)
    with open(args.manifest) as source:
        manifest = json.load(source)
    library = build_signature_library(manifest, args.name)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as target:
        json.dump(library, target, indent=1)
    print(f"{len(library['Models'])} signature-only models -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
