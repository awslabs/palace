#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Spatial coupon geometry from a version-2 `SpatialEdgeCluster` signature (the v2 cluster
contract of SURFACE-RESPONSE-IDENTIFICATION.md (a) / (d), decision 91 lane 2).

The identification describes a cluster by its canonical signature alone: the claimed
perimeter portions `P` = [x0, y0, x1, y1] in units of R in the canonical frame (origin at the
length-weighted centroid, z = process normal), each with its in-plane gap direction (toward
the non-metal side), its canonical conductor label, interface types and boundary law, plus
the cluster's vertices (corners / endpoints / junctions). No mesh facets are exported. The
coupon is therefore built as a pure function of the signature, in the canonical frame:

* the edge rows of the spatial generator (`generate_spatial_response.py` Edges: Point,
  GapDirection, ProcessNormal, Interval, Conductor, InterfaceSlot) are the portions scaled
  by R; a portion end that touches no vertex and no other portion is a claim cut - the chain
  continues straight there (a vertex within 2R of the cluster's claimed perimeter would have
  joined the cluster) - and the row is lengthened so that the generator's box rule extends it
  to the coupon box; a portion end at a vertex or at another portion keeps its length (the
  box rule may still extend a long row there: the metal comes from the mask below, so the
  overshoot only pads the box);
* two claim cuts of the same chain FACING each other (collinear free ends of the same gap
  direction, conductor, interfaces and law, less than 2R apart - the piece between them is
  claimed by another feature, e.g. a translational stack piece beyond the cluster's event
  reach between two claims of the cluster; S1p's 41-edge loop end, stage 1) are an interior
  cut: the metal edge continues straight across the piece, so the mask gets a bridging chain
  segment between the two ends and neither end is lengthened; the rows and the model's Edges
  stay the claimed portions (the piece is not claimed by this model);
* the plan-view mask is the set of faces of the planar arrangement of the extended chains
  inside the coupon box (the generator's own box, `coupon_bounds`) that lie on the metal
  side of their bounding portions (the gap direction points away from the metal), one facet
  polygon set per conductor, ear-clipped into triangles; the plan-view boundary follows by
  `prepare_surface_response_coupons.canonical_plan_view_boundary`.

Every geometric decision is checked and fails closed: chains crossing inside the box, a face
whose bounding portions disagree on the metal side or the conductor, a face bounded by the
box alone, a portion end whose interface types match no slot of the record. Chain ends that
coincide within COINCIDENCE_OVER_R (the connection tolerance of end_states: signature
coordinates are rounded to the 1e-6 R grid, so an arc's chord end and the straight portion
it meets can differ by one quantum) are snapped to one point before the arrangement is built.

Spatial-support contract v3 (USER decision 281, decisions 282 / 285 / 286): a signature
carrying ``Box`` + ``Context`` (the identification keyed it with the device plan clipped to the
box: non-legacy context or a grown box) is built in its OWN frame (the canonical frame, M =
identity) inside the signature's Box, and its metal is the arrangement of the claims plus the
Context pieces (the continuation chains, ``Chain`` true, and the foreign edges): no straight
extension, no fictitious metal (the D2 / D3-C defect), no interior-cut bridge (the piece
between two facing claim cuts is a Context piece claimed by the other feature). The coupon then carries
``Geometry.SupportBox`` (the generator's box), ``Geometry.Edges`` = the claims rows plus the
context rows (``Context`` true, ``Chain``: the mesher's owner lookup and attributes for every
conductor of the plan) and ``Geometry.ForeignEdges`` (the foreign pieces as 3D segments: the
config's ``EdgeExcludeSegments``, rule B5). A claims-only signature takes the legacy path
unchanged (byte-identical generator inputs, test_cluster_signature_geometry).
"""
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
for path in (str(HERE), str(HERE.parents[1] / "cpw2d"), str(HERE.parents[1] / "surface_response_identification")):
    if path not in sys.path:
        sys.path.insert(0, path)
import generate_spatial_response as spatial_generator  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402
import signature_library  # noqa: E402

PROCESS_NORMAL = (0.0, 0.0, 1.0)
# Signature coordinates live on the 1e-6 R grid (SignatureLengthQuantumOverR); two portion
# ends or a portion end and a vertex coincide when they agree to within a few quanta (an arc
# end computed from its centre and radius and the rounded straight end it meets differ by up
# to one quantum). The one tolerance of "the same point" in this module: end_states, the
# interior-cut bridges and the arrangement's node snapping all use it.
COINCIDENCE_OVER_R = 1.0e-5
# Block (b) step 0, F0-a (curved-clusters-20261005/DESIGN.md section 0): the serialised Gap
# components and portion ends are rounded to the signature quantum (1e-6 R), so a short
# OBLIQUE straight portion of length L (units of R) reads |tangent . gap| up to
# sqrt(2) / 2 q + 2 q / L from the rounding alone (the axis-aligned gaps of rectilinear
# clusters read exactly 0). The perpendicularity test admits twice that bound; EVERY straight
# row with 0 < |tangent . gap| <= bound has its gap RE-DERIVED as the exact perpendicular of
# its chord with the serialised sign (decision 317 MAJOR-1 option (a): no threshold inside
# the bound), so the generator's and Palace's frames stay exactly orthogonal; a row with
# |tangent . gap| == 0 exactly (every rectilinear row built so far) is kept bitwise - it
# re-derives to itself. An arc CHORD's gap is the arc's radial direction at the chord's
# middle, perpendicular by construction up to round-off: it is not a serialised number and
# is kept as computed.
SIGNATURE_QUANTUM_OVER_R = 1.0e-6
GAP_PERPENDICULARITY_MARGIN = 2.0
# An interior cut (two facing claim cuts of one chain) spans less than the cluster's event
# reach: the piece between two claims of the same cluster is a translational remainder shorter
# than 2R (a longer piece would carry its own events and belong to the cluster).
INTERIOR_CUT_REACH_OVER_R = 2.0
MASK_REGULARIZATION = {"Version": 1, "PhysicalBoundary": "TaperAndRound", "ContinuationBoundary": "Vertical"}


class SignatureGeometryError(ValueError):
    """The signature admits no unambiguous coupon geometry (the reason is named)."""


# Block (b) DESIGN section 2 (b) + A1 (decision 303): the ends of an ARC entry (claim or
# context) are fixed FIRST - the serialised ends, each possibly replaced by a joint snap onto
# the end of a neighbouring piece of the same conductor within ARC_JOINT_SNAP_OVER_R (the
# identification's arc-fit tolerance kArcFitToleranceOverRadius 1e-3 R + 2 quanta, read
# inclusive) or by a face snap onto a box face within one quantum - and the circle is REBUILT
# through them (centre = the point of the perpendicular bisector nearest the serialised
# centre), so every chord vertex is concyclic to double precision, the straight neighbours
# are untouched and the arrangement closes at the joints. The rebuilt circle deviates from
# the signature's by at most the snap (asserted <= ARC_REBUILD_TOLERANCE_OVER_R, inclusive).
# Straight / straight joints keep COINCIDENCE_OVER_R.
ARC_FIT_TOLERANCE_OVER_R = 1.0e-3
ARC_JOINT_SNAP_OVER_R = ARC_FIT_TOLERANCE_OVER_R + 2.0 * SIGNATURE_QUANTUM_OVER_R
ARC_REBUILD_TOLERANCE_OVER_R = ARC_FIT_TOLERANCE_OVER_R + 2.0 * SIGNATURE_QUANTUM_OVER_R
ARC_FACE_SNAP_OVER_R = SIGNATURE_QUANTUM_OVER_R
HALF_QUANTUM_OVER_R = 0.5 * SIGNATURE_QUANTUM_OVER_R


def _arc_sweep(a, b, c, m):
    """The signed sweep from a to b about c through m (signature_library's rule); 2 pi for a
    closed circle (equal ends)."""
    r = math.hypot(a[0] - c[0], a[1] - c[1])
    if math.hypot(a[0] - b[0], a[1] - b[1]) <= 1.0e-9 * max(r, 1.0):
        return 2.0 * math.pi, True
    angle = lambda q: math.atan2(q[1] - c[1], q[0] - c[0])  # noqa: E731
    ta, tb, tm = angle(a), angle(b), angle(m)
    ccw = (tb - ta) % (2.0 * math.pi)
    return (ccw if ((tm - ta) % (2.0 * math.pi)) <= ccw + 1.0e-12 else ccw - 2.0 * math.pi), False


def rebuilt_arc(a, b, c, m, radius, step_degrees=signature_library.CLUSTER_ARC_CHORD_STEP_DEGREES,
                max_chord_over_R=signature_library.CLUSTER_ARC_CHORD_MAX_LENGTH_OVER_R):
    """The chord vertices of an arc whose ENDS a, b (mesh units; the fixed, possibly snapped
    ends) are kept exactly and whose circle is rebuilt through them: centre c' = the point of
    the perpendicular bisector of ab nearest the serialised centre c, radius r' = |a - c'|; the
    side of the arc is the side of the serialised midpoint m; interior vertices at equal
    angular steps (n = max(ceil(sweep / step), ceil(r' sweep / (max_chord R)), 1)). A closed
    circle (equal ends) keeps the serialised circle. Returns (vertices, centre, radius, sweep)
    with vertices[0] is a and vertices[-1] is b exactly."""
    a, b, c, m = (np.asarray(v, dtype=float) for v in (a, b, c, m))
    sweep, closed = _arc_sweep(a, b, c, m)
    if closed:
        centre, r = c, float(np.linalg.norm(a - c))
    else:
        midpoint = 0.5 * (a + b)
        chord = b - a
        bisector = np.asarray([-chord[1], chord[0]]) / float(np.linalg.norm(chord))
        centre = midpoint + float(np.dot(c - midpoint, bisector)) * bisector
        r = float(np.linalg.norm(a - centre))
        sweep, _ = _arc_sweep(a, b, centre, m)
    ta = math.atan2(a[1] - centre[1], a[0] - centre[0])
    n = max(int(math.ceil(abs(sweep) / math.radians(step_degrees) - 1.0e-9)),
            int(math.ceil(r * abs(sweep) / (max_chord_over_R * radius) - 1.0e-9)), 1)
    vertices = [a.copy()]
    for k in range(1, n):
        t = ta + sweep * k / n
        vertices.append(np.asarray([centre[0] + r * math.cos(t), centre[1] + r * math.sin(t)]))
    vertices.append(a.copy() if closed else b.copy())
    return vertices, centre, r, sweep


def _serialised_entries(signature, radius):
    """Every entry of the signature (claims then context) in mesh units: {Kind, Index, A, B,
    Arc (centre, midpoint) or None, Conductor, Entry}."""
    entries = []
    for kind, key in (("Claim", "Portions"), ("Context", "Context")):
        for index, entry in enumerate(signature.get(key, [])):
            p = [float(v) * radius for v in entry["P"]]
            arc = None if "Arc" not in entry else [float(v) * radius for v in entry["Arc"]]
            entries.append({"Kind": kind, "Index": index, "A": np.asarray(p[:2]), "B": np.asarray(p[2:]),
                            "Arc": arc, "Conductor": int(entry["Conductor"]), "Entry": entry})
    return entries


def arc_end_snaps(signature, radius):
    """The fixed ends of every ARC entry of the signature (block (b) DESIGN A1 (1), 2 (b)):
    per (Kind, Index) the pair (a', b') in mesh units, and the JointSnaps records. An arc end
    is snapped onto a box face it lies within one quantum of (inclusive), then onto the end of
    a neighbouring piece of the same conductor within ARC_JOINT_SNAP_OVER_R (inclusive): the
    straight neighbour's end (the device vertex) wins, between two arcs the end of the arc of
    smaller serialised radius; the pieces must be each other's nearest such end (two candidate
    points within the tolerance, or a non-mutual pair, fail closed). Closed circles are not
    snapped."""
    entries = _serialised_entries(signature, radius)
    box = support_box(signature, radius)
    joint_tolerance = (ARC_JOINT_SNAP_OVER_R + HALF_QUANTUM_OVER_R) * radius
    face_tolerance = (ARC_FACE_SNAP_OVER_R + HALF_QUANTUM_OVER_R) * radius
    coincidence = COINCIDENCE_OVER_R * radius

    def radius_of(entry):
        return float(np.linalg.norm(entry["A"] - np.asarray(entry["Arc"][:2]))) if entry["Arc"] else math.inf

    def is_closed(entry):
        return entry["Arc"] is not None and float(np.linalg.norm(entry["A"] - entry["B"])) <= 1.0e-9 * max(radius_of(entry), 1.0)

    def candidates(i, point):
        """Candidate target points for the end `point` of entry i: the ends of other pieces of
        the same conductor within the joint tolerance, grouped by coincident point."""
        groups = []
        for j, other in enumerate(entries):
            if j == i or other["Conductor"] != entries[i]["Conductor"]:
                continue
            for end_index, q in enumerate((other["A"], other["B"])):
                distance = float(np.linalg.norm(point - q))
                if distance > joint_tolerance:
                    continue
                for group in groups:
                    if float(np.linalg.norm(group["Point"] - q)) <= coincidence:
                        group["Members"].append((j, end_index, distance))
                        break
                else:
                    groups.append({"Point": q, "Members": [(j, end_index, distance)]})
        return groups

    fixed = {}
    records = []
    for i, entry in enumerate(entries):
        if entry["Arc"] is None or is_closed(entry):
            continue
        ends = []
        for end_index, point in enumerate((entry["A"], entry["B"])):
            snapped = point.copy()
            if box is not None:
                for axis, faces in ((0, (box[0], box[2])), (1, (box[1], box[3]))):
                    for face_index, face in enumerate(faces):
                        if snapped[axis] != face and abs(snapped[axis] - face) <= face_tolerance:
                            records.append({"Piece": [entry["Kind"], entry["Index"]], "End": end_index, "Class": "Face",
                                            "Face": 2 * face_index + axis, "DistanceOverR": abs(snapped[axis] - face) / radius})
                            snapped[axis] = face
            groups = candidates(i, snapped)
            if len(groups) > 1:
                raise SignatureGeometryError(f"{entry['Kind'].lower()} {entry['Index']} end {end_index} lies within the arc joint "
                                             f"tolerance {ARC_JOINT_SNAP_OVER_R:g} R of {len(groups)} distinct piece ends (ambiguous joint)")
            if groups:
                group = groups[0]
                # The target: a straight end (the device vertex) wins; between arcs the end
                # of the arc of smaller serialised radius (this arc keeps its own end when it
                # is the smaller one).
                straight = [m for m in group["Members"] if entries[m[0]]["Arc"] is None]
                members = straight or group["Members"]
                j, k, distance = min(members, key=lambda item: (radius_of(entries[item[0]]), item[2], item[0], item[1]))
                target = entries[j]["A"] if k == 0 else entries[j]["B"]
                # Mutual: the target end's own candidates are this end's point alone.
                back = candidates(j, target)
                if len(back) != 1 or not any(m[0] == i for m in back[0]["Members"]):
                    raise SignatureGeometryError(f"{entry['Kind'].lower()} {entry['Index']} end {end_index} and "
                                                 f"{entries[j]['Kind'].lower()} {entries[j]['Index']} end {k} are not each "
                                                 f"other's nearest free end (fail closed)")
                if not straight and radius_of(entry) <= radius_of(entries[j]):
                    target = snapped  # the smaller (or equal) arc keeps its serialised end
                if float(np.linalg.norm(target - snapped)) > 0.0:
                    records.append({"Piece": [entry["Kind"], entry["Index"]], "End": end_index,
                                    "Class": "ArcJoint" if straight else "ArcArcJoint",
                                    "To": [entries[j]["Kind"], entries[j]["Index"], k],
                                    "DistanceOverR": float(np.linalg.norm(target - snapped)) / radius})
                    snapped = target.copy()
            ends.append(snapped)
        fixed[(entry["Kind"], entry["Index"])] = (ends[0], ends[1])
    return fixed, records


# Block (b) DESIGN A3 (1): the turn at a joint of two entries (the angle between their
# tangents at the shared end, read on the SIGNATURE, pre-snap) classifies the joint SMOOTH
# (|turn| <= JUNCTION_TANGENT_ANGLE, the curved-offset constant of the mesher: no corner, a
# shared tube section) or a CORNER (the ball, the caps and the clearance rule). One constant,
# one concept; the mesher reads the classification from the plan-view boundary tags.
JUNCTION_TANGENT_ANGLE = 1.0e-4


def _entry_tangent_at(entry, point, radius):
    """The unit tangent (travel direction A -> B) of a serialised entry at its end `point`
    (A or B): the chord direction of a straight entry, the circle tangent of an arc entry."""
    a, b = entry["A"], entry["B"]
    if entry["Arc"] is None:
        d = b - a
        return d / float(np.linalg.norm(d))
    c = np.asarray(entry["Arc"][:2])
    m = np.asarray(entry["Arc"][2:])
    sweep, closed = _arc_sweep(a, b, c, m)
    if closed:
        raise SignatureGeometryError("a closed circle has no joint")
    rad = point - c
    tangent = np.asarray([-rad[1], rad[0]]) / float(np.linalg.norm(rad))
    return tangent if sweep > 0.0 else -tangent


def rebuilt_arcs(signature, radius):
    """Every ARC entry of the signature (claims then context) rebuilt through its fixed ends
    (arc_end_snaps + rebuilt_arc), in mesh units of the canonical frame: a list of
    {ArcId (1-based, in entry order), Kind, Index, Vertices (the chord vertices, ends exact),
    Centre, Radius, Sweep (signed, the travel A -> B), Sign (GapRadial), Conductor,
    Joints [turn at A, turn at B]} with the JointSnaps records (arc_end_snaps' plus one ``Arc``
    record per rebuilt arc {Piece, ArcDeviationOverR, CentreShiftOverR, RadiusShiftOverR,
    Chords, Centre, RadiusOverR}); a rebuilt arc farther than ARC_REBUILD_TOLERANCE_OVER_R
    from the signature's circle fails closed.  The joint turn at an end is the angle between
    the arc's serialised tangent there and the serialised tangent of the neighbouring entry
    of the same conductor whose end lies within ARC_JOINT_SNAP_OVER_R (design A3 (1): read on
    the signature, before any snap); None at an end without a neighbour (a box face, a free
    end).  A closed circle has no joints and keeps its serialised circle."""
    fixed, records = arc_end_snaps(signature, radius)
    entries = _serialised_entries(signature, radius)
    tolerance = (ARC_REBUILD_TOLERANCE_OVER_R + HALF_QUANTUM_OVER_R) * radius
    joint_tolerance = (ARC_JOINT_SNAP_OVER_R + HALF_QUANTUM_OVER_R) * radius
    arcs = []
    for i, entry in enumerate(entries):
        if entry["Arc"] is None:
            continue
        kind, index = entry["Kind"], entry["Index"]
        a, b = entry["A"], entry["B"]
        c, m = np.asarray(entry["Arc"][:2]), np.asarray(entry["Arc"][2:])
        a_fixed, b_fixed = fixed.get((kind, index), (a, b))
        vertices, centre, r, sweep = rebuilt_arc(a_fixed, b_fixed, c, m, radius)
        r_signature = float(np.linalg.norm(a - c))
        deviation = max(abs(float(np.linalg.norm(q - c)) - r_signature) for q in vertices)
        if deviation > tolerance:
            raise SignatureGeometryError(f"{kind.lower()} {index}: arc rebuild outside the fit tolerance (the rebuilt chord "
                                         f"vertices deviate {deviation / radius:.3e} R from the signature's circle, tolerance "
                                         f"{ARC_REBUILD_TOLERANCE_OVER_R:g} R: a snap bent a short arc too far)")
        records.append({"Piece": [kind, index], "Class": "Arc", "ArcDeviationOverR": deviation / radius,
                        "CentreShiftOverR": float(np.linalg.norm(centre - c)) / radius,
                        "RadiusShiftOverR": abs(r - r_signature) / radius, "Chords": len(vertices) - 1,
                        "Centre": [float(centre[0]) / radius, float(centre[1]) / radius], "RadiusOverR": r / radius})
        closed = float(np.linalg.norm(a - b)) <= 1.0e-9 * max(r_signature, 1.0)
        joints = [None, None]
        if not closed:
            for end_index, point in enumerate((a, b)):
                own = _entry_tangent_at(entry, point, radius)
                for j, other in enumerate(entries):
                    if j == i or other["Conductor"] != entry["Conductor"]:
                        continue
                    for other_point in (other["A"], other["B"]):
                        if float(np.linalg.norm(point - other_point)) > joint_tolerance:
                            continue
                        if other["Arc"] is not None and float(np.linalg.norm(other["A"] - other["B"])) <= 1.0e-9 * radius:
                            continue
                        neighbour = _entry_tangent_at(other, other_point, radius)
                        # The two travel directions meet head-on or tail-to-head at the joint:
                        # the turn is the angle between the lines of travel through it.
                        cosine = abs(float(np.dot(own, neighbour)))
                        turn = math.acos(min(1.0, cosine))
                        joints[end_index] = turn if joints[end_index] is None else min(joints[end_index], turn)
        arcs.append({"ArcId": len(arcs) + 1, "Kind": kind, "Index": index, "Vertices": vertices, "Centre": centre,
                     "Radius": r, "Sweep": sweep, "Sign": int(entry["Entry"]["GapRadial"]), "Conductor": entry["Conductor"],
                     "Joints": joints, "Closed": closed})
    return arcs, records


def chorded_entries(signature, radius, include_context=False):
    """The builder's plan-view edges: signature_library.cluster_plan_view_edges with every
    arc entry rebuilt through its fixed ends (arc_end_snaps + rebuilt_arc), so the chord ends
    are the serialised / snapped ends exactly and every chord vertex is concyclic. Returns
    (edges, records): the edges in cluster_plan_view_edges' encoding (``Chord`` / ``Chords``
    on arc chords, ``Context`` / ``Chain`` on context entries), the JointSnaps records with
    one ``Arc`` record per rebuilt arc {Piece, ArcDeviationOverR, CentreShiftOverR,
    RadiusShiftOverR, Chords}; a rebuilt arc farther than ARC_REBUILD_TOLERANCE_OVER_R from
    the signature's circle fails closed."""
    arcs, records = rebuilt_arcs(signature, radius)
    rebuilt = {(arc["Kind"], arc["Index"]): arc for arc in arcs}
    edges = []
    entries = [(index, portion, False) for index, portion in enumerate(signature["Portions"])]
    if include_context:
        entries += [(index, portion, True) for index, portion in enumerate(signature.get("Context", []))]
    for index, portion, context in entries:
        common = {"Conductor": int(portion["Conductor"]), "Interfaces": portion.get("Interfaces", []), "Law": portion.get("Law"), "Portion": index}
        if context:
            common.update({"Context": True, "Chain": bool(portion.get("Chain", False))})
        p = [float(v) * radius for v in portion["P"]]
        if "Arc" not in portion:
            edges.append(dict(common, P0=(p[0], p[1]), P1=(p[2], p[3]), Gap=(float(portion["Gap"][0]), float(portion["Gap"][1]))))
            continue
        arc = rebuilt[("Context" if context else "Claim", index)]
        vertices, centre, sweep, sign = arc["Vertices"], arc["Centre"], arc["Sweep"], arc["Sign"]
        n = len(vertices) - 1
        ta = math.atan2(vertices[0][1] - centre[1], vertices[0][0] - centre[0])
        for k in range(n):
            tmid = ta + sweep * (k + 0.5) / n
            edges.append(dict(common, P0=(float(vertices[k][0]), float(vertices[k][1])), P1=(float(vertices[k + 1][0]), float(vertices[k + 1][1])),
                              Gap=(sign * math.cos(tmid), sign * math.sin(tmid)), Chord=k, Chords=n))
    return edges, records


def context_from_signature(signature, radius):
    """The Context pieces of a contract-v3 SpatialEdgeCluster signature (decision 282 rule B1 /
    B3: the device plan clipped to the Box, minus the claims) in mesh units (canonical frame),
    chorded like the portions, each with its ``Chain`` flag (True = the coupon's own
    continuation chain, False = a FOREIGN edge: present in the geometry, excluded from the
    within-R accounting by ``EdgeExcludeSegments``, rule B5). Empty for a claims-only
    (contract-2) signature."""
    pieces = []
    if "Context" not in signature:
        return pieces
    edges, _ = chorded_entries(signature, radius, include_context=True)
    for edge in edges:
        if not edge.get("Context"):
            continue
        p0, p1 = np.asarray(edge["P0"], dtype=float), np.asarray(edge["P1"], dtype=float)
        label = f"context piece {edge['Portion']}"
        gap, rederived, _ = perpendicular_gap(p0, p1, np.asarray(edge["Gap"], dtype=float), radius, label,
                                             serialised=edge.get("Chord") is None)
        length = float(np.linalg.norm(p1 - p0))
        pieces.append({"P0": p0, "P1": p1, "Gap": gap, "Length": length, "Conductor": int(edge["Conductor"]),
                       "Interfaces": sorted(edge.get("Interfaces") or []), "Law": edge.get("Law") or '{"Type":"PEC"}',
                       "Portion": edge["Portion"], "Chain": bool(edge.get("Chain", False)), "GapRederived": rederived})
    return pieces


def support_box(signature, radius):
    """The signature's support box [x0, y0, x1, y1] in mesh units (canonical frame): the
    serialised ``Box`` of a contract-v3 signature (rule B2 + the face rules T1-T3, grown by the
    identification); None for a claims-only signature (the generator's own ``coupon_bounds``
    then applies, as before)."""
    if "Box" not in signature:
        return None
    box = [float(v) * radius for v in signature["Box"]]
    if len(box) != 4 or box[2] <= box[0] or box[3] <= box[1]:
        raise SignatureGeometryError(f"the signature's Box {signature['Box']} is not a box")
    return box


def gap_perpendicularity_bound(length_over_R):
    """The quantisation bound on |tangent . gap| of a serialised straight row of length
    ``length_over_R`` (units of R) with the F0-a margin: GAP_PERPENDICULARITY_MARGIN x
    (sqrt(2) / 2 q + 2 q / L), q = SIGNATURE_QUANTUM_OVER_R."""
    q = SIGNATURE_QUANTUM_OVER_R
    return GAP_PERPENDICULARITY_MARGIN * (math.sqrt(2.0) * 0.5 * q + 2.0 * q / length_over_R)


def perpendicular_gap(p0, p1, gap, radius, label, serialised=True):
    """The unit gap of a straight row (F0-a, decision 317 MAJOR-1 option (a)): the serialised
    gap bitwise when |tangent . gap| == 0 exactly, else - for every deviation within the
    quantisation bound - the exact perpendicular of the chord with the serialised sign, else
    a fail-closed refusal naming the bound. An arc chord (``serialised`` False: its gap is the
    arc's radial direction at the chord's middle, not a serialised number) is kept as computed
    within the same bound. Returns (unit gap, re-derived flag, deviation)."""
    norm = float(np.linalg.norm(gap))
    if norm <= 0.0:
        raise SignatureGeometryError(f"{label} has an invalid P or Gap")
    length = float(np.linalg.norm(p1 - p0))
    if length <= 0.0:
        raise SignatureGeometryError(f"{label} has zero length")
    tangent = (p1 - p0) / length
    gap = gap / norm
    deviation = abs(float(np.dot(tangent, gap)))
    bound = gap_perpendicularity_bound(length / radius)
    if deviation > bound:
        raise SignatureGeometryError(f"{label}: Gap is not perpendicular to the portion (|tangent . gap| = "
                                     f"{deviation:.3e} > the quantisation bound {bound:.3e} for a {length / radius:.6f} R row)")
    if deviation == 0.0 or not serialised:
        return gap, False, deviation
    perpendicular = np.asarray([tangent[1], -tangent[0]])
    sign = 1.0 if float(np.dot(perpendicular, gap)) >= 0.0 else -1.0
    return sign * perpendicular + 0.0, True, deviation  # + 0.0: no negative zero in the rows


def portions_from_signature(signature, radius):
    """The portions of a SpatialEdgeCluster signature in mesh units (canonical frame). An arc
    portion (option A: ``Arc`` = centre + midpoint, ``GapRadial``) is chorded at the canonical
    step (signature_library.cluster_plan_view_edges: 5 deg / 0.25 R), each chord a straight
    portion whose gap direction is the arc's radial direction at the chord's middle; the
    chords are what the coupon's plan view and the model's Edges carry (Palace places a model
    carrying its Signature with the identity map and verifies the chords against the arc).
    A straight portion's gap is tested against the quantisation bound and re-derived when
    the serialised rounding tilts it (perpendicular_gap, F0-a); ``GapRederived`` marks it."""
    if signature.get("Type") != "SpatialEdgeCluster":
        raise SignatureGeometryError(f"not a SpatialEdgeCluster signature: {signature.get('Type')!r}")
    portions = []
    for index, entry in enumerate(signature["Portions"]):
        if "Arc" not in entry and ("Gap" not in entry or len(entry["P"]) != 4):
            raise SignatureGeometryError(f"portion {index} has an invalid P or Gap")
    edges, _ = chorded_entries(signature, radius)
    for edge in edges:
        index = edge["Portion"]
        p0, p1 = np.asarray(edge["P0"], dtype=float), np.asarray(edge["P1"], dtype=float)
        gap, rederived, _ = perpendicular_gap(p0, p1, np.asarray(edge["Gap"], dtype=float), radius, f"portion {index}",
                                             serialised=edge.get("Chord") is None)
        length = float(np.linalg.norm(p1 - p0))
        portions.append({"P0": p0, "P1": p1, "Gap": gap, "Length": length, "Conductor": int(edge["Conductor"]),
                         "Interfaces": sorted(edge.get("Interfaces") or []), "Law": edge.get("Law") or '{"Type":"PEC"}',
                         "Portion": index, "GapRederived": rederived})
    return portions


def vertex_points(signature, radius):
    return [np.asarray([float(v) * radius for v in vertex["P"]]) for vertex in signature.get("Vertices", [])]


def slot_of(interfaces, record_interfaces):
    """The interface slot whose type set equals the portion's interface types."""
    by_slot = {}
    for entry in record_interfaces:
        by_slot.setdefault(int(entry.get("Slot", 0)), set()).add(entry["Type"])
    matches = [slot for slot, types in by_slot.items() if sorted(types) == sorted(interfaces)]
    if len(matches) != 1:
        raise SignatureGeometryError(f"portion interfaces {interfaces} match {len(matches)} slots of {sorted(by_slot)}")
    return matches[0]


def end_states(portions, vertices, radius):
    """Per portion: (begin free, end free): an end touching a vertex or another portion's end
    is connected; every other end is a claim cut (free)."""
    tolerance = COINCIDENCE_OVER_R * radius
    ends = [(p["P0"], p["P1"]) for p in portions]
    states = []
    for i, (a, b) in enumerate(ends):
        free = []
        for point in (a, b):
            connected = any(np.linalg.norm(point - v) <= tolerance for v in vertices)
            connected = connected or any(
                np.linalg.norm(point - other) <= tolerance for j, pair in enumerate(ends) if j != i for other in pair)
            free.append(not connected)
        states.append(tuple(free))
    return states


def interior_bridges(portions, states, radius):
    """The interior cuts of the chains: pairs of free ends (portion index, end index: 0 = P0,
    1 = P1) that face each other on one line - the second end lies on the ray leaving the
    first portion at its free end, within COINCIDENCE_OVER_R x R of the line, less than
    INTERIOR_CUT_REACH_OVER_R x R away, with the same gap direction, conductor, interfaces and
    law, and its own ray leading back to the first. Returns (bridges, states): the bridges as
    (i, end_i, j, end_j) with i < j and the states with the bridged ends connected. Fails
    closed when a free end faces more than one candidate (no unambiguous continuation)."""
    tolerance = COINCIDENCE_OVER_R * radius
    reach = INTERIOR_CUT_REACH_OVER_R * radius
    free_ends = [(i, end) for i, state in enumerate(states) for end, free in enumerate(state) if free]

    def ray(i, end):
        portion = portions[i]
        point = portion["P1"] if end else portion["P0"]
        direction = (portion["P1"] - portion["P0"]) / portion["Length"]
        return point, (direction if end else -direction)

    def same_chain(a, b):
        return (a["Conductor"] == b["Conductor"] and sorted(a["Interfaces"]) == sorted(b["Interfaces"])
                and a["Law"] == b["Law"] and float(np.linalg.norm(a["Gap"] - b["Gap"])) <= 1.0e-6)

    def facing(i, end_i, j, end_j):
        if i == j or not same_chain(portions[i], portions[j]):
            return False
        point_i, direction_i = ray(i, end_i)
        point_j, direction_j = ray(j, end_j)
        offset = point_j - point_i
        along = float(np.dot(offset, direction_i))
        if not tolerance < along < reach or float(np.linalg.norm(offset - along * direction_i)) > tolerance:
            return False
        return float(np.dot(direction_i, direction_j)) < -1.0 + 1.0e-9   # anti-parallel: j's ray leads back

    partner = {}
    for i, end_i in free_ends:
        candidates = [(j, end_j) for j, end_j in free_ends if facing(i, end_i, j, end_j)]
        if len(candidates) > 1:
            raise SignatureGeometryError(f"portion {i} end {end_i} faces {len(candidates)} collinear free ends")
        if candidates:
            partner[(i, end_i)] = candidates[0]
    bridges = []
    states = [list(state) for state in states]
    for (i, end_i), (j, end_j) in sorted(partner.items()):
        if partner.get((j, end_j)) != (i, end_i):
            raise SignatureGeometryError(f"portion {i} end {end_i} faces portion {j} end {end_j} but not the reverse")
        if i < j:
            bridges.append((i, end_i, j, end_j))
            states[i][end_i] = False
            states[j][end_j] = False
    return bridges, [tuple(state) for state in states]


def bridge_segments(portions, bridges):
    """The mask's chain segments across the interior cuts (the metal edge continues straight
    between the two facing ends; the gap direction and conductor of the chain)."""
    segments = []
    for i, end_i, j, end_j in bridges:
        p0 = portions[i]["P1"] if end_i else portions[i]["P0"]
        p1 = portions[j]["P1"] if end_j else portions[j]["P0"]
        segments.append({"P0": p0.copy(), "P1": p1.copy(), "Gap": portions[i]["Gap"], "Conductor": portions[i]["Conductor"]})
    return segments


def edge_rows(portions, states, radius, record_interfaces, boundary_condition):
    """The generator's Edges (canonical frame): Point on the portion, Interval along the
    generator's tangent (gap x normal), free ends lengthened to at least R from Point so the
    box rule (`extended_interval`) carries them 2R further, to the coupon box."""
    rows = []
    for portion, (begin_free, end_free) in zip(portions, states):
        gap = portion["Gap"]
        tangent = np.asarray([gap[1], -gap[0]])   # np.cross(gap, normal) for normal +z
        p0, p1, length = portion["P0"], portion["P1"], portion["Length"]
        midpoint = 0.5 * (p0 + p1)
        # The end at +tangent from the midpoint is the interval end, the other the begin.
        forward_is_p1 = float(np.dot(p1 - midpoint, tangent)) > 0.0
        end_is_free = end_free if forward_is_p1 else begin_free
        begin_is_free = begin_free if forward_is_p1 else end_free
        half = 0.5 * length
        begin = -max(half, radius) if begin_is_free else -half
        end = max(half, radius) if end_is_free else half
        law = portion["Law"]
        rows.append({"Point": [float(midpoint[0]), float(midpoint[1]), 0.0],
                     "GapDirection": [float(gap[0]), float(gap[1]), 0.0],
                     "ProcessNormal": list(PROCESS_NORMAL), "Interval": [float(begin), float(end)],
                     "Conductor": portion["Conductor"], "InterfaceSlot": slot_of(portion["Interfaces"], record_interfaces),
                     "BoundaryCondition": json.loads(law) if isinstance(law, str) else dict(law or boundary_condition)})
    return rows


def context_slot(piece, pieces_and_portions, record_interfaces):
    """The interface slot of a context piece: its own interface types' slot; a piece without
    interface types (a perimeter segment excluded before the identification — a port cut has
    no target interface) takes the slot of the nearest claim or context piece of its conductor
    that has one (its edge surfaces are built like its neighbours')."""
    if piece["Interfaces"]:
        return slot_of(piece["Interfaces"], record_interfaces)
    best = None
    for other in pieces_and_portions:
        if other is piece or not other["Interfaces"] or other["Conductor"] != piece["Conductor"]:
            continue
        distance = min(np.linalg.norm(a - b) for a in (piece["P0"], piece["P1"]) for b in (other["P0"], other["P1"]))
        if best is None or distance < best[0]:
            best = (distance, other)
    if best is None:
        raise SignatureGeometryError("a context piece without interface types has no neighbour of its conductor "
                                     "to take its interface slot from")
    return slot_of(best[1]["Interfaces"], record_interfaces)


def context_rows(pieces, record_interfaces, boundary_condition, portions=()):
    """The generator's Edges rows of the context pieces (canonical frame): Point at the piece's
    midpoint, the exact Interval [-L/2, L/2] along the generator's tangent (no lengthening: a
    context piece already reaches the face it crosses, the box is the signature's), flagged
    ``Context`` (the trace basis ignores them: knot columns come from the mask vertices on the
    faces, interior cap hats stay within R of the CLAIMS) with the ``Chain`` class."""
    rows = []
    for piece in pieces:
        gap = piece["Gap"]
        midpoint = 0.5 * (piece["P0"] + piece["P1"])
        half = 0.5 * piece["Length"]
        law = piece["Law"]
        rows.append({"Point": [float(midpoint[0]), float(midpoint[1]), 0.0],
                     "GapDirection": [float(gap[0]), float(gap[1]), 0.0],
                     "ProcessNormal": list(PROCESS_NORMAL), "Interval": [-half, half],
                     "Conductor": piece["Conductor"],
                     "InterfaceSlot": context_slot(piece, list(portions) + list(pieces), record_interfaces),
                     "BoundaryCondition": json.loads(law) if isinstance(law, str) else dict(law or boundary_condition),
                     "Context": True, "Chain": bool(piece["Chain"])})
    return rows


def foreign_edge_segments(pieces):
    """The FOREIGN context pieces (Chain false) as 3D segments [x0, y0, z, x1, y1, z] in mesh
    units of the canonical frame (z = 0, the process plane): the coupon config's
    ``EdgeExcludeSegments`` (rule B5), so that the within-R accounting covers the coupon's own
    edges only."""
    return [[float(p["P0"][0]), float(p["P0"][1]), 0.0, float(p["P1"][0]), float(p["P1"][1]), 0.0]
            for p in pieces if not p["Chain"]]


def exact_portion_edges(portions, rows):
    """The model's Edges: the claimed portions exactly (Point = P0, Interval [0, L] along the
    generator's tangent, or [-L, 0] when the tangent points from P1 to P0) - what
    ModelClusterSignature canonicalises to the feature's own signature."""
    edges = []
    for portion, row in zip(portions, rows):
        gap = portion["Gap"]
        tangent = np.asarray([gap[1], -gap[0]])
        forward = float(np.dot(portion["P1"] - portion["P0"], tangent)) > 0.0
        interval = [0.0, portion["Length"]] if forward else [-portion["Length"], 0.0]
        edges.append({**row, "Point": [float(portion["P0"][0]), float(portion["P0"][1]), 0.0], "Interval": interval})
    return edges


def _quantize(point, quantum):
    return (int(math.floor(point[0] / quantum + 0.5)), int(math.floor(point[1] / quantum + 0.5)))


def _segment_intersection(a0, a1, b0, b1, tolerance):
    """The proper intersection point of two segments (interior of both), or None."""
    d1, d2 = a1 - a0, b1 - b0
    denominator = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denominator) <= 1.0e-14 * (np.linalg.norm(d1) * np.linalg.norm(d2)):
        return None
    r = b0 - a0
    t = (r[0] * d2[1] - r[1] * d2[0]) / denominator
    u = (r[0] * d1[1] - r[1] * d1[0]) / denominator
    la, lb = np.linalg.norm(d1), np.linalg.norm(d2)
    if tolerance / la < t < 1.0 - tolerance / la and tolerance / lb < u < 1.0 - tolerance / lb:
        return a0 + t * d1
    return None


def _clip_to_box(point, direction, box):
    """Farthest parameter s >= 0 with point + s direction inside the box (ray exit)."""
    (x0, y0), (x1, y1) = box
    s = math.inf
    for k, (lo, hi) in enumerate(((x0, x1), (y0, y1))):
        if abs(direction[k]) > 1.0e-15:
            candidates = [(lo - point[k]) / direction[k], (hi - point[k]) / direction[k]]
            s = min(s, max(candidates))
    if not math.isfinite(s) or s < 0.0:
        raise SignatureGeometryError("a chain end cannot be extended to the coupon box")
    return s


def extended_chain_segments(portions, states, box, radius):
    """The portions with their free ends extended to the box boundary (in the frame of `box`);
    every portion must lie inside the box."""
    (x0, y0), (x1, y1) = box
    segments = []
    for portion, (begin_free, end_free) in zip(portions, states):
        p0, p1 = portion["P0"].copy(), portion["P1"].copy()
        for point in (p0, p1):
            if not (x0 - 1.0e-9 * radius <= point[0] <= x1 + 1.0e-9 * radius and
                    y0 - 1.0e-9 * radius <= point[1] <= y1 + 1.0e-9 * radius):
                raise SignatureGeometryError("a claimed portion lies outside the coupon box")
        direction = (p1 - p0) / portion["Length"]
        if begin_free:
            p0 = p0 - _clip_to_box(p0, -direction, box) * direction
        if end_free:
            p1 = p1 + _clip_to_box(p1, direction, box) * direction
        segments.append({"P0": p0, "P1": p1, "Gap": portion["Gap"], "Conductor": portion["Conductor"]})
    return segments


def snap_chain_ends(segments, radius):
    """The segments with every end that coincides with an earlier end within
    COINCIDENCE_OVER_R x R moved onto it (the connection tolerance of end_states; the
    arrangement below resolves nodes far more finely and would otherwise see a one-quantum
    rounding difference as a gap in the chain)."""
    tolerance = COINCIDENCE_OVER_R * radius
    representatives = []

    def snapped(point):
        for representative in representatives:
            if float(np.linalg.norm(point - representative)) <= tolerance:
                return representative.copy()
        representatives.append(np.asarray(point, dtype=float).copy())
        return representatives[-1].copy()

    return [{**s, "P0": snapped(s["P0"]), "P1": snapped(s["P1"])} for s in segments]


def plan_view_faces(segments, box, radius):
    """Faces of the planar arrangement of the chain segments and the box boundary: a list of
    (polygon points ccw, metal flag, conductor). Fails closed on crossing chains, on a face
    whose bounding chain segments disagree, and on a face bounded by the box alone."""
    quantum = 1.0e-9 * radius
    tolerance = 1.0e-7 * radius
    (x0, y0), (x1, y1) = box
    segments = snap_chain_ends(segments, radius)
    # Proper crossings between chain segments are a geometry the signature cannot describe.
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            if _segment_intersection(segments[i]["P0"], segments[i]["P1"], segments[j]["P0"], segments[j]["P1"],
                                     tolerance) is not None:
                raise SignatureGeometryError(f"chain segments {i} and {j} cross inside the coupon box")
    # Nodes: quantized points; the box corners; chain ends; T-junctions split the segment they touch.
    nodes = {}

    def node(point):
        key = _quantize(point, quantum)
        if key not in nodes:
            nodes[key] = np.asarray([key[0] * quantum, key[1] * quantum])
        return key

    raw = [(np.asarray((x0, y0)), np.asarray((x1, y0)), None), (np.asarray((x1, y0)), np.asarray((x1, y1)), None),
           (np.asarray((x1, y1)), np.asarray((x0, y1)), None), (np.asarray((x0, y1)), np.asarray((x0, y0)), None)]
    raw += [(s["P0"], s["P1"], index) for index, s in enumerate(segments)]
    points = [p for a, b, _ in raw for p in (a, b)]
    edges = []
    for a, b, owner in raw:
        direction = b - a
        length = float(np.linalg.norm(direction))
        if length <= tolerance:
            raise SignatureGeometryError("a zero-length chain segment")
        split = [0.0, 1.0]
        for p in points:
            offset = p - a
            t = float(np.dot(offset, direction)) / length ** 2
            if tolerance / length < t < 1.0 - tolerance / length:
                distance = abs(direction[0] * offset[1] - direction[1] * offset[0]) / length
                if distance <= tolerance:
                    split.append(t)
        split = sorted(set(split))
        for t0, t1 in zip(split, split[1:]):
            key0, key1 = node(a + t0 * direction), node(a + t1 * direction)
            if key0 != key1:
                edges.append((key0, key1, owner))
    # Half-edge structure: outgoing half-edges per node sorted by angle.
    outgoing = {}
    half_edges = []
    for key0, key1, owner in edges:
        for start, stop in ((key0, key1), (key1, key0)):
            half_edges.append((start, stop, owner))
            outgoing.setdefault(start, []).append(len(half_edges) - 1)

    def angle(index):
        start, stop, _ = half_edges[index]
        d = nodes[stop] - nodes[start]
        return math.atan2(d[1], d[0])

    for start in outgoing:
        outgoing[start].sort(key=angle)
    twin = {}
    for index, (start, stop, owner) in enumerate(half_edges):
        twin[index] = index ^ 1
    # Face tracing: from a half-edge, the next is the outgoing edge at its head that comes
    # right after the twin in clockwise order (the face on the left of every half-edge).
    visited = [False] * len(half_edges)
    faces = []
    for seed in range(len(half_edges)):
        if visited[seed]:
            continue
        cycle = []
        current = seed
        while not visited[current]:
            visited[current] = True
            cycle.append(current)
            start, stop, _ = half_edges[current]
            candidates = outgoing[stop]
            position = candidates.index(twin[current])
            current = candidates[(position - 1) % len(candidates)]
        if current != seed:
            raise SignatureGeometryError("the plan-view arrangement is not a closed cell complex")
        polygon = [nodes[half_edges[index][0]] for index in cycle]
        area = 0.5 * sum(p[0] * q[1] - p[1] * q[0] for p, q in zip(polygon, polygon[1:] + polygon[:1]))
        if area <= 0.0:
            continue   # the outer face (clockwise)
        metal, conductors = set(), set()
        for index in cycle:
            start, stop, owner = half_edges[index]
            if owner is None:
                continue
            gap = segments[owner]["Gap"]
            d = nodes[stop] - nodes[start]
            left = np.asarray((-d[1], d[0]))
            metal.add(float(np.dot(left, gap)) < 0.0)
            conductors.add(segments[owner]["Conductor"])
        if not metal:
            raise SignatureGeometryError("a face of the coupon box is bounded by no chain segment")
        if len(metal) != 1:
            raise SignatureGeometryError("the chain segments bounding a face disagree on its metal side")
        is_metal = metal.pop()
        if is_metal and len(conductors) != 1:
            raise SignatureGeometryError("a metal face is bounded by portions of different conductors")
        faces.append((polygon, is_metal, conductors.pop() if is_metal else None))
    if not faces:
        raise SignatureGeometryError("the plan-view arrangement has no faces")
    return faces


def triangulate(polygon):
    """Ear clipping of a simple counter-clockwise polygon (list of 2D points)."""
    points = [np.asarray(p, dtype=float) for p in polygon]
    # Drop repeated / collinear vertices first.
    changed = True
    while changed and len(points) > 3:
        changed = False
        for i in range(len(points)):
            a, b, c = points[i - 1], points[i], points[(i + 1) % len(points)]
            cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if abs(cross) <= 1.0e-18 or np.linalg.norm(b - a) <= 0.0:
                points.pop(i)
                changed = True
                break
    triangles = []
    indices = list(range(len(points)))

    def inside(p, a, b, c):
        def side(u, v, w):
            return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])
        return side(a, b, p) >= -1.0e-15 and side(b, c, p) >= -1.0e-15 and side(c, a, p) >= -1.0e-15

    guard = 0
    while len(indices) > 3:
        guard += 1
        if guard > 10 * len(polygon) ** 2:
            raise SignatureGeometryError("ear clipping did not converge (non-simple polygon)")
        clipped = False
        for k in range(len(indices)):
            i0, i1, i2 = indices[k - 1], indices[k], indices[(k + 1) % len(indices)]
            a, b, c = points[i0], points[i1], points[i2]
            cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
            if cross <= 1.0e-18:
                continue   # reflex or degenerate vertex
            if any(inside(points[j], a, b, c) for j in indices if j not in (i0, i1, i2)):
                continue
            triangles.append([a.tolist(), b.tolist(), c.tolist()])
            indices.pop(k)
            clipped = True
            break
        if not clipped:
            raise SignatureGeometryError("ear clipping found no ear (non-simple polygon)")
    triangles.append([points[indices[0]].tolist(), points[indices[1]].tolist(), points[indices[2]].tolist()])
    return triangles


def model_edges(record, radius):
    """The exact-portion Edges of a version-2 record's model (see exact_portion_edges)."""
    signature = record["Signature"] if "Signature" in record else record["Geometry"]["Signature"]
    portions = portions_from_signature(signature, radius)
    states = end_states(portions, vertex_points(signature, radius), radius)
    rows = edge_rows(portions, states, radius, record["Interfaces"], record.get("BoundaryCondition", {"Type": "PEC"}))
    return exact_portion_edges(portions, rows)


def cluster_coupon(record, radius, metal_thickness, overetch):
    """The spatial generator's coupon (Topology, Geometry {Edges, EdgeCount, PlanViewFacets,
    PlanViewBoundary, MaskRegularization, Signature}, Interfaces, BoundaryCondition) of a
    version-2 SpatialEdgeCluster requirement record, plus the model's exact-portion Edges:
    (coupon, model_edges)."""
    signature = record["Signature"] if "Signature" in record else record["Geometry"]["Signature"]
    portions = portions_from_signature(signature, radius)
    vertices = vertex_points(signature, radius)
    states = end_states(portions, vertices, radius)
    box_from_signature = support_box(signature, radius)
    bridges = []
    if box_from_signature is None:
        # Interior cuts (two facing claim cuts of one chain, stage 1's loop end) are bridged in
        # the legacy mask only: under contract v3 the piece between the cuts is a Context piece
        # of the device plan (claimed by the other feature, Chain false) and is drawn from the
        # signature itself - a straight bridge would duplicate it.
        bridges, states = interior_bridges(portions, states, radius)
    boundary_condition = record.get("BoundaryCondition", {"Type": "PEC"})
    rows = edge_rows(portions, states, radius, record["Interfaces"], boundary_condition)
    exact_rows = rows
    coupon = {"Topology": "SpatialEdgeCluster",
              "Geometry": {"EdgeCount": len(rows), "Edges": rows, "Signature": signature},
              "Interfaces": record["Interfaces"], "BoundaryCondition": boundary_condition}
    if box_from_signature is None:
        # The legacy contract (contract 2, decision 236; every context piece a straight
        # continuation of a claim, no growth): the generator's own frame and box (a rotation
        # about the process normal of the canonical frame), the free ends extended straight
        # to the box; the mask is built inside that box and returned in the canonical frame.
        # Byte-identical to the pre-v3 builder.
        frame, local_edges, _ = spatial_generator.normalize_geometry(coupon, radius)
        lower, upper = spatial_generator.coupon_bounds(local_edges, radius, metal_thickness, overetch)
        box = ((float(lower[0]), float(lower[1])), (float(upper[0]), float(upper[1])))
        rotation = np.asarray(frame)[:2, :2]

        def to_local(point):
            return rotation @ np.asarray(point)

        local_portions = [{**p, "P0": to_local(p["P0"]), "P1": to_local(p["P1"]), "Gap": to_local(p["Gap"])} for p in portions]
        segments = extended_chain_segments(local_portions, states, box, radius) + bridge_segments(local_portions, bridges)
        conductors = {p["Conductor"] for p in portions}
    else:
        # Contract 3 (decision 282 rules B1-B3, the R1a ruling MAJOR-2): the box frame IS the
        # canonical frame (M = identity), the box the signature's grown Box, the metal the
        # device plan: the claims plus the Context pieces (own continuation chains and foreign
        # edges), cut at the faces by the identification; nothing is continued straight past a
        # device vertex and no fictitious metal boundary exists.
        frame = np.identity(3)
        pieces = context_from_signature(signature, radius)
        box = ((box_from_signature[0], box_from_signature[1]), (box_from_signature[2], box_from_signature[3]))
        (x0, y0), (x1, y1) = box
        for item in portions + pieces:
            for point in (item["P0"], item["P1"]):
                if not (x0 - 1.0e-6 * radius <= point[0] <= x1 + 1.0e-6 * radius and
                        y0 - 1.0e-6 * radius <= point[1] <= y1 + 1.0e-6 * radius):
                    raise SignatureGeometryError("a claimed portion or context piece lies outside the signature's Box")
        segments = [{"P0": p["P0"], "P1": p["P1"], "Gap": p["Gap"], "Conductor": p["Conductor"]} for p in portions + pieces]
        # The claims rows exactly (no lengthening: the continuation is in the context) and
        # the context rows.
        rows = exact_portion_edges(portions, rows) + context_rows(pieces, record["Interfaces"], boundary_condition, portions)
        coupon["Geometry"]["Edges"] = rows
        coupon["Geometry"]["SupportBox"] = list(box_from_signature)
        coupon["Geometry"]["ForeignEdges"] = foreign_edge_segments(pieces)
        coupon["Geometry"]["ContextEdgeCount"] = len(pieces)
        conductors = {p["Conductor"] for p in portions} | {p["Conductor"] for p in pieces}
    faces = plan_view_faces(segments, box, radius)
    facets = []
    for polygon, is_metal, conductor in faces:
        if not is_metal:
            continue
        for triangle in triangulate(polygon):
            points = [(np.asarray(frame).T @ np.asarray([x, y, 0.0])).tolist() for x, y in triangle]
            facets.append({"Conductor": conductor, "Points": points})
    if {f["Conductor"] for f in facets} != conductors:
        raise SignatureGeometryError("the plan-view mask does not cover every conductor of the signature")
    coupon["Geometry"]["PlanViewFacets"] = facets
    coupon["Geometry"]["PlanViewBoundary"] = planner.canonical_plan_view_boundary(facets, radius, 2)
    coupon["Geometry"]["MaskRegularization"] = dict(MASK_REGULARIZATION)
    _, joint_snaps = chorded_entries(signature, radius, include_context=True)
    if joint_snaps:
        # Block (b) DESIGN A1 (3) / 2 (b): the arc end snaps (face / joint) and the rebuilt
        # circles' deviations from the signature's (recorded; a legacy straight coupon has
        # none, so its coupon.json is unchanged).
        coupon["Geometry"]["JointSnaps"] = joint_snaps
    if bridges:
        coupon["Geometry"]["InteriorCuts"] = [
            {"Portions": [portions[i]["Portion"], portions[j]["Portion"]],
             "Ends": [[float(v) for v in (portions[i]["P1"] if end_i else portions[i]["P0"])],
                      [float(v) for v in (portions[j]["P1"] if end_j else portions[j]["P0"])]],
             "LengthOverR": float(np.linalg.norm((portions[j]["P1"] if end_j else portions[j]["P0"])
                                                 - (portions[i]["P1"] if end_i else portions[i]["P0"]))) / radius}
            for i, end_i, j, end_j in bridges]
    return coupon, exact_portion_edges(portions, exact_rows)
