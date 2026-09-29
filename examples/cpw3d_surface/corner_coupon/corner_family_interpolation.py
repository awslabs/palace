#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Angle interpolation of the corner family (USER decision 121 (C); corner-qualification block
2026-09-29): the Python mirror of the runtime's CornerBasisEvents / SelectCornerFamilyStencil
(palace/models/cornertracebasis.cpp, used by MatchCornerFamily).

Geometric EVENTS of the trace basis rule: a knot of a ring that meets the metal passes a vertex
of the fixed layout — a box corner or a side midpoint. At every event the perimeter-ordered
band triangulation flips a quad diagonal and the hats JUMP (measured O(1); the corner slave does
not help: the knot swaps order with it). Coupons built with a segment connectivity (TraceBasis
ConnectivityAngleDegrees: the band merge keyed by the rule's layout at that angle) have one
triangulation over the whole segment, which removes the midpoint flips; a knot-corner passage
remains a jump between two segments, so the family carries one coupon per side there.

STENCIL rule: coupons sharing a connectivity angle form a segment whose node angles must lie in
one corner-event-free interval with the connectivity angle (fail closed); a device angle
strictly inside a segment's node range is interpolated by Lagrange on the segment's nodes
nearest to it (cubic on four, else quadratic / linear: never across a corner event); an angle
equal to a node (within the tolerance) is exact — with several coupons at that angle the legacy
one (the tie triangulation at its own angle, e.g. the recorded 90 / 135 coupons) is preferred,
else the coupon of the lower-angle segment; legacy coupons (no connectivity record) are never
interpolated; outside the node range refused (no extrapolation). The node tolerance is closed on
the exact side and open on the interior side with the same floating-point difference
(|node - angle| <= tol exact, > tol interior): an angle at the boundary is exact or interior,
never in no segment (qualification review 2026-09-29 m1). The segment structure (a segment
across a knot-corner passage, two coupons at one angle, overlapping segments) is checked by
check_segments = the runtime's CheckCornerFamilySegments, which fails closed at library load."""
import math

import generate_corner_response as generator

SIGNATURE_ANGLE_TOLERANCE_DEGREES = 1.0e-2  # kSignatureAngleToleranceDegrees
FIRST_CROSSING_FRACTION = 0.5  # the +x arm meets the box at (R, 0)


def role_coefficients(topology, ring_size=8):
    """(a, b) per knot role: the role's unwrapped perimeter fraction is a s2 + b with s2 the
    second crossing's unwrapped fraction in (0.5, 1.5] (metal_ring_layout)."""
    if ring_size != 2 + generator.METAL_INTERIOR_KNOTS + generator.FREE_KNOTS:
        raise ValueError("ring size inconsistent with the trace basis rule")
    first = FIRST_CROSSING_FRACTION
    roles = {"crossing2": (1.0, 0.0)}
    for m in range(1, generator.METAL_INTERIOR_KNOTS + 1):
        t = m / (generator.METAL_INTERIOR_KNOTS + 1)
        roles[f"metal{m}"] = (t, first * (1.0 - t)) if topology == "convex" else (
            1.0 - t, (first + 1.0) * t)
    for k in range(1, generator.FREE_KNOTS + 1):
        t = k / (generator.FREE_KNOTS + 1)
        roles[f"free{k}"] = (1.0 - t, (first + 1.0) * t) if topology == "convex" else (
            t, first * (1.0 - t))
    return roles


def basis_events(topology, ring_size=8):
    """Every event in (0, 180] degrees: dicts angle_degrees / role / fraction (the fixed-layout
    vertex passed) / corner (a box corner, else a side midpoint), sorted by angle."""
    first = FIRST_CROSSING_FRACTION
    events = []
    for role, (a, b) in role_coefficients(topology, ring_size).items():
        for j in range(ring_size):
            c = j / ring_size
            for n in range(3):
                s2 = (c + n - b) / a
                if s2 <= first + generator.KNOT_COINCIDENCE_FRACTION or s2 > first + 0.5 + 1e-12:
                    continue
                x, y, _ = generator.square_perimeter_point(1.0, 0.0, s2 % 1.0)
                angle = math.degrees(math.atan2(y, x))
                events.append({
                    "angle_degrees": angle + 360.0 if angle <= 0.0 else angle,
                    "role": role,
                    "fraction": c,
                    "corner": abs((c * 4.0) % 1.0 - 0.5) <= 1e-12,
                })
    return sorted(events, key=lambda event: (event["angle_degrees"], event["role"]))


def corner_event_angles(topology, ring_size=8, tolerance=SIGNATURE_ANGLE_TOLERANCE_DEGREES):
    """The distinct knot-corner passage angles (the segment boundaries)."""
    angles = []
    for event in basis_events(topology, ring_size):
        if event["corner"] and (not angles or event["angle_degrees"] - angles[-1] > tolerance):
            angles.append(event["angle_degrees"])
    return angles


def lagrange_weights(abscissae, x):
    return [
        math.prod((x - abscissae[j]) / (abscissae[i] - abscissae[j])
                  for j in range(len(abscissae)) if j != i)
        for i in range(len(abscissae))
    ]


def _group_segments(nodes, topology, ring_size, tolerance):
    """(boundaries, {connectivity angle: members sorted by angle}, legacy present)."""
    boundaries = corner_event_angles(topology, ring_size, tolerance)
    segments = {}
    legacy = False
    for node in nodes:
        if node[1] is None:
            legacy = True
            continue
        for key in segments:
            if abs(key - node[1]) <= tolerance:
                segments[key].append(node)
                break
        else:
            segments[node[1]] = [node]
    for members in segments.values():
        members.sort(key=lambda node: node[0])
    return boundaries, segments, legacy


def check_segments(nodes, topology, ring_size=8, tolerance=SIGNATURE_ANGLE_TOLERANCE_DEGREES):
    """The runtime's CheckCornerFamilySegments (fail closed at library load, ReadProcessLibrary):
    the reason string, empty when the segment structure is consistent — every segment's
    connectivity angle off the knot-corner passages and its nodes in one event-free interval with
    it, one coupon per angle within a segment, segments overlapping at shared node angles only."""
    boundaries, segments, _ = _group_segments(nodes, topology, ring_size, tolerance)

    def interval(value):
        return sum(1 for boundary in boundaries if boundary < value - tolerance)

    def on_boundary(value):
        return any(abs(value - boundary) <= tolerance for boundary in boundaries)

    ranges = []
    for key, members in segments.items():
        if on_boundary(key):
            return f"segment connectivity angle {key:g} deg lies on a knot-corner passage"
        for member in members:
            if not (abs(member[0] - key) <= tolerance or interval(member[0]) == interval(key)
                    or (on_boundary(member[0]) and interval(member[0]) + 1 == interval(key))):
                return (f"segment with connectivity angle {key:g} deg has a node at {member[0]:g} "
                        "deg across a knot-corner passage of the trace basis: split the segment")
        for previous, following in zip(members, members[1:]):
            if not following[0] - previous[0] > tolerance:
                return f"segment with connectivity angle {key:g} deg has two coupons at {following[0]:g} deg"
        ranges.append((members[0][0], members[-1][0]))
    ranges.sort()
    for previous, following in zip(ranges, ranges[1:]):
        if not following[0] >= previous[1] - tolerance:
            return (f"segments [{previous[0]:g}, {previous[1]:g}] and [{following[0]:g}, "
                    f"{following[1]:g}] deg overlap beyond a shared node angle")
    return ""


def select_stencil(nodes, angle, topology, ring_size=8,
                   tolerance=SIGNATURE_ANGLE_TOLERANCE_DEGREES):
    """The runtime's stencil for `nodes` = [(angle_degrees, connectivity_angle_degrees or
    None, index), ...]. Returns a dict with `nodes` = [(index, weight), ...], `rule`
    (exact / linear / quadratic / cubic), `connectivity_angle_degrees`, `base` (the nearest
    node), or `reason` when refused. Raises ValueError where the runtime fails closed (the
    segment structure of check_segments, verified at library load)."""
    if not nodes:
        raise ValueError("a corner family stencil needs nodes")

    def beyond(first, second):
        return abs(first - second) > tolerance

    exact = None
    for node in nodes:
        if beyond(node[0], angle):
            continue
        if (exact is None or (node[1] is None and exact[1] is not None)
                or (node[1] is not None and exact[1] is not None and node[1] < exact[1])):
            exact = node
    if exact is not None:
        return {"nodes": [(exact[2], 1.0)], "rule": "exact",
                "connectivity_angle_degrees": exact[1], "base": exact[2]}
    low, high = min(n[0] for n in nodes), max(n[0] for n in nodes)
    convexity = topology
    if angle < low:
        return {"reason": f"corner angle {angle:g} deg is sharper than the smallest {convexity} "
                          f"coupon angle {low:g} deg (no extrapolation)"}
    if angle > high:
        return {"reason": f"corner angle {angle:g} deg is wider than the widest {convexity} "
                          f"coupon angle {high:g} deg (no extrapolation; the first-order regime "
                          "needs the straight anchor)"}
    reason = check_segments(nodes, topology, ring_size, tolerance)
    if reason:
        raise ValueError(reason)
    _, segments, legacy = _group_segments(nodes, topology, ring_size, tolerance)
    segment = None
    for key, members in segments.items():
        if (members[0][0] < angle < members[-1][0] and beyond(members[0][0], angle)
                and beyond(members[-1][0], angle)):
            segment, segment_key = members, key
    if segment is None:
        if legacy:
            return {"reason": f"corner angle {angle:g} deg needs interpolation but the corner family "
                              "is built without segment connectivity records (TraceBasis "
                              "ConnectivityAngleDegrees): its trace bases jump at every knot "
                              "passage of a box vertex; rebuild the family on segment connectivity "
                              "(corner-qualification block 2026-09-29)"}
        return {"reason": f"corner angle {angle:g} deg lies in no segment of the corner family "
                          "(segments end at the knot-corner passages of the trace basis: a node "
                          "is needed on each side of every passage)"}
    interval_index = 0
    while interval_index + 1 < len(segment) and segment[interval_index + 1][0] < angle:
        interval_index += 1
    begin, end = 0, len(segment)
    if len(segment) > 4:
        begin = min(interval_index - 1 if interval_index > 0 else 0, len(segment) - 4)
        end = begin + 4
    window = segment[begin:end]
    weights = lagrange_weights([node[0] for node in window], angle)
    base = min(window, key=lambda node: abs(node[0] - angle))[2]
    return {"nodes": [(node[2], weight) for node, weight in zip(window, weights)],
            "rule": {4: "cubic", 3: "quadratic", 2: "linear"}[len(window)],
            "connectivity_angle_degrees": segment_key, "base": base}


def segment_reference_angle(low, high):
    """The recorded choice of a segment's connectivity angle: its midpoint."""
    return 0.5 * (low + high)


if __name__ == "__main__":
    for topology in ("convex", "concave"):
        print(topology)
        for event in basis_events(topology):
            print(f"  {event['angle_degrees']:14.9f}  {event['role']:10s} fraction {event['fraction']:.3f} "
                  f"{'CORNER' if event['corner'] else 'midpoint'}")
