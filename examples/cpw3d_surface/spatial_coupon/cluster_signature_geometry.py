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
* the plan-view mask is the set of faces of the planar arrangement of the extended chains
  inside the coupon box (the generator's own box, `coupon_bounds`) that lie on the metal
  side of their bounding portions (the gap direction points away from the metal), one facet
  polygon set per conductor, ear-clipped into triangles; the plan-view boundary follows by
  `prepare_surface_response_coupons.canonical_plan_view_boundary`.

Every geometric decision is checked and fails closed: chains crossing inside the box, a face
whose bounding portions disagree on the metal side or the conductor, a face bounded by the
box alone, a portion end whose interface types match no slot of the record.
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
# ends or a portion end and a vertex coincide when they agree well within that quantum.
COINCIDENCE_OVER_R = 1.0e-5
MASK_REGULARIZATION = {"Version": 1, "PhysicalBoundary": "TaperAndRound", "ContinuationBoundary": "Vertical"}


class SignatureGeometryError(ValueError):
    """The signature admits no unambiguous coupon geometry (the reason is named)."""


def portions_from_signature(signature, radius):
    """The portions of a SpatialEdgeCluster signature in mesh units (canonical frame). An arc
    portion (option A: ``Arc`` = centre + midpoint, ``GapRadial``) is chorded at the canonical
    step (signature_library.cluster_plan_view_edges: 5 deg / 0.25 R), each chord a straight
    portion whose gap direction is the arc's radial direction at the chord's middle; the
    chords are what the coupon's plan view and the model's Edges carry (Palace places a model
    carrying its Signature with the identity map and verifies the chords against the arc)."""
    if signature.get("Type") != "SpatialEdgeCluster":
        raise SignatureGeometryError(f"not a SpatialEdgeCluster signature: {signature.get('Type')!r}")
    portions = []
    for index, entry in enumerate(signature["Portions"]):
        if "Arc" not in entry and ("Gap" not in entry or len(entry["P"]) != 4):
            raise SignatureGeometryError(f"portion {index} has an invalid P or Gap")
    for edge in signature_library.cluster_plan_view_edges(signature, radius):
        index = edge["Portion"]
        p0, p1 = np.asarray(edge["P0"], dtype=float), np.asarray(edge["P1"], dtype=float)
        gap = np.asarray(edge["Gap"], dtype=float)
        norm = np.linalg.norm(gap)
        if norm <= 0.0:
            raise SignatureGeometryError(f"portion {index} has an invalid P or Gap")
        length = float(np.linalg.norm(p1 - p0))
        if length <= 0.0:
            raise SignatureGeometryError(f"portion {index} has zero length")
        tangent = (p1 - p0) / length
        gap = gap / norm
        if abs(float(np.dot(tangent, gap))) > 1.0e-6:
            raise SignatureGeometryError(f"portion {index}: Gap is not perpendicular to the portion")
        portions.append({"P0": p0, "P1": p1, "Gap": gap, "Length": length, "Conductor": int(edge["Conductor"]),
                         "Interfaces": sorted(edge.get("Interfaces") or []), "Law": edge.get("Law") or '{"Type":"PEC"}',
                         "Portion": index})
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


def plan_view_faces(segments, box, radius):
    """Faces of the planar arrangement of the chain segments and the box boundary: a list of
    (polygon points ccw, metal flag, conductor). Fails closed on crossing chains, on a face
    whose bounding chain segments disagree, and on a face bounded by the box alone."""
    quantum = 1.0e-9 * radius
    tolerance = 1.0e-7 * radius
    (x0, y0), (x1, y1) = box
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
    boundary_condition = record.get("BoundaryCondition", {"Type": "PEC"})
    rows = edge_rows(portions, states, radius, record["Interfaces"], boundary_condition)
    coupon = {"Topology": "SpatialEdgeCluster",
              "Geometry": {"EdgeCount": len(rows), "Edges": rows, "Signature": signature},
              "Interfaces": record["Interfaces"], "BoundaryCondition": boundary_condition}
    # The generator's own frame and box (a rotation about the process normal of the canonical
    # frame): the mask is built inside that box and returned in the canonical frame.
    frame, local_edges, _ = spatial_generator.normalize_geometry(coupon, radius)
    lower, upper = spatial_generator.coupon_bounds(local_edges, radius, metal_thickness, overetch)
    box = ((float(lower[0]), float(lower[1])), (float(upper[0]), float(upper[1])))
    rotation = np.asarray(frame)[:2, :2]

    def to_local(point):
        return rotation @ np.asarray(point)

    local_portions = [{**p, "P0": to_local(p["P0"]), "P1": to_local(p["P1"]), "Gap": to_local(p["Gap"])} for p in portions]
    segments = extended_chain_segments(local_portions, states, box, radius)
    faces = plan_view_faces(segments, box, radius)
    facets = []
    for polygon, is_metal, conductor in faces:
        if not is_metal:
            continue
        for triangle in triangulate(polygon):
            points = [(np.asarray(frame).T @ np.asarray([x, y, 0.0])).tolist() for x, y in triangle]
            facets.append({"Conductor": conductor, "Points": points})
    conductors = {p["Conductor"] for p in portions}
    if {f["Conductor"] for f in facets} != conductors:
        raise SignatureGeometryError("the plan-view mask does not cover every conductor of the signature")
    coupon["Geometry"]["PlanViewFacets"] = facets
    coupon["Geometry"]["PlanViewBoundary"] = planner.canonical_plan_view_boundary(facets, radius, 2)
    coupon["Geometry"]["MaskRegularization"] = dict(MASK_REGULARIZATION)
    return coupon, exact_portion_edges(portions, rows)
