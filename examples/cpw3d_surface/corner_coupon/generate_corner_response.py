#!/usr/bin/env python3

"""Generate a spatial trace basis and Palace configs for a 3D corner coupon."""

import argparse
import json
from pathlib import Path

import numpy as np


INTERFACES = {
    "SA": (4.0, 2.0e-3),
    "MS": (11.47, 3.0e-4),
    "MA": (10.0, 3.0e-2),
}

# Boundary attributes of mesh_corner_coupon.jl (thin and fabricated): the matching box and
# the substrate-air surface. The metal edge lines of every interface are the perimeter of
# the SA surface (the metal outline of the thin sheet; the foot of the fabricated slab where
# the SA plane / trench wall meets the MS face) minus its edges on the matching box. This is
# the line the legacy automatic extraction returned; the version-2 perimeter classification
# (identification phase 3) classifies every edge of the fabricated slab as a fold between
# non-coplanar metal faces (MS / sidewall, top / sidewall) and retains no one-sided metal
# perimeter, so `AutomaticEdges` finds nothing on a fabricated corner coupon.
MATCHING_SURFACE_ATTRIBUTE = 1
SA_ATTRIBUTE = 3


def square_perimeter_point(half_width, z, fraction):
    """Point of the square |x|, |y| <= half_width at the perimeter fraction from (-half_width,
    0) counterclockwise (down the left side first), the traversal of square_ring."""
    coordinate = (fraction % 1.0) * 8.0 * half_width
    if coordinate < half_width:
        return (-half_width, -coordinate, z)
    if coordinate < 3.0 * half_width:
        return (-half_width + coordinate - half_width, -half_width, z)
    if coordinate < 5.0 * half_width:
        return (half_width, -half_width + coordinate - 3.0 * half_width, z)
    if coordinate < 7.0 * half_width:
        return (half_width - (coordinate - 5.0 * half_width), half_width, z)
    return (-half_width, half_width - (coordinate - 7.0 * half_width), z)


def square_perimeter_fraction(half_width, point):
    """Inverse of square_perimeter_point for a point on the square's perimeter."""
    x, y = float(point[0]), float(point[1])
    tolerance = 1.0e-12 * half_width
    if abs(x + half_width) <= tolerance and y <= 0.0:
        coordinate = -y
    elif abs(y + half_width) <= tolerance:
        coordinate = half_width + (x + half_width)
    elif abs(x - half_width) <= tolerance:
        coordinate = 3.0 * half_width + (y + half_width)
    elif abs(y - half_width) <= tolerance:
        coordinate = 5.0 * half_width + (half_width - x)
    elif abs(x + half_width) <= tolerance:
        coordinate = 7.0 * half_width + (half_width - y)
    else:
        raise ValueError("point does not lie on the square perimeter")
    return (coordinate / (8.0 * half_width)) % 1.0


def square_ring(half_width, z, size):
    if size < 8 or size % 8:
        raise ValueError("ring size must be a multiple of eight and at least eight")
    perimeter_coordinate = np.arange(size, dtype=float) * (8.0 * half_width / size)
    # Start at (-half_width, 0) so the Maxwell anchor path lies in the gap
    # instead of crossing the metal quadrant.
    perimeter_coordinate = (perimeter_coordinate + 7.0 * half_width) % (
        8.0 * half_width
    )
    points = []
    for coordinate in perimeter_coordinate:
        if coordinate < 2.0 * half_width:
            x = -half_width + coordinate
            y = -half_width
        elif coordinate < 4.0 * half_width:
            x = half_width
            y = -half_width + coordinate - 2.0 * half_width
        elif coordinate < 6.0 * half_width:
            x = half_width - coordinate + 4.0 * half_width
            y = half_width
        else:
            x = -half_width
            y = half_width - coordinate + 6.0 * half_width
        points.append((x, y, z))
    points = np.asarray(points)
    start = np.argmin(np.linalg.norm(points - np.array([-half_width, 0.0, z]), axis=1))
    return np.roll(points, -start, axis=0)


# Trace basis rule of the corner family (supervisor decision on the corner-family review,
# 2026-09-29). Every box ring that meets the metal (z = 0: the thin sheet and the slab foot;
# z = MetalThickness: the slab top) carries KNOTS whose semantics do not depend on the corner
# angle: the two crossings of the metal arms with the ring (PEC, zero set), METAL_INTERIOR_KNOTS
# knots at equal perimeter-arc-length fractions of the metal arc between them (PEC, zero set)
# and FREE_KNOTS knots at equal fractions of the free arc. Every node of the family therefore
# has the same knot count, the same zero set and like-to-like free knots whose positions vary
# smoothly with the angle, so the entrywise (Lagrange) blend of the nodes' matrices is well
# posed; with METAL_INTERIOR_KNOTS = 1 and FREE_KNOTS = 5 the 90-degree convex node reproduces
# the lane-2 8-knot layout byte-identically (0 / 45 / 90 deg PEC; 135 ... 315 deg free). The
# box corners that are not knots stay vertices of the trace triangulation (the trace surface
# lies on the box): SLAVE vertices whose hat values are the linear interpolation in the
# perimeter fraction between the two neighbouring knots (trace-vertices.csv columns parent_a /
# parent_b / weight_a; basis 0). Rings that do not meet the metal keep the fixed layout of
# square_ring (knots at the fractions k / RingSize; every corner a knot). No free hat has
# support on the PEC part of the box contour (checked by free_hat_pec_support). The straight
# anchor (180 degrees) follows the same rule (its second crossing is (-R, 0)).
METAL_INTERIOR_KNOTS = 1
FREE_KNOTS = 5
FRACTION_PARAMETRISATION = "PerimeterArcLength"
# Knot coincidence in perimeter fraction (corner-basis-fix review m3, 2026-09-29; the same
# value as kKnotCoincidenceFraction in palace/models/cornertracebasis.cpp): a free or
# metal-interior knot within this fraction of a fixed-layout fraction k / RingSize (a box
# corner or a side midpoint) takes that fraction exactly (snap_fraction), and a box corner
# within it of any knot gets no slave vertex. The value is set by the surface mortar's
# degenerate-triangle threshold (surfaceresponseoperator.cpp: area <= 1e-14 max(1, L^2)
# in mesh units): the nearest a slave corner can come to a knot is 8e-6 R along the
# perimeter (the perimeter is 8R), so the smallest triangle of a constructed basis has an
# edge of 8e-6 R and an area >= 4e-6 R h with h the smallest ring spacing (R = 1.9 um,
# h >= 0.05 um: 3.8e-7 um^2, six orders above the threshold and an aspect ratio <= 1e4),
# while a snap moves a knot by at most 8e-6 R = 15 pm at R = 1.9 um, below the coupon mesh
# resolution. The former value 1e-9 was a floating-point identity tolerance only and left
# slivers of edge 8e-9 .. 8e-6 R (not degenerate, but arbitrarily thin). Crossing knots are
# never moved (they must lie on the arm: the gate's PEC knot at the crossing); a crossing
# within the fraction of a corner suppresses that corner's slave instead, so the trace
# surface cuts the corner by at most 8e-6 R.
KNOT_COINCIDENCE_FRACTION = 1.0e-6
# Floating-point identity of a fraction with k / RingSize for the reuse of square_ring's
# coordinates in ring_points (byte identity of the recorded fixed layout).
FIXED_FRACTION_IDENTITY = 1.0e-12


def trace_basis_rule(ring_size, connectivity_angle_degrees=None):
    """The TraceBasis record. With a connectivity angle the coupon is a node of an
    interpolation SEGMENT of the corner family (corner-qualification block 2026-09-29): the
    bands next to its metal rings are triangulated in the merge order of the rule's layout at
    that angle (connectivity_keys), the same for every coupon of the segment, so the family's
    hats do not jump between the segment's nodes (they do at every knot passage of a
    fixed-layout vertex under the perimeter-ordered merge; corner_family_interpolation.py
    lists the events). Without one the coupon is a legacy node: exact matches only."""
    record = {
        "Rule": (
            "rings at z = 0 and z = MetalThickness: knots = the two metal-arm crossings "
            "(PEC), MetalInteriorKnots at equal fractions of the metal arc (PEC) and "
            "FreeKnots at equal fractions of the free arc (free), fractions = perimeter arc "
            "length; box corners that are no knot are slave trace vertices (linear in the "
            "fraction between the neighbouring knots); other rings: knots at the fractions "
            "k / RingSize"
        ),
        "RingSize": int(ring_size),
        "MetalInteriorKnots": METAL_INTERIOR_KNOTS,
        "FreeKnots": FREE_KNOTS,
        "Fractions": FRACTION_PARAMETRISATION,
        "MetalLevels": ["0", "MetalThickness"],
    }
    if connectivity_angle_degrees is not None:
        record["ConnectivityAngleDegrees"] = float(connectivity_angle_degrees)
        record["Connectivity"] = (
            "the bands next to the metal rings are merged in the perimeter order of the "
            "rule's layout at ConnectivityAngleDegrees (knots by role, box corners by their "
            "fraction): one triangulation for every coupon of the interpolation segment"
        )
    return record


def arm_crossing_fractions(radius, angle_degrees):
    """Perimeter fractions where the two metal arms of a sharp or rounded corner (the arm
    lines through the apex: the first along +x, the second at the corner angle
    counterclockwise; a rounded corner's tangency points lie inside the box, so the arms
    are straight where they cross it) meet the matching box |x|, |y| <= radius. The
    straight anchor's second arm is the -x axis. Returned in counterclockwise order from
    the first arm: (first, second) with first < second <= first + 1."""
    angle, _, _, _ = corner_frame(angle_degrees)
    fractions = []
    for direction in (np.array([1.0, 0.0]), np.array([np.cos(angle), np.sin(angle)])):
        scale = radius / np.max(np.abs(direction))
        point = scale * direction
        point = np.where(
            np.abs(np.abs(point) - radius) <= 1.0e-12 * radius,
            np.sign(point) * radius,
            point,
        )
        fractions.append(square_perimeter_fraction(radius, point))
    first, second = fractions
    if second <= first + KNOT_COINCIDENCE_FRACTION:
        second += 1.0
    return first, second


def metal_arc_fractions(radius, angle_degrees, topology):
    """The metal and the free arc of a ring that meets the metal as unwrapped perimeter
    fraction intervals (start, end) with start < end <= start + 1, counterclockwise. The
    metal arc of a convex corner runs from the first arm crossing to the second (the sector
    0 .. angle); the concave corner's metal is the complement."""
    first, second = arm_crossing_fractions(radius, angle_degrees)
    if topology == "convex":
        return (first, second), (second, first + 1.0)
    return (second, first + 1.0), (first, second)


def heldout_ring_layout(radius, angle_degrees, topology, ring_size):
    """Vertices of a ring that meets the metal on the fine held-out surface (USER decision
    149 (6), option (c)): the fixed layout k / ring_size plus the two metal-arm crossings
    (a crossing within KNOT_COINCIDENCE_FRACTION of a fixed vertex is that vertex), so the
    piecewise-linear held-out trace, zero at every vertex of the metal arc, vanishes exactly
    on the PEC part of the ring and nowhere else. Returns (fraction, kind, slot) triples in
    perimeter order with kind "zero" on the closed metal arc and "free" elsewhere; slots
    are sequential; no slave vertices (every box corner is a fixed vertex)."""
    metal, _ = metal_arc_fractions(radius, angle_degrees, topology)
    fractions = [k / ring_size for k in range(ring_size)]
    for crossing in metal:
        crossing %= 1.0
        if all(
            min(abs(crossing - fraction), 1.0 - abs(crossing - fraction))
            > KNOT_COINCIDENCE_FRACTION
            for fraction in fractions
        ):
            fractions.append(crossing)
    fractions.sort()
    layout = []
    for slot, fraction in enumerate(fractions):
        kind = "zero" if perimeter_arc_distance(fraction, metal) == 0.0 else "free"
        layout.append((fraction, kind, slot))
    return layout


def perimeter_arc_distance(fraction, arc):
    """Perimeter-fraction distance from `fraction` (any real) to the closed arc
    (start, end), start < end <= start + 1 unwrapped; zero on the arc (a fraction within
    FIXED_FRACTION_IDENTITY of an end point is on the arc, so a vertex placed at a crossing
    counts as PEC whatever the rounding of its fraction)."""
    start, end = arc
    shifted = start + (fraction - start) % 1.0
    if (
        shifted <= end + FIXED_FRACTION_IDENTITY
        or shifted >= start + 1.0 - FIXED_FRACTION_IDENTITY
    ):
        return 0.0
    return min(shifted - end, start + 1.0 - shifted)


def metal_ring_layout(radius, angle_degrees, topology, ring_size):
    """Knots and slave vertices of a ring that meets the metal, by the trace basis rule.
    Returns a list of vertices in perimeter (counterclockwise) order, each (fraction, kind,
    slot) with kind "zero" (PEC knot), "free" (free knot) or "slave" (box corner between two
    knots, slot None); the fraction is in [0, 1). The slot is the knot's basis position
    within the ring, fixed by ROLE so that every node of the family has the same zero set
    and like-to-like free knots at the same indices: convex rings are ordered free 2, free
    3, free 4, free 5, first crossing, metal interior, second crossing, free 1 (free k = the
    k-th free knot counterclockwise from the second crossing), which is the lane-2 order of
    the 90-degree node (its free 2 lies at (-R, 0), the ring start); concave rings are
    ordered first crossing, free 1 .. free 5 (counterclockwise from the first crossing),
    second crossing, metal interior."""
    if ring_size != 2 + METAL_INTERIOR_KNOTS + FREE_KNOTS:
        raise ValueError(
            "the trace basis rule needs RingSize = 2 crossings + MetalInteriorKnots + "
            f"FreeKnots = {2 + METAL_INTERIOR_KNOTS + FREE_KNOTS}"
        )
    first, second = arm_crossing_fractions(radius, angle_degrees)
    metal, free = metal_arc_fractions(radius, angle_degrees, topology)
    roles = {"crossing1": first, "crossing2": second % 1.0}
    for m in range(1, METAL_INTERIOR_KNOTS + 1):
        roles[f"metal{m}"] = snap_fraction(
            (metal[0] + (metal[1] - metal[0]) * m / (METAL_INTERIOR_KNOTS + 1)) % 1.0,
            ring_size,
        )
    for k in range(1, FREE_KNOTS + 1):
        roles[f"free{k}"] = snap_fraction(
            (free[0] + (free[1] - free[0]) * k / (FREE_KNOTS + 1)) % 1.0, ring_size
        )
    order = ring_role_order(topology)
    knots = []
    for slot, role in enumerate(order):
        kind = "free" if role.startswith("free") else "zero"
        knots.append((roles[role], kind, slot))
    knots.sort(key=lambda knot: knot[0])
    for previous, following in zip(knots, knots[1:]):
        if following[0] - previous[0] <= KNOT_COINCIDENCE_FRACTION:
            raise ValueError("two knots of the trace basis rule coincide")
    vertices = list(knots)
    for corner in (0.125, 0.375, 0.625, 0.875):
        separation = min(
            min(abs(corner - fraction), 1.0 - abs(corner - fraction))
            for fraction, _, _ in knots
        )
        if separation > KNOT_COINCIDENCE_FRACTION:
            vertices.append((corner, "slave", None))
    vertices.sort(key=lambda vertex: vertex[0])
    return vertices


def ring_role_order(topology):
    """Basis order of the knot roles within a ring that meets the metal (see
    metal_ring_layout)."""
    free = [f"free{k}" for k in range(1, FREE_KNOTS + 1)]
    metal = [f"metal{m}" for m in range(1, METAL_INTERIOR_KNOTS + 1)]
    if topology == "convex":
        return free[1:] + ["crossing1"] + metal + ["crossing2"] + free[:1]
    return ["crossing1"] + free + ["crossing2"] + metal


def zero_slots(topology):
    return [
        slot for slot, role in enumerate(ring_role_order(topology))
        if not role.startswith("free")
    ]


def fixed_ring_layout(ring_size):
    return [(k / ring_size, "free", k) for k in range(ring_size)]


def snap_fraction(fraction, ring_size):
    """A fraction within KNOT_COINCIDENCE_FRACTION of a fixed knot k / ring_size (a box
    corner or a side midpoint) is that fraction exactly (the C++ CornerMetalRingLayout
    applies the same snap; see KNOT_COINCIDENCE_FRACTION)."""
    nearest = round(fraction * ring_size)
    if abs(fraction - nearest / ring_size) <= KNOT_COINCIDENCE_FRACTION:
        return (nearest % ring_size) / ring_size
    return fraction


def ring_points(half_width, z, layout, ring_size):
    """Points of a ring layout: square_perimeter_point at every fraction, except that a
    fraction equal to k / ring_size (within FIXED_FRACTION_IDENTITY) reuses square_ring's
    coordinates so the recorded fixed layout stays byte-identical (the two arithmetics
    leave different rounding residues of ~1e-15 R). The C++ runtime (SquarePerimeterPoint in
    cornertracebasis.cpp, the same formula as square_perimeter_point) places every knot by
    the fraction alone: the two agree to floating-point rounding (~1e-15 R) at the fixed
    fractions and to FIXED_FRACTION_IDENTITY x 8 R elsewhere; MatchCornerFamily checks the
    recorded files against the C++ rule at 1e-9 R (test_ring_points_agree_with_
    square_perimeter_point)."""
    fixed = square_ring(half_width, z, ring_size)
    points = []
    for fraction, _, _ in layout:
        index = round(fraction * ring_size)
        if abs(fraction - index / ring_size) <= FIXED_FRACTION_IDENTITY:
            points.append(tuple(fixed[index % ring_size]))
        else:
            points.append(square_perimeter_point(half_width, z, fraction))
    return np.asarray(points)


def connect_rings_by_fraction(
    triangles, first, first_fractions, second, second_fractions
):
    """Triangulate the band between two closed rings that share the perimeter parameter
    (both start at fraction 0 and run counterclockwise) but may carry different vertices.
    Equal vertex sets reproduce the lane-2 connect_rings exactly (same triangles, same
    order); a vertex of one ring strictly between two vertices of the other adds one
    triangle to the quad. `first` / `second` are the offsets of the rings' vertices."""
    first_count, second_count = len(first_fractions), len(second_fractions)
    tolerance = 1.0e-12
    i = j = 0
    while i < first_count or j < second_count:
        next_first = (
            first_fractions[i + 1] if i + 1 < first_count else 1.0
        ) if i < first_count else np.inf
        next_second = (
            second_fractions[j + 1] if j + 1 < second_count else 1.0
        ) if j < second_count else np.inf
        first_vertex = first + (i % first_count)
        second_vertex = second + (j % second_count)
        if next_first < next_second - tolerance:
            triangles.append((first_vertex, first + (i + 1) % first_count, second_vertex))
            i += 1
        elif next_second < next_first - tolerance:
            triangles.append((first_vertex, second + (j + 1) % second_count, second_vertex))
            j += 1
        else:
            next_first_vertex = first + (i + 1) % first_count
            next_second_vertex = second + (j + 1) % second_count
            triangles.append((first_vertex, next_first_vertex, next_second_vertex))
            triangles.append((first_vertex, next_second_vertex, second_vertex))
            i += 1
            j += 1


def cap_ring(triangles, offset, size, reverse):
    for index in range(1, size - 1):
        triangle = (offset, offset + index, offset + index + 1)
        triangles.append(tuple(reversed(triangle)) if reverse else triangle)


def metal_levels(metal_thickness):
    """The box rings that meet the metal: the sheet / slab foot at z = 0 and the slab top."""
    return (0.0, metal_thickness)


class TraceSurface:
    """The matching-box trace triangulation: vertices (knots = basis functions, in the
    ring order; slave vertices appended after every knot), triangles, per-ring knot counts
    (ContourGroups), the zero set (0-based knot indices) and the slave table."""

    def __init__(self):
        self.knot_points = []
        self.knot_zero = []
        self.contour_groups = []
        self.slaves = []  # (point, parent_a, parent_b, weight_a) with 0-based parents
        self.triangles = []

    @property
    def basis_size(self):
        return len(self.knot_points)

    def vertex_points(self):
        return np.vstack(
            [np.asarray(self.knot_points)]
            + ([np.asarray([slave[0] for slave in self.slaves])] if self.slaves else [])
        )

    def hat_values(self, basis):
        """Values of the hat of knot `basis` (0-based) at every vertex."""
        values = np.zeros(self.basis_size + len(self.slaves))
        values[basis] = 1.0
        for offset, (_, parent_a, parent_b, weight_a) in enumerate(self.slaves):
            if parent_a == basis:
                values[self.basis_size + offset] += weight_a
            if parent_b == basis:
                values[self.basis_size + offset] += 1.0 - weight_a
        return values

    def zero_trace_indices(self):
        return [index for index, zero in enumerate(self.knot_zero) if zero]


def connectivity_keys(layout, key_layout):
    """Merge keys of a metal ring's vertices for connect_rings_by_fraction: a knot's key is
    the perimeter fraction of the SAME ROLE (slot) in `key_layout` (the rule's layout at the
    segment's connectivity angle), a slave's key is its corner fraction. The band
    triangulation is then constant over the angles sharing the connectivity angle (the
    fraction-ordered merge flips a quad's diagonal whenever a knot passes a fixed-ring vertex,
    a jump of the hats; see CornerTraceBasisRule ConnectivityAngleDegrees). Raises when the
    key order is not the ring's perimeter order (a knot-corner passage lies between the two
    angles: the triangles would fold)."""
    key_by_slot = {slot: fraction for fraction, kind, slot in key_layout if kind != "slave"}
    keyed = []
    for fraction, kind, slot in layout:
        key = fraction if kind == "slave" else key_by_slot[slot]
        keyed.append((key, fraction, kind, slot))
    keyed.sort(key=lambda vertex: vertex[0])
    fractions = [vertex[1] for vertex in keyed]
    descents = sum(
        1 for previous, following in zip(fractions, fractions[1:] + fractions[:1])
        if following < previous
    )
    if descents > 1:
        raise ValueError(
            "the connectivity angle lies across a knot-corner passage from the ring's angle"
        )
    return keyed


def build_surface(
    radius,
    ring_size,
    metal_thickness,
    overetch_depth,
    angle_degrees=None,
    topology="convex",
    cap_centers=False,
    connectivity_angle_degrees=None,
    crossing_vertices=False,
):
    """The trace surface of a corner coupon (TraceSurface). With an angle, the rings at
    metal_levels follow the trace basis rule (metal_ring_layout); without one (the probe
    surfaces, whose traces vanish on the metal band) every ring is the fixed layout and the
    result carries no slave vertices. With an angle and crossing_vertices the rings at
    metal_levels are the fixed layout plus the two metal-arm crossings (heldout_ring_layout:
    the fine held-out surface, whose trace vanishes on the PEC part of the rings only). With
    a connectivity angle the bands next to the metal rings are triangulated in the merge
    order of the rule's layout at THAT angle (connectivity_keys; the knot positions stay
    those of angle_degrees), else in the perimeter order at angle_degrees."""
    if crossing_vertices and (angle_degrees is None or connectivity_angle_degrees is not None):
        raise ValueError(
            "crossing vertices need an angle and take no connectivity angle (the held-out "
            "surface is not a basis)"
        )
    # Keep the trace triangulation conforming to every fabrication plane that
    # reaches the matching surface. The coupon mesh resolves these intersections
    # so the narrow trace hats across the process zone have active boundary DOFs.
    levels = sorted(
        {
            -radius,
            -radius / 3.0,
            -overetch_depth,
            0.0,
            metal_thickness,
            radius / 3.0,
            radius,
        }
    )
    tolerance = 1.0e-12 * radius
    surface = TraceSurface()
    rings = []  # per ring: list of (fraction, vertex index) in perimeter order
    for half_width, level in (
        [(radius, level) for level in levels]
        + [(radius / 3.0, radius), (radius / 3.0, -radius)]
    ):
        meets_metal = (
            angle_degrees is not None
            and half_width == radius
            and any(abs(level - metal) <= tolerance for metal in metal_levels(metal_thickness))
        )
        layout = (
            (
                heldout_ring_layout(radius, angle_degrees, topology, ring_size)
                if crossing_vertices
                else metal_ring_layout(radius, angle_degrees, topology, ring_size)
            )
            if meets_metal
            else fixed_ring_layout(ring_size)
        )
        if meets_metal and connectivity_angle_degrees is not None:
            keyed = connectivity_keys(
                layout,
                metal_ring_layout(radius, connectivity_angle_degrees, topology, ring_size),
            )
            key_of = {(fraction, slot): key for key, fraction, _, slot in keyed}
        else:
            key_of = {(fraction, slot): fraction for fraction, _, slot in layout}
        points = ring_points(half_width, level, layout, ring_size)
        knot_offset = surface.basis_size
        knot_count = sum(1 for _, kind, _ in layout if kind != "slave")
        surface.knot_points.extend([None] * knot_count)
        surface.knot_zero.extend([None] * knot_count)
        ring = []  # (fraction, knot index or None, merge key) in perimeter order
        pending_slaves = []
        for (fraction, kind, slot), point in zip(layout, points):
            if kind == "slave":
                pending_slaves.append((fraction, point, len(ring)))
                ring.append((fraction, None, key_of[(fraction, slot)]))
                continue
            surface.knot_points[knot_offset + slot] = tuple(point)
            surface.knot_zero[knot_offset + slot] = kind == "zero"
            ring.append((fraction, knot_offset + slot, key_of[(fraction, slot)]))
        if any(point is None for point in surface.knot_points[knot_offset:]):
            raise ValueError("ring layout does not fill every basis slot")
        surface.contour_groups.append(knot_count)
        knot_positions = [
            position for position, (_, index, _) in enumerate(ring) if index is not None
        ]
        for fraction, point, position in pending_slaves:
            before = [p for p in knot_positions if p < position]
            after = [p for p in knot_positions if p > position]
            previous = max(before) if before else max(knot_positions)
            following = min(after) if after else min(knot_positions)
            f_previous = ring[previous][0] - (1.0 if previous > position else 0.0)
            f_following = ring[following][0] + (1.0 if following < position else 0.0)
            weight_a = (f_following - fraction) / (f_following - f_previous)
            surface.slaves.append(
                (tuple(point), ring[previous][1], ring[following][1], weight_a)
            )
        rings.append(ring)
    # Vertex numbering: knots first (basis order), then the slaves in ring (perimeter) order.
    # The bands are merged in KEY order (= the perimeter order without a connectivity angle).
    slave_counter = [surface.basis_size]

    def ring_vertices(ring):
        keyed = []
        for fraction, index, key in ring:
            if index is None:
                index = slave_counter[0]
                slave_counter[0] += 1
            keyed.append((key, index))
        keyed.sort(key=lambda vertex: vertex[0])
        return np.asarray([key for key, _ in keyed]), [index for _, index in keyed]

    resolved_rings = [ring_vertices(ring) for ring in rings]

    def connect(first, second):
        first_fractions, first_indices = resolved_rings[first]
        second_fractions, second_indices = resolved_rings[second]
        local = []
        connect_rings_by_fraction(
            local, 0, first_fractions, 1000000, second_fractions
        )
        for triangle in local:
            surface.triangles.append(
                tuple(
                    second_indices[vertex - 1000000] if vertex >= 1000000 else first_indices[vertex]
                    for vertex in triangle
                )
            )

    outer = len(levels)
    for ring in range(outer - 1):
        connect(ring, ring + 1)
    top_inner, bottom_inner = outer, outer + 1
    connect(outer - 1, top_inner)
    connect(bottom_inner, 0)
    top_indices = resolved_rings[top_inner][1]
    bottom_indices = resolved_rings[bottom_inner][1]
    if cap_centers:
        top_center = surface.basis_size + len(surface.slaves)
        bottom_center = top_center + 1
        surface.slaves.append(((0.0, 0.0, radius), None, None, None))
        surface.slaves.append(((0.0, 0.0, -radius), None, None, None))
        for index in range(ring_size):
            next_index = (index + 1) % ring_size
            surface.triangles.append(
                (top_indices[index], top_indices[next_index], top_center)
            )
            surface.triangles.append(
                (bottom_indices[next_index], bottom_indices[index], bottom_center)
            )
    else:
        local = []
        cap_ring(local, 0, ring_size, False)
        surface.triangles.extend(tuple(top_indices[v] for v in t) for t in local)
        local = []
        cap_ring(local, 0, ring_size, True)
        surface.triangles.extend(tuple(bottom_indices[v] for v in t) for t in local)
    surface.triangles = np.asarray(surface.triangles, dtype=int)
    return surface


def free_hat_pec_support(points, contour_groups, zero_trace_indices, in_metal, slaves=()):
    """The basis gate of the box contour (corner-family review 2026-09-29): along every closed
    ring, the metal footprint's boundary must cross the contour at PEC-constrained knots
    only, so that no free hat (nonzero between its knot and the neighbouring knots) has
    support on the PEC part of the contour. `in_metal(points)` is the closed footprint
    mask at the ring's level. Slave vertices (point, parent_a, parent_b, weight) are
    checked as part of the segment between their parents. Returns the offending (free
    knot, neighbour, metal fraction of the segment) triples as 0-based indices; empty when
    the basis passes. A segment is sampled at 64 interior points (the footprints are
    polygonal / circular: the sampling finds every crossing wider than 1/64 of a segment;
    the exact crossing of a straight arm is placed by arm_crossing_fractions, which this
    check audits)."""
    constrained = set(int(index) for index in zero_trace_indices)
    offending = []
    offset = 0
    samples = (np.arange(64) + 0.5) / 64.0
    slave_between = {}
    for point, parent_a, parent_b, _ in slaves:
        if parent_a is None:
            continue
        slave_between.setdefault((parent_a, parent_b), []).append(np.asarray(point))
        slave_between.setdefault((parent_b, parent_a), []).append(np.asarray(point))

    def polyline(knot, following):
        # The ring segment from knot to following through the slave corners between them.
        corners = slave_between.get((knot, following), [])
        vertices = [points[knot]] + corners + [points[following]]
        # Order the interior corners along the way from knot to following.
        if len(corners) > 1:
            corners = sorted(
                corners, key=lambda c: np.linalg.norm(c - points[knot])
            )
            vertices = [points[knot]] + corners + [points[following]]
        return vertices

    for size in contour_groups:
        for local in range(size):
            knot = offset + local
            following = offset + (local + 1) % size
            if knot in constrained and following in constrained:
                continue
            inside = []
            vertices = polyline(knot, following)
            for a, b in zip(vertices, vertices[1:]):
                segment = a + samples[:, None] * (b - a)
                inside.append(in_metal(segment))
            inside = np.concatenate(inside)
            if np.any(inside):
                free = knot if knot not in constrained else following
                other = following if free == knot else knot
                offending.append((free, other, float(np.mean(inside))))
        offset += size
    return offending


def write_trace_mesh(output, surface):
    """trace-vertices.csv / trace-triangles.csv of a TraceSurface: the knots (basis = the
    1-based knot index; conductor 1 = PEC constrained) followed by the slave vertices (basis
    0, conductor 0, parents parent_a / parent_b 1-based with the weight of parent_a)."""
    points = surface.vertex_points()
    vertex_rows = []
    for index, point in enumerate(points, start=1):
        if index <= surface.basis_size:
            zero = surface.knot_zero[index - 1]
            vertex_rows.append((index, *point, index, 1 if zero else 0, 0, 0, 0.0))
        else:
            _, parent_a, parent_b, weight_a = surface.slaves[index - surface.basis_size - 1]
            vertex_rows.append(
                (index, *point, 0, 0, parent_a + 1, parent_b + 1, weight_a)
            )
    np.savetxt(
        output / "trace-vertices.csv",
        np.asarray(vertex_rows),
        delimiter=",",
        header="vertex,x,y,z,basis,conductor,parent_a,parent_b,weight_a",
        comments="",
        fmt=("%d", "%.16e", "%.16e", "%.16e", "%d", "%d", "%d", "%d", "%.16e"),
    )
    triangle_rows = [
        (index, *(np.asarray(triangle, dtype=int) + 1))
        for index, triangle in enumerate(surface.triangles, start=1)
    ]
    np.savetxt(
        output / "trace-triangles.csv",
        np.asarray(triangle_rows),
        delimiter=",",
        header="triangle,vertex_i,vertex_j,vertex_k",
        comments="",
        fmt="%d",
    )


def write_basis(output, surface):
    """basis-points.csv (the knots) and one hat trace per knot on every triangle of the
    trace surface (slave vertices carry the interpolated hat value)."""
    np.savetxt(
        output / "basis-points.csv",
        np.asarray(surface.knot_points),
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )
    trace_directory = output / "traces"
    trace_directory.mkdir(parents=True, exist_ok=True)
    for stale_trace in trace_directory.glob("basis-*.csv"):
        stale_trace.unlink()
    points = surface.vertex_points()
    paths = []
    for basis in range(surface.basis_size):
        path = trace_directory / f"basis-{basis + 1:03d}.csv"
        write_surface_trace(path, points, surface.triangles, surface.hat_values(basis))
        paths.append(path)
    return paths


def write_surface_trace(path, points, triangles, values):
    rows = []
    for triangle_index, triangle in enumerate(triangles, start=1):
        for vertex in triangle:
            rows.append((*points[vertex], values[vertex], triangle_index))
    np.savetxt(
        path,
        np.asarray(rows),
        delimiter=",",
        header="x,y,z,V,triangle",
        comments="",
        fmt=("%.16e", "%.16e", "%.16e", "%.16e", "%d"),
    )


def metal_band_cutoff(points, radius, metal_thickness):
    """Smoothly suppress a trace on both thin and fabricated PEC cuts (the convergence
    probes: the whole metal band, PEC or not)."""
    transition = radius / 3.0
    distance = np.maximum(-points[:, 2], points[:, 2] - metal_thickness)
    coordinate = np.clip(distance / transition, 0.0, 1.0)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def pec_contour_distance(points, radius, angle_degrees, topology, metal_thickness):
    """Distance (in length units) from box-surface points to the PEC part of the matching
    box, the metal band z in [0, MetalThickness] over the metal arc of the perimeter (the
    swept footprint of both coupons: the thin sheet's cut and the fabricated slab's foot and
    top edge cross the rings at the same arm crossings; a sloped sidewall lies on the metal
    side of them): the hypotenuse of the perimeter-arc distance to the metal arc
    (perimeter_arc_distance x 8 R) and the vertical distance to the band. Points off the
    perimeter (the cap rings and centres at z = -/+ R) take the vertical distance alone."""
    metal, _ = metal_arc_fractions(radius, angle_degrees, topology)
    tolerance = 1.0e-12 * radius
    lateral = np.zeros(len(points))
    for index, point in enumerate(points):
        if max(abs(point[0]), abs(point[1])) >= radius - tolerance:
            lateral[index] = 8.0 * radius * perimeter_arc_distance(
                square_perimeter_fraction(radius, point), metal
            )
    vertical = np.maximum.reduce(
        (-points[:, 2], points[:, 2] - metal_thickness, np.zeros(len(points)))
    )
    return np.hypot(lateral, vertical)


def heldout_cutoff(points, radius, angle_degrees, topology, metal_thickness):
    """The held-out trace's cutoff (USER decision 149 (6), option (c)): the smoothstep over
    R / 3 of the distance to the PEC part of the box, so it vanishes exactly on the PEC part
    of the metal rings and nowhere else — every free knot of the trace basis rule, on the
    metal rings included, is excited and the self-check judges the basis (the former cutoff
    suppressed the whole metal band, so the held-out coefficients were exactly zero on both
    metal rings and the check never saw the free knots next to the metal). The same form as
    the 2D generators' distance-to-the-cut cutoff (cpw2d/generate_edge_response.py)."""
    distance = pec_contour_distance(points, radius, angle_degrees, topology, metal_thickness)
    coordinate = np.clip(distance / (radius / 3.0), 0.0, 1.0)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def heldout_polynomial(points, radius):
    """The held-out trace before its cutoff: a fixed low-order polynomial of the box
    coordinates over R (the same polynomial under the former band cutoff and option (c),
    so a recorded coefficient file can be classified by its cutoff alone)."""
    x = points[:, 0] / radius
    y = points[:, 1] / radius
    z = points[:, 2] / radius
    return 0.35 + 0.20 * x - 0.15 * y + 0.10 * z + 0.08 * x * y + 0.06 * z * z


def heldout_potential(points, radius, metal_thickness, angle_degrees, topology):
    # Both coupons have conductor cuts in the matching surface at z = 0, while the
    # fabricated cut extends to z = metal_thickness. One smooth trace compatible with both
    # cuts (zero on the swept PEC part of the box) avoids an order-dependent Dirichlet jump
    # where the matching and grounded boundaries meet.
    return heldout_cutoff(
        points, radius, angle_degrees, topology, metal_thickness
    ) * heldout_polynomial(points, radius)


def convergence_probe_potentials(points, radius, metal_thickness):
    """Return a small low-order trace space for response convergence tests."""
    x = points[:, 0] / radius
    y = points[:, 1] / radius
    z = points[:, 2] / radius
    cutoff = metal_band_cutoff(points, radius, metal_thickness)
    probes = (
        ("common", "cutoff", np.ones_like(x)),
        ("x-linear", "cutoff*x/R", x),
        ("y-linear", "cutoff*y/R", y),
        ("z-linear", "cutoff*z/R", z),
        ("xy-mixed", "cutoff*x*y/R^2", x * y),
        ("z-quadratic", "cutoff*z^2/R^2", z * z),
    )
    values = np.column_stack([cutoff * probe[2] for probe in probes])
    if np.linalg.matrix_rank(values) != len(probes):
        raise ValueError("corner convergence probes are linearly dependent")
    return [
        (name, expression, values[:, index])
        for index, (name, expression, _) in enumerate(probes)
    ]


# The angle-interpolated corner family's straight anchor (USER decision 121 (C)): a sharp
# "corner" of exactly 180 degrees, the straight edge through the corner box on the same basis.
STRAIGHT_ANCHOR_TOLERANCE_DEGREES = 1.0e-9


def is_straight_anchor(angle_degrees):
    return abs(angle_degrees - 180.0) <= STRAIGHT_ANCHOR_TOLERANCE_DEGREES


def corner_frame(angle_degrees):
    angle = np.deg2rad(angle_degrees)
    if not 0.0 < angle <= np.pi + np.deg2rad(STRAIGHT_ANCHOR_TOLERANCE_DEGREES):
        raise ValueError("corner angle must lie in (0, 180] degrees")
    first_normal = np.array([0.0, 1.0])
    second_normal = np.array([np.sin(angle), -np.cos(angle)])
    bisector = np.array([np.cos(0.5 * angle), np.sin(0.5 * angle)])
    return angle, first_normal, second_normal, bisector


def corner_center(angle_degrees, corner_radius):
    angle, _, _, bisector = corner_frame(angle_degrees)
    return corner_radius * bisector / np.sin(0.5 * angle)


def metal_footprint_mask(
    points,
    radius,
    angle_degrees,
    corner_radius,
    offset,
    topology,
):
    tolerance = 1.0e-12 * radius
    _, first_normal, second_normal, _ = corner_frame(angle_degrees)
    coordinates = points[:, :2]
    first_distance = coordinates @ first_normal
    second_distance = coordinates @ second_normal
    in_wedge = (first_distance >= offset - tolerance) & (
        second_distance >= offset - tolerance
    )
    if corner_radius > 0.0:
        arc_radius = corner_radius - offset
        if arc_radius <= 0.0:
            raise ValueError(
                "corner offset must be smaller than the plan-view corner radius"
            )
        center = corner_center(angle_degrees, corner_radius)
        rounded_wedge = (
            (first_distance >= corner_radius - tolerance)
            | (second_distance >= corner_radius - tolerance)
            | (
                np.sum((coordinates - center) ** 2, axis=1)
                <= arc_radius**2 + tolerance**2
            )
        )
        in_wedge &= rounded_wedge
    if topology == "convex":
        in_metal = in_wedge
    else:
        on_corner_boundary = (
            (np.abs(first_distance - offset) <= tolerance)
            & (second_distance >= corner_radius - tolerance)
        ) | (
            (np.abs(second_distance - offset) <= tolerance)
            & (first_distance >= corner_radius - tolerance)
        )
        in_metal = ~in_wedge | on_corner_boundary
    return in_metal


def thin_metal_mask(points, radius, angle_degrees, corner_radius, topology):
    tolerance = 1.0e-12 * radius
    return (np.abs(points[:, 2]) <= tolerance) & metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )


def pec_trace_mask(
    points,
    radius,
    angle_degrees,
    corner_radius,
    metal_thickness,
    sidewall_angle,
    topology,
):
    tolerance = 1.0e-12 * radius
    bottom = (np.abs(points[:, 2]) <= tolerance) & metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )
    pullback = metal_thickness / np.tan(np.deg2rad(sidewall_angle))
    top_offset = pullback if topology == "convex" else -pullback
    top_footprint = metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        top_offset,
        topology,
    )
    swept_footprint = top_footprint | metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )
    top = (
        np.abs(points[:, 2] - metal_thickness) <= tolerance
    ) & swept_footprint
    return bottom | top


def dielectric(
    index, attributes, interface_type, radius, thickness, permittivity
):
    _, loss_tangent = INTERFACES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "LocalizeEdgeEnergy": True,
        "SaveLocalEdgeEnergy": False,
        "EdgeAttributes": [SA_ATTRIBUTE],
        "EdgeExcludeAttributes": [MATCHING_SURFACE_ATTRIBUTE],
        "EdgeFrameNormal": [0.0, 0.0, 1.0],
        "EdgeDistances": [radius],
    }


def make_config(
    output,
    name,
    mesh,
    traces,
    radius,
    order,
    fabricated,
    substrate_permittivity,
    interface_layers,
):
    interfaces = (
        [
            dielectric(1, [3], "SA", radius, *interface_layers["SA"]),
            dielectric(2, [2], "MS", radius, *interface_layers["MS"]),
            dielectric(3, [4], "MA", radius, *interface_layers["MA"]),
        ]
        if fabricated
        else [
            dielectric(1, [3], "SA", radius, *interface_layers["SA"]),
            dielectric(2, [2], "MS", radius, *interface_layers["MS"]),
            dielectric(3, [2], "MA", radius, *interface_layers["MA"]),
        ]
    )
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
        },
        "Model": {
            "Mesh": str(mesh),
            "L0": 1.0e-6,
            "Refinement": {"MaxIts": 0},
        },
        "Domains": {
            "Materials": [
                {"Attributes": [1], "Permittivity": substrate_permittivity},
                {"Attributes": [2], "Permittivity": 1.0},
            ],
            "Postprocessing": {
                "Energy": [
                    {"Index": 1, "Attributes": [1]},
                    {"Index": 2, "Attributes": [2]},
                ]
            },
        },
        "Boundaries": {
            "Ground": {"Attributes": [2, 4] if fabricated else [2]},
            "PrescribedPotential": [
                {
                    "Index": index,
                    "Attributes": [1],
                    "DataFile": str(trace),
                }
                for index, trace in enumerate(traces, start=1)
            ],
            "Postprocessing": {"Dielectric": interfaces},
        },
        "Solver": {
            "Order": order,
            "Electrostatic": {
                "Save": 0,
                "ResponseMatrix": True,
                "AggregateResponseMatrix": True,
            },
            "Linear": {
                "Type": "BoomerAMG",
                "KSPType": "CG",
                "Tol": 1.0e-10,
                "MaxIts": 1000,
                # These generated coupon solves never adapt the mesh.
                "EstimatorTol": 5.0e-1,
                "EstimatorMaxIts": 5,
                "EstimatorMG": True,
            },
        },
    }


def write_library(
    output,
    radius,
    angle_degrees,
    corner_radius,
    contour_groups,
    zero_trace_indices,
    topology,
    metal_thickness,
    overetch_depth,
    sidewall_angle,
    top_rounding,
    trench_rounding,
    substrate_permittivity,
    interface_layers,
    trace_basis=None,
):
    topology_name = f"{topology.capitalize()}Corner"
    model_name = f"{topology}-corner-{angle_degrees:g}deg"
    if corner_radius > 0.0:
        model_name += f"-r{corner_radius:g}um"
    if trace_basis is not None and "ConnectivityAngleDegrees" in trace_basis:
        # A segment node (one coupon per side of a knot-corner passage shares the angle):
        # the name carries the segment's connectivity angle (model names are unique).
        model_name += f"-c{trace_basis['ConnectivityAngleDegrees']:g}"
    reference = (
        [*corner_center(angle_degrees, corner_radius), 0.0]
        if topology == "convex" and corner_radius > 0.0
        else [0.0, 0.0, 0.0]
    )
    corner_radius_tolerance = max(
        0.02 * corner_radius,
        1.0e-3 * radius,
    )
    model = {
        "Name": model_name,
        "Topology": topology_name,
        "Angle": angle_degrees,
        # Corner family records (USER decision 121 (C)): the angle and the convexity by name
        # (the same values as Angle / Topology; the family interpolates sharp coupons of one
        # Convexity in the turn 180 - AngleDegrees; Angle 180 = the family's straight anchor).
        "AngleDegrees": angle_degrees,
        "Convexity": topology.capitalize(),
        "AngleTolerance": 2.0,
        "CornerRadius": corner_radius,
        "CornerRadiusTolerance": corner_radius_tolerance,
        "Reference": reference,
        "FabricatedMatrix": "postpro/fabricated/domain-response-matrix.csv",
        "ThinMatrix": "postpro/thin/domain-response-matrix.csv",
        "FabricatedSurfaceMatrix":
            "postpro/fabricated/surface-response-matrix-aggregate.csv",
        "ThinSurfaceMatrix":
            "postpro/thin/surface-response-matrix-aggregate.csv",
        "BasisPoints": "basis-points.csv",
        "TraceMesh": {
            "Vertices": "trace-vertices.csv",
            "Triangles": "trace-triangles.csv",
        },
        "ContourGroups": contour_groups,
        "Interfaces": [
            {"Type": "SA", "Coupon": 1},
            {"Type": "MS", "Coupon": 2},
            {"Type": "MA", "Coupon": 3},
        ],
    }
    if zero_trace_indices:
        model["ZeroTraceIndices"] = zero_trace_indices
    # The trace basis rule (corner-family review 2026-09-29 and the supervisor decision on
    # it): angle-independent knot semantics on the rings that meet the metal, recorded for
    # the runtime's like-to-like check of the family (MatchCornerFamily) and the audit.
    if trace_basis is not None:
        model["TraceBasis"] = trace_basis
    library = {
        "Version": 3,
        "TraceLiftVersion": 2,
        "Name": (
            "100nm-metal-50nm-overetch-10nm-rounding-"
            f"{topology}-corner-r{corner_radius:g}um-prototype"
        ),
        "MatchingRadius": radius,
        "Fabrication": {
            "LengthUnit": "um",
            "MetalThickness": metal_thickness,
            "OveretchDepth": overetch_depth,
            "SidewallAngle": sidewall_angle,
            "TopRounding": top_rounding,
            "TrenchRounding": trench_rounding,
            "SubstratePermittivity": substrate_permittivity,
            "InterfaceLayers": {
                interface_type: {
                    "Thickness": layer[0],
                    "Permittivity": layer[1],
                }
                for interface_type, layer in interface_layers.items()
            },
        },
        "Models": [model],
    }
    path = output / "process-library.json"
    path.write_text(json.dumps(library, indent=2) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--angle", type=float, default=90.0)
    parser.add_argument(
        "--connectivity-angle",
        type=float,
        default=None,
        help="segment connectivity angle of a corner family node (degrees; the band merge "
        "order of the rule's layout at this angle; recorded as TraceBasis "
        "ConnectivityAngleDegrees); default: the perimeter order at --angle (legacy node)",
    )
    parser.add_argument("--corner-radius", type=float, default=0.0)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch-depth", type=float, default=0.05)
    parser.add_argument("--sidewall-angle", type=float, default=80.0)
    parser.add_argument("--top-rounding", type=float, default=0.01)
    parser.add_argument("--trench-rounding", type=float, default=0.01)
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    parser.add_argument("--sa-thickness", type=float, default=0.002)
    parser.add_argument("--sa-permittivity", type=float, default=4.0)
    parser.add_argument("--ms-thickness", type=float, default=0.002)
    parser.add_argument("--ms-permittivity", type=float, default=11.47)
    parser.add_argument("--ma-thickness", type=float, default=0.002)
    parser.add_argument("--ma-permittivity", type=float, default=10.0)
    parser.add_argument(
        "--topology", choices=("convex", "concave"), default="convex"
    )
    args = parser.parse_args()
    if args.radius <= 0.0:
        parser.error("--radius must be positive")
    if not (0.0 < args.angle < 180.0 or is_straight_anchor(args.angle)):
        parser.error("--angle must lie strictly between zero and 180 degrees (180 = the straight anchor)")
    if is_straight_anchor(args.angle) and args.corner_radius > 0.0:
        parser.error("the straight anchor (--angle 180) has no corner radius")
    if not 0.0 <= args.corner_radius < args.radius:
        parser.error("--corner-radius must lie in [0, radius)")
    if args.connectivity_angle is not None and (
        not 0.0 < args.connectivity_angle < 180.0 or args.corner_radius > 0.0
    ):
        parser.error(
            "--connectivity-angle must lie strictly between zero and 180 degrees and is a "
            "sharp corner family option"
        )
    tangent_distance = (
        0.0 if is_straight_anchor(args.angle)
        else args.corner_radius / np.tan(0.5 * np.deg2rad(args.angle))
    )
    if args.corner_radius > 0.0 and tangent_distance >= args.radius:
        parser.error(
            "rounded-corner tangency points must lie inside the matching box"
        )
    if args.order < 1:
        parser.error("--order must be positive")
    if not 0.0 < args.metal_thickness < args.radius:
        parser.error("--metal-thickness must lie between zero and the radius")
    if not 0.0 <= args.overetch_depth < args.radius:
        parser.error("--overetch-depth must be nonnegative and smaller than the radius")
    if not 0.0 < args.sidewall_angle <= 90.0:
        parser.error("--sidewall-angle must lie in (0, 90]")
    if not 0.0 <= args.top_rounding < args.metal_thickness:
        parser.error(
            "--top-rounding must be nonnegative and smaller than the metal thickness"
        )
    if not 0.0 <= args.trench_rounding <= args.overetch_depth:
        parser.error(
            "--trench-rounding must lie between zero and the overetch depth"
        )
    pullback = args.metal_thickness / np.tan(
        np.deg2rad(args.sidewall_angle)
    )
    if (
        args.corner_radius > 0.0
        and args.topology == "convex"
        and pullback >= args.corner_radius
    ):
        parser.error(
            "--sidewall-angle gives a pullback larger than the corner radius"
        )
    trench_pullback = args.overetch_depth / np.tan(
        np.deg2rad(args.sidewall_angle)
    )
    if (
        args.corner_radius > 0.0
        and args.topology == "concave"
        and trench_pullback >= args.corner_radius
    ):
        parser.error(
            "--sidewall-angle gives a trench pullback larger than the corner radius"
        )
    material_values = (
        args.substrate_permittivity,
        args.sa_thickness,
        args.sa_permittivity,
        args.ms_thickness,
        args.ms_permittivity,
        args.ma_thickness,
        args.ma_permittivity,
    )
    if any(not np.isfinite(value) or value <= 0.0 for value in material_values):
        parser.error("substrate and interface-layer properties must be finite and positive")
    interface_layers = {
        "SA": (args.sa_thickness, args.sa_permittivity),
        "MS": (args.ms_thickness, args.ms_permittivity),
        "MA": (args.ma_thickness, args.ma_permittivity),
    }

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    thin_mesh = args.thin_mesh.expanduser().resolve()
    fabricated_mesh = args.fabricated_mesh.expanduser().resolve()
    for mesh in (thin_mesh, fabricated_mesh):
        if not mesh.is_file():
            raise FileNotFoundError(mesh)

    surface = build_surface(
        args.radius,
        args.ring_size,
        args.metal_thickness,
        args.overetch_depth,
        angle_degrees=args.angle,
        topology=args.topology,
        connectivity_angle_degrees=args.connectivity_angle,
    )
    points = np.asarray(surface.knot_points)
    contour_groups = surface.contour_groups
    traces = write_basis(output, surface)

    def pec_mask(sample_points):
        return pec_trace_mask(
            sample_points,
            args.radius,
            args.angle,
            args.corner_radius,
            args.metal_thickness,
            args.sidewall_angle,
            args.topology,
        )

    # The zero set of the rule (crossings + metal-interior knots) must be exactly the set of
    # knots on the PEC footprint of both coupons: the rule and the geometry agree.
    pec_knots = pec_mask(points)
    rule_zero = np.asarray(surface.knot_zero, dtype=bool)
    if not np.array_equal(pec_knots, rule_zero):
        raise ValueError(
            "trace basis rule and PEC footprint disagree on the zero set: rule "
            f"{(np.flatnonzero(rule_zero) + 1).tolist()} vs footprint "
            f"{(np.flatnonzero(pec_knots) + 1).tolist()}"
        )
    zero_trace_indices = (np.flatnonzero(pec_knots) + 1).tolist()
    if not zero_trace_indices:
        raise ValueError("corner coupon has no PEC-constrained trace knots")
    # The basis gate: every crossing of the metal edge with a box ring is a PEC knot, so no
    # free hat has support on the PEC part of the contour (the corner-family review's root
    # cause: a free hat across the metal edge imposes conflicting Dirichlet data and puts a
    # near-singular field into the MS / MA layers of the fabricated coupon).
    offending = free_hat_pec_support(
        points, contour_groups, np.flatnonzero(pec_knots), pec_mask, surface.slaves
    )
    if offending:
        raise ValueError(
            "free trace hats with support on the PEC part of the box contour (free knot, "
            f"neighbour, metal fraction of the segment; 1-based): "
            f"{[(free + 1, fixed + 1, round(fraction, 3)) for free, fixed, fraction in offending]}"
        )
    write_trace_mesh(output, surface)
    for name, mesh, fabricated in (
        ("thin", thin_mesh, False),
        ("fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            traces,
            args.radius,
            args.order,
            fabricated,
            args.substrate_permittivity,
            interface_layers,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")
    library = write_library(
        output,
        args.radius,
        args.angle,
        args.corner_radius,
        contour_groups,
        zero_trace_indices,
        args.topology,
        args.metal_thickness,
        args.overetch_depth,
        args.sidewall_angle,
        args.top_rounding,
        args.trench_rounding,
        args.substrate_permittivity,
        interface_layers,
        trace_basis_rule(args.ring_size, args.connectivity_angle),
    )

    fine_surface = build_surface(
        args.radius,
        max(16, 4 * args.ring_size),
        args.metal_thickness,
        args.overetch_depth,
        angle_degrees=args.angle,
        topology=args.topology,
        cap_centers=True,
        crossing_vertices=True,
    )
    fine_points = fine_surface.vertex_points()
    fine_triangles = fine_surface.triangles
    fine_values = heldout_potential(
        fine_points,
        args.radius,
        args.metal_thickness,
        args.angle,
        args.topology,
    )
    # Option (c) of USER decision 149 (6): the fine trace is zero at every PEC vertex (the
    # crossings are vertices, so it vanishes on the whole PEC part of the box) and nonzero
    # at every vertex off the PEC part, the free knots of the metal rings included.
    fine_pec = pec_mask(fine_points)
    fine_cutoff = heldout_cutoff(
        fine_points, args.radius, args.angle, args.topology, args.metal_thickness
    )
    if np.any(fine_values[fine_pec] != 0.0) or np.any(fine_cutoff[~fine_pec] <= 0.0):
        raise ValueError(
            "held-out trace does not vanish exactly on the PEC part of the box and nowhere "
            "else"
        )
    heldout_trace = output / "heldout-trace.csv"
    write_surface_trace(heldout_trace, fine_points, fine_triangles, fine_values)
    heldout_coefficients = heldout_potential(
        points,
        args.radius,
        args.metal_thickness,
        args.angle,
        args.topology,
    )
    coarse_cutoff = heldout_cutoff(
        points, args.radius, args.angle, args.topology, args.metal_thickness
    )
    if np.any(heldout_coefficients[pec_knots] != 0.0) or np.any(
        coarse_cutoff[~pec_knots] <= 0.0
    ):
        raise ValueError(
            "held-out coefficients do not vanish exactly on the zero set and excite every "
            "free knot"
        )
    np.savetxt(
        output / "heldout-coefficients.csv",
        heldout_coefficients,
        delimiter=",",
        header="coefficient_V",
        comments="",
        fmt="%.16e",
    )
    for name, mesh, fabricated in (
        ("heldout-thin", thin_mesh, False),
        ("heldout-fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            [heldout_trace],
            args.radius,
            args.order,
            fabricated,
            args.substrate_permittivity,
            interface_layers,
        )
        config["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        config["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = False
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")

    probe_directory = output / "convergence-probes"
    probe_directory.mkdir(parents=True, exist_ok=True)
    for stale_probe in probe_directory.glob("probe-*.csv"):
        stale_probe.unlink()
    convergence_probes = convergence_probe_potentials(
        fine_points,
        args.radius,
        args.metal_thickness,
    )
    probe_paths = []
    probe_metadata = []
    for index, (name, expression, values) in enumerate(
        convergence_probes, start=1
    ):
        path = probe_directory / f"probe-{index:02d}-{name}.csv"
        write_surface_trace(path, fine_points, fine_triangles, values)
        probe_paths.append(path)
        probe_metadata.append(
            {
                "Index": index,
                "Name": name,
                "Expression": expression,
            }
        )
    (output / "probe-manifest.json").write_text(
        json.dumps({"Version": 1, "Probes": probe_metadata}, indent=2) + "\n"
    )
    for name, mesh, fabricated in (
        ("probe-thin", thin_mesh, False),
        ("probe-fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            probe_paths,
            args.radius,
            args.order,
            fabricated,
            args.substrate_permittivity,
            interface_layers,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")

    print(f"Generated {len(points)} spatial basis traces")
    print(f"Generated {len(convergence_probes)} convergence probe traces")
    print(f"ContourGroups: {contour_groups}")
    if zero_trace_indices:
        print(f"ZeroTraceIndices: {zero_trace_indices}")
    print(output / "thin.json")
    print(output / "fabricated.json")
    print(library)
    print(output / "heldout-thin.json")
    print(output / "heldout-fabricated.json")
    print(output / "probe-thin.json")
    print(output / "probe-fabricated.json")


if __name__ == "__main__":
    main()
