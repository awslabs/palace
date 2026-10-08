#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the identification audit on a synthetic tiny mesh with known answers."""

import contextlib
import io
import math
import json
import os
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import audit, manifest as M, perimeter as P  # noqa: E402
from surface_response_identification.msh2 import read_msh2  # noqa: E402


def write_msh2(path, nodes, elements, names, binary):
    """nodes: list of (x, y, z) (tags 1..N); elements: list of (type, physical, node tags)."""
    if binary:
        out = bytearray(b"$MeshFormat\n2.2 1 8\n" + struct.pack("<i", 1) + b"\n$EndMeshFormat\n")
    else:
        out = bytearray(b"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    out += f"$PhysicalNames\n{len(names)}\n".encode()
    for dimension, tag, name in names:
        out += f'{dimension} {tag} "{name}"\n'.encode()
    out += b"$EndPhysicalNames\n"
    out += f"$Nodes\n{len(nodes)}\n".encode()
    for tag, (x, y, z) in enumerate(nodes, start=1):
        out += struct.pack("<iddd", tag, x, y, z) if binary else f"{tag} {x!r} {y!r} {z!r}\n".encode()
    if binary:
        out += b"\n"
    out += b"$EndNodes\n"
    out += f"$Elements\n{len(elements)}\n".encode()
    if binary:
        by_type = {}
        for element_type, physical, tags in elements:
            by_type.setdefault(element_type, []).append((physical, tags))
        number = 1
        for element_type, entries in by_type.items():
            out += struct.pack("<3i", element_type, len(entries), 2)
            for physical, tags in entries:
                out += struct.pack(f"<{3 + len(tags)}i", number, physical, physical, *tags)
                number += 1
        out += b"\n"
    else:
        for number, (element_type, physical, tags) in enumerate(elements, start=1):
            out += f"{number} {element_type} 2 {physical} {physical} {' '.join(map(str, tags))}\n".encode()
    out += b"$EndElements\n"
    with open(path, "wb") as target:
        target.write(bytes(out))


def island_mesh():
    """A 4 x 2 metal rectangle (attribute 5) made of 4 triangles on the z = 0 plane, a
    second metal rectangle 2 units away at x = 6..8 (a strip pair at separation exactly 2),
    SA faces (8) around them, and a vertical metal wall (attribute 5) standing on the first
    rectangle's interior line x = 2 with a span at z = 1 (non-planar wall edges, a fold at the
    wall top, a nonmanifold foot line, and with R = 0.5 cross-layer zones on the parts of the
    first rectangle's long edges within 2R = 1 of the wall and on the span's short edges).
    One tetrahedron references the SA face so it counts as exterior; the 'outer' face 3 is
    exterior and coincides with the first rectangle's left edge (truncation). A far node
    keeps the span off the bounding box of the mesh (a PEC simulation box is not metal)."""
    nodes = [
        (0, 0, 0), (2, 0, 0), (4, 0, 0),  # 1 2 3
        (0, 2, 0), (2, 2, 0), (4, 2, 0),  # 4 5 6
        (6, 0, 0), (8, 0, 0), (6, 2, 0), (8, 2, 0),  # 7 8 9 10
        (2, 0, 1), (2, 2, 1), (3, 0, 1), (3, 2, 1),  # 11 12 13 14  wall top + span
        (0, 0, -3), (0, 2, -3),  # 15 16 outer box below the left edge
        (4, -2, 0), (6, -2, 0),  # 17 18 SA below the gap
        (5, 1, -1),  # 19 apex of the tetrahedra
        (4, 0, 3),  # 20 unused: keeps the span off the bounding box
    ]
    elements = [
        (2, 5, (1, 2, 5)), (2, 5, (1, 5, 4)), (2, 5, (2, 3, 6)), (2, 5, (2, 6, 5)),  # island A
        (2, 5, (7, 8, 10)), (2, 5, (7, 10, 9)),  # island B
        (2, 5, (2, 5, 12)), (2, 5, (2, 12, 11)),  # wall x = 2, z 0..1
        (2, 5, (11, 12, 14)), (2, 5, (11, 14, 13)),  # span z = 1
        (2, 8, (3, 7, 9)), (2, 8, (3, 9, 6)),  # SA in the gap
        (2, 8, (3, 17, 18)), (2, 8, (3, 18, 7)),  # SA below the gap
        (2, 3, (1, 4, 16)), (2, 3, (1, 16, 15)),  # outer, exterior, on the left edge
        (4, 1, (3, 7, 9, 19)), (4, 1, (3, 9, 6, 19)),  # tets making the SA faces exterior
        (4, 1, (1, 4, 16, 19)),  # tet making the outer face exterior
    ]
    names = [(2, 3, "outer"), (2, 5, "metal"), (2, 8, "substrate_air"), (3, 1, "substrate")]
    return nodes, elements, names


CONFIG = {
    "Model": {"L0": 1.0e-6},
    "Boundaries": {
        "PEC": {"Attributes": [5]},
        "Postprocessing": {
            "Dielectric": [
                {"Index": 1, "Attributes": [8], "Type": "SA", "EdgeDistances": [0.5], "AutomaticEdges": True},
                {"Index": 2, "Attributes": [5], "Type": "MS", "EdgeDistances": [0.5], "AutomaticEdges": True},
            ]
        },
    },
    "Solver": {"SurfaceResponseCorrection": {"Library": "x", "TargetInterfaces": [1, 2]}},
}


def make_manifest(requirements, metal_segments=None, chains=None):
    manifest = {
        "Version": 1,
        "Complete": all(r.get("Status") != "Missing" for r in requirements),
        "LengthUnit": "mesh",
        "Library": {"Path": "/x/lib.json", "Name": "x", "MatchingRadius": 0.5, "DecisionQuantization": {"LengthRelativeToMatchingRadius": 1e-8, "Direction": 1e-12}},
        "MeshDimension": 3,
        "Maxwell": True,
        "Summary": {},
        "Requirements": requirements,
    }
    if metal_segments is not None:
        manifest["Statistics"] = {"Geometry": {"MetalSegments": metal_segments, "PhysicalChains": chains}}
    return manifest


def requirement(topology, count, length, geometry, status="Exact", **extra):
    entry = {"Dimension": 3, "Topology": topology, "Status": status, "Geometry": geometry, "Interfaces": [{"Slot": 0, "Type": "SA", "Target": 1}], "BoundaryCondition": {"Type": "PEC"}, "Count": count, "TotalEdgeLength": length}
    entry.update(extra)
    return entry


def subdivide_closed_polygon(points, spacing):
    """Collinear vertices inserted on every edge of a closed polygon: pieces of at most `spacing`."""
    out = []
    for i, a in enumerate(points):
        b = points[(i + 1) % len(points)]
        pieces = max(1, int(math.ceil(np.linalg.norm(b - a) / spacing - 1.0e-9)))
        out.extend(a + (b - a) * k / pieces for k in range(pieces))
    return out


def closed_metal_perimeter(points, R):
    """One closed PEC / MA chain in the plane z = 0 through `points` (counter-clockwise, metal inside), vertices classified."""
    vertices = [P.PerimeterVertex(point=np.array([p[0], p[1], 0.0])) for p in points]
    edges = []
    n = len(points)
    for i in range(n):
        j = (i + 1) % n
        d = vertices[j].point - vertices[i].point
        edges.append(P.PerimeterEdge(vertices=(i, j), length=float(np.linalg.norm(d)), attributes=(5,), kind="PHYSICAL", conductors=("PEC",), interfaces=((0, "MA"),), inward=np.array([-d[1], d[0], 0.0]) / np.linalg.norm(d), chain=0))
        vertices[i].edges.append(i)
        vertices[j].edges.append(i)
    perimeter = P.Perimeter(vertices=vertices, edges=edges, process_normal=np.array([0.0, 0.0, 1.0]), planes=[0.0], chains=1)
    P.classify_vertices(perimeter, R)
    return perimeter


def arc_reading(points, R):
    """The mirror's arc set of the closed polygon as sorted (radius, turn, joints, rounded) tuples."""
    arcs = P.arc_groups(closed_metal_perimeter(points, R), R)
    return sorted((round(a["Radius"], 6), round(a["TurnDegrees"], 6), len(a["Joints"]), a["Rounded"]) for a in arcs)


class PerimeterTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.nodes, self.elements, self.names = island_mesh()

    def tearDown(self):
        self.directory.cleanup()

    def mesh(self, binary):
        path = os.path.join(self.directory.name, "island.msh2")
        write_msh2(path, self.nodes, self.elements, self.names, binary)
        return read_msh2(path)

    def test_reader_binary_and_ascii_agree(self):
        a, b = self.mesh(True), self.mesh(False)
        np.testing.assert_array_equal(a.coordinates, b.coordinates)
        np.testing.assert_array_equal(a.corner_indices(2), b.corner_indices(2))
        np.testing.assert_array_equal(a.physical_tags(4), b.physical_tags(4))
        self.assertEqual(a.physical_names[(2, 5)], "metal")

    def test_perimeter_classes_and_lengths(self):
        perimeter = P.extract_perimeter(self.mesh(True), CONFIG, radius=0.5)
        np.testing.assert_allclose(np.abs(perimeter.process_normal), [0, 0, 1])
        self.assertEqual(perimeter.primary_plane, 0)
        self.assertEqual(sorted(perimeter.planes), [0.0, 1.0])
        by_kind = {}
        for e in perimeter.edges:
            by_kind.setdefault(e.kind, []).append(e)
        # Island A: bottom 4 + top 4 + right 2 one-sided planar (5 edges), left edge 2
        # (truncation), interior line x = 2 (nonmanifold: two island faces + the wall). The
        # span at z = 1 is planar metal on its own plane: far edge 2 + two short edges 1
        # (3 PHYSICAL edges). Island B: 4 edges, 8.
        self.assertAlmostEqual(sum(e.length for e in by_kind["PHYSICAL"]), 10.0 + 4.0 + 8.0)
        self.assertEqual(len(by_kind["PHYSICAL"]), 12)
        self.assertAlmostEqual(sum(e.length for e in by_kind["TRUNCATION"]), 2.0)
        self.assertAlmostEqual(sum(e.length for e in by_kind["NONMANIFOLD"]), 2.0)
        # Wall vertical edges 2 x 1 (one-sided faces off the process plane); the wall top is
        # a fold with the span.
        self.assertAlmostEqual(sum(e.length for e in by_kind["NONPLANAR"]), 2.0)
        self.assertAlmostEqual(sum(e.length for e in by_kind["FOLD"]), 2.0)
        # Cross-layer zones (R = 0.5, reach 1): the parts of A's long edges within 1 of the
        # wall (x in (1, 3): 2 + 2) and the span's short edges (within 1 of the wall: 1 + 1);
        # A's right edge (2 away) and the span's far edge (exactly 1 away) are not zones.
        self.assertAlmostEqual(perimeter.length("CROSS_LAYER"), 6.0)
        self.assertEqual(sum(1 for e in by_kind["PHYSICAL"] if e.cross_layer), 6)
        # Conductors by edge connectivity: A with its wall and span, B.
        self.assertEqual(perimeter.components, 2)
        self.assertEqual({e.component for e in by_kind["PHYSICAL"]}, {0, 1})
        # Interfaces: the gap-facing edges of A and B carry SA + MS, the others MS only.
        signatures = sorted({e.interfaces for e in by_kind["PHYSICAL"]})
        self.assertIn(((1, "SA"), (2, "MS")), signatures)
        # Vertices: A's left joints (0,0),(0,2) have one physical edge (ENDPOINT); the wall
        # foot vertices (2,0),(2,2) join two island segments and the wall's vertical edge
        # (JUNCTION); A's right corners, the wall-top / span corners and B's four corners
        # are 90 degree CORNERs.
        kinds = {}
        for v in perimeter.vertices:
            kinds[v.physical_kind] = kinds.get(v.physical_kind, 0) + 1
        self.assertEqual(kinds["CORNER"], 10)
        self.assertEqual(kinds["ENDPOINT"], 2)
        self.assertEqual(kinds["JUNCTION"], 2)
        self.assertNotIn("REGULAR", kinds)
        # Excluded vertices: the two junctions and the two wall-top corners sit on the wall
        # (off-plane metal within reach); the span's far corners are exactly 1 away.
        self.assertEqual(len(perimeter.excluded_vertices), 4)
        # Physical chains: A bottom / top split at the junctions (2 + 2), A right, the span's
        # three edges, B's loop broken at 4 corners -> 12 chains.
        self.assertEqual(perimeter.chains, 12)

    def test_interactions_report_the_strip_at_exactly_R(self):
        perimeter = P.extract_perimeter(self.mesh(True), CONFIG, radius=2.0)
        interactions = P.edge_interactions(perimeter, 2.0)
        # A's right edge x = 4 faces B's left edge x = 6 across the SA gap: a parallel pair at
        # separation exactly R = 2 (the transmon shield-strip case).
        parallel = [(d, c) for _, _, d, c in interactions if c >= 1.0 - 1e-8]
        self.assertIn(2.0, [round(d, 12) for d, _ in parallel])
        self.assertTrue(all(d <= 4.0 + 1e-9 for _, _, d, _ in interactions))

    def test_corner_snap_rule_matches_classifier(self):
        # The geometric joint noise rule (metaledge.hpp kJointNoiseSagittaOverRadius = 0.05,
        # USER decision 121 (B)): a joint between two unit pieces turning by t is REGULAR iff
        # (1 / 2) tan(t / 4) < 0.05 R. At R = 1 the threshold turn is 4 atan(0.1) = 22.83 deg:
        # 22.8 deg is REGULAR, 22.9 deg a CORNER, and so is a 1.001 deg joint at R = 0.0025
        # (threshold 4 atan(0.00025) = 0.0573 deg) where the former 1 deg angular threshold
        # read it; a 29.999 deg joint is a corner at R = 1 unless an arc absorbs it.
        for turn, radius, expected in ((22.8, 1.0, "REGULAR"), (22.9, 1.0, "CORNER"), (1.001, 0.0025, "CORNER"), (1.001, 1.0, "REGULAR"), (29.999, 1.0, "CORNER")):
            angle = np.radians(turn)
            nodes = [(0, 0, 0), (1, 0, 0), (1 + np.cos(angle), np.sin(angle), 0), (0, 1, 0), (1, 1, 0), (1 + np.cos(angle), 1 + np.sin(angle), 0)]
            elements = [(2, 5, (1, 2, 5)), (2, 5, (1, 5, 4)), (2, 5, (2, 3, 6)), (2, 5, (2, 6, 5))]
            path = os.path.join(self.directory.name, "turn.msh2")
            write_msh2(path, nodes, elements, [(2, 5, "metal")], True)
            perimeter = P.extract_perimeter(read_msh2(path), CONFIG, radius=radius)
            vertex = perimeter.vertices[[i for i, v in enumerate(perimeter.vertices) if np.allclose(v.point, (1, 0, 0))][0]]
            self.assertEqual(vertex.physical_kind, expected, msg=f"turn {turn} R {radius}")
            self.assertEqual(P.joint_is_noise(angle, 1.0, radius), expected == "REGULAR")

    def test_rounded_runs_follow_the_classifier(self):
        # One open chain: arm along -x -> 90 deg fillet (radius 1 = R / 2, 8 chords) -> arm
        # along +y; the mesh puts a collinear vertex in the middle of one fillet chord and
        # another one on each arm (refinement midpoints), and a second bend with one-chord
        # "arms" (a polyline arc of radius 5 = 2.5 R, no fillet) sits on the +y arm.
        R = 2.0
        radius = 1.0
        arc = [np.array([-radius + radius * math.cos(t), radius + radius * math.sin(t), 0.0]) for t in np.linspace(-math.pi / 2, 0.0, 9)]
        points = [np.array([-6.0, 0.0, 0.0]), np.array([-3.5, 0.0, 0.0])]
        for i, p in enumerate(arc):
            if i == 4:
                points.append(0.5 * (arc[3] + arc[4]))  # collinear midpoint inside the fillet
            points.append(p)
        points += [np.array([0.0, 3.0, 0.0]), np.array([0.0, 6.0, 0.0])]
        # Bend of radius 5 (> R): 45 deg in 20 chords of 0.196 with a collinear split after
        # every 5 chords (the audit must not read the sub-arcs as fillets with one-chord arms).
        big = 5.0
        for j in range(1, 21):
            t = math.radians(45.0 * j / 20)
            q = np.array([big - big * math.cos(t), 6.0 + big * math.sin(t), 0.0])
            if j % 5 == 0 and j < 20:
                points.append(0.5 * (points[-1] + q))
            points.append(q)
        points.append(points[-1] + 4.0 * np.array([math.sin(math.radians(45.0)), math.cos(math.radians(45.0)), 0.0]))  # tangent arm
        vertices = [P.PerimeterVertex(point=p, physical_kind="REGULAR") for p in points]
        vertices[0].physical_kind = vertices[-1].physical_kind = "ENDPOINT"
        edges = []
        for i in range(len(points) - 1):
            d = points[i + 1] - points[i]
            edges.append(P.PerimeterEdge(vertices=(i, i + 1), length=float(np.linalg.norm(d)), attributes=(5,), kind="PHYSICAL", conductors=("PEC",), interfaces=((0, "MA"),), inward=np.array([-d[1], d[0], 0.0]) / np.linalg.norm(d), chain=0))
            vertices[i].edges.append(i)
            vertices[i + 1].edges.append(i)
        perimeter = P.Perimeter(vertices=vertices, edges=edges, process_normal=np.array([0.0, 0.0, 1.0]), planes=[0.0], chains=1)
        runs = P.rounded_runs(perimeter, R)
        rounded = [r for r in runs if r["Rounded"]]
        self.assertEqual(len(rounded), 1, msg=str(runs))
        self.assertAlmostEqual(rounded[0]["Radius"], radius, places=9)
        self.assertAlmostEqual(rounded[0]["AngleDegrees"], 90.0, places=9)
        self.assertAlmostEqual(rounded[0]["TotalTurnDegrees"], 90.0, places=9)
        # The radius-5 bend is one arc sequence (its collinear splits merge), not a fillet.
        self.assertEqual(len(runs), 2, msg=str(runs))
        self.assertAlmostEqual([r for r in runs if not r["Rounded"]][0]["TotalTurnDegrees"], 45.0, places=9)


    def test_arc_groups_resolve_every_chord(self):
        # The every-point clause of the arc rule (USER decision 184 (2), C++ ChordsResolved;
        # the audit mirror lacked it: lane-I follow-up, decision 203 MAJOR-2). A 500 um
        # straight island edge whose ends each carry two 1.2 / 2.4 deg noise joints on 12 um
        # pieces (exactly concyclic by mirror symmetry; the E8-3 unit geometry) is NO arc: the
        # 500 um chord between two noise joints bows 7.8 um off the ~4 mm circle. The 90 deg
        # fillet of radius 1.0 um = R / 1.9 in two chords (points on the circle; its 22.5 /
        # 45 deg joints on 0.77 um chords are noise too) that follows the far noise joint is a rounded corner of its
        # own: before, the joints-only concyclicity read the straight edge and its noise
        # joints as a false 4 mm bend (DS-OSC-003 / DS-SCT-002 false clusters) and, at the
        # DS-CTX-003 loop ends, swallowed such a fillet into a false 26.6 R bend, so the audit
        # counted 8,486 rounded corners against the classifier's 8,490.
        R = 1.9
        half, piece, deg = 250.0, 12.0, math.pi / 180.0
        m0, m1 = np.array([-half, 0.0, 0.0]), np.array([half, 0.0, 0.0])
        q1 = m1 + piece * np.array([math.cos(1.2 * deg), math.sin(1.2 * deg), 0.0])
        q2 = q1 + piece * np.array([math.cos(3.6 * deg), math.sin(3.6 * deg), 0.0])
        p1 = np.array([-q1[0], q1[1], 0.0])
        p2 = np.array([-q2[0], q2[1], 0.0])
        points = [np.array([p2[0] - 20.0, p2[1] - 20.0, 0.0]), p2, p1, m0, m1, q1, q2]
        # A straight lead of 6 um along the last piece's direction, then the fillet (radius
        # 1.0, 90 deg toward +y in 2 chords) and a 6 um arm.
        t = (q2 - q1) / np.linalg.norm(q2 - q1)
        n = np.array([-t[1], t[0], 0.0])
        lead_end = q2 + 6.0 * t
        points.append(lead_end)
        centre = lead_end + 1.0 * n
        for k in (1, 2):
            theta = 0.5 * math.pi * k / 2
            points.append(centre - 1.0 * n * math.cos(theta) + 1.0 * t * math.sin(theta))
        points.append(points[-1] + 6.0 * n)
        vertices = [P.PerimeterVertex(point=p) for p in points]
        edges = []
        for i in range(len(points) - 1):
            d = points[i + 1] - points[i]
            edges.append(P.PerimeterEdge(vertices=(i, i + 1), length=float(np.linalg.norm(d)), attributes=(5,), kind="PHYSICAL", conductors=("PEC",), interfaces=((0, "MA"),), inward=np.array([-d[1], d[0], 0.0]) / np.linalg.norm(d), chain=0))
            vertices[i].edges.append(i)
            vertices[i + 1].edges.append(i)
        perimeter = P.Perimeter(vertices=vertices, edges=edges, process_normal=np.array([0.0, 0.0, 1.0]), planes=[0.0], chains=1)
        P.classify_vertices(perimeter, R)
        kinds = [v.physical_kind for v in vertices]
        # The 1.2 / 2.4 deg joints on 12 um pieces and the fillet's 22.5 / 45 deg joints on
        # 0.77 um chords are all noise under the geometric joint rule (implied sagitta below
        # 0.05 R): one chain, regular vertices; the arc rule tells them apart.
        self.assertEqual(kinds[2:10], ["REGULAR"] * 8, msg=str(kinds))
        arcs = P.arc_groups(perimeter, R)
        bends = [a for a in arcs if not a["Rounded"]]
        rounded = [a for a in arcs if a["Rounded"]]
        self.assertEqual(bends, [], msg=str(arcs))
        self.assertEqual(len(rounded), 1, msg=str(arcs))
        self.assertAlmostEqual(rounded[0]["Radius"], 1.0, places=6)
        self.assertAlmostEqual(rounded[0]["TurnDegrees"], 90.0, places=6)


    def test_arc_groups_collinear_subdivision_and_start_rule(self):
        # Decisions 212 / 213 (VALIDATION-PLAN (h)-8), the audit mirror of the C++ arc test:
        # (A) the end-joint test of a least-squares bend reads the arm's far vertex off the
        # arm's straight PIECE (the rigid-run joint), not the adjacent mesh edge, so inserting
        # collinear vertices on a chord (a thin mesh at LC 4 um) does not change the reading;
        # (B') a chord arm at the range's first joint whose far joint is itself absorbable is no
        # arm, so a closed-loop scan starting inside an arc (its longest piece a chord of the
        # arc) finds the whole arc from its real start instead of chopping it there. Geometry
        # (the unit section "a closed loop whose longest piece is a chord inside the arc"): an
        # 8 um bar whose 6 um leads meet a 120 um bend of ten unequal chords (6-10 deg, the
        # longest mid-arc) at a 20 deg kink; one 9-joint bend per side (radii 116 / 124), the
        # kink joints corners, identical joint-only, subdivided at 4 um and for every start
        # vertex of the loop.
        R = 1.9
        deg = math.pi / 180.0

        def kinked_bar(width, radius, chord_degrees, lead, kink):
            angles = [0.0]
            for c in chord_degrees:
                angles.append(angles[-1] + c * deg)
            a_in, a_out = -kink * deg, angles[-1] + kink * deg
            centre_line = [np.array([-lead * math.cos(a_in), -lead * math.sin(a_in)])]
            normals = [np.array([-math.sin(a_in), math.cos(a_in)])]
            for a in angles:
                centre_line.append(np.array([radius * math.sin(a), radius - radius * math.cos(a)]))
                normals.append(np.array([-math.sin(a), math.cos(a)]))
            end = centre_line[-1]
            centre_line.append(end + lead * np.array([math.cos(a_out), math.sin(a_out)]))
            normals.append(np.array([-math.sin(a_out), math.cos(a_out)]))
            h = 0.5 * width
            right = [p - h * n for p, n in zip(centre_line, normals)]
            left = [p + h * n for p, n in zip(centre_line, normals)]
            return right + left[::-1]  # counter-clockwise

        def reading(points):
            return arc_reading(points, R)

        design = kinked_bar(8.0, 120.0, [6.0, 7.0, 8.0, 7.0, 6.0, 10.0, 6.0, 7.0, 8.0, 7.0], 6.0, 20.0)
        plain = reading(design)
        self.assertEqual([(r, j, rounded) for r, _, j, rounded in plain], [(116.0, 9, False), (124.0, 9, False)], msg=str(plain))
        self.assertEqual(reading(subdivide_closed_polygon(design, 4.0)), plain)
        for start in range(1, len(design)):
            rotated = design[start:] + design[:start]
            self.assertEqual(reading(rotated), plain, msg="start vertex %d" % start)
            self.assertEqual(reading(subdivide_closed_polygon(rotated, 4.0)), plain, msg="start vertex %d subdivided" % start)

    def test_arc_groups_coarse_pad_collinear_subdivision(self):
        # Decision 212 fix A alone (the C4 class of the unit section "coarse polygon arcs read
        # alike joint-only and with subdivided chords"): a round pad of `chords` equal chords
        # with a 5 um lead attached at the top; the attach joints turn above the 50 deg cap, so
        # the bend is the least-squares fit over the interior joints and its end arms are the
        # first and last chords of the circle. Joint-only the arm's far vertex (the attach
        # joint) is on the circle; subdivided at 4 um the mirror read the far end of the mesh
        # edge instead, 0.8 um off the circle, and the pads lost their bends (main: [] for the
        # r 23 / 36 / 80 pads). The reading is the <= 180 deg cutting of the current rule:
        # r 23 um x 8 chords -> one 4-joint bend (DS-CTX-003 C4), the ground's r 36 x 9 -> two
        # 4-joint bends; the finely chorded r 95.5 x 60 pad (tangent bends) is the control.
        R = 1.9

        def round_pad_with_lead(rho, chords, w, L):
            theta0 = math.atan2(math.sqrt(rho * rho - w * w), w)
            y0 = math.sqrt(rho * rho - w * w)
            a0, a1 = math.pi - theta0, 2.0 * math.pi + theta0
            points = [np.array([w, y0 + L]), np.array([-w, y0 + L])]
            for k in range(chords + 1):
                angle = a0 + (a1 - a0) * k / chords
                points.append(np.array([rho * math.cos(angle), rho * math.sin(angle)]))
            return points  # counter-clockwise

        for rho, chords, bends in ((23.0, 8, [(23.0, 4)]), (36.0, 9, [(36.0, 4), (36.0, 4)]), (80.0, 8, [(80.0, 4)]), (95.5, 60, [(95.5, 29), (95.5, 30)])):
            design = round_pad_with_lead(rho, chords, 2.5, 60.0)
            plain = arc_reading(design, R)
            self.assertEqual([(r, j) for r, _, j, rounded in plain if not rounded], bends, msg="pad %g x %d: %s" % (rho, chords, plain))
            self.assertEqual(arc_reading(subdivide_closed_polygon(design, 4.0), R), plain, msg="pad %g x %d subdivided" % (rho, chords))
            self.assertEqual(arc_reading(subdivide_closed_polygon(design, 1.0), R), plain, msg="pad %g x %d subdivided at 1 um" % (rho, chords))


class ManifestTest(unittest.TestCase):
    def test_digest_ignores_status_and_library(self):
        a = make_manifest([requirement("IsolatedEdge", 3, 6.0, {"EdgeCount": 1}, "Exact", SelectedModels=[{"Name": "iso", "Topology": "IsolatedEdge", "Weight": 1.0}], NormalizedLibraryDistance=0.0)])
        b = make_manifest([requirement("IsolatedEdge", 3, 6.0, {"EdgeCount": 1}, "Missing", Reason="No model")])
        b["Library"]["Path"] = "/elsewhere/lib.json"
        b["Statistics"] = {"Geometry": {"MetalSegments": 99}}
        self.assertEqual(M.canonical_digest(a), M.canonical_digest(b))
        self.assertTrue(M.diff_manifests(a, b)["Identical"])

    def test_digest_and_diff_see_geometry_changes(self):
        a = make_manifest([requirement("ConvexCorner", 4, 16.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        b = make_manifest([requirement("ConvexCorner", 4, 16.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0}), requirement("SameConductorStrip", 2, 4.0, {"EdgeCount": 2, "Separation": 2.0})])
        c = make_manifest([requirement("ConvexCorner", 3, 12.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        self.assertNotEqual(M.canonical_digest(a), M.canonical_digest(b))
        diff = M.diff_manifests(a, b)
        self.assertEqual(len(diff["Added"]), 1)
        self.assertEqual(diff["Added"][0]["Topology"], "SameConductorStrip")
        self.assertEqual(len(M.diff_manifests(a, c)["Changed"]), 1)

    def test_geometry_only_digest_is_refinement_invariant_for_translational(self):
        a = make_manifest([requirement("IsolatedEdge", 3, 6.0, {"EdgeCount": 1}), requirement("ConvexCorner", 4, 16.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        b = make_manifest([requirement("IsolatedEdge", 6, 6.0, {"EdgeCount": 1}), requirement("ConvexCorner", 4, 16.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        self.assertNotEqual(M.canonical_digest(a), M.canonical_digest(b))
        self.assertEqual(M.canonical_digest(a, geometry_only_counts=True), M.canonical_digest(b, geometry_only_counts=True))
        c = make_manifest([requirement("IsolatedEdge", 6, 6.0, {"EdgeCount": 1}), requirement("ConvexCorner", 5, 20.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        self.assertNotEqual(M.canonical_digest(a, geometry_only_counts=True), M.canonical_digest(c, geometry_only_counts=True))

    def test_weight_defects(self):
        m = make_manifest([requirement("IsolatedEdge", 1, 2.0, {"EdgeCount": 1}, SelectedModels=[{"Name": "a", "Topology": "IsolatedEdge", "Weight": 0.6}, {"Name": "b", "Topology": "IsolatedEdge", "Weight": 0.3}])])
        self.assertEqual(len(M.summarize(m)["WeightDefects"]), 1)

    def test_log_parser(self):
        with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as log:
            log.write("Running with 4 MPI processes\nAdded 12 elements in 1 iterations of local bisection for under-resolved interior boundaries\n")
            log.write("\x1b[33mOmitting 2 of 2 three-dimensional target edge segments which are within 2R of a physical metal edge with a different interface mapping.\x1b[0m\n")
            log.write("Omitting 67 of 3074 three-dimensional target edge segments in unsupported local interaction neighborhoods (nonparallel: 2, incompatible process normal: 0, process-normal offset: 0, unclassified topology: 0, missing library model: 64, multi-edge: 0). Unclassified pairs by conductor ownership: same = 0, different = 0.\n")
            log.write(" Matched physical edge segments: 3007\n Matched corner patches: 34\nThe selected three-dimensional metal perimeter has 6 unmatched corner, endpoint, or junction vertices.\n")
            name = log.name
        parsed = M.parse_palace_log(name)
        os.unlink(name)
        self.assertEqual(parsed["Ranks"], 4)
        self.assertEqual(parsed["Bisection"], {"AddedElements": 12, "Iterations": 1})
        self.assertEqual(parsed["OmittedSegments"], 2 + 67 - 64)
        self.assertEqual(parsed["OmittedUnsupported"][0]["Nonparallel"], 2)
        self.assertEqual(parsed["MatchedSegments"], 3007)
        self.assertEqual(parsed["UnmatchedVertices"], 6)


class AuditGateTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        nodes, elements, names = island_mesh()
        self.mesh_path = os.path.join(self.directory.name, "island.msh2")
        write_msh2(self.mesh_path, nodes, elements, names, True)
        self.config_path = os.path.join(self.directory.name, "config.json")
        with open(self.config_path, "w") as target:
            json.dump(CONFIG, target)

    def tearDown(self):
        self.directory.cleanup()

    def run_gates(self, manifest, log=None):
        manifest_path = os.path.join(self.directory.name, "req.json")
        with open(manifest_path, "w") as target:
            json.dump(manifest, target)
        argv = ["--mesh", self.mesh_path, "--config", self.config_path, "--manifest", manifest_path, "--output-prefix", os.path.join(self.directory.name, "out", "audit")]
        if log:
            log_path = os.path.join(self.directory.name, "palace.log")
            with open(log_path, "w") as target:
                target.write(log)
            argv += ["--log", log_path]
        with contextlib.redirect_stdout(io.StringIO()):
            code = audit.main(argv)
        with open(os.path.join(self.directory.name, "out", "audit.json")) as source:
            result = json.load(source)
        result["ByGate"] = {g["Gate"]: g["Detail"] for g in result["Gates"]}
        return code, {g["Gate"]: g["Status"] for g in result["Gates"]}, result

    def test_complete_partition_passes_except_the_silent_exclusions(self):
        # Version-1 (aggregate) manifest: 22 of physical length in 12 segments of which 6 lie
        # in cross-layer zones (16 targeted), 8 feature corners; the 2 endpoints of A's open
        # chains lie on the truncation edge and are simulation cuts, not features.
        manifest = make_manifest(
            [
                requirement("IsolatedEdge", 10, 12.0, {"EdgeCount": 1}),
                requirement("SameConductorStrip", 2, 4.0, {"EdgeCount": 2, "Separation": 2.0}),
                requirement("ConvexCorner", 8, 32.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0}),
            ],
            metal_segments=13,
            chains=12,
        )
        code, gates, result = self.run_gates(manifest)
        self.assertEqual(gates["perimeter-agreement"], "PASS")
        self.assertEqual(gates["chain-agreement"], "PASS")
        self.assertEqual(gates["A1-length-partition"], "PASS")
        self.assertEqual(gates["A1-count-partition"], "PASS")
        self.assertEqual(gates["A1-multiplicity"], "PASS")
        self.assertEqual(gates["A1-vertex-census"], "PASS")
        self.assertEqual(result["ByGate"]["A1-vertex-census"]["AuditTruncationCuts"], 2)
        self.assertEqual(gates["A2-cluster-balls"], "PASS")
        # The wall / span / zones are never reported by a version-1 manifest: the exclusion
        # gate fails (excluded length: nonplanar 2 + fold 2 + nonmanifold 2 + zones 6).
        self.assertEqual(gates["A1-exclusions-recorded"], "FAIL")
        self.assertEqual(code, 1)
        self.assertAlmostEqual(result["GapBound"]["ExcludedLength"], 12.0)

    def test_gap_and_double_count_fail(self):
        manifest = make_manifest([requirement("IsolatedEdge", 11, 14.0, {"EdgeCount": 1}), requirement("ConvexCorner", 8, 32.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        code, gates, result = self.run_gates(manifest, log="Omitting 1 of 12 three-dimensional target edge segments which are within 2R of a physical metal edge with a different interface mapping.\n")
        self.assertEqual(gates["A1-length-partition"], "FAIL")
        self.assertEqual(gates["A1-count-partition"], "FAIL")
        self.assertEqual(gates["A1-multiplicity"], "NOT-EVALUABLE")
        self.assertAlmostEqual(result["GapBound"]["OmittedLength"], 2.0)
        self.assertTrue(result["ByGate"]["A1-count-partition"]["Reconciled"])
        manifest = make_manifest([requirement("IsolatedEdge", 10, 20.0, {"EdgeCount": 1}), requirement("ConvexCorner", 8, 32.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        code, gates, _ = self.run_gates(manifest)
        self.assertEqual(gates["A1-count-partition"], "FAIL")
        self.assertEqual(code, 1)

    def test_vertex_census_residual_and_compare(self):
        manifest = make_manifest([requirement("IsolatedEdge", 12, 16.0, {"EdgeCount": 1}), requirement("ConvexCorner", 3, 12.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        code, gates, result = self.run_gates(manifest)
        # 8 audit feature corners, 3 in the manifest, no clusters: a silent drop of 5 -> FAIL.
        self.assertEqual(gates["A1-vertex-census"], "FAIL")
        self.assertEqual(result["ByGate"]["A1-vertex-census"]["Residual"], 5)
        other = os.path.join(self.directory.name, "other.json")
        with open(other, "w") as target:
            json.dump(make_manifest([requirement("IsolatedEdge", 12, 16.0, {"EdgeCount": 1}), requirement("ConvexCorner", 3, 12.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0}, "Missing", Reason="x")]), target)
        manifest_path = os.path.join(self.directory.name, "req.json")
        with contextlib.redirect_stdout(io.StringIO()):
            code = audit.main(["--mesh", self.mesh_path, "--config", self.config_path, "--manifest", manifest_path, "--compare", other, "--output-prefix", os.path.join(self.directory.name, "out", "cmp")])
        with open(os.path.join(self.directory.name, "out", "cmp.json")) as source:
            result = json.load(source)
        self.assertTrue(result["Compare"]["Diff"]["Identical"])
        self.assertEqual({g["Gate"]: g["Status"] for g in result["Gates"]}["A3/A5-set-identity"], "PASS")
        # With a cluster record present the residual is not evaluable (absorbed or dropped).
        clustered = make_manifest([requirement("IsolatedEdge", 12, 16.0, {"EdgeCount": 1}), requirement("ConvexCorner", 3, 12.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0}), requirement("SpatialEdgeCluster", 1, 2.0, {"EdgeCount": 3, "Edges": []})])
        _, gates_clustered, _ = self.run_gates(clustered)
        self.assertEqual(gates_clustered["A1-vertex-census"], "NOT-EVALUABLE")


class IdentificationGateTest(AuditGateTest):
    """Version-2 manifests: the gates read the per-segment assignment, vertex and exclusion
    tables (palace/models/SURFACE-RESPONSE-IDENTIFICATION.md)."""

    EXCLUSION_CLASS = {"TRUNCATION": "TruncationCut", "NONPLANAR": "NonPlanar", "FOLD": "NonPlanar", "NONMANIFOLD": "NonManifold", "EMBEDDED": "UndeterminedProcessSide", "BOX": "SimulationBoundary"}

    def identification(self):
        """A complete identification of the island mesh built from the audit's own perimeter:
        one IsolatedEdge feature per chain claiming its segments outside the cross-layer zones
        (recorded as ExcludedPortions of the CrossLayer record), eight corners, four excluded
        vertices, the wall / span and truncation edges recorded as exclusions."""
        from .msh2 import read_msh2

        perimeter = P.extract_perimeter(read_msh2(self.mesh_path), CONFIG, radius=0.5)
        features, segments, vertices, exclusions = [], [], [], {}
        order = []
        chain_feature = {}

        def record(cls, count, length):
            if cls not in exclusions:
                exclusions[cls] = [0, 0.0]
                order.append(cls)
            exclusions[cls][0] += count
            exclusions[cls][1] += length

        for edge in perimeter.edges:
            p0, p1 = perimeter.edge_points(edge)
            key = [list(map(float, p0)), list(map(float, p1))]
            if key[1] < key[0]:
                key.reverse()
            entry = {"Key": key, "Length": edge.length, "Chain": edge.chain}
            if edge.kind != "PHYSICAL":
                cls = self.EXCLUSION_CLASS[edge.kind]
                entry["Exclusion"] = {"Class": cls, "Reason": "test"}
                record(cls, 1, edge.length)
            else:
                if edge.chain not in chain_feature:
                    chain_feature[edge.chain] = len(features)
                    features.append({"Id": len(features), "Type": "IsolatedEdge", "Signature": {"Type": "IsolatedEdge", "Interfaces": ["MS", "SA"], "Law": "{\"Type\":\"PEC\"}"}, "Hash": f"h{edge.chain}", "Chirality": 1, "Length": 0.0, "Portions": [], "Vertices": [], "Frame": {"Origin": [0, 0, 0], "Axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}, "Match": {"Status": "Missing"}})
                f = features[chain_feature[edge.chain]]
                # Portions along the canonical key direction; zones are in edge order.
                forward = list(map(float, p0)) == key[0]
                zones = sorted(z if forward else (edge.length - z[1], edge.length - z[0]) for z in edge.cross_layer)
                cursor = 0.0
                portions = []
                for a, b in zones:
                    if a > cursor + 1e-12:
                        portions.append([cursor, a, f["Id"]])
                    cursor = b
                if cursor < edge.length - 1e-12:
                    portions.append([cursor, edge.length, f["Id"]])
                for a, b, _ in portions:
                    f["Portions"].append([len(segments), a, b])
                    f["Length"] += b - a
                entry["Portions"] = portions
                if zones:
                    record("CrossLayer", len(zones), sum(b - a for a, b in zones))
                    entry["ExcludedPortions"] = [[a, b, order.index("CrossLayer")] for a, b in zones]
            segments.append(entry)
        for index, vertex in enumerate(perimeter.vertices):
            if vertex.physical_kind in (None, "REGULAR"):
                continue
            if vertex.physical_kind == "ENDPOINT":
                vertices.append({"Vertex": index, "Type": "TruncationCut"})
            elif index in perimeter.excluded_vertices:
                vertices.append({"Vertex": index, "Type": "Excluded"})
            else:
                vertices.append({"Vertex": index, "Type": "ConvexCorner", "TurnDegrees": 90.0, "Feature": 0})
        assigned = sum(f["Length"] for f in features)
        excluded = sum(v[1] for v in exclusions.values())
        return {
            "Version": 2,
            "MatchingRadius": 0.5,
            "Conventions": {"JointNoiseSagittaOverR": 0.05, "ArcMaxJointTurnDegrees": 50.0},
            "ReferenceProcessNormal": [0.0, 0.0, 1.0],
            "Features": features,
            "Segments": segments,
            "Vertices": vertices,
            "Exclusions": [{"Class": k, "Reason": "test", "Count": exclusions[k][0], "Length": exclusions[k][1]} for k in order],
            "Totals": {"PerimeterLength": assigned + excluded, "AssignedLength": assigned, "ExcludedLength": excluded},
            # The classifier always writes Diagnostics (gate A9 reads SamePriorityClaimOverlaps).
            "Diagnostics": {"SamePriorityClaimOverlaps": 0, "StackGeometricOffsetIntervals": 0, "StackCompositionCapHits": 0, "StackCompositionCap": 64,
                            "ClusterExtension": {"Passes": 1, "AbsorbedPortions": 0, "AbsorbedLength": 0.0, "VertexFeaturesJoined": 0}, "StackEndThirdBodyLength": 0.0},
            "GeometryDigest": "digest",
        }

    def v2_manifest(self, identification):
        manifest = make_manifest([requirement("IsolatedEdge", 12, 16.0, {}), requirement("ConvexCorner", 8, 32.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        manifest["Version"] = 2
        manifest["Identification"] = identification
        return manifest

    def test_complete_identification_passes(self):
        code, gates, result = self.run_gates(self.v2_manifest(self.identification()))
        for name in ("perimeter-agreement", "A1-length-partition", "A1-count-partition", "A1-multiplicity", "A1-weights", "A1-vertex-census", "A1-exclusions-recorded", "A2-cluster-balls", "A9-same-priority-claim-overlaps"):
            self.assertEqual(gates[name], "PASS", name)
        self.assertEqual(code, 0)
        self.assertEqual(result["Identification"]["Features"]["IsolatedEdge"], 12)
        self.assertAlmostEqual(result["Identification"]["Totals"]["AssignedLength"], 16.0)
        self.assertEqual(result["ByGate"]["A1-vertex-census"]["ManifestExcludedVertices"], 4)

    def test_a9_needs_the_diagnostics(self):
        """A manifest without Diagnostics cannot be gated on the claim overlaps (NOT-EVALUABLE
        counts as a failure), a manifest reporting an overlap fails A9 (review fix-3 M1)."""
        ident = self.identification()
        del ident["Diagnostics"]
        code, gates, _ = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A9-same-priority-claim-overlaps"], "NOT-EVALUABLE")
        self.assertEqual(code, 1)
        ident = self.identification()
        ident["Diagnostics"]["SamePriorityClaimOverlaps"] = 2
        code, gates, _ = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A9-same-priority-claim-overlaps"], "FAIL")
        self.assertEqual(code, 1)

    def test_defects_fail_the_exact_gates(self):
        ident = self.identification()
        first = next(s for s in ident["Segments"] if "Portions" in s)
        first["Portions"][0][1] *= 0.5  # a gap on one segment
        _, gates, result = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A1-multiplicity"], "FAIL")
        self.assertEqual(result["ByGate"]["A1-multiplicity"]["Defects"], 1)
        ident = self.identification()
        ident["Exclusions"] = [e for e in ident["Exclusions"] if e["Class"] != "NonPlanar"]
        _, gates, _ = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A1-exclusions-recorded"], "FAIL")
        ident = self.identification()
        next(v for v in ident["Vertices"] if v["Type"] == "ConvexCorner")["Feature"] = -1
        _, gates, _ = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A1-vertex-census"], "FAIL")
        ident = self.identification()
        ident["Vertices"].remove(next(v for v in ident["Vertices"] if v["Type"] == "ConvexCorner"))
        _, gates, _ = self.run_gates(self.v2_manifest(ident))
        self.assertEqual(gates["A1-vertex-census"], "FAIL")

    def test_compare_uses_the_geometry_digest(self):
        ident = self.identification()
        other_path = os.path.join(self.directory.name, "other.json")
        other = self.v2_manifest(self.identification())
        other["Identification"]["GeometryDigest"] = "different"
        with open(other_path, "w") as target:
            json.dump(other, target)
        manifest_path = os.path.join(self.directory.name, "req.json")
        with open(manifest_path, "w") as target:
            json.dump(self.v2_manifest(ident), target)
        with contextlib.redirect_stdout(io.StringIO()):
            audit.main(["--mesh", self.mesh_path, "--config", self.config_path, "--manifest", manifest_path, "--compare", other_path, "--output-prefix", os.path.join(self.directory.name, "out", "cmp")])
        with open(os.path.join(self.directory.name, "out", "cmp.json")) as source:
            result = json.load(source)
        self.assertEqual({g["Gate"]: g["Status"] for g in result["Gates"]}["A3/A5-set-identity"], "FAIL")
        self.assertFalse(result["Compare"]["GeometryDigestIdentical"])


class PatchGateTest(unittest.TestCase):
    """Gates A7 on the patch dry run (phase 4): the patched set is the matched set, every
    portion of a matched longitudinal feature is one quadrature interval with weights summing
    to one, vertex / cluster features carry one patch, nothing on unmatched / excluded."""

    R = 2.0

    def identification(self):
        # Segments 0-1: an isolated edge (feature 0, matched); 2-3: a strip pair on two
        # chains (feature 1, matched); 4: a corner window + the corner (feature 2, matched);
        # 5: unmatched isolated edge (feature 3); 6: excluded (truncation).
        features = [
            {"Id": 0, "Type": "IsolatedEdge", "Length": 5.0, "Portions": [[0, 0.0, 3.0], [1, 0.0, 2.0]], "Vertices": [], "Match": {"Status": "Matched", "Model": "iso"}},
            # Both sides of the strip lie on one chain (a slot): the sides come from Sides.
            {"Id": 1, "Type": "SameConductorStrip", "Length": 4.0, "Portions": [[2, 0.0, 2.0], [3, 0.0, 2.0]], "Sides": [0, 1], "Vertices": [], "Match": {"Status": "Matched", "Model": "strip"}},
            {"Id": 2, "Type": "ConvexCorner", "Length": 2.0, "Portions": [[4, 0.0, 2.0]], "Vertices": [7], "Match": {"Status": "Matched", "Model": "corner"}},
            {"Id": 3, "Type": "IsolatedEdge", "Length": 1.0, "Portions": [[5, 0.0, 1.0]], "Vertices": [], "Match": {"Status": "Missing"}},
        ]
        segments = [
            {"Key": [[0, 0, 0], [3, 0, 0]], "Length": 3.0, "Chain": 0, "Portions": [[0.0, 3.0, 0]]},
            {"Key": [[3, 0, 0], [5, 0, 0]], "Length": 2.0, "Chain": 0, "Portions": [[0.0, 2.0, 0]]},
            {"Key": [[0, 1, 0], [2, 1, 0]], "Length": 2.0, "Chain": 1, "Portions": [[0.0, 2.0, 1]]},
            {"Key": [[0, 2, 0], [2, 2, 0]], "Length": 2.0, "Chain": 1, "Portions": [[0.0, 2.0, 1]]},
            {"Key": [[0, 3, 0], [2, 3, 0]], "Length": 2.0, "Chain": 3, "Portions": [[0.0, 2.0, 2]]},
            {"Key": [[0, 4, 0], [1, 4, 0]], "Length": 1.0, "Chain": 4, "Portions": [[0.0, 1.0, 3]]},
            {"Key": [[0, 5, 0], [1, 5, 0]], "Length": 1.0, "Chain": 5, "Exclusion": {"Class": "TruncationCut", "Reason": "test"}},
        ]
        return {"Version": 2, "MatchingRadius": self.R, "Features": features, "Segments": segments, "Vertices": [], "Exclusions": [{"Class": "TruncationCut", "Reason": "test", "Count": 1, "Length": 1.0}], "Totals": {"PerimeterLength": 13.0, "AssignedLength": 12.0, "ExcludedLength": 1.0}, "GeometryDigest": "d"}

    def patches(self):
        rows = []

        def longitudinal(feature, topology, segment, s0, s1, side, depth):
            for q, weight in ((0.2113248654051871, 0.5), (0.7886751345948129, 0.5)):
                rows.append({"Patch": len(rows), "Feature": feature, "Topology": topology, "Model": topology, "ModelIndex": 1, "Weight": weight * (s1 - s0) * side / depth, "ModelWeight": 1.0, "QuadratureWeight": weight, "SideFactor": side, "CouponDepth": depth, "Segment": segment, "S0": s0, "S1": s1})

        longitudinal(0, "isolated edge", 0, 0.0, 3.0, 1.0, 2.0)
        longitudinal(0, "isolated edge", 1, 0.0, 2.0, 1.0, 2.0)
        longitudinal(1, "same-conductor strip", 2, 0.0, 2.0, 0.5, 2.0)
        longitudinal(1, "same-conductor strip", 3, 0.0, 2.0, 0.5, 2.0)
        rows.append({"Patch": len(rows), "Feature": 2, "Topology": "convex corner", "Model": "corner", "ModelIndex": 2, "Weight": 1.0, "ModelWeight": 1.0, "QuadratureWeight": 1.0, "SideFactor": 1.0, "CouponDepth": 0.0, "Segment": -1, "S0": 0.0, "S1": 0.0})
        return rows

    def gates(self, ident, rows):
        gates, summary = audit.patch_gates(ident, rows, self.R)
        return {g["Gate"]: g["Status"] for g in gates}, {g["Gate"]: g["Detail"] for g in gates}, summary

    def test_complete_dry_run_passes(self):
        status, detail, summary = self.gates(self.identification(), self.patches())
        self.assertEqual(status, {"A7-patch-features": "PASS", "A7-patch-exclusions": "PASS", "A7-patch-coverage": "PASS", "A7-patch-weights": "PASS"})
        self.assertAlmostEqual(summary["CoveredLength"], 11.0)  # the unmatched edge (1.0) is not covered
        self.assertAlmostEqual(summary["CoveredFractionOfAssigned"], 11.0 / 12.0)

    def test_defects_fail(self):
        ident = self.identification()
        rows = self.patches()
        status, detail, _ = self.gates(ident, rows[:-1])  # the corner without a patch
        self.assertEqual(status["A7-patch-features"], "FAIL")
        self.assertEqual(detail["A7-patch-features"]["MatchedNotPatched"], [2])
        rows = self.patches()
        rows[0]["Feature"] = 3  # a patch on the unmatched feature
        status, _, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-exclusions"], "FAIL")
        rows = self.patches()
        for r in rows[:2]:
            r["Segment"] = 6  # a patch on the excluded segment
        status, _, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-exclusions"], "FAIL")
        rows = self.patches()
        for r in rows[:2]:
            r["S1"] = 2.0  # a shorter interval: the portion is not covered exactly
        status, _, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-coverage"], "FAIL")
        rows = self.patches()
        rows[0]["QuadratureWeight"] = 0.4  # quadrature weights no longer sum to one
        status, _, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        rows = self.patches()
        for r in rows[4:8]:
            r["SideFactor"] = 1.0  # a pair integrated once per side would double count
        status, _, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")

    def test_continuation_ownership_scales_the_quadrature_weights(self):
        # Decision 236 (2): a cell owned by a spatial coupon keeps its portion [S0, S1) with
        # its quadrature weight scaled by the kept fraction; the manifest's
        # ContinuationOwnership.OwnedCells reconciles the portion's sum (1 - owned / length).
        ident = self.identification()
        rows = self.patches()
        # Feature 0, segment 0 (portion [0, 3], two half cells of 1.5): the first cell wholly
        # owned (weight 0), the second clipped to its outer 0.9 (kept fraction 0.6).
        rows[0]["QuadratureWeight"] = 0.0
        rows[0]["Weight"] = 0.0
        rows[1]["QuadratureWeight"] = 0.5 * 0.6
        rows[1]["Weight"] = 0.5 * 0.6 * 3.0 / 2.0
        ident["Diagnostics"] = {"ContinuationOwnership": {"OwnedCells": [
            {"Patch": 0, "Feature": 0, "Segment": 0, "S0": 0.0, "S1": 3.0, "CellLength": 1.5, "OwnedLength": 1.5},
            {"Patch": 1, "Feature": 0, "Segment": 0, "S0": 0.0, "S1": 3.0, "CellLength": 1.5, "OwnedLength": 0.6},
        ]}}
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "PASS", detail["A7-patch-weights"])
        self.assertEqual(status["A7-patch-coverage"], "PASS")
        # Without the record the scaled weights are a defect (fail-before).
        del ident["Diagnostics"]
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        self.assertEqual(detail["A7-patch-weights"]["Examples"][0]["Defect"], "quadrature x model weights do not sum to 1 - owned / portion")

    def test_corner_arm_extension_scales_the_quadrature_weights(self):
        # Decision 511 O2 (i): the arm cell of a matched rounded corner beginning at the claim
        # end is extended back to the square exit; its quadrature weight scales by new / old
        # and the manifest's CornerArmExtension.Cells reconciles the portion's sum
        # (1 + extended / length), the mirror image of the F1 trim's record.
        ident = self.identification()
        rows = self.patches()
        # Feature 0, segment 1 (portion [0, 2], two half cells of 1.0): the first cell gains
        # 0.5 (new / old = 1.5).
        rows[2]["QuadratureWeight"] = 0.5 * 1.5
        rows[2]["Weight"] = 0.5 * 1.5 * 2.0 / 2.0
        ident["Diagnostics"] = {"CornerArmExtension": {"Cells": [
            {"Patch": 2, "Corner": 2, "Feature": 0, "Segment": 1, "ExtendedLength": 0.5},
        ]}}
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "PASS", detail["A7-patch-weights"])
        self.assertEqual(status["A7-patch-coverage"], "PASS")
        # Without the record the grown weights are a defect (fail-before).
        del ident["Diagnostics"]
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        self.assertEqual(detail["A7-patch-weights"]["Examples"][0]["Defect"], "quadrature x model weights do not sum to 1 - owned / portion")

    def test_domain_boundary_exclusion_keeps_the_portion_tiled(self):
        # Decision 258: a patch whose placed coupon section leaves the device mesh is written
        # with Weight 0 and its unscaled quadrature weight (the portion stays tiled), listed
        # in Diagnostics.DomainBoundaryExclusions with its cell length.
        ident = self.identification()
        rows = self.patches()
        rows[0]["Weight"] = 0.0
        ident["Diagnostics"] = {"DomainBoundaryExclusions": {"Count": 1, "Patches": [
            {"Patch": 0, "Feature": 0, "Segment": 0, "S0": 0.0, "S1": 3.0, "CellLength": 1.5, "PortionLength": 3.0, "OutsidePoints": 3},
        ]}}
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "PASS", detail["A7-patch-weights"])
        self.assertEqual(status["A7-patch-coverage"], "PASS")
        self.assertEqual(detail["A7-patch-weights"]["DomainBoundaryExcludedPatches"], 1)
        self.assertAlmostEqual(detail["A7-patch-weights"]["DomainBoundaryExcludedLength"], 1.5)
        # The record's cell length must be the patch's quadrature x portion.
        ident["Diagnostics"]["DomainBoundaryExclusions"]["Patches"][0]["CellLength"] = 3.0
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        self.assertEqual(detail["A7-patch-weights"]["Examples"][0]["Defect"], "domain-boundary record cell length is not the patch cell")
        # Without the record a zero weight is the weight-formula defect (fail-before); with
        # the record a nonzero weight on an excluded patch is a defect.
        del ident["Diagnostics"]
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        self.assertEqual(detail["A7-patch-weights"]["Examples"][0]["Defect"], "weight formula")
        rows[0]["Weight"] = 0.5 * 3.0 / 2.0
        ident["Diagnostics"] = {"DomainBoundaryExclusions": {"Count": 1, "Patches": [
            {"Patch": 0, "Feature": 0, "Segment": 0, "S0": 0.0, "S1": 3.0, "CellLength": 1.5, "PortionLength": 3.0, "OutsidePoints": 3},
        ]}}
        status, detail, _ = self.gates(ident, rows)
        self.assertEqual(status["A7-patch-weights"], "FAIL")
        self.assertEqual(detail["A7-patch-weights"]["Examples"][0]["Defect"], "domain-boundary excluded patch with a nonzero weight")


class SignatureLibraryTest(unittest.TestCase):
    LAW = "{\"Type\":\"PEC\"}"

    def stack(self, feature_id, offsets, gaps=(1, -1, 1, -1), conductors=(1, 2, 2, 1)):
        edges = [{"OffsetOverR": o, "GapSide": g, "Conductor": c, "Interfaces": ["SA"], "Law": self.LAW}
                 for o, g, c in zip(offsets, gaps, conductors)]
        return {"Id": feature_id, "Type": "ParallelEdgeCluster", "Hash": f"{feature_id:064x}",
                "Signature": {"Type": "ParallelEdgeCluster", "Edges": edges}, "Length": 1.0}

    def test_one_model_per_coupon_with_version1_parameters(self):
        from .signature_library import build_signature_library

        manifest = {
            "Identification": {
                "MatchingRadius": 2.0,
                "Features": [
                    {"Id": 0, "Type": "IsolatedEdge", "Hash": "a" * 64, "Signature": {"Type": "IsolatedEdge", "Interfaces": ["SA"], "Law": self.LAW}},
                    {"Id": 1, "Type": "IsolatedEdge", "Hash": "a" * 64, "Signature": {"Type": "IsolatedEdge", "Interfaces": ["SA"], "Law": self.LAW}},
                    {"Id": 2, "Type": "ConcaveCorner", "Hash": "b" * 64, "Signature": {"Type": "ConcaveCorner", "AngleDegrees": 90.0, "CornerRadiusOverR": 0.25, "Interfaces": ["SA"], "Law": self.LAW}},
                    {"Id": 3, "Type": "DifferentConductorGap", "Hash": "c" * 64, "Signature": {"Type": "DifferentConductorGap", "SeparationOverR": 0.95,
                     "Edges": [{"OffsetOverR": 0.0, "GapSide": 1, "Conductor": 1, "Interfaces": ["SA"], "Law": self.LAW}, {"OffsetOverR": 0.95, "GapSide": -1, "Conductor": 2, "Interfaces": ["SA"], "Law": self.LAW}]}},
                    {"Id": 4, "Type": "SpatialEdgeCluster", "Hash": "d" * 64, "Signature": {"Type": "SpatialEdgeCluster", "EdgeCount": 3, "Portions": [{"Conductor": 1}, {"Conductor": 2}, {"Conductor": 1}]}},
                    {"Id": 5, "Type": "UnclassifiedParallelPair", "Hash": "e" * 64, "Signature": {"Type": "UnclassifiedParallelPair"}},
                ],
            }
        }
        library = build_signature_library(manifest)
        by_topology = {}
        for m in library["Models"]:
            by_topology.setdefault(m["Topology"], []).append(m)
        self.assertEqual(len(library["Models"]), 4)  # one per coupon; the unclassified pair has no model
        self.assertEqual(library["MatchingRadius"], 2.0)
        self.assertEqual(by_topology["IsolatedEdge"][0]["Instances"], 2)
        self.assertEqual(by_topology["ConcaveCorner"][0]["Angle"], 90.0)
        self.assertEqual(by_topology["ConcaveCorner"][0]["CornerRadius"], 0.5)
        self.assertEqual(by_topology["DifferentConductorGap"][0]["Separation"], 1.9)
        self.assertEqual(len(by_topology["DifferentConductorGap"][0]["ConductorReferences"]), 2)
        self.assertEqual(len(by_topology["SpatialEdgeCluster"][0]["ConductorReferences"]), 2)
        self.assertNotIn("Edges", by_topology["SpatialEdgeCluster"][0])
        self.assertEqual(by_topology["IsolatedEdge"][0]["CouponDepth"], 2.0)
        for m in library["Models"]:
            self.assertEqual(m["Signature"]["Type"], m["Topology"])
            self.assertEqual(m["ParameterSpread"], 0.0)

    def test_instances_within_the_tolerance_are_one_coupon(self):
        """Decision 85(1): three 4-edge stacks whose offsets agree within 1e-3 R (one in the mirror
        orientation) are one coupon at the midpoint representative; a fourth 2e-3 R away is a
        second coupon; the grouping does not depend on the feature order."""
        from .signature_library import build_signature_library, signature_deviation, mirror_translational

        a = self.stack(0, [0.0, 1.0004, 2.0006, 3.0008])
        b = self.stack(1, [0.0, 0.9996, 1.9998, 3.0002])
        c = mirror_translational(self.stack(2, [0.0, 1.0002, 2.0004, 3.0006])["Signature"])
        c = {"Id": 2, "Type": "ParallelEdgeCluster", "Hash": "2" * 64, "Signature": c, "Length": 1.0}
        d = self.stack(3, [0.0, 1.0030, 2.0040, 3.0050])
        for order in ([a, b, c, d], [d, c, b, a], [b, d, a, c]):
            library = build_signature_library({"Identification": {"MatchingRadius": 2.0, "Features": order}})
            self.assertEqual(len(library["Models"]), 2, [m["Signature"]["Edges"] for m in library["Models"]])
            coupon = max(library["Models"], key=lambda m: m["Instances"])
            self.assertEqual(coupon["Instances"], 3)
            offsets = [e["OffsetOverR"] for e in coupon["Signature"]["Edges"]]
            # The lead is b (the smallest serialisation); a is nearer to it in its mirror
            # orientation (0 / 1.0002 / 2.0004 / 3.0008): midpoints of the aligned ranges.
            self.assertEqual([round(o, 6) for o in offsets], [0.0, 0.9999, 2.0001, 3.0005])
            self.assertLessEqual(coupon["ParameterSpread"], 1.0)
            for member in (a, b, c):
                self.assertLessEqual(signature_deviation(coupon["Signature"], member["Signature"]), 1.0)
            self.assertGreater(signature_deviation(coupon["Signature"], d["Signature"]), 1.0)
        self.assertIsNone(signature_deviation(a["Signature"], self.stack(9, [0.0, 1.0, 2.0, 3.0], conductors=(1, 2, 3, 1))["Signature"]))


class FacingCheckTest(unittest.TestCase):
    """facing_check.py on hand-built version-2 identifications: an isolated chain at y = 0 facing
    cluster metal at y = 1 (R = 2, 2R = 4) so that every isolated sample faces; the
    `SubToleranceFeature` exemption follows the narrowed rule of decision 89."""

    R = 2.0
    TINY = 1.0e-3 * R * 0.5  # 1 nm at R = 2 um (the tolerance is 2 nm)
    LAW = json.dumps({"Type": "PEC"})

    def segment(self, x0, x1, y, chain, portions):
        return {"Key": [[x0, y, 0.0], [x1, y, 0.0]], "Length": x1 - x0, "Chain": chain, "Portions": [[a, b, f] for a, b, f in portions]}

    def feature(self, feature_id, ftype, portions):
        frame = {"Origin": [0.0, 0.0, 0.0], "Axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}
        return {"Id": feature_id, "Type": ftype, "Signature": {"Type": ftype, "Interfaces": ["SA"], "Law": self.LAW}, "Hash": f"{feature_id:064x}",
                "Chirality": 1, "Length": sum(b - a for _, a, b in portions), "Portions": [[s, a, b] for s, a, b in portions], "Vertices": [], "Frame": frame,
                "Match": {"Status": "Missing"}}

    def cluster_side(self):
        """The facing metal: one SpatialEdgeCluster segment (feature 1) at y = 1 over x in [-2, 12];
        clusters are not sampled by the facing check, so only the isolated side is gated."""
        return self.segment(-2.0, 12.0, 1.0, 9, [(0.0, 14.0, 1)])

    def run_check(self, segments, features):
        from .facing_check import facing_check

        identification = {"Version": 2, "MatchingRadius": self.R, "Features": features, "Segments": segments, "Vertices": [], "Exclusions": [],
                          "Diagnostics": {"SamePriorityClaimOverlaps": 0}, "GeometryDigest": "digest"}
        return facing_check({"Version": 2, "MatchingRadius": self.R, "Identification": identification}, spacing=0.5)

    def test_sub_tolerance_mesh_segment_of_a_long_isolated_edge_is_not_exempt(self):
        """A 1 nm mesh segment inside a long isolated edge (both run neighbours claimed by the same
        feature) is an ordinary portion: it faces the cluster metal and is NOT exempt (DS-SCT-001
        at R 2.0 had 896 such segments; review fix-4 m-A)."""
        tiny = self.TINY
        segments = [self.segment(0.0, 4.0, 0.0, 0, [(0.0, 4.0, 0)]),
                    self.segment(4.0, 4.0 + tiny, 0.0, 0, [(0.0, tiny, 0)]),
                    self.segment(4.0 + tiny, 10.0, 0.0, 0, [(0.0, 6.0 - tiny, 0)]), self.cluster_side()]
        features = [self.feature(0, "IsolatedEdge", [(0, 0.0, 4.0), (1, 0.0, tiny), (2, 0.0, 6.0 - tiny)]),
                    self.feature(1, "SpatialEdgeCluster", [(3, 0.0, 14.0)])]
        result = self.run_check(segments, features)
        self.assertFalse(result["Gates"]["IsolatedFacing"])
        self.assertEqual(result["SubToleranceFeatures"]["Count"], {})
        self.assertEqual(result["SubToleranceFeatures"]["ExemptLength"], 0.0)
        self.assertNotIn("SubToleranceFeature", result["Isolated"]["ExcludedLength"])
        self.assertAlmostEqual(result["Isolated"]["UnexcludedFacingLength"], 10.0, places=9)

    def test_sub_tolerance_remainder_between_other_claims_is_exempt(self):
        """A sub-tolerance isolated portion bounded on both sides along its run by other features'
        claims (within one segment, and as a whole segment between two claimed neighbours) is the
        unclaimed remainder of the assignment: exempt as SubToleranceFeature, bounded by the
        tolerance per portion, counted."""
        tiny = self.TINY
        segments = [self.segment(0.0, 4.0, 0.0, 0, [(0.0, 2.0, 1), (2.0, 2.0 + tiny, 0), (2.0 + tiny, 4.0, 1)]),
                    self.segment(4.0, 6.0, 0.0, 0, [(0.0, 2.0, 1)]),
                    self.segment(6.0, 6.0 + tiny, 0.0, 0, [(0.0, tiny, 0)]),
                    self.segment(6.0 + tiny, 10.0, 0.0, 0, [(0.0, 4.0 - tiny, 1)]),
                    # The rest of the isolated feature far away (nothing to face), so that the
                    # feature as a whole is not sub-tolerance.
                    self.segment(100.0, 110.0, 0.0, 5, [(0.0, 10.0, 0)]), self.cluster_side()]
        features = [self.feature(0, "IsolatedEdge", [(0, 2.0, 2.0 + tiny), (2, 0.0, tiny), (4, 0.0, 10.0)]),
                    self.feature(1, "SpatialEdgeCluster", [(0, 0.0, 2.0), (0, 2.0 + tiny, 4.0), (1, 0.0, 2.0), (3, 0.0, 4.0 - tiny), (5, 0.0, 14.0)])]
        result = self.run_check(segments, features)
        self.assertTrue(result["Gates"]["IsolatedFacing"])
        self.assertEqual(result["SubToleranceFeatures"]["Count"], {"IsolatedEdge": 2})
        self.assertAlmostEqual(result["SubToleranceFeatures"]["ExemptLength"], 2 * tiny, places=12)
        self.assertAlmostEqual(result["Isolated"]["ExcludedLength"]["SubToleranceFeature"], 2 * tiny, places=12)
        self.assertEqual(result["Isolated"]["UnexcludedFacingLength"], 0.0)

    def test_whole_sub_tolerance_feature_is_exempt(self):
        """An isolated feature that is as a whole shorter than the tolerance (a 1 nm chain with no
        claimed neighbour along its run) is exempt; the same sliver with a longer sibling portion
        elsewhere is not (its run neighbours are not other features' claims)."""
        tiny = self.TINY
        segments = [self.segment(5.0, 5.0 + tiny, 0.0, 0, [(0.0, tiny, 0)]), self.cluster_side()]
        features = [self.feature(0, "IsolatedEdge", [(0, 0.0, tiny)]), self.feature(1, "SpatialEdgeCluster", [(1, 0.0, 14.0)])]
        result = self.run_check(segments, features)
        self.assertTrue(result["Gates"]["IsolatedFacing"])
        self.assertEqual(result["SubToleranceFeatures"]["Count"], {"IsolatedEdge": 1})
        self.assertAlmostEqual(result["SubToleranceFeatures"]["ExemptLength"], tiny, places=12)
        segments = [self.segment(5.0, 5.0 + tiny, 0.0, 0, [(0.0, tiny, 0)]), self.segment(100.0, 110.0, 0.0, 5, [(0.0, 10.0, 0)]), self.cluster_side()]
        features = [self.feature(0, "IsolatedEdge", [(0, 0.0, tiny), (1, 0.0, 10.0)]), self.feature(1, "SpatialEdgeCluster", [(2, 0.0, 14.0)])]
        result = self.run_check(segments, features)
        self.assertFalse(result["Gates"]["IsolatedFacing"])
        self.assertEqual(result["SubToleranceFeatures"]["Count"], {})
        self.assertAlmostEqual(result["Isolated"]["UnexcludedFacingLength"], tiny, places=12)


if __name__ == "__main__":
    unittest.main()
