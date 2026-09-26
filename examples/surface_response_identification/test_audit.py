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
        perimeter = P.extract_perimeter(self.mesh(True), CONFIG)
        interactions = P.edge_interactions(perimeter, 2.0)
        # A's right edge x = 4 faces B's left edge x = 6 across the SA gap: a parallel pair at
        # separation exactly R = 2 (the transmon shield-strip case).
        parallel = [(d, c) for _, _, d, c in interactions if c >= 1.0 - 1e-8]
        self.assertIn(2.0, [round(d, 12) for d, _ in parallel])
        self.assertTrue(all(d <= 4.0 + 1e-9 for _, _, d, _ in interactions))

    def test_corner_snap_rule_matches_classifier(self):
        # A polyline vertex turning by 29.999 deg is REGULAR, one turning by 30.001 is a
        # CORNER (metaledge.cpp: 30 degree tolerance on the 1e-12 direction grid).
        for turn, expected in ((29.999, "REGULAR"), (30.001, "CORNER"), (30.0, "REGULAR")):
            angle = np.radians(turn)
            nodes = [(0, 0, 0), (1, 0, 0), (1 + np.cos(angle), np.sin(angle), 0), (0, 1, 0), (1, 1, 0), (1 + np.cos(angle), 1 + np.sin(angle), 0)]
            elements = [(2, 5, (1, 2, 5)), (2, 5, (1, 5, 4)), (2, 5, (2, 3, 6)), (2, 5, (2, 6, 5))]
            path = os.path.join(self.directory.name, "turn.msh2")
            write_msh2(path, nodes, elements, [(2, 5, "metal")], True)
            perimeter = P.extract_perimeter(read_msh2(path), CONFIG)
            vertex = perimeter.vertices[[i for i, v in enumerate(perimeter.vertices) if np.allclose(v.point, (1, 0, 0))][0]]
            self.assertEqual(vertex.physical_kind, expected, msg=f"turn {turn}")

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
            "Conventions": {"CornerTurnToleranceDegrees": 30.0},
            "ReferenceProcessNormal": [0.0, 0.0, 1.0],
            "Features": features,
            "Segments": segments,
            "Vertices": vertices,
            "Exclusions": [{"Class": k, "Reason": "test", "Count": exclusions[k][0], "Length": exclusions[k][1]} for k in order],
            "Totals": {"PerimeterLength": assigned + excluded, "AssignedLength": assigned, "ExcludedLength": excluded},
            "GeometryDigest": "digest",
        }

    def v2_manifest(self, identification):
        manifest = make_manifest([requirement("IsolatedEdge", 12, 16.0, {}), requirement("ConvexCorner", 8, 32.0, {"AngleDegrees": 90.0, "CornerRadius": 0.0})])
        manifest["Version"] = 2
        manifest["Identification"] = identification
        return manifest

    def test_complete_identification_passes(self):
        code, gates, result = self.run_gates(self.v2_manifest(self.identification()))
        for name in ("perimeter-agreement", "A1-length-partition", "A1-count-partition", "A1-multiplicity", "A1-weights", "A1-vertex-census", "A1-exclusions-recorded", "A2-cluster-balls"):
            self.assertEqual(gates[name], "PASS", name)
        self.assertEqual(code, 0)
        self.assertEqual(result["Identification"]["Features"]["IsolatedEdge"], 12)
        self.assertAlmostEqual(result["Identification"]["Totals"]["AssignedLength"], 16.0)
        self.assertEqual(result["ByGate"]["A1-vertex-census"]["ManifestExcludedVertices"], 4)

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


if __name__ == "__main__":
    unittest.main()
