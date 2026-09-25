#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the identification audit on a synthetic tiny mesh with known answers."""

import contextlib
import io
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


if __name__ == "__main__":
    unittest.main()
