# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""derive_semantic_contract reproduces the frozen four-/ten-edge contracts from their
immutable inputs and a label census, and fails closed on inconsistent inputs."""
import copy
import csv
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest

from derive_semantic_contract import coupon_radius, derive, semantic_corners
from semantic_mesh_contract import (CORNER_FREE_RULE, INVARIANT_CORNER_RULE, box_face_cut_end,
                                    boundary_semantic_corners, contract_is_corner_free, invariant_corner,
                                    plan_view_quantum, quantised_side_dot, validate_semantic_contract)

HERE = Path(__file__).resolve().parent
TEN_EDGE = HERE / "testdata" / "ten-edge-6791f1c84123"
FOUR_EDGE = HERE / "testdata" / "four-edge-9d2cb9bbb3fe"
# The plan-view quantum of the Radius-1.9 production coupons (1.9e-9) and of the Radius-2
# repository fixtures.
QUANTUM_1_9 = plan_view_quantum(1.9)
QUANTUM_2 = plan_view_quantum(2.0)


def census_file(directory, labels):
    path = Path(directory) / "build-census.json"
    path.write_text(json.dumps({"InterfaceAreas": [{"Attribute": label, "Area": 1.0}
                                                   for label in labels]}))
    return path


class DeriveSemanticContractTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_reproduces_the_frozen_ten_edge_contract(self):
        frozen = json.loads((TEN_EDGE / "semantic-contract.json").read_text())
        labels = [item["Attribute"] for item in frozen["BoundaryLabels"]]
        derived = derive(TEN_EDGE, census_file(self.tmp, labels))
        self.assertEqual(list(derived), list(frozen))
        for key in frozen:
            if key != "Derivation":
                self.assertEqual(derived[key], frozen[key], key)
        self.assertEqual(derived["Derivation"]["ProcessLibrarySHA256"],
                         frozen["Derivation"]["ProcessLibrarySHA256"])
        self.assertEqual(derived["Derivation"]["PlanViewBoundarySHA256"],
                         frozen["Derivation"]["PlanViewBoundarySHA256"])
        self.assertEqual(derived["Derivation"]["BuildCensusInterfaceLabels"], sorted(labels))
        self.assertNotIn("Provisional", derived["Derivation"])

    def test_reproduces_the_four_edge_attributes_with_the_unetched_plane(self):
        # The four-edge device footprint leaves an un-etched plane (3000); the frozen
        # contract's older role strings differ, every attribute, adjacency, corner,
        # support count and topology is reproduced.
        frozen = json.loads((FOUR_EDGE / "semantic-contract.json").read_text())
        labels = [item["Attribute"] for item in frozen["BoundaryLabels"]]
        derived = derive(FOUR_EDGE, census_file(self.tmp, labels))
        self.assertEqual([(item["Attribute"], item["AdjacentMaterials"], item.get("Protected"))
                          for item in derived["BoundaryLabels"]],
                         [(item["Attribute"], item["AdjacentMaterials"], item.get("Protected"))
                          for item in frozen["BoundaryLabels"]])
        self.assertEqual(derived["SemanticCorners"], frozen["SemanticCorners"])
        self.assertEqual(derived["FeatureTopology"], frozen["FeatureTopology"])
        self.assertEqual(derived["MetricSurfaceRoles"], frozen["MetricSurfaceRoles"])
        self.assertEqual(derived["BoundaryLabels"][1]["Role"], "un-etched-substrate-vacuum-slot-0")

    def test_provisional_contract_omits_the_unetched_plane_and_says_so(self):
        derived = derive(FOUR_EDGE)
        self.assertEqual([item["Attribute"] for item in derived["BoundaryLabels"]],
                         [1, 3100, 5001, 6001])
        self.assertIn("Provisional", derived["Derivation"])

    def test_census_labels_outside_the_families_or_missing_fail_closed(self):
        with self.assertRaisesRegex(ValueError, r"outside the derived families \[4001\]"):
            derive(FOUR_EDGE, census_file(self.tmp, [1, 3100, 4001, 5001, 6001]))
        with self.assertRaisesRegex(ValueError, r"required labels missing \[6001\]"):
            derive(FOUR_EDGE, census_file(self.tmp, [1, 3100, 5001]))
        with self.assertRaisesRegex(ValueError, "records no InterfaceAreas"):
            derive(FOUR_EDGE, census_file(self.tmp, []))

    def test_fixture_without_process_library_takes_its_pairs_from_the_signature(self):
        # The synthetic six-edge fixture binds named files in the shared testdata
        # directory and no process library: two slots x two conductors from the
        # signature, ten Physical corners, the source recorded, provisional without a census.
        derived = derive(HERE / "testdata", signature=HERE / "testdata" / "six-edge-cluster-signature.csv",
                         boundary=HERE / "testdata" / "six-edge-cluster-boundary.csv", radius=2.0)
        self.assertEqual([item["Attribute"] for item in derived["BoundaryLabels"]],
                         [1, 3100, 3101, 5001, 5002, 5101, 5102, 6001, 6002, 6101, 6102])
        self.assertEqual(len(derived["SemanticCorners"]), 10)
        self.assertEqual(derived["FeatureTopology"]["BoundaryPhysicalVertexCount"], 10)
        self.assertIsNone(derived["Derivation"]["ProcessLibrarySHA256"])
        self.assertIn("no process library", derived["Derivation"]["SlotConductorSource"])
        self.assertIn("Provisional", derived["Derivation"])
        frozen = json.loads((HERE / "testdata" / "six-edge-semantic.json").read_text())
        self.assertEqual({k: v for k, v in derived.items() if k != "Derivation"},
                         {k: v for k, v in frozen.items() if k != "Derivation"})
        with self.assertRaisesRegex(FileNotFoundError, "mesh-signature.csv"):
            derive(HERE / "testdata", radius=2.0)
        # A source without process.toml needs the radius (the plan-view quantum 1e-9 R).
        with self.assertRaisesRegex(ValueError, "process.toml"):
            derive(HERE / "testdata", signature=HERE / "testdata" / "six-edge-cluster-signature.csv",
                   boundary=HERE / "testdata" / "six-edge-cluster-boundary.csv")

    def test_coupon_radius_comes_from_process_toml_or_the_explicit_radius(self):
        # The frozen fixtures: process.toml Radius 2.0 (the number the mesher receives as
        # --radius); an explicit radius must agree with it.
        self.assertEqual(coupon_radius(FOUR_EDGE, None), 2.0)
        self.assertEqual(coupon_radius(FOUR_EDGE, 2.0), 2.0)
        with self.assertRaisesRegex(ValueError, "differs from"):
            coupon_radius(FOUR_EDGE, 1.9)
        self.assertEqual(coupon_radius(self.tmp, 1.9), 1.9)
        with self.assertRaisesRegex(ValueError, "process.toml"):
            coupon_radius(self.tmp, None)
        (self.tmp / "process.toml").write_text('Units = "um"\nRadius = 0.0\n')
        with self.assertRaisesRegex(ValueError, "Radius must be a positive number"):
            coupon_radius(self.tmp, None)

    def test_multi_model_or_non_cluster_library_is_a_contract_error(self):
        broken = self.tmp / "source"
        shutil.copytree(FOUR_EDGE, broken)
        library_path = broken / "process-library.json"
        library = json.loads(library_path.read_text())
        library_path.write_text(json.dumps({**library, "Models": library["Models"] * 2}))
        with self.assertRaisesRegex(ValueError, "exactly one model .* found 2"):
            derive(broken)
        arms = json.loads(json.dumps(library))
        arms["Models"][0]["Topology"] = "Arms"
        library_path.write_text(json.dumps(arms))
        with self.assertRaisesRegex(ValueError, "topology 'Arms' is not SpatialEdgeCluster"):
            derive(broken)
        bare = json.loads(json.dumps(library))
        del bare["Models"][0]["Edges"][0]["InterfaceSlot"]
        library_path.write_text(json.dumps(bare))
        with self.assertRaisesRegex(ValueError, "must carry InterfaceSlot and Conductor"):
            derive(broken)

    def test_signature_and_process_library_pairs_must_agree(self):
        broken = self.tmp / "source"
        shutil.copytree(FOUR_EDGE, broken)
        signature = (broken / "mesh-signature.csv").read_text().splitlines()
        signature[1] = signature[1].replace(",0,1,", ",1,1,", 1)
        (broken / "mesh-signature.csv").write_text("\n".join(signature) + "\n")
        with self.assertRaisesRegex(ValueError, "differ from the signature"):
            derive(broken)


if __name__ == "__main__":
    unittest.main()


def write_boundary(path, loops):
    """loops: [[(class, x, y), ...], ...] in loop order (a vertex carries its outgoing side's
    class, as the plan-view builder writes it)."""
    lines = ["Loop,Vertex,Conductor,Plane,Hole,Class,X,Y"]
    for loop_index, loop in enumerate(loops, start=1):
        for vertex, (cls, x, y) in enumerate(loop, start=1):
            lines.append(f"{loop_index},{vertex},1,0.0,0,{cls},{x!r},{y!r}")
    Path(path).write_text("\n".join(lines) + "\n")
    return path


class BoxFaceCutEndTest(unittest.TestCase):
    """Supervisor decision 320 (block (b) design A2 / A6): a Physical vertex whose incoming
    side is a Continuation (box) side has a single metal side there; it is a box-face cut
    end - not a semantic corner - unless its side is exactly perpendicular to the face."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_oblique_box_vertex_is_a_cut_end_and_a_perpendicular_one_a_legacy_corner(self):
        # The O1 2-edge geometry (coupons-a probe source): vertex 1 (5.9475643, -10.9262198) on
        # the face y0 after the Continuation side 5 -> 1, its side leaving at 22.5 degrees.
        o1 = [("Physical", 5.9475643, -10.9262198), ("Physical", 3.2135498, -4.3257243),
              ("Physical", -0.5864502, 4.848287), ("Physical", -1.3135498, 6.603659),
              ("Continuation", -1.3135498, -10.9262198)]
        self.assertTrue(box_face_cut_end((o1[4][1], o1[4][2]), "Continuation", (o1[0][1], o1[0][2]),
                                         "Physical", (o1[1][1], o1[1][2])))
        corners, cut_ends, invariant = boundary_semantic_corners(
            [{"Loop": "1", "Class": cls, "X": repr(x), "Y": repr(y), "Plane": "0.0"} for cls, x, y in o1],
            QUANTUM_1_9)
        self.assertEqual(cut_ends, [[5.9475643, -10.9262198, 0.0]])
        self.assertEqual([c[:2] for c in corners],
                         [[3.2135498, -4.3257243], [-0.5864502, 4.848287], [-1.3135498, 6.603659]])
        # Design round 2 F5-A: the two oblique kinks are invariant corners; the vertex 4 where the
        # oblique side meets the box side x = -1.3135498 is not exactly perpendicular either.
        self.assertEqual(invariant, corners)
        # The C2 19-edge loop-2 vertices 14 / 15: the side leaving (-1.3574569, 12.9453859) runs
        # along x = -1.3574569 exactly -> theta == 0 -> the legacy corner (decision 311's
        # "SemanticCorner on face y1"), bitwise unchanged.
        self.assertFalse(box_face_cut_end((-0.9484572, 12.9453859), "Continuation",
                                          (-1.3574569, 12.9453859), "Physical", (-1.3574569, 12.4703859)))
        # Two Physical sides meeting at a box vertex stay a corner at any angle; a Continuation
        # vertex is never a corner.
        self.assertFalse(box_face_cut_end((1.0, 2.0), "Physical", (4.0, 0.0), "Physical", (1.0, -2.0)))
        self.assertFalse(box_face_cut_end((1.0, 0.0), "Physical", (4.0, 0.0), "Continuation", (4.0, 4.0)))
        # A Continuation side must run along one box face.
        with self.assertRaisesRegex(ValueError, "along one box face"):
            box_face_cut_end((1.0, 1.0), "Continuation", (4.0, 0.0), "Physical", (1.0, -2.0))

    def test_derived_contract_classifies_the_oblique_box_vertex_as_a_cut_endpoint(self):
        # Two 45-degree metal tips at the origin inside the box [-4, 4]^2 whose oblique side
        # runs from the tip to the box corner-free point (4, 4) of the face x = 4 ... the same
        # oblique side with the metal on either side of it: the vertex class at (4, 4) is that
        # of the OUTGOING side, so one orientation reads Continuation there and the other
        # Physical; both are box-face cut ends of the topology, never corners.
        source = self.tmp / "tip"
        source.mkdir()
        # (a) metal = the triangle (0, 0) -> (4, 4) -> (0, 4): side A leaves at 45 degrees
        #     (Continuation class at (4, 4)), side B along +y leaves (0, 4) perpendicularly
        #     (Physical class after the box side y = 4: the legacy corner, theta == 0).
        rows = ["Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm",
                "1,0,1,2,2,0,0.70710678118654757,-0.70710678118654757,0,0.70710678118654757,0.70710678118654757,0,1,-2.8284271247461903,2.8284271247461903,0",
                "2,0,1,0,2,0,-1,0,0,0,-1,0,1,-2,2,0"]
        (source / "mesh-signature.csv").write_text("\n".join(rows) + "\n")
        boundary = write_boundary(source / "plan-view-boundary.csv",
                                  [[("Physical", 0.0, 0.0), ("Continuation", 4.0, 4.0), ("Physical", 0.0, 4.0)]])
        corners, cut_ends, invariant = semantic_corners(boundary, QUANTUM_2)
        self.assertEqual(corners, [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(cut_ends, [])
        # Design round 2 F5-A: the 45-degree tip is an invariant corner (sides (4, 4) and
        # (0, 4) away from it: dot product 16 != 0); the theta-0 box vertex (0, 4) is legacy
        # (sides (0, -4) and (4, 0): exactly 0). Recorded under Derivation only where it acts.
        self.assertEqual(invariant, [[0.0, 0.0, 0.0]])
        derived = derive(source, radius=2.0)
        self.assertEqual(derived["SemanticCorners"], [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(derived["FeatureTopology"]["CutEndpoints"], [[4.0, 4.0, 0.0]])
        self.assertNotIn("BoxFaceCutEnds", derived["Derivation"])
        self.assertEqual(derived["Derivation"]["InvariantCorners"]["Points"], [[0.0, 0.0, 0.0]])
        self.assertIn("kappa_reg", derived["Derivation"]["InvariantCorners"]["Rule"])
        # (b) metal = the triangle (0, 0) -> (4, 0) -> (4, 4): the box side x = 4 ARRIVES at
        #     (4, 4), whose outgoing side is the oblique one -> Physical class, a single metal
        #     side, theta 45 > 0: a box-face cut end by decision 320 (the legacy rule would have
        #     made it a corner with h_K + a ball on the face), recorded under Derivation.
        rows = ["Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm",
                "1,0,1,2,0,0,0,-1,0,1,0,0,1,-2,2,0",
                "2,0,1,2,2,0,-0.70710678118654757,0.70710678118654757,0,-0.70710678118654757,-0.70710678118654757,0,1,-2.8284271247461903,2.8284271247461903,0"]
        (source / "mesh-signature.csv").write_text("\n".join(rows) + "\n")
        boundary = write_boundary(source / "plan-view-boundary.csv",
                                  [[("Physical", 0.0, 0.0), ("Continuation", 4.0, 0.0), ("Physical", 4.0, 4.0)]])
        corners, cut_ends, invariant = semantic_corners(boundary, QUANTUM_2)
        self.assertEqual(corners, [[0.0, 0.0, 0.0]])
        self.assertEqual(cut_ends, [[4.0, 4.0, 0.0]])
        self.assertEqual(invariant, [[0.0, 0.0, 0.0]])
        derived = derive(source, radius=2.0)
        self.assertEqual(derived["SemanticCorners"], [[0.0, 0.0, 0.0]])
        self.assertEqual(derived["FeatureTopology"]["CutEndpoints"], [[4.0, 0.0, 0.0], [4.0, 4.0, 0.0]])
        self.assertEqual(derived["Derivation"]["BoxFaceCutEnds"]["Points"], [[4.0, 4.0, 0.0]])
        self.assertIn("decision 320", derived["Derivation"]["BoxFaceCutEnds"]["Rule"])
        # (c) the same triangle with the oblique side exactly perpendicular instead - the
        #     legacy rectilinear case: (4, 4) -> (0, 4) is not possible here, so take the
        #     rectangle (0, 0) -> (4, 0) -> (4, 4) -> (0, 4): its Physical-class box vertex (0, 4)
        #     after the box side y = 4 leaves along x = 0 exactly -> a corner, nothing recorded.
        rows = ["Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm",
                "1,0,1,2,0,0,0,-1,0,1,0,0,1,-2,2,0",
                "2,0,1,0,2,0,-1,0,0,0,-1,0,1,-2,2,0"]
        (source / "mesh-signature.csv").write_text("\n".join(rows) + "\n")
        boundary = write_boundary(source / "plan-view-boundary.csv",
                                  [[("Physical", 0.0, 0.0), ("Continuation", 4.0, 0.0),
                                    ("Continuation", 4.0, 4.0), ("Physical", 0.0, 4.0)]])
        corners, cut_ends, invariant = semantic_corners(boundary, QUANTUM_2)
        self.assertEqual(corners, [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(cut_ends, [])
        # Rectilinear: no invariant corner, nothing recorded (every built contract unchanged).
        self.assertEqual(invariant, [])
        derived = derive(source, radius=2.0)
        self.assertNotIn("BoxFaceCutEnds", derived["Derivation"])
        self.assertNotIn("InvariantCorners", derived["Derivation"])
        self.assertEqual(derived["FeatureTopology"]["CutEndpoints"], [[4.0, 0.0, 0.0]])


# The seven semantic corners of the SCT loop end 1b26671c9080 (sct002-S1p, Radius 1.9; the
# merge-preparation registration PBS 57005) with their two plan-view loop neighbours as the
# plan-view boundary CSV spells them: the two 135-degree corners, the box corner, a rigidly
# ROTATED perpendicular corner (both sides tilted 9.5e-7 rad: (2.0000008, -0.0000019) and
# (0.0000019, 2.0000008), exactly 1000 quanta off the axes) and three corners whose sides
# are tilted by real ~1e-6-rad amounts (supervisor decision 416's erratum to 410).
LOOP_END_CORNERS = [
    ((10.8001909, -14.6336328), (14.6001928, -10.833638500000001), (14.600213700000001, 11.1663608), True),
    ((14.6001928, -10.833638500000001), (14.600213700000001, 11.1663608), (10.8002023, 14.966381700000001), True),
    ((-19.730892, 14.966381700000001), (-19.730892, 1.1663834), (-6.8997949, 1.1663834), False),
    ((-6.899793, 3.1663823), (-8.8997938, 3.1663842), (-8.8997919, 5.166385), False),
    ((-8.8997938, 3.1663842), (-8.8997919, 5.166385), (-7.8997915, 5.1663831), True),
    ((-8.8997919, 5.166385), (-7.8997915, 5.1663831), (-7.8997915, 6.1663835), True),
    ((-7.8998029, -5.8336156), (-7.899801, -4.8336171000000006), (-19.730892, -4.8336171000000006), True),
]


class ExactQuantisedPredicateTest(unittest.TestCase):
    """Supervisor decision 416: the invariant-corner predicate is integer arithmetic on the
    plan-view quantum counts, so a perpendicular corner is legacy whatever its rotation and
    a truly tilted corner invariant; the float dot product was exact for axis-aligned sides
    only."""

    def test_loop_end_corners_rotated_perpendicular_legacy_real_tilts_invariant(self):
        for previous, point, following, expected in LOOP_END_CORNERS:
            # Every coordinate is an integer count of the 1.9e-9 quantum.
            for value in (*previous, *point, *following):
                self.assertLess(abs(value / QUANTUM_1_9 - round(value / QUANTUM_1_9)), 1.0e-5)
            self.assertEqual(invariant_corner(previous, point, following, QUANTUM_1_9), expected, point)
        self.assertEqual(sum(expected for *_, expected in LOOP_END_CORNERS), 5)
        # The rotated perpendicular corner: the float dot product of the float side vectors
        # is -1.776e-15 (the coordinates' rounding), the quantised one exactly 0.
        previous, point, following, _ = LOOP_END_CORNERS[3]
        a = (previous[0] - point[0], previous[1] - point[1])
        b = (following[0] - point[0], following[1] - point[1])
        self.assertNotEqual(a[0] * b[0] + a[1] * b[1], 0.0)
        self.assertEqual(quantised_side_dot(previous, point, following, QUANTUM_1_9), 0)
        # The three real tilts: 526316000 x 1000 (twice) and 6226890000 x 1000 quanta squared.
        self.assertEqual([quantised_side_dot(p, q, f, QUANTUM_1_9) for p, q, f, _ in LOOP_END_CORNERS[4:]],
                         [526316000 * 1000, 526316000 * 1000, 6226890000 * 1000])
        # Rigid integer rotations on the quantum grid (3-4-5 / 5-12-13 triangles) stay legacy.
        for c, s in ((3, 4), (4, 3), (-4, 3), (5, 12)):
            q = 100_000_000
            point = (0.38, -0.76)
            previous = (point[0] + QUANTUM_1_9 * c * q, point[1] + QUANTUM_1_9 * s * q)
            following = (point[0] - QUANTUM_1_9 * s * q, point[1] + QUANTUM_1_9 * c * q)
            self.assertFalse(invariant_corner(previous, point, following, QUANTUM_1_9))

    def test_rows_classify_the_loop_end_corners_and_the_rule_names_the_quantum(self):
        # One three-vertex loop per corner (the triangle previous -> corner -> following; the
        # triangle's other two vertices are corners of their own sides, so each loop is read
        # alone and judged at its middle vertex).
        for previous, point, following, expected in LOOP_END_CORNERS:
            rows = [{"Loop": "1", "Vertex": str(vertex), "Class": "Physical", "X": repr(x), "Y": repr(y),
                     "Plane": "0.0"} for vertex, (x, y) in enumerate((previous, point, following), start=1)]
            corners, cut_ends, invariant = boundary_semantic_corners(rows, QUANTUM_1_9)
            self.assertEqual(len(corners), 3)
            self.assertEqual(cut_ends, [])
            self.assertEqual([*point, 0.0] in invariant, expected, point)
        self.assertIn("integer arithmetic", INVARIANT_CORNER_RULE)
        self.assertIn("1e-9 R", INVARIANT_CORNER_RULE)
        self.assertIn("rigidly rotated perpendicular corner", INVARIANT_CORNER_RULE)


def write_bent_band_source(directory, *, ground=False, tilt_degrees=2.5, chord_degrees=5.0):
    """The corner-free synthetics of mesher design round 3 class (1) (G.1.6), the Python twin of
    test_arc_tubes.jl write_bent_band_inputs: a metal band of width 0.5 (R 0.5) with one exactly
    tangent 90-degree bend (outer convex arc 0.75, inner CONCAVE arc 0.25 about (0, 0.75)) entering
    through the x0 face and leaving through the y1 face, rotated by `tilt_degrees` so that every
    straight leg crosses its face obliquely (a box-face cut end, decision 320); `ground` adds the
    region beyond a 0.25 gap whose edge carries a second concave arc of radius 1.0 (S-G1-b).
    Writes signature.csv / boundary.csv; returns (signature, boundary, chords)."""
    centre, inner, outer, ground_radius, sweep = (0.0, 0.75), 0.25, 0.75, 1.0, 90.0
    n = math.ceil(sweep / chord_degrees - 1e-9)
    c, s = math.cos(math.radians(tilt_degrees)), math.sin(math.radians(tilt_degrees))
    rotate = lambda p: (centre[0] + c * (p[0] - centre[0]) - s * (p[1] - centre[1]),
                        centre[1] + s * (p[0] - centre[0]) + c * (p[1] - centre[1]))
    west, north = (-c, -s), (-s, c)
    rows = []

    def straight(a, b):
        d = (b[0] - a[0], b[1] - a[1]); L = math.hypot(*d); t = (d[0] / L, d[1] / L)
        rows.append({"Slot": 0, "Conductor": 1, "Px": 0.5 * (a[0] + b[0]), "Py": 0.5 * (a[1] + b[1]), "Pz": 0.0,
                     "Gx": t[1], "Gy": -t[0], "Gz": 0.0, "Tx": t[0], "Ty": t[1], "Tz": 0.0, "Nz": 1,
                     "S0": -0.5 * L, "S1": 0.5 * L, "VertexArm": 0})

    def arc_points(rho):
        return [rotate((centre[0] + rho * math.cos(math.radians(-90.0 + sweep * k / n)),
                        centre[1] + rho * math.sin(math.radians(-90.0 + sweep * k / n)))) for k in range(n + 1)]
    outer_points, inner_points = arc_points(outer), arc_points(inner)[::-1]
    ground_points = arc_points(ground_radius)[::-1]
    along = lambda p, d, length: (p[0] + length * d[0], p[1] + length * d[1])
    chains = [(outer_points, west, north, outer, 1, 1), (inner_points, north, west, inner, -1, 2)]
    if ground:
        chains.append((ground_points, north, west, ground_radius, -1, 3))
    for points, before, after, _, _, _ in chains:
        straight(along(points[0], before, 1.0), points[0])
        for k in range(n):
            straight(points[k], points[k + 1])
        straight(points[-1], along(points[-1], after, 1.0))
    # The box of the mesher's row rule (row_coupon_bounds: every straight claim of length 2 R
    # extended by 2 R at both ends, padded by R transversally, then R): the legs reach the faces.
    xs, ys = [], []
    for row in rows:
        for end in (row["S0"] - 1.0 if row["S0"] <= -0.5 else row["S0"], row["S1"] + 1.0 if row["S1"] >= 0.5 else row["S1"]):
            for side in (-0.5, 0.5):
                xs.append(row["Px"] + end * row["Tx"] + side * row["Gx"])
                ys.append(row["Py"] + end * row["Ty"] + side * row["Gy"])
    lower, upper = (min(xs) - 0.5, min(ys) - 0.5), (max(xs) + 0.5, max(ys) + 0.5)
    x0_face = lambda p: (lower[0], p[1] + (lower[0] - p[0]) * west[1] / west[0])
    y1_face = lambda p: (p[0] + (upper[1] - p[1]) * north[0] / north[1], upper[1])
    band = [x0_face(outer_points[0]), *outer_points, y1_face(outer_points[-1]), y1_face(inner_points[0]),
            *inner_points, x0_face(inner_points[-1])]
    tags = {1: {**{2 + k: (1, outer, 1) for k in range(n)}, **{n + 5 + k: (2, inner, -1) for k in range(n)}}}
    joints = {1: {2, n + 2, n + 5, 2 * n + 5}}
    loops = [band]
    if ground:
        loops.append([y1_face(ground_points[0]), *ground_points, x0_face(ground_points[-1]), lower,
                      (upper[0], lower[1]), upper])
        tags[2] = {2 + k: (3, ground_radius, -1) for k in range(n)}
        joints[2] = {2, n + 2}
    on_face = lambda p, q: any(abs(p[d] - bound[d]) <= 1e-9 and abs(q[d] - bound[d]) <= 1e-9
                               for d in range(2) for bound in (lower, upper))
    signature = Path(directory) / "signature.csv"
    with signature.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["Index", *rows[0]])
        writer.writeheader()
        for index, row in enumerate(rows, start=1):
            writer.writerow({"Index": index, **row})
    boundary = Path(directory) / "boundary.csv"
    with boundary.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["Loop", "Vertex", "Conductor", "Plane", "Hole", "Class", "X", "Y", "ArcId", "ArcCx", "ArcCy",
                         "ArcR", "ArcSign", "JointTurn", "JointSmooth"])
        for loop_index, polygon in enumerate(loops, start=1):
            m = len(polygon)
            for i, point in enumerate(polygon, start=1):
                cls = "Continuation" if on_face(point, polygon[i % m]) else "Physical"
                tag = tags[loop_index].get(i)
                arc = [tag[0], repr(centre[0]), repr(centre[1]), repr(tag[1]), tag[2]] if tag else [""] * 5
                joint = ["0.0", "1"] if i in joints[loop_index] else ["", ""]
                writer.writerow([loop_index, i, 1, "0.0", 0, cls, repr(point[0]), repr(point[1]), *arc, *joint])
    return signature, boundary, n


class CornerFreeContractTest(unittest.TestCase):
    """Mesher design round 3 class (1) (decisions 491 / 510; DESIGN-part-G G.1.3 Option A): a
    boundary whose every Physical vertex is an arc vertex or a box-face cut end derives an EMPTY
    corner list with Derivation.CornerFree (integer counts only); the validator admits the empty
    list iff the record is present; a boundary with no Physical vertex at all keeps the legacy
    stop (S-G1-c)."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_corner_free_band_derives_an_empty_corner_list_with_the_record(self):
        for ground in (False, True):
            signature, boundary, chords = write_bent_band_source(self.tmp, ground=ground)
            contract = derive(self.tmp, signature=signature, boundary=boundary, radius=0.5)
            arcs = 3 if ground else 2
            self.assertEqual(contract["SemanticCorners"], [])
            self.assertEqual(contract["Derivation"]["CornerFree"],
                             {"Rule": CORNER_FREE_RULE, "ArcInteriorVertices": arcs * (chords - 1),
                              "SmoothJoints": 2 * arcs, "BoxFaceCutEnds": arcs})
            # One cut end per face crossing: the vertex LEAVING the face (its incoming side the
            # box side; the arriving vertex carries the Continuation class of its outgoing side).
            self.assertTrue(contract_is_corner_free(contract))
            self.assertEqual(len(contract["Derivation"]["BoxFaceCutEnds"]["Points"]), arcs)
            self.assertNotIn("InvariantCorners", contract["Derivation"])
            topology = contract["FeatureTopology"]
            self.assertEqual(topology["ArcFeatureCount"], arcs)
            self.assertEqual(topology["PhysicalFeatureCount"], 3 * arcs)
            self.assertEqual(len(topology["CutEndpoints"]), 2 * arcs)
            self.assertEqual(len(topology["CADSubdivisionEndpoints"]), 2 * arcs)
            self.assertEqual(topology["BoundaryContinuationVertexCount"], 2 if not ground else 6)
            # The whole record is integers and the rule: nothing platform-dependent (class (10)).
            self.assertTrue(all(isinstance(v, int) for k, v in contract["Derivation"]["CornerFree"].items()
                                if k != "Rule"))
        # The same outline crossing the faces exactly perpendicularly (tilt 0) keeps decision 320's
        # legacy corners at the face vertices: not corner-free, no record.
        signature, boundary, _ = write_bent_band_source(self.tmp, tilt_degrees=0.0)
        perpendicular = derive(self.tmp, signature=signature, boundary=boundary, radius=0.5)
        self.assertEqual(len(perpendicular["SemanticCorners"]), 2)
        self.assertNotIn("CornerFree", perpendicular["Derivation"])
        self.assertFalse(contract_is_corner_free(perpendicular))

    def test_validator_admits_the_empty_list_only_under_the_record_and_never_beside_corners(self):
        signature, boundary, _ = write_bent_band_source(self.tmp)
        contract = derive(self.tmp, signature=signature, boundary=boundary, radius=0.5)
        self.assertIs(validate_semantic_contract(contract), contract)
        legacy = copy.deepcopy(contract)
        del legacy["Derivation"]["CornerFree"]
        with self.assertRaisesRegex(ValueError, "at least one 3D point"):
            validate_semantic_contract(legacy)
        mixed = copy.deepcopy(contract)
        mixed["SemanticCorners"] = [[0.0, 0.0, 0.0]]
        with self.assertRaisesRegex(ValueError, "CornerFree is recorded on a contract with semantic corners"):
            validate_semantic_contract(mixed)
        for name, value in (("ArcInteriorVertices", -1), ("SmoothJoints", 2.0), ("BoxFaceCutEnds", True)):
            broken = copy.deepcopy(contract)
            broken["Derivation"]["CornerFree"][name] = value
            with self.assertRaisesRegex(ValueError, "CornerFree must record"):
                contract_is_corner_free(broken)
        no_arc = copy.deepcopy(contract)
        no_arc["Derivation"]["CornerFree"].update(ArcInteriorVertices=0, SmoothJoints=0)
        with self.assertRaisesRegex(ValueError, "at least one arc vertex"):
            contract_is_corner_free(no_arc)
        wrong_rule = copy.deepcopy(contract)
        wrong_rule["Derivation"]["CornerFree"]["Rule"] = "another rule"
        with self.assertRaisesRegex(ValueError, "CornerFree must record"):
            contract_is_corner_free(wrong_rule)
        # Every frozen fixture contract is corner-bearing: no record, the rule inert.
        for fixture in (TEN_EDGE, FOUR_EDGE):
            frozen = json.loads((fixture / "semantic-contract.json").read_text())
            self.assertFalse(contract_is_corner_free(frozen))
            self.assertIs(validate_semantic_contract(frozen), frozen)

    def test_boundary_without_any_physical_vertex_keeps_the_legacy_stop(self):
        # S-G1-c: a boundary of box sides only (no Physical vertex, no arc vertex) is malformed;
        # so is a straight oblique strip whose only Physical vertices are cut ends - the empty
        # set is admitted only when the ARC classification produced it.
        header = ["Loop", "Vertex", "Conductor", "Plane", "Hole", "Class", "X", "Y"]
        def write(rows, name):
            path = self.tmp / name
            with path.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(header)
                writer.writerows(rows)
            return path
        box_only = write([[1, i + 1, 1, "0.0", 0, "Continuation", repr(x), repr(y)]
                          for i, (x, y) in enumerate(((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)))],
                         "box-only.csv")
        with self.assertRaisesRegex(ValueError, "classifies no Physical vertex"):
            semantic_corners(box_only, QUANTUM_2)
        oblique = write([[1, 1, 1, "0.0", 0, "Physical", "-1.0", "-0.5"], [1, 2, 1, "0.0", 0, "Continuation", "1.0", "-0.3"],
                         [1, 3, 1, "0.0", 0, "Physical", "1.0", "0.3"], [1, 4, 1, "0.0", 0, "Continuation", "-1.0", "0.5"]],
                        "oblique-strip.csv")
        with self.assertRaisesRegex(ValueError, "classifies no Physical vertex"):
            semantic_corners(oblique, QUANTUM_2)
