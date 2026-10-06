# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""derive_semantic_contract reproduces the frozen four-/ten-edge contracts from their
immutable inputs and a label census, and fails closed on inconsistent inputs."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from derive_semantic_contract import coupon_radius, derive, semantic_corners
from semantic_mesh_contract import (INVARIANT_CORNER_RULE, box_face_cut_end, boundary_semantic_corners,
                                    invariant_corner, plan_view_quantum, quantised_side_dot)

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
