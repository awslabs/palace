# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""derive_semantic_contract reproduces the frozen four-/ten-edge contracts from their
immutable inputs and a label census, and fails closed on inconsistent inputs."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from derive_semantic_contract import derive, semantic_corners
from semantic_mesh_contract import box_face_cut_end, boundary_semantic_corners

HERE = Path(__file__).resolve().parent
TEN_EDGE = HERE / "testdata" / "ten-edge-6791f1c84123"
FOUR_EDGE = HERE / "testdata" / "four-edge-9d2cb9bbb3fe"


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
                         boundary=HERE / "testdata" / "six-edge-cluster-boundary.csv")
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
            derive(HERE / "testdata")

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
        corners, cut_ends = boundary_semantic_corners(
            [{"Loop": "1", "Class": cls, "X": repr(x), "Y": repr(y), "Plane": "0.0"} for cls, x, y in o1])
        self.assertEqual(cut_ends, [[5.9475643, -10.9262198, 0.0]])
        self.assertEqual([c[:2] for c in corners],
                         [[3.2135498, -4.3257243], [-0.5864502, 4.848287], [-1.3135498, 6.603659]])
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
        corners, cut_ends = semantic_corners(boundary)
        self.assertEqual(corners, [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(cut_ends, [])
        derived = derive(source)
        self.assertEqual(derived["SemanticCorners"], [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(derived["FeatureTopology"]["CutEndpoints"], [[4.0, 4.0, 0.0]])
        self.assertNotIn("BoxFaceCutEnds", derived["Derivation"])
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
        corners, cut_ends = semantic_corners(boundary)
        self.assertEqual(corners, [[0.0, 0.0, 0.0]])
        self.assertEqual(cut_ends, [[4.0, 4.0, 0.0]])
        derived = derive(source)
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
        corners, cut_ends = semantic_corners(boundary)
        self.assertEqual(corners, [[0.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
        self.assertEqual(cut_ends, [])
        derived = derive(source)
        self.assertNotIn("BoxFaceCutEnds", derived["Derivation"])
        self.assertEqual(derived["FeatureTopology"]["CutEndpoints"], [[4.0, 0.0, 0.0]])
