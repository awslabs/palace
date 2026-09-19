# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""derive_semantic_contract reproduces the frozen four-/ten-edge contracts from their
immutable inputs and a label census, and fails closed on inconsistent inputs."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from derive_semantic_contract import derive

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
