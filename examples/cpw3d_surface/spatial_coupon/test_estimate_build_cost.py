# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""estimate_build_cost: the pre-build element estimate follows the recipe's size laws
from the frozen inputs, reproduces its recorded calibration and fails closed against
the element cap before any build (headroom gate)."""
import copy
import json
from pathlib import Path
import unittest

from estimate_build_cost import (actual_counts, case_paths, coupon_box, edge_chains, estimate, gate,
                                 read_edges, tube_rings)
from general_mesh_manifest import (BUILD_COST_ESTIMATE_KEY, preflight_build_cost, validate_build_cost_estimate_model,
                                   validate_manifest)

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "geometry-independence-suite.json"
# Recorded tetrahedron counts of the verified calibration roots (manifest
# ProductionRecipe.BuildCostEstimate.Calibration).
CALIBRATION_TETRAHEDRA = {"four-edge-9d2cb9bbb3fe": 1731140, "ten-edge-6791f1c84123": 3208149,
                          "three-edge-419576fdab24": 1659373, "two-edge-8dd4bc70f183": 425916,
                          "two-edge-3f8992613e95": 1322026, "three-edge-current-calibration": 1311526,
                          "concave-multislot": 1423816}


class EstimateBuildCostTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads(MANIFEST.read_text())
        cls.options = cls.manifest["ProductionRecipe"]["BuildCommandOptions"]
        cls.model = cls.manifest["ProductionRecipe"][BUILD_COST_ESTIMATE_KEY]

    def case(self, case_id):
        return next(case for case in self.manifest["Cases"] if case["Id"] == case_id)

    def test_box_and_tube_recipe_follow_the_mesher(self):
        paths = case_paths(self.manifest, MANIFEST, self.case("four-edge-9d2cb9bbb3fe"))
        lower, upper = coupon_box(read_edges(paths["Signature"]), 2.0, 0.1, 0.05)
        self.assertEqual((lower.tolist(), upper.tolist()), ([-6.0, -8.0, -2.05], [10.0, 8.0, 2.1]))
        sizes, radius = tube_rings(0.00025, 2.0, 0.05, 0.1, 0.1)
        self.assertEqual(len(sizes), 7)
        self.assertAlmostEqual(radius, 0.03175)
        result = estimate(paths, self.options, self.model["TetrahedraPerCubicSize"])
        tubes = result["Tubes"]
        self.assertEqual((tubes["Sides"], tubes["Rings"], tubes["Sectors"]), (4, 7, 9))
        self.assertEqual(tubes["Prisms"], tubes["Layers"] * 117)
        self.assertEqual(tubes["Pyramids"], tubes["Layers"] * 9)
        self.assertEqual(result["Corners"], 4)
        self.assertEqual(result["Sizes"]["TangentialSize"], 0.05)

    def test_cad_subdivided_edge_keeps_the_unsubdivided_box(self):
        # Decision 47: the two collinear half rows of one-edge-cad-subdivided chain into
        # the straight edge's union, so both fixtures share the box and the estimate.
        straight = estimate(case_paths(self.manifest, MANIFEST, self.case("one-edge-straight")), self.options,
                            self.model["TetrahedraPerCubicSize"])
        subdivided = estimate(case_paths(self.manifest, MANIFEST, self.case("one-edge-cad-subdivided")),
                              self.options, self.model["TetrahedraPerCubicSize"])
        self.assertEqual(subdivided["Box"], straight["Box"])
        self.assertEqual(straight["Box"]["Lower"][:2], [-4.0, -8.0])
        self.assertAlmostEqual(subdivided["EstimatedElements"], straight["EstimatedElements"])
        ten = read_edges(case_paths(self.manifest, MANIFEST, self.case("ten-edge-6791f1c84123"))["Signature"])
        self.assertEqual(edge_chains(ten), {})      # collinear rows separated by a slot never chain

    def test_size_bound_enters_the_estimate(self):
        result = estimate(case_paths(self.manifest, MANIFEST, self.case("concave-multislot")),
                          self.options, self.model["TetrahedraPerCubicSize"])
        self.assertEqual(result["Sizes"]["FarSize"], 0.04)
        self.assertEqual(result["Sizes"]["TangentialSize"], 0.04)
        self.assertEqual(result["Sizes"]["RequestedTangentialSize"], 0.05)
        self.assertIsNone(result["TraceBasis"])

    def test_reproduces_the_recorded_calibration_within_ten_percent(self):
        ratios = []
        for case_id, tetrahedra in CALIBRATION_TETRAHEDRA.items():
            result = estimate(case_paths(self.manifest, MANIFEST, self.case(case_id)), self.options,
                              self.model["TetrahedraPerCubicSize"])
            ratio = result["EstimatedTetrahedra"] / tetrahedra
            ratios.append(ratio)
            self.assertGreaterEqual(ratio, 0.999, case_id)      # never below a calibrated build
            self.assertLessEqual(ratio, 1.11, case_id)
        self.assertAlmostEqual(min(ratios), 1.0, places=2)
        recorded = self.model["Calibration"]["IntegralOverActualTetrahedra"]
        self.assertAlmostEqual(min(recorded.values()), self.model["TetrahedraPerCubicSize"], places=6)

    def test_needles_dominate_a_needle_heavy_basis(self):
        ten = estimate(case_paths(self.manifest, MANIFEST, self.case("ten-edge-6791f1c84123")), self.options,
                       self.model["TetrahedraPerCubicSize"])
        self.assertGreater(ten["Integrals"]["TraceBasis"], 0.5 * ten["Integrals"]["FarField"])
        self.assertLess(ten["TraceBasis"]["MinimumRequestedSize"], 0.003)

    def test_gate_fails_closed_against_the_cap_before_building(self):
        result = gate(self.manifest, MANIFEST, self.case("ten-edge-6791f1c84123"))
        self.assertTrue(result["Passed"])
        self.assertLess(result["EstimateOverCap"], 1.0)
        tight = copy.deepcopy(self.manifest)
        tight["Gates"]["MaximumElements"] = 3000000
        result = gate(tight, MANIFEST, self.case("ten-edge-6791f1c84123"))
        self.assertFalse(result["Passed"])
        self.assertGreater(result["EstimateOverCap"], 1.0)
        self.assertFalse(preflight_build_cost(tight, MANIFEST, self.case("ten-edge-6791f1c84123"))["Passed"])
        self.assertIsNone(preflight_build_cost({"Gates": tight["Gates"]}, MANIFEST, self.case("ten-edge-6791f1c84123")))

    def test_manifest_requires_the_model(self):
        validate_manifest(self.manifest, MANIFEST)
        broken = copy.deepcopy(self.manifest)
        del broken["ProductionRecipe"][BUILD_COST_ESTIMATE_KEY]
        with self.assertRaisesRegex(ValueError, "BuildCostEstimate model"):
            validate_manifest(broken, MANIFEST)
        with self.assertRaisesRegex(ValueError, "BuildCostEstimate model"):
            validate_build_cost_estimate_model({BUILD_COST_ESTIMATE_KEY: {**self.model, "TetrahedraPerCubicSize": 0.0}})

    def test_actual_counts_read_the_census_quality(self):
        census = HERE / "testdata" / "does-not-exist.json"
        with self.assertRaises(OSError):
            actual_counts(census)


if __name__ == "__main__":
    unittest.main()
