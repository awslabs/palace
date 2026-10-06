# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""estimate_build_cost: the pre-build element estimate follows the recipe's size laws
from the frozen inputs, reproduces its recorded calibration and fails closed against
the element cap before any build (headroom gate)."""
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import tomllib
import unittest

from estimate_build_cost import (actual_counts, case_paths, coupon_box, edge_chains, estimate, face_end_layers,
                                 gate, metal_sides, read_edges, read_loops, read_support_box, row_coupon_box,
                                 semantic_corner_count, tube_rings)
from general_mesh_manifest import (BUILD_COST_ESTIMATE_KEY, preflight_build_cost, thin_build_options,
                                   validate_build_cost_estimate_model, validate_manifest)

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "geometry-independence-suite.json"
# Recorded tetrahedron counts of the verified calibration roots (manifest
# ProductionRecipe.BuildCostEstimate.Calibration).
CALIBRATION_TETRAHEDRA = {"four-edge-9d2cb9bbb3fe": 1731140, "ten-edge-6791f1c84123": 3208149,
                          "three-edge-419576fdab24": 1659373, "two-edge-8dd4bc70f183": 425916,
                          "two-edge-3f8992613e95": 1322026, "three-edge-current-calibration": 1311526,
                          "concave-multislot": 1423816}
# The two S1p device-plan (contract 3, decision 282) 19-edge coupons whose R2b registration
# failed closed at the 6 M cap under the legacy box rule (decision 290: the estimator
# ignored the SupportBox): their recorded probe inputs (SHA-256 of the register probe
# manifest's Source.Files; the mesh recipe is testdata/generality-mesh-recipe.json), the
# SupportBox, the R2b SupportBox estimates they must reproduce and the legacy-rule
# estimates they must no longer produce (elements, rounded).
DEVICE_PLAN_FIXTURES = {
    "spatial-19-edge-39ab2ffd68ec": {
        "Files": {"mesh-signature.csv": "10a7caf154c759ea676cb51ed7a7fa94d96ee8a95a3e7877a071d5dc2b4c4b71",
                  "plan-view-boundary.csv": "449ff67ab123a31329934f8f8a5247d9e4064f41b498c2f13ae835a5aa3eefbc",
                  "process.toml": "6d235486b76eb537b8569f97e3040324861204c08c3dd04aff3f27b4d972a569",
                  "process-library.json": "f5dfeccee521b32cae64e5c0c4130b0449eef0780ea53fa30d4fc85937e9af42",
                  "trace-vertices.csv": "c5ab9b7979a9139939dbf46a60392b0e5c2e93dfbab1d936a37e2dc780a49401",
                  "trace-triangles.csv": "2c4df741415df94c03efab5d1a5437a8a05594a0ad813876571a5441c275f2f1"},
        "SupportBox": [-11.5472367, -12.6033897, 13.1277631, 12.971609699999998],
        "Elements": {"fabricated": 5070435, "thin": 4660940},
        "LegacyRuleElements": {"fabricated": 8074017, "thin": 7575173}},
    "spatial-19-edge-c5952e69af66": {
        "Files": {"mesh-signature.csv": "d830f8911c6c839d25046c55a57a72985a183172687c712fe4b6e53f2165e140",
                  "plan-view-boundary.csv": "4fb9cc3661e5d7edb3829f3392fc928885bb677af090a2f11661201fe345a346",
                  "process.toml": "6d235486b76eb537b8569f97e3040324861204c08c3dd04aff3f27b4d972a569",
                  "process-library.json": "d0409fa77cfab1e0b671d5780226c06c0998ed5a0fdbacaaaf86ac05d0ac16d6",
                  "trace-vertices.csv": "67803d288f3855cb597f40a995e3b2ae29ac3c994606d4f300b32e15bc1566d1",
                  "trace-triangles.csv": "f2c5d28007940fa6737d3865f045b439a64f8c627f3be490c6c817a4eb068a55"},
        "SupportBox": [-12.080608499999999, -12.7173916, 12.5943913, 12.857607799999998],
        "Elements": {"fabricated": 5025642, "thin": 4627115},
        "LegacyRuleElements": {"fabricated": 8092891, "thin": 7588639}}}
DEVICE_PLAN_ROLES = {"Signature": "mesh-signature.csv", "Boundary": "plan-view-boundary.csv", "Process": "process.toml",
                     "ProcessLibrary": "process-library.json", "TraceVertices": "trace-vertices.csv",
                     "TraceTriangles": "trace-triangles.csv"}


def device_plan_paths(case_id):
    paths = {role: HERE / "testdata" / case_id / name for role, name in DEVICE_PLAN_ROLES.items()}
    paths["MeshRecipe"] = HERE / "testdata" / "generality-mesh-recipe.json"
    return paths


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

    def test_box_pads_each_layer_by_its_process_normal(self):
        # Decision 48: Overetch on the substrate side (-Nz), MetalThickness on the metal
        # side (+Nz); the opposed-layers fixture (upward at 0, downward at 0.6) is padded
        # by the trench above its downward layer, the hole fixture (upward only) as before.
        paths = case_paths(self.manifest, MANIFEST, self.case("opposed-layers"))
        lower, upper = coupon_box(read_edges(paths["Signature"]), 0.5, 0.06, 0.02)
        self.assertAlmostEqual(lower[2], -0.52)
        self.assertAlmostEqual(upper[2], 0.6 + 0.5 + 0.02)
        paths = case_paths(self.manifest, MANIFEST, self.case("hole"))
        lower, upper = coupon_box(read_edges(paths["Signature"]), 0.5, 0.08, 0.03)
        self.assertAlmostEqual(lower[2], -0.53)
        self.assertAlmostEqual(upper[2], 0.58)
        result = estimate(paths, self.options, self.model["TetrahedraPerCubicSize"])
        self.assertEqual(result["Tubes"]["Sides"], 8)      # the hole's four sides carry tubes too

    def test_device_plan_coupon_box_is_the_process_library_support_box(self):
        """Decision 290 (a): a device-plan coupon (contract 3) is estimated over the mesher's
        box (coupon_bounds with the process library's SupportBox: the plane box the
        SupportBox, the z extent the rows'), and a plan side closing the coupon on a
        SupportBox face is not a tubed metal side; the two S1p 19-edge coupons reproduce
        the R2b SupportBox estimates from their recorded inputs (fail-before: the legacy
        rule's 8.07 / 8.09 M over the 6 M cap, the R2b registration failure)."""
        cap = self.manifest["Gates"]["MaximumElements"]
        thin_options = thin_build_options(self.manifest["ProductionRecipe"])
        for case_id, recorded in DEVICE_PLAN_FIXTURES.items():
            paths = device_plan_paths(case_id)
            for name, digest in recorded["Files"].items():
                self.assertEqual(hashlib.sha256((HERE / "testdata" / case_id / name).read_bytes()).hexdigest(), digest,
                                 f"{case_id}/{name} is not the recorded probe input")
            process = tomllib.loads(paths["Process"].read_text())
            radius, thickness, overetch = process["Radius"], process["MetalThickness"], process["Overetch"]
            support_box = read_support_box(paths["ProcessLibrary"])
            self.assertEqual(support_box, recorded["SupportBox"])
            edges = read_edges(paths["Signature"])
            self.assertEqual((len(edges), sum(edge["context"] for edge in edges)), (32, 13))
            lower, upper = coupon_box(edges, radius, thickness, overetch, support_box=support_box)
            row_lower, row_upper = row_coupon_box(edges, radius, thickness, overetch)
            self.assertEqual(lower[:2].tolist() + upper[:2].tolist(), support_box)
            self.assertEqual((lower[2], upper[2]), (row_lower[2], row_upper[2]))
            self.assertAlmostEqual(lower[2], -radius - overetch)
            self.assertAlmostEqual(upper[2], radius + thickness)
            # The grown signature box lies inside the legacy rule's (equal on the face y1 here).
            self.assertTrue(all(lower[:2] >= row_lower[:2]) and all(upper[:2] <= row_upper[:2]))
            self.assertLess(float((upper - lower).prod()), float((row_upper - row_lower).prod()))
            loops = read_loops(paths["Boundary"])
            # 31 straight plan sides; the 5 closing the plan on the SupportBox faces carry no tube
            # (under the legacy box they all did).
            self.assertEqual(len(metal_sides(loops, lower, upper, 1e-8 * radius)[0]), 26)
            self.assertEqual(len(metal_sides(loops, row_lower, row_upper, 1e-8 * radius)[0]), 31)
            for kind, options in (("fabricated", self.options), ("thin", thin_options)):
                result = estimate(paths, options, self.model["TetrahedraPerCubicSize"], kind=kind)
                self.assertEqual(result["Box"]["Lower"][:2] + result["Box"]["Upper"][:2], support_box)
                self.assertAlmostEqual(result["Box"]["Volume"], float((upper - lower).prod()))
                self.assertAlmostEqual(result["Box"]["Volume"], 2492.699265066, places=6)
                self.assertEqual(result["Tubes"]["Sides"], 26)
                self.assertEqual(round(result["EstimatedElements"]), recorded["Elements"][kind], (case_id, kind))
                self.assertLess(result["EstimatedElements"], cap)
                self.assertGreater(recorded["LegacyRuleElements"][kind], cap)

    def test_context_rows_without_a_support_box_fail_closed(self):
        """The mesher's rule: device-plan rows (Context) need the SupportBox; a library of
        several models or a degenerate SupportBox is refused."""
        paths = device_plan_paths("spatial-19-edge-39ab2ffd68ec")
        edges = read_edges(paths["Signature"])
        with self.assertRaisesRegex(ValueError, "SupportBox"):
            coupon_box(edges, 1.9, 0.1, 0.05)
        with self.assertRaisesRegex(ValueError, "SupportBox"):
            estimate({k: v for k, v in paths.items() if k != "ProcessLibrary"}, self.options,
                     self.model["TetrahedraPerCubicSize"])
        library = json.loads(paths["ProcessLibrary"].read_text())
        box = library["Models"][0]["SupportBox"]
        with tempfile.TemporaryDirectory() as directory:
            broken = Path(directory) / "process-library.json"
            broken.write_text(json.dumps({**library, "Models": library["Models"] * 2}))
            with self.assertRaisesRegex(ValueError, "exactly one model"):
                read_support_box(broken)
            broken.write_text(json.dumps({**library, "Models": [{**library["Models"][0],
                                                                  "SupportBox": [box[2], box[1], box[0], box[3]]}]}))
            with self.assertRaisesRegex(ValueError, "not a box"):
                read_support_box(broken)

    def test_legacy_coupon_estimate_is_unchanged_by_the_support_box_rule(self):
        """A legacy (contract 2) coupon's process library binds no SupportBox: its box is
        the rows' rule with or without the library, and its estimate is byte-identical
        to the estimate without a process library."""
        self.assertIsNone(read_support_box(None))
        for case_id in ("four-edge-9d2cb9bbb3fe", "ten-edge-6791f1c84123"):
            paths = case_paths(self.manifest, MANIFEST, self.case(case_id))
            self.assertIn("ProcessLibrary", paths)
            self.assertIsNone(read_support_box(paths["ProcessLibrary"]))
            edges = read_edges(paths["Signature"])
            self.assertFalse(any(edge["context"] for edge in edges))
            lower, upper = coupon_box(edges, 2.0, 0.1, 0.05)
            row_lower, row_upper = row_coupon_box(edges, 2.0, 0.1, 0.05)
            self.assertEqual((lower.tolist(), upper.tolist()), (row_lower.tolist(), row_upper.tolist()))
            with_library = estimate(paths, self.options, self.model["TetrahedraPerCubicSize"])
            without = estimate({k: v for k, v in paths.items() if k != "ProcessLibrary"}, self.options,
                               self.model["TetrahedraPerCubicSize"])
            self.assertEqual(json.dumps(with_library, sort_keys=True), json.dumps(without, sort_keys=True))

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

    def test_relabel_calibration_case_is_estimated_at_its_parent_options(self):
        """A label-only Relabel case (decision 56) declares no BuildCommandOptions: its
        preflight estimate is the parent mesh's (ProductionValues alone), under the cap."""
        sizing_path = HERE / "geometry-independence-calibration-sizing.json"
        sizing = json.loads(sizing_path.read_text())
        relabel = next(case for case in sizing["Cases"] if case["Calibration"].get("Relabel") is not None)
        self.assertNotIn("BuildCommandOptions", relabel["Calibration"])
        result = preflight_build_cost(sizing, sizing_path, relabel)
        self.assertTrue(result["Passed"])
        self.assertEqual(result["Options"], relabel["Calibration"]["ProductionValues"])
        self.assertIn("Relabel", result["Origin"]["Options"])
        parent = gate(self.manifest, MANIFEST, self.case(relabel["Calibration"]["BaseCase"]))
        self.assertEqual(result["EstimatedElements"], parent["EstimatedElements"])

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


class FaceEndEstimateTest(unittest.TestCase):
    """Block (b) design A2 / decision 320 in the estimator: a side end on a box face with a
    single metal side and a tilt theta > 0 adds its sheared end block of m layers; a
    Physical-class box vertex is a corner only when its side is exactly perpendicular."""

    def test_oblique_face_ends_add_end_block_layers_and_lose_their_corner(self):
        lower, upper = (-4.0, -4.0), (4.0, 4.0)
        # The O4-like tip: side A leaves x = 4 at 45 degrees (Continuation class at the exit),
        # side B along +y leaves y = 4 perpendicularly (a Physical-class legacy box corner).
        legacy = [[((0.0, 0.0), "Physical", 0.0), ((4.0, 4.0), "Continuation", 0.0),
                   ((0.0, 4.0), "Physical", 0.0)]]
        sides, caps, tilts = metal_sides(legacy, lower, upper, 1e-9)
        self.assertEqual(len(sides), 2)
        self.assertEqual(caps, 4)                       # the tip only: two tubes x two side ends
        self.assertEqual(tilts, [1.0])                   # tan 45 at the A exit; the B exit is exact
        self.assertEqual(semantic_corner_count(legacy), 2)
        # The box side arriving at the oblique exit (Physical class there): a cut end, no corner.
        oblique = [[((0.0, 0.0), "Physical", 0.0), ((4.0, 0.0), "Continuation", 0.0),
                    ((4.0, 4.0), "Physical", 0.0)]]
        sides, caps, tilts = metal_sides(oblique, lower, upper, 1e-9)
        self.assertEqual(len(sides), 2)
        self.assertEqual(tilts, [1.0])
        self.assertEqual(semantic_corner_count(oblique), 1)
        # Two metal sides meeting at a box vertex: a corner, no face end (an island touching x = 4).
        island = [[((4.0, 0.0), "Physical", 0.0), ((2.0, 1.0), "Physical", 0.0), ((2.0, -1.0), "Physical", 0.0)]]
        sides, caps, tilts = metal_sides(island, lower, upper, 1e-9)
        self.assertEqual((len(sides), tilts), (3, []))
        self.assertEqual(semantic_corner_count(island), 3)
        # The end-block layer counts of the production fabricated tube (R 31.75 nm, h_K 16 nm,
        # TangentialSize 50 nm): 45 degrees -> 2 layers at lc_end 50 nm; 74.3 degrees -> 3 at
        # lc_end 114 nm; 12.6 degrees -> 1 (design A2 (4)); thin (R 62, h_K 32) at 45 -> 3.
        ring_sizes, radius = tube_rings(0.00025, 2.0, 0.05, 0.1, 0.1)
        self.assertEqual(face_end_layers([1.0], ring_sizes, radius, 0.05), 2)
        self.assertEqual(face_end_layers([math.tan(math.radians(74.3))], ring_sizes, radius, 0.05), 3)
        self.assertEqual(face_end_layers([math.tan(math.radians(12.6))], ring_sizes, radius, 0.05), 1)
        self.assertEqual(face_end_layers([1.0, 1.0], ring_sizes, radius, 0.05), 4)
        thin_sizes, thin_radius = tube_rings(0.002, 2.0, 0.05, 0.1, 0.1, thin=True)
        self.assertEqual(face_end_layers([1.0], thin_sizes, thin_radius, 0.05), 3)
        self.assertEqual(face_end_layers([], ring_sizes, radius, 0.05), 0)
