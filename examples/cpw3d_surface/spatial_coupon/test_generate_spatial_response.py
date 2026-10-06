# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""generate_spatial_response: the conductor labelling of the trace knots against the plan-view
mask (decision 360 (c*): a tolerance of 4 x the 1e-9 R plan quantum, inclusive, fail-closed to
PEC), the labelling validator (no FREE knot within the tolerance of a metal outline; the near-
outline inspection record) and the mask-frame check (design round 2 review MINOR-5)."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cluster_signature_geometry as csg  # noqa: E402
import generate_spatial_response as gsr  # noqa: E402

RADIUS = 1.9
THICKNESS = 0.1
# Conductor 1 covers x in [-1, 0], y in [-1, 1] on the plane z = 0 (two mask facets).
FACETS = [{"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, -1.0], [0.0, 1.0]]},
          {"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, 1.0], [-1.0, 1.0]]}]
EDGES = [{"Conductor": 1, "Point": [-1.0, 0.0, 0.0], "ProcessNormal": [0.0, 0.0, 1.0], "Tangent": [0.0, 1.0, 0.0],
          "GapDirection": [-1.0, 0.0, 0.0], "Interval": [-1.0, 1.0], "VertexArm": False}]
QUANTUM = 1.0e-9 * RADIUS


def labels_of(points):
    return gsr.conductor_at_points(np.asarray(points, dtype=float), EDGES, RADIUS, THICKNESS, 90.0, FACETS)


class ConductorLabelToleranceTest(unittest.TestCase):
    def test_tolerance_is_four_plan_quanta_inclusive(self):
        self.assertEqual(gsr.CONDUCTOR_LABEL_TOLERANCE_OVER_R, 4.0e-9)
        inside = [-0.5, 0.0, 0.0]
        on_outline = [-1.0, 0.5, 0.0]
        rounding = [-1.0 - 0.3 * QUANTUM, -1.0 - 0.3 * QUANTUM, THICKNESS]  # a corner knot ~0.4 quanta outside (M2)
        three_quanta = [-1.0 - 3.0 * QUANTUM, 0.5, 0.0]
        four_quanta = [-1.0 - 4.0 * QUANTUM, 0.5, THICKNESS]                 # exactly the tolerance: inclusive
        five_quanta = [-1.0 - 5.0 * QUANTUM, 0.5, 0.0]
        above_metal = [-0.5, 0.0, 0.5]
        labels = labels_of([inside, on_outline, rounding, three_quanta, four_quanta, five_quanta, above_metal])
        self.assertEqual(labels.tolist(), [1, 1, 1, 1, 1, 0, 0])

    def test_edge_only_signature_without_facets_is_unchanged(self):
        # The edge-only signature puts the metal on the side opposite the gap direction (x >= -1 here).
        points = np.array([[-1.0 - 3.0 * QUANTUM, 0.5, 0.0], [-0.5, 0.0, 0.0], [-1.0 - 0.5 * 1.0e-10 * RADIUS, 0.0, 0.0]])
        labels = gsr.conductor_at_points(points, EDGES, RADIUS, THICKNESS, 90.0, [])
        self.assertEqual(labels.tolist(), [0, 1, 1])  # tolerance 1e-10 R, no plan-view quantum involved


class NearOutlineFreeKnotsTest(unittest.TestCase):
    def test_validator_binds_the_labels_and_records_the_inspection_band(self):
        points = np.array([[-1.0 - 3.0 * QUANTUM, 0.5, 0.0],        # within the tolerance: must be labelled
                           [-1.0 - 5.0 * QUANTUM, 0.5, THICKNESS],  # free, inside the 1e-3 R inspection band
                           [-1.0 - 2.0e-3 * RADIUS, 0.5, 0.0],      # free, beyond the band
                           [-1.0 - 3.0 * QUANTUM, 0.5, 0.5],        # above the metal: not judged
                           [-0.5, 0.0, 0.0]])                       # in the metal
        labels = labels_of(points)
        self.assertEqual(labels.tolist(), [1, 0, 0, 0, 1])
        violations, inspect = gsr.near_outline_free_knots(points, labels, EDGES, RADIUS, THICKNESS, FACETS)
        self.assertEqual(violations, [])
        self.assertEqual([k["Vertex"] for k in inspect], [2])
        self.assertAlmostEqual(inspect[0]["Distance"], 5.0 * QUANTUM, places=15)
        # A label set that leaves a knot within the tolerance free is a violation (fail closed).
        stale = labels.copy()
        stale[0] = 0
        violations, inspect = gsr.near_outline_free_knots(points, stale, EDGES, RADIUS, THICKNESS, FACETS)
        self.assertEqual([k["Vertex"] for k in violations], [1])
        self.assertEqual([k["Vertex"] for k in inspect], [2])
        self.assertEqual(gsr.near_outline_free_knots(points, stale, EDGES, RADIUS, THICKNESS, []), ([], []))


class MaskFrameTest(unittest.TestCase):
    def test_mask_inside_the_box_passes_and_a_foreign_frame_fails_closed(self):
        lower, upper = np.array([-1.0, -1.0, -2.0]), np.array([3.0, 1.0, 2.0])
        self.assertEqual(gsr.validate_mask_frame(FACETS, lower, upper, RADIUS), [-1.0, -1.0, 0.0, 1.0])
        self.assertIsNone(gsr.validate_mask_frame([], lower, upper, RADIUS))
        slack = gsr.CONDUCTOR_LABEL_TOLERANCE_OVER_R * RADIUS
        self.assertIsNotNone(gsr.validate_mask_frame(FACETS, lower + [0.5 * slack, 0.0, 0.0], upper, RADIUS))
        rotated = [{"Conductor": 1, "Plane": 0.0, "Points": [[2.0, 3.0], [12.0, 3.0], [12.0, 9.0]]}]
        with self.assertRaises(ValueError) as caught:
            gsr.validate_mask_frame(rotated, lower, upper, RADIUS)
        self.assertIn("not in the mesh frame", str(caught.exception))
        with self.assertRaises(ValueError):
            gsr.validate_mask_frame(FACETS, lower + [2.0 * slack, 0.0, 0.0], upper, RADIUS)

    def test_generator_stops_before_any_output_on_a_foreign_frame(self):
        """Decision 369 review MINOR-2: the frame check runs before the first output file, so a
        mis-framed mask leaves generation-failure.json (Stage PlanViewMaskFrame) and nothing
        else: no partial mesh-signature.csv next to it."""
        law = '{"Type":"PEC"}'
        portions = [{"Conductor": 1, "Gap": [0.0, 1.0], "Interfaces": ["MA", "MS", "SA"], "Law": law, "P": [-2.0, 0.0, 0.0, 0.0]},
                    {"Conductor": 1, "Gap": [-1.0, 0.0], "Interfaces": ["MA", "MS", "SA"], "Law": law, "P": [0.0, 0.0, 0.0, 2.0]}]
        signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 2, "Portions": portions,
                     "Vertices": [{"P": [0.0, 0.0], "TurnDegrees": 90.0, "Type": "ConcaveCorner"}]}
        cluster = {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 2, "Signature": signature}, "Signature": signature,
                   "Interfaces": [{"Slot": 0, "Type": t, "Target": n} for n, t in ((3, "MA"), (2, "MS"), (1, "SA"))],
                   "BoundaryCondition": {"Type": "PEC"}}
        coupon, _ = csg.cluster_coupon(cluster, RADIUS, THICKNESS, 0.05)
        for facet in coupon["Geometry"]["PlanViewFacets"]:  # the mask of another frame: shifted by 1 um
            facet["Points"] = [[p[0] + 1.0, p[1]] + list(p[2:]) for p in facet["Points"]]
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "coupon.json").write_text(json.dumps(coupon))
            result = subprocess.run([sys.executable, str(HERE / "generate_spatial_response.py"), str(Path(tmp) / "coupon.json"),
                                     "--output", str(Path(tmp) / "out"), "--radius", str(RADIUS), "--metal-thickness",
                                     str(THICKNESS), "--overetch-depth", "0.05", "--sidewall-angle", "90", "--top-rounding", "0",
                                     "--trench-rounding", "0", "--model-name", "m", "--basis-only"],
                                    capture_output=True, text=True, timeout=120)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not in the mesh frame", result.stderr)
            failure = json.loads((Path(tmp) / "out" / "generation-failure.json").read_text())
            self.assertEqual(failure["Stage"], "PlanViewMaskFrame")
            self.assertEqual(sorted(p.name for p in (Path(tmp) / "out").iterdir()), ["generation-failure.json"])


if __name__ == "__main__":
    unittest.main()
