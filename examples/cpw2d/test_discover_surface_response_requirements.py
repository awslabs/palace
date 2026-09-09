#!/usr/bin/env python3

import json
import sys
import unittest
from pathlib import Path

CPW2D = Path(__file__).parent
SPATIAL = CPW2D.parent / "cpw3d_surface" / "spatial_coupon"
sys.path.insert(0, str(CPW2D))
sys.path.insert(0, str(SPATIAL))
import discover_surface_response_requirements as DISCOVERY  # noqa: E402
import generate_spatial_response as SPATIAL_RESPONSE  # noqa: E402


class DiscoverSurfaceResponseRequirementsTest(unittest.TestCase):
    def requirement(self, topology, geometry, edges=None):
        if edges is not None:
            geometry = {**geometry, "Edges": edges}
        return {
            "Topology": topology,
            "Geometry": geometry,
            "BoundaryCondition": {"Type": "PEC"},
            "Interfaces": [
                {"Slot": 0, "Target": 1, "Type": "SA"},
                {"Slot": 0, "Target": 2, "Type": "MS"},
            ],
        }

    def test_single_and_multiple_conductor_placeholders(self):
        single = self.requirement(
            "SpatialEdgeCluster",
            {"EdgeCount": 2},
            [
                {
                    "Conductor": 1,
                    "Point": [0.0, 0.0, 0.0],
                    "GapDirection": [1.0, 0.0, 0.0],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-1.0, 1.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
                {
                    "Conductor": 1,
                    "Point": [1.0, 0.0, 0.0],
                    "GapDirection": [-1.0, 0.0, 0.0],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-1.0, 1.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
            ],
        )
        single["Geometry"]["PlanViewBoundary"] = [
            {
                "Conductor": 1,
                "Segments": [[[0, 0, 0], [1, 0, 0]]],
                "ContinuationSegments": [],
            }
        ]
        single["Geometry"]["MaskRegularization"] = {"Version": 1}
        _, model = DISCOVERY.placeholder_model(
            single,
            2.0,
            {"MetalThickness": 0.1, "OveretchDepth": 0.05},
        )
        self.assertIn("Reference", model)
        self.assertNotIn("ConductorReferences", model)
        self.assertEqual(
            model["PlanViewBoundary"], single["Geometry"]["PlanViewBoundary"]
        )
        self.assertEqual(
            model["MaskRegularization"], single["Geometry"]["MaskRegularization"]
        )
        self.assertEqual(len(model["SupportPoints"]), 8)

        different = self.requirement(
            "DifferentConductorGap", {"EdgeCount": 2, "Separation": 1.0}
        )
        _, model = DISCOVERY.placeholder_model(different, 2.0)
        self.assertEqual(len(model["ConductorReferences"]), 2)
        self.assertNotIn("Reference", model)

    def test_virtual_and_production_spatial_support_are_identical(self):
        direction = 2.0**-0.5
        geometry = {
            "EdgeCount": 2,
            "Edges": [
                {
                    "Conductor": 1,
                    "Point": [0.123456789, 0.0, -0.987654321],
                    "GapDirection": [direction, 0.0, direction],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-2.0, 0.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
                {
                    "Conductor": 1,
                    "Point": [1.234567891, 0.0, 0.345678912],
                    "GapDirection": [-direction, 0.0, -direction],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [0.0, 2.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
            ],
        }
        fabrication = {"MetalThickness": 0.1, "OveretchDepth": 0.05}
        virtual = DISCOVERY.spatial_support_points(geometry, 2.0, fabrication)
        coupon = {
            "Topology": "SpatialEdgeCluster",
            "Geometry": geometry,
            "BoundaryCondition": {"Type": "PEC"},
        }
        frame, edges, _ = SPATIAL_RESPONSE.normalize_geometry(coupon, 2.0)
        lower, upper = SPATIAL_RESPONSE.coupon_bounds(edges, 2.0, 0.1, 0.05)
        production = SPATIAL_RESPONSE.matching_support_points(
            lower, upper, frame, 2.0
        )
        self.assertEqual(virtual, production)

    def test_restore_source_status_marks_only_placeholders_missing(self):
        manifest = {
            "Complete": True,
            "Library": {"Name": "closure", "Path": "/tmp/closure.json"},
            "Summary": {},
            "Requirements": [
                {
                    "Status": "Exact",
                    "Count": 2,
                    "TotalEdgeLength": 3.0,
                    "SelectedModels": [{"Name": "real", "Weight": 1.0}],
                },
                {
                    "Status": "Exact",
                    "Count": 4,
                    "TotalEdgeLength": 5.0,
                    "SelectedModels": [
                        {"Name": "__preflight_placeholder_deadbeef", "Weight": 1.0}
                    ],
                    "NormalizedLibraryDistance": 0.0,
                },
            ],
        }
        source = {"Name": "source", "__SourcePath": Path("/tmp/source.json")}
        result = DISCOVERY.restore_source_status(
            manifest,
            source,
            {
                "__preflight_placeholder_deadbeef": {
                    "Topology": "IsolatedEdge",
                    "Geometry": {"EdgeCount": 1},
                    "BoundaryCondition": {"Type": "PEC"},
                    "Interfaces": [],
                }
            },
        )
        self.assertFalse(result["Complete"])
        self.assertEqual(
            result["Summary"],
            {
                "Counts": {"Exact": 2, "Interpolated": 0, "Missing": 4},
                "TotalEdgeLengths": {
                    "Exact": 3.0,
                    "Interpolated": 0.0,
                    "Missing": 5.0,
                },
            },
        )
        missing = result["Requirements"][1]
        self.assertEqual(missing["Status"], "Missing")
        self.assertNotIn("SelectedModels", missing)
        self.assertNotIn("NormalizedLibraryDistance", missing)


if __name__ == "__main__":
    unittest.main()
