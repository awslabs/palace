#!/usr/bin/env python3

import json
import sys
import unittest
from pathlib import Path

CPW2D = Path(__file__).parent
sys.path.insert(0, str(CPW2D))
import discover_surface_response_requirements as DISCOVERY  # noqa: E402


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
        _, model = DISCOVERY.placeholder_model(single, 2.0)
        self.assertIn("Reference", model)
        self.assertNotIn("ConductorReferences", model)

        different = self.requirement(
            "DifferentConductorGap", {"EdgeCount": 2, "Separation": 1.0}
        )
        _, model = DISCOVERY.placeholder_model(different, 2.0)
        self.assertEqual(len(model["ConductorReferences"]), 2)
        self.assertNotIn("Reference", model)

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
