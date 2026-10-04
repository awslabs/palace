#!/usr/bin/env python3

import json
import os
import stat
import sys
import tempfile
import textwrap
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

    # Two version-2 Missing requirements of a fake device: an IsolatedEdge (Hash aaaa...) and a
    # 3-edge ParallelEdgeCluster (Hash bbbb...); the fake Palace below reads every pass's
    # library and reports a requirement Exact once a placeholder named after its Hash exists.
    FAKE_REQUIREMENTS = [
        {"Topology": "IsolatedEdge", "Hash": "a" * 64, "Count": 4, "Instances": 2, "TotalEdgeLength": 10.0,
         "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
         "Signature": {"Type": "IsolatedEdge", "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}"}},
        {"Topology": "ParallelEdgeCluster", "Hash": "b" * 64, "Count": 7, "Instances": 2, "TotalEdgeLength": 8.2,
         "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
         "Geometry": {"EdgeCount": 3},
         "Signature": {"Type": "ParallelEdgeCluster", "Edges": [
             {"Conductor": 1, "GapSide": -1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 0.0},
             {"Conductor": 1, "GapSide": 1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 1.052632},
             {"Conductor": 1, "GapSide": -1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 2.105263}]}},
    ]

    def fake_device(self, tmp):
        """A device config + source library + an executable fake Palace (geometry preflight only)."""
        requirements = tmp / "fake-requirements.json"
        requirements.write_text(json.dumps(self.FAKE_REQUIREMENTS))
        library = tmp / "process-library.json"
        library.write_text(json.dumps({"Version": 3, "Name": "fake-seed", "MatchingRadius": 1.9, "Models": []}))
        config = tmp / "device.json"
        config.write_text(json.dumps({
            "Problem": {"Type": "Electrostatic", "Output": str(tmp / "postpro")},
            "Solver": {"Electrostatic": {"ResponseCorrection": {"Library": str(library), "UnmatchedPolicy": "Warn"}}}}))
        palace = tmp / "fake-palace.py"
        palace.write_text(textwrap.dedent(f"""\
            #!{sys.executable}
            import json, pathlib, sys
            assert sys.argv[1] == "--surface-response-preflight"
            config = json.load(open(sys.argv[2]))
            response = config["Solver"]["Electrostatic"]["ResponseCorrection"]
            library = json.load(open(response["Library"]))
            names = {{model["Name"] for model in library["Models"]}}
            requirements = json.load(open({str(requirements)!r}))
            counts = {{"Exact": 0, "Interpolated": 0, "Missing": 0}}
            lengths = {{"Exact": 0.0, "Interpolated": 0.0, "Missing": 0.0}}
            for requirement in requirements:
                placeholder = "__preflight_placeholder_" + requirement["Hash"][:16]
                if placeholder in names:
                    requirement["Status"] = "Exact"
                    requirement["SelectedModels"] = [{{"Name": placeholder, "Weight": 1.0}}]
                else:
                    requirement["Status"] = "Missing"
                counts[requirement["Status"]] += requirement["Count"]
                lengths[requirement["Status"]] += requirement["TotalEdgeLength"]
            output = pathlib.Path(config["Problem"]["Output"]); output.mkdir(parents=True, exist_ok=True)
            manifest = {{"Version": 2, "Complete": counts["Missing"] == 0, "Requirements": requirements,
                        "Library": {{"Name": library["Name"], "Path": response["Library"],
                                    "MatchingRadius": library["MatchingRadius"]}},
                        "Summary": {{"Counts": counts, "TotalEdgeLengths": lengths}}}}
            (output / "surface-response-requirements.json").write_text(json.dumps(manifest))
            """))
        palace.chmod(palace.stat().st_mode | stat.S_IXUSR)
        return config, palace

    def test_omitted_requirement_gets_no_placeholder_and_stays_missing(self):
        """`--omit-requirement HASH_PREFIX`: the matching Missing requirement gets no placeholder in
        any pass, does not stall the closure, stays Missing in the final manifest (counted in Summary.Missing, Complete false) and is listed under
        OmittedRequirements of the manifest and the closure history; the other requirement
        closes as before. A prefix matching nothing fails closed."""
        with tempfile.TemporaryDirectory(prefix="discover-omit-") as tmp:
            tmp = Path(tmp)
            config, palace = self.fake_device(tmp)
            # Without an omission the fake device closes in two passes with two placeholders.
            manifest = DISCOVERY.discover(config, tmp / "full", palace)
            self.assertEqual(manifest["Summary"]["Counts"], {"Exact": 0, "Interpolated": 0, "Missing": 11})
            self.assertEqual(manifest["OmittedRequirements"], [])
            history = json.loads((tmp / "full" / "closure-history.json").read_text())
            self.assertEqual(history["PlaceholderCount"], 2)
            self.assertEqual(len(history["Passes"]), 2)
            self.assertEqual(sorted(p["Topology"] for p in history["Passes"][0]["AddedPlaceholders"]),
                             ["IsolatedEdge", "ParallelEdgeCluster"])
            # The 3-edge cluster omitted: one placeholder, the omitted requirement Missing throughout.
            manifest = DISCOVERY.discover(config, tmp / "omit", palace, omit_requirements=["bbbbbbbbbbbb"])
            self.assertFalse(manifest["Complete"])
            self.assertEqual(manifest["Summary"]["Counts"], {"Exact": 0, "Interpolated": 0, "Missing": 11})
            by_hash = {requirement["Hash"]: requirement for requirement in manifest["Requirements"]}
            self.assertEqual(by_hash["a" * 64]["Status"], "Missing")   # the restored placeholder
            self.assertEqual(by_hash["a" * 64]["Reason"], "Missing from source library after exhaustive geometry discovery")
            self.assertEqual(by_hash["b" * 64]["Status"], "Missing")   # never had a placeholder
            self.assertNotIn("Reason", by_hash["b" * 64])
            self.assertEqual(manifest["OmittedRequirements"],
                             [{"Hash": "b" * 64, "Prefix": "bbbbbbbbbbbb", "Topology": "ParallelEdgeCluster", "Count": 7,
                               "Instances": 2, "TotalEdgeLength": 8.2}])
            history = json.loads((tmp / "omit" / "closure-history.json").read_text())
            self.assertEqual(history["PlaceholderCount"], 1)
            self.assertEqual(history["OmittedRequirements"], manifest["OmittedRequirements"])
            self.assertEqual(len(history["Passes"]), 2)
            self.assertEqual([p["Topology"] for p in history["Passes"][0]["AddedPlaceholders"]], ["IsolatedEdge"])
            self.assertEqual(history["Passes"][1]["AddedPlaceholders"], [])
            for pass_dir in ("pass-01", "pass-02"):
                library = json.loads((tmp / "omit" / pass_dir / "process-library.json").read_text())
                self.assertNotIn("__preflight_placeholder_" + "b" * 16, {m["Name"] for m in library["Models"]})
            # A prefix that matches no Missing requirement is a mistake, not a no-op.
            with self.assertRaisesRegex(RuntimeError, "matching no Missing requirement.*cccccccc"):
                DISCOVERY.discover(config, tmp / "typo", palace, omit_requirements=["cccccccc"])
            with self.assertRaises(ValueError):
                DISCOVERY.discover(config, tmp / "empty", palace, omit_requirements=[""])

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
