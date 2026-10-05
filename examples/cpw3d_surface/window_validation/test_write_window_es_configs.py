#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the window electrostatic configuration writer (`write_window_es_configs.py`).

    python3 -m unittest test_write_window_es_configs
"""

import copy
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from write_window_es_configs import LAYERS, fabricated_config, main, terminal_labels, thin_config  # noqa: E402

SQUARE = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
POLYGON_SET = {
    "Version": 1,
    "Name": "W",
    "Box": {"X": [0.0, 10.0], "Y": [0.0, 10.0]},
    "Planes": [
        {"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "SubstrateThickness": 525.0,
         "Polygons": [{"Conductor": "ground", "Outer": SQUARE, "Holes": [[[2.0, 2.0], [8.0, 2.0], [8.0, 8.0], [2.0, 8.0]]]},
                      {"Conductor": "trace_7", "Outer": [[3.0, 3.0], [7.0, 3.0], [7.0, 7.0], [3.0, 7.0]]}]},
        {"Name": "L2", "SurfaceZ": 4.8, "Facing": "down", "SubstrateThickness": 525.0,
         "Polygons": [{"Conductor": "ground", "Outer": SQUARE, "Holes": [[[1.0, 1.0], [9.0, 1.0], [9.0, 9.0], [1.0, 9.0]]]},
                      {"Conductor": "island_3", "Outer": [[2.0, 2.0], [8.0, 2.0], [8.0, 8.0], [2.0, 8.0]]}]},
    ],
    "Bumps": [{"Conductor": "ground", "Footprint": [[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8]]}],
    "Vacuum": {"Below": 0.0, "Above": 0.0},
    "Terminals": ["island_3", "trace_7"],
}

# The mesher's manifest (`attributes` lists substrate_backside even when no face carries it).
FABRICATED_MANIFEST = {
    "attributes": {
        "substrate": {"attribute": 1, "dimension": 3}, "vacuum": {"attribute": 2, "dimension": 3},
        "exterior_boundary": {"attribute": 3, "dimension": 2}, "ground_air": {"attribute": 4, "dimension": 2},
        "ground_substrate": {"attribute": 5, "dimension": 2}, "substrate_air": {"attribute": 6, "dimension": 2},
        "island_3_air": {"attribute": 7, "dimension": 2}, "island_3_substrate": {"attribute": 8, "dimension": 2},
        "substrate_backside": {"attribute": 9, "dimension": 2},
        "trace_7_air": {"attribute": 10, "dimension": 2}, "trace_7_substrate": {"attribute": 11, "dimension": 2},
    },
    "surface_attribute_counts": {"3": 100, "4": 50, "5": 40, "6": 60, "7": 20, "8": 10, "10": 12, "11": 6},
}

THIN_MANIFEST = {
    "Attributes": {
        "bump_surface": 153, "exterior_boundary": 154, "gap_L1": 145, "gap_L2": 28, "ground_L1": 114, "ground_L2": 126,
        "island_3_L2": 1002, "trace_7_L1": 1001, "substrate_l1": 2, "substrate_l2": 1, "vacuum": 3,
    },
    "Planes": [{"Name": "L1", "Facing": "up", "SurfaceZ": 0.0}, {"Name": "L2", "Facing": "down", "SurfaceZ": 4.8}],
    "Terminals": ["island_3", "trace_7"],
}


def by_type(config, kind):
    return [d for d in config["Boundaries"]["Postprocessing"]["Dielectric"] if d["Type"] == kind]


class TerminalLabelsTest(unittest.TestCase):
    def test_order_from_terminals(self):
        self.assertEqual(terminal_labels(POLYGON_SET), ["island_3", "trace_7"])

    def test_default_sorted(self):
        polygon_set = copy.deepcopy(POLYGON_SET)
        del polygon_set["Terminals"]
        self.assertEqual(terminal_labels(polygon_set), ["island_3", "trace_7"])

    def test_incomplete_terminals_refused(self):
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Terminals"] = ["island_3"]
        with self.assertRaises(ValueError):
            terminal_labels(polygon_set)


class FabricatedConfigTest(unittest.TestCase):
    def test_transmon_conventions(self):
        config = fabricated_config(POLYGON_SET, FABRICATED_MANIFEST, "/m/W.msh2", "/m/postpro")
        self.assertEqual(config["Problem"]["Type"], "Electrostatic")
        self.assertEqual(config["Model"], {"Mesh": "/m/W.msh2", "L0": 1.0e-6, "CrackInternalBoundaryElements": False, "Refinement": {"MaxIts": 0}})
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [4, 5])
        self.assertEqual(config["Boundaries"]["Terminal"], [{"Index": 1, "Attributes": [7, 8]}, {"Index": 2, "Attributes": [10, 11]}])
        sa, ms, ma = by_type(config, "SA"), by_type(config, "MS"), by_type(config, "MA")
        # No attribute 9 with Vacuum 0 (review d198 m5): the SA integral is substrate_air only.
        self.assertEqual([d["Attributes"] for d in sa], [[6]])
        self.assertEqual([d["Attributes"] for d in ms], [[5, 8, 11]])
        self.assertEqual([d["Attributes"] for d in ma], [[4, 7, 10]])
        for d in sa + ms + ma:
            self.assertEqual(d["Thickness"], 0.002)
            self.assertNotIn("AutomaticEdges", d)
        self.assertEqual((sa[0]["Permittivity"], ms[0]["Permittivity"], ma[0]["Permittivity"]), (4.0, 11.45, 10.0))
        self.assertEqual([d["Index"] for d in config["Boundaries"]["Postprocessing"]["Dielectric"]], [1, 2, 3])
        self.assertEqual(config["Solver"], {"Order": 4, "Device": "CPU", "Electrostatic": {"Save": 0}, "Linear": {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1.0e-10, "MaxIts": 5000}})
        self.assertEqual(config["Domains"]["Materials"][0], {"Attributes": [1], "Permittivity": 11.45})
        self.assertEqual(config["Domains"]["Postprocessing"]["Energy"], [{"Index": 1, "Attributes": [1]}, {"Index": 2, "Attributes": [2]}])

    def test_order_and_cheap_estimator(self):
        config = fabricated_config(POLYGON_SET, FABRICATED_MANIFEST, "m", "p", order=5, estimator_cheap=True)
        self.assertEqual(config["Solver"]["Order"], 5)
        self.assertEqual(config["Solver"]["Linear"]["EstimatorTol"], 0.1)
        self.assertEqual(config["Solver"]["Linear"]["EstimatorMaxIts"], 100)
        self.assertEqual(config["Model"]["Refinement"], {"MaxIts": 0})

    def test_backside_sa_only_with_vacuum(self):
        # The transmon box (Vacuum > 0): substrate_backside carries faces -> a separate SA index 4.
        manifest = copy.deepcopy(FABRICATED_MANIFEST)
        manifest["surface_attribute_counts"]["9"] = 30
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Vacuum"] = {"Below": 475.0, "Above": 1000.0}
        config = fabricated_config(polygon_set, manifest, "m", "p")
        self.assertEqual([(d["Index"], d["Attributes"]) for d in by_type(config, "SA")], [(1, [6]), (4, [9])])
        # Faces on 9 with Vacuum 0, or Vacuum > 0 without faces on 9: inconsistent input, refused.
        with self.assertRaises(ValueError):
            fabricated_config(POLYGON_SET, manifest, "m", "p")
        with self.assertRaises(ValueError):
            fabricated_config(polygon_set, FABRICATED_MANIFEST, "m", "p")

    def test_single_up_plane_vacuum_above_has_no_backside(self):
        # One `up` plane with Vacuum.Above > 0 and Below 0 (the OSC windows): the vacuum lies over
        # the metal, the backside is the box wall -> no substrate_backside faces, accepted
        # (graded-sweep-mesher M1.5 item 5); Below > 0 on the same plane needs faces on 9.
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Planes"] = [polygon_set["Planes"][0]]
        polygon_set["Terminals"] = ["trace_7"]
        polygon_set["Vacuum"] = {"Below": 0.0, "Above": 1000.0}
        manifest = copy.deepcopy(FABRICATED_MANIFEST)
        del manifest["attributes"]["island_3_air"], manifest["attributes"]["island_3_substrate"]
        config = fabricated_config(polygon_set, manifest, "m", "p")
        self.assertEqual([(d["Index"], d["Attributes"]) for d in by_type(config, "SA")], [(1, [6])])
        polygon_set["Vacuum"] = {"Below": 475.0, "Above": 1000.0}
        with self.assertRaises(ValueError):
            fabricated_config(polygon_set, manifest, "m", "p")
        manifest["surface_attribute_counts"]["9"] = 30
        config = fabricated_config(polygon_set, manifest, "m", "p")
        self.assertEqual([(d["Index"], d["Attributes"]) for d in by_type(config, "SA")], [(1, [6]), (4, [9])])

    def test_missing_or_empty_terminal_attribute_refused(self):
        manifest = copy.deepcopy(FABRICATED_MANIFEST)
        del manifest["attributes"]["trace_7_air"]
        with self.assertRaises(ValueError):
            fabricated_config(POLYGON_SET, manifest, "m", "p")
        manifest = copy.deepcopy(FABRICATED_MANIFEST)
        manifest["surface_attribute_counts"]["11"] = 0
        with self.assertRaises(ValueError):
            fabricated_config(POLYGON_SET, manifest, "m", "p")


class ThinConfigTest(unittest.TestCase):
    def test_path_T(self):
        config = thin_config(POLYGON_SET, THIN_MANIFEST, "/m/thin.msh2", "/lib/process-library.json", "T", "/m/postpro")
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [114, 126, 153])
        self.assertEqual(config["Boundaries"]["Terminal"], [{"Index": 1, "Attributes": [1002]}, {"Index": 2, "Attributes": [1001]}])
        dielectrics = config["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertEqual([(d["Index"], d["Type"], d["Attributes"]) for d in dielectrics],
                         [(1, "SA", [28, 145]), (2, "MS", [114, 1001]), (3, "MA", [114, 1001]), (4, "MS", [126, 1002]), (5, "MA", [126, 1002])])
        self.assertEqual([d.get("EdgeFrameNormal") for d in dielectrics], [None, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])
        for d in dielectrics:
            self.assertTrue(d["AutomaticEdges"])
            self.assertEqual(d["EdgeDistances"], [1.9])
        self.assertEqual(config["Domains"]["Materials"][0], {"Attributes": [1, 2], "Permittivity": 11.45})
        self.assertEqual(config["Domains"]["Materials"][1], {"Attributes": [3], "Permittivity": 1.0})
        self.assertEqual(config["Model"]["Refinement"], {"MaxIts": 11, "UniformLevels": 0, "SerialUniformLevels": 0, "Tol": 0.01, "UpdateFraction": 0.7, "MaximumImbalance": 1.15,
                                                         "MaxNCLevels": 8, "Nonconformal": True, "SaveAdaptIterations": True, "SaveAdaptMesh": False, "MaxSize": 400000000})
        self.assertEqual(config["Solver"]["Order"], 5)
        response = config["Solver"]["Electrostatic"]["ResponseCorrection"]
        self.assertEqual(response["Library"], "/lib/process-library.json")
        self.assertEqual(response["TargetInterfaces"], [1, 2, 3, 4, 5])
        self.assertEqual((response["TraceCoupling"], response["MortarOversampling"], response["SolveTol"], response["CorrectionMode"], response["TranslationalDomainCorrection"], response["PatchConstruction"]),
                         ("SurfaceMortar", 2, 1.0e-6, "Both", "FixedTrace", "Features"))
        self.assertEqual(config["Solver"]["Electrostatic"]["Save"], 1)
        self.assertEqual(config["Solver"]["Linear"], {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1.0e-10, "MaxIts": 1000, "EstimatorTol": 0.1, "EstimatorMG": True})

    def test_path_P4_and_preflight(self):
        p4 = thin_config(POLYGON_SET, THIN_MANIFEST, "m", "lib", "P4", "p")
        self.assertEqual((p4["Solver"]["Order"], p4["Model"]["Refinement"]["MaxIts"], p4["Model"]["Refinement"]["MaxSize"]), (4, 15, 100000000))
        preflight = thin_config(POLYGON_SET, THIN_MANIFEST, "m", "lib", "preflight", "p", radius=2.1)
        self.assertEqual(preflight["Solver"]["Order"], 2)
        self.assertEqual(preflight["Model"]["Refinement"], {"MaxIts": 0, "UniformLevels": 0, "SerialUniformLevels": 0})
        self.assertEqual(preflight["Solver"]["Electrostatic"]["ResponseCorrection"]["TraceCoupling"], "Collocated")
        self.assertEqual(preflight["Solver"]["Electrostatic"]["Save"], 0)
        self.assertEqual(preflight["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"], [2.1])
        with self.assertRaises(ValueError):
            thin_config(POLYGON_SET, THIN_MANIFEST, "m", "lib", "p3", "p")

    def test_plane_without_gap_keeps_metal_interfaces(self):
        # Decisions 329 / 341: a plane that is one ground sheet over the whole window (S5 / S6's
        # L2) has no gap sheet -> no SA on that plane, but its metal sheet keeps MS / MA as PLAIN
        # interfaces (no AutomaticEdges / EdgeDistances / EdgeFrameNormal, the fabricated
        # writer's form) outside TargetInterfaces: a response target needs a metal perimeter the
        # edgeless sheet has not, and its raw energy is exact. The earlier rule dropped all three.
        manifest = copy.deepcopy(THIN_MANIFEST)
        del manifest["Attributes"]["gap_L2"]
        del manifest["Attributes"]["island_3_L2"]
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Planes"][1]["Polygons"] = [{"Conductor": "ground", "Outer": SQUARE}]
        polygon_set["Terminals"] = ["trace_7"]
        config = thin_config(polygon_set, manifest, "m", "lib", "T", "p")
        dielectrics = config["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertEqual([(d["Index"], d["Type"], d["Attributes"]) for d in dielectrics],
                         [(1, "SA", [145]), (2, "MS", [114, 1001]), (3, "MA", [114, 1001]), (4, "MS", [126]), (5, "MA", [126])])
        self.assertEqual([d.get("EdgeFrameNormal") for d in dielectrics], [None, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], None, None])
        self.assertEqual([d.get("AutomaticEdges", False) for d in dielectrics], [True, True, True, False, False])
        self.assertEqual(["EdgeDistances" in d for d in dielectrics], [True, True, True, False, False])
        for d in dielectrics[3:]:
            self.assertEqual((d["Thickness"], d["Permittivity"], d["LossTan"]), (dielectrics[1]["Thickness"], LAYERS[d["Type"]]["Permittivity"], LAYERS[d["Type"]]["LossTan"]))
        self.assertEqual(config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"], [1, 2, 3])
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [114, 126, 153])
        # The same on the other plane (an `up` gapless L1 under an L2 with edges).
        manifest = copy.deepcopy(THIN_MANIFEST)
        del manifest["Attributes"]["gap_L1"]
        del manifest["Attributes"]["trace_7_L1"]
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Planes"][0]["Polygons"] = [{"Conductor": "ground", "Outer": SQUARE}]
        polygon_set["Terminals"] = ["island_3"]
        config = thin_config(polygon_set, manifest, "m", "lib", "T", "p")
        self.assertEqual([(d["Type"], d["Attributes"], d.get("EdgeFrameNormal"), d.get("AutomaticEdges", False)) for d in config["Boundaries"]["Postprocessing"]["Dielectric"]],
                         [("SA", [28], None, True), ("MS", [114], None, False), ("MA", [114], None, False), ("MS", [126, 1002], [0.0, 0.0, -1.0], True), ("MA", [126, 1002], [0.0, 0.0, -1.0], True)])
        self.assertEqual(config["Solver"]["Electrostatic"]["ResponseCorrection"]["TargetInterfaces"], [1, 4, 5])
        # No gap sheet on any plane: no metal edge anywhere, refused.
        del manifest["Attributes"]["gap_L2"]
        with self.assertRaises(ValueError):
            thin_config(polygon_set, manifest, "m", "lib", "T", "p")

    def test_plane_without_sheet_refused(self):
        # A plane of the polygon set with no metal sheet in the mesh is inconsistent input.
        manifest = copy.deepcopy(THIN_MANIFEST)
        del manifest["Attributes"]["gap_L2"], manifest["Attributes"]["ground_L2"], manifest["Attributes"]["island_3_L2"]
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Planes"][1]["Polygons"] = [{"Conductor": "ground", "Outer": SQUARE}]
        polygon_set["Terminals"] = ["trace_7"]
        with self.assertRaises(ValueError):
            thin_config(polygon_set, manifest, "m", "lib", "T", "p")

    def test_terminal_without_sheet_refused(self):
        manifest = copy.deepcopy(THIN_MANIFEST)
        del manifest["Attributes"]["trace_7_L1"]
        with self.assertRaises(ValueError):
            thin_config(POLYGON_SET, manifest, "m", "lib", "T", "p")

    def test_terminal_bump_shells_in_its_terminal(self):
        # Decision 321: the thin mesher groups bump shells per conductor; `bump_<label>` joins the
        # label's Terminal (C4's bump 6 carries the terminal), `bump_surface` (ground) stays in Ground.
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Bumps"].append({"Conductor": "island_3", "Footprint": [[4.0, 4.0], [6.0, 4.0], [6.0, 6.0], [4.0, 6.0]]})
        manifest = copy.deepcopy(THIN_MANIFEST)
        manifest["Attributes"]["bump_island_3"] = 1003
        config = thin_config(polygon_set, manifest, "m", "lib", "T", "p")
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [114, 126, 153])
        self.assertEqual(config["Boundaries"]["Terminal"], [{"Index": 1, "Attributes": [1002, 1003]}, {"Index": 2, "Attributes": [1001]}])
        # The bump shells are no interface target: the dielectric sheets are unchanged.
        self.assertEqual([d["Attributes"] for d in config["Boundaries"]["Postprocessing"]["Dielectric"]], [[28, 145], [114, 1001], [114, 1001], [126, 1002], [126, 1002]])
        # Only terminal bumps: no bump_surface -> Ground is the ground sheets alone.
        del manifest["Attributes"]["bump_surface"]
        polygon_set["Bumps"] = polygon_set["Bumps"][1:]
        config = thin_config(polygon_set, manifest, "m", "lib", "T", "p")
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [114, 126])
        self.assertEqual(config["Boundaries"]["Terminal"][0]["Attributes"], [1002, 1003])

    def test_bump_group_without_terminal_refused(self):
        # Fail closed: a bump_<label> the Terminals do not name would otherwise be silently dropped.
        manifest = copy.deepcopy(THIN_MANIFEST)
        manifest["Attributes"]["bump_pad_9"] = 1013
        with self.assertRaises(ValueError):
            thin_config(POLYGON_SET, manifest, "m", "lib", "T", "p")
        manifest = copy.deepcopy(THIN_MANIFEST)
        manifest["Attributes"]["surface_L1"] = 1011
        polygon_set = copy.deepcopy(POLYGON_SET)
        polygon_set["Planes"][0]["Polygons"].append({"Conductor": "surface", "Outer": [[8.5, 8.5], [9.0, 8.5], [9.0, 9.0], [8.5, 9.0]]})
        polygon_set["Terminals"] = ["island_3", "trace_7", "surface"]
        with self.assertRaises(ValueError):
            thin_config(polygon_set, manifest, "m", "lib", "T", "p")


class MainTest(unittest.TestCase):
    def test_writes_both_kinds(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = {}
            for name, data in (("set.json", POLYGON_SET), ("fab.json", FABRICATED_MANIFEST), ("thin.json", THIN_MANIFEST), ("lib.json", {"MatchingRadius": 2.5})):
                paths[name] = os.path.join(directory, name)
                with open(paths[name], "w") as target:
                    json.dump(data, target)
            out = os.path.join(directory, "out", "ref.json")
            main(["fabricated", "--polygon-set", paths["set.json"], "--manifest", paths["fab.json"], "--mesh", "W.msh2", "--postpro", "pp", "--output", out, "--order", "5", "--estimator-cheap"])
            with open(out) as source:
                config = json.load(source)
            self.assertEqual(config["Solver"]["Order"], 5)
            self.assertEqual(config["Solver"]["Linear"]["EstimatorMaxIts"], 100)
            out = os.path.join(directory, "out", "thin.json")
            main(["thin", "--polygon-set", paths["set.json"], "--manifest", paths["thin.json"], "--mesh", "t.msh2", "--postpro", "pp", "--output", out, "--library", paths["lib.json"], "--path", "P4"])
            with open(out) as source:
                config = json.load(source)
            self.assertEqual(config["Solver"]["Order"], 4)
            # The library's MatchingRadius is authoritative for the edge patches.
            self.assertEqual(config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"], [2.5])
            self.assertEqual(config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"], paths["lib.json"])


if __name__ == "__main__":
    unittest.main()
