#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the generic geometry-only preflight configuration writer."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import perimeter as P  # noqa: E402
from surface_response_identification.preflight_config import main, preflight_config  # noqa: E402

LIBRARY = {
    "Version": 3,
    "MatchingRadius": 2.5,
    "LengthUnit": "um",
    "Models": [],
    "Fabrication": {
        "SubstratePermittivity": 11.7,
        "InterfaceLayers": {
            "SA": {"Thickness": 0.003, "Permittivity": 4.5},
            "MS": {"Thickness": 0.002, "Permittivity": 11.47},
            "MA": {"Thickness": 0.002, "Permittivity": 10.0},
        },
    },
}


class PreflightConfigTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.library = os.path.join(self.directory.name, "library.json")
        with open(self.library, "w") as target:
            json.dump(LIBRARY, target)

    def tearDown(self):
        self.directory.cleanup()

    def test_metal_interfaces_and_targets_follow_the_arguments(self):
        config = preflight_config("m.msh2", ground=[5, 6], terminals=[[9], [10, 11]], sa=[8], library=self.library, output_directory="out")
        self.assertEqual(config["Boundaries"]["Ground"]["Attributes"], [5, 6])
        self.assertEqual(config["Boundaries"]["Terminal"], [{"Index": 1, "Attributes": [9]}, {"Index": 2, "Attributes": [10, 11]}])
        # The audit's metal set (metaledge.cpp union) sees exactly ground + terminals.
        self.assertEqual(sorted(P.metal_attributes(config)), [5, 6, 9, 10, 11])
        self.assertEqual(P.conductor_of_attribute(config, 11), 2)
        interfaces = P.interface_attributes(config)
        self.assertEqual(interfaces[(1, "SA")], [8])
        self.assertEqual(interfaces[(2, "MS")], [5, 6, 9, 10, 11])
        self.assertEqual(interfaces[(3, "MA")], [5, 6, 9, 10, 11])
        self.assertEqual(P.target_interfaces(config), [1, 2, 3])
        # Layer data and the radius come from the library, not from the defaults.
        for entry in config["Boundaries"]["Postprocessing"]["Dielectric"]:
            self.assertEqual(entry["EdgeDistances"], [2.5])
            self.assertTrue(entry["AutomaticEdges"])
        sa = config["Boundaries"]["Postprocessing"]["Dielectric"][0]
        self.assertEqual((sa["Thickness"], sa["Permittivity"]), (0.003, 4.5))
        self.assertEqual(config["Domains"]["Materials"][0]["Permittivity"], 11.7)
        self.assertEqual(config["Problem"]["Type"], "Electrostatic")
        self.assertEqual(config["Model"]["Refinement"]["UniformLevels"], 0)
        self.assertNotIn("CrackInternalBoundaryElements", config["Model"])

    def test_without_sa_only_ms_and_ma_are_targeted(self):
        config = preflight_config("m.msh2", ground=[4], terminals=[], sa=[], library=self.library, output_directory="out", ms=[4], ma=[4], uniform_levels=1, crack=False)
        self.assertNotIn("Terminal", config["Boundaries"])
        self.assertEqual([e["Type"] for e in config["Boundaries"]["Postprocessing"]["Dielectric"]], ["MS", "MA"])
        self.assertEqual(P.target_interfaces(config), [2, 3])
        self.assertEqual(config["Model"]["Refinement"]["UniformLevels"], 1)
        self.assertIs(config["Model"]["CrackInternalBoundaryElements"], False)

    def test_command_line_writes_the_file(self):
        output = os.path.join(self.directory.name, "cfg", "config.json")
        main(["--mesh", "m.msh2", "--ground", "5", "--terminal", "9", "--sa", "8", "--library", self.library, "--output", output, "--no-crack", "--l0", "1e-3"])
        with open(output) as source:
            config = json.load(source)
        self.assertEqual(config["Model"]["L0"], 1.0e-3)
        self.assertEqual(config["Problem"]["Output"], os.path.join(self.directory.name, "cfg", "postpro"))
        self.assertIs(config["Model"]["CrackInternalBoundaryElements"], False)

    def test_planes_give_each_metal_plane_its_own_targets_and_frame_normal(self):
        config = preflight_config("m.msh2", ground=[114, 126, 10], terminals=[], sa=[28, 145], library=self.library, output_directory="out",
                                  substrate_attributes=[1, 2], vacuum_attributes=[3, 4],
                                  planes=[([114], [0, 0, 1]), ([126, 10], [0, 0, -1])])
        dielectrics = config["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertEqual([(d["Index"], d["Type"], d["Attributes"], d.get("EdgeFrameNormal")) for d in dielectrics], [
            (1, "SA", [28, 145], None),
            (2, "MS", [114], [0.0, 0.0, 1.0]), (3, "MA", [114], [0.0, 0.0, 1.0]),
            (4, "MS", [10, 126], [0.0, 0.0, -1.0]), (5, "MA", [10, 126], [0.0, 0.0, -1.0]),
        ])
        self.assertEqual(P.target_interfaces(config), [1, 2, 3, 4, 5])
        self.assertEqual(config["Domains"]["Materials"][0]["Attributes"], [1, 2])
        self.assertEqual(config["Domains"]["Materials"][1]["Attributes"], [3, 4])

    def test_plane_option_takes_disjoint_metal_subsets(self):
        output = os.path.join(self.directory.name, "cfg.json")
        with self.assertRaises(SystemExit):  # 153 is not metal
            main(["--mesh", "m.msh2", "--ground", "114", "126", "--library", self.library, "--output", output,
                  "--plane", "114,153@0,0,1"])
        with self.assertRaises(SystemExit):  # 114 twice
            main(["--mesh", "m.msh2", "--ground", "114", "126", "--library", self.library, "--output", output,
                  "--plane", "114@0,0,1", "--plane", "114,126@0,0,-1"])
        # Metal outside every plane (153: bumps) stays untargeted.
        main(["--mesh", "m.msh2", "--ground", "114", "126", "153", "--library", self.library, "--output", output,
              "--plane", "114@0,0,1", "--plane", "126@0,0,-1"])
        with open(output) as source:
            config = json.load(source)
        normals = [d.get("EdgeFrameNormal") for d in config["Boundaries"]["Postprocessing"]["Dielectric"]]
        self.assertEqual(normals, [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])


if __name__ == "__main__":
    unittest.main()
