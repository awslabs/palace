# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The seed's edge layer construction (mesh_spatial_coupon.jl) follows the recorded
EdgeSize / GrowthRatio / EdgeLayerAspect rules the metric stage and the stage
contract bind."""
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


class EdgeLayerSeedJuliaTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        cls.project = REPO / "test" / "examples"

    def test_row_offsets_subdivisions_and_row_nodes_follow_the_bound_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            j = json.dumps
            code = f"""
include(joinpath({j(str(HERE))}, "mesh_spatial_coupon.jl"))
span = [0.0 0.05 0.10 0.15; 0.0 0.0 0.0 0.0; 0.1 0.1 0.1 0.1]
nodes16 = edge_layer_row_nodes(span, 0.001, [0.0, -1.0, 0.0], 16)
nodes1 = edge_layer_row_nodes(span, 0.031, [0.0, -1.0, 0.0], 1)
record = Dict{{String, Any}}(
    "offsets2" => edge_layer_row_offsets(0.001, 2.0, 0.025),
    "offsets15" => edge_layer_row_offsets(0.001, 1.5, 0.025),
    "offsetsNone" => edge_layer_row_offsets(0.025, 2.0, 0.025),
    "subdivisions" => [tangential_subdivision(0.05, 0.004), tangential_subdivision(0.05, 0.008),
                       tangential_subdivision(0.05, 0.016), tangential_subdivision(0.05, 0.064)],
    "nodes16" => [collect(nodes16[:, i]) for i in axes(nodes16, 2)],
    "nodes1" => [collect(nodes1[:, i]) for i in axes(nodes1, 2)],
    "zigzag" => EDGE_LAYER_ROW_ZIGZAG,
    "families" => collect(METAL_SURFACE_FAMILIES))
open({j(str(root / 'record.json'))}, "w") do io; write_json(io, record); end
"""
            result = subprocess.run([self.julia, "--startup-file=no", f"--project={self.project}",
                                     "-e", code], cwd=REPO, capture_output=True, text=True,
                                    check=False, timeout=600)
            self.assertEqual(result.returncode, 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
            record = json.loads((root / "record.json").read_text())
        # Geometric layers strictly below NormalSize: 1, 2, 4, 8, 16 nm (ratio 2) at the
        # cumulative distances 1, 3, 7, 15, 31 nm; eight layers for ratio 1.5; none when
        # EdgeSize equals NormalSize.  The contract recomputes the same sequence.
        np.testing.assert_allclose(record["offsets2"], [.001, .003, .007, .015, .031], rtol=1e-12)
        self.assertEqual(len(record["offsets15"]), 8)
        np.testing.assert_allclose(record["offsets15"][:3], [.001, .0025, .00475], rtol=1e-12)
        self.assertEqual(record["offsetsNone"], [])
        # Nested powers of two: 3.125 / 6.25 / 12.5 / 50 nm spacing for targets 4 / 8 / 16 / 64 nm.
        self.assertEqual(record["subdivisions"], [16, 8, 4, 1])
        # Row nodes: the span grid subdivided (endpoints excluded), offset by the row
        # distance with every second node RowZigzag farther from the ridge.
        nodes16 = np.asarray(record["nodes16"])
        self.assertEqual(len(nodes16), 3 * 16 - 1)
        np.testing.assert_allclose(nodes16[:, 0], np.arange(1, 48) * .05 / 16, rtol=1e-12)
        np.testing.assert_allclose(nodes16[::2, 1], -.001 * (1. + record["zigzag"]), rtol=1e-12)
        np.testing.assert_allclose(nodes16[1::2, 1], -.001, rtol=1e-12)
        np.testing.assert_allclose(nodes16[:, 2], .1)
        nodes1 = np.asarray(record["nodes1"])
        np.testing.assert_allclose(nodes1[:, 0], [.05, .10], rtol=1e-12)
        self.assertEqual(record["families"], [4, 5, 6])


if __name__ == "__main__":
    unittest.main()
