# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""compare_case_records: two builds of one case are identical when only timings, resources,
tool digests, paths and commands differ; any other field or mesh byte is reported."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_case_records  # noqa: E402


def write_case(root, *, seconds, value, mesh, producer):
    root.mkdir(parents=True)
    (root / "identity-mesh-topology-quality.json").write_text(json.dumps({
        "Kind": "mesh-topology-quality", "Producer": {"Name": "p.py", "SHA256": producer},
        "Command": ["python3", str(root / "x")], "Measurements": {"MeshQuality": {"MinimumScaledJacobian": value}},
        "Nested": [{"Seconds": seconds, "Value": 1.5}]}))
    (root / "gmsh-build.log.json").write_text(json.dumps({"Seconds": seconds, "PeakProcessTreeRSSBytes": seconds * 10,
                                                          "Stage": "gmsh-build"}))
    (root / "identity.msh").write_bytes(mesh)


class CompareCaseRecordsTest(unittest.TestCase):
    def test_timings_tools_and_paths_are_ignored_values_and_bytes_are_not(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_case(root / "before", seconds=10.0, value=0.25, mesh=b"mesh-a", producer="aa")
            write_case(root / "same", seconds=12.5, value=0.25, mesh=b"mesh-a", producer="bb")
            write_case(root / "value", seconds=10.0, value=0.2500000001, mesh=b"mesh-a", producer="aa")
            write_case(root / "mesh", seconds=10.0, value=0.25, mesh=b"mesh-b", producer="aa")
            report = compare_case_records.compare(root / "before", root / "same")
            self.assertTrue(report["Identical"], report)
            self.assertEqual(report["ByteIdentical"], ["identity.msh"])
            self.assertEqual(sorted(report["Compared"]), ["gmsh-build.log.json", "identity-mesh-topology-quality.json"])
            report = compare_case_records.compare(root / "before", root / "value")
            self.assertFalse(report["Identical"])
            self.assertEqual([item["Path"] for item in report["Differences"]],
                             ["identity-mesh-topology-quality.json/Measurements/MeshQuality/MinimumScaledJacobian"])
            report = compare_case_records.compare(root / "before", root / "mesh")
            self.assertEqual((report["Identical"], report["ByteDifferent"]), (False, ["identity.msh"]))
            self.assertEqual(compare_case_records.main([str(root / "before"), str(root / "same")]), 0)
            self.assertEqual(compare_case_records.main([str(root / "before"), str(root / "value")]), 1)


if __name__ == "__main__":
    unittest.main()
