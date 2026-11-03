#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""attach_magnetic_leaves.py (mx-sc-magnetic DESIGN-v2 section 3; decision 1229): the vacuum twin's leaves are attached
to a COPY of a Version-3 library with the MagneticSource provenance and the MagneticRules; a twin on another basis is
refused (r1); the r3 same-mesh statistic and the two-level r4 PSD statistic are recorded with the verdict the loader
will reach; the input library is not rewritten."""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

CPW2D = Path(__file__).resolve().parent
sys.path.insert(0, str(CPW2D))
import attach_magnetic_leaves as AML  # noqa: E402


def write_matrix(path, matrix):
    with open(path, "w") as handle:
        handle.write("  basis_i,  basis_j,                   Q_ij (J)\n")
        n = matrix.shape[0]
        for i in range(n):
            for j in range(i, n):
                handle.write(" %.2e, %.2e, %+.12e\n" % (i + 1, j + 1, matrix[i, j]))


def write_points(path, count, shift=0.0):
    with open(path, "w") as handle:
        handle.write("x,y,z\n")
        for k in range(count):
            handle.write("%.16e,%.16e,0.0\n" % (0.1 * k + shift, -0.2))


class AttachMagneticLeavesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.root = root
        rng = np.random.default_rng(3)
        n = 5  # 4 hats + 1 conductor state
        a = rng.normal(size=(n, n))
        self.thin = a @ a.T * 1e-13
        self.fab = self.thin * 0.9
        write_points(root / "basis-points.csv", 4)
        write_matrix(root / "thin.csv", self.thin)
        write_matrix(root / "fab.csv", self.fab)
        self.library = {
            "Version": 3, "TraceLiftVersion": 3, "Name": "unit", "MatchingRadius": 1.9,
            "Fabrication": {"LengthUnit": "um", "MetalThickness": 0.1, "OveretchDepth": 0.05,
                            "SidewallAngleDegrees": 90.0, "TopRoundingRadius": 0.0, "BottomRoundingRadius": 0.0,
                            "SubstratePermittivity": 11.45,
                            "InterfaceLayers": {"SA": {"Thickness": 0.002, "Permittivity": 4.0}}},
            "Models": [{"Name": "pair", "Topology": "DifferentConductorGap", "Separation": 2.0,
                        "FabricatedMatrix": "fab.csv", "ThinMatrix": "thin.csv", "BasisPoints": "basis-points.csv",
                        "ConductorReferences": [[-0.95, 0.0, 0.0], [2.95, 0.0, 0.0]],
                        "Interfaces": [{"Type": "SA", "Coupon": 1}]}]}
        with open(root / "library.json", "w") as handle:
            json.dump(self.library, handle)
        self.library_sha = AML.sha256_of(root / "library.json")
        # The vacuum twin: S_fab - S_thin PSD (the duality's sign); the regenerated electric pair = the record's.
        b = rng.normal(size=(n, n))
        self.s_thin = b @ b.T * 1e-14
        c = rng.normal(size=(n, n))
        self.s_fab = self.s_thin + c @ c.T * 1e-15
        write_matrix(root / "vac-thin.csv", self.s_thin)
        write_matrix(root / "vac-fab.csv", self.s_fab)
        write_points(root / "vac-basis-points.csv", 4)
        write_points(root / "other-basis-points.csv", 4, shift=1.0e-3)

    def tearDown(self):
        self.tmp.cleanup()

    def manifest(self, **overrides):
        entry = {"Model": "pair", "VacuumThin": str(self.root / "vac-thin.csv"),
                 "VacuumFabricated": str(self.root / "vac-fab.csv"),
                 "VacuumBasisPoints": str(self.root / "vac-basis-points.csv"),
                 "RegeneratedElectricThin": str(self.root / "thin.csv"),
                 "RegeneratedElectricFabricated": str(self.root / "fab.csv"),
                 "MeshSha256": {"Thin": "aa", "Fabricated": "bb"}, "Toolchain": "julia 1.10.4 / gmsh 4.13.1"}
        entry.update(overrides)
        path = self.root / "attach.json"
        with open(path, "w") as handle:
            json.dump([entry], handle)
        return str(path)

    def test_reader_symmetrises_the_upper_triangle(self):
        matrix = AML.read_domain_response_matrix(self.root / "thin.csv")
        self.assertEqual(matrix.shape, (5, 5))
        self.assertTrue(np.allclose(matrix, matrix.T))
        self.assertTrue(np.allclose(matrix, self.thin, rtol=1e-12, atol=0.0))

    def test_attach_writes_leaves_rules_and_provenance(self):
        out_dir = self.root / "out"
        rc = AML.main(["--library", str(self.root / "library.json"), "--attach", self.manifest(),
                       "--output", str(out_dir / "library-mag.json")])
        self.assertEqual(rc, 0)
        with open(out_dir / "library-mag.json") as handle:
            out = json.load(handle)
        model = out["Models"][0]
        self.assertEqual(out["MagneticRules"], AML.DEFAULT_RULES)
        self.assertEqual(out["Name"], "unit-mag")
        # The leaf paths are RELATIVE to the output library's directory (relocatable).
        self.assertFalse(os.path.isabs(model["MagneticFabricatedMatrix"]))
        self.assertEqual(model["MagneticFabricatedMatrix"],
                         os.path.join("models", "pair", "magnetic-fabricated-domain-response-matrix.csv"))
        self.assertTrue(os.path.exists(out_dir / model["MagneticThinMatrix"]))
        source = model["MagneticSource"]
        self.assertEqual(source["Version"], 1)
        self.assertEqual(source["VacuumPermittivity"], 1.0)
        self.assertEqual(source["Fabrication"],
                         {"MetalThickness": 0.1, "OveretchDepth": 0.05, "SidewallAngle": 90.0,
                          "TopRoundingRadius": 0.0, "BottomRoundingRadius": 0.0, "SubstratePermittivity": 11.45})
        self.assertEqual(source["ElectricRegenerationMaxRelativeDeviation"], 0.0)
        self.assertEqual((source["BasisPointCount"], source["ConductorStateCount"]), (4, 1))
        self.assertTrue(source["BasisPointsIdentical"])
        self.assertEqual(source["NegativeEigenvalues"], 0)
        self.assertLess(source["PsdRatio"], 1e-12)
        self.assertEqual(source["MeshSha256"], {"Thin": "aa", "Fabricated": "bb"})
        # The electric entries and the input file are untouched.
        self.assertEqual(model["FabricatedMatrix"], "fab.csv")
        self.assertEqual(AML.sha256_of(self.root / "library.json"), self.library_sha)
        with open(out_dir / "library-mag-magnetic-records.json") as handle:
            records = json.load(handle)
        self.assertEqual(records["Models"][0]["Verdict"], "Applied")
        # The copied leaf is the twin file byte for byte.
        self.assertEqual(AML.sha256_of(out_dir / model["MagneticThinMatrix"]),
                         AML.sha256_of(self.root / "vac-thin.csv"))

    def test_r1_refuses_a_twin_on_another_basis(self):
        with self.assertRaises(ValueError) as context:
            AML.main(["--library", str(self.root / "library.json"),
                      "--attach", self.manifest(VacuumBasisPoints=str(self.root / "other-basis-points.csv")),
                      "--output", str(self.root / "out-r1" / "library-mag.json")])
        self.assertIn("r1", str(context.exception))

    def test_r3_and_r4_statistics_and_verdicts(self):
        # A regenerated pair that differs: the deviation is recorded and the verdict reads Refused r3.
        write_matrix(self.root / "regen-thin.csv", self.thin * (1.0 + 8.1e-5))
        out_dir = self.root / "out-r3"
        AML.main(["--library", str(self.root / "library.json"),
                  "--attach", self.manifest(RegeneratedElectricThin=str(self.root / "regen-thin.csv")),
                  "--output", str(out_dir / "library-mag.json")])
        with open(out_dir / "library-mag-magnetic-records.json") as handle:
            record = json.load(handle)["Models"][0]
        self.assertAlmostEqual(record["ElectricRegenerationMaxRelativeDeviation"], 8.1e-5, places=9)
        self.assertEqual(record["Verdict"], "Refused r3")
        # Without a regenerated pair the deviation is -1 (not recorded): Refused r3 by the default tolerance.
        out_dir = self.root / "out-norec"
        AML.main(["--library", str(self.root / "library.json"),
                  "--attach", self.manifest(RegeneratedElectricThin=None, RegeneratedElectricFabricated=None),
                  "--output", str(out_dir / "library-mag.json")])
        with open(out_dir / "library-mag-magnetic-records.json") as handle:
            self.assertEqual(json.load(handle)["Models"][0]["ElectricRegenerationMaxRelativeDeviation"], -1.0)
        # r4 two-level: a defect with a negative eigenvalue at 5 % (Finding) and at 50 % (Refused r4), and the
        # library-level ratio moving the line.
        for ratio, verdict, extra in ((0.05, "Finding", []), (0.5, "Refused r4", []),
                                      (0.5, "Finding", ["--psd-refusal-ratio", "0.9"])):
            s_fab = self.s_thin + np.diag([1.0, 1.0, 1.0, 1.0, -ratio]) * 1e-15
            write_matrix(self.root / "vac-fab-neg.csv", s_fab)
            out_dir = self.root / ("out-r4-%s-%d" % (verdict.replace(" ", ""), len(extra)))
            AML.main(["--library", str(self.root / "library.json"),
                      "--attach", self.manifest(VacuumFabricated=str(self.root / "vac-fab-neg.csv")),
                      "--output", str(out_dir / "library-mag.json")] + extra)
            with open(out_dir / "library-mag-magnetic-records.json") as handle:
                record = json.load(handle)["Models"][0]
            self.assertEqual(record["Verdict"], verdict, record)
            self.assertEqual(record["NegativeEigenvalues"], 1)
            self.assertAlmostEqual(record["PsdRatio"], ratio, places=12)

    def test_refuses_an_existing_output_and_a_blend_entry(self):
        with self.assertRaises(SystemExit):
            AML.main(["--library", str(self.root / "library.json"), "--attach", self.manifest(),
                      "--output", str(self.root / "library.json")])
        library = dict(self.library)
        library["Models"] = [dict(self.library["Models"][0], Blend=[{"Weight": 1.0, "FabricatedMatrix": "fab.csv",
                                                                       "ThinMatrix": "thin.csv"}])]
        with open(self.root / "blend.json", "w") as handle:
            json.dump(library, handle)
        with self.assertRaises(ValueError):
            AML.main(["--library", str(self.root / "blend.json"), "--attach", self.manifest(),
                      "--output", str(self.root / "out-blend" / "library-mag.json")])


if __name__ == "__main__":
    unittest.main()
