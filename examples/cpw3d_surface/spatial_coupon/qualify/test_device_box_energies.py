# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""device_box_energies: the device-side domain reading of the (F) twin-consistency check
(decision 299 (1a)) on a synthetic mesh - the placed-box membership, the inside / straddle
relabelling with overlapping boxes, the config extension and the energy record."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import device_box_energies as db  # noqa: E402


def cube_tetrahedra(nx, ny, nz, spacing=1.0):
    """A structured tetrahedral mesh of [0, nx] x [0, ny] x [0, nz] (5 tets per cell is not
    needed: 6 tets per cell through the Kuhn split); materials 1 below z = nz / 2, 3 above."""
    points = np.array([[i, j, k] for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)], dtype=float) * spacing

    def index(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)
    tets, materials = [], []
    kuhn = [(0, 1, 3, 7), (0, 1, 5, 7), (0, 2, 3, 7), (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 6, 7)]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corners = [index(i + a, j + b, k + c) for c in (0, 1) for b in (0, 1) for a in (0, 1)]
                for tet in kuhn:
                    tets.append([corners[v] for v in tet])
                    materials.append(1 if k < nz / 2 else 3)
    return points, np.array(tets), np.array(materials)


class RelabelTest(unittest.TestCase):
    def setUp(self):
        self.points, self.tets, self.materials = cube_tetrahedra(4, 4, 4)
        # A box of the model's canonical support [-1, 1] x [-1, 1] x [-1, 1] placed with its origin
        # at (1.5, 1.5, 2) and axes u = +x, v = -y, w = -z (the flip-chip frame of the S1p boxes).
        support = [[sx, sy, sz] for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)]
        self.box = {"Model": "m", "Origin": [1.5, 1.5, 2.0], "Axes": [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
                    "SupportPoints": support}
        self.other = {**self.box, "Model": "n", "Origin": [2.5, 1.5, 2.0]}  # overlaps the first box in x

    def test_membership_and_single_box(self):
        inside = db.box_membership(self.points, self.box)
        xyz = self.points
        expected = (np.abs(xyz[:, 0] - 1.5) <= 1 + 1e-9) & (np.abs(xyz[:, 1] - 1.5) <= 1 + 1e-9) & (np.abs(xyz[:, 2] - 2.0) <= 1 + 1e-9)
        np.testing.assert_array_equal(inside, expected)
        attributes, records, combinations = db.relabel_tetrahedra(self.points, self.tets, self.materials, [self.box])
        record = records[0]
        self.assertEqual(record["Index"], 1)
        # The box [0.5, 2.5] x [0.5, 2.5] x [1, 3] contains the unit cells [1, 2]^2 x [1, 3] fully: 2 cells x 6 tets.
        self.assertEqual(record["Tetrahedra"][db.STATUS_INSIDE], 12)
        self.assertAlmostEqual(record["Volume"][db.STATUS_INSIDE], 2.0)
        self.assertGreater(record["Tetrahedra"][db.STATUS_STRADDLE], 0)
        inside_attributes = {a for values in record["Attributes"][db.STATUS_INSIDE].values() for a in values}
        straddle_attributes = {a for values in record["Attributes"][db.STATUS_STRADDLE].values() for a in values}
        self.assertTrue(inside_attributes.isdisjoint(straddle_attributes))
        self.assertTrue(all(a > self.materials.max() for a in inside_attributes | straddle_attributes))
        untouched = ~np.isin(attributes, list(inside_attributes | straddle_attributes))
        np.testing.assert_array_equal(attributes[untouched], self.materials[untouched])
        # Every combination keeps its original material and the inside cells straddle z = 2 (both materials).
        self.assertEqual(sorted(record["Attributes"][db.STATUS_INSIDE]), ["1", "3"])
        self.assertEqual(sum(c["Tetrahedra"] for c in combinations), int((attributes > self.materials.max()).sum()))

    def test_overlapping_boxes_are_read_exactly(self):
        attributes, records, combinations = db.relabel_tetrahedra(self.points, self.tets, self.materials, [self.box, self.other])
        first, second = records
        shared = [c for c in combinations if all(code != 0 for code in c["BoxStatus"])]
        self.assertTrue(shared, "the overlapping boxes share tetrahedra")
        for record, column in ((first, 0), (second, 1)):
            for status, code in ((db.STATUS_INSIDE, 1), (db.STATUS_STRADDLE, 2)):
                listed = {a for values in record["Attributes"][status].values() for a in values}
                expected = {c["Attribute"] for c in combinations if c["BoxStatus"][column] == code}
                self.assertEqual(listed, expected)
                self.assertEqual(record["Tetrahedra"][status], sum(c["Tetrahedra"] for c in combinations
                                                                   if c["BoxStatus"][column] == code))
        self.assertEqual(first["Tetrahedra"][db.STATUS_INSIDE], second["Tetrahedra"][db.STATUS_INSIDE])
        # One attribute per (status combination, material): the attribute numbers continue the volume range.
        self.assertEqual([c["Attribute"] for c in combinations], list(range(4, 4 + len(combinations))))

    def test_config_and_record(self):
        _, records, combinations = db.relabel_tetrahedra(self.points, self.tets, self.materials, [self.box, self.other])
        boxes_record = {"Mesh": "boxes.msh", "Boxes": records, "Combinations": combinations}
        reference = {"Model": {"Mesh": "device.msh"}, "Problem": {"Output": "out"}, "Solver": {"Order": 3},
                     "Domains": {"Materials": [{"Attributes": [1, 2], "Permittivity": 11.45}, {"Attributes": [3], "Permittivity": 1.0}],
                                 "Postprocessing": {"Energy": [{"Index": 1, "Attributes": [1, 2]}, {"Index": 2, "Attributes": [3]}]}}}
        config = db.box_config(reference, boxes_record, "/cluster/boxes.msh", 5, "/cluster/p5/postpro")
        self.assertEqual((config["Model"]["Mesh"], config["Solver"]["Order"], config["Problem"]["Output"]),
                         ("/cluster/boxes.msh", 5, "/cluster/p5/postpro"))
        substrate = {c["Attribute"] for c in combinations if c["Material"] == 1}
        vacuum = {c["Attribute"] for c in combinations if c["Material"] == 3}
        self.assertEqual(set(config["Domains"]["Materials"][0]["Attributes"]), {1, 2} | substrate)
        self.assertEqual(set(config["Domains"]["Materials"][1]["Attributes"]), {3} | vacuum)
        energies = {e["Index"]: set(e["Attributes"]) for e in config["Domains"]["Postprocessing"]["Energy"]}
        self.assertEqual(energies[1], {1, 2} | substrate)
        self.assertEqual(energies[2], {3} | vacuum)
        self.assertEqual(sorted(energies), [1, 2, 1001, 1002, 2001, 2002])
        self.assertEqual(energies[1001], {a for v in records[0]["Attributes"][db.STATUS_INSIDE].values() for a in v})
        self.assertEqual(energies[2002], {a for v in records[1]["Attributes"][db.STATUS_STRADDLE].values() for a in v})
        self.assertEqual(reference["Domains"]["Materials"][0]["Attributes"], [1, 2])  # the reference is untouched
        with tempfile.TemporaryDirectory() as tmp:
            runs = {}
            for order, (e_in, e_straddle) in (("p4", (1.0, 0.5)), ("p5", (1.02, 0.5))):
                postpro = Path(tmp) / order
                postpro.mkdir()
                header = "i,E_elec (J),E_mag (J),E_elec[1] (J),E_elec[1001] (J),E_elec[1002] (J),E_elec[2001] (J),E_elec[2002] (J)"
                (postpro / "domain-E.csv").write_text(header + "\n" + ",".join(
                    f"{v:+.12e}" for v in (1, 10.0, 0.0, 6.0, e_in, e_straddle, 2 * e_in, e_straddle)) + "\n")
                runs[order] = str(postpro)
            twin = {"m": {"p4": 1.21, "p5": 1.215}}
            record = db.box_energy_record(boxes_record, runs, twin)
            self.assertEqual(record["Orders"], ["p4", "p5"])
            first = record["Boxes"][0]
            self.assertEqual(first["Device"]["p4"]["Bracket"], [1.0, 1.5])
            self.assertAlmostEqual(first["Device"]["p4"]["Central"], 1.25)
            self.assertAlmostEqual(first["DeviceStep"]["p4->p5"]["Central"], 0.02)
            self.assertAlmostEqual(first["DeviceStep"]["p4->p5"]["Inside"], 0.02)
            self.assertAlmostEqual(first["Twin"]["Step"]["p4->p5"], 0.005)
            self.assertAlmostEqual(first["Twin"]["TwinOverDevice"]["p4"]["Central"], 1.21 / 1.25)
            self.assertEqual(first["Twin"]["TwinOverDevice"]["p4"]["Bracket"], sorted([1.21 / 1.0, 1.21 / 1.5]))
            self.assertIsNone(record["Boxes"][1]["Twin"])  # no twin energies for model n
            self.assertAlmostEqual(record["Boxes"][1]["Device"]["p5"]["Inside"], 2.04)
            self.assertIn("INFORMATION ONLY", record["Rule"])
            with self.assertRaises(db.DeviceBoxError):
                db.domain_energies(runs["p4"], source=2)


if __name__ == "__main__":
    unittest.main()
