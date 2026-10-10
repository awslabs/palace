# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The BOUND RING in mesh_stage_contract (round 4 F-3 redesign, decisions 629 / 630 (4) (i)): a tube row
carrying Tubes[].BoundRing reads its node envelope from BoundRing.EnvelopeRadius (the farthest bound node,
never below Radius + PyramidHeight) in the arc face-end crossing slope / node-circle margin / end-block
layer count and in the arc curvature term; the ring's radii are re-derived from |plane| / |sin theta|
within the reach; a row without the key keeps the plain envelope bitwise."""
import copy
import math
import unittest

from mesh_stage_contract import (arc_crossing_slope, tube_row_envelope, validate_arc_tubes,
                                 validate_tube_face_ends)
from test_mesh_stage_contract_arcs import arc_census, tagged_loop


def production_bound_ring(envelope=0.05 / math.sin(math.radians(60.0))):
    r60 = 0.05 / math.sin(math.radians(60.0))
    return {"Reach": 0.03175 + 2.0 * 0.016, "Rays": [240.0, 270.0, 300.0, 420.0, 450.0],
            "Radii": [r60, 0.05, r60, r60, 0.05], "Planes": [-0.05, -0.05, -0.05, 0.05, 0.05],
            "CoveredSectors": [[240.0, 270.0], [270.0, 300.0], [420.0, 450.0]],
            "ShoulderSectors": [[210.0, 240.0], [300.0, 330.0], [390.0, 420.0]],
            "Pivots": [240.0, 240.0, 300.0, 300.0, 420.0, 420.0], "UnpairedSectors": [],
            "MaximumRadius": r60, "EnvelopeRadius": envelope}


class TubeRowEnvelopeTest(unittest.TestCase):
    def test_plain_row_keeps_the_plain_envelope(self):
        self.assertEqual(tube_row_envelope({"Length": 1.0}, 0.03975), 0.03975)

    def test_bound_ring_row_reads_its_envelope(self):
        row = {"BoundRing": production_bound_ring()}
        self.assertAlmostEqual(tube_row_envelope(row, 0.03975), 0.05 / math.sin(math.radians(60.0)), places=15)

    def test_bound_ring_fails_closed(self):
        def rejected(mutate, message):
            ring = production_bound_ring()
            mutate(ring)
            with self.assertRaisesRegex(ValueError, message):
                tube_row_envelope({"BoundRing": ring}, 0.03975)

        rejected(lambda r: r.__setitem__("EnvelopeRadius", 0.039), "below the section's Radius")
        rejected(lambda r: r.pop("EnvelopeRadius"), "EnvelopeRadius")
        rejected(lambda r: r.pop("CoveredSectors"), "lacks its CoveredSectors")
        rejected(lambda r: r["Radii"].__setitem__(1, 0.051), "does not follow")
        rejected(lambda r: r["Radii"].pop(), "differ in length")
        rejected(lambda r: r.__setitem__("Reach", 0.05), "does not follow")
        rejected(lambda r: r["Planes"].__setitem__(0, 0.05), "does not follow")
        with self.assertRaisesRegex(ValueError, "not a record"):
            tube_row_envelope({"BoundRing": 3}, 0.03975)


class BoundRingFaceEndTest(unittest.TestCase):
    """An arc face end on a bound-ring tube: every formula reads the bound envelope (the plain-envelope
    record of the same geometry is REJECTED by name, and vice versa)."""

    def face_end_row(self, bound, theta=45.0):
        radius, h_pyr, lc = 0.03175, 0.0079375, 0.05
        envelope = 0.05 / math.sin(math.radians(60.0)) if bound else radius + h_pyr
        rho, centre = 1.5 / math.sin(math.radians(theta)), [0.0, 0.0]
        box = ([-2.0, -1.0, -1.0], [1.5, 3.0, 1.0])
        slope = arc_crossing_slope(rho, envelope, 1.5)
        lc_end = max(lc, 4.0 * h_pyr * slope)
        m = max(1, math.ceil(2.0 * envelope * slope / lc_end * (1.0 - 1e-9)))
        shear = envelope * slope
        margin = (rho - 1.5) / envelope
        record = {"Face": "x1", "End": "end", "ThetaDegrees": theta, "Layers": m, "EndSpacing": lc_end,
                  "EnvelopeShear": shear, "LayerThicknessRange": [lc_end - shear / m, lc_end + shear / m],
                  "OverLength": shear + lc, "Kappa": [0.0, 0.0], "CrossingSlope": slope, "NodeCircleMargin": margin}
        row = {"Arc": {"ArcId": 1, "Centre": centre, "Radius": rho, "Sign": 1, "Part": 1, "Parts": 1},
               "FaceEnds": [record]}
        if bound:
            row["BoundRing"] = production_bound_ring()
            row["StationLaw"] = "SideShared(bottom axis)"
        return row, lc, {"Radius": radius, "PyramidHeight": h_pyr}, box, lc_end, margin

    def test_bound_row_validates_with_the_bound_envelope(self):
        row, lc, section, box, lc_end, margin = self.face_end_row(True)
        self.assertEqual(validate_tube_face_ends(row, lc, section, box=box, fabricated=True), lc_end)
        # 15.6 plain envelopes -> 10.76 bound envelopes at 45 degrees: still dominated (fe75p5r13p3 10.66).
        self.assertAlmostEqual(margin, 10.76, places=1)
        # The layer count of the bound end block exceeds the plain row's (a larger envelope shear).
        plain, *_ = self.face_end_row(False)
        self.assertGreater(row["FaceEnds"][0]["Layers"], plain["FaceEnds"][0]["Layers"] - 1)
        self.assertGreater(row["FaceEnds"][0]["EnvelopeShear"], plain["FaceEnds"][0]["EnvelopeShear"])

    def test_mixed_records_are_rejected(self):
        row, lc, section, box, _, _ = self.face_end_row(True)
        plain, *_ = self.face_end_row(False)
        # the plain geometry's record on a bound-ring row, and the bound record on a plain row
        mixed = copy.deepcopy(row)
        mixed["FaceEnds"] = plain["FaceEnds"]
        with self.assertRaisesRegex(ValueError, "NodeCircleMargin does not follow"):
            validate_tube_face_ends(mixed, lc, section, box=box)
        mixed = copy.deepcopy(plain)
        mixed["FaceEnds"] = row["FaceEnds"]
        with self.assertRaisesRegex(ValueError, "NodeCircleMargin does not follow"):
            validate_tube_face_ends(mixed, lc, section, box=box)

    def test_arc_curvature_reads_the_bound_envelope(self):
        rows = tagged_loop(90.0, 27)
        tubes, tube_rows = arc_census(rows, 1, 90.0)
        tubes["Tubes"][0]["BoundRing"] = production_bound_ring()
        self.assertEqual(validate_arc_tubes(tubes, tubes["Tubes"], rows), 1)
        # an envelope above 0.25 x the arc radius (rho 1.0 here) fails closed by the same message
        tubes["Tubes"][0]["BoundRing"]["EnvelopeRadius"] = 0.26
        with self.assertRaisesRegex(ValueError, "does not follow its tagged arc"):
            validate_arc_tubes(tubes, tubes["Tubes"], rows)


if __name__ == "__main__":
    unittest.main()
