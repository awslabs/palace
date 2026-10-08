#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Block (b) step 4.3 (T2 arcs): the manifest contract's arc bindings - boundary_arc_runs /
arc_part_count / metal_loop_side_points / metal_loop_arc_parts on tagged plan-view boundary rows,
scope_classes ArcSides, validate_arc_tubes positives and negatives."""
import copy
import math
from pathlib import Path
import re
import unittest

from mesh_stage_contract import (ARC_CORNER_JOINT_TURN_RANGE_RADIANS, ARC_FACE_END_TILT_RANGE_DEGREES,
                                 ARC_JOINT_TURN_BOUND_RADIANS, ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS,
                                 FABRICATED_FACE_END_TILT_BOUND_DEGREES, FACE_END_DERIVED_APEX_ABOVE_DEGREES,
                                 RECIPE_SCOPE_GUARDS, THIN_FACE_END_TILT_BOUND_DEGREES, arc_crossing_slope,
                                 arc_part_count, boundary_arc_runs, metal_loop_arc_parts, metal_loop_side_points,
                                 scope_classes, scope_guard_in_text, validate_arc_tubes, validate_tube_face_ends)

HERE = Path(__file__).resolve().parent

HEADER = ("Loop", "Vertex", "Conductor", "Plane", "Hole", "Class", "X", "Y", "ArcId", "ArcCx", "ArcCy", "ArcR", "ArcSign",
          "JointTurn", "JointSmooth")


def tagged_loop(sweep_degrees, chords, radius=1.0):
    """A loop: a straight side from the x0 face (-3, 0) to (0, 0), an arc of `sweep_degrees` about
    (0, radius) from angle -90 upwards (convex, sign +1) chorded `chords` times, a straight side
    back to the box corner region and the box side closing it."""
    centre = (0.0, radius)
    points = [(-3.0, 0.0)]
    for k in range(chords + 1):
        theta = math.radians(-90.0 + sweep_degrees * k / chords)
        points.append((centre[0] + radius * math.cos(theta), centre[1] + radius * math.sin(theta)))
    points.append((3.0, 4.0))
    points.append((-3.0, 4.0))
    rows = []
    n = len(points)
    for i, (x, y) in enumerate(points):
        arc = 1 <= i <= chords
        row = {"Loop": "1", "Vertex": str(i + 1), "Conductor": "1", "Plane": "0.0", "Hole": "0",
               "Class": "Continuation" if i == n - 1 else "Physical", "X": repr(x), "Y": repr(y),
               "ArcId": "1" if arc else "", "ArcCx": repr(centre[0]) if arc else "", "ArcCy": repr(centre[1]) if arc else "",
               "ArcR": repr(radius) if arc else "", "ArcSign": "1" if arc else "",
               "JointTurn": "0.0" if i in (1, chords + 1) else "", "JointSmooth": "1" if i in (1, chords + 1) else ""}
        rows.append(row)
    return rows


class ArcRunsTest(unittest.TestCase):
    def test_part_count_rule(self):
        self.assertEqual(arc_part_count(math.radians(90.0)), 1)
        self.assertEqual(arc_part_count(math.radians(90.000001)), 2)
        self.assertEqual(arc_part_count(math.radians(-180.0)), 2)
        self.assertEqual(arc_part_count(math.radians(180.000001)), 3)
        self.assertEqual(arc_part_count(math.radians(10.0)), 1)

    def test_tagged_runs_sweep_and_parts_do_not_depend_on_the_chord_count(self):
        for chords in (18, 36):
            runs = boundary_arc_runs(tagged_loop(90.0, chords))
            self.assertEqual(len(runs), 1)
            self.assertEqual(len(runs[0]), 1)
            run = runs[0][0]
            self.assertEqual(run["ArcId"], 1)
            self.assertEqual(run["Chords"], chords)
            self.assertAlmostEqual(run["Sweep"], math.pi / 2, places=12)
            self.assertEqual(run["Parts"], 1)
            self.assertEqual(run["Sign"], 1)
        run = boundary_arc_runs(tagged_loop(135.0, 27))[0][0]
        self.assertAlmostEqual(math.degrees(run["Sweep"]), 135.0, places=9)
        self.assertEqual(run["Parts"], 2)

    def test_sides_and_parts(self):
        rows = tagged_loop(90.0, 18)
        lower, upper = (-3.0, -1.0), (3.0, 4.0)
        sides = metal_loop_side_points(rows, lower, upper, 1e-9)
        self.assertEqual(len(sides), 1)
        # The two straight sides (the x0 face side is a box side; the chords are not sides).
        self.assertEqual(len(sides[0]), 2)
        self.assertEqual(metal_loop_arc_parts(rows), [1])
        self.assertEqual(metal_loop_arc_parts(tagged_loop(135.0, 27)), [2])
        legacy = [{k: v for k, v in row.items() if k not in ("ArcId", "ArcCx", "ArcCy", "ArcR", "ArcSign", "JointTurn", "JointSmooth")}
                  for row in rows]
        self.assertEqual(boundary_arc_runs(legacy), [[]])
        self.assertEqual(len(metal_loop_side_points(legacy, lower, upper, 1e-9)[0]), 2 + 18)

    def test_scope_class(self):
        signature = [{"Slot": "0", "Conductor": "1", "Pz": "0.0", "Nz": "1"}]
        rows = tagged_loop(90.0, 18)
        self.assertIn("ArcSides", scope_classes(signature, rows, overetch=0.05))
        legacy = [{k: v for k, v in row.items() if k not in ("ArcId", "ArcCx", "ArcCy", "ArcR", "ArcSign", "JointTurn", "JointSmooth")}
                  for row in rows]
        self.assertNotIn("ArcSides", scope_classes(signature, legacy, overetch=0.05))

    def test_non_consecutive_run_fails_closed(self):
        rows = tagged_loop(90.0, 18)
        rows[5]["ArcId"] = ""
        with self.assertRaisesRegex(ValueError, "not one consecutive run"):
            boundary_arc_runs(rows)


def arc_census(rows, parts, sweep_degrees):
    section = {"Radius": 0.03175, "PyramidHeight": 0.0079375}
    tube_rows = []
    for part in range(1, parts + 1):
        tube_rows.append({"Arc": {"ArcId": 1, "Centre": [0.0, 1.0], "Radius": 1.0, "Sign": 1, "Part": part, "Parts": parts,
                                  "SweepDegrees": sweep_degrees / parts}, "Length": math.radians(sweep_degrees) / parts})
    straight = [{"Length": 3.0, "Joints": [{"End": "end", "TiltRadians": 0.0, "PlaneCut": False}]},
                {"Length": 3.0, "Joints": [{"End": "start", "TiltRadians": 1.0e-6, "PlaneCut": True}]}]
    tube_rows += straight
    tubes = {"Section": section, "Tubes": tube_rows,
             "ArcTubes": {"Count": parts, "JointEnds": 2, "PartSplits": parts - 1, "SharedSections": 2 + parts - 1,
                          "TotalArcLength": math.radians(sweep_degrees), "Rule": "arc", "SmoothJointRule": "joint"}}
    return tubes, tube_rows


class ValidateArcTubesTest(unittest.TestCase):
    def test_positive(self):
        for sweep, parts in ((90.0, 1), (135.0, 2)):
            rows = tagged_loop(sweep, 27)
            tubes, tube_rows = arc_census(rows, parts, sweep)
            self.assertEqual(validate_arc_tubes(tubes, tube_rows, rows), parts)

    def test_no_arcs(self):
        rows = tagged_loop(90.0, 18)
        tubes = {"Section": {"Radius": 0.03, "PyramidHeight": 0.01}, "Tubes": [{"Length": 1.0}]}
        self.assertEqual(validate_arc_tubes(tubes, tubes["Tubes"], rows), 0)
        tubes["ArcTubes"] = {"Count": 0}
        with self.assertRaisesRegex(ValueError, "without arc rows"):
            validate_arc_tubes(tubes, tubes["Tubes"], rows)

    def test_negatives(self):
        rows = tagged_loop(135.0, 27)
        tubes, tube_rows = arc_census(rows, 2, 135.0)

        def rejected(mutate, message):
            copied = copy.deepcopy(tubes)
            mutate(copied)
            with self.assertRaisesRegex(ValueError, message):
                validate_arc_tubes(copied, copied["Tubes"], rows)

        rejected(lambda c: c["Tubes"][0]["Arc"].__setitem__("Radius", 1.001), "does not follow its tagged arc")
        rejected(lambda c: c["Tubes"][0]["Arc"].__setitem__("Parts", 3), "does not follow its tagged arc")
        rejected(lambda c: c["Tubes"][0]["Arc"].__setitem__("Sign", -1), "does not follow its tagged arc")
        rejected(lambda c: c["Tubes"][0]["Arc"].__setitem__("ArcId", 7), "does not carry")
        rejected(lambda c: c["Tubes"][0].__setitem__("Length", 5.0), "does not follow its tagged arc")
        rejected(lambda c: c["Section"].__setitem__("Radius", 0.3), "does not follow its tagged arc")
        rejected(lambda c: c["ArcTubes"].__setitem__("SharedSections", 7), "arc summary does not match")
        rejected(lambda c: c["ArcTubes"].__setitem__("PartSplits", 0), "arc summary does not match")
        rejected(lambda c: c["ArcTubes"].__setitem__("Count", 1), "arc summary does not match")
        rejected(lambda c: c["Tubes"][2]["Joints"][0].__setitem__("PlaneCut", True), "joint record")
        rejected(lambda c: c.pop("ArcTubes"), "arc summary does not match")


class ArcFaceEndCrossingSlopeTest(unittest.TestCase):
    """Mesher design round 3 class (9) 9H (part M 3.3 Fact 2; decisions 491 / 510 O9): an arc row's
    face-end formulas read the arc's crossing slope max_u |s'(u)| = rho h / (r sqrt(r^2 - h^2)) at the
    inner node circle r = rho - (R + h_pyr) in place of |tan theta|; the record's CrossingSlope is
    bound to it, recomputed from the row's Arc and the census CouponBox."""

    def test_slope_formula(self):
        for theta in (2.1, 45.0, 74.3):
            for rho in (1.56, 58.5):
                h = rho * math.sin(math.radians(theta))
                self.assertAlmostEqual(arc_crossing_slope(rho, 0.0, h) / math.tan(math.radians(theta)), 1.0, places=12)
                self.assertGreater(arc_crossing_slope(rho, 0.04, h), math.tan(math.radians(theta)))
                self.assertAlmostEqual(arc_crossing_slope(1e9 * rho, 0.04, 1e9 * h) / math.tan(math.radians(theta)), 1.0, places=8)
        # the mesher's value on the 32dc558f4810 geometry (rho 58.5, h 56.324, envelope 0.04)
        self.assertAlmostEqual(arc_crossing_slope(58.5, 0.04, 56.324169260983396), 3.5997093489397574, places=12)
        with self.assertRaisesRegex(ValueError, "does not reach the face plane"):
            arc_crossing_slope(1.0, 0.04, 0.97)
        with self.assertRaisesRegex(ValueError, "envelope"):
            arc_crossing_slope(1.0, 1.0, 0.5)

    def face_end_row(self, slope_scale=1.0, drop_slope=False, theta=45.0):
        radius, h_pyr, lc = 0.03175, 0.0079375, 0.05
        rho, centre = 1.5 / math.sin(math.radians(theta)), [0.0, 0.0]
        box = ([-2.0, -1.0, -1.0], [1.5, 3.0, 1.0])      # the x1 face at 1.5 = rho sin theta from the centre
        slope = arc_crossing_slope(rho, radius + h_pyr, 1.5)
        lc_end = max(lc, 4.0 * h_pyr * slope)
        m = max(1, math.ceil(2.0 * (radius + h_pyr) * slope / lc_end * (1.0 - 1e-9)))
        shear = (radius + h_pyr) * slope
        record = {"Face": "x1", "End": "end", "ThetaDegrees": theta, "Layers": m, "EndSpacing": lc_end,
                  "EnvelopeShear": shear, "LayerThicknessRange": [lc_end - shear / m, lc_end + shear / m],
                  "OverLength": shear + lc, "Kappa": [0.0, 0.0], "CrossingSlope": slope * slope_scale}
        if drop_slope:
            record.pop("CrossingSlope")
        row = {"Arc": {"ArcId": 1, "Centre": centre, "Radius": rho, "Sign": 1, "Part": 1, "Parts": 1},
               "FaceEnds": [record]}
        return row, lc, {"Radius": radius, "PyramidHeight": h_pyr}, box, lc_end, slope

    def test_arc_row_binds_the_recomputed_slope(self):
        row, lc, section, box, lc_end, slope = self.face_end_row()
        self.assertEqual(validate_tube_face_ends(row, lc, section, box=box), lc_end)
        self.assertGreater(slope, math.tan(math.radians(45.0)))
        with self.assertRaisesRegex(ValueError, "CrossingSlope does not follow"):
            validate_tube_face_ends(self.face_end_row(slope_scale=1.0 + 1e-9)[0], lc, section, box=box)
        with self.assertRaisesRegex(ValueError, "lacks its CrossingSlope"):
            validate_tube_face_ends(self.face_end_row(drop_slope=True)[0], lc, section, box=box)
        with self.assertRaisesRegex(ValueError, "needs the census CouponBox"):
            validate_tube_face_ends(row, lc, section)
        # the straight law (tan theta) is NOT the arc's: at 74.3 degrees (the 32dc558f4810 tilt) the
        # arc record's spacing / layers do not follow a straight row's formulas
        steep, lc, section, box, lc_end, slope = self.face_end_row(theta=74.3)
        self.assertEqual(validate_tube_face_ends(steep, lc, section, box=box), lc_end)
        self.assertGreater(slope, 1.5 * math.tan(math.radians(74.3)))
        steep.pop("Arc")
        steep["FaceEnds"][0]["CrossingSlope"] = math.tan(math.radians(74.3))
        with self.assertRaisesRegex(ValueError, "face-end rule"):
            validate_tube_face_ends(steep, lc, section, box=box)

    def test_straight_row_slope_is_tan_theta(self):
        radius, h_pyr, lc, theta = 0.03175, 0.0079375, 0.05, 45.0
        slope = math.tan(math.radians(theta))
        lc_end = max(lc, 4.0 * h_pyr * slope)
        m = max(1, math.ceil(2.0 * (radius + h_pyr) * slope / lc_end * (1.0 - 1e-9)))
        shear = (radius + h_pyr) * slope
        record = {"Face": "x1", "End": "end", "ThetaDegrees": theta, "Layers": m, "EndSpacing": lc_end,
                  "EnvelopeShear": shear, "LayerThicknessRange": [lc_end - shear / m, lc_end + shear / m],
                  "OverLength": shear + lc, "Kappa": [-slope, 0.0], "CrossingSlope": slope}
        section = {"Radius": radius, "PyramidHeight": h_pyr}
        self.assertEqual(validate_tube_face_ends({"FaceEnds": [record]}, lc, section), lc_end)
        record.pop("CrossingSlope")          # a pre-9H straight record carries none
        self.assertEqual(validate_tube_face_ends({"FaceEnds": [record]}, lc, section), lc_end)
        record["CrossingSlope"] = 1.01 * slope
        with self.assertRaisesRegex(ValueError, "CrossingSlope does not follow"):
            validate_tube_face_ends({"FaceEnds": [record]}, lc, section)


class FaceEndDerivedApexTest(unittest.TestCase):
    """Mesher design round 3 class (8) 8A-bitwise (part M 5.2; decisions 510 O8 / 563): above the kind's largest
    built face-end tilt the block's pyramids take h_pyr_end = min(h_pyr, TangentialSize / (4 s)); the
    record's PyramidHeight is bound, and the apex / regime formulas read it. Production thin section:
    EdgeSize 2 nm, 5 rings -> R 62 nm, h_pyr 16 nm; fabricated: R 31.75 nm, h_pyr 8 nm; lc 50 nm."""

    def straight_record(self, theta, radius, h_pyr, lc, derived, drop=False):
        slope = math.tan(math.radians(theta))
        apex_height = min(h_pyr, lc / (4.0 * slope)) if derived else h_pyr
        lc_end = max(lc, 4.0 * apex_height * slope)
        m = max(1, math.ceil(2.0 * (radius + h_pyr) * slope / lc_end * (1.0 - 1e-9)))
        shear = (radius + h_pyr) * slope
        record = {"Face": "x1", "End": "end", "ThetaDegrees": theta, "Layers": m, "EndSpacing": lc_end,
                  "EnvelopeShear": shear, "LayerThicknessRange": [lc_end - shear / m, lc_end + shear / m],
                  "OverLength": shear + lc, "Kappa": [-slope, 0.0], "CrossingSlope": slope, "PyramidHeight": apex_height}
        if drop:
            record.pop("PyramidHeight")
        return {"FaceEnds": [record]}, lc_end, apex_height

    def test_thin_above_70_takes_the_derived_height(self):
        radius, h_pyr, lc = 0.062, 0.016, 0.05
        section = {"Radius": radius, "PyramidHeight": h_pyr}
        self.assertEqual(FACE_END_DERIVED_APEX_ABOVE_DEGREES, {"fabricated": 70.0, "thin": 70.0})
        # thin 74.3 degrees (the 32dc558f4810 crossing): h_pyr_end 3.49 nm, lc_end 50 nm, m 12 (E3's block
        # read lc_end 247 nm, m 3)
        row, lc_end, apex_height = self.straight_record(74.3, radius, h_pyr, lc, derived=True)
        self.assertEqual(validate_tube_face_ends(row, lc, section, fabricated=False), lc_end)
        self.assertEqual(lc_end, lc)
        self.assertAlmostEqual(apex_height, lc / (4.0 * math.tan(math.radians(74.3))))
        self.assertEqual(row["FaceEnds"][0]["Layers"], 12)
        # the A2 (4) block (h_pyr) is refused there ...
        with self.assertRaisesRegex(ValueError, "PyramidHeight does not follow"):
            validate_tube_face_ends(self.straight_record(74.3, radius, h_pyr, lc, derived=False)[0], lc, section,
                                    fabricated=False)
        with self.assertRaisesRegex(ValueError, "lacks its PyramidHeight"):
            validate_tube_face_ends(self.straight_record(74.3, radius, h_pyr, lc, derived=True, drop=True)[0], lc,
                                    section, fabricated=False)
        # ... and at 70 degrees (the V10 thin 70, the thin admission top), 45 (the stored O4 thin) and
        # below, the A2 (4) block holds, with or without the field (decision 563: no built block changes)
        for theta in (70.0, 45.0, 20.0):
            row, lc_end, apex_height = self.straight_record(theta, radius, h_pyr, lc, derived=False)
            self.assertEqual(apex_height, h_pyr)
            self.assertEqual(validate_tube_face_ends(row, lc, section, fabricated=False), lc_end)
            row, lc_end, _ = self.straight_record(theta, radius, h_pyr, lc, derived=False, drop=True)
            self.assertEqual(validate_tube_face_ends(row, lc, section, fabricated=False), lc_end)
        # the FABRICATED kind keeps the A2 (4) block at 70 and derives above it
        fab = {"Radius": 0.03175, "PyramidHeight": 0.008}
        row, lc_end, _ = self.straight_record(70.0, 0.03175, 0.008, lc, derived=False)
        self.assertEqual(validate_tube_face_ends(row, lc, fab, fabricated=True), lc_end)
        row, lc_end, apex_height = self.straight_record(75.5, 0.03175, 0.008, lc, derived=True)
        self.assertEqual(validate_tube_face_ends(row, lc, fab, fabricated=True), lc_end)
        self.assertEqual(row["FaceEnds"][0]["Layers"], 7)
        self.assertLess(apex_height, 0.008)
        # without the kind no derivation is assumed (a record's PyramidHeight = h_pyr is still bound)
        row, lc_end, _ = self.straight_record(75.5, 0.03175, 0.008, lc, derived=False)
        self.assertEqual(validate_tube_face_ends(row, lc, fab), lc_end)


class ArcScopeGuardsTest(unittest.TestCase):
    """Decision 391 MAJOR-2 (ii): the arc face-end and arc-joint-tilt guards are in the contract's
    scope list, with the mesher's list spelled identically (ids, order, detection origin) and the
    loop end's tested turn bound."""

    def julia_guards(self):
        source = (HERE / "mesh_spatial_coupon.jl").read_text()
        block = re.search(r"const RECIPE_SCOPE_GUARDS = \[(.*?)\n(?:const|\S)", source, re.S).group(1)
        return [(m.group(1), m.group(2)) for m in re.finditer(r'\("([A-Za-z]+)", "(inputs|build)",', block)]

    def test_guards_present_and_parsed(self):
        # Round 3 B2 (decision 510): the interim guards ArcArcJoint (class 11) and CollarFaceEnd (class 5)
        # join the list; SteepFaceCrossing also names the thin tested-range rule (class 8, 8B).
        for guard in ("ArcFaceEnds", "ArcJointTilt", "ArcArcJoint", "CollarFaceEnd", "SteepFaceCrossing"):
            self.assertEqual(RECIPE_SCOPE_GUARDS[guard], "build")
            self.assertEqual(scope_guard_in_text(f"ERROR: ScopeGuard[{guard}]: an arc ...; arc 1 part 1"), guard)
        self.assertEqual(ARC_JOINT_TURN_BOUND_RADIANS, 1.6e-6)

    def test_lists_identical_with_the_mesher(self):
        julia = self.julia_guards()
        self.assertEqual([guard for guard, _ in julia], list(RECIPE_SCOPE_GUARDS))
        self.assertEqual(dict(julia), dict(RECIPE_SCOPE_GUARDS))
        source = (HERE / "mesh_spatial_coupon.jl").read_text()
        bound = re.search(r"const ARC_JOINT_TURN_BOUND = ([0-9.e+-]+)", source).group(1)
        self.assertEqual(float(bound), ARC_JOINT_TURN_BOUND_RADIANS)
        # Round 2b (decision 437 (3)): the tested ranges that lift the guards, parsed from the mesher.
        smooth = re.search(r"const ARC_SMOOTH_JOINT_TURN_BOUND = ([0-9.e+-]+)", source).group(1)
        self.assertEqual(float(smooth), ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS)
        # Round 3 class (6) (decision 556): the corner range is the BUILT range PER KIND, one named
        # constant with a fabricated and a thin field in both lists.
        corner = re.search(r"const ARC_CORNER_JOINT_TURN_RANGE = \(fabricated=\(([0-9.e+-]+), deg2rad\(([0-9.]+)\)\), "
                           r"thin=\(([0-9.e+-]+), deg2rad\(([0-9.]+)\)\)\)", source)
        self.assertIsNotNone(corner)
        self.assertEqual(float(corner.group(1)), ARC_CORNER_JOINT_TURN_RANGE_RADIANS["fabricated"][0])
        self.assertAlmostEqual(math.radians(float(corner.group(2))), ARC_CORNER_JOINT_TURN_RANGE_RADIANS["fabricated"][1])
        self.assertEqual(float(corner.group(3)), ARC_CORNER_JOINT_TURN_RANGE_RADIANS["thin"][0])
        self.assertAlmostEqual(math.radians(float(corner.group(4))), ARC_CORNER_JOINT_TURN_RANGE_RADIANS["thin"][1])
        self.assertEqual(set(ARC_CORNER_JOINT_TURN_RANGE_RADIANS), {"fabricated", "thin"})
        self.assertEqual((float(corner.group(2)), float(corner.group(4))), (15.0, 30.0))
        self.assertLess(ARC_CORNER_JOINT_TURN_RANGE_RADIANS["fabricated"][1], ARC_CORNER_JOINT_TURN_RANGE_RADIANS["thin"][1])
        # Round 3 (decisions 497 / 510 MAJOR-1 / O6): ONE range constant, (lowest built, largest built), in
        # both lists; the thin tested-range bound of 8B likewise.
        tilt = re.search(r"const ARC_FACE_END_TILT_RANGE = \(deg2rad\(([0-9.]+)\), deg2rad\(([0-9.]+)\)\)", source)
        self.assertEqual((float(tilt.group(1)), float(tilt.group(2))), ARC_FACE_END_TILT_RANGE_DEGREES)
        self.assertEqual(ARC_FACE_END_TILT_RANGE_DEGREES, (0.1, 70.0))
        thin = re.search(r"const THIN_FACE_END_TILT_BOUND = deg2rad\(([0-9.]+)\)", source).group(1)
        self.assertEqual(float(thin), THIN_FACE_END_TILT_BOUND_DEGREES)
        # Round 3 8A-bitwise (decision 510 O8): the per-kind thresholds of the derived face-end pyramid height.
        # Round 3 B4 (decision 566): the fabricated tested-range bound, the kind's largest built straight tilt.
        fab_bound = re.search(r"const FABRICATED_FACE_END_TILT_BOUND = deg2rad\(([0-9.]+)\)", source).group(1)
        self.assertEqual(float(fab_bound), FABRICATED_FACE_END_TILT_BOUND_DEGREES)
        self.assertEqual(FABRICATED_FACE_END_TILT_BOUND_DEGREES, 70.0)
        self.assertEqual(FACE_END_DERIVED_APEX_ABOVE_DEGREES["fabricated"], FABRICATED_FACE_END_TILT_BOUND_DEGREES)
        # (decision 563: each kind's largest BUILT tilt; the thin one IS the thin admission top)
        apex = re.search(r"const FACE_END_DERIVED_APEX_ABOVE = \(fabricated=deg2rad\(([0-9.]+)\), thin=THIN_FACE_END_TILT_BOUND\)",
                         source)
        self.assertIsNotNone(apex)
        self.assertEqual({"fabricated": float(apex.group(1)), "thin": float(thin)}, FACE_END_DERIVED_APEX_ABOVE_DEGREES)
        self.assertEqual(FACE_END_DERIVED_APEX_ABOVE_DEGREES["fabricated"], ARC_FACE_END_TILT_RANGE_DEGREES[1])
        self.assertEqual(FACE_END_DERIVED_APEX_ABOVE_DEGREES["thin"], THIN_FACE_END_TILT_BOUND_DEGREES)
        self.assertLessEqual(ARC_FACE_END_TILT_RANGE_DEGREES[1], THIN_FACE_END_TILT_BOUND_DEGREES)
        self.assertLess(ARC_JOINT_TURN_BOUND_RADIANS, ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS)
        for kind in ("fabricated", "thin"):
            self.assertLess(ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS, ARC_CORNER_JOINT_TURN_RANGE_RADIANS[kind][0])


if __name__ == "__main__":
    unittest.main()
