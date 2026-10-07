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

from mesh_stage_contract import (ARC_CORNER_JOINT_TURN_RANGE_RADIANS, ARC_FACE_END_TILT_BOUND_DEGREES,
                                 ARC_JOINT_TURN_BOUND_RADIANS, ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS,
                                 RECIPE_SCOPE_GUARDS, arc_part_count, boundary_arc_runs,
                                 metal_loop_arc_parts, metal_loop_side_points, scope_classes, scope_guard_in_text,
                                 validate_arc_tubes)

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


class ArcScopeGuardsTest(unittest.TestCase):
    """Decision 391 MAJOR-2 (ii): the arc face-end and arc-joint-tilt guards are in the contract's
    scope list, with the mesher's list spelled identically (ids, order, detection origin) and the
    loop end's tested turn bound."""

    def julia_guards(self):
        source = (HERE / "mesh_spatial_coupon.jl").read_text()
        block = re.search(r"const RECIPE_SCOPE_GUARDS = \[(.*?)\n(?:const|\S)", source, re.S).group(1)
        return [(m.group(1), m.group(2)) for m in re.finditer(r'\("([A-Za-z]+)", "(inputs|build)",', block)]

    def test_guards_present_and_parsed(self):
        for guard in ("ArcFaceEnds", "ArcJointTilt"):
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
        corner = re.search(r"const ARC_CORNER_JOINT_TURN_RANGE = \(([0-9.e+-]+), deg2rad\(([0-9.]+)\)\)", source)
        self.assertEqual(float(corner.group(1)), ARC_CORNER_JOINT_TURN_RANGE_RADIANS[0])
        self.assertAlmostEqual(math.radians(float(corner.group(2))), ARC_CORNER_JOINT_TURN_RANGE_RADIANS[1])
        tilt = re.search(r"const ARC_FACE_END_TILT_BOUND = deg2rad\(([0-9.]+)\)", source).group(1)
        self.assertEqual(float(tilt), ARC_FACE_END_TILT_BOUND_DEGREES)
        self.assertLess(ARC_JOINT_TURN_BOUND_RADIANS, ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS)
        self.assertLess(ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS, ARC_CORNER_JOINT_TURN_RANGE_RADIANS[0])


if __name__ == "__main__":
    unittest.main()
