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
                                 RECIPE_SCOPE_GUARDS, THIN_FACE_END_TILT_BOUND_DEGREES, arc_part_count,
                                 boundary_arc_arc_joints, boundary_arc_runs, metal_loop_arc_parts,
                                 metal_loop_side_points, read_csv_rows, scope_classes, scope_guard_in_text,
                                 validate_arc_tubes)
import json

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
    def test_arc_chain_column_is_optional_for_every_reader(self):
        # Round 3 R7 (decision 510 MINOR-6): a tagged boundary WITHOUT the ArcChain column (every
        # stored arc coupon) reads Chain 0; one WITH it reads the chain; the contract reader
        # (semantic_mesh_contract.boundary_arc_tags) likewise; a straight boundary has no column.
        import semantic_mesh_contract
        without = tagged_loop(90.0, 18)
        self.assertNotIn("ArcChain", without[0])
        self.assertEqual([run["Chain"] for run in boundary_arc_runs(without)[0]], [0])
        tags, _ = semantic_mesh_contract.boundary_arc_tags(without)
        self.assertEqual({tag["ArcChain"] for tag in tags if tag is not None}, {0})
        with_chain = [dict(row, ArcChain="3" if row["ArcId"] else "") for row in without]
        self.assertEqual([run["Chain"] for run in boundary_arc_runs(with_chain)[0]], [3])
        tags, _ = semantic_mesh_contract.boundary_arc_tags(with_chain)
        self.assertEqual({tag["ArcChain"] for tag in tags if tag is not None}, {3})
        # The stored rounded-strip boundary (round 2b, no arc column at all) reads no runs.
        import csv
        with (HERE / "testdata" / "rounded-strip" / "plan-view-boundary.csv").open(newline="") as stream:
            stored = list(csv.DictReader(stream))
        self.assertNotIn("ArcId", stored[0])
        self.assertEqual(boundary_arc_runs(stored), [[]])

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


ARC_ARC = HERE / "testdata" / "arc-arc-joints"


def excerpt(name):
    """A reduced REAL excerpt of a stored build census (its Provenance names the full census's sha256)."""
    return json.loads((ARC_ARC / f"{name}-census.excerpt.json").read_text())


class ArcArcJointsTest(unittest.TestCase):
    """Round 3 class (11), fix 11 (decisions 491 / 510 / 577): the stage contract's arc-arc shared-section term.
    Under fix 11 two DISTINCT tagged arc runs of ONE circle meeting at a smooth joint share one tube section per
    placement (owned by the earlier tube, no JointEnd row), so the census identity reads SharedSections =
    JointEnds + PartSplits + TubesPerSide x ArcArcJoints; the joints are DERIVED from the bound boundary
    (consecutive tagged runs, JointSmooth 1, one circle within the arc-fit tolerance, one ArcSign; the ArcChain
    column, where present, must agree), so the STORED census of R1 (32b0083dad90, PBS 59700 on main
    931f8be05a: the census that stopped the per-entry verification at the pre-B3 identity, decision 577)
    validates with NO mesher change; an explicit ArcTubes.ArcArcJoints / ArcArcSections field is optional."""

    R1_JOINTS = [(1, 24, 14, 42), (1, 13, 23, 180), (1, 22, 12, 213), (1, 11, 21, 433)]

    def validate(self, census, boundary):
        tubes = census["PrismTubes"]
        return validate_arc_tubes(tubes, tubes["Tubes"], boundary, census["CouponBox"]["Radius"])

    def test_r1_stored_censuses_validate_with_the_derived_arc_arc_term(self):
        boundary = read_csv_rows(ARC_ARC / "r1-32b0083dad90-plan-view-boundary.csv")
        self.assertEqual(boundary[0].get("ArcChain"), "")
        self.assertEqual(boundary_arc_arc_joints(boundary), self.R1_JOINTS)      # the four claim / context chains
        for kind, expected in (("thin", (46, 40, 2, 1, 28)), ("fab", (92, 80, 4, 2, 56))):
            census = excerpt(f"r1-32b0083dad90-{kind}")
            summary = census["PrismTubes"]["ArcTubes"]
            shared, joint_ends, part_splits, per_side, arc_rows = expected
            self.assertEqual((summary["SharedSections"], summary["JointEnds"], summary["PartSplits"],
                              census["PrismTubes"]["Section"]["TubesPerSide"]), (shared, joint_ends, part_splits, per_side))
            self.assertNotIn("ArcArcJoints", summary)                           # the stored census predates the field
            self.assertEqual(self.validate(census, boundary), arc_rows)
            self.assertEqual(shared, joint_ends + part_splits + per_side * len(self.R1_JOINTS))
            # The pre-B3 identity (SharedSections = JointEnds + PartSplits) is what stopped R1: a census
            # written to it would now be refused - the term is derived, not assumed.
            narrowed = copy.deepcopy(census)
            narrowed["PrismTubes"]["ArcTubes"]["SharedSections"] = joint_ends + part_splits
            with self.assertRaisesRegex(ValueError, "arc summary does not match"):
                self.validate(narrowed, boundary)
        # Without the ArcChain column (a tagged boundary of the B2 era) the joints derive from the geometry alone.
        legacy = [{k: v for k, v in row.items() if k != "ArcChain"} for row in boundary]
        self.assertEqual(boundary_arc_arc_joints(legacy), self.R1_JOINTS)
        self.assertEqual(self.validate(excerpt("r1-32b0083dad90-thin"), legacy), 28)

    def test_b3_synthetics_s1_s3(self):
        # S1: the strip arc split at -55 degrees into ids 1 / 2 (one joint): fab 6 = 4 + 0 + 2 x 1, thin 3 = 2 + 0 + 1 x 1.
        # S3: three members (two joints), ids ordered (1, 2, 3) and permuted (3, 1, 2): fab 8 = 4 + 0 + 2 x 2 either way.
        for name, joints, shared in (("s1-fab", [(1, 1, 2, 9)], 6), ("s1-thin", [(1, 1, 2, 9)], 3),
                                     ("s3-fab", [(1, 1, 2, 8), (1, 2, 3, 14)], 8), ("s3p-fab", [(1, 3, 1, 8), (1, 1, 2, 14)], 8)):
            boundary = read_csv_rows(ARC_ARC / f"{name}-boundary.csv")
            self.assertNotIn("ArcChain", boundary[0])                           # the Julia fixtures write the 7 columns
            self.assertEqual(boundary_arc_arc_joints(boundary), joints, name)
            census = excerpt(name)
            self.assertEqual(census["PrismTubes"]["ArcTubes"]["SharedSections"], shared, name)
            self.assertEqual(self.validate(census, boundary), census["PrismTubes"]["ArcTubes"]["Count"], name)

    def test_optional_explicit_fields(self):
        boundary = read_csv_rows(ARC_ARC / "s1-fab-boundary.csv")
        census = excerpt("s1-fab")
        census["PrismTubes"]["ArcTubes"].update({"ArcArcJoints": 1, "ArcArcSections": 2})
        self.assertEqual(self.validate(census, boundary), 4)
        for name, wrong in (("ArcArcJoints", 2), ("ArcArcSections", 1)):
            bad = copy.deepcopy(census)
            bad["PrismTubes"]["ArcTubes"][name] = wrong
            with self.assertRaisesRegex(ValueError, f"{name} {wrong} differs from the .* derived"):
                self.validate(bad, boundary)

    def test_fail_closed(self):
        boundary = read_csv_rows(ARC_ARC / "r1-32b0083dad90-plan-view-boundary.csv")
        census = excerpt("r1-32b0083dad90-thin")
        joint_row = next(row for row in boundary if int(row["Vertex"]) == 42)
        self.assertEqual((joint_row["ArcId"], joint_row["ArcChain"], joint_row["JointSmooth"]), ("14", "4", "1"))
        # (a) ArcChain disagreeing with the geometry (one run of the pair un-chained).
        broken = copy.deepcopy(boundary)
        for row in broken:
            if row["ArcId"] == "14":
                row["ArcChain"] = "0"
        with self.assertRaisesRegex(ValueError, "ArcChain columns disagree .* cannot be derived"):
            boundary_arc_arc_joints(broken)
        # (b) a chain pair that is not one circle (the second run's radius off by 1e-3 relative).
        broken = copy.deepcopy(boundary)
        for row in broken:
            if row["ArcId"] == "14":
                row["ArcR"] = repr(float(row["ArcR"]) * (1.0 + 1.0e-3))
        with self.assertRaisesRegex(ValueError, "share ArcChain 4 .* not one circle .* cannot be derived"):
            boundary_arc_arc_joints(broken)
        # (c) the same circle but opposite ArcSign: not a joint by geometry, and the chain claims one -> fail closed.
        broken = copy.deepcopy(boundary)
        for row in broken:
            if row["ArcId"] == "14":
                row["ArcSign"] = "-1"
        with self.assertRaisesRegex(ValueError, "cannot be derived"):
            boundary_arc_arc_joints(broken)
        # (d) a joint vertex not tagged smooth is no shared section: the identity then wants 45, the census says 46.
        unsmooth = copy.deepcopy(boundary)
        next(row for row in unsmooth if int(row["Vertex"]) == 42)["JointSmooth"] = "0"
        self.assertEqual(len(boundary_arc_arc_joints(unsmooth)), 3)
        with self.assertRaisesRegex(ValueError, "arc summary does not match"):
            self.validate(census, unsmooth)
        # (e) a section without TubesPerSide cannot carry the term.
        no_side = copy.deepcopy(census)
        del no_side["PrismTubes"]["Section"]["TubesPerSide"]
        with self.assertRaisesRegex(ValueError, "lacks TubesPerSide"):
            self.validate(no_side, boundary)
        # (f) a derived joint whose arc is not built.
        unbuilt = copy.deepcopy(census)
        unbuilt["PrismTubes"]["Tubes"] = [row for row in unbuilt["PrismTubes"]["Tubes"]
                                          if not ("Arc" in row and row["Arc"]["ArcId"] == 14)]
        with self.assertRaisesRegex(ValueError, "not all built as arc tubes"):
            self.validate(unbuilt, boundary)

    def test_pre_b3_identity_unchanged_without_arc_arc_joints(self):
        # No consecutive runs (the one-arc fixtures): the term is 0 and the pre-B3 identity stands bitwise; two
        # consecutive runs whose joint vertex is NOT smooth-tagged, or of distinct circles, add nothing.
        rows = tagged_loop(135.0, 27)
        self.assertEqual(boundary_arc_arc_joints(rows), [])
        tubes, tube_rows = arc_census(rows, 2, 135.0)
        self.assertEqual(validate_arc_tubes(tubes, tube_rows, rows), 2)
        self.assertEqual(validate_arc_tubes(tubes, tube_rows, rows, 1.0), 2)
        boundary = read_csv_rows(ARC_ARC / "s1-fab-boundary.csv")
        corner = [dict(row, JointSmooth="0" if int(row["Vertex"]) == 9 else row["JointSmooth"]) for row in boundary]
        self.assertEqual(boundary_arc_arc_joints(corner), [])
        distinct = [dict(row, ArcR=repr(1.001) if row["ArcId"] == "2" else row["ArcR"]) for row in boundary]
        self.assertEqual(boundary_arc_arc_joints(distinct), [])


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
        corner = re.search(r"const ARC_CORNER_JOINT_TURN_RANGE = \(([0-9.e+-]+), deg2rad\(([0-9.]+)\)\)", source)
        self.assertEqual(float(corner.group(1)), ARC_CORNER_JOINT_TURN_RANGE_RADIANS[0])
        self.assertAlmostEqual(math.radians(float(corner.group(2))), ARC_CORNER_JOINT_TURN_RANGE_RADIANS[1])
        # Round 3 (decisions 497 / 510 MAJOR-1 / O6): ONE range constant, (lowest built, largest built), in
        # both lists; the thin tested-range bound of 8B likewise.
        tilt = re.search(r"const ARC_FACE_END_TILT_RANGE = \(deg2rad\(([0-9.]+)\), deg2rad\(([0-9.]+)\)\)", source)
        self.assertEqual((float(tilt.group(1)), float(tilt.group(2))), ARC_FACE_END_TILT_RANGE_DEGREES)
        self.assertEqual(ARC_FACE_END_TILT_RANGE_DEGREES, (0.1, 70.0))
        thin = re.search(r"const THIN_FACE_END_TILT_BOUND = deg2rad\(([0-9.]+)\)", source).group(1)
        self.assertEqual(float(thin), THIN_FACE_END_TILT_BOUND_DEGREES)
        self.assertLessEqual(ARC_FACE_END_TILT_RANGE_DEGREES[1], THIN_FACE_END_TILT_BOUND_DEGREES)
        self.assertLess(ARC_JOINT_TURN_BOUND_RADIANS, ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS)
        self.assertLess(ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS, ARC_CORNER_JOINT_TURN_RANGE_RADIANS[0])


if __name__ == "__main__":
    unittest.main()
