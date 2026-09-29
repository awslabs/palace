#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner family's angle-interpolation rule (corner-qualification block 2026-09-29):
the trace basis events (knot passages of fixed-layout vertices) and the stencil rule that
never straddles a knot-corner passage — the Python mirror of CornerBasisEvents /
SelectCornerFamilyStencil (cornertracebasis.cpp), pinned to the same numbers as the C++ unit
test test-cornerfamily.cpp."""

import importlib.util
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FAMILY = load("corner_family_interpolation")
GENERATOR = load("generate_corner_response")

# The qualified family's node set (both convexities): (angle, connectivity angle) per
# segment; the recorded legacy tie coupons at 90 / 135 / 180 are kept beside them.
CONVEX_NODES = [
    (75.0, 82.5), (80.0, 82.5), (85.0, 82.5), (90.0, 82.5),
    (90.0, 112.5), (105.0, 112.5), (120.0, 112.5), (135.0, 112.5),
    (135.0, 144.2), (144.0, 144.2), (150.0, 144.2), (153.434948822922, 144.2),
    (153.434948822922, 166.7), (159.0, 166.7), (165.0, 166.7), (180.0, 166.7),
]
CONCAVE_NODES = [
    (75.0, 82.5), (80.0, 82.5), (85.0, 82.5), (90.0, 82.5),
    (90.0, 112.5), (105.0, 112.5), (120.0, 112.5), (135.0, 112.5),
    (135.0, 146.6), (143.0, 146.6), (150.0, 146.6), (158.19859051364818, 146.6),
    (158.19859051364818, 169.1), (165.0, 169.1), (172.0, 169.1), (180.0, 169.1),
]


def indexed(nodes):
    return [(angle, connectivity, index) for index, (angle, connectivity) in enumerate(nodes)]


class CornerBasisEventsTest(unittest.TestCase):
    def test_events_of_the_rule(self):
        """Knot-corner passages in the family range [75, 180]: convex 90 (free 1 / 3 / 5 and
        metal 1 at the four corners: the fixed layout), 135 (the second crossing at (-R, R)),
        153.435 = 180 - atan(1/2) (free 2 at (-R, -R)); concave 90 (free 3 (R, R), metal 1
        (-R, -R)), 135 (free 2 (R, R) and the crossing), 158.199 = 180 - atan(2/5) (free 5 at
        (-R, R)). Side-midpoint passages (removed by the segment connectivity): convex
        141.340 = 180 - atan(4/5) (free 1 at (-R, 0)), concave 111.801 = 180 - atan(5/2)
        (free 5 at (0, R)); the anchor's crossing reaches (-R, 0) at 180."""
        convex = FAMILY.corner_event_angles("convex")
        concave = FAMILY.corner_event_angles("concave")
        for actual, expected in (
            ([a for a in convex if a >= 75.0], [90.0, 135.0, 180.0 - math.degrees(math.atan(0.5))]),
            ([a for a in concave if a >= 75.0], [90.0, 135.0, 180.0 - math.degrees(math.atan(0.4))]),
        ):
            self.assertEqual(len(actual), len(expected))
            for a, e in zip(actual, expected):
                self.assertAlmostEqual(a, e, places=9)
        self.assertAlmostEqual(convex[-1], 153.434948822922, places=9)
        self.assertAlmostEqual(concave[-1], 158.19859051364818, places=9)
        # (Midpoint passages at a knot-corner passage angle, e.g. the fixed layout at 90, are
        # absorbed by that event.)
        midpoints = {
            topology: [
                e for e in FAMILY.basis_events(topology)
                if not e["corner"] and 75.0 < e["angle_degrees"] < 180.0
                and all(abs(e["angle_degrees"] - c) > 1e-6 for c in FAMILY.corner_event_angles(topology))
            ]
            for topology in ("convex", "concave")
        }
        self.assertEqual([(e["role"], e["fraction"]) for e in midpoints["convex"]], [("free1", 0.0)])
        self.assertAlmostEqual(midpoints["convex"][0]["angle_degrees"], 180.0 - math.degrees(math.atan(0.8)), places=9)
        self.assertEqual([(e["role"], e["fraction"]) for e in midpoints["concave"]], [("free5", 0.75)])
        self.assertAlmostEqual(midpoints["concave"][0]["angle_degrees"], 180.0 - math.degrees(math.atan(2.5)), places=9)
        # Every event is a real passage: the knot of that role sits on the fixed vertex.
        for topology in ("convex", "concave"):
            for event in FAMILY.basis_events(topology):
                if not 75.0 <= event["angle_degrees"] <= 180.0:
                    continue
                layout = GENERATOR.metal_ring_layout(1.9, event["angle_degrees"], topology, 8)
                order = GENERATOR.ring_role_order(topology)
                fractions = {order[slot]: fraction for fraction, kind, slot in layout if kind != "slave"}
                self.assertAlmostEqual(fractions[event["role"]], event["fraction"], places=9, msg=str(event))

    def test_events_are_shared_with_the_cpp_rule(self):
        """The C++ unit test pins the same event angles (test-cornerfamily.cpp)."""
        source = (ROOT / "../../../test/unit/test-cornerfamily.cpp").resolve()
        self.assertTrue(source.is_file(), f"C++ unit test missing: {source}")
        text = source.read_text()
        for value in ("153.434948822922", "158.198590513648", "141.340191745910", "111.801409486352"):
            self.assertIn(value, text)


class StencilTest(unittest.TestCase):
    def test_qualified_family_is_cubic_in_every_segment(self):
        for topology, nodes in (("convex", CONVEX_NODES), ("concave", CONCAVE_NODES)):
            family = indexed(nodes)
            keys = sorted({connectivity for _, connectivity in nodes})
            for angle in (78.0, 82.5, 100.0, 112.5, 127.5, 142.5, 147.0, 172.5, 176.0):
                stencil = FAMILY.select_stencil(family, angle, topology)
                self.assertNotIn("reason", stencil, (topology, angle))
                self.assertEqual(stencil["rule"], "cubic", (topology, angle))
                self.assertAlmostEqual(sum(w for _, w in stencil["nodes"]), 1.0, places=12)
                # Every node of the stencil shares the segment's connectivity angle and lies on
                # the same side of every knot-corner passage as the angle.
                for index, _ in stencil["nodes"]:
                    self.assertEqual(family[index][1], stencil["connectivity_angle_degrees"])
                    for boundary in FAMILY.corner_event_angles(topology):
                        self.assertFalse(
                            min(angle, family[index][0]) + 1e-2 < boundary < max(angle, family[index][0]) - 1e-2,
                            (topology, angle, family[index], boundary),
                        )
                self.assertIn(stencil["connectivity_angle_degrees"], keys)
        # The recorded held-out 82.5 (convex): linear on 75 / 90 before the block, now cubic
        # on 75 / 80 / 85 / 90- of the [75, 90] segment.
        stencil = FAMILY.select_stencil(indexed(CONVEX_NODES), 82.5, "convex")
        self.assertEqual([CONVEX_NODES[i][0] for i, _ in stencil["nodes"]], [75.0, 80.0, 85.0, 90.0])
        self.assertEqual(stencil["connectivity_angle_degrees"], 82.5)

    def test_exact_node_preference(self):
        """At an event angle with two per-side coupons and a legacy tie coupon the legacy one
        is exact; without it the lower-angle segment's coupon."""
        family = indexed(CONVEX_NODES) + [(90.0, None, 99)]
        self.assertEqual(FAMILY.select_stencil(family, 90.0, "convex")["nodes"], [(99, 1.0)])
        stencil = FAMILY.select_stencil(indexed(CONVEX_NODES), 90.0, "convex")
        self.assertEqual(stencil["rule"], "exact")
        self.assertEqual(CONVEX_NODES[stencil["nodes"][0][0]], (90.0, 82.5))
        stencil = FAMILY.select_stencil(indexed(CONVEX_NODES), 135.0, "convex")
        self.assertEqual(CONVEX_NODES[stencil["nodes"][0][0]], (135.0, 112.5))

    def test_legacy_family_is_exact_only(self):
        legacy = [(a, None, i) for i, a in enumerate((75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0))]
        self.assertEqual(FAMILY.select_stencil(legacy, 120.0, "convex")["rule"], "exact")
        reason = FAMILY.select_stencil(legacy, 112.5, "convex")["reason"]
        self.assertIn("without segment connectivity records", reason)
        # A legacy coupon beside keyed segments does not extend them.
        mixed = indexed(CONVEX_NODES) + [(60.0, None, 99)]
        self.assertIn("without segment connectivity", FAMILY.select_stencil(mixed, 70.0, "convex")["reason"])

    def test_refusals_and_fail_closed(self):
        family = indexed(CONVEX_NODES)
        self.assertIn("sharper than", FAMILY.select_stencil(family, 70.0, "convex")["reason"])
        no_anchor = [n for n in family if n[0] < 180.0]
        self.assertIn("wider than", FAMILY.select_stencil(no_anchor, 175.0, "convex")["reason"])
        # A segment across a knot-corner passage (nodes 120 and 150 keyed at 112.5).
        with self.assertRaises(ValueError):
            FAMILY.select_stencil([(120.0, 112.5, 0), (150.0, 112.5, 1)], 130.0, "convex")
        # A connectivity angle on a passage.
        with self.assertRaises(ValueError):
            FAMILY.select_stencil([(120.0, 135.0, 0), (130.0, 135.0, 1)], 125.0, "convex")
        # A gap between segments (no coupon on the upper side of 135).
        gap = [n for n in family if n[1] != 144.2]
        self.assertIn("no segment", FAMILY.select_stencil(gap, 142.5, "convex")["reason"])
        # Lower orders where a segment has fewer nodes.
        short = [(135.0, 144.2, 0), (150.0, 144.2, 1), (153.434948822922, 144.2, 2)]
        self.assertEqual(FAMILY.select_stencil(short, 142.5, "convex")["rule"], "quadratic")
        self.assertEqual(FAMILY.select_stencil(short[:2], 142.5, "convex")["rule"], "linear")
        self.assertIn("no segment", FAMILY.select_stencil(short[:1] + [(75.0, 82.5, 5)], 100.0, "convex")["reason"])


if __name__ == "__main__":
    unittest.main()
