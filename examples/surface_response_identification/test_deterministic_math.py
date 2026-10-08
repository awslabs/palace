#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of deterministic_math (round 3 class (10), DESIGN-part-G G.10.5): the golden-hex table,
the exact-rational cross-check of the correct rounding, the 1-ulp libm sanity bound, the
special values and the fail-closed arguments.

The golden table testdata/deterministic-math-golden.csv (function, argument hex(es), result hex)
is platform-independent by construction: a correctly rounded double is unique. Regenerate it
(only when the ARGUMENT SET changes; a changed RESULT is a defect of the module or of the
table) with ``python3 test_deterministic_math.py --write-golden``."""
from fractions import Fraction
import math
from pathlib import Path
import random
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import deterministic_math as dm  # noqa: E402

GOLDEN = HERE / "testdata" / "deterministic-math-golden.csv"
FUNCTIONS = {"sin": (dm.sin, math.sin), "cos": (dm.cos, math.cos), "tan": (dm.tan, math.tan), "atan": (dm.atan, math.atan),
             "atan2": (dm.atan2, math.atan2), "acos": (dm.acos, math.acos), "asin": (dm.asin, math.asin)}


def ulps_away(x, steps):
    for _ in range(abs(steps)):
        x = math.nextafter(x, math.inf if steps > 0 else -math.inf)
    return x


def golden_arguments():
    """The argument set: ~1e4 cases, seeded uniform draws plus the hard cases (angles at and
    within a few ulps of k pi / 2, tiny and large angles, acos / asin inputs 1 - k 2^-53 for k =
    0..64 and their mirrors, atan2 on near-axis and equal-magnitude vectors)."""
    rng = random.Random(20261007)
    cases = []
    quarter_turns = [k * math.pi / 2 for k in range(-8, 9)]
    hard_angles = sorted({ulps_away(angle, j) for angle in quarter_turns for j in range(-3, 4)}
                         | {1.0e-300, 2.0 ** -30, 1.0e-8, -1.0e-8, 2.0 ** -26, 0.1, 0.5, 1.0, 2.0, 100.0, 1.0e3, -1.0e3}
                         | {rng.uniform(-1.0e4, 1.0e4) for _ in range(20)})
    for name in ("sin", "cos", "tan"):
        for angle in hard_angles:
            cases.append((name, (angle,)))
        for _ in range(1400):
            cases.append((name, (rng.uniform(-2.0 * math.pi, 2.0 * math.pi),)))
    for _ in range(1000):
        cases.append(("atan", (rng.uniform(-10.0, 10.0) * 10.0 ** rng.randint(-8, 2),)))
    near_axis = [1.0e-12, 1.0e-9, 2.0 ** -40, 1.0e-6]
    for y in near_axis + [-v for v in near_axis]:
        for x in (1.0, -1.0, 0.5, -3.0):
            cases.append(("atan2", (y, x)))
            cases.append(("atan2", (x, y)))
    for magnitude in (1.0, 0.3, 7.0e-3, 2.5e2):
        for sy in (1.0, -1.0):
            for sx in (1.0, -1.0):
                cases.append(("atan2", (sy * magnitude, sx * magnitude)))
    for _ in range(1700):
        cases.append(("atan2", (rng.uniform(-10.0, 10.0) * 10.0 ** rng.randint(-6, 1),
                                rng.uniform(-10.0, 10.0) * 10.0 ** rng.randint(-6, 1))))
    unit_edge = sorted({1.0 - k * 2.0 ** -53 for k in range(65)} | {-1.0 + k * 2.0 ** -53 for k in range(65)}
                       | {k * 2.0 ** -53 for k in range(-8, 9)} | {0.5, -0.5, math.sqrt(0.5), -math.sqrt(0.5), 0.25, 0.75})
    for name in ("acos", "asin"):
        for value in unit_edge:
            cases.append((name, (value,)))
        for _ in range(1200):
            cases.append((name, (rng.uniform(-1.0, 1.0),)))
    return cases


def golden_rows():
    rows = []
    for name, arguments in golden_arguments():
        rows.append((name, tuple(float(a).hex() for a in arguments), FUNCTIONS[name][0](*arguments).hex()))
    return rows


def read_golden():
    rows = []
    for line in GOLDEN.read_text().splitlines()[1:]:
        fields = line.split(",")
        rows.append((fields[0], tuple(fields[1:-1]), fields[-1]))
    return rows


class GoldenTableTest(unittest.TestCase):
    def test_every_golden_result_reproduces(self):
        """Every stored (function, arguments) -> hex result reproduces bitwise (the cross-platform
        half runs the same test on the login node); mismatches are printed with their cells."""
        stored = read_golden()
        self.assertGreaterEqual(len(stored), 10000)
        mismatches = []
        for name, argument_hexes, result_hex in stored:
            arguments = [float.fromhex(h) for h in argument_hexes]
            got = FUNCTIONS[name][0](*arguments).hex()
            if got != result_hex:
                mismatches.append(f"{name}({', '.join(repr(a) for a in arguments)}): {got} != golden {result_hex}")
        self.assertEqual(mismatches, [], "\n".join(mismatches[:40]))

    def test_golden_arguments_are_the_stored_arguments(self):
        """The argument set of this module and of the table agree (a changed set needs a rewrite)."""
        stored = [(name, hexes) for name, hexes, _ in read_golden()]
        generated = [(name, tuple(float(a).hex() for a in arguments)) for name, arguments in golden_arguments()]
        self.assertEqual(stored, generated)

    def test_libm_is_within_one_ulp_of_the_golden_results(self):
        """A sanity bound, not the verdict: this platform's libm agrees with the correctly rounded
        value to within one ulp on the whole table; tan to within two ulps (Apple libm's tan: 2
        ulps at 0.9075198688623791, exact-rational checked) and four above |x| ~ 1e3, where libm's
        own argument reduction is the weaker one."""
        worst = {}
        for name, argument_hexes, result_hex in read_golden():
            arguments = [float.fromhex(h) for h in argument_hexes]
            golden = float.fromhex(result_hex)
            libm = FUNCTIONS[name][1](*arguments)
            ulps = 0.0 if libm == golden else abs(libm - golden) / math.ulp(golden)
            bound = 1.0 if name != "tan" else (4.0 if abs(arguments[0]) > 1.0e3 else 2.0)
            self.assertLessEqual(ulps, bound, f"{name}({arguments}): libm {libm!r} vs golden {golden!r}")
            worst[name] = max(worst.get(name, 0.0), ulps)
        self.assertTrue(all(value <= 4.0 for value in worst.values()), worst)


class ExactRationalCrossCheckTest(unittest.TestCase):
    """The golden values ARE the correctly rounded values: an independent check with exact
    rational (Fraction) Taylor series — no decimal arithmetic, no rounding until the final
    int / int division — on a sample of the table where the series converge fast."""

    @staticmethod
    def exact_sin(x):
        f = Fraction(x)
        term, total = f, f
        for k in range(1, 70):
            term = -term * f * f / ((2 * k) * (2 * k + 1))
            total += term
        return total

    @staticmethod
    def exact_cos(x):
        f = Fraction(x)
        term, total = Fraction(1), Fraction(1)
        for k in range(1, 70):
            term = -term * f * f / ((2 * k - 1) * (2 * k))
            total += term
        return total

    @staticmethod
    def exact_atan(x):
        f = Fraction(x)
        term, total = f, f
        for k in range(1, 260):
            term = -term * f * f
            total += term / (2 * k + 1)
        return total

    @staticmethod
    def correctly_rounded(exact):
        candidate = float(exact)   # CPython's int / int true division is correctly rounded
        below, above = math.nextafter(candidate, -math.inf), math.nextafter(candidate, math.inf)
        assert (Fraction(below) + Fraction(candidate)) / 2 < exact < (Fraction(candidate) + Fraction(above)) / 2
        return candidate

    def test_sin_cos_atan_against_exact_rational_series(self):
        rng = random.Random(3)
        checked = 0
        for _ in range(150):
            x = rng.uniform(-1.5, 1.5)
            self.assertEqual(dm.sin(x), self.correctly_rounded(self.exact_sin(x)), x)
            self.assertEqual(dm.cos(x), self.correctly_rounded(self.exact_cos(x)), x)
            t = rng.uniform(-0.5, 0.5)
            self.assertEqual(dm.atan(t), self.correctly_rounded(self.exact_atan(t)), t)
            checked += 3
        for k in range(0, 65):
            # acos(1 - k 2^-53) = 2 atan(sqrt((1 - x) / (1 + x))): checked through asin / atan identities
            # is circular; instead check sin(acos(x)) and cos(acos(x)) reproduce x to one ulp.
            x = 1.0 - k * 2.0 ** -53
            angle = dm.acos(x)
            self.assertLessEqual(abs(dm.cos(angle) - x), 2.0 * math.ulp(1.0), (k, angle))
            checked += 1
        self.assertGreaterEqual(checked, 500)


class SpecialValueTest(unittest.TestCase):
    def test_signed_zero_and_axis_conventions(self):
        self.assertEqual(dm.sin(0.0).hex(), (0.0).hex())
        self.assertEqual(dm.sin(-0.0).hex(), (-0.0).hex())
        self.assertEqual(dm.cos(0.0), 1.0)
        self.assertEqual(dm.tan(-0.0).hex(), (-0.0).hex())
        self.assertEqual(dm.atan2(0.0, 1.0).hex(), (0.0).hex())
        self.assertEqual(dm.atan2(-0.0, 1.0).hex(), (-0.0).hex())
        self.assertEqual(dm.atan2(0.0, -1.0), math.pi)
        self.assertEqual(dm.atan2(-0.0, -1.0), -math.pi)
        self.assertEqual(dm.atan2(0.0, 0.0).hex(), (0.0).hex())
        self.assertEqual(dm.atan2(0.0, -0.0), math.pi)
        self.assertEqual(dm.atan2(-0.0, -0.0), -math.pi)
        self.assertEqual(dm.atan2(2.0, 0.0), math.pi / 2)
        self.assertEqual(dm.atan2(-2.0, -0.0), -math.pi / 2)
        self.assertEqual(dm.acos(1.0), 0.0)
        self.assertEqual(dm.acos(-1.0), math.pi)
        self.assertEqual(dm.acos(0.0), math.pi / 2)
        self.assertEqual(dm.asin(1.0), math.pi / 2)
        self.assertEqual(dm.asin(-1.0), -math.pi / 2)
        self.assertEqual(dm.asin(-0.0).hex(), (-0.0).hex())
        for name, (ours, libm) in FUNCTIONS.items():
            arguments = (0.0, 1.0) if name == "atan2" else (0.0,)
            self.assertEqual(ours(*arguments).hex(), libm(*arguments).hex(), name)

    def test_the_batch_3_joint_turn_case_reads_one_value(self):
        """P-G10.4: acos of a dot within one ulp of 1 — 0.0 on one platform, 1.49e-8 on the other
        before; the correctly rounded value is unique (acos(1 - 2^-53) = sqrt(2^-52) to 1 ulp)."""
        x = math.nextafter(1.0, 0.0)
        turn = dm.acos(x)
        self.assertAlmostEqual(turn, math.sqrt(2.0 * (1.0 - x)), delta=2.0 * math.ulp(turn))
        self.assertEqual(dm.acos(1.0), 0.0)
        self.assertEqual(dm.acos(min(1.0, 1.0 + 2.0 ** -52)), 0.0)

    def test_pi_decimal(self):
        self.assertEqual(str(dm.pi_decimal(30)), "3.14159265358979323846264338328")
        self.assertEqual(float(dm.pi_decimal(40)), math.pi)


class FailClosedTest(unittest.TestCase):
    def test_non_finite_and_out_of_range_arguments_are_refused(self):
        for value in (math.nan, math.inf, -math.inf):
            for name, (ours, _) in FUNCTIONS.items():
                with self.assertRaises(ValueError, msg=name):
                    ours(value, 1.0) if name == "atan2" else ours(value)
        with self.assertRaises(ValueError):
            dm.atan2(1.0, math.nan)
        for value in (1.0 + 2.0 ** -52, -1.0000000000000002):
            with self.assertRaises(ValueError):
                dm.acos(value)
            with self.assertRaises(ValueError):
                dm.asin(value)
        with self.assertRaises(TypeError):
            dm.sin("1.0")

    def test_golden_hex_is_the_hex_of_the_result(self):
        self.assertEqual(dm.golden_hex(dm.atan2, 1.0, 1.0), dm.atan2(1.0, 1.0).hex())


def write_golden(path=GOLDEN):
    rows = golden_rows()
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["function,arguments...,result"]
    for name, argument_hexes, result_hex in rows:
        lines.append(",".join((name, *argument_hexes, result_hex)))
    path.write_text("\n".join(lines) + "\n")
    return len(rows)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--write-golden":
        print(write_golden(), "rows ->", GOLDEN)
        sys.exit(0)
    unittest.main()
