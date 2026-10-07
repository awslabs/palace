#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Correctly rounded transcendental functions for the content-hashed source path.

Round-3 class (10) (supervisor decisions 492 / 493 / 510, DESIGN-part-G G.10 Option B): the
case id of a coupon is the SHA-256 of its generated source files, which are written at full
double precision. IEEE 754 makes + - x / sqrt and the integer conversions correctly rounded,
hence identical on every platform; it does NOT cover sin / cos / tan / atan2 / acos / asin,
which glibc, Apple libm and numpy's SIMD loops round differently in a fraction of a percent
of arguments (the batch-2 / batch-3 case-id variants: <= 5.1e-14 relative in the arc chord
directions, 0.0 vs 1.49e-8 rad from acos of a dot within one ulp of 1). The correctly rounded
value of a function at a double argument is unique, so a function that returns it is the
same on every platform by definition — no quantum, no knife-edge.

Backend (decision 510 open choice O4): mpmath is NOT in either Python environment of record
(the local /opt/homebrew python 3.14.2 nor the cluster venv coupon-v2 python 3.13.12 of the
manifest's runtime digests), so this module evaluates with the standard-library ``decimal``
module: its own argument reduction and Taylor / arctangent series at a working precision of
WORKING_DIGITS significant digits, rounded ONCE to a double by exact rational comparison
(``float(Fraction(decimal))``: CPython's int / int true division is correctly rounded). A
Ziv-style guard (``_round_to_double``) re-evaluates at twice the precision whenever the
approximation lies within its error bound of a rounding boundary (the midpoint of two
adjacent doubles), so the returned double is the correctly rounded value for every argument
that is not an exact midpoint — impossible for a transcendental value of a nonzero rational
argument; the exactly representable results (sin 0, cos 0, atan2(0, x), acos 1, asin 0, ...)
are returned directly.

THE SCALAR-ARITHMETIC RULE (G.10.3) that this module is one half of: every float written into
a content-hashed source is produced by CPython scalar float arithmetic (+ - x / sqrt, int,
``math.hypot``: CPython's own correctly scaled implementation since 3.10, no libm) on
serialised inputs, and by this module for the transcendental functions; numpy arrays carry
data on that path but do no arithmetic on it (numpy's 2-vector ``dot`` / ``linalg.norm`` go
through BLAS: Accelerate contracts a0 b0 + a1 b1 into an FMA, OpenBLAS does not — the probe
of impl-B0/tools/probe_numpy_scalar_agreement.py); the frame rotation is applied with scalar
arithmetic. provenance.json records ``FloatSerialisation`` = {Rule FLOAT_SERIALISATION_RULE,
Module sha256 of this file}.

API (all arguments and results are Python floats; non-finite arguments fail closed):
    sin(x), cos(x), tan(x), atan2(y, x), atan(x), acos(x), asin(x)
    pi_decimal(digits): pi to ``digits`` significant digits (cached)
    golden_hex(function, *arguments) -> ``float.hex()`` of the result (the tests' key)
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import math

FLOAT_SERIALISATION_RULE = "deterministic-trig-v1"
# Working precision of the first evaluation (significant decimal digits). A double carries
# ~16 digits; the series below are truncated when a term falls below 10^-(prec + 2), and the
# accumulated rounding of a few hundred context operations stays far inside GUARD_DIGITS.
WORKING_DIGITS = 48
# The error bound claimed for an evaluation at precision ``prec``: 10^-(prec - GUARD_DIGITS)
# relative to the result (the series truncation and <= ~10^3 context roundings of 10^-prec
# each are covered by 10^8 of slack).
GUARD_DIGITS = 8
# Precision doublings allowed by the Ziv guard before failing closed (an exact midpoint).
MAX_PRECISION_DOUBLINGS = 3

_PI_CACHE = {}


def pi_decimal(digits):
    """pi to ``digits`` significant digits (the decimal module's documented recipe), cached."""
    if digits in _PI_CACHE:
        return _PI_CACHE[digits]
    with localcontext() as context:
        context.prec = digits + 4
        three = Decimal(3)
        lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            n, na = n + na, na + 8
            d, da = d + da, da + 32
            t = (t * n) / d
            s += t
    with localcontext() as context:
        context.prec = digits
        value = +s
    _PI_CACHE[digits] = value
    return value


def _check_finite(name, *values):
    for value in values:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise TypeError(f"deterministic_math.{name}: a real number is required, not {type(value).__name__}")
        if not math.isfinite(value):
            raise ValueError(f"deterministic_math.{name}: a finite argument is required, not {value!r}")


def _round_to_double(evaluate, name):
    """``evaluate(prec)`` returns the Decimal approximation of the exact result at ``prec``
    significant digits with relative error <= 10^-(prec - GUARD_DIGITS). Returns the
    correctly rounded double, re-evaluating at twice the precision while the approximation
    lies within its error bound of a midpoint between two adjacent doubles."""
    prec = WORKING_DIGITS
    for _ in range(MAX_PRECISION_DOUBLINGS + 1):
        approximation = Fraction(evaluate(prec))
        if approximation == 0:
            return 0.0
        error = abs(approximation) * Fraction(10) ** (-(prec - GUARD_DIGITS))
        candidate = float(approximation)              # correctly rounded int / int division
        if not math.isfinite(candidate):
            raise OverflowError(f"deterministic_math.{name}: the result overflows a double")
        below, above = math.nextafter(candidate, -math.inf), math.nextafter(candidate, math.inf)
        lower_midpoint = (Fraction(below) + Fraction(candidate)) / 2
        upper_midpoint = (Fraction(candidate) + Fraction(above)) / 2
        if (approximation - lower_midpoint > error) and (upper_midpoint - approximation > error):
            return candidate
        prec *= 2
    raise ArithmeticError(f"deterministic_math.{name}: the result lies on a rounding boundary at {prec // 2} digits "
                          "(an exact midpoint; not a transcendental value)")


def _reduce_half_pi(x, prec):
    """(r, quadrant): x = quadrant (pi / 2) + r with r in [-pi/4, pi/4] (approximately: the
    reduction is exact to the working precision plus the digits lost to |x|)."""
    with localcontext() as context:
        context.prec = prec + 20 + max(0, len(str(int(abs(x)))))
        half_pi = pi_decimal(context.prec) / 2
        value = Decimal(x)                            # exact
        quadrant = int((value / half_pi).to_integral_value())   # nearest
        r = value - quadrant * half_pi
    return r, quadrant % 4


def _sin_series(r, prec):
    """sin r by its Taylor series, |r| <= pi / 4 + epsilon."""
    with localcontext() as context:
        context.prec = prec + 4
        r = +r
        r2 = r * r
        term, total, k = r, r, 1
        threshold = Decimal(10) ** (-(prec + 6))
        while True:
            term = -term * r2 / ((2 * k) * (2 * k + 1))
            total += term
            k += 1
            if abs(term) < threshold * abs(total) or term == 0:
                return total


def _cos_series(r, prec):
    """cos r by its Taylor series, |r| <= pi / 4 + epsilon."""
    with localcontext() as context:
        context.prec = prec + 4
        r = +r
        r2 = r * r
        term, total, k = Decimal(1), Decimal(1), 1
        threshold = Decimal(10) ** (-(prec + 6))
        while True:
            term = -term * r2 / ((2 * k - 1) * (2 * k))
            total += term
            k += 1
            if abs(term) < threshold * abs(total) or term == 0:
                return total


def _sin_decimal(x, prec):
    r, quadrant = _reduce_half_pi(x, prec)
    if quadrant == 0:
        return _sin_series(r, prec)
    if quadrant == 1:
        return _cos_series(r, prec)
    if quadrant == 2:
        return -_sin_series(r, prec)
    return -_cos_series(r, prec)


def _cos_decimal(x, prec):
    r, quadrant = _reduce_half_pi(x, prec)
    if quadrant == 0:
        return _cos_series(r, prec)
    if quadrant == 1:
        return -_sin_series(r, prec)
    if quadrant == 2:
        return -_cos_series(r, prec)
    return _sin_series(r, prec)


def sin(x):
    """The correctly rounded sine of a finite double."""
    _check_finite("sin", x)
    x = float(x)
    if x == 0.0:
        return x                                      # sin(+-0) = +-0
    return _round_to_double(lambda prec: _sin_decimal(x, prec), "sin")


def cos(x):
    """The correctly rounded cosine of a finite double."""
    _check_finite("cos", x)
    x = float(x)
    if x == 0.0:
        return 1.0
    return _round_to_double(lambda prec: _cos_decimal(x, prec), "cos")


def tan(x):
    """The correctly rounded tangent of a finite double."""
    _check_finite("tan", x)
    x = float(x)
    if x == 0.0:
        return x

    def evaluate(prec):
        with localcontext() as context:
            context.prec = prec + 4
            return _sin_decimal(x, prec + 4) / _cos_decimal(x, prec + 4)
    return _round_to_double(evaluate, "tan")


def _atan_decimal(t, prec):
    """atan t for a Decimal (or Fraction-exact) t >= 0 by argument halving (atan t = 2 atan
    (t / (1 + sqrt(1 + t^2)))) until t < 1/8, then the alternating series."""
    with localcontext() as context:
        context.prec = prec + 6
        t = +t
        if t < 0:
            return -_atan_decimal(-t, prec)
        halvings = 0
        eighth = Decimal(1) / 8
        while t > eighth:
            t = t / (1 + (1 + t * t).sqrt())
            halvings += 1
        t2 = t * t
        term, total, k = t, t, 1
        threshold = Decimal(10) ** (-(prec + 8))
        while True:
            term = -term * t2
            contribution = term / (2 * k + 1)
            total += contribution
            k += 1
            if abs(contribution) < threshold * abs(total) or contribution == 0:
                break
        return total * (2 ** halvings)


def atan(x):
    """The correctly rounded arctangent of a finite double."""
    _check_finite("atan", x)
    x = float(x)
    if x == 0.0:
        return x
    return _round_to_double(lambda prec: _atan_decimal(Decimal(x), prec), "atan")


def atan2(y, x):
    """The correctly rounded atan2(y, x) of finite doubles, with the C99 / IEEE signed-zero
    conventions: atan2(+-0, x > 0) = +-0, atan2(+-0, x < 0) = +-pi, atan2(+-0, +0) = +-0,
    atan2(+-0, -0) = +-pi, atan2(y > 0, +-0) = pi / 2, atan2(y < 0, +-0) = -pi / 2."""
    _check_finite("atan2", y, x)
    y, x = float(y), float(x)
    if y == 0.0:
        if x > 0.0 or (x == 0.0 and not math.copysign(1.0, x) < 0.0):
            return y                                  # +-0
        return math.copysign(math.pi, y)              # math.pi is the correctly rounded pi
    if x == 0.0:
        return math.copysign(math.pi / 2, y)          # pi / 2: an exact halving of math.pi

    def evaluate(prec):
        with localcontext() as context:
            context.prec = prec + 8
            ratio = Decimal(y) / Decimal(x)           # both exact; one rounding at prec + 8
            if abs(ratio) <= 1:
                angle = _atan_decimal(abs(ratio), prec)
            else:
                angle = pi_decimal(prec + 8) / 2 - _atan_decimal(1 / abs(ratio), prec)
            if x < 0.0:
                angle = pi_decimal(prec + 8) - angle
            return angle if y > 0.0 else -angle
    return _round_to_double(evaluate, "atan2")


def acos(x):
    """The correctly rounded arccosine of a double in [-1, 1] (acos 1 = 0, acos -1 = pi, acos 0
    = pi / 2 exactly as rounded constants; elsewhere 2 atan(sqrt((1 - x) / (1 + x))) with the
    difference 1 - x formed exactly)."""
    _check_finite("acos", x)
    x = float(x)
    if not -1.0 <= x <= 1.0:
        raise ValueError(f"deterministic_math.acos: {x!r} is outside [-1, 1]")
    if x == 1.0:
        return 0.0
    if x == -1.0:
        return math.pi
    if x == 0.0:
        return math.pi / 2

    def evaluate(prec):
        with localcontext() as context:
            context.prec = prec + 8
            exact = Fraction(x)
            ratio = Fraction(1) - exact
            ratio /= Fraction(1) + exact              # exact rational (1 - x) / (1 + x)
            root = (Decimal(ratio.numerator) / Decimal(ratio.denominator)).sqrt()
            return 2 * _atan_decimal(root, prec)
    return _round_to_double(evaluate, "acos")


def asin(x):
    """The correctly rounded arcsine of a double in [-1, 1] (asin 0 = +-0, asin +-1 = +-pi / 2
    as rounded constants; elsewhere atan(x / sqrt(1 - x^2)) with 1 - x^2 formed exactly)."""
    _check_finite("asin", x)
    x = float(x)
    if not -1.0 <= x <= 1.0:
        raise ValueError(f"deterministic_math.asin: {x!r} is outside [-1, 1]")
    if x == 0.0:
        return x
    if abs(x) == 1.0:
        return math.copysign(math.pi / 2, x)

    def evaluate(prec):
        with localcontext() as context:
            context.prec = prec + 8
            exact = Fraction(x)
            complement = Fraction(1) - exact * exact  # exact rational 1 - x^2
            root = (Decimal(complement.numerator) / Decimal(complement.denominator)).sqrt()
            ratio = Decimal(abs(x)) / root
            angle = _atan_decimal(ratio, prec)
            return angle if x > 0.0 else -angle
    return _round_to_double(evaluate, "asin")


def golden_hex(function, *arguments):
    """``float.hex()`` of ``function(*arguments)``: the platform-independent key of the golden
    tables in test_deterministic_math.py."""
    return function(*arguments).hex()
