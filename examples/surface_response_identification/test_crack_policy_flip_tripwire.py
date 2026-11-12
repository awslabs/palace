#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tripwire of ruling 1586 (upstream integration P2; SL9 review m-2): CRACK_POLICY_DEFAULT may leave "legacy"
only once the crack policy reaches the library-build orchestrators AND their cache keys (the Tool /
Generator fingerprints of the coupon specs). Until then a cache built on cracked coupons would be reused
unchanged for an uncracked build, which must never happen. The plumbing is a P3 item and a precondition of
the flip ruling (decision 1570 (2))."""

import os
from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import crack_policy as C  # noqa: E402

EXAMPLES = Path(__file__).resolve().parents[1]
# The library-build orchestrators the policy must reach (SL9 review m-2).
ORCHESTRATORS = (
    "cpw2d/prepare_surface_response_coupons.py",
    "cpw2d/build_surface_process_library.py",
    "cpw3d_surface/spatial_coupon/device_coupons.py",
    "cpw3d_surface/spatial_coupon/qualify/qualify_library.py",
    "cpw3d_surface/spatial_coupon/qualify/case_inputs.py",
)
FINGERPRINT_CALL = re.compile(r"(?<!def )\b(tool_fingerprint|generator_fingerprint)\(")


def call_arguments(text, start):
    """The text between the parenthesis opened at text[start - 1] and its match."""
    depth, i = 1, start
    while depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return text[start:i - 1]


def missing_plumbing():
    """What the flip still lacks: an orchestrator that passes no crack policy, or a fingerprint (cache key)
    that does not include crack_policy.py."""
    out = []
    for rel in ORCHESTRATORS:
        text = (EXAMPLES / rel).read_text()
        if "crack_policy=" not in text:
            out.append(f"{rel}: passes no crack policy (crack_policy=...)")
        for match in FINGERPRINT_CALL.finditer(text):
            args = call_arguments(text, match.end())
            if "crack_policy" not in args.lower():
                line = text.count("\n", 0, match.start()) + 1
                out.append(f"{rel}:{line}: {match.group(1)}(...) does not fingerprint crack_policy.py")
    return out


def flip_problems(default):
    """Empty when the default may stand: "legacy", or any default once the plumbing exists."""
    return [] if default == "legacy" else missing_plumbing()


class CrackPolicyFlipTripwireTest(unittest.TestCase):
    def test_default_leaves_legacy_only_with_the_plumbing(self):
        problems = flip_problems(C.CRACK_POLICY_DEFAULT)
        self.assertEqual(problems, [], "ruling 1586: CRACK_POLICY_DEFAULT = "
                         f"{C.CRACK_POLICY_DEFAULT!r} needs the crack policy in every library-build orchestrator "
                         "and in every Tool / Generator fingerprint first (a cache built cracked must never be "
                         "reused for an uncracked build): " + "; ".join(problems))

    def test_tripwire_reads_every_orchestrator(self):
        self.assertEqual(flip_problems("legacy"), [])
        self.assertEqual(flip_problems("default"), missing_plumbing())
        calls = sum(len(FINGERPRINT_CALL.findall((EXAMPLES / rel).read_text())) for rel in ORCHESTRATORS)
        self.assertGreaterEqual(calls, 6)  # 4 ToolFingerprint + 2 GeneratorFingerprint specs at P2


if __name__ == "__main__":
    unittest.main()
