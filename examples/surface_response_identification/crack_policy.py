#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Model.CrackInternalBoundaryElements policy of the fresh Palace config writers (decision 1570).

Upstream Palace PR 994 changed what an UNSET key means: unset = upstream's default (PEC / Ground /
Terminal sheets are not cracked; PMC / ZeroCharge and absorbing (farfield) boundaries are always cracked;
Conductivity, Impedance and RationalImpedance boundaries are cracked only where marked "Crack");
true = the legacy behaviour (every interior boundary with a boundary condition except
lumped ports is cracked; deprecated); false = nothing is cracked. Upstream's default has no explicit
spelling, so the two policies are:

  legacy   writes Model.CrackInternalBoundaryElements: true (the behaviour of an unset key on the
           pre-integration binary, whose ModelData::crack_bdr_elements defaults to true);
  default  omits the key (upstream's default; per-boundary "Crack" entries are the caller's business).

Every fresh writer takes --crack-policy with CRACK_POLICY_DEFAULT as its default value. The library-build
orchestrators (the coupon planner, the straight-family builder, the device adapter, the spatial
qualification) take it too, pass it to the writers they run, and key their caches by it (decision 1586).
"""

from pathlib import Path

CRACK_KEY = "CrackInternalBoundaryElements"
CRACK_POLICIES = ("legacy", "default")
# This module's own path: the orchestrators' Tool / Generator fingerprints include its bytes, so a change of
# the policy rules (or of CRACK_POLICY_DEFAULT) changes their cache keys.
CRACK_POLICY_SOURCE = Path(__file__).resolve()

# The library models and every record of the stage-F freeze were made on cracked meshes, so the
# writers default to "legacy" (decision 1570; USER decision 1566 (5): no library rebuild now).
# Upstream's uncracked default is the TARGET (USER decision 1566 (3)): flipping this one value to
# "default" is a later ruling, after the uncracked test tiers pass.
CRACK_POLICY_DEFAULT = "legacy"


def add_crack_policy_argument(parser):
    """Add --crack-policy {legacy,default} to an argparse parser."""
    parser.add_argument(
        "--crack-policy", choices=CRACK_POLICIES, default=CRACK_POLICY_DEFAULT,
        help=f"Model.{CRACK_KEY}: 'legacy' writes true, 'default' omits the key "
             f"(upstream's uncracked default) (default: {CRACK_POLICY_DEFAULT})")


def check_crack_policy(policy):
    """Return `policy` if it is one of CRACK_POLICIES; refuse anything else."""
    if policy not in CRACK_POLICIES:
        raise ValueError(f"unknown crack policy {policy!r} (expected one of {CRACK_POLICIES})")
    return policy


def crack_policy_of(args):
    """The crack policy of parsed orchestrator options: args.crack_policy, else CRACK_POLICY_DEFAULT (the
    orchestrators' spec functions are also called with partial option namespaces)."""
    return check_crack_policy(getattr(args, "crack_policy", CRACK_POLICY_DEFAULT))


def crack_policy_arguments(crack_policy):
    """The command-line arguments that pass `crack_policy` to a writer or orchestrator run as a subprocess."""
    return ["--crack-policy", check_crack_policy(crack_policy)]


def apply_crack_policy(model, policy):
    """Set Model.CrackInternalBoundaryElements of a fresh config's Model dict per the policy; returns
    the dict. A fresh writer's Model never carries the key already: that is refused."""
    check_crack_policy(policy)
    if CRACK_KEY in model:
        raise ValueError(f"Model.{CRACK_KEY} is already set; a crack policy applies to a fresh Model only")
    if policy == "legacy":
        model[CRACK_KEY] = True
    return model
