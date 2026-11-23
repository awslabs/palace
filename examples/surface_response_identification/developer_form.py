#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The two forms of the response-correction block a config writer of the examples writes (the v1 cut, decisions
2377 / 2379; v1/PARKED.md).

Palace corrects by post-processing only by default; a parked option of a response-correction block is refused by
name unless the same block sets "Developer": true (configfile.cpp RequireDeveloperOption). A writer builds its
RECORDED block (the settings of its records) and writes one of two forms:

  v1 form                  (the default) the recorded block less the parked options it selects: each option's
                           selecting setting set to its v1 value and the option's own keys dropped (PARKED_OPTIONS);
                           a device run with a library correction takes the crack policy V1_CRACK_POLICY.
  recorded developer form  (--developer) the recorded block with the defaults of record it relied on stated
                           (RECORD_DEFAULTS: Palace's defaults before the cut) and "Developer": true appended, with the
                           recorded crack policy RECORDED_CRACK_POLICY (the records ran on cracked meshes).

Un-parking one option later is ONE line here: delete its PARKED_OPTIONS entry (with its RequireDeveloperOption line
in Palace); the v1 form of every writer then keeps that option as recorded. Values from the v1 form are not
comparable with records made before this cut until the benchmark's collocated twins are read.
"""

import copy

DEVELOPER_KEY = "Developer"
V1_FORM = "v1 form"
DEVELOPER_FORM = "recorded developer form"
# The crack policies (crack_policy.py) of a device run with a library correction: the v1 form's (Palace refuses the
# global crack key over an interior perfect conductor without the switch) and the records'.
V1_CRACK_POLICY = "default"
RECORDED_CRACK_POLICY = "legacy"
# The defaults a recorded block relied on (Palace's before the cut), stated by the recorded developer form.
RECORD_DEFAULTS = {"CorrectionMode": "Both"}

# The parked options a recorded block may select: name -> (the selecting setting, the values that select it, the
# option's own keys (Palace reads them with the option only), the selecting setting's v1 value).
PARKED_OPTIONS = {
    "self-consistent solve": ("CorrectionMode", ("SelfConsistent", "Both"),
                              ("SolveTol", "TranslationalDomainCorrection"), "PostprocessOnly"),
    "mortar trace coupling": ("TraceCoupling", ("SurfaceMortar",), ("MortarOversampling",), "Collocated"),
}
# The other parked settings of a block, which no writer selects: (key path, the parked value; None = any value).
PARKED_SETTINGS = (
    (("PatchConstruction",), "Legacy"),
    (("ConductorConsistency", "Exclusion"), "Sticky"),
    (("AllowElectricOnlyFallback",), True),
    (("Magnetic",), True),
    (("ConductorConsistency", "ControlStateAnchorSetbackOverR"), None),
)
CORRECTION_PATHS = (("Solver", "Electrostatic", "ResponseCorrection"), ("Solver", "SurfaceResponseCorrection"))


def add_developer_argument(parser):
    parser.add_argument("--developer", action="store_true",
                        help=f"write the {DEVELOPER_FORM} (the recorded parked settings with \"{DEVELOPER_KEY}\": "
                             f"true) instead of the {V1_FORM}")


def form_name(developer):
    return DEVELOPER_FORM if developer else V1_FORM


def _node(tree, path):
    for key in path:
        tree = tree.get(key) if isinstance(tree, dict) else None
    return tree


def parked_block_settings(block):
    """The parked settings a response-correction block states, as Palace names them; an option's own keys count
    where the block does not select the option (Palace does not read them there)."""
    found = []
    for setting, values, keys, _ in PARKED_OPTIONS.values():
        if block.get(setting) in values:
            found.append(f"{setting} {block[setting]}")
        found.extend(key for key in keys if key in block)
    for path, value in PARKED_SETTINGS:
        parent = _node(block, path[:-1])
        if isinstance(parent, dict) and path[-1] in parent and (value is None or parent[path[-1]] == value):
            found.append(".".join(path) + ("" if value is None else f" {value}"))
    if "Library" not in block:
        found.append("the hand-placed form (no Library)")
    return found


def parked_config_settings(config):
    """The parked settings of a whole config: those of its correction blocks and, in a correction run, the global
    crack key."""
    found = []
    for path in CORRECTION_PATHS:
        block = _node(config, path)
        if isinstance(block, dict):
            found.extend(f"{'.'.join(path)}: {setting}" for setting in parked_block_settings(block))
            if block.get(DEVELOPER_KEY) is not None:
                found.append(f"{'.'.join(path)}: {DEVELOPER_KEY}")
    if any(isinstance(_node(config, path), dict) for path in CORRECTION_PATHS):
        if config.get("Model", {}).get("CrackInternalBoundaryElements") is True:
            found.append("Model.CrackInternalBoundaryElements true")
    return found


def device_run_crack_policy(crack_policy, developer):
    """The crack policy of a device run with a library correction: `crack_policy` when given, else the form's
    (V1_CRACK_POLICY; RECORDED_CRACK_POLICY with `developer`); the legacy policy, which cracks the interior perfect
    conductors, is a parked option there and needs `developer`."""
    policy = crack_policy or (RECORDED_CRACK_POLICY if developer else V1_CRACK_POLICY)
    if policy == "legacy" and not developer:
        raise ValueError("the crack policy legacy cracks the interior perfect conductors, a parked developer option of a "
                         "library correction run: write the recorded developer form (--developer)")
    return policy


def v1_block(recorded):
    """The v1 form of a recorded block; refuses a block that keeps a parked setting no PARKED_OPTIONS entry removes."""
    block = copy.deepcopy(recorded)
    for setting, values, keys, v1_value in PARKED_OPTIONS.values():
        if block.get(setting) in values:
            block[setting] = v1_value
        for key in keys:
            block.pop(key, None)
    left = parked_block_settings(block) + ([DEVELOPER_KEY] if DEVELOPER_KEY in block else [])
    if left:
        raise ValueError(f"the recorded block keeps the parked settings {left}: it has no {V1_FORM}")
    return block


def with_switch(block):
    """`block` with "Developer": true appended (a copy)."""
    block = copy.deepcopy(block)
    block.pop(DEVELOPER_KEY, None)
    block[DEVELOPER_KEY] = True
    return block


def developer_block(recorded):
    """The recorded developer form of a recorded block: its defaults of record stated, then the switch appended."""
    block = copy.deepcopy(recorded)
    for key, value in RECORD_DEFAULTS.items():
        block.setdefault(key, value)
    return with_switch(block)
