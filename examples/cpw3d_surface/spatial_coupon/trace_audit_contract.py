#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Semantic trace-contract seam for future three/six/ten-edge physics audits.

This validates excitation intent and source-bank membership only. It makes no
field, response, accuracy, or physics qualification claim.
"""
import math


REQUIRED_SINGLETON_KINDS = {"smooth-global", "localized-matching",
                            "deterministic-combination"}
ALLOWED_KINDS = REQUIRED_SINGLETON_KINDS | {"conductor-interface"}


def _conductor_interface(item):
    conductor = item.get("Conductor")
    interface = item.get("Interface")
    if (not isinstance(conductor, int) or conductor <= 0 or
            not isinstance(interface, str) or not interface.strip()):
        raise ValueError("Conductor/interface semantics must be explicit")
    return conductor, interface


def trace_groups(contract):
    if not isinstance(contract, dict) or contract.get("Version") != 1:
        raise ValueError("Unsupported trace audit contract")
    sources = contract.get("SourceIndices")
    constrained = contract.get("ConstrainedIndices")
    classes = contract.get("SemanticExcitationClasses")
    required_interfaces = contract.get("RequiredConductorInterfaces")
    if (not isinstance(sources, list) or not sources or
            any(not isinstance(value, int) or value <= 0 for value in sources) or
            len(set(sources)) != len(sources)):
        raise ValueError("SourceIndices must be nonempty, positive, and unique")
    if (not isinstance(constrained, list) or len(set(constrained)) != len(constrained) or
            not set(constrained) <= set(sources)):
        raise ValueError("ConstrainedIndices must be a unique source subset")
    if not isinstance(required_interfaces, list) or not required_interfaces:
        raise ValueError("RequiredConductorInterfaces must be nonempty")
    required_pairs = [_conductor_interface(item) if isinstance(item, dict) else None
                      for item in required_interfaces]
    if None in required_pairs or len(set(required_pairs)) != len(required_pairs):
        raise ValueError("Required conductor/interface classes must be unique")
    if not isinstance(classes, list) or not classes:
        raise ValueError("SemanticExcitationClasses must be nonempty")

    names = []
    kinds = []
    semantic = {}
    actual_pairs = []
    for item in classes:
        if (not isinstance(item, dict) or not isinstance(item.get("Name"), str) or
                not item["Name"].strip() or item.get("Kind") not in ALLOWED_KINDS or
                not isinstance(item.get("Indices"), list) or not item["Indices"] or
                any(not isinstance(value, int) for value in item["Indices"]) or
                len(set(item["Indices"])) != len(item["Indices"]) or
                not set(item["Indices"]) <= set(sources)):
            raise ValueError("Every semantic excitation needs a unique name, explicit kind, "
                             "and nonempty unique source subset")
        names.append(item["Name"])
        kinds.append(item["Kind"])
        if item["Kind"] == "conductor-interface":
            actual_pairs.append(_conductor_interface(item))
        elif "Conductor" in item or "Interface" in item:
            raise ValueError("Only conductor-interface excitations name conductor semantics")
        if item["Kind"] == "deterministic-combination":
            coefficients = item.get("Coefficients")
            if (not isinstance(coefficients, list) or len(coefficients) != len(item["Indices"]) or
                    any(not isinstance(value, (int, float)) or isinstance(value, bool) or
                        not math.isfinite(value) or value == 0 for value in coefficients)):
                raise ValueError("Deterministic combinations need finite nonzero coefficients")
        elif "Coefficients" in item:
            raise ValueError("Coefficients are reserved for deterministic combinations")
        semantic[item["Name"]] = {"kind": item["Kind"], "indices": item["Indices"]}
    if len(set(names)) != len(names):
        raise ValueError("Semantic excitation names must be unique")
    if not REQUIRED_SINGLETON_KINDS <= set(kinds):
        raise ValueError("Smooth/global, localized matching, and deterministic roles are required")
    if len(actual_pairs) != len(set(actual_pairs)) or set(actual_pairs) != set(required_pairs):
        raise ValueError("Every required conductor/interface class needs exactly one excitation")

    constrained_set = set(constrained)
    return {"all": sources,
            "free": [value for value in sources if value not in constrained_set],
            "constrained": constrained,
            "semantic": semantic,
            "conductor_interfaces": required_pairs,
            "PhysicsQualified": False}
