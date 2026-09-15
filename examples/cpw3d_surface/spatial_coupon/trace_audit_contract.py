#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Geometry-independent seam for future coupon trace audits.

The existing 135-column corrected-model diagnostic is intentionally not wired
through this seam and remains outside the geometry-independence release gate.
New three/six/ten-edge producers must supply this semantic contract instead of
copying the historical 40/95 partition.
"""


def trace_groups(contract):
    if not isinstance(contract, dict) or contract.get("Version") != 1:
        raise ValueError("Unsupported trace audit contract")
    sources = contract.get("SourceIndices")
    constrained = contract.get("ConstrainedIndices")
    classes = contract.get("SemanticExcitationClasses")
    if (not isinstance(sources, list) or not sources or
            any(not isinstance(value, int) or value <= 0 for value in sources) or
            len(set(sources)) != len(sources)):
        raise ValueError("SourceIndices must be nonempty, positive, and unique")
    if (not isinstance(constrained, list) or len(set(constrained)) != len(constrained) or
            not set(constrained) <= set(sources)):
        raise ValueError("ConstrainedIndices must be a unique source subset")
    if (not isinstance(classes, list) or not classes or
            any(not isinstance(item, dict) or not item.get("Name") or
                not isinstance(item.get("Indices"), list) or not item["Indices"] or
                not set(item["Indices"]) <= set(sources) for item in classes)):
        raise ValueError("Semantic excitation classes must be nonempty source subsets")
    return {"all": sources,
            "free": [value for value in sources if value not in set(constrained)],
            "constrained": constrained,
            "semantic": {item["Name"]: item["Indices"] for item in classes}}
