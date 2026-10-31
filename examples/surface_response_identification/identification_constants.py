#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The identification's matching-band constants shared BY NAME with the C++ (identification fix round
DESIGN sections 4 / 5 and A.5; decisions 1062 / 1075 / 1079 / 1097 Q2): the Python readers / discovery /
library tools use these values, and test_identification_constants.py asserts each equals its C++ definition
in palace/models/surfaceresponseidentification.hpp (the constant's name is the link; a change on either side
fails the test).

- GDS_SCATTER_OVER_RADIUS (kGdsScatterOverRadius): the translational pair families' near-match band
  max(SeparationTolerance, GDS_SCATTER_OVER_RADIUS x R) -> status NearMatched (Exact keeps precedence).
- GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS (kGeometricNearMatchClaimsOverRadius): the claims band of the
  geometric near-match (every ParallelEdgeCluster offset / SpatialEdgeCluster claims number) -> status
  GeometricNearMatched.
- GEOMETRIC_NEAR_MATCH_CONTEXT_OVER_RADIUS / _ANGLE_DEGREES / _GAP_COSINE: the context band, the angle band
  and the gap-direction cosine band of the same rule (SpatialEdgeCluster keys; the C++'s one-to-one context
  assignment is the rule, decision 1079 Q2).
- GEOMETRIC_NEAR_PAIR_CORRECTION_BAR: the per-Type correction-difference bar (2 %) of a library pair inside
  the band (DESIGN A.4 (c) / A.5 (b): a calibration bar set after the measurement; a pair above it is RECORDED
  and listed, never excluded).
"""

import os
import re

CXX_HEADER = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "palace", "models",
                          "surfaceresponseidentification.hpp")

GDS_SCATTER_OVER_RADIUS = 5.0e-3
GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS = 1.0e-2
GEOMETRIC_NEAR_MATCH_CONTEXT_OVER_RADIUS = 1.0e-1
GEOMETRIC_NEAR_MATCH_ANGLE_DEGREES = 0.1
GEOMETRIC_NEAR_MATCH_GAP_COSINE = 2.0e-2
GEOMETRIC_NEAR_PAIR_CORRECTION_BAR = 0.02

# Python name -> the C++ constexpr it mirrors.
CXX_NAMES = {
    "GDS_SCATTER_OVER_RADIUS": "kGdsScatterOverRadius",
    "GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS": "kGeometricNearMatchClaimsOverRadius",
    "GEOMETRIC_NEAR_MATCH_CONTEXT_OVER_RADIUS": "kGeometricNearMatchContextOverRadius",
    "GEOMETRIC_NEAR_MATCH_ANGLE_DEGREES": "kGeometricNearMatchAngleDegrees",
    "GEOMETRIC_NEAR_MATCH_GAP_COSINE": "kGeometricNearMatchGapCosine",
}

# The statuses the readers count beside Exact / Interpolated / Missing (DESIGN sections 4 / 5; decision
# 1075 Q1 / 1079): a key with one of these is CORRECTED (not Missing: no placeholder, no coupon).
NEAR_MATCHED_STATUSES = ("NearMatched", "GeometricNearMatched")
# The translational pair topologies the GDS scatter band applies to (the C++ FindLibraryModel).
TRANSLATIONAL_PAIR_TOPOLOGIES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip")


def cxx_constant(name, header=CXX_HEADER):
    """The value of `constexpr double <name> = <literal>;` in the C++ header (ValueError when absent)."""
    with open(header) as f:
        text = f.read()
    m = re.search(r"constexpr\s+double\s+" + re.escape(name) + r"\s*=\s*([-+0-9.eE]+)\s*;", text)
    if m is None:
        raise ValueError(f"{name} is not a numeric constexpr double of {header}")
    return float(m.group(1))


def pair_near_match_band(separation_tolerance, radius):
    """The band (mesh / library length units) a translational pair key is near-matched within."""
    return max(float(separation_tolerance), GDS_SCATTER_OVER_RADIUS * float(radius))
