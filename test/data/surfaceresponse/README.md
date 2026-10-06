<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Surface-response golden files

- `corner-arm-trim-90deg-trim-disabled-patches.csv` — the surface-response patch dry run
  (`surface-response-patches.csv`) of the 90-degree two-corner PEC lead of
  `test/unit/test-cornerarmtrim.cpp` ("SurfaceResponseOperator corner-arm trim", section
  "a 90-degree corner keeps the legacy layout"), written by a build in which the
  `ApplyCornerArmTrim` call at the end of `BuildFeaturePatches` was removed (the corner-arm
  trim of decision 394 F1 disabled; branch `simlapointe/corner-overlap-uncovered-energy`,
  2026-10-06). Serial and 2-rank runs wrote identical bytes (sha256
  aa7703a0215b21dc2490484cd04c62aa21430fe7b56ea3e04a5890081b75cf25). The test compares the
  trimming build's dry run against it byte for byte (decision 399 MINOR-6): at 90 degrees
  the trim must leave the legacy layout untouched.
