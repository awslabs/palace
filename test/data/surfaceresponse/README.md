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
    the trim must leave the legacy layout untouched. The `DomainWeight` column (decision 559
    MINOR-3 / rerun-3 plan O-14, the mirror-formed contract lane `simlapointe/mirror-formed-
    contract`, 2026-10-08) was appended to every row of the golden (`,DomainWeight` in the
    header, `,1` per patch: the dry run's domain weights are all 1 here) when the patch dry
    run gained it; no other byte changed.

  - `corner-metal-graded-rule/{concave,convex}-{90,60,48.75}deg-basis-points.csv` + `manifest.json` —
    the Python generator's basis points of the METAL-GRADED corner trace basis rule
    (`generate_corner_response.METAL_GRADED_RULE`, `all-rings-follow-metal-metal-graded`: the refined
    rule with MetalKnotGrading (1/9, 1/3), MetalInteriorKnots 9; 11 rings x 20 knots = 220; R = 1.9,
    MetalThickness 0.1, OveretchDepth 0.05; floats written with `repr`), at 90 / 60 / 48.75 degrees for
    both convexities, with each file's zero set (22), ContourGroups and sha256 in `manifest.json`
    (the concave-MS decide step, branch `simlapointe/concave-ms-decide`, 2026-10-09; decisions 590 /
    629 O-9; the (B) 2.4 (iii) test form of decision 565 MAJOR-2). `test-cornerbasisrefinement.cpp`
    `CornerMetalGradedRulePythonParity` reads them and requires the C++ rule's knots
    (`MakeCornerBoxSeed` + `BuildCornerTraceBasis`) within 1e-9 R, the same zero set and the same
    ring counts — the Python / C++ layout parity pin on real generator output (the 48.75-degree
    concave node exercises the scaled FREE grading, the 48.75-degree convex node the scaled METAL
    grading: its metal arc 0.93 R is below the reference 1.4226 R).
