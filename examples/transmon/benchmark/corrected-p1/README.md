<!--
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
-->

# Fixed-mesh corrected single-transmon benchmark (p=1)

This benchmark exercises the complete Maxwell surface-response path on the checked-in
coarse single-transmon mesh:

- automatic geometry matching;
- distributed response patches and Maxwell contour traces;
- fixed-trace and fixed-flux postprocessing;
- self-consistent corrected eigenmodes;
- corrected divergence-free projection;
- mode pairing and confidence output.

The device and coupon FEM order is `p=1`. AMR, uniform refinement, `EdgeRefinement`, field
export, and localized per-segment output are disabled. The run requests one mode, although
the eigensolver may converge and report an additional mode from its invariant subspace.

## Benchmark-only process library

`library/` contains an unqualified, intentionally coarse software benchmark library with:

- 12 requested isolated-edge knots, with one PEC-constrained knot removed to leave an
  effective 11-function trace basis;
- a 12-function same-conductor-strip basis at 2 um separation;
- p=1 coupon solves on quadratic coupon meshes.

The exact generation settings and command are recorded in `generation.json`. Each model
is explicitly marked `BoundaryLawQualification.Status: Unqualified`; this library is
**not** a physically qualified process library and must not be used for design results.
Unsupported nonparallel, spatial, and corner neighborhoods use the documented `Warn`
behavior; their exact coverage fractions are part of the numerical baseline.

## Running

From the repository root:

```text
python3 examples/transmon/benchmark/corrected-p1/run_corrected_eigenmode.py \
  --output /tmp/transmon-corrected-p1 \
  --palace build/bin/palace \
  --launcher "$(command -v mpiexec)" \
  --ranks 1 2 4
```

The default is one run at each rank count. On the development Apple arm64 host, observed
wall times were approximately 153, 84, and 60 seconds for one, two, and four ranks. Peak
node memory was approximately 3.7--4.0 GiB. These values are observational, not portable
pass/fail limits.

The harness verifies:

- exact fixture hashes;
- p=1 and fixed-mesh safety invariants;
- the requested MPI size;
- exact patch/trace/line/point/stencil workload counts;
- raw, fixed-trace, fixed-flux, and self-consistent CSV values against rank-specific
  numerical baselines.

To create a reviewed candidate after an intentional numerical or fixture change:

```text
python3 examples/transmon/benchmark/corrected-p1/run_corrected_eigenmode.py \
  --output /tmp/transmon-corrected-p1-candidate-run \
  --ranks 1 2 4 \
  --write-baseline-candidate /tmp/transmon-corrected-p1-candidate.json
```

Review the numerical and workload diff before replacing
`corrected-eigenmode-baseline.json`.
