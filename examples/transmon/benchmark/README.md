<!--
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
-->

# Surface-response preflight benchmark

This benchmark uses the checked-in coarse single-transmon mesh to measure automatic
surface-response geometry classification without assembling or solving a field problem.
It disables AMR, uniform refinement, `EdgeRefinement`, field output, and localized
per-segment energy output. The ordinary deterministic crack/interface preprocessing used
by the practical solve mesh remains enabled and is part of the baseline.

The benchmark checks that the canonical preflight manifest is invariant across repetitions
and MPI rank counts. Timing and memory measurements are observational and are not used as
portable pass/fail thresholds.

From the repository root:

```text
python3 examples/transmon/benchmark/run_surface_response_preflight.py \
  --output /tmp/transmon-surface-preflight-benchmark \
  --palace build/bin/palace \
  --launcher "$(command -v mpiexec)" \
  --ranks 1 2 4 \
  --repeats 3
```

The command writes one configuration, log, raw manifest, and canonical manifest per run,
plus `benchmark-results.json`. Absolute process-library paths and the optional
observability-only `Statistics` object are excluded from the correctness hash.

The checked-in baseline was recorded at Palace commit `da467f864` on an Apple arm64 CPU
with `OMP_NUM_THREADS=1`. Median external wall times were 4.59, 3.87, and 3.52 seconds for
one, two, and four MPI ranks, respectively. These timings describe one machine only.

To create a candidate after an intentional fixture or manifest change:

```text
python3 examples/transmon/benchmark/run_surface_response_preflight.py \
  --output /tmp/transmon-surface-preflight-candidate-run \
  --ranks 1 2 4 \
  --repeats 1 \
  --write-baseline-candidate /tmp/transmon-surface-preflight-candidate.json
```

Review the complete canonical-manifest diff and fixture hashes, then explicitly copy the
candidate to `surface-response-preflight-baseline.json`. Candidate generation requires
exactly the 1/2/4-rank set and refuses to overwrite the checked-in baseline directly.

Do not add `EdgeRefinement`, uniform refinement, or AMR to this benchmark. Refinement is
validated separately on small meshes with explicit resource bounds.
