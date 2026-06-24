# Autoresearch: Palace ParaView streaming for large SingleTransmon AMR

## Objective
Reduce the **ParaView** timer for the large SingleTransmon AMR GPU workload on p4d, starting from branch:

```text
origin/hughcars/libceed-output-functionals-dev-paraview-stream
```

This branch already replaced the worst `-auto2` behavior of materializing derived domain ParaView fields as intermediate `ParGridFunction`s. It adds `CeedParaViewDataCollection` and a `DomainFieldEvaluator` VTU point-buffer mode. Your job is to keep pushing that direction and reduce the large AMR SingleTransmon ParaView time further.

## Workload
The measurement script runs:

- `examples/transmon/transmon_amr.json`
- Eigenmode `N=2`, `Save=2`
- `Solver.Order=3`
- `Model.Refinement.MaxIts=2`
- `Solver.Device="GPU"`
- `OutputFormats.Paraview=true`
- `OutputFormats.GridFunction=true`
- current p4d autoresearch objective uses 2 GPUs: `PALACE_AR_GPU_NP=2` (keep fixed during the search; run all-GPU validation separately before finalizing)

This is the same large nonconforming AMR case that exposed the `-auto2` ParaView regression:

- `origin/main` reference observed earlier: ParaView ~220 s, GridFunction ~188 s, total ~983 s
- `-auto2` before this branch: ParaView ~524 s, GridFunction ~6.8 s, total ~1117 s
- expected target for this branch: keep GridFunction fast while driving ParaView toward/below the old coefficient-streaming path

## Primary metric
- **transmon_amr_paraview_seconds** (seconds, lower is better): Avg. value from Palace's `Elapsed Time Report`, row `Paraview`.

Keep a result only when this metric improves and scalar-output checks pass. Re-run close wins if the improvement is within noise.

## Secondary metrics
The measure script also emits:

- `transmon_amr_total_seconds`
- `transmon_amr_postprocessing_seconds`
- `transmon_amr_gridfunction_seconds`
- `transmon_amr_operator_construction_seconds`
- `scalar_max_rel`, `scalar_max_abs`, `scalar_compared_files`

Do not trade a large GridFunction or total-runtime regression for a tiny ParaView win. The main win must be in output/postprocessing, not the solver.


## Build loop
`.auto/measure.sh` keeps the workflow Spack-mediated but avoids a full Spack reinstall when possible:

1. Find the kept Spack Palace CMake/Ninja stage for this worktree.
2. Run `spack -e "$ROOT" build-env palace -- ninja -C "$build_dir" -j "$PALACE_AR_BUILD_JOBS" install`.
3. If no kept stage exists, or incremental Ninja fails, fall back to `spack install --only package --overwrite --keep-stage -j "$PALACE_AR_BUILD_JOBS"` so the next iteration has a reusable stage.

Do not bypass Spack with an ad-hoc CMake build; use the measured harness path. Set `PALACE_AR_INCREMENTAL_BUILD=0` only to force the conservative fallback.

## Correctness guard
`.auto/measure.sh` creates a scalar CSV baseline on the first unmodified run and compares later runs against it. `.auto/checks.sh` fails if scalar drift exceeds the threshold recorded by the measure script.

The scalar guard is intentionally local to this transmon workload. Before finalizing any kept changes, run the broader p4d GPU+CPU validation suite from the project playbook.

## Files in scope
Primary files:

- `palace/utils/ceedparaviewdatacollection.hpp`
- `palace/utils/ceedparaviewdatacollection.cpp`
- `palace/fem/surfacefunctional.hpp` (the embedded `DomainFieldEvaluator`)
- `palace/fem/surfacefunctional.cpp` (the embedded `DomainFieldEvaluator`)
- `palace/models/postoperator.hpp`
- `palace/models/postoperator.cpp`

Possible supporting files:

- `palace/utils/CMakeLists.txt`
- small test files under `test/unit/` if you add targeted coverage

## Off limits
- Do **not** optimize or tune eigensolver/linear solver behavior.
- Do **not** drop required output fields or disable ParaView/GridFunction output to win the metric.
- Do **not** remove scalar correctness checks.
- Do **not** run direct CMake/make builds. Use Spack only.
- Do **not** use more than `PALACE_AR_BUILD_JOBS` build jobs; default here is 32 on p4d.
- Do **not** change transient, BoundaryMode, or 2D fallback behavior unless required to keep compilation working.

## Starting implementation notes
Current branch state:

- `CeedParaViewDataCollection` subclasses MFEM `ParaViewDataCollection` and adds domain/boundary point-field maps.
- `DomainFieldEvaluator` can fill either an interpolatory L2 `ParGridFunction` or a full VTU point buffer in MFEM's element/refined-point traversal order.
- `PostOperator` registers derived domain fields (`U_e`, `U_m`, `S`) as direct domain point fields for ParaView and still uses GridFunction evaluators for MFEM grid-function output.
- Boundary libCEED visualization buffers are registered as direct boundary point fields.

Promising directions:

1. Avoid full-size temporary point buffers where possible. Current streaming still allocates complete point buffers before writing; a chunked writer could evaluate one geometry/rank chunk, copy/write it, then reuse scratch.
2. Reduce duplicate host/device transfers and repeated `HostRead()` synchronization during VTU writing.
3. Avoid repeated refined-geometry and header work when MFEM already traverses the same element/order lattice.
4. Split the timer inside ParaView save if needed, but remove instrumentation before keeping unless it remains useful and low-noise.
5. Preserve the direct VTU point-field API shape: it should be easy to reason about and should not depend on fuzzy floating-point point lookup.

## How to run

```bash
./.auto/measure.sh
```

This builds with Spack and runs the large transmon AMR case. It prints `METRIC` lines for autoresearch.

## Loop rules
1. First action: run `.auto/measure.sh` unmodified and log the baseline before editing code.
2. Inspect source and logs deeply before each change; this workload is expensive.
3. Keep atomic commits with clear metric deltas.
4. Discard regressions or correctness failures.
5. Update this file or `.auto/ideas.md` when you learn something that should survive context compaction.
Operational storage policy: preserve scalar CSV baselines and an optional origin/main ParaView/GridFunction reference tree at PALACE_AR_REFERENCE_POST_ROOT. Candidate runs must write large output to PALACE_AR_POST_ROOT on /home/ubuntu/workspace and overwrite it each iteration, keeping only the most recent candidate ParaView output. Do not preserve per-commit/per-run ParaView directories. The origin/main reference output files are expected to remain valid because output filenames and values should not change across writer-performance experiments.
