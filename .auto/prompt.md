# Autoresearch: Palace ParaView writer proxy optimization

## Objective
Reduce the **ParaView writer time** reported by the hidden nonconforming writer proxy benchmark on p4d.

The benchmark is solve-free: it builds deterministic high-order scalar/vector fields on a refined transmon mesh and writes both ParaView VTU and MFEM GridFunction outputs. Use it to isolate writer/materialization overhead without spending iterations on eigensolver runtime.

## Workload
The measurement script runs the hidden Catch2 benchmark:

```text
palace-unit-tests "[paraview-writer-proxy]" --skip-benchmarks
```

through the Spack environment, under MPI with:

```text
PALACE_PROXY_MPI_NP=2
PALACE_PROXY_ORDER=3
PALACE_PROXY_REFINE_STEPS=2
PALACE_PROXY_REFINE_PROB=0.3
PALACE_PROXY_SEED=20260624
PALACE_PROXY_COMPRESSION=1
```

Keep these defaults fixed during the search unless a human explicitly changes the proxy target. The deterministic seed and refinement controls are part of the benchmark definition.

## Primary metric
- **proxy_paraview_seconds** (seconds, lower is better): wall time for `mfem::ParaViewDataCollection::Save()` in the proxy test.

Keep a result only when this metric improves and `.auto/checks.sh` passes. Re-run close wins if the improvement is plausibly noise.

## Secondary metrics
The measure script also emits:

- `proxy_gridfunction_seconds`
- `proxy_paraview_bytes`
- `proxy_gridfunction_bytes`
- `proxy_mpi_ranks`
- `proxy_elements`
- `proxy_scalar_true_dofs`
- `proxy_vector_true_dofs`
- `proxy_ok`

Do not trade a large GridFunction regression for a tiny ParaView win. The useful win is a writer/materialization improvement, not a change in field count, problem size, solver behavior, or MPI rank count.

## Build loop
`.auto/measure.sh` keeps the workflow Spack-mediated:

1. Ensure the environment is concretized with Palace test dependencies (`catch2`).
2. Run `spack -e "$ROOT" install --test=root ... --keep-stage` so `palace-unit-tests` and test data are installed.
3. Run the hidden proxy benchmark through `spack build-env palace`.

Do not bypass Spack with an ad-hoc CMake or make build. Use the measured harness path.

## Correctness guard
The proxy benchmark is a hidden Catch2 test with deterministic mesh refinement/data and internal assertions. `.auto/checks.sh` requires the last proxy run to have completed successfully.

Before finalizing any kept code changes, run the broader p4d GPU+CPU validation suite from the project playbook and representative full-output cases separately.

## Files in scope
Primary files:

- `palace/utils/ceedparaviewdatacollection.hpp`
- `palace/utils/ceedparaviewdatacollection.cpp`
- `palace/fem/surfacefunctional.hpp` (the embedded `DomainFieldEvaluator`)
- `palace/fem/surfacefunctional.cpp` (the embedded `DomainFieldEvaluator`)
- `palace/models/postoperator.hpp`
- `palace/models/postoperator.cpp`
- `test/unit/test-paraview-writer-proxy.cpp`

Possible supporting files:

- `palace/utils/CMakeLists.txt`
- `test/unit/CMakeLists.txt`
- small targeted tests under `test/unit/`

## Off limits
- Do **not** optimize or tune eigensolver/linear solver behavior.
- Do **not** drop required output fields or disable ParaView/GridFunction output to win the metric.
- Do **not** change the proxy MPI rank count or deterministic proxy size controls unless a human changes the benchmark target.
- Do **not** run direct CMake/make builds. Use Spack only.
- Do **not** use more than `PALACE_AR_BUILD_JOBS` build jobs; default here is 32 on p4d.
- Do **not** change transient, BoundaryMode, or 2D fallback behavior unless required to keep compilation working.

## Starting implementation notes
Current branch state:

- `CeedParaViewDataCollection` subclasses MFEM `ParaViewDataCollection` and adds domain/boundary point-field maps.
- `DomainFieldEvaluator` can fill either an interpolatory L2 `ParGridFunction` or a full VTU point buffer in MFEM's element/refined-point traversal order.
- `PostOperator` registers derived fields as direct point fields for ParaView and still uses GridFunction evaluators for MFEM grid-function output.
- The proxy benchmark writes synthetic high-order scalar/vector fields and prints comparable ParaView/GridFunction times and byte counts.

Promising directions:

1. Avoid full-size temporary point buffers where possible. Current streaming still allocates complete point buffers before writing; a chunked writer could evaluate one geometry/rank chunk, write it, then reuse scratch.
2. Reduce duplicate host/device transfers and repeated `HostRead()` synchronization during VTU writing.
3. Avoid repeated refined-geometry and header work when MFEM already traverses the same element/order lattice.
4. Compare MFEM GridFunction save mechanics against VTU binary/base64/compression paths to isolate representation overhead.
5. Preserve the direct VTU point-field API shape: it should be easy to reason about and should not depend on fuzzy floating-point point lookup.

## How to run

```bash
./.auto/measure.sh
```

## Loop rules
1. First action: run `.auto/measure.sh` unmodified and log the proxy baseline before editing code.
2. Inspect source and benchmark logs before each change.
3. Keep atomic commits with clear metric deltas.
4. Discard regressions or correctness failures.
5. Update this file or `.auto/ideas.md` when you learn something that should survive context compaction.
