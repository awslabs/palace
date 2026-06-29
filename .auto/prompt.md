# Autoresearch: Transmon conforming-AMR GPU postprocessing OOM

## Objective
Make the GPU/libCEED postprocessing path complete `examples/transmon/transmon_amr.json` when the generated run uses conforming AMR:

```json
"Model": { "Refinement": { "Nonconformal": false, "MaxIts": 1 } }
```

The run must keep GPU execution and both ParaView and GridFunction output enabled. The current failure happens after one conforming refinement, during boundary point-field postprocessing in `fem::ApplyAddGroupOperators`, with libCEED CUDA OOM from `CeedSetDeviceGenericArray_Cuda`.

## Metrics
- **Primary**: `transmon_ncfalse_penalty_s` (s, lower is better). Successful runs report the Palace total time; OOM runs report `999999` as a failure penalty.
- **Secondary**:
  - `transmon_ncfalse_oom`: 1 for CUDA OOM, 0 for successful completion.
  - `transmon_ncfalse_exit`: Palace process exit code.
  - `transmon_ncfalse_initial_elements`: initial mesh element count.
  - `transmon_ncfalse_refined_elements`: post-conforming-refinement element count.
  - `transmon_ncfalse_bdr_elems_max`: largest boundary SurfaceFunctional element count reached before completion/failure.
  - `transmon_ncfalse_initial_paraview_s`, `transmon_ncfalse_initial_gridfunction_s`: initial-mesh output timings if available.
  - `transmon_ncfalse_nvrtc`, `transmon_ncfalse_cuModuleLoadData`: CUDA/NVRTC counts if the process exits cleanly enough for the preload destructor.

## How to Run

```bash
./.auto/measure.sh
```

The script rebuilds Palace through Spack (`-j4`), generates `/tmp/palace_auto_transmon_ncfalse_<tag>.json`, runs the installed binary with `PALACE_SURFACE_PROFILE=1`, and prints `METRIC name=value` lines.

For bootstrapping from an already captured log, use:

```bash
PALACE_MEASURE_EXISTING_LOG=/tmp/transmon_amr_branch_ncfalse_maxits1_full.log ./.auto/measure.sh
```

## Files in Scope

Primary Palace-side files likely involved in memory ownership or peak device memory:

- `palace/fem/ceed_group_operator.{hpp,cpp}` — shared libCEED group apply helper; current OOM reports here. Good place for optional diagnostics, output-vector lifetime changes, or apply-time memory reduction.
- `palace/fem/output_functionals.{hpp,cpp}` — `SurfaceFunctional` assembly, AtPoints boundary groups, retained group data, local output buffers.
- `palace/fem/point_field_evaluator.{hpp,cpp}` — boundary visualization point-field wrapper; currently can release/rebuild boundary evaluators to reduce retained AtPoints operators.
- `palace/models/postoperator.{hpp,cpp}` — controls when field-output evaluators are assembled and retained.
- `palace/drivers/eigensolver.cpp` — already releases `ksp`, `P`, and `A` after eigensolve to reduce postprocessing peak memory.
- `palace/fem/domain_point_field_evaluator.cpp`, `palace/fem/facenbrexchange.cpp`, `palace/fem/interpolator.cpp` — related users of `CeedGroupOperator` and cleanup helper.

Diagnostics may also temporarily touch libCEED CUDA allocation files under `/home/ubuntu/libceed-atpoints-poc`, but do not keep libCEED diagnostic patches in final Palace commits unless explicitly justified.

## Off Limits

- Do not disable ParaView, GridFunction, GPU execution, or boundary-field output to make the benchmark pass.
- Do not change the benchmark physics or silently reduce mesh/refinement beyond the specified `Nonconformal=false, MaxIts=1` reproducer.
- Do not broaden into unrelated BoundaryMode/transient/2D behavior.
- Do not bypass Spack builds; always build/test via `source /home/ubuntu/spack/share/spack/setup-env.sh` and `spack -e /home/ubuntu/palace install -j4` or overwrite equivalent.
- Never build with more than `-j4`.
- Do not leave broad temporary libCEED hacks or diagnostic prints enabled by default.

## Current Reproducer and Evidence

Known failing config/log:

- `/tmp/transmon_amr_branch_ncfalse_maxits1.json`
- `/tmp/transmon_amr_branch_ncfalse_maxits1_full.log`

Key facts from that log:

```text
Marked 1814/62950 elements for refinement (70.00% of the error, θ = 0.70)
Conforming mesh refinement added 148645 elements (initial = 62950, final = 211595)
...
SurfaceFunctional profile kind=BDR_ENERGY_E groups=5 elems=37364
SurfaceFunctional profile kind=BDR_FIELD_E groups=5 elems=37364
SurfaceFunctional profile kind=BDR_FLUX_Q groups=4 elems=37364
SurfaceFunctional profile kind=BDR_ENERGY_M groups=5 elems=37364
SurfaceFunctional profile kind=BDR_FIELD_B groups=5 elems=37364
SurfaceFunctional profile kind=BDR_CURRENT_J groups=4 elems=37364
SurfaceFunctional profile kind=BDR_POYNTING groups=5 elems=37364
MFEM abort: ... CeedSetDeviceGenericArray_Cuda(): CUDA_ERROR_OUT_OF_MEMORY
 ... in palace::fem::ApplyAddGroupOperators(...)
 ... palace/fem/ceed_group_operator.cpp:75
```

`MaxIts=1` still solves/postprocesses the refined mesh; it does not avoid the failure.

`Nonconformal=true, MaxIts=2` completed after adding the eigensolver `ksp/P/A` release, but `Nonconformal=false, MaxIts=1` still OOMs because the first conforming refinement grows the mesh to 211595 elements and boundary point-field groups to 37364 boundary elements.

## What's Been Tried

- Fixed shared error-estimator MAGMA failure by avoiding fake basis reduction in `AssembleCeedElementErrorIntegrator`; uses duplicate-index `CEED_EVAL_NONE` restriction instead.
- Added current-style `DestroyGroupOperators(...)` and fixed a concrete `CeedOperatorFieldGetVector(...)` reference leak in `ApplyAddGroupOperators`.
- Avoid retaining unused non-FARFIELD `CeedQFunctionContext` handles.
- Made local `CeedAssemblyScratch` structs non-copyable.
- Added lazy field-output initialization in `PostOperator` so field-output libCEED objects are not created in the constructor.
- Added boundary `PointFieldEvaluator` release/rebuild logic; this reduces retained boundary evaluator objects but increases repeated JIT/rebuild work.
- Released eigensolver `ksp`, `P`, and `A` after eigenvector rescaling. This allowed the nonconformal full AMR transmon case to complete, but not the larger conforming-refinement case.
- Tried retaining boundary evaluators again for performance; conforming AMR still OOMed. Restored release-after-use.

## Immediate Investigation Plan

1. Attribute device memory around `fem::ApplyAddGroupOperators` on the conforming refined mesh:
   - before each boundary kind/group apply,
   - before/after `CeedVectorCreate`, `CeedVectorSetArray`, and `CeedOperatorApplyAdd`,
   - output vector sizes and group labels.
2. Identify whether peak memory is dominated by:
   - retained libCEED group operators/restrictions/bases,
   - per-kind full boundary output buffers,
   - MFEM vector device aliases created by `out.ReadWrite(true)`,
   - libCEED internal active vector/basis work arrays,
   - duplicate boundary evaluators across modes/real-imag writes.
3. Prefer fixes that lower peak memory without disabling outputs, e.g. destroy per-kind boundary evaluators immediately after their last use, process boundary fields one kind/mode at a time, avoid retained `out_vec` wrappers if they hold owned device storage unexpectedly, or stream smaller output chunks if feasible.

## Loop Rules for This Session

- Baseline failing/OOM metric is `999999`; a successful full `Nonconformal=false, MaxIts=1` run is a major improvement and should be kept if outputs are sane.
- If an experiment still OOMs, discard it unless it provides diagnostic-only evidence; record the exact attribution in ASI before rollback.
- Keep changes focused and upstreamable. Prefer manual cleanup helpers matching existing style over broad RAII refactors.
