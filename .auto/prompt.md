# Autoresearch: Nsight-guided Palace AMR1 volume ParaView GPU postprocessing

## Objective
Optimize Palace/libCEED **volume/domain** ParaView postprocessing for the large SingleTransmon AMR1 workload. Keep the full AMR1 run and full output enabled, but use an env-gated profile line to measure the domain `paraview->Save()` calls separately from boundary `paraview_bdr->Save()`.

The previous segment fixed the no-fallback nonconforming boundary trace explosion, taking full cumulative ParaView time from 1151s to a best observed 48.35s. This new segment should focus on volume/domain postprocessing work: base domain fields (`E`, `B`) and derived domain point fields (`U_e`, `U_m`, `S`), not boundary trace-map grouping.

## Metrics
- **Primary**: `volume_paraview_seconds` (seconds, lower is better) — sum of max-over-ranks domain `paraview->Save()` times from `VolumeProfile step=... domain_save_max=...` lines. This excludes boundary save time and final rank/indicator save.
- **Secondary**:
  - `paraview_seconds` — Palace cumulative full ParaView timer; do not improve volume by making full output worse dramatically.
  - `boundary_paraview_seconds` — boundary save time from profile lines, to detect accidental boundary regressions.
  - `final_volume_paraview_seconds` — final domain rank/indicator save time.
  - `total_seconds`, `postprocessing_seconds`, `operator_construction_seconds`.
  - `paraview_memory_hwm_gb`, `output_bytes` — don't win by dropping fields/output.
  - `nsys_report_bytes`, `nsys_stats_ready`, `cuda_log`, `ceed_gpu`, `run_rc`.

## How to Run

```bash
PALACE_CEED_NONCONFORMING_PARAVIEW=1 PALACE_SURFACE_PROFILE=1 ./.auto/measure.sh
```

The measure script:
1. Builds Palace through Spack package-only install.
2. Creates a deterministic AMR1 config from `examples/transmon/transmon_amr.json` with:
   - `Solver.Order = 3`
   - `Model.Refinement.MaxIts = 1`
   - ParaView enabled
   - GridFunction disabled
   - GPU device
3. Runs `mpirun -n 8` with rank 0 under Nsight Systems and other ranks unprofiled.
4. Enables `PALACE_VOLUME_PROFILE=1` to print domain/boundary save split timings.
5. Emits `METRIC ...=...` lines.
6. Writes latest logs/profiles under `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/`.

Important files to inspect after each run:
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/palace.log`
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/nsys_stats.txt`
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/profile/allranks*.nsys-rep`

## Profiling workflow
Use NVIDIA profiling tools as the guide. Nsight Systems is collected every measured iteration. Use Nsight Compute selectively only after Nsight Systems points to suspicious kernels/memops.

CPU sampling may be unavailable; CUDA API/kernel/memcpy, OS runtime/file access, MPI, and the env-gated `VolumeProfile` split are the main evidence.

## Files in Scope
Primary volume/domain files:
- `palace/utils/ceedparaviewdatacollection.cpp/.hpp` — custom ParaView point-field writer, appended payload packing, domain point evaluator registration.
- `palace/fem/domain_point_field_evaluator.cpp/.hpp` — libCEED domain derived point-field evaluator for `U_e`, `U_m`, `S`.
- `palace/fem/point_field_evaluator.cpp/.hpp` — facade for domain/boundary point fields.
- `palace/models/postoperator.cpp/.hpp` — output registration policy and profile split instrumentation.
- `palace/fem/ceed_group_operator.cpp/.hpp` — libCEED apply helpers if profiler points there.

Secondary/reference:
- `palace/fem/output_functionals.cpp/.hpp` — boundary path; avoid touching unless a volume change needs boundary interaction.
- `test/unit/test-postoperator-boundary-viz.cpp`, `test/unit/test-surfacefunctional.cpp`.

Autoresearch/session files:
- `.auto/measure.sh`, `.auto/prompt.md`, `.auto/ideas.md`, `.auto/log.jsonl`.

## Off Limits / hard constraints
- Use Spack only. Do not run direct CMake/make outside Spack.
- Do not optimize by dropping domain output fields, reducing AMR/refinement/order, changing rank count, or disabling ParaView data.
- Keep full output enabled during measurement. The primary metric isolates domain save time via profiling, not by skipping boundary output.
- Do not switch to `PALACE_LEGACY_SURFACE_POSTPRO=1` as a final answer.
- Do not trust runs unless logs contain:
  - `Device configuration: cuda,cpu`
  - `libCEED backend: /gpu/cuda/magma`
- Preserve candidate output by overwriting latest output each iteration; do not accumulate large VTU outputs/profiles beyond latest artifacts.

## Previous segment lessons
Kept boundary innovations:
- Finite MFEM/NCMesh NC trace-map grouping collapsed per-element mapped boundary operators.
- Split local/ghost NC traces so local side used AtPoints.
- Grouped ghost-only local accumulation independently of remote trace-map identity.
- Used AtPoints export evaluators for face-neighbor ghost field exchange.
- Evaluated continuous boundary `E_t`/`B_n` from MFEM boundary/face trace DOFs.
- Cached marked boundary trace indices.

Failed/deferred ideas:
- Metadata-only writer refined-point-count caching regressed.
- Trace zero-fill avoidance and trace IR pointer caching regressed or were noise.
- Thread-local payload buffer reuse regressed.
- Larger VTU stream buffers reduced write syscall count but regressed primary time; file writes were not the bottleneck.
- Skipping unused domain GridFunction evaluator assembly improved `operator_construction_seconds`/HWM but regressed the old full ParaView metric; revisit only if construction becomes the target.

## Current hypothesis backlog for volume/domain
- First establish a volume split baseline with `PALACE_VOLUME_PROFILE=1` and inspect `volume_paraview_seconds` vs boundary/final times.
- If domain derived fields dominate: explore fusing `U_e`, `U_m`, and `S` into fewer libCEED point evaluation passes or a multi-field point evaluator/writer API, so E/B/geometry/material data are loaded once.
- If base `RegisterField` output dominates: inspect MFEM domain GridFunction VTU writing and consider a pointstream/alternate writer for base E/B fields, preserving output fields and layout.
- If CUDA memcpy/UM dominates: make host/device ownership explicit around point buffers and `HostRead`; consider host staging/pinned buffers only with profiler evidence.
- If OS writes dominate despite prior stream-buffer failure: inspect actual domain vs boundary write sizes/counts before changing writer code.
- If operator construction dominates for volume metric: revisit buffer-only domain evaluator assembly, but require a primary `volume_paraview_seconds` win and clean implementation.

## What not to repeat
- Do not repeat boundary trace-map work for volume; volume elements do not have face-neighbor trace-map identity.
- Do not repeat metadata-only writer caching, thread-local payload reuse, or ofstream buffer enlargement without new profiling evidence.
- Do not disable boundary output to define a fake win; use the profile split while preserving full output.
