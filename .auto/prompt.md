# Autoresearch: Nsight-guided Palace AMR1 ParaView GPU postprocessing

## Objective
Optimize Palace/libCEED GPU postprocessing output for the large SingleTransmon AMR1 ParaView workload. Use NVIDIA profiling tools as the primary guide: Nsight Systems every measured iteration, and Nsight Compute selectively for suspicious kernels/memops.

The current default libCEED pointstream path is catastrophically slow for AMR output even though the legacy coefficient ParaView path is fast. The goal is to preserve GPU/libCEED postprocessing wins while reducing the AMR1 ParaView time and memory blowup.

## Metrics
- **Primary**: `paraview_seconds` (seconds, lower is better) — cumulative final AMR1 ParaView time from Palace's elapsed-time report under Nsight Systems.
- **Secondary**:
  - `total_seconds` — final cumulative Palace total time.
  - `postprocessing_seconds` — final cumulative postprocessing timer.
  - `operator_construction_seconds` — setup overhead; watch for excessive libCEED operator proliferation.
  - `paraview_memory_hwm_gb` — final high-water memory associated with ParaView phase.
  - `output_bytes` — don't win by dropping fields/output.
  - `nsys_report_bytes` / `nsys_stats_ready` — profiler health.
  - `cuda_log` / `ceed_gpu` — must both be 1.

## How to Run

```bash
./.auto/measure.sh
```

The measure script:
1. Builds Palace through Spack package-only install.
2. Creates a deterministic AMR1 config from `examples/transmon/transmon_amr.json` with:
   - `Solver.Order = 3`
   - `Model.Refinement.MaxIts = 1`
   - ParaView enabled
   - GridFunction disabled
   - GPU device
3. Runs `mpirun -n 8` with rank 0 under Nsight Systems and the other ranks unprofiled.
4. Emits `METRIC ...=...` lines.
5. Writes latest logs/profiles under `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/`.

Important files to inspect after each run:
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/palace.log`
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/nsys_stats.txt`
- `/home/ubuntu/scratch/palace_nsight_autoresearch/latest/profile/rank0*.nsys-rep`

## Profiling workflow

Use profilers before guessing. Do not add manual timing instrumentation as the first response to a performance mystery.

### Nsight Systems
The measurement already collects:

```bash
nsys profile --trace=cuda,nvtx,osrt,mpi --mpi-impl=mpich --osrt-file-access=true ...
```

Use `nsys stats` summaries to decide whether the bottleneck is:
- GPU kernels / libCEED point evaluation,
- CUDA memcpy / Unified Memory migration / synchronization,
- OS runtime/file writes,
- MPI waiting/imbalance,
- or CPU-side processing not visible as GPU work.

CPU IP/backtrace sampling may be disabled by kernel settings on this instance. That is okay: CUDA API/kernel/memcpy, OS runtime/file access, and MPI stats are still useful.

### Nsight Compute
Use `ncu` selectively only after Nsight Systems identifies suspicious kernels. Do not run full ncu on every experiment; it is too slow. Capture one rank / one or a few kernels when needed.

## Baselines and facts

Full AMR2 prior measurements:
- `origin/main`: ParaView ~220s, GridFunction ~188s, total ~983s, HWM ~11.6G.
- `auto2`: ParaView ~524s, GridFunction ~6.8s, total ~1117s, HWM ~27.2G.
- bad remote interleave: ParaView ~3658s, HWM ~135G.

AMR1/AMR2 pointstream sweep on this local branch before profiler-driven work:
- `local_pointstream_pv_only`: AMR0 PV ~18s, AMR1 PV ~1053s, AMR2 PV ~3636s, HWM ~149G.
- `local_legacy_gate_pv_gf` (`PALACE_LEGACY_SURFACE_POSTPRO=1`): AMR0 PV ~43s, AMR1 PV ~95s, AMR2 PV ~194s, HWM ~11.8G.
- Bulk contiguous `os.write` packing was tried for point fields. AMR1 remained ~1088s, so scalar-at-a-time writes were not the dominant bottleneck.

This means the bottleneck likely occurs before or around whole-domain pointstream evaluation/transfer/memory behavior, not simply at the final file write.

## Files in Scope

Primary:
- `palace/utils/ceedparaviewdatacollection.cpp/.hpp` — custom ParaView point-field writer; domain/boundary point evaluator registration and payload writing.
- `palace/fem/domain_point_field_evaluator.cpp/.hpp` — libCEED domain derived point-field evaluator (`U_e`, `U_m`, `S`).
- `palace/fem/point_field_evaluator.cpp/.hpp` — facade for domain/boundary visualization point fields.
- `palace/fem/ceed_group_operator.cpp/.hpp` — libCEED operator application helpers.
- `palace/models/postoperator.cpp/.hpp` — registration policy for ParaView/GridFunction output and backend selection.
- `palace/fem/output_functionals.cpp/.hpp` — boundary point/reduction backend; touch only if profiler points there.

Secondary/test files:
- `test/unit/test-postoperator-boundary-viz.cpp`
- `test/unit/test-surfacefunctional.cpp`

Autoresearch/session files:
- `.auto/measure.sh`
- `.auto/prompt.md`
- `.auto/ideas.md`

## Off Limits / hard constraints

- Use Spack only. Do not run direct CMake/make outside Spack.
- Do not optimize by dropping output fields, reducing AMR/refinement/order, changing rank count, or disabling ParaView data.
- Do not switch to `PALACE_LEGACY_SURFACE_POSTPRO=1` as the final answer unless explicitly documenting it as a fallback/gating strategy; the optimization target is GPU/libCEED postprocessing.
- Do not trust runs unless logs contain:
  - `Device configuration: cuda,cpu`
  - `libCEED backend: /gpu/cuda/magma`
- Preserve candidate output by overwriting latest output each iteration; do not accumulate large VTU outputs/profiles beyond the latest artifacts.
- p4d root EBS is persistent, but `/home/ubuntu/scratch` is instance-store scratch and can be wiped on stop.

## Current hypothesis backlog

Start with profiler evidence, then try targeted changes. Promising areas:

- Determine whether AMR1 time is GPU eval, CUDA memcpy/sync, UM migration, OS writes, or MPI imbalance using `nsys_stats.txt`.
- If GPU eval dominates: inspect libCEED AtPoints/buffer operator structure and reduce duplicated full-domain operator work.
- If memcpy/UM dominates: make host/device ownership explicit around `HostRead`, consider host staging, pinned buffers, or avoiding UM migration patterns.
- If OS writes dominate: keep bulk contiguous writes, check write sizes/counts in `osrt_sum`, avoid per-scalar or small-block writes.
- If full-domain scratch memory dominates: consider chunk-specific libCEED operators, but only after proving whole-domain buffers are the bottleneck. libCEED has no apply-time element range API; chunking requires chunk-specific restrictions/operators.
- If base `RegisterField` output dominates: explore pointstream/alternate writer for base E/B fields or hybrid policy.
- Hybrid fallback idea: use legacy coefficient ParaView for AMR domain derived fields while retaining libCEED GridFunction/reduction wins. Treat this as a possible guarded fallback, not the primary GPU pointstream optimization.

## What not to repeat

- Boundary-only interleaving does not solve the AMR domain ParaView regression.
- Whole-domain pointstream without understanding profiler output reproduces the ~1000s AMR1 / ~3600s AMR2 disaster.
- Bulk final `os.write` packing alone did not improve AMR1 pointstream timing.
