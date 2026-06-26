# Ideas backlog

- Use Nsight Systems `osrt_sum` / syscall summaries to verify whether writes are still many/small after bulk payload packing. If file IO still dominates, inspect actual write sizes and consider bypassing iostream buffering or using POSIX `write`/`writev` for appended payloads.
- If CUDA memcpy/sync/UM dominates, make pointfield host staging explicit: evaluate to device, copy to pinned host buffer, reshape there; avoid accidental Umpire/UM migration of huge vectors.
- If GPU kernels dominate, inspect libCEED AtPoints kernels with Nsight Compute on rank 0 for the AMR1 output phase; check occupancy, memory throughput, and whether evaluating full-domain U_e/U_m/S separately repeats expensive geometry/material work.
- If operator construction dominates, lazily assemble only the buffer path for ParaView-only runs; avoid assembling both GridFunction and point-buffer operators when GridFunction output is disabled.
- Explore fusing domain derived point fields (`U_e`, `U_m`, `S`) into one libCEED point evaluation pass so E/B/material/geometry are loaded once and multiple outputs are written together. This needs VTK multi-array handling or temporary tuple buffers.
- If whole-domain scratch memory is still the issue, design chunk-specific libCEED restrictions/operators. libCEED has no apply-time element subset/range parameter; chunking requires chunk-specific restrictions/operators.
- Hybrid fallback: gate AMR/nonconforming ParaView domain derived fields back to legacy coefficient output while retaining libCEED for GridFunction, reductions, CPW, and boundary cases. Use if profiling shows pointstream cannot be fixed quickly.
