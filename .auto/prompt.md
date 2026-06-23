# Autoresearch: Palace GPU SurfaceFunctional reduction performance

## Objective
Improve Palace GPU postprocessing performance by moving/keeping 3D statics + driven surface reductions on efficient libCEED kernel paths with much less setup/JIT/operator proliferation.

The current bottleneck is not the arithmetic in the surface integral QFunctions. Profiling showed the integral kernels take only milliseconds while postprocessing takes tens of seconds. The issue is that mapped volume-reference face quadrature coordinates are baked into `SurfaceFunctional` grouping and basis/operator specialization, creating many libCEED operators and NVRTC compiles for equivalent trace configurations.

The first target is local `SURFACE_FLUX` reductions in `palace/fem/surfacefunctional.cpp`. Once that is working and measured, reuse the same local low-level context for `INTERFACE_EPR`, then `FARFIELD`.

## Branch
Work on this branch only:

```text
hughcars/libceed-output-functionals-dev-auto
```

It was created from the current tip of:

```text
hughcars/libceed-output-functionals-dev
```

Commit each kept experiment to the `-auto` branch so the improvement history is visible. The non-auto branch is the comparison baseline.

## Primary metric
- **surface_score_seconds** (seconds, lower is better):
  `spheres_postprocessing + rings_postprocessing + cpw_postprocessing`

## Secondary metrics
Always report/consider:

- `spheres_postprocessing`, `rings_postprocessing`, `cpw_postprocessing`
- `cpw_farfields`
- `spheres_nvrtc`, `rings_nvrtc`, `cpw_nvrtc`
- `spheres_cuModuleLoadData`, `rings_cuModuleLoadData`, `cpw_cuModuleLoadData`
- `spheres_surface_flux_groups`, `rings_surface_flux_groups`, `cpw_surface_flux_groups`
- `cpw_interface_epr_groups`, `cpw_farfield_groups`
- `spheres_total`, `rings_total`, `cpw_total`

## Baseline cases
The measurement script runs these on the remote GPU host:

1. `examples/spheres/spheres.json` — electrostatic `SURFACE_FLUX`
2. `examples/rings/rings.json` — magnetostatic `SURFACE_FLUX`
3. `examples/cpw/cpw_lumped_uniform.json` — driven rich postprocessing, one driven point at `17 GHz`, field output saved

Starting reference from prior measurements:

- `spheres`: postprocessing ~3.5 s, `nvrtcCompileProgram` ~254, `SURFACE_FLUX` groups `7` and `6`
- `rings`: postprocessing ~33 s, `nvrtcCompileProgram` ~1830, `SURFACE_FLUX` groups `37` and `53`
- `cpw_lumped`: postprocessing ~23 s, far fields ~3 s, `nvrtcCompileProgram` ~1610
- Temporary grouping instrumentation showed these collapse to one coarse group when mapped coordinates are ignored:
  - `spheres SURFACE_FLUX`: `7/6 -> 1/1`
  - `rings SURFACE_FLUX`: `37/53 -> 1/1`
  - `cpw_lumped SURFACE_FLUX`, `INTERFACE_EPR`, and `FARFIELD`: each line -> `1`

## Files in scope
Primary low-level scope:

- `palace/fem/surfacefunctional.cpp`
- `palace/fem/surfacefunctional.hpp`
- `palace/fem/qfunctions/32/surf_32_qf.h`
- `palace/fem/libceed/functional.cpp`
- `palace/fem/libceed/functional.hpp`

Tests/guards:

- `test/unit/CMakeLists.txt`
- `test/unit/test-*.cpp`

Later batching stage only, after low-level grouping works:

- `palace/models/surfacepostoperator.cpp`
- `palace/models/surfacepostoperator.hpp`
- `palace/models/lumpedportoperator.cpp`
- `palace/models/waveportoperator.cpp`

## Off limits
Do not optimize or edit:

- linear solver/preconditioner/coarse solver paths
- BoundaryMode, transient, or 2D fallback paths unless required to keep compilation working
- example meshes
- Spack environment files, generated profiling JSONs/logs, `postpro/`, build directories
- libCEED, unless a Palace-side change proves a missing libCEED primitive is strictly required

## Hard constraints
- Use Spack only. Never run direct CMake/make.
- Do not exceed `-j4`.
- Remote build command is:
  ```bash
  source /home/ubuntu/spack/share/spack/setup-env.sh
  spack -e /home/ubuntu/palace install -j4
  ```
- Preserve correctness: capacitance/inductance/S-parameters/participations/far-field scalar outputs should remain consistent with baseline logs.
- Keep changes atomic. Commit each kept experiment with a clear metric-oriented message.
- Prefer reducing group/JIT/operator proliferation over moving cost into heavier kernels.

## Technical strategy hints
The current fine grouping includes mapped volume-reference quadrature coordinates. The implementation should make mapped trace points runtime data rather than JIT key data while preserving 2D surface integral semantics.

When considering an “AtPoints-backed reduction,” treat it only as an implementation mechanism for runtime trace-point basis evaluation. The operation remains a 2D surface integral reduction, not ParaView/gridfunction point output.

Surface-only quantities should not participate in volume trace specialization:

- surface-only: `dS`, unit normal `n`, oriented weighted normal `wN`, physical coordinate `x`, normal flip, `x0` outward orientation
- still requires adjacent cell trace: H(curl) `E`, H(div) `B`, Piola transforms, material side attributes, two-sided/ghost semantics

## Suggested order
1. Collapse local `SURFACE_FLUX` mapped-reference-coordinate grouping.
2. Extend the same local collapsed path to `INTERFACE_EPR`.
3. Extend the same local collapsed path to `FARFIELD`.
4. Then handle ghost/two-sided parallel interfaces.
5. Then consider batching same-kind surface postprocessing requests in `models/`.

Surface geometry precompute (`w`, `n`, `wN`, `x` EVAL_NONE arrays) is useful cleanup and may be folded into the grouping work if it simplifies implementation, but it is not by itself the main performance milestone.

## Keep/discard policy
Keep only if:

- `.auto/measure.sh` exits 0,
- `.auto/checks.sh` passes,
- `surface_score_seconds` improves, or a targeted submetric improves substantially with no meaningful regression,
- and the direction matches the current stage goal (e.g. group count collapse for `SURFACE_FLUX`).

Discard if:

- correctness changes,
- group counts do not improve for the targeted functional,
- `nvrtcCompileProgram` or postprocessing time regress materially,
- or the code becomes broad/fragile without a metric win.

## Resuming
On every resume:

1. Read this file.
2. Read `AGENTS.md`.
3. Skim `surface_functional_kernel_proliferation_report.md` if present.
4. Read `.auto/log.jsonl` tail and `git log --oneline --decorate -20`.
5. Continue with the next compact experiment; do not restart broad surveys.

## What's been tried
Prior profiling showed the surface integral QFunction GPU time is tiny: rings flux kernels were ~2.3 ms while postprocessing was ~33 s; CPW EPR/flux/farfield kernels were under ~10 ms while postprocessing was ~23 s. Therefore optimize setup/JIT/operator proliferation, not arithmetic throughput.

## Manual baseline before autoresearch loop
A manual `.auto/measure.sh` run during setup succeeded with:

```text
METRIC surface_score_seconds=76.767000
METRIC spheres_postprocessing=6.085000
METRIC spheres_nvrtc=254
METRIC spheres_surface_flux_groups=13
METRIC rings_postprocessing=39.405000
METRIC rings_nvrtc=1830
METRIC rings_surface_flux_groups=90
METRIC cpw_postprocessing=31.277000
METRIC cpw_farfields=3.478000
METRIC cpw_nvrtc=1610
METRIC cpw_surface_flux_groups=45
METRIC cpw_interface_epr_groups=23
METRIC cpw_farfield_groups=8
```

Use the next `run_experiment` baseline as the official autoresearch baseline, but these values are a sanity check for whether measurement is working.
