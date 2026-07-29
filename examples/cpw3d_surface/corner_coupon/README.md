# Three-dimensional corner response coupon

This directory generates local electrostatic response models for convex and
concave corner angles strictly between zero and 180 degrees. The coupon frame
is the frame used by Palace's automatic corner placement:

- the corner is at the origin;
- the first incident metal-edge arm points along positive `x`;
- the second arm is counterclockwise at the requested corner angle; and
- positive `z` points from the substrate into vacuum.

The fabricated prototype uses 100 nm metal, 50 nm substrate overetch, and 80
degree sidewalls. The thin and fabricated geometries can use either a sharp
plan-view corner or a circular fillet with a specified mask radius.

## Generate the prototype

From `examples/cpw3d_surface`:

```sh
mkdir -p corner_coupon/mesh
julia corner_coupon/mesh_corner_coupon.jl \
  convex-thin corner_coupon/mesh/corner_thin.msh 0.5
julia corner_coupon/mesh_corner_coupon.jl \
  convex-fabricated corner_coupon/mesh/corner_fabricated.msh 0.5

python3 corner_coupon/generate_corner_response.py \
  --output corner_coupon/calibration \
  --thin-mesh corner_coupon/mesh/corner_thin.msh \
  --fabricated-mesh corner_coupon/mesh/corner_fabricated.msh \
  --corner-radius 0.5 \
  --ring-size 8 \
  --order 1

../../build/bin/palace --serial corner_coupon/calibration/thin.json
../../build/bin/palace --serial corner_coupon/calibration/fabricated.json

python3 corner_coupon/finalize_corner_response.py \
  corner_coupon/calibration

../../build/bin/palace --serial \
  corner_coupon/calibration/heldout-thin.json
../../build/bin/palace --serial \
  corner_coupon/calibration/heldout-fabricated.json

python3 corner_coupon/finalize_corner_response.py \
  corner_coupon/calibration
```

The mesh generator imports `Gmsh.jl`. Use an instantiated Julia environment
containing that package. Omit the final mesh-generator argument and
`--corner-radius` to generate a sharp corner. Replace `convex-thin` and
`convex-fabricated` with `concave-thin` and `concave-fabricated`, and pass
`--topology concave` to `generate_corner_response.py`, to generate the
corresponding reentrant-corner entry. Pass the same `--angle ANGLE_DEGREES` to
both the mesh generator and response generator for a non-default angle. The
automatic coupon planner does this from Palace's `AngleDegrees` requirement
and uses at least quadratic mesh geometry so circular fillets remain curved.
A rounded coupon also requires `rho cot(theta / 2) < R`, so both fillet
tangency points lie inside the matching box, as required by runtime corner
recognition.

For the default process, the matching surface has nine closed square contours with
eight trace knots per contour. The resulting 72 coefficients include every corner
of the box matching surface and contours at the trench bottom, thin-metal plane,
and metal top. The coupon mesh resolves the intersections of these process planes
with the matching surface. Pass
`--metal-thickness` and `--overetch-depth` when those process dimensions differ
from the defaults. Each prescribed
potential is stored as a triangulated `x,y,z,V,triangle` surface trace.

`finalize_corner_response.py` checks the domain and surface matrices, aggregates
the per-edge response localized to the union of the physical-edge radius-`R`
tubes into compact matrices, and evaluates a smooth held-out boundary
excitation defined on a finer matching-surface triangulation. The held-out
trace vanishes smoothly throughout the metal-thickness band. It is therefore
compatible with both the thin and fabricated conductor cuts and does not
introduce an artificial, order-dependent Dirichlet singularity where the
matching surface meets the grounded metal.

Generated response configs enable `AggregateResponseMatrix`, so Palace performs the
physical-edge sum before storing `surface-response-matrix.csv`. The finalizer still maps
the radius-`R` core columns into the compact library format, but no longer needs a
hundreds-of-megabytes intermediate per-edge CSV.

For a rounded convex corner with angle `theta`, the generated library uses the
fillet center `(rho cot(theta / 2), rho, 0)` as the local PEC reference. This
reduces to `(rho, rho, 0)` at 90 degrees. The virtual sharp corner at the origin
is outside the rounded conductor and must not be used to gauge the coupon trace.
For a rounded concave corner, the origin remains on the conductor side of the
fillet and is the reference. The generator derives `ZeroTraceIndices` from the
swept envelope of the thin-plane and fabricated conductor cuts for sharp and
rounded corners of both topologies, so Maxwell contour knots on either PEC cut
are fixed to the reference-conductor potential. This conservative envelope
keeps a tapered cut compatible even when its sub-mesh-scale pullback lies
between trace knots.

At runtime Palace joins connected curved perimeter facets and estimates the
radius from the two tangent-arm setbacks and the corner angle,
`rho = mean(setback) tan(angle / 2)`. It rejects asymmetric joins and invalid
total turns before matching both corner angle and radius. This estimate is
insensitive to the number and chord lengths of the facets used to mesh a
circular fillet. A matched corner replaces the fillet and the remaining
distance `R - rho` along both straight arms.

When no model matches within its declared radius tolerance, Palace searches the
library's `CornerRadiusInterpolation` records for a held-out-qualified span with
the same topology and a compatible angle. Merely including lower- and
upper-radius models does not enable interpolation. Both radii must be positive;
the sharp coupon is not used as a rounded-radius endpoint, and Palace never
extrapolates. The response and conductor reference use convex interpolation
weights. The bracket width divided by `R` contributes to the reported maximum
library distance. If a sub-`R` fillet has neither an exact model nor a qualified
span, `UnmatchedPolicy: "Error"` terminates matching. With `Warn`, Palace retains
straight-edge response through that neighborhood and reports the unsupported
corner.

## Calibration cost and checks

There are no AMR iterations, and BoomerAMG is used for the electrostatic
solves. The error-estimator tolerance is intentionally relaxed because its
result is not used without AMR.

The automatic planner requires both FEM-order and mesh-resolution convergence
before merging a corner model. By default it compares p2 to p3 on the requested
mesh, then compares p3 responses at `2 * corner_lc_fine` and
`corner_lc_fine`. `--corner-h-factors` can add intermediate coarse-to-fine
levels but must contain at least two factors ending at one.

For `R = 2 um` and `rho = 0.5 um`, the refined convex thin/fabricated meshes
contain 331,606/397,672 tetrahedra; the corresponding concave meshes contain
332,156/395,665. Direct fabricated held-out energies changed as follows:

| Topology | Order change | Domain | SA | MS | MA |
|:--|:--|--:|--:|--:|--:|
| Convex | p1 to p2 | `-1.31%` | `+4.14%` | `+1.19%` | `+10.16%` |
| Convex | p2 to p3 | `-0.14%` | `+0.38%` | `-0.23%` | `+2.49%` |
| Concave | p1 to p2 | `-1.16%` | `+5.04%` | `+0.22%` | `+4.06%` |
| Concave | p2 to p3 | `-0.04%` | `+0.59%` | `-0.17%` | `+1.22%` |

The p1 convex 72-knot response matrix reproduced the independent compatible
trace with domain, SA, MS, and MA errors of `+2.07%`, `+0.07%`, `+5.31%`, and
`-4.49%`, respectively. This separates trace-basis error from FEM-order error.
Use `compare_coupon_convergence.py` with ordered `--case NAME=CALIBRATION`
arguments to reproduce both comparisons. Matrix conditioning and convergence
are evaluated only on the free trace subspace; PEC-constrained rows are not
physical response modes.

These tests validate the spatial trace representation and show convergence of
the fabrication-resolved surface energies for one smooth excitation. Final
process matrices should use p2 or higher. The tests do not by themselves
validate the physical corner correction in a device.

The generator also writes `probe-thin.json` and `probe-fabricated.json` for an
economical FEM-order convergence test. These configs use six smooth,
PEC-compatible outer traces instead of the 72 trace hats. Run both configs at
each order and compare the resulting projected response matrices with:

```sh
python3 corner_coupon/compare_probe_convergence.py \
  --case p1=/path/to/p1-calibration \
  --case p2=/path/to/p2-calibration \
  --case p3=/path/to/p3-calibration
```

The comparison reports the thin, fabricated, and fabricated-minus-thin domain
matrices, together with the thin and fabricated SA/MS/MA matrices. The matrix
change is a relative Frobenius norm. When a matrix is positive definite, the
worst-energy column is the maximum relative quadratic-energy change over all
linear combinations of the six probes. The domain defect is important because
cancellation can make normalization and the self-consistent correction converge
more slowly than either coupon domain response.

Corrected participation does not use a fabricated-minus-thin surface matrix.
Palace replaces the measured device energy inside `R` with the complete
fabricated surface response. The thin surface matrices are retained as coupon
diagnostics and library-format compatibility data, but their order sensitivity
does not enter the fixed-trace, fixed-flux, or self-consistent corrected surface
energy.

For the refined `R = 2 um`, `rho = 0.5 um` meshes above, the p2-to-p3 projected
matrix changes were:

| Topology | Matrix | Domain | SA | MS | MA |
|:--|:--|--:|--:|--:|--:|
| Convex | Thin | `0.04%` | `3.99%` | `5.05%` | `6.83%` |
| Convex | Fabricated | `0.09%` | `0.53%` | `0.27%` | `2.72%` |
| Convex | Domain defect | `1.59%` | - | - | - |
| Concave | Thin | `0.04%` | `5.25%` | `2.40%` | `3.19%` |
| Concave | Fabricated | `0.05%` | `0.55%` | `0.15%` | `1.39%` |
| Concave | Domain defect | `2.39%` | - | - | - |

The worst p2-to-p3 quadratic-energy changes over the probe space for the
fabricated responses were `0.25/1.03/2.32/4.92%` for convex
domain/SA/MS/MA and `0.16/1.14/1.24/4.99%` for concave. Thus p2 gives
percent-level dominant fabricated responses, while p3 is still preferable for
final corner libraries, particularly for weak MA combinations and the
domain-response defect. The six-source Palace totals were approximately
`30 s`, `3--4 min`, and `12--14 min` per coupon at p1, p2, and p3,
respectively, on the local Apple M3 Pro.

This projected test is suitable for local convergence studies; final process
libraries still require the full trace basis.

`run_probe_convergence.py` prepares and optionally runs these probe configs at
multiple FEM orders, then writes `probe-convergence.json`. By default the lower
order transitions are diagnostic and the acceptance gates apply to the final,
highest-order transition. Pass `--gate-all-transitions` to require every
adjacent-order comparison to pass.

```sh
python3 corner_coupon/run_probe_convergence.py \
  /path/to/calibration \
  --output /tmp/convex-r0p5-convergence \
  --palace ../../build/bin/palace \
  --orders 1 2 3
```

The default gates are 5% for fabricated response-matrix change, 10% for the
worst fabricated quadratic-energy change, and 5% for the
fabricated-minus-thin domain-defect matrix. Thin surface matrices are reported
but not gated because corrected participation uses the complete fabricated
surface response.

### Rounded-radius interpolation

Independent `rho = 0.25 um` and `rho = 0.75 um` calibrations were linearly
interpolated to the held-out `rho = 0.5 um` coupon. Relative errors for the
smooth held-out trace, compared with an independently calibrated
`rho = 0.5 um` response, were:

| Interface | Fixed trace | Fixed flux |
|:--|--:|--:|
| SA | `+0.31%` | `-0.15%` |
| MS | `+2.16%` | `+2.70%` |
| MA | `+1.36%` | `+1.05%` |

The individual response-matrix norm errors were `1.3--3.9%`. The interpolated
fabricated-minus-thin domain defect differed by `+16.5%` for this trace because
it subtracts two substantially larger domain responses. This cancellation
makes normalization-energy interpolation less robust than direct surface-energy
interpolation and remains visible through the fixed-trace/fixed-flux closure
diagnostic. Consequently this historical bracket does not pass the current
default 10% energy gate and must not be embedded as a qualified interpolation
span.

For three independently completed calibration directories, create the
held-out qualification report with:

```text
python3 corner_coupon/qualify_corner_interpolation.py \
  --lower corner_coupon/radius-0.25 \
  --heldout corner_coupon/radius-0.5 \
  --upper corner_coupon/radius-0.75 \
  --output corner_coupon/radius-interpolation.json
```

The qualifier checks thin and fabricated domain matrices, fabricated SA/MS/MA
matrices, the fabricated-minus-thin domain defect, and fixed-trace and
per-endpoint fixed-flux energies. Its default limits are 5% for matrix errors
and 10% for energy errors. It exits unsuccessfully and records `Passed: false`
when any checked quantity fails.

## Rectangular-island validation

> **Warning**
> The prototype device numbers in the validation sections below predate the
> swept-cut `ZeroTraceIndices` and compatible held-out trace. They are retained
> only as historical closure diagnostics and must be rerun before they can
> validate a library generated by the current scripts.

The scripts in `../corner_validation` compare a corrected zero-thickness
rounded rectangular metal island with an independently meshed
fabrication-resolved island. The process library combines the `rho = 0.5 um`
corner entry with the existing straight-edge entry. For a second-order device
solve, Palace detected all four fillets automatically. The participation-ratio
errors relative to the fabricated reference were:

| Interface | Raw | Fixed trace | Fixed flux | Self-consistent |
|:--|--:|--:|--:|--:|
| SA | `-1.99%` | `-6.09%` | `-1.19%` | `+2.24%` |
| MS | `-2.10%` | `-6.13%` | `+2.97%` | `-10.97%` |
| MA | `-53.03%` | `+1.72%` | `-13.09%` | `-17.16%` |

The fixed-trace/fixed-flux closure spread is 9.0--12.6%. Maxwell correction is
postprocessing-only, so those two columns, together with the reported closure
spread, are the relevant current estimate and confidence diagnostic.

### Driven-Maxwell validation

The same rounded island can be driven capacitively from its right edge. Generate matching
thin and fabricated meshes with:

```sh
julia corner_validation/mesh_corner_validation.jl \
  thin corner_validation/corner_maxwell_thin.msh 0.5 maxwell
julia corner_validation/mesh_corner_validation.jl \
  fabricated corner_validation/corner_maxwell_fabricated.msh 0.5 maxwell

python3 corner_validation/prepare_maxwell_validation.py \
  --output corner_validation/maxwell-p2 \
  --library corner_validation/process-library/process-library.json \
  --thin-mesh corner_validation/corner_maxwell_thin.msh \
  --fabricated-mesh corner_validation/corner_maxwell_fabricated.msh \
  --order 2 \
  --frequency 50
```

The fabricated mesh assigns the overetched SA footprint directly below the artificial
lumped-port sheet to boundary attribute 6. This footprint is excluded from both
validation metrics because the port replaces the corresponding SA patch in the thin
mesh and its metal-edge segment is excluded from correction.

For the second-order driven solve, the fabricated reference participation ratios and
relative thin-metal errors were:

| Interface | Fabricated reference | Raw | Fixed trace | Fixed flux | Self-consistent |
|:--|--:|--:|--:|--:|--:|
| SA | `5.2655e-4` | `-9.89%` | `-14.58%` | `-10.01%` | `-8.49%` |
| MS | `1.0980e-3` | `-6.30%` | `-10.24%` | `-1.21%` | `-16.85%` |
| MA | `1.9729e-5` | `-54.54%` | `-2.03%` | `-15.90%` | `-19.84%` |

The local Maxwell field passes the electrical-size and loop-closure checks, and all
straight and rounded edge neighborhoods match the process library. The maximum
fixed-trace/fixed-flux closure spread is `12.6%`, however, so Palace correctly reports a
failed confidence diagnostic. The result demonstrates accurate MA reconstruction with
fixed trace and MS reconstruction with fixed flux, but does not justify selecting one
closure universally or fitting an empirical blend to this geometry.

Replacing the exact `rho = 0.5 um` corner entry with the interpolated
`rho = 0.25/0.75 um` entries gave:

| Interface | Exact trace | Interpolated trace | Exact flux | Interpolated flux |
|:--|--:|--:|--:|--:|
| SA | `-14.58%` | `-13.90%` | `-10.01%` | `-9.00%` |
| MS | `-10.24%` | `-8.69%` | `-1.21%` | `+0.68%` |
| MA | `-2.03%` | `-0.05%` | `-15.90%` | `-13.85%` |

Thus radius interpolation changes the corrected device energies by about
`0.6--2.6%`, consistent with the direct coupon test. The remaining device error
is dominated by response-closure limitations rather than radius interpolation:
the interpolated run reports a `13.0%` closure spread and correctly fails the
confidence diagnostic.

## Rounded-aperture validation

The same validation generator can complement the island into a rounded
rectangular aperture. The fabricated reference uses a tapered opening through
the metal sheet and a substrate overetch opening that narrows toward the trench
bottom:

```sh
julia corner_validation/mesh_corner_validation.jl \
  thin corner_validation/aperture_thin.msh 0.5 aperture
julia corner_validation/mesh_corner_validation.jl \
  fabricated corner_validation/aperture_fabricated.msh 0.5 aperture
```

A process library containing only the generic isolated-edge model and the
independently calibrated rounded 90 degree concave-corner model automatically
matched all four aperture corners. No device boundary splitting or explicit
edge selection was used. The remote top and bottom planes are grounded while
the sidewalls use zero charge, avoiding an artificial Dirichlet discontinuity
where the metal sheet reaches the truncated device boundary.

For a second-order device solve with one AMR iteration, using the fine
175k- and 210k-element thin and fabricated corner coupons, the
participation-ratio errors were:

| Interface | Raw | Fixed trace | Fixed flux | Self-consistent |
|:--|--:|--:|--:|--:|
| SA | `-5.28%` | `-5.52%` | `+1.32%` | `+6.23%` |
| MS | `-0.65%` | `+0.09%` | `+2.39%` | `+0.25%` |
| MA | `-17.43%` | `+9.61%` | `-1.11%` | `-2.47%` |

The corrected thin participations changed by at most `0.03%` from AMR 0 to
AMR 1. The fabricated reference changed by `0.20%`, `0.16%`, and `0.03%` for
SA, MS, and MA, respectively. A first-order reference converged much more
slowly: after four AMR iterations its SA and MA values remained `4.7%` and
`6.6%` below the second-order AMR-1 result even though its bulk electric energy
was within `0.02%`. This is a surface-trace discretization effect, so the
second-order result is used as the reference above.

The maximum fixed-trace/fixed-flux closure spread is `21.5%`, on MA. Fixed flux
is the most uniform estimate for this aperture, but the closure diagnostic
correctly prevents promoting that observation into a geometry-independent
choice.

In a separate first-order AMR-2 ablation, removing the concave model increased
the fixed-trace SA error from `+0.81%` to `+26.34%`, the fixed-flux MA error
from `+4.88%` to `+12.62%`, and the self-consistent SA/MA errors from
`+7.02%/+9.52%` to `+29.75%/+22.12%`. The concave model is therefore necessary
for this geometry even though the remaining closure spread prevents a
production accuracy claim.

The generator accepts `--substrate-permittivity` and per-interface
`--sa-thickness`, `--sa-permittivity`, `--ms-thickness`,
`--ms-permittivity`, `--ma-thickness`, and `--ma-permittivity` options.
These values drive the coupon solves and are recorded in the version-3
`Fabrication.InterfaceLayers` metadata. They must match the straight-edge
coupons and the runtime dielectric postprocessing configuration.

### Rounded-corner probe convergence

For 100 nm metal, 50 nm overetch, and a `0.5 um` plan-view corner radius,
six-field p1/p2/p3 probe studies gave the following p2-to-p3 changes:

| Topology | Quantity | Matrix change | Worst energy change |
|----------|----------|--------------:|--------------------:|
| Convex | Fabricated domain | 0.09% | 0.24% |
| Convex | Fabricated SA | 0.53% | 1.03% |
| Convex | Fabricated MS | 0.27% | 2.32% |
| Convex | Fabricated MA | 2.72% | 4.92% |
| Concave | Fabricated domain | 0.05% | 0.16% |
| Concave | Fabricated SA | 0.55% | 1.14% |
| Concave | Fabricated MS | 0.15% | 1.24% |
| Concave | Fabricated MA | 1.39% | 4.99% |

The fabricated-minus-thin domain defect changed by 1.59% for the convex
coupon and 2.39% for the concave coupon. Both probe studies therefore pass
the default 5% matrix and 10% worst-energy convergence gates.

Full p2 matrices regenerated after adding all metal-interior trace knots to
`ZeroTraceIndices` pass the default 10% held-out reconstruction gate. The
errors below compare direct held-out solves against the quadratic forms
evaluated from the response matrices:

| Topology | Coupon | Domain | SA | MS | MA |
|----------|--------|-------:|---:|---:|---:|
| Convex | Thin | `+1.88%` | `-0.10%` | `+4.45%` | `-2.56%` |
| Convex | Fabricated | `+2.02%` | `+0.09%` | `+5.16%` | `-4.06%` |
| Concave | Thin | `+0.47%` | `+5.40%` | `-0.24%` | `+1.14%` |
| Concave | Fabricated | `+0.40%` | `+8.37%` | `-0.44%` | `+1.21%` |

The free-trace domain matrices are positive definite, with condition numbers
of 83 or less, and all interface matrices are positive semidefinite to the
numerical tolerance used by the qualifier. These tests qualify the rounded
90 degree, `0.5 um`-radius p2 coupon matrices for this trace basis and process
description; they do not replace the radius, angle, and nearby-edge coverage
still listed below.

### Acute and obtuse probe convergence

Rounded corner coupons with `R = 2 um`, `rho = 0.5 um`, and the default
process were tested at `lc_fine = 0.02 um`. Passing p2-to-p3 cases were also
tested at fixed p3 from `lc_fine = 0.04 um` to `0.02 um`:

| Topology | Study | Result | Domain defect | Fabricated domain | SA | MS | MA | Max. worst energy |
|----------|-------|:------:|--------------:|------------------:|---:|---:|---:|------------------:|
| 45 degree convex | p2 to p3 | Pass | `0.66%` | `0.07%` | `0.52%` | `0.26%` | `2.41%` | `5.16%` |
| 45 degree convex | h2 to h1 | Pass | `0.90%` | `0.03%` | `0.38%` | `0.71%` | `2.01%` | `5.70%` |
| 45 degree concave | p2 to p3 | **Fail** | `5.09%` | `0.05%` | `0.68%` | `0.17%` | `0.82%` | `3.45%` |
| 135 degree convex | p2 to p3 | Pass | `1.67%` | `0.09%` | `0.58%` | `0.15%` | `1.56%` | `3.97%` |
| 135 degree convex | h2 to h1 | Pass | `1.25%` | `0.02%` | `0.40%` | `0.50%` | `1.00%` | `3.22%` |
| 135 degree concave | p2 to p3 | Pass | `3.35%` | `0.07%` | `0.60%` | `0.15%` | `1.17%` | `2.79%` |
| 135 degree concave | h2 to h1 | Pass | `1.95%` | `0.01%` | `0.47%` | `0.34%` | `0.92%` | `2.90%` |

The matrix columns report relative Frobenius changes; the last column is the
largest fabricated worst-probe energy change across domain, SA, MS, and MA.
For the 45 degree concave case, the defect norm is only `0.666%` of the full
fabricated domain matrix. Its p2-to-p3 change is `0.034%` of that full matrix,
but remains `5.09%` of the correction itself and therefore correctly fails the
5% defect gate. A p3-to-p4 check and subsequent h-study remain required.

The passing studies validate the geometry and six-probe response convergence.
Full 72-trace matrices and independent held-out excitations have not yet been
run, so these angle families are not yet qualified process-library entries.

## Remaining validation

Before adding these matrices to a production process library:

1. Add any missing cross-sectional process rounding and fabrication details.
2. Check matching-radius and trace-basis convergence.
3. Generate at least two positive-radius entries, sweep held-out plan-view
   radii, and validate the interpolated response against fabricated references.
4. Converge the fabrication-resolved concave-corner device reference.
5. Combine the convex- and concave-corner entries with the process's isolated-
   and paired-edge entries.
6. Resolve the 45 degree concave p-order failure and converge additional corner
   angles and nearby-edge configurations.

The generic geometry and trace-mask paths also have valid high-order meshes for
sharp acute and obtuse coupons. Those sharp families still require the same p,
h, full-matrix, and held-out studies.

Each generated `process-library.json` contains only the selected corner
topology and radius. Combine the required entries before using the result as a
practical process library. A passed held-out report can be attached while
combining:

```text
python3 corner_coupon/combine_process_libraries.py \
  --output corner_coupon/process-library \
  --corner-interpolation-qualification \
    corner_coupon/radius-interpolation.json \
  corner_coupon/radius-0.25/process-library.json \
  corner_coupon/radius-0.75/process-library.json
```

The combiner rejects failed reports, missing endpoint models, duplicate spans,
and interpolation metadata in pre-version-3 libraries.
