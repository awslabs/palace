<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Graded coupon meshes: generality and acceptance

These tools are research/qualification tools. Generating a valid mesh, passing a
probe comparison, or preserving CAD measures does **not** automatically install
or qualify a response library.

## Accuracy targets are an error budget, not a universal 0.1% rule

`compare_coupon_probe.py` now has two explicit profiles:

| Profile | Domain energy | Fabricated interface energies |
|---|---:|---:|
| `engineering` (default, provisional) | 0.05% | 0.5% |
| `strict` (diagnostic) | 0.01% | 0.1% |

`--domain-tol` and `--surface-tol` override those values (relative fractions, so
`--surface-tol 0.01` means 1%). Reports store the actual tolerances and always
state `LibraryQualified: false`. The strict profile is useful only when the
reference and observable justify it. The retained three-edge reference's p4/p5
six-probe study still changed by percent-level amounts for some fabricated
interface quantities; it is not an established 0.1%-accurate continuum reference.

The engineering targets are starting points, not an accumulated device-error
bound. Required checks remain:

- thin domain response, **not edge-inclusive raw thin-metal SPR**;
- fabricated interface response for several independent excitations;
- original basis and conductor-state probes, including difficult localized modes;
- the domain defect `F-T`, its global embedding/stability, and final observables;
- independent mesh/order checks when reference uncertainty is appreciable.

A cancellation-sensitive defect may require a different absolute/relative scale
than either domain matrix. Tolerance relaxation does not justify dropping modes,
interfaces, overlap checks, or convergence checks.

## Geometry is a separate input

Pass `--process graded-process-example.toml` to `mesh_graded_tet_experiment.jl`.
The file explicitly declares units, radius, metal thickness, overetch, wall angle,
and rounding. Dimensions are currently in micrometres; unknown settings, missing
file units, and non-finite inputs fail. The resolved process is saved beside the
mesh as `.process.toml`. Meshing sizes and FEM solution order are separate from
these dimensions.

The supported implementation is not based on model names or edge counts:

- polygonal masks, concave boundaries, multiple conductors and slots;
- arbitrary **in-plane** orientations and translations;
- holes, with material-side-aware offsets independent of input winding;
- parallel/opposed process planes, including repeated conductor labels;
- straight CAD curves even when OCC represents them as Bezier curves;
- exact circle distance sizing and validated conic models for generated fillets.

For a conic size model, the polygonal distance approximation is biased toward
refinement using its interpolation-error bound (at most 1% of the requested
minimum size). It does not replace the CAD curve with a polygon. Trimmed conics
are checked against sampled positions and derivatives; this is not a general
B-spline certification algorithm.

### Explicit limits

- Substrate half-spaces and fabricated layer bands must not overlap.
- General tilted/nonparallel planes are not supported by this coupon generator.
- Nonconic nonlinear curves fail rather than being replaced with straight chords.
- Nonzero offsets of curved plan-view boundaries and shrinking nonconvex holes
  require additional topology-aware support and fail explicitly.
- A hole that disappears between loft levels is rejected. A convex hole that
  fully disappears at both offset levels is handled as a collapsed void, not a
  reflected polygon.
- Current attribute encoding permits slots 0–9 and conductor IDs 1–99; these are
  label-format limits, not a limit on the number of edges.
- General multi-slot surface partitions still need energy-weighted convergence
  checks. Centroid ownership is exclusive, but its sampled ambiguity is only a
  diagnostic. Do not confuse a valid label set with an accurate partition. The
  opt-in quadrature partition pilot is described in
  [quadrature ownership results](quadrature-ownership-results-20260910.md); it
  resolves polygonal ownership at integration points without further PDE solves.

## Bugs exposed by the broader cases

- A CAD bounding box's padding could cause a requested submicron fillet to be
  silently skipped. Curve coordinates are now used to select process-plane edges;
  a nonzero requested fillet with no eligible edges fails explicitly. Tests also
  verify that rounding actually changes the solid, rather than only trusting the
  request metadata.
- Hole offset direction depended on point winding, and an over-eroded hole could
  reopen as a reflected polygon. Material-side orientation and convex half-plane
  clipping now handle those cases. An independent trench-volume test verifies
  that etching reaches the hole while retaining the metal footprint.
- Whole CAD-face labels cannot provide the expected multi-slot partition. The
  scout suite uses a separately derived expected-label contract and checks every
  expected label after element-wise assignment.
- Nearest-edge ties now prefer stable physical labels instead of CSV row order.
- Layer ownership is selected before edge/slot ownership, preventing a nearby
  edge on the opposing layer from capturing a surface patch.

## Repeatable generality scouts

Build the optional audit tools against an existing Palace/MFEM installation:

```sh
make -f examples/cpw3d_surface/spatial_coupon/audit-tools.mk \
  PALACE_BUILD=/path/to/palace/build AUDIT_OUTPUT=/tmp/coupon-audit
```

Run the suite with a Julia environment that contains Gmsh:

```sh
python3 examples/cpw3d_surface/spatial_coupon/run_general_mesh_suite.py \
  --root /tmp/new-general-coupon-study \
  --julia-project /path/to/julia/environment \
  --audit-bin /tmp/coupon-audit/audit_coupon_mesh \
  --measures-bin /tmp/coupon-audit/audit_mesh_measures
```

It tests strips, an L/concave mask, its rotated/translated form, a T junction,
two and three conductors, a hole, opposing layers, rounding, and sloped walls.
Each mesher is separately time/RSS bounded, with one mesher active at a time.
Tools/binaries are snapshotted and hashed; a new output root is required.

Checks include independent metal-area formulas, physical group coverage,
exterior/material-interface face coverage, positive quadrature Jacobians,
integrated mesh/CAD volumes and areas, and quadrature-order comparison. Curved
area agreement has a 0.1% scout target; affine geometry uses a 0.0001% target.
These are geometric bookkeeping checks, not local normal-error or response-energy
bounds. The quality auditor's legacy process-distance histograms are not used for
cases with different process dimensions.

`run_mesh_rotation_probe.py` adds a small p5 multi-conductor solve before/after a
rigid transform of identical connectivity. Its tight covariance tolerance tests
coordinate invariance, not mesh convergence. It does not require an HPC allocation.

The BAMG surface-metric and fixed-shell hybrid prototypes remain experimental.
Their failures are recorded in `mesh-followup-20260909.md`; neither is silently
selected as a production fallback.
