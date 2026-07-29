```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Thin-metal surface-response roadmap

## Scope

This document defines the work required to make fabrication-response-corrected
SA, MS, and MA participation ratios a production feature for thin-metal simulations
of superconducting quantum circuits. It covers electrostatic, boundary-mode, driven,
and eigenmode simulations on two- and three-dimensional meshes.

The correction replaces the nonconvergent interface energy near an ideal zero-thickness
metal edge with a response computed from a local fabrication-resolved coupon. The device
mesh remains thin. A process library records metal thickness, overetch, taper, rounding,
interface films, and the qualified response matrices needed by the local geometries found
in the device.

This roadmap separates three levels of readiness:

  - **Implemented:** The code path and focused tests exist.
  - **Qualified:** Held-out coupon or device tests establish an accuracy envelope.
  - **Production:** The complete workflow satisfies the acceptance gates below and fails
    closed outside its qualified envelope.

The feature is currently experimental. Core operators and geometry matching are
implemented, but the complete practical-geometry validation matrix has not passed.

## Definition of feature complete

The feature is complete for production when all of the following hold:

  1. A user supplies the same thin-metal device mesh and SA/MS/MA interface attributes
     used by raw postprocessing, plus a fabrication-process description or process
     library. No correction-specific edge or corner attributes are required.
  2. Geometry preflight classifies every correction neighborhood before the field solve.
     It automatically retrieves or generates qualified coupons, or reports an exact
     unsupported signature and exits according to `UnmatchedPolicy`.
  3. The supported solver and boundary-condition matrix below is complete. Raw fields,
     raw participation, and the raw AMR sequence remain unchanged and available.
  4. Isolated edges, close parallel edges, corners, smooth bends, trace caps and branch
     intersections, and nonparallel or cross-layer interactions are either corrected by
     a qualified local model or explicitly rejected. No geometric fallback silently
     invents a fabricated conductor.
  5. Every reported corrected value carries coverage, trace/flux closure, electrical-size,
     boundary-law, curvature, library-distance, and solver-specific confidence data.
  6. Held-out coupon and independently fabricated device references satisfy the accuracy,
     convergence, and performance gates in this document.
  7. Process libraries are reproducible, content-addressed, portable to remote runs, and
     retain the generator version and complete qualification evidence.

## Capability matrix

### Solver support

| Problem | Raw | Fixed trace/flux | Self-consistent | Remaining work |
|:--|:--:|:--:|:--:|:--|
| 2D electrostatic | Implemented | Implemented | Implemented | Production validation |
| 3D electrostatic | Implemented | Implemented | Implemented | Practical-device validation |
| 2D boundary mode | Implemented | Implemented | Not planned initially | Validate quasi-TEM scope |
| 3D uniform driven | Implemented | Implemented | Implemented | High-order/AMR validation |
| 3D adaptive PROM driven | Implemented | Implemented | Missing | Corrected reduced model |
| 3D linear undamped eigenmode | Implemented | Implemented | Implemented | High-order/AMR validation |
| 3D damped or frequency-dependent eigenmode | Implemented | Implemented | Missing | Nonlinear corrected EVP |

Standard GridFunction and ParaView output is currently the raw thin-metal field. Export
of the self-consistent corrected field is missing.

### Metal boundary support

| Boundary law | Runtime matching | Automatic coupon generation | Physical qualification |
|:--|:--:|:--:|:--:|
| PEC | Implemented | Implemented | Partial |
| Conductivity | Implemented | Missing | Missing |
| Impedance | Implemented | Missing | Missing |
| Rational impedance | Implemented | Missing | Missing |

Non-PEC libraries match numerical boundary parameters exactly. It remains to establish
whether a quasi-electrostatic coupon is independent of frequency and surface impedance.
If not, process libraries must index or interpolate response matrices in frequency and
boundary-law parameters.

### Local geometry support

| Geometry | Classification | Coupon generation | Qualification |
|:--|:--:|:--:|:--:|
| Isolated straight edge | Implemented | Implemented | Qualified in CPW studies |
| Parallel pair | Implemented | Implemented | Qualified over selected separations |
| Three or more parallel edges | Implemented | Implemented | Incomplete |
| Sharp 90 degree convex/concave corner | Implemented | Implemented | Partial |
| Rounded 90 degree convex/concave corner | Implemented | Implemented | Partial |
| Arbitrary corner angle | Implemented | Implemented | Partial |
| Smooth curved edge | Piecewise linear | Reuses edge/corner models | Missing |
| Manifold trace cap or T/cross mask | Corners, bends, and spatial clusters | Exact manifold mask | Incomplete |
| Open-graph endpoint | Implemented matching | Requires an exact closed mask | Incomplete |
| Nonmanifold junction | Implemented matching | Requires process regularization | Unsupported automatically |
| Nonparallel local cluster | Exact signature | Implemented | Incomplete |
| Cross-layer local cluster | Exact multi-slot signature | Implemented | Incomplete |
| Exact plan-view conductor mask | Implemented | Classified manifold mask | Incomplete |

Multiple independent fabrication layers and per-interface diagnostics are implemented.
Spatial coupon generation currently requires parallel or antiparallel fabrication planes.

## Current implementation

The following infrastructure is present:

  - Automatic two- and three-dimensional metal-edge extraction and conductor ownership.
  - Isolated, paired, parallel-cluster, corner, endpoint, junction, and exact spatial
    neighborhood classification.
  - Exact clipped plan-view masks for ambiguous spatial neighborhoods, canonicalized
    independently of surface triangulation.
  - Coupon-plan aggregation by canonical plan-view boundary, so equivalent source-facet
    triangulations retain one generated coupon and accumulate their occurrence count.
  - Classified physical and artificial continuation boundaries for exact endpoint,
    junction, and spatial masks. Fabricated coupons loft the complete manifold footprint,
    taper and round physical boundaries, and apply overetch only along those boundaries.
    Open and nonmanifold masks fail closed.
  - Circular chains in classified masks are reconstructed as smooth OCC arcs. Spatial
    trench geometry uses continuous offset-mask shells, automatic preparation requires
    at least quadratic geometry, and nonpositive curved-element Jacobians fail closed.
  - Spatial qualification separately gates FEM-order and local h-convergence. The h-study
    refines both thin and fabricated meshes at the final FEM order and reuses the finest
    completed p-study response.
  - Separation interpolation and held-out-qualification-gated rounded-corner radius
    interpolation, without extrapolation. Adjacent corner models alone do not authorize
    a runtime interpolation span.
  - Geometry-only preflight and a metadata-only version-3 process-library seed.
  - PEC coupon planning, generation, probe convergence, held-out qualification,
    content-addressed caching, partial merging, and explicit failure manifests.
  - Separate raw, fixed-trace, fixed-flux, and supported self-consistent output.
  - Correction at every AMR iteration, with cached geometry and distributed compact
    response traces.
  - Per-interface matched fractions and confidence diagnostics which prevent aggregation
    across large and small interfaces from hiding a local failure.

## Ordered implementation plan

### P0: Fabricated spatial geometry

**Goal:** A generated spatial coupon must represent the fabrication of every physical
plan-view boundary in its exact mask.

**Implementation status:** Complete for classified manifold masks. Spatial-cluster
preflight exports exact facets and classified boundaries; generation applies one
explicit regularization policy and rejects open or nonmanifold masks. A planar trace
cap or T/cross branch is represented by its manifold perimeter and therefore appears as
corners, smooth chains, and possibly an exact spatial cluster, not as a degree-one or
degree-three perimeter vertex. Physical qualification of these coupon families remains
part of P1.

  1. Represent mask boundaries as classified physical edges rather than vertical clipping
     planes.
  2. Apply process taper and top/bottom rounding to every exterior mask boundary.
  3. Derive trace caps from the exact device mask rather than inventing a closure from an
     open edge arm.
  4. Apply the same fabrication transform to every branch and reentrant boundary of a
     manifold T/cross mask.
  5. Reject nonmanifold masks, unresolved features, and incompatible cross-layer masks
     before meshing.
  6. Include the regularization policy in the coupon cache key and process-library model.

**Exit gate:** Trace-cap, T/cross-mask, and overlapping nonparallel-cluster coupons
generate thin and fabricated meshes without invented boundaries. Probe and held-out
qualification either pass or produce a stable, signature-specific failure. A graph
endpoint must carry an exact closed mask. Masks whose components touch only at a point
remain unsupported until a process-specific nonmanifold regularization policy is
defined.

### P1: Complete local coupon families

**Implementation status:** Geometry generation is available for the listed exact
manifold-mask families and arbitrary convex/concave corner angles, but physical
qualification is incomplete. Arbitrary-angle corner candidates use exact wedge and
circular-fillet geometry, quadratic mesh geometry, and a positive-Jacobian gate. Only
the recorded 90-degree rounded families have full-matrix and held-out qualification.
Rounded 45-degree convex, 135-degree convex, and 135-degree concave candidates pass
the p2-to-p3 and fixed-p3 h-refinement probe gates, but their full 72-trace and
held-out studies remain open. The rounded 45-degree concave candidate narrowly fails
the domain-defect p2-to-p3 gate (`5.09%` versus `5%`), despite a `0.05%` fabricated
domain change and sub-`1%` fabricated interface changes. Its defect is `0.666%` of
the fabricated matrix and changes by `0.034%` of that full response, but the
correction-specific gate remains fail-closed pending a p3-to-p4 study. Corner and
spatial candidates must pass both p-order and independent fixed-p h-refinement gates.
A rounded endpoint has valid quadratic geometry, yet its p=3 response from
`lc_fine=0.025 um` to `0.0125 um` still changes by `13.14%` for MA, `2.39%` for MS
(`16.01%` worst probe energy), and `8.23%` for SA. Automatic qualification records
this as an h-convergence failure instead of accepting the coupon on p-convergence
alone.

  1. Complete the 45-degree concave p3-to-p4 and h-refinement checks, then run
     full-matrix and held-out qualification for the passing acute and obtuse families.
  2. Run and qualify radius interpolation studies for rounded corners and smooth bends.
     The qualifier, version-3 metadata, fail-closed runtime selection, and combiner path
     are implemented. The historical 0.25/0.75-to-0.5 micron rounded-corner study fails
     the current strict defect-energy gate and is not a qualified span.
  3. Qualify manifold trace caps and T/cross masks over cap radius, branch angle, and
     branch width using corner, curved-edge, and exact spatial-cluster models.
  4. Qualify parallel clusters with three or more edges over separation and conductor
     ownership.
  5. Qualify nonparallel and cross-layer clusters over angle, offset, and layer spacing.
  6. Replace piecewise-linear curved-edge classification with curvature-aware segmentation
     and a discretization-independent signature.

**Exit gate:** The held-out local suite below passes at two matching radii and at least
two fabrication processes without geometry-specific tuning.

### P2: Non-PEC process response

  1. Generate conductivity, impedance, and rational-impedance coupon pairs.
  2. Sweep frequency and boundary-law parameters against fabrication-resolved Maxwell
     coupons.
  3. Decide whether matrices are universal, indexed, or interpolated in those parameters.
  4. Encode the result in process-library compatibility and cache keys.

**Exit gate:** Runtime matching cannot select a coupon outside its qualified frequency
and boundary-law range, and held-out Maxwell coupon tests pass the accuracy gates.

### P3: Solver completeness

  1. Add self-consistent correction to adaptive PROM driven sweeps.
  2. Add the corrected nonlinear operator required by damped and frequency-dependent
     eigenmodes.
  3. Add optional corrected-field and corrected-flux GridFunction export.
  4. Determine whether raw-field AMR resolves every response contour. If not, add a
     correction-trace indicator without changing raw output semantics.

**Exit gate:** Corrected results agree between direct and reduced driven solves, mode
pairing is robust near crossings, corrected fields reproduce reported energies, and AMR
convergence is independent of the initial mesh.

### P4: Practical 3D validation

Run independent thin/raw, thin/corrected, and fabrication-resolved references for:

  1. Extruded CPW, for electrostatic, driven, and eigenmode baselines.
  2. Sharp and rounded convex and concave islands or apertures.
  3. Straight and curved resonators.
  4. A single transmon with readout and resonator.
  5. A capacitive or inductive coupler.
  6. Flip-chip transmon and coupler layouts with separate L1/L2 SA/MS/MA.
  7. At least one finite-impedance Maxwell case after P2.

Use order refinement and AMR until the fabricated reference and corrected thin result
independently satisfy the convergence gate. Report each physical interface separately;
global L1+L2 or SA+MS+MA aggregates are not acceptance metrics.

### P5: Confidence, performance, and release workflow

  1. Fit confidence limits using held-out passes and failures rather than hand-selected
     thresholds.
  2. Stage process libraries and matrix files automatically for remote execution.
  3. Benchmark one- and multi-node transmon, coupler, and flip-chip jobs.
  4. Add solver-level regression cases for preflight, correction output, and raw-output
     invariance.
  5. Freeze the process-library schema and document migration/versioning.

## Validation and acceptance gates

These are release targets. A study may use stricter gates.

### Coupon convergence

  - The final two FEM orders change every material fabricated-response matrix by at most
    5% in relative Frobenius norm.
  - The worst quadratic-energy change over the probe space is at most 10%.
  - The fabricated-minus-thin domain-response defect changes by at most 5%.
  - Held-out traces not used to construct the basis reproduce each material interface
    energy within 5% when that channel carries at least 1% of the coupon response.
  - Weaker channels report both relative and response-normalized absolute error and must
    not dominate a confidence decision through an unstable relative denominator.

### Device accuracy

  - Fabricated-reference participation changes by at most 2% over the final two
    order/refinement levels.
  - Corrected thin participation changes by at most 2% over the final two levels.
  - Corrected participation is within 5% of the fabricated reference for each material
    interface carrying at least 1% of the device's maximum interface participation.
  - No accepted confidence result may exceed 10% error on a material interface in the
    held-out device suite.
  - Raw thin participation is recorded at every level but is not required to converge.

### Coverage and diagnostics

  - Validation runs use `UnmatchedPolicy: "Error"` and achieve complete geometric
    coverage.
  - `kR`, boundary-law verification, loop residual, local trace closure, curvature,
    library distance, and corrected-mode overlap meet their documented confidence limits.
  - A failed confidence gate must not suppress raw or corrected output, but the output
    must be clearly marked unqualified.

### Performance

For practical multi-node cases, relative to the same raw solve:

  - Postprocessing-only correction targets at most 10% additional wall time and memory.
  - Self-consistent correction targets at most 25% additional wall time and 15% peak
    memory.
  - Preflight and coupon preparation are accounted separately from the device solve.
  - No rank may retain a replicated device-scale surface trace or mesh representation.

## Validation artifacts

Every completed validation case must preserve:

  - thin and fabrication-resolved mesh-generation inputs;
  - resolved Palace configurations and Palace commit;
  - process-library JSON, matrices, cache key, and qualification manifest;
  - order/AMR convergence tables for raw, fixed-trace, fixed-flux, and self-consistent
    values;
  - per-interface errors and confidence diagnostics;
  - elapsed-time and peak-memory reports; and
  - a script which reproduces comparison tables without hand-edited values.

## Explicitly out of scope

This roadmap does not require a surface BEM, a fabrication-resolved full-device mesh, or
a microscopic model of Josephson-junction barriers. Junction leads and nearby planar
metal edges are in scope; the junction barrier remains represented by Palace's lumped
element model unless a separate feature changes that contract.

## Maintenance

Update this document in the same commit whenever solver support, coupon-generation
coverage, confidence gates, or validation status changes. A code path is not marked
qualified solely because it compiles or passes a synthetic unit test.
