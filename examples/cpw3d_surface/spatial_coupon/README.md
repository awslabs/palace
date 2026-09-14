<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Spatial coupon meshers

## Exact blockwise response matrices

Large high-order coupons can exceed node memory because the ordinary response-matrix
path retains every potential and recovered-flux field. The default fast path streams the
complete source basis through one Palace process, retaining only one source field while
reusing the mesh, assembled operator, AMG hierarchy, and preconditioner:

```sh
python3 prepare_blockwise_response.py spatial_fabricated.json \
  --output /tmp/fabricated-stream \
  --execution-mode streaming \
  --reduction-block-size 3 \
  --ranks 6 \
  --local

PALACE=/path/to/palace-arm64.bin /tmp/fabricated-stream/run-blocks.sh
PALACE=/path/to/palace-arm64.bin /tmp/fabricated-stream/run-reducer.sh
```

For a generated coupon directory containing both thin and fabricated configs, prepare the
complete exact two-command workflow with:

```sh
python3 prepare_fast_coupon_response.py /path/to/coupon-output \
  --ranks 6 --reduction-block-size 6 --local
PALACE=/path/to/palace-arm64.bin \
  /path/to/coupon-output/fast-response/run-all.sh
```

The second command streams and reduces both cases, then installs the matrices at the
paths already referenced by `process-library.json`.

`PALACE_RESPONSE_ARCHIVE_ONLY=1` skips every ordinary per-source table, energy
measurement, estimator update, and within-run response matrix. It archives the potential
and recovered flux and immediately releases the volume field. This is exact: a 72-source
local benchmark produced bitwise-identical archives to source blocks while reducing the
worker phase from 78.6 s to 5.4 s (14.4 times). Increasing the reducer block from 6 to 24
reduced its time from 25.5 s to 9.8 s; choose the largest block that fits memory.

Source blocks remain available for machines that cannot hold one coupon operator or for
independent cluster scheduling:

```sh
python3 prepare_blockwise_response.py spatial_fabricated.json \
  --output /tmp/fabricated-blockwise \
  --execution-mode blocks \
  --sources-per-block 3 \
  --reduction-block-size 3 \
  --ranks 192
```

The command preserves every original `PrescribedPotential.Index`, writes source configs,
a full-source reducer config, a manifest, and local- or PBS-friendly shell commands.
Archive workers use:

```sh
PALACE_RESPONSE_ARCHIVE_DIR=/shared/archive mpirun ... palace block-0000.json
```

Each MPI rank writes deterministic `V` and recovered-`D` true-dof files with validated
binary headers. After every block succeeds, run the reducer on the same rank count and
mesh partition:

```sh
PALACE_RESPONSE_ARCHIVE_DIR=/shared/archive \
PALACE_RESPONSE_REDUCE_ONLY=1 \
PALACE_RESPONSE_BLOCK_SIZE=3 \
  mpirun ... palace reducer.json
```

The reducer loads at most two blocks, uses Palace's existing batched interface Gram
assembly, and computes tiled domain products. The resulting CSV matrices retain the full
nodal basis and are not a reduced-order approximation. All archive and reducer runs must
use the same MPI rank count, Palace build, mesh, and partitioning. Archive files can be
deleted after reduction.

`mesh_spatial_coupon_swept.jl` is the explicit swept-mesh implementation for
fabrication-resolved spatial response coupons. It is separate from
`mesh_spatial_coupon.jl`, whose field-based anisotropic tetrahedral mode remains
experimental. See [graded mesh qualification](graded-mesh-qualification.md) for
configurable process geometry, generality scouts, supported/unsupported cases,
and the provisional engineering versus strict accuracy profiles. The focused
[matched volume-growth benchmark](matched-volume-growth-benchmark.md) and its
[measured mesh results](mesh-growth-results-20260910.md) track DOFs, conditioning,
source timings, memory, and accuracy without changing the integration method.

`frozen_volume_study.jl` generates independent volume candidates from one fixed
surface mesh. For multislot coupons, `relabel_frozen_interface_mesh.jl` creates a
derived **Gmsh 2.2** mesh after volume generation. It applies exact exclusive
ownership, retains the frozen input unchanged, and verifies in-memory IDs plus
serialized geometry/connectivity, quality, grouped physical-family areas,
expected attributes, and SHA-bound partition certificates. MSH 2.2 is
intentional: it is the format consumed by the Palace/MFEM reader used by these
campaigns.

The [graded-library campaign notes](graded-library-campaign-20260910.md) record
full-scale source-contract checks and remaining qualification gates. Retained
artificial etch collars can be exported with `export_etch_footprint.jl` and
supplied explicitly using `--etch-boundary`; this preserves their actual planar
footprint rather than changing process geometry during a mesh comparison.
`audit_trace_continuity.py` detects side/cap T-junction jumps, and
`repair_triangle_trace_caps.py` writes separate corrected traces while recording
that the source functions changed. Such repairs require matched controls and
must not be presented as byte-identical source migration.

The [corrected full-generation workflow](graded-library-generation-20260911.md)
completes missing box corners without dropping retained basis degrees of freedom.
It freezes the source/config contracts, gates mesh jobs before response jobs,
reserves shared archive space, and validates complete matrices before publication.
`palace --mesh-statistics CONFIG` reports H1 size after actual mesh preprocessing
(including internal-boundary cracking), so archive budgets do not rely on raw
mesh counts.

## Trace-sizing safety policy (experimental tetrahedral mesher)

Trace **source data**, CAD break-line constraints, and scalar mesh sizing are separate
contracts. None of the controls below removes or rewrites any source `TraceMesh` triangle,
node, basis DOF, or prescribed value. They do not change `p`, tolerances, materials, or
physical geometry. Keeping a CSV byte-identical does **not** imply its piecewise function
is represented exactly by the resulting FEM trace space.

`mesh_graded_tet_experiment.jl` defaults to `TET_TRACE_SIZE=0` and
`TET_TRACE_SIZE_SCOPE=off`: the physical-interface distance sizing remains available
without extra trace-driven scalar refinement. A positive size is rejected unless both
scope and `TET_TRACE_CONSTRAINT_MODE` are explicitly selected. Lengths are in micrometres.

| `TET_TRACE_SIZE_SCOPE` | Where the extra scalar trace size applies |
| --- | --- |
| `off` | Nowhere; requires size zero and relative size zero. |
| `matching` | Only dim=2 queries whose CAD entity tags are in the matching-surface set. |
| `matching-and-volume` | Those matching surfaces plus dim=3 queries, explicitly opting into volume grading. |
| `legacy-global` | Historical replay only: all dimensions/entities, including physical surfaces, from **every** triangle edge regardless of CAD mode. |

In the two scoped modes, dim=1 curves (including matching curves), dim=0 points,
unknown-dimension queries, and physical-interface surfaces receive no trace sizing.
The existing physical-interface field and any separately selected non-trace controls
remain in effect. Shared conforming boundaries may still affect the generated mesh;
unchanged callback values are not a promise of identical mesh connectivity or counts.
No anisotropy is introduced by these scalar policies.

The explicit shared CAD/sizing mode is `all` (all triangle edges on box surfaces), `sides`
(side-face edges only), or `levels` (interior horizontal box rings only). These are
meshing subsets, never source-basis reductions. `levels` does not retain every prescribed
break line, and uses the maximum trace target on each ring. `none` creates no CAD trace
lines and is permitted only with scalar sizing off; the supplied trace geometry is
still validated, without creating/fragmenting CAD entities for it. Local triangle
validity is checked before selection in every mode; complete box coverage and source
continuity remain separate gates.

A **physical-only diagnostic** can use the following environment with the usual
mesher arguments, including `--matching-trace /path/to/unchanged/basis-0001.csv`:

```sh
export TET_TRACE_SIZE=0 TET_TRACE_SIZE_SCOPE=off TET_TRACE_RELATIVE_SIZE=0
export TET_TRACE_CONSTRAINT_MODE=none
```

For an explicitly unqualified matching-only scout, select e.g.
`TET_TRACE_SIZE=0.01 TET_TRACE_SIZE_SCOPE=matching TET_TRACE_CONSTRAINT_MODE=levels`.
`TET_TRACE_SURFACE_GROWTH` defaults to 4 for scoped policies; volume growth uses the
existing trace volume profile (default 0.5). `TET_TRACE_RELATIVE_SIZE=0` requests the
fixed maximum target. Positive relative sizing in `all`/`sides` caps each edge target
by that edge's opposite altitude times the relative factor. This is only a geometric
heuristic: long edges of skinny triangles can still demand tiny isotropic sizes.
Per-edge altitude alone neither resolves the accuracy question nor cures the cost.

To replay the old global behavior, explicitly select `legacy-global`, a positive size,
and the historical CAD mode. It intentionally ignores that mode when choosing scalar
segments, assigns each triangle's **minimum** altitude to all three edges for relative
sizing, and uses growth 0.5 in all dimensions (independent of volume trace profiles).
A different `TET_TRACE_SURFACE_GROWTH` is rejected for this replay policy. This preserves
the old pathological far-cap/diagonal refinement; it is not an accuracy recommendation.

New `prepare_graded_library_campaign.py` invocations require `--trace-size-scope`
(`--trace-size` is also needed for a positive policy), in addition to `--trace-mode`.
Optional `--trace-relative-size` and `--trace-surface-growth` are recorded. An old
manifest with positive `TraceSize` but no `TraceSizeScope` fails **before cache reuse or
jobs**; it is not silently migrated. Use a fresh, explicitly chosen recipe, never edit
retained campaign artifacts to bypass this failure. All trace controls and helper code
are recipe-bound; the runner excludes inherited `TET_*` overrides. The separate
`prepare_graded_library_meshes.py` retains its no-scalar-sizing baseline explicitly in
its plan/runner, with an explicit CAD-mode choice. Choosing a diagnostic policy does not
qualify it for production.

Outstanding gates, before accepting any diagnostic/scoped mesh for a library:
complete immutable source contract and box coverage/continuity; actual FEM boundary
projection/interpolation error for the full trace basis; achieved interface-normal
resolution and acceptable Jacobians; processed H1 DOFs; and complete matched domain
and per-interface response/energy matrices at unchanged order, tolerance, materials,
and geometry. Mesh size, unchanged source bytes, and local altitude are not substitutes.
No new PDE accuracy result is claimed here.

## Straight-edge swept milestone

The straight-edge branch accepts one continuing edge without a plan-view mask.
It:

- constructs explicit graded coordinates in the edge-normal cross-section;
- constructs longitudinal stations independently from `--lc-tangent`;
- triangulates the conforming substrate/vacuum cross-section and sweeps it into
  quadratic prism elements;
- represents fabricated sidewalls and top/trench fillets in the cross-section;
- preserves the Palace domain, matching-surface, SA, thin-metal, MS, and MA
  physical attributes;
- rejects non-prism, nonmanifold, over-budget, or nonpositive-Jacobian meshes;
- reopens the serialized MSH file and records independently checked evidence in
  `<mesh>.metadata.json`.

Generate a test coupon from the repository root with an instantiated Julia
Gmsh environment:

```sh
julia --project=examples \
  examples/cpw3d_surface/spatial_coupon/mesh_spatial_coupon_swept.jl \
  examples/cpw3d_surface/spatial_coupon/testdata/straight-edge-signature.csv \
  fabricated /tmp/straight-edge-fabricated.msh \
  --radius 1.0 \
  --metal-thickness 0.2 \
  --overetch 0.1 \
  --top-radius 0.02 \
  --bottom-radius 0.02 \
  --lc-normal 0.1 \
  --lc-tangent 1.0 \
  --lc-far 0.4 \
  --mesh-order 2
```

`--lc-fine` is retained as an alias for `--lc-normal` so the existing coupon
orchestrator can invoke either mesher. Select this implementation explicitly in
that orchestrator with `--spatial-mesher swept`. Mesher identity and the normal
growth ratio are included in the cache fingerprint.

Run the executable mesh checks with:

```sh
python3 \
  examples/cpw3d_surface/spatial_coupon/test_mesh_spatial_coupon_swept.py
```

The tests generate thin and fabricated meshes and check normal/tangent
independence, quadratic prism topology, physical attributes, serialized counts,
positive scaled Jacobians, face conformity, and budget failure cleanup.

## Masked edge-cluster transition mesher

A second branch accepts masked `SpatialEdgeCluster` inputs with any positive
edge count. It uses the classified mask boundary as geometry evidence and
constructs geometrically spaced offset rows around every active finite segment.
OCC fragments all row intersections in one planar mesh, so orthogonal, parallel,
nearby-edge, and multi-edge junction transition regions are conforming. The
planar triangles are then swept through independently graded
fabrication-normal coordinates into quadratic prisms.

The masked mesher is intentionally limited to 90 degree sidewalls and zero
rounding, which is the process in the checked-in single-transmon preflight seed.
Nonzero taper or rounding fails closed. All ten consolidated transmon spatial
families—seven 2-edge masks and one each with 3, 4, and 6 edges—generate both
thin and fabricated meshes with positive scaled Jacobians and the expected
multi-slot and multi-conductor attributes. At 20 nm normal and 300 nm tangent
spacing, the remaining fabricated transmon families contain:

| Edge count | Quadratic nodes | Prisms |
| ---------: | --------------: | -----: |
| 3          |       1,030,114 | 248,064 |
| 4          |       1,365,159 | 331,034 |
| 6          |         967,387 | 232,692 |

The six-edge, two-conductor/two-slot mesh passes a Palace held-out solve at p=1
and exposes distinct conductor/slot MS and MA attributes.

The representative orthogonal two-edge coupon scales as follows with fixed
`lc_tangent = lc_far = 0.3 um`, a `0.2 um` collar, growth ratio 1.4, and
quadratic geometry:

| First normal spacing | Thin nodes | Fabricated nodes | Fabricated prisms |
| -------------------: | ---------: | ---------------: | -----------------:|
| 20 nm                |    734,461 |          959,804 |            230,892 |
| 10 nm                |    918,213 |        1,259,570 |            304,544 |
| 5 nm                 |  1,120,405 |        1,591,768 |            386,300 |
| 2 nm                 |  1,371,889 |        2,156,767 |            525,360 |
| 1 nm                 |  1,643,151 |        2,672,949 |            652,676 |

Thus a 20-fold reduction in first spacing increases fabricated quadratic nodes
by 2.8 times, rather than by a quadratic or cubic factor. Mesh generation for
each level took less than 15 seconds on the development Apple M3 Pro.

The spatial probe traces now use a smooth cutoff over the complete thin/fabricated
metal-height band. This removes the artificial order-dependent Dirichlet jump at
the grounded-metal/matching-surface intersection. At fixed p=1, the representative
fabricated response changed from 2 nm to 1 nm by:

| Quantity      | Matrix norm change | Worst probe-energy change |
|:--------------| -----------------:| ------------------------:|
| Domain        |             0.001% |                     0.41% |
| Domain defect |             0.099% |                         - |
| MA            |             1.424% |                     2.46% |
| MS            |             0.068% |                     1.27% |
| SA            |             0.427% |                     0.54% |

At p=2, the 5 nm to 2 nm transition changed the fabricated domain matrix by
`0.0004%` and the fabricated MA/MS/SA matrices by
`1.203%/0.023%/0.440%`. Worst probe-energy changes were
`0.029%/2.085%/0.230%/0.577%` for domain/MA/MS/SA. Thus the 2 nm mesh passes the
existing fabricated-response gates for this representative coupon; the p=1
2 nm to 1 nm comparison independently confirms the final trend. Final library
qualification still requires the independent p-order study on the selected
mesh.

### Multiconductor trace lift

Multiconductor terminal probes originally failed to converge even though their
smooth contour blocks converged. The exact mask labeling had been intersected
with finite edge strips, so conductor continuation regions on the matching
surface were incorrectly left as free zero-potential trace nodes. A terminal
source then imposed one volt on the complete metal boundary against zero on
those mislabeled contact nodes, producing a mesh-dependent Dirichlet jump.

For exact masks, conductor labels now come from the complete conductor facets.
Each nonreference conductor source also carries a contact lift which is one on
its excluded matching-surface contact nodes and zero on every free trace knot.
The lift therefore changes no runtime contour coefficient, while making the
terminal/matching intersection compatible. Multiconductor libraries use
`TraceLiftVersion = 3`.

With the corrected lift, the six-edge family passes p3 to p4 at 2 nm:

| Quantity | Fabricated matrix change | Worst probe-energy change |
|:---------| ------------------------:| ------------------------:|
| Domain   |                    0.295% |                    0.866% |
| MA       |                    2.404% |                    8.467% |
| MS       |                    1.459% |                    2.140% |
| SA       |                    0.884% |                    2.396% |

Its 2 nm to 1 nm p3 transition also passes, with a `2.46%` domain-defect change
and fabricated domain/MA/MS/SA matrix changes of
`0.005%/0.611%/0.035%/0.256%`. The six-edge library model should therefore use
a 1 nm normal mesh and p4; the other nine spatial models use 2 nm and p4.

## Remaining transition milestones

General tapered/rounded masked fabrication and endpoint/junction families still
fail closed. They do not fall back to the field-based mesher. The checked-in
single-transmon spatial edge-cluster requirement set is covered.

`replay_device_trace.py --max-relative-error VALUE` provides the pass/fail gate
for the later replay step. It writes per-domain and per-interface relative
errors, the maximum error, the requested limit, and `Passed` to
`device-trace-replay.json`, and exits unsuccessfully when the limit is exceeded.
