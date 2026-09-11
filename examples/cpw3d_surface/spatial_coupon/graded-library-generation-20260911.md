<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Corrected graded-library generation campaign — 2026-09-11

Campaign: `/data/home/simlap/transmon_library_graded_v2_20260911`.
Retained reference: `/data/home/simlap/transmon_library_speed_20260909`.

This supersedes the blocked preparation described in
[the previous campaign notes](graded-library-campaign-20260910.md). It is an
isolated candidate build, not automatic production-library promotion.

## Corrections

- Box corners are explicit at every trace level. No source triangle cuts across
  a corner of the actual matching boundary.
- Cap triangulation preserves every collinear boundary node. Source functions
  are continuous across the complete closed matching surface.
- All retained nodal degrees of freedom survive. New corner nodes are added;
  `basis-contract.json` records old-to-new indices, bounded coordinate snapping,
  source hashes, conductor states, and the deliberate source-definition change.
- Per-kind conductor `TerminalAttributes`, original p-order/linear settings,
  materials, interface settings, and process dimensions are retained.
- The legacy finite etch footprint is preserved explicitly where required.
- Every Palace-consumed mesh is **MSH 2.2 binary**. Relabeling has serialized
  geometry/connectivity checks and remapped partition certificates.

| Spatial case | Original sources per kind | Corrected sources per kind | Order |
|---|---:|---:|---:|
| Two-edge `3f8992613e95` | 115 | 125 | 4 |
| Three-edge `419576fdab24` | 130 | 135 | 4 |
| Four-edge `9d2cb9bbb3fe` | 80 | 80 | 4 |
| Three-edge `7f03270dca8e` | 180 | 185 | 5 |
| Ten-edge `6791f1c84123` | 161 | 161 | 5 |
| Two-edge `8dd4bc70f183` | 68 | 78 | 5 |

Including the eight unchanged non-spatial thin/fabricated cases, the complete
campaign has **20 cases and 2,198 sources**, compared with 2,138 previously.
The three geometry-discovery descriptors remain present but are not evaluated.
Strict metadata preflight passed: **3,120 exact matches, zero interpolated or
missing requirements**.

## Mesh policy and accuracy scope

The selected candidate policy uses level constraints plus local grading around
all declared trace edges, not full cap-diagonal CAD constraints. The latter
produced much worse element conditioning. Physical surface minima are 2 nm
(thin) and 0.5 nm (fabricated); the volume minimum is independently 2 nm. Trace
sizing is 10 nm, with rapid 3D growth away from both process and trace features.
Surface meshing uses Delaunay (algorithm 5); volume meshing uses HXT without the
additional long-running Netgen optimization.

The corrected three-edge nine-source thin diagnostic showed:

- full CAD trace constraints: 3,489,692 tets, maximum kappa 1743;
- level constraints with local trace grading: 3,357,707 tets, maximum kappa 53.9;
- p4 pilot time 605.4 s versus 311.3 s; p5 1880.3 s versus 815.6 s;
- p4-to-p5 matrix-norm change on the level mesh about 0.0191%; the sharpest
  sampled direction still changed by about 0.914%.

These are diagnostics, not a full-basis or continuum qualification. Full
matrices, fabricated MA/MS/SA response, slot ownership, domain defects, and
corrected-device observables remain validation outputs of the larger campaign.
Thin edge-inclusive raw SPR is not used as an accuracy gate.

## Pipeline integrity gates

`prepare_graded_library_campaign.py` freezes the model/config/source contracts
and prepares independent mesh and response PBS jobs.
`run_graded_library_mesh.py` verifies physical geometry, complete expected
attributes, tetrahedral validity, serialized format, and an exact recipe binding
all relevant input/tool hashes and grading settings. Unbound pilot-mesh reuse
is rejected.

Archive admission must use the **processed** H1 space. Palace may crack internal
PEC boundaries and refine crack-adjacent elements; raw mesh topology is not a
safe count. `palace --mesh-statistics CONFIG` runs the normal preprocessing path
and reports H1 size without assembling a PDE operator. A regression demonstrates
14 raw p2 DOFs becoming 20 after cracking, matching the real solve. A 42,822,378-
DOF production-scale check likewise matched the frozen solver across 1 and 192
ranks. Cracking is not disabled.

The mesh-statistics helper is separate from the frozen response executable.
Response jobs still use binary SHA256
`e1d14c903c6b0b41ce65221a57c569ab18c9c804e9b0ac0a6e1dfcffd8cf6cbd`.
Probe costs are recorded separately from worker/reducer timing.

`archive_storage.py` serializes shared-space reservations and accounts for
outstanding writers and sealed archives. Source/mesh/config/executable hashes,
complete archived source/rank sets, and convergence are checked.
`run_graded_library_case.py` requires mandatory energy columns, exact interface/
edge/radius groups, complete new-basis pairs, finite values, and energy-matrix
sanity before copying results. Unchanged controls must also reproduce their old
matrices. Corner core-energy conversion explicitly writes LF endings.

The manifest is published last, only after every numerical case completes, as
`GeneratedCandidateNotAccuracyQualified`. This label is intentional.

## Execution

- Mesh jobs: one c8g.48xlarge each, bounded single-threaded meshing.
- Response jobs: original per-case r8g.48xlarge node/rank allocations.
- Project `DS-EM-FEM`; subnet `subnet-0c98d793bbcebb39a`.
- `submit_graded_library_campaign.py` enforces fewer than or equal to 40 user
  jobs and binds response jobs to successful mesh jobs with `afterok`.
- The eight unchanged control cases completed successfully before the full
  spatial launch. Spatial mesh jobs are 41679–41690; their response jobs are
  41691–41702.

Inspect `submissions.json`, per-case `status.json`, per-mesh `mesh-state.json`,
`archive-storage.json`, and `summary.json` under the campaign directory.
Run `summarize_graded_library_campaign.py ROOT` to refresh the summary.

The historical 9 h 19 min generation is not a clean mesh-only speedup baseline:
the corrected basis adds sources and changes previously invalid cap functions.
Report source counts, worker/reducer times, meshing time, setup and resource use
separately. Existing reference libraries and meshes remain unchanged.
