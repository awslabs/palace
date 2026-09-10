<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Local all-tet mesher debugging (2026-09-09)

No new SOCA jobs were submitted. Only small existing input/reference files were fetched.
Gmsh version on the local machine: `4.15.0-git`; SOCA had used `4.13.1`.

## Geometry corrections

- Thin exact-mask geometry now uses the full mask, not the intersection of finite edge
  strips with the mask. The three-edge thin coupon's metal area changed from the
  incorrect 63.63636 um^2 to 122.54545 um^2, matching the retained prism reference.
- Etch collars subtract the union of **all** retained conductor masks. This prevents one
  conductor's collar etching underneath another. The ten-edge fabricated coupon's
  substrate-volume deficit of 0.5 um^3 disappeared.
- Three-edge thin/fabricated CAD areas and volumes match their references to roughly
  3e-11 relative. Ten-edge total physical family areas/volumes also match, but per-slot
  surface assignment is still incorrect: assigning a whole CAD face from one sample
  point cannot represent the reference's multiple slots on that face.

The experiment requires `--reference-measures CSV` before full meshing. Missing physical
attributes or area/volume differences >1e-5 fail **before** mesh generation. A geometry-only
mode records CAD measures without generating elements. This is a necessary check, not a
complete proof of geometry equivalence or field accuracy. The ten-edge case currently
fails this strict gate and is not meshed for physics comparison.

## Local results (not qualified response libraries)

Outputs: `/tmp/local-tet-debug`.

| Trial | Elements | Runtime | Peak process-tree RSS | Max kappa |
|---|---:|---:|---:|---:|
| Fully isotropic, 2 nm target, thin | No mesh | stopped at 180 s | 0.75 GiB | — |
| Fully isotropic, 20 nm scout, thin | 277,724 | 37.2 s | 1.83 GiB | not audited |
| 2 nm face/volume target, 100 nm edge-tangent cap, growth 1, thin | 140,416 | 24.5 s | 1.31 GiB | 13.67 |
| Same, fabricated | 138,761 | 29.0 s | 2.54 GiB | 18.84 |
| Same, growth 0.5, thin | 490,799 | 100.9 s | 3.22 GiB | 10.48 |

The isotropic 2 nm run was active in Gmsh's initial 2D surface triangulation at timeout,
not out of memory. The 20 nm scout only tests the meshing pipeline; it is not an accuracy
candidate.

The successful 2 nm-target trials use **dimension-dependent boundary sizing**:
`TET_EDGE_TANGENT_SIZE=0.1` caps horizontal feature-curve spacing at 100 nm, with an
endpoint ramp down to the minimum target. Face/volume sizing still uses exact 3D
point-to-segment distance. This is NOT uniformly isotropic 2 nm meshing near edges, nor
an enforced 2 nm first-normal layer. The actual near-edge discretization and p5 energies
must be qualified before any speedup/accuracy claim. Metadata explicitly records
`EdgeCurveTangentialSize` and `NormalLayerEnforced: false`.

The successful meshes are all tetrahedral, geometry order 2. Their measured material
volumes and boundary areas agree with the retained three-edge references to about
3e-11 relative. They have no nonpositive element-center Jacobians. Post-mesh max center
versus corner volume discrepancies are 1.5e-7 (thin) and 3e-7 (fabricated); aggregate
geometry measures nevertheless agree closely. A center-only audit is not an all-point
Jacobian certificate; Gmsh also checks positive minSJ during generation.

A p5 held-out solve on the small thin candidate assembled a 3,083,795-DOF H1 space, but
both six-rank and two-rank local attempts hit the conservative 8 GiB process-tree RSS
cap during setup. They were terminated in about 25 seconds, and no participation or
accuracy result is available. No larger-memory retry or remote solve was launched.

## Reproduce locally

```sh
TET_EDGE_TANGENT_SIZE=0.1 python3 run_bounded_mesher.py \
  --seconds 120 --memory-gib 8 --log /tmp/trial.log -- \
  julia --startup-file=no --project=/path/to/palace/test/examples \
  mesh_graded_tet_experiment.jl /path/to/signature-dir thin /tmp/trial.msh \
  1.0 0.002 0.5 --reference-measures /path/to/reference-measures.csv
```

The reference CSV has `dimension,attribute,measure` columns, with dimension 2 for areas
and 3 for volumes, using the same original mesh units. Use `--geometry-only` to stop at
CAD validation. Existing mesh outputs/logs are not overwritten.

Regression tests in `test_graded_tet_geometry.jl` passed 11 assertions covering exact-mask
membership beyond the finite edge strips, conductor-union etch exclusion, the named
ARM-compatible C callback, exact segment-interior sizing, and dimension-dependent sizing.
The bounded runner terminates the complete local process group on a time/memory limit.

Remaining work: correct multi-slot surface partitioning, measure actual near-wall
resolution, and run p5 selected-source/held-out energy comparisons on a host with enough
memory. Do not substitute these trial meshes or matrices into the reference library yet.
