<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Coupon mesh-quality study, 2026-09-09

Read-only audit of six existing coupon meshes. No mesh, response library, solver setting,
or running generation job was changed. SOCA green job **41282** completed successfully
in 165 seconds on one c8g.48xlarge node, using a single analysis process. The largest
mesh analysis took 57 seconds and about 6 GiB RSS. No high-order solution space or PDE
operator was assembled. Each case had a 300-second timeout.

## Reproduction and scope

- Audit source: `audit_coupon_mesh.cpp` in this directory.
- Analytic-prism test: `test_audit_coupon_mesh.py` (set `PALACE_MESH_AUDIT_BINARY`).
- Remote working directory: `/data/home/simlap/coupon_mesh_quality_20260909`.
- Full JSON results: `reports/{06,08,09}-{thin,fabricated}.json`.
- Source input paths/hashes and timing: `status.json`, `logs/*.time`.
- Actual link command: `build-command.txt`; uses the installed MFEM from the benchmark.

The three model identifiers are:

- 06: `spatialedgecluster_edgecount-3_419576fdab24`, p4.
- 08: `spatialedgecluster_edgecount-3_7f03270dca8e`, p5.
- 09: `spatialedgecluster_edgecount-10_6791f1c84123`, p5.

Kappa is exactly MFEM's element-center `GetElementJacobian` singular-value ratio,
including its perfect-reference-element normalization. Extrema reproduce the existing
Palace logs. Statistics are element-count percentiles, not volume-weighted percentiles.
Volume fractions use the center Jacobian; comparison with straight-prism corner volumes
agrees within 2.3e-12 relative and maximum XY drift along the sweep is zero. All center
Jacobians are positive. This is not an all-quadrature-point Jacobian certificate.

Coordinates are the original mesh coordinates, in micrometres. Distance statistics use
centroids: vertical distance to the process planes and XY distance to `Physical` segments
of the original `plan-view-boundary.csv`, excluding `Continuation` edges. The fabricated
plane offsets are the study process's 100 nm thickness and 50 nm overetch. These distances
are not exact distances to the finite MA/MS/SA surface patches.

## Element condition numbers

| Model/case | Elements | Median kappa | P90 | P99 | P99.9 | Maximum | Count >1000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 06 thin | 394,416 | 19.92 | 123.72 | 173.43 | 290.48 | 1419.30 | 12 |
| 06 fabricated | 620,832 | 42.23 | 147.06 | 173.43 | 240.49 | 1174.78 | 8 |
| 08 thin | 1,504,260 | 13.02 | 61.86 | 234.37 | 328.12 | 988.30 | 0 |
| 08 fabricated | 2,350,416 | 18.22 | 61.86 | 240.23 | 297.19 | 818.03 | 0 |
| 09 thin | 2,902,680 | 12.99 | 57.75 | 234.37 | 377.23 | 5606.45 | 260 |
| 09 fabricated | 4,393,024 | 16.10 | 57.74 | 212.28 | 305.61 | 4640.54 | 256 |

The very large extrema are rare, not characteristic of the entire mesh. Every element
with kappa above 1000 in these six meshes is more than 0.2 um from a process plane.
The condition numbers alone do not establish their contribution to PCG runtime.

## In-plane triangle shape

Thin/fabricated cases share their plan triangulation.

| Model | Unique plan triangles | Median plan kappa | Maximum plan kappa | Minimum angle | Minimum altitude |
|---|---:|---:|---:|---:|---:|
| 06 | 10,956 | 1.41 | 416.83 | 0.16565 deg | 0.3764 nm |
| 08 | 41,785 | 2.74 | 288.68 | 0.22918 deg | 0.5764 nm |
| 09 | 80,630 | 2.99 | 430.93 | 0.15394 deg | 0.1010 nm |

The raw `TrianglesWithEdgeBelow2nm` counter includes floating-point deviations from
nominal 2 nm spacing; it must not be read as a count of genuinely undersized features.
Minimum altitude, rather than shortest edge, exposes the clearest issue here.

For model 09, one worst plan triangle has vertices (um):

```
(5.498000000000, 0.400000000000)
(5.500000000000, 0.400000000000)
(5.537545787546, 0.398000000000)
```

Its shortest edge is 2 nm, but staggering produces a 0.101 nm altitude and 0.154-degree
angle. This is not a 0.1 nm fabrication feature. It is a skinny triangulation of a region
near a 90-degree corner. The thin mesh repeats that triangle through all 36 layers.
A worst prism lies at approximately (5.51185, 0.39933, 1.80978) um, where a 580.43 nm
vertical interval combines with this tiny in-plane altitude to give kappa 5606.45.
The fabricated counterpart has a 480.43 nm interval and kappa 4640.54.

## Refinement propagation

The plan mesh is repeated through 36 thin or 64 fabricated vertical intervals, excluding
metal volumes. The p5 meshes retain 2 nm first spacing, 100 nm near-edge tangential
spacing, 500 nm far spacing, 800 nm process-core width, and growth ratio 1.4.

For all thin cases:

- 38.89% of elements are >0.2 um from the process plane, representing 90.48% of volume.
- Fine in-plane rows remain present at these distant heights.

For the fabricated p5 cases:

- 24.89% (model 08) / 25.70% (model 09) of elements are >0.2 um from process planes,
  representing 87.54% / 87.87% of volume.
- Elements within 0.2 um of process planes but >0.8 um from physical edges account for
  6.69% / 12.02% of all elements. Thin layers are propagated into these regions too.

However, the mesh is NOT uniformly over-refined everywhere: the region both >0.2 um from
planes and >0.8 um from edges contains only 2.08% of elements but 58.34% of volume in
model 08 fabricated, and 4.63% of elements but 68.35% of volume in model 09 fabricated.
Simple far-field size inflation alone therefore cannot be assumed to yield a large gain.

## Recommended controlled experiments

Keep the current benchmark untouched and freeze p-order, physical geometry, trace basis,
quadrature, material data, and solver tolerances for the experiments.

1. **Plan-shape repair:** align neighboring graded-row nodes or locally retriangulate near
   corners. Preserve physical boundaries and the 2 nm requested normal spacing. Measure
   minimum angle/altitude, condition-number tails, DOFs, and source solve times. Repair may
   add a few elements while removing an artificial sub-nm altitude; that is acceptable.
2. **Height-dependent plan coarsening:** stop extruding the same fine plan mesh to the
   outer top/bottom. Use validated conforming transitions or supported nonconforming
   refinement. Do not coarsen the physical boundary layers or matching-surface trace
   representation blindly.
3. **Distance-dependent tangential grading:** compare 100 nm near-edge spacing with graded
   larger spacing away from corners/interaction regions while retaining near-wall normal
   resolution. The full library includes localized boundary excitations, so smooth-field
   convergence alone is insufficient.

Use representative hard/easy original nodal sources, independent held-out excitations,
and eventually full domain and per-interface matrix comparisons. Report wall time at the
same rank count, not just iteration count. No speedup or retained accuracy from these mesh
changes has yet been demonstrated by this audit.
