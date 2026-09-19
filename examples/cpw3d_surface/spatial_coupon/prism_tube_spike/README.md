# Prism edge tube feasibility SPIKE

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

**Status: feasibility spike (supervisor decision 37, 2026-09-18). Not production. Not part
of the canonical pipeline, its manifests, gates or audits.** The production mesher
`../mesh_spatial_coupon.jl`, the metric/adaptation stages and their defaults are untouched;
this directory only `include`s the production CAD helpers.

## What it builds

A four-edge spatial coupon whose straight metal edges (top edge `z = thickness` and bottom
edge `z = 0` of every metal segment) are each surrounded, on the dielectric side, by a
**localized prism tube**:

- 2D cross-section meshed once: geometric rings of radial size 0.25/0.5/1/2/4/8/16 nm
  (ratio 2, cumulative radius 31.75 nm), 9 angular sectors of 30 degrees over the 270 degree
  dielectric side; 71 nodes / 117 triangles per section. The top edge sees vacuum only; the
  bottom edge sees 3 substrate sectors under the metal (bounded by the metal bottom face,
  label 5001, and the trench wall, label 3100) and 6 vacuum sectors (trench and above).
- Extruded along the edge at the production tangential spacing 50 nm -> prisms
  (117 per layer). The tube stops at the 4 nm-graded corner balls (0.1 um from every
  semantic corner, tets there as in production) and runs to the outer box at continuation
  vertices.
- Tube/tet interface: every lateral quadrangle of the outermost prism ring carries an
  **explicit pyramid** (apex 8 nm outside the quad on the sector bisector), so the surface
  presented to Gmsh's 3D Delaunay mesher is all triangles and it meshes pure tetrahedra.
  This is the robust choice: letting Gmsh close the quadrangles with its own pyramids
  (`Mesh.Algorithm3D = 1` does support that) was tried first and Gmsh's tetrahedral optimizer
  then produced duplicated / overlapping tetrahedra next to the pyramids (2 duplicates,
  4 faces shared by 3-4 elements in the tiny case); with the explicit pyramids the mesh is
  conforming with the optimizer on.
- Everything outside the tubes: the production size fields (anisotropic curve attractor
  NormalSize 25 nm / TangentialSize 50 nm / FarSize 160 nm, the 4 nm corner grading law,
  the trace-basis size callback) and the production labels (1 / 3000 / 3100 / 5001 / 6001,
  materials 1 substrate / 2 vacuum), Gmsh 2.2 binary output. Gmsh only, no MMG.

Tube radius 31.75 nm rather than the 50 nm protected band: the trench floor lies 50 nm
below the bottom edge, so a 50 nm tube (or its pyramids) would touch it; 7 rings of the
requested sizes end at 31.75 nm and leave 10 nm between the pyramid apexes and the floor.

## How the explicit mesh survives Gmsh's generator

Gmsh's `generate()` deletes every face mesh before its 1D pass and every volume mesh before
its 2D/3D passes, `Mesh.MeshOnlyEmpty` protects only the dimension being meshed, and
`Mesh.Renumber` (default on) renumbers nodes after each pass. The tube mesh is therefore
installed in three phases (`../prism_edge_tubes.jl`): points and curves before `generate(2)`; the
faces (cap triangles, radial quadrangles on the metal / trench-wall strips, pyramid
triangles on the lateral faces) after `generate(2)`; the OCC tube volumes are removed
(non-recursively) before `generate(3)` and replaced afterwards by discrete volumes carrying
the interior nodes, prisms and pyramids. `Mesh.Renumber = 0` keeps the explicit node tags
valid. Tube CAD entities are matched to the fragmented OCC entities by centroid and the
match must be one-to-one (a split tube entity is an error).

## Files

| file | role |
|---|---|
| `../prism_edge_tubes.jl` (moved out of the spike for the Gmsh-only production mesher, decision 38) | cross-section, tube frames, OCC tube volumes, entity matching, three-phase explicit mesh installation |
| `hybrid_census.jl` | element counts / quality per type (Gmsh SICN, gamma), label areas, duplicate removal, JSON writer |
| `mesh_tiny_hybrid.jl` | Step 0: small box + one metal ridge on an etched substrate, `hybrid` or pure `tets`, Gmsh 2.2 |
| `mesh_prism_tube_coupon.jl` | Step 1: the four-edge hybrid coupon from the frozen inputs of a case directory |
| `check_hybrid_mesh.py` | independent check of the written file (meshio): counts, orientation, conformity, coverage, areas, corner scaled Jacobian / condition per type, SHA256 |
| `test_prism_tube.jl`, `test_check_hybrid_mesh.py` | targeted tests |

```
julia --startup-file=no --project=test/examples examples/cpw3d_surface/spatial_coupon/prism_tube_spike/mesh_prism_tube_coupon.jl \
    examples/cpw3d_surface/spatial_coupon/testdata/four-edge-9d2cb9bbb3fe out.msh --census census.json
python3 examples/cpw3d_surface/spatial_coupon/prism_tube_spike/check_hybrid_mesh.py out.msh report.json \
    --expect-area 3000=11.925867628990442 --expect-area 3100=210.56267181021096 --expect-area 5001=36 --expect-area 6001=39.6
```

## Results (2026-09-18, HEAD 026f42a8a, Gmsh 4.15.0-git, local Palace build v0.17.0-370-gc0f09d47e)

Step 0 (tiny case, terminal = top/bottom box faces, ground = metal): the hybrid mesh
(20,291 tets + 9,360 prisms + 720 pyramids) loads in Palace with `Tetrahedron / Prism /
Pyramid` geometries, converges, and the surface participations are finite with quadrangle
boundary elements on 6001 / 5001 / 3100. p_MA at p4: hybrid 7.837e-5, pure tets with the
same size field (25 nm at the edge) 7.278e-5 (-7.1 %), pure tets refined to 8 nm at the
edge 7.486e-5 (-4.5 %); the hybrid's p-sequence is the most converged (p2 -> p4: +1.2 %
vs +4.9 % / +2.8 %). PCG 26 vs 9-10 iterations.

Step 1 (four-edge): 1,588,715 elements = 1,410,299 tets + 165,672 prisms + 12,744
pyramids, 352,379 nodes, Gmsh 2D + 3D in 25 s (mesher < 4 min single-threaded including
Julia compilation). Conforming (every face shared by exactly 1 or 2 elements; every
quadrangle by prism+prism or prism+pyramid; every labeled surface element a face of the
right number of volume elements), every element positively oriented, no duplicates.
Label areas equal the reference invariants exactly (1: 777.2, 3000: 11.925867628990439,
3100: 210.56267181021096, 5001: 36.0, 6001: 39.6 um^2; volumes 514.396 / 544.404 um^3).
Quality per type (corner scaled Jacobian / Jacobian condition; report only): prisms
0.228 min / 590 max (the 380:1 longitudinal anisotropy of the innermost ring; Gmsh SICN
0.007), pyramids 0.29 / 8.7, tets P1 0.106 / P50 3.6, min 1.7e-5 / max 1774 - all 702 tets
below 0.01 (total volume 3.9e-6 um^3) and the 5 above condition 1000 are the corner-ball
tets glued to the 0.25 nm cap triangles at the tube ends.

Local Palace solve of the four-edge hybrid mesh, source 7 (edge-interior control), p2,
6 ranks: 2.72M H1 unknowns, PCG 23 iterations, 95 s, 10.2 GB. Against the graded_v2
reference (8.74M tets, p4): E_elec +0.19 %, E_MA -0.39 % (p_MA -0.58 %), E_MS -0.14 %,
E_SA(total) -4.5 %. For comparison at the same source the tet meshes gave p_MA
EL4c p3 -9.16 % / p4 -7.16 %, EL1c p3 -6.89 % / p4 -5.52 % / p5 -4.33 %.

## Known limitations / follow-ups

- Tube caps at the corner balls glue 0.25 nm triangles to 25 nm tets: ~700 sliver tets per
  coupon (condition up to 1774). A graded cap fan or extending the tube into the ball would
  remove them; not done in the spike.
- Prism anisotropy 380:1 at the innermost ring: PCG iterations ~2.5x the tet mesh at p4 on
  the tiny case; 23 at p2 on the four-edge coupon.
- The four-edge coupon only (single slot, single conductor, sharp vertical fabricated
  geometry, explicit etch footprint carrying every metal edge); the tube runs to the box at
  continuation vertices and stops 0.1 um before semantic corners.
- Step 2 (2D cross-section MA per unit length with the cpw2d machinery) was not done.
