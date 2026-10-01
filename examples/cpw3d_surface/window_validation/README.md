<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Fabricated reference mesher for validation windows (polygon sets)

The single transmon's fabricated reference (`single-transmon-finitemetal-anisotropic-20260824/generate_shared_plan_mesh.jl`: plan mesh with edge-normal boundary layers and a fixed tangent spacing on every metal edge, swept through explicit z levels, 100-nm vertical metal excluded, 50-nm overetch, conforming metal-air / metal-substrate / substrate-air shells) generalised to the validation plan's windows (`VALIDATION-PLAN.md` section (c), USER decision 184): a polygon set instead of the DeviceLayout geometry, any number of labelled conductors, one or two metal levels, bumps, a truncated box.

| file | role |
|---|---|
| `PolygonWindowMesh.jl` | the mesher (module): schema reader, plan fragmentation and classification, per-partition one-sided boundary-layer meshing, welding, z-level stack, sweep, attribute writing, manifest |
| `mesh_polygon_window.jl` | CLI: `julia --project=. mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUT.msh2 [--plan-only] [--exact-band-thickness]` |
| `validate_window_mesh.jl` | independent validation by physical NAME (adjacency rules, first-order simplices, no shared faces, minSICN, manifest counts, per-attribute areas / volumes) |
| `export_transmon_polygon_set.jl` | the transmon footprint as a polygon set (SingleTransmon environment; curved edges as their transfinite chords at the tangent target) |
| `synthetic_two_level_window.jl` | tiny two-level window (two facing CPW stubs + a bump) for the tests |
| `compare_mesh_validations.py` | two validation manifests of one geometry side by side (counts, areas, volumes) |
| `test_polygon_window_mesh.jl` | `julia --project=. test_polygon_window_mesh.jl` |
| `SCHEMA.md` | the polygon-set input schema (fields, units, orientation and clipping rules) for the writer (lane W) |
| `trial_window_extract_to_polygon_set.jl` | TRIAL-only converter of a lane-W `window_extract.py` extract (chip-wide loops, clipped at the box here by OCC) to a polygon set, for sizing / robustness trials of a window mesh before lane W's schema'd output exists (supervisor decision 187); never a stage-1 reference input |

## Polygon-set JSON

Micrometres, plan-view loops per plane with holes, conductor labels, bump footprints (the full schema: `SCHEMA.md`):

```json
{"Version": 1, "Name": "S1p",
 "Box": {"X": [x0, x1], "Y": [y0, y1]},
 "Process": {"MetalThickness": 0.1, "Overetch": 0.05},
 "Planes": [
   {"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "SubstrateThickness": 525.0,
    "Polygons": [{"Conductor": "ground", "Outer": [[x, y], ...], "Holes": [[[x, y], ...]]},
                 {"Conductor": "island_331", "Outer": [[x, y], ...], "Holes": []}]},
   {"Name": "L2", "SurfaceZ": 4.8, "Facing": "down", "SubstrateThickness": 300.0, "Polygons": [...]}],
 "Bumps": [{"Conductor": "ground", "Footprint": [[x, y], ...]}],
 "Vacuum": {"Below": 0.0, "Above": 0.0},
 "Terminals": ["island_331"]}
```

  - A plane facing `up` has its substrate below `SurfaceZ` and its metal in `[SurfaceZ, SurfaceZ + MetalThickness]`; facing `down` is the flip-chip L2 (substrate above, metal below). Two planes must be a lower `up` and an upper `down` plane.
  - `ground` is shared by every plane and bump (attributes 4 / 5); every other label is a terminal with its own pair (`Terminals` fixes the order; default sorted). Polygons of one plane must not overlap; two conductors must not touch. Vertices may lie on the box wall (wall edges get no boundary layer).
  - Bumps are excluded metal columns from the L1 metal top to the L2 metal bottom; their footprint must lie on metal of both planes.
  - `Vacuum.Below` / `Above`: vacuum beyond the outermost substrate backsides (0 = the backside is the box wall, homogeneous Neumann). A single `up` plane needs `Above` > 0 (the transmon: 1000 above, 475 below).

## Process, resolution, attributes

Metal 0.1 um with vertical sidewalls, overetch 0.05 um of the exposed substrate, metal volume excluded. Boundary layers of first height r (growth 2, band ~1.55 um) on every metal edge of every plane and bump in ONE plan mesh, tangent spacing t; metal and overetch bands resolved at r in z (`max(2, ceil(0.1 / r))` and `max(1, ceil(0.05 / r))` layers); far-field z levels at 0.1 / 0.2 / 0.5 / 1 / 2 / 5 / 10 / 20 / 50 / 100 um into each substrate and 0.15 / 0.2 / 0.3 / 0.5 / 1 / 2 / 5 / ... / 200 um into the vacuum (the transmon reference's stack; between two planes each side grades to the midpoint of the gap). Tetrahedra = plan triangles x (levels - 1) x 3 minus the metal.

| attribute | dim | group |
|---:|---|---|
| 1 / 2 | 3 | substrate / vacuum |
| 3 | 2 | exterior_boundary |
| 4 / 5 | 2 | ground_air / ground_substrate |
| 6 | 2 | substrate_air (top surfaces and overetch steps, both planes) |
| 7 / 8 | 2 | `<terminal 1>_air` / `_substrate` |
| 9 | 2 | substrate_backside (a substrate backside against vacuum, when `Vacuum` > 0) |
| 10 / 11, 12 / 13, ... | 2 | further terminals |

SA / MS / MA participations are surface integrals on these shells (Palace `Postprocessing.Dielectric` with `Thickness` 0.002), as in the transmon reference configs; the manifest (`OUT.json`) records the table, counts, per-attribute areas / volumes, first-layer heights, tangent statistics and z levels.

### Boundary-layer Thickness: a deliberate deviation from the recorded transmon generator

The recorded generator passes the Gmsh BoundaryLayer `Thickness` as the exact geometric sum r (2^n - 1) of the n layers, so floating-point rounding decides whether a column gets n or n - 1 rows (~8 % of the columns short at r10, ~20 % at r50, also in the recorded meshes; on Linux that mix fails Gmsh's edge recovery at r10, supervisor decision 188). This mesher passes the sum x (1 + 1e-6) by DEFAULT (every column exactly n rows, deterministic across platforms; manifest `radial_band_thickness_mode` = `geometric_sum_x_1p000001`); `--exact-band-thickness` (`exact_geometric_sum`) reproduces the recorded meshes. The r50 regeneration of the transmon footprint is recorded under both settings (`fabricated-window-mesher-20261002/REPORT.md`); the accepted transmon reference, computed on the recorded meshes, is unaffected.

Memory: the transmon at r50 / t5 needs ~9 GB (16.5 M tets), r10 ~15 GB; window-sized sets are small. Run full-size generations on a cluster node.
