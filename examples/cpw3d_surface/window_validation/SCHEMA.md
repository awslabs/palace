<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Polygon-set input of the fabricated window mesher (`Version` 1)

The input of `mesh_polygon_window.jl` (read by `PolygonWindowMesh.read_polygon_set`): one JSON
object describing a validation window in plan view. Lane W (window extraction) writes it; the
mesher builds the fabricated three-dimensional mesh from it. Lane W owns the geometry decisions
that happen BEFORE this file (the window box, the clipping of the chip's metal at the box,
termination of cut traces, which bodies are terminals, bump halos); the mesher only refuses
inconsistent input.

## Units and coordinates

  - Micrometres everywhere (lengths, coordinates, z).
  - Plan coordinates `[x, y]` are absolute chip coordinates (the identification manifest's).
  - z is the chip normal: `SurfaceZ` of a plane is the substrate surface of that plane's chip;
    a plane facing `up` has its substrate below `SurfaceZ` and its metal in
    `[SurfaceZ, SurfaceZ + MetalThickness]`; a plane facing `down` (the flip-chip L2) has its
    substrate above `SurfaceZ` and its metal in `[SurfaceZ - MetalThickness, SurfaceZ]`.

## Top-level fields

| field | type | required | meaning |
|---|---|---|---|
| `Version` | integer | no (default 1) | schema version; only 1 is accepted |
| `Name` | string | no (default `window`) | the window name (manifest, Gmsh model name) |
| `Box` | `{"X": [x0, x1], "Y": [y0, y1]}` | yes | the plan rectangle of the truncated box, `x0 < x1`, `y0 < y1`; every polygon must lie inside it (vertices may lie ON the wall; a polygon leaving the box is refused) |
| `Process` | `{"MetalThickness": 0.1, "Overetch": 0.05}` | no (defaults shown) | metal thickness and overetch recess of the exposed substrate, both planes |
| `Planes` | array of 1 or 2 plane objects | yes | the metal levels (below) |
| `Bumps` | array of bump objects | no (default none) | excluded metal columns between the two planes (below); need two planes |
| `Vacuum` | `{"Below": b, "Above": a}` | no (default 0 / 0) | vacuum beyond the outermost substrate backsides (0 = the backside is the box wall, homogeneous Neumann; > 0 = the backside is a `substrate_backside` face against vacuum). A single `up` plane needs `Above` > the metal thickness (the transmon: `Below` 475, `Above` 1000) |
| `Terminals` | array of strings | no (default: every non-ground label sorted) | the attribute order of the terminals; must list every non-`ground` conductor label exactly once |
| `MatchingRadius` | number | with two planes (unless the mesher is run with `--cross-plane-snap-um`) | the identification's matching radius R (um, e.g. 1.9); sets the cross-plane snap distance delta = 0.05 R (below) |

## Plane object

| field | type | required | meaning |
|---|---|---|---|
| `Name` | string | yes | e.g. `L1`, `L2` (used in the manifest's per-plane statistics) |
| `SurfaceZ` | number | yes | substrate surface z |
| `Facing` | `"up"` or `"down"` | yes | see units / coordinates; with two planes the lower one must face `up` and the upper one `down`, and their `SurfaceZ` must differ by more than twice the metal thickness |
| `SubstrateThickness` | number | yes | thickness of this plane's substrate (its backside is at `SurfaceZ -/+ SubstrateThickness` for `up` / `down`) |
| `Polygons` | array of polygon objects | yes (at least one) | the metal footprint of this plane |

## Polygon object (a metal body of one plane)

| field | type | required | meaning |
|---|---|---|---|
| `Conductor` | string | yes | conductor label: `ground` is the shared ground (every plane, every bump; attributes 4 / 5); every other label is a terminal with its own metal-air / metal-substrate attribute pair (7 / 8, 10 / 11, ...). The same label on two planes or in several polygons is ONE conductor (e.g. an L1 ground and an L2 ground, or a body that a bump joins) |
| `Outer` | `[[x, y], ...]` | yes (>= 3 vertices) | the outer loop, straight edges, not self-intersecting, the last vertex is NOT repeated |
| `Holes` | `[[[x, y], ...], ...]` | no (default none) | hole loops (>= 3 vertices each), strictly inside the outer loop, pairwise disjoint; the metal is `Outer` minus the holes |

Rules:

  - Loop orientation is free (either winding for `Outer` and for every hole): the mesher
    orients every loop counterclockwise itself before passing it to Gmsh OCC (which expects
    the hole wires in the OUTER wire's orientation). Do not rely on a winding convention.
  - Curved metal edges are polygonised by the writer at (or below) the tangent target t that
    the mesh will be generated with (the transmon export uses the transfinite chords at t = 5):
    every polygon edge becomes a transfinite curve of `ceil(length / t)` segments, so a chord
    longer than t is subdivided and a chord shorter than t stays one segment.
  - Polygons of one plane must not overlap (a shared edge between two bodies is also refused:
    two conductors must not touch). Polygons of different planes may overlap freely (that is the
    flip-chip stack).
  - Metal that reaches the window wall is represented by vertices ON the wall (a ground plane
    `Outer` equal to the box rectangle with holes is the common case). Edges lying on the wall
    are not metal edges (no boundary layer; the wall is `exterior_boundary`, attribute 3).
  - Every polygon must lie inside the box (clip at the box BEFORE writing; the mesher refuses a
    polygon that leaves the plan rectangle).
  - Cross-plane reconciliation (supervisor decision on the S1p trial): edges of the two planes
    that are nominally coincident in plan (aligned ground edges, the rounded corners of both
    chips) usually arrive with different vertex samplings and produce sliver partitions. The
    mesher snaps every vertex of the UPPER plane within delta = 0.05 x `MatchingRadius`
    (0.095 um at R = 1.9, the identification's joint-noise resolution; cross-plane edges are
    separated vertically by the chip gap, so smaller plan offsets are physically irrelevant)
    onto the lower plane (vertices first, then segments; ties by
    coordinates) and inserts the cross vertices on both chains, so coincident runs become
    identical point sequences; the lower plane's geometry never moves. Within ONE plane
    nothing is snapped: two polygons closer than delta (a sub-delta slot, touching metal) are
    REFUSED — the writer must resolve them. The manifest records the moved vertices, the
    maximum displacement, the inserted vertices and the coincident run length
    (`cross_plane_reconciliation`): the reference differs from the (unreconciled) thin
    geometry by at most delta on those runs.

## Bump object

| field | type | required | meaning |
|---|---|---|---|
| `Conductor` | string | no (default `ground`) | the conductor of the column (and thus of the bodies it joins on both planes) |
| `Footprint` | `[[x, y], ...]` | yes (>= 3 vertices) | the plan footprint; it must lie on metal of the SAME conductor on both planes (a footprint not on metal of both planes is refused) and footprints must not overlap |

The bump is an excluded metal column from the L1 metal top (`SurfaceZ_1 + MetalThickness`) to the
L2 metal bottom (`SurfaceZ_2 - MetalThickness`); its sidewall is a metal-air face of its
conductor. The bump height is therefore implied by the two planes (no `Height` field).

## Example

```json
{"Version": 1, "Name": "S1p",
 "Box": {"X": [290.0, 765.0], "Y": [-200.0, 80.0]},
 "Process": {"MetalThickness": 0.1, "Overetch": 0.05},
 "Planes": [
   {"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "SubstrateThickness": 525.0,
    "Polygons": [{"Conductor": "ground",
                  "Outer": [[290.0, -200.0], [765.0, -200.0], [765.0, 80.0], [290.0, 80.0]],
                  "Holes": [[[300.0, -10.0], [400.0, -10.0], [400.0, 10.0], [300.0, 10.0]]]},
                 {"Conductor": "trace_l1",
                  "Outer": [[300.0, -2.0], [400.0, -2.0], [400.0, 2.0], [300.0, 2.0]]}]},
   {"Name": "L2", "SurfaceZ": 4.8, "Facing": "down", "SubstrateThickness": 300.0,
    "Polygons": [{"Conductor": "ground",
                  "Outer": [[290.0, -200.0], [765.0, -200.0], [765.0, 80.0], [290.0, 80.0]],
                  "Holes": [[[500.0, -10.0], [600.0, -10.0], [600.0, 10.0], [500.0, 10.0]]]},
                 {"Conductor": "island_331",
                  "Outer": [[500.0, -2.0], [600.0, -2.0], [600.0, 2.0], [500.0, 2.0]]}]}],
 "Bumps": [{"Conductor": "ground",
            "Footprint": [[541.5, -191.2], [561.5, -191.2], [561.5, -171.2], [541.5, -171.2]]}],
 "Vacuum": {"Below": 0.0, "Above": 0.0},
 "Terminals": ["trace_l1", "island_331"]}
```

The resulting attribute table is written to the manifest (`OUT.json`, key `attributes`): 3D 1
substrate / 2 vacuum; 2D 3 `exterior_boundary`, 4 / 5 `ground_air` / `ground_substrate`,
6 `substrate_air`, 7 / 8 `trace_l1_air` / `_substrate`, 9 `substrate_backside`, 10 / 11
`island_331_air` / `_substrate`.
