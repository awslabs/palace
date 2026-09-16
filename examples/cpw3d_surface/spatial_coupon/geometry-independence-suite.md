<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is a fail-closed inventory and gate contract
for the mesh-only phase of the coupon-library acceleration plan. It does **not**
qualify response fields, trace accuracy, a production library, or all release
gates.

## Frozen inputs and required matrix

Every runnable case requires SHA-256-frozen `Signature`, `Boundary`, `Mask`,
`Process`, `SemanticContract`, and `MeshRecipe` roles. Expected materials,
labels/adjacency, corners, and protected supports come only from that contract.
The device etch footprint is a recorded choice, never an omission: a case
declares either a SHA-256-frozen `RetainedEtch` file (`retained-etch.csv`, the
same plan-view loop format the graded_v2 producer consumed through
`--etch-boundary`) or `EtchFootprint: "producer-default"`; a case with neither,
both, or another word fails preflight. When a `RetainedEtch` is declared, seed
generation binds it as the `source-retained-etch` input and its command must
pass exactly that path to `--etch-boundary`; without one, `--etch-boundary` is
forbidden and the seed census records `EtchBoundary: "producer-default"`. The
census also records the footprint's SHA-256 and the per-label interface areas
(`InterfaceAreas`, um^2) measured from the seed, so the footprint is asserted
from the mesh (reported, not gated; no expected-areas input exists). Only remote
campaign input `07` (four-edge) carries a `retained-etch.csv`
(`49fe5072...86d6`; graded_v2 `campaign.json` binds it for `07-fabricated`
only); `05`, `06`, `08`, `09` (ten-edge) and `10` have none and their graded_v2
references used the producer default, which the ten-edge case therefore records
explicitly. With the bound footprint the four-edge seed's areas are 3000 =
11.9259 um^2 and 3100 = 210.5627 um^2 against the graded_v2 reference 11.926 /
210.563 (the producer default gave a 2 x 2 square and 218.0). Exporting a
producer-default footprint as a bound file is a P2 follow-up.
Every etch footprint polygon - each device loop, or each producer-default collar
and strip - is simplified before any CAD face is created from it
(`simplify_footprint_polygon`): consecutive edges are merged while every vertex
between their outer endpoints lies within `FOOTPRINT_COLLINEAR_TOLERANCE` (the
same 1e-6 as `COPLANAR_TOLERANCE`; the stage contract requires equality) times
the merged edge's length of the merged edge, so there is one CAD face per
genuine facet and no near-coplanar sliver walls (the four-edge device footprint
carries a CSV-precision chain, loop 2 vertices 8-10, whose walls differed by
2e-8 to 1.6e-6 in normal). A kept vertex bends the wall by more than the
tolerance, so the metric stage always sees it as a dihedral. The census records
every simplified polygon (`FootprintPolygons`: conductor, plane, hole, points,
removed vertex indices, maximum deviation with its local scale) and a summary
(`FootprintSimplification`); the contract requires the deviation to stay within
tolerance x local scale (reported otherwise, not gated).
Modeling statement (supervisor decision 18(c)): etch footprint edges - the
trench wall/floor and wall/surface junctions of the device retained etch or of
the producer-default collars - are physical dielectric step edges. They are
legitimate feature segments alongside the signature (metal) edges, and they
legitimately carry the `NormalSize` band (the seed-derived `PhysicalSegments`
already give it to them; SA sensitivity at trench/cut junctions per the physics
pilot). No size is attached to them separately. The metric stage binds the seed
census and records the simplified footprint edges in the recipe as
`FootprintSegments` (provenance = the bound retained-etch SHA-256 or
`producer-default`, tolerance, polygon count, `[x0, y0, z, x1, y1, z]` segments
on the process plane); the stage contract requires them to equal the census's
polygon edges. The trace-diagonal detector treats a band aligned with a
signature edge or a footprint edge as a feature band and flags only bands
aligned with neither (`AlignedWithFeature`; `FeatureSegments` records both
counts and `FootprintSegmentProvenance` the source). The four-edge device
footprint's tilted trench walls, which the detector previously counted as 42
diagonal bands, are footprint edges.
Modeling statement (four-edge physics main run, 80 sources vs the graded_v2
reference): the lines where the Dirichlet cut surface meets a dielectric step or
material interface - the trench floor and walls and the un-etched
substrate-vacuum plane meeting the coupon box - are the recipe's dominant
under-resolved feature (13 wide hats on the z = -0.05 / 0 rings off by 0.4-4% in
energy and up to 25% in SA participation; one ordinary AMR cycle put 100% of its
235 marks within 0.3 um of the cut/trench junction line and repaired most of it).
They are therefore feature lines with the process-edge band. The seed makes
every curve of a material-interface surface (substrate volume on one side,
vacuum on the other) that lies on the outer box a feature curve
(`JunctionCurves` in the census: count and total length), so the frozen cut
surface carries the anisotropic band; the metric stage derives the same lines
from the seed's shared edges (a non-coplanar edge with a cut-surface triangle
and a material-interface triangle; roles from the contract's `CutSurfaceRoles`
and two-material `AdjacentMaterialSets`, `material_interface_attributes`),
chains them like the physical graph and records them as `JunctionSegments`
(segments, count, total length, label sets, provenance, rule). They receive
exactly the `PhysicalSegments` band law (NormalSize transverse band with the
ProtectedDistance/FarGrowth grading, TangentialSize along the line,
SurfaceProtectionRadius freeze); `PhysicalSegments` are unchanged and precede
them in the metric intersection order (`BandSegmentOrder`). The stage contract
requires the recipe's segments to be self-consistent, to name the contract's
cut and interface labels, and to total the census's junction-curve length
within `COPLANAR_TOLERANCE`. The trace-diagonal detector treats a band aligned
with a junction line as a feature band (`AlignedWithJunctionSegment`;
`FeatureSegments.Junction`). Cut/cut box edges and cut/conductor edges are not
junctions (the AMR probe marked nothing there). No size was added: the four-edge
junction set is 22 segments totalling 60.5 um on the box faces.
Trace-aware cut-surface sizing (same physics run): nine 22-48 nm-wide hats at
x ~ 1.93-2.0 on the y = 8 side are unresolved by the ~0.13-0.16 um cut-surface
elements (energy up to -21%), and the z = -0.05 and z = 0 source rings are only
0.05 um apart. The cut surface is frozen from the seed, so the seed sizes it from
the bound trace basis: a case that freezes `BasisContract`, `TraceVertices`,
`TraceTriangles` and `ProcessLibrary` (all four or none; `TRACE_BASIS_ROLES`)
binds them to seed generation and metric preparation as
`source-basis-contract` / `source-trace-vertices` / `source-trace-triangles` /
`source-process-library` (`--trace-basis-contract`, `--trace-vertices`,
`--trace-triangles`, `--process-library`, required together with the semantic
contract; forbidden otherwise). The basis vertices are in the process-library
frame and are placed in the mesh frame exactly as the campaign producer does
(`trace_basis.process_frame`: local z = process normal, local x = gap direction
of the first edge); the basis box must equal the coupon box and every vertex
must lie on the seed's cut-surface planes. Rule (`TraceBasisSizeRatio`, the
only new parameter, dimensionless, default 1.0 = at least one element per basis
edge; `--trace-basis-size-ratio`, passed identically by both stages): on the
cut surface the element size must not exceed the ratio times the shortest edge
of the basis triangle containing the point - a per-triangle rule, because the
hat of a basis vertex varies linearly over the whole incident triangle, so its
support is resolved where the hat varies only when the whole triangle is
discretized at that scale (the per-edge alternative resolves only the edges).
The seed applies it through the Gmsh size callback on top of the scalar
corner-isotropy background (`min(background, ratio x shortest edge + slope x
distance to the triangle)` with the process-band grading slope, capped by
`lc_far`; `Mesh.MeshSizeMin` follows the smallest requested size, the only
sub-`lc_fine` request any field can make) and records `TraceBasisSizing` in the
census (ratio, rule, frame, input digests, box, counts, basis edges below the
far size, minimum requested size, mesh-frame triangles, mesh size minimum,
slope). The metric stage binds the same files, requires the census record to
match (digests, ratio, triangles), caps the volume metric isotropically by
`min(FarSize, ratio x shortest edge + FarGrowth x distance)` so the existing
far/grading law is kept away from narrow hats (MMG's hmin = NormalSize still
floors the volume metric; the frozen cut triangles keep the seed's sizes), and
records `TraceBasisSizing` in the recipe with the statistics: unique basis
edges, edges below the far size, triangles below it, minimum requested size,
and `CutSurfaceSize` (min/median/max longest edge of the seed's cut triangles;
per narrow basis triangle the largest extent of the cut triangles centred in it
along its shortest-edge direction over the requested size - reported, not
gated). The stage contract requires both stages to bind the same basis or
none, the same ratio in both commands and both records, equal input digests,
equal mesh-frame triangles and the statistics; the manifest requires the seed
and metric inputs to be exactly the case's frozen files. Without a bound basis
nothing is sized, recorded or passed. The four-edge basis has 234 unique edges,
46 below the 0.16 um far size (minimum 0.0217 um, 76 of 156 triangles); the
ten-edge basis 774 edges, 146 below (minimum 0.0492 um, 244 of 516).
Supervisor decision 21 (trace-basis bands in the diagonal detector): the
four-edge narrow-hat basis triangles are slivers (0.022-0.048 um base, 11.4 um
long, fanning from (1.93..2, 8) to (-6, 0) on the bottom and top faces) and the
per-triangle rule puts a 0.02-0.05 um band along their long edges, which the
detector counted as two diagonal over-refinements (first build under the rule:
identity failed `trace-diagonal-overrefinement`, rotate-z passed only because
three edges crossed the 0.05 threshold at roundoff and split the component).
The detector exists to reject ARBITRARY trace-diagonal refinement (the
swept-prism era artifact); a band that lies on a bound trace-basis edge is
source-driven: the hat gradient across the sliver is part of the imposed
Dirichlet data along the whole sliver, resolving it is required to represent
the source (the -21%/-18% energies of sources 10/11/74/75 are exactly that
under-representation), and the graded_v2 reference resolved the same slivers
with CAD trace diagonals. The metric stage records the bound basis's unique
edges as `TraceBasisEdges` (input digests, count, metric-frame segments; the
contract requires them to be exactly the recorded basis triangles' edges); the
detector treats a line-like band as source-driven only when it lies ON one of
them - direction aligned AND both band endpoints within 2 x ShortEdgeThreshold
of the edge segment, not its line (`OnTraceBasisEdge`; position-aware, stricter
than the direction-only signature/footprint/junction checks) - and reports such
bands separately (`TraceBasisEdgeBands`: count, total span, maximum RMS width;
`LineLikeBandsAlignedWith` per feature class) so the arbitrary-diagonal guard
stays visible. Unaligned bands must still be zero (gate unchanged); with no
bound basis the behaviour is unchanged. P2 (recorded, not fixed): the
detector's components are formed from edges at or below the threshold, so a
band whose edges sit at the threshold can split into sub-half-diameter pieces
under a rigid motion (three edges of 16,819 crossed 0.05 by roundoff in the
rotate-z variant); a tolerance-consistent component rule is a follow-up.
Every case has identity and `rotate-z-0.63` variants with explicit transforms
and a fixed comparison pair. Concave/multislot, hole, rounded/filleted, and
opposed-layer controls are ordinary required cases. Feature-scaling and
CAD-subdivision-sensitivity comparisons are declared in the manifest.

The matrix contains 12 cases and 24 required case/variant entries. All 12 cases
now have hash-frozen local source contracts. A bounded read-only assessment on
`soca-green` copied only the approved source-contract files from campaign inputs
`07` and `09`; it copied no mesh, field, solution, response matrix, or source
bank. Campaign input `07` maps to
`spatialedgecluster_edgecount-4_9d2cb9bbb3fe` and has exactly four signature
rows. Input `09` maps to `spatialedgecluster_edgecount-10_6791f1c84123` and has
exactly ten. The process-library model names, local/remote SHA-256 values, trace
mesh references, slots, conductors, topology, and target-repository Apache-2.0
license were reviewed. Their semantic contracts derive label families and
slot/conductor coverage from the source process model and semantic corners from
plan-view vertices explicitly classified `Physical`.

## Preflight

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --preflight-only \
  --root /tmp/coupon-generality-preflight
```

This now succeeds for all 12 source cases and verifies the four-/ten-edge row
counts as four and ten. `--input CASE=DIRECTORY` cannot bypass a mismatched
hash.

## Production evidence chain

A version-3 normalized record is assembled from five separate records:

1. `bounded-run`: `run_bounded_mesher.py` records the exact command, its
   absolute working directory (so relative argv paths resolve unambiguously),
   constrained environment, launcher digest, resource measurements, and output
   mesh digest.
2. `mesh-topology-quality`: the real Gmsh mesh is parsed and audited for
   topology, labels/material adjacency, corners, protected labels, achieved
   directional widths, diagonal bands, and Jacobian quality.
3. `mesh-complexity`: H1 DOFs are counted from mesh topology/order and semantic
   feature/CAD-subdivision counts are derived from the signature.
4. `mesh-invariants`: material volumes and labeled surface areas are integrated
   from the mesh.
5. `variant-transform`: candidate points, connectivity, and labels are checked
   against the exact 4x4 transform of the independently audited identity mesh.

Every record binds the same mesh digest, every manifest source digest, exact
transform and digest, case/variant, command/environment, and frozen producer
digest. The bounded record consumes five separately executed reports forming a
digest-linked DAG: seed generation, metric preparation, native adapter/MMG,
label restoration, and final Gmsh publication. A tool is accepted only when its
frozen path occurs in that stage's executed argv. The DAG binds the seed mesh,
seed corner census, MMG seed, tensor metric, pins, fixed triangles, restoration
recipe, adapted mesh, restored mesh, final candidate, and ownership partition. Merely declaring
an adaptor/MMG file cannot pass.

The topology audit reads Gmsh physical volume names, derives ownership from the
publisher's partition report, compares protected surface measures to the seed,
measures corner/subdivision/cut neighborhoods, and detects connected short-edge
surface bands rather than asserting constants. Feature and CAD-subdivision
counts are frozen in `FeatureTopology` and independently reproduced from finite
oriented signature segments plus boundary `Physical`/`Continuation` classes.
Records and meshes must be content-distinct by SHA-256 across matrix entries.

Run `general_mesh_audit_producer.py` once for each required kind, passing all
five `--stage-report STAGE=PATH` bindings to `bounded-run` and
`mesh-topology-quality`. Normalize only those five audit outputs:

```sh
python3 normalize_general_mesh_evidence.py \
  geometry-independence-suite.json CASE VARIANT \
  /large/run/CASE--VARIANT/candidate.msh \
  /large/audits/CASE--VARIANT.json \
  --audit-record bounded-run=/large/run/CASE--VARIANT/bounded.json \
  --audit-record mesh-topology-quality=/large/run/CASE--VARIANT/topology.json \
  --audit-record mesh-complexity=/large/run/CASE--VARIANT/complexity.json \
  --audit-record mesh-invariants=/large/run/CASE--VARIANT/invariants.json \
  --audit-record variant-transform=/large/run/CASE--VARIANT/transform.json
```

One case can be verified independently of the matrix with
`verify_canonical_case_entries.py MANIFEST AUDIT_ROOT OUTPUT.json CASE` (a frozen
tool): it judges every manifest variant of the case with the unchanged
production functions (manifest validation, immutable-input hashes, evidence
gates, mesh readability, bound audit/stage records, exact canonical-build reuse)
and evaluates the frozen covariance comparison whenever both compared variants
have evidence - even when a variant fails - so the report lists every failure
(`Failures`, per-entry `GateFailures`/`Error`, `TransformComparisonFailures`,
`CanonicalReuseFailures`). The report is written whether or not the case
passed; the exit status is nonzero unless everything passed. It is a
single-case check, not the matrix, physics or release qualification.

After all 24 records exist:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --audit-root /large/audits \
  --root /large/coupon-generality-summary
```

No release claim is permitted until all production records and scaling
comparisons pass.

## Trace tooling boundary

The historical 135-column/40-constrained/95-free tools remain explicitly
outside this gate. `trace_audit_contract.py` requires unique excitation names,
explicit kinds, smooth/global and localized-matching roles, exact coverage of
every declared conductor/interface class, and a fixed finite-coefficient
combination. It remains only a future arbitrary-source seam; no three/six/ten-
edge physics claim is made.

## Canonical build and mandatory placement publication

Every proper-rigid placement now uses the seven-stage contract in
`mesh_stage_contract.py`. The first six stages are one reusable source-local
canonical build: identity source validation, seed generation, metric preparation,
native MMG adaptation, label restoration, and canonical Gmsh publication. The
seventh stage, `proper-rigid-publication`, is mandatory even for identity.

`canonical_mesh_build.py` computes the cache key from all immutable source,
process, semantic, recipe, and gate hashes plus the exact canonical tool-role
set derived from the six stages (runtime, validator, mesher, metric preparer,
adaptation wrapper, adapter, runtime-resolved MMG library, restorer, publisher).
Its build record is derived from the six validated stage reports: it binds the
exact sixteen canonical artifact roles (every stage output, including the seed
corner census, the adaptation receipt and the source-local restored mesh) and
each stage report's SHA-256, and rejects missing or extra roles. A placement may share the six
canonical reports only when the cache key, build hash, every canonical artifact
hash and every stage-report hash match. Final meshes, transform receipts,
ownership outputs, and audit records remain per-placement and content-distinct.
Canonical and placement seconds/RSS are reported separately.

Seed generation honors the same isotropic corner ball that metric preparation
prescribes. `mesh_spatial_coupon.jl` takes the canonical semantic contract
(`--semantic-contract`, bound to the `canonical-semantic-contract` input of the
`canonical-source-validation` stage) and `--corner-isotropy-radius`, and around
every contract `SemanticCorners` point (pulled back through the contract's
`RigidTransform`, which must equal the seed's `--rigid-transform`) it imposes the
isotropic size `--lc-fine` inside the radius, graded to `--lc-far` with the same
linear slope as the process edge band (`(lc_far - lc_fine)` over the attractor
transition width). The band field itself is evaluated unchanged through a
`MathEval` `min(F1, ...)` (a `Min`/`MinAniso` wrapper re-derives the anisotropic
attractor's size and was measured to change the band mesh, so it is not used).
Longitudinal feature curves whose transfinite `lc_tangent` spacing would cross a
corner ball are meshed explicitly instead (`Mesh.MeshOnlyEmpty`), so the
corner-adjacent ridge edges start at the isotropic size; curves out of reach keep
the ordinary transfinite spacing. The explicit placement keeps the transfinite
`lc_tangent` grid wherever the corner size law equals `lc_tangent` (beyond
`CornerLawReach` = radius + (`lc_tangent` - `lc_fine`) / slope from every corner)
and fills only the gaps between kept grid nodes and the curve ends by
size-weighted arclength quadrature of the law. Grid preservation is required, not
cosmetic: the interior node row of a ridge-to-ridge face (a metal sidewall) exists
only because the ridge nodes of its two longitudinal curves are exactly aligned
(0.1 x 0.1 squares whose 0.141 diagonals the face mesher splits); rows that are
misaligned by a different node count or phase produce 0.1 x 0.11 triangles that
need no interior node, leaving full-height slivers along the whole face. A 1D
mesh driven by the size field or a size callback cannot reproduce the exact grid
(the scalar anisotropic attractor evaluates to sqrt(`lc_tangent` x `lc_fine`) on
the curve, and integral inversion shifts the phase), so it was measured to remove
the sidewall interior row exactly as the first, non-grid-preserving explicit
placement did. No new sizes are introduced. The stage contract requires the seed's
`--lc-fine` and `--corner-isotropy-radius` to equal the recipe's `NormalSize` and
`CornerIsotropyRadius`, and the seed's recorded `seed-corner-census` artifact
(`--corner-census`) to name the recipe's `TruePhysicalCorners` with the same size
and radius and to carry a `LongitudinalFaces` census; a seed stage without these
bindings, or with different values, is rejected. The census (edge
minimum/median/maximum inside each corner ball, the count and fraction of ball
edges longer than sqrt(2) x `lc_fine`, the corner-incident maximum aspect and
ring radii; per face bounded by at least two longitudinal curves, the interior
node count and along-edge histogram away from the corners and the count of
full-height triangles, i.e. triangles whose vertices all lie on the bounding
curves and span two longitudinal curves) is reported so the property is asserted
rather than assumed; it is **not** a pass/fail gate. Gmsh's Delaunay volume mesh
still leaves a tail of interior ball edges up to about twice the size (the
surface and ridge edges sit at the target), so MMG may still split interior edges
in the balls.

Metric preparation classifies seed surface features by geometry, not by label:
`edge_volume_metric.surface_features` makes a shared triangle edge a feature
(MMG ridge, pin at graph turns/junctions, anisotropic `PhysicalSegment`,
protected band) only when its incident triangles are not coplanar within the
fixed dimensionless `COPLANAR_TOLERANCE` (sine of the normal angle and plane
offset relative to the local triangle diameter). A reference change between
coplanar triangles, such as a slot-ownership split of one plane, remains a
plain reference boundary without ridge, metric or protected-band constraints;
per-support relabeling and projection in the restorer are unchanged. MMG,
however, reconstructs every reference boundary as a discrete geometric curve,
and with the 1e-8 Hausdorff bound an unpinned seam bend is refined without
limit (a ten-edge attempt reached 23.6M tetrahedra and aborted). Pins are
therefore the turns and junctions of the complete reference-boundary graph
(edges whose incident triangles differ in reference or are not coplanar):
`geometric-corner` pins come from the feature graph alone, `reference-turn`
pins only from coplanar seams; straight seam vertices are not pinned. The
recipe records `PinnedVertexKinds`, `PinnedVertices` and `PinPolicy`.
Straightness (`_continues_straight`, also used for chain merging in
`feature_chains`) admits a deviation of at most `COPLANAR_TOLERANCE` radians
(sine of the angle, 1e-6), which is *stricter* than the former
`cos > -1 + 1e-10` collinearity test (about 1.4e-5 rad): a seed whose
straight-edge vertices carry noise between 1e-6 and 1.4e-5 rad would now
produce extra chain breaks (spurious feature corners, pins and segments), not
merged chains. The canonical seeds are exact planes far below 1e-6 and their
segments, pins and corners are unchanged.
Reference-turn pins are MMG required vertices only; they never enter the
semantic-corner set, the protected supports, or the ownership/covariance audits
as geometry.

The same tolerance is the only rule for deciding that two triangles lie on one
plane anywhere in the chain: `edge_volume_metric.plane_deviation` (the larger of
the sine of the normal angle and the offset of a point from the plane relative
to max(local triangle diameter, distance to the reference point) - dimensionless,
origin-independent, rigid-covariant), `cluster_coplanar_triangles` (first-fit
clustering per label) and `match_equivalent_planes` (one-to-one pairing of the
clusters of two meshes). The metric stage's `PlanarSupports`, the planar-patch
areas of `audit_edge_metric_mesh.analyze`, the protected-surface audit patches
and the diagonal-band detector's planes are all these equivalence classes; no
plane is keyed by rounded coordinates any more (8-digit rounding fragmented the
ten-edge seed's `x = 9.8333...` matching plane into seven keys from 1e-8
normal noise). The metric recipe records `PlanarSupportEquivalence`.

Protected labels that occupy one plane and are joined by coplanar label seams -
edges shared by exactly two surface triangles of one plane class whose labels
differ while their contract roles differ only in the slot index
(`etched-substrate-vacuum-slot-0`/`slot-1`, `conductor-k-slot-i-ms`/`-ma`
across slots) - are audited as their union on that plane (supervisor decision
20): the per-triangle slot partition of one physical surface is label-only,
response ownership is certified point-wise by the quadrature ownership audit,
and another triangulation cannot reproduce a per-triangle partition to 1e-8
(ten-edge trench floor: 3100 = 6.4855 vs 6.5087 um^2 per label, union 193.0
both). Roles that differ otherwise are never unioned, and an edge carrying a
third, non-coplanar surface triangle (a metal footprint edge with its sidewall)
is a feature line, not a seam. Measure, boundary and topology thresholds are
unchanged; the per-label areas and seam length of every union are recorded as
`CoplanarSlotUnion` diagnostics.
Metric preparation also protects the seed's corner balls: every seed surface
triangle with at least one vertex within `CornerIsotropyRadius` of a contract
semantic corner is a fixed (MMG required) triangle, in addition to the cut
surfaces and the protected edge bands. The seed's isotropic Gmsh corner ball at
`NormalSize` (asserted by the seed census) is the intended corner
discretization; MMG's anisotropic adaptation adds value along the edges, not at
the corners, and when the balls were left free MMG re-meshed the corner surface
and volume with ring edges down to `NormalSize` / 2 and corner-incident aspects
of 7-15 that the bounded restoration repair cannot reach (measured on the
grid-preserving seeds: four-edge corner 2 at 10.59, ten-edge 5 of 10 corners).
The volume inside a ball is still adapted to the same isotropic metric. Freezing
the balls does not by itself fix the corners: the frozen surface ring is
preserved (four-edge ring radii 0.0243-0.0259, seed corner cells of aspect
3.4-7.8) but MMG inserts free interior vertices at about `NormalSize` / 2
(0.012-0.017, below its own `hmin`) next to the required corner vertex, and the
resulting cells have aspects 7-12 (four-edge after MMG [7.19, 11.90, 10.59,
5.19]; ten-edge [10.99, 11.95, 5.28, 5.29, 10.19, 9.18, 5.94, 9.89, 7.38,
10.50]) that the bounded 0.75 x `NormalSize` smoothing cannot undo (four-edge
repair reaches 4.57 / 5.72 / 6.68 / 3.80). The edge bands, however, improve
further (four-edge restoration components 7 -> 5; ten-edge sub-0.02
scaled-Jacobian cells 23 -> 6), so the protection is kept. The
recipe records `ProtectedCornerBalls` (radius, rule, total and per-corner frozen
triangle counts); the stage contract requires it, with the recipe's
`CornerIsotropyRadius`, the recipe's semantic corners and counts covered by the
fixed-triangle total (reported, not gated). The far-field budget policy depends
only on the seed element count and is unaffected; the ownership and covariance
audits see the same frozen cut/matching surface everywhere, since the corner
balls on a cut were the only previously unfrozen cut triangles.

Label restoration first removes MMG's sub-`hmin` corner insertions
(`_collapse_corner_ball_vertices`): inside every contract corner ball
(`CornerIsotropyRadius`), a free interior vertex (no planar support, not pinned)
whose shortest incident edge is below the recipe `NormalSize` (the `hmin` MMG was
given) is collapsed onto the corner vertex, a non-free neighbor, or a free ring
vertex of the corner at or beyond `hmin`, whichever valid cavity has the best
worst cell; a collapse stands only when every remapped cavity cell keeps a
positive orientation and a scaled Jacobian of at least the 0.02 quality target
(new cavity cells have no original floor) and its worst cell is no worse than
the cells it replaces (a collapse onto the corner itself was measured to be
valid yet leave an aspect-54 sliver where a ring vertex gave 3.4), a corner's
collapses are committed together and rolled back if its corner-incident aspect
did not improve, no vertex moves and no boundary triangle changes. The only threshold
is the recipe `hmin` already handed to MMG; the cells that survive keep their
original floors and the unchanged 0.75 x `NormalSize` smoothing, 4.0 corner gate
and 0.01 scaled-Jacobian gate follow. The restoration report records
`CornerBallCollapses` (per corner: collapsed vertices, aspect before and after
the collapse, rolled-back count), `CollapsedCornerVertices`,
`CornerAspectsBeforeCollapse` and `CollapseMinimumSize`. Measured on the
frozen-ball WIP adaptations: four-edge corners [7.19, 11.90, 10.59, 5.19] ->
after collapse [5.13, 4.89, 4.77, 5.01] (27/12/35/20 collapsed vertices);
ten-edge [10.99, 11.95, 5.28, 5.29, 10.19, 9.18, 5.94, 9.89, 7.38, 10.50] ->
[6.18, 4.96, 5.28, 4.05, 4.46, 4.71, 4.42, 4.89, 4.77, 5.43] (15-25 collapsed
vertices per corner; corner 2 rolled back: its worst cell's free vertex sits at
0.034, above `hmin`). The smoothing then proceeds from these ranges. Label restoration compares repair displacements to the bound with a
roundoff-only relative allowance (`DISPLACEMENT_ROUNDOFF_TOLERANCE`, 1e-12) so
a vertex clamped onto the displacement ball is not misreported as an overshoot;
a true overshoot is a rejected component, never an exception. Semantic corners
are repaired in alternating one-ring/two-ring passes over the non-fixed
vertices inside the recipe's corner ball (`CornerIsotropyRadius`; both rings
are filtered by the same radius, so the two-ring pass is a superset of the
one-ring pass), chained on one candidate; the *best* chained candidate is
committed as a single transaction only when it satisfies the cumulative
displacement bound, the original incident floors and the `MaximumCornerAspect`
gate, and `AchievedAspect` reports that committed candidate. The initial point
of every chained pass is clipped onto the optimizer bounds so a parameter
saturated by the previous pass cannot be rejected by SciPy as infeasible on
roundoff. 0.95 x gate stays the optimizer objective, results in (target, gate]
are committed and counted as `CornerRepairsTargetMissedGateSatisfied`, results
above the gate are `RejectedCornerRepairs`, and a corner or component whose
vertices are all fixed is a recorded rejection rather than an exception. Pinned
vertices (`PinnedVertices`, matched on the native adapted coordinates) and
matching-surface vertices never move. The
final scaled-Jacobian, corner-aspect and displacement gates are unchanged.

Native adaptation must be launched through `run_native_mmg_adaptation.py`.
There is no caller-provided hmax: the wrapper reads
`FarFieldBudgetPolicy.EffectiveFarSize` from the bound restoration recipe,
requires it to equal the recipe's `FarSize`, passes that exact value to the
reviewed adapter, and emits an adaptation receipt containing the executed argv.
The wrapper also resolves the adapter's single MMG3D shared library from its
link table and rpath entries (`otool` on macOS, `readelf` elsewhere), requires
it to be the `--mmg-library` bound as the stage's `mmg-library` tool role,
removes dynamic-loader search overrides from the adapter environment, and
records the library path, link reference and SHA-256 in the receipt. The stage
contract requires the receipt's adapter and library hashes to equal the bound
tools, so replacing the loaded library changes `CanonicalBuildId`. Preflight
resolution:

```sh
python3 -c 'from run_native_mmg_adaptation import resolve_mmg_library; \
  print(resolve_mmg_library("$EDGE_METRIC_ADAPTER_EXE"))'
```

The normal, tangent, corner, protected-band targets and 4,000,000-element cap
remain unchanged.

Every stage command must name exactly its bound inputs and outputs: seed
generation binds the source signature/boundary/mask, the canonical semantic
contract (`--semantic-contract`) and its corner census (`--corner-census`),
metric preparation binds the seed mesh, the canonical semantic contract and
supports (`--semantic-contract/--transformed-supports`) and the seed census
(`--seed-census`), canonical publication binds
the source process/signature/boundary (`--process/--signature/--boundary`), and
proper-rigid publication binds source semantic/signature/boundary/mask/process,
the canonical build record and every ownership/semantic/support output. Named
options must equal the bindings and positionally consumed paths must occur in
argv; `run_bounded_mesher.py` refuses to launch otherwise. Source validation
and canonical publication still take a source directory as a positional
argument, but the directory binds nothing: the files those stages consume are
named only by their bound `--semantic-input/--signature/--boundary/--mask` and
`--process/--signature/--boundary` options.

`publish_rigid_coupon_mesh.py` requires a fresh output and a finite orthonormal
homogeneous transform with determinant +1. It changes only Gmsh 2.2 node
coordinates, preserving node tags, cell-block order, connectivity, physical and
geometrical tags, physical names, and point/cell/field data. It reconstructs the
transformed semantic and support contracts from immutable source, checks exact
`R*x+t`, positive orientation, quality and invariant measures, and invokes
`audit_rigid_coupon_ownership.jl` with explicit `--process`, `--signature` and
`--boundary` files. That audit classifies quadrature points in the inverse
source-local frame and rejects any final owner-label mismatch or
nonpositive/nonclosing owner partition. The receipt records the ownership argv,
runtime/auditor hashes, output hashes and every source-input hash; the stage
contract requires the argv options to equal the placement bindings and the
hashes to equal the bound inputs/outputs. Identity therefore gets a fresh path
and receipt; reflection remains rejected.

A complete four-edge and ten-edge evidence run still requires an executable
reviewed build of `adapt_edge_metric.cpp` through `EDGE_METRIC_ADAPTER_EXE` and
its resolved MMG3D library bound as `mmg-library`. The fixture tests build a
tiny native adapter linked against a fixture `libmmg3d`
(`testdata/build_tiny_native_adapter.py`), so they require a C compiler.
Use fresh directories, 1800 seconds and 8 GiB for each bounded stage (16 GiB for
the aggregate audit), and the unchanged 4M element cap. Build the four-edge
canonical candidate once, publish identity and rotate-z through stage seven,
then build only ten-edge identity. Do not run ten-edge rotation, reflection,
physics, p4/p5, library, or HPC work in this mesh-only gate.
