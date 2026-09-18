<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is a fail-closed inventory and gate contract
for the mesh-only phase of the coupon-library acceleration plan. It does **not**
qualify response fields, trace accuracy, a production library, or all release
gates.

## Frozen inputs and required matrix

`geometry-independence-suite.json` is the only production manifest.

## Production recipe (supervisor decision 34B, 2026-09-17)

The production recipe is the EL4c calibration recipe, recorded in the manifest's
`ProductionRecipe` block and executed by every production stage command:
seed `--lc-tangent 0.05` (V1: the ridge grid at 2 x NormalSize), metric
`--far-growth 0.5` (V2: half the growth rate of the normal size beyond the
protected distance, widening the graded layer around every metal edge and
junction line from 0.135 to 0.22 um), the seeded transverse edge layer
`--edge-size 0.004 --edge-growth-ratio 2 --edge-layer-aspect 4` on seed and
metric (rows 4/12/28 nm below NormalSize 0.025 with nested tangential rows
12.5/25/50 nm on the 50 nm ridge grid; RequiredReach 0.0334; the adapter
`--hmin 0.004` = EdgeSize), corner grading `--corner-size 0.004` on seed and
metric (shells 4/8/16 nm at 4/12/28 nm inside every 0.1 um corner ball, rows to
the ball boundary, 0.1 um of un-layered edge per corner), the corner balls and
the layer as MMG required tetrahedra (decision 30) with the seed-side
required-region gates (`--maximum-corner-aspect 4 --minimum-scaled-jacobian .01
--maximum-jacobian-condition 1000 --maximum-quality-displacement-over-normal
.75`, equal on seed and restorer). Every other parameter is unchanged (NormalSize
0.25 x thickness, TangentialSize / CornerIsotropyRadius 4 x NormalSize, FarSize
0.08 x radius, protected distance and surface radius 2 x NormalSize,
TraceBasisSizeRatio 1, `--hgrad 1.15 --hausd 1e-8`, restoration bounds 0.25 /
0.75). Before 34B production executed `--lc-tangent 0.1`, `--far-growth 1.0`,
no edge layer (`--edge-size 0`), `--hmin` NormalSize and no corner grading;
those values are recorded as `ProductionRecipe.ValuesBefore34B` and as the
calibration cases' `ProductionValuesBefore34B`.
Binding: `verify_canonical_case_entries.py` and `run_general_mesh_suite.py`
(`general_mesh_manifest.validate_production_recipe_commands`) require the
recorded seed-generation / metric-preparation / native-adaptation-mmg commands
of every production case to execute each `SeedCommandOptions` /
`MetricCommandOptions` / `AdaptationCommandOptions` pair exactly once at its
value; the recipe's `EdgeLayer`, `CornerGrading` and `RequiredTetrahedra`
records are bound by the stage contract as before. The canonical cache key
does not encode recipe options (unchanged, documented limitation). A
calibration manifest never carries `ProductionRecipe` (`validate_manifest`
rejects it) and production cases still never carry `Calibration` or `EdgeLayer`
blocks.
Physical gates unchanged: MinimumScaledJacobian 0.01, MaximumJacobianCondition
1000, MaximumCornerAspect 4.0, MaximumProtectedMeasureError 1e-8, orientation,
ownership closure, MaximumElements 4,000,000, 1800 s / 8 GiB per stage (audits
16 GiB), trace-diagonal 0. The calibration-only EdgeLayerQualityRule and the
5,000,000 element cap stay in the calibration manifest only.
Evidence (four-edge coupon vs the graded_v2 reference, 8.74M tets; p4, 80
sources, 60 free-view): production before 34B (physics-02) E max 4.1%, SA max
15%, p_MA median -10.08%, p_MS -1.46%, 0.453 node-h; V1/V2 (physics-03) E 60/60
within 1% (max 0.95%), SA max 8.7%, p_MA -9.11 / -8.99%, 0.454 node-h; EL4
(physics-04) p_MA -7.05%, p_MS -0.81%, E/SA no regression, PCG 22.4, 0.749
node-h; EL1 1 nm x 50 nm (physics-05) p_MA -10.94% (worse), E/SA unchanged,
0.836 node-h, rejected; EL4c (physics-06, PBS 44717) E 60/60 within 1% (max
0.81%), SA max 8.56%, p_MS median -0.75%, p_MA median -6.63% / magnitude-
weighted -5.00% / strongest-20 -5.05%, PCG mean 22.6 / max 33, 0.752 node-h per
80-source p4 coupon, 3,471,507 tets (2.5x fewer than the reference). E and SA
never regressed across V1/V2/EL4/EL4c.
Design gate (supervisor decision 35, 2026-09-18): `achieved-anisotropy`
(MinimumAchievedAspect 1.5 and MaximumNormalFactor 2.0, values unchanged) judges
the band within one NormalSize of the physical segments outside the recorded
edge layer (`audit_edge_metric_mesh.achieved_anisotropy`). When the recorded
layer covers every cell of that band the gate is not applicable by construction
(`AchievedAnisotropy.Gate` = `not-applicable: layer-covered band`, `Samples` 0;
`general_mesh_manifest.layer_covered_band` accepts only the complete record):
the layer's design statement is the bound EdgeLayer aspect rule
(`AchievedAnisotropy.EdgeLayer`), and the layer-adjacent band - the band cells
within three NormalSize outside the layer - is recorded informationally
(`AchievedAnisotropy.LayerAdjacentBand`: cells, TangentialP50, transverse P90s,
`TransverseP90OverNormalSize` - the MaximumNormalFactor-equivalent - and the
nearest-vertex distance range to a span), never gated. Any one-NormalSize band
sample with cells outside a recorded layer (no layer, or an edge without one)
is judged by 1.5 and the normal factor as before. Rationale: the layer-adjacent
band is the transition shell between the 4 nm layer and the 25 nm band, graded
by the layer rows and the frozen 50 nm seed grid; its anisotropy is not a
design intent (EL4c identity: the one-NormalSize band is 26,870 layer cells and
0 others; the three-NormalSize band outside the layer is 5,955 cells with
nearest vertices 33-83 nm from a span, TangentialP50 50.0 nm, transverse P90
40.5/66.5 nm, 2.66 x NormalSize; the seed grid caps the band tangential at
50 nm, so 1.5 is unreachable by any sample of a lc-tangent-0.05 recipe - V2
measured 0.0571 / 0.0499 = 1.14); E/SA never regressed across V1/V2/EL4/EL4c
(physics-03..06). The covariance comparison compares the layer-adjacent band
of both placements when the gate is not applicable and requires both to agree
on whether it applied.
Trace-diagonal detector alignment (supervisor decision 36, 2026-09-18): the
edge layer's nested 12.5 nm rows make every metal edge a line-like short-edge
band on the metal faces (pre-34B the frozen 50 nm grid had none). On the
ten-edge coupon's 1.0 um edge x = -1 (y in [-0.6, 0.4]) two of the five
per-face bands (5001 at z = 0, 6001 at z = 0.1; RMS width 15-19 nm) had their
SVD principal axis tilted by 2.2-2.9 mrad from the signature direction - a
one-sided 1 um x 0.02 um point cloud resolves its direction only to about
width / span = 1.5e-2 - and the fixed cosine tolerance 1e-6 (1.4 mrad) counted
them as diagonal over-refinements while their siblings passed at 0.5-1.1 mrad.
Direction alignment with a signature, footprint or junction segment is now
judged within the band's own resolvability: sin(angle) <= RMSWidth / Span,
with the former cosine floor 1e-6 kept for degenerate widths
(`_direction_aligned`; the position-aware trace-basis rule is unchanged). The
gate is unchanged (unaligned bands must be 0); a 45-degree diagonal (0.785 rad)
or any band off every feature direction by more than its own aspect stays
flagged, and the four-edge/ten-edge 5.5-10 um bands (angles 1e-5 to 3e-4 rad,
resolvability 1e-3 to 8e-3) are unaffected. Every band records
`AlignmentAngles` (signature/footprint/junction), `DirectionResolvability`
and the verdict, and `FeatureSegments.Alignment` states the rule. Should a
resolvability-accepted band ever turn out to lie off its feature, the
position-aware alternative (both band endpoints on the signature segment
lifted to the band's plane by the metal thickness, as the trace-basis rule
does) remains available. Unit-tested: a 1 um band tilted 5 mrad accepted, the
same band at 20 mrad rejected, a 45-degree band rejected, a 10 um band at
0 accepted and at 5 mrad (above its 1.1e-3 resolvability) rejected.
Requalification under the production recipe (commits 70dc4c369 / 3fa5ee380,
adapter 72f741e3..., MMG 5.6 a97d9580...; roots
`/tmp/coupon-canonical-four-edge-9d2cb9bbb3fe-70dc4c369-20260917-233153` and
`/tmp/coupon-canonical-ten-edge-6791f1c84123-70dc4c369-20260918-001452`; the
audits/verification were regenerated on the same roots under the decision-36
producer, CanonicalBuildIds unchanged): four-edge identity + rotate-z-0.63
VERIFIED - 3,471,480 tets (seed 1,507,719; 200,348 required), minimum scaled
Jacobian 0.0200, maximum Jacobian condition 638.6, corners 3.76/3.29/3.72/3.46
== seed, 0 repairs, protected 1.04e-10, ownership closure 8.0e-13, 0 diagonal
bands, far-field pressure 1.0, design gate not applicable (layer 26,870 cells,
tangential P50 12.9 nm, transverse P90 17.1/35.4 nm; layer-adjacent band 6,022
cells, 50.0 nm vs 40.4/66.4 nm, 2.66 x NormalSize, vertices 33-71 nm from a
span), canonical 219 s / 5.08 GiB, placement 150 s / 5.28 GiB, audits 269 /
395 s, verification 1,170 s; identity SHA256 45ab1d37..., CanonicalBuildId
1b532338.... Ten-edge identity + rotate-z-0.63 VERIFIED - 3,570,533 tets
(estimate before building 3.65M; seed 2,125,258; 310,643 required; far-field
budget policy pressure 1.143, effective far 0.1828 / growth 0.571, as
designed), minimum scaled Jacobian 0.0200, maximum condition 507.7, ten
corners 3.37-3.75 == seed, 0 repairs, protected 1.50e-10, closure 6.6e-13, 0
diagonal bands (2 before decision 36), design gate not applicable (layer
54,937 cells; adjacent band 11,825 cells, 50.1 nm vs 41.3/68.2 nm, 2.73 x
NormalSize), canonical 284 s / 5.51 GiB, placement 159 s / 5.80 GiB, audits
308 / 592 s, verification 1,430 s; identity SHA256 12f485e7...,
CanonicalBuildId 548c28fe.... Layer census on both: transverse tet-edge P50
per shell 0-2/2-5/5-10/10-25/25-50 nm = 4.0/4.2/8.6/26.9/33.7 nm, 0 layer
cells below 0.02, corner balls 725-1,506 cells each.
The other ten manifest cases were run stages-only under the same driver and
none is buildable with its frozen inputs; every failure precedes the recipe:
one-edge-straight and one-edge-cad-subdivided (metric: 2 contract
SemanticCorners vs 1 boundary `Physical` vertex, "Transformed physical boundary
differs from semantic corners"), three-edge-current-calibration (metric: seed
boundary labels differ from the frozen semantic contract), two-edge-multislot
(seed: semantic corner (0, 0, 0) absent from the seed CAD), two-edge-transition
and six-edge-cluster (seed: the required-region gates fail after optimization -
corner aspect 4.11, 63 / 558 cells below 0.01, 16 / 175 above condition 1000 -
on Radius 12.5 fixtures whose contracts also disagree with their boundaries: 1
vs 5 and 4 vs 10 Physical vertices), concave-multislot, hole, rounded-strip and
opposed-layers (seed: "tangential mesh size must lie between fine and far
sizes" - their Radius 0.5 process gives FarSize 0.04 below the recipe tangential
0.05, and below the pre-34B 0.1 as well). Only the four-edge and ten-edge
contracts have SemanticCorners equal to the boundary's Physical vertices; the
fixture contracts were authored for the preflight feature-topology counts and
never passed through the seven stages. Repairing them means rewriting frozen
immutable inputs, a separate decision.

`geometry-independence-calibration-ma.json` is a separately labeled CALIBRATION
manifest (supervisor decision 22: MA/MS metal-edge-layer h-study) whose
cases re-mesh the four-edge inputs with the recipe parameters listed in their
`Calibration` blocks (V1 seed `--lc-tangent .05`; V2 additionally metric
`--far-growth 0.5`; EL4/EL1/EL4c/EL1c below); it mirrors the production tools
and gates except the labeled anisotropy-design gate `MinimumAchievedAspect 0.9`,
which can never be used by the production suite (asserted by
`test_general_mesh_manifest.py`). The per-case verifier binds each calibration
label to its build: the recorded `seed-generation` / `metric-preparation` /
`native-adaptation-mmg` commands must execute exactly every
`SeedCommandOptions` / `MetricCommandOptions` / `AdaptationCommandOptions` pair
(numerically, exactly once) and none of the `ProductionValuesBefore34B` of those
options (the production values at the time of the study), and an undeclared
production option may appear only at its pre-34B production value - the
canonical cache key does not encode recipe options, so a root labeled V2 but
built with the V1 options is rejected. The EL4c case is labeled
`AdoptedAsProductionRecipe`; the study records stay bound to the pre-34B values. `refreeze_manifest_tools.py` recomputes the
repository-tool digests of the production manifest and mirrors `Tools` /
`StageToolSHA256` into the calibration manifest in one step (`--check` reports
stale digests; runtimes, adapter and MMG library are never recomputed).

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
ten-edge basis 774 edges, 151 below its effective far size 0.1799 um (the
far-budget policy raised the requested 0.16; minimum 0.0492 um, 250 of 516).
The cut-surface rule is honoured to within ~2x on the seed cut surface
(`MaximumExtentOverRequested` 1.88 four-edge, 1.82 ten-edge; reported, not
gated) because the smallest requests (0.0217 um four-edge) fall below the
adapter `hmin` = NormalSize 0.025; the narrowest slivers are resolved at about
twice the request.
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
bound basis the behaviour is unchanged. Rigid-motion covariance of the detector
(review P2, supervisor decision 23; measured on the preserved four-edge
b8d323f5a root): (i) 681 seed-grid edges sit at exactly 2 x NormalSize = 0.05
with construction roundoff up to ~7e-12 (clusters at +6.3e-13, +4e-12,
+6.7e-12), far above the former 64-eps tolerance, so identity classified
17,609 short edges and rotate-z 17,606; the short-edge tolerance is now
`COPLANAR_TOLERANCE x ShortEdgeThreshold` (the shared dimensionless family,
5e-8 um here) and both placements classify exactly the same 16,837 interior
short edges. (ii) The "span > half the surface diameter" rule used the
axis-aligned bounding-box diagonal of the patch, which a rotation inflates
(bottom/top faces 22.63 -> 31.61 um for rotate-z-0.63 while the band span
11.44 um is invariant): the two basis bands vanished under rotation. The
diameter is now the in-plane convex-hull diameter of the patch
(`_planar_diameter`), a rigid-motion invariant that equals the former value on
every patch of the four-edge and ten-edge identity meshes (all rectangles or
corner-to-corner footprint spans), so identity records are unchanged; the
rotate-z four-edge variant now reports identity's 2/41/2/2 line-like bands
with 2 basis bands (span 22.8757 um), and the ten-edge variants stay at 0
bands (12,309 short edges both). Unit tests cover construction roundoff at the
threshold and the rotated-square band.
Supervisor decisions 27-29 (seeded transverse edge layer, 2026-09-16): the
metal-air participation is edge-singular (|E|^2 ~ r^-2/3) and NormalSize 0.025
is 25x too coarse right at the metal edges, so the seed can carry a geometric
transverse layer (`--edge-size` EdgeSize, `--edge-growth-ratio` GrowthRatio,
`--edge-layer-aspect` EdgeLayerAspect; production passes none, i.e. EdgeSize 0)
on every face bounding a metal edge (longitudinal feature curve of a metal
surface family, junction curves excluded): one embedded explicit node row per
layer of size EdgeSize x GrowthRatio^(k-1) below NormalSize at the cumulative
layer distance from the ridge, inside the ridge span on the `lc_tangent` grid
outside the corner size law. The scaled Jacobian of a tetrahedron whose corner
has three tangential edges is (hn/ht)^2, so `MinimumScaledJacobian` 0.01
(repair target 0.02) bounds the layer anisotropy to roughly 5 - nm resolution
has to come from a thin layer with nested tangential refinement, not from 50
nm-long slivers (prisms would not have this limit; not pursued): inside the
span the ridge grid is subdivided by the smallest power of two with spacing <=
Aspect x EdgeSize, row k by the nested power of two with spacing <= Aspect x
its size, the subdivision halves interval by interval towards the span ends
(the taper, no rows there), and every second row node is 5% of the offset
farther out (perfectly aligned rows form Delaunay-degenerate rectangles that
left zero-volume tets in the boundary recovery). The census records
`EdgeLayer` (sizes, rows, subdivisions, taper, corner taper offset, spans per
curve, row nodes). The metric stage binds the census layer (`--edge-size`,
`--edge-growth-ratio`, `--edge-layer-aspect` equal; LayerThickness within
`SurfaceProtectionRadius`, so the frozen band covers the layer; every span on
a band segment) and prescribes the continuous form of the same layers along
the spans: hn(r) = EdgeSize + (GrowthRatio - 1) r up to the reach (NormalSize
- EdgeSize) / (GrowthRatio - 1), then the ordinary band law; tangential size
capped at Aspect x hn blending into the band's; spans intersected after the
band segments; seed cells within the protection radius of a span are excluded
from the far-field budget policy's seed load (recorded) so the far field is
not coarsened by the transient layer seed; the adapter hmin is EdgeSize. The
restoration's CAD-correction, repair-displacement and corner-collapse bounds
are relative to the local prescribed size at each vertex capped at NormalSize
(`local_bound_size`: production bounds unchanged, EdgeSize-based inside the
layer; per-vertex bound statistics in the report) and the layer's frozen
surface vertices never move. `mesh_stage_contract.validate_edge_layer` binds
the seed command, the metric command, the census, the recipe and the adapter
hmin to one layer or none (`validate_canonical_dag`). Measured on four-edge
(V2 + layer, `edge_layer_census.py`): at EdgeSize 1 nm / Aspect 4 the layer
is surface-driven (rows 3.125/6.25/12.5/25/50 nm) and costs ~1.3M tets
(4.70M total, over the cap) with 1,955 cells below 0.02 in 1,522 repair
components; at 4 nm (rows 12.5/25/50 nm at 4/12/28 nm) 3.57M tets, layer
285k, 4 cells below 0.01 / 135 below 0.02, the metric law followed
(transverse P50 4.0/4.2/8.0/17/29 nm by shell); the 4 nm build failed the
label-restoration corner gate at (0, 8, 0) (aspect 5.39 after MMG, the
collapse rolled back and the bounded move repair cannot fix a corner cell with
three constrained vertices; identical with the previous restorer), so no
edge-layer case was qualified by that build.
Supervisor decision 30 (required tetrahedra, 2026-09-16): MMG's corner output
had blocked four campaigns, so the near-corner/near-edge region is now
deterministic - the seed defines it and MMG must not touch it. The metric
stage emits the bound artifact `required-tetrahedra.txt` (1-based seed
tetrahedron indices; `edge_volume_metric.required_tetrahedra`): every seed
cell whose centroid lies within `CornerIsotropyRadius` of a semantic corner
(also without an edge layer) and, with a layer, every seed cell with a vertex
within LayerThickness x (1 + RowZigzag) + EdgeSize of a recorded span (the
rows and the cells touching them, so the layer is kept without holes); the
recipe `RequiredTetrahedra` record carries the rule, the radius, the reach and
one count per corner and per span. The reviewed adapter takes the list as the
flag `--required-tetrahedra FILE` (removed from the positional arguments, so
every existing command shape is unchanged) and calls
`MMG3D_Set_requiredTetrahedron` plus `MMG3D_Set_requiredVertex` for each
listed cell (MMG 5.6 also tags the vertices itself, `MMG3D_set_reqBoundaries`),
failing closed on out-of-range or malformed indices; the wrapper requires the
list, checks it against the recipe record and binds it by SHA-256 in the
receipt (`RequiredTetrahedraSHA256`, count); the stage contract binds it as a
metric output and an adaptation input (`--required-tetrahedra`) and
`validate_required_region` checks the record, the list and the seed gates
below. MMG keeps the listed cells verbatim (vertex order included: the
corner-incident aspects after adaptation equal the seed census's); MMG writes
them under the binary Medit `RequiredTetrahedra` keyword (libMeshb code 12),
which meshio mis-parses, so `mesh_array_io.read_medit_binary` reads the native
output by keyword positions and flags the kept cells (`medit:required`); the
restorer freezes their vertices (no collapse, no repair move; `RequiredTetrahedra`,
`RequiredVertices`, `MaximumRequiredVertexCorrectionUm` in the report). Because
the required region's quality after adaptation is the seed's, the seed stage
satisfies the gates itself: `--maximum-corner-aspect`, `--minimum-scaled-jacobian`
and `--maximum-quality-displacement-over-normal` (the restorer's values, bound
equal by the contract) turn on `optimize_required_seed_region!` - the
restorer's bounded rule on the seed (surface vertices in the null space of
their triangle normals, three planes fix a vertex, ratio x the local
prescribed size, touched cells keep min(original, 2 x gate) scaled Jacobian),
as greedy coordinate descent (26 directions for a free vertex, 8 in a plane,
2 on a line; halving steps down to 1/128 of the bound) on the corner-incident
aspects through their 32-norm (a smooth proxy of the maximum that keeps
descending where several cells tie for the worst; target 0.95 x the gate on
the true maximum) and on components of required cells below 2 x the gate; the
census records `SeedQualityOptimization` (before/after per corner, required
minimum scaled Jacobian, moved vertices, bound usage) and the seed fails closed
when a corner exceeds the gate or a required cell stays below it. Measured on
the four-edge 4 nm seed: corners 3.41/4.63/3.42/8.69 -> 3.41/3.70/3.42/3.50,
the 88 required cells below 0.02 (needles: 4-6 nm surface triangles joined to
an interior vertex 48 nm away, minimum 0.0092) all lifted to >= 0.02 with 358
vertices moved by at most 18.75 nm (the bound), +10 s of seed time; on the
production seeds the corner balls alone are required (four-edge 2,882 cells,
corners 3.27/4.89/6.25/8.16 -> 3.27/3.39/3.67/3.56; ten-edge 8,562 cells, worst
corner 8.45 -> 3.75, all ten <= 3.86). MMG kept every required cell verbatim
(the preserved 4 nm seed: 199,572 of 199,572, 3,450,791 tets in 46 s, every
cell below 0.02 a seed cell) and the label restoration passes with 0 collapses
and no corner repair (four-edge 222 s, ten-edge 181 s, 4 nm 35 s). The adapter
build is recorded in `testdata/adapter-build.json` (command, compiler, flags,
rpath, source/executable/dylib SHA-256); its manifest digest is refrozen only
through `refreeze_manifest_tools.py --adapter-mmg PATH`, which accepts an
executable only when the record names its digest and the repository's
`adapt_edge_metric.cpp` digest. The machine's Julia launcher is refrozen the
same way (`--julia-runtime PATH`, the three Julia runtime roles together).
Supervisor decision 31 (2026-09-17): the achieved-anisotropy design gate
(TangentialP50 >= MinimumAchievedAspect x transverse P90 over the band cells)
is a statement about the metric-driven band, which a seeded edge layer
contradicts by construction (decision 28 caps the layer's tangential size at
Aspect x hn: 12.5 nm at the 4 nm rows against a 35 nm transverse P90 of the
band sample). The audit producer therefore computes the band statistics over
band cells outside a recorded edge layer (`directional_widths`: cells whose
centroid lies within LayerThickness + EdgeSize of a restoration-recipe
`EdgeLayer` span are excluded, counted as `ExcludedEdgeLayerCells`, and
reported as `AchievedAnisotropy.EdgeLayer` with their own percentiles); the
rule is keyed on the recipe's EdgeLayer record, never on a case, gate values
are unchanged (production 1.5, calibration 0.9), and the layer's design
statement remains the bound EdgeLayer aspect rule.
Decision 35 (above, "Production recipe") completes this: with the production
layer covering the whole one-NormalSize band the gate is not applicable by
construction and the layer-adjacent band is recorded, never gated.
Review of 836a1c19f (P1/P2, 2026-09-17): the seed optimizer gated the required
set computed before its vertex moves (4 nm root: 199,572 gated, 199,992 listed
by the metric stage on the moved seed), so `optimize_required_region!` now
recomputes the set on the moved positions, runs one more scaled-Jacobian pass
when the membership changed and judges the gates on the final set (census
`RequiredTetrahedra`, with `RequiredTetrahedraBeforeMoves` and
`RequiredSetRecomputations` reported); `validate_required_region` requires
census count == recipe `Count` == the label restorer's `RequiredTetrahedra`
(its `.projection.json` report; the restorer itself fails closed when the
adapted mesh carries a different count). There is one recorded layer reach,
`EdgeLayer.RequiredReach` = LayerThickness x (1 + RowZigzag) + EdgeSize
(`edge_volume_metric.EDGE_LAYER_CELL_RULE`: a tetrahedron with a vertex within
it is a layer cell), used by the required region, the restorer's frozen surface
vertices and the audits' layer exclusion alike (the audit previously used the
centroid within LayerThickness + EdgeSize). `read_medit_binary` rejects Medit
version 4 (64-bit counts) instead of mis-parsing it.
Supervisor decision 32 (2026-09-17, user directive): the 1 nm x 50 nm layer
(EdgeSize 0.001, growth 2, rows 1/3/7/15/31 nm, tangential = lc_tangent 0.05
with no nested subdivision, EdgeLayerAspect = lc_tangent / EdgeSize = 50) has
layer cells whose scaled Jacobian is (EdgeSize / lc_tangent)^2 ~ 4e-4 by
construction, so MinimumScaledJacobian 0.01 cannot judge it. The calibration
manifest, and only it, carries the layer-local quality rule
`Gates.EdgeLayerQualityRule` (`edge_volume_metric.EDGE_LAYER_QUALITY_RULE`):
inside the recorded edge layer (`EDGE_LAYER_CELL_RULE`) a cell passes when its
scaled Jacobian exceeds the roundoff floor 1e-12 (positive orientation) and its
longest edge over its shortest height (`tetrahedron_edge_aspect`) is at most
`MaximumEdgeAspect` = 2 x lc_tangent / EdgeSize = 100 (design aspect: the
ridge-grid diagonal sqrt(2) x 50 nm over the 1 nm first slab = 70.7; factor 2
for the row zigzag and the Delaunay slab split); the layer minimum scaled
Jacobian and its cells per decade are reported as diagnostics; every other gate
(MinimumScaledJacobian, MaximumJacobianCondition, corners, protected surfaces,
ownership, diagonal, resources) judges every cell outside the layer at its
production value. The rule is applied consistently: the seed optimizer
(`--edge-layer-maximum-aspect`, bounded descent on the edge aspect of layer
components above 0.95 x the bound, scaled-Jacobian passes on the other
required cells, fails closed), the label restorer (same option: layer cells are
never repair targets, gated by `edge_layer_quality`, and every layer cell must
be an MMG required cell), the audit producer (`MeshQuality.EdgeLayer` /
`MeshQuality.OutsideEdgeLayer`), `general_mesh_manifest.audit_manifest_evidence`
(`edge-layer-quality`; `mesh-quality-jacobian` on the outside statistics when
the manifest carries the rule, on the whole mesh otherwise) and the verifier
(`validate_edge_layer_quality_rule_binding`: the seed and restorer of a case
declaring `Calibration.EdgeLayerQualityRule` execute exactly the manifest bound;
every other case neither). Production is protected structurally:
`validate_manifest` refuses `Gates.EdgeLayerQualityRule` in a manifest without a
`Calibration` block and refuses any production case with a `Calibration` or
`EdgeLayer` block (preflight fails closed); `test_general_mesh_manifest` asserts
both and the gate negatives (a layer cell with negative orientation or an
aspect above the bound fails; a non-layer cell below 0.01 still fails).
Measured on the four-edge 1 nm seed (`four-edge-calib-ma-el1nm-50`, root
`/tmp/coupon-calibration-ma-el1nm-50-dd9f57120-20260917-110949`): 81,821 layer
cells with edge aspects P50 10.6 / P90 70.7 (the design value) / P99 104.5 /
maximum 1,154 before repair; the tail is Gmsh's volume split (interior vertices
0.3 nm from a row node under the face plane, 13 cells at 1,150 that no bounded
move can fix), so the seed optimizer first collapses free interior layer
vertices with an incident edge below EdgeSize onto the neighbour whose cavity
has the best worst edge aspect (9 vertices, 44 cells removed, 56 remapped;
`collapse_short_layer_edges!`) and then descends on the remaining components
with the layer cells' scaled-Jacobian floors replaced by the orientation floor
(their vertex-0 scaled Jacobian is not a quality measure and was measured to
block every aspect repair): maximum 1,154 -> 827 -> 95.0, 0 cells above 100,
layer minimum scaled Jacobian 3.0e-4 (diagnostic), +12 s of seed time. MMG kept
all 84,653 required cells (census == recipe == restoration), adaptation
3,137,539 tets in 49 s, restoration 0 components / 0 collapses in 29 s, corners
3.71/3.80/3.42/3.84, outside-layer minimum scaled Jacobian 0.0428, protected
1.0e-10, ownership closed, 0 diagonal bands, identity and rotate-z covariant;
transverse tet-edge P50 per shell 0-2/2-5/5-10/10-25 nm = 1.05/2.1/4.2/16.8 nm
at 50 nm tangential. Every physical gate passes on both variants; the
calibration design gate `achieved-anisotropy` fails on the band next to the
layer (TangentialP50 55.8 nm vs transverse P90 46/77 nm), the EL4 outcome that
decision 31 leaves to physics.
The rule is case-keyed (`general_mesh_manifest.case_gates`, used by the verifier
and the suite runner alike): a case that does not declare
`Calibration.EdgeLayerQualityRule` is never judged by it, since its seed and
restorer never executed the bound. The 4 nm layer case
(`four-edge-calib-ma-edge-layer-4nm`, 3,450,085 tets, identity SHA 22d6d926...)
therefore stays under its scaled-Jacobian-gated record ("physical gates pass
under the SJ gate", design gate pending physics, decision 31/32 rulings): it is
judged by MinimumScaledJacobian 0.01 on its whole mesh, exactly as when it was
verified, and it is not re-verified under the layer rule (its layer set has
slivers at edge aspect 761 with vertex-0 scaled Jacobian >= 0.02 that the rule
would fail and a rebuild with the seed collapse/descent would be needed to
pass). The manifest Gates, rule included, remain the canonical cache key of
every build of the manifest. The restorer computes the adapter's required flags
before the corner-ball collapse and compacts them with the cells, so the layer
rule's "every layer cell is a required cell" check stays aligned when a collapse
removes cells; the seed's Gmsh connectivity mutation
(`apply_seed_cell_collapse!`) is pinned by a tiny-box round-trip test (written
points == used points, written cells == census).
Supervisor decision 33 (2026-09-17, user directive): physics-05's per-segment MA
located the residual of control 1 (63% of its MA within 0.5 um of the outer
corner (0, -2)) and of control 23 (77% within 0.5 um of the termination
(10, 0)) in the un-layered corner regions - the 25 nm isotropic balls plus the
0.137-0.187 um taper before the layer rows. Corner grading (`CornerGrading`,
seed and metric option `--corner-size`, `mesh_stage_contract.validate_corner_grading`):
inside every frozen corner ball (CornerIsotropyRadius 0.1 around each semantic
corner) the uniform NormalSize is replaced by a geometric isotropic grading from
CornerSize at the corner point growing by the edge-layer ratio to NormalSize -
shells of size CornerSize x GrowthRatio^(k-1) ending at the cumulative radii
CornerSize (GrowthRatio^k - 1) / (GrowthRatio - 1) (4 nm: sizes 4/8/16 nm to
4/12/28 nm, NormalSize from the reach 0.021 to the radius), the edge layer's
rows with CornerSize for EdgeSize. The seed carries the shells (Gmsh MathEval
step field; the ridge nodes through a ball fall on the shell radii and a node
sits on the ball boundary; the required-region optimizer's bounds are 0.75 x
the shell size), the metric prescribes the continuous form min(NormalSize,
CornerSize + (GrowthRatio - 1) d) (`edge_volume_metric.corner_ball_size`,
which the shells never exceed) and the restorer's local size follows it; the
corner-ball cells are MMG required tetrahedra, so the adapted corners are the
seed's, gated by the seed (corner aspect <= 4, scaled Jacobian >= 0.01). With
a layer the rows start at the ridge node on the ball boundary with no taper
(census `EdgeLayer.LayerReachesCornerBall`, `UnlayeredEdgeLengthPerCorner`: per
corner the distance to the nearest span end of each layered edge, 0.1 here
against 0.187 before); a ball-boundary crossing within half NormalSize of a
CAD vertex is not a node (a vertical corner edge of the metal thickness 0.1
ends exactly on the radius; measured: a node at 0.9 nm from the top corner
made six flat cells). `--corner-size 0` (absent) is the production ball
(uniform NormalSize), which the EL4 and EL1 records keep. Seed census per
corner: `Shells` (edges by midpoint distance, P50/P90 against the shell size,
cells by centroid); `edge_layer_census.py` reports the same on the final mesh
(`CornerBalls`). Case `four-edge-calib-ma-el4c` = EL4 + CornerSize 0.004:
seed 1,507,725 tets (EL4 1,508,461), 200,354 required (EL4 199,992), shell
edge P50 per corner 3.8-5.8 / 6.0-7.8 / 12.1-17.7 / 27.2-29.5 nm against
4/8/16/25 nm, corners 3.76/3.29/3.72/3.23, un-layered length 0.1 at every
layered edge. Case `four-edge-calib-ma-el1c` = EL4c with EdgeSize 0.001
(rows 1/3/7/15/31 nm, aspect-4 nested rows 3.125-50 nm) at EL4c's corner
shells (CornerSize 0.004, unchanged, so EL4c -> EL1c changes exactly the edge
layer; CornerSize 0.001 failed the seed's gates: corner aspect 4.81 from a
0.41 nm Gmsh surface edge next to the 1 nm shell node, one layer/ball-junction
cell at scaled Jacobian 0.0073 - the recorded P2 test cases for a corner-ball
collapse) declares the labeled
calibration-only element cap `Calibration.MaximumElements` 5,000,000
(supervisor decision 33): `general_mesh_manifest.validate_case_element_cap`
requires the manifest's `Calibration.GateDeviations.MaximumElements` to name
exactly the declaring cases with the production value 4,000,000, and
`case_gates` judges only those cases by the cap; `Gates.MaximumElements` stays
4,000,000 in both manifests (the build cache key, and the gate of every other
case), no production case may declare a cap, and the stage bounds are unchanged.
Measured outcomes of decision 33 (roots
`/tmp/coupon-calibration-ma-el4c-4f9946870-20260917-155730` and
`/tmp/coupon-calibration-ma-el1c-cb3f64d81-20260917-160517`, adapter 72f741e3...,
MMG 5.6 a97d9580...): EL4c PASSES every physical gate on identity and
rotate-z-0.63 - 3,471,507 tets (EL4 3,450,085), minimum scaled Jacobian 0.0200,
maximum Jacobian condition 898.6, corners 3.76/3.29/3.72/3.23 == seed, 0
restoration repairs, protected 1.0e-10, ownership closed, identity SHA256
5fb3a5c8..., CanonicalBuildSHA256 819f50cf... - and fails only the calibration
design gate `achieved-anisotropy` on the band next to the layer (TangentialP50
50.0 nm vs transverse P90 40.5/66.5 nm; decision 31: physics adjudicates;
physics-06 ran on it). EL1c as built at cb3f64d81 FAILED the physical gate
`mesh-quality-jacobian`: MaximumJacobianCondition 5115.5 > 1000 on 56
edge-interior layer cells on the metal faces (209 > 500), each with a 0.02-0.24 nm
edge between a face vertex and a 1 nm row node (minimum scaled Jacobian 0.01696
passes; 4,479,207 tets < the 5,000,000 cap; corners 3.55/3.24/3.42/3.41 == seed;
identity SHA256 9e9e1567...); its identity variant was audited, normalized and
verified, its rotate-z-0.63 variant audit producer exceeded the 1800 s audit
bound twice (once alone: the rotated topology/protected-surface pass takes
20 min on 4.48 M tets) - recorded as "identity verified; rotate-z audit exceeded
the audit bound", no bound change; not run in physics. Both cases share
`CanonicalBuildId` 40d6db94... (the cache key is source + process + contract +
recipe + manifest `Gates` + tools; neither the case `Calibration` options nor
the case-level element cap enter it), so `run_general_mesh_suite` over the whole
calibration manifest would refuse the second case ("shared canonical stages
require exact cache key and hashes"): calibration cases are verified per case
only (`verify_canonical_case_entries.py`), a documented limitation.
Decision 34 (2026-09-17) located the EL1c failure on the seed itself: the 56
cells (and the CornerSize 1 nm probe's 0.41 nm corner edge) are made by the
seed's bounded descent, which moved face vertices to 0.013-0.24 nm of a row or
ridge node while repairing scaled Jacobians and corners - the Gmsh seed has no
edge below half EdgeSize (the sub-size collapse finds 0 candidates below the
threshold). The seed now (a) collapses sub-size edges (`collapse_short_edges!`,
threshold `SEED_COLLAPSE_SIZE_FRACTION` 0.5 x the smallest prescribed size:
EdgeSize or CornerSize; a threshold at the size collapsed 14,667 legitimate layer
face vertices between the 1 and 3 nm rows and 2,677 row nodes) for interior and
surface vertices alike, a surface vertex only along its own surface (onto a
surface neighbour in every plane of the vertex carrying every line/triangle
support of the vertex: a face vertex within its face, a row or ridge node along
its curve, CAD points never; remapped triangles keep their normal; the cavity is
no worse in maximum edge aspect and Jacobian condition and better in one, and
no cell falls below the scaled-Jacobian gate where none was; census
`SubSizeEdgeCollapse` with every collapse's position and quality), with and
without the layer quality rule; (b) guards the descent so that every touched
cell keeps max(original, 0.95 x MaximumJacobianCondition) Jacobian condition and
a moved vertex in the collapse region keeps its edges at or above min(original,
threshold); (c) gates the scaled-Jacobian-gated required cells by
`--maximum-jacobian-condition` (the manifest MaximumJacobianCondition 1000,
carried equally by the seed and the restorer, `REQUIRED_REGION_GATE_OPTIONS`;
census `RequiredMaximumJacobianCondition{Before,After}`,
`RequiredCellsAboveConditionAfter`; restorer `RequiredMaximumJacobianCondition`),
failing closed at the seed instead of after the hour-long build. EL1c seed under
(a)-(c): 2,147,264 tets (unchanged connectivity), 834,434 required, corners
3.55/3.24/3.42/3.41, required maximum condition 221.7 before the moves and 924.2
after (0 above 1000; without the guards 8,525 on 43 cells), minimum scaled
Jacobian 0.0163, 188 s. Required region of the seed, old -> new build: maximum
condition 5115.5 -> 924.2, cells above 1000/500/200: 56/209/983 -> 0/106/741,
shortest edge 0.024 nm -> 0.500 nm (1,228 -> 0 edges below 0.5 nm), minimum
scaled Jacobian 0.01696 -> 0.01630 (one cell in [0.01, 0.02) in both, the
layer/ball-junction cell of the amendment).
EL1c rebuilt under decision 34 (root
`/tmp/coupon-calibration-ma-el1c-cb3f64d81-20260917-205035`, tools of commit
3c6b4c7fe, adapter 72f741e3..., MMG a97d9580...): 4,483,816 tets (< the
5,000,000 cap; 890,977 H1 DOFs), MMG kept all 834,434 required tetrahedra
(recipe == receipt == restorer), restoration 1 repair component rejected (the
frozen junction cell, 0.0163 >= 0.01), 0 corner-ball collapses. Identity
variant: every PHYSICAL gate passes - minimum scaled Jacobian 0.01630, maximum
Jacobian condition 924.2 (outside the layer 221.7), corners 3.55/3.24/3.42/3.41
== seed, protected surfaces 1.04e-10, ownership closure 3.1e-13 (0 unmatched, 0
overlapping), 0 diagonal bands, canonical build 444.8 s / 7.47 GiB and
placement 201.1 s / 7.39 GiB within 1800 s / 8 GiB (the publication stages now
peak at 7.4-7.5 GiB, 0.5 GiB under the bound) - and the calibration design gate
`achieved-anisotropy` fails on the layer-adjacent band (TangentialP50 41.2 nm vs
transverse P90 35.4/60.6 nm), as for EL4/EL4c (decision 31). Layer census:
transverse tet-edge P50 per shell 0-2/2-5/5-10/10-25/25-50 nm =
1.05/3.5/13.4/20.0/31.2 nm at tangential 3.1/3.5/7.7/15.4/32.8 nm; corner balls
779/1156/773/1522 cells, 0 below 0.02. The rotate-z-0.63 variant audit producer
again exceeded the 1800 s audit bound alone (bounded-run, topology/quality,
complexity and invariants written at +24 min, the variant-transform record not
reached; its written records equal the identity's: scaled Jacobian 0.01630,
condition 924.2, corners, protected 1.04e-10, ownership closed), so EL1c is
recorded as "identity verified on every physical gate; rotate-z audit exceeded
the 1800 s audit bound" - no bound change (`per-entry-verification.json`:
identity failures = [achieved-anisotropy] only; rotate-z evidence missing).
identity.msh == canonical.msh SHA256 2a3de0d9..., rotate-z-0.63.msh 64bd3c3b...,
CanonicalBuildId 1742cc0c... (new tool digests), CanonicalBuildSHA256 30854ab8....
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
digest. The bounded record consumes seven separately executed reports forming a
digest-linked DAG: canonical source validation, seed generation, metric
preparation, native adapter/MMG, label restoration, final Gmsh publication, and
proper rigid publication. A tool is accepted only when its
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
seven `--stage-report STAGE=PATH` bindings to `bounded-run` and
`mesh-topology-quality`, or run its `variant-audits` kind once per variant: one
bounded process (the same `run_bounded_mesher.py` limits, 1800 s / 16 GiB, with
the five records declared as `--artifact`s) reads the mesh once and writes the
five records `<prefix><kind>.json` through the same per-kind producer functions,
so every record equals its standalone twin except the recorded `Command`
(unit-tested); on the four-edge identity variant this replaced 434 s of five
separate processes (each re-reading the 137 MB mesh; peak 3.1 GiB) by one
226 s run (peak 3.0 GiB).
Normalize only those five audit outputs:

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
