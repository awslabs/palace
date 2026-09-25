# Geometry identification: audit tooling and baseline

Tooling for the geometry-identification block of the surface-response (SPR) correction
(decisions 72 / 73 in `coupon-accuracy-assessment-20260913/SUPERVISOR-DECISIONS.md`, plan in
`GEOMETRY-IDENTIFICATION-PLAN.md`). Everything here is identification only: geometry-only
preflights (`palace --surface-response-preflight`), the audit, layout oracles and synthetic
layouts. No field solves, no library builds. The classifier (`palace/`) is not modified by
this block: every disagreement below is a recorded defect for the fix block.

| Module | Purpose |
|---|---|
| `msh2.py` | binary / ASCII MSH 2.2 reader (first- and second-order simplices) |
| `perimeter.py` | metal perimeter E recomputed from the mesh with the classifier's definitions (metal attribute union, face-direction classes on the distinct geometric faces, truncation, 30 deg corner rule on the 1e-12 direction grid, physical chains, edge-connected components) plus the excluded classes of decision 73(3): NONPLANAR / FOLD, CROSS_LAYER zones within 2R of off-plane metal, NONMANIFOLD, EMBEDDED, BOX; `rounded_runs` reads fillet arcs with the classifier's rounded-corner rule |
| `manifest.py` | canonical digest / set diff / aggregate summary of `surface-response-requirements.json` (drops Status, SelectedModels, NormalizedLibraryDistance, Reason, Library.Path, Statistics); preflight log parser |
| `audit.py` | gates A1 (perimeter agreement, length and count partition, multiplicity proxy, weights, vertex census, recorded exclusions) and A2 (cluster balls), B1 gap bound; JSON + Markdown; exit 1 on any failure; `--compare` for A3 / A5 set identity |
| `preflight_config.py` | generic geometry-only preflight configuration (Electrostatic; Ground + Terminal metal for conductor identity; typed SA / MS / MA interfaces with the library radius) |
| `preflight_matrix.py` | ranks x libraries x UniformLevels x CrackInternalBoundaryElements cells with the audit and cross-cell digests (A3 / A4 / A5) |
| `synthetic_layouts.py` + `synthetic_layouts.jl` | 62 synthetic stress layouts with exact oracles (Gmsh / OCC via Julia), runner and oracle comparison (A6) |
| `tag_metal_components.py` | N-island conductor tagger for production meshes; repairs: drop metal triangles duplicating port faces, drop attributes, add a material-interface (substrate_air) group |
| `refine_msh2.py` | uniform 1 -> 8 / 1 -> 4 refinement of a first-order MSH 2.2 mesh outside Palace, so the refined mesh exists for the audit (real A5 test) |

Tests (`python3 -m unittest discover -s examples/surface_response_identification -p 'test_*.py' -t examples` from the repository
root, or per module): `test_audit.py` (18, incl. 3 version-2 gate tests), `test_preflight_config.py` (3),
`test_synthetic_layouts.py` (15), `test_tag_metal_components.py` (5), `test_refine_msh2.py` (1).

Typical use:

```bash
cd examples
python3 -m surface_response_identification.preflight_matrix --mesh M.msh2 --ground 5 6 7 \
    --terminal 9 --sa 8 --library seed=transmon/benchmark/transmon_surface_process_seed.json \
    --ranks 1 2 4 6 --uniform-levels 0 1 --crack true false --output OUT
python3 -m surface_response_identification.synthetic_layouts --output OUT/synthetic --ranks 1 2 4 --uniform-levels 0 1
python3 -m surface_response_identification.audit --mesh M.msh2 --config cfg.json \
    --manifest postpro/surface-response-requirements.json --log palace.log --output-prefix OUT/audit
```

## Manifest version 2 (identification fix block, phase 1)

The classifier now runs the pure-geometry identification of
`palace/models/SURFACE-RESPONSE-IDENTIFICATION.md` (`palace/models/surfaceresponseidentification.cpp`)
in the preflight: `surface-response-requirements.json` carries `Version: 2`, the version-1
`Requirements` derived from the features (aggregated by signature, `Hash` per record, `Status`
from the key-based matching pass), the legacy per-pass records under `LegacyRequirements`
(comparison only) and the contract under `Identification` (`Features` with canonical
`Signature` / `Hash` / `Chirality` / `Portions` / `Vertices` / `Frame` / `Match`, `Segments` with
`Key` and `Portions` or `Exclusion`, `Vertices`, `Exclusions`, `Totals`, `GeometryDigest`,
`Conventions`). The audit reads the version-2 contract when present (exact partition, vertex
census, exclusions by class, cluster disjointness, `GeometryDigest` + per-signature lengths
for A3 / A5) and falls back to the aggregate version-1 reading otherwise; `preflight_matrix`
and `synthetic_layouts` compare the `GeometryDigest` across cells. A library model may carry
the feature's `Signature` object (`"Signature": {...}`) and is then matched by key. The
synthetic oracle follows the design rules (sampled event cores, corners join a cluster when a
core lies within 3R, parallel pairs survive outside cluster regions and corner windows).
Phase-1 results: `coupon-accuracy-assessment-20260913/geometry-identification-fix-20260924/phase1/REPORT.md`.

Phase 2 (curved-edge chain rule, design item (b) 7): signatures are rounded on the recorded
`SignatureLengthQuantumOverR` = 1e-6 / `SignatureAngleQuantumDegrees` = 1e-6 grids; chains are
polylines whose sub-corner joints carry a windowed curvature (window R, turns spread over the
adjacent half-chords); portions with a windowed bend radius below `StraightBendRadiusOverR` = 10
(decision 75; phase 2 used 20) are the curved classes `CurvedEdge` / `CurvedSameConductorStrip` / `CurvedSameConductorGap` /
`CurvedDifferentConductorGap` (`RadiusOverR` in the signature), straight-like portions keep the
straight classes with a `BendRadiusOverR` annotation on every feature; two chains pair along a bend
when their closest-point separation is constant within `PairSeparationToleranceRelative` = 0.05
(concentric arcs and CPW gaps along bends are pairs, never clusters). The synthetic oracle applies
the same rules (`design_bent_pairs`; arc-bar layouts carry their design bend for the expected
classes; checks `A6-parallel-pairs` with the pair tolerance and `A6-bent-pair-classes`).
Phase-2 results: `coupon-accuracy-assessment-20260913/geometry-identification-fix-20260924/phase2/REPORT.md`.

Phase 3 (decision 74 step 4: conductor identity, crack-independence, decision-73(3) exclusions,
embedded sheets, mesh preparation):

* **Perimeter from the distinct geometric metal faces** (`palace/utils/metaledge.cpp`): a face
  edge is classified by the in-plane inward directions of the distinct faces owning it — one
  direction = one-sided perimeter (PHYSICAL / TRUNCATION), two opposite = interior, two
  non-coplanar = FOLD, three or more = NONMANIFOLD. Coincident crack copies count once, so
  `CrackInternalBoundaryElements` true and false give the same perimeter and identification
  (same `GeometryDigest` and lengths), and sheets with one material on both sides (airbridge
  spans, embedded metal) are no longer cancelled. `perimeter.py` applies the same rule
  (kinds PHYSICAL, TRUNCATION, EMBEDDED, NONPLANAR, FOLD, NONMANIFOLD, BOX).
* **Conductor identity = edge-connected metal component**, independent of the attribute
  numbering and of the problem type (electrostatic Terminal / Ground labels and Maxwell PEC
  attributes alike); the labels only produce a warning when one component carries distinct
  labels. Consequences: two ground planes cut by the simulation box are two conductors (a CPW
  cross-section is a three-conductor four-edge cluster); two sheets touching at a single vertex
  are two conductors (no galvanic connection) and the vertex is a `PointContact` record + warning;
  the `Junction` signature carries `ArmConductors`. The audit's `perimeter.py` reports the
  edge-connected component of every edge; `tag_metal_components.py` is no longer needed for
  conductor identity (it remains a mesh-repair tool).
* **Exclusions reported, never silent** (`Identification.Exclusions`, each with count and
  length; the segment table carries `Exclusion` for whole segments and `ExcludedPortions`
  `[s0, s1, exclusion index]` for analytic zones): `SimulationBoundary` (metal faces on the
  mesh bounding box, a PEC box), `NonPlanar` (folds; one-sided faces not parallel to the layer
  normal: walls, staple legs, vias), `NonManifold` (a wall standing on a sheet),
  `UndeterminedProcessSide` (same material on both sides and no `EdgeFrameNormal` on the target
  interface — the configured process layers decide the side, else the record), `CrossLayer`
  (the parts of a planar run within 2R of metal off its own plane — a facing layer across a
  gap, a wall, a staple — solved analytically on the distance to every such face, so the zones
  are refinement invariant; `Conventions.CrossLayerReachOverR` = 2), `Untargeted`,
  `TruncationCut`. Vertices: `Excluded` (every incident run excluded or within 2R of off-plane
  metal), `ExclusionCut` (a chain cut by an excluded segment; no endpoint feature), `PointContact`
  flag. Every other metal plane is identified in its own right (a second chip is not excluded
  unless it faces metal within 2R). The legacy classification omits exactly the segments the
  identification excludes.
* **Layer normal** = area-weighted principal direction of the metal face normals (faces on the
  mesh bounding box do not vote); planarity tolerance = the parallelism tolerance 1e-8.
* **Mesh preparation**: coincident metal boundary elements of one attribute are tolerated
  (deduplicated); coincident boundary elements of different metal attributes abort naming both
  attributes and the face location; Palace already refuses a face with two boundary elements at
  load (`geodata.cpp GetFaceToBdrElementMap`, attributes named) — `tag_metal_components.py
  --drop-metal-duplicates` is the deterministic repair. A missing `substrate_air` group is not a
  blocker (MS / MA targets only).
* **Audit**: split-tolerant perimeter agreement (cracking bisects elements next to under-resolved
  sheets, second-order edges are sampled), partition with `ExcludedPortions`, vertex census with
  excluded vertices and point contacts, exclusion lengths per class compared with the mesh
  (incl. the analytic CrossLayer zones); `synthetic_layouts.py --crack-false` adds a
  `CrackInternalBoundaryElements=false` cell per library and compares its digest; the oracle's
  conductor is the sheet (connectivity), `gap-same-*` are slots in one U-shaped sheet, and the
  off-plane exclusions (facing sheets, walls) are computed analytically (`A6-excluded-classes`).
* `StraightBendRadiusOverR` = 10 (decision 75; phase 2 used 20).
* **Pairs along bends decided on the curve separation** (`PairSeparationEstimate`): the
  separation of two constant-separation chains is the smaller of the two directional maxima of
  the sampled closest-point distance (an inscribed polyline meets its curve at the vertices; an
  exact offset polyline keeps corresponding chords at the design separation), the pair
  interacts iff that separation is below 2R on the quantized grid — the straight-pair answer at
  every discretisation — the candidate facing region is 2R (1 + 0.05) (`PairCandidateReachOverR`),
  and the cross-chord interactions of a constant-separation pair are never event cores whether
  or not it interacts. Found on DS-SCT-001: its 4 um CPW gaps (= 2R) dipped to 3.9999 mid-chord
  along the 250 um bends and formed three clusters of 874-1,184 edges (1.5-3.2 mm) claiming
  every strip. Synthetic `gap-bend-r{50,250}-{2R,2Rminus,2Rplus}-step{1,5,15}` (18 layouts: two
  concentric 8 um bars with a gap of 2R and 2R +/- 1e-3 R): oracle = the straight-pair answer
  (2R and 2R + 1e-3 R: isolated edges, no cluster, no pair; 2R - 1e-3 R: one
  `DifferentConductorGap` and the two corner pairs across the gap as clusters); the oracle's
  corner pairs are strictly within 2R (the classifier's quantized decision); `arc_bar(...,
  centre_y=...)` builds concentric bars.
* **Legacy pair construction** (`surfaceresponseoperator.cpp` `VerifyParallelOverlap`): the
  parallel-overlap verification tolerance follows the parallel class (cosine deficit 1e-8 =
  sqrt(2e-8) rad times the projected lengths) instead of 1e-10 relative; DS-SCT-001 aborted on
  190 CPW chord pairs 3e-5 to 1e-4 rad apart.
* **Mesh preparation, DS-SCT-001**: `tag_metal_components.py --drop-metal-duplicates` (16 port
  triangles duplicating metal faces; ground attr 4 + island attr 9), then `preflight_matrix
  --frame-normal 0 0 1` because the L1 metal has the substrate volume (attr 1, z in [-525, 5] um)
  on both sides (without it every sheet segment is an `UndeterminedProcessSide` record).
* **Survey findings recorded (PENDING)**: `CanonicalClusterSignature` serialises a cluster once
  per distinct portion direction (~290 s on the 1,184-edge DS-SCT-001 cluster before the fix;
  chip-scale clusters need a bounded frame rule); DS-SCT-001 A1 vertex census: the audit reads
  29 rounded runs on the BSpline chords, the manifest 18 `RoundedCorner` vertices (the audit's
  rounded-run reading and the classifier's fillet rule disagree on 11 sub-R chord runs with
  total turns 11-80 deg).
Phase-3 results: `coupon-accuracy-assessment-20260913/geometry-identification-fix-20260924/phase3/REPORT.md`;
block summary: `coupon-accuracy-assessment-20260913/geometry-identification-fix-20260924/SUMMARY.md`.

## Geometry identification: baseline audit (2026-09-24, executable 9ef5256b / v0.17.0-572-g5876402f7)

Evidence: `coupon-accuracy-assessment-20260913/geometry-identification-baseline-20260924/`
(`transmon/matrix`, `synthetic/`, `synthetic-nosa/`, `survey/chain2`, `survey/sct001`,
`REPORT.md`). Libraries: `seed` = `transmon/benchmark/transmon_surface_process_seed.json`
(Models []), `isolated` = `transmon/benchmark/corrected-p1/library/process-library.json`
(isolated-edge + 2 um strip), `full` = the 11-model device library JSON of the
classification-fix block (preflight-only: no matrices are read by the preflight).

### Transmon (island-tagged coarse mesh, 3,146 perimeter segments, 35.30 mm targeted length)

| Invariant | Result |
|---|---|
| A4 rank determinism (1, 2, 4, 6) | PASS for every library and level: identical canonical digests |
| A3 library independence (seed vs full, set diff) | FAIL, two components. (i) Manifest-format artefact: IsolatedEdge 1 record / 6 um (seed: "No compatible isolated-edge model ... correction is disabled", one representative segment) vs 3,007 / 34,853 um (full, Exact) — the geometry is the same, the record is not. (ii) Genuine library dependence: SameConductorStrip (2 um) 170 segments / 668 um (seed) vs 64 / 250 um (full); nonparallel omissions 12 vs 2; clusters: seed has 3-edge (12 um), 4-edge (16 um), 4-edge (16 um) records, full has 3-edge (8 um) and 4-edge (12 um) — the 2 cluster models matched by the full library absorb 106 strip segments (418 um), 10 nonparallel pairs and one 4-edge cluster that the no-model run reports separately. The 64 strips at exactly R are common to both. |
| A5 refinement, Palace UniformLevels 1 (digest only) | FAIL: seed 6 cluster records replaced by 7 (+1 strip record at separation 1.0 um); full: ConvexCorner 12 -> 13, 7 clusters replaced by 5, the IsolatedEdge record removed |
| A5 refinement, external 1 -> 8 split (`refine_msh2.py`, 36.8 MB, audited against the refined mesh; `transmon/matrix-r1`) | FAIL, and the digest-only reading hid the meaning: with the full library the refined mesh logs `Nearby three-dimensional metal edges are not parallel; correction is disabled for this interface group!` — Exact 37 (corners + 2 clusters, 160 um) / Missing 131, Matched physical edge segments 0, i.e. the whole 34.9 mm isolated-edge correction is switched off by one refinement level (coarse: 3,043 Exact). Perimeter agreement is exact (6,292 = 6,292 segments, 48 = 48 chains, no bisection); nonparallel omissions 2 -> 18, cross-interface 14 -> 24; ConvexCorner 12 -> 13 (a corner freed from a vanished 4-edge cluster). External and Palace-internal refinement agree on every non-cluster record; their cluster records differ (mesh-order-dependent representative events). |
| A1 length partition | FAIL: full library assigns 35,103 of 35,297 um (deficit 193.6 um = 0.55 %); seed 674 um only (correction disabled: no isolated-edge model) |
| A1 count partition | FAIL: 3,071 assigned + 17 omitted (2 + 12 cross-interface, 2 nonparallel) = 3,088 reconciled; omitted != 0 |
| A1 vertex census | NOT-EVALUABLE: 48 audit corners, 34 manifest corners; 14 absorbed by 5-6 spatial clusters or dropped (the manifest does not enumerate them) |
| A1 exclusions recorded | FAIL: 160 um non-planar + 256 um cross-layer (plane z = 10 um, 1,280 um^2) + 80 um non-manifold metal in attribute 5 are silently absent (classifier: 0 incompatible-process-normal / 0 process-normal-offset omissions) |
| A2 cluster balls | PASS (audit-side vertex clusters; minimum centre distance 6.0 >= 2R) |
| CrackInternalBoundaryElements false | ABORT (exit 134) in all 16 transmon cells, all 4 chain2 cells and all 4 DS-SCT-001 cells: `Verification failed: (norm_squared > 1.0e-20) is false: --> Unable to infer the metal-to-gap direction for an automatic edge segment! ... in function: BuildMetalEdgeGapDirections ... palace/utils/metaledge.cpp:1208` — the uncracked mesh has both materials' faces on one boundary element, so the metal-side face sum cancels; the classifier only works on cracked meshes |
| Knife edges in the layout | 180 parallel pairs at exactly R = 2 um, 4 at exactly 2R (mesh census); the 64 same-conductor-strip segments at exactly R (250 um) have no model in the full library |

### Synthetic stress layouts (62 layouts, ranks 1 / 2 / 4, UniformLevels 0 / 1, seed and isolated libraries)

A4 holds on all 62 layouts (identical digests over ranks). Oracle agreement (A6): mesh perimeter
length and corner list agree with the polygon oracle on all 62; the manifest agrees with the
oracle's expected classes on 60 (exceptions below). Observed classifier rules:

* Corner rule: a vertex is a corner iff its turn exceeds 30 deg (interior angle < 150 deg);
  `corner-150` (turn exactly 30 deg) is straight: 4 corners, `corner-135` 6. Concave corners
  report the gap-side angle (`ConcaveCorner@30` for the 30 deg bent bar).
* Interaction threshold: pairs at 3.9 um are SameConductorGap / DifferentConductorGap /
  SameConductorStrip; pairs at exactly 2R = 4.0 um and at 4.1 um are isolated (quantized
  decision: exactly 2R is not within). Same/different conductor gaps are distinguished only
  through Terminal attributes.
* Separation exactly R = 2 um (`strip-2`): the strip's four 90 deg corners are reported as
  ConvexCorner (4) while at 1.95 and 2.05 um they are absorbed into the spatial cluster; the
  transmon's 64 strips at exactly R are the same knife edge.
* Every corner within 2R of another feature becomes a spatial cluster: gaps <= 3.9 um make one
  6-edge cluster spanning the whole gap and omit 16-24 um (15-23 % of the perimeter) as
  "nonparallel"; T with a 3 um stem: 4-edge cluster, 20 % omitted; X with 3 um arms: 2 clusters,
  44 % omitted; 3 x 3 um island and 1 / 3 um apertures: one 4-edge cluster, 50 % of the length
  not assigned; 10 um aperture: 4 concave corners, complete.
* Acute corners: at 30 and 45 deg the arms' segments near the corner are within 2R of each
  other and form a cluster in addition to the corner record (double description: corner-30 has
  1 concave + 3 convex corners + clusters of 2 and 4 edges, 36 % omitted); at 60 deg and above
  none. Taper 80 deg corners become 2-edge clusters (the 100 deg ones stay corners).
* Curved edges: a fillet run (consecutive sub-tolerance turns) is a rounded corner iff its
  tangent distances from the virtual corner are < R and equal within 5 % and the fillet radius
  is in (0, R): the 0.5 um fillets of the 8 x 6 um island / aperture are
  `ConvexCorner@90 radius 0.5` / `ConcaveCorner@90 radius 0.5`. Polyline arcs of radius 5-250
  um (1-20 deg per vertex) never produce corners; the two edges of a 3 um wide arc bar are
  strips per parallel chord pair and the adjacent non-parallel chords within 2R are omitted:
  7-39 % of the perimeter (r5: 27 / 30 / 16 %, r20: 39 / 38 / 8 %, r50: 39 / 12 / 6 %,
  r250: 17 / 8 % for 1 / 5 / 20 deg per vertex) — the omission grows as the discretisation
  gets finer, so the "straight edge" reading of decision 73(1)(b) is not what the classifier
  does today for bend radii up to 125 R.
* A5 (one uniform refinement) holds only on layouts without clusters and without fillets:
  every spatial-cluster record changes under refinement (Edges intervals / points follow the
  mesh segments), and the rounded corners of the 8 x 6 um island DISAPPEAR after refinement
  (the collinear midpoints inserted on the chords split the curved run: `MatchRoundedRun`
  requires consecutive curved vertices).
* A3 (seed vs isolated): the requirement set differs on every layout with an isolated edge:
  without an isolated-edge model the classifier records one representative segment (Count 1,
  the segment's length) and disables the group.
* Excluded classes: the facing sheet 3 um above the process plane and the vertical wall leave
  the manifest identical to the bare rectangle (no omission counter, no record): silent.
* Pair records carry one edge length (TotalEdgeLength = Count x segment length): in
  `gap-same-2` 12 pairs / 12 um plus 72 um isolated of 104 um with 20 um omitted, so the
  partner side of a pair is counted as isolated (or double-counted); without a per-segment
  assignment in the manifest the audit cannot tell (A1-multiplicity NOT-EVALUABLE).

### Survey geometries (local)

* Two-transmon chain (`DeviceLayout.jl/examples/SingleTransmon/transmon_chain_2.msh2`, 413k
  nodes, second order): N-island tagger finds 2 islands (139 triangles each); 588 MB RSS at 1
  rank, 12-17 s; A4 holds at 1 vs 6 ranks; crack=false aborts; strip set depends on the
  library (340 / 1,337 um seed vs 308 / 613 um isolated); 44 omitted segments (80 nonparallel
  + missing models); 88 audit vertices vs 52 manifest corners; DifferentConductorGap at 1.0 um
  (2 segments) appears with the two Terminal islands. No substrate_air group: MS / MA only.
* DS-SCT-001 bump-Q subgraph (`mre_port_bug_fixed.msh2`): Palace rejects the mesh (16 port
  triangles duplicate metal faces: `geodata.cpp:3707 A non-periodic face cannot have multiple
  boundary elements`); after `--drop-metal-duplicates` the metal splits into ground + 1 island
  (313 triangles). The classifier then finds 207 segments (the island) of the audit's 7,679:
  the whole ground perimeter (27.1 mm, 48 corners, 939 parallel pairs at exactly 2.0 um and
  768 at exactly 4.0 um) is invisible; with the ground alone it aborts "found no metal
  perimeter". Root cause in the mesh: every L1 metal face has material attribute 1 on both
  sides (the volume attributes put the sheet inside one material), so the cracked sheet's two
  copies fall in one side component and cancel in the odd-incidence perimeter — the same
  mechanism silently drops the transmon's airbridge (80 faces with vacuum on both sides).
  Adding a substrate_air group (`--add-interface-group`) does not help (the interface is not
  adjacent to the metal); crack=false gives the same abort.
* Blockers hit: conductor identity (single PEC attribute -> conductor 0; N-island tagging
  works when islands touch no ground-adjacent attribute), duplicate boundary elements,
  metal embedded in one material, missing SA group (harmless for MS / MA targets: the
  `synthetic-nosa/` cells reproduce the SA cells' requirement sets), second-order meshes
  (classifier bisects: 12,113 segments vs 7,143 corner-edge segments).

### Consolidated findings for the fix block (ranked)

1. Canonical contract: the manifest aggregates by geometry and carries no per-segment
   assignment, no feature positions, no cluster membership, no excluded-class record; the
   requirement set depends on the library (isolated-edge representative, strip records with
   / without a strip model); cluster records depend on the mesh (representative event / tie
   break); pair records cover one side. Gates A1-multiplicity / vertex census are therefore
   NOT-EVALUABLE by construction (decision 73(4) asks for gates).
2. Whole-geometry sites: clusters absorb every corner within 2R and omit the neighbouring
   segments as "nonparallel" (15-50 % of small features); acute corners are described twice
   (corner + cluster); knife-edge behaviour at exactly R (strip-2 corners kept, 1.95 / 2.05
   absorbed).
3. Curved-edge rule: fillets with radius < R are rounded corners but vanish under uniform
   refinement (run continuity); arcs with bend radius >= R are chords whose non-parallel
   neighbours are omitted (7-39 %): decision 73(1)(b) needs a chain-level (accumulated turn)
   rule and a recorded bend-radius threshold.
4. Conductor identity: PEC -> conductor 0 for every attribute; the tagger is a workaround;
   `pec_attribute_conductors` plumbing or component labelling in the classifier.
5. Perimeter extraction on sheets with one material on both sides (airbridges, embedded
   metal, TSV walls) cancels silently; excluded classes (non-planar, cross-layer,
   non-manifold) are dropped without a record (decision 73(3) requires a reported gap).
6. CrackInternalBoundaryElements=false aborts on every real mesh (`metaledge.cpp:1208` gap direction inference); one uniform refinement disables the whole isolated-edge correction on the transmon ("not parallel ... correction is disabled for this interface group").
7. Duplicate boundary elements / missing SA group are mesh-preparation blockers (tools here).
