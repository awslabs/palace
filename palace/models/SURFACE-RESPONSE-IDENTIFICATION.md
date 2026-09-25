# Surface-response geometry identification: the canonical contract (design, 2026-09-24)

Scope: the geometry-identification step of the fabrication-process surface-response
correction (decisions 69-74 of `coupon-accuracy-assessment-20260913/SUPERVISOR-DECISIONS.md`,
invariants A1-A6 of `GEOMETRY-IDENTIFICATION-PLAN.md`). This document is the contract that the
implementation in `surfaceresponseidentification.{hpp,cpp}` follows, that the preflight
manifest (`surface-response-requirements.json`, version 2) serialises, and that the audit tool
(`examples/surface_response_identification/audit.py`) gates. Phase 1 implements (a), (b), (c),
(d) below for the manifest and the library-matching pass; the solve-path patch construction
consuming the same feature list is the first item of phase 2 unless it fits in phase 1.

## (a) The contract

Identification is a **pure function of (mesh, process parameters, R)** — never of the library.
Its output is:

1. **Segment table.** Every metal perimeter segment of the mesh (canonical global key: the two
   endpoint coordinates on the decision grid, lexicographically ordered) is mapped to a list of
   *portions* `[s0, s1) -> feature id` whose lengths sum to the segment length, or to exactly
   one *exclusion record* `{class, reason}`. A segment is never both assigned and excluded, and
   no portion of an assigned segment is claimed by two features (A1: partition by length and by
   count, multiplicity 1).
2. **Vertex table.** Every perimeter vertex whose turning angle exceeds the recorded corner
   threshold (`CornerTurnToleranceDegrees` = 30, `metaledge.cpp`: `corner_angle_tolerance_degrees`,
   compared on the 1e-12 direction grid) and every endpoint / junction vertex that is not a
   simulation cut is mapped to exactly one vertex feature (corner / endpoint / junction) or to
   exactly one cluster (membership). Vertices on the truncation boundary are cuts, not
   features, and are listed with `Class: "TruncationCut"`.
3. **Feature list** with canonical **signatures**: `Type` plus dimensionless parameters only —
   separations / R, offsets / R, angles in degrees, corner radii / R, gap sides, interface
   types, conductor topology (same / different conductor by connectivity, relabelled by first
   appearance in the canonical order), metal boundary law. The signature is translation and
   rotation invariant; mirror images produce the same signature and differ only in the
   `Chirality` flag (+1 / -1), which is *not* part of the signature hash. Every feature carries
   `Hash` = SHA-256 of its canonical signature string, its claimed length, and its members
   (segment portions, vertices) in mesh coordinates for the audit.
4. **Exclusions** with class, reason, count and length, so that
   `sum(assigned) + sum(excluded) = total perimeter length` to roundoff.
5. **GeometryDigest**: SHA-256 over the sorted feature signatures (with multiplicity) and the
   exclusion classes. Identical for any rank count (A4), for any library (A3), and for
   refinements of the same layout up to the documented segment split (A5). Feature lengths are
   continuous quantities whose sums differ at roundoff between meshes of the same layout; they
   are not hashed but compared per feature with a tolerance by the audit.

**Library matching is a separate pass.** It looks features up by signature (the library models'
signatures are computed from their stored geometry — separation, angle, corner radius, cluster
`Edges` / sites — by the same canonicalisation) and can only mark a feature `Matched` /
`Missing`. An unmatched feature disables the correction for **that feature only**; no interface
group, chain, or neighbour is affected. Existing library models were discovered under the old
classification and are expected to be `Missing` (decisions 70 / 71: regeneration deferred).

## (b) Feature definitions from the whole geometry

Notation: R = matching radius (mesh units); all thresholds are multiples of R; all decisions
are quantized as in (c). A *chain* is a maximal perimeter path through regular vertices
(`physical_chain`, `metaledge.cpp`): a polyline of straight runs joined at sub-corner vertices
(item 7, the curved-edge chain rule). Chain tangents are canonical: from the lexicographically
smaller endpoint.

1. **Through-vertex pairs.** Two points p, q on two chains that meet at a vertex v are
   *through-vertex* when both are within 2R of v. Such pairs are described by the vertex
   feature (or by the cluster the vertex belongs to), never by an interaction event. This is
   the existing convention (`ConnectedNearVertex`) restated on points instead of mesh
   segments so that it is mesh independent; it makes corners with interior angle >= 60 deg
   plain corners and lets the arms of sharper corners interact beyond 2R from the vertex
   (2a sin(theta/2) < 2R for a >= 2R has solutions iff theta < 60 deg).
2. **Parallel pairs.** Two single-run chains whose tangents are parallel (|t1 . t2| >= 1 - 1e-8
   on the direction grid; chains with joints pair through item 7) with in-plane separation d < 2R and overlapping longitudinal extent form a
   *parallel interaction* over the overlap interval. Gap directions facing each other -> `Gap`
   (same or different conductor from connectivity); pointing away from each other -> `Strip`
   (the interval between them is one metal strip); the same direction -> impossible in one
   metal plane, recorded as exclusion `IncompatibleParallelPair`. Parallel interactions do not
   create events. Along a common longitudinal axis the chains are split at every projected
   chain endpoint; between consecutive split points the set of mutually interacting parallel
   chains (connected through d < 2R links) is constant and forms one *translational feature*
   over that interval: 2 chains -> `SameConductorGap` / `DifferentConductorGap` /
   `SameConductorStrip` with `Separation / R`; >= 3 chains -> `ParallelEdgeCluster` with sorted
   offsets / R, gap sides and conductor labels. Adjacent intervals with the same chain set and
   the same signature merge. The feature claims the interval on **every** participating chain.
3. **Interaction events and clusters.** An *event* is a pair of points (p, q) on distinct,
   non-parallel chains, not through-vertex, with |p - q| < 2R. The event set of a chain pair
   is the sub-interval of each chain whose distance to the other chain's sub-interval is
   < 2R (computed on the chains' straight geometry, not per mesh segment). The *cluster
   region* is the connected union of radius-R balls around all event points: two event cores
   belong to the same cluster when their distance is < 2R (union-find in canonical order;
   the result does not depend on the order). Every chain portion inside the region (distance
   to an event core < R: an interval on each chain, solved analytically) and every vertex
   inside the region belongs to that cluster. A vertex feature whose through-vertex zone (2R)
   reaches into a cluster region (distance from the vertex to an event core < 3R) joins the
   cluster together with its own window, so an acute corner and the cluster its arms form are
   **one description** (no corner record next to the cluster). Portions are split at the
   region boundary canonically (the analytic interval endpoints on the chain, on the decision
   grid). Two vertex features closer than 2R have overlapping radius-R windows (invariant A2)
   and are an event of their own (both vertices are degenerate event cores), so the corners of
   a strip end or of an aperture narrower than 2R form one cluster instead of overlapping
   corner descriptions. Cluster portions have priority over vertex windows and parallel
   features.
4. **Vertex features.** A corner (2 chains, turn > 30 deg), endpoint (1 chain) or junction
   (>= 3 chains) that is not inside a cluster claims a *window* of length R along each of its
   chains, shortened to half the chain length when the chain ends at another vertex feature
   within 2R (two corners R apart on a strip end each claim half of the end edge). Corner
   signature: `ConvexCorner` / `ConcaveCorner` (from the gap side), interior angle on the gap
   side in degrees, `CornerRadius / R` (0 for a sharp corner; the rounded-run rule of the
   legacy classifier is retained for fillets on runs instead of mesh vertices: a maximal
   sequence of straight runs shorter than R with sub-threshold turns at both ends, bounded by
   two longer arm runs, whose tangent distances from the virtual corner are < R and equal
   within 5 % and whose fillet radius is in (0, R); refinement cannot break the sequence
   because collinear mesh vertices merge into one run), interface types and boundary law. Endpoint: interface types, law. Junction: sorted arm angles.
5. **Isolated edges.** Every chain portion not claimed by a cluster, a vertex window or a
   translational feature is an `IsolatedEdge` portion; one feature per straight run
   (chain), signature = interface types + boundary law.
7. **Curved-edge chain rule** (decision 73(1); phase 2). A chain is defined by its
   significant vertices only: collinear splits (refinement midpoints, second-order mid-edge
   nodes) merge into one run, and the sub-corner joints between runs (turn <= 30 deg) are the
   polyline's bends. The turn at every joint (except the joints of a detected fillet, which
   the rounded corner accounts for) is spread over the two adjacent half-chords — the
   discrete curvature density, exact for a polygon inscribed in a circle at any chord length —
   and the *windowed curvature* at a chain point is the mean density over a window of length
   `CurvatureWindowOverR` = 1 x R centred on it (the response at a point integrates the geometry
   within ~R; the window is clipped at the ends of an open chain and periodic on a closed one).
   The windowed bend radius is its inverse. A chain point is *curved* when the windowed bend
   radius is below `StraightBendRadiusOverR` = 10 x R (quantized strict less), otherwise
   *straight-like*; the crossings are solved on the piecewise-linear windowed curvature so
   that they do not depend on the mesh. Rationale for 10 (decision 75, 2026-09-24): the
   first-order curvature correction to an edge response scales as R / radius (the response
   integrates the field over distances <= R from the edge and an in-plane bend perturbs that
   geometry at relative order R / radius), so at 10 R it is <= 10 % of the local edge
   correction — ~1e-3 of the corrected edge participation for edge corrections of a few per
   cent; the transmon's 38.9 um = 19.4 R CPW bends are straight-like, and curved coupons for
   1 < radius / R < 10 are deferred to the library regeneration (phase 2 used 20, i.e. 5 %; no
   synthetic class changes between the two). Consequences:
   * a straight-like chain portion is described by the straight features (isolated edge,
     pair) with a `BendRadiusOverR` annotation on every feature (the tightest windowed radius
     over its portions; not hashed; null on straight chains);
   * an unpaired curved portion is a `CurvedEdge` feature (one per curved chain section)
     with `RadiusOverR` = the section's tightest windowed radius in the signature;
   * a fillet at a corner (radius < R; the run-based rounded-corner rule of item 4, computed
     from the arms' accumulated turn and the tangent distances, hence refinement-invariant)
     is a rounded corner and takes no part in the curvature;
   * **pairs along bends**: two chains that are not both single straight runs are a
     constant-separation pair when their closest-point separation over the mutually paired
     intervals (points within the candidate reach `PairCandidateReachOverR` = 2 (1 + 0.05) R
     of the other chain, not beyond either end of it — the half-plane past a chain end along
     its outward tangent — and outside the 2R zones of shared vertices) is constant within
     `PairSeparationToleranceRelative` = 0.05 of the minimum (a polyline of sub-corner turns
     at constant width varies by at most 1 / cos(15 deg) - 1 = 3.5 %; the pair response
     sensitivity d dR/dd is O(1)). **Whether the pair interacts is decided at the chain
     level on the separation of the underlying curves** (`PairSeparationEstimate`), so that
     the discretisation never changes the classification: the chords of a polyline inscribed
     in a curve lie inside it (sagitta c^2 / (8 rho); two concentric inscribed polylines are
     w cos(turn / 2) apart mid-chord and exactly w apart at their vertices), while the exact
     offset polyline of a bent path keeps corresponding chords at the design separation and
     the outer side's samples near the joints project onto the inner vertices at up to
     w / cos(turn / 2). In both constructions the sampled closest-point distance from one
     chain to the other reaches the curve separation w as its maximum on the side whose
     maximum is smaller: separation = min over the two directions of the maximum sampled
     distance (the samples include the run ends, i.e. the vertices; 16 per piece). The pair
     interacts iff separation < 2R on the quantized grid — the same strict-less decision as
     a straight parallel pair at that separation (item (c)): a CPW gap of exactly 2R along a
     bend is isolated edges like a straight one (DS-SCT-001's 4 um gaps at R = 2 um dipped to
     3.9999 mid-chord along its 250 um bends and had formed three clusters of 874-1,184 edges
     claiming every strip). An interacting pair claims both chains and is split by curvature
     class: the straight-like pieces form a `SameConductorStrip` / `SameConductorGap` /
     `DifferentConductorGap` (separation = the pair separation above on the signature grid,
     one value per chain pair), the curved pieces a `CurvedSameConductorStrip` /
     `CurvedSameConductorGap` / `CurvedDifferentConductorGap` with `RadiusOverR` = the
     tightest windowed radius of the two sides (the inner side of concentric arcs). The
     cross-chord interactions of a constant-separation pair are never events, whether or not
     it interacts, so the concentric chords of a bend never form a spatial cluster and
     nothing is omitted as "nonparallel". Chains whose separation is not constant (acute
     corner arms, tapers) keep the event rule of item 3. Two single straight runs keep the
     translational rule of item 2 (exactly parallel or events).
   Limitations recorded: a slowly tapering pair (separation drift > 5 % over the paired
   interval) is a cluster, as before, and so is a pair whose separation is constant over most
   of its length but not all (the constancy test is per chain pair; along its bends such a
   pair at exactly 2R would again dip into events); the joint between a straight lead and a coarse polyline
   arc is intrinsically ambiguous (its turn is spread over the adjacent half-chords, so up to
   half a lead may join the curved class at coarse discretisations — classes are
   discretisation-independent, lengths are not); a chain pair with two bends of different
   radii yields one curved feature with the tighter radius.
8. **Exclusions** (decision 73(3), recorded with length): `TruncationCut` (segments on the
   simulation boundary), `Untargeted` (no target interface on the segment),
   `NonPlanar` (segment process normal not parallel to the reference process normal: walls,
   staples), `CrossLayer` (planar segment whose offset along the process normal differs from
   the primary metal plane — the plane carrying the largest perimeter length), and
   `IncompatibleParallelPair`. Perimeter that never reaches the classifier (metal sheets with
   the same material on both sides cancel in the odd-incidence perimeter, see the baseline
   report) is a phase-4 item and is *not* covered by this table.

**Signature grids.** Every continuous signature quantity (portion endpoints / R, offsets and
separations / R, corner radii / R, gap-direction cosines) is rounded to the recorded grid
`SignatureLengthQuantumOverR` = 1e-6 and every angle to `SignatureAngleQuantumDegrees` = 1e-6
deg. Rationale: the coordinates entering a signature are chip-scale mesh coordinates and
bisection results whose roundoff is ~ulp(|p|) (1e-12 um at 10 mm), so the grid must be far
above it for translated / rotated copies of one feature to hash identically (flip probability
per coordinate ~ roundoff / grid ~ 1e-6 at R = 2 um), and it must be far below any resolution the
response can depend on (the response varies on the scale R): 1e-6 R = 2 pm at R = 2 um. The
decision grid 1e-8 R of (c) is for geometric decisions (interaction / membership), not for
hashed values. Version-1 records derived from a signature (Separation, CornerRadius) inherit
the signature grid.

**Cluster frame and signature.** Origin = length-weighted centroid of the cluster's claimed
portions; z = process normal. The in-plane x axis is chosen among the finite candidate set
{chain tangents of the cluster, both signs, and their in-plane perpendiculars}; for each
candidate and each handedness (y = z x x or y = -(z x x)) the cluster is serialised (portions
as `[x0, y0, x1, y1] / R` on the signature grid, gap side, interface types, conductor label by
first appearance, vertices as `[x, y] / R` + type + turn) and sorted; the lexicographically
smallest serialisation is the signature, and the handedness that produced it is the
`Chirality` (+1 / -1; 0 when both handedness values reach the minimal serialisation, i.e. the
cluster is its own mirror image). A mirror image yields the same signature with the opposite
chirality; a rotated or translated copy yields the same signature and chirality. Symmetric
clusters tie between equivalent frames and produce the same string. Cost recorded (phase 3,
PENDING): the frame search serialises the cluster once per distinct portion direction — O(directions x
portions log portions) JSON serialisations — so a cluster of ~1,000 edges with ~1,000 distinct
directions (DS-SCT-001 before the pair-along-bend fix) took ~290 s; chip-scale clusters need a bounded
frame rule (e.g. candidates restricted to the frames minimising the first serialised entry, an exact
pruning of the same lexicographic rule). This replaces the representative-event site
(first closest pair in candidate order), the exhaustive spatial closure (4 x 2R merging) and the
"nonparallel" omission of the legacy classifier.

## (c) Behaviour at exactly R and 2R

All length decisions use the decision-69 quantizer: lengths are rounded to the grid
`1e-8 R` (`kDecisionLengthQuantumRelativeToMatchingRadius`), direction cosines to `1e-12`
(`kDecisionDirectionQuantum`), and every "within" test is a strict `<` on the quantized values.
Consequences (recorded in the manifest under `Library.DecisionQuantization` and
`Identification.Conventions`):

* two parallel edges at exactly 2R do **not** interact (both isolated edges); at 2R - quantum
  they form a pair;
* an event needs |p - q| < 2R; through-vertex exclusion applies when either point is < 2R from
  the shared vertex; a vertex joins a cluster when its distance to an event core is < 3R;
* a chain portion belongs to a cluster when its distance to an event core is < R;
* the strip at exactly R (`strip-2`, the transmon's 64 segments) is a `SameConductorStrip` at
  `Separation / R = 1` on both sides and its corners are plain corners: no knife edge between
  1.95 / 2 / 2.05 um beyond the separation value itself.

## (d) Manifest version 2

`surface-response-requirements.json` keeps every version-1 field (`Requirements` with
`Topology`, `Geometry`, `Interfaces`, `BoundaryCondition`, `Status`, `Count`,
`TotalEdgeLength`, `Summary`, `Library`, `Statistics`) so that `manifest.py`,
`preflight_matrix.py` and the qualify tooling keep working; `Version` becomes 2 and the
version-1 `Requirements` are **derived from the version-2 features** (one record per distinct
signature with `Count` = features or mesh segments as before and `Status` from the key-based
matching pass). The new top-level `Identification` object carries the contract:

```
"Identification": {
  "Version": 2,
  "MatchingRadius": R,
  "Conventions": {"CornerTurnToleranceDegrees": 30, "InteractionDistanceOverR": 2,
                  "ThroughVertexZoneOverR": 2, "ClusterBallOverR": 1,
                  "VertexJoinsClusterOverR": 3, "VertexWindowOverR": 1,
                  "ParallelCosineTolerance": 1e-8, "RoundedCornerTangentTolerance": 0.05,
                  "SignatureLengthQuantumOverR": 1e-6, "SignatureAngleQuantumDegrees": 1e-6,
                  "StraightBendRadiusOverR": 10, "CurvatureWindowOverR": 1,
                  "PairSeparationToleranceRelative": 0.05, "PairSeparationSamplesPerInterval": 16,
                  "PairSeparationEstimate": "min over the two chains of the maximum sampled closest-point distance to the other chain (the curve separation at the vertices); interacting iff < 2R",
                  "PairCandidateReachOverR": 2.1,
                  "Comparison": "strict less on the quantized grid"},
  "ReferenceProcessNormal": [nx, ny, nz],
  "Features": [ {"Id": k, "Type": "...", "Signature": {...}, "Hash": "sha256", "Chirality": +-1,
                 "BendRadiusOverR": r | null,
                 "Length": L, "Portions": [[segment, s0, s1], ...], "Vertices": [v, ...],
                 "Frame": {"Origin": [...], "Axes": [[...],[...],[...]]},
                 "Match": {"Status": "Matched" | "Missing", "Model": "name"} } ],
  "Segments":  [ {"Key": [[x0,y0,z0],[x1,y1,z1]], "Length": L, "Chain": c,
                  "Portions": [[s0, s1, feature], ...] } | {"Key": ..., "Length": L,
                  "Exclusion": {"Class": "...", "Reason": "..."}} ],
  "Vertices":  [ {"Point": [...], "Type": "Corner|Endpoint|Junction", "TurnDegrees": t,
                  "Feature": k} | {"Point": ..., "Class": "TruncationCut"} ],
  "Exclusions": [ {"Class": "...", "Reason": "...", "Count": n, "Length": L} ],
  "Totals": {"PerimeterLength": L, "AssignedLength": La, "ExcludedLength": Le},
  "GeometryDigest": "sha256"
}
```

Lengths and coordinates are in mesh units on the `ScaleLength` grid of the version-1 manifest.
The audit reads `Identification` when present (A1 partition from `Segments`, vertex census from
`Vertices`, exclusions from `Exclusions`, A2 from the cluster frames, A3 / A5 from
`GeometryDigest` and the feature signatures) and falls back to the version-1 aggregate reading
otherwise.
