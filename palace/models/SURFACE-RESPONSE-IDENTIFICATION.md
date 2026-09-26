# Surface-response geometry identification: the canonical contract (design, 2026-09-24)

Scope: the geometry-identification step of the fabrication-process surface-response
correction (decisions 69-74 of `coupon-accuracy-assessment-20260913/SUPERVISOR-DECISIONS.md`,
invariants A1-A6 of `GEOMETRY-IDENTIFICATION-PLAN.md`). This document is the contract that the
implementation in `surfaceresponseidentification.{hpp,cpp}` follows, that the preflight
manifest (`surface-response-requirements.json`, version 2) serialises, that the audit tool
(`examples/surface_response_identification/audit.py`) gates, and that the solve-path patch
construction (`BuildFeaturePatches`, `surfaceresponseoperator.cpp`) consumes verbatim: (a)-(d)
the identification and the manifest (phases 1-3), (e) the patches built from the features
(phase 4; the legacy classification stays behind `PatchConstruction = "Legacy"`).

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

**Numbering (decision 82 infrastructure, 2026-09-25).** The point, face, segment and vertex
numbering of the perimeter is a function of the set of distinct points only: the canonical
points are renumbered by their sorted quantized coordinates (the 1e-10 x extent grid, ties by
the raw coordinates; the representative of a merged point is its lexicographically smallest
copy), the distinct faces by their sorted canonical vertex sets, the segments by their point
pairs and the vertices by first appearance in segment order. The whole manifest apart from
`Statistics` is therefore byte-identical at any rank count and for any gathered order of the
crack copies (the first-encounter numbering had permuted 2,647 segment-table entries between
np8 and np192 on DS-SCT-002). The metal faces are gathered on the root only, the perimeter and
the identification are computed on the root and the compact results (segments, vertices; the
feature list, segment and vertex tables, exclusions) are broadcast; every rank receives only
its own retained facets. No rank other than the root holds the gathered faces, the global
faces or the identification state.

**Vertex census rule.** The vertex table (item 2) is the census of the vertices of the
classifier's one-sided (PHYSICAL-type) segments — non-planar faces and one-sided box edges
included, folds and non-manifold edges not. A vertex all of whose incident edges are folds /
non-manifold edges (the corner of a PEC box where three box faces meet, the base corners of a
bump) is a vertex of no one-sided perimeter and has no record on either side; the audit's
`physical_kind` counts the same set (a `BOX` edge shared by two box faces is a fold).

**Library matching is a separate pass.** It looks features up by signature (the library models'
signatures are computed from their stored geometry — separation, angle, corner radius, cluster
`Edges` / sites — by the same canonicalisation): a model of the feature's topology whose
parameters lie within the signature tolerance (decision 85(1), item 7 below; the nearest one),
and can only mark a feature `Matched` / `Missing`. An unmatched feature disables the correction for **that feature only**; no interface
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
   inside the region belongs to that cluster. A vertex feature within the interaction
   distance of an event core (distance from the vertex to a core < 2R: its radius-R window
   overlaps the core's radius-R ball) joins the cluster together with its own window
   (`VertexJoinsClusterOverR` = 2; decision 82(1), 2026-09-25: the former 3R join was the one
   distance differing from 2R and joined the L2 junction regions of DS-SCT-002 to the L1
   flux-loop ends 2.4R below; an acute corner whose arms interact from 2R along the arms is
   now a corner feature next to a separate arm cluster, abutting at R along each arm).
   **Arc rule (decision 82(3), 2026-09-25; replaces the run-based rounded-corner rule of
   item 4 and the half-chord curvature of fitted arcs in item 7).** Along every perimeter
   path through vertices with exactly two path segments (regular or corner; a path stops at
   endpoints, junctions and cuts) the joints are the non-collinear vertices. From every
   unconsumed joint, the largest range of at least three following joints that are
   connected by pieces shorter than the interaction distance 2R, turn the same way, total
   at most 180 deg and fit ONE circle with the two arm tangents (`ArcFitToleranceRelative`
   = 0.05: tangent lengths from the virtual corner equal within 5 %, every joint within 5 %
   of the radius from the centre; antiparallel arms: the radius is half their separation
   and the tangent points face each other) is an arc; a closed loop starts its scan after
   its longest piece. **Amendments of decision 85(1) (2026-09-26):** (i) a piece of 2R or
   longer stops the arc only when either of its joints is a corner-class turn (> 30 deg: a
   polygonal design stays corners); between two sub-corner joints a long piece is a chord
   of a smooth polyline bend (DS-SCT-001's 367 / 371 / 373 um route bends are meshed with
   4 um chords = exactly 2R at R = 2 um and the exact stack separations are their radius
   differences), and in that long-chord regime the scan stops at the first failed fit (the
   short-chord scan keeps its full range so that a perturbed joint inside a fillet is still
   absorbed by a larger range); (ii) a BEND's circle (radius >= R) is the algebraic
   least-squares fit of its joint vertices — exact for an inscribed polyline and a set
   function of the joints; the tangent-length construction assumes arms tangent to the
   circle, which a route bend preceded by a spline piece is not (0.5 % radius error and a
   0.8 um centre offset on DS-SCT-001's 370 um bends) — while a rounded corner keeps the
   tangent-length radius (its arms are tangent by construction; the fillet gate is
   unchanged); (iii) a closed path whose joints all turn one way through 360 deg on one
   circle (a round pad, hole or via) is ONE arc of total turn 2 pi and a bend of exact
   radius whatever its radius: a circle of radius < R is one `CurvedEdge` with
   `RadiusOverR` < 1, never two 180 deg "rounded corners" split at a numbering-dependent
   joint (review m6; every joint must be sub-corner: a square hole is a 4-corner polygon);
   (iv) both traversal directions of every path are scanned and the arc set absorbing more
   joints is applied (ties: more arcs, then the smaller serialisation) — the greedy scan
   from the first unconsumed joint depended on the input orientation, and a mirrored mesh
   now gives the mirrored arcs. Gate (`permute_msh2.py`): DS-SCT-001 and a synthetic stack
   renumbered (seeded node / element permutation) give identical manifest content, mirrored
   in x the same signature multiset, lengths and digest with the cluster chirality negated.
   **Self-pairing point-wise (2026-09-26):** the partner search of a chain facing itself and
   the self events exclude the part of a run within pi R of arc length of the POINT (the
   former run-level exclusion dropped a whole 28 um leg of a hairpin for every point of its
   fold, so that the fold end facing the far leg at 1.6R had no partner and no event; it is
   a cluster now). Two-joint polylines are never arcs (a chamfer, and a square strip end
   — the diameter chord of a semicircle — cannot be told from a one-chord arc: they stay
   corners). An arc of radius < R whose total turn exceeds the corner threshold is ONE
   rounded corner: `ConvexCorner` / `ConcaveCorner` by the side of the centre (convex when
   the centre lies on the metal side of the first arm; well defined for a U-turn),
   `AngleDegrees` = 180 - total turn, `CornerRadiusOverR` = radius / R from the tangent
   lengths (exact for an inscribed polygon at any chord count), claiming the arc runs and
   R along each arm; it is a chain of its own between its tangent points and its two arms
   meet THROUGH it (the through-vertex zone of item 1 is taken within 2R of the arc's runs),
   so the event and pair rules read a filleted corner like a sharp one and a U-turned
   narrow strip keeps its strip pair; like a sharp corner it separates its arms into
   distinct chains (one isolated edge per arm). An arc of radius >= R is a bend of exactly
   that radius: it stays inside its chain, its joints (corner vertices included) contribute
   no half-chord density, the density over the arc is 1 / radius, a curved section on it
   reads `RadiusOverR` = radius / R exactly, and the chains a corner vertex inside it
   separated are merged (a coarse bend with a super-threshold joint is the same chain as a
   fine one). Corner vertices absorbed by an arc are `RoundedCornerVertex` (feature = the
   rounded corner) / `BendVertex` records of the vertex table, never corner features. The
   description therefore does not depend on the number of chords as long as every chord is
   shorter than 2R (a coarser polyline is a different geometry at the scale of R: its kinks
   are real corners); gate: synthetic fillets at rho / R in {0.1 ... 20} x turns {45, 90,
   135, 180} x {2, 4, 8, 16} chords + a seeded chord perturbation. Recorded limitations: a
   chain never pairs with itself (a bend of radius >= R folding an edge back onto itself
   within 2R is not paired — a wide U-turn is beyond 2R anyway); the polyline's inscribed
   circle differs from a design radius by the discretisation (an offset polyline at 20 deg
   per vertex reads 3 % off; within the fit tolerance). Portions are split at the
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
   (chain), signature = interface types + boundary law. Claims are resolved per run in the
   order cluster > vertex window > translational; a claimed or unclaimed piece shorter than
   `SignatureLengthQuantumOverR` x R is roundoff between two claim boundaries (a cluster ball
   cutting a pair piece next to a run end) and joins the adjacent portion on its run — a
   feature whose whole length is below the signature grid is not a feature (DS-SCT-001 had
   three `CurvedSameConductorStrip` records of 1.6e-7 um in total).
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
   * **pairs along bends** (phase 3: local constancy): two chains that are not two exactly
     parallel straight runs (those keep the translational rule of item 2) are sampled over
     their candidate facing regions — the points within `PairCandidateReachOverR` = 2 (1 + 0.05)
     R of the other chain, not beyond either end of it (the half-plane past a chain end along
     its outward tangent), outside the 2R zones of shared vertices, cut at the curved
     boundaries, at most `PairSampleSpacingOverR` = 0.5 R apart. A sample is **locally
     constant** when the sampled distances within `PairConstancyWindowOverR` = 1 R of it along
     its own chain vary by at most `PairSeparationToleranceRelative` = 0.05 of their minimum (a
     polyline of sub-corner turns at constant width varies by at most 1 / cos(15 deg) - 1 =
     3.5 %; the pair response sensitivity d dR/dd is O(1)). The constant portions are the pair;
     the portions that are not (divergence at tees and port ends, fast tapers, acute corner
     arms) keep the event rule of item 3 — a slow taper is a pair, a tee is a cluster. **Whether
     a constant portion interacts is decided on the separation of the underlying curves**
     (`PairSeparationEstimate`), so that the discretisation never changes the classification:
     the chords of a polyline inscribed in a curve lie inside it (two concentric inscribed
     polylines are w cos(turn / 2) apart mid-chord and exactly w apart at their vertices),
     while the exact offset polyline of a bent path keeps corresponding chords at the design
     separation and the outer side's samples near the joints project onto the inner vertices
     at up to w / cos(turn / 2). In both constructions the sampled closest-point distance from
     one chain to the other reaches the curve separation w as its maximum on the side whose
     maximum is smaller: a sample's separation = min over the two chains of the maximum
     sampled distance within a window of half-width max(R, local chord) where the chain bends
     (windowed curvature > 0: the window holds a vertex of an inscribed polyline and a full
     chord of an offset polyline) and R on straight runs (no dip; a taper is read locally),
     about the sample on its own chain and about its foot on the other chain. That chord
     reading C is exact for an offset polyline; two polylines inscribed in the curves at
     aligned angles are C = w cos(turn / 2) apart everywhere (chords and vertex-to-polyline
     alike) for a curve separation w, and the polyline pair alone cannot tell the two
     constructions apart (they differ at order turn^2): the inscribed reading is
     C / cos(turn / 2) with turn = the larger local joint turn of the two chains. **The portion
     interacts iff both readings are below 2R** on the quantized grid — the same strict-less
     decision as a straight parallel pair (item (c)), taken on the non-interacting side of the
     recorded discretisation ambiguity w (1 / cos(turn / 2) - 1) (4e-5 w at 1 deg per vertex,
     1e-3 w at 5 deg, 9e-3 w at 15 deg): a CPW gap of exactly 2R along a bend is isolated
     edges like a straight one, and a gap of 2R - 1e-3 R interacts only where the joint turn
     is below 2 acos(1 - 5e-4) = 3.6 deg (DS-SCT-001's 4 um gaps at R = 2 um read 3.9998
     mid-chord along its 250 um bends (1.1 deg joints) and had formed three clusters of
     874-1,184 edges claiming every strip; its pairs read 3.999-4.200 over the whole facing
     region because of the tee / port ends, which a global constancy test rejected). An interacting portion set claims both
     chains and is split by curvature class: the straight-like pieces form a
     `SameConductorStrip` / `SameConductorGap` / `DifferentConductorGap` (separation: one value
     per chain pair and class by the EXACT-parameter rule of decision 85(1) below — formerly
     the sample-count-weighted mean chord reading, which carried the sagitta bias of the
     coarser chords of a polyline bend, c^2 / (8 rho) = 31 nm = 1.5 % of a 2 um gap for
     7.9 um chords on a 370 um bend, so that one design cross-section hashed to several
     keys: review B1), the curved pieces a `CurvedSameConductorStrip` /
     `CurvedSameConductorGap` / `CurvedDifferentConductorGap` with `RadiusOverR` = the tightest
     windowed radius of the two sides (the inner side of concentric arcs). The cross-chord
     interactions of a locally constant portion are never events, whether or not it
     interacts, so the concentric chords of a bend never form a spatial cluster and nothing
     is omitted as "nonparallel".
   Phase 4 refinements found by the patch construction on DS-SCT-001 (2 um trace between
   2 um gaps: trace + gap = 2R): (i) "not beyond either end of the other chain" is local —
   past the end plane along the outward tangent AND within the reach of that end (the
   half-plane alone is right for a straight partner but past the end plane of a bent route
   it covered points facing the partner's interior, so the bent ground edges lost their
   whole facing region and the pairs were described on one side only); (ii) the constant,
   interacting pieces of a chain pair are grouped by separation (split where consecutive
   mean separations differ by more than the 5 % pair tolerance) and each group is its own
   pair with its own frame and class — a closed trace loop faces the same ground chain at
   the gap (near side) and across the strip (far side, 2R); (iii) a pair or parallel cluster
   whose claim on one of its sides is lost entirely to an earlier claim of the same priority
   is not a pair: its surviving pieces return to the chain's isolated / curved edge (claim
   resolution, `Sides` in the manifest = the signature's edge order for chirality +1). Still
   open (PENDING, supervisor): a bent CPW narrower than 2R is a multi-chain neighbourhood
   (ground - gap - trace - gap - ground within 2R) that the pairwise bent rule cannot
   express and the translational component rule only covers for straight runs; on
   DS-SCT-001 this leaves knife-edge pieces at exactly 2R across the trace whose two sides
   read differently (unequal chords: 7.5 um ground chords against 4 um trace chords) — the
   patch construction refuses such a pair (its sides do not face each other at its
   separation within 10 %, twice the pair tolerance) and reports it — RESOLVED by the
   stack rule below (2026-09-25).
   **Stack rule (decision 82(2), 2026-09-25; replaces the PENDING item above and the
   decision-78 tie-break).** The locally constant, interacting facing relations of the
   bent-pair rule (one *link* per chain pair and separation group) and the translational
   spans of the rigid parallel runs are the pairwise facing relations of the perimeter.
   Links and spans sharing a run over a common interval (union-find in canonical order) are
   ONE cross-section component; a lone two-edge link or span is the pair feature above; every
   other component is assembled per cross-section: along every member chain the pieces are
   cut at the piece ends of every link on it (per curvature class), at the ends of the
   higher-priority claims on it (cluster portions, vertex windows) and at the images of the
   other members' cuts through the links (the foot on the partner chain), and the composition
   at the middle of every elementary interval is read off the active links: from the chain,
   the partner chains reached through the links (breadth first, each at the foot of the
   previous point; a chain facing itself takes its foot outside the pi R neighbourhood),
   ordered along the in-plane normal. k >= 3 edges = a `ParallelEdgeCluster`
   (`CurvedParallelEdgeCluster` where a member is curved: `RadiusOverR` = the tightest
   windowed radius over the members), k = 2 the pair classes; the signature offsets are the
   sums of the CONSECUTIVE links' separations (a non-consecutive link within 2R, e.g. the
   outer edges of a 4-edge stack at 4 um and R = 2.1 um, does not enter), the gap sides and
   conductors are read at the members; sides in the canonical signature order (chirality +1
   for every asymmetric cross-section; a symmetric one, chirality 0, puts the member with the
   smallest (chain, position) on side 0); one feature per (type, signature, curvature class,
   member chains) — a straight stack on the two leads of a bend is ONE feature, the bend a
   second. **Stack-end rule:** a member taken by a cluster portion or a vertex window is not
   part of the cross-section and is not traversed: the stack ends at the claim boundary and
   the remaining members are recomposed there (a smaller stack, a pair, or nothing — the
   member alone returns to the chain's isolated / curved edge as a recorded
   `ClusterNeighbour` / `VertexNeighbour` — superseded by decision 85(2) below: the member
   returns to the single-edge remainder, where the cluster extension tests it, and every
   component, a lone two-edge link or span included, is assembled per cross-section with
   this taken rule since 2026-09-26). The pairwise candidates inside a stack are
   superseded by it: no two claims of the pair priority ever overlap and claim resolution never
   decides by feature id (`Diagnostics.SamePriorityClaimOverlaps` = 0 is gated). The chord
   reading of a locally constant sample takes the window maximum over the locally CONSTANT
   samples only (the window of a full chord crossing a taper kink read 2.15 - 2.7 um for a
   2 um gap on DS-SCT-001 and formed a separation group of its own). **Self-pairing:** a chain
   folding back onto itself pairs with itself: every point of a run is a candidate and its
   partner is the closest point of the chain at least `SelfPairNeighbourhoodOverR` = pi R of
   arc length away (on the tightest bend, radius R, the chord reaches 2R after half a turn =
   pi R of arc; Schur's comparison: two points closer than pi R along a chain of curvature
   <= 1 / R are within 2R of each other along any bend of radius >= R, so only points farther
   apart along the chain can face each other across a fold). Recorded geometric fact: a smooth
   chain of curvature <= 1 / R has its legs >= 2R apart after a 180 deg turn (displacement =
   int sin(theta) / kappa dtheta >= 2R), so a semicircular hairpin never faces itself within
   2R (legs closer than 2R have an inner fold below R = a rounded corner whose arms are
   distinct chains); the rule applies to non-circular folds of sub-corner joints (no arc
   fits) and to convergent legs. Tried and rejected (2026-09-25): making a rounded corner's
   foreign concentric neighbour event-eligible turned each DS-SCT-001 flux-loop end into one
   160 um cluster (the event set of a chain facing an arc reaches sqrt(3) R past the tangent
   points, the R balls another R, and the cores merged with the trace-junction clusters 1.5 R
   away); a rounded corner stays a vertex whose concentric neighbour is a constant non-event,
   and the pair / stack side facing its arc or window is a recorded `VertexNeighbour`.
   **Sampling margins (supervisor addition (b), verified):** the bent-pair candidate reach
   2R (1 + 0.05) enters only the candidate gathering (`RunIntervalWithin`,
   `SegmentSegmentDistance`); the interaction decision is `quantizer.Less(upper, 2R)` (strict,
   3D); the mutual-sides test at the pair separation x 1.05 trims a side whose partner piece
   is gone (a facing test); the assembly accepts a foot within the reach only along a link
   already decided interacting. **Facing gates (audit A8, `facing_check.py`):** the length of
   isolated / curved edge portions facing another edge of the plane within 2R and of pair /
   stack sides facing a third edge within 2R must be zero apart from the recorded exclusions
   `AtExactly2R` (strict rule), `ClusterNeighbour` (the facing portion is a cluster's: claim
   radius R), `VertexNeighbour` (corner window, rounded-corner arc, endpoint / junction
   window), `ThroughVertex` (sample and foot within 2R of one vertex feature),
   `SelfNeighbourhood` (own chain within pi R of arc); the own sides of a pair / stack are its
   member chains (superseded by decision 85(2) below: every facing segment within 2R is
   tested, the own sides are the portions plus the member chains within the feature's
   reach, and the `ClusterNeighbour` / `VertexNeighbour` exemptions are gone). Gate on synthetic stacks (straight k = 3..6 symmetric / asymmetric / around
   2R / wide ground, curved k = 3..6 at rho / R = 3, 8, 30, taper, U-ring, hairpins, kinked
   hairpin: 33 / 35 layouts, 2 meshes not loadable), DS-SCT-001 (three 4-edge stacks of 2.2 /
   1.2 / 3.2 mm, isolated facing 676 -> 0 um unexcluded), transmon, two-transmon chain.
   Limitations recorded: a taper faster than 5 % per 2R is events (a cluster), a slower one
   is a pair described by its mean separation (a slow taper crossing 2R is split at the
   crossing sample); the joint between a straight lead and a coarse polyline
   arc is intrinsically ambiguous (its turn is spread over the adjacent half-chords, so up to
   half a lead may join the curved class at coarse discretisations — classes are
   discretisation-independent, lengths are not); a chain pair with two bends of different
   radii yields one curved feature with the tighter radius.
   **Exact parameters and the signature tolerance (decision 85(1), 2026-09-26; review B1).**
   Every continuous parameter of a signature is an exact geometric reading where the
   geometry allows one: corner angles and junction arms from straight arm directions,
   rounded-corner radii from the tangent lengths, bend radii from the least-squares circle
   of the joints, and the pair / stack separations per link and curvature class as (a) the
   perpendicular distance of two exactly parallel straight runs neither of which bends
   within R of the sample or its foot (the leads of a route: constant and exact for the
   design lines; a chord of a coarse polyline bend has a joint within R and is not read
   this way — aligned inscribed chords of a 370 um bend are 31 nm inside the circle), or
   (b) the radius difference of two fitted bend arcs whose centres coincide within the
   parameter tolerance (the centre offset bounds the error of the difference; two local
   circles of a spline fitted piecewise are not concentric). A class with any exact piece
   takes the length-weighted mean of its exact pieces; otherwise it takes the
   length-weighted mean chord reading of the bent-pair rule (a slow taper, a spline route
   with chords beyond 2R) and the feature carries `ExactParameters` false (an annotation,
   not hashed). Sample-count weighting is gone. The stack offsets are the sums of the
   consecutive links' separations of the cross-section's class (a straight stack on the
   leads no longer carries the bend's readings). The chord reading keeps its recorded
   discretisation ambiguity (order c^2 / (8 rho) of the coarser chords), which the
   tolerance below absorbs where it is below 1e-3 R.
   **Signature tolerance and library grouping.** `SignatureParameterToleranceOverR` = 1e-3
   (offsets, separations, radii in R) and `SignatureAngleToleranceDegrees` = 1e-2 (corner
   angles, junction arms). Rationale: 1e-3 R is far above the chord ambiguity of sub-2-degree
   polyline joints (w (1 / cos(turn / 2) - 1) < 1e-4 R at w = 4R) and the signature grid
   (1e-6 R), and far below any separation difference the response resolves (pair sensitivity
   O(1) per R: 1e-3 of the pair correction); 1e-2 deg is a displacement of 1.7e-4 R at
   distance R, below the length tolerance, and the angles are exact anyway. The *topology
   key* of a signature is the signature with every continuous parameter (`OffsetOverR`,
   `SeparationOverR`, `RadiusOverR`, `CornerRadiusOverR`; `AngleDegrees`,
   `ArmAnglesDegrees`) removed. Two signatures of one type agree when their topology keys are
   equal — trying both orientations of a translational signature, since two instances of one
   asymmetric stack can canonicalise to mirror orientations when their offsets differ by less
   than the tolerance — and every parameter differs by at most its tolerance
   (`SignatureDeviation` = max |difference| / tolerance <= 1). The *library* groups feature
   instances of one type and topology whose parameters agree within the tolerance into ONE
   coupon: single linkage over the distinct signatures (union-find; the result does not
   depend on the instance order) at the representative signature `RepresentativeSignature`
   — every parameter the midpoint of its range over the group, each instance taken in the
   orientation nearest to the lexicographically smallest instance, substituted into that
   instance (a function of the set alone) — with the group's `Instances` and
   `ParameterSpread` (max normalised deviation of a member from the representative; a
   single-linkage chain wider than twice the tolerance would leave members unmatched and is
   visible here) recorded in the version-1 `Requirements` record and in the signature-only
   library (`signature_library.py`, the same rule in Python). The *matcher*
   (`LibrarySignatureIndex`, `surfaceresponseoperator.cpp`) indexes the library models by
   topology key (both orientations) and gives a feature the model of its topology within the
   tolerance with the smallest deviation, ties by model name — deterministic and independent
   of the feature / model order; `Match.Deviation` is recorded per feature. A
   `SpatialEdgeCluster` signature (a whole plan-view geometry in a canonical frame whose
   lexicographically minimal frame is discontinuous in the coordinates) has no tolerance: it
   matches exactly. Gate (A5 set identity, re-mesh): the stack suite, the bent-pair cases
   and the fillet suite meshed at >= 3 mesh sizes plus a seeded perturbation give the same
   topology keys and parameters within the tolerance; DS-SCT-001's three identical 2 / 2 / 2
   um routes give ONE key (0 / 1 / 2 / 3 R exactly, chirality 0; before: three keys 0.04-0.6 %
   apart).
   **Cluster extension (decision 85(2), 2026-09-26; review M4).** Every *single-edge portion*
   — the remainder of a run outside every cluster portion, vertex window, pair / stack claim
   and CrossLayer zone, i.e. what would become an `IsolatedEdge` / `CurvedEdge` — whose 3D
   distance to a cluster's claimed perimeter on another chain (or on its own chain beyond the
   pi R self-pair neighbourhood) is strictly below 2R, faced ACROSS (the perpendicular
   projection of the point onto the claimed piece falls inside the piece; at an interior
   joint of a claimed chain interval the piece's domain is extended by 2R tan(turn), the
   width of the wedge between consecutive perpendicular domains; never past the end of the
   claimed interval — a diagonal reach would creep along a sub-2R strip, whereas the
   partner beyond the claim end pairs with the free continuation), outside the through-vertex
   zones of the shared vertices that are NOT members of that cluster, joins the cluster over
   that sub-interval. The same for a vertex feature outside every cluster: a single-edge
   portion faced across within 2R of its window (outside its own through-vertex zone — the
   arms of a corner meet through it) makes that vertex a cluster together with the portion.
   A portion within 2R of several owners merges them (union-find over clusters and joined
   vertex features; order independent; the merged cluster is numbered by its smallest
   member). Pairs and stacks are joint descriptions and are never absorbed. The extension
   iterates to closure: an enlarged claim moves the stack ends (the stack-end rule
   recomposes the members on the new claims) and the recomposed stacks leave new single-edge
   portions to test; the loop stops when a pass absorbs less than the signature parameter
   tolerance 1e-3 R in total — below the resolution of every parameter and gate (DS-SCT-001:
   2 passes, 127 portions, 44.6 um; transmon 12.8 um; two-transmon chain 25.6 um — the
   former `ClusterNeighbour` lengths plus the recomposed pair-side pieces; DS-SCT-002 at
   R 2.1: pass 1 absorbed everything, passes 3-15 had re-cut 1 nm slivers for 8 s each). The pair / stack length within 2R of a cluster's claimed
   perimeter (across, outside the through-vertex zones) is the *stack-end third body*,
   reported as `Diagnostics.StackEndThirdBodyLength` (DS-SCT-001: 148 um) and by the facing
   gate as the recorded pair-side class `StackEndThirdBody`. Facing gates (audit A8,
   `facing_check.py`) after this decision: every facing segment within 2R is tested (the
   nearest alone masked a second edge behind an excluded one); the own sides of a pair /
   stack are the segments of its portions plus the member chains within the feature's reach
   (its lateral span x 1.05) of them, not the whole chains; the remaining exemptions are
   `AtExactly2R`, `ThroughVertex` (design item 1: the arms of one vertex feature meet through
   it — the identification's own through-vertex convention, kept), `SelfNeighbourhood`,
   and for pair / stack sides only `StackEndThirdBody` and `StackEndRecomposition` (facing
   another pair / stack that shares a member chain: the route's cross-section recomposed at a
   stack end or class boundary, the partner's foot across the cut; sub-chord lengths,
   DS-SCT-001 2.1 um), each with its length. Diagnostics also record the stack assembly's
   `StackCompositionCap` (64 members: far above any physical stack — the chip's widest has
   12 edges — bounding a runaway traversal through inconsistent links; hits counted,
   review m2) and `StackGeometricOffsetIntervals` (elementary intervals whose consecutive
   members had no active link, offsets from the geometric lateral distance). The
   mutual-sides facing test at separation x 1.05 lets a pair side overhang a truncated
   partner piece by at most sqrt(1.05^2 - 1) s = 0.32 s (`MutualSidesOverhangOverSeparation`,
   review m3).
8. **Plane rule (decision 82(1)).** Every run belongs to one metal plane (the distinct
   offsets of the run midpoints along the reference process normal on the decision grid).
   Events, core merging, site-site cores, the vertex-site join, bent pairs and the
   translational direction classes are taken between runs of ONE plane only: no pair,
   cluster or vertex join ever spans two planes. Metal of another plane within the
   interaction distance (2R, 3D, strict) is the `CrossLayer` exclusion of item 9 (the
   face-based zones cover every cross-plane edge pair within 2R, since the other plane's edge
   is the edge of a metal face). Distances recorded: interaction 2R (events, through-vertex
   zone, core merge, site-site, vertex join, CrossLayer reach); the cluster ball R is the
   claim radius of the region, not an interaction decision; the bent-pair candidate reach
   2R (1 + 0.05) is the SAMPLING margin of the constancy test (its interaction decision is
   the same strict < 2R on the curve separation) and the mutual-sides facing test uses the
   pair's own separation x 1.05 — both recorded here as the two non-2R constants that
   remain, neither decides an interaction.
   **Port rule (decision 82(5)).** Ports are not metal: a one-sided perimeter segment
   coincident with an edge of a LumpedPort / WavePort boundary face of the configuration
   (attributes from `Boundaries.LumpedPort[].Attributes` / `Elements[].Attributes` and
   `Boundaries.WavePort[].Attributes`, metal attributes excluded; never from names) is the
   `Port` exclusion (segment type PORT in `metaledge.cpp`, checked before the truncation
   test); the physical-edge graph stops there like at a truncation, so a vertex ending a
   port-bordering segment is a `PortCut` (never a corner / endpoint feature) and the metal
   edges running into the port are chains ending at the cut. A lumped port bridging the gap
   between two leads therefore leaves the lead ends excluded (2 x the lead width) with four
   PortCut vertices and the leads' long edges as isolated edges to the cut.
9. **Exclusions** (decision 73(3), recorded with length): `TruncationCut` (segments on the
   simulation boundary), `Port` (item 8), `Untargeted` (no target interface on the segment),
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
  the shared vertex; a vertex joins a cluster when its distance to an event core is < 2R;
* a chain portion belongs to a cluster when its distance to an event core is < R;
* the strip at exactly R (`strip-2`, the transmon's 64 segments) is a `SameConductorStrip` at
  `Separation / R = 1` on both sides and its corners are plain corners: no knife edge between
  1.95 / 2 / 2.05 um beyond the separation value itself.

**Knife-edge census (decision 82(4), 2026-09-25).** The strict rule is kept and the matching
radius is chosen off the common layout dimensions (the default of the identification tools is
R = 2.1 um: thresholds R 2.1, 2R 4.2, 10R 21 um; the library's `MatchingRadius` stays
authoritative and R = 2 um libraries are unchanged). So that any design can check its R, the
manifest reports `Identification.KnifeEdgeCensus`: the perimeter length with another perimeter
point (3D; the same chain beyond the self-pair neighbourhood; runs sharing a vertex excluded)
at a distance within `KnifeEdgeBandRelative` = 0.01 of R and of 2R, the chain length whose
windowed bend radius lies within 1 % of `StraightBendRadiusOverR` R, and the vertices whose
turn lies within 1 % of the corner threshold, each split into the below / above sides
(samples every 0.5 R). DS-SCT-001 at R = 2 um: 7,247 um within 1 % of R and 7,175 um within
1 % of 2R (the 2 / 2 / 2 um flux lines). The library continuity gate
(`coupon_library.py continuity`, `qualification-gates.json` LibraryContinuity: a pair / stack
model whose consecutive separations are all >= 2R (1 - 0.01) responds per edge and per unit
length like the isolated-edge model within 1 %) runs on every written process library.

## (e) Patch construction from the features (phase 4; solve path and patch dry run)

The three-dimensional correction patches are built from the feature list (`ResponseCorrection.
PatchConstruction = "Features"`, the default; `"Legacy"` keeps the former per-interface-group
classification for comparison only). Matching is the key-based pass of (a); then

* a feature **matched** by signature becomes its patches, built from its own portions,
  vertices and canonical frame (below); a feature **unmatched** (no model with its signature,
  or a type the library cannot model such as `UnclassifiedParallelPair`) is **omitted alone**:
  no interface group, chain, pair or neighbour loses its correction; under `UnmatchedPolicy =
  Error` the solve aborts with the count of unmatched features (never in the preflight);
* **excluded segments** carry no feature and are never corrected; analytic exclusion zones
  are outside every portion by construction of the assignment.

Patch per class (weights in mesh units; `CouponDepth` = the model's longitudinal depth):

| Feature | Patches | Frame (u, v, w) | Weight |
|---|---|---|---|
| `IsolatedEdge`, `CurvedEdge` | one per quadrature point of every portion (`2 x order` Gauss points) | u = gap direction, v = process normal of the segment | `(s1 - s0) x w_q / CouponDepth` |
| pairs (`SameConductorGap`, `DifferentConductorGap`, `SameConductorStrip`, `Curved*`) | quadrature on **both** sides, side factor 1 / 2 (the longitudinal measure is the mean of the two sides: exact for a straight pair, the centreline for concentric arcs); at a sample p its foot q on the partner's portions | origin (e1 + e2) / 2, u from the model's first edge e1 toward e2, v = mean process normal; the first edge is the lower side along the feature's lateral axis `Frame.Axes[1]` (the higher one for `Chirality` -1: the canonical orientation is the mirror) | `(s1 - s0) x w_q x 1/2 / CouponDepth` |
| `ParallelEdgeCluster` | quadrature on every side, side factor 1 / n; origin on the canonical first edge at the sample's longitudinal coordinate; anchors on the first edge of every conductor label | u = lateral axis toward increasing canonical offsets, v = mean process normal | `(s1 - s0) x w_q / n / CouponDepth` |
| `ConvexCorner`, `ConcaveCorner` (sharp or rounded), `Endpoint`, `Junction` | one patch at `Frame.Origin` (the vertex or the virtual corner of a fillet) | `Frame.Axes` (below) | model weight (1) |
| `SpatialEdgeCluster` | one patch | the model's canonical frame composed with the feature's: a model-frame point m maps to `F.origin + F.axes^T M.axes (m - M.origin)`, M from `CanonicalClusterSignature` of the model's stored edges (identity for a model keyed by its `Signature` alone, which is built in the canonical frame) | model weight (1) |

**Vertex-feature frames** (`Frame` of the manifest, shared by the library builder): corner:
x = the first arm away from the (virtual) corner, the arms ordered so that the second is
counterclockwise about the process normal (a corner is its own mirror image), y = n x x;
endpoint: x = the arm, y = +-(n x x) toward the gap; junction: x = the canonical first arm of
`CanonicalJunctionSignature`, y = +-(n x x) so that the canonical arm order proceeds
counterclockwise in (x, y). A legacy junction model (absolute `ArmAngles`) is mapped by its own
canonical order (first arm angle theta, orientation): u = cos(theta) D - sigma sin(theta) (n x D),
v = sigma (n x u), sigma = +1 when both orientations agree.

**Library contract.** A model keyed by its `Signature` (the feature's canonical object, `Type`
included; `Signature.Type` must equal `Topology`) needs no version-1 geometry parameters of its
own; a model matches every feature of its topology whose parameters lie within the signature
tolerance (decision 85(1): 1e-3 R for lengths, 1e-2 deg for angles; the nearest model wins,
ties by name), and the library builder groups the feature instances of one topology within
the tolerance into one coupon at their representative signature (the version-1 `Requirements`
record carries `Signature`, `Instances`, `ParameterSpread`, `ExactParameters`; the
`ParallelEdgeCluster` / `CurvedParallelEdgeCluster` `Geometry` carries `Edges`, `EdgeCount`
and, for the curved class, `BendRadius`); the curved classes (`CurvedEdge`, `CurvedSameConductorGap`, `CurvedDifferentConductorGap`,
`CurvedSameConductorStrip`) exist only as signature-keyed models and are patched like their
straight analogues along the curved portions; a cluster model keyed by its signature needs no
`Edges`. A model's interface types are part of its key (a model mapping MA + MS + SA never
matches an SA-only feature). Runtime models are one per (library model, target interfaces by
slot; slot k = the k-th distinct target map of the feature's portions in sorted order).

**Patch dry run.** `palace --surface-response-preflight` builds the same patches without a
field solve and writes `surface-response-patches.csv` next to the manifest: `Patch, Feature,
Topology, Model, ModelIndex, Weight, ModelWeight, QuadratureWeight, SideFactor, CouponDepth,
Segment, S0, S1, Origin, AxisU, AxisV, AxisW` (manifest units; `Segment` = the manifest
segment index, `[S0, S1)` the portion from the segment's canonical key origin; vertex and
cluster patches carry `Segment` -1 and `CouponDepth` 0). The audit's gates A7: the patched
feature set equals the matched set; every portion of a matched longitudinal feature is exactly
one quadrature interval (sum of quadrature x model weights = 1) and `Weight = ModelWeight x
QuadratureWeight x (S1 - S0) x SideFactor / CouponDepth` with `SideFactor` = 1 / claimed chains;
vertex / cluster features carry patches without a portion whose model weights sum to 1; no
patch on an unmatched feature, an excluded segment or an excluded portion. With the
signature-only library built from the manifest's own features
(`examples/surface_response_identification/signature_library.py`) the covered length equals
`Totals.AssignedLength`: the whole perimeter minus the recorded exclusions.

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
                  "VertexJoinsClusterOverR": 2, "VertexWindowOverR": 1,
                  "ParallelCosineTolerance": 1e-8, "RoundedCornerTangentTolerance": 0.05, "ArcFitToleranceRelative": 0.05,
                  "SignatureLengthQuantumOverR": 1e-6, "SignatureAngleQuantumDegrees": 1e-6,
                  "StraightBendRadiusOverR": 10, "CurvatureWindowOverR": 1,
                  "PairSeparationToleranceRelative": 0.05, "PairSeparationSamplesPerInterval": 16,
                  "PairSeparationEstimate": "per sample: chord reading C = min over the two chains of the maximum sampled closest-point distance within max(R, local chord) of the sample / its foot where the chain bends, R on straight runs; inscribed reading C / cos(turn / 2) with the larger local joint turn; interacting iff both < 2R; feature separation = mean C",
                  "PairConstancyWindowOverR": 1, "PairSampleSpacingOverR": 0.5,
                  "PairCandidateReachOverR": 2.1, "SamplingMargins": "...", "CrossLayerReachOverR": 2,
                  "SelfPairNeighbourhoodOverR": 3.14159, "StackRule": "...",
                  "StackCompositionCap": 64, "ClusterExtensionRule": "...",
                  "SignatureParameterToleranceOverR": 1e-3, "SignatureAngleToleranceDegrees": 1e-2,
                  "SignatureMatching": "...", "MutualSidesOverhangOverSeparation": 0.32,
                  "KnifeEdgeBandRelative": 0.01, "FacingGateExclusions": "...",
                  "PlaneRule": "...", "PortRule": "...",
                  "Comparison": "strict less on the quantized grid"},
  "ReferenceProcessNormal": [nx, ny, nz],
  "Features": [ {"Id": k, "Type": "...", "Signature": {...}, "Hash": "sha256", "Chirality": +-1,
                 "BendRadiusOverR": r | null, "ExactParameters": true | false,
                 "Length": L, "Portions": [[segment, s0, s1], ...], "Vertices": [v, ...],
                 "Frame": {"Origin": [...], "Axes": [[...],[...],[...]]},
                 "Match": {"Status": "Matched" | "Missing", "Model": "name", "Deviation": d} } ],
  "Segments":  [ {"Key": [[x0,y0,z0],[x1,y1,z1]], "Length": L, "Chain": c,
                  "Portions": [[s0, s1, feature], ...] } | {"Key": ..., "Length": L,
                  "Exclusion": {"Class": "...", "Reason": "..."}} ],
  "Vertices":  [ {"Point": [...], "Type": "Corner|Endpoint|Junction", "TurnDegrees": t,
                  "Feature": k} | {"Point": ..., "Class": "TruncationCut"} ],
  "Exclusions": [ {"Class": "...", "Reason": "...", "Count": n, "Length": L} ],
  "Totals": {"PerimeterLength": L, "AssignedLength": La, "ExcludedLength": Le},
  "Diagnostics": {"SamePriorityClaimOverlaps": 0, "Rule": "...",
                  "StackGeometricOffsetIntervals": 0, "StackCompositionCapHits": 0, "StackCompositionCap": 64,
                  "ClusterExtension": {"Passes": n, "AbsorbedPortions": n, "AbsorbedLength": L,
                                       "VertexFeaturesJoined": n, "Rule": "..."},
                  "StackEndThirdBodyLength": L},
  "KnifeEdgeCensus": {"BandRelative": 0.01, "SampleSpacingOverR": 0.5, "SampledLength": L,
                      "Distance": {"R": {"Below": L, "Above": L, "Total": L}, "2R": {...}},
                      "BendRadius": {"10R": {...}}, "CornerTurnDegrees": {"30": {...}}},
  "GeometryDigest": "sha256"
}
```

Lengths and coordinates are in mesh units on the `ScaleLength` grid of the version-1 manifest.
The audit reads `Identification` when present (A1 partition from `Segments`, vertex census from
`Vertices`, exclusions from `Exclusions`, A2 from the cluster frames, A3 / A5 from
`GeometryDigest` and the feature signatures) and falls back to the version-1 aggregate reading
otherwise.
