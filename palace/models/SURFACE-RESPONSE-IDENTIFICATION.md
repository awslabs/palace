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
2. **Vertex table.** Every perimeter vertex that is a joint (two segments, not collinear on
   the 1e-12 direction grid) and NOT noise under the geometric joint noise rule
   (`JointNoiseSagittaOverR` = 0.05, `metaledge.hpp`: `kJointNoiseSagittaOverRadius` /
   `JointIsNoise`; USER decision 121 (B), 2026-09-28, replacing the 1 deg angular threshold
   of decision 117(4) and the 30 deg corner-class threshold of decision 73: a joint turning
   by t between two straight pieces — collinear mesh segments merged — is a straight
   continuation when the sagitta (c / 2) tan(t / 4) it implies on the SHORTER adjacent piece
   c is below 0.05 R, item 7) and that no fitted arc absorbs (item 3, arc rule), and every
   endpoint / junction vertex that is not a
   simulation cut is mapped to exactly one vertex feature (corner / endpoint / junction) or to
   exactly one cluster (membership); a joint absorbed by an arc is a `RoundedCornerVertex` /
   `BendVertex` record of the arc. Vertices on the truncation boundary are cuts, not
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
(`physical_chain`, `metaledge.cpp`): a polyline of straight runs joined at sub-noise vertices
(the geometric joint noise rule, `JointNoiseSagittaOverR` = 0.05: implied sagitta on the
shorter adjacent piece below 0.05 R; items 1 and 7, the curved-edge chain rule), the
chains a fitted bend arc's absorbed corners separated being merged by the identification
(item 3) and a rounded corner's arc taking a chain of its own. Chain tangents are canonical: from the lexicographically
smaller endpoint.

1. **Through-vertex pairs.** Two points p, q on two chains that meet at a vertex v are
   *through-vertex* when both are within 2R of v. Such pairs are described by the vertex
   feature (or by the cluster the vertex belongs to), never by an interaction event. This is
   the existing convention (`ConnectedNearVertex`) restated on points instead of mesh
   segments so that it is mesh independent; it makes corners with interior angle >= 60 deg
   plain corners and lets the arms of sharper corners interact beyond 2R from the vertex
   (2a sin(theta/2) < 2R for a >= 2R has solutions iff theta < 60 deg).
   **Joint noise rule (USER decision 121 (B), 2026-09-28; `metaledge.hpp`
   `kJointNoiseSagittaOverRadius` = 0.05, `JointIsNoise`, mirrored by `perimeter.joint_is_noise`
   and the oracle's `_design_joint_is_noise`).** A two-segment vertex turning by t between two
   straight pieces (the maximal collinear runs on either side, so that a refinement midpoint
   or a second-order mid-edge node never changes the pieces: A5) is a straight continuation of
   its chain — a REGULAR joint whose turn feeds the windowed curvature (item 7) — iff the
   deviation from straight it implies at the resolution of the correction is below
   `JointNoiseSagittaOverR` x R = 0.05 R: the sagitta (c / 2) tan(t / 4) of a chord of length
   c = the SHORTER adjacent piece read as one chord of a circle turning t per chord
   (sagitta = rho (1 - cos(t / 2)) with c = 2 rho sin(t / 2)); compared on a 1e-9 relative
   grid. Every other joint is a corner vertex and a chain break unless the arc rule absorbs
   it. Why this quantity and not the vertex distance from the line through its neighbours:
   (i) it is the very quantity the arc rule records per chord (`Arcs[].MaxChordSagittaOverR`,
   the mesh-coarseness diagnostic at the same 0.05 R), so a sub-noise joint is a chord joint
   of a curve the correction cannot resolve — one constant for the noise threshold, the arc
   diagnostic, the extraction, the identification's end-joint test, the census and the
   audit mirror; (ii) sub-nm mesh slivers (the 1-2 nm flux-line segments of DS-SCT-002 whose
   roundoff directions turned by more than 1 deg and fragmented the stacks under the angular
   threshold, A8 2.1 mm) are bounded by c / 2 whatever their turn, a U-turn included
   (tan(pi / 4) = 1); (iii) spline steps of 1-6 deg on 1-5 um chords (0.5-2.6 R at R 1.9)
   imply 2-65 nm = 0.001-0.034 R and read as the smooth curves they discretise, where the
   vertex distance c sin(t) (4x larger at small t) would keep 6 deg joints on 1 um chords as
   174 deg corners; (iv) a real corner is never noise: 5 um arms turning 20 deg imply 0.23 R,
   the chips' longest chords (5.3 R) turning 4.3 deg imply 0.05 R (the knife edge, in the
   census as `JointNoiseSagittaOverR.0.05` on the shorter run). Consequences: the transmon's
   sampled fillets and the chips' 174-179 deg "corners" of the angular threshold are chains
   or arcs; a joint of 7.5 deg on 6 um leads at R 2 um implies 0.049 R (noise) — the census
   band records such cases. The extraction receives the threshold in mesh units
   (`MetalSurfaceExtraction::joint_noise_sagitta` = 0.05 x the library's matching radius; the
   post-processing edge tools pass their largest edge distance); the legacy per-group
   classifier (comparison only) keeps its angular 30 deg class
   (`corner_turn_tolerance_degrees`).
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
   **Arc rule — the CONCYCLICITY form (USER decisions 121 / 122, 2026-09-28; amends the
   sagitta form of 117(4) below and supersedes the joint-turn / piece-length form of decisions
   82(3) and 85(1), whose fixed points are kept where noted).** Along every perimeter path
   through vertices with exactly two path segments (regular or corner; a path stops at
   endpoints, junctions and cuts) the joints are the non-collinear vertices. From every
   unconsumed joint, the largest range of at least three following joints that turn the same
   way, **each by less than `ArcMaxJointTurnDegrees` = 50**, and total at most 180 deg is ONE
   arc **iff its joint vertices lie on one circle within the recorded fit tolerance
   `ArcFitToleranceOverR` = 1e-3 (the signature parameter tolerance)** — whatever the chord
   sagitta. "We do not want to artificially classify curves as corners": a coarsely meshed
   design curve is the curve it discretises (the transmon's 38.9 um CPW bends at 11-16 deg
   per chord, sagitta 0.09-0.20 R, are bends again; DS-OSC-003's 4 / 6-chord meanders are
   bends), and its coarseness is REPORTED: the largest chord sagitta rho (1 - cos(central
   angle / 2)) of every arc is recorded (`Arcs[].MaxChordSagittaOverR`, on the 1e-6 R
   signature grid) and the manifest's `MeshCoarsenessWarning` lists the arcs whose recorded
   sagitta is at or above `SagittaOverR` = 0.05 R (count, arc length, worst sagitta, arc
   indices; the identification log prints the same) — a mesh-coarseness diagnostic, not a
   membership test (decision-117 review m2: the former strict-below membership test is gone;
   the diagnostic's sense is "at or above 0.05 on the recorded grid"). The joint-turn cap
   keeps regular polygons corners (a square turns 90 deg per joint, a hexagon 60; an octagon
   at 45 deg per joint is a circle at the resolution of the correction — recorded knife edge,
   in the census as `ArcMaxJointTurnDegrees.50`). The rule applies to bends (radius >= R) and
   rounded corners (radius < R) alike, with no piece-length rule. The circle is (a) the one
   tangent to both arms at the end joints (centre at Ta + rho na, rho from the tangent
   lengths; antiparallel arms: half their separation) when every joint — the last one
   included, which unequal tangent lengths put off the circle — lies on it within the fit
   tolerance; else (b), **for a radius >= R over at least FOUR joints** (decision-117 review
   m3, revisited under the concyclicity rule and KEPT: three points are concyclic whatever
   they are, so without the sagitta test a 3-joint least-squares "bend" would absorb ANY three
   consecutive same-sign joints — a lead-end corner, an arc's first joint and its neighbour
   passed as a bend of radius 350 um under the former rule; with three joints only the
   tangent construction, whose circle is fixed by the arms, is a test), the algebraic
   least-squares circle of the joints (exact for an inscribed polyline, a set function of the
   joints), with each arm meeting the circle's tangent at its end joint within the geometric
   joint noise rule on the shorter of the arm piece and the first chord (a spline piece
   joining a bend tangent-continuously; read from the fitted circle, whose tangent moves by
   1e-3 R / rho for a 1e-3 R position error) or lying on the circle as a chord (an arc
   starting at a corner that lies on its circle — the sharp end of a rounded slot: the corner
   stays a corner, the arc starts at its next joint; a lead meeting a circular arc at a 90 deg
   corner is neither, and the corner is never absorbed). **Collinear-subdivision invariance
   (supervisor decisions 212 / 213, 2026-10-02; VALIDATION-PLAN (h)-8).** The arm whose far
   vertex the chord test reads is the arm's straight PIECE — from the end joint to the
   neighbouring joint or the path end, the rigid-run joint — never the adjacent mesh
   segment: inserting collinear vertices on a chord is a geometric no-op and the reading must
   not depend on it. Until 2026-10-02 the test read the far end of the mesh segment, so a
   coarse round pad whose lead attaches through joints above the turn cap (DS-CTX-003 C4:
   r 23 um, eight 43-deg chords, sagitta 0.86 R; the ground's r 36 hole) was a bend over
   its interior joints on the chip mesh (one mesh edge per design chord: the attach joint
   is on the circle) and sharp corners on the thin window mesh subdividing the same chords
   at 4 um (the first sub-vertex 0.8 um off the circle): 168 um of E1 class D over S4 / C3 /
   C4 / O1 / O3 / O4. Every other quantity of the arc test (joints, pieces, the every-point
   clause, the tangent construction) was already a rigid-run quantity: a joint is a
   non-collinear vertex at the direction quantum (1e-12 on the cosine), the same test that
   bounds the rigid runs and the extraction's pieces. **First-joint absorbability (decision
   213).** A chord arm at the range's FIRST joint whose far joint p is itself absorbable —
   unconsumed, the same sign, below the turn cap, and its own arm consistent with the circle
   (tangent-noise or a chord; p lies on the circle by the chord-arm clause itself) — is no
   arm: the range is not maximal at its start and the fit is refused; the scan reaches the
   arc's real start later (cyclically on a loop). This ends the recorded exposure of the
   closed-loop start rule below for bends of four or more joints: when the loop's longest
   piece is a chord INSIDE such a bend (short leads, long chords) the scan started inside it
   and, on a joint-only mesh, accepted the sub-range anchored at the loop start (its first
   arm, the previous chord, ends on the circle) and chopped the arc there, while on a
   subdivided mesh the sub-vertex failed the former mesh-segment test and the whole arc was
   found from its real start — the two discretisations disagreed, and the two-direction scan
   rescued the joint-only reading only where one direction's sub-range had fewer than four
   joints. A corner on the circle whose own arm is neither tangent nor a chord (a lead
   meeting a round pad at a sub-cap angle) is not absorbable: the arc still starts at its
   next joint and the corner stays a corner. The LAST joint's chord-arm clause is unchanged:
   a polyline turning more than 180 deg (the C4 pad) is still cut into <= 180-deg bends
   there under the current rule — which reading is physically right for such coarse polygons
   (corners or arc) stays the deferred rule question of decision 184 / D6. Tests: unit case
   `SurfaceResponseIdentificationCollinearSubdivision` (coarse pads of 10-60 um chords and
   136-174 deg joints joint-only vs subdivided at 4 um, irregularly and by the mesher; straight
   edges and sharp corners with subdivided arms; the start-rule loops incl. one whose longest
   chord lies mid-arc, with every start vertex), the python mirror's
   `test_arc_groups_collinear_subdivision_and_start_rule`, and the generic
   "collinear-subdivision invariance" gate (`subdivision_gate.py`: every synthetic layout of
   the existing gates re-identified with every mesh edge subdivided at two spacings and every
   closed loop with a rotated start vertex -> identical GeometryDigest / features). A rounded corner keeps the tangent
   construction only (its arms are tangent by construction and its site — virtual corner, arm
   directions — has no meaning otherwise). **Every-point clause (USER decision 184 (2),
   2026-10-01; the E8-3 / E8-4 false arcs of stage 0).** The concyclicity holds at EVERY
   point of the range, not only at its joint vertices: the interior points of the chord
   between two consecutive joints lie off the circle by the chord's sagitta on it, which is
   admissible only where it is below the joint noise resolution (`JointNoiseSagittaOverR` R
   = 0.05 R, the deviation from straight the noise rule itself cannot resolve: a finely
   chorded curve) or where one of the chord's end joints is a REAL joint (not noise under
   the geometric rule on its shorter piece: the design chord of a coarsely discretised
   curve, which decision 122 reads as the curve whatever its sagitta). ONE real end joint
   suffices (not both; recorded choice, 2026-10-02): a concyclic long chord between one real
   joint and one noise joint stays an arc chord — a straight lead leaving a circle is tangent
   to it, not a chord, so the configuration is geometrically unlikely, and the choice is the
   one consistent with decision 122's coarse-polyline exception. A chord deviating
   from the circle by more than the resolution between two NOISE joints is a straight
   edge and its range is no arc (`ChordsResolved`, applied to the tangent circle, the
   least-squares bend and the closed circle alike; the python audit mirror `perimeter.
   arc_groups` carries the same clause since 2026-10-02 — without it the audit's A1 vertex
   census still fitted the false bends and, at the DS-CTX-003 loop ends, read a 2-chord
   1.0 um fillet into one, 8,486 vs the classifier's 8,490 rounded corners): a 5,500 um straight trace edge whose
   ends carry two 1.2 / 2.4 deg joints on 12 um pieces — the first chords of the bends the
   run enters, exactly concyclic by mirror symmetry — read as a 132 mm bend bowing 15 R off
   the metal, and the arc-aware cluster geometry found that bow within 2R of the parallel
   ground edge 16 um away (DS-OSC-003: 3.9 mm of false 2-edge clusters; DS-SCT-002: 2.1 mm).
   For an inscribed polyline of equal chords the joint's implied sagitta and the chord
   sagitta are one quantity, so a uniformly chorded arc passes iff its joints are real or
   its chords are below the resolution: no new threshold enters. On a closed path the
   piece lengths wrap at every straddling joint pair (`CyclicGap`; formerly only the pair of
   the list's first and last joints wrapped, so the pieces of the pair straddling the
   path's start read negative and the noise rule accepted any turn there). Every joint that no arc absorbs and that is not
   noise under the joint noise rule (item 1) is a corner feature — a two-joint chamfer, a
   square strip end (the diameter chord of a semicircle cannot be told from a one-chord
   arc), a polygon of 50 deg joints, the 90 deg lead-end corner of a rounded slot. Recorded
   consequences: a spline-discretised bend whose vertices lie on one circle within 1e-3 R
   only over a few joints is chopped into exact-fit arcs (as before) and its leftover joints
   are corners only when they are not noise (a 1-6 deg spline step on chords below ~5 R is
   noise: the 174-179 deg corner classes of the angular threshold are gone); a mitred offset
   polyline of an arc (vertices at h / cos(step / 2) from the centreline, end vertices at
   h / cos(step / 4)) is an arc where its joints lie on the circle tangent to both leads
   (a two-chord mitre always; wider bars at coarser steps read corners at the polyline ends
   where a CAD offset of an arc path — concentric arcs at the same angles — reads one bend).
   The fillet gate and the arc-cluster gate hold over the discretisations whose joints all
   turn less than the cap (`Fillet.Resolved` / `ArcCluster.Resolved` in `synthetic_layouts`:
   a 2-chord 135 / 180 deg fillet turns 67.5 / 90 deg at its middle joint and is corners, a
   4-chord 180 deg fillet at 45 deg is an arc; a closed circle needs more than 360 / 50
   chords: an octagon), the coarse chords among them carry the diagnostic (`Coarse`); the
   polygon cases are reported with their corner readings, not gated.
   **Sagitta form of USER decision 117(4) (2026-09-28, amended above):** membership required
   every chord's sagitta strictly below `SagittaOverR` R on the 1e-6 R grid and no joint-turn
   threshold; it read the transmon's CPW bends as 387 corners of 164-175 deg (no library
   coupons), DS-OSC-003's meanders as 135 / 150 deg corners and, with its 1 deg angular noise
   threshold, the chips' spline joints as 24k / 15k corners of 150-180 deg (840 / 763 distinct
   angles) and sub-nm slivers as stack-fragmenting corners (user-decisions-117-sagitta):
   the reason for decisions 121 / 122.
   **Former form (decision 82(3), 2026-09-25; replaced the run-based rounded-corner rule of
   item 4 and the half-chord curvature of fitted arcs in item 7):** the largest range of at
   least three following joints connected by pieces shorter than 2R, turning the same way,
   totalling at most 180 deg and fitting ONE circle with the two arm tangents
   (`ArcFitToleranceRelative` = 0.05: tangent lengths equal within 5 %, every joint within
   5 % of the radius) was an arc; a closed loop starts its scan after its longest piece
   (kept). **Amendments of decision 85(1) (2026-09-26):** (i) a piece of 2R or longer
   stopped the arc only when either of its joints was a corner-class turn (> 30 deg), a long
   piece between sub-corner joints being a chord of a smooth polyline bend (DS-SCT-001's
   367 / 371 / 373 um route bends meshed with 4 um chords = exactly 2R at R = 2 um), and in
   that long-chord regime the scan stopped at the first failed fit — superseded by the
   sagitta cap; (ii) a BEND's circle (radius >= R) is the algebraic least-squares fit of its
   joint vertices where the arms are not tangent (0.5 % radius error and a 0.8 um centre
   offset of the tangent construction on DS-SCT-001's 370 um bends), a rounded corner keeps
   the tangent-length radius — kept, with the four-joint minimum and the end-joint rule
   above; (iii) a closed path whose joints all turn one way through 360 deg on one circle (a
   round pad, hole or via) is ONE arc of total turn 2 pi and a bend of exact radius whatever
   its radius: a circle of radius < R is one `CurvedEdge` with `RadiusOverR` < 1, never two
   180 deg "rounded corners" split at a numbering-dependent joint (review m6) — kept; the
   former "every joint sub-corner" condition (a square hole is a 4-corner polygon) is
   replaced by the sagitta cap (a square hole with a side above 0.24 R is four corners);
   (iv) both traversal directions of every path are scanned and the arc set absorbing more
   joints is applied (ties: FEWER arcs, then the smaller serialisation of the set, per arc
   (radius / R, turn, joint count, centre distance from the path's vertex centroid / R,
   first joint's arc-length distance from the nearer path end / R) on the signature grid
   — radius, turn, joint count and centre distance are invariant under translation,
   rotation, mirroring and path reversal; the first joint's position is read in the scan
   direction, so the per-arc key is not a set function of the geometric arc, but the PAIR of
   serialisations {forward, backward} is orientation invariant (a mirrored or reversed path
   scanned the other way gives key(forward') = key(backward)) and the comparison picks the
   congruent set whatever the input orientation (restated, review fix-4 m-C); the
   former serialisation of absolute centre coordinates was not, and on DS-SCT-001 mirrored
   it chose the other chopping of a spline route bend into exact-fit arcs, which moved a
   curved-class boundary by 0.036 um, a cluster core by 0.12 um and a cluster boundary by
   0.03 um (review fix-3 m2; fixed 2026-09-26, the mirrored manifest is now byte-identical
   apart from the chirality signs)) — the greedy scan from the first unconsumed joint
   depended on the input orientation, and a mirrored mesh gives the mirrored arcs. Recorded
   limitation: a spline-discretised bend (vertices within the 1e-3 R exact-fit tolerance of
   a circle only over ~10 joints) is chopped into consecutive exact-fit arcs whose cut
   positions depend on the scan direction; the tie-break makes the choice orientation-free
   but the chopping itself is a property of the polyline, not of the design curve. Only a
   path whose two direction scans give congruent arc sets with different joints (a path
   that is its own mirror image along its length) remains ambiguous. **Closed-loop start
   rule (recorded, review fix-4 m-C):** a closed path has no ends; both direction scans start
   at the joint following the path's LONGEST piece (in the scan direction), so that no arc is
   split by the loop start (a fillet, rounded corner or bend never contains the longest piece
   of a loop with straight sides; a loop that is one circle is the single 360 deg arc of
   (iii)); when several pieces tie for the longest, the first from the path's seed vertex is
   taken — the seed is the canonical (coordinate-sorted) numbering's first segment of the
   loop and therefore coordinate-dependent — and the tie-break's first-joint position is 0 for
   every arc of a loop. A different start point of a loop is not covered by the two-direction
   scan; the exposure was a loop whose tied longest pieces lie inside arcs (all-equal chords of
   a polygon are corners or one circle) — closed for bends of four or more joints by the
   first-joint absorbability clause above (decision 213), which refuses a sub-arc anchored at
   the loop start whatever the start vertex. Gate: the mirror gate's rotation variant identifies
   closed filleted loops rotated by 37 deg (and re-numbered) and requires identical content.
   Former recorded
   constants of the arc rule, retired by the sagitta form (`ArcFitToleranceRelative` = 0.05,
   `ArcFitAbsoluteToleranceOverR` = 0.05, `ArcInscribedAngleToleranceRelative` = 0.05,
   `ArcTangentLengthPreferenceOverTolerance` = 10, `RoundedCornerTangentTolerance` = 0.05):
   a joint of an arc lay within min(5 % of the radius, 0.05 R) of the circle (the relative
   tolerance alone let a 200 um route bend absorb an adjoining spline joint 10 um off its
   circle); in the long-chord regime every interior joint's turn had to equal half the sum of
   the central angles of its two chords within 5 % (implied, to the slop of the fit
   tolerance, by joints on one circle turning one way — and far stricter than that slop on
   short chords, where it chopped exact 1 deg polylines: retired); a bend kept the
   tangent-length circle when every joint lay on it within 10x the parameter tolerance, else
   the least-squares circle (now: the tangent circle when the joints lie on it within the fit
   tolerance itself, else the least-squares bend). Gates (`permute_msh2.py`): DS-SCT-001 and a synthetic stack
   renumbered (seeded node / element permutation) give identical manifest content; mirrored
   in x, DS-SCT-001 and the stack suite give identical content with every
   `SpatialEdgeCluster` chirality negated (the mirror gate, decision 88); rotated by 37 deg
   about the plan-view normal (`--rotate-degrees`), the 42 closed filleted loops of the
   fillet gate give identical coordinate-free content and digest (the closed-loop rotation
   variant of the mirror gate, review fix-4 m-C).
   **Self-pairing point-wise (2026-09-26):** the partner search of a chain facing itself and
   the self events exclude the part of a run within pi R of arc length of the POINT (the
   former run-level exclusion dropped a whole 28 um leg of a hairpin for every point of its
   fold, so that the fold end facing the far leg at 1.6R had no partner and no event; it is
   a cluster now). Two-joint polylines are never arcs (a chamfer, and a square strip end
   — the diameter chord of a semicircle — cannot be told from a one-chord arc: they stay
   corners). An arc of radius < R (tangent arms; any total turn — a rounded corner of small
   turn is a corner of large angle) is ONE rounded corner: `ConvexCorner` / `ConcaveCorner` by the side of the centre (convex when
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
   **Curve-aware cluster geometry (option A, decision 91(1), 2026-09-27; fix-4 residual 1).**
   Every run on a fitted arc (a rounded corner's or a bend's chord) carries the piece of the
   fitted circle between its end angles (`BuildRunArcGeometry`; the run parameter maps
   linearly onto the angle), and the cluster machinery evaluates it ON THE ARC: the event
   sublevel sets (`CoresOnRun`), the core-core merge, the vertex join, the radius-R claims
   (`RunIntervalWithinPiece`), the through zones it uses, the extension's across rule (the
   domain of an arc piece is the radial projection inside its angular range; the wedge turn at
   a piece end is read between arc tangents, so two chords of one circle or a chord meeting
   its tangent arm turn by nothing) and the stack-end images of cluster / window ends on
   concentric partner arcs (`ArcAwareImage`: a taken end on a bend cuts the concentric partner
   radially). The distance along an arc is not convex: its sublevel sets are bracketed by
   samples at most `ArcSampleSpacingOverR` = 0.05 R apart (at least 8 per piece) and the
   crossings bisected (`SampledSublevelIntervals`), so the result depends on the geometry alone
   to the bisection precision; two chords keep the former exact formulas, so a geometry without
   arcs in its clusters is unchanged. The cluster's claims on one arc are ONE signature portion
   (`ClusterSignaturePortions`), serialised in the frame as `{P: sorted ends, Arc: [centre,
   midpoint], GapRadial: +1 metal inside the circle / -1 outside}` (a closed circle claimed
   whole has no ends of its own: it is serialised as centre + radius, with `P` the circle's
   point on the frame's +x side, repeated, and the midpoint at its antipode, so the mesh's
   first joint never enters); the frame candidates add the arc's end tangents and end radial
   directions (none for a closed circle; a cluster of closed circles only takes the
   centroid-to-centre directions, or a fixed in-plane axis when they are concentric), the
   origin is the length-weighted centroid with the arcs' analytic centroids, `EdgeCount`
   counts an arc portion once. Known behaviour of the arc fit (recorded from the synthetic
   arc set): two tangent fillets joined by a straight end edge shorter than the arc-fit
   tolerance (e.g. two 0.75 R fillets of a 1.5 R + 0.3 R finger end) lie on one circle within
   the tolerance and are absorbed into ONE semicircular arc — a straight edge shorter than the
   fit tolerance between tangent fillets does not survive as a straight portion (the synthetic
   set uses 0.25 / 0.5 / 0.7 R fillets for that reason). The manifest
   records the fitted arcs (`Arcs`: centre, radius, turn, kind RoundedCorner | Bend, joints,
   segments) and every chord segment's arc (`Segments[].Arc`); the coupon builder chords a
   signature arc at `ClusterArcChordStepDegrees` = 5 deg, finer so that no chord exceeds
   `ClusterArcChordMaxLengthOverR` = 0.25 R (`signature_library.cluster_plan_view_edges`), so
   the coupon geometry is a function of the signature and not of the device mesh. The facing
   gate reads hits involving arc segments on the arcs too (`facing_check.ArcGeometry`). Result:
   the re-mesh gate 74 / 74 (fix-4: 71 / 74; hairpin 12.08-13.40 um and u-ring 73.12-73.36 um
   were the chord-based extents), 222 / 222 digests identical across the mesh variants; the
   hairpin whose strip is exactly 2R wide (`hairpin-rho1p5-g1`) reads a rounded corner plus a
   curved edge on every mesh (the concentric design circles are 2R apart: no event, the strict
   rule), where the inscribed chords used to dip below 2R and form a mesh-dependent cluster
   (the recorded knife edge, resolved on the design geometry); DS-SCT-001's three CPW
   termination clusters go from 337 / 337 / 379 chord edges (three hashes) to 29 / 29 / 33
   edges with 13-14 arc portions (two of them one hash: the same design box).
4. **Vertex features.** A corner (2 chains, a joint that is not noise under the geometric
   joint noise rule `JointNoiseSagittaOverR` = 0.05 and not absorbed by a fitted arc), endpoint (1 chain) or junction
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
   order cluster > vertex window > translational. **Sliver rule (supervisor decision 222,
   2026-10-02; `Conventions.SliverRule`, `PortionMinimumLengthOverR`):** NO portion shorter
   than the signature parameter tolerance `SignatureParameterToleranceOverR` x R = 1e-3 R
   exists, a portion being a maximal contiguous stretch of one feature side along a chain
   (across the chain's runs and mesh segments: a sub-tolerance MESH SEGMENT inside a long edge
   continues that edge's portion and is not one). A shorter stretch is roundoff between two
   claim boundaries — a cluster ball cutting a pair piece, a claim ending next to a run end,
   the perpendicular foot of a neighbour's cut end on a near-parallel member — that no
   signature parameter resolves (every length is matched within the tolerance) and no coupon
   models: it joins its adjacent portion on the chain, the LONGER of its two neighbours (ties:
   the one before it along the chain), taking that neighbour's feature and side; the join is
   decided by the chain order and the neighbours' lengths, never by a feature id. A pair /
   stack side whose surviving length is below the tolerance is no side (the feature
   dissolves, its pieces return to the run). A stretch with no adjacent portion on its chain
   (a whole chain, or a piece between two `CrossLayer` zones, shorter than the tolerance) has
   nothing to join and stays, counted. Reported as `Diagnostics.SubTolerancePortionsJoined`
   (`Count`, `Length`, `MaxLength`, `Isolated`). The cluster signatures (`EmitClusters`,
   the `Portions` of a cluster's key) are emitted BEFORE this rule runs in `Assign`, so a
   cluster's signature and its assigned portions (hence the spatial patch `Claims`) can
   differ by one sub-1e-3 R sliver joined to or from the cluster — within the matching
   tolerance by construction, and no worse than the former 1e-6 R join. Why the tolerance
   and not the signature grid:
   the former rule joined only pieces at or below `SignatureLengthQuantumOverR` x R = 1e-6 R
   (DS-SCT-001's three `CurvedSameConductorStrip` records of 1.6e-7 um in total), a knife-edge
   the stage-1 S1p window missed by 2.3 %: its 4-edge flux-line stack's claim on two members
   started at s = 1.9435e-6 um = 2 um x the 9.7e-7 rad near-parallel tilt of the window set
   (the foot of the neighbour edge's cut end at the y = -200 truncation cut), the 1.0229e-6 R
   remainder on one member belonged to the 3-edge stack whose real end lay 64 um away (a
   sample placed on it pointed its lateral axis along the segment and the
   `ParallelEdgeCluster` placement failed closed, `LongitudinalCellOffsets` |cos| 0.0627) and
   its twin read as a 2e-6 um `IsolatedEdge`; the same class sat in the E1 windows S1 / S1p /
   S4 and a 3.9e-4 um cluster portion of CTX C1. The tolerance is the length below which the
   signature contract already treats two readings as one ("Signature tolerance and library
   grouping" under item 7), so nothing a
   coupon could distinguish is merged; the placement's `MFEM_VERIFY` stays as the fail-closed
   net. Unit test `SurfaceResponseIdentificationSubTolerancePortions` (a near-parallel stack
   whose members end on a truncation cut with a 1.5e-6 rad tilt; a 1 nm mesh segment inside a
   long edge).
7. **Curved-edge chain rule** (decision 73(1); phase 2). A chain is defined by its
   significant vertices only: collinear splits (refinement midpoints, second-order mid-edge
   nodes) merge into one run, and the joints between runs that remain inside a chain — the
   joints of fitted bend arcs and the sub-noise joints (implied sagitta on the shorter
   adjacent piece below `JointNoiseSagittaOverR` R = 0.05 R; USER decision 121 (B): every
   other joint is a corner and a chain break) — are the
   polyline's bends. The turn at every sub-noise joint (the joints of a fitted arc are
   accounted for by the arc: a rounded corner's turn, a bend's exact density 1 / radius) is
   spread over the two adjacent half-chords — the
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
     over its portions; not hashed; null on straight chains) and, on such a feature, every
     portion's **signed windowed turn toward the metal** (`PortionTurns`, radians: the
     integral over the portion of the *signed* windowed curvature, positive where the edge
     bends around its metal — a disk-like edge — negative around its gap; not hashed; the
     sum `TurnTowardMetal` for a one-sided feature). The turn is the weight of the
     first-order curvature term (design (e): C4, decision 108);
   * an unpaired curved portion is a `CurvedEdge` feature (one per curved chain section)
     with `RadiusOverR` = the section's tightest windowed radius and `Convexity` in the
     signature: `Convex` when the signed windowed curvature over the section bends around
     the metal (the edge of a disk; the curved coupon of a disk edge models it), `Concave`
     around the gap (a hole edge), `Mixed` when both senses reach the curved regime inside
     one section (an S-bend whose two arcs lie within the window of each other: reported
     in the signature and never read as one convexity — such a feature is unmatched with
     the reason; a straight-like wobble of the minority sense is not a bend of the class).
     Recorded knife edge (review J/C4 m4): a section is `CurvedEdge` on the UNSIGNED
     windowed curvature while `Mixed` needs BOTH signed senses in the curved regime, so a
     wobble whose two senses each stay below the threshold while their unsigned sum crosses
     it (two opposite joints of ~5.3-5.7 deg within one window W = R, or three alternating
     ~4 deg joints) is read as one convexity at the unsigned kappa and corrected by that
     family; the band is ~0.5 deg wide and above the ~1 deg noise floor of real polylines.
     The sign convention is the coupons': the curvature family (`Kappa`, `Convexity`
     records) of a disk edge is `Convex`;
   * a fillet at a corner (radius < R; the run-based rounded-corner rule of item 4, computed
     from the arms' accumulated turn and the tangent distances, hence refinement-invariant)
     is a rounded corner and takes no part in the curvature;
   * **pairs along bends** (phase 3: local constancy): two chains that are not two
     parallel straight runs (those keep the translational rule of item 2; **one parallel
     relation for every rule**, USER decision 184 (1), 2026-10-01: two rigid runs of one
     plane are parallel iff |cos| of their tangents exceeds 1 - `ParallelCosineTolerance`
     = 1 - 1e-8, i.e. within 1.4e-4 rad — dimensionless, hence independent of the length
     unit and of R — and the *parallel classes* are the connected components of that
     relation per plane, built from the runs' 1e-9 `DirectionKey` buckets sorted by their
     in-plane angle with consecutive buckets within the tolerance linked, the 0 / pi wrap
     included: a canonical, numbering-independent partition (`BuildDirectionClasses`,
     `ParallelRigidRuns`). Parallel rigid runs pair and stack by the translational rule
     ONLY; every other rigid pair goes through the bent-pair and event rules. Before this
     rule the translational stage grouped the runs on the 1e-9 grid while the two other
     stages skipped rigid pairs within the cosine tolerance, so runs tilted by between
     1e-9 and 1.4e-4 rad — a mesh-noise wobble of 1.8e-4 um over a 186 um chip edge next
     to an exactly axis-aligned partner — were paired by NEITHER: the DS-SCT-002 flux
     lines read as a 3-edge stack plus an isolated fourth edge facing it at 1.05 R over
     12 x 175.7 um (stage-0 E8-1)) are sampled over
     their candidate facing regions — the points within `PairCandidateReachOverR` = 2 (1 + 0.05)
     R of the other chain, not beyond either end of it (the half-plane past a chain end along
     its outward tangent), outside the 2R zones of shared vertices, cut at the curved
     boundaries, at most `PairSampleSpacingOverR` = 0.5 R apart. A sample is **locally
     constant** when the sampled distances within `PairConstancyWindowOverR` = 1 R of it along
     its own chain vary by at most `PairSeparationToleranceRelative` = 0.05 of their minimum (a
     polyline of sub-noise turns at constant width varies by less than 1 / cos(0.5 deg) - 1, a
     bend arc's chords of sagitta <= 0.05 R at constant width by up to 0.05 R / separation =
     2.5 % at 2R; the pair response sensitivity d dR/dd is O(1)). The constant portions are the pair;
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
     sampled distance within a window of half-width max(R, min(local chord,
     `PairChordWindowCapOverR` = 2 x R)) where the chain bends within `PairBendProximityOverR`
     = 1 x R of the sample / the foot along the chain (a joint with a turn or a fitted bend
     arc: the window holds a vertex of an inscribed polyline and a full chord of an offset
     polyline) and R elsewhere (no dip; a taper is read locally), about the sample on its own
     chain and about its foot on the other chain. Both constants are dimensionless multiples
     of R. The locality is NOT the windowed curvature (implementation corrected 2026-10-01,
     USER decision 184 (3) as approved): the curvature rule spreads every joint's turn over
     its two half-runs, so a 0.1 deg taper kink between a 200 um lead and a 400 um taper run
     makes kappa > 0 over 100 um of the lead; the former "bends where kappa > 0" test with a
     window of the whole run length read the taper from the lead (the DS-CTX-003 532 um
     1 / 2 / 1 um stacks read 2.2 / 2.7 / 2.0 um once the sub-piece rule below stopped
     averaging the exact samples only; the straight leads of the re-mesh gate's 10 R curved
     stacks were mesh-dependent). The half-run spreading itself is unchanged: it is the
     curvature rule's density (the curved class, the bend radius annotation, the signed turn
     of a portion, the BendRadius census band). That chord
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
     windowed radius of the two sides (the inner side of concentric arcs) and `Convexity` =
     that of the signature's FIRST edge (side 0 for chirality >= 0, the last side for -1:
     the edge the library model's first edge e1 is placed on), read on that side's pieces
     or, when it carries no curvature of its own, as the opposite of the far side's
     (concentric edges bend in opposite senses relative to their gaps: for a gap the inner
     edge is convex, for a strip the inner edge is concave; for a cross-section that is its
     own mirror image, chirality 0, side 0 is the outermost edge of the bend — decision 214
     (i) — so the field is DERIVED there: a symmetric curved gap reads Concave, a symmetric
     curved strip Convex, a symmetric curved stack the convexity of its outermost edge, Convex
     when that edge's `GapSide` is -1 (gap away from the stack) and Concave when +1; the field
     is kept in every curved signature, chirality 0 included, for one schema, and carries
     independent information for chiral cross-sections only, whose side 0 is fixed by the
     edge pattern whatever the bend sense); curved
     stacks (`CurvedParallelEdgeCluster`) record it the same way, and when neither outer side of
     a k > 2 stack carries curvature (the bent side is an interior one) the first interior
     side that does is read, related to the first edge through the signature's `GapSide`s
     (same convexity when the gaps point the same way, opposite otherwise; review J/C4 m5 —
     formerly an abort). The cross-chord
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
   for every asymmetric cross-section; a symmetric one, chirality 0, has no observable
   orientation when straight and puts the member with the smallest (chain, position) on side
   0 — but a CURVED symmetric cross-section has one, since its `Convexity` is read on side 0
   and the chain ids follow the coordinate-sorted canonical numbering: that key made the
   Convexity of the two- and three-trace curved stacks a function of the frame (decision 214
   (i), found by the subdivision gate's rotation variant: Convex <-> Concave at 37 deg; no
   chip carries a symmetric curved pair or stack). The bend sense orients a curved symmetric
   cross-section instead: the lateral axis points toward the centre of curvature (the sum
   over the members of the signed windowed curvature toward the metal times the metal side
   along the lateral axis), so side 0 is the OUTERMOST edge of the bend — the key decides
   only where that sum vanishes; a rotation, a translation and any numbering of the perimeter
   then read the same signature, and the Frame's lateral axis of such a feature is canonical
   too); one feature per (type, signature, curvature class,
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
   distinct chains); the rule applies to non-circular folds of sub-noise joints (no arc
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
   is gone (a facing test; the separation is max(the feature's mean, the piece's LOCAL
   separation to the partner chains) — the local separation samples the piece's two ends and
   its midpoint only, so it is a lower bound of the piece's true maximum separation, exact
   for a monotone taper (the recorded case) and for a constant pair: a slow taper's wider end
   is still facing — with the
   mean alone DS-CTX-003's flux-line launchers, 1 -> 6 um over 560 um, lost the inner gap's
   pair where the outer links ended and read four isolated edges at 3.7 um, decision 103 /
   review m11; the synthetic `stack-taper-through-2R` is the reproducer); the assembly
   accepts a foot within the reach only along a link already decided interacting. **Facing gates (audit A8, `facing_check.py`):** the length of
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
   within `PairBendProximityOverR` = 1 x R of the sample or its foot along its chain — a
   joint with a turn or a fitted bend arc, not the windowed curvature — (the leads of a
   route: constant and exact for the design lines; a chord of a coarse polyline bend has a
   joint within R and is not read this way — aligned inscribed chords of a 370 um bend are
   31 nm inside the circle), or
   (b) the radius difference of two fitted bend arcs whose centres coincide within the
   parameter tolerance (the centre offset bounds the error of the difference; two local
   circles of a spline fitted piecewise are not concentric). **Separation of a sub-piece and
   of a link (decision 85(1) as AMENDED by USER decision 184 (3), 2026-10-01; the stage-0
   E8-7 finding: DS-CTX-003 launcher strips 3.68-3.79 um wide keyed as the 2.0 um strip,
   `ExactParameters` true, over 36 x ~46 um).** A piece's separation reflects its OWN
   geometry. The straight reading (a) is exact at a sample only where the sample's foot is
   the perpendicular projection onto the partner run (closest-point distance = line
   distance within the parameter tolerance; a foot clamped at a run end, opposite a
   diverging piece, is no exact reading). A sub-piece (the consecutive samples of one
   constancy / interaction status along a run) is EXACT only when it has an exact sample
   and every sample's OWN local separation — its exact reading where it has one, else, away
   from any bend within `PairBendProximityOverR` R of the sample / its foot, its
   closest-point distance to the partner chain — lies within the parameter tolerance
   (1e-3 R) of the exact mean, which it then reads. The chord-read samples next to a bend or
   a taper kink are not judged: they have no exact reading of their own, their window-max
   chord reading (reaching up to `PairChordWindowCapOverR` R past the joint) holds the
   neighbouring piece's distances, and the end sample of a lead at the joint reads the
   partner's first chord (2.85 cos(2.5 deg) for a fine bend: 1.4e-3 R off the line) — they
   are the junction's geometry, read off the chords like the bend. The chord reading stays
   the pair-level VALUE estimator for polyline bends; no new constant (the proximity is
   `PairBendProximityOverR`). Otherwise the sub-piece reads the length-weighted mean of the
   chord readings over the sub-piece (equally spaced samples; unchanged by the 2026-10-01
   correction) and is not exact.
   **Exact stretches (USER decision 203, 2026-10-02; the review's MAJOR-1: exactness
   all-or-nothing per run).** The exactness is LOCAL along the run, like the readings. Within
   a sub-piece of one constancy / interaction status, an EXACT STRETCH is a maximal
   contiguous stretch of samples that HAVE an exact reading (the straight reading (a) with
   its perpendicular-foot guard, or (b)) whose readings all lie within the parameter
   tolerance (1e-3 R) of the stretch's exact mean (a reading that would take the mean more
   than the tolerance from any member starts a new stretch). An exact stretch splits off as
   its own EXACT sub-piece, reading that mean, only when its exact samples span at least
   `ExactStretchMinLengthOverR` = 1 x R (dimensionless); a shorter one merges into the
   adjacent non-exact stretch(es). A sub-piece that is ONE exact stretch stays exact whatever
   its length, as before (the sub-R exact stack ends of the chips are not re-keyed by this
   clause: the threshold governs splitting only). The judged samples WITHOUT an exact reading
   (away from any bend within `PairBendProximityOverR` R) form NON-EXACT stretches by their
   status — they are no longer compared with the exact mean: a judged sample lacks an exact
   reading only where its partner is not parallel within the cosine tolerance or its foot is
   not perpendicular, so a coincidental agreement of its distance with the mean must not make
   it exact. The near-bend samples (not judged) never form, split or decide a stretch on
   their own: they join an adjacent stretch within `PairBendProximityOverR` R along the run of
   that stretch's last exact / judged sample — measured between the samples' CELLS, i.e. R
   plus one sample spacing (at most R / 2): the exact samples end just before (joint - R),
   so the near-bend sample at the joint lies R + up to one spacing from the last exact
   sample — (between an exact and a non-exact stretch the first R goes to the exact one; a
   short gap of them between two exact stretches whose readings agree is transparent, one
   exact stretch across the partner's noise joint); contiguous near-bend samples farther than
   R from every exact / judged sample form a NON-EXACT stretch of their own (the bend
   exemption is a 1 R zone, not a licence for an unbounded unjudged region). Consecutive
   non-exact stretches are one non-exact stretch reading the length-weighted mean of the
   chord readings as before. **Partner-run cut.** A non-exact stretch is additionally cut
   where its samples' foot crosses a run boundary of the partner chain, so that each
   non-exact sub-piece faces ONE partner run and reads that run's local chord mean — the
   value the partner's own piece reads; exact stretches are not cut (their exact group holds
   every agreeing piece). Reason: the two sides of a pair were discretised differently (a
   taper side per run, a straight side as one stretch over a hundred partner runs), so one
   long non-exact stretch over a changing partner (mean 2.58 um over a 2.0 -> 3.7 um taper)
   could not pair with the partner's per-run pieces and both sides read as IsolatedEdges
   over the taper start. This exposes a consequence of the existing 184 (3) link rule, now
   written down: chord pieces within the 5 % pair tolerance of an exact group join it, so in
   the reproducer the first ~73 um of the taper (2.0 -> 2.1 um) key WITH the exact 2.0 um
   strip, on both sides alike (the established pair tolerance, not a new error).
   **Invariant of the link grouping: mutually facing pieces belong to the same group.** Two
   non-exact pieces P and Q on the two sides face each other mutually when every foot of P
   lies on Q's run, every foot of Q on P's run, and each is the only non-exact piece of its
   run facing the other's run. The two sides sample one geometry with different sample sets
   and agree to ~1e-3 only, so at a group threshold (the 5 % boundary of a taper, w = 2.1 um
   for the 2.0 um exact group) one side's piece joined the exact group and its facing piece
   the chord group, leaving both unpaired: a strip edge read as an IsolatedEdge sliver of one
   chord — a misidentification. Where a mutual pair straddles a threshold, the piece sampled
   over exactly its own run (its interval is the whole run) decides and the other joins its
   group; when both or neither span their run, the piece on the lower run index decides.
   Each piece is in at most one mutual pair and every move is decided on the original
   grouping: deterministic and order-independent. The 5 % grouping itself (the post-pilot
   taper-reading question D6) is unchanged; the stage log counts the moved pieces with their
   length and sites.
   Consequences: a straight
   run facing a partner that is parallel over part of its length keys that part exact and the
   rest non-exact (the reviewer's reproducer, a 450 um run facing 150 um of straight 2.0 um
   strip and a quadratic one-sided taper to 3.7 um: before, one judged sample at the far end
   made the whole run non-exact at its 450 um chord mean, the 2.0 um exact group had no
   partner piece and the 0.53 R strip read as two IsolatedEdges; with every taper point
   within R of a noise joint (chords < 2R) and no judged sample at all, the twin mode: the
   whole run read exact at the lead's value); an extremely slow taper whose chords are
   parallel within the cosine tolerance and pass the foot guard keys a staircase of exact
   stretches, one per 2e-3 R of separation change (the parameter-tolerance semantics; the
   stage log counts the sub-pieces cut into stretches and the exact stretches per chip); an
   exact stretch absorbs up to R of a slowly changing partner (a mis-keying of at most the
   taper's slope x R: 0.006 R in the reproducer); the translational mortar-strip unit test's
   pure-shear "bend" shape (both strip sides sheared by 8 deg at mid-length) reads TWO exact
   strips, 0.25 and the perpendicular width 0.25 cos(8 deg) 1 % apart (two exact groups under
   decision 85(1)'s 1e-3 R, since 2026-10-01; formerly one chord-chained strip at a mean), and
   since the partner-run cut the clamped-foot region of the sheared side facing the straight
   partner across the kink is its own portion of the 0.25 strip (both its Gauss patches on a
   tilted frame; formerly one). Recorded census (asked for with the rule):
   the stage log lists the sub-R non-exact stretches that sit between two exact stretches of
   the same mean (fragmentation), with their sites. Considered alternative, not taken because
   it changes the judging globally: judging the near-bend samples too at the tolerance
   widened by their discretisation ambiguity d (1 - cos(turn)).
   The pieces of a chain pair are grouped into links so that a
   link holds ONE separation: the exact pieces within the parameter tolerance of each other
   form an exact group (value = their length-weighted mean), a chord piece joins the exact
   group within the pair tolerance (`PairSeparationToleranceRelative` = 0.05) of its own
   separation (the nearest when several: a chord reading of that design separation), and
   the remaining chord pieces form chord groups by consecutive steps within the pair
   tolerance (value = the length-weighted mean chord reading, `ExactParameters` false). A
   piece is never keyed by the exact readings of pieces beyond its own tolerance: formerly
   the pieces were chained by the 5 % step alone and "a class with any exact piece takes
   the mean of its exact pieces", so a slow taper linked an exact 2 um lead to a 3.7 um
   strip through its intermediate readings and the whole link took the lead's 2.0 um. A
   slow taper is still a pair described by its mean separation (the recorded limitation
   above); the taper READING rule (where the pair ends, the taper joints) is a USER
   question after the pilots, not settled here. Superseded text kept for the record: "A
   class with any exact piece takes the length-weighted mean of its exact pieces; otherwise
   it takes the length-weighted mean chord reading of the bent-pair rule (a slow taper, a
   spline route with chords beyond 2R) and the feature carries `ExactParameters` false (an
   annotation, not hashed)." Sample-count weighting is gone. The stack offsets are the sums
   of the
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
   instance (a function of the set alone) — with the group's `Instances` (the number of
   FEATURE instances in the group; `DistinctSignatures` the number of distinct signatures
   among them; the version-1 `Count` stays the mesh-segment / feature count of the record
   as before — one semantics in the version-1 record and in `signature_library.py`, review
   fix-3 m5) and `ParameterSpread` (max normalised deviation of a member from the representative; a
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
   and the fillet suite meshed at >= 3 mesh sizes plus a seeded node perturbation (decision
   88(3): the perimeter nodes moved radially by at most the signature parameter tolerance
   1e-3 R and the interior nodes placed by Gmsh's Delaunay algorithm, as the re-mesh gate's
   Delaunay variant; the former 1 %-of-chord vertex noise, 1.0-1.6e-2 R on the long chords
   of the q = 5-9 bends, was a different geometry under the exact-parameter contract) give the same
   topology keys and parameters within the tolerance; DS-SCT-001's three identical 2 / 2 / 2
   um routes give ONE key (0 / 1 / 2 / 3 R exactly, chirality 0; before: three keys 0.04-0.6 %
   apart).
   **Cluster extension (decision 85(2), 2026-09-26; review M4; the "across" wording ratified
   by decision 88(2)).** Every *single-edge portion*
   — the remainder of a run outside every cluster portion, vertex window, pair / stack claim
   and CrossLayer zone, i.e. what would become an `IsolatedEdge` / `CurvedEdge` — whose 3D
   distance to a cluster's claimed perimeter on another chain (or on its own chain beyond the
   pi R self-pair neighbourhood) is strictly below 2R, faced ACROSS (the perpendicular
   projection of the point onto the claimed piece falls inside the piece; at an interior
   joint of a claimed chain interval the piece's domain is extended by 2R tan(turn), the
   width of the wedge between consecutive perpendicular domains, with the turn capped at
   `ClusterExtensionWedgeCapDegrees` = 30 deg (formerly the corner threshold itself; kept at
   30 deg when the corner threshold became the 1 deg noise threshold so that the extension is
   unchanged — a piece end at a sub-noise joint or at an arc joint turns by less anyway, and
   a corner site's window end reaches at most this wedge; the cap bounds the wedge at
   2R tan(30 deg) = 1.15 R); never past the end of the claimed interval — a
   diagonal reach would creep along a sub-2R strip, whereas the partner beyond the claim
   end pairs with the free continuation; a claimed piece ending where its CHAIN ends — a
   strip end, the tangent point of a rounded corner's arc chain — faces the metal around
   that end within the full 2R ball, since nothing continues beyond it), outside the through-vertex
   zones of the shared vertices that are NOT members of that cluster, joins the cluster over
   that sub-interval. The same for a vertex feature outside every cluster: a single-edge
   portion faced across within 2R of its window (outside its own through-vertex zone — the
   arms of a corner meet through it) makes that vertex a cluster together with the portion.
   A portion within 2R of several owners merges them (union-find over clusters and joined
   vertex features; order independent; the merged cluster is numbered by its smallest
   member). Pairs and stacks are joint descriptions and are never absorbed — with ONE
   exception (**supervisor decision 224, 2026-10-02**; `Conventions.ClusterExtensionRule`,
   `Diagnostics.ClusterExtension.TranslationalPiecesAbsorbed`): a pair / stack stretch that
   exists ONLY because a cluster's claim boundary cut it is absorbed by that cluster. A
   stretch is a maximal contiguous run of one cross-section's claims (one signature key; the
   stack assembly can emit one feature per mesh segment) along a chain; it must lie ENTIRELY
   within the cluster ball radius `ClusterBallOverR` x R of the cluster's claimed pieces (the
   3D distance sublevel sets of the cluster's own claims, arcs included; the claims AT THAT
   PASS — the ball grows with the absorptions, so a k-member stack end sheds up to k - 1
   recomposition pieces over successive passes, e.g. stack-k4-1-1p5-3: 0.959 + 0.268 +
   0.268 um per end over three passes) — shorter than the
   interaction distance from the cluster's own claims, so the cluster's coupon describes it
   while the stack's translational patches would correct the same surface a second time —
   and be one of two classes, counted separately (`TwoSided` / `StackEndRecomposition`):
   TWO-SIDED, bounded at both ends (within the decision quantum; the chain's closed wrap
   counts) by claimed intervals of the SAME cluster = a piece shorter than 2R between two of
   its claims (the S1p stage-1 window's 41-edge loop end, key fcdbab7d58ad: its three
   vertical leads, the loop wire's edges and the ground edge at offsets 0 / 2 / 4 um, were
   each split into two 5.131 um cluster portions by a 1.738 um piece of the 3-edge stack
   69ca648cc16f, the six free ends being claim cuts inside the cluster); STACK-END
   RECOMPOSITION, adjacent to the cluster at one end and continuing, at the other, a pair /
   stack claim of a cross-section with strictly more members containing its own = the
   members a cluster takes later continue past the member it takes first (the members are
   compared as sets of CHAINS, `feature_members`; the device path splits physical chains at
   corners, so a pair's two edges are two chains and the member test is exact — two edges of
   one pair sharing a chain would not be told apart) (the S1p 3-edge
   stack end, 3 x 1.0 um at y = -136.3 between the 4-edge stack and the loop end's claims on
   the leads, inside the loop end's hull). Absorbed in the same pass as the single-edge
   absorptions (the next pass recomposes the stacks around the enlarged claims; the stretch
   is cut into its runs' intervals). NEVER between two DIFFERENT clusters (the placement's
   spatial-vs-spatial overlap check covers overlapping coupon boxes), never a stretch whose
   far end is free, a bend or a smaller cross-section: strips and gaps alongside a cluster,
   the halves of a bent strip and the gap middle between two pad-end clusters are genuine
   translational features the 2D models describe. On S1p the rule removes 69ca648cc16f
   entirely (8.21 um -> 0: the two-sided 1.738 um piece and the one-sided 1.0 um stack end;
   the other leads of each follow as single-edge remainders) and re-keys the loop end to
   5ed91f8890c0 (38 edges, 197.8757 um). Two broader forms were tried the same day and
   rejected: the two-sided class alone (left the stack-end piece inside the loop end's
   coupon volume, which the ownership record below lists) and the "ball only" form (any
   stretch adjacent to a cluster and within its ball: it ate 10 identity cells — DS-SCT-001
   +43.6 um of cluster, the arc bar of syn-arc-r5-step20, the gap middles of the
   syn-gap-different cells — and the bent-strip halves of the translational mortar strip
   test). A translational stretch the rule leaves inside a spatial coupon's volume is
   RECORDED by the placement's ownership record (decision 236, never an abort): a
   translational STRETCH (the identification's portion unit, carried as
   `IdentifiedPortion::stretch` and the patch provenance `Stretch` / geometry cache version
   6) whose every longitudinal cell lies strictly inside one spatial support's box is
   listed under `Diagnostics.TranslationalStretchesInsideSpatialSupport` (feature, stretch,
   model, length, box, class) with a warning (`FindTranslationalStretchInsideSpatialSupport`;
   judged per stretch, never per cell). The spatial support's box is the bounding box of the
   model's basis points placed by the patch frame = the cluster's claims + 3R past every
   claim-cut end along its edge (the coupon generator's `coupon_bounds`: every claim-cut end
   continued STRAIGHT by 2R, then R of padding; transversely R + R of padding = 2R; a coupon
   exists only when the plan span of its box is <= 16 R — not "about R"), so the first cells
   of every stack portion adjacent to a cluster lie inside its box legitimately.
   The coupon WAS calibrated on its claims AND their straight continuations to the box face
   (the decision-236 geometry; SUPERSEDED for the coupon metal by the spatial-support
   contract v3 below, under which the coupon metal is the device plan clipped to the box and
   a continuation follows the device chain through its vertices — the ownership record's
   Continuation class is unchanged: on the device the straight continuation and the chain
   coincide wherever a stretch continues a claim straight), so a recorded stretch is classed
   Continuation when it continues a claimed portion of that
   cluster through its claim cut BY THE SIDE'S OWN EDGE (decision 252): one of its cells lies
   on the claim's mesh segment (the claim boundary cut that segment), or it runs parallel to
   the claim (within the signature angle tolerance), one of its cell ends ON ITS OWN EDGE —
   the cell end shifted by the patch provenance `EdgeOffset` along AxisU: 0 for a single
   edge, -/+ half the separation for the two sides of a pair (whose cells sit on the
   midline), the side's offset from the first side for a stack (whose cells sit on the first
   side); carried by geometry cache version 6 with the mesh `Segment` — abuts a claim end
   along the chain within the signature parameter tolerance 1e-3 R and within the same
   tolerance transversely (a claim boundary snapped onto a mesh vertex of the same edge) AND
   it extends beyond that claim end (every cell end on the outward side of the abutting end,
   away from the claim's other end, within the tolerance; a parallel stretch lying alongside
   the claim over the claim's own range, one end aligned with a claim end, is Foreign — it
   does not go through the claim cut). The side of a pair or stack whose own edge the cluster
   does not claim is Foreign whatever its cells' proximity to the claim end (before decision
   252 the transverse reach was R on the midline cells, so the unclaimed edge's half of a
   pair correction was removed wherever its cells abutted a claimed neighbour's cut) — the
   double count the continuation ownership at placement
   removes (the accepted transmon library `transmon-r1p9-folded-refined-corners`: 6 records,
   all Continuation, 13.28 um, 0 Foreign, GeometryDigest 9ada660bf6e4 unchanged — the two
   sides of each 2-um strip (features 0 / 1, 4 x 0.869 um from the 3-edge clusters' claim
   cut at y = 540.131 to the feature split at y = 541.0, x = +-43) and the two 4.9-um
   isolated edges (features 70 / 73, y = 12) continuing the 10-edge cluster's claims from
   x = +-5.2, ending 0.8 um inside its box face at +-10.9 = 5.2 + 3R) — and Foreign
   otherwise (metal absent from the coupon's twins: a second-order model mismatch, not a
   double count). The record runs in the operator constructor (metadata
   `SurfaceResponse.Diagnostics`, with the warning) and in the preflight on the spatial
   models whose basis points the library provides (manifest `Identification.Diagnostics`
   and the same warning; a signature placeholder of a Missing key has none and is listed
   under `SpatialSupportsWithoutBasisPoints`). HOW THE COUPON'S MARGIN AND THE
   TRANSLATIONAL PATCHES MEET (continuation ownership at placement, decision 236 (2) / 244,
   2026-10-02, `ApplyContinuationOwnership`): the cells of every translational stretch that
   continues a cluster's claim (the Continuation criterion above, judged on the whole
   stretch whether or not it lies inside the box) are OWNED by that coupon inside its box —
   the coupon's twins already carry the straight continuation of the claim to the box
   face — and CLIPPED exactly at the face: the kept part of a cell is the part outside
   EVERY box whose claims its stretch continues (a cell on the continuations of two coupons
   is removed once; no first-wins), the patch weight and its provenance quadrature weight
   scale by kept / cell (a cell wholly inside keeps weight 0 and the operator skips it;
   the portion [S0, S1) stays, so a portion's quadrature weights sum to 1 - owned /
   portion length, which the A7 audit reconciles from the record), the origin moves to the
   kept interval's midpoint with the cell symmetric about it (the clipped cell equals a cell
   built directly on the kept interval), the Maxwell anchors move with it. Foreign cells,
   the cells of a continuing stretch outside the box and the stack-end cells of a stretch
   that continues no claim are untouched; a curved cell never continues a claim (parallel
   within the signature angle tolerance), so an arc continuing an arc claim keeps its
   patches — a residual double count the record lengths show; the cell of one SIDE is the
   quantum: a pair's or stack's side is owned only where its own edge continues the claim
   (the segment branch or the own-edge abutment above), the cells of an unclaimed side stay
   (unit case: a pair with one edge claimed — that side owned, the other untouched; the
   transmon's owned strips have both edges claimed and every owned stretch is classed by the
   segment branch, so its record is unchanged). `Diagnostics.ContinuationOwnership` lists every owned
   cell (patch, feature, stretch, portion, owned length, owners with the attributed length —
   a shared cell split at the midpoint between the two continued claim ends) and the owned
   length per spatial support; `TranslationalStretchesInsideSpatialSupport` carries the
   owned length per record and in total. The accepted transmon library: 128 cells / 112.40 um
   owned (110 wholly, 18 clipped; 44.0 um of isolated-edge cells on the 3-edge and 10-edge
   coupons' continuations + 68.4 um of 2-um strip midline cells whose both edges the 3-edge
   and 4-edge clusters claim), 9.5e-4 of the translational patch weight. COUPON VS COUPON
   (decision 244, `FindSpatialSupportMarginOverlaps`): two cluster boxes overlapping in
   their interiors no longer abort when the overlap is MARGINS ONLY — no claim SEGMENT of
   either enters the bounding box of the other's claims by a positive length beyond the
   tolerance (tested on the segment, so a claim crossing the hull with both ends outside
   counts; decision 252); a claim inside the other's claims still aborts. Each coupon continues every claim CUT end (an end no other
   claim of the same cluster shares within 1e-3 R) straight to its own box face; the length
   of those continuations lying on the other coupon's claims (margin over claims) or on the
   other coupon's continuations (margin over margin, the bridging piece between two claim
   cuts) is corrected by both coupons — a double count the placement cannot remove (the
   dense coupon operator is not clippable), recorded under
   `Diagnostics.SpatialSupportMarginOverlaps` per pair with a warning (the S1p 19-edge
   pair 5f0ccfcb3b8b / bd43654a77c6: boxes x [568.55, 592.75] and [583.05, 607.25], margin
   over claims 2.80 + 3.80 um, margin over margin 2.90 um = 9.50 um, 7.8 % of the two
   coupons' 122.2 um of claims, 0.33 % of the device's corrected edge length). FOLLOW-UP
   (option D of the design note): shrinking the continuation from 2R to R reduces every
   margin double count, with library rebuilds. Unit tests
   `SurfaceResponseIdentificationStackPieceInsideCluster`,
   `SurfaceResponseOperatorTranslationalStretchOwnership`,
   `SurfaceResponseOperatorContinuationOwnership`,
   `SurfaceResponseOperatorSpatialSupportMarginOverlaps`. The extension
   iterates to closure: an enlarged claim moves the stack ends (the stack-end rule
   recomposes the members on the new claims) and the recomposed stacks leave new single-edge
   portions to test; the loop stops when a pass absorbs at most the signature parameter
   tolerance `ClusterExtensionClosureOverR` = 1e-3 R in total — below the resolution of
   every parameter and gate. That last pass IS applied and the pairs / stacks are not
   recomposed again (stated, review fix-3 m8): its pieces come from the unclaimed remainder,
   so no pair / stack claim overlaps them and the partition gates hold; the stack cut images
   on the other member chains then differ from the final claims by less than the tolerance
   (a sub-tolerance inconsistency of an exact-match cluster signature's last piece). The
   alternative — testing the closure before applying — was tried and rejected (2026-09-26):
   it left the sub-tolerance unclaimed remainder as isolated slivers facing the cluster metal
   (DS-SCT-001: a 0.24 nm `IsolatedEdge` failing the A8 facing gate), and recomposing after
   the last application re-creates such slivers at the moved cut images (the 1 nm slivers
   of DS-SCT-002's passes 3-15). (DS-SCT-001: 2 passes, 127 portions, 44.7 um; transmon
   12.8 um; two-transmon chain 25.6 um — the former `ClusterNeighbour` lengths plus the
   recomposed pair-side pieces; DS-SCT-002 at R 2.1: pass 1 absorbed everything, passes
   3-15 had re-cut 1 nm slivers for 8 s each). The pair / stack claimed length satisfying
   the SAME across rule (within 2R of a cluster's claimed perimeter or a free vertex
   feature's window, faced across, outside the through-vertex zones of non-member vertices;
   `ExtendClusters` evaluated on the pair / stack claims without absorbing) is the
   *stack-end third body*, reported as `Diagnostics.StackEndThirdBodyLength`
   (`Diagnostics.StackEndThirdBodyRule`; DS-SCT-001: 155 um) and sampled by the facing gate
   as the pair-side class `StackEndThirdBody` (0.5 R samples, across = the facing direction
   within 60 deg of the sample normal, the through-vertex zones excluded first: the two
   readings are one definition, review fix-3 m5). **Why the two readings differ on DS-SCT-001
   (reconciled interval by interval, review fix-4; `PALACE_IDENTIFICATION_DEBUG_EXTENSION`
   dumps the measured intervals, `facing_check --excluded-samples` every excluded sample):**
   every gate sample read as `StackEndThirdBody` or `ThroughVertex` lies INSIDE an analytic
   interval (R 2.0: 34.9 + 92.9 = 127.8 um of the 155.1 um; R 1.9: 21.0 + 51.3 = 72.3 of
   100.8 um) — the gate's `ThroughVertex` class is a finer split of the same length, since
   the analytic rule excludes only the through-vertex zones of vertices shared by the two
   chains that are not the owner's members; nothing the gate reads as third body lies
   outside the analytic intervals (`StackEndRecomposition`, pair facing pair, is by
   definition not third body). The remainder (30.6 um at R 2.0, 29.4 um at R 1.9, minus
   0.9-3.3 um of 0.5 um sample discretisation at the interval ends) is stack claim length
   whose facing cluster or corner-window claim lies ON THE STACK'S OWN MEMBER CHAINS within
   its lateral reach — the flux-line ends, where a member chain continues into the
   termination corner cluster — which the gate's own-side rule (portions + member chains
   within reach, review fix-3 M5) masks before any facing test (28.6 um of it is across in
   the gate's sense, 2.0 um only within the analytic perpendicular domain / wedge). One
   definition, two instruments: the analytic reading is the authoritative
   `StackEndThirdBodyLength`; the gate's number is a lower bound that equals it wherever no
   member chain carries cluster metal within the stack's reach (the chip: 45.9 = 45.9 um).
   Reported metric (decision 88(2), not a
   gate): `facing_check.ClusterProximityNotAcross` — the isolated / curved edge length with
   cluster or vertex-feature metal within 2R in ANY direction and no across hit, i.e. what
   the across rule leaves single-edge, split into the `ThroughVertex` / `SelfNeighbourhood`
   exemptions and the residual. Facing gates (audit A8,
   `facing_check.py`) after this decision: every facing segment within 2R is tested (the
   nearest alone masked a second edge behind an excluded one); the own sides of a pair /
   stack are the segments of its portions plus the member chains within the feature's reach
   (its lateral span x 1.05) of them, not the whole chains; the remaining exemptions are
   `AtExactly2R`, `ThroughVertex` (design item 1: the arms of one vertex feature meet through
   it — the identification's own through-vertex convention, kept), `SelfNeighbourhood`,
   and for pair / stack sides only `StackEndThirdBody` and `StackEndRecomposition` (facing
   another pair / stack that shares a member chain: the route's cross-section recomposed at a
   stack end or class boundary, the partner's foot across the cut; sub-chord lengths,
   DS-SCT-001 2.1 um), each with its length. The `StackEndRecomposition` exemption is
   BOUNDED (review fix-3 m4): a recomposition is a sub-chord event, so one site (one ordered
   pair of features) may exempt at most `STACK_END_RECOMPOSITION_SITE_CAP_OVER_R` = 1 R and
   the manifest at most 1 R per pair / stack feature involved in total; beyond either cap
   the length is not exempt (two overlapping two-edge pairs across a sub-2R trace — the
   decision-78 defect — would run along the whole route and fail the gate); the site count,
   total and largest site are reported in R. The per-site cap is the OPERATIVE bound (review
   fix-4 m-D): with every site at most 1 R the per-feature total can only bind when one
   feature carries more than one site per involved feature, so the total is reported against
   its cap but the site cap is what refuses a route-long recomposition (DS-SCT-001: 4 / 5
   sites, largest 0.38 / 0.43 R, total 1.05 / 1.30 R at R 2.0 / 1.9). One further recorded
   exemption of the isolated gate (2026-09-26; NARROWED and ratified by decision 89):
   `SubToleranceFeature` — an isolated / curved PORTION shorter than the signature parameter
   tolerance 1e-3 R that is an unclaimed remainder, i.e. bounded on both sides along its run by
   OTHER features' claims (DS-SCT-001 at R 1.9: one 0.18 nm `IsolatedEdge` piece between a
   `SameConductorGap` claim and a `SpatialEdgeCluster` claim on one segment, facing the
   neighbouring route at 2.0 um), or a WHOLE isolated / curved feature shorter than 1e-3 R; it
   has no resolvable parameter, the exemption is bounded to 1e-3 R per such portion and their
   count is reported. A sub-tolerance MESH SEGMENT inside a long isolated edge — its run
   neighbours claimed by the same feature — is an ordinary portion of that feature and is NOT
   exempt (the former per-segment predicate exempted 896 whole 1-2 nm segments of one 2,898 um
   `IsolatedEdge` on DS-SCT-001 at R 2.0, none of them facing; under the narrowed rule 0 exempt
   portions at R 2.0 and 1 at R 1.9, both gates PASS). Since the sliver rule of item (b) 5
   (decision 222) such sub-tolerance remainders join their adjacent portion in the
   identification itself, so the exemption is expected to count 0 on every manifest of this
   version (it stays in the gate as the record of the class and for older manifests). Diagnostics also record the stack assembly's
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
   pair's own separation x 1.05 (the mean, or the piece's local separation when larger) —
   both recorded here as the two non-2R constants that remain, neither decides an
   interaction. Recorded taper residual: the chord reading is the window MAXIMUM over R (a
   full chord on a curved chain), so along a slow taper a pair ends where that maximum
   reaches 2R — up to slope x window before the strict local crossing (DS-CTX-003 launchers:
   0.03-0.13 um = 4-16 um along the edge; audit A8 81.8 um at 5 sites after the fix, 529.5 um
   at 10 sites before); a rule question for the USER list, not a threshold to add here.
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
portions; z = the cluster's OWN signed process normal (substrate -> vacuum: the reference
direction oriented by the claimed-length weighted mean of its runs' signed normals,
`FeatureProcessNormal`, fail closed naming the cluster when its runs disagree in sign — see the
frame rule below; never the sign-canonical device reference normal n_ref itself). The in-plane
x axis is chosen among the finite candidate set
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

**Frame rule: every spatial / vertex feature is framed with its OWN signed process normal
(decision 266, 2026-10-02; flip-chip devices).** The device-global `ReferenceProcessNormal`
n_ref is sign-canonical (`SignCanonical`: the first nonzero component positive), so on a
flip-chip device the upright bottom chip (normal +z) and the flipped top chip (normal -z, its
substrate ABOVE its metal) share n_ref = +z; it is a planarity / parallelism reference only
(|dot| tests for `NonPlanar`, the plane numbering, the in-plane direction classes, the
translational lateral order) and is never a frame normal. A coupon is built upright (support
w in [-1.95, +2.0] R: substrate at w < 0, vacuum at w > 0), so a frame whose w is n_ref on a
flipped plane places the coupon's substrate in the vacuum gap and its vacuum in the top chip's
substrate — the media swapped for every spatial feature of the flipped plane (the S1p defect:
its two 19-edge clusters and its L2 corners, AxisW = +z at z = 4.8, while the translational /
curved / stack patches, built from the segment frames with the segment's own signed normal
(AxisV = -z), were right). The rule: `SpatialEdgeCluster` z = the signed normal of its runs
(claimed-length weighted mean, `EmitClusters`), `ConvexCorner` / `ConcaveCorner` / `Endpoint` /
`Junction` z = the signed normal of the site's incident runs (length weighted,
`SiteProcessNormal`; the junction's arm angles are measured about the same normal), the
isolated / curved representative axes (informative; their patches come from the segment
frames) (t, n x t, n) with the run's own signed n; a feature whose runs disagree in sign has no
substrate -> vacuum side and the identification fails closed naming the feature
(`FeatureProcessNormal`). The direction returned is +-n_ref oriented by the feature's own runs
(`SignedReferenceNormal`: every framed run is parallel to n_ref within the NonPlanar tolerance
1e-8, so this IS the feature's normal, and single-plane frames stay bit-identical — the transmon
patch CSV is byte-identical). The frames are right-handed (x, n x x, n) where the feature is its
own mirror image or its canonical chirality is +1; the mirror frame (x, -(n x x), n) of a
chirality -1 cluster (and the gap-side choice of an endpoint, the arm-order choice of a
junction) is a PLAN-VIEW mirror, which keeps w on the vacuum side and maps the model onto its
mirror-image feature exactly (the response is reflection invariant). Key changes of flipped-chip
spatial features: NONE of the hashes — a cluster signature is serialised in both handedness
values and its hash is handedness-invariant (a plan-view mirror image has the same hash with the
opposite `Chirality`; corner and junction signatures are mirror-invariant by construction), so
the fix changes, for every spatial / vertex feature of a flipped plane, the frame's w (-z instead
of +z), the cluster's `Chirality` (its sign flips: the old reading was the mirror image's) and,
for corners, the in-plane axes (x = the other arm, the frame rotated by 180 deg about the
corner's bisector); the in-plane placement map of clusters, junctions and endpoints is
unchanged (y = chirality x (n x x) is the same vector for both readings), so a Signature-keyed
coupon built for the old key fits the new frame without a rebuild and only its w side turns
over. Single-plane devices (the transmon) and the upright plane of a flip-chip device are
unchanged (n = n_ref there): the transmon preflight digest 9ada660bf6e4, records and patch CSV
are byte-identical. `GeometryDigest` hashes signature keys only and does not change. Unit test
`SurfaceResponseIdentificationFlippedPlaneFrames`: two planes carrying congruent copies of one
asymmetric cluster + corners (the flipped plane = the upright one rotated by 180 deg about an
in-plane axis): on the flipped plane w = -z, the coupon's w < 0 half-space lands in that plane's
substrate, (hash, chirality) and the rigid image of the upright frame are reproduced, the old
reading (the same geometry framed with +z) is the mirror (same hash, opposite chirality, w = +z),
the upright plane reads as alone, and a mixed-sign cluster aborts naming the cluster.

**Spatial-support contract v3: the coupon metal is the DEVICE PLAN clipped to the support box
(USER decision 281, supervisor decision 282, 2026-10-03; design note
`device-plan-coupons-20261003/DESIGN.md` B1-B3 / T1-T3 with the review rulings MAJOR-1 (i) and
MAJOR-2).** Rationale (blocks D2 / D3-C, decisions 277-278): under the decision-236 contract a
spatial coupon's metal was the claims + every claim-cut end continued STRAIGHT 3R to the box
face, whatever the device does there. On the S1p 19-edge junction-lead coupons the island
protrusion's edge was continued straight past the device's 90-deg corner 1.16 R beyond the cut:
26.7 / 19.9 um^2 of FICTITIOUS island metal at the island potential where the device has gap
at 10-19 V; the window trace imposed the device potential on the face knot rows 50 nm under
that metal, an MS-selective corner field that the decisive replay (D3-C, PBS 55513 / 55514)
measured as 99.8 % of the +5.79e-16 J MS surplus (x1.92 of the reference band, 39 points of the
window MS); the neighbour's leads (foreign metal) were absent from both twins; C2p shows the
opposite form (continuations running INSIDE device metal: a fictitious EDGE the conductor
consistency gate cannot see). The contract therefore changes at the identification, where the
key is formed: for every `SpatialEdgeCluster` the identification computes, in the cluster's
canonical frame and in units of R, (1) the claims-derived SUPPORT BOX (rule B2 = the coupon
generator's `coupon_bounds` / `edge_rows` / `extended_interval`, ported as
`SupportBoxFromSignature` on the serialised `{Portions, Vertices}` of the frame — the same
numbers the Python builder reads; `signature_library.cluster_support_box` is its Python
mirror, bit-identical on the transmon / S1p / C2p clusters): every claimed portion is a row
along its own tangent, a claim-cut end (no vertex and no other portion end within
`SupportEndCoincidenceOverR` = 1e-5) is lengthened to at least R from the row midpoint, every
row end at or beyond R continues by `SupportContinuationOverR` = 2R, the rows are widened by R
on both sides and the bounding box padded by `SupportPaddingOverR` = R (a claim cut is 3R from
its face, the claims 2R from the lateral faces); an arc portion is chorded at
`ClusterArcChordStepDegrees` / `ClusterArcChordMaxLengthOverR` as the builder chords it; (2)
the CONTEXT: every perimeter run of the cluster's plane (straight, or on its fitted arc —
the identification's `Runs` / `Arcs`, never the raw mesh polyline) clipped to the box, the
cluster's own claims removed, in the portion encoding `{P | Arc, Gap | GapRadial, Conductor,
Interfaces, Law}` plus `Chain`: true on the pieces connected to the claims inside the box
through run ends and device vertices (the OWN continuation chains, rule B3: the within-R
accounting and the continuation ownership follow them) and false on FOREIGN edges (present in
both coupon twins so that the fields are consistent; excluded from the within-R accounting;
their own patches untouched). A piece lying on ANOTHER spatial cluster's claimed interval is
foreign by definition (decision 285 (1), the R1a review's MAJOR-1: the other coupon owns it):
the chain stops where such a claim begins and never enters it, the context pieces are split
at the other clusters' claim cuts so that the boundary is representable, and the manifest
census lists those pieces under `Context.ClaimedByOtherFeature` (count, length, the owning
cluster and feature — not hashed). Chains along device edges claimed by TRANSLATIONAL
features continue as before (continuation ownership, decisions 236 / 243-252). Conductor
labels by first appearance over the sorted `Portions` THEN the sorted `Context` (a foreign
conductor touching no claim takes the next label); (3)
the FACE RULES, dimensionless: T1 a piece end within `SupportFaceSnapOverR` = 1e-3 R of a face
is ON the face (a crossing) and a context piece shorter than that is dropped (the sliver
quantum of decision 222); T2 every device edge inside the box keeps
`SupportFaceClearanceOverR` = 0.25 R (the trace basis' finest knot scale) from every face it
does not cross, every piece end not on a face (device vertex, claim end) and every crossing
keep that clearance from every other face (a crossing near a box corner), every crossing meets
its face at sin(theta) >= 0.25 (theta >= 14.5 deg), and two crossings of one face closer than
0.25 R may bound METAL (a narrow lead: allowed, `FaceRules.MinCrossSectionOverR` /
`NarrowCrossSections` recorded) but not gap (a channel the trace basis cannot resolve). T2 is
TWO-SIDED (decision 285 (2), the R1a review's MAJOR-2 — the face trace cannot represent a 1/r
edge field within 0.25 R of the face on either side): the plan is also clipped to the box
dilated by the clearance, and in that shell a run end lying on neither boundary (a device
vertex or joint within 0.25 R OUTSIDE the face) or a piece entering and leaving the shell
without crossing the face (an edge running within 0.25 R outside the face) fails that face
(`FaceRules.ExteriorVertices` / `ExteriorEdges`, summed over every T2 pass; their distances
enter `MinClearanceOverR` and the threshold band); the exterior tail of a face crossing, from
the face to the shell boundary, is the crossing itself. An end that exists only because
ANOTHER cluster's claim cut split the run there is no device vertex (R1 final review MINOR-1,
decision 288 (3)): it is exempt from the end tests (still a crossing when it lies on a face),
and the shell clip is not split at cuts, so a key never depends on a neighbour's claim EXTENT
beyond the Chain / ClaimedByOtherFeature classification. QUANTISED COMPARISONS (decision 287
(b), the R1 final review's MAJOR-1 — the knife edge): the box faces sit on the 1e-6 R
signature grid (`SupportComparisonQuantumOverR`) and a growth step moves them by a grid
multiple, so a device edge lying ON a claims-box face (the S1p / C2p 19-edge pattern: the
neighbour lead's claimed edge on the face) reads EXACTLY 0.25 R from the moved face up to the
quantisation residual (half a quantum) and mesh / float noise; the former
`distance + 1e-9 < clearance` decided both production keys on that noise. Every length
threshold of the face rules is now read on the grid: a value within half a quantum of the
threshold is AT the threshold and takes the rule's own inclusive side — the vertex /
crossing clearance and the piece clearance pass when `distance >= 0.25 - 0.5e-6` R, the
crossing separation bounds a gap only below that, a piece end is on a face when
`distance < 1e-3 + 0.5e-6` R (on-face, touches, the shell boundary, the beyond-side test) and a
piece is a sliver only below `1e-3 - 0.5e-6` R. Left as they were, with the reason: the
crossing sine (dimensionless, on no length grid), the span cap (compares box coordinates that
are already quantised, 1e-12 relative) and the shell membership (the clip reads with the snap).
Consequence: keys are geometry-precise to the 1e-6 R quantum — a one-quantum move of a device
edge is a real geometric change that re-keys (and may legitimately change the growth: an edge
one quantum inside the clearance fails it) — and sub-quantum noise (a regenerated mesh of the
same device, float round-off, the quantisation residual) cannot move a key or a growth
sequence (unit scene 11: +-4e-7 R perturbations of the corner and of the frame-defining lead
give the same growth and the same key; +-1e-6 R re-keys). T3 every face failing T2 moves
outward by `SupportFaceGrowthStepOverR` = 0.25 R per step — all failing faces of a step
together (review MINOR-2: one order) — and T1 / T2 are re-evaluated on the grown box (new
edges enter), at most `SupportFaceGrowthMaxSteps` = 12 steps per face and never past the plan
span cap `SupportSpanCapOverR` = 16 R (`Growth.Steps` counts the applied steps,
`AttemptedSteps` also the refused one, `StepReasons` the first failure behind every applied
step, `Growth.Passes` EVERY T2 pass with its box, failing faces, every failure — R1 final
review MINOR-8: a step growing two faces lists both reasons — and its band / exterior
readings; `FaceRules.ThresholdBandHits` sums the band readings of every pass (decision 287
(a): a key whose growth was decided at a threshold is flagged whatever pass read it;
`FinalPassThresholdBandHits` the last pass alone), the identification log and the preflight
warn on every cluster with a non-zero count, `Diagnostics.SpatialSupport.ThresholdBandHits`
sums the clusters); a claims-derived
box already beyond the cap is keyed and recorded `ExceedsSpanCap` (the builder's cap and its
per-case override, decision 244 (i), decide as before) but may not grow (decision 285 (3): the
cap bites through growth only); a cluster no box satisfies is an `UnboxableFeature`: its
signature carries `"Unboxable": true` (a Missing placeholder no builder makes — never a
knife-edge coupon) and the record names the face and the reason (`Unboxable` text,
`UnboxableReason` = `GrowthStepsExhausted` | `SpanCapRefusedGrowth`). The two thresholds of
the box rule itself (the 1e-5 end coincidence, the 2R continuation of a connected row end at
R) are read within the knife-edge band too (`FaceRules.BoxRuleThresholdBandHits`). The
perimeter segments inside the box excluded before the identification (a port cut, an
undetermined process side, a non-manifold or non-planar edge — never a window truncation
cut) ARE plan pieces (rule B1, "nothing inside the box is omitted": a lead ending on a lumped
port keeps its end edge in the coupon's metal boundary): straight, foreign (`Chain: false`,
never entered by the chain, excluded from the within-R accounting as the device's own
perimeter excludes them), hashed with the context and listed under
`Context.ExcludedSegments` with their exclusion class. KEY
RULE (ruling MAJOR-1, option (i)): `Box` + `Context` enter the hashed signature ONLY when the
context is not what the decision-236 contract already draws, i.e. unless every context piece
is a straight chain piece abutting a claim-cut end, collinear with that claim and reaching a
face of the box (`LegacyEquivalent`) with no growth — then the device plan clipped to the box
IS the legacy coupon geometry and the claims-only key stands byte for byte (contract 2:
unchanged coupons keep their keys, no library migration; the transmon's 10-edge JJ and 4-edge
models and every corner). Otherwise (contract 3) the signature is the lexicographically
smallest `{Box, Context, Portions, Vertices}` over the candidate frames of the CLAIMS (the
portion tangents and perpendiculars, both signs and handedness values, as before) with the box
and the context recomputed in every frame (`CanonicalClusterSignatureWithSupport`; the box
frame IS the canonical frame, ruling MAJOR-2 — for axis-aligned clusters every candidate frame
gives the same geometric box, for an arc cluster the never-built loop ends change), so that
mirror images keep ONE key with opposite `Chirality` and rotated / translated / re-meshed
copies the same key (unit test `SurfaceResponseIdentificationSpatialSupportContract`);
`EdgeCount` stays the claimed count and the model `Edges` / the claims / the A10 placement
check / the continuation ownership at placement are unchanged. Vertex features (corners,
junctions, endpoints) are unchanged in signature and key: a corner coupon's box is [-R, R]^2
about the vertex and any other perimeter within 2R joins a cluster, so the device plan
clipped to a corner box is the two arms (56 / 56 production corners unchanged, R0 census); a
vertex feature lying inside a cluster's box on a chain is listed under
`Context.ChainVertices` with its distance from the nearest face for the placement's vertex
ownership (rule B4: a corner closer than R to a face has an arm partly outside the box, review
MINOR-1 (a), and a vertex inside two boxes appears in both records, MINOR-1 (b)), the others
under `ForeignVertices`; every entry names its owner (`Feature` = the vertex feature's id for
a free site, or the id of the cluster whose claims hold the site, with that `Cluster` index
— a site of another cluster is never a chain vertex), so that the placement never owns
another cluster's claimed vertex. The per-feature record plus the `Diagnostics.SpatialSupport`
summary stand in for a `KnifeEdgeCensus` band on T2 (review MINOR-7: a vertex shared by two
pieces is read once per piece end). The manifest records the whole evaluation per feature
(`Features[].SpatialSupport`: the claims box and the grown box, the growth steps per face, the
face-rule readings with their threshold-band hits at `KnifeEdgeBandRelative`, the context
census, the legacy contract's straight continuation and its FICTITIOUS part lying on no device
edge — 7.58 R / 4.95 R = 14.4 / 9.4 um on the S1p 19-edge coupons, the D2 census figures —
and the window truncation inside the box) and a summary under `Diagnostics.SpatialSupport`.
Census on the production devices with the first binary (R1a, 2026-10-03): transmon 1 of 3
cluster keys changes (the 3-edge `da0179c64591` -> contract 3: one foreign ground piece
4.38 R of edge with a foreign convex corner, no growth; 10-edge / 4-edge contract 2 and the
34 corners byte-identical, so the accepted library's models keep their keys except the
3-edge pair); S1p 3 / 3 change (the two 19-edge coupons: fictitious continuation gone, the
neighbour's leads foreign context, the boxes grown by 0.5 R / 0.25 R on the face the
neighbour lead's corner sits on; the loop end: contract 3 at 17.96 R > the cap, as unbuildable
as before), 10 / 10 corners unchanged; C2p 3 / 3 change (the 45-edge loop end UNBOXABLE: its
19 R box may not grow), 12 / 12 corners unchanged. Re-census with the decision-285 rules
(R1b): S1p 48e28baa8abb (grown 4 steps on x1) / 9f11f35e955e (1 step), the loop end
UNBOXABLE (SpanCapRefusedGrowth: a device vertex 0.105 R outside face x1); C2p 3d077d2793f3 /
0b6c2bd06426 (grown [0, 2, 0, 4]), the loop end unboxable; transmon 3-edge 35b2900ae07f and
10-edge unchanged, the 4-edge re-keyed 22c658057893 (decision 286 alias); every corner
unchanged (34 / 10 / 12).

The coupon BUILDER (R1b, `cluster_signature_geometry.cluster_coupon`) consumes `Box` +
`Context` (`signature_library.cluster_plan_view_edges(..., include_context=True)`): a
signature with a `Box` is built in its OWN frame (the canonical frame, M = identity:
`generate_spatial_response.frame_from_geometry` / `trace_basis.process_frame` / the mesher's
`process_frame` read the model's `SupportBox`), inside the signature's Box (the generator's and
the mesher's `coupon_bounds` take it: `Geometry.SupportBox`, the library model's `SupportBox`),
and its metal is the arrangement of the claims plus the context pieces cut at the faces (the
metal side = -Gap; `plan_view_faces` fails closed on crossings, disagreeing metal sides and
box-only faces): no straight extension, no `interior_bridges`, no fictitious metal. The
DESIGN's assertion `coupon_bounds == Box` is not a run-time check: the signature's `Box`
REPLACES the generator's box (`coupon_bounds(..., support_box)`), the C++ record being the
authority, and the two-language identity is the evidence (`tools/box_identity.log`: Python
`cluster_support_box` vs the C++ `ClaimsBox` 0.0 R on 11 / 11 production instances; R1 final
review MINOR-6). Open for CONTRACT-2 keys (R1a MINOR-6): a legacy-equivalent cluster with
mutually oblique portions is judged on the canonical-frame box and built on
`frame_from_geometry`'s (first gap direction) box; every production legacy cluster is
axis-aligned. The
generator's rows are the exact claims plus the context rows (`Context`, `Chain` columns of
mesh-signature.csv; the trace basis ignores them — the face knots of every conductor
cross-section follow from the mask vertices on the faces, B6, and the interior cap hats stay
within R of the claims; the mesher's owner lookup and attributes read them, so every
conductor of the plan, foreign ones included, gets its metal / SA / MS / MA surfaces); the
model carries `ContextEdges` (the run config's further conductors: terminal attributes,
conductor states) and `ForeignEdges` (the `Chain: false` pieces as 3D segments) which
`case_inputs.derive` writes into every Dielectric entry as `EdgeExcludeSegments` (rule B5: the
within-R accounting covers the coupon's own edges only; `EdgeExcludeSegmentTolerance` defaults
to 1e-3 R). A claims-only signature takes the legacy path unchanged — byte-identical generator
inputs (`test_cluster_signature_geometry.LegacyByteIdentityTest` against the a79b6af748 fixture).
Device-plan coupons are LARGER than their legacy twins (the plan inside the grown box): the
pre-build element estimate may exceed a suite's cap (the fixture transmon's JJ coupon at 6.4 M
> 6 M), a legitimate fail-closed outcome of the headroom gate; `test_device_coupons` therefore
asserts `PreflightPassed` only for the cases that pass the gate (R1 final review MINOR-5: the
former "every fixture case passes" is gone by design), and R2's S1p rebuild must re-estimate
its device-plan coupons against the 6 M cap.

**(F) qualification upgrade (decision 282 with the MAJOR-3 ruling, 285 (5);
`qualify/spatial_qualification.py`, `coupon_library.py spatial-qualify`).** `traces` writes
the dense held-out traces of a coupon source directory — T2: the C - 1 conductor states and
the potential of a unit line charge 5 R outside each in-plane face at the process plane with
every conductor GROUNDED (the DESIGN's "superposed with the states read at the references"
needs conductor potentials other than 0 / 1 V in one excitation, which PrescribedPotential's
one-volt `TerminalAttributes` cannot impose for two states at once: a trace with one nonzero
state is scaled to the 1 V terminal and its energies by s^2, a trace with two distinct nonzero
states is recorded `Unsupported`, never approximated — a `TerminalPotential` field of Palace
would lift it; this deviation from DESIGN section 6 is the (F) contract for R2 unless the
supervisor rules otherwise, R1 final review MINOR-4; the S1p 19-edge coupons have one
conductor state and are supported); T1: the device trace of a registration run's surface-response-traces.csv
(coefficients relative to the reference conductor, the states after the contour) — and the
fabricated / thin solve configs at the library order p4 and the control order p5 (the run
config's sources replaced, the response matrix off). `evaluate` reads the four runs' within-R
energies per class (SA, MS, MA raw — "MA sharp" when the run carries radial shells — and the
domain: `p_surf E_elec - E_out` at the largest R; domain-E.csv, surface-Q.csv,
surface-Q-edge.csv) and the p4 basis runs' matrices (domain-response-matrix.csv,
surface-response-matrix.csv, the whole-interface group), and judges, per class and trace, the
closure `|E_thin,p4 + t^T (Q_fab - Q_thin) t - E_fab,p5| / E_fab,p5 <= 0.02`, the DOMAIN
twin-consistency `|dE_dom,p5 - dE_dom,p4| / E_fab,dom,p5 <= 0.02` (decision 299, below) and the
matrix identity `|E_fab,p4 - t^T Q_fab t| / E_fab,p4 <= 1e-6`; the thin twin's SURFACE p-steps
are recorded as information; the conductor-consistency gate of decision 277 is the rebuild
acceptance, resolution-aware (decision 299 (2), below: the model's records on every supplied
registration-device solve all probed with MaxRatio <= 1e-4, Count 0, the solves at the
production orders p4 and p5 on the c0 mesh; a model without a record is untestable and never
qualified); the reference box integral of a
sub-tagged window reference ((b), D2a) is read as `ft_model / REF_A` with REF_A = E_in +
E_straddle / 2, the bracket [ft / (E_in + E_straddle), ft / E_in] and the decision-218 marker
(0.05 validated class / 0.10 new). Statuses stamped into process-library.json
(`QualificationStatus`, `SpatialQualification`, `LibraryQualified`): PendingQualification ->
Qualified ((a) on every dense trace, the identity, the gate) -> WindowValidated ((b) inside the
marker on >= 1 window); a failing criterion -> Failed. Smoke (R1b, local, 2 ranks, no
qualification): the synthetic device-plan coupon of `test_cluster_signature_geometry`
(3 conductors, 72 sources, coarse fabricated / thin meshes) through the whole pipeline at
p1 / p2 — the identity holds on every trace and class (<= 3e-8; the readers' conventions),
closure / p-stability fail by 1-40 % at these orders (expected: a 10 % p1 -> p2 energy change)
and the status is Failed; with the transmon's gate record the model is untestable -> Failed.

**Which coupon matrices act in the device correction (decision 299, 2026-10-04; from
`SurfaceResponseOperator::GetElectrostaticResponse` and `ElectrostaticSolver` ApplyResponse).**
The fixed-trace (ft) interface energy of a target interface is `E_out(R) + sum_patches w_p
t_p^T Q_fab,surf t_p`: the device's own raw thin energy within R of its edges is DROPPED and
replaced by the FABRICATED coupon's within-R energy under the device trace — the thin twin's
surface matrix does not enter. The thin twin acts through its DOMAIN matrix alone: the
fixed-trace domain correction `1/2 w t^T (Q_fab,dom - Q_thin,dom) t`, the self-consistent (sc)
operator `K + P^T (Q_fab,dom - Q_thin,dom) P` (ApplyUneliminated), and the fixed-flux (ff)
transform `F_thin^+ ...` built from the two DOMAIN matrices (the ff surface energy is
`t^T Q_fab,surf t` evaluated on the fixed-flux trace). The surface defects `Q_fab,surf -
Q_thin,surf` are assembled (`surface_defects`) but their only consumer
`SurfaceResponseOperator::GetEnergyCorrection` has NO caller in any driver — recorded cleanup
follow-up: remove it or give it a test; not changed in the decision-299 branch. Consequences:
(i) the decision-282 thin-side (F) criterion (the correction's surface p-stability, which
failed on both S1p v3 coupons by 3-13 % of E_fab, decision 293) judged a quantity with no
path into ft / ff / sc — the decision-293 rationale "the thin p-dependence is meant to cancel
the device's" is corrected to "unused": the thin SURFACE p-dependence is the decision-66
log-divergence at the recorded 2-nm sheet cutoff and is recorded as information
(`ThinSurfacePSteps`); (ii) the (F) thin-side criterion of record is the DOMAIN
twin-consistency `|dE_dom,p5 - dE_dom,p4| / E_fab,dom,p5 <= 0.02` per dense trace (the
device's box domain energy is finite and p-convergent, so the domain defect must be; the
residual is the thin-side error of a p5 device run corrected with the p4 library, in the
closure's unit; both S1p v3 coupons read <= 1e-4); (iii) a device-side domain reading is
recorded as information (`qualify/device_box_energies.py`: the registration device's thin
mesh with the placed boxes sub-tagged in the volume — one attribute per (per-box inside /
straddle status, material) combination so overlapping boxes are read exactly —, one solve at
p4 and p5 on that one mesh, the box domain energy bracket `[E_in, E_in + E_straddle]` and its
p-step next to the twin's domain energy under the device trace; not a criterion: the device
mesh resolves the sheet edges at micrometres against the twin's 2 nm, and a sheet-edge
field's domain energy converges like h^1; on the S1p c0 mesh no tetrahedron lies inside the
3.95-um-thick box, one uniform split gives 5 % of the box volume inside).

**Legacy-contract aliases (USER decision 283, 2026-10-03).** The transmon's 3-edge cluster
`spatialedgecluster_edgecount-3_da0179c64591` (claims-only key 7c4b31a894f9, two mirror
instances) is kept as a LEGACY-CONTRACT key: its accepted coupon (the decision-236
straight-continuation geometry) stays in the accepted library `transmon-r1p9-folded-refined-
corners` until the next transmon rebuild, its v3 context (16.9 um^2 of foreign ground beyond R,
0 fictitious; closure shift expected <= 0.1 point) recorded, not rebuilt. Under v3 its key is
35b2900ae07f (contract 3: one foreign ground piece), so the accepted library would read it
Missing; instead the library lists an EXPLICIT alias on the legacy model
(`LegacyContractAliases: [{Key, ContextDigest, Reason, Context}]`, `Key` = the v3 hash,
`ContextDigest` = `SpatialSupportContextDigest` = sha256 of the serialised `{Box, Context}`,
recorded per contract-3 feature as `SpatialSupport.ContextDigest`;
`signature_library.legacy_contract_alias(feature, reason)` writes the entry from a manifest).
The library load refuses an alias on a model that is not a Signature-keyed
`SpatialEdgeCluster`, a key / digest that is not 64 hex, an alias of the model's own key or one
key listed twice. The matching pass (`RunGeometryIdentification`) looks an alias up by the
feature's EXACT hash only after the normal key match fails, then `ResolveLegacyContractAlias`
fails closed unless the feature's context digest equals the alias's AND the feature's
claims-only key (`SpatialSupport.ClaimsKey`: the claims canonicalised alone, recorded by the
identification — the Portions of a contract-3 signature are serialised in the frame minimising
Box + Context and cannot be re-keyed by removing those members) equals the legacy model's
Signature key; the feature is then PLACED in its claims-only canonical frame
(`Features[].ClaimsFrame`: the legacy model's Edges and basis points live there; the v3 frame
may differ by a rotation / reflection — on the transmon it did, and the continuation-ownership
clipping with it) and the match is recorded: `Match.LegacyContract` (model, key, digest,
claims key, reason, the recorded Box + Context, the placement frame) and `Match.Note` on the
feature, `LegacyContract` on the version-1 requirement record (Status Exact: the legacy coupon
is applied), `Summary.LegacyContract` + `Summary.Counts.LegacyContract` /
`TotalEdgeLengths.LegacyContract` in the preflight inventory (the features are counted in Exact
as well), `SurfaceResponse.Diagnostics.LegacyContract` in the operator record (through
`ResponseCorrectionData.legacy_contract`, carried by the geometry cache under
`LegacyContract`, optional) with a warning. Never a fallback: a feature whose hash is not a
listed key is Missing as before. Verified on the transmon with the accepted library + the alias
(`device-plan-coupons-20261003/r1/transmon/library-alias`): the 3-edge features 8 / 12 resolve
to `spatialedgecluster_edgecount-3_da0179c64591` with `LegacyContract`, every feature matched
to the same model in the same placement frame as the accepted preflight, Counts Exact 3,099 /
35,296.81 um unchanged, the patches CSV BYTE-IDENTICAL (b362840fcc34) and the
ContinuationOwnership record identical; a wrong digest aborts the preflight. Unit tests
`SurfaceResponseIdentificationLegacyContractAlias` (resolution, the three abort paths, the
claims frame, the broadcast / manifest record) and the cache round trip in
`SurfaceResponseOperatorContinuationOwnership`. Decision 286 extended the alias into a POLICY:
an accepted-library model that v3 re-keys is kept through an explicit, recorded alias until
that library's next rebuild — the transmon's 4-edge `spatialedgecluster_edgecount-4_78e4c00ca560`
(its claims box has a ground concave corner 0.158 R outside the x0 / y0 faces; the two-sided
T2 of decision 285 (2) grows the box 0.5 R on both faces and keys 19.8 R of foreign ground
edge) is the second alias; the transmon stop condition then reads Counts Exact 3,099 /
LegacyContract 4 / Missing 0 with the patches CSV byte-identical.

**Block (b) phase 1 (decision 303; `curved-clusters-20261005/DESIGN.md` sections 2-4 + AMENDMENT
1 A1 / A4, 2026-10-05) — the quantum near-match, per-case span-cap allowances, arc context.**

*Quantum near-match of cluster keys (DESIGN section 4; ruling: max quanta 4).* A
`SpatialEdgeCluster` signature's TOPOLOGY KEY is the signature with every number of
`Portions[].P / Arc / Gap`, `Context[].P / Arc / Gap`, `Box` and `Vertices[].P / TurnDegrees`
replaced by null (entry ORDER preserved: the canonical serialisation sorts the entries, a
permutation is another key) and its PARAMETERS are those numbers in traversal order — lengths on
the 1e-6 R grid (the unit Gap components and the Box included) and angles on the 1e-6 deg grid
(`SplitSignatureParameters`, with the JSON path of every number). Two clusters of one topology
key whose numbers agree within `kClusterQuantumNearMatchMaxQuanta` = 4 quanta describe ONE
geometry at the resolution of the grid: the same design cell at different chip positions rounds
~3-6 % of its coordinates differently by sub-quantum float noise (the three stage-2 loop-end keys
284d6c2b5b66 / 9e103a0f291c / 20ac3e14a128 differ by exactly one quantum in 17 / 8 of their 290
numbers; the nearest real near-key is >= 2,000 quanta away). Every quantum threshold is read
HALF-QUANTUM INCLUSIVE (decisions 287 (b) / 288 (2); the phase-1 review's MINOR-1, decision
317): an on-grid difference of k quanta computes to k +- 1e-9 in floating point, so a feature
matches iff max |delta| <= 4 + 1/2 quanta (`WithinClusterQuantumNearMatch`) and two models or
two allowances are one geometry iff max |delta| <= 8 + 1/2 (`ClusterQuantumDuplicate`;
`kClusterQuantumInclusiveMargin` = 0.5, recorded as `Conventions.ClusterQuantumInclusiveMargin`
and `QuantumNearMatch.InclusiveMargin`); exactly 4 / 8 on-grid quanta match / are refused,
5 / 9 are not (unit tests through the matcher, the allowance resolver and the Python mirror).
`SignatureDeviation` of two clusters = max |delta| / (4.5 q) (<= 1 matches; no mirror
orientation: the chirality is folded into the canonical key), `ClusterSignatureQuantumDifference`
gives the largest |delta| in quanta and the differing paths. The matching pass (`LibrarySignatureIndex::Match`, by topology
key, the nearest model, ties by rank then name) records a near-match as
`Features[].Match.QuantumNearMatch {ModelKey, FeatureKey, MaxDeltaQuanta, DifferingNumbers
{Count, Paths}, MaxQuanta, Rule}` with a `Match.Note`; the status stays Matched / Exact (the
model's coupon is applied in the feature's own canonical frame with M = identity, the A10 /
A10-extended checks unchanged); `Summary.Counts.QuantumNearMatched` /
`TotalEdgeLengths.QuantumNearMatched` and `Summary.QuantumNearMatch {Features, Length,
MaxQuanta, Keys [{Model, ModelKey, FeatureKey, MaxDeltaQuanta, DifferingNumbers, Features}],
Rule}` in the inventory (the FeatureKey -> ModelKey map the thin-run guard and the library
tooling read; the guard itself reads the requirement records' Status, Exact for a near-matched
key), `SurfaceResponse.Diagnostics.QuantumNearMatch` in the operator record (through
`ResponseCorrectionData.quantum_near_match`, carried by the geometry cache under
`QuantumNearMatch`; cache version 7 -> 8, a stale cache is refused). Library side: two cluster
models of one topology within 2 x 4 = 8 quanta are refused at load ("two models, one geometry"),
so no feature can be within 4 quanta of two models; the requirement records group near-matching
features into one record whose Hash is the lexicographically smallest member's
(`RepresentativeSignature` of clusters: the lexicographically smallest SERIALISATION — the
design wrote "smallest key"; deterministic, mirrored by `representative_signature`; review
MINOR-7) and list the other members' keys under `NearKeys` (the grouping is single-linkage over
`SignatureDeviation` <= 1, as for the translational types: three variants in one window could
chain to 8 quanta in principle — none on the census; review MINOR-8); `signature_library.py`
mirrors the comparator (`split_cluster_parameters` / `cluster_quantum_difference` /
`signature_deviation` / `within_cluster_quantum_near_match` / `cluster_quantum_duplicate`) and
writes `NearKeys` on its models. The near-match is carried by the MANIFEST (`Features[].Match.
QuantumNearMatch`, `Summary.QuantumNearMatch.Keys`) — the patches CSV keeps its `Model` /
`Feature` columns and stays byte-identical for every exact-matched library (the transmon stop
condition); design MINOR-3's "ModelKey + FeatureKey in the CSV" is NOT implemented, ratified
by decision 317 (review MINOR-2). The library-side grouping is per window (the planner runs
one device config): `NearKeys` is populated only when two variants fall in ONE window; across
windows the build list names the representative (the loop end 1b26671c9080 serves
e52df843d197; a census-level grouping tool over `cluster_quantum_difference` is a follow-up,
review MINOR-3 / decision 317). Legacy-contract aliases stay EXACT-hash (decision 283,
unchanged); a permuted entry order stays Missing (the residual knife-edge, recorded).
`Conventions.ClusterQuantumNearMatchMaxQuanta / ClusterQuantumInclusiveMargin /
ClusterQuantumNearMatch`.

*Per-case span-cap allowances (DESIGN section 3 (a) + A4).* `Solver.Electrostatic.
ResponseCorrection.SpatialSupport.SpanCapAllowances[]` / `Solver.SurfaceResponseCorrection.
SpatialSupport.SpanCapAllowances[]` (beside `Library` / `PatchConstruction`; the design note's
`Boundaries.SurfaceResponse` path does not exist — review MINOR-6) `{ClaimsSignature,
SpanCapOverR, Label, Reason, Approval}` carries the plan
span cap of ONE approved closed feature that cannot be split (decision 244 (i)), keyed by the
EMBEDDED claims-only signature of the feature (verbatim as the inventory exports it:
`Features[].SpatialSupport.ClaimsSignature` of a refused cluster — SpanCapRefusedGrowth /
UnboxableFeature —, the Missing placeholder's claims + `"Unboxable": true` accepted as well);
the Label is the recorded hash prefix, never compared. Resolution ORDER: the cluster's
claims-only signature -> the allowance (the quantum near-match comparator, the nearest within 4
quanta, ties by Label) -> box growth with span_cap = the allowance -> the contract-3 key — the
allowance cannot depend on the boxed key it produces and one entry serves every window instance
of the feature. Records: `SpatialSupport.SpanCapAllowance {Label, SpanCapOverR, Reason,
Approval, MatchedQuanta, DifferingNumbers, Rule}` per feature (and `SpanCapOverR` /
`ExceedsSpanCap` read against the allowance), `SpanCapAllowance` on the requirement record (the
planner copies it onto the plan coupon; `device_coupons` passes its SpanCapOverR as the
generator's `--support-span-cap` unless a CLI option names the hash), `UnusedSpanCapAllowances`
(Labels matching no cluster of the run: a warning, never an abort). Fail closed at the start of
the identification: SpanCapOverR below `kSupportSpanCapOverRadius` (an allowance never lowers
the cap), two allowances within 8 (+ 1/2) quanta ("two allowances, one geometry"), a
ClaimsSignature with a Box / Context or without Portions. The generator's claim radius
(`normalize_geometry`) is half the case's recorded span cap (default 16 R -> 8 R,
byte-identical). Measured on the stage-2 census (`curved-clusters-20261005/phase1/recensus`):
with the loop end's allowance (20 R) and 005bec161f6d's (21 R) exactly the four Unboxable keys
box (MatchedQuanta 0 / 0 / 1 / 1 / 1 on the five SCT windows), everything else unchanged; the
loop end boxes to two contract-3 keys of one topology (S1 / S1p 1b26671c9080, S2p / S3p / S4
e52df843d197, 8 numbers one quantum apart) that one model resolves (Exact / QuantumNearMatch).
The allowance VALUES are ratified by decision 317 / STAGE2-PLAN AMENDMENT 3 (the recorded
`Approval`) on the measured need (`phase1/fixes/allowance-sweep`, caps 16..24 R): the loop end's
claims box spans 17.96 R and first boxes at cap 19 R (Box span 18.46 R), 005bec161f6d's spans
19.00 R and first boxes at cap 20 R (19.50 R); both allowances are NON-BINDING — Box and key
unchanged at allowance + 2 R and at every cap from the first boxing one to 24 R. An evidence
config that names its `Output` inside the evidence tree is re-run ONLY from a copy with a
scratch `Output` (review MINOR-4): a re-run in place overwrites the stored inventory.

*Arc context (DESIGN section 2 + A1).* `DevicePerimeterDistance` (the A10 check extended to the
context) reads a device segment lying on a fitted arc (`Segments[].Arc`) on the ARC of its
circle (`ArcChordDistance`: the point's in-plane projection inside the chord's angular interval
reads the radial residual + out-of-plane offset, outside it the nearer chord end; a chord
collinear with the centre falls back to the chord), exact and chord-independent — an arc
context entry's end cut at a box face lies on the fitted circle, up to the chord sagitta
(1.6e-3..5.4e-2 R on the census) from the device polyline; `placement_check.py` mirrors it
(`device_perimeter_distance`) and gates every placed context entry end (`A10-context`).
`VerifySpatialEdgesInSignatureFrame` compares a model's chord rows with its signature's arc at
`kArcFitToleranceOverRadius` + 2 quanta, inclusive (A1 (5)). The builder
(`cluster_signature_geometry`) fixes an arc entry's ENDS first — the serialised ends, each
snapped onto a box face within one quantum (inclusive) and then onto the end of a neighbouring
piece of the same conductor within the arc-fit tolerance + 2 quanta (inclusive; the straight
neighbour's device vertex wins, between two arcs the end of the smaller radius; the two must be
each other's nearest such end, two candidate points fail closed; closed circles are not
snapped) — and REBUILDS the circle through them (`rebuilt_arc`: centre = the point of the
perpendicular bisector nearest the serialised centre, radius |a' - c'|, the side of the
serialised midpoint, equal angular steps), so every chord vertex is concyclic to double
precision, the straight neighbours are untouched and the arrangement closes; the rebuilt
circle's deviation from the signature's is asserted <= the fit tolerance + 2 quanta
(`chorded_entries`) and recorded with every snap in `coupon.json Geometry.JointSnaps[]`
({Piece, End, Class Face | ArcJoint | ArcArcJoint, To, DistanceOverR} and one `Arc` record per
arc {ArcDeviationOverR, CentreShiftOverR, RadiusShiftOverR, Chords, Centre, RadiusOverR});
straight coupons are unchanged (legacy byte identity). Step 0 (F0-a) of the same block: the
builder's gap perpendicularity test admits the quantisation bound 2 x (sqrt(2) / 2 q + 2 q / L)
of a serialised straight row and (decision 317 MAJOR-1, option (a) = the design's letter)
re-derives EVERY straight row with 0 < |tangent . gap| <= bound as the exact perpendicular of
the chord with the serialised sign (`GapRederived`; no threshold inside the bound — the S1p
loop end's portion 30 at |tangent . gap| = 9.999999999995e-07 is re-derived like its 17 other
oblique-by-rounding rows), keeps a row with |tangent . gap| == 0 exactly bitwise (every row of
the eight group-A coupons built by stage 2.2: their sources regenerate byte-identically) and
leaves an arc CHORD's analytic radial gap as computed (not a serialised number); beyond the
bound it fails closed. Near-collinear line / line joints (|turn| <= `JUNCTION_TANGENT_ANGLE`
1e-4 rad, same class) merge in both plan-view boundaries (A3 (1)).

**Placement of a contract-3 (device-plan) model (decision 285 (4), R1b).** The cluster patch of
a model whose Signature carries `Box` + `Context` records in its provenance the model's
`SupportBox` and its chain pieces (the `Chain: true` context entries, arcs chorded) in the
patch's local frame in units of R (geometry cache version 7), and the placement verifies A10
EXTENDED TO THE CONTEXT: the model's `Context` must equal the matched feature's and the two
ends of every context entry, placed by the patch frame, must lie within 1e-3 R of a device
edge (`DevicePerimeterDistance` on the identification's segments; an arc entry's ends are
device vertices, its chords lie on the fitted circle); a mis-keyed library or a placement
frame defect aborts, and the count of verified ends is printed with the matching summary.
LIMITATION (R1 final review MINOR-3): the ends are tested against the RAW mesh polyline, so an
arc context piece cut at a box FACE ends on the fitted circle, up to the chord sagitta off the
polyline (the S1p loop end: 0.054 R), and a legitimately keyed arc-crossing coupon would be
refused at placement — fail closed and loud; today every arc cluster is unboxable, so nothing
in production reaches it; lifting it means testing such ends against the fitted arc. The
VERTEX OWNERSHIP of rule B4 (`ApplyContinuationOwnership`, after the translational cells): a
vertex feature's patch (corner / junction / endpoint coupon: coupon depth 0, no claims) whose
vertex lies within 1e-3 R of a chain piece END of a contract-3 support, inside that support's
box (the signature box placed by the patch frame; the same plane) is owned by the coupon —
the device-plan coupon contains the real corner, so the corner coupon's band on the chain
arms would be counted twice: the patch keeps weight 0 (once; a vertex inside two boxes lists
both owners, review MINOR-1 (b)), its face distance is recorded and a vertex closer than R to
a face is flagged `ArmOutsideBox` (an arm partly outside the box, MINOR-1 (a)) with the lost
length recorded, `LostArmLengthOverR` = max(0, 1 - FaceDistanceOverR) of the first owner's box
(R1 final review MINOR-2: the part of the corner's R window beyond the box that no coupon
corrects once the vertex patch has weight 0 — C2p feature 10: 0.66 R of one arm); a vertex on
another cluster's claims is never on a chain (the chain stops there, decision 285 (1)) and a
legacy (claims-only) model owns no vertex. `ApplyContinuationOwnership` takes R explicitly
(`matching_radius`, MINOR-7) next to the continuation tolerance. Recorded under
`Diagnostics.ContinuationOwnership.Vertices` {Count, Shared, Records [{Kind: Vertex, Patch,
Feature, Model, Origin, FaceDistanceOverR, ChainEndDistanceOverR, ArmOutsideBox,
LostArmLengthOverR, Owners}]}
in the preflight and the operator (the owned patches are skipped like wholly owned cells). A
contract-3 placeholder without basis points (a signature-only library: the preflight of a
Missing key) takes its support bounds from the Signature's box, R above and below the plane
(`BySupport[].FromSignatureBox`), so the dry run's ownership records are complete; a library
model keyed by an `Unboxable` signature is refused at load (no coupon exists for it). Census
(R1b, signature-only placeholders of the v3 keys): S1p — the 4 L2 90-deg corners (features 7,
8, 10, 12) owned by the two 19-edge boxes, corner 10 by both (as rule B4 predicted), 56
context ends verified; C2p — 3 corners owned (one shared, one at 0.342 R from a face with an
arm outside), 54 ends verified; transmon — no v3 model, nothing owned, the accepted patches
byte-identical. Unit test `SurfaceResponseOperatorContinuationOwnership` (vertex ownership:
the chain corner owned once by two supports, interior chain points / other planes / other
boxes / legacy models untouched, idempotent, the Diagnostics entry; the v7 cache round trip).

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
  1.95 / 2 / 2.05 um beyond the separation value itself;
* **recorded knife-edge residual (synthetic `gap-bend-r{50,250}-2Rminus-step{5,15}`, review
  fix-4):** two concentric 8 um bars at the design gap 2R - 1e-3 R along a polyline bend whose
  offset sides are MITRED at the joints (vertex displacement g / (2 cos(step / 2)), an exact
  offset path) carry two separations at once: the chords are parallel at exactly g = 3.998 um
  (< 2R), while the joint vertices lie on circles whose radii differ by g / cos(step / 2) =
  4.033 um at 15 deg and 4.002 um at 5 deg per vertex (> 2R; 3.99815 um at 1 deg, < 2R). The
  identification decides the pair from the fitted-circle separation (the curve reading, so
  that mid-chord dips of a coarse polyline never create events — DS-SCT-001's 4 um CPW gaps at
  R = 2 um) and reads the 5 / 15 deg bends as isolated edges; the facing gate measures the
  chord distance 3.998 um < 2R and flags 79-139 um (A8 FAIL; the 1 deg variants and every
  2R / 2R + 1e-3 R variant pass, the exact-2R ones as `AtExactly2R`). The two readings of one
  polyline differ by g (1 / cos(step / 2) - 1) = 8.6e-3 g at 15 deg = 17x the 1e-3 R margin
  of the design value below 2R, so the polyline is not one geometry at the scale of the
  knife edge; no consistent reading resolves it and neither is wrong. Not in any gate set
  (the identity set records the manifests); a design at 2R - 1e-3 R on a coarse polyline
  bend is exactly what the knife-edge census flags (perimeter within 1 % of 2R).

**Knife-edge census (decision 82(4), 2026-09-25; R = 1.9 um by decision 88(1), 2026-09-26).**
The strict rule is kept and the matching radius is chosen off the common layout dimensions:
the default of the identification tools (`preflight_config.py`, seed
`seeds/preflight_seed_r1p9.json`; the R 2.1 seed of decision 82(4) is kept for comparison
runs) is R = 1.9 um — thresholds R 1.9, 2R 3.8, 10R 19 um; the 2D locality study
(2026-09-26) found R = 1.9 / 2.0 / 2.1 um equivalent within 0.2 % at the knife edge and far,
so the choice is one of cost: at 1.9 um the 4 um separations stay isolated edges (DS-SCT-002:
895 um within 1 % of 2R = 3.8 um against 577.6 mm at R = 2.0 um and 7.0 mm at 2.1 um; 2
stack coupons); the library's `MatchingRadius` stays authoritative and R = 2 um libraries are
unchanged. So that any design can check its R, the
manifest reports `Identification.KnifeEdgeCensus`: the perimeter length with another perimeter
point (3D; the same chain beyond the self-pair neighbourhood; runs sharing a vertex excluded)
at a distance within `KnifeEdgeBandRelative` = 0.01 of R and of 2R, the chain length whose
windowed bend radius lies within 1 % of `StraightBendRadiusOverR` R, the two-run vertices
whose implied sagitta (c / 2) tan(turn / 4) on the SHORTER run lies within 1 % of
`JointNoiseSagittaOverR` R (the joint noise rule, `JointNoiseSagittaOverR`), whose turn lies
within 1 % of `ArcMaxJointTurnDegrees` (the arc rule's cap, `ArcMaxJointTurnDegrees`) and
whose LONGER run, read as a chord at the joint's turn, has a sagitta within 1 % of
`SagittaOverR` R (the mesh-coarseness diagnostic, `ArcSagittaOverR`), each split into the
below / above sides (samples every 0.5 R); and (USER decision 184 (4), 2026-10-01; the
stage-0 E6 finding: the DS-CTX-003 loop-end clusters went from 42 to 48 edges at R 1.85 while
every distance band stayed quiet) the **cluster-composition band**
(`KnifeEdgeCensus.ClusterComposition`): a cluster's membership is the outcome of the whole
cluster machinery (events, cores, joins, claims, the extension), not of one distance, so the
band is measured the way the E6 study measured it — the same perimeter (the input as
received, the vertex classification at R kept) is identified again at R (1 - 0.01) and at
R (1 + 0.01), silently (no census of its own, root rank only), every `SpatialEdgeCluster` at
R is matched to the cluster at the other radius sharing the most claimed perimeter with it,
and its composition is unchanged iff the match has the same `EdgeCount` and the same set of
member vertices (corner / endpoint / junction mesh vertices). The claimed length at R of the
changed clusters is reported on the Below / Above side (Total = both) with the counts
`ClustersBelow` / `ClustersAbove` of `Clusters`, and in the identification log. Clusters
that exist only at the other radius are not counted (their perimeter's interaction
distances are the distance bands' business). DS-SCT-001 at R = 2 um: 7,247 um within 1 % of R
and 7,175 um within 1 % of 2R (the 2 / 2 / 2 um flux lines). The library continuity gate
(`coupon_library.py continuity`, `qualification-gates.json` LibraryContinuity Version 3, USER
decisions 117(3) and 2026-09-28 on the decision-117 review MAJOR-2: a pair / stack model
whose consecutive separations are all >= 2R (1 - 0.01) responds like the isolated-edge model
per basis function at matched local positions — gated per edge on the energy-weighted
matched offset within 1 % (`MaximumRelativeOffset`) AND per matched basis function within
5 % (`PerBasisFunctionLimit`: a single hat carries the mesh noise of two independently meshed
coupons, +-1-2 %, so the per-hat limit is looser than the aggregate; the recorded 1.985 R
stack / pair pass with worst hats -1.8 / -2.3 / -1.6 and -1.3 %); an edge with zero matched
basis functions is Failed, never NotApplicable) runs on every written process library.

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
| `ParallelEdgeCluster` | quadrature on every side, side factor 1 / n; origin = the sample's foot on the canonical first side (the sample itself on that side); anchors = the feet on the first side of every conductor label | u = from the origin to its foot on the last side (the local lateral: a stack following straight-like bends turns with them — a feature-wide frame placed the meandering DS-SCT-001 flux-line stacks up to 49 deg off, found by the A10 placement audit), v = mean process normal | `(s1 - s0) x w_q / n / CouponDepth` |
| `ConvexCorner`, `ConcaveCorner` (sharp or rounded), `Endpoint`, `Junction` | one patch at `Frame.Origin` (the vertex or the virtual corner of a fillet) | `Frame.Axes` (below; w = the site's own signed process normal, decision 266) | model weight (1) |
| `SpatialEdgeCluster` | one patch | a model carrying its `Signature` is built in that Signature's canonical frame and is placed with the identity map (m maps to `F.origin + F.axes^T m`; `F.axes[2]` = the cluster's own signed process normal, decision 266, so the model's w > 0 is the feature's vacuum on a flipped plane too; its stored `Edges`, when present, are verified at library load to lie on the Signature's portions — straight portions as segments, arc portions on their circle — within the signature tolerance, each edge on ONE portion whose relabelled `Conductor` is the edge's and whose `Interfaces` set is the set mapped to the edge's `InterfaceSlot` (a mirror-symmetric geometry with asymmetric labels lands on the portions but is refused), `VerifySpatialEdgesInSignatureFrame`, fail closed); a legacy model without a `Signature` maps a model-frame point m to `F.origin + F.axes^T M.axes (m - M.origin)`, M from `CanonicalClusterSignature` of its stored straight edges | model weight (1) |

**Longitudinal cell of a translational patch (`SurfaceMortar` strip; decision 152).** The
weight is the patch's dimensionless MEASURE (the fraction of the model's `CouponDepth` it
integrates, times the side and model factors); the strip over which the surface mortar
projects the device trace is a LENGTH, carried separately as the patch's longitudinal cell:
the interval of its portion that its quadrature point integrates, as offsets along the patch
AxisW from the origin in mesh units, ordered begin <= end (`longitudinal_cell`; `StripBegin`,
`StripEnd` of the dry run in manifest units; `LongitudinalCell` of the geometry cache, version
3; version 4 adds the patch's `Feature` / `Stretch` provenance for the ownership check of
decision 224, version 5 the spatial cluster patch's `Claims` and the library's
`MatchingRadius` for the ownership record's continuation class of decision 236, version 6
the patch's mesh `Segment` and own-edge `EdgeOffset` for the per-side continuation
ownership of decision 252 — a version-5 cache dropped the segment, so cached patches were
classified by the abutment branch alone; older caches refused; version 7 (decision 285 (4))
the contract-3 cluster patch's `SupportBox` and chain pieces `Chain` for the vertex
ownership of rule B4). The cells of one portion tile it exactly — cell q = [cumulative
Gauss weight before q, + w_q] in order of increasing position, so every cell contains its
point (`LongitudinalQuadratureCells`; Gauss cells centred on their points would not tile) —
and every patch built from a quadrature point carries the full cell whatever its weight
factors (both sides of a pair, every side of a stack, the two co-located patches of a
first-order split). The mortar samples the strip in ceil(cell length / mortar resolution)
slices at the slice midpoints, each slice a transverse projection at the contour resolution
(oversampled by `MortarOversampling`), and averages them with equal weights (lift and
transpose alike): the patch coefficient is the trace's LENGTH-AVERAGE over its cell, so the
patch energy `weight x c^T Q c` is the Jensen lower bound of the strip's energy and the sum
over a portion is a surface functional of the field over the whole matched perimeter. A
patch without a cell ({0, 0}: 2D, spatial and vertex models, explicit configuration syntax)
is one cross-section at its origin. A pair or stack patch's AxisW follows the PARTNER side
(AxisU points at the sample's closest foot on the other side), so on the slow tapers and
sub-noise polyline bends the identification classifies as pairs AxisW is not parallel to the
sample's own segment (1 - |cos| of 1e-6..1e-3): the arc cell is projected onto AxisW by
Dot(tangent, AxisW) — the sign orients it, the magnitude shortens it — because the
cross-section perpendicular to AxisW at the projected offset contains the segment point at that
arc offset (the foot lies in the plane perpendicular to AxisW through the origin), so the slices
sweep exactly the sample's cell of its own segment (`LongitudinalCellOffsets`); a frame below
|cos| = 0.95 (`kLongitudinalAxisCosineTolerance`, the paired-edge topology's facing threshold)
is not a translational frame and fails closed. Before this fix the dimensionless weight was used as the
strip length in the solver's nondimensional mesh coordinates (`mortar_longitudinal_subdivisions`,
`longitudinal_coordinate`), i.e. a strip of l_cell x Lc / `CouponDepth` centred on the Gauss
point (Lc the mesh's characteristic length): 3.8 cells on the transmon (Lc 4 mm, `CouponDepth`
1055 um; ~3 slices overlapping the neighbouring cells), one cross-section where that length is
below the mortar resolution (small islands); the fixed-trace and self-consistent trace maps both
change with the fix wherever the strips were not already single slices.

**Domain-boundary exclusion (decision 258, 2026-10-02; `FindDomainBoundaryExclusions`).** A
coupon's trace coupling is undefined beyond the device domain: its basis points sample the
device field on the coupon contour, and a contour point outside the mesh has no value (the
mortar-resolution lookup at the model's first basis point and the point location both fail
closed — the S1p thin runs 3a / 3b aborted in the constructor on two isolated-edge patches,
216 / 222, where the L1 / L2 leads meet the window's left cut x = 290 at ~13 deg to its
normal). This happens wherever a metal edge meets an ARTIFICIAL domain cut — a window cut,
a chip outline (decision 193 coverage item) — within ~R |sin theta| of it, theta the angle
between the edge and the cut's normal: the cross-section at the edge end sticks out of the
cut by R |sin theta| on one side (0.42 um on S1p; a perpendicular lead, theta = 0, sticks out
nothing and its first cell ends exactly on the cut). RULE: a library-placed 3D patch any of
whose PLACED COUPON POINTS lies outside the device mesh is NOT applied and is recorded as a
`DomainBoundary` exclusion. The placed points are the model's basis (contour) points and its
conductor references in the patch frame, at the patch origin cross-section and, for a
translational patch with a longitudinal cell, at both cell ends moved
`kSignatureParameterToleranceOverRadius` x R = 1e-3 R inward along AxisW (the only
tolerance, dimensionless in R: a cut lead's first cell ends exactly ON the cut and the
S1p bottom-wall cluster strips end on the wall within 8e-6 um of a 1e-6-rad frame tilt;
their sample slices lie at least half a slice inside; without the end sections a longer
cut portion whose origin section is inside would pass the preflight and abort on its
strip slices — the unit case's second cell). ONE containment test decides for the
preflight and the solve: every rank locates every tested point in its local mesh with the
operator's own `ElementPointLocator` (its reference-space tolerance 1e-9 on linear
simplices, the inverse transformation otherwise, the routing box tolerance) and the found
flags are OR-reduced over the communicator, so the decision is the partition's union —
identical for any rank count and for both of the operator's locator paths
(`ElementPointLocator` below 64 ranks, `FindPointsGSLIB` at 64 and above or with
`PALACE_RESPONSE_USE_GSLIB_POINTS`), which afterwards locate the APPLIED patches' points
only. An excluded patch keeps weight 0 (the operator skips it like a wholly owned cell; the
dry run writes it with Weight 0, its unscaled QuadratureWeight and cell, no new column) and
its portion stays tiled (the A7 identity holds; the audit reads the record). RECORD
(`Identification.Diagnostics.DomainBoundaryExclusions` of the preflight manifest, the
operator's `Diagnostics` in the metadata, one summary line in both logs with the test's
wall time): per patch the feature, model, origin, CELL (begin / end offsets and
`CellLength` = the length left uncorrected: quadrature weight x portion), PORTION
(`Segment`, `S0`, `S1`, `PortionLength`), tested / outside point counts, the outside point
nearest the mesh with `NearestDistance` = its distance to the nearest element bounding box
(exact under an axis-aligned cut, a lower bound under a tilted face, 0 when the point
lies inside a sheared element's box); totals `CellLength`, `PortionLength`, `ByFeature`.
INVENTORY: `Summary.Counts.DomainBoundary` (patches) and `Summary.TotalEdgeLengths.
DomainBoundary` (the CELL length) next to Missing — it is not a library gap, the matched
totals are untouched — so the B1 gap bound can include it; `Summary.DomainBoundary` adds
`Features` and the portion sum `PortionLength` (information: on S1p the two excluded
cells are 1.0145 + 0.9637 = 1.978 um of 2 x 3 patches on portions of 3.652 + 3.469 =
7.121 um; the "~7.12 um" of decision 258 is that portion sum). The excluded cells REMAIN
inside `Exact` (they are matched; S1p: 1.978 of Exact 2,364.778 um): `Exact` +
`Interpolated` + `Missing` is the identified length and `DomainBoundary` is the part of
`Exact` left uncorrected — a consumer (the B1 bound) adds `DomainBoundary` to `Missing`
and must NEVER sum Exact + Interpolated + Missing + DomainBoundary. FAIL CLOSED (decision
260, review MAJOR-1): the exclusion is for a cut THROUGH a placed coupon; a candidate whose
FIRST conductor reference at the origin section (the metal-edge point, which lies on a mesh
face for every correctly placed patch, a metal edge on a chip-outline face included) is
not located, or none of whose tested points is, is a misplaced or mis-scaled coupon (a
library in the wrong units, a radius beyond the domain) and the test aborts naming the
patch (0-based, as in the record and the dry run), model, point and coordinates; after the
exclusion at least one applied patch must remain when there were candidates
(`MFEM_VERIFY`). Explicit (configured) 2D patches are never excluded: a point of theirs
outside the mesh aborts as before, and every point of an applied patch the operator cannot
locate still fails closed, the message naming the point, its coordinates and the patch
(model). The candidate set is the placed patches after the continuation ownership; the
operator additionally omits the vertex patches subordinated to an overlapping cluster box
before the test — a KNOWN LIMIT of the preflight (S1p tests 2494 patches in the
preflight, 2490 in the operator, the same 2 excluded): that subordination lives in the
operator's pairwise spatial-support loop with its overlap abort, and the preflight does
not reproduce it. The totals cannot differ: a subordinated vertex patch is spatial, its
`longitudinal_cell` is {0, 0} and it carries no portion (`Segment` -1), so it contributes
0 to `CellLength` and `PortionLength`; only `Count` could show one more entry in the
preflight record. The log summary prints patch indices 1-based like the ownership
summaries and names the base; the record and the dry run are 0-based. Cost: all ranks test all points (3 sections x (basis + references)
per translational patch): transmon 10,366 patches / 2.98 M points in 1.8 s on 1 rank, S1p
2,494 / 763 k in 0.24 s on 2 ranks; the transmon preflight digest 9ada660bf6e4, record and
dry run are unchanged (no exclusion: its domain is far from the metal). Unit test
`SurfaceResponseOperator domain-boundary exclusion` (`test-domainboundary.cpp`): a lead cut
by a domain face tilted out of the metal plane (theta 24.2 deg) on 1 and 2 ranks, the
default and the forced-GSLIB locator path, the notch fail-closed case, and the misplaced
coupon (reference off the mesh, every point off the mesh, no applied patch left) aborts.

**Conductor-consistency gate (decision 277 (A), 2026-10-03;
`SurfaceResponseOperator::ApplyConductorConsistencyGate`, SOLVE TIME ONLY).** A spatial
coupon holds its metal cross-sections on the box faces at the conductor potentials (the
trace mesh's conductor vertices, fixed in the surface mortar) while the device trace is
imposed on the free knots around them. Where the coupon's metal is not the device's — the
S1p 19-edge junction-lead coupons continue the island protrusion's edge straight 3R past
the device's corner, inventing 20-27 um^2 of island metal over device gap at 10-19 V
(decisions 271-277: an MS-specific surplus of x1.92 the reference, 39 points of the window
MS) — the device potential at the coupon's metal differs from the conductor's and every
energy of the patch under the device trace is wrong. The gate PROBES the device potential
at every conductor vertex of the trace mesh on the PROCESS PLANE (w = 0: the coupon's metal
bottom, which lies on the device's metal sheet, a Dirichlet surface for real metal) — the
probe points are appended to the patch's sampled point list after the conductor references
and the trace walks stop before them — and forms, per applied spatial surface-mortar patch
and per excitation, MaxRatio = max |V_device(knot) - V_device(conductor reference)| /
Normalization with Normalization = max(Amplitude, `kConductorConsistencyAmplitudeFloor` =
1e-3 x ExcitationPotential), Amplitude = max |trace coefficient| of the patch incl. its
conductor states (= |state| for a two-conductor coupon without Gibbs overshoots;
dimensionless, defined for single-conductor coupons and for excitations whose conductors
sit at one potential) and ExcitationPotential = max |V| of the excitation (its largest
terminal potential). The floor (decision 279, MINOR-3) keeps a patch whose trace amplitude
is a near-zero fraction of the excitation (noise over noise) from being excluded
spuriously: such a patch carries at most 1e-6 of a unit-amplitude patch's energy (energy
~ amplitude^2), a defective coupon is still excluded at the excitation that fields it (the
exclusion is sticky), and the record marks it `FloorApplied` (the transmon's far corner
patches read amplitudes 1e-6 to 1e-3 V under a 1-V excitation; the S1p flagged coupons
~1.1 V).
MaxRatio > `kConductorConsistencyTolerance` = 0.02 EXCLUDES the patch like a DomainBoundary
cell: weight 0 from that excitation on (fixed-trace energies, the self-consistent operator
and its right-hand side, ModelCatalog weights), `surface-response-patches.csv` rewritten,
the record `SurfaceResponse.Diagnostics.ConductorConsistency` of palace.json (every tested
patch of every excitation: Amplitude, State, ExcitationPotential, Normalization,
FloorApplied, MaxDeviation, MaxRatio, MaxRatioOverState, the worst knot's vertex /
conductor / device point, OffPlaneMaxRatio, AdjacentMaxRatio, ClaimLength, CellLength,
Excluded; `ExcludedPatches`, `Count`, `ClaimLength` = the excluded clusters' claimed portion
length left uncorrected, `CellLength`; lengths and coordinates in mesh-file units = the
device coordinates, like the DomainBoundary record) and a log summary. The record is
complete on the solver side; the B1 COVERAGE BOUND IS NOT APPLIED BY THE SOLVER: the B1
consumer (the analysis tabulation, e.g. `stage1-20261002/thin/tools/tabulate_s1p.py` for
S1p) must add `ClaimLength` to the Missing and DomainBoundary lengths it bounds (decision
279, MINOR-1). CALIBRATION
(`coupon-accuracy-assessment-20260913/conductor-consistency-20261003/gate/`): real metal
reads exactly the Dirichlet value — S1p c0 field at the 31 real-metal plane knots of the two
19-edge coupons <= 7e-7, the operator on the S1p thin smoke <= 1.5e-9 on 6 corner patches,
0 on the transmon-like unit case — while the fictitious blocks read 0.20-0.43 of the
amplitude (0.22-0.48 of the state; 3 + 4 knots); 0.02 is x10 below the weakest flagged knot
and >= 3e4 above the real-metal maximum, and a 5 % potential mismatch would already move
the 15 %-residual MS closure by several points (so the tolerance is not larger). THE
REAL-METAL RESIDUAL IS RESOLUTION-DEPENDENT (decision 299 (2), from the A6 record and its
erratum): the probe point sits a roundoff distance delta off the device sheet edge (the
placement's frame arithmetic: S1p box 3's worst point reads (582.5750002, -116.4999995) for
the nominal (582.575, -116.5), 0.2-0.5 nm; the transmon's 4-edge end knots 0.05-0.75 nm) and
the FE potential there differs from the Dirichlet value by delta x the FE gradient at the
sheet-edge singularity, which grows with the resolved singularity: ~p^2 in the order (the
transmon 1.9e-8 -> 7.7e-8 from Order 1 to 2; S1p box 3 2.2e-7 -> 6.3e-7 from Order 3 to 5 on
c0) and x1.3-1.4 per AMR cycle (box 3 6.3e-7 at c0 -> 1.09e-6 at c2 -> 5.4e-6 at c7 on T, 4.8e-7
-> 4.3e-6 on P4; box 1 1.3e-7 -> 7.8e-7; the corners <= 9e-9). The solve-time exclusion 0.02
is unchanged (4e3 above the most refined real-metal reading); the (F) qualification's
acceptance (`spatial_qualification.GATE_MAX_RATIO`) is 1e-4 on EVERY supplied device solve,
to be supplied at the production orders p4 and p5 on the registration (c0) mesh (the
decision-282 value 1e-6 was a c0 / Order-3 statement, exceeded by box 3 from c2 on): 1e-4 x
2^(k/2) stays below 0.02 for k < 15 halvings (every production sequence to date <= 8 cycles,
>= 12x margin) and lies >= 2e3 below fictitious metal. RECORDED
ONLY, never gated: the conductor vertices OFF the plane (the metal top rows, 0.1 um into
the device gap for a thin device: the normal field x the thickness, 1.3-2.8 % on S1p real
metal) and the ADJACENT free knots (sharing a trace-triangle edge with a plane conductor
vertex in its column, the 50-nm trench-floor row; their mortar coefficient against the
conductor value reads the near-edge field x 50 nm: up to 0.27 of the amplitude on the
transmon's real metal (10-edge junction coupon), 0.11-0.17 on its corners — they cannot
discriminate real metal from the S1p blocks' 0.21-0.36 and are information). The plane
row is load-bearing: the OffPlaneMaxRatio of real metal (1.3-4.8 % on S1p) already sits
above the tolerance, so a coupon frame whose w = 0 did not coincide with the device's
metal sheet (a thick-metal device keyed on its metal top) would be excluded wholesale —
loudly and recorded, not silently. FAIL CLOSED: a conductor of the trace mesh without a
vertex on the process plane aborts at construction naming the model and conductor (its
cross-section cannot be probed). RECORDED UNTESTABLE (`UnprobedModels` [{Model, ModelIndex,
Reason}] in the record + a warning naming the model; decision 279, MINOR-2): a spatial
model applied collocated, and a spatial surface-mortar trace mesh without any conductor
vertex — legitimate cases exist (a finite-impedance coupon's metal knots are free by
construction; the ring-path corner models without `ZeroTraceIndices` of the unit libraries),
and such a model has no metal cross-section at a fixed potential for the gate to compare,
so it is not gated; neither production library has one (every applied transmon / S1p
spatial patch is probed). The PREFLIGHT cannot
evaluate the gate (no device trace): the manifest states so under
`Summary.ConductorConsistency` (`Evaluated` false, the tolerance, where the record lives);
its digest, inventory and dry run are unchanged. Transmon (MEASURED, decision 279: the
T-cont configuration on the same mesh and library at Order 1 / 2, MaxIts 0, 2 ranks;
`conductor-consistency-20261003/gate/transmon/local-solve-*`): Count 0 of 33 patches / 250
plane knots; 31 patches read MaxRatio <= 3e-13 at both orders (exact Dirichlet values),
the two 4-edge clusters 7.7e-8 at Order 2 / 1.9e-8 at Order 1 (their cross-section end
knots lie 0.05-0.75 nm outside the device lead edges: the library's 1e-7-um coordinate
rounding, read as the local gap gradient x the offset, >= 2.6e5 below the tolerance); 12 /
13 far corner patches with amplitudes under 1e-3 V are FloorApplied; the offline census of
the 250 plane knots against the mesh's metal plan finds 242 on device metal of the same
conductor as their reference and the 8 four-edge end knots within 0.75 nm of it (on no
other conductor); every output CSV of the Order-2 run is byte-identical to main's binary
(e8bc64b3bf) and palace.json differs only in the new record — the preflight digest
9ada660bf6e4, record and CSV unchanged;
S1p thin smoke: exactly the two
19-edge patches excluded (0.411 / 0.395 of the amplitude at the fictitious-block face
knots), 122.2 um of claims left uncorrected, raw / C identical, every other model's
fixed-trace energy identical. Unit test `SurfaceResponseOperator conductor-consistency
gate` (`test-conductorconsistency.cpp`): a hand-placed two-conductor coupon across the gap
between two islands (consistent: applied, MaxRatio 0; a conductor vertex over the gap:
excluded, recorded, sc operator 0, sticky; the same under a 1e-4-scaled excitation; the
same coupon under an excitation it barely sees (amplitude 1e-6 of the 1-V potential):
FloorApplied, not excluded; conductor off the plane: aborts; no conductor vertex: recorded
untestable, applied), 1 and 2 ranks.

**Vertex-feature frames** (`Frame` of the manifest, shared by the library builder): n = the
site's OWN signed process normal (substrate -> vacuum, oriented by the length-weighted mean of
its incident runs' normals; decision 266, frame rule of (b) — on a flipped plane n = -z, never
the sign-canonical n_ref as such); corner: x = the first arm away from the (virtual) corner, the arms
ordered so that the second is counterclockwise about n (a corner is its own mirror image),
y = n x x (right-handed); endpoint: x = the arm, y = +-(n x x) toward the gap; junction: x = the
canonical first arm of `CanonicalJunctionSignature`, the arm angles measured about n, y =
+-(n x x) so that the canonical arm order proceeds counterclockwise in (x, y). A legacy junction model (absolute `ArmAngles`) is mapped by its own
canonical order (first arm angle theta, orientation): u = cos(theta) D - sigma sin(theta) (n x D),
v = sigma (n x u), sigma = +1 when both orientations agree.

**Curvature (C4, decision 108; the two-dimensional axisymmetric path uses the same rules).**
A curved feature (`CurvedEdge`, `CurvedSameConductorGap`, `CurvedDifferentConductorGap`,
`CurvedSameConductorStrip`) that no Signature-keyed model matches exactly is modelled by the
**curvature family** of its straight analogue: the anchor (the straight model at the
feature's separation) and the library's coupons of the curved topology with `Kappa` =
R / rho (rho the edge radius, for a pair the INNER edge's) and `Convexity` (that of the
model's first edge e1: metal inside the bend), at kappa = 1 / `RadiusOverR` and the
feature's `Convexity`, with `FindCurvedLibraryModel`'s recorded rule: an exact node is that
coupon, kappa <= 1 / `StraightBendRadiusOverR` the linear combination of the anchor and the
node AT that kappa, otherwise the cubic Lagrange interpolant on the four nearest nodes; the
combination is one runtime model `<anchor>@<convexity>-kappa<k>-<rule>` whose matrices are
the weighted sum of the nodes' matrices rescaled to the anchor's `CouponDepth` (Lagrange
weights may be negative), patched like the straight analogue (`Match.Note` records the rule,
the version-1 record `CurvatureFamily` the anchor / nodes / weights, `Status`
`Interpolated`). Never silently straight: a curved feature the family cannot model (no
anchor, no coupons of its convexity, kappa above the largest node, an anchor mapping other
interfaces, `Mixed` convexity) is unmatched with the reason in `Match.Note`; a curved pair
whose (class, separation) group has no family is unmatched the same way. **First-order term
of a straight-like feature** (`BendRadiusOverR` >= `StraightBendRadiusOverR` on an
`IsolatedEdge` / straight pair): the recorded linear rule, evaluated per portion where the
curvature is — a straight-like feature spans bends of both senses and straight legs (a
meander), so one feature-wide kappa would carry the wrong sign and magnitude. Every
quadrature point of a portion with a nonzero signed turn splits between the anchor (model
weight 1 - a) and the family node at kappa 0.1 of the turn's convexity (weight a =
kappa_local / 0.1, kappa_local = R |turn| / length, clamped to [0, 1]; for a pair the inner
edge's kappa: an outer-side sample reads rho_inner = rho - s, and the node's convexity is
e1's — the far side turns in the opposite sense): two co-located positive patches whose
assembled matrices are the linear interpolant (identical to blending the matrices for the
fixed-trace and self-consistent closures, equal to second order in the node - anchor
difference for the fixed-flux transform, which is per model). A missing first-order node
leaves the anchor alone on that portion, counted and warned (`Curvature:` summary line;
never silent); the A7 audit's per-interval rule (sum of quadrature x model weights = 1) and
the per-patch weight formula hold with two models per interval. The A10 placement audit
resolves a family runtime model to its anchor and checks, for every curved pair patch (blend
or first-order node), that e1 sits on the side of the coupon's recorded `Convexity` (the
inner circle of a gap / the outer circle of a strip for `Convex`, read on the claimed fitted
arc under e1).

**Angle-interpolated corner family (USER decision 121 (C), 2026-09-28; `MatchCornerFamily`).**
A SHARP corner feature (`ConvexCorner` / `ConcaveCorner`, `CornerRadiusOverR` 0) that no
Signature-keyed model matches within the angle tolerance is modelled by the **corner family**
of its convexity: the library's sharp corner coupons of the same topology, interfaces and law
(records `Angle` / `AngleDegrees`, `Convexity`; one box basis for the whole family, checked at
match time) interpolated in the TURN t = 180 - `AngleDegrees`. The interpolation variable is
the turn because the coupon response is smooth in the arm direction and t = 0 is the straight
edge, where the first-order regime is anchored: the family's anchor is the straight edge
through the corner box (`Angle` 180, built by the corner generator on the same basis — the
"180-deg anchor"; no device corner has that angle, a corner being a joint that is not noise).
Rule (mirrors the curvature family's): an exact node within `SignatureAngleToleranceDegrees`
-> that coupon; t at or below the smallest coupon turn (the largest coupon angle, e.g. 165
deg -> 15 deg) -> linear between the anchor and that node (first order in the turn: the
corner excess of a small turn is O(t), as the curvature family's first-order regime in
kappa); otherwise cubic Lagrange on the four nearest abscissae (anchor included); t above the
largest coupon turn (an angle sharper than the sharpest coupon) -> unmatched with the reason;
a first-order corner without the anchor -> unmatched with the reason; a rounded corner ->
refused (per-radius coupons and `CornerRadiusInterpolation` only, no angle family) — never
silently straight, never the nearest coupon. The runtime model is the BLEND of the nodes'
matrices, entry by entry (Lagrange weights may be negative; corner coupons have no coupon
depth: the weights apply as they are), named `<base>@corner-angle<deg>-<rule>` on the NEAREST
node (its conductor references), patched once with weight one in the FEATURE's frame at the
feature's angle like an exact coupon; `Match.Note` records the rule, the version-1 record
`CornerFamily` the base / angle / turn / convexity / rule / nodes and weights, `Status`
`Interpolated`. The A10 placement audit resolves the runtime model to its base coupon at the
INTERPOLATED angle (`placement_check.resolve_model`), so the arms checked are the feature's.

**Corner trace basis rule (corner-family review 2026-09-29, supervisor decision 132;
`generate_corner_response.py` record `TraceBasis`, `palace/models/cornertracebasis.cpp`).**
The review's root cause of the non-90-degree MS / MA over-correction: on the lane-2 layout
(8 knots per box ring at the square's corners and side midpoints, angle-independent) a metal
arm at 75 / 105 / 120 / 150 / 165 deg crosses a box ring BETWEEN a free knot and a PEC knot,
so the free hat is nonzero on the PEC part of the box contour — conflicting Dirichlet data
in the coupon solve and a near-singular field in the 2 nm MS / MA layers (the fabricated
Q_MS of that knot 6-7 orders above its clean value), added as the absolute fabricated
surface form at runtime while the domain defect (fab - thin) cancelled it. The rule: on every
box ring that meets the metal (z = 0: the thin sheet and the slab foot; z = MetalThickness:
the slab top) the knots are the two crossings of the metal arms with the ring (PEC, zero set),
`MetalInteriorKnots` = 1 knot at equal perimeter-arc-length fractions of the metal arc between
them (PEC) and `FreeKnots` = 5 knots at equal fractions of the free arc (free); the knots are
ordered by ROLE within the ring (convex: free 2 .. 5, first crossing, metal interior, second
crossing, free 1 — the lane-2 order of the 90-degree CONVEX node, which the rule reproduces
byte-identically; concave: first crossing, free 1 .. 5, second crossing, metal interior — the
concave 90-degree node differs from lane 2's by construction, 5 free knots on the free
quadrant arc instead of 1), so
every node of a family has the same knot count (72), the same zero set (convex 1-based
29-31 / 37-39, concave 25 / 31-33 / 39 / 40) and like-to-like free knots whose POSITIONS
vary smoothly with the angle: the entrywise blend is well posed. Fractions are perimeter arc
length from (-R, 0) counterclockwise (`Fractions` `PerimeterArcLength`). The box corners that
are no knot are SLAVE vertices of the trace triangulation (`trace-vertices.csv` columns
`parent_a` / `parent_b` / `weight_a`: the linear interpolation in the fraction between the
neighbouring knots; `MortarVertex::ForEachBasis` in the surface mortar; the Maxwell path
refuses them), so the trace surface stays on the box; rings that do not meet the metal keep
the fixed layout. No free hat has support on the PEC part of the box contour
(`free_hat_pec_support`; the library gate `library_basis_gate.py` for every coupon type).
Runtime: the library load FAILS CLOSED on any spatial corner coupon whose metal arm crosses
a box ring at no PEC knot (`CheckCornerBasisCrossings`: the lane-2 layout at 120 / 165 deg is
refused; 90 / 135 / 180 pass); `MatchCornerFamily` requires the rule on every node (a lane-2
coupon is unmatched with the reason), checks one knot semantics (rule parameters, basis
size, `ContourGroups`, `ZeroTraceIndices`, the fixed rings, every node's files against the
rule at its angle) and CONSTRUCTS the runtime basis at the feature's angle
(`BuildCornerTraceBasis`: the base's fixed rings + the metal rings by the rule + slave
corners, carried by the runtime model as `constructed_basis_points` / trace vertices /
triangles and by the geometry cache). `ZeroTraceIndices` semantics on the electrostatic
path (review m5): a PEC knot's trace is zero whatever the device potential at its point (a
slab-top knot lies in the thin device's air); the surface mortar and, since this fix, the
collocated lift enforce it (RATIFIED, supervisor decision 134 (a), 2026-09-29: physically
right and consistent with the mortar; quoted device numbers use `TraceCoupling`
`SurfaceMortar` with `MortarOversampling` 2 — the protocol's lift — and none of the
corner-basis-fix block's quoted numbers was produced with the collocated lift). Knot
coincidence (review m3): `KNOT_COINCIDENCE_FRACTION` / `kKnotCoincidenceFraction` = 1e-6 of
the perimeter — a free or metal-interior knot within it of a fixed-layout fraction
k / RingSize (a box corner or side midpoint) takes that fraction exactly and a corner within
it of any knot gets no slave vertex, so the smallest triangle a constructed basis can
contain has an edge of 8e-6 R (area >= 4e-6 R h, six orders above the mortar's
degenerate-triangle threshold `area <= 1e-14 max(1, L^2)`), while a snap moves a knot by at
most 8e-6 R (15 pm at R = 1.9 um); crossing knots are never moved (they stay on the arm for
the gates; a crossing within the band suppresses the corner's slave). The former 1e-9 was a
floating-point identity tolerance and left arbitrarily thin (non-degenerate) slivers. The
ring selection of `CheckCornerBasisCrossings` (commit 5bd5089a7): a ring is judged when it
has PEC knots AND every knot lies on the perimeter of the centred square |x|, |y| <= its
half width (a corner-coupon box ring); the commit message's motivation ("the concave family's
metal rings reach neither the left nor the bottom side") does not describe the test — the
code is right, the message is not, recorded here rather than rewritten. Recorded limitation
of the rule (held-out interpolation,
corner-basis-fix-20260929): a free knot passes a box corner at some angle (its hat straddles
the corner on one side only), a kink of the matrix entries in the angle; the cubic stencil
75 / 90 / 105 / 120 of the convex 82.5-deg held-out angle contains the passing of free 1 (the
knot next to the second crossing) through (-R, R) and its fabricated MA interpolates to
-6.3 %, the linear first-order regime at 172.5 deg to -2.7 / +3.1 % (convex / concave);
SA / MS interpolate to <= 0.2 % everywhere (participation-referenced).

**Kink-aware interpolation: trace basis EVENTS, segment connectivity, per-side nodes
(corner-qualification block 2026-09-29, supervisor decisions 134 (c) / 135; `CornerBasisEvents`,
`SelectCornerFamilyStencil`, `BuildCornerTraceBasis` connectivity, Python mirror
`corner_family_interpolation.py`).** Measured on the generator's trace surfaces (hats sampled
on the bands next to both metal rings at theta -/+ 1e-3 deg; a no-event angle gives 4e-5): the
"kink" above is a JUMP. Every angle at which a knot of a metal ring passes a vertex of the
fixed layout on the neighbouring rings — a box corner OR a side midpoint — flips the diagonal
of a band quad in the perimeter-ordered merge (`connect_rings_by_fraction`) and the hats
change by O(1) (0.42 at 90, 0.50 at 135 / 141.34 / 153.43 / 158.20, 0.63 at 111.80); a box
corner's slave vertex does not help because the knot swaps its perimeter order with the slave.
The general rule is therefore: **a stencil never straddles a geometric event of the trace
basis.** The events of the rule (M = 1, F = 5), computed from the affine relation between
every knot's fraction and the second crossing's (`CornerBasisEvents`; identical in
`corner_family_interpolation.basis_events`), in the family range [75, 180]: CONVEX knot-corner
passages 90 (free 1 / 3 / 5 and metal 1 at the four corners: the fixed layout), 135 (the
second crossing at (-R, R)), 153.434948822922 = 180 - atan(1/2) (free 2 at (-R, -R));
side-midpoint passages 141.340191745910 = 180 - atan(4/5) (free 1 at (-R, 0)) and 180 (the
anchor's crossing at (-R, 0)); CONCAVE corner passages 90 (free 3 at (R, R), metal 1 at
(-R, -R)), 135 (free 2 at (R, R) and the crossing), 158.198590513648 = 180 - atan(2/5) (free 5
at (-R, R)); midpoint passage 111.801409486352 = 180 - atan(5/2) (free 5 at (0, R)) and 180.
(1) **Segment connectivity** (`TraceBasis.ConnectivityAngleDegrees`, generator option
`--connectivity-angle`, requirement `Geometry.ConnectivityAngleDegrees`, part of the coupon
Id): the bands next to the metal rings are merged in the perimeter order of the rule's layout
at the connectivity angle (a knot's key is its ROLE's fraction there, a slave's its corner;
`connectivity_keys`), while the knots sit at the coupon's own angle. The triangulation is then
one and the same for every coupon of a SEGMENT (the coupons sharing the connectivity angle),
which removes the side-midpoint flips inside a segment (hat difference 6e-5 across 141.34 with
keys 144.2), and it can never fold: every band triangle lies between two consecutive rungs on
one box side because the corner rungs are ties, and the key order equals the perimeter order
whenever no knot-corner passage lies between the two angles (refused otherwise, both
languages). A knot-corner passage cannot be removed: it is a jump between the two
neighbouring segments' triangulations, and the recorded coupon at such an angle (the
perimeter-ordered tie merge) equals NEITHER side (checked at 90 / 135 / 180). (2) **Per-side
nodes**: the family carries one coupon per side at every knot-corner passage (the same knots,
the segment's connectivity: e.g. `90-` with keys 82.5 for [75, 90] and `90+` with keys 112.5
for [90, 135]) and the anchor rebuilt with its segment's keys; the connectivity angle of a
segment is recorded as its midpoint (convex 82.5 / 112.5 / 144.2 / 166.7, concave 82.5 /
112.5 / 146.6 / 169.1). (3) **Node density** (supervisor): every event-free segment has at
least four nodes at most 15 deg apart, so the interpolation is cubic everywhere — convex 75,
80, 85, 90- | 90+, 105, 120, 135- | 135+, 144, 150, 153.435- | 153.435+, 159, 165, 180;
concave 75, 80, 85, 90- | 90+, 105, 120, 135- | 135+, 143, 150, 158.199- | 158.199+, 165, 172,
180 (16 coupons per convexity). (4) **Stencil rule** (`SelectCornerFamilyStencil`, used by
`MatchCornerFamily`): the segment whose node range contains the angle strictly inside gives
the Lagrange stencil on its nodes nearest to the angle (cubic on four, else the highest order
the segment supports: quadratic / linear); an angle within `SignatureAngleToleranceDegrees` of a
node is EXACT — with several coupons at that angle the legacy one (no connectivity record: the
tie triangulation at its own angle, e.g. the recorded 90 / 135 coupons, the lane-2-identical
90) is preferred, else the coupon of the lower-angle segment (the Signature index ranks the
candidates legacy 0 / lower segment 1 / upper segment 2 before the name tie-break, so the
key-based exact match agrees); the runtime basis of an interpolated corner is constructed
with the segment's connectivity; outside the node range refused (no extrapolation; the
first-order regime is the anchor's segment, cubic like every other). The node tolerance is
closed on the exact side and open on the interior side by the SAME floating-point difference
(|node - angle| <= tol exact, > tol interior; the range refusals follow the exact test), so an
angle at the tolerance boundary of a node is exact or interior, never "in no segment"
(qualification review 2026-09-29 m1; identical in the Python mirror). FAIL CLOSED AT LIBRARY
LOAD (`ReadProcessLibrary`, `CheckCornerFamilySegments`; decision 137 (1), review m3): for every
corner family of the library (the sharp trace-basis-rule coupons of one topology, interface
set and boundary law) a segment whose connectivity angle or nodes lie across a knot-corner
passage, two coupons at one angle in one segment, or overlapping segments refuse the library
before any corner is matched (`SelectCornerFamilyStencil` asserts the same precondition); the
per-node check that a segment node's recorded trace triangulation is the rule's for its
connectivity angle stays at match time (`MatchCornerFamily`, it reads the trace meshes). A
LEGACY family (coupons without connectivity records — the recorded corner-basis-fix libraries,
the lane-2 90-degree coupons of the transmon libraries) is never interpolated (reason: "corner
family built without segment connectivity records … rebuild") while its exact matches stay
valid (a coupon at the device angle is self-consistent), so the verified 90-degree-only device
libraries keep working. The version-1 record `CornerFamily` carries
`ConnectivityAngleDegrees`. The per-side coupons serve only as stencil nodes; the held-out
check of the qualified family (`corner-qualification-20260929/`) judges the rule in both forms.
**Production libraries carry the legacy tie coupons at 90 / 135 / 180 beside the per-side
nodes** (supervisor decision 140 (1), qualification review M1): the exact-node preference then
keeps the recorded corner-basis-fix 90-degree coupons for every exact 90-degree corner (the
convex one = the transmon-verified model, byte-equal matrices; the concave legacy tie 90 is
the corner-basis-fix rule node with 5 free knots on the free arc, which differs from the
lane-2 transmon concave 90 by construction), while the per-side coupons remain the stencil
nodes of their segments. Measured on the held-out trace
of the recorded matrices, the three convex 90-degree coupons (identical knots, three band
triangulations) differ by **MA_fab +10.85 % (90-, keys 82.5, vs the verified legacy 90; MA
defect +14.25 %)** and by at most 0.7 % on every energy for 90+ (keys 112.5; MA_fab -0.00 %):
this is the **MA representation sensitivity of the coarse box basis** (8-knot rings, R/3 ring
spacing: the MA response depends at the 10 % level on how the band right above the metal top
edge is triangulated), recorded as an UNCERTAINTY of the corner family's MA — an open item
(a finer band basis: an extra ring near the metal top or 16 metal-ring knots, decision 140 (3)),
not a defect of the rule (every stencil is consistent with one triangulation).

**Refined trace basis: `RingLayout` `AllRingsFollowMetal` (corner-basis refinement
2026-09-30, USER decision 161; `generate_corner_response.TraceBasisRule` / `REFINED_RULE`,
C++ `CornerRingLayout` / `RefinedCornerTraceBasisRule`; evidence
`corner-basis-refinement-20260930/REPORT.md`).** Under the option-(c) held-out trace (USER
decision 149 (6)) 37 of the 40 recorded corner coupons failed the 10 % self-check (decision
154): the matrices were exact, the coarse basis could not represent a trace nonzero on the
metal rings. Measured cause (Phase 1, 7 recorded coupons, the coarse interpolant of the (c)
coefficients solved directly against the fine held-out solve): (i) the fraction MISMATCH
between a rule ring and its fixed neighbours (spurious band gradients in the 0.05-um band
next to the metal: the concave thin SA +195 %); (ii) the lateral RESOLUTION of the R / 3 ramp
next to each crossing with free knots 0.67-1 R apart (the convex 90 node = the fixed layout,
no mismatch: thin SA -14 / MA -12 %, fab MA -13 % are pure resolution); (iii) the VERTICAL
structure right above the metal top over the metal arc and below the trench, which the
recorded self-check could not see because its fine reference shared the basis's seven
z-levels (the recorded reference was off by up to 16 % on the fabricated MA and 7 % on MS,
measured by refining its levels). The refined rule: EVERY ring of the box — the outer rings
at -R, -R/3, -OveretchDepth, 0, MetalThickness, MetalThickness + k OveretchDepth for k in
`ExtraLevelsAboveOverOveretch` = [1, 4] (0.15 and 0.30 um at the recorded process: the mirror
of the trench ring, and the ring that resolves the trace right above the metal top over the
metal arc, without which the concave family's fabricated MA read +7 %), R/3, R, and the two
inner cap rings — carries the SAME angle-dependent knot fractions: the two crossings,
`MetalInteriorKnots` = 5 at equal fractions of the metal arc, `FreeKnots` = 9 with
`FreeKnotGrading` [1/3, 2/3] (a knot at R/3 and one at 2R/3 along the perimeter from each
crossing on the free side, 5 at equal fractions between; ACUTE CONCAVE nodes — block (b)
DESIGN A9 family 4, decision 318: a concave wedge's free arc is (2 - cot theta) R and the
unscaled grading exhausts it below 56.3 deg, so on a free arc shorter than
`FreeKnotGradingReferenceFreeArcOverR` = 2 - cot 60 deg = 1.4226497308103743 (the concave
60-degree node's, the family's sharpest qualified node) the two graded distances are scaled
by FreeArc / Reference — the sharper node keeps the 60-degree node's layout proportions with
the same slots, like-to-like indices and zero set; the generator writes the key only on
records where the scaling is active and an absent key reads as that default, so every record
written before (free arcs at or above the reference) is unchanged and the family stays on one
rule; `FreeKnotGradingScale`, `kFreeKnotGradingReferenceFreeArcOverR`; STATUS, decisions 325 /
328, 2026-10-05: the scaling rule stands; the 48.75-degree held-out coupon's failure to
finalize — two scaled free hats of the z = -R cap ring had no active boundary DOF in the
fabricated p4 solve — is a COUPON MESH defect: the coupon prescribes hats by nodal
interpolation, and the cap rings' inner knots (9.4 nm apart at 60 degrees, 7.4 at 48.75) were
unresolved by the 300-nm far mesh at every concave node <= 60 degrees, the published 60 node
included; round 2 adds the generator-side trace resolvability gate (`trace_resolvability.py`,
5 p^2 active boundary nodes per free hat) and the knot-gap mesh sizing of
`mesh_corner_coupon.jl --trace-mesh`, a corner-coupon RECIPE change (new cache keys; the
published caches are not rebuilt); the DEVICE runtime is unaffected: the SurfaceMortar lift is
an L2 projection on the coupon's trace triangles with the device field sampled by
FindPointsGSLIB, independent of the device mesh's DOF layout
(test-cornerbasisrefinement.cpp CornerRefinedRuleAcuteConcaveDeviceMortar); the branch is NOT
for merge until round 2 is reviewed); PEC = crossings +
metal-interior knots on the two metal rings only (14 of 176 knots); the box corners are slaves
on every ring; each cap is a fan from a centre slave at the mean of the cap ring's two crossing
knots
(the recorded fixed cap's fan diagonal gave the centre that value too; a fan from a ring
vertex is degenerate as soon as two consecutive dense knots share the apex's side); every
band is a regular column grid with one diagonal orientation. Consequences: NO events (no
fixed vertex is ever passed; a knot passing a box-corner slave on every ring at once leaves
the interpolant continuous — measured 1-3e-5 at the MetalRingsOnly event angles, the
smooth-angle level; `CornerBasisEvents` returns the empty list), ONE segment per convexity
(no `ConnectivityAngleDegrees` — a coupon carrying one is refused at library load; no
per-side coupons, no legacy tie nodes: exact 90 = the single 90 node, the recorded 10.85 %
per-side MA spread is 0 by construction), the stencil the cubic sliding window on the four
nearest nodes (one-sided at the range ends: the family's nodes 75 / 80 / 90 / 105 / 120 /
135 / 150 / 165 / 180 — the 80 node added after the review of decision 167 so the 78 and 82.5
held-out angles are bracketed, held-out 78 / 82.5 / 112.5 / 142.5 / 172.5 / 176). Measured with
the 80 node (combined-verification-20260930/fixes): the gating option-(c) interpolation gate
passes both convexities, the 78-deg residuals fall clear (convex MA +0.449 -> -0.121 %), but the
convex MA at 82.5 reads -0.379 % (0.12 points under the 0.5 % gate; +0.362 before): a coupon MA
RESOLUTION floor, not the node spacing — the convex fabricated MA matrices change 2.7-3.3 %
between the coupon mesh factors 2 and 1 at 75-82.5 deg (SA 0.7, MS 0.2 %), i.e. the family's
MA is h-converged only to ~3 % at the acute end and each coupon's MA carries a few tenths of a
per cent of coupon-to-coupon scatter that no interpolant between nodes can remove. Open
follow-up (not blocking): a coupon-MA-resolution study (finer `lc_fine` on the corner coupons).
The transmon's two 90-deg corner models carry 1.5-1.6 % of its fabricated MA surface energy
(T p5; the concave 90 ~0), so a 3 % corner MA h-error moves the transmon MA by ~0.05 % — below
its error bars. The
planner (`prepare_surface_response_coupons.py --corner-trace-basis`) builds the refined rule
by default on SHARP corners only; a rounded corner (`CornerRadius` > 0) keeps the legacy rule
(the refined rule is qualified on sharp corners only; a refined request there is refused). The library load checks
every AllRingsFollowMetal coupon against the rule at its own angle (`CheckCornerRuleCouponFiles`:
the outer ring levels, every basis point, every knot's zero flag against `ZeroTraceIndices`,
the trace mesh's vertices / slave parents /
triangle set; fail closed; skipped for rounded corners, which carry no refined rule). The runtime constructs the basis of an interpolated corner by the
rule at the device angle as before (all rings, the centre slaves; the SurfaceMortar lift
through `MortarVertex::ForEachBasis`). Cost: 176 knots (72 before); the coupon's trace solves
are one operator with many right-hand sides. Recorded Phase-1 measurements (worst |error| of
the (c) self-check over the 7 coupons against the converged reference): the recorded rule
195 %, neighbour rings following alone 15.6 % (it uncovers the resolution error), all rings
with the recorded knots (M1 F5) 16.8 %, M3 F11 graded 7-level 30 % (fab MA against the
converged reference), the chosen rule 6.3 % — the thin DOMAIN of a concave coupon, a KNOWN
representation error common to every all-rings candidate (the coarse trace sits 2-4 % low in
the band -R/3 .. -OveretchDepth, linear in z against the hypot ramp next to the crossings,
and above MetalThickness + R/3, where no knot sits at the level the cutoff reaches 1; not
device-gated); every SA / MS / MA within 5.3 %; denser lateral grading (R/6 .. 2R/3, F 11-15)
and more far-side knots do not move it, more rings below the trench worsen MS (the graded
knots interpolate the hypot ramp linearly on those rings). **The held-out self-check now
judges the VERTICAL representation too (supervisor decision 2026-09-30, option (iii)):** the
reference surface (`heldout_reference_levels`, `HELDOUT_REFERENCE_RING_SIZE` 64) is
decoupled from the basis levels — the standard levels plus rings every OveretchDepth across
both R/3 ramps (from the metal top up to and including MetalThickness + R/3, where the (c)
cutoff reaches 1, and from the trench floor down to -R/3): measured convergence on the 7
recorded coupons S0 (7 levels) -> S1 (+ t + d): fab MA -5..-16 %; -> S2 (+ t + 2d, -2d): MS
-3.5..-6, MA -3..-6; -> S3 (+ t + 4d, -4d): MA -1.2..-2.9; -> S4 (+ 0.5, -0.4): domain +2..+3,
MS +4.8..+7.4 (not monotone: the linear-in-z interpolation misses the smoothstep); -> S5
(+ 0.4, -0.3): < 0.2 %; -> S6 (0.05-um spacing + t + R/3): domain +1.5..+2.3, MA +1..+2; -> S7
(0.025-um spacing): < 0.4 % on every energy — the recorded spacing is S6's (within 0.4 % of
S7). The reference is the same for every rule: a `--trace-basis legacy` rebuild reproduces the
recorded coupon's basis, trace mesh and zero set byte-identically but is judged against the
converged reference, not the recorded 7-level one. The option-(c) traces GATE the family's interpolation check
(`qualification-gates.json` Version 5, `GatingTrace`; USER decision 161 (2)); the band-trace
verdict is reported alongside.

**Library contract.** A model keyed by its `Signature` (the feature's canonical object, `Type`
included; `Signature.Type` must equal `Topology`) needs no version-1 geometry parameters of its
own; a model matches every feature of its topology whose parameters lie within the signature
tolerance (decision 85(1): 1e-3 R for lengths, 1e-2 deg for angles; the nearest model wins,
ties by name), and the library builder groups the feature instances of one topology within
the tolerance into one coupon at their representative signature (the version-1 `Requirements`
record carries `Signature`, `Instances` (feature instances), `DistinctSignatures`,
`ParameterSpread`, `ExactParameters`; the
`ParallelEdgeCluster` / `CurvedParallelEdgeCluster` `Geometry` carries `Edges`, `EdgeCount`
and, for the curved class, `BendRadius`); the curved classes (`CurvedEdge`, `CurvedSameConductorGap`, `CurvedDifferentConductorGap`,
`CurvedSameConductorStrip`) exist only as signature-keyed models and are patched like their
straight analogues along the curved portions; a cluster model keyed by its signature needs no
`Edges`, and when it stores them (the library builder: the Signature's straight portions as
edges and every arc portion chorded at `ClusterArcChordStepDegrees` / `ClusterArcChordMaxLengthOverR`,
`cluster_signature_geometry.portions_from_signature`) they must be the Signature's portions in
its canonical frame (Point = P x R, Interval along gap x normal, process normal +z): the
placement is then the identity, exact for arc clusters, whose chords could not be
re-canonicalised into arc portions (a chorded model re-canonicalised from its chords landed in
another frame: the 2394fdb0c failure class; the re-canonicalisation also chose, for the transmon's
180-deg-symmetric 10-edge JJ cluster, the frame rotated by 180 deg — geometry-identical, so the
placement audit passed, but the two leads are different conductors: the coupon's conductor-1 edges
sat on the conductor-2 lead, self-consistent SA -0.10 % on the device; the identity placement with
the label check above removes this class). A model's interface types are part of its key (a model mapping MA + MS + SA never
matches an SA-only feature). Runtime models are one per (library model, target interfaces by
slot; slot k = the k-th distinct target map of the feature's portions in sorted order).

**Patch dry run.** `palace --surface-response-preflight` builds the same patches without a
field solve and writes `surface-response-patches.csv` next to the manifest: `Patch, Feature,
Topology, Model, ModelIndex, Weight, ModelWeight, QuadratureWeight, SideFactor, CouponDepth,
Segment, S0, S1, Origin, AxisU, AxisV, AxisW, StripBegin, StripEnd` (manifest units;
`Segment` = the manifest segment index, `[S0, S1)` the portion from the segment's canonical
key origin; vertex and cluster patches carry `Segment` -1 and `CouponDepth` 0;
`[StripBegin, StripEnd]` the longitudinal cell as offsets along AxisW from the origin, the
cells of one portion tiling it, 0 / 0 for patches without a strip). The audit's gates A7: the patched
feature set equals the matched set; every portion of a matched longitudinal feature is exactly
one quadrature interval (sum of quadrature x model weights = 1) and `Weight = ModelWeight x
QuadratureWeight x (S1 - S0) x SideFactor / CouponDepth` with `SideFactor` = 1 / claimed chains;
vertex / cluster features carry patches without a portion whose model weights sum to 1; no
patch on an unmatched feature, an excluded segment or an excluded portion. With the
signature-only library built from the manifest's own features
(`examples/surface_response_identification/signature_library.py`) the covered length equals
`Totals.AssignedLength`: the whole perimeter minus the recorded exclusions.

**Placement audit (gates A10, `placement_check.py`; geometry only).** The coverage gates see
where a coupon is applied, not how it is oriented: a misplaced coupon (a version-2 cluster
model canonicalised to another frame, 2394fdb0c: SA +19 / MS +527 %) passes every A1 / A7 gate.
With the library the preflight ran with (`Library.Path` of the manifest or `--library`), every
matched cluster / corner / stack / pair model is mapped through the patch frame the dry run
wrote and must land on the feature's claimed portions within the signature parameter
tolerance (1e-3 R): a `SpatialEdgeCluster` model's `Edges` (or the chorded portions of a
Signature-only model) — every edge endpoint on a claimed straight sub-segment or on a claimed
arc (radially on the fitted circle, within the claimed angular range), and every claimed
straight endpoint / arc end on the mapped model polyline; a corner model's arms (+u and the
model `Angle` counterclockwise; a rounded corner's claimed chords on the fillet circle of the
model `CornerRadius`); a stack model's `Offset`s along u from the patch origin on the claimed
side of that offset (sides reversed for `Chirality` -1); a pair model's `Separation` about the
patch origin. Model lengths are library units times R_mesh / R_library. `A10-placement-*` per
class (features, checks, worst deviation / R, defects); a model the gate cannot evaluate is
listed under `A10-placement-evaluable`, never silently passed. On the lane-2 transmon preflight
the gate fails the pre-2394fdb0c run (worst 3.39 R on the 3-edge cluster) and passes the fixed
one (5.1e-7 R). Sensitivity per class (recorded): clusters and corners are read at the
signature tolerance (1e-3 R); stack and pair patches are read on the mesh chords with the
chord sagitta and the pair rule's own 5 % of the offset as allowances, so a stack / pair frame
defect smaller than 5 % of the offset (0.1 R at a 2 R offset) is invisible to the gate by
construction (the two stack defects it found were 0.7-2.3 R), and a slow taper wider than 5 %
about its mean fails by construction (the coupon sits at the mean separation).

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
  "Conventions": {"JointNoiseSagittaOverR": 0.05, "CornerRule": "...", "InteractionDistanceOverR": 2,
                  "ThroughVertexZoneOverR": 2, "ClusterBallOverR": 1,
                  "VertexJoinsClusterOverR": 2, "VertexWindowOverR": 1,
                  "ParallelCosineTolerance": 1e-8, "ParallelClassRule": "...",
                  "ArcFitToleranceOverR": 1e-3,
                  "ArcMaxJointTurnDegrees": 50, "SagittaOverR": 0.05,
                  "ArcRule": "...",
                  "ArcSampleSpacingOverR": 0.05, "ClusterArcChordStepDegrees": 5,
                  "ClusterArcChordMaxLengthOverR": 0.25, "ClusterGeometryRule": "...",
                  "SignatureLengthQuantumOverR": 1e-6, "SignatureAngleQuantumDegrees": 1e-6,
                  "StraightBendRadiusOverR": 10, "CurvatureWindowOverR": 1,
                  "PairSeparationToleranceRelative": 0.05, "PairSeparationSamplesPerInterval": 16,
                  "PairSeparationEstimate": "per sample: chord reading C = min over the two chains of the maximum sampled closest-point distance within max(R, min(local chord, PairChordWindowCapOverR R)) of the sample / its foot where the chain bends within PairBendProximityOverR R of it along the chain, R elsewhere; inscribed reading C / cos(turn / 2) with the larger local joint turn; interacting iff both < 2R; feature separation = mean C",
                  "PairConstancyWindowOverR": 1, "PairSampleSpacingOverR": 0.5,
                  "PairBendProximityOverR": 1, "PairChordWindowCapOverR": 2,
                  "PairCandidateReachOverR": 2.1, "SamplingMargins": "...", "CrossLayerReachOverR": 2,
                  "SelfPairNeighbourhoodOverR": 3.14159, "StackRule": "...",
                  "StackCompositionCap": 64, "ClusterExtensionWedgeCapDegrees": 30,
                  "ClusterExtensionClosureOverR": 1e-3, "ClusterExtensionRule": "...",
                  "SignatureParameterToleranceOverR": 1e-3, "SignatureAngleToleranceDegrees": 1e-2,
                  "SignatureMatching": "...", "MutualSidesOverhangOverSeparation": 0.32,
                  "KnifeEdgeBandRelative": 0.01, "FacingGateExclusions": "...",
                  "PlaneRule": "...", "PortRule": "...",
                  "Comparison": "strict less on the quantized grid"},
  "ReferenceProcessNormal": [nx, ny, nz],
  "Features": [ {"Id": k, "Type": "...", "Signature": {...}, "Hash": "sha256", "Chirality": +-1,
                 "BendRadiusOverR": r | null, "ExactParameters": true | false,
                 "Length": L, "Portions": [[segment, s0, s1], ...], "Vertices": [v, ...],
                 "PortionTurns": [t, ...] (features with a bend: signed turn toward the metal per portion, radians),
                 "TurnTowardMetal": t (one-sided features with a bend),
                 "Frame": {"Origin": [...], "Axes": [[...],[...],[...]]},
                 "ClaimsFrame": {"Origin": [...], "Axes": [[...],[...],[...]], "Chirality": c} (SpatialEdgeCluster: the claims-only canonical frame),
                 "SpatialSupport": {"Contract": 2 | 3 | 0, "ClaimsKey": "sha256", "ContextDigest": "sha256" | null,
                                    "ClaimsBox": [x0, y0, x1, y1], "Box": [...], "SpanOverR": s,
                                    "SpanCapOverR": 16, "ExceedsSpanCap": b, "LegacyEquivalent": b,
                                    "Growth": {"StepOverR": 0.25, "MaxSteps": 12, "Steps": [n_x0, n_y0, n_x1, n_y1], "AttemptedSteps": [...], "StepReasons": ["..."], "Grown": b},
                                    "FaceRules": {"SnapOverR": 1e-3, "ClearanceOverR": 0.25, "SliversDropped": n, "Crossings": [4 counts],
                                                  "MinClearanceOverR": d | null, "MinCrossingSine": s | null, "MinCrossSectionOverR": w | null,
                                                  "NarrowCrossSections": n, "ExteriorVertices": n, "ExteriorEdges": n,
                                                  "ThresholdBandRelative": 0.01, "ThresholdBandHits": n, "BoxRuleThresholdBandHits": n},
                                    "Context": {"Pieces": n, "ChainPieces": n, "ForeignPieces": n, "ChainLengthOverR": L, "ForeignLengthOverR": L,
                                                "ForeignConductors": n,
                                                "ClaimedByOtherFeature": {"Pieces": n, "LengthOverR": L, "Entries": [{"P": [x0, y0, x1, y1], "LengthOverR": L, "Cluster": c, "Feature": k}, ...]},
                                                "ExcludedSegments": {"Pieces": n, "LengthOverR": L, "Entries": [{"P": [...], "LengthOverR": L, "Class": "Port" | ..., "Reason": "..."}, ...]},
                                                "ChainVertices": [{"P": [x, y], "Type": "...", "FaceDistanceOverR": d, "Site": s, "Cluster": c | null, "Feature": k | null}, ...],
                                                "ForeignVertices": [...]},
                                    "LegacyContinuation": {"StraightContinuationLengthOverR": L, "FictitiousContinuationLengthOverR": L},
                                    "Truncation": {"Segments": n, "LengthOverR": L},
                                    "Unboxable": null | "reason", "UnboxableReason": null | "GrowthStepsExhausted" | "SpanCapRefusedGrowth",
                                    "SpanCapAllowance": {"Label", "SpanCapOverR", "Reason", "Approval", "MatchedQuanta", "DifferingNumbers", "Rule"} (resolved, block (b)),
                                    "ClaimsSignature": {...} (a refused cluster's claims-only signature: the SpanCapAllowances entry's object)}
                                   (SpatialEdgeCluster only; units of R in the feature's frame; the spatial-support contract v3 below),
                 "Match": {"Status": "Matched" | "Missing", "Model": "name", "Deviation": d,
                           "Note": "curvature family: <rule> at kappa k (<convexity>)" | "curvature family: <refusal reason>",
                           "LegacyContract": {"Model", "Key", "ContextDigest", "ClaimsKey", "Reason", "Context": {"Box", "Context"},
                                              "PlacementFrame", "Rule"} (matched through a library alias, USER decision 283),
                           "QuantumNearMatch": {"ModelKey", "FeatureKey", "MaxDeltaQuanta", "DifferingNumbers": {"Count", "Paths"},
                                                "MaxQuanta", "InclusiveMargin", "Rule"} (a SpatialEdgeCluster matched within 4 (+ 1/2) signature quanta, block (b))} } ],
  "Segments":  [ {"Key": [[x0,y0,z0],[x1,y1,z1]], "Length": L, "Chain": c, "Arc": a (chord of Arcs[a]; absent otherwise),
                  "Portions": [[s0, s1, feature], ...] } | {"Key": ..., "Length": L,
                  "Exclusion": {"Class": "...", "Reason": "..."}} ],
  "Arcs":      [ {"Center": [...], "Radius": r, "RadiusOverR": r / R, "TurnDegrees": t,
                  "MaxChordSagittaOverR": s, "Kind": "RoundedCorner" | "Bend", "Joints": n, "Segments": n} ],
  "MeshCoarsenessWarning": {"SagittaOverR": 0.05, "Rule": "...", "Count": n, "Length": L,
                            "WorstMaxChordSagittaOverR": s, "Arcs": [arc indices with s >= 0.05]},
  "Vertices":  [ {"Point": [...], "Type": "Corner|Endpoint|Junction", "TurnDegrees": t,
                  "Feature": k} | {"Point": ..., "Class": "TruncationCut"} ],
  "Exclusions": [ {"Class": "...", "Reason": "...", "Count": n, "Length": L} ],
  "Totals": {"PerimeterLength": L, "AssignedLength": La, "ExcludedLength": Le},
  "Diagnostics": {"SamePriorityClaimOverlaps": 0, "Rule": "...",
                  "StackGeometricOffsetIntervals": 0, "StackCompositionCapHits": 0, "StackCompositionCap": 64,
                  "ClusterExtension": {"Passes": n, "AbsorbedPortions": n, "AbsorbedLength": L,
                                       "VertexFeaturesJoined": n, "Rule": "..."},
                  "SpatialSupport": {"Clusters": n, "ClaimsKeyed": n, "ContextKeyed": n, "Grown": n, "Unboxable": n,
                                     "ExceedingSpanCap": n, "ChainLength": L, "ForeignLength": L,
                                     "FictitiousContinuationLength": L, "ThresholdBandHits": n, "Rule": "..."},
                  "StackEndThirdBodyLength": L, "StackEndThirdBodyRule": "..."},
  "KnifeEdgeCensus": {"BandRelative": 0.01, "SampleSpacingOverR": 0.5, "SampledLength": L,
                      "Distance": {"R": {"Below": L, "Above": L, "Total": L}, "2R": {...}},
                      "BendRadius": {"10R": {...}}, "JointNoiseSagittaOverR": {"0.05": {...}},
                      "ArcMaxJointTurnDegrees": {"50": {...}}, "ArcSagittaOverR": {"0.05": {...}},
                      "ClusterComposition": {"BandRelative": 0.01, "Clusters": n, "ClustersBelow": n,
                                             "ClustersAbove": n, "Below": L, "Above": L, "Total": L,
                                             "Rule": "..."}},
  "UnusedSpanCapAllowances": ["label", ...] (SpanCapAllowances entries no cluster of the run resolved; block (b)),
  "GeometryDigest": "sha256"
}
```

Lengths and coordinates are in mesh units on the `ScaleLength` grid of the version-1 manifest.
The audit reads `Identification` when present (A1 partition from `Segments`, vertex census from
`Vertices`, exclusions from `Exclusions`, A2 from the cluster frames, A3 / A5 from
`GeometryDigest` and the feature signatures) and falls back to the version-1 aggregate reading
otherwise.
