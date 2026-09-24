# Partition-independent metal-edge classification (decision 69, 2026-09-23)

Commits 84bea9490 (Fix 1: canonical perimeter vertex coordinates — the nodal grid function instead of the
partition-averaged `GetVertex`, distinct crack copies averaged in sorted order), bdfd57542 (Fix 2: canonical-order
reductions of the process normals / gap directions + the regression test "Automatic metal edge classification is
partition independent"), cdfe5d4be (Fix 3: gathered plan-view facets sorted), 5876402f7 (Fix 4: quantized decisions —
distances compared on the 1e-8 R grid, direction cosines on the 1e-12 grid, recorded under Library.DecisionQuantization).
Linux executable palace-archive-estimate-9ef5256bc9ca9abc109954fb97b096879b3b472f18483f867b8de0b73a832058.bin (PBS 48697,
source freeze 222a3d33... of 5876402f7) is the qualify default (build_plan.DEFAULT_FROZEN_BINARY_SHA256); 170439c4...
stays recorded as PREVIOUS_ / COST_MODEL_FROZEN_BINARY_SHA256 (the solver stages are unchanged).
Evidence: /Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/classification-fix-20260923/ (REPORT.md).

## Acceptance
- Local (1, 2, 3, 4, 6 ranks) and cluster (192 / 384 / 768 ranks, PBS 48700 / 48701 / 48702) geometry-only preflight of
  the island-tagged coarse transmon (mesh c4fe532f...) against the Version-3 device library (42b19686...):
  surface-response-requirements.json byte-identical apart from Statistics / Path (sha256 6aa56d6c905e5e70): 3,007 isolated
  segments, 34 corners, 2 spatial patches; the 384-rank abort (metaledge.cpp:557) is gone.
- AMR p4 (10 cycles, SurfaceMortar 2) at 768 ranks (PBS 48702): the recorded 48257 element path and raw values reproduced
  (8,501,504 elements, 64.6M DOFs; raw SA -0.35 / MS +10.65 / MA -49.05 / C +0.46 %). At 192 ranks (48700, 1 node: cycles
  1-10, final solve OOM as the recorded 48250; 48872 on 2 nodes, see summary) the classification is identical and the
  corrected values agree to 3.8e-7 relative on the identical meshes of cycles 1-2; the AMR loop's own partition dependence
  (Dorfler marking on solver-tolerance indicators, the same two paths as the recorded 192 / 768 runs) separates the meshes
  from cycle 3 (<= 2.8e-3 relative) — expected, record only.
- Corrected values, 768 ranks, final: SA -4.74 / MS -3.25 / MA -2.22 / C +0.15 % vs r5nm-t5um p5 (recorded 170439c4:
  -2.97 / -0.80 / +0.28 / +0.14). The gate "equal to the recorded values" is NOT MET BY DESIGN: 318 um of the 35.3 mm
  perimeter (4 spatial-cluster sites + 64 two-micron same-conductor strip segments) have no model in a library whose
  cluster models were discovered with the old, partition-dependent classification (below). Decision 70 (deferred): the
  library and the device run must use the same classification executable (re-discovery + rebuild).

## The four changed sites (event dumps of both builds at 1 rank, PALACE_DEBUG_SPATIAL_EVENTS, not committed)
The spatial-cluster machinery keeps ONE interaction event per pair of perimeter chains — the closest segment pair, the
first found in candidate order — and the site point of every chain in a cluster is the mean of its events; events merge
into a cluster when their centers lie within 16 um (4 x 2R, exhaustive spatial closure). For parallel edges at a uniform
distance every overlapping segment pair ties, so the representative was the first pair in Allgather rank order in the old
build and is the first pair in the canonical (sorted-perimeter) order now. Both builds find the same 44 cross-interface
chain pairs and the same same-signature chain pairs away from the junction; only the tied representatives differ.
The relevant layout (plan-view, um): a 1-um JJ lead per conductor (x = +-0.5) meeting across a 1-um gap at the junction;
two 2-um-wide vertical strips at x = -44..-42 and 42..44 rising from a pad top edge at y = 535 to a 2-um horizontal bar
(edges y = 662 / 664, x = -44..44) — a rectangular frame with four 2-um corners.
1. 10-edge cross-interface junction cluster (40 um): model 6791... has the c2 lead-edge site at 1.8333, the new site is
   1.5714. The chain pair (30, 34) = the two long edges of the c1 lead (1 um apart, y 6.5-12) is represented by the segment
   pair y 6.5-8.33 / 6.5-8.33 (canonical) instead of 6.5-8.33 / 8.33-10.17 (rank order); the event point moves along the
   lead and the site mean moves from -11/6 to -11/7 um. Both are representatives of the same uniform 1-um gap: the old
   value was an ordering artefact, the new one is a function of the mesh. The model encodes the artefact -> Missing.
2. Right strip base, 3-edge cluster (12 um) Missing; the old 7f03... (3 edges at z 0 / 7.82 / 15.64) was an artefact:
   the strip's chain pair (39, 40) was represented by the segment pair y 546.9-550.9 / 550.6-554.5 — 15.6 um UP the
   strip, rank order — so the "3-edge cluster" was the strip base (event (39, 41) at y 535) plus the strip body 15 um
   away; the new build represents (39, 40) at the base (y 531-535, canonical) and reads a clean strip-base corner: two
   parallel edges at separation 2 plus the perpendicular pad edge, the mirror image of the left base (model 4195...,
   Exact). The library has only the left handedness -> Missing.
3. Horizontal bar: the old 2-edge model 3f89... (two parallel edges, 4-um window, nothing else) was the bar's chain pair
   (21, 24) represented at x = -20..-16 / -18..-14 (rank order), > 16 um from both corner event groups (x = +-42), so it
   formed a standalone "cluster" whose exhaustive closure then covered the bar body. The new build represents (21, 24)
   at x = -44..-40 / -42..-38 (the bar's left end, canonical); the event joins the top-left corner cluster (the 4-edge
   model 9d2c..., Exact) and the bar body becomes what it is geometrically: a SameConductorStrip family at separation
   2.0 um (64 segments, 250.18 um) — Missing because the library has only same-conductor-strip-1um.
4. 4-edge cross-conductor cluster (16 um, two split edges 12 um apart) Missing, in place of the old 8-edge (32 um) and
   2-edge cross-conductor (8 um) Missing junction descriptions: the 14 junction segments within 2R of an edge with a
   different interface mapping are omitted from the same-signature path in the new build (log "Omitting 2 of 2" / "12 of
   3086") because their cross-interface cluster is Missing; in the old build the matched 10-edge model kept them in that
   path, which described the junction a second time (11 extra junction events in the old per-group dump, none in the new).
   Which chains form the new 16-um cluster is UNSETTLED from the data in hand (it is not the junction and not a frame
   corner); every non-junction per-group event is identical in both builds except the two tied representatives above.
Tie rule: replacing every quantized "<" by "<=" gives a byte-identical classification — no decision sits on the quantized
threshold in this device; the changes come from canonical ordering (Fixes 1-3), not from the tie convention (Fix 4).

## 192-rank AMR p4 on 2 nodes (PBS 48872, --map-by ppr:96:node, 15:38)
Unused for the acceptance (geometry classification is accepted on the preflight identity; AMR solves are not part of it).
For the record: classification 3,007 / 34 / 2, preflight sha 6aa56d6c905e5e70 (identical), final 63,741,334 DOFs (the
recorded 48256 path), SA -4.72 / MS -3.24 / MA -2.20 / C +0.15 % (768 ranks: -4.74 / -3.25 / -2.22 / +0.15).
