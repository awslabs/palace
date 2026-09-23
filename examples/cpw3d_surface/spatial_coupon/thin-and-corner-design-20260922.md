# Thin spatial coupons and corners under the device library — design draft (decision 65, part 3)

Draft section, not merged into the main doc (spatial_coupon/README / DEVICE-LIBRARY): the scoping the user
asked for before choosing.  Branch simlapointe/device-verification, 2026-09-22.

## (a) THIN spatial coupons under the prism-tube recipe

What the device correction needs: `BuildDomainResponseMatrices(fabricated, thin, ...)` forms the defect
fabricated − thin, so every SpatialEdgeCluster model needs the thin (zero-thickness metal) response on the
same trace basis.  Today `ThinMetal` is an `inputs`-detected scope guard of the prism-tube recipe
(mesh_spatial_coupon.jl 4981 `fabricated || scope_error("ThinMetal", "kind thin")`; mesh_stage_contract
RECIPE_SCOPE_GUARDS) and the qualify writer records the six device models NotLoadable for that reason.

Zero-thickness geometry per edge (thin kind): the metal is a PEC sheet in the process plane z = P_z; there
is no sidewall, no top edge, no trench (overetch does not apply), so the sheet edge is ONE line at z = P_z
where four faces meet: metal sheet behind the edge (MS on its substrate side, MA on its vacuum side —
one physical surface, Palace distinguishes the sides by the interface Type, exactly as the 2D
`thin_metal` tag 2 carries both MS and MA and as the coarse transmon config puts MS and MA on
attribute 5), and the exposed substrate ahead of the edge (SA).  Changes to the recipe, component by component:

| component | fabricated (today) | thin (to do) |
|---|---|---|
| tubes per edge | top tube (z = P_z + t) and bottom tube (z = P_z, trench lip) | one tube at z = P_z; `Layers`/rings/sectors laws unchanged; NormalSize = 0.25 x thickness has no thickness -> take the fabricated value (0.25 nm inner ring: the same MA cutoff as the fabricated coupon, required for the defect to be consistent) |
| tube section | sectors split by the metal top plane (top tube) and by the sidewall / trench floor (bottom tube) | sectors split by the plane z = P_z on the FULL diameter (metal sheet behind, SA floor ahead): the top-tube quarter-plane split reused, the sidewall/trench split dropped |
| trench / overetch | trench volume, SA_side, SA_floor at −overetch, junction curves of the etch cut surfaces | none: SA is the substrate plane outside the plan-view mask; no junction curves (the `JunctionProxyLength` term of estimate_build_cost drops) |
| MA / MS faces | MA = metal top + sidewall, MS = metal foot | the two sides of ONE sheet surface (one physical group, two interface entries) — ownership/relabel must emit one surface label per (conductor, plane) instead of top/foot/side |
| corner balls | one ball per semantic corner and tube cap (two planes) | one ball per semantic corner and cap at z = P_z only |
| census / validator | interface_ownership.jl bands (MS below, MA above the metal slab) | bands collapse to the plane: label by the plan-view mask (inside footprint -> sheet, outside -> SA); audit_rigid_coupon_ownership.jl and the labels-only census need the thin rule; ScopeGuard ThinMetal removed, `ThinMetal` moves to RECIPE_SCOPE_SUPPORTED_CLASSES |
| estimate_build_cost | two tubes + junction proxy per edge | one tube, no junction proxy: roughly half the tube prisms; the trace-basis and far-field terms unchanged |
| manifest | general_mesh_manifest frozen tool digests | one refreeze (mesher, ownership auditor, census) + the calibration manifests |

Effort: **M** (mesher `if fabricated` branches at 5075 / 5154 / 5332 plus the tube-section split, the
ownership labeler, the census/validator rule, estimate_build_cost, tests, one refreeze; no new physics; no
Palace change).  Not S: the tube/ownership code is fabricated-only throughout; not L: the thin section is a
strict simplification of the top tube.

Interim alternative — the legacy graded-tet thin path (`generate_spatial_response.py` spatial_thin +
mesh_spatial_coupon.jl thin kind without tubes) is available today.  Recorded graded_v2 thin costs at
1,536 ranks (8 x 192-rank nodes): 05 1,686 + 2,765 s, 06 2,168 + 3,826 s, 07 1,100 + 1,231 s,
10 1,391 + 1,530 s = 0.65-1.67 wall-h = **5.2-13.3 node-h per thin coupon** (the "0.8-1.7 node-h" of the
decision-64 scout are wall-hours, not node-hours).  Six device coupons -> ~40-80 node-h against the whole
fabricated device library at 3.4-8.3 node-h per run; and a graded-tet thin coupon at 0.5 nm edge size has
a different MA cutoff than the 0.25 nm prism rings, so the defect would mix two cutoffs (the eps^(1/3)
modelled-reference machinery of ma_tail.py would have to be applied to the thin side).  Recommendation:
do the M-sized recipe extension (consistent cutoff, ~0.3-0.8 node-h per thin coupon at the fabricated
p4 cost or less), do NOT run the legacy thin path for the device library.

## (b) CORNERS: CornerCoupon tooling vs a single-corner spatial coupon

Device census: ConcaveCorner 18 + ConvexCorner 12 occurrences, all 90 deg, CornerRadius 0 ("No compatible
sharp-corner model"); two models (convex-90, concave-90) x (thin, fabricated) cover them.

Existing CornerCoupon (`corner_coupon/`): mesher `mesh_corner_coupon.jl` (3D graded tets, Threshold field
SizeMin lc_fine = 0.02 um at DistMin 2 lc_fine, SizeMax 0.3 um at 0.5 R; sharp or filleted corner, angle
option, thin and fabricated kinds), basis = 9 closed square contours x 8 knots = 72 coefficients (its own
contour basis, `generate_corner_response.py`; NOT the device trace basis), finalizer aggregates the tube
energies into the compact matrices, held-out qualification + p2->p3 + h convergence in the planner
(`prepare_surface_response_coupons.py build_corner`).  Library: Version 3 with Angle/CornerRadius
tolerances, ConductorReferences, ZeroTraceIndices — `combine_process_libraries.py` accepts it next to
SpatialEdgeCluster models and Palace matches corners by topology/angle/radius (`ReadProcessLibrary`
corner branch 1690-1702), so it is compatible with the qualify writer's Version-3 file (a `--merge-into`
of a combined file or a second merge step).  Cost (README, Apple M3 Pro): 6-probe totals 30 s / 3-4 min /
12-14 min at p1 / p2 / p3; the 72-knot full response ~12x -> ~40 min (p2) to ~3 h (p3) per coupon kind
locally, 2 models x 2 kinds -> 3-12 h local or one small cluster job.  MA: lc_fine 20 nm at the edge is
80x the prism ring (0.25 nm); by the eps^(1/3) tail law the MA cutoff deficit is ~4.3x the fabricated
coupons' 1.9-2.7% -> ~8-12% MA deficit (README: MA +10% p1->p2, +2.5% p2->p3, not converged).  Either
reduce `--corner-lc-fine` (cost grows ~h^-2 in the edge band) or record MA_raw with the modelled deficit.

Single-corner spatial coupon (prism-tube recipe): a 2-edge SpatialEdgeCluster whose edges meet at one
semantic corner.  The mesher/basis producer supports it today — semantic corners with corner balls are
part of every device cluster (10-edge: 10 corners; 3-edge 5d3b: 5) and the 2-edge class is built
(spatial-2-edge-8ce3fb213a2d: 1,215,960 elements, 1,118 s local build, p4 0.31 node-h) — so the coupon
cost is ~1.2M elements / ~0.3 node-h fabricated + the thin counterpart once (a) exists, with the device
trace basis and the 0.25 nm ring cutoff (consistent MA).  What does NOT exist: Palace's requirement
classifier emits corners as ConvexCorner/ConcaveCorner requirements (CornerCoupon route), and
`FindSpatialClusterLibraryModel` matches SpatialEdgeCluster models only to sites the classifier grouped
as a spatial cluster; a corner-only spatial model would need the classifier to emit the two arms of an
isolated corner as a 2-edge cluster requirement (with the arm lengths inside R) — a Palace change
(classification + the plan-view boundary export for the arms), M/L, plus the discovery census and
`device_coupons.py` routing.  estimate_build_cost needs a frozen case to run; the 2-edge device case is
the measured proxy above.

Recommendation: **CornerCoupon now** (existing, Version-3 compatible, no Palace change, 2 models cover
all 30 corners), with the MA cutoff caveat recorded (run `--corner-lc-fine` 0.02 -> 0.005 as the h
convergence and record the eps^(1/3) modelled deficit next to MA_raw); the single-corner spatial coupon
is the consistent-cutoff successor once (a) is in and the classifier change is approved.

## Writer workaround to remove later

Palace reads `ThinMatrix` unconditionally as a string, also in `--surface-response-preflight`
(surfaceresponseoperator.cpp:1756); the qualify writer therefore emits `process-library-preflight.json`
(geometry-only; NotLoadable models point at the declared thin copy paths) next to the honest
`process-library.json` (ThinMatrix null).  A Palace change accepting null in geometry-only mode would make
the second file unnecessary — noted, not done.
