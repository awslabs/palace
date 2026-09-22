<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is a fail-closed inventory and gate contract
for the mesh-only phase of the coupon-library acceleration plan. It does **not**
qualify response fields, trace accuracy, a production library, or all release
gates.

## Frozen inputs and required matrix

`geometry-independence-suite.json` is the only production manifest. It freezes the
Gmsh-only pipeline (`Pipeline: gmsh-only`, decision 38, below); the labeled
calibration manifest `geometry-independence-calibration-ma.json` freezes the legacy
MMG pipeline (`Pipeline: legacy-mmg`) so its calibration cases and evidence remain
verifiable, and the labeled calibration manifest
`geometry-independence-calibration-sizing.json` (supervisor decision 41) freezes the
Gmsh-only pipeline for the sizing calibration of the four-edge case (below).

### Resource bounds are guidelines (supervisor decision 61c / user decision 60(3), 2026-09-21)

The three build-machine bounds of the manifests' `Gates` are guidelines for the 36 GiB
workstation (at most two concurrent builds), not restrictions on the physics:
`MaximumElements` 4,000,000 -> **6,000,000**, `MaximumSeconds` 1800 -> **3600** per
bounded stage, `MaximumRSSGiB` 8 -> **12** (two 12 GiB builds and the operating system
fit 36 GiB). The production manifest carries the new values, mirrored into both
calibration manifests; every manifest records the previous values, the decision and
the reason under `Gates.PreviousValues` (the `Gates` block is the canonical build cache
key, so a root built under the previous bounds carries its own recorded `Gates`).
Reason: the transmon device library (decision 58) built its largest coupon,
`spatial-3-edge-5d3b5e644745`, at 3,037,912 elements with a verification stage at
1692 s (0.94 of the 1800 s bound, the run's one headroom flag) - the device coupons
approach bounds calibrated on the gallery cases. Nothing else moves: every physical
gate (positive orientation, SJ 0.01, condition 1000, corner aspect 4.0, protected
1e-8, closure 1e-12) is unchanged; the estimate gate
(`estimate_build_cost.preflight_build_cost` against `Gates.MaximumElements`), the
mesher's `--max-elements` / `--max-nodes` budget (a fail-closed gate on the generated
mesh: `FarFieldBudgetPolicy` Pressure 1, the size field never depends on the cap), the
bounded-stage limits (`run_bounded_mesher.py` at `Gates.MaximumSeconds` /
`Gates.MaximumRSSGiB`, audits included) and the 0.9 headroom flags
(`run_gmsh_only_matrix.HEADROOM_FRACTION`) keep failing closed at the NEW values
(`test_general_mesh_manifest.test_resource_bounds_are_guidelines_raised_by_decision_61c`).
The decision-33 calibration-only element cap (`Calibration.MaximumElements` 5,000,000
of `four-edge-calib-ma-el1c`, recorded build 4,483,816 tets) is retired
(`Calibration.RetiredGateDeviations.MaximumElements` of the MA calibration manifest;
the case is judged by the 6,000,000 gate like every other case; the per-case
`UnchangedProductionParameters.seed["--max-elements"] 4000000` texts of the
calibration cases record what those builds executed and are unchanged). Evidence
recorded before 2026-09-21 quotes the bounds of its time (4M / 1800 s / 8 GiB).

### Production recipe (supervisor decision 42, 2026-09-19): Gmsh-only, TraceBasisSizeRatio 0.5

The production recipe is the Gmsh-only build of decisions 38-42
(`ProductionRecipe` of `geometry-independence-suite.json`): prism edge tubes on every
straight metal edge (inner ring 0.25 nm, ratio 2, layers following the size field on
the axis with the surface layer at on-box ends, explicit pyramids), isotropic corner
balls graded to the tube inner size, the etch footprint, junction-line and
footprint bands, the decision-39 volume laws growing with FarGrowth 0.5, the
trace-basis cut-surface / volume law at **TraceBasisSizeRatio 0.5** (decision 42:
size <= 0.5 x the minimum altitude of the nearest bound-basis triangle + FarGrowth x
the distance to it; dimensionless, contract-derived, 1.0 before; the measure was the
shortest edge until decision 43, below) and the far field by size fields. `ProductionRecipe.BuildCommandOptions` = `--lc-tangent 0.05
--edge-size 0.00025 --edge-growth-ratio 2.0 --corner-size 0.00025 --far-growth 0.5
--trace-basis-size-ratio 0.5`; `validate_production_recipe_commands` requires the
recorded `gmsh-build` command of every production case to execute each option
exactly once at its value - the ratio exactly once at 0.5 when the case binds a
trace basis and never otherwise -, `run_gmsh_only_case.py` takes the ratio from the
manifest options alone (no default: a basis-binding case under a manifest without
the option fails closed), and `mesh_stage_contract.validate_trace_basis_sizing`
binds the census `TraceBasisSizing.Ratio` to the executed value. Every physical
gate is unchanged (positive orientation and Jacobian condition <= 1000 for every
element type, tetrahedra SJ >= 0.01, corner aspect 4.0, protected 1e-8, closure
1e-12; the resource guidelines 4M elements, 1800 s / 8 GiB of the time - 6M / 3600 s /
12 GiB since decision 61c). **The MMG path (metric preparation, native
MMG adaptation, label restoration, required-region optimizer, tetrahedral edge
layer) is legacy calibration only**: it survives under the labeled
`geometry-independence-calibration-ma.json` (`Pipeline: legacy-mmg`) so its study
records stay verifiable, and it is never executed by a production case.

Physics evidence of the recipe on the four-edge coupon (graded_v2 reference, p4, 80
sources / 60 free-view; the assessment directory
`coupon-accuracy-assessment-20260913`):

| run | mesh (identity SHA-256, elements, H1 p4) | recipe state | E at 74 / 10 | E within 1% (worst) | p_SA 48 / 47 / 35 | p_SA within 1 / 2 / 5% | p_MS median; worst | p_MA median / weighted / strongest-20 (within 2%) | node-h per coupon |
|---|---|---|---|---|---|---|---|---|---|
| physics-08 (PBS 44988, decision 37 spike) | 7c17487f..., 1,588,715, 21.5M | prism tubes, ratio 1.0, no volume grading | +10.2 / +9.4% | 48 | -13.1 / +6.9 / +8.3% | - | -1.06%; +/-6-10% at 8 sources | +0.47 / +0.52 / +0.76% (49) | 0.55 |
| physics-09 (PBS 45076, decision 38 DAG) | e57ce087..., 1,680,422, 22.2M | production DAG, 50 nm junction layer | +10.3 / +9.8% | 48 | -13.1 / +7 / +7% | - | +/-6-10% | +0.58 / +0.70 / +0.86% (54) | 0.516 |
| physics-10 (PBS 45085, decisions 39-40) | 647b2079..., 1,742,434, 23.0M | volume laws, variable tube layers, ratio 1.0 | +6.5 / +7.3% | 54 (+7.3%) | -6.1 / +2.8 / +2.1% | 29 / 39 / 54 | -0.69%; +7.2% (30) | +0.67 / +0.72 / +0.80% (56) | 0.520 |
| **physics-11 (PBS 45145, V-a = decision 42 production)** | **80966c7d..., 1,869,209, 24.8M** | **ratio 0.5, surface layer** | **+0.03 / +0.15%** | **60 (+0.43%)** | **-5.2 / +1.8 / +1.9%** | **35 / 44 / 58** | **-0.60%; -1.70% (47)** | **+0.76 / +0.76 / +0.76% (55; 20/20 strongest)** | **0.508** |
| physics-12 (PBS 45146, V-b, rejected) | 7db360be..., 2,119,019, 31.4M | TangentialSize 25 nm, ratio 1.0 | +7.1 / +7.2% | 54 (+7.2%) | -6.1 / +2.9 / +2.1% | 29 / 39 / 54 | -0.75%; +8.1% (30) | +0.61 / +0.71 / +0.82% (55) | 0.593 |
| EL4c (physics-06, PBS 44717; retired MMG production) | 3,471,507 tets | legacy MMG recipe | +0.5 / -0.2% | 60 (+0.81%) | -8.6 / +0.4 / +0.7% | 35 / 40 / 56 | -0.75%; -9.3% (26) | -6.63 / -5.00 / -5.05% (6) | 0.752 |

V-a reaches the EL4c level or better on E, SA and MS while MA stays converged, at
0.98x the physics-10 cost and 0.68x EL4c's; V-b moves nothing (<= 0.53 points) at
1.14x cost. Decision 42 adopts the V-a ratio: the production four-edge root rebuilt
under `geometry-independence-suite.json` reproduces the V-a identity mesh
(`80966c7db44dabc49ac7bb068bbaee0e0828e413b115cee8c066c886866b6108`), so physics-11
binds to the production root (evidence below). Residuals recorded by physics-11: MA
movers 47 -4.2% / 74 +3.2% and the alternating MA at 35; p_MS at 26 alternates with
p on every mesh; the reference is a p4 anchor.

### Trace rule measure and composed curve spacing (supervisor decision 43, 2026-09-19)

Physics on the real gallery case `three-edge-419576fdab24` (input 06; evidence
below) transferred the recipe at EL4c level or better on every class except four
isolated 0.05 um slivers (sources 25 / 26 / 133 / 134 at (-8, -0.5) / (-8, -0.55) on
the bottom / top box edges: E +4..+8%, MA / MS / SA -3..-10%), while the same sliver
at z = -0.05 (52 / 53) was accurate. Class label used here and in the physics
RESULTS (gallery-physics-06 / -06b): the **12-source narrow-hat class** of case 06 =
hats narrower than 0.1 um, split into the **8 junction-adjacent slivers** (14 / 16 /
41 / 43 / 52 / 53 / 122 / 124, 33-90 nm wide, at most 0.09 um from a junction, covered
by the junction / band rules) and the **4 isolated box-edge slivers** (25 / 26 / 133 /
134, 50 nm wide, 0.5 um from any junction, outside every tube / junction / corner
rule). Only the isolated four failed under decision 42; the whole class transfers
under decision 44. The diagnosis on the verified root refuted the
first hypothesis (box-edge curve spacing): the outer-box edges are not explicitly
1D-meshed (they are excluded from the feature curves and Gmsh meshes them from the
composed callback field), and at the four supports the box-edge nodes were at
24.8 nm against 24.9 nm prescribed (achieved / prescribed P50 0.99 within 0.3 um),
the cut triangles at mean-edge / prescribed 0.90 on both faces (0.96 at 52 / 53)
and the tetrahedra within 0.1 um at longest-edge / prescribed 1.49 (1.45 at 52 /
53) - the realization was the same at the bad and the good supports. The cause was
the prescription's measure: the four sources are the endpoints of the 49.8 nm
shortest edge of two needle basis triangles on the bottom / top faces
([(-8, -0.55), (-8, -0.5), (-6, 8)]: edges 0.0498 / 8.73 / 8.78 um) whose minimum
altitude is 11.3 nm, so their hats vary at 1 / 11.3 nm across an 8.7 um needle
while the rule "ratio x shortest edge" prescribed 24.9 nm = 2.2x that scale. Exactly
the four free sources whose own hat altitude is below 0.6 x their shortest adjacent
edge were the four outliers (every other source has an altitude >= 0.11 um); the
four-edge basis has no such needle (its 22 nm slivers 74 / 10 are right triangles,
altitude = shortest edge, and are accurate).

Decision 43 (no new parameter; the recorded rules change):

- **Trace rule measure = the basis triangle's minimum altitude** (2 x area /
  longest edge = 1 / the largest gradient of its three vertex hats), ratio 0.5
  unchanged: the quantity the mesh must resolve is the Dirichlet datum's variation,
  and the shortest edge is a proxy that coincides with the altitude for right
  slivers (where 0.5 was calibrated) and overstates it for needles. The census
  `TraceBasisSizing` records `SizeMeasure`, `MinimumBasisAltitude`,
  `MinimumRequestedSize` = Ratio x it, and the report-only needle counts
  `NeedleTriangles` / `NeedleTrianglesBelowFarSize` (altitude < 0.6 x shortest edge,
  `NeedleRule`; the threshold `trace_basis.NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE` is
  defined once in Python - the stage contract and the fixture producer import it, the
  Julia census constant mirrors it and the contract requires the recorded value to
  equal it - and is a report-only classification threshold that never enters a size:
  altitude / shortest edge = (b / c) sin C is 1 for the right slivers the ratio was
  calibrated on, 0.866 equilateral, 0.707 right isosceles, so 0.6 lies below every
  well-shaped triangle and counts exactly the triangles whose shortest-edge proxy
  overstated the hat scale by more than 1.67x - case 06's needles are at 0.23);
  `mesh_stage_contract.validate_trace_basis_size_measure` (called by
  `validate_gmsh_build_census`) recomputes all of them from the recorded mesh-frame
  triangles. The Python metric rule (`trace_basis.requested_sizes`,
  `basis_statistics`, `cut_surface_size_report` - now measured across the minimum
  altitude, `MaximumExtentAcrossMinimumAltitude`) and the Julia seed rule are the
  same measure (`test_trace_basis.py`). Needle basis triangles are a property of
  the reference basis triangulation, not of the coupon mesh: avoiding them in the
  library's own basis construction is a future design item, not part of this
  change. Basis statistics under the altitude measure: four-edge minimum altitude
  15.3 nm (shortest edge 21.7 nm; 4 needles, none below the far size), ten-edge
  5.5 nm (49.2 nm; 50 needles, 34 below the far size, the worst a 19.4 um needle),
  three-edge 06 11.3 nm (32.7 nm; 18 / 6), two-edge 10 48.2 nm (50.0 nm; 12 / 2).
- **Every explicitly 1D-meshed curve follows the composed size field** (the
  decision-40 principle applied to the curves, `CURVE_SPACING_RULE`): the band
  curves (NormalSize) and the un-tubed metal ridge parts (TangentialSize) take
  min(their explicit spacing, the composed field on the curve: corner law of the
  graded points, trace rule, band rule, corner exterior rule), sampled along the
  curve and gradient-limited to (GrowthRatio - 1) / GrowthRatio exactly as the tube
  axes are (`graded_tube_stations`); the spacing grid is kept on every grid interval
  where the limited law equals the spacing (ridge alignment across faces), the
  remaining gaps are equidistributed in the arclength integral of the reciprocal
  law with ceil(integral) intervals (every node interval <= the size it spans). The
  census `CurveSpacing` records the rule, the growth cap and one row per curve
  (kind junction / band / metal, segment, length, spacing, graded flag, grid
  intervals kept, interior nodes, node spacing min / P50 / max, prescribed minimum,
  achieved-over-prescribed min / P50 / max); `validate_curve_spacing` binds the rows
  (spacing = NormalSize on band / junction rows and TangentialSize on metal rows,
  statistics ordered and within the spacing, achieved-over-prescribed within the
  growth cap - an interval is bounded by the law over its span, which steps down by
  at most GrowthRatio inside it at a corner-ball shell -, kept grid within the grid,
  the band curves covered). A law marginally below the spacing on one grid interval
  splits it in two (ceil), so a band curve can show a 12.5 nm minimum interval. Before this the
  band curves ignored the trace rule along their length (a basis sliver crossing a
  junction line saw 25 nm curve nodes at a 11 nm request; the surface mesher cannot
  refine a curve's nodes), the most plausible remaining cause of the junction-source
  SA residual of physics-11 (47 / 35 / 31 at +1.8 / +1.9 / +0.8%, 0.6-1.4 points
  above EL4c) - recorded as such, to be adjudicated by the next physics run.
  `corner_isotropic_curve_nodes` (the corner-law-only placement of the prism-tube
  spike and the face-census tests) is now the composed rule without the trace and
  band laws.

Cost of the altitude measure (probe builds before adoption, identical otherwise):
case 06 1,693,002 -> 1,843,365 elements (+8.9%; the bottom-face needle at source 25
from 116 triangles at 36 nm mean edge within 0.3 um to 868 at 6.9 nm, the
tetrahedra within 0.1 um from 136 to 940); ten-edge 2,673,691 -> 3,495,791 (+30.7%,
87% of the unchanged 4M cap, accepted by the supervisor; the cap fails closed).
The production roots are rebuilt below.

### Gmsh-only sizing calibration manifest (supervisor decision 41, 2026-09-19)

`geometry-independence-calibration-sizing.json` is a CALIBRATION-ONLY manifest of
the Gmsh-only pipeline: it carries a `Calibration` block (no `ProductionRecipe`),
`GateDeviations` naming the one labeled per-case deviation (decision 53 below; every
other gate at its production value), and `Tools` /
`StageToolSHA256` mirrored from production by `refreeze_manifest_tools.py`
(`CALIBRATION_MANIFESTS`; `--check` reports a stale mirror of either calibration
manifest). Its cases clone the production four-edge case (same immutable inputs,
hashes, semantic contract, variants and covariance comparison; `InventoryStatus`
Calibration, `Calibration.BaseCase`) and declare `Calibration.BuildCommandOptions` -
the build options that differ from the baseline, any subset - against
`Calibration.ProductionValues`, **the production `BuildCommandOptions` at the time
of the study (before decision 42: `--trace-basis-size-ratio 1.0`), a historical
baseline exactly like the legacy manifest's `ProductionValuesBefore34B`; it is not
updated after an adoption** (the production manifest now carries
`--trace-basis-size-ratio 0.5`): `four-edge-calib-sizing-tbr-0.5` (V-a,
`--trace-basis-size-ratio 0.5`, labeled `AdoptedAsProductionRecipe` by decision 42 -
the calibration entry proving the adoption, the production four-edge root reproduces
its identity mesh SHA-256 `80966c7d...`) and `four-edge-calib-sizing-lct-0.025` (V-b,
`--lc-tangent 0.025`, rejected by physics-12). `validate_manifest`
(`validate_calibration_case_options`) requires the label, non-empty options and
finite baseline values differing from the declared ones; `verify_canonical_case_entries.
validate_calibration_commands` is pipeline-aware (`PIPELINE_CALIBRATION_STAGE_OPTIONS`:
the legacy seed / metric / adaptation / restoration blocks against
`ProductionValuesBefore34B`, the Gmsh-only `gmsh-build` block against
`ProductionValues`) and requires the recorded build command to execute every declared
option exactly once at its value, never at its baseline value, and every undeclared
baseline option only at its baseline value - the canonical cache key does not
encode recipe options, so the recorded command is the binding of the label to the
build (a root labeled V-b but built with the baseline options, or with V-a's, is
rejected). `run_gmsh_only_case.py CASE --manifest geometry-independence-calibration-sizing.json`
builds a calibration case from `ProductionValues` overridden by `BuildCommandOptions`
(`case_build_options`; `--trace-basis-size-ratio` is one of the options, taken from
the options alone and passed with the bound trace basis; a production case executes
the production recipe's 0.5) and writes `CALIBRATION.txt` instead of
`PRODUCTION.txt` in the root.

**Decision 53 (2026-09-20): `two-edge-calib-tube-rings-0.125nm`.** The third case of
the sizing manifest clones the production two-edge-8dd4bc70f183 case and declares
`--edge-size 0.000125` / `--corner-size 0.000125` against the CURRENT production values
(decision 42, `--trace-basis-size-ratio 0.5`; `ProductionValuesRule`): the innermost
tube ring halved at the unchanged growth 2, so the ring rule yields one more ring (8
rings to 31.875 nm; production 7 to 31.75 nm) and the corner balls grade from the same
size - the h-refinement lever for the two-edge p_MA strongest-20 excess at 53 / 58.
Halving the inner ring at the production tube layer spacing doubles the innermost
prisms' Jacobian condition by construction (ring / layer aspect): the production
two-edge worst prism is 586.30, the refined ring set's 1172.6 = 2 x 586.30, above the
production bound 1000, so the case declares `Calibration.MaximumJacobianCondition 1200`
labeled under `Calibration.GateDeviations.MaximumJacobianCondition` exactly on the
decision-33 element-cap precedent (`Production` 1000 recorded, `Cases` naming it,
`ProductionUse FORBIDDEN`; `general_mesh_manifest.validate_case_jacobian_condition`,
applied by `case_gates` in the mesher command of `run_gmsh_only_case.py` - the manifest
`Gates` stay the canonical cache key, the judged gates are written to `case-gates.json`
- and in the verification). The pre-build estimate gate now judges calibration cases
too (`estimate_build_cost.build_options_and_model`: the case's own labeled options with
the production manifest's model), and `run_gmsh_only_matrix.py` accepts a labeled
Gmsh-only calibration manifest (records `Calibration` per case and
`Library.Manifest.Kind calibration`). Build (2026-09-20): estimate 574,787 (0.144 of
the cap; production two-edge 560,105), actual 537,069 elements (427,629 tets + 102,600
prisms + 6,840 pyramids; production 521,676 = 425,916 + 88,920 + 6,840: +15,393 = +2.95%,
the prisms +15.4% = the eighth ring), H1 p4 8,445,107 (production 7,973,827, +5.9%), 760
tube layers on 35.17 um of tube (unchanged), spacing 15.8-49.7 nm (unchanged), maximum
prism edge aspect 759.7 (production 379.9), prism Jacobian condition max 1172.56 (11,772
cells above 1000, 0 nonpositive, min scaled Jacobian 0.0305), tetrahedra max condition
124.4 / min scaled Jacobian 0.0200, pyramids 8.69 / 0.292, cap regions 124.4 / 0.0200 -
every gate but the labeled prism bound at its production value; gmsh-build 41 s /
2.45 GiB; identity + rotate-z verified (Passed). The record for `qualify` is the
calibration root's `library-build.json` (`Library.Manifest.Kind calibration`;
`qualify` reads the run parameters from `Calibration.ProductionManifest`); the physics
run of this case must report PCG iteration counts and kappa next to the production
two-edge run (supervisor: a material PCG degradation bounds the ring lever by
conditioning - itself a finding).

**Outcome of decision 53 (2026-09-21, PBS 46094, `qualify` at f46d34324 unchanged; record
`qualify/he-check-20260921/HE-CHECK.md`, root `/tmp/library-he-check-01`).** The ring-refined
variant through the frozen gates (p4 + p5 main, p3 / p5 controls including 53 / 58, gated
p5 vs the case-10 p5 anchor): E 78/78 within 1%, p_MS 77/78/78, p-sequence controls
passed, verdict Failed on p_MA - median +1.72%, weighted +1.65%, strongest-20 16/20 with
53 +4.49%, 58 +4.47%, 69 +2.03%, 71 +2.41% (production p5: +1.25 / +1.12%, 18/20, 53 +3.50%,
58 +3.46%). **Our MA is NOT stable under the ring**: at equal p5 the ring set moves p_MA by
+0.96 / +0.97% at 53 / 58, +0.47% at the strongest-20 median and +0.54% at the all-78
median, upward at every source, while E and p_MS move < 0.1%; the p-step p4 -> p5 is +0.7-0.9%
at 53 / 58 on both meshes (ring p3 / p4 / p5 at 53: +2.99 / +3.78 / +4.49%, Aitken r 0.91).
The decision-53 "moves" branch applies (ring-set sensitivity about -1.0% p_MA per halving
of the inner ring at 53 / 58, -0.5% at the median source), with the qualification that every
refinement moves p_MA AWAY from the p5 anchor: the anchor is at least as under-resolved in
the same direction, so the ring lever alone cannot close the 2% statement against this
reference and a "reference-limited" annotation is not supported by this evidence either.
Conditioning: Palace kappa max 831 vs 416, PCG mean 21.6 / 20.4 (p4 / p5) vs 11.9 / 10.8
(Tol 1e-10 vs the acceptance's 1e-8 explains about 1.25x; the ring about 1.5x at equal Tol),
seconds per iteration unchanged; 0.968 node-h (1.67x the production run at equal p).
Remote archives (16.4 GB) deleted after digest verification, 22 MB left. No code change.

**Localization and the sharp-edge tail (decision 55, same run).** The per-segment MA table of
the p4 local-edge stage (controls incl. 53 / 58 on the ring mesh; the production run has no
53 / 58 there and no run has a p5 local-edge stage): 53 = the hat on the cut face x = 4 at
(4, -0.5, 0.1), 0.5 um from the metal-edge / cut junction (4, 0) along conductor 2's y = 0 edge,
4.03 um from the nearest 90-degree corner (0, 0); 58 its mirror at (-5, 1.5, 0.1). Their MA is
99.2% on that single edge line - 80% straight edge, 20% junction ball (r < 0.1 um), 0.00% corner
balls; 52% within 0.25 um of the junction. The facing conductor ends (x = -1 / x = 0) are 1.0 um
apart vs 2 x (R_tube 0.032 + pyramid 0.008 + band 0.05-0.1) = 0.2-0.3 um: no tube / band
interaction. Ring - production per shell at the 8 common controls (p4): total +0.60% = straight
edge +0.51%, corner balls +0.08%, junction balls +0.03%, uniform along the edges at +0.5-0.7% of
each 0.25-um bin - the +1% ring step and the +0.8% p-step at 53 / 58 live along the straight
edge next to the junction: the lever is the tube ring set / tangential spacing, not the
corner-ball grading. MODEL (labeled): for a sharp 90-degree edge E ~ r^(-1/3), the MA integrand
~ r^(-2/3), the edge MA inside eps ~ eps^(1/3), so each ring halving recovers 1 - 2^(-1/3) = 0.206
of the remaining edge deficit (increment ratio 0.79). The two increments we have - two-edge
0.25 -> 0.125 nm +0.5-0.6% along the edges, four-edge 0.5 nm-tet reference -> 0.25 nm ring
+0.62 / +0.79% at the edge-interior controls 80 / 23 (median +0.82%) - have the ratio 0.55 / 0.70 =
0.79, consistent within the noise of two coupons; the implied remaining edge-MA deficit (hedged
extrapolation) is 3.4 / 2.7 / 2.1% at 0.5 / 0.25 / 0.125 nm on both chains (53 / 58: 4.7 / 3.7% at
0.25 / 0.125 nm). Conclusion: the strongest-20 2% statement vs a 0.5 nm-anchored reference is not
decidable by refinement - both sequences share the same slowly converging edge tail; the
resolutions (analytic shell extrapolation of the edge remainder from per-ring MA labels, or a
rounded-edge process model with its own references) are physics / definition decisions for the
supervisor / user, not mesh levers (`qualify/he-check-20260921/HE-CHECK.md`).

### Radial MA shells: measuring the edge tail directly (supervisor decision 56, 2026-09-21)

Before the MA-definition decision, the tail is measured rather than modelled: the
PRODUCTION two-edge identity mesh (SHA256 `5d01204e...`, the acceptance-run mesh of PBS 46023)
is relabeled, label-only, into per-tube-ring radial shells of its MA surfaces and run through
`qualify`, so the MA of every source is recorded per ring. Tool `relabel_radial_ma_shells.py`
(byte-level on the MSH 2.2 binary; nothing is built): a surface element of an MA label
(6000 + 100 slot + conductor: metal top face and sidewalls) goes to shell k when the distance
of its centroid to the nearest metal edge line - the Physical plan-view boundary segments at
the metal top level (plane + Nz MetalThickness, the sharp 90-degree "top" edges) and at the
plane (the metal / trench "bottom" edges) - lies in [r_(k-1), r_k) of the production ring set
(7 rings: 0.25, 0.75, 1.75, 3.75, 7.75, 15.75, 31.75 nm), the far shell beyond the tube radius;
label = 10000 x ordinal + parent (ordinal 1 far, 1 + k top ring k, 1 + K + k bottom ring k),
each (shell, conductor) its own boundary label of the same interface type MA. Assertions:
the $Nodes block byte-identical and every element identical in type / nodes / tag count with
only the (physical, elementary) pair of the 14,575 MA elements changed (29,150 integers;
the parent names leave $PhysicalNames, the 30 shell names enter); the shells sum to the
parent areas (4.9 um^2 each, 1.9e-14 relative vs the parent's ownership partition record);
every relabeled element's parent equals the parent ownership certificate's attribute; the
radial partition audited point-wise (degree-4 positive rules; closure 0; straddling measure
0.11% of the MA, in the graded corner-ball / cap triangles - tube quadrangles are
ring-aligned). The parent's mixed-element, protected-surface and ownership-quadrature
audits hold by identity (geometry, connectivity, element order and every other label
unchanged; Gmsh reopens the file with 34 surface groups and the parent's min SICN). Census
`identity.msh.radial-shells.json` (ring radii, edge lines, per-shell label / parent / kind /
ring / radii / element counts / area / straddling), build record `library-build.json`
(`Status relabeled`, the parent's entity counts, the `Relabel` binding) under
`/tmp/coupon-calibration-radial-ma-20260921`.

Manifest: case `two-edge-calib-radial-ma-shells` of the sizing calibration manifest
(`Calibration.Relabel` - Kind, Tool, BaseCase, ParentMeshSHA256, RingRadii - instead of
`BuildCommandOptions`; `general_mesh_manifest.validate_calibration_relabel`;
`run_gmsh_only_case.case_build_options` refuses to build it) with the labeled physics-run
deviation `Calibration.PhysicsRun {LinearTol 1e-8}` (`GateDeviations.LinearTol`, Production
1e-10, ProductionUse FORBIDDEN, `validate_case_linear_tol`;
`qualify/case_inputs.physics_run_parameters(manifest, case)` applies it only when labeled):
the acceptance run PBS 46023 ran the producer default 1e-8, and the label-only proof at the
physics level is the summed shells reproducing that run at the deterministic-solve level,
which the recipe's 1e-10 would mask by the recorded sub-4e-4% Tol effect (immaterial for
the per-shell profile). `qualify` on a relabeled case: the base config is derived against
the parent labels (it equals the acceptance p5 worker config apart from paths / Order) and
compared with the reference; `case_inputs.expand_radial_shells` then replaces the parent
labels by the shell labels in every attribute list (Ground, TerminalAttributes, edge lists)
and the one MA Dielectric entry by one entry per shell ordinal (indices after the base
entries' maximum: 15 MA entries 3-17 next to MS 2), so a participation of the type sums the
shells (`compare_matrices`, `p_sequence`, `ma_ms_offsets` label the reference matrices by the
reference's own interface map, the run's by the shell map; per-entry rows cover the common
indices) and the per-shell matrices are recorded in the fetched `surface-response-matrix.csv`.
Analysis `qualify/radial_ma_profile.py`: per source / order / edge kind the ring energies
Q_k, the power-law fit Q_k = c L (r_k^(1+alpha) - r_(k-1)^(1+alpha)) / (1 + alpha) over rings
2-7 (the innermost ring holds the singularity a polynomial element cannot represent; alpha
with its standard error), the remainder inside ring 1 = model ring-1 energy - resolved
ring-1 energy (fitted alpha and the theoretical -2/3), the extrapolated sharp-edge MA and the
deficit, and a clean-power-law flag (log residual RMS < 0.05, local slopes within 0.25 of
alpha). Tests: `test_relabel_radial_ma_shells.py` (synthetic strip relabel, the two-edge edge
lines, the shell expansion, the Tol deviation, the two-map comparison, the fit),
`test_general_mesh_manifest.py` (the case and both deviations). Evidence:
`qualify/radial-ma-20260921/RADIAL-MA.md`.

**Outcome of decision 56 (2026-09-21, PBS 46164, 0.657 node-h, root `/tmp/library-radial-ma-01`).**
Label-only proof: the shells summed per source reproduce the production acceptance run PBS
46023 at max |rel| 4.4e-13 (Q_MA), 1.1e-13 (E), 2.8e-13 (Q_MS) at p4 and p5 - the workers are
the acceptance's to the second and PCG count (322 / 1092 s, 11.90 / 10.82 iterations), the
reducers +50 / +163 s for 16 vs 2 surface integrals - and the gated p5 verdict is the
acceptance's to the digit (Failed on p_MA, 53 +3.50%). Profile (78 free sources; every number
from `qualify/radial_ma_profile.py`, `qualify/radial-ma-20260921/radial-ma-profile.md`): the sharp
90-degree top edge follows a r^(-2/3) tail in its inner rings - the ring 2-3 slope (0.25-1.75 nm)
is within 0.05 of -2/3 at every source, the 3-4 slope (1.75-3.75 nm) within 0.10 at 56 / 78 (49 /
78 within 0.05 on both); median local slopes -0.68 / -0.64 / -0.67 between rings 2-3 / 3-4 / 4-5,
fitted alpha over rings 2-4 -0.648 (quartiles -0.657 / -0.613, se 0.006), 53 / 58 -0.672 +/- 0.011;
the outer rings (4-32 nm) bend as the finite geometry takes over (flatter at the cut-face hats,
steeper at 53 / 58), which is why a fit over rings 2-7 (-0.62 median, -0.70 at 53 / 58) is a mixture
and 22 / 78 top profiles are not "clean" there (13 near-junction hats, 4 junction rings, 4
cut-face wide hats, 1 junction column; their inner slopes are still -0.68 / -0.64 ... -0.51). The
innermost ring (0-0.25 nm, 5.7% of the MA at the median source, 10.9% at 53 / 58) holds 0.706
(p5; range 0.687-0.710 over all sources) / 0.665 (p4) of the energy the -2/3 law anchored on
ring 2 predicts - a discretization property, not a source property; the p4 -> p5 MA step (+0.79%
median) is +6.2% in top ring 1, <= 0.13% in rings 2-7 (p-converged), +0.9% in the far shell (84%
ring 1 / 20% far at 53 / 58). The bottom (metal / trench) edge is not a -2/3 edge: alpha -0.33
(Fit2-4 over the 45 sources whose bottom rings carry > 20% of the MA), ring-1 p-step +1.4%,
remainder 0.1% of the MA. Implied deficit of the production value below its sharp-edge limit
(p5): 1.9-2.7% at the median source and 4.4-4.8% at 53 / 58 across the tool's estimators; the
named `Consistent` estimator (top edge -2/3 anchored on ring 2 + bottom edge at its own fitted
law) gives +2.68% (quartiles +2.22 / +3.21%), strongest-20 +2.65%, 53 / 58 +4.44% (p4: +3.07% /
+5.11%), Fit2-4 +1.89% / +4.75%. `Consistent` is preferred because 0.206 x its deficit
reproduces the ring-refined run's ring-halving step (0.50% / 0.91% vs +0.54% / +0.96%; Fit2-4
0.39%) - one check, not two: the HE-CHECK model's 2.7% / 4.7% were derived from that same step,
so its "confirmation" IS the ring-halving match. Evidence for the options (no decision): (i)
analytic extrapolation is computable per source from the run's own per-ring labels with a
source-independent ring-1 factor, but it removes the ring-1 deficit only (the far shell keeps
20-46% of the p-step), rests on the -2/3 assumption inside 0.25 nm (estimator spread 0.8 points
at the median: the exponent -0.648 vs -2/3 and the amplitude anchor) and must be applied to the
0.5 nm-anchored reference too (~5.6% at 53 / 58 under the same law); (ii) a rounding radius of
nm order acts on 5.7 / 10.9% (inside 0.25 nm) to 9.3 / 17.6% (inside 0.75 nm) of the MA at the
median / 53 - a process quantity of several percent (`qualify/radial-ma-20260921/RADIAL-MA.md`).
Review follow-up (decision 58 notes, same day): the headline estimator is the tool's
`Consistent` (COMBINED_ESTIMATORS; `HeadlineEstimator` in the JSON) and the bottom-share
floor a recorded rule (`BottomShareFloor` 0.20); `relabel_radial_ma_shells.closure_check`
closes the quadrature weights against the elements' cross-product areas (`element_area`;
1.1e-15 on the production mesh) instead of the same sum on both sides; `analyze()` /
`markdown()` guard every None fit and are tested on a synthetic shell run
(`test_relabel_radial_ma_shells.py`); `analyze_basis_needles.delaunay_faces` flips every
all-boundary box face where the producer flips the two caps only - identical on the
three-level boxes of every coupon, recorded as a difference for two-level boxes.

### Radial MA shells on every production coupon (supervisor decision 61a, 2026-09-21)

User decision 60(1) adopts the analytic sharp-edge tail extrapolation, per source from the
per-ring MA labels; decision 61a puts the labels into the recipe: every production coupon is
published with the per-ring MA shells of `relabel_radial_ma_shells.py`, applied label-only at
the placement stage. `publish_rigid_coupon_mesh.py` (proper-rigid-publication) transforms the
canonical mesh, asserts the exact structure against it as before (`ParentLabeledMeshSHA256` in
the receipt), then relabels the MA surfaces of the published bytes into the shells: the ring
radii are the gmsh-build census's `PrismTubes.Section.RingRadii` bound through the canonical
build record (a canonical build without a gmsh-build census - the legacy MMG pipeline -
publishes no shell and records why; a gmsh-build census without a ring set fails closed), the
metal edge lines come from the bound boundary / signature / process, and the centroids and
quadrature points are classified in source-local coordinates through the inverse rigid map
(so the rotated variant carries the identity's shells element by element; asserted in
`test_relabel_radial_ma_shells.ProductionRadialShellsTest`). The census
`<variant>.msh.radial-shells.json` (ring radii, edge lines, per-shell label / parent / kind /
ring / radii / counts / area / straddling, the byte-level `LabelOnly` assertions, the shells
summed against the canonical ownership partition's parent areas, the quadrature closure) is
bound by SHA-256 in the transform receipt (`RadialShells`), and the publisher re-reads the
shelled file to assert that its parent view equals the parent-labeled publication in every
cell block, connectivity and physical label and in every material / label measure. The stage
contract (inputs, artifacts, tool roles) is unchanged; the tool hashes are refrozen and
`interface_ownership.jl` / `relabel_radial_ma_shells.py` join the manifest `Tools`.

The auditors run on the shell labels, they do not skip them: the family rule of
`interface_ownership.jl` reads a surface attribute >= 10000 as a shell (`shell_parent` =
attribute mod 10000, `shell_prefix` = the rest), classifies and certifies the parent (its
family, conductor, slot) and returns the shell prefix plus the re-derived parent, so
`audit_rigid_coupon_ownership.jl` reproduces the published shell labels element by element
(the decision-56 mesh: 30 shell labels, 39,858 / 39,858 certified) and its partition /
quadrature CSVs are per shell label; `general_mesh_audit_producer._ownership_report` sums
the shell rows per parent against the semantic contract and records the shell rows
(`OwnershipClosure.RadialShells`). Every Python audit reads the mesh through
`mixed_mesh.parent_label_view` (`read_audit_mesh`; the verifier's recomputation too), so the
contract, adjacency, protected-surface, invariants and covariance audits judge the parent
surfaces as before - the shell partition is label-only and is audited by its own census. The
library build record binds the identity census (`RadialShells`: Kind, Shells {Path, SHA256},
ShellCount, RingRadii, MAParents, ParentMesh) and `coupon-library qualify` consumes it
exactly as the calibration `Relabel` binding (one MA interface per shell, per-type sums, the
base config for the reference comparison; a shelled mesh without a bound census fails closed
in `case_inputs.derive`).

**The sharp-edge MA in qualify (user decision 60(1); `qualify/ma_tail.py`).** After the
reduction qualify computes, per source and per main order, `MA_raw` = the sum of the per-ring
MA shells, `MA_tail` = the remainder of the sharp-edge law inside the innermost ring and
`MA_sharp = MA_raw + MA_tail`. The tail is the `Consistent` estimator of the decision-56
evidence (`radial_ma_profile.py`, RADIAL-MA.md): the top (sharp 90-degree) edge's r^(-2/3) law
anchored on ring 2 - Q_1 model = Q_2 r_1^(1/3) / (r_2^(1/3) - r_1^(1/3)), remainder = model
minus the resolved ring 1, the ring-1 factor (resolved / model) measured on THIS run at each
order and never carried over - plus the bottom (metal / trench) edge's own law fitted over
rings 2-4. Per source the record keeps alpha of the top edge fitted over rings 2-4 with its
standard error, the ring-1 factor, MA_raw / MA_tail / MA_sharp and the estimator spread (the
Fit2-4 and Theory@2 deficits next to `Consistent`: the honest uncertainty of the tail;
`comparison/ma-tail.{json,md}`, `ma-sharp-<stage>.csv`, the library record's `MATail` and
`Qualification.MASharp`, the process-library model's `MA {MA_raw, MA_sharp, p_MA_raw,
p_MA_sharp}` per source). The reference is extrapolated by the same rule only when its ring /
edge sizing is recorded (a shelled run: its own shells); the graded_v2 tet references record
none, so they are 'reference unextrapolated' and their sharp-edge deficit is MODELLED per
source from the edge size given as `--reference-edge-size-nm` (0.5 for graded_v2) via the
eps^(1/3) law of decision 55 - deficit_ref = deficit_run(`Consistent`) x (eps_ref / eps_run)^(1/3)
with eps_run the run's innermost ring radius - and the comparison records raw-vs-raw
(`p_MA_rel`) next to sharp-vs-modelled (`p_MA_sharp_rel`; ma-ms-offsets and key-sources carry
both). The gate table (`qualification-gates.json`, thresholds unchanged) names the quantity:
`p_MA.Quantity = p_MA_sharp` with the recorded `QuantityRule`, and the p-sequence control's
`MAObservable = p_MA_sharp` (|d45| < 5% on the extrapolated MA replaces the raw gate, whose
step was dominated by the innermost ring's resolution); the gate record carries `Quantity`,
`QuantityRule`, the raw free distribution (`RawFree`) and `MAObservable`; a coupon without
shells is gated on the raw p_MA and says so. On the stored decision-56 run the analysis
reproduces RADIAL-MA.md (p5 `Consistent` median 2.68%, 53 / 58 4.44%; Fit2-4 1.89% / 4.75%;
alpha -0.648; ring-1 factor 0.706, range 0.687-0.710; p4 3.07% / 5.11%) and the raw gate
numbers of the acceptance (median +1.25%, 53 +3.50%); sharp-vs-modelled (0.5 nm reference,
modelled deficit median 3.38%) gives median +0.55%, 51 / 76 / 78 within 1 / 2 / 5%, strongest-20
18 / 20, 53 / 58 +2.37% (`test_ma_tail.StoredRadialRunAnalysisTest`).

**Outcome (2026-09-22, `qualify/ma-sharp-20260922/MA-SHARP.md`).** The five production roots were
rebuilt with the shells and verified - two-edge 10 `5d01204e` -> `179a405d` (byte for byte the
decision-56 relabel), two-edge 05 `9f011ed3` -> `86cfad22`, four-edge `1d536c44` -> `f36d2d83`,
ten-edge `8a871f22` -> `e1a68075` (60 shells, four MA parents), three-edge 06 `95837ed5` ->
`4e2adbe1` - each labels-only against the previous identity (`verify_label_only_republication.py`:
$Nodes and every element identical apart from the MA elements' label pair, the receipt's
ParentLabeledMeshSHA256 = the previous identity). The first publication of 06 failed on a
`metal_edge_lines` defect (the conductor's process-normal set consumed by the first plan-view loop),
fixed forward with its test. Acceptance (PBS 46729 / 46730, 1.76 node-h): four-edge MA_raw = PBS
46024 at 8.5e-13 (labels only); two-edge 10 at 3.3e-4 vs PBS 46023 because 46023 ran Tol 1e-8 (its
mesh is PBS 46164's, which reproduced 46023 at 4e-13). Four-edge PASSES on MA_sharp (free median
+0.14% vs +0.79% raw, 44 / 53 / 60 within 1 / 2 / 5%, strongest-20 median +0.18%, 20 / 20 within 1%;
deficit `Consistent` 2.86% median, alpha -0.655, ring-1 factor 0.666 at p4 = two-edge's 0.665);
two-edge 10 still FAILS the strongest-20 at 53 / 58: +2.37% / +2.33% sharp-vs-modelled (raw +3.50% /
+3.46%; median +0.55%, 51 / 76 / 78). The estimator spread moves the 53 / 58 excess by 0.07 points
(Fit2-4 +2.30%, Theory@2 +2.35%): it is not explained by the tail's uncertainty; the modelled 0.5 nm
reference deficit (5.59% at 53) would have to be 5.98% for the 2% statement, i.e. the tet reference
resolving its innermost 0.5 nm 7% worse than a prism ring - its own ring-1 factor is not measured, the
excess stays reference-limited at 2.3-2.4%. Every p-sequence control on p_MA_sharp passed (max
+1.67% / +1.52%). Full `unittest discover` 342 tests OK (33 skipped); preflight passes on the three
manifests.

## Gmsh-only production pipeline (supervisor decision 38, 2026-09-18)

The coupon mesh is generated by Gmsh alone. Decision 37's prism-tube spike showed
that a localized prism tube around every straight metal edge (rings from 0.25 nm
with ratio 2, extruded along the edge at the tangential spacing - since decision 40
in layers following the size field on the axis -, explicit pyramids
on the lateral quadrangles) converges MA at 21.5M H1 DOFs (p_MA vs the 0.5 nm
reference median +0.47% / weighted +0.52% / strongest-20 +0.76%, p3/p4/p5
increments <= ~1 point, 0.57x the production DOFs), where every tetrahedral edge
layer under MMG stalled at 23% of the residual per 4x step (EL4c -6.63%, EL1c
-5.35%). The spike's E/SA regressions at the narrow hats (74/10) and the
cut/trench-junction sources were attributed to its coarse Gmsh volume (no far
growth, no junction-line bands, no tube-adjacent band): the trace-basis cut-surface
callback was present. The metric-preparation, native-adaptation-mmg and
label-restoration stages, the adapter/MMG tool roles, the required-region
optimizer and the tetrahedral edge layer are retired from production; the code and
the calibration manifest keep them (labeled `legacy-mmg`) as history.

The production canonical DAG is `canonical-source-validation -> gmsh-build ->
canonical-gmsh-publication` (one reusable source-local canonical build, cache key
from the immutable sources, gates and the three stages' tool roles), followed by
the mandatory `proper-rigid-publication` per placement, the consolidated audits,
normalization and verification. `mesh_stage_contract.pipeline_of` identifies a
pipeline by its exact stage set; a manifest's `StageToolSHA256` freezes exactly one
pipeline and an explicit `Pipeline` key must agree.

`gmsh-build` runs `mesh_spatial_coupon.jl --prism-tubes true` (the seed-generation
mesher, unified with the spike: `prism_edge_tubes.jl`) with the recipe
`ProductionRecipe.BuildCommandOptions` and writes the labeled mixed-element Gmsh 2.2
binary mesh (`gmsh-mesh`: materials 1/2, matching surface, the semantic interface
labels; multi-slot coupons through the mesher's ownership postprocessor) and the
build census (`build-census`, the bound build report):

- prism edge tubes on the top and bottom edge of every straight metal segment
  (every `Physical` side of the plan-view loops, excluding box sides): ring k has
  size EdgeSize x GrowthRatio^(k-1), the ring count is the largest K with r_K + h_K
  <= min(Overetch, MetalThickness / 2, CornerIsotropyRadius) (7 rings, radius
  31.75 nm on the four/ten-edge process), 30-degree sectors on the dielectric side
  (vacuum above the top edge; substrate / vacuum split at the trench wall below the
  bottom edge; the etch footprint must carry the metal edge), extruded in layers
  whose thickness follows the composed size field on the tube axis (supervisor
  decision 40, below; every layer <= TangentialSize; the largest layer is the
  recorded per-tube Spacing), lateral quadrangles closed by explicit pyramids of
  height 0.5 x the outermost ring size; tubes end at the outer box and
  R / tan(phi / 2) + h_K before a semantic corner (phi the in-plane angle of the
  metal edges meeting there; recorded);
- isotropic corner balls graded to the tube inner size: CornerSize == EdgeSize is
  required (one graded law; shells 0.25/0.5/1/2/4/8/16 nm to NormalSize inside the
  0.1 um ball), and every tube cap centre before a corner is a graded point of the
  same law, so the un-tubed edge part and the cap region are tetrahedra graded from
  EdgeSize - the spike's cap slivers (702 tets < SJ 0.01) are gone: the census
  `CapRegions` (tetrahedra with a vertex within the tube radius of a cap centre)
  measure min SJ 0.043 / max condition 40 on the four-edge build; the corner balls
  keep the seed-side bounded descent with the production gates (corner aspect <=
  4.0, target 3.8) with every tube node fixed;
- explicit volume size laws in the Gmsh size callback (recorded under
  `PrismTubes.SizeLaws`): tube band `size = min(FarSize, NormalSize + FarGrowth x
  max(d_axis - (R + pyramid height), 0))`; feature-curve band (junction lines,
  footprint edges, un-tubed edge parts) `size = min(FarSize, NormalSize +
  RadialGrowth x min(r, 2 NormalSize) + FarGrowth x max(r - 2 NormalSize, 0))`
  (the metric stage's band law; RadialGrowth 1 and ProtectedDistance 2 x
  NormalSize are mesher constants bound by the census validator); the trace-basis
  cut-surface rule (TraceBasisSizeRatio x local basis edge) composed by `min` with
  the background attractor / graded-point fields; the far field FarSize with the
  fail-closed element cap (no far-budget pressure: the requested far size is used
  as is and the build fails above `--max-elements`);
- volume size laws of supervisor decision 39 (physics-09 localized the E/SA/MS
  regressions of the Gmsh-only mesh to volume sizing next to the narrow-hat apexes,
  the corner balls and the junction lines): (a) the trace rule is a volume law,
  size <= TraceBasisSizeRatio x the local basis edge + FarGrowth x the distance to
  the basis triangle (`TraceBasisSizing.GradingSlope` = FarGrowth; before: the
  process-band slope 0.675); (b) the corner-ball exterior, size = NormalSize +
  FarGrowth x the distance beyond CornerIsotropyRadius from the nearest graded
  point (`SizeLaws.CornerExteriorRule`); (c) the junction lines carry the band law
  throughout the volume (`SizeLaws.JunctionVolumeRule`, the BandRule) and are
  1D-meshed at NormalSize (below). All three compose by `min` in the size callback
  and are measured in NormalSize shells (`SizeLaws.Achieved`: mean / longest edge
  percentiles and achieved-over-prescribed around the trace apexes within
  CornerIsotropyRadius, outside the corner balls, around the junction lines;
  report, not gate);
- tube layers following the size field (supervisor decision 40, 2026-09-19; the E
  +10% at the narrow-hat sources 74/10 sits where a tube terminates on the cut
  surface at (2,8), its extrusion normal to the Dirichlet surface, and the 50 nm
  prism layers under-resolved the 20-50 nm decay of the hats into the volume): the
  layer boundaries of every tube equidistribute the arclength integral of the
  reciprocal of the composed size field evaluated on the tube axis - `min(
  TangentialSize, the corner-ball law of the semantic corners along the edge, the
  trace-basis volume rule, the band rule, the corner-exterior rule)`
  (`PrismTubes.TubeAxisSizeLaw`) - gradient-limited along the axis to
  (GrowthRatio - 1) / GrowthRatio, with ceil(integral) layers, so a layer is the
  size at its midpoint, never above TangentialSize, and consecutive layers differ
  by at most GrowthRatio (`PrismTubes.LayerRule`; `LayerGrowthCap` = GrowthRatio;
  the mesher fails closed above it). No new parameter: the same laws as everywhere
  else; the cross-section rings, the pyramids and the tetrahedral interface are
  unchanged. Excluded from the axis law: the tube rule (the tube's exterior, which
  reads NormalSize on the axis) and the tube cap centres as ball-law points -
  grading the layers to CornerSize at a cap makes the lateral pyramid faces
  CornerSize x outer-arc slivers by construction, and the corner-ball tetrahedra
  against them fail the scaled-Jacobian gate (four-edge probe: 38 cells below 0.01,
  minimum 0.0060); the cap centres act on the axis through the corner-exterior rule
  only (NormalSize within CornerIsotropyRadius of a cap). The census records per
  tube (`PrismTubes.Tubes[].LayerThickness`) the layer count, thickness minimum /
  P50 / maximum, the thickness and the prescribed size at both ends, the
  achieved-over-prescribed range (layer over the gradient-limited size at its
  midpoint) and the neighbour ratio, the tube end points with `EndsOnBox`, and over
  all tubes `PrismTubes.LayerThickness` (with
  `LayersBelowTangentialSizeOverGrowthRatio`, the layers the axis field refines
  below TangentialSize / GrowthRatio); `validate_gmsh_build_census` requires the
  rule strings, `LayerGrowthCap` = `--edge-growth-ratio`, every neighbour ratio
  within it, the thickness order Minimum <= P50 <= Maximum <= TangentialSize with
  `SpacingMinimum` / `SpacingMaximum` equal to the extremes, the below-count within
  the layer count, and every row's record with its Maximum equal to the row's
  Spacing and, at an end on the box, the end layer within the prescribed size.
  Decision 41 (surface layer): at a tube end on the outer box the end layer is the
  field evaluated at the surface (at most the gradient-limited field over the layer
  span, iterated from the surface value), not the size at the layer midpoint, so
  that a field growing away from the cut surface (the trace rule at FarGrowth) is
  resolved from the surface; the rest of the tube is equidistributed as before. No
  new parameter. Four-edge (probe, 083874fee): layers at the
  (2,8) box end 28.2 nm (prescribed 21.7 nm by the trace rule, growing at
  FarGrowth; 21.7 nm under decision 41), 25 nm at the corner ends, 16 nm at the bottom-tube ends 16 nm before
  the on-box corners (10,0) / (0,8), 49.9 nm mid-tube; 1436 -> 1466 layers;
- band curves 1D-meshed at the band law (review P1 of the phase-3 evidence): the
  non-metal longitudinal feature curves - the cut-surface / material-interface
  junction lines and the footprint edges parallel to a metal edge - are placed at
  NormalSize spacing (the band law on the line, composed with the corner law), so
  the first cell layer against a junction is NormalSize transversally; on the
  lc_tangent ridge grid it was TangentialSize (four-edge cut surface near junctions
  P50 49.8 nm under commit 1a0289994), because the surface mesher cannot refine a
  curve's nodes. The metal ridges keep the tangential grid (they carry the tubes).
  The census records the band curves (`PrismTubes.BandCurves`: count, length,
  Spacing = NormalSize, segments) and the achieved first-layer transverse size
  against the junction lines (`PrismTubes.Bands.JunctionFirstLayer`: P10/P50/P90/
  maximum and achieved-over-prescribed for the matching-surface elements and the
  tetrahedra with a node on a line; report, not gate);
- the census computed by the mesher itself (Gmsh's `minSJ` is the high-order
  mapping Jacobian, identically 1 for order-1 elements): per element type the
  corner-frame scaled Jacobian, Jacobian condition and orientation (tetrahedra,
  prisms, pyramids; gated fail closed: positive orientation and condition <= 1000
  for every type, SJ >= 0.01 for tetrahedra), tube rings / prisms / pyramids /
  layers / spacing / prism edge aspect, cap-region quality, cut-surface edge
  statistics (all, near junctions, achieved over requested against the trace rule:
  measured only on the narrow basis triangles that contain a matching-element
  centroid, the unmeasured count reported), junction / footprint band tetrahedron
  sizes (the footprint statistic over the footprint sides not on the outer box -
  the feature curves - with the box sides reported separately as
  `FootprintOnBox`, a far-field sample), the junction first-layer sizes, size
  laws, corner shells, interface areas with quadrangles, footprint simplification,
  straight junction segments, element counts per type and the element-budget
  record.

`mesh_stage_contract.validate_gmsh_build_census` binds the census to the build
command and the canonical semantic contract (InnerSize = `--edge-size` =
CornerSize, GrowthRatio, TangentialSize = `--lc-tangent`, NormalSize = `--lc-fine`,
FarSize = `--lc-far`, FarGrowth, `--prism-tubes true`, rings geometric, the sector
angle = `--tube-sector-degrees` or its default 30, PyramidHeight = 0.5 x the
outermost ring (PyramidHeightOverOuterRing), the band law's RadialGrowth 1 and
ProtectedDistance 2 x NormalSize, the decision-39 volume laws growing with
FarGrowth (trace GradingSlope, CornerExteriorGrowth, CornerIsotropyRadius, the
achieved shells recorded), band curves at NormalSize spacing with the
junction first-layer records present, every tube
spacing within the tangential size, the decision-40 layer records (rule strings,
`LayerGrowthCap` = GrowthRatio, neighbour ratios within it, thickness extremes =
`SpacingMinimum` / `SpacingMaximum`), per-type quality within the command's gates,
cap regions within the tetrahedral gates, corner aspects within the corner gate,
element count within `--max-elements`, footprint provenance, junction segments
straight, trace basis iff bound, no edge layer). `canonical-gmsh-publication`
(`relabel_frozen_interface_mesh.jl`) applies the exact ownership partition to the
Gmsh mesh: `interface_ownership.jl` / `label_interface_patches.jl` accept linear
quadrangles (the tube's radial faces on the metal and trench-wall strips) with the
positive-weight Gauss4 rule of each element type and four-node certificates.

Audits on mixed meshes (`mixed_mesh.py`): the topology, adjacency, material volume,
label area, protected-surface and trace-diagonal audits run on a conforming
simplicial view (prisms -> 3, pyramids -> 2 tetrahedra, quadrangles -> 2 triangles
by the minimum-vertex diagonal rule; exact measures for the planar-faced tube
cells); element quality is measured on the native cells per type
(`MeshQuality.ByType`; top level: orientation and condition over every type, scaled
Jacobian over the tetrahedra; the gate `mesh-quality-jacobian` judges both); H1
counts and element counts are native. The achieved-anisotropy design gate is
replaced by the prism tube design statement (`AchievedAnisotropy.Gate =
"not-applicable: prism edge tubes"`: the recorded rings / extrusion / prisms /
pyramids / spacing / cap and cut-surface statistics, the mesh's prism and pyramid
counts required to equal the census; informational). The trace-diagonal detector
judges a band on a trace-basis edge within the band's own direction resolvability
(decision 36) and its endpoints within 2 x ShortEdgeThreshold + RMSWidth of the
edge segment (the four-edge narrow-hat bands on the top/bottom cut faces: 11.4 um
on the basis edges (2,8)-(-6,0), RMS width 33 nm, tilt ~3 mrad). Physical gates
unchanged: positive orientation and Jacobian condition <= 1000 for every element
type, scaled Jacobian >= 0.01 for tetrahedra, MaximumCornerAspect 4.0,
protected-surface tolerance 1e-8, ownership closure 1e-12, trace-diagonal 0; the
resource guidelines were MaximumElements 4,000,000, 1800 s / 8 GiB per stage (audits
16 GiB) at the time (6,000,000 / 3600 s / 12 GiB since decision 61c).

`run_gmsh_only_case.py CASE` drives a production case through the DAG, the rigid
placements, the consolidated audits, normalization and per-entry verification from
the manifest case alone (roots `/tmp/coupon-gmsh-only-<case>-<commit>-<timestamp>`).

### Evidence (Gmsh-only four-edge and ten-edge, commit 083874fee, 2026-09-19)

Both production cases were rebuilt under the committed tools by
`run_gmsh_only_case.py` after the band-curve fix (review P1 of the 1a0289994
evidence: junction curves 1D-meshed at NormalSize, commit 8022ade03), the recipe
bindings (review P2s, 3852dbd47), the decision-39 volume size laws (f85c30932)
and the decision-40 tube layers following the size field on the axis (083874fee),
and verified identity + rotate-z with empty failure lists
(`per-entry-verification.json`: `Passed true`, `Failures []`,
`TransformComparisonFailures []`, one shared CanonicalBuildId per case).

| case | elements (tets + prisms + pyramids) | nodes = H1 p1 (H1 p4) | tubes / layers / layer thickness min-P50-max | rings | cap regions (min SJ / max cond) | tets min SJ / max cond | prisms max cond | pyramids max cond | corners (10 / 4) | protected | closure | diagonal bands | stage s / GiB (build; publication) | audits s / GiB | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| four-edge-9d2cb9bbb3fe | **1,742,434** = 1,557,718 + 171,522 + 13,194 (+0.43% vs f85c30932's 1,734,992; +3.7% vs 1a0289994's 1,680,422; 0.50x the 34B production 3,471,480) | 383,805 (23.0M) | 8 / 1,466 / 15.9 - 49.86 - 49.99 nm (f85c30932: 1,436 uniform layers of 49.70-49.93 nm) | 7 (0.25 ... 16 nm, R 31.75 nm) | 12: 0.0548 / 20.7 | 0.0442 / 42.6 (f85c30932: 0.0270 / 86.0) | 589.4 | 8.74 | 2.74 / 3.59 / 2.84 / 3.25 | 0 (support vertices 1e-15) | 8.4e-14 (Gauss4, 633,714 points) | 0 (8 signature-, 30 footprint-, 8 junction-aligned; 0 on trace-basis edges: ShortEdgeThreshold 0.032, below the narrow-hat band width) | 137.4 / 4.13 (gmsh-build), 31.7 / 3.20 (publication); 122-129 / 3.9 per placement | 291-298 / 4.1 | 877 / 2.8 |
| ten-edge-6791f1c84123 | **2,465,185** = 2,175,133 + 269,334 + 20,718 (+0.89% vs 2,443,456; 0.69x the 34B production 3,570,533) | 556,661 (33.2M) | 20 / 2,302 / 15.9 - 49.82 - 49.94 nm (before 2,224 uniform layers of 47.6-49.9 nm) | 7 | 36: 0.0401 / 43.4 | 0.0217 / 101.2 (before 0.0218 / 93.4) | 588.8 | 8.73 | 3.34-3.78 (10 corners) | 0 (1.8e-15) | 7.8e-14 (970,074 points; 10 owner labels: 3100/3101, 5001/5002, 5101/5102, 6001/6002, 6101/6102) | 0 (10 signature-, 10 footprint-, 6 junction-aligned) | 230.5 / 4.53 (gmsh-build), 50.6 / 4.26 (publication); 175-180 / 5.0 per placement | 457-461 / 5.0 | 1,111 / 3.8 |

Decision-40 tube layers (census `PrismTubes.LayerThickness`, per tube
`Tubes[].LayerThickness`; the layer is the gradient-limited axis size at its
midpoint, achieved-over-prescribed P50 0.994-1.000 on every tube; the largest
neighbour ratio 1.59 (four) / 1.53 (ten) against the cap GrowthRatio 2): four-edge
- the two tubes ending on the cut surface at (2,8) where the narrow-hat sources
74/10 live: first layer 28.2 nm against the trace rule's 21.7 nm on the surface
(the layer grows at FarGrowth 0.5 from the surface: 21.7 + 0.5 x 14 nm at its
midpoint), then 25 nm inside the corner ball at the (2,0) end; the two tubes ending
on the box at (10,-2): 49.9 nm there (no narrow hat; the far size), 25 nm at the
(0,-2) corner end; the (2,0)->(10,0) and (0,-2)->(0,8) tubes: 25 nm at the
corner-clearance ends, 15.9 nm at the bottom-tube ends 16 nm before the on-box
corners (10,0) / (0,8) (the corner-ball shell size there; the top tubes are
sqrt(16^2 + 100^2) nm from those corners, 25 nm); 49.9 nm mid-tube; per tube 162-205
layers (before 159-200). Ten-edge: 40 tube ends - 4 on the box at 49.8 nm (50 nm
prescribed), 34 corner ends at 24.9-25.0 nm, 2 at 15.9 nm (16 nm); per tube 23-165
layers. (The census field `LayersBelowTangentialSize` of these builds counted every
layer - every layer is strictly below the cap by construction - and was replaced by
`LayersBelowTangentialSizeOverGrowthRatio`, the field-refined layers.)
Trace-apex shells (below) tightened: four-edge 25-50 / 50-100 nm longest edge P50
38.6 -> 32.6 nm / 49.3 -> 42.7 nm, ten-edge 49.0 -> 35.8 nm / 62.9 -> 54.0 nm; the
four-edge tetrahedra minimum scaled Jacobian rose from 0.027 to 0.044 (25 nm layers
against the cap-graded tetrahedra instead of 50 nm). Not viable (probe under the
same commit's tools, recorded here): grading the layers to CornerSize at the cap
centres - 38 corner-ball tetrahedra below the SJ gate (minimum 0.0060) after the
seed optimization; hence the axis law excludes the cap centres as ball-law points.

H1 p1 is `Measurements.Complexity.H1DOFs` of the identity `mesh-complexity`
record (= the node count); H1 p4 is `mixed_mesh.h1_dofs(mesh, 4)` on `identity.msh`
(not a recorded field). The table below and the decision-44 evidence table were
computed before `mixed_mesh.h1_dofs` adopted Palace's pyramid interior count
(`(p - 1)^3`, the Fuentes H1 pyramid of Palace's MFEM build; the earlier Bergot count
`(p - 1)(p - 2)(2p - 3)/6` undercounts by 22 x pyramids at p4): the Palace-printed H1
of the three-edge root is 24,844,050 = 24,551,802 + 22 x 13,284, and
`h1_dofs_from_counts` now reproduces the printed p3 / p4 / p5 counts exactly (test).
The `coupon-library build` record stores the H1 entity counts so `qualify` estimates
every order from them.

Junction first layer (census `PrismTubes.Bands.JunctionFirstLayer`, prescribed
NormalSize 25 nm; before = the 1a0289994 builds, whose junction curves sat on the
50 nm ridge grid): four-edge cut-surface elements with a node on a junction line,
transverse P50 / P90 25.1 / 28.8 nm (achieved over prescribed 1.005 / 1.15;
before: cut surface near junctions edge P50 49.8 nm), tetrahedra 27.9 / 29.8 nm;
ten-edge 25.1 / 28.4 nm (before P50 49.4 nm), tetrahedra 27.9 / 29.4 nm. Band
curves: four-edge 24 (98.4 um: the 12 horizontal junction lines and the 12
axis-parallel footprint edges), ten-edge 12 (80.0 um).

Decision-39 volume laws (census `PrismTubes.SizeLaws.Achieved`, NormalSize shells,
longest edge P50 / achieved-over-prescribed P50 per shell): trace apexes
(shortest-edge endpoints of the narrow basis triangles; the physics-09 observable
was the median longest edge within 0.1 um of the narrow-hat apexes, 57 nm on the
1a0289994 mesh vs 42 nm on the MMG mesh) - four-edge 54 apexes, 0-25 / 25-50 /
50-100 nm shells: 4.9 / 32.6 / 42.7 nm (2,517 / 1,437 / 4,002 cells; ratios 0.05 /
0.45 / 0.45; f85c30932: 4.8 / 38.6 / 49.3 nm; the innermost shell sits inside the
tube of the metal edge the hats touch), ten-edge 180 apexes: 4.9 / 35.8 / 54.0 nm
(ratios 0.06 / 0.45 / 0.53; before 4.7 / 49.0 / 62.9 nm); corner exterior (R = 0.1
um; shells R + 0-50 / 50-100 / 100-200 nm) four-edge 33.5 / 45.0 / 67.4 nm (before
41.6 / 51.1 / 69.5; the anisotropic attractor field is finer than the exterior law,
so the law is a recorded bound the mesh already satisfies), ten-edge 33.5 / 44.0 /
66.6 nm; junction lines (0-25 / 25-50 / 50-100 / 100-200 nm) four-edge 41.0 / 75.2 /
104 / 141 nm, ten-edge 40.0 / 75.0 / 102 / 121 nm. `TraceBasisSizing.GradingSlope`
is FarGrowth 0.5 (before 0.675).

Cut-surface sizes (`PrismTubes.CutSurface`): four-edge matching surface edge P50
0.160 (far); the trace rule is measured on the narrow basis triangles that contain
a matching-element centroid: four-edge 60 of 76 measured (16 unmeasured),
achieved/requested P50 0.64, maximum 1.12; ten-edge 147 of 244 measured (97
unmeasured), P50 1.00, maximum 1.20 - some narrow triangles are met only within
~20% (recorded, not gated). Footprint band tetrahedra (`Bands.Footprint`, the
footprint sides not on the outer box; the box sides separately as
`Bands.FootprintOnBox`): the four-edge interior sides carry the band (mean edge
P50 ~23 nm); the ten-edge producer-default footprint's interior sides are
conductor-strip boundaries at z = 0 without a material step (the ten-edge coupon
has no 3000/3001 un-etched label), so they carry no band (P50 ~178 nm; the 228 nm
maximum of the 1a0289994 statistic came from these sides, not from the box sides).
The tube band at the tube surface starts at NormalSize and grows with FarGrowth
0.5 to FarSize 0.16. Roots (binaries kept):
`/tmp/coupon-gmsh-only-four-edge-9d2cb9bbb3fe-083874fee-20260918-225243`
(identity.msh SHA-256
`647b2079786f02784a1e855ae197d874df19598e1f33cae40fe25b6d357e4c6f`, rotate-z
`c0a2b136...`, CanonicalBuildId
`1675018f3cac4a80194c0464a0648178214f80880f9c010eef63d6bed5011db7`),
`/tmp/coupon-gmsh-only-ten-edge-6791f1c84123-083874fee-20260918-225246`
(identity.msh SHA-256
`eb9986a40f6bb4e97a169fe478aa2fe647bf363ba2a1cf5f82306f5a1c3a5f73`, rotate-z
`6864b616...`, CanonicalBuildId
`15fb423f9f8a1f418f20fd9581dfbcef20605204609348bd9c9139f692c25d61`). The superseded
1a0289994 roots (four-edge identity `e57ce087...`, 1,680,422 elements, the mesh
of physics-09; ten-edge `7ff972ab...`), the intermediate 3852dbd47 roots (four
`dfb58639...`, 1,714,641; ten `9ff5ee0d...`, 2,422,099; band curves without the
volume laws) and the f85c30932 roots (four `a94983da...`, 1,734,992, CanonicalBuildId
`e8931bb1...`; ten `897c46bc...`, 2,443,456, `27ee76e3...`; uniform tube layers)
keep their records; their mesh binaries were removed once these roots verified. Physics-09 (on the 1a0289994 four-edge mesh) adjudicated the
junction SA against the 50 nm first layer and localized the E/SA/MS regressions to
the volume sizing these meshes now carry.

### Evidence (decision 41 sizing calibration variants V-a / V-b, commit efd393f38, 2026-09-19)

Physics-10 on the decision-40 four-edge mesh (`647b2079...`) kept MA converged but
left E +6.5/+7.3% at the narrow-hat sources 74/10 (cells within 0.1 um of the apexes:
146/134 at median longest edge 53/57 nm against EL1c's 232/241 at 42/44 nm; the first
prism layer at the (2,8) cut 28.2 nm against 21.7 nm prescribed), junction p_SA
+2..+3% and the z = 0 junction ring p_MS +5..+7%. Decision 41 calibrates the two
recorded recipe parameters that control these volumes, each together with the
surface layer at on-box tube ends (f7c23957b, above): V-a `--trace-basis-size-ratio
0.5` and V-b `--lc-tangent 0.025`, as the labeled cases of
`geometry-independence-calibration-sizing.json`. Mesher-only probes (f7c23957b) first
bounded the element counts under the 4M cap (V-a 1,869,209, 121 s / 3.9 GB; V-b
2,119,019, 126 s / 4.2 GB), then both cases were built through
`run_gmsh_only_case.py` under the calibration manifest (identity + rotate-z-0.63,
audits, per-entry verification: `Passed true`, `Failures []`,
`TransformComparisonFailures []`, `CanonicalReuseFailures []`; every physical gate at
its production value). The production manifest was unchanged at efd393f38 (ratio 1.0 / 50 nm); decision 42 adopted the V-a ratio afterwards (above).

| case | elements (tets + prisms + pyramids) | nodes = H1 p1 (H1 p3 / p4 / p5, exact hybrid entity sums) | tubes / layers / thickness min-P50-max (nm) / max neighbour ratio / below TangentialSize / 2 | first layer at the (2,8) cut (top, bottom) | trace-apex cells within 0.1 um of 74 / 10 (median longest edge) | junction first layer cut-surface P50 / P90 (nm) | tets min SJ / max cond | prisms / pyramids max cond | caps (12) min SJ / max cond | corners (4) | protected / closure / diagonal bands | stage s / GiB (build; publication; placement) | audits s / GiB (identity; rotate-z) | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| four-edge-calib-sizing-tbr-0.5 (V-a) | **1,869,209** = 1,683,863 + 172,107 + 13,239 (+7.3% vs 083874fee's 1,742,434) | 414,758 (10.55M / 24.83M / 48.29M; p4 +6.5%) | 8 / 1,471 / 10.85 - 49.75 - 50.0 / 1.93 (cap 2) / 52 | **10.85 nm** = the trace rule at ratio 0.5 on the surface (prescribed 10.85) | **575 / 533 at 29.5 / 31.6 nm** (r < 0.05: 170 / 161 at 25.0 / 26.2; r < 0.25: 1,918 / 1,824 at 36.2 / 35.8); cut-surface triangles within 0.1 um 257 / 249 at 12.7 / 13.5 nm | 21.6 / 27.9 (A/P 0.87: the trace rule refines the cut under the junctions too) | 0.0296 / 118.1 | 668.6 / 8.74 | 0.0305 / 67.6 | 3.13 / 3.44 / 3.47 / 3.54 | 0 (9.9e-16) / 2.5e-13 (653,577 points) / 0 | 115.1 / 4.27; 27.8 / 3.18; 104-106 / 3.82 | 265 / 2.82; 322 / 2.86 | 972 / 2.92 |
| four-edge-calib-sizing-lct-0.025 (V-b) | **2,119,019** = 1,757,651 + 335,556 + 25,812 (+21.6%) | 514,194 (13.31M / 31.43M / 61.25M; p4 +34.8%) | 8 / 2,868 / 16.0 - 24.94 - 25.0 / 1.34 / 0 | **21.7 nm** = the trace rule at ratio 1.0 on the surface (prescribed 21.7; 28.2 before f7c23957b) | 146 / 134 at 52.1 / 56.8 nm (= the 083874fee mesh: the apexes at z = +/-2.1 are not on a tube); source 26 (z = -0.05, on the tube): 907 at 24.9 nm (before 613 at 44 nm) | 25.1 / 28.8 | 0.0267 / 87.8 | 299.4 / 4.48 | 0.0318 / 41.7 | 3.11 / 3.76 / 2.89 / 3.58 | 0 (9.9e-16) / 1.4e-13 (922,506 points) / 0 | 122.3 / 4.76; 29.6 / 3.74; 117-123 / 4.47 | 289 / 3.22; 378 / 3.20 | 1,085 / 3.25 |

Trace-apex census (the physics-10 statistic: volume cells with centroid within r of
the source apex, `matching_surface_local_sizing.py` of the assessment on the identity
meshes): V-a multiplies the cells within 0.1 um of 74 / 10 by 3.9 / 4.0 and halves
the median longest edge (53 / 57 -> 29.5 / 31.6 nm; EL1c 42 / 44 nm at 232 / 241
cells), the census trace-apex shells 25-50 / 50-100 nm move 32.6 -> 32.1 / 42.7 ->
41.3 nm (binned over all 54 apexes); V-b leaves the apex volumes unchanged and
halves the layers along every tube (2,868 layers of 25 nm; prism max condition 589 ->
299, pyramid 8.74 -> 4.48). Tube ends: both variants place the surface layer at every
on-box end (V-a (2,8) 10.85 nm, (10,-2) 50.0 / 25.0 nm; V-b 21.7 nm, 25.0 nm), the
corner-clearance ends at 24.9-25.0 nm and the bottom-tube ends before the on-box
corners at 15.9-16.0 nm; V-a's largest neighbour ratio 1.93 is the surface layer
against the trace rule growing at FarGrowth 0.5 (the analytic bound 1.95 at
GrowthRatio 2, prism_edge_tubes.jl). Both cases share CanonicalBuildId
`41e4f24391f99725716443a2eea578c00f419f0471c3af075c9d806990a55b1e` (the cache key
does not encode recipe options, the documented calibration limitation: the label is
bound to the build by the recorded gmsh-build command). Roots (binaries kept for
physics-11/12): `/tmp/coupon-gmsh-only-four-edge-calib-sizing-tbr-0.5-efd393f38-20260919-021146`
(identity.msh SHA-256
`80966c7db44dabc49ac7bb068bbaee0e0828e413b115cee8c066c886866b6108`, rotate-z
`97bfb753...`) and
`/tmp/coupon-gmsh-only-four-edge-calib-sizing-lct-0.025-efd393f38-20260919-021148`
(identity.msh `7db360bed45f2c7f94bf6857c603ed212ef6a1699ea0dce6e26677e2e4f98615`,
rotate-z `5b000119...`); the mesher-only probe meshes were deleted. Physics-11/12
adjudicated E at 74/10, SA at 48/47/35 and MS against physics-10 (table above):
V-a adopted (decision 42), V-b rejected.

### Evidence (decision 42 production roots: four-edge and ten-edge at TraceBasisSizeRatio 0.5, commit f6efe6367, 2026-09-19)

Both production cases were rebuilt by `run_gmsh_only_case.py` under the adopted
recipe (production manifest with `--trace-basis-size-ratio 0.5`) and verified
identity + rotate-z-0.63 with empty failure lists (`per-entry-verification.json`:
`Passed true`, `Failures []`, `TransformComparisonFailures []`,
`CanonicalReuseFailures []`). The roots were launched from the f6efe6367 working
tree before that commit was recorded, so their directory names carry the parent
21cf534ff; the evidence is bound to f6efe6367 by content: the verification
reports record the manifest SHA-256 `04f8df72434eb08c5089805510918c504a7f5a74
32762bae240cd19f958c00b5` (= `geometry-independence-suite.json` at f6efe6367) and
the refrozen tool digests (`general_mesh_manifest.py` `212bcf44...`). The
four-edge identity mesh is byte-identical to the V-a calibration mesh
(SHA-256 `80966c7db44dabc49ac7bb068bbaee0e0828e413b115cee8c066c886866b6108`,
the mesh physics-11 ran on; rotate-z `97bfb753...` identical too, CanonicalBuildId
`41e4f243...` shared with the calibration roots because the cache key does not
encode recipe options), so physics-11 binds to the production root.

| case | elements (tets + prisms + pyramids) | nodes = H1 p1 (H1 p4) | tubes / layers / thickness min-P50-max (nm) / max neighbour ratio | caps min SJ / max cond | tets min SJ / max cond | prisms / pyramids max cond | corners | protected (measure / vertex) | closure (points) | diagonal bands | trace-apex shells 0-25 / 25-50 / 50-100 nm (cells; longest edge P50 nm) | junction first layer cut P50 / tets P50 (nm) | stage s / GiB (build; publication; placement) | audits s / GiB (identity; rotate-z) | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| four-edge-9d2cb9bbb3fe | **1,869,209** = 1,683,863 + 172,107 + 13,239 (= V-a; +7.3% vs 083874fee's 1,742,434; 0.54x the 34B production 3,471,480) | 414,758 (24,827,674) | 8 / 1,471 / 10.85 - 49.75 - 50.0 / 1.93 | 12: 0.0305 / 67.6 | 0.0296 / 118.1 | 668.6 / 8.74 | 3.13 / 3.44 / 3.47 / 3.54 | 0 / 9.9e-16 | 2.5e-13 (653,577; owners 3000/3100/5001/6001) | 0 | 54 apexes: 2,802 / 2,260 / 5,949 cells at 6.4 / 32.1 / 41.3 (083874fee 4.9 / 32.6 / 42.7) | 21.6 / 25.5 | 115.6 / 4.25; 29.7 / 3.34; 96.9-111.8 / 3.8 | 284.6 / 2.88; 346.1 / 2.85 | 991 / 2.95 |
| ten-edge-6791f1c84123 | **2,673,691** = 2,383,387 + 269,568 + 20,736 (+8.5% vs 083874fee's 2,465,185; 0.75x the 34B production 3,570,533) | 611,846 (35,774,961) | 20 / 2,304 / 15.9 - 49.82 - 50.0 / 1.83 | 36: 0.0274 / 82.8 | 0.0200 / 173.9 | 589.5 / 8.74 | 3.34-3.80 (10 corners) | 0 / 1.3e-15 | 1.1e-13 (1,031,142; 10 owner labels 3100/3101, 5001/5002, 5101/5102, 6001/6002, 6101/6102) | 0 | 204 apexes (180 at ratio 1.0: more basis triangles fall below FarSize at 0.5): 3,155 / 3,788 / 10,360 cells at 8.9 / 46.1 / 56.2 | 21.7 / 25.3 | 194.3 / 5.04; 36.4 / 4.34; 148.7-151.8 / 5.26 | 344.0 / 3.85; 496.1 / 3.83 | 1,396 / 4.35 |

Every physical gate passed at its production value (positive orientation and
condition <= 1000 for every element type; tetrahedra SJ >= 0.01; corner aspect <=
4.0; protected <= 1e-8; closure <= 1e-12; <= 4M elements; 1800 s / 8 GiB per
stage). The ten-edge tetrahedra minimum scaled Jacobian 0.0200 sits at 2x the gate
(083874fee: 0.0217); the trace-apex shell statistic of the ten-edge is not
comparable with the ratio-1.0 build because the apex set grew from 180 to 204
(the wider basis triangles now requesting sizes below FarSize join the set). Roots
(binaries kept; the physics-11 mesh is the four-edge identity):
`/tmp/coupon-gmsh-only-four-edge-9d2cb9bbb3fe-21cf534ff-20260919-050014`
(identity.msh `80966c7db44dabc49ac7bb068bbaee0e0828e413b115cee8c066c886866b6108`,
rotate-z `97bfb753c7148c421ebb1fe7be89811e7d88ebbd6309d8a3a24a42e4cedd94bc`,
CanonicalBuildId `41e4f24391f99725716443a2eea578c00f419f0471c3af075c9d806990a55b1e`),
`/tmp/coupon-gmsh-only-ten-edge-6791f1c84123-21cf534ff-20260919-050014`
(identity.msh `082cf7d432c983ad3a489a1c3d4d9c13c80a9c117194d65c03c342b54899a1ee`,
rotate-z `d2424ec8c3aeb10f205c90633112512d181bc5184e3b679123d0e346e842fbbb`,
CanonicalBuildId `2c47ef1224af3660a96fbaf800752df4fc928e463ae72f0aa52f0b36d01f37c7`).
The superseded 083874fee roots (four `647b2079...`, the physics-10 mesh; ten
`eb9986a4...`) and the V-b calibration root (`7db360be...`, physics-12) keep their
records; their mesh binaries were removed once these roots verified. The V-a
calibration root keeps its records and its (identical) binaries were removed in
favour of the production root above.

### Evidence (real gallery cases three-edge-419576fdab24 and two-edge-8dd4bc70f183, commit a22b471c1 tools, 2026-09-19)

The two graded_v2 gallery inputs with completed references (`06`, `10`; above)
were built through `run_gmsh_only_case.py` under the production recipe
(TraceBasisSizeRatio 0.5, every option and gate at its production value, no
case-specific constant) and verified identity + rotate-z-0.63 with empty failure
lists (`Passed true`, `Failures []`, `TransformComparisonFailures []`,
`CanonicalReuseFailures []`; the verification reports bind the manifest SHA-256
`1c99aced64fe7235dd1a87b81a2d8d558ebf43f11fd09b19689deba9625e24ca`, the manifest
carrying both cases). The semantic contracts were derived by
`derive_semantic_contract.py` from the copied inputs and the census of a
production-option probe build of the same inputs (`Derivation.BuildCensusSHA256`);
the probe meshes were byte-identical to the production builds' `gmsh-build.msh`
(same inputs, options and corners; the contract enters the build through its
SemanticCorners only).

| case | elements (tets + prisms + pyramids) | nodes = H1 p1 (H1 p4) | tubes / layers / thickness min-P50-max (nm) / max neighbour ratio | caps min SJ / max cond | tets min SJ / max cond | prisms / pyramids max cond | corners | protected (measure / vertex) | closure (points; owners) | diagonal bands | interface areas um^2 | trace-apex shells 0-25 / 25-50 / 50-100 nm (cells; longest edge P50 nm) | trace rule achieved/requested P50 / max (narrow triangles) | junction first layer cut / tets P50 (nm) | stage s / GiB (build; publication; placement) | audits s / GiB | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| three-edge-419576fdab24 (input 06, 135 traces) | **1,693,002** = 1,507,026 + 172,692 + 13,284 | 390,906 (22,737,354) | 10 / 1,476 / 15.65 - 49.85 - 50.0 / 1.94 | 14: 0.0412 / 27.7 | 0.0259 / 71.5 | 589.5 / 8.74 | 3.56 / 2.71 / 3.44 / 3.41 / 3.56 (5 corners, two on the box at the island loop) | 0 / 1.3e-15 | 2.7e-13 (637,170; 3100 / 5001 / 6001) | 0 | 1: 694.6, 3100: 125.8, 5001: 100.0, 6001: 103.6 | 97 apexes: 4,084 / 3,106 / 8,439 at 5.8 / 35.6 / 45.8 | 0.98 / 1.15 (132) | 21.7 / 25.2 | 108.0 / 4.32; 26.4 / 3.15; 102.5-107.5 / 3.7 | 223 / 2.67; 289 / 2.69 | 775 / 2.70 |
| two-edge-8dd4bc70f183 (input 10, 77 traces) | **516,662** = 420,902 + 88,920 + 6,840 | 133,763 (7,764,541) | 12 / 760 / 15.82 - 49.43 - 50.0 / 1.82 | 20: 0.0483 / 19.6 | 0.0276 / 70.5 | 589.5 / 8.74 | 3.70 / 3.49 / 3.39 / 3.04 / 3.04 / 3.64 (6 corners, two at the strip ends inside the box) | 0 / 4.9e-16 | 4.7e-14 (277,578; 3100 / 5001 / 5002 / 6001 / 6002) | 0 | 1: 206.0, 3100: 37.9, 5001: 4.0, 5002: 4.0, 6001: 4.9, 6002: 4.9 | 75 apexes: 3,147 / 2,026 / 5,797 at 4.3 / 36.9 / 45.5 | 0.86 / 1.08 (112) | 21.7 / 25.6 | 45.5 / 2.37; 14.3 / 1.39; 34.6-35.5 / 1.6 | 65 / 1.00; 85 / 0.92 | 265 / 0.92 |

Both are far under the caps (4M elements, 1800 s, 8 GiB) and pass every physical
gate at its production value; no gate, option or rule was changed or tuned for
them. The two-edge case exercises two conductors in one slot with finite metal
strips (their far ends at x = -5 and x = 4 are Physical metal edges without a
signature row: they get tubes and graded corners like every Physical side) and
the three-edge case a metal island loop whose Physical vertices lie on the coupon
box. Recorded risk (both): the etch footprint is the producer default (no
`retained-etch.csv` exists for these inputs); the `InterfaceAreas` above are the
mesh's per-label areas to be checked against the reference meshes by the physics
preflight's area invariants before any accuracy comparison. Roots (binaries kept
for the physics preflight):
`/tmp/coupon-gmsh-only-three-edge-419576fdab24-a22b471c1-20260919-053602`
(identity.msh `f55095d80c42c7e2709ff22f82e4ea4eef57cc908f723e7dc3ef29aba2e51728`,
rotate-z `9e78545ab4582e2cc5c06d6d2538d8c142ed4d8df7694dbf24545eef001a3af3`,
CanonicalBuildId `6d6a2e8e4e75b4cbbdc63eb6fb869418ae1323784bb5a47bf994c415eadead66`),
`/tmp/coupon-gmsh-only-two-edge-8dd4bc70f183-a22b471c1-20260919-053602`
(identity.msh `15e9386589150c5dbfc3428063a6bee186898f4c5ecb9000af315aa0e9a74233`,
rotate-z `36a640bfef0382ae775f1a71e3e5b0517f05d53c79290337a9b60e30cc00c29c`,
CanonicalBuildId `61d7a9564b5c46e7fc36103bb43d0c6651e8a3522ead0fa81dee86d0d71bb672`).
The probe roots keep their censuses (bound by the contracts' Derivation); their
mesh binaries were removed. Validation of the decision-42 state (commits f6efe6367
.. 2565fedf0): `python3 -m unittest discover -s . -p "test_*.py"` ran 250 tests,
OK (33 skipped); `run_general_mesh_suite.py --preflight-only` passes for the
production manifest (14 cases), the MA calibration manifest (6) and the sizing
calibration manifest (2); `refreeze_manifest_tools.py --check` current.

### Evidence (decision 43 production roots: three-edge 06, four-edge, ten-edge, two-edge 10; commits b6378f173 / 72185ce89, 2026-09-19)

All four production cases were rebuilt through `run_gmsh_only_case.py` under the
unchanged manifest options and gates (identity + rotate-z-0.63, audits,
verification `Passed true`, `Failures []`, `TransformComparisonFailures []`,
`CanonicalReuseFailures []`; the verification reports bind the manifest SHA-256
`452fd4c32c0174d4...`). The three-edge and four-edge roots were launched from
b6378f173 and verified under 72185ce89 (the fix-forward of the curve-spacing
validator bound; the mesher is identical in both commits) - their first audit
attempts, rejected by the too-strict bound, are kept under
`superseded-audit-attempts/` / `failed-audits-b6378f173-validator/` in the roots.

| case | elements (tets + prisms + pyramids); delta vs decision 42 | nodes (H1 p4) | tets min SJ / max cond | prisms / pyramids max cond | caps min SJ / max cond | corners | protected (measure / vertex) | closure | diagonal | trace-apex shells 0-25 / 25-50 / 50-100 nm (cells; longest edge P50 nm) | trace rule achieved/requested P50 / max (narrow) | curves (graded) / min node spacing nm | build s / GiB | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| three-edge-419576fdab24 (06) | **1,845,349** = 1,659,373 + 172,692 + 13,284; **+152,347 tets (+9.0%)**, prisms / pyramids unchanged | 428,366 (24,551,802) | 0.0201 / 79.2 | 589.4 / 8.74 | 14: 0.0201 / 79.2 | 3.45 / 3.42 / 3.72 / 3.10 / 3.62 | 0 / 3.5e-18 | 3.9e-13 | 0 | 103 apexes: 4,844 / 3,706 / 9,685 at 6.7 / 33.6 / 43.2 | 0.97 / 1.06 (136) | 24 (24) / 0.23 (metal ridge parts at the corners; band curves 10.5, junction 17.7) | 117 / 4.6 | 926 / 2.9 |
| four-edge-9d2cb9bbb3fe | **1,916,486** = 1,731,140 + 172,107 + 13,239; **+47,277 tets (+2.5%)**, prisms / pyramids unchanged | 426,520 (25,100,892) | 0.0213 / 97.7 | 588.1 / 8.72 | 12: 0.0213 / 97.7 | 3.10 / 3.78 / 3.71 / 3.43 | 0 / 9.9e-16 | 1.5e-13 | 0 | 58 apexes: 3,244 / 2,488 / 6,342 at 5.6 / 29.4 / 40.1 | 0.92 / 1.15 (78) | 36 (28) / 0.23 (band 11.0, junction 12.9) | 123 / 4.6 | 1,040 / 2.9 |
| ten-edge-6791f1c84123 | **3,498,453** = 3,208,149 + 269,568 + 20,736; **+824,762 tets (+30.8%)**, prisms / pyramids unchanged; **87.5% of the 4M cap (501,547 headroom; the cap fails closed)** | 814,867 (45,596,995) | 0.0358 / 89.3 | 588.8 / 8.73 | 36: 0.0369 / 20.1 | 2.87 - 3.77 (10) | 0 / 1.3e-15 | 7.7e-14 (10 owner labels) | 0 | 224 apexes: 5,877 / 6,694 / 16,648 at 8.2 / 33.5 / 43.4 | 0.98 / 1.09 (302) | 48 (42) / 0.23 (band 12.5, junction 21.6) | 311 / **7.66 of 8** | 1,607 of 1,800 / 5.3 (concurrent with another audit) |
| two-edge-8dd4bc70f183 (10) | **521,676** = 425,916 + 88,920 + 6,840; **+5,014 tets (+1.2%)**, prisms / pyramids unchanged | 134,927 (7,823,347) | 0.0202 / 99.3 | 586.3 / 8.69 | 20: 0.0319 / 55.7 | 3.30 / 3.54 / 3.73 / 3.25 / 3.74 / 3.38 | 0 / 0 | 1.2e-13 | 0 | 75 apexes: 3,185 / 2,063 / 6,075 at 4.7 / 37.2 / 45.5 | 0.86 / 1.08 (112) | 14 (12) / 0.23 (junction 21.6, band 25.0) | 48 / 2.4 | 261 / 0.9 |

Roots and digests (binaries kept for the physics preflight; the decision-42 roots'
mesh binaries were deleted, their censuses and reports kept):
`/tmp/coupon-gmsh-only-three-edge-419576fdab24-b6378f173-20260919-162000`
(identity.msh `95837ed5aca05bd52a571b27ba935014a4b40f16d1bd1c9cb9437a3bb5205a70`,
rotate-z `bc352df932df76609159d18800ddc610c189cbb5e87a44bd057535b776a75cc3`,
CanonicalBuildId `1a2a44781a334b67ff57f3f50eb744f6235af4262c429c71d857eef4472c46ae`);
`/tmp/coupon-gmsh-only-four-edge-9d2cb9bbb3fe-b6378f173-20260919-162000`
(identity.msh `1d536c448885cda2d8af7441c63a9b7b591c18b0b127366c1958ba65030a9ba5`,
rotate-z `5443f2d8357159abda8d6425814598f859966ada2f66947317bb2b6360c4d72e`,
CanonicalBuildId `6257fa3ce62225b11e1cb49419afc136206f9a83f3737c443236f16829b2dfd2`);
`/tmp/coupon-gmsh-only-ten-edge-6791f1c84123-72185ce89-20260919-170018`
(identity.msh `8a871f222da9370da81df84585e9a3c772a62aa67188d5934a239f697a2ca75a`,
rotate-z `b1c15ecce0971df92c88543f32a240351b8bfdbfa392420c0120b88d8ecd729d`,
CanonicalBuildId `bb8bf518c44b298a7ede5d3ccdd3712690c77336547ec2e6d938d66f6bf7a4b2`);
`/tmp/coupon-gmsh-only-two-edge-8dd4bc70f183-72185ce89-20260919-170018`
(identity.msh `5d01204e3396744cdbcf6f53fd1ff85d9f8c218eec1f0b9e11fc6a58aa469b04`,
rotate-z `ebb790ad1b89a2b264780fc8f9d8e450bb47071fe9f6f89b797b1647c504139a`,
CanonicalBuildId `5b42e72caf710458bbf13dae80845d802e7a41b7694b2ae5f8d68ffb3a6014a9`).

**The four-edge production root no longer equals the physics-11 mesh**
(`80966c7d...`, 1,869,209 elements): +47,277 tetrahedra (+2.5%), located where the
four-edge basis has triangles whose altitude is below the shortest edge (its
narrowest 21.7 nm slivers have a 15.3 nm altitude) and along the band curves the
trace rule now grades (the junction slivers 74 / 10: top / bottom-face triangles
within 0.3 um 412 -> 658 / 650 at mean edge 14.7 -> 10.5 nm, side-face triangles
unchanged in size, tetrahedra within 0.1 um 575 / 533 -> 760 / 747 at longest edge
P50 29.5 / 31.6 -> 25.4 / 26.1 nm; box-edge curve nodes 10.8-28 -> 7.6-21 nm). The
physics-11 evidence of decision 42 stays bound to its own mesh SHA; the production
root is the decision-43 mesh.

Sliver-support statistics of case 06 before / after (the diagnosis measure, within
0.3 um of the apex on each face and 0.1 um in the volume): sources 25 / 26 / 133 /
134 - box-edge curve node spacing 24.8-37 -> 5.7-11.5 nm; side-face (x = -8) cut
triangles 186-195 at mean edge 30-31 nm -> 305-315 at 25 nm; bottom / top-face
(needle) triangles 104-124 at 35-36 nm -> 732-868 at 6.9-7.0 nm; tetrahedra 115-136
at longest edge P50 52-57 nm -> 606-948 at 17.2-18.3 nm. The twin 52 / 53 (no
needle) is nearly unchanged (431 -> 456 triangles at 26.9 -> 26.0 nm; tetrahedra
210 -> 257 at 51 -> 49 nm). Every band / junction curve of case 06 is graded (the
trace rule at the trench rows prescribes 24.9 nm < NormalSize along the junction
lines, and the basis rows lower the footprint band curves locally); the un-tubed
metal ridge parts inside the corner clearance grade to CornerSize (0.23-0.25 nm) as
before. Physics on the new 06 root (gallery-physics-06b) closes the class.

Validation of the decision-43 state (b6378f173, 72185ce89): `python3 -m unittest
discover -s . -p "test_*.py"` ran 251 tests, OK (33 skipped);
`run_general_mesh_suite.py --preflight-only` passes for the production manifest
(14 cases), the MA calibration manifest (6) and the sizing calibration manifest (2);
`refreeze_manifest_tools.py --check` current. Disk hygiene: the 19 superseded mesh
binaries (1.70 GB: the four decision-42 roots and the three probe meshes) are listed
in `/tmp/coupon-sliver-fix-20260919/deleted-binaries.txt`.

### Coupon-scale size bound (decision 45(b), 2026-09-19): TangentialSize = min(--lc-tangent, FarSize)

The recipe fixes the tangential spacing in absolute units (`--lc-tangent` 0.05 um:
the tube extrusion spacing and the metal ridge grid) while `FarSize =
FarSizeOverRadius x Radius` scales with the coupon. On the Radius-0.5 synthetic
coupons FarSize is 0.04 um and the mesher used to fail closed ("tangential mesh
size must lie between fine and far sizes"). The rule, implemented in
`mesh_spatial_coupon.jl` (`SIZE_BOUND_RULE`) and nowhere else: **every coarsening
size prescription that exceeds the coupon-scale FarSize is bounded by it** -
`TangentialSize = min(--lc-tangent, FarSize)`. Rationale: FarSize is the coarsest
size the coupon admits (the size prescribed at its matching surface), so the
along-edge spacing of the metal-edge tubes can never legitimately exceed it; the
bound is dimensionless (it acts exactly when `Radius < --lc-tangent /
FarSizeOverRadius` = 0.625 um at the production values) and is the identity for
every production coupon (Radius 2, FarSize 0.16). The resolution sizes - NormalSize
(0.25 x MetalThickness) and EdgeSize = CornerSize - are never bounded: a fine size
above FarSize is a contradictory recipe and still fails closed. No new parameter and
no case constant. Recorded: census `SizeBounds` (`Rule`, `FarSize`,
`RequestedTangentialSize` = `--lc-tangent`, `TangentialSize`,
`TangentialSizeBoundByFarSize`); bound: `mesh_stage_contract.validate_size_bounds`
(called by `validate_gmsh_build_census`) requires the record to follow the command
exactly and the tube record's TangentialSize to equal the bound value; the recipe
binding still requires the command to execute `--lc-tangent` at 0.05 (the request is
the recipe, the bound is the coupon's). Manifest `ProductionRecipe.Parameters.
TangentialSize` states the rule. Tests: `GmshOnlyPipelineTest.
test_gmsh_build_census_contract_negatives` (record present, flag, request, bound value,
tube TangentialSize; a bounded command accepted); the fixture producer records the same
rule.

### Pre-build element estimate: the headroom gate (milestone review of f6efe6367..bcb9e18af, P1; 2026-09-19)

The ten-edge root sits at 87.5% of the 4M cap and the cap fails closed only after a
5-minute build. `estimate_build_cost.py` now estimates a production case's element
count BEFORE anything is built, from the frozen inputs alone, and the case fails
closed when the estimate exceeds the unchanged `Gates.MaximumElements`: in
`run_general_mesh_suite.py --preflight-only` (per-case `BuildCostEstimate` record;
the preflight fails) and in `run_gmsh_only_case.py` before its first stage (the root
receives `build-cost-estimate.json` and nothing else; a tight-cap dry run on case 05
refuses in under a second). The estimate is the size-field integral of the recipe's
own laws, N = (1 / TetrahedraPerCubicSize) x sum over components of the integral of
dV / h^3, plus the recipe's prisms and pyramids: far field (box volume / FarSize^3),
tube band (NormalSize at the tube surface growing with FarGrowth, per tube length on
the dielectric side), corner balls and tube caps (graded shells from CornerSize, then
the corner exterior law), junction / band curves on a recorded proxy length (the
metal-loop perimeter: the cut-surface junction lines are a producer outcome), and the
trace basis - every basis triangle whose requested size s = TraceBasisSizeRatio x
minimum altitude lies below FarSize contributes the half-space Steiner shell integral
of s + FarGrowth x d over its area and perimeter, so a needle (small altitude, long
perimeter) is charged what it costs. Tubes: Layers = length / TangentialSize (bounded
by FarSize), prisms per layer Sectors + 2 Sectors (Rings - 1) = 117, pyramids Sectors
= 9 at the production process. `TetrahedraPerCubicSize` = 0.216 is the one
dimensionless model constant (tetrahedra Gmsh realizes per h^3 of the prescribed
field; an equilateral tetrahedron has volume 0.1179 h^3): calibrated as the MINIMUM
integral / actual ratio over the seven verified production roots below, so the
estimate never falls below any calibrated build (fail-closed direction); recorded with
the per-root ratios in `ProductionRecipe.BuildCostEstimate` and validated by
`validate_manifest`. It changes no mesh and no gate. Estimate vs actual (the estimate
of case 05 was evaluated on its inputs; its build had already run when the gate was
added, so the comparison is post hoc for 05 as for the others):

| case | estimated elements | actual | ratio | tets est / actual (ratio) | prisms est / actual; pyramids | integral shares far / tube / corners / junction / trace | estimate / cap |
|---|---|---|---|---|---|---|---|
| 05 | 1,486,171 | 1,485,070 | 1.001 | 1,324,891 / 1,322,026 (1.002) | 149,760 / 151,398; 11,520 / 11,646 | 0.57 / 0.13 / 0.00 / 0.02 / 0.29 | 0.372 |
| 06 | 1,886,489 | 1,845,349 | 1.022 | 1,705,049 / 1,659,373 (1.028) | 168,480 / 172,692; 12,960 / 13,284 | 0.62 / 0.11 / 0.02 / 0.02 / 0.23 | 0.472 |
| four-edge | 1,966,879 | 1,916,486 | 1.026 | 1,785,439 / 1,731,140 (1.031) | 168,480 / 172,107; 12,960 / 13,239 | 0.67 / 0.10 / 0.02 / 0.02 / 0.19 | 0.492 |
| ten-edge | 3,561,695 | 3,498,453 | 1.018 | 3,278,951 / 3,208,149 (1.022) | 262,548 / 269,568; 20,196 / 20,736 | 0.52 / 0.09 / 0.04 / 0.02 / 0.33 | 0.890 |
| 10 | 560,105 | 521,676 | 1.074 | 469,385 / 425,916 (1.102) | 84,240 / 88,920; 6,480 / 6,840 | 0.45 / 0.20 / 0.13 / 0.03 / 0.19 | 0.140 |
| three-edge-current-calibration | 1,496,561 | 1,496,998 | 1.000 | 1,315,121 / 1,311,526 (1.003) | 168,480 / 172,224; 12,960 / 13,248 | 0.80 / 0.14 / 0.03 / 0.03 / 0.00 | 0.374 |
| concave-multislot | 1,448,136 | 1,454,056 | 0.996 | 1,423,944 / 1,423,816 (1.000) | 22,176 / 27,720; 2,016 / 2,520 | 0.95 / 0.01 / 0.04 / 0.00 / 0.00 | 0.362 |

The needle-heavy ten-edge basis is where the gate matters: its trace-basis integral
(232k, min requested size 2.7 nm) is 63% of its far field and the estimate 3.56M is
0.89 of the cap; evaluated with the pre-decision-44 shortest-edge measure instead of
the altitude, the same model predicts the decision-44 tetrahedra growth as +26.9% /
+7.1% / +2.7% on ten-edge / 06 / four-edge against the measured +34.6% / +10.1% /
+2.8%, and the decision-42 ten-edge root (2,383,387 tets) at 1.084x. The tetrahedra
ratios lie in 1.000-1.102 and the element ratios in 0.996-1.074: the prisms are
under-counted where the layers follow the size field (20% on the 1 um concave coupon
whose short tubes are mostly inside corner balls, 2-3% elsewhere), which the
tetrahedra over-count covers except on that coupon (-0.4%). Calibration manifests
(labeled experiments) are not gated by the estimate.

### Synthetic fixtures regenerated on the Radius-2 process (supervisor decision 46) and the coupon box rule for CAD-subdivided edges (decision 47), 2026-09-19

Decision 46: the five Radius-12.5 fixtures are synthetic generality tests, not gallery
references, so a self-contradictory fixture has no evidentiary value and is
regenerated as a consistent one: each is re-bound to the existing Radius-2 sharp
process (`testdata/four-edge-9d2cb9bbb3fe/process.toml`, SHA `883b8b34...`: identical to
`generality-sharp-process.toml` except Radius 12.5 -> 2.0 and the license header - the
box every loop was authored on), its contract re-derived by
`derive_semantic_contract.py` with a probe census, and its source hashes re-frozen as
fixture version 2 (`FixtureVersion: 2`); the version-1 bindings (process SHA
`703ec878...`, the pre-45(b) contract SHAs and the provisional ones of 869465f32) are
kept in the manifest's `RetiredFixtures` record with the contradiction that retired
them - never edited silently. The contradiction, per fixture (box the loop was
authored on vs the Radius-12.5 box the version-1 binding produced; every Continuation
vertex lies on the former): one-edge-straight [-4, 4] x [-8, 8] vs [-25, 25] x
[-14.5, 14.5] (loop vertices (-4, +-8), (0, 8)); two-edge-transition [-8, 6] x [-8, 8]
vs [-23, 27] x [-25, 25]; two-edge-multislot [-8, 8] x [-4, 8] vs [-22.5, 27.5] x
[-25, 25]; six-edge-cluster [-9.9167, 8.9167] x [-4.5, 4.5] vs [-26, 25] x [-25.5,
25.5]; one-edge-cad-subdivided shares one-edge-straight's loop but its two collinear
half rows (P = (0, -1) / (0, 1), S in [-1, 1]) never reached Radius from their own
reference points, so even at Radius 2 the pre-47 mesher shrank its box to [-4, 4] x
[-4, 4] and the loop lay outside it (`RetiredFixtures` carries these numbers; probe
logs under `/tmp/coupon-matrix-case05-20260919/probe-root-<case>-v1` /
`probe-root-<case>`). Version-2 probe builds (production options): one-edge-straight
796,062 elements (areas 1: 452.8, 3100: 64.8, 5001: 64.0, 6001: 65.6),
two-edge-transition 1,496,338 (694.6 / 125.8 / 100.0 / 103.6 - the 06 plan at Radius
2), two-edge-multislot 1,117,402 (7 labels, 3000 / 3001 absent), six-edge-cluster
1,329,382 (11 labels), one-edge-cad-subdivided 796,062 with `gmsh-build.msh`
byte-identical to one-edge-straight's (`1d19bb1c...`).

Decision 47 (the fifth fixture): a coupon box that depends on how the source CAD
subdivides a straight edge is a generality defect the `cad-subdivision-sensitivity`
comparison exists to catch, so it is fixed in the mesher, gated by an inertness
probe. Rule (`mesh_spatial_coupon.jl` `register_edge_chains!` / `extended_interval`,
`COUPON_BOX_RULE`): collinear rows of one metal edge that touch end to end (same slot,
conductor, plane, tangent and gap, no vertex arm) form a chain, and the 2 x Radius
extension is decided on the chain's union interval about its midpoint (extended at
both ends iff the union reaches Radius from the midpoint, i.e. total length >= 2 R -
the single symmetric row rule applied to the union); every other row keeps the
single-row rule exactly. Rows that merely share a line but do not touch (two edges
separated by a slot) never chain. Recorded in the census `CouponBox` (rule, radius,
bounds, `EdgeChains` with rows / union / union length, `ChainedRows`,
`ExtendedChains`) and bound by `mesh_stage_contract.validate_coupon_box` (rule text,
consistent bounds, chains with >= 2 distinct rows, recomputed counts; negatives in
`GmshOnlyPipelineTest`); Julia test `test_edge_chains.jl` (a subdivided edge reproduces
the unsubdivided box and extension, a chain shorter than 2 R is not extended, gapped /
other-conductor / vertex-arm rows do not chain); the estimator's box replica carries
the same rule (`estimate_build_cost.edge_chains`, tested). Inertness probe
(`decision47-box-probe.txt`, the Python box replica cross-checked equal to the
mesher's recorded box on the five trace-basis cases): the coupon box of every
registered case of the three manifests is UNCHANGED by the rule - the only touching
chain in the inventory is one-edge-cad-subdivided's (union length 4.0); ten-edge's four
and six-edge's two collinear row pairs are separated by slots (e.g. ten-edge (-6.5,
-0.6) S [0, 2] and (-6.5, 0.4) S [-2, 0]: a 1 um gap) and do not chain. Per case
(before .. after, x-y):

```
suite             one-edge-straight                  R=2.0   before [-4.0, -8.0]..[4.0, 8.0] after [-4.0, -8.0]..[4.0, 8.0] chains=[] UNCHANGED
suite             one-edge-cad-subdivided            R=12.5  before [-25.0, -14.5]..[25.0, 14.5] after [-25.0, -14.5]..[25.0, 14.5] chains=[(2, np.float64(4.0))] UNCHANGED
suite             two-edge-transition                R=2.0   before [-8.0, -8.0]..[6.0, 8.0] after [-8.0, -8.0]..[6.0, 8.0] chains=[] UNCHANGED
suite             two-edge-multislot                 R=2.0   before [-8.0, -4.0]..[8.0, 8.0] after [-8.0, -4.0]..[8.0, 8.0] chains=[] UNCHANGED
suite             three-edge-current-calibration     R=2.0   before [-8.0, -8.0]..[6.0, 8.0] after [-8.0, -8.0]..[6.0, 8.0] chains=[] UNCHANGED
suite             four-edge-9d2cb9bbb3fe             R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
suite             six-edge-cluster                   R=2.0   before [-9.9167, -4.5]..[8.9167, 4.5] after [-9.9167, -4.5]..[8.9167, 4.5] chains=[] UNCHANGED
suite             ten-edge-6791f1c84123              R=2.0   before [-11.75, -8.6]..[9.8333, 8.4] after [-11.75, -8.6]..[9.8333, 8.4] chains=[] UNCHANGED
suite             three-edge-419576fdab24            R=2.0   before [-8.0, -8.0]..[6.0, 8.0] after [-8.0, -8.0]..[6.0, 8.0] chains=[] UNCHANGED
suite             two-edge-8dd4bc70f183              R=2.0   before [-5.0, -2.0]..[4.0, 3.0] after [-5.0, -2.0]..[4.0, 3.0] chains=[] UNCHANGED
suite             two-edge-3f8992613e95              R=2.0   before [-4.0, -8.0]..[6.0, 8.0] after [-4.0, -8.0]..[6.0, 8.0] chains=[] UNCHANGED
suite             concave-multislot                  R=0.5   before [-1.5, -1.5]..[2.7, 2.5] after [-1.5, -1.5]..[2.7, 2.5] chains=[] UNCHANGED
suite             hole                               R=0.5   before [-2.3, -2.3]..[2.3, 2.3] after [-2.3, -2.3]..[2.3, 2.3] chains=[] UNCHANGED
suite             rounded-strip                      R=0.5   before [-2.1, -1.2]..[2.1, 1.2] after [-2.1, -1.2]..[2.1, 1.2] chains=[] UNCHANGED
suite             opposed-layers                     R=0.5   before [-2.1, -2.1]..[2.1, 2.1] after [-2.1, -2.1]..[2.1, 2.1] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-v1              R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-v2              R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-edge-layer-4nm  R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-el1nm-50        R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-el4c            R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-ma    four-edge-calib-ma-el1c            R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-sizingfour-edge-calib-sizing-tbr-0.5     R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
calibration-sizingfour-edge-calib-sizing-lct-0.025   R=2.0   before [-6.0, -8.0]..[10.0, 8.0] after [-6.0, -8.0]..[10.0, 8.0] chains=[] UNCHANGED
```

The unchanged real roots therefore stay valid: identical box and identical mesher
behaviour on unchained rows; verified by rebuilding the two-edge 10 case under the new
mesher and comparing its `gmsh-build.msh` SHA with the decision-43 root (evidence
below).

### Evidence (decisions 48 / 49: hole, opposed-layers and the concave-multislot relabel rebuild; commits 42f997c62 / 97ed84952 / bdfe87890, 2026-09-20)

Built through `run_gmsh_only_case.py` under the unchanged production recipe (every
option and gate at its production value; audits at 8 GiB), identity + rotate-z-0.63,
verified `Passed true` with empty `Failures` / `TransformComparisonFailures`; every
root passed the headroom gate first and its `build-summary.json` records
`Status built`. Roots `/tmp/coupon-gmsh-only-<case>-<commit>-*`:

| case | elements = tets + prisms + pyramids | H1 DOFs | estimate (ratio) | tets SJ / cond | prisms / pyramids cond | corners (n) | protected err / vertex | closure (points) | diag | interface areas (um^2) | build s / GiB | verification s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hole (v2; bdfe87890, identity `6ab22454...`, CanonicalBuildId `a16919cb...`) | **1,852,882** = 1,795,858 + 52,272 + 4,752 | 345,684 | 1,854,486 (1.001) | 0.0245 / 45.9 | 465.2 / 13.7 | 3.42 - 3.76 (8) | 0 / 0 | 1.6e-13 (364,020) | 0 | 1: 62.744, 3100: 19.224, 5001: 2.200, 6001: 2.904 | 114 / 4.00 | 600 |
| opposed-layers (v2; 97ed84952, identity `cc6e04f7...`, CanonicalBuildId `cc4e4c44...`) | **2,289,512** = 2,240,552 + 44,064 + 4,896 | 414,330 | 2,252,546 (0.984) | 0.0130 / 142.2 | 464.4 / 27.4 | 2.77 - 3.70 (8) | 0 / 2.2e-16 | 1.8e-13 (551,952) | 0 | 1: 62.832, 3000 / 3001: 1.680 / 1.680, 3100 / 3101: 15.248 / 15.248, 5001 / 5101: 0.961 / 0.959, 6001 / 6101: 1.202 / 1.198 | 185 / 4.33 | 871 |
| concave-multislot (v2; 97ed84952, identity `3c86b72a...`, CanonicalBuildId `d255af28...`) | **1,454,056** = 1,423,816 + 27,720 + 2,520 | 265,608 | 1,448,136 (0.996) | 0.0207 / 139.4 | 469.7 / 13.9 | 3.23 - 3.73 (6) | 0 / 2.2e-16 | 2.6e-13 (252,990) | 0 | 1: 51.804, 3000 / 3001: 0.4533 / 0.1767, 3100 / 3101: 8.950 / 6.830, 5001 / 5101: 0.325 / 0.245, 6001 / 6101: 0.517 / 0.405 | 113 / 3.47 | - |

- hole: the hole's four sides carry tubes pointing into it (16 tubes = 2 x 8 sides of
  both loops); the metal areas are exact (5001 = 1.6^2 - 0.6^2 = 2.2; 6001 = 2.2 +
  8.8 x 0.08 = 2.904; 3100 = trench floor 18.96 + walls 8.8 x 0.03 = 19.224); the
  producer-default collars etch the whole Radius-0.5 coupon (no 3000 plane). The
  first hole root (42f997c62 + the extension-1 mesher) and the rebuild under the final
  mesher of bdfe87890 give byte-identical `gmsh-build.msh` (`49258267...`) and
  `identity.msh` (`6ab22454...`): the extension-2 labeling change is inert on it. The
  earlier root's binaries were deleted (evidence kept).
- opposed-layers: upward layer at z = 0 and downward layer at z = 0.6, two slots on
  both; box z in [-0.52, 1.12] (per-sign padding); vacuum gap 0.48 um against the
  facing reach 2 x 0.03975; the un-etched strips of both planes (3000 / 3001 = 2 x 4.2
  x 0.2 = 1.68 each) and the trenches (floor 15.0 + walls 0.08 + collar walls 0.168 =
  15.248 each) are exact; the multi-slot ownership postprocessor certified both facing
  layers (closure 1.8e-13 over 551,952 points). The tetrahedral minimum scaled
  Jacobian 0.0130 is the smallest margin of the matrix over the 0.01 gate.
- concave-multislot: element count identical to the decision-45(b) build (1,454,056;
  the mesh differs only by the un-etched labels 3000 / 3001 = 0.4533 / 0.1767 um^2 the
  fixed labeling restores); the superseded 869465f32 root's binaries were deleted.
- Probe roots (`/tmp/coupon-scope-20260920/probe-root-*`, the mirror-covariance strips
  and the mislabeled opposed-layers probe) keep their censuses / logs; their meshes are
  listed in `/tmp/coupon-scope-20260920/deleted-binaries.txt`.

**Matrix status after decision 48 (15 cases): built and verified 14** (the twelve of
decisions 46 / 47 with concave-multislot rebuilt at version 2, plus hole and
opposed-layers); **unsupported class 1**: rounded-strip (`TopRounding`, recorded by the
preflight as `UnsupportedClass`, not a failure; the rounded-edge extension stays out of
scope per decision 48(3): a process-model / reference decision first). Validation sweep
at bdfe87890: `unittest discover` **263 OK (skipped=33)**; preflights production 15
cases / 14 passed / 1 unsupported (`rounded-strip: TopRounding`; max estimate 0.89 of
the cap), calibration-ma 6 / 6, calibration-sizing 2 / 2.

### Evidence (decisions 46 / 47: the five regenerated fixtures and the two-edge 10 inertness rebuild; commit 8b0057dbd, 2026-09-19)

All six built through `run_gmsh_only_case.py` under the production recipe (options and
gates unchanged; audits at 8 GiB), identity + rotate-z-0.63, verified `Passed true`
with empty failure lists; the verification reports bind the manifest SHA-256
`7120d18c...` (commit 8b0057dbd). Every root passed the headroom gate first
(`build-cost-estimate.json`).

| case (fixture v2) | elements (tets + prisms + pyramids) | nodes | estimate (ratio) | tets min SJ / max cond | prisms / pyramids max cond | tubes / layers | corners | protected | closure (points; owners) | interface areas um^2 | box; chains | build s / GiB | verification s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| one-edge-straight | **796,062** = 714,666 + 75,582 + 5,814 | 176,133 | 794,118 (0.998) | 0.0396 / 23.8 | 589.5 / 8.74 | 2 / 646 | 3.67 | 0 / 3.5e-18 | 3.1e-13 (287,781; 3100 / 5001 / 6001) | 1: 452.8, 3100: 64.8, 5001: 64.0, 6001: 65.6 | [-4, 4] x [-8, 8]; 0 | 44.4 / 2.5 | 331 |
| one-edge-cad-subdivided | **796,062** (identical) | 176,133 | 794,118 (0.998) | 0.0396 / 23.8 | 589.5 / 8.74 | 2 / 646 | 3.67 | 0 / 3.5e-18 | 3.1e-13 (287,781) | identical | [-4, 4] x [-8, 8]; **1 chain (rows 1-2, union 4.0, extended)** | 44.5 / 2.5 | 333 |
| two-edge-transition | **1,496,338** = 1,310,866 + 172,224 + 13,248 | 339,728 | 1,496,561 (1.000) | 0.0257 / 63.1 | 589.5 / 8.74 | 10 / 1,472 | 3.05 - 3.50 (5) | 0 / 3.5e-18 | 5.9e-13 (593,391) | 1: 694.6, 3100: 125.8, 5001: 100.0, 6001: 103.6 | [-8, 6] x [-8, 8]; 0 | 76.8 / 3.9 | 675 |
| two-edge-multislot | **1,117,402** = 1,029,328 + 81,783 + 6,291 | 236,752 | 1,136,205 (1.017) | 0.0351 / 100.8 | 589.5 / 8.74 | 6 / 699 | 3.30 / 3.52 / 3.20 | 0 / 3.5e-18 | 3.3e-13 (377,358; 3100 / 3101 / 5001 / 5101 / 6001 / 6101) | 1: 613.6, 3100: 29.53, 3101: 57.32, 5001: 31.84, 5101: 74.16, 6001: 32.62, 6101: 75.08 | [-8, 8] x [-4, 8]; 0 | 60.9 / 3.3 | 457 |
| six-edge-cluster | **1,329,382** = 1,119,970 + 194,454 + 14,958 | 318,937 | 1,371,715 (1.032) | 0.0201 / 126.4 | 589.5 / 8.74 | 20 / 1,662 | 3.23 - 3.79 (10) | 0 / 1.4e-17 | 5.2e-13 (626,208; 10 owners) | 1: 566.85, 3100: 8.97, 3101: 90.03, 5001 / 5002: 0.254 / 0.255, 5101 / 5102: 36.00 / 36.00, 6001 / 6002: 0.351 / 0.356, 6101 / 6102: 37.90 / 37.89 | [-9.917, 8.917] x [-4.5, 4.5]; 0 | 82.7 / 3.7 | 548 |
| two-edge-8dd4bc70f183 (10, inertness rebuild) | **521,676** = 425,916 + 88,920 + 6,840 | 134,927 | 560,105 (1.074) | 0.0202 / 99.3 | 586.3 / 8.69 | 12 / 760 | 3.25 - 3.74 (6) | 0 / 0 | 1.2e-13 (279,090) | 1: 206.0, 3100: 37.9, 5001 / 5002: 4.0, 6001 / 6002: 4.9 | [-5, 4] x [-2, 3]; 0 | 43.9 / 2.4 | 238 |

The decision-47 rule is inert on the real cases as the probe predicted: the two-edge
10 rebuild under the new mesher gives `gmsh-build.msh` `10afee1f...` and identity.msh
`5d01204e3396744cdbcf6f53fd1ff85d9f8c218eec1f0b9e11fc6a58aa469b04`, byte-identical to
the decision-43 root (`...-72185ce89-20260919-170018`), whose binaries are kept while
the rebuild's were deleted (its census / reports are the evidence). The
`cad-subdivision-sensitivity` comparison is live again and trivially satisfied: the
subdivided and the straight fixture produce byte-identical identity meshes
(`ff9898a5e001f64e96f6226e25828d67455aeee69bbc6c8249e58500bf5be00f`; normalized DOF
ratio 1.0 against the 1.25 bound), the only difference being the census `CouponBox`
chain record. Roots: `/tmp/coupon-gmsh-only-<case>-8b0057dbd-2026091[6-9]-*` (identity
SHA-256s: one-edge-straight / -cad-subdivided `ff9898a5...`; two-edge-transition
`455f0dc754e27cfd7b75f91b5887c3751630dc53c6583dd04aa6d4349006d7a1`; two-edge-multislot
`45ded77f1dea16af56740d34bf88e2d80d1422bc004df798bdd9ea04c06c9eb6`; six-edge-cluster
`5f749dbaab96cc9d0af0e22337cf9ad89385af040eb1fcc24660dd65f013d40c`; CanonicalBuildIds
1fd47f70... / 722caab8... / e3f5fda2... / 7e91cc6f... / f38a201c...). The first
decision-46 attempts of one-edge-straight and two-edge-transition (roots
`...-7abfd08ff-20260919-155819-superseded-mesher-edit`) were superseded when the
decision-47 mesher edit landed during their audits (tool digest changed); their and
the probe meshes are listed in `deleted-binaries.txt`.

**Matrix status after decisions 46 / 47 (15 cases; superseded by decision 48 above: 14 built / 1 unsupported class): built and verified 12**
(four-edge, ten-edge, three-edge 06, two-edge 10, two-edge 05,
three-edge-current-calibration, concave-multislot, one-edge-straight,
one-edge-cad-subdivided, two-edge-transition, two-edge-multislot, six-edge-cluster);
**unbuildable 3**: hole, rounded-strip, opposed-layers - outside the prism-tube
recipe's stated scope (hole loop, rounded edge, downward layer), a producer feature
decision, not a repair.

### Recipe scope as a recorded statement (supervisor decision 48, 2026-09-20)

Every fail-closed guard of the prism-tube recipe is now a recorded, machine-readable
statement, so that a library run distinguishes "unsupported class" from a bug:

- **Classes.** `mesh_spatial_coupon.jl` `RECIPE_SCOPE_SUPPORTED_CLASSES` are the input
  classes the recipe builds (`ContinuationVertices`, `DeviceFootprint`, `ExteriorLoops`,
  `MultipleConductors`, `MultipleLayers`, `MultipleSlots`, `TraceBasis`);
  `RECIPE_SCOPE_GUARDS` are the classes it fails closed on, each with a stable id, a
  statement and its detection origin: visible in the frozen inputs (`HoleLoops`,
  `DownwardLayers`, `TopRounding`, `TrenchRounding`, `SlopedSidewalls`, `ThinMetal`,
  `NoTrench`) or only in a derived quantity during the build (`ShallowTrench` - the
  pyramids would reach the trench floor; `NarrowTransverseBound` - no ring fits
  min(Overetch, MetalThickness / 2, CornerIsotropyRadius); `FreeEdgeEnds` - an edge end
  neither a semantic corner nor on the box; `ShortEdges` - no tube interval remains
  after the corner clearances; `FootprintWithoutEdge` - an explicit footprint without
  the metal edge; `FootprintTopology` - a producer-default collar whose region is not one
  simple polygon, decision 54a below). `mesh_stage_contract.py` spells the same two lists
  (`RECIPE_SCOPE_SUPPORTED_CLASSES`, `RECIPE_SCOPE_GUARDS`).
- **Guard messages.** A guard fails with `ScopeGuard[<id>]: <statement>; <detail>`
  (`scope_error`), replacing the former prose messages ("Prism edge tubes support
  exterior conductor loops only", "Prism tubes require sharp vertical fabricated
  geometry", "upward process layers only", "The tube pyramids would reach the trench
  floor", "neither a semantic corner nor on the box", ...).
- **Census.** The build census records a `Scope` block: `Recipe`, `SupportedClasses`,
  `GuardedClasses`, `Guards[]` (id, origin, statement), `ExhibitedClasses` (the classes
  this input exhibits, from the loops, the signature layers and the process options)
  and `MetalLoops[]` (per plan-view loop: conductor, plane, hole flag, vertices and the
  straight sides not on the outer box). The former descriptive `Scope` string is now
  `Purpose`. `validate_gmsh_build_census` -> `validate_recipe_scope` binds the block: the
  lists equal the contract's, the exhibited classes equal the classification recomputed
  from the bound signature / boundary inputs and the command's process options (none of
  them guarded - a guarded class never reaches a census), the loop sides equal the
  recount from the bound boundary against the census `CouponBox`, and
  `PrismTubes.TubeCount = 2 x` the sides of all loops (negatives in
  `test_general_mesh_manifest.py`; Julia `test_prism_tube_build.jl`).
- **Drivers.** `run_gmsh_only_case.py` writes `build-summary.json` in the root with
  `Status` `built` / `unsupported-class` / `failed`: a mesher stop whose log carries
  `ScopeGuard[<id>]` is recorded with the id (`UNSUPPORTED_CLASS <id>` on stderr),
  distinctly from any other failure. The manifest preflight (`run_general_mesh_suite.py
  --preflight-only`) classifies every Gmsh-only case from its frozen inputs
  (`Cases[].Scope.ExhibitedClasses / UnsupportedClasses`), records a case outside the
  scope as `UnsupportedClass` with the error `unsupported class <id>` (summary
  `UnsupportedClassCases`), and does not count it as a preflight failure of the matrix:
  the case is never built or judged, so a full matrix run reports it as not passed with
  its class, not as a bug.

#### Interior conductor loops (holes) in scope (decision 48, extension 1)

`HoleLoops` moved from the guarded to the supported classes. In
`metal_edge_segments` the tube normal points away from the metal: out of an exterior
loop and into a hole (the polygon interior of a hole is dielectric: the same
`loop.hole` inversion the tetrahedral path applies in `offset_loop_points`); every
other component was loop-agnostic by inspection (corner clearance from the in-plane
angle between the two tube edges, hole corners = contract corners from the boundary's
`Physical` vertices, producer-default collars via `offset_hole_points`, junction and
band curves from the fragmented CAD, ownership by conductor and z-band). Added:

- `NarrowHoles` guard (build-detected): a hole must be wider than twice the tube reach
  `Radius + PyramidHeight + ProtectedDistance` (2 x NormalSize, the band law's protected
  distance) between any two of its non-adjacent sides (`hole_facing_width`), so the
  tubes facing each other across it keep disjoint bands; recorded as
  `PrismTubes.Section.FacingReach / FacingRule` (hole fixture: 0.6 um against 2 x 0.05975).
- The etch-footprint check `assert_etch_carries_edge` now runs over the hole sides too
  (a device footprint must carry them; negative in `test_prism_tube_build.jl`).
- Census tube rows carry `Hole`; `validate_recipe_scope` requires `TubeCount = 2 x` the
  straight sides of ALL loops and the hole flags of `MetalLoops` to agree with the
  exhibited `HoleLoops` class.
- Julia tests: hole tube normals point towards the hole centre, a hole coupon has 2 x
  sides of all loops tubes (16 for the square-with-square-hole), right-angle hole
  corners take the exterior right-angle clearance, facing width / segment distance.
- Fixture `hole` re-frozen as `FixtureVersion 2`: the version-1 contract listed 2 of the
  8 `Physical` vertices as corners (retired in `RetiredFixtures` with the
  `ScopeGuard[FreeEdgeEnds]` evidence); the contract is re-derived by
  `derive_semantic_contract.py` with the probe census (labels 1 / 3100 / 5001 / 6001:
  the Radius-0.5 coupon is fully etched under the producer-default collars, no 3000
  plane). The production build's `gmsh-build.msh` equals the probe's byte for byte.

#### Downward process layers (Nz = -1, flip-chip) in scope (decision 48, extension 2)

`DownwardLayers` moved from the guarded to the supported classes. In
`build_edge_tubes!` the tube frame follows the layer: `b = (0, 0, Nz)`, `e = n x b`
(the extrusion sense mirrors; both `along` senses were already handled), the top tube
lies on the metal top face at `plane + Nz x MetalThickness`, the bottom tube on the
plane; the sections are defined in the (n, b) frame, so "vacuum above / substrate
below" mirrors with b (`TubeFrameRule` in the census Section). Added:

- `NarrowLayerGap` guard (build-detected): the vacuum gap between the metal top faces
  of an upward layer and the downward layer above it (`layer_groups` admits no other
  two-layer configuration) must exceed twice the tube reach `Radius + PyramidHeight +
  ProtectedDistance`, the same rule as `NarrowHoles` (opposed-layers fixture: gap
  0.48 um against 2 x 0.03975).
- Per-layer-sign vertical box padding in `coupon_bounds` and
  `estimate_build_cost.coupon_box`: Overetch on the substrate side (-Nz),
  MetalThickness on the metal side (+Nz) of every row's plane (`COUPON_BOX_RULE`);
  identity for upward-only coupons, opposed-layers box z in [-0.52, 1.12].
- Census tube rows carry `Layer` (Nz); `validate_tube_layers` binds every row to the
  signature's Nz on its plane and to `Origin[3] = Plane + Layer x --metal-thickness`
  (top) / `Plane` (bottom).
- **Interface labeling by the surface's own layer (generic defect fixed; supervisor
  decision 49, 2026-09-20: layer-band selection of the interface label + the box
  tolerance, the mislabeling condition Radius < 1 um, concave-multislot FixtureVersion
  2).** The CAD
  interface classification took the nearest signature edge over ALL layers to decide
  the un-etched plane (3000 + s) vs the etched trench (3100 + s) and the metal
  surface's slot, and compared the surface's z-range with that edge's plane at the
  source tolerance 1e-7 x Radius - below the 1e-7 padding of the OCC bounding box, so
  for Radius < 1 um every un-etched plane was labeled 3100 + s (Radius-2 cases: 2e-7 >
  1e-7, unaffected). Now `surface_process_layer` selects the layer whose band
  `[plane - Nz x Overetch, plane + Nz x MetalThickness]` contains the surface (fail
  closed otherwise), its edges own the surface, and the flat-plane test uses the box
  tolerance like every other bounding-box comparison. Effect: opposed-layers records
  3000 / 3001 (1.68 um^2 each) and 3100 / 3101 (15.248 each) instead of 3100 / 3101
  (16.928 each) with an identical mesh; concave-multislot (Radius 0.5) gains 3000 /
  3001 (0.4533 / 0.1767 um^2) and is re-frozen as `FixtureVersion 2` with its version-1
  contract retired; hole (fully etched) is unchanged. Inertness on production: all five
  gallery cases (four-edge 07, ten-edge 09, three-edge 06, two-edge 10, two-edge 05)
  and the six Radius-2 fixtures bind Radius 2 um, where the source tolerance 2e-7
  exceeds the 1e-7 CAD padding, so the old rule already labeled their planes correctly
  by construction; the cheapest real case, two-edge 10 (`two-edge-8dd4bc70f183`), was
  rebuilt under the fixed mesher (bdfe87890): `gmsh-build.msh` `10afee1f...` and
  `identity.msh` `5d01204e...` byte-identical to the kept decision-43 root
  (`/tmp/coupon-gmsh-only-two-edge-8dd4bc70f183-72185ce89-20260919-170018`; the
  rebuild's binaries were deleted, its logs kept under
  `/tmp/coupon-scope-20260920/two-edge-10-inertness`). Regression guard: the Julia
  testset "un-etched plane of a Radius-0.5 coupon is labeled 3000" builds a
  single-slot L-shaped Radius-0.5 coupon whose producer-default collar leaves the
  0.63 um^2 notch un-etched and requires label 3000 with that area (and 3100 = 15.90);
  under the old rule it fails (3000 absent, 3100 = 16.53 - verified on a patched copy).
- Julia tests (`test_prism_tube_build.jl`): per-sign box padding; downward tube frames
  are the mirror of the upward ones (b, origin z, e; same intervals, sections and
  materials); the `NarrowLayerGap` guard at the threshold; mirror covariance of a
  downward-only strip coupon against the upward one through the full mesher: the
  interface areas, tube rows, prism / pyramid counts and the tube end cross-sections
  mirror to roundoff, every gate passes on both, the tube node sets mirror within one
  axial sampling step (`TangentialSize / TUBE_LAYER_SAMPLES_PER_SIZE`: the layer
  stations follow the axis field sampled from `s_start`, so reversing the extrusion
  sense moves a station by < 1.5e-3 um here - the direction dependence two oppositely
  traversed exterior edges already have) and the tetrahedra are reported, not asserted
  equal (Gmsh's Delaunay kernel is not reflection-covariant: 27,699 vs 27,732 tets on
  the probe).
- Fixture `opposed-layers` re-frozen as `FixtureVersion 2` (contract re-derived: 8
  corners, 9 labels); the multi-slot ownership postprocessor runs on two facing
  layers (first production case with two slots on different planes).

### Synthetic matrix under the Gmsh-only recipe (decision 45(b), 2026-09-19)

The 2026-09-17 record (below, under the retired recipe) found ten of the twelve
then-registered cases unbuildable. Under the Gmsh-only production recipe every one
was re-examined; contracts are now derived, never authored:

- **Contracts regenerated with `derive_semantic_contract.py`** (no hand-edited
  number): `one-edge-semantic.json`, `one-edge-subdivided-semantic.json`,
  `two-edge-transition-semantic.json`, `two-edge-multislot-semantic.json`,
  `six-edge-semantic.json`, `three-edge-semantic.json` (the six whose corners
  contradicted their boundaries) and `concave-multislot/semantic-contract.json`.
  The tool now accepts the case's bound file names (`--signature`, `--boundary`,
  `--process-library`) and, when no process library is bound (every synthetic
  fixture), takes the slot / conductor pairs from the signature and records
  `Derivation.SlotConductorSource`. Semantic corners are now exactly the boundary's
  `Physical` vertices (1 / 1 / 5 / 3 / 10 / 5 / 6), so the canonical-source
  validation ("Transformed physical boundary differs from semantic corners") passes
  for all; label families are the producer's (`etched-substrate-vacuum-slot-s`,
  `conductor-c-slot-s-ms/-ma`; the un-etched 3000 + s plane only where a build
  census shows it). `three-edge-semantic.json` and `concave-multislot` are bound to
  the gmsh-build census of a production-option probe build (`Derivation.
  BuildCensusSHA256`); the five Radius-12.5 fixtures cannot be built (below), so
  their contracts carry `Derivation.Provisional` (un-etched plane unconfirmed) and
  their required label set. Preflight passes for all 15 cases; the six-edge
  fixture's contract still drives the rigid-transform and multislot-census tests.
- **The two "seed-gate failures" (two-edge-transition, six-edge-cluster)** were
  MMG-era seed-side failures (required-region optimization gates); under the
  Gmsh-only builder those gates do not exist as a separate stage, and both cases
  fail earlier, for the same frozen-input reason as the other Radius-12.5 fixtures
  (next item). No generic mesher bug was found: the mesher's fail-closed message is
  correct for the inputs.
- **Five Radius-12.5 fixtures are unbuildable with their frozen inputs**
  (`one-edge-straight`, `one-edge-cad-subdivided`, `two-edge-transition`,
  `two-edge-multislot`, `six-edge-cluster`): gmsh-build stops at "Metal edge end
  (x, y) is neither a semantic corner nor on the box". Their plan-view boundary
  loops were authored on the Radius-2 coupon box - every `Continuation` vertex lies
  exactly on the box the signature spans at Radius 2 (one-edge: box x +-4 / y +-8,
  loop vertices (-4, +-8), (0, 8); two-edge-transition: box [-8, 6] x [-8, 8];
  two-edge-multislot: [-8, 8] x [-4, 8]; six-edge: [-9.9167, 8.9167] x [-4.5, 4.5];
  one-edge-cad-subdivided binds the same loop to a signature whose Radius-2 box is
  only +-4, so it is inconsistent at either radius) - but the cases bind
  `generality-sharp-process.toml` (Radius 12.5: boxes +-25 x +-14.5 and larger), so
  the metal loops close in the coupon interior without a Physical edge. A
  Continuation vertex means "the metal continues past the box"; a loop vertex inside
  the coupon that is not a semantic corner is a contradiction between two frozen
  inputs (process vs boundary), not a contract or mesher defect. Repairing it means
  re-binding those cases to a Radius-2 process (or re-authoring the loops): a
  frozen-input rewrite, left as a decision. Evidence: the probe roots
  `/tmp/coupon-matrix-case05-20260919/probe-root-<case>/gmsh-build.log` and the
  box / loop numbers above.
- **Three Radius-0.5 fixtures are outside the production recipe's stated scope**:
  `hole` ("Prism edge tubes support exterior conductor loops only" - its loop is a
  hole), `rounded-strip` ("Prism tubes require sharp vertical fabricated geometry" -
  TopRounding 0.005), `opposed-layers` ("Prism tubes support upward process layers
  only" - four of its eight edges have Nz = -1). Each is the prism-tube recipe's
  own fail-closed scope statement (decision 38); extending the tubes to hole loops,
  rounded edges or downward layers is producer feature work, not a repair. The size
  bound above did act on all three before they stopped (their probe logs show no
  tangential-size error), so nothing else hides behind these messages. Their
  contracts were not in the contradictory six and are unchanged.
- **Built**: `three-edge-current-calibration` (Radius 2, the 06 inputs without a
  trace basis) and `concave-multislot` (Radius 0.5, the first production build
  with `TangentialSizeBoundByFarSize true`: TangentialSize 0.05 -> 0.04), plus the
  new gallery case `two-edge-3f8992613e95` (input 05): evidence below.

### Evidence (decision 45(b): gallery case 05, three-edge-current-calibration, concave-multislot; commit 869465f32, 2026-09-19)

Built through `run_gmsh_only_case.py` under the production recipe (every option and
gate at its production value; audits at 8 GiB), identity + rotate-z-0.63, verified
`Passed true` with empty failure lists (`Failures`, `TransformComparisonFailures`,
`CanonicalReuseFailures`); the verification reports bind the manifest SHA-256
`cef6ba510fb2863be5f50d380b8fbe2dc565af8f441f70b907c3a285c7644f21` (commit
869465f32). The probe builds that fed `derive_semantic_contract.py` produced
`gmsh-build.msh` files byte-identical to the production builds' (case 05 also to a
build with the pre-size-bound mesher bcb9e18af: the bound is the identity at Radius
2), so the contract enters the build through its SemanticCorners only.

| case | elements (tets + prisms + pyramids) | nodes | tubes / layers / thickness min-P50-max (nm) / max neighbour ratio | caps min SJ / max cond | tets min SJ / max cond | prisms / pyramids max cond | corners | protected (measure / vertex) | closure (points; owners) | diagonal | interface areas um^2 | SizeBounds | build s / GiB | audits s (identity / rotate) | verification s / GiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| two-edge-3f8992613e95 (input 05, 125 traces; trace ratio 0.5, 246 basis triangles, min altitude 14.9 nm, 18 needles / 10 below FarSize) | **1,485,070** = 1,322,026 + 151,398 + 11,646 | 351,848 | 4 / 1,294 / 15.98 - 49.88 - 49.91 / 1.93 | 4: 0.0714 / 18.7 | 0.0485 / 28.0 | 588.5 / 8.72 | 3.58 / 3.39 (both on the box) | 0 / 3.5e-18 (rot. 8.5e-16) | 1.38e-13 (488,406; 3100 / 5001 / 6001) | 0 | 1: 535.4, 3100: 129.6, 5001: 32.0, 6001: 35.2 | FarSize 0.16, request 0.05 -> 0.05, bound false | 89.8 / 4.25 | 176 / 239 | 705 / 2.5 |
| three-edge-current-calibration (the 06 inputs without a trace basis) | **1,496,998** = 1,311,526 + 172,224 + 13,248 | 339,844 | 10 / 1,472 / 15.65 - 49.85 - 50.0 / 1.55 | 14: 0.0276 / 72.7 | 0.0276 / 72.7 | 589.5 / 8.74 | 3.14 / 3.27 / 3.52 / 3.80 / 3.57 | 0 / 3.5e-18 (rot. 8.9e-16) | 6.57e-13 (593,385; 3100 / 5001 / 6001) | 0 | 1: 694.6, 3100: 125.8, 5001: 100.0, 6001: 103.6 (= case 06) | FarSize 0.16, 0.05 -> 0.05, false | 79.7 / 3.8 | 172 / 239 | 694 / 2.5 |
| concave-multislot (Radius 0.5) | **1,454,056** = 1,423,816 + 27,720 + 2,520 | 265,608 | 12 / 280 / 16.56 - 29.58 - 39.83 / 1.28 | 24: 0.0572 / 22.7 | 0.0207 / 139.4 | 469.7 / 13.9 | 3.23 / 3.67 / 3.73 / 3.49 / 3.26 / 3.31 | 0 / 2.2e-16 (rot. 4.4e-16) | 2.66e-13 (252,990; 3100 / 3101 / 5001 / 5101 / 6001 / 6101) | 0 | 1: 51.80, 3100: 9.403, 3101: 7.007, 5001: 0.325, 5101: 0.245, 6001: 0.517, 6101: 0.405 | **FarSize 0.04, request 0.05 -> 0.04, bound true** (tube spacing max 39.8 nm) | 75.8 / 3.4 | 137 / 140 | 515 / 2.2 |

Case 06's mesh (1,845,349) exceeds the same inputs without a trace basis by 348,351
tetrahedra: the trace-basis rule's cost on that geometry. The concave-multislot
coupon (1 x 1 um box at FarSize 0.04 = 2 x NormalSize) is nearly uniform and
1.45M elements dense; it passes every physical gate unchanged. Roots (binaries kept):
`/tmp/coupon-gmsh-only-two-edge-3f8992613e95-869465f32-20260919-134938` (identity.msh
`9f011ed3a76e93f1245f3f893ad17dc0b924791161efc0e26b2e794fa9e51b6c`, rotate-z
`d520f9ed0ad61288c9347be39b6776bbfc2a935ee48711865e4f794b21fff52e`, CanonicalBuildId
`ac79dcb458e1104e46d9113fc116db23a18179b0c668de35e6a5b2e32b948ead`);
`/tmp/coupon-gmsh-only-three-edge-current-calibration-869465f32-20260919-134938`
(`5e48f2f7ee23f1f6f61710249e9eb833e3c3b7cd5b36bfd94011eebcc6ec37e3`,
`0b0a46232dbb78b9c48f86b7588f3a432feb24817f7a58f14526face2c900f49`,
`a03ce80a67dd8147d1b9ff2b25a1ef4723edb037dc04ca82ee42f3c4add2d65e`);
`/tmp/coupon-gmsh-only-concave-multislot-869465f32-20260919-141300`
(`7e7ca8c26c6baa701923154a4ba32dcb48a587937e05ff663950cf1e15653e9c`,
`734fca04f8e26fda8650b19d553d3de11363a27a5d20d612ff604ceb7d83e501`,
`9eff740479ddf11bb545414f40332af8aa3a69f6b925a5acf67b21952e4c8e7b`). Probe roots
(censuses kept, meshes deleted) and the unbuildable cases' logs:
`/tmp/coupon-matrix-case05-20260919/probe-root-<case>`; the 17 superseded mesh
binaries (1.24 GB: the case-05 build launched before the size-bound commit, whose
rotate-z audit failed on the changed mesher digest, and the probe meshes) are listed
in `deleted-binaries.txt` there.

**Matrix status under the Gmsh-only recipe before decisions 46 / 47 (15 cases): built and verified 7**
(four-edge, ten-edge, three-edge 06, two-edge 10, two-edge 05,
three-edge-current-calibration, concave-multislot); **unbuildable 8**: five
Radius-12.5 fixtures whose frozen boundary loops were authored on the Radius-2 box
(process vs boundary contradiction, a frozen-input rewrite decision) and three
Radius-0.5 fixtures outside the prism-tube recipe's stated scope (hole loop, rounded
edge, downward layer). Of the ten cases the 2026-09-17 record found unbuildable, two
now build (three-edge-current-calibration: its contract contradicted its labels and
corners; concave-multislot: FarSize below TangentialSize) and eight remain so for
the two reasons above, none of them a contract or size-rule question any more.

Validation of the decision-45(b) state (6e4700cec, 869465f32): `python3 -m unittest
discover -s . -p "test_*.py"` ran 252 tests, OK (33 skipped);
`run_general_mesh_suite.py --preflight-only` passes for the production manifest (15
cases), the MA calibration manifest (6) and the sizing calibration manifest (2);
`refreeze_manifest_tools.py --check` current.
After the milestone-review scope addition (71fbe4d9d: headroom gate, P2 items) the
same sweep ran 260 tests, OK (33 skipped); preflights 15 (every case with its
`BuildCostEstimate`, maximum 0.89 of the cap on ten-edge) / 6 / 2 pass;
`refreeze_manifest_tools.py --check` current. The 869465f32 roots stay bound to their
commit's manifest (`cef6ba51...`); nothing was rebuilt.
After decisions 46 / 47 (8b0057dbd) the sweep ran 261 tests, OK (33 skipped);
preflights 15 / 6 / 2 pass (maximum estimate 0.89 of the cap); `refreeze_manifest_tools.py
--check` current; the Julia `test_edge_chains.jl` testsets pass (17 assertions).

## Retired production recipe (supervisor decision 34B, 2026-09-17; legacy MMG pipeline)

Retired from production by decision 38; recorded as
`ProductionRecipe.RetiredLegacyRecipe` and kept verifiable under the calibration
manifest. The 34B production recipe was the EL4c calibration recipe, recorded in the manifest's
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
cut surface the element size must not exceed the ratio times the minimum
altitude of the basis triangle containing the point (the shortest edge until
decision 43, above) - a per-triangle rule, because the
hat of a basis vertex varies linearly over the whole incident triangle, so its
support is resolved where the hat varies only when the whole triangle is
discretized at that scale (the per-edge alternative resolves only the edges).
The seed applies it through the Gmsh size callback on top of the scalar
corner-isotropy background (`min(background, ratio x minimum altitude + slope x
distance to the triangle)` with the process-band grading slope, capped by
`lc_far`; `Mesh.MeshSizeMin` follows the smallest requested size, the only
sub-`lc_fine` request any field can make) and records `TraceBasisSizing` in the
census (ratio, rule, frame, input digests, box, counts, basis edges below the
far size, minimum requested size, mesh-frame triangles, mesh size minimum,
slope). The metric stage binds the same files, requires the census record to
match (digests, ratio, triangles), caps the volume metric isotropically by
`min(FarSize, ratio x minimum altitude + FarGrowth x distance)` so the existing
far/grading law is kept away from narrow hats (MMG's hmin = NormalSize still
floors the volume metric; the frozen cut triangles keep the seed's sizes), and
records `TraceBasisSizing` in the recipe with the statistics: unique basis
edges, edges below the far size, triangles below it, minimum requested size,
and `CutSurfaceSize` (min/median/max longest edge of the seed's cut triangles;
per narrow basis triangle the largest extent of the cut triangles centred in it
across its minimum altitude over the requested size - reported, not
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
(Decision 61c, 2026-09-21: the cap deviation is RETIRED - `Gates.MaximumElements` is
6,000,000 in every manifest and the EL1c case declares no cap; the record is kept under
`Calibration.RetiredGateDeviations.MaximumElements`.)
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

The matrix contains 15 cases and 30 required case/variant entries. All 15 cases
have hash-frozen local source contracts. A bounded read-only assessment on
`soca-green` copied only the approved source-contract files from campaign inputs
`07` and `09` (2026-09-14) and, under decision 42 (2026-09-19, `soca-green-job`),
from the graded_v2 gallery inputs `06` and `10`, which have completed references
(`reference-library/models/007-...-419576fdab24` and `010-...-8dd4bc70f183`), and
under decision 45(b) from input `05` (`two-edge-3f8992613e95`:
`spatialedgecluster_edgecount-2_3f8992613e95`, one 2 um wide metal strip of one
conductor spanning the coupon, two signature rows, two Physical vertices on the box,
125 trace files, Order-4 reference on `production-meshes/05-fabricated/coupon.msh`
`591e7183...`; no `retained-etch.csv`, producer-default footprint recorded as the same
risk); it
copied no mesh, field, solution, response matrix, or source bank (the per-source
`traces/` directories are referenced by count only). Campaign input `07` maps to
`spatialedgecluster_edgecount-4_9d2cb9bbb3fe` and has exactly four signature
rows. Input `09` maps to `spatialedgecluster_edgecount-10_6791f1c84123` and has
exactly ten. Input `06` maps to `spatialedgecluster_edgecount-3_419576fdab24`
(`three-edge-419576fdab24`: three signature rows, one slot / one conductor, a metal
island loop in the plan view whose two Physical vertices sit on the coupon box;
135 trace files) and input `10` to `spatialedgecluster_edgecount-2_8dd4bc70f183`
(`two-edge-8dd4bc70f183`: two signature rows, one slot / two conductors - two
finite metal strips whose far ends are Physical metal edges without a signature
row; 77 trace files). The process-library model names, local/remote SHA-256
values (remote `sha256sum` = local for every copied file, recorded in each
`provenance.json` with the copied `basis-points.csv`, `zero-trace.csv` and
`spatial_fabricated.json`), trace mesh references, slots, conductors, topology,
and target-repository Apache-2.0 license were reviewed. Neither `06` nor `10`
carries a `retained-etch.csv`, so both declare `EtchFootprint:
"producer-default"` - recorded as a risk: the producer-default collars are a
producer outcome, not a bound device footprint, and the physics preflight must
check the interface-area invariants against the reference mesh before any
accuracy comparison.

Semantic contracts are derived, never authored: `derive_semantic_contract.py
SOURCE_DIR OUTPUT [--build-census CENSUS]` derives label families and
slot/conductor coverage from `Models[0].Edges` (checked against the signature),
semantic corners from the plan-view vertices classified `Physical`, and
`FeatureTopology` from the finite signature segments. Whether a slot's un-etched
plane (3000 + slot) exists is a producer outcome of the etch footprint and the
coupon box, so the tool takes the label set from the `InterfaceAreas` of a
gmsh-build census of the same inputs (a production-option probe built with the
provisional contract the tool writes without a census), requires every census
label to belong to a derived family and every non-optional family to be present,
and records the census digest and labels under `Derivation`. It reproduces the
frozen ten-edge contract exactly and the four-edge contract up to the two older
substrate-vacuum role strings (`test_derive_semantic_contract.py`). The `06` /
`10` probes recorded labels 1 / 3100 / 5001 / 6001 and 1 / 3100 / 5001 / 5002 /
6001 / 6002: no un-etched plane under the producer-default footprint. The workflow is
two-pass by design: (1) derive without a census (PROVISIONAL: the required label set
only), (2) a census-only probe build of the same inputs under the production options
with that provisional contract (`run_gmsh_only_case.py --census-only`: headroom gate,
source validation, gmsh-build with its census, canonical publication - not the
placements: their ownership check judges the mesh against the contract and the
provisional one lacks exactly the un-etched plane labels the census confirms, as the
device four-edge coupon under the producer-default footprint showed), then derive
again with `--build-census PROBE/build-census.json`. The probe census is not validated by the stage contract
(only its `InterfaceAreas` labels are read, and only labels inside the derived
families are accepted); the production build that follows binds the final contract
and validates its own census, and since the contract enters the build through its
SemanticCorners only, the probe's `gmsh-build.msh` equals the production one (verified
byte-identical on 05 / 06 / 10 / three-edge-current-calibration). A process library
with more than one model, another topology than `SpatialEdgeCluster` (e.g. an Arms
model) or edges without `InterfaceSlot` / `Conductor` is a contract error with that
message, not a lookup failure.

## Two commands (supervisor decision 48, 2026-09-20)

The library is two commands of `coupon_library.py`: `build` (the mesh path:
registration, job-pool build, `library-build.json`) and `qualify` (the physics path:
the per-run physics scripts of `coupon-accuracy-assessment-20260913` moved into
`qualify/` as parametrized modules, plan / estimate / submit under the 40-job
accounting, fetch / validate / hygiene, machine-readable class gates,
`PendingQualification` without a reference, `library-qualification.json` and the
process-library entries).

```sh
python3 coupon_library.py build \
  --register <CASE_ID>=<SOURCE_DIR> --footprint producer-default \
  --inventory-status RepositoryAssessmentFixture \
  [--mesh-recipe examples/cpw3d_surface/spatial_coupon/testdata/generality-mesh-recipe.json] \
  [--case ID ...] [--jobs 2] [--root /tmp/coupon-library-build-<commit>-<ts>]
```

1. **Registration** (`register_case.py`, one call per `--register`): a source
   directory (`mesh-signature.csv`, `plan-view-boundary.csv`, `plan-view-mask.csv`,
   `process.toml`; optionally `process-library.json` with the trace basis
   `basis-contract.json` / `trace-vertices.csv` / `trace-triangles.csv` - all four or
   none -, `provenance.json`, `retained-etch.csv`; its own `mesh-recipe.json` or
   `--mesh-recipe`) becomes a manifest case: SHA-256 of every source file, the recipe
   scope classes of the inputs (a guarded class stops the registration as
   `unsupported-class` with the guard id, no probe built), the two-pass contract
   derivation orchestrated automatically (provisional contract -> census-only probe
   build of a staging copy under the production recipe, headroom gate included ->
   `derive_semantic_contract.py --build-census`), the production `Variants` (identity
   and rotate-z) / `TransformComparison` / `SignatureColumns` shared by every existing
   case (fail closed when they differ), `InventoryStatus`, `Features` (default: the
   exhibited scope classes), `FixtureVersion` and a `Provenance` statement (commit,
   footprint declaration, probe root), then `refreeze_manifest_tools.py`. The etch
   footprint is an explicit declaration - `--footprint bound` freezes the directory's
   `retained-etch.csv`, `--footprint producer-default` requires its absence - and the
   command fails closed without it. Idempotent by content: recorded source digests
   equal to the directory's reuse the case (nothing written); any changed source
   becomes the next `FixtureVersion` and the previous entry's changed bindings go to
   `RetiredFixtures.Entries` (`Reason`, `Evidence`), never a silent edit. The outcome
   is `WORK/register-case.json` (`Status` registered / reused / unsupported-class /
   failed, `StoppedBy`, `Scope`, `SourceSHA256`, `ContractSHA256`, `ProbeRoot`).
2. **Build** (`run_gmsh_only_matrix.py`): the selected cases (default: every case,
   the registered ones included) run as a pool of `--jobs` drivers (default 2), each
   the unchanged `run_gmsh_only_case.py` under the manifest bounds for every stage
   (`--audit-memory-gib` = `Gates.MaximumRSSGiB`, so audits and verification share
   the stage bounds - 3600 s / 12 GiB since decision 61c, 1800 s / 8 GiB before):
   headroom gate -> canonical DAG -> rigid
   placements -> audits -> `verify_canonical_case_entries.py`. A case fails closed on
   its own and the others continue; the exit status is nonzero unless every case
   passed.
3. **Record** `ROOT/library-build.json`: per case `Scope` (exhibited / unsupported
   classes), `Status` (built / unsupported-class / failed) and `StoppedBy` (the exact
   `ScopeGuard` id, `HeadroomGate` MaximumElements, the stopped `Stage` with its
   `StopReason`, or `Verification` with the failure list), `CanonicalBuildId`,
   identity / rotate-z mesh SHA-256 and paths, elements by type, `H1` DOFs at
   `--h1-order` (default 4), `Estimate` (estimated vs actual elements,
   `EstimateOverActual`, `EstimateOverCap`), every bounded stage's wall seconds / peak
   GiB / limits, the verification verdict and `HeadroomFlags` (any measure at or above
   0.9 of its bound: elements or estimate vs the cap, a stage's seconds or peak RSS vs
   its limit - the margin rule the throughput plan lacked); library totals: cases
   attempted / built / passed / unsupported / failed, `FlaggedCases`, wall clock,
   jobs, bounds, commit and manifest digest. Nothing in the record is measured by the
   command itself: every number is read from the per-case root's stage reports,
   census, estimate, audits and verification report.

4. **`build --device` (supervisor decision 52; `device_coupons.py`).** A device layout
   is the only input: `--device <Palace config> --palace <executable>` runs the
   discovery closure (`examples/cpw2d/discover_surface_response_requirements.py`:
   geometry preflights against the config's process seed - the version-3 `Fabrication`
   metadata of `Solver.*ResponseCorrection.Library` -, never a mesh or a solve), routes
   every requirement with the planner (`prepare_surface_response_coupons.plan_from_manifest`)
   and turns every `SpatialCoupon` requirement into a source directory:
   `generate_spatial_response.py --basis-only` with the planner's canonical plan-view
   boundary / mask regularization and the seed's fabrication writes the signature
   files, the trace basis (`basis-contract.json` with the geometry report,
   `FrameFitResidual` and the digests of every source trace, `trace-vertices.csv`,
   `trace-triangles.csv`, `basis-points.csv`, `zero-trace.csv`, `conductor-N.csv`) and
   the model's `process-library.json`; `process.toml` from the fabrication; a
   `provenance.json` naming the device config, the seed, the closure manifest, the
   requirement, the generator command, the ring size (the planner's default 16 - the
   basis every gallery case was produced with) and the cap triangulation. **The device
   basis triangulates the two matching-box caps with Delaunay flips by default**
   (`--cap-triangulation delaunay`, supervisor decision 57: the ear-clipped caps
   produced the needle triangles that drove 20-40% of the mesh cost and the only
   remaining device fail-closed cases, decision 54b below; recorded in
   `basis-contract.json` `CapTriangulation`, the device record's `TraceBasis` and the
   provenance, so every device case id / hat set changed with the default - they are
   new cases; `--cap-triangulation ear-clipping` selects the gallery producer's caps
   explicitly and the gallery reference cases keep their producer's basis for
   comparability). The directory is named by the content
   hash of its bound source files (`spatial-<edges>-edge-<hash12>`), so the same device
   geometry maps to the same case and `register_case.py` reuses it by content; every
   directory is registered (footprint `producer-default` - the device path binds no
   `retained-etch.csv` -, `InventoryStatus DeviceDerived`, the mesh recipe every
   trace-basis case of the manifest binds or `--mesh-recipe`). Corner and straight-edge
   requirements are recorded `OutOfScope` with their builder (their own families), never
   dropped. `--build-limit N` builds the N smallest registered cases by the pre-build
   estimate and records the others `registered-unbuilt` with their estimates
   (`Library.BuildLimit.Ranked`). On the transmon example (`examples/transmon/
   transmon_surface_coarse.json` with the benchmark process seed) the closure yields six
   spatial coupons - the four-edge `9d2cb9bbb3fe`, three-edge `419576fdab24`, two-edge
   `3f8992613e95` and ten-edge `6791f1c84123` gallery models among them, whose
   `mesh-signature.csv` / `plan-view-boundary.csv` / `process.toml` the device path
   reproduces BYTE FOR BYTE from the device config (the `plan-view-mask.csv` facet
   tessellation follows the device mesh: byte-identical for the four-edge and
   three-edge, the same footprint area per conductor to 1e-14 for the two-edge and
   ten-edge) - plus three out-of-scope families (`test_device_coupons.py`; the checked-in seed's MS
   permittivity 11.45 differs from the config's 11.47 and Palace refuses the mismatch,
   so the test binds the config's interface layers into a copy of the seed). The device
   trace basis is the producer's own (`build_matching_surface`: 120 vertices on the
   four-edge box against the gallery's 80 retained ones), so device coupons are new
   cases with their own sources, not the gallery references' basis (decision 44's open
   design item stands).

Nothing case-specific is hard-coded: the shared case fields, the recipe, the gates,
the bounds and the scope vocabulary come from the manifest and the stage contract.
Tests: `test_register_case.py` (a temporary copy of a repository source directory:
hash reuse, changed-source versioning with the retired binding, footprint / recipe /
inventory fail-closed, unsupported class from the inputs and from the probe),
`test_run_gmsh_only_matrix.py` (the record schema, stop attribution and headroom
flags from synthetic roots; with Julia, `coupon-library build` end to end: the
smallest trace-basis gallery case registered as a temporary copy through the real
probe - its derived contract equals the frozen one - and built next to
`one-edge-straight` as a pool of two, root kept under `/tmp/coupon-library-e2e-*`).

### `coupon-library qualify`

```sh
python3 coupon_library.py qualify \
  --build-record /tmp/coupon-library-build-<commit>-<ts>/library-build.json \
  --reference <graded_v2 campaign dir or none> \
  --remote soca-green-job:/data/home/simlap/coupon_accuracy_assessment_20260913 \
  [--orders p5] --controls p3,p5 --max-jobs 40 [--job-policy speed|frugal|fixed [--fixed-jobs N]] \
  [--frozen-binary-sha256 HEX] [--reducer-block-size N] \
  [--case ID ...] [--stage-prefix NAME] [--control-source I ...] [--root DIR] [--dry-run] [--resume]
```

The frozen executable defaults to the manifest's `ProductionRecipe.PhysicsRun.FrozenExecutable`
(`170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b`, the streaming one-pass Gram
build of decision 62(4); `PreviousSHA256 b28f089a…` with its provenance; decision 63), else
`build_plan.DEFAULT_FROZEN_BINARY_SHA256` (the same digest); `--frozen-binary-sha256` overrides
and the origin is recorded per coupon (`FrozenExecutable {SHA256, Origin, Rule}`). The reducer
block size defaults to `PhysicsRun.ReducerBlockSize` (48; `--reducer-block-size` overrides).

Every constant the physics runs edited per copy is an argument or a record field:
the mesh path and SHA-256 (the build record's identity variant, re-hashed before
anything is written), the source count and the trace files (every
`PrescribedPotential` entry of the reference config, digests pinned), the
ZeroTrace set (the basis contract's `ZeroTraceIndices`), the control sources (by
geometric class, or `--control-source`), the remote host / root (`--remote`), the
stage prefix (`--stage-prefix`, default the case id), the orders, `Solver.Order` /
`Linear.Tol` / materials / interfaces (byte-for-byte the reference config), the
frozen executable hash (`--frozen-binary-sha256`; the executable is
`<root>/palace-archive-estimate-<sha>.bin`), the cluster facts
(`qualify/cluster-profile.json`: queue, project, instance, 192 ranks, 6 h walltime,
20,700 s runner deadline, modules, `mpiexec_bound.sh`, the 40-job user cap) and the
measured cost rates (`qualify/cost-model.json`: the physics-11 V-a worker / reducer /
local-edge rates with the V-a entity counts, whose closed-form H1 must reproduce the
measured counts on load).

1. **Run inputs from the case itself (supervisor decision 52).** The Palace config, the
   source traces and the zero-trace set are derived from the manifest case's own frozen
   sources by `qualify/case_inputs.py`: `process-library.json` (substrate permittivity,
   interface layer thickness / permittivity, matching radius, the model's edges - slot,
   conductor, process normal placed in the mesh frame by the bound process frame - and
   its `Interfaces[].Coupon` index -> type map), the trace basis (`basis-contract.json`,
   `trace-vertices.csv`, `trace-triangles.csv`) regenerated as the producer wrote it
   (`generate_spatial_response.write_surface_trace`: `basis-NNNN.csv` = the hat of basis
   vertex NNNN, `conductor-N.csv` = the lift of every conductor but the first,
   `zero-trace.csv`; the coordinates are within the contract's `FrameFitResidual` of the
   producer's files - the canonical-frame round trip -, V and triangle columns
   identical, digests recorded under `Inputs.Sources`), the identity mesh's
   `$PhysicalNames` (the attribute candidates are filtered by them exactly as the
   producer's `make_config` did on the reference mesh; every attribute the config names
   must exist - `AttributeCheck`, before any submission) and the recipe's
   `ProductionRecipe.PhysicsRun` (`Order` 4, `LinearTol` 1e-10 with their calibration
   provenance: the four-edge reference and every recorded qualification ran them; bound
   by `general_mesh_manifest.validate_physics_run`, used by every case). The five
   gallery cases' derived configs equal the configs their graded_v2 references ran apart
   from `Model.Mesh`, `Problem.Output`, the DataFile directory and - for the p5 / Tol 1e-8
   references (gallery 10, ten-edge) - `Solver.Order` / `Linear.Tol` only
   (`test_qualify_dry_run.test_run_config_derived_from_the_case_equals_every_gallery_reference`).
   A case whose signature has a downward layer (`Nz = -1`) or more than one layer stops
   with `StoppedBy ScopeGuard` (`DownwardLayers` / `MultipleLayers`,
   `locate_sources.check_layers`: the z-level role assignment covers one upward layer)
   until the roles are assigned per layer band.
   **Reference by content.** `--reference DIR` has the layout of the physics runs'
   `reference/` trees: `inputs-<key>/` (`basis-contract.json`, the producer's
   `spatial_fabricated.json`) and `case-<key>-fabricated/` (`worker.json` = the config
   the reference ran, preferred over the producer's;
   `reducer/{domain,surface}-response-matrix.csv`). A coupon is bound to `inputs-<key>`
   whose `basis-contract.json` digest equals the manifest's `BasisContract` digest
   (`qualify/reference_campaign.py`); the reference's config must equal the derived one
   apart from the path fields and Order / Tol (`StoppedBy Reference` with the
   differences otherwise); no match = `StoppedBy Reference` (skipped: pass
   `--reference none` to run on the case's own inputs); inputs without reducer matrices
   = the coupon runs and is `PendingQualification`. `--reference none` (mandatory
   spelling: `--reference` is required) runs every coupon on its own inputs: the
   p-sequence controls alone are evaluated - `PendingQualification` when they pass,
   `Failed` when one fails, never `Passed`.
2. **Sources and controls.** `locate_sources.py` (box from the trace vertices, z levels
   from the apex heights, metal loops and junctions from the plan-view boundary, the
   3000 / 3100 adjacency from a bound `retained-etch.csv`) and `classify_sources.py`
   (ZeroTrace; junction rings; junction columns; narrow hats next to a junction;
   near-junction hats; isolated narrow hats; box 3D corners; wide hats bottom / top,
   metal-top ring, substrate / trench rings; conductor terminals = `TerminalAttributes`
   sources). The 8 controls (`--control-count`) are one source per class in that
   priority order, cycling, lowest index first (`choose_controls`), unless
   `--control-source` names them (the recorded campaigns' supervisor-specified sets).
3. **Configs and plan.** The main orders of a coupon are the recipe's `PhysicsRun`
   order, then `--orders`, then the reference's own `Solver.Order` when it differs
   (gallery case 10: reference p5 -> main stages p4 and p5, the recorded gallery-10
   layout; the recipe order
   stays the library order: cost coupon, local-edge stage, p-sequence main; recorded
   per coupon as `Orders`). `build_configs.py` derives worker / reducer at every main
   order on all sources, at every `--controls` order (highest first) on the controls,
   and the ordinary-path local-edge `config.json` at the main order on the controls
   (`SaveLocalEdgeEnergy` true); only `Model.Mesh`, `Problem.Output`, the trace
   directory, the source subset and `Solver.Order` differ from the reference config.
   `estimate_stages.py` scales the cost model by the exact H1 ratio
   (`mixed_mesh.h1_dofs_from_counts` on the build record's `H1.EntityCounts`) at 1 /
   1.5 / 2x the measured PCG counts; a coupon whose 2x total with the 35% + 300 s
   preflight margin exceeds the walltime is not one job; whose Palace peak exceeds
   0.75 of the node is `StoppedBy Estimate` before any plan. **Job policy - the
   per-coupon source split (supervisor decision 61b, user decision 60(2)):**
   `qualify/job_split.py` partitions the coupon's main-order source set into N
   contiguous blocks run as N independent worker jobs (the same mesh and configs apart
   from the `PrescribedPotential` subset, `worker-block<k>.json` per block with its own
   `Problem.Output`), every block archiving into the stage's ONE archive directory
   (Palace archives one file per source, rank and field - `source-NNNNNN-rank-NNNNNN-
   {V,D}.bin`, header-checked -, so the union of the N archives is that directory), and
   the main-order reducer runs once, in its own job, on the union after every worker
   job completed (the reduction is linear in the archive union: the matrices equal a
   single job's to roundoff - acceptance below). The p3 / p5 controls and the
   local-edge stage stay in the first worker job (`ControlsJob worker-1`), whose block
   is shortened by their estimated share so every worker job ends together; when they
   leave no room for a block the first job carries them alone (`separate`, recorded).
   Policy `{Mode: speed | frugal | fixed, MaxJobs, WalltimeSeconds, FixedJobs}`, every
   job estimated from the cost model at 2x the measured PCG counts with the 35% +
   300 s preflight margin: `frugal` = the smallest N whose every job fits the
   walltime; `speed` = the N minimizing the estimated critical path (the longest
   worker job, then the reducer job), the smallest such N on a tie; `fixed` = N given.
   N <= `--max-jobs` (the run's concurrency), the user job cap and the source count; a
   coupon fails closed only when even the maximal split does not fit (`StoppedBy
   Estimate` with the candidates table). N = 1 is the single job of the recorded
   campaigns, byte-identical (`main/plan.json`, `job.pbs`); a split coupon has
   `main/jobs/<worker-k | reducer>/{plan.json, job.pbs}` (the runner's status and logs
   per job directory). `--job-policy` / `--fixed-jobs` select it; the manifest's
   `ProductionRecipe.PhysicsRun.JobPolicy` is the recorded default (`frugal` until the
   user chooses speed vs frugality); the policy used, the candidates (N, longest job,
   critical path, node seconds), the blocks and per-job estimates are recorded per
   coupon (`JobPolicy`, `Split`, `Jobs`). `build_plan.py`: pinned SHA-256 of
   the mesh, every config and every trace; `CapSeconds` = 2 x the stage's 2x-PCG
   estimate rounded up to 300 s and bounded by the deadline, `MinimumSeconds` = the
   1x estimate rounded up (a block worker at its block's source count); `job.pbs` from
   the cluster profile; `run_stages.py` (the unchanged bounded runner, executable /
   hash / MPI wrapper read from the plan).
4. **Submission and results** (`qualify/remote.py`; not under `--dry-run`): up to
   `--max-jobs` of the run's jobs are queued / running at once (the library run's
   concurrency; the 40-job user cap is checked at every `qsub`): a coupon's worker
   jobs, then - once every worker job's `status.json` is complete and the archive
   union holds sources x ranks potential files (`archive-union.json`, fail closed) -
   its reducer job; a single-job coupon is one job. Per coupon: rsync of mesh / traces /
   `main/` to `<root>/<run>/<case>/`, `qsub` after a read-only `qstat` count of the
   user's jobs against the cap (`submission.json`, `submission-<job>.json`); every
   active job is polled read-only once per interval (job state and the runner's
   `status.json`; a job is in the queue while `qstat` reports Q / R / E / H / W / T / S / B -
   `remote.IN_QUEUE_STATES` -, a held job (H: the SOCA dispatcher holds a job whose
   compute-node stack failed, `error_message CF:ROLLBACK_COMPLETE:retry=N`, and releases
   it itself at `retry_eligible_after`) stays active with its reason logged); a worker
   job that left the queue incomplete stops the coupon (its
   other jobs finish on their own, never qdel'd, recorded `JobsLeftRunning`); a coupon
   whose last job left the queue is
   fetched while the others run - rsync of `main/` without the archives, `sha256sum`
   of every fetched CSV against the remote (`result-csv-sha256.json`),
   `run_graded_library_case.validate_matrix` on every reducer matrix (complete,
   symmetric, nonnegative), then `du` + `rm -rf` of the response archives
   (`remote-archive-deletion.json`) - and analyzed, and the next pending coupon takes
   the freed slot. Any stage not `complete`, a PCG non-convergence, a digest mismatch
   or an invalid matrix is a recorded stop. The library's `CriticalPathSeconds` is
   measured from the first submission to the last fetch; a coupon's own
   `Cost.CriticalPathSeconds` from its first submission to its fetch (the reducer
   job's queue wait included), its `Cost.Jobs` every job's estimate, actual seconds and
   node-h, its `Cost.JobNodeHours` the sum over its jobs (the merged main stage reads
   the blocks' per-source timings in block order: `summarize_cost.merge_split_statuses`);
   the library totals keep their meaning (node-h summed over every job, jobs counted at
   every `qsub`, `Splits` per coupon). `--resume` on the same
   `--root` adopts the job ids a previous driver recorded (`<case>/submission*.json`;
   the re-derived plans must be byte-identical - both sides read in path order -, else a
   recorded `Resume` stop) and
   monitors / fetches / analyzes from there, so a lost driver (VPN drop, a stalled
   process) never causes a second submission; the monitor wait is sliced against
   the wall clock (a single 90 s sleep of the idle driver was observed not to return
   on macOS during the acceptance run). A poll whose ssh round trip fails (no
   `POLL_MARKER` back: login host unreachable, banner timeout) is a recorded
   `Monitor.TransportFailures` count that keeps the job active and spends one poll of
   the budget - never "the job left the queue" -, and a fetch whose rsync fails is a
   recorded `Fetch` stop (`--resume` fetches later), not a crash of the driver (the
   2026-09-21 device run lost its login host for 30+ min with two jobs running).
5. **Qualification.** `compare_matrices.py` (main vs reference, controls vs reference,
   main vs the higher control, the lower control vs main), `classify_sources.py`
   class statistics, `ma_ms_offsets.py` (distributions, reference-p_MA-weighted view,
   strongest-20), `p_sequence.py` (d_low / d_high / r / Aitken limit at the controls),
   `key_sources.py`, `summarize_cost.py` (per-source PCG and seconds, node-h = wall x
   nodes, the full-coupon extrapolation) and `gates.py` on the frozen
   `qualify/qualification-gates.json` (its SHA-256 in every record; free view = the
   reference's sources minus the ZeroTrace knots and zero-energy sources; anchor
   "vs p<reference order>"):

   | gate | statement |
   |---|---|
   | E | every source of the four wide classes within 1% (the all-free count and the worst source reported) |
   | p_MA | free signed median within 1%; reference-p_MA-weighted mean within 1%; every one of the 20 strongest reference-MA sources (all free sources when fewer) within 2%; every free source within 5% |
   | p_MS | free |median| within 1%; every free source within 5% |
   | p_SA | >= 2/3 of the free sources within 2% and >= 90% within 5% (the EL4c level 40 / 56 of 60) |
   | p-sequence controls | every control's step to the higher order d_high = (p_high - p_main)/|p_high| within 1% for E and 5% for p_MA / p_MS / p_SA |

   Anchor rule: when a main stage was solved at the reference order, the gates
   evaluate that same-order comparison and the other main orders are recorded as
   `Informational`; otherwise the first main order is gated against the named anchor
   (`GatedOrder`, `GatedComparison`). A participation whose interface the reference
   config does not postprocess (no `Postprocessing.Dielectric` entry of that type;
   case 10 declares MA and MS only) is `NotApplicable` for its gate and its p-sequence
   observable - recorded with the declared interfaces, never a failure.
   The interface index -> type map of the response matrices is read from the run
   config's `Postprocessing.Dielectric` entries (never `{1 MA, 2 MS, 3 SA}`); a
   participation p_X sums every postprocessed interface of type X (the ten-edge
   postprocesses MA / MS per slot). The class thresholds (narrow hat width 0.1 um,
   junction reach 0.6 um) are the frozen gate table's `SourceClasses` block, covered by
   its digest. Verdict `Passed` only when every gate passes; `Failed` otherwise, a
   failed p-sequence control included; `PendingQualification` only when there are no
   reference matrices AND every p-sequence control passed (never `Passed`). On the stored CSVs: physics-11 passes with
   E 60 / 60, p_SA 35 / 44 / 58, p_MS 51 / 60 / 60, p_MA 32 / 55 / 60 (strongest-20 15
   / 20 within 1%, 20 / 20 within 2%), 0.508 node-h (0.14x the reference's 3.64);
   gallery-06b passes with E 93 / 95 (both misses in the narrow class), p_SA 47 / 69 /
   95, p_MS 77 / 94 / 95, p_MA 68 / 91 / 95 (19 / 20, 20 / 20), 0.974 node-h;
   gallery-06 (before decision 44) fails p_MA (25 / 133 beyond 5%) - the table
   reproduces the RESULTS.md class counts (tests); gallery-10 (reference p5, gated at
   its p5 main stage, p_SA not applicable) reproduces its RESULTS.md p5 row - E 78 / 78,
   p_MA 33 / 59 / 78 (median +1.28%), p_MS 75 / 78 / 78 - and FAILS the p_MA
   strongest-20 statement at the two z = 0.1 near-junction hats 53 / 58 (+2.8 / +3.5%),
   the finding RESULTS.md reports as the systematic far-surface MA offset at equal p.
6. **Records.** `ROOT/library-qualification.json`: per coupon `Status` (qualified /
   pending-qualification / failed / planned / skipped) and `StoppedBy` (Build, Manifest,
   Reference, Mesh, Estimate, JobBudget, Monitor, Fetch, Stages, Verification,
   MatrixValidation), the reference binding, sources / classes / controls, stage layout,
   estimate (H1 by order, job seconds by PCG factor, peak GB, node-h of the main stage),
   plan (pins, caps), remote layout, submission / monitor / fetch / digest / matrix /
   deletion records, `Qualification` (verdict, gates passed, not-applicable gates, anchor,
   gated stage / order / comparison, class statistics, offsets, weighted p_MA, the other
   main orders as `Informational`), `Cost` (library-order main-stage H1 / PCG / seconds /
   node-h, every main stage under `MainStages`, job node-h, the reference's node-h from
   its `status.json`, the ratio); library totals: coupons by status, stopped coupons with
   the reason, node-h, critical-path seconds from the first submission to the last fetch
   (jobs overlap up to `--max-jobs`), per-job wall seconds, jobs submitted vs `--max-jobs`
   and the cap (every qsub counted at submission, a stop after it included), orders,
   binary hash, profile. `ROOT/qualification-gates.json` (the table used),
   `ROOT/process-library.json` (each coupon's model from its own `process-library.json`
   with the fetched matrices, `CouponMesh`, `Qualification` and `LibraryQualified` only
   when Passed; `--merge-into PATH` = a previous run's `process-library.json`, read only:
   its models this run did not qualify are kept ahead of this run's, a model of the same
   `Name` is replaced, `MergedFrom {Path, SHA256, Root, Kept, Replaced}` recorded - a
   coupon qualified later joins its device library without re-running the others).

`--dry-run` writes steps 1-3 and the gate table without contacting anything; the
recorded campaigns are its fixtures: on the four-edge case with `--stage-prefix va`
and the physics-11 controls, every generated `worker.json` / `reducer.json` /
`config.json` equals `four-edge-physics-11/main/*/` apart from paths and the plan
equals its `plan.json` in stages (names, config files, environment, dependencies,
order) and pins the same 80 trace files (their digests are the regenerated traces',
recorded; the mesh pin is the build record's); the same on the gallery-06 case against
`gallery-physics-06b` (135 traces). Tests:
`test_qualify_gates.py` (the gate evaluation, the recorded class table and locations,
the estimator against the recorded 06b `stage-estimate.json`, node-h, the cap rule,
the control choice), `test_qualify_dry_run.py` (the dry runs above, the analysis of
the recorded results through the same records - physics-11 / 06b Passed, gallery-10
gated at its p5 main stage with p_SA not applicable -, the concurrent scheduler against
a fake remote replaying the recorded trees, PendingQualification, the fail-closed
stops, `--reference none`, the derived-config equality on the five gallery cases, the
altered-reference stop). Nothing four-edge-specific remains hard-coded: no source
count, control index, remote path, interface index map or mesh digest is in the code.

**Still manual between a device layout and a qualified library.** (1) The remote
prerequisites: the frozen executable `palace-archive-estimate-<sha>.bin` and
`mpiexec_bound.sh` must pre-exist under `--remote ROOT`, and another cluster means
editing `qualify/cluster-profile.json`. (2) The etch footprint declaration of a
registered case (`--footprint bound` with a `retained-etch.csv`, or `producer-default`)
is a deliberate per-case statement (decision 16), never inferred. (3) A device whose
discovery yields corner or straight-edge requirements gets those families from their
own builders (`corner_coupon/`, `cpw2d/`): `build --device` registers the
SpatialEdgeCluster coupons and records the others as out of this library's scope. (4)
A qualification without a reference is `PendingQualification` at best: an accuracy
statement needs a graded_v2 reference (or the decision-53 calibration path) and is not
produced by the command alone.

### Acceptance of the source split (supervisor decision 61b / 61d, 2026-09-21)

Two live runs, records under `qualify/split-acceptance-20260921/` and
`qualify/device-library-20260921/spatial-3-edge-5d3b5e644745/` (7f03).

**1. Split correctness - the two-edge 10 under `--job-policy fixed --fixed-jobs 2`** (78 sources,
mesh `5d01204e…`, `--reference none --controls p3,p5`; PBS 46337 worker-1 = block 1 sources 1-13 with
the p3 / p5 controls and the local-edge stage, 46338 worker-2 = block 2 sources 14-78, 46339 the
reducer on the union of 14,976 = 78 x 192 potentials). `qualify/compare_split_matrices.py` against the
single job 46023 of the 2026-09-20 acceptance (same mesh and p4 configs): domain Q_ij 3,081 entries,
3,042 bit-for-bit, largest per-entry relative difference 8.65e-13 (source 13) = 2.3e-14 of the largest
entry; surface Q_ij / normal / total 6,162 entries, 6,111 bit-for-bit, 8.47e-13 (source 53) = 3.4e-14
of the largest entry; R and the tangential columns exact - every difference one unit in the 12th printed
digit (the union's file order changes the summation order of the linear reduction, nothing else):
**equal to roundoff**. Cost: 339 + 272 + 102 = 713 s = 0.198 node-h over three jobs against the single
job's 696 s for the same p4 stages (+2.4 %, three job start-ups); compute path 441 s (worker-1 then the
reducer) against 696 s; measured `CriticalPathSeconds` 2,346 because the reducer job was held three
times by the SOCA dispatcher for r8g.48xlarge capacity (1,689 s of queue). The acceptance surfaced two
driver defects, fixed forward: a held job (PBS H) was read as "left the queue" and stopped the coupon
(`remote.IN_QUEUE_STATES`, commit 9cf2b837f), and `--resume` of a split coupon compared byte-identical
plans in two orders (commit 415ed0428); the third driver resumed and finished the run.

**2. 7f03 under `--job-policy speed --max-jobs 4` (decision 61d) - RUNNING at the time of this
commit.** `spatial-3-edge-5d3b5e644745` (225 sources, 3.04M elements, the coupon that did not fit
one 6 h job in the decision-58 run), `--reference none --orders p4 --controls p3,p5 --merge-into
/tmp/library-device-transmon-01/process-library.json`, root `/tmp/library-device-transmon-02`. Dry
run (identical to the recorded one): N = 4 of 4, blocks [17, 70, 69, 69], controls + local-edge in
worker-1, per-job 2x-PCG estimates 91 / 92 / 91 / 91 min and the reducer 149 min (the longest job),
critical path 241 min, node time 8.55 h (N = 1 does not fit at 29,272 s; N = 2 fits at 10,551 s
longest / 19,492 s path; N = 3 8,941 / 16,133 s). Live: PBS 46341-46344 submitted 21:47Z, all
running 21:53Z without a capacity hold; worker-1 2,026 s (block 1 500 s, p5 controls 786 + 117 s,
p3 controls 135 + 61 s, local-edge), worker-2 1,777 s, worker-3 1,598 s, worker-4 1,589 s (every
worker 2.3-3.4x under its 2x estimate); archive union 43,200 = 225 x 192 potentials; reducer job
46350 submitted 22:27Z and running since 22:28Z. The driver is detached (`nohup`, log
`/tmp/library-device-transmon-02/qualify.log`) and will fetch, verify the digests, delete the
archives, evaluate the p-sequence gates and write `process-library.json` merged with the five
decision-58 models; its outcome, cost (jobs, node-h, critical path) and records are the next
worker's evidence commit (resume with the same command plus `--resume` if the driver is lost).

### Live acceptance of the two commands (supervisor decision 51, 2026-09-20)

`qualify/acceptance-20260920/` holds the records of the live run that qualifies the two
commands (`ACCEPTANCE.md`, `library-build.json`, `library-qualification.json`,
`process-library.json`, per-coupon `qualification.json` / `cost-summary.json` / submission,
digest, matrix-validation and archive-deletion records; no CSV, mesh or log). `build`
rebuilt four-edge (identity byte-identical to the physics-13 mesh `1d536c44…`) and two-edge
(the decision-44 production mesh `5d01204e…`, not gallery-10's pre-decision-43/44 mesh) in
1913 s; `qualify` ran both coupons concurrently (PBS 46023 / 46024, 34:32 and 48:50,
critical path 3141 s, 1.389 node-h, archives deleted after the digest check). Four-edge:
Passed, CSVs bit-identical to physics-13 (E 60 / 60, p_MA 30 / 53 / 60 with 20 / 20 strongest,
p_MS 47 / 59 / 60, p_SA 34 / 46 / 58, 0.532 node-h). Two-edge, gated at its p5 main stage
against the p5 reference with p_SA not applicable: Failed on the p_MA strongest-20 at 53 / 58
exactly as the recorded gallery-10 CSVs do through the same gates (E 78 / 78, p_MA
34 / 60 / 78 vs 33 / 59 / 78 recorded, p_MS 77 / 78 / 78 vs 75 / 78 / 78; per-source
differences ≤ 2.1% = the recipe step 43 / 44 on the two-edge mesh; 0.118 / 0.382 node-h at
p4 / p5). Defects fixed forward during the run: not-applicable interfaces, sequential
submission, global `--orders`, `rsync --mkpath`, the runner uploaded as a tree, a lost driver
without `--resume`.

### Device-library blockers (supervisor decision 54, 2026-09-21): the 5-edge collar-union footprint (54a) and the needle cost of the device basis (54b)

**54a - the 5-edge mesher defect.** The transmon 5-edge coupon (`spatial-5-edge-ea1fbd054c1e`,
7 physical sides: two 6 um wide notches x in (-6, 0) and (10, 16) below y = 2.5 flanking a
10 um tooth) stopped its registration probe at `tube tool (3, 10) has 3 volume
descendants` (prism_edge_tubes.jl `tube_volume_after_fragment`). Reproduced on the source
directory (17 s, geometry only) and localized: tool (3, 10) is the bottom tube's vacuum
sectors of the side (-6, 2.5) -> (-6, -8) (the left notch wall, y from 2.452 to -7.984).
The cause is not a tube interaction but the producer-default etch footprint:
`boundary_strips` offsets the metal loop by the collar width 3R = 6 um with
`offset_loop_points`, a miter offset that shifts every Physical side by 6 um and joins
consecutive shifted lines. Two Physical sides facing each other across a dielectric gap
narrower than twice the collar fold the miter polygon back through the gap: here the left
notch wall x = -6 offsets to x = 0 and the notch's other wall x = 0 offsets to x = -6, so
the polygon runs (0, -8) -> (0, -3.5) -> (-6, -3.5) -> (-6, -8.5) -> (16, -8.5) -> ...: its
side (-6, -3.5) -> (-6, -8.5) lies ON the metal wall x = -6, the axis of that wall's tubes,
and it crosses the box side y = -8 (self-intersections at (-6, -8) and (10, -8); the OCC
warnings `BOPAlgo_AlertAcquiredSelfIntersection` / `FaceBuilderUnusedEdges`). The face OCC
builds from the self-intersecting wire acquires arbitrary seams, and the fragment splits the
tube volume along them. A survey of every plan-view loop of the manifest (five gallery cases,
the calibration input, the six Radius-2 fixtures, the four Radius-0.5 assessment fixtures
and the six device coupons) finds exactly one self-intersecting collar polygon: the 5-edge's
(`/tmp/coupon-five-edge-54a-20260921/survey`).

Fix (`mesh_spatial_coupon.jl`, generic): `collar_loop_points` builds the loft polygon of an
exterior loop: the miter polygon while it is simple (`polygon_is_simple`: no side shorter
than the tolerance or folding back onto its predecessor, no two non-adjacent sides within
the tolerance; construction `MiterOffset`), otherwise the outer boundary of the union of
`collar_pieces` - the loop, one rectangle per Physical side and one miter kite per convex
metal corner (the miter point of `offset_loop_points`; a right-angle box junction gives a
degenerate kite), every piece clipped to the coupon box by Sutherland-Hodgman
(`clip_polygon_to_box`) - traced by `polygon_union_boundary` (every piece edge split where
another meets it, a sub-segment kept when exactly one of its sides is inside some piece,
chained with the union on the left from the lexicographically smallest vertex;
construction `CollarUnion`). The two constructions describe the same region wherever the
miter polygon is simple, so the union form is entered only where the old one was invalid.
Fail-closed: a union boundary that is not one simple loop (an un-etched island inside a gap
narrower than twice the collar - a keyhole - or a pinch point) stops at the new recorded
scope guard `FootprintTopology` (`RECIPE_SCOPE_GUARDS`, `mesh_stage_contract.py`); an
inward (positive) self-intersecting offset has no union form and errors. Census: every
`FootprintPolygons[]` record carries `Construction` (`MiterOffset` / `CollarUnion` /
`HoleOffset` / `EdgeStrip`), `FootprintSimplification` the `ConstructionRule` and
`CollarUnionPolygons`. On the 5-edge the union is the whole box: every point of both
notches is within 6 um of a wall and every other exposed point within 6 um of the tooth or
the legs, so the footprint polygon is the box, the trench covers the whole exposed plane
(no 3000 label; 3100 117.8 + 3101 65.8 um^2), and every tube tool has one descendant.

Byte-identity of the verified roots (the union path is never entered): two-edge 10
(`two-edge-8dd4bc70f183`) and four-edge rebuilt under the fixed mesher with
`run_gmsh_only_case.py --stages-only` (`/tmp/coupon-five-edge-54a-20260921/`, binaries
deleted after hashing, `deleted-mesh-sha256.txt`): `gmsh-build.msh` `10afee1f...` (two-edge:
the recorded gmsh-build of the production identity `5d01204e...`) and `96645749...`
(four-edge: the recorded gmsh-build of the physics-13 / acceptance identity `1d536c44...`),
both BYTE-IDENTICAL; census `FootprintPolygons` all `MiterOffset`, `CollarUnionPolygons 0`.
Every other case has a simple miter polygon (survey above) and takes the unchanged path.
Regression tests (`test_prism_tube_build.jl`, testset "producer-default collar of facing
sides closer than twice the collar"): a 2D notch 8 um wide under a 6 um collar (miter
polygon self-intersecting, union = the half-plane y >= -6 of the box; pieces counted and
inside the box; the 16 um notch keeps `MiterOffset` with the identical polygon; zero offset
is the loop; inward self-intersection and a missing box fail closed; a keyhole stops at
`ScopeGuard[FootprintTopology]`), and a Radius-0.5 coupon build whose notch is exactly one
collar (1.5 um) wide - the device failure mode: verified to stop with `tube tool (3, 5) has
2 volume descendants` under the old construction on a scratch copy - now builds with
`CollarUnion`, the box as footprint, 16 tubes / 21 matched tube volumes, no 3000 plane.

Device 5-edge under the production recipe (`/tmp/coupon-five-edge-54a-20260921/device`, a
fresh copy of the production manifest; `register_case.py` census-only probe, then
`run_gmsh_only_case.py`): estimate 3,452,743 (0.863 of the cap; est / actual 1.126); built
**3,066,661** elements = 2,793,997 tets + 253,188 prisms + 19,476 pyramids, gmsh-build 195 s
/ 6.54 GiB (the probe 177 s / 7.09 GiB - the 8 GiB bound is close), 14 tubes / 21 tube
volumes, min SJ 0.0319, tet condition 126.5, prism condition 589.1, corner aspects <= 3.72,
interfaces 1 1174.4 / 3100 117.805 / 3101 65.845 / 5001 199.555 / 5101 35.445 / 6001 203.174
/ 6101 37.126 um^2, basis 28 needles (14 below FarSize, min altitude 11.8 nm) - the ear-
clipped device basis; identity `identity.msh` `6b94299e...` (gmsh-build `3220b214...`,
rotate-z `0663a207...`; every mesh gate passed; CanonicalBuildId `53dc9f55...`). The case
then FAILED CLOSED at per-entry verification (1,601 s of the 1,800 s bound): both
variants `gate failure: trace-diagonal-overrefinement` - two global short-edge bands of
27,613 / 27,647 vertices along the box edge y = -8 on the matching-surface bottom and top
faces (z = -2.05 / 2.1), spanning the 26 um box, RMS width 24 nm, tilted 0.0057 rad from
the edge against the band's own resolvability 0.0009 rad, aligned with no feature within
that resolvability. They are the needle bands of the ear-clipped device basis: its cap fan
joins the y = -8 ring vertices (2-5 um apart) to an apex 0.153 um off the line at the
opposite corner (altitudes 11.8 / 19.1 / 38.3 nm over 26 / 24 / 21 um long sides; tilt 0.153
/ 26 = 0.0059 rad), so the trace rule requests 6 nm along the whole edge - the needle cost
of 54b, here failing a production gate rather than the cap. Recorded fail-closed
(`build-summary.json`, `per-entry-verification.json`); the root's mesh binaries were deleted
after hashing (`deleted-mesh-sha256.txt`). The same coupon with the Delaunay device basis
(54b below) is `spatial-5-edge-9bfa8265e2e7` (same sources, `trace-triangles.csv` and
`basis-contract.json` the only differing files; 8 needles, 2 below FarSize, min altitude
46.8 nm): estimate 2,905,085 (0.726 of the cap; est / actual 1.077), **BUILT AND VERIFIED** -
2,696,939 elements = 2,424,275 tets + 253,188 prisms + 19,476 pyramids (0.88 x the
ear-clipped build), gmsh-build 190 s / 5.55 GiB, verification 1,335 s, min SJ 0.0229, tet
condition 100.6, prism 589.1, corner aspects <= 3.66, 0 global diagonal bands (12 / 12 / 4
aligned with signature / footprint / junction), interfaces 1 1174.4 / 3100 117.861 / 3101
65.789 / 5001 199.56 / 5101 35.44 / 6001 203.179 / 6101 37.121 um^2, identity `identity.msh`
`ddf1081b694fba74f604e79442e5dc82ce52b1821b05273ca8b24e2281452895` (gmsh-build `db66deed...`,
rotate-z `d209228f...`), CanonicalBuildId `c136a66e...`, both variants Passed, transform
comparison Passed.

**54b - needle cost (`analyze_basis_needles.py`, analysis; record
`/tmp/coupon-five-edge-54a-20260921/needle/*.json`).** The tool decomposes
`estimate_build_cost.estimate` into its components and splits the trace-basis law by basis
triangle class - needle: minimum altitude < 0.6 x shortest edge (trace_basis.py); right
sliver: any other triangle with requested size 0.5 x altitude < FarSize; wide: the rest - and,
given a built root, attributes every tetrahedron to the size law governing at its centroid
(tube band, corner / cap ball, junction band, trace-basis per class, FarSize; ties at FarSize
to the far field). Where the needles come from: the producer's `cap_ring` ear-clips the two
box caps; a ring with clustered vertices (the ten-edge: three vertices 56 nm apart on the
-y side) leaves ears that join the short edge to a vertex far along the same side (19.4 um
away at 5.6 deg: altitude 56 nm x sin 5.6 deg = 5.5 nm), and the trace rule then requests
2.7 nm over a 39 um perimeter. The side strips (56 nm x 2 um) are right slivers (altitude =
the ring spacing), the class TraceBasisSizeRatio 0.5 was calibrated on.

| case | sources | estimate (cap 4M) | needles / right slivers / wide | min altitude needle / sliver | needle tets (share) | right-sliver tets | far / tube band / prisms+pyramids / balls / junction |
|---|---|---|---|---|---|---|---|
| three-edge 7f03 (`spatial-3-edge-d96d52b95001`) | 225 | 5,025,593 (1.256) | 38 / 224 / 184 | 1.0 nm / 15.1 nm | 1,977,952 (**39.4%**) | 310,611 (6.2%) | 1,814,861 / 409,725 / 397,908 / 36,732 / 77,804 |
| ten-edge gallery (`ten-edge-6791f1c84123`) | 260 | 3,561,695 (0.890) | 34 / 268 / 214 | 5.5 nm / 35.1 nm | 744,253 (**20.9%**) | 330,017 (9.3%) | 1,721,083 / 290,773 / 282,744 / 118,673 / 74,151 |
| ten-edge production root `8a871f22...` (actual, by governing law) | 260 | 3,498,453 built | - | - | 921,337 (**26.3%**) | 225,603 (6.4%) | far 1,620,702 / tube band 363,132 / prisms+pyramids 290,304 / balls 77,317 / junction 58 |
| device ten-edge (`spatial-10-edge-1ce0c327a734`, 8.03 GiB probe) | 300 | 3,625,496 (0.906) | 40 / 304 / 252 | 5.5 nm / - | 794,975 (21.9%) | 343,097 | - |
| device 5-edge | 185 | 3,452,743 (0.863) | 14 / 182 / 170 | 11.8 nm / - | 538,214 (15.6%) | 255,192 | 1,951,317 / 275,196 / 267,120 / 87,592 / 78,113 |
| four-edge gallery | 80 | 1,966,879 | 0 / 78 / 78 | - / 15.3 nm | 0 | - | - |

The needle share of the estimate is 39.4% on 7f03 and 20.9% on the ten-edge (26.3% of the
built ten-edge by governing law; the estimator charges the needles' Steiner shells at
their perimeter and under-attributes them relative to the built mesh, whose trace-law tets
also include the shell overlap with the far field).

Basis constructions for DEVICE coupons (the library's own basis; never the gallery
references):

(a) **Delaunay caps** (`generate_spatial_response.py --cap-triangulation delaunay`;
`delaunay_flip_cap`, Lawson edge flips of the ear-clipped caps to the max-min-angle
triangulation of the SAME vertices; ~40 lines; default `ear-clipping` = the frozen
producer, so every gallery and every existing device directory is unchanged; plumbed
through `device_coupons.py` / `coupon-library build --device --cap-triangulation`,
recorded in `basis-contract.json` `CapTriangulation` {Method, MinimumAltitude,
NeedleTriangles}, the device record and the provenance; the content hash names a new
case). The short cluster edges are joined across the face, so every narrow altitude is a
ring spacing times the sine of the join angle: 7f03 min altitude 1.0 -> 12.8 nm, needles
38 -> 2 (two 0.18 um corner triangles), estimate 5,025,593 -> **3,135,224** (0.624 x; 0.784
of the cap); ten-edge 5.5 -> 32.0 nm, 34 -> 8 (0.107 um cluster edges joined at 35 deg,
altitude 62 nm), 3,561,695 -> 2,795,156 (0.785 x); device ten-edge 3,625,496 -> 2,837,463
(0.783 x); device 5-edge 3,452,743 -> 2,905,085 (0.841 x). Sources unchanged, so the
physics cost index (sources x elements) falls by the same factor. `analyze_basis_needles.py`
reproduces the producer's result exactly: the estimate of the directory generated with
`--cap-triangulation delaunay` equals the tool's re-triangulated estimate (3,135,224).
Tests: `test_box_trace.py` (the gallery ten-edge testdata basis: 34 -> <= 8 needles below
0.32 um, min altitude 5.5 -> > 30 nm, sides untouched, closed oriented surface, fixed point
of the flips; `build_matching_surface` with either method: same vertices and closure, an
unknown method fails closed).

(b) **Graded refinement** (modelled, not built): a basis whose local edge on every face is
the ring spacing graded at slope BasisGrowth (h = min_i spacing_i + BasisGrowth x |x - v_i|),
interior sources 2 / (sqrt 3 h^2) per area, trace cost the volume-law shell with s = 0.5 x
(sqrt 3 / 2) h. The refinement cannot raise an altitude above the ring spacing (the hats
of clustered ring vertices vary over that spacing whatever the interior), it only shortens
the narrow hats' supports, and every added vertex is a Palace solve:

| case | BasisGrowth | sources | estimate | elements x | sources x | physics cost (sources x elements) |
|---|---|---|---|---|---|---|
| 7f03 | 0.5 / 1.0 / 2.0 | 3,092 / 1,308 / 625 | 2,763,814 / 2,745,369 / 2,739,686 | 0.550 / 0.546 / 0.545 | 13.7 / 5.8 / 2.8 | **7.6 / 3.2 / 1.5 x** |
| ten-edge | 0.5 / 1.0 / 2.0 | 3,281 / 1,444 / 711 | 2,515,615 / 2,496,585 / 2,490,453 | 0.706 / 0.701 / 0.699 | 12.6 / 5.6 / 2.7 | **8.9 / 3.9 / 1.9 x** |

Against Delaunay caps (0.624 / 0.785 x at 1.0 x sources) the refinement buys a further
0.08-0.09 x of elements for 2.7-14 x the sources: the trade is lost at every slope.
**Proposal:** adopt (a) for device coupons - `--cap-triangulation delaunay` as the device
basis rule - and do not pursue (b). **Adopted (supervisor decision 57, 2026-09-21):**
`delaunay` is the `build --device` default (`device_coupons.DEFAULT_CAP_TRIANGULATION`;
`generate_spatial_response.py` itself keeps `ear-clipping` as its default because it is
the gallery producer); every device case id and hat set changed with it (new cases);
`--cap-triangulation ear-clipping` remains an explicit option; no gallery reference
case is touched (`test_device_coupons.py`). The remaining
right-sliver cost (6-9% of the estimate) is the ring spacing itself (16-56 nm clusters
inherited from the device mesh) and would need ring coarsening, a source-resolution
decision, not a triangulation one. Note the ten-edge's 8.03 GiB probe is a memory-bound
question the Delaunay basis addresses only through the element count (0.785 x).

Build of 7f03 with the Delaunay basis (`spatial-3-edge-5d3b5e644745`: the recorded
`coupon.json` of the 2026-09-20 device run regenerated with `--cap-triangulation delaunay`;
every bound file byte-identical to `spatial-3-edge-d96d52b95001` apart from
`trace-triangles.csv` and `basis-contract.json`; registered into the same manifest copy):
estimate gate 3,135,224 (0.784 of the cap) PASSED where the ear-clipped basis failed closed
at 1.256; probe 3,037,912 elements in 194 s / 6.35 GiB (est / actual 1.032). **BUILT AND
VERIFIED** (`spatial-3-edge-5d3b5e644745-attempt2`): 3,037,912 elements = 2,635,846 tets +
373,347 prisms + 28,719 pyramids, gmsh-build 193 s / 6.51 GiB, min SJ 0.0244, tet condition
108.6, prism 589.4, corner aspects <= 3.79, 0 global diagonal bands, interfaces 1 1115.836 /
3100 268.309 / 5001 122.545 / 6001 130.436 um^2, identity `identity.msh`
`e038c5efbb50145aacc18f768317344f03b8af1a47e0482a0296d9ab1ffe62fc` (gmsh-build `25fe87f3...`,
rotate-z `ed1ea942...`), CanonicalBuildId `a9a0eaf4...`, both variants Passed; verification
1,561 s of the 1,800 s bound. A first attempt of the same case produced byte-identical
meshes and stopped at the verification TIME bound (1,800 s, `StopReason timeout`) while the
full unit-test sweep and a second verification ran concurrently
(`spatial-3-edge-5d3b5e644745/build-summary.first-attempt-timeout.json`, meshes deleted after
hashing): the bound is a machine resource guard (decision 54c) and a 3M-element coupon's
verification sits within 15% of it - the ten-edge's 1,607 s again.


### First complete device library run (supervisor decision 58, 2026-09-21)

`qualify/device-library-20260921/` (`DEVICE-LIBRARY.md`, `library-build.json`,
`device-coupons.json`, the dry-run `library-qualification.dry-run.json`, the live
`library-qualification.json` / `process-library.json` / `qualification-gates.json`, per coupon
`submission.json` / `stage-estimate.json` / `qualification.json` / `cost-summary.json` /
`matrix-validation.json` / `remote-archive-deletion.json` / `result-csv-sha256.json` /
`qstat-xf.txt` / `comparison/p-sequence-controls.{json,md}` / `comparison/class-statistics.json`;
no CSV, mesh or log) records the transmon layout run
through both commands with the decision-57 Delaunay device basis. **Build** (`build --device`,
fresh production-manifest copy, no `--build-limit`, pool of 2, one workstation): discovery +
basis generation + registration of all SIX spatial coupons in 990 s (every census-only probe
under the bounds; the ten-edge probe 177 s / 5.94 GiB where the ear-clipped basis hit 8.03 GiB),
build pool 7,338 s, end to end 8,493 s; **6 / 6 built and verified** (identity + rotate-z
Passed, every gate at its production value): 2-edge 1,215,960 (H1 p4 17.1M), 3-edge `9cd9`
1,680,987 (22.9M), 4-edge 1,895,487 (25.2M), ten-edge 2,643,905 (35.9M; 0.709 of the cap,
gmsh 5.92 GiB), 5-edge 2,696,939 (35.9M), 7f03 3,037,912 (42.7M; 0.784 of the cap; the one
headroom flag: verification 1,692 s = 0.94 of the 1,800 s bound under the concurrent 5-edge
build, decision 54c); estimate / actual 1.03-1.08 everywhere. Nothing fails closed in the build
any more: the three 2026-09-20 stops (7f03 1.256 x cap, ten-edge memory, 5-edge tube tool) are
removed by 54a / 54b. `build` does not reuse roots by content: the 5-edge and 7f03 were rebuilt
and came out BYTE-IDENTICAL to the decision-54 identities `ddf1081b…` / `e038c5ef…`.
**Qualify** (`--reference none --orders p4 --controls p3,p5 --max-jobs 2`): the dry run plans
five coupons in one 6 h job each (runner minutes at the measured PCG counts / at 2.0x with
preflight and margin: 2-edge 59 / 130, 4-edge 69 / 154, 3-edge 94 / 201, 5-edge 158 / 332,
ten-edge 164 / 342; main-stage node-h 0.79 / 0.86 / 1.30 / 2.21 / 2.31 = 7.47 at 1.0x; Palace
peaks 135-285 GB of 1,485) and fails 7f03 closed at the estimate: 225 sources x 42.7M H1 (p4) /
83.1M (p5) = 21,460 s at 2.0x PCG, 29,272 s with preflight and margin > the 21,600 s walltime -
a job-size limit of the one-job-per-coupon plan (the mesh is built and verified), resolvable
only by a two-job split of the source set (unsupported) or a longer walltime (user decision).
The live run (the five, largest first) submitted PBS 46219 (ten-edge) and 46220 (5-edge) at
11:28Z; from 11:38Z the cluster login host was unreachable for more than two hours (the two jobs
kept running), the driver read the empty poll as "left the queue", its fetch failed and it
crashed - fixed forward (transport-failure polls and fetch stops above, commit `5a25d034e`) and
the run resumed with `--resume` on the same root at 15:48Z: both submissions adopted
(`Monitor.Resumed`, no second qsub), fetched at 15:52 / 15:55Z; 46221 (3-edge) and 46222 (4-edge)
submitted into the freed slots, a second login-host outage 16:21-16:40Z recorded as 8
`Monitor.TransportFailures` on each (ssh rc 255, jobs kept active, 50 / 38 of 240 polls spent),
both fetched; 46255 (2-edge) last, fetched 18:11:38Z. **Outcome: 5 / 5 PendingQualification, 0
Failed, 0 Passed** - every p-sequence gate passed at the gated order p4 (8 control sources per
coupon; max |d45| of E 0.09-0.49% against 1%, of the participations 1.21-2.68% against 5%: ten-edge
p_MA 1.22%, 5-edge p_MS 1.69%, 3-edge p_MA 2.68%, 4-edge p_SA 1.21%, 2-edge p_MS 1.63%; every
stage complete, every matrix validated, 14 digests per coupon, all response archives deleted on
the cluster); `process-library.json` carries the five models with `LibraryQualified false` - the
rule "never Passed without a reference" holds, so nothing is qualified until a reference exists.
Cost: PBS job wall 8,618 / 8,501 / 5,241 / 3,984 / 3,380 s (ten-edge, 5-edge, 3-edge, 4-edge,
2-edge; PCG mean 14-18, max 39-40), **8.257 node-h** (main stages 6.82; the 1.0x-PCG estimate
4-14% above the job, the 2.0x-with-margin gate value 2.3-2.4x); `CriticalPathSeconds` 24,313
(11:26:25Z first submission -> 18:11:38Z last fetch), of which 7,158 s is the first outage's dead
time (46219 / 46220 ended 13:53 / 13:51Z), 17,155 s of job + queue chain; 5 jobs submitted, 2
concurrent, user cap 40. Every MA value of this run is the sharp-edge value at the recipe's 0.25
nm cutoff (decisions 55 / 56). Manual after the run: the 7f03 job split or walltime, the MA
definition, the reference for the device geometries (the accuracy statement), the corner and
isolated-edge requirements outside the spatial scope.

### Performance (user decision 62, 2026-09-21): reducer block size, labels-only probe, streaming Gram

Records under `qualify/perf-20260921/` (`step1-block-size/`, `step2-labels-only/`,
`step4-streaming-gram/`, `linux-build/`; no CSV, mesh, archive or binary). Every acceptance is the
smallest sufficient run: the two-edge-8dd4bc70f183 p3 control (8 sources) / p4 stage (78 sources) on
the cluster, the six transmon coupons locally. Physical gates, quadrature rule / order / weights /
sample set / ownership are unchanged; the rotate-z policy and the quadrature order were not touched.

**Step 1 - reducer block size 6 -> 48 (`PALACE_RESPONSE_BLOCK_SIZE`).** The reducer reads,
differentiates and evaluates every archived source once per block pair its block takes part in
(N x ceil(N / b) evaluations): b = 6 evaluated every source 32 times at 191 sources. Recorded as
`ProductionRecipe.PhysicsRun.ReducerBlockSize {Value 48, PreviousValue 6, Rule, Provenance}` (the
memory rationale: 2b resident fields x ~14 MB per rank at p4, b = 48 -> est. 350-450 GB of 1,485 GiB;
the `MinimumMemAvailableBytes` admission guard unchanged), `build_plan.DEFAULT_REDUCER_BLOCK_SIZE`
/ `reducer_environment(b)`, the plan's `ReducerBlockSize` (PLAN_VERSION 3), `qualify
--reducer-block-size` (the origin recorded per coupon: command line / manifest / built-in default);
`estimate_stages` takes the block size (cost model `MeasuredBlockSize 6`, the block-pair part split
into the evaluation part N x ceil(N / b) and the pair-scaled Gram by `ReducerEvaluationFraction`;
the reducer peak grows by `ReducerResidentFieldGBPerMillionH1` per resident field); `summarize_cost`
reads the block size the reducer stage ran at. Acceptance (PBS 46685, frozen b28, the p3 control and
the p4 stage archived once each and reduced at b = 6 and b = 48 from the same archive): every
surface- / domain-response-matrix VALUE bit-identical (p4: 36,972 surface + 3,081 domain entries;
p3: 432 + 36); the CSV row ORDER differs (block-pair emission order; every consumer groups rows by
key). p4 reducer wall 99.2 -> 28.4 s (Operator Construction 86.9 -> 15.3 s; 91 -> 3 block pairs),
Palace peak 45.6 -> 91.1 GB (node used 90.8 -> 136.5 GiB); p3 12.2 -> 11.8 s (setup-bound).
Calibration from the pair: 0.0835 s per source evaluation, 2.2 s of pair-scaled Gram -> evaluation
fraction 0.974 (recorded 0.97, the Gram remainder rounded up); 0.542 GB of peak per resident field
at H1 7.97M (0.068 GB per million H1 DOFs). Projected library reducers (the decision-58 device
library, b28 executable, from the recorded b = 6 times with the 0.974 split): 10-edge p4 3,969 ->
636 s, 5-edge 3,956 -> 647, 3-edge 2,020 -> 344, 2-edge 992 -> 203, 4-edge 1,017 -> 207; controls
53-116 -> 51-99 s; all reducer stages 12,592 -> 2,607 s (-2.77 node-h of the 8.26 node-h run).

**Step 2 - labels-only registration probe and parallel registration.** The registration learned
the label set (the un-etched plane 3000 + slot per slot) from a full census-only production build
of every coupon (947 s of the 990 s stage). `mesh_spatial_coupon.jl --labels-only OUT.json` stops
right after the CAD-entity labelling - after every `occ.fragment`, before any `mesh.generate` -
and writes `{InterfaceAreas: [{Attribute, Name, CADSurfaces}]}` with the recipe scope record, the
semantic contract identity, the corners and the box; `run_gmsh_only_case.py --labels-only` runs
the headroom gate, the source validation and the gmsh-build mesher command plus `--labels-only`
under the stage bounds (outside the frozen gmsh-build stage contract: no mesh artifact);
`register_case.run_probe_build` uses it, `derive_semantic_contract` consumes the census unchanged
(`BuildCensusInterfaceLabels`), the production build still fails closed on a label mismatch
(`mesh_stage_contract.validate_gmsh_build_census`). Registration is split into
`prepare_registration` (digests, scope, probe, derivation; the manifest read only) and
`commit_registration` (the serial manifest append; the version decided at commit);
`device_coupons.register_device_sources` prepares coupons in a pool (`--register-jobs`, default 2;
`coupon-library build` too) and commits in coupon order (`RegistrationPool` recorded). Acceptance
(local): the six transmon coupons re-registered from copies of their recorded source directories
into a scratch manifest (pool of 3): identical case ids, identical contracts (BoundaryLabels,
SemanticCorners, FeatureTopology, ProtectedSupports, CutSurfaceRoles, MetricSurfaceRoles,
UnmatchedPolicy, VolumeMaterials; the only difference `Derivation.BuildCensusSHA256`, the digest
of the census file itself), identical source role digests, FixtureVersion 1; labels-only pass
14.5-15.4 s / 0.9 GiB per coupon; registration wall 990 s -> 32.9 s. Existing registrations are
untouched (reuse by content); a coupon registered anew carries a new contract digest (its
`Derivation.BuildCensusSHA256`), hence a new CanonicalBuildId.

**Step 3 (audit dedupe / caching)** belongs to the main tree (not part of this block).

**Step 4 - streaming one-pass Gram in Palace (a new frozen executable).** In the reduce-only path
`SurfacePostOperator::CacheInterfaceResponseSamples` caches, per localized interface, the
boundary elements, quadrature rule, points, physical weights, edge distances and ownership
selection of the batched traversal once; `EvaluateInterfaceResponseRow` evaluates one field once
into its sqrt(0.5 w)-scaled 4-vector amplitude row; `AssembleInterfaceResponseMatrices` forms S = F
W F^T per interface (the total and the inside-radius window variants) with a cache-blocked
symmetric rank-k update over contiguous rows and reduces over the ranks.
`ElectrostaticSolver::PostprocessArchivedResponseMatrix` streams every archived source once (read,
gradient, face-neighbour exchange, evaluation into its rows, one `M_elec` matvec against the
resident potentials for the domain Gram), assembles at the end and emits the CSV rows in (i <= j)
order; the block size no longer changes the work (recorded in the log); memory = base + N x (sum
over interfaces of 4 Q_local + L_local) x 8 bytes (the resident potentials, one L-vector per
source, are what keep the domain Gram at one matvec per source); reduce-only skips the stiffness
assembly and `ksp.SetOperators` (the AMG setup); BlockTimer sections `Archive Reduction` /
`Archive Read` / `Sample Evaluation` / `Gram Assembly` / `Domain Gram`; per-rank sample counts
(min / max / total per interface) logged and parsed by `run_stages` (`ReductionSamples`,
`StreamedSourcesProgress`). Unit test: the streamed matrices equal the batched
`GetInterfaceElectricFieldEnergyMatrices` to 1e-12 (ownership-partitioned interfaces, total and
inside variants, serial and 2 ranks). Local smoke (the device 2-edge coupon, p2, 8 sources, 6
ranks, the main-tree binary as OLD): max per-entry relative difference 7.1e-13 (domain) / 4.1e-13
(surface) - the CSV print resolution; reducer 25.5 -> 16.5 s, peak 10.5 -> 6.8 GB. Linux frozen
executable: `qualify/perf-20260921/linux-build/` (the recorded freeze / build procedure re-pointed
at `source-freeze-perf62`, HEAD b1e7e9e9d, tar SHA256 b3103728…; PBS 46717, 177 s, six jobs, 1,856
provenance inputs unchanged) -> `palace-archive-estimate-170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b.bin`
under the remote root, SHA256 `170439c4…`; the recorded archive-estimate synthetic tests pass with
it (PBS 46720, 12 tests). **Acceptance (PBS 46718, `step4-streaming-gram/acceptance.json`, block
size 48 everywhere):** the new executable reduced the PBS 46685 archives - p3 control (8 sources)
bit-identical to the b28 reduction (468 entries), p4 (78 sources) max per-entry relative difference
8.4e-13 (domain) / 3.9e-13 (surface), 3,041 / 3,081 and 36,956 / 36,972 entries bit-identical (the
differences are the last printed digit of the 13-digit CSV); a fresh p4 worker archive reduced by
b28 and by the new executable in the same job gives the same numbers, and its b28 reduction equals
PBS 46685's bit for bit (the worker is deterministic across jobs). Reducer wall old (b28, b = 48)
-> new: p4 28.1 -> 15.5 s (Palace 26.3 -> 13.8 s; Operator Construction 15.6 -> 0.33 s - the
stiffness / AMG setup gone - the reduction itself ~5 s: archive read max 4.2 s, sample evaluation
max 4.2 s, Gram 0.2 s, domain Gram 1.1 s; the remaining ~9 s are the mesh preprocessing, the
initialization and the disk IO of any Palace run), Palace peak 91.2 -> 33.3 GB (node used 136 ->
76 GiB); p3 12.0 s / 22.0 GB. Per-rank samples: min 0, max 8,889 of 438,064 total (the interface
faces sit on a subset of the ranks - the sample evaluation is 0.23 s on the lightest rank and 4.2
s on the heaviest: the next lever, if any, is the surface-face balance of the partition). Against
the b = 6 baseline of the same stage (PBS 46685: 99.2 s, 45.6 GB): 6.4x faster, 0.73x the peak.
**Projected library reducers with the new executable** (the decision-58 device library, from the
recorded b = 6 times: the setup ~ the worker's non-source time minus the stiffness / AMG setup, one
evaluation per source at the measured per-evaluation rate, the pair-scaled Gram unchanged - an upper
bound, since the Gram kernel is also faster): 10-edge p4 3,969 -> ~280 s, 5-edge 3,956 -> ~280,
3-edge 2,020 -> ~150, 4-edge 1,017 -> ~110, 2-edge 992 -> ~90; all reducer stages 12,592 -> ~1,480
s (-3.1 node-h of the 8.26 node-h run; b = 48 alone: -2.77). The reducer stops being the deadline
floor of the source split (`job_split`): a coupon's job is the worker.

**Integration (decision 63, 2026-09-22; merge commit 4fb1f26d7 of simlapointe/perf-reducer
af52cb0d8..79a2500b9 into the MA-shell tree 4ef0a21ac).** (1) The frozen executable of every
`qualify` stage is now `170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b`:
recorded as `ProductionRecipe.PhysicsRun.FrozenExecutable {SHA256, PreviousSHA256
b28f089ae12c25863493566b2b8ca11af2c8ffb0e273e7aa67a2b42046eacf27, Rule, Provenance}` (validated by
`general_mesh_manifest.validate_physics_run` and `case_inputs.physics_run_parameters`),
`build_plan.DEFAULT_FROZEN_BINARY_SHA256` / `PREVIOUS_FROZEN_BINARY_SHA256` / `FROZEN_BINARY_RULE`;
`--frozen-binary-sha256` is optional and overrides (origin recorded per coupon and in the library
record). (2) The reducer block-size rule records that with the streaming executable the work is
independent of b and b = N (one block; N x (4 Q_local + L_local) x 8 bytes resident) is permitted;
the default stays 48. (3) The cost model's reducer rates remain the b28 / b = 6 measurements (a
conservative upper bound for the streaming reducer; re-measuring them is not part of this
integration). (4) The worker path is unchanged in the new executable; the merged qualify's live
acceptance (below) compares its p3 control matrices to PBS 46718's. The three manifests were
re-frozen after the merge (`general_mesh_manifest.py` ff44e45b8d2b -> be9ed95f8f2c) and again with the
`FrozenExecutable` validation (-> 9bb1b89600c1).

**Step 3 - audit dedupe / caching (decision 62(3); commits b88394b4e..296c403f9, refreeze 5590074f8;
records under `qualify/perf-20260921/step3-audit-dedupe/`).** Seven commits, one per proposal of the
performance review, each leaving every record bit-identical: (4) `mesh_array_io.read_mesh` parses a
MSH 2.2 binary file directly into the meshio.Mesh `meshio.read` returns (Gmsh writes one header per
element; meshio issues one `numpy.fromfile` per element: 10.4 s per read of the 602,862-element
two-edge mesh, 37-40 s of the 2.1M-element four-edge / three-edge meshes, against 0.01-0.06 s;
checked field by field on 27 recorded production meshes and by `test_mesh_array_io`); (1) the
producer memoizes `simplicial_view` / `volume_quality` / `analyze` / `_mesh_invariants` / the
corner-frame aspects per mesh object, the identity covariance comparison uses the canonical mesh
object itself once the inverse-transformed candidate equals it bit for bit, and bit-identical
protected patches normalize once; (2) verification reads the audited mesh once per variant and
passes it to every check (the canonical mesh cached across the case's variants); (5) the collinear
boundary-vertex removal runs on a min-heap of candidates (the same removal sequence as the restarted
scan), the material adjacency is judged per distinct (low, high) pair, the prism / pyramid /
quadrangle splits are vectorized (checked equal cell by cell) - the edge lengths of
`_global_diagonal_bands` stay scalar: the 1-D BLAS norm differs from a vectorized norm in the last
bit on this machine and `ShortEdgeThreshold` is their median; (8) the H1 entity counts are computed
once per volume-connectivity digest (`<case>/h1-entity-counts.json` for the producers of one case;
the verifier never reads it); (6) the Julia ownership audit snapshots the mesh as tag-sorted arrays
instead of one Dict entry per element (old and new agree element by element; the same negatives;
the mesher include, 2.1 s, stays); (7) the rigid publication proves the $Nodes-only change at the
byte level and parses the output once (the canonical mesh is that parse with the canonical
coordinates), receipts unchanged. Acceptance: two-edge-8dd4bc70f183 rebuilt before (f03c7eabd) and
after, compared with `compare_case_records.py` (41 records field by field, 18 meshes / CSVs byte
for byte; only timings, resources, tool digests and the ids derived from them, paths, commands and
environments are skipped): **0 differences**.

| two-edge 10 stage | before (s) | after (s) | peak GiB before -> after |
|---|---|---|---|
| gmsh-build | 44.4 | 41.6 | 2.51 -> 2.33 |
| canonical-publish | 14.3 | 13.6 | 1.49 -> 1.55 |
| identity publication | 54.7 | 16.5 | 1.73 -> 1.56 |
| rotate-z publication | 52.5 | 16.4 | 1.70 -> 1.56 |
| identity variant audits | 55.6 | 12.3 | 1.02 -> 0.86 |
| rotate-z variant audits | 70.6 | 12.2 | 0.92 -> 0.86 |
| per-entry verification | 213.8 | 21.7 | 0.95 -> 0.88 |
| **build wall** | **523** | **138** | |

**Integrated before / after and the projected device library.** Registration: 990 -> 32.9 s (step
2, measured on the six transmon coupons). Build: applying the measured two-edge ratios
(publications 0.31, audits 0.19, verification 0.10; gmsh-build and canonical-publish unchanged) to
the decision-58 device library's recorded stage times gives the serial case sum 13,145 -> ~2,940 s
(10-edge 2,654 -> ~610, 5-edge 2,560 -> ~580, 3-edge 5d3b 3,377 -> ~720, 3-edge 9cd9 1,661 ->
~370, 4-edge 1,780 -> ~405, 2-edge 1,114 -> ~255) and the pool-of-2 build wall 7,338 -> ~1,500 s
(the large coupons' meshio reads grew linearly with the element count, so the ratio is, if
anything, conservative there). Qualify: the reducer stages 12,592 -> ~1,480 s with the new
executable at b = 48 (step 4), i.e. the 8.26 node-h device-library run -> ~5.2 node-h; the workers
are unchanged. Every record of the two rebuilt roots is bit-identical apart from timings; the
physical gates, the quadrature rule / order and the rotate-z policy are untouched.

**Live acceptance of the merged qualify (PBS 46987, 2026-09-22; records under
`qualify/integration-20260922/`).** two-edge-8dd4bc70f183 from the step-3 AFTER build record
(identity 179a405d), `--reference none --controls p3`, the default executable 170439c4… (origin: the
manifest's `FrozenExecutable`) and block size 48 (origin: the manifest's `ReducerBlockSize`) on both
reducers: p4 worker 514.7 s + reducer 15.8 s (Palace peak 33.8 GB), p3 control worker 32.2 s +
reducer 10.4 s (22.3 GB), local-edge; 0.188 node-h; the job was held / re-queued eleven times by the
dispatcher (CF:ROLLBACK_COMPLETE, r8g.48xlarge capacity) before running. The verdict is `Failed` on
`PSequenceControls` alone: with the p3 control only there is no step to the higher order - the run
design of the smallest sufficient acceptance, not a physics result. Comparison target: **PBS 46730**
(the decision-61a MA-sharp acceptance: the same qualify inputs - the shelled mesh 179a405d, Tol
1e-10, the same eight controls - with b28f089a… at b = 6). PBS 46718 is not a valid target: it
reduced the PBS 46685 archives of the pre-shell mesh 5d01204e at Tol 1e-8 (no MA shell interfaces,
another solve tolerance), so bit-identity with it is impossible for any implementation (supervisor
decision, 2026-09-22). Result: the p3 control's 36 domain and 3,456 surface entries are **bit-identical**
to PBS 46730's (row order differs: (i <= j) vs block-pair order; every consumer is key-based); the
p4 stage 3,021 / 3,081 domain and 295,746 / 295,776 surface entries bit-identical, max per-entry
relative difference 9.7e-13 / 6.7e-13 (the CSV print resolution, as PBS 46718's 8.4e-13); the p4
local-edge output bit-identical including its row order (the worker / direct-solve path of the new
executable equals b28's). Reducer walls vs PBS 46730: p4 148.6 -> 15.8 s, p3 12.4 -> 10.4 s; Palace
peak 45.5 -> 33.8 GB and 29.9 -> 22.3 GB. One defect fixed forward on the way: `ma_tail.markdown`
formatted a `None` strongest-source median with `:.3f` for a run without a reference (the driver
stopped after the fetch / verification / archive deletion; `f3` prints n/a, test added, the analysis
completed under `--resume` without a second job).

## Preflight

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --preflight-only \
  --root /tmp/coupon-generality-preflight
```

This now succeeds for all 15 source cases and verifies the four-/ten-edge row
counts as four and ten. `--input CASE=DIRECTORY` cannot bypass a mismatched
hash. A case outside the recipe scope is recorded as `UnsupportedClass` (summary
`UnsupportedClassCases`, decision 48) and excluded from the build, not counted as a
preflight failure.

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
digest. The bounded record consumes the pipeline's separately executed reports
forming a digest-linked DAG: under the Gmsh-only production pipeline four
(canonical source validation, gmsh-build, final Gmsh publication, proper rigid
publication; the bounded record carries `Pipeline`); under the legacy MMG pipeline
seven (canonical source validation, seed generation, metric
preparation, native adapter/MMG, label restoration, final Gmsh publication, and
proper rigid publication). A tool is accepted only when its
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

Every proper-rigid placement uses the stage contract in `mesh_stage_contract.py`.
Under the Gmsh-only production pipeline the canonical build is three stages
(identity source validation, gmsh-build, canonical Gmsh publication; see above);
this section describes the legacy MMG pipeline's seven-stage contract, retired
from production by decision 38 and kept under the calibration manifest. Its first
six stages are one reusable source-local canonical build: identity source
validation, seed generation, metric preparation, native MMG adaptation, label
restoration, and canonical Gmsh publication. The seventh stage,
`proper-rigid-publication`, is mandatory even for identity (in both pipelines).

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
