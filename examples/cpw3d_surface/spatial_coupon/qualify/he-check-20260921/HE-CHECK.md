# Decision 53 — tube h_e check on two-edge 10 (ring-refined calibration variant), 2026-09-21

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

Question (decisions 52/53): is the two-edge p_MA strongest-20 excess at 53 / 58 (+2.8 / +3.5% vs the graded_v2 case-10 p5
anchor, verdict Failed under the frozen gates) the reference's under-resolution (our MA stable under h and p refinement)
or ours (our MA moves)? **Answer: our MA moves — not converged at the tube — and it moves AWAY from the reference under
both h (ring) and p refinement (E and p_MS do not move); the movement lives along the straight metal edges (the sharp-edge
MA tail), not in the corner balls (localization + model below).**

## Run (one job, `coupon_library.py qualify` at f46d34324, tool unchanged)
Build record `/tmp/coupon-calibration-tube-rings-20260920/library-build.json` (case `two-edge-calib-tube-rings-0.125nm`, mesh
`682b9117…`, 537,069 elements, labeled prism-condition deviation 1200); reference mirror `/tmp/library-acceptance-01-reference`
(case 10, p5 anchor, MA/MS); `--orders p4,p5 --controls p3,p5 --control-source 1 2 3 5 6 20 25 49 53 58 --max-jobs 1`,
binary `b28f089a…`, `--remote soca-green-job:/data/home/simlap/coupon_accuracy_assessment_20260913`, root `/tmp/library-he-check-01`
(dry run first: configs equal the acceptance p5 worker apart from paths and `Solver.Linear.Tol` 1e-10 = the decision-52 recipe
value; acceptance / gallery-10 ran 1e-8, a sub-4e-4% energy effect).
PBS 46094: qsub 02:18:33Z, R 02:27:20Z (r8g.48xlarge, 192 ranks, normal-g, DS-EM-FEM, efa), `job_state F`, `Exit_status 0`, walltime
00:58:04; 0 non-convergences; 16 CSVs hash-verified; matrices validated; p-sequence controls passed; archives deleted after
verification (16.4 GB, `Remaining ''`; re-verified read-only: 0 `archive/`, 22 MB left); read-only 90 s monitoring; job 46013
untouched; no qdel; no recorded campaign touched. Cost 3483 s = 0.968 node-h (p4 0.184, p5 0.636; production run 0.575 =
0.118 + 0.382): the ring set costs 1.67x at equal p, 0.050x the reference at p4.

## Verdict and class counts (free view, 78 sources, gated p5 vs the p5 anchor)
| run (mesh) | E <1% | p_MA median / weighted | p_MA <1/2/5% | strongest-20 <2% (worst) | p_MS <1/2/5% | verdict |
|---|---|---|---|---|---|---|
| production p4 (informational) | 76/78 | +0.10 / +0.29% | 56/76/78 | 18/20 (+2.67% at 53) | 73/77/78 | |
| production p5 (acceptance, PBS 46023) | 78/78 | +1.25 / +1.12% | 34/60/78 | 18/20 (+3.50% at 53) | 77/78/78 | Failed (53, 58) |
| gallery-10 p5 (recorded, recipe step aside) | 78/78 | +1.28 / +1.14% | 33/59/78 | 18/20 (+3.48% at 58) | 75/78/78 | Failed (53, 58) |
| **ring-refined p4** (informational) | 76/78 | +0.61 / +0.92% | 47/66/78 | 18/20 (+3.78% at 53) | 71/77/78 | |
| **ring-refined p5** (gated, PBS 46094) | **78/78** | **+1.72 / +1.65%** | **22/48/78** | **16/20 (+4.49% at 53)** | **77/78/78** | **Failed (53, 58, 69, 71)** |

## Ring (h) step at equal p, ring-refined vs production (same tool and recipe HEAD)
| quantity | p4: median all 78 / max | p5: median all 78 / max | p5 strongest-20 median / max | p5 at 53 / 58 |
|---|---|---|---|---|
| E | +0.00% / −0.17% (74) | +0.00% / +0.10% (65) | −0.00% / +0.06% | −0.01% / −0.02% |
| p_MA | **+0.67% / +1.37% (74)** | **+0.54% / +0.97% (58)** | **+0.47% / +0.97%** | **+0.96% / +0.97%** |
| p_MS | −0.02% / −2.00% (36) | +0.03% / +0.64% (36) | +0.03% / −0.12% | +0.03% / +0.10% |
p_MA at p5 moves upward at every one of the 78 sources; vs gallery-10 p5 the same (median +0.50%, 53 +0.98%, 58 +0.96%).
p-sequence vs the anchor — production 53: p4 +2.67 → p5 +3.50% (step +0.81%), 58: +2.66 → +3.46%; ring 53: p3 +2.99 → p4 +3.78 →
p5 +4.49% (steps +0.76 / +0.69%, Aitken r 0.91, unreliable), 58: +3.38 → +3.57 → +4.47%; strongest-20 p-step medians +0.73% /
+0.59% (production / ring); E steps −0.03…−0.04% (converged).

## Where 53 / 58 sit and where their MA lives (per-segment MA table of the p4 local-edge stage, manual edge 3100, R = 2 um)
Geometry: two collinear 1-um-wide strips at the metal-top level z = 0.1 — conductor 1 x ∈ [−5, −1], conductor 2 x ∈ [0, 4],
y ∈ [0, 1]; semantic 90° corners (balls r 0.1 um) at (−1, 0), (−1, 1), (0, 0), (0, 1); metal-edge / cut junctions (no ball)
at (−5, 0), (−5, 1), (4, 0), (4, 1). Source 53 = hat on the cut face x = 4 at (4, −0.5, 0.1): 0.5 um from the junction (4, 0)
along the y = 0 edge of conductor 2, 4.03 um from the nearest corner (0, 0), 5.0 um from conductor 1's end; 58 is its mirror
at (−5, 1.5, 0.1) (junction (−5, 1), corner (−1, 1)). The facing conductor ends x = −1 / x = 0 are 1.0 um apart vs
2 x (R_tube 0.032 + pyramid 0.008 + band 0.05–0.1) ≈ 0.2–0.3 um: the tubes / bands do not interact.
Ring mesh, p4 (616 segments, 18.0 um of edge): the MA of 53 is 99.2% on the single nearest edge line y = 0 (58: 99.2% on
y = 1), **80.3% straight edge, 19.7% junction ball (r < 0.1 um of (4, 0)), 0.00% corner balls**, 52% within 0.25 um and
81% within 0.5 um of the junction; "far surface" (> 2 um from every edge) is empty — the strips are 1 um wide, so the
table localizes along the edge only, not transversally. The production run has no per-segment table for 53 / 58 (not in
its control set) and no run has a p5 local-edge stage: the h-step at 53 / 58 and the p-step are not localizable directly.
Ring − production per shell at the 8 common controls (1, 2, 3, 5, 6, 20, 25, 49; p4; % of the production E_MA): total
+0.60% mean = **straight edge +0.51%**, corner balls +0.08% (their share 0–19%), junction balls +0.03%, far 0. Per 0.25-um
bin along the edges the step is uniform, +0.5–0.7% of each bin's energy (corner-ball bins ≈ +1.0%, junction-end bins at 25
+0.30–0.37%); the y = 0 line carries +0.55% of the mean +0.60%. **Reading: the +1% ring step at 53 / 58 and the +0.8%
p-step live along the straight edge next to the junction (0% of their MA is in a corner ball): the lever is the tube ring
set / tangential spacing along the edge, not the corner-ball grading** (their per-unit step, +0.96% of a 100%-edge MA, is
~1.7x the common controls' +0.5–0.6%: the junction-end field on the cut face has a heavier singular tail).

## MODEL (labeled; not a measurement): the sharp-edge MA tail
For a sharp 90° metal edge the exterior field scales as r^(−1/3), the MA integrand as r^(−2/3), the edge MA inside a radius
eps as eps^(1/3): the deficit D(eps) of an edge integral cut off at eps satisfies D(eps/2) = 2^(−1/3) D(eps) = 0.794 D(eps) —
each halving of the innermost ring recovers 1 − 2^(−1/3) = 0.206 of the remaining deficit, successive increments have
the constant ratio 0.79. Test on the two increments we have: two-edge production 0.25 nm → ring 0.125 nm: +0.5–0.6% of
the edge MA (straight-edge +0.51%, bins +0.5–0.7%; 53 / 58 +0.96%); four-edge 0.5 nm-tet reference (p4) → production
0.25 nm ring (p4): +0.62 / +0.79% at the edge-interior controls 80 / 23, median +0.82% (physics-13 / acceptance;
confounded by tet-vs-ring discretization at equal nominal size). Ratio 0.55 / 0.70 ≈ 0.79 — **consistent with the model
within the noise of two different coupons.** Implied remaining edge-MA deficit (hedged model extrapolation, I = 0.206 D):
four-edge chain D(0.5) ≈ 3.4%, D(0.25) ≈ 2.7%, D(0.125) ≈ 2.1%; two-edge chain D(0.25) ≈ 2.7%, D(0.125) ≈ 2.1% (the chains
agree); at 53 / 58 (I = 0.96%) D(0.25) ≈ 4.7%, D(0.125) ≈ 3.7%. A 0.5 nm-anchored p5 reference is then itself ≈ 3.4% low in
its edge MA and our 0.25 nm ring ≈ 2.7% low: expected equal-p gap ≈ +0.7% (observed all-78 median +1.25% at p5), larger
where the edge tail is heavier (53 / 58 +3.5%).

## Adjudication
1. **Not stable**: one ring halving moves p_MA at 53 / 58 by +0.96 / +0.97% at p5 (+1.08 / +0.88% at p4), the strongest-20
   by +0.47% and every free source by +0.54% (medians) — one p-step's worth. Decision-53 "moves" branch: our MA is not
   h-converged at the production ring set (≈ −1.0% p_MA per halving at 53 / 58, ≈ −0.5% at the median source; E, p_MS < 0.1%).
2. **Direction**: every refinement (h at the tube, p on either mesh) increases p_MA and widens the gap to the p5 anchor
   (53: +2.67 → +3.50 → +3.78 → +4.49%): the anchor is at least as under-resolved in the same direction; refining the ring
   alone cannot close the 2% statement (it widens it), and "our value is stable" (the reference-limited branch) is false.
3. **Conclusion**: the strongest-20 2% statement vs a 0.5 nm-anchored reference is not decidable by refinement — both
   sequences share the same slowly converging sharp-edge tail (0.79 per halving under the model). The two candidate
   resolutions — (i) analytic shell extrapolation of the edge remainder from per-ring MA labels (tube rings as label seams,
   MA per radial shell, Richardson in eps^(1/3)), or (ii) a rounded-edge process model with its own references — are
   physics / definition decisions for the supervisor / user, not mesh levers.
4. **Conditioning (requested)**: Palace kappa max 831 vs 416, min h 5.9e-6 vs 1.0e-5; PCG mean 21.6 / 20.4 (p4 / p5, max
   46 / 53) vs 11.9 / 10.8 (max 29 / 33); median avg reduction factor 0.327 / 0.294 vs 0.183 / 0.122 — Tol 1e-8 at the
   ring's own factors would take ~16.5 / 15.0 iterations, so Tol explains ~1.25x of the 1.8–1.9x and the ring's doubled
   prism condition ~1.5x at equal Tol; seconds per iteration unchanged (0.309 / 1.218 vs 0.318 / 1.260).

## Records (this root; no CSV committed; no code change)
`library-qualification.json`, `qualification-gates.json` (SHA `aea5d5f9…`), `process-library.json`, per-case qualification /
cost-summary / submission / result-csv-sha256 / matrix-validation / remote-archive-deletion JSON, `comparison/`, `results/main/`
(CSVs, logs, `qstat -xf`, p4 `surface-Q-edge-local.csv`), `he_check_compare.py` / `he_check_shells.py` + their JSON, `qualify.driver.log`.
HE_CHECK=MOVES;JOB=46094.ip-192-168-54-24.us-west-2.compute.internal;RECORD=/tmp/library-he-check-01/library-qualification.json
