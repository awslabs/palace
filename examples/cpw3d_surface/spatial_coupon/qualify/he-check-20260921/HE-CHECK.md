# Decision 53 — tube h_e check on two-edge 10 (ring-refined calibration variant), 2026-09-21

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

Question (decisions 52/53): is the two-edge p_MA strongest-20 excess at 53 / 58 (+2.8 / +3.5% vs the graded_v2
case-10 p5 anchor, verdict Failed under the frozen gates) the reference's under-resolution (our MA stable under
h and p refinement) or ours (our MA moves)? **Answer: our MA moves — it is not converged at the tube — and it
moves AWAY from the reference under both h (ring) and p refinement (E and p_MS do not move).**

## Run (one job, `coupon_library.py qualify` at f46d34324, tool unchanged)
Build record `/tmp/coupon-calibration-tube-rings-20260920/library-build.json` (case `two-edge-calib-tube-rings-0.125nm`,
identity mesh `682b9117…`, 537,069 elements, labeled prism-condition deviation 1200); reference mirror
`/tmp/library-acceptance-01-reference` (case 10, p5 anchor, MA/MS); `--orders p4,p5 --controls p3,p5`,
`--control-source 1 2 3 5 6 20 25 49 53 58` (the acceptance's class set + 53/58), `--max-jobs 1`, frozen binary
`b28f089a…`, `--remote soca-green-job:/data/home/simlap/coupon_accuracy_assessment_20260913`, local root
`/tmp/library-he-check-01` (dry run first at `…-01-dryrun`: configs equal the acceptance p5 worker apart from paths and
`Solver.Linear.Tol` 1e-10 = the recipe value of decision 52; the acceptance / gallery-10 ran the producer default 1e-8,
a sub-4e-4% energy effect per the manifest). PBS 46094: qsub 02:18:33Z, R 02:27:20Z on ip-192-168-37-130
(r8g.48xlarge, 192 ranks, normal-g, DS-EM-FEM, efa, one instance), `job_state F`, `Exit_status 0`, walltime 00:58:04;
0 non-convergences; 16 result CSVs hash-verified vs the remote; matrices validated (complete / symmetric / PSD);
p-sequence controls all passed. Remote archives deleted after verification (5.1G + 9.7G + 1.3G + 282M, `Remaining ''`,
verified read-only afterwards: 0 `archive/` directories, 22 MB left under `library-he-check-01/`). Monitoring
read-only at 90 s; the unrelated running job 46013 untouched; no qdel; no recorded campaign directory touched.
Cost: job 3483 s = 0.968 node-h (p4 stage 0.184, p5 0.636; production acceptance run 2071 s = 0.575 node-h,
p4 0.118, p5 0.382): the ring set costs 1.67x at equal p, 0.050x the reference's 3.68 node-h at p4.

## Verdict and class counts (free view, 78 sources, gated p5 vs the p5 anchor)
| run (mesh) | E <1% | p_MA median / weighted | p_MA <1/2/5% | strongest-20 <2% (worst) | p_MS <1/2/5% | verdict |
|---|---|---|---|---|---|---|
| production p4 (informational) | 76/78 | +0.10 / +0.29% | 56/76/78 | 18/20 (+2.67% at 53) | 73/77/78 | |
| production p5 (acceptance, PBS 46023) | 78/78 | +1.25 / +1.12% | 34/60/78 | 18/20 (+3.50% at 53) | 77/78/78 | Failed (53, 58) |
| gallery-10 p5 (recorded, recipe step aside) | 78/78 | +1.28 / +1.14% | 33/59/78 | 18/20 (+3.48% at 58) | 75/78/78 | Failed (53, 58) |
| **ring-refined p4** (informational) | 76/78 | +0.61 / +0.92% | 47/66/78 | 18/20 (+3.78% at 53) | 71/77/78 | |
| **ring-refined p5** (gated, PBS 46094) | **78/78** | **+1.72 / +1.65%** | **22/48/78** | **16/20 (+4.49% at 53)** | **77/78/78** | **Failed (53, 58, 69, 71)** |

## Ring (h) step at equal p: ring-refined vs production, same tool and recipe HEAD
| quantity | p4: median all 78 / max | p5: median all 78 / max | p5 strongest-20 median / max | p5 at 53 / 58 |
|---|---|---|---|---|
| E | +0.00% / −0.17% (74) | +0.00% / +0.10% (65) | −0.00% / +0.06% | −0.01% / −0.02% |
| p_MA | **+0.67% / +1.37% (74)** | **+0.54% / +0.97% (58)** | **+0.47% / +0.97%** | **+0.96% / +0.97%** |
| p_MS | −0.02% / −2.00% (36) | +0.03% / +0.64% (36) | +0.03% / −0.12% | +0.03% / +0.10% |
p_MA at p5 moves upward at every one of the 78 sources (31/78 within 0.5%, 78/78 within 1%); vs gallery-10 p5 the same
picture (median +0.50%, 53 +0.98%, 58 +0.96%). E and p_MS are ring-invariant: only the MA layer integral responds.

## p-sequence of p_MA at 53 / 58 (offset vs the p5 anchor; values in he-check-tables.json)
| mesh | p3 | p4 | p5 | p-step (p5−p4)/p5 | Aitken r (p3/p4/p5) |
|---|---|---|---|---|---|
| production, 53 | n/a (not a control) | +2.67% | +3.50% | +0.81% | n/a |
| production, 58 | n/a | +2.66% | +3.46% | +0.78% | n/a |
| gallery-10, 53 / 58 | n/a | +2.76 / +2.69% | +3.48 / +3.48% | +0.70 / +0.77% | n/a |
| ring-refined, 53 | +2.99% | +3.78% | +4.49% | +0.69% | 0.91 (p_inf +11%, unreliable at r→1) |
| ring-refined, 58 | +3.38% | +3.57% | +4.47% | +0.86% | 4.8 (non-geometric) |
Strongest-20 p-step medians: production +0.73%, ring +0.59%; all-78 medians +0.90% / +0.79%. On both meshes p_MA
rises monotonically with p; E steps are −0.03…−0.04% (converged), p_MS steps alternate within ±0.4%.

## Adjudication
1. **Not stable**: the ring halving (0.25 → 0.125 nm inner ring, one more ring, corners graded) moves p_MA at 53 / 58 by
   +0.96 / +0.97% at p5 (+1.08 / +0.88% at p4), the strongest-20 by +0.47% (median) and every free source by +0.54%
   (median) — the same size as one p-step (+0.7…0.9%). Under the decision-53 criterion this is the "moves" branch: our
   MA at the metal-edge junction is not h-converged at the production ring set, and the ring set is the lever (ring-set
   sensitivity: d p_MA / d log2(h_inner) ≈ −1.0% per halving at 53 / 58, ≈ −0.5% per halving at the median source;
   E and p_MS: < 0.1%).
2. **Direction**: every refinement we can apply (h at the tube, p on either mesh) increases p_MA and enlarges the gap
   to the p5 anchor (53: +2.67 → +3.50 → +3.78 → +4.49%). If the anchor were the converged value, refinement would
   move us toward it; it does not. The anchor is therefore at least as under-resolved as our production mesh, in the
   same direction (its own p-sequence is not available — the decision-53 caveat). The 2% strongest-20 statement at
   53 / 58 compares two unconverged MA values; refining the ring set alone cannot close it against this reference
   (it widens it), and a "reference-limited" annotation is NOT supported by this run's evidence either (our value is not
   stable). Neither branch's wording is recommended as is; the finding is reported for the supervisor's decision.
3. **Conditioning (requested)**: Palace mesh kappa max 831 vs 416 (production), min h 5.9e-6 vs 1.0e-5; PCG mean 21.6 /
   20.4 (p4 / p5, max 46 / 53) vs 11.9 / 10.8 (max 29 / 33); median avg reduction factor 0.327 / 0.294 vs 0.183 /
   0.122. At the ring run's own reduction factors, Tol 1e-8 would take ~16.5 / 15.0 iterations (the Tol difference
   explains ~1.25x of the 1.8–1.9x), so the ring doubles the prism condition and costs ~1.5x PCG iterations at equal
   Tol — a material but bounded degradation (s / iteration unchanged: 0.309 / 1.218 vs 0.318 / 1.260).

## Records (this root; no CSV committed)
`library-qualification.json`, `qualification-gates.json` (SHA `aea5d5f9…`: the acceptance's thresholds unchanged, plus the
3cd7a3aee SourceClasses / InterfaceTypes text blocks), `process-library.json`,
`two-edge-calib-tube-rings-0.125nm/{qualification,cost-summary,submission,result-csv-sha256,matrix-validation,
remote-archive-deletion}.json`, `comparison/` (per-source CSV/MD vs the anchor, p-sequence controls incl. 53 / 58),
`results/main/` (CSVs, logs, `qstat -xf`), `he_check_compare.py` + `he-check-tables.json` (this adjudication),
`qualify.driver.log`. Repository: no code change; this file copied under `qualify/he-check-20260921/` and a doc paragraph.

HE_CHECK=MOVES;JOB=46094.ip-192-168-54-24.us-west-2.compute.internal;RECORD=/tmp/library-he-check-01/library-qualification.json
