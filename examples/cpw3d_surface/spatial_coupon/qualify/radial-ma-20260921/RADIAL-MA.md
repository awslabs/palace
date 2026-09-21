# Decision 56 — radial MA profile of two-edge 10 (label-only shell relabel of the production mesh), 2026-09-21

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

Question (decisions 55 / 56): measure the sharp-edge MA tail directly instead of modelling it — MA per tube-ring shell per
source on the PRODUCTION two-edge mesh, the power law of the tail, the remainder inside the innermost ring, the implied
deficit of the production value. **Answer: the sharp 90° edge follows the r^(−2/3) law in every resolved inner ring at
every source (median local slopes −0.68 / −0.64 / −0.67 between rings 2–3 / 3–4 / 4–5, fitted alpha −0.648 (quartiles
−0.657 / −0.613, se 0.006), 53 / 58 −0.672 ± 0.011); the innermost ring (0–0.25 nm) holds only 0.706 (p5) / 0.665 (p4) of
the energy the −2/3 law anchored on ring 2 predicts — a source-independent factor (range 0.687–0.710 at p5) — so the
production MA is 2.7% (median, p5; quartiles 2.2 / 3.2%) and 4.4% at 53 / 58 below its sharp-edge limit, consistent with
the HE-CHECK model (2.7% / 4.7%) and reproducing the observed ring-halving step (0.206 × deficit: 0.50% median, 0.91% at
53 vs the measured +0.54% / +0.96%).**

## Run (one job, `coupon_library.py qualify`, root `/tmp/library-radial-ma-01`)
Mesh: the production identity mesh `5d01204e…` (PBS 46023's) relabeled by `relabel_radial_ma_shells.py` into 30 MA labels
(15 per conductor: far, top rings 1–7, bottom rings 1–7; label 10000 × ordinal + parent), SHA256 `179a405d…`, 28,360,844
bytes: $Nodes byte-identical, 602,862 elements identical apart from the (physical, elementary) pair of the 14,575 MA
elements (29,150 integers), shells sum to the parent areas 4.9 um² at 1.9e-14, parent ownership certificate reproduced
element by element, radial quadrature closure 0 (straddling 0.11%, corner-ball / cap triangles). Case `two-edge-calib-
radial-ma-shells` (sizing calibration manifest; `Calibration.Relabel`, labeled `PhysicsRun {LinearTol 1e-8}` =
PBS 46023's Tol so the reproduction is at the deterministic-solve level; supervisor option A), build record
`/tmp/coupon-calibration-radial-ma-20260921/library-build.json`. `qualify --reference /tmp/library-acceptance-01-reference
--orders p4,p5 --controls p3,p5 --control-source 1 2 3 5 6 20 25 49 53 58 --max-jobs 1`, binary `b28f089a…`, dry run first
(base config == the acceptance p5 worker config apart from paths / Order; 15 MA Dielectric entries 3–17 next to MS 2).
PBS 46164: qsub 07:20:55Z, R 07:23:26Z (r8g.48xlarge, 192 ranks, normal-g, DS-EM-FEM), `job_state F`, `Exit_status 0`,
walltime 00:39:27; 0 non-convergences; CSVs hash-verified; matrices validated (every shell PSD); archives deleted after
verification (4.8G + 9.2G + 1.2G + 267M, `Remaining ''`; re-verified read-only: 0 `archive/`, 52 MB left); 90 s read-only
polls; job 46013 untouched; no qdel. Cost 2366 s = 0.657 node-h (p4 0.131, p5 0.426; acceptance 0.575): the workers are
identical (322 / 1092 s vs 325 / 1095 s, PCG 11.90 / 10.82 = the acceptance's), the reducers +50 / +163 s for 16 vs 2
surface integrals.

## Label-only proof at the physics level (shells summed per source vs PBS 46023, same mesh geometry, same Tol)
p4: 78 sources, max |rel| E 1.1e-13, Q_MA 4.4e-13, Q_MS 1.1e-13; p5: E 0, Q_MA 3.6e-13, Q_MS 2.8e-13 — the run reproduces
the production acceptance run at 1e-10 per source (1e-12, in fact); the gated p5 verdict is the acceptance's to the digit
(Failed on p_MA: median +1.25%, weighted +1.12%, 34/60/78 within 1/2/5%, strongest-20 18/20, worst 53 +3.50%; E 78/78,
p_MS 77/78/78; p-sequence controls passed). Shell-sum identity: the 15 shells of each source sum to the whole MA by
construction (the surface integral is additive over the partition); the Tol effect (1e-8 vs the recipe's 1e-10, sub-4e-4%)
is immaterial for the per-shell profile (ring ratios of 1.5–2) and for the equal-p comparison with the ring run (0.5–1%).

## MA per shell (top edge = the sharp 90° edge at z = 0.1; bottom edge = metal / trench at z = 0; 78 free sources)
Shares (median, p5): top-edge rings 41.9%, bottom-edge rings 28.0%, far shell (> 31.75 nm from every edge line) 26.9%;
top ring 1 (0–0.25 nm) alone 5.7% of the MA (53 / 58: 10.9%), bottom ring 1 1.0%. Top-edge Q_k / Q_MA at 53, rings 1..7
(p5): 0.109 0.068 0.072 0.084 0.100 0.120 0.136 — ring 1 exceeds ring 2 (as a −2/3 law demands: Q_1 / Q_2 = 2.26 for
the law, 1.60 resolved). Local slopes d log(Q_k / h_k) / d log r (median over sources), rings 1–2 .. 6–7, p5: top −0.84,
−0.68, −0.64, −0.67, −0.62, −0.57; bottom −0.33, −0.28, −0.34, −0.37, −0.38, −0.42. The p4 profiles are the same to
< 0.01 in every slope but the first.

| estimator (top edge unless noted) | alpha median (quartiles) | se | strongest-20 | 53 / 58 | deficit median (q) p5 | strongest-20 | 53 / 58 | p4 median / 53 |
|---|---|---:|---:|---:|---|---:|---:|---|
| Fit2-4 (rings 2–4, 0.25–3.75 nm) | −0.648 (−0.657 / −0.613) | 0.006 | −0.645 | −0.672 | top part +1.01%; total +1.89% (+0.96 / +2.85%) | +1.76% | +4.75% | +2.21% / +5.40% |
| Fit2-K (rings 2–7) | −0.622 (−0.654 / −0.593) | 0.012 | −0.615 | −0.704 ± 0.011 | total +1.97% (+1.35 / +3.04%) | +1.80% | +7.26% | +2.31% / +7.95% |
| Theory@2 (−2/3 anchored on ring 2) | −2/3 | – | – | – | top part +2.41%; total (−2/3 on both edges) +4.44% | +4.44% | +4.52% | +2.76% / +5.20% |
| bottom edge Fit2-4 (45 sources with bottom share > 20%) | −0.332 (−0.392 / −0.301) | 0.013 | | 53 / 58: no power law (3% share) | bottom part (own law) +0.10% median | | | |
| **consistent: top −2/3 @ ring 2 + bottom own law** | | | | | **+2.68% (+2.22 / +3.21%)** | **+2.65%** | **+4.44%** | **+3.07% / +5.11%** |

Ring 1 resolved / model (−2/3 anchored on ring 2): p5 median 0.706 (0.687–0.710 over all 78 sources), p4 0.665 (0.648–
0.668) — the innermost ring's resolution factor is a property of the discretization (0.25 nm inner ring, p), not of the
source. p-step p4 → p5 per ring (median relative change of Q_k): top rings 1..7 +6.2%, +0.06%, +0.01%, +0.13%, +0.08%,
+0.07%, −0.03%; bottom +1.4% then ≤ 0.07%; far +0.9%; Q_MA +0.79% (E −0.07%). At 53 / 58 the +0.76% MA step is 84% top
ring 1 and 20% far shell (rings 2–7 −5%); median source 40% ring 1 / 46% far. Rings 2–7 are p-converged; the far shell
(corner balls, junction, top-face interior) carries a p-dependence of its own, outside the tube.

Strongest-20 detail (p5; class; shares top / bottom / far %; slopes 2–3, 3–4, 4–5; alpha Fit2-4; Fit2-K (se); clean; ring-1
resolved / model; deficit top −2/3; consistent total):
78 terminal 43/25/32 −0.66 −0.53 −0.71 −0.594 −0.623 (0.009) yes 0.697 +2.49% +2.63% · 77 junction column 34/1/65 −0.67 −0.59 −0.69
−0.620 −0.651 (0.008) yes 0.700 +2.05% +2.05% · **58 / 53 near-junction hat 69/3/28 −0.69 −0.68 −0.69 −0.672 −0.704 (0.011) yes 0.710
+4.44% +4.44%** · 29 / 38 junction ring 18/64/18 −0.68 −0.65 −0.60 −0.651 −0.554 (0.036) no 0.707 +0.93% +3.25 / +3.34% · 25 / 35
junction ring 11/75/14 −0.68 −0.64 −0.58 −0.645 −0.517 (0.047) no 0.706 +0.52% +3.40% · 76 near-junction hat 5/0/94 −0.69 −0.67 −0.65
−0.662 −0.633 yes 0.708 +0.32% · 69 box corner 43/2/55 −0.68 −0.62 −0.69 −0.635 −0.664 yes 0.702 +2.66% · 19 junction column 37/39/25
−0.66 −0.49 −0.69 −0.568 −0.590 yes 0.693 +2.04% +2.25% · 68 junction column 34/1/65 −0.68 −0.64 −0.67 −0.646 −0.650 (0.002) yes 0.704
+2.03% · 48 / 43 near-junction hat 34/26/40 −0.63 −0.51 −0.32 −0.562 −0.306 (0.073) no 0.691 +1.08% +1.05% · 24 / 34 near-junction hat
37/39/24 −0.68 −0.66 −0.63 −0.658 −0.602 no 0.708 +2.04% +2.21% · 63 near-junction hat 40/2/58 −0.68 −0.64 −0.68 −0.646 −0.664 yes
0.705 +2.46% · 71 wide hat 44/3/53 −0.68 −0.63 −0.67 −0.645 −0.665 yes 0.705 +2.72% · 31 / 21 wide hat 45/28/27 −0.66 −0.52 −0.68
−0.586 −0.616 yes 0.695 +2.62% +2.71 / +2.73%.

Clean power law (Fit2-K residual RMS < 0.05 and every fitted local slope within 0.25 of alpha): top edge 56 / 78 (p5; p4
58 / 78), bottom edge 26 / 78. The 22 non-clean top profiles are 13 near-junction hats, 4 junction rings, 4 wide hats
next to the cut faces and 1 junction column: their inner slopes 2–3 / 3–4 are still −0.68 / −0.64 … −0.51 (49 sources
within 0.05 of −2/3 on both, 56 within 0.10), the departure is in rings 5–7 (4–32 nm) where the profile bends — flatter
(−0.3 … +0.14 at 43 / 48, the near-junction hats on the cut faces x = ±4 / −5) or steeper (−0.73 / −0.80 at 53 / 58) as the
finite geometry (0.1 um metal thickness, the junction with the cut face) takes over from the edge asymptotics. The
bottom edge is not a −2/3 edge: alpha −0.33 (E ~ r^(−0.17); the conductor corner on the substrate / trench interface),
its ring 1 moves +1.4% with p, its remainder is 0.1% of the MA. The Fit2-K alpha (−0.62 median, −0.70 at 53 / 58) is a
mixture of the edge exponent and the outer bend and gives the least stable extrapolation (deficit 1.4–7.3% between
sources with the same inner profile); Fit2-4 and Theory@2 agree on the inner law.

## What the data says (evidence only, no decision)
1. **The sharp-edge tail is real and universal**: r^(−2/3) holds in every source's inner rings (0.25–3.75 nm) including
   53 / 58 and the junction-ring sources; the resolved ring-1 energy is 0.71 of the law's at p5 (0.67 at p4) at every
   source. The production MA (0.25 nm inner ring, p5) is therefore ~2.4% (top edge, median; 53 / 58 4.4%) below the
   sharp-edge limit; with the bottom edge's own law ~2.7% (quartiles 2.2 / 3.2%; strongest-20 2.7%). The HE-CHECK model's
   2.7% / 4.7% at 0.25 nm is confirmed by measurement, and the model's ring-halving prediction (0.206 × deficit = 0.50% /
   0.91%) matches the ring-refined run's +0.54% / +0.96%.
2. **(i) Analytic extrapolation** is supported: the per-ring MA labels make the remainder computable per source from the
   run itself (Theory@2: −2/3 anchored on the nearest resolved ring; no fit needed), the correction factor is
   source-independent (0.706 ± 0.01 at p5) so a per-ring rule would be one number per (ring set, order), and the outer
   rings are p-converged. Caveats the data shows: the far shell carries 20–46% of the remaining p-step (corner balls /
   junction / top-face interior, outside the tube), so the edge extrapolation removes the ring-1 deficit but not the
   whole p-dependence; the extrapolation is only as good as the −2/3 assumption inside 0.25 nm (the fitted inner exponent
   −0.648 ± 0.006 vs −0.667 changes the ring-1 model by ~5%, the deficit by ~0.5 points); the reference (0.5 nm tets)
   would carry its own, larger deficit (~5.6% at 53 / 58 under the same law), so extrapolating both sides is required
   before any 2% statement.
3. **(ii) Rounded edges**: the measured profile says where a rounding radius rho acts — 5.7% (median) / 10.9% (53 / 58)
   of the MA sits inside 0.25 nm and 9.3% / 17.6% inside 0.75 nm, so any rho of nm order changes the MA by several
   percent and the value becomes a process quantity; the −2/3 law would then be cut off physically rather than
   numerically, and the far-shell p-dependence remains either way. The data does not favour one option; it quantifies
   both: the tail below 0.25 nm is 2.4–4.4% of the MA and its shape is known.

## Records
`/tmp/library-radial-ma-01`: `library-qualification.json`, `qualification-gates.json`, `process-library.json`, per-case
records (submission / result-csv-sha256 / matrix-validation / remote-archive-deletion JSON, `comparison/`, `results/main/`
CSVs incl. the per-shell `surface-response-matrix.csv` of every stage, logs, `qstat -xf`), `radial-ma-profile.{json,md}`
(`qualify/radial_ma_profile.py`), `qualify.driver.log`; mesh + census + build record under
`/tmp/coupon-calibration-radial-ma-20260921`; dry run `/tmp/library-radial-ma-01-dryrun`.
RADIAL=-0.648;2.68%;JOB=46164.ip-192-168-54-24.us-west-2.compute.internal
