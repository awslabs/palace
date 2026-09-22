# Decision 61a acceptance — the sharp-edge MA on the production four-edge and two-edge 10 coupons, 2026-09-22

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

Question (user decision 60(1), supervisor decision 61a/d): with the per-ring MA shells on every production coupon and
`qualify` computing MA_sharp = MA_raw + MA_tail per source, re-evaluate the four-edge and two-edge 10 records - MA_raw must
equal the acceptance runs (PBS 46024 / 46023) to roundoff (labels only), and the strongest-20 statement is re-read on the
sharp values against the 0.5 nm graded_v2 references (unextrapolated, modelled by the eps^(1/3) law). **Answer: the
shells change nothing but labels - four-edge MA_raw = PBS 46024 at 8.5e-13 (E 3.3e-13, Q_MS 5.3e-13, 80 sources); the
two-edge production mesh is byte for byte the decision-56 relabel (`179a405d…`) whose run 46164 reproduced PBS 46023 at
4.4e-13, and this run (recipe Tol 1e-10 vs 46023's 1e-8) agrees with 46023 at 3.3e-4 = the Tol effect. Four-edge
PASSES on MA_sharp (free median +0.14% vs +0.79% raw; 44 / 53 / 60 within 1 / 2 / 5% vs 30 / 53 / 60; strongest-20
median +0.18%, 20 / 20 within 1%); two-edge 10 still FAILS on the strongest-20 at 53 / 58: +2.37% / +2.33% sharp-vs-
modelled (raw +3.50% / +3.46%), the median at +0.55% (raw +1.25%), 51 / 76 / 78 within 1 / 2 / 5% (raw 34 / 60 / 78).
The 53 / 58 excess is NOT explained by the estimator uncertainty: Fit2-4 / Theory@2 / Consistent give +2.30 / +2.35 /
+2.37% (spread 0.07 points); the modelled 0.5 nm reference deficit at 53 is 5.59% and it would have to be 5.98% for the
2% statement (8.1% for equality) - 0.4 points more than the eps^(1/3) scaling of the run's own deficit gives, i.e. the
tet reference resolving its innermost 0.5 nm somewhat worse than a prism ring of the same size does (its ring-1 factor
is not measured). Evidence only; no gate or threshold changed.**

## Production rebuild (decision 61a step 1; root `/tmp/coupon-library-build-73dd25a22-20260922-013409`, jobs 2)
Every production coupon is published with the per-ring radial MA shells (`publish_rigid_coupon_mesh.apply_radial_shells`
at the placement stage, label-only, census bound in the transform receipt; the Julia ownership auditor and every Python
audit run on the shelled mesh through the family rule / parent view). Rebuilt and VERIFIED (per-entry verification
passed, 7 rings, inner ring 0.25 nm); `verify_label_only_republication.py` against the previous production / newest
local identity: $Nodes byte-identical, every element identical apart from the (physical, elementary) pair of the MA
elements, the receipt's ParentLabeledMeshSHA256 = the previous identity (records under `label-only/`):

| coupon | previous identity (= parent-labeled) | shelled identity | shells (parents) | relabeled / elements | straddling | parent-area closure |
|---|---|---|---|---:|---:|---:|
| two-edge 10 (8dd4bc70f183) | `5d01204e` (PBS 46023's) | `179a405d` (= the decision-56 relabel) | 30 (6001, 6002) | 14,575 / 602,862 | 1.1e-3 | 3e-14 |
| two-edge 05 (3f8992613e95) | `9f011ed3` (Sep-19 root) | `86cfad22` | 15 (6001) | 22,256 / 1,725,520 | 1.6e-5 | 4e-14 |
| four-edge (9d2cb9bbb3fe) | `1d536c44` (PBS 46024's) | `f36d2d83` | 15 (6001) | 27,033 / 2,169,612 | 1.7e-4 | 1e-13 |
| ten-edge (6791f1c84123) | `8a871f22` (Sep-19 root) | `e1a68075` | 60 (6001, 6002, 6101, 6102) | 59,209 / 4,101,387 | 1.3e-4 | 3e-14 |
| three-edge 06 (419576fdab24) | `95837ed5` (Sep-19 root) | `4e2adbe1` | 15 (6001) | 34,377 / 2,126,224 | 6.4e-5 | 7e-14 |

Three-edge 06 failed its first publication: `relabel_radial_ma_shells.metal_edge_lines` consumed (`set.pop`) the
conductor's process-normal set shared by every plan-view loop of that conductor, so the second loop of conductor 1 read
an empty set. Fixed forward (commit 79df1f416, test `test_metal_edge_lines_of_two_loops_of_one_conductor`, tool hash
refrozen) and rebuilt in `/tmp/coupon-library-build-79df1f416-20260922-023049` (`library-build/library-build-root2-three-edge.json`).
The two-edge 05 / 10 censuses record the pre-fix tool hash `98ce14c1…` (published before the fix; the fix does not change
the output where the old code did not fail - two-edge 10 is byte-identical to the decision-56 mesh). Build walls: 563 /
1523 / 2281 / 4012 s (two-edge 10 / 05 / four-edge / ten-edge), 1982 s (06); the mesher is deterministic (every
parent-labeled publication equals the previous identity byte for byte).

## Run (`coupon_library.py qualify`, root `/tmp/library-ma-sharp-01`, remote `coupon_accuracy_assessment_20260913/library-ma-sharp-01`)
`--build-record <root 1>/library-build.json --case four-edge-9d2cb9bbb3fe --case two-edge-8dd4bc70f183 --reference
/tmp/library-acceptance-01-reference --frozen-binary-sha256 b28f089a… --max-jobs 2 --monitor-interval 90
--reference-edge-size-nm 0.5`; dry run first (`/tmp/library-ma-sharp-01-dryrun`: controls = the acceptance's by class,
base config == the reference config apart from paths / Order / Tol, 15 MA shell entries per config, two-edge adds the
p5 main). Policy frugal, N = 1 per coupon (one instance each, 192 ranks, r8g.48xlarge, normal-g, DS-EM-FEM):
four-edge PBS 46729 (R 03:15:15Z, F Exit 0, walltime 00:49:28, host 192.168.50.255; p4 worker 1379 s / reducer 562 s,
PCG 20.05 mean / 40 max; controls p5 532 + 85, p3 87 + 42; 0.824 node-h), two-edge 10 PBS 46730 (R 03:15:16Z, F Exit
0, walltime 00:56:21, host 192.168.51.62; p4 508 + 149 s, p5 1875 + 443 s; 0.939 node-h). Total 1.763 node-h, critical
path 3631 s; CSVs hash-verified, matrices validated (every shell PSD), archives deleted after verification (0
`archive/` left, 140 MB remaining, re-verified read-only); 90 s read-only polls; jobs 46013 / 46251 / 46350 / 46716 of
other work untouched; no qdel.

## MA_raw vs the acceptance runs (shells summed per source; `radial-ma-profile.json` Reproduction)
Four-edge vs PBS 46024 (p4, 80 sources, both Tol 1e-10): max |rel| E 3.3e-13, Q_MA 8.5e-13, Q_MS 5.3e-13 - roundoff.
Two-edge 10 vs PBS 46023 (78 sources): E 5.0e-12 / 1.2e-11 (p4 / p5) but Q_MA 2.5e-4 / 3.3e-4 and Q_MS 3.8e-4 / 1.5e-4:
46023 ran the producer default Tol 1e-8 (the acceptance tool was fixed forward after its two-edge job), this run the
recipe's 1e-10 - the CG stopping point differs, sub-4e-4% on the surface integrals (RADIAL-MA.md's Tol statement). The
roundoff identity for two-edge holds through the decision-56 run: this production identity `179a405d…` IS the mesh of
PBS 46164 (Tol 1e-8), which reproduced 46023 at E 0 / Q_MA 3.6e-13 / Q_MS 2.8e-13 (p5).

## MA_sharp (`comparison/ma-tail.{json,md}`, `Qualification.MASharp`; the `Consistent` estimator, ring-1 factor measured on each run)
| coupon / order | alpha top rings 2-4 median (quartiles; se) | ring-1 factor median (range) | deficit Fit2-4 | Theory@2 | **Consistent** (quartiles) | strongest-20 | modelled reference deficit (0.5 nm) |
|---|---|---|---:|---:|---|---:|---:|
| four-edge p4 (gated) | -0.655 (-0.661 / -0.639; 0.003) | 0.666 (0.644-0.668) | +2.27% | +4.68% | **+2.86%** (+2.5 / +3.3%) | +2.14% | 3.64% median |
| two-edge 10 p4 | -0.648 (-0.657 / -0.613; 0.006) | 0.665 (0.648-0.668) | +2.21% | +4.88% | **+3.07%** | +3.02% | - |
| two-edge 10 p5 (gated) | -0.648 (-0.658 / -0.609; 0.006) | 0.706 (0.687-0.710) | +1.89% | +4.44% | **+2.68%** (+2.2 / +3.2%) | +2.65% | 3.38% median |

The two-edge numbers equal RADIAL-MA.md's to the digit (same mesh, Tol immaterial for the profile; 53 / 58: alpha
-0.672 ± 0.002, ring-1 0.710, deficit Consistent 4.44%, Fit2-4 4.75%, Theory@2 4.52%). The four-edge p4 ring-1 factor
0.666 equals the two-edge p4 0.665 - the factor is a property of the discretization (0.25 nm inner ring, p), not of the
coupon; the four-edge alpha -0.655 is 0.007 steeper. p-sequence on p_MA_sharp (|d45| < 5%): four-edge max +1.67%
(source 4), two-edge max +1.52% (source 3); every control passed on both coupons.

## Gates (p_MA on p_MA_sharp vs the modelled reference; `qualification.json`)
| coupon | verdict | free median raw -> sharp | within 1 / 2 / 5% raw -> sharp | weighted raw -> sharp | strongest-20 median raw -> sharp | strongest within 2% | worst |
|---|---|---|---|---|---|---|---|
| four-edge (vs p4) | **Passed** | +0.79% -> +0.14% | 30 / 53 / 60 -> 44 / 53 / 60 | - -> +0.24% | +0.79% -> +0.18% | 20 / 20 (20 / 20 within 1%) | 74: +3.92% -> +3.29% (free bound 5%) |
| two-edge 10 (vs p5) | **Failed** (p_MA strongest) | +1.25% -> +0.55% | 34 / 60 / 78 -> 51 / 76 / 78 | +1.12% -> +0.42% | +0.98% -> +0.46% | 18 / 20 (53, 58 failing) | 53: +3.50% -> +2.37%, 58: +3.46% -> +2.33% |

Strongest-20 raw-vs-raw and sharp-vs-modelled (`ma-ms-offsets.json` `p_MA_rel:strongest20` / `p_MA_sharp_rel:strongest20`):
four-edge raw median +0.79% (range -0.26 … +1.54%, 14 / 20 within 1%), sharp +0.18% (-0.51 … +0.81%, 20 / 20 within 1%);
two-edge raw +0.98% (-0.06 … +3.50%, 10 within 1%), sharp +0.46% (-0.72 … +2.37%, 15 within 1%).

## Is the 53 / 58 excess explained within the estimator uncertainty? (evidence only)
No. At 53 the raw excess +3.50% becomes +2.37% with `Consistent` (run deficit 4.44%, modelled reference deficit 4.44% x
2^(1/3) = 5.59%), +2.30% with Fit2-4 on both sides (4.75% / 5.98%) and +2.35% with Theory@2 (4.52% / 5.69%): the
estimator spread moves the sharp excess by 0.07 points, all above the 2% bound. The 2% statement would need a reference
deficit of 5.98% (equality 8.1%) against the modelled 5.59% - the eps^(1/3) law scales the RUN's ring-1 deficit to the
reference's 0.5 nm edge size and thereby assumes the tet reference resolves its innermost 0.5 nm as well as a prism ring
does (the same ring-1 factor); a tet band of nominal size 0.5 nm resolving the r^(-2/3) singularity 7% worse would close
the gap. The near-junction hats 53 / 58 carry the largest top-edge share (69%) and ring-1 share (10.9%) of any source,
so they are the most sensitive to the reference's own edge resolution - which is not recorded and cannot be
extrapolated; the excess stays classified reference-limited (decision 55's annotation), now at 2.3-2.4% instead of 3.5%.
The four-edge worst free source 74 (+3.92% raw -> +3.29% sharp; deficit Consistent 2.42%, Fit2-4 -0.19%: its bottom
edge dominates and has no -2/3 law) is under the free 5% bound in both readings.

## Records
`library-qualification.json`, `process-library.json` (the four-edge model LibraryQualified with `MA {MA_raw, MA_sharp,
p_MA_raw, p_MA_sharp}` per source; the two-edge model not qualified), `qualification-gates.json` (the table used, with
`p_MA.Quantity` / `PSequenceControls.MAObservable` = p_MA_sharp), per coupon `qualification.json`, `cost-summary.json`,
`submission.json`, `remote-archive-deletion.json`, `result-csv-sha256.json`, `matrix-validation.json`,
`radial-ma-profile.{json,md}` (the profile and the Reproduction vs the acceptance run), `comparison/{ma-tail.json,
ma-tail.md, ma-ms-offsets.json, ma-ms-offsets.md, p-sequence-controls.json, key-sources.json}`; `label-only/` (the five
labels-only records), `library-build/` (the two library-build records), `qualify.driver.log`, `qstat-xf-46729-46730.txt`. The
per-shell response CSVs stay under `/tmp/library-ma-sharp-01/<case>/results/main/<stage>/reducer/`.
