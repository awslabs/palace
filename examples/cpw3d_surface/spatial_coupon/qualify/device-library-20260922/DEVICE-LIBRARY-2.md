# Device library end to end with the integrated tooling (user decision 64a, 2026-09-22)

The same transmon layout as the first complete run (decision 58, 2026-09-21: `examples/transmon/transmon_surface_coarse.json`
with the interface-layer seed copy, Delaunay device basis) rebuilt and re-qualified with the decision-62 / 63 tooling: labels-only
registration probes, the audit dedupe / caching, per-ring MA shells at publication (MA_sharp), the streaming one-pass Gram
executable `170439c4…` at reducer block size 48 and the decision-61b source split (`--job-policy speed --max-jobs 6`). Build tree
`307086ae7` (root `/tmp/coupon-device-transmon-20260922`, manifest SHA `4173b409…`), qualify `ToolCommit 7c7bb5fe5` (root
`/tmp/library-device-transmon-03`, remote `coupon_accuracy_assessment_20260913/library-device-transmon-03`, build record SHA
`378879f4…`, cost model `22ee226d…` = the physics-11 / b28 model the jobs were planned with). Records: `qualify/device-library-20260922/`.

## Build (local, 6 cores, pool of 3; no cluster node-hours)

Discovery 08:42:07 -> 08:42:37Z (30 s; 3,120 requirements, 35,297 um, 9 planned coupons, 6 spatial). Registration of the six coupons
with the labels-only probes, in parallel: 08:42:39 -> 08:43:12Z = **33 s** (was 990 s). Pool wall **1,462 s** (was 7,338 s); end to
end 08:42:07 -> 09:07:34Z = **1,527 s** (was 8,493 s). Every case id, element count and H1 count identical to 2026-09-21.
**Label-only check** (`label-only-vs-20260921/<case>.json`, 6/6 `labels-only`): each 2026-09-22 identity mesh equals the 2026-09-21
Delaunay identity in the `$Nodes` block and in every element apart from the (physical, elementary) pair of the MA elements moved
into the ring shells; the publication receipt's `ParentLabeledMeshSHA256` is the 2026-09-21 identity digest for every coupon.

| coupon (sources) | elements | H1 p4 | gmsh-build s | publication id / rot s | audits id / rot s | verification s | case wall s 09-21 -> 09-22 | MA elements relabeled |
|---|---|---|---|---|---|---|---|---|
| spatial-10-edge-65450ff47b9b (191) | 2,643,905 | 35,861,283 | 246 | 70 / 71 | 107 / 108 | 161 | 2,659 -> 814 | 59,237 (parents 6001/6002/6101/6102) |
| spatial-2-edge-8ce3fb213a2d (150) | 1,215,960 | 17,099,910 | 82 | 26 / 28 | 43 / 43 | 69 | 1,118 -> 322 | 22,218 |
| spatial-3-edge-9cd906b81cbe (175) | 1,680,987 | 22,886,084 | 122 | 40 / 41 | 59 / 62 | 110 | 1,665 -> 471 | 34,414 |
| spatial-3-edge-5d3b5e644745 (225) | 3,037,912 | 42,696,378 | 254 | 63 / 64 | 104 / 116 | 177 | 3,383 -> 835 | 62,207 |
| spatial-4-edge-e51d7380245e (120) | 1,895,487 | 25,167,936 | 160 | 37 / 34 | 60 / 59 | 113 | 1,784 -> 501 | 27,028 |
| spatial-5-edge-9bfa8265e2e7 (185) | 2,696,939 | 35,941,049 | 198 | 53 / 54 | 77 / 83 | 140 | 2,566 -> 648 | 62,410 (6001/6101) |

Stage sums over the six, 2026-09-21 -> 2026-09-22 (s): gmsh-build 832 -> 1,062 and canonical-publish 193 -> 221 (unchanged code; a
pool of 3 instead of 2 on the same 6 cores), identity publication 740 -> 291, rotate-z publication 743 -> 293, identity audits 1,593 ->
450, rotate-z audits 2,544 -> 471, per-entry verification 6,500 -> 771; serial sum of the case walls 13,175 -> 3,591 s (0.27; the
decision-62 projection from the two-edge ratios was ~2,940 s). 0 failures, 0 headroom flags (bounds 6M / 3,600 s / 12 GiB; the 7f03
verification that sat at 0.94 of the old 1,800 s bound took 177 s). Out of scope as before: corners 30 -> CornerCoupon, isolated
edges 3,082 -> StraightEdgeBuilder.

## Qualify (PBS 47080-47467, one r8g.48xlarge / 192 ranks per job, up to 6 concurrent; `library-qualification.json`)

Six coupons under `--job-policy speed --max-jobs 6`: N contiguous source blocks as N worker jobs (job 1 = the p3 / p5 controls and the
local-edge stage alone), one reducer job per coupon on the archive union. 7f03, fail-closed as one job on 2026-09-21 (29,272 s), ran as
6 + 1 jobs. p-sequence gate per control source (8 per coupon): |d45| of E within 1%, of every participation within 5%.

| coupon | PBS | H1 p4 | src | jobs W+R (block sizes) | node-h job / main | PCG mean / max (s per it) | worker s (sum of blocks) / reducer s | critical path s | max E abs d45 | max participation abs d45 | MA_sharp deficit median [quartiles] (ring-1 factor) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10-edge 65450ff47b9b | 47080-47108 | 35,861,283 | 191 | 5+1 (0/48/48/48/47) | 1.37 / 1.03 | 14.7 / 40 (1.15) | 3,565 / 133 | 6,661 | 0.19% (src 4) | p_MA 1.22% (src 126) | 2.39% [0.24, 2.46] (0.666) | PendingQualification |
| 2-edge 8ce3fb213a2d | 47086-47218 | 17,099,910 | 150 | 4+1 (0/50/50/50) | 0.69 / 0.49 | 18.4 / 39 (0.57) | 1,701 / 52 | 10,211 | 0.10% (src 7) | p_MS 1.63% (src 33) | 3.33% [2.89, 3.55] (0.667) | PendingQualification |
| 3-edge 9cd906b81cbe | 47109-47289 | 22,886,084 | 175 | 5+1 (0/44/44/44/43) | 0.94 / 0.72 | 18.0 / 40 (0.74) | 2,522 / 78 | 12,086 | 0.32% (src 2) | p_MA 2.68% (src 2) | 2.63% [2.28, 3.35] (0.667) | PendingQualification |
| 3-edge 5d3b5e644745 (7f03) | 47214-47300 | 42,696,378 | 225 | 6+1 (0/45/45/45/45/45) | 2.00 / 1.58 | 17.7 / 40 (1.28) | 5,527 / 158 | 8,945 | 0.19% (src 6) | p_SA 1.02% (src 46) | 2.51% [2.11, 3.23] (0.667) | PendingQualification |
| 4-edge e51d7380245e | 47290-47358 | 25,167,936 | 120 | 4+1 (0/40/40/40) | 0.86 / 0.59 | 18.3 / 40 (0.86) | 2,049 / 71 | 9,469 | 0.49% (src 15) | p_SA 1.21% (src 15) | 3.15% [2.83, 3.39] (0.667) | PendingQualification |
| 5-edge 9bfa8265e2e7 | 47335-47467 | 35,941,049 | 185 | 5+1 (0/47/46/46/46) | 1.34 / 0.98 | 14.2 / 40 (1.18) | 3,383 / 135 | 11,843 | 0.09% (src 1) | p_MS 1.69% (src 38) | 0.01% [0.00, 0.40] (0.666) | PendingQualification |

Every p-sequence gate passed (E <= 0.49% of 1%, participations <= 2.68% of 5% - the same maxima and sources as 2026-09-21: the worker
path is unchanged); verdict `PendingQualification` on all six ("no reference matrices: only the p-sequence controls were evaluated
(passed); never Passed") - **0 Failed, 0 Passed, nothing LibraryQualified**; `process-library.json` carries the six models with their
fabricated matrices (18,336 / 605,088 domain / surface rows for the ten-edge, 36 / 1,188 per control). PCG counts per coupon equal the
2026-09-21 ones. MA_sharp (decision 61a: MA_raw + the `Consistent` tail, ring-1 factor measured per run at 0.666-0.667 on every coupon):
the tail is 2.4-3.3% of MA_raw at the median source of five coupons and 0.01% [0.00, 0.40] on the 5-edge (tail summary over 119 of its 185
sources: the tail is negligible against MA_raw for most of them - as measured, not interpreted here). 7f03 is new: 2.00 node-h, 1.58 of them in the p4 main stage.
Every stage `complete`, every matrix validated, 18-20 result digests per coupon verified, every archive (51 / 20 / 30 / 71 / 23 / 50 GB
p4 + controls) deleted on the cluster after the fetch (`Remaining` empty); no qdel, no other job touched.

**Cost.** Job seconds 4,939 / 2,471 / 3,397 / 7,188 / 3,104 / 4,832 = **25,932 s = 7.203 node-h** (main stages 5.38); the five coupons
of 2026-09-21 alone 5.21 node-h vs 8.257 (the decision-62 projection was ~5.2). **Reducer share**: 626 s = **2.4%** of the job seconds
(2026-09-21: 11,964 s = 40.3%; per coupon 3,972 -> 133 s ten-edge, 994 -> 52 s two-edge). **Estimate vs actual** (planned with the
physics-11 / b28 model at b = 48): every job's 2.0x-with-margin estimate 2.5-4.4x the worker job and 8.2-11.6x the reducer job; the
per-coupon node-time estimate at 2.0x 7,832-24,226 s vs 2,471-7,188 s actual. **Jobs**: 35 submitted (5 x 5+1, 4+1 x 2, 6+1), never
more than 6 running (`--max-jobs 6`, cap 40). **Critical path** (first submission 16:10:46Z -> last fetch 00:00:13Z) **28,167 s**
(7 h 49 min) for 25,932 s of node time: queue-dominated - 119,824 job-seconds of queue wait (per job 3-133 min, median ~47) from the
dispatcher's r8g.48xlarge capacity holds (`CF:ROLLBACK_COMPLETE`, 650 held polls, up to retry 15; the 5-edge's two first worker jobs
waited 2 h 13 min, its reducer 38 min). The planner's estimated critical path per coupon (77-96 min at 2.0x) assumed immediate starts;
the measured per-coupon paths were 1.9-3.4 h. Transport: the login host unreachable 16:17-17:00Z and 22:27-22:43Z (ssh rc 255; 126
transport-failure polls, every job kept active, nothing recomputed) and one 14-min ssh stall 23:20-23:34Z; the driver ran unattended
16:09:39Z -> 00:00:13Z (WallClockSeconds 28,259) and exited on its own. Two reporting workers were lost to outages; none of the driver.

## BEFORE / AFTER (2026-09-21 -> 2026-09-22)

| quantity | 2026-09-21 (decision 58) | 2026-09-22 (decision 64a) |
|---|---|---|
| registration of the six coupons | 990 s (census probes, serial) | 33 s (labels-only probes, parallel) |
| build pool wall / end to end | 7,338 s / 8,493 s (pool of 2) | 1,462 s / 1,527 s (pool of 3) |
| serial sum of the case walls | 13,175 s | 3,591 s |
| qualify node-h | 8.257 (5 coupons; 7f03 fail-closed) | 7.203 (6 coupons; the same five 5.21) |
| reducer seconds / share of job seconds | 11,964 s / 40.3% | 626 s / 2.4% |
| jobs / concurrency | 5 / 2 | 35 / 6 (split) |
| critical path | 24,313 s (17,155 s without the 7,158 s outage dead time) | 28,167 s (queue-dominated: 119,824 job-s of capacity holds) |
| cost model | physics-11 / b28 / b = 6 (`22ee226d…`) | refit from this run (`qualify/cost-model.json`; see below) |

**Cost-model refit** (`qualify/refit_cost_model.py`, record `qualify/device-library-20260922/cost-model-refit.json`): every rate the
largest the six coupons imply on the 7f03 reference mesh (p4 1.455 s per PCG iteration, 31.3 s per source, worker non-source 78 s per
job, reducer 260 s vs 157.5 measured / 1,050 under the previous model; p3 / p5 / local-edge likewise; evaluation fraction and resident-
field rate KEPT); the previous model kept byte for byte as `cost-model-physics11-b28.json`. Self-check per job: minimum stage estimate /
measured wall 1.002 (the previous model's 0.835 on the two-edge p5 control worker), 1x job estimate / stage wall 1.02-1.24, worst-PCG
with margin over the coupon's node seconds 3.18-3.76 (previous 3.02-3.62): conservative. 7f03 still does not fit one job (22,681 s at
2.0x + margin) and the speed split reproduces the recorded N = 6.

## What remains manual

The seed copy with the config's interface layers; reading PendingQualification (no reference for device geometries: user decision 59
accepts p-sequence-only entries, `LibraryQualified` stays false in the tool); the MA definition for consumers (MA_raw vs MA_sharp;
Palace reads MA_raw); thin coupons, corners (30) and isolated edges (3,082) - decision 65 / 66 on the device-verification branch;
the r8g capacity holds (a cluster condition, not the tool's: the planner's critical path is a lower bound); the recorded queue waits.

## Validation of the recording commits

`8ebe7340a` (records), `05f664772` (refit + tests): test_refit_cost_model 7 (RecordedRefitTest un-skipped), test_qualify_gates 22,
test_qualify_dry_run 26 OK; full `unittest discover` and the three preflights (`preflight-final/`: suite 15 PreflightPassed with the
rounded-strip UnsupportedClass TopRounding as recorded, calibration-sizing 4 / 4, calibration-ma 6 / 6), `refreeze_manifest_tools.py
--check` current - results in the final commit's message.
