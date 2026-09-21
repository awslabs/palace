# First complete device library run (supervisor decision 58, 2026-09-21)

What a device layout -> coupon library run took end to end: `coupon_library.py build --device` on
`examples/transmon/transmon_surface_coarse.json` (the benchmark process seed with the config's
interface layers bound into a copy - the checked-in seed's MS permittivity 11.45 differs from the
config's 11.47 and Palace refuses the mismatch; recorded in the 2026-09-20 device run and repeated
byte-identically here) with the decision-57 Delaunay device basis, then `qualify --reference none
--orders p4 --controls p3,p5 --max-jobs 2` on every verified root. Repo commits `4538c7fc7` (decision
57 default; the build record's commit), `5a25d034e` (monitor transport fix), `8759f6b59` (the qualify
record's ToolCommit). Build root `/tmp/coupon-device-transmon-delaunay-20260921` (fresh production-
manifest copy, SHA `4caf40dd…`), qualify root `/tmp/library-device-transmon-01` (remote
`coupon_accuracy_assessment_20260913/library-device-transmon-01`; the build record SHA `1c4f5eaf…`).
Every MA value of this run is the sharp-edge value at the recipe's 0.25 nm cutoff (decisions 55 / 56:
~1.9-2.7% below its sharp-edge limit at the median source, 4.4-4.8% at 53/58 of the two-edge 10);
the MA definition is the user's pending decision and applies to both sides of any later comparison.

## Build (local machine, 6 cores, pool of 2; no cluster node-hours)

Discovery (3,120 requirements, 35,297 um of edge, all Missing -> 9 planned coupons, 6 spatial) +
basis generation + registration of all six spatial coupons (census-only probes, the ten-edge probe
177 s / 5.94 GiB under the 8 GiB bound): 02:01:40 -> 02:18:10 (990 s). Build pool 7,338 s. End to
end 8,493 s (2 h 22 min) of one workstation. Every device case id changed with the Delaunay default
(new cases; the 5-edge `9bfa8265e2e7` and 7f03 `5d3b5e644745` are the decision-54 ids and were rebuilt
BYTE-IDENTICAL to their recorded identities `ddf1081b…` / `e038c5ef…` - `build` does not reuse roots by
content). Cap-triangulation needles left (min altitude): 2-edge 12 (49 nm), ten-edge 20 (26 nm) -
the ring spacing, not needles.

| coupon (sources; model) | elements = tet + prism + pyr | H1 p4 | estimate / actual (x cap) | gmsh-build s / GiB | verification s / GiB | wall s | flags | identity |
|---|---|---|---|---|---|---|---|---|
| spatial-2-edge-8ce3fb213a2d (150; 3f8992613e95) | 1,215,960 = 1,052,916 + 151,398 + 11,646 | 17,099,910 | 1,254,183 / 1,215,960 = 1.031 (0.314) | 71 / 3.55 | 547 / 2.05 | 1,118 | - | `7839a4f6341a` |
| spatial-3-edge-9cd906b81cbe (175; 419576fdab24) | 1,680,987 = 1,495,011 + 172,692 + 13,284 | 22,886,084 | 1,763,366 / 1,680,987 = 1.049 (0.441) | 105 / 4.24 | 826 / 2.81 | 1,665 | - | `9b2bd97cbe74` |
| spatial-4-edge-e51d7380245e (120; 9d2cb9bbb3fe) | 1,895,487 = 1,710,141 + 172,107 + 13,239 | 25,167,936 | 1,968,010 / 1,895,487 = 1.038 (0.492) | 115 / 4.15 | 877 / 2.97 | 1,784 | - | `957e3129b495` |
| spatial-10-edge-65450ff47b9b (191; 6791f1c84123) | 2,643,905 = 2,353,601 + 269,568 + 20,736 | 35,861,283 | 2,837,463 / 2,643,905 = 1.073 (0.709) | 195 / 5.92 | 1,330 / 4.05 | 2,659 | - | `316d7a1898c2` |
| spatial-5-edge-9bfa8265e2e7 (185; 938a30f63a3d) | 2,696,939 = 2,424,275 + 253,188 + 19,476 | 35,941,049 | 2,905,085 / 2,696,939 = 1.077 (0.726) | 163 / 5.59 | 1,228 / 3.84 | 2,566 | - | `ddf1081b694f` |
| spatial-3-edge-5d3b5e644745 (225; 7f03270dca8e) | 3,037,912 = 2,635,846 + 373,347 + 28,719 | 42,696,378 | 3,135,224 / 3,037,912 = 1.032 (0.784) | 182 / 6.63 | 1,692 / 4.93 | 3,383 | verification 1692 s >= 0.9 x 1800 s | `e038c5efbb50` |

All six built, every variant (identity + rotate-z) verified Passed, 0 failures; every gate at its
production value (4M / 8 GiB / 1800 s, SJ 0.01, condition 1000, corner 4.0). Unbuildable: none -
the three 2026-09-20 fail-closed coupons (7f03 1.256 x cap, ten-edge 8.03 GiB probe, 5-edge tube
tool) are all gone with the Delaunay basis (54b) and the collar-union footprint (54a). One headroom
flag: 7f03's verification at 0.94 of the 1800 s bound under a concurrent build (54c). Out of scope,
recorded: ConcaveCorner 18 / ConvexCorner 12 occurrences -> CornerCoupon; IsolatedEdge 3,082 ->
StraightEdgeBuilder.

## Qualify - estimate gate (dry run, `/tmp/library-device-transmon-01-dry`; `library-qualification.dry-run.json`)

| coupon | runner minutes at measured PCG / at 2.0x (+35% and preflight) | node-h main stage (1.0x) | Palace peak GB of 1485 | fits one 6 h job |
|---|---|---|---|---|
| 2-edge | 59 / 130 | 0.785 | 135 | yes |
| 4-edge | 69 / 154 | 0.856 | 199 | yes |
| 3-edge 9cd9 | 94 / 201 | 1.303 | 181 | yes |
| 5-edge | 158 / 332 | 2.215 | 285 | yes |
| ten-edge | 164 / 342 | 2.314 | 284 | yes |
| 3-edge 5d3b (7f03) | 240 / 488 | 3.502 | 338 | **NO: fail closed, not submitted** |

7f03 (225 sources x 42,696,378 H1 at p4, 83,063,960 at p5): 14,391 s at the measured PCG counts,
21,460 s at 2.0x PCG, 29,272 s with preflight and margin > the 21,600 s walltime (`StoppedBy.Kind
Estimate`: "does NOT fit one job of 6 h at 2.0x PCG (488 min; largest Palace peak 338 GB of 1485
GiB): fail closed, not submitted"). A JOB-SIZE limit of the one-job-per-coupon plan, not a recipe or
mesh failure (the mesh built and verified): either a two-job split of the source set (`qualify` does
not support it) or a longer walltime - both user decisions. Excluded from the live run with `--case`.

## Qualify - live run (PBS, one r8g.48xlarge / 192 ranks per job, two concurrent; `library-qualification.json`)

Five coupons, largest first; the p-sequence gate (qualification-gates.json `aea5d5f9…`): for every
control source (8 per coupon) |d45| = |(p5 - p4) / p5| of the domain energy E within 1% and of every
participation p_MA / p_MS / p_SA within 5% (d34 = (p4 - p3) / p5 recorded). Max over the 8 sources:

| coupon | PBS | elements | H1 p4 | sources | node-h job / main stage | PCG mean / max (s per iteration) | max E |d45| | max participation |d45| | max E |d34| / participation |d34| | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| spatial-10-edge-65450ff47b9b | 46219 | 2,643,905 | 35,861,283 | 191 | 2.39 / 2.04 | 14.7 / 40 (1.15) | 0.19% (src 4) | p_MA 1.22% (src 126) | 0.49% / p_SA 1.82% | PendingQualification |
| spatial-5-edge-9bfa8265e2e7 | 46220 | 2,696,939 | 35,941,049 | 185 | 2.36 / 1.99 | 14.2 / 40 (1.18) | 0.09% (src 1) | p_MS 1.69% (src 38) | 1.13% / p_SA 2.86% | PendingQualification |
| spatial-3-edge-9cd906b81cbe | 46221 | 1,680,987 | 22,886,084 | 175 | 1.46 / 1.23 | 18.0 / 40 (0.74) | 0.32% (src 2) | p_MA 2.68% (src 2) | 0.61% / p_MA 2.83% | PendingQualification |
| spatial-4-edge-e51d7380245e | 46222 | 1,895,487 | 25,167,936 | 120 | 1.11 / 0.83 | 18.3 / 40 (0.86) | 0.49% (src 15) | p_SA 1.21% (src 15) | 1.54% / p_SA 1.84% | PendingQualification |
| spatial-2-edge-8ce3fb213a2d | 46255 | 1,215,960 | 17,099,910 | 150 | 0.94 / 0.73 | 18.4 / 39 (0.57) | 0.10% (src 7) | p_MS 1.63% (src 33) | 0.47% / p_SA 1.54% | PendingQualification |

Every p-sequence gate passed (E <= 0.49% of 1%, participations <= 2.68% of 5%, 0 failing observables);
verdict PendingQualification on all five ("no reference matrices: only the p-sequence controls were
evaluated (passed); never Passed") - **0 Failed, 0 Passed, nothing is qualified**: `process-library.json`
carries the five models (`spatialedgecluster_edgecount-{10,5,3,4,2}_…`) with their fabricated matrices,
`LibraryQualified false` on every entry. Every stage `complete`, every reducer matrix validated
(domain / surface rows 18,336 / 91,680 ten-edge … 11,325 / 33,975 two-edge; 36 / 108-180 controls), 14
result digests per coupon, the three response archives of every coupon deleted on the cluster after the
fetch (`remote-archive-deletion.json`, Remaining empty; 51 / 50 / 30 / 23 / 20 GB before deletion).

**Cost.** Job wall (PBS): 8,618 / 8,501 / 5,241 / 3,984 / 3,380 s = 29,724 s; **node-h 8.257** (main
stages 6.82; 192 ranks each, the controls and local-edge stages the rest). Estimate at 1.0x PCG vs
job: 9,828 / 8,618, 9,475 / 8,501, 5,646 / 5,241, 4,132 / 3,984, 3,540 / 3,380 - the measured PCG
estimate over by 4-14%, the 2.0x-with-margin gate value 2.3-2.4x the job. **Critical path** (first submission
11:26:25Z -> last fetch 18:11:38Z) **24,313 s** (6 h 45 min), of which 7,158 s is dead time of the
first outage (46219 / 46220 ended 13:53 / 13:51Z, fetched 15:52 / 15:55Z after the resume); the job +
queue chain alone 17,155 s (4 h 46 min: 10-edge 2:23:39 -> 4-edge 1:06:24 -> 2-edge 0:56:21 in one
slot, 5-edge 2:21:42 -> 3-edge 1:27:22 in the other). **Jobs:** 5 submitted, 2 concurrent (`--max-jobs
2`), user cap 40 (4 user jobs before the last qsub). Driver wall of the resumed run 8,603 s.

**Transport failures recorded during the outages.** (1) 11:38Z: the login host unreachable for more
than two hours with 46219 / 46220 running; the pre-fix driver read the empty poll as "left the queue",
its fetch failed and it crashed (no second submission; the jobs completed on the cluster). Fixed forward
(`5a25d034e`); `--resume` on the same root at 15:48Z adopted both submissions (`Monitor.Resumed true`,
1 poll each) and fetched them. (2) 16:21:40Z-16:40:20Z with 46221 / 46222 running: ssh rc 255 on 8
consecutive polls -> `Monitor.TransportFailures 8` on each job, jobs kept active, the polls spent from the
budget (50 / 38 of 240), no stop; both fetched normally. (3) Two workers of the reporting side were
lost to the same outages; the driver ran through unattended and finished at 18:11:38Z.

## What remains manual

The seed copy with the config's interface layers (a fixture mismatch, not the command's); the 7f03
job split or walltime; the MA definition (decisions 55 / 56 - every MA here is the 0.25 nm sharp-edge
value); reading the p-sequence verdicts: PendingQualification is not qualified - the accuracy statement
needs a reference (`--reference` on a gallery-anchored case; there is none for these device geometries),
so registering these five entries as LibraryQualified is a user decision the tool refuses to take; the
corner (30 occurrences) and isolated-edge (3,082) requirements outside the spatial-coupon scope.
