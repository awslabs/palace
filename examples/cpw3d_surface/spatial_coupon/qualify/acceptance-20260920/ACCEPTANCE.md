# coupon-library qualify — live acceptance run (supervisor decision 51, 2026-09-20)

<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

Repo `palace2`, build at `96ec539c4`; qualify tool fixed forward during the run to `526c282da`
(commits `9a159cc5d` concurrency / reference-order main stage / NotApplicable interfaces,
`bf7053855` tool commit in the record, `c6d9fe439` portable rsync upload, `526c282da` `--resume`
+ sliced wait). Local root `/tmp/library-acceptance-01`; remote
`soca-green-job:/data/home/simlap/coupon_accuracy_assessment_20260913/library-acceptance-01/`
(the command's remote ROOT holds the frozen binary `b28f089a…` and `mpiexec_bound.sh`; the run
directory is the local root's name). Reference mirror `/tmp/library-acceptance-01-reference`
(195 files fetched read-only from the remote graded_v2 `inputs/07`, `inputs/10`,
`cases/{07,10}-fabricated/{worker.json,status.json,reducer/*.csv}`; every SHA-256 equals the
remote's and the recorded campaigns' `reference/` copies; `local-mirror-sha256.txt`).

## Step 1 — `coupon_library.py build` (root `/tmp/coupon-library-build-96ec539c4-20260920-125056`)

| case | status | identity SHA-256 | elements (tet/prism/pyr) | H1 p4 | estimate/actual | wall |
|---|---|---|---|---|---|---|
| four-edge-9d2cb9bbb3fe | built, passed | `1d536c44…9ba5` = physics-13 mesh (byte-identical) | 1,916,486 (1,731,140/172,107/13,239) | 25,392,150 | 1.026 | 1913 s (pool of 2) |
| two-edge-8dd4bc70f183 | built, passed | `5d01204e…9b04` = the production root since decision 44; NOT gallery-10's `15e93865…` (a22b471c1, 516,662 el.) | 521,676 (425,916/88,920/6,840) | 7,973,827 | 1.074 | |

`build` has no cross-root cache: both cases were rebuilt from source under the manifest bounds
(every stage ≤ 121 s / 4.6 GiB of 1800 s / 8 GiB); reuse is by content — the four-edge identity
mesh reproduces the recorded physics-13 mesh byte for byte, the two-edge mesh reproduces the
decision-44 production root (the gallery-10 campaign predates decisions 43/44, so its mesh is
not reproducible from HEAD). CanonicalBuildIds are new (`54e800bc…`, `a9499778…`: the frozen
tool set changed since the recorded roots).

## Step 2 — dry run (`/tmp/library-acceptance-01-dryrun`, `dry-run-diff-vs-recorded.txt`)

Main-stage `worker.json` / `reducer.json` equal `four-edge-physics-13/main/pr-p4/*` and
`gallery-physics-10/main/g10-p4/*`, `g10-p5/*` apart from paths; plan stage lists (names,
config, environment, dependencies, order) equal the recorded plans (prefix aside; two-edge gets
the p5 main stage because its reference is p5); every trace pin equal (80 / 78); the four-edge
mesh pin equal, the two-edge mesh pin the production mesh. Differences: the control subsets
(class-based defaults `[1,2,4,5,10,17,20,49]` / `[1,2,3,5,6,20,25,49]` vs the recorded
supervisor-specified `[1,7,23,26,31,34,35,47,48,80]` / `[1,7,21,25,26,35,43,78]`), and the
estimate-derived caps (e.g. four-edge p4-worker 5400 s vs the hand-set 9000 s; measured 1380 s).

## Step 3 — live run: PBS 46023 (two-edge) and 46024 (four-edge), concurrent

Submitted 20:35:31 / 20:35:41 UTC (1 unrelated user job running, cap 40), both started 20:38:26
UTC on one r8g.48xlarge each (192 ranks), `job_state F`, `Exit_status 0`, walltime 00:34:32 /
00:48:50; in-job preflight passed (pins, binary `b28f089a…`); 0 PCG non-convergences; 16 / 14
result CSVs hash-verified against the remote; every reducer matrix validated; the remote
`archive/` directories deleted after verification (two-edge 4.8G + 9.2G + 961M + 214M; four-edge
16G + 3.0G + 665M; `Remaining ''`, remote `main/` now 5.3M / 9.8M).

### Reproduction table

| coupon | gated comparison | verdict | E within 1% | p_MA within 1/2/5% (median, weighted, strongest-20 within 2%) | p_MS within 1/2/5% (median) | p_SA within 1/2/5% | p-seq controls |
|---|---|---|---|---|---|---|---|
| four-edge, RECORDED physics-13 | p4 vs p4 ref | Passed | 60/60 (worst −0.99% at 10) | 30/53/60 (+0.79%, +0.77%, 20/20) | 47/59/60 (−0.63%) | 34/46/58 | pass |
| **four-edge, LIVE** | p4 vs p4 ref | **Passed** | **60/60 (−0.99% at 10)** | **30/53/60 (+0.7915%, +0.7675%, 20/20)** | **47/59/60 (−0.6275%)** | **34/46/58** | pass (8 class controls) |
| two-edge, RECORDED gallery-10 (through the command's gates) | p5 vs p5 ref | Failed (p_MA strongest-20: 53, 58) | 78/78 (+0.92% at 18) | 33/59/78 (+1.28%, +1.14%, 18/20) | 75/78/78 (+0.08%) | n/a (no SA interface) | pass |
| **two-edge, LIVE (production mesh)** | p5 vs p5 ref | **Failed (p_MA strongest-20: 53, 58)** | **78/78 (+0.68% at 7)** | **34/60/78 (+1.25%, +1.12%, 18/20)** | **77/78/78 (+0.07%)** | n/a | pass |
| two-edge, informational p4 vs p5 ref: recorded → live | | | 70/78 → 70/78 | 57/76/78 → 56/76/78 | 66/78/78 → 66/78/78 | | |

### Live CSVs vs recorded CSVs (per-source maximum relative difference of every observable)

- **four-edge (identical mesh `1d536c44…`)**: E, Q_MA, Q_MS, Q_SA, p_MA, p_MS, p_SA — **0.000e+00
  at every one of the 80 sources**; the whole domain matrix (3,240 entries) and surface matrix
  (58,320 entries) are **bit-identical** to `four-edge-physics-13/results/main/pr-p4/reducer/*`
  (the run is deterministic on the same node type and rank count; 12-digit CSVs).
- **two-edge (production mesh `5d01204e…` vs gallery-10's `15e93865…`; the decision-43/44 recipe
  step, never measured on two-edge before)**: p5 stage — E max 5.1e-3 (18), Q_MA 1.9e-2 (66),
  Q_MS 2.1e-2 (35), p_MA 1.8e-2 (66), p_MS 2.1e-2 (35), signed medians |≤ 4.6e-4|; p4 stage —
  E 9.1e-3 (18), p_MA 1.2e-2 (17), p_MS 1.7e-2 (35). Off-diagonal maxima (0.24 / 0.099 at
  (9,29)) are on entries ~1e-4 of the diagonal scale. Explained: 521,676 vs 516,662 elements
  (H1 p5 15,510,556 vs 15,397,771) from the altitude-measured trace rule and the composed curve
  spacing; the class counts move by ≤ 2 sources and the gate outcome is identical.

### Cost accounting (`library-qualification.json` → `Library`, per-case `Cost`)

| | two-edge p4 | two-edge p5 | four-edge p4 | recorded |
|---|---|---|---|---|
| node-h per coupon | 0.1178 (424.0 s) | 0.3818 (1374.4 s) | 0.5318 (1914.6 s) | 0.118 / 0.383 / 0.532 |
| PCG mean / max | 11.90 / 29 | 10.82 / 33 | 20.05 / 40 | 11.9/30, 11.0/34, 20.1/40 |
| vs reference node-h | 0.032× (3.680) | | 0.146× (3.639) | 0.032× / 0.146× |

Jobs 2 submitted / `--max-jobs` 2 / cap 40; job wall 2071 s + 2929 s = **1.389 node-h**;
**critical path 3141 s** (first submission 20:35:31 → last fetch 21:27:52 UTC; the jobs
overlapped); controls: two-edge p5 137 + 23 s, p3 27 + 10 s, local-edge 75 s; four-edge p5
531 + 84 s, p3 87 + 42 s, local-edge 269 s. Process-library entries written
(`process-library.json`: four-edge `LibraryQualified true`, two-edge `false`, verdicts bound).

## Defects of the command found and fixed forward (each its own commit; never by editing the recorded campaigns)

1. p_SA gate and p-sequence p_SA observable failed on a reference without an SA interface
   (gallery 10) → `NotApplicable`, recorded, never a failure (`9a159cc5d`).
2. Coupons ran sequentially in one invocation → up to `--max-jobs` concurrent jobs, round-robin
   read-only monitoring, fetch / analyze on completion, measured critical path (`9a159cc5d`).
3. `--orders` was global → the reference order is added as a main stage per coupon and the
   same-order comparison is gated (others informational) (`9a159cc5d`).
4. `rsync --mkpath` is rejected by the macOS openrsync client; the runner was uploaded as a
   directory → ssh `mkdir -p` + plain rsync, `upload_file` (`c6d9fe439`; first live upload).
5. The first driver's single `time.sleep(90)` never returned (31 min observed; a probe in the
   same conditions slept exactly 90 s) and the "fetch later with the recorded job id" stop had
   no mechanism → `--resume` (adopts `<case>/submission.json` when the re-derived plan is
   byte-identical; used for this run at 21:07 UTC, no second submission) and a sliced wait
   (`526c282da`). The first driver's checkpoint is kept as
   `library-qualification.first-driver-checkpoint.json`.

## Verdict

ACCEPTANCE = **PASS**: the four-edge coupon reproduces the recorded physics-13 qualification
exactly (bit-identical CSVs, verdict Passed, every class count and the cost); the two-edge coupon
reproduces the recorded gallery-10 gate outcome (Failed on the p_MA strongest-20 at 53 / 58 —
a genuine two-edge finding at equal p, not a command defect) with class counts within ≤ 2
sources on the production mesh; cost accounting and process-library entries filled.
