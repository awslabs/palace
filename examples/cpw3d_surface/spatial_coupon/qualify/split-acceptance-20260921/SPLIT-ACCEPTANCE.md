# Acceptance of the per-coupon source split (supervisor decision 61b, user decision 60(2); 2026-09-21)

`coupon-library qualify --reference none --controls p3,p5 --max-jobs 2 --job-policy fixed --fixed-jobs 2`
on the two-edge 10 coupon (`two-edge-8dd4bc70f183`, mesh `5d01204e3396…`, 78 sources, the decision-51
build `/tmp/decision-61-split/build-two-edge-tol8`), root `/tmp/library-split-acceptance-01` (remote
`coupon_accuracy_assessment_20260913/library-split-acceptance-01`), frozen binary `b28f089a…`, tool
commits 7045a37d6 (plans, submission) / 9cf2b837f + 415ed0428 (the two fix-forward commits below; the
resumed driver). Compared with the single-job acceptance record `/tmp/library-acceptance-01` (PBS 46023,
2026-09-20, the same mesh and p4 configs; its p4 reducer matrices).

## Run

| job | PBS | stages (s) | wall s | estimate s (2x PCG, +35% + 300 s) | node-h |
|---|---|---|---|---|---|
| worker-1 (block 1 = sources 1-13; controls + local-edge) | 46337 | p4-worker-block1 66, p5-control-worker 138, p5-control-reducer 23, p3-control-worker 25, p3-control-reducer 12, p4-local-edge 74 | 339 | 1206 | 0.094 |
| worker-2 (block 2 = sources 14-78) | 46338 | p4-worker-block2 271 | 272 | 1205 | 0.076 |
| reducer (union 14,976 = 78 x 192 potentials, `archive-union.json`) | 46339 | p4-reducer 101 | 102 | 512 | 0.028 |

Total 713 s = **0.198 node-h** over 3 jobs; `CriticalPathSeconds` 2,346 (first `qsub` 21:04:24Z -> fetch
21:43:30Z) of which 1,689 s is the reducer job's queue time under three SOCA capacity holds (`Hold_Types
u`, `Resource_List.error_message CF:ROLLBACK_COMPLETE:retry=1..3`, released by the dispatcher at
`retry_eligible_after`; ran 21:41:38Z); the compute path is worker-1 339 s + reducer 102 s = 441 s. The
single job (46023) ran the same p4 stages in 325 + 99 + 197 + 75 = 696 s (its 2,071 s / 0.575 node-h
include the p5 main stage 1,375 s that this run did not request): the split's overhead is +17 s (+2.4 %,
three job start-ups), the p4 main stage 338 vs 325 s worker and 101 vs 99 s reducer. Verdict
PendingQualification (no reference; the p-sequence controls passed - `comparison/p-sequence-controls.md`);
the fetched CSV digests matched the remote (`result-csv-sha256.json`), every reducer matrix complete /
symmetric / nonnegative (`matrix-validation.json`), the three archives (4.8G p4, 961M p5-control, 214M
p3-control) deleted after the digest check (`remote-archive-deletion.json`, `Remaining` empty).

## Matrices: split reduction = single job to roundoff (PASS)

`qualify/compare_split_matrices.py` on the two p4 reducer outputs (`split-vs-single-job-matrices.json`;
rows matched by key columns, 12 printed significant digits):

| matrix / column | entries | bit-for-bit | max per-entry relative difference (source) | max absolute | relative to the largest entry |
|---|---|---|---|---|---|
| domain Q_ij (J) | 3,081 | 3,042 | 8.649e-13 (13; then 53: 8.43e-13, 59: 6.76e-13; 26 / 78 sources differ at all) | 1.0e-29 | 2.27e-14 |
| surface Q_ij, Q_ij normal, Q_total_ij, Q_total_ij normal (J) | 6,162 each | 6,111 | 8.470e-13 (53; then 58: 7.88e-13, 17: 7.35e-13; 34 / 78 sources) | 1.0e-31 | 3.42e-14 |
| surface R (m), Q_ij tangential, Q_total_ij tangential | 6,162 each | 6,162 | 0 | 0 | 0 |

`EQUAL True` at the 1e-9 tolerance: every difference is one unit in the 12th printed digit (the union's
file order changes only the summation order of the linear reduction). Per-source maxima for every source:
`PerSourceMaxRelativeDifference` in the record.

## Defects found by this acceptance, fixed forward

1. **A held PBS job was read as "left the queue"** (commit 9cf2b837f). The reducer job was held by the
   SOCA dispatcher 4 min after its submission (capacity: the compute-node stack rolled back) and the
   driver - `poll_job` / `remote.monitor` tested `not in ("Q", "R", "E")` - stopped the coupon with
   `Fetch: no status.json fetched for reducer` (`library-qualification.first-driver.json`). Fix:
   `remote.IN_QUEUE_STATES` = Q R E H W T S B (`in_queue`), the poll's `qstat` grep carries `Hold_Types`
   and `error_message` and the held reason is logged; `test_held_job_is_still_in_the_queue`.
2. **`--resume` of a split coupon compared byte-identical plans in two orders** (commit 415ed0428): the
   recorded plans by sorted glob (reducer, worker-1, worker-2) against the derived plans in job order
   (worker-1, worker-2, reducer) - `Resume: the plans derived now differ` on identical files
   (`library-qualification.second-driver.json`). Fix: `plans_text` concatenates both sides in path
   order; the split fake-remote test resumes the coupon (three submissions adopted, no upload / qsub).
   The third driver (`--resume`, identical command) adopted 46337 / 46338 / 46339 and finished the run.

Records here: the library records (three drivers), gates, `process-library.json`, per coupon the three
plans / submissions / `status.json` / `qstat -xf`, cost, qualification, archive union, matrix
validation, digests, archive deletion, the p-sequence controls and the matrix comparison. No CSV, mesh,
log or archive.
