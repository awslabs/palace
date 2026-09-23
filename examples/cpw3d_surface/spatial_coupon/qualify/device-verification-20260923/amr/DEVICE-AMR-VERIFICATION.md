# AMR device-level verification of the transmon process library (decision 68, 2026-09-23)

Branch `simlapointe/device-verification` (worktree /tmp/palace2-devver-20260922). Evidence: /tmp/library-device-amr-01 (cluster/:
configs, job scripts, fetched CSV/JSON/logs per run, tabulate_amr.py, amr-tabulation.txt, amr-summary.json; smoke-p2/,
feasibility-edgeref/) and /tmp/device-prism-era-20260923 (prism-era copies + prism-era-tabulation.txt). Committed here: this record,
device-amr-summary.json (per cycle: DOFs, elements, raw / corrected values, closures, confidence, min h, final-mesh near-edge sizes,
wall, memory) and the p5 job script. Follows the coarse-mesh record ../DEVICE-VERIFICATION.md (decision 66).

## Setup
- `prepare_surface_response_electrostatic.py --amr-max-its 10 --trace-coupling SurfaceMortar --mortar-oversampling 2
  --save-adapt-mesh --amr-max-size <guard>` (commit 2661ba44a: the Refinement block of the recorded prism-library device runs, Tol
  0.01, UpdateFraction 0.7, MaximumImbalance 1.15, MaxNCLevels 8, Nonconformal, SaveAdaptIterations; default still MaxIts 0; the
  checked-in transmon_surface_amr.json EdgeRefinement tubes are opt-in only: measured x3.5 elements per pass, > 1e9 elements if
  applied) + the two recorded decision-66 edits (UnmatchedPolicy Warn, MS layer 11.47). Substrate 11.45, CorrectionMode Both,
  FixedTrace domain correction, SolveTol 1e-6. MaxSize guards 20M / 40M / 80M DOFs (p3 / p4 / p5) never triggered.
- Mesh device-run/mesh/transmon_surface_p1_island9.msh2 (sha256 c4fe532f...), library /tmp/library-device-full-01
  (process-library.json sha256 42b19686..., 11 models), frozen executable 170439c4..., m8g.48xlarge, remote root
  soca:/data/home/simlap/coupon_accuracy_assessment_20260913/device-verification-01/amr/.
- Jobs: 48246 p3 (1 node, 192 ranks), 48248 p5 (4 nodes, 768 ranks; USER: p5 on 4 m8g, MaxIts 10), 48257 p4 (4 nodes, 768),
  48256 p4 (2 nodes, 192 ranks = 96/node, --map-by ppr:96:node), 48405 p3 (4 nodes, 768). Multi-node because the final solves need
  the memory (p5 99M DOFs: 1.3 TB total; p4 64M: 0.76-0.95 TB) and because 768 ranks give the complete edge matching (below). Failed attempts: 48247 p4 at 384 ranks aborted in 35 s
  (ExtractMetalEdgeGeometry "Unable to associate a metal perimeter edge with one supporting metal surface", metaledge.cpp:557:
  partition-dependent; 192 and 768 ranks pass); 48250 p4 "192 ranks on 2 nodes" was packed on one node by Open MPI (hostfile slots
  ignored; --map-by ppr:96:node fixes it) and its 63.7M-DOF final solve was OOM-killed (cycles 1-10 kept as p4-attempt2-1node).
- Per-cycle near-edge diagnostic: min h per solved mesh from the Palace log (x Lc = 4 mm); the adapted mesh through
  examples/transmon/amr_edge_resolution.py (commit 3490b3449: leaf tetrahedra within 2 / 0.5 / 0.1 um of the metal edges, min / P05 /
  median longest edge). Palace quirk found: SaveAdaptMesh writes through the iterationNN symlink left by SaveAdaptIterations, so only
  the final mesh survives a run (the iteration02 copy is the same file) -> the per-cycle sizes come from the log's min h.

## AMR path (identical at every order to 1 %; p5 reproduces the prism-era element counts to 0.1 %)
Elements 138,247 -> 138.9k -> 140.9k -> 146.2k -> 159.0k -> 191.2k -> 271.4k -> 468.4k -> 998.1k -> 2.40M -> 6.29M (p5; p3 7.34M, p4
8.50M): the refinement is Dorfler-marked at the metal edges (min h 0.952 um for cycles 1-3, then 0.602, 0.301, 0.286, 0.143, 0.072,
0.036, 0.018, 0.0090 um). H1 DOFs: p3 0.689M -> 19.6M, p4 1.60M -> 64.6M, p5 3.07M -> 99.4M. Final mesh within 0.1 um of the edges
(p5): 502k elements, longest edge min 15 nm / P05 45 nm / median 110 nm (shortest element edge anywhere 5.5 nm): the device mesh never
refines below the 2 nm thin cutoff of the library; within 2 um: 2.40M elements, median 0.25 um.

## Results per cycle vs paladin_runs single-transmon-es-r5nm-t5um-p5 (SA 9.9135e-5, MS 2.1274e-4, MA 4.0002e-6, C 8.70481e-14 F)
sc = self-consistent corrected; MA also x 1.0852 (MA_sharp / MA_raw of the 2D isolated-edge p3 model = 99.4 % of the length);
conf = postprocessing confidence flag (0 = fail), spread = max trace-closure spread. Full tables: device-amr-summary.json.

p5, 768 ranks (PBS 48248), matched 3,071 segments / 34 corner / 6 spatial patches:
| cycle | DOFs | SA raw / sc | MS raw / sc | MA raw / sc / sharp | C raw / corr | conf, spread |
|---|---|---|---|---|---|---|
| 1 (coarse) | 3.07M | -46.9 / -24.9 % | -52.0 / -27.1 % | -77.5 / -27.9 / -21.8 % | +3.3 / +3.0 % | 0, 0.075 |
| 4 | 3.19M | -40.1 / -14.9 | -43.6 / -14.5 | -73.2 / -15.6 / -8.4 | +1.9 / +1.7 | 0, 0.081 |
| 6 | 3.86M | -30.4 / -7.9 | -31.0 / -6.8 | -68.3 / -6.2 / +1.8 | +1.0 / +0.8 | 0, 0.084 |
| 8 | 8.06M | -17.5 / -4.1 | -13.7 / -2.4 | -60.5 / -1.1 / +7.3 | +0.6 / +0.3 | 0, 0.085 |
| 9 | 16.1M | -11.0 / -3.4 | -4.8 / -1.4 | -56.4 / -0.3 / +8.2 | +0.5 / +0.2 | 0, 0.085 |
| 10 | 38.0M | -3.9 / -3.0 | +4.9 / -0.8 | -52.0 / +0.2 / +8.8 | +0.5 / +0.2 | 0, 0.085 |
| 11 (final) | 99.4M | +3.1 / **-2.7 %** | +14.4 / **-0.5 %** | -47.6 / **+0.5 / +9.0 %** | +0.4 / **+0.1 %** | 0, 0.086 |
Final absolute: SA 9.650e-5 (fixed-trace 9.418e-5, fixed-flux 1.005e-4), MS 2.117e-4, MA 4.019e-6 (sharp-scaled 4.361e-6), C 8.716e-14 F.
Against r5nm p4 / r10nm p5: SA -2.3 / -2.3, MS -0.1 / -0.4, MA +2.7 / +0.7 (sharp +11.4 / +9.3), C +0.1 / +0.1 %.

p4, 768 ranks (48257), same matching: cycle 1 -29.9 / -30.2 / -35.6 / +4.6 %; cycle 6 -9.2 / -8.7 / -9.2 / +1.0; cycle 10 -3.2 / -1.2 /
-0.1 / +0.2; final (64.6M DOFs) **SA -3.0 / MS -0.8 / MA +0.3 (sharp +8.8) / C +0.1 %** (raw -0.4 / +10.6 / -49.1 / +0.5).
p3, 768 ranks (48405), same matching: cycle 1 -44.4 / -44.3 / -49.8 / +7.7 %; cycle 6 -17.6 / -16.4 / -20.2 / +2.2; cycle 10 -4.4 / -2.6 /
-0.9 / +0.3; final (19.9M DOFs) **SA -3.7 / MS -1.7 / MA -0.4 (sharp +8.1) / C +0.2 %** (raw -10.1 / -0.7 / -53.7 / +0.5).
Order series at 768 ranks, final cycle: SA -3.7 / -3.0 / -2.7, MS -1.7 / -0.8 / -0.5, MA -0.4 / +0.3 / +0.5 % (p3 / p4 / p5): monotonic,
p4 -> p5 <= 0.3 points.
p3, 192 ranks (48246; 3,019 segments / 4 spatial patches, 69 omitted): cycle 1 -40.5 / -39.4 / -43.3 / +7.6 %; cycle 6 -13.2 / -10.8 /
-13.1 / +2.2; cycle 10 +0.1 / +3.4 / +6.4 / +0.3; final (19.6M DOFs) **SA +0.8 / MS +4.2 / MA +6.9 (sharp +16.1) / C +0.2 %**.
p4, 192 ranks on 2 nodes (48256; 48250 identical through cycle 10): cycle 10 (23.6M) +1.3 / +4.7 / +7.3 / +0.2; final (63.7M DOFs)
**+1.6 / +5.1 / +7.6 (sharp +16.8) / +0.1 %** - vs the 768-rank p4 final -3.0 / -0.8 / +0.3 / +0.1 on the same raw values (-0.4 / +10.6 / -49.1).
Raw values are rank-independent (p4 final: -0.4 / +10.6 / -49.1 % at both counts); the 4.6 / 5.9 / 7.3-point difference in the
corrected SA / MS / MA is entirely the partition-dependent classification (192 ranks drop 52 segments = 42 "missing library model" + 26
nonparallel, and 2 of the 6 spatial patches; 768 ranks match 3,071 / 6, the 6-rank local preflight 3,068 / 6). On the coarse mesh
this was <= 0.3 points (decision 66); at the AMR operating point it is the largest single effect after the mesh itself.

## Context and side-by-side
- Coarse mesh, single solve, Collocated (decision 66): p2 -48.3 / -60.4 / -61.2 / +16.9 %; p5 -19.1 / -19.8 / -15.4 (sharp -8.2) / +3.0;
  p8 -5.5 / -2.7 / +3.9 (sharp +12.8) / +1.3 (192 ranks). The AMR p5 cycle 1 (-24.9 / -27.1 / -27.9) differs from the coarse p5 by the
  SurfaceMortar coupling and the 768-rank matching.
- Prism era (transmon_final_library_device_20260830, 2026-09-06..08; library final-library-exhaustive-p5 = prism coupons, the same
  coarse mesh without attribute 9, substrate 11.47, p5, the same Refinement block, SurfaceMortar oversampling 2, 4 r8g / 2 m8g):
  best agreement at cycle 5 (SA -3.6 / MS +0.8 / MA +1.2 / C +1.2 %), then a monotonic drift to +8.4 / +14.3 / +13.4 / +0.2 % at cycle
  12 (284M DOFs; amr10 on 2 m8g: +4.5 / +7.0 / +30.0 at 98M DOFs); the confidence flag failed at every cycle (spread 0.09-0.10, the
  Collocated variants 0.5-0.65). Same mesh path, same settings, full library: no drift (cycles 9 -> 11 move SA +0.7, MS +0.9, MA +0.8
  points and flatten); the prism-era drift was a library effect, not the device mesh refining below the thin cutoff (it does not: 15 nm
  min longest edge at the edges vs 2 nm).

## Confidence flag and coverage
The flag fails at every cycle and order (max trace-closure spread 0.073 -> 0.086, growing slightly with refinement); the three
closures (fixed-trace / fixed-flux / self-consistent) bracket the reference at the final cycle: p5 SA -5.0 / +1.4 / -2.7 %. Coverage
as in decision 66: Exact 3,107 requirements / 99.89 % of the metal-edge length; 768 ranks match 3,071 segments (12 within 2R of
another interface + 3 unsupported omitted), 34 corner and 6 spatial patches, 1 unmatched interface group (Warn).

## Cost (Palace wall; m8g.48xlarge)
p5: 586 s on 4 nodes = 0.65 node-h (per solve 23 -> 285 -> 586 s; peak 332 GB / node), p4: 451 s on 4 nodes = 0.50 node-h (peak 244
GB / node; 1,095 s = 0.61 node-h at 192 ranks on 2 nodes, peak 385 GB), p3: 477 s on 1 node = 0.13 node-h (peak 407 GB; 296 s on 4
nodes = 0.33); the near-edge diagnostic ~70 s per mesh. Reference r5nm-t5um p5: 38.1M
elements, 802M DOFs, 9,079 s on 3,840 ranks (20 nodes) = 50.4 node-h, peak 4.8 TB -> the corrected AMR device run at 99M DOFs costs
1.3 % of the reference; the seven AMR jobs (two failed attempts included) ~2.7 node-h. The prism-era p5 amr10 on 2 m8g took 1,124 s.

## Where the agreement stands and what limits it now
1. At the intended operating point (p4-p5, 10 AMR cycles, SurfaceMortar 2, full matching) the corrected device observables agree with
   the brute-force r5nm-t5um p5 reference to **SA -2.7 to -3.0 %, MS -0.5 to -0.8 %, MA +0.3 to +0.5 % (raw-cutoff) / +8.8 to +9.0 %
   (MA_sharp-scaled), C +0.1 %**, order-converged (p4 = p5 within 0.3 points) and mesh-flattening (cycle 10 -> 11: <= 0.3 points at
   p5); the reference spread (r5nm p4, r10nm p5) is 0.4 points on SA / MS and 2 points on MA.
2. Limits, in order: (a) the partition-dependent edge classification of the frozen Palace executable (192 ranks: 4.5-7 points; 384
   ranks: abort) - a Palace fix, not a library item; the reported numbers use the complete 768-rank matching; (b) the MA convention
   (raw-cutoff +0.5 % vs sharp-extrapolated +9 %; the reference is itself a ~5 nm-cutoff, ~2 %-unconverged value); (c) the confidence
   flag / closure spread 8.6 % (fixed-trace -5.0 %, fixed-flux +1.4 % on SA at p5): the closure-independent local response is not
   established by the flag even where the corrected value lands within 3 %; (d) the residual SA deficit (-2.7 %) exceeds MS (-0.5 %)
   and is stable from cycle 8 on: a library-side term (2D isolated-edge SA response, unqualified 2D / corner models, 0.11 % uncorrected
   length, the 11.47 MS layer -0.17 %) rather than device discretization.
3. Not done: a fix of the Palace classifier, the 8-edge / 2-edge spatial coupons, the r2p5nm reference for the MA convention.
