# Device-level verification of the transmon process library (decision 66, 2026-09-23)

Branch `simlapointe/device-verification` (worktree /tmp/palace2-devver-20260922); evidence root /tmp/library-device-full-01
(`device-run/`: configs, preflight, cluster-4789{8,9}/ and cluster-47906/ fetched outputs, local-p{2..6}/, compare.py,
comparison.txt, device-verification-summary.json). Committed: qualify/device-verification-20260923/ (this document + the
summary JSON), examples/transmon/tag_island_conductor.py (+ test), the writer fix in qualify/qualify_library.py (+ tests).

## Setup
- Device: examples/transmon/prepare_surface_response_electrostatic.py on the coarse thin transmon mesh, Electrostatic, island =
  Terminal 1 (attr 9), ground = [5, 6, 7], SA attr 8, MS / MA on [5, 6, 7, 9] at 2 nm; `--correction-mode Both`,
  TranslationalDomainCorrection FixedTrace, Collocated traces, SolveTol 1e-6, BoomerAMG / CG 1e-10, MaxIts 0 refinement.
- Two edits of the produced JSON, recorded: `UnmatchedPolicy` Error -> Warn (the script hard-codes Error); MS layer permittivity
  11.45 -> 11.47. Permittivity choice: substrate `--substrate-permittivity 11.45` (the references' and the coupons'
  `Fabrication.SubstratePermittivity`); the library `InterfaceLayers.MS` is 11.47 (the coupons' MS layer) and Palace checks the device
  layer against it to 1e-10, while the script ties MS := substrate -> the MS edit. Bias vs the 11.45 references: p_MS ~ -0.17 %.
- Mesh: the checked-in `transmon_surface_p1.msh2` has NO attribute 9 (metal = one physical surface 5); PBS 47890 failed on the
  terminal. New tool `tag_island_conductor.py` (edge-connected components of surface 5; the component touching no port patch 6 / 7
  becomes 9 "island": 139 of 13,086 triangles, a 24 x 625.5 um strip with a 30 um gap; the ground plane and the port-joined feedline
  stay 5) -> `device-run/mesh/transmon_surface_p1_island9.msh2` (sha256 c4fe532f...). Reference terminal area 14,885 um^2 agrees.
- Library: /tmp/library-device-full-01/process-library.json (Version 3, R 2.0, 11 models). Its six fabricated SURFACE matrices were
  the reducer's shelled files (MA shells 4.., no MA base index): the geometry-only preflight loads them but Palace aborts the solve
  ("interface mapping refers to a missing coupon surface response"; local + PBS 47897). Cause: `write_process_library --previous` kept
  the -03 Version-1 models as-is; the shell collapse ran only for the run's own cases. Fix (749c94899): kept models are collapsed
  through the previous root's `library-qualification.json` RadialShells, and every written model's `Interfaces[].Coupon` is checked
  against both surface matrices (fail closed). Library re-assembled (assemble.py unchanged): 376 MB (was 801), all 22 checks pass.
- Executable: the frozen `palace-archive-estimate-170439c4...bin` (source b1e7e9e9d; `git log b1e7e9e9d..HEAD -- palace/` is empty on
  the branch -> repository-consistent; ResponseCorrection present). Local dry run: `palace -dry-run` + `-serial
  --surface-response-preflight` (3.4 s) on the exact config before the cluster. Instance m8g.48xlarge (USER decision 67; the run
  peaks at 0.3-0.7 GB / rank), 192 ranks; remote root .../coupon_accuracy_assessment_20260913/device-verification-01.
- Jobs: 47898 = the deliverable (p2 + p3, 192 ranks); 47899 (p2 at 6 / 48 ranks; p4-p6 at 192) and 47906 (p6 at 48; p7, p8 at 192)
  extend it because every solve costs seconds. Local 6-rank checks (main-tree build v0.17.0-370) at p2-p6 reproduce the cluster
  numbers to 0.4-1.8 % (the difference is the rank-dependent matching below, not the build).

## Coverage (Warn; `preflight-p2/postpro/surface-response-requirements.json`, the same in every 6-rank run)
Exact 3,107 requirements / 35,282.2 um = **99.89 % of the metal-edge length**; Missing 2 / 40.0 um, both SpatialEdgeCluster: an
8-edge (32 um) and one 2-edge (8 um) cluster without a coupon. Matched by class: IsolatedEdge 3,066 segments / 35,069.8 um (99.4 % of
the length -> the 2D `isolated-edge` model carries the correction), corners 34 patches / 136 um (concave 22 / convex 12; the two
CornerCoupon models), SpatialEdgeCluster 6 patches / 76.5 um (10-edge 26.45, 2-edge 8, 3-edge 8 + 10, 4-edge 2 x 12; the 5-edge model
unused); 20 nonparallel segments omitted, 4 unmatched vertices.
**Rank dependence (same binary):** 6 ranks 3,068 segments / 6 spatial patches; 48 ranks 3,065 / 6; 192 ranks 3,019 / 4 (69 omitted:
26 nonparallel + 42 "missing library model"). Effect on the corrected values at p6: <= 0.3 points (48 vs 192 ranks). A Palace
classifier defect (partition-dependent segment classification), recorded, not fixed here.

## Results vs paladin_runs single-transmon-es-r5nm-t5um-p5 (SA 9.9135e-5, MS 2.1274e-4, MA 4.0002e-6, C 8.70481e-14 F)
Per interface p_surf: raw thin sheet | self-consistent corrected (sc) | MA also sc x 1.0852 (MA_sharp / MA_raw of the isolated-edge
p3 model); relative differences to the reference. C corrected = C_raw x E_corr / E_raw (terminal-C.csv is raw only).

| run (ranks) | SA raw / sc | MS raw / sc | MA raw / sc / sc x sharp | C raw / corr |
|---|---|---|---|---|
| 47898 p2 (192) | -63.2 / **-48.3 %** | -78.0 / **-60.4 %** | -89.6 / **-61.2 / -57.9 %** | +17.0 / **+16.9 %** |
| 47898 p3 (192) | -56.5 / -39.5 % | -65.6 / -38.1 % | -83.8 / -40.0 / -35.0 % | +7.8 / +7.6 % |
| 47899 p4 (192) | -51.2 / -24.5 % | -57.9 / -23.4 % | -80.2 / -24.3 / -17.8 % | +4.7 / +4.5 % |
| 47899 p5 (192) | -46.9 / -19.1 % | -52.0 / -19.8 % | -77.5 / -15.4 / -8.2 % | +3.3 / +3.0 % |
| 47899 p6 (192) | -43.4 / -12.6 % | -47.1 / -11.2 % | -75.3 / -6.1 / +1.9 % | +2.4 / +2.2 % |
| 47906 p7 (192) | -40.4 / -9.0 % | -43.0 / -4.0 % | -73.3 / +0.8 / +9.4 % | +1.9 / +1.6 % |
| 47906 p8 (192) | -37.8 / **-5.5 %** | -39.3 / **-2.7 %** | -71.6 / **+3.9 / +12.8 %** | +1.6 / **+1.3 %** |

Absolute values at p8: SA 6.164e-5 raw -> 9.370e-5 sc (fixed-trace 9.201e-5, fixed-flux 9.734e-5); MS 1.291e-4 -> 2.069e-4
(2.023e-4 / 2.105e-4); MA 1.136e-6 -> 4.156e-6 (4.243e-6 / 4.039e-6); C 8.843e-14 -> 8.817e-14 F. At p2: SA 3.649e-5 -> 5.130e-5,
MS 4.682e-5 -> 8.420e-5, MA 4.180e-7 -> 1.552e-6, C 1.018e-13 -> 1.017e-13 F. Fixed-trace / fixed-flux / self-consistent closures
agree within 1-6 % at every order; the postprocessing confidence flag FAILS at every order (max trace-closure spread 0.070-0.077,
93 % of the response above the 5 % limit). Against the other references (r5nm p4, r10nm p5) the p8 differences shift by <= 1 point
(SA / MS) and 2-3 points (MA, the references' own ~2 % floor).

## Cost
Every device solve is seconds: PBS 47898 p2 15.0 s wall (Palace 11.3 s, Response Setup 4.2 s), p3 12.6 s, 285-305 MB / rank;
p4 15 s, p5 21 s, p6 37 s, p7 64 s, p8 139 s (717 MB / rank) at 192 ranks; p6 at 48 ranks 87 s. The five jobs (two failed attempts
included) used ~0.3 node-h; the library upload is 376 MB. Locally (6 ranks): p2 11 s, p4 51 s, p6 405 s / 2.9 GB.

## What limits the agreement
1. **The coarse device mesh, not the library.** 138k tets; the island's triangle edges are 14.9 um median / 46 um max, the SA plane
   6.7 um median: the raw thin-sheet field trace at R = 2 um from the edges is unresolved at low order. Raw p_surf move 30-70 % and
   the corrected values 40-55 points between p2 and p8, all three interfaces at the same relative deficit at each order (a global
   under-resolution of the trace amplitude); the corrected values approach the reference monotonically (p7 -> p8: SA +3.5, MS +1.3,
   MA +3.1 points) and are not converged at p8. C itself is +17 % at p2 and +1.3 % at p8 (raw +1.6 %; the thick-metal correction
   of C is -0.3 %). The library correction is validated to the extent the raw field allows: at p8 SA -5.5 %, MS -2.7 %, MA +3.9 %.
2. MA definition: MA_raw (2 nm cutoff of the 2D isolated edge that carries 99.4 % of the length) lands at +3.9 % of the reference;
   the sharp-edge extrapolation x 1.0852 at +12.8 %. The r5nm reference resolves the edge at 5 nm and is itself unconverged by ~2 %
   (p4 -> p5), i.e. closer to a 5 nm-cutoff value than to the sharp limit; both are reported, neither is tuned.
3. Uncorrected geometry: 0.11 % of the edge length (two spatial clusters) plus, at 192 ranks, 1.6 % of the segments dropped by the
   partition-dependent classifier (<= 0.3 points at p6); the 2D isolated-edge / corner models are unqualified (part-1 gates), the thin
   coupons' PSequenceControls verdict is Failed at the fixed 2 nm cutoff (convention recorded, not converged).
4. Layer permittivity 11.47 vs 11.45 (-0.17 % on p_MS) and the confidence flag (closure spread 7 %: the raw field does not define a
   closure-independent local response on this mesh; the three closures still agree within 1-6 %).

## Next steps
- Re-run on a device mesh that resolves the trace ring (edge-adjacent elements <= 1-2 um, e.g. the transmon_surface_amr EdgeRefinement
  blocks or a graded remesh) at p2-p3: the correction should then reproduce the p8 agreement at a fraction of the dofs.
- Fix the partition-dependent edge-segment classification in Palace (3,019 vs 3,068 segments, 4 vs 6 spatial patches at 192 ranks).
- Add the 8-edge and 2-edge spatial coupons (40 um), qualify the 2D / corner models, and read the thin PSequenceControls verdict.
- Decide the MA reporting convention against a reference converged in its own edge cutoff (r2p5nm p5 did not complete).
