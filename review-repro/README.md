# Review reproducers — PR 979 (FluxLoopExcitation)

Failing cases for the review of `Boundaries.FluxLoopExcitation` and
`Solver.Magnetostatic.FluxLoopMatrixSweep`. All configs run from
`examples/circular_hole/` (they use the existing `mesh/sheet_w_two_holes.msh`
and the `double_circular_hole_multi_flux.json` setup). Results below were
observed at the PR head on linux-aarch64, 1 rank; every case was run from the
committed files.

The feature itself works as advertised: `ctrl_excitation.json` (sweep on, one
`[+1, -1]` excitation) recovers the imposed fluxoids exactly (±1.000000000e+00),
writes `terminal-fluxexc.csv` with the stored energy
E = 1.130431842459e-13 J, and exits 0.

## m1_savefield_overridden.json — `SaveField: false` is ignored under the `Save` count

- **Setup:** sweep off, `Solver.Magnetostatic.Save: 1`, one excitation with
  `SaveField: false`.
- **Expected:** no ParaView field for the excitation (docs and schema promise
  `SaveField` is "independent of the `Save` count").
- **Observed:** the excitation is written anyway. With the sweep off it lands on
  step 0, and `ShouldWriteFields(step)` also accepts `step < output_n_post`, so
  the field is written regardless of `SaveField`:
  `postpro/out_m1/paraview/magnetostatic/Cycle000002/` contains the full
  A, B, U_m field set, followed by the final save in `Cycle000003/`.

## m2_final_save_overwrite.json — the final ParaView save clobbers a swept loop's directory

- **Setup:** swept flux loops with indices 1 and 3, one excitation with
  `Index: 2` (driving loop 1), `Save: 2`.
- **Writes:** loop 1 → `Cycle000001`, excitation → `Cycle000002`,
  loop 3 → `Cycle000003`. The final save then writes to `GetCycle() + 1 = 3`,
  over `Cycle000003`.
- **Observed:** `Cycle000003/` holds only `Rank` and `Indicator` — the swept
  loop 3 field data (A, B, U_m, still present in `Cycle000001/` for loop 1) is
  gone. The user loses the field of a swept loop without any warning.
- **Note:** the final save always writes to the last written cycle + 1, so any
  source whose index equals the last written index + 1 is at risk; this also
  pre-exists on main for mixed current-port/flux-loop setups where the indices
  align. The excitation path adds an easy new way to hit it.

## m4_zero_fluxoid_nan.json — a zero fluxoid produces NaN and aborts the run

- **Setup:** sweep off, no current ports, one excitation with
  `FluxAmounts: [0.0]` — a zero right-hand side (a "ground state" row that a
  parameter sweep script can plausibly request).
- **Observed:** the CG convergence test is `res < max(rel_tol * 0, abs_tol)`,
  with the default `abs_tol = 0` this is `0 < 0` and fails, and the solver
  divides `0/0`:
  - `KSP residual norm ||r||_B = nan` from iteration 1 onward,
  - `Linear solver did not converge, norm(Ax-b)/norm(b) = nan (norm(b) = 0.000e+00)!`
  - `FluxLoopExcitation 2 combined solve did not converge; skipping!`
  - with no other states the run then aborts (exit 134):
    `Verification failed: (solve_converged_) is false: --> Initial solve did
    not converge; no converged result to fall back on!`
- **Expected:** either A = 0 with E = 0 (a legitimate zero state) or a clear
  validation error at config parse time.
