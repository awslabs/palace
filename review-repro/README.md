# Reproducers for the PR 992 review (TwoSided London sheets)

These configs come from review of `dnpham23/london-two-port-sheet` at
`6c71c655f8`. They are review artifacts, not regression tests, and are not meant
to be merged.

Run everything from the **repository root**: mesh and output paths are relative
to it. Outputs go to `review-repro/postpro/<name>/` (git-ignored). All numbers
below were measured on darwin-m1 with a spack build of the PR head.

```bash
palace -np 4 review-repro/<config>.json
```

## Meshes

- `washer_*`: `test/data/regression/input/circular_hole_london/mesh/circular_hole.msh`,
  the existing regression mesh. Nothing to generate.
- `pec_*`: `review-repro/pec.msh` is not committed (about 3 MB). Generate it with
  gmsh's Python API (4.15.2), which gives deterministic output:
  ```bash
  python3 review-repro/make_pec_mesh.py   # writes review-repro/pec.msh
  ```
  The geometry is a 10 µm PEC box (attribute 2, mirror-symmetric about z = 0).
  An interior film on z = 0 (attribute 8) spans the whole cross-section, so its
  outer edge lies on the PEC walls. It has a radius-1.5 µm hole whose cap is
  attribute 9. Both sides of the film are volume, so the "film must be fully
  interior" requirement is met. This is the "shared ground plane" case the docs
  advertise.

Throughout, λ = 0.1 µm, and d = 0.2 µm unless stated. The `*_sym_reference`
configs replace the two-sided film with a single sheet of
`KineticInductance = ½ μ0 λ coth(d/2λ)`. Both geometries are mirror-symmetric,
so the field is the same on both faces and the two-sided result must match that
reference.

## Configs

### `pec_two_sided.json`: blocker, C not eliminated on essential DOFs
- **Shows:** the cross-face coupling C is added to the operator after K's PEC
  (essential) DOFs are eliminated. The rows and columns of C on the film's PEC
  edge are therefore never eliminated.
- **Observed (np 4):** `PCG solver did NOT converge in 200 iterations`
  (stalls at a relative residual of 9.15e-08 against `Tol` 1e-10). The run
  aborts with "London solve ... did not converge". No inductance is written.
- **Expected:** converges in about 20 iterations to L ≈ 1.6886e-13 H, matching
  `pec_sym_reference.json`.

### `pec_two_sided_loose_tol.json`: the same blocker, silently wrong answer
- **Shows:** `pec_two_sided.json` with `Tol` 1e-6, so CG reports convergence.
- **Observed (np 4):** converged in 88 iterations, **L = 1.257406507288e-13 H**.
  That is 25.5% below the correct value, with no warning.
- **Expected:** L ≈ 1.6886e-13 H.

### `pec_sym_reference.json`: reference for the two configs above
- **Shows:** the equivalent single sheet with `KineticInductance` 8.250044e-14 H/sq.
- **Observed (np 4):** 20 iterations, **L = 1.688602505288e-13 H**.
- With the reference fix on this branch (commit 2), `pec_two_sided.json`
  converges in 22 iterations to L = 1.68880e-13 H. That is within 1.2e-4 of
  this reference; the mesh is not exactly mirror-symmetric.

### `washer_two_sided_nc_amr.json`: TwoSided + nonconformal AMR aborts
- **Shows:** `BuildTwoPortCoupling` maps local DOFs to global true DOFs with
  `ParFiniteElementSpace::GetGlobalTDofNumber`. On a nonconforming mesh that
  only works for owned, non-slave DOFs.
- **Observed:**
  - np 2 aborts on AMR iteration 0, before any refinement:
    `Verification failed: (ldof_ltdof[ldof] >= 0) ... ldof 422 not a true DOF.`
    The mesh is already converted to NC in `Partition`, and the film has shared DOFs.
  - np 1 solves iteration 0 (L as in the regression case), then aborts on
    iteration 2 with `ldof 18140 not a true DOF` (slave DOFs).
- **Expected:** either a correct solve, with C assembled through the
  prolongation P (PᵀC P), or a clear config-time rejection of
  `TwoSided` + `Refinement.MaxIts > 0`.

### `washer_two_sided_conforming_amr.json`: TwoSided + conforming AMR aborts
- **Shows:** the two sides of the crack are marked and refined independently,
  so after refinement the film triangles on the two faces no longer match.
- **Observed:** iteration 1 solves. After refinement (494 elements marked,
  73750 total), iteration 2 aborts with
  `Two-sided (two-port) sheet has 384 unpaired face(s)` at np 2. The count is
  192 at np 1 and 768 at np 4, because the unpaired count is also multiplied by
  the number of ranks. The message ("the film must be fully interior") is
  misleading here.
- **Expected:** a correct solve, or a clear config-time rejection.
- For comparison, uniform refinement (`Refinement.UniformLevels: 1`) works:
  L = 3.587255e-12 H.

### `washer_two_sided_thin_d0.005.json`: preconditioner degrades as d/λ → 0
- **Shows:** a thin film (d = 0.005 µm, d/λ = 0.05) on the regression washer,
  with `MaxIts` raised to 2000. The operator K + C is preconditioned by K
  without C, which misjudges the soft common mode by about 2λ²/d².
- **Observed (np 4):** correct L = 1.781158867432e-11 H, but **310 CG
  iterations**. The reference takes 84. At d = 0.02 µm it is 133 against 61.
  Palace's default `MaxIts` of 100 would fail.
- **Expected:** iterations comparable to the single-sheet reference. C is
  already an assembled `HypreParMatrix`, so it could be added to the finest
  assembled preconditioner matrix. Alternatively, warn when d/λ is small.

### `washer_sym_reference_thin_d0.005.json`: reference for the thin case
- **Observed (np 4):** 84 iterations, L = 1.781158824456e-11 H.

## Not reproduced here (found by reading the code)

Screened current-port steps create a temporary `step_op` (K_step + C) for each
step. The driver's `set_operator` decides whether to rebind by comparing
addresses. If a freed `step_op` address is reused by the next step with a
different short key, rebinding is skipped and that step runs with the previous
key's preconditioner. Those screened steps also do not eliminate C for the extra
shorted-port DOFs.
