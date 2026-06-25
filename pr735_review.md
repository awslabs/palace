# Review: PR #735 — "Move the regression suite from Julia into Catch2"

**Author:** hughcars
**Branch:** `hughcars/regression-test-refactor` → `main` (stacked on `hughcars/boundarymode-2d`)
**Size:** 27 files, +1657 / −1547
**Head reviewed:** `6fa1bb46` · merge-base `5b947d9d`
**URL:** https://github.com/awslabs/palace/pull/735

Replaces the Julia regression harness (`test/examples/runtests.jl` + helpers) with
Catch2 cases built into `palace-unit-tests`, and lifts a `palace::Run` library entry
point out of `main.cpp` so production `main()` and the test binary share one code path.

Reviewed across four dimensions (separation/layering, code duplication, dead code,
clean code), each verified empirically against the PR head. One inter-reviewer
contradiction (the `MakeOutputFolder` double-call) was resolved by direct inspection —
see Mj1.

## Verdict: Request changes

Strong, well-structured refactor: the `palace::Run` extraction has correct layering
(test → library, never inverted), clean MPI ownership, a single source of truth for the
case list (a real improvement over the Julia harness), and faithful tolerance
preservation. But a regression *test suite* is held to a higher bar — "does it actually
compare data?" — and two findings fail it (B1, B2). There is also a verified
production-behavior regression (Mj1) and broken GitLab CI (Mj2).

* * *

## BLOCKERS — the suite can report PASS without testing

### B1. `TestFarfield` reads hardcoded column indices that match neither reference file

`test/unit/regression/cases.cpp:34-58` computes `|E|` from columns `t[3]..t[8]`, but the
actual reference files have different layouts:

  - `test/examples/ref/antenna/antenna_short_dipole/farfield-rE.csv`: `r*Re{E_x}` is at
    **index 4**, not 3 — off by one, folds the `phi` angle into the "magnitude."
  - `test/examples/ref/cpw/lumped_eigen/farfield-rE.csv`: `r*Re{E_x}` is at **index 6** —
    the check reads entirely wrong columns (`exc, theta, phi, …`).

Both farfield magnitude checks therefore compare meaningless, partition-independent
quantities and **likely pass spuriously**, defeating the one check meant to verify the
radiated-field magnitude.

**Fix:** look up columns by `header_text` (as `CompareComplexMagnitudes` already does),
or at minimum assert the header names at those indices.

### B2. No minimum-row guard → header-only / empty output passes vacuously

`regression_helpers.cpp` `CompareCSVFiles` loops `n_rows = min(actual.n_rows(), reference.n_rows())`. `ValidateCSVTables` only checks row-count equality when
`!opts.skip_rowcount`. For the **12 cases with `skip_rowcount=true`**, if `palace::Run`
writes a CSV with correct headers but **zero data rows** (e.g. an eigensolve converging 0
modes, or truncated output), `min=0` → the per-row loop never executes → the case
**PASSES with nothing ever compared**. `CheckSetEqual` only confirms the file *exists*.
The custom-check path (`CompareComplexMagnitudes`/`TestFarfield`) has the same hole.

**Fix:** after validation, `CHECK(reference.n_rows() > 0)` and require
`actual.n_rows() >= reference.n_rows()` (or `==`) even when `skip_rowcount` is set.

* * *

## MAJOR

### Mj1. `MakeOutputFolder` called twice in production (verified — resolves an inter-reviewer contradiction)

PR head `palace/main.cpp:228` calls `MakeOutputFolder(iodata, world_comm)`, then
`palace::Run` calls it again at `palace/driver.cpp:30` (invoked from `main.cpp:269`). The
old code called it exactly once (`5b947d9d:palace/main.cpp:239`).

Between the two calls, `main.cpp:234` writes `<output>/<stem>_resolved.json` into the
directory, so the second call (inside `Run`) finds the directory non-empty and emits a
spurious **"Output folder is not empty; program will overwrite content!"** warning on
*every* production run, plus redundant path re-validation and a re-run write-probe.
Undocumented behavior change (CHANGELOG mentions only the test-harness move).

**Fix:** single owner — drop the `main.cpp:228` call and let `Run` own output-dir
creation, reconciling the `WriteResolvedConfig` ordering (it needs the dir to exist).

### Mj2. `.gitlab-ci.yml` still runs the deleted Julia harness

`.gitlab-ci.yml` is **not touched by this PR**, yet 12 build/test jobs each run:

```
julia --project=test/examples -e 'using Pkg; Pkg.instantiate()'
julia --project=test/examples --color=yes test/examples/runtests.jl
```

(lines 52-53, 74-75, 92-93, 108-109, 121-122, 134-135, 148-149, 163-164, 181-182,
196-197, 212-213, 233-234). Both `test/examples/Project.toml` and
`test/examples/runtests.jl` are deleted on this branch, so every GitLab test job breaks
at `Pkg.instantiate()` / file-not-found.

**Action:** port `.gitlab-ci.yml` to the ctest model, or confirm the GitLab pipeline is
retired (if dormant, downgrade to MINOR but still an incomplete-removal loose end).

### Mj3. Staging failures deadlock all ranks until the CTest timeout

`regression_helpers.cpp` `ScopedExampleStage` ctor performs throwing `create_directories`
/ `directory_iterator` / `create_symlink` on rank 0 *before* the `Mpi::Barrier(comm)`. If
rank 0 throws (permissions, a stale root-owned `palace-regression`, ENOSPC), it unwinds
out of the ctor while non-root ranks block forever at the barrier → the job hangs until
the 600s/1800s CTest `TIMEOUT` instead of failing fast. `remove_all(stage_, ec)` also
ignores `ec`.

**Fix:** stage with `error_code` overloads, broadcast success/failure to all ranks, and
`REQUIRE` symmetrically so a staging error fails cleanly on every rank.

### Mj4. Two live execution paths will drift

There are two independent reimplementations of "run the regression suite":

  - **CTest path** (`ctest -L "^regression$"`, used by linux/macos workflows + `make tests`): rank count from `PALACE_REGRESSION_NUMPROC`, case set from the CMake
    `TEST_SPEC`.
  - **Direct-binary path** (`.github/actions/run-regression-tests/action.yml`): rank count
    from inline `NUM_PROC_TEST`, case set re-derived via
    `palace-unit-tests "[Regression]~[Long]" --list-tests | awk '/^  [a-z]/ {print $1}'`.

They agree today but nothing enforces it; adding a case / changing a tag / changing rank
count requires edits in multiple places. The awk scrape is also fragile — a future
Catch2 list-format change (or an uppercase/digit-leading case name) yields **zero cases
with exit 0** (silent pass), with no "ran 0 cases" guard.

This divergence is *not* forced by the spack context. The `run-regression-tests` action
can run from the build tree (it already has a `testing-on-build-runner` input for the
build-machine case), so `ctest -L "^regression$"` and the generated
`run_regression_test.sh` wrapper are available to it just as they are to the
linux/macos workflows. The direct `--list-tests | awk | mpirun` loop reimplements — more
fragilely — what the CTest registration already provides (per-case process isolation,
rank count via `PALACE_REGRESSION_NUMPROC`, per-case `TIMEOUT`, `SKIP_RETURN_CODE`
handling). The only thing the action does that the CTest wrapper doesn't is inject the
runtime solver/device override flags (`--palace-linear-solver`, `--palace-device GPU`);
that can ride an env var the wrapper appends, exactly like `PALACE_REGRESSION_NUMPROC`
already does for rank count.

**Fix:** collapse the action onto CTest — have `run_regression_test.sh.in` append an
optional `$PALACE_SOLVER_OVERRIDE`, then replace the action's `translate_cases` /
`run_cases` loop with `PALACE_REGRESSION_NUMPROC=$NUM_PROC_TEST ctest -L "^regression$"`
(CUDA → `PALACE_REGRESSION_NUMPROC=1` + the device override). All platforms then drive
the suite identically and the brittle scrape is removed.

* * *

## MINOR / NIT

  - **NaN/Inf and malformed CSV cells** *(clean-code, verified-latent)*: tablecsv skips
    empty/`"NULL"` cells, under-populating a column → unchecked `data[r]` out-of-bounds
    read (UB) rather than a clean failure; a non-numeric cell triggers `MFEM_ABORT` (MFEM
    built with `MFEM_USE_EXCEPTIONS=NO`) which aborts the process instead of producing a
    Catch2 failure. No current fixture triggers it.
  - **Tolerances / excluded-columns scattered as inline magic literals** *(duplication,
    MINOR)*: `opts.rtol`/`opts.atol` repeated ~30× each; `{"Maximum","Minimum"}` repeated
    11×; the eigen exclusion list hand-rolled 5× *despite* `kEigenExcluded` existing. The
    deleted Julia harness cascaded shared `reltol`/`abstol` variables, so this is a DRY
    regression. Add a tolerance-tier table + named column-set constants
    (`kMinMaxExcluded`, `kPortColumns`) — the single highest-value cleanup.
  - **`TestFarfield` is a subset of `CompareComplexMagnitudes`** *(duplication, MINOR)*:
    both implement the `sqrt(Re²+Im²)` magnitude compare (~25 duplicated lines); the
    general one could subsume the farfield one.
  - **Shared-build library path** *(CI, MINOR-suspected)*: the regression CI steps don't
    set `LD_LIBRARY_PATH` / `DYLD_LIBRARY_PATH` for the build-tree binary, while the
    sibling unit-test step does. The old code ran the *installed* binary; could fail to
    dlopen on the `shared` matrix entry.
  - **`BlockTimer::Reset()` comment overstated** *(layering, MINOR)*: peak memory comes
    from `getrusage(...).ru_maxrss` (monotonic, unresettable), so absolute peak *does*
    bleed across in-process cases — only per-phase deltas are correct. Tighten the comment.
  - **Misleading `-j` fan-out comments** *(CI, NIT)*: both workflows' comments describe a
    CTest concurrency model the commands don't use (no `-j` flag).
  - **Unused `<fmt/format.h>`** in `test/unit/regression_helpers.cpp:16` *(dead-code,
    NIT)*.
  - **`LoadConfigJson` bypasses `ParseAndValidate`** *(duplication, NIT)*: re-derives
    `PreprocessFile` + `json::parse`, skipping the duplicate-key / schema validation.
    Intentional (needs the mutable json to inject overrides) but undocumented.
  - **`Run` returns `void`** *(layering, MINOR)*: signals failure only via `MFEM_ABORT`
    (process abort), so a bad case kills the whole test binary rather than failing one
    case. Acceptable given upstream validation, but a CLI-shaped assumption baked into a
    library API.

* * *

## PR-description discrepancies to correct

  - Body says **27 cases**; the branch has **30** `TEST_CASE`s.
  - Body says **transmon not ported**; `transmon_coarse` + `transmon_amr` *are* registered
    (`[Regression][Long]`, `cases.cpp:373,388`) with full fixtures. Only `adapter/slp` is
    genuinely omitted (correctly — no fixture exists).

* * *

## Verified-clean (do not re-litigate)

  - **MPI lifecycle**: single owner. `MPI_Init` only in each `main` (production
    `main.cpp:158`, test `test/unit/main.cpp:28`, both guarded); `MPI_Finalize` via the
    `Mpi` singleton destructor (idempotent). `palace::Run` never touches Init/Finalize —
    correct for repeated in-process calls.
  - **No production→test dependency inversion**: `driver.cpp` includes only `palace/`
    headers; CMake puts it in the library target, `main.cpp` in the executable; both the
    executable and test binary link the library. Flow is test → library.
  - **Single source of truth for the case list**: only the `TEST_CASE` declarations; CI and
    CMake derive dynamically. Eliminates the Julia harness's duplicated `TEST_CASES` array.
  - **`make tests`, `docs/src/developer/testing.md`, `scripts/measure-test-coverage`, and
    the GitHub workflows/actions** are all correctly rewired (no stale Julia references on
    the GitHub side).
  - **Antenna/Julia dependency** is genuinely still needed on exactly one Linux matrix entry
    (x86 + gcc + aocl) for `plot_farfield.jl`; its staged path matches the harness output
    and the stage is intentionally not torn down. macOS correctly drops `setup-julia`.
  - **`[Regression]` excluded from fast PR runs** everywhere (catch_discover_tests specs,
    bare-invocation auto-filter, `ctest -LE`), so PRs don't pay ~39 min unless the
    regression step runs. Build-before-test ordering correct; timeouts have headroom.
  - **`run_regression_test.sh.in`** is used (configured + set as `TEST_EXECUTOR`), not dead
    scaffolding.
  - **`main.cpp` ↔ `driver.cpp` extraction**: the lifted block was deleted from `main.cpp`,
    not duplicated; includes trimmed; no residual forked setup.
  - **`timer.hpp` `start_time` const→non-const**: required for `Reset()`'s `timer = Timer()`; correct and minimal.

* * *

## Priority to address

 1. **B1 + B2** — false passes; the entire purpose of the suite.
 2. **Mj1** — production behavior regression (spurious warning every run).
 3. **Mj2** — broken GitLab CI (confirm whether the pipeline is active).
 4. **Mj3 + Mj4** — robustness and execution-path consolidation.
 5. The rest is cleanup.
