<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# Palace unit tests

This directory houses the Catch2-based test binary `palace-unit-tests`.
It serves two roles:

 1. **Unit tests** — fast, in-process checks of operators, solvers, and
    support code. Tagged `[Serial]` or `[Parallel]` (the latter is
    auto-applied when the binary is launched under MPI with more than
    one rank). Registered with CTest via `catch_discover_tests` and run
    as part of the default `ctest` sweep at hardcoded rank counts
    (1 for serial, 2 for parallel).

 2. **Regression cases** — full-mesh, full-solve runs against the
    `examples/` directory, with outputs diffed against
    `test/examples/ref/`. Tagged `[Regression]` (exclusively — no
    `[Serial]` / `[Parallel]` / `[GPU]` tag; the rank count is
    chosen at invocation time rather than baked into the case). These
    replace the Julia suite that used to live in `test/examples/`.

## Running

### Default — unit tests only

```
./palace-unit-tests                # serial unit tests
mpirun -n 2 ./palace-unit-tests    # parallel unit tests
```

`[Regression]` is excluded by default, both from the binary's own
filter and from the two `catch_discover_tests` registrations. No
opt-out needed.

### Regression cases

Regression cases are not registered with CTest. Invoke them
explicitly:

```
mpirun -n $NUM_PROC_TEST ./palace-unit-tests "[Regression]"
```

`NUM_PROC_TEST` is the same rank count the old Julia harness used
(computed in `.github/workflows/build-and-test-*.yml` as
`min(nproc, NUM_PROC_TEST_MAX)`, halved under coverage). Regression
reference CSVs were generated at that rank count, so running at a
different one will produce spurious mismatches on cases sensitive to
partition layout.

The regression cases are also registered as individual `regression-*`
CTest entries with a `regression` label, so they can be fanned out via
ctest. Each ctest entry launches one fresh `palace-unit-tests` process
through the `run_regression_test.sh` wrapper generated in the build
directory; the wrapper reads `$PALACE_REGRESSION_NUMPROC` for the rank
count (default 2) and `ctest -j N` controls how many cases run
concurrently. Total ranks in flight = `PALACE_REGRESSION_NUMPROC * N`.

```
# 8 cases at 2 ranks each, in parallel:
PALACE_REGRESSION_NUMPROC=2 ctest -L "^regression$" -j 8 --output-on-failure

# 1 case at 8 ranks at a time:
PALACE_REGRESSION_NUMPROC=8 ctest -L "^regression$" -j 1 --output-on-failure

# Just the cpw eigen cases:
PALACE_REGRESSION_NUMPROC=4 ctest -L "^regression$" -R cpw_.*_eigen

# Long regression cases (transmon eigenmodes, ~10 min each):
PALACE_REGRESSION_NUMPROC=2 ctest -L "^long$" -j 1 --output-on-failure
```

Tag combinations work as you'd expect. Because regression cases are
tagged `[Regression]` only, you filter them with the test name or a
wildcard spec rather than a secondary tag:

```
mpirun -n 4 ./palace-unit-tests "[Regression]"         # every regression
mpirun -n 2 ./palace-unit-tests "[Regression]~*cpw*"   # skip cpw cases
./palace-unit-tests rings                              # single case, 1 rank
```

## Where the cases live

All cases live in a single TU:

```
test/unit/regression/cases.cpp
```

Each `TEST_CASE` calls `palace::test::RunRegressionCase(...)` with the
case directory, config JSON, postpro subdirectory, and case-specific
tolerances / excluded-column list / custom per-CSV checks (for
phase-invariant eigenmode and farfield magnitude comparisons).

Each regression `TEST_CASE` carries `[Serial][Parallel][GPU][Regression]` —
`[Regression]` is the category tag, `[Serial]/[Parallel]/[GPU]` say the
case is valid under any rank count and on either CPU or GPU. The rank
count itself comes from the launch command (either the `mpirun -n N`
invocation for a direct selector run, or `$PALACE_REGRESSION_NUMPROC`
for the ctest-driven fan-out path). The two transmon eigen cases
additionally carry `[Long]`; they're registered as `long-*` ctest
entries (label `long`) and kept out of the default sweep.

## Overrides

The regression machinery reads the source-tree `examples/` and
`test/examples/ref/` roots via a two-level override chain:

| Precedence | Mechanism                                                                                                                                                                                 |
|:---------- |:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1 (high)   | `--examples-dir` / `--regression-ref-dir` / `--regression-run-dir` CLI flags                                                                                                              |
| 2 (low)    | Compile-time `PALACE_EXAMPLES_DIR_DEFAULT` / `PALACE_REGRESSION_REF_DIR_DEFAULT` (wired from CMake); run-dir falls back to `std::filesystem::temp_directory_path() / "palace-regression"` |

Solver knobs mirroring the old Julia flags:

| CLI flag                 | Effect                                     |
|:------------------------ |:------------------------------------------ |
| `--palace-linear-solver` | Overrides `Solver.Linear.Type` per case    |
| `--palace-eigensolver`   | Overrides `Solver.Eigenmode.Type` per case |
| `--palace-device`        | Overrides `Solver.Device` per case         |

Omitting `--palace-linear-solver` / `--palace-eigensolver` uses the old Julia harness default of `Default`. Individual cases preserve the Julia harness policy: most consume the global linear-solver knob, while cases that used `linear_solver="Default"` continue to force `Default` even in solver-matrix CI jobs.

## Adding a new regression case

 1. Drop the config under `examples/<name>/` and the reference postpro
    tree under `test/examples/ref/<name>/<subdir>/`.
 2. Add a `TEST_CASE("<name>", "[Serial][Parallel][GPU][Regression]")`
    to `test/unit/regression/cases.cpp`. Tack on `[Long]` if the case
    runs longer than the default `regression-` ctest TIMEOUT.
    to set `rtol`, `atol`, `excluded_columns`, `skip_rowcount`, and
    any `custom_checks` callbacks.
 3. Build `unit-tests`, then run
    `mpirun -n $NUM_PROC_TEST ./palace-unit-tests "<name>"`
    to validate.
