# Palace Example Tests

This directory contains regression tests for the Palace examples. The tests
validate the correctness of Palace simulations by comparing output results
against reference data.

## Running Tests

### Prerequisites

  - Julia
  - Palace executable in PATH or specified via environment variable

### Setup

First, instantiate the Julia environment:

```bash
julia --project -e "using Pkg; Pkg.instantiate()"
```

You need to do this step only the very first time.

### Environment Variables

  - `PALACE_TEST`: Path to Palace executable and optional arguments (default: "`palace`")
  - `NUM_PROC_TEST`: Number of MPI processes (default: number of physical cores)
  - `OMP_NUM_THREADS`: Number of OpenMP threads (default: 1)
  - `TEST_CASES`: Space-separated list of test cases to run (default: all examples)

### Execution

Run all tests:

```bash
julia --project runtests.jl
```

Run specific test cases:

```bash
TEST_CASES="spheres rings" julia --project runtests.jl
```

Run with custom Palace executable:

```bash
PALACE_TEST="../../build/bin/palace" julia --project runtests.jl
```

Run with custom Palace executable and custom MPI launcher:

```bash
PALACE_TEST="../../build/bin/palace --launcher mpiexec" julia --project runtests.jl
```

Each test case runs Palace simulations and compares generated CSV files against
reference data using configurable tolerances. When Palace behavior changes
legitimately (e.g., algorithm improvements), reference data can be updated using
the baseline script:

```bash
./baseline                    # Update all reference data
./baseline -e spheres         # Update specific example
./baseline --dry-run          # Test without updating files
./baseline -np 4              # Use 4 MPI processes
```
