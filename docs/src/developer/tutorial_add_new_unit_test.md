```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Tutorial: Adding a New Unit Test

This tutorial demonstrates how to add unit tests to the *Palace* test suite.
Prerequisites: familiarity with compiling and running unit tests (see [Running
unit tests](testing.md#Running-unit-tests)).

## Overview

As a motivating example, we'll test the `Sum(MPI_Comm comm, Vector& vec)`
function. This function has certain characteristics that make it interesting to
test:

  - `Vector`s may have different sizes across MPI processes
  - `Vector` data may be stored on CPU or GPU
  - The function performs MPI communication to compute global sums

These characteristics require testing across multiple configurations to ensure
correctness.

Create a new file `test-vector-sum.cpp` in `test/unit/`:

```cpp
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_test_macros.hpp>

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
using namespace Catch;

TEST_CASE("Vector Sum - Basic", "[myvector][Serial]")
{
  Vector v(2);
  v(0) = 1.0;
  v(1) = 2.0;

  double sum = linalg::Sum(Mpi::World(), v);
  REQUIRE_THAT(sum, Catch::Matchers::WithinRel(3.0));
}

}  // namespace palace
```

This defines a new `Catch2` test case. The key components are:

  - Test name (`Vector Sum - Basic`): Must be unique across the test suite
  - Tags: `[myvector]` (arbitrary, used for filtering) and `[Serial]` (special tag, more on this later)
  - `WithinRel()`: Handles floating-point comparison tolerances

To compile our test, we need to add it to the list of sources in the `CMakeLists.txt` in `test/unit`:

```cmake
add_executable(unit-tests
  # ... existing files ...
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector-sum.cpp
)
```

Then, build and run:

```bash
# Build tests in the build directory
make palace-tests

# Run this specific test
bin/palace-unit-tests "[myvector]"
```

We should see that `All tests passed`.

The test that we wrote is a reasonable first test, but it does not check that
the feature works with GPUs or with multiple processes.

In particular, the test fails with multiple MPI processes because each process
has a copy of `v` and contributes to the sum. Similarly, if we were to run this
on a GPU we would find that the test passes, but we would also observe no
activity on the device.

Let us extend this case so that we can write a more comprehensive test that is
also meaningful on GPU and with MPI:

```cpp
TEST_CASE("MyTest Vector Sum", "[myvector][Serial][Parallel][GPU]")
{
  Vector v(2);
  v.UseDevice(true);
  auto d_v = v.Write();

  mfem::forall(v.Size(), [=] MFEM_HOST_DEVICE (int i){
    d_v[i] = rank + 1.0;
  });

  double sum = linalg::Sum(Mpi::World(), v);
  double expected = Mpi::Size(Mpi::World()) * 3.0;
  REQUIRE_THAT(value, Catch::Matchers::WithinRel(expected))
}
```

We added GPU compatibility by (see, [MFEM documentation](https://mfem.org/gpu-support/)):

 1. Adding `v.UseDevice(true);`, which defines our intent to use `v` for computations on the device.
 2. Defining `auto d_v = v.Write();`, a pointer to area of memory on the device.
 3. Using `forall` to execute the execute the function on the device.

When we add GPU code, we need to make sure that the file is compiled with the
correct compiler. To do so, we add the file to `TARGET_SOURCES_DEVICE` in the
`CMakeLists.txt` file, it is not already there.

```cmake
set(TARGET_SOURCES_DEVICE
  # ... existing files ...
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector-sum.cpp
)
```

Note that this code still works for CPU, when the device is not a GPU, allowing
us to test both devices with the same code. Note also that here we are working
with a vector with only two elements and this is highly inefficient on GPUs (but
we are not concerned with performance here, only correctness).

To add MPI compatibility, we changed the `expected` value to account for how
many copies of the vector `v` there are.

With these additions, this test case is a meaningful and interesting test in all
the possible configurations, so we we added the `[Parallel]` and `[GPU]` tags.

*Palace* uses three special tags for execution control:

  - `[Serial]`: Runs only with single MPI process
  - `[Parallel]`: Runs only with multiple MPI processes
  - `[GPU]`: Runs only when GPU devices are available

To understand why we need this, let us add a test that checks that the sum is
correct when vectors have different lengths on different MPI processes (for
simplicity let us ignore GPU compatibility):

```cpp
TEST_CASE("MyTest Vector Sum - Different Lengths", "[myvector][Parallel]")
{
  Vector v;
  
  if (Mpi::Root(Mpi::World())){
    v.SetSize(2);
    v(0) = 10;
    v(1) = 20;
  } else {
    v.SetSize(1);
    v(0) = 3;
  }
  
  double sum = linalg::Sum(Mpi::World(), v);
  double expected = 3 * (Mpi::Size(Mpi::World()) - 1) + 30;
  REQUIRE_THAT(value, Catch::Matchers::WithinRel(expected))
}
```

This test is useful because it checks that `Sum` is not implemented making
assumptions on the length of the vector. This test is also meaningless when run
with less than 2 MPI processes, so we removed the `[Serial]` tag.

Sometimes, tests need to write to the filesystem. In this case, it is often best
to create temporary working directories. This can be accomplished with the
`PerRankTempDir` fixture, which gives each MPI rank its own directory:

```cpp
#include "fixtures.hpp"

TEST_CASE_METHOD(palace::test::PerRankTempDir, "MyTest Print", "[myvector][Serial]") {
  // temp_dir is available and will be cleaned up automatically.
  auto file_path = temp_dir / "vector.txt";
  Vector v;
  v = 1;
  {
     std::ofstream file(file_path);
     v.Print(file);
  }
  CHECK(std::filesystem::exists(file_path));
}
```

For tests where all ranks need to share the same directory, use `SharedTempDir`:

```cpp
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "MyTest Shared", "[myvector][Parallel]") {
  // All ranks share temp_dir.
}
```

Suppose you want to compare the result of some operations with some pre-existing
file `expected_vector.txt`. To do this, we first need to save the file in
`test/unit/data`. Then, we can access it as

```cpp
auto path_expected_vector = fs::path(PALACE_TEST_DATA_DIR) / "expected_vector.txt"
```

`path_expected_vector` points to the `expected_vector.txt` file.
