# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# CMake configuration for Palace unit tests
#

# Add Catch2
if(PALACE_BUILD_EXTERNAL_DEPS)
  Include(FetchContent)
  FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )
  FetchContent_MakeAvailable(Catch2)
else()
  find_package(Catch2 3 QUIET)
  if(NOT Catch2_FOUND)
    message(WARNING "Catch2 not found, unit-tests target will not be available")
    return()
  endif()
endif()

# Add executable target
add_executable(unit-tests
  ${CMAKE_CURRENT_SOURCE_DIR}/fixtures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-2d-submesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-basesolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-boundarymodeoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-coefficient.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-config.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-constants.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-domainpostoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-floquetport.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-geodata.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-helpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-libceed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-lumpedportintegration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-materialoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-memoryreporting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-nondimensionalize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-orthog.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-postoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-postoperatorcsv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-rap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-romoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-schema.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-strattonchu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-surfaceimpedanceoperator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-tablecsv.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-waveportimpedance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/regression_helpers.cpp
)

# Set output name for installed binary
set_target_properties(unit-tests PROPERTIES OUTPUT_NAME palace-unit-tests)

# Let regression/*.cpp find regression_helpers.hpp via a bare include.
target_include_directories(unit-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Below we link using WHOLE_ARCHIVE with coverage, otherwise LLVM won't report coverage of
# all palace files. The actual linker commands to do this by hand is delicate. (see
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_LINK_LIBRARY_USING_FEATURE.html#loading-a-whole-static-library)
target_link_libraries(unit-tests
  PRIVATE Catch2::Catch2
          $<$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>:coverage_flags>
          $<$<BOOL:${PALACE_BUILD_WITH_SANITIZERS}>:sanitizer_flags>
          $<IF:$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>,$<LINK_LIBRARY:WHOLE_ARCHIVE,${LIB_TARGET_NAME}>,${LIB_TARGET_NAME}>
          )

# Handle device source code
set(TARGET_SOURCES_DEVICE
  ${CMAKE_SOURCE_DIR}/fem/libceed/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-orthog.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test-vector.cpp
)
if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  if(PALACE_WITH_CUDA)
    set(LANGUAGE_PROPERTY CUDA)
  elseif(PALACE_WITH_HIP)
    set(LANGUAGE_PROPERTY HIP)
  endif()
  set(COMPILE_OPTIONS_PROPERTY "-Wno-pedantic")
  if(PALACE_WITH_OPENMP)
    set(COMPILE_OPTIONS_PROPERTY "${COMPILE_OPTIONS_PROPERTY} ${OpenMP_CXX_FLAGS}")
  endif()
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    PROPERTY LANGUAGE ${LANGUAGE_PROPERTY}
  )
  set_property(
    SOURCE ${TARGET_SOURCES_DEVICE}
    APPEND PROPERTY COMPILE_OPTIONS "${COMPILE_OPTIONS_PROPERTY}"
  )
endif()

# Resolve linker error with CUDA and static linkage
if(PALACE_WITH_CUDA)
  set_target_properties(unit-tests PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

# Add JIT source file path definition for libCEED
set_property(
  SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  APPEND PROPERTY COMPILE_DEFINITIONS "PALACE_LIBCEED_JIT_SOURCE_DIR=\"${CMAKE_INSTALL_PREFIX}/include/palace/\""
)

# Add path for test data files. 
target_compile_definitions(
  unit-tests
  PRIVATE "PALACE_TEST_DATA_DIR=\"${CMAKE_INSTALL_PREFIX}/share/palace/test/data\""
  # Source-tree roots used by the regression suite. Baking these into the
  # binary keeps the common path ("run the built binary") free of env
  # plumbing. RunRegressionCase honours a PALACE_EXAMPLES_DIR /
  # PALACE_REGRESSION_REF_DIR environment override and a Catch2 CLI flag
  # for relocated installs; these are the last-resort defaults.
  PRIVATE "PALACE_EXAMPLES_DIR_DEFAULT=\"${CMAKE_SOURCE_DIR}/../examples\""
  PRIVATE "PALACE_REGRESSION_REF_DIR_DEFAULT=\"${CMAKE_CURRENT_SOURCE_DIR}/../examples/ref\""
  PRIVATE $<TARGET_PROPERTY:${LIB_TARGET_NAME},COMPILE_DEFINITIONS>
          $<$<AND:$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:IntelLLVM>>>:PALACE_WITH_COVERAGE="LLVM">
          $<$<AND:$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>,$<CXX_COMPILER_ID:GNU>>:PALACE_WITH_COVERAGE="GCC">
)

# Status messages for test target
message(STATUS "Configured unit-tests target for unit tests using Catch2")

# Install Catch2 libraries
if(PALACE_BUILD_EXTERNAL_DEPS)
  install(
    TARGETS Catch2
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endif()

# Install unit-tests executable to bin directory
install(
  TARGETS unit-tests
  RUNTIME DESTINATION bin
)

# Install test data files to share directory
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
  DESTINATION share/palace/test/
)

# Install jsons from the all the examples (used in some tests)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/
  DESTINATION share/palace/test/data/examples
  FILES_MATCHING
    PATTERN "*.json"
    PATTERN "*/postpro/*" EXCLUDE     # Skip palace.json files
)

# Enable CTest for coverage support
include(CTest)

# Discover and register tests with CTest using Catch2's built-in function
include(Catch)

# Common coverage environment variable setup
set(SHOULD_SET_LLVM_PROFILE_FILE "$<AND:$<BOOL:${PALACE_BUILD_WITH_COVERAGE}>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:IntelLLVM>>>")

# Serial tests
catch_discover_tests(unit-tests
  TEST_PREFIX "serial-"
  # Exclude [Regression] so the fast unit-test sweep stays fast; regression
  # cases are invoked directly via `mpirun -np N palace-unit-tests
  # "[Regression]"` from the workflow (see
  # .github/workflows/build-and-test-*.yml).
  TEST_SPEC "[Serial]~[Regression]"
  EXTRA_ARGS --skip-benchmarks
  PROPERTIES
    SKIP_RETURN_CODE 4
    # When building with coverage, set LLVM_PROFILE_FILE for parallel execution
    ENVIRONMENT "$<${SHOULD_SET_LLVM_PROFILE_FILE}:LLVM_PROFILE_FILE=${CMAKE_BINARY_DIR}/coverage/serial/coverage-%p.profraw>"
)

# MPI tests  
catch_discover_tests(unit-tests
  TEST_PREFIX "mpi-"
  TEST_SPEC "[Parallel]~[Regression]"
  TEST_EXECUTOR ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
  EXTRA_ARGS ${MPIEXEC_POSTFLAGS} --skip-benchmarks
  PROPERTIES
    PROCESSORS 2
    SKIP_RETURN_CODE 4
    # When building with coverage, set LLVM_PROFILE_FILE for parallel execution
    ENVIRONMENT "$<${SHOULD_SET_LLVM_PROFILE_FILE}:LLVM_PROFILE_FILE=${CMAKE_BINARY_DIR}/coverage/mpi/coverage-%p.profraw>"
)

# Regression tests: full-mesh Palace solves diffed against test/examples/ref/.
# Regression tests: full-mesh Palace solves diffed against
# test/examples/ref/. Slow enough to deserve their own ctest prefix
# (`ctest -L "^regression$"` to select, `-LE "^(regression|long)$"` to
# skip both regression flavors).
#
# Rank count is a *runtime* knob, not a build option: each ctest entry
# launches via the generated `run_regression_test.sh` wrapper, which
# reads `$PALACE_REGRESSION_NUMPROC` (default below) and exec's mpirun.
# Use `PALACE_REGRESSION_NUMPROC=8 ctest -L "^regression$" -j 1` to fan
# out at 8 ranks per case; `PALACE_REGRESSION_NUMPROC=4 ctest -L
# "^regression$" -j 2` to run two 4-rank cases in parallel; etc.
set(PALACE_REGRESSION_NUMPROC_DEFAULT 2 CACHE STRING
    "Default MPI rank count for the [Regression] ctest registration. Overridden at run time by the PALACE_REGRESSION_NUMPROC environment variable.")

list(JOIN MPIEXEC_PREFLAGS " " MPIEXEC_PREFLAGS_STR)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/run_regression_test.sh.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_regression_test.sh
  @ONLY)
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/run_regression_test.sh
     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                 GROUP_READ GROUP_EXECUTE
                 WORLD_READ WORLD_EXECUTE)

catch_discover_tests(unit-tests
  TEST_PREFIX "regression-"
  TEST_SPEC "[Regression]~[Long]"
  TEST_EXECUTOR ${CMAKE_CURRENT_BINARY_DIR}/run_regression_test.sh
  EXTRA_ARGS ${MPIEXEC_POSTFLAGS}
  PROPERTIES
    LABELS regression
    SKIP_RETURN_CODE 4
    # Cylinder eigen cases can take ~50s on 2 ranks; give a generous ceiling.
    TIMEOUT 600
    ENVIRONMENT "$<${SHOULD_SET_LLVM_PROFILE_FILE}:LLVM_PROFILE_FILE=${CMAKE_BINARY_DIR}/coverage/regression/coverage-%p.profraw>"
)

# Long regression tests: full-mesh Palace solves that take 10+ minutes
# even at modest rank counts (transmon eigenmodes). Tagged
# `[Regression][Long]` and registered as `long-*` ctest entries with
# label `long`. Default `palace-unit-tests` and
# `ctest -L "^regression$"` runs skip them; the long-tests CI workflow
# (triggered by the `trigger-long-tests` PR label) selects them via
# `ctest -L "^long$"`.
catch_discover_tests(unit-tests
  TEST_PREFIX "long-"
  TEST_SPEC "[Long]"
  TEST_EXECUTOR ${CMAKE_CURRENT_BINARY_DIR}/run_regression_test.sh
  EXTRA_ARGS ${MPIEXEC_POSTFLAGS}
  PROPERTIES
    LABELS long
    SKIP_RETURN_CODE 4
    # Transmon eigen cases take ~10 minutes at np=2; allow 30 min headroom.
    TIMEOUT 1800
    ENVIRONMENT "$<${SHOULD_SET_LLVM_PROFILE_FILE}:LLVM_PROFILE_FILE=${CMAKE_BINARY_DIR}/coverage/regression/coverage-%p.profraw>"
)

# GPU tests (if available)
if(PALACE_WITH_CUDA)
  set(PALACE_GPU_DEVICE cuda)
elseif(PALACE_WITH_HIP)
  set(PALACE_GPU_DEVICE hip)
endif()

if(PALACE_WITH_CUDA OR PALACE_WITH_HIP)
  catch_discover_tests(unit-tests
    TEST_PREFIX "gpu-"
    # Exclude [Regression] so the GPU sweep stays focused on unit
    # tests; regression cases (including the [Long] transmon ones,
    # which also carry [Regression]) have their own registration and
    # are exercised on GPU via .github/actions/run-regression-tests.
    TEST_SPEC "[GPU]~[Regression]"
    EXTRA_ARGS --device ${PALACE_GPU_DEVICE} --skip-benchmarks
    PROPERTIES
      # Don't allow more than one GPU test at a time.
      # In theory, we could allow one per GPU, but that's much harder to implement.
      RESOURCE_LOCK gpu
      SKIP_RETURN_CODE 4

      # When building with coverage, set LLVM_PROFILE_FILE for parallel execution
      ENVIRONMENT "$<${SHOULD_SET_LLVM_PROFILE_FILE}:LLVM_PROFILE_FILE=${CMAKE_BINARY_DIR}/coverage/gpu/coverage-%p.profraw>"
  )

  # TODO: ADD MPI+GPU tests. This requires checking that we have enough GPUs
  # (one per rank).
endif()

# Enable testing
enable_testing()

