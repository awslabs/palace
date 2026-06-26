# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Relocatable CTest entry point for an *installed* palace-unit-tests (no build
# tree). `catch_discover_tests` only writes build-tree CTestTestfiles, so the
# install-only Spack test jobs have nothing for `ctest` to run. This file is
# installed verbatim next to the test data as CTestTestfile.cmake, so those
# jobs can:
#
#   ctest --test-dir <prefix>/share/palace/test -L "^regression$"
#
# It is fully generic -- nothing is baked at build time, so it never needs
# regeneration. Everything is taken at `ctest` time:
#
#   * the binary is found on PATH (e.g. via `spack load palace`);
#   * the case list is discovered from the binary itself;
#   * ranks / threads / solver / device come from the environment.
#
# Knobs (all optional):
#   PALACE_TEST_NUMPROC        MPI ranks per case          (default 2)
#   PALACE_TEST_OMP_THREADS    OpenMP threads per process  (default 1)
#   PALACE_TEST_LINEAR_SOLVER  Solver.Linear.Type override (default: case JSON)
#   PALACE_TEST_DEVICE         MFEM device, e.g. cuda      (default: cpu)
#
# Each case reserves PALACE_TEST_NUMPROC * PALACE_TEST_OMP_THREADS CTest slots
# (PROCESSORS) and runs with OMP_NUM_THREADS pinned to the same thread count, so
# the slot accounting is derived from a single source and cannot drift.

find_program(_palace_bin NAMES palace-unit-tests)
if(NOT _palace_bin)
  set(_palace_bin palace-unit-tests)
endif()
find_program(_palace_mpiexec NAMES mpirun mpiexec)
if(NOT _palace_mpiexec)
  set(_palace_mpiexec mpirun)
endif()

set(_palace_np "$ENV{PALACE_TEST_NUMPROC}")
if(NOT _palace_np)
  set(_palace_np 2)
endif()
set(_palace_omp "$ENV{PALACE_TEST_OMP_THREADS}")
if(NOT _palace_omp)
  set(_palace_omp 1)
endif()
math(EXPR _palace_procs "${_palace_np} * ${_palace_omp}")

# Per-run solver/device, appended to each case command only (not to discovery,
# so listing never needs a GPU).
set(_palace_args "")
if(NOT "$ENV{PALACE_TEST_LINEAR_SOLVER}" STREQUAL "")
  list(APPEND _palace_args "--palace-linear-solver" "$ENV{PALACE_TEST_LINEAR_SOLVER}")
endif()
if(NOT "$ENV{PALACE_TEST_DEVICE}" STREQUAL "")
  list(APPEND _palace_args "--device" "$ENV{PALACE_TEST_DEVICE}")
endif()

set(_palace_tmp "$ENV{TMPDIR}")
if(NOT _palace_tmp)
  set(_palace_tmp "/tmp")
endif()

# On a device (GPU) run, only register cases tagged [GPU]: palace-unit-tests
# intersects [GPU] when --device is set, so a non-GPU case would match nothing
# and exit non-zero. Mirror that filter at discovery so it isn't registered.
set(_palace_tag_filter "")
if(NOT "$ENV{PALACE_TEST_DEVICE}" STREQUAL "")
  set(_palace_tag_filter "[GPU]")
endif()

function(_palace_register spec label)
  # Discover via a file, not OUTPUT_VARIABLE: ctest's script interpreter
  # truncates large execute_process captures, but file output is complete.
  set(_xml_file "${_palace_tmp}/palace-ctest-${label}.xml")
  execute_process(
    COMMAND "${_palace_bin}" "${spec}" --list-tests --reporter xml
    OUTPUT_FILE "${_xml_file}"
    RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "Failed to list '${spec}' cases from ${_palace_bin} (${_rc})")
  endif()
  file(READ "${_xml_file}" _xml)
  string(REGEX MATCHALL "<Name>[^<]*</Name>" _names "${_xml}")
  foreach(_tag IN LISTS _names)
    string(REGEX REPLACE "<Name>(.*)</Name>" "\\1" _case "${_tag}")
    add_test("${label}-${_case}"
      "${_palace_mpiexec}" "-n" "${_palace_np}" "${_palace_bin}" "${_case}" ${_palace_args})
    set_tests_properties("${label}-${_case}" PROPERTIES
      LABELS "${label}"
      PROCESSORS "${_palace_procs}"
      SKIP_RETURN_CODE 4
      ENVIRONMENT "OMP_NUM_THREADS=${_palace_omp}")
  endforeach()
endfunction()

_palace_register("[Regression]${_palace_tag_filter}~[Long]" "regression")
_palace_register("[Long]${_palace_tag_filter}" "long")
