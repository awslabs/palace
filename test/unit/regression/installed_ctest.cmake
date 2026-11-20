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
# The "long" label also runs the unit tests on MFEM's debug device, on 1 and 2 processes
# (see the end of this file).
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
# Each case reserves its number of processes (PALACE_TEST_NUMPROC for the regression
# and long cases) times PALACE_TEST_OMP_THREADS CTest slots (PROCESSORS) and runs with
# OMP_NUM_THREADS pinned to the same thread count, so the slot accounting is derived
# from a single source and cannot drift.

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

# Per-run solver/device, appended to each case command only (not to discovery,
# so listing never needs a GPU).
set(_palace_solver_args "")
if(NOT "$ENV{PALACE_TEST_LINEAR_SOLVER}" STREQUAL "")
  list(APPEND _palace_solver_args "--palace-linear-solver" "$ENV{PALACE_TEST_LINEAR_SOLVER}")
endif()
set(_palace_args ${_palace_solver_args})
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

# Register the cases matching SPEC as "<PREFIX>-<case>" tests with label LABEL, each run on
# NUMPROC processes with the extra ARGS.
function(_palace_register)
  cmake_parse_arguments(PARSE_ARGV 0 arg "" "SPEC;LABEL;PREFIX;NUMPROC" "ARGS")
  math(EXPR _procs "${arg_NUMPROC} * ${_palace_omp}")
  # Discover on the processes the cases run on, as catch_discover_tests does: the binary
  # selects [Serial] cases on one process and [Parallel] ones on more. Discover via a file,
  # not OUTPUT_VARIABLE: ctest's script interpreter truncates large execute_process
  # captures, but file output is complete.
  set(_xml_file "${_palace_tmp}/palace-ctest-${arg_PREFIX}.xml")
  execute_process(
    COMMAND "${_palace_mpiexec}" "-n" "${arg_NUMPROC}" "${_palace_bin}" "${arg_SPEC}"
      --list-tests --reporter xml
    OUTPUT_FILE "${_xml_file}"
    RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "Failed to list '${arg_SPEC}' cases from ${_palace_bin} (${_rc})")
  endif()
  file(READ "${_xml_file}" _xml)
  string(REGEX MATCHALL "<Name>[^<]*</Name>" _names "${_xml}")
  set(_seen "")
  foreach(_tag IN LISTS _names)
    string(REGEX REPLACE "<Name>(.*)</Name>" "\\1" _case "${_tag}")
    # Cases with the same name run together, as one test. (No IN_LIST: ctest reads this file
    # without a policy version.)
    list(FIND _seen "${_case}" _index)
    if(NOT _index EQUAL -1)
      continue()
    endif()
    list(APPEND _seen "${_case}")
    # Escape the characters that Catch2 parses in a test spec, as CatchAddTests.cmake does
    # (backslash first).
    set(_spec "${_case}")
    foreach(_char "\\" "," "[" "]")
      string(REPLACE "${_char}" "\\${_char}" _spec "${_spec}")
    endforeach()
    add_test("${arg_PREFIX}-${_case}"
      "${_palace_mpiexec}" "-n" "${arg_NUMPROC}" "${_palace_bin}" "${_spec}" ${arg_ARGS})
    set_tests_properties("${arg_PREFIX}-${_case}" PROPERTIES
      LABELS "${arg_LABEL}"
      PROCESSORS "${_procs}"
      SKIP_RETURN_CODE 4
      ENVIRONMENT "OMP_NUM_THREADS=${_palace_omp}")
  endforeach()
endfunction()

_palace_register(SPEC "[Regression]${_palace_tag_filter}~[Long]" LABEL regression
  PREFIX regression NUMPROC ${_palace_np} ARGS ${_palace_args})
_palace_register(SPEC "[Long]${_palace_tag_filter}" LABEL long
  PREFIX long NUMPROC ${_palace_np} ARGS ${_palace_args})

# Run the unit tests on MFEM's debug device with the long cases, always on the CPU
# (PALACE_TEST_DEVICE does not apply). Without a GPU, the debug device fails on a host access
# to data whose valid copy is on the device (operator[] and operator() do not synchronize),
# which on a GPU would read stale data, and it runs libCEED on its CPU backend.
foreach(_palace_debug_np 1 2)
  _palace_register(SPEC "~[Regression]~[Long]" LABEL long
    PREFIX debug-np${_palace_debug_np} NUMPROC ${_palace_debug_np}
    ARGS ${_palace_solver_args} --device debug --skip-benchmarks)
endforeach()
