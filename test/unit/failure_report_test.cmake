# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Runs the hidden "Failure report on abort" test case, which aborts on purpose, and checks
# the palace.json it leaves. Arguments (-D): EXE, MPIEXEC, NP_FLAG, NP, PREFLAGS, POSTFLAGS,
# MODE (<root|last|root-last|all>[-dtor]-<abort|throw>) and DIR. The test case waits for
# 1 s for all processes to fail.

file(REMOVE_RECURSE "${DIR}")
file(MAKE_DIRECTORY "${DIR}")
string(TIMESTAMP t0 "%s" UTC)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E env PALACE_TEST_FAILURE=${MODE} PALACE_TEST_FAILURE_DIR=${DIR}
          ${MPIEXEC} ${NP_FLAG} ${NP} ${PREFLAGS} ${EXE} ${POSTFLAGS}
          "Failure report on abort"
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE out
  TIMEOUT 300
)
string(TIMESTAMP t1 "%s" UTC)
math(EXPR elapsed "${t1} - ${t0}")
if(NOT rc MATCHES "^[0-9]+$" OR rc EQUAL 0)
  message(FATAL_ERROR "The run did not abort (${rc}):\n${out}")
endif()
if(NOT EXISTS "${DIR}/palace.json")
  message(FATAL_ERROR "No palace.json (exit code ${rc}):\n${out}")
endif()
file(READ "${DIR}/palace.json" meta)

function(expect path expected)
  string(REPLACE "/" ";" keys "${path}")
  string(JSON value ERROR_VARIABLE err GET "${meta}" ${keys})
  if(NOT "${value}" STREQUAL "${expected}")
    message(FATAL_ERROR "${path} is \"${value}\" (${err}), expected \"${expected}\":\n${meta}\n${out}")
  endif()
endfunction()

if(MODE MATCHES "^all" OR (MODE MATCHES "^root-last" AND NP EQUAL 2))
  set(failed all)
elseif(MODE MATCHES "^last")
  set(failed last)
else()
  set(failed root)
endif()
if(failed STREQUAL "last")
  # Only a non-root process fails: palace.json is left as it was, and the process waits
  # for 3 s (for the root to fail too) before aborting.
  string(JSON value ERROR_VARIABLE err GET "${meta}" Complete)
  if(NOT err)
    message(FATAL_ERROR "Unexpected report:\n${meta}")
  endif()
  expect(Problem/MPISize ${NP})
  if(elapsed LESS 3)
    message(FATAL_ERROR "Aborted after ${elapsed} s, before the end of the wait:\n${out}")
  endif()
else()
  expect(Complete OFF)
  expect(Error/Phase EPS)
  if(MODE MATCHES "throw$")
    expect(Error/Message "Test failure")
  else()
    expect(Error/Message "Test failure on rank 0")
  endif()
  if(failed STREQUAL "all")
    expect(Error/Collective ON)
    expect(Test/Ranks ${NP})
    expect(ElapsedTime/Counts/EPS 1)
  else()
    expect(Error/Collective OFF)
    expect(RootProcess/Test/Ranks 1)
    expect(RootProcess/ElapsedTime/Counts/EPS 1)
  endif()
endif()
