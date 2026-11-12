# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Runs the hidden "Failure report on abort" test case, which aborts on purpose, and checks
# the palace.json it leaves. Arguments (-D): EXE, MPIEXEC, NP_FLAG, NP, PREFLAGS, POSTFLAGS,
# MODE (root-abort: the root process alone fails, through an MFEM error message and
# MPI_Abort; all-throw: all processes fail, through an uncaught exception) and DIR. Either
# way, the root writes its own report.

file(REMOVE_RECURSE "${DIR}")
file(MAKE_DIRECTORY "${DIR}")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E env PALACE_TEST_FAILURE=${MODE} PALACE_TEST_FAILURE_DIR=${DIR}
          ${MPIEXEC} ${NP_FLAG} ${NP} ${PREFLAGS} ${EXE} ${POSTFLAGS}
          "Failure report on abort"
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE out
  TIMEOUT 300
)
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

expect(Complete OFF)
expect(Error/Phase EPS)
if(MODE MATCHES "throw$")
  expect(Error/Message "Test failure")
else()
  expect(Error/Message "Test failure on rank 0")
endif()
expect(RootProcess/Test/Ranks 1)
expect(RootProcess/ElapsedTime/Counts/EPS 1)
