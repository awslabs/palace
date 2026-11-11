#!/bin/sh
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# CTest launcher for the MPI test cases: runs the command (`mpiexec -n N
# palace-unit-tests <case>`) with a private Open MPI shared-memory backing
# directory and removes it afterwards.
#
# Open MPI's btl/sm (btl/vader in 4.x) backs each rank with the file
# <dir>/sm_segment.<host>.<uid>.<jobid>.<rank>, <dir> = /dev/shm by default,
# and <jobid> keeps only a 16-bit hash of the launcher's PID. With dozens of
# concurrent mpiexec under `ctest -j`, two jobs occasionally get the same file
# names and share segments: their ranks crash, hang in the first collectives,
# or fail in MPI_Init with "open(2) ... No such file or directory".
base=${OMPI_MCA_btl_sm_backing_directory:-${OMPI_MCA_btl_vader_backing_directory:-/dev/shm}}
if [ -d "$base" ] && [ -w "$base" ] && dir=$(mktemp -d "$base/palace-test.XXXXXX"); then
  trap 'rm -rf "$dir"' EXIT
  trap 'exit 129' HUP
  trap 'exit 130' INT
  trap 'exit 143' TERM
  OMPI_MCA_btl_sm_backing_directory=$dir
  OMPI_MCA_btl_vader_backing_directory=$dir
  export OMPI_MCA_btl_sm_backing_directory OMPI_MCA_btl_vader_backing_directory
  "$@"
  exit $?
fi
exec "$@"
