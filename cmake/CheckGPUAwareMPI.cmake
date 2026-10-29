# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Check whether MPI can communicate CUDA or HIP device buffers, querying Open MPI with
# MPIX_Query_{cuda,rocm}_support and MPICH (4.0 or later) with MPIX_GPU_query_support. Other
# MPI implementations are taken as not GPU aware.
#

function(check_gpu_aware_mpi VAR)
  if(PALACE_WITH_CUDA)
    set(_gpu CUDA)
  else()
    set(_gpu HIP)
  endif()
  try_run(_run_result _compile_result
    SOURCE_FROM_CONTENT check_gpu_aware_mpi.c [=[
#include <stdio.h>
#include <mpi.h>
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif

int main(int argc, char *argv[])
{
  int aware = 0;
  MPI_Init(&argc, &argv);
#if defined(PALACE_GPU_CUDA) && defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  aware = MPIX_Query_cuda_support();
#elif defined(PALACE_GPU_HIP) && defined(MPIX_ROCM_AWARE_SUPPORT) && MPIX_ROCM_AWARE_SUPPORT
  aware = MPIX_Query_rocm_support();
#elif defined(PALACE_GPU_CUDA) && defined(MPICH_NUMVERSION) && MPICH_NUMVERSION >= 40000000
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA, &aware);
#elif defined(PALACE_GPU_HIP) && defined(MPICH_NUMVERSION) && MPICH_NUMVERSION >= 40100000
  MPIX_GPU_query_support(MPIX_GPU_SUPPORT_HIP, &aware);
#endif
  MPI_Finalize();
  printf("PALACE_GPU_AWARE_MPI=%d\n", aware);
  return 0;
}
]=]
    COMPILE_DEFINITIONS -DPALACE_GPU_${_gpu}
    LINK_LIBRARIES MPI::MPI_C
    RUN_OUTPUT_VARIABLE _output
  )
  if(_compile_result AND _run_result EQUAL 0 AND _output MATCHES "PALACE_GPU_AWARE_MPI=1")
    set(${VAR} ON PARENT_SCOPE)
    message(STATUS "MPI is ${_gpu} aware")
  else()
    set(${VAR} OFF PARENT_SCOPE)
    message(STATUS "MPI is not ${_gpu} aware")
  endif()
endfunction()
