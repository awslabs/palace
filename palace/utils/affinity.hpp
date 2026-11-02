// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_AFFINITY_HPP
#define PALACE_UTILS_AFFINITY_HPP

#include <mpi.h>

namespace palace::utils
{

// If the MPI launcher packed the ranks of a node onto its first cores (one core per rank,
// local rank i on core i, as OpenMPI's default core mapping does) and this leaves NUMA
// domains with fewer ranks than others, rebind each rank so the ranks of a node are spread
// evenly across the NUMA domains (and their memory bandwidth). Call before MPI
// initialization so that memory is first touched on the new domain. Disabled with
// PALACE_DISABLE_REBIND=1.
void SpreadPackedProcessBinding();

// After MPI initialization: undo the rebinding on a node if not all of its ranks agree or
// cores would be shared, and report it.
void CheckProcessBinding(MPI_Comm comm);

}  // namespace palace::utils

#endif  // PALACE_UTILS_AFFINITY_HPP
