// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <mfem.hpp>

namespace palace
{

// Parallel-capable interface identification in true-DOF space. Given the parent finite
// element space and the region/environment domain attribute sets, mark each parent true DOF
// touched by region element support, by environment support, and by both (the shared
// interface). Works on any MPI partition: L-vector attribute marks are reduced to true DOFs
// across ranks via the parallel prolongation. Returns the global interface true-DOF count.
int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker,
                          mfem::Array<int> &interface_marker);

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
