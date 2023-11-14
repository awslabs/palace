// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_DORFLER_HPP
#define PALACE_UTILS_DORFLER_HPP

#include <array>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace mfem
{

class ParMesh;

}  // namespace mfem

namespace palace::utils
{

// Given a vector of estimates, e, and a fraction, compute a partition value E, such that
// that the set of all estimates with value greater than E, K_E, is the smallest set to
// achieve sum_{K_E} e² >= fraction * sum e². Namely the smallest set of elements that
// will mark the top fraction of the sum of the squared error. Returns as the second element
// in the pair the actual fraction of the total error.
// Reference: Dörfler, A convergent adaptive algorithm for Poisson’s equation, SIAM J.
//            Numer. Anal. (1996).
std::array<double, 2> ComputeDorflerThreshold(MPI_Comm comm, const Vector &e,
                                              double fraction);

// Given a nonconforming mesh, target fraction and error estimates, compute a threshold
// value and actual fraction that will mark the largest number of elements that make up the
// specified fraction of the total coarsening opportunities. This is analogous to
// ComputeDorflerThreshold, but operates only the list of available derefinement
// opportunities within the mesh.
std::array<double, 2> ComputeDorflerCoarseningThreshold(const mfem::ParMesh &mesh,
                                                        const Vector &e, double fraction);

}  // namespace palace::utils

#endif  // PALACE_UTILS_DORFLER_HPP
