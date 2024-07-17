// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_COEFFICIENT_HPP
#define PALACE_LIBCEED_COEFFICIENT_HPP

#include <vector>
#include <mfem.hpp>

union CeedIntScalar;

namespace palace
{

class MaterialPropertyCoefficient;

namespace ceed
{

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim, const MaterialPropertyCoefficient *Q, double a = 1.0);

std::vector<CeedIntScalar>
PopulateCoefficientContext(int dim_mass, const MaterialPropertyCoefficient *Q_mass, int dim,
                           const MaterialPropertyCoefficient *Q, double a_mass = 1.0,
                           double a = 1.0);

struct QuadratureCoefficient
{
  int ncomp;
  mfem::Vector data;
  QuadratureCoefficient(int ncomp, int size) : ncomp(ncomp), data(size) {}
};

QuadratureCoefficient InitCoefficient(mfem::Coefficient &Q, mfem::ParMesh &mesh,
                                      const mfem::IntegrationRule &ir,
                                      const std::vector<int> &indices, bool use_bdr);

QuadratureCoefficient InitCoefficient(mfem::VectorCoefficient &VQ, mfem::ParMesh &mesh,
                                      const mfem::IntegrationRule &ir,
                                      const std::vector<int> &indices, bool use_bdr);

QuadratureCoefficient InitCoefficient(mfem::MatrixCoefficient &MQ, mfem::ParMesh &mesh,
                                      const mfem::IntegrationRule &ir,
                                      const std::vector<int> &indices, bool use_bdr);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_COEFFICIENT_HPP
