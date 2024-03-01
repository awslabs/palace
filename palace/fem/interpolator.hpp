// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTERPOLATOR_HPP
#define PALACE_FEM_INTERPOLATOR_HPP

#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class GridFunction;
class IoData;

//
// A class which wraps MFEM's GSLIB interface for high-order field interpolation.
//
class InterpolationOperator
{
private:
#if defined(MFEM_USE_GSLIB)
  mfem::FindPointsGSLIB op;
#endif
  std::vector<int> op_idx;

  std::vector<double> ProbeField(const mfem::ParGridFunction &U);

public:
  InterpolationOperator(const IoData &iodata, mfem::ParMesh &mesh);

  const auto &GetProbes() const { return op_idx; }

  std::vector<std::complex<double>> ProbeField(const GridFunction &U);
};

}  // namespace palace

#endif  // PALACE_FEM_INTERPOLATOR_HPP
