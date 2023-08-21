// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_ERROR_INDICATORS_HPP
#define PALACE_FEM_ERROR_INDICATORS_HPP

#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{
// Storage for error estimation results from the solve. Required in the AMR loop. This is
// richer than the IndicatorsAndNormalization because it stores derived quantities, and a
// communicator for use in adaptation.
class ErrorIndicator
{
public:
  ErrorIndicator(Vector &&local, double normalization)
    : local(std::move(local)), normalization(normalization), n(1)
  {
  }
  ErrorIndicator(const std::vector<ErrorIndicator> &ind)
  {
    for (const auto &i : ind)
    {
      AddIndicator(i);
    }
  }
  ErrorIndicator() = default;

  // Return the local error indicators.
  const auto &Local() const { return local; }
  // Return the global error indicator.
  inline auto Norml2(MPI_Comm comm) const { return linalg::Norml2(comm, local); }
  // Return the largest local error indicator.
  inline auto Max(MPI_Comm comm) const
  {
    double max = local.Max();
    Mpi::GlobalMax(1, &max, comm);
    return max;
  }
  // Return the smallest local error indicator.
  inline auto Min(MPI_Comm comm) const
  {
    double min = local.Min();
    Mpi::GlobalMin(1, &min, comm);
    return min;
  }
  // Return the mean local error indicator.
  inline auto Mean(MPI_Comm comm) const
  {
    int size = local.Size();
    Mpi::GlobalSum(1, &size, comm);
    auto sum = local.Sum();
    Mpi::GlobalSum(1, &sum, comm);
    return sum / size;
  }
  // Return the normalization constant for the absolute error.
  inline auto GetNormalization() const { return normalization; }
  // Add a set of indicators to the running totals.
  void AddIndicator(const ErrorIndicator &indicators);
  // Return a vector of postprocess data.
  std::array<double, 5> GetPostprocessData(MPI_Comm comm) const
  {
    return {Norml2(comm), Min(comm), Max(comm), Mean(comm), GetNormalization()};
  }

protected:
  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local;
  // Normalization constant.
  double normalization;
  // Number of samples. Mutability required to guarantee operation.
  int n = 0;
};

}  // namespace palace

#endif  // PALACE_FEM_ERROR_INDICATORS_HPP
