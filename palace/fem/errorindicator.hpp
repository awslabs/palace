// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_ERROR_INDICATORS_HPP
#define PALACE_FEM_ERROR_INDICATORS_HPP

#include <array>
#include <vector>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

//
// Storage for error estimation results from a simulation which involves one or more solves,
// required in the AMR loop.
//
class ErrorIndicator
{
protected:
  // Elemental localized error indicators. Used for marking elements for
  // refinement and coarsening.
  Vector local;

  // Number of samples.
  int n;

public:
  ErrorIndicator(Vector &&local) : local(std::move(local)), n(1)
  {
    this->local.UseDevice(true);
  }
  ErrorIndicator() : n(0) { local.UseDevice(true); }

  // Add an indicator to the running total.
  void AddIndicator(const Vector &indicator);

  // Return the local error indicator.
  const auto &Local() const { return local; }

  // Return the global error indicator.
  auto Norml2(MPI_Comm comm) const { return linalg::Norml2(comm, local); }

  // Return the largest local error indicator.
  auto Max(MPI_Comm comm) const
  {
    auto max = local.Max();
    Mpi::GlobalMax(1, &max, comm);
    return max;
  }

  // Return the smallest local error indicator.
  auto Min(MPI_Comm comm) const
  {
    auto min = local.Min();
    Mpi::GlobalMin(1, &min, comm);
    return min;
  }

  // Return the mean local error indicator.
  auto Mean(MPI_Comm comm) const { return linalg::Mean(comm, local); }

  struct SummaryStatistics
  {
    double norm;
    double min;
    double max;
    double mean;
  };

  SummaryStatistics GetSummaryStatistics(MPI_Comm comm) const
  {
    return {Norml2(comm), Min(comm), Max(comm), Mean(comm)};
  }
};

}  // namespace palace

#endif  // PALACE_FEM_ERROR_INDICATORS_HPP
