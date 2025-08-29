// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacefluxoperator.hpp"
#include <memory>
#include "models/surfacecurlsolver.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceFluxData::SurfaceFluxData(const config::FluxLoopData &data)
  : hole_attributes(data.hole_attributes), flux_amounts(data.flux_amounts),
    fluxloop_pec(data.fluxloop_pec),
    direction(data.direction.begin(), data.direction.end()),
    regularization(data.regularization)
{
}

double SurfaceFluxData::GetExcitationFlux() const
{
  return flux_amounts.empty() ? 0.0 : flux_amounts[0];
}

SurfaceFluxOperator::SurfaceFluxOperator(const IoData &iodata)
  : solver_config_(std::make_unique<IoData>(iodata))
{
  SetUpBoundaryProperties(iodata);
  PrintBoundaryInfo(iodata);
}

void SurfaceFluxOperator::SetUpBoundaryProperties(const IoData &iodata)
{
  for (const auto &[idx, data] : iodata.boundaries.fluxloop)
  {
    sources.try_emplace(idx, data);
  }
}

void SurfaceFluxOperator::PrintBoundaryInfo(const IoData &iodata)
{
  if (sources.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring flux loop excitation sources:\n");
  for (const auto &[idx, data] : sources)
  {
    Mpi::Print(" Index = {:d}, holes = [{}]\n", idx, fmt::join(data.hole_attributes, ","));
  }
}

const SurfaceFluxData &SurfaceFluxOperator::GetSource(int idx) const
{
  auto it = sources.find(idx);
  MFEM_VERIFY(it != sources.end(), "Unknown flux source index requested!");
  return it->second;
}

palace::Vector
SurfaceFluxOperator::SolveSurfaceCurlProblem(int idx, const Mesh &mesh,
                                             const FiniteElementSpace &nd_fespace) const
{
  const auto &data = GetSource(idx);
  auto result =
      palace::SolveSurfaceCurlProblem(data, *solver_config_, mesh, nd_fespace, idx);
  return *result;
}

}  // namespace palace