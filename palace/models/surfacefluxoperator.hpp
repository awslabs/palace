// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class Mesh;
class FiniteElementSpace;

namespace config
{
struct FluxLoopData;
}  // namespace config

class SurfaceFluxData
{
public:
  std::vector<int> hole_attributes;
  std::vector<double> flux_amounts;
  std::vector<int> fluxloop_pec;
  std::vector<double> direction;
  double regularization;

public:
  SurfaceFluxData(const config::FluxLoopData &data);
  double GetExcitationFlux() const;
};

class SurfaceFluxOperator
{
private:
  std::map<int, SurfaceFluxData> sources;
  std::unique_ptr<IoData> solver_config_;  // Store only solver configuration

  void SetUpBoundaryProperties(const IoData &iodata);
  void PrintBoundaryInfo(const IoData &iodata);

public:
  SurfaceFluxOperator(const IoData &iodata);

  const SurfaceFluxData &GetSource(int idx) const;
  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto Size() const { return sources.size(); }

  // Solve surface curl problem for given flux loop index
  palace::Vector SolveSurfaceCurlProblem(int idx, const Mesh &mesh,
                                         const FiniteElementSpace &nd_fespace) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP