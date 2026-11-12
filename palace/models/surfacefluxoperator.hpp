// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP

#include <algorithm>
#include <map>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class IoData;
class Mesh;
class FiniteElementSpace;
template <ProblemType T>
class PostOperator;

namespace config
{
struct FluxLoopData;
}  // namespace config

class SurfaceFluxData
{
public:
  std::vector<int> hole_attributes;
  std::vector<double> flux_amounts;
  std::vector<int> film_attributes;
  std::vector<double> direction;
  // Excitation this loop is driven in (0: held at zero fluxoid in every excitation).
  int excitation;

public:
  SurfaceFluxData(const config::FluxLoopData &data);
  double GetExcitationFlux() const;
};

class SurfaceFluxOperator
{
private:
  std::map<int, SurfaceFluxData> sources;

  // Excitation index -> the flux loops driven together in it, each at its own fluxoid.
  // Each loop is its own excitation unless the configuration groups them, in which case
  // the solve gives the superposed state of the group.
  std::map<int, std::vector<int>> excitations;

  void SetUpBoundaryProperties(const IoData &iodata);
  void PrintBoundaryInfo(const IoData &iodata);

public:
  SurfaceFluxOperator() = default;
  SurfaceFluxOperator(const IoData &iodata);

  const SurfaceFluxData &GetSource(int idx) const;
  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto Size() const { return sources.size(); }

  const auto &GetExcitations() const { return excitations; }

  // Whether every excitation drives a single loop, so the solves span an inductance matrix
  // (the excitation and loop indices then agree, see config::Validate).
  bool HasSingleLoopExcitations() const
  {
    return std::ranges::all_of(excitations,
                               [](const auto &x) { return x.second.size() == 1; });
  }

  // Solve surface curl problem for given flux loop index
  palace::Vector
  SolveSurfaceCurlProblem(int idx, const Mesh &mesh, const FiniteElementSpace &nd_fespace,
                          PostOperator<ProblemType::MAGNETOSTATIC> &post_op) const;
  void SolveSurfaceCurlProblem(int idx, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                               palace::Vector &result) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_FLUX_OPERATOR_HPP
