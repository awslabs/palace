// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CURRENT_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CURRENT_OPERATOR_HPP

#include <map>
#include <memory>
#include <optional>
#include <vector>
#include <mfem.hpp>
#include "fem/lumpedelement.hpp"

namespace palace
{

class IoData;
class SumVectorCoefficient;

namespace config
{

struct SurfaceCurrentData;

}  // namespace config

//
// Helper classes for surface current boundaries in a model.
//
struct SurfaceCurrentAperture
{
  // Boundary attributes tiling the spanning surface, and a normalized Cartesian reference
  // direction defining its positive oriented normal.
  std::vector<int> attributes;
  mfem::Vector direction;
};

struct SurfaceCurrentElement
{
  std::unique_ptr<LumpedElementData> source;
  std::optional<SurfaceCurrentAperture> aperture;

  // Fraction of the total port current carried by this parallel element. This same weight
  // must be used when combining the element flux linkages into the generalized port flux.
  double current_fraction;
};

class SurfaceCurrentData
{
public:
  // A multielement source consists of parallel elements which add to deliver the total
  // source current. Each element owns its corresponding aperture when mixed current/flux
  // inductance extraction is configured.
  std::vector<SurfaceCurrentElement> elements;

  SurfaceCurrentData(const config::SurfaceCurrentData &data, const mfem::ParMesh &mesh);

  double GetExcitationCurrent() const;
};

//
// A class handling surface current boundaries.
//
class SurfaceCurrentOperator
{
private:
  // Mapping from source index to data structure containing source surface current
  // information.
  std::map<int, SurfaceCurrentData> sources;

  void SetUpBoundaryProperties(const std::map<int, config::SurfaceCurrentData> &current,
                               const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const mfem::ParMesh &mesh);

public:
  SurfaceCurrentOperator(const std::map<int, config::SurfaceCurrentData> &current,
                         const mfem::ParMesh &mesh);
  SurfaceCurrentOperator(const IoData &iodata, const mfem::ParMesh &mesh);

  // Access data structures for the surface current source with the given index.
  const SurfaceCurrentData &GetSource(int idx) const;
  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto rbegin() const { return sources.rbegin(); }
  auto rend() const { return sources.rend(); }
  auto Size() const { return sources.size(); }

  // Returns array of surface current source attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to the right-hand side source term vector for a surface current
  // excitation at the specified boundaries, -J_inc for the real version (versus the
  // full -iω J_inc for the complex one).
  void AddExcitationBdrCoefficients(SumVectorCoefficient &fb);
  void AddExcitationBdrCoefficients(int idx, SumVectorCoefficient &fb);
  void AddExcitationBdrCoefficients(const SurfaceCurrentData &data,
                                    SumVectorCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CURRENT_OPERATOR_HPP
