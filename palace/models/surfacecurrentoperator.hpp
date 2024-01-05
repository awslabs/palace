// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CURRENT_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_CURRENT_OPERATOR_HPP

#include <map>
#include <memory>
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
// Helper class for surface current boundaries in a model.
//
class SurfaceCurrentData
{
public:
  // To accomodate multielement surface current sources, a current source may be made up
  // of elements with different attributes and directions which add to deliver the same
  // total source current.
  std::vector<std::unique_ptr<LumpedElementData>> elems;

public:
  SurfaceCurrentData(const config::SurfaceCurrentData &data,
                     mfem::ParFiniteElementSpace &h1_fespace);

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

  void SetUpBoundaryProperties(const IoData &iodata,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  SurfaceCurrentOperator(const IoData &iodata, mfem::ParFiniteElementSpace &h1_fespace);

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
