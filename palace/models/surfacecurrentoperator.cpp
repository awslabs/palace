// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurrentoperator.hpp"

#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceCurrentData::SurfaceCurrentData(const config::SurfaceCurrentData &data,
                                       const mfem::ParMesh &mesh)
{
  MFEM_VERIFY(!data.elements.empty(),
              "Surface current source must have at least one element!");
  const double current_fraction = 1.0 / static_cast<double>(data.elements.size());

  // Construct the parallel source elements together with their corresponding apertures.
  for (const auto &elem : data.elements)
  {
    mfem::Array<int> attr_list;
    attr_list.Append(elem.attributes.data(), static_cast<int>(elem.attributes.size()));
    std::unique_ptr<LumpedElementData> source;
    switch (elem.coordinate_system)
    {
      case CoordinateSystem::CYLINDRICAL:
        source = std::make_unique<CoaxialElementData>(elem.direction, attr_list, mesh);
        break;
      case CoordinateSystem::CARTESIAN:
        source = std::make_unique<UniformElementData>(elem.direction, attr_list, mesh);
        break;
    }

    std::optional<SurfaceCurrentAperture> aperture;
    if (elem.aperture)
    {
      aperture.emplace();
      aperture->attributes = elem.aperture->attributes;
      aperture->direction.SetSize(static_cast<int>(elem.aperture->direction.size()));
      for (int i = 0; i < aperture->direction.Size(); i++)
      {
        aperture->direction[i] = elem.aperture->direction[i];
      }
    }
    elements.push_back({std::move(source), std::move(aperture), current_fraction});
  }
}

double SurfaceCurrentData::GetExcitationCurrent() const
{
  // Ideal unit current source for each index.
  return 1.0;
}

SurfaceCurrentOperator::SurfaceCurrentOperator(const IoData &iodata,
                                               const mfem::ParMesh &mesh)
{
  SetUpBoundaryProperties(iodata.boundaries.current, mesh);
  PrintBoundaryInfo(mesh);
}

SurfaceCurrentOperator::SurfaceCurrentOperator(
    const std::map<int, config::SurfaceCurrentData> &current, const mfem::ParMesh &mesh)
{
  SetUpBoundaryProperties(current, mesh);
  PrintBoundaryInfo(mesh);
}

void SurfaceCurrentOperator::SetUpBoundaryProperties(
    const std::map<int, config::SurfaceCurrentData> &current, const mfem::ParMesh &mesh)
{
  // Check that surface current boundary attributes have been specified correctly.
  if (!current.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), source_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    source_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : current)
    {
      for (const auto &elem : data.elements)
      {
        for (auto attr : elem.attributes)
        {
          MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                      "Surface current boundary attribute tags must be positive and "
                      "correspond to boundaries in the mesh!");
          MFEM_VERIFY(bdr_attr_marker[attr - 1],
                      "Unknown surface current boundary attribute " << attr << "!");
          MFEM_VERIFY(
              !source_marker[attr - 1],
              "Boundary attribute is assigned to more than one surface current source!");
          source_marker[attr - 1] = 1;
        }
        if (elem.aperture)
        {
          MFEM_VERIFY(!elem.aperture->attributes.empty(),
                      "SurfaceCurrent Aperture must contain at least one boundary "
                      "attribute for SurfaceCurrent index "
                          << idx << "!");
          for (auto attr : elem.aperture->attributes)
          {
            MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max && bdr_attr_marker[attr - 1],
                        "Unknown SurfaceCurrent Aperture boundary attribute "
                            << attr << " for SurfaceCurrent index " << idx << "!");
          }
        }
      }
    }
  }

  // Set up surface current data structures.
  for (const auto &[idx, data] : current)
  {
    sources.try_emplace(idx, data, mesh);
  }
}

void SurfaceCurrentOperator::PrintBoundaryInfo(const mfem::ParMesh &mesh)
{
  if (sources.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring surface current excitation source term at attributes:\n");
  for (const auto &[idx, data] : sources)
  {
    for (const auto &elem : data.elements)
    {
      for (auto attr : elem.source->GetAttrList())
      {
        Mpi::Print(" {:d}: Index = {:d}, n = ({:+.1f})\n", attr, idx,
                   fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
      }
    }
  }
}

const SurfaceCurrentData &SurfaceCurrentOperator::GetSource(int idx) const
{
  auto it = sources.find(idx);
  MFEM_VERIFY(it != sources.end(), "Unknown current source index requested!");
  return it->second;
}

mfem::Array<int> SurfaceCurrentOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : sources)
  {
    for (const auto &elem : data.elements)
    {
      attr_list.Append(elem.source->GetAttrList());
    }
  }
  return attr_list;
}

void SurfaceCurrentOperator::AddExcitationBdrCoefficients(SumVectorCoefficient &fb)
{
  // Construct the RHS source term for surface current boundaries, which looks like
  // -iω J_inc for a surface current boundary. The chosen surface current J_inc corresponds
  // to a unit current excitation. Note: The real RHS returned here does not yet have the
  // factor of (iω) included, so works for time domain simulations requiring RHS -J_inc
  // (t).
  for (const auto &[idx, data] : sources)
  {
    AddExcitationBdrCoefficients(data, fb);
  }
}

void SurfaceCurrentOperator::AddExcitationBdrCoefficients(int idx, SumVectorCoefficient &fb)
{
  // Construct the RHS source term for a single surface current boundary index.
  AddExcitationBdrCoefficients(GetSource(idx), fb);
}

void SurfaceCurrentOperator::AddExcitationBdrCoefficients(const SurfaceCurrentData &data,
                                                          SumVectorCoefficient &fb)
{
  // Add excited boundaries to the linear form, with a unit current distributed across
  // all elements of the current source in parallel.
  for (const auto &elem : data.elements)
  {
    const double Jinc = elem.current_fraction / elem.source->GetGeometryWidth();
    fb.AddCoefficient(elem.source->GetModeCoefficient(-Jinc));
  }
}

}  // namespace palace
