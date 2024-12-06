// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurrentoperator.hpp"

#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceCurrentData::SurfaceCurrentData(const config::SurfaceCurrentData &data,
                                       const mfem::ParMesh &mesh)
{
  // Construct the source elements allowing for a possible multielement surface current
  // sources.
  for (const auto &elem : data.elements)
  {
    mfem::Array<int> attr_list;
    attr_list.Append(elem.attributes.data(), elem.attributes.size());
    switch (elem.coordinate_system)
    {
      case config::internal::ElementData::CoordinateSystem::CYLINDRICAL:
        elems.push_back(
            std::make_unique<CoaxialElementData>(elem.direction, attr_list, mesh));
        break;
      case config::internal::ElementData::CoordinateSystem::CARTESIAN:
        elems.push_back(
            std::make_unique<UniformElementData>(elem.direction, attr_list, mesh));
        break;
    }
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
  // Set up surface current source boundaries.
  SetUpBoundaryProperties(iodata, mesh);
  PrintBoundaryInfo(iodata, mesh);
}

void SurfaceCurrentOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                     const mfem::ParMesh &mesh)
{
  // Check that surface current boundary attributes have been specified correctly.
  if (!iodata.boundaries.current.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), source_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    source_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : iodata.boundaries.current)
    {
      for (const auto &elem : data.elements)
      {
        for (auto attr : elem.attributes)
        {
          MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                      "Surface current boundary attribute tags must be non-negative and "
                      "correspond to boundaries in the mesh!");
          MFEM_VERIFY(bdr_attr_marker[attr - 1],
                      "Unknown surface current boundary attribute " << attr << "!");
          MFEM_VERIFY(
              !source_marker[attr - 1],
              "Boundary attribute is assigned to more than one surface current source!");
          source_marker[attr - 1] = 1;
        }
      }
    }
  }

  // Set up surface current data structures.
  for (const auto &[idx, data] : iodata.boundaries.current)
  {
    sources.try_emplace(idx, data, mesh);
  }
}

void SurfaceCurrentOperator::PrintBoundaryInfo(const IoData &iodata,
                                               const mfem::ParMesh &mesh)
{
  if (sources.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring surface current excitation source term at attributes:\n");
  for (const auto &[idx, data] : sources)
  {
    for (const auto &elem : data.elems)
    {
      for (auto attr : elem->GetAttrList())
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
    for (const auto &elem : data.elems)
    {
      attr_list.Append(elem->GetAttrList());
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
  for (const auto &elem : data.elems)
  {
    const double Jinc = 1.0 / (elem->GetGeometryWidth() * data.elems.size());
    fb.AddCoefficient(elem->GetModeCoefficient(-Jinc));
  }
}

}  // namespace palace
