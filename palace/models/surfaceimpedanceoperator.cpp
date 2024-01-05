// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceimpedanceoperator.hpp"

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceImpedanceOperator::SurfaceImpedanceOperator(const IoData &iodata,
                                                   const MaterialOperator &mat_op,
                                                   const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  // Print out BC info for all impedance boundary attributes.
  SetUpBoundaryProperties(iodata, mesh);
  PrintBoundaryInfo(iodata, mesh);
}

void SurfaceImpedanceOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                       const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  if (!iodata.boundaries.impedance.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), impedance_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    impedance_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &data : iodata.boundaries.impedance)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Impedance boundary attribute tags must be non-negative and correspond "
                    "to attributes in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown impedance boundary attribute " << attr << "!");
        MFEM_VERIFY(
            !impedance_marker[attr - 1],
            "Multiple definitions of impedance boundary properties for boundary attribute "
                << attr << "!");
        impedance_marker[attr - 1] = 1;
      }
    }
  }

  // Impedance boundaries are defined using the user provided impedance per square.
  boundaries.reserve(iodata.boundaries.impedance.size());
  for (const auto &data : iodata.boundaries.impedance)
  {
    MFEM_VERIFY(std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0,
                "Impedance boundary has no Rs, Ls, or Cs defined!");
    auto &bdr = boundaries.emplace_back();
    bdr.Rs = data.Rs;
    bdr.Ls = data.Ls;
    bdr.Cs = data.Cs;
    bdr.attr_list.Append(data.attributes.data(), data.attributes.size());
  }
}

void SurfaceImpedanceOperator::PrintBoundaryInfo(const IoData &iodata,
                                                 const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin impedance BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      mfem::Vector normal = mesh::GetSurfaceNormal(mesh, attr);
      bool comma = false;
      Mpi::Print(" {:d}:", attr);
      if (std::abs(bdr.Rs) > 0.0)
      {
        Mpi::Print(" Rs = {:.3e} Ω/sq",
                   iodata.DimensionalizeValue(IoData::ValueType::IMPEDANCE, bdr.Rs));
        comma = true;
      }
      if (std::abs(bdr.Ls) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(" Ls = {:.3e} H/sq",
                   iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, bdr.Ls));
        comma = true;
      }
      if (std::abs(bdr.Cs) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(" Cs = {:.3e} F/sq",
                   iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, bdr.Cs));
        comma = true;
      }
      if (comma)
      {
        Mpi::Print(",");
      }
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(" n = ({:+.1f}, {:+.1f}, {:+.1f})", normal(0), normal(1), normal(2));
      }
      else
      {
        Mpi::Print(" n = ({:+.1f}, {:+.1f})", normal(0), normal(1));
      }
      Mpi::Print("\n");
    }
  }
}

mfem::Array<int> SurfaceImpedanceOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

mfem::Array<int> SurfaceImpedanceOperator::GetRsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Rs) > 0.0)
    {
      attr_list.Append(bdr.attr_list);
    }
  }
  return attr_list;
}

mfem::Array<int> SurfaceImpedanceOperator::GetLsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Ls) > 0.0)
    {
      attr_list.Append(bdr.attr_list);
    }
  }
  return attr_list;
}

mfem::Array<int> SurfaceImpedanceOperator::GetCsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Cs) > 0.0)
    {
      attr_list.Append(bdr.attr_list);
    }
  }
  return attr_list;
}

void SurfaceImpedanceOperator::AddStiffnessBdrCoefficients(double coef,
                                                           MaterialPropertyCoefficient &fb)
{
  // Lumped inductor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Ls) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetBdrAttributeGlobalToLocal(bdr.attr_list),
                             coef / bdr.Ls);
    }
  }
}

void SurfaceImpedanceOperator::AddDampingBdrCoefficients(double coef,
                                                         MaterialPropertyCoefficient &fb)
{
  // Lumped resistor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Rs) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetBdrAttributeGlobalToLocal(bdr.attr_list),
                             coef / bdr.Rs);
    }
  }
}

void SurfaceImpedanceOperator::AddMassBdrCoefficients(double coef,
                                                      MaterialPropertyCoefficient &fb)
{
  // Lumped capacitor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Cs) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetBdrAttributeGlobalToLocal(bdr.attr_list),
                             coef * bdr.Cs);
    }
  }
}

}  // namespace palace
