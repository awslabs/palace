// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceimpedanceoperator.hpp"

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceImpedanceOperator::SurfaceImpedanceOperator(const IoData &iodata,
                                                   const MaterialOperator &mat_op,
                                                   mfem::ParMesh &mesh)
  : mat_op(mat_op), impedance_attr(SetUpBoundaryProperties(iodata, mesh))
{
  // Print out BC info for all impedance boundary attributes.
  PrintBoundaryInfo(iodata, mesh);
}

void SurfaceImpedanceOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                       const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!iodata.boundaries.impedance.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
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
      }
    }
  }

  // Impedance boundaries are defined using the user provided impedance per square.
  mfem::Array<bool> marker(bdr_attr_max);
  marker = false;
  impedance_data.reserve(iodata.boundaries.impedance.size());
  for (const auto &data : iodata.boundaries.impedance)
  {
    MFEM_VERIFY(std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0,
                "Impedance boundary has no Rs, Ls, or Cs defined!");
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          !marker[attr - 1],
          "Multiple definitions of impedance boundary properties for boundary attribute "
              << attr << "!");
      marker[attr - 1] = true;
    }
    auto &Z = impedance_data.emplace_back();
    Z.Rsinv = (std::abs(data.Rs) > 0.0) ? 1.0 / data.Rs : 0.0;
    Z.Lsinv = (std::abs(data.Ls) > 0.0) ? 1.0 / data.Ls : 0.0;
    Z.Cs = (std::abs(data.Cs) > 0.0) ? data.Cs : 0.0;
    Z.attr.Append(data.attributes.data(), data.attributes.size());
  }

  // Mark selected boundary attributes from the mesh as impedance.
  mfem::Array<int> impedance_bcs;
  for (const auto &data : iodata.boundaries.impedance)
  {
    impedance_bcs.Append(data.attributes.data(), data.attributes.size());
    if (std::abs(data.Rs) > 0.0)
    {
      impedance_Rs_attr.Append(data.attributes.data(), data.attributes.size());
    }
    if (std::abs(data.Ls) > 0.0)
    {
      impedance_Ls_attr.Append(data.attributes.data(), data.attributes.size());
    }
    if (std::abs(data.Cs) > 0.0)
    {
      impedance_Cs_attr.Append(data.attributes.data(), data.attributes.size());
    }
  }
  return impedance_bcs;
}

void SurfaceImpedanceOperator::PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh)
{

  // XX TODO MARKER...

  if (impedance_data.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin impedance BC at attributes:\n");
  for (const auto &Z : impedance_data)
  {
    for (auto attr : Z.attr)
    {
      mfem::Vector nor;
      mesh::GetSurfaceNormal(mesh, attr, nor);
      bool comma = false;
      Mpi::Print(" {:d}:", attr);
      if (std::abs(Z.Rsinv(i)) > 0.0)
      {
        Mpi::Print(
            " Rs = {:.3e} Ω/sq",
            iodata.DimensionalizeValue(IoData::ValueType::IMPEDANCE, 1.0 / Z.Rsinv(i)));
        comma = true;
      }
      if (std::abs(Z.Lsinv(i)) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(
            " Ls = {:.3e} H/sq",
            iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, 1.0 / Z.Lsinv(i)));
        comma = true;
      }
      if (std::abs(Z.Cs(i)) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(" Cs = {:.3e} F/sq",
                   iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, Z.Cs(i)));
        comma = true;
      }
      if (comma)
      {
        Mpi::Print(",");
      }
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(" n = ({:+.1f}, {:+.1f}, {:+.1f})", nor(0), nor(1), nor(2));
      }
      else
      {
        Mpi::Print(" n = ({:+.1f}, {:+.1f})", nor(0), nor(1));
      }
      Mpi::Print("\n");
    }
  }
}

void SurfaceImpedanceOperator::AddStiffnessBdrCoefficients(double coef,
                                                           MaterialPropertyCoefficient &fb)
{
  // Lumped inductor boundaries.
  for (const auto &Z : impedance_data)
  {
    if (std::abs(Z.Lsinv) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetAttributeGlobalToLocal(Z.attr), Z.Lsinv, coef);
    }
  }
}

void SurfaceImpedanceOperator::AddDampingBdrCoefficients(double coef,
                                                         MaterialPropertyCoefficient &fb)
{
  // Lumped resistor boundaries.
  for (const auto &Z : impedance_data)
  {
    if (std::abs(Z.Rsinv) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetAttributeGlobalToLocal(Z.attr), Z.Rsinv, coef);
    }
  }
}

void SurfaceImpedanceOperator::AddMassBdrCoefficients(double coef,
                                                      MaterialPropertyCoefficient &fb)
{
  // Lumped capacitor boundaries.
  for (const auto &Z : impedance_data)
  {
    if (std::abs(Z.Cs) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetAttributeGlobalToLocal(Z.attr), Z.Cs, coef);
    }
  }
}

}  // namespace palace
