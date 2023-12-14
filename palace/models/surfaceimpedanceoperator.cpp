// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceimpedanceoperator.hpp"

#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

SurfaceImpedanceOperator::SurfaceImpedanceOperator(const IoData &iodata,
                                                   const mfem::ParMesh &mesh)
{
  // Set up impedance boundary conditions.
  SetUpBoundaryProperties(iodata, mesh);
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
  Z_Rsinv.SetSize(bdr_attr_max);
  Z_Lsinv.SetSize(bdr_attr_max);
  Z_Cs.SetSize(bdr_attr_max);
  Z_Rsinv = 0.0;
  Z_Lsinv = 0.0;
  Z_Cs = 0.0;
  for (const auto &data : iodata.boundaries.impedance)
  {
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          Z_Rsinv(attr - 1) == 0.0 && Z_Lsinv(attr - 1) == 0.0 && Z_Cs(attr - 1) == 0.0,
          "Multiple definitions of impedance boundary properties for boundary attribute "
              << attr << "!");
      Z_Rsinv(attr - 1) = (std::abs(data.Rs) > 0.0) ? 1.0 / data.Rs : 0.0;
      Z_Lsinv(attr - 1) = (std::abs(data.Ls) > 0.0) ? 1.0 / data.Ls : 0.0;
      Z_Cs(attr - 1) = (std::abs(data.Cs) > 0.0) ? data.Cs : 0.0;
      MFEM_VERIFY(std::abs(Z_Rsinv(attr - 1)) + std::abs(Z_Lsinv(attr - 1)) +
                          std::abs(Z_Cs(attr - 1)) >
                      0.0,
                  "Impedance boundary has no Rs, Ls, or Cs defined!");
    }
  }

  // Mark selected boundary attributes from the mesh as impedance.
  mfem::Array<int> impedance_bcs, impedance_Rs_bcs, impedance_Ls_bcs, impedance_Cs_bcs;
  for (const auto &data : iodata.boundaries.impedance)
  {
    for (auto attr : data.attributes)
    {
      impedance_bcs.Append(attr);
      if (std::abs(Z_Rsinv(attr - 1)) > 0.0)
      {
        impedance_Rs_bcs.Append(attr);
      }
      if (std::abs(Z_Lsinv(attr - 1)) > 0.0)
      {
        impedance_Ls_bcs.Append(attr);
      }
      if (std::abs(Z_Cs(attr - 1)) > 0.0)
      {
        impedance_Cs_bcs.Append(attr);
      }
    }
  }
  mesh::AttrToMarker(bdr_attr_max, impedance_bcs, impedance_marker);
  mesh::AttrToMarker(bdr_attr_max, impedance_Rs_bcs, impedance_Rs_marker);
  mesh::AttrToMarker(bdr_attr_max, impedance_Ls_bcs, impedance_Ls_marker);
  mesh::AttrToMarker(bdr_attr_max, impedance_Cs_bcs, impedance_Cs_marker);
}

void SurfaceImpedanceOperator::PrintBoundaryInfo(const IoData &iodata,
                                                 const mfem::ParMesh &mesh)
{
  if (impedance_marker.Size() && impedance_marker.Max() == 0)
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin impedance BC at attributes:\n");
  for (int i = 0; i < impedance_marker.Size(); i++)
  {
    if (impedance_marker[i])
    {
      const int attr = i + 1;
      mfem::Vector normal = mesh::GetSurfaceNormal(mesh, attr);
      bool comma = false;
      Mpi::Print(" {:d}:", attr);
      if (std::abs(Z_Rsinv(i)) > 0.0)
      {
        Mpi::Print(
            " Rs = {:.3e} Ω/sq",
            iodata.DimensionalizeValue(IoData::ValueType::IMPEDANCE, 1.0 / Z_Rsinv(i)));
        comma = true;
      }
      if (std::abs(Z_Lsinv(i)) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(
            " Ls = {:.3e} H/sq",
            iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, 1.0 / Z_Lsinv(i)));
        comma = true;
      }
      if (std::abs(Z_Cs(i)) > 0.0)
      {
        if (comma)
        {
          Mpi::Print(",");
        }
        Mpi::Print(" Cs = {:.3e} F/sq",
                   iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, Z_Cs(i)));
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

void SurfaceImpedanceOperator::AddStiffnessBdrCoefficients(double coef,
                                                           SumMatrixCoefficient &fb)
{
  // Lumped inductor boundaries.
  if (impedance_Ls_marker.Size() && impedance_Ls_marker.Max() > 0)
  {
    mfem::Vector v(Z_Lsinv);
    v *= coef;
    auto f = std::make_unique<mfem::PWConstCoefficient>(v);
    fb.AddCoefficient(std::make_unique<mfem::PWConstCoefficient>(v), impedance_Ls_marker);
  }
}

void SurfaceImpedanceOperator::AddMassBdrCoefficients(double coef, SumMatrixCoefficient &fb)
{
  // Lumped capacitor boundaries.
  if (impedance_Cs_marker.Size() && impedance_Cs_marker.Max() > 0)
  {
    mfem::Vector v(Z_Cs);
    v *= coef;
    fb.AddCoefficient(std::make_unique<mfem::PWConstCoefficient>(v), impedance_Cs_marker);
  }
}

void SurfaceImpedanceOperator::AddDampingBdrCoefficients(double coef,
                                                         SumMatrixCoefficient &fb)
{
  // Lumped resistor boundaries.
  if (impedance_Rs_marker.Size() && impedance_Rs_marker.Max() > 0)
  {
    mfem::Vector v(Z_Rsinv);
    v *= coef;
    fb.AddCoefficient(std::make_unique<mfem::PWConstCoefficient>(v), impedance_Rs_marker);
  }
}

}  // namespace palace
