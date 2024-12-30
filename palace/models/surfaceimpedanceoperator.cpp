// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceimpedanceoperator.hpp"

#include <set>
#include <fmt/ranges.h>
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!iodata.boundaries.impedance.empty())
  {
    mfem::Array<int> impedance_marker(bdr_attr_max);
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    impedance_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (const auto &data : iodata.boundaries.impedance)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(
            !impedance_marker[attr - 1],
            "Multiple definitions of impedance boundary properties for boundary attribute "
                << attr << "!");
        impedance_marker[attr - 1] = 1;
        // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
        //             "Impedance boundary attribute tags must be non-negative and
        //             correspond " "to attributes in the mesh!");
        // MFEM_VERIFY(bdr_attr_marker[attr - 1],
        //             "Unknown impedance boundary attribute " << attr << "!");
        if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
        {
          bdr_warn_list.insert(attr);
        }
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown impedance boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
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
    bdr.attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
    }
  }
}

void SurfaceImpedanceOperator::PrintBoundaryInfo(const IoData &iodata,
                                                 const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buf{};  // Output buffer & buffer append lambda for cleaner code
  auto to = [&buf](auto fmt, auto &&...args)
  { fmt::format_to(std::back_inserter(buf), fmt, std::forward<decltype(args)>(args)...); };

  using VT = IoData::ValueType;

  to("\nConfiguring Robin impedance BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      to(" {:d}:", attr);
      if (std::abs(bdr.Rs) > 0.0)
      {
        to(" Rs = {:.3e} Ω/sq,", iodata.DimensionalizeValue(VT::IMPEDANCE, bdr.Rs));
      }
      if (std::abs(bdr.Ls) > 0.0)
      {
        to(" Ls = {:.3e} H/sq,", iodata.DimensionalizeValue(VT::INDUCTANCE, bdr.Ls));
      }
      if (std::abs(bdr.Cs) > 0.0)
      {
        to(" Cs = {:.3e} F/sq,", iodata.DimensionalizeValue(VT::CAPACITANCE, bdr.Cs));
      }
      to(" n = ({:+.1f})\n", fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buf));
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

void SurfaceImpedanceOperator::AddStiffnessBdrCoefficients(double coeff,
                                                           MaterialPropertyCoefficient &fb)
{
  // Lumped inductor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Ls) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(bdr.attr_list), coeff / bdr.Ls);
    }
  }
}

void SurfaceImpedanceOperator::AddDampingBdrCoefficients(double coeff,
                                                         MaterialPropertyCoefficient &fb)
{
  // Lumped resistor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Rs) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(bdr.attr_list), coeff / bdr.Rs);
    }
  }
}

void SurfaceImpedanceOperator::AddMassBdrCoefficients(double coeff,
                                                      MaterialPropertyCoefficient &fb)
{
  // Lumped capacitor boundaries.
  for (const auto &bdr : boundaries)
  {
    if (std::abs(bdr.Cs) > 0.0)
    {
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(bdr.attr_list), coeff * bdr.Cs);
    }
  }
}

}  // namespace palace
