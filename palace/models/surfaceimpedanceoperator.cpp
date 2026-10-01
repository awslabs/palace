// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceimpedanceoperator.hpp"

#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

SurfaceImpedanceOperator::SurfaceImpedanceOperator(
    const std::vector<config::ImpedanceData> &impedance,
    const std::unordered_set<int> &cracked_attributes, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(impedance, cracked_attributes, mesh);
  PrintBoundaryInfo(units, mesh);
}

SurfaceImpedanceOperator::SurfaceImpedanceOperator(const IoData &iodata,
                                                   const MaterialOperator &mat_op,
                                                   const mfem::ParMesh &mesh)
  : SurfaceImpedanceOperator(iodata.boundaries.impedance,
                             iodata.boundaries.cracked_attributes, iodata.units, mat_op,
                             mesh)
{
}

void SurfaceImpedanceOperator::SetUpBoundaryProperties(
    const std::vector<config::ImpedanceData> &impedance,
    const std::unordered_set<int> &cracked_attributes, const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!impedance.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, impedance, "impedance");
  }

  // Impedance boundaries are defined using the user provided impedance per square.
  boundaries.reserve(impedance.size());
  for (const auto &data : impedance)
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
      if (!mesh_attrs.Contains(attr))
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
      // Per-attribute scaling to account for increased area when using mesh cracking.
      bdr.attr_scaling[attr] =
          (cracked_attributes.find(attr) != cracked_attributes.end()) ? 2.0 : 1.0;
    }
  }
}

void SurfaceImpedanceOperator::PrintBoundaryInfo(const Units &units,
                                                 const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  fmt::format_to(out, "\nConfiguring Robin impedance BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      fmt::format_to(out, " {:d}:", attr);
      if (std::abs(bdr.Rs) > 0.0)
      {
        fmt::format_to(out, " Rs = {:.3e} Ω/sq,",
                       units.Dimensionalize<VT::IMPEDANCE>(bdr.Rs));
      }
      if (std::abs(bdr.Ls) > 0.0)
      {
        fmt::format_to(out, " Ls = {:.3e} H/sq,",
                       units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls));
      }
      if (std::abs(bdr.Cs) > 0.0)
      {
        fmt::format_to(out, " Cs = {:.3e} F/sq,",
                       units.Dimensionalize<VT::CAPACITANCE>(bdr.Cs));
      }
      fmt::format_to(out, " n = ({:+.1f})\n",
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buffer));
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
      for (auto attr : bdr.attr_list)
      {
        const double s = bdr.attr_scaling.at(attr);
        fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Ls * s));
      }
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
      for (auto attr : bdr.attr_list)
      {
        const double s = bdr.attr_scaling.at(attr);
        fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Rs * s));
      }
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
      for (auto attr : bdr.attr_list)
      {
        const double s = bdr.attr_scaling.at(attr);
        fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff * bdr.Cs / s);
      }
    }
  }
}

}  // namespace palace
