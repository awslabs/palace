// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedelementoperator.hpp"

#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/units.hpp"

namespace palace
{

LumpedElementData::LumpedElementData(const config::LumpedElementData &data,
                                     const MaterialOperator &mat_op,
                                     const mfem::ParMesh &mesh)
  : mat_op(mat_op), topology(data.topology)
{
  bool has_circ = (std::abs(data.R) + std::abs(data.L) + std::abs(data.C) > 0.0);
  bool has_surf = (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0);
  MFEM_VERIFY(has_circ || has_surf,
              "Lumped element boundary has no R/L/C or Rs/Ls/Cs defined!");
  MFEM_VERIFY(!(has_circ && has_surf),
              "Lumped element boundary has both R/L/C and Rs/Ls/Cs defined!");

  // Single element only — no multielement support.
  MFEM_VERIFY(data.elements.size() == 1,
              "Lumped element boundary must have exactly one element!");
  const auto &e = data.elements[0];
  mfem::Array<int> attr_list;
  attr_list.Append(e.attributes.data(), e.attributes.size());
  switch (e.coordinate_system)
  {
    case CoordinateSystem::CYLINDRICAL:
      elem = std::make_unique<CoaxialLumpedGeometry>(e.direction, attr_list, mesh);
      break;
    case CoordinateSystem::CARTESIAN:
      elem = std::make_unique<UniformLumpedGeometry>(e.direction, attr_list, mesh);
      break;
  }

  if (has_circ)
  {
    R = data.R;
    L = data.L;
    C = data.C;
  }
  else
  {
    const double sq = elem->GetGeometryWidth() / elem->GetGeometryLength();
    if (topology == config::LumpedElementTopology::PARALLEL)
    {
      // Parallel: admittances add → Y_tot = Σ(sq/Rs)
      R = (std::abs(data.Rs) > 0.0) ? data.Rs / sq : 0.0;
      L = (std::abs(data.Ls) > 0.0) ? data.Ls / sq : 0.0;
      C = (std::abs(data.Cs) > 0.0) ? data.Cs * sq : 0.0;
    }
    else
    {
      // Series: impedances add → Z_tot = Σ(Rs*sq)
      R = (std::abs(data.Rs) > 0.0) ? data.Rs * sq : 0.0;
      L = (std::abs(data.Ls) > 0.0) ? data.Ls * sq : 0.0;
      C = (std::abs(data.Cs) > 0.0) ? data.Cs / sq : 0.0;
    }
  }
}

LumpedElementOperator::LumpedElementOperator(
    const std::map<int, config::LumpedElementData> &lumpedelement, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
{
  SetUpBoundaryProperties(lumpedelement, mat_op, mesh);
  PrintBoundaryInfo(units, mesh);
}

LumpedElementOperator::LumpedElementOperator(const IoData &iodata,
                                             const MaterialOperator &mat_op,
                                             const mfem::ParMesh &mesh)
  : LumpedElementOperator(iodata.boundaries.lumpedelement, iodata.units, mat_op, mesh)
{
}

void LumpedElementOperator::SetUpBoundaryProperties(
    const std::map<int, config::LumpedElementData> &lumpedelement,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
{
  if (lumpedelement.empty()) return;

  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker(bdr_attr_max), elem_marker(bdr_attr_max);
  bdr_attr_marker = 0;
  elem_marker = 0;
  for (auto attr : mesh.bdr_attributes)
  {
    bdr_attr_marker[attr - 1] = 1;
  }
  for (const auto &[idx, data] : lumpedelement)
  {
    for (const auto &elem : data.elements)
    {
      for (auto attr : elem.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Lumped element boundary attribute tags must be non-negative and "
                    "correspond to boundaries in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown lumped element boundary attribute " << attr << "!");
        MFEM_VERIFY(!elem_marker[attr - 1],
                    "Boundary attribute is assigned to more than one lumped element!");
        elem_marker[attr - 1] = 1;
      }
    }
  }

  for (const auto &[idx, data] : lumpedelement)
  {
    elements.try_emplace(idx, data, mat_op, mesh);
  }
}

void LumpedElementOperator::PrintBoundaryInfo(const Units &units,
                                              const mfem::ParMesh &mesh)
{
  using VT = Units::ValueType;
  if (elements.empty()) return;

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};

  fmt::format_to(out, "\nConfiguring Robin impedance BC for lumped elements at attributes:\n");
  for (const auto &[idx, data] : elements)
  {
    for (auto attr : data.elem->GetAttrList())
    {
      fmt::format_to(out, " {:d}:", attr);
      if (std::abs(data.R) > 0.0)
      {
        double Rs = data.R * data.GetToSquare();
        fmt::format_to(out, " Rs = {:.3e} Ω/sq,",
                       units.Dimensionalize<VT::IMPEDANCE>(Rs));
      }
      if (std::abs(data.L) > 0.0)
      {
        double Ls = data.L * data.GetToSquare();
        fmt::format_to(out, " Ls = {:.3e} H/sq,",
                       units.Dimensionalize<VT::INDUCTANCE>(Ls));
      }
      if (std::abs(data.C) > 0.0)
      {
        double Cs = data.C / data.GetToSquare();
        fmt::format_to(out, " Cs = {:.3e} F/sq,",
                       units.Dimensionalize<VT::CAPACITANCE>(Cs));
      }
      fmt::format_to(out, " n = ({:+.1f})\n",
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }

  fmt::format_to(out, "\nConfiguring lumped element circuit properties:\n");
  for (const auto &[idx, data] : elements)
  {
    std::vector<std::string> props;
    if (std::abs(data.R) > 0.0)
      props.emplace_back(fmt::format(
          "R = {:.3e} Ω", units.Dimensionalize<VT::IMPEDANCE>(data.R)));
    if (std::abs(data.L) > 0.0)
      props.emplace_back(fmt::format(
          "L = {:.3e} H", units.Dimensionalize<VT::INDUCTANCE>(data.L)));
    if (std::abs(data.C) > 0.0)
      props.emplace_back(fmt::format(
          "C = {:.3e} F", units.Dimensionalize<VT::CAPACITANCE>(data.C)));
    fmt::format_to(out, " Index = {:d}: {}\n", idx, fmt::join(props, ", "));
  }

  Mpi::Print("{}", fmt::to_string(buffer));
}

const LumpedElementData &LumpedElementOperator::GetElement(int idx) const
{
  auto it = elements.find(idx);
  MFEM_VERIFY(it != elements.end(), "Unknown lumped element index requested!");
  return it->second;
}

mfem::Array<int> LumpedElementOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : elements)
    attr_list.Append(data.elem->GetAttrList());
  return attr_list;
}

mfem::Array<int> LumpedElementOperator::GetRsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : elements)
    if (std::abs(data.R) > 0.0)
      attr_list.Append(data.elem->GetAttrList());
  return attr_list;
}

mfem::Array<int> LumpedElementOperator::GetLsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : elements)
    if (std::abs(data.L) > 0.0)
      attr_list.Append(data.elem->GetAttrList());
  return attr_list;
}

mfem::Array<int> LumpedElementOperator::GetCsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : elements)
    if (std::abs(data.C) > 0.0)
      attr_list.Append(data.elem->GetAttrList());
  return attr_list;
}

void LumpedElementOperator::AddStiffnessBdrCoefficients(double coeff,
                                                        MaterialPropertyCoefficient &fb)
{
  for (const auto &[idx, data] : elements)
    if (std::abs(data.L) > 0.0)
    {
      const double Ls = data.L * data.GetToSquare();
      fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(data.elem->GetAttrList()),
                             coeff / Ls);
    }
}

void LumpedElementOperator::AddDampingBdrCoefficients(double coeff,
                                                      MaterialPropertyCoefficient &fb)
{
  for (const auto &[idx, data] : elements)
    if (std::abs(data.R) > 0.0)
    {
      const double Rs = data.R * data.GetToSquare();
      fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(data.elem->GetAttrList()),
                             coeff / Rs);
    }
}

void LumpedElementOperator::AddMassBdrCoefficients(double coeff,
                                                   MaterialPropertyCoefficient &fb)
{
  for (const auto &[idx, data] : elements)
    if (std::abs(data.C) > 0.0)
    {
      const double Cs = data.C / data.GetToSquare();
      fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(data.elem->GetAttrList()),
                             coeff * Cs);
    }
}

}  // namespace palace
