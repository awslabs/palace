// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superconductorsheetoperator.hpp"

#include <cmath>
#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

SuperconductorSheetOperator::SuperconductorSheetOperator(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(superconductor, cracked_attributes, mesh);
  PrintBoundaryInfo(units, mesh);
}

SuperconductorSheetOperator::SuperconductorSheetOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         const mfem::ParMesh &mesh)
  : SuperconductorSheetOperator(iodata.boundaries.superconductor,
                                iodata.boundaries.cracked_attributes, iodata.units, mat_op,
                                mesh)
{
}

void SuperconductorSheetOperator::SetUpBoundaryProperties(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const mfem::ParMesh &mesh)
{
  // Check that superconductor sheet boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!superconductor.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, superconductor, "superconductor sheet");
  }

  // Kinetic sheet inductance L_ksq [H/sq]: supplied directly, or from (lambda, d) via the
  // finite-thickness London sheet L_ksq = lambda*coth(d/lambda) (nondimensional, mu0
  // absorbed). Reduces to the thin-film Pearl limit lambda^2/d for d << lambda and
  // saturates at lambda for d >> lambda, capturing the through-thickness current profile
  // without resolving the film thickness.
  boundaries.reserve(superconductor.size());
  for (const auto &data : superconductor)
  {
    const double Ls =
        (data.Ls > 0.0) ? data.Ls : KineticSheetInductance(data.lambda_L, data.thickness);
    MFEM_VERIFY(Ls > 0.0,
                "Superconductor sheet has non-positive kinetic sheet inductance!");
    auto &bdr = boundaries.emplace_back();
    bdr.Ls = Ls;
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

void SuperconductorSheetOperator::PrintBoundaryInfo(const Units &units,
                                                    const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  fmt::format_to(out, "\nConfiguring superconductor sheet BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      fmt::format_to(out, " {:d}: L_ksq = {:.3e} H/sq, n = ({:+.1f})\n", attr,
                     units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buffer));
}

double SuperconductorSheetOperator::KineticSheetInductance(double lambda_L,
                                                           double thickness)
{
  // Finite-thickness London sheet inductance L_ksq = lambda * coth(d/lambda), with
  // coth(x) = 1/tanh(x). Reduces to lambda^2/d for d << lambda; saturates at lambda for
  // d >> lambda.
  return lambda_L / std::tanh(thickness / lambda_L);
}

mfem::Array<int> SuperconductorSheetOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

void SuperconductorSheetOperator::AddStiffnessBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fb) const
{
  // Kinetic sheet inductance boundaries: add (coeff / L_ksq) as a tangential surface mass.
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      const double s = bdr.attr_scaling.at(attr);
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Ls * s));
    }
  }
}

}  // namespace palace
