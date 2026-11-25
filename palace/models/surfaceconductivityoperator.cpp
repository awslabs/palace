// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceconductivityoperator.hpp"

#include <algorithm>
#include <set>
#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/permittivitymodel.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

using namespace std::complex_literals;

SurfaceConductivityOperator::SurfaceConductivityOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         const mfem::ParMesh &mesh)
  : SurfaceConductivityOperator(iodata.boundaries.conductivity, iodata.problem.type,
                                iodata.units, mat_op, mesh)
{
}

SurfaceConductivityOperator::SurfaceConductivityOperator(
    const std::vector<config::ConductivityData> &conductivity, ProblemType problem_type,
    const Units &units, const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(conductivity, problem_type, mesh);
  PrintBoundaryInfo(units, mesh);
}

void SurfaceConductivityOperator::SetUpBoundaryProperties(
    const std::vector<config::ConductivityData> &conductivity, ProblemType problem_type,
    const mfem::ParMesh &mesh)
{
  // Check that conductivity boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!conductivity.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, conductivity, "conductivity");
  }

  // Conductivity boundaries on the interior of the domain (not cracked) are conducting
  // sheets with a surface on either side. Attributes with both interior and exterior
  // boundary elements are treated as exterior.
  mfem::Array<int> int_marker, any_int_marker;
  if (!conductivity.empty() && mesh_attrs.Max() > 0)
  {
    int_marker.SetSize(mesh_attrs.Max());
    any_int_marker.SetSize(mesh_attrs.Max());
    int_marker = 1;
    any_int_marker = 1;
    mesh.UnmarkInternalBoundaries(int_marker, true);
    mesh.UnmarkInternalBoundaries(any_int_marker, false);
  }
  std::set<int> mixed_warn_list, external_warn_list;

  // Finite conductivity boundaries are defined using the user provided surface conductivity
  // and optionally conductor thickness.
  boundaries.reserve(conductivity.size());
  for (const auto &data : conductivity)
  {
    MFEM_VERIFY((data.sigma > 0.0 || data.epsilon_r > 0.0) && data.mu_r > 0.0,
                "Conductivity boundary has no conductivity or permittivity, or no "
                "permeability defined!");
    MFEM_VERIFY(data.h >= 0.0, "Conductivity boundary should have non-negative thickness!");
    mfem::Array<int> attr_list, ext_attr_list, int_attr_list;
    for (auto attr : data.attributes)
    {
      if (!mesh_attrs.Contains(attr))
      {
        continue;  // Can just ignore if wrong
      }
      attr_list.Append(attr);
      if (!int_marker[attr - 1])
      {
        int_attr_list.Append(attr);
        if (data.external)
        {
          external_warn_list.insert(attr);
        }
      }
      else
      {
        ext_attr_list.Append(attr);
        if (!any_int_marker[attr - 1])
        {
          mixed_warn_list.insert(attr);
        }
      }
    }

    // External surfaces have twice the effective thickness since the BC is applied at one
    // side. This does not apply to the interior boundaries (conducting sheets with two
    // surfaces), so with a finite thickness, these are separated into a group of their own,
    // with a single surface impedance per group.
    auto AddGroup = [&](const mfem::Array<int> &attrs, const mfem::Array<int> &ext_attrs,
                        const mfem::Array<int> &int_attrs, double h)
    {
      auto &bdr = boundaries.emplace_back();
      bdr.sigma = data.sigma;
      bdr.eps = data.epsilon_r;
      bdr.terms = data.permittivity_terms;
      bdr.mu = data.mu_r;
      bdr.h = h;
      bdr.attr_list = attrs;
      bdr.ext_attr_list = ext_attrs;
      bdr.int_attr_list = int_attrs;
    };
    const mfem::Array<int> none;
    if (data.external && data.h > 0.0 && ext_attr_list.Size() > 0 &&
        int_attr_list.Size() > 0)
    {
      AddGroup(ext_attr_list, ext_attr_list, none, 2.0 * data.h);
      AddGroup(int_attr_list, none, int_attr_list, data.h);
    }
    else
    {
      AddGroup(attr_list, ext_attr_list, int_attr_list,
               (data.external && ext_attr_list.Size() > 0) ? 2.0 * data.h : data.h);
    }
  }
  if (!external_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning("Conductivity boundary attributes with \"External\": true are interior "
                 "boundaries (conducting sheets with two surfaces), the thickness "
                 "correction for exterior surfaces does not apply!");
    utils::PrettyPrint(external_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  if (!mixed_warn_list.empty())
  {
    Mpi::Print("\n");
    Mpi::Warning(
        "Conductivity boundary attributes with both interior and exterior boundary "
        "elements are treated as exterior boundaries (one-sided), consider "
        "separating them into different attributes!");
    utils::PrettyPrint(mixed_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  MFEM_VERIFY(boundaries.empty() || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE ||
                  problem_type == ProblemType::BOUNDARYMODE,
              "Finite conductivity boundaries are only available for frequency "
              "domain simulations!");
}

void SurfaceConductivityOperator::PrintBoundaryInfo(const Units &units,
                                                    const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin finite conductivity BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      Mpi::Print(" {:d}: σ = {:.3e} S/m", attr,
                 units.Dimensionalize<Units::ValueType::CONDUCTIVITY>(bdr.sigma));
      if (bdr.eps > 0.0)
      {
        Mpi::Print(", ε = {:.3e}{}", bdr.eps,
                   bdr.terms.empty() ? "" : " (frequency-dependent)");
      }
      if (bdr.h > 0.0)
      {
        Mpi::Print(", h = {:.3e} m", units.Dimensionalize<Units::ValueType::LENGTH>(bdr.h));
      }
      Mpi::Print(", n = ({:+.1f})", fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
      const bool two_sided = std::find(bdr.int_attr_list.begin(), bdr.int_attr_list.end(),
                                       attr) != bdr.int_attr_list.end();
      Mpi::Print("{}\n", two_sided ? " (interior, two-sided)" : "");
    }
  }
}

mfem::Array<int> SurfaceConductivityOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

bool SurfaceConductivityOperator::HasRealFrequencyTerms() const
{
  return std::ranges::any_of(
      boundaries, [](const auto &bdr)
      { return std::ranges::any_of(bdr.terms, permittivity::IsRealFrequencyTerm); });
}

std::complex<double>
SurfaceConductivityOperator::EvaluateScalar(std::size_t group_idx,
                                            std::complex<double> omega) const
{
  // Surface impedance Z = iωμ/k of the conductor, k² = iωμσ(ω) with Re k > 0: (1 + i)/(σδ)
  // for a constant conductivity with skin depth δ, and Z coth(kh/2) for a thickness h
  // (Ansys HFSS user manual, "Surface Impedance Boundary Condition for Metal Traces of
  // Finite Thickness"). σ(ω) includes the displacement and polarization currents
  // iω(ε∞ + χ(iω)) of a given permittivity. The BC coefficient iω/Z = (k/μ) tanh(kh/2) is
  // analytic in ω, so a complex ω continues its real-ω value.
  MFEM_ASSERT(group_idx < boundaries.size(),
              "SurfaceConductivityOperator::EvaluateScalar: group_idx out of range!");
  const auto &bdr = boundaries[group_idx];
  const std::complex<double> s = 1i * omega;
  std::complex<double> sigma = bdr.sigma;
  if (bdr.eps > 0.0)
  {
    sigma += s * (bdr.eps + permittivity::Susceptibility(bdr.terms, s));
  }
  const std::complex<double> k = std::sqrt(s * bdr.mu * sigma);
  return (bdr.h > 0.0) ? k / bdr.mu * std::tanh(0.5 * k * bdr.h) : k / bdr.mu;
}

void SurfaceConductivityOperator::AddBoundaryMassBdrCoefficients(
    std::size_t g, MaterialPropertyCoefficient &fb, double coeff) const
{
  // ω-independent unit-coefficient boundary mass for group g (the operator that the scalar
  // i·ω/Z_g(ω) multiplies), times coeff. Interior conducting sheets have two surfaces, each
  // with surface admittance 1/Z. No contribution for an inactive group.
  if (IsActive(g))
  {
    const auto &bdr = boundaries[g];
    fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(bdr.ext_attr_list), coeff);
    fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(bdr.int_attr_list), 2.0 * coeff);
  }
}

}  // namespace palace
