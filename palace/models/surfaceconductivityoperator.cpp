// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceconductivityoperator.hpp"

#include <set>
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!conductivity.empty())
  {
    mfem::Array<int> conductivity_marker(bdr_attr_max);
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    conductivity_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (const auto &data : conductivity)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(!conductivity_marker[attr - 1],
                    "Multiple definitions of conductivity boundary properties for boundary "
                    "attribute "
                        << attr << "!");
        conductivity_marker[attr - 1] = 1;
        if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
        {
          bdr_warn_list.insert(attr);
        }
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning(
          "Unknown conductivity boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // Finite conductivity boundaries are defined using the user provided surface conductivity
  // and optionally conductor thickness.
  boundaries.reserve(conductivity.size());
  for (const auto &data : conductivity)
  {
    MFEM_VERIFY(data.sigma > 0.0 && data.mu_r > 0.0,
                "Conductivity boundary has no conductivity or no "
                "permeability defined!");
    MFEM_VERIFY(data.h >= 0.0, "Conductivity boundary should have non-negative thickness!");
    auto &bdr = boundaries.emplace_back();
    bdr.sigma = data.sigma;
    bdr.mu = data.mu_r;
    bdr.h = data.h;
    if (data.external)
    {
      // External surfaces have twice the effective thickness since the BC is applied at one
      // side.
      bdr.h *= 2.0;
    }
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
      if (bdr.h > 0.0)
      {
        Mpi::Print(", h = {:.3e} m", units.Dimensionalize<Units::ValueType::LENGTH>(bdr.h));
      }
      Mpi::Print(", n = ({:+.1f})\n", fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
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

void SurfaceConductivityOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  // Real-ω stamping: per-group complex coefficient i·ω/Z(ω) split into (fbr, fbi). The
  // formula and finite-thickness correction live in EvaluateScalar — see that method
  // for the closed-form skin-depth derivation and the analytic continuation that the
  // complex-λ path reuses.
  for (std::size_t g = 0; g < boundaries.size(); g++)
  {
    if (!(std::abs(boundaries[g].sigma) > 0.0))
    {
      continue;
    }
    const std::complex<double> s = EvaluateScalar(g, std::complex<double>{omega, 0.0});
    fbr.AddMaterialProperty(mat_op.GetCeedBdrAttributes(boundaries[g].attr_list), s.real());
    fbi.AddMaterialProperty(mat_op.GetCeedBdrAttributes(boundaries[g].attr_list), s.imag());
  }
}

void SurfaceConductivityOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi)
{
  // Complex-ω stamping: identical to the real-ω overload but evaluates the per-group
  // i·ω/Z(ω) at a genuinely complex ω (analytic continuation, ω = -i·λ). The skin-depth
  // sqrt and finite-thickness correction in EvaluateScalar already handle complex ω.
  for (std::size_t g = 0; g < boundaries.size(); g++)
  {
    if (!(std::abs(boundaries[g].sigma) > 0.0))
    {
      continue;
    }
    const std::complex<double> s = EvaluateScalar(g, omega);
    fbr.AddMaterialProperty(mat_op.GetCeedBdrAttributes(boundaries[g].attr_list), s.real());
    fbi.AddMaterialProperty(mat_op.GetCeedBdrAttributes(boundaries[g].attr_list), s.imag());
  }
}

void SurfaceConductivityOperator::AddBoundaryMassBdrCoefficients(
    std::size_t group_idx, MaterialPropertyCoefficient &fb, double coeff) const
{
  // Per-group ω-independent boundary ND-mass with optional scalar scaling. For complex-λ
  // analytic continuation we factor the coefficient out (assemble M_σ_g once) and apply
  // the per-group complex scalar i·ω/Z(ω) with ω = -i·λ at runtime via EvaluateScalar.
  if (group_idx >= boundaries.size() || !(std::abs(boundaries[group_idx].sigma) > 0.0))
  {
    return;
  }
  fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(boundaries[group_idx].attr_list),
                         coeff);
}

std::complex<double>
SurfaceConductivityOperator::EvaluateScalar(std::size_t group_idx,
                                            std::complex<double> omega) const
{
  // Surface-impedance Leontovich BC; complex-ω overload of the formula in
  // AddExtraSystemBdrCoefficients. With ω complex (analytic continuation under
  // ω = -i·λ for upper-half-plane λ), std::sqrt's principal branch maps the BC onto
  // the physical decaying sheet (Re(δ) > 0, Re(Z) > 0). Finite-thickness correction
  // uses std::complex overloads of cosh/cos/sinh/sin (entire functions in the
  // sqrt(ω) variable, no extra branches introduced).
  MFEM_VERIFY(group_idx < boundaries.size(),
              "SurfaceConductivityOperator::EvaluateScalar: group_idx out of range.");
  const auto &bdr = boundaries[group_idx];
  if (!(std::abs(bdr.sigma) > 0.0))
  {
    return {0.0, 0.0};
  }
  // Skin depth, complex when ω is complex. Principal branch of std::sqrt has a cut on
  // the negative real axis; for upper-half-plane λ → ω = -i·λ in the right half-plane,
  // so Re(δ) > 0 and the BC sits on the physical (decaying) sheet.
  const std::complex<double> delta = std::sqrt(2.0 / (bdr.mu * bdr.sigma * omega));
  const std::complex<double> base = 1.0 / (bdr.sigma * delta);
  // Half-space limit: Z = (1+i)·base.
  std::complex<double> Z = std::complex<double>{1.0, 1.0} * base;
  if (bdr.h > 0.0)
  {
    // Finite-thickness HFSS correction. The real-ω implementation applies two distinct
    // real ratios to Re(Z) and Im(Z) separately; analytically continued to complex ω
    // that becomes a uniform complex correction whose real and imaginary parts at
    // ω real reduce to the original ratios. Equivalent form:
    //   Z = base · [(sinh ν + sin ν) + i·(sinh ν − sin ν)] / (cosh ν − cos ν)
    // with ν = h/δ. cosh, cos, sinh, sin of complex ν are entire functions, so no new
    // branch points are introduced.
    const std::complex<double> nu = bdr.h / delta;
    const std::complex<double> den = std::cosh(nu) - std::cos(nu);
    const std::complex<double> shp = std::sinh(nu) + std::sin(nu);
    const std::complex<double> shm = std::sinh(nu) - std::sin(nu);
    Z = base * (shp + std::complex<double>{0.0, 1.0} * shm) / den;
  }
  return std::complex<double>{0.0, 1.0} * omega / Z;
}

}  // namespace palace
