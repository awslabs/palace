// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "dispersiveconductivityoperator.hpp"

#include <cmath>
#include <numbers>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/enum_string.hpp"
#include "utils/units.hpp"

namespace palace
{

using namespace std::complex_literals;

DispersiveConductivityOperator::DispersiveConductivityOperator(
    const std::vector<config::MaterialData> &materials, ProblemType problem_type,
    const Units &units, const MaterialOperator &mat_op)
  : mat_op(mat_op)
{
  // Groups are built from all materials, so that they are the same on all processes.
  for (const auto &data : materials)
  {
    if (!data.complex_sigma || !data.complex_sigma->IsDispersive())
    {
      continue;
    }
    ComplexConductivityModel model(*data.complex_sigma);
    const double omega_ref = data.complex_sigma->reference_freq;
    const auto sigma_ref = model(omega_ref);
    bool active = mat_op.GetCeedAttributes(data.attributes).Size() > 0;
    Mpi::GlobalOr(1, &active, mat_op.GetMesh().GetComm());
    groups.push_back({std::move(model), data.attributes, omega_ref, sigma_ref, active});
  }
  MFEM_VERIFY(groups.empty() || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::BOUNDARYMODE,
              "Dispersive complex conductivity domains are only available for driven and "
              "boundary mode simulations!");
  PrintDomainInfo(units);
}

void DispersiveConductivityOperator::PrintDomainInfo(const Units &units) const
{
  if (groups.empty())
  {
    return;
  }
  auto freq = [&units](double omega)
  {
    return units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
           (2 * std::numbers::pi);
  };
  Mpi::Print("\nConfiguring dispersive complex conductivity at domain attributes:\n");
  for (const auto &group : groups)
  {
    Mpi::Print(
        " {}: {} model, f_ref = {:.3e} GHz, σ1 = {:.3e} S/m",
        fmt::join(group.attributes, ","), ToString(group.model.GetType()),
        freq(group.omega_ref),
        units.Dimensionalize<Units::ValueType::CONDUCTIVITY>(group.sigma_ref.real()));
    if (const double omega_sigma_2 = -group.omega_ref * group.sigma_ref.imag();
        omega_sigma_2 > 0.0)
    {
      // ω σ2 = 1/(μ0 λ²) defines the London depth at the reference frequency.
      Mpi::Print(", λ = {:.3e} m", units.Dimensionalize<Units::ValueType::LENGTH>(
                                       1.0 / std::sqrt(omega_sigma_2)));
    }
    if (group.model.GetType() == ComplexConductivityType::MATTIS_BARDEEN)
    {
      Mpi::Print(", 2Δ/h = {:.3e} GHz", freq(2.0 * group.model.GetGapFrequency()));
    }
    Mpi::Print("\n");
  }
}

std::complex<double> DispersiveConductivityOperator::EvaluateScalar(std::size_t g,
                                                                    double omega) const
{
  const auto &group = groups[g];
  return 1i * omega * (group.model(omega) - group.sigma_ref.real()) +
         group.omega_ref * group.sigma_ref.imag();
}

double DispersiveConductivityOperator::EvaluateKineticScalar(std::size_t g,
                                                             double omega) const
{
  // -dσ2/dω = Im dσ/dω, less ω_ref σ2,ref / ω² of the London anchor.
  const auto &group = groups[g];
  return group.model.Derivative(omega).imag() +
         group.omega_ref * group.sigma_ref.imag() / (omega * omega);
}

double DispersiveConductivityOperator::EvaluateLossScalar(std::size_t g, double omega) const
{
  const auto &group = groups[g];
  return group.model(omega).real() - group.sigma_ref.real();
}

void DispersiveConductivityOperator::AddExtraSystemCoefficients(
    double omega, MaterialPropertyCoefficient &fr, MaterialPropertyCoefficient &fi) const
{
  for (std::size_t g = 0; g < groups.size(); g++)
  {
    const auto attr_list = mat_op.GetCeedAttributes(groups[g].attributes);
    if (attr_list.Size() > 0)
    {
      const auto s = EvaluateScalar(g, omega);
      fr.AddMaterialProperty(attr_list, s.real());
      fi.AddMaterialProperty(attr_list, s.imag());
    }
  }
}

void DispersiveConductivityOperator::AddUnitCoefficients(
    std::size_t g, MaterialPropertyCoefficient &f) const
{
  const auto attr_list = mat_op.GetCeedAttributes(groups[g].attributes);
  if (attr_list.Size() > 0)
  {
    f.AddMaterialProperty(attr_list, 1.0);
  }
}

}  // namespace palace
