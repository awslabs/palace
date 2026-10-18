// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
#define PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP

#include <functional>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class DispersiveConductivityOperator;
class GridFunction;
class FiniteElementSpace;
class IoData;
class MaterialOperator;

namespace config
{

struct DomainPostData;

}  // namespace config

//
// Class to handle domain energy postprocessing. We use a leading factor of 1/2 instead of
// 1/4 even though the eigenmodes are peak phasors and not RMS normalized because the same
// peak phasors are used to compute the voltages/currents which are 2x the time-averaged
// values. This correctly yields an EPR of 1 in cases where expected.
//
class DomainPostOperator
{
public:
  // Temporary vectors for inner product calculations.
  mutable Vector D, H;

  // Bilinear forms for computing field energy integrals over domains.
  std::unique_ptr<Operator> M_elec, M_mag;
  std::map<int, std::pair<std::unique_ptr<Operator>, std::unique_ptr<Operator>>> M_i;

  // Mass operators whose weighted quadratic forms give the kinetic energy 1/2 ∫ -dσ2/dω
  // |E|² and conduction loss 1/2 ∫ σ1 |E|²: the London and conductivity coefficients, and
  // the unit mass of each dispersive group weighted by its residuals.
  struct ConductionForm
  {
    std::unique_ptr<Operator> M;
    std::map<int, std::unique_ptr<Operator>> M_i;
    std::function<std::pair<double, double>(double)> weights;  // (Kinetic, loss) at ω.
  };
  std::vector<ConductionForm> cond_forms;
  bool has_kinetic = false, has_loss = false;

  DomainPostOperator() = default;
  // The kinetic energy and conduction loss are postprocessed when disp_op is given.
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace,
                     const DispersiveConductivityOperator *disp_op = nullptr);
  DomainPostOperator(const config::DomainPostData &postpro, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace,
                     const FiniteElementSpace &rt_fespace);
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &fespace);

  // Mode analysis: ND space for electric field energy only (no B-field).
  // Uses VectorFEMassIntegrator with ε (not CurlCurlIntegrator).
  DomainPostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const FiniteElementSpace &nd_fespace, bool electric_energy_only);

  // Get volume integrals computing the electric or magnetic field energy in the entire
  // domain.
  double GetElectricFieldEnergy(const GridFunction &E) const;
  double GetMagneticFieldEnergy(const GridFunction &B) const;

  // Get volume integrals for the electric or magnetic field energy in a portion of the
  // domain.
  double GetDomainElectricFieldEnergy(int idx, const GridFunction &E) const;
  double GetDomainMagneticFieldEnergy(int idx, const GridFunction &B) const;

  // Get the kinetic energy and conduction loss at angular frequency ω in the entire domain
  // or in a portion of the domain.
  std::pair<double, double> GetKineticEnergyAndLoss(const GridFunction &E,
                                                    double omega) const;
  std::pair<double, double> GetDomainKineticEnergyAndLoss(int idx, const GridFunction &E,
                                                          double omega) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DOMAIN_POST_OPERATOR_HPP
