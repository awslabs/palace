// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP
#define PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP

#include <complex>
#include <map>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class DrivenSubstructure;
struct DrivenSubstructureModel;
class IoData;
struct SignatureMap;

namespace config
{

struct BoundaryData;
struct DomainData;
struct SubstructuringData;

}  // namespace config

//
// Helpers of the model of driven substructuring, for its offline and online sweeps.
//

// The lumped ports on the environment side (and its wave ports, as -index), and an abort
// for lumped ports touching Γ (their excitation and voltage would depend on both sides).
std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::map<int, std::vector<char>> &port_dofs);

// The parameters of the environment's ports that enter S_E (see
// DrivenSubstructureModel::port_data).
std::vector<double> EnvironmentPortData(const IoData &iodata,
                                        const std::vector<int> &env_ports);

// The full S_E (n x n, column-major) of a record's lower triangle (by columns), mapped to
// the current interface as M^-T S_E M^-1 with map (rows[g]: row g of M^-T), on one
// triangle.
std::vector<std::complex<double>> FullSchur(const std::complex<double> *lower, int n,
                                            const SignatureMap *map);

// The index of x in y (relative tolerance tol), -1 if none.
int FindValue(const std::vector<double> &y, double x, double tol = 1.0e-10);

// From a model, the field is known in the region (and on Γ) only: the domain and boundary
// postprocessing restricted to energies of region domains, without totals or participation
// ratios, and to the surface fluxes and probes the region holds, with a warning of what is
// skipped (and that the environment's ports, if any, have no power). mode names the run in
// the warning. Collective.
void RestrictPostprocessingToRegion(const mfem::ParMesh &mesh,
                                    const config::SubstructuringData &sub,
                                    const DrivenSubstructure &ds,
                                    const DrivenSubstructureModel &model, bool env_ports,
                                    const char *mode, config::DomainData &domains,
                                    config::BoundaryData &boundaries);

}  // namespace palace

#endif  // PALACE_MODELS_CONDENSED_ENVIRONMENT_HPP
