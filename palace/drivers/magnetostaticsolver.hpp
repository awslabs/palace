// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
#define PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP

#include <memory>
#include <vector>
#include "drivers/basesolver.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class CurlCurlOperator;
class ErrorIndicator;
class Mesh;
template <ProblemType>
class PostOperator;
class SurfaceCurrentOperator;
class SurfaceFluxOperator;

//
// Driver class for magnetostatic simulations.
//
class MagnetostaticSolver : public BaseSolver
{
private:
  void PostprocessTerminals(PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                            const SurfaceCurrentOperator &surf_j_op,
                            const SurfaceFluxOperator &surf_flux_op,
                            const std::vector<Vector> &A, const std::vector<double> &I_inc,
                            const std::vector<double> &Phi_inc,
                            const mfem::DenseMatrix &linked_flux,
                            const std::vector<Vector> &london_ah,
                            const std::vector<Vector> &london_ms_shifted,
                            const std::vector<Vector> &london_ml) const;

  // Report the stored energy and imposed loop fluxoids of each excitation, when grouped
  // flux loop excitations leave no inductance matrix.
  void PostprocessExcitations(const CurlCurlOperator &curlcurl_op,
                              const std::vector<Vector> &london_ah,
                              const std::vector<double> &domain_energy,
                              const std::vector<double> &london_kinetic) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
