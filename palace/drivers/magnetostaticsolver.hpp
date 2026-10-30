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
                            const std::vector<Vector> &london_ms_shifted) const;

  // Inductance, reluctance and mutual matrices from the cross-energies of the current and
  // flux-loop excitations (currents first), written to terminal-M.csv, terminal-Minv.csv
  // and terminal-Mm.csv with the excitations.
  void PostprocessInductance(const std::vector<int> &current_idx,
                             const std::vector<int> &flux_idx,
                             const mfem::DenseMatrix &cross_energy,
                             const std::vector<double> &I_inc,
                             const std::vector<double> &Phi_inc,
                             const mfem::DenseMatrix &linked_flux) const;

  // Whether each excitation column has reciprocal mutuals (current ports Open when
  // inactive, and every flux loop).
  std::vector<bool> ReciprocalColumns(const std::vector<int> &current_idx, int n) const;

  // Report the stored energy and imposed loop fluxoids of each excitation, when grouped
  // flux loop excitations leave no inductance matrix.
  void PostprocessExcitations(const CurlCurlOperator &curlcurl_op,
                              const std::vector<Vector> &london_ah,
                              const std::vector<double> &domain_energy,
                              const std::vector<double> &london_kinetic) const;

  // Write terminal-Phi.csv with the stored energy of each excitation (currents first) and
  // the fluxoid it imposes on each flux loop.
  void PrintExcitations(const std::vector<int> &excitation_idx,
                        const std::vector<int> &loop_idx, const std::vector<double> &energy,
                        const mfem::DenseMatrix &fluxoid) const;

  std::pair<ErrorIndicator, long long int>
  Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const override;

public:
  using BaseSolver::BaseSolver;

  // An online substructured run takes the reference length of its saved model, so that a
  // mesh of the region alone is nondimensionalized as the model was.
  void Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                  MPI_Comm comm) const override;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_MAGNETOSTATIC_SOLVER_HPP
