// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
#define PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class MumpsSchurSolver;
class SpaceOperator;

//
// Exact per-frequency substructuring of a driven problem: the environment operator
// A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's operators restricted
// to the environment, is condensed onto the interface Γ by a partial MUMPS factorization
// of its real symmetric form [[Ar, Ai], [Ai, -Ar]], whose Schur complement is the same
// form of S_E(ω). The region operator plus S_E(ω) on Γ is then factored in the same form.
//
class DrivenSubstructure
{
public:
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes);
  ~DrivenSubstructure();

  // Condense the environment and factor the region against it at the angular frequency ω
  // (nondimensional). The analyses of the first frequency are reused at the next ones.
  // Collective.
  void Condense(double omega);

  // Solve the coupled problem at the frequency of the last condensation for a batch of
  // right-hand sides (full true-DOF vectors): the region with the interface loads of the
  // environment sources, then the environment interior. Collective.
  void Solve(const std::vector<const ComplexVector *> &rhs, std::vector<ComplexVector> &u);

  // S_E(ω) of the last condensation on rank 0 (|Γ| x |Γ|, column-major, complex
  // symmetric), in interface order.
  const std::vector<std::complex<double>> &Schur() const { return S; }

  // Global true DOFs of the interface in interface order (replicated), and its size.
  const std::vector<HYPRE_BigInt> &InterfaceTrueDofs() const { return gamma_tdofs; }
  int InterfaceSize() const { return static_cast<int>(gamma_tdofs.size()); }

  // Local true-DOF classification: environment interior, interface (excluding the
  // Dirichlet DOFs).
  const std::vector<char> &EnvironmentInterior() const { return is_env_int; }
  const std::vector<char> &Interface() const { return is_gamma; }

private:
  SpaceOperator &space_op;
  std::vector<int> region_attrs, env_attrs;
  std::vector<char> is_env_int, is_gamma, is_region_free;
  // Local true DOFs pinned in the environment operator (outside its interior and Γ) and in
  // the region operator (outside the region-free DOFs, which include Γ).
  mfem::Array<int> other_env, other_region;
  std::vector<HYPRE_BigInt> gamma_tdofs;
  int gamma_off = 0, gamma_nloc = 0;
  std::unique_ptr<ComplexOperator> K_env, C_env, M_env, K_region, C_region, M_region;
  std::unique_ptr<MumpsSchurSolver> env_schur, region_lu;
  std::unique_ptr<mfem::HypreParMatrix> env_op;  // real form of A_E(ω), last frequency
  std::vector<std::complex<double>> S, S_rows;   // S_E on rank 0, and this rank's rows

  // The real form of K + iω C - ω² M + A2(ω) for one side (+ X on Γ), with the given DOFs
  // pinned.
  std::unique_ptr<mfem::HypreParMatrix>
  BlockOperator(double omega, const std::vector<int> &attrs, const ComplexOperator *K,
                const ComplexOperator *C, const ComplexOperator *M,
                const mfem::Array<int> &pinned, const mfem::HypreParMatrix *Xr = nullptr,
                const mfem::HypreParMatrix *Xi = nullptr);

  // S_E on Γ as parent-space matrices (dense Γ x Γ block): real and imaginary parts.
  std::pair<std::unique_ptr<mfem::HypreParMatrix>, std::unique_ptr<mfem::HypreParMatrix>>
  InterfaceMatrices() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
