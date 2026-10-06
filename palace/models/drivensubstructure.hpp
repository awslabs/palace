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

template <typename T>
class MumpsSchurSolverT;
class SpaceOperator;

//
// Exact per-frequency substructuring of a driven problem: the environment operator
// A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's operators restricted
// to the environment, is condensed onto the interface Γ by a partial complex symmetric
// MUMPS factorization, which gives S_E(ω). The region is condensed onto Γ in the same way,
// and the interface system S_R(ω) + S_E(ω) is factored densely.
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
  std::vector<char> is_env_int, is_gamma, is_region_int;
  // Local true DOFs pinned in the environment operator (outside its interior and Γ) and in
  // the region operator (outside the region interior and Γ).
  mfem::Array<int> other_env, other_region;
  std::vector<HYPRE_BigInt> gamma_tdofs;
  std::vector<int> gamma_cnt, gamma_disp;  // interface DOFs per rank
  std::unique_ptr<ComplexOperator> K_env, C_env, M_env, K_region, C_region, M_region;
  std::unique_ptr<MumpsSchurSolverT<std::complex<double>>> env_schur, region_schur;
  // Real and imaginary parts of a side operator, A_E(ω) or A_R(ω) (until factored).
  struct Parts
  {
    std::unique_ptr<mfem::HypreParMatrix> real, imag;
  };
  Parts env_op, region_op;
  // On rank 0: S_E, and the factored interface system S_R + S_E with its pivots.
  std::vector<std::complex<double>> S, T;
  std::vector<int> T_piv;

  // K + iω C - ω² M + A2(ω) for one side, with the given DOFs pinned.
  Parts SideOperator(double omega, const std::vector<int> &attrs, const ComplexOperator *K,
                     const ComplexOperator *C, const ComplexOperator *M,
                     const mfem::Array<int> &pinned);
};

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
