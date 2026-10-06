// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
#define PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace
{

class MumpsSchurSolver;
class SpaceOperator;

//
// Exact per-frequency condensation of the environment of a driven problem: the
// environment operator A_E(ω) = K_E + iω C_E - ω² M_E + A2_E(ω), assembled with Palace's
// operators restricted to the environment, is condensed onto the interface Γ by a partial
// MUMPS factorization of its real symmetric form [[Ar, Ai], [Ai, -Ar]], whose Schur
// complement is the same form of S_E(ω).
//
class DrivenSubstructure
{
public:
  DrivenSubstructure(SpaceOperator &space_op, const std::vector<int> &region_attributes,
                     const std::vector<int> &environment_attributes);
  ~DrivenSubstructure();

  // Condense the environment at the angular frequency ω (nondimensional). The analysis of
  // the first condensation is reused at the next frequencies. Collective.
  void Condense(double omega);

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
  std::vector<int> env_attrs;
  std::vector<char> is_env_int, is_gamma;
  mfem::Array<int> other;  // local true DOFs pinned in the environment operator
  std::vector<HYPRE_BigInt> gamma_tdofs;
  std::unique_ptr<ComplexOperator> K, C, M;
  std::unique_ptr<MumpsSchurSolver> schur;
  std::vector<std::complex<double>> S;

  // The real form of A_E(ω), with the DOFs outside the environment interior and Γ pinned.
  std::unique_ptr<mfem::HypreParMatrix> BlockOperator(double omega);
};

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
