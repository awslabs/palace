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
  std::vector<char> is_env_int, is_gamma, is_region_int;
  std::vector<HYPRE_BigInt> gamma_tdofs;
  std::vector<int> gamma_cnt, gamma_disp;  // interface DOFs per rank

  // The operator K + iω C - ω² M + A2(ω) of one side (environment or region), with its
  // pinned DOFs, factored by MUMPS: the lower-triangle pattern, the entries of the
  // frequency-independent parts in pattern order (the assembled operators are not kept),
  // and the positions of the unit diagonal of the pinned DOFs.
  struct Side
  {
    std::vector<int> attrs;
    // Local true DOFs pinned in the side's operator: outside its interior and Γ.
    mfem::Array<int> pinned;
    std::vector<int> row_ptr, unit;
    std::vector<std::vector<double>> parts;  // Kr, Ki, Cr, Ci, Mr, Mi (empty if absent)
    bool extra = false;                      // A2(ω) is present
    // The local entries (1-based COO) at the first frequency, until factored.
    std::vector<int> irn, jcn;
    std::vector<std::complex<double>> val;
    std::unique_ptr<MumpsSchurSolverT<std::complex<double>>> schur;
  };
  Side env, region;
  // On rank 0: S_E, and the factored interface system S_R + S_E with its pivots.
  std::vector<std::complex<double>> S, T;
  std::vector<int> T_piv;

  // The assembled parts of A2(ω) of a side, with the pinned DOFs eliminated.
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> ExtraParts(const Side &side,
                                                                double omega);

  // The pattern and entries of a side at the first frequency.
  void Setup(Side &side, double omega);

  // Factor a side (set up at ω), or refactor it at another frequency.
  void Factor(Side &side, double omega);

  // The entries of a side at ω on its pattern (columns jcn), given its A2(ω).
  void Fill(const Side &side, double omega,
            const std::vector<std::unique_ptr<mfem::HypreParMatrix>> &extra,
            const std::vector<int> &jcn, std::vector<std::complex<double>> &val) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_DRIVEN_SUBSTRUCTURE_HPP
