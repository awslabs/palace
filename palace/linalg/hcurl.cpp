// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hcurl.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a,
                           const FiniteElementSpace &fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&a,
                                     const FiniteElementSpace &fespace)
{
  return std::make_unique<ParOperator>(std::move(a), fespace);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&a,
                                            const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(a), nullptr, fespace);
}

}  // namespace

template <typename VecType>
WeightedHCurlNormSolver<VecType>::WeightedHCurlNormSolver(
    const MaterialOperator &mat_op, FiniteElementSpaceHierarchy &nd_fespaces,
    FiniteElementSpaceHierarchy &h1_fespaces,
    const std::vector<mfem::Array<int>> &nd_dbc_tdof_lists,
    const std::vector<mfem::Array<int>> &h1_dbc_tdof_lists, double tol, int max_it,
    int print)
{
  MFEM_VERIFY(h1_fespaces.GetNumLevels() == nd_fespaces.GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");
  const auto n_levels = nd_fespaces.GetNumLevels();
  {
    constexpr bool skip_zeros = false;
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm a(nd_fespaces.GetFinestFESpace()), a_aux(h1_fespaces.GetFinestFESpace());
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_func, epsilon_func);
    a_aux.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
    // a.AssembleQuadratureData();
    // a_aux.AssembleQuadratureData();
    auto a_vec = a.Assemble(nd_fespaces, skip_zeros);
    auto a_aux_vec = a_aux.Assemble(h1_fespaces, skip_zeros);
    auto A_mg = std::make_unique<BaseMultigridOperator<OperType>>(n_levels);
    for (bool aux : {false, true})
    {
      for (std::size_t l = 0; l < n_levels; l++)
      {
        const auto &fespace_l =
            aux ? h1_fespaces.GetFESpaceAtLevel(l) : nd_fespaces.GetFESpaceAtLevel(l);
        const auto &dbc_tdof_lists_l = aux ? h1_dbc_tdof_lists[l] : nd_dbc_tdof_lists[l];
        auto A_l = BuildLevelParOperator<OperType>(std::move(aux ? a_aux_vec[l] : a_vec[l]),
                                                   fespace_l);
        A_l->SetEssentialTrueDofs(dbc_tdof_lists_l, Operator::DiagonalPolicy::DIAG_ONE);
        if (aux)
        {
          A_mg->AddAuxiliaryOperator(std::move(A_l));
        }
        else
        {
          A_mg->AddOperator(std::move(A_l));
        }
      }
    }
    A = std::move(A_mg);
  }

  // The system matrix K + M is real and SPD. We use Hypre's AMS solver as the coarse-level
  // multigrid solve.
  auto ams = std::make_unique<MfemWrapperSolver<OperType>>(std::make_unique<HypreAmsSolver>(
      nd_fespaces.GetFESpaceAtLevel(0), h1_fespaces.GetFESpaceAtLevel(0), 1, 1, false, true,
      false, 0));
  std::unique_ptr<Solver<OperType>> pc;
  if (n_levels > 1)
  {
    const auto G = nd_fespaces.GetDiscreteInterpolators(h1_fespaces);
    const int mg_smooth_order =
        std::max(nd_fespaces.GetFinestFESpace().GetMaxElementOrder(), 2);
    pc = std::make_unique<GeometricMultigridSolver<OperType>>(
        nd_fespaces.GetFinestFESpace().GetComm(), std::move(ams),
        nd_fespaces.GetProlongationOperators(), &G, 1, 1, mg_smooth_order, 1.0, 0.0, true);
  }
  else
  {
    pc = std::move(ams);
  }

  auto pcg =
      std::make_unique<CgSolver<OperType>>(nd_fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(pc));
  ksp->SetOperators(*A, *A);
}

template class WeightedHCurlNormSolver<Vector>;
template class WeightedHCurlNormSolver<ComplexVector>;

}  // namespace palace
