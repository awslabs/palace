// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/omp.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a,
                           const FiniteElementSpace &trial_fespace,
                           const FiniteElementSpace &test_fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&a,
                                     const FiniteElementSpace &trial_fespace,
                                     const FiniteElementSpace &test_fespace)
{
  return std::make_unique<ParOperator>(std::move(a), trial_fespace, test_fespace, false);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&a,
                                            const FiniteElementSpace &trial_fespace,
                                            const FiniteElementSpace &test_fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(a), nullptr, trial_fespace,
                                              test_fespace, false);
}

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a, const FiniteElementSpace &fespace)
{
  return BuildLevelParOperator<OperType>(std::move(a), fespace, fespace);
}

template <typename OperType>
std::unique_ptr<OperType> GetMassMatrix(const FiniteElementSpaceHierarchy &fespaces,
                                        bool use_mg)
{
  constexpr bool skip_zeros = false;
  BilinearForm m(fespaces.GetFinestFESpace());
  m.AddDomainIntegrator<VectorFEMassIntegrator>();
  if (!use_mg)
  {
    return BuildLevelParOperator<OperType>(m.Assemble(skip_zeros),
                                           fespaces.GetFinestFESpace());
  }
  else
  {
    auto m_vec = m.Assemble(fespaces, skip_zeros);
    auto M_mg = std::make_unique<BaseMultigridOperator<OperType>>(fespaces.GetNumLevels());
    for (std::size_t l = 0; l < fespaces.GetNumLevels(); l++)
    {
      const auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
      M_mg->AddOperator(BuildLevelParOperator<OperType>(std::move(m_vec[l]), fespace_l));
    }
    return M_mg;
  }
}

template <typename OperType>
auto ConfigureLinearSolver(const FiniteElementSpaceHierarchy &fespaces, double tol,
                           int max_it, int print, bool use_mg)
{
  // The system matrix for the projection is real, SPD and diagonally dominant.
  std::unique_ptr<Solver<OperType>> pc;
  if (!use_mg)
  {
    pc = std::make_unique<JacobiSmoother<OperType>>();
  }
  else
  {
    auto amg = std::make_unique<BoomerAmgSolver>(1, 1, 0);
    amg->SetStrengthThresh(0.8);  // More coarsening to save memory
    const int mg_smooth_order =
        std::max(fespaces.GetFinestFESpace().GetMaxElementOrder(), 2) +
        (dynamic_cast<const mfem::RT_FECollection *>(
             &fespaces.GetFinestFESpace().GetFEColl()) != nullptr);
    pc = std::make_unique<GeometricMultigridSolver<OperType>>(
        fespaces.GetFinestFESpace().GetComm(),
        std::make_unique<MfemWrapperSolver<OperType>>(std::move(amg)),
        fespaces.GetProlongationOperators(), nullptr, 1, 1, mg_smooth_order, 1.0, 0.0,
        true);
  }

  auto pcg =
      std::make_unique<CgSolver<OperType>>(fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  return std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(pc));
}

}  // namespace

template <typename VecType>
FluxProjector<VecType>::FluxProjector(const MaterialOperator &mat_op,
                                      const FiniteElementSpaceHierarchy &nd_fespaces,
                                      double tol, int max_it, int print, bool use_mg)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  use_mg = use_mg && (nd_fespaces.GetNumLevels() > 1);
  const auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    BilinearForm flux(nd_fespace);
    flux.AddDomainIntegrator<MixedVectorCurlIntegrator>(muinv_func);
    Flux = BuildLevelParOperator<OperType>(flux.PartialAssemble(), nd_fespace);
  }
  M = GetMassMatrix<OperType>(nd_fespaces, use_mg);
  ksp = ConfigureLinearSolver<OperType>(nd_fespaces, tol, max_it, print, use_mg);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(nd_fespace.GetTrueVSize());
  rhs.UseDevice(true);
}

template <typename VecType>
FluxProjector<VecType>::FluxProjector(const MaterialOperator &mat_op,
                                      const FiniteElementSpace &h1_fespace,
                                      const FiniteElementSpaceHierarchy &rt_fespaces,
                                      double tol, int max_it, int print, bool use_mg)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  use_mg = use_mg && (rt_fespaces.GetNumLevels() > 1);
  const auto &rt_fespace = rt_fespaces.GetFinestFESpace();
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm flux(h1_fespace, rt_fespace);
    flux.AddDomainIntegrator<MixedVectorGradientIntegrator>(epsilon_func);
    Flux = BuildLevelParOperator<OperType>(flux.PartialAssemble(), h1_fespace, rt_fespace);
  }
  M = GetMassMatrix<OperType>(rt_fespaces, use_mg);
  ksp = ConfigureLinearSolver<OperType>(rt_fespaces, tol, max_it, print, use_mg);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(rt_fespace.GetTrueVSize());
  rhs.UseDevice(true);
}

template <typename VecType>
void FluxProjector<VecType>::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::SOLVE_ESTIMATOR);
  MFEM_ASSERT(y.Size() == rhs.Size(), "Invalid vector dimensions for FluxProjector::Mult!");
  Flux->Mult(x, rhs);
  // Mpi::Print(" Computing smooth flux projection for error estimation\n");
  ksp->Mult(rhs, y);
}

template <typename VecType>
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpaceHierarchy &nd_fespaces, double tol,
    int max_it, int print, bool use_mg)
  : mat_op(mat_op), nd_fespace(nd_fespaces.GetFinestFESpace()),
    projector(mat_op, nd_fespaces, tol, max_it, print, use_mg), U_gf(nd_fespace.GetVSize()),
    F(nd_fespace.GetTrueVSize()), F_gf(nd_fespace.GetVSize())
{
  U_gf.UseDevice(true);
  F.UseDevice(true);
  F_gf.UseDevice(true);
}

template <typename VecType>
void CurlFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &U,
                                                        ErrorIndicator &indicator) const
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(U, F);
  if constexpr (std::is_same<VecType, ComplexVector>::value)
  {
    nd_fespace.GetProlongationMatrix()->Mult(U.Real(), U_gf.Real());
    nd_fespace.GetProlongationMatrix()->Mult(U.Imag(), U_gf.Imag());
    nd_fespace.GetProlongationMatrix()->Mult(F.Real(), F_gf.Real());
    nd_fespace.GetProlongationMatrix()->Mult(F.Imag(), F_gf.Imag());
    U_gf.Real().HostRead();
    U_gf.Imag().HostRead();
    F_gf.Real().HostRead();
    F_gf.Imag().HostRead();
  }
  else
  {
    nd_fespace.GetProlongationMatrix()->Mult(U, U_gf);
    nd_fespace.GetProlongationMatrix()->Mult(F, F_gf);
    U_gf.HostRead();
    F_gf.HostRead();
  }

  // Loop over elements and accumulate the estimates from this component. The discontinuous
  // flux is μ⁻¹ ∇ × U.
  const auto &mesh = nd_fespace.GetParMesh();
  Vector estimates(mesh.GetNE());
  auto *h_estimates = estimates.HostWrite();
  double norm2 = 0.0;
  PalacePragmaOmp(parallel reduction(+ : norm2))
  {
    // Assuming dim == space_dim == curl_dim
    mfem::IsoparametricTransformation T;
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    mfem::Vector V_ip(mesh.SpaceDimension()), V_smooth(mesh.SpaceDimension()),
        V_tmp(mesh.SpaceDimension()), loc_gf;
    mfem::DenseMatrix Interp, Curl;

    double loc_norm2 = 0.0;
    PalacePragmaOmp(for schedule(static))
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const mfem::FiniteElement &fe = *nd_fespace.Get().GetFE(e);
      mesh.GetElementTransformation(e, &T);
      nd_fespace.Get().GetElementDofs(e, dofs, dof_trans);
      Interp.SetSize(fe.GetDof(), V_ip.Size());
      Curl.SetSize(fe.GetDof(), V_ip.Size());
      const int q_order = fem::DefaultIntegrationOrder::Get(T);
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(mesh.GetElementGeometry(e), q_order);

      double elem_err = 0.0;
      for (int i = 0; i < ir.GetNPoints(); i++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        T.SetIntPoint(&ip);
        fe.CalcVShape(ip, Interp);
        fe.CalcCurlShape(ip, Curl);
        const double w = ip.weight * T.Weight();

        auto AccumulateError = [&](const Vector &U_gf_, const Vector &F_gf_)
        {
          // μ⁻¹ ∇ × U -- μ⁻¹ |J|⁻¹ J ∇ × Uᵣ
          U_gf_.GetSubVector(dofs, loc_gf);
          if (dof_trans.GetDofTransformation())
          {
            dof_trans.InvTransformPrimal(loc_gf);
          }
          Curl.MultTranspose(loc_gf, V_ip);
          T.Jacobian().Mult(V_ip, V_smooth);
          mat_op.GetInvPermeability(T.Attribute).Mult(V_smooth, V_ip);
          V_ip *= 1.0 / T.Weight();

          // Smooth flux -- J⁻ᵀ Fᵣ
          F_gf_.GetSubVector(dofs, loc_gf);
          if (dof_trans.GetDofTransformation())
          {
            dof_trans.InvTransformPrimal(loc_gf);
          }
          Interp.MultTranspose(loc_gf, V_tmp);
          T.InverseJacobian().MultTranspose(V_tmp, V_smooth);

          V_smooth -= V_ip;
          elem_err += w * (V_smooth * V_smooth);
          return w * (V_ip * V_ip);
        };
        if constexpr (std::is_same<VecType, ComplexVector>::value)
        {
          loc_norm2 += AccumulateError(U_gf.Real(), F_gf.Real());
          loc_norm2 += AccumulateError(U_gf.Imag(), F_gf.Imag());
        }
        else
        {
          loc_norm2 += AccumulateError(U_gf, F_gf);
        }
      }
      h_estimates[e] = std::sqrt(elem_err);
    }
    norm2 += loc_norm2;
  }
  estimates.UseDevice(true);

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &norm2, mesh.GetComm());
  if (norm2 > 0.0)
  {
    estimates *= 1.0 / std::sqrt(norm2);
  }
  indicator.AddIndicator(estimates);
}

GradFluxErrorEstimator::GradFluxErrorEstimator(const MaterialOperator &mat_op,
                                               FiniteElementSpace &h1_fespace,
                                               FiniteElementSpaceHierarchy &rt_fespaces,
                                               double tol, int max_it, int print,
                                               bool use_mg)
  : mat_op(mat_op), h1_fespace(h1_fespace), rt_fespace(rt_fespaces.GetFinestFESpace()),
    projector(mat_op, h1_fespace, rt_fespaces, tol, max_it, print, use_mg),
    U_gf(h1_fespace.GetVSize()), F(rt_fespace.GetTrueVSize()), F_gf(rt_fespace.GetVSize())
{
  U_gf.UseDevice(true);
  F.UseDevice(true);
  F_gf.UseDevice(true);
}

void GradFluxErrorEstimator::AddErrorIndicator(const Vector &U,
                                               ErrorIndicator &indicator) const
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(U, F);
  h1_fespace.GetProlongationMatrix()->Mult(U, U_gf);
  rt_fespace.GetProlongationMatrix()->Mult(F, F_gf);
  U_gf.HostRead();
  F_gf.HostRead();

  // Loop over elements and accumulate the estimates from this component. The discontinuous
  // flux is ε ∇U.
  const auto &mesh = h1_fespace.GetParMesh();
  Vector estimates(mesh.GetNE());
  auto *h_estimates = estimates.HostWrite();
  double norm2 = 0.0;
  PalacePragmaOmp(parallel reduction(+ : norm2))
  {
    // Assuming dim == space_dim
    mfem::IsoparametricTransformation T;
    mfem::Array<int> h1_dofs, rt_dofs;
    mfem::Vector V_ip(mesh.SpaceDimension()), V_smooth(mesh.SpaceDimension()), loc_gf,
        V_tmp(mesh.SpaceDimension());
    mfem::DenseMatrix Interp, Grad;

    double loc_norm2 = 0.0;
    PalacePragmaOmp(for schedule(static))
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const mfem::FiniteElement &h1_fe = *h1_fespace.Get().GetFE(e);
      const mfem::FiniteElement &rt_fe = *rt_fespace.Get().GetFE(e);
      mesh.GetElementTransformation(e, &T);
      h1_fespace.Get().GetElementDofs(e, h1_dofs);
      rt_fespace.Get().GetElementDofs(e, rt_dofs);
      Interp.SetSize(rt_fe.GetDof(), V_ip.Size());
      Grad.SetSize(h1_fe.GetDof(), V_ip.Size());
      const int q_order = fem::DefaultIntegrationOrder::Get(T);
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(mesh.GetElementGeometry(e), q_order);

      double elem_err = 0.0;
      for (int i = 0; i < ir.GetNPoints(); i++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        T.SetIntPoint(&ip);
        rt_fe.CalcVShape(ip, Interp);
        h1_fe.CalcDShape(ip, Grad);
        const double w = ip.weight * T.Weight();

        // ε ∇U -- ε J⁻ᵀ ∇Uᵣ
        U_gf.GetSubVector(h1_dofs, loc_gf);
        Grad.MultTranspose(loc_gf, V_ip);
        T.InverseJacobian().MultTranspose(V_ip, V_smooth);
        mat_op.GetPermittivityReal(T.Attribute).Mult(V_smooth, V_ip);

        // Smooth flux -- |J|⁻¹ J Fᵣ
        F_gf.GetSubVector(rt_dofs, loc_gf);
        Interp.MultTranspose(loc_gf, V_tmp);
        T.Jacobian().Mult(V_tmp, V_smooth);
        V_smooth /= T.Weight();

        V_smooth -= V_ip;
        elem_err += w * (V_smooth * V_smooth);
        loc_norm2 += w * (V_ip * V_ip);
      }
      h_estimates[e] = std::sqrt(elem_err);
    }
    norm2 += loc_norm2;
  }
  estimates.UseDevice(true);

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &norm2, mesh.GetComm());
  if (norm2 > 0.0)
  {
    estimates *= 1.0 / std::sqrt(norm2);
  }
  indicator.AddIndicator(estimates);
}

template class FluxProjector<Vector>;
template class FluxProjector<ComplexVector>;
template class CurlFluxErrorEstimator<Vector>;
template class CurlFluxErrorEstimator<ComplexVector>;

}  // namespace palace
