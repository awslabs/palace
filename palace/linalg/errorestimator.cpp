// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
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
auto WrapOperator(std::unique_ptr<Operator> &&op);

template <>
auto WrapOperator<Operator>(std::unique_ptr<Operator> &&op)
{
  return std::move(op);
}

template <>
auto WrapOperator<ComplexOperator>(std::unique_ptr<Operator> &&op)
{
  return std::make_unique<ComplexWrapperOperator>(std::move(op), nullptr);
}

auto GetMassMatrix(const FiniteElementSpace &fespace)
{
  constexpr bool skip_zeros = false;
  BilinearForm m(fespace);
  m.AddDomainIntegrator<VectorFEMassIntegrator>();
  return std::make_unique<ParOperator>(m.Assemble(skip_zeros), fespace);
}

template <typename OperType>
auto ConfigureLinearSolver(MPI_Comm comm, double tol, int max_it, int print)
{
  // The system matrix for the projection is real, SPD and diagonally dominant.
  auto pc = std::make_unique<JacobiSmoother<OperType>>();
  auto pcg = std::make_unique<CgSolver<OperType>>(comm, print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  return std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(pc));
}

}  // namespace

template <typename VecType>
FluxProjector<VecType>::FluxProjector(const MaterialOperator &mat_op,
                                      const FiniteElementSpace &nd_fespace, double tol,
                                      int max_it, int print)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    BilinearForm flux(nd_fespace);
    flux.AddDomainIntegrator<MixedVectorCurlIntegrator>(muinv_func);
    Flux = WrapOperator<OperType>(
        std::make_unique<ParOperator>(flux.PartialAssemble(), nd_fespace));
  }
  M = WrapOperator<OperType>(GetMassMatrix(nd_fespace));

  ksp = ConfigureLinearSolver<OperType>(nd_fespace.GetComm(), tol, max_it, print);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(nd_fespace.GetTrueVSize());
  rhs.UseDevice(true);
}

template <typename VecType>
FluxProjector<VecType>::FluxProjector(const MaterialOperator &mat_op,
                                      const FiniteElementSpace &h1_fespace,
                                      const FiniteElementSpace &rt_fespace, double tol,
                                      int max_it, int print)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm flux(h1_fespace, rt_fespace);
    flux.AddDomainIntegrator<MixedVectorGradientIntegrator>(epsilon_func);
    Flux = WrapOperator<OperType>(std::make_unique<ParOperator>(
        flux.PartialAssemble(), h1_fespace, rt_fespace, false));
  }
  M = WrapOperator<OperType>(GetMassMatrix(rt_fespace));

  ksp = ConfigureLinearSolver<OperType>(h1_fespace.GetComm(), tol, max_it, print);
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
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(const MaterialOperator &mat_op,
                                                        FiniteElementSpace &nd_fespace,
                                                        double tol, int max_it, int print)
  : mat_op(mat_op), nd_fespace(nd_fespace),
    projector(mat_op, nd_fespace, tol, max_it, print), F(nd_fespace.GetTrueVSize()),
    F_gf(nd_fespace.GetVSize()), U_gf(nd_fespace.GetVSize())
{
  F.UseDevice(true);
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
          // μ⁻¹ ∇ × U -- μ⁻¹ |J|⁻¹J ∇ × Uᵣ
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
                                               FiniteElementSpace &h1_fespace, double tol,
                                               int max_it, int print)
  : mat_op(mat_op), h1_fespace(h1_fespace),
    rt_fec(std::make_unique<mfem::RT_FECollection>(h1_fespace.GetFEColl().GetOrder() - 1,
                                                   h1_fespace.SpaceDimension())),
    rt_fespace(std::make_unique<FiniteElementSpace>(h1_fespace.GetMesh(), rt_fec.get())),
    projector(mat_op, h1_fespace, *rt_fespace, tol, max_it, print),
    F(rt_fespace->GetTrueVSize()), F_gf(rt_fespace->GetVSize()), U_gf(h1_fespace.GetVSize())
{
  F.UseDevice(true);
}

void GradFluxErrorEstimator::AddErrorIndicator(const Vector &U,
                                               ErrorIndicator &indicator) const
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(U, F);
  h1_fespace.GetProlongationMatrix()->Mult(U, U_gf);
  rt_fespace->GetProlongationMatrix()->Mult(F, F_gf);
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
      const mfem::FiniteElement &rt_fe = *rt_fespace->Get().GetFE(e);
      mesh.GetElementTransformation(e, &T);
      h1_fespace.Get().GetElementDofs(e, h1_dofs);
      rt_fespace->Get().GetElementDofs(e, rt_dofs);
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

        // ε ∇U -- ε J⁻¹ ∇Uᵣ
        U_gf.GetSubVector(h1_dofs, loc_gf);
        Grad.MultTranspose(loc_gf, V_ip);
        T.InverseJacobian().MultTranspose(V_ip, V_smooth);
        mat_op.GetPermittivityReal(T.Attribute).Mult(V_smooth, V_ip);

        // Smooth flux -- |J|⁻¹J Fᵣ
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
