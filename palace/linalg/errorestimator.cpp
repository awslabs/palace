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

auto GetMassMatrix(const FiniteElementSpace &fespace)
{
  constexpr bool skip_zeros = false;
  const int dim = fespace.Dimension();
  const auto type = fespace.GetFEColl().GetRangeType(dim);
  BilinearForm m(fespace);
  if (type == mfem::FiniteElement::SCALAR)
  {
    MFEM_VERIFY(fespace.GetVDim() == 1,
                "Scalar mass matrix hierarchy assumes a component-wise solve.");
    m.AddDomainIntegrator<MassIntegrator>();
  }
  else
  {
    m.AddDomainIntegrator<VectorFEMassIntegrator>();
  }
  return std::make_unique<ParOperator>(m.Assemble(skip_zeros), fespace);
}

auto ConfigureLinearSolver(MPI_Comm comm, double tol, int max_it, int print)
{
  // The system matrix for the projection is real, SPD and diagonally dominant.
  auto pc = std::make_unique<JacobiSmoother<Operator>>();
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  return std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
}

}  // namespace

FluxProjector::FluxProjector(const MaterialOperator &mat_op,
                             const FiniteElementSpace &nd_fespace, double tol, int max_it,
                             int print)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    BilinearForm flux(nd_fespace);
    flux.AddDomainIntegrator<MixedVectorCurlIntegrator>(muinv_func);
    Flux = std::make_unique<ParOperator>(flux.PartialAssemble(), nd_fespace);
  }
  M = GetMassMatrix(nd_fespace);

  ksp = ConfigureLinearSolver(nd_fespace.GetComm(), tol, max_it, print);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(nd_fespace.GetTrueVSize());
}

FluxProjector::FluxProjector(const MaterialOperator &mat_op,
                             const FiniteElementSpace &h1_fespace,
                             const FiniteElementSpace &h1d_fespace, double tol, int max_it,
                             int print)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  {
    // Flux operator is always partially assembled.
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm flux(h1_fespace, h1d_fespace);
    flux.AddDomainIntegrator<GradientIntegrator>(epsilon_func);
    Flux = std::make_unique<ParOperator>(flux.PartialAssemble(), h1_fespace, h1d_fespace,
                                         false);
  }
  M = GetMassMatrix(h1_fespace);

  ksp = ConfigureLinearSolver(h1_fespace.GetComm(), tol, max_it, print);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(h1d_fespace.GetTrueVSize());
}

template <typename VecType>
void FluxProjector::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::SOLVEESTIMATOR);
  MFEM_ASSERT(y.Size() == rhs.Size(), "Invalid vector dimensions for FluxProjector::Mult!");
  MFEM_ASSERT(
      y.Size() % x.Size() == 0,
      "Invalid vector dimension for FluxProjector::Mult, does not yield even blocking!");
  auto MultImpl = [this](const Vector &x_, Vector &y_)
  {
    const int vdim = y_.Size() / x_.Size();
    Flux->Mult(x_, rhs);
    if (vdim == 1)
    {
      // Mpi::Print(" Computing smooth flux projection for error estimation\n");
      ksp->Mult(rhs, y_);
    }
    else
    {
      for (int i = 0; i < vdim; i++)
      {
        // Mpi::Print(" Computing smooth flux projection of flux component {:d}/{:d} for "
        //            "error estimation\n",
        //            i + 1, vdim);
        const Vector rhsb(rhs, i * x_.Size(), x_.Size());
        Vector yb(y_, i * x_.Size(), x_.Size());
        ksp->Mult(rhsb, yb);
      }
    }
  };
  if constexpr (std::is_same<VecType, ComplexVector>::value)
  {
    MultImpl(x.Real(), y.Real());
    MultImpl(x.Imag(), y.Imag());
  }
  else
  {
    MultImpl(x, y);
  }
}

template <typename VecType>
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(const MaterialOperator &mat_op,
                                                        FiniteElementSpace &nd_fespace,
                                                        double tol, int max_it, int print)
  : mat_op(mat_op), nd_fespace(nd_fespace),
    projector(mat_op, nd_fespace, tol, max_it, print), F(nd_fespace.GetTrueVSize()),
    F_gf(&nd_fespace.Get()), U_gf(&nd_fespace.Get())
{
}

template <typename VecType>
ErrorIndicator CurlFluxErrorEstimator<VecType>::ComputeIndicators(const VecType &U) const
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(U, F);
  if constexpr (std::is_same<VecType, ComplexVector>::value)
  {
    F_gf.real().SetFromTrueDofs(F.Real());
    F_gf.imag().SetFromTrueDofs(F.Imag());
    U_gf.real().SetFromTrueDofs(U.Real());
    U_gf.imag().SetFromTrueDofs(U.Imag());
  }
  else
  {
    F_gf.SetFromTrueDofs(F);
    U_gf.SetFromTrueDofs(U);
  }

  // Loop over elements and accumulate the estimates from this component. The discontinuous
  // flux is μ⁻¹ ∇ × U.
  auto &mesh = nd_fespace.GetParMesh();
  Vector estimates(mesh.GetNE());
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

        auto AccumulateError =
            [&](const mfem::ParGridFunction &U_gf_, const mfem::ParGridFunction &F_gf_)
        {
          // μ⁻¹ ∇ × U
          U_gf_.GetSubVector(dofs, loc_gf);
          if (dof_trans.GetDofTransformation())
          {
            dof_trans.InvTransformPrimal(loc_gf);
          }
          Curl.MultTranspose(loc_gf, V_ip);
          T.Jacobian().Mult(V_ip, V_smooth);
          mat_op.GetInvPermeability(T.Attribute).Mult(V_smooth, V_ip);
          V_ip *= 1.0 / T.Weight();

          // Smooth flux
          F_gf_.GetSubVector(dofs, loc_gf);
          if (dof_trans.GetDofTransformation())
          {
            dof_trans.InvTransformPrimal(loc_gf);
          }
          Interp.MultTranspose(loc_gf, V_tmp);
          T.InverseJacobian().MultTranspose(V_tmp, V_smooth);

          V_smooth -= V_ip;
          elem_err += w * (V_smooth * V_smooth);
          loc_norm2 += w * (V_ip * V_ip);
        };
        if constexpr (std::is_same<VecType, ComplexVector>::value)
        {
          AccumulateError(U_gf.real(), F_gf.real());
          AccumulateError(U_gf.imag(), F_gf.imag());
        }
        else
        {
          AccumulateError(U_gf, F_gf);
        }
      }
      estimates[e] = std::sqrt(elem_err);
    }
    norm2 += loc_norm2;
  }

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &norm2, mesh.GetComm());
  if (norm2 > 0.0)
  {
    estimates *= 1.0 / std::sqrt(norm2);
  }
  return ErrorIndicator(std::move(estimates));
}

GradFluxErrorEstimator::GradFluxErrorEstimator(const MaterialOperator &mat_op,
                                               FiniteElementSpace &h1_fespace, double tol,
                                               int max_it, int print)
  : mat_op(mat_op), h1_fespace(h1_fespace),
    h1d_fespace(std::make_unique<FiniteElementSpace>(
        h1_fespace.GetMesh(), &h1_fespace.GetFEColl(), h1_fespace.SpaceDimension(),
        mfem::Ordering::byNODES)),
    projector(mat_op, h1_fespace, *h1d_fespace, tol, max_it, print),
    F(h1d_fespace->GetTrueVSize()), F_gf(&h1d_fespace->Get()), U_gf(&h1_fespace.Get())
{
}

ErrorIndicator GradFluxErrorEstimator::ComputeIndicators(const Vector &U) const
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(U, F);
  F_gf.SetFromTrueDofs(F);
  U_gf.SetFromTrueDofs(U);

  // Loop over elements and accumulate the estimates from this component. The discontinuous
  // flux is ε ∇U.
  auto &mesh = h1_fespace.GetParMesh();
  Vector estimates(mesh.GetNE());
  double norm2 = 0.0;
  PalacePragmaOmp(parallel reduction(+ : norm2))
  {
    // Assuming dim == space_dim
    mfem::IsoparametricTransformation T;
    mfem::Array<int> dofs, vdofs;
    mfem::Vector V_ip(h1d_fespace->GetVDim()), V_smooth(h1d_fespace->GetVDim()), loc_gf;
    mfem::Vector Interp;
    mfem::DenseMatrix Grad;

    double loc_norm2 = 0.0;
    PalacePragmaOmp(for schedule(static))
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const mfem::FiniteElement &fe = *h1d_fespace->Get().GetFE(e);
      mesh.GetElementTransformation(e, &T);
      h1_fespace.Get().GetElementDofs(e, dofs);
      vdofs = dofs;
      h1d_fespace->Get().DofsToVDofs(vdofs);
      Interp.SetSize(fe.GetDof());
      Grad.SetSize(fe.GetDof(), V_ip.Size());
      const int q_order = fem::DefaultIntegrationOrder::Get(T);
      const mfem::IntegrationRule &ir =
          mfem::IntRules.Get(mesh.GetElementGeometry(e), q_order);

      double elem_err = 0.0;
      for (int i = 0; i < ir.GetNPoints(); i++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(i);
        T.SetIntPoint(&ip);
        fe.CalcShape(ip, Interp);
        fe.CalcDShape(ip, Grad);
        const double w = ip.weight * T.Weight();

        // ε ∇U
        U_gf.GetSubVector(dofs, loc_gf);
        Grad.MultTranspose(loc_gf, V_ip);
        T.InverseJacobian().MultTranspose(V_ip, V_smooth);
        mat_op.GetPermittivityReal(T.Attribute).Mult(V_smooth, V_ip);

        // Smooth flux
        F_gf.GetSubVector(vdofs, loc_gf);
        for (int k = 0; k < h1d_fespace->GetVDim(); k++)
        {
          V_smooth(k) = Interp * (&loc_gf(Interp.Size() * k));
        }

        V_smooth -= V_ip;
        elem_err += w * (V_smooth * V_smooth);
        loc_norm2 += w * (V_ip * V_ip);
      }
      estimates[e] = std::sqrt(elem_err);
    }
    norm2 += loc_norm2;
  }

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &norm2, mesh.GetComm());
  if (norm2 > 0.0)
  {
    estimates *= 1.0 / std::sqrt(norm2);
  }
  return ErrorIndicator(std::move(estimates));
}

template void FluxProjector::Mult(const Vector &, Vector &) const;
template void FluxProjector::Mult(const ComplexVector &, ComplexVector &) const;

template class CurlFluxErrorEstimator<Vector>;
template class CurlFluxErrorEstimator<ComplexVector>;

}  // namespace palace
