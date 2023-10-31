// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/omp.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

std::unique_ptr<Operator> GetMassMatrix(const FiniteElementSpaceHierarchy &fespaces,
                                        int pa_order_threshold)
{
  constexpr bool skip_zeros = false;
  const int dim = fespaces.GetFinestFESpace().GetParMesh()->Dimension();
  const auto type = fespaces.GetFinestFESpace().FEColl()->GetRangeType(dim);
  auto M = std::make_unique<MultigridOperator>(fespaces.GetNumLevels());
  for (std::size_t l = 0; l < fespaces.GetNumLevels(); l++)
  {
    // Force coarse level operator to be fully assembled always.
    const auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
    BilinearForm m(fespace_l);
    if (type == mfem::FiniteElement::SCALAR)
    {
      MFEM_VERIFY(fespace_l.GetVDim() == 1,
                  "Scalar mass matrix hierarchy assumes a component-wise solve.");
      m.AddDomainIntegrator<MassIntegrator>();
    }
    else
    {
      m.AddDomainIntegrator<VectorFEMassIntegrator>();
    }
    auto M_l = std::make_unique<ParOperator>(
        m.Assemble((l > 0) ? pa_order_threshold : 99, skip_zeros), fespace_l);
    M->AddOperator(std::move(M_l));
  }
  return M;
}

std::unique_ptr<KspSolver>
ConfigureLinearSolver(const FiniteElementSpaceHierarchy &fespaces, double tol, int max_it,
                      int print)
{
  // The system matrix for the projection is real and SPD.
  auto amg =
      std::make_unique<WrapperSolver<Operator>>(std::make_unique<BoomerAmgSolver>(1, 1, 0));
  std::unique_ptr<Solver<Operator>> pc;
  if (fespaces.GetNumLevels() > 1)
  {
    const int mg_smooth_order =
        std::max(fespaces.GetFinestFESpace().GetMaxElementOrder(), 2);
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        std::move(amg), fespaces.GetProlongationOperators(), nullptr, 1, 1, mg_smooth_order,
        1.0, 0.0, true);
  }
  else
  {
    pc = std::move(amg);
  }

  auto pcg =
      std::make_unique<CgSolver<Operator>>(fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  return std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
}

}  // namespace

FluxProjector::FluxProjector(const MaterialOperator &mat_op,
                             const FiniteElementSpaceHierarchy &nd_fespaces, double tol,
                             int max_it, int print, int pa_order_threshold)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  {
    // Flux operator is always partially assembled.
    constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
    MaterialPropertyCoefficient<MatType> muinv_func(mat_op);
    BilinearForm flux(nd_fespaces.GetFinestFESpace());
    flux.AddDomainIntegrator<MixedVectorCurlIntegrator>(muinv_func);
    Flux = std::make_unique<ParOperator>(flux.Assemble(), nd_fespaces.GetFinestFESpace());
  }
  M = GetMassMatrix(nd_fespaces, pa_order_threshold);

  ksp = ConfigureLinearSolver(nd_fespaces, tol, max_it, print);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(nd_fespaces.GetFinestFESpace().GetTrueVSize());
}

FluxProjector::FluxProjector(const MaterialOperator &mat_op,
                             const FiniteElementSpaceHierarchy &h1_fespaces,
                             const FiniteElementSpace &h1d_fespace, double tol, int max_it,
                             int print, int pa_order_threshold)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  {
    // Flux operator is always partially assembled.
    constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_REAL;
    MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
    BilinearForm flux(h1_fespaces.GetFinestFESpace(), h1d_fespace);
    flux.AddDomainIntegrator<GradientIntegrator>(epsilon_func);
    Flux = std::make_unique<ParOperator>(flux.Assemble(), h1_fespaces.GetFinestFESpace(),
                                         h1d_fespace, false);
  }
  M = GetMassMatrix(h1_fespaces, pa_order_threshold);

  ksp = ConfigureLinearSolver(h1_fespaces, tol, max_it, print);
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
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(
    const MaterialOperator &mat_op, const FiniteElementSpaceHierarchy &nd_fespaces,
    double tol, int max_it, int print, int pa_order_threshold)
  : mat_op(mat_op), nd_fespace(nd_fespaces.GetFinestFESpace()),
    projector(mat_op, nd_fespaces, tol, max_it, print, pa_order_threshold),
    F(nd_fespace.GetTrueVSize()), F_gf(const_cast<FiniteElementSpace *>(&nd_fespace)),
    U_gf(const_cast<FiniteElementSpace *>(&nd_fespace))
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
  auto &mesh = *nd_fespace.GetParMesh();
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
      const mfem::FiniteElement &fe = *nd_fespace.GetFE(e);
      mesh.GetElementTransformation(e, &T);
      nd_fespace.GetElementDofs(e, dofs, dof_trans);
      Interp.SetSize(fe.GetDof(), V_ip.Size());
      Curl.SetSize(fe.GetDof(), V_ip.Size());
      const int q_order = fem::GetDefaultIntegrationOrder(fe, fe, T);
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

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const MaterialOperator &mat_op, const FiniteElementSpaceHierarchy &h1_fespaces,
    double tol, int max_it, int print, int pa_order_threshold)
  : mat_op(mat_op), h1_fespace(h1_fespaces.GetFinestFESpace()),
    h1d_fespace(std::make_unique<FiniteElementSpace>(
        h1_fespace.GetParMesh(), h1_fespace.FEColl(),
        h1_fespace.GetParMesh()->SpaceDimension(), mfem::Ordering::byNODES)),
    projector(mat_op, h1_fespaces, *h1d_fespace, tol, max_it, print, pa_order_threshold),
    F(h1d_fespace->GetTrueVSize()), F_gf(h1d_fespace.get()),
    U_gf(const_cast<FiniteElementSpace *>(&h1_fespace))
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
  auto &mesh = *h1_fespace.GetParMesh();
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
      const mfem::FiniteElement &fe = *h1d_fespace->GetFE(e);
      mesh.GetElementTransformation(e, &T);
      h1_fespace.GetElementDofs(e, dofs);
      vdofs = dofs;
      h1d_fespace->DofsToVDofs(vdofs);
      Interp.SetSize(fe.GetDof());
      Grad.SetSize(fe.GetDof(), V_ip.Size());
      const int q_order = fem::GetDefaultIntegrationOrder(fe, fe, T);
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
