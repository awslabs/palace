// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <limits>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

std::unique_ptr<Operator> GetMassMatrix(mfem::ParFiniteElementSpaceHierarchy &fespaces,
                                        int pa_order_threshold, int skip_zeros)
{
  const int dim = fespaces.GetFinestFESpace().GetParMesh()->Dimension();
  const auto type = fespaces.GetFinestFESpace().FEColl()->GetRangeType(dim);
  auto M = std::make_unique<MultigridOperator>(fespaces.GetNumLevels());
  for (int l = 0; l < fespaces.GetNumLevels(); l++)
  {
    // Force coarse level operator to be fully assembled always.
    auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
    auto m = std::make_unique<mfem::SymmetricBilinearForm>(&fespace_l);
    if (type == mfem::FiniteElement::SCALAR)
    {
      MFEM_VERIFY(fespace_l.GetVDim() == 1,
                  "Scalar mass matrix hierarchy assumes a component-wise solve.");
      m->AddDomainIntegrator(new mfem::MassIntegrator);
    }
    else
    {
      m->AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
    }
    auto M_l = std::make_unique<ParOperator>(
        fem::AssembleOperator(std::move(m), true, (l > 0) ? pa_order_threshold : 100,
                              skip_zeros),
        fespace_l);
    M->AddOperator(std::move(M_l));
  }
  return M;
}

std::unique_ptr<KspSolver>
ConfigureLinearSolver(mfem::ParFiniteElementSpaceHierarchy &fespaces, double tol,
                      int max_it, int print, int pa_order_threshold)
{
  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg =
      std::make_unique<WrapperSolver<Operator>>(std::make_unique<BoomerAmgSolver>(1, 1, 0));
  std::unique_ptr<Solver<Operator>> pc;
  if (fespaces.GetNumLevels() > 1)
  {
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        std::move(amg), fespaces, nullptr, 1, 1, 2, 1.0, 0.0, true, pa_order_threshold);
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
                             mfem::ParFiniteElementSpaceHierarchy &nd_fespaces, double tol,
                             int max_it, int print, int pa_order_threshold)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  constexpr int skip_zeros = 0;
  {
    // XX TODO: No partial assembly is available yet for this operator.
    constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
    MaterialPropertyCoefficient<MatType> muinv_func(mat_op);
    auto flux = std::make_unique<mfem::MixedBilinearForm>(&nd_fespaces.GetFinestFESpace(),
                                                          &nd_fespaces.GetFinestFESpace());
    flux->AddDomainIntegrator(new mfem::MixedVectorCurlIntegrator(muinv_func));
    flux->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    flux->Assemble(skip_zeros);
    flux->Finalize(skip_zeros);
    Flux = std::make_unique<ParOperator>(std::move(flux), nd_fespaces.GetFinestFESpace(),
                                         nd_fespaces.GetFinestFESpace(), false);
  }
  M = GetMassMatrix(nd_fespaces, pa_order_threshold, skip_zeros);

  ksp = ConfigureLinearSolver(nd_fespaces, tol, max_it, print, pa_order_threshold);
  ksp->SetOperators(*M, *M);

  rhs.SetSize(nd_fespaces.GetFinestFESpace().GetTrueVSize());
}

FluxProjector::FluxProjector(const MaterialOperator &mat_op,
                             mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
                             mfem::ParFiniteElementSpace &h1d_fespace, double tol,
                             int max_it, int print, int pa_order_threshold)
{
  BlockTimer bt(Timer::CONSTRUCTESTIMATOR);
  constexpr int skip_zeros = 0;
  {
    // XX TODO: Matrix coefficient support for GradientIntegrator.
    // XX TODO: No partial assembly is available yet for this operator.
    //  constexpr auto MatType = MaterialPropertyType::PERMITTIVITY_REAL;
    //  MaterialPropertyCoefficient<MatType> epsilon_func(mat_op);
    auto flux = std::make_unique<mfem::MixedBilinearForm>(&h1_fespaces.GetFinestFESpace(),
                                                          &h1d_fespace);
    flux->AddDomainIntegrator(new mfem::GradientIntegrator);
    flux->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
    flux->Assemble(skip_zeros);
    flux->Finalize(skip_zeros);
    Flux = std::make_unique<ParOperator>(std::move(flux), h1_fespaces.GetFinestFESpace(),
                                         h1d_fespace, false);
  }
  M = GetMassMatrix(h1_fespaces, pa_order_threshold, skip_zeros);

  ksp = ConfigureLinearSolver(h1_fespaces, tol, max_it, print, pa_order_threshold);
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
    const MaterialOperator &mat_op, mfem::ParFiniteElementSpaceHierarchy &nd_fespaces,
    double tol, int max_it, int print_level, int pa_order_threshold)
  : mat_op(mat_op), nd_fespaces(nd_fespaces),
    projector(mat_op, nd_fespaces, tol, max_it, print_level, pa_order_threshold),
    F(nd_fespaces.GetFinestFESpace().GetTrueVSize()), F_gf(&nd_fespaces.GetFinestFESpace()),
    U_gf(&nd_fespaces.GetFinestFESpace())
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
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  auto &mesh = *nd_fespace.GetParMesh();
  Vector estimates(mesh.GetNE()), V_ip(mesh.SpaceDimension()), V_smooth;
  double normalization = 0.0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const mfem::FiniteElement &fe = *nd_fespace.GetFE(e);
    mfem::ElementTransformation &T = *mesh.GetElementTransformation(e);
    const int q_order = 2 * fe.GetOrder() + T.OrderW();
    const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), q_order);
    double err = 0.0;
    for (int i = 0; i < ir.GetNPoints(); i++)
    {
      const mfem::IntegrationPoint &ip = ir.IntPoint(i);
      T.SetIntPoint(&ip);
      const double w = ip.weight * T.Weight();
      if constexpr (std::is_same<VecType, ComplexVector>::value)
      {
        // Real part
        U_gf.real().GetCurl(T, V_smooth);
        mat_op.GetInvPermeability(T.Attribute).Mult(V_smooth, V_ip);
        F_gf.real().GetVectorValue(T, ip, V_smooth);
        V_smooth -= V_ip;
        err += w * (V_smooth * V_smooth);
        normalization += w * (V_ip * V_ip);

        // Imaginary part
        U_gf.imag().GetCurl(T, V_smooth);
        mat_op.GetInvPermeability(T.Attribute).Mult(V_smooth, V_ip);
        F_gf.imag().GetVectorValue(T, ip, V_smooth);
        V_smooth -= V_ip;
        err += w * (V_smooth * V_smooth);
        normalization += w * (V_ip * V_ip);
      }
      else
      {
        U_gf.GetCurl(T, V_smooth);
        mat_op.GetInvPermeability(T.Attribute).Mult(V_smooth, V_ip);
        F_gf.GetVectorValue(T, ip, V_smooth);
        V_smooth -= V_ip;
        err += w * (V_smooth * V_smooth);
        normalization += w * (V_ip * V_ip);
      }
    }
    estimates[e] = std::sqrt(err);
  }

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &normalization, mesh.GetComm());
  normalization = std::sqrt(normalization);
  if (normalization > 0.0)
  {
    estimates *= 1.0 / normalization;
  }
  return ErrorIndicator(std::move(estimates), normalization);
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const MaterialOperator &mat_op, mfem::ParFiniteElementSpaceHierarchy &h1_fespaces,
    double tol, int max_it, int print_level, int pa_order_threshold)
  : mat_op(mat_op), h1_fespaces(h1_fespaces),
    h1d_fespace(h1_fespaces.GetFinestFESpace().GetParMesh(),
                h1_fespaces.GetFinestFESpace().FEColl(),
                h1_fespaces.GetFinestFESpace().GetParMesh()->SpaceDimension(),
                mfem::Ordering::byNODES),
    projector(mat_op, h1_fespaces, h1d_fespace, tol, max_it, print_level,
              pa_order_threshold),
    F(h1d_fespace.GetTrueVSize()), F_gf(&h1d_fespace), U_gf(&h1_fespaces.GetFinestFESpace())
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
  auto &h1_fespace = h1_fespaces.GetFinestFESpace();
  auto &mesh = *h1_fespace.GetParMesh();
  Vector estimates(mesh.GetNE()), V_ip(mesh.SpaceDimension()), V_smooth;
  double normalization = 0.0;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const mfem::FiniteElement &fe = *h1d_fespace.GetFE(e);
    mfem::ElementTransformation &T = *mesh.GetElementTransformation(e);
    const int q_order = 2 * fe.GetOrder() + T.OrderW();
    const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), q_order);
    double err = 0.0;
    for (int i = 0; i < ir.GetNPoints(); i++)
    {
      // XX TODO: For now the flux is just ∇U since the coefficient support for
      // matrix-valued
      //          permittivity is not yet there (coming soon).
      const mfem::IntegrationPoint &ip = ir.IntPoint(i);
      T.SetIntPoint(&ip);
      const double w = ip.weight * T.Weight();
      U_gf.GetGradient(T, V_smooth);
      // mat_op.GetPermittivityReal(T.Attribute).Mult(V_smooth, V_ip);
      V_ip = V_smooth;
      F_gf.GetVectorValue(T, ip, V_smooth);
      V_smooth -= V_ip;
      err += w * (V_smooth * V_smooth);
      normalization += w * (V_ip * V_ip);
    }
    estimates[e] = std::sqrt(err);
  }

  // Finalize the element-wise error estimates.
  Mpi::GlobalSum(1, &normalization, mesh.GetComm());
  normalization = std::sqrt(normalization);
  if (normalization > 0.0)
  {
    estimates *= 1.0 / normalization;
  }
  return ErrorIndicator(std::move(estimates), normalization);
}

template void FluxProjector::Mult(const Vector &, Vector &) const;
template void FluxProjector::Mult(const ComplexVector &, ComplexVector &) const;

template class CurlFluxErrorEstimator<Vector>;
template class CurlFluxErrorEstimator<ComplexVector>;

}  // namespace palace
