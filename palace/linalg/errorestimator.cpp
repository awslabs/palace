// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include <limits>
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace fem;

// Given a finite element space hierarchy, construct a vector of mass matrix
// operators corresponding to each level.
template <typename SmoothFluxFiniteElementCollection>
std::unique_ptr<Operator> BuildMassMatrixOperator(mfem::ParFiniteElementSpaceHierarchy &h,
                                                  int pa_order_threshold)
{
  constexpr int skip_zeros = 0;

  constexpr bool ScalarFESpace =
      std::is_same<SmoothFluxFiniteElementCollection, mfem::H1_FECollection>::value ||
      std::is_same<SmoothFluxFiniteElementCollection, mfem::L2_FECollection>::value;

  // Assemble the bilinear form operator
  auto M = std::make_unique<MultigridOperator>(h.GetNumLevels());
  for (int l = 0; l < h.GetNumLevels(); l++)
  {
    auto &h_l = h.GetFESpaceAtLevel(l);
    auto m = std::make_unique<mfem::SymmetricBilinearForm>(&h_l);

    if constexpr (ScalarFESpace)
    {
      MFEM_VERIFY(h_l.GetVDim() == 1,
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
        h_l);

    // Set the essential dofs (none).
    M->AddOperator(std::move(M_l));
  }
  return M;
}

template <typename SmoothFluxFiniteElementCollection>
FluxProjector<SmoothFluxFiniteElementCollection>::FluxProjector(
    mfem::ParFiniteElementSpaceHierarchy &fespaces, double tol, int max_it, int print,
    int pa_order_threshold)
  : M(BuildMassMatrixOperator<SmoothFluxFiniteElementCollection>(fespaces,
                                                                 pa_order_threshold))
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

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));
  ksp->SetOperators(*M, *M);
}

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    mfem::ParFiniteElementSpaceHierarchy &nd_fespaces)
  : mat_op(mat_op), nd_fespaces(nd_fespaces),
    smooth_projector(nd_fespaces, iodata.solver.linear.estimator_tol,
                     iodata.solver.linear.estimator_max_it, iodata.problem.verbose,
                     iodata.solver.pa_order_threshold),
    smooth_flux(nd_fespaces.GetFinestFESpace().GetTrueVSize()),
    flux_rhs(nd_fespaces.GetFinestFESpace().GetTrueVSize()),
    field_gf(&nd_fespaces.GetFinestFESpace()),
    smooth_flux_gf(&nd_fespaces.GetFinestFESpace())
{
}

template <>
ErrorIndicator CurlFluxErrorEstimator::ComputeIndicators(const ComplexVector &v) const
{
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  const int nelem = nd_fespace.GetNE();

  Vector smooth_vec, coarse_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;
  for (bool real : {true, false})
  {
    field_gf.SetFromTrueDofs(real ? v.Real() : v.Imag());

    // Coefficients for computing the discontinuous flux component, i.e. (W, μ⁻¹∇ × V).
    CurlFluxCoefficient coef(field_gf, mat_op);
    {
      // Given the RHS vector of non-smooth flux, construct a flux projector and perform
      // mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
      Mpi::Print("Computing smooth flux approximation of {} component\n",
                 real ? "real" : "imaginary");
      BlockTimer bt(Timer::ESTSOLVE);
      mfem::LinearForm rhs(&nd_fespace);
      rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
      rhs.UseFastAssembly(false);
      rhs.Assemble();
      nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, flux_rhs);
      smooth_projector.Mult(flux_rhs, smooth_flux);
    }
    smooth_flux_gf.SetFromTrueDofs(smooth_flux);
    smooth_flux_gf.ExchangeFaceNbrData();

    // Loop over elements and accumulate the estimates from this component
    for (int e = 0; e < nd_fespace.GetNE(); e++)
    {
      auto &T = *nd_fespace.GetElementTransformation(e);
      // integration order 2p + q
      const auto &ir = mfem::IntRules.Get(T.GetGeometryType(),
                                          2 * nd_fespace.GetFE(e)->GetOrder() + T.Order());
      for (const auto &ip : ir)
      {
        T.SetIntPoint(&ip);

        smooth_flux_gf.GetVectorValue(e, ip, smooth_vec);
        coef.Eval(coarse_vec, T, ip);
        double w_i = ip.weight * T.Weight();
        constexpr int vdim = 3;
        for (int c = 0; c < vdim; c++)
        {
          estimates[e] += w_i * std::pow(smooth_vec[c] - coarse_vec[c], 2.0);
          normalization += w_i * coarse_vec[c] * coarse_vec[c];
        }
      }
    }
  }
  linalg::Sqrt(estimates);

  Mpi::GlobalSum(1, &normalization, nd_fespace.GetComm());
  normalization = std::sqrt(normalization);
  if (normalization > 0)
  {
    estimates /= normalization;
  }
  return ErrorIndicator(std::move(estimates), normalization);
}

template <>
ErrorIndicator CurlFluxErrorEstimator::ComputeIndicators(const Vector &v) const
{
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  field_gf.SetFromTrueDofs(v);
  const int nelem = nd_fespace.GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient coef(field_gf, mat_op);
  {
    // Given the RHS vector of non-smooth flux, construct a flux projector and perform mass
    // matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
    Mpi::Print("Computing smooth flux approximation\n");
    BlockTimer bt(Timer::ESTSOLVE);
    mfem::LinearForm rhs(&nd_fespace);
    rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, flux_rhs);
    smooth_projector.Mult(flux_rhs, smooth_flux);
  }

  // Given a complex solution represented with a ComplexVector, build a ComplexGridFunction
  // for evaluation.
  smooth_flux_gf.SetFromTrueDofs(smooth_flux);
  smooth_flux_gf.ExchangeFaceNbrData();

  Vector smooth_vec, coarse_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;
  for (int e = 0; e < nd_fespace.GetNE(); e++)
  {
    auto &T = *nd_fespace.GetElementTransformation(e);
    // integration order 2p + q
    const auto &ir = mfem::IntRules.Get(T.GetGeometryType(),
                                        2 * nd_fespace.GetFE(e)->GetOrder() + T.Order());
    for (const auto &ip : ir)
    {
      T.SetIntPoint(&ip);
      smooth_flux_gf.GetVectorValue(e, ip, smooth_vec);
      coef.Eval(coarse_vec, T, ip);
      const double w_i = ip.weight * T.Weight();
      constexpr int vdim = 3;
      for (int c = 0; c < vdim; c++)
      {
        estimates[e] += w_i * std::pow(smooth_vec[c] - coarse_vec[c], 2.0);
        normalization += w_i * coarse_vec[c] * coarse_vec[c];
      }
    }
    estimates[e] = std::sqrt(estimates[e]);
  }

  Mpi::GlobalSum(1, &normalization, nd_fespace.GetComm());
  normalization = std::sqrt(normalization);
  if (normalization > 0)
  {
    estimates /= normalization;
  }
  return ErrorIndicator(std::move(estimates), normalization);
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    mfem::ParFiniteElementSpaceHierarchy &h1_fespaces)
  : mat_op(mat_op), h1_fespaces(h1_fespaces),
    smooth_projector(h1_fespaces, iodata.solver.linear.estimator_tol,
                     iodata.solver.linear.estimator_max_it, iodata.problem.verbose,
                     iodata.solver.pa_order_threshold),
    smooth_flux(h1_fespaces.GetFinestFESpace().GetTrueVSize()),
    flux_rhs(h1_fespaces.GetFinestFESpace().GetTrueVSize()),
    field_gf(&h1_fespaces.GetFinestFESpace()),
    smooth_flux_gf(&h1_fespaces.GetFinestFESpace())
{
}

ErrorIndicator GradFluxErrorEstimator::ComputeIndicators(const Vector &v) const
{
  auto &h1_fespace = h1_fespaces.GetFinestFESpace();
  const int sdim = h1_fespace.GetMesh()->SpaceDimension();
  field_gf.SetFromTrueDofs(v);
  const int nelem = h1_fespace.GetNE();

  Vector estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;

  // Coefficient for computing the discontinuous flux., i.e. (Vᵢ, (ϵ ∇ ϕ)ᵢ).
  GradFluxCoefficient coef(field_gf, mat_op);
  for (int c = 0; c < sdim; c++)
  {
    coef.SetComponent(c);
    {
      // Given the RHS vector of non-smooth flux component, compute fᵢ = M⁻¹ f̂ᵢ.
      Mpi::Print("Computing smooth flux approximation of component {}\n", c);
      BlockTimer bt0(Timer::ESTSOLVE);
      mfem::LinearForm rhs(&h1_fespace);
      rhs.AddDomainIntegrator(new mfem::DomainLFIntegrator(coef));
      rhs.UseFastAssembly(false);
      rhs.Assemble();
      h1_fespace.GetProlongationMatrix()->MultTranspose(rhs, flux_rhs);
      smooth_projector.Mult(flux_rhs, smooth_flux);
    }
    smooth_flux_gf.SetFromTrueDofs(smooth_flux);
    smooth_flux_gf.ExchangeFaceNbrData();
    for (int e = 0; e < h1_fespace.GetNE(); e++)
    {
      auto &T = *h1_fespace.GetElementTransformation(e);
      // integration order 2p + q
      const auto &ir = mfem::IntRules.Get(T.GetGeometryType(),
                                          2 * h1_fespace.GetFE(e)->GetOrder() + T.Order());
      for (const auto &ip : ir)
      {
        T.SetIntPoint(&ip);
        double smooth_val = smooth_flux_gf.GetValue(e, ip);
        double coarse_val = coef.Eval(T, ip);
        const double w_i = ip.weight * T.Weight();
        estimates[e] += w_i * std::pow(smooth_val - coarse_val, 2.0);
        normalization += w_i * std::pow(coarse_val, 2.0);
      }
    }
    if constexpr (false)
    {
      // Debugging branch generates some intermediate fields for paraview.
      mfem::ParaViewDataCollection paraview("debug_coeff" + std::to_string(c),
                                            h1_fespace.GetParMesh());
      paraview.RegisterCoeffField("Flux", &coef);
      paraview.RegisterField("SmoothFlux", &smooth_flux_gf);
      paraview.Save();
    }
  }
  linalg::Sqrt(estimates);

  Mpi::GlobalSum(1, &normalization, h1_fespace.GetComm());
  normalization = std::sqrt(normalization);

  if (normalization > 0)
  {
    estimates /= normalization;
  }
  return ErrorIndicator(std::move(estimates), normalization);
}

}  // namespace palace
