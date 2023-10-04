// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <limits>
#include "errorestimator.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/errorindicators.hpp"
#include "utils/iodata.hpp"
#include "linalg/amg.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/rap.hpp"
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

  constexpr bool ScalarFESpace = std::is_same<SmoothFluxFiniteElementCollection, mfem::H1_FECollection>::value
    || std::is_same<SmoothFluxFiniteElementCollection, mfem::L2_FECollection>::value;

  // Assemble the bilinear form operator
  auto M = std::make_unique<MultigridOperator>(h.GetNumLevels());
  for (int l = 0; l < h.GetNumLevels(); l++)
  {
    auto &h_l = h.GetFESpaceAtLevel(l);
    auto m = std::make_unique<mfem::SymmetricBilinearForm>(&h_l);

    if constexpr (ScalarFESpace)
    {
      MFEM_ASSERT(h_l.GetVDim() == 1,
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
FluxProjector<SmoothFluxFiniteElementCollection>::FluxProjector(mfem::ParFiniteElementSpaceHierarchy &smooth_flux_fespace, double tol, int max_it, int print, int pa_order_threshold)
  : M(BuildMassMatrixOperator<SmoothFluxFiniteElementCollection>(smooth_flux_fespace, pa_order_threshold))
{
  // The system matrix for the projection is real and SPD. For the coarse-level AMG solve,
  // we don't use an exact solve on the coarsest level.
  auto amg =
      std::make_unique<WrapperSolver<Operator>>(std::make_unique<BoomerAmgSolver>(1, 1, 0));
  auto gmg = std::make_unique<GeometricMultigridSolver<Operator>>(
      std::move(amg), smooth_flux_fespace, nullptr, 1, 1, 2, 1.0, 0.0, true,
      pa_order_threshold);

  auto pcg = std::make_unique<CgSolver<Operator>>(
      smooth_flux_fespace.GetFinestFESpace().GetComm(), print);

  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);

  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(gmg));
  ksp->SetOperators(*M, *M);

  tmp.SetSize(smooth_flux_fespace.GetFinestFESpace().GetTrueVSize());
}

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    mfem::ParFiniteElementSpaceHierarchy &nd_fespaces)
  : mat_op(mat_op), nd_fespaces(nd_fespaces),
    smooth_projector(nd_fespaces, iodata.solver.linear.tol, 200, 0,
                     iodata.solver.pa_order_threshold),
    smooth_flux(nd_fespaces.GetFinestFESpace().GetTrueVSize()),
    flux_rhs(nd_fespaces.GetFinestFESpace().GetTrueVSize())
{
}

template <>
IndicatorsAndNormalization CurlFluxErrorEstimator::ComputeIndicators(
  const ComplexVector& v, bool normalize)
{
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  mfem::ParComplexGridFunction field(&nd_fespace);
  field.real().SetFromTrueDofs(v.Real());
  field.imag().SetFromTrueDofs(v.Imag());
  const int nelem = nd_fespace.GetNE();

  Vector smooth_vec, coarse_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;

  mfem::ParGridFunction flux_func(&nd_fespaces.GetFinestFESpace());
  for (bool real : {true, false})
  {
    auto &field_component = real ? field.real() : field.imag();

    // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
    CurlFluxCoefficient coef(field_component, mat_op);
    {
      mfem::LinearForm rhs(&nd_fespace);
      rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
      rhs.UseFastAssembly(false);
      rhs.Assemble();
      nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, flux_rhs);
    }

    // Given the RHS vector of non-smooth flux, construct a flux projector and perform mass
    // matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
    smooth_projector.Mult(flux_rhs, smooth_flux);
    flux_func.SetFromTrueDofs(smooth_flux);
    flux_func.ExchangeFaceNbrData();

    // Loop over elements and accumulate the estimates from this component
    for (int e = 0; e < nd_fespace.GetNE(); e++)
    {
      auto &T = *nd_fespace.GetElementTransformation(e);
      // integration order 2p + q
      const auto &ir = mfem::IntRules.Get(T.GetGeometryType(), 2 * nd_fespace.GetFE(e)->GetOrder() + T.Order());
      for (const auto &ip : ir)
      {
        T.SetIntPoint(&ip);

        flux_func.GetVectorValue(e, ip, smooth_vec);
        coef.Eval(coarse_vec, T, ip);
        double w_i = ip.weight * T.Weight();
        for (int c = 0; c < 3; c++)
        {
          estimates[e] += w_i * std::pow(smooth_vec[c] - coarse_vec[c], 2.0);
          normalization += w_i * coarse_vec[c] * coarse_vec[c];
        }
      }
    }
  }
  for (auto &e : estimates)
  {
    e = std::sqrt(e);
  }

  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);
  if (normalize)
  {
    std::for_each(estimates.begin(), estimates.end(),
                  [&normalization](auto &x) { x /= normalization; });
  }
  return {estimates, normalization};
}

template <>
IndicatorsAndNormalization
CurlFluxErrorEstimator::ComputeIndicators(const Vector &v, bool normalize)
{
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  mfem::ParGridFunction field(&nd_fespace);
  field.SetFromTrueDofs(v);
  const int nelem = nd_fespaces.GetFinestFESpace().GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient coef(field, mat_op);
  {
    mfem::LinearForm rhs(&nd_fespace);
    rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, flux_rhs);
  }

  // Given the RHS vector of non-smooth flux, construct a flux projector and perform mass
  // matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  smooth_projector.Mult(flux_rhs, smooth_flux);

  // Given a complex solution represented with a ComplexVector, build a ComplexGridFunction
  // for evaluation.
  auto build_func = [](const Vector &f, mfem::ParFiniteElementSpace &fespace)
  {
    mfem::ParGridFunction flux(&fespace);
    flux.SetFromTrueDofs(f);
    flux.ExchangeFaceNbrData();
    return flux;
  };
  auto smooth_flux_func = build_func(smooth_flux, nd_fespaces.GetFinestFESpace());

  Vector smooth_vec, coarse_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;
  for (int e = 0; e < nd_fespace.GetNE(); e++)
  {
    auto &T = *nd_fespace.GetElementTransformation(e);
    // integration order 2p + q
    const auto &ir = mfem::IntRules.Get(T.GetGeometryType(), 2 * nd_fespace.GetFE(e)->GetOrder() + T.Order());
    for (const auto &ip : ir)
    {
      T.SetIntPoint(&ip);
      smooth_flux_func.GetVectorValue(e, ip, smooth_vec);
      coef.Eval(coarse_vec, T, ip);
      for (int c = 0; c < 3; c++)
      {
        estimates[e] += ip.weight * T.Weight() * std::pow(smooth_vec[c] - coarse_vec[c], 2.0);
        normalization += ip.weight * T.Weight() * coarse_vec[c] * coarse_vec[c];
      }
    }
    estimates[e] = std::sqrt(estimates[e]);
  }

  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);
  if (normalize)
  {
    std::for_each(estimates.begin(), estimates.end(),
                  [&normalization](auto &x) { x /= normalization; });
  }
  return {estimates, normalization};
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fespace)
  : mat_op(mat_op), fespace(fespace),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    smooth_flux_component_fespace(ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        iodata.solver.pa_order_threshold, mesh, smooth_flux_fecs)),
    smooth_flux_fespace(mesh.back().get(), smooth_flux_fecs.back().get(),
                    mesh.back()->Dimension()),
    smooth_projector(smooth_flux_component_fespace, iodata.solver.linear.tol, 200, 0,
                     iodata.solver.pa_order_threshold),
    coarse_flux_fec(iodata.solver.order, mesh.back()->Dimension(),
                    mfem::BasisType::GaussLobatto),
    coarse_flux_fespace(mesh.back().get(), &coarse_flux_fec, mesh.back()->Dimension()),
    scalar_mass_matrices(fespace.GetNE()), smooth_to_coarse_embed(fespace.GetNE())
{
  mfem::MassIntegrator mass_integrator;

  for (int e = 0; e < fespace.GetNE(); e++)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto &coarse_fe = *coarse_flux_fespace.GetFE(e);
    auto &T = *fespace.GetElementTransformation(e);
    mass_integrator.AssembleElementMatrix(coarse_fe, T, scalar_mass_matrices[e]);

    const auto &smooth_fe = *smooth_flux_component_fespace.GetFinestFESpace().GetFE(e);
    coarse_fe.Project(smooth_fe, T, smooth_to_coarse_embed[e]);
  }
}

IndicatorsAndNormalization GradFluxErrorEstimator::ComputeIndicators(const Vector &v,
                                                              bool normalize) const
{
  const int sdim = fespace.GetMesh()->SpaceDimension();
  mfem::ParGridFunction field(&fespace);
  field.SetFromTrueDofs(v);
  const int nelem = smooth_flux_fespace.GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (V, ϵ ∇ ϕ).
  GradFluxCoefficient coef(field, mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &fespace, auto &coef)
  {
    Vector RHS(fespace.GetTrueVSize());

    mfem::LinearForm rhs(&fespace);
    rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    fespace.GetProlongationMatrix()->MultTranspose(rhs, RHS);

    return RHS;
  };
  auto smooth_flux_rhs = rhs_from_coef(smooth_flux_fespace, coef);

  // Given the RHS vector of non-smooth flux, construct a flux projector and perform
  // component wise mass matrix inversion in the appropriate space, giving fᵢ = M⁻¹ f̂ᵢ.
  auto build_flux = [sdim](const FluxProjector<mfem::H1_FECollection> &proj, Vector &rhs)
  {
    // Use a copy construction to match appropriate size.
    Vector flux(rhs.Size());
    flux = 0.0;

    // Apply the flux projector component wise.
    const int ndof = flux.Size();
    const int stride = ndof / sdim;
    MFEM_ASSERT(ndof % 3 == 0, "!");

    Vector flux_comp, rhs_comp;
    for (int i = 0; i < 3; i++)
    {
      flux_comp.MakeRef(flux, i * stride, stride);
      rhs_comp.MakeRef(rhs, i * stride, stride);
      proj.Mult(rhs_comp, flux_comp);
    }
    return flux;
  };
  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);

  // Given a solution represented with a Vector, build a GridFunction for evaluation.
  auto build_func = [](const Vector &f, mfem::ParFiniteElementSpace &fespace)
  {
    mfem::ParGridFunction flux(&fespace);
    flux.SetFromTrueDofs(f);
    flux.ExchangeFaceNbrData();
    return flux;
  };
  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fespace);

  mfem::ParGridFunction coarse_flux(&coarse_flux_fespace);
  coarse_flux.ProjectCoefficient(coef);

  Vector coarse_vec, smooth_vec, coarse_sub_vec, smooth_sub_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;
  for (int e = 0; e < fespace.GetNE(); e++)
  {
    coarse_flux.GetElementDofValues(e, coarse_vec);
    smooth_flux_func.GetElementDofValues(e, smooth_vec);

    MFEM_ASSERT(coarse_vec.Size() == smooth_vec.Size() && coarse_vec.Size() % 3 == 0,
                "H1d and L2d spaces should have the same number of elemental dofs, "
                "and should be exactly divisible by 3: "
                    << coarse_vec.Size() << " " << smooth_vec.Size());

    const int ndof = coarse_vec.Size() / sdim;
    coarse_sub_vec.SetSize(ndof);
    smooth_sub_vec.SetSize(ndof);
    for (int c = 0; c < 3; c++)
    {
      coarse_sub_vec.MakeRef(coarse_vec, c * ndof);
      smooth_sub_vec.MakeRef(smooth_vec, c * ndof);

      normalization += scalar_mass_matrices[e].InnerProduct(coarse_sub_vec, coarse_sub_vec);
      smooth_to_coarse_embed[e].AddMult(smooth_sub_vec, coarse_sub_vec, -1.0);  // Embed
      estimates[e] += scalar_mass_matrices[e].InnerProduct(coarse_sub_vec,
                                                           coarse_sub_vec);  // Integrate
    }
    estimates[e] = std::sqrt(estimates[e]);
  }
  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);

  if constexpr (false)
  {
    // Debugging branch generates some intermediate fields for paraview.
    mfem::ParaViewDataCollection paraview("debug", fespace.GetParMesh());
    paraview.RegisterVCoeffField("Flux", &coef);

    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fespace);
    paraview.RegisterField("SmoothFlux", &smooth_flux_func);

    mfem::L2_FECollection est_fec(0, 3);
    mfem::ParFiniteElementSpace est_fespace(fespace.GetParMesh(), &est_fec);
    mfem::ParGridFunction est_field(&est_fespace);
    est_field.SetFromTrueDofs(estimates);

    paraview.RegisterField("ErrorIndicator", &est_field);
    paraview.Save();
  }
  if (normalize)
  {
    std::for_each(estimates.begin(), estimates.end(),
                  [&normalization](auto &x) { x /= normalization; });
  }
  return {estimates, normalization};
}


}  // namespace palace
