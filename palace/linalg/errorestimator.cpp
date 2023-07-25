// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/errorindicators.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace utils;

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)), // TODO: pc_lor?
    smooth_flux_fes(ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
        iodata.solver.pa_order_threshold, mesh, smooth_flux_fecs)),
    smooth_projector(smooth_flux_fes, iodata.solver.linear.tol, 200, 0, iodata.solver.pa_order_threshold),
    coarse_flux_fec(iodata.solver.order, mesh.back()->Dimension(),
                    mfem::BasisType::GaussLobatto),
    coarse_flux_fes(mesh.back().get(), &coarse_flux_fec, mesh.back()->Dimension()),
    scalar_mass_matrices(fes.GetNE()), smooth_to_coarse_embed(fes.GetNE())
{
  mfem::MassIntegrator mass_integrator;
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto &coarse_fe = *coarse_flux_fes.GetFE(e);
    auto &T = *coarse_flux_fes.GetElementTransformation(e);
    mass_integrator.AssembleElementMatrix(coarse_fe, T, scalar_mass_matrices[e]);

    const auto &smooth_fe = *smooth_flux_fes.GetFinestFESpace().GetFE(e);
    coarse_fe.Project(smooth_fe, T, smooth_to_coarse_embed[e]);
  }
}

Vector CurlFluxErrorEstimator::operator()(const ComplexVector &v) const
{
  Timer local_timer;
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(v.Real());
  field.imag().SetFromTrueDofs(v.Imag());

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &flux_fes, auto &coef)
  {
    Vector RHS(flux_fes.GetTrueVSize());

    mfem::LinearForm rhs(&flux_fes);
    rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    flux_fes.GetProlongationMatrix()->MultTranspose(rhs, RHS);

    return RHS;
  };

  const auto smooth_flux_rhs =
      ComplexVector(rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
                    rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), imag_coef));

  local_timer.construct_time += local_timer.Lap();

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_flux = [](const FluxProjector &proj, const ComplexVector &flux_coef)
  {
    // Use a copy construction to match appropriate size.
    ComplexVector flux(flux_coef);
    proj.Mult(flux_coef, flux);
    return flux;
  };
  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);

  local_timer.solve_time += local_timer.Lap();

  // Given a complex solution represented with a PetscParVector, build a
  // ComplexGridFunction for evaluation.
  auto build_func = [](const ComplexVector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParComplexGridFunction flux(&fes);
    flux.real().SetFromTrueDofs(f.Real());
    flux.imag().SetFromTrueDofs(f.Imag());
    flux.real().ExchangeFaceNbrData();
    flux.imag().ExchangeFaceNbrData();
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

  mfem::ParComplexGridFunction coarse_flux_func(&coarse_flux_fes);

  local_timer.construct_time += local_timer.Lap();
  coarse_flux_func.real().ProjectCoefficient(real_coef);
  coarse_flux_func.imag().ProjectCoefficient(imag_coef);

  local_timer.solve_time += local_timer.Lap();

  // Loop over elements, embed the smooth flux into the coarse flux space, then
  // compute squared integral using a component-wise mass matrix.
  Vector smooth_vec, coarse_vec, sub_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;

  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // real
    smooth_flux_func.real().GetElementDofValues(e, smooth_vec);
    coarse_flux_func.real().GetElementDofValues(e, coarse_vec);

    const int ndof = coarse_vec.Size() / 3;
    sub_vec.SetSize(ndof);
    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      normalization += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    smooth_to_coarse_embed[e].AddMult(smooth_vec, coarse_vec, -1.0);
    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      estimates[e] += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    // imag
    smooth_flux_func.imag().GetElementDofValues(e, smooth_vec);
    coarse_flux_func.imag().GetElementDofValues(e, coarse_vec);

    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      normalization += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    smooth_to_coarse_embed[e].AddMult(smooth_vec, coarse_vec, -1.0);
    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      estimates[e] += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    estimates[e] = std::sqrt(estimates[e]);
  }
  local_timer.est_solve_time += local_timer.Lap();

  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);

  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  local_timer.postpro_time += local_timer.Lap();

  if constexpr (false)
  {
    Mpi::Print("CurlFluxErrorEstimation Profiling");
    local_timer.Reduce(Mpi::World());

    using namespace std::chrono;
    // clang-format off
    constexpr int w = 5;  // Total column width
    Mpi::Print("\n"
               "Estimation Time Report (ms)\n"
               "============================{}\n"
               "RHS Construction            {:{}}\n"
               "Smooth Flux Solve           {:{}}\n"
               "GridFunction construction   {:{}}\n"
               "Estimator Evaluation        {:{}}\n"
               "Normalization               {:{}}\n"
               "----------------------------{}\n"
               "Total Estimation            {:{}}\n",
               std::string(w, '='),
               duration_cast<milliseconds>(local_timer.construct_time).count(), w,
               duration_cast<milliseconds>(local_timer.solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_construction_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.postpro_time).count(), w,
               std::string(w, '-'),
               duration_cast<milliseconds>(local_timer.total_time).count(), w);
    // clang-format on
  }

  return estimates;
}

Vector CurlFluxErrorEstimator::operator()(const Vector &v) const
{
  Timer local_timer;
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient coef(field, mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &flux_fes, auto &coef)
  {
    Vector RHS(flux_fes.GetTrueVSize());

    mfem::LinearForm rhs(&flux_fes);
    rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    flux_fes.GetProlongationMatrix()->MultTranspose(rhs, RHS);

    return RHS;
  };

  const auto smooth_flux_rhs = rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), coef);
  local_timer.construct_time += local_timer.Lap();

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_flux = [](const FluxProjector &proj, const Vector &flux_coef)
  {
    // Use a copy construction to match appropriate size.
    Vector flux(flux_coef);
    proj.Mult(flux_coef, flux);
    return flux;
  };
  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);

  local_timer.solve_time += local_timer.Lap();

  // Given a complex solution represented with a PetscParVector, build a
  // ComplexGridFunction for evaluation.
  auto build_func = [](const Vector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction flux(&fes);
    flux.SetFromTrueDofs(f);
    flux.ExchangeFaceNbrData();
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

  mfem::ParGridFunction coarse_flux_func(&coarse_flux_fes);

  local_timer.construct_time += local_timer.Lap();
  coarse_flux_func.ProjectCoefficient(coef);

  local_timer.solve_time += local_timer.Lap();

  // Loop over elements, embed the smooth flux into the coarse flux space, then
  // compute squared integral using a component-wise mass matrix.
  Vector smooth_vec, coarse_vec, sub_vec, estimates(nelem);
  estimates = 0.0;
  double normalization = 0.0;

  for (int e = 0; e < fes.GetNE(); ++e)
  {
    smooth_flux_func.GetElementDofValues(e, smooth_vec);
    coarse_flux_func.GetElementDofValues(e, coarse_vec);

    const int ndof = coarse_vec.Size() / 3;
    sub_vec.SetSize(ndof);
    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      normalization += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    smooth_to_coarse_embed[e].AddMult(smooth_vec, coarse_vec, -1.0);
    for (int c = 0; c < 3; ++c)
    {
      sub_vec.MakeRef(coarse_vec, c * ndof);
      estimates[e] += scalar_mass_matrices[e].InnerProduct(sub_vec, sub_vec);
    }

    estimates[e] = std::sqrt(estimates[e]);
  }
  local_timer.est_solve_time += local_timer.Lap();

  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);

  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  local_timer.postpro_time += local_timer.Lap();

  if constexpr (false)
  {
    Mpi::Print("CurlFluxErrorEstimation Profiling");
    local_timer.Reduce(Mpi::World());

    using namespace std::chrono;
    // clang-format off
    constexpr int w = 5;  // Total column width
    Mpi::Print("\n"
               "Estimation Time Report (ms)\n"
               "============================{}\n"
               "RHS Construction            {:{}}\n"
               "Smooth Flux Solve           {:{}}\n"
               "GridFunction construction   {:{}}\n"
               "Estimator Evaluation        {:{}}\n"
               "Normalization               {:{}}\n"
               "----------------------------{}\n"
               "Total Estimation            {:{}}\n",
               std::string(w, '='),
               duration_cast<milliseconds>(local_timer.construct_time).count(), w,
               duration_cast<milliseconds>(local_timer.solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_construction_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.postpro_time).count(), w,
               std::string(w, '-'),
               duration_cast<milliseconds>(local_timer.total_time).count(), w);
    // clang-format on
  }

  return estimates;
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
      iodata.solver.order, mesh.back()->Dimension(), iodata.solver.linear.mg_max_levels,
        iodata.solver.linear.mg_coarsen_type, false)),
    smooth_flux_component_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        iodata.solver.linear.mg_max_levels, iodata.solver.linear.mg_legacy_transfer,
      iodata.solver.pa_order_threshold, mesh, smooth_flux_fecs)),
    smooth_flux_fes(mesh.back().get(), smooth_flux_fecs.back().get(), mesh.back()->Dimension()),
    smooth_projector(smooth_flux_component_fes, iodata.solver.linear.tol, 200, 0, iodata.solver.pa_order_threshold),
    coarse_flux_fec(iodata.solver.order, mesh.back()->Dimension(),
                    mfem::BasisType::GaussLobatto),
    coarse_flux_fes(mesh.back().get(), &coarse_flux_fec, mesh.back()->Dimension()),
    scalar_mass_matrices(fes.GetNE()), smooth_to_coarse_embed(fes.GetNE())
{
  mfem::MassIntegrator mass_integrator;

  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto &coarse_fe = *coarse_flux_fes.GetFE(e);
    auto &T = *fes.GetElementTransformation(e);
    mass_integrator.AssembleElementMatrix(coarse_fe, T, scalar_mass_matrices[e]);

    const auto &smooth_fe = *smooth_flux_component_fes.GetFinestFESpace().GetFE(e);
    coarse_fe.Project(smooth_fe, T, smooth_to_coarse_embed[e]);
  }
}

Vector GradFluxErrorEstimator::operator()(const Vector &v) const
{
  Timer local_timer;
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  const int nelem = smooth_flux_fes.GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (V, ϵ ∇ ϕ).
  GradFluxCoefficient coef(field, mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &fes, auto &coef)
  {
    Vector RHS(fes.GetTrueVSize());

    mfem::LinearForm rhs(&fes);
    rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
    rhs.UseFastAssembly(false);
    rhs.Assemble();
    fes.GetProlongationMatrix()->MultTranspose(rhs, RHS);

    return RHS;
  };

  auto smooth_flux_rhs = rhs_from_coef(smooth_flux_fes, coef);
  local_timer.construct_time += local_timer.Lap();

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform component wise mass matrix inversion in the appropriate space,
  // giving fᵢ = M⁻¹ f̂ᵢ.
  auto build_flux = [](const FluxProjector &proj, Vector &rhs)
  {
    // Use a copy construction to match appropriate size.
    Vector flux(rhs.Size());
    flux = 0.0;

    // Apply the flux projector component wise.
    const int ndof = flux.Size();
    const int stride = ndof / 3;
    MFEM_ASSERT(ndof % 3 == 0, "!");

    Vector flux_comp, rhs_comp;
    for (std::size_t i = 0; i < 3; ++i)
    {
      flux_comp.MakeRef(flux, i * stride, stride);
      rhs_comp.MakeRef(rhs, i * stride, stride);
      proj.Mult(rhs_comp, flux_comp);
    }

    return flux;
  };

  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);

  local_timer.solve_time += local_timer.Lap();

  // Given a solution represented with a Vector, build a GridFunction for evaluation.
  auto build_func = [](const Vector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction flux(&fes);
    flux.SetFromTrueDofs(f);
    flux.ExchangeFaceNbrData();
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes);

  mfem::ParGridFunction coarse_flux(&coarse_flux_fes);
  local_timer.construct_time += local_timer.Lap();

  coarse_flux.ProjectCoefficient(coef);

  local_timer.solve_time += local_timer.Lap();

  Vector coarse_vec, smooth_vec, coarse_sub_vec, smooth_sub_vec, estimates(nelem);
  estimates = 0.0;

  double normalization = 0.0;
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    coarse_flux.GetElementDofValues(e, coarse_vec);
    smooth_flux_func.GetElementDofValues(e, smooth_vec);

    MFEM_ASSERT(coarse_vec.Size() == smooth_vec.Size() && coarse_vec.Size() % 3 == 0,
                "H1d and L2d spaces should have the same number of elemental dofs, "
                "and should be exactly divisible by 3: "
                    << coarse_vec.Size() << " " << smooth_vec.Size());

    const int ndof = coarse_vec.Size() / 3;
    coarse_sub_vec.SetSize(ndof);
    smooth_sub_vec.SetSize(ndof);

    for (int c = 0; c < 3; ++c)
    {
      coarse_sub_vec.MakeRef(coarse_vec, c * ndof);
      smooth_sub_vec.MakeRef(smooth_vec, c * ndof);

      normalization += scalar_mass_matrices[e].InnerProduct(coarse_sub_vec, coarse_sub_vec);

      // Embed
      smooth_to_coarse_embed[e].AddMult(smooth_sub_vec, coarse_sub_vec, -1.0);

      // Integrate
      estimates[e] += scalar_mass_matrices[e].InnerProduct(coarse_sub_vec, coarse_sub_vec);
    }

    estimates[e] = std::sqrt(estimates[e]);
  }

  Mpi::GlobalSum(1, &normalization, field.ParFESpace()->GetComm());
  normalization = std::sqrt(normalization);

  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  local_timer.postpro_time += local_timer.Lap();

  if constexpr (false)
  {
    // Debugging branch generates some intermediate fields for paraview.
    // Also reports the local timing results from above.
    mfem::ParaViewDataCollection paraview("debug", fes.GetParMesh());
    paraview.RegisterVCoeffField("Flux", &coef);

    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes);
    paraview.RegisterField("SmoothFlux", &smooth_flux_func);

    mfem::L2_FECollection est_fec(0, 3);
    mfem::ParFiniteElementSpace est_fes(fes.GetParMesh(), &est_fec);
    mfem::ParGridFunction est_field(&est_fes);
    est_field.SetFromTrueDofs(estimates);

    paraview.RegisterField("ErrorIndicator", &est_field);

    paraview.Save();
  }

  if constexpr (false)
  {
    Mpi::Print("GradFluxErrorEstimation Profiling");
    local_timer.Reduce(Mpi::World());

    using namespace std::chrono;
    // clang-format off
    constexpr int w = 5;  // Total column width
    Mpi::Print("\n"
               "Estimation Time Report (ms)\n"
               "============================{}\n"
               "RHS Construction            {:{}}\n"
               "Smooth Flux Solve           {:{}}\n"
               "GridFunction construction   {:{}}\n"
               "Estimator Evaluation        {:{}}\n"
               "Normalization               {:{}}\n"
               "----------------------------{}\n"
               "Total Estimation            {:{}}\n",
               std::string(w, '='),
               duration_cast<milliseconds>(local_timer.construct_time).count(), w,
               duration_cast<milliseconds>(local_timer.solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_construction_time).count(), w,
               duration_cast<milliseconds>(local_timer.est_solve_time).count(), w,
               duration_cast<milliseconds>(local_timer.postpro_time).count(), w,
               std::string(w, '-'),
               duration_cast<milliseconds>(local_timer.total_time).count(), w);
    // clang-format on
  }

  return estimates;
}
}  // namespace palace
