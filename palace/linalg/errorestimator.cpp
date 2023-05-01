// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include "fem/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/errorindicators.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"
#include "utils/mfemintegrators.hpp"
#include "utils/multigrid.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Given a grid function defining a vector solution, compute the error relative
// to a vector coefficient.
mfem::Vector ComputeElementL2Errors(const mfem::ParGridFunction &sol,
                                    mfem::VectorCoefficient &exsol, int quad_order_increment = 0);

// Given a grid function defining a vector/scalar solution, compute the L2 norm of the
// solution.
double ComputeVectorL2Norm(const mfem::ParGridFunction &sol);
double ComputeScalarL2Norm(const mfem::ParGridFunction &sol);

}  // namespace

using namespace utils;

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    smooth_projector(smooth_flux_fes, 1e-10, 200, iodata.problem.verbose),
    mass(&fes),
    coarse_flux_fecs(ConstructFECollections<mfem::L2_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    coarse_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(
        mesh, coarse_flux_fecs, mesh.back()->Dimension())),
    coarse_projector(coarse_flux_fes, 1e-10, 200, iodata.problem.verbose),
    smooth_to_coarse(&smooth_flux_fes.GetFinestFESpace(), &coarse_flux_fes.GetFinestFESpace())
{
  mfem::MassIntegrator mass_integrator;

  scalar_mass_matrices.resize(fes.GetNE());
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto * const fe = coarse_flux_fes.GetFinestFESpace().GetFE(e);
    auto * const T = coarse_flux_fes.GetFinestFESpace().GetElementTransformation(e);
    mass_integrator.AssembleElementMatrix(*fe, *T, scalar_mass_matrices[e]);
  }

  mass.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
  mass.Assemble();
  mass.Finalize();

  smooth_to_coarse.AddDomainIntegrator(new mfem::IdentityInterpolator);
  smooth_to_coarse.Assemble();
  smooth_to_coarse.Finalize();
}

CurlFluxErrorEstimator::CurlFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    smooth_projector(smooth_flux_fes, 1e-10, 200, iodata.problem.verbose),
    mass(&fes),
    coarse_flux_fecs(ConstructFECollections<mfem::L2_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh->Dimension())),
    coarse_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(
        mesh, coarse_flux_fecs, mesh->Dimension())),
    coarse_projector(coarse_flux_fes, 1e-10, 200, iodata.problem.verbose),
    smooth_to_coarse(&smooth_flux_fes.GetFinestFESpace(), &coarse_flux_fes.GetFinestFESpace())
{
  mfem::MassIntegrator mass_integrator;

  scalar_mass_matrices.resize(fes.GetNE());
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto * const fe = coarse_flux_fes.GetFinestFESpace().GetFE(e);
    auto * const T = coarse_flux_fes.GetFinestFESpace().GetElementTransformation(e);
    mass_integrator.AssembleElementMatrix(*fe, *T, scalar_mass_matrices[e]);
  }

  mass.AddDomainIntegrator(new mfem::VectorFEMassIntegrator);
  mass.Assemble();
  mass.Finalize();

  smooth_to_coarse.AddDomainIntegrator(new mfem::IdentityInterpolator);
  smooth_to_coarse.Assemble();
  smooth_to_coarse.Finalize();
}

mfem::Vector CurlFluxErrorEstimator::operator()(const petsc::PetscParVector &v) const
{
  const auto cv = v.GetToVectors();
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(cv.real);
  field.imag().SetFromTrueDofs(cv.imag);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &flux_fes, auto &coef)
  {
    mfem::Vector RHS(flux_fes.GetTrueVSize());

    mfem::ParLinearForm rhs(&flux_fes);
    if (flux_fes.FEColl()->GetContType() == mfem::FiniteElementCollection::TANGENTIAL)
    {
      rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    }
    else
    {
      rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
    }
    rhs.UseFastAssembly(true);
    rhs.Assemble();
    rhs.ParallelAssemble(RHS);

    return RHS;
  };

  const auto smooth_flux_rhs = petsc::PetscParVector(
      v.GetComm(), rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
      rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), imag_coef));

  const auto coarse_flux_rhs = petsc::PetscParVector(
      v.GetComm(), rhs_from_coef(coarse_flux_fes.GetFinestFESpace(), real_coef),
      rhs_from_coef(coarse_flux_fes.GetFinestFESpace(), imag_coef));

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_flux = [this](const FluxProjector &proj, const petsc::PetscParVector &flux_coef)
  {
    // Use a copy construction to match appropriate size.
    petsc::PetscParVector flux(flux_coef);
    proj.Mult(flux_coef, flux);
    return flux;
  };
  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);
  auto coarse_flux = build_flux(coarse_projector, coarse_flux_rhs);

  auto smooth_flux_components = smooth_flux.GetToVectors();
  auto coarse_flux_components = coarse_flux.GetToVectors();

  mfem::Array<const mfem::Vector*> smooth_comps(2);
  mfem::Array<mfem::Vector*> coarse_comps(2);
  smooth_comps[0] = &smooth_flux_components.real;
  smooth_comps[1] = &smooth_flux_components.imag;
  coarse_comps[0] = &coarse_flux_components.real;
  coarse_comps[1] = &coarse_flux_components.imag;

  std::cout << "coarse_flux_components " << coarse_flux_components.real.Norml2() << " " << coarse_flux_components.imag.Norml2() << '\n';
  smooth_to_coarse.ArrayAddMult(smooth_comps, coarse_comps, -1.0);
  std::cout << "coarse_flux_components " << coarse_flux_components.real.Norml2() << " " << coarse_flux_components.imag.Norml2() << '\n';

  // Compute using the more optimized approach.

  mfem::Array<int> dofs;
  mfem::Vector sub_vec;
  mfem::Vector real_error2(nelem), imag_error2(nelem), estimates2(nelem); // TODO: combine these
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    coarse_flux_fes.GetFinestFESpace().GetElementVDofs(e, dofs);

    coarse_flux_components.real.GetSubVector(dofs, sub_vec);
    const int ndof = sub_vec.Size() / 3;
    real_error2[e] = scalar_mass_matrices[e].InnerProduct(sub_vec.GetData(), sub_vec.GetData());
    for (int c = 1; c < 3; ++c)
    {
      double *slice = sub_vec.GetData() + c * ndof;
      real_error2[e] += scalar_mass_matrices[e].InnerProduct(slice, slice);
    }

    coarse_flux_components.imag.GetSubVector(dofs, sub_vec);
    imag_error2[e] = scalar_mass_matrices[e].InnerProduct(sub_vec.GetData(), sub_vec.GetData());
    for (int c = 1; c < 3; ++c)
    {
      double *slice = sub_vec.GetData() + c * ndof;
      imag_error2[e] += scalar_mass_matrices[e].InnerProduct(slice, slice);
    }

    estimates2[e] = std::sqrt(real_error2[e] + imag_error2[e]);
    real_error2[e] = std::sqrt(real_error2[e]);
    imag_error2[e] = std::sqrt(imag_error2[e]);
  }

  // Compute using the current expensive process.

  // Given a complex solution represented with a PetscParVector, build a
  // ComplexGridFunction for evaluation.
  auto build_func = [](const petsc::PetscParVector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParComplexGridFunction flux(&fes);
    const auto fi = f.GetToVectors();
    flux.real().SetFromTrueDofs(fi.real);
    flux.imag().SetFromTrueDofs(fi.imag);
    flux.real().ExchangeFaceNbrData();
    flux.imag().ExchangeFaceNbrData();
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());
  auto coarse_flux_func = build_func(coarse_flux, coarse_flux_fes.GetFinestFESpace());

  auto real_error = ComputeElementL2Errors(smooth_flux_func.real(), real_coef);
  auto imag_error = ComputeElementL2Errors(smooth_flux_func.imag(), imag_coef);

  // Compute the magnitude of the complex valued error.
  auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
  mfem::Vector estimates(real_error.Size());
  std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                 estimates.begin(), magnitude);


  std::cout << "Delta sum real " << real_error2.Sum() - real_error.Sum() << '\n';
  std::cout << "Delta sum imag " << imag_error2.Sum() - imag_error.Sum() << '\n';
  std::cout << "Delta sum estimates " << estimates2.Sum() - estimates.Sum() << '\n';

  std::cout << "Sum real2 " << real_error2.Sum() << '\n';
  std::cout << "Sum imag2 " << imag_error2.Sum() << '\n';
  std::cout << "Sum estimates2 " << estimates2.Sum() << '\n';

  std::cout << "Sum real " << real_error.Sum() << '\n';
  std::cout << "Sum imag " << imag_error.Sum() << '\n';
  std::cout << "Sum estimates " << estimates.Sum() << '\n';

  std::cout << "Delta real " << (real_error2 -= real_error).Normlinf() << '\n';
  std::cout << "Delta imag " << (imag_error2 -= imag_error).Normlinf() << '\n';

  std::cout << "Delta estimates " << (estimates2 -= estimates).Normlinf() << '\n';

  // Normalize the error by the solution L2 norm.
  const auto normalization = std::sqrt(std::pow(ComputeVectorL2Norm(field.real()), 2.0) +
                                       std::pow(ComputeVectorL2Norm(field.imag()), 2.0));

  double normalization2 = mass.InnerProduct(cv.real, cv.real) + mass.InnerProduct(cv.imag, cv.imag);
  Mpi::GlobalSum(1, &normalization2, field.ParFESpace()->GetComm());
  normalization2 = std::sqrt(normalization2);

  std::cout << "normalization " << normalization << " normalization2 " << normalization2 << '\n';

  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  return estimates;
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh.back()->Dimension())),
    smooth_projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose),
    mass(&fes),
    coarse_flux_fecs(ConstructFECollections<mfem::L2_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    coarse_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(
        mesh, coarse_flux_fecs, mesh.back()->Dimension())),
    coarse_projector(coarse_flux_fes, iodata.solver.linear.tol, 200, iodata.problem.verbose),
    smooth_to_coarse(smooth_flux_fes.GetFinestFESpace(), coarse_flux_fes.GetFinestFESpace())
{
  auto mass_integrator = new mfem::MassIntegrator;

  scalar_mass_matrices.resize(fes.GetNE());
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto * const fe = fes.GetFE(e);
    auto * const T = fes.GetElementTransformation(e);
    mass_integrator->AssembleElementMatrix(*fe, *T, scalar_mass_matrices[e]);
  }

  mass.AddDomainIntegrator(mass_integrator);
  mass.Assemble();
  mass.Finalize();
}

GradFluxErrorEstimator::GradFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh->Dimension())),
    smooth_projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose),
    mass(&fes),
    coarse_flux_fecs(ConstructFECollections<mfem::L2_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh->Dimension())),
    coarse_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(
        mesh, coarse_flux_fecs, mesh->Dimension())),
    coarse_projector(coarse_flux_fes, iodata.solver.linear.tol, 200, iodata.problem.verbose),
    smooth_to_coarse(smooth_flux_fes.GetFinestFESpace(), coarse_flux_fes.GetFinestFESpace())
{
  auto mass_integrator = new mfem::MassIntegrator;

  scalar_mass_matrices.resize(fes.GetNE());
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    // Loop over each element, and save an elemental mass matrix.
    // Will exploit the fact that vector L2 mass matrix components are independent.
    const auto * const fe = fes.GetFE(e);
    auto * const T = fes.GetElementTransformation(e);
    mass_integrator->AssembleElementMatrix(*fe, *T, scalar_mass_matrices[e]);
  }

  mass.AddDomainIntegrator(mass_integrator);
  mass.Assemble();
  mass.Finalize();
}

mfem::Vector GradFluxErrorEstimator::operator()(const mfem::Vector &v) const
{
  Timer local_timer;
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  // Coefficients for computing the discontinuous flux., i.e. (V, ϵ ∇ ϕ).
  GradFluxCoefficient coef(field, mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &fes, auto &coef)
  {
    mfem::Vector RHS(fes.GetTrueVSize());

    mfem::ParLinearForm rhs(&fes);
    rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
    rhs.UseFastAssembly(true);
    rhs.Assemble();
    rhs.ParallelAssemble(RHS);

    return RHS;
  };

  const auto smooth_flux_rhs = rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), coef);
  local_timer.construct_time += local_timer.Lap();

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_flux = [this](const FluxProjector &proj, const mfem::Vector &flux_coef)
  {
    // Use a copy construction to match appropriate size.
    mfem::Vector flux(flux_coef);
    proj.Mult(flux_coef, flux);
    return flux;
  };
  auto smooth_flux = build_flux(smooth_projector, smooth_flux_rhs);

  local_timer.solve_time += local_timer.Lap();

  // Given a solution represented with a Vector, build a GridFunction for evaluation.
  auto build_func = [](const mfem::Vector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction flux(&fes);
    flux.SetFromTrueDofs(f);
    flux.ExchangeFaceNbrData();
    return flux;
  };

  const auto coarse_flux_rhs = rhs_from_coef(coarse_flux_fes.GetFinestFESpace(), coef);
  local_timer.construct_time += local_timer.Lap();
  auto coarse_flux = build_flux(coarse_projector, coarse_flux_rhs);
  local_timer.solve_time += local_timer.Lap();

  smooth_to_coarse.AddMult(smooth_flux, coarse_flux, -1.0);

  mfem::Array<int> dofs;
  mfem::Vector coarse_sub_vec;
  mfem::Vector estimates(nelem);
  for (int e = 0; e < fes.GetNE(); ++e)
  {
    coarse_flux_fes.GetFinestFESpace().GetElementVDofs(e, dofs);
    coarse_flux.GetSubVector(dofs, coarse_sub_vec);

    const int ndof = coarse_sub_vec.Size() / 3;
    estimates[e] = scalar_mass_matrices[e].InnerProduct(coarse_sub_vec.GetData(), coarse_sub_vec.GetData());
    for (int c = 1; c < 3; ++c)
    {
      double *slice = coarse_sub_vec.GetData() + c * ndof;
      estimates[e] += scalar_mass_matrices[e].InnerProduct(slice, slice);
    }

    estimates[e] = std::sqrt(estimates[e]);
  }

  double normalization = mass.InnerProduct(field, field);
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

    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());
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

namespace
{
// mfem::Vector ComputeElementL2Errors(const mfem::ParGridFunction &sol,
//                                     mfem::VectorCoefficient &exsol, int quad_order_increment)
// {
//   using namespace std::chrono;
//   std::vector<double> timers(7, 0.0);
//   auto &fes = *sol.ParFESpace();

//   mfem::Vector error(fes.GetNE());
//   error = 0.0;

//   mfem::DenseMatrix vals, exact_vals;

//   auto clock = steady_clock::now();
//   auto lap = steady_clock::now();
//   for (int i = 0; i < fes.GetNE(); ++i)
//   {
//     clock = steady_clock::now();
//     const auto &fe = *fes.GetFE(i);
//     auto &T = *fes.GetElementTransformation(i);
//     // Can drop one order as the coefficient is the gradient of a P-order field.
//     const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
//     lap = steady_clock::now();
//     timers[0] += duration_cast<microseconds>(lap - clock).count();
//     clock = steady_clock::now();

//     sol.GetVectorValues(T, ir, vals);

//     lap = steady_clock::now();
//     timers[1] += duration_cast<microseconds>(lap - clock).count();
//     clock = steady_clock::now();
//     exsol.Eval(exact_vals, T, ir);

//     lap = steady_clock::now();
//     timers[2] += duration_cast<microseconds>(lap - clock).count();
//     clock = steady_clock::now();

//     vals -= exact_vals;

//     lap = steady_clock::now();
//     timers[3] += duration_cast<microseconds>(lap - clock).count();
//     clock = steady_clock::now();

//     for (int j = 0; j < ir.GetNPoints(); ++j)
//     {
//       const auto &ip = ir.IntPoint(j);
//       T.SetIntPoint(&ip);

//       lap = steady_clock::now();
//       timers[4] += duration_cast<microseconds>(lap - clock).count();
//       clock = steady_clock::now();

//       error[i] += ip.weight * T.Weight() *
//         (vals(0, j)*vals(0, j) + vals(1, j)*vals(1, j) + vals(2, j)*vals(2, j));

//       lap = steady_clock::now();
//       timers[5] += duration_cast<microseconds>(lap - clock).count();
//       clock = steady_clock::now();
//     }
//     error[i] = std::sqrt(std::abs(error[i]));

//     lap = steady_clock::now();
//     timers[6] += duration_cast<microseconds>(lap - clock).count();
//     clock = steady_clock::now();
//   }

//   if constexpr (true)
//   {
//     Mpi::GlobalSum(7, timers.data(), Mpi::World());
//     if (Mpi::Root(Mpi::World()))
//     {
//       std::cout << "Timers within ComputeElementL2Errors\n";
//       std::cout << "GetIRandTransforms  " << timers[0] << " μs\n";
//       std::cout << "GetVectorValues     " << timers[1] << " μs\n";
//       std::cout << "EvaluateCoefficient " << timers[2] << " μs\n";
//       std::cout << "TakeDifference      " << timers[3] << " μs\n";
//       std::cout << "SetIntPoint         " << timers[4] << " μs\n";
//       std::cout << "WeightedSummation   " << timers[5] << " μs\n";
//       std::cout << "Sqrt                " << timers[6] << " μs\n";
//     }
//   }

//   return error;
// }


mfem::Vector ComputeElementL2Errors(const mfem::ParGridFunction &sol,
                                    mfem::VectorCoefficient &exsol, int quad_order_increment)
{
  using namespace std::chrono;
  std::vector<double> timers(7, 0.0);
  auto &fes = *sol.ParFESpace();

  mfem::Vector error(fes.GetNE());
  error = 0.0;

  mfem::Vector vals, exact_vals;

  auto clock = steady_clock::now();
  auto lap = steady_clock::now();
  for (int i = 0; i < fes.GetNE(); ++i)
  {
    clock = steady_clock::now();
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    // Can drop one order as the coefficient is the gradient of a P-order field.
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
    lap = steady_clock::now();
    timers[0] += duration_cast<microseconds>(lap - clock).count();
    clock = steady_clock::now();

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      lap = steady_clock::now();
      timers[1] += duration_cast<microseconds>(lap - clock).count();
      clock = steady_clock::now();

      sol.GetVectorValue(i, ip, vals);

      lap = steady_clock::now();
      timers[2] += duration_cast<microseconds>(lap - clock).count();
      clock = steady_clock::now();

      exsol.Eval(exact_vals, T, ip);

      lap = steady_clock::now();
      timers[3] += duration_cast<microseconds>(lap - clock).count();
      clock = steady_clock::now();

      vals -= exact_vals;

      lap = steady_clock::now();
      timers[4] += duration_cast<microseconds>(lap - clock).count();
      clock = steady_clock::now();

      error[i] += ip.weight * T.Weight() * (vals * vals);

      lap = steady_clock::now();
      timers[5] += duration_cast<microseconds>(lap - clock).count();
      clock = steady_clock::now();
    }
    error[i] = std::sqrt(std::abs(error[i]));

    lap = steady_clock::now();
    timers[6] += duration_cast<microseconds>(lap - clock).count();
    clock = steady_clock::now();
  }

  if constexpr (true)
  {
    Mpi::GlobalSum(7, timers.data(), Mpi::World());
    if (Mpi::Root(Mpi::World()))
    {
      std::cout << "Timers within ComputeElementL2Errors\n";
      std::cout << "GetIRandTransforms  " << timers[0] << " μs\n";
      std::cout << "SetIntPoint         " << timers[1] << " μs\n";
      std::cout << "GetVectorValues     " << timers[2] << " μs\n";
      std::cout << "EvaluateCoefficient " << timers[3] << " μs\n";
      std::cout << "TakeDifference      " << timers[4] << " μs\n";
      std::cout << "WeightedSummation   " << timers[5] << " μs\n";
      std::cout << "Sqrt                " << timers[6] << " μs\n";
    }
  }

  return error;
}

double ComputeVectorL2Norm(const mfem::ParGridFunction &sol)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  mfem::DenseMatrix vals;

  double norm2 = 0.0;
  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T);
    sol.GetVectorValues(T, ir, vals);

    MFEM_ASSERT(vals.Height() == 3, "!");

    double elem_norm2 = 0;
    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      // Combine components
      double component = 0.0;
      for (int c = 0; c < 3; ++c)
      {
        component += std::pow(std::abs(vals(c, j)), 2.0);
      }

      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      elem_norm2 += ip.weight * T.Weight() *
        (vals(0,j)*vals(0,j) + vals(1,j)*vals(1,j) + vals(2,j)*vals(2,j));
    }

    // Negative quadrature weights might case elemental norms to have been
    // negative, correct this before accumulation.
    norm2 += std::abs(elem_norm2);
  }

  Mpi::GlobalSum(1, &norm2, Mpi::World());

  return std::sqrt(norm2);
}

double ComputeScalarL2Norm(const mfem::ParGridFunction &sol)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  mfem::Vector vals;

  double norm2 = 0.0;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T);
    sol.GetValues(T, ir, vals);

    double elem_norm2 = 0.0;
    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);
      elem_norm2 += ip.weight * T.Weight() * vals(j) * vals(j);
    }

    // Negative quadrature weights might cause elemental norms to have been
    // negative.
    norm2 += std::abs(elem_norm2);
  }

  Mpi::GlobalSum(1, &norm2, Mpi::World());

  return std::sqrt(norm2);
}
}  // namespace

}  // namespace palace
