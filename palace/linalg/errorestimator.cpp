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
                                    mfem::VectorCoefficient &exsol);

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
    projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose)
{
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
    projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose)
{
}

mfem::Vector CurlFluxErrorEstimator::operator()(const petsc::PetscParVector &v) const
{
  const auto cv = v.GetToVectors();
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(cv.real);
  field.imag().SetFromTrueDofs(cv.imag);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();
  mfem::Vector real_error(nelem), imag_error(nelem);

  constexpr int normp = 2;  // 2 norm ensures no under integration.

  // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
  CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &smooth_flux_fes, auto &coef)
  {
    mfem::Vector RHS(smooth_flux_fes.GetTrueVSize());

    mfem::ParLinearForm rhs(&smooth_flux_fes);
    rhs.AddDomainIntegrator(new VectorFEDomainLFIntegrator(coef));
    rhs.UseFastAssembly(true);
    rhs.Assemble();
    rhs.ParallelAssemble(RHS);

    return RHS;
  };

  const auto pflux = petsc::PetscParVector(
      v.GetComm(), rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
      rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), imag_coef));

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_smooth_flux = [this](const petsc::PetscParVector &flux)
  {
    // Use a copy construction to match appropriate size.
    petsc::PetscParVector smooth_flux(flux);
    projector.Mult(flux, smooth_flux);
    return smooth_flux;
  };
  auto smooth_flux = build_smooth_flux(pflux);

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

  real_error = ComputeElementL2Errors(smooth_flux_func.real(), real_coef);
  imag_error = ComputeElementL2Errors(smooth_flux_func.imag(), imag_coef);

  // Compute the magnitude of the complex valued error.
  auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
  mfem::Vector estimates(real_error.Size());
  std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                 estimates.begin(), magnitude);

  // Normalize the error by the solution L2 norm.
  const auto normalization = std::sqrt(std::pow(ComputeVectorL2Norm(field.real()), 2.0) +
                                       std::pow(ComputeVectorL2Norm(field.imag()), 2.0));

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
    projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose)
{
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
    projector(smooth_flux_fes, iodata.solver.linear.tol * 10, 200, iodata.problem.verbose)
{
}

mfem::Vector GradFluxErrorEstimator::operator()(const mfem::Vector &v) const
{
  // Timer local_timer;
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();

  constexpr int normp = 2;  // 2 norm ensures no under integration.

  // Coefficients for computing the discontinuous flux., i.e. (V, ϵ ∇ ϕ).
  GradFluxCoefficient coef(field, mat_op);
  auto rhs_from_coef = [](mfem::ParFiniteElementSpace &smooth_flux_fes, auto &coef)
  {
    mfem::Vector RHS(smooth_flux_fes.GetTrueVSize());

    mfem::ParLinearForm rhs(&smooth_flux_fes);
    rhs.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(coef));
    rhs.UseFastAssembly(true);
    rhs.Assemble();
    rhs.ParallelAssemble(RHS);

    return RHS;
  };

  // mfem::ParaViewDataCollection paraview("debug", fes.GetParMesh());
  // paraview.SetLengthScale(1e-2);
  // paraview.RegisterVCoeffField("Flux", &coef);

  const auto flux = rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), coef);
  // local_timer.construct_time += local_timer.Lap();

  // Given the RHS vector of non-smooth flux, construct a flux projector and
  // perform mass matrix inversion in the appropriate space, giving f = M⁻¹ f̂.
  auto build_smooth_flux = [this](const mfem::Vector &flux)
  {
    // Use a copy construction to match appropriate size.
    mfem::Vector smooth_flux(flux);
    projector.Mult(flux, smooth_flux);
    return smooth_flux;
  };
  auto smooth_flux = build_smooth_flux(flux);

  // local_timer.solve_time += local_timer.Lap();

  // Given a solution represented with a Vector, build a GridFunction for evaluation.
  auto build_func = [](const mfem::Vector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction flux(&fes);
    flux.SetFromTrueDofs(f);
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

  // paraview.RegisterField("SmoothFlux", &smooth_flux_func);

  // local_timer.est_construction_time += local_timer.Lap();

  auto estimates = ComputeElementL2Errors(smooth_flux_func, coef);

  // mfem::L2_FECollection est_fec(0, 3);
  // mfem::ParFiniteElementSpace est_fes(fes.GetParMesh(), &est_fec);
  // mfem::ParGridFunction est_field(&est_fes);
  // est_field.SetFromTrueDofs(estimates);

  // paraview.RegisterField("ErrorIndicator", &est_field);

  // local_timer.est_solve_time += local_timer.Lap();

  // Normalize the error by the solution L2 norm, ensures reductions are well scaled.
  const auto normalization = ComputeScalarL2Norm(field);
  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  auto max_loc = std::max_element(estimates.begin(), estimates.end());
  // local_timer.postpro_time += local_timer.Lap();

  // mfem::ParGridFunction norm_est_field(&est_fes);
  // norm_est_field.SetFromTrueDofs(estimates);
  // paraview.RegisterField("NormalizedErrorIndicator", &norm_est_field);

  // paraview.Save();

  // Mpi::Print("GradFluxErrorEstimationProfiling");
  // local_timer.Reduce(Mpi::World());

  // local_timer.Print(Mpi::World());

  return estimates;
}

namespace
{
mfem::Vector ComputeElementL2Errors(const mfem::ParGridFunction &sol,
                                    mfem::VectorCoefficient &exsol)
{
  auto &fes = *sol.ParFESpace();

  mfem::Vector error(fes.GetNE());
  error = 0.0;

  mfem::DenseMatrix vals, exact_vals;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T);

    sol.GetVectorValues(T, ir, vals);
    exsol.Eval(exact_vals, T, ir);

    vals -= exact_vals;

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      error[i] += ip.weight * T.Weight() *
        (vals(0, j)*vals(0, j) + vals(1, j)*vals(1, j) + vals(2, j)*vals(2, j));
    }
    error[i] = std::sqrt(std::abs(error[i]));
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
        (vals(0, j)*vals(0, j) + vals(1, j)*vals(1, j) + vals(2, j)*vals(2, j));
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
