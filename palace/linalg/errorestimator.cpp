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

namespace palace
{

namespace
{

// Given a grid function defining a vector solution, compute the error relative
// to a vector coefficient. The default quadrature rule exactly integrates 2p +
// q polynomials, but the quadrature order can be increased or decreased via
// quad_order_increment.
mfem::Vector ComputeElementLpErrors(const mfem::ParGridFunction &sol, double p,
                                    mfem::VectorCoefficient &exsol,
                                    int quad_order_increment = 0);

// Given a grid function defining a vector/scalar solution, compute the Lp norm of the
// solution. The default quadrature rule exactly integrates 2p +
// q polynomials, but the quadrature order can be increased or decreased via
// quad_order_increment for non L2 norms.
double ComputeVectorLpNorm(const mfem::ParGridFunction &sol, double p = 2,
                           int quad_order_increment = 0);
double ComputeScalarLpNorm(const mfem::ParGridFunction &sol, double p = 2,
                           int quad_order_increment = 0);

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
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
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
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
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

  const auto flux =
      ComplexVector(rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
                    rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), imag_coef));

  const auto pflux = petsc::PetscParVector(v.GetComm(), flux);

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

  real_error = ComputeElementLpErrors(smooth_flux_func.real(), normp, real_coef);
  imag_error = ComputeElementLpErrors(smooth_flux_func.imag(), normp, imag_coef);

  // Compute the magnitude of the complex valued error.
  auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
  mfem::Vector estimates(real_error.Size());
  std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                 estimates.begin(), magnitude);

  // Normalize the error by the solution L2 norm.
  const auto normalization = std::sqrt(std::pow(ComputeVectorLpNorm(field.real(), 2), 2.0) +
                                       std::pow(ComputeVectorLpNorm(field.imag(), 2), 2.0));

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
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
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
    projector(smooth_flux_fes, iodata.solver.linear.tol / 10, 200, iodata.problem.verbose)
{
}

mfem::Vector GradFluxErrorEstimator::operator()(const mfem::Vector &v) const
{
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

  const auto flux = rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), coef);

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

  // Given a solution represented with a Vector, build a GridFunction for evaluation.
  auto build_func = [](const mfem::Vector &f, mfem::ParFiniteElementSpace &fes)
  {
    mfem::ParGridFunction flux(&fes);
    flux.SetFromTrueDofs(f);
    return flux;
  };

  auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

  auto estimates = ComputeElementLpErrors(smooth_flux_func, normp, coef);

  // Normalize the error by the solution L2 norm, ensures reductions are well scaled.
  const auto normalization = ComputeScalarLpNorm(field, 2);
  std::for_each(estimates.begin(), estimates.end(),
                [&normalization](auto &x) { x /= normalization; });

  return estimates;
}

namespace
{
mfem::Vector ComputeElementLpErrors(const mfem::ParGridFunction &sol, double p,
                                    mfem::VectorCoefficient &exsol,
                                    int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::Vector error(fes.GetNE());
  error = 0.0;

  mfem::DenseMatrix vals, exact_vals;
  mfem::Vector loc_errs;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);

    sol.GetVectorValues(T, ir, vals);
    exsol.Eval(exact_vals, T, ir);

    vals -= exact_vals;

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      // Combine components
      double component = 0.0;
      for (int c = 0; c < vals.Height(); ++c)
      {
        component += std::pow(std::abs(vals(c, j)), p);
      }

      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      error[i] += ip.weight * T.Weight() * component;
    }
    error[i] = std::pow(std::abs(error[i]), 1 / p);
  }
  return error;
}

double ComputeVectorLpNorm(const mfem::ParGridFunction &sol, double p,
                           int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::DenseMatrix vals;
  mfem::Vector loc, elem_norm(nelem);

  elem_norm = 0.0;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
    sol.GetVectorValues(T, ir, vals);

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      // Combine components
      double component = 0.0;
      for (int c = 0; c < vals.Height(); ++c)
      {
        component += std::pow(std::abs(vals(c, j)), p);
      }

      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);

      elem_norm[i] += ip.weight * T.Weight() * component;
    }

    // Negative quadrature weights might case elemental norms to have been
    // negative, correct this before accumulation.
    elem_norm[i] = std::abs(elem_norm[i]);
  }

  // (∑ₑ ∫ₑ|u|^p)^(1/p)
  auto norm = std::accumulate(elem_norm.begin(), elem_norm.end(), 0.0);

  Mpi::GlobalSum(1, &norm, Mpi::World());

  norm = std::pow(norm, 1.0 / p);

  return norm;
}

double ComputeScalarLpNorm(const mfem::ParGridFunction &sol, double p,
                           int quad_order_increment)
{
  auto &fes = *sol.ParFESpace();
  const int nelem = fes.GetNE();

  MFEM_VERIFY(p < mfem::infinity(), "Do not use this routine to compute Linf norm\n");

  mfem::Vector vals, elem_norm(nelem);

  elem_norm = 0.0;

  for (int i = 0; i < fes.GetNE(); ++i)
  {
    const auto &fe = *fes.GetFE(i);
    auto &T = *fes.GetElementTransformation(i);
    const auto &ir = *utils::GetDefaultRule(fe, T, quad_order_increment);
    sol.GetValues(T, ir, vals);

    for (int j = 0; j < ir.GetNPoints(); ++j)
    {
      const auto &ip = ir.IntPoint(j);
      T.SetIntPoint(&ip);
      elem_norm[i] += ip.weight * T.Weight() * std::pow(std::abs(vals(j)), p);
    }

    // Negative quadrature weights might cause elemental norms to have been
    // negative.
    elem_norm[i] = std::abs(elem_norm[i]);
  }

  // (∑ₑ ∫ₑ|u|^p)^(1/p)
  auto norm = std::accumulate(elem_norm.begin(), elem_norm.end(), 0.0);

  Mpi::GlobalSum(1, &norm, Mpi::World());

  norm = std::pow(norm, 1.0 / p);

  return norm;
}
}  // namespace

}  // namespace palace
