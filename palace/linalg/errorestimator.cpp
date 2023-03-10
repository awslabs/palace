// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include "fem/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"
#include "utils/multigrid.hpp"

namespace palace
{

namespace
{
template <typename SmoothFluxFECollection>
struct PropertySelector;

template <>
struct PropertySelector<mfem::H1_FECollection>
{
  static constexpr auto mat = MaterialPropertyType::PERMITTIVITY_ABS;
  static constexpr auto dim = 3;
  using FluxFEC = mfem::L2_FECollection;
  using FluxCoefficient = GradFluxCoefficient;
  using Integrator = mfem::DiffusionIntegrator;
};

template <>
struct PropertySelector<mfem::ND_FECollection>
{
  static constexpr auto mat = MaterialPropertyType::INV_PERMEABILITY;
  static constexpr auto dim = 1;
  using FluxFEC = mfem::RT_FECollection;
  using FluxCoefficient = CurlFluxCoefficient;
  using Integrator = mfem::CurlCurlIntegrator;
};

}  // namespace

using namespace utils;

CurlFluxErrorEstimator::CurlFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    flux_fes(mesh.back().get(), &flux_fec),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    projector(smooth_flux_fes, iodata.solver.linear.tol, iodata.solver.linear.max_it, 0)
{
}

CurlFluxErrorEstimator::CurlFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh->Dimension()),
    flux_fes(mesh.get(), &flux_fec),
    smooth_flux_fecs(ConstructFECollections<mfem::ND_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        mesh, smooth_flux_fecs)),
    projector(smooth_flux_fes, iodata.solver.linear.tol, iodata.solver.linear.max_it, 0)
{
}

mfem::Vector CurlFluxErrorEstimator::operator()(const petsc::PetscParVector &v,
                                                bool use_mfem) const
{
  const auto cv = v.GetToVectors();
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(cv.real);
  field.imag().SetFromTrueDofs(cv.imag);

  constexpr int normp = 1;
  if (use_mfem)
  {
    const int order = fes.GetElementOrder(0);  // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> coef(mat_op);
    mfem::CurlCurlIntegrator flux_integrator(coef);  // L2ZZErrorEstimator ignores

    mfem::Vector real_error;
    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.real(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               real_error, normp));

    mfem::Vector imag_error;
    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.imag(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               imag_error, normp));

    // Compute the magnitude of the complex valued error.
    auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
    mfem::Vector estimates(real_error.Size());
    std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                   estimates.begin(), magnitude);

    return estimates;
  }
  else
  {
    // Home rolled L2ZZErrorEstimator.

    // Compute the non-smooth flux RHS., i.e. (W, μ⁻¹∇ × V)
    const auto flux = [&]()
    {
      CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);
      const auto ndof = smooth_flux_fes.GetFinestFESpace().GetTrueVSize();

      ComplexVector flux(ndof);
      {
        mfem::ParLinearForm rhs(&flux_fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(real_coef));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.real);
      }

      {
        mfem::ParLinearForm rhs(&flux_fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(imag_coef));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.imag);
      }

      return petsc::PetscParVector(v.GetComm(), flux);
    }();

    // Given the RHS vector of non-smooth flux, construct a flux projector and
    // perform mass matrix inversion in the appropriate space.
    const auto smooth_flux = [&]()
    {
      // use copy for convenience
      petsc::PetscParVector smooth_flux(flux);
      projector.Mult(flux, smooth_flux);
      return smooth_flux;
    }();

    // Given the two fluxes, create grid functions in order to allow integration
    // over each element.
    auto build_func = [this](const petsc::PetscParVector &f, auto &flux_fes)
    {
      mfem::ParComplexGridFunction flux(&flux_fes);
      const auto fi = f.GetToVectors();
      flux.real().SetFromTrueDofs(fi.real);
      flux.imag().SetFromTrueDofs(fi.imag);
      flux.real().ExchangeFaceNbrData();
      flux.imag().ExchangeFaceNbrData();
      return flux;
    };

    auto flux_func = build_func(flux, flux_fes);
    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

    // Given the flux and the smooth flux, compute the norm of the difference
    // between the two on each element.
    const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();
    std::vector<std::complex<double>> estimate_components;
    estimate_components.reserve(nelem);

    for (int i = 0; i < nelem; i++)
    {
      const auto real_error = mfem::ComputeElementLpDistance(normp, i, flux_func.real(),
                                                             smooth_flux_func.real());

      const auto imag_error = mfem::ComputeElementLpDistance(normp, i, flux_func.imag(),
                                                             smooth_flux_func.imag());

      estimate_components.emplace_back(real_error, imag_error);
    }

    // Compute the magnitude of the complex valued error.
    auto magnitude = [](const auto &x)
    { return std::sqrt(x.real() * x.real() + x.imag() * x.imag()); };
    mfem::Vector estimates(estimate_components.size());
    std::transform(estimate_components.begin(), estimate_components.end(),
                   estimates.begin(), magnitude);

    return estimates;
  }
}

GradFluxErrorEstimator::GradFluxErrorEstimator(
    const IoData &iodata, const MaterialOperator &mat_op,
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh.back()->Dimension()),
    flux_fes(mesh.back().get(), &flux_fec, mesh.back()->Dimension()),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order,
        mesh.back()->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh.back()->Dimension())),
    projector(smooth_flux_fes, iodata.solver.linear.tol, iodata.solver.linear.max_it, 0)
{
}

GradFluxErrorEstimator::GradFluxErrorEstimator(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               std::unique_ptr<mfem::ParMesh> &mesh,
                                               mfem::ParFiniteElementSpace &fes)
  : mat_op(mat_op), fes(fes), flux_fec(iodata.solver.order - 1, mesh->Dimension()),
    flux_fes(mesh.get(), &flux_fec, mesh->Dimension()),
    smooth_flux_fecs(ConstructFECollections<mfem::H1_FECollection>(
        iodata.solver.linear.mat_gmg, false, iodata.solver.order, mesh->Dimension())),
    smooth_flux_fes(utils::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        mesh, smooth_flux_fecs, mesh->Dimension())),
    projector(smooth_flux_fes, iodata.solver.linear.tol, iodata.solver.linear.max_it, 0)
{
}

mfem::Vector GradFluxErrorEstimator::operator()(const mfem::Vector &v, bool use_mfem) const
{
  mfem::ParGridFunction field(&fes);
  field.SetFromTrueDofs(v);

  constexpr int normp = 1;
  if (use_mfem)
  {
    const int order = fes.GetElementOrder(0);  // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_ABS> coef(mat_op);
    mfem::DiffusionIntegrator flux_integrator(coef);

    mfem::Vector estimates;
    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field,
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               estimates, normp));

    return estimates;
  }
  else
  {
    // Home rolled L2ZZErrorEstimator.

    // Compute the non-smooth flux RHS, i.e. (w, ϵ ∇ ϕ)
    const auto flux = [&]()
    {
      CurlFluxCoefficient coef(field, mat_op);
      const auto ndof = smooth_flux_fes.GetFinestFESpace().GetTrueVSize();

      mfem::Vector flux(ndof);

      mfem::ParLinearForm rhs(&flux_fes);
      rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(coef));
      rhs.UseFastAssembly(true);
      rhs.Assemble();
      rhs.ParallelAssemble(flux);

      return flux;
    }();

    // Given the RHS vector of non-smooth flux, construct a flux projector and
    // perform mass matrix inversion in the appropriate space.
    const auto smooth_flux = [&]()
    {
      mfem::Vector smooth_flux(flux.Size());
      projector.Mult(flux, smooth_flux);
      return smooth_flux;
    }();

    // Given the two fluxes, create grid functions in order to allow integration
    // over each element.
    auto build_func = [this](const mfem::Vector &f, auto &flux_fes)
    {
      mfem::ParGridFunction flux(&fes);
      flux.SetFromTrueDofs(f);
      flux.ExchangeFaceNbrData();
      return flux;
    };

    auto flux_func = build_func(flux, flux_fes);
    auto smooth_flux_func = build_func(smooth_flux, smooth_flux_fes.GetFinestFESpace());

    // Given the flux and the smooth flux, compute the norm of the difference
    // between the two on each element.
    const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();
    mfem::Vector estimates(nelem);

    for (int i = 0; i < nelem; i++)
    {
      estimates[i] = mfem::ComputeElementLpDistance(normp, i, flux_func, smooth_flux_func);
    }

    return estimates;
  }
}

}  // namespace palace
