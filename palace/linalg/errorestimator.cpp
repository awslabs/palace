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
    projector(smooth_flux_fes, 1e-12, 200)
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
    projector(smooth_flux_fes, 1e-12, 200)
{
}

mfem::Vector CurlFluxErrorEstimator::operator()(const petsc::PetscParVector &v,
                                                bool use_mfem) const
{
  const auto cv = v.GetToVectors();
  mfem::ParComplexGridFunction field(&fes);
  field.real().SetFromTrueDofs(cv.real);
  field.imag().SetFromTrueDofs(cv.imag);

  const int nelem = smooth_flux_fes.GetFinestFESpace().GetNE();
  mfem::Vector real_error(nelem), imag_error(nelem);

  constexpr int normp = 1;
  if (use_mfem)
  {
    const int order = fes.GetElementOrder(0);  // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> coef(mat_op);
    mfem::CurlCurlIntegrator flux_integrator(coef);  // L2ZZErrorEstimator ignores

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.real(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               real_error, normp));

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.imag(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               imag_error, normp));

  }
  else
  {
    // Home rolled L2ZZErrorEstimator.

    // Compute the non-smooth flux RHS., i.e. (W, μ⁻¹∇ × V)
    CurlFluxCoefficient real_coef(field.real(), mat_op), imag_coef(field.imag(), mat_op);

    // This lambda computes the flux as achieved within mfem. Use this to bench against.
    auto mfem_flux = [&]()
    {
      mfem::ParComplexGridFunction flux_func(&flux_fes);
      flux_func.real() = 0.0;
      flux_func.imag() = 0.0;

      mfem::Array<int> xdofs, fdofs;
      mfem::Vector el_x, el_f, alt_el_f;
      mfem::CurlCurlIntegrator curl;

      for (int i = 0; i < fes.GetNE(); i++)
      {
        auto *T = fes.GetElementTransformation(i);

        fes.GetElementVDofs(i, xdofs);
        flux_fes.GetElementVDofs(i, fdofs);

        field.real().GetSubVector(xdofs, el_x);
        curl.ComputeElementFlux(*field.real().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);

        flux_func.real().AddElementVector(fdofs, el_f);

        field.imag().GetSubVector(xdofs, el_x);
        curl.ComputeElementFlux(*field.imag().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);
        flux_func.imag().AddElementVector(fdofs, el_f);
      }

      const auto ndof = smooth_flux_fes.GetFinestFESpace().GetTrueVSize();

      ComplexVector flux(ndof);
      {
        mfem::VectorGridFunctionCoefficient f_real(&flux_func.real());
        mfem::ParLinearForm rhs(&smooth_flux_fes.GetFinestFESpace());
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_real));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.real);
      }

      {
        mfem::VectorGridFunctionCoefficient f_imag(&flux_func.imag());
        mfem::ParLinearForm rhs(&smooth_flux_fes.GetFinestFESpace());
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_imag));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.imag);
      }

      return petsc::PetscParVector(v.GetComm(), flux);
    };

    auto palace_flux = [&]()
    {
      // TODO: Make a drop in here that uses muinv properly, but also recovers
      // the mfem_flux for muinv of identity.
    };


    const auto flux = mfem_flux();


    // Given the RHS vector of non-smooth flux, construct a flux projector and
    // perform mass matrix inversion in the appropriate space.
    // f = M⁻¹ f̂
    auto smooth_flux = [&]()
    {
      // use copy for convenience
      petsc::PetscParVector smooth_flux(flux);
      projector.Mult(flux, smooth_flux);
      return smooth_flux;
    }();

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

    auto flux_func = build_func(flux, smooth_flux_fes.GetFinestFESpace());

    for (int i = 0; i < smooth_flux_fes.GetFinestFESpace().GetNE(); i++)
    {
      real_error(i) = ComputeElementLpDistance(normp, i, smooth_flux_func.real(), flux_func.real());
      imag_error(i) = ComputeElementLpDistance(normp, i, smooth_flux_func.imag(), flux_func.imag());
    }
  }

  // Compute the magnitude of the complex valued error.
  auto magnitude = [](const auto &r, const auto &i) { return std::sqrt(r * r + i * i); };
  mfem::Vector estimates(real_error.Size());
  std::transform(real_error.begin(), real_error.end(), imag_error.begin(),
                  estimates.begin(), magnitude);

  return estimates;
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
    projector(smooth_flux_fes, 1e-12, 200)
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
    projector(smooth_flux_fes, 1e-12, 200)
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
