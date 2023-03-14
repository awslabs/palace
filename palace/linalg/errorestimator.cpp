// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"
#include "fem/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"
#include "utils/mfemintegrators.hpp"
#include "utils/multigrid.hpp"

namespace palace
{

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

  constexpr int normp = 2; // 2 norm ensures no under integration.
  if (use_mfem)
  {
    // This is used for comparison purposes only.
    // TODO: Delete this once the alternative direct construction is ready.
    // NB: There are bugs within L2ZZErrorEstimator presently, that ND simplices
    // with p > 1 will not work, and also a flux fec of RT elements will also
    // give incorrect answers.
    const int order = fes.GetElementOrder(0); // Assumes no mixed p.
    auto *const pmesh = fes.GetParMesh();

    MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> coef(mat_op);
    mfem::CurlCurlIntegrator flux_integrator(coef);  // L2ZZErrorEstimator ignores coef

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.real(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               real_error, normp));

    static_cast<void>(mfem::L2ZZErrorEstimator(flux_integrator, field.imag(),
                                               smooth_flux_fes.GetFinestFESpace(), flux_fes,
                                               imag_error, normp));
  }
  else
  {
    // This lambda computes the flux as achieved within mfem. Use this to bench against.
    auto mfem_flux_func = [this, &field]()
    {
      mfem::ParComplexGridFunction flux_func(&flux_fes);
      flux_func.real() = 0.0;
      flux_func.imag() = 0.0;

      mfem::Array<int> xdofs, fdofs;
      mfem::Vector el_x, el_f;
      mfem::CurlCurlIntegrator curl;

      for (int i = 0; i < fes.GetNE(); ++i)
      {
        const auto * const xdoftrans = fes.GetElementVDofs(i, xdofs);
        field.real().GetSubVector(xdofs, el_x);
        if (xdoftrans)
        {
          xdoftrans->InvTransformPrimal(el_x);
        }

        auto *T = field.real().ParFESpace()->GetElementTransformation(i);
        curl.ComputeElementFlux(*field.real().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);

        const auto * const fdoftrans = flux_fes.GetElementVDofs(i, fdofs);
        if (fdoftrans)
        {
          fdoftrans->TransformPrimal(el_f);
        }

        flux_func.real().SetSubVector(fdofs, el_f);

        field.imag().GetSubVector(xdofs, el_x);
        if (xdoftrans)
        {
          xdoftrans->InvTransformPrimal(el_x);
        }
        curl.ComputeElementFlux(*field.imag().ParFESpace()->GetFE(i), *T, el_x,
                                *flux_fes.GetFE(i), el_f, false);
        if (fdoftrans)
        {
          fdoftrans->TransformPrimal(el_f);
        }
        flux_func.imag().SetSubVector(fdofs, el_f);
      }

      return flux_func;
    };

    // This lambda computes the curl flux function by forming a projector. In
    // theory same as the above, but still doesn't allow for coefficients. This
    // is another benching device, helpful for debugging.
    auto projector_flux_func = [this, &cv, &field]()
    {
      // Interpolate the weighted curl of the field in fes onto the discrete
      // flux space.
      mfem::ParDiscreteLinearOperator curl(&fes, &flux_fes); // (domain, range)

      auto muinv = MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY>(mat_op);
      curl.AddDomainInterpolator(new mfem::CurlInterpolator);
      curl.Assemble();
      curl.Finalize();
      auto hyprecurlop = std::unique_ptr<mfem::HypreParMatrix>(curl.ParallelAssemble());

      ComplexVector flux(flux_fes.GetTrueVSize());

      hyprecurlop->Mult(cv.real, flux.real);
      hyprecurlop->Mult(cv.imag, flux.imag);

      mfem::ParComplexGridFunction flux_func(&flux_fes);
      flux_func.real().SetFromTrueDofs(flux.real);
      flux_func.imag().SetFromTrueDofs(flux.imag);

      return flux_func;
    };

    // Given a flux function built elsewhere, construct the linear operator RHS.
    auto flux_rhs_from_func = [&v](auto &fes, mfem::ParComplexGridFunction &flux_func)
    {
      const auto ndof = fes.GetTrueVSize();

      ComplexVector flux(ndof);
      {
        mfem::VectorGridFunctionCoefficient f_real(&flux_func.real());
        mfem::ParLinearForm rhs(&fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_real));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.real);
      }

      {
        mfem::VectorGridFunctionCoefficient f_imag(&flux_func.imag());
        mfem::ParLinearForm rhs(&fes);
        rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(f_imag));
        rhs.UseFastAssembly(true);
        rhs.Assemble();
        rhs.ParallelAssemble(flux.imag);
      }

      return flux;
    };

    // Compare the two flux functions, should be identical.
    auto comp = [this](const auto &mflux, const auto &pflux)
    {
      bool match = true;

      mfem::Array<int> flux_dofs;
      mfem::Vector mflux_val, pflux_val;
      for (int e = 0; e < flux_fes.GetNE(); e++)
      {
        flux_fes.GetElementVDofs(e, flux_dofs);

        pflux.GetSubVector(flux_dofs, pflux_val);
        mflux.GetSubVector(flux_dofs, mflux_val);

        MFEM_ASSERT(mflux_val.Size() == pflux_val.Size(), "Must match size.");
        constexpr double tol = 1e-6;
        for (int i = 0; i < mflux_val.Size(); ++i)
        {
          auto diff = std::abs(mflux_val(i) - pflux_val(i));
          if (diff > tol)
          {
            std::cout << "Mismatch on e " << e << " i " << i << ": " << mflux_val(i) << " " << pflux_val(i) << '\n';
            match = false;
          }
        }
      }

      return match;
    };

    MFEM_ASSERT(comp(mfem_flux_func().real(), projector_flux_func().real()), "Mismatch between projector and L2ZZ construction real values");
    MFEM_ASSERT(comp(mfem_flux_func().imag(), projector_flux_func().imag()), "Mismatch between projector and L2ZZ construction imag values");

    // Coefficients for computing the discontinuous flux., i.e. (W, μ⁻¹∇ × V).
    // The code from here down will ultimately be the way to calculate the flux.
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

    // Switching between these two gives identical.
    // auto flux_func = projector_flux_func();
    // const auto flux = flux_rhs_from_func(smooth_flux_fes.GetFinestFESpace(), flux_func);
    const auto flux = ComplexVector(rhs_from_coef(smooth_flux_fes.GetFinestFESpace(), real_coef),
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

    // Specify a minimum integration rule. Using normp = 2 above ensures the
    // element integrands are polynomials of O(2p + q), thus we can exactly
    // integrate using this minimum rule.
    std::vector<const mfem::IntegrationRule *> irs;irs.reserve(nelem);
    for (int i = 0; i < smooth_flux_fes.GetFinestFESpace().GetNE(); ++i)
    {
      const int order = 2 * smooth_flux_fes.GetFinestFESpace().GetOrder(i)
                      + smooth_flux_fes.GetFinestFESpace().GetParMesh()->GetElementTransformation(i)->OrderW();
      irs.push_back(&mfem::IntRules.Get(smooth_flux_fes.GetFinestFESpace().GetParMesh()->GetElementGeometry(i), order));
    }

    smooth_flux_func.real().ComputeElementLpErrors(normp, real_coef, real_error, nullptr, nullptr, irs.data());
    smooth_flux_func.imag().ComputeElementLpErrors(normp, imag_coef, imag_error, nullptr, nullptr, irs.data());
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

    for (int i = 0; i < nelem; ++i)
    {
      estimates[i] = mfem::ComputeElementLpDistance(normp, i, flux_func, smooth_flux_func);
    }

    return estimates;
  }
}

}  // namespace palace
