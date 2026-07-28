// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/multigrid.hpp"
#include "fem/singulardofs.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

std::unique_ptr<mfem::HypreParMatrix> AddMatrices(double alpha,
                                                  const mfem::HypreParMatrix &a,
                                                  double beta,
                                                  const mfem::HypreParMatrix &b)
{
  return std::unique_ptr<mfem::HypreParMatrix>(mfem::Add(alpha, a, beta, b));
}

fem::singular::ParallelSparseOperatorBlocks
CombineOperatorBlocks(double alpha, const fem::singular::ParallelSparseOperatorBlocks &a,
                      double beta, const fem::singular::ParallelSparseOperatorBlocks &b)
{
  MFEM_VERIFY(a.enrichment_enrichment && a.standard_enrichment && a.enrichment_standard &&
                  b.enrichment_enrichment && b.standard_enrichment && b.enrichment_standard,
              "Cannot combine incomplete singular operator blocks!");
  fem::singular::ParallelSparseOperatorBlocks result;
  result.enrichment_enrichment =
      AddMatrices(alpha, *a.enrichment_enrichment, beta, *b.enrichment_enrichment);
  result.standard_enrichment =
      AddMatrices(alpha, *a.standard_enrichment, beta, *b.standard_enrichment);
  result.enrichment_standard =
      AddMatrices(alpha, *a.enrichment_standard, beta, *b.enrichment_standard);
  return result;
}

fem::singular::ParallelSparseOperatorBlocks
CombineOperatorBlocks(double alpha, const fem::singular::ParallelSparseOperatorBlocks &a,
                      double beta, const fem::singular::ParallelSparseOperatorBlocks &b,
                      double gamma, const fem::singular::ParallelSparseOperatorBlocks &c)
{
  auto ab = CombineOperatorBlocks(alpha, a, beta, b);
  return CombineOperatorBlocks(1.0, ab, gamma, c);
}

}  // namespace

BoundaryModeOperator::BoundaryModeOperator(
    const IoData &iodata_, const std::vector<std::unique_ptr<Mesh>> &mesh,
    const MaterialOperator &mat_op_in,
    const fem::singular::TriangleFeatureTopology *singular_features_in,
    const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids_in)
  : iodata(iodata_), solver_order(iodata_.solver.order), solve_mesh(mesh.back().get()),
    singular_features(singular_features_in), source_vertex_ids(source_vertex_ids_in),
    mat_op(mat_op_in)
{
  MFEM_VERIFY(solve_mesh->Dimension() == 2,
              "BoundaryMode solver requires a 2D mesh (waveguide cross-section). When "
              "\"Attributes\" is set the driver's Preprocess hook extracts it from "
              "the 3D parent before this operator is constructed.");

  // FE spaces.
  SetUpFESpaces(mesh);

  // Boundary operators.
  auto &pmesh = solve_mesh->Get();
  surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, mat_op, pmesh);
  farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, mat_op, pmesh);
  surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(iodata, mat_op, pmesh);
  surf_rz_op = std::make_unique<SurfaceRationalImpedanceOperator>(iodata, mat_op, pmesh);

  if (HasSingularEnrichment())
  {
    SetUpSingularEnrichment();
  }

  // Frequency-independent block matrices. Atn is the ND/H1 gradient coupling,
  // Btn = -Atn^T, Btt is the ND mass.
  if (!HasSingularEnrichment())
  {
    std::tie(Atnr, Atni) = mode_assembly::AssembleAtn(GetNDSpace(), GetH1Space(), mat_op);
  }
  auto [bttr_tmp, btti_tmp] = mode_assembly::AssembleBtt(GetNDSpace(), mat_op);
  MFEM_VERIFY(!btti_tmp, "BoundaryMode inverse-permeability mass must be real!");
  if (HasSingularEnrichment())
  {
    Bttr = fem::singular::BuildParallelEnrichedOperator(*bttr_tmp,
                                                        singular_mu_matrices.nd_mass);
    Atnr.reset(mfem::ParMult(Bttr.get(), singular_gradient.get(), true));
    MFEM_VERIFY(Atnr, "Failed to multiply the singular ND mass and exact gradient!");
    *Atnr *= -1.0;
  }
  else
  {
    Bttr = std::move(bttr_tmp);
  }
  Btnr.reset(Atnr->Transpose());
  *Btnr *= -1.0;
  if (Atni)
  {
    Btni.reset(Atni->Transpose());
    *Btni *= -1.0;
  }
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             GetNDGlobalTrueVSize(), GetH1GlobalTrueVSize(),
             GetNDGlobalTrueVSize() + GetH1GlobalTrueVSize());
}

BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAtt(std::complex<double> omega, double sigma) const
{
  auto [standard_real, standard_imag] =
      mode_assembly::AssembleAtt(GetNDSpace(), mat_op, nullptr, *surf_z_op, *farfield_op,
                                 *surf_sigma_op, *surf_rz_op, omega, sigma);
  if (!HasSingularEnrichment())
  {
    return {std::move(standard_real), std::move(standard_imag)};
  }
  MFEM_VERIFY(omega.imag() == 0.0,
              "Triangular BoundaryMode singular enrichment currently supports only "
              "real excitation frequencies!");
  const double omega_squared = omega.real() * omega.real();
  auto enrichment_real = CombineOperatorBlocks(
      1.0, singular_mu_matrices.nd_curl_curl, -omega_squared,
      singular_epsilon_matrices.nd_mass, -sigma, singular_mu_matrices.nd_mass);
  auto combined_real =
      fem::singular::BuildParallelEnrichedOperator(*standard_real, enrichment_real);
  std::unique_ptr<mfem::HypreParMatrix> combined_imaginary;
  if (standard_imag)
  {
    MFEM_VERIFY(singular_epsilon_imag_matrices.nd_mass.enrichment_enrichment,
                "Triangular BoundaryMode imaginary permittivity blocks were not "
                "assembled!");
    auto enrichment_imaginary =
        CombineOperatorBlocks(-omega_squared, singular_epsilon_imag_matrices.nd_mass, 0.0,
                              singular_epsilon_imag_matrices.nd_mass);
    combined_imaginary =
        fem::singular::BuildParallelEnrichedOperator(*standard_imag, enrichment_imaginary);
  }
  return {std::move(combined_real), std::move(combined_imaginary)};
}

BoundaryModeOperator::ComplexHypreParMatrix
BoundaryModeOperator::AssembleAnn(std::complex<double> omega) const
{
  auto [standard_real, standard_imag] =
      mode_assembly::AssembleAnn(GetH1Space(), mat_op, nullptr, *surf_z_op, *farfield_op,
                                 *surf_sigma_op, *surf_rz_op, omega);
  if (!HasSingularEnrichment())
  {
    return {std::move(standard_real), std::move(standard_imag)};
  }
  MFEM_VERIFY(omega.imag() == 0.0,
              "Triangular BoundaryMode singular enrichment currently supports only "
              "real excitation frequencies!");
  const double omega_squared = omega.real() * omega.real();
  auto enrichment_real =
      CombineOperatorBlocks(-1.0, singular_mu_matrices.h1_diffusion, omega_squared,
                            singular_epsilon_matrices.h1_mass);
  auto combined_real =
      fem::singular::BuildParallelEnrichedOperator(*standard_real, enrichment_real);
  std::unique_ptr<mfem::HypreParMatrix> combined_imaginary;
  if (standard_imag)
  {
    MFEM_VERIFY(singular_epsilon_imag_matrices.h1_mass.enrichment_enrichment,
                "Triangular BoundaryMode imaginary permittivity blocks were not "
                "assembled!");
    auto enrichment_imaginary =
        CombineOperatorBlocks(omega_squared, singular_epsilon_imag_matrices.h1_mass, 0.0,
                              singular_epsilon_imag_matrices.h1_mass);
    combined_imaginary =
        fem::singular::BuildParallelEnrichedOperator(*standard_imag, enrichment_imaginary);
  }
  return {std::move(combined_real), std::move(combined_imaginary)};
}

void BoundaryModeOperator::ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn,
                                                ComplexVector &et, ComplexVector &en) const
{
  mode_assembly::ApplyVDBackTransform(e0, kn, GetNDTrueVSize(), GetH1TrueVSize(), et, en);
}

std::complex<double>
BoundaryModeOperator::ComputePoyntingPower(double omega, std::complex<double> kn,
                                           const ComplexVector &et,
                                           const ComplexVector &en) const
{
  if (!Bttr)
  {
    return 0.0;
  }
  auto comm = GetComm();
  std::complex<double> P = 0.5 * std::conj(kn) / omega * linalg::Dot(comm, et, *Bttr, et);
  if (Atnr && en.Size() == GetH1TrueVSize())
  {
    ComplexWrapperOperator Atn(const_cast<mfem::HypreParMatrix *>(Atnr.get()),
                               const_cast<mfem::HypreParMatrix *>(Atni.get()));
    P += std::complex<double>(0.0, 1.0) / (2.0 * omega) * linalg::Dot(comm, en, Atn, et);
  }
  return P;
}

BoundaryModeOperator::SingularCoefficientNorms
BoundaryModeOperator::ComputeSingularCoefficientNorms(const ComplexVector &et,
                                                      const ComplexVector &en) const
{
  MFEM_VERIFY(HasSingularEnrichment() && singular_dofs && singular_numbering &&
                  et.Size() == GetNDTrueVSize() && en.Size() == GetH1TrueVSize(),
              "Singular BoundaryMode coefficient diagnostics received inconsistent "
              "topology or vector dimensions!");
  const int rank = Mpi::Rank(GetComm());
  const int standard_nd_size = GetNDSpace().GetTrueVSize();
  const int standard_h1_size = GetH1Space().GetTrueVSize();
  std::array<double, 3> squared_norms{0.0, 0.0, 0.0};
  for (std::size_t local = 0; local < singular_dofs->nd_dofs.size(); local++)
  {
    if (singular_numbering->nd.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt owned =
        singular_numbering->nd.local_to_true[local] - singular_numbering->nd.owned_offset;
    MFEM_VERIFY(owned >= 0 && owned < singular_numbering->nd.owned_size,
                "Owned singular ND true DOF is outside the local owner block!");
    const int index = standard_nd_size + static_cast<int>(owned);
    const double magnitude_squared =
        et.Real()[index] * et.Real()[index] + et.Imag()[index] * et.Imag()[index];
    const auto family = singular_dofs->nd_dofs[local].family;
    if (family == fem::singular::HigherOrderBasisFamily::NODE_GRADIENT)
    {
      squared_norms[0] += magnitude_squared;
    }
    else
    {
      MFEM_VERIFY(family == fem::singular::HigherOrderBasisFamily::NODE_ROTATIONAL,
                  "Triangular singular ND diagnostics found an unsupported basis "
                  "family!");
      squared_norms[1] += magnitude_squared;
    }
  }
  for (std::size_t local = 0; local < singular_dofs->h1_dofs.size(); local++)
  {
    if (singular_numbering->h1.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt owned =
        singular_numbering->h1.local_to_true[local] - singular_numbering->h1.owned_offset;
    MFEM_VERIFY(owned >= 0 && owned < singular_numbering->h1.owned_size &&
                    singular_dofs->h1_dofs[local].family ==
                        fem::singular::HigherOrderBasisFamily::NODE_GRADIENT,
                "Owned singular H1 true DOF is inconsistent with the triangular "
                "gradient topology!");
    const int index = standard_h1_size + static_cast<int>(owned);
    squared_norms[2] +=
        en.Real()[index] * en.Real()[index] + en.Imag()[index] * en.Imag()[index];
  }
  Mpi::GlobalSum(static_cast<int>(squared_norms.size()), squared_norms.data(), GetComm());
  return {std::sqrt(squared_norms[0]), std::sqrt(squared_norms[1]),
          std::sqrt(squared_norms[2])};
}

BoundaryModeOperator::SingularFieldEnergies
BoundaryModeOperator::ComputeSingularFieldEnergies(double omega, std::complex<double> kn,
                                                   const ComplexVector &et,
                                                   const ComplexVector &en) const
{
  MFEM_VERIFY(HasSingularEnrichment() && singular_epsilon_nd_mass &&
                  singular_epsilon_h1_mass && singular_mu_nd_curl_curl &&
                  std::isfinite(omega) && omega > 0.0 && std::isfinite(kn.real()) &&
                  std::isfinite(kn.imag()) && et.Size() == GetNDTrueVSize() &&
                  en.Size() == GetH1TrueVSize(),
              "Singular BoundaryMode energy calculation received inconsistent "
              "operators, frequency, or vector dimensions!");

  const auto quadratic_energy =
      [this](const ComplexVector &field, const mfem::HypreParMatrix &matrix, double scale)
  {
    Vector work(matrix.Height());
    work.UseDevice(true);
    matrix.Mult(field.Real(), work);
    const double real_product = linalg::Dot(GetComm(), field.Real(), work);
    matrix.Mult(field.Imag(), work);
    const double imaginary_product = linalg::Dot(GetComm(), field.Imag(), work);
    const double tolerance =
        1.0e-10 * std::max({1.0, std::abs(real_product), std::abs(imaginary_product)});
    MFEM_VERIFY(std::isfinite(real_product) && std::isfinite(imaginary_product) &&
                    real_product >= -tolerance && imaginary_product >= -tolerance &&
                    std::isfinite(scale) && scale > 0.0,
                "Singular BoundaryMode energy quadratic form is not finite and "
                "nonnegative!");
    return scale * (std::max(0.0, real_product) + std::max(0.0, imaginary_product));
  };

  ComplexVector grad_en(GetNDTrueVSize());
  singular_gradient->Mult(en.Real(), grad_en.Real());
  singular_gradient->Mult(en.Imag(), grad_en.Imag());
  ComplexVector transverse_magnetic_coordinate(et);
  // Bt = R(-kn Et + i grad(En))/omega, where R is a ninety-degree
  // rotation. Isotropic mu^{-1} makes R an isometry for the Btt inner product.
  transverse_magnetic_coordinate *= -kn;
  transverse_magnetic_coordinate.AXPY(std::complex<double>(0.0, 1.0), grad_en);

  return {quadratic_energy(et, *singular_epsilon_nd_mass, 0.5),
          quadratic_energy(en, *singular_epsilon_h1_mass, 0.5),
          quadratic_energy(transverse_magnetic_coordinate, *Bttr, 0.5 / (omega * omega)),
          quadratic_energy(et, *singular_mu_nd_curl_curl, 0.5 / (omega * omega))};
}

void BoundaryModeOperator::SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  const auto &mg = iodata.solver.linear;
  const int dim = solve_mesh->Dimension();

  // Dirichlet bdr attrs. The solve mesh already carries PEC on every edge that needs it:
  // for direct 2D, the user defined it; for submesh extraction, the driver folded the
  // inherited 3D boundary conditions into 2D attributes before this runs.
  {
    const auto &pmesh = solve_mesh->Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    for (auto attr : iodata.boundaries.auxpec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    dbc_bcs.Sort();
    dbc_bcs.Unique();
  }

  // FE collections.
  nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_aux_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  // RT collection for estimator flux recovery (order = solver_order - 1 since the
  // recovered flux lives one polynomial degree below the primary ND space).
  const int rt_mg_max_levels = mg.estimator_mg ? mg.mg_max_levels : 1;
  rt_fecs = fem::ConstructFECollections<mfem::RT_FECollection>(
      solver_order - 1, dim, rt_mg_max_levels, mg.mg_coarsening, false);

  nd_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
      mg.mg_max_levels, mesh, nd_fecs, &dbc_bcs, &nd_dbc_tdof_lists);
  h1_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, mesh, h1_fecs, &dbc_bcs, &h1_dbc_tdof_lists);
  h1_aux_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, mesh, h1_aux_fecs, &dbc_bcs, &h1_aux_dbc_tdof_lists);
  rt_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
      rt_mg_max_levels, mesh, rt_fecs);

  // L2 curl space for 2D B-field.
  l2_curl_fec = std::make_unique<mfem::L2_FECollection>(
      solver_order - 1, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL);
  l2_curl_fespace = std::make_unique<FiniteElementSpace>(*solve_mesh, l2_curl_fec.get());

  const int standard_nd_size = GetNDSpace().GetTrueVSize();
  combined_dbc_tdof_list.Append(nd_dbc_tdof_lists.back());
  for (int dof : h1_dbc_tdof_lists.back())
  {
    combined_dbc_tdof_list.Append(standard_nd_size + dof);
  }
}

void BoundaryModeOperator::SetUpSingularEnrichment()
{
  MFEM_VERIFY(singular_features && source_vertex_ids &&
                  source_vertex_ids->size() ==
                      static_cast<std::size_t>(solve_mesh->Get().GetNV()) &&
                  iodata.solver.singular_elements.order == 1 &&
                  iodata.solver.linear.mg_max_levels == 1,
              "BoundaryMode singular enrichment requires complete source topology, "
              "singular order one, and one finite-element level!");
  MFEM_VERIFY(!mat_op.HasConductivity() && !mat_op.HasLondonDepth(),
              "BoundaryMode singular enrichment does not support bulk conductivity or "
              "London penetration depth!");

  singular_dofs = std::make_unique<fem::singular::TriangleDofTopology>(
      fem::singular::BuildLocalTriangleDofTopology(solve_mesh->Get(), *singular_features,
                                                   *source_vertex_ids,
                                                   iodata.solver.singular_elements.order));
  singular_numbering = std::make_unique<fem::singular::ParallelDofNumbering>(
      fem::singular::BuildParallelDofNumbering(GetComm(), *singular_dofs));
  MFEM_VERIFY(singular_numbering->h1.owned_size <= std::numeric_limits<int>::max() &&
                  singular_numbering->nd.owned_size <= std::numeric_limits<int>::max(),
              "BoundaryMode singular local true-DOF count exceeds integer limits!");

  std::vector<fem::singular::IsotropicMaterialCoefficients> mu_materials(
      solve_mesh->GetNE(), {1.0, 1.0});
  std::vector<fem::singular::IsotropicMaterialCoefficients> epsilon_materials(
      solve_mesh->GetNE(), {1.0, 1.0});
  std::vector<fem::singular::IsotropicMaterialCoefficients> epsilon_imag_materials(
      solve_mesh->GetNE(), {0.0, 1.0});
  for (int element = 0; element < solve_mesh->GetNE(); element++)
  {
    const int attribute = solve_mesh->Get().GetAttribute(element);
    MFEM_VERIFY(mat_op.IsIsotropic(attribute),
                "Triangular BoundaryMode singular enrichment requires isotropic "
                "permittivity and permeability in every triangle! Domain attribute: "
                    << attribute);
    if (singular_features->elements[element].nodes.empty())
    {
      continue;
    }
    const double inverse_permeability = mat_op.GetInvPermeability(attribute)(0, 0);
    const double permittivity = mat_op.GetPermittivityReal(attribute)(0, 0);
    mu_materials[element] = {inverse_permeability, inverse_permeability};
    epsilon_materials[element] = {permittivity, 1.0};
    epsilon_imag_materials[element] = {mat_op.GetPermittivityImag(attribute)(0, 0), 1.0};
  }
  const fem::singular::AdaptiveAssemblyOptions options{
      iodata.solver.singular_elements.quadrature_order,
      iodata.solver.singular_elements.abs_tol, iodata.solver.singular_elements.rel_tol,
      iodata.solver.singular_elements.max_subdivisions};
  const auto local_mu = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      *singular_dofs, GetH1Space().Get(), GetNDSpace().Get(), mu_materials, options);
  const auto local_epsilon = fem::singular::AssembleLocalSparseEnrichmentMatrices(
      *singular_dofs, GetH1Space().Get(), GetNDSpace().Get(), epsilon_materials, options);
  singular_mu_matrices = fem::singular::AssembleParallelSparseEnrichmentMatrices(
      local_mu, *singular_numbering, GetH1Space().Get(), GetNDSpace().Get());
  singular_epsilon_matrices = fem::singular::AssembleParallelSparseEnrichmentMatrices(
      local_epsilon, *singular_numbering, GetH1Space().Get(), GetNDSpace().Get());
  if (mat_op.HasLossTangent())
  {
    const auto local_epsilon_imag = fem::singular::AssembleLocalSparseEnrichmentMatrices(
        *singular_dofs, GetH1Space().Get(), GetNDSpace().Get(), epsilon_imag_materials,
        options);
    singular_epsilon_imag_matrices =
        fem::singular::AssembleParallelSparseEnrichmentMatrices(
            local_epsilon_imag, *singular_numbering, GetH1Space().Get(),
            GetNDSpace().Get());
  }

  {
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal());
    BilinearForm mass(GetNDSpace());
    mass.AddDomainIntegrator<VectorFEMassIntegrator>(epsilon_func);
    auto standard =
        ParOperator(mass.FullAssemble(false), GetNDSpace()).StealParallelAssemble();
    singular_epsilon_nd_mass = fem::singular::BuildParallelEnrichedOperator(
        *standard, singular_epsilon_matrices.nd_mass);
  }
  {
    MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityScalar());
    BilinearForm mass(GetH1Space());
    mass.AddDomainIntegrator<MassIntegrator>(epsilon_func);
    auto standard =
        ParOperator(mass.FullAssemble(false), GetH1Space()).StealParallelAssemble();
    singular_epsilon_h1_mass = fem::singular::BuildParallelEnrichedOperator(
        *standard, singular_epsilon_matrices.h1_mass);
  }
  {
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetCurlCurlInvPermeability());
    BilinearForm stiffness(GetNDSpace());
    stiffness.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
    auto standard =
        ParOperator(stiffness.FullAssemble(false), GetNDSpace()).StealParallelAssemble();
    singular_mu_nd_curl_curl = fem::singular::BuildParallelEnrichedOperator(
        *standard, singular_mu_matrices.nd_curl_curl);
  }

  auto enrichment_gradient =
      fem::singular::BuildParallelEnrichmentGradient(GetComm(), *singular_numbering);
  const auto &standard_gradient_operator =
      GetNDSpace().GetDiscreteInterpolator(GetH1Space());
  const auto *standard_gradient =
      dynamic_cast<const ParOperator *>(&standard_gradient_operator);
  MFEM_VERIFY(standard_gradient,
              "BoundaryMode singular enrichment requires an assembled standard "
              "gradient!");
  singular_gradient = fem::singular::BuildParallelEnrichedGradient(
      standard_gradient->ParallelAssemble(), *enrichment_gradient);

  singular_h1_essential_true_dofs = fem::singular::GetEssentialTriangleH1TrueDofs(
      GetComm(), *singular_features, *singular_dofs, *singular_numbering);
  singular_nd_essential_true_dofs = fem::singular::GetEssentialTriangleNDTrueDofs(
      GetComm(), *singular_features, *singular_dofs, *singular_numbering);

  const int standard_nd_size = GetNDSpace().GetTrueVSize();
  const int combined_nd_size = GetNDTrueVSize();
  combined_dbc_tdof_list.DeleteAll();
  combined_dbc_tdof_list.Append(nd_dbc_tdof_lists.back());
  for (int dof : singular_nd_essential_true_dofs)
  {
    combined_dbc_tdof_list.Append(standard_nd_size + dof);
  }
  for (int dof : h1_dbc_tdof_lists.back())
  {
    combined_dbc_tdof_list.Append(combined_nd_size + dof);
  }
  const int standard_h1_size = GetH1Space().GetTrueVSize();
  for (int dof : singular_h1_essential_true_dofs)
  {
    combined_dbc_tdof_list.Append(combined_nd_size + standard_h1_size + dof);
  }
  combined_dbc_tdof_list.Sort();
  MFEM_VERIFY(
      std::adjacent_find(combined_dbc_tdof_list.begin(), combined_dbc_tdof_list.end()) ==
          combined_dbc_tdof_list.end(),
      "BoundaryMode singular essential true DOFs are not unique!");

  Mpi::Print(" Singular BoundaryMode enrichment: {:d} ND + {:d} H1 global true "
             "DOFs\n",
             singular_numbering->nd.global_size, singular_numbering->h1.global_size);
}

}  // namespace palace
