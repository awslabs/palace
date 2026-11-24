// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeoperatorassembly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceterms.hpp"
#include "utils/communication.hpp"

namespace palace::mode_assembly
{

namespace
{
constexpr bool skip_zeros = false;

void EliminateEssentialDofs(mfem::HypreParMatrix &A, const mfem::Array<int> &dbc_tdof_list,
                            Operator::DiagonalPolicy diag_policy)
{
  // HypreParMatrix::EliminateBC only finds the eliminated column entries of a row through
  // the transposed entry, so it assumes a structurally symmetric matrix. A component with
  // only one of the off-diagonal blocks (Atn without -Btn, or the reverse) is not, so
  // eliminate the rows and columns for a general sparsity pattern first and then set the
  // diagonal.
  MFEM_ASSERT(diag_policy != Operator::DIAG_KEEP,
              "Unsupported diagonal policy for block system matrix elimination!");
  std::unique_ptr<mfem::HypreParMatrix> Ae(A.EliminateRowsCols(dbc_tdof_list));
  A.EliminateBC(dbc_tdof_list, diag_policy);
}

void AddScaled(std::unique_ptr<mfem::HypreParMatrix> &sum, const mfem::HypreParMatrix *term,
               double coefficient)
{
  if (!term || coefficient == 0.0)
  {
    return;
  }
  if (!sum)
  {
    sum = std::make_unique<mfem::HypreParMatrix>(*term);
    *sum *= coefficient;
  }
  else
  {
    sum.reset(mfem::Add(1.0, *sum, coefficient, *term));
  }
}

void AddScaledComplex(std::unique_ptr<mfem::HypreParMatrix> &real,
                      std::unique_ptr<mfem::HypreParMatrix> &imag,
                      const mfem::HypreParMatrix *term_real,
                      const mfem::HypreParMatrix *term_imag,
                      std::complex<double> coefficient)
{
  AddScaled(real, term_real, coefficient.real());
  AddScaled(real, term_imag, -coefficient.imag());
  AddScaled(imag, term_real, coefficient.imag());
  AddScaled(imag, term_imag, coefficient.real());
}
}  // namespace

mfem::DenseMatrix GetTransverseInvPermeability(const mfem::DenseMatrix &muinv,
                                               const mfem::Vector *normal)
{
  const int dim = muinv.Height();
  mfem::DenseMatrix nu(dim);
  if (!normal)
  {
    // Rᵀ μ⁻¹ R for R = [0, -1; 1, 0].
    MFEM_VERIFY(dim == 2 && muinv.Width() == 2,
                "The transverse inverse permeability of a 2D cross-section needs the 2 x 2 "
                "in-plane inverse permeability!");
    nu(0, 0) = muinv(1, 1);
    nu(0, 1) = -muinv(1, 0);
    nu(1, 0) = -muinv(0, 1);
    nu(1, 1) = muinv(0, 0);
    return nu;
  }
  MFEM_VERIFY(
      dim == 3 && muinv.Width() == 3 && normal->Size() == 3,
      "The transverse inverse permeability of a cross-section in 3D needs the 3 x 3 "
      "inverse permeability and a normal!");
  mfem::Vector n(*normal);
  n /= n.Norml2();
  mfem::DenseMatrix N(3), M(3);  // N u = n × u
  N = 0.0;
  N(0, 1) = -n(2);
  N(0, 2) = n(1);
  N(1, 0) = n(2);
  N(1, 2) = -n(0);
  N(2, 0) = -n(1);
  N(2, 1) = n(0);
  mfem::Mult(muinv, N, M);
  mfem::MultAtB(N, M, nu);
  return nu;
}

mfem::DenseTensor GetTransverseInvPermeability(const mfem::DenseTensor &muinv,
                                               const mfem::Vector *normal)
{
  mfem::DenseTensor nu(muinv.SizeI(), muinv.SizeJ(), muinv.SizeK());
  for (int k = 0; k < muinv.SizeK(); k++)
  {
    nu(k) = GetTransverseInvPermeability(muinv(k), normal);
  }
  return nu;
}

ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal)
{
  MaterialPropertyCoefficient nu_func(
      mat_op.GetAttributeToMaterial(),
      GetTransverseInvPermeability(mat_op.GetInvPermeability(), normal), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(nu_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal)
{
  MaterialPropertyCoefficient nu_func(
      mat_op.GetAttributeToMaterial(),
      GetTransverseInvPermeability(mat_op.GetInvPermeability(), normal));
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(nu_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

OperatorComponent::OperatorComponent(CoefficientType type, int index,
                                     std::unique_ptr<mfem::HypreParMatrix> &&Ar,
                                     std::unique_ptr<mfem::HypreParMatrix> &&Ai)
  : type(type), index(index), Ar(std::move(Ar)), Ai(std::move(Ai)),
    op(std::make_unique<ComplexWrapperOperator>(this->Ar.get(), this->Ai.get()))
{
}

ModeOperatorModel::ModeOperatorModel(
    const FiniteElementSpace &nd_fespace, const FiniteElementSpace &h1_fespace,
    const MaterialOperator &mat_op, const mfem::Vector *normal, SurfaceTerms &surface_terms,
    const mfem::HypreParMatrix &Bttr, const mfem::HypreParMatrix *Atnr,
    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Btnr,
    const mfem::Array<int> &dbc_tdof_list)
  : mat_op(mat_op), surface_terms(surface_terms)
{
  auto assemble_nd = [&](BilinearForm &form)
  {
    return ParOperator(form.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
  };
  auto assemble_h1 = [&](BilinearForm &form)
  {
    return ParOperator(form.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
  };
  auto add_component = [&](CoefficientType type, int index,
                           std::unique_ptr<mfem::HypreParMatrix> Attr,
                           std::unique_ptr<mfem::HypreParMatrix> Atti,
                           std::unique_ptr<mfem::HypreParMatrix> Annr,
                           std::unique_ptr<mfem::HypreParMatrix> Anni)
  {
    const bool constant = type == CoefficientType::CONSTANT;
    const bool shift = type == CoefficientType::SHIFT;
    std::unique_ptr<mfem::HypreParMatrix> lower_left;
    if (shift)
    {
      lower_left = std::make_unique<mfem::HypreParMatrix>(*Btnr);
      *lower_left *= -1.0;
    }
    auto [Ar, Ai] = BuildSystemMatrixA(
        nd_fespace, h1_fespace, dbc_tdof_list, Attr.get(), Atti.get(),
        constant ? Atnr : nullptr, constant ? Atni : nullptr, Annr.get(), Anni.get(),
        lower_left.get(), constant ? Operator::DIAG_ONE : Operator::DIAG_ZERO);
    if (Ar || Ai)
    {
      components.emplace_back(type, index, std::move(Ar), std::move(Ai));
    }
  };

  const auto &attr_to_mat = mat_op.GetAttributeToMaterial();
  const int n_attr = attr_to_mat.Size();
  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();

  // Constant: curl-curl/negative diffusion, London mass, and surface inductance.
  {
    MaterialPropertyCoefficient muinv_cc(attr_to_mat,
                                         normal ? mat_op.GetInvPermeability()
                                                : mat_op.GetCurlCurlInvPermeability());
    if (normal)
    {
      muinv_cc.NormalProjectedCoefficient(*normal);
    }
    BilinearForm att(nd_fespace);
    att.AddDomainIntegrator<CurlCurlIntegrator>(muinv_cc);
    MaterialPropertyCoefficient london_t(n_attr);
    if (mat_op.HasLondonDepth())
    {
      london_t.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepth());
      att.AddDomainIntegrator<VectorFEMassIntegrator>(london_t);
    }
    MaterialPropertyCoefficient surf_k_t(max_bdr_attr);
    surface_terms.AddStiffnessBdrCoefficients(1.0, surf_k_t);
    if (!surf_k_t.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(surf_k_t);
    }
    auto Attr = assemble_nd(att);

    // The diffusion of the out-of-plane component comes from the in-plane part of the curl,
    // with ν, not the out-of-plane component of μ⁻¹.
    MaterialPropertyCoefficient neg_nu(
        attr_to_mat, GetTransverseInvPermeability(mat_op.GetInvPermeability(), normal),
        -1.0);
    BilinearForm ann(h1_fespace);
    ann.AddDomainIntegrator<DiffusionIntegrator>(neg_nu);
    // The H1 block is assembled with the opposite sign (negative diffusion), so the London
    // mass, like the surface inductance, enters with a negative coefficient.
    MaterialPropertyCoefficient london_n(n_attr);
    if (mat_op.HasLondonDepth())
    {
      if (!normal)
      {
        london_n.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepthScalar(), -1.0);
      }
      else
      {
        london_n.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepth(), -1.0);
        london_n.NormalProjectedCoefficient(*normal);
      }
      ann.AddDomainIntegrator<MassIntegrator>(london_n);
    }
    MaterialPropertyCoefficient surf_k_n(max_bdr_attr);
    surface_terms.AddStiffnessBdrCoefficients(-1.0, surf_k_n);
    if (!surf_k_n.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(surf_k_n);
    }
    add_component(CoefficientType::CONSTANT, -1, std::move(Attr), nullptr, assemble_h1(ann),
                  nullptr);
  }

  // Linear omega: i times bulk conductivity and resistive/farfield damping.
  if (mat_op.HasConductivity() || surface_terms.HasDamping())
  {
    BilinearForm atti(nd_fespace);
    MaterialPropertyCoefficient cond_t(n_attr);
    if (mat_op.HasConductivity())
    {
      cond_t.AddCoefficient(attr_to_mat, mat_op.GetConductivity());
      atti.AddDomainIntegrator<VectorFEMassIntegrator>(cond_t);
    }
    MaterialPropertyCoefficient damp_t(max_bdr_attr);
    surface_terms.AddDampingBdrCoefficients(1.0, damp_t);
    if (!damp_t.empty())
    {
      atti.AddBoundaryIntegrator<VectorFEMassIntegrator>(damp_t);
    }

    BilinearForm anni(h1_fespace);
    MaterialPropertyCoefficient cond_n(n_attr);
    if (mat_op.HasConductivity())
    {
      cond_n.AddCoefficient(
          attr_to_mat, normal ? mat_op.GetConductivity() : mat_op.GetConductivityScalar(),
          -1.0);
      if (normal)
      {
        cond_n.NormalProjectedCoefficient(*normal);
      }
      anni.AddDomainIntegrator<MassIntegrator>(cond_n);
    }
    MaterialPropertyCoefficient damp_n(max_bdr_attr);
    surface_terms.AddScalarDampingBdrCoefficients(-1.0, damp_n);
    if (!damp_n.empty())
    {
      anni.AddBoundaryIntegrator<MassIntegrator>(damp_n);
    }
    add_component(CoefficientType::OMEGA, -1, nullptr, assemble_nd(atti), nullptr,
                  assemble_h1(anni));
  }

  // Quadratic omega: signed permittivity/loss and surface capacitance.
  {
    BilinearForm attr(nd_fespace), atti(nd_fespace);
    MaterialPropertyCoefficient neg_eps(attr_to_mat, mat_op.GetPermittivityReal(), -1.0);
    attr.AddDomainIntegrator<VectorFEMassIntegrator>(neg_eps);
    MaterialPropertyCoefficient surf_m_t(max_bdr_attr);
    surface_terms.AddMassBdrCoefficients(-1.0, surf_m_t);
    if (!surf_m_t.empty())
    {
      attr.AddBoundaryIntegrator<VectorFEMassIntegrator>(surf_m_t);
    }
    std::unique_ptr<mfem::HypreParMatrix> Atti;
    if (mat_op.HasLossTangent())
    {
      MaterialPropertyCoefficient neg_eps_i(attr_to_mat, mat_op.GetPermittivityImag(),
                                            -1.0);
      atti.AddDomainIntegrator<VectorFEMassIntegrator>(neg_eps_i);
      Atti = assemble_nd(atti);
    }

    BilinearForm annr(h1_fespace), anni(h1_fespace);
    MaterialPropertyCoefficient eps_n(attr_to_mat, normal ? mat_op.GetPermittivityReal()
                                                          : mat_op.GetPermittivityScalar());
    if (normal)
    {
      eps_n.NormalProjectedCoefficient(*normal);
    }
    annr.AddDomainIntegrator<MassIntegrator>(eps_n);
    MaterialPropertyCoefficient surf_m_n(max_bdr_attr);
    surface_terms.AddMassBdrCoefficients(1.0, surf_m_n);
    if (!surf_m_n.empty())
    {
      annr.AddBoundaryIntegrator<MassIntegrator>(surf_m_n);
    }
    std::unique_ptr<mfem::HypreParMatrix> Anni;
    if (mat_op.HasLossTangent())
    {
      MaterialPropertyCoefficient eps_i_n(attr_to_mat,
                                          normal ? mat_op.GetPermittivityImag()
                                                 : mat_op.GetPermittivityImagScalar());
      if (normal)
      {
        eps_i_n.NormalProjectedCoefficient(*normal);
      }
      anni.AddDomainIntegrator<MassIntegrator>(eps_i_n);
      Anni = assemble_h1(anni);
    }
    add_component(CoefficientType::OMEGA_SQUARED, -1, assemble_nd(attr), std::move(Atti),
                  assemble_h1(annr), std::move(Anni));
  }

  // Frequency-dependent permittivity beyond ε∞: the ω² ε terms above gain ω² χ(iω) = −g(iω)
  // on each dispersive material, with g = s² χ(s) the nonlinear 3D system coefficient. One
  // unit mass per material, like the 3D A2, since its terms share one coefficient. The
  // material is isotropic, so the normal projection of χ in the H1 block is χ itself.
  for (std::size_t i = 0; i < mat_op.NumFrequencyDependentPermittivityMaterials(); i++)
  {
    if (!mat_op.HasFrequencyDependentPermittivityA2(i) ||
        !mat_op.HasFrequencyDependentPermittivitySupport(i))
    {
      continue;
    }
    const auto attrs =
        mat_op.GetCeedAttributes(mat_op.GetFrequencyDependentPermittivityAttributes(i));
    MaterialPropertyCoefficient unit_t(n_attr), unit_n(n_attr);
    unit_t.AddMaterialProperty(attrs, 1.0);
    unit_n.AddMaterialProperty(attrs, -1.0);
    BilinearForm att(nd_fespace), ann(h1_fespace);
    if (!unit_t.empty())
    {
      att.AddDomainIntegrator<VectorFEMassIntegrator>(unit_t);
      ann.AddDomainIntegrator<MassIntegrator>(unit_n);
    }
    add_component(CoefficientType::FREQUENCY_DEPENDENT_PERMITTIVITY, static_cast<int>(i),
                  assemble_nd(att), nullptr, assemble_h1(ann), nullptr);
  }

  auto shift = std::make_unique<mfem::HypreParMatrix>(Bttr);
  *shift *= -1.0;
  add_component(CoefficientType::SHIFT, -1, std::move(shift), nullptr, nullptr, nullptr);

  auto add_boundary_component = [&](CoefficientType type, int index,
                                    MaterialPropertyCoefficient &unit_t,
                                    MaterialPropertyCoefficient &unit_n)
  {
    // Processes without elements on the boundary have an empty coefficient and assemble
    // their (zero) part of the component without the integrators.
    BilinearForm att(nd_fespace), ann(h1_fespace);
    if (!unit_t.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(unit_t);
    }
    auto Attr = assemble_nd(att);
    unit_n *= -1.0;
    if (!unit_n.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(unit_n);
    }
    add_component(type, index, std::move(Attr), nullptr, assemble_h1(ann), nullptr);
  };
  for (int i = 0; i < surface_terms.GetNumExtraSystemBdrMassTerms(); i++)
  {
    MaterialPropertyCoefficient unit_t(max_bdr_attr), unit_n(max_bdr_attr);
    surface_terms.AddExtraSystemBdrMassTermCoefficients(i, unit_t);
    surface_terms.AddScalarExtraSystemBdrMassTermCoefficients(i, unit_n);
    add_boundary_component(CoefficientType::SURFACE_TERM, i, unit_t, unit_n);
  }

  // Second-order absorbing boundary terms, the components -[0, G; 0, D] and -[M, 0; Gᵀ, 0]
  // (see ModeOperatorModel), for the terms with boundary elements on any process.
  for (int i = 0; i < surface_terms.GetNumSecondOrderFarfieldTerms(); i++)
  {
    MaterialPropertyCoefficient unit(max_bdr_attr);
    surface_terms.AddSecondOrderFarfieldTermCoefficients(i, unit, -1.0);
    bool any = !unit.empty();
    Mpi::GlobalOr(1, &any, nd_fespace.GetComm());
    if (!any)
    {
      continue;
    }
    BilinearForm m(nd_fespace), d(h1_fespace), g(h1_fespace, nd_fespace);
    if (!unit.empty())
    {
      m.AddBoundaryIntegrator<VectorFEMassIntegrator>(unit);
      d.AddBoundaryIntegrator<DiffusionIntegrator>(unit);
      g.AddBoundaryIntegrator<MixedVectorGradientIntegrator>(unit);
    }
    auto M = assemble_nd(m);
    auto D = assemble_h1(d);
    auto G = ParOperator(g.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
                 .StealParallelAssemble();
    std::unique_ptr<mfem::HypreParMatrix> Gt(G->Transpose());
    auto [Cr, Ci] =
        BuildSystemMatrixA(nd_fespace, h1_fespace, dbc_tdof_list, nullptr, nullptr, G.get(),
                           nullptr, D.get(), nullptr, nullptr, Operator::DIAG_ZERO);
    components.emplace_back(CoefficientType::SECOND_ORDER_FARFIELD, i, std::move(Cr),
                            std::move(Ci));
    auto [Sr, Si] =
        BuildSystemMatrixA(nd_fespace, h1_fespace, dbc_tdof_list, M.get(), nullptr, nullptr,
                           nullptr, nullptr, nullptr, Gt.get(), Operator::DIAG_ZERO);
    components.emplace_back(CoefficientType::SECOND_ORDER_FARFIELD_SHIFT, i, std::move(Sr),
                            std::move(Si));
    farfield_lower_blocks.emplace_back(i, std::move(Gt));
  }
}

void ModeOperatorModel::AddSecondOrderFarfieldLowerBlock(
    std::complex<double> omega, double sigma, std::unique_ptr<mfem::HypreParMatrix> &real,
    std::unique_ptr<mfem::HypreParMatrix> &imag) const
{
  for (const auto &[i, block] : farfield_lower_blocks)
  {
    AddScaledComplex(real, imag, block.get(), nullptr,
                     sigma * surface_terms.GetSecondOrderFarfieldTermScale(i, omega));
  }
}

std::complex<double>
ModeOperatorModel::EvaluateCoefficient(const OperatorComponent &component,
                                       std::complex<double> omega, double sigma) const
{
  switch (component.type)
  {
    case CoefficientType::CONSTANT:
      return 1.0;
    case CoefficientType::OMEGA:
      return omega;
    case CoefficientType::OMEGA_SQUARED:
      return omega * omega;
    case CoefficientType::SHIFT:
      return sigma;
    case CoefficientType::SURFACE_TERM:
      return surface_terms.GetExtraSystemBdrMassTermScale(component.index, omega);
    case CoefficientType::FREQUENCY_DEPENDENT_PERMITTIVITY:
      return mat_op.EvaluateFrequencyDependentPermittivityA2(
          static_cast<std::size_t>(component.index),
          std::complex<double>(0.0, 1.0) * omega);
    case CoefficientType::SECOND_ORDER_FARFIELD:
      return surface_terms.GetSecondOrderFarfieldTermScale(component.index, omega);
    case CoefficientType::SECOND_ORDER_FARFIELD_SHIFT:
      return sigma * surface_terms.GetSecondOrderFarfieldTermScale(component.index, omega);
  }
  MFEM_ABORT("Unknown mode-operator coefficient type!");
  return 0.0;
}

std::complex<double>
ModeOperatorModel::EvaluateBCoefficient(const OperatorComponent &component,
                                        std::complex<double> omega) const
{
  return (component.type == CoefficientType::SECOND_ORDER_FARFIELD_SHIFT)
             ? -surface_terms.GetSecondOrderFarfieldTermScale(component.index, omega)
             : 0.0;
}

bool ModeOperatorModel::HasFrequencyDependentB() const
{
  return std::ranges::any_of(
      components, [](const auto &component)
      { return component.type == CoefficientType::SECOND_ORDER_FARFIELD_SHIFT; });
}

std::unique_ptr<ComplexOperator>
ModeOperatorModel::AssembleB(const mfem::HypreParMatrix &B0,
                             std::complex<double> omega) const
{
  auto real = std::make_unique<mfem::HypreParMatrix>(B0);
  std::unique_ptr<mfem::HypreParMatrix> imag;
  for (const auto &component : components)
  {
    const auto coeff = EvaluateBCoefficient(component, omega);
    if (coeff != 0.0)
    {
      AddScaledComplex(real, imag, component.Ar.get(), component.Ai.get(), coeff);
    }
  }
  return std::make_unique<ComplexWrapperOperator>(std::move(real), std::move(imag));
}

std::unique_ptr<ComplexOperator> ModeOperatorModel::Assemble(std::complex<double> omega,
                                                             double sigma) const
{
  std::unique_ptr<mfem::HypreParMatrix> real, imag;
  for (const auto &component : components)
  {
    AddScaledComplex(real, imag, component.Ar.get(), component.Ai.get(),
                     EvaluateCoefficient(component, omega, sigma));
  }
  return std::make_unique<ComplexWrapperOperator>(std::move(real), std::move(imag));
}

ComplexHypreParMatrix BuildSystemMatrixA(
    const FiniteElementSpace &nd_fespace, const FiniteElementSpace &h1_fespace,
    const mfem::Array<int> &dbc_tdof_list, const mfem::HypreParMatrix *Attr,
    const mfem::HypreParMatrix *Atti, const mfem::HypreParMatrix *Atnr,
    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Annr,
    const mfem::HypreParMatrix *Anni, const mfem::HypreParMatrix *shifted_Btnr,
    Operator::DiagonalPolicy diag_policy)
{
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  Vector dtt(nd_fespace.GetTrueVSize()), dnn(h1_fespace.GetTrueVSize());
  dtt.UseDevice(false);
  dnn.UseDevice(false);
  dtt = 0.0;
  dnn = 0.0;
  auto diag_tt = std::make_unique<mfem::SparseMatrix>(dtt);
  auto diag_nn = std::make_unique<mfem::SparseMatrix>(dnn);
  auto Dtt_zero = std::make_unique<mfem::HypreParMatrix>(
      nd_fespace.Get().GetComm(), nd_fespace.Get().GlobalTrueVSize(),
      nd_fespace.Get().GetTrueDofOffsets(), diag_tt.get());
  auto Dnn_zero = std::make_unique<mfem::HypreParMatrix>(
      h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
      h1_fespace.Get().GetTrueDofOffsets(), diag_nn.get());

  std::unique_ptr<mfem::HypreParMatrix> Ar;
  if (Attr || Atnr || Annr || shifted_Btnr)
  {
    blocks(0, 0) = Attr ? Attr : Dtt_zero.get();
    blocks(0, 1) = Atnr;
    blocks(1, 0) = shifted_Btnr;
    blocks(1, 1) = Annr ? Annr : Dnn_zero.get();
    Ar.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti || Atni || Anni)
  {
    blocks(0, 0) = Atti ? Atti : Dtt_zero.get();
    blocks(0, 1) = Atni;
    blocks(1, 0) = nullptr;
    blocks(1, 1) = Anni ? Anni : Dnn_zero.get();
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }
  if (Ar)
  {
    EliminateEssentialDofs(*Ar, dbc_tdof_list, diag_policy);
  }
  if (Ai)
  {
    EliminateEssentialDofs(*Ai, dbc_tdof_list, Operator::DIAG_ZERO);
  }
  return {std::move(Ar), std::move(Ai)};
}

void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en)
{
  et.Real().MakeRef(e0.Real(), 0, nd_size);
  et.Imag().MakeRef(e0.Imag(), 0, nd_size);
  en.Real().MakeRef(e0.Real(), nd_size, h1_size);
  en.Imag().MakeRef(e0.Imag(), nd_size, h1_size);
  const auto ikn_inv = 1.0 / (std::complex<double>(0.0, 1.0) * kn);
  ComplexVector::AXPBY(ikn_inv, en.Real(), en.Imag(), 0.0, en.Real(), en.Imag());
}

std::complex<double> PropagationConstant(std::complex<double> kn2,
                                         std::complex<double> omega)
{
  // Real kₙ² at the real frequency Re{ω}, for kₙ² = n²ω² − k_c² with real n² and k_c².
  const double kn2_real = (omega.imag() != 0.0 && omega.real() != 0.0)
                              ? kn2.real() + 0.5 * kn2.imag() * omega.imag() / omega.real()
                              : kn2.real();
  return (kn2_real < 0.0) ? std::complex<double>(0.0, -1.0) * std::sqrt(-kn2)
                          : std::sqrt(kn2);
}

double TargetDistance(std::complex<double> kn, double kn_target)
{
  const double distance = std::abs(kn - kn_target);
  return std::isnan(distance) ? std::numeric_limits<double>::infinity() : distance;
}

}  // namespace palace::mode_assembly
