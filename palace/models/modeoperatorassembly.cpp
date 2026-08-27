// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeoperatorassembly.hpp"

#include <cmath>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"

namespace palace::mode_assembly
{

namespace
{
constexpr bool skip_zeros = false;

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

ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op)
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op)
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
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
    const MaterialOperator &mat_op, const mfem::Vector *normal,
    SurfaceImpedanceOperator &surf_z_op, FarfieldBoundaryOperator &farfield_op,
    SurfaceConductivityOperator &surf_sigma_op,
    SurfaceRationalImpedanceOperator &surf_rz_op, const mfem::HypreParMatrix &Bttr,
    const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
    const mfem::HypreParMatrix *Btnr, const mfem::Array<int> &dbc_tdof_list)
  : surf_sigma_op(surf_sigma_op), surf_rz_op(surf_rz_op)
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
    surf_z_op.AddStiffnessBdrCoefficients(1.0, surf_k_t);
    if (!surf_k_t.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(surf_k_t);
    }
    auto Attr = assemble_nd(att);

    MaterialPropertyCoefficient neg_muinv(attr_to_mat, mat_op.GetInvPermeability(), -1.0);
    if (normal)
    {
      neg_muinv.NormalProjectedCoefficient(*normal);
    }
    BilinearForm ann(h1_fespace);
    ann.AddDomainIntegrator<DiffusionIntegrator>(neg_muinv);
    MaterialPropertyCoefficient london_n(n_attr);
    if (mat_op.HasLondonDepth())
    {
      if (!normal)
      {
        london_n.AddCoefficient(attr_to_mat, mat_op.GetInvLondonDepthScalar());
      }
      else
      {
        const auto &ild = mat_op.GetInvLondonDepth();
        mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
        for (int k = 0; k < ild.SizeK(); k++)
        {
          ild_scalar(0, 0, k) = ild(0, 0, k);
        }
        london_n.AddCoefficient(attr_to_mat, ild_scalar);
        london_n.NormalProjectedCoefficient(*normal);
      }
      ann.AddDomainIntegrator<MassIntegrator>(london_n);
    }
    MaterialPropertyCoefficient surf_k_n(max_bdr_attr);
    surf_z_op.AddStiffnessBdrCoefficients(-1.0, surf_k_n);
    if (!surf_k_n.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(surf_k_n);
    }
    add_component(CoefficientType::CONSTANT, -1, std::move(Attr), nullptr, assemble_h1(ann),
                  nullptr);
  }

  // Linear omega: i times bulk conductivity and resistive/farfield damping.
  if (mat_op.HasConductivity() || surf_z_op.GetRsAttrList().Size() > 0 ||
      farfield_op.GetAttrList().Size() > 0)
  {
    BilinearForm atti(nd_fespace);
    MaterialPropertyCoefficient cond_t(n_attr);
    if (mat_op.HasConductivity())
    {
      cond_t.AddCoefficient(attr_to_mat, mat_op.GetConductivity());
      atti.AddDomainIntegrator<VectorFEMassIntegrator>(cond_t);
    }
    MaterialPropertyCoefficient damp_t(max_bdr_attr);
    surf_z_op.AddDampingBdrCoefficients(1.0, damp_t);
    farfield_op.AddDampingBdrCoefficients(1.0, damp_t);
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
    surf_z_op.AddDampingBdrCoefficients(-1.0, damp_n);
    if (farfield_op.GetAttrList().Size() > 0)
    {
      const auto &inv_z = mat_op.GetInvImpedance();
      const auto &bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
      for (auto attr : farfield_op.GetAttrList())
      {
        const int mat_idx =
            (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
        const double inv_z0 = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
        auto ceed_attrs = mat_op.GetCeedBdrAttributes(attr);
        if (ceed_attrs.Size() > 0)
        {
          damp_n.AddMaterialProperty(ceed_attrs, inv_z0, -1.0);
        }
      }
    }
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
    surf_z_op.AddMassBdrCoefficients(-1.0, surf_m_t);
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
    surf_z_op.AddMassBdrCoefficients(1.0, surf_m_n);
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

  auto shift = std::make_unique<mfem::HypreParMatrix>(Bttr);
  *shift *= -1.0;
  add_component(CoefficientType::SHIFT, -1, std::move(shift), nullptr, nullptr, nullptr);

  auto add_boundary_component =
      [&](CoefficientType type, int index, MaterialPropertyCoefficient &unit)
  {
    BilinearForm att(nd_fespace), ann(h1_fespace);
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(unit);
    auto Attr = assemble_nd(att);
    unit *= -1.0;
    ann.AddBoundaryIntegrator<MassIntegrator>(unit);
    add_component(type, index, std::move(Attr), nullptr, assemble_h1(ann), nullptr);
  };
  for (std::size_t g = 0; g < surf_sigma_op.Size(); g++)
  {
    if (surf_sigma_op.IsActive(g))
    {
      MaterialPropertyCoefficient unit(max_bdr_attr);
      surf_sigma_op.AddBoundaryMassBdrCoefficients(g, unit);
      add_boundary_component(CoefficientType::SURFACE_CONDUCTIVITY, static_cast<int>(g),
                             unit);
    }
  }
  for (int b = 0; b < surf_rz_op.GetNumBoundaries(); b++)
  {
    MaterialPropertyCoefficient unit(max_bdr_attr);
    surf_rz_op.AddUnitBdrCoefficients(b, unit);
    add_boundary_component(CoefficientType::RATIONAL_IMPEDANCE, b, unit);
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
    case CoefficientType::SURFACE_CONDUCTIVITY:
      return surf_sigma_op.EvaluateScalar(static_cast<std::size_t>(component.index), omega);
    case CoefficientType::RATIONAL_IMPEDANCE:
      return surf_rz_op.EvalRobinCoefficient(component.index,
                                             std::complex<double>(0.0, 1.0) * omega);
  }
  MFEM_ABORT("Unknown mode-operator coefficient type!");
  return 0.0;
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
    Ar->EliminateBC(dbc_tdof_list, diag_policy);
  }
  if (Ai)
  {
    Ai->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
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

}  // namespace palace::mode_assembly
