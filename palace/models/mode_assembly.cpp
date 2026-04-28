// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "mode_assembly.hpp"

#include <mfem.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"

namespace palace::mode_assembly
{

namespace
{
constexpr bool skip_zeros = false;
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

ComplexHypreParMatrix
AssembleAtt(const FiniteElementSpace &nd_fespace, const MaterialOperator &mat_op,
            const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
            FarfieldBoundaryOperator &farfield_op,
            SurfaceConductivityOperator &surf_sigma_op, double omega, double sigma)
{
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient eps_shifted_func(
      mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), -omega * omega);
  eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvPermeability(), -sigma);
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepth(), 1.0);
  }

  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);
  surf_z_op.AddStiffnessBdrCoefficients(1.0, fbr);
  surf_z_op.AddDampingBdrCoefficients(omega, fbi);
  surf_z_op.AddMassBdrCoefficients(-omega * omega, fbr);
  farfield_op.AddDampingBdrCoefficients(omega, fbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_func);
  if (!fbr.empty())
  {
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }
  auto Attr_assembled =
      ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();

  std::unique_ptr<mfem::HypreParMatrix> Atti_assembled;
  {
    const bool has_imag =
        mat_op.HasLossTangent() || mat_op.HasConductivity() || !fbi.empty();
    if (has_imag)
    {
      // Coefficients must outlive the BilinearForm (integrators hold raw pointers).
      const int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient negepstandelta_func(n_attr);
      MaterialPropertyCoefficient fi_domain(n_attr);
      if (mat_op.HasLossTangent())
      {
        negepstandelta_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityImag(), -omega * omega);
      }
      if (mat_op.HasConductivity())
      {
        fi_domain.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(),
                                 omega);
      }
      BilinearForm atti(nd_fespace);
      if (!negepstandelta_func.empty())
      {
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
      }
      if (!fi_domain.empty())
      {
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(fi_domain);
      }
      if (!fbi.empty())
      {
        atti.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbi);
      }
      Atti_assembled =
          ParOperator(atti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
    }
  }
  return {std::move(Attr_assembled), std::move(Atti_assembled)};
}

ComplexHypreParMatrix AssembleAnn(const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal,
                                  SurfaceImpedanceOperator &surf_z_op,
                                  FarfieldBoundaryOperator &farfield_op,
                                  SurfaceConductivityOperator &surf_sigma_op, double omega)
{
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             normal ? mat_op.GetPermittivityReal()
                                                    : mat_op.GetPermittivityScalar(),
                                             omega * omega);
  if (normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*normal);
  }
  if (mat_op.HasLondonDepth())
  {
    if (!normal)
    {
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepthScalar());
    }
    else
    {
      const auto &ild = mat_op.GetInvLondonDepth();
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(), ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*normal);
    }
  }

  const int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr), nn_fbi(max_bdr_attr);
  surf_z_op.AddStiffnessBdrCoefficients(-1.0, nn_fbr);
  surf_z_op.AddDampingBdrCoefficients(-omega, nn_fbi);
  surf_z_op.AddMassBdrCoefficients(omega * omega, nn_fbr);
  if (farfield_op.GetAttrList().Size() > 0)
  {
    // Farfield boundary: scalar inverse impedance for the H1 mass integrator.
    const auto &farfield_attrs = farfield_op.GetAttrList();
    const auto &inv_z = mat_op.GetInvImpedance();
    const auto &bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
    for (auto attr : farfield_attrs)
    {
      int mat_idx =
          (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
      double inv_z0_scalar = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
      auto ceed_attrs = mat_op.GetCeedBdrAttributes(attr);
      if (ceed_attrs.Size() > 0)
      {
        nn_fbi.AddMaterialProperty(ceed_attrs, inv_z0_scalar, -omega);
      }
    }
  }
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr), cond_i(max_bdr_attr);
    surf_sigma_op.AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
    if (!cond_i.empty())
    {
      cond_i *= -1.0;
      nn_fbi.AddCoefficient(cond_i.GetAttributeToMaterial(),
                            cond_i.GetMaterialProperties());
    }
  }

  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
  if (!nn_fbr.empty())
  {
    annr.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
  }
  auto Annr_assembled =
      ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  std::unique_ptr<mfem::HypreParMatrix> Anni_assembled;
  {
    const bool has_imag = mat_op.HasLossTangent() || !nn_fbi.empty();
    if (has_imag)
    {
      const int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient posepsi_h1_func(n_attr);
      if (mat_op.HasLossTangent())
      {
        if (normal)
        {
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImag(), omega * omega);
          posepsi_h1_func.NormalProjectedCoefficient(*normal);
        }
        else
        {
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImagScalar(), omega * omega);
        }
      }
      BilinearForm anni(h1_fespace);
      if (!posepsi_h1_func.empty())
      {
        anni.AddDomainIntegrator<MassIntegrator>(posepsi_h1_func);
      }
      if (!nn_fbi.empty())
      {
        anni.AddBoundaryIntegrator<MassIntegrator>(nn_fbi);
      }
      Anni_assembled =
          ParOperator(anni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
    }
  }
  return {std::move(Annr_assembled), std::move(Anni_assembled)};
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
