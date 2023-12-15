// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_COEFFICIENT_HPP
#define PALACE_FEM_COEFFICIENT_HPP

#include <complex>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "models/materialoperator.hpp"

namespace palace
{

//
// Derived coefficients which compute single values on internal boundaries where a possibly
// discontinuous function is given as an input grid function. These are all cheap to
// construct by design. All methods assume the provided grid function is ready for parallel
// comm on shared faces after a call to ExchangeFaceNbrData.
//

// Base class for coefficients which need to evaluate a GridFunction in a domain element
// attached to a boundary element, or both domain elements on either side for internal
// boundaries.
class BdrGridFunctionCoefficient
{
protected:
  // XX TODO: For thread-safety (multiple threads evaluating a coefficient simultaneously),
  //          the FET, FET.Elem1, and FET.Elem2 objects cannot be shared
  const mfem::ParMesh &mesh;
  const std::unordered_map<int, int> &local_to_shared;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2, TF;

  void GetBdrElementNeighborTransformations(mfem::ElementTransformation &T,
                                            const mfem::IntegrationPoint &ip,
                                            mfem::Vector *C1 = nullptr);

public:
  BdrGridFunctionCoefficient(const mfem::ParMesh &mesh,
                             const std::unordered_map<int, int> &local_to_shared)
    : mesh(mesh), local_to_shared(local_to_shared)
  {
  }

  // For a boundary element, return the element transformation objects for the neighboring
  // domain elements. FET.Elem2 may be nullptr if the boundary is a true one-sided boundary,
  // but if it is shared with another subdomain then it will be populated. Expects
  // ParMesh::ExchangeFaceNbrData has been called already.
  static void GetBdrElementNeighborTransformations(
      int i, const mfem::ParMesh &mesh, const std::unordered_map<int, int> &local_to_shared,
      mfem::FaceElementTransformations &FET, mfem::IsoparametricTransformation &T1,
      mfem::IsoparametricTransformation &T2, const mfem::IntegrationPoint *ip = nullptr);

  // Return normal vector to the boundary element at an integration point (it is assumed
  // that the element transformation has already been configured at the integration point of
  // interest).
  static void GetNormal(mfem::ElementTransformation &T, mfem::Vector &normal)
  {
    normal.SetSize(T.GetSpaceDim());
    mfem::CalcOrtho(T.Jacobian(), normal);
    normal /= normal.Norml2();
  }
};

// Computes surface current J_s = n x H on boundaries from B as a vector grid function
// where n is an inward normal (computes -n x H for outward normal n). For a two-sided
// internal boundary, the contributions from both sides add.
class BdrCurrentVectorCoefficient : public mfem::VectorCoefficient,
                                    public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &B;
  const MaterialOperator &mat_op;
  mfem::Vector C1, W, VU, VL, nor;

public:
  BdrCurrentVectorCoefficient(const mfem::ParGridFunction &gf,
                              const MaterialOperator &mat_op)
    : mfem::VectorCoefficient(mat_op.SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                 mat_op.GetLocalToSharedFaceMap()),
      B(gf), mat_op(mat_op), C1(gf.VectorDim()), W(gf.VectorDim()), VU(gf.VectorDim()),
      VL(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(vdim == 3, "BdrJVectorCoefficient expects a mesh in 3D space!");
    GetBdrElementNeighborTransformations(T, ip, &C1);

    // For interior faces, compute J_s = -n x H = -n x μ⁻¹(B1 - B2), where B1 (B2) is B in
    // el1 (el2) and n points out from el1.
    B.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), W);
    mat_op.GetInvPermeability(FET.Elem1->Attribute).Mult(W, VU);
    if (FET.Elem2)
    {
      // Double-sided, not a true boundary.
      B.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
      mat_op.GetInvPermeability(FET.Elem2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
    GetNormal(T, nor);
    V.SetSize(vdim);
    if (C1 * nor < 0.0)
    {
      V[0] = -nor[1] * VU[2] + nor[2] * VU[1];
      V[1] = -nor[2] * VU[0] + nor[0] * VU[2];
      V[2] = -nor[0] * VU[1] + nor[1] * VU[0];
    }
    else
    {
      V[0] = nor[1] * VU[2] - nor[2] * VU[1];
      V[1] = nor[2] * VU[0] - nor[0] * VU[2];
      V[2] = nor[0] * VU[1] - nor[1] * VU[0];
    }
  }
};

// Computes a single-valued surface charge ρ_s = D ⋅ n on boundaries from E given as a
// vector grid function. For a two-sided internal boundary, the contributions from both
// sides add.
class BdrChargeCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &E;
  const MaterialOperator &mat_op;
  mfem::Vector C1, W, VU, VL, nor;

public:
  BdrChargeCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &mat_op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      mat_op.GetLocalToSharedFaceMap()),
      E(gf), mat_op(mat_op), C1(gf.VectorDim()), W(gf.VectorDim()), VU(gf.VectorDim()),
      VL(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    GetBdrElementNeighborTransformations(T, ip, &C1);

    // For interior faces, compute D ⋅ n = ε (E1 - E2) ⋅ n, where E1 (E2) is E in el1 (el2)
    // to get a single-valued function.
    E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), W);
    mat_op.GetPermittivityReal(FET.Elem1->Attribute).Mult(W, VU);
    if (FET.Elem2)
    {
      E.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
      mat_op.GetPermittivityReal(FET.Elem2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
    GetNormal(T, nor);
    return (C1 * nor < 0.0) ? -(VU * nor) : VU * nor;
  }
};

// Computes the flux Φ_s = B ⋅ n on interior boundary elements using the user specified
// normal direction. Manually implements InnerProductCoefficient and
// VectorGridFunctionCoefficient to allow for evaluating the flux on internal boundaries.
class BdrFluxCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &B;
  const mfem::Vector dir;
  mfem::Vector V, VL, nor;

public:
  BdrFluxCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &mat_op,
                     const mfem::Vector &d)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      mat_op.GetLocalToSharedFaceMap()),
      B(gf), dir(d), V(gf.VectorDim()), VL(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    GetBdrElementNeighborTransformations(T, ip);

    // For interior faces, compute the average value. Since this is only used for
    // continuous (normal or tangential) values, we don't care that we average out the
    // discontinuous (tangential or normal) parts.
    B.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
    if (FET.Elem2)
    {
      B.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), VL);
      V += VL;
      V *= 0.5;
    }

    // Orient sign with the global direction.
    GetNormal(T, nor);
    return (dir * nor < 0.0) ? -(V * nor) : V * nor;
  }
};

// Helper for DielectricInterfaceCoefficient.
enum class DielectricInterfaceType
{
  DEFAULT,
  MA,
  MS,
  SA
};

// Computes a single-valued α Eᵀ E on boundaries from E given as a vector grid function.
// Uses the neighbor element on a user specified side to compute a single-sided value for
// potentially discontinuous solutions for an interior boundary element. The four cases
// correspond to a generic interface vs. specializations for metal-air, metal-substrate,
// and subtrate-air interfaces following:
//   J. Wenner et al., Surface loss simulations of superconducting coplanar waveguide
//     resonators, Appl. Phys. Lett. (2011).
template <DielectricInterfaceType Type>
class DielectricInterfaceCoefficient : public mfem::Coefficient,
                                       public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &E;
  const MaterialOperator &mat_op;
  const double ts, epsilon;
  const mfem::Vector side;
  mfem::Vector C1, V, nor;

  int Initialize(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                 mfem::Vector &V)
  {
    // Get neighboring elements.
    GetBdrElementNeighborTransformations(T, ip, &C1);

    // Get the single-sided solution.
    if (!FET.Elem2)
    {
      // Ignore side, solution is single-valued.
      E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
      return FET.Elem1->Attribute;
    }
    if (!side.Size())
    {
      // With no side specified, try to take the solution from the element which corresponds
      // to the vacuum domain, or at least the one with the higher speed of light.
      if (mat_op.GetLightSpeedMin(FET.Elem2->Attribute) >
          mat_op.GetLightSpeedMax(FET.Elem1->Attribute))
      {
        E.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), V);
        return FET.Elem2->Attribute;
      }
      E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
      return FET.Elem1->Attribute;
    }
    if (C1 * side < 0.0)
    {
      // Get solution in el2.
      E.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), V);
      return FET.Elem2->Attribute;
    }
    // Get solution in el1.
    E.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
    return FET.Elem1->Attribute;
  }

public:
  DielectricInterfaceCoefficient(const mfem::ParGridFunction &gf,
                                 const MaterialOperator &mat_op, double ti, double ei,
                                 const mfem::Vector &s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      mat_op.GetLocalToSharedFaceMap()),
      E(gf), mat_op(mat_op), ts(ti), epsilon(ei), side(s), C1(gf.VectorDim()),
      V(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    MFEM_ABORT("DielectricInterfaceCoefficient::Eval() is not implemented for this "
               "interface type!");
    return 0.0;
  }
};

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::MA>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  Initialize(T, ip, V);
  GetNormal(T, nor);

  // Metal-air interface: 0.5 * t / ε_MA * |E_n|² .
  double Vn = V * nor;
  return 0.5 * ts / epsilon * (Vn * Vn);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::MS>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  int attr = Initialize(T, ip, V);
  GetNormal(T, nor);

  // Metal-substrate interface: 0.5 * t * (ε_S)² / ε_MS * |E_n|² .
  const double Vn = V * nor;
  const double epsilon_S = mat_op.GetPermittivityReal(attr).InnerProduct(nor, nor);
  return 0.5 * ts * std::pow(epsilon_S, 2) / epsilon * (Vn * Vn);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::SA>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  Initialize(T, ip, V);
  GetNormal(T, nor);

  // Substrate-air interface: 0.5 * t * (ε_SA * |E_t|² + 1 / ε_MS * |E_n|²) .
  double Vn = V * nor;
  V.Add(-Vn, nor);
  return 0.5 * ts * (epsilon * (V * V) + (Vn * Vn) / epsilon);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::DEFAULT>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  Initialize(T, ip, V);

  // No specific interface, use full field evaluation: 0.5 * t * ε * |E|² .
  return 0.5 * ts * epsilon * (V * V);
}

// Helper for EnergyDensityCoefficient.
enum class EnergyDensityType
{
  ELECTRIC,
  MAGNETIC
};

// Returns the local energy density evaluated as 1/2 Dᴴ E or 1/2 Bᴴ H for real-valued
// material coefficients. For internal boundary elements, the solution is taken on the side
// of the element with the larger-valued material property (permittivity or permeability).
template <EnergyDensityType Type, typename GridFunctionType>
class EnergyDensityCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const GridFunctionType &U;
  const MaterialOperator &mat_op;
  mfem::Vector V;

  double GetLocalEnergyDensity(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip, int attr);

public:
  EnergyDensityCoefficient(const GridFunctionType &gf, const MaterialOperator &mat_op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      mat_op.GetLocalToSharedFaceMap()),
      U(gf), mat_op(mat_op), V(mat_op.SpaceDimension())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      return GetLocalEnergyDensity(T, ip, T.Attribute);
    }
    if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT)
    {
      // Get neighboring elements.
      GetBdrElementNeighborTransformations(T, ip);

      // For interior faces, compute the value on the side where the speed of light is
      // smaller (typically should choose the non-vacuum side).
      if (FET.Elem2 && mat_op.GetLightSpeedMax(FET.Elem2->Attribute) <
                           mat_op.GetLightSpeedMin(FET.Elem1->Attribute))
      {
        return GetLocalEnergyDensity(*FET.Elem2, FET.Elem2->GetIntPoint(),
                                     FET.Elem2->Attribute);
      }
      else
      {
        return GetLocalEnergyDensity(*FET.Elem1, FET.Elem1->GetIntPoint(),
                                     FET.Elem1->Attribute);
      }
    }
    MFEM_ABORT("Unsupported element type in EnergyDensityCoefficient!");
    return 0.0;
  }
};

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, mfem::ParComplexGridFunction>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  // Only the real part of the permittivity contributes to the energy (imaginary part
  // cancels out in the inner product due to symmetry).
  U.real().GetVectorValue(T, ip, V);
  double res = mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
  U.imag().GetVectorValue(T, ip, V);
  res += mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
  return 0.5 * res;
}

template <>
inline double EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, mfem::ParGridFunction>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  U.GetVectorValue(T, ip, V);
  return 0.5 * mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
}

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, mfem::ParComplexGridFunction>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  U.real().GetVectorValue(T, ip, V);
  double res = mat_op.GetInvPermeability(attr).InnerProduct(V, V);
  U.imag().GetVectorValue(T, ip, V);
  res += mat_op.GetInvPermeability(attr).InnerProduct(V, V);
  return 0.5 * res;
}

template <>
inline double EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, mfem::ParGridFunction>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  U.GetVectorValue(T, ip, V);
  return 0.5 * mat_op.GetInvPermeability(attr).InnerProduct(V, V);
}

// Returns the local field evaluated on a boundary element. For internal boundary elements,
// the solution is taken on the side of the element with the larger-valued material
// property (permittivity or permeability).
class BdrFieldVectorCoefficient : public mfem::VectorCoefficient,
                                  public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &U;
  const MaterialOperator &mat_op;

public:
  BdrFieldVectorCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &mat_op)
    : mfem::VectorCoefficient(mat_op.SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                 mat_op.GetLocalToSharedFaceMap()),
      U(gf), mat_op(mat_op)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    GetBdrElementNeighborTransformations(T, ip);

    // For interior faces, compute the value on the side where the speed of light is
    // smaller (typically should choose the non-vacuum side).
    if (FET.Elem2 && mat_op.GetLightSpeedMax(FET.Elem2->Attribute) <
                         mat_op.GetLightSpeedMin(FET.Elem1->Attribute))
    {
      U.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), V);
    }
    else
    {
      U.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
    }
  }
};

class BdrFieldCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &U;
  const MaterialOperator &mat_op;

public:
  BdrFieldCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &mat_op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      mat_op.GetLocalToSharedFaceMap()),
      U(gf), mat_op(mat_op)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    GetBdrElementNeighborTransformations(T, ip);

    // For interior faces, compute the value on the side where the speed of light is
    // smaller (typically should choose the non-vacuum side).
    if (FET.Elem2 && mat_op.GetLightSpeedMax(FET.Elem2->Attribute) <
                         mat_op.GetLightSpeedMin(FET.Elem1->Attribute))
    {
      return U.GetValue(*FET.Elem2, FET.Elem2->GetIntPoint());
    }
    else
    {
      return U.GetValue(*FET.Elem1, FET.Elem1->GetIntPoint());
    }
  }
};

//
// More helpful coefficient types. Wrapper coefficients allow additions of scalar and vector
// or matrix coefficients. Restricted coefficients only compute the coefficient if for the
// given list of attributes. Sum coefficients own a list of coefficients to add.
//

class VectorWrappedCoefficient : public mfem::VectorCoefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coeff;

public:
  VectorWrappedCoefficient(int dim, std::unique_ptr<mfem::Coefficient> &&coeff)
    : mfem::VectorCoefficient(dim), coeff(std::move(coeff))
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    V.SetSize(vdim);
    V = coeff->Eval(T, ip);
  }
};

class MatrixWrappedCoefficient : public mfem::MatrixCoefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coeff;

public:
  MatrixWrappedCoefficient(int dim, std::unique_ptr<mfem::Coefficient> &&coeff)
    : mfem::MatrixCoefficient(dim), coeff(std::move(coeff))
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    K.Diag(coeff->Eval(T, ip), height);
  }
};

class RestrictedCoefficient : public mfem::Coefficient
{
private:
  std::unique_ptr<mfem::Coefficient> coeff;
  const mfem::Array<int> &attr;

public:
  RestrictedCoefficient(std::unique_ptr<mfem::Coefficient> &&coeff,
                        const mfem::Array<int> &attr)
    : mfem::Coefficient(), coeff(std::move(coeff)), attr(attr)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    return (attr.Find(T.Attribute) < 0) ? 0.0 : coeff->Eval(T, ip);
  }
};

class RestrictedVectorCoefficient : public mfem::VectorCoefficient
{
private:
  std::unique_ptr<mfem::VectorCoefficient> coeff;
  const mfem::Array<int> &attr;

public:
  RestrictedVectorCoefficient(std::unique_ptr<mfem::VectorCoefficient> &&coeff,
                              const mfem::Array<int> &attr)
    : mfem::VectorCoefficient(coeff->GetVDim()), coeff(std::move(coeff)), attr(attr)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    if (attr.Find(T.Attribute) < 0)
    {
      V.SetSize(vdim);
      V = 0.0;
    }
    else
    {
      coeff->Eval(V, T, ip);
    }
  }
};

class RestrictedMatrixCoefficient : public mfem::MatrixCoefficient
{
private:
  std::unique_ptr<mfem::MatrixCoefficient> coeff;
  const mfem::Array<int> &attr;

public:
  RestrictedMatrixCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coeff,
                              const mfem::Array<int> &attr)
    : mfem::MatrixCoefficient(coeff->GetHeight(), coeff->GetWidth()),
      coeff(std::move(coeff)), attr(attr)
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    if (attr.Find(T.Attribute) < 0)
    {
      K.SetSize(height, width);
      K = 0.0;
    }
    else
    {
      coeff->Eval(K, T, ip);
    }
  }
};

class SumCoefficient : public mfem::Coefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::Coefficient>, double>> c;

public:
  SumCoefficient() : mfem::Coefficient() {}

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coeff, double a = 1.0)
  {
    c.emplace_back(std::move(coeff), a);
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    double val = 0.0;
    for (auto &[coeff, a] : c)
    {
      val += a * coeff->Eval(T, ip);
    }
    return val;
  }
};

class SumVectorCoefficient : public mfem::VectorCoefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::VectorCoefficient>, double>> c;

public:
  SumVectorCoefficient(int d) : mfem::VectorCoefficient(d) {}

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::VectorCoefficient> &&coeff, double a = 1.0)
  {
    MFEM_VERIFY(coeff->GetVDim() == vdim,
                "Invalid VectorCoefficient dimensions for SumVectorCoefficient!");
    c.emplace_back(std::move(coeff), a);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coeff, double a = 1.0)
  {
    c.emplace_back(std::make_unique<VectorWrappedCoefficient>(vdim, std::move(coeff)), a);
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector U(vdim);
    V.SetSize(vdim);
    V = 0.0;
    for (auto &[coeff, a] : c)
    {
      coeff->Eval(U, T, ip);
      V.Add(a, U);
    }
  }
};

class SumMatrixCoefficient : public mfem::MatrixCoefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::MatrixCoefficient>, double>> c;

public:
  SumMatrixCoefficient(int d) : mfem::MatrixCoefficient(d) {}
  SumMatrixCoefficient(int h, int w) : mfem::MatrixCoefficient(h, w) {}

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coeff, double a)
  {
    MFEM_VERIFY(coeff->GetHeight() == height && coeff->GetWidth() == width,
                "Invalid MatrixCoefficient dimensions for SumMatrixCoefficient!");
    c.emplace_back(std::move(coeff), a);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coeff, double a)
  {
    MFEM_VERIFY(width == height, "MatrixWrappedCoefficient can only be constructed for "
                                 "square MatrixCoefficient objects!");
    c.emplace_back(std::make_unique<MatrixWrappedCoefficient>(height, std::move(coeff)), a);
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::DenseMatrix M(height, width);
    K.SetSize(height, width);
    K = 0.0;
    for (auto &[coeff, a] : c)
    {
      coeff->Eval(M, T, ip);
      K.Add(a, M);
    }
  }
};

}  // namespace palace

#endif  // PALACE_FEM_COEFFICIENT_HPP
