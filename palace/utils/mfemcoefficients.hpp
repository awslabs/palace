// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MFEM_COEFFICIENTS_HPP
#define PALACE_MFEM_COEFFICIENTS_HPP

#include <complex>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "fem/materialoperator.hpp"

namespace palace
{

//
// Derived coefficients which compute single values on internal boundaries where a possibly
// discontinuous function is given as an input grid function. These are all cheap to
// construct by design. All methods assume the provided grid function is ready for parallel
// comm on shared faces after a call to ExchangeFaceNbrData.
//

class BdrGridFunctionCoefficient
{
private:
  // From mfem::GridFunction::GetVectorValue.
  static mfem::IntegrationPoint be_to_bfe(mfem::Geometry::Type geom, int o,
                                          const mfem::IntegrationPoint &ip);

protected:
  mfem::ParMesh &mesh;
  const std::map<int, int> &local_to_shared;

  void GetElementTransformations(mfem::ElementTransformation &T,
                                 const mfem::IntegrationPoint &ip,
                                 mfem::ElementTransformation *&T1,
                                 mfem::ElementTransformation *&T2,
                                 mfem::Vector *C1 = nullptr);

public:
  BdrGridFunctionCoefficient(mfem::ParMesh &msh, const std::map<int, int> &l2s)
    : mesh(msh), local_to_shared(l2s)
  {
  }

  static void GetNormal(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                        mfem::Vector &normal);
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

public:
  BdrCurrentVectorCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op,
                              const std::map<int, int> &l2s)
    : mfem::VectorCoefficient(gf.ParFESpace()->GetParMesh()->SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s), B(gf), mat_op(op)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_VERIFY(vdim == 3, "BdrJVectorCoefficient expects a mesh in 3D space!");
    mfem::ElementTransformation *T1, *T2;
    mfem::Vector C1;
    GetElementTransformations(T, ip, T1, T2, &C1);

    // For interior faces, compute J_s = -n x H = -n x μ⁻¹(B1 - B2), where B1 (B2) is B in
    // el1 (el2) and n points out from el1.
    mfem::Vector W, VU(vdim);
    B.GetVectorValue(*T1, T1->GetIntPoint(), W);
    mat_op.GetInvPermeability(T1->Attribute).Mult(W, VU);
    if (T2)
    {
      // Double-sided, not a true boundary.
      mfem::Vector VL(vdim);
      B.GetVectorValue(*T2, T2->GetIntPoint(), W);
      mat_op.GetInvPermeability(T2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
    mfem::Vector nor;
    GetNormal(T, ip, nor);
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

public:
  BdrChargeCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op,
                       const std::map<int, int> &l2s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s),
      E(gf), mat_op(op)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    mfem::Vector C1;
    GetElementTransformations(T, ip, T1, T2, &C1);

    // For interior faces, compute D ⋅ n = ε (E1 - E2) ⋅ n, where E1 (E2) is E in el1 (el2)
    // to get a single-valued function.
    mfem::Vector W, VU(E.VectorDim());
    E.GetVectorValue(*T1, T1->GetIntPoint(), W);
    mat_op.GetPermittivityReal(T1->Attribute).Mult(W, VU);
    if (T2)
    {
      mfem::Vector VL(E.VectorDim());
      E.GetVectorValue(*T2, T2->GetIntPoint(), W);
      mat_op.GetPermittivityReal(T2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
    mfem::Vector nor;
    GetNormal(T, ip, nor);
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

public:
  BdrFluxCoefficient(const mfem::ParGridFunction &gf, mfem::Vector d,
                     const std::map<int, int> &l2s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s),
      B(gf), dir(std::move(d))
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    GetElementTransformations(T, ip, T1, T2);

    // For interior faces, compute the average value. Since this is only used for
    // continuous (normal or tangential) values, we don't care that we average out the
    // discontinuous (tangential or normal) parts.
    mfem::Vector V;
    B.GetVectorValue(*T1, T1->GetIntPoint(), V);
    if (T2)
    {
      mfem::Vector VL;
      B.GetVectorValue(*T2, T2->GetIntPoint(), VL);
      V += VL;
      V *= 0.5;
    }

    // Orient sign with the global direction.
    mfem::Vector nor;
    GetNormal(T, ip, nor);
    return (dir * nor < 0.0) ? -(V * nor) : V * nor;
  }
};

// Computes a single-valued α Eᵀ E on boundaries from E given as a vector grid function.
// Uses the neighbor element on a user specified side to compute a single-sided value for
// potentially discontinuous solutions for an interior boundary element. The four cases
// correspond to a generic interface vs. specializations for metal-air, metal-substrate,
// and subtrate-air interfaces following:
//   J. Wenner et al., Surface loss simulations of superconducting coplanar waveguide
//     resonators, Appl. Phys. Lett. (2011).
enum class DielectricInterfaceType
{
  DEFAULT,
  MA,
  MS,
  SA,
  INVALID = -1
};

template <DielectricInterfaceType Type>
class DielectricInterfaceCoefficient : public mfem::Coefficient,
                                       public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &E;
  const MaterialOperator &mat_op;
  const double ts, epsilon;
  const mfem::Vector side;

  int Initialize(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                 mfem::Vector &V)
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    mfem::Vector C1;
    GetElementTransformations(T, ip, T1, T2, &C1);

    // Get the single-sided solution.
    if (!T2)
    {
      // Ignore side, solution is single-valued.
      E.GetVectorValue(*T1, T1->GetIntPoint(), V);
      return T1->Attribute;
    }
    if (!side.Size())
    {
      // With no side specified, try to take the solution from the element which corresponds
      // to the vacuum domain, or at least the one with the higher speed of light.
      if (mat_op.GetLightSpeedMin(T2->Attribute) > mat_op.GetLightSpeedMax(T1->Attribute))
      {
        E.GetVectorValue(*T2, T2->GetIntPoint(), V);
        return T2->Attribute;
      }
      E.GetVectorValue(*T1, T1->GetIntPoint(), V);
      return T1->Attribute;
    }
    if (C1 * side < 0.0)
    {
      // Get solution in el2.
      E.GetVectorValue(*T2, T2->GetIntPoint(), V);
      return T2->Attribute;
    }
    // Get solution in el1.
    E.GetVectorValue(*T1, T1->GetIntPoint(), V);
    return T1->Attribute;
  }

public:
  DielectricInterfaceCoefficient(const mfem::ParGridFunction &gf,
                                 const MaterialOperator &op, double ti, double ei,
                                 mfem::Vector s, const std::map<int, int> &l2s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s),
      E(gf), mat_op(op), ts(ti), epsilon(ei), side(std::move(s))
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
  mfem::Vector V, nor;
  Initialize(T, ip, V);
  GetNormal(T, ip, nor);

  // Metal-air interface: 0.5 * t / ϵ_MA * |E_n|² .
  double Vn = V * nor;
  return 0.5 * ts / epsilon * (Vn * Vn);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::MS>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  mfem::Vector V, nor;
  int attr = Initialize(T, ip, V);
  GetNormal(T, ip, nor);

  // Metal-substrate interface: 0.5 * t * (ϵ_S)² / ϵ_MS * |E_n|² .
  const double Vn = V * nor;
  const double epsilon_S = mat_op.GetPermittivityReal(attr).InnerProduct(nor, nor);
  return 0.5 * ts * std::pow(epsilon_S, 2) / epsilon * (Vn * Vn);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::SA>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  mfem::Vector V, nor;
  Initialize(T, ip, V);
  GetNormal(T, ip, nor);

  // Substrate-air interface: 0.5 * t * (ϵ_SA * |E_t|² + 1 / ϵ_MS * |E_n|²) .
  double Vn = V * nor;
  mfem::Vector Vt(V.Size());
  add(V, -Vn, nor, Vt);
  return 0.5 * ts * (epsilon * (Vt * Vt) + (Vn * Vn) / epsilon);
}

template <>
inline double DielectricInterfaceCoefficient<DielectricInterfaceType::DEFAULT>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution and neighboring element attribute.
  mfem::Vector V;
  Initialize(T, ip, V);

  // No specific interface, use full field evaluation: 0.5 * t * ϵ * |E|² .
  return 0.5 * ts * epsilon * (V * V);
}

// Returns the property value of the material for the given index. When called on a
// boundary element (external only), uses the material property of the neighboring domain
// element.
enum class MaterialPropertyType
{
  INV_PERMEABILITY,
  PERMITTIVITY_REAL,
  PERMITTIVITY_IMAG,
  PERMITTIVITY_ABS,
  CONDUCTIVITY,
  INV_LONDON_DEPTH,
  INV_Z0,
  INV_PERMEABILITY_C0,
  INVALID = -1
};

template <MaterialPropertyType Type>
class MaterialPropertyCoefficient : public mfem::MatrixCoefficient
{
private:
  const MaterialOperator &mat_op;
  const double coef;

  static int GetAttribute(mfem::ElementTransformation &T)
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      return T.mesh->GetAttribute(T.ElementNo);
    }
    if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT)
    {
      int i, o;
      int iel1, iel2, info1, info2;
      const mfem::Mesh &mesh = *T.mesh;
      mesh.GetBdrElementFace(T.ElementNo, &i, &o);
      mesh.GetFaceElements(i, &iel1, &iel2);
      mesh.GetFaceInfos(i, &info1, &info2);
      MFEM_VERIFY(info2 < 0, "MaterialPropertyCoefficient should only be used for exterior "
                             "(single-sided) boundaries!");
      return mesh.GetAttribute(iel1);
    }
    MFEM_ABORT("Unsupported element type in MaterialPropertyCoefficient!");
    return 0;
  }

public:
  MaterialPropertyCoefficient(const MaterialOperator &op, double c = 1.0)
    : mfem::MatrixCoefficient(op.SpaceDimension()), mat_op(op), coef(c)
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    MFEM_ABORT("MaterialPropertyCoefficient::Eval() is not implemented for this "
               "material property type!");
  }
};

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetInvPermeability(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetPermittivityReal(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_IMAG>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetPermittivityImag(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_ABS>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetPermittivityAbs(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::CONDUCTIVITY>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetConductivity(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::INV_LONDON_DEPTH>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetInvLondonDepth(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::INV_Z0>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  K = mat_op.GetInvImpedance(GetAttribute(T));
  K *= coef;
}

template <>
inline void MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY_C0>::Eval(
    mfem::DenseMatrix &K, mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  const int attr = GetAttribute(T);
  Mult(mat_op.GetInvPermeability(attr), mat_op.GetLightSpeed(attr), K);
  K *= coef;
}

// Returns the local energy density evaluated as 1/2 Dᴴ E or 1/2 Bᴴ H for real-valued
// material coefficients. For internal boundary elements, the solution is taken on the side
// of the element with the larger-valued material property (permittivity or permeability).
enum class EnergyDensityType
{
  ELECTRIC,
  ELECTRIC_LOSS,
  MAGNETIC,
  INVALID = -1
};
enum class EnergyDensityValueType
{
  COMPLEX,
  REAL,
  INVALID = -1
};

template <EnergyDensityType Type, EnergyDensityValueType ValueType>
class EnergyDensityCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  typedef
      typename std::conditional<ValueType == EnergyDensityValueType::COMPLEX,
                                mfem::ParComplexGridFunction, mfem::ParGridFunction>::type
          GridFunctionType;
  const GridFunctionType &U;
  const MaterialOperator &mat_op;

  double GetLocalEnergyDensity(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip, int attr)
  {
    MFEM_ABORT(
        "EnergyDensityCoefficient::GetLocalEnergyDensity() is not implemented for this "
        "value type!");
    return 0.0;
  }

public:
  EnergyDensityCoefficient(const GridFunctionType &gf, const MaterialOperator &op,
                           const std::map<int, int> &l2s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s),
      U(gf), mat_op(op)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      T.SetIntPoint(&ip);
      return GetLocalEnergyDensity(T, ip, mesh.GetAttribute(T.ElementNo));
    }
    if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT)
    {
      // Get neighboring elements.
      mfem::ElementTransformation *T1, *T2;
      GetElementTransformations(T, ip, T1, T2);

      // For interior faces, compute the value on the side where the material property is
      // larger (typically should choose the non-vacuum side).
      if (T2 &&
          mat_op.GetLightSpeedMax(T2->Attribute) < mat_op.GetLightSpeedMin(T1->Attribute))
      {
        return GetLocalEnergyDensity(*T2, T2->GetIntPoint(), T2->Attribute);
      }
      else
      {
        return GetLocalEnergyDensity(*T1, T1->GetIntPoint(), T1->Attribute);
      }
    }
    MFEM_ABORT("Unsupported element type in EnergyDensityCoefficient!");
    return 0.0;
  }
};

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::COMPLEX>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  // Only the real part of the permittivity contributes to the energy (imaginary part
  // cancels out in the inner product due to symmetry).
  mfem::Vector V;
  U.real().GetVectorValue(T, ip, V);
  double res = mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
  U.imag().GetVectorValue(T, ip, V);
  res += mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
  return 0.5 * res;
}

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::REAL>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  mfem::Vector V;
  U.GetVectorValue(T, ip, V);
  return 0.5 * mat_op.GetPermittivityReal(attr).InnerProduct(V, V);
}

template <>
inline double EnergyDensityCoefficient<EnergyDensityType::ELECTRIC_LOSS,
                                       EnergyDensityValueType::COMPLEX>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  mfem::Vector V;
  U.real().GetVectorValue(T, ip, V);
  double res = mat_op.GetPermittivityImag(attr).InnerProduct(V, V);
  U.imag().GetVectorValue(T, ip, V);
  res += mat_op.GetPermittivityImag(attr).InnerProduct(V, V);
  return -0.5 * res;
}

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::ELECTRIC_LOSS, EnergyDensityValueType::REAL>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  mfem::Vector V;
  U.GetVectorValue(T, ip, V);
  return -0.5 * mat_op.GetPermittivityImag(attr).InnerProduct(V, V);
}

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, EnergyDensityValueType::COMPLEX>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  mfem::Vector V;
  U.real().GetVectorValue(T, ip, V);
  double res = mat_op.GetInvPermeability(attr).InnerProduct(V, V);
  U.imag().GetVectorValue(T, ip, V);
  res += mat_op.GetInvPermeability(attr).InnerProduct(V, V);
  return 0.5 * res;
}

template <>
inline double
EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, EnergyDensityValueType::REAL>::
    GetLocalEnergyDensity(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                          int attr)
{
  mfem::Vector V;
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
  BdrFieldVectorCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op,
                            const std::map<int, int> &l2s)
    : mfem::VectorCoefficient(gf.ParFESpace()->GetParMesh()->SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s), U(gf), mat_op(op)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    GetElementTransformations(T, ip, T1, T2);

    // For interior faces, compute the value on the side where the material property is
    // larger (typically should choose the non-vacuum side).
    if (T2 &&
        mat_op.GetLightSpeedMax(T2->Attribute) < mat_op.GetLightSpeedMin(T1->Attribute))
    {
      U.GetVectorValue(*T2, T2->GetIntPoint(), V);
    }
    else
    {
      U.GetVectorValue(*T1, T1->GetIntPoint(), V);
    }
  }
};

class BdrFieldCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &U;
  const MaterialOperator &mat_op;

public:
  BdrFieldCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op,
                      const std::map<int, int> &l2s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), l2s),
      U(gf), mat_op(op)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    GetElementTransformations(T, ip, T1, T2);

    // For interior faces, compute the value on the side where the material property is
    // larger (typically should choose the non-vacuum side).
    if (T2 &&
        mat_op.GetLightSpeedMax(T2->Attribute) < mat_op.GetLightSpeedMin(T1->Attribute))
    {
      return U.GetValue(*T2, T2->GetIntPoint());
    }
    else
    {
      return U.GetValue(*T1, T1->GetIntPoint());
    }
  }
};

// From mfem::GridFunction::GetVectorValue.
inline mfem::IntegrationPoint
BdrGridFunctionCoefficient::be_to_bfe(mfem::Geometry::Type geom, int o,
                                      const mfem::IntegrationPoint &ip)
{
  mfem::IntegrationPoint fip = {};
  if (geom == mfem::Geometry::TRIANGLE)
  {
    if (o == 2)
    {
      fip.x = 1.0 - ip.x - ip.y;
      fip.y = ip.x;
    }
    else if (o == 4)
    {
      fip.x = ip.y;
      fip.y = 1.0 - ip.x - ip.y;
    }
    else
    {
      fip.x = ip.x;
      fip.y = ip.y;
    }
    fip.z = ip.z;
  }
  else
  {
    if (o == 2)
    {
      fip.x = ip.y;
      fip.y = 1.0 - ip.x;
    }
    else if (o == 4)
    {
      fip.x = 1.0 - ip.x;
      fip.y = 1.0 - ip.y;
    }
    else if (o == 6)
    {
      fip.x = 1.0 - ip.y;
      fip.y = ip.x;
    }
    else
    {
      fip.x = ip.x;
      fip.y = ip.y;
    }
    fip.z = ip.z;
  }
  fip.weight = ip.weight;
  fip.index = ip.index;
  return fip;
}

inline void BdrGridFunctionCoefficient::GetElementTransformations(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
    mfem::ElementTransformation *&T1, mfem::ElementTransformation *&T2, mfem::Vector *C1)
{
  // Return transformations for elements attached to boundary element T. T1 always exists
  // but T2 may not if the element is truly a single-sided boundary.
  MFEM_VERIFY(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
              "Unexpected element type in BdrGridFunctionCoefficient!");
  MFEM_VERIFY(&mesh == T.mesh, "Invalid mesh for BdrGridFunctionCoefficient!");
  int i, o;
  int iel1, iel2, info1, info2;
  mesh.GetBdrElementFace(T.ElementNo, &i, &o);
  mesh.GetFaceElements(i, &iel1, &iel2);
  mesh.GetFaceInfos(i, &info1, &info2);

  mfem::FaceElementTransformations *FET;
  if (info2 >= 0 && iel2 < 0)
  {
    // Face is shared with another subdomain.
    const int &ishared = local_to_shared.at(i);
    FET = mesh.GetSharedFaceTransformations(ishared);
  }
  else
  {
    // Face is either internal to the subdomain, or a true one-sided boundary.
    FET = mesh.GetFaceElementTransformations(i);
  }

  // Boundary elements and boundary faces may have different orientations so adjust the
  // integration point if necessary. See mfem::GridFunction::GetValue and GetVectorValue.
  mfem::IntegrationPoint fip = be_to_bfe(FET->GetGeometryType(), o, ip);
  FET->SetAllIntPoints(&fip);
  T1 = &FET->GetElement1Transformation();
  T2 = (info2 >= 0) ? &FET->GetElement2Transformation() : nullptr;

  // If desired, get vector pointing from center of boundary element into element 1 for
  // orientations.
  if (C1)
  {
    mfem::Vector CF(T.GetSpaceDim());
    mfem::ElementTransformation &TF = *mesh.GetFaceTransformation(i);
    TF.Transform(mfem::Geometries.GetCenter(mesh.GetFaceGeometry(i)), CF);

    C1->SetSize(T.GetSpaceDim());
    T1->Transform(mfem::Geometries.GetCenter(T1->GetGeometryType()), *C1);
    *C1 -= CF;  // Points into element 1 from the face
  }
}

inline void BdrGridFunctionCoefficient::GetNormal(mfem::ElementTransformation &T,
                                                  const mfem::IntegrationPoint &ip,
                                                  mfem::Vector &normal)
{
  // Return normal vector to the boundary element at the provided integration point.
  normal.SetSize(T.GetSpaceDim());
  T.SetIntPoint(&ip);
  mfem::CalcOrtho(T.Jacobian(), normal);
  normal /= normal.Norml2();
}

// Wraps a mfem::MatrixCoefficient to compute a scalar coefficient as nᵀ M n. Only works
// for square matrix coefficients of size equal to the spatial dimension.
class NormalProjectedCoefficient : public mfem::Coefficient
{
  std::unique_ptr<mfem::MatrixCoefficient> c;

public:
  NormalProjectedCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coef)
    : mfem::Coefficient(), c(std::move(coef))
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector nor;
    mfem::DenseMatrix K(c->GetHeight(), c->GetWidth());
    c->Eval(K, T, ip);
    BdrGridFunctionCoefficient::GetNormal(T, ip, nor);
    return K.InnerProduct(nor, nor);
  }
};

class VectorWrappedCoefficient : public mfem::VectorCoefficient
{
private:
  std::unique_ptr<mfem::Coefficient> c;

public:
  VectorWrappedCoefficient(int d, std::unique_ptr<mfem::Coefficient> &&coef)
    : mfem::VectorCoefficient(d), c(std::move(coef))
  {
  }

  void SetTime(double t) override
  {
    mfem::VectorCoefficient::SetTime(t);
    c->SetTime(t);
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector U(vdim);
    V.SetSize(vdim);
    V = c->Eval(T, ip);
  }
};

class MatrixWrappedCoefficient : public mfem::MatrixCoefficient
{
private:
  std::unique_ptr<mfem::Coefficient> c;

public:
  MatrixWrappedCoefficient(int d, std::unique_ptr<mfem::Coefficient> &&coef)
    : mfem::MatrixCoefficient(d), c(std::move(coef))
  {
  }

  void SetTime(double t) override
  {
    mfem::MatrixCoefficient::SetTime(t);
    c->SetTime(t);
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    K.Diag(c->Eval(T, ip), height);
  }
};

class SumCoefficient : public mfem::Coefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::Coefficient>, const mfem::Array<int> *>> c;

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> *marker)
  {
    c.emplace_back(std::move(coef), marker);
  }

public:

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef)
  {
    AddCoefficient(std::move(coef), nullptr);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> &marker)
  {
    AddCoefficient(std::move(coef), &marker);
  }

  void SetTime(double t) override
  {
    mfem::Coefficient::SetTime(t);
    for (auto &[coef, marker] : c)
    {
      coef->SetTime(t);
    }
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    double val = 0.0;
    for (auto &[coef, marker] : c)
    {
      if (!marker || (*marker)[T.Attribute - 1])
      {
        val += coef->Eval(T, ip);
      }
    }
    return val;
  }
};

class SumVectorCoefficient : public mfem::VectorCoefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::VectorCoefficient>, const mfem::Array<int> *>>
      c;

  void AddCoefficient(std::unique_ptr<mfem::VectorCoefficient> &&coef,
                      const mfem::Array<int> *marker)
  {
    MFEM_VERIFY(coef->GetVDim() == vdim,
                "Invalid VectorCoefficient dimensions for SumVectorCoefficient!");
    c.emplace_back(std::move(coef), marker);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> *marker)
  {
    c.emplace_back(std::make_unique<VectorWrappedCoefficient>(vdim, std::move(coef)),
                   marker);
  }

public:
  SumVectorCoefficient(int d) : mfem::VectorCoefficient(d) {}

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::VectorCoefficient> &&coef)
  {
    AddCoefficient(std::move(coef), nullptr);
  }

  void AddCoefficient(std::unique_ptr<mfem::VectorCoefficient> &&coef,
                      const mfem::Array<int> &marker)
  {
    AddCoefficient(std::move(coef), &marker);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef)
  {
    AddCoefficient(std::move(coef), nullptr);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> &marker)
  {
    AddCoefficient(std::move(coef), &marker);
  }

  void SetTime(double t) override
  {
    mfem::VectorCoefficient::SetTime(t);
    for (auto &[coef, marker] : c)
    {
      coef->SetTime(t);
    }
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::Vector U(vdim);
    V.SetSize(vdim);
    V = 0.0;
    for (auto &[coef, marker] : c)
    {
      if (!marker || (*marker)[T.Attribute - 1])
      {
        coef->Eval(U, T, ip);
        V += U;
      }
    }
  }
};

class SumMatrixCoefficient : public mfem::MatrixCoefficient
{
private:
  std::vector<std::pair<std::unique_ptr<mfem::MatrixCoefficient>, const mfem::Array<int> *>>
      c;

  void AddCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coef,
                      const mfem::Array<int> *marker)
  {
    MFEM_VERIFY(coef->GetHeight() == height && coef->GetWidth() == width,
                "Invalid MatrixCoefficient dimensions for SumMatrixCoefficient!");
    c.emplace_back(std::move(coef), marker);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> *marker)
  {
    MFEM_VERIFY(width == height, "MatrixWrappedCoefficient can only be constructed for "
                                 "square MatrixCoefficient objects!");
    c.emplace_back(std::make_unique<MatrixWrappedCoefficient>(height, std::move(coef)),
                   marker);
  }

public:
  SumMatrixCoefficient(int d) : mfem::MatrixCoefficient(d) {}
  SumMatrixCoefficient(int h, int w) : mfem::MatrixCoefficient(h, w) {}

  bool empty() const { return c.empty(); }

  void AddCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coef)
  {
    AddCoefficient(std::move(coef), nullptr);
  }

  void AddCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coef,
                      const mfem::Array<int> &marker)
  {
    AddCoefficient(std::move(coef), &marker);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef)
  {
    AddCoefficient(std::move(coef), nullptr);
  }

  void AddCoefficient(std::unique_ptr<mfem::Coefficient> &&coef,
                      const mfem::Array<int> &marker)
  {
    AddCoefficient(std::move(coef), &marker);
  }

  void SetTime(double t) override
  {
    mfem::MatrixCoefficient::SetTime(t);
    for (auto &[coef, marker] : c)
    {
      coef->SetTime(t);
    }
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::DenseMatrix M(height, width);
    K.SetSize(height, width);
    K = 0.0;
    for (auto &[coef, marker] : c)
    {
      if (!marker || (*marker)[T.Attribute - 1])
      {
        coef->Eval(M, T, ip);
        K += M;
      }
    }
  }
};

}  // namespace palace

#endif  // PALACE_MFEM_COEFFICIENTS_HPP
