// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_COEFFICIENT_HPP
#define PALACE_FEM_COEFFICIENT_HPP

#include <complex>
#include <map>
#include <memory>
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

enum class MaterialPropertyType
{
  INV_PERMEABILITY,
  PERMITTIVITY_REAL,
  PERMITTIVITY_IMAG,
  PERMITTIVITY_ABS,
  CONDUCTIVITY,
  INV_LONDON_DEPTH,
  INV_Z0,
  INV_PERMEABILITY_C0
};

enum class MeshElementType
{
  ELEMENT,
  BDR_ELEMENT,
  SUBMESH,
  BDR_SUBMESH
};

// Returns the property value of the material for the given index. Two separate classes for
// domain element access and boundary element access, which returns the material property of
// the neighboring domain element.
template <MaterialPropertyType MatType, MeshElementType ElemType = MeshElementType::ELEMENT>
class MaterialPropertyCoefficient : public mfem::MatrixCoefficient
{
private:
  const MaterialOperator &mat_op;
  const double coef;

  static int GetAttribute(mfem::ElementTransformation &T)
  {
    if constexpr (ElemType == MeshElementType::SUBMESH ||
                  ElemType == MeshElementType::BDR_SUBMESH)
    {
      MFEM_ASSERT(
          T.ElementType == mfem::ElementTransformation::ELEMENT,
          "Invalid usage of MaterialPropertyCoefficient for given MeshElementType!");
      const mfem::ParSubMesh &submesh = *static_cast<mfem::ParSubMesh *>(T.mesh);
      const mfem::ParMesh &mesh = *submesh.GetParent();
      if constexpr (ElemType == MeshElementType::SUBMESH)
      {
        MFEM_ASSERT(
            const_cast<mfem::ParSubMesh &>(submesh).GetFrom() ==
                mfem::SubMesh::From::Domain,
            "Invalid usage of MaterialPropertyCoefficient for given MeshElementType!");
        return mesh.GetAttribute(submesh.GetParentElementIDMap()[T.ElementNo]);
      }
      else if constexpr (ElemType == MeshElementType::BDR_SUBMESH)
      {
        MFEM_ASSERT(
            const_cast<mfem::ParSubMesh &>(submesh).GetFrom() ==
                mfem::SubMesh::From::Boundary,
            "Invalid usage of MaterialPropertyCoefficient for given MeshElementType!");
        int i, o, iel1, iel2;
        mesh.GetBdrElementFace(submesh.GetParentElementIDMap()[T.ElementNo], &i, &o);
        mesh.GetFaceElements(i, &iel1, &iel2);
#ifdef MFEM_DEBUG
        int info1, info2, nc;
        mesh.GetFaceInfos(i, &info1, &info2, &nc);
        MFEM_VERIFY(nc == -1 && iel2 < 0 && info2 < 0,
                    "MaterialPropertyCoefficient should only be used for exterior "
                    "(single-sided) boundaries!");
#endif
        return mesh.GetAttribute(iel1);
      }
    }
    else if constexpr (ElemType == MeshElementType::ELEMENT)
    {
      MFEM_ASSERT(
          T.ElementType == mfem::ElementTransformation::ELEMENT,
          "Invalid usage of MaterialPropertyCoefficient for given MeshElementType!");
      return T.Attribute;
    }
    else if constexpr (ElemType == MeshElementType::BDR_ELEMENT)
    {
      MFEM_ASSERT(
          T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
          "Invalid usage of MaterialPropertyCoefficient for given MeshElementType!");
      int i, o, iel1, iel2;
      const mfem::Mesh &mesh = *T.mesh;
      mesh.GetBdrElementFace(T.ElementNo, &i, &o);
      mesh.GetFaceElements(i, &iel1, &iel2);
#ifdef MFEM_DEBUG
      int info1, info2, nc;
      mesh.GetFaceInfos(i, &info1, &info2, &nc);
      MFEM_VERIFY(nc == -1 && iel2 < 0 && info2 < 0,
                  "MaterialPropertyCoefficient should only be used for exterior "
                  "(single-sided) boundaries!");
#endif
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
    if constexpr (MatType == MaterialPropertyType::INV_PERMEABILITY)
    {
      K = mat_op.GetInvPermeability(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::PERMITTIVITY_REAL)
    {
      K = mat_op.GetPermittivityReal(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::PERMITTIVITY_IMAG)
    {
      K = mat_op.GetPermittivityImag(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::PERMITTIVITY_ABS)
    {
      K = mat_op.GetPermittivityAbs(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::CONDUCTIVITY)
    {
      K = mat_op.GetConductivity(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::INV_LONDON_DEPTH)
    {
      K = mat_op.GetInvLondonDepth(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::INV_Z0)
    {
      K = mat_op.GetInvImpedance(GetAttribute(T));
    }
    else if constexpr (MatType == MaterialPropertyType::INV_PERMEABILITY_C0)
    {
      const int attr = GetAttribute(T);
      K.SetSize(height, width);
      Mult(mat_op.GetInvPermeability(attr), mat_op.GetLightSpeed(attr), K);
    }
    else
    {
      MFEM_ABORT("MaterialPropertyCoefficient::Eval() is not implemented for this "
                 "material property type!");
    }
    K *= coef;
  }
};

// Base class for coefficients which need to evaluate a GridFunction in a domain element
// attached to a boundary element, or both domain elements on either side for internal
// boundaries.
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
  BdrGridFunctionCoefficient(mfem::ParMesh &mesh, const std::map<int, int> &local_to_shared)
    : mesh(mesh), local_to_shared(local_to_shared)
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
  mutable mfem::Vector C1, W, VU, VL, nor;

public:
  BdrCurrentVectorCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op)
    : mfem::VectorCoefficient(gf.ParFESpace()->GetParMesh()->SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                 op.GetLocalToSharedFaceMap()),
      B(gf), mat_op(op), C1(gf.VectorDim()), W(gf.VectorDim()), VU(gf.VectorDim()),
      VL(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(vdim == 3, "BdrJVectorCoefficient expects a mesh in 3D space!");
    mfem::ElementTransformation *T1, *T2;
    GetElementTransformations(T, ip, T1, T2, &C1);

    // For interior faces, compute J_s = -n x H = -n x μ⁻¹(B1 - B2), where B1 (B2) is B in
    // el1 (el2) and n points out from el1.
    B.GetVectorValue(*T1, T1->GetIntPoint(), W);
    mat_op.GetInvPermeability(T1->Attribute).Mult(W, VU);
    if (T2)
    {
      // Double-sided, not a true boundary.
      B.GetVectorValue(*T2, T2->GetIntPoint(), W);
      mat_op.GetInvPermeability(T2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
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
  mutable mfem::Vector C1, W, VU, VL, nor;

public:
  BdrChargeCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      op.GetLocalToSharedFaceMap()),
      E(gf), mat_op(op), C1(gf.VectorDim()), W(gf.VectorDim()), VU(gf.VectorDim()),
      VL(gf.VectorDim()), nor(gf.VectorDim())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
    GetElementTransformations(T, ip, T1, T2, &C1);

    // For interior faces, compute D ⋅ n = ε (E1 - E2) ⋅ n, where E1 (E2) is E in el1 (el2)
    // to get a single-valued function.
    E.GetVectorValue(*T1, T1->GetIntPoint(), W);
    mat_op.GetPermittivityReal(T1->Attribute).Mult(W, VU);
    if (T2)
    {
      E.GetVectorValue(*T2, T2->GetIntPoint(), W);
      mat_op.GetPermittivityReal(T2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into el1.
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
  mutable mfem::Vector V, VL, nor;

public:
  BdrFluxCoefficient(const mfem::ParGridFunction &gf, mfem::Vector d,
                     const std::map<int, int> &local_to_shared)
    : mfem::Coefficient(),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(), local_to_shared), B(gf),
      dir(std::move(d)), V(gf.VectorDim()), VL(gf.VectorDim()), nor(gf.VectorDim())
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
    B.GetVectorValue(*T1, T1->GetIntPoint(), V);
    if (T2)
    {
      B.GetVectorValue(*T2, T2->GetIntPoint(), VL);
      V += VL;
      V *= 0.5;
    }

    // Orient sign with the global direction.
    GetNormal(T, ip, nor);
    return (dir * nor < 0.0) ? -(V * nor) : V * nor;
  }
};

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
  mutable mfem::Vector C1, V, nor;

  int Initialize(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                 mfem::Vector &V)
  {
    // Get neighboring elements.
    mfem::ElementTransformation *T1, *T2;
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
                                 mfem::Vector s)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      op.GetLocalToSharedFaceMap()),
      E(gf), mat_op(op), ts(ti), epsilon(ei), side(std::move(s)), C1(gf.VectorDim()),
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
  Initialize(T, ip, V);
  GetNormal(T, ip, nor);

  // Substrate-air interface: 0.5 * t * (ϵ_SA * |E_t|² + 1 / ϵ_MS * |E_n|²) .
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

  // No specific interface, use full field evaluation: 0.5 * t * ϵ * |E|² .
  return 0.5 * ts * epsilon * (V * V);
}

// Computes the flux, μ⁻¹ ∇ × X, of a field, X, where X can be the electric field E, or the
// magnetic vector potential A.
class CurlFluxCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParGridFunction &X;
  MaterialPropertyCoefficient<MaterialPropertyType::INV_PERMEABILITY> coef;

  mfem::DenseMatrix muinv;
  mfem::Vector curl;

public:
  CurlFluxCoefficient(const mfem::ParGridFunction &pgf, const MaterialOperator &op)
    : mfem::VectorCoefficient(pgf.ParFESpace()->GetParMesh()->SpaceDimension()), X(pgf),
      coef(op, 1.0), muinv(3), curl(3)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    V.SetSize(3);
    coef.Eval(muinv, T, ip);
    X.GetCurl(T, curl);
    muinv.Mult(curl, V);
  }
};

// Computes the flux, ϵ ∇ ϕ, of the electrostatic potential ϕ.
class GradFluxCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParGridFunction &phi;
  MaterialPropertyCoefficient<MaterialPropertyType::PERMITTIVITY_REAL> coef;

  mfem::Vector grad;
  mfem::DenseMatrix eps;

public:
  GradFluxCoefficient(const mfem::ParGridFunction &pgf, const MaterialOperator &op)
    : mfem::VectorCoefficient(pgf.ParFESpace()->GetParMesh()->SpaceDimension()), phi(pgf),
      coef(op, 1.0), grad(3), eps(3)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    V.SetSize(3);
    coef.Eval(eps, T, ip);
    phi.GetGradient(T, grad);
    eps.Mult(grad, V);
  }
};

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
  mutable mfem::Vector V;

  double GetLocalEnergyDensity(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip, int attr);

public:
  EnergyDensityCoefficient(const GridFunctionType &gf, const MaterialOperator &op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      op.GetLocalToSharedFaceMap()),
      U(gf), mat_op(op), V(gf.ParFESpace()->GetParMesh()->SpaceDimension())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      T.SetIntPoint(&ip);
      return GetLocalEnergyDensity(T, ip, T.Attribute);
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
  BdrFieldVectorCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op)
    : mfem::VectorCoefficient(gf.ParFESpace()->GetParMesh()->SpaceDimension()),
      BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                 op.GetLocalToSharedFaceMap()),
      U(gf), mat_op(op)
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
  BdrFieldCoefficient(const mfem::ParGridFunction &gf, const MaterialOperator &op)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*gf.ParFESpace()->GetParMesh(),
                                                      op.GetLocalToSharedFaceMap()),
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

// Wraps a mfem::MatrixCoefficient to compute a scalar coefficient as nᵀ M n. Only works
// for square matrix coefficients of size equal to the spatial dimension.
class NormalProjectedCoefficient : public mfem::Coefficient
{
  std::unique_ptr<mfem::MatrixCoefficient> c;
  mutable mfem::DenseMatrix K;
  mutable mfem::Vector nor;

public:
  NormalProjectedCoefficient(std::unique_ptr<mfem::MatrixCoefficient> &&coef)
    : mfem::Coefficient(), c(std::move(coef)), K(c->GetHeight(), c->GetWidth()),
      nor(c->GetHeight())
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
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
  mutable mfem::Vector U;

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
  SumVectorCoefficient(int d) : mfem::VectorCoefficient(d), U(d) {}

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
  mutable mfem::DenseMatrix M;

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
  SumMatrixCoefficient(int d) : mfem::MatrixCoefficient(d), M(d) {}
  SumMatrixCoefficient(int h, int w) : mfem::MatrixCoefficient(h, w), M(h, w) {}

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
  }
  fip.z = ip.z;
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
  MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
              "Unexpected element type in BdrGridFunctionCoefficient!");
  MFEM_ASSERT(&mesh == T.mesh, "Invalid mesh for BdrGridFunctionCoefficient!");
  int i, o;
  int iel1, iel2, info1, info2;
  mesh.GetBdrElementFace(T.ElementNo, &i, &o);
  mesh.GetFaceElements(i, &iel1, &iel2);
  mesh.GetFaceInfos(i, &info1, &info2);  // XX TODO: Nonconforming support

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

}  // namespace palace

#endif  // PALACE_FEM_COEFFICIENT_HPP
