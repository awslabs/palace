// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_COEFFICIENT_HPP
#define PALACE_FEM_COEFFICIENT_HPP

#include <complex>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "fem/gridfunction.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/geodata.hpp"
#include "utils/labels.hpp"

// XX TODO: Add bulk element Eval() overrides to speed up postprocessing (also needed in
//          mfem::DataCollection classes.

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
// boundaries, with optional scaling factor.
class BdrGridFunctionCoefficient
{
protected:
  // XX TODO: For thread-safety (multiple threads evaluating a coefficient simultaneously),
  //          the FET, FET.Elem1, and FET.Elem2 objects cannot be shared.
  const mfem::ParMesh &mesh;
  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  const double scaling;  // scaling factor used for unit conversions

  bool GetBdrElementNeighborTransformations(int i, const mfem::IntegrationPoint &ip)
  {
    // Get the element transformations neighboring the element, and optionally set the
    // integration point too.
    return GetBdrElementNeighborTransformations(i, mesh, FET, T1, T2, &ip);
  }

public:
  BdrGridFunctionCoefficient(const mfem::ParMesh &mesh, double scaling = 1.0)
    : mesh(mesh), scaling(scaling)
  {
  }

  // For a boundary element, return the element transformation objects for the neighboring
  // domain elements. FET.Elem2 may be nullptr if the boundary is a true one-sided boundary,
  // but if it is shared with another subdomain then it will be populated. Expects
  // ParMesh::ExchangeFaceNbrData has been called already.
  static bool GetBdrElementNeighborTransformations(
      int i, const mfem::ParMesh &mesh, mfem::FaceElementTransformations &FET,
      mfem::IsoparametricTransformation &T1, mfem::IsoparametricTransformation &T2,
      const mfem::IntegrationPoint *ip = nullptr);

  // Return normal vector to the boundary element at an integration point. For a face
  // element, the normal points out of the element (from element 1 into element 2, if it
  // exists). This convention can be flipped with the optional parameter. It is assumed
  // that the element transformation has already been configured at the integration point
  // of interest.
  static void GetNormal(mfem::ElementTransformation &T, mfem::Vector &normal,
                        bool invert = false)
  {
    MFEM_ASSERT(normal.Size() == T.GetSpaceDim(),
                "Size mismatch for normal vector (space dimension = " << T.GetSpaceDim()
                                                                      << ")!");
    mfem::CalcOrtho(T.Jacobian(), normal);
    normal /= invert ? -normal.Norml2() : normal.Norml2();
  }
};

// Computes surface current Jₛ = n x H = n x μ⁻¹ B on boundaries from B as a vector grid
// function where n is an inward normal (computes -n x H for outward normal n). For a
// two-sided internal boundary, the contributions from both sides add.
class BdrSurfaceCurrentVectorCoefficient : public mfem::VectorCoefficient,
                                           public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &B;
  const MaterialOperator &mat_op;

public:
  BdrSurfaceCurrentVectorCoefficient(const mfem::ParGridFunction &B,
                                     const MaterialOperator &mat_op, double scaling = 1.0)
    : mfem::VectorCoefficient(B.VectorDim()),
      BdrGridFunctionCoefficient(*B.ParFESpace()->GetParMesh(), scaling), B(B),
      mat_op(mat_op)
  {
  }

  using mfem::VectorCoefficient::Eval;
  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in BdrSurfaceCurrentVectorCoefficient!");
    bool ori = GetBdrElementNeighborTransformations(T.ElementNo, ip);

    // For interior faces, compute Jₛ = n x H = n x μ⁻¹ (B1 - B2), where B1 (B2) is B in
    // element 1 (element 2) and n points into element 1.
    double W_data[3], VU_data[3];
    mfem::Vector W(W_data, vdim), VU(VU_data, vdim);
    B.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), W);
    mat_op.GetInvPermeability(FET.Elem1->Attribute).Mult(W, VU);
    if (FET.Elem2)
    {
      // Double-sided, not a true boundary. Add result with opposite normal.
      double VL_data[3];
      mfem::Vector VL(VL_data, vdim);
      B.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
      mat_op.GetInvPermeability(FET.Elem2->Attribute).Mult(W, VL);
      VU -= VL;
    }

    // Orient with normal pointing into element 1.
    double normal_data[3];
    mfem::Vector normal(normal_data, vdim);
    GetNormal(T, normal, ori);
    V.SetSize(vdim);
    linalg::Cross3(normal, VU, V);

    V *= scaling;
  }
};

// Computes the flux Φₛ = F ⋅ n with F = B or ε E on interior boundary elements using B or
// E given as a vector grid function. For a two-sided internal boundary, the contributions
// from both sides can either add or be averaged.
template <SurfaceFlux Type>
class BdrSurfaceFluxCoefficient : public mfem::Coefficient,
                                  public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction *E, *B;
  const MaterialOperator &mat_op;
  bool two_sided;
  const mfem::Vector &x0;
  void GetLocalFlux(mfem::ElementTransformation &T, mfem::Vector &V) const;

public:
  BdrSurfaceFluxCoefficient(const mfem::ParGridFunction *E, const mfem::ParGridFunction *B,
                            const MaterialOperator &mat_op, bool two_sided,
                            const mfem::Vector &x0, double scaling = 1.0)
    : mfem::Coefficient(),
      BdrGridFunctionCoefficient(
          E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh(), scaling),
      E(E), B(B), mat_op(mat_op), two_sided(two_sided), x0(x0)
  {
    MFEM_VERIFY((E || (Type != SurfaceFlux::ELECTRIC && Type != SurfaceFlux::POWER)) &&
                    (B || (Type != SurfaceFlux::MAGNETIC && Type != SurfaceFlux::POWER)),
                "Missing E or B field grid function for surface flux coefficient!");
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in BdrSurfaceFluxCoefficient!");
    bool ori = GetBdrElementNeighborTransformations(T.ElementNo, ip);

    // For interior faces, compute either F ⋅ n as the average or by adding the
    // contributions from opposite sides with opposite normals.
    const int vdim = T.GetSpaceDim();
    double VU_data[3];
    mfem::Vector VU(VU_data, vdim);
    GetLocalFlux(*FET.Elem1, VU);
    if (FET.Elem2)
    {
      // Double-sided, not a true boundary.
      double VL_data[3];
      mfem::Vector VL(VL_data, vdim);
      GetLocalFlux(*FET.Elem2, VL);
      if (two_sided)
      {
        // Add result with opposite normal. This only happens when crack_bdr_elements =
        // false (two_sided = true doesn't make sense for an internal boundary without an
        // associated BC).
        VU -= VL;
      }
      else
      {
        // Take the average of the values on both sides.
        add(0.5, VU, VL, VU);
      }
    }

    // Dot with normal direction and assign appropriate sign. The normal is oriented to
    // point into element 1.
    double normal_data[3];
    mfem::Vector normal(normal_data, vdim);
    GetNormal(T, normal, ori);
    double flux = VU * normal;
    if (two_sided)
    {
      return flux;
    }
    else
    {
      // Orient outward from the surface with the given center.
      double x_data[3];
      mfem::Vector x(x_data, vdim);
      T.Transform(ip, x);
      x -= x0;
      return (x * normal < 0.0) ? -flux : flux;
    }
  }
};

template <>
inline void BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>::GetLocalFlux(
    mfem::ElementTransformation &T, mfem::Vector &V) const
{
  // Flux D.
  double W_data[3];
  mfem::Vector W(W_data, T.GetSpaceDim());
  E->GetVectorValue(T, T.GetIntPoint(), W);
  mat_op.GetPermittivityReal(T.Attribute).Mult(W, V);
  V *= scaling;
}

template <>
inline void BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>::GetLocalFlux(
    mfem::ElementTransformation &T, mfem::Vector &V) const
{
  // Flux B.
  B->GetVectorValue(T, T.GetIntPoint(), V);
  V *= scaling;
}

template <>
inline void
BdrSurfaceFluxCoefficient<SurfaceFlux::POWER>::GetLocalFlux(mfem::ElementTransformation &T,
                                                            mfem::Vector &V) const
{
  // Flux E x H = E x μ⁻¹ B.
  double W1_data[3], W2_data[3];
  mfem::Vector W1(W1_data, T.GetSpaceDim()), W2(W2_data, T.GetSpaceDim());
  B->GetVectorValue(T, T.GetIntPoint(), W1);
  mat_op.GetInvPermeability(T.Attribute).Mult(W1, W2);
  E->GetVectorValue(T, T.GetIntPoint(), W1);
  V.SetSize(W1.Size());
  linalg::Cross3(W1, W2, V);
  V *= scaling;
}

// Computes a single-valued α Eᵀ E on boundaries from E given as a vector grid function.
// Uses the neighbor element on a user specified side to compute a single-sided value for
// potentially discontinuous solutions for an interior boundary element. The four cases
// correspond to a generic interface vs. specializations for metal-air, metal-substrate,
// and substrate-air interfaces following:
//   J. Wenner et al., Surface loss simulations of superconducting coplanar waveguide
//     resonators, Appl. Phys. Lett. (2011).
template <InterfaceDielectric Type>
class InterfaceDielectricCoefficient : public mfem::Coefficient,
                                       public BdrGridFunctionCoefficient
{
private:
  const GridFunction &E;
  const MaterialOperator &mat_op;
  const double t_i, epsilon_i;

  void Initialize(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
                  mfem::Vector *normal)
  {
    // Get neighboring elements and the normal vector, oriented to point into element 1.
    MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in InterfaceDielectricCoefficient!");
    bool ori = GetBdrElementNeighborTransformations(T.ElementNo, ip);
    if (normal)
    {
      GetNormal(T, *normal, ori);
    }
  }

  int GetLocalVectorValue(const mfem::ParGridFunction &U, mfem::Vector &V,
                          bool vacuum_side) const
  {
    constexpr double threshold = 1.0 - 1.0e-6;
    const bool use_elem1 =
        ((vacuum_side && mat_op.GetLightSpeedMax(FET.Elem1->Attribute) >= threshold) ||
         (!vacuum_side && mat_op.GetLightSpeedMax(FET.Elem1->Attribute) < threshold));
    const bool use_elem2 =
        (FET.Elem2 &&
         ((vacuum_side && mat_op.GetLightSpeedMax(FET.Elem2->Attribute) >= threshold) ||
          (!vacuum_side && mat_op.GetLightSpeedMax(FET.Elem2->Attribute) < threshold)));
    if (use_elem1)
    {
      U.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
      if (use_elem2)
      {
        // Double-sided, not a true boundary. Just average the solution from both sides.
        double W_data[3];
        mfem::Vector W(W_data, V.Size());
        U.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
        add(0.5, V, W, V);
      }
      return FET.Elem1->Attribute;
    }
    else if (use_elem2)
    {
      U.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), V);
      return FET.Elem2->Attribute;
    }
    else
    {
      return 0;
    }
  }

public:
  InterfaceDielectricCoefficient(const GridFunction &E, const MaterialOperator &mat_op,
                                 double t_i, double epsilon_i)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*E.ParFESpace()->GetParMesh()), E(E),
      mat_op(mat_op), t_i(t_i), epsilon_i(epsilon_i)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;
};

template <>
inline double InterfaceDielectricCoefficient<InterfaceDielectric::DEFAULT>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution. Don't use lightspeed detection for differentiating side.
  auto GetLocalVectorValueDefault = [this](const mfem::ParGridFunction &U, mfem::Vector &V)
  {
    U.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
    if (FET.Elem2)
    {
      // Double-sided, not a true boundary. Just average the field solution from both sides.
      double W_data[3];
      mfem::Vector W(W_data, V.Size());
      U.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
      add(0.5, V, W, V);
    }
  };
  double V_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim());
  Initialize(T, ip, nullptr);
  GetLocalVectorValueDefault(E.Real(), V);
  double V2 = V * V;
  if (E.HasImag())
  {
    GetLocalVectorValueDefault(E.Imag(), V);
    V2 += V * V;
  }

  // No specific interface, use full field evaluation: 0.5 * t * ε * |E|² .
  return 0.5 * t_i * epsilon_i * V2;
}

template <>
inline double InterfaceDielectricCoefficient<InterfaceDielectric::MA>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution on air (vacuum) side and neighboring element attribute.
  double V_data[3], normal_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim()), normal(normal_data, T.GetSpaceDim());
  Initialize(T, ip, &normal);
  int attr = GetLocalVectorValue(E.Real(), V, true);
  if (attr <= 0)
  {
    return 0.0;
  }
  double Vn = V * normal;
  double Vn2 = Vn * Vn;
  if (E.HasImag())
  {
    GetLocalVectorValue(E.Imag(), V, true);
    Vn = V * normal;
    Vn2 += Vn * Vn;
  }

  // Metal-air interface: 0.5 * t / ε_MA * |E_n|² .
  return 0.5 * (t_i / epsilon_i) * Vn2;
}

template <>
inline double InterfaceDielectricCoefficient<InterfaceDielectric::MS>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution on substrate side and neighboring element attribute.
  double V_data[3], W_data[3], normal_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim()), W(W_data, T.GetSpaceDim()),
      normal(normal_data, T.GetSpaceDim());
  Initialize(T, ip, &normal);
  int attr = GetLocalVectorValue(E.Real(), V, false);
  if (attr <= 0)
  {
    return 0.0;
  }
  mat_op.GetPermittivityReal(attr).Mult(V, W);
  double Vn = W * normal;
  double Vn2 = Vn * Vn;
  if (E.HasImag())
  {
    GetLocalVectorValue(E.Imag(), V, false);
    mat_op.GetPermittivityReal(attr).Mult(V, W);
    Vn = W * normal;
    Vn2 += Vn * Vn;
  }

  // Metal-substrate interface: 0.5 * t / ε_MS * |(ε_S E)_n|² .
  return 0.5 * (t_i / epsilon_i) * Vn2;
}

template <>
inline double InterfaceDielectricCoefficient<InterfaceDielectric::SA>::Eval(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
  // Get single-sided solution on air side and neighboring element attribute.
  double V_data[3], normal_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim()), normal(normal_data, T.GetSpaceDim());
  Initialize(T, ip, &normal);
  int attr = GetLocalVectorValue(E.Real(), V, true);
  if (attr <= 0)
  {
    return 0.0;
  }
  double Vn = V * normal;
  V.Add(-Vn, normal);
  double Vn2 = Vn * Vn;
  double Vt2 = V * V;
  if (E.HasImag())
  {
    GetLocalVectorValue(E.Imag(), V, true);
    Vn = V * normal;
    V.Add(-Vn, normal);
    Vn2 += Vn * Vn;
    Vt2 += V * V;
  }

  // Substrate-air interface: 0.5 * t * (ε_SA * |E_t|² + 1 / ε_SA * |E_n|²) .
  return 0.5 * t_i * ((epsilon_i * Vt2) + (Vn2 / epsilon_i));
}

// Helper for EnergyDensityCoefficient.
enum class EnergyDensityType
{
  ELECTRIC,
  MAGNETIC
};

// Returns the local energy density evaluated as 1/2 Dᴴ E or 1/2 Hᴴ B for real-valued
// material coefficients. For internal boundary elements, the solution is averaged across
// the interface.
template <EnergyDensityType Type>
class EnergyDensityCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const GridFunction &U;
  const MaterialOperator &mat_op;
  double GetLocalEnergyDensity(mfem::ElementTransformation &T) const;

public:
  EnergyDensityCoefficient(const GridFunction &U, const MaterialOperator &mat_op,
                           double scaling = 1.0)
    : mfem::Coefficient(),
      BdrGridFunctionCoefficient(*U.ParFESpace()->GetParMesh(), scaling), U(U),
      mat_op(mat_op)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      return GetLocalEnergyDensity(T);
    }
    else if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT)
    {
      // Get neighboring elements.
      GetBdrElementNeighborTransformations(T.ElementNo, ip);

      // For interior faces, compute the average value.
      if (FET.Elem2)
      {
        return 0.5 *
               (GetLocalEnergyDensity(*FET.Elem1) + GetLocalEnergyDensity(*FET.Elem2));
      }
      else
      {
        return GetLocalEnergyDensity(*FET.Elem1);
      }
    }
    MFEM_ABORT("Unsupported element type in EnergyDensityCoefficient!");
    return 0.0;
  }
};

template <>
inline double EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>::GetLocalEnergyDensity(
    mfem::ElementTransformation &T) const
{
  // Only the real part of the permittivity contributes to the energy (imaginary part
  // cancels out in the inner product due to symmetry).
  double V_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim());
  U.Real().GetVectorValue(T, T.GetIntPoint(), V);
  double dot = mat_op.GetPermittivityReal(T.Attribute).InnerProduct(V, V);
  if (U.HasImag())
  {
    U.Imag().GetVectorValue(T, T.GetIntPoint(), V);
    dot += mat_op.GetPermittivityReal(T.Attribute).InnerProduct(V, V);
  }
  return 0.5 * dot * scaling;
}

template <>
inline double EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>::GetLocalEnergyDensity(
    mfem::ElementTransformation &T) const
{
  double V_data[3];
  mfem::Vector V(V_data, T.GetSpaceDim());
  U.Real().GetVectorValue(T, T.GetIntPoint(), V);
  double dot = mat_op.GetInvPermeability(T.Attribute).InnerProduct(V, V);
  if (U.HasImag())
  {
    U.Imag().GetVectorValue(T, T.GetIntPoint(), V);
    dot += mat_op.GetInvPermeability(T.Attribute).InnerProduct(V, V);
  }
  return 0.5 * dot * scaling;
}

// Compute time-averaged Poynting vector Re{E x H⋆}, without the typical factor of 1/2. For
// internal boundary elements, the solution is taken as the average.
class PoyntingVectorCoefficient : public mfem::VectorCoefficient,
                                  public BdrGridFunctionCoefficient
{
private:
  const GridFunction &E, &B;
  const MaterialOperator &mat_op;

  void GetLocalPower(mfem::ElementTransformation &T, mfem::Vector &V) const
  {
    double W1_data[3], W2_data[3];
    mfem::Vector W1(W1_data, T.GetSpaceDim()), W2(W2_data, T.GetSpaceDim());
    B.Real().GetVectorValue(T, T.GetIntPoint(), W1);
    mat_op.GetInvPermeability(T.Attribute).Mult(W1, W2);
    E.Real().GetVectorValue(T, T.GetIntPoint(), W1);
    V.SetSize(vdim);
    linalg::Cross3(W1, W2, V);
    if (E.HasImag())
    {
      B.Imag().GetVectorValue(T, T.GetIntPoint(), W1);
      mat_op.GetInvPermeability(T.Attribute).Mult(W1, W2);
      E.Imag().GetVectorValue(T, T.GetIntPoint(), W1);
      linalg::Cross3(W1, W2, V, true);
    }
    V *= scaling;
  }

public:
  PoyntingVectorCoefficient(const GridFunction &E, const GridFunction &B,
                            const MaterialOperator &mat_op, double scaling = 1.0)
    : mfem::VectorCoefficient(E.VectorDim()),
      BdrGridFunctionCoefficient(*E.ParFESpace()->GetParMesh(), scaling), E(E), B(B),
      mat_op(mat_op)
  {
  }

  using mfem::VectorCoefficient::Eval;
  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    if (T.ElementType == mfem::ElementTransformation::ELEMENT)
    {
      GetLocalPower(T, V);
      return;
    }
    else if (T.ElementType == mfem::ElementTransformation::BDR_ELEMENT)
    {
      // Get neighboring elements.
      GetBdrElementNeighborTransformations(T.ElementNo, ip);

      // For interior faces, compute the value on the desired side.
      GetLocalPower(*FET.Elem1, V);
      if (FET.Elem2)
      {
        double W_data[3];
        mfem::Vector W(W_data, V.Size());
        GetLocalPower(*FET.Elem2, W);
        add(0.5, V, W, V);
      }
      return;
    }
    MFEM_ABORT("Unsupported element type in PoyntingVectorCoefficient!");
  }
};

// Returns the local vector field evaluated on a boundary element. For internal boundary
// elements the solution is the average.
class BdrFieldVectorCoefficient : public mfem::VectorCoefficient,
                                  public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &U;

public:
  BdrFieldVectorCoefficient(const mfem::ParGridFunction &U)
    : mfem::VectorCoefficient(U.VectorDim()),
      BdrGridFunctionCoefficient(*U.ParFESpace()->GetParMesh()), U(U)
  {
  }

  using mfem::VectorCoefficient::Eval;
  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in BdrFieldVectorCoefficient!");
    GetBdrElementNeighborTransformations(T.ElementNo, ip);

    // For interior faces, compute the average.
    U.GetVectorValue(*FET.Elem1, FET.Elem1->GetIntPoint(), V);
    if (FET.Elem2)
    {
      double W_data[3];
      mfem::Vector W(W_data, V.Size());
      U.GetVectorValue(*FET.Elem2, FET.Elem2->GetIntPoint(), W);
      add(0.5, V, W, V);
    }
  }
};

// Returns the local scalar field evaluated on a boundary element. For internal boundary
// elements the solution is the average.
class BdrFieldCoefficient : public mfem::Coefficient, public BdrGridFunctionCoefficient
{
private:
  const mfem::ParGridFunction &U;

public:
  BdrFieldCoefficient(const mfem::ParGridFunction &U)
    : mfem::Coefficient(), BdrGridFunctionCoefficient(*U.ParFESpace()->GetParMesh()), U(U)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    // Get neighboring elements.
    MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in BdrFieldCoefficient!");
    GetBdrElementNeighborTransformations(T.ElementNo, ip);

    // For interior faces, compute the average.
    if (FET.Elem2)
    {
      return 0.5 * (U.GetValue(*FET.Elem1, FET.Elem1->GetIntPoint()),
                    U.GetValue(*FET.Elem2, FET.Elem2->GetIntPoint()));
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

  using mfem::VectorCoefficient::Eval;
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

template <typename Coefficient>
class RestrictedCoefficient : public Coefficient
{
private:
  mfem::Array<int> attr_marker;

public:
  template <typename... T>
  RestrictedCoefficient(const mfem::Array<int> &attr_list, T &&...args)
    : Coefficient(std::forward<T>(args)...),
      attr_marker(mesh::AttrToMarker(attr_list.Size() ? attr_list.Max() : 0, attr_list))
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    return (T.Attribute > attr_marker.Size() || !attr_marker[T.Attribute - 1])
               ? 0.0
               : Coefficient::Eval(T, ip);
  }
};

template <typename Coefficient>
class RestrictedVectorCoefficient : public Coefficient
{
private:
  mfem::Array<int> attr_marker;

public:
  template <typename... T>
  RestrictedVectorCoefficient(const mfem::Array<int> &attr_list, T &&...args)
    : Coefficient(std::forward<T>(args)...),
      attr_marker(mesh::AttrToMarker(attr_list.Size() ? attr_list.Max() : 0, attr_list))
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    if (T.Attribute > attr_marker.Size() || !attr_marker[T.Attribute - 1])
    {
      V.SetSize(this->vdim);
      V = 0.0;
    }
    else
    {
      Coefficient::Eval(V, T, ip);
    }
  }
};

template <typename Coefficient>
class RestrictedMatrixCoefficient : public Coefficient
{
private:
  mfem::Array<int> attr_marker;

public:
  template <typename... T>
  RestrictedMatrixCoefficient(const mfem::Array<int> &attr_list, T &&...args)
    : Coefficient(std::forward<T>(args)...),
      attr_marker(mesh::AttrToMarker(attr_list.Size() ? attr_list.Max() : 0, attr_list))
  {
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    if (T.Attribute > attr_marker.Size() || !attr_marker[T.Attribute - 1])
    {
      K.SetSize(this->height, this->width);
      K = 0.0;
    }
    else
    {
      Coefficient::Eval(K, T, ip);
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

  using mfem::VectorCoefficient::Eval;
  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    double U_data[3];
    mfem::Vector U(U_data, vdim);
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
    double M_data[9];
    mfem::DenseMatrix M(M_data, height, width);
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
