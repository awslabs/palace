// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Wrapper for MFEM's ParFiniteElementSpace class, with extensions for Palace.
//
class FiniteElementSpace
{
private:
  // Underlying MFEM object.
  mfem::ParFiniteElementSpace fespace;

  // Reference to the underlying mesh object (not owned).
  Mesh &mesh;

  // Members for constructing libCEED operators.
  mutable ceed::CeedObjectMap<CeedBasis> basis;
  mutable ceed::CeedObjectMap<CeedElemRestriction> restr, interp_restr, interp_range_restr;

  // Temporary storage for operator applications.
  mutable ComplexVector tx, lx, ly;

  // Members for discrete interpolators from an auxiliary space to a primal space.
  mutable const FiniteElementSpace *aux_fespace;
  mutable std::unique_ptr<Operator> G;

  bool HasUniqueInterpRestriction(const mfem::FiniteElement &fe) const
  {
    // For interpolation operators and tensor-product elements, we need native (not
    // lexicographic) ordering.
    const mfem::TensorBasisElement *tfe =
        dynamic_cast<const mfem::TensorBasisElement *>(&fe);
    return (tfe && tfe->GetDofMap().Size() > 0 &&
            fe.GetRangeType() != mfem::FiniteElement::VECTOR);
  }

  bool HasUniqueInterpRangeRestriction(const mfem::FiniteElement &fe) const
  {
    // The range restriction for interpolation operators needs to use a special
    // DofTransformation (not equal to the transpose of the domain restriction).
    if (mesh.Dimension() < 3)
    {
      return false;
    }
    const auto geom = fe.GetGeomType();
    const auto *dof_trans = fespace.FEColl()->DofTransformationForGeometry(geom);
    return (dof_trans && !dof_trans->IsIdentity());
  }

  const Operator &BuildDiscreteInterpolator() const;

public:
  template <typename... T>
  FiniteElementSpace(Mesh &mesh, T &&...args)
    : fespace(&mesh.Get(), std::forward<T>(args)...), mesh(mesh), aux_fespace(nullptr)
  {
    ResetCeedObjects();
    tx.UseDevice(true);
    lx.UseDevice(true);
    ly.UseDevice(true);
  }
  virtual ~FiniteElementSpace() { ResetCeedObjects(); }

  const auto &Get() const { return fespace; }
  auto &Get() { return fespace; }

  operator const mfem::ParFiniteElementSpace &() const { return Get(); }
  operator mfem::ParFiniteElementSpace &() { return Get(); }

  const auto &GetFEColl() const { return *Get().FEColl(); }
  auto &GetFEColl() { return *Get().FEColl(); }

  const auto &GetMesh() const { return mesh; }
  auto &GetMesh() { return mesh; }

  const auto &GetParMesh() const { return mesh.Get(); }
  auto &GetParMesh() { return mesh.Get(); }

  auto GetVDim() const { return Get().GetVDim(); }
  auto GetVSize() const { return Get().GetVSize(); }
  auto GlobalVSize() const { return Get().GlobalVSize(); }
  auto GetTrueVSize() const { return Get().GetTrueVSize(); }
  auto GlobalTrueVSize() const { return Get().GlobalTrueVSize(); }
  auto Dimension() const { return mesh.Get().Dimension(); }
  auto SpaceDimension() const { return mesh.Get().SpaceDimension(); }
  auto GetMaxElementOrder() const { return Get().GetMaxElementOrder(); }

  const auto *GetProlongationMatrix() const { return Get().GetProlongationMatrix(); }
  const auto *GetRestrictionMatrix() const { return Get().GetRestrictionMatrix(); }

  // Return the discrete gradient, curl, or divergence matrix interpolating from the
  // auxiliary to the primal space, constructing it on the fly as necessary.
  const auto &GetDiscreteInterpolator(const FiniteElementSpace &aux_fespace_) const
  {
    if (&aux_fespace_ != aux_fespace)
    {
      G.reset();
      aux_fespace = &aux_fespace_;
    }
    return G ? *G : BuildDiscreteInterpolator();
  }

  // Return the basis object for elements of the given element geometry type.
  CeedBasis GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const;

  // Return the element restriction object for the given element set (all with the same
  // geometry type).
  CeedElemRestriction GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                             const std::vector<int> &indices) const;

  // If the space has a special element restriction for discrete interpolators, return that.
  // Otherwise return the same restriction as given by GetCeedElemRestriction.
  CeedElemRestriction GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                   const std::vector<int> &indices) const;

  // If the space has a special element restriction for the range space of discrete
  // interpolators, return that. Otherwise return the same restriction as given by
  // GetCeedElemRestriction.
  CeedElemRestriction
  GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                    const std::vector<int> &indices) const;

  // Clear the cached basis and element restriction objects owned by the finite element
  // space.
  void ResetCeedObjects();

  void Update() { ResetCeedObjects(); }

  static CeedBasis BuildCeedBasis(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                                  mfem::Geometry::Type geom);
  static CeedElemRestriction
  BuildCeedElemRestriction(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                           mfem::Geometry::Type geom, const std::vector<int> &indices,
                           bool is_interp = false, bool is_interp_range = false);

  template <typename VecType>
  auto &GetTVector() const
  {
    tx.SetSize(GetTrueVSize());
    if constexpr (std::is_same<VecType, ComplexVector>::value)
    {
      return tx;
    }
    else
    {
      return tx.Real();
    }
  }

  template <typename VecType>
  auto &GetLVector() const
  {
    lx.SetSize(GetVSize());
    if constexpr (std::is_same<VecType, ComplexVector>::value)
    {
      return lx;
    }
    else
    {
      return lx.Real();
    }
  }

  template <typename VecType>
  auto &GetLVector2() const
  {
    ly.SetSize(GetVSize());
    if constexpr (std::is_same<VecType, ComplexVector>::value)
    {
      return ly;
    }
    else
    {
      return ly.Real();
    }
  }

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return fespace.GetComm(); }
};

//
// A collection of FiniteElementSpace objects constructed on the same mesh with the ability
// to construct the prolongation operators between them as needed.
//
class FiniteElementSpaceHierarchy
{
protected:
  std::vector<std::unique_ptr<FiniteElementSpace>> fespaces;
  mutable std::vector<std::unique_ptr<Operator>> P;

  const Operator &BuildProlongationAtLevel(std::size_t l) const;

public:
  FiniteElementSpaceHierarchy() = default;
  FiniteElementSpaceHierarchy(std::unique_ptr<FiniteElementSpace> &&fespace)
  {
    AddLevel(std::move(fespace));
  }

  auto GetNumLevels() const { return fespaces.size(); }

  void AddLevel(std::unique_ptr<FiniteElementSpace> &&fespace)
  {
    fespaces.push_back(std::move(fespace));
    P.push_back(nullptr);
  }

  auto &GetFESpaceAtLevel(std::size_t l)
  {
    MFEM_ASSERT(l < GetNumLevels(),
                "Out of bounds request for finite element space at level " << l << "!");
    return *fespaces[l];
  }
  const auto &GetFESpaceAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l < GetNumLevels(),
                "Out of bounds request for finite element space at level " << l << "!");
    return *fespaces[l];
  }

  auto &GetFinestFESpace()
  {
    MFEM_ASSERT(GetNumLevels() > 0,
                "Out of bounds request for finite element space at level 0!");
    return *fespaces.back();
  }
  const auto &GetFinestFESpace() const
  {
    MFEM_ASSERT(GetNumLevels() > 0,
                "Out of bounds request for finite element space at level 0!");
    return *fespaces.back();
  }

  const auto &GetProlongationAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l + 1 < GetNumLevels(),
                "Out of bounds request for finite element space prolongation at level "
                    << l << "!");
    return P[l] ? *P[l] : BuildProlongationAtLevel(l);
  }

  std::vector<const Operator *> GetProlongationOperators() const
  {
    MFEM_ASSERT(GetNumLevels() > 1,
                "Out of bounds request for finite element space prolongation at level 0!");
    std::vector<const Operator *> P_(GetNumLevels() - 1);
    for (std::size_t l = 0; l < P_.size(); l++)
    {
      P_[l] = &GetProlongationAtLevel(l);
    }
    return P_;
  }

  const auto &GetDiscreteInterpolatorAtLevel(std::size_t l,
                                             const FiniteElementSpace &aux_fespace) const
  {
    return GetFESpaceAtLevel(l).GetDiscreteInterpolator(aux_fespace);
  }

  std::vector<const Operator *>
  GetDiscreteInterpolators(const FiniteElementSpaceHierarchy &aux_fespaces) const
  {
    std::vector<const Operator *> G_(GetNumLevels());
    G_[0] = nullptr;  // No discrete interpolator for coarsest level
    for (std::size_t l = 1; l < G_.size(); l++)
    {
      G_[l] = &GetDiscreteInterpolatorAtLevel(l, aux_fespaces.GetFESpaceAtLevel(l));
    }
    return G_;
  }
};

}  // namespace palace

#endif  // PALACE_FEM_FESPACE_HPP
