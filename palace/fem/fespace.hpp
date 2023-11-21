// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"
#include "linalg/operator.hpp"

namespace palace
{

//
// Wrapper for MFEM's ParFiniteElementSpace class, where the finite element space object
// is constructed with a unique ID associated with it. This is useful for defining equality
// operations between spaces (either different spaces on the same mesh, or the same space
// type on different meshes).
//
class FiniteElementSpace : public mfem::ParFiniteElementSpace
{
private:
  // Members used to define equality between two spaces.
  static std::size_t global_id;
  mutable std::size_t id;
  mutable long int prev_sequence;
  mutable bool init = false;

  // Members for constructing libCEED operators.
  mutable ceed::CeedObjectMap<CeedBasis> basis;
  mutable ceed::CeedObjectMap<CeedElemRestriction> restr, interp_restr, interp_range_restr;

  CeedBasis BuildCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const;
  CeedElemRestriction BuildCeedElemRestriction(Ceed ceed, const std::vector<int> &indices,
                                               bool use_bdr, bool is_interp = false,
                                               bool is_interp_range = false) const;

  bool HasUniqueInterpRestriction(const mfem::FiniteElement &fe)
  {
    // For interpolation operators and tensor-product elements, we need native (not
    // lexicographic) ordering.
    const mfem::TensorBasisElement *tfe =
        dynamic_cast<const mfem::TensorBasisElement *>(&fe);
    return (tfe && tfe->GetDofMap().Size() > 0 &&
            fe.GetRangeType() != mfem::FiniteElement::VECTOR);
  }

  bool HasUniqueInterpRangeRestriction(const mfem::FiniteElement &fe)
  {
    // The range restriction for interpolation operators needs to use a special
    // DofTransformation (not equal to the transpose of the domain restriction).
    const auto geom = fe.GetGeomType();
    return (DoFTransArray[geom] && !DoFTransArray[geom]->IsIdentity());
  }

public:
  using mfem::ParFiniteElementSpace::ParFiniteElementSpace;
  FiniteElementSpace(const mfem::ParFiniteElementSpace &fespace)
    : mfem::ParFiniteElementSpace(fespace)
  {
  }
  ~FiniteElementSpace();

  // Get the ID associated with the instance of this class. If the underlying sequence has
  // changed (due to a mesh update, for example), regenerate the ID.
  std::size_t GetId() const;

  // Operator overload for equality comparisons between two spaces.
  bool operator==(const FiniteElementSpace &fespace) const
  {
    return GetId() == fespace.GetId();
  }

  // Return the basis object for elements of the given element geometry type.
  const auto GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
  {
    const auto it = basis.find(std::make_pair(ceed, geom));
    return (it != basis.end()) ? it->second : BuildCeedBasis(ceed, geom);
  }

  // Return the element restriction object for the given element set (all with the same
  // geometry type).
  const auto GetCeedElemRestriction(Ceed ceed, const std::vector<int> &indices) const
  {
    MFEM_ASSERT(!indices.empty(),
                "Cannot create CeedElemRestriction for an empty mesh partition!");
    const auto geom = GetParMesh()->GetElementGeometry(indices[0]);
    const auto it = restr.find(std::make_pair(ceed, geom));
    return (it != restr.end()) ? it->second
                               : BuildCeedElemRestriction(ceed, indices, false);
  }

  // Return the element restriction object for the given boundary element set (all with the
  // same geometry type).
  const auto GetBdrCeedElemRestriction(Ceed ceed, const std::vector<int> &indices) const
  {
    MFEM_ASSERT(!indices.empty(),
                "Cannot create boundary CeedElemRestriction for an empty mesh partition!");
    const auto geom = GetParMesh()->GetBdrElementGeometry(indices[0]);
    const auto it = restr.find(std::make_pair(ceed, geom));
    return (it != restr.end()) ? it->second : BuildCeedElemRestriction(ceed, indices, true);
  }

  // If the space has a special element restriction for discrete interpolators, return that.
  // Otherwise return the same restiction as given by GetCeedElemRestriction.
  const auto GetInterpCeedElemRestriction(Ceed ceed, const std::vector<int> &indices) const
  {
    MFEM_ASSERT(!indices.empty(),
                "Cannot create boundary CeedElemRestriction for an empty mesh partition!");
    const auto geom = GetParMesh()->GetElementGeometry(indices[0]);
    const mfem::FiniteElement &fe = *FEColl()->FiniteElementForGeometry(geom);
    if (!HasUniqueInterpRestriction(fe))
    {
      return GetCeedElemRestriction(ceed, indices);
    }
    const auto it = interp_restr.find(std::make_pair(ceed, geom));
    return (it != interp_restr.end())
               ? it->second
               : BuildCeedElemRestriction(ceed, indices, true, true, false);
  }

  // If the space has a special element restriction for the range space of discrete
  // interpolators, return that. Otherwise return the same restiction as given by
  // GetCeedElemRestriction.
  const auto GetInterpRangeCeedElemRestriction(Ceed ceed,
                                               const std::vector<int> &indices) const
  {
    MFEM_ASSERT(!indices.empty(),
                "Cannot create boundary CeedElemRestriction for an empty mesh partition!");
    const auto geom = GetParMesh()->GetElementGeometry(indices[0]);
    const mfem::FiniteElement &fe = *FEColl()->FiniteElementForGeometry(geom);
    if (!HasUniqueInterpRangeRestriction(fe))
    {
      return GetInterpCeedElemRestriction(ceed, indices);
    }
    const auto it = interp_range_restr.find(std::make_pair(ceed, geom));
    return (it != interp_range_restr.end())
               ? it->second
               : BuildCeedElemRestriction(ceed, indices, true, true, true);
  }
};

//
// An AuxiliaryFiniteElement space is a FiniteElementSpace which allows for lazy
// construction of the interpolation operator (discrete gradient or curl) from the primal
// space to this one.
//
class AuxiliaryFiniteElementSpace : public FiniteElementSpace
{
private:
  const FiniteElementSpace &primal_fespace;
  mutable std::unique_ptr<Operator> G;

  const Operator &BuildDiscreteInterpolator() const;

public:
  template <typename... T>
  AuxiliaryFiniteElementSpace(const FiniteElementSpace &primal_fespace, T &&...args)
    : FiniteElementSpace(std::forward<T>(args)...), primal_fespace(primal_fespace)
  {
  }

  // Return the discrete gradient or discrete curl matrix interpolating from the auxiliary
  // to the primal space, constructing it on the fly as necessary.
  const auto &GetDiscreteInterpolator() const
  {
    return G ? *G : BuildDiscreteInterpolator();
  }
};

//
// A collection of FiniteElementSpace objects constructed on the same mesh with the ability
// to construct the prolongation operators between them as needed.
//
template <typename FESpace>
class BaseFiniteElementSpaceHierarchy
{
  static_assert(std::is_base_of<FiniteElementSpace, FESpace>::value,
                "A space hierarchy can only be constructed of FiniteElementSpace objects!");

protected:
  std::vector<std::unique_ptr<FESpace>> fespaces;
  mutable std::vector<std::unique_ptr<Operator>> P;

  const Operator &BuildProlongationAtLevel(std::size_t l) const;

public:
  BaseFiniteElementSpaceHierarchy() = default;
  BaseFiniteElementSpaceHierarchy(std::unique_ptr<FESpace> &&fespace)
  {
    AddLevel(std::move(fespace));
  }

  auto GetNumLevels() const { return fespaces.size(); }

  void AddLevel(std::unique_ptr<FESpace> &&fespace)
  {
    fespaces.push_back(std::move(fespace));
    P.push_back(nullptr);
  }

  auto &GetFESpaceAtLevel(std::size_t l)
  {
    MFEM_ASSERT(l >= 0 && l < GetNumLevels(),
                "Out of bounds request for finite element space at level " << l << "!");
    return *fespaces[l];
  }
  const auto &GetFESpaceAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l >= 0 && l < GetNumLevels(),
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
    MFEM_ASSERT(l >= 0 && l < GetNumLevels() - 1,
                "Out of bounds request for finite element space prolongation at level "
                    << l << "!");
    return P[l] ? *P[l] : BuildProlongationAtLevel(l);
  }

  std::vector<const Operator *> GetProlongationOperators() const
  {
    std::vector<const Operator *> P_(GetNumLevels() - 1);
    for (std::size_t l = 0; l < P_.size(); l++)
    {
      P_[l] = &GetProlongationAtLevel(l);
    }
    return P_;
  }
};

class FiniteElementSpaceHierarchy
  : public BaseFiniteElementSpaceHierarchy<FiniteElementSpace>
{
public:
  using BaseFiniteElementSpaceHierarchy<
      FiniteElementSpace>::BaseFiniteElementSpaceHierarchy;
};

//
// A special type of FiniteElementSpaceHierarchy where all members are auxiliary finite
// element spaces.
//
class AuxiliaryFiniteElementSpaceHierarchy
  : public BaseFiniteElementSpaceHierarchy<AuxiliaryFiniteElementSpace>
{
public:
  using BaseFiniteElementSpaceHierarchy<
      AuxiliaryFiniteElementSpace>::BaseFiniteElementSpaceHierarchy;

  const auto &GetDiscreteInterpolatorAtLevel(std::size_t l) const
  {
    return GetFESpaceAtLevel(l).GetDiscreteInterpolator();
  }

  std::vector<const Operator *> GetDiscreteInterpolators() const
  {
    std::vector<const Operator *> G_(GetNumLevels());
    for (std::size_t l = 0; l < G_.size(); l++)
    {
      G_[l] = &GetDiscreteInterpolatorAtLevel(l);
    }
    return G_;
  }
};

}  // namespace palace

#endif  // PALACE_FEM_FESPACE_HPP
