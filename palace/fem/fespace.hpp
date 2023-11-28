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

namespace fem
{

// Construct mesh data structures for assembling libCEED operators on a (mixed) mesh:
//   - Mesh element indices for threads and element geometry types.
//   - Geometry factor quadrature point data (w |J|, adj(J)^T / |J|, J / |J|) for domain
//     and boundary elements.
//   - Attributes for domain and boundary elements. The attributes are not the same as the
//     mesh element attributes as they map to a compressed (1-based) list of used
//     attributes on this MPI process.
ceed::CeedObjectMap<ceed::CeedGeomFactorData> SetUpCeedGeomFactorData(mfem::ParMesh &mesh);

// Regenerate mesh geometry factor data after a previous call to SetUpCeedGeomFactorData.
// Can be used when changing to a new quadrature rule, for example.
void UpdateCeedGeomFactorData(const mfem::ParMesh &mesh,
                              ceed::CeedObjectMap<ceed::CeedGeomFactorData> &geom_data);

}  // namespace fem

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
  const ceed::CeedObjectMap<ceed::CeedGeomFactorData> *geom_data = nullptr;

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
    const auto geom = fe.GetGeomType();
    return (DoFTransArray[geom] && !DoFTransArray[geom]->IsIdentity());
  }

public:
  using mfem::ParFiniteElementSpace::ParFiniteElementSpace;
  FiniteElementSpace(const mfem::ParFiniteElementSpace &fespace)
    : mfem::ParFiniteElementSpace(fespace)
  {
  }
  ~FiniteElementSpace() { DestroyCeedObjects(); }

  // Get the ID associated with the instance of this class. If the underlying sequence has
  // changed (due to a mesh update, for example), regenerate the ID.
  std::size_t GetId() const;

  // Operator overload for equality comparisons between two spaces.
  bool operator==(const FiniteElementSpace &fespace) const
  {
    return GetId() == fespace.GetId();
  }

  // Clear the cached basis and element restriction objects owned by the finite element
  // space.
  void DestroyCeedObjects();

  // Return the basis object for elements of the given element geometry type.
  const CeedBasis GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const;

  // Return the element restriction object for the given element set (all with the same
  // geometry type).
  const CeedElemRestriction GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                   const std::vector<int> &indices) const;

  // If the space has a special element restriction for discrete interpolators, return that.
  // Otherwise return the same restiction as given by GetCeedElemRestriction.
  const CeedElemRestriction
  GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                               const std::vector<int> &indices) const;

  // If the space has a special element restriction for the range space of discrete
  // interpolators, return that. Otherwise return the same restiction as given by
  // GetCeedElemRestriction.
  const CeedElemRestriction
  GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                    const std::vector<int> &indices) const;

  static CeedBasis BuildCeedBasis(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                                  mfem::Geometry::Type geom);
  static CeedElemRestriction
  BuildCeedElemRestriction(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                           mfem::Geometry::Type geom, const std::vector<int> &indices,
                           bool is_interp = false, bool is_interp_range = false);

  // Set and access the (not owned) geometry factor data associated with the underlying mesh
  // object for the space.
  void SetCeedGeomFactorData(const ceed::CeedObjectMap<ceed::CeedGeomFactorData> &data)
  {
    geom_data = &data;
  }
  const auto &GetCeedGeomFactorData() const
  {
    MFEM_ASSERT(
        geom_data,
        "Must call SetCeedGeomFactorData before accessing with GetCeedGeomFactorData!");
    return *geom_data;
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

  void SetCeedGeomFactorData(const ceed::CeedObjectMap<ceed::CeedGeomFactorData> &data)
  {
    for (auto &fespace : fespaces)
    {
      fespace->SetCeedGeomFactorData(data);
    }
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
