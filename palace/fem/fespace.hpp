// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/brokenspace.hpp"
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
  // Underlying MFEM object, owned unless this is a broken view of another space.
  std::unique_ptr<mfem::ParFiniteElementSpace> owned_fespace;
  mfem::ParFiniteElementSpace &fespace;

  // Reference to the underlying mesh object (not owned).
  Mesh &mesh;

  // Members for constructing libCEED operators. active_restr caches the single shared
  // lossy-domain or active-boundary subset for each geometry type.
  mutable ceed::CeedObjectMap<CeedBasis> basis;
  mutable ceed::CeedObjectMap<CeedElemRestriction> restr, active_restr, interp_restr,
      interp_range_restr;

  // Temporary storage for operator applications.
  mutable ComplexVector tx, lx, ly;

  // Members for discrete interpolators from an auxiliary space to a primal space.
  mutable const FiniteElementSpace *aux_fespace;
  mutable std::unique_ptr<Operator> G;

  // Data for a space which is broken (discontinuous) across interior boundaries, see
  // fem/brokenspace.hpp: the space it is a view of and the sides it was broken with,
  // element DOF overrides for the element restriction (for local element e, entries
  // [override_offsets[e], override_offsets[e + 1]) give native local DOF indices and their
  // L-vector indices in the copy block), sizes, and the broken prolongation matrix. For the
  // interpolators between broken spaces: the side label of (the elements reading) each
  // local true DOF, the original ones then the copies, and the version of each owned true
  // DOF of the underlying space read by the elements with each side label, as sorted
  // (global true DOF, side label, version) triplets.
  struct BrokenData
  {
    const FiniteElementSpace *unbroken = nullptr;
    const CrackSides *sides = nullptr;
    std::vector<int> override_offsets, override_local, override_ldof;
    int vsize = 0, tsize = 0;
    HYPRE_BigInt global_vsize = 0, global_tsize = 0;
    fem::BrokenProlongationMatrix P;
    std::vector<std::int64_t> version_label;
    std::vector<std::array<std::int64_t, 3>> resolved;
    std::vector<HYPRE_BigInt> col_starts;
  };
  std::unique_ptr<BrokenData> broken;

  void InitBroken(const CrackSides &sides);

  // For the given (global true DOF of the underlying space, side label) pairs, return the
  // global true DOF of this broken space read by the elements with that side label (the
  // version of the DOF on their side), or -1 when no element with that label reads the
  // DOF. Collective.
  std::vector<HYPRE_BigInt>
  GetVersionedTrueDofs(const std::vector<std::array<std::int64_t, 2>> &keys) const;

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
    : owned_fespace(std::make_unique<mfem::ParFiniteElementSpace>(
          &mesh.Get(), std::forward<T>(args)...)),
      fespace(*owned_fespace), mesh(mesh), aux_fespace(nullptr)
  {
    ResetCeedObjects();
    tx.UseDevice(true);
    lx.UseDevice(true);
    ly.UseDevice(true);
  }

  // Construct a view of the given space which is broken (discontinuous) across interior
  // boundaries, given the interior boundary sides of each local element (see
  // fem/brokenspace.hpp). The L-vector of the view is the L-vector of the given space
  // extended with a block of copied L-DOFs, one for each version of the true DOFs of its
  // row read by the elements on the sides of split entities other than the first, and the
  // true DOF vector is extended with the copies of the true DOFs of split entities (one per
  // side other than the first). The view shares the underlying MFEM space with the given
  // space, which must outlive the view and cannot be updated while the view exists. Only
  // element (domain) restrictions and the prolongation (a HypreParMatrix, with the global
  // offsets of the L-DOFs and true DOFs) are available for a broken space, which cannot be
  // used with MFEM assembly or for MFEM grid functions. The discrete interpolators between
  // two spaces broken with the same sides are assembled matrices (see
  // BuildBrokenInterpolator). Collective.
  FiniteElementSpace(FiniteElementSpace &fespace, const CrackSides &sides);

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
  auto GetVSize() const { return broken ? broken->vsize : Get().GetVSize(); }
  auto GlobalVSize() const { return broken ? broken->global_vsize : Get().GlobalVSize(); }
  auto GetTrueVSize() const { return broken ? broken->tsize : Get().GetTrueVSize(); }
  auto GlobalTrueVSize() const
  {
    return broken ? broken->global_tsize : Get().GlobalTrueVSize();
  }
  auto Dimension() const { return mesh.Get().Dimension(); }
  auto SpaceDimension() const { return mesh.Get().SpaceDimension(); }
  auto GetMaxElementOrder() const { return Get().GetMaxElementOrder(); }

  const Operator *GetProlongationMatrix() const
  {
    return broken ? broken->P.P.get() : Get().GetProlongationMatrix();
  }

  // The prolongation as a HypreParMatrix, and the global offsets of the L-DOFs and of the
  // true DOFs, as for mfem::ParFiniteElementSpace (for the parallel assembly of operators).
  const mfem::HypreParMatrix *Dof_TrueDof_Matrix() const
  {
    return broken ? broken->P.P.get() : Get().Dof_TrueDof_Matrix();
  }
  HYPRE_BigInt *GetDofOffsets() const
  {
    return broken ? broken->P.row_offsets.GetData() : Get().GetDofOffsets();
  }
  HYPRE_BigInt *GetTrueDofOffsets() const
  {
    return broken ? broken->P.col_offsets.GetData() : Get().GetTrueDofOffsets();
  }
  const mfem::SparseMatrix *GetRestrictionMatrix() const
  {
    MFEM_VERIFY(!broken, "Restriction matrix is not available for a broken space!");
    return Get().GetRestrictionMatrix();
  }

  bool IsBroken() const { return broken != nullptr; }

  // The space a broken space is a view of (this space when it is not broken).
  const FiniteElementSpace &GetUnbroken() const
  {
    return broken ? *broken->unbroken : *this;
  }

  // Signed element DOFs into the L-vector (see mfem::FiniteElementSpace::GetElementDofs),
  // which for a broken space read the copies of the side of the element.
  void GetElementDofs(int e, mfem::Array<int> &dofs,
                      mfem::DofTransformation &dof_trans) const;

  // Signed boundary element DOFs into the L-vector (see
  // mfem::FiniteElementSpace::GetBdrElementDofs), which for a broken space are those read
  // by the element adjacent to the boundary element: a boundary standing on one side of an
  // interior boundary reads the versions of that side.
  void GetBdrElementDofs(int be, mfem::Array<int> &dofs,
                         mfem::DofTransformation &dof_trans) const;

  // Essential true DOFs of the boundary elements with the marked attributes (see
  // mfem::FiniteElementSpace::GetEssentialTrueDofs). On a broken space, a boundary element
  // makes essential the versions read by the element adjacent to it, so a boundary standing
  // on one side of an interior boundary leaves the other side free. Collective.
  void GetEssentialTrueDofs(const mfem::Array<int> &bdr_marker,
                            mfem::Array<int> &ess_tdof_list) const;

  // Number of versions of the local true DOF t of the underlying space (1 for an unbroken
  // space), and the local true DOF of its version v (t for v = 0).
  int GetNumVersions(int t) const { return broken ? broken->P.num_versions[t] : 1; }
  int GetVersionTrueDof(int t, int v) const
  {
    return broken ? broken->P.GetVersionTrueDof(t, v) : t;
  }

  // Map from the true DOFs of a broken space to those of the underlying space, local to
  // each process: FIRST takes version 0 (the original DOFs, exact away from the interior
  // boundaries), ALL sums all versions (its transpose copies a field of the underlying
  // space to every version), and SIDES sums the traces of the two sides of a two-sided
  // interior boundary (all versions, and twice the original for true DOFs which are not
  // split, which both sides read).
  enum class VersionMap
  {
    FIRST,
    ALL,
    SIDES
  };
  std::unique_ptr<mfem::HypreParMatrix> BuildVersionMap(VersionMap type) const;

  // Build the interpolator between two spaces broken with the same sides from the assembled
  // interpolator G between the underlying spaces (on their true DOFs): for every entry
  // G(t, s) and every version v of the range DOF t, G_b((t, v), (s, w)) = G(t, s), where w
  // is the version of s read by the elements reading version v of t (the elements around an
  // entity on one side of the interior boundaries are connected through faces containing
  // the entity, so they all read the same version of each DOF of its closure, which are the
  // only DOFs an interpolated DOF depends on). This is exact also for the copies read only
  // through the constraints of hanging entities, whose values are those of a trace in the
  // master's trace space. On two meshes (the h-prolongation of a hierarchy), the side
  // labels of the fine mesh must be those of the coarse mesh (Mesh::InheritCrackSides).
  // Collective.
  static std::unique_ptr<mfem::HypreParMatrix>
  BuildBrokenInterpolator(const mfem::HypreParMatrix &G,
                          const FiniteElementSpace &trial_fespace,
                          const FiniteElementSpace &test_fespace);

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

  // Return the element restriction object for the full element index list or the shared
  // active subset of the given geometry factor data.
  CeedElemRestriction GetCeedElemRestriction(
      Ceed ceed, mfem::Geometry::Type geom, const ceed::CeedGeomFactorData &data,
      ceed::CeedElementSubset subset = ceed::CeedElementSubset::Full) const;

  // If the space has a special element restriction for discrete interpolators, return that.
  // Otherwise return the same restriction as given by GetCeedElemRestriction.
  CeedElemRestriction
  GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                               const ceed::CeedGeomFactorData &data) const;

  // If the space has a special element restriction for the range space of discrete
  // interpolators, return that. Otherwise return the same restriction as given by
  // GetCeedElemRestriction.
  CeedElemRestriction
  GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                    const ceed::CeedGeomFactorData &data) const;

  // Clear the cached basis and element restriction objects owned by the finite element
  // space.
  void ResetCeedObjects();

  void Update()
  {
    MFEM_VERIFY(!broken, "Broken finite element spaces cannot be updated!");
    ResetCeedObjects();
  }

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
    if constexpr (std::is_same_v<VecType, ComplexVector>)
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
    if constexpr (std::is_same_v<VecType, ComplexVector>)
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
    if constexpr (std::is_same_v<VecType, ComplexVector>)
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
