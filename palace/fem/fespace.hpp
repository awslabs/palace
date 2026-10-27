// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <initializer_list>
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

// One owner for the FE metadata, actual element/trace maps, and parallel constraints.
// An MFEM space is deliberately not an access path to a broken field.
class FiniteElementSpace
{
private:
  friend class FiniteElementSpaceHierarchy;
  mfem::ParFiniteElementSpace fespace;
  Mesh &mesh;
  std::vector<int> break_attributes;
  std::unique_ptr<fem::BrokenDofLayout> broken;
  mutable ceed::CeedObjectMap<CeedBasis> basis;
  mutable ceed::CeedObjectMap<CeedElemRestriction> restr, active_restr, interp_restr,
      interp_range_restr;
  mutable ComplexVector tx, lx, ly;
  mutable const FiniteElementSpace *aux_fespace = nullptr;
  mutable std::unique_ptr<Operator> G;

  void InitBroken();
  void RemapDofs(int e, mfem::Array<int> &dofs) const;
  CeedElemRestriction BuildSpaceRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                            const std::vector<int> &indices,
                                            bool interp = false, bool range = false) const;
  const Operator &BuildDiscreteInterpolator() const;
  static std::unique_ptr<mfem::HypreParMatrix>
  AssembleInterpolator(const FiniteElementSpace &trial, const FiniteElementSpace &test,
                       mfem::DiscreteInterpolator *interpolator, bool refinement = false);
  bool HasUniqueInterpRestriction(const mfem::FiniteElement &fe) const
  {
    const auto *tfe = dynamic_cast<const mfem::TensorBasisElement *>(&fe);
    return tfe && tfe->GetDofMap().Size() &&
           fe.GetRangeType() != mfem::FiniteElement::VECTOR;
  }
  bool HasUniqueInterpRangeRestriction(const mfem::FiniteElement &fe) const
  {
    const auto *dt = fespace.FEColl()->DofTransformationForGeometry(fe.GetGeomType());
    return mesh.Dimension() == 3 && dt && !dt->IsIdentity();
  }

public:
  FiniteElementSpace(Mesh &mesh, const mfem::FiniteElementCollection *fec, int vdim = 1,
                     int ordering = mfem::Ordering::byNODES);
  // Missing attributes are ignored. Empty list is the ordinary conforming path.
  // Nonempty lists require uniform order and vdim=1. Construction is collective.
  FiniteElementSpace(Mesh &mesh, const mfem::FiniteElementCollection *fec,
                     const std::vector<int> &break_attributes);
  FiniteElementSpace(Mesh &mesh, const mfem::FiniteElementCollection *fec,
                     std::initializer_list<int> attributes)
    : FiniteElementSpace(mesh, fec, std::vector<int>(attributes))
  {
  }
  virtual ~FiniteElementSpace() { ResetCeedObjects(); }

  const auto &Get() const
  {
    MFEM_VERIFY(!broken, "A broken space is not an MFEM space; use its DOF API!");
    return fespace;
  }
  auto &Get()
  {
    MFEM_VERIFY(!broken, "A broken space is not an MFEM space; use its DOF API!");
    return fespace;
  }
  operator const mfem::ParFiniteElementSpace &() const { return Get(); }
  operator mfem::ParFiniteElementSpace &() { return Get(); }

  const auto &GetFEColl() const { return *fespace.FEColl(); }
  auto &GetFEColl() { return *fespace.FEColl(); }
  const auto *GetFE(int e) const { return fespace.GetFE(e); }
  const auto *GetBE(int be) const { return fespace.GetBE(be); }
  const auto &GetBreakAttributes() const { return break_attributes; }
  bool IsBroken() const { return broken != nullptr; }
  const auto &GetMesh() const { return mesh; }
  auto &GetMesh() { return mesh; }
  const auto &GetParMesh() const { return mesh.Get(); }
  auto &GetParMesh() { return mesh.Get(); }
  int Dimension() const { return mesh.Dimension(); }
  int SpaceDimension() const { return mesh.SpaceDimension(); }
  int GetVDim() const { return fespace.GetVDim(); }
  int GetVectorDim() const { return fespace.GetVectorDim(); }
  int GetMaxElementOrder() const { return fespace.GetMaxElementOrder(); }
  int GetVSize() const { return broken ? broken->P->Height() : fespace.GetVSize(); }
  int GetTrueVSize() const { return broken ? broken->P->Width() : fespace.GetTrueVSize(); }
  auto GlobalVSize() const
  {
    return broken ? broken->P->GetGlobalNumRows() : fespace.GlobalVSize();
  }
  auto GlobalTrueVSize() const
  {
    return broken ? broken->P->GetGlobalNumCols() : fespace.GlobalTrueVSize();
  }
  MPI_Comm GetComm() const { return fespace.GetComm(); }
  const Operator *GetProlongationMatrix() const
  {
    return broken ? broken->P.get() : fespace.GetProlongationMatrix();
  }
  const mfem::HypreParMatrix *Dof_TrueDof_Matrix() const
  {
    return broken ? broken->P.get() : fespace.Dof_TrueDof_Matrix();
  }
  HYPRE_BigInt *GetDofOffsets() const
  {
    return broken ? broken->row_offsets.GetData() : fespace.GetDofOffsets();
  }
  HYPRE_BigInt *GetTrueDofOffsets() const
  {
    return broken ? broken->col_offsets.GetData() : fespace.GetTrueDofOffsets();
  }
  const mfem::SparseMatrix *GetRestrictionMatrix() const
  {
    return broken ? broken->R.get() : fespace.GetRestrictionMatrix();
  }

  void GetElementDofs(int e, mfem::Array<int> &dofs, mfem::DofTransformation &dt) const;
  void GetElementVDofs(int e, mfem::Array<int> &dofs, mfem::DofTransformation &dt) const
  {
    if (broken)
    {
      GetElementDofs(e, dofs, dt);
    }
    else
    {
      fespace.GetElementVDofs(e, dofs, dt);
    }
  }
  void GetBdrElementDofs(int be, mfem::Array<int> &dofs, mfem::DofTransformation &dt) const;
  // DOFs in an entity's closure, in the copies read by this incident element.
  void GetElementEntityDofs(int e, int dimension, int local_entity,
                            mfem::Array<int> &dofs) const;
  // Element 0 is local; element 1 may be a face neighbour. Both return volume DOFs
  // in their native chart; the face transformations supply the common/subface chart.
  void GetInteriorFaceDofs(int face, mfem::Array<int> &first, mfem::DofTransformation &dt1,
                           mfem::Array<int> &second, mfem::DofTransformation &dt2) const;
  const mfem::FiniteElement *GetFaceNbrFE(int e) const { return fespace.GetFaceNbrFE(e); }
  void GetFaceNbrElementDofs(int e, mfem::Array<int> &dofs,
                             mfem::DofTransformation &dt) const;
  // Collective: each selected physical interior patch is returned on exactly one
  // rank holding its boundary representation, including master/slave and shared faces.
  std::vector<int> GetInteriorBoundaryFaces(const std::vector<int> &attributes) const;
  void GetEssentialTrueDofs(const mfem::Array<int> &marker, mfem::Array<int> &dofs) const;
  void CopyFromConforming(const Vector &x, Vector &y) const;
  void GetTrueDofCoordinates(int component, Vector &coordinates) const;

  const auto &GetDiscreteInterpolator(const FiniteElementSpace &aux) const
  {
    if (&aux != aux_fespace)
    {
      G.reset();
      aux_fespace = &aux;
    }
    return G ? *G : BuildDiscreteInterpolator();
  }
  CeedBasis GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const;
  CeedElemRestriction GetCeedElemRestriction(
      Ceed ceed, mfem::Geometry::Type geom, const ceed::CeedGeomFactorData &data,
      ceed::CeedElementSubset subset = ceed::CeedElementSubset::Full) const;
  CeedElemRestriction
  GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                               const ceed::CeedGeomFactorData &data) const;
  CeedElemRestriction
  GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                    const ceed::CeedGeomFactorData &data) const;
  void ResetCeedObjects();
  // Recompute the representation after a mesh update. Old vectors/operators are invalid.
  void Update();
  static CeedBasis BuildCeedBasis(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                                  mfem::Geometry::Type geom);
  static CeedElemRestriction
  BuildCeedElemRestriction(const mfem::FiniteElementSpace &fespace, Ceed ceed,
                           mfem::Geometry::Type geom, const std::vector<int> &indices,
                           bool interp = false, bool range = false);

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
};

class FiniteElementSpaceHierarchy
{
protected:
  std::vector<std::unique_ptr<FiniteElementSpace>> fespaces;
  mutable std::vector<std::unique_ptr<Operator>> P;
  const Operator &BuildProlongationAtLevel(std::size_t l) const;

public:
  FiniteElementSpaceHierarchy() = default;
  FiniteElementSpaceHierarchy(std::unique_ptr<FiniteElementSpace> &&space)
  {
    AddLevel(std::move(space));
  }
  auto GetNumLevels() const { return fespaces.size(); }
  void AddLevel(std::unique_ptr<FiniteElementSpace> &&space)
  {
    MFEM_VERIFY(fespaces.empty() ||
                    fespaces.back()->GetBreakAttributes() == space->GetBreakAttributes(),
                "Hierarchy levels must have the same break attributes!");
    fespaces.push_back(std::move(space));
    P.push_back(nullptr);
  }
  auto &GetFESpaceAtLevel(std::size_t l)
  {
    MFEM_ASSERT(l < GetNumLevels(), "Invalid finite element space level!");
    return *fespaces[l];
  }
  const auto &GetFESpaceAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l < GetNumLevels(), "Invalid finite element space level!");
    return *fespaces[l];
  }
  auto &GetFinestFESpace() { return GetFESpaceAtLevel(GetNumLevels() - 1); }
  const auto &GetFinestFESpace() const { return GetFESpaceAtLevel(GetNumLevels() - 1); }
  const auto &GetProlongationAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l + 1 < GetNumLevels(), "Invalid prolongation level!");
    return P[l] ? *P[l] : BuildProlongationAtLevel(l);
  }
  std::vector<const Operator *> GetProlongationOperators() const
  {
    std::vector<const Operator *> result(GetNumLevels() - 1);
    for (std::size_t l = 0; l < result.size(); l++)
    {
      result[l] = &GetProlongationAtLevel(l);
    }
    return result;
  }
  const auto &GetDiscreteInterpolatorAtLevel(std::size_t l,
                                             const FiniteElementSpace &aux) const
  {
    return GetFESpaceAtLevel(l).GetDiscreteInterpolator(aux);
  }
  std::vector<const Operator *>
  GetDiscreteInterpolators(const FiniteElementSpaceHierarchy &aux) const
  {
    std::vector<const Operator *> result(GetNumLevels(), nullptr);
    for (std::size_t l = 1; l < result.size(); l++)
    {
      result[l] = &GetDiscreteInterpolatorAtLevel(l, aux.GetFESpaceAtLevel(l));
    }
    return result;
  }
};

}  // namespace palace
#endif  // PALACE_FEM_FESPACE_HPP
