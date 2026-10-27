// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include <algorithm>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/restriction.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

FiniteElementSpace::FiniteElementSpace(Mesh &mesh, const mfem::FiniteElementCollection *fec,
                                       int vdim, int ordering)
  : fespace(&mesh.Get(), fec, vdim, ordering), mesh(mesh)
{
  ResetCeedObjects();
  tx.UseDevice(true);
  lx.UseDevice(true);
  ly.UseDevice(true);
}

FiniteElementSpace::FiniteElementSpace(Mesh &mesh, const mfem::FiniteElementCollection *fec,
                                       const std::vector<int> &attributes)
  : FiniteElementSpace(mesh, fec)
{
  break_attributes = attributes;
  std::sort(break_attributes.begin(), break_attributes.end());
  break_attributes.erase(std::unique(break_attributes.begin(), break_attributes.end()),
                         break_attributes.end());
  InitBroken();
}

void FiniteElementSpace::InitBroken()
{
  if (std::any_of(break_attributes.begin(), break_attributes.end(),
                  [this](int attr) { return mesh.Get().bdr_attributes.Find(attr) >= 0; }))
  {
    broken = fem::BuildBrokenDofLayout(fespace, break_attributes);
    // Complete construction collectively, before threaded assembly or applications.
    broken->R = fem::BuildBrokenRestriction(*this);
  }
}

void FiniteElementSpace::Update()
{
  ResetCeedObjects();
  G.reset();
  aux_fespace = nullptr;
  broken.reset();
  fespace.Update(false);
  InitBroken();
}

void FiniteElementSpace::GetElementDofs(int e, mfem::Array<int> &dofs,
                                        mfem::DofTransformation &dt) const
{
  fespace.GetElementDofs(e, dofs, dt);
  if (broken)
  {
    const int start = broken->offsets[e];
    for (int j = 0; j < dofs.Size(); j++)
    {
      dofs[j] = broken->dofs[start + j];
    }
  }
}

void FiniteElementSpace::RemapDofs(int e, mfem::Array<int> &dofs) const
{
  if (!broken)
  {
    return;
  }
  mfem::Array<int> original;
  mfem::DofTransformation dt;
  fespace.GetElementDofs(e, original, dt);
  for (auto &d : dofs)
  {
    const int l = d >= 0 ? d : -1 - d;
    int j = 0;
    while (j < original.Size() && (original[j] >= 0 ? original[j] : -1 - original[j]) != l)
    {
      j++;
    }
    MFEM_VERIFY(j < original.Size(), "Trace DOF missing from its incident element!");
    const int b = broken->dofs[broken->offsets[e] + j];
    const int lb = b >= 0 ? b : -1 - b;
    d = d >= 0 ? lb : -1 - lb;
  }
}

void FiniteElementSpace::GetBdrElementDofs(int be, mfem::Array<int> &dofs,
                                           mfem::DofTransformation &dt) const
{
  fespace.GetBdrElementDofs(be, dofs, dt);
  int e, info;
  mesh.Get().GetBdrElementAdjacentElement(be, e, info);
  RemapDofs(e, dofs);
}

void FiniteElementSpace::GetElementEntityDofs(int e, int dimension, int local_entity,
                                              mfem::Array<int> &dofs) const
{
  mfem::Array<int> entities, orientation;
  if (dimension == 1)
  {
    mesh.Get().GetElementEdges(e, entities, orientation);
    fespace.GetEdgeDofs(entities[local_entity], dofs);
  }
  else
  {
    mesh.Get().GetElementFaces(e, entities, orientation);
    fespace.GetFaceDofs(entities[local_entity], dofs);
  }
  RemapDofs(e, dofs);
}

void FiniteElementSpace::GetInteriorFaceDofs(int f, mfem::Array<int> &first,
                                             mfem::DofTransformation &dt1,
                                             mfem::Array<int> &second,
                                             mfem::DofTransformation &dt2) const
{
  const auto face = mesh.Get().GetFaceInformation(f);
  MFEM_VERIFY(face.IsLocal() || face.IsShared(), "Expected an interior face patch!");
  GetElementDofs(face.element[0].index, first, dt1);
  if (face.IsLocal())
  {
    GetElementDofs(face.element[1].index, second, dt2);
  }
  else
  {
    const int e = face.element[1].index;
    fespace.GetFaceNbrElementVDofs(e, second, dt2);
    MFEM_VERIFY(broken, "Neighbour trace rows require an owned broken space!");
    for (int j = 0; j < second.Size(); j++)
    {
      second[j] = broken->nbr_dofs[broken->nbr_offsets[e] + j];
    }
  }
}

std::vector<int>
FiniteElementSpace::GetInteriorBoundaryFaces(const std::vector<int> &attributes) const
{
  auto &pmesh = const_cast<mfem::ParMesh &>(mesh.Get());
  const auto comm = GetComm();
  const int rank = Mpi::Rank(comm);
  std::vector<int> selected(attributes);
  std::sort(selected.begin(), selected.end());
  mfem::L2_FECollection p0(0, Dimension());
  mfem::ParFiniteElementSpace header_space(&pmesh, &p0);
  mfem::ParGridFunction header(&header_space);
  std::vector<int> marked(pmesh.GetNE(), 0);
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (std::binary_search(selected.begin(), selected.end(), pmesh.GetBdrAttribute(be)))
    {
      int e, info;
      pmesh.GetBdrElementAdjacentElement(be, e, info);
      marked[e] |= 1 << pmesh.DecodeFaceInfoLocalIndex(info);
    }
  }
  mfem::Array<int> dofs;
  mfem::DofTransformation dt;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    header_space.GetElementDofs(e, dofs);
    header(dofs[0]) = marked[e];
  }
  header.ExchangeFaceNbrData();
  std::vector<int> neighbour(pmesh.GetNFaceNeighborElements());
  for (int e = 0; e < static_cast<int>(neighbour.size()); e++)
  {
    header_space.GetFaceNbrElementVDofs(e, dofs, dt);
    neighbour[e] = static_cast<int>(header.FaceNbrData()(dofs[0]));
  }
  std::vector<int> result;
  for (int f = 0; f < pmesh.GetNumFacesWithGhost(); f++)
  {
    const auto face = pmesh.GetFaceInformation(f);
    if (!face.IsLocal() && !face.IsShared())
    {
      continue;
    }
    const int e = face.element[0].index, other = face.element[1].index;
    const bool here = (marked[e] >> face.element[0].local_face_id) & 1;
    const bool there = ((face.IsShared() ? neighbour[other] : marked[other]) >>
                        face.element[1].local_face_id) &
                       1;
    if (face.IsLocal())
    {
      if (here || there)
      {
        result.push_back(f);
      }
    }
    else if (here)
    {
      const auto *offset = pmesh.face_nbr_elements_offset.GetData();
      const int fn =
          std::upper_bound(offset, offset + pmesh.GetNFaceNeighbors() + 1, other) - offset -
          1;
      if (!there || rank < pmesh.GetFaceNbrRank(fn))
      {
        result.push_back(f);
      }
    }
  }
  return result;
}

void FiniteElementSpace::GetEssentialTrueDofs(const mfem::Array<int> &marker,
                                              mfem::Array<int> &dofs) const
{
  if (!broken)
  {
    fespace.GetEssentialTrueDofs(marker, dofs);
    return;
  }
  Vector local(GetVSize()), owned(GetTrueVSize());
  local = 0.0;
  mfem::Array<int> boundary;
  mfem::DofTransformation dt;
  for (int be = 0; be < mesh.Get().GetNBE(); be++)
  {
    const int attr = mesh.Get().GetBdrAttribute(be);
    if (attr <= marker.Size() && marker[attr - 1])
    {
      GetBdrElementDofs(be, boundary, dt);
      for (auto d : boundary)
      {
        local(d >= 0 ? d : -1 - d) = 1.0;
      }
    }
  }
  broken->P->AbsMultTranspose(1.0, local, 0.0, owned);
  dofs.SetSize(0);
  for (int t = 0; t < owned.Size(); t++)
  {
    if (owned(t) > 0.0)
    {
      dofs.Append(t);
    }
  }
}

void FiniteElementSpace::CopyFromConforming(const Vector &x, Vector &y) const
{
  MFEM_VERIFY(x.Size() == fespace.GetTrueVSize(), "Invalid conforming field size!");
  if (!broken)
  {
    y = x;
    return;
  }
  Vector original(fespace.GetVSize()), local(GetVSize());
  fespace.GetProlongationMatrix()->Mult(x, original);
  local = 0.0;
  for (int l = 0; l < local.Size(); l++)
  {
    const int source = broken->row_origin[l];
    if (source >= 0)
    {
      local(l) = original(source);
    }
  }
  y.SetSize(GetTrueVSize());
  broken->R->Mult(local, y);
}

void FiniteElementSpace::GetTrueDofCoordinates(int component, Vector &coordinates) const
{
  MFEM_VERIFY(GetMaxElementOrder() == 1 && GetVDim() == 1 &&
                  GetFEColl().GetMapType(Dimension()) == mfem::FiniteElement::VALUE &&
                  component >= 0 && component < SpaceDimension(),
              "True DOF coordinates require order-1 scalar H1!");
  mfem::ParGridFunction local(const_cast<mfem::ParFiniteElementSpace *>(&fespace));
  MFEM_VERIFY(local.Size() == mesh.Get().GetNV(), "Unexpected H1 coordinate layout!");
  for (int v = 0; v < local.Size(); v++)
  {
    local(v) = mesh.Get().GetVertex(v)[component];
  }
  Vector true_coords(fespace.GetTrueVSize());
  local.ParallelProject(true_coords);
  CopyFromConforming(true_coords, coordinates);
}

CeedBasis FiniteElementSpace::GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  auto &map = basis.at(ceed);
  auto it = map.find(geom);
  return it != map.end()
             ? it->second
             : map.emplace(geom, BuildCeedBasis(fespace, ceed, geom)).first->second;
}

CeedElemRestriction
FiniteElementSpace::BuildSpaceRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                          const std::vector<int> &indices, bool interp,
                                          bool range) const
{
  if (!broken)
  {
    return BuildCeedElemRestriction(fespace, ceed, geom, indices, interp, range);
  }
  ceed::ElementDofRemap remap;
  remap.l_size = GetVSize();
  remap.apply = [this](int index, bool boundary, mfem::Array<int> &dofs)
  {
    if (boundary)
    {
      int e, info;
      mesh.Get().GetBdrElementAdjacentElement(index, e, info);
      RemapDofs(e, dofs);
    }
    else
    {
      for (int j = 0; j < dofs.Size(); j++)
      {
        dofs[j] = broken->dofs[broken->offsets[index] + j];
      }
    }
  };
  CeedElemRestriction result;
  ceed::InitRestriction(fespace, indices, mfem::Geometry::Dimension[geom] != Dimension(),
                        interp, range, ceed, &result, &remap);
  return result;
}

CeedElemRestriction
FiniteElementSpace::GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                           const ceed::CeedGeomFactorData &data,
                                           ceed::CeedElementSubset subset) const
{
  const bool active = subset == ceed::CeedElementSubset::Active &&
                      data.active_indices.size() < data.indices.size();
  auto &map = (active ? active_restr : restr).at(ceed);
  auto it = map.find(geom);
  if (it != map.end())
  {
    return it->second;
  }
  std::vector<int> selected;
  if (active)
  {
    MFEM_VERIFY(!data.active_indices.empty(), "Empty active restriction subset!");
    for (auto i : data.active_indices)
    {
      selected.push_back(data.indices[i]);
    }
  }
  return map
      .emplace(geom, BuildSpaceRestriction(ceed, geom, active ? selected : data.indices))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                 const ceed::CeedGeomFactorData &data) const
{
  if (!HasUniqueInterpRestriction(*GetFEColl().FiniteElementForGeometry(geom)))
  {
    return GetCeedElemRestriction(ceed, geom, data);
  }
  auto &map = interp_restr.at(ceed);
  auto it = map.find(geom);
  return it != map.end()
             ? it->second
             : map.emplace(geom, BuildSpaceRestriction(ceed, geom, data.indices, true))
                   .first->second;
}

CeedElemRestriction FiniteElementSpace::GetInterpRangeCeedElemRestriction(
    Ceed ceed, mfem::Geometry::Type geom, const ceed::CeedGeomFactorData &data) const
{
  if (!HasUniqueInterpRangeRestriction(*GetFEColl().FiniteElementForGeometry(geom)))
  {
    return GetInterpCeedElemRestriction(ceed, geom, data);
  }
  auto &map = interp_range_restr.at(ceed);
  auto it = map.find(geom);
  return it != map.end() ? it->second
                         : map.emplace(geom, BuildSpaceRestriction(ceed, geom, data.indices,
                                                                   true, true))
                               .first->second;
}

void FiniteElementSpace::ResetCeedObjects()
{
  for (auto &[ceed, map] : basis)
  {
    for (auto &[geom, val] : map)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&val));
    }
  }
  auto Destroy = [](auto &objects)
  {
    for (auto &[ceed, map] : objects)
    {
      for (auto &[geom, val] : map)
      {
        PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
      }
    }
    objects.clear();
  };
  Destroy(restr);
  Destroy(active_restr);
  Destroy(interp_restr);
  Destroy(interp_range_restr);
  basis.clear();
  for (auto ceed : ceed::internal::GetCeedObjects())
  {
    basis.emplace(ceed, ceed::GeometryObjectMap<CeedBasis>());
    restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    active_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_range_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
  }
}

CeedBasis FiniteElementSpace::BuildCeedBasis(const mfem::FiniteElementSpace &space,
                                             Ceed ceed, mfem::Geometry::Type geom)
{
  mfem::IsoparametricTransformation T;
  const auto *nodal =
      space.GetMesh()->GetNodalFESpace()->FEColl()->FiniteElementForGeometry(geom);
  if (!nodal)
  {
    nodal =
        space.GetMesh()->GetNodalFESpace()->FEColl()->TraceFiniteElementForGeometry(geom);
  }
  T.SetFE(nodal);
  const auto &ir = mfem::IntRules.Get(geom, fem::DefaultIntegrationOrder::Get(T));
  const auto *fe = space.FEColl()->FiniteElementForGeometry(geom);
  if (!fe)
  {
    fe = space.FEColl()->TraceFiniteElementForGeometry(geom);
  }
  CeedBasis result;
  ceed::InitBasis(*fe, ir, space.GetVDim(), ceed, &result);
  return result;
}

CeedElemRestriction FiniteElementSpace::BuildCeedElemRestriction(
    const mfem::FiniteElementSpace &space, Ceed ceed, mfem::Geometry::Type geom,
    const std::vector<int> &indices, bool interp, bool range)
{
  CeedElemRestriction result;
  ceed::InitRestriction(space, indices,
                        mfem::Geometry::Dimension[geom] != space.GetMesh()->Dimension(),
                        interp, range, ceed, &result);
  return result;
}

std::unique_ptr<mfem::HypreParMatrix> FiniteElementSpace::AssembleInterpolator(
    const FiniteElementSpace &trial, const FiniteElementSpace &test,
    mfem::DiscreteInterpolator *interpolator, bool refinement)
{
  MFEM_VERIFY(trial.GetBreakAttributes() == test.GetBreakAttributes() &&
                  trial.GetVDim() == 1 && test.GetVDim() == 1,
              "Interpolators require scalar layouts with the same break list!");
  auto &mesh = const_cast<mfem::ParMesh &>(test.GetParMesh());
  const mfem::CoarseFineTransformations *cf =
      refinement ? &mesh.GetRefinementTransforms() : nullptr;
  if (refinement)
  {
    MFEM_VERIFY(
        mesh.GetLastOperation() == mfem::Mesh::REFINE &&
            cf->embeddings.Size() == mesh.GetNE() &&
            &trial.GetFEColl() == &test.GetFEColl(),
        "H-transfer requires the same collection and parent refinement embeddings!");
  }
  auto local = std::make_unique<mfem::SparseMatrix>(test.GetVSize(), trial.GetVSize());
  Vector count(test.GetVSize());
  count = 0.0;
  mfem::Array<int> in, out;
  mfem::DofTransformation idt, odt;
  mfem::DenseMatrix native, matrix;
  mfem::IsoparametricTransformation embedding;
  Vector v, w;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int parent = refinement ? cf->embeddings[e].parent : e;
    trial.GetElementDofs(parent, in, idt);
    test.GetElementDofs(e, out, odt);
    if (refinement)
    {
      // Palace-assembled h-transfer: MFEM's assembled RefinementMatrix misses
      // ND tetrahedral DOF transformations. Retain both primal transformations.
      const auto geom = mesh.GetElementGeometry(e);
      embedding.SetIdentityTransformation(geom);
      embedding.SetPointMat(cf->point_matrices[geom](cf->embeddings[e].matrix));
      test.GetFE(e)->GetLocalInterpolation(embedding, native);
    }
    else
    {
      interpolator->AssembleElementMatrix2(*trial.GetFE(e), *test.GetFE(e),
                                           *mesh.GetElementTransformation(e), native);
    }
    matrix.SetSize(out.Size(), in.Size());
    v.SetSize(in.Size());
    w.SetSize(out.Size());
    for (int j = 0; j < in.Size(); j++)
    {
      v = 0.0;
      v(j) = 1.0;
      if (!idt.IsIdentity())
      {
        idt.InvTransformPrimal(v);
      }
      native.Mult(v, w);
      if (!odt.IsIdentity())
      {
        odt.TransformPrimal(w);
      }
      matrix.SetCol(j, w);
    }
    local->AddSubMatrix(out, in, matrix, 0);
    for (auto d : out)
    {
      count(d >= 0 ? d : -1 - d) += 1.0;
    }
  }
  local->Finalize();
  // Interpolation writes a coefficient once. This normalization is in I, never R.
  for (int l = 0; l < local->Height(); l++)
  {
    for (int k = local->GetI()[l]; k < local->GetI()[l + 1]; k++)
    {
      local->GetData()[k] /= count(l);
    }
  }
  ParOperator result(std::move(local), trial, test, true);
  return result.StealParallelAssemble();
}

const Operator &FiniteElementSpace::BuildDiscreteInterpolator() const
{
  const auto &trial = *aux_fespace;
  MFEM_VERIFY(trial.GetBreakAttributes() == GetBreakAttributes(),
              "Discrete interpolators require the same break attributes!");
  const int dim = Dimension();
  const auto in = trial.GetFEColl().GetMapType(dim), out = GetFEColl().GetMapType(dim);
  if (!broken && !trial.broken && in == mfem::FiniteElement::H_CURL &&
      out == mfem::FiniteElement::INTEGRAL)
  {
    // Preserve today's conforming scalar-curl assembly and its write convention.
    mfem::DiscreteLinearOperator interp(
        const_cast<mfem::ParFiniteElementSpace *>(&trial.fespace),
        const_cast<mfem::ParFiniteElementSpace *>(&fespace));
    interp.AddDomainInterpolator(new mfem::CurlInterpolator);
    interp.Assemble();
    interp.Finalize();
    G = std::make_unique<ParOperator>(std::unique_ptr<mfem::SparseMatrix>(interp.LoseMat()),
                                      trial, *this, true);
  }
  else if (broken || trial.broken)
  {
    std::unique_ptr<mfem::DiscreteInterpolator> interp;
    if (in == mfem::FiniteElement::VALUE && out == mfem::FiniteElement::H_CURL)
    {
      interp = std::make_unique<mfem::GradientInterpolator>();
    }
    else if (in == mfem::FiniteElement::H_CURL &&
             (out == mfem::FiniteElement::H_DIV || out == mfem::FiniteElement::INTEGRAL))
    {
      interp = std::make_unique<mfem::CurlInterpolator>();
    }
    else if (in == mfem::FiniteElement::H_DIV && out == mfem::FiniteElement::INTEGRAL)
    {
      interp = std::make_unique<mfem::DivergenceInterpolator>();
    }
    MFEM_VERIFY(interp, "Unsupported discrete interpolator spaces!");
    G = AssembleInterpolator(trial, *this, interp.get());
  }
  else
  {
    DiscreteLinearOperator interp(trial, *this);
    if (in == mfem::FiniteElement::VALUE && out == mfem::FiniteElement::H_CURL)
    {
      interp.AddDomainInterpolator<GradientInterpolator>();
    }
    else if (in == mfem::FiniteElement::H_CURL && out == mfem::FiniteElement::H_DIV)
    {
      interp.AddDomainInterpolator<CurlInterpolator>();
    }
    else if (in == mfem::FiniteElement::H_DIV && out == mfem::FiniteElement::INTEGRAL)
    {
      interp.AddDomainInterpolator<DivergenceInterpolator>();
    }
    else
    {
      MFEM_ABORT("Unsupported discrete interpolator spaces!");
    }
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial, *this, true);
  }
  return *G;
}

const Operator &FiniteElementSpaceHierarchy::BuildProlongationAtLevel(std::size_t l) const
{
  const auto &coarse = *fespaces[l], &fine = *fespaces[l + 1];
  const bool refinement = &coarse.GetMesh() != &fine.GetMesh();
  if (coarse.IsBroken() || fine.IsBroken())
  {
    mfem::IdentityInterpolator interp;
    P[l] = FiniteElementSpace::AssembleInterpolator(coarse, fine, &interp, refinement);
  }
  else if (refinement)
  {
    P[l] = std::make_unique<ParOperator>(
        std::make_unique<mfem::TransferOperator>(coarse, fine), coarse, fine, true);
  }
  else
  {
    DiscreteLinearOperator interp(coarse, fine);
    interp.AddDomainInterpolator<IdentityInterpolator>();
    P[l] = std::make_unique<ParOperator>(interp.PartialAssemble(), coarse, fine, true);
  }
  return *P[l];
}

}  // namespace palace
