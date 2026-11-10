// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <typeinfo>
#include <utility>
#include <ceed/backend.h>
#include <mfem.hpp>
#include <mfem/general/forall.hpp>
#include "fem/fespace.hpp"
#include "fem/qfunctions/apply/complex_apply_qf.h"
#include "linalg/hypre.hpp"
#include "utils/omp.hpp"

namespace palace::ceed
{

Operator::Operator(int h, int w) : palace::Operator(h, w)
{
  const std::size_t nt = internal::GetCeedObjects().size();
  op.resize(nt, nullptr);
  op_t.resize(nt, nullptr);
  u.resize(nt, nullptr);
  v.resize(nt, nullptr);
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    CeedOperator loc_op, loc_op_t;
    CeedVector loc_u, loc_v;
    PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &loc_op));
    PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &loc_op_t));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, width, &loc_u));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, height, &loc_v));
    op[id] = loc_op;
    op_t[id] = loc_op_t;
    u[id] = loc_u;
    v[id] = loc_v;
  }
  temp.UseDevice(true);
}

void Operator::AddSubOperator(CeedOperator sub_op, CeedOperator sub_op_t)
{
  // This should be called from within a OpenMP parallel region.
  const int id = utils::GetThreadNum();
  MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
              "Out of bounds access for thread number " << id << "!");
  Ceed ceed;
  PalaceCeedCallBackend(CeedOperatorGetCeed(sub_op, &ceed));
  CeedSize l_in, l_out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(sub_op, &l_in, &l_out));
  MFEM_VERIFY((l_in < 0 || mfem::internal::to_int(l_in) == width) &&
                  (l_out < 0 || mfem::internal::to_int(l_out) == height),
              "Dimensions mismatch for CeedOperator!");
  PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(op[id], sub_op));
  PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
  if (sub_op_t)
  {
    Ceed ceed_t;
    PalaceCeedCallBackend(CeedOperatorGetCeed(sub_op_t, &ceed_t));
    MFEM_VERIFY(ceed_t == ceed, "Ceed context mismatch for transpose CeedOperator!");
    CeedSize l_in_t, l_out_t;
    PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(sub_op_t, &l_in_t, &l_out_t));
    MFEM_VERIFY(l_in_t == l_out && l_out_t == l_in,
                "Dimensions mismatch for transpose CeedOperator!");
    PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(op_t[id], sub_op_t));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
  }
}

void Operator::Finalize()
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(op[id]));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(op_t[id]));
  }
}

void Operator::DestroyAssemblyData() const
{
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));
  }
}

void Operator::AssembleDiagonal(Vector &diag) const
{
  Ceed ceed;
  MFEM_VERIFY(diag.Size() == height, "Invalid size for diagonal vector!");
  diag = 0.0;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  const CeedMemType mem = GetUsableMemType(ceed);
  auto *diag_data = diag.ReadWrite(mem == CEED_MEM_DEVICE);

  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(v[id], mem, CEED_USE_POINTER, diag_data));
    PalaceCeedCall(
        ceed, CeedOperatorLinearAssembleAddDiagonal(op[id], v[id], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[id], mem, nullptr));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));
  }
}

namespace
{

inline void CeedAddMult(const std::vector<CeedOperator> &op,
                        const std::vector<CeedVector> &u, const std::vector<CeedVector> &v,
                        const Vector &x, Vector &y)
{
  Ceed ceed;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  const CeedMemType mem = GetUsableMemType(ceed);
  const auto *x_data = x.Read(mem == CEED_MEM_DEVICE);
  auto *y_data = y.ReadWrite(mem == CEED_MEM_DEVICE);

  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(u[id], mem, CEED_USE_POINTER,
                                            const_cast<CeedScalar *>(x_data)));
    PalaceCeedCall(ceed, CeedVectorSetArray(v[id], mem, CEED_USE_POINTER, y_data));
    PalaceCeedCall(ceed,
                   CeedOperatorApplyAdd(op[id], u[id], v[id], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(u[id], mem, nullptr));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[id], mem, nullptr));
  }
}

}  // namespace

void Operator::Mult(const Vector &x, Vector &y) const
{
  if (ea)
  {
    ElementAssemblyMult(x, y, false);
    return;
  }
  y = 0.0;
  CeedAddMult(op, u, v, x, y);
  if (dof_multiplicity.Size() > 0)
  {
    y *= dof_multiplicity;
  }
}

void Operator::AddMult(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0, "ceed::Operator::AddMult only supports coefficient = 1.0!");
  if (ea)
  {
    ElementAssemblyMult(x, y, true);
    return;
  }
  if (dof_multiplicity.Size() > 0)
  {
    temp.SetSize(height);
    temp = 0.0;
    CeedAddMult(op, u, v, x, temp);
    {
      const auto *d_dof_multiplicity = dof_multiplicity.Read();
      const auto *d_temp = temp.Read();
      auto *d_y = y.ReadWrite();
      mfem::forall(height, [=] MFEM_HOST_DEVICE(int i)
                   { d_y[i] += d_dof_multiplicity[i] * d_temp[i]; });
    }
  }
  else
  {
    CeedAddMult(op, u, v, x, y);
  }
}

void Operator::MultTranspose(const Vector &x, Vector &y) const
{
  y = 0.0;
  AddMultTranspose(x, y);
}

void Operator::AddMultTranspose(const Vector &x, Vector &y, const double a) const
{
  MFEM_VERIFY(a == 1.0,
              "ceed::Operator::AddMultTranspose only supports coefficient = 1.0!");
  if (dof_multiplicity.Size() > 0)
  {
    temp.SetSize(height);
    {
      const auto *d_dof_multiplicity = dof_multiplicity.Read();
      const auto *d_x = x.Read();
      auto *d_temp = temp.Write();
      mfem::forall(height, [=] MFEM_HOST_DEVICE(int i)
                   { d_temp[i] = d_dof_multiplicity[i] * d_x[i]; });
    }
    CeedAddMult(op_t, v, u, temp, y);
  }
  else
  {
    CeedAddMult(op_t, v, u, x, y);
  }
}

namespace
{

// Scope the references acquired while inspecting or constructing CEED objects across
// every rejection path. Fields and composites retain their own references when objects
// are installed, so releasing these temporary handles leaves the originals intact.
template <typename T, int (*Destroy)(T *)>
struct ScopedCeed
{
  T value = nullptr;
  ScopedCeed() = default;
  ScopedCeed(const ScopedCeed &) = delete;
  ScopedCeed &operator=(const ScopedCeed &) = delete;
  ~ScopedCeed() { PalaceCeedCallBackend(Destroy(&value)); }
};

using ScopedContext = ScopedCeed<Ceed, CeedDestroy>;
using ScopedBasis = ScopedCeed<CeedBasis, CeedBasisDestroy>;
using ScopedRestriction = ScopedCeed<CeedElemRestriction, CeedElemRestrictionDestroy>;
using ScopedVector = ScopedCeed<CeedVector, CeedVectorDestroy>;
using ScopedQFunction = ScopedCeed<CeedQFunction, CeedQFunctionDestroy>;
using ScopedOperator = ScopedCeed<CeedOperator, CeedOperatorDestroy>;

bool SameContext(Ceed ceed, Ceed object_ceed)
{
  ScopedContext parent;
  PalaceCeedCall(ceed, CeedGetParent(object_ceed, &parent.value));
  return parent.value == ceed;
}

struct PackedField
{
  ScopedRestriction restriction;
  ScopedBasis basis;
  ScopedVector vector;

  void Read(Ceed ceed, CeedOperatorField field)
  {
    PalaceCeedCall(ceed, CeedOperatorFieldGetData(field, nullptr, &restriction.value,
                                                  &basis.value, &vector.value));
  }
};

struct PackedLeaf
{
  PackedField active, qdata;
  CeedInt elements = 0, nodes = 0, qpts = 0;
  CeedEvalMode mode = CEED_EVAL_INTERP;
  bool curl = false;
};

bool FieldMatches(Ceed ceed, CeedQFunctionField field, CeedInt size, CeedEvalMode mode)
{
  CeedInt actual_size;
  CeedEvalMode actual_mode;
  PalaceCeedCall(ceed,
                 CeedQFunctionFieldGetData(field, nullptr, &actual_size, &actual_mode));
  return size == actual_size && mode == actual_mode;
}

std::unique_ptr<PackedLeaf> InspectPackedLeaf(Ceed ceed, CeedOperator op, int size)
{
  // Recognize QData-assembled, non-tensor 3D H(curl)/H(div) volume kernels and H1
  // diffusion kernels. Names identify candidates; the field and restriction checks below
  // establish the layout assumed by the packed kernels. Unsupported leaves remain in the
  // original action.
  bool composite;
  PalaceCeedCall(ceed, CeedOperatorIsComposite(op, &composite));
  if (composite)
  {
    return {};
  }
  bool at_points;
  PalaceCeedCall(ceed, CeedOperatorIsAtPoints(op, &at_points));
  if (at_points || !SameContext(ceed, CeedOperatorReturnCeed(op)))
  {
    return {};
  }
  ScopedQFunction qf;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(op, &qf.value));
  const char *name;
  PalaceCeedCall(ceed, CeedQFunctionGetKernelName(qf.value, &name));
  if (!name)
  {
    return {};
  }
  const bool curl = std::strcmp(name, "f_apply_33") == 0;
  if (!curl && std::strcmp(name, "f_apply_3") != 0)
  {
    return {};
  }
  CeedInt ni, no, qni, qno;
  CeedOperatorField *inputs, *outputs;
  CeedQFunctionField *qinputs, *qoutputs;
  PalaceCeedCall(ceed, CeedOperatorGetFields(op, &ni, &inputs, &no, &outputs));
  PalaceCeedCall(ceed, CeedQFunctionGetFields(qf.value, &qni, &qinputs, &qno, &qoutputs));
  if (ni != (curl ? 3 : 2) || no != (curl ? 2 : 1) || qni != ni || qno != no ||
      !FieldMatches(ceed, qinputs[0], curl ? 18 : 9, CEED_EVAL_NONE))
  {
    return {};
  }
  auto leaf = std::make_unique<PackedLeaf>();
  leaf->curl = curl;
  leaf->active.Read(ceed, inputs[1]);
  auto basis = leaf->active.basis.value;
  auto restriction = leaf->active.restriction.value;
  if (!basis || basis == CEED_BASIS_NONE || !restriction ||
      restriction == CEED_ELEMRESTRICTION_NONE ||
      leaf->active.vector.value != CEED_VECTOR_ACTIVE ||
      !SameContext(ceed, CeedBasisReturnCeed(basis)) ||
      !SameContext(ceed, CeedElemRestrictionReturnCeed(restriction)))
  {
    return {};
  }
  // Vector spaces interpolate their three-component fields; H1 diffusion acts on
  // gradients of a scalar field.
  {
    CeedFESpace space;
    PalaceCeedCall(ceed, CeedBasisGetFESpace(basis, &space));
    leaf->mode = (space == CEED_FE_SPACE_H1) ? CEED_EVAL_GRAD : CEED_EVAL_INTERP;
  }
  // Validate all active field modes/sizes and the actual trial/test handles. Names
  // alone cannot distinguish mass from diffusion, mixed spaces, or boundary operators.
  for (CeedInt j = 0; j < no; j++)
  {
    const auto mode = j == 0 ? leaf->mode : CEED_EVAL_CURL;
    if (!FieldMatches(ceed, qinputs[j + 1], 3, mode) ||
        !FieldMatches(ceed, qoutputs[j], 3, mode))
    {
      return {};
    }
    PackedField in, out;
    in.Read(ceed, inputs[j + 1]);
    out.Read(ceed, outputs[j]);
    if (in.vector.value != CEED_VECTOR_ACTIVE || out.vector.value != CEED_VECTOR_ACTIVE ||
        in.basis.value != basis || out.basis.value != basis ||
        in.restriction.value != restriction || out.restriction.value != restriction)
    {
      return {};
    }
  }
  // One FE field already interpolates to three physical components. Packing adds a
  // second copy for the imaginary coefficients. The restriction clone paths below
  // support unblocked offset maps, including signs and H(curl) orientation transforms.
  bool tensor;
  CeedFESpace space;
  CeedInt dim, components, qcomponents, rcomponents, block;
  CeedSize length, in_length, out_length;
  CeedRestrictionType type;
  PalaceCeedCall(ceed, CeedBasisIsTensor(basis, &tensor));
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &dim));
  PalaceCeedCall(ceed, CeedBasisGetFESpace(basis, &space));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(basis, &components));
  PalaceCeedCall(ceed, CeedBasisGetNumNodes(basis, &leaf->nodes));
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(basis, &leaf->qpts));
  PalaceCeedCall(ceed,
                 CeedBasisGetNumQuadratureComponents(basis, leaf->mode, &qcomponents));
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(restriction, &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(restriction, &rcomponents));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(restriction, &leaf->elements));
  PalaceCeedCall(ceed, CeedElemRestrictionGetBlockSize(restriction, &block));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(restriction, &length));
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op, &in_length, &out_length));
  if (tensor || dim != 3 || components != 1 || qcomponents != 3 || rcomponents != 1 ||
      (space != CEED_FE_SPACE_HCURL && space != CEED_FE_SPACE_HDIV &&
       space != CEED_FE_SPACE_H1) ||
      (curl && space != CEED_FE_SPACE_HCURL) || block != 1 || leaf->elements <= 0 ||
      leaf->qpts <= 0 || length != size || in_length != size || out_length != size ||
      (type != CEED_RESTRICTION_STANDARD && type != CEED_RESTRICTION_ORIENTED &&
       type != CEED_RESTRICTION_CURL_ORIENTED))
  {
    return {};
  }
  if (curl)
  {
    PalaceCeedCall(
        ceed, CeedBasisGetNumQuadratureComponents(basis, CEED_EVAL_CURL, &qcomponents));
    if (qcomponents != 3)
    {
      return {};
    }
  }
  leaf->qdata.Read(ceed, inputs[0]);
  auto qr = leaf->qdata.restriction.value;
  auto qv = leaf->qdata.vector.value;
  if (!qr || qr == CEED_ELEMRESTRICTION_NONE || !qv || qv == CEED_VECTOR_ACTIVE ||
      qv == CEED_VECTOR_NONE || leaf->qdata.basis.value != CEED_BASIS_NONE ||
      !SameContext(ceed, CeedElemRestrictionReturnCeed(qr)) ||
      !SameContext(ceed, CeedVectorReturnCeed(qv)))
  {
    return {};
  }
  // Finalization of the original operator is a caller precondition and validates
  // basis/restriction dimensions and passive QData storage. Reuse the QData vector
  // and its restriction, so in-place value updates remain visible to both operators.
  return leaf;
}

// The second component starts size entries after the first in the L-vector, giving
// [Re(x) | Im(x)] with the same element offsets and orientations for both parts.
// Copy those maps before restoring the borrowed arrays from the original restriction.
void ClonePackedRestriction(Ceed ceed, CeedElemRestriction original, CeedInt elements,
                            CeedInt nodes, int size, CeedElemRestriction *packed)
{
  CeedRestrictionType type;
  const CeedInt *offsets;
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(original, &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetOffsets(original, CEED_MEM_HOST, &offsets));
  const CeedSize length = 2 * static_cast<CeedSize>(size);
  if (type == CEED_RESTRICTION_CURL_ORIENTED)
  {
    const CeedInt8 *orientations;
    PalaceCeedCall(ceed, CeedElemRestrictionGetCurlOrientations(original, CEED_MEM_HOST,
                                                                &orientations));
    PalaceCeedCall(ceed, CeedElemRestrictionCreateCurlOriented(
                             ceed, elements, nodes, 2, size, length, CEED_MEM_HOST,
                             CEED_COPY_VALUES, offsets, orientations, packed));
    PalaceCeedCall(ceed,
                   CeedElemRestrictionRestoreCurlOrientations(original, &orientations));
  }
  else if (type == CEED_RESTRICTION_ORIENTED)
  {
    const bool *orientations;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetOrientations(original, CEED_MEM_HOST, &orientations));
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, elements, nodes, 2, size, length, CEED_MEM_HOST,
                             CEED_COPY_VALUES, offsets, orientations, packed));
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreOrientations(original, &orientations));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, elements, nodes, 2, size, length,
                                                   CEED_MEM_HOST, CEED_COPY_VALUES, offsets,
                                                   packed));
  }
  PalaceCeedCall(ceed, CeedElemRestrictionRestoreOffsets(original, &offsets));
}

// Two copies of the vector basis interpolate the real/imaginary coefficient blocks
// into [physical component][Re/Im][Q] order. This gives six Q-length blocks:
// Re(u_x), Im(u_x), Re(u_y), Im(u_y), Re(u_z), Im(u_z).
void ClonePackedBasis(Ceed ceed, const PackedLeaf &leaf, CeedBasis *packed)
{
  auto original = leaf.active.basis.value;
  CeedFESpace space;
  CeedElemTopology topology;
  const CeedScalar *interp, *derivative, *qref, *weights;
  PalaceCeedCall(ceed, CeedBasisGetFESpace(original, &space));
  PalaceCeedCall(ceed, CeedBasisGetTopology(original, &topology));
  PalaceCeedCall(ceed, CeedBasisGetInterp(original, &interp));
  PalaceCeedCall(ceed, CeedBasisGetQRef(original, &qref));
  PalaceCeedCall(ceed, CeedBasisGetQWeights(original, &weights));
  if (space == CEED_FE_SPACE_HCURL)
  {
    PalaceCeedCall(ceed, CeedBasisGetCurl(original, &derivative));
    PalaceCeedCall(ceed, CeedBasisCreateHcurl(ceed, topology, 2, leaf.nodes, leaf.qpts,
                                              interp, derivative, qref, weights, packed));
  }
  else if (space == CEED_FE_SPACE_H1)
  {
    PalaceCeedCall(ceed, CeedBasisGetGrad(original, &derivative));
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, topology, 2, leaf.nodes, leaf.qpts, interp,
                                           derivative, qref, weights, packed));
  }
  else
  {
    PalaceCeedCall(ceed, CeedBasisGetDiv(original, &derivative));
    PalaceCeedCall(ceed, CeedBasisCreateHdiv(ceed, topology, 2, leaf.nodes, leaf.qpts,
                                             interp, derivative, qref, weights, packed));
  }
}

// The active subset of a finite element space's element list, when the imaginary leaf of a
// pair covers it and the real leaf covers the full list. The geometry data positions the
// subset inside the full list in mesh order.
struct PackedSubset
{
  const CeedGeomFactorData *data = nullptr;
  mfem::Geometry::Type geom = mfem::Geometry::INVALID;
};

// Identify the shared element subset covered by an imaginary leaf whose element restriction
// differs from that of the real leaf it would pair with. The restriction handles cached by
// the finite element space tie the two leaves to its full element list and to the subset of
// that list, in the same mesh order and with the same orientation maps, so the recorded
// positions address the real leaf's elements. Anything else stays unmatched.
PackedSubset MatchPackedSubset(Ceed ceed, const FiniteElementSpace &fespace,
                               const PackedLeaf &real, const PackedLeaf &imag)
{
  CeedElemTopology topology;
  PalaceCeedCall(ceed, CeedBasisGetTopology(real.active.basis.value, &topology));
  const auto geom = GetMfemTopology(topology);
  const auto &geom_data_map = fespace.GetMesh().GetCeedGeomFactorData(ceed);
  const auto it = geom_data_map.find(geom);
  if (it == geom_data_map.end())
  {
    return {};
  }
  const auto &data = it->second;
  if (data.indices.size() != static_cast<std::size_t>(real.elements) ||
      data.active_indices.size() != static_cast<std::size_t>(imag.elements) ||
      data.active_indices.empty() || data.active_indices.size() >= data.indices.size())
  {
    return {};
  }
  // The subset views of the real leaf's QData are derived from its strided layout.
  bool strided;
  PalaceCeedCall(ceed,
                 CeedElemRestrictionIsStrided(real.qdata.restriction.value, &strided));
  if (!strided)
  {
    return {};
  }
  if (real.active.restriction.value !=
          fespace.GetCeedElemRestriction(ceed, geom, data, CeedElementSubset::Full) ||
      imag.active.restriction.value !=
          fespace.GetCeedElemRestriction(ceed, geom, data, CeedElementSubset::Active))
  {
    return {};
  }
  return {&data, geom};
}

// Read the QData of a subset of a leaf's elements without copying any of it. Each offset
// addresses the first component at one quadrature point of a subset element in the leaf's
// full QData vector; the component stride then selects the remaining components at that
// point. The vector layout is the backend's, queried from the strided restriction (as
// [nodes, components, elements] strides) rather than assumed, as for the shared geometry
// factor vector in Mesh::GetCeedGeomFactorData.
void CreateSubsetQDataRestriction(Ceed ceed, CeedElemRestriction full,
                                  const std::vector<int> &positions,
                                  CeedElemRestriction *subset)
{
  bool strided;
  CeedInt num_qpts, num_comp, layout[3];
  CeedSize length;
  PalaceCeedCall(ceed, CeedElemRestrictionIsStrided(full, &strided));
  MFEM_ASSERT(strided, "Packed QData subset views require a strided QData layout!");
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(full, &num_qpts));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(full, &num_comp));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(full, &length));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLLayout(full, layout));
  std::vector<CeedInt> offsets(positions.size() * num_qpts);
  for (std::size_t k = 0; k < positions.size(); k++)
  {
    const CeedInt elem_offset = positions[k] * layout[2];
    for (CeedInt q = 0; q < num_qpts; q++)
    {
      offsets[k * num_qpts + q] = elem_offset + q * layout[0];
    }
  }
  PalaceCeedCall(
      ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(positions.size()),
                                      num_qpts, num_comp, layout[1], length, CEED_MEM_HOST,
                                      CEED_COPY_VALUES, offsets.data(), subset));
}

void AddPackedPair(Ceed ceed, const PackedLeaf &real, const PackedLeaf &imag, int size,
                   const PackedSubset &subset, Operator &packed)
{
  ScopedRestriction restriction, qdata;
  ScopedBasis basis;
  ScopedQFunction qf;
  ScopedOperator op;
  // The pair is applied over the imaginary leaf's elements: either the same element list as
  // the real leaf, or the active subset of it. Both leaves share the basis, so the node
  // count matches, and the cloned orientation maps coincide with the real leaf's on the
  // subset elements because both restrictions derive from the same space's element list.
  ClonePackedRestriction(ceed, imag.active.restriction.value, imag.elements, imag.nodes,
                         size, &restriction.value);
  ClonePackedBasis(ceed, real, &basis.value);
  if (subset.data)
  {
    // The real coefficients of the subset elements, in place in the real leaf's QData.
    CreateSubsetQDataRestriction(ceed, real.qdata.restriction.value,
                                 subset.data->active_indices, &qdata.value);
  }
  // Both kernels implement complex mass. Only f_apply_complex_33 adds curl-curl,
  // with a real tensor acting independently on Re(x) and Im(x). The pairing loop
  // must therefore leave imaginary curl leaves in the remainder.
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                           ceed, 1, real.curl ? f_apply_complex_33 : f_apply_complex_3,
                           real.curl ? PalaceQFunctionRelativePath(f_apply_complex_33_loc)
                                     : PalaceQFunctionRelativePath(f_apply_complex_3_loc),
                           &qf.value));
  PalaceCeedCall(ceed,
                 CeedQFunctionAddInput(qf.value, "qr", real.curl ? 18 : 9, CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "qi", 9, CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "u", 6, real.mode));
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "v", 6, real.mode));
  if (real.curl)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "curl_u", 6, CEED_EVAL_CURL));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "curl_v", 6, CEED_EVAL_CURL));
  }
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &op.value));
  PalaceCeedCall(
      ceed, CeedOperatorSetField(op.value, "qr",
                                 qdata.value ? qdata.value : real.qdata.restriction.value,
                                 CEED_BASIS_NONE, real.qdata.vector.value));
  PalaceCeedCall(ceed, CeedOperatorSetField(op.value, "qi", imag.qdata.restriction.value,
                                            CEED_BASIS_NONE, imag.qdata.vector.value));
  for (const auto *field : {"u", "v", "curl_u", "curl_v"})
  {
    if (real.curl || field[0] != 'c')
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(op.value, field, restriction.value,
                                                basis.value, CEED_VECTOR_ACTIVE));
    }
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(op.value));
  packed.AddSubOperator(std::exchange(op.value, nullptr));
}

// Pack a leaf that has no partner in the other part: a real leaf acts identically on
// [Re(x) | Im(x)], an imaginary (mass) leaf maps them to [-Im | Re] parts of i A x. The
// leaf's QData vector is reused, either over all of its elements or over a subset of them
// given by an element restriction and the matching QData positions.
void AddPackedSingle(Ceed ceed, const PackedLeaf &leaf, bool imag, int size,
                     CeedElemRestriction subset_restriction, CeedInt subset_elements,
                     CeedElemRestriction subset_qdata, Operator &packed)
{
  ScopedRestriction restriction;
  ScopedBasis basis;
  ScopedQFunction qf;
  ScopedOperator op;
  MFEM_ASSERT(!imag || !leaf.curl, "Imaginary curl leaves are not packed!");
  ClonePackedRestriction(
      ceed, subset_restriction ? subset_restriction : leaf.active.restriction.value,
      subset_restriction ? subset_elements : leaf.elements, leaf.nodes, size,
      &restriction.value);
  ClonePackedBasis(ceed, leaf, &basis.value);
  if (imag)
  {
    PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                             ceed, 1, f_apply_imag_3,
                             PalaceQFunctionRelativePath(f_apply_imag_3_loc), &qf.value));
  }
  else
  {
    PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                             ceed, 1, leaf.curl ? f_apply_real_33 : f_apply_real_3,
                             leaf.curl ? PalaceQFunctionRelativePath(f_apply_real_33_loc)
                                       : PalaceQFunctionRelativePath(f_apply_real_3_loc),
                             &qf.value));
  }
  PalaceCeedCall(ceed,
                 CeedQFunctionAddInput(qf.value, "q", leaf.curl ? 18 : 9, CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "u", 6, leaf.mode));
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "v", 6, leaf.mode));
  if (leaf.curl)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "curl_u", 6, CEED_EVAL_CURL));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "curl_v", 6, CEED_EVAL_CURL));
  }
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &op.value));
  PalaceCeedCall(
      ceed, CeedOperatorSetField(op.value, "q",
                                 subset_qdata ? subset_qdata : leaf.qdata.restriction.value,
                                 CEED_BASIS_NONE, leaf.qdata.vector.value));
  for (const auto *field : {"u", "v", "curl_u", "curl_v"})
  {
    if (leaf.curl || field[0] != 'c')
    {
      PalaceCeedCall(ceed, CeedOperatorSetField(op.value, field, restriction.value,
                                                basis.value, CEED_VECTOR_ACTIVE));
    }
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(op.value));
  packed.AddSubOperator(std::exchange(op.value, nullptr));
}

// Pack a real leaf's action over the elements outside the subset its pair was packed over,
// so that the two together cover each of its elements exactly once. The leaf's bases and
// QData vector are reused; only the element restrictions are those of the complement.
void AddPackedComplement(Ceed ceed, const FiniteElementSpace &fespace,
                         const PackedSubset &subset, const PackedLeaf &leaf, int size,
                         Operator &packed)
{
  // Derive the complement only for this CPU packing path rather than storing it on every
  // mesh and finite element space.
  std::vector<int> positions, elements;
  positions.reserve(subset.data->indices.size() - subset.data->active_indices.size());
  elements.reserve(positions.capacity());
  auto active = subset.data->active_indices.begin();
  for (std::size_t i = 0; i < subset.data->indices.size(); i++)
  {
    if (active != subset.data->active_indices.end() && *active == static_cast<int>(i))
    {
      ++active;
    }
    else
    {
      positions.push_back(static_cast<int>(i));
      elements.push_back(subset.data->indices[i]);
    }
  }
  MFEM_ASSERT(!positions.empty(), "Empty packed-operator complement!");

  ScopedRestriction restriction, qdata;
  restriction.value =
      FiniteElementSpace::BuildCeedElemRestriction(fespace, ceed, subset.geom, elements);
  CreateSubsetQDataRestriction(ceed, leaf.qdata.restriction.value, positions, &qdata.value);
  AddPackedSingle(ceed, leaf, false, size, restriction.value,
                  static_cast<CeedInt>(positions.size()), qdata.value, packed);
}

// Element assembly (EA) of leaves for CPU forward applications. The dense element
// matrices of a leaf, including any orientation transforms of its restriction, are stored
// in L-vector coordinates: an application gathers each element's coefficients, multiplies
// them by the element matrix, and scatters the result. Real matrices act on both
// [Re(x) | Im(x)]; imaginary ones contribute i A x.
struct ElementMatrices
{
  CeedInt elements = 0, nodes = 0;
  bool symmetric = false;
  std::vector<CeedInt> indices;
  std::vector<CeedScalar> values;
};

// Element assembly replaces the streaming of passive QData and the basis contractions of
// partial assembly by the streaming of the element matrices, with far fewer floating point
// operations. It is considered for single-component leaves of non-tensor (simplex) H1 and
// boundary bases with up to 45 nodes per element (H1 p <= 4 on tetrahedra); larger and
// tensor-product elements keep partial assembly.
bool IsElementAssemblyLeaf(Ceed ceed, CeedOperator op, int size, int num_rhs)
{
  constexpr CeedInt max_nodes = 45;
  bool composite, at_points;
  PalaceCeedCall(ceed, CeedOperatorIsComposite(op, &composite));
  if (composite)
  {
    return false;
  }
  PalaceCeedCall(ceed, CeedOperatorIsAtPoints(op, &at_points));
  if (at_points || !SameContext(ceed, CeedOperatorReturnCeed(op)))
  {
    return false;
  }
  ScopedRestriction in, out;
  ScopedBasis basis_in, basis_out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestrictions(op, &in.value, &out.value));
  PalaceCeedCall(ceed, CeedOperatorGetActiveBases(op, &basis_in.value, &basis_out.value));
  if (!in.value || in.value != out.value || !basis_in.value ||
      basis_in.value == CEED_BASIS_NONE || basis_in.value != basis_out.value)
  {
    return false;
  }
  bool tensor;
  CeedInt nodes, components;
  CeedSize length;
  PalaceCeedCall(ceed, CeedBasisIsTensor(basis_in.value, &tensor));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(in.value, &nodes));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(in.value, &components));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(in.value, &length));
  CeedFESpace space;
  CeedInt dim;
  PalaceCeedCall(ceed, CeedBasisGetFESpace(basis_in.value, &space));
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis_in.value, &dim));
  if (tensor || components != 1 || length != size || nodes <= 0 || nodes > max_nodes)
  {
    return false;
  }
  // Leaves with a curl or a divergence (below) keep partial assembly: element matrices of
  // curl-curl terms do not preserve the discrete curl-gradient cancellation of the factored
  // action, which perturbs ill-conditioned (e.g. near-resonant driven) solves beyond
  // rounding. So do all volume H(curl) and H(div) leaves in 3D, whose mass terms are mostly
  // bandwidth bound either way.
  if (space != CEED_FE_SPACE_H1 && dim == 3)
  {
    return false;
  }
  // Per element, partial assembly streams the passive QData and the element dofs and
  // performs the basis contractions for each right hand side, while element assembly
  // streams the (symmetric) element matrix once for all right hand sides. Choose element
  // assembly when it moves fewer bytes, counting contraction flops at the flop-per-byte
  // balance of a fully loaded CPU socket.
  constexpr double flops_per_byte = 4.0;
  CeedInt num_qpts;
  PalaceCeedCall(ceed, CeedOperatorGetNumQuadraturePoints(op, &num_qpts));
  ScopedQFunction qf;
  CeedInt num_in, num_out;
  CeedQFunctionField *qf_in, *qf_out;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(op, &qf.value));
  PalaceCeedCall(ceed,
                 CeedQFunctionGetFields(qf.value, &num_in, &qf_in, &num_out, &qf_out));
  double pa_bytes = 4.0 * nodes, pa_flops = 0.0;
  for (CeedInt i = 0; i < num_in + num_out; i++)
  {
    CeedEvalMode mode;
    CeedInt field_size;
    PalaceCeedCall(ceed,
                   CeedQFunctionFieldGetData(i < num_in ? qf_in[i] : qf_out[i - num_in],
                                             nullptr, &field_size, &mode));
    if (mode == CEED_EVAL_CURL || mode == CEED_EVAL_DIV)
    {
      return false;
    }
    if (mode == CEED_EVAL_NONE)
    {
      pa_bytes += 8.0 * field_size * num_qpts;
    }
    else if (mode != CEED_EVAL_WEIGHT)
    {
      pa_flops += 2.0 * field_size * num_qpts * nodes;
    }
  }
  const double ea_bytes = 4.0 * nodes * (nodes + 1) + 4.0 * nodes;
  return ea_bytes <= pa_bytes + num_rhs * pa_flops / flops_per_byte;
}

// A new operator with the QFunction and fields of a leaf. Assembling it instead of the leaf
// leaves the leaf without backend data: the CPU backends allocate element vectors for all
// fields on first use, including a copy of the passive inputs, which an element assembled
// leaf does not need unless libCEED applies or assembles it later.
CeedOperator CopyLeaf(CeedOperator leaf)
{
  Ceed ceed = CeedOperatorReturnCeed(leaf);
  ScopedQFunction qf;
  ScopedOperator copy;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(leaf, &qf.value));
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &copy.value));
  CeedInt num_in, num_out;
  CeedOperatorField *in, *out;
  PalaceCeedCall(ceed, CeedOperatorGetFields(leaf, &num_in, &in, &num_out, &out));
  for (CeedInt i = 0; i < num_in + num_out; i++)
  {
    const char *name;
    ScopedRestriction restriction;
    ScopedBasis basis;
    ScopedVector vector;
    PalaceCeedCall(ceed, CeedOperatorFieldGetData(i < num_in ? in[i] : out[i - num_in],
                                                  &name, &restriction.value, &basis.value,
                                                  &vector.value));
    PalaceCeedCall(ceed, CeedOperatorSetField(copy.value, name, restriction.value,
                                              basis.value, vector.value));
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(copy.value));
  return std::exchange(copy.value, nullptr);
}

std::unique_ptr<ElementMatrices> AssembleElementMatrices(Ceed ceed, CeedOperator op,
                                                         int size)
{
  ScopedRestriction in, out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestrictions(op, &in.value, &out.value));
  auto em = std::make_unique<ElementMatrices>();
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(in.value, &em->elements));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(in.value, &em->nodes));
  const std::size_t P = em->nodes, ne = em->elements;
  CeedSize nnz;
  CeedInt *rows, *cols;
  ScopedVector vals;
  {
    ScopedOperator leaf;
    leaf.value = CopyLeaf(op);
    PalaceCeedCall(ceed,
                   CeedOperatorLinearAssembleSymbolic(leaf.value, &nnz, &rows, &cols));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, nnz, &vals.value));
    PalaceCeedCall(ceed, CeedOperatorLinearAssemble(leaf.value, vals.value));
  }
  // The entries come in element blocks, each with the element's row coefficients in order
  // and the same column coefficients, since the leaf uses one active restriction.
  bool valid = (static_cast<std::size_t>(nnz) == ne * P * P);
  em->indices.resize(ne * P);
  for (std::size_t e = 0; valid && e < ne; e++)
  {
    const CeedInt *r = rows + e * P * P, *c = cols + e * P * P;
    for (std::size_t i = 0; valid && i < P; i++)
    {
      em->indices[e * P + i] = r[i * P];
      if (r[i * P] < 0 || r[i * P] >= size)
      {
        valid = false;
        break;
      }
      for (std::size_t j = 0; j < P; j++)
      {
        if (r[i * P + j] != r[i * P] || c[i * P + j] != r[j * P])
        {
          valid = false;
          break;
        }
      }
    }
  }
  // The sparsity arrays are allocated by libCEED with calloc.
  std::free(rows);
  std::free(cols);
  if (!valid)
  {
    return {};
  }
  const CeedScalar *v;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals.value, CEED_MEM_HOST, &v));
  // Store the upper triangle when every element matrix is symmetric up to rounding.
  bool symmetric = true;
  for (std::size_t e = 0; symmetric && e < ne; e++)
  {
    const CeedScalar *a = v + e * P * P;
    CeedScalar scale = 0.0;
    for (std::size_t k = 0; k < P * P; k++)
    {
      scale = std::max(scale, std::abs(a[k]));
    }
    const CeedScalar tol = 1.0e-13 * scale;
    for (std::size_t i = 0; symmetric && i < P; i++)
    {
      for (std::size_t j = i + 1; j < P; j++)
      {
        if (std::abs(a[i * P + j] - a[j * P + i]) > tol)
        {
          symmetric = false;
          break;
        }
      }
    }
  }
  em->symmetric = symmetric;
  if (symmetric)
  {
    em->values.resize(ne * P * (P + 1) / 2);
    CeedScalar *u = em->values.data();
    for (std::size_t e = 0; e < ne; e++)
    {
      const CeedScalar *a = v + e * P * P;
      for (std::size_t i = 0; i < P; i++)
      {
        *u++ = a[i * P + i];
        for (std::size_t j = i + 1; j < P; j++)
        {
          *u++ = 0.5 * (a[i * P + j] + a[j * P + i]);
        }
      }
    }
  }
  else
  {
    em->values.assign(v, v + nnz);
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals.value, &v));
  return em;
}

// y += A x (real) or y += i A x (imaginary) for x = [xr | xi] and y = [yr | yi].
template <bool Symmetric, bool Imag>
void ApplyElementMatrices(const ElementMatrices &em, const CeedScalar *__restrict__ xr,
                          const CeedScalar *__restrict__ xi, CeedScalar *__restrict__ yr,
                          CeedScalar *__restrict__ yi, CeedScalar *__restrict__ work)
{
  const std::size_t P = em.nodes;
  const std::size_t S = Symmetric ? P * (P + 1) / 2 : P * P;
  CeedScalar *__restrict__ xe = work, *__restrict__ ye = work + 2 * P;
  for (std::size_t e = 0; e < static_cast<std::size_t>(em.elements); e++)
  {
    const CeedInt *id = em.indices.data() + e * P;
    const CeedScalar *__restrict__ a = em.values.data() + e * S;
    for (std::size_t j = 0; j < P; j++)
    {
      xe[2 * j] = xr[id[j]];
      xe[2 * j + 1] = xi[id[j]];
    }
    if constexpr (Symmetric)
    {
      for (std::size_t j = 0; j < 2 * P; j++)
      {
        ye[j] = 0.0;
      }
      for (std::size_t i = 0; i < P; i++)
      {
        const CeedScalar x0 = xe[2 * i], x1 = xe[2 * i + 1];
        CeedScalar s0 = a[0] * x0, s1 = a[0] * x1;
        a++;
        for (std::size_t j = i + 1; j < P; j++)
        {
          const CeedScalar aij = a[j - i - 1];
          s0 += aij * xe[2 * j];
          s1 += aij * xe[2 * j + 1];
          ye[2 * j] += aij * x0;
          ye[2 * j + 1] += aij * x1;
        }
        a += P - i - 1;
        ye[2 * i] += s0;
        ye[2 * i + 1] += s1;
      }
    }
    else
    {
      for (std::size_t i = 0; i < P; i++)
      {
        CeedScalar s0 = 0.0, s1 = 0.0;
        for (std::size_t j = 0; j < P; j++)
        {
          s0 += a[i * P + j] * xe[2 * j];
          s1 += a[i * P + j] * xe[2 * j + 1];
        }
        ye[2 * i] = s0;
        ye[2 * i + 1] = s1;
      }
    }
    for (std::size_t i = 0; i < P; i++)
    {
      if constexpr (Imag)
      {
        yr[id[i]] -= ye[2 * i + 1];
        yi[id[i]] += ye[2 * i];
      }
      else
      {
        yr[id[i]] += ye[2 * i];
        yi[id[i]] += ye[2 * i + 1];
      }
    }
  }
}

// y += A x for a single real vector.
template <bool Symmetric>
void ApplyElementMatrices(const ElementMatrices &em, const CeedScalar *__restrict__ x,
                          CeedScalar *__restrict__ y, CeedScalar *__restrict__ work)
{
  const std::size_t P = em.nodes;
  const std::size_t S = Symmetric ? P * (P + 1) / 2 : P * P;
  CeedScalar *__restrict__ xe = work, *__restrict__ ye = work + P;
  for (std::size_t e = 0; e < static_cast<std::size_t>(em.elements); e++)
  {
    const CeedInt *id = em.indices.data() + e * P;
    const CeedScalar *__restrict__ a = em.values.data() + e * S;
    for (std::size_t j = 0; j < P; j++)
    {
      xe[j] = x[id[j]];
    }
    if constexpr (Symmetric)
    {
      for (std::size_t j = 0; j < P; j++)
      {
        ye[j] = 0.0;
      }
      for (std::size_t i = 0; i < P; i++)
      {
        const CeedScalar xi = xe[i];
        CeedScalar s = a[0] * xi;
        a++;
        for (std::size_t j = i + 1; j < P; j++)
        {
          const CeedScalar aij = a[j - i - 1];
          s += aij * xe[j];
          ye[j] += aij * xi;
        }
        a += P - i - 1;
        ye[i] += s;
      }
    }
    else
    {
      for (std::size_t i = 0; i < P; i++)
      {
        CeedScalar s = 0.0;
        for (std::size_t j = 0; j < P; j++)
        {
          s += a[i * P + j] * xe[j];
        }
        ye[i] = s;
      }
    }
    for (std::size_t i = 0; i < P; i++)
    {
      y[id[i]] += ye[i];
    }
  }
}

}  // namespace

struct Operator::ElementAssembly
{
  // Element-assembled children of the composite, and the composite of the others.
  std::vector<CeedInt> children;
  std::vector<std::shared_ptr<const ElementMatrices>> leaves;
  CeedOperator remainder = nullptr;
  mutable Vector work;

  const std::shared_ptr<const ElementMatrices> *Find(CeedInt child) const
  {
    for (std::size_t k = 0; k < children.size(); k++)
    {
      if (children[k] == child)
      {
        return &leaves[k];
      }
    }
    return nullptr;
  }
};

namespace
{

// An element-assembled leaf of the real or imaginary part of a complex operator.
struct ElementLeaf
{
  std::shared_ptr<const ElementMatrices> matrices;
  bool imag;
};

void AddMultElementMatrices(const std::vector<ElementLeaf> &ea, const CeedScalar *xr,
                            const CeedScalar *xi, CeedScalar *yr, CeedScalar *yi,
                            Vector &work)
{
  for (const auto &[em, imag] : ea)
  {
    if (work.Size() < 4 * em->nodes)
    {
      work.SetSize(4 * em->nodes);
    }
    CeedScalar *w = work.HostWrite();
    if (em->symmetric)
    {
      imag ? ApplyElementMatrices<true, true>(*em, xr, xi, yr, yi, w)
           : ApplyElementMatrices<true, false>(*em, xr, xi, yr, yi, w);
    }
    else
    {
      imag ? ApplyElementMatrices<false, true>(*em, xr, xi, yr, yi, w)
           : ApplyElementMatrices<false, false>(*em, xr, xi, yr, yi, w);
    }
  }
}

// The base owns the finalized originals for diagonal, transpose, full assembly, and
// coarsening. Only forward applications use the packed composite and its remainder.
class PackedComplexOperator final : public ComplexWrapperOperator
{
private:
  std::unique_ptr<Operator> packed;
  std::unique_ptr<ComplexWrapperOperator> remainder;
  std::vector<ElementLeaf> ea;
  // Persistent L-vectors back the real/imaginary views. Forward calls overwrite this
  // scratch, so applications on the same instance must not overlap.
  mutable Vector input, output, work;
  mutable ComplexVector input_view, output_view;

  const ComplexVector &Apply(const ComplexVector &x) const
  {
    // Stage x before touching the caller's output. Mult/AddMult consume the private
    // result only after all actions finish, including when their x and y alias.
    input_view = x;
    // The remainder writes the caller-facing output directly; the packed and element
    // assembled actions then accumulate on top with their supported +1 coefficient. All
    // read the private staged input, so there is no output zeroing or trailing AXPY.
    if (remainder)
    {
      remainder->Mult(input_view, output_view);
      if (packed)
      {
        packed->AddMult(input, output);
      }
    }
    else if (packed)
    {
      packed->Mult(input, output);
    }
    else
    {
      output = 0.0;
    }
    if (!ea.empty())
    {
      const CeedScalar *xr = input.HostRead();
      CeedScalar *yr = output.HostReadWrite();
      AddMultElementMatrices(ea, xr, xr + width, yr, yr + height, work);
    }
    return output_view;
  }

  // Element matrices and the remainder read the parts of x in place, unless they alias
  // the output.
  static bool Overlaps(const Vector &a, const Vector &b)
  {
    const auto *pa = a.HostRead(), *pb = b.HostRead();
    return pa < pb + b.Size() && pb < pa + a.Size();
  }

  bool ApplyInPlace(const ComplexVector &x, ComplexVector &y) const
  {
    if (packed || Overlaps(x.Real(), y.Real()) || Overlaps(x.Real(), y.Imag()) ||
        Overlaps(x.Imag(), y.Real()) || Overlaps(x.Imag(), y.Imag()))
    {
      return false;
    }
    if (remainder)
    {
      remainder->Mult(x, y);
    }
    else
    {
      y = 0.0;
    }
    AddMultElementMatrices(ea, x.Real().HostRead(), x.Imag().HostRead(),
                           y.Real().HostReadWrite(), y.Imag().HostReadWrite(), work);
    return true;
  }

public:
  PackedComplexOperator(std::unique_ptr<palace::Operator> &&Ar,
                        std::unique_ptr<palace::Operator> &&Ai,
                        std::unique_ptr<Operator> &&packed,
                        std::unique_ptr<ComplexWrapperOperator> &&remainder,
                        std::vector<ElementLeaf> &&ea)
    : ComplexWrapperOperator(std::move(Ar), std::move(Ai)), packed(std::move(packed)),
      remainder(std::move(remainder)), ea(std::move(ea)), input(2 * width),
      output(2 * height)
  {
    // Set the owning buffers' memory policy before MakeRef allocates storage and
    // creates the real/imaginary aliases.
    input.UseDevice(true);
    output.UseDevice(true);
    input_view.MakeRef(input, 0, width);
    output_view.MakeRef(output, 0, height);
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    MFEM_ASSERT(x.Size() == width && y.Size() == height,
                "Invalid dimensions for packed complex application!");
    if (!ApplyInPlace(x, y))
    {
      y = Apply(x);
    }
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               std::complex<double> a) const override
  {
    MFEM_ASSERT(x.Size() == width && y.Size() == height,
                "Invalid dimensions for packed complex application!");
    if (a != std::complex<double>{0.0})
    {
      // LocalOperator() and CreateComplexOperator() expose AddMult to local callers.
      // The inherited real-scale path requests unsupported negative CEED coefficients.
      y.AXPY(a, Apply(x));
    }
  }
};

// Apply compatible terms of (Ar + i Ai)x through shared two-component restriction
// and basis work for the real and imaginary parts, avoiding separate wrapper
// applications. The packed QFunctions preserve the original mass/curl algebra.
// The finite element space the leaves were assembled from identifies its own element
// restrictions, so that a pair over one of its element subsets can be recognized.
std::unique_ptr<ComplexWrapperOperator>
PackComplexOperator(std::unique_ptr<palace::Operator> &Ar,
                    std::unique_ptr<palace::Operator> &Ai,
                    const FiniteElementSpace &fespace)
{
  // Require one host CEED context and one OpenMP thread per MPI process;
  // utils::InParallel() checks for an active OpenMP region.
  // Equal square spaces let the remainder reuse input storage; doubled dimensions
  // must fit int. Multiplicity scaling acts outside the CEED children and would be
  // lost by packing.
  const auto *real = dynamic_cast<const Operator *>(Ar.get());
  const auto *imag = dynamic_cast<const Operator *>(Ai.get());
  if ((!real && !imag) || (Ar && !real) || (Ai && !imag) || internal::NumCeeds() != 1 ||
      utils::GetMaxThreads() > 1 || utils::InParallel() ||
      mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    return {};
  }
  const auto *first = real ? real : imag;
  for (const auto *part : {real, imag})
  {
    // Do not bypass an unknown subclass's application semantics.
    if (part &&
        (part->Size() != 1 || part->HasDofMultiplicity() ||
         part->Width() != part->Height() || part->Width() != first->Width() ||
         (typeid(*part) != typeid(Operator) && typeid(*part) != typeid(SymmetricOperator))))
    {
      return {};
    }
  }
  if (first->Width() <= 0 || first->Width() > std::numeric_limits<int>::max() / 2)
  {
    return {};
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  CeedMemType mem;
  const char *resource;
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  if (mem != CEED_MEM_HOST || std::strncmp(resource, "/cpu/", 5) != 0 ||
      (real && !SameContext(ceed, CeedOperatorReturnCeed((*real)[0]))) ||
      (imag && !SameContext(ceed, CeedOperatorReturnCeed((*imag)[0]))))
  {
    return {};
  }
  const int size = first->Width();
  // Pair only direct children of Palace's composites. Keep unmatched nested
  // composites intact in the remainder rather than flattening them, preserving
  // their action and libCEED's composite size limits.
  CeedInt count[2] = {0, 0};
  CeedOperator *children[2] = {nullptr, nullptr};
  for (int part = 0; part < 2; part++)
  {
    const auto *part_op = part == 0 ? real : imag;
    if (part_op)
    {
      auto op = (*part_op)[0];
      PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op, &count[part]));
      PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op, &children[part]));
    }
  }
  // The original parts apply their element-assembled leaves themselves, so that their own
  // applications (used for example in projections onto reduced bases) and the complex
  // action share the same element matrices. The choice of leaves is made for the complex
  // action, which applies each part to two vectors.
  const Operator::ElementAssembly *part_ea[2] = {nullptr, nullptr};
  for (int part = 0; part < 2; part++)
  {
    auto *part_op = dynamic_cast<Operator *>((part == 0 ? Ar : Ai).get());
    if (part_op && part_op->UseElementAssembly(2))
    {
      part_ea[part] = part_op->GetElementAssembly();
    }
  }
  std::vector<std::unique_ptr<PackedLeaf>> fields[2];
  std::vector<bool> matched[2];
  std::vector<ElementLeaf> ea;
  for (int part = 0; part < 2; part++)
  {
    matched[part].resize(count[part], false);
    for (CeedInt j = 0; j < count[part]; j++)
    {
      const auto *em = part_ea[part] ? part_ea[part]->Find(j) : nullptr;
      if (em)
      {
        ea.push_back({*em, part == 1});
        matched[part][j] = true;
        fields[part].emplace_back();
        continue;
      }
      fields[part].push_back(InspectPackedLeaf(ceed, children[part][j], size));
    }
  }
  // A FiniteElementSpace caches basis/restriction handles per context and geometry.
  // Requiring both identities ties each pair to the same quadrature, element ordering,
  // and orientation maps. Separately constructed equivalents remain unmatched.
  std::unique_ptr<Operator> packed;
  std::vector<std::pair<std::size_t, PackedSubset>> complements;
  for (std::size_t r = 0; r < fields[0].size(); r++)
  {
    if (!fields[0][r] || matched[0][r])
    {
      continue;
    }
    for (std::size_t i = 0; i < fields[1].size(); i++)
    {
      if (matched[1][i] || !fields[1][i] || fields[1][i]->curl ||
          fields[0][r]->active.basis.value != fields[1][i]->active.basis.value)
      {
        continue;
      }
      // The two leaves either cover the same elements, or the imaginary leaf covers the
      // space's shared active subset of the real leaf's full element list. In the latter
      // case the pair is applied over the subset and the real leaf's remaining elements
      // are applied on their own below.
      PackedSubset subset;
      if (fields[0][r]->active.restriction.value != fields[1][i]->active.restriction.value)
      {
        subset = MatchPackedSubset(ceed, fespace, *fields[0][r], *fields[1][i]);
        if (!subset.data)
        {
          continue;
        }
      }
      if (!packed)
      {
        packed = std::make_unique<Operator>(2 * size, 2 * size);
      }
      AddPackedPair(ceed, *fields[0][r], *fields[1][i], size, subset, *packed);
      if (subset.data)
      {
        complements.emplace_back(r, subset);
      }
      matched[0][r] = matched[1][i] = true;
      break;
    }
  }
  // A real leaf packed over a subset still owes the elements that subset leaves out, and
  // recognized leaves without a partner act on both parts of x on their own. Imaginary
  // curl leaves stay in the remainder, like for pairing.
  for (const auto &[r, subset] : complements)
  {
    AddPackedComplement(ceed, fespace, subset, *fields[0][r], size, *packed);
  }
  for (int part = 0; part < 2; part++)
  {
    for (std::size_t j = 0; j < fields[part].size(); j++)
    {
      if (matched[part][j] || !fields[part][j] || (part == 1 && fields[part][j]->curl))
      {
        continue;
      }
      if (!packed)
      {
        packed = std::make_unique<Operator>(2 * size, 2 * size);
      }
      AddPackedSingle(ceed, *fields[part][j], part == 1, size, nullptr, 0, nullptr,
                      *packed);
      matched[part][j] = true;
    }
  }
  if (!packed && ea.empty())
  {
    return {};  // This also handles ranks with no local volume terms.
  }
  // Keep every unmatched child exactly once. In particular, differing real/imaginary
  // boundary terms must not disable fusion of compatible domain contributions.
  // These private remainders are only applied forward, so their transpose composites
  // intentionally stay empty. Transposes use the retained original operators.
  std::unique_ptr<palace::Operator> remainder_parts[2];
  for (int part = 0; part < 2; part++)
  {
    std::unique_ptr<Operator> remainder;
    for (CeedInt j = 0; j < count[part]; j++)
    {
      if (!matched[part][j])
      {
        if (!remainder)
        {
          remainder = std::make_unique<Operator>(size, size);
        }
        CeedOperator copy = nullptr;
        PalaceCeedCall(ceed, CeedOperatorReferenceCopy(children[part][j], &copy));
        remainder->AddSubOperator(copy);
      }
    }
    if (remainder)
    {
      remainder->Finalize();
      remainder_parts[part] = std::move(remainder);
    }
  }
  std::unique_ptr<ComplexWrapperOperator> remainder;
  if (remainder_parts[0] || remainder_parts[1])
  {
    remainder = std::make_unique<ComplexWrapperOperator>(std::move(remainder_parts[0]),
                                                         std::move(remainder_parts[1]));
  }
  if (packed)
  {
    packed->Finalize();
  }
  return std::make_unique<PackedComplexOperator>(
      std::move(Ar), std::move(Ai), std::move(packed), std::move(remainder), std::move(ea));
}

}  // namespace

Operator::~Operator()
{
  if (ea && ea->remainder)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(ea->remainder, &ceed));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&ea->remainder));
  }
  PalacePragmaOmp(parallel if (op.size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op[id]));
    PalaceCeedCall(ceed, CeedOperatorDestroy(&op_t[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&u[id]));
    PalaceCeedCall(ceed, CeedVectorDestroy(&v[id]));
  }
}

bool Operator::UseElementAssembly(int num_rhs)
{
  if (ea)
  {
    return true;
  }
  if (op.size() != 1 || internal::NumCeeds() != 1 || utils::GetMaxThreads() > 1 ||
      utils::InParallel() || mfem::Device::Allows(mfem::Backend::DEVICE_MASK) ||
      HasDofMultiplicity() || width != height || width <= 0)
  {
    return false;
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  CeedMemType mem;
  const char *resource;
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  if (mem != CEED_MEM_HOST || std::strncmp(resource, "/cpu/", 5) != 0 ||
      !SameContext(ceed, CeedOperatorReturnCeed(op[0])))
  {
    return false;
  }
  CeedInt count;
  CeedOperator *children;
  PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[0], &count));
  PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[0], &children));
  auto data = std::make_unique<ElementAssembly>();
  std::vector<CeedOperator> others;
  for (CeedInt j = 0; j < count; j++)
  {
    std::unique_ptr<ElementMatrices> em;
    if (IsElementAssemblyLeaf(ceed, children[j], width, num_rhs))
    {
      em = AssembleElementMatrices(ceed, children[j], width);
    }
    if (em)
    {
      data->children.push_back(j);
      data->leaves.push_back(std::move(em));
    }
    else
    {
      others.push_back(children[j]);
    }
  }
  if (data->leaves.empty())
  {
    return false;
  }
  if (!others.empty())
  {
    PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &data->remainder));
    for (auto child : others)
    {
      PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(data->remainder, child));
    }
    PalaceCeedCall(ceed, CeedOperatorCheckReady(data->remainder));
  }
  ea = std::move(data);
  return true;
}

void Operator::ElementAssemblyMult(const Vector &x, Vector &y, bool add) const
{
  if (!add)
  {
    y = 0.0;
  }
  if (ea->remainder)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(ea->remainder, &ceed));
    PalaceCeedCall(ceed, CeedVectorSetArray(u[0], CEED_MEM_HOST, CEED_USE_POINTER,
                                            const_cast<CeedScalar *>(x.HostRead())));
    PalaceCeedCall(
        ceed, CeedVectorSetArray(v[0], CEED_MEM_HOST, CEED_USE_POINTER, y.HostReadWrite()));
    PalaceCeedCall(ceed,
                   CeedOperatorApplyAdd(ea->remainder, u[0], v[0], CEED_REQUEST_IMMEDIATE));
    PalaceCeedCall(ceed, CeedVectorTakeArray(u[0], CEED_MEM_HOST, nullptr));
    PalaceCeedCall(ceed, CeedVectorTakeArray(v[0], CEED_MEM_HOST, nullptr));
  }
  const CeedScalar *xd = x.HostRead();
  CeedScalar *yd = y.HostReadWrite();
  for (const auto &em : ea->leaves)
  {
    if (ea->work.Size() < 2 * em->nodes)
    {
      ea->work.SetSize(2 * em->nodes);
    }
    CeedScalar *w = ea->work.HostWrite();
    em->symmetric ? ApplyElementMatrices<true>(*em, xd, yd, w)
                  : ApplyElementMatrices<false>(*em, xd, yd, w);
  }
}

// Construction fixes the child pairing and clones active basis/restriction data.
// Ownership transfer requires finalized operators whose structure and outer scaling
// stay unchanged; values in the shared passive QData vectors may still be updated.
// A rejected packing attempt leaves ownership intact for the ordinary wrapper.
std::unique_ptr<ComplexWrapperOperator>
CreateComplexOperator(std::unique_ptr<palace::Operator> &&Ar,
                      std::unique_ptr<palace::Operator> &&Ai,
                      const FiniteElementSpace &fespace)
{
  auto packed = PackComplexOperator(Ar, Ai, fespace);
  if (packed)
  {
    return packed;
  }
  return std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

namespace
{

int CeedInternalCallocArray(size_t n, size_t unit, void *p)
{
  *(void **)p = calloc(n, unit);
  MFEM_ASSERT(!n || !unit || *(void **)p,
              "calloc failed to allocate " << n << " members of size " << unit << "!");
  return 0;
}

int CeedInternalFree(void *p)
{
  free(*(void **)p);
  *(void **)p = nullptr;
  return 0;
}

#define CeedInternalCalloc(n, p) CeedInternalCallocArray((n), sizeof(**(p)), p)

void CeedOperatorAssembleCOO(Ceed ceed, CeedOperator op, bool skip_zeros, CeedSize *nnz,
                             CeedInt **rows, CeedInt **cols, CeedVector *vals,
                             CeedMemType *mem)
{
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, mem));

  // Assemble sparsity pattern (rows, cols are always host pointers).
  PalaceCeedCall(ceed, CeedOperatorLinearAssembleSymbolic(op, nnz, rows, cols));

  // Assemble values.
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, vals));
  PalaceCeedCall(ceed, CeedOperatorLinearAssemble(op, *vals));

  // Filter out zero entries. For now, eliminating zeros happens all on the host.
  // std::cout << "  Operator full assembly (COO) has " << *nnz << " NNZ";
  if (skip_zeros && *nnz > 0)
  {
    // XX TODO: Use Thrust for this (thrust::copy_if and thrust::zip_iterator)
    CeedInt *new_rows, *new_cols;
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_rows));
    PalaceCeedCall(ceed, CeedInternalCalloc(*nnz, &new_cols));

    CeedVector new_vals;
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, *nnz, &new_vals));

    CeedSize q = 0;
    const CeedScalar *vals_array;
    CeedScalar *new_vals_array;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(*vals, CEED_MEM_HOST, &vals_array));
    PalaceCeedCall(ceed, CeedVectorGetArrayWrite(new_vals, CEED_MEM_HOST, &new_vals_array));
    for (CeedSize k = 0; k < *nnz; k++)
    {
      if (vals_array[k] != 0.0)
      {
        new_rows[q] = (*rows)[k];
        new_cols[q] = (*cols)[k];
        new_vals_array[q] = vals_array[k];
        q++;
      }
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(*vals, &vals_array));
    PalaceCeedCall(ceed, CeedVectorRestoreArray(new_vals, &new_vals_array));

    PalaceCeedCall(ceed, CeedInternalFree(rows));
    PalaceCeedCall(ceed, CeedInternalFree(cols));
    PalaceCeedCall(ceed, CeedVectorDestroy(vals));

    *nnz = q;
    *rows = new_rows;
    *cols = new_cols;
    *vals = new_vals;

    // std::cout << " (new NNZ after removal: " << *nnz << ")";
  }
  // std::cout << "\n";
}

std::unique_ptr<hypre::HypreCSRMatrix> OperatorCOOtoCSR(Ceed ceed, CeedInt m, CeedInt n,
                                                        CeedSize nnz, CeedInt *rows,
                                                        CeedInt *cols, CeedVector vals,
                                                        CeedMemType mem, bool set)
{
  // Preallocate CSR memory on host (like PETSc's MatSetValuesCOO). Check for overflow for
  // large nonzero counts.
  const int nnz_int = mfem::internal::to_int(nnz);
  mfem::Array<int> I(m + 1), J(nnz_int), perm(nnz_int), Jmap(nnz_int + 1);
  I = 0;
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(),
            [&](const int &i, const int &j) { return (rows[i] < rows[j]); });

  int q = -1;  // True nnz index
  for (int k = 0; k < nnz_int;)
  {
    // Sort column entries in the row.
    const int row = rows[perm[k]];
    const int start = k;
    while (k < nnz_int && rows[perm[k]] == row)
    {
      k++;
    }
    std::sort(perm.begin() + start, perm.begin() + k,
              [&](const int &i, const int &j) { return (cols[i] < cols[j]); });

    q++;
    I[row + 1] = 1;
    J[q] = cols[perm[start]];
    Jmap[q + 1] = 1;
    for (int p = start + 1; p < k; p++)
    {
      if (cols[perm[p]] != cols[perm[p - 1]])
      {
        // New nonzero.
        q++;
        I[row + 1]++;
        J[q] = cols[perm[p]];
        Jmap[q + 1] = 1;
      }
      else
      {
        Jmap[q + 1]++;
      }
    }
  }
  PalaceCeedCall(ceed, CeedInternalFree(&rows));
  PalaceCeedCall(ceed, CeedInternalFree(&cols));

  // Finalize I, Jmap.
  const int nnz_new = q + 1;
  I[0] = 0;
  for (int i = 0; i < m; i++)
  {
    I[i + 1] += I[i];
  }
  Jmap[0] = 0;
  for (int k = 0; k < nnz_new; k++)
  {
    Jmap[k + 1] += Jmap[k];
  }

  // Construct and fill the final CSR matrix. On GPU, MFEM and Hypre share the same memory
  // space. On CPU, the inner nested OpenMP loop (if enabled in MFEM) should be ignored.
  auto mat = std::make_unique<hypre::HypreCSRMatrix>(m, n, nnz_new);
  {
    const auto *d_I_old = I.Read();
    auto *d_I = mat->GetI();
    mfem::forall(m + 1, [=] MFEM_HOST_DEVICE(int i) { d_I[i] = d_I_old[i]; });
  }
  {
    const auto *d_J_old = J.Read();
    auto *d_J = mat->GetJ();
    mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k) { d_J[k] = d_J_old[k]; });
  }
  {
    auto FillValues = [&](const double *vals_array)
    {
      const auto *d_perm = perm.Read();
      const auto *d_Jmap = Jmap.Read();
      auto *d_A = mat->GetData();
      if (set)
      {
        mfem::forall(nnz_new, [=] MFEM_HOST_DEVICE(int k)
                     { d_A[k] = vals_array[d_perm[d_Jmap[k]]]; });
      }
      else
      {
        mfem::forall(nnz_new,
                     [=] MFEM_HOST_DEVICE(int k)
                     {
                       double sum = 0.0;
                       for (int p = d_Jmap[k]; p < d_Jmap[k + 1]; p++)
                       {
                         sum += vals_array[d_perm[p]];
                       }
                       d_A[k] = sum;
                     });
      }
    };
    Ceed ceed;
    const CeedScalar *vals_array;
    PalaceCeedCallBackend(CeedVectorGetCeed(vals, &ceed));
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals, mem, &vals_array));
    if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem != CEED_MEM_DEVICE)
    {
      // Copy values to device before filling.
      Vector d_vals(nnz_int);
      {
        auto *d_vals_array = d_vals.HostWrite();
        PalacePragmaOmp(parallel for schedule(static))
        for (int k = 0; k < nnz_int; k++)
        {
          d_vals_array[k] = vals_array[k];
        }
      }
      FillValues(d_vals.Read());
    }
    else
    {
      // No copy required.
      FillValues(vals_array);
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals, &vals_array));
    PalaceCeedCall(ceed, CeedVectorDestroy(&vals));
  }

  return mat;
}

}  // namespace

std::unique_ptr<hypre::HypreCSRMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                                bool skip_zeros, bool set)
{
  // Assemble operators on each thread.
  std::vector<std::unique_ptr<hypre::HypreCSRMatrix>> loc_mat(op.Size());
  PalacePragmaOmp(parallel if (op.Size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op.Size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op[id], &ceed));

    // Check if the operator is empty, otherwise assemble.
    CeedInt nsub_ops;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &nsub_ops));
    if (nsub_ops == 0)
    {
      loc_mat[id] = std::make_unique<hypre::HypreCSRMatrix>(op.Height(), op.Width(), 0);
    }
    else
    {
      // First, get matrix on master thread in COO format, with rows/cols always on host
      // and vals potentially on the device. Process skipping zeros if desired.
      CeedSize nnz;
      CeedInt *rows, *cols;
      CeedVector vals;
      CeedMemType mem;
      CeedOperatorAssembleCOO(ceed, op[id], skip_zeros, &nnz, &rows, &cols, &vals, &mem);
      PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));

      // Convert COO to CSR (on each thread). The COO memory is free'd internally.
      loc_mat[id] =
          OperatorCOOtoCSR(ceed, op.Height(), op.Width(), nnz, rows, cols, vals, mem, set);
    }
  }

  // Add CSR matrix objects from each thread (HYPRE's hypre_CSRMatrixAdd uses threads
  // internally as available). We have to scale the duplicated nonzeros when set = true.
  auto mat = std::move(loc_mat[0]);
  std::unique_ptr<hypre::HypreCSRMatrix> b_mat;
  if (set && op.Size() > 1)
  {
    b_mat = std::make_unique<hypre::HypreCSRMatrix>(hypre_CSRMatrixClone(*mat, 0));
    hypre_CSRMatrixSetConstantValues(*b_mat, 1.0);
    for (std::size_t id = 1; id < op.Size(); id++)
    {
      hypre_CSRMatrix *b_loc_mat = hypre_CSRMatrixClone(*loc_mat[id], 0);
      hypre_CSRMatrixSetConstantValues(b_loc_mat, 1.0);
      b_mat = std::make_unique<hypre::HypreCSRMatrix>(
          hypre_CSRMatrixAdd(1.0, *b_mat, 1.0, b_loc_mat));
      hypre_CSRMatrixDestroy(b_loc_mat);
    }
  }
  for (std::size_t id = 1; id < op.Size(); id++)
  {
    mat = std::make_unique<hypre::HypreCSRMatrix>(
        hypre_CSRMatrixAdd(1.0, *mat, 1.0, *loc_mat[id]));
  }
  if (set && op.Size() > 1)
  {
    const auto *d_b_data = b_mat->GetData();
    auto *d_data = mat->GetData();
    mfem::forall(mat->NNZ(),
                 [=] MFEM_HOST_DEVICE(int i) { d_data[i] *= 1.0 / d_b_data[i]; });
  }

  return mat;
}

std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse)
{
  auto SingleOperatorCoarsen =
      [&fespace_coarse](Ceed ceed, CeedOperator op_fine, CeedOperator *op_coarse)
  {
    CeedBasis basis_fine;
    CeedElemTopology geom;
    PalaceCeedCall(ceed, CeedOperatorGetActiveBasis(op_fine, &basis_fine));
    PalaceCeedCall(ceed, CeedBasisGetTopology(basis_fine, &geom));

    // The coarse operator reuses the fine operator's quadrature data, so its finite element
    // restriction must consume the identical element subset in the same mesh order.
    const auto mfem_geom = GetMfemTopology(geom);
    const auto &geom_data =
        fespace_coarse.GetMesh().GetCeedGeomFactorData(ceed).at(mfem_geom);
    CeedElemRestriction restr_fine;
    CeedInt num_elem;
    PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestriction(op_fine, &restr_fine));
    PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(restr_fine, &num_elem));
    CeedElementSubset subset = CeedElementSubset::Full;
    if (num_elem != static_cast<CeedInt>(geom_data.indices.size()))
    {
      MFEM_VERIFY(num_elem == static_cast<CeedInt>(geom_data.active_indices.size()),
                  "Fine operator uses an unknown element subset in CeedOperatorCoarsen!");
      subset = CeedElementSubset::Active;
    }
    CeedElemRestriction restr_coarse =
        fespace_coarse.GetCeedElemRestriction(ceed, mfem_geom, geom_data, subset);
    CeedBasis basis_coarse = fespace_coarse.GetCeedBasis(ceed, GetMfemTopology(geom));

    PalaceCeedCall(ceed, CeedOperatorMultigridLevelCreate(op_fine, nullptr, restr_coarse,
                                                          basis_coarse, op_coarse, nullptr,
                                                          nullptr));
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(*op_coarse));
  };

  // Initialize the coarse operator.
  auto op_coarse = std::make_unique<SymmetricOperator>(fespace_coarse.GetVSize(),
                                                       fespace_coarse.GetVSize());

  // Assemble the coarse operator by coarsening each sub-operator (over threads, geometry
  // types, integrators) of the original fine operator.
  PalacePragmaOmp(parallel if (op_fine.Size() > 1))
  {
    const int id = utils::GetThreadNum();
    MFEM_ASSERT(static_cast<std::size_t>(id) < op_fine.Size(),
                "Out of bounds access for thread number " << id << "!");
    Ceed ceed;
    PalaceCeedCallBackend(CeedOperatorGetCeed(op_fine[id], &ceed));
    {
      Ceed ceed_parent;
      PalaceCeedCall(ceed, CeedGetParent(ceed, &ceed_parent));
      if (ceed_parent)
      {
        ceed = ceed_parent;
      }
    }
    CeedInt nsub_ops_fine;
    CeedOperator *sub_ops_fine;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op_fine[id], &nsub_ops_fine));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op_fine[id], &sub_ops_fine));
    for (CeedInt k = 0; k < nsub_ops_fine; k++)
    {
      CeedOperator sub_op_coarse;
      SingleOperatorCoarsen(ceed, sub_ops_fine[k], &sub_op_coarse);
      // Sub-operator owned by ceed::Operator.
      op_coarse->AddSubOperator(sub_op_coarse);
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  op_coarse->Finalize();

  return op_coarse;
}

}  // namespace palace::ceed
