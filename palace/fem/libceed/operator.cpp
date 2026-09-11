// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <cstring>
#include <limits>
#include <numeric>
#include <typeinfo>
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

Operator::~Operator()
{
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
  application_revision.fetch_add(1, std::memory_order_relaxed);
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
  CeedMemType mem;
  MFEM_VERIFY(diag.Size() == height, "Invalid size for diagonal vector!");
  diag = 0.0;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    mem = CEED_MEM_HOST;
  }
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
  CeedMemType mem;
  PalaceCeedCallBackend(CeedOperatorGetCeed(op[0], &ceed));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    mem = CEED_MEM_HOST;
  }
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

// Reference-counted metadata must be destroyed even when a match is rejected.
// Keep ownership local while inspecting each operator.
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
  // Selection relies on the names and field layouts of Palace's assembled apply kernels.
  bool at_points, composite;
  PalaceCeedCall(ceed, CeedOperatorIsComposite(op, &composite));
  PalaceCeedCall(ceed, CeedOperatorIsAtPoints(op, &at_points));
  if (composite || at_points || !SameContext(ceed, CeedOperatorReturnCeed(op)))
  {
    return {};
  }
  ScopedQFunction qf;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(op, &qf.value));
  const char *name;
  PalaceCeedCall(ceed, CeedQFunctionGetKernelName(qf.value, &name));
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
  // Validate all active field modes/sizes and the actual trial/test handles. Names
  // alone cannot distinguish mass from diffusion, mixed spaces, or boundary operators.
  for (CeedInt j = 0; j < no; j++)
  {
    const auto mode = j == 0 ? CEED_EVAL_INTERP : CEED_EVAL_CURL;
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
  bool tensor;
  CeedFESpace space;
  CeedInt dim, components, qcomponents, rcomponents, rnodes, block;
  CeedSize length, in_length, out_length;
  CeedRestrictionType type;
  PalaceCeedCall(ceed, CeedBasisIsTensor(basis, &tensor));
  PalaceCeedCall(ceed, CeedBasisGetDimension(basis, &dim));
  PalaceCeedCall(ceed, CeedBasisGetFESpace(basis, &space));
  PalaceCeedCall(ceed, CeedBasisGetNumComponents(basis, &components));
  PalaceCeedCall(ceed, CeedBasisGetNumNodes(basis, &leaf->nodes));
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(basis, &leaf->qpts));
  PalaceCeedCall(
      ceed, CeedBasisGetNumQuadratureComponents(basis, CEED_EVAL_INTERP, &qcomponents));
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(restriction, &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(restriction, &rcomponents));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(restriction, &rnodes));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(restriction, &leaf->elements));
  PalaceCeedCall(ceed, CeedElemRestrictionGetBlockSize(restriction, &block));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(restriction, &length));
  PalaceCeedCall(ceed, CeedOperatorGetActiveVectorLengths(op, &in_length, &out_length));
  if (tensor || dim != 3 || components != 1 || qcomponents != 3 || rcomponents != 1 ||
      (space != CEED_FE_SPACE_HCURL && space != CEED_FE_SPACE_HDIV) ||
      (curl && space != CEED_FE_SPACE_HCURL) || rnodes != leaf->nodes || block != 1 ||
      leaf->elements <= 0 || leaf->qpts <= 0 || length != size || in_length != size ||
      out_length != size ||
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
  CeedInt qe, qn, qc;
  CeedSize qlength, vlength;
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(qr, &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(qr, &qe));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(qr, &qn));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(qr, &qc));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(qr, &qlength));
  PalaceCeedCall(ceed, CeedVectorGetLength(qv, &vlength));
  if (type != CEED_RESTRICTION_STRIDED || qe != leaf->elements || qn != leaf->qpts ||
      qc != (curl ? 18 : 9) || qlength != static_cast<CeedSize>(qe) * qn * qc ||
      vlength != qlength)
  {
    return {};
  }
  return leaf;
}

bool CollectPackedLeaves(Ceed ceed, CeedOperator op, std::vector<CeedOperator> &leaves)
{
  bool composite;
  PalaceCeedCall(ceed, CeedOperatorIsComposite(op, &composite));
  if (!composite)
  {
    if (leaves.size() == CEED_COMPOSITE_MAX)
    {
      return false;
    }
    leaves.push_back(op);  // Borrowed until the original wrapper takes ownership.
    return true;
  }
  CeedInt count;
  CeedOperator *sub;
  PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op, &count));
  PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op, &sub));
  for (CeedInt j = 0; j < count; j++)
  {
    if (!CollectPackedLeaves(ceed, sub[j], leaves))
    {
      return false;
    }
  }
  return true;
}

void ClonePackedRestriction(Ceed ceed, const PackedLeaf &leaf, int size,
                            CeedElemRestriction *packed)
{
  auto original = leaf.active.restriction.value;
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
    PalaceCeedCall(ceed,
                   CeedElemRestrictionCreateCurlOriented(
                       ceed, leaf.elements, leaf.nodes, 2, size, length, CEED_MEM_HOST,
                       CEED_COPY_VALUES, offsets, orientations, packed));
    PalaceCeedCall(ceed,
                   CeedElemRestrictionRestoreCurlOrientations(original, &orientations));
  }
  else if (type == CEED_RESTRICTION_ORIENTED)
  {
    const bool *orientations;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetOrientations(original, CEED_MEM_HOST, &orientations));
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(ceed, leaf.elements, leaf.nodes,
                                                           2, size, length, CEED_MEM_HOST,
                                                           CEED_COPY_VALUES, offsets,
                                                           orientations, packed));
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreOrientations(original, &orientations));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, leaf.elements, leaf.nodes, 2, size,
                                                   length, CEED_MEM_HOST, CEED_COPY_VALUES,
                                                   offsets, packed));
  }
  PalaceCeedCall(ceed, CeedElemRestrictionRestoreOffsets(original, &offsets));
}

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
  else
  {
    PalaceCeedCall(ceed, CeedBasisGetDiv(original, &derivative));
    PalaceCeedCall(ceed, CeedBasisCreateHdiv(ceed, topology, 2, leaf.nodes, leaf.qpts,
                                             interp, derivative, qref, weights, packed));
  }
}

void AddPackedPair(Ceed ceed, const PackedLeaf &real, const PackedLeaf &imag, int size,
                   CeedOperator composite)
{
  ScopedRestriction restriction;
  ScopedBasis basis;
  ScopedQFunction qf;
  ScopedOperator op;
  ClonePackedRestriction(ceed, real, size, &restriction.value);
  ClonePackedBasis(ceed, real, &basis.value);
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                           ceed, 1, real.curl ? f_apply_complex_33 : f_apply_complex_3,
                           real.curl ? PalaceQFunctionRelativePath(f_apply_complex_33_loc)
                                     : PalaceQFunctionRelativePath(f_apply_complex_3_loc),
                           &qf.value));
  PalaceCeedCall(ceed,
                 CeedQFunctionAddInput(qf.value, "qr", real.curl ? 18 : 9, CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "qi", 9, CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "u", 6, CEED_EVAL_INTERP));
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "v", 6, CEED_EVAL_INTERP));
  if (real.curl)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "curl_u", 6, CEED_EVAL_CURL));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "curl_v", 6, CEED_EVAL_CURL));
  }
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &op.value));
  PalaceCeedCall(ceed, CeedOperatorSetField(op.value, "qr", real.qdata.restriction.value,
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
  PalaceCeedCall(ceed, CeedOperatorCompositeAddSub(composite, op.value));
}

}  // namespace

struct PackedComplexOperator::Data
{
  ScopedContext ceed;
  ScopedOperator packed;
  ScopedVector x, y;
  std::unique_ptr<ComplexWrapperOperator> remainder;
  mutable ComplexVector temp;
  std::size_t pairs = 0, remainder_terms = 0;
  const Operator *source[2] = {};
  std::size_t revision[2] = {};
};

std::unique_ptr<PackedComplexOperator::Data>
PackedComplexOperator::Build(const palace::Operator *Ar, const palace::Operator *Ai)
{
  const auto *real = dynamic_cast<const Operator *>(Ar);
  const auto *imag = dynamic_cast<const Operator *>(Ai);
  if (!real || !imag || internal::NumCeeds() != 1 || utils::GetMaxThreads() > 1 ||
      utils::InParallel() || mfem::Device::Allows(mfem::Backend::DEVICE_MASK) ||
      real->Size() != 1 || imag->Size() != 1 || real->HasDofMultiplicity() ||
      imag->HasDofMultiplicity() || real->Width() != real->Height() ||
      real->Width() != imag->Width() || real->Height() != imag->Height() ||
      real->Width() <= 0 || real->Width() > std::numeric_limits<int>::max() / 2)
  {
    return {};
  }
  // Do not bypass an unknown subclass's application semantics.
  if ((typeid(*real) != typeid(Operator) && typeid(*real) != typeid(SymmetricOperator)) ||
      (typeid(*imag) != typeid(Operator) && typeid(*imag) != typeid(SymmetricOperator)))
  {
    return {};
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  CeedMemType mem;
  const char *resource;
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  if (mem != CEED_MEM_HOST || std::strncmp(resource, "/cpu/", 5) != 0 ||
      !SameContext(ceed, CeedOperatorReturnCeed((*real)[0])) ||
      !SameContext(ceed, CeedOperatorReturnCeed((*imag)[0])))
  {
    return {};
  }
  const int size = real->Width();
  std::vector<CeedOperator> leaves[2];
  if (!CollectPackedLeaves(ceed, (*real)[0], leaves[0]) ||
      !CollectPackedLeaves(ceed, (*imag)[0], leaves[1]))
  {
    return {};  // Do not exceed the cap when flattening either packed or remainder terms.
  }
  std::vector<std::unique_ptr<PackedLeaf>> fields[2];
  std::vector<bool> matched[2];
  for (int part = 0; part < 2; part++)
  {
    matched[part].resize(leaves[part].size(), false);
    for (auto leaf : leaves[part])
    {
      fields[part].push_back(InspectPackedLeaf(ceed, leaf, size));
    }
  }
  std::unique_ptr<Data> data;
  for (std::size_t r = 0; r < fields[0].size(); r++)
  {
    if (!fields[0][r])
    {
      continue;
    }
    for (std::size_t i = 0; i < fields[1].size(); i++)
    {
      if (matched[1][i] || !fields[1][i] || fields[1][i]->curl ||
          fields[0][r]->active.basis.value != fields[1][i]->active.basis.value ||
          fields[0][r]->active.restriction.value != fields[1][i]->active.restriction.value)
      {
        continue;
      }
      if (!data)
      {
        data = std::make_unique<Data>();
        PalaceCeedCall(ceed, CeedReferenceCopy(ceed, &data->ceed.value));
        PalaceCeedCall(ceed, CeedOperatorCreateComposite(ceed, &data->packed.value));
      }
      AddPackedPair(ceed, *fields[0][r], *fields[1][i], size, data->packed.value);
      matched[0][r] = matched[1][i] = true;
      data->pairs++;
      break;
    }
  }
  if (!data)
  {
    return {};  // This also handles ranks with no local volume terms.
  }
  // Keep every unmatched leaf exactly once. In particular, differing real/imaginary
  // boundary terms must not disable fusion of compatible domain contributions.
  std::unique_ptr<palace::Operator> remainder_parts[2];
  for (int part = 0; part < 2; part++)
  {
    std::unique_ptr<Operator> remainder;
    for (std::size_t j = 0; j < leaves[part].size(); j++)
    {
      if (!matched[part][j])
      {
        if (!remainder)
        {
          remainder = std::make_unique<Operator>(size, size);
        }
        CeedOperator copy = nullptr;
        PalaceCeedCall(ceed, CeedOperatorReferenceCopy(leaves[part][j], &copy));
        remainder->AddSubOperator(copy);
        data->remainder_terms++;
      }
    }
    if (remainder)
    {
      remainder->Finalize();
      remainder_parts[part] = std::move(remainder);
    }
  }
  if (remainder_parts[0] || remainder_parts[1])
  {
    data->remainder = std::make_unique<ComplexWrapperOperator>(
        std::move(remainder_parts[0]), std::move(remainder_parts[1]));
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(data->packed.value));
  PalaceCeedCall(ceed,
                 CeedVectorCreate(ceed, 2 * static_cast<CeedSize>(size), &data->x.value));
  PalaceCeedCall(ceed,
                 CeedVectorCreate(ceed, 2 * static_cast<CeedSize>(size), &data->y.value));
  data->source[0] = real;
  data->source[1] = imag;
  data->revision[0] = real->ApplicationRevision();
  data->revision[1] = imag->ApplicationRevision();
  return data;
}

PackedComplexOperator::PackedComplexOperator(std::unique_ptr<palace::Operator> &&Ar,
                                             std::unique_ptr<palace::Operator> &&Ai,
                                             std::unique_ptr<Data> &&data)
  : ComplexWrapperOperator(std::move(Ar), std::move(Ai)), data(std::move(data))
{
}

PackedComplexOperator::PackedComplexOperator(const palace::Operator *Ar,
                                             const palace::Operator *Ai,
                                             std::unique_ptr<Data> &&data)
  : ComplexWrapperOperator(Ar, Ai), data(std::move(data))
{
}

PackedComplexOperator::~PackedComplexOperator() = default;

bool PackedComplexOperator::CanApplyPacked() const
{
  return data->revision[0] == data->source[0]->ApplicationRevision() &&
         data->revision[1] == data->source[1]->ApplicationRevision();
}

std::size_t PackedComplexOperator::NumFusedPairs() const
{
  return CanApplyPacked() ? data->pairs : 0;
}

std::size_t PackedComplexOperator::NumRemainderTerms() const
{
  return data->remainder_terms;
}

void PackedComplexOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  MFEM_ASSERT(x.Size() == width && y.Size() == height,
              "Invalid dimensions for PackedComplexOperator::Mult!");
  if (!CanApplyPacked())
  {
    ComplexWrapperOperator::Mult(x, y);
    return;
  }
  Ceed ceed = data->ceed.value;
  CeedScalar *packed_x;
  const std::size_t bytes = static_cast<std::size_t>(width) * sizeof(CeedScalar);
  PalaceCeedCall(ceed, CeedVectorGetArrayWrite(data->x.value, CEED_MEM_HOST, &packed_x));
  std::memcpy(packed_x, x.Real().HostRead(), bytes);
  std::memcpy(packed_x + width, x.Imag().HostRead(), bytes);
  PalaceCeedCall(ceed, CeedVectorRestoreArray(data->x.value, &packed_x));
  PalaceCeedCall(ceed, CeedVectorSetValue(data->y.value, 0.0));
  PalaceCeedCall(ceed, CeedOperatorApplyAdd(data->packed.value, data->x.value,
                                            data->y.value, CEED_REQUEST_IMMEDIATE));
  if (data->remainder)
  {
    data->remainder->Mult(x, y);
  }
  const CeedScalar *packed_y;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(data->y.value, CEED_MEM_HOST, &packed_y));
  if (data->remainder)
  {
    auto *yr = y.Real().HostReadWrite();
    auto *yi = y.Imag().HostReadWrite();
    for (int j = 0; j < height; j++)
    {
      yr[j] += packed_y[j];
      yi[j] += packed_y[j + height];
    }
  }
  else
  {
    std::memcpy(y.Real().HostWrite(), packed_y, bytes);
    std::memcpy(y.Imag().HostWrite(), packed_y + height, bytes);
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(data->y.value, &packed_y));
}

void PackedComplexOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                    std::complex<double> a) const
{
  if (a != std::complex<double>{0.0})
  {
    // Stay on the packed path, also for non-unit complex coefficients. The real CEED
    // operators only support AddMult with coefficient one.
    data->temp.SetSize(height);
    Mult(x, data->temp);
    y.AXPY(a, data->temp);
  }
}

std::unique_ptr<ComplexWrapperOperator>
CreateComplexOperator(std::unique_ptr<palace::Operator> &&Ar,
                      std::unique_ptr<palace::Operator> &&Ai)
{
  auto data = PackedComplexOperator::Build(Ar.get(), Ai.get());
  if (data)
  {
    return std::unique_ptr<ComplexWrapperOperator>(
        new PackedComplexOperator(std::move(Ar), std::move(Ai), std::move(data)));
  }
  return std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

std::unique_ptr<ComplexWrapperOperator> CreateComplexOperator(const palace::Operator *Ar,
                                                              const palace::Operator *Ai)
{
  auto data = PackedComplexOperator::Build(Ar, Ai);
  if (data)
  {
    return std::unique_ptr<ComplexWrapperOperator>(
        new PackedComplexOperator(Ar, Ai, std::move(data)));
  }
  return std::make_unique<ComplexWrapperOperator>(Ar, Ai);
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
  for (int k = 0; k < nnz_int; k++)
  {
    perm[k] = k;
  }
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

    const auto &geom_data =
        fespace_coarse.GetMesh().GetCeedGeomFactorData(ceed).at(GetMfemTopology(geom));
    CeedElemRestriction restr_coarse = fespace_coarse.GetCeedElemRestriction(
        ceed, GetMfemTopology(geom), geom_data.indices);
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
      op_coarse->AddSubOperator(sub_op_coarse);  // Sub-operator owned by ceed::Operator
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  op_coarse->Finalize();

  return op_coarse;
}

}  // namespace palace::ceed
