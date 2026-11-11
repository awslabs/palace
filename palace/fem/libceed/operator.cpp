// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "operator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
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

  // Whether the QData is the same at all points of each element, and for curl leaves the
  // number of leading points with nonzero curl interpolation rows.
  bool element_qdata = false;
  CeedInt curl_points = 0;

  // Whether the QData is stored once per element and interpolated to the points by a basis
  // with one node (see ceed::UsePerElementQuadratureData), and for element QData the factor
  // of each point (the interpolation values of that basis, or one).
  bool stored_per_element = false;
  std::vector<CeedScalar> point_scale;
};

bool FieldMatches(Ceed ceed, CeedQFunctionField field, CeedInt size, CeedEvalMode mode)
{
  CeedInt actual_size;
  CeedEvalMode actual_mode;
  PalaceCeedCall(ceed,
                 CeedQFunctionFieldGetData(field, nullptr, &actual_size, &actual_mode));
  return size == actual_size && mode == actual_mode;
}

// Bases with unit weights for affine elements (see ceed::InitReducedBasis) give QData
// which is the same at all points of an element, since the geometry factors and
// coefficients are. Such QData is read once per element instead of once per point.
bool HasElementQData(Ceed ceed, const PackedLeaf &leaf, CeedInt num_comp)
{
  bool strided;
  CeedInt layout[3];
  PalaceCeedCall(ceed,
                 CeedElemRestrictionIsStrided(leaf.qdata.restriction.value, &strided));
  if (!strided || leaf.qpts < 2)
  {
    return false;
  }
  PalaceCeedCall(ceed, CeedElemRestrictionGetLLayout(leaf.qdata.restriction.value, layout));
  const CeedScalar *qd;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(leaf.qdata.vector.value, CEED_MEM_HOST, &qd));
  bool constant = true;
  for (CeedInt e = 0; constant && e < leaf.elements; e++)
  {
    for (CeedInt c = 0; constant && c < num_comp; c++)
    {
      const CeedSize offset =
          static_cast<CeedSize>(e) * layout[2] + static_cast<CeedSize>(c) * layout[1];
      for (CeedInt q = 1; q < leaf.qpts; q++)
      {
        if (qd[offset + static_cast<CeedSize>(q) * layout[0]] != qd[offset])
        {
          constant = false;
          break;
        }
      }
    }
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(leaf.qdata.vector.value, &qd));
  return constant;
}

// Number of leading points at which the curl interpolation of a basis has nonzero rows
// (for each component). For the hierarchical bases of ceed::InitReducedBasis, the curl is
// represented on fewer points than the values.
CeedInt CurlPoints(Ceed ceed, const PackedLeaf &leaf)
{
  const CeedScalar *curl;
  PalaceCeedCall(ceed, CeedBasisGetCurl(leaf.active.basis.value, &curl));
  CeedInt num_points = 0;
  for (CeedInt d = 0; d < 3; d++)
  {
    for (CeedInt q = num_points; q < leaf.qpts; q++)
    {
      const CeedScalar *row =
          curl + (static_cast<std::size_t>(d) * leaf.qpts + q) * leaf.nodes;
      if (std::any_of(row, row + leaf.nodes, [](CeedScalar b) { return b != 0.0; }))
      {
        num_points = q + 1;
      }
    }
  }
  return num_points;
}

// Element restriction reading the QData of each (subset) element at its first point, with
// each node selecting a group of components starting at the given offset.
void CreateElementQDataRestriction(Ceed ceed, CeedElemRestriction full,
                                   const std::vector<int> *positions, CeedInt num_elem,
                                   const std::vector<CeedInt> &comp_offsets,
                                   CeedInt num_comp, CeedElemRestriction *restr)
{
  CeedInt layout[3];
  CeedSize length;
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(full, &length));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLLayout(full, layout));
  const CeedInt num_nodes = static_cast<CeedInt>(comp_offsets.size());
  std::vector<CeedInt> offsets(static_cast<std::size_t>(num_elem) * num_nodes);
  for (CeedInt k = 0; k < num_elem; k++)
  {
    const CeedInt e = positions ? (*positions)[k] : k;
    for (CeedInt j = 0; j < num_nodes; j++)
    {
      offsets[k * num_nodes + j] = e * layout[2] + comp_offsets[j] * layout[1];
    }
  }
  PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, num_elem, num_nodes, num_comp,
                                                 layout[1], length, CEED_MEM_HOST,
                                                 CEED_COPY_VALUES, offsets.data(), restr));
}

// Basis copying the values at its nodes to the points: point q gets the value of node
// node_of_point[q] times point_scale[q], or zero if the node is negative.
void CreateSelectorBasis(Ceed ceed, CeedInt num_comp, CeedInt num_nodes,
                         const std::vector<CeedInt> &node_of_point,
                         const std::vector<CeedScalar> &point_scale, CeedBasis *basis)
{
  const CeedInt num_qpts = static_cast<CeedInt>(node_of_point.size());
  MFEM_ASSERT(point_scale.size() == node_of_point.size(),
              "Invalid point scaling for a selector basis!");
  std::vector<CeedScalar> interp(static_cast<std::size_t>(num_qpts) * num_nodes, 0.0),
      grad(interp.size(), 0.0), q_ref(num_qpts, 0.0), q_weight(num_qpts, 1.0);
  for (CeedInt q = 0; q < num_qpts; q++)
  {
    if (node_of_point[q] >= 0)
    {
      interp[q * num_nodes + node_of_point[q]] = point_scale[q];
    }
  }
  PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, num_comp, num_nodes,
                                         num_qpts, interp.data(), grad.data(), q_ref.data(),
                                         q_weight.data(), basis));
}

// Factors of the points of a selector basis reading a leaf's element QData: those of the
// leaf's points, followed by those of its leading curl points for a stacked basis.
std::vector<CeedScalar> SelectorScale(const PackedLeaf &leaf, CeedInt num_qpts)
{
  std::vector<CeedScalar> scale(num_qpts);
  for (CeedInt q = 0; q < num_qpts; q++)
  {
    scale[q] = leaf.point_scale[(q < leaf.qpts) ? q : q - leaf.qpts];
  }
  return scale;
}

// For a curl leaf with element QData and fewer curl points than points, one interpolation
// gives the values at all points followed by the curl at the curl points (two copies for
// the real and imaginary coefficient blocks, see ClonePackedBasis). The mass and curl-curl
// tensors then act on the respective points, with the same result as the original pair of
// interpolations, since the omitted curl rows are zero and the QData is constant.
void CreateStackedBasis(Ceed ceed, const PackedLeaf &leaf, CeedBasis *stacked)
{
  CeedElemTopology topology;
  const CeedScalar *interp, *curl;
  auto original = leaf.active.basis.value;
  PalaceCeedCall(ceed, CeedBasisGetTopology(original, &topology));
  PalaceCeedCall(ceed, CeedBasisGetInterp(original, &interp));
  PalaceCeedCall(ceed, CeedBasisGetCurl(original, &curl));
  const CeedInt Q = leaf.qpts, Qs = leaf.qpts + leaf.curl_points, P = leaf.nodes;
  std::vector<CeedScalar> s_interp(static_cast<std::size_t>(3) * Qs * P),
      s_curl(s_interp.size(), 0.0), q_ref(static_cast<std::size_t>(3) * Qs, 0.0),
      q_weight(Qs, 1.0);
  for (CeedInt d = 0; d < 3; d++)
  {
    std::copy_n(interp + static_cast<std::size_t>(d) * Q * P, Q * P,
                s_interp.begin() + static_cast<std::size_t>(d) * Qs * P);
    std::copy_n(curl + static_cast<std::size_t>(d) * Q * P, leaf.curl_points * P,
                s_interp.begin() + (static_cast<std::size_t>(d) * Qs + Q) * P);
  }
  PalaceCeedCall(ceed, CeedBasisCreateHcurl(ceed, topology, 2, P, Qs, s_interp.data(),
                                            s_curl.data(), q_ref.data(), q_weight.data(),
                                            stacked));
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
  const CeedInt q_size = curl ? 18 : 9;
  if (ni != (curl ? 3 : 2) || no != (curl ? 2 : 1) || qni != ni || qno != no ||
      (!FieldMatches(ceed, qinputs[0], q_size, CEED_EVAL_NONE) &&
       !FieldMatches(ceed, qinputs[0], q_size, CEED_EVAL_INTERP)))
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
  // gradients of a scalar field and H(curl) curl-curl (without mass) on the curl.
  {
    CeedFESpace space;
    PalaceCeedCall(ceed, CeedBasisGetFESpace(basis, &space));
    leaf->mode = (space == CEED_FE_SPACE_H1) ? CEED_EVAL_GRAD : CEED_EVAL_INTERP;
    CeedEvalMode mode;
    PalaceCeedCall(ceed, CeedQFunctionFieldGetEvalMode(qinputs[1], &mode));
    if (space == CEED_FE_SPACE_HCURL && !curl && mode == CEED_EVAL_CURL)
    {
      leaf->mode = CEED_EVAL_CURL;
    }
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
  auto qb = leaf->qdata.basis.value;
  if (!qr || qr == CEED_ELEMRESTRICTION_NONE || !qv || qv == CEED_VECTOR_ACTIVE ||
      qv == CEED_VECTOR_NONE || !qb ||
      !SameContext(ceed, CeedElemRestrictionReturnCeed(qr)) ||
      !SameContext(ceed, CeedVectorReturnCeed(qv)))
  {
    return {};
  }
  // Finalization of the original operator is a caller precondition and validates
  // basis/restriction dimensions and passive QData storage. Reuse the QData vector
  // and its restriction, so in-place value updates remain visible to both operators
  // (element QData stays constant on each element under coefficient updates).
  if (qb != CEED_BASIS_NONE)
  {
    // QData stored once per element (strided, one node per element), interpolated to the
    // points by a basis with one node whose values are the factors of the points.
    bool strided, q_tensor;
    CeedInt q_nodes, q_qpts, q_comp, r_nodes, r_comp;
    PalaceCeedCall(ceed, CeedElemRestrictionIsStrided(qr, &strided));
    PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(qr, &r_nodes));
    PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(qr, &r_comp));
    PalaceCeedCall(ceed, CeedBasisIsTensor(qb, &q_tensor));
    PalaceCeedCall(ceed, CeedBasisGetNumNodes(qb, &q_nodes));
    PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(qb, &q_qpts));
    PalaceCeedCall(ceed, CeedBasisGetNumComponents(qb, &q_comp));
    if (!FieldMatches(ceed, qinputs[0], q_size, CEED_EVAL_INTERP) ||
        !SameContext(ceed, CeedBasisReturnCeed(qb)) || !strided || q_tensor ||
        q_nodes != 1 || q_qpts != leaf->qpts || q_comp != q_size || r_nodes != 1 ||
        r_comp != q_size)
    {
      return {};
    }
    const CeedScalar *interp;
    PalaceCeedCall(ceed, CeedBasisGetInterp(qb, &interp));
    leaf->stored_per_element = leaf->element_qdata = true;
    leaf->point_scale.assign(interp, interp + leaf->qpts);
  }
  else
  {
    if (!FieldMatches(ceed, qinputs[0], q_size, CEED_EVAL_NONE))
    {
      return {};
    }
    leaf->element_qdata = HasElementQData(ceed, *leaf, q_size);
    leaf->point_scale.assign(leaf->qpts, 1.0);
  }
  leaf->curl_points = curl ? CurlPoints(ceed, *leaf) : 0;
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
  // The subset views of the real leaf's QData are derived from its strided layout.
  bool strided;
  PalaceCeedCall(ceed,
                 CeedElemRestrictionIsStrided(real.qdata.restriction.value, &strided));
  if (!strided)
  {
    return {};
  }
  // The leaves may also be assembled over a part of the element list (see
  // CeedGeomFactorData), with its own active subset.
  std::vector<const CeedGeomFactorData *> candidates = {&it->second};
  for (const auto &part : it->second.parts)
  {
    candidates.push_back(&part);
  }
  for (const auto *candidate : candidates)
  {
    const auto &data = *candidate;
    if (data.indices.size() != static_cast<std::size_t>(real.elements) ||
        data.active_indices.size() != static_cast<std::size_t>(imag.elements) ||
        data.active_indices.empty() || data.active_indices.size() >= data.indices.size())
    {
      continue;
    }
    if (real.active.restriction.value ==
            fespace.GetCeedElemRestriction(ceed, geom, data, CeedElementSubset::Full) &&
        imag.active.restriction.value ==
            fespace.GetCeedElemRestriction(ceed, geom, data, CeedElementSubset::Active))
    {
      return {&data, geom};
    }
  }
  return {};
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
  ScopedRestriction restriction, qdata, qdata_imag;
  ScopedBasis basis, qdata_basis, qdata_imag_basis;
  ScopedQFunction qf;
  ScopedOperator op;
  // The pair is applied over the imaginary leaf's elements: either the same element list as
  // the real leaf, or the active subset of it. Both leaves share the basis, so the node
  // count matches, and the cloned orientation maps coincide with the real leaf's on the
  // subset elements because both restrictions derive from the same space's element list.
  ClonePackedRestriction(ceed, imag.active.restriction.value, imag.elements, imag.nodes,
                         size, &restriction.value);
  // With element QData for both leaves, it is read once per element, and a curl leaf with
  // fewer curl points than points interpolates its values and curls together.
  const bool element_qdata = real.element_qdata && imag.element_qdata;
  const bool stacked = element_qdata && real.curl && real.curl_points < real.qpts;
  const bool curl = real.curl && !stacked;
  const CeedInt qr_size = curl ? 18 : 9;
  if (stacked)
  {
    CreateStackedBasis(ceed, real, &basis.value);
  }
  else
  {
    ClonePackedBasis(ceed, real, &basis.value);
  }
  if (element_qdata)
  {
    // The real coefficients of the subset elements, in place in the real leaf's QData. The
    // imaginary (mass) tensor does not act at the curl points of a stacked pair.
    const CeedInt num_qpts = stacked ? real.qpts + real.curl_points : real.qpts;
    std::vector<CeedInt> real_nodes(num_qpts, 0), imag_nodes(num_qpts, 0);
    for (CeedInt q = real.qpts; q < num_qpts; q++)
    {
      real_nodes[q] = 1;
      imag_nodes[q] = -1;
    }
    CreateElementQDataRestriction(
        ceed, real.qdata.restriction.value,
        subset.data ? &subset.data->active_indices : nullptr, imag.elements,
        stacked ? std::vector<CeedInt>{0, 9} : std::vector<CeedInt>{0}, qr_size,
        &qdata.value);
    CreateElementQDataRestriction(ceed, imag.qdata.restriction.value, nullptr,
                                  imag.elements, {0}, 9, &qdata_imag.value);
    CreateSelectorBasis(ceed, qr_size, stacked ? 2 : 1, real_nodes,
                        SelectorScale(real, num_qpts), &qdata_basis.value);
    CreateSelectorBasis(ceed, 9, 1, imag_nodes, SelectorScale(imag, num_qpts),
                        &qdata_imag_basis.value);
  }
  else if (subset.data)
  {
    // The real coefficients of the subset elements, in place in the real leaf's QData.
    CreateSubsetQDataRestriction(ceed, real.qdata.restriction.value,
                                 subset.data->active_indices, &qdata.value);
  }
  // Both kernels implement complex mass. Only f_apply_complex_33 adds curl-curl,
  // with a real tensor acting independently on Re(x) and Im(x). The pairing loop
  // must therefore leave imaginary curl leaves in the remainder.
  PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                           ceed, 1, curl ? f_apply_complex_33 : f_apply_complex_3,
                           curl ? PalaceQFunctionRelativePath(f_apply_complex_33_loc)
                                : PalaceQFunctionRelativePath(f_apply_complex_3_loc),
                           &qf.value));
  const CeedEvalMode qdata_mode = element_qdata ? CEED_EVAL_INTERP : CEED_EVAL_NONE;
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "qr", qr_size, qdata_mode));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "qi", 9, qdata_mode));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "u", 6, real.mode));
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "v", 6, real.mode));
  if (curl)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "curl_u", 6, CEED_EVAL_CURL));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "curl_v", 6, CEED_EVAL_CURL));
  }
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &op.value));
  PalaceCeedCall(
      ceed, CeedOperatorSetField(op.value, "qr",
                                 qdata.value ? qdata.value : real.qdata.restriction.value,
                                 element_qdata ? qdata_basis.value : CEED_BASIS_NONE,
                                 real.qdata.vector.value));
  PalaceCeedCall(ceed, CeedOperatorSetField(
                           op.value, "qi",
                           element_qdata ? qdata_imag.value : imag.qdata.restriction.value,
                           element_qdata ? qdata_imag_basis.value : CEED_BASIS_NONE,
                           imag.qdata.vector.value));
  for (const auto *field : {"u", "v", "curl_u", "curl_v"})
  {
    if (curl || field[0] != 'c')
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
                     CeedElemRestriction subset_restriction,
                     const std::vector<int> *positions, Operator &packed)
{
  ScopedRestriction restriction, qdata;
  ScopedBasis basis, qdata_basis;
  ScopedQFunction qf;
  ScopedOperator op;
  MFEM_ASSERT(!imag || !leaf.curl, "Imaginary curl leaves are not packed!");
  const CeedInt num_elem =
      positions ? static_cast<CeedInt>(positions->size()) : leaf.elements;
  ClonePackedRestriction(
      ceed, subset_restriction ? subset_restriction : leaf.active.restriction.value,
      num_elem, leaf.nodes, size, &restriction.value);
  const bool stacked = leaf.element_qdata && leaf.curl && leaf.curl_points < leaf.qpts;
  const bool curl = leaf.curl && !stacked;
  const CeedInt q_size = curl ? 18 : 9;
  if (stacked)
  {
    CreateStackedBasis(ceed, leaf, &basis.value);
  }
  else
  {
    ClonePackedBasis(ceed, leaf, &basis.value);
  }
  if (leaf.element_qdata)
  {
    const CeedInt num_qpts = stacked ? leaf.qpts + leaf.curl_points : leaf.qpts;
    std::vector<CeedInt> nodes(num_qpts, 0);
    for (CeedInt q = leaf.qpts; q < num_qpts; q++)
    {
      nodes[q] = 1;
    }
    CreateElementQDataRestriction(ceed, leaf.qdata.restriction.value, positions, num_elem,
                                  stacked ? std::vector<CeedInt>{0, 9}
                                          : std::vector<CeedInt>{0},
                                  q_size, &qdata.value);
    CreateSelectorBasis(ceed, q_size, stacked ? 2 : 1, nodes, SelectorScale(leaf, num_qpts),
                        &qdata_basis.value);
  }
  else if (positions)
  {
    CreateSubsetQDataRestriction(ceed, leaf.qdata.restriction.value, *positions,
                                 &qdata.value);
  }
  if (imag)
  {
    PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                             ceed, 1, f_apply_imag_3,
                             PalaceQFunctionRelativePath(f_apply_imag_3_loc), &qf.value));
  }
  else
  {
    PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                             ceed, 1, curl ? f_apply_real_33 : f_apply_real_3,
                             curl ? PalaceQFunctionRelativePath(f_apply_real_33_loc)
                                  : PalaceQFunctionRelativePath(f_apply_real_3_loc),
                             &qf.value));
  }
  PalaceCeedCall(
      ceed, CeedQFunctionAddInput(qf.value, "q", q_size,
                                  leaf.element_qdata ? CEED_EVAL_INTERP : CEED_EVAL_NONE));
  PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "u", 6, leaf.mode));
  PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "v", 6, leaf.mode));
  if (curl)
  {
    PalaceCeedCall(ceed, CeedQFunctionAddInput(qf.value, "curl_u", 6, CEED_EVAL_CURL));
    PalaceCeedCall(ceed, CeedQFunctionAddOutput(qf.value, "curl_v", 6, CEED_EVAL_CURL));
  }
  PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf.value, CEED_QFUNCTION_NONE,
                                          CEED_QFUNCTION_NONE, &op.value));
  PalaceCeedCall(
      ceed, CeedOperatorSetField(op.value, "q",
                                 qdata.value ? qdata.value : leaf.qdata.restriction.value,
                                 leaf.element_qdata ? qdata_basis.value : CEED_BASIS_NONE,
                                 leaf.qdata.vector.value));
  for (const auto *field : {"u", "v", "curl_u", "curl_v"})
  {
    if (curl || field[0] != 'c')
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

  ScopedRestriction restriction;
  restriction.value =
      FiniteElementSpace::BuildCeedElemRestriction(fespace, ceed, subset.geom, elements);
  AddPackedSingle(ceed, leaf, false, size, restriction.value, &positions, packed);
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
  // Whether the values are stored in single precision (in single_values) instead of
  // double precision (in values).
  bool single = false;
  std::vector<CeedInt> indices;
  std::vector<CeedScalar> values;
  std::vector<float> single_values;

  // The indices mapped by a local dof map (see ceed::AddMultMapped) and its key.
  mutable std::vector<CeedInt> mapped_indices;
  mutable std::uint64_t mapped_key = 0;

  // The leaf and its passive inputs (QData) and QFunction context, with their states when
  // the values were assembled. The values of passive inputs may be updated in place after
  // construction, in which case the element matrices are reassembled before use.
  CeedOperator leaf = nullptr;
  std::vector<CeedVector> passive;
  std::vector<uint64_t> states;
  CeedQFunctionContext context = nullptr;
  uint64_t context_state = 0;
  // Copy of the context data the values were assembled with. A context can be shared by
  // several leaves (the coarse operators of a p-multigrid hierarchy reuse the QFunction of
  // the fine operator), and applying or assembling any of them advances the shared state
  // without changing the data.
  std::vector<char> context_data;

  ElementMatrices() = default;
  ElementMatrices(const ElementMatrices &) = delete;
  ElementMatrices &operator=(const ElementMatrices &) = delete;

  std::size_t Bytes() const
  {
    return indices.size() * sizeof(CeedInt) + values.size() * sizeof(CeedScalar) +
           single_values.size() * sizeof(float);
  }

  ~ElementMatrices()
  {
    for (auto &vec : passive)
    {
      PalaceCeedCallBackend(CeedVectorDestroy(&vec));
    }
    PalaceCeedCallBackend(CeedQFunctionContextDestroy(&context));
    PalaceCeedCallBackend(CeedOperatorDestroy(&leaf));
  }
};

}  // namespace

// A leaf of the real (part = 0) or imaginary (part = 1) part of a term, with a reference to
// its active restriction and its element matrices (null if it has none).
struct ElementMatrixTerms::Leaf
{
  int term = 0, part = 0;
  CeedElemRestriction restr = nullptr;
  std::shared_ptr<ElementMatrices> em;
};

namespace
{

// Element assembly replaces the streaming of passive QData and the basis contractions of
// partial assembly by the streaming of the element matrices, with far fewer floating point
// operations. It is considered for single-component leaves of non-tensor (simplex) H1 and
// boundary bases with up to 45 nodes per element (H1 p <= 4 on tetrahedra); larger and
// tensor-product elements keep partial assembly. For the single-precision approximation of
// an operator (see CreateSinglePrecisionOperator), volume H(curl) and H(div) leaves are
// also considered, with up to 84 nodes per element (H(curl) p <= 4 on tetrahedra). With
// structural, only whether the leaf would be considered is checked, for the leaves of the
// terms of sums which may be element-assembled (see ElementMatrixTerms): the choice
// depends on the QFunction of the sum.
bool IsElementAssemblyLeaf(Ceed ceed, CeedOperator op, int size, int num_rhs,
                           bool single = false, bool approximation = false,
                           bool structural = false)
{
  const CeedInt max_nodes = approximation ? 84 : 45;
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
  if (structural)
  {
    // See the exclusion of volume H(curl) and H(div) leaves below.
    return space == CEED_FE_SPACE_H1 || dim != 3 || approximation;
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
  CeedOperatorField *op_in;
  PalaceCeedCall(ceed, CeedOperatorGetQFunction(op, &qf.value));
  PalaceCeedCall(ceed,
                 CeedQFunctionGetFields(qf.value, &num_in, &qf_in, &num_out, &qf_out));
  PalaceCeedCall(ceed, CeedOperatorGetFields(op, nullptr, &op_in, nullptr, nullptr));
  double pa_bytes = 4.0 * nodes, pa_flops = 0.0;
  for (CeedInt i = 0; i < num_in + num_out; i++)
  {
    CeedEvalMode mode;
    CeedInt field_size;
    PalaceCeedCall(ceed,
                   CeedQFunctionFieldGetData(i < num_in ? qf_in[i] : qf_out[i - num_in],
                                             nullptr, &field_size, &mode));
    // Passive inputs with a basis hold QData stored once per element (see
    // ceed::UsePerElementQuadratureData). They count as the QData at every point which the
    // criterion was calibrated with, so that the choice does not depend on the storage.
    bool passive = false;
    if (i < num_in && mode != CEED_EVAL_NONE && mode != CEED_EVAL_WEIGHT)
    {
      ScopedVector vec;
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(op_in[i], &vec.value));
      passive =
          vec.value && vec.value != CEED_VECTOR_ACTIVE && vec.value != CEED_VECTOR_NONE;
    }
    if (mode == CEED_EVAL_NONE || passive)
    {
      pa_bytes += 8.0 * field_size * num_qpts;
    }
    else if (mode != CEED_EVAL_WEIGHT)
    {
      pa_flops += 2.0 * field_size * num_qpts * nodes;
    }
  }
  // Volume H(curl)/H(div) leaves keep partial assembly: element matrices of curl-curl terms
  // do not preserve the discrete curl-gradient cancellation of the factored action, which
  // perturbs ill-conditioned (e.g. near-resonant driven) solves beyond rounding, and their
  // mass terms are mostly bandwidth bound either way. Applications which tolerate the
  // single-precision rounding of the element matrices tolerate this perturbation too.
  if (space != CEED_FE_SPACE_H1 && dim == 3 && !approximation)
  {
    return false;
  }
  const double ea_bytes = (single ? 2.0 : 4.0) * nodes * (nodes + 1) + 4.0 * nodes;
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

// The restriction of a leaf field to the elements [begin, end), or nullptr if the type of
// the restriction is not supported.
CeedElemRestriction RestrictElements(Ceed ceed, CeedElemRestriction full, CeedInt begin,
                                     CeedInt end)
{
  CeedElemRestriction subset = nullptr;
  CeedRestrictionType type;
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(full, &type));
  if (type == CEED_RESTRICTION_STRIDED)
  {
    std::vector<int> positions(end - begin);
    std::iota(positions.begin(), positions.end(), begin);
    CreateSubsetQDataRestriction(ceed, full, positions, &subset);
    return subset;
  }
  if (type != CEED_RESTRICTION_STANDARD && type != CEED_RESTRICTION_ORIENTED &&
      type != CEED_RESTRICTION_CURL_ORIENTED)
  {
    return nullptr;
  }
  CeedInt nodes, num_comp, comp_stride;
  CeedSize length;
  const CeedInt *offsets;
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(full, &nodes));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(full, &num_comp));
  PalaceCeedCall(ceed, CeedElemRestrictionGetCompStride(full, &comp_stride));
  PalaceCeedCall(ceed, CeedElemRestrictionGetLVectorSize(full, &length));
  PalaceCeedCall(ceed, CeedElemRestrictionGetOffsets(full, CEED_MEM_HOST, &offsets));
  const std::size_t first = static_cast<std::size_t>(begin) * nodes;
  if (type == CEED_RESTRICTION_CURL_ORIENTED)
  {
    // Tridiagonal orientation transforms, three entries per node.
    const CeedInt8 *orientations;
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetCurlOrientations(full, CEED_MEM_HOST, &orientations));
    PalaceCeedCall(ceed, CeedElemRestrictionCreateCurlOriented(
                             ceed, end - begin, nodes, num_comp, comp_stride, length,
                             CEED_MEM_HOST, CEED_COPY_VALUES, offsets + first,
                             orientations + 3 * first, &subset));
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreCurlOrientations(full, &orientations));
  }
  else if (type == CEED_RESTRICTION_ORIENTED)
  {
    const bool *orientations;
    PalaceCeedCall(ceed,
                   CeedElemRestrictionGetOrientations(full, CEED_MEM_HOST, &orientations));
    PalaceCeedCall(ceed, CeedElemRestrictionCreateOriented(
                             ceed, end - begin, nodes, num_comp, comp_stride, length,
                             CEED_MEM_HOST, CEED_COPY_VALUES, offsets + first,
                             orientations + first, &subset));
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreOrientations(full, &orientations));
  }
  else
  {
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                             ceed, end - begin, nodes, num_comp, comp_stride, length,
                             CEED_MEM_HOST, CEED_COPY_VALUES, offsets + first, &subset));
  }
  PalaceCeedCall(ceed, CeedElemRestrictionRestoreOffsets(full, &offsets));
  return subset;
}

// A new operator with the QFunction and fields of a leaf over its elements [begin, end), or
// nullptr if a field restriction cannot be restricted to them.
CeedOperator CopyLeafElements(CeedOperator leaf, CeedInt begin, CeedInt end)
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
  // Fields sharing a restriction (the active input and output) share its restriction to
  // the elements.
  std::vector<std::pair<CeedElemRestriction, CeedElemRestriction>> subsets;
  bool valid = true;
  for (CeedInt i = 0; valid && i < num_in + num_out; i++)
  {
    const char *name;
    ScopedRestriction restriction;
    ScopedBasis basis;
    ScopedVector vector;
    PalaceCeedCall(ceed, CeedOperatorFieldGetData(i < num_in ? in[i] : out[i - num_in],
                                                  &name, &restriction.value, &basis.value,
                                                  &vector.value));
    CeedElemRestriction subset = restriction.value;
    if (restriction.value != CEED_ELEMRESTRICTION_NONE)
    {
      auto it = std::find_if(subsets.begin(), subsets.end(), [&](const auto &entry)
                             { return entry.first == restriction.value; });
      if (it == subsets.end())
      {
        subset = RestrictElements(ceed, restriction.value, begin, end);
        valid = (subset != nullptr);
        if (valid)
        {
          subsets.emplace_back(restriction.value, subset);
        }
      }
      else
      {
        subset = it->second;
      }
    }
    if (valid)
    {
      PalaceCeedCall(
          ceed, CeedOperatorSetField(copy.value, name, subset, basis.value, vector.value));
    }
  }
  for (auto &[full, subset] : subsets)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&subset));
  }
  if (!valid)
  {
    return nullptr;
  }
  PalaceCeedCall(ceed, CeedOperatorCheckReady(copy.value));
  return std::exchange(copy.value, nullptr);
}

// Record the states of the passive inputs and QFunction context of the leaf which the
// values of em correspond to.
void RecordElementMatricesState(Ceed ceed, ElementMatrices &em)
{
  for (std::size_t k = 0; k < em.passive.size(); k++)
  {
    PalaceCeedCall(ceed, CeedVectorGetState(em.passive[k], &em.states[k]));
  }
  if (em.context)
  {
    PalaceCeedCall(ceed, CeedQFunctionContextGetState(em.context, &em.context_state));
    std::size_t context_size;
    PalaceCeedCall(ceed, CeedQFunctionContextGetContextSize(em.context, &context_size));
    const void *data;
    PalaceCeedCall(ceed, CeedQFunctionContextGetDataRead(em.context, CEED_MEM_HOST, &data));
    em.context_data.assign(static_cast<const char *>(data),
                           static_cast<const char *>(data) + context_size);
    PalaceCeedCall(ceed, CeedQFunctionContextRestoreDataRead(em.context, &data));
  }
}

// Assemble the element matrices of a leaf into em. With structure, also extract and
// validate its element structure, returning false if the assembled entries do not have the
// expected element block layout. The leaf is assembled through copies over chunks of its
// elements, which bounds the temporary memory of the assembly.
bool FillElementMatrices(Ceed ceed, ElementMatrices &em, int size, bool structure)
{
  const std::size_t P = em.nodes, ne = em.elements, S = P * (P + 1) / 2;
  constexpr std::size_t chunk_bytes = std::size_t(64) << 20;
  std::size_t chunk = std::max(std::size_t(1), chunk_bytes / (sizeof(CeedScalar) * P * P));
  if (chunk < ne)
  {
    ScopedOperator test;
    test.value = CopyLeafElements(em.leaf, 0, 1);
    if (!test.value)
    {
      chunk = ne;
    }
  }
  if (structure)
  {
    em.indices.resize(ne * P);
  }
  // Values are stored as the upper triangles of symmetric matrices until an element matrix
  // is not symmetric up to rounding, from which point all are stored in full.
  bool symmetric = true;
  auto Resize = [&](auto &values) { values.assign(ne * S, 0); };
  auto Expand = [&](auto &values, std::size_t elements)
  {
    std::remove_reference_t<decltype(values)> full(ne * P * P);
    for (std::size_t e = 0; e < elements; e++)
    {
      const auto *u = values.data() + e * S;
      auto *a = full.data() + e * P * P;
      for (std::size_t i = 0; i < P; i++)
      {
        a[i * P + i] = *u++;
        for (std::size_t j = i + 1; j < P; j++)
        {
          a[i * P + j] = a[j * P + i] = *u++;
        }
      }
    }
    values.swap(full);
  };
  auto Store = [&](auto &values, const CeedScalar *v, std::size_t begin, std::size_t count)
  {
    using T = typename std::decay_t<decltype(values)>::value_type;
    if (symmetric)
    {
      T *u = values.data() + begin * S;
      for (std::size_t e = 0; e < count; e++)
      {
        const CeedScalar *a = v + e * P * P;
        for (std::size_t i = 0; i < P; i++)
        {
          *u++ = static_cast<T>(a[i * P + i]);
          for (std::size_t j = i + 1; j < P; j++)
          {
            *u++ = static_cast<T>(0.5 * (a[i * P + j] + a[j * P + i]));
          }
        }
      }
    }
    else
    {
      std::transform(v, v + count * P * P, values.begin() + begin * P * P,
                     [](CeedScalar a) { return static_cast<T>(a); });
    }
  };
  if (em.single)
  {
    std::vector<CeedScalar>().swap(em.values);
    Resize(em.single_values);
  }
  else
  {
    std::vector<float>().swap(em.single_values);
    Resize(em.values);
  }
  for (std::size_t begin = 0; begin < ne; begin += chunk)
  {
    const std::size_t count = std::min(chunk, ne - begin);
    ScopedOperator leaf;
    leaf.value = (count == ne) ? CopyLeaf(em.leaf)
                               : CopyLeafElements(em.leaf, static_cast<CeedInt>(begin),
                                                  static_cast<CeedInt>(begin + count));
    if (structure)
    {
      // The entries come in element blocks, each with the element's row coefficients in
      // order and the same column coefficients, since the leaf uses one active
      // restriction. The sparsity is released before the values are assembled.
      CeedSize nnz;
      CeedInt *rows, *cols;
      PalaceCeedCall(ceed,
                     CeedOperatorLinearAssembleSymbolic(leaf.value, &nnz, &rows, &cols));
      bool valid = (nnz == static_cast<CeedSize>(count * P * P));
      for (std::size_t e = 0; valid && e < count; e++)
      {
        const CeedInt *r = rows + e * P * P, *c = cols + e * P * P;
        for (std::size_t i = 0; valid && i < P; i++)
        {
          em.indices[(begin + e) * P + i] = r[i * P];
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
        return false;
      }
    }
    ScopedVector vals;
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, count * P * P, &vals.value));
    PalaceCeedCall(ceed, CeedOperatorLinearAssemble(leaf.value, vals.value));
    // Release the backend data of the assembled copy before storing the values.
    PalaceCeedCall(ceed, CeedOperatorDestroy(&leaf.value));
    const CeedScalar *v;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(vals.value, CEED_MEM_HOST, &v));
    for (std::size_t e = 0; symmetric && e < count; e++)
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
      if (!symmetric)
      {
        em.single ? Expand(em.single_values, begin) : Expand(em.values, begin);
      }
    }
    em.single ? Store(em.single_values, v, begin, count)
              : Store(em.values, v, begin, count);
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(vals.value, &v));
  }
  em.symmetric = symmetric;
  RecordElementMatricesState(ceed, em);
  return true;
}

// Element matrices (without values) for a leaf, with references to the leaf and to its
// passive inputs and QFunction context.
std::unique_ptr<ElementMatrices> NewElementMatrices(Ceed ceed, CeedOperator op, bool single)
{
  ScopedRestriction in, out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestrictions(op, &in.value, &out.value));
  auto em = std::make_unique<ElementMatrices>();
  em->single = single;
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(in.value, &em->elements));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(in.value, &em->nodes));
  PalaceCeedCall(ceed, CeedOperatorReferenceCopy(op, &em->leaf));
  CeedInt num_in;
  CeedOperatorField *inputs;
  PalaceCeedCall(ceed, CeedOperatorGetFields(op, &num_in, &inputs, nullptr, nullptr));
  for (CeedInt i = 0; i < num_in; i++)
  {
    CeedVector vec = nullptr;
    PalaceCeedCall(ceed, CeedOperatorFieldGetVector(inputs[i], &vec));
    if (vec && vec != CEED_VECTOR_ACTIVE && vec != CEED_VECTOR_NONE)
    {
      em->passive.push_back(vec);
      em->states.push_back(0);
    }
    else
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&vec));
    }
  }
  PalaceCeedCall(ceed, CeedOperatorGetContext(op, &em->context));
  return em;
}

std::unique_ptr<ElementMatrices> AssembleElementMatrices(Ceed ceed, CeedOperator op,
                                                         int size, bool single = false)
{
  auto em = NewElementMatrices(ceed, op, single);
  if (!FillElementMatrices(ceed, *em, size, true))
  {
    return {};
  }
  return em;
}

// Reassemble the element matrices if any passive input or the context of the leaf changed.
void RefreshElementMatrices(ElementMatrices &em)
{
  Ceed ceed = CeedOperatorReturnCeed(em.leaf);
  bool current = true;
  for (std::size_t k = 0; current && k < em.passive.size(); k++)
  {
    uint64_t state;
    PalaceCeedCall(ceed, CeedVectorGetState(em.passive[k], &state));
    current = (state == em.states[k]);
  }
  if (current && em.context)
  {
    uint64_t state;
    PalaceCeedCall(ceed, CeedQFunctionContextGetState(em.context, &state));
    if (state != em.context_state)
    {
      // Compare the data itself: the state of a shared context also advances when another
      // leaf is applied or assembled. Reassembling here would then repeat at every switch
      // between the levels of a multigrid hierarchy.
      std::size_t context_size;
      PalaceCeedCall(ceed, CeedQFunctionContextGetContextSize(em.context, &context_size));
      const void *data;
      PalaceCeedCall(ceed,
                     CeedQFunctionContextGetDataRead(em.context, CEED_MEM_HOST, &data));
      current = (context_size == em.context_data.size() &&
                 std::memcmp(data, em.context_data.data(), context_size) == 0);
      PalaceCeedCall(ceed, CeedQFunctionContextRestoreDataRead(em.context, &data));
      if (current)
      {
        em.context_state = state;
      }
    }
  }
  if (!current)
  {
    MFEM_VERIFY(FillElementMatrices(ceed, em, std::numeric_limits<int>::max(), false),
                "Failed to reassemble element matrices after an update of their inputs!");
  }
}

// Release the references of em to its leaf and to the inputs of the leaf, for element
// matrices which are not refreshed.
void ReleaseElementMatricesLeaf(ElementMatrices &em)
{
  for (auto &vec : em.passive)
  {
    PalaceCeedCallBackend(CeedVectorDestroy(&vec));
  }
  em.passive.clear();
  em.states.clear();
  PalaceCeedCallBackend(CeedQFunctionContextDestroy(&em.context));
  std::vector<char>().swap(em.context_data);
  PalaceCeedCallBackend(CeedOperatorDestroy(&em.leaf));
}

// The weight of a leaf of a term T = Tr + i Ti with coefficient c in the real part
// Re(c) Tr - Im(c) Ti, or with imag in the imaginary part Im(c) Tr + Re(c) Ti, of a
// combination.
double CombinationWeight(const ElementMatrixCombination &comb,
                         const ElementMatrixTerms::Leaf &leaf, bool imag)
{
  const auto c = comb.coeff[leaf.term];
  return imag ? (leaf.part ? c.real() : c.imag()) : (leaf.part ? -c.imag() : c.real());
}

// Whether each child of the real part (or with imag, the imaginary part) of a combination
// of terms may have its element matrices formed from those of the terms (see
// ElementMatrixTerms): its restriction is that of no other child, and every weighted leaf
// of the terms has the restriction of a child.
std::vector<bool> CombinableLeaves(Ceed ceed, const ElementMatrixCombination &comb,
                                   bool imag, bool approximation, int num_rhs,
                                   CeedInt count, CeedOperator *children)
{
  std::vector<bool> combinable(count, false);
  const auto *terms = comb.terms.get();
  if (!terms || !terms->IsValid() || terms->IsApproximation() != approximation ||
      terms->NumRHS() != num_rhs ||
      static_cast<int>(comb.coeff.size()) != terms->NumTerms())
  {
    return combinable;
  }
  std::vector<CeedElemRestriction> restr(count, nullptr);
  for (CeedInt j = 0; j < count; j++)
  {
    ScopedRestriction in, out;
    PalaceCeedCall(
        ceed, CeedOperatorGetActiveElemRestrictions(children[j], &in.value, &out.value));
    if (in.value && in.value == out.value)
    {
      restr[j] = in.value;
    }
  }
  for (const auto &leaf : terms->GetLeaves())
  {
    if (CombinationWeight(comb, leaf, imag) != 0.0 &&
        std::find(restr.begin(), restr.end(), leaf.restr) == restr.end())
    {
      return combinable;
    }
  }
  for (CeedInt j = 0; j < count; j++)
  {
    combinable[j] = restr[j] && std::count(restr.begin(), restr.end(), restr[j]) == 1;
  }
  return combinable;
}

// Single-precision element matrices of a child of an operator part (see CombinableLeaves),
// formed from those of the leaves of the terms over its restriction, with the arithmetic in
// double precision. Returns null if the terms have no weighted leaf over the restriction or
// a leaf over it without symmetric single-precision element matrices.
std::unique_ptr<ElementMatrices>
CombineElementMatrices(Ceed ceed, CeedOperator op, const ElementMatrixCombination &comb,
                       bool imag)
{
  ScopedRestriction in, out;
  PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestrictions(op, &in.value, &out.value));
  CeedInt ne, P;
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(in.value, &ne));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(in.value, &P));
  std::vector<std::pair<double, const float *>> sum;
  const std::vector<CeedInt> *indices = nullptr;
  for (const auto &leaf : comb.terms->GetLeaves())
  {
    if (leaf.restr != in.value)
    {
      continue;
    }
    const auto *t = leaf.em.get();
    if (!t || !t->single || !t->symmetric || t->elements != ne || t->nodes != P)
    {
      return {};
    }
    const double w = CombinationWeight(comb, leaf, imag);
    if (w != 0.0)
    {
      sum.emplace_back(w, t->single_values.data());
      indices = &t->indices;
    }
  }
  if (sum.empty())
  {
    return {};
  }
  auto em = NewElementMatrices(ceed, op, true);
  em->symmetric = true;
  em->indices = *indices;
  const std::size_t n = static_cast<std::size_t>(ne) * P * (P + 1) / 2;
  em->single_values.resize(n);
  constexpr std::size_t block = 512;
  double acc[block];
  for (std::size_t b = 0; b < n; b += block)
  {
    const std::size_t m = std::min(block, n - b);
    const auto [w0, t0] = sum.front();
    for (std::size_t i = 0; i < m; i++)
    {
      acc[i] = w0 * t0[b + i];
    }
    for (std::size_t k = 1; k < sum.size(); k++)
    {
      const auto [w, t] = sum[k];
      for (std::size_t i = 0; i < m; i++)
      {
        acc[i] += w * t[b + i];
      }
    }
    float *v = em->single_values.data() + b;
    for (std::size_t i = 0; i < m; i++)
    {
      v[i] = static_cast<float>(acc[i]);
    }
  }
  RecordElementMatricesState(ceed, *em);
  return em;
}

// Products y += A x of an element matrix with the element coefficients of one, or of two
// interleaved (x[2 j + c], c = 0, 1 for the real and imaginary parts), right hand sides,
// with the arithmetic in double precision for values stored in either precision. For the
// packed upper triangle of a symmetric matrix (row i holds the columns j >= i), the
// entries above the diagonal contribute to the rows of both of their indices. Rows are
// processed four at a time, so that the accumulations of the four rows along the columns
// are independent and not bound by the latency of the additions.
template <typename T>
inline void SymmetricElementMult2(const std::size_t P, const T *__restrict__ a,
                                  const CeedScalar *__restrict__ x,
                                  CeedScalar *__restrict__ y)
{
  std::size_t i = 0;
  for (; i + 4 <= P; i += 4)
  {
    const T *__restrict__ a0 = a;
    const T *__restrict__ a1 = a0 + (P - i);
    const T *__restrict__ a2 = a1 + (P - i - 1);
    const T *__restrict__ a3 = a2 + (P - i - 2);
    a = a3 + (P - i - 3);
    const CeedScalar x0r = x[2 * i], x0i = x[2 * i + 1], x1r = x[2 * i + 2],
                     x1i = x[2 * i + 3], x2r = x[2 * i + 4], x2i = x[2 * i + 5],
                     x3r = x[2 * i + 6], x3i = x[2 * i + 7];
    // The diagonal block of rows and columns i, ..., i + 3.
    const CeedScalar d00 = a0[0], d01 = a0[1], d02 = a0[2], d03 = a0[3], d11 = a1[0],
                     d12 = a1[1], d13 = a1[2], d22 = a2[0], d23 = a2[1], d33 = a3[0];
    CeedScalar s0r = d00 * x0r + d01 * x1r + d02 * x2r + d03 * x3r;
    CeedScalar s0i = d00 * x0i + d01 * x1i + d02 * x2i + d03 * x3i;
    CeedScalar s1r = d01 * x0r + d11 * x1r + d12 * x2r + d13 * x3r;
    CeedScalar s1i = d01 * x0i + d11 * x1i + d12 * x2i + d13 * x3i;
    CeedScalar s2r = d02 * x0r + d12 * x1r + d22 * x2r + d23 * x3r;
    CeedScalar s2i = d02 * x0i + d12 * x1i + d22 * x2i + d23 * x3i;
    CeedScalar s3r = d03 * x0r + d13 * x1r + d23 * x2r + d33 * x3r;
    CeedScalar s3i = d03 * x0i + d13 * x1i + d23 * x2i + d33 * x3i;
    // The columns j >= i + 4 of the four rows, b_k[j - i - 4] = A(i + k, j).
    const T *__restrict__ b0 = a0 + 4;
    const T *__restrict__ b1 = a1 + 3;
    const T *__restrict__ b2 = a2 + 2;
    const T *__restrict__ b3 = a3 + 1;
    const CeedScalar *__restrict__ xj = x + 2 * (i + 4);
    CeedScalar *__restrict__ yj = y + 2 * (i + 4);
    for (std::size_t m = 0; m < P - i - 4; m++)
    {
      const CeedScalar c0 = b0[m], c1 = b1[m], c2 = b2[m], c3 = b3[m];
      const CeedScalar xr = xj[2 * m], xi = xj[2 * m + 1];
      s0r += c0 * xr;
      s0i += c0 * xi;
      s1r += c1 * xr;
      s1i += c1 * xi;
      s2r += c2 * xr;
      s2i += c2 * xi;
      s3r += c3 * xr;
      s3i += c3 * xi;
      yj[2 * m] += c0 * x0r + c1 * x1r + c2 * x2r + c3 * x3r;
      yj[2 * m + 1] += c0 * x0i + c1 * x1i + c2 * x2i + c3 * x3i;
    }
    y[2 * i] += s0r;
    y[2 * i + 1] += s0i;
    y[2 * i + 2] += s1r;
    y[2 * i + 3] += s1i;
    y[2 * i + 4] += s2r;
    y[2 * i + 5] += s2i;
    y[2 * i + 6] += s3r;
    y[2 * i + 7] += s3i;
  }
  for (; i < P; i++)
  {
    const T *__restrict__ a0 = a;
    a += P - i;
    const CeedScalar x0r = x[2 * i], x0i = x[2 * i + 1];
    CeedScalar s0r = a0[0] * x0r, s0i = a0[0] * x0i;
    for (std::size_t j = i + 1; j < P; j++)
    {
      const CeedScalar c0 = a0[j - i];
      s0r += c0 * x[2 * j];
      s0i += c0 * x[2 * j + 1];
      y[2 * j] += c0 * x0r;
      y[2 * j + 1] += c0 * x0i;
    }
    y[2 * i] += s0r;
    y[2 * i + 1] += s0i;
  }
}

template <typename T>
inline void FullElementMult2(const std::size_t P, const T *__restrict__ a,
                             const CeedScalar *__restrict__ x, CeedScalar *__restrict__ y)
{
  std::size_t i = 0;
  for (; i + 4 <= P; i += 4)
  {
    const T *__restrict__ a0 = a + i * P;
    const T *__restrict__ a1 = a0 + P;
    const T *__restrict__ a2 = a1 + P;
    const T *__restrict__ a3 = a2 + P;
    CeedScalar s0r = 0.0, s0i = 0.0, s1r = 0.0, s1i = 0.0, s2r = 0.0, s2i = 0.0, s3r = 0.0,
               s3i = 0.0;
    for (std::size_t j = 0; j < P; j++)
    {
      const CeedScalar c0 = a0[j], c1 = a1[j], c2 = a2[j], c3 = a3[j];
      const CeedScalar xr = x[2 * j], xi = x[2 * j + 1];
      s0r += c0 * xr;
      s0i += c0 * xi;
      s1r += c1 * xr;
      s1i += c1 * xi;
      s2r += c2 * xr;
      s2i += c2 * xi;
      s3r += c3 * xr;
      s3i += c3 * xi;
    }
    y[2 * i] += s0r;
    y[2 * i + 1] += s0i;
    y[2 * i + 2] += s1r;
    y[2 * i + 3] += s1i;
    y[2 * i + 4] += s2r;
    y[2 * i + 5] += s2i;
    y[2 * i + 6] += s3r;
    y[2 * i + 7] += s3i;
  }
  for (; i < P; i++)
  {
    const T *__restrict__ a0 = a + i * P;
    CeedScalar s0r = 0.0, s0i = 0.0;
    for (std::size_t j = 0; j < P; j++)
    {
      s0r += a0[j] * x[2 * j];
      s0i += a0[j] * x[2 * j + 1];
    }
    y[2 * i] += s0r;
    y[2 * i + 1] += s0i;
  }
}

template <typename T>
inline void SymmetricElementMult1(const std::size_t P, const T *__restrict__ a,
                                  const CeedScalar *__restrict__ x,
                                  CeedScalar *__restrict__ y)
{
  std::size_t i = 0;
  for (; i + 4 <= P; i += 4)
  {
    const T *__restrict__ a0 = a;
    const T *__restrict__ a1 = a0 + (P - i);
    const T *__restrict__ a2 = a1 + (P - i - 1);
    const T *__restrict__ a3 = a2 + (P - i - 2);
    a = a3 + (P - i - 3);
    const CeedScalar x0 = x[i], x1 = x[i + 1], x2 = x[i + 2], x3 = x[i + 3];
    const CeedScalar d00 = a0[0], d01 = a0[1], d02 = a0[2], d03 = a0[3], d11 = a1[0],
                     d12 = a1[1], d13 = a1[2], d22 = a2[0], d23 = a2[1], d33 = a3[0];
    CeedScalar s0 = d00 * x0 + d01 * x1 + d02 * x2 + d03 * x3;
    CeedScalar s1 = d01 * x0 + d11 * x1 + d12 * x2 + d13 * x3;
    CeedScalar s2 = d02 * x0 + d12 * x1 + d22 * x2 + d23 * x3;
    CeedScalar s3 = d03 * x0 + d13 * x1 + d23 * x2 + d33 * x3;
    const T *__restrict__ b0 = a0 + 4;
    const T *__restrict__ b1 = a1 + 3;
    const T *__restrict__ b2 = a2 + 2;
    const T *__restrict__ b3 = a3 + 1;
    const CeedScalar *__restrict__ xj = x + i + 4;
    CeedScalar *__restrict__ yj = y + i + 4;
    for (std::size_t m = 0; m < P - i - 4; m++)
    {
      const CeedScalar c0 = b0[m], c1 = b1[m], c2 = b2[m], c3 = b3[m];
      const CeedScalar xm = xj[m];
      s0 += c0 * xm;
      s1 += c1 * xm;
      s2 += c2 * xm;
      s3 += c3 * xm;
      yj[m] += c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3;
    }
    y[i] += s0;
    y[i + 1] += s1;
    y[i + 2] += s2;
    y[i + 3] += s3;
  }
  for (; i < P; i++)
  {
    const T *__restrict__ a0 = a;
    a += P - i;
    const CeedScalar x0 = x[i];
    CeedScalar s0 = a0[0] * x0;
    for (std::size_t j = i + 1; j < P; j++)
    {
      const CeedScalar c0 = a0[j - i];
      s0 += c0 * x[j];
      y[j] += c0 * x0;
    }
    y[i] += s0;
  }
}

template <typename T>
inline void FullElementMult1(const std::size_t P, const T *__restrict__ a,
                             const CeedScalar *__restrict__ x, CeedScalar *__restrict__ y)
{
  std::size_t i = 0;
  for (; i + 4 <= P; i += 4)
  {
    const T *__restrict__ a0 = a + i * P;
    const T *__restrict__ a1 = a0 + P;
    const T *__restrict__ a2 = a1 + P;
    const T *__restrict__ a3 = a2 + P;
    CeedScalar s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (std::size_t j = 0; j < P; j++)
    {
      s0 += a0[j] * x[j];
      s1 += a1[j] * x[j];
      s2 += a2[j] * x[j];
      s3 += a3[j] * x[j];
    }
    y[i] += s0;
    y[i + 1] += s1;
    y[i + 2] += s2;
    y[i + 3] += s3;
  }
  for (; i < P; i++)
  {
    const T *__restrict__ a0 = a + i * P;
    CeedScalar s0 = 0.0;
    for (std::size_t j = 0; j < P; j++)
    {
      s0 += a0[j] * x[j];
    }
    y[i] += s0;
  }
}

template <typename T>
const T *ElementValues(const ElementMatrices &em)
{
  if constexpr (std::is_same_v<T, float>)
  {
    return em.single_values.data();
  }
  else
  {
    return em.values.data();
  }
}

// y += (Ar + i Ai) x for x = [xr | xi] and y = [yr | yi], with the element matrices ar of
// the real part and ai of the imaginary part (one of them may be null) over the same
// elements: each element's coefficients are gathered once, interleaved, and the result
// scattered once. The work array holds 6 P values.
template <typename Tr, typename Ti>
void AddMultElements2(const ElementMatrices *ar, const ElementMatrices *ai,
                      const CeedScalar *__restrict__ xr, const CeedScalar *__restrict__ xi,
                      CeedScalar *__restrict__ yr, CeedScalar *__restrict__ yi,
                      CeedScalar *__restrict__ work)
{
  const ElementMatrices &em = ar ? *ar : *ai;
  const std::size_t P = em.nodes, ne = em.elements;
  const std::size_t Sr = ar ? (ar->symmetric ? P * (P + 1) / 2 : P * P) : 0;
  const std::size_t Si = ai ? (ai->symmetric ? P * (P + 1) / 2 : P * P) : 0;
  const Tr *va = ar ? ElementValues<Tr>(*ar) : nullptr;
  const Ti *vb = ai ? ElementValues<Ti>(*ai) : nullptr;
  CeedScalar *__restrict__ xe = work;
  CeedScalar *__restrict__ xs = work + 2 * P;
  CeedScalar *__restrict__ ye = work + 4 * P;
  for (std::size_t e = 0; e < ne; e++)
  {
    const CeedInt *id = em.indices.data() + e * P;
    for (std::size_t j = 0; j < P; j++)
    {
      xe[2 * j] = xr[id[j]];
      xe[2 * j + 1] = xi[id[j]];
      ye[2 * j] = 0.0;
      ye[2 * j + 1] = 0.0;
    }
    if (ar)
    {
      ar->symmetric ? SymmetricElementMult2(P, va + e * Sr, xe, ye)
                    : FullElementMult2(P, va + e * Sr, xe, ye);
    }
    if (ai)
    {
      // The imaginary part acts on i x.
      for (std::size_t j = 0; j < P; j++)
      {
        xs[2 * j] = -xe[2 * j + 1];
        xs[2 * j + 1] = xe[2 * j];
      }
      ai->symmetric ? SymmetricElementMult2(P, vb + e * Si, xs, ye)
                    : FullElementMult2(P, vb + e * Si, xs, ye);
    }
    for (std::size_t i = 0; i < P; i++)
    {
      yr[id[i]] += ye[2 * i];
      yi[id[i]] += ye[2 * i + 1];
    }
  }
}

// y += A x for a single real vector. The work array holds 2 P values.
template <typename T>
void AddMultElements1(const ElementMatrices &em, const CeedScalar *__restrict__ x,
                      CeedScalar *__restrict__ y, CeedScalar *__restrict__ work)
{
  const std::size_t P = em.nodes, ne = em.elements;
  const std::size_t S = em.symmetric ? P * (P + 1) / 2 : P * P;
  const T *va = ElementValues<T>(em);
  CeedScalar *__restrict__ xe = work;
  CeedScalar *__restrict__ ye = work + P;
  for (std::size_t e = 0; e < ne; e++)
  {
    const CeedInt *id = em.indices.data() + e * P;
    for (std::size_t j = 0; j < P; j++)
    {
      xe[j] = x[id[j]];
      ye[j] = 0.0;
    }
    em.symmetric ? SymmetricElementMult1(P, va + e * S, xe, ye)
                 : FullElementMult1(P, va + e * S, xe, ye);
    for (std::size_t i = 0; i < P; i++)
    {
      y[id[i]] += ye[i];
    }
  }
}

// y += (Ar + i Ai) x (see AddMultElements2), after reassembling element matrices whose
// inputs changed.
void AddMultElementPair(ElementMatrices *ar, ElementMatrices *ai, const CeedScalar *xr,
                        const CeedScalar *xi, CeedScalar *yr, CeedScalar *yi, Vector &work)
{
  for (auto *em : {ar, ai})
  {
    if (em)
    {
      RefreshElementMatrices(*em);
    }
  }
  const std::size_t P = (ar ? ar : ai)->nodes;
  if (static_cast<std::size_t>(work.Size()) < 6 * P)
  {
    work.SetSize(6 * P);
  }
  CeedScalar *w = work.HostWrite();
  const bool sr = ar ? ar->single : ai->single, si = ai ? ai->single : sr;
  if (sr && si)
  {
    AddMultElements2<float, float>(ar, ai, xr, xi, yr, yi, w);
  }
  else if (sr)
  {
    AddMultElements2<float, CeedScalar>(ar, ai, xr, xi, yr, yi, w);
  }
  else if (si)
  {
    AddMultElements2<CeedScalar, float>(ar, ai, xr, xi, yr, yi, w);
  }
  else
  {
    AddMultElements2<CeedScalar, CeedScalar>(ar, ai, xr, xi, yr, yi, w);
  }
}

// The element indices of em mapped by a local dof map (built at the first use of the map).
const CeedInt *MappedIndices(const ElementMatrices &em, const LocalDofMap &map)
{
  if (em.mapped_key != map.key)
  {
    em.mapped_indices.resize(em.indices.size());
    for (std::size_t k = 0; k < em.indices.size(); k++)
    {
      em.mapped_indices[k] = map.map[em.indices[k]];
    }
    em.mapped_key = map.key;
  }
  return em.mapped_indices.data();
}

// Entry m of a local dof map (see LocalDofMap): true dof m ≥ 0, or extra dof -1 - m.
inline std::size_t MappedEntry(CeedInt m)
{
  return static_cast<std::size_t>(m ^ (m >> 31));
}

// AddMultElements2 with the element coefficients read from, and the results added to, the
// true-dof vectors (xtr, xti, ytr, yti) or the extra-dof vectors (xxr, ...) of a local dof
// map (mapped indices mid).
template <typename Tr, typename Ti>
void AddMultElements2Mapped(const ElementMatrices *ar, const ElementMatrices *ai,
                            const CeedInt *mid, const CeedScalar *xtr,
                            const CeedScalar *xti, const CeedScalar *xxr,
                            const CeedScalar *xxi, CeedScalar *ytr, CeedScalar *yti,
                            CeedScalar *yxr, CeedScalar *yxi, CeedScalar *__restrict__ work)
{
  const ElementMatrices &em = ar ? *ar : *ai;
  const std::size_t P = em.nodes, ne = em.elements;
  const std::size_t Sr = ar ? (ar->symmetric ? P * (P + 1) / 2 : P * P) : 0;
  const std::size_t Si = ai ? (ai->symmetric ? P * (P + 1) / 2 : P * P) : 0;
  const Tr *va = ar ? ElementValues<Tr>(*ar) : nullptr;
  const Ti *vb = ai ? ElementValues<Ti>(*ai) : nullptr;
  CeedScalar *__restrict__ xe = work;
  CeedScalar *__restrict__ xs = work + 2 * P;
  CeedScalar *__restrict__ ye = work + 4 * P;
  for (std::size_t e = 0; e < ne; e++)
  {
    const CeedInt *id = mid + e * P;
    for (std::size_t j = 0; j < P; j++)
    {
      const bool t = id[j] >= 0;
      const std::size_t i = MappedEntry(id[j]);
      xe[2 * j] = (t ? xtr : xxr)[i];
      xe[2 * j + 1] = (t ? xti : xxi)[i];
      ye[2 * j] = 0.0;
      ye[2 * j + 1] = 0.0;
    }
    if (ar)
    {
      ar->symmetric ? SymmetricElementMult2(P, va + e * Sr, xe, ye)
                    : FullElementMult2(P, va + e * Sr, xe, ye);
    }
    if (ai)
    {
      for (std::size_t j = 0; j < P; j++)
      {
        xs[2 * j] = -xe[2 * j + 1];
        xs[2 * j + 1] = xe[2 * j];
      }
      ai->symmetric ? SymmetricElementMult2(P, vb + e * Si, xs, ye)
                    : FullElementMult2(P, vb + e * Si, xs, ye);
    }
    for (std::size_t j = 0; j < P; j++)
    {
      const bool t = id[j] >= 0;
      const std::size_t i = MappedEntry(id[j]);
      (t ? ytr : yxr)[i] += ye[2 * j];
      (t ? yti : yxi)[i] += ye[2 * j + 1];
    }
  }
}

template <typename T>
void AddMultElements1Mapped(const ElementMatrices &em, const CeedInt *mid,
                            const CeedScalar *xt, const CeedScalar *xx, CeedScalar *yt,
                            CeedScalar *yx, CeedScalar *__restrict__ work)
{
  const std::size_t P = em.nodes, ne = em.elements;
  const std::size_t S = em.symmetric ? P * (P + 1) / 2 : P * P;
  const T *va = ElementValues<T>(em);
  CeedScalar *__restrict__ xe = work;
  CeedScalar *__restrict__ ye = work + P;
  for (std::size_t e = 0; e < ne; e++)
  {
    const CeedInt *id = mid + e * P;
    for (std::size_t j = 0; j < P; j++)
    {
      xe[j] = (id[j] >= 0 ? xt : xx)[MappedEntry(id[j])];
      ye[j] = 0.0;
    }
    em.symmetric ? SymmetricElementMult1(P, va + e * S, xe, ye)
                 : FullElementMult1(P, va + e * S, xe, ye);
    for (std::size_t j = 0; j < P; j++)
    {
      (id[j] >= 0 ? yt : yx)[MappedEntry(id[j])] += ye[j];
    }
  }
}

void AddMultElementPairMapped(ElementMatrices *ar, ElementMatrices *ai,
                              const LocalDofMap &map, const CeedScalar *xtr,
                              const CeedScalar *xti, const CeedScalar *xxr,
                              const CeedScalar *xxi, CeedScalar *ytr, CeedScalar *yti,
                              CeedScalar *yxr, CeedScalar *yxi, Vector &work)
{
  for (auto *em : {ar, ai})
  {
    if (em)
    {
      RefreshElementMatrices(*em);
    }
  }
  const ElementMatrices &em = ar ? *ar : *ai;
  const std::size_t P = em.nodes;
  if (static_cast<std::size_t>(work.Size()) < 6 * P)
  {
    work.SetSize(6 * P);
  }
  CeedScalar *w = work.HostWrite();
  const CeedInt *mid = MappedIndices(em, map);
  const bool sr = ar ? ar->single : ai->single, si = ai ? ai->single : sr;
  if (sr && si)
  {
    AddMultElements2Mapped<float, float>(ar, ai, mid, xtr, xti, xxr, xxi, ytr, yti, yxr,
                                         yxi, w);
  }
  else if (sr)
  {
    AddMultElements2Mapped<float, CeedScalar>(ar, ai, mid, xtr, xti, xxr, xxi, ytr, yti,
                                              yxr, yxi, w);
  }
  else if (si)
  {
    AddMultElements2Mapped<CeedScalar, float>(ar, ai, mid, xtr, xti, xxr, xxi, ytr, yti,
                                              yxr, yxi, w);
  }
  else
  {
    AddMultElements2Mapped<CeedScalar, CeedScalar>(ar, ai, mid, xtr, xti, xxr, xxi, ytr,
                                                   yti, yxr, yxi, w);
  }
}

void AddMultElementMatricesMapped(ElementMatrices &em, const LocalDofMap &map,
                                  const CeedScalar *xt, const CeedScalar *xx,
                                  CeedScalar *yt, CeedScalar *yx, Vector &work)
{
  RefreshElementMatrices(em);
  if (static_cast<std::size_t>(work.Size()) < 2 * static_cast<std::size_t>(em.nodes))
  {
    work.SetSize(2 * em.nodes);
  }
  CeedScalar *w = work.HostWrite();
  const CeedInt *mid = MappedIndices(em, map);
  em.single ? AddMultElements1Mapped<float>(em, mid, xt, xx, yt, yx, w)
            : AddMultElements1Mapped<CeedScalar>(em, mid, xt, xx, yt, yx, w);
}

// y += A x for a single real vector, after reassembling element matrices whose inputs
// changed.
void AddMultElementMatrices(ElementMatrices &em, const CeedScalar *x, CeedScalar *y,
                            Vector &work)
{
  RefreshElementMatrices(em);
  if (static_cast<std::size_t>(work.Size()) < 2 * static_cast<std::size_t>(em.nodes))
  {
    work.SetSize(2 * em.nodes);
  }
  CeedScalar *w = work.HostWrite();
  em.single ? AddMultElements1<float>(em, x, y, w)
            : AddMultElements1<CeedScalar>(em, x, y, w);
}

}  // namespace

struct Operator::ElementAssembly
{
  // Element-assembled children of the composite, and the composite of the others.
  std::vector<CeedInt> children;
  std::vector<std::shared_ptr<ElementMatrices>> leaves;
  CeedOperator remainder = nullptr;
  mutable Vector work;

  const std::shared_ptr<ElementMatrices> *Find(CeedInt child) const
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
  std::shared_ptr<ElementMatrices> matrices;
  bool imag;
};

void AddMultElementMatrices(const std::vector<ElementLeaf> &ea, const CeedScalar *xr,
                            const CeedScalar *xi, CeedScalar *yr, CeedScalar *yi,
                            Vector &work)
{
  for (const auto &[em, imag] : ea)
  {
    AddMultElementPair(imag ? nullptr : em.get(), imag ? em.get() : nullptr, xr, xi, yr, yi,
                       work);
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

  // Only element-assembled terms (as ApplyInPlace).
  bool SupportsMapped() const { return !packed && !remainder && !ea.empty(); }

  void AddMultMapped(const LocalDofMap &map, const CeedScalar *xtr, const CeedScalar *xti,
                     const CeedScalar *xxr, const CeedScalar *xxi, CeedScalar *ytr,
                     CeedScalar *yti, CeedScalar *yxr, CeedScalar *yxi) const
  {
    for (const auto &[em, imag] : ea)
    {
      AddMultElementPairMapped(imag ? nullptr : em.get(), imag ? em.get() : nullptr, map,
                               xtr, xti, xxr, xxi, ytr, yti, yxr, yxi, work);
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
  // action share the same element matrices.
  const Operator::ElementAssembly *part_ea[2] = {nullptr, nullptr};
  for (int part = 0; part < 2; part++)
  {
    auto *part_op = const_cast<Operator *>(part == 0 ? real : imag);
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
          fields[0][r]->mode != fields[1][i]->mode ||
          fields[0][r]->active.basis.value != fields[1][i]->active.basis.value)
      {
        continue;
      }
      // A pair reads the QData of both leaves per element or both per point, and QData
      // stored once per element can only be read per element. Such leaves without a
      // matching partner are packed on their own below.
      if (fields[0][r]->element_qdata != fields[1][i]->element_qdata &&
          (fields[0][r]->stored_per_element || fields[1][i]->stored_per_element))
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
      AddPackedSingle(ceed, *fields[part][j], part == 1, size, nullptr, nullptr, *packed);
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

// Largest magnitude of the entries of each row of the element matrices, maximized into
// scale over the elements sharing an L-vector coefficient. Returns false if a row of an
// element exceeds tol times the given scale of its coefficient, when tol > 0.
bool ElementRowMaxima(const ElementMatrices &em, std::vector<double> &scale, double tol)
{
  const std::size_t P = em.nodes, ne = em.elements;
  std::vector<double> row(P);
  auto Rows = [&](const auto *v)
  {
    for (std::size_t e = 0; e < ne; e++)
    {
      std::fill(row.begin(), row.end(), 0.0);
      if (em.symmetric)
      {
        const auto *a = v + e * P * (P + 1) / 2;
        for (std::size_t i = 0; i < P; i++)
        {
          for (std::size_t j = i; j < P; j++)
          {
            const double aij = std::abs(static_cast<double>(*a++));
            row[i] = std::max(row[i], aij);
            row[j] = std::max(row[j], aij);
          }
        }
      }
      else
      {
        const auto *a = v + e * P * P;
        for (std::size_t i = 0; i < P; i++)
        {
          for (std::size_t j = 0; j < P; j++)
          {
            row[i] = std::max(row[i], std::abs(static_cast<double>(a[i * P + j])));
          }
        }
      }
      const CeedInt *id = em.indices.data() + e * P;
      for (std::size_t i = 0; i < P; i++)
      {
        if (tol > 0.0)
        {
          if (row[i] > tol * scale[id[i]])
          {
            return false;
          }
        }
        else
        {
          scale[id[i]] = std::max(scale[id[i]], row[i]);
        }
      }
    }
    return true;
  };
  return em.single ? Rows(em.single_values.data()) : Rows(em.values.data());
}

// Whether a real or complex operator can be approximated in single precision, and its
// composite children.
bool SinglePrecisionParts(const std::array<const Operator *, 2> &parts, int &size,
                          std::array<CeedInt, 2> &count,
                          std::array<CeedOperator *, 2> &children)
{
  if ((!parts[0] && !parts[1]) || internal::NumCeeds() != 1 || utils::GetMaxThreads() > 1 ||
      utils::InParallel() || mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    return false;
  }
  const auto *first = parts[0] ? parts[0] : parts[1];
  for (const auto *part : parts)
  {
    if (part &&
        (part->Size() != 1 || part->HasDofMultiplicity() ||
         part->Width() != part->Height() || part->Width() != first->Width() ||
         (typeid(*part) != typeid(Operator) && typeid(*part) != typeid(SymmetricOperator))))
    {
      return false;
    }
  }
  size = first->Width();
  if (size <= 0)
  {
    return false;
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  CeedMemType mem;
  const char *resource;
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  if (mem != CEED_MEM_HOST || std::strncmp(resource, "/cpu/", 5) != 0)
  {
    return false;
  }
  for (int k = 0; k < 2; k++)
  {
    count[k] = 0;
    children[k] = nullptr;
    if (parts[k])
    {
      if (!SameContext(ceed, CeedOperatorReturnCeed((*parts[k])[0])))
      {
        return false;
      }
      PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub((*parts[k])[0], &count[k]));
      PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList((*parts[k])[0], &children[k]));
    }
  }
  return true;
}

// The element matrices of the children of a part for its single-precision approximation:
// those of the part's own element assembly, or single-precision ones for the eligible
// children (formed from those of the terms of a combination of which the part is the real
// part, or with imag the imaginary part, where possible). The other children are left null.
// Returns the bytes of the new matrices.
std::size_t SinglePrecisionLeaves(Ceed ceed, const Operator &part, CeedInt count,
                                  CeedOperator *children, int size, int num_rhs,
                                  std::vector<std::shared_ptr<ElementMatrices>> &leaves,
                                  bool &any_single,
                                  const ElementMatrixCombination *comb = nullptr,
                                  bool imag = false)
{
  std::size_t bytes = 0;
  const auto *ea = part.GetElementAssembly();
  leaves.assign(count, nullptr);
  const auto combinable =
      comb ? CombinableLeaves(ceed, *comb, imag, true, num_rhs, count, children)
           : std::vector<bool>(count, false);
  for (CeedInt j = 0; j < count; j++)
  {
    if (const auto *em = ea ? ea->Find(j) : nullptr)
    {
      leaves[j] = *em;
      continue;
    }
    if (IsElementAssemblyLeaf(ceed, children[j], size, num_rhs, true, true))
    {
      std::unique_ptr<ElementMatrices> em;
      if (combinable[j])
      {
        em = CombineElementMatrices(ceed, children[j], *comb, imag);
      }
      if (!em)
      {
        em = AssembleElementMatrices(ceed, children[j], size, true);
      }
      if (em)
      {
        bytes += em->Bytes();
        leaves[j] = std::move(em);
        any_single = true;
      }
    }
  }
  return bytes;
}

// The children of a part without element matrices, as a new finalized composite operator.
std::unique_ptr<Operator>
SinglePrecisionRemainder(Ceed ceed, CeedInt count, CeedOperator *children, int size,
                         const std::vector<std::shared_ptr<ElementMatrices>> &leaves,
                         const std::vector<bool> &dropped)
{
  std::unique_ptr<Operator> remainder;
  for (CeedInt j = 0; j < count; j++)
  {
    if (!leaves[j] && !dropped[j])
    {
      if (!remainder)
      {
        remainder = std::make_unique<Operator>(size, size);
      }
      CeedOperator copy = nullptr;
      PalaceCeedCall(ceed, CeedOperatorReferenceCopy(children[j], &copy));
      remainder->AddSubOperator(copy);
    }
  }
  if (remainder)
  {
    remainder->Finalize();
  }
  return remainder;
}

// Single-precision approximation of a real operator (see CreateSinglePrecisionOperator).
class SinglePrecisionRealOperator final : public palace::Operator
{
private:
  std::vector<std::shared_ptr<ElementMatrices>> leaves;
  std::unique_ptr<Operator> remainder;
  mutable Vector work, temp;

public:
  SinglePrecisionRealOperator(int size,
                              std::vector<std::shared_ptr<ElementMatrices>> &&leaves,
                              std::unique_ptr<Operator> &&remainder)
    : palace::Operator(size), leaves(std::move(leaves)), remainder(std::move(remainder))
  {
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    if (remainder)
    {
      remainder->Mult(x, y);
    }
    else
    {
      y = 0.0;
    }
    const CeedScalar *xd = x.HostRead();
    CeedScalar *yd = y.HostReadWrite();
    for (const auto &em : leaves)
    {
      AddMultElementMatrices(*em, xd, yd, work);
    }
  }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    temp.SetSize(height);
    Mult(x, temp);
    y.Add(a, temp);
  }

  bool SupportsMapped() const { return !remainder; }

  void AddMultMapped(const LocalDofMap &map, const CeedScalar *xt, const CeedScalar *xx,
                     CeedScalar *yt, CeedScalar *yx) const
  {
    for (const auto &em : leaves)
    {
      AddMultElementMatricesMapped(*em, map, xt, xx, yt, yx, work);
    }
  }
};

// Single-precision approximation of a complex operator (see CreateSinglePrecisionOperator):
// terms with the element matrices of a real leaf, of an imaginary leaf, or of both over the
// same elements, and the leaves applied with their own action.
class SinglePrecisionComplexOperator final : public ComplexOperator
{
public:
  struct Term
  {
    std::shared_ptr<ElementMatrices> real, imag;
  };

private:
  std::vector<Term> terms;
  std::unique_ptr<ComplexWrapperOperator> remainder;
  mutable Vector work;
  mutable ComplexVector temp;

public:
  SinglePrecisionComplexOperator(int size, std::vector<Term> &&terms,
                                 std::unique_ptr<ComplexWrapperOperator> &&remainder)
    : ComplexOperator(size), terms(std::move(terms)), remainder(std::move(remainder))
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    if (remainder)
    {
      remainder->Mult(x, y);
    }
    else
    {
      y = 0.0;
    }
    const CeedScalar *xr = x.Real().HostRead(), *xi = x.Imag().HostRead();
    CeedScalar *yr = y.Real().HostReadWrite(), *yi = y.Imag().HostReadWrite();
    for (const auto &[real, imag] : terms)
    {
      AddMultElementPair(real.get(), imag.get(), xr, xi, yr, yi, work);
    }
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    temp.SetSize(height);
    Mult(x, temp);
    y.AXPY(a, temp);
  }

  bool SupportsMapped() const { return !remainder; }

  void AddMultMapped(const LocalDofMap &map, const CeedScalar *xtr, const CeedScalar *xti,
                     const CeedScalar *xxr, const CeedScalar *xxi, CeedScalar *ytr,
                     CeedScalar *yti, CeedScalar *yxr, CeedScalar *yxi) const
  {
    for (const auto &[real, imag] : terms)
    {
      AddMultElementPairMapped(real.get(), imag.get(), map, xtr, xti, xxr, xxi, ytr, yti,
                               yxr, yxi, work);
    }
  }
};

}  // namespace

bool SupportsMapped(const palace::Operator &op)
{
  if (const auto *sop = dynamic_cast<const SinglePrecisionRealOperator *>(&op))
  {
    return sop->SupportsMapped();
  }
  if (const auto *cop = dynamic_cast<const Operator *>(&op))
  {
    const auto *ea = cop->GetElementAssembly();
    return ea && !ea->remainder && !cop->HasDofMultiplicity();
  }
  return false;
}

bool SupportsMapped(const ComplexOperator &op)
{
  if (const auto *sop = dynamic_cast<const SinglePrecisionComplexOperator *>(&op))
  {
    return sop->SupportsMapped();
  }
  if (const auto *pop = dynamic_cast<const PackedComplexOperator *>(&op))
  {
    return pop->SupportsMapped();
  }
  return false;
}

void AddMultMapped(const palace::Operator &op, const LocalDofMap &map, const CeedScalar *xt,
                   const CeedScalar *xx, CeedScalar *yt, CeedScalar *yx)
{
  if (const auto *sop = dynamic_cast<const SinglePrecisionRealOperator *>(&op))
  {
    sop->AddMultMapped(map, xt, xx, yt, yx);
    return;
  }
  const auto *cop = dynamic_cast<const Operator *>(&op);
  const auto *ea = cop ? cop->GetElementAssembly() : nullptr;
  MFEM_VERIFY(ea, "Unsupported operator for AddMultMapped!");
  for (const auto &em : ea->leaves)
  {
    AddMultElementMatricesMapped(*em, map, xt, xx, yt, yx, ea->work);
  }
}

void AddMultMapped(const ComplexOperator &op, const LocalDofMap &map, const CeedScalar *xtr,
                   const CeedScalar *xti, const CeedScalar *xxr, const CeedScalar *xxi,
                   CeedScalar *ytr, CeedScalar *yti, CeedScalar *yxr, CeedScalar *yxi)
{
  if (const auto *sop = dynamic_cast<const SinglePrecisionComplexOperator *>(&op))
  {
    sop->AddMultMapped(map, xtr, xti, xxr, xxi, ytr, yti, yxr, yxi);
    return;
  }
  const auto *pop = dynamic_cast<const PackedComplexOperator *>(&op);
  MFEM_VERIFY(pop, "Unsupported operator for AddMultMapped!");
  pop->AddMultMapped(map, xtr, xti, xxr, xxi, ytr, yti, yxr, yxi);
}

ElementMatrixTerms::ElementMatrixTerms(bool approximation, int num_rhs)
  : approximation(approximation), num_rhs(num_rhs)
{
}

ElementMatrixTerms::~ElementMatrixTerms()
{
  for (auto &leaf : leaves)
  {
    PalaceCeedCallBackend(CeedElemRestrictionDestroy(&leaf.restr));
  }
}

void ElementMatrixTerms::AddTerm(const palace::Operator *Tr, const palace::Operator *Ti)
{
  const int term = num_terms++;
  if (!valid || (!Tr && !Ti))
  {
    return;
  }
  const std::array<const Operator *, 2> parts = {dynamic_cast<const Operator *>(Tr),
                                                 dynamic_cast<const Operator *>(Ti)};
  int size;
  std::array<CeedInt, 2> count;
  std::array<CeedOperator *, 2> children;
  if ((Tr && !parts[0]) || (Ti && !parts[1]) ||
      !SinglePrecisionParts(parts, size, count, children))
  {
    valid = false;
    return;
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  for (int k = 0; k < 2; k++)
  {
    for (CeedInt j = 0; j < count[k]; j++)
    {
      ScopedRestriction in, out;
      PalaceCeedCall(ceed, CeedOperatorGetActiveElemRestrictions(children[k][j], &in.value,
                                                                 &out.value));
      if (!in.value || in.value != out.value)
      {
        valid = false;
        return;
      }
      Leaf leaf;
      leaf.term = term;
      leaf.part = k;
      PalaceCeedCall(ceed, CeedElemRestrictionReferenceCopy(in.value, &leaf.restr));
      // The leaves which the approximations of the sums (see SinglePrecisionLeaves), or
      // the sums themselves (see Operator::UseElementAssembly), may element-assemble.
      if (IsElementAssemblyLeaf(ceed, children[k][j], size, num_rhs, approximation,
                                approximation, true))
      {
        auto em = AssembleElementMatrices(ceed, children[k][j], size, true);
        if (em && em->symmetric)
        {
          // The element matrices of the terms are only combined, never refreshed: release
          // the QData of the term.
          ReleaseElementMatricesLeaf(*em);
          leaf.em = std::move(em);
        }
      }
      leaves.push_back(std::move(leaf));
    }
  }
}

std::size_t ElementMatrixTerms::Bytes() const
{
  std::size_t bytes = 0;
  for (const auto &leaf : leaves)
  {
    bytes += leaf.em ? leaf.em->Bytes() : 0;
  }
  return bytes;
}

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
  const auto combinable = ea_single ? CombinableLeaves(ceed, ea_comb, ea_comb_imag, false,
                                                       num_rhs, count, children)
                                    : std::vector<bool>(count, false);
  for (CeedInt j = 0; j < count; j++)
  {
    std::unique_ptr<ElementMatrices> em;
    // Single-precision storage keeps the choice of the leaves: for small elements, partial
    // assembly of a real and imaginary pair in one pass beats two element-assembled passes
    // even with the smaller element matrices.
    if (IsElementAssemblyLeaf(ceed, children[j], width, num_rhs))
    {
      if (combinable[j])
      {
        em = CombineElementMatrices(ceed, children[j], ea_comb, ea_comb_imag);
      }
      if (!em)
      {
        em = AssembleElementMatrices(ceed, children[j], width, ea_single);
      }
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
  data->work.SetSize(0);
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
    AddMultElementMatrices(*em, xd, yd, ea->work);
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

std::unique_ptr<palace::Operator> CreateSinglePrecisionOperator(const palace::Operator &A,
                                                                std::size_t *bytes)
{
  if (bytes)
  {
    *bytes = 0;
  }
  int size;
  std::array<CeedInt, 2> count;
  std::array<CeedOperator *, 2> children;
  const auto *op = dynamic_cast<const Operator *>(&A);
  if (!op || !SinglePrecisionParts({op, nullptr}, size, count, children))
  {
    return {};
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  bool any_single = false;
  std::vector<std::shared_ptr<ElementMatrices>> leaves;
  const std::size_t new_bytes =
      SinglePrecisionLeaves(ceed, *op, count[0], children[0], size, 1, leaves, any_single);
  if (!any_single)
  {
    return {};
  }
  auto remainder = SinglePrecisionRemainder(ceed, count[0], children[0], size, leaves,
                                            std::vector<bool>(count[0], false));
  std::vector<std::shared_ptr<ElementMatrices>> terms;
  for (auto &em : leaves)
  {
    if (em)
    {
      terms.push_back(std::move(em));
    }
  }
  if (bytes)
  {
    *bytes = new_bytes;
  }
  return std::make_unique<SinglePrecisionRealOperator>(size, std::move(terms),
                                                       std::move(remainder));
}

std::unique_ptr<ComplexOperator>
CreateSinglePrecisionOperator(const palace::Operator *Ar, const palace::Operator *Ai,
                              const FiniteElementSpace &fespace, std::size_t *bytes,
                              int *dropped, const ElementMatrixCombination *comb)
{
  if (bytes)
  {
    *bytes = 0;
  }
  if (dropped)
  {
    *dropped = 0;
  }
  const std::array<const Operator *, 2> parts = {dynamic_cast<const Operator *>(Ar),
                                                 dynamic_cast<const Operator *>(Ai)};
  int size;
  std::array<CeedInt, 2> count;
  std::array<CeedOperator *, 2> children;
  if ((Ar && !parts[0]) || (Ai && !parts[1]) ||
      !SinglePrecisionParts(parts, size, count, children))
  {
    return {};
  }
  Ceed ceed = internal::GetCeedObjects()[0];
  bool any_single = false;
  std::size_t new_bytes = 0;
  std::array<std::vector<std::shared_ptr<ElementMatrices>>, 2> leaves;
  for (int k = 0; k < 2; k++)
  {
    if (parts[k])
    {
      new_bytes += SinglePrecisionLeaves(ceed, *parts[k], count[k], children[k], size, 2,
                                         leaves[k], any_single, comb, k == 1);
    }
  }
  if (!any_single)
  {
    return {};
  }

  // Leave out the imaginary volume leaves which are negligible in single precision
  // relative to the element-assembled real part, row by row. Boundary terms are kept.
  constexpr double drop_tol = 1.0e-6;
  std::array<std::vector<bool>, 2> skip = {std::vector<bool>(count[0], false),
                                           std::vector<bool>(count[1], false)};
  int num_dropped = 0;
  if (parts[0] && parts[1])
  {
    auto LeafDimension = [&](CeedOperator leaf)
    {
      ScopedBasis basis_in, basis_out;
      PalaceCeedCall(ceed,
                     CeedOperatorGetActiveBases(leaf, &basis_in.value, &basis_out.value));
      CeedInt dim = 0;
      if (basis_in.value && basis_in.value != CEED_BASIS_NONE)
      {
        PalaceCeedCall(ceed, CeedBasisGetDimension(basis_in.value, &dim));
      }
      return dim;
    };
    CeedInt volume_dim = 0;
    for (int k = 0; k < 2; k++)
    {
      for (CeedInt j = 0; j < count[k]; j++)
      {
        volume_dim = std::max(volume_dim, LeafDimension(children[k][j]));
      }
    }
    std::vector<double> scale(size, 0.0);
    for (const auto &em : leaves[0])
    {
      if (em)
      {
        ElementRowMaxima(*em, scale, 0.0);
      }
    }
    for (CeedInt j = 0; j < count[1]; j++)
    {
      auto &em = leaves[1][j];
      if (em && LeafDimension(children[1][j]) == volume_dim &&
          ElementRowMaxima(*em, scale, drop_tol))
      {
        if (em->single && em.use_count() == 1)
        {
          new_bytes -= em->Bytes();
        }
        em.reset();
        skip[1][j] = true;
        num_dropped++;
      }
    }
  }

  // The children without element matrices keep their own action.
  std::array<std::unique_ptr<palace::Operator>, 2> remainder_parts;
  for (int k = 0; k < 2; k++)
  {
    if (parts[k])
    {
      remainder_parts[k] =
          SinglePrecisionRemainder(ceed, count[k], children[k], size, leaves[k], skip[k]);
    }
  }
  std::unique_ptr<ComplexWrapperOperator> remainder;
  if (remainder_parts[0] || remainder_parts[1])
  {
    remainder = CreateComplexOperator(std::move(remainder_parts[0]),
                                      std::move(remainder_parts[1]), fespace);
  }

  // Pair the single-precision real and imaginary leaves over the same elements, so that
  // each element is gathered and scattered once for both.
  std::vector<SinglePrecisionComplexOperator::Term> terms;
  std::vector<bool> paired(count[1], false);
  for (auto &real : leaves[0])
  {
    if (!real)
    {
      continue;
    }
    SinglePrecisionComplexOperator::Term term{std::move(real), nullptr};
    for (CeedInt j = 0; term.real->single && j < count[1]; j++)
    {
      auto &imag = leaves[1][j];
      if (!paired[j] && imag && imag->single && imag.use_count() == 1 &&
          imag->elements == term.real->elements && imag->nodes == term.real->nodes &&
          imag->indices == term.real->indices)
      {
        // The real leaf's indices serve both.
        new_bytes -= imag->indices.size() * sizeof(CeedInt);
        std::vector<CeedInt>().swap(imag->indices);
        term.imag = imag;
        paired[j] = true;
        break;
      }
    }
    terms.push_back(std::move(term));
  }
  for (CeedInt j = 0; j < count[1]; j++)
  {
    if (leaves[1][j] && !paired[j])
    {
      terms.push_back({nullptr, std::move(leaves[1][j])});
    }
  }

  if (bytes)
  {
    *bytes = new_bytes;
  }
  if (dropped)
  {
    *dropped = num_dropped;
  }
  return std::make_unique<SinglePrecisionComplexOperator>(size, std::move(terms),
                                                          std::move(remainder));
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

namespace
{

// Whether AssembleExactDiagonal supports a leaf: a non-composite operator, not at points,
// with one active input and one active output basis sharing one restriction with standard,
// oriented or curl-oriented transforms. Empty leaves are supported.
bool ExactDiagonalLeaf(Ceed ceed, CeedOperator leaf)
{
  bool composite, at_points;
  PalaceCeedCall(ceed, CeedOperatorIsComposite(leaf, &composite));
  if (composite)
  {
    return false;
  }
  PalaceCeedCall(ceed, CeedOperatorIsAtPoints(leaf, &at_points));
  if (at_points)
  {
    return false;
  }
  CeedInt num_elem;
  PalaceCeedCall(ceed, CeedOperatorGetNumElements(leaf, &num_elem));
  if (num_elem == 0)
  {
    return true;
  }
  CeedOperatorAssemblyData data;
  PalaceCeedCall(ceed, CeedOperatorGetOperatorAssemblyData(leaf, &data));
  CeedInt num_bases_in, num_bases_out, *num_modes_in, *num_modes_out;
  PalaceCeedCall(ceed, CeedOperatorAssemblyDataGetEvalModes(
                           data, &num_bases_in, &num_modes_in, nullptr, nullptr,
                           &num_bases_out, &num_modes_out, nullptr, nullptr, nullptr));
  if (num_bases_in != 1 || num_bases_out != 1 || num_modes_in[0] <= 0 ||
      num_modes_out[0] <= 0)
  {
    return false;
  }
  CeedElemRestriction *rstrs_in, *rstrs_out;
  PalaceCeedCall(ceed, CeedOperatorAssemblyDataGetElemRestrictions(data, nullptr, &rstrs_in,
                                                                   nullptr, &rstrs_out));
  if (rstrs_in[0] != rstrs_out[0])
  {
    return false;
  }
  CeedRestrictionType type;
  CeedInt rstr_elem;
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(rstrs_in[0], &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumElements(rstrs_in[0], &rstr_elem));
  return (rstr_elem == num_elem &&
          (type == CEED_RESTRICTION_STANDARD || type == CEED_RESTRICTION_ORIENTED ||
           type == CEED_RESTRICTION_CURL_ORIENTED));
}

// Add the diagonal of a supported leaf in L-vector coordinates to diag (if exact_local),
// and pass its entries between distinct coupled L-dofs to add_entry. The element matrices
// are formed from the assembled QFunction like CeedOperatorLinearAssemble: A_e = Bᵀ D_e B
// in the element basis, then Tᵀ A_e T with the orientation transforms T of the restriction.
// Only the entries needed are computed: the diagonal of Tᵀ A_e T (with tridiagonal curl
// orientation transforms, the entries of A_e between neighboring nodes), and, for the
// elements with at least two coupled L-dofs, its entries between those. With all_entries,
// the entries (including diagonal ones) between all marked L-dofs of every element with an
// L-dof marked 2 or more are passed instead (see AssembleMarkedEntries).
void AddExactLeafDiagonal(Ceed ceed, CeedOperator leaf, const char *coupled,
                          bool exact_local, CeedScalar *diag,
                          const std::function<void(int, int, double)> &add_entry,
                          bool all_entries = false)
{
  CeedInt num_elem;
  PalaceCeedCall(ceed, CeedOperatorGetNumElements(leaf, &num_elem));
  if (num_elem == 0)
  {
    return;
  }

  // Active restriction, its transforms and the element structure. With an approximate local
  // diagonal, only elements with coupled L-dofs need any entries.
  CeedOperatorAssemblyData data;
  PalaceCeedCall(ceed, CeedOperatorGetOperatorAssemblyData(leaf, &data));
  CeedElemRestriction *rstrs_in;
  PalaceCeedCall(ceed, CeedOperatorAssemblyDataGetElemRestrictions(data, nullptr, &rstrs_in,
                                                                   nullptr, nullptr));
  CeedElemRestriction rstr = rstrs_in[0];
  CeedRestrictionType type;
  CeedInt P, num_comp, comp_stride;
  PalaceCeedCall(ceed, CeedElemRestrictionGetType(rstr, &type));
  PalaceCeedCall(ceed, CeedElemRestrictionGetElementSize(rstr, &P));
  PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(rstr, &num_comp));
  PalaceCeedCall(ceed, CeedElemRestrictionGetCompStride(rstr, &comp_stride));
  const CeedInt *offsets;
  const bool *orients = nullptr;
  const CeedInt8 *curl_orients = nullptr;
  PalaceCeedCall(ceed, CeedElemRestrictionGetOffsets(rstr, CEED_MEM_HOST, &offsets));
  if (type == CEED_RESTRICTION_ORIENTED)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionGetOrientations(rstr, CEED_MEM_HOST, &orients));
  }
  else if (type == CEED_RESTRICTION_CURL_ORIENTED)
  {
    PalaceCeedCall(
        ceed, CeedElemRestrictionGetCurlOrientations(rstr, CEED_MEM_HOST, &curl_orients));
  }
  auto Index = [&](CeedInt e, CeedInt i, CeedInt c)
  { return offsets[static_cast<std::size_t>(e) * P + i] + c * comp_stride; };
  std::vector<char> full(num_elem, 0);
  bool any_full = false;
  if (coupled)
  {
    for (CeedInt e = 0; e < num_elem; e++)
    {
      int count = 0, count_shared = 0;
      for (CeedInt i = 0; i < P; i++)
      {
        for (CeedInt c = 0; c < num_comp; c++)
        {
          const char m = coupled[Index(e, i, c)];
          count += (m != 0);
          count_shared += (m > 1);
        }
      }
      full[e] = all_entries ? (count_shared > 0) : (count >= 2 && count_shared > 0);
      any_full = any_full || full[e];
    }
  }
  if (exact_local || any_full)
  {
    // Assembled QFunction, eval modes and basis matrices (B[(q * E + m) * P + n]).
    CeedVector qf_vec = nullptr;
    CeedElemRestriction qf_rstr = nullptr;
    PalaceCeedCall(ceed, CeedOperatorLinearAssembleQFunctionBuildOrUpdate(
                             leaf, &qf_vec, &qf_rstr, CEED_REQUEST_IMMEDIATE));
    CeedInt layout[3];
    PalaceCeedCall(ceed, CeedElemRestrictionGetELayout(qf_rstr, layout));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&qf_rstr));
    const CeedScalar *qf;
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(qf_vec, CEED_MEM_HOST, &qf));
    CeedInt *num_modes_in, *num_modes_out;
    CeedSize **mode_offsets_in, **mode_offsets_out, num_out_comp;
    PalaceCeedCall(ceed,
                   CeedOperatorAssemblyDataGetEvalModes(
                       data, nullptr, &num_modes_in, nullptr, &mode_offsets_in, nullptr,
                       &num_modes_out, nullptr, &mode_offsets_out, &num_out_comp));
    CeedBasis *bases_in;
    const CeedScalar **B_ins, **B_outs;
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataGetBases(data, nullptr, &bases_in, &B_ins,
                                                          nullptr, nullptr, &B_outs));
    CeedInt Q = P;
    if (bases_in[0] != CEED_BASIS_NONE)
    {
      PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(bases_in[0], &Q));
    }
    const CeedInt E_in = num_modes_in[0], E_out = num_modes_out[0], K = Q * E_in;
    const CeedScalar *B_in = B_ins[0], *B_out = B_outs[0];

    // Bᵀ (nodes by points and input modes) for contiguous products, and work arrays.
    std::vector<CeedScalar> Bt(static_cast<std::size_t>(P) * K), BtD(Bt.size()),
        D(static_cast<std::size_t>(E_out) * E_in), A;
    for (CeedInt k = 0; k < K; k++)
    {
      for (CeedInt n = 0; n < P; n++)
      {
        Bt[static_cast<std::size_t>(n) * K + k] = B_in[static_cast<std::size_t>(k) * P + n];
      }
    }
    auto Entry = [&](CeedInt a, CeedInt b)
    {
      const CeedScalar *x = BtD.data() + static_cast<std::size_t>(a) * K,
                       *y = Bt.data() + static_cast<std::size_t>(b) * K;
      CeedScalar sum = 0.0;
      for (CeedInt k = 0; k < K; k++)
      {
        sum += x[k] * y[k];
      }
      return sum;
    };
    // Coefficients of the transform of node i: (node, weight) pairs of Tᵀ row i.
    auto Transform = [&](CeedInt e, CeedInt i, CeedInt(&nodes)[3], CeedScalar(&w)[3])
    {
      int count = 0;
      if (curl_orients)
      {
        const CeedInt8 *co = curl_orients + 3 * (static_cast<std::size_t>(e) * P);
        if (i > 0 && co[3 * i - 1] != 0)
        {
          nodes[count] = i - 1;
          w[count++] = co[3 * i - 1];
        }
        if (co[3 * i + 1] != 0)
        {
          nodes[count] = i;
          w[count++] = co[3 * i + 1];
        }
        if (i < P - 1 && co[3 * i + 3] != 0)
        {
          nodes[count] = i + 1;
          w[count++] = co[3 * i + 3];
        }
      }
      else
      {
        nodes[count] = i;
        w[count++] = (orients && orients[static_cast<std::size_t>(e) * P + i]) ? -1.0 : 1.0;
      }
      return count;
    };
    // For the entries between coupled L-dofs: the coupled nodes of the element for the
    // output and input components, the nodes of A_e their transforms involve, and the block
    // of A_e between those.
    std::vector<CeedInt> coupled_out, coupled_in, needed_out, needed_in, pos_out(P),
        pos_in(P);
    auto Needed = [&](CeedInt e, CeedInt comp, std::vector<CeedInt> &nodes_c,
                      std::vector<CeedInt> &needed, std::vector<CeedInt> &pos)
    {
      nodes_c.clear();
      needed.clear();
      std::fill(pos.begin(), pos.end(), -1);
      for (CeedInt i = 0; i < P; i++)
      {
        if (!coupled[Index(e, i, comp)])
        {
          continue;
        }
        nodes_c.push_back(i);
        CeedInt nodes[3];
        CeedScalar w[3];
        const int count = Transform(e, i, nodes, w);
        for (int a = 0; a < count; a++)
        {
          if (pos[nodes[a]] < 0)
          {
            pos[nodes[a]] = static_cast<CeedInt>(needed.size());
            needed.push_back(nodes[a]);
          }
        }
      }
    };
    for (CeedInt e = 0; e < num_elem; e++)
    {
      if (!exact_local && !full[e])
      {
        continue;
      }
      for (CeedInt comp_in = 0; comp_in < num_comp; comp_in++)
      {
        for (CeedInt comp_out = 0; comp_out < num_comp; comp_out++)
        {
          const bool diagonal = exact_local && comp_in == comp_out;
          if (!diagonal && !full[e])
          {
            continue;
          }
          // BtD[n][q * E_in + m] = Σ_m' B_out[(q * E_out + m') * P + n] D_q[m'][m].
          std::fill(BtD.begin(), BtD.end(), 0.0);
          for (CeedInt q = 0; q < Q; q++)
          {
            for (CeedInt mo = 0; mo < E_out; mo++)
            {
              for (CeedInt mi = 0; mi < E_in; mi++)
              {
                const CeedSize index = (mode_offsets_in[0][mi] + comp_in) * num_out_comp +
                                       mode_offsets_out[0][mo] + comp_out;
                D[mo * E_in + mi] = qf[q * layout[0] + index * layout[1] +
                                       static_cast<CeedSize>(e) * layout[2]];
              }
            }
            for (CeedInt n = 0; n < P; n++)
            {
              CeedScalar *row = BtD.data() + static_cast<std::size_t>(n) * K + q * E_in;
              for (CeedInt mo = 0; mo < E_out; mo++)
              {
                const CeedScalar b =
                    B_out[(static_cast<std::size_t>(q) * E_out + mo) * P + n];
                if (b == 0.0)
                {
                  continue;
                }
                for (CeedInt mi = 0; mi < E_in; mi++)
                {
                  row[mi] += b * D[mo * E_in + mi];
                }
              }
            }
          }
          if (diagonal)
          {
            // Diagonal of Tᵀ A_e T.
            for (CeedInt i = 0; i < P; i++)
            {
              CeedInt nodes[3];
              CeedScalar w[3];
              const int count = Transform(e, i, nodes, w);
              CeedScalar sum = 0.0;
              for (int a = 0; a < count; a++)
              {
                for (int b = 0; b < count; b++)
                {
                  sum += w[a] * w[b] * Entry(nodes[a], nodes[b]);
                }
              }
              diag[Index(e, i, comp_out)] += sum;
            }
          }
          if (!full[e])
          {
            continue;
          }
          // Entries (Tᵀ A_e T)_ij between distinct coupled L-dofs.
          Needed(e, comp_out, coupled_out, needed_out, pos_out);
          Needed(e, comp_in, coupled_in, needed_in, pos_in);
          const std::size_t n_out = needed_out.size(), n_in = needed_in.size();
          A.resize(n_out * n_in);
          for (std::size_t a = 0; a < n_out; a++)
          {
            for (std::size_t b = 0; b < n_in; b++)
            {
              A[a * n_in + b] = Entry(needed_out[a], needed_in[b]);
            }
          }
          for (const CeedInt i : coupled_out)
          {
            const CeedInt r = Index(e, i, comp_out);
            CeedInt nodes_i[3];
            CeedScalar w_i[3];
            const int count_i = Transform(e, i, nodes_i, w_i);
            for (const CeedInt j : coupled_in)
            {
              const CeedInt c = Index(e, j, comp_in);
              if (r == c && !all_entries)
              {
                continue;
              }
              CeedInt nodes_j[3];
              CeedScalar w_j[3];
              const int count_j = Transform(e, j, nodes_j, w_j);
              CeedScalar sum = 0.0;
              for (int a = 0; a < count_i; a++)
              {
                for (int b = 0; b < count_j; b++)
                {
                  sum +=
                      w_i[a] * w_j[b] * A[pos_out[nodes_i[a]] * n_in + pos_in[nodes_j[b]]];
                }
              }
              add_entry(r, c, sum);
            }
          }
        }
      }
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(qf_vec, &qf));
    PalaceCeedCall(ceed, CeedVectorDestroy(&qf_vec));
  }
  if (orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreOrientations(rstr, &orients));
  }
  if (curl_orients)
  {
    PalaceCeedCall(ceed, CeedElemRestrictionRestoreCurlOrientations(rstr, &curl_orients));
  }
  PalaceCeedCall(ceed, CeedElemRestrictionRestoreOffsets(rstr, &offsets));
}

// Whether AssembleExactDiagonal supports all leaves of an operator.
bool ExactDiagonalSupported(const Operator &op)
{
  bool supported = true;
  for (std::size_t id = 0; id < op.Size() && supported; id++)
  {
    Ceed ceed = CeedOperatorReturnCeed(op[id]);
    CeedInt count;
    CeedOperator *children;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &count));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[id], &children));
    for (CeedInt k = 0; k < count && supported; k++)
    {
      supported = ExactDiagonalLeaf(ceed, children[k]);
    }
    PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(op[id]));
  }
  return supported;
}

}  // namespace

bool AssembleExactDiagonal(const Operator &op, const std::vector<char> *coupled,
                           bool exact_local, Vector &diag,
                           const std::function<void(int, int, double)> &add_entry)
{
  MFEM_VERIFY(op.Height() == op.Width(),
              "Exact diagonal assembly requires a square operator!");
  MFEM_VERIFY(!coupled || static_cast<int>(coupled->size()) == op.Height(),
              "Invalid size for the coupled L-dof marker!");
  // Check all leaves first, so that an unsupported operator leaves no partial result.
  if (!ExactDiagonalSupported(op))
  {
    return false;
  }
  diag.SetSize(op.Height());
  if (exact_local)
  {
    diag = 0.0;
  }
  else
  {
    // The approximate diagonal of libCEED (absolute values of the orientation transforms).
    op.AssembleDiagonal(diag);
  }
  CeedScalar *d = diag.HostReadWrite();
  const char *marker = coupled ? coupled->data() : nullptr;
  // The leaves of all threads are processed by the calling thread, which accumulates into
  // the same diagonal.
  for (std::size_t id = 0; id < op.Size(); id++)
  {
    Ceed ceed = CeedOperatorReturnCeed(op[id]);
    CeedInt count;
    CeedOperator *children;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &count));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[id], &children));
    for (CeedInt k = 0; k < count; k++)
    {
      AddExactLeafDiagonal(ceed, children[k], marker, exact_local, d, add_entry);
      PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(children[k]));
    }
  }
  return true;
}

bool AssembleMarkedEntries(const Operator &op, const std::vector<char> &marker,
                           const std::function<void(int, int, double)> &add_entry)
{
  MFEM_VERIFY(op.Height() == op.Width(),
              "Marked entry assembly requires a square operator!");
  MFEM_VERIFY(static_cast<int>(marker.size()) == op.Height(),
              "Invalid size for the L-dof marker!");
  if (!ExactDiagonalSupported(op))
  {
    return false;
  }
  for (std::size_t id = 0; id < op.Size(); id++)
  {
    Ceed ceed = CeedOperatorReturnCeed(op[id]);
    CeedInt count;
    CeedOperator *children;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[id], &count));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[id], &children));
    for (CeedInt k = 0; k < count; k++)
    {
      AddExactLeafDiagonal(ceed, children[k], marker.data(), false, nullptr, add_entry,
                           true);
      PalaceCeedCall(ceed, CeedOperatorAssemblyDataStrip(children[k]));
    }
  }
  return true;
}

}  // namespace palace::ceed
