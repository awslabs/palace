// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basis.hpp"

#include <cstdint>
#include <cstring>
#include <map>
#include <mutex>
#include <vector>
#include <mfem.hpp>
#include "utils/diagnostic.hpp"

namespace palace::ceed
{

namespace
{

using BasisKey = std::vector<std::uint64_t>;
using BasisCache = std::map<BasisKey, CeedBasis>;

BasisCache &GetBasisCache()
{
  static auto *cache = new BasisCache();
  return *cache;
}

std::mutex &GetBasisCacheMutex()
{
  static auto *cache_mutex = new std::mutex();
  return *cache_mutex;
}

void InitTensorBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                     CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  // The x-coordinates of the first `Q` points of the integration rule are the points of
  // the corresponding 1D rule. We also scale the weights accordingly.
  const mfem::DofToQuad &maps = fe.GetDofToQuad(ir, mfem::DofToQuad::TENSOR);
  const int dim = fe.GetDim();
  const int P = maps.ndof;
  const int Q = maps.nqpt;
  mfem::Vector qX(Q), qW(Q);
  double w_sum = 0.0;
  for (int i = 0; i < Q; i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    qX(i) = ip.x;
    qW(i) = ip.weight;
    w_sum += ip.weight;
  }
  qW *= 1.0 / w_sum;

  PalaceCeedCall(ceed, CeedBasisCreateTensorH1(ceed, dim, num_comp, P, Q, maps.Bt.GetData(),
                                               maps.Gt.GetData(), qX.GetData(),
                                               qW.GetData(), basis));
}

void InitNonTensorBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                        CeedInt num_comp, Ceed ceed, CeedBasis *basis, bool cache_fixed)
{
  const mfem::DofToQuad &maps = fe.GetDofToQuad(ir, mfem::DofToQuad::FULL);
  const int dim = fe.GetDim();
  const int P = maps.ndof;
  const int Q = maps.nqpt;
  // libCEED expects reference coordinates in component-major layout,
  // q_ref[d * Q + q].  MFEM DenseMatrix is column-major, so store the logical
  // transpose (points as rows, dimensions as columns): qX(q, d) has exactly the
  // layout libCEED expects while keeping matrix-style indexing.
  mfem::DenseMatrix qX(Q, dim);
  mfem::Vector qW(Q);
  for (int i = 0; i < Q; i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    qX(i, 0) = ip.x;
    if (dim > 1)
    {
      qX(i, 1) = ip.y;
    }
    if (dim > 2)
    {
      qX(i, 2) = ip.z;
    }
    qW(i) = ip.weight;
  }

  auto Create = [&](CeedBasis *created)
  {
    if (fe.GetMapType() == mfem::FiniteElement::H_DIV)
    {
      PalaceCeedCall(ceed,
                     CeedBasisCreateHdiv(ceed, GetCeedTopology(fe.GetGeomType()), num_comp,
                                         P, Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                         qX.GetData(), qW.GetData(), created));
    }
    else if (fe.GetMapType() == mfem::FiniteElement::H_CURL)
    {
      PalaceCeedCall(ceed,
                     CeedBasisCreateHcurl(ceed, GetCeedTopology(fe.GetGeomType()), num_comp,
                                          P, Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                          qX.GetData(), qW.GetData(), created));
    }
    else
    {
      PalaceCeedCall(ceed,
                     CeedBasisCreateH1(ceed, GetCeedTopology(fe.GetGeomType()), num_comp, P,
                                       Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                       qX.GetData(), qW.GetData(), created));
    }
  };
  if (!cache_fixed)
  {
    Create(basis);
    return;
  }

  // Fixed face/domain point clouds recur across output fields. Reuse the immutable
  // ordinary libCEED basis so GPU backend setup and tabulation storage are paid once per
  // exact FE/rule pair. The key contains the full MFEM tabulation, not FE object
  // addresses, so finite-element collection destruction/address reuse cannot alias it.
  BasisKey key;
  key.reserve(16 + maps.Bt.Size() + maps.Gt.Size() + qX.Height() * qX.Width() + qW.Size());
  auto AppendInt = [&](std::uint64_t value) { key.push_back(value); };
  auto AppendDouble = [&](double value)
  {
    std::uint64_t bits;
    static_assert(sizeof(bits) == sizeof(value));
    std::memcpy(&bits, &value, sizeof(bits));
    key.push_back(bits);
  };
  auto AppendArray = [&](const auto &array)
  {
    AppendInt(static_cast<std::uint64_t>(array.Size()));
    for (int i = 0; i < array.Size(); i++)
    {
      AppendDouble(array.GetData()[i]);
    }
  };
  auto AppendMatrix = [&](const mfem::DenseMatrix &mat)
  {
    AppendInt(static_cast<std::uint64_t>(mat.Height()));
    AppendInt(static_cast<std::uint64_t>(mat.Width()));
    for (int i = 0; i < mat.Height() * mat.Width(); i++)
    {
      AppendDouble(mat.GetData()[i]);
    }
  };
  AppendInt(static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(ceed)));
  AppendInt(static_cast<std::uint64_t>(fe.GetGeomType()));
  AppendInt(static_cast<std::uint64_t>(fe.GetMapType()));
  AppendInt(static_cast<std::uint64_t>(fe.GetRangeType()));
  AppendInt(static_cast<std::uint64_t>(dim));
  AppendInt(static_cast<std::uint64_t>(num_comp));
  AppendInt(static_cast<std::uint64_t>(P));
  AppendInt(static_cast<std::uint64_t>(Q));
  AppendArray(maps.Bt);
  AppendArray(maps.Gt);
  AppendMatrix(qX);
  AppendInt(static_cast<std::uint64_t>(qW.Size()));
  for (int i = 0; i < qW.Size(); i++)
  {
    AppendDouble(qW(i));
  }

  // The cache lives until explicit libCEED finalization, matching registered
  // IntegrationRules retained for MFEM's pointer-keyed DofToQuad cache. Stored handles
  // keep one reference; callers receive normal reference-counted copies owned by their
  // operators.
  auto &cache = GetBasisCache();
  std::lock_guard<std::mutex> lock(GetBasisCacheMutex());
  auto it = cache.find(key);
  if (it == cache.end())
  {
    CeedBasis cached_basis;
    Create(&cached_basis);
    it = cache.emplace(std::move(key), cached_basis).first;
  }
  *basis = nullptr;
  PalaceCeedCall(ceed, CeedBasisReferenceCopy(it->second, basis));
}

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

void InitCeedInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                               const mfem::FiniteElement &test_fe, CeedInt trial_num_comp,
                               CeedInt test_num_comp, Ceed ceed, CeedBasis *basis)
{
  // Basis projection operator using libCEED.
  CeedBasis trial_basis, test_basis;
  const int P = std::max(trial_fe.GetDof(), test_fe.GetDof()), ir_order_max = 100;
  int ir_order = std::max(trial_fe.GetOrder(), test_fe.GetOrder());
  for (; ir_order < ir_order_max; ir_order++)
  {
    if (mfem::IntRules.Get(trial_fe.GetGeomType(), ir_order).GetNPoints() >= P)
    {
      break;
    }
  }
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(trial_fe.GetGeomType(), ir_order);

  InitBasis(trial_fe, ir, trial_num_comp, ceed, &trial_basis);
  InitBasis(test_fe, ir, test_num_comp, ceed, &test_basis);
  PalaceCeedCall(ceed, CeedBasisCreateProjection(trial_basis, test_basis, basis));
  PalaceCeedCall(ceed, CeedBasisDestroy(&trial_basis));
  PalaceCeedCall(ceed, CeedBasisDestroy(&test_basis));
}

PalacePragmaDiagnosticPop

void InitMfemInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                               const mfem::FiniteElement &test_fe, CeedInt trial_num_comp,
                               CeedInt test_num_comp, Ceed ceed, CeedBasis *basis)
{
  MFEM_VERIFY(trial_num_comp == test_num_comp && trial_num_comp == 1,
              "libCEED discrete linear operator requires same vdim = 1 for trial and test "
              "FE spaces!");
  const int trial_P = trial_fe.GetDof();
  const int test_P = test_fe.GetDof();
  mfem::DenseMatrix Bt, Gt(trial_P, test_P);
  mfem::Vector qX(test_P), qW(test_P);
  mfem::IsoparametricTransformation dummy;
  dummy.SetIdentityTransformation(trial_fe.GetGeomType());
  if (trial_fe.GetMapType() == test_fe.GetMapType())
  {
    // Prolongation.
    test_fe.GetTransferMatrix(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::VALUE &&
           test_fe.GetMapType() == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator.
    test_fe.ProjectGrad(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::H_CURL &&
           test_fe.GetMapType() == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator.
    test_fe.ProjectCurl(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::H_DIV &&
           test_fe.GetMapType() == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator.
    test_fe.ProjectDiv(trial_fe, dummy, Bt);
  }
  else
  {
    MFEM_ABORT("Unsupported trial/test FE spaces for libCEED discrete linear operator!");
  }
  Bt.Transpose();
  Gt = 0.0;
  qX = 0.0;
  qW = 0.0;

  // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
  PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, trial_num_comp, trial_P,
                                         test_P, Bt.GetData(), Gt.GetData(), qX.GetData(),
                                         qW.GetData(), basis));
}

}  // namespace

namespace internal
{

void FinalizeBasisCache()
{
  auto &cache = GetBasisCache();
  std::lock_guard<std::mutex> lock(GetBasisCacheMutex());
  for (auto &[key, basis] : cache)
  {
    (void)key;
    Ceed ceed = CeedBasisReturnCeed(basis);
    PalaceCeedCall(ceed, CeedBasisDestroy(&basis));
  }
  cache.clear();
}

}  // namespace internal

void InitBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
               CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  if constexpr (false)
  {
    std::cout << "New basis (" << ceed << ", " << &fe << ", " << &ir << ")\n";
  }
  const bool tensor = dynamic_cast<const mfem::TensorBasisElement *>(&fe) != nullptr;
  const bool vector = fe.GetRangeType() == mfem::FiniteElement::VECTOR;
  if (tensor && !vector)
  {
    InitTensorBasis(fe, ir, num_comp, ceed, basis);
  }
  else
  {
    InitNonTensorBasis(fe, ir, num_comp, ceed, basis, false);
  }
}

void InitBasisFromRule(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                       CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  // Fixed mapped rules require FULL tabulation: their points need not form a tensor
  // product. The caller keeps the rule alive because MFEM caches DofToQuad by its
  // pointer inside the shared finite element.
  InitNonTensorBasis(fe, ir, num_comp, ceed, basis, false);
}

void InitCachedBasisFromRule(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                             CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  InitNonTensorBasis(fe, ir, num_comp, ceed, basis, true);
}

void InitInterpolatorBasis(const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, CeedInt trial_num_comp,
                           CeedInt test_num_comp, Ceed ceed, CeedBasis *basis)
{
  if constexpr (false)
  {
    std::cout << "New interpolator basis (" << ceed << ", " << &trial_fe << ", " << &test_fe
              << ", " << (trial_fe.GetMapType() == test_fe.GetMapType()) << ")\n";
  }
  if constexpr (false)
  {
    if (trial_fe.GetMapType() == test_fe.GetMapType())
    {
      InitCeedInterpolatorBasis(trial_fe, test_fe, trial_num_comp, test_num_comp, ceed,
                                basis);
    }
    else
    {
      InitMfemInterpolatorBasis(trial_fe, test_fe, trial_num_comp, test_num_comp, ceed,
                                basis);
    }
  }
  else
  {
    InitMfemInterpolatorBasis(trial_fe, test_fe, trial_num_comp, test_num_comp, ceed,
                              basis);
  }
}

}  // namespace palace::ceed
