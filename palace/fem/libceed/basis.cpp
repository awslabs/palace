// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basis.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <mfem.hpp>
#include "utils/diagnostic.hpp"

namespace palace::ceed
{

namespace
{

int SimplexNumModes(mfem::Geometry::Type geom, int degree)
{
  if (geom == mfem::Geometry::TRIANGLE)
  {
    return (degree + 1) * (degree + 2) / 2;
  }
  MFEM_VERIFY(geom == mfem::Geometry::TETRAHEDRON,
              "AtPoints lattice requested for a non-simplex element!");
  return (degree + 1) * (degree + 2) * (degree + 3) / 6;
}

const mfem::IntegrationRule &GetRegisteredSimplexLatticeRule(mfem::Geometry::Type geom,
                                                             int degree)
{
  using RuleKey = std::pair<int, int>;
  static std::map<RuleKey, std::unique_ptr<mfem::IntegrationRule>> registry;
  static std::mutex registry_mutex;
  const RuleKey key{static_cast<int>(geom), degree};
  std::lock_guard<std::mutex> lock(registry_mutex);
  auto it = registry.find(key);
  if (it == registry.end())
  {
    auto ir = std::make_unique<mfem::IntegrationRule>(SimplexNumModes(geom, degree));
    int q = 0;
    if (degree == 0)
    {
      if (geom == mfem::Geometry::TRIANGLE)
      {
        ir->IntPoint(q).Set2(1.0 / 3.0, 1.0 / 3.0);
      }
      else
      {
        ir->IntPoint(q).Set3(0.25, 0.25, 0.25);
      }
      ir->IntPoint(q++).weight = 1.0;
    }
    else if (geom == mfem::Geometry::TRIANGLE)
    {
      for (int total = 0; total <= degree; total++)
      {
        for (int i = 0; i <= total; i++)
        {
          const int j = total - i;
          ir->IntPoint(q).Set2(static_cast<double>(i) / degree,
                               static_cast<double>(j) / degree);
          ir->IntPoint(q++).weight = 1.0;
        }
      }
    }
    else
    {
      for (int total = 0; total <= degree; total++)
      {
        for (int i = 0; i <= total; i++)
        {
          for (int j = 0; j <= total - i; j++)
          {
            const int k = total - i - j;
            ir->IntPoint(q).Set3(static_cast<double>(i) / degree,
                                 static_cast<double>(j) / degree,
                                 static_cast<double>(k) / degree);
            ir->IntPoint(q++).weight = 1.0;
          }
        }
      }
    }
    MFEM_ASSERT(q == ir->GetNPoints(), "Invalid simplex lattice rule size!");
    it = registry.emplace(key, std::move(ir)).first;
  }
  return *it->second;
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
                        CeedInt num_comp, Ceed ceed, CeedBasis *basis)
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

  if (fe.GetMapType() == mfem::FiniteElement::H_DIV)
  {
    PalaceCeedCall(ceed,
                   CeedBasisCreateHdiv(ceed, GetCeedTopology(fe.GetGeomType()), num_comp, P,
                                       Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                       qX.GetData(), qW.GetData(), basis));
  }
  else if (fe.GetMapType() == mfem::FiniteElement::H_CURL)
  {
    PalaceCeedCall(ceed,
                   CeedBasisCreateHcurl(ceed, GetCeedTopology(fe.GetGeomType()), num_comp,
                                        P, Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                        qX.GetData(), qW.GetData(), basis));
  }
  else
  {
    PalaceCeedCall(ceed,
                   CeedBasisCreateH1(ceed, GetCeedTopology(fe.GetGeomType()), num_comp, P,
                                     Q, maps.Bt.GetData(), maps.Gt.GetData(), qX.GetData(),
                                     qW.GetData(), basis));
  }
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
    InitNonTensorBasis(fe, ir, num_comp, ceed, basis);
  }
}

void InitBasisAtPoints(const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                       CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  // Always use full tabulation: the integration rule points may be arbitrary (not a
  // tensor-product rule), which the tensor basis construction cannot represent. The
  // corresponding element restriction must use native dof ordering. The caller must
  // keep ir alive while the finite element can retain its pointer-keyed DofToQuad cache.
  InitNonTensorBasis(fe, ir, num_comp, ceed, basis);
}

void InitSimplexBasisAtPoints(const mfem::FiniteElement &fe, bool grad_only,
                              CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  const auto geom = fe.GetGeomType();
  MFEM_VERIFY(geom == mfem::Geometry::TRIANGLE || geom == mfem::Geometry::TETRAHEDRON,
              "Simplex AtPoints basis requested for a non-simplex element!");
  // MAGMA's non-tensor AtPoints basis construction requires tabulation points that
  // overdetermine the complete polynomial space. One extra lattice degree provides a
  // deterministic overdetermined reconstruction. The registered rule has application
  // lifetime to satisfy MFEM's pointer-keyed DofToQuad cache.
  const int degree = std::max(0, fe.GetOrder() - (grad_only ? 1 : 0) + 1);
  InitBasisAtPoints(fe, GetRegisteredSimplexLatticeRule(geom, degree), num_comp, ceed,
                    basis);
}

void InitTetBasisAtPoints(const mfem::FiniteElement &fe, bool grad_only, CeedInt num_comp,
                          Ceed ceed, CeedBasis *basis)
{
  MFEM_VERIFY(fe.GetGeomType() == mfem::Geometry::TETRAHEDRON,
              "Tetrahedral AtPoints basis requested for a non-tetrahedral element!");
  InitSimplexBasisAtPoints(fe, grad_only, num_comp, ceed, basis);
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
