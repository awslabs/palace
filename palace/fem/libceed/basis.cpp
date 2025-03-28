// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basis.hpp"

#include <mfem.hpp>
#include "utils/diagnostic.hpp"

namespace palace::ceed
{

namespace
{

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
  mfem::DenseMatrix qX(dim, Q);
  mfem::Vector qW(Q);
  for (int i = 0; i < Q; i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    qX(0, i) = ip.x;
    if (dim > 1)
    {
      qX(1, i) = ip.y;
    }
    if (dim > 2)
    {
      qX(2, i) = ip.z;
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
