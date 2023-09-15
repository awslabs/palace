// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basis.hpp"

namespace palace::ceed
{

namespace internal
{

std::unordered_map<BasisKey, CeedBasis, BasisHash> basis_map;

}  // namespace internal

namespace
{

inline CeedElemTopology GetCeedTopology(mfem::Geometry::Type geom)
{
  switch (geom)
  {
    case mfem::Geometry::SEGMENT:
      return CEED_TOPOLOGY_LINE;
    case mfem::Geometry::TRIANGLE:
      return CEED_TOPOLOGY_TRIANGLE;
    case mfem::Geometry::SQUARE:
      return CEED_TOPOLOGY_QUAD;
    case mfem::Geometry::TETRAHEDRON:
      return CEED_TOPOLOGY_TET;
    case mfem::Geometry::CUBE:
      return CEED_TOPOLOGY_HEX;
    case mfem::Geometry::PRISM:
      return CEED_TOPOLOGY_PRISM;
    case mfem::Geometry::PYRAMID:
      return CEED_TOPOLOGY_PYRAMID;
    default:
      MFEM_ABORT("This type of element is not supported!");
      return CEED_TOPOLOGY_LINE;  // Silence compiler warning
  }
}

void InitTensorBasis(const mfem::FiniteElementSpace &fespace, const mfem::FiniteElement &fe,
                     const mfem::IntegrationRule &ir, Ceed ceed, CeedBasis *basis)
{
  const mfem::DofToQuad &maps = fe.GetDofToQuad(ir, mfem::DofToQuad::TENSOR);
  const int dim = fe.GetDim();
  const int ncomp = fespace.GetVDim();
  const int P = maps.ndof;
  const int Q = maps.nqpt;
  mfem::Vector qX(Q), qW(Q);
  // The x-coordinates of the first `Q` points of the integration rule are the points of
  // the corresponding 1D rule. We also scale the weights accordingly.
  double w_sum = 0.0;
  for (int i = 0; i < Q; i++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(i);
    qX(i) = ip.x;
    qW(i) = ip.weight;
    w_sum += ip.weight;
  }
  qW *= 1.0 / w_sum;
  PalaceCeedCall(ceed, CeedBasisCreateTensorH1(ceed, dim, ncomp, P, Q, maps.Bt.GetData(),
                                               maps.Gt.GetData(), qX.GetData(),
                                               qW.GetData(), basis));
}

void InitNonTensorBasis(const mfem::FiniteElementSpace &fespace,
                        const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir,
                        Ceed ceed, CeedBasis *basis)
{
  const mfem::DofToQuad &maps = fe.GetDofToQuad(ir, mfem::DofToQuad::FULL);
  const int dim = fe.GetDim();
  const int ncomp = fespace.GetVDim();
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
    PalaceCeedCall(ceed, CeedBasisCreateHdiv(ceed, GetCeedTopology(fe.GetGeomType()), ncomp,
                                             P, Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                             qX.GetData(), qW.GetData(), basis));
  }
  else if (fe.GetMapType() == mfem::FiniteElement::H_CURL)
  {
    PalaceCeedCall(ceed,
                   CeedBasisCreateHcurl(ceed, GetCeedTopology(fe.GetGeomType()), ncomp, P,
                                        Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                        qX.GetData(), qW.GetData(), basis));
  }
  else
  {
    PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, GetCeedTopology(fe.GetGeomType()), ncomp,
                                           P, Q, maps.Bt.GetData(), maps.Gt.GetData(),
                                           qX.GetData(), qW.GetData(), basis));
  }
}

#if 0
void InitCeedInterpolatorBasis(const mfem::FiniteElementSpace &trial_fespace,
                                      const mfem::FiniteElementSpace &test_fespace,
                                      const mfem::FiniteElement &trial_fe,
                                      const mfem::FiniteElement &test_fe,
                                      Ceed ceed,
                                      CeedBasis *basis)
{
   // Basis projection operator using libCEED
   CeedBasis trial_basis, test_basis;
   const int P = std::max(trial_fe.GetDof(), test_fe.GetDof()), ir_order_max = 100;
   int ir_order = std::max(trial_fe.GetOrder(), test_fe.GetOrder());
   for (; ir_order < ir_order_max; ir_order++)
   {
      if (IntRules.Get(trial_fe.GetGeomType(), ir_order).GetNPoints() >= P) { break; }
   }
   const mfem::IntegrationRule &ir = IntRules.Get(trial_fe.GetGeomType(), ir_order);
   InitBasis(trial_fespace, trial_fe, ir, ceed, &trial_basis);
   InitBasis(test_fespace, test_fe, ir, ceed, &test_basis);
   PalaceCeedCall(ceed, CeedBasisCreateProjection(trial_basis, test_basis, basis));
}
#endif

void InitMFEMInterpolatorBasis(const mfem::FiniteElementSpace &trial_fespace,
                               const mfem::FiniteElementSpace &test_fespace,
                               const mfem::FiniteElement &trial_fe,
                               const mfem::FiniteElement &test_fe, Ceed ceed,
                               CeedBasis *basis)
{
  MFEM_VERIFY(
      trial_fespace.GetVDim() == test_fespace.GetVDim(),
      "libCEED discrete linear operator requires same vdim for trial and test FE spaces!");
  const int dim = trial_fe.GetDim();
  const int ncomp = trial_fespace.GetVDim();
  const int trial_P = trial_fe.GetDof();
  const int test_P = test_fe.GetDof();
  mfem::DenseMatrix qX(dim, test_P), Gt(trial_P, test_P * dim), Bt;
  mfem::Vector qW(test_P);
  mfem::IsoparametricTransformation dummy;
  dummy.SetIdentityTransformation(trial_fe.GetGeomType());
  if (trial_fe.GetMapType() == test_fe.GetMapType())
  {
    // Prolongation
    test_fe.GetTransferMatrix(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::VALUE &&
           test_fe.GetMapType() == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator
    test_fe.ProjectGrad(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::H_CURL &&
           test_fe.GetMapType() == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator
    test_fe.ProjectCurl(trial_fe, dummy, Bt);
  }
  else if (trial_fe.GetMapType() == mfem::FiniteElement::H_DIV &&
           test_fe.GetMapType() == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator
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
  PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, GetCeedTopology(trial_fe.GetGeomType()),
                                         ncomp, trial_P, test_P, Bt.GetData(), Gt.GetData(),
                                         qX.GetData(), qW.GetData(), basis));
}

}  // namespace

void InitBasis(const mfem::FiniteElementSpace &fespace, const mfem::FiniteElement &fe,
               const mfem::IntegrationRule &ir, Ceed ceed, CeedBasis *basis)
{
  // Check for fespace -> basis in hash table.
  const int ncomp = fespace.GetVDim();
  internal::BasisKey key(ceed, (void *)&fe, (void *)&ir, ncomp);

  // Initialize or retrieve key values.
  auto basis_itr = internal::basis_map.find(key);
  if (basis_itr == internal::basis_map.end())
  {
    const bool tensor = dynamic_cast<const mfem::TensorBasisElement *>(&fe) != nullptr;
    const bool vector = fe.GetRangeType() == mfem::FiniteElement::VECTOR;
    if (tensor && !vector)
    {
      InitTensorBasis(fespace, fe, ir, ceed, basis);
    }
    else
    {
      InitNonTensorBasis(fespace, fe, ir, ceed, basis);
    }
    PalacePragmaOmp(critical(InitBasis))
    {
      internal::basis_map[key] = *basis;
    }
  }
  else
  {
    *basis = basis_itr->second;
  }
}

void InitInterpolatorBasis(const mfem::FiniteElementSpace &trial_fespace,
                           const mfem::FiniteElementSpace &test_fespace,
                           const mfem::FiniteElement &trial_fe,
                           const mfem::FiniteElement &test_fe, Ceed ceed, CeedBasis *basis)
{
  // Check for fespace -> basis in hash table.
  const int ncomp = trial_fespace.GetVDim();  // Assumed same as test_fespace
  internal::BasisKey key(ceed, (void *)&trial_fe, (void *)&test_fe, ncomp);

  // Initialize or retrieve key values.
  auto basis_itr = internal::basis_map.find(key);
  if (basis_itr == internal::basis_map.end())
  {
#if 0
       if (trial_fe.GetMapType() == test_fe.GetMapType())
       {
          InitCeedInterpolatorBasis(trial_fespace, test_fespace, trial_fe, test_fe, ceed, basis);
       }
       else
#endif
    {
      InitMFEMInterpolatorBasis(trial_fespace, test_fespace, trial_fe, test_fe, ceed,
                                basis);
    }
    PalacePragmaOmp(critical(InitBasis))
    {
      internal::basis_map[key] = *basis;
    }
  }
  else
  {
    *basis = basis_itr->second;
  }
}

}  // namespace palace::ceed
