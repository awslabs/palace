// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/restriction.hpp"
#include "linalg/rap.hpp"

namespace palace
{

CeedBasis FiniteElementSpace::GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  auto it = basis.find(ceed);
  MFEM_ASSERT(it != basis.end(), "Unknown Ceed context in GetCeedBasis!");
  auto &basis_map = it->second;
  auto basis_it = basis_map.find(geom);
  if (basis_it != basis_map.end())
  {
    return basis_it->second;
  }
  return basis_map.emplace(geom, BuildCeedBasis(*this, ceed, geom)).first->second;
}

CeedElemRestriction
FiniteElementSpace::GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                           const std::vector<int> &indices) const
{
  auto it = restr.find(ceed);
  MFEM_ASSERT(it != restr.end(), "Unknown Ceed context in GetCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  return restr_map.emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                 const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRestriction(fe))
  {
    return GetCeedElemRestriction(ceed, geom, indices);
  }
  auto it = interp_restr.find(ceed);
  MFEM_ASSERT(it != interp_restr.end(),
              "Unknown Ceed context in GetInterpCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  return restr_map
      .emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices, true, false))
      .first->second;
}

CeedElemRestriction
FiniteElementSpace::GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                      const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRangeRestriction(fe))
  {
    return GetInterpCeedElemRestriction(ceed, geom, indices);
  }
  auto it = interp_range_restr.find(ceed);
  MFEM_ASSERT(it != interp_range_restr.end(),
              "Unknown Ceed context in GetInterpRangeCeedElemRestriction!");
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  return restr_map
      .emplace(geom, BuildCeedElemRestriction(*this, ceed, geom, indices, true, true))
      .first->second;
}

void FiniteElementSpace::ResetCeedObjects()
{
  for (auto &[ceed, basis_map] : basis)
  {
    for (auto &[key, val] : basis_map)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&val));
    }
  }
  for (auto &[ceed, restr_map] : restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  for (auto &[ceed, restr_map] : interp_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  for (auto &[ceed, restr_map] : interp_range_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  basis.clear();
  restr.clear();
  interp_restr.clear();
  interp_range_restr.clear();
  for (std::size_t i = 0; i < ceed::internal::GetCeedObjects().size(); i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    basis.emplace(ceed, ceed::GeometryObjectMap<CeedBasis>());
    restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
    interp_range_restr.emplace(ceed, ceed::GeometryObjectMap<CeedElemRestriction>());
  }
}

CeedBasis FiniteElementSpace::BuildCeedBasis(const mfem::FiniteElementSpace &fespace,
                                             Ceed ceed, mfem::Geometry::Type geom)
{
  // Find the appropriate integration rule for the element.
  mfem::IsoparametricTransformation T;
  const mfem::FiniteElement *fe_nodal =
      fespace.GetMesh()->GetNodalFESpace()->FEColl()->FiniteElementForGeometry(geom);
  if (!fe_nodal)
  {
    fe_nodal =
        fespace.GetMesh()->GetNodalFESpace()->FEColl()->TraceFiniteElementForGeometry(geom);
  }
  T.SetFE(fe_nodal);
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);

  // Build the libCEED basis.
  CeedBasis val;
  const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
  if (!fe)
  {
    fe = fespace.FEColl()->TraceFiniteElementForGeometry(geom);
  }
  const int vdim = fespace.GetVDim();
  ceed::InitBasis(*fe, ir, vdim, ceed, &val);
  return val;
}

CeedElemRestriction FiniteElementSpace::BuildCeedElemRestriction(
    const mfem::FiniteElementSpace &fespace, Ceed ceed, mfem::Geometry::Type geom,
    const std::vector<int> &indices, bool is_interp, bool is_interp_range)
{
  // Construct the libCEED element restriction for this element type.
  CeedElemRestriction val;
  const bool use_bdr = (mfem::Geometry::Dimension[geom] != fespace.GetMesh()->Dimension());
  ceed::InitRestriction(fespace, indices, use_bdr, is_interp, is_interp_range, ceed, &val);
  return val;
}

const Operator &FiniteElementSpace::BuildDiscreteInterpolator() const
{
  // Allow finite element spaces to be swapped in their order (intended as deriv(aux) ->
  // primal). G is always partially assembled.
  const int dim = Dimension();
  const bool swap =
      (aux_fespace->GetFEColl().GetMapType(dim) == GetFEColl().GetDerivMapType(dim));
  MFEM_VERIFY(!swap, "Incorrect order for primal/auxiliary (test/trial) spaces in discrete "
                     "interpolator construction!");
  MFEM_VERIFY(
      GetFEColl().GetMapType(dim) == aux_fespace->GetFEColl().GetDerivMapType(dim),
      "Unsupported trial/test FE spaces for FiniteElementSpace discrete interpolator!");
  const FiniteElementSpace &trial_fespace = !swap ? *aux_fespace : *this;
  const FiniteElementSpace &test_fespace = !swap ? *this : *aux_fespace;
  const auto aux_map_type = trial_fespace.GetFEColl().GetMapType(dim);
  const auto primal_map_type = test_fespace.GetFEColl().GetMapType(dim);
  if (aux_map_type == mfem::FiniteElement::VALUE &&
      primal_map_type == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator.
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<GradientInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_CURL &&
           primal_map_type == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator.
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<CurlInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_DIV &&
           primal_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator.
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<DivergenceInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else
  {
    MFEM_ABORT(
        "Unsupported trial/test FE spaces for FiniteElementSpace discrete interpolator!");
  }

  return *G;
}

const Operator &FiniteElementSpaceHierarchy::BuildProlongationAtLevel(std::size_t l) const
{
  // P is always partially assembled.
  MFEM_VERIFY(l + 1 < GetNumLevels(),
              "Can only construct a finite element space prolongation with more than one "
              "space in the hierarchy!");
  Mpi::Print("BuildProlongation at level: {}\n", l);
  if (&fespaces[l]->GetMesh() != &fespaces[l + 1]->GetMesh())
  {
    /*
    Mpi::Print("fespace.cpp L233 different meshes get refine transfer operator?\n");
    refine_op[l] = std::make_unique<mfem::TransferOperator>(*fespaces[l], *fespaces[l + 1]);
    std::cout << "refine_op rank: " << Mpi::Rank(fespaces[l + 1]->GetComm()) << " width/height: " << refine_op[l]->Width() << " " << refine_op[l]->Height() << "\n";
    Mpi::Print("fespace.cpp L235 fespaces[l + 1]->GetParMesh().Rebalance()!\n");
    fespaces[l + 1]->GetParMesh().Rebalance();
    fespaces[l + 1]->GetMesh().Update(); //not sure if needed
    fespaces[l + 1]->Update(); // not sure if needed...
    Mpi::Print("fespace.cpp L239 GetUpdateOperator\n");
    rebalance_op[l].reset(const_cast<Operator*>(fespaces[l + 1]->Get().GetUpdateOperator()));
    std::cout << "rebalance_op rank: " << Mpi::Rank(fespaces[l + 1]->GetComm()) << " width/height: " << rebalance_op[l]->Width() << " " << rebalance_op[l]->Height() << "\n";
    Mpi::Print("fespace.cpp L241 combine refine and rebalance operators\n");
    transfer_op[l] = std::make_unique<ProductOperator>(*rebalance_op[l], *refine_op[l]);
    std::cout << "transfer_op rank: " << Mpi::Rank(fespaces[l + 1]->GetComm()) << " width/height: " << transfer_op[l]->Width() << " " << transfer_op[l]->Height() << "\n";
    // reconstruct fespace?
    //FiniteElementSpace h1d_fespace(h1_fespace.GetMesh(), &h1_fespace.GetFEColl(), space_dim, mfem::Ordering::byVDIM);
    auto mesh = fespaces[l + 1]->GetMesh();
    auto fecoll = &fespaces[l + 1]->GetFEColl();
    //FiniteElementSpace test_space(mesh, fecoll);
    //fespaces[l + 1] = std::make_unique<FiniteElementSpace>(mesh, fecoll);//
    fespaces[l + 1].reset(new FiniteElementSpace(mesh, fecoll));//
    Mpi::Print("fespace.cpp L243 create P[l]\n");
    */
    P[l] = std::make_unique<ParOperator>(
        std::make_unique<mfem::TransferOperator>(*fespaces[l], *fespaces[l + 1]), // transferoperator is problematic if a rebalance happened
        //std::make_unique<ProductOperator>(*rebalance_op[l], *refine_op[l]),
        *fespaces[l], *fespaces[l + 1], true);
    Mpi::Print("fespace.cpp L248\n");
  }
  else
  {
    Mpi::Print("fespace.cpp L255 same mesh, calling update?\n");
    fespaces[l + 1]->GetMesh().Update(); // not sure if needed?
    fespaces[l + 1]->Update(); // not sure if needed...
    fespaces[l + 1]->Get().Update(); // not sure if needed...
    fespaces[l]->GetMesh().Update(); // not sure if needed?
    fespaces[l]->Update(); // not sure if needed...
    fespaces[l]->Get().Update(); // not sure if needed...
    Mpi::Print("fespace.cpp L258 discrete linear operator\n");
    DiscreteLinearOperator p(*fespaces[l], *fespaces[l + 1]);
    Mpi::Print("fespace.cpp L260\n");
    p.AddDomainInterpolator<IdentityInterpolator>();
    Mpi::Print("fespace.cpp L262\n");
    P[l] = std::make_unique<ParOperator>(p.PartialAssemble(), *fespaces[l],
                                         *fespaces[l + 1], true);
    Mpi::Print("fespace.cpp L265\n");
  }

  return *P[l];
}

}  // namespace palace
