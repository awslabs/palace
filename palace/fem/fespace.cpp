// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/restriction.hpp"
#include "linalg/rap.hpp"
#include "utils/omp.hpp"

namespace palace
{

const CeedBasis FiniteElementSpace::GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  // No two threads should ever be calling this simultaneously with the same Ceed context.
  auto it = basis.find(ceed);
  if (it == basis.end())
  {
    PalacePragmaOmp(critical(InitBasis))
    {
      it = basis.emplace(ceed, ceed::CeedGeomObjectMap<CeedBasis>()).first;
    }
  }
  auto &basis_map = it->second;
  auto basis_it = basis_map.find(geom);
  if (basis_it != basis_map.end())
  {
    return basis_it->second;
  }
  auto val = BuildCeedBasis(*this, ceed, geom);
  basis_map.emplace(geom, val);
  return val;
}

const CeedElemRestriction
FiniteElementSpace::GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                           const std::vector<int> &indices) const
{
  // No two threads should ever be calling this simultaneously with the same Ceed context.
  auto it = restr.find(ceed);
  if (it == restr.end())
  {
    PalacePragmaOmp(critical(InitRestriction))
    {
      it = restr.emplace(ceed, ceed::CeedGeomObjectMap<CeedElemRestriction>()).first;
    }
  }
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices);
  restr_map.emplace(geom, val);
  return val;
}

const CeedElemRestriction
FiniteElementSpace::GetInterpCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                 const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRestriction(fe))
  {
    return GetCeedElemRestriction(ceed, geom, indices);
  }
  // No two threads should ever be calling this simultaneously with the same Ceed context.
  auto it = interp_restr.find(ceed);
  if (it == interp_restr.end())
  {
    PalacePragmaOmp(critical(InitInterpRestriction))
    {
      it = interp_restr.emplace(ceed, ceed::CeedGeomObjectMap<CeedElemRestriction>()).first;
    }
  }
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices, true, false);
  restr_map.emplace(geom, val);
  return val;
}

const CeedElemRestriction
FiniteElementSpace::GetInterpRangeCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                                      const std::vector<int> &indices) const
{
  const mfem::FiniteElement &fe = *GetFEColl().FiniteElementForGeometry(geom);
  if (!HasUniqueInterpRangeRestriction(fe))
  {
    return GetInterpCeedElemRestriction(ceed, geom, indices);
  }
  // No two threads should ever be calling this simultaneously with the same Ceed context.
  auto it = interp_range_restr.find(ceed);
  if (it == interp_range_restr.end())
  {
    PalacePragmaOmp(critical(InitInterpRangeRestriction))
    {
      it = interp_range_restr.emplace(ceed, ceed::CeedGeomObjectMap<CeedElemRestriction>())
               .first;
    }
  }
  auto &restr_map = it->second;
  auto restr_it = restr_map.find(geom);
  if (restr_it != restr_map.end())
  {
    return restr_it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices, true, true);
  restr_map.emplace(geom, val);
  return val;
}

void FiniteElementSpace::DestroyCeedObjects()
{
  for (auto &[ceed, basis_map] : basis)
  {
    for (auto &[key, val] : basis_map)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&val));
    }
  }
  basis.clear();
  for (auto &[ceed, restr_map] : restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  restr.clear();
  for (auto &[ceed, restr_map] : interp_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  interp_restr.clear();
  for (auto &[ceed, restr_map] : interp_range_restr)
  {
    for (auto &[key, val] : restr_map)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
    }
  }
  interp_range_restr.clear();
}

CeedBasis FiniteElementSpace::BuildCeedBasis(const mfem::FiniteElementSpace &fespace,
                                             Ceed ceed, mfem::Geometry::Type geom)
{
  // Find the appropriate integration rule for the element.
  mfem::IsoparametricTransformation T;
  T.SetFE(fespace.GetMesh()->GetNodalFESpace()->FEColl()->FiniteElementForGeometry(geom));
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);

  // Build the libCEED basis.
  CeedBasis val;
  const mfem::FiniteElement &fe = *fespace.FEColl()->FiniteElementForGeometry(geom);
  const int vdim = fespace.GetVDim();
  ceed::InitBasis(fe, ir, vdim, ceed, &val);
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

const Operator &AuxiliaryFiniteElementSpace::BuildDiscreteInterpolator() const
{
  // Allow finite element spaces to be swapped in their order (intended as deriv(aux) ->
  // primal). G is always partially assembled.
  const int dim = Dimension();
  const bool swap =
      (GetFEColl().GetMapType(dim) == primal_fespace.GetFEColl().GetDerivMapType(dim));
  const FiniteElementSpace &trial_fespace = swap ? primal_fespace : *this;
  const FiniteElementSpace &test_fespace = swap ? *this : primal_fespace;
  const auto aux_map_type = trial_fespace.GetFEColl().GetMapType(dim);
  const auto primal_map_type = test_fespace.GetFEColl().GetMapType(dim);
  if (aux_map_type == mfem::FiniteElement::VALUE &&
      primal_map_type == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<GradientInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_CURL &&
           primal_map_type == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<CurlInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_DIV &&
           primal_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator
    DiscreteLinearOperator interp(trial_fespace, test_fespace);
    interp.AddDomainInterpolator<DivergenceInterpolator>();
    G = std::make_unique<ParOperator>(interp.PartialAssemble(), trial_fespace, test_fespace,
                                      true);
  }
  else
  {
    MFEM_ABORT("Unsupported trial/test FE spaces for AuxiliaryFiniteElementSpace discrete "
               "interpolator!");
  }

  return *G;
}

template <typename FESpace>
const Operator &
BaseFiniteElementSpaceHierarchy<FESpace>::BuildProlongationAtLevel(std::size_t l) const
{
  // P is always partially assembled.
  MFEM_VERIFY(l >= 0 && l < GetNumLevels() - 1,
              "Can only construct a finite element space prolongation with more than one "
              "space in the hierarchy!");
  if (&fespaces[l]->GetMesh() != &fespaces[l + 1]->GetMesh())
  {
    P[l] = std::make_unique<ParOperator>(
        std::make_unique<mfem::TransferOperator>(*fespaces[l], *fespaces[l + 1]),
        *fespaces[l], *fespaces[l + 1], true);
  }
  else
  {
    DiscreteLinearOperator p(*fespaces[l], *fespaces[l + 1]);
    p.AddDomainInterpolator<IdentityInterpolator>();
    P[l] = std::make_unique<ParOperator>(p.PartialAssemble(), *fespaces[l],
                                         *fespaces[l + 1], true);
  }

  return *P[l];
}

template class BaseFiniteElementSpaceHierarchy<FiniteElementSpace>;
template class BaseFiniteElementSpaceHierarchy<AuxiliaryFiniteElementSpace>;

}  // namespace palace
