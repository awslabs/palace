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

std::size_t FiniteElementSpace::global_id = 0;

std::size_t FiniteElementSpace::GetGlobalId()
{
  std::size_t id;
  PalacePragmaOmp(critical(GetGlobalId))
  {
    id = global_id++;
  }
  return id;
}

std::size_t FiniteElementSpace::GetId() const
{
  PalacePragmaOmp(critical(GetId))
  {
    if (sequence != fespace->GetSequence())
    {
      id = GetGlobalId();
      sequence = fespace->GetSequence();
    }
  }
  return id;
}

void FiniteElementSpace::DestroyCeedObjects()
{
  for (auto [key, val] : basis)
  {
    Ceed ceed = key.first;
    PalaceCeedCall(ceed, CeedBasisDestroy(&val));
  }
  basis.clear();
  for (auto [key, val] : restr)
  {
    Ceed ceed = key.first;
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
  }
  restr.clear();
  for (auto [key, val] : interp_restr)
  {
    Ceed ceed = key.first;
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
  }
  interp_restr.clear();
  for (auto [key, val] : interp_range_restr)
  {
    Ceed ceed = key.first;
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&val));
  }
  interp_range_restr.clear();
}

const CeedBasis FiniteElementSpace::GetCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  const auto key = std::make_pair(ceed, geom);
  const auto it = basis.find(key);
  if (it != basis.end())
  {
    return it->second;
  }
  auto val = BuildCeedBasis(*this, ceed, geom);
  PalacePragmaOmp(critical(InitBasis))
  {
    basis.emplace(key, val);
  }
  return val;
}

const CeedElemRestriction
FiniteElementSpace::GetCeedElemRestriction(Ceed ceed, mfem::Geometry::Type geom,
                                           const std::vector<int> &indices) const
{
  const auto key = std::make_pair(ceed, geom);
  const auto it = restr.find(key);
  if (it != restr.end())
  {
    return it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices);
  PalacePragmaOmp(critical(InitRestriction))
  {
    restr.emplace(key, val);
  }
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
  const auto key = std::make_pair(ceed, geom);
  const auto it = interp_restr.find(key);
  if (it != interp_restr.end())
  {
    return it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices, true, false);
  PalacePragmaOmp(critical(InitInterpRestriction))
  {
    interp_restr.emplace(key, val);
  }
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
  const auto key = std::make_pair(ceed, geom);
  const auto it = interp_range_restr.find(key);
  if (it != interp_range_restr.end())
  {
    return it->second;
  }
  auto val = BuildCeedElemRestriction(*this, ceed, geom, indices, true, true);
  PalacePragmaOmp(critical(InitInterpRangeRestriction))
  {
    interp_range_restr.emplace(key, val);
  }
  return val;
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
