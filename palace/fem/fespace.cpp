// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/libceed/utils.hpp"
#include "linalg/rap.hpp"
#include "utils/geodata.hpp"
#include "utils/omp.hpp"

namespace palace
{

std::size_t FiniteElementSpace::global_id = 0;

std::size_t FiniteElementSpace::GetId() const
{
  PalacePragmaOmp(critical(GetId))
  {
    if (!init || GetSequence() != prev_sequence)
    {
      id = global_id++;
      prev_sequence = GetSequence();
      init = true;
    }
  }
  return id;
}

FiniteElementSpace::~FiniteElementSpace()
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

CeedBasis FiniteElementSpace::BuildCeedBasis(Ceed ceed, mfem::Geometry::Type geom) const
{
  // Find the appropriate integration rule for the element.
  mfem::IsoparametricTransformation T;
  T->SetFE(GetParMesh()->GetNodes()->FEColl()->FiniteElementForGeometry(geom));
  const int q_order = fem::DefaultIntegrationOrder::Get(T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom, q_order);

  // Build the libCEED basis.
  const mfem::FiniteElement &fe = FEColl()->FiniteElementForGeometry(geom);
  const int vdim = GetVDim();
  const auto key = std::make_pair(ceed, geom);
  CeedBasis val;
  ceed::InitBasis(fe, ir, vdim, ceed, &val);
  PalacePragmaOmp(critical(InitBasis))
  {
    basis.emplace(key, val);
  }
  return val;
}

CeedElemRestriction FiniteElementSpace::BuildCeedElemRestriction(
    Ceed ceed, const std::vector<std::size_t> &indices, bool use_bdr, bool is_interp,
    bool is_interp_range) const
{
  const auto geom = use_bdr ? GetParMesh()->GetBdrElementGeometry(indices[0])
                            : GetParMesh()->GetElementGeometry(indices[0]);
  const auto key = std::make_pair(ceed, geom);
  CeedElementRestriction val;
  ceed::InitRestriction(*this, indices, use_bdr, is_interp, is_interp_range, ceed, &val);
  PalacePragmaOmp(critical(InitRestriction))
  {
    if (!is_interp && !is_interp_range)
    {
      restr.emplace(key, val);
    }
    else if (is_interp)
    {
      interp_restr.emplace(key, val);
    }
    else if (is_interp_range)
    {
      interp_range_restr.emplace(key, val);
    }
  }
  return val;
}

const Operator &AuxiliaryFiniteElementSpace::BuildDiscreteInterpolator() const
{
  // Allow finite element spaces to be swapped in their order (intended as deriv(aux) ->
  // primal). G is always partially assembled.
  const int dim = GetParMesh()->Dimension();
  const bool swap =
      (FEColl()->GetMapType(dim) == primal_fespace.FEColl()->GetDerivMapType(dim));
  const FiniteElementSpace &trial_fespace = swap ? primal_fespace : *this;
  const FiniteElementSpace &test_fespace = swap ? *this : primal_fespace;
  const auto aux_map_type = trial_fespace.FEColl()->GetMapType(dim);
  const auto primal_map_type = test_fespace.FEColl()->GetMapType(dim);
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
  if (fespaces[l]->GetParMesh() != fespaces[l + 1]->GetParMesh())
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
