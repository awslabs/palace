// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "utils/omp.hpp"

namespace palace
{

std::size_t FiniteElementSpace::GetGlobalId()
{
  static std::size_t global_id = 0;
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
    if (sequence != fespace.GetSequence())
    {
      id = GetGlobalId();
      sequence = fespace.GetSequence();
    }
  }
  return id;
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
