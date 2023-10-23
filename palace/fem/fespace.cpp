// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
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

const Operator &AuxiliaryFiniteElementSpace::BuildDiscreteInterpolator() const
{
  // G is always partially assembled.
  const int dim = GetParMesh()->Dimension();
  const auto aux_map_type = FEColl()->GetMapType(dim);
  const auto primal_map_type = primal_fespace.FEColl()->GetMapType(dim);
  if (aux_map_type == mfem::FiniteElement::VALUE &&
      primal_map_type == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator
    DiscreteLinearOperator interp(*this, primal_fespace);
    interp.AddDomainInterpolator<GradientInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), *this, primal_fespace, true);
  }
  else if (primal_map_type == mfem::FiniteElement::VALUE &&
           aux_map_type == mfem::FiniteElement::H_CURL)
  {
    // Discrete gradient interpolator (spaces reversed)
    DiscreteLinearOperator interp(primal_fespace, *this);
    interp.AddDomainInterpolator<GradientInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), primal_fespace, *this, true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_CURL &&
           primal_map_type == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator
    DiscreteLinearOperator interp(*this, primal_fespace);
    interp.AddDomainInterpolator<CurlInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), *this, primal_fespace, true);
  }
  else if (primal_map_type == mfem::FiniteElement::H_CURL &&
           aux_map_type == mfem::FiniteElement::H_DIV)
  {
    // Discrete curl interpolator (spaces reversed)
    DiscreteLinearOperator interp(primal_fespace, *this);
    interp.AddDomainInterpolator<CurlInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), primal_fespace, *this, true);
  }
  else if (aux_map_type == mfem::FiniteElement::H_DIV &&
           primal_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator
    DiscreteLinearOperator interp(*this, primal_fespace);
    interp.AddDomainInterpolator<DivergenceInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), *this, primal_fespace, true);
  }
  else if (primal_map_type == mfem::FiniteElement::H_DIV &&
           aux_map_type == mfem::FiniteElement::INTEGRAL)
  {
    // Discrete divergence interpolator (spaces reversed)
    DiscreteLinearOperator interp(primal_fespace, *this);
    interp.AddDomainInterpolator<DivergenceInterpolator>();
    G = std::make_unique<ParOperator>(interp.Assemble(), primal_fespace, *this, true);
  }
  else
  {
    MFEM_ABORT("Unsupported trial/test FE spaces for AuxiliaryFiniteElementSpace discrete "
               "interpolator!");
  }

  return *G;
}

const Operator &FiniteElementSpaceHierarchy::BuildProlongationAtLevel(std::size_t l) const
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
    P[l] =
        std::make_unique<ParOperator>(p.Assemble(), *fespaces[l], *fespaces[l + 1], true);
  }

  return *P[l];
}

}  // namespace palace
