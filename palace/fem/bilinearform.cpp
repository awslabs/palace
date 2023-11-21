// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bilinearform.hpp"

#include "fem/libceed/basis.hpp"
#include "fem/libceed/utils.hpp"
#include "utils/omp.hpp"

namespace palace
{

std::unique_ptr<ceed::Operator> BilinearForm::PartialAssemble() const
{
  MFEM_VERIFY(trial_fespace.GetParMesh() == test_fespace.GetParMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  const mfem::ParMesh &mesh = *trial_fespace.GetParMesh();

  // Initialize the operator.
  std::unique_ptr<ceed::Operator> op;
  if (&trial_fespace == &test_fespace)
  {
    op = std::make_unique<ceed::SymmetricOperator>(test_fespace.GetVSize(),
                                                   trial_fespace.GetVSize());
  }
  else
  {
    op =
        std::make_unique<ceed::Operator>(test_fespace.GetVSize(), trial_fespace.GetVSize());
  }

  // Assemble the libCEED operator in parallel, each thread builds a composite operator.
  // This should work fine if some threads create an empty operator (no elements or boundary
  // elements).
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];

    // Initialize the composite operator on each thread.
    CeedOperator loc_op;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));

    for (const auto &[key, val] : mat_op.GetElementIndices())
    {
      if (key.first != ceed)
      {
        continue;
      }
      const auto geom = key.second;
      const auto &indices = val;
      const auto &geom_data = mat_op.GetGeomFactorData(ceed, geom);
      CeedElementRestriction trial_restr =
          trial_fespace.GetCeedElemRestriction(ceed, geom, indices);
      CeedElementRestriction test_restr =
          test_fespace.GetCeedElemRestriction(ceed, geom, indices);
      CeedBasis trial_basis = trial_fespace.GetCeedBasis(ceed, geom);
      CeedBasis test_basis = test_fespace.GetCeedBasis(ceed, geom);

      if (mfem::Geometry::Dimension[geom] == mesh.Dimension())
      {
        // Assemble domain integrators on this element geometry type.
        for (const auto &integ : domain_integs)
        {
          CeedOperator sub_op;
          integ->Assemble(geom_data, ceed, trial_restr, test_restr, trial_basis, test_basis,
                          &sub_op);
          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
        }
      }
      else
      {
        // Assemble boundary integrators on this element geometry type.
        for (const auto &integ : boundary_integs)
        {
          CeedOperator sub_op;
          integ->AssembleBoundary(geom_data, ceed, trial_restr, test_restr, trial_basis,
                                  test_basis, &sub_op);
          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
        }
      }
    }
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    op->AddOper(loc_op);  // Thread-safe
  }

  return op;
}

std::unique_ptr<mfem::SparseMatrix> BilinearForm::FullAssemble(const ceed::Operator &op,
                                                               bool skip_zeros, bool set)
{
  return ceed::CeedOperatorFullAssemble(op, skip_zeros, set);
}

std::unique_ptr<ceed::Operator> DiscreteLinearOperator::PartialAssemble() const
{
  MFEM_VERIFY(trial_fespace.GetParMesh() == test_fespace.GetParMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  const mfem::ParMesh &mesh = *trial_fespace.GetParMesh();

  // Initialize the operator.
  auto op =
      std::make_unique<ceed::Operator>(test_fespace.GetVSize(), trial_fespace.GetVSize());

  // Assemble the libCEED operator in parallel, each thread builds a composite operator.
  // This should work fine if some threads create an empty operator (no elements or bounday
  // elements).
  MFEM_VERIFY(domain_integs.size() == 1,
              "DiscreteLinearOperator should only have a single domain interpolator!");
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];

    // Initialize the composite operators for each thread.
    CeedOperator loc_op, loc_op_t;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op_t));

    for (const auto &[key, val] : mat_op.GetIndices(ceed))
    {
      const auto geom = key;
      const std::vector<int> &indices = val;

      if (mfem::Geometry::Dimension[geom] == mesh.Dimension())
      {
        // Assemble domain interpolators on this element geometry type.
        CeedElementRestriction trial_restr =
            trial_fespace.GetInterpCeedElemRestriction(ceed, geom, indices);
        CeedElementRestriction test_restr =
            test_fespace.GetInterpRangeCeedElemRestriction(ceed, geom, indices);

        // Construct the interpolator basis.
        CeedBasis interp_basis;
        const mfem::FiniteElement &trial_fe =
            trial_fespace.FEColl()->FiniteElementForGeometry(geom);
        const mfem::FiniteElement &test_fe =
            test_fespace.FEColl()->FiniteElementForGeometry(geom);
        const int trial_vdim = trial_fespace.GetVDim();
        const int test_vdim = test_fespace.GetVDim();
        ceed::InitInterpolatorBasis(trial_fe, test_fe, trial_vdim, test_vdim, ceed,
                                    &interp_basis);

        for (const auto &interp : domain_interps)
        {
          CeedOperator sub_op, sub_op_t;
          interp->Assemble(ceed, trial_restr, test_restr, interp_basis, &sub_op, &sub_op_t);
          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op_t, sub_op_t));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
        }

        // Basis is owned by the operator.
        PalaceCeedCall(ceed, CeedBasisDestroy(&interp_basis));
      }
    }
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op_t));
    op->AddOper(loc_op, loc_op_t);  // Thread-safe
  }

  // Construct dof multiplicity vector for scaling to account for dofs shared between
  // elements (on host, then copy to device).
  Vector test_multiplicity(test_fespace.GetVSize());
  test_multiplicity = 0.0;
  mfem::Array<int> dofs;
  auto *h_mult = test_multiplicity.HostReadWrite();
  for (int i = 0; i < test_fespace.GetNE(); i++)
  {
    test_fespace.GetElementVDofs(i, dofs);
    for (int j = 0; j < dofs.Size(); j++)
    {
      const int k = dofs[j];
      h_mult[(k >= 0) ? k : -1 - k] += 1.0;
    }
  }
  test_multiplicity.UseDevice(true);
  test_multiplicity.Reciprocal();
  op->SetDofMultiplicity(std::move(test_multiplicity));

  return op;
}

}  // namespace palace
