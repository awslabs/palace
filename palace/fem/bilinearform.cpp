// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bilinearform.hpp"

#include "fem/fespace.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "utils/omp.hpp"

namespace palace
{

void BilinearForm::AssembleQuadratureData()
{
  for (auto &integ : domain_integs)
  {
    integ->AssembleQuadratureData();
  }
  for (auto &integ : boundary_integs)
  {
    integ->AssembleQuadratureData();
  }
}

std::unique_ptr<ceed::Operator>
BilinearForm::PartialAssemble(const FiniteElementSpace &trial_fespace,
                              const FiniteElementSpace &test_fespace) const
{
  MFEM_VERIFY(&trial_fespace.GetMesh() == &test_fespace.GetMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  const auto &mesh = trial_fespace.GetMesh();

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
  PalacePragmaOmp(parallel if (nt > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      const auto trial_map_type =
          trial_fespace.GetFEColl().GetMapType(mfem::Geometry::Dimension[geom]);
      const auto test_map_type =
          test_fespace.GetFEColl().GetMapType(mfem::Geometry::Dimension[geom]);

      if (mfem::Geometry::Dimension[geom] == mesh.Dimension() && !domain_integs.empty())
      {
        // Assemble domain integrators on this element geometry type.
        CeedElemRestriction trial_restr =
            trial_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
        CeedElemRestriction test_restr =
            test_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
        CeedBasis trial_basis = trial_fespace.GetCeedBasis(ceed, geom);
        CeedBasis test_basis = test_fespace.GetCeedBasis(ceed, geom);

        for (const auto &integ : domain_integs)
        {
          CeedOperator sub_op;
          integ->SetMapTypes(trial_map_type, test_map_type);
          integ->Assemble(ceed, trial_restr, test_restr, trial_basis, test_basis,
                          data.geom_data, data.geom_data_restr, &sub_op);
          op->AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator
        }
      }
      else if (mfem::Geometry::Dimension[geom] == mesh.Dimension() - 1 &&
               !boundary_integs.empty())
      {
        // Assemble boundary integrators on this element geometry type.
        CeedElemRestriction trial_restr =
            trial_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
        CeedElemRestriction test_restr =
            test_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
        CeedBasis trial_basis = trial_fespace.GetCeedBasis(ceed, geom);
        CeedBasis test_basis = test_fespace.GetCeedBasis(ceed, geom);

        for (const auto &integ : boundary_integs)
        {
          CeedOperator sub_op;
          integ->SetMapTypes(trial_map_type, test_map_type);
          integ->Assemble(ceed, trial_restr, test_restr, trial_basis, test_basis,
                          data.geom_data, data.geom_data_restr, &sub_op);
          op->AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator
        }
      }
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  op->Finalize();

  return op;
}

std::unique_ptr<hypre::HypreCSRMatrix> BilinearForm::FullAssemble(const ceed::Operator &op,
                                                                  bool skip_zeros, bool set)
{
  return ceed::CeedOperatorFullAssemble(op, skip_zeros, set);
}

namespace
{

bool UseFullAssembly(const FiniteElementSpace &trial_fespace,
                     const FiniteElementSpace &test_fespace, int pa_order_threshold)
{
  // Returns order such that the miniumum for all element types is 1. MFEM's
  // RT_FECollection actually already returns order + 1 for GetOrder() for historical
  // reasons.
  const auto &trial_fec = trial_fespace.GetFEColl();
  const auto &test_fec = test_fespace.GetFEColl();
  int max_order = std::max(
      dynamic_cast<const mfem::L2_FECollection *>(&trial_fec) ? trial_fec.GetOrder() + 1
                                                              : trial_fec.GetOrder(),
      dynamic_cast<const mfem::L2_FECollection *>(&test_fec) ? test_fec.GetOrder() + 1
                                                             : test_fec.GetOrder());
  return (max_order < pa_order_threshold);
}

bool UseFullAssembly(const FiniteElementSpace &fespace, int pa_order_threshold)
{
  return UseFullAssembly(fespace, fespace, pa_order_threshold);
}

}  // namespace

std::unique_ptr<Operator> BilinearForm::Assemble(bool skip_zeros) const
{
  if (UseFullAssembly(trial_fespace, test_fespace, pa_order_threshold))
  {
    return FullAssemble(skip_zeros);
  }
  else
  {
    return PartialAssemble();
  }
}

std::vector<std::unique_ptr<Operator>>
BilinearForm::Assemble(const FiniteElementSpaceHierarchy &fespaces, bool skip_zeros,
                       std::size_t l0) const
{
  // Only available for square operators (same test and trial spaces).
  MFEM_VERIFY(&trial_fespace == &test_fespace &&
                  &fespaces.GetFinestFESpace() == &trial_fespace,
              "Assembly on a FiniteElementSpaceHierarchy should have the same BilinearForm "
              "spaces and fine space of the hierarchy!");

  // First partially assemble all of the operators.
  MFEM_VERIFY(l0 < fespaces.GetNumLevels(),
              "No levels available for operator coarsening (l0 = " << l0 << ")!");
  std::vector<std::unique_ptr<ceed::Operator>> pa_ops;
  pa_ops.reserve(fespaces.GetNumLevels() - l0);
  for (std::size_t l = l0; l < fespaces.GetNumLevels(); l++)
  {
    if (l > l0 && &fespaces.GetFESpaceAtLevel(l).GetMesh() ==
                      &fespaces.GetFESpaceAtLevel(l - 1).GetMesh())
    {
      pa_ops.push_back(
          ceed::CeedOperatorCoarsen(*pa_ops.back(), fespaces.GetFESpaceAtLevel(l)));
    }
    else
    {
      pa_ops.push_back(
          PartialAssemble(fespaces.GetFESpaceAtLevel(l), fespaces.GetFESpaceAtLevel(l)));
    }
  }

  // Construct the final operators using full or partial assemble as needed. We do not
  // force the coarse-level operator to be fully assembled always, it will be only assembled
  // as needed for parallel assembly.
  std::vector<std::unique_ptr<Operator>> ops;
  ops.reserve(fespaces.GetNumLevels() - l0);
  for (std::size_t l = l0; l < fespaces.GetNumLevels(); l++)
  {
    if (UseFullAssembly(fespaces.GetFESpaceAtLevel(l), pa_order_threshold))
    {
      ops.push_back(FullAssemble(*pa_ops[l - l0], skip_zeros));
    }
    else
    {
      ops.push_back(std::move(pa_ops[l - l0]));
    }
  }

  return ops;
}

std::unique_ptr<ceed::Operator> DiscreteLinearOperator::PartialAssemble() const
{
  MFEM_VERIFY(&trial_fespace.GetMesh() == &test_fespace.GetMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  const auto &mesh = trial_fespace.GetMesh();

  // Initialize the operator.
  auto op =
      std::make_unique<ceed::Operator>(test_fespace.GetVSize(), trial_fespace.GetVSize());

  // Assemble the libCEED operator in parallel, each thread builds a composite operator.
  // This should work fine if some threads create an empty operator (no elements or boundary
  // elements).
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel if (nt > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      if (mfem::Geometry::Dimension[geom] == mesh.Dimension() && !domain_interps.empty())
      {
        // Assemble domain interpolators on this element geometry type.
        CeedElemRestriction trial_restr =
            trial_fespace.GetInterpCeedElemRestriction(ceed, geom, data.indices);
        CeedElemRestriction test_restr =
            test_fespace.GetInterpRangeCeedElemRestriction(ceed, geom, data.indices);

        // Construct the interpolator basis.
        CeedBasis interp_basis;
        const mfem::FiniteElement &trial_fe =
            *trial_fespace.GetFEColl().FiniteElementForGeometry(geom);
        const mfem::FiniteElement &test_fe =
            *test_fespace.GetFEColl().FiniteElementForGeometry(geom);
        const int trial_vdim = trial_fespace.GetVDim();
        const int test_vdim = test_fespace.GetVDim();
        ceed::InitInterpolatorBasis(trial_fe, test_fe, trial_vdim, test_vdim, ceed,
                                    &interp_basis);

        for (const auto &interp : domain_interps)
        {
          CeedOperator sub_op, sub_op_t;
          interp->Assemble(ceed, trial_restr, test_restr, interp_basis, &sub_op, &sub_op_t);
          op->AddSubOperator(sub_op, sub_op_t);  // Sub-operator owned by ceed::Operator
        }

        // Basis is owned by the operator.
        PalaceCeedCall(ceed, CeedBasisDestroy(&interp_basis));
      }
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  op->Finalize();

  // Construct dof multiplicity vector for scaling to account for dofs shared between
  // elements (on host, then copy to device).
  Vector test_multiplicity(test_fespace.GetVSize());
  test_multiplicity = 0.0;
  auto *h_mult = test_multiplicity.HostReadWrite();
  PalacePragmaOmp(parallel)
  {
    mfem::Array<int> dofs;
    mfem::DofTransformation dof_trans;
    PalacePragmaOmp(for schedule(static))
    for (int i = 0; i < test_fespace.GetMesh().GetNE(); i++)
    {
      test_fespace.Get().GetElementVDofs(i, dofs, dof_trans);
      for (int j = 0; j < dofs.Size(); j++)
      {
        const int k = dofs[j];
        PalacePragmaOmp(atomic update)
        h_mult[(k >= 0) ? k : -1 - k] += 1.0;
      }
    }
  }
  test_multiplicity.UseDevice(true);
  test_multiplicity.Reciprocal();
  op->SetDofMultiplicity(std::move(test_multiplicity));

  return op;
}

}  // namespace palace
