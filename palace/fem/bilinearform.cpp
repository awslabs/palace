// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bilinearform.hpp"

#include <unordered_map>
#include <ceed.h>
#include "fem/fespace.hpp"
#include "fem/libceed/hash.hpp"
#include "fem/libceed/operator.hpp"
#include "fem/libceed/utils.hpp"
#include "utils/omp.hpp"

namespace palace
{

namespace
{

using ceed::internal::FiniteElementKey;
using ceed::internal::FiniteElementPairHash;
using ceed::internal::FiniteElementPairKey;

// Count the number of elements of each type in the local mesh.
std::unordered_map<FiniteElementPairKey, std::vector<int>, FiniteElementPairHash>
GetElementIndices(const mfem::ParFiniteElementSpace &trial_fespace,
                  const mfem::ParFiniteElementSpace &test_fespace, bool use_bdr, int start,
                  int stop)
{
  std::unordered_map<FiniteElementPairKey, int, FiniteElementPairHash> counts, offsets;
  std::unordered_map<FiniteElementPairKey, std::vector<int>, FiniteElementPairHash>
      element_indices;

  // Count the number of elements of each type and order.
  for (int i = start; i < stop; i++)
  {
    const mfem::FiniteElement &trial_fe =
        use_bdr ? *trial_fespace.GetBE(i) : *trial_fespace.GetFE(i);
    const mfem::FiniteElement &test_fe =
        use_bdr ? *test_fespace.GetBE(i) : *test_fespace.GetFE(i);
    FiniteElementPairKey key =
        std::make_pair(FiniteElementKey(trial_fe), FiniteElementKey(test_fe));
    auto value = counts.find(key);
    if (value == counts.end())
    {
      counts[key] = 1;
    }
    else
    {
      value->second++;
    }
  }

  // Populate the indices arrays for each element type.
  for (const auto &value : counts)
  {
    offsets[value.first] = 0;
    element_indices[value.first] = std::vector<int>(value.second);
  }
  for (int i = start; i < stop; i++)
  {
    const mfem::FiniteElement &trial_fe =
        use_bdr ? *trial_fespace.GetBE(i) : *trial_fespace.GetFE(i);
    const mfem::FiniteElement &test_fe =
        use_bdr ? *test_fespace.GetBE(i) : *test_fespace.GetFE(i);
    FiniteElementPairKey key =
        std::make_pair(FiniteElementKey(trial_fe), FiniteElementKey(test_fe));
    int &offset = offsets[key];
    std::vector<int> &indices = element_indices[key];
    indices[offset++] = i;
  }

  return element_indices;
}

}  // namespace

std::unique_ptr<ceed::Operator> BilinearForm::Assemble() const
{
  MFEM_VERIFY(trial_fespace.GetParMesh() == test_fespace.GetParMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  mfem::ParMesh &mesh = *trial_fespace.GetParMesh();
  {
    // In the following, we copy the mesh FE space for the nodes as a
    // palace::FiniteElementSpace and replace it in the nodal grid function. Unfortunately
    // mfem::ParFiniteElementSpace does not have a move constructor to make this more
    // efficient, but it's only done once for the lifetime of the mesh.
    mesh.EnsureNodes();
    mfem::GridFunction *mesh_nodes = mesh.GetNodes();
    mfem::FiniteElementSpace *mesh_fespace = mesh_nodes->FESpace();
    MFEM_VERIFY(dynamic_cast<mfem::ParFiniteElementSpace *>(mesh_fespace),
                "Unexpected non-parallel FiniteElementSpace for mesh nodes!");
    if (!dynamic_cast<FiniteElementSpace *>(mesh_fespace))
    {
      // Ensure the FiniteElementCollection associated with the original nodes is not
      // deleted.
      auto *new_mesh_fespace =
          new FiniteElementSpace(*static_cast<mfem::ParFiniteElementSpace *>(mesh_fespace));
      mfem::FiniteElementCollection *mesh_fec = mesh_nodes->OwnFEC();
      MFEM_VERIFY(mesh_fec, "Replacing the FiniteElementSpace for mesh nodes is only "
                            "possible when it owns its fec/fes members!");
      mesh_nodes->MakeOwner(nullptr);
      mesh.SetNodalFESpace(new_mesh_fespace);
      mfem::GridFunction *new_mesh_nodes = mesh.GetNodes();
      new_mesh_nodes->MakeOwner(mesh_fec);
      delete mesh_fespace;
    }
  }

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
  // This should work fine if some threads create an empty operator (no elements or bounday
  // elements).
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    CeedOperator loc_op, loc_op_t;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op_t));

    // Domain integrators first.
    if (!domain_integs.empty())
    {
      const int ne = mesh.GetNE();
      const int stride = (ne + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, ne);
      const bool use_bdr = false;

      const auto element_indices =
          GetElementIndices(trial_fespace, test_fespace, use_bdr, start, stop);

      for (const auto &value : element_indices)
      {
        const std::vector<int> &indices = value.second;
        const int q_order = fem::GetDefaultIntegrationOrder(
            trial_fespace, test_fespace, indices, use_bdr, q_extra_pk, q_extra_qk);
        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(mesh.GetElementGeometry(indices[0]), q_order);

        for (const auto &integ : domain_integs)
        {
          CeedOperator sub_op, sub_op_t;
          integ->Assemble(trial_fespace, test_fespace, ir, indices, ceed, &sub_op,
                          &sub_op_t);

          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
          if (sub_op_t)
          {
            PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op_t, sub_op_t));
            PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
          }
        }
      }
    }

    // Boundary integrators next.
    if (!boundary_integs.empty())
    {
      const int nbe = mesh.GetNBE();
      const int stride = (nbe + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, nbe);
      const bool use_bdr = true;

      const auto element_indices =
          GetElementIndices(trial_fespace, test_fespace, use_bdr, start, stop);

      for (const auto &value : element_indices)
      {
        const std::vector<int> &indices = value.second;
        const int q_order = fem::GetDefaultIntegrationOrder(
            trial_fespace, test_fespace, indices, use_bdr, q_extra_pk, q_extra_qk);
        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(mesh.GetBdrElementGeometry(indices[0]), q_order);

        for (const auto &integ : boundary_integs)
        {
          CeedOperator sub_op, sub_op_t;
          integ->AssembleBoundary(trial_fespace, test_fespace, ir, indices, ceed, &sub_op,
                                  &sub_op_t);

          PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op, sub_op));
          PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op));
          if (sub_op_t)
          {
            PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op_t, sub_op_t));
            PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
          }
        }
      }
    }

    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op_t));
    op->AddOper(loc_op, loc_op_t);  // Thread-safe
  }

  return op;
}

std::unique_ptr<mfem::SparseMatrix> BilinearForm::FullAssemble(const ceed::Operator &op,
                                                               bool skip_zeros) const
{
  return ceed::CeedOperatorFullAssemble(op, skip_zeros, false);
}

std::unique_ptr<ceed::Operator> DiscreteLinearOperator::Assemble() const
{
  // Construct dof multiplicity vector for scaling to account for dofs shared between
  // elements (on host, then copy to device).
  const auto &test_fespace = a.GetTestSpace();
  Vector test_multiplicity(test_fespace.GetVSize());
  test_multiplicity = 0.0;
  mfem::Array<int> dofs;
  auto h_mult = test_multiplicity.HostReadWrite();
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

  auto op = a.Assemble();
  op->SetDofMultiplicity(std::move(test_multiplicity));
  return op;
}

std::unique_ptr<mfem::SparseMatrix>
DiscreteLinearOperator::FullAssemble(const ceed::Operator &op, bool skip_zeros) const
{
  return ceed::CeedOperatorFullAssemble(op, skip_zeros, true);
}

}  // namespace palace
