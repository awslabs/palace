// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "bilinearform.hpp"

#include <array>
#include <unordered_map>
#include <ceed.h>
#include "fem/libceed/operator.hpp"
#include "fem/libceed/utils.hpp"
#include "linalg/vector.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif

namespace palace
{

namespace
{

using ElementKey = std::array<int, 3>;

struct ElementHash
{
  std::size_t operator()(const ElementKey &k) const
  {
    return ceed::CeedHashCombine(
        ceed::CeedHashCombine(ceed::CeedHash(k[0]), ceed::CeedHash(k[1])),
        ceed::CeedHash(k[2]));
  }
};

// Count the number of elements of each type in the local mesh.
std::unordered_map<ElementKey, std::vector<int>, ElementHash>
GetElementIndices(const mfem::FiniteElementSpace &trial_fespace,
                  const mfem::FiniteElementSpace &test_fespace, bool use_bdr, int start,
                  int stop)
{
  mfem::Mesh &mesh = *trial_fespace.GetMesh();
  std::unordered_map<ElementKey, int, ElementHash> counts, offsets;
  std::unordered_map<ElementKey, std::vector<int>, ElementHash> element_indices;

  // Count the number of elements of each type and order.
  for (int i = start; i < stop; i++)
  {
    const mfem::FiniteElement &trial_fe =
        use_bdr ? *trial_fespace.GetBE(i) : *trial_fespace.GetFE(i);
    const mfem::FiniteElement &test_fe =
        use_bdr ? *test_fespace.GetBE(i) : *test_fespace.GetFE(i);
    mfem::Element::Type type = use_bdr ? mesh.GetBdrElementType(i) : mesh.GetElementType(i);
    ElementKey key = {type, trial_fe.GetOrder(), test_fe.GetOrder()};
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
    mfem::Element::Type type = use_bdr ? mesh.GetBdrElementType(i) : mesh.GetElementType(i);
    ElementKey key = {type, trial_fe.GetOrder(), test_fe.GetOrder()};
    int &offset = offsets[key];
    std::vector<int> &indices = element_indices[key];
    indices[offset++] = i;
  }

  return element_indices;
}

}  // namespace

std::unique_ptr<Operator> BilinearForm::Assemble() const
{
  MFEM_VERIFY(trial_fespace.GetMesh() == test_fespace.GetMesh(),
              "Trial and test finite element spaces must correspond to the same mesh!");
  mfem::Mesh &mesh = *trial_fespace.GetMesh();
  mesh.EnsureNodes();
  if (q_order < 0)
  {
    q_order = fem::GetDefaultIntegrationOrder(trial_fespace, test_fespace);
  }
  std::unique_ptr<ceed::Operator> op;
  if (&trial_fespace == &test_fespace)
  {
    op = std::make_unique<ceed::SymmetricOperator>();
  }
  else
  {
    op = std::make_unique<ceed::Operator>();
  }

  // Assemble the libCEED operator in parallel, each thread builds a composite operator.
  // This should work fine if some threads create an empty operator (no elements or bounday
  // elements).
  PalacePragmaOmp(parallel)
  {
#if defined(MFEM_USE_OPENMP)
    const int nt = omp_get_num_threads();
    const int tid = omp_get_thread_num();
#else
    const int nt = 1;
    const int tid = 0;
#endif
    Ceed ceed = ceed::internal::ceed[tid];
    CeedOperator loc_op, loc_op_t = nullptr;
    PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op));

    // Domain integrators first.
    if (!domain_integs.empty())
    {
      const int ne = mesh.GetNE();
      const int stride = (ne + nt - 1) / nt;
      const int start = tid * stride;
      const int stop = std::min(start + stride, ne);
      const bool use_bdr = false;

      const auto element_indices =
          GetElementIndices(trial_fespace, test_fespace, use_bdr, start, stop);

      for (const auto &value : element_indices)
      {
        const std::vector<int> &indices = value.second;
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
            if (!loc_op_t)
            {
              PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op_t));
            }
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
      const int start = tid * stride;
      const int stop = std::min(start + stride, nbe);
      const bool use_bdr = true;

      const auto element_indices =
          GetElementIndices(trial_fespace, test_fespace, use_bdr, start, stop);

      for (const auto &value : element_indices)
      {
        const std::vector<int> &indices = value.second;
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
            if (!loc_op_t)
            {
              PalaceCeedCall(ceed, CeedCompositeOperatorCreate(ceed, &loc_op_t));
            }
            PalaceCeedCall(ceed, CeedCompositeOperatorAddSub(loc_op_t, sub_op_t));
            PalaceCeedCall(ceed, CeedOperatorDestroy(&sub_op_t));
          }
        }
      }
    }

    PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op));
    if (loc_op_t)
    {
      PalaceCeedCall(ceed, CeedOperatorCheckReady(loc_op_t));
    }
    op->AddOper(loc_op, loc_op_t);  // Thread-safe
  }

  return op;
}

std::unique_ptr<mfem::SparseMatrix> BilinearForm::FullAssemble(bool skip_zeros) const
{
  auto op = Assemble();
  return ceed::CeedOperatorFullAssemble(*static_cast<ceed::Operator *>(op.get()),
                                        skip_zeros, false);
}

std::unique_ptr<Operator> DiscreteLinearOperator::Assemble() const
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
  static_cast<ceed::Operator *>(op.get())->SetDofMultiplicity(std::move(test_multiplicity));
  return op;
}

std::unique_ptr<mfem::SparseMatrix>
DiscreteLinearOperator::FullAssemble(bool skip_zeros) const
{
  auto op = a.Assemble();
  return ceed::CeedOperatorFullAssemble(*static_cast<ceed::Operator *>(op.get()),
                                        skip_zeros, true);
}

}  // namespace palace
