// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "feutils.hpp"
#include "utils/mfemoperators.hpp"

namespace palace
{

template <typename FECollection>
std::vector<std::unique_ptr<FECollection>> ConstructFECollections(bool pc_pmg, bool pc_lor,
                                                                  int p, int dim)
{
  // If the solver will use a LOR preconditioner, we need to construct with a specific basis
  // type.
  MFEM_VERIFY(p > 0, "Nedelec space order must be positive!");
  int b1 = mfem::BasisType::GaussLobatto, b2 = mfem::BasisType::GaussLegendre;
  if (pc_lor)
  {
    b2 = mfem::BasisType::IntegratedGLL;
  }
  std::vector<std::unique_ptr<FECollection>> fecs;

  constexpr bool is_vector_fec = std::is_same<FECollection, mfem::ND_FECollection>::value ||
                                 std::is_same<FECollection, mfem::RT_FECollection>::value;

  if (pc_pmg)
  {
    fecs.reserve(p);
    for (int o = 1; o <= p; o++)
    {
      if constexpr (is_vector_fec)
      {
        fecs.push_back(std::make_unique<FECollection>(o, dim, b1, b2));
      }
      else
      {
        fecs.push_back(std::make_unique<FECollection>(o, dim, b1));
      }
    }
  }
  else
  {
    fecs.reserve(1);
    if constexpr (is_vector_fec)
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1, b2));
    }
    else
    {
      fecs.push_back(std::make_unique<FECollection>(p, dim, b1));
    }
  }
  return fecs;
}

template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
    const std::vector<std::unique_ptr<FECollection>> &fecs,
    const mfem::Array<int> &dbc_marker)
{
  MFEM_VERIFY(!mesh.empty() && !fecs.empty(),
              "Empty mesh or FE collection for FE space construction!");
  auto *fespace = new mfem::ParFiniteElementSpace(mesh[0].get(), fecs[0].get());
  mfem::ParFiniteElementSpaceHierarchy fespaces(mesh[0].get(), fespace, false, true);
  // h-refinement
  for (std::size_t l = 1; l < mesh.size(); l++)
  {
    fespace = new mfem::ParFiniteElementSpace(mesh[l].get(), fecs[0].get());
    auto *P =
        new ZeroWrapTransferOperator(fespaces.GetFinestFESpace(), *fespace, dbc_marker);
    fespaces.AddLevel(mesh[l].get(), fespace, P, false, true, true);
  }
  // p-refinement
  for (std::size_t l = 1; l < fecs.size(); l++)
  {
    fespace = new mfem::ParFiniteElementSpace(mesh.back().get(), fecs[l].get());
    auto *P =
        new ZeroWrapTransferOperator(fespaces.GetFinestFESpace(), *fespace, dbc_marker);
    fespaces.AddLevel(mesh.back().get(), fespace, P, false, true, true);
  }
  return fespaces;
}

template <typename FECollection>
mfem::ParFiniteElementSpaceHierarchy
ConstructFiniteElementSpaceHierarchy(mfem::ParMesh &mesh, const FECollection &fec)
{
  auto *fespace = new mfem::ParFiniteElementSpace(&mesh, &fec);
  return mfem::ParFiniteElementSpaceHierarchy(&mesh, fespace, false, true);
}

// explicit instantiations
template std::vector<std::unique_ptr<mfem::ND_FECollection>>
ConstructFECollections(bool, bool, int, int);
template std::vector<std::unique_ptr<mfem::H1_FECollection>>
ConstructFECollections(bool, bool, int, int);

template mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    const std::vector<std::unique_ptr<mfem::ParMesh>> &,
    const std::vector<std::unique_ptr<mfem::ND_FECollection>> &, const mfem::Array<int> &);
template mfem::ParFiniteElementSpaceHierarchy ConstructFiniteElementSpaceHierarchy(
    const std::vector<std::unique_ptr<mfem::ParMesh>> &,
    const std::vector<std::unique_ptr<mfem::H1_FECollection>> &, const mfem::Array<int> &);

template mfem::ParFiniteElementSpaceHierarchy
ConstructFiniteElementSpaceHierarchy(mfem::ParMesh &, const mfem::ND_FECollection &);
template mfem::ParFiniteElementSpaceHierarchy
ConstructFiniteElementSpaceHierarchy(mfem::ParMesh &, const mfem::H1_FECollection &);

}  // namespace palace
