// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fespace.hpp"

#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/libceed/utils.hpp"
#include "linalg/rap.hpp"
#include "utils/omp.hpp"

namespace palace
{

namespace
{

// Count the number of elements of each type in the local mesh.
std::unordered_map<mfem::Geometry::Type, std::vector<int>>
GetElementIndices(const mfem::ParMesh &mesh, bool use_bdr, int start, int stop)
{
  std::unordered_map<mfem::Geometry::Type, int> counts, offsets;
  std::unordered_map<mfem::Geometry::Type, std::vector<int>> element_indices;

  for (int i = start; i < stop; i++)
  {
    const auto key = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto it = counts.find(key);
    if (it == counts.end())
    {
      counts[key] = 1;
    }
    else
    {
      it->second++;
    }
  }

  // Populate the indices arrays for each element geometry.
  for (const auto &[key, val] : counts)
  {
    offsets[key] = 0;
    element_indices[key] = std::vector<int>(val);
  }
  for (int i = start; i < stop; i++)
  {
    const auto key = use_bdr ? mesh.GetBdrElementGeometry(i) : mesh.GetElementGeometry(i);
    auto &offset = offsets[key];
    auto &indices = element_indices[key];
    indices[offset++] = i;
  }

  return element_indices;
}

ceed::CeedGeomFactorData AssembleGeometryData(const mfem::GridFunction &mesh_nodes,
                                              Ceed ceed, mfem::Geometry::Type geom,
                                              std::vector<int> &indices)
{
  const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
  const mfem::Mesh &mesh = *mesh_fespace.GetMesh();

  // XX TODO: In practice we do not need to compute all of the geometry data for every
  //          simulation type.
  auto data = ceed::CeedGeomFactorDataCreate(ceed);
  data->dim = mfem::Geometry::Dimension[geom];
  data->space_dim = mesh.SpaceDimension();

  // Compute the required geometry factors at quadrature points.
  constexpr auto info = ceed::GeomFactorInfo::Determinant | ceed::GeomFactorInfo::Adjugate |
                        ceed::GeomFactorInfo::Jacobian;
  CeedElemRestriction mesh_restr =
      FiniteElementSpace::BuildCeedElemRestriction(mesh_fespace, ceed, geom, indices);
  CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(mesh_fespace, ceed, geom);
  CeedInt nqpts;
  PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &nqpts));
  ceed::AssembleCeedGeometryData(info, ceed, mesh_restr, mesh_basis, mesh_nodes, data);
  PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
  PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));

  // Compute element attribute quadrature data. This should ideally be a single scalar per
  // element but all fields associated with a CeedOperator require the same number of
  // quadrature points.

  // XX TODO COMPRESS TO ONLY USED ATTRIBUTES FROM INDICES (this rank, this thread?)

  {
    const auto ne = indices.size();
    const bool use_bdr = (data->dim != mesh.Dimension());
    data->attr.SetSize(ne * nqpts);
    for (std::size_t j = 0; j < ne; j++)
    {
      const int e = indices[j];
      const int attr = use_bdr ? mesh.GetBdrAttribute(e) : mesh.GetAttribute(e);
      for (CeedInt q = 0; q < nqpts; q++)
      {
        data->attr[j * nqpts + q] = attr;
      }
    }
    PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(ceed, ne, nqpts, 1, ne * nqpts,
                                                          CEED_STRIDES_BACKEND,
                                                          &data->attr_restr));
    ceed::InitCeedVector(data->attr, ceed, &data->attr_vec);
  }

  // Save mesh element indices.
  data->indices = std::move(indices);

  return data;
}

}  // namespace

namespace fem
{

ceed::CeedObjectMap<ceed::CeedGeomFactorData> SetUpCeedGeomFactorData(mfem::ParMesh &mesh)
{
  mesh.EnsureNodes();
  const mfem::GridFunction &mesh_nodes = *mesh.GetNodes();

  // The geometry factor object to be constructed.
  ceed::CeedObjectMap<ceed::CeedGeomFactorData> geom_data;

  // Create a list of the element indices in the mesh corresponding to a given thread and
  // element geometry type. libCEED operators will be constructed in parallel over threads,
  // where each thread builds a composite operator with sub-operators for each geometry.
  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];

    // First domain elements.
    {
      const int ne = mesh.GetNE();
      const int stride = (ne + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, ne);
      constexpr bool use_bdr = false;
      auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data =
            AssembleGeometryData(mesh_nodes, ceed, geom, indices);
        PalacePragmaOmp(critical(GeomFactorData))
        {
          geom_data.emplace(std::make_pair(ceed, geom), std::move(data));
        }
      }
    }

    // Then boundary elements.
    {
      const int nbe = mesh.GetNBE();
      const int stride = (nbe + nt - 1) / nt;
      const int start = i * stride;
      const int stop = std::min(start + stride, nbe);
      constexpr bool use_bdr = true;
      auto element_indices = GetElementIndices(mesh, use_bdr, start, stop);
      for (auto &[geom, indices] : element_indices)
      {
        ceed::CeedGeomFactorData data =
            AssembleGeometryData(mesh_nodes, ceed, geom, indices);
        PalacePragmaOmp(critical(GeomFactorData))
        {
          geom_data.emplace(std::make_pair(ceed, geom), std::move(data));
        }
      }
    }
  }

  return geom_data;
}

void UpdateCeedGeomFactorData(const mfem::ParMesh &mesh,
                              ceed::CeedObjectMap<ceed::CeedGeomFactorData> &geom_data)
{
  // Regenerate the geometry factor data for each element type and replace the entries in
  // the table.
  const mfem::GridFunction &mesh_nodes = *mesh.GetNodes();

  const std::size_t nt = ceed::internal::GetCeedObjects().size();
  PalacePragmaOmp(parallel for schedule(static))
  for (std::size_t i = 0; i < nt; i++)
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[i];

    for (auto &[key, val] : geom_data)
    {
      if (key.first != ceed)
      {
        continue;
      }
      const auto geom = key.second;
      const auto &orig_data = val;

      geom_data[key] = AssembleGeometryData(mesh_nodes, ceed, geom, orig_data->indices);
    }
  }
}

}  // namespace fem

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
  const mfem::FiniteElement &fe = *FEColl()->FiniteElementForGeometry(geom);
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
  const mfem::FiniteElement &fe = *FEColl()->FiniteElementForGeometry(geom);
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
  const int dim = GetMesh()->Dimension();
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
