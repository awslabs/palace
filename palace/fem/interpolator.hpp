// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTERPOLATOR_HPP
#define PALACE_FEM_INTERPOLATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "utils/configfile.hpp"

namespace palace
{

class GridFunction;
class IoData;
class Units;
class FiniteElementSpace;

//
// A class which wraps MFEM's GSLIB interface for high-order field interpolation.
//
class InterpolationOperator
{
private:
#if defined(MFEM_USE_GSLIB)
  mfem::FindPointsGSLIB op;
#endif
  std::vector<int> op_idx;

  int v_dim_fes;  // dimension of interpolated vector from NDSpace

  // With domains, the submesh of their elements, in which the probes are found and the
  // fields interpolated, and for each space of a field, its counterpart on the submesh with
  // a field there and the transfer to it.
  std::unique_ptr<mfem::ParSubMesh> submesh;
  struct SubmeshField
  {
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::ParGridFunction> field;
    std::unique_ptr<mfem::ParTransferMap> transfer;
  };
  std::map<const mfem::ParFiniteElementSpace *, SubmeshField> submesh_fields;

  std::vector<double> ProbeField(const mfem::ParGridFunction &U);

public:
  // With domains (attributes), the probes are found, and the fields interpolated, in the
  // elements of the domains only: a probe on a face of one of them is found in one of them.
  InterpolationOperator(const std::map<int, config::ProbeData> &probe, const Units &units,
                        FiniteElementSpace &nd_space, const std::vector<int> &domains = {});
  InterpolationOperator(const IoData &iodata, FiniteElementSpace &nd_space);

  auto GetVDim() const { return v_dim_fes; }
  const auto &GetProbes() const { return op_idx; }

  std::vector<std::complex<double>> ProbeField(const GridFunction &U);
};

namespace fem
{

// The indices of the probes in the elements with the given attributes, as an
// InterpolationOperator with these domains finds them (also those on a face of one of these
// elements). Collective.
std::vector<int> ProbesInDomains(const std::map<int, config::ProbeData> &probe,
                                 const mfem::ParMesh &mesh,
                                 const std::vector<int> &attributes);

// Interpolate a function on a serial or parallel mesh to a different mesh, using GSLIB.
// Similar to MFEM's field-interp miniapp.
void InterpolateFunction(const mfem::GridFunction &U, mfem::GridFunction &V);

// Interpolate a function at a specific list of points, specified using the provided
// ordering. The output vector values are always arranged byVDIM.
void InterpolateFunction(const mfem::Vector &xyz, const mfem::GridFunction &U,
                         mfem::Vector &V,
                         mfem::Ordering::Type ordering = mfem::Ordering::byNODES);

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_INTERPOLATOR_HPP
