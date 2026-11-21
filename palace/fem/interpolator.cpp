// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "interpolator.hpp"

#include <algorithm>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

constexpr auto GSLIB_BB_TOL = 0.01;  // MFEM defaults, slightly reduced bounding box
constexpr auto GSLIB_NEWTON_TOL = 1.0e-12;

#if defined(MFEM_USE_GSLIB)
// Find the probes with GSLIB, all processes searching for all points: their indices, in the
// order of the points.
std::vector<int> FindProbes(mfem::FindPointsGSLIB &op,
                            const std::map<int, config::ProbeData> &probe,
                            mfem::ParMesh &mesh)
{
  const int dim = mesh.SpaceDimension();
  MFEM_VERIFY(
      mesh.Dimension() == dim,
      "Probe postprocessing functionality requires mesh dimension == space dimension!");
  const int npts = static_cast<int>(probe.size());
  mfem::Vector xyz(npts * dim);
  std::vector<int> idx(npts);
  int i = 0;
  for (const auto &[index, data] : probe)
  {
    for (int d = 0; d < dim; d++)
    {
      // Use default ordering byNODES.
      xyz(d * npts + i) = data.center[d];
    }
    idx[i++] = index;
  }
  op.Setup(mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);
  op.FindPoints(xyz, mfem::Ordering::byNODES);
  return idx;
}

// The submesh of the elements with the given attributes, with nodes for the point search.
mfem::ParSubMesh DomainSubmesh(const mfem::ParMesh &mesh,
                               const std::vector<int> &attributes)
{
  mfem::Array<int> attr(const_cast<int *>(attributes.data()),
                        static_cast<int>(attributes.size()));
  auto submesh = mfem::ParSubMesh::CreateFromDomain(mesh, attr);
  submesh.EnsureNodes();
  return submesh;
}
#endif

}  // namespace

InterpolationOperator::InterpolationOperator(const std::map<int, config::ProbeData> &probe,
                                             const Units &units,
                                             FiniteElementSpace &nd_space,
                                             const std::vector<int> &domains)
#if defined(MFEM_USE_GSLIB)
  : op(nd_space.GetParMesh().GetComm()), v_dim_fes(nd_space.Get().GetVectorDim())
{
  // Set up probes interpolation. All processes search for all points.
  if (probe.empty())
  {
    return;
  }
  if (!domains.empty())
  {
    submesh =
        std::make_unique<mfem::ParSubMesh>(DomainSubmesh(nd_space.GetParMesh(), domains));
  }
  op_idx = FindProbes(op, probe, submesh ? *submesh : nd_space.GetParMesh());
  op.SetDefaultInterpolationValue(0.0);
  int i = 0;
  for (const auto &[idx, data] : probe)
  {
    if (op.GetCode()[i++] == 2)
    {
      Mpi::Warning(
          "Probe {:d} at ({:.3e}) m could not be found!\n Using default value 0.0!\n", idx,
          fmt::join(units.Dimensionalize<Units::ValueType::LENGTH>(data.center), ", "));
    }
  }
}
#else
{
  MFEM_CONTRACT_VAR(GSLIB_BB_TOL);
  MFEM_CONTRACT_VAR(GSLIB_NEWTON_TOL);
  MFEM_CONTRACT_VAR(units);
  MFEM_CONTRACT_VAR(nd_space);
  MFEM_CONTRACT_VAR(domains);
  MFEM_VERIFY(probe.empty(), "InterpolationOperator class requires MFEM_USE_GSLIB!");
}
#endif

InterpolationOperator::InterpolationOperator(const IoData &iodata,
                                             FiniteElementSpace &nd_space)
  : InterpolationOperator(iodata.domains.postpro.probe, iodata.units, nd_space)
{
}

std::vector<double> InterpolationOperator::ProbeField(const mfem::ParGridFunction &U)
{
#if defined(MFEM_USE_GSLIB)
  // On the submesh of the domains, the field transferred there.
  const mfem::ParGridFunction *V = &U;
  if (submesh)
  {
    MFEM_ASSERT(U.ParFESpace()->GetParMesh() == submesh->GetParent(),
                "A probed field must be on the mesh of the probes!");
    auto &s = submesh_fields[U.ParFESpace()];
    if (!s.fespace)
    {
      const auto &fes = *U.ParFESpace();
      s.fespace = std::make_unique<mfem::ParFiniteElementSpace>(
          submesh.get(), fes.FEColl(), fes.GetVDim(), fes.GetOrdering());
      s.field = std::make_unique<mfem::ParGridFunction>(s.fespace.get());
      s.transfer = std::make_unique<mfem::ParTransferMap>(fes, *s.fespace);
    }
    s.transfer->Transfer(U, *s.field);
    V = s.field.get();
  }
  // Interpolate with byVDIM output ordering. GSLIB moves the result to the device when the
  // field is on the device, so read it back on the host.
  const int npts = op.GetCode().Size();
  const int vdim = V->VectorDim();
  mfem::Vector v(npts * vdim);
  op.Interpolate(*V, v, mfem::Ordering::byVDIM);
  const double *h_v = v.HostRead();
  return {h_v, h_v + npts * vdim};
#else
  MFEM_ABORT("InterpolationOperator class requires MFEM_USE_GSLIB!");
  return {};
#endif
}

std::vector<std::complex<double>> InterpolationOperator::ProbeField(const GridFunction &U)
{
  std::vector<double> vr = ProbeField(U.Real());
  if (U.HasImag())
  {
    std::vector<double> vi = ProbeField(U.Imag());
    std::vector<std::complex<double>> vals(vr.size());
    std::ranges::transform(vr, vi, vals.begin(), [](double xr, double xi)
                           { return std::complex<double>(xr, xi); });
    return vals;
  }
  else
  {
    return {vr.begin(), vr.end()};
  }
}

namespace fem
{

std::vector<int> ProbesInDomains(const std::map<int, config::ProbeData> &probe,
                                 const mfem::ParMesh &mesh,
                                 const std::vector<int> &attributes)
{
  std::vector<int> in_domains;
  if (probe.empty() || attributes.empty())
  {
    return in_domains;
  }
#if defined(MFEM_USE_GSLIB)
  MPI_Comm comm = mesh.GetComm();
  auto submesh = DomainSubmesh(mesh, attributes);
  mfem::FindPointsGSLIB op(comm);
  const auto idx = FindProbes(op, probe, submesh);
  // Each process searches for all points, but the probes are those of the root's search.
  const int n = static_cast<int>(idx.size());
  std::vector<int> found(n, 0);
  for (int i = 0; Mpi::Root(comm) && i < n; i++)
  {
    found[i] = (op.GetCode()[i] != 2);
  }
  Mpi::Broadcast(n, found.data(), 0, comm);
  for (int i = 0; i < n; i++)
  {
    if (found[i])
    {
      in_domains.push_back(idx[i]);
    }
  }
#else
  MFEM_CONTRACT_VAR(mesh);
  MFEM_CONTRACT_VAR(attributes);
  MFEM_ABORT("Probes require MFEM_USE_GSLIB!");
#endif
  return in_domains;
}

void InterpolateFunction(const mfem::GridFunction &U, mfem::GridFunction &V)
{
#if defined(MFEM_USE_GSLIB)
  // Generate list of points where the grid function will be evaluated. If the grid function
  // to interpolate is an H1 space of the same order as the mesh nodes, we can use the
  // mesh node points directly. Otherwise, for a different basis order or type, we generate
  // the interpolation points in the physical space manually.
  auto &dest_mesh = *V.FESpace()->GetMesh();
  MFEM_VERIFY(dest_mesh.GetNodes(), "Destination mesh has no nodal FE space!");
  const int dim = dest_mesh.SpaceDimension();
  mfem::Vector xyz;
  mfem::Ordering::Type ordering;
  const auto *dest_fec_h1 =
      dynamic_cast<const mfem::H1_FECollection *>(V.FESpace()->FEColl());
  const auto *dest_nodes_h1 = dynamic_cast<const mfem::H1_FECollection *>(
      dest_mesh.GetNodes()->FESpace()->FEColl());
  int dest_fespace_order = V.FESpace()->GetMaxElementOrder();
  int dest_nodes_order = dest_mesh.GetNodes()->FESpace()->GetMaxElementOrder();
  dest_mesh.GetNodes()->HostRead();
  if (dest_fec_h1 && dest_nodes_h1 && dest_fespace_order == dest_nodes_order)
  {
    xyz.MakeRef(*dest_mesh.GetNodes(), 0, dest_mesh.GetNodes()->Size());
    ordering = dest_mesh.GetNodes()->FESpace()->GetOrdering();
  }
  else
  {
    int npts = 0, offset = 0;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      npts += V.FESpace()->GetFE(i)->GetNodes().GetNPoints();
    }
    xyz.SetSize(npts * dim);
    mfem::DenseMatrix pointmat;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
      mfem::ElementTransformation &T = *dest_mesh.GetElementTransformation(i);
      T.Transform(fe.GetNodes(), pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        for (int d = 0; d < dim; d++)
        {
          // Use default ordering byNODES.
          xyz(d * npts + offset + j) = pointmat(d, j);
        }
      }
      offset += pointmat.Width();
    }
    ordering = mfem::Ordering::byNODES;
  }
  const int npts = xyz.Size() / dim;

  // Set up the interpolator.
  auto &src_mesh = *U.FESpace()->GetMesh();
  MFEM_VERIFY(src_mesh.GetNodes(), "Source mesh has no nodal FE space!");
  auto *src_pmesh = dynamic_cast<mfem::ParMesh *>(&src_mesh);
  MPI_Comm comm = (src_pmesh) ? src_pmesh->GetComm() : MPI_COMM_SELF;
  mfem::FindPointsGSLIB op(comm);
  op.Setup(src_mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);

  // Perform the interpolation and fill the target GridFunction (see MFEM's field-interp
  // miniapp).
  const int vdim = U.VectorDim();
  mfem::Vector vals(npts * vdim);
  op.SetDefaultInterpolationValue(0.0);
  op.SetL2AvgType(mfem::FindPointsGSLIB::NONE);
  op.Interpolate(xyz, U, vals, ordering);
  const auto *dest_fec_l2 =
      dynamic_cast<const mfem::L2_FECollection *>(V.FESpace()->FEColl());
  if (dest_fec_h1 || dest_fec_l2)
  {
    if (dest_fec_h1 && dest_fespace_order != dest_nodes_order)
    {
      // H1 with order != mesh order needs to handle duplicated interpolation points.
      mfem::Vector elem_vals;
      mfem::Array<int> vdofs;
      int offset = 0;
      for (int i = 0; i < dest_mesh.GetNE(); i++)
      {
        const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
        const int elem_npts = fe.GetNodes().GetNPoints();
        elem_vals.SetSize(elem_npts * vdim);
        for (int d = 0; d < vdim; d++)
        {
          for (int j = 0; j < elem_npts; j++)
          {
            // Arrange element values byNODES to align with GetElementVDofs.
            int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                          ? d * npts + offset + j
                          : (offset + j) * vdim + d;
            elem_vals(d * elem_npts + j) = vals(idx);
          }
        }
        const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
        if (dof_trans)
        {
          dof_trans->TransformPrimal(elem_vals);
        }
        V.SetSubVector(vdofs, elem_vals);
        offset += elem_npts;
      }
    }
    else
    {
      // Otherwise, H1 and L2 copy interpolated values to vdofs.
      MFEM_ASSERT(V.Size() == vals.Size(),
                  "Unexpected size mismatch for interpolated values and grid function!");
      V = vals;
    }
  }
  else
  {
    // H(div) or H(curl) use ProjectFromNodes.
    mfem::Vector elem_vals, v;
    mfem::Array<int> vdofs;
    int offset = 0;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
      mfem::ElementTransformation &T = *dest_mesh.GetElementTransformation(i);
      const int elem_npts = fe.GetNodes().GetNPoints();
      elem_vals.SetSize(elem_npts * vdim);
      for (int d = 0; d < vdim; d++)
      {
        for (int j = 0; j < elem_npts; j++)
        {
          // Arrange element values byVDIM for ProjectFromNodes.
          int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                        ? d * npts + offset + j
                        : (offset + j) * vdim + d;
          elem_vals(j * vdim + d) = vals(idx);
        }
      }
      const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
      v.SetSize(vdofs.Size());
      fe.ProjectFromNodes(elem_vals, T, v);
      if (dof_trans)
      {
        dof_trans->TransformPrimal(v);
      }
      V.SetSubVector(vdofs, v);
      offset += elem_npts;
    }
  }
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

void InterpolateFunction(const mfem::Vector &xyz, const mfem::GridFunction &U,
                         mfem::Vector &vals, mfem::Ordering::Type ordering)
{
#if defined(MFEM_USE_GSLIB)
  // Set up the interpolator.
  auto &src_mesh = *U.FESpace()->GetMesh();
  MFEM_VERIFY(src_mesh.GetNodes(), "Source mesh has no nodal FE space!");
  const int dim = src_mesh.SpaceDimension();
  const int npts = xyz.Size() / dim;
  auto *src_pmesh = dynamic_cast<mfem::ParMesh *>(&src_mesh);
  MPI_Comm comm = (src_pmesh) ? src_pmesh->GetComm() : MPI_COMM_SELF;
  mfem::FindPointsGSLIB op(comm);
  op.Setup(src_mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);

  // Perform the interpolation, with the ordering of the returned values matching the
  // ordering of the source grid function.
  const int vdim = U.VectorDim();
  MFEM_VERIFY(vals.Size() == npts * vdim, "Incorrect size for interpolated values vector!");
  op.SetDefaultInterpolationValue(0.0);
  op.SetL2AvgType(mfem::FindPointsGSLIB::NONE);
  op.Interpolate(xyz, U, vals, ordering);
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

}  // namespace fem

}  // namespace palace
