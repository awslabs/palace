// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "interpolator.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

constexpr auto GSLIB_BB_TOL = 0.01;  // MFEM defaults, slightly reduced bounding box
constexpr auto GSLIB_NEWTON_TOL = 1.0e-12;

}  // namespace

InterpolationOperator::InterpolationOperator(const std::map<int, config::ProbeData> &probe,
                                             const Units &units,
                                             FiniteElementSpace &nd_space)
#if defined(MFEM_USE_GSLIB)
  : op(nd_space.GetParMesh().GetComm()), v_dim_fes(nd_space.Get().GetVectorDim())
{
  auto &mesh = nd_space.GetParMesh();
  // Set up probes interpolation. All processes search for all points.
  if (probe.empty())
  {
    return;
  }
  const int dim = mesh.SpaceDimension();
  MFEM_VERIFY(
      mesh.Dimension() == dim,
      "Probe postprocessing functionality requires mesh dimension == space dimension!");
  const int npts = static_cast<int>(probe.size());
  mfem::Vector xyz(npts * dim);
  op_idx.resize(npts);
  int i = 0;
  for (const auto &[idx, data] : probe)
  {
    for (int d = 0; d < dim; d++)
    {
      // Use default ordering byNODES.
      xyz(d * npts + i) = data.center[d];
    }
    op_idx[i++] = idx;
  }
  op.Setup(mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);
  op.FindPoints(xyz, mfem::Ordering::byNODES);
  op.SetDefaultInterpolationValue(0.0);
  i = 0;
  for (const auto &[idx, data] : probe)
  {
    if (op.GetCode()[i++] == 2)
    {
      Mpi::Warning(
          "Probe {:d} at ({:.3e}) m could not be found!\n Using default value 0.0!\n", idx,
          fmt::join(units.Dimensionalize<Units::ValueType::LENGTH>(data.center), ", "));
    }
  }

  // With reordered elements, the element which GSLIB selects for a point on the boundary
  // between elements (the first in the local element order) depends on the reordering and
  // on the partition: choose instead the element first in the order of the input mesh
  // among those containing the point (within a reference-space tolerance), on all
  // processes. Points which no element contains (found by GSLIB on the boundary of the
  // mesh) keep the GSLIB value.
  const auto keys = mesh::GetElementInputKeys(mesh);
  int use_keys = !keys.empty();
  Mpi::GlobalMin(1, &use_keys, mesh.GetComm());
  if (!use_keys)
  {
    return;
  }
  probe_mesh = &mesh;
  const int ne = mesh.GetNE();
  std::vector<double> bb(2 * dim * static_cast<std::size_t>(ne));
  {
    mfem::Array<int> v;
    for (int e = 0; e < ne; e++)
    {
      double *lo = &bb[2 * dim * e], *hi = lo + dim;
      std::fill(lo, lo + dim, mfem::infinity());
      std::fill(hi, hi + dim, -mfem::infinity());
      mesh.GetElementVertices(e, v);
      for (const int j : v)
      {
        const double *x = mesh.GetVertex(j);
        for (int d = 0; d < dim; d++)
        {
          lo[d] = std::min(lo[d], x[d]);
          hi[d] = std::max(hi[d], x[d]);
        }
      }
      // Inflated, since curved elements can extend beyond their vertices.
      for (int d = 0; d < dim; d++)
      {
        const double m = 0.1 * (hi[d] - lo[d]) + 1.0e-12 * std::abs(hi[d]);
        lo[d] -= m;
        hi[d] += m;
      }
    }
  }
  constexpr auto none = std::numeric_limits<std::int64_t>::max();
  std::vector<std::int64_t> k0(npts, none), k1(npts, none);
  std::vector<int> best_elem(npts, -1);
  std::vector<mfem::IntegrationPoint> best_ip(npts);
  mfem::Vector pt(dim);
  for (int p = 0; p < npts; p++)
  {
    if (op.GetCode()[p] == 2)
    {
      continue;
    }
    for (int d = 0; d < dim; d++)
    {
      pt(d) = xyz(d * npts + p);
    }
    for (int e = 0; e < ne; e++)
    {
      const double *lo = &bb[2 * dim * e], *hi = lo + dim;
      bool in = true;
      for (int d = 0; d < dim && in; d++)
      {
        in = (pt(d) >= lo[d] && pt(d) <= hi[d]);
      }
      if (!in || (keys[e][0] > k0[p] || (keys[e][0] == k0[p] && keys[e][1] >= k1[p])))
      {
        continue;
      }
      mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
      mfem::InverseElementTransformation inv(T);
      inv.SetSolverType(mfem::InverseElementTransformation::NewtonElementProject);
      inv.SetPhysicalRelTol(1.0e-14);
      inv.SetElementTol(1.0e-10);
      mfem::IntegrationPoint ip;
      if (inv.Transform(pt, ip) == mfem::InverseElementTransformation::Inside)
      {
        k0[p] = keys[e][0];
        k1[p] = keys[e][1];
        best_elem[p] = e;
        best_ip[p] = ip;
      }
    }
  }
  std::vector<std::int64_t> g0(k0), g1(npts);
  MPI_Allreduce(MPI_IN_PLACE, g0.data(), npts, MPI_INT64_T, MPI_MIN, mesh.GetComm());
  for (int p = 0; p < npts; p++)
  {
    g1[p] = (k0[p] == g0[p]) ? k1[p] : none;
  }
  MPI_Allreduce(MPI_IN_PLACE, g1.data(), npts, MPI_INT64_T, MPI_MIN, mesh.GetComm());
  is_canonical.assign(npts, 0);
  for (int p = 0; p < npts; p++)
  {
    if (g0[p] == none)
    {
      continue;
    }
    is_canonical[p] = 1;
    if (best_elem[p] >= 0 && k0[p] == g0[p] && k1[p] == g1[p])
    {
      canonical.push_back({p, best_elem[p], best_ip[p]});
    }
  }
}
#else
{
  MFEM_CONTRACT_VAR(GSLIB_BB_TOL);
  MFEM_CONTRACT_VAR(GSLIB_NEWTON_TOL);
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
  // Interpolated vector values are returned from GSLIB interpolator with the same ordering
  // as the source grid function, which we transform to byVDIM for output.
  const int npts = op.GetCode().Size();
  const int vdim = U.VectorDim();
  std::vector<double> vals(npts * vdim);
  if (U.FESpace()->GetOrdering() == mfem::Ordering::byVDIM)
  {
    mfem::Vector v(vals.data(), npts * vdim);
    op.Interpolate(U, v);
  }
  else
  {
    mfem::Vector v(npts * vdim);
    op.Interpolate(U, v);
    for (int d = 0; d < vdim; d++)
    {
      for (int i = 0; i < npts; i++)
      {
        vals[i * vdim + d] = v(d * npts + i);
      }
    }
  }
  if (!is_canonical.empty())
  {
    // The values of the probes on element boundaries in their canonical element, from the
    // process holding it.
    std::vector<double> cvals(npts * vdim, 0.0);
    mfem::Vector val(vdim);
    for (const auto &c : canonical)
    {
      mfem::ElementTransformation *T =
          const_cast<mfem::ParMesh *>(probe_mesh)->GetElementTransformation(c.elem);
      T->SetIntPoint(&c.ip);
      U.GetVectorValue(*T, c.ip, val);
      for (int d = 0; d < vdim; d++)
      {
        cvals[c.index * vdim + d] = val(d);
      }
    }
    Mpi::GlobalSum(npts * vdim, cvals.data(), probe_mesh->GetComm());
    for (int i = 0; i < npts; i++)
    {
      if (is_canonical[i])
      {
        std::copy_n(cvals.data() + i * vdim, vdim, vals.data() + i * vdim);
      }
    }
  }
  return vals;
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
