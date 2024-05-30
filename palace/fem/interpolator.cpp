// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "interpolator.hpp"

#include <algorithm>
#include "fem/gridfunction.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

#if defined(MFEM_USE_GSLIB)
InterpolationOperator::InterpolationOperator(const IoData &iodata, mfem::ParMesh &mesh)
  : op(mesh.GetComm())
#else
InterpolationOperator::InterpolationOperator(const IoData &iodata, mfem::ParMesh &mesh)
#endif
{
#if defined(MFEM_USE_GSLIB)
  // Set up probes interpolation. All processes search for all points.
  if (iodata.domains.postpro.probe.empty())
  {
    return;
  }
  constexpr double bb_t = 0.1;  // MFEM defaults
  constexpr double newton_tol = 1.0e-12;
  const int dim = mesh.SpaceDimension();
  MFEM_VERIFY(
      mesh.Dimension() == dim,
      "Probe postprocessing functionality requires mesh dimension == space dimension!");
  const int npts = static_cast<int>(iodata.domains.postpro.probe.size());
  mfem::Vector xyz(npts * dim);
  op_idx.resize(npts);
  int i = 0;
  for (const auto &[idx, data] : iodata.domains.postpro.probe)
  {
    // Use default ordering byNODES.
    xyz(i) = data.center[0];
    xyz(npts + i) = data.center[1];
    if (dim == 3)
    {
      xyz(2 * npts + i) = data.center[2];
    }
    op_idx[i++] = idx;
  }
  op.Setup(mesh, bb_t, newton_tol, npts);
  op.FindPoints(xyz, mfem::Ordering::byNODES);
  op.SetDefaultInterpolationValue(0.0);
  i = 0;
  for (const auto &[idx, data] : iodata.domains.postpro.probe)
  {
    if (op.GetCode()[i++] == 2)
    {
      Mpi::Warning("Probe {:d} at ({:.3e}, {:.3e}, {:.3e}) m could not be found!\n"
                   "Using default value 0.0!\n",
                   idx,
                   iodata.DimensionalizeValue(IoData::ValueType::LENGTH, data.center[0]),
                   iodata.DimensionalizeValue(IoData::ValueType::LENGTH, data.center[1]),
                   iodata.DimensionalizeValue(IoData::ValueType::LENGTH, data.center[2]));
    }
  }
#else
  MFEM_VERIFY(iodata.domains.postpro.probe.empty(),
              "InterpolationOperator class requires MFEM_USE_GSLIB!");
#endif
}

std::vector<double> InterpolationOperator::ProbeField(const mfem::ParGridFunction &U)
{
#if defined(MFEM_USE_GSLIB)
  // Interpolated vector values are returned from GSLIB interpolator byNODES, which we
  // transform to byVDIM for output.
  const int npts = op.GetCode().Size();
  const int vdim = U.VectorDim();
  std::vector<double> vals(npts * vdim);
  mfem::Vector v(npts * vdim);
  op.Interpolate(U, v);
  for (int d = 0; d < vdim; d++)
  {
    for (int i = 0; i < npts; i++)
    {
      vals[i * vdim + d] = v(d * npts + i);
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
    std::transform(vr.begin(), vr.end(), vi.begin(), vals.begin(),
                   [](double xr, double xi) { return std::complex<double>(xr, xi); });
    return vals;
  }
  else
  {
    return std::vector<std::complex<double>>(vr.begin(), vr.end());
  }
}

namespace fem
{

void InterpolateFunction(const mfem::GridFunction &U, mfem::GridFunction &V)
{
#if defined(MFEM_USE_GSLIB)
  // Generate list of points where the grid function will be evaluated.
  auto &dest_mesh = *V.FESpace()->GetMesh();
  MFEM_VERIFY(dest_mesh.GetNodes(), "Destination mesh has no nodal FE space!");
  const int dim = dest_mesh.SpaceDimension();
  mfem::Vector xyz;
  mfem::Ordering::Type ordering;
  const auto *dest_fec_h1 =
      dynamic_cast<const mfem::H1_FECollection *>(V.FESpace()->FEColl());
  const auto *dest_fec_l2 =
      dynamic_cast<const mfem::L2_FECollection *>(V.FESpace()->FEColl());
  const auto *dest_nodes_h1 = dynamic_cast<const mfem::H1_FECollection *>(
      dest_mesh.GetNodes()->FESpace()->FEColl());
  const auto *dest_nodes_l2 = dynamic_cast<const mfem::L2_FECollection *>(
      dest_mesh.GetNodes()->FESpace()->FEColl());
  int dest_fespace_order = V.FESpace()->GetMaxElementOrder();
  int dest_nodes_order = dest_mesh.GetNodes()->FESpace()->GetMaxElementOrder();
  dest_mesh.GetNodes()->HostRead();
  if (((dest_fec_h1 && dest_nodes_h1) || (dest_fec_l2 && dest_nodes_l2)) &&
      dest_fespace_order == dest_nodes_order)
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
      mfem::ElementTransformation &T = *V.FESpace()->GetElementTransformation(i);
      T.Transform(fe.GetNodes(), pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        // Fill with ordering byNODES.
        xyz(offset + j) = pointmat(0, j);
        xyz(npts + offset + j) = pointmat(1, j);
        if (dim == 3)
        {
          xyz(2 * npts + offset + j) = pointmat(2, j);
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
  constexpr double bb_t = 0.1;  // MFEM defaults
  constexpr double newton_tol = 1.0e-12;
  op.Setup(src_mesh, bb_t, newton_tol, npts);

  // Perform the interpolation and fill the target GridFunction (see MFEM's field-interp
  // miniapp).
  const int vdim = V.VectorDim();
  mfem::Vector vals(npts * vdim);
  op.SetDefaultInterpolationValue(0.0);
  op.SetL2AvgType(mfem::FindPointsGSLIB::NONE);
  op.Interpolate(xyz, U, vals, ordering);
  if (dest_fec_h1 || dest_fec_l2)
  {
    if ((dest_fec_h1 && dest_fespace_order == dest_nodes_order) || dest_fec_l2)
    {
      // H1 and L2 (copy interpolated values to vdofs).
      MFEM_ASSERT(V.Size() == vals.Size(),
                  "Unexpected size mismatch for interpolated values and grid function!");
      V = vals;
    }
    else
    {
      // H1 with order != mesh order needs to handle duplicated points.
      mfem::Vector elem_vals;
      mfem::Array<int> vdofs;
      int offset = 0;
      for (int i = 0; i < dest_mesh.GetNE(); i++)
      {
        const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
        const int elem_npts = fe.GetNodes().GetNPoints();
        elem_vals.SetSize(elem_npts * vdim);

        for (int j = 0; j < elem_npts; j++)
        {
          for (int d = 0; d < vdim; d++)
          {
            // Arrange values byNODES to align with GetElementVDofs.
            int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                          ? d * npts + offset + j
                          : (offset + j) * dim + d;
            elem_vals(d * elem_npts + j) = vals(idx);
          }
        }
        offset += elem_npts;
        const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
        if (dof_trans)
        {
          dof_trans->TransformPrimal(elem_vals);
        }
        V.SetSubVector(vdofs, elem_vals);
      }
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
      mfem::ElementTransformation &T = *V.FESpace()->GetElementTransformation(i);
      const int elem_npts = fe.GetNodes().GetNPoints();
      elem_vals.SetSize(elem_npts * vdim);
      for (int j = 0; j < elem_npts; j++)
      {
        for (int d = 0; d < vdim; d++)
        {
          // Arrange values byVDIM.
          int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                        ? d * npts + offset + j
                        : (offset + j) * dim + d;
          elem_vals(j * vdim + d) = vals(idx);
        }
      }
      offset += elem_npts;
      const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
      v.SetSize(vdofs.Size());
      fe.ProjectFromNodes(elem_vals, T, v);
      if (dof_trans)
      {
        dof_trans->TransformPrimal(v);
      }
      V.SetSubVector(vdofs, v);
    }
  }
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

}  // namespace fem

}  // namespace palace
