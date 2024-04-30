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
  const double bb_t = 0.1;  // MFEM defaults
  const double newton_tol = 1.0e-12;
  const int npts = static_cast<int>(iodata.domains.postpro.probe.size());
  MFEM_VERIFY(
      mesh.Dimension() == mesh.SpaceDimension(),
      "Probe postprocessing functionality requires mesh dimension == space dimension!");
  mfem::Vector xyz(npts * mesh.SpaceDimension());
  op_idx.resize(npts);
  int i = 0;
  for (const auto &[idx, data] : iodata.domains.postpro.probe)
  {
    // Use default ordering byNODES.
    xyz(i) = data.center[0];
    xyz(npts + i) = data.center[1];
    if (mesh.SpaceDimension() == 3)
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
  const int dim = U.VectorDim();
  std::vector<double> vals(npts * dim);
  mfem::Vector v(npts * dim);
  op.Interpolate(U, v);
  for (int d = 0; d < dim; d++)
  {
    for (int i = 0; i < npts; i++)
    {
      vals[i * dim + d] = v(d * npts + i);
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

}  // namespace palace
