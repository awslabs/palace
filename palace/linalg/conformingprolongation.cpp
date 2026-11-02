// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "conformingprolongation.hpp"

namespace palace
{

ConformingProlongation::ConformingProlongation(
    const mfem::ConformingProlongationOperator &P)
  : Operator(P.Height(), P.Width()), P(P),
    gc2(P.GetGroupCommunicator().GetGroupTopology(), mfem::GroupCommunicator::byNeighbor)
{
  const auto &gc = P.GetGroupCommunicator();
  const mfem::Table &group_ldof = gc.GroupLDofTable();
  const int ngroups = group_ldof.Size();
  mfem::Table &group_ldof2 = gc2.GroupLDofTable();
  group_ldof2.MakeI(ngroups);
  for (int g = 0; g < ngroups; g++)
  {
    group_ldof2.AddColumnsInRow(g, 2 * group_ldof.RowSize(g));
  }
  group_ldof2.MakeJ();
  for (int g = 0; g < ngroups; g++)
  {
    const int *row = group_ldof.GetRow(g);
    for (int k = 0; k < group_ldof.RowSize(g); k++)
    {
      group_ldof2.AddConnection(g, 2 * row[k]);
      group_ldof2.AddConnection(g, 2 * row[k] + 1);
    }
  }
  group_ldof2.ShiftUpI();
  gc2.Finalize();

  // As in the mfem::ConformingProlongationOperator constructor.
  for (int g = 1; g < ngroups; g++)
  {
    if (!gc.GetGroupTopology().IAmMaster(g))
    {
      external_ldofs.Append(group_ldof.GetRow(g), group_ldof.RowSize(g));
    }
  }
  external_ldofs.Sort();
  MFEM_VERIFY(external_ldofs.Size() == Height() - Width(),
              "Unexpected number of external L-DOFs for conforming prolongation!");
}

void ConformingProlongation::Mult(const Vector &xr, const Vector &xi, Vector &yr,
                                  Vector &yi) const
{
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostWrite(), *YI = yi.HostWrite();
  const int m = external_ldofs.Size(), n = Height();
  buf.resize(2 * static_cast<std::size_t>(n));
  double *B = buf.data();

  // Owned L-DOFs from the true DOFs (true DOF j - i for L-DOF j after i external ones).
  int j = 0;
  for (int i = 0; i <= m; i++)
  {
    const int end = (i < m) ? external_ldofs[i] : n;
    for (; j < end; j++)
    {
      B[2 * j] = XR[j - i];
      B[2 * j + 1] = XI[j - i];
    }
    j = end + 1;
  }
  gc2.BcastBegin(B, 0);
  gc2.BcastEnd(B, 0);
  for (int k = 0; k < n; k++)
  {
    YR[k] = B[2 * k];
    YI[k] = B[2 * k + 1];
  }
}

void ConformingProlongation::MultTranspose(const Vector &xr, const Vector &xi, Vector &yr,
                                           Vector &yi) const
{
  const double *XR = xr.HostRead(), *XI = xi.HostRead();
  double *YR = yr.HostWrite(), *YI = yi.HostWrite();
  const int m = external_ldofs.Size(), n = Height();
  buf.resize(2 * static_cast<std::size_t>(n));
  double *B = buf.data();
  for (int k = 0; k < n; k++)
  {
    B[2 * k] = XR[k];
    B[2 * k + 1] = XI[k];
  }
  gc2.ReduceBegin(B);
  gc2.ReduceEnd(B, 0, mfem::GroupCommunicator::Sum<double>);
  int j = 0;
  for (int i = 0; i <= m; i++)
  {
    const int end = (i < m) ? external_ldofs[i] : n;
    for (; j < end; j++)
    {
      YR[j - i] = B[2 * j];
      YI[j - i] = B[2 * j + 1];
    }
    j = end + 1;
  }
}

}  // namespace palace
