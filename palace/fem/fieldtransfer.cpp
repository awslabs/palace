// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fieldtransfer.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace palace
{

VectorFieldTransfer::VectorFieldTransfer(mfem::ParMesh &mesh, int order)
  : fec(order, mesh.Dimension(), mfem::BasisType::GaussLegendre),
    fespace(&mesh, &fec, mesh.SpaceDimension(), mfem::Ordering::byNODES),
    mesh_sequence(mesh.GetSequence())
{
}

void VectorFieldTransfer::AddField(mfem::ParFiniteElementSpace &field_fespace,
                                   const Vector &tv)
{
  MFEM_VERIFY(field_fespace.GetParMesh() == fespace.GetParMesh(),
              "Invalid space for VectorFieldTransfer::AddField!");
  mfem::ParGridFunction gf(&field_fespace);
  gf.SetFromTrueDofs(tv);
  mfem::VectorGridFunctionCoefficient coeff(&gf);
  MFEM_VERIFY(coeff.GetVDim() == fespace.GetVDim(),
              "Invalid field dimension for VectorFieldTransfer::AddField!");
  auto &field = fields.emplace_back(std::make_unique<mfem::ParGridFunction>(&fespace));
  field->ProjectCoefficient(coeff);
}

void VectorFieldTransfer::Update()
{
  fespace.Update();
  for (auto &field : fields)
  {
    field->Update();
  }
  fespace.UpdatesFinished();
  mesh_sequence = fespace.GetParMesh()->GetSequence();
}

void VectorFieldTransfer::Remark(const std::vector<int> &old_vertices, int marking)
{
  markings = marking;
  // The marking changes nothing but the vertex order of elements: the dofs of the L2 space
  // keep their places, and only the values of an element move between its nodes. The nodes
  // of an L2 simplex are symmetric under a reordering of its vertices, so the move is a
  // permutation of the element's nodes, the same for every element with the same
  // reordering.
  fespace.Update(false);
  for (auto &field : fields)
  {
    field->Update();
  }
  fespace.UpdatesFinished();
  const mfem::ParMesh &mesh = *fespace.GetParMesh();
  mesh_sequence = mesh.GetSequence();
  const int dim = mesh.Dimension(), nv = dim + 1, vdim = fespace.GetVDim();
  MFEM_VERIFY(old_vertices.size() == static_cast<std::size_t>(nv * mesh.GetNE()),
              "Invalid vertex lists for VectorFieldTransfer::Remark!");
  auto Barycentric = [dim](const mfem::IntegrationPoint &ip, double *l)
  {
    l[1] = ip.x;
    l[2] = ip.y;
    l[3] = (dim == 3) ? ip.z : 0.0;
    l[0] = 1.0 - l[1] - l[2] - l[3];
  };
  std::map<std::vector<int>, std::vector<int>> perms;  // By the position of each new vertex
  mfem::Array<int> vdofs;
  std::vector<double> values;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int *vn = mesh.GetElement(e)->GetVertices(), *vo = &old_vertices[nv * e];
    std::vector<int> pos(nv);  // vn[a] = vo[pos[a]]
    for (int a = 0; a < nv; a++)
    {
      pos[a] = static_cast<int>(std::find(vo, vo + nv, vn[a]) - vo);
      MFEM_VERIFY(pos[a] < nv, "Element vertices changed in VectorFieldTransfer::Remark!");
    }
    if (std::is_sorted(pos.begin(), pos.end()))
    {
      continue;  // Unchanged order
    }
    const auto &nodes = fespace.GetFE(e)->GetNodes();
    const int n = nodes.GetNPoints();
    auto it = perms.find(pos);
    if (it == perms.end())
    {
      // New node i sits where the old node perm[i] did.
      std::vector<int> perm(n, -1);
      for (int i = 0; i < n; i++)
      {
        double ln[4], lo[4], lj[4];
        Barycentric(nodes.IntPoint(i), ln);
        for (int a = 0; a < nv; a++)
        {
          lo[pos[a]] = ln[a];
        }
        for (int j = 0; j < n && perm[i] < 0; j++)
        {
          Barycentric(nodes.IntPoint(j), lj);
          bool same = true;
          for (int a = 0; a < nv; a++)
          {
            same = same && std::abs(lj[a] - lo[a]) < 1.0e-12;
          }
          perm[i] = same ? j : -1;
        }
        MFEM_VERIFY(perm[i] >= 0, "Nodes not symmetric in VectorFieldTransfer::Remark!");
      }
      it = perms.emplace(pos, std::move(perm)).first;
    }
    const auto &perm = it->second;
    fespace.GetElementVDofs(e, vdofs);
    for (auto &field : fields)
    {
      double *data = field->HostReadWrite();
      values.resize(vdofs.Size());
      for (int k = 0; k < vdofs.Size(); k++)
      {
        values[k] = data[vdofs[k]];
      }
      for (int c = 0; c < vdim; c++)
      {
        for (int i = 0; i < n; i++)
        {
          data[vdofs[c * n + i]] = values[c * n + perm[i]];
        }
      }
    }
  }
}

bool VectorFieldTransfer::IsCurrent(const mfem::ParMesh &mesh) const
{
  return fespace.GetParMesh() == &mesh && mesh_sequence == mesh.GetSequence();
}

void VectorFieldTransfer::GetField(std::size_t i,
                                   mfem::ParFiniteElementSpace &field_fespace,
                                   Vector &tv) const
{
  MFEM_VERIFY(i < fields.size() && field_fespace.GetParMesh() == fespace.GetParMesh(),
              "Invalid field or space for VectorFieldTransfer::GetField!");
  mfem::ParGridFunction gf(&field_fespace);
  mfem::VectorGridFunctionCoefficient coeff(fields[i].get());
  gf.ProjectCoefficient(coeff);
  tv.SetSize(field_fespace.GetTrueVSize());
  gf.GetTrueDofs(tv);
}

}  // namespace palace
