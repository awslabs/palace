// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FIELD_TRANSFER_HPP
#define PALACE_FEM_FIELD_TRANSFER_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

//
// Vector fields of a finite element space (for example, H(curl)) kept through the
// refinement and rebalancing of their mesh. They are stored in a discontinuous vector L2
// space of the same order, which represents them exactly on affine elements and has no
// oriented degrees of freedom: MFEM's transfer of H(curl) and H(div) grid functions
// through ParMesh::Rebalance mishandles their signs (mfem/mfem#5125).
//
class VectorFieldTransfer
{
private:
  mfem::L2_FECollection fec;
  mfem::ParFiniteElementSpace fespace;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> fields;

  // Number of the last marking the fields were put back in place for (Remark).
  int markings = 0;

  // Sequence number of the mesh followed by the fields.
  long mesh_sequence;

public:
  VectorFieldTransfer(mfem::ParMesh &mesh, int order);

  // Add a field given by its true dofs in a space on the same mesh.
  void AddField(mfem::ParFiniteElementSpace &field_fespace, const Vector &tv);

  // Follow a modification of the mesh (refinement or rebalancing), to be called after each
  // one.
  void Update();

  // Follow a new marking of a conforming simplex mesh, which reordered the vertices of its
  // elements: old_vertices holds the vertices of each element before it (dim + 1 each), and
  // marking counts the markings, which IsCurrent checks, so that a missed one is not
  // silent.
  void Remark(const std::vector<int> &old_vertices, int marking = 0);

  // Mesh of the fields (not owned).
  const mfem::ParMesh *GetParMesh() const { return fespace.GetParMesh(); }

  // Whether the fields follow the current state of the given mesh object.
  bool IsCurrent(const mfem::ParMesh &mesh) const;
  bool IsCurrent(const mfem::ParMesh &mesh, int marking) const
  {
    return IsCurrent(mesh) && marking == markings;
  }

  // Number of stored fields.
  auto Size() const { return fields.size(); }

  // Interpolate the i-th field into a space on the same mesh, as true dofs.
  void GetField(std::size_t i, mfem::ParFiniteElementSpace &field_fespace,
                Vector &tv) const;
};

}  // namespace palace

#endif  // PALACE_FEM_FIELD_TRANSFER_HPP
