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
// Vector fields of a finite element space (for example, H(curl)) kept through refinement
// and rebalancing of their mesh, for interpolation into a space on the modified mesh. The
// fields are stored in a discontinuous vector L2 space of the same order, which represents
// them exactly (for affine elements) and has no oriented degrees of freedom: MFEM's
// transfer of H(curl) and H(div) grid functions through ParMesh::Rebalance does not
// account for the signs of their degrees of freedom.
//
class VectorFieldTransfer
{
private:
  mfem::L2_FECollection fec;
  mfem::ParFiniteElementSpace fespace;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> fields;

  // Sequence number of the mesh followed by the fields.
  long mesh_sequence;

public:
  VectorFieldTransfer(mfem::ParMesh &mesh, int order);

  // Add a field given by its true dofs in a space on the same mesh.
  void AddField(mfem::ParFiniteElementSpace &field_fespace, const Vector &tv);

  // Follow a modification of the mesh (refinement or rebalancing), to be called after each
  // one.
  void Update();

  // Mesh of the fields (not owned).
  const mfem::ParMesh *GetParMesh() const { return fespace.GetParMesh(); }

  // Whether the fields follow the current state of the given mesh object.
  bool IsCurrent(const mfem::ParMesh &mesh) const;

  // Number of stored fields.
  auto Size() const { return fields.size(); }

  // Interpolate the i-th field into a space on the same mesh, as true dofs.
  void GetField(std::size_t i, mfem::ParFiniteElementSpace &field_fespace,
                Vector &tv) const;
};

}  // namespace palace

#endif  // PALACE_FEM_FIELD_TRANSFER_HPP
