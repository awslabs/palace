// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fieldtransfer.hpp"

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
