// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "coefficient.hpp"

namespace palace
{

void BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
    int i, mfem::ParMesh &mesh, const std::unordered_map<int, int> &local_to_shared,
    mfem::ElementTransformation *&T1, mfem::ElementTransformation *&T2,
    const mfem::IntegrationPoint *ip)
{
  // Return transformations for elements attached to the given boundary element. T1 always
  // exists but T2 may not if the element is truly a single-sided boundary.
  int f, o;
  int iel1, iel2, info1, info2;
  mesh.GetBdrElementFace(i, &f, &o);
  mesh.GetFaceElements(f, &iel1, &iel2);
  mesh.GetFaceInfos(f, &info1, &info2);

  // XX TODO THREAD SAFETY: THIS IS NOT THREAD SAFE, SEE MFEM PR AND UPGRADE!

  // Master faces can never be boundary elements, thus only need to check for the state of
  // info2 and el2, and do not need to access the ncface numbering. See mfem::Mesh::FaceInfo
  // for details.
  mfem::FaceElementTransformations *FET;
  if (info2 >= 0 && iel2 < 0)
  {
    // Face is shared with another subdomain.
    const int &ishared = local_to_shared.at(f);
    FET = mesh.GetSharedFaceTransformations(ishared);
  }
  else
  {
    // Face is either internal to the subdomain, or a true one-sided boundary.
    FET = mesh.GetFaceElementTransformations(f);
  }
  T1 = &FET->GetElement1Transformation();
  T2 = (info2 >= 0) ? &FET->GetElement2Transformation() : nullptr;

  // Boundary elements and boundary faces may have different orientations so adjust the
  // integration point if necessary. See mfem::GridFunction::GetValue and GetVectorValue.
  if (ip)
  {
    mfem::IntegrationPoint fip =
        mfem::Mesh::TransformBdrElementToFace(FET->GetGeometryType(), o, *ip);
    FET->SetAllIntPoints(&fip);
  }
}

void BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
    mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip,
    mfem::ElementTransformation *&T1, mfem::ElementTransformation *&T2, mfem::Vector *C1)
{
  // Get the element transformations neighboring the element, and set the integration point
  // too.
  MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
              "Unexpected element type in BdrGridFunctionCoefficient!");
  GetBdrElementNeighborTransformations(T.ElementNo, mesh, local_to_shared, T1, T2, &ip);

  // If desired, get vector pointing from center of boundary element into element 1 for
  // orientations.
  if (C1)
  {
    mfem::Vector CF(T.GetSpaceDim());
    mfem::ElementTransformation &TF = *mesh.GetFaceTransformation(T.ElementNo);
    TF.Transform(mfem::Geometries.GetCenter(mesh.GetFaceGeometry(T.ElementNo)), CF);

    C1->SetSize(T.GetSpaceDim());
    T1->Transform(mfem::Geometries.GetCenter(T1->GetGeometryType()), *C1);
    *C1 -= CF;  // Points into element 1 from the face
  }
}

}  // namespace palace
