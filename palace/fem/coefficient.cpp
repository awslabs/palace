// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "coefficient.hpp"

namespace palace
{

bool BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(
    int i, const mfem::ParMesh &mesh, mfem::FaceElementTransformations &FET,
    mfem::IsoparametricTransformation &T1, mfem::IsoparametricTransformation &T2,
    const mfem::IntegrationPoint *ip)
{
  // Return transformations for elements attached to the given boundary element. FET.Elem1
  // always exists but FET.Elem2 may not if the element is truly a single-sided boundary.
  int f, o;
  int iel1, iel2, info1, info2;
  mesh.GetBdrElementFace(i, &f, &o);
  mesh.GetFaceElements(f, &iel1, &iel2);
  mesh.GetFaceInfos(f, &info1, &info2);

  // Master faces can never be boundary elements, thus only need to check for the state of
  // info2 and el2, and do not need to access the ncface numbering. See mfem::Mesh::FaceInfo
  // for details.
  if (info2 >= 0 && iel2 < 0)
  {
    // Face is shared with another subdomain.
    mesh.GetSharedFaceTransformationsByLocalIndex(f, FET, T1, T2);
  }
  else
  {
    // Face is either internal to the subdomain, or a true one-sided boundary.
    mesh.GetFaceElementTransformations(f, FET, T1, T2);
  }

  // Boundary elements and boundary faces may have different orientations so adjust the
  // integration point if necessary. See mfem::GridFunction::GetValue and GetVectorValue.
  if (ip)
  {
    mfem::IntegrationPoint fip =
        mfem::Mesh::TransformBdrElementToFace(FET.GetGeometryType(), o, *ip);
    FET.SetAllIntPoints(&fip);
  }

  // Return whether or not the boundary element and face share the same orientations.
  return (o % 2 == 0);
}

}  // namespace palace
