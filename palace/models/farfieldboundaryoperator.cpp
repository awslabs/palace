// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "farfieldboundaryoperator.hpp"

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

FarfieldBoundaryOperator::FarfieldBoundaryOperator(const IoData &iodata,
                                                   const MaterialOperator &mat_op,
                                                   const mfem::ParMesh &mesh)
  : mat_op(mat_op), farfield_attr(SetUpBoundaryProperties(iodata, mesh))
{
  // Print out BC info for all farfield attributes.
  if (farfield_attr.Size())
  {
    Mpi::Print("\nConfiguring Robin absorbing BC (order {:d}) at attributes:\n", order);
    std::sort(farfield_attr.begin(), farfield_attr.end());
    utils::PrettyPrint(farfield_attr);
  }
}

mfem::Array<int>
FarfieldBoundaryOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                  const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  if (!iodata.boundaries.farfield.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (auto attr : iodata.boundaries.farfield.attributes)
    {
      MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                  "Absorbing boundary attribute tags must be non-negative and correspond "
                  "to attributes in the mesh!");
      MFEM_VERIFY(bdr_attr_marker[attr - 1],
                  "Unknown absorbing boundary attribute " << attr << "!");
    }
  }

  // Set the order of the farfield boundary condition.
  order = iodata.boundaries.farfield.order;

  // Mark selected boundary attributes from the mesh as farfield.
  mfem::Array<int> farfield_bcs;
  farfield_bcs.Append(iodata.boundaries.farfield.attributes.data(),
                      iodata.boundaries.farfield.attributes.size());
  MFEM_VERIFY(farfield_bcs.Size() == 0 || order < 2 ||
                  iodata.problem.type == config::ProblemData::Type::DRIVEN,
              "Second-order farfield boundaries are only available for frequency "
              "domain driven simulations!");
  return farfield_bcs;
}

void FarfieldBoundaryOperator::AddDampingBdrCoefficients(double coef,
                                                         MaterialPropertyCoefficient &fb)
{
  // First-order absorbing boundary condition.
  if (farfield_attr.Size())
  {
    MaterialPropertyCoefficient invz0_func(mat_op, mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvImpedance());
    invz0_func.RestrictCoefficient(mat_op.GetBdrAttributeGlobalToLocal(farfield_attr));
    fb.AddCoefficient(invz0_func.GetAttributeToMaterial(),
                      invz0_func.GetMaterialProperties(), coef);
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi)
{
  // Contribution for second-order absorbing BC. See Jin Section 9.3 for reference. The β
  // coefficient for the second-order ABC is 1/(2ik+2/r). Taking the radius of curvature as
  // infinity (plane wave scattering), the r-dependence vanishes and the contribution is
  // purely imaginary. Multiplying through by μ⁻¹ we get the material coefficient to ω as
  // 1 / (μ √με). Also, this implementation ignores the divergence term ∇⋅Eₜ, as COMSOL
  // does as well.
  if (farfield_attr.Size() && order > 1)
  {
    mfem::DenseTensor muinvc0(mat_op.GetLightSpeed());
    for (int k = 0; k < muinvc0.SizeK(); k++)
    {
      Mult(mat_op.GetInvPermeability()(k), mat_op.GetLightSpeed()(k), muinvc0(k));
    }
    MaterialPropertyCoefficient muinvc0_func(mat_op, mat_op.GetBdrAttributeToMaterial(),
                                             muinvc0);
    muinvc0_func.RestrictCoefficient(mat_op.GetBdrAttributeGlobalToLocal(farfield_attr));

    // Instead getting the correct normal of farfield boundary elements, just pick the
    // the first element normal. This is fine as long as the farfield material properties
    // are not anisotropic.
    mfem::Vector normal(mat_op.SpaceDimension());
    normal = 0.0;
    normal(0) = 1.0;
    muinvc0_func.NormalProjectedCoefficient(normal);

    dfbi.AddCoefficient(muinvc0_func.GetAttributeToMaterial(),
                        muinvc0_func.GetMaterialProperties(), 0.5 / omega);
  }
}

}  // namespace palace
