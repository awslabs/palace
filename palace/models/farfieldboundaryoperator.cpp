// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "farfieldboundaryoperator.hpp"

#include <set>
#include "linalg/densematrix.hpp"
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
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!iodata.boundaries.farfield.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : iodata.boundaries.farfield.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "Absorbing boundary attribute tags must be non-negative and correspond
      //             " "to attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr - 1],
      //             "Unknown absorbing boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
      if (!bdr_warn_list.empty())
      {
        Mpi::Print("\n");
        Mpi::Warning(
            "Unknown absorbing boundary attributes!\nSolver will just ignore them!");
        utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
        Mpi::Print("\n");
      }
    }
  }

  // Set the order of the farfield boundary condition.
  order = iodata.boundaries.farfield.order;

  // Mark selected boundary attributes from the mesh as farfield.
  mfem::Array<int> farfield_bcs;
  farfield_bcs.Reserve(static_cast<int>(iodata.boundaries.farfield.attributes.size()));
  for (auto attr : iodata.boundaries.farfield.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;  // Can just ignore if wrong
    }
    farfield_bcs.Append(attr);
  }
  MFEM_VERIFY(farfield_bcs.Size() == 0 || order < 2 ||
                  iodata.problem.type == ProblemType::DRIVEN ||
                  iodata.problem.type == ProblemType::EIGENMODE,
              "Second-order farfield boundaries are only available for frequency "
              "domain simulations!");
  return farfield_bcs;
}

void FarfieldBoundaryOperator::AddDampingBdrCoefficients(double coeff,
                                                         MaterialPropertyCoefficient &fb)
{
  // First-order absorbing boundary condition.
  if (farfield_attr.Size())
  {
    MaterialPropertyCoefficient invz0_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvImpedance());
    invz0_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(farfield_attr));
    fb.AddCoefficient(invz0_func.GetAttributeToMaterial(),
                      invz0_func.GetMaterialProperties(), coeff);
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi)
{
  // Contribution for second-order absorbing BC. See Jin Section 9.3 for reference. The β
  // coefficient for the second-order ABC is 1/(2ik+2/r). Taking the radius of curvature
  // as infinity (plane wave scattering), the r-dependence vanishes and the contribution
  // is purely imaginary. Multiplying through by μ⁻¹ we get the material coefficient to ω
  // as 1 / (μ √(με)). Also, this implementation ignores the divergence term ∇⋅Eₜ, as
  // COMSOL does as well.
  if (farfield_attr.Size() && order > 1)
  {
    mfem::DenseTensor muinvc0 =
        linalg::Mult(mat_op.GetInvPermeability(), mat_op.GetLightSpeed());
    MaterialPropertyCoefficient muinvc0_func(mat_op.GetBdrAttributeToMaterial(), muinvc0);
    muinvc0_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(farfield_attr));

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
