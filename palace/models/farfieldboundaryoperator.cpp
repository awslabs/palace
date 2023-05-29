// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "farfieldboundaryoperator.hpp"

#include "fem/coefficient.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

FarfieldBoundaryOperator::FarfieldBoundaryOperator(const IoData &iodata,
                                                   const MaterialOperator &mat,
                                                   const mfem::ParMesh &mesh)
  : mat_op(mat)
{
  // Set up impedance boundary conditions.
  SetUpBoundaryProperties(iodata, mesh);

  // Print out BC info for all farfield attributes.
  if (farfield_marker.Max() > 0)
  {
    Mpi::Print("\nConfiguring Robin absorbing BC (order {:d}) at attributes:\n", order);
    utils::PrettyPrintMarker(farfield_marker);
  }
}

void FarfieldBoundaryOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                       const mfem::ParMesh &mesh)
{
  // Check that impedance boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Max();
  if (!iodata.boundaries.farfield.empty())
  {
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
  MFEM_VERIFY(iodata.boundaries.farfield.attributes.empty() || order < 2 ||
                  iodata.problem.type == config::ProblemData::Type::DRIVEN,
              "Second-order farfield boundaries are only available for frequency "
              "domain driven simulations!");
  mesh::AttrToMarker(bdr_attr_max, iodata.boundaries.farfield.attributes, farfield_marker);
}

void FarfieldBoundaryOperator::AddDampingBdrCoefficients(double coef,
                                                         SumMatrixCoefficient &fb)
{
  // First-order absorbing boundary condition.
  if (farfield_marker.Max() > 0)
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_Z0;
    fb.AddCoefficient(
        std::make_unique<BdrMaterialPropertyCoefficient<MatType>>(mat_op, coef),
        farfield_marker);
  }
}

void FarfieldBoundaryOperator::AddExtraSystemBdrCoefficients(double omega,
                                                             SumCoefficient &dfbr,
                                                             SumCoefficient &dfbi)
{
  // Contribution for second-order absorbing BC. See Jin Section 9.3 for reference. The β
  // coefficient for the second-order ABC is 1/(2ik+2/r). Taking the radius of curvature as
  // infinity (plane wave scattering), the r-dependence vanishes and the contribution is
  // purely imaginary. Multiplying through by μ⁻¹ we get the material coefficient to ω as
  // 1 / (μ √με). Also, this implementation ignores the divergence term ∇⋅Eₜ, as COMSOL
  // does as well.
  if (farfield_marker.Max() > 0 && order > 1)
  {
    constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY_C0;
    dfbi.AddCoefficient(
        std::make_unique<NormalProjectedCoefficient>(
            std::make_unique<BdrMaterialPropertyCoefficient<MatType>>(mat_op, 0.5 / omega)),
        farfield_marker);
  }
}

}  // namespace palace
