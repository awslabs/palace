// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceconductivityoperator.hpp"

#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

SurfaceConductivityOperator::SurfaceConductivityOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         mfem::ParMesh &mesh)
  : mat_op(mat_op), conductivity_attr(SetUpBoundaryProperties(iodata, mesh))
{
  // Print out BC info for all finite conductivity boundary attributes.
  PrintBoundaryInfo(iodata, mesh);
}

void SurfaceConductivityOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                          const mfem::ParMesh &mesh)
{
  // Check that conductivity boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!iodata.boundaries.conductivity.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &data : iodata.boundaries.conductivity)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Conductivity boundary attribute tags must be non-negative and "
                    "correspond to attributes in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown conductivity boundary attribute " << attr << "!");
      }
    }
  }

  // Finite conductivity boundaries are defined using the user provided surface conductivity
  // and optionally conductor thickness.
  mfem::Array<bool> marker(bdr_attr_max);
  marker = false;
  conductivity_data.reserve(iodata.boundaries.conductivity.size());
  for (const auto &data : iodata.boundaries.conductivity)
  {
    MFEM_VERIFY(data.sigma > 0.0 && data.mu_r > 0.0,
                "Conductivity boundary has no conductivity or no "
                "permeability defined!");
    MFEM_VERIFY(data.h >= 0.0, "Conductivity boundary should have non-negative thickness!");
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          !marker[attr - 1],
          "Multiple definitions of conductivity boundary properties for boundary attribute "
              << attr << "!");
      marker[attr - 1] = true;
    }
    auto &bdr = conductivity_data.emplace_back();
    bdr.sigma = data.sigma;
    bdr.mu = data.mu_r;
    bdr.h = data.h;
    if (data.external)
    {
      // External surfaces have twice the effective thickness since the BC is applied at one
      // side.
      bdr.h *= 2.0;
    }
    bdr.attr.Append(data.attributes.data(), data.attributes.size());
  }

  // Mark selected boundary attributes from the mesh as finite conductivity.
  mfem::Array<int> conductivity_bcs;
  for (const auto &data : iodata.boundaries.conductivity)
  {
    conductivity_bcs.Append(data.attributes.data(), data.attributes.size());
  }
  MFEM_VERIFY(conductivity_bcs.Size() == 0 ||
                  iodata.problem.type == config::ProblemData::Type::DRIVEN,
              "Finite conductivity boundaries are only available for frequency "
              "domain driven simulations!");
  return conductivity_bcs;
}

void SurfaceConductivityOperator::PrintBoundaryInfo(const IoData &iodata,
                                                    mfem::ParMesh &mesh)
{
  if (conductivity_data.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin finite conductivity BC at attributes:\n");
  for (const auto &bdr : conductivity_data)
  {
    for (auto attr : bdr.attr)
    {
      mfem::Vector nor;
      mesh::GetSurfaceNormal(mesh, attr, nor);
      Mpi::Print(" {:d}: σ = {:.3e} S/m", attr,
                 iodata.DimensionalizeValue(IoData::ValueType::CONDUCTIVITY, bdr.sigma(i)));
      if (bdr.h(i) > 0.0)
      {
        Mpi::Print(", h = {:.3e} m",
                   iodata.DimensionalizeValue(IoData::ValueType::LENGTH, bdr.h(i)));
      }
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f}, {:+.1f})", nor(0), nor(1), nor(2));
      }
      else
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f})", nor(0), nor(1));
      }
      Mpi::Print("\n");
    }
  }
}

void SurfaceConductivityOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  // If the provided conductor thickness is empty (zero), prescribe a surface impedance
  // (1+i)/σδ, where δ is the skin depth. If it is nonzero, use a finite thickness
  // modification which correctly produces the DC limit when h << δ. See the Ansys HFSS
  // user manual section titled "Surface Impedance Boundary Condition for Metal Traces of
  // Finite Thickness."
  for (const auto &bdr : conductivity_data)
  {
    if (std::abs(bdr.sigma) > 0.0)
    {
      double delta = std::sqrt(2.0 / (bdr.mu(i) * bdr.sigma(i) * omega));
      std::complex<double> Z = 1.0 / (bdr.sigma(i) * delta);
      Z.imag(Z.real());
      if (bdr.h(i) > 0.0)
      {
        double nu = bdr.h(i) / delta;
        double den = std::cosh(nu) - std::cos(nu);
        Z.real(Z.real() * (std::sinh(nu) + std::sin(nu)) / den);
        Z.imag(Z.imag() * (std::sinh(nu) - std::sin(nu)) / den);
      }
      // The BC term has coefficient iω/Z (like for standard lumped surface impedance).
      std::complex<double> s(1i * omega / Z);
      fbr.AddMaterialProperty(mat_op.GetAttributeGlobalToLocal(bdr.attr), s.real());
      fbi.AddMaterialProperty(mat_op.GetAttributeGlobalToLocal(bdr.attr), s.imag());
    }
  }
}

}  // namespace palace
