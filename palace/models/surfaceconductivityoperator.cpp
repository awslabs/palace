// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfaceconductivityoperator.hpp"

#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

SurfaceConductivityOperator::SurfaceConductivityOperator(const IoData &iodata,
                                                         const mfem::ParMesh &mesh)
{
  // Set up finite conductivity boundary conditions.
  SetUpBoundaryProperties(iodata, mesh);
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
  bdr_sigma.SetSize(bdr_attr_max);
  bdr_mu.SetSize(bdr_attr_max);
  bdr_h.SetSize(bdr_attr_max);
  bdr_sigma = 0.0;
  bdr_mu = 0.0;
  bdr_h = 0.0;
  for (const auto &data : iodata.boundaries.conductivity)
  {
    MFEM_VERIFY(data.sigma > 0.0 && data.mu_r > 0.0,
                "Conductivity boundary has no conductivity or no "
                "permeability defined!");
    MFEM_VERIFY(data.h >= 0.0, "Conductivity boundary should have non-negative thickness!");
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          bdr_sigma(attr - 1) == 0.0 && bdr_mu(attr - 1) == 0.0 && bdr_h(attr - 1) == 0.0,
          "Multiple definitions of conductivity boundary properties for boundary attribute "
              << attr << "!");
      bdr_sigma(attr - 1) = data.sigma;
      bdr_mu(attr - 1) = data.mu_r;
      bdr_h(attr - 1) = data.h;
      if (data.external)
      {
        // External surfaces have twice the effective thickness since the BC is applied at
        // one side.
        bdr_h(attr - 1) *= 2.0;
      }
    }
  }

  // Mark selected boundary attributes from the mesh as finite conductivity.
  mfem::Array<int> conductivity_bcs;
  for (const auto &data : iodata.boundaries.conductivity)
  {
    for (auto attr : data.attributes)
    {
      conductivity_bcs.Append(attr);
    }
  }
  MFEM_VERIFY(conductivity_bcs.Size() == 0 ||
                  iodata.problem.type == config::ProblemData::Type::DRIVEN,
              "Finite conductivity boundaries are only available for frequency "
              "domain driven simulations!");
  mesh::AttrToMarker(bdr_attr_max, conductivity_bcs, conductivity_marker);
}

void SurfaceConductivityOperator::PrintBoundaryInfo(const IoData &iodata,
                                                    const mfem::ParMesh &mesh)
{
  if (conductivity_marker.Size() && conductivity_marker.Max() == 0)
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin finite conductivity BC at attributes:\n");
  for (int i = 0; i < conductivity_marker.Size(); i++)
  {
    if (conductivity_marker[i])
    {
      const int attr = i + 1;
      mfem::Vector normal = mesh::GetSurfaceNormal(mesh, attr);
      Mpi::Print(" {:d}: σ = {:.3e} S/m", attr,
                 iodata.DimensionalizeValue(IoData::ValueType::CONDUCTIVITY, bdr_sigma(i)));
      if (bdr_h(i) > 0.0)
      {
        Mpi::Print(", h = {:.3e} m",
                   iodata.DimensionalizeValue(IoData::ValueType::LENGTH, bdr_h(i)));
      }
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f}, {:+.1f})", normal(0), normal(1), normal(2));
      }
      else
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f})", normal(0), normal(1));
      }
      Mpi::Print("\n");
    }
  }
}

void SurfaceConductivityOperator::AddExtraSystemBdrCoefficients(double omega,
                                                                SumMatrixCoefficient &fbr,
                                                                SumMatrixCoefficient &fbi)
{
  if (conductivity_marker.Size() && conductivity_marker.Max() > 0)
  {
    // If the provided conductor thickness is empty (zero), prescribe a surface impedance
    // (1+i)/σδ, where δ is the skin depth. If it is nonzero, use a finite thickness
    // modification which correctly produces the DC limit when h << δ. See the Ansys HFSS
    // user manual section titled "Surface Impedance Boundary Condition for Metal Traces of
    // Finite Thickness."
    mfem::Vector vr(bdr_sigma.Size()), vi(bdr_sigma.Size());
    for (int i = 0; i < bdr_sigma.Size(); i++)
    {
      if (bdr_sigma(i) > 0.0)
      {
        double delta = std::sqrt(2.0 / (bdr_mu(i) * bdr_sigma(i) * omega));
        std::complex<double> Z = 1.0 / (bdr_sigma(i) * delta);
        Z.imag(Z.real());
        if (bdr_h(i) > 0.0)
        {
          double nu = bdr_h(i) / delta;
          double den = std::cosh(nu) - std::cos(nu);
          Z.real(Z.real() * (std::sinh(nu) + std::sin(nu)) / den);
          Z.imag(Z.imag() * (std::sinh(nu) - std::sin(nu)) / den);
        }
        // The BC term has coefficient iω/Z (like for standard lumped surface impedance).
        std::complex<double> s(1i * omega / Z);
        vr(i) = s.real();
        vi(i) = s.imag();
      }
      else
      {
        vr(i) = vi(i) = 0.0;  // Not a conductivity boundary
      }
    }
    fbr.AddCoefficient(std::make_unique<mfem::PWConstCoefficient>(vr), conductivity_marker);
    fbi.AddCoefficient(std::make_unique<mfem::PWConstCoefficient>(vi), conductivity_marker);
  }
}

}  // namespace palace
