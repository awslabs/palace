// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "linalg/operator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

void FullWaveSingularFeatures::Preprocess(const IoData &iodata,
                                          const std::unique_ptr<mfem::Mesh> &serial_mesh,
                                          MPI_Comm comm)
{
  serial_sheet_features = {};
  local_sheet_features = {};
  serial_line_features = {};
  local_line_features = {};
  source_vertex_ids.clear();
  source_element_ids.clear();
  dimension = 0;
  if (!iodata.solver.singular_elements.Enabled())
  {
    return;
  }

  if (Mpi::Root(comm))
  {
    MFEM_VERIFY(serial_mesh,
                "Root rank has no serial mesh for full-wave singular feature extraction!");
    dimension = serial_mesh->Dimension();
    if (dimension == 3)
    {
      serial_sheet_features = fem::singular::ExtractSerialSheetFeatures(
          *serial_mesh, iodata.solver.singular_elements.attributes);
    }
    else
    {
      MFEM_VERIFY(dimension == 2 && serial_mesh->SpaceDimension() == 2,
                  "Full-wave singular enrichment requires a 2D triangular or 3D "
                  "tetrahedral mesh!");
      MFEM_VERIFY(iodata.solver.singular_elements.order == 1,
                  "Triangular full-wave singular enrichment currently supports only "
                  "Solver.SingularElements.Order = 1!");
      serial_line_features = fem::singular::ExtractSerialLineTipFeatures(
          *serial_mesh, iodata.solver.singular_elements.attributes);
    }
  }
  Mpi::Broadcast(1, &dimension, 0, comm);
  if (dimension == 2)
  {
    fem::singular::BroadcastSerialLineTipFeatures(serial_line_features, comm);
    MFEM_VERIFY(!serial_line_features.Empty(),
                "Full-wave singular extraction produced no internal PEC line tips!");
  }
  else
  {
    MFEM_VERIFY(dimension == 3, "Full-wave singular enrichment requires a 2D or 3D mesh!");
    fem::singular::BroadcastSerialSheetFeatures(serial_sheet_features, comm);
    MFEM_VERIFY(!serial_sheet_features.Empty(),
                "Full-wave singular extraction produced no PEC sheet-edge features!");
  }
}

void FullWaveSingularFeatures::ProcessPartitionedMesh(
    const IoData &iodata, const mfem::ParMesh &parallel_mesh,
    const mesh::PartitionMetadata &metadata)
{
  MFEM_VERIFY(iodata.solver.singular_elements.Enabled() &&
                  parallel_mesh.Dimension() == dimension,
              "Unexpected or inconsistent full-wave singular partition metadata!");
  if (dimension == 2)
  {
    local_line_features = fem::singular::DistributeSerialLineTipFeatures(
        serial_line_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  else
  {
    local_sheet_features = fem::singular::DistributeSerialSheetFeatures(
        serial_sheet_features, parallel_mesh, metadata.source_vertex_ids,
        metadata.source_element_ids);
  }
  source_vertex_ids = metadata.source_vertex_ids;
  source_element_ids = metadata.source_element_ids;
  MFEM_VERIFY(source_vertex_ids.size() == static_cast<std::size_t>(parallel_mesh.GetNV()) &&
                  source_element_ids.size() ==
                      static_cast<std::size_t>(parallel_mesh.GetNE()),
              "Full-wave singular source-entity metadata is incomplete!");
}

const fem::singular::FeatureTopology *FullWaveSingularFeatures::GetSheetFeatures() const
{
  return dimension == 3 && !local_sheet_features.Empty() ? &local_sheet_features : nullptr;
}

const fem::singular::TriangleFeatureTopology *
FullWaveSingularFeatures::GetLineFeatures() const
{
  return dimension == 2 && !local_line_features.Empty() ? &local_line_features : nullptr;
}

const std::vector<fem::singular::GlobalVertexId> *
FullWaveSingularFeatures::GetSourceVertexIds() const
{
  return source_vertex_ids.empty() ? nullptr : &source_vertex_ids;
}

double SingularComplexQuadraticForm(MPI_Comm comm, const mfem::Operator &op,
                                    const ComplexVector &x)
{
  MFEM_VERIFY(op.Height() == op.Width() && op.Width() == x.Size(),
              "Invalid dimensions for a singular full-wave quadratic form!");
  Vector work(op.Height());
  work.UseDevice(true);
  op.Mult(x.Real(), work);
  double value = linalg::Dot(comm, x.Real(), work);
  op.Mult(x.Imag(), work);
  value += linalg::Dot(comm, x.Imag(), work);
  MFEM_VERIFY(std::isfinite(value),
              "Singular full-wave quadratic form produced a nonfinite value!");
  return value;
}

SingularFullWaveEnergy MeasureSingularFullWaveEnergy(MPI_Comm comm,
                                                     const mfem::Operator &mass,
                                                     const mfem::Operator &stiffness,
                                                     const ComplexVector &electric_field,
                                                     std::complex<double> omega)
{
  const double omega_squared = std::norm(omega);
  MFEM_VERIFY(omega_squared > 0.0 && std::isfinite(omega_squared),
              "Singular full-wave energy requires a finite nonzero frequency!");
  const double electric_form = SingularComplexQuadraticForm(comm, mass, electric_field);
  const double magnetic_form =
      SingularComplexQuadraticForm(comm, stiffness, electric_field) / omega_squared;
  const double scale = std::max({1.0, std::abs(electric_form), std::abs(magnetic_form)});
  constexpr double positivity_tolerance = 256.0 * std::numeric_limits<double>::epsilon();
  MFEM_VERIFY(electric_form >= -positivity_tolerance * scale &&
                  magnetic_form >= -positivity_tolerance * scale,
              "Singular full-wave energy is negative beyond floating-point roundoff!");
  return {0.5 * std::max(0.0, electric_form), 0.5 * std::max(0.0, magnetic_form)};
}

}  // namespace palace
