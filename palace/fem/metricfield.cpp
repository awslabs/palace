// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "metricfield.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::metric
{

namespace
{

// Compute the volume of a single element via numerical quadrature.
double ElementVolume(mfem::Mesh &mesh, int elem_idx)
{
  mfem::IsoparametricTransformation T;
  mesh.GetElementTransformation(elem_idx, &T);
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), T.OrderJ());
  double vol = 0.0;
  for (int j = 0; j < ir.GetNPoints(); j++)
  {
    const mfem::IntegrationPoint &ip = ir.IntPoint(j);
    T.SetIntPoint(&ip);
    vol += ip.weight * T.Weight();
  }
  return vol;
}

// Compute the characteristic element size from its volume.
// For dimension d: h = vol^{1/d}.
double CharacteristicSize(double volume, int dim)
{
  return std::pow(std::abs(volume), 1.0 / dim);
}

}  // namespace

mfem::Vector ComputeIsotropicMetric(mfem::ParMesh &mesh, const ErrorIndicator &indicator,
                                    int order, double h_min, double h_max)
{
  MPI_Comm comm = mesh.GetComm();
  const int dim = mesh.Dimension();

  // Gather the error indicator vector from all ranks to rank 0. The local indicator is
  // sized by the number of local elements on each process.
  const auto &local_err = indicator.Local();
  const int local_ne = local_err.Size();

  // Collect element counts from all ranks.
  int num_procs = Mpi::Size(comm);
  std::vector<int> recv_counts(num_procs);
  MPI_Gather(&local_ne, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, comm);

  std::vector<int> displs(num_procs, 0);
  int global_ne = 0;
  if (Mpi::Root(comm))
  {
    for (int i = 0; i < num_procs; i++)
    {
      displs[i] = global_ne;
      global_ne += recv_counts[i];
    }
  }

  // Gather all error indicators to rank 0. The error indicator vector lives on the host
  // for this operation.
  std::vector<double> global_err(global_ne);
  {
    const double *local_data = local_err.HostRead();
    MPI_Gatherv(local_data, local_ne, MPI_DOUBLE, global_err.data(), recv_counts.data(),
                displs.data(), MPI_DOUBLE, 0, comm);
  }

  // Get the serial mesh on rank 0 for vertex/element connectivity.
  // GetSerialMesh(0) gathers the full mesh to rank 0; other ranks get an empty mesh.
  // The element ordering in the serial mesh matches the Gatherv ordering above because
  // GetSerialMesh groups elements by rank in rank order, which is the same order as our
  // MPI_Gatherv.
  mfem::Mesh smesh = mesh.GetSerialMesh(0);

  // Only rank 0 computes the metric; other ranks return an empty vector.
  if (!Mpi::Root(comm))
  {
    return mfem::Vector();
  }

  const int nv = smesh.GetNV();
  const int ne = smesh.GetNE();
  MFEM_VERIFY(ne == global_ne,
              "Mismatch between gathered error indicators (" << global_ne
                                                             << ") and serial mesh elements ("
                                                             << ne << ")!");

  // Step 1: Average element errors to vertices, weighted by element volume.
  // Also compute per-element volumes and per-element characteristic sizes.
  std::vector<double> elem_vol(ne);
  std::vector<double> elem_h(ne);
  std::vector<double> vertex_err(nv, 0.0);
  std::vector<double> vertex_weight(nv, 0.0);

  for (int i = 0; i < ne; i++)
  {
    elem_vol[i] = ElementVolume(smesh, i);
    elem_h[i] = CharacteristicSize(elem_vol[i], dim);

    mfem::Array<int> verts;
    smesh.GetElementVertices(i, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      vertex_err[verts[j]] += global_err[i] * elem_vol[i];
      vertex_weight[verts[j]] += elem_vol[i];
    }
  }

  for (int i = 0; i < nv; i++)
  {
    if (vertex_weight[i] > 0.0)
    {
      vertex_err[i] /= vertex_weight[i];
    }
  }

  // Compute average mesh size and average vertex error for the equidistribution formula.
  double h_avg = 0.0;
  double total_vol = 0.0;
  for (int i = 0; i < ne; i++)
  {
    h_avg += elem_h[i] * elem_vol[i];
    total_vol += elem_vol[i];
  }
  h_avg /= total_vol;

  double eta_avg = 0.0;
  double total_weight = 0.0;
  for (int i = 0; i < nv; i++)
  {
    eta_avg += vertex_err[i] * vertex_weight[i];
    total_weight += vertex_weight[i];
  }
  eta_avg /= total_weight;

  // Avoid division by zero when all errors are identical.
  if (eta_avg <= 0.0)
  {
    eta_avg = 1.0;
  }

  // Step 2: Compute desired element size per vertex using equidistribution.
  //   h_k = h_avg * (eta_k / eta_avg)^{-2/(2p+d)}
  const double exponent = -2.0 / (2.0 * order + dim);

  // If h_min/h_max are not set, use reasonable defaults.
  if (h_min <= 0.0)
  {
    h_min = h_avg * 1.0e-3;
  }
  if (h_max <= 0.0)
  {
    h_max = h_avg * 1.0e3;
  }

  // Step 3: Build per-vertex isotropic metric tensor.
  // Upper-triangle storage: 3 components for 2D (M11, M12, M22),
  //                         6 components for 3D (M11, M12, M13, M22, M23, M33).
  const int ncomp = (dim == 2) ? 3 : 6;
  mfem::Vector metric(nv * ncomp);
  metric = 0.0;

  for (int i = 0; i < nv; i++)
  {
    double ratio = vertex_err[i] / eta_avg;
    double h_k = h_avg * std::pow(ratio, exponent);

    // Clamp to bounds.
    h_k = std::max(h_min, std::min(h_max, h_k));

    double m = 1.0 / (h_k * h_k);

    if (dim == 2)
    {
      // M11, M12, M22
      metric[i * 3 + 0] = m;
      metric[i * 3 + 1] = 0.0;
      metric[i * 3 + 2] = m;
    }
    else
    {
      // M11, M12, M13, M22, M23, M33
      metric[i * 6 + 0] = m;
      metric[i * 6 + 1] = 0.0;
      metric[i * 6 + 2] = 0.0;
      metric[i * 6 + 3] = m;
      metric[i * 6 + 4] = 0.0;
      metric[i * 6 + 5] = m;
    }
  }

  return metric;
}

void WriteSolb(const std::string &filename, const mfem::Vector &metric, int num_vertices,
               int dim)
{
  const int ncomp = (dim == 2) ? 3 : 6;
  MFEM_VERIFY(metric.Size() == num_vertices * ncomp,
              "Metric vector size (" << metric.Size() << ") does not match expected ("
                                     << num_vertices * ncomp << ")!");

  std::ofstream fo(filename);
  MFEM_VERIFY(fo.is_open(), "Unable to open metric file \"" << filename << "\"!");

  fo << "MeshVersionFormatted 2\n";
  fo << "Dimension " << dim << "\n\n";
  // SolAtVertices: count, num_fields, field_type (3 = symmetric matrix).
  fo << "SolAtVertices\n";
  fo << num_vertices << " 1 3\n";

  fo << std::scientific;
  fo.precision(15);

  const double *data = metric.HostRead();
  for (int i = 0; i < num_vertices; i++)
  {
    for (int j = 0; j < ncomp; j++)
    {
      if (j > 0)
      {
        fo << ' ';
      }
      fo << data[i * ncomp + j];
    }
    fo << '\n';
  }

  fo << "\nEnd\n";
  fo.close();
}

void WriteMeshMeshb(const std::string &filename, mfem::ParMesh &mesh)
{
  MPI_Comm comm = mesh.GetComm();

  // Gather the full mesh to rank 0.
  mfem::Mesh smesh = mesh.GetSerialMesh(0);

  // Only rank 0 writes.
  if (!Mpi::Root(comm))
  {
    return;
  }

  const int sdim = smesh.SpaceDimension();
  const int nv = smesh.GetNV();

  std::ofstream fo(filename);
  MFEM_VERIFY(fo.is_open(), "Unable to open mesh file \"" << filename << "\"!");

  fo << "MeshVersionFormatted 2\n";
  fo << "Dimension " << sdim << "\n\n";

  // Write vertices.
  fo << "Vertices\n";
  fo << nv << '\n';
  fo << std::scientific;
  fo.precision(15);
  for (int i = 0; i < nv; i++)
  {
    const double *coord = smesh.GetVertex(i);
    for (int j = 0; j < sdim; j++)
    {
      fo << coord[j] << ' ';
    }
    // Reference ID (0 for all vertices).
    fo << "0\n";
  }

  // Collect elements by type.
  std::vector<std::vector<int>> triangles, tetrahedra, quadrilaterals, hexahedra;

  for (int i = 0; i < smesh.GetNE(); i++)
  {
    mfem::Array<int> verts;
    smesh.GetElementVertices(i, verts);
    int attr = smesh.GetAttribute(i);
    auto geom = smesh.GetElementGeometry(i);

    if (geom == mfem::Geometry::TRIANGLE)
    {
      std::vector<int> elem;
      for (int k = 0; k < verts.Size(); k++)
      {
        elem.push_back(verts[k] + 1);  // 1-based indexing
      }
      elem.push_back(attr);
      triangles.push_back(std::move(elem));
    }
    else if (geom == mfem::Geometry::TETRAHEDRON)
    {
      std::vector<int> elem;
      for (int k = 0; k < verts.Size(); k++)
      {
        elem.push_back(verts[k] + 1);
      }
      elem.push_back(attr);
      tetrahedra.push_back(std::move(elem));
    }
    else if (geom == mfem::Geometry::SQUARE)
    {
      std::vector<int> elem;
      for (int k = 0; k < verts.Size(); k++)
      {
        elem.push_back(verts[k] + 1);
      }
      elem.push_back(attr);
      quadrilaterals.push_back(std::move(elem));
    }
    else if (geom == mfem::Geometry::CUBE)
    {
      std::vector<int> elem;
      for (int k = 0; k < verts.Size(); k++)
      {
        elem.push_back(verts[k] + 1);
      }
      elem.push_back(attr);
      hexahedra.push_back(std::move(elem));
    }
  }

  // Write element sections.
  auto WriteSection = [&fo](const std::string &name,
                            const std::vector<std::vector<int>> &elems)
  {
    if (elems.empty())
    {
      return;
    }
    fo << '\n' << name << '\n';
    fo << elems.size() << '\n';
    for (const auto &elem : elems)
    {
      for (std::size_t k = 0; k < elem.size(); k++)
      {
        if (k > 0)
        {
          fo << ' ';
        }
        fo << elem[k];
      }
      fo << '\n';
    }
  };

  WriteSection("Triangles", triangles);
  WriteSection("Quadrilaterals", quadrilaterals);
  WriteSection("Tetrahedra", tetrahedra);
  WriteSection("Hexahedra", hexahedra);

  // Write boundary elements.
  std::vector<std::vector<int>> bdr_triangles, bdr_edges, bdr_quads;
  for (int i = 0; i < smesh.GetNBE(); i++)
  {
    mfem::Array<int> verts;
    smesh.GetBdrElementVertices(i, verts);
    int attr = smesh.GetBdrAttribute(i);
    auto geom = smesh.GetBdrElementGeometry(i);

    std::vector<int> elem;
    for (int k = 0; k < verts.Size(); k++)
    {
      elem.push_back(verts[k] + 1);
    }
    elem.push_back(attr);

    if (geom == mfem::Geometry::SEGMENT)
    {
      bdr_edges.push_back(std::move(elem));
    }
    else if (geom == mfem::Geometry::TRIANGLE)
    {
      bdr_triangles.push_back(std::move(elem));
    }
    else if (geom == mfem::Geometry::SQUARE)
    {
      bdr_quads.push_back(std::move(elem));
    }
  }

  WriteSection("Edges", bdr_edges);
  if (!bdr_triangles.empty())
  {
    // Boundary triangles are distinct from volume triangles in .meshb format; write them
    // under the same "Triangles" keyword only if there are no volume triangles (2D mesh).
    // Otherwise they appear as boundary faces sharing the keyword. For simplicity in this
    // initial implementation, only write boundary faces if they are the sole face type.
    if (triangles.empty())
    {
      WriteSection("Triangles", bdr_triangles);
    }
  }
  if (!bdr_quads.empty() && quadrilaterals.empty())
  {
    WriteSection("Quadrilaterals", bdr_quads);
  }

  fo << "\nEnd\n";
  fo.close();
}

}  // namespace palace::metric
