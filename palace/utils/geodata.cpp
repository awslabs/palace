// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"

#include <array>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Floating point precision for mesh IO. This precision is important, make sure nothing is
// lost!
const auto MSH_FLT_PRECISION = std::numeric_limits<double>::max_digits10;

// Load the serial mesh from disk.
std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &);

// Optionally reorder mesh elements based on MFEM's internal reordeing tools for improved
// cache usage.
void ReorderMesh(mfem::Mesh &);

// Generate element-based mesh partitioning, using either a provided file or METIS.
std::unique_ptr<int[]> GetMeshPartitioning(mfem::Mesh &, int, const std::string &);

// Cleanup the provided serial mesh by removing unnecessary domain and elements, adding
// boundary elements for material interfaces and exterior boundaries, and adding boundary
// elements for subdomain interfaces.
std::map<int, std::array<int, 2>> CheckMesh(std::unique_ptr<mfem::Mesh> &,
                                            const std::unique_ptr<int[]> &, const IoData &,
                                            bool, bool, bool);

// Given a serial mesh on the root processor and element partitioning, create a parallel
// mesh oer the given communicator.
std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm, std::unique_ptr<mfem::Mesh> &,
                                              std::unique_ptr<int[]> &);

// Get list of domain and boundary attribute markers used in configuration file for mesh
// cleaning.
void GetUsedAttributeMarkers(const IoData &, int, int, mfem::Array<int> &,
                             mfem::Array<int> &);

}  // namespace

namespace mesh
{

std::unique_ptr<mfem::ParMesh> ReadMesh(MPI_Comm comm, const IoData &iodata, bool reorder,
                                        bool clean, bool add_bdr, bool unassembled,
                                        Timer &timer)
{
  // On root, read the serial mesh (converting format if necessary), and do all necessary
  // serial preprocessing. When finished, distribute the mesh to all processes. Count disk
  // I/O time separately for the mesh read from file.
  std::unique_ptr<mfem::Mesh> smesh;
  auto t0 = timer.Now();
  if (Mpi::Root(comm))
  {
    // Optionally reorder elements (and vertices) based on spatial location after loading
    // the serial mesh.
    smesh = LoadMesh(iodata.model.mesh);
    if (reorder)
    {
      ReorderMesh(*smesh);
    }
  }
  Mpi::Barrier(comm);
  timer.io_time += timer.Now() - t0;

  std::unique_ptr<int[]> partitioning;
  if (Mpi::Root(comm))
  {
    // Generate the parallel mesh partitioning on the root process.
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), iodata.model.partition);

    // Clean up unused domain elements from the mesh, add new boundary elements for material
    // interfaces if not present, and optionally (when running unassembled) add subdomain
    // interface boundary elements.
    std::map<int, std::array<int, 2>> attr_map =
        CheckMesh(smesh, partitioning, iodata, clean, add_bdr, unassembled);
  }

  // Construct the parallel mesh data structure by distributing the serial mesh from the
  // root process. The serial mesh and partitioning are deleted inside.
  std::unique_ptr<mfem::ParMesh> mesh = DistributeMesh(comm, smesh, partitioning);

#if 0
  {
    std::string tmp = iodata.problem.output;
    if (tmp.back() != '/')
    {
      tmp += '/';
    }
    tmp += "tmp/";
    if (Mpi::Root(comm) && !std::filesystem::exists(tmp))
    {
      std::filesystem::create_directories(tmp);
    }
    int width = 1 + static_cast<int>(std::log10(Mpi::Size(comm)-1));
    std::unique_ptr<mfem::Mesh> gsmesh = LoadMesh(iodata.model.mesh);
    std::unique_ptr<int[]> gpartitioning = GetMeshPartitioning(*gsmesh, Mpi::Size(comm));
    mfem::ParMesh gpmesh(comm, *gsmesh, gpartitioning.get(), 0);
    {
      std::string pfile = mfem::MakeParFilename(tmp + "part.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
      fo.precision(MSH_FLT_PRECISION);
      gpmesh.ParPrint(fo);
    }
    {
      std::string pfile = mfem::MakeParFilename(tmp + "final.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
      fo.precision(MSH_FLT_PRECISION);
      mesh->ParPrint(fo);
    }
  }
#endif

  return mesh;
}

void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
{
  // Prepare for uniform and region-based refinement.
  MFEM_VERIFY(mesh.size() == 1,
              "Input mesh vector before refinement has more than a single mesh!");
  int uniform_ref_levels = iodata.model.refinement.uniform_ref_levels;
  int max_region_ref_levels = 0;
  for (const auto &box : iodata.model.refinement.GetBoxes())
  {
    if (max_region_ref_levels < box.ref_levels)
    {
      max_region_ref_levels = box.ref_levels;
    }
  }
  for (const auto &sphere : iodata.model.refinement.GetSpheres())
  {
    if (max_region_ref_levels < sphere.ref_levels)
    {
      max_region_ref_levels = sphere.ref_levels;
    }
  }
  if (iodata.solver.linear.mg_max_levels > 1)
  {
    mesh.reserve(1 + uniform_ref_levels + max_region_ref_levels);
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the coarse mesh so that the refinements are still true refinements of
  // the original mesh (required for geometric multigrid). Otherwise, it happens after
  // refinement.
  if (iodata.model.reorient_tet && mesh.capacity() > 1)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    mesh[0]->ReorientTetMesh();
#pragma GCC diagnostic pop
  }

  // Uniformly refine the mesh further in parallel, saving the level meshes for geometric
  // coarsening later on if desired.
  for (int l = 0; l < uniform_ref_levels; l++)
  {
    if (mesh.capacity() > 1)
    {
      mesh.push_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
    }
    mesh.back()->UniformRefinement();
  }

  // Proceed with region-based refinement, level-by-level for all regions. Currently support
  // box and sphere region shapes. Any overlap between regions is ignored (take the union,
  // don't double-refine).
  if (max_region_ref_levels > 0 &&
      (mesh[0]->MeshGenerator() & 2 || mesh[0]->MeshGenerator() & 4 ||
       mesh[0]->MeshGenerator() & 8))
  {
    // XX TODO: Region-based refinement won't work if the ParMesh has been constructed from
    //          a conforming mesh, but nonconforming refinement is needed. Unclear if the
    //          current mesh distribution scheme will work even for a conforming serial mesh
    //          which is a NCMesh after Mesh::EnsureNCMesh is called.
    MFEM_ABORT("Region-based refinement is currently only supported for simplex meshes!");
  }
  int region_ref_level = 0;
  bool use_nodes = (mesh.back()->GetNodes() != nullptr);
  int ref = use_nodes ? mesh.back()->GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  int dim = mesh.back()->SpaceDimension();
  while (region_ref_level < max_region_ref_levels)
  {
    // Mark elements for refinement in all regions. An element is marked for refinement if
    // any of its vertices are inside any refinement region for the given level.
    mfem::Array<mfem::Refinement> refs;
    for (int i = 0; i < mesh.back()->GetNE(); i++)
    {
      bool refine = false;
      mfem::DenseMatrix pointmat;
      if (use_nodes)
      {
        mfem::ElementTransformation *T = mesh.back()->GetElementTransformation(i);
        mfem::Geometry::Type geo = mesh.back()->GetElementGeometry(i);
        mfem::RefinedGeometry *RefG = mfem::GlobGeometryRefiner.Refine(geo, ref);
        T->Transform(RefG->RefPts, pointmat);
      }
      else
      {
        mfem::Array<int> verts;
        mesh.back()->GetElementVertices(i, verts);
        pointmat.SetSize(dim, verts.Size());
        for (int j = 0; j < verts.Size(); j++)
        {
          const double *coord = mesh.back()->GetVertex(verts[j]);
          for (int d = 0; d < dim; d++)
          {
            pointmat(d, j) = coord[d];
          }
        }
      }
      for (const auto &box : iodata.model.refinement.GetBoxes())
      {
        if (region_ref_level < box.ref_levels)
        {
          for (int j = 0; j < pointmat.Width(); j++)
          {
            // Check if the point is inside the box.
            int d = 0;
            for (; d < pointmat.Height(); d++)
            {
              if (pointmat(d, j) < box.bbmin[d] || pointmat(d, j) > box.bbmax[d])
              {
                break;
              }
            }
            if (d == dim)
            {
              refine = true;
              break;
            }
          }
          if (refine)
          {
            break;
          }
        }
      }
      if (refine)
      {
        refs.Append(mfem::Refinement(i));
        continue;
      }
      for (const auto &sphere : iodata.model.refinement.GetSpheres())
      {
        if (region_ref_level < sphere.ref_levels)
        {
          for (int j = 0; j < pointmat.Width(); j++)
          {
            // Check if the point is inside the sphere.
            double dist = 0.0;
            for (int d = 0; d < pointmat.Height(); d++)
            {
              double s = pointmat(d, j) - sphere.center[d];
              dist += s * s;
            }
            if (dist <= sphere.r * sphere.r)
            {
              refine = true;
              break;
            }
          }
          if (refine)
          {
            break;
          }
        }
      }
      if (refine)
      {
        refs.Append(mfem::Refinement(i));
      }
    }

    // Do the refinement. For tensor element meshes, this may make the mesh nonconforming
    // (adds hanging nodes).
    if (mesh.capacity() > 1)
    {
      mesh.push_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
    }
    mesh.back()->GeneralRefinement(refs, -1);
    region_ref_level++;
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the mesh after refinement if there is a single mesh (doesn't work with
  // h-refinement geometric multigrid).
  if (iodata.model.reorient_tet && mesh.size() == 1)
  {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    mesh[0]->ReorientTetMesh();
#pragma GCC diagnostic pop
  }

  // Print some mesh information.
  mfem::Vector bbmin, bbmax;
  mesh[0]->GetBoundingBox(bbmin, bbmax);
  const double Lc = iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0);
  Mpi::Print(mesh[0]->GetComm(), "\nMesh curvature order: {:d}\nMesh bounding box:\n",
             mesh[0]->GetNodes() ? mesh[0]->GetNodes()->FESpace()->GetMaxElementOrder()
                                 : 1);
  if (mesh[0]->SpaceDimension() == 3)
  {
    Mpi::Print(mesh[0]->GetComm(),
               " (Xmin, Ymin, Zmin) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n"
               " (Xmax, Ymax, Zmax) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n",
               bbmin[0] * Lc, bbmin[1] * Lc, bbmin[2] * Lc, bbmax[0] * Lc, bbmax[1] * Lc,
               bbmax[2] * Lc);
  }
  else
  {
    Mpi::Print(mesh[0]->GetComm(),
               " (Xmin, Ymin) = ({:+.3e}, {:+.3e}) m\n"
               " (Xmax, Ymax) = ({:+.3e}, {:+.3e}) m\n",
               bbmin[0] * Lc, bbmin[1] * Lc, bbmax[0] * Lc, bbmax[1] * Lc);
  }
  Mpi::Print(mesh[0]->GetComm(), "\n{}", (mesh.size() > 1) ? "Coarse " : "");
  mesh[0]->PrintInfo();
  if (mesh.size() > 1)
  {
    Mpi::Print(mesh[0]->GetComm(), "\nRefined ");
    mesh.back()->PrintInfo();
  }
}

void AttrToMarker(int max_attr, const mfem::Array<int> &attrs, mfem::Array<int> &marker)
{
  MFEM_VERIFY(attrs.Size() == 0 || attrs.Max() <= max_attr,
              "Invalid attribute number present (" << attrs.Max() << ")!");
  marker.SetSize(max_attr);
  if (attrs.Size() == 1 && attrs[0] == -1)
  {
    marker = 1;
  }
  else
  {
    marker = 0;
    for (auto attr : attrs)
    {
      MFEM_VERIFY(attr > 0, "Attribute number less than one!");
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeate attribute in attribute list!");
      marker[attr - 1] = 1;
    }
  }
}

void AttrToMarker(int max_attr, const std::vector<int> &attrs, mfem::Array<int> &marker)
{
  MFEM_VERIFY(attrs.empty() || *std::max_element(attrs.begin(), attrs.end()) <= max_attr,
              "Invalid attribute number present ("
                  << *std::max_element(attrs.begin(), attrs.end()) << ")!");
  marker.SetSize(max_attr);
  if (attrs.size() == 1 && attrs[0] == -1)
  {
    marker = 1;
  }
  else
  {
    marker = 0;
    for (auto attr : attrs)
    {
      MFEM_VERIFY(attr > 0, "Attribute number less than one!");
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeate attribute in attribute list!");
      marker[attr - 1] = 1;
    }
  }
}

namespace
{

// Compute the cross product between two 3-vectors expressed as arrays.
std::array<double, 3> Cross(const std::array<double, 3> &n_1,
                            const std::array<double, 3> &n_2)
{
  return {n_1[1] * n_2[2] - n_1[2] * n_2[1], -n_1[0] * n_2[2] + n_1[2] * n_2[0],
          n_1[0] * n_2[1] - n_1[1] * n_2[0]};
};

// Helper for converting a container with begin(), end() to a string for printing.
template <typename T>
std::string String(const T &v)
{
  std::stringstream msg;
  for (auto x : v)
  {
    msg << x << ' ';
  }
  return msg.str();
};
template <typename T, std::size_t N>
std::string String(const std::vector<std::array<T, N>> &v)
{
  std::stringstream msg;
  for (const auto &x : v)
  {
    msg << "{";
    for (auto y : x)
      msg << y << ' ';
    msg << "} ";
  }
  return msg.str();
};

}  // namespace

double BoundingBox::Area() const
{
  auto cross = Cross(normals[0], normals[1]);
  return 4 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
}

double BoundingBox::Volume() const
{
  // Compute area in plane, then extrude above and below to get volume.
  return planar
             ? 0
             : 2 *
                   std::sqrt(normals[2][0] * normals[2][0] + normals[2][1] * normals[2][1] +
                             normals[2][2] * normals[2][2]) *
                   Area();
}

BoundingBox BoundingBoxFromPointCloud(MPI_Comm comm,
                                      std::vector<std::array<double, 3>> vertices)
{
  using Array3 = std::array<double, 3>;

  const auto num_vertices = int(vertices.size());
  const int dominant_rank = Mpi::FindMax<int>(num_vertices, comm).rank;

  // dominant_rank will perform the calculation
  std::vector<int> recv_counts(Mpi::Size(comm)), displacements;
  MPI_Gather(&num_vertices, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, dominant_rank,
             comm);

  std::vector<std::array<double, 3>> collected_vertices;
  if (Mpi::Rank(comm) == dominant_rank)
  {
    // First displacement is zero, then after is the partial sum of recv_counts.
    displacements.resize(Mpi::Size(comm));
    displacements[0] = 0;
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);

    // Add on slots at the end of vertices for the incoming data.
    collected_vertices.resize(std::accumulate(recv_counts.begin(), recv_counts.end(), 0));
  }

  // Gather the data to the dominant rank.
  const auto mpi_array_3_type = Mpi::GetArrayType<double, 3>();
  MPI_Gatherv(vertices.data(), num_vertices, mpi_array_3_type, collected_vertices.data(),
              recv_counts.data(), displacements.data(), mpi_array_3_type, dominant_rank,
              comm);

  BoundingBox box;
  if (Mpi::Rank(comm) == dominant_rank)
  {
    vertices = std::move(collected_vertices);
    // dominant_rank now has a fully stocked vector of vertices. All other ranks
    // will wait for results.
    // Deduplicate vertices. Given floating point precision, need a tolerance.
    auto vertex_equality = [](const auto &x, const auto &y)
    {
      constexpr double tolerance = 10 * std::numeric_limits<double>::epsilon();
      return std::abs(x[0] - y[0]) < tolerance && std::abs(x[1] - y[1]) < tolerance &&
             std::abs(x[2] - y[2]) < tolerance;
    };

    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end(), vertex_equality),
                   vertices.end());

    MFEM_VERIFY(vertices.size() >= 4,
                "A bounding box requires a minimum of four vertices for this algorithm");

    // Define a diagonal of the ASSUMED cuboid bounding box. Store references as
    // this is useful for checking pointers later.
    const auto &v_000 = *std::min_element(vertices.begin(), vertices.end());
    const auto &v_111 = *std::max_element(vertices.begin(), vertices.end());

    MFEM_VERIFY(&v_000 != &v_111, "Minimum and Maximum extents cannot be identical");

    const auto origin = v_000;  // Save a copy to make undoing the transform later easier.

    // Return y <- a * x + y
    auto Inplace = [](Array3 &y, const Array3 &x, double a = 1.0)
    {
      y[0] += a * x[0];
      y[1] += a * x[1];
      y[2] += a * x[2];
    };

    // Dot product
    auto Dot = [](const Array3 &x_1, const Array3 &x_2)
    { return x_1[0] * x_2[0] + x_1[1] * x_2[1] + x_1[2] * x_2[2]; };

    // Scale to unit magnitude
    auto Normalize = [Dot](Array3 &v)
    {
      auto mag = std::sqrt(Dot(v, v));
      MFEM_ASSERT(mag >= 0, "!");
      v[0] /= mag;
      v[1] /= mag;
      v[2] /= mag;
    };

    // Given a pair of vertices, define a normalized vector from x_1 to x_2
    auto Normal = [Inplace, Normalize](const Array3 &x_1, Array3 x_2)
    {
      Inplace(x_2, x_1, -1);
      Normalize(x_2);
      return x_2;
    };

    // v_1 -= (v_1 . v_2) v_2. ASSUMES v_2 normalized
    auto Orthogonalize = [Inplace, Dot](Array3 &v_1, const Array3 &v_2)
    { Inplace(v_1, v_2, -Dot(v_1, v_2)); };

    // Orient all vertex positions relative to origin. Makes projection and
    // orientation logic simpler.
    for (auto &v : vertices)
    {
      Inplace(v, origin, -1.0);
    }

    const auto n_1 = Normal(v_000, v_111);

    // Compute the distance from the normal axis. Note: everything has been
    // oriented relative to v_000 == (0,0,0).
    auto PerpendicularDistance = [&n_1, Dot, Inplace](Array3 v)
    {
      Inplace(v, n_1, -Dot(v, n_1));
      return std::sqrt(Dot(v, v));
    };

    // Find the vertex furthest from the diagonal axis. We cannot know yet if
    // this defines (001) or (011).
    const auto &t_0 =
        *std::max_element(vertices.begin(), vertices.end(),
                          [PerpendicularDistance](const Array3 &x, const Array3 &y)
                          { return PerpendicularDistance(x) < PerpendicularDistance(y); });

    MFEM_VERIFY(&t_0 != &v_000, "Vertices are degenerate!");
    MFEM_VERIFY(&t_0 != &v_111, "Vertices are degenerate!");

    // Use the discovered vertex define a second normal direction, and thus a plane.
    auto n_2 = t_0;
    Orthogonalize(n_2, n_1);
    Normalize(n_2);

    // n_1 and n_2 now define a planar coordinate system intersecting the main
    // diagonal, and two opposite edges of the cuboid.

    // Now look for a component that maximizes distance from the planar system:
    // complete the axes with a cross, then use a dot product to pick the
    // greatest deviation.
    const auto n_3 = Cross(n_1, n_2);

    auto OutOfPlaneComp = [&n_1, &n_2, &n_3, Dot](const Array3 &v_1, const Array3 &v_2)
    {
      // Precedence of directions is in reverse order of discovery of the
      // directions. The most important deciding feature is the distance in the
      // out of plane direction, then in the first off-diagonal, then finally in
      // the diagonal direction.
      const Array3 dist_1{{Dot(v_1, n_3), Dot(v_1, n_2), Dot(v_1, n_1)}};
      const Array3 dist_2{{Dot(v_2, n_3), Dot(v_2, n_2), Dot(v_2, n_1)}};
      return dist_1 < dist_2;
    };
    const auto &t_1 = *std::max_element(vertices.begin(), vertices.end(), OutOfPlaneComp);

    // There is a degeneration if the final point is within the plane defined by
    // (v_000, v_001, v_111).
    constexpr double planar_tolerance = 1e-9;
    box.planar = std::abs(Dot(n_3, t_0)) < planar_tolerance &&
                 std::abs(Dot(n_3, t_1)) < planar_tolerance;

    if (!box.planar)
    {
      MFEM_VERIFY(&t_0 != &t_1, "Degenerate coordinates");

      MFEM_VERIFY(&t_0 != &v_000, "Degenerate coordinates");
      MFEM_VERIFY(&t_0 != &v_111, "Degenerate coordinates");

      MFEM_VERIFY(&t_1 != &v_000, "Degenerate coordinates");
      MFEM_VERIFY(&t_1 != &v_111, "Degenerate coordinates");
    }

    // If t_1 points to v_000, t_0 or v_111, then the data is coplanar.
    // Establish if t_0 is a diagonal or not (using Pythagoras).
    // Only pick t_1 for v_001 if the points are non planar, and t_0 is longer.
    bool t_0_gt_t_1 = Dot(t_0, t_0) >= Dot(t_1, t_1);
    const auto &v_001 = !box.planar && t_0_gt_t_1 ? t_1 : t_0;
    const auto &v_011 = !box.planar && t_0_gt_t_1 ? t_0 : t_1;

    // In the v_000 coordinate system, can easily compute the center as halfway
    // along the main diagonal.
    box.center = {origin[0] + v_111[0] / 2, origin[1] + v_111[1] / 2,
                  origin[2] + v_111[2] / 2};

    // The length in each direction is then given by traversing the edges of the
    // cuboid in turn
    box.normals[0] = Array3{(v_001[0] - v_000[0]) / 2, (v_001[1] - v_000[1]) / 2,
                            (v_001[2] - v_000[2]) / 2};
    box.normals[1] = box.planar
                         ? Array3{(v_111[0] - v_001[0]) / 2, (v_111[1] - v_001[1]) / 2,
                                  (v_111[2] - v_001[2]) / 2}
                         : Array3{(v_011[0] - v_001[0]) / 2, (v_011[1] - v_001[1]) / 2,
                                  (v_011[2] - v_001[2]) / 2};
    box.normals[2] = box.planar
                         ? Array3{0, 0, 0}
                         : Array3{(v_111[0] - v_011[0]) / 2, (v_111[1] - v_011[1]) / 2,
                                  (v_111[2] - v_011[2]) / 2};

    // Make sure the longest dimension comes first
    std::sort(box.normals.begin(), box.normals.end(),
              [Dot](const Array3 &x, const Array3 &y) { return Dot(x, x) > Dot(y, y); });
  }

  Mpi::Broadcast(3, box.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, box.normals.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &box.planar, dominant_rank, comm);

  return box;
}

BoundingBox GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr)
{
  // Collect set of vertices to check
  std::set<int> vertex_indices;

  using Array3 = std::array<double, 3>;
  std::vector<Array3> vertices;

  if (mesh.GetNodes() == nullptr)
  {
    // Linear mesh, work with element vertices directly.
    mfem::Array<int> v;
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); ++i)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetBdrElementVertices(i, v);
        vertex_indices.insert(v.begin(), v.end());
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetElementVertices(i, v);
        vertex_indices.insert(v.begin(), v.end());
      }
    }

    for (auto i : vertex_indices)
    {
      const auto &vx = mesh.GetVertex(i);
      vertices.push_back({vx[0], vx[1], vx[2]});
    }
  }
  else
  {
    // Nonlinear mesh, need to process point matrices.
    const int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    mfem::DenseMatrix pointmat;  // 3 x N
    mfem::ElementTransformation *T;

    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); ++i)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }

        T = mesh.GetBdrElementTransformation(i);
        T->Transform(
            mfem::GlobGeometryRefiner.Refine(mesh.GetBdrElementGeometry(i), ref)->RefPts,
            pointmat);
        for (int j = 0; j < pointmat.Width(); ++j)
        {
          vertices.push_back({pointmat(0, j), pointmat(1, j), pointmat(2, j)});
        }
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); ++i)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }

        T = mesh.GetElementTransformation(i);
        T->Transform(
            mfem::GlobGeometryRefiner.Refine(mesh.GetElementGeometry(i), ref)->RefPts,
            pointmat);
        for (int j = 0; j < pointmat.Width(); ++j)
        {
          vertices.push_back({pointmat(0, j), pointmat(1, j), pointmat(2, j)});
        }
      }
    }
  }

  return BoundingBoxFromPointCloud(mesh.GetComm(), std::move(vertices));
}

BoundingBox GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;

  return GetBoundingBox(mesh, marker, bdr);
}

BoundingBall BoundingBallFromPointCloud(MPI_Comm comm,
                                        std::vector<std::array<double, 3>> vertices)
{
  using Array3 = std::array<double, 3>;

  auto num_vertices = int(vertices.size());
  const int dominant_rank = Mpi::FindMax<int>(num_vertices, comm).rank;

  const bool is_root = dominant_rank == Mpi::Rank(comm);

  // dominant_rank will perform the calculation
  std::vector<int> recv_counts(Mpi::Size(comm)), displacements;
  MPI_Gather(&num_vertices, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, dominant_rank,
             comm);

  std::vector<Array3> collected_vertices;
  if (is_root)
  {
    // First displacement is zero, then after is the partial sum of recv_counts.
    displacements.resize(Mpi::Size(comm));
    displacements[0] = 0;
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);

    // Add on slots at the end of vertices for the incoming data.
    collected_vertices.resize(std::accumulate(recv_counts.begin(), recv_counts.end(), 0));
  }

  // Gather the data to the dominant rank.
  const auto mpi_array_3_type = Mpi::GetArrayType<double, 3>();
  MPI_Gatherv(vertices.data(), vertices.size(), mpi_array_3_type, collected_vertices.data(),
              recv_counts.data(), displacements.data(), mpi_array_3_type, dominant_rank,
              comm);

  BoundingBall ball;
  if (Mpi::Rank(comm) == dominant_rank)
  {
    vertices = std::move(collected_vertices);
    // dominant_rank now has a fully stocked vector of vertices.
    // Use lexicographic sort to get min and max vertices. These will be
    // separated by the largest L1 distance.

    const auto &min = *std::min_element(vertices.begin(), vertices.end());
    const auto &max = *std::max_element(vertices.begin(), vertices.end());

    Array3 delta{{max[0] - min[0], max[1] - min[1], max[2] - min[2]}};

    ball.radius =
        std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]) / 2;
    ball.center = {min[0] + delta[0] / 2, min[1] + delta[1] / 2, min[2] + delta[2] / 2};

    // Project onto this candidate diameter, and pick a vertex furthest away.
    // Check that this resulting distance is less than or equal to the radius,
    // and use the resulting direction to compute another in plane vector.

    // Assumes all delta are normalized, and applies a common origin as part of the
    // projection.
    auto PerpendicularDistance = [min](const std::initializer_list<Array3> &deltas, auto v)
    {
      v[0] -= min[0];
      v[1] -= min[1];
      v[2] -= min[2];
      for (const auto &d : deltas)
      {
        double dot = d[0] * v[0] + d[1] * v[1] + d[2] * v[2];
        for (auto i : {0, 1, 2})
        {
          v[i] -= dot * d[i];
        }
      }
      return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    };

    // Normalize the diagonal.
    delta[0] /= (ball.radius * 2);
    delta[1] /= (ball.radius * 2);
    delta[2] /= (ball.radius * 2);

    const auto &perp = *std::max_element(
        vertices.begin(), vertices.end(),
        [&delta, PerpendicularDistance](const auto &x, const auto &y)
        { return PerpendicularDistance({delta}, x) < PerpendicularDistance({delta}, y); });

    constexpr double radius_tol = 1e-6;
    MFEM_VERIFY(std::abs(PerpendicularDistance({delta}, perp) - ball.radius) /
                        ball.radius <=
                    radius_tol,
                "A point perpendicular must be contained in the ball: "
                    << PerpendicularDistance({delta}, perp) << " vs " << ball.radius);

    auto n_radial = perp;
    n_radial[0] -= ball.center[0];
    n_radial[1] -= ball.center[1];
    n_radial[2] -= ball.center[2];

    auto radius_2 = std::sqrt(n_radial[0] * n_radial[0] + n_radial[1] * n_radial[1] +
                              n_radial[2] * n_radial[2]);
    n_radial[0] /= radius_2;
    n_radial[1] /= radius_2;
    n_radial[2] /= radius_2;

    // Compute a perpendicular to the circle using the cross product
    ball.planar_normal[0] = delta[1] * n_radial[2] - delta[2] * n_radial[1];
    ball.planar_normal[1] = delta[2] * n_radial[0] - delta[0] * n_radial[2];
    ball.planar_normal[2] = delta[0] * n_radial[1] - delta[1] * n_radial[0];

    // Compute the point furthest out of the plane discovered. If below
    // tolerance, this means the circle is 2D.

    const auto &out_of_plane = *std::max_element(
        vertices.begin(), vertices.end(),
        [&delta, &n_radial, PerpendicularDistance](const auto &x, const auto &y)
        {
          return PerpendicularDistance({delta, n_radial}, x) <
                 PerpendicularDistance({delta, n_radial}, y);
        });

    constexpr double planar_tolerance = 1e-9;
    if (PerpendicularDistance({delta, n_radial}, out_of_plane) < planar_tolerance)
    {
      // The points are functionally coplanar, zero out the normal
      ball.planar_normal[0] *= 0;
      ball.planar_normal[1] *= 0;
      ball.planar_normal[2] *= 0;
    }
    else
    {
      MFEM_VERIFY(PerpendicularDistance({delta, n_radial}, out_of_plane) <= ball.radius,
                  "A point perpendicular must be contained in the ball");
    }
  }

  Mpi::Broadcast(3, ball.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3, ball.planar_normal.data(), dominant_rank, comm);
  Mpi::Broadcast(1, &ball.radius, dominant_rank, comm);

  return ball;
}

BoundingBall GetBoundingBall(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr)
{
  // Collect set of vertices to check.
  std::set<int> vertex_indices;
  mfem::Array<int> v;
  if (bdr)
  {
    for (int i = 0; i < mesh.GetNBE(); ++i)
    {
      if (!marker[mesh.GetBdrAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetBdrElementVertices(i, v);
      vertex_indices.insert(v.begin(), v.end());
    }
  }
  else
  {
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      if (!marker[mesh.GetAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetElementVertices(i, v);
      vertex_indices.insert(v.begin(), v.end());
    }
  }

  // Extract the coordinates for each vertex
  using Array3 = std::array<double, 3>;

  std::vector<Array3> vertices;
  for (auto i : vertex_indices)
  {
    const auto &vx = mesh.GetVertex(i);
    vertices.push_back({vx[0], vx[1], vx[2]});
  }

  return BoundingBallFromPointCloud(mesh.GetComm(), vertices);
}

BoundingBall GetBoundingBall(mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;

  return GetBoundingBall(mesh, marker, bdr);
}

void GetCartesianBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                             mfem::Vector &max)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetCartesianBoundingBox(mesh, marker, bdr, min, max);
}

void GetCartesianBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                             mfem::Vector &min, mfem::Vector &max)
{
  int dim = mesh.SpaceDimension();
  min.SetSize(dim);
  max.SetSize(dim);
  for (int d = 0; d < dim; d++)
  {
    min(d) = mfem::infinity();
    max(d) = -mfem::infinity();
  }
  if (mesh.GetNodes() == nullptr)
  {
    auto BBUpdate = [&mesh, &dim, &min, &max](mfem::Array<int> &verts) -> void
    {
      for (int j = 0; j < verts.Size(); j++)
      {
        const double *coord = mesh.GetVertex(verts[j]);
        for (int d = 0; d < dim; d++)
        {
          if (coord[d] < min(d))
          {
            min(d) = coord[d];
          }
          if (coord[d] > max(d))
          {
            max(d) = coord[d];
          }
        }
      }
    };
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mfem::Array<int> verts;
        mesh.GetBdrElementVertices(i, verts);
        BBUpdate(verts);
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mfem::Array<int> verts;
        mesh.GetElementVertices(i, verts);
        BBUpdate(verts);
      }
    }
  }
  else
  {
    int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    auto BBUpdate = [&ref, &min, &max](mfem::ElementTransformation *T,
                                       mfem::Geometry::Type &geo) -> void
    {
      mfem::DenseMatrix pointmat;
      mfem::RefinedGeometry *RefG = mfem::GlobGeometryRefiner.Refine(geo, ref);
      T->Transform(RefG->RefPts, pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        for (int d = 0; d < pointmat.Height(); d++)
        {
          if (pointmat(d, j) < min(d))
          {
            min(d) = pointmat(d, j);
          }
          if (pointmat(d, j) > max(d))
          {
            max(d) = pointmat(d, j);
          }
        }
      }
    };
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
        mfem::Geometry::Type geo = mesh.GetBdrElementGeometry(i);
        BBUpdate(T, geo);
      }
    }
    else
    {
      for (int i = 0; i < mesh.GetNE(); i++)
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
        mfem::Geometry::Type geo = mesh.GetElementGeometry(i);
        BBUpdate(T, geo);
      }
    }
  }
  auto *Min = min.HostReadWrite();
  auto *Max = max.HostReadWrite();
  Mpi::GlobalMin(dim, Min, mesh.GetComm());
  Mpi::GlobalMax(dim, Max, mesh.GetComm());
}

void GetSurfaceNormal(mfem::ParMesh &mesh, int attr, mfem::Vector &normal)
{
  mfem::Array<int> marker(mesh.bdr_attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetSurfaceNormal(mesh, marker, normal);
}

void GetSurfaceNormal(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                      mfem::Vector &normal)
{
  int dim = mesh.SpaceDimension();
  mfem::Vector nor(dim);
  normal.SetSize(dim);
  normal = 0.0;
  bool init = false;
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    if (!marker[mesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
    const mfem::IntegrationPoint &ip =
        mfem::Geometries.GetCenter(mesh.GetBdrElementGeometry(i));
    T->SetIntPoint(&ip);
    mfem::CalcOrtho(T->Jacobian(), nor);
    if (!init)
    {
      normal = nor;
      init = true;
    }
    else
    {
      // Check orientation and make sure consistent on this process. If a boundary has
      // conflicting normal definitions, use the first value.
      if (nor * normal < 0.0)
      {
        normal -= nor;
      }
      else
      {
        normal += nor;
      }
    }
  }
  // If different processors have different normal orientations, take that from the lowest
  // rank processor.
  MPI_Comm comm = mesh.GetComm();
  int rank = Mpi::Size(comm);
  mfem::Vector glob_normal(dim);
  if (init)
  {
    rank = Mpi::Rank(comm);
  }
  Mpi::GlobalMin(1, &rank, comm);
  if (rank == Mpi::Size(comm))
  {
    // No boundary elements of attribute attr.
    normal = 0.0;
    return;
  }
  if (rank == Mpi::Rank(comm))
  {
    glob_normal = normal;
  }
  {
    auto *GlobNormal = glob_normal.HostReadWrite();
    Mpi::Broadcast(dim, GlobNormal, rank, comm);
  }
  if (init && normal * glob_normal < 0.0)
  {
    normal.Neg();
  }
  {
    auto *Normal = normal.HostReadWrite();
    Mpi::GlobalSum(dim, Normal, comm);
  }
  normal /= normal.Norml2();
  // if (dim == 3)
  // {
  //   Mpi::Print(comm, " Surface normal {:d} = ({:+.3e}, {:+.3e}, {:+.3e})", attr,
  //   normal(0),
  //              normal(1), normal(2));
  // }
  // else
  // {
  //   Mpi::Print(comm, " Surface normal {:d} = ({:+.3e}, {:+.3e})", attr, normal(0),
  //              normal(1));
  // }
}

}  // namespace mesh

namespace
{

std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &path)
{
  // Read the (serial) mesh from the given mesh file. Handle preparation for refinement and
  // orientations here to avoid possible reorientations and reordering later on. MFEM
  // supports a native mesh format (.mesh), VTK/VTU, Gmsh, as well as some others. We use
  // built-in converters for the types we know, otherwise rely on MFEM to do the conversion
  // or error out if not supported.
  std::filesystem::path mfile(path);
  if (mfile.extension() == ".mphtxt" || mfile.extension() == ".mphbin" ||
      mfile.extension() == ".nas" || mfile.extension() == ".bdf")
  {
    // Put translated mesh in temporary string buffer.
    std::stringstream fi(std::stringstream::in | std::stringstream::out);
    // fi << std::fixed;
    fi << std::scientific;
    fi.precision(MSH_FLT_PRECISION);

#if 0
    // Put translated mesh in temporary storage (directory is created and destroyed in
    // calling function).
    std::string tmp = iodata.problem.output;
    if (tmp.back() != '/')
    {
      tmp += '/';
    }
    tmp += "tmp/serial.msh";
    std::ofstream fo(tmp);
    // mfem::ofgzstream fo(tmp, true);  // Use zlib compression if available
    // fo << std::fixed;
    fo << std::scientific;
    fo.precision(MSH_FLT_PRECISION);
#endif

    if (mfile.extension() == ".mphtxt" || mfile.extension() == ".mphbin")
    {
      mesh::ConvertMeshComsol(path, fi);
      // mesh::ConvertMeshComsol(path, fo);
    }
    else
    {
      mesh::ConvertMeshNastran(path, fi);
      // mesh::ConvertMeshNastran(path, fo);
    }

#if 0
    std::ifstream fi(tmp);
    // mfem::ifgzstream fi(tmp);
    if (!fi.good())
    {
      MFEM_ABORT("Unable to open translated mesh file \"" << tmp << "\"!");
    }
#endif
    return std::make_unique<mfem::Mesh>(fi, 1, 1, true);
  }
  // Otherwise, just rely on MFEM load the mesh.
  std::ifstream fi(path);
  if (!fi.good())
  {
    MFEM_ABORT("Unable to open mesh file \"" << path << "\"!");
  }
  std::unique_ptr mesh = std::make_unique<mfem::Mesh>(fi, 1, 1, true);
  mesh->EnsureNodes();
  return mesh;
}

void ReorderMesh(mfem::Mesh &mesh)
{
  mfem::Array<int> ordering;

#if 0
  // Gecko reordering.
  Mpi::Print(\n);
  mfem::Array<int> tentative;
  int outer = 3, inner = 3, window = 4, period = 2;
  double best_cost = mfem::infinity();
  for (int i = 0; i < outer; i++)
  {
    int seed = i+1;
    double cost = mesh.GetGeckoElementOrdering(tentative, inner, window, eriod, seed, true);
    if (cost < best_cost)
    {
      ordering = tentative;
      best_cost = cost;
    }
  }
  Mpi::Print("Final cost: {:e}\n", best_cost);
#endif

  // (Faster) Hilbert reordering.
  mesh.GetHilbertElementOrdering(ordering);
  mesh.ReorderElements(ordering);
}

std::unique_ptr<int[]> GetMeshPartitioning(mfem::Mesh &mesh, int size,
                                           const std::string &partition)
{
  MFEM_VERIFY(size <= mesh.GetNE(), "Mesh partitioning must have parts <= mesh elements ("
                                        << size << " vs. " << mesh.GetNE() << ")!");
  if (partition.length() == 0)
  {
    const int part_method = 1;
    std::unique_ptr<int[]> partitioning(mesh.GeneratePartitioning(size, part_method));
    Mpi::Print("Finished partitioning mesh into {:d} subdomain{}\n", size,
               (size > 1) ? "s" : "");
    return partitioning;
  }
  // User can optionally specify a mesh partitioning file as generated from the MFEM
  // mesh-explorer miniapp, for example. It has the format:
  //
  //   number_of_elements <NE>
  //   number_of_processors <NPART>
  //   <part[0]>
  //     ...
  //   <part[NE-1]>
  //
  int nel, np;
  std::ifstream part_ifs(partition);
  part_ifs.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
  part_ifs >> nel;
  if (nel != mesh.GetNE())
  {
    MFEM_ABORT("Invalid partitioning file (number of elements)!");
  }
  part_ifs.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
  part_ifs >> np;
  if (np != size)
  {
    MFEM_ABORT("Invalid partitioning file (number of processors)!");
  }
  auto partitioning = std::make_unique<int[]>(mesh.GetNE());
  int i = 0;
  while (i < mesh.GetNE())
  {
    part_ifs >> partitioning[i++];
  }
  Mpi::Print("Read mesh partitioning into {:d} subdomain{} from disk\n", size,
             (size > 1) ? "s" : "");
  return partitioning;
}

std::map<int, std::array<int, 2>> CheckMesh(std::unique_ptr<mfem::Mesh> &orig_mesh,
                                            const std::unique_ptr<int[]> &partitioning,
                                            const IoData &iodata, bool clean_elem,
                                            bool add_bdr, bool add_subdomain)
{
  // - Check that all external boundaries of the mesh have a corresponding boundary
  //   condition.
  // - If desired, create a new mesh which has added boundary elements for all material
  //   interfaces if these elements do not yet exist.
  // - If desired, create a new mesh which has removed all domain elements which do not have
  //   an associated material property specified in the input file.
  MFEM_VERIFY(orig_mesh->Dimension() == 3 && !orig_mesh->Nonconforming(),
              "Nonconforming or 2D meshes have not been tested yet!");
  mfem::Array<int> mat_marker, bdr_marker;
  GetUsedAttributeMarkers(iodata, orig_mesh->attributes.Max(),
                          orig_mesh->bdr_attributes.Max(), mat_marker, bdr_marker);
  bool warn = false;
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    int attr = orig_mesh->GetBdrAttribute(be);
    if (!bdr_marker[attr - 1])
    {
      int f, o, e1, e2;
      orig_mesh->GetBdrElementFace(be, &f, &o);
      orig_mesh->GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)  // Internal boundary elements are allowed to have no BC
      {
        warn = true;
        break;
      }
    }
  }
  if (warn)
  {
    Mpi::Warning("One or more external boundary attributes has no associated boundary "
                 "condition!\n\"PMC\"/\"ZeroCharge\" condition is assumed!\n");
  }

  // Mapping from new interface boundary attribute tags to vector of neighboring domain
  // attributes (when adding new boundary elements).
  std::map<int, std::array<int, 2>> new_attr_map;
  if (!clean_elem && !add_bdr && !add_subdomain)
  {
    return new_attr_map;
  }

  // Count deleted or added domain and boundary elements.
  int new_ne = orig_mesh->GetNE();
  int new_nbdr = orig_mesh->GetNBE();
  mfem::Array<bool> elem_delete, bdr_delete;
  mfem::Array<int> orig_bdr_faces, add_bdr_faces;
  elem_delete.SetSize(orig_mesh->GetNE(), false);
  bdr_delete.SetSize(orig_mesh->GetNBE(), false);
  orig_bdr_faces.SetSize(orig_mesh->GetNumFaces(), -1);
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    int f, o;
    orig_mesh->GetBdrElementFace(be, &f, &o);
    MFEM_VERIFY(orig_bdr_faces[f] < 0,
                "Mesh should not define boundary elements multiple times!");
    orig_bdr_faces[f] = be;
  }
  if (add_bdr || add_subdomain)
  {
    add_bdr_faces.SetSize(orig_mesh->GetNumFaces(), -1);
  }

  if (clean_elem)
  {
    // Delete domain and boundary elements which have no associated material or BC attribute
    // from the mesh.
    for (int e = 0; e < orig_mesh->GetNE(); e++)
    {
      int attr = orig_mesh->GetAttribute(e);
      if (!mat_marker[attr - 1])
      {
        elem_delete[e] = true;
        new_ne--;
      }
    }

    // Make sure to remove any boundary elements which are no longer attached to elements in
    // the domain.
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be >= 0)
      {
        int e1, e2;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        if ((e1 < 0 || elem_delete[e1]) && (e2 < 0 || elem_delete[e2]))
        {
          // Mpi::Print("Deleting an unattached boundary element!\n");
          bdr_delete[be] = true;
          new_nbdr--;
        }
      }
    }
    if (new_ne < orig_mesh->GetNE())
    {
      Mpi::Print("Removed {:d} unmarked domain elements from the mesh\n",
                 orig_mesh->GetNE() - new_ne);
    }
    if (new_nbdr < orig_mesh->GetNBE())
    {
      Mpi::Print("Removed {:d} unattached boundary elements from the mesh\n",
                 orig_mesh->GetNBE() - new_nbdr);
    }
  }
  int new_ne_step1 = new_ne;
  int new_nbdr_step1 = new_nbdr;

  if (add_bdr)
  {
    // Add new boundary elements at material interfaces or on the exterior boundary of the
    // simulation domain, if there is not already a boundary element present.
    MFEM_VERIFY(!orig_mesh->Nonconforming(), "Adding material interface boundary elements "
                                             "is not supported for nonconforming meshes!");
    int add_bdr_ext = 0, add_bdr_int = 0;
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        bool no_e1 = (e1 < 0 || elem_delete[e1]);
        bool no_e2 = (e2 < 0 || elem_delete[e2]);
        if ((no_e1 || no_e2) && !(no_e1 && no_e2))
        {
          // Mpi::Print("Adding exterior boundary element!\n");
          add_bdr_faces[f] = 1;
          add_bdr_ext++;
        }
        else if (orig_mesh->GetAttribute(e1) != orig_mesh->GetAttribute(e2))
        {
          // Add new boundary element at material interface between two domains.
          // Mpi::Print("Adding material interface boundary element!\n");
          add_bdr_faces[f] = 1;
          add_bdr_int++;
        }
      }
    }
    new_nbdr += (add_bdr_ext + add_bdr_int);
    if (add_bdr_ext > 0)
    {
      Mpi::Print("Added {:d} boundary elements for exterior boundaries to the mesh\n",
                 add_bdr_ext);
    }
    if (add_bdr_int > 0)
    {
      Mpi::Print("Added {:d} boundary elements for material interfaces to the mesh\n",
                 add_bdr_int);
    }
  }
  int new_ne_step2 = new_ne;
  int new_nbdr_step2 = new_nbdr;

  if (add_subdomain)
  {
    // Add new boundary elements at interfaces between elements beloning to different
    // subdomains. This uses similar code to mfem::Mesh::PrintWithPartitioning.
    MFEM_VERIFY(partitioning, "Cannot add subdomain interface boundary elements without "
                              "supplied mesh partitioning!");
    MFEM_VERIFY(!orig_mesh->Nonconforming(), "Adding subdomain interface boundary elements "
                                             "is not supported for nonconforming meshes!");
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        bool no_e1 = (e1 < 0 || elem_delete[e1]);
        bool no_e2 = (e2 < 0 || elem_delete[e2]);
        if (!no_e1 && !no_e2 && partitioning[e1] != partitioning[e2])
        {
          // Internal face is connected to two elements belonging to different subdomains
          // (this works for conforming meshes).
          add_bdr_faces[f] = 2;
          new_nbdr += 2;
        }
      }
      // else
      // {
      //   // This face is attached to a boundary element. We could define a new boundary
      //   // element with opposite orientation to ensure both subdomains in the distributed
      //   // ParMesh have the boundary element.
      // }
    }
    if (new_nbdr > new_nbdr_step2)
    {
      Mpi::Print("Added boundary elements for subdomain interfaces to the mesh\n",
                 new_nbdr - new_nbdr_step2);
    }
  }

  // Create the new mesh.
  if (new_ne == new_ne_step1 && new_ne_step1 == new_ne_step2 &&
      new_ne_step2 == orig_mesh->GetNE() && new_nbdr == new_nbdr_step1 &&
      new_nbdr_step1 == new_nbdr_step2 && new_nbdr_step2 == orig_mesh->GetNBE())
  {
    return new_attr_map;
  }
  std::unique_ptr<mfem::Mesh> new_mesh =
      std::make_unique<mfem::Mesh>(orig_mesh->Dimension(), orig_mesh->GetNV(), new_ne,
                                   new_nbdr, orig_mesh->SpaceDimension());

  // Copy vertices and non-deleted domain and boundary elements.
  for (int v = 0; v < orig_mesh->GetNV(); v++)
  {
    new_mesh->AddVertex(orig_mesh->GetVertex(v));
  }
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    if (!elem_delete[e])
    {
      mfem::Element *ne = orig_mesh->GetElement(e)->Duplicate(new_mesh.get());
      new_mesh->AddElement(ne);
    }
  }
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    if (!bdr_delete[be])
    {
      mfem::Element *ne = orig_mesh->GetBdrElement(be)->Duplicate(new_mesh.get());
      new_mesh->AddBdrElement(ne);
    }
  }

  // Add new boundary elements.
  if (add_bdr || add_subdomain)
  {
    auto FlipVertices = [](mfem::Element *e)
    {
      mfem::Array<int> v;
      e->GetVertices(v);
      int start = 0, end = v.Size() - 1;
      while (start < end)
      {
        int t = v[start];
        v[start] = v[end];
        v[end] = t;
        start++;
        end--;
      }
      e->SetVertices(v.HostRead());
    };

    // 1-based, some boundary attributes may be empty since they were removed from the
    // original mesh, but to keep indices the same as config file we don't compact the
    // list.
    int max_bdr_attr = orig_mesh->bdr_attributes.Max();
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      if (add_bdr_faces[f] > 0)
      {
        // Assign new unique attribute based on attached elements (we want the material
        // properties on the face to average those on the elements). This is used later on
        // when integrating the transmission condition on the subdomain interface. Save the
        // inverse so that the attributes of e1 and e2 can be easily referenced using the
        // new attribute. Since attributes are in 1-based indexing, a, b > 0.
        int e1, e2, a = 0, b = 0;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        bool no_e1 = (e1 < 0 || elem_delete[e1]);
        bool no_e2 = (e2 < 0 || elem_delete[e2]);
        if (!no_e1 && !no_e2)
        {
          a = std::max(orig_mesh->GetAttribute(e1), orig_mesh->GetAttribute(e2));
          b = (a == orig_mesh->GetAttribute(e1)) ? orig_mesh->GetAttribute(e2)
                                                 : orig_mesh->GetAttribute(e1);
        }
        else if (!no_e1)
        {
          a = orig_mesh->GetAttribute(e1);
          b = 0;
        }
        else if (!no_e2)
        {
          a = orig_mesh->GetAttribute(e2);
          b = 0;
        }
        MFEM_VERIFY(a + b > 0, "Invalid new boundary element attribute!");
        int new_attr = max_bdr_attr + (a * (a - 1)) / 2 + b;  // At least max_bdr_attr+1
        if (new_attr_map.find(new_attr) == new_attr_map.end())
        {
          new_attr_map.emplace(new_attr, std::array<int, 2>{a, b});
        }

        // Add the boundary elements with the new boundary attribute.
        mfem::Element *ne = orig_mesh->GetFace(f)->Duplicate(new_mesh.get());
        ne->SetAttribute(new_attr);
        new_mesh->AddBdrElement(ne);
        if (add_bdr_faces[f] > 1)
        {
          // Flip order of vertices to reverse normal direction of second added element.
          ne = orig_mesh->GetFace(f)->Duplicate(new_mesh.get());
          FlipVertices(ne);
          ne->SetAttribute(new_attr);
          new_mesh->AddBdrElement(ne);
          // Mpi::Print("Adding two BE with attr {:d} from elements {:d} and {:d}\n",
          //            new_attr, a, b);
        }
      }
    }
  }

  // Finalize new mesh and replace the old one. If a curved mesh, set up the new mesh by
  // projecting nodes onto the new mesh for the non-trimmed vdofs (accounts for new
  // boundary elements too since no new dofs are added). See the MFEM trimmer miniapp for
  // reference.
  new_mesh->FinalizeTopology();
  new_mesh->Finalize();
  new_mesh->RemoveUnusedVertices();
  if (orig_mesh->GetNodes())
  {
    const mfem::GridFunction *nodes = orig_mesh->GetNodes();
    const mfem::FiniteElementSpace *fespace = nodes->FESpace();

    mfem::Ordering::Type ordering = fespace->GetOrdering();
    int order = fespace->GetMaxElementOrder();
    int sdim = orig_mesh->SpaceDimension();
    bool discont =
        dynamic_cast<const mfem::L2_FECollection *>(fespace->FEColl()) != nullptr;

    new_mesh->SetCurvature(order, discont, sdim, ordering);
    mfem::GridFunction *new_nodes = new_mesh->GetNodes();
    const mfem::FiniteElementSpace *new_fespace = new_nodes->FESpace();

    // The element loop works because we know the mapping from old_mesh to new_mesh element
    // indices from the insertion order.
    mfem::Array<int> vdofs, new_vdofs;
    mfem::Vector loc_vec;
    int te = 0;
    for (int e = 0; e < orig_mesh->GetNE(); e++)
    {
      if (!elem_delete[e])
      {
        fespace->GetElementVDofs(e, vdofs);
        nodes->GetSubVector(vdofs, loc_vec);
        new_fespace->GetElementVDofs(te, new_vdofs);
        new_nodes->SetSubVector(new_vdofs, loc_vec);
        te++;
      }
    }
  }
  orig_mesh = std::move(new_mesh);
  return new_attr_map;
}

std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm comm,
                                              std::unique_ptr<mfem::Mesh> &smesh,
                                              std::unique_ptr<int[]> &partitioning)
{
  // Take a serial mesh and partitioning on the root process and construct the global
  // parallel mesh. For now, prefer the MPI-based version.
#if 0
  {
    // Write each processor's component to file.
    std::string tmp = iodata.problem.output;
    if (tmp.back() != '/')
    {
      tmp += '/';
    }
    tmp += "tmp/";
    int width = 1 + static_cast<int>(std::log10(Mpi::Size(comm) - 1));
    if (Mpi::Root(comm))
    {
      if (!std::filesystem::exists(tmp))
      {
        std::filesystem::create_directories(tmp);
      }
      mfem::MeshPartitioner partitioner(*smesh, Mpi::Size(comm), partitioning.get());
      for (int i = 0; i < Mpi::Size(comm); i++)
      {
        mfem::MeshPart part;
        partitioner.ExtractPart(i, part);
        std::string pfile = mfem::MakeParFilename(tmp + "part.", i, ".mesh", width);
        std::ofstream fo(pfile);
        // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
        // fo << std::fixed;
        fo << std::scientific;
        fo.precision(MSH_FLT_PRECISION);
        part.Print(fo);
      }
    }

    // Each process loads its own partitioned mesh file and constructs the parallel mesh.
    std::string pfile =
        mfem::MakeParFilename(tmp + "part.", Mpi::Rank(comm), ".mesh", width);
    int exists = 0;
    while (!exists)  // Wait for root to finish writing all files
    {
      exists = std::filesystem::exists(pfile);
      Mpi::GlobalMax(1, &exists, comm);
    }
    std::ifstream fi(pfile);
    // mfem::ifgzstream fi(pfile);
    if (!fi.good())
    {
      MFEM_ABORT("Unable to open partitioned mesh file \"" << pfile << "\"!");
    }
    auto pmesh = std::make_unique<mfem::ParMesh>(comm, fi);
    Mpi::Barrier(comm);
    if (Mpi::Root(comm))
    {
      std::filesystem::remove_all(tmp);  // Remove the temporary directory
    }
    return pmesh;
  }
#endif
  {
    // Send each processor's component as a byte string.
    std::vector<std::string> so;
    if (Mpi::Root(comm))
    {
      mfem::MeshPartitioner partitioner(*smesh, Mpi::Size(comm), partitioning.get());
      so.reserve(Mpi::Size(comm));
      for (int i = 0; i < Mpi::Size(comm); i++)
      {
        mfem::MeshPart part;
        partitioner.ExtractPart(i, part);
        std::ostringstream fo(std::stringstream::out);
        // fo << std::fixed;
        fo << std::scientific;
        fo.precision(MSH_FLT_PRECISION);
        part.Print(fo);
        so.push_back(fo.str());
        // so.push_back((i > 0) ? zlib::CompressString(fo.str()) : fo.str());
      }
    }

    // Scatter the partitioned mesh files and generate the parallel mesh.
    if (Mpi::Root(comm))
    {
      std::vector<MPI_Request> send_requests(Mpi::Size(comm) - 1, MPI_REQUEST_NULL);
      for (int i = 1; i < Mpi::Size(comm); i++)
      {
        int ilen = static_cast<int>(so[i].length());
        MFEM_VERIFY(so[i].length() == (std::size_t)ilen,
                    "Overflow error distributing parallel mesh!");
        MPI_Isend(so[i].c_str(), ilen, MPI_CHAR, i, i, comm, &send_requests[i - 1]);
      }
      std::istringstream fi(so[0]);  // This is never compressed
      auto pmesh = std::make_unique<mfem::ParMesh>(comm, fi);
      MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(),
                  MPI_STATUSES_IGNORE);
      return pmesh;
    }
    int rlen;
    MPI_Status status;
    MPI_Probe(0, Mpi::Rank(comm), comm, &status);
    MPI_Get_count(&status, MPI_CHAR, &rlen);

    std::string si;
    si.resize(rlen);
    MPI_Recv(si.data(), rlen, MPI_CHAR, 0, Mpi::Rank(comm), comm, MPI_STATUS_IGNORE);
    std::istringstream fi(si);
    // std::istringstream fi(zlib::DecompressString(si));
    return std::make_unique<mfem::ParMesh>(comm, fi);
  }
}

void GetUsedAttributeMarkers(const IoData &iodata, int n_mat, int n_bdr,
                             mfem::Array<int> &mat_marker, mfem::Array<int> &bdr_marker)
{
  mfem::Array<int> mat_attr, bdr_attr;
  mat_attr.Reserve(static_cast<int>(iodata.domains.attributes.size()));
  for (auto attr : iodata.domains.attributes)
  {
    mat_attr.Append(attr);
  }
  bdr_attr.Reserve(static_cast<int>(iodata.boundaries.attributes.size()));
  for (auto attr : iodata.boundaries.attributes)
  {
    bdr_attr.Append(attr);
  }
  mesh::AttrToMarker(n_mat, mat_attr, mat_marker);
  mesh::AttrToMarker(n_bdr, bdr_attr, bdr_marker);
}

}  // namespace

}  // namespace palace
