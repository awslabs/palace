// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/timer.hpp"

namespace palace
{

using Vector3dMap = Eigen::Map<Eigen::Vector3d>;
using CVector3dMap = Eigen::Map<const Eigen::Vector3d>;

namespace
{

// Floating point precision for mesh IO. This precision is important, make sure nothing is
// lost!
constexpr auto MSH_FLT_PRECISION = std::numeric_limits<double>::max_digits10;

// Load the serial mesh from disk.
std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &, bool);

// Optionally reorder mesh elements based on MFEM's internal reordeing tools for improved
// cache usage.
void ReorderMesh(mfem::Mesh &);

// Generate element-based mesh partitioning, using either a provided file or METIS.
std::unique_ptr<int[]> GetMeshPartitioning(mfem::Mesh &, int, const std::string & = "");

// Cleanup the provided serial mesh by removing unnecessary domain and elements, adding
// boundary elements for material interfaces and exterior boundaries, and adding boundary
// elements for subdomain interfaces.
std::map<int, std::array<int, 2>> CheckMesh(mfem::Mesh &, const std::unique_ptr<int[]> &,
                                            const IoData &, bool, bool, bool);

// Given a serial mesh on the root processor and element partitioning, create a parallel
// mesh over the given communicator.
std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm, const std::unique_ptr<mfem::Mesh> &,
                                              const std::unique_ptr<int[]> & = nullptr,
                                              const std::string & = "");

// Get list of domain and boundary attribute markers used in configuration file for mesh
// cleaning.
void GetUsedAttributeMarkers(const IoData &, int, int, mfem::Array<int> &,
                             mfem::Array<int> &);

// Rebalance a conformal mesh across processor ranks, using the MeshPartitioner. Gathers the
// mesh onto the root rank before scattering the partitioned mesh.
void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &, double, const std::string &);

struct ElementTypeInfo
{
  bool has_simplices;
  bool has_tensors;
  bool has_wedges;
  bool has_pyramids;
};

// Simplified helper for describing the element types in a mesh.
ElementTypeInfo CheckElements(mfem::Mesh &mesh)
{
  // MeshGenerator is reduced over the communicator. This checks for geometries on any
  // processor.
  auto meshgen = mesh.MeshGenerator();
  return {bool(meshgen & 1), bool(meshgen & 2), bool(meshgen & 4), bool(meshgen & 8)};
}

}  // namespace

namespace mesh
{

std::unique_ptr<mfem::ParMesh> ReadMesh(MPI_Comm comm, const IoData &iodata, bool reorder,
                                        bool clean, bool add_bdr, bool unassembled)
{
  // If possible on root, read the serial mesh (converting format if necessary), and do all
  // necessary serial preprocessing. When finished, distribute the mesh to all processes.
  // Count disk I/O time separately for the mesh read from file.

  // If not adapting, or performing conformal adaptation, can use the mesh partitioner.
  std::unique_ptr<mfem::Mesh> smesh;
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = refinement.max_it > 0;
  const bool use_mesh_partitioner = !use_amr || !refinement.nonconformal;
  {
    BlockTimer bt(Timer::IO);
    if (Mpi::Root(comm) || !use_mesh_partitioner)
    {
      // Optionally reorder elements (and vertices) based on spatial location after loading
      // the serial mesh.
      smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature);
      if (reorder)
      {
        smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature);
        if (reorder)
        {
          ReorderMesh(*smesh);
        }
      }
      MFEM_VERIFY(!(smesh->Nonconforming() && use_mesh_partitioner),
                  "Cannot use mesh partitioner on a nonconforming mesh");
    }
    Mpi::Barrier(comm);
  }

  std::unique_ptr<int[]> partitioning;
  if (Mpi::Root(comm) || !use_mesh_partitioner)
  {
    // Check the the AMR specification and the mesh elements are compatible.
    const auto element_types = CheckElements(*smesh);
    MFEM_VERIFY(!use_amr || !element_types.has_tensors || refinement.nonconformal,
                "If there are tensor elements, AMR must be nonconformal");
    MFEM_VERIFY(!use_amr || !element_types.has_pyramids || refinement.nonconformal,
                "If there are pyramid elements, AMR must be nonconformal");
    MFEM_VERIFY(!use_amr || !element_types.has_wedges || refinement.nonconformal,
                "If there are wedge elements, AMR must be nonconformal");

    // Generate the mesh partitioning.
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), iodata.model.partition);

    // Clean up unused domain elements from the mesh, add new boundary elements for material
    // interfaces if not present, and optionally (when running unassembled) add subdomain
    // interface boundary elements. Can only clean up conforming meshes, assumes that any
    // nonconformal mesh was generated by adaptation and thus does not need checking.
    if (smesh->Conforming())
    {
      static_cast<void>(
          CheckMesh(*smesh, partitioning, iodata, clean, add_bdr, unassembled));
    }
  }

  std::unique_ptr<mfem::ParMesh> pmesh;
  if (use_mesh_partitioner)
  {
    pmesh = DistributeMesh(comm, smesh, partitioning, iodata.problem.output);
  }
  else
  {
    if (refinement.nonconformal && use_amr)
    {
      smesh->EnsureNCMesh(true);
    }
    pmesh = std::make_unique<mfem::ParMesh>(comm, *smesh, partitioning.get());
  }

  if constexpr (false)
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
    int width = 1 + static_cast<int>(std::log10(Mpi::Size(comm) - 1));
    std::unique_ptr<mfem::Mesh> gsmesh =
        LoadMesh(iodata.model.mesh, iodata.model.remove_curvature);
    std::unique_ptr<int[]> gpartitioning = GetMeshPartitioning(*gsmesh, Mpi::Size(comm));
    mfem::ParMesh gpmesh(comm, *gsmesh, gpartitioning.get(), 0);
    {
      std::string pfile =
          mfem::MakeParFilename(tmp + "part.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
      fo.precision(MSH_FLT_PRECISION);
      gpmesh.ParPrint(fo);
    }
    {
      std::string pfile =
          mfem::MakeParFilename(tmp + "final.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
      fo.precision(MSH_FLT_PRECISION);
      pmesh->ParPrint(fo);
    }
  }

  return pmesh;
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
    PalacePragmaDiagnosticPush
    PalacePragmaDiagnosticDisableDeprecated
    mesh[0]->ReorientTetMesh();
    PalacePragmaDiagnosticPop
  }

  // Uniformly refine the mesh further in parallel, saving the level meshes for geometric
  // coarsening later on if desired.
  for (int l = 0; l < uniform_ref_levels; l++)
  {
    if (mesh.capacity() > 1)
    {
      mesh.emplace_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
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
  const bool use_nodes = (mesh.back()->GetNodes() != nullptr);
  const int ref = use_nodes ? mesh.back()->GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  const int dim = mesh.back()->SpaceDimension();
  int region_ref_level = 0;
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
      mesh.emplace_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
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
    PalacePragmaDiagnosticPush
    PalacePragmaDiagnosticDisableDeprecated
    mesh[0]->ReorientTetMesh();
    PalacePragmaDiagnosticPop
  }

  // Print some mesh information.
  mfem::Vector bbmin, bbmax;
  mesh[0]->GetBoundingBox(bbmin, bbmax);
  const double Lc = iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0);
  Mpi::Print(mesh[0]->GetComm(), "\nMesh curvature order: {}\nMesh bounding box:\n",
             mesh[0]->GetNodes()
                 ? std::to_string(mesh[0]->GetNodes()->FESpace()->GetMaxElementOrder())
                 : "None");
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

namespace
{

void ScaleMesh(mfem::Mesh &mesh, double L)
{
  for (int i = 0; i < mesh.GetNV(); i++)
  {
    double *v = mesh.GetVertex(i);
    std::transform(v, v + mesh.SpaceDimension(), v, [L](double val) { return val * L; });
  }
  if (mesh.GetNodes())
  {
    *mesh.GetNodes() *= L;
  }
}

}  // namespace

void DimensionalizeMesh(mfem::Mesh &mesh, double L)
{
  ScaleMesh(mesh, L);
}

void NondimensionalizeMesh(mfem::Mesh &mesh, double L)
{
  ScaleMesh(mesh, 1.0 / L);
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

void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr, mfem::Vector &min,
                               mfem::Vector &max)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  GetAxisAlignedBoundingBox(mesh, marker, bdr, min, max);
}

void GetAxisAlignedBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                               bool bdr, mfem::Vector &min, mfem::Vector &max)
{
  int dim = mesh.SpaceDimension();
  min.SetSize(dim);
  max.SetSize(dim);
  for (int d = 0; d < dim; d++)
  {
    min(d) = mfem::infinity();
    max(d) = -mfem::infinity();
  }
  if (!mesh.GetNodes())
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
    const int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
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

double BoundingBox::Area() const
{
  return 4.0 *
         CVector3dMap(normals[0].data()).cross(CVector3dMap(normals[1].data())).norm();
}

double BoundingBox::Volume() const
{
  return planar ? 0.0 : 2.0 * CVector3dMap(normals[2].data()).norm() * Area();
}

std::array<double, 3> BoundingBox::Lengths() const
{
  return {2.0 * CVector3dMap(normals[0].data()).norm(),
          2.0 * CVector3dMap(normals[1].data()).norm(),
          2.0 * CVector3dMap(normals[2].data()).norm()};
}

std::array<double, 3> BoundingBox::Deviation(const std::array<double, 3> &direction) const
{
  const auto eig_dir = CVector3dMap(direction.data());
  std::array<double, 3> deviation_deg;
  for (std::size_t i = 0; i < 3; i++)
  {
    deviation_deg[i] =
        std::acos(std::min(1.0, std::abs(eig_dir.normalized().dot(
                                    CVector3dMap(normals[i].data()).normalized())))) *
        (180.0 / M_PI);
  }
  return deviation_deg;
}

namespace
{

// Compute a lexicographic comparison of Eigen Vector3d.
bool EigenLE(const Eigen::Vector3d &x, const Eigen::Vector3d &y)
{
  return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
};

// Helper for collecting a point cloud from a mesh, used in calculating bounding boxes and
// bounding balls. Returns the dominant rank, for which the vertices argument will be
// filled, while all other ranks will have an empty vector. Vertices are de-duplicated to a
// certain floating point precision.
int CollectPointCloudOnRoot(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                            std::vector<Eigen::Vector3d> &vertices)
{
  std::set<int> vertex_indices;
  if (!mesh.GetNodes())
  {
    // Linear mesh, work with element vertices directly.
    mfem::Array<int> v;
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
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
    if (bdr)
    {
      for (int i = 0; i < mesh.GetNBE(); i++)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mfem::ElementTransformation *T = mesh.GetBdrElementTransformation(i);
        T->Transform(
            mfem::GlobGeometryRefiner.Refine(mesh.GetBdrElementGeometry(i), ref)->RefPts,
            pointmat);
        for (int j = 0; j < pointmat.Width(); j++)
        {
          vertices.push_back({pointmat(0, j), pointmat(1, j), pointmat(2, j)});
        }
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
        T->Transform(
            mfem::GlobGeometryRefiner.Refine(mesh.GetElementGeometry(i), ref)->RefPts,
            pointmat);
        for (int j = 0; j < pointmat.Width(); j++)
        {
          vertices.push_back({pointmat(0, j), pointmat(1, j), pointmat(2, j)});
        }
      }
    }
  }

  // dominant_rank will perform the calculation.
  MPI_Comm comm = mesh.GetComm();
  const auto num_vertices = int(vertices.size());
  const int dominant_rank = [&]()
  {
    int vert = num_vertices, rank = Mpi::Rank(comm);
    Mpi::GlobalMaxLoc(1, &vert, &rank, comm);
    return rank;
  }();
  std::vector<int> recv_counts(Mpi::Size(comm)), displacements;
  std::vector<Eigen::Vector3d> collected_vertices;
  MPI_Gather(&num_vertices, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, dominant_rank,
             comm);
  if (dominant_rank == Mpi::Rank(comm))
  {
    // First displacement is zero, then after is the partial sum of recv_counts.
    displacements.resize(Mpi::Size(comm));
    displacements[0] = 0;
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);

    // Add on slots at the end of vertices for the incoming data.
    collected_vertices.resize(std::accumulate(recv_counts.begin(), recv_counts.end(), 0));

    // MPI transfer will be done with MPI_DOUBLE, so duplicate all these values.
    for (auto &x : displacements)
    {
      x *= 3;
    }
    for (auto &x : recv_counts)
    {
      x *= 3;
    }
  }

  // Gather the data to the dominant rank.
  static_assert(sizeof(Eigen::Vector3d) == 3 * sizeof(double));
  MPI_Gatherv(vertices.data(), 3 * num_vertices, MPI_DOUBLE, collected_vertices.data(),
              recv_counts.data(), displacements.data(), MPI_DOUBLE, dominant_rank, comm);

  // Deduplicate vertices. Given floating point precision, need a tolerance.
  if (dominant_rank == Mpi::Rank(comm))
  {
    auto vertex_equality = [](const auto &x, const auto &y)
    {
      constexpr double tolerance = 10.0 * std::numeric_limits<double>::epsilon();
      return std::abs(x[0] - y[0]) < tolerance && std::abs(x[1] - y[1]) < tolerance &&
             std::abs(x[2] - y[2]) < tolerance;
    };
    vertices = std::move(collected_vertices);
    std::sort(vertices.begin(), vertices.end(), EigenLE);
    vertices.erase(std::unique(vertices.begin(), vertices.end(), vertex_equality),
                   vertices.end());
  }
  else
  {
    vertices.clear();
  }

  return dominant_rank;
}

// Calculates a bounding box from a point cloud, result is broadcast across all processes.
BoundingBox BoundingBoxFromPointCloud(MPI_Comm comm,
                                      const std::vector<Eigen::Vector3d> &vertices,
                                      int dominant_rank)
{
  BoundingBox box;
  if (dominant_rank == Mpi::Rank(comm))
  {
    // Pick a candidate 000 vertex using lexicographic sort. This can be vulnerable to
    // floating point precision if the box is axis aligned, but not floating point exact.
    // Pick candidate 111 as the furthest from this candidate, then reassign 000 as the
    // furthest from 111. Such a pair has to form the diagonal for a point cloud defining a
    // box. Verify that p_111 is also the maximum distance from p_000 -> a diagonal is
    // found.
    MFEM_VERIFY(vertices.size() >= 4,
                "A bounding box requires a minimum of four vertices for this algorithm!");
    auto p_000 = std::min_element(vertices.begin(), vertices.end(), EigenLE);
    auto DistFromP_000 = [&p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
    { return (x - *p_000).norm() < (y - *p_000).norm(); };
    auto p_111 = std::max_element(vertices.begin(), vertices.end(), DistFromP_000);
    auto DistFromP_111 = [&p_111](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
    { return (x - *p_111).norm() < (y - *p_111).norm(); };
    p_000 = std::max_element(vertices.begin(), vertices.end(), DistFromP_111);
    MFEM_VERIFY(std::max_element(vertices.begin(), vertices.end(), DistFromP_000) == p_111,
                "p_000 and p_111 must be mutually opposing points!");

    // Define a diagonal of the ASSUMED cuboid bounding box. Store references as this is
    // useful for checking pointers later.
    const auto &v_000 = *p_000;
    const auto &v_111 = *p_111;
    MFEM_VERIFY(&v_000 != &v_111, "Minimum and maximum extents cannot be identical!");
    const auto origin = v_000;
    const Eigen::Vector3d n_1 = (v_111 - v_000).normalized();

    // Compute the distance from the normal axis. Note: everything has been oriented
    // relative to v_000 == (0,0,0).
    auto PerpendicularDistance = [&n_1, &origin](const Eigen::Vector3d &v)
    { return ((v - origin) - (v - origin).dot(n_1) * n_1).norm(); };

    // Find the vertex furthest from the diagonal axis. We cannot know yet if this defines
    // (001) or (011).
    const auto &t_0 =
        *std::max_element(vertices.begin(), vertices.end(),
                          [PerpendicularDistance](const auto &x, const auto &y)
                          { return PerpendicularDistance(x) < PerpendicularDistance(y); });
    MFEM_VERIFY(&t_0 != &v_000, "Vertices are degenerate!");
    MFEM_VERIFY(&t_0 != &v_111, "Vertices are degenerate!");

    // Use the discovered vertex to define a second direction and thus a plane.
    const Eigen::Vector3d n_2 =
        ((t_0 - origin) - (t_0 - origin).dot(n_1) * n_1).normalized();

    // n_1 and n_2 now define a planar coordinate system intersecting the main diagonal, and
    // two opposite edges of the cuboid. Now look for a component that maximizes distance
    // from the planar system: complete the axes with a cross, then use a dot product to
    // pick the greatest deviation.
    auto OutOfPlaneDistance = [&n_1, &n_2, &origin](const Eigen::Vector3d &v)
    {
      return ((v - origin) - (v - origin).dot(n_1) * n_1 - (v - origin).dot(n_2) * n_2)
          .norm();
    };

    // Collect the furthest point from the plane.
    auto max_distance = OutOfPlaneDistance(
        *std::max_element(vertices.begin(), vertices.end(),
                          [OutOfPlaneDistance](const auto &x, const auto &y)
                          { return OutOfPlaneDistance(x) < OutOfPlaneDistance(y); }));

    constexpr double rel_tol = 1e-6;
    box.planar = max_distance < (rel_tol * (v_111 - v_000).norm());

    // Given numerical tolerance, collect other points with an almost matching distance.
    std::vector<Eigen::Vector3d> vertices_out_of_plane;
    const double cooincident_tolerance = rel_tol * max_distance;
    std::copy_if(
        vertices.begin(), vertices.end(), std::back_inserter(vertices_out_of_plane),
        [OutOfPlaneDistance, cooincident_tolerance, max_distance](const auto &v)
        { return std::abs(OutOfPlaneDistance(v) - max_distance) < cooincident_tolerance; });

    // Given candidates t_0 and t_1, the closer to origin defines v_001.
    const auto &t_1 = box.planar
                          ? t_0
                          : *std::min_element(vertices_out_of_plane.begin(),
                                              vertices_out_of_plane.end(), DistFromP_000);
    const bool t_0_gt_t_1 =
        (t_0 - origin).norm() > (t_1 - origin).norm();  // If planar t_1 == t_0
    const auto &v_001 = t_0_gt_t_1 ? t_1 : t_0;
    const auto &v_011 = box.planar ? v_111 : (t_0_gt_t_1 ? t_0 : t_1);

    // Compute the center as halfway along the main diagonal.
    Vector3dMap(box.center.data()) = 0.5 * (v_000 + v_111);

    // The length in each direction is given by traversing the edges of the cuboid in turn.
    Vector3dMap(box.normals[0].data()) = 0.5 * (v_001 - v_000);
    Vector3dMap(box.normals[1].data()) = 0.5 * (v_011 - v_001);
    Vector3dMap(box.normals[2].data()) = 0.5 * (v_111 - v_011);

    // Make sure the longest dimension comes first.
    std::sort(box.normals.begin(), box.normals.end(),
              [](const auto &x, const auto &y)
              { return CVector3dMap(x.data()).norm() > CVector3dMap(y.data()).norm(); });
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, box.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, box.normals.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &box.planar, dominant_rank, comm);

  return box;
}

// Calculates a bounding ball from a point cloud, result is broadcast across all processes.
BoundingBall BoundingBallFromPointCloud(MPI_Comm comm,
                                        const std::vector<Eigen::Vector3d> &vertices,
                                        int dominant_rank)
{
  BoundingBall ball;
  if (dominant_rank == Mpi::Rank(comm))
  {
    // Pick a candidate 000 vertex using lexicographic sort. This can be vulnerable to
    // floating point precision if there is no directly opposed vertex.
    // Pick candidate 111 as the furthest from this candidate, then reassign 000 as the
    // furthest from 111. Such a pair has to form the diagonal for a point cloud defining a
    // ball. Verify that p_111 is also the maximum distance from p_000 -> a diagonal is
    // found.
    MFEM_VERIFY(vertices.size() >= 3,
                "A bounding ball requires a minimum of three vertices for this algorithm!");
    auto p_000 = std::min_element(vertices.begin(), vertices.end(), EigenLE);
    auto DistFromP_000 = [&p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
    { return (x - *p_000).norm() < (y - *p_000).norm(); };
    auto p_111 = std::max_element(vertices.begin(), vertices.end(), DistFromP_000);
    auto DistFromP_111 = [&p_111](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
    { return (x - *p_111).norm() < (y - *p_111).norm(); };
    p_000 = std::max_element(vertices.begin(), vertices.end(), DistFromP_111);
    MFEM_VERIFY(std::max_element(vertices.begin(), vertices.end(), DistFromP_000) == p_111,
                "p_000 and p_111 must be mutually opposing points!");

    const auto &min = *p_000;
    const auto &max = *p_111;
    Eigen::Vector3d delta = max - min;
    ball.radius = 0.5 * delta.norm();
    Vector3dMap(ball.center.data()) = 0.5 * (min + max);

    // Project onto this candidate diameter, and pick a vertex furthest away. Check that
    // this resulting distance is less than or equal to the radius, and use the resulting
    // direction to compute another in plane vector. Assumes all delta are normalized, and
    // applies a common origin as part of the projection.
    auto PerpendicularDistance = [min](const std::initializer_list<Eigen::Vector3d> &deltas,
                                       const Eigen::Vector3d &vin)
    {
      Eigen::Vector3d v = vin - min;
      for (const auto &d : deltas)
      {
        v -= d.dot(v) * d;
      }
      return v.norm();
    };

    delta.normalize();
    const auto &perp = *std::max_element(
        vertices.begin(), vertices.end(),
        [&delta, PerpendicularDistance](const auto &x, const auto &y)
        { return PerpendicularDistance({delta}, x) < PerpendicularDistance({delta}, y); });
    constexpr double rel_tol = 1.0e-6;
    MFEM_VERIFY(std::abs(PerpendicularDistance({delta}, perp) - ball.radius) <=
                    rel_tol * ball.radius,
                "Furthest point perpendicular must be on the exterior of the ball: "
                    << PerpendicularDistance({delta}, perp) << " vs. " << ball.radius
                    << "!");

    // Compute a perpendicular to the circle using the cross product.
    const Eigen::Vector3d n_radial = (perp - CVector3dMap(ball.center.data())).normalized();
    Vector3dMap(ball.planar_normal.data()) = delta.cross(n_radial).normalized();

    // Compute the point furthest out of the plane discovered. If below tolerance, this
    // means the ball is 2D.
    const auto &out_of_plane = *std::max_element(
        vertices.begin(), vertices.end(),
        [&delta, &n_radial, PerpendicularDistance](const auto &x, const auto &y)
        {
          return PerpendicularDistance({delta, n_radial}, x) <
                 PerpendicularDistance({delta, n_radial}, y);
        });

    ball.planar =
        PerpendicularDistance({delta, n_radial}, out_of_plane) / ball.radius < rel_tol;
    if (!ball.planar)
    {
      // The points are not functionally coplanar, zero out the normal.
      MFEM_VERIFY(std::abs(PerpendicularDistance({delta}, perp) - ball.radius) <=
                      rel_tol * ball.radius,
                  "Furthest point perpendicular must be on the exterior of the sphere!");
      Vector3dMap(ball.planar_normal.data()) *= 0;
    }
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, ball.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3, ball.planar_normal.data(), dominant_rank, comm);
  Mpi::Broadcast(1, &ball.radius, dominant_rank, comm);
  Mpi::Broadcast(1, &ball.planar, dominant_rank, comm);

  return ball;
}

double LengthFromPointCloud(MPI_Comm comm, const std::vector<Eigen::Vector3d> &vertices,
                            int dominant_rank, const std::array<double, 3> &dir)
{
  double length;
  if (dominant_rank == Mpi::Rank(comm))
  {
    CVector3dMap direction(dir.data());

    auto Dot = [&](const auto &x, const auto &y)
    { return direction.dot(x) < direction.dot(y); };
    auto p_min = std::min_element(vertices.begin(), vertices.end(), Dot);
    auto p_max = std::max_element(vertices.begin(), vertices.end(), Dot);

    length = (*p_max - *p_min).dot(direction.normalized());
  }
  Mpi::Broadcast(1, &length, dominant_rank, comm);
  return length;
}

}  // namespace

double GetProjectedLength(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr,
                          const std::array<double, 3> &dir)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  return LengthFromPointCloud(mesh.GetComm(), vertices, dominant_rank, dir);
}

double GetProjectedLength(mfem::ParMesh &mesh, int attr, bool bdr,
                          const std::array<double, 3> &dir)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetProjectedLength(mesh, marker, bdr, dir);
}

BoundingBox GetBoundingBox(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  return BoundingBoxFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
}

BoundingBox GetBoundingBox(mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetBoundingBox(mesh, marker, bdr);
}

BoundingBall GetBoundingBall(mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  return BoundingBallFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
}

BoundingBall GetBoundingBall(mfem::ParMesh &mesh, int attr, bool bdr)
{
  mfem::Array<int> marker(bdr ? mesh.bdr_attributes.Max() : mesh.attributes.Max());
  marker = 0;
  marker[attr - 1] = 1;
  return GetBoundingBall(mesh, marker, bdr);
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

double RebalanceMesh(const IoData &iodata, std::unique_ptr<mfem::ParMesh> &mesh, double tol)
{
  BlockTimer bt0(Timer::REBALANCE);
  const bool save_adapt_mesh = iodata.model.refinement.save_adapt_mesh;
  std::string serial_mesh_file;
  if (save_adapt_mesh)
  {
    serial_mesh_file = iodata.problem.output;
    if (serial_mesh_file.back() != '/')
    {
      serial_mesh_file += '/';
    }
    serial_mesh_file += "serial.mesh";
  }

  MPI_Comm comm = mesh->GetComm();
  if (Mpi::Size(comm) == 1)
  {
    if (save_adapt_mesh)
    {
      BlockTimer bt1(Timer::IO);
      std::ofstream fo(serial_mesh_file);
      fo.precision(MSH_FLT_PRECISION);
      mesh::DimensionalizeMesh(*mesh, iodata.GetLengthScale());
      mesh->mfem::Mesh::Print(fo);
      mesh::NondimensionalizeMesh(*mesh, iodata.GetLengthScale());
    }
    return 1.0;
  }

  // If there is more than one processor, may perform rebalancing.
  mesh->ExchangeFaceNbrData();
  int min_elem, max_elem;
  min_elem = max_elem = mesh->GetNE();
  Mpi::GlobalMin(1, &min_elem, comm);
  Mpi::GlobalMax(1, &max_elem, comm);
  const double ratio = double(max_elem) / min_elem;
  if constexpr (false)
  {
    Mpi::Print("Rebalancing: max/min elements per processor = {:d}/{:d} (ratio = {:.3e}, "
               "tol = {:.3e})\n",
               max_elem, min_elem, ratio, tol);
  }
  if (ratio > tol)
  {
    if (mesh->Nonconforming() && save_adapt_mesh)
    {
      // Do not need to duplicate the mesh, as rebalancing will undo this.
      mfem::Array<int> serial_partition(mesh->GetNE());
      serial_partition = 0;
      mesh->Rebalance(serial_partition);
      BlockTimer bt1(Timer::IO);
      if (Mpi::Root(comm))
      {
        std::ofstream fo(serial_mesh_file);
        fo.precision(MSH_FLT_PRECISION);
        mesh::DimensionalizeMesh(*mesh, iodata.GetLengthScale());
        mesh->Mesh::Print(fo);
        mesh::NondimensionalizeMesh(*mesh, iodata.GetLengthScale());
      }
      Mpi::Barrier(comm);
    }
    if (mesh->Nonconforming())
    {
      mesh->Rebalance();
    }
    else
    {
      // Without access to a refinement tree, partitioning must be done on the root
      // processor and then redistributed.
      RebalanceConformalMesh(mesh, iodata.GetLengthScale(), serial_mesh_file);
    }
  }
  else if (save_adapt_mesh)
  {
    // Given no rebalancing will be done, need to handle the serial write more carefully.
    // This requires creating a separate serial mesh.
    if (mesh->Nonconforming())
    {
      mfem::ParMesh smesh(*mesh);
      mfem::Array<int> serial_partition(mesh->GetNE());
      serial_partition = 0;
      smesh.Rebalance(serial_partition);
      BlockTimer bt1(Timer::IO);
      if (Mpi::Root(comm))
      {
        std::ofstream fo(serial_mesh_file);
        fo.precision(MSH_FLT_PRECISION);
        mesh::DimensionalizeMesh(smesh, iodata.GetLengthScale());
        smesh.Mesh::Print(fo);  // Do not need to nondimensionalize the temporary mesh
      }
      Mpi::Barrier(comm);
    }
    else
    {
      auto smesh = std::make_unique<mfem::Mesh>(mesh->GetSerialMesh(0));
      BlockTimer bt1(Timer::IO);
      if (Mpi::Rank(comm) == 0)
      {
        std::ofstream fo(serial_mesh_file);
        fo.precision(MSH_FLT_PRECISION);
        mesh::DimensionalizeMesh(*smesh, iodata.GetLengthScale());
        smesh->Print(fo);  // Do not need to nondimensionalize the temporary mesh
      }
      Mpi::Barrier(comm);
    }
  }
  mesh->ExchangeFaceNbrData();
  return ratio;
}

}  // namespace mesh

namespace
{

std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &path, bool remove_curvature)
{
  // Read the (serial) mesh from the given mesh file. Handle preparation for refinement and
  // orientations here to avoid possible reorientations and reordering later on. MFEM
  // supports a native mesh format (.mesh), VTK/VTU, Gmsh, as well as some others. We use
  // built-in converters for the types we know, otherwise rely on MFEM to do the conversion
  // or error out if not supported.
  constexpr bool generate_edges = true, refine = true, fix_orientation = true;
  std::unique_ptr<mfem::Mesh> mesh;
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

    mesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
  }
  else
  {
    // Otherwise, just rely on MFEM load the mesh.
    std::ifstream fi(path);
    if (!fi.good())
    {
      MFEM_ABORT("Unable to open mesh file \"" << path << "\"!");
    }
    mesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
  }
  if (remove_curvature)
  {
    if (mesh->GetNodes())
    {
      mfem::GridFunction *nodes = nullptr;
      int own_nodes = true;
      mesh->SwapNodes(nodes, own_nodes);
      if (own_nodes)
      {
        delete nodes;
      }
    }
  }
  else
  {
    mesh->EnsureNodes();
  }
  return mesh;
}

void ReorderMesh(mfem::Mesh &mesh)
{
  mfem::Array<int> ordering;

  if constexpr (false)
  {
    // Gecko reordering.
    mfem::Array<int> tentative;
    int outer = 3, inner = 3, window = 4, period = 2;
    double best_cost = mfem::infinity();
    for (int i = 0; i < outer; i++)
    {
      int seed = i + 1;
      double cost =
          mesh.GetGeckoElementOrdering(tentative, inner, window, period, seed, true);
      if (cost < best_cost)
      {
        ordering = tentative;
        best_cost = cost;
      }
    }
    Mpi::Print("Final cost: {:e}\n", best_cost);
  }

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
  int ne, np;
  std::ifstream part_ifs(partition);
  part_ifs.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
  part_ifs >> ne;
  if (ne != mesh.GetNE())
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

std::map<int, std::array<int, 2>> CheckMesh(mfem::Mesh &orig_mesh,
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
  MFEM_VERIFY(orig_mesh.Dimension() == 3 && !orig_mesh.Nonconforming(),
              "Nonconforming or 2D meshes have not been tested yet!");
  MFEM_VERIFY(dynamic_cast<mfem::ParMesh *>(&orig_mesh) == nullptr,
              "This function does not work for ParMesh");
  mfem::Array<int> mat_marker, bdr_marker;
  GetUsedAttributeMarkers(
      iodata, orig_mesh.attributes.Size() ? orig_mesh.attributes.Max() : 0,
      orig_mesh.bdr_attributes.Size() ? orig_mesh.bdr_attributes.Max() : 0, mat_marker,
      bdr_marker);
  bool warn = false;
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    int attr = orig_mesh.GetBdrAttribute(be);
    if (!bdr_marker[attr - 1])
    {
      int f, o, e1, e2;
      orig_mesh.GetBdrElementFace(be, &f, &o);
      orig_mesh.GetFaceElements(f, &e1, &e2);
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
  int new_ne = orig_mesh.GetNE();
  int new_nbdr = orig_mesh.GetNBE();
  mfem::Array<bool> elem_delete, bdr_delete;
  mfem::Array<int> orig_bdr_faces, add_bdr_faces;
  elem_delete.SetSize(orig_mesh.GetNE(), false);
  bdr_delete.SetSize(orig_mesh.GetNBE(), false);
  orig_bdr_faces.SetSize(orig_mesh.GetNumFaces(), -1);
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    int f, o;
    orig_mesh.GetBdrElementFace(be, &f, &o);
    MFEM_VERIFY(orig_bdr_faces[f] < 0,
                "Mesh should not define boundary elements multiple times!");
    orig_bdr_faces[f] = be;
  }
  if (add_bdr || add_subdomain)
  {
    add_bdr_faces.SetSize(orig_mesh.GetNumFaces(), -1);
  }

  if (clean_elem)
  {
    // Delete domain and boundary elements which have no associated material or BC attribute
    // from the mesh.
    for (int e = 0; e < orig_mesh.GetNE(); e++)
    {
      int attr = orig_mesh.GetAttribute(e);
      if (!mat_marker[attr - 1])
      {
        elem_delete[e] = true;
        new_ne--;
      }
    }

    // Make sure to remove any boundary elements which are no longer attached to elements in
    // the domain.
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be >= 0)
      {
        int e1, e2;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        if ((e1 < 0 || elem_delete[e1]) && (e2 < 0 || elem_delete[e2]))
        {
          // Mpi::Print("Deleting an unattached boundary element!\n");
          bdr_delete[be] = true;
          new_nbdr--;
        }
      }
    }
    if (new_ne < orig_mesh.GetNE())
    {
      Mpi::Print("Removed {:d} unmarked domain elements from the mesh\n",
                 orig_mesh.GetNE() - new_ne);
    }
    if (new_nbdr < orig_mesh.GetNBE())
    {
      Mpi::Print("Removed {:d} unattached boundary elements from the mesh\n",
                 orig_mesh.GetNBE() - new_nbdr);
    }
  }
  int new_ne_step1 = new_ne;
  int new_nbdr_step1 = new_nbdr;

  if (add_bdr)
  {
    // Add new boundary elements at material interfaces or on the exterior boundary of the
    // simulation domain, if there is not already a boundary element present.
    MFEM_VERIFY(!orig_mesh.Nonconforming(), "Adding material interface boundary elements "
                                            "is not supported for nonconforming meshes!");
    int add_bdr_ext = 0, add_bdr_int = 0;
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;
        orig_mesh.GetFaceElements(f, &e1, &e2);

        bool no_e1 = (e1 < 0 || elem_delete[e1]);
        bool no_e2 = (e2 < 0 || elem_delete[e2]);
        if ((no_e1 || no_e2) && !(no_e1 && no_e2))
        {
          // Mpi::Print("Adding exterior boundary element!\n");
          add_bdr_faces[f] = 1;
          add_bdr_ext++;
        }
        else if (orig_mesh.GetAttribute(e1) != orig_mesh.GetAttribute(e2))
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
    // Add new boundary elements at interfaces between elements belonging to different
    // subdomains. This uses similar code to mfem::Mesh::PrintWithPartitioning.
    MFEM_VERIFY(partitioning, "Cannot add subdomain interface boundary elements without "
                              "supplied mesh partitioning!");
    MFEM_VERIFY(!orig_mesh.Nonconforming(), "Adding subdomain interface boundary elements "
                                            "is not supported for nonconforming meshes!");
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      const int &be = orig_bdr_faces[f];
      if (be < 0 && add_bdr_faces[f] < 0)
      {
        int e1, e2;

        orig_mesh.GetFaceElements(f, &e1, &e2);
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
      Mpi::Print("Added {:d} boundary elements for subdomain interfaces to the mesh\n",
                 new_nbdr - new_nbdr_step2);
    }
  }

  // Create the new mesh.
  if (new_ne == new_ne_step1 && new_ne_step1 == new_ne_step2 &&
      new_ne_step2 == orig_mesh.GetNE() && new_nbdr == new_nbdr_step1 &&
      new_nbdr_step1 == new_nbdr_step2 && new_nbdr_step2 == orig_mesh.GetNBE())
  {
    return new_attr_map;
  }
  mfem::Mesh new_mesh(orig_mesh.Dimension(), orig_mesh.GetNV(), new_ne, new_nbdr,
                      orig_mesh.SpaceDimension());

  // Copy vertices and non-deleted domain and boundary elements.
  for (int v = 0; v < orig_mesh.GetNV(); v++)
  {
    new_mesh.AddVertex(orig_mesh.GetVertex(v));
  }
  for (int e = 0; e < orig_mesh.GetNE(); e++)
  {
    if (!elem_delete[e])
    {
      mfem::Element *el = orig_mesh.GetElement(e)->Duplicate(&new_mesh);
      new_mesh.AddElement(el);
    }
  }
  for (int be = 0; be < orig_mesh.GetNBE(); be++)
  {
    if (!bdr_delete[be])
    {
      mfem::Element *el = orig_mesh.GetBdrElement(be)->Duplicate(&new_mesh);
      new_mesh.AddBdrElement(el);
    }
  }

  // Add new boundary elements.
  if (add_bdr || add_subdomain)
  {
    auto FlipVertices = [](mfem::Element *el)
    {
      mfem::Array<int> v;
      el->GetVertices(v);
      std::reverse(v.begin(), v.end());
      el->SetVertices(v.HostRead());
    };

    // 1-based, some boundary attributes may be empty since they were removed from the
    // original mesh, but to keep indices the same as config file we don't compact the
    // list.
    int max_bdr_attr = orig_mesh.bdr_attributes.Size() ? orig_mesh.bdr_attributes.Max() : 0;
    for (int f = 0; f < orig_mesh.GetNumFaces(); f++)
    {
      if (add_bdr_faces[f] > 0)
      {
        // Assign new unique attribute based on attached elements (we want the material
        // properties on the face to average those on the elements). This is used later on
        // when integrating the transmission condition on the subdomain interface. Save the
        // inverse so that the attributes of e1 and e2 can be easily referenced using the
        // new attribute. Since attributes are in 1-based indexing, a, b > 0.
        int e1, e2, a = 0, b = 0;
        orig_mesh.GetFaceElements(f, &e1, &e2);
        bool no_e1 = (e1 < 0 || elem_delete[e1]);
        bool no_e2 = (e2 < 0 || elem_delete[e2]);
        if (!no_e1 && !no_e2)
        {
          a = std::max(orig_mesh.GetAttribute(e1), orig_mesh.GetAttribute(e2));
          b = (a == orig_mesh.GetAttribute(e1)) ? orig_mesh.GetAttribute(e2)
                                                : orig_mesh.GetAttribute(e1);
        }
        else if (!no_e1)
        {
          a = orig_mesh.GetAttribute(e1);
          b = 0;
        }
        else if (!no_e2)
        {
          a = orig_mesh.GetAttribute(e2);
          b = 0;
        }
        MFEM_VERIFY(a + b > 0, "Invalid new boundary element attribute!");
        int new_attr = max_bdr_attr +
                       (b > 0 ? (a * (a - 1)) / 2 + b : a);  // At least max_bdr_attr + 1
        if (new_attr_map.find(new_attr) == new_attr_map.end())
        {
          new_attr_map.emplace(new_attr, std::array<int, 2>{a, b});
        }

        // Add the boundary elements with the new boundary attribute.
        mfem::Element *el = orig_mesh.GetFace(f)->Duplicate(&new_mesh);
        el->SetAttribute(new_attr);
        new_mesh.AddBdrElement(el);
        if (add_bdr_faces[f] > 1)
        {
          // Flip order of vertices to reverse normal direction of second added element.
          el = orig_mesh.GetFace(f)->Duplicate(&new_mesh);
          FlipVertices(el);
          el->SetAttribute(new_attr);
          new_mesh.AddBdrElement(el);
          // Mpi::Print("Adding two BE with attr {:d} from elements {:d} and {:d}\n",
          //            new_attr, a, b);
        }
      }
    }
  }

  // Finalize new mesh and replace the old one. If a curved mesh, set up the new mesh by
  // projecting nodes onto the new mesh for the non-trimmed vdofs (accounts for new
  // boundary elements too since no new dofs are added). See the MFEM trimmer miniapp for
  // reference. After we have copied the high-order nodes information, topological changes
  // in Mesh::Finalize are OK (with refine = true).
  constexpr bool generate_bdr = false, refine = true, fix_orientation = true;
  new_mesh.FinalizeTopology(generate_bdr);
  new_mesh.RemoveUnusedVertices();
  if (orig_mesh.GetNodes())
  {
    const mfem::GridFunction *nodes = orig_mesh.GetNodes();
    const mfem::FiniteElementSpace *fespace = nodes->FESpace();

    mfem::Ordering::Type ordering = fespace->GetOrdering();
    int order = fespace->GetMaxElementOrder();
    int sdim = orig_mesh.SpaceDimension();
    bool discont =
        dynamic_cast<const mfem::L2_FECollection *>(fespace->FEColl()) != nullptr;

    new_mesh.SetCurvature(order, discont, sdim, ordering);
    mfem::GridFunction *new_nodes = new_mesh.GetNodes();
    const mfem::FiniteElementSpace *new_fespace = new_nodes->FESpace();

    // The element loop works because we know the mapping from old_mesh to new_mesh element
    // indices from the insertion order.
    mfem::Array<int> vdofs, new_vdofs;
    mfem::Vector loc_vec;
    int te = 0;
    for (int e = 0; e < orig_mesh.GetNE(); e++)
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
  new_mesh.Finalize(refine, fix_orientation);
  orig_mesh = std::move(new_mesh);
  return new_attr_map;
}

std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm comm,
                                              const std::unique_ptr<mfem::Mesh> &smesh,
                                              const std::unique_ptr<int[]> &partitioning,
                                              const std::string &output_dir)
{
  // Take a serial mesh and partitioning on the root process and construct the global
  // parallel mesh. For now, prefer the MPI-based version. When constructing the ParMesh, we
  // pass arguments to ensure no topological changes (this isn't required since the serial
  // mesh was marked for refinement).
  constexpr bool generate_edges = true, refine = true, fix_orientation = true;
  if constexpr (false)
  {
    // Write each processor's component to file.
    std::string tmp = output_dir;
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
    auto pmesh =
        std::make_unique<mfem::ParMesh>(comm, fi, generate_edges, refine, fix_orientation);
    Mpi::Barrier(comm);
    if (Mpi::Root(comm))
    {
      std::filesystem::remove_all(tmp);  // Remove the temporary directory
    }
    return pmesh;
  }
  else
  {
    // Send each processor's component as a byte string.
    std::unique_ptr<mfem::ParMesh> pmesh;
    if (Mpi::Root(comm))
    {
      mfem::MeshPartitioner partitioner(*smesh, Mpi::Size(comm), partitioning.get());
      std::vector<MPI_Request> send_requests(Mpi::Size(comm) - 1, MPI_REQUEST_NULL);
      std::vector<std::string> so;
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
        if (i > 0)
        {
          int slen = static_cast<int>(so[i].length());
          MFEM_VERIFY(so[i].length() == (std::size_t)slen,
                      "Overflow error distributing parallel mesh!");
          MPI_Isend(so[i].data(), slen, MPI_CHAR, i, i, comm, &send_requests[i - 1]);
        }
      }
      std::istringstream fi(so[0]);  // This is never compressed
      pmesh = std::make_unique<mfem::ParMesh>(comm, fi, generate_edges, refine,
                                              fix_orientation);
      MPI_Waitall(static_cast<int>(send_requests.size()), send_requests.data(),
                  MPI_STATUSES_IGNORE);
    }
    else
    {
      MPI_Status status;
      int rlen;
      std::string si;
      MPI_Probe(0, Mpi::Rank(comm), comm, &status);
      MPI_Get_count(&status, MPI_CHAR, &rlen);
      si.resize(rlen);
      MPI_Recv(si.data(), rlen, MPI_CHAR, 0, Mpi::Rank(comm), comm, MPI_STATUS_IGNORE);
      std::istringstream fi(si);
      // std::istringstream fi(zlib::DecompressString(si));
      pmesh = std::make_unique<mfem::ParMesh>(comm, fi, generate_edges, refine,
                                              fix_orientation);
    }
    return pmesh;
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

void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &pmesh, double length_scale,
                            const std::string &serial_mesh_file)
{
  // Write the parallel mesh to a stream as a serial mesh, then read back in and partition
  // using METIS.
  MPI_Comm comm = pmesh->GetComm();
  constexpr bool generate_edges = true, refine = true, fix_orientation = true,
                 generate_bdr = false;
  std::unique_ptr<mfem::Mesh> smesh;
  std::unique_ptr<int[]> partitioning;
  if constexpr (false)
  {
    // Write the serial mesh to a stream and read that through the Mesh constructor.
    std::stringstream fo;
    fo.precision(MSH_FLT_PRECISION);
    pmesh->PrintAsSerial(fo);
    pmesh.reset();
    if (Mpi::Root(comm))
    {
      smesh = std::make_unique<mfem::Mesh>(fo, generate_edges, refine, fix_orientation);
    }
  }
  else
  {
    // Directly ingest the generated Mesh and release the no longer needed memory.
    smesh = std::make_unique<mfem::Mesh>(pmesh->GetSerialMesh(0));
    pmesh.reset();
    if (!Mpi::Root(comm))
    {
      smesh.reset();
    }
  }
  if (Mpi::Root(comm))
  {
    smesh->FinalizeTopology(generate_bdr);
    smesh->Finalize(refine, fix_orientation);
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm));
  }

  // Construct the parallel mesh.
  pmesh = DistributeMesh(comm, smesh, partitioning);
  if (!serial_mesh_file.empty())
  {
    BlockTimer bt(Timer::IO);
    if (Mpi::Root(comm))
    {
      std::ofstream fo(serial_mesh_file);
      fo.precision(MSH_FLT_PRECISION);
      mesh::DimensionalizeMesh(*smesh, length_scale);
      smesh->Print(fo);  // Do not need to nondimensionalize the temporary mesh
    }
    Mpi::Barrier(comm);
  }
}

}  // namespace

}  // namespace palace
