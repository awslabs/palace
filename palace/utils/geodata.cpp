// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Eigen/Dense>
#include "fem/interpolator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"
#include "utils/meshio.hpp"
#include "utils/omp.hpp"
#include "utils/prettyprint.hpp"
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
std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &, bool,
                                     const config::BoundaryData &);

// Clean the provided serial mesh by removing unused domain and boundary elements.
void CleanMesh(std::unique_ptr<mfem::Mesh> &, const std::vector<int> &, bool = true);

// Create a new mesh by splitting all elements of the mesh into simplices or hexes
// (using tet-to-hex). Optionally preserves curvature of the original mesh by interpolating
// the high-order nodes with GSLIB.
void SplitMeshElements(std::unique_ptr<mfem::Mesh> &, bool, bool, bool = true);

// Optionally reorder mesh elements based on MFEM's internal reordeing tools for improved
// cache usage.
void ReorderMeshElements(mfem::Mesh &, bool = true);

// Check that mesh boundary conditions are given for external boundaries.
std::unordered_map<int, int> CheckMesh(const mfem::Mesh &, const config::BoundaryData &);

// Adding boundary elements for material interfaces and exterior boundaries, and "crack"
// desired internal boundary elements to disconnect the elements on either side.
int AddInterfaceBdrElements(const IoData &, std::unique_ptr<mfem::Mesh> &,
                            std::unordered_map<int, int> &);

// Generate element-based mesh partitioning, using either a provided file or METIS.
std::unique_ptr<int[]> GetMeshPartitioning(const mfem::Mesh &, int,
                                           const std::string & = "", bool = true);

// Given a serial mesh on the root processor and element partitioning, create a parallel
// mesh over the given communicator. The serial mesh is destroyed when no longer needed.
std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm, std::unique_ptr<mfem::Mesh> &,
                                              const int *, const std::string & = "");

// Rebalance a conformal mesh across processor ranks, using the MeshPartitioner. Gathers the
// mesh onto the root rank before scattering the partitioned mesh.
void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &);

}  // namespace

namespace mesh
{

std::unique_ptr<mfem::ParMesh> ReadMesh(const IoData &iodata, MPI_Comm comm)
{
  // If possible on root, read the serial mesh (converting format if necessary), and do all
  // necessary serial preprocessing. When finished, distribute the mesh to all processes.
  // Count disk I/O time separately for the mesh read from file.
  BlockTimer bt0(Timer::MESH_PREPROCESS);

  // If not doing any local adaptation, or performing conformal adaptation, we can use the
  // mesh partitioner.
  std::unique_ptr<mfem::Mesh> smesh;
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = (refinement.max_it > 0) || [&refinement]()
  {
    for (const auto &box : refinement.GetBoxes())
    {
      if (box.ref_levels > 0)
      {
        return true;
      }
    }
    for (const auto &sphere : refinement.GetSpheres())
    {
      if (sphere.ref_levels > 0)
      {
        return true;
      }
    }
    return false;
  }();
  const bool use_mesh_partitioner = !use_amr || !refinement.nonconformal;
  MPI_Comm node_comm;
  if (!use_mesh_partitioner)
  {
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
  }

  // Only one process per node reads the serial mesh.
  {
    BlockTimer bt1(Timer::IO);
    if ((use_mesh_partitioner && Mpi::Root(comm)) ||
        (!use_mesh_partitioner && Mpi::Root(node_comm)))
    {
      smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature, iodata.boundaries);
      MFEM_VERIFY(!(smesh->Nonconforming() && use_mesh_partitioner),
                  "Cannot use mesh partitioner on a nonconforming mesh!");
    }
    Mpi::Barrier(comm);
  }

  // Do some mesh preprocessing, and generate the partitioning.
  std::unique_ptr<int[]> partitioning;
  if (smesh)
  {
    // Check the the AMR specification and the mesh elements are compatible.
    const auto element_types = CheckElements(*smesh);
    MFEM_VERIFY(!use_amr || iodata.model.make_simplex || !element_types.has_hexahedra ||
                    refinement.nonconformal,
                "If there are tensor elements, AMR must be nonconformal!");
    MFEM_VERIFY(!use_amr || iodata.model.make_simplex || !element_types.has_prisms ||
                    refinement.nonconformal,
                "If there are wedge elements, AMR must be nonconformal!");
    MFEM_VERIFY(!use_amr || iodata.model.make_simplex || !element_types.has_pyramids ||
                    refinement.nonconformal,
                "If there are pyramid elements, AMR must be nonconformal!");

    // Clean up unused domain elements from the mesh.
    if (iodata.model.clean_unused_elements)
    {
      std::vector<int> attr_list;
      std::merge(iodata.domains.attributes.begin(), iodata.domains.attributes.end(),
                 iodata.domains.postpro.attributes.begin(),
                 iodata.domains.postpro.attributes.end(), std::back_inserter(attr_list));
      attr_list.erase(std::unique(attr_list.begin(), attr_list.end()), attr_list.end());
      CleanMesh(smesh, attr_list);
    }

    // Optionally convert mesh elements to simplices, for example in order to enable
    // conformal mesh refinement, or hexes.
    if (iodata.model.make_simplex || iodata.model.make_hex)
    {
      SplitMeshElements(smesh, iodata.model.make_simplex, iodata.model.make_hex);
    }

    // Optionally reorder elements (and vertices) based on spatial location after loading
    // the serial mesh.
    if (iodata.model.reorder_elements)
    {
      ReorderMeshElements(*smesh);
    }

    // Refine the serial mesh (not typically used, prefer parallel uniform refinement
    // instead).
    {
      int ne = smesh->GetNE();
      for (int l = 0; l < iodata.model.refinement.ser_uniform_ref_levels; l++)
      {
        smesh->UniformRefinement();
      }
      if (iodata.model.refinement.ser_uniform_ref_levels > 0)
      {
        Mpi::Print("Serial uniform mesh refinement levels added {:d} elements (initial = "
                   "{:d}, final = {:d})\n",
                   smesh->GetNE() - ne, ne, smesh->GetNE());
      }
    }

    // Check the final mesh, throwing warnings if there are exterior boundaries with no
    // associated boundary condition.
    auto face_to_be = CheckMesh(*smesh, iodata.boundaries);

    // Add new boundary elements for material interfaces if not present (with new unique
    // boundary attributes). Also duplicate internal boundary elements associated with
    // cracks if desired.
    if (iodata.model.crack_bdr_elements || iodata.model.add_bdr_elements)
    {
      // Split all internal (non periodic) boundary elements for boundary attributes where
      // BC are applied (not just postprocessing).
      while (AddInterfaceBdrElements(iodata, smesh, face_to_be) != 1)
      {
        // May require multiple calls due to early exit/retry approach.
      }
    }

    // Finally, finalize the serial mesh. Mark tetrahedral meshes for refinement. There
    // should be no need to fix orientation as this was done during initial mesh loading
    // from disk.
    constexpr bool refine = true, fix_orientation = false;
    smesh->Finalize(refine, fix_orientation);

    // Generate the mesh partitioning.
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), iodata.model.partitioning);
  }

  // Distribute the mesh.
  std::unique_ptr<mfem::ParMesh> pmesh;
  if (use_mesh_partitioner)
  {
    pmesh = DistributeMesh(comm, smesh, partitioning.get(), iodata.problem.output);
  }
  else
  {
    // Send the preprocessed serial mesh and partitioning as a byte string.
    constexpr bool generate_edges = false, refine = true, fix_orientation = false;
    std::string so;
    int slen = 0;
    if (smesh)
    {
      std::ostringstream fo(std::stringstream::out);
      // fo << std::fixed;
      fo << std::scientific;
      fo.precision(MSH_FLT_PRECISION);
      smesh->Print(fo);
      smesh.reset();  // Root process needs to rebuild the mesh to ensure consistency with
                      // the saved serial mesh (refinement marking, for example)
      so = fo.str();
      // so = zlib::CompressString(fo.str());
      slen = static_cast<int>(so.size());
      MFEM_VERIFY(so.size() == (std::size_t)slen, "Overflow in stringbuffer size!");
    }
    Mpi::Broadcast(1, &slen, 0, node_comm);
    if (so.empty())
    {
      so.resize(slen);
    }
    Mpi::Broadcast(slen, so.data(), 0, node_comm);
    {
      std::istringstream fi(so);
      // std::istringstream fi(zlib::DecompressString(so));
      smesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
      so.clear();
    }
    if (refinement.nonconformal && use_amr)
    {
      smesh->EnsureNCMesh(true);
    }
    if (!partitioning)
    {
      partitioning = std::make_unique<int[]>(smesh->GetNE());
    }
    Mpi::Broadcast(smesh->GetNE(), partitioning.get(), 0, node_comm);
    MPI_Comm_free(&node_comm);
    pmesh = std::make_unique<mfem::ParMesh>(comm, *smesh, partitioning.get());
    smesh.reset();
  }

  if constexpr (false)
  {
    auto tmp = fs::path(iodata.problem.output) / "tmp";
    if (Mpi::Root(comm) && !fs::exists(tmp))
    {
      fs::create_directories(tmp);
    }
    int width = 1 + static_cast<int>(std::log10(Mpi::Size(comm) - 1));
    std::unique_ptr<mfem::Mesh> gsmesh =
        LoadMesh(iodata.model.mesh, iodata.model.remove_curvature, iodata.boundaries);
    std::unique_ptr<int[]> gpartitioning = GetMeshPartitioning(*gsmesh, Mpi::Size(comm));
    mfem::ParMesh gpmesh(comm, *gsmesh, gpartitioning.get(), 0);
    {
      std::string pfile =
          mfem::MakeParFilename(tmp.string() + "part.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      // mfem::ofgzstream fo(pfile, true);  // Use zlib compression if available
      fo.precision(MSH_FLT_PRECISION);
      gpmesh.ParPrint(fo);
    }
    {
      std::string pfile =
          mfem::MakeParFilename(tmp.string() + "final.", Mpi::Rank(comm), ".mesh", width);
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
  if (iodata.solver.linear.mg_use_mesh && iodata.solver.linear.mg_max_levels > 1)
  {
    mesh.reserve(1 + uniform_ref_levels + max_region_ref_levels);
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the coarse mesh so that the refinements are still true refinements of
  // the original mesh (required for geometric multigrid). Otherwise, it happens after
  // refinement.
  if (iodata.model.reorient_tet_mesh && mesh.capacity() > 1)
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

  // Simplex meshes need to be re-finalized in order to use local refinement (see
  // the docstring for mfem::Mesh::UniformRefinement).
  const auto element_types = mesh::CheckElements(*mesh.back());
  if (element_types.has_simplices && uniform_ref_levels > 0 &&
      (max_region_ref_levels > 0 || iodata.model.refinement.max_it > 0))
  {
    constexpr bool refine = true, fix_orientation = false;
    Mpi::Print("\nFlattening mesh sequence:\n Local mesh refinement will start from the "
               "final uniformly-refined mesh\n");
    mesh.erase(mesh.begin(), mesh.end() - 1);
    mesh.back()->Finalize(refine, fix_orientation);
  }

  // Proceed with region-based refinement, level-by-level for all regions. Currently support
  // box and sphere region shapes. Any overlap between regions is ignored (take the union,
  // don't double-refine).
  MFEM_VERIFY(
      max_region_ref_levels == 0 ||
          !(element_types.has_hexahedra || element_types.has_prisms ||
            element_types.has_pyramids) ||
          mesh.back()->Nonconforming(),
      "Region-based refinement for non-simplex meshes requires a nonconformal mesh!");
  const bool use_nodes = (mesh.back()->GetNodes() != nullptr);
  const int ref = use_nodes ? mesh.back()->GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  const int dim = mesh.back()->SpaceDimension();
  int region_ref_level = 0;
  while (region_ref_level < max_region_ref_levels)
  {
    // Mark elements for refinement in all regions. An element is marked for refinement if
    // any of its vertices are inside any refinement region for the given level.
    mfem::Array<mfem::Refinement> refinements;
    for (int i = 0; i < mesh.back()->GetNE(); i++)
    {
      bool refine = false;
      mfem::DenseMatrix pointmat;
      if (use_nodes)
      {
        mfem::ElementTransformation &T = *mesh.back()->GetElementTransformation(i);
        mfem::RefinedGeometry *RefG =
            mfem::GlobGeometryRefiner.Refine(T.GetGeometryType(), ref);
        T.Transform(RefG->RefPts, pointmat);
      }
      else
      {
        const int *verts = mesh.back()->GetElement(i)->GetVertices();
        const int nv = mesh.back()->GetElement(i)->GetNVertices();
        pointmat.SetSize(dim, nv);
        for (int j = 0; j < nv; j++)
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
        refinements.Append(mfem::Refinement(i));
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
        refinements.Append(mfem::Refinement(i));
      }
    }

    // Do the refinement. For tensor element meshes, this may make the mesh nonconforming
    // (adds hanging nodes).
    if (mesh.capacity() > 1)
    {
      mesh.emplace_back(std::make_unique<mfem::ParMesh>(*mesh.back()));
    }
    mesh.back()->GeneralRefinement(refinements, -1);
    region_ref_level++;
  }
  if (max_region_ref_levels > 0 && mesh.capacity() == 1)
  {
    RebalanceMesh(iodata, mesh[0]);
  }

  // Prior to MFEM's PR #1046, the tetrahedral mesh required reorientation after all mesh
  // refinement in order to define higher-order Nedelec spaces on it. This is technically
  // not required after MFEM's PR #1046, but in case you want to be absolutely sure, we
  // reorient only the mesh after refinement if there is a single mesh (doesn't work with
  // h-refinement geometric multigrid).
  if (iodata.model.reorient_tet_mesh && mesh.capacity() == 1)
  {
    PalacePragmaDiagnosticPush
    PalacePragmaDiagnosticDisableDeprecated
    mesh[0]->ReorientTetMesh();
    PalacePragmaDiagnosticPop
  }

  // Print some mesh information.
  mfem::Vector bbmin, bbmax;
  GetAxisAlignedBoundingBox(*mesh[0], bbmin, bbmax);
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
  PalacePragmaOmp(parallel for schedule(static))
  for (int i = 0; i < mesh.GetNV(); i++)
  {
    double *v = mesh.GetVertex(i);
    std::transform(v, v + mesh.SpaceDimension(), v, [L](double val) { return val * L; });
  }
  if (auto *pmesh = dynamic_cast<mfem::ParMesh *>(&mesh))
  {
    PalacePragmaOmp(parallel for schedule(static))
    for (int i = 0; i < pmesh->face_nbr_vertices.Size(); i++)
    {
      double *v = pmesh->face_nbr_vertices[i]();
      std::transform(v, v + mesh.SpaceDimension(), v, [L](double val) { return val * L; });
    }
  }
  if (mesh.GetNodes())
  {
    *mesh.GetNodes() *= L;
    if (auto *pnodes = dynamic_cast<mfem::ParGridFunction *>(mesh.GetNodes()))
    {
      pnodes->FaceNbrData() *= L;
    }
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

std::vector<mfem::Geometry::Type> ElementTypeInfo::GetGeomTypes() const
{
  std::vector<mfem::Geometry::Type> geom_types;
  if (has_simplices)
  {
    geom_types.push_back(mfem::Geometry::TETRAHEDRON);
  }
  if (has_hexahedra)
  {
    geom_types.push_back(mfem::Geometry::CUBE);
  }
  if (has_prisms)
  {
    geom_types.push_back(mfem::Geometry::PRISM);
  }
  if (has_pyramids)
  {
    geom_types.push_back(mfem::Geometry::PYRAMID);
  }
  return geom_types;
}

ElementTypeInfo CheckElements(const mfem::Mesh &mesh)
{
  // MeshGenerator is reduced over the communicator. This checks for geometries on any
  // processor.
  auto meshgen = mesh.MeshGenerator();
  return {bool(meshgen & 1), bool(meshgen & 2), bool(meshgen & 4), bool(meshgen & 8)};
}

bool CheckRefinementFlags(const mfem::Mesh &mesh)
{
  bool marked = true;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const mfem::Element *el = mesh.GetElement(e);
    const int geom = el->GetGeometryType();
    if (geom == mfem::Geometry::TETRAHEDRON)
    {
      const mfem::Tetrahedron *tet = static_cast<const mfem::Tetrahedron *>(el);
      if (tet->GetRefinementFlag() == 0)
      {
        marked = false;
        break;
      }
    }
  }
  if (const auto *pmesh = dynamic_cast<const mfem::ParMesh *>(&mesh))
  {
    Mpi::GlobalAnd(1, &marked, pmesh->GetComm());
  }
  return marked;
}

void AttrToMarker(int max_attr, const int *attr_list, int attr_list_size,
                  mfem::Array<int> &marker, bool skip_invalid)
{
  MFEM_VERIFY(skip_invalid || attr_list_size == 0 ||
                  *std::max_element(attr_list, attr_list + attr_list_size) <= max_attr,
              "Invalid attribute number present ("
                  << *std::max_element(attr_list, attr_list + attr_list_size) << ")!");
  marker.SetSize(max_attr);
  if (attr_list_size == 1 && attr_list[0] == -1)
  {
    marker = 1;
  }
  else
  {
    marker = 0;
    for (int i = 0; i < attr_list_size; i++)
    {
      int attr = attr_list[i];
      if ((attr <= 0 || attr > max_attr) && skip_invalid)
      {
        continue;
      }
      MFEM_VERIFY(attr > 0, "Attribute number less than one!");
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeate attribute in attribute list!");
      marker[attr - 1] = 1;
    }
  }
}

void GetAxisAlignedBoundingBox(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
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
    auto BBUpdate =
        [&mesh, &dim](const int *v, int nv, mfem::Vector &min, mfem::Vector &max)
    {
      for (int j = 0; j < nv; j++)
      {
        const double *coord = mesh.GetVertex(v[j]);
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
    PalacePragmaOmp(parallel)
    {
      mfem::Vector loc_min(dim), loc_max(dim);
      for (int d = 0; d < dim; d++)
      {
        loc_min(d) = mfem::infinity();
        loc_max(d) = -mfem::infinity();
      }
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetBdrElement(i)->GetVertices();
          BBUpdate(verts, mesh.GetBdrElement(i)->GetNVertices(), loc_min, loc_max);
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetElement(i)->GetVertices();
          BBUpdate(verts, mesh.GetElement(i)->GetNVertices(), loc_min, loc_max);
        }
      }
      PalacePragmaOmp(critical(BBUpdate))
      {
        for (int d = 0; d < dim; d++)
        {
          min(d) = std::min(min(d), loc_min(d));
          max(d) = std::max(max(d), loc_max(d));
        }
      }
    }
  }
  else
  {
    mesh.GetNodes()->HostRead();
    const int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    auto BBUpdate = [&ref](mfem::GeometryRefiner &refiner, mfem::ElementTransformation &T,
                           mfem::DenseMatrix &pointmat, mfem::Vector &min,
                           mfem::Vector &max)
    {
      mfem::RefinedGeometry *RefG = refiner.Refine(T.GetGeometryType(), ref);
      T.Transform(RefG->RefPts, pointmat);
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
    PalacePragmaOmp(parallel)
    {
      mfem::Vector loc_min(dim), loc_max(dim);
      for (int d = 0; d < dim; d++)
      {
        loc_min(d) = mfem::infinity();
        loc_max(d) = -mfem::infinity();
      }
      mfem::GeometryRefiner refiner;
      mfem::IsoparametricTransformation T;
      mfem::DenseMatrix pointmat;
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetBdrElementTransformation(i, &T);
          BBUpdate(refiner, T, pointmat, loc_min, loc_max);
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetElementTransformation(i, &T);
          BBUpdate(refiner, T, pointmat, loc_min, loc_max);
        }
      }
      PalacePragmaOmp(critical(BBUpdate))
      {
        for (int d = 0; d < dim; d++)
        {
          min(d) = std::min(min(d), loc_min(d));
          max(d) = std::max(max(d), loc_max(d));
        }
      }
    }
  }
  Mpi::GlobalMin(dim, min.HostReadWrite(), mesh.GetComm());
  Mpi::GlobalMax(dim, max.HostReadWrite(), mesh.GetComm());
}

double BoundingBox::Area() const
{
  return 4.0 * CVector3dMap(axes[0].data()).cross(CVector3dMap(axes[1].data())).norm();
}

double BoundingBox::Volume() const
{
  return planar ? 0.0 : 2.0 * CVector3dMap(axes[2].data()).norm() * Area();
}

std::array<std::array<double, 3>, 3> BoundingBox::Normals() const
{
  std::array<std::array<double, 3>, 3> normals = {axes[0], axes[1], axes[2]};
  Vector3dMap(normals[0].data()).normalize();
  Vector3dMap(normals[1].data()).normalize();
  Vector3dMap(normals[2].data()).normalize();
  return normals;
}

std::array<double, 3> BoundingBox::Lengths() const
{
  return {2.0 * CVector3dMap(axes[0].data()).norm(),
          2.0 * CVector3dMap(axes[1].data()).norm(),
          2.0 * CVector3dMap(axes[2].data()).norm()};
}

std::array<double, 3> BoundingBox::Deviations(const std::array<double, 3> &direction) const
{
  const auto eig_dir = CVector3dMap(direction.data());
  std::array<double, 3> deviation_deg;
  for (std::size_t i = 0; i < 3; i++)
  {
    deviation_deg[i] =
        std::acos(std::min(1.0, std::abs(eig_dir.normalized().dot(
                                    CVector3dMap(axes[i].data()).normalized())))) *
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
}

// Helper for collecting a point cloud from a mesh, used in calculating bounding boxes and
// bounding balls. Returns the dominant rank, for which the vertices argument will be
// filled, while all other ranks will have an empty vector. Vertices are de-duplicated to a
// certain floating point precision.
int CollectPointCloudOnRoot(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr, std::vector<Eigen::Vector3d> &vertices)
{
  if (!mesh.GetNodes())
  {
    // Linear mesh, work with element vertices directly.
    PalacePragmaOmp(parallel)
    {
      std::unordered_set<int> vertex_indices;
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetBdrElement(i)->GetVertices();
          vertex_indices.insert(verts, verts + mesh.GetBdrElement(i)->GetNVertices());
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetElement(i)->GetVertices();
          vertex_indices.insert(verts, verts + mesh.GetElement(i)->GetNVertices());
        }
      }
      PalacePragmaOmp(critical(PointCloud))
      {
        for (auto i : vertex_indices)
        {
          const auto &vx = mesh.GetVertex(i);
          vertices.emplace_back(vx[0], vx[1], vx[2]);
        }
      }
    }
  }
  else
  {
    // Curved mesh, need to process point matrices.
    const int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    auto AddPoints = [&](mfem::GeometryRefiner &refiner, mfem::ElementTransformation &T,
                         mfem::DenseMatrix &pointmat,
                         std::vector<Eigen::Vector3d> &loc_vertices)
    {
      mfem::RefinedGeometry *RefG = refiner.Refine(T.GetGeometryType(), ref);
      T.Transform(RefG->RefPts, pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        loc_vertices.emplace_back(pointmat(0, j), pointmat(1, j), pointmat(2, j));
      }
    };
    PalacePragmaOmp(parallel)
    {
      mfem::GeometryRefiner refiner;
      mfem::IsoparametricTransformation T;
      mfem::DenseMatrix pointmat;  // 3 x N
      std::vector<Eigen::Vector3d> loc_vertices;
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetBdrElementTransformation(i, &T);
          AddPoints(refiner, T, pointmat, loc_vertices);
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetElementTransformation(i, &T);
          AddPoints(refiner, T, pointmat, loc_vertices);
        }
      }
      PalacePragmaOmp(critical(PointCloud))
      {
        for (const auto &v : loc_vertices)
        {
          vertices.push_back(v);
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

// Compute the distance from a point orthogonal to the list of normal axes, relative to
// the given origin.
auto PerpendicularDistance(const std::initializer_list<Eigen::Vector3d> &normals,
                           const Eigen::Vector3d &origin, const Eigen::Vector3d &v)
{
  Eigen::Vector3d v0 = v - origin;
  for (const auto &n : normals)
  {
    v0 -= n.dot(v0) * n;
  }
  return v0.norm();
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
    auto p_111 =
        std::max_element(vertices.begin(), vertices.end(),
                         [p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                         { return (x - *p_000).norm() < (y - *p_000).norm(); });
    p_000 = std::max_element(vertices.begin(), vertices.end(),
                             [p_111](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                             { return (x - *p_111).norm() < (y - *p_111).norm(); });
    MFEM_ASSERT(std::max_element(vertices.begin(), vertices.end(),
                                 [p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                 { return (x - *p_000).norm() < (y - *p_000).norm(); }) ==
                    p_111,
                "p_000 and p_111 must be mutually opposing points!");

    // Define a diagonal of the ASSUMED cuboid bounding box.
    const auto &v_000 = *p_000;
    const auto &v_111 = *p_111;
    MFEM_VERIFY(&v_000 != &v_111, "Minimum and maximum extents cannot be identical!");
    const Eigen::Vector3d origin = v_000;
    const Eigen::Vector3d n_1 = (v_111 - v_000).normalized();

    // Find the vertex furthest from the diagonal axis. We cannot know yet if this defines
    // (001) or (011).
    const auto &t_0 = *std::max_element(vertices.begin(), vertices.end(),
                                        [&](const auto &x, const auto &y)
                                        {
                                          return PerpendicularDistance({n_1}, origin, x) <
                                                 PerpendicularDistance({n_1}, origin, y);
                                        });
    MFEM_VERIFY(&t_0 != &v_000 && &t_0 != &v_111, "Vertices are degenerate!");

    // Use the discovered vertex to define a second direction and thus a plane. n_1 and n_2
    // now define a planar coordinate system intersecting the main diagonal, and two
    // opposite edges of the cuboid.
    const Eigen::Vector3d n_2 =
        ((t_0 - origin) - ((t_0 - origin).dot(n_1) * n_1)).normalized();

    // Collect the furthest point from the plane to determine if the box is planar. Look for
    // a component that maximizes distance from the planar system: complete the axes with a
    // cross, then use a dot product to pick the greatest deviation.
    constexpr double rel_tol = 1.0e-6;
    auto max_distance = PerpendicularDistance(
        {n_1, n_2}, origin,
        *std::max_element(vertices.begin(), vertices.end(),
                          [&](const auto &x, const auto &y)
                          {
                            return PerpendicularDistance({n_1, n_2}, origin, x) <
                                   PerpendicularDistance({n_1, n_2}, origin, y);
                          }));
    box.planar = (max_distance < rel_tol * (v_111 - v_000).norm());

    // For the non-planar case, collect points furthest from the plane and choose the one
    // closest to the origin as the next vertex which might be (001) or (011).
    const auto &t_1 = [&]()
    {
      if (box.planar)
      {
        return t_0;
      }
      std::vector<Eigen::Vector3d> vertices_out_of_plane;
      std::copy_if(vertices.begin(), vertices.end(),
                   std::back_inserter(vertices_out_of_plane),
                   [&](const auto &v)
                   {
                     return std::abs(PerpendicularDistance({n_1, n_2}, origin, v) -
                                     max_distance) < rel_tol * max_distance;
                   });
      return *std::min_element(vertices_out_of_plane.begin(), vertices_out_of_plane.end(),
                               [&](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                               { return (x - origin).norm() < (y - origin).norm(); });
    }();

    // Given candidates t_0 and t_1, the closer to origin defines v_001.
    const bool t_0_gt_t_1 = (t_0 - origin).norm() > (t_1 - origin).norm();
    const auto &v_001 = t_0_gt_t_1 ? t_1 : t_0;
    const auto &v_011 = box.planar ? v_111 : (t_0_gt_t_1 ? t_0 : t_1);

    // Compute the center as halfway along the main diagonal.
    Vector3dMap(box.center.data()) = 0.5 * (v_000 + v_111);

    // Compute the box axes. Using the 4 extremal points, we find the first two axes as the
    // edges which are closest to perpendicular. For a perfect rectangular prism point
    // cloud, we could instead compute the axes and length in each direction using the
    // found edges of the cuboid, but this does not work for non-rectangular prism
    // cross-sections or pyramid shapes.
    {
      const auto [e_0, e_1] = [&v_000, &v_001, &v_011, &v_111]()
      {
        std::array<const Eigen::Vector3d *, 4> verts = {&v_000, &v_001, &v_011, &v_111};
        Eigen::Vector3d e_0 = Eigen::Vector3d::Zero(), e_1 = Eigen::Vector3d::Zero();
        double dot_min = mfem::infinity();
        for (int i_0 = 0; i_0 < 4; i_0++)
        {
          for (int j_0 = i_0 + 1; j_0 < 4; j_0++)
          {
            for (int i_1 = 0; i_1 < 4; i_1++)
            {
              for (int j_1 = i_1 + 1; j_1 < 4; j_1++)
              {
                if (i_1 == i_0 && j_1 == j_0)
                {
                  continue;
                }
                const auto e_ij_0 = (*verts[j_0] - *verts[i_0]).normalized();
                const auto e_ij_1 = (*verts[j_1] - *verts[i_1]).normalized();
                const auto dot = std::abs(e_ij_0.dot(e_ij_1));
                if (dot < dot_min)
                {
                  dot_min = dot;
                  e_0 = e_ij_0;
                  e_1 = e_ij_1;
                  if (dot_min < rel_tol)
                  {
                    return std::make_pair(e_0, e_1);
                  }
                }
              }
            }
          }
        }
        return std::make_pair(e_0, e_1);
      }();
      Vector3dMap(box.axes[0].data()) = e_0;
      Vector3dMap(box.axes[1].data()) = e_1;
      Vector3dMap(box.axes[2].data()) =
          box.planar ? Eigen::Vector3d::Zero() : e_0.cross(e_1);
    }

    // Scale axes by length of the box in each direction.
    std::array<double, 3> l = {0.0};
    for (const auto &v : {v_000, v_001, v_011, v_111})
    {
      const auto v_0 = v - Vector3dMap(box.center.data());
      l[0] = std::max(l[0], std::abs(v_0.dot(Vector3dMap(box.axes[0].data()))));
      l[1] = std::max(l[1], std::abs(v_0.dot(Vector3dMap(box.axes[1].data()))));
      l[2] = std::max(l[2], std::abs(v_0.dot(Vector3dMap(box.axes[2].data()))));
    }
    Vector3dMap(box.axes[0].data()) *= l[0];
    Vector3dMap(box.axes[1].data()) *= l[1];
    Vector3dMap(box.axes[2].data()) *= l[2];

    // Make sure the longest dimension comes first.
    std::sort(box.axes.begin(), box.axes.end(), [](const auto &x, const auto &y)
              { return CVector3dMap(x.data()).norm() > CVector3dMap(y.data()).norm(); });
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, box.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, box.axes.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &box.planar, dominant_rank, comm);

  return box;
}

// For the public interface, we can use a BoundingBox as a generalization of a BoundingBall.
// Internally, however, it's nice to work with a specific ball data type.
struct BoundingBall
{
  Eigen::Vector3d origin;
  double radius;
  bool planar;
};

// Use 4 points to define a sphere in 3D. If the points are coplanar, 3 of them are used to
// define a circle which is interpreted as the equator of the sphere. We assume the points
// are unique and not collinear.
BoundingBall SphereFromPoints(const std::vector<std::size_t> &indices,
                              const std::vector<Eigen::Vector3d> &vertices)
{
  // Given 0 or 1 points, just return a radius of 0.
  MFEM_VERIFY(
      indices.size() <= 4,
      "Determining a sphere in 3D requires 4 points (and a circle requires 3 points)!");
  BoundingBall ball;
  ball.planar = (indices.size() < 4);
  if (indices.size() < 2)
  {
    ball.origin = Eigen::Vector3d::Zero();
    ball.radius = 0.0;
    return ball;
  }

  // For two points, construct a circle with the segment as its diameter. This could also
  // handle the collinear case for more than 2 points.
  if (indices.size() == 2)
  {
    ball.origin = 0.5 * (vertices[indices[0]] + vertices[indices[1]]);
    ball.radius = (vertices[indices[0]] - ball.origin).norm();
    return ball;
  }

  // Check for coplanarity.
  constexpr double rel_tol = 1.0e-6;
  const Eigen::Vector3d AB = vertices[indices[1]] - vertices[indices[0]];
  const Eigen::Vector3d AC = vertices[indices[2]] - vertices[indices[0]];
  const Eigen::Vector3d ABAC = AB.cross(AC);
  Eigen::Vector3d AD = Eigen::Vector3d::Zero();
  if (!ball.planar)
  {
    AD = vertices[indices[3]] - vertices[indices[0]];
    ball.planar = (std::abs(AD.dot(ABAC)) < rel_tol * AD.norm() * ABAC.norm());
  }

  // Construct a circle passing through 3 points.
  // See: https://en.wikipedia.org/wiki/Circumcircle#Higher_dimensions.
  if (ball.planar)
  {
    ball.origin = (0.5 / ABAC.squaredNorm()) *
                  ((AB.squaredNorm() * AC) - (AC.squaredNorm() * AB)).cross(ABAC);
    ball.radius = ball.origin.norm();
    ball.origin += vertices[indices[0]];
#if defined(MFEM_DEBUG)
    const auto r1 = (vertices[indices[1]] - ball.origin).norm();
    const auto r2 = (vertices[indices[2]] - ball.origin).norm();
    MFEM_VERIFY((1.0 - rel_tol) * ball.radius < r1 && r1 < (1.0 + rel_tol) * ball.radius &&
                    (1.0 - rel_tol) * ball.radius < r2 &&
                    r2 < (1.0 + rel_tol) * ball.radius,
                "Invalid circle calculated from 3 points!");
#endif
    return ball;
  }

  // Construct a sphere passing through 4 points.
  // See: https://steve.hollasch.net/cgindex/geometry/sphere4pts.html.
  Eigen::Matrix3d C;
  Eigen::Vector3d d;
  const auto s = vertices[indices[0]].squaredNorm();
  C.row(0) = AB.transpose();
  C.row(1) = AC.transpose();
  C.row(2) = AD.transpose();
  d(0) = 0.5 * (vertices[indices[1]].squaredNorm() - s);
  d(1) = 0.5 * (vertices[indices[2]].squaredNorm() - s);
  d(2) = 0.5 * (vertices[indices[3]].squaredNorm() - s);
  ball.origin = C.inverse() * d;  // 3x3 matrix inverse might be faster than general LU
                                  // if Eigen uses the explicit closed-form solution
  ball.radius = (vertices[indices[0]] - ball.origin).norm();
#if defined(MFEM_DEBUG)
  const auto r1 = (vertices[indices[1]] - ball.origin).norm();
  const auto r2 = (vertices[indices[2]] - ball.origin).norm();
  const auto r3 = (vertices[indices[3]] - ball.origin).norm();
  MFEM_VERIFY((1.0 - rel_tol) * ball.radius < r1 && r1 < (1.0 + rel_tol) * ball.radius &&
                  (1.0 - rel_tol) * ball.radius < r2 &&
                  r2 < (1.0 + rel_tol) * ball.radius &&
                  (1.0 - rel_tol) * ball.radius < r3 && r3 < (1.0 + rel_tol) * ball.radius,
              "Invalid sphere calculated from 3 points!");
#endif
  return ball;
}

BoundingBall Welzl(std::vector<std::size_t> P, std::vector<std::size_t> R,
                   const std::vector<Eigen::Vector3d> &vertices)
{
  // Base case.
  if (R.size() == 4 || P.empty())
  {
    return SphereFromPoints(R, vertices);
  }

  // Choose a p ∈ P randomly, and recurse for (P \ {p}, R). The set P has already been
  // randomized on input.
  const std::size_t p = P.back();
  P.pop_back();
  BoundingBall D = Welzl(P, R, vertices);

  // If p is outside the sphere, recurse for (P \ {p}, R U {p}).
  constexpr double rel_tol = 1.0e-6;
  if ((vertices[p] - D.origin).norm() >= (1.0 + rel_tol) * D.radius)
  {
    R.push_back(p);
    D = Welzl(P, R, vertices);
  }

  return D;
}

// Calculates a bounding ball from a point cloud using Welzl's algorithm, result is
// broadcast across all processes. We don't operate on the convex hull, since the number of
// points should be small enough that operating on the full set should be OK. If only three
// points are provided, the bounding circle is computed (likewise for if the points are
// coplanar).
BoundingBox BoundingBallFromPointCloud(MPI_Comm comm,
                                       const std::vector<Eigen::Vector3d> &vertices,
                                       int dominant_rank)
{
  BoundingBox ball;
  if (dominant_rank == Mpi::Rank(comm))
  {
    MFEM_VERIFY(vertices.size() >= 3,
                "A bounding ball requires a minimum of three vertices for this algorithm!");
    std::vector<std::size_t> indices(vertices.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Acceleration from https://informatica.vu.lt/journal/INFORMATICA/article/1251. Allow
    // for duplicate points and just add the 4 points to the end of the indicies list to be
    // considered first. The two points are not necessarily the maximizer of the distance
    // between all pairs, but they should be a good estimate.
    {
      auto p_1 = std::min_element(vertices.begin(), vertices.end(), EigenLE);
      auto p_2 = std::max_element(vertices.begin(), vertices.end(),
                                  [p_1](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                  { return (x - *p_1).norm() < (y - *p_1).norm(); });
      p_1 = std::max_element(vertices.begin(), vertices.end(),
                             [p_2](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                             { return (x - *p_2).norm() < (y - *p_2).norm(); });

      // Find the next point as the vertex furthest from the initial axis.
      const Eigen::Vector3d n_1 = (*p_2 - *p_1).normalized();
      auto p_3 = std::max_element(vertices.begin(), vertices.end(),
                                  [&](const auto &x, const auto &y)
                                  {
                                    return PerpendicularDistance({n_1}, *p_1, x) <
                                           PerpendicularDistance({n_1}, *p_1, y);
                                  });
      auto p_4 = std::max_element(vertices.begin(), vertices.end(),
                                  [p_3](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                  { return (x - *p_3).norm() < (y - *p_3).norm(); });
      MFEM_VERIFY(p_3 != p_1 && p_3 != p_2 && p_4 != p_1 && p_4 != p_2,
                  "Vertices are degenerate!");

      // Start search with these points, which should be roughly extremal. With the search
      // for p_3 done in an orthogonal direction, p_1, p_2, p_3, and p_4 should all be
      // unique.
      std::swap(indices[indices.size() - 1], indices[p_1 - vertices.begin()]);
      std::swap(indices[indices.size() - 2], indices[p_2 - vertices.begin()]);
      std::swap(indices[indices.size() - 3], indices[p_3 - vertices.begin()]);
      std::swap(indices[indices.size() - 4], indices[p_4 - vertices.begin()]);
    }

    // Randomly permute the point set.
    {
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(indices.begin(), indices.end() - 4, g);
    }

    // Compute the bounding ball.
    BoundingBall min_ball = Welzl(indices, {}, vertices);
    Vector3dMap(ball.center.data()) = min_ball.origin;
    Vector3dMap(ball.axes[0].data()) = Eigen::Vector3d(min_ball.radius, 0.0, 0.0);
    Vector3dMap(ball.axes[1].data()) = Eigen::Vector3d(0.0, min_ball.radius, 0.0);
    Vector3dMap(ball.axes[2].data()) = Eigen::Vector3d(0.0, 0.0, min_ball.radius);
    ball.planar = min_ball.planar;
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, ball.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, ball.axes.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &ball.planar, dominant_rank, comm);

  return ball;
}

}  // namespace

BoundingBox GetBoundingBox(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                           bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  return BoundingBoxFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
}

BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  return BoundingBallFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
}

double GetProjectedLength(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                          bool bdr, const std::array<double, 3> &dir)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  double length;
  if (dominant_rank == Mpi::Rank(mesh.GetComm()))
  {
    CVector3dMap direction(dir.data());
    auto Dot = [&](const auto &x, const auto &y)
    { return direction.dot(x) < direction.dot(y); };
    auto p_min = std::min_element(vertices.begin(), vertices.end(), Dot);
    auto p_max = std::max_element(vertices.begin(), vertices.end(), Dot);
    length = (*p_max - *p_min).dot(direction.normalized());
  }
  Mpi::Broadcast(1, &length, dominant_rank, mesh.GetComm());
  return length;
}

double GetDistanceFromPoint(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr, const std::array<double, 3> &origin, bool max)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  double dist;
  if (dominant_rank == Mpi::Rank(mesh.GetComm()))
  {
    CVector3dMap x0(origin.data());
    auto p =
        max ? std::max_element(vertices.begin(), vertices.end(),
                               [&x0](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                               { return (x - x0).norm() < (y - x0).norm(); })
            : std::min_element(vertices.begin(), vertices.end(),
                               [&x0](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                               { return (x - x0).norm() < (y - x0).norm(); });
    dist = (*p - x0).norm();
  }
  Mpi::Broadcast(1, &dist, dominant_rank, mesh.GetComm());
  return dist;
}

mfem::Vector GetSurfaceNormal(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                              bool average)
{
  int dim = mesh.SpaceDimension();
  mfem::IsoparametricTransformation T;
  mfem::Vector loc_normal(dim), normal(dim);
  normal = 0.0;
  bool init = false;
  auto UpdateNormal = [&](mfem::ElementTransformation &T)
  {
    const mfem::IntegrationPoint &ip = mfem::Geometries.GetCenter(T.GetGeometryType());
    T.SetIntPoint(&ip);
    mfem::CalcOrtho(T.Jacobian(), loc_normal);
    if (!init)
    {
      normal = loc_normal;
      init = true;
    }
    else
    {
      // Check orientation and make sure consistent on this process. If a boundary has
      // conflicting normal definitions, use the first value.
      if (loc_normal * normal < 0.0)
      {
        normal -= loc_normal;
      }
      else
      {
        normal += loc_normal;
      }
    }
  };
  if (mesh.Dimension() == mesh.SpaceDimension())
  {
    // Loop over boundary elements.
    for (int i = 0; i < mesh.GetNBE(); i++)
    {
      if (!marker[mesh.GetBdrAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetBdrElementTransformation(i, &T);
      UpdateNormal(T);
      if (!average)
      {
        break;
      }
    }
  }
  else
  {
    // Loop over domain elements.
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      if (!marker[mesh.GetAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetElementTransformation(i, &T);
      UpdateNormal(T);
      if (!average)
      {
        break;
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
    // No boundary elements are marked.
    normal = 0.0;
    return normal;
  }
  if (rank == Mpi::Rank(comm))
  {
    glob_normal = normal;
  }
  Mpi::Broadcast(dim, glob_normal.HostReadWrite(), rank, comm);
  if (average)
  {
    if (init && normal * glob_normal < 0.0)
    {
      normal.Neg();
    }
    Mpi::GlobalSum(dim, normal.HostReadWrite(), comm);
  }
  else
  {
    normal = glob_normal;
  }
  normal /= normal.Norml2();

  if constexpr (false)
  {
    Mpi::Print(comm, " Surface normal = ({:+.3e})", fmt::join(normal, ", "));
  }

  return normal;
}

double GetSurfaceArea(const mfem::ParMesh &mesh, const mfem::Array<int> &marker)
{
  double area = 0.0;
  PalacePragmaOmp(parallel reduction(+ : area))
  {
    mfem::IsoparametricTransformation T;
    PalacePragmaOmp(for schedule(static))
    for (int i = 0; i < mesh.GetNBE(); i++)
    {
      if (!marker[mesh.GetBdrAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetBdrElementTransformation(i, &T);
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), T.OrderJ());
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(j);
        T.SetIntPoint(&ip);
        area += ip.weight * T.Weight();
      }
    }
  }
  Mpi::GlobalSum(1, &area, mesh.GetComm());
  return area;
}

double GetVolume(const mfem::ParMesh &mesh, const mfem::Array<int> &marker)
{
  double volume = 0.0;
  PalacePragmaOmp(parallel reduction(+ : volume))
  {
    mfem::IsoparametricTransformation T;
    PalacePragmaOmp(for schedule(static))
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      if (!marker[mesh.GetAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetElementTransformation(i, &T);
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), T.OrderJ());
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(j);
        T.SetIntPoint(&ip);
        volume += ip.weight * T.Weight();
      }
    }
  }
  Mpi::GlobalSum(1, &volume, mesh.GetComm());
  return volume;
}

double RebalanceMesh(const IoData &iodata, std::unique_ptr<mfem::ParMesh> &mesh)
{
  BlockTimer bt0(Timer::REBALANCE);
  MPI_Comm comm = mesh->GetComm();
  if (iodata.model.refinement.save_adapt_mesh)
  {
    // Create a separate serial mesh to write to disk.
    auto sfile = fs::path(iodata.problem.output) / fs::path(iodata.model.mesh).stem();
    sfile += ".mesh";

    auto PrintSerial = [&](mfem::Mesh &smesh)
    {
      BlockTimer bt1(Timer::IO);
      if (Mpi::Root(comm))
      {
        std::ofstream fo(sfile);
        // mfem::ofgzstream fo(sfile, true);  // Use zlib compression if available
        // fo << std::fixed;
        fo << std::scientific;
        fo.precision(MSH_FLT_PRECISION);
        mesh::DimensionalizeMesh(smesh, iodata.GetMeshLengthScale());
        smesh.Mesh::Print(fo);  // Do not need to nondimensionalize the temporary mesh
      }
      Mpi::Barrier(comm);
    };

    if (mesh->Nonconforming())
    {
      mfem::ParMesh smesh(*mesh);
      mfem::Array<int> serial_partition(mesh->GetNE());
      serial_partition = 0;
      smesh.Rebalance(serial_partition);
      PrintSerial(smesh);
    }
    else
    {
      mfem::Mesh smesh = mesh->GetSerialMesh(0);
      PrintSerial(smesh);
    }
  }

  // If there is more than one processor, may perform rebalancing.
  if (Mpi::Size(comm) == 1)
  {
    return 1.0;
  }
  int min_elem, max_elem;
  min_elem = max_elem = mesh->GetNE();
  Mpi::GlobalMin(1, &min_elem, comm);
  Mpi::GlobalMax(1, &max_elem, comm);
  const double ratio = double(max_elem) / min_elem;
  const double tol = iodata.model.refinement.maximum_imbalance;
  if constexpr (false)
  {
    Mpi::Print("Rebalancing: max/min elements per processor = {:d}/{:d} (ratio = {:.3e}, "
               "tol = {:.3e})\n",
               max_elem, min_elem, ratio, tol);
  }
  if (ratio > tol)
  {
    if (mesh->Nonconforming())
    {
      mesh->Rebalance();
    }
    else
    {
      // Without access to a refinement tree, partitioning must be done on the root
      // processor and then redistributed.
      RebalanceConformalMesh(mesh);
    }
  }
  return ratio;
}

}  // namespace mesh

namespace
{

std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &mesh_file, bool remove_curvature,
                                     const config::BoundaryData &boundaries)
{
  // Read the (serial) mesh from the given mesh file. Handle preparation for refinement and
  // orientations here to avoid possible reorientations and reordering later on. MFEM
  // supports a native mesh format (.mesh), VTK/VTU, Gmsh, as well as some others. We use
  // built-in converters for the types we know, otherwise rely on MFEM to do the conversion
  // or error out if not supported.
  constexpr bool generate_edges = false, refine = false, fix_orientation = true;
  std::unique_ptr<mfem::Mesh> mesh;
  fs::path mesh_path(mesh_file);
  if (mesh_path.extension() == ".mphtxt" || mesh_path.extension() == ".mphbin" ||
      mesh_path.extension() == ".nas" || mesh_path.extension() == ".bdf")
  {
    // Put translated mesh in temporary string buffer.
    std::stringstream fi(std::stringstream::in | std::stringstream::out);
    // fi << std::fixed;
    fi << std::scientific;
    fi.precision(MSH_FLT_PRECISION);
    if (mesh_path.extension() == ".mphtxt" || mesh_path.extension() == ".mphbin")
    {
      mesh::ConvertMeshComsol(mesh_file, fi, remove_curvature);
      // mesh::ConvertMeshComsol(mesh_file, fo, remove_curvature);
    }
    else
    {
      mesh::ConvertMeshNastran(mesh_file, fi, remove_curvature);
      // mesh::ConvertMeshNastran(mesh_file, fo, remove_curvature);
    }
    mesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
  }
  else
  {
    // Otherwise, just rely on MFEM load the mesh.
    std::ifstream fi(mesh_file);
    if (!fi.good())
    {
      MFEM_ABORT("Unable to open mesh file \"" << mesh_file << "\"!");
    }
    mesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
  }
  if (remove_curvature)
  {
    mesh->SetCurvature(-1);
  }
  else
  {
    mesh->EnsureNodes();
  }
  if (!boundaries.periodic.empty())
  {
    auto periodic_mesh = std::move(mesh);
    for (const auto &data : boundaries.periodic)
    {
      mfem::Vector translation(data.translation.size());
      std::copy(data.translation.begin(), data.translation.end(), translation.GetData());
      auto periodic_mapping =
          periodic_mesh->CreatePeriodicVertexMapping({translation}, 1E-6);
      auto p_mesh = std::make_unique<mfem::Mesh>(
          mfem::Mesh::MakePeriodic(*periodic_mesh, periodic_mapping));
      periodic_mesh = std::move(p_mesh);
    }
    mesh = std::move(periodic_mesh);
  }
  return mesh;
}

template <typename T = std::vector<int>>
void TransferHighOrderNodes(const mfem::Mesh &orig_mesh, mfem::Mesh &new_mesh,
                            const T *elem_delete_map = nullptr)
{
  // This accounts for new boundary elements too since no new dofs are added. See the MFEM
  // trimmer miniapp for reference.
  MFEM_VERIFY(orig_mesh.GetNodes(), "No high-order nodes information to transfer!");
  const mfem::GridFunction *nodes = orig_mesh.GetNodes();
  const mfem::FiniteElementSpace *fespace = nodes->FESpace();
  mfem::Ordering::Type ordering = fespace->GetOrdering();
  int order = fespace->GetMaxElementOrder();
  int sdim = orig_mesh.SpaceDimension();
  bool discont =
      (dynamic_cast<const mfem::L2_FECollection *>(fespace->FEColl()) != nullptr);
  new_mesh.SetCurvature(order, discont, sdim, ordering);
  mfem::GridFunction *new_nodes = new_mesh.GetNodes();
  const mfem::FiniteElementSpace *new_fespace = new_nodes->FESpace();

  // Transfer dofs from the old mesh to the new ones. Either consider all elements (works
  // for orientation or numbering changes), or use the prescribed old to new element index
  // map.
  mfem::Array<int> vdofs;
  mfem::Vector loc_vec;
  for (int e = 0; e < orig_mesh.GetNE(); e++)
  {
    if (!elem_delete_map || (*elem_delete_map)[e] >= 0)
    {
      // No need for DofTransformation here since spaces are H1 or L2.
      fespace->GetElementVDofs(e, vdofs);
      nodes->GetSubVector(vdofs, loc_vec);
      new_fespace->GetElementVDofs(!elem_delete_map ? e : (*elem_delete_map)[e], vdofs);
      new_nodes->SetSubVector(vdofs, loc_vec);
    }
  }
}

void CleanMesh(std::unique_ptr<mfem::Mesh> &orig_mesh,
               const std::vector<int> &mat_attr_list, bool preserve_curvature)
{
  auto mat_marker = mesh::AttrToMarker(
      orig_mesh->attributes.Size() ? orig_mesh->attributes.Max() : 0, mat_attr_list, true);
  std::vector<int> elem_delete_map(orig_mesh->GetNE(), -1),
      bdr_elem_delete_map(orig_mesh->GetNBE(), -1);

  // Delete domain and boundary elements which have no associated material or BC attribute
  // from the mesh.
  int new_ne = 0;
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    if (mat_marker[orig_mesh->GetAttribute(e) - 1])
    {
      elem_delete_map[e] = new_ne++;
    }
  }

  // Make sure to remove any boundary elements which are no longer attached to elements in
  // the domain.
  int new_nbe = 0;
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    int f, o, e1, e2;
    orig_mesh->GetBdrElementFace(be, &f, &o);
    orig_mesh->GetFaceElements(f, &e1, &e2);
    bool no_e1 = (e1 < 0 || elem_delete_map[e1] < 0);
    bool no_e2 = (e2 < 0 || elem_delete_map[e2] < 0);
    if (!no_e1 || !no_e2)
    {
      bdr_elem_delete_map[be] = new_nbe++;
    }
    else if constexpr (false)
    {
      Mpi::Print("Deleting an unattached boundary element!\n");
    }
  }
  if (new_ne < orig_mesh->GetNE())
  {
    Mpi::Print("Removed {:d} unmarked domain elements from the mesh\n",
               orig_mesh->GetNE() - new_ne);
  }
  if (new_nbe < orig_mesh->GetNBE())
  {
    Mpi::Print("Removed {:d} unattached boundary elements from the mesh\n",
               orig_mesh->GetNBE() - new_nbe);
  }

  // Create the new mesh.
  if (new_ne == orig_mesh->GetNE() && new_nbe == orig_mesh->GetNBE())
  {
    return;
  }
  MFEM_VERIFY(!orig_mesh->Nonconforming(),
              "Mesh element cleaning is not supported for nonconforming meshes!");
  auto new_mesh =
      std::make_unique<mfem::Mesh>(orig_mesh->Dimension(), orig_mesh->GetNV(), new_ne,
                                   new_nbe, orig_mesh->SpaceDimension());

  // Copy vertices and non-deleted domain and boundary elements.
  for (int v = 0; v < orig_mesh->GetNV(); v++)
  {
    new_mesh->AddVertex(orig_mesh->GetVertex(v));
  }
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    if (elem_delete_map[e] >= 0)
    {
      mfem::Element *el = orig_mesh->GetElement(e)->Duplicate(new_mesh.get());
      new_mesh->AddElement(el);
    }
  }
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    if (bdr_elem_delete_map[be] >= 0)
    {
      mfem::Element *bdr_el = orig_mesh->GetBdrElement(be)->Duplicate(new_mesh.get());
      new_mesh->AddBdrElement(bdr_el);
    }
  }

  // Finalize the new mesh topology and replace the old mesh. If a curved mesh, set up the
  // new mesh by projecting nodes onto the new mesh for the non-trimmed vdofs. No need to
  // mark for refinement or fix orientations, since everything is copied from the previous
  // mesh.
  constexpr bool generate_bdr = false;
  new_mesh->FinalizeTopology(generate_bdr);
  new_mesh->RemoveUnusedVertices();  // Remove vertices from the deleted elements
  if (orig_mesh->GetNodes())
  {
    TransferHighOrderNodes(*orig_mesh, *new_mesh, &elem_delete_map);
  }
  orig_mesh = std::move(new_mesh);
}

mfem::Mesh MeshTetToHex(const mfem::Mesh &orig_mesh)
{
  // Courtesy of https://gist.github.com/pazner/e9376f77055c0918d7c43e034e9e5888, only
  // supports tetrahedral elements for now. Eventually should be expanded to support prism
  // and pyramid elements but this mixed mesh support requires a bit more work.
  MFEM_VERIFY(orig_mesh.Dimension() == 3, "Tet-to-hex conversion only supports 3D meshes!");
  {
    // This checks the local mesh on each process, but the assertion failing on any single
    // process will terminate the program.
    mfem::Array<mfem::Geometry::Type> geoms;
    orig_mesh.GetGeometries(3, geoms);
    MFEM_VERIFY(geoms.Size() == 1 && geoms[0] == mfem::Geometry::TETRAHEDRON,
                "Tet-to-hex conversion only works for pure tetrahedral meshes!");
  }

  // Add new vertices in every edge, face, and volume. Each tet is subdivided into 4 hexes,
  // and each triangular face subdivided into 3 quads.
  const int nv_tet = orig_mesh.GetNV();
  const int nedge_tet = orig_mesh.GetNEdges();
  const int nface_tet = orig_mesh.GetNFaces();
  const int ne_tet = orig_mesh.GetNE();
  const int nbe_tet = orig_mesh.GetNBE();
  const int nv = nv_tet + nedge_tet + nface_tet + ne_tet;
  const int ne = 4 * ne_tet;
  const int nbe = 3 * nbe_tet;
  mfem::Mesh hex_mesh(orig_mesh.Dimension(), nv, ne, nbe, orig_mesh.SpaceDimension());

  // Add original vertices.
  for (int v = 0; v < nv_tet; v++)
  {
    hex_mesh.AddVertex(orig_mesh.GetVertex(v));
  }

  // Add midpoints of edges, faces, and elements.
  auto AddCentroid = [&orig_mesh, &hex_mesh](const int *verts, int nv)
  {
    double coord[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < nv; i++)
    {
      for (int d = 0; d < orig_mesh.SpaceDimension(); d++)
      {
        coord[d] += orig_mesh.GetVertex(verts[i])[d] / nv;
      }
    }
    hex_mesh.AddVertex(coord);
  };
  {
    mfem::Array<int> verts;
    for (int e = 0; e < nedge_tet; ++e)
    {
      orig_mesh.GetEdgeVertices(e, verts);
      AddCentroid(verts.GetData(), verts.Size());
    }
  }
  for (int f = 0; f < nface_tet; ++f)
  {
    AddCentroid(orig_mesh.GetFace(f)->GetVertices(), orig_mesh.GetFace(f)->GetNVertices());
  }
  for (int e = 0; e < ne_tet; ++e)
  {
    AddCentroid(orig_mesh.GetElement(e)->GetVertices(),
                orig_mesh.GetElement(e)->GetNVertices());
  }

  // Connectivity of tetrahedron vertices to the edges.
  constexpr int tet_vertex_edge_map[4 * 3] = {0, 1, 2, 3, 0, 4, 1, 3, 5, 5, 4, 2};
  constexpr int tet_vertex_face_map[4 * 3] = {3, 2, 1, 3, 0, 2, 3, 1, 0, 0, 1, 2};
  constexpr int tri_vertex_edge_map[3 * 2] = {0, 2, 1, 0, 2, 1};

  // Add four hexahedra for each tetrahedron.
  {
    mfem::Array<int> edges, faces, orients;
    for (int e = 0; e < ne_tet; ++e)
    {
      const int *verts = orig_mesh.GetElement(e)->GetVertices();
      orig_mesh.GetElementEdges(e, edges, orients);
      orig_mesh.GetElementFaces(e, faces, orients);

      // One hex for each vertex of the tet.
      for (int i = 0; i < 4; ++i)
      {
        int hex_v[8];
        hex_v[0] = verts[i];
        hex_v[1] = nv_tet + edges[tet_vertex_edge_map[3 * i + 0]];
        hex_v[2] = nv_tet + nedge_tet + faces[tet_vertex_face_map[3 * i + 0]];
        hex_v[3] = nv_tet + edges[tet_vertex_edge_map[3 * i + 1]];
        hex_v[4] = nv_tet + edges[tet_vertex_edge_map[3 * i + 2]];
        hex_v[5] = nv_tet + nedge_tet + faces[tet_vertex_face_map[3 * i + 1]];
        hex_v[6] = nv_tet + nedge_tet + nface_tet + e;
        hex_v[7] = nv_tet + nedge_tet + faces[tet_vertex_face_map[3 * i + 2]];
        hex_mesh.AddHex(hex_v, orig_mesh.GetAttribute(e));
      }
    }
  }

  // Add the boundary elements.
  {
    mfem::Array<int> edges, orients;
    for (int be = 0; be < nbe_tet; ++be)
    {
      int f, o;
      const int *verts = orig_mesh.GetBdrElement(be)->GetVertices();
      orig_mesh.GetBdrElementEdges(be, edges, orients);
      orig_mesh.GetBdrElementFace(be, &f, &o);

      // One quad for each vertex of the tri.
      for (int i = 0; i < 3; ++i)
      {
        int quad_v[4];
        quad_v[0] = verts[i];
        quad_v[1] = nv_tet + edges[tri_vertex_edge_map[2 * i + 0]];
        quad_v[2] = nv_tet + nedge_tet + f;
        quad_v[3] = nv_tet + edges[tri_vertex_edge_map[2 * i + 1]];
        hex_mesh.AddBdrQuad(quad_v, orig_mesh.GetBdrAttribute(be));
      }
    }
  }

  // Finalize the hex mesh topology. The mesh will be marked for refinement later on.
  constexpr bool generate_bdr = false;
  hex_mesh.FinalizeTopology(generate_bdr);
  return hex_mesh;
}

void SplitMeshElements(std::unique_ptr<mfem::Mesh> &orig_mesh, bool make_simplex,
                       bool make_hex, bool preserve_curvature)
{
  if (!make_simplex && !make_hex)
  {
    return;
  }
  mfem::Mesh *mesh = orig_mesh.get();
  mfem::Mesh new_mesh;

  // In order to track the mapping from parent mesh element to split mesh elements for
  // high-order node interpolation, we use a unique attribute per parent mesh element. No
  // need to call mfem::Mesh::SetAttributes at the end, since we won't use the attribute
  // list data structure.
  mfem::Array<int> orig_mesh_attr;
  auto BackUpAttributes = [&orig_mesh_attr](mfem::Mesh &mesh)
  {
    orig_mesh_attr.SetSize(mesh.GetNE());
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      orig_mesh_attr[e] = mesh.GetAttribute(e);
      mesh.SetAttribute(e, 1 + e);
    }
  };

  // Convert all element types to simplices.
  if (make_simplex)
  {
    const auto element_types = mesh::CheckElements(*mesh);
    if (element_types.has_hexahedra || element_types.has_prisms ||
        element_types.has_pyramids)
    {
      MFEM_VERIFY(!mesh->Nonconforming(),
                  "Mesh element splitting is not supported for nonconforming meshes!");
      MFEM_VERIFY(
          !element_types.has_pyramids,
          "Splitting mesh elements to simplices does not support pyramid elements yet!");
      if (preserve_curvature && mesh->GetNodes() && !orig_mesh_attr.Size())
      {
        BackUpAttributes(*mesh);
      }
      int ne = mesh->GetNE();
      new_mesh = mfem::Mesh::MakeSimplicial(*mesh);
      Mpi::Print("Added {:d} elements to the mesh during conversion to simplices\n",
                 new_mesh.GetNE() - ne);
      mesh = &new_mesh;
    }
  }

  // Convert all element types to hexahedra (currently only tet-to-hex).
  if (make_hex)
  {
    const auto element_types = mesh::CheckElements(*mesh);
    if (element_types.has_simplices || element_types.has_prisms ||
        element_types.has_pyramids)
    {
      MFEM_VERIFY(!mesh->Nonconforming(),
                  "Mesh element splitting is not supported for nonconforming meshes!");
      MFEM_VERIFY(!element_types.has_prisms && !element_types.has_pyramids,
                  "Splitting mesh elements to hexahedra only supports simplex elements "
                  "(tetrahedra) for now!");
      if (preserve_curvature && mesh->GetNodes() && !orig_mesh_attr.Size())
      {
        BackUpAttributes(*mesh);
      }
      int ne = mesh->GetNE();
      new_mesh = MeshTetToHex(*mesh);
      Mpi::Print("Added {:d} elements to the mesh during conversion to hexahedra\n",
                 new_mesh.GetNE() - ne);
      mesh = &new_mesh;
    }
  }

  // Return if no modifications were made.
  if (mesh == orig_mesh.get())
  {
    return;
  }

  // The previous splitting functions remove curvature information from the new mesh. So, if
  // needed, we interpolate it onto the new mesh with GSLIB. The topology of the mesh must
  // be finalized for this to work (assumes FinalizeTopology has been called on the new
  // mesh).
  if (preserve_curvature && orig_mesh->GetNodes())
  {
    // Prepare to interpolate the grid function for high-order nodes from the old mesh to
    // the new one. This first sets up the new mesh as a linear mesh and constructs a
    // separate high-order grid function for storing the interpolated nodes.
    new_mesh.EnsureNodes();
    const mfem::GridFunction *nodes = orig_mesh->GetNodes();
    const mfem::FiniteElementSpace *fespace = nodes->FESpace();
    mfem::Ordering::Type ordering = fespace->GetOrdering();
    int order = fespace->GetMaxElementOrder();
    int sdim = orig_mesh->SpaceDimension();
    bool discont =
        (dynamic_cast<const mfem::L2_FECollection *>(fespace->FEColl()) != nullptr);
    mfem::FiniteElementSpace new_fespace(&new_mesh, fespace->FEColl(), sdim, ordering);
    mfem::GridFunction new_nodes(&new_fespace);

    mfem::Array<int> vdofs;
    mfem::Vector vals, elem_vals, xyz;
    mfem::DenseMatrix pointmat;
    int start = 0;
    for (int e = 0; e < orig_mesh->GetNE(); e++)
    {
      // Get the high-order nodes restricted to this parent element, always returned
      // byNODES.
      fespace->GetElementVDofs(e, vdofs);
      nodes->GetSubVector(vdofs, vals);

      // Find all child elements of this parent in the split mesh and restore their correct
      // original attribute. This works because child elements are added in contiguous
      // batches in the same order as the parents.
      int attr = new_mesh.GetAttribute(start);
      int end = start;
      while (end < new_mesh.GetNE() && new_mesh.GetAttribute(end) == attr)
      {
        new_mesh.SetAttribute(end++, orig_mesh_attr[e]);
      }

      // Special case: Parent element is unsplit, no interpolation is needed.
      if (end == start + 1)
      {
        new_fespace.GetElementVDofs(start, vdofs);
        new_nodes.SetSubVector(vdofs, vals);
        start = end;
        continue;
      }

      // Create a list of points at which to interpolate the high order node information.
      int npts = 0, offset = 0;
      for (int i = start; i < end; i++)
      {
        npts += new_fespace.GetFE(i)->GetNodes().GetNPoints();
      }
      xyz.SetSize(npts * sdim);
      for (int i = start; i < end; i++)
      {
        const mfem::FiniteElement &fe = *new_fespace.GetFE(i);
        mfem::ElementTransformation &T = *new_mesh.GetElementTransformation(i);
        T.Transform(fe.GetNodes(), pointmat);
        for (int d = 0; d < sdim; d++)
        {
          for (int j = 0; j < pointmat.Width(); j++)
          {
            // Use default ordering byNODES.
            xyz(d * npts + offset + j) = pointmat(d, j);
          }
        }
        offset += pointmat.Width();
      }

      // Create a single element linear mesh for interpolation.
      const mfem::Element *el = orig_mesh->GetElement(e);
      mfem::Mesh parent_mesh(orig_mesh->Dimension(), el->GetNVertices(), 1, 0,
                             orig_mesh->SpaceDimension());
      parent_mesh.AddElement(el->Duplicate(&parent_mesh));
      for (int i = 0; i < el->GetNVertices(); i++)
      {
        int v = parent_mesh.AddVertex(orig_mesh->GetVertex(el->GetVertices()[i]));
        parent_mesh.GetElement(0)->GetVertices()[i] = v;
      }
      constexpr bool generate_bdr = false;
      parent_mesh.FinalizeTopology(generate_bdr);

      // Create a grid function on the single element parent mesh with the high-order nodal
      // grid function to interpolate.
      parent_mesh.EnsureNodes();
      mfem::FiniteElementSpace parent_fespace(&parent_mesh, fespace->FEColl(), sdim,
                                              mfem::Ordering::byNODES);
      mfem::GridFunction parent_nodes(&parent_fespace);
      MFEM_ASSERT(parent_nodes.Size() == vals.Size(),
                  "Unexpected size mismatch for high-order node interpolation!");
      parent_nodes = vals;

      // Interpolate the high-order nodes grid function and copy into the new mesh nodes
      // grid function.
      vals.SetSize(npts * sdim);
      fem::InterpolateFunction(xyz, parent_nodes, vals, mfem::Ordering::byNODES);
      offset = 0;
      for (int i = start; i < end; i++)
      {
        const int elem_npts = new_fespace.GetFE(i)->GetNodes().GetNPoints();
        elem_vals.SetSize(elem_npts * sdim);
        for (int d = 0; d < sdim; d++)
        {
          for (int j = 0; j < elem_npts; j++)
          {
            // Arrange element values byNODES to align with GetElementVDofs.
            elem_vals(d * elem_npts + j) = vals(d * npts + offset + j);
          }
        }
        new_fespace.GetElementVDofs(i, vdofs);
        MFEM_ASSERT(vdofs.Size() == elem_vals.Size(),
                    "Unexpected size mismatch for high-order node interpolation!");
        new_nodes.SetSubVector(vdofs, elem_vals);
        offset += elem_npts;
      }

      // Prepare for next parent element.
      start = end;
    }
    MFEM_VERIFY(start == new_mesh.GetNE(),
                "Premature abort for high-order curvature in mesh splitting!");

    // Finally, copy the nodal grid function to the new mesh.
    new_mesh.SetCurvature(order, discont, sdim, ordering);
    MFEM_VERIFY(new_mesh.GetNodes()->Size() == new_nodes.Size(),
                "Unexpected size mismatch for nodes!");
    new_mesh.SetNodes(new_nodes);
  }
  orig_mesh = std::make_unique<mfem::Mesh>(std::move(new_mesh));  // Call move constructor
}

void ReorderMeshElements(mfem::Mesh &mesh, bool print)
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
    if (print)
    {
      Mpi::Print("Final cost: {:e}\n", best_cost);
    }
  }
  else
  {
    // (Faster) Hilbert reordering.
    mesh.GetHilbertElementOrdering(ordering);
    mesh.ReorderElements(ordering);
  }
}

std::unordered_map<int, int> GetFaceToBdrElementMap(const mfem::Mesh &mesh,
                                                    const config::BoundaryData &boundaries)
{
  std::unordered_map<int, int> face_to_be;
  face_to_be.reserve(mesh.GetNBE());
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    int f, o, e1 = -1, e2 = -1;
    mesh.GetBdrElementFace(be, &f, &o);
    if (!boundaries.periodic.empty())
    {
      int attr = mesh.GetBdrAttribute(be);
      for (const auto &data : boundaries.periodic)
      {
        const auto &da = data.donor_attributes, &ra = data.receiver_attributes;
        auto donor = std::find(da.begin(), da.end(), attr) != da.end();
        auto receiver = std::find(ra.begin(), ra.end(), attr) != ra.end();
        if (donor || receiver)
        {
          mesh.GetFaceElements(f, &e1, &e2);
          MFEM_VERIFY(e1 >= 0 && e2 >= 0,
                      "Mesh is not periodic on attribute " << attr << "!");
        }
      }
    }
    MFEM_VERIFY((e1 >= 0 && e2 >= 0) || face_to_be.find(f) == face_to_be.end(),
                "A non-periodic face (" << f << ") cannot have multiple boundary elements");
    face_to_be[f] = be;
  }
  return face_to_be;
}

std::unordered_map<int, int> CheckMesh(const mfem::Mesh &mesh,
                                       const config::BoundaryData &boundaries)
{
  // Check for:
  //   (1) Boundary elements with no prescribed boundary condition, and
  //   (2) Boundary faces which have no boundary element.
  auto bdr_marker =
      mesh::AttrToMarker(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0,
                         boundaries.attributes, true);
  std::unordered_map<int, int> face_to_be = GetFaceToBdrElementMap(mesh, boundaries);
  std::unordered_set<int> bdr_warn_list;
  int bdr_face_warn = 0;
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    int e1, e2;
    mesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0)
    {
      continue;  // Only consider true exterior faces
    }
    auto it = face_to_be.find(f);
    if (it != face_to_be.end())
    {
      int attr = mesh.GetBdrAttribute(it->second);
      if (!bdr_marker[attr - 1])
      {
        // Boundary element with no prescribed boundary condition.
        bdr_warn_list.insert(attr);
      }
    }
    else
    {
      // Boundary face with no attached boundary element.
      bdr_face_warn++;
    }
  }
  if (!bdr_warn_list.empty())
  {
    Mpi::Warning("One or more external boundary attributes has no associated boundary "
                 "condition!\n\"PMC\"/\"ZeroCharge\" condition is assumed!\n");
    utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
    Mpi::Print("\n");
  }
  if (bdr_face_warn)
  {
    Mpi::Warning("{:d} mesh faces with no associated boundary element exist on the domain "
                 "boundary!\n",
                 bdr_face_warn);
  }
  return face_to_be;
}

template <typename T>
class EdgeRefinementMesh : public mfem::Mesh
{
private:
  // Container with keys being pairs of vertex indicies (match row/column indices of v_to_v)
  // of edges which we desire to be refined.
  const T &refinement_edges;

  void MarkTetMeshForRefinement(const mfem::DSTable &v_to_v) override
  {
    // The standard tetrahedral refinement in MFEM is to mark the longest edge of each tet
    // for the first refinement. We hack this marking here in order to prioritize refinement
    // of edges which are part of internal boundary faces being marked for refinement. This
    // should hopefully limit the amount of extra refinement required to ensure conformity
    // after the marked elements are refined. Marking will then discover only the longest
    // edges, which are those within the boundary to be cracked.
    mfem::Array<mfem::real_t> lengths;
    GetEdgeLengths2(v_to_v, lengths);
    const auto min_length = 0.01 * lengths.Min();
    for (int i = 0; i < v_to_v.NumberOfRows(); i++)
    {
      for (mfem::DSTable::RowIterator it(v_to_v, i); !it; ++it)
      {
        int j = it.Column();
        if (refinement_edges.find({i, j}) == refinement_edges.end())
        {
          // "Zero" the edge lengths which do not connect vertices on the interface. Avoid
          // zero-length edges just in case.
          lengths[it.Index()] = min_length;
        }
      }
    }

    // Finish marking (see mfem::Mesh::MarkTetMeshForRefinement).
    mfem::Array<int> indices(NumOfEdges);
    std::iota(indices.begin(), indices.end(), 0);
    for (int i = 0; i < NumOfElements; i++)
    {
      if (elements[i]->GetType() == mfem::Element::TETRAHEDRON)
      {
        MFEM_ASSERT(dynamic_cast<mfem::Tetrahedron *>(elements[i]),
                    "Unexpected non-Tetrahedron element type!");
        static_cast<mfem::Tetrahedron *>(elements[i])->MarkEdge(v_to_v, lengths, indices);
      }
    }
    for (int i = 0; i < NumOfBdrElements; i++)
    {
      if (boundary[i]->GetType() == mfem::Element::TRIANGLE)
      {
        MFEM_ASSERT(dynamic_cast<mfem::Triangle *>(boundary[i]),
                    "Unexpected non-Triangle element type!");
        static_cast<mfem::Triangle *>(boundary[i])->MarkEdge(v_to_v, lengths, indices);
      }
    }
  }

public:
  EdgeRefinementMesh(mfem::Mesh &&mesh, const T &refinement_edges)
    : mfem::Mesh(std::move(mesh)), refinement_edges(refinement_edges)
  {
  }
};

template <typename T>
struct UnorderedPair
{
  T first, second;
  UnorderedPair(T first, T second) : first(first), second(second) {}
  bool operator==(const UnorderedPair &v) const
  {
    return ((v.first == first && v.second == second) ||
            (v.first == second && v.second == first));
  }
};

template <typename T>
struct UnorderedPairHasher
{
  std::size_t operator()(const UnorderedPair<T> &v) const
  {
    // Simple hash function for a pair, see https://tinyurl.com/2k4phapb.
    return std::hash<T>()(std::min(v.first, v.second)) ^
           std::hash<T>()(std::max(v.first, v.second)) << 1;
  }
};

int AddInterfaceBdrElements(const IoData &iodata, std::unique_ptr<mfem::Mesh> &orig_mesh,
                            std::unordered_map<int, int> &face_to_be)
{
  // Return if nothing to do. Otherwise, count vertices and boundary elements to add.
  if (iodata.boundaries.attributes.empty() && !iodata.model.add_bdr_elements)
  {
    return 1;  // Success
  }
  if (face_to_be.size() != static_cast<std::size_t>(orig_mesh->GetNBE()))
  {
    face_to_be = GetFaceToBdrElementMap(*orig_mesh, iodata.boundaries);
  }
  int new_nv = orig_mesh->GetNV();
  int new_nbe = orig_mesh->GetNBE();

  // Duplicate internal boundary elements from the given boundary attribute list, cracking
  // the mesh such that domain elements on either side are no longer coupled. Correctly
  // handles cracks where more than two domains intersect as well as seams where a crack
  // ends and the vertices are not duplicated.
  std::unordered_set<int> crack_bdr_elem;
  std::unordered_map<int, std::vector<std::pair<int, std::unordered_set<int>>>>
      crack_vert_duplicates;
  std::unique_ptr<mfem::Table> vert_to_elem;
  if (!iodata.boundaries.attributes.empty() && iodata.model.crack_bdr_elements)
  {
    auto crack_bdr_marker = mesh::AttrToMarker(
        orig_mesh->bdr_attributes.Size() ? orig_mesh->bdr_attributes.Max() : 0,
        iodata.boundaries.attributes, true);
    for (int be = 0; be < orig_mesh->GetNBE(); be++)
    {
      if (crack_bdr_marker[orig_mesh->GetBdrAttribute(be) - 1])
      {
        int f, o, e1, e2;
        orig_mesh->GetBdrElementFace(be, &f, &o);
        orig_mesh->GetFaceElements(f, &e1, &e2);
        if (e1 >= 0 && e2 >= 0)
        {
          crack_bdr_elem.insert(be);
        }
      }
    }
    MFEM_VERIFY(crack_bdr_elem.empty() || !orig_mesh->Nonconforming(),
                "Duplicating internal boundary elements for interior boundaries is not "
                "supported for nonconforming meshes!");

    vert_to_elem.reset(orig_mesh->GetVertexToElementTable());  // Owned by caller
    const mfem::Table &elem_to_face = orig_mesh->ElementToFaceTable();
    int new_nv_dups = 0;
    for (auto be : crack_bdr_elem)
    {
      const mfem::Element *bdr_el = orig_mesh->GetBdrElement(be);
      const int *verts = bdr_el->GetVertices();
      for (int i = 0; i < bdr_el->GetNVertices(); i++)
      {
        // Skip vertices we have already processed.
        const auto v = verts[i];
        if (crack_vert_duplicates.find(v) != crack_vert_duplicates.end())
        {
          continue;
        }

        // Collect connected components of elements connected to the vertex. Perform BFS
        // on graph of all elements connected to this vertex, where adjacencies are
        // determined by face connectivity excluding the crack faces.
        std::vector<std::unordered_set<int>> components;
        const int *elems = vert_to_elem->GetRow(v);
        std::unordered_set<int> unvisited(elems, elems + vert_to_elem->RowSize(v));
        while (!unvisited.empty())
        {
          auto &component = components.emplace_back();
          component.reserve(unvisited.size());
          std::queue<int> que;
          que.push(*unvisited.begin());
          unvisited.erase(unvisited.begin());
          while (!que.empty())
          {
            // Process the current node.
            int e = que.front();
            que.pop();
            component.insert(e);

            // Add neighbors.
            const int *faces = elem_to_face.GetRow(e);
            for (int j = 0; j < elem_to_face.RowSize(e); j++)
            {
              const auto f = faces[j];
              {
                auto it = face_to_be.find(f);
                if (it != face_to_be.end() &&
                    crack_bdr_elem.find(it->second) != crack_bdr_elem.end())
                {
                  // Skip element-element connectivities which cross the crack.
                  continue;
                }
              }
              int e1, e2;
              orig_mesh->GetFaceElements(f, &e1, &e2);
              MFEM_VERIFY(
                  e == e1 || e == e2,
                  "Unexpected face-element connectivity in internal boundary cracking!");
              int nbr = (e == e1) ? e2 : e1;
              if (nbr >= 0)
              {
                auto it = unvisited.find(nbr);
                if (it != unvisited.end())
                {
                  que.push(nbr);
                  unvisited.erase(it);
                }
              }
            }
          }
        }
        MFEM_VERIFY(
            !components.empty(),
            "No connected components found for elements adjacent to a crack vertex!");
#if defined(MFEM_DEBUG)
        {
          std::size_t visited_size = 0;
          for (const auto &component : components)
          {
            visited_size += component.size();
          }
          MFEM_VERIFY(visited_size == static_cast<std::size_t>(vert_to_elem->RowSize(v)),
                      "Failed to visit all elements in neighborhood of vertex when "
                      "counting connected components!");
        }
#endif

        // Save mapping from original vertex to duplicate vertices, and the corresponding
        // element groupings requiring renumbering. The first group doesn't need
        // renumbering so is not saved. We still keep entries for non-duplicated crack
        // vertices in the set to track them as processed, however.
        auto &vert_components = crack_vert_duplicates.try_emplace(v).first->second;
        for (auto it = components.begin() + 1; it != components.end(); ++it)
        {
          vert_components.emplace_back(-1, std::move(*it));
          new_nv_dups++;
        }
      }
    }

    // After processing all boundary elements, check if there are any elements which need
    // refinement in order to successfully decouple both sides. This happens if we have an
    // edge interior to the crack which connects to seam vertices (non-duplicated vertices
    // attached to crack boundary elements). A previous implementation of refinement
    // considered for refinement just all boundary elements with all attached vertices
    // lying on the seam, but this doesn't catch all required cases.
    if (iodata.model.refine_crack_elements)
    {
      std::unordered_map<UnorderedPair<int>, std::vector<int>, UnorderedPairHasher<int>>
          coarse_crack_edge_to_be;
      for (auto be : crack_bdr_elem)
      {
        const mfem::Element *bdr_el = orig_mesh->GetBdrElement(be);
        const int *verts = bdr_el->GetVertices();
        for (int i = 0; i < bdr_el->GetNEdges(); i++)
        {
          auto v0 = verts[bdr_el->GetEdgeVertices(i)[0]],
               v1 = verts[bdr_el->GetEdgeVertices(i)[1]];
          MFEM_ASSERT(crack_vert_duplicates.find(v0) != crack_vert_duplicates.end() &&
                          crack_vert_duplicates.find(v1) != crack_vert_duplicates.end(),
                      "Unable to locate crack vertices for an interior boundary element!");
          if (crack_vert_duplicates[v0].empty() && crack_vert_duplicates[v1].empty())
          {
            // This is a seam edge, so add the attached boundary element to a list. The
            // check for the edge being interior to the crack is indicated by visiting more
            // than once.
            auto it = coarse_crack_edge_to_be.find({v0, v1});
            auto &adjacent_be =
                (it == coarse_crack_edge_to_be.end())
                    ? coarse_crack_edge_to_be.try_emplace({v0, v1}).first->second
                    : it->second;
            adjacent_be.push_back(be);
          }
        }
      }
      for (auto it = coarse_crack_edge_to_be.begin(); it != coarse_crack_edge_to_be.end();)
      {
        // Remove all seam edges which are on the "outside" of the crack (visited only
        // once).
        if (it->second.size() == 1)
        {
          it = coarse_crack_edge_to_be.erase(it);
        }
        else
        {
          ++it;
        }
      }
      // Static reporting variables so can persist across retries.
      static int new_ne_ref = 0;
      static int new_ref_its = 0;
      if (!coarse_crack_edge_to_be.empty())
      {
        // Locally refine the mesh using conformal refinement. If necessary, convert the
        // mesh to simplices first to enable conforming refinement (this will do nothing
        // if the mesh is already a simplex mesh).
        // Note: Eventually we can implement manual conforming face refinement of pairs of
        // elements sharing a face for all element types (insert a vertex at the boundary
        // element center and connect it to all other element vertices). For now, this adds
        // complexity and making use of conformal simplex refinement seems good enough for
        // most use cases.
        int ne = orig_mesh->GetNE();
        SplitMeshElements(orig_mesh, true, false);
        if (ne != orig_mesh->GetNE())
        {
          face_to_be.clear();
          return 0;  // Mesh was converted to simplices, start over
        }
        std::unordered_map<int, int> elem_to_refine;
        for (const auto &[edge, adjacent_be] : coarse_crack_edge_to_be)
        {
          for (auto be : adjacent_be)
          {
            int f, o, e1, e2;
            orig_mesh->GetBdrElementFace(be, &f, &o);
            orig_mesh->GetFaceElements(f, &e1, &e2);
            MFEM_ASSERT(e1 >= 0 && e2 >= 0,
                        "Invalid internal boundary element connectivity!");
            elem_to_refine[e1]++;  // Value-initialized to 0 at first access
            elem_to_refine[e2]++;
          }
        }
        mfem::Array<mfem::Refinement> refinements;
        refinements.Reserve(elem_to_refine.size());
        for (const auto &[e, count] : elem_to_refine)
        {
          // Tetrahedral bisection (vs. default octasection) will result in fewer added
          // elements at the cost of a potential minor mesh quality degredation.
          refinements.Append(mfem::Refinement(e, mfem::Refinement::X));
          // refinements.Append(mfem::Refinement(e, (count > 1) ? mfem::Refinement::XY
          //                                                    : mfem::Refinement::X));
        }
        if (mesh::CheckElements(*orig_mesh).has_simplices)
        {
          // Mark tetrahedral mesh for refinement before doing local refinement. This is a
          // bit of a strange pattern to override the standard conforming refinement of the
          // mfem::Mesh class. We want to implement our own edge marking of the tetrahedra,
          // so we move the mesh to a constructed derived class object, mark it, and then
          // move assign it to the original base class object before refining. All of these
          // moves should be cheap without any extra memory allocation. Also, we mark the
          // mesh every time to ensure multiple rounds of refinement target the interior
          // boundary (we don't care about preserving the refinement hierarchy).
          constexpr bool refine = true, fix_orientation = false;
          EdgeRefinementMesh ref_mesh(std::move(*orig_mesh), coarse_crack_edge_to_be);
          ref_mesh.Finalize(refine, fix_orientation);
          *orig_mesh = std::move(ref_mesh);
        }
        orig_mesh->GeneralRefinement(refinements, 0);
        new_ne_ref += orig_mesh->GetNE() - ne;
        new_ref_its++;
        face_to_be.clear();
        return 0;  // Mesh was refined (conformally), start over
      }
      else if (new_ne_ref > 0)
      {
        Mpi::Print(
            "Added {:d} elements in {:d} iterations of local bisection for under-resolved "
            "interior boundaries\n",
            new_ne_ref, new_ref_its);
      }
    }

    new_nv += new_nv_dups;
    new_nbe += crack_bdr_elem.size();
    if (crack_bdr_elem.size() > 0)
    {
      Mpi::Print("Added {:d} duplicate vertices for interior boundaries in the mesh\n",
                 new_nv_dups);
      Mpi::Print(
          "Added {:d} duplicate boundary elements for interior boundaries in the mesh\n",
          crack_bdr_elem.size());
    }
  }

  // Add new boundary elements at material interfaces or on the exterior boundary of the
  // simulation domain, if there is not already a boundary element present.
  std::unordered_map<int, int> new_face_bdr_elem;
  if (iodata.model.add_bdr_elements)
  {
    int new_nbe_ext = 0, new_nbe_int = 0;
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      // Skip all faces which already have an associated boundary element (this includes
      // any boundary elements which were duplicated during cracking in the previous step).
      if (face_to_be.find(f) != face_to_be.end())
      {
        continue;
      }
      int e1, e2;
      orig_mesh->GetFaceElements(f, &e1, &e2);
      if (e1 < 0 || e2 < 0)
      {
        if constexpr (false)
        {
          Mpi::Print("Adding exterior boundary element!\n");
        }
        new_face_bdr_elem[f] = 1;
        new_nbe_ext++;
      }
      else if (orig_mesh->GetAttribute(e1) != orig_mesh->GetAttribute(e2))
      {
        if constexpr (false)
        {
          Mpi::Print("Adding material interface boundary element!\n");
        }
        new_face_bdr_elem[f] = 1;
        new_nbe_int++;
      }
    }
    MFEM_VERIFY(new_nbe_ext + new_nbe_int == 0 || !orig_mesh->Nonconforming(),
                "Adding material interface boundary elements is not supported for "
                "nonconforming meshes!");

    new_nbe += (new_nbe_ext + new_nbe_int);
    if (new_nbe_ext > 0)
    {
      Mpi::Print("Added {:d} boundary elements for exterior boundaries to the mesh\n",
                 new_nbe_ext);
    }
    if (new_nbe_int > 0)
    {
      Mpi::Print("Added {:d} boundary elements for material interfaces to the mesh\n",
                 new_nbe_int);
    }
  }

  // Create the new mesh. We can't just add the new vertices and boundary elements to the
  // original mesh, because we need to keep it around in order to transfer the high-order
  // nodes information to the new mesh.
  if (new_nv == orig_mesh->GetNV() && new_nbe == orig_mesh->GetNBE())
  {
    return 1;  // Success
  }
  auto new_mesh =
      std::make_unique<mfem::Mesh>(orig_mesh->Dimension(), new_nv, orig_mesh->GetNE(),
                                   new_nbe, orig_mesh->SpaceDimension());

  // Copy vertices and domain and boundary elements.
  for (int v = 0; v < orig_mesh->GetNV(); v++)
  {
    new_mesh->AddVertex(orig_mesh->GetVertex(v));
  }
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    mfem::Element *el = orig_mesh->GetElement(e)->Duplicate(new_mesh.get());
    new_mesh->AddElement(el);
  }
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    mfem::Element *bdr_el = orig_mesh->GetBdrElement(be)->Duplicate(new_mesh.get());
    new_mesh->AddBdrElement(bdr_el);
  }

  // Add duplicated vertices from interior boundary cracking, renumber the vertices of
  // domain and boundary elements to tear the mesh, and add new crack boundary elements.
  if (!iodata.boundaries.attributes.empty() && !crack_bdr_elem.empty())
  {
    // Add duplicate vertices. We assign the vertex number of the duplicated vertex in order
    // to update the element connectivities in the next step.
    for (auto &[orig_v, vert_components] : crack_vert_duplicates)
    {
      for (auto &[dup_v, component] : vert_components)
      {
        dup_v = new_mesh->AddVertex(orig_mesh->GetVertex(orig_v));
      }
    }

    // Renumber the duplicated vertex in the domain elements.
    for (const auto &[orig_v, vert_components] : crack_vert_duplicates)
    {
      if (vert_components.empty())
      {
        continue;  // Can skip vertices which were not duplicated
      }
      const int *elems = vert_to_elem->GetRow(orig_v);
      for (int i = 0; i < vert_to_elem->RowSize(orig_v); i++)
      {
        // Find vertex in the element.
        const auto e = elems[i];
        mfem::Element *el = new_mesh->GetElement(e);
        int *verts = el->GetVertices(), j;
        for (j = 0; j < el->GetNVertices(); j++)
        {
          if (verts[j] == orig_v)
          {
            break;
          }
        }
        MFEM_VERIFY(j < el->GetNVertices(), "Unable to locate vertex in element!");

        // Find the correct duplicate for this vertex. It's OK if the element is not in
        // any of the connected components, this indicates that it keeps the original
        // vertex and its connectivity is unmodified.
        for (const auto &[dup_v, component] : vert_components)
        {
          if (component.find(e) != component.end())
          {
            verts[j] = dup_v;
            break;
          }
        }
      }
    }

    // Finally, we insert the new duplicate boundary elements for the crack interface and
    // also renumber the original boundary elements. To renumber the original boundary
    // elements in the mesh, we use the updated vertex connectivity from the torn elements
    // in the new mesh (done above).
    const mfem::Table &elem_to_face = orig_mesh->ElementToFaceTable();
    for (int be = 0; be < orig_mesh->GetNBE(); be++)
    {
      // Whether on the crack or not, we renumber the boundary element vertices as needed
      // based on the neighboring element. For non-crack boundary elements, both
      // neighboring elements must be part of the same connected component. First we find
      // the index of the face in the old element, which should match the new element.
      int f, o, e1, e2;
      orig_mesh->GetBdrElementFace(be, &f, &o);
      orig_mesh->GetFaceElements(f, &e1, &e2);
      MFEM_VERIFY(e1 >= 0, "Boundary element with no attached elements!");
      const int *faces = elem_to_face.GetRow(e1);
      int i;
      for (i = 0; i < elem_to_face.RowSize(e1); i++)
      {
        if (faces[i] == f)
        {
          break;
        }
      }
      MFEM_VERIFY(i < elem_to_face.RowSize(e1), "Unable to locate face in element!");

      // Update the boundary element vertices.
      mfem::Element *bdr_el = new_mesh->GetBdrElement(be);
      const mfem::Element *el = new_mesh->GetElement(e1);
      for (int j = 0; j < bdr_el->GetNVertices(); j++)
      {
        bdr_el->GetVertices()[j] = el->GetVertices()[el->GetFaceVertices(i)[j]];
      }

      // Add the duplicate boundary element for boundary elements on the crack.
      if (crack_bdr_elem.find(be) != crack_bdr_elem.end())
      {
        faces = elem_to_face.GetRow(e2);
        for (i = 0; i < elem_to_face.RowSize(e2); i++)
        {
          if (faces[i] == f)
          {
            break;
          }
        }
        MFEM_VERIFY(i < elem_to_face.RowSize(e2), "Unable to locate face in element!");

        // Add the interface boundary element attached to element 2 (the other part of the
        // pair has been attached to element 1 in the previous step).
        bdr_el = bdr_el->Duplicate(new_mesh.get());
        el = new_mesh->GetElement(e2);
        for (int j = 0; j < bdr_el->GetNVertices(); j++)
        {
          bdr_el->GetVertices()[j] = el->GetVertices()[el->GetFaceVertices(i)[j]];
        }
        new_mesh->AddBdrElement(bdr_el);
      }
    }
  }

  // Add new boundary elements.
  if (iodata.model.add_bdr_elements && !new_face_bdr_elem.empty())
  {
    // Some (1-based) boundary attributes may be empty since they were removed from the
    // original mesh, but to keep attributes the same as config file we don't compress the
    // list.
    const mfem::Table &elem_to_face = orig_mesh->ElementToFaceTable();
    int bdr_attr_max =
        orig_mesh->bdr_attributes.Size() ? orig_mesh->bdr_attributes.Max() : 0;
    for (int f = 0; f < orig_mesh->GetNumFaces(); f++)
    {
      if (new_face_bdr_elem[f] > 0)
      {
        // Assign new unique attribute based on attached elements. Save so that the
        // attributes of e1 and e2 can be easily referenced using the new attribute. Since
        // attributes are in 1-based indexing, a, b > 0. See also
        // https://en.wikipedia.org/wiki/Pairing_function.
        int e1, e2, a = 0, b = 0;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        if (e1 >= 0 && e2 >= 0)
        {
          a = std::max(orig_mesh->GetAttribute(e1), orig_mesh->GetAttribute(e2));
          b = (a == orig_mesh->GetAttribute(e1)) ? orig_mesh->GetAttribute(e2)
                                                 : orig_mesh->GetAttribute(e1);
        }
        else  // e1 >= 0
        {
          a = orig_mesh->GetAttribute(e1);
          b = 0;
        }
        MFEM_VERIFY(a + b > 0, "Invalid new boundary element attribute!");
        long long int l_new_attr =
            bdr_attr_max + (((a + b) * (long long int)(a + b + 1)) / 2) + a;
        int new_attr = mfem::internal::to_int(l_new_attr);  // At least bdr_attr_max + 1

        // Add the boundary elements with the new boundary attribute. The element vertices
        // may have been renumbered in the new mesh, so the new face is not necessarily
        // just a duplicate of the old one. First we find the index of the face in the old
        // element, which should match the new element.
        const int *faces = elem_to_face.GetRow(e1);
        int i;
        for (i = 0; i < elem_to_face.RowSize(e1); i++)
        {
          if (faces[i] == f)
          {
            break;
          }
        }
        MFEM_VERIFY(i < elem_to_face.RowSize(e1), "Unable to locate face in element!");

        // Now add the boundary element(s).
        mfem::Element *bdr_el = orig_mesh->GetFace(f)->Duplicate(new_mesh.get());
        bdr_el->SetAttribute(new_attr);
        const mfem::Element *el = new_mesh->GetElement(e1);
        for (int j = 0; j < bdr_el->GetNVertices(); j++)
        {
          bdr_el->GetVertices()[j] = el->GetVertices()[el->GetFaceVertices(i)[j]];
        }
        new_mesh->AddBdrElement(bdr_el);
        if constexpr (false)
        {
          Mpi::Print(
              "Adding boundary element with attribute {:d} from elements {:d} and {:d}\n",
              new_attr, a, b);
        }
        if (new_face_bdr_elem[f] > 1)
        {
          // Flip order of vertices to reverse normal direction of the second added element.
          bdr_el = bdr_el->Duplicate(new_mesh.get());
          std::reverse(bdr_el->GetVertices(),
                       bdr_el->GetVertices() + bdr_el->GetNVertices());
          new_mesh->AddBdrElement(bdr_el);
          if constexpr (false)
          {
            Mpi::Print("Adding second boundary element with attribute {:d} from elements "
                       "{:d} and {:d}\n",
                       new_attr, a, b);
          }
        }
      }
    }
  }

  // Finalize the new mesh topology, and copy mesh curvature information if needed. This
  // copies the nodes over correctly accounting for the element topology changes (the number
  // of elements in the mesh has not changed, just their connectivity has).
  constexpr bool generate_bdr = false;
  new_mesh->FinalizeTopology(generate_bdr);
  if (orig_mesh->GetNodes())
  {
    TransferHighOrderNodes(*orig_mesh, *new_mesh);
  }

  // If we have added cracks for interior boundary elements, apply a very very small
  // perturbation to separate the duplicated boundary elements on either side and prevent
  // them from lying exactly on top of each other. This is mostly just for visualization
  // and can be increased in magnitude for debugging.
  if (!iodata.boundaries.attributes.empty() && !crack_bdr_elem.empty() &&
      iodata.model.crack_displ_factor > 0.0)
  {
    // mfem::Mesh::MoveNodes expects byNODES ordering when using vertices.
    mfem::GridFunction *nodes = new_mesh->GetNodes();
    mfem::Ordering::Type ordering =
        nodes ? nodes->FESpace()->GetOrdering() : mfem::Ordering::byNODES;
    int sdim = new_mesh->SpaceDimension();
    int nv = nodes ? nodes->Size() / sdim : new_mesh->GetNV();
    auto Index = [ordering, sdim, nv](int v, int d)
    { return (ordering == mfem::Ordering::byVDIM) ? sdim * v + d : d * nv + v; };
    mfem::Vector normal(sdim);
    mfem::IsoparametricTransformation T;
    mfem::Array<int> dofs;

    // Compute the displacement as the average normal of the attached boundary elements.
    mfem::Vector displacements(nv * sdim);
    displacements = 0.0;
    double h_min = mfem::infinity();
    const mfem::Table &elem_to_face = orig_mesh->ElementToFaceTable();
    const mfem::Table &new_elem_to_face = new_mesh->ElementToFaceTable();
    for (auto be : crack_bdr_elem)
    {
      // Get the neighboring elements (same indices in the old and new mesh).
      int f, o, e1, e2;
      orig_mesh->GetBdrElementFace(be, &f, &o);
      orig_mesh->GetFaceElements(f, &e1, &e2);

      // Perturb both new boundary elements in opposite directions.
      for (auto e : {e1, e2})
      {
        // Find the index of the face in the old element, which matches the new element, so
        // we can get the list of all vertices or nodes to perturb.
        const int *faces = elem_to_face.GetRow(e);
        int i;
        for (i = 0; i < elem_to_face.RowSize(e); i++)
        {
          if (faces[i] == f)
          {
            break;
          }
        }
        MFEM_VERIFY(i < elem_to_face.RowSize(e), "Unable to locate face in element!");

        // Compute the element normal, oriented to point outward from element 1 initially.
        int new_f = new_elem_to_face.GetRow(e)[i];
        if (e == e1)
        {
          new_mesh->GetFaceTransformation(new_f, &T);
          const mfem::IntegrationPoint &ip =
              mfem::Geometries.GetCenter(T.GetGeometryType());
          T.SetIntPoint(&ip);
          mfem::CalcOrtho(T.Jacobian(), normal);
          double s = normal.Norml2();
          h_min = std::min(h_min, std::sqrt(s));
          normal /= -s;  // We could also area-weight the average normal
        }
        else  // e == e2
        {
          normal *= -1.0;
        }

        // For all "nodes" associated with this crack face, update the direction of their
        // displacements.
        auto NodeUpdate = [&](int v)
        {
          for (int d = 0; d < sdim; d++)
          {
            const int idx = Index(v, d);
            displacements(idx) += normal(d);
          }
        };
        if (nodes)
        {
          nodes->FESpace()->GetFaceDofs(new_f, dofs);
          for (int j = 0; j < dofs.Size(); j++)
          {
            NodeUpdate(dofs[j]);
          }
        }
        else
        {
          const mfem::Element *el = new_mesh->GetElement(e);
          for (int j = 0; j < el->GetNFaceVertices(i); j++)
          {
            NodeUpdate(el->GetVertices()[el->GetFaceVertices(i)[j]]);
          }
        }
      }
    }
    for (int v = 0; v < nv; v++)
    {
      double s = 0.0;
      for (int d = 0; d < sdim; d++)
      {
        const int idx = Index(v, d);
        s += displacements(idx) * displacements(idx);
      }
      if (s > 0.0)
      {
        s = std::sqrt(s);
        for (int d = 0; d < sdim; d++)
        {
          const int idx = Index(v, d);
          displacements(idx) /= s;
        }
      }
    }

    // Scale and apply the displacements. We don't need to do anything special to constrain
    // the displacements at seam vertices (and associated high-order nodes on seam edges) to
    // to zero, because the normals from both sides will average out to zero.
    displacements *= (iodata.model.crack_displ_factor * h_min /
                      (nodes ? nodes->FESpace()->GetMaxElementOrder() : 1));
    new_mesh->MoveNodes(displacements);
  }

  orig_mesh = std::move(new_mesh);
  return 1;  // Success
}

std::unique_ptr<int[]> GetMeshPartitioning(const mfem::Mesh &mesh, int size,
                                           const std::string &part_file, bool print)
{
  MFEM_VERIFY(size <= mesh.GetNE(), "Mesh partitioning must have parts <= mesh elements ("
                                        << size << " vs. " << mesh.GetNE() << ")!");
  if (part_file.length() == 0)
  {
    const int part_method = 1;
    std::unique_ptr<int[]> partitioning(
        const_cast<mfem::Mesh &>(mesh).GeneratePartitioning(size, part_method));
    if (print)
    {
      Mpi::Print("Finished partitioning mesh into {:d} subdomain{}\n", size,
                 (size > 1) ? "s" : "");
    }
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
  std::ifstream part_ifs(part_file);
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
  if (print)
  {
    Mpi::Print("Read mesh partitioning into {:d} subdomain{} from disk\n", size,
               (size > 1) ? "s" : "");
  }
  return partitioning;
}

std::unique_ptr<mfem::ParMesh> DistributeMesh(MPI_Comm comm,
                                              std::unique_ptr<mfem::Mesh> &smesh,
                                              const int *partitioning,
                                              const std::string &output_dir)
{
  // Take a serial mesh and partitioning on the root process and construct the global
  // parallel mesh. For now, prefer the MPI-based version to the file IO one. When
  // constructing the ParMesh, we mark for refinement since refinement flags are not copied
  // from the serial mesh. Beware that mfem::ParMesh constructor argument order is not the
  // same as mfem::Mesh! Each processor's component gets sent as a byte string.
  constexpr bool generate_edges = false, refine = true, fix_orientation = false;
  std::unique_ptr<mfem::ParMesh> pmesh;
  if (Mpi::Root(comm))
  {
    mfem::MeshPartitioner partitioner(*smesh, Mpi::Size(comm),
                                      const_cast<int *>(partitioning));
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
    smesh.reset();
    std::istringstream fi(so[0]);  // This is never compressed
    pmesh =
        std::make_unique<mfem::ParMesh>(comm, fi, refine, generate_edges, fix_orientation);
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
    pmesh =
        std::make_unique<mfem::ParMesh>(comm, fi, refine, generate_edges, fix_orientation);
  }
  return pmesh;
}

void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &pmesh)
{
  // Write the parallel mesh to a stream as a serial mesh, then read back in and partition
  // using METIS.
  MPI_Comm comm = pmesh->GetComm();
  constexpr bool generate_edges = false, generate_bdr = false, refine = true,
                 fix_orientation = false;
  std::unique_ptr<mfem::Mesh> smesh;
  if constexpr (false)
  {
    // Write the serial mesh to a stream and read that through the Mesh constructor.
    std::ostringstream fo(std::stringstream::out);
    // fo << std::fixed;
    fo << std::scientific;
    fo.precision(MSH_FLT_PRECISION);
    pmesh->PrintAsSerial(fo);
    if (Mpi::Root(comm))
    {
      smesh = std::make_unique<mfem::Mesh>(fo, generate_edges, refine, fix_orientation);
    }
  }
  else
  {
    // Directly ingest the generated Mesh and release the no longer needed memory.
    smesh = std::make_unique<mfem::Mesh>(pmesh->GetSerialMesh(0));
    if (Mpi::Root(comm))
    {
      smesh->FinalizeTopology(generate_bdr);
      smesh->Finalize(refine, fix_orientation);
    }
    else
    {
      smesh.reset();
    }
  }

  // (Re)-construct the parallel mesh.
  std::unique_ptr<int[]> partitioning;
  if (Mpi::Root(comm))
  {
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), "", false);
  }
  pmesh = DistributeMesh(comm, smesh, partitioning.get());
}

}  // namespace

}  // namespace palace
