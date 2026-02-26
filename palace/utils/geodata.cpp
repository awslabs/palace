// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"
#include "geodata_impl.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Eigen/Dense>
#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "fem/interpolator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"
#include "utils/filesystem.hpp"
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
void CleanMesh(std::unique_ptr<mfem::Mesh> &, const std::vector<int> &);

// Create a new mesh by splitting all elements of the mesh into simplices or hexes
// (using tet-to-hex). Optionally preserves curvature of the original mesh by interpolating
// the high-order nodes with GSLIB.
void SplitMeshElements(std::unique_ptr<mfem::Mesh> &, bool, bool);

// Optionally reorder mesh elements based on MFEM's internal reordering tools for improved
// cache usage.
void ReorderMeshElements(mfem::Mesh &, bool = true);

// Check that mesh boundary conditions are given for external boundaries.
std::unordered_map<int, int> CheckMesh(const mfem::Mesh &, const config::BoundaryData &);

// Adding boundary elements for material interfaces and exterior boundaries, and "crack"
// desired internal boundary elements to disconnect the elements on either side.
int AddInterfaceBdrElements(IoData &, std::unique_ptr<mfem::Mesh> &,
                            std::unordered_map<int, int> &, MPI_Comm comm);

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

std::unique_ptr<mfem::ParMesh> ReadMesh(IoData &iodata, MPI_Comm comm)
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

  const bool use_mesh_partitioner = [&]()
  {
    // Root must load the mesh to discover if nonconformal, as a previously adapted mesh
    // might be reused for nonadaptive simulations.
    BlockTimer bt(Timer::IO);
    bool use_mesh_partitioner = !use_amr || !refinement.nonconformal;
    if (Mpi::Root(comm))
    {
      smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature, iodata.boundaries);
      use_mesh_partitioner &= smesh->Conforming();  // The initial mesh must be conformal
    }
    Mpi::Broadcast(1, &use_mesh_partitioner, 0, comm);
    return use_mesh_partitioner;
  }();

  MPI_Comm node_comm;
  if (!use_mesh_partitioner)
  {
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
  }

  {
    BlockTimer bt1(Timer::IO);
    if (!use_mesh_partitioner && Mpi::Root(node_comm) && !Mpi::Root(comm))
    {
      // Only one process per node reads the serial mesh, if not using mesh partitioner.
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
    MFEM_VERIFY(
        smesh->Conforming() || !use_amr || refinement.nonconformal,
        "The provided mesh is nonconformal, only nonconformal AMR can be performed!");

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
    if (smesh->Conforming())
    {
      auto face_to_be = CheckMesh(*smesh, iodata.boundaries);

      // Add new boundary elements for material interfaces if not present (with new unique
      // boundary attributes). Also duplicate internal boundary elements associated with
      // cracks if desired.
      if ((iodata.model.crack_bdr_elements || iodata.model.add_bdr_elements))
      {
        // Split all internal (non periodic) boundary elements for boundary attributes where
        // BC are applied (not just postprocessing).
        while (AddInterfaceBdrElements(iodata, smesh, face_to_be, comm) != 1)
        {
          // May require multiple calls due to early exit/retry approach.
        }
      }
    }
    else
    {
      Mpi::Warning("{} is a nonconformal mesh, assumed from previous AMR iteration.\n"
                   "Skipping mesh modification preprocessing steps!\n\n",
                   iodata.model.mesh);
    }

    // Finally, finalize the serial mesh. Mark tetrahedral meshes for refinement. There
    // should be no need to fix orientation as this was done during initial mesh loading
    // from disk.
    constexpr bool refine = true, fix_orientation = false;
    smesh->Finalize(refine, fix_orientation);

    // Generate the mesh partitioning.
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), iodata.model.partitioning);
  }

  // Broadcast cracked boundary attributes to other ranks.
  if ((iodata.model.crack_bdr_elements || iodata.model.add_bdr_elements))
  {
    int size = iodata.boundaries.cracked_attributes.size();
    Mpi::Broadcast(1, &size, 0, comm);
    std::vector<int> data;
    if (Mpi::Root(comm))
    {
      data.assign(iodata.boundaries.cracked_attributes.begin(),
                  iodata.boundaries.cracked_attributes.end());
    }
    else
    {
      data.resize(size);
    }
    Mpi::Broadcast(size, data.data(), 0, comm);

    if (!Mpi::Root(comm))
    {
      iodata.boundaries.cracked_attributes.clear();
      iodata.boundaries.cracked_attributes.insert(data.begin(), data.end());
    }
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
  const double Lc = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
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
  const int ne = 4 * ne_tet;    // 4 hex per tet
  const int nbe = 3 * nbe_tet;  // 3 square per tri
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

  // Connectivity of tetrahedron vertices to the edges. The vertices of the new mesh are
  // ordered so that the original tet vertices are first, then the vertices splitting each
  // edge, then the vertices at the center of each triangle face, then the center of the
  // tet. Thus the edge/face/element numbers can be used to index into the new array of
  // vertices, and the element local edge/face can be used to extract the global edge/face
  // index, and thus the corresponding vertex.
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

  // All elements have now been added, can construct the higher order field.
  if (orig_mesh.GetNodes())
  {
    hex_mesh.EnsureNodes();
    // Higher order associated to vertices are unchanged, and those for
    // previously existing edges. DOFs associated to new elements need to be set.
    const int sdim = orig_mesh.SpaceDimension();
    auto *orig_fespace = orig_mesh.GetNodes()->FESpace();
    hex_mesh.SetCurvature(orig_fespace->GetMaxElementOrder(), orig_fespace->IsDGSpace(),
                          orig_mesh.SpaceDimension(), orig_fespace->GetOrdering());

    // Need to convert the hexahedra local coordinate system into the parent tetrahedra
    // system. Each hexahedra spans a different set of the tet's reference coordinates. To
    // convert between, define the reference coordinate locations of each of the vertices
    // the hexahedra will use, then perform trilinear interpolation in the reference space.

    auto [vert_loc, edge_loc, face_loc] = []()
    {
      std::array<std::array<double, 3>, 4> vert_loc{};
      vert_loc[0] = {0.0, 0.0, 0.0};
      vert_loc[1] = {1.0, 0.0, 0.0};
      vert_loc[2] = {0.0, 1.0, 0.0};
      vert_loc[3] = {0.0, 0.0, 1.0};
      std::array<std::array<double, 3>, 6> edge_loc{};
      edge_loc[0] = {0.5, 0.0, 0.0};
      edge_loc[1] = {0.0, 0.5, 0.0};
      edge_loc[2] = {0.0, 0.0, 0.5};
      edge_loc[3] = {0.5, 0.5, 0.0};
      edge_loc[4] = {0.5, 0.0, 0.5};
      edge_loc[5] = {0.0, 0.5, 0.5};
      std::array<std::array<double, 3>, 6> face_loc{};
      face_loc[0] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
      face_loc[1] = {0.0, 1.0 / 3, 1.0 / 3};
      face_loc[2] = {1.0 / 3, 0.0, 1.0 / 3};
      face_loc[3] = {1.0 / 3, 1.0 / 3, 0.0};
      return std::make_tuple(vert_loc, edge_loc, face_loc);
    }();
    std::array<double, 3> centroid{{0.25, 0.25, 0.25}};

    // We assume the Nodes field is of a single order, and there is a single tet originally.
    // The nodes within the reference hex and parent tet are always the same, so we use the
    // typical FE. We then exploit the fact the map between reference spaces is always
    // linear, and construct the transformation explicitly.
    const auto *orig_FE = orig_mesh.GetNodes()->FESpace()->GetTypicalFE();
    const auto *child_FE = hex_mesh.GetNodes()->FESpace()->GetTypicalFE();
    // Original shape function (i), at new element nodes (j), for each new element (k).
    mfem::DenseTensor shape(orig_FE->GetDof(), child_FE->GetDof(), 4);
    mfem::Vector col;  // For slicing into matrices within shape
    for (int i = 0; i < 4; i++)
    {
      // Collect the vertices of the new hex within the tet.
      std::array<std::array<double, 3>, 8> hex_verts;
      hex_verts[0] = vert_loc[i];
      hex_verts[1] = edge_loc[tet_vertex_edge_map[3 * i + 0]];
      hex_verts[2] = face_loc[tet_vertex_face_map[3 * i + 0]];
      hex_verts[3] = edge_loc[tet_vertex_edge_map[3 * i + 1]];
      hex_verts[4] = edge_loc[tet_vertex_edge_map[3 * i + 2]];
      hex_verts[5] = face_loc[tet_vertex_face_map[3 * i + 1]];
      hex_verts[6] = centroid;
      hex_verts[7] = face_loc[tet_vertex_face_map[3 * i + 2]];
      for (int j = 0; j < child_FE->GetNodes().Size(); j++)
      {
        const auto &cn = child_FE->GetNodes()[j];
        mfem::IntegrationPoint cn_in_orig;

        // Perform trilinear interpolation from (u,v,w) the unit ref coords in the new hex,
        // and the corresponding nodes in the containing tet.
        // clang-format off
        // x component
        cn_in_orig.x =
          hex_verts[0][0] * (1-cn.x) * (1-cn.y) * (1-cn.z) +
          hex_verts[1][0] * cn.x     * (1-cn.y) * (1-cn.z) +
          hex_verts[2][0] * cn.x     * cn.y     * (1-cn.z) +
          hex_verts[3][0] * (1-cn.x) * cn.y     * (1-cn.z) +
          hex_verts[4][0] * (1-cn.x) * (1-cn.y) * cn.z     +
          hex_verts[5][0] * cn.x     * (1-cn.y) * cn.z     +
          hex_verts[6][0] * cn.x     * cn.y     * cn.z     +
          hex_verts[7][0] * (1-cn.x) * cn.y     * cn.z;

        // y component
        cn_in_orig.y =
          hex_verts[0][1] * (1-cn.x) * (1-cn.y) * (1-cn.z) +
          hex_verts[1][1] * cn.x     * (1-cn.y) * (1-cn.z) +
          hex_verts[2][1] * cn.x     * cn.y     * (1-cn.z) +
          hex_verts[3][1] * (1-cn.x) * cn.y     * (1-cn.z) +
          hex_verts[4][1] * (1-cn.x) * (1-cn.y) * cn.z     +
          hex_verts[5][1] * cn.x     * (1-cn.y) * cn.z     +
          hex_verts[6][1] * cn.x     * cn.y     * cn.z     +
          hex_verts[7][1] * (1-cn.x) * cn.y     * cn.z;

        // z component
        cn_in_orig.z =
          hex_verts[0][2] * (1-cn.x) * (1-cn.y) * (1-cn.z) +
          hex_verts[1][2] * cn.x     * (1-cn.y) * (1-cn.z) +
          hex_verts[2][2] * cn.x     * cn.y     * (1-cn.z) +
          hex_verts[3][2] * (1-cn.x) * cn.y     * (1-cn.z) +
          hex_verts[4][2] * (1-cn.x) * (1-cn.y) * cn.z     +
          hex_verts[5][2] * cn.x     * (1-cn.y) * cn.z     +
          hex_verts[6][2] * cn.x     * cn.y     * cn.z     +
          hex_verts[7][2] * (1-cn.x) * cn.y     * cn.z;
        // clang-format on
        shape(i).GetColumnReference(j, col);
        orig_FE->CalcShape(cn_in_orig, col);
      }
    }

    // Each submatrix of shape tensor now encodes the reference coordinates of each hex
    // within the containing tet. Extracting the specific element dof values, and applying
    // to the correct shape slice will now give the requisite higher order dofs evaluated at
    // the refined elements nodes.
    mfem::Array<int> hex_dofs;
    mfem::DenseMatrix point_matrix(child_FE->GetDof(), sdim);  // nnode_child x sdim
    mfem::Vector dof_vals(orig_FE->GetDof() * sdim);
    mfem::DenseMatrix dof_vals_mat(dof_vals.GetData(), orig_FE->GetDof(), sdim);
    for (int e = 0; e < ne_tet; ++e)
    {
      // Returns byNODES no matter what, because FiniteElementSpace::GetElementVDofs does.
      // Matches the GetElementVDofs call below, which similarly always uses byNODES.
      orig_mesh.GetNodes()->GetElementDofValues(e, dof_vals);
      for (int i = 0; i < 4; i++)
      {
        // shape(i) : orig_FE->GetDof() x hex_FE->GetDof()
        // dof_vals_mat : orig_FE->GetDof() x sdim
        // point_matrix : child_FE->GetDof() x sdim
        MultAtB(shape(i), dof_vals_mat, point_matrix);
        hex_mesh.GetNodes()->FESpace()->GetElementVDofs(4 * e + i, hex_dofs);
        hex_mesh.GetNodes()->SetSubVector(hex_dofs, point_matrix.GetData());
      }
    }
  }

  return hex_mesh;
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

std::vector<mfem::Geometry::Type> ElementTypeInfo::GetGeomTypes(int dim) const
{
  std::vector<mfem::Geometry::Type> geom_types;
  if (dim == 2)
  {
    if (has_simplices)
    {
      geom_types.push_back(mfem::Geometry::TRIANGLE);
    }
    if (has_hexahedra)
    {
      geom_types.push_back(mfem::Geometry::SQUARE);
    }
  }
  else
  {
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
      MFEM_VERIFY(marker[attr - 1] == 0, "Repeated attribute in attribute list!");
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
  const int dim = Dim();
  if (dim == 3)
  {
    return 4.0 *
           CVector3dMap(axes.GetColumn(0)).cross(CVector3dMap(axes.GetColumn(1))).norm();
  }
  // 2D: area of the parallelogram spanned by the two axis vectors (2D cross product).
  const double *a0 = axes.GetColumn(0);
  const double *a1 = axes.GetColumn(1);
  return 4.0 * std::abs(a0[0] * a1[1] - a0[1] * a1[0]);
}

double BoundingBox::Volume() const
{
  if (Dim() < 3 || planar)
  {
    return 0.0;
  }
  return 2.0 * CVector3dMap(axes.GetColumn(2)).norm() * Area();
}

mfem::DenseMatrix BoundingBox::Normals() const
{
  const int dim = Dim();
  mfem::DenseMatrix normals(axes);
  for (int i = 0; i < dim; i++)
  {
    mfem::Vector col(normals.GetColumn(i), normals.Height());
    double nrm = col.Norml2();
    if (nrm > 0.0)
    {
      col /= nrm;
    }
  }
  return normals;
}

mfem::Vector BoundingBox::Lengths() const
{
  const int dim = Dim();
  const int h = axes.Height();
  mfem::Vector lengths(dim);
  for (int i = 0; i < dim; i++)
  {
    const double *col = axes.GetColumn(i);
    double nrm = 0.0;
    for (int j = 0; j < h; j++)
    {
      nrm += col[j] * col[j];
    }
    lengths(i) = 2.0 * std::sqrt(nrm);
  }
  return lengths;
}

mfem::Vector BoundingBox::Deviations(const mfem::Vector &direction) const
{
  const int dim = Dim();
  const int h = axes.Height();
  MFEM_VERIFY(direction.Size() == dim,
              "Direction vector size must match bounding box dimension!");
  double dir_norm = direction.Norml2();
  mfem::Vector deviation_deg(dim);
  for (int i = 0; i < dim; i++)
  {
    const double *col = axes.GetColumn(i);
    double ax_norm = 0.0, dot = 0.0;
    for (int j = 0; j < h; j++)
    {
      ax_norm += col[j] * col[j];
      dot += direction(j) * col[j];
    }
    ax_norm = std::sqrt(ax_norm);
    double cosine = (dir_norm > 0.0 && ax_norm > 0.0) ? dot / (dir_norm * ax_norm) : 0.0;
    deviation_deg(i) = std::acos(std::min(1.0, std::abs(cosine))) * (180.0 / M_PI);
  }
  return deviation_deg;
}

namespace
{

// Trim a 3D BoundingBox to the given spatial dimension by extracting the leading
// dim x dim submatrix from axes and first dim entries of center.
void TrimBoundingBox(BoundingBox &box, int sdim)
{
  if (sdim < 3 && box.Dim() == 3)
  {
    mfem::Vector c(sdim);
    for (int i = 0; i < sdim; i++)
    {
      c(i) = box.center(i);
    }
    box.center = c;

    mfem::DenseMatrix a(sdim, sdim);
    for (int j = 0; j < sdim; j++)
    {
      for (int i = 0; i < sdim; i++)
      {
        a(i, j) = box.axes(i, j);
      }
    }
    box.axes = a;
  }
}

}  // namespace

BoundingBox GetBoundingBox(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                           bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  auto box = BoundingBoxFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
  TrimBoundingBox(box, mesh.SpaceDimension());
  return box;
}

BoundingBox GetBoundingBall(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  auto ball = BoundingBallFromPointCloud(mesh.GetComm(), vertices, dominant_rank);
  TrimBoundingBox(ball, mesh.SpaceDimension());
  return ball;
}

double GetProjectedLength(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                          bool bdr, const mfem::Vector &dir)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  double length;
  if (dominant_rank == Mpi::Rank(mesh.GetComm()))
  {
    Eigen::Vector3d direction = Eigen::Vector3d::Zero();
    for (int i = 0; i < dir.Size(); i++)
    {
      direction(i) = dir(i);
    }
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
                            bool bdr, const mfem::Vector &origin, bool max)
{
  std::vector<Eigen::Vector3d> vertices;
  int dominant_rank = CollectPointCloudOnRoot(mesh, marker, bdr, vertices);
  double dist;
  if (dominant_rank == Mpi::Rank(mesh.GetComm()))
  {
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
    for (int i = 0; i < origin.Size(); i++)
    {
      x0(i) = origin(i);
    }
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

// Given a mesh and boundary attribute marker array, compute a normal for the surface. If
// not averaging, use the first entry.
mfem::Vector GetSurfaceNormal(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                              bool average)
{
  int dim = mesh.SpaceDimension();
  mfem::IsoparametricTransformation T;
  mfem::Vector loc_normal(dim), normal(dim);
  normal = 0.0;
  if (mesh.Dimension() == mesh.SpaceDimension())
  {
    // Loop over boundary elements. Exit early if not averaging and non-zero normal.
    for (int i = 0; i < mesh.GetNBE() && !(!average && normal.Norml2() > 0.0); i++)
    {
      if (!marker[mesh.GetBdrAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetBdrElementTransformation(i, &T);
      mesh::Normal(T, loc_normal, &normal);
      normal += loc_normal;
    }
  }
  else
  {
    // Loop over domain elements. Exit early if not averaging and non-zero normal.
    for (int i = 0; i < mesh.GetNE() && !(!average && normal.Norml2() > 0.0); i++)
    {
      if (!marker[mesh.GetAttribute(i) - 1])
      {
        continue;
      }
      mesh.GetElementTransformation(i, &T);
      mesh::Normal(T, loc_normal, &normal);
      normal += loc_normal;
    }
  }

  // If different processors have different normal orientations, take that from the lowest
  // rank processor.
  MPI_Comm comm = mesh.GetComm();
  int rank = Mpi::Size(comm);
  mfem::Vector glob_normal(dim);
  if (normal.Norml2() > 0.0)
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
    if (normal * glob_normal < 0.0)
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

std::unique_ptr<mfem::Mesh> ExtractStandalone2DSubmesh(
    const mfem::ParMesh &parent_mesh, const mfem::Array<int> &surface_attrs,
    const std::vector<int> &pec_bdr_attrs, mfem::Vector &surface_normal,
    mfem::Vector &centroid, mfem::Vector &e1, mfem::Vector &e2)
{
  MPI_Comm comm = parent_mesh.GetComm();
  const int nprocs = Mpi::Size(comm);

  // Step 1: Extract ParSubMesh from 3D boundary.
  auto par_submesh = mfem::ParSubMesh::CreateFromBoundary(parent_mesh, surface_attrs);

  // Step 2: Compute surface normal (MPI collective).
  surface_normal = GetSurfaceNormal(par_submesh);
  Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", surface_normal(0),
             surface_normal(1), surface_normal(2));

  // Step 3: Remap domain and boundary attributes.
  RemapSubMeshAttributes(par_submesh);
  RemapSubMeshBdrAttributes(par_submesh, surface_attrs);

  // Step 4: Collect PEC internal edges via global vertex matching across all ranks.
  // PEC boundary faces in the 3D parent are distributed; each rank contributes its
  // local PEC edges, gathered to rank 0 for matching against submesh edges.
  std::vector<std::array<int, 3>> pec_internal_edges;  // {v1, v2, attr}
  {
    // Global vertex numbering (collective MPI operation).
    mfem::Array<HYPRE_BigInt> pvert_gi;
    parent_mesh.GetGlobalVertexIndices(pvert_gi);

    // Each rank: collect PEC edges as global vertex pairs.
    std::vector<int> local_pec_flat;
    {
      mfem::Array<int> edges, orientations, ev;
      for (int be = 0; be < parent_mesh.GetNBE(); be++)
      {
        int attr = parent_mesh.GetBdrAttribute(be);
        if (std::find(pec_bdr_attrs.begin(), pec_bdr_attrs.end(), attr) ==
            pec_bdr_attrs.end())
        {
          continue;
        }
        parent_mesh.GetBdrElementEdges(be, edges, orientations);
        for (int j = 0; j < edges.Size(); j++)
        {
          parent_mesh.GetEdgeVertices(edges[j], ev);
          int gv0 = static_cast<int>(pvert_gi[ev[0]]);
          int gv1 = static_cast<int>(pvert_gi[ev[1]]);
          local_pec_flat.insert(local_pec_flat.end(),
                                {std::min(gv0, gv1), std::max(gv0, gv1), attr});
        }
      }
    }

    // Each rank: collect submesh edge → global vertex pair mappings.
    std::vector<int> local_sub_flat;
    {
      const auto &parent_edge_map = par_submesh.GetParentEdgeIDMap();
      for (int i = 0; i < parent_edge_map.Size(); i++)
      {
        mfem::Array<int> pev, sev;
        parent_mesh.GetEdgeVertices(parent_edge_map[i], pev);
        par_submesh.GetEdgeVertices(i, sev);
        int gv0 = static_cast<int>(pvert_gi[pev[0]]);
        int gv1 = static_cast<int>(pvert_gi[pev[1]]);
        local_sub_flat.insert(local_sub_flat.end(),
                              {std::min(gv0, gv1), std::max(gv0, gv1), sev[0], sev[1]});
      }
    }

    // Gather both sets on rank 0.
    auto GatherFlat = [&](const std::vector<int> &local, int stride) -> std::vector<int>
    {
      int local_n = static_cast<int>(local.size()) / stride;
      std::vector<int> counts(nprocs);
      MPI_Gather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);

      std::vector<int> fcounts(nprocs), fdispls(nprocs);
      int total = 0;
      for (int i = 0; i < nprocs; i++)
      {
        fdispls[i] = total * stride;
        fcounts[i] = counts[i] * stride;
        total += counts[i];
      }
      std::vector<int> result;
      if (Mpi::Root(comm))
      {
        result.resize(total * stride);
      }
      MPI_Gatherv(local.data(), local_n * stride, MPI_INT, result.data(), fcounts.data(),
                  fdispls.data(), MPI_INT, 0, comm);
      return result;
    };

    auto all_pec = GatherFlat(local_pec_flat, 3);
    auto all_sub = GatherFlat(local_sub_flat, 4);

    // Rank 0: match PEC edges against submesh edges.
    if (Mpi::Root(comm))
    {
      std::map<std::pair<int, int>, std::pair<int, int>> gvpair_to_sub;
      for (std::size_t i = 0; i < all_sub.size() / 4; i++)
      {
        gvpair_to_sub[{all_sub[4 * i], all_sub[4 * i + 1]}] = {all_sub[4 * i + 2],
                                                               all_sub[4 * i + 3]};
      }
      std::set<std::pair<int, int>> seen;
      for (std::size_t i = 0; i < all_pec.size() / 3; i++)
      {
        auto key = std::make_pair(all_pec[3 * i], all_pec[3 * i + 1]);
        auto it = gvpair_to_sub.find(key);
        if (it != gvpair_to_sub.end())
        {
          auto [sv0, sv1] = it->second;
          if (seen.insert({std::min(sv0, sv1), std::max(sv0, sv1)}).second)
          {
            pec_internal_edges.push_back({sv0, sv1, all_pec[3 * i + 2]});
          }
        }
      }
    }

    // Broadcast PEC edges to all ranks.
    int n_pec = static_cast<int>(pec_internal_edges.size());
    MPI_Bcast(&n_pec, 1, MPI_INT, 0, comm);
    pec_internal_edges.resize(n_pec);
    MPI_Bcast(pec_internal_edges.data(), n_pec * 3, MPI_INT, 0, comm);

    Mpi::Print(" Found {:d} PEC internal edges on the cross-section\n", n_pec);
  }

  // Step 5: Gather parallel submesh as serial mesh. PrintAsOne replaces domain element
  // attributes with (rank + 1), so we gather the real attributes separately.
  std::vector<int> global_elem_attrs, global_bdr_attrs;
  {
    auto GatherAttrs = [&](int n_local, auto get_attr) -> std::vector<int>
    {
      std::vector<int> local(n_local);
      for (int i = 0; i < n_local; i++)
      {
        local[i] = get_attr(i);
      }
      std::vector<int> counts(nprocs);
      MPI_Gather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
      int total = 0;
      std::vector<int> displs(nprocs);
      for (int i = 0; i < nprocs; i++)
      {
        displs[i] = total;
        total += counts[i];
      }
      std::vector<int> result(total);
      MPI_Gatherv(local.data(), n_local, MPI_INT, result.data(), counts.data(),
                  displs.data(), MPI_INT, 0, comm);
      // Broadcast to all ranks.
      MPI_Bcast(&total, 1, MPI_INT, 0, comm);
      result.resize(total);
      MPI_Bcast(result.data(), total, MPI_INT, 0, comm);
      return result;
    };

    global_elem_attrs = GatherAttrs(par_submesh.GetNE(),
                                    [&](int i) { return par_submesh.GetAttribute(i); });
    global_bdr_attrs = GatherAttrs(par_submesh.GetNBE(),
                                   [&](int i) { return par_submesh.GetBdrAttribute(i); });
  }

  // PrintAsOne (collective) + broadcast serial mesh string.
  std::string serial_str;
  {
    std::ostringstream oss;
    par_submesh.PrintAsOne(oss);
    serial_str = oss.str();
  }
  int str_size = static_cast<int>(serial_str.size());
  MPI_Bcast(&str_size, 1, MPI_INT, 0, comm);
  serial_str.resize(str_size);
  MPI_Bcast(serial_str.data(), str_size, MPI_CHAR, 0, comm);

  // Read serial mesh and restore correct attributes.
  std::istringstream iss(serial_str);
  auto serial_mesh = std::make_unique<mfem::Mesh>(iss);
  int ne_total = static_cast<int>(global_elem_attrs.size());
  MFEM_VERIFY(serial_mesh->GetNE() == ne_total, "PrintAsOne element count mismatch!");
  for (int i = 0; i < ne_total; i++)
  {
    serial_mesh->SetAttribute(i, global_elem_attrs[i]);
  }
  int nbe_total = static_cast<int>(global_bdr_attrs.size());
  for (int i = 0; i < std::min(nbe_total, serial_mesh->GetNBE()); i++)
  {
    serial_mesh->SetBdrAttribute(i, global_bdr_attrs[i]);
  }
  serial_mesh->SetAttributes();

  // Step 6: Add PEC internal edges as boundary elements.
  if (!pec_internal_edges.empty())
  {
    int dim_s = serial_mesh->Dimension();
    int sdim_s = serial_mesh->SpaceDimension();
    int nv_s = serial_mesh->GetNV();
    int ne_s = serial_mesh->GetNE();
    int nbe_s = serial_mesh->GetNBE();
    int n_pec = static_cast<int>(pec_internal_edges.size());

    auto new_mesh = std::make_unique<mfem::Mesh>(dim_s, nv_s, ne_s, nbe_s + n_pec, sdim_s);
    for (int v = 0; v < nv_s; v++)
    {
      new_mesh->AddVertex(serial_mesh->GetVertex(v));
    }
    for (int e = 0; e < ne_s; e++)
    {
      new_mesh->AddElement(serial_mesh->GetElement(e)->Duplicate(new_mesh.get()));
    }
    for (int be = 0; be < nbe_s; be++)
    {
      new_mesh->AddBdrElement(serial_mesh->GetBdrElement(be)->Duplicate(new_mesh.get()));
    }
    for (const auto &edge : pec_internal_edges)
    {
      int vs[2] = {edge[0], edge[1]};
      new_mesh->AddBdrSegment(vs, edge[2]);
    }
    new_mesh->FinalizeTopology();
    new_mesh->EnsureNodes();
    serial_mesh = std::move(new_mesh);
  }

  // Step 7: Project from 3D to 2D coordinates.
  {
    // Build orthonormal tangent frame from the surface normal.
    e1.SetSize(3);
    e2.SetSize(3);
    {
      int min_idx = 0;
      double min_val = std::abs(surface_normal(0));
      for (int d = 1; d < 3; d++)
      {
        if (std::abs(surface_normal(d)) < min_val)
        {
          min_val = std::abs(surface_normal(d));
          min_idx = d;
        }
      }
      mfem::Vector axis(3);
      axis = 0.0;
      axis(min_idx) = 1.0;
      e1(0) = axis(1) * surface_normal(2) - axis(2) * surface_normal(1);
      e1(1) = axis(2) * surface_normal(0) - axis(0) * surface_normal(2);
      e1(2) = axis(0) * surface_normal(1) - axis(1) * surface_normal(0);
      e1 /= e1.Norml2();
      e2(0) = surface_normal(1) * e1(2) - surface_normal(2) * e1(1);
      e2(1) = surface_normal(2) * e1(0) - surface_normal(0) * e1(2);
      e2(2) = surface_normal(0) * e1(1) - surface_normal(1) * e1(0);
      e2 /= e2.Norml2();
    }

    // Compute centroid.
    centroid.SetSize(3);
    centroid = 0.0;
    int nv = serial_mesh->GetNV();
    for (int i = 0; i < nv; i++)
    {
      const double *v = serial_mesh->GetVertex(i);
      for (int d = 0; d < 3; d++)
      {
        centroid(d) += v[d];
      }
    }
    if (nv > 0)
    {
      centroid /= nv;
    }

    // Project node coordinates to 2D.
    int mesh_order = 1;
    if (serial_mesh->GetNodes())
    {
      mesh_order = serial_mesh->GetNodes()->FESpace()->GetMaxElementOrder();
    }
    mfem::GridFunction *nodes3d = serial_mesh->GetNodes();
    const int sdim = serial_mesh->SpaceDimension();
    if (nodes3d)
    {
      const int npts = nodes3d->Size() / sdim;
      std::vector<std::array<double, 2>> projected(npts);
      for (int i = 0; i < npts; i++)
      {
        double coords[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < sdim; d++)
        {
          int idx = (nodes3d->FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                        ? d * npts + i
                        : i * sdim + d;
          coords[d] = (*nodes3d)(idx);
        }
        double dx = coords[0] - centroid(0);
        double dy = coords[1] - centroid(1);
        double dz = coords[2] - centroid(2);
        projected[i][0] = dx * e1(0) + dy * e1(1) + dz * e1(2);
        projected[i][1] = dx * e2(0) + dy * e2(1) + dz * e2(2);
      }
      serial_mesh->SetCurvature(mesh_order, false, 2, mfem::Ordering::byNODES);
      mfem::GridFunction *nodes2d = serial_mesh->GetNodes();
      for (int i = 0; i < npts; i++)
      {
        (*nodes2d)(0 * npts + i) = projected[i][0];
        (*nodes2d)(1 * npts + i) = projected[i][1];
      }
    }
  }

  return serial_mesh;
}

void RemapSubMeshAttributes(mfem::ParSubMesh &submesh)
{
  MFEM_VERIFY(submesh.GetFrom() == mfem::SubMesh::From::Boundary,
              "RemapSubMeshAttributes requires a boundary ParSubMesh!");

  const auto &parent = *submesh.GetParent();
  const mfem::Array<int> &parent_elem_map = submesh.GetParentElementIDMap();

  mfem::FaceElementTransformations FET;
  mfem::IsoparametricTransformation T1, T2;
  for (int i = 0; i < submesh.GetNE(); i++)
  {
    int parent_bdr_elem = parent_elem_map[i];
    BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(parent_bdr_elem,
                                                                     parent, FET, T1, T2);
    // Use the neighboring domain element's attribute. For internal boundaries, use the
    // element with the lower attribute (same convention as Palace's mesh.cpp).
    int nbr_attr = FET.Elem1->Attribute;
    if (FET.Elem2 && FET.Elem2->Attribute < nbr_attr)
    {
      nbr_attr = FET.Elem2->Attribute;
    }
    submesh.SetAttribute(i, nbr_attr);
  }
}

void RemapSubMeshBdrAttributes(mfem::ParSubMesh &submesh,
                               const mfem::Array<int> &surface_attrs)
{
  MFEM_VERIFY(submesh.GetFrom() == mfem::SubMesh::From::Boundary,
              "RemapSubMeshBdrAttributes requires a boundary ParSubMesh!");

  const auto &parent = *submesh.GetParent();

  // Build a set of surface attributes for quick lookup.
  std::unordered_set<int> surface_attr_set;
  for (int i = 0; i < surface_attrs.Size(); i++)
  {
    surface_attr_set.insert(surface_attrs[i]);
  }

  // Build a map from parent edge index to parent boundary face attribute. For edges shared
  // by multiple parent boundary faces, prefer the face that is NOT part of the mode
  // analysis surface (since the submesh boundary lies at the intersection of the surface
  // with adjacent boundaries, and we want the adjacent boundary's attribute).
  std::unordered_map<int, int> edge_to_bdr_attr;
  mfem::Array<int> edges, orientations;
  for (int be = 0; be < parent.GetNBE(); be++)
  {
    int attr = parent.GetBdrAttribute(be);
    parent.GetBdrElementEdges(be, edges, orientations);
    bool is_surface = (surface_attr_set.count(attr) > 0);
    for (int j = 0; j < edges.Size(); j++)
    {
      auto it = edge_to_bdr_attr.find(edges[j]);
      if (it == edge_to_bdr_attr.end())
      {
        edge_to_bdr_attr[edges[j]] = attr;
      }
      else if (!is_surface)
      {
        // Prefer the non-surface attribute (adjacent boundary's attribute).
        it->second = attr;
      }
    }
  }

  // For each submesh boundary element (a 1D edge on the 2D submesh), trace back to the
  // parent edge and assign the corresponding boundary attribute.
  const mfem::Array<int> &parent_edge_map = submesh.GetParentEdgeIDMap();
  for (int sbe = 0; sbe < submesh.GetNBE(); sbe++)
  {
    int submesh_edge = submesh.GetBdrElementFaceIndex(sbe);
    MFEM_ASSERT(submesh_edge >= 0 && submesh_edge < parent_edge_map.Size(),
                "Submesh boundary edge index out of parent edge map range!");
    int parent_edge = parent_edge_map[submesh_edge];
    auto it = edge_to_bdr_attr.find(parent_edge);
    if (it != edge_to_bdr_attr.end())
    {
      submesh.SetBdrAttribute(sbe, it->second);
    }
    // If not found, the edge has no parent boundary face (unlikely for a boundary submesh
    // but possible for internal edges). Leave the default attribute.
  }

  // Rebuild the bdr_attributes array to reflect the new attribute values.
  submesh.SetAttributes();
}

mfem::Vector ProjectSubmeshTo2D(mfem::ParMesh &submesh, mfem::Vector *out_centroid,
                                mfem::Vector *out_e1, mfem::Vector *out_e2)
{
  MFEM_VERIFY(submesh.Dimension() == 2 && submesh.SpaceDimension() == 3,
              "ProjectSubmeshTo2D requires a 2D submesh with 3D ambient coordinates!");

  // Compute the surface normal from all domain elements of the submesh.
  mfem::Vector normal = GetSurfaceNormal(submesh);

  // Build an orthonormal tangent frame (e1, e2) perpendicular to the normal.
  // Choose e1 as the cross product of normal with the axis most perpendicular to it.
  mfem::Vector e1(3), e2(3);
  {
    // Pick the coordinate axis least aligned with the normal.
    int min_idx = 0;
    double min_val = std::abs(normal(0));
    for (int d = 1; d < 3; d++)
    {
      if (std::abs(normal(d)) < min_val)
      {
        min_val = std::abs(normal(d));
        min_idx = d;
      }
    }
    mfem::Vector axis(3);
    axis = 0.0;
    axis(min_idx) = 1.0;

    // e1 = normalize(axis x normal)
    e1(0) = axis(1) * normal(2) - axis(2) * normal(1);
    e1(1) = axis(2) * normal(0) - axis(0) * normal(2);
    e1(2) = axis(0) * normal(1) - axis(1) * normal(0);
    e1 /= e1.Norml2();

    // e2 = normal x e1
    e2(0) = normal(1) * e1(2) - normal(2) * e1(1);
    e2(1) = normal(2) * e1(0) - normal(0) * e1(2);
    e2(2) = normal(0) * e1(1) - normal(1) * e1(0);
    e2 /= e2.Norml2();
  }

  // Compute the centroid of all submesh vertices (using local vertices, then averaging
  // across MPI ranks).
  mfem::Vector centroid(3);
  centroid = 0.0;
  int nv_local = submesh.GetNV();
  for (int i = 0; i < nv_local; i++)
  {
    const double *v = submesh.GetVertex(i);
    for (int d = 0; d < 3; d++)
    {
      centroid(d) += v[d];
    }
  }
  int nv_global = nv_local;
  MPI_Comm comm = submesh.GetComm();
  Mpi::GlobalSum(3, centroid.HostReadWrite(), comm);
  Mpi::GlobalSum(1, &nv_global, comm);
  if (nv_global > 0)
  {
    centroid /= nv_global;
  }

  // Read the 3D node coordinates, project to 2D, then rebuild the mesh with 2D nodes.
  // We read all coordinates first, then use SetCurvature to rebuild with SpaceDim=2.
  {
    // Determine the mesh curvature order from the existing nodes. On empty partitions
    // (zero elements), GetMaxElementOrder() may fail, so use MPI reduction.
    int mesh_order = 0;
    if (submesh.GetNodes() && submesh.GetNE() > 0)
    {
      mesh_order = submesh.GetNodes()->FESpace()->GetMaxElementOrder();
    }
    Mpi::GlobalMax(1, &mesh_order, comm);
    if (mesh_order < 1)
    {
      mesh_order = 1;
    }

    // Read all 3D node coordinates. For high-order meshes, nodes includes interior DOFs;
    // for linear meshes with nodes, it's just vertices.
    mfem::GridFunction *nodes3d = submesh.GetNodes();
    const int sdim = submesh.SpaceDimension();
    std::vector<std::array<double, 2>> projected;

    if (nodes3d)
    {
      const int npts = nodes3d->Size() / sdim;
      projected.resize(npts);
      for (int i = 0; i < npts; i++)
      {
        double coords[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < sdim; d++)
        {
          // Handle both byNODES and byVDIM ordering.
          int idx = (nodes3d->FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                        ? d * npts + i
                        : i * sdim + d;
          coords[d] = (*nodes3d)(idx);
        }
        double dx = coords[0] - centroid(0);
        double dy = coords[1] - centroid(1);
        double dz = coords[2] - centroid(2);
        projected[i][0] = dx * e1(0) + dy * e1(1) + dz * e1(2);
        projected[i][1] = dx * e2(0) + dy * e2(1) + dz * e2(2);
      }
    }
    else
    {
      projected.resize(nv_local);
      for (int i = 0; i < nv_local; i++)
      {
        const double *v = submesh.GetVertex(i);
        double dx = v[0] - centroid(0);
        double dy = v[1] - centroid(1);
        double dz = v[2] - centroid(2);
        projected[i][0] = dx * e1(0) + dy * e1(1) + dz * e1(2);
        projected[i][1] = dx * e2(0) + dy * e2(1) + dz * e2(2);
      }
    }

    // Rebuild the mesh with 2D coordinates. SetCurvature creates a new Nodes
    // GridFunction with the specified vdim (SpaceDim), replacing the old 3D one.
    if (submesh.GetNE() > 0)
    {
      submesh.SetCurvature(mesh_order, false, 2, mfem::Ordering::byNODES);
      mfem::GridFunction *nodes2d = submesh.GetNodes();
      const int npts = static_cast<int>(projected.size());
      MFEM_VERIFY(nodes2d->Size() == 2 * npts,
                  "ProjectSubmeshTo2D: mismatch between projected points ("
                      << npts << ") and new nodes size (" << nodes2d->Size() << ")!");
      for (int i = 0; i < npts; i++)
      {
        (*nodes2d)(0 * npts + i) = projected[i][0];
        (*nodes2d)(1 * npts + i) = projected[i][1];
      }
    }
    else
    {
      // Empty partition: SetCurvature calls GetTypicalElementGeometry() which fails with
      // zero elements. Manually create an empty Nodes GridFunction with vdim=2 to set
      // SpaceDimension to 2.
      auto *fec = new mfem::H1_FECollection(mesh_order, submesh.Dimension());
      auto *fes = new mfem::FiniteElementSpace(&submesh, fec, 2, mfem::Ordering::byNODES);
      auto *nodes = new mfem::GridFunction(fes);
      nodes->MakeOwner(fec);
      submesh.NewNodes(*nodes, true);
    }
  }

  MFEM_VERIFY(submesh.SpaceDimension() == 2,
              "ProjectSubmeshTo2D: SpaceDimension should be 2 after projection!");

  // Output the transform parameters for projecting additional coordinates.
  if (out_centroid)
  {
    *out_centroid = centroid;
  }
  if (out_e1)
  {
    *out_e1 = e1;
  }
  if (out_e2)
  {
    *out_e2 = e2;
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
        mesh::DimensionalizeMesh(smesh, iodata.units.GetMeshLengthRelativeScale());
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
  if (!boundaries.periodic.boundary_pairs.empty())
  {
    auto periodic_mesh = std::move(mesh);

    for (const auto &data : boundaries.periodic.boundary_pairs)
    {
      auto periodic_mapping = mesh::DeterminePeriodicVertexMapping(periodic_mesh, data);
      if (!periodic_mapping.empty())
      {
        auto p_mesh = std::make_unique<mfem::Mesh>(
            mfem::Mesh::MakePeriodic(*periodic_mesh, periodic_mapping));
        periodic_mesh = std::move(p_mesh);
      }
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
               const std::vector<int> &mat_attr_list)
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

void SplitMeshElements(std::unique_ptr<mfem::Mesh> &orig_mesh, bool make_simplex,
                       bool make_hex)
{
  if (!make_simplex && !make_hex)
  {
    return;
  }
  mfem::Mesh *mesh = orig_mesh.get();
  mfem::Mesh new_mesh;

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
      int ne = mesh->GetNE();
      new_mesh = mesh::MeshTetToHex(*mesh);
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
  orig_mesh = std::make_unique<mfem::Mesh>(std::move(new_mesh));  // Call move constructor
  orig_mesh->FinalizeTopology();
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
    int attr = mesh.GetBdrAttribute(be);
    if (!boundaries.periodic.boundary_pairs.empty())
    {
      for (const auto &data : boundaries.periodic.boundary_pairs)
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
                "A non-periodic face ("
                    << f << ") cannot have multiple boundary elements! Attributes: " << attr
                    << ' ' << mesh.GetBdrAttribute(face_to_be[f]));
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
  // Container with keys being pairs of vertex indices (match row/column indices of v_to_v)
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

int AddInterfaceBdrElements(IoData &iodata, std::unique_ptr<mfem::Mesh> &orig_mesh,
                            std::unordered_map<int, int> &face_to_be, MPI_Comm comm)
{
  // Exclude some internal boundary conditions for which cracking would give invalid
  // results: lumpedports in particular.
  const auto crack_boundary_attributes = [&iodata]()
  {
    auto cba = iodata.boundaries.attributes;
    // Remove lumped port attributes.
    for (const auto &[idx, data] : iodata.boundaries.lumpedport)
    {
      for (const auto &e : data.elements)
      {
        auto attr_in_elem = [&](auto x)
        {
          return std::find(e.attributes.begin(), e.attributes.end(), x) !=
                 e.attributes.end();
        };
        cba.erase(std::remove_if(cba.begin(), cba.end(), attr_in_elem), cba.end());
      }
    }
    return cba;
  }();

  // Return if nothing to do. Otherwise, count vertices and boundary elements to add.
  if (crack_boundary_attributes.empty() && !iodata.model.add_bdr_elements)
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
  if (!crack_boundary_attributes.empty() && iodata.model.crack_bdr_elements)
  {
    auto crack_bdr_marker = mesh::AttrToMarker(
        orig_mesh->bdr_attributes.Size() ? orig_mesh->bdr_attributes.Max() : 0,
        crack_boundary_attributes, true);
    std::unordered_set<int> external_attributes;
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
          iodata.boundaries.cracked_attributes.insert(orig_mesh->GetBdrAttribute(be));
        }
        else
        {
          external_attributes.insert(orig_mesh->GetBdrAttribute(be));
        }
      }
    }
    MFEM_VERIFY(crack_bdr_elem.empty() || !orig_mesh->Nonconforming(),
                "Duplicating internal boundary elements for interior boundaries is not "
                "supported for nonconforming meshes!");
    std::vector<int> mixed_attributes;
    std::set_intersection(iodata.boundaries.cracked_attributes.begin(),
                          iodata.boundaries.cracked_attributes.end(),
                          external_attributes.begin(), external_attributes.end(),
                          std::back_inserter(mixed_attributes));
    if (!mixed_attributes.empty())
    {
      MFEM_WARNING("Found boundary attribute with internal and external boundary elements: "
                   << fmt::format("{}", fmt::join(mixed_attributes, " "))
                   << ". Impedance boundary conditions for these attributes will give "
                      "erroneous results, consider separating into different attributes!");
    }
    vert_to_elem.reset(orig_mesh->GetVertexToElementTable());  // Owned by caller
    const mfem::Table &elem_to_face =
        (orig_mesh->Dimension() == 2 ? orig_mesh->ElementToEdgeTable()
                                     : orig_mesh->ElementToFaceTable());
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
          // elements at the cost of a potential minor mesh quality degradation.
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

  // Export mesh after pre-processing, before cracking boundary elements.
  if (iodata.model.export_prerefined_mesh && Mpi::Root(comm))
  {
    auto pos = iodata.model.mesh.find_last_of(".");
    std::string meshfile = iodata.model.mesh.substr(0, pos) + "_preprocessed.mesh";
    std::ofstream fo(meshfile);
    fo.precision(MSH_FLT_PRECISION);
    orig_mesh->Print(fo);
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
  if (!crack_boundary_attributes.empty() && !crack_bdr_elem.empty())
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
    const mfem::Table &elem_to_face =
        (orig_mesh->Dimension() == 2 ? orig_mesh->ElementToEdgeTable()
                                     : orig_mesh->ElementToFaceTable());
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
        bdr_el->GetVertices()[j] =
            el->GetVertices()[(orig_mesh->Dimension() == 2 ? el->GetEdgeVertices(i)
                                                           : el->GetFaceVertices(i))[j]];
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
          bdr_el->GetVertices()[j] =
              el->GetVertices()[(orig_mesh->Dimension() == 2 ? el->GetEdgeVertices(i)
                                                             : el->GetFaceVertices(i))[j]];
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
    const mfem::Table &elem_to_face =
        (orig_mesh->Dimension() == 2 ? orig_mesh->ElementToEdgeTable()
                                     : orig_mesh->ElementToFaceTable());
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
          bdr_el->GetVertices()[j] =
              el->GetVertices()[(orig_mesh->Dimension() == 2 ? el->GetEdgeVertices(i)
                                                             : el->GetFaceVertices(i))[j]];
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
  if (!crack_boundary_attributes.empty() && !crack_bdr_elem.empty() &&
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
    const mfem::Table &elem_to_face =
        (orig_mesh->Dimension() == 2 ? orig_mesh->ElementToEdgeTable()
                                     : orig_mesh->ElementToFaceTable());
    const mfem::Table &new_elem_to_face =
        (new_mesh->Dimension() == 2 ? new_mesh->ElementToEdgeTable()
                                    : new_mesh->ElementToFaceTable());
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
            NodeUpdate(el->GetVertices()[(orig_mesh->Dimension() == 2
                                              ? el->GetEdgeVertices(i)
                                              : el->GetFaceVertices(i))[j]]);
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
