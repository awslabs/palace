// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "geodata.hpp"
#include "geodata_impl.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory_resource>
#include <mutex>
#include <numbers>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#endif
#include <Eigen/Dense>
#include <mpi.h>
#if defined(MFEM_USE_ZLIB)
#include <zlib.h>
#endif
#include <ptscotch.h>
#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "fem/interpolator.hpp"
#include "linalg/sliverpatch.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"
#include "utils/filesystem.hpp"
#include "utils/gmshreader.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/meshio.hpp"
#include "utils/ncmeshdist.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"
#include "utils/perfstats.hpp"
#include "utils/prettyprint.hpp"
#include "utils/timer.hpp"
#include "utils/units.hpp"

namespace palace
{

using Vector3dMap = Eigen::Map<Eigen::Vector3d>;
using CVector3dMap = Eigen::Map<const Eigen::Vector3d>;

namespace
{

// Floating point precision for mesh IO. This precision is important, make sure nothing is
// lost!
constexpr auto MSH_FLT_PRECISION = std::numeric_limits<double>::max_digits10;

void UnpinThread(std::thread &t);

// Wait for all processes without spinning, leaving the CPUs to the threads of the processes
// working meanwhile.
void SleepyBarrier(MPI_Comm comm)
{
  MPI_Request req;
  MPI_Ibarrier(comm, &req);
  for (int flag = 0;;)
  {
    MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
    if (flag)
    {
      break;
    }
    std::this_thread::sleep_for(std::chrono::microseconds(200));
  }
}

// Threads for the serial mesh preprocessing of the processes holding the serial mesh: one
// per process of the node, the others waiting without spinning (set by mesh::Load).
int serial_mesh_threads = 1;

// Call f(k, t) for k in [0, n) from nt threads t (the calling one is t = 0).
template <typename Func>
void ParallelForThreads(std::size_t n, int nt, Func &&f)
{
  std::atomic<std::size_t> next = 0;
  auto Work = [&](int t)
  {
    for (std::size_t k; (k = next.fetch_add(1)) < n;)
    {
      f(k, t);
    }
  };
  std::vector<std::thread> threads;
  for (int t = 1; t < nt; t++)
  {
    threads.emplace_back(Work, t);
    UnpinThread(threads.back());
  }
  Work(0);
  for (auto &t : threads)
  {
    t.join();
  }
}

// mfem::Mesh::GetEdgeLengths2 for the vertices of a mesh without nodes (the squared length
// of each edge of v_to_v, by edge index), with nt threads over the rows of v_to_v.
void EdgeLengths2(const mfem::DSTable &v_to_v, const mfem::Array<mfem::Vertex> &vertices,
                  int sdim, int nedges, int nt, mfem::Array<mfem::real_t> &lengths)
{
  lengths.SetSize(nedges);
  constexpr int block = 4096;
  const int nrows = v_to_v.NumberOfRows();
  ParallelForThreads((nrows + block - 1) / block, nt,
                     [&](std::size_t b, int)
                     {
                       const int end = std::min(nrows, static_cast<int>(b + 1) * block);
                       for (int i = static_cast<int>(b) * block; i < end; i++)
                       {
                         for (mfem::DSTable::RowIterator it(v_to_v, i); !it; ++it)
                         {
                           const mfem::real_t *vi = vertices[i]();
                           const mfem::real_t *vj = vertices[it.Column()]();
                           mfem::real_t l = 0.;
                           for (int k = 0; k < sdim; k++)
                           {
                             l += (vi[k] - vj[k]) * (vi[k] - vj[k]);
                           }
                           lengths[it.Index()] = l;
                         }
                       }
                     });
}

// Load the serial mesh from disk.
std::unique_ptr<mfem::Mesh> LoadMesh(const std::string &, bool,
                                     const config::BoundaryData &, bool = false,
                                     bool * = nullptr);

// Read a serial mesh from a stream like the mfem::Mesh stream constructor, but removing the
// boundary elements which are not faces of any element of the mesh.
std::unique_ptr<mfem::Mesh> ReadMeshStream(std::istream &, bool, bool, bool,
                                           bool * = nullptr);
std::unique_ptr<mfem::Mesh> FinishReadMesh(std::unique_ptr<mfem::Mesh>, bool, bool,
                                           bool * = nullptr);
bool HasTopology(const mfem::Mesh &);

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
                                              const int *, const std::string & = "",
                                              bool = false);

// Rebalance a conformal mesh across processor ranks, using the MeshPartitioner. Gathers the
// mesh onto the root rank before scattering the partitioned mesh.
void RebalanceConformalMesh(std::unique_ptr<mfem::ParMesh> &);

// Apply box/sphere region-based refinement to a serial mesh in-place. For tensor-element
// meshes the refinement is non-conforming and the mesh must already be nonconforming.
// Refinement levels are applied one at a time; at each level every element whose geometry
// intersects at least one active (level-appropriate) refinement region is refined.
void RegionRefine(const config::RefinementData &refinement, mfem::Mesh &mesh);

}  // namespace

namespace mesh
{

namespace
{

// True if the input file requests any h-refinement — adaptive, box, or sphere.
bool UseAmr(const config::RefinementData &refinement)
{
  if (refinement.max_it > 0)
  {
    return true;
  }
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
}

}  // namespace

namespace
{

// A mesh marking its tetrahedra for refinement as mfem::Mesh::MarkTetMeshForRefinement does
// (the longest edge of each tetrahedron, ties broken by the edge index), with threads for
// the elements, whose marking changes only themselves (see serial_mesh_threads).
class ThreadedMarkMesh : public mfem::Mesh
{
  void MarkTetMeshForRefinement(const mfem::DSTable &v_to_v) override
  {
    mfem::Array<mfem::real_t> lengths;
    if (Nodes)
    {
      GetEdgeLengths2(v_to_v, lengths);
    }
    else
    {
      EdgeLengths2(v_to_v, vertices, spaceDim, NumOfEdges, serial_mesh_threads, lengths);
    }
    mfem::Array<int> idx(NumOfEdges);
    std::iota(idx.begin(), idx.end(), 0);
    constexpr int block = 4096;
    const int nt = serial_mesh_threads;
    ParallelForThreads(
        (NumOfElements + block - 1) / block, nt,
        [&](std::size_t b, int)
        {
          const int end = std::min(NumOfElements, static_cast<int>(b + 1) * block);
          for (int i = static_cast<int>(b) * block; i < end; i++)
          {
            if (elements[i]->GetType() == mfem::Element::TETRAHEDRON)
            {
              MFEM_ASSERT(dynamic_cast<mfem::Tetrahedron *>(elements[i]),
                          "Unexpected non-Tetrahedron element type!");
              static_cast<mfem::Tetrahedron *>(elements[i])->MarkEdge(v_to_v, lengths, idx);
            }
          }
        });
    ParallelForThreads(
        (NumOfBdrElements + block - 1) / block, nt,
        [&](std::size_t b, int)
        {
          const int end = std::min(NumOfBdrElements, static_cast<int>(b + 1) * block);
          for (int i = static_cast<int>(b) * block; i < end; i++)
          {
            if (boundary[i]->GetType() == mfem::Element::TRIANGLE)
            {
              MFEM_ASSERT(dynamic_cast<mfem::Triangle *>(boundary[i]),
                          "Unexpected non-Triangle element type!");
              static_cast<mfem::Triangle *>(boundary[i])->MarkEdge(v_to_v, lengths, idx);
            }
          }
        });
  }

public:
  explicit ThreadedMarkMesh(mfem::Mesh &&mesh) : mfem::Mesh(std::move(mesh)) {}
};

}  // namespace

std::unique_ptr<mfem::Mesh> Load(IoData &iodata, MPI_Comm comm)
{
  BlockTimer bt0(Timer::MESH_PREPROCESS);

  // The processes holding the serial mesh preprocess it with one thread per process of
  // their node (PALACE_SETUP_THREADS overrides), while the others wait without spinning
  // (until the end of Load).
  {
    MPI_Comm node_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
    serial_mesh_threads = Mpi::Size(node_comm);
    MPI_Comm_free(&node_comm);
    if (const char *env = std::getenv("PALACE_SETUP_THREADS"))
    {
      serial_mesh_threads = std::max(1, std::atoi(env));
    }
  }

  const auto &refinement = iodata.model.refinement;
  const bool use_amr = UseAmr(refinement);

  // Unless an element splitting, reordering, or serial refinement (which depend on the
  // element orientations) or periodicity requires the mesh finalized once read, the
  // orientation fix after reading is done with the refinement marking below.
  const bool defer_finalize = !iodata.model.make_simplex && !iodata.model.make_hex &&
                              !iodata.model.reorder_elements &&
                              refinement.ser_uniform_ref_levels == 0 &&
                              iodata.boundaries.periodic.boundary_pairs.empty();
  bool finalize_deferred = false;

  // Root loads first so conformality of the input mesh can be factored into the
  // distribution-path decision. Non-conformal AMR on an already-nonconforming mesh
  // forces the byte-string broadcast path, since MFEM's MeshPartitioner doesn't accept
  // nonconforming meshes. When the configuration alone excludes the partitioner, one
  // process per node reads the serial mesh, all of them at once.
  std::unique_ptr<mfem::Mesh> smesh;
  bool use_mesh_partitioner = !use_amr || !refinement.nonconformal;
  const bool load_per_node = !use_mesh_partitioner;
  if (!load_per_node)
  {
    BlockTimer bt(Timer::IO);
    if (Mpi::Root(comm))
    {
      perf::ScopedTally t(perf::RegisterTally("Startup/MeshLoad/Root", Timer::IO));
      smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature, iodata.boundaries,
                       true, defer_finalize ? &finalize_deferred : nullptr);
      use_mesh_partitioner &= smesh->Conforming();
    }
    perf::ScopedTally t(perf::RegisterTally("Wait/IO/MeshLoad", Timer::IO));
    Mpi::Broadcast(1, &use_mesh_partitioner, 0, comm);
  }

  MPI_Comm node_comm = MPI_COMM_NULL;
  if (!use_mesh_partitioner)
  {
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
  }
  {
    BlockTimer bt1(Timer::IO);
    if (!use_mesh_partitioner && Mpi::Root(node_comm) &&
        (load_per_node || !Mpi::Root(comm)))
    {
      // Only one process per node reads the serial mesh when not using the partitioner.
      perf::ScopedTally t(perf::RegisterTally("Startup/MeshLoad/NodeRoot", Timer::IO));
      smesh = LoadMesh(iodata.model.mesh, iodata.model.remove_curvature, iodata.boundaries,
                       true, defer_finalize ? &finalize_deferred : nullptr);
      MFEM_VERIFY(!(smesh->Nonconforming() && use_mesh_partitioner),
                  "Cannot use mesh partitioner on a nonconforming mesh!");
    }
    perf::ScopedTally t(perf::RegisterTally("Wait/IO/MeshLoad", Timer::IO));
    Mpi::Barrier(comm);
  }
  if (node_comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&node_comm);
  }

  if (!smesh)
  {
    SleepyBarrier(comm);
    return smesh;
  }

  // The mesh is read without nodes when it has none (a linear mesh): first-order nodes,
  // which every later step updates along with the vertices, are added only once the mesh
  // topology is final, or before a step which handles the nodes of a mesh differently than
  // its vertices. The vertices, and the nodes once added, are the same as with nodes from
  // the start.
  auto EnsureNodes = [&]()
  {
    if (!iodata.model.remove_curvature && !smesh->GetNodes())
    {
      smesh->EnsureNodes();
    }
  };

  // AMR / element compatibility.
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
  MFEM_VERIFY(smesh->Conforming() || !use_amr || refinement.nonconformal,
              "The provided mesh is nonconformal, only nonconformal AMR can be performed!");
  MFEM_VERIFY(iodata.boundaries.periodic.boundary_pairs.empty() || !use_amr ||
                  refinement.nonconformal,
              "Conformal mesh refinement is not supported on periodic meshes; set "
              "\"Refinement\"/\"Nonconformal\": true to refine a periodic mesh!");

  // Clean up unused domain elements from the mesh.
  if (iodata.model.clean_unused_elements)
  {
    std::vector<int> attr_list;
    std::merge(iodata.domains.attributes.begin(), iodata.domains.attributes.end(),
               iodata.domains.postpro.attributes.begin(),
               iodata.domains.postpro.attributes.end(), std::back_inserter(attr_list));
    attr_list.erase(std::ranges::unique(attr_list).begin(), attr_list.end());
    CleanMesh(smesh, attr_list);
  }

  // Optionally convert mesh elements to simplices or hexes.
  if (iodata.model.make_simplex || iodata.model.make_hex)
  {
    EnsureNodes();
    SplitMeshElements(smesh, iodata.model.make_simplex, iodata.model.make_hex);
  }

  // Optionally reorder elements/vertices by spatial location for cache locality.
  if (iodata.model.reorder_elements)
  {
    EnsureNodes();
    ReorderMeshElements(*smesh);
  }

  // Serial uniform refinement (seldom used; prefer parallel uniform refinement in
  // RefineMesh).
  {
    int ne = smesh->GetNE();
    if (refinement.ser_uniform_ref_levels > 0)
    {
      EnsureNodes();
    }
    for (int l = 0; l < refinement.ser_uniform_ref_levels; l++)
    {
      smesh->UniformRefinement();
    }
    if (refinement.ser_uniform_ref_levels > 0)
    {
      Mpi::Print("Serial uniform mesh refinement levels added {:d} elements (initial = "
                 "{:d}, final = {:d})\n",
                 smesh->GetNE() - ne, ne, smesh->GetNE());
    }
  }

  // Finalize before any further refinement: GeneralRefinement (used by RegionRefine)
  // needs tet refinement flags set, which Finalize(refine=true) provides. Initial
  // orientation was already fixed at load.
  {
    constexpr bool refine = true;
    const bool fix_orientation = finalize_deferred;
    ThreadedMarkMesh mesh(std::move(*smesh));
    mesh.Finalize(refine, fix_orientation);
    *smesh = std::move(mesh);
  }

  // Box / sphere region refinement on the serial mesh, before partitioning. Done here
  // (rather than in parallel RefineMesh) so the user-facing 3D box / sphere geometry
  // stays in sync with the mesh the problem actually solves on — BoundaryMode's
  // Preprocess may extract a 2D submesh from this refined 3D mesh. A nonconformal mesh
  // cannot be cracked along interior boundaries: a nonconformal input mesh, for example
  // from a previous adaptive simulation, must already be cracked along the interior
  // boundaries which are to be cracked, and region refinement makes a mesh with
  // hexahedra, prisms, or pyramids nonconformal. Both are an error when interior boundaries
  // are to be cracked, since the cracking is part of the model.
  const bool conforming = smesh->Conforming();
  std::set<int> interior_crack_attrs;
  if (iodata.model.crack_bdr_elements.value_or(true))
  {
    const auto crack_attrs =
        iodata.boundaries.GetMeshCrackAttributes(iodata.model.crack_bdr_elements);
    const int bdr_attr_max = smesh->bdr_attributes.Size() ? smesh->bdr_attributes.Max() : 0;
    const auto crack_marker = mesh::AttrToMarker(bdr_attr_max, crack_attrs, true);
    for (int be = 0; be < smesh->GetNBE(); be++)
    {
      const int attr = smesh->GetBdrAttribute(be);
      int e1, e2;
      smesh->GetFaceElements(smesh->GetBdrElementFaceIndex(be), &e1, &e2);
      if (attr > 0 && attr <= bdr_attr_max && crack_marker[attr - 1] && e1 >= 0 && e2 >= 0)
      {
        interior_crack_attrs.insert(attr);
      }
    }
  }
  MFEM_VERIFY(interior_crack_attrs.empty() || conforming,
              "The nonconformal mesh "
                  << iodata.model.mesh
                  << " is not cracked along the interior boundaries with "
                     "attributes "
                  << fmt::format("{}", fmt::join(interior_crack_attrs, ", "))
                  << " (with \"Crack\": true, or with boundary conditions "
                     "applying to each side separately), and cannot be "
                     "cracked! Use a mesh which is cracked along them, for "
                     "example from a previous adaptive simulation with the "
                     "same boundary conditions.");
  if (!refinement.GetBoxes().empty() || !refinement.GetSpheres().empty())
  {
    EnsureNodes();
  }
  RegionRefine(refinement, *smesh);
  MFEM_VERIFY(
      interior_crack_attrs.empty() || smesh->Conforming(),
      "Refinement boxes and spheres make a mesh with hexahedra, prisms, or pyramids "
      "nonconformal, which cannot be cracked along the interior boundaries with "
      "attributes "
          << fmt::format("{}", fmt::join(interior_crack_attrs, ", "))
          << " (with \"Crack\": true, or with boundary conditions applying to "
             "each side separately)! Remove the refinement regions, or use a "
             "simplex mesh.");

  // Exterior-boundary check and optional material-interface / crack boundary element
  // insertion. Only meaningful on an initial conformal mesh.
  if (smesh->Conforming())
  {
    auto face_to_be = CheckMesh(*smesh, iodata.boundaries);
    if (iodata.model.crack_bdr_elements.value_or(true) || iodata.model.add_bdr_elements)
    {
      while (AddInterfaceBdrElements(iodata, smesh, face_to_be, comm) != 1)
      {
        // May require multiple calls due to early exit/retry approach.
      }
    }
  }
  else if (conforming)
  {
    Mpi::Warning("Refinement boxes and spheres make the mesh nonconformal, skipping mesh "
                 "modification preprocessing steps!\n\n");
  }
  else
  {
    Mpi::Warning("{} is a nonconformal mesh, assumed from previous AMR iteration.\n"
                 "Skipping mesh modification preprocessing steps!\n\n",
                 iodata.model.mesh);
  }

  EnsureNodes();
  SleepyBarrier(comm);
  return smesh;
}

std::unique_ptr<mfem::ParMesh> Partition(IoData &iodata, std::unique_ptr<mfem::Mesh> smesh,
                                         MPI_Comm comm)
{
  BlockTimer bt0(Timer::MESH_PREPROCESS);
  perf::ScopedTally t(perf::RegisterTally("Startup/Partition", Timer::MESH_PREPROCESS));
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = UseAmr(refinement);

  // Re-derive the distribution-path decision that Load already used. Both halves of the
  // pipeline run on the same iodata + the same initial mesh, so the result is identical.
  bool use_mesh_partitioner = !use_amr || !refinement.nonconformal;
  if (Mpi::Root(comm) && smesh)
  {
    use_mesh_partitioner &= smesh->Conforming();
  }
  Mpi::Broadcast(1, &use_mesh_partitioner, 0, comm);

  // Whether the serial mesh is conforming before the conversion to NC below (an NC input
  // mesh, for example from a previous adaptive simulation, carries a refinement tree).
  bool input_conforming = false;
  if (Mpi::Root(comm) && smesh)
  {
    input_conforming = smesh->Conforming();
  }
  Mpi::Broadcast(1, &input_conforming, 0, comm);

  // For nonconformal AMR, each loading rank converts its own serial copy to NC here
  // (after Load has done all preprocessing and EnsureNodes), before partitioning. The
  // resulting ParMesh inherits the NC structure with Nodes intact.
  if (use_amr && refinement.nonconformal && smesh)
  {
    smesh->EnsureNCMesh(true);
  }

  // Broadcast cracked boundary attributes from root to all ranks.
  if (iodata.model.crack_bdr_elements.value_or(true) || iodata.model.add_bdr_elements)
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
  if (iodata.model.crack_bdr_elements.value_or(false))
  {
    iodata.boundaries.split_impedance_attributes = iodata.boundaries.cracked_attributes;
  }

  // Generate partitioning from the serial mesh on loading ranks.
  std::unique_ptr<int[]> partitioning;
  if (smesh)
  {
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), iodata.model.partitioning);
  }

  std::unique_ptr<mfem::ParMesh> pmesh;
  if (use_mesh_partitioner)
  {
    // A conforming mesh which is refined after its distribution (conformal adaptive or
    // parallel uniform refinement) keeps the element order of the serial mesh: the tie
    // breaks of the marking of the tetrahedra for refinement follow the local edge order.
    pmesh = DistributeMesh(comm, smesh, partitioning.get(), iodata.problem.output,
                           iodata.model.local_reorder && !use_amr &&
                               refinement.uniform_ref_levels == 0);
  }
  else
  {
    auto DistributeSerialMesh = [&]()
    {
      // Send the preprocessed serial mesh and partitioning as a byte string. The
      // serialized mesh can exceed INT_MAX bytes for large meshes, so use a 64-bit
      // length and a chunked broadcast.
      std::unique_ptr<mfem::ParMesh> spmesh;
      MPI_Comm node_comm;
      MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                          &node_comm);
      constexpr bool generate_edges = false, refine = true, fix_orientation = false;
      std::string so;
      std::int64_t slen = 0;
      if (smesh)
      {
        std::ostringstream fo(std::stringstream::out);
        fo << std::scientific;
        fo.precision(MSH_FLT_PRECISION);
        smesh->Print(fo);
        smesh.reset();  // Root process needs to rebuild the mesh to ensure consistency
                        // with the saved serial mesh (refinement marking, for example)
        so = fo.str();
        slen = static_cast<std::int64_t>(so.size());
      }
      Mpi::Broadcast(1, &slen, 0, node_comm);
      if (so.empty())
      {
        so.resize(slen);
      }
      Mpi::BroadcastLarge(slen, so.data(), 0, node_comm);
      {
        std::istringstream fi(so);
        smesh = std::make_unique<mfem::Mesh>(fi, generate_edges, refine, fix_orientation);
        so.clear();
      }
      // !use_mesh_partitioner implies use_amr && nonconformal. EnsureNCMesh has already
      // been called on each loading rank's serial copy above; the rebuilt mesh below
      // needs the same treatment.
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
      spmesh = std::make_unique<mfem::ParMesh>(comm, *smesh, partitioning.get());
      smesh.reset();
      SetRootElementOrigin(*spmesh, {});  // Root elements in the input order
      return spmesh;
    };

    // Every process builds only the NCMesh of the root elements (which mfem::ParNCMesh
    // replicates) instead of the full serial mesh; a nonconforming input mesh has its
    // refinement tree replayed in parallel. A refinement tree which cannot be replayed is
    // sent as the serialized serial mesh.
    pmesh = input_conforming ? DistributeNCMesh(comm, smesh, partitioning.get(),
                                                iodata.model.local_reorder)
                             : DistributeNCMeshTree(comm, smesh, partitioning.get());
    if (!pmesh)
    {
      pmesh = DistributeSerialMesh();
    }
    smesh.reset();
  }
  memory_reporting::TrimHeap();

  // Debug: optionally write the partitioned and final meshes to tmp/ for inspection.
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
      fo.precision(MSH_FLT_PRECISION);
      gpmesh.ParPrint(fo);
    }
    {
      std::string pfile =
          mfem::MakeParFilename(tmp.string() + "final.", Mpi::Rank(comm), ".mesh", width);
      std::ofstream fo(pfile);
      fo.precision(MSH_FLT_PRECISION);
      pmesh->ParPrint(fo);
    }
  }

  return pmesh;
}

std::unique_ptr<mfem::ParMesh> ReadMesh(IoData &iodata, MPI_Comm comm)
{
  return Partition(iodata, Load(iodata, comm), comm);
}

namespace
{

// mfem::Mesh::GetBoundingBox of a serial mesh. For first-order nodes on a simplex mesh the
// box is that of the nodal values: the points of mfem::Mesh::GetBoundingBox are the element
// vertices and edge midpoints, whose coordinates are the node values or lie between two of
// them (also in floating point), and the box is found without element transformations.
void GetSerialBoundingBox(const mfem::Mesh &mesh, mfem::Vector &bbmin, mfem::Vector &bbmax)
{
  const mfem::GridFunction *nodes = mesh.GetNodes();
  const bool simplex_mesh = [&mesh]()
  {
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const auto geom = mesh.GetElementBaseGeometry(e);
      if (geom != mfem::Geometry::TETRAHEDRON && geom != mfem::Geometry::TRIANGLE)
      {
        return false;
      }
    }
    return true;
  }();
  const auto *fespace = nodes ? nodes->FESpace() : nullptr;
  if (!nodes || !simplex_mesh || fespace->GetMaxElementOrder() != 1 ||
      !dynamic_cast<const mfem::H1_FECollection *>(fespace->FEColl()))
  {
    const_cast<mfem::Mesh &>(mesh).GetBoundingBox(bbmin, bbmax);
    return;
  }
  const int sdim = mesh.SpaceDimension();
  bbmin.SetSize(sdim);
  bbmax.SetSize(sdim);
  bbmin = mfem::infinity();
  bbmax = -mfem::infinity();
  for (int i = 0; i < fespace->GetNDofs(); i++)
  {
    for (int d = 0; d < sdim; d++)
    {
      const double x = (*nodes)(fespace->DofToVDof(i, d));
      bbmin(d) = std::min(bbmin(d), x);
      bbmax(d) = std::max(bbmax(d), x);
    }
  }
}

}  // namespace

double ComputeReferenceLength(const std::unique_ptr<mfem::Mesh> &mesh, MPI_Comm comm)
{
  double Lc_local = 0.0;
  if (mesh)
  {
    mfem::Vector bbmin, bbmax;
    GetSerialBoundingBox(*mesh, bbmin, bbmax);
    bbmax -= bbmin;
    Lc_local = *std::ranges::max_element(bbmax);
  }
  Mpi::GlobalMax(1, &Lc_local, comm);
  return Lc_local;
}

void RefineMesh(const IoData &iodata, std::vector<std::unique_ptr<mfem::ParMesh>> &mesh)
{
  // Parallel uniform refinement only. Box / sphere region refinement happens in Load.
  MFEM_VERIFY(mesh.size() == 1,
              "Input mesh vector before refinement has more than a single mesh!");
  int uniform_ref_levels = iodata.model.refinement.uniform_ref_levels;
  if (iodata.solver.linear.mg_use_mesh && iodata.solver.linear.mg_max_levels > 1)
  {
    mesh.reserve(1 + uniform_ref_levels);
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
      if (const std::vector<int> *origin = GetRootElementOrigin(*mesh[mesh.size() - 2]))
      {
        SetRootElementOrigin(*mesh.back(), std::vector<int>(*origin));
      }
    }
    mesh.back()->UniformRefinement();
  }

  // Simplex meshes need to be re-finalized in order to use local refinement (see
  // the docstring for mfem::Mesh::UniformRefinement).
  const auto element_types = mesh::CheckElements(*mesh.back());
  if (element_types.has_simplices && uniform_ref_levels > 0 &&
      iodata.model.refinement.max_it > 0)
  {
    constexpr bool refine = true, fix_orientation = false;
    Mpi::Print("\nFlattening mesh sequence:\n Local mesh refinement will start from the "
               "final uniformly-refined mesh\n");
    mesh.erase(mesh.begin(), mesh.end() - 1);
    mesh.back()->Finalize(refine, fix_orientation);
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
  PrintMeshInfo(*mesh[0], iodata, /*full=*/true, (mesh.size() > 1) ? "Coarse " : "");
  if (mesh.size() > 1)
  {
    PrintMeshInfo(*mesh.back(), iodata, /*full=*/false, "Refined ");
  }
}

void PrintMeshInfo(mfem::ParMesh &mesh, const IoData &iodata, bool full,
                   const std::string &prefix)
{
  if (full)
  {
    mfem::Vector bbmin, bbmax;
    GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
    const double Lc = iodata.units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
    Mpi::Print(mesh.GetComm(), "\nMesh curvature order: {}\nMesh bounding box:\n",
               mesh.GetNodes()
                   ? std::to_string(mesh.GetNodes()->FESpace()->GetMaxElementOrder())
                   : "None");
    if (mesh.SpaceDimension() == 3)
    {
      Mpi::Print(mesh.GetComm(),
                 " (Xmin, Ymin, Zmin) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n"
                 " (Xmax, Ymax, Zmax) = ({:+.3e}, {:+.3e}, {:+.3e}) m\n",
                 bbmin[0] * Lc, bbmin[1] * Lc, bbmin[2] * Lc, bbmax[0] * Lc, bbmax[1] * Lc,
                 bbmax[2] * Lc);
    }
    else
    {
      Mpi::Print(mesh.GetComm(),
                 " (Xmin, Ymin) = ({:+.3e}, {:+.3e}) m\n"
                 " (Xmax, Ymax) = ({:+.3e}, {:+.3e}) m\n",
                 bbmin[0] * Lc, bbmin[1] * Lc, bbmax[0] * Lc, bbmax[1] * Lc);
    }
  }
  Mpi::Print(mesh.GetComm(), "\n{}", prefix);
  mesh.PrintInfo();

  // Degenerate (sliver) tetrahedra can slow down the convergence of the iterative linear
  // solvers, see config["Solver"]["Linear"]["MGSliverAngle"].
  SliverStatistics stats;
  FindSliverElements(mesh, 0.0, &stats);
  if (stats.num_below_5 > 0)
  {
    Mpi::Print(mesh.GetComm(),
               " Tetrahedra with a minimum dihedral angle below 5 degrees: {:d} ({:d} "
               "below 2, {:d} below 1, smallest {:.3g})\n",
               stats.num_below_5, stats.num_below_2, stats.num_below_1, stats.min_angle);
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

void Nondimensionalize(const Units &units, mfem::Mesh &mesh)
{
  NondimensionalizeMesh(mesh, units.GetMeshLengthRelativeScale());
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

std::vector<int> HilbertOrder(const double *x, int n, int sdim)
{
  std::array<double, 3> lo, hi;
  lo.fill(mfem::infinity());
  hi.fill(-mfem::infinity());
  for (int i = 0; i < n; i++)
  {
    for (int d = 0; d < sdim; d++)
    {
      lo[d] = std::min(lo[d], x[i * sdim + d]);
      hi[d] = std::max(hi[d], x[i * sdim + d]);
    }
  }
  double len = 0.0;
  for (int d = 0; d < sdim; d++)
  {
    len = std::max(len, hi[d] - lo[d]);
  }
  return HilbertOrder(x, n, sdim, lo.data(), len);
}

std::vector<int> HilbertOrder(const double *x, int n, int sdim, const double *lo,
                              double len)
{
  const int b = (sdim == 3) ? 21 : 31;
  const double nb = std::ldexp(1.0, b);
  std::vector<std::pair<std::uint64_t, int>> key(n);
  for (int i = 0; i < n; i++)
  {
    std::array<std::uint32_t, 3> X = {0, 0, 0};
    for (int d = 0; d < sdim; d++)
    {
      const double t = (len > 0.0) ? (x[i * sdim + d] - lo[d]) / len : 0.0;
      X[d] = static_cast<std::uint32_t>(std::clamp(std::floor(t * nb), 0.0, nb - 1.0));
    }
    key[i] = {HilbertIndex(X, sdim, b), i};
  }
  std::sort(key.begin(), key.end());
  std::vector<int> order(n);
  for (int i = 0; i < n; i++)
  {
    order[i] = key[i].second;
  }
  return order;
}

namespace
{

// Root element origins of the reordered nonconforming meshes, by their mfem::ParNCMesh
// (which an mfem::ParMesh keeps through refinement and rebalancing), with the number of
// root elements as a check.
std::map<const mfem::NCMesh *, std::vector<int>> &RootElementOrigins()
{
  static std::map<const mfem::NCMesh *, std::vector<int>> origins;
  return origins;
}

// Element origins of the reordered conforming meshes, by their mfem::ParMesh, with the
// number of elements as a check.
std::map<const mfem::Mesh *, std::vector<std::int64_t>> &ElementOrigins()
{
  static std::map<const mfem::Mesh *, std::vector<std::int64_t>> origins;
  return origins;
}

// Access to the leaf elements of an mfem::NCMesh (the tree index of each mesh element).
struct NCMeshLeaves : public mfem::ParNCMesh
{
  static const mfem::Array<int> &Get(const mfem::NCMesh &ncmesh)
  {
    return ncmesh.*(&NCMeshLeaves::leaf_elements);
  }
};

}  // namespace

void SetRootElementOrigin(const mfem::ParMesh &mesh, std::vector<int> &&origin)
{
  if (!mesh.pncmesh)
  {
    return;
  }
  if (origin.empty())
  {
    RootElementOrigins().erase(mesh.pncmesh);
  }
  else
  {
    RootElementOrigins()[mesh.pncmesh] = std::move(origin);
  }
}

const std::vector<int> *GetRootElementOrigin(const mfem::ParMesh &mesh)
{
  if (!mesh.pncmesh)
  {
    return nullptr;
  }
  const auto it = RootElementOrigins().find(mesh.pncmesh);
  if (it == RootElementOrigins().end() ||
      static_cast<int>(it->second.size()) !=
          const_cast<mfem::ParNCMesh *>(mesh.pncmesh)->GetNumRootElements())
  {
    return nullptr;
  }
  return &it->second;
}

void SetElementOrigin(const mfem::ParMesh &mesh, std::vector<std::int64_t> &&origin)
{
  if (origin.empty())
  {
    ElementOrigins().erase(&mesh);
  }
  else
  {
    ElementOrigins()[&mesh] = std::move(origin);
  }
}

std::vector<std::array<std::int64_t, 2>> GetElementInputKeys(const mfem::ParMesh &mesh)
{
  std::vector<std::array<std::int64_t, 2>> keys;
  if (const std::vector<int> *origin = GetRootElementOrigin(mesh))
  {
    // The input index of the root element, then the path from the root element (the child
    // index at each level, 4 bits per level from the most significant ones, offset by 1),
    // whose order is the depth-first order of the leaves of a tetrahedral tree.
    const mfem::NCMesh &ncmesh = *mesh.ncmesh;
    const mfem::Array<int> &leaves = NCMeshLeaves::Get(ncmesh);
    keys.resize(mesh.GetNE());
    std::vector<int> path;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      int elem = leaves[e];
      path.clear();
      while (ncmesh.GetElement(elem).parent >= 0)
      {
        const int parent = ncmesh.GetElement(elem).parent;
        const auto &pa = ncmesh.GetElement(parent);
        int c = 0;
        while (pa.child[c] != elem)
        {
          c++;
        }
        path.push_back(c);
        elem = parent;
      }
      std::uint64_t key = 0;
      int shift = 60;
      for (auto it = path.rbegin(); it != path.rend() && shift >= 0; ++it, shift -= 4)
      {
        key |= static_cast<std::uint64_t>(*it + 1) << shift;
      }
      keys[e] = {(*origin)[elem], static_cast<std::int64_t>(key >> 1)};
    }
  }
  else if (const auto it = ElementOrigins().find(&mesh);
           !mesh.Nonconforming() && it != ElementOrigins().end() &&
           static_cast<int>(it->second.size()) == mesh.GetNE())
  {
    keys.resize(mesh.GetNE());
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      keys[e] = {it->second[e], 0};
    }
  }
  return keys;
}

int GetMaxBdrAttribute(const mfem::ParMesh &mesh)
{
  int bdr_attr_max = 0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    bdr_attr_max = std::max(bdr_attr_max, mesh.GetBdrAttribute(be));
  }
  Mpi::GlobalMax(1, &bdr_attr_max, mesh.GetComm());
  return bdr_attr_max;
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
    deviation_deg(i) =
        std::acos(std::min(1.0, std::abs(cosine))) * (180.0 / std::numbers::pi);
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
    auto p_min = std::ranges::min_element(vertices, Dot);
    auto p_max = std::ranges::max_element(vertices, Dot);
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
    auto p = max ? std::ranges::max_element(
                       vertices, [&x0](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                       { return (x - x0).norm() < (y - x0).norm(); })
                 : std::ranges::min_element(
                       vertices, [&x0](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
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

mfem::Vector GetSurfaceNormal(const mfem::Mesh &mesh, const mfem::Array<int> &marker,
                              bool average)
{
  // Serial counterpart of the ParMesh overload; no MPI. Used by ProjectSubmeshTo2D on
  // the pre-partition serial mesh.
  const int dim = mesh.SpaceDimension();
  mfem::IsoparametricTransformation T;
  mfem::Vector loc_normal(dim), normal(dim);
  normal = 0.0;
  if (mesh.Dimension() == mesh.SpaceDimension())
  {
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
  const double norm = normal.Norml2();
  if (norm > 0.0)
  {
    normal /= norm;
  }
  return normal;
}

namespace
{

// Trait helpers dispatching between serial (mfem::SubMesh / mfem::Mesh) and parallel
// (mfem::ParSubMesh / mfem::ParMesh) submesh primitives so a single templated body
// serves both. The serial versions return MPI_COMM_SELF; collectives over MPI_COMM_SELF
// degenerate to no-ops.
inline MPI_Comm GetSubMeshComm(const mfem::SubMesh &)
{
  return MPI_COMM_SELF;
}
inline MPI_Comm GetSubMeshComm(const mfem::ParSubMesh &s)
{
  return s.GetComm();
}

// Parent neighbor-transformation lookup. ParMesh delegates to the Palace helper that
// handles cross-rank shared faces; serial Mesh uses MFEM's face-transformation API
// directly (no cross-rank neighbors).
inline void GetBdrElementNeighborTransforms(int i, const mfem::Mesh &mesh,
                                            mfem::FaceElementTransformations &FET,
                                            mfem::IsoparametricTransformation &T1,
                                            mfem::IsoparametricTransformation &T2)
{
  int f, o;
  mesh.GetBdrElementFace(i, &f, &o);
  mesh.GetFaceElementTransformations(f, FET, T1, T2);
}
inline void GetBdrElementNeighborTransforms(int i, const mfem::ParMesh &mesh,
                                            mfem::FaceElementTransformations &FET,
                                            mfem::IsoparametricTransformation &T1,
                                            mfem::IsoparametricTransformation &T2)
{
  BdrGridFunctionCoefficient::GetBdrElementNeighborTransformations(i, mesh, FET, T1, T2);
}

}  // namespace

template <class SubMeshT>
void RemapSubMeshAttributes(SubMeshT &submesh)
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
    GetBdrElementNeighborTransforms(parent_bdr_elem, parent, FET, T1, T2);
    // Use the neighboring domain element's attribute. For internal boundaries, use the
    // element with the lower attribute (same convention as Palace's mesh.cpp).
    int nbr_attr = FET.Elem1->Attribute;
    if (FET.Elem2 && FET.Elem2->Attribute < nbr_attr)
    {
      nbr_attr = FET.Elem2->Attribute;
    }
    submesh.SetAttribute(i, nbr_attr);
  }
  // Refresh the cached distinct-attribute array to match the rewritten per-element
  // attributes.
  submesh.SetAttributes();
}

namespace
{

// The submesh boundary edges and parent boundary face edges matched below can be owned by
// different processes, and on a nonconforming mesh one can contain the other rather than
// equal it. Neither the ParMesh vertex or edge numbering, nor the NCMesh node ids, nor the
// slave/master classification of edges (which only reflects the elements a process sees)
// identify an edge consistently across processes after parallel refinement and
// rebalancing, so edges are matched geometrically: an edge contains another if the
// endpoints of the other lie on it. This is tested on the process which owns the
// containing edge, where its exact (possibly curved) geometry is available, against the
// endpoints of the edges of all processes. Each matching is done in both directions, once
// with each set of edges gathered.

// Tolerance of the containment test, relative to the edge lengths. The two sides of a
// cracked internal boundary are only separated by CrackDisplacementFactor (by default
// 1.0e-12) times the mesh size, so they match as one.
constexpr double EDGE_REL_TOL = 1.0e-6;

// The exact geometry of a mesh edge, for testing whether points lie on it.
class EdgeCurve
{
private:
  int sdim, order;
  mutable mfem::IsoparametricTransformation T;
  mutable mfem::IntegrationPoint ip;
  // Points at the parameters k / (2 order), k = 0, ..., 2 order.
  std::vector<double> samples;

  void Eval(double t, double *x) const
  {
    ip.x = t;
    mfem::Vector xv(x, sdim);
    T.Transform(ip, xv);
  }

public:
  std::array<double, 3> bb_min = {0.0, 0.0, 0.0}, bb_max = {0.0, 0.0, 0.0};
  double length = 0.0;  // Distance between the endpoints

  EdgeCurve(const mfem::Mesh &mesh, int edge) : sdim(mesh.SpaceDimension())
  {
    // Use the edge transformation rather than the vertex coordinates, which are not kept
    // up to date with the nodes of a mesh with a nodal grid function.
    mesh.GetEdgeTransformation(edge, &T);
    order = std::max(T.GetFE()->GetOrder(), 1);
    const int n = 2 * order + 1;
    samples.resize(static_cast<std::size_t>(n) * sdim);
    for (int k = 0; k < n; k++)
    {
      Eval(static_cast<double>(k) / (n - 1), samples.data() + k * sdim);
    }
    double l2 = 0.0;
    for (int d = 0; d < sdim; d++)
    {
      const double c = samples[(n - 1) * sdim + d] - samples[d];
      l2 += c * c;
    }
    length = std::sqrt(l2);
    // The bounding box of the samples, enlarged to contain the curve between them.
    const double margin = (order > 1 ? 0.25 : EDGE_REL_TOL) * length;
    for (int d = 0; d < sdim; d++)
    {
      bb_min[d] = mfem::infinity();
      bb_max[d] = -mfem::infinity();
      for (int k = 0; k < n; k++)
      {
        bb_min[d] = std::min(bb_min[d], samples[k * sdim + d]);
        bb_max[d] = std::max(bb_max[d], samples[k * sdim + d]);
      }
      bb_min[d] -= margin;
      bb_max[d] += margin;
    }
  }

  const double *Point(int k) const { return samples.data() + k * sdim; }
  const double *Endpoint(int k) const { return Point(k == 0 ? 0 : 2 * order); }

  // Distance from a point to the curve: closest sample, then Gauss-Newton on the edge
  // parameter, which converges in one step for straight edges.
  double Distance(const double *p) const
  {
    const int n = 2 * order + 1;
    auto Dist2 = [this, p](const double *x)
    {
      double d2 = 0.0;
      for (int d = 0; d < sdim; d++)
      {
        d2 += (x[d] - p[d]) * (x[d] - p[d]);
      }
      return d2;
    };
    int k0 = 0;
    for (int k = 1; k < n; k++)
    {
      if (Dist2(Point(k)) < Dist2(Point(k0)))
      {
        k0 = k;
      }
    }
    double t = static_cast<double>(k0) / (n - 1);
    std::array<double, 3> x;
    for (int it = 0; it < 20; it++)
    {
      Eval(t, x.data());
      T.SetIntPoint(&ip);
      const auto &J = T.Jacobian();
      double g = 0.0, h = 0.0;
      for (int d = 0; d < sdim; d++)
      {
        g += J(d, 0) * (x[d] - p[d]);
        h += J(d, 0) * J(d, 0);
      }
      if (h <= 0.0)
      {
        break;
      }
      const double t_new = std::clamp(t - g / h, 0.0, 1.0);
      const bool converged = std::abs(t_new - t) < 1.0e-14;
      t = t_new;
      if (converged)
      {
        break;
      }
    }
    Eval(t, x.data());
    return std::sqrt(Dist2(x.data()));
  }

  // Whether the edge with the given endpoints and length lies on this edge.
  bool Contains(const double *p0, const double *p1, double other_length) const
  {
    const double tol = EDGE_REL_TOL * std::max(length, other_length);
    return Distance(p0) <= tol && Distance(p1) <= tol;
  }
};

// A set of edges given by their endpoints (two points of the space dimension each), with
// a search for the edges whose bounding box overlaps a given box. The bounding boxes are
// those of the endpoints, unless given (min and max corners of each edge in turn).
class EdgeSearch
{
private:
  int sdim, axis = 0;
  std::vector<double> pts, bb_min, bb_max, length, sorted_min;
  std::vector<int> order;
  double max_extent = 0.0;

public:
  // Global bounding box of the endpoints, and a tolerance for box tests relative to its
  // size.
  std::vector<double> glob_min, glob_max;
  double bb_tol = 0.0;

  EdgeSearch(std::vector<double> &&points, int sdim, const std::vector<double> &boxes = {})
    : sdim(sdim), pts(std::move(points)), glob_min(sdim, mfem::infinity()),
      glob_max(sdim, -mfem::infinity())
  {
    const int n = Size();
    MFEM_VERIFY(boxes.empty() || boxes.size() == static_cast<std::size_t>(n) * 2 * sdim,
                "Invalid edge bounding boxes!");
    bb_min.resize(static_cast<std::size_t>(n) * sdim);
    bb_max.resize(static_cast<std::size_t>(n) * sdim);
    length.resize(n);
    for (int i = 0; i < n; i++)
    {
      const double *p = Endpoint(i, 0), *q = Endpoint(i, 1);
      double l2 = 0.0;
      for (int d = 0; d < sdim; d++)
      {
        bb_min[i * sdim + d] =
            boxes.empty() ? std::min(p[d], q[d]) : boxes[2 * i * sdim + d];
        bb_max[i * sdim + d] =
            boxes.empty() ? std::max(p[d], q[d]) : boxes[(2 * i + 1) * sdim + d];
        glob_min[d] = std::min(glob_min[d], bb_min[i * sdim + d]);
        glob_max[d] = std::max(glob_max[d], bb_max[i * sdim + d]);
        l2 += (q[d] - p[d]) * (q[d] - p[d]);
      }
      length[i] = std::sqrt(l2);
    }
    double extent = 0.0;
    for (int d = 0; d < sdim && n > 0; d++)
    {
      if (glob_max[d] - glob_min[d] > extent)
      {
        extent = glob_max[d] - glob_min[d];
        axis = d;
      }
    }
    bb_tol = 1.0e-8 * extent;
    order.resize(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int a, int b)
              { return bb_min[a * this->sdim + axis] < bb_min[b * this->sdim + axis]; });
    sorted_min.resize(n);
    for (int j = 0; j < n; j++)
    {
      sorted_min[j] = bb_min[order[j] * sdim + axis];
      max_extent = std::max(max_extent, bb_max[order[j] * sdim + axis] - sorted_min[j]);
    }
  }

  int Size() const { return static_cast<int>(pts.size()) / (2 * sdim); }
  const double *Endpoint(int i, int k) const
  {
    return pts.data() + (static_cast<std::size_t>(i) * 2 + k) * sdim;
  }
  double Length(int i) const { return length[i]; }

  // Whether a box overlaps the global bounding box.
  bool Overlaps(const double *q_min, const double *q_max) const
  {
    for (int d = 0; d < sdim; d++)
    {
      if (q_max[d] < glob_min[d] - bb_tol || q_min[d] > glob_max[d] + bb_tol)
      {
        return false;
      }
    }
    return Size() > 0;
  }

  // Call f(i) for the edges whose endpoints' bounding box overlaps the given box.
  template <typename F>
  void ForEachCandidate(const double *q_min, const double *q_max, F &&f) const
  {
    auto it = std::lower_bound(sorted_min.begin(), sorted_min.end(),
                               q_min[axis] - max_extent - bb_tol);
    for (auto j = static_cast<int>(it - sorted_min.begin());
         j < Size() && sorted_min[j] <= q_max[axis] + bb_tol; j++)
    {
      const int i = order[j];
      bool overlap = true;
      for (int d = 0; d < sdim && overlap; d++)
      {
        overlap = bb_max[i * sdim + d] >= q_min[d] - bb_tol &&
                  bb_min[i * sdim + d] <= q_max[d] + bb_tol;
      }
      if (overlap)
      {
        f(i);
      }
    }
  }
};

// Gather the given values from all processes, returning them concatenated in order of rank
// and the offsets of the values of each rank.
template <typename T>
std::vector<T> GatherAll(const std::vector<T> &local, MPI_Comm comm,
                         std::vector<int> *offsets = nullptr)
{
  const int size = Mpi::Size(comm);
  const int count = static_cast<int>(local.size());
  std::vector<int> counts(size), displs(size);
  Mpi::Allgather(1, &count, counts.data(), comm);
  int total = 0;
  for (int r = 0; r < size; r++)
  {
    displs[r] = total;
    total += counts[r];
  }
  std::vector<T> all(total);
  Mpi::Allgatherv(count, local.data(), all.data(), counts.data(), displs.data(), comm);
  if (offsets)
  {
    offsets->assign(displs.begin(), displs.end());
    offsets->push_back(total);
  }
  return all;
}

// The curve (in parent coordinates) of a submesh edge: the submesh coordinates may already
// have been projected to 2D.
template <class SubMeshT>
EdgeCurve GetSubMeshEdgeCurve(const SubMeshT &submesh, int submesh_edge)
{
  const mfem::Array<int> &parent_edge_map = submesh.GetParentEdgeIDMap();
  MFEM_VERIFY(submesh_edge >= 0 && submesh_edge < parent_edge_map.Size() &&
                  parent_edge_map[submesh_edge] >= 0,
              "Submesh edge has no parent edge!");
  return EdgeCurve(*submesh.GetParent(), parent_edge_map[submesh_edge]);
}

void AppendEndpoints(const EdgeCurve &curve, int sdim, std::vector<double> &pts)
{
  for (int k = 0; k < 2; k++)
  {
    pts.insert(pts.end(), curve.Endpoint(k), curve.Endpoint(k) + sdim);
  }
}

void AppendBoundingBox(const EdgeCurve &curve, int sdim, std::vector<double> &boxes)
{
  boxes.insert(boxes.end(), curve.bb_min.begin(), curve.bb_min.begin() + sdim);
  boxes.insert(boxes.end(), curve.bb_max.begin(), curve.bb_max.begin() + sdim);
}

// Whether both endpoints of an edge lie in a box, within a tolerance.
bool EndpointsInside(const EdgeCurve &curve, int sdim, const std::vector<double> &q_min,
                     const std::vector<double> &q_max, double tol)
{
  for (int k = 0; k < 2; k++)
  {
    for (int d = 0; d < sdim; d++)
    {
      const double x = curve.Endpoint(k)[d];
      if (x < q_min[d] - tol || x > q_max[d] + tol)
      {
        return false;
      }
    }
  }
  return true;
}

// The local edges of the parent boundary faces with the given attributes, as (edge,
// attribute) pairs with each pair listed once.
template <typename AttrFilter>
std::vector<std::pair<int, int>> GetParentBoundaryEdges(const mfem::Mesh &parent,
                                                        AttrFilter &&keep)
{
  std::set<std::pair<int, int>> edge_attrs;
  mfem::Array<int> edges, orientations;
  for (int be = 0; be < parent.GetNBE(); be++)
  {
    const int attr = parent.GetBdrAttribute(be);
    if (!keep(attr))
    {
      continue;
    }
    parent.GetBdrElementEdges(be, edges, orientations);
    for (int edge : edges)
    {
      edge_attrs.insert({edge, attr});
    }
  }
  return {edge_attrs.begin(), edge_attrs.end()};
}

}  // namespace

template <class SubMeshT>
void RemapSubMeshBdrAttributes(SubMeshT &submesh, const mfem::Array<int> &surface_attrs)
{
  MFEM_VERIFY(submesh.GetFrom() == mfem::SubMesh::From::Boundary,
              "RemapSubMeshBdrAttributes requires a boundary submesh!");

  // Give each submesh boundary (perimeter) edge the attribute of the parent boundary face
  // it lies on, other than the selected surface itself. If several physical boundaries
  // meet at the same edge, choose the smallest attribute deterministically; an edge with
  // no adjacent non-surface parent boundary keeps MFEM's generated attribute (it can be an
  // artificial NC submesh boundary). The edges are matched geometrically, see above.
  const auto &parent = *submesh.GetParent();
  MPI_Comm comm = GetSubMeshComm(submesh);
  const int sdim = parent.SpaceDimension();
  std::unordered_set<int> surface_attr_set(surface_attrs.begin(), surface_attrs.end());

  // Gather the endpoints and bounding boxes of the submesh boundary edges.
  const int nbe = submesh.GetNBE();
  std::vector<double> local_sub_pts, local_sub_boxes;
  local_sub_pts.reserve(static_cast<std::size_t>(nbe) * 2 * sdim);
  local_sub_boxes.reserve(static_cast<std::size_t>(nbe) * 2 * sdim);
  for (int sbe = 0; sbe < nbe; sbe++)
  {
    const EdgeCurve sub = GetSubMeshEdgeCurve(submesh, submesh.GetBdrElementFaceIndex(sbe));
    AppendEndpoints(sub, sdim, local_sub_pts);
    AppendBoundingBox(sub, sdim, local_sub_boxes);
  }
  std::vector<int> sub_offsets;
  const EdgeSearch sub_edges(GatherAll(local_sub_pts, comm, &sub_offsets), sdim,
                             GatherAll(local_sub_boxes, comm));
  const int nedges = sub_edges.Size();
  if (nedges == 0)
  {
    return;
  }
  const int offset = sub_offsets[Mpi::Rank(comm)] / (2 * sdim);
  std::vector<int> attr(nedges, std::numeric_limits<int>::max());

  // Submesh edges lying on local parent boundary face edges, and gather the endpoints of
  // the parent boundary face edges lying on the surface for the converse.
  std::vector<double> local_wall_pts;
  std::vector<int> local_wall_attrs;
  for (const auto &edge_attr :
       GetParentBoundaryEdges(parent, [&](int a) { return !surface_attr_set.contains(a); }))
  {
    // Not a structured binding, which OpenMP does not support capturing in the lambda.
    const int be_attr = edge_attr.second;
    const EdgeCurve wall(parent, edge_attr.first);
    if (!sub_edges.Overlaps(wall.bb_min.data(), wall.bb_max.data()))
    {
      continue;
    }
    sub_edges.ForEachCandidate(wall.bb_min.data(), wall.bb_max.data(),
                               [&](int i)
                               {
                                 if (be_attr < attr[i] &&
                                     wall.Contains(sub_edges.Endpoint(i, 0),
                                                   sub_edges.Endpoint(i, 1),
                                                   sub_edges.Length(i)))
                                 {
                                   attr[i] = be_attr;
                                 }
                               });
    if (EndpointsInside(wall, sdim, sub_edges.glob_min, sub_edges.glob_max,
                        sub_edges.bb_tol))
    {
      AppendEndpoints(wall, sdim, local_wall_pts);
      local_wall_attrs.push_back(be_attr);
    }
  }

  // Parent boundary face edges lying on local submesh boundary edges.
  const std::vector<int> wall_attrs = GatherAll(local_wall_attrs, comm);
  const EdgeSearch walls(GatherAll(local_wall_pts, comm), sdim);
  for (int sbe = 0; sbe < nbe && walls.Size() > 0; sbe++)
  {
    const EdgeCurve sub = GetSubMeshEdgeCurve(submesh, submesh.GetBdrElementFaceIndex(sbe));
    int &a = attr[offset + sbe];
    walls.ForEachCandidate(
        sub.bb_min.data(), sub.bb_max.data(),
        [&](int i)
        {
          if (wall_attrs[i] < a &&
              sub.Contains(walls.Endpoint(i, 0), walls.Endpoint(i, 1), walls.Length(i)))
          {
            a = wall_attrs[i];
          }
        });
  }
  Mpi::GlobalMin(nedges, attr.data(), comm);

  for (int sbe = 0; sbe < nbe; sbe++)
  {
    if (attr[offset + sbe] != std::numeric_limits<int>::max())
    {
      submesh.SetBdrAttribute(sbe, attr[offset + sbe]);
    }
  }

  // Do not touch submesh.bdr_attributes directly: modifying it on a ParSubMesh can
  // corrupt internal MFEM state. RebuildCeedAttributes reads per-element attributes.
}

template <class SubMeshT>
void AddSubMeshInternalBoundaryElements(SubMeshT &submesh,
                                        const mfem::Array<int> &surface_attrs,
                                        const std::vector<int> &internal_bdr_attrs)
{
  MFEM_VERIFY(submesh.GetFrom() == mfem::SubMesh::From::Boundary,
              "AddSubMeshInternalBoundaryElements requires a boundary submesh!");
  if (internal_bdr_attrs.empty())
  {
    return;
  }

  // Add a submesh boundary element on each interior submesh edge lying on a parent
  // boundary face with one of the given attributes (for example an uncracked PEC sheet
  // crossing the surface), with the smallest such attribute. The edges are matched
  // geometrically, see above. On a nonconforming submesh, boundary elements are only
  // added on conforming and slave edges, not on master edges (whose essential degrees of
  // freedom follow from their slaves). An edge shared by several processes is added by
  // the lowest of them, so that it is only integrated over once.
  const auto &parent = *submesh.GetParent();
  MPI_Comm comm = GetSubMeshComm(submesh);
  const int rank = Mpi::Rank(comm);
  const int sdim = parent.SpaceDimension();
  std::unordered_set<int> surface_attr_set(surface_attrs.begin(), surface_attrs.end());
  std::unordered_set<int> internal_attr_set;
  for (int a : internal_bdr_attrs)
  {
    if (!surface_attr_set.contains(a))
    {
      internal_attr_set.insert(a);
    }
  }

  // Candidate submesh edges: interior (not already boundary), non-master edges.
  std::unordered_set<int> bdr_edges;
  for (int sbe = 0; sbe < submesh.GetNBE(); sbe++)
  {
    bdr_edges.insert(submesh.GetBdrElementFaceIndex(sbe));
  }
  std::vector<int> cand_edges;
  std::vector<double> glob_min(sdim, mfem::infinity()), glob_max(sdim, -mfem::infinity());
  for (int e = 0; e < submesh.GetNEdges(); e++)
  {
    const EdgeCurve curve = GetSubMeshEdgeCurve(submesh, e);
    for (int d = 0; d < sdim; d++)
    {
      glob_min[d] = std::min(glob_min[d], curve.bb_min[d]);
      glob_max[d] = std::max(glob_max[d], curve.bb_max[d]);
    }
    if (bdr_edges.contains(e) ||
        (submesh.Nonconforming() && submesh.GetFaceInformation(e).IsNonconformingCoarse()))
    {
      continue;
    }
    cand_edges.push_back(e);
  }
  Mpi::GlobalMin(sdim, glob_min.data(), comm);
  Mpi::GlobalMax(sdim, glob_max.data(), comm);
  double extent = 0.0;
  for (int d = 0; d < sdim; d++)
  {
    extent = std::max(extent, glob_max[d] - glob_min[d]);
  }
  const double bb_tol = 1.0e-8 * extent;

  // Gather the endpoints of the internal boundary face edges lying on the surface.
  const auto internal_edges =
      GetParentBoundaryEdges(parent, [&](int a) { return internal_attr_set.contains(a); });
  std::vector<double> local_wall_pts, local_wall_boxes;
  std::vector<int> local_wall_attrs;
  for (const auto &[edge, be_attr] : internal_edges)
  {
    const EdgeCurve wall(parent, edge);
    if (EndpointsInside(wall, sdim, glob_min, glob_max, bb_tol))
    {
      AppendEndpoints(wall, sdim, local_wall_pts);
      AppendBoundingBox(wall, sdim, local_wall_boxes);
      local_wall_attrs.push_back(be_attr);
    }
  }
  const std::vector<int> wall_attrs = GatherAll(local_wall_attrs, comm);
  const EdgeSearch walls(GatherAll(local_wall_pts, comm), sdim,
                         GatherAll(local_wall_boxes, comm));
  if (walls.Size() == 0)
  {
    return;
  }

  // Internal boundary edges lying on local candidate edges, and gather the endpoints of the
  // candidate edges near internal boundary edges (whose bounding boxes overlap) for the
  // converse.
  std::vector<int> match_attr(cand_edges.size(), std::numeric_limits<int>::max());
  std::vector<double> local_near_pts;
  std::vector<int> local_near;
  std::vector<double> q_min(sdim), q_max(sdim);
  for (std::size_t j = 0; j < cand_edges.size(); j++)
  {
    const EdgeCurve cand = GetSubMeshEdgeCurve(submesh, cand_edges[j]);
    bool near = false;
    walls.ForEachCandidate(
        cand.bb_min.data(), cand.bb_max.data(),
        [&](int i)
        {
          near = true;
          if (wall_attrs[i] < match_attr[j] &&
              cand.Contains(walls.Endpoint(i, 0), walls.Endpoint(i, 1), walls.Length(i)))
          {
            match_attr[j] = wall_attrs[i];
          }
        });
    if (near)
    {
      AppendEndpoints(cand, sdim, local_near_pts);
      local_near.push_back(static_cast<int>(j));
    }
  }

  // Candidate edges lying on local internal boundary edges.
  std::vector<int> near_offsets;
  const EdgeSearch near_edges(GatherAll(local_near_pts, comm, &near_offsets), sdim);
  std::vector<int> near_attr(near_edges.Size(), std::numeric_limits<int>::max());
  for (const auto &edge_attr : internal_edges)
  {
    // Not a structured binding, which OpenMP does not support capturing in the lambda.
    const int be_attr = edge_attr.second;
    const EdgeCurve wall(parent, edge_attr.first);
    if (!near_edges.Overlaps(wall.bb_min.data(), wall.bb_max.data()))
    {
      continue;
    }
    near_edges.ForEachCandidate(wall.bb_min.data(), wall.bb_max.data(),
                                [&](int i)
                                {
                                  if (be_attr < near_attr[i] &&
                                      wall.Contains(near_edges.Endpoint(i, 0),
                                                    near_edges.Endpoint(i, 1),
                                                    near_edges.Length(i)))
                                  {
                                    near_attr[i] = be_attr;
                                  }
                                });
  }
  Mpi::GlobalMin(static_cast<int>(near_attr.size()), near_attr.data(), comm);
  const int near_offset = near_offsets[rank] / (2 * sdim);
  for (std::size_t k = 0; k < local_near.size(); k++)
  {
    match_attr[local_near[k]] =
        std::min(match_attr[local_near[k]], near_attr[near_offset + k]);
  }

  // Each edge is added by the lowest process on which it matched.
  std::vector<double> local_match_pts;
  std::vector<int> local_matches;
  for (std::size_t j = 0; j < cand_edges.size(); j++)
  {
    if (match_attr[j] != std::numeric_limits<int>::max())
    {
      local_matches.push_back(static_cast<int>(j));
      AppendEndpoints(GetSubMeshEdgeCurve(submesh, cand_edges[j]), sdim, local_match_pts);
    }
  }
  std::vector<int> match_offsets;
  const EdgeSearch all_matches(GatherAll(local_match_pts, comm, &match_offsets), sdim);
  std::vector<int> match_rank(all_matches.Size());
  for (int r = 0; r < Mpi::Size(comm); r++)
  {
    for (int i = match_offsets[r] / (2 * sdim); i < match_offsets[r + 1] / (2 * sdim); i++)
    {
      match_rank[i] = r;
    }
  }
  mfem::Array<mfem::Element *> new_bdr_elems;
  mfem::Array<int> new_be_to_face;
  mfem::Array<int> edge_verts;
  const int match_offset = match_offsets[rank] / (2 * sdim);
  for (std::size_t k = 0; k < local_matches.size(); k++)
  {
    const int j = local_matches[k];
    const int self = match_offset + static_cast<int>(k);
    const double *p0 = all_matches.Endpoint(self, 0), *p1 = all_matches.Endpoint(self, 1);
    const double tol = EDGE_REL_TOL * all_matches.Length(self);
    for (int d = 0; d < sdim; d++)
    {
      q_min[d] = std::min(p0[d], p1[d]);
      q_max[d] = std::max(p0[d], p1[d]);
    }
    bool owner = true;
    all_matches.ForEachCandidate(
        q_min.data(), q_max.data(),
        [&](int i)
        {
          if (match_rank[i] >= rank || !owner)
          {
            return;
          }
          // The same edge, with its endpoints in either order.
          auto Close = [&](const double *a, const double *b)
          {
            double d2 = 0.0;
            for (int d = 0; d < sdim; d++)
            {
              d2 += (a[d] - b[d]) * (a[d] - b[d]);
            }
            return std::sqrt(d2) <= tol;
          };
          const double *q0 = all_matches.Endpoint(i, 0), *q1 = all_matches.Endpoint(i, 1);
          owner = !((Close(p0, q0) && Close(p1, q1)) || (Close(p0, q1) && Close(p1, q0)));
        });
    if (!owner)
    {
      continue;
    }
    submesh.GetEdgeVertices(cand_edges[j], edge_verts);
    new_bdr_elems.Append(new mfem::Segment(edge_verts[0], edge_verts[1], match_attr[j]));
    new_be_to_face.Append(cand_edges[j]);
  }
  if (new_bdr_elems.Size() > 0)
  {
    // AddBdrElements adds the elements with their topology mapping, avoiding
    // FinalizeTopology (which deadlocks on a ParSubMesh).
    submesh.AddBdrElements(new_bdr_elems, new_be_to_face);
  }
}

mfem::Vector ProjectSubmeshTo2D(mfem::Mesh &submesh, mfem::Vector &centroid,
                                mfem::Vector &e1, mfem::Vector &e2)
{
  MFEM_VERIFY(submesh.Dimension() == 2 && submesh.SpaceDimension() == 3,
              "ProjectSubmeshTo2D requires a 2D submesh with 3D ambient coordinates!");

  // Serial operation on the pre-partitioned mesh: compute the normal locally (no MPI
  // reductions), build the tangent frame into the caller-owned outputs (centroid, e1,
  // e2), then rewrite node coordinates in place.
  mfem::Vector normal = GetSurfaceNormal(submesh);

  // Build an orthonormal tangent frame (e1, e2) perpendicular to the normal.
  // Choose e1 as the cross product of normal with the axis most perpendicular to it.
  e1.SetSize(3);
  e2.SetSize(3);
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

  // Compute the centroid of all submesh vertices. Serial — no MPI reduction needed.
  centroid.SetSize(3);
  centroid = 0.0;
  const int nv = submesh.GetNV();
  for (int i = 0; i < nv; i++)
  {
    const double *v = submesh.GetVertex(i);
    for (int d = 0; d < 3; d++)
    {
      centroid(d) += v[d];
    }
  }
  if (nv > 0)
  {
    centroid /= nv;
  }

  // Read the 3D node coordinates, project to 2D, then rebuild the mesh with 2D nodes.
  // We read all coordinates first, then use SetCurvature to rebuild with SpaceDim=2.
  {
    // Determine the mesh curvature order from the existing nodes.
    int mesh_order = 1;
    if (submesh.GetNodes() && submesh.GetNE() > 0)
    {
      mesh_order = submesh.GetNodes()->FESpace()->GetMaxElementOrder();
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
      projected.resize(nv);
      for (int i = 0; i < nv; i++)
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
    // GridFunction with the specified vdim (SpaceDim). Extraction runs on a rank that
    // holds the serial mesh, so submesh.GetNE() > 0 is required.
    MFEM_VERIFY(submesh.GetNE() > 0,
                "ProjectSubmeshTo2D called on an empty mesh — extraction must run on a "
                "rank that holds the serial mesh!");
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

  MFEM_VERIFY(submesh.SpaceDimension() == 2,
              "ProjectSubmeshTo2D: SpaceDimension should be 2 after projection!");
  return normal;
}

std::unique_ptr<mfem::Mesh> RebuildNonconformingMesh(const mfem::Mesh &orig)
{
  MFEM_VERIFY(orig.Nonconforming() && orig.Dimension() == 2 &&
                  !dynamic_cast<const mfem::ParMesh *>(&orig),
              "RebuildNonconformingMesh requires a serial 2D nonconforming mesh!");
  const int sdim = orig.SpaceDimension();
  auto mesh =
      std::make_unique<mfem::Mesh>(2, orig.GetNV(), orig.GetNE(), orig.GetNBE(), sdim);

  // Vertex coordinates from the nodes: the vertex coordinates of a nonconforming SubMesh
  // do not follow its vertex numbering.
  for (int i = 0; i < orig.GetNV(); i++)
  {
    double x[3] = {0.0, 0.0, 0.0};
    orig.GetNode(i, x);
    mesh->AddVertex(x);
  }
  for (int i = 0; i < orig.GetNE(); i++)
  {
    mesh->AddElement(orig.GetElement(i)->Duplicate(mesh.get()));
  }
  for (int i = 0; i < orig.GetNBE(); i++)
  {
    mesh->AddBdrElement(orig.GetBdrElement(i)->Duplicate(mesh.get()));
  }

  // MFEM keeps the attribute of a boundary edge refined on one side on its master edge
  // (which has no boundary element), and passes it on to the children of the master when
  // the elements on its unrefined side are refined. Add the missing boundary elements of
  // the master and slave edges of a master with a boundary element on any of them: the
  // nonconforming mesh keeps their attributes, but only the boundary elements of the
  // slaves.
  auto &edge_list = orig.ncmesh->GetEdgeList();
  mfem::Array<int> v;
  {
    std::unordered_map<int, int> edge_attr;
    for (int i = 0; i < orig.GetNBE(); i++)
    {
      edge_attr[orig.GetBdrElementFaceIndex(i)] = orig.GetBdrAttribute(i);
    }
    for (const auto &master : edge_list.masters)
    {
      int attr = 0;
      auto Visit = [&](int edge)
      {
        auto it = edge_attr.find(edge);
        if (it != edge_attr.end())
        {
          MFEM_VERIFY(attr == 0 || attr == it->second,
                      "Different boundary attributes on the edges of a nonconforming "
                      "master edge!");
          attr = it->second;
        }
      };
      Visit(master.index);
      for (int s = master.slaves_begin; s < master.slaves_end; s++)
      {
        Visit(edge_list.slaves[s].index);
      }
      if (attr == 0)
      {
        continue;
      }
      auto Add = [&](int edge)
      {
        if (edge_attr.count(edge) == 0)
        {
          orig.GetEdgeVertices(edge, v);
          mesh->AddBdrSegment(v[0], v[1], attr);
        }
      };
      Add(master.index);
      for (int s = master.slaves_begin; s < master.slaves_end; s++)
      {
        Add(edge_list.slaves[s].index);
      }
    }
  }

  // The hanging vertices are the vertices of the slaves of a master edge in its interior.
  // Their positions along the master edge follow from the reference coordinates of the
  // slaves (exact dyadic numbers), starting from the endpoints of the master in the
  // orientation for which they are consistent. A vertex at k / 2^l (k odd) bisects the
  // edge between the vertices at (k - 1) / 2^l and (k + 1) / 2^l.
  constexpr int max_level = 30;
  constexpr double tol = 1.0e-12;
  const auto &point_matrices = edge_list.point_matrices[mfem::Geometry::SEGMENT];
  std::vector<std::array<int, 4>> parents;  // Level, vertex, parent vertices
  for (const auto &master : edge_list.masters)
  {
    orig.GetEdgeVertices(master.index, v);
    const int a = v[0], b = v[1];
    std::unordered_map<int, double> pos;
    bool consistent = false;
    for (int orientation = 0; orientation < 2 && !consistent; orientation++)
    {
      pos = {{a, orientation ? 1.0 : 0.0}, {b, orientation ? 0.0 : 1.0}};
      consistent = true;
      for (bool changed = true; changed && consistent;)
      {
        changed = false;
        for (int s = master.slaves_begin; s < master.slaves_end && consistent; s++)
        {
          const auto &slave = edge_list.slaves[s];
          const auto &pm = *point_matrices[slave.matrix];
          orig.GetEdgeVertices(slave.index, v);
          auto it0 = pos.find(v[0]), it1 = pos.find(v[1]);
          if (it0 == pos.end() && it1 == pos.end())
          {
            continue;
          }
          const int known = (it0 != pos.end()) ? 0 : 1;
          const double t = (known == 0) ? it0->second : it1->second;
          const double t_other = (std::abs(t - pm(0, 0)) < tol)   ? pm(0, 1)
                                 : (std::abs(t - pm(0, 1)) < tol) ? pm(0, 0)
                                                                  : -1.0;
          if (t_other < 0.0)
          {
            consistent = false;
          }
          else if (it0 == pos.end() || it1 == pos.end())
          {
            pos[v[1 - known]] = t_other;
            changed = true;
          }
          else if (std::abs(pos[v[1 - known]] - t_other) >= tol)
          {
            consistent = false;
          }
        }
      }
    }
    MFEM_VERIFY(consistent, "Inconsistent slave edges of a nonconforming master edge!");

    std::unordered_map<long long, int> vertex_at;
    for (const auto &[i, t] : pos)
    {
      vertex_at[std::llround(t * (1 << max_level))] = i;
    }
    for (const auto &[i, t] : pos)
    {
      if (i == a || i == b)
      {
        continue;
      }
      const long long k = std::llround(t * (1 << max_level));
      MFEM_VERIFY(std::abs(t * (1 << max_level) - k) < 1.0e-6 && k > 0,
                  "Nonconforming master edge is not bisected!");
      int level = max_level;
      for (long long j = k; j % 2 == 0; j /= 2)
      {
        level--;
      }
      const long long dk = 1LL << (max_level - level);
      auto p1 = vertex_at.find(k - dk), p2 = vertex_at.find(k + dk);
      MFEM_VERIFY(p1 != vertex_at.end() && p2 != vertex_at.end(),
                  "Missing parent vertex for a hanging vertex of a nonconforming mesh!");
      parents.push_back({level, i, p1->second, p2->second});
    }
  }
  std::sort(parents.begin(), parents.end());
  for (const auto &[level, i, p1, p2] : parents)
  {
    mesh->AddVertexParents(i, p1, p2);
  }

  // The nonconforming mesh is built from the elements, boundary elements and hanging
  // vertices. It keeps the order of the elements and of their vertices, so the nodes can be
  // copied element by element.
  mesh->FinalizeTopology();
  mesh->Finalize();
  if (const auto *nodes = orig.GetNodes())
  {
    const auto &fes = *nodes->FESpace();
    mesh->SetCurvature(fes.GetMaxElementOrder(), fes.IsDGSpace(), sdim, fes.GetOrdering());
    auto &new_nodes = *mesh->GetNodes();
    mfem::Array<int> vdofs, new_vdofs;
    mfem::Vector vals;
    for (int i = 0; i < orig.GetNE(); i++)
    {
      fes.GetElementVDofs(i, vdofs);
      new_nodes.FESpace()->GetElementVDofs(i, new_vdofs);
      nodes->GetSubVector(vdofs, vals);
      new_nodes.SetSubVector(new_vdofs, vals);
    }
  }
  MFEM_VERIFY(mesh->GetNE() == orig.GetNE(),
              "Unexpected number of elements after rebuilding a nonconforming mesh!");
  for (int i = 0; i < orig.GetNE(); i++)
  {
    mfem::Array<int> ov, nv;
    orig.GetElementVertices(i, ov);
    mesh->GetElementVertices(i, nv);
    MFEM_VERIFY(ov.Size() == nv.Size(),
                "Unexpected element after rebuilding a nonconforming mesh!");
    for (int j = 0; j < ov.Size(); j++)
    {
      double xo[3] = {0.0, 0.0, 0.0}, xn[3] = {0.0, 0.0, 0.0};
      orig.GetNode(ov[j], xo);
      mesh->GetNode(nv[j], xn);
      MFEM_VERIFY(std::abs(xo[0] - xn[0]) + std::abs(xo[1] - xn[1]) +
                          std::abs(xo[2] - xn[2]) <=
                      1.0e-12 * (1.0 + std::abs(xo[0]) + std::abs(xo[1]) + std::abs(xo[2])),
                  "Unexpected element vertices after rebuilding a nonconforming mesh!");
    }
  }
  return mesh;
}

// Explicit instantiations for the submesh helpers. Serial (mfem::SubMesh) is used by
// BoundaryModeSolver on the pre-partition mesh; parallel (mfem::ParSubMesh) by WavePort
// after partitioning.
template void RemapSubMeshAttributes<mfem::SubMesh>(mfem::SubMesh &);
template void RemapSubMeshAttributes<mfem::ParSubMesh>(mfem::ParSubMesh &);
template void RemapSubMeshBdrAttributes<mfem::SubMesh>(mfem::SubMesh &,
                                                       const mfem::Array<int> &);
template void RemapSubMeshBdrAttributes<mfem::ParSubMesh>(mfem::ParSubMesh &,
                                                          const mfem::Array<int> &);
template void AddSubMeshInternalBoundaryElements<mfem::SubMesh>(mfem::SubMesh &,
                                                                const mfem::Array<int> &,
                                                                const std::vector<int> &);
template void AddSubMeshInternalBoundaryElements<mfem::ParSubMesh>(
    mfem::ParSubMesh &, const mfem::Array<int> &, const std::vector<int> &);

namespace
{

// Surface area (bdr == true) or volume of the mesh regions selected by marker.
double GetMeasure(const mfem::ParMesh &mesh, const mfem::Array<int> &marker, bool bdr)
{
  double measure = 0.0;
  PalacePragmaOmp(parallel reduction(+ : measure))
  {
    mfem::IsoparametricTransformation T;
    const int ne = bdr ? mesh.GetNBE() : mesh.GetNE();
    PalacePragmaOmp(for schedule(static))
    for (int i = 0; i < ne; i++)
    {
      if (bdr)
      {
        if (!marker[mesh.GetBdrAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetBdrElementTransformation(i, &T);
      }
      else
      {
        if (!marker[mesh.GetAttribute(i) - 1])
        {
          continue;
        }
        mesh.GetElementTransformation(i, &T);
      }
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(T.GetGeometryType(), T.OrderJ());
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(j);
        T.SetIntPoint(&ip);
        measure += ip.weight * T.Weight();
      }
    }
  }
  Mpi::GlobalSum(1, &measure, mesh.GetComm());
  return measure;
}

}  // namespace

double GetSurfaceArea(const mfem::ParMesh &mesh, const mfem::Array<int> &marker)
{
  return GetMeasure(mesh, marker, true);
}

double GetVolume(const mfem::ParMesh &mesh, const mfem::Array<int> &marker)
{
  return GetMeasure(mesh, marker, false);
}

std::unique_ptr<mfem::ParMesh> DistributeSerialMesh(MPI_Comm comm,
                                                    std::unique_ptr<mfem::Mesh> &smesh)
{
  // Generate METIS partitioning on root and distribute using the MeshPartitioner pipeline,
  // which correctly handles shared entity topology and edge orientations.
  std::unique_ptr<int[]> partitioning;
  if (Mpi::Root(comm))
  {
    partitioning = GetMeshPartitioning(*smesh, Mpi::Size(comm), "", false);
  }
  return DistributeMesh(comm, smesh, partitioning.get());
}

namespace
{

// Fill geometry-resolved topological entity counts from a gathered serial mesh. `counts` is
// reset on entry, so this is safe to call repeatedly with the same struct.
//
// Cells are a raw per-geometry count (interior DOFs are never constrained, so raw == true).
// On a conforming mesh, raw vertex, edge, and face counts are also true counts. On a
// nonconforming mesh, CompleteMeshEntityCounts supplies true vertices and edges
// collectively from distributed H1(1) and ND(1) spaces. True faces are counted here by
// dropping hanging SLAVE faces from the NC face list; this list construction is
// non-collective in MFEM.
void FillMeshEntityCounts(mfem::Mesh &smesh, MeshEntityCounts &counts)
{
  // Reset first so the accumulating members (cells, true_faces) never sum across AMR
  // iterations that reuse the same struct.
  counts = MeshEntityCounts{};

  const int dim = smesh.Dimension();
  counts.dim = dim;

  for (int i = 0; i < smesh.GetNE(); i++)
  {
    counts.cells[smesh.GetElementBaseGeometry(i)]++;
  }

  const bool nonconforming = smesh.Nonconforming();
  if (!nonconforming)
  {
    counts.true_vertices = smesh.GetNV();
    counts.true_edges = (dim >= 2) ? smesh.GetNEdges() : 0;
  }

  if (dim == 3)
  {
    // Per-geometry TRUE-face count (GetNFaces() == 0 in 2D). Conforming: every leaf face is
    // true. Nonconforming: count non-SLAVE leaf faces from the NC face list (== RT(0) true
    // size; see the block comment above). Bucketed by GetFaceGeometry either way.
    const mfem::NCMesh::NCList *face_list =
        nonconforming ? &smesh.ncmesh->GetFaceList() : nullptr;
    for (int f = 0; f < smesh.GetNFaces(); f++)
    {
      if (face_list == nullptr ||
          face_list->GetMeshIdType(f) != mfem::NCMesh::NCList::MeshIdType::SLAVE)
      {
        counts.true_faces[smesh.GetFaceGeometry(f)]++;
      }
    }
  }

  counts.domain_attributes =
      std::vector<int>(smesh.attributes.begin(), smesh.attributes.end());
  counts.boundary_attributes =
      std::vector<int>(smesh.bdr_attributes.begin(), smesh.bdr_attributes.end());
  counts.valid = !nonconforming;
}

// Partition the elements of a nonconforming mesh with PT-Scotch, on the element adjacency
// graph through all conforming and nonconforming faces. MFEM's default partition splits the
// sequence of refinement tree leaves into equal chunks, which fragments the subdomains over
// the AMR iterations. The parts are numbered so that the ranks of a node hold neighboring
// subdomains (see NumberPartsOnNodes). Returns false if PT-Scotch fails or would leave a
// process without elements.
bool GetNonconformingRebalancePartitioning(mfem::ParMesh &mesh,
                                           mfem::Array<int> &partitioning)
{
  MPI_Comm comm = mesh.GetComm();
  const int size = Mpi::Size(comm), rank = Mpi::Rank(comm);
  const int ne = mesh.GetNE(), dim = mesh.SpaceDimension();
  std::vector<SCOTCH_Num> vtxdist(size + 1, 0);
  {
    long long ne_ll = ne;
    std::vector<long long> counts(size);
    MPI_Allgather(&ne_ll, 1, MPI_LONG_LONG, counts.data(), 1, MPI_LONG_LONG, comm);
    for (int i = 0; i < size; i++)
    {
      vtxdist[i + 1] = vtxdist[i] + static_cast<SCOTCH_Num>(counts[i]);
    }
    if (std::ranges::min(counts) == 0)
    {
      return false;  // Leave the empty processes to MFEM's partition
    }
  }

  // The elements are numbered in the order of the input mesh: for a mesh with reordered
  // root elements (see SetRootElementOrigin), at[k] is the local element at position
  // pos[at[k]] = k in the order of their root elements in the input mesh (the leaves of a
  // root element are consecutive and in the same order), so that the partition does not
  // depend on the reordering.
  std::vector<int> pos(ne), at(ne);
  std::iota(at.begin(), at.end(), 0);
  if (const std::vector<int> *origin = GetRootElementOrigin(mesh))
  {
    mfem::CoarseFineTransformations transforms;
    mesh.ncmesh->GetRootTransforms(transforms);
    std::ranges::stable_sort(at, {}, [&](int e)
                             { return (*origin)[transforms.embeddings[e].parent]; });
  }
  for (int k = 0; k < ne; k++)
  {
    pos[at[k]] = k;
  }

  // Calls f(e1, e2, v2) for each face between elements e1 and e2 != e1, where v2 is
  // value(e2) exchanged in a piecewise constant field, and e2 < 0 for a face neighbor.
  mfem::L2_FECollection fec(0, mesh.Dimension());
  mfem::ParFiniteElementSpace fes(&mesh, &fec);
  mfem::ParGridFunction field(&fes);
  mfem::Array<int> vdofs;
  auto ForEachFace = [&](auto &&value, auto &&f)
  {
    double *h_field = field.HostWrite();
    for (int e = 0; e < ne; e++)
    {
      fes.GetElementVDofs(e, vdofs);
      h_field[vdofs[0]] = static_cast<double>(value(e));
    }
    field.ExchangeFaceNbrData();
    const double *nbr = field.FaceNbrData().HostRead();
    for (int i = 0; i < mesh.GetNumFacesWithGhost(); i++)
    {
      const auto info = mesh.GetFaceInformation(i);
      const int e1 = info.element[0].index, e2 = info.element[1].index;
      if (info.IsLocal() && e1 != e2)  // A periodic face can connect an element to itself
      {
        f(e1, e2, static_cast<double>(value(e2)));
      }
      else if (info.IsShared())
      {
        fes.GetFaceNbrElementVDofs(e2, vdofs);
        f(e1, -1, nbr[vdofs[0]]);
      }
    }
  };
  std::vector<std::vector<SCOTCH_Num>> adj(ne);
  ForEachFace([&](int e) { return vtxdist[rank] + pos[e]; },
              [&](int e1, int e2, double g2)
              {
                adj[e1].push_back(static_cast<SCOTCH_Num>(g2));
                if (e2 >= 0)
                {
                  adj[e2].push_back(vtxdist[rank] + pos[e1]);
                }
              });
  std::vector<SCOTCH_Num> deg(ne), adjncy;
  for (int k = 0; k < ne; k++)
  {
    auto &a = adj[at[k]];
    std::ranges::sort(a);
    const auto dup = std::ranges::unique(a);
    a.erase(dup.begin(), dup.end());
    deg[k] = static_cast<SCOTCH_Num>(a.size());
    adjncy.insert(adjncy.end(), a.begin(), a.end());
  }
  const auto nedge = static_cast<SCOTCH_Num>(adjncy.size());

  // PT-Scotch's cost grows with the number of processes, as its band graphs connect all
  // pairs of processes. So only the first process of each group of consecutive ranks
  // partitions, with the rows of its group, a contiguous range of the global numbering:
  // about 100000 elements per process and at least 8 processes. SPEED and SAFETY skip the
  // second bipartition attempt and the diffusion. Whatever its build options, the context
  // runs PT-Scotch on one thread per process, with deterministic communication and a fixed
  // seed.
  const int nproc = static_cast<int>(
      std::clamp<SCOTCH_Num>((vtxdist[size] + 99999) / 100000, std::min(size, 8), size));
  MPI_Comm group_comm, part_comm;
  MPI_Comm_split(comm, static_cast<int>(static_cast<long long>(rank) * nproc / size), rank,
                 &group_comm);
  const int group_rank = Mpi::Rank(group_comm), group_size = Mpi::Size(group_comm);
  const bool leader = (group_rank == 0);
  MPI_Comm_split(comm, leader ? 0 : MPI_UNDEFINED, rank, &part_comm);
  std::vector<int> group_counts(2 * group_size), vcounts(group_size), ecounts(group_size),
      vdispls(group_size + 1, 0), edispls(group_size + 1, 0);
  {
    const int loc[2] = {ne, static_cast<int>(nedge)};
    MPI_Gather(loc, 2, MPI_INT, group_counts.data(), 2, MPI_INT, 0, group_comm);
  }
  for (int i = 0; i < group_size; i++)
  {
    vcounts[i] = group_counts[2 * i];
    ecounts[i] = group_counts[2 * i + 1];
    vdispls[i + 1] = vdispls[i] + vcounts[i];
    edispls[i + 1] = edispls[i] + ecounts[i];
  }
  const SCOTCH_Num gne = leader ? vdispls[group_size] : 0,
                   gnedge = leader ? edispls[group_size] : 0;
  std::vector<SCOTCH_Num> gxadj(gne + 1, 0), gadjncy(gnedge + 1), gpart(gne + 1);
  const MPI_Datatype num_type = mpi::DataType<SCOTCH_Num>();
  MPI_Gatherv(deg.data(), ne, num_type, gxadj.data() + 1, vcounts.data(), vdispls.data(),
              num_type, 0, group_comm);
  MPI_Gatherv(adjncy.data(), static_cast<int>(nedge), num_type, gadjncy.data(),
              ecounts.data(), edispls.data(), num_type, 0, group_comm);
  std::partial_sum(gxadj.begin(), gxadj.end(), gxadj.begin());
  int ok = 1;
  if (leader)
  {
    SCOTCH_Context context;
    SCOTCH_Dgraph graph, context_graph;
    SCOTCH_Strat strat;
    SCOTCH_contextInit(&context);
    SCOTCH_dgraphInit(&graph, part_comm);
    SCOTCH_stratInit(&strat);
    ok = (SCOTCH_contextOptionSetNum(&context, SCOTCH_OPTIONNUMDETERMINISTIC, 1) == 0 &&
          SCOTCH_contextRandomClone(&context) == 0 &&
          SCOTCH_contextThreadSpawn(&context, 1, nullptr) == 0);
    if (ok)
    {
      SCOTCH_contextRandomSeed(&context, 1);
    }
    Mpi::GlobalMin(1, &ok, part_comm);
    ok = ok &&
         (SCOTCH_dgraphBuild(&graph, 0, gne, gne, gxadj.data(), nullptr, nullptr, nullptr,
                             gnedge, gnedge, gadjncy.data(), nullptr, nullptr) == 0);
    ok = ok && (SCOTCH_contextBindDgraph(&context, &graph, &context_graph) == 0);
    ok = ok && (SCOTCH_stratDgraphMapBuild(&strat, SCOTCH_STRATSPEED | SCOTCH_STRATSAFETY,
                                           Mpi::Size(part_comm), size, 0.02) == 0);
    ok = ok && (SCOTCH_dgraphPart(&context_graph, size, &strat, gpart.data()) == 0);
    SCOTCH_stratExit(&strat);
    SCOTCH_dgraphExit(&graph);
    SCOTCH_contextExit(&context);
    MPI_Comm_free(&part_comm);
  }
  std::vector<SCOTCH_Num> part(ne);
  {
    std::vector<SCOTCH_Num> kpart(ne);
    MPI_Scatterv(gpart.data(), vcounts.data(), vdispls.data(), num_type, kpart.data(), ne,
                 num_type, 0, group_comm);
    for (int e = 0; e < ne; e++)
    {
      part[e] = kpart[pos[e]];
    }
  }
  MPI_Comm_free(&group_comm);
  Mpi::GlobalMin(1, &ok, comm);
  if (!ok)
  {
    return false;
  }

  // Centroid, number of elements, and mean global element index of each part.
  std::vector<double> xc(static_cast<std::size_t>(size) * dim, 0.0), cnt(size, 0.0),
      gidx(size, 0.0);
  mfem::Vector center;
  for (int k = 0; k < ne; k++)
  {
    const int e = at[k];
    mesh.GetElementCenter(e, center);  // From the nodes, correct on periodic meshes
    for (int d = 0; d < dim; d++)
    {
      xc[part[e] * dim + d] += center(d);
    }
    cnt[part[e]] += 1.0;
    gidx[part[e]] += static_cast<double>(vtxdist[rank] + k);
  }
  Mpi::GlobalSum(static_cast<int>(xc.size()), xc.data(), comm);
  Mpi::GlobalSum(size, cnt.data(), comm);
  Mpi::GlobalSum(size, gidx.data(), comm);
  if (std::find(cnt.begin(), cnt.end(), 0.0) != cnt.end())
  {
    return false;  // Do not leave a process without elements
  }

  // Gather on the root the faces between the parts, with weight 2 for a local face and 1 on
  // each side of a shared face, and the node of each rank, labeled by its lowest rank.
  std::map<std::pair<int, int>, int> faces;
  ForEachFace([&](int e) { return part[e]; },
              [&](int e1, int e2, double p2)
              {
                const int p = static_cast<int>(part[e1]), q = static_cast<int>(p2);
                if (p != q)
                {
                  faces[{std::min(p, q), std::max(p, q)}] += (e2 >= 0) ? 2 : 1;
                }
              });
  static_assert(sizeof(std::array<int, 3>) == 3 * sizeof(int));
  std::vector<std::array<int, 3>> loc, pairs;
  for (const auto &[pq, w] : faces)
  {
    loc.push_back({pq.first, pq.second, w});
  }
  std::vector<int> counts(size), displs(size + 1, 0);
  const int nloc = 3 * static_cast<int>(loc.size());
  MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  for (int i = 0; i < size; i++)
  {
    displs[i + 1] = displs[i] + counts[i];
  }
  pairs.resize(rank == 0 ? displs[size] / 3 : 0);  // Repeated pairs add up their weights
  MPI_Gatherv(loc.data(), nloc, MPI_INT, pairs.data(), counts.data(), displs.data(),
              MPI_INT, 0, comm);
  std::vector<int> rank_node(size);
  {
    MPI_Comm node_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node_comm);
    int node = rank;
    Mpi::Broadcast(1, &node, 0, node_comm);
    MPI_Comm_free(&node_comm);
    MPI_Gather(&node, 1, MPI_INT, rank_node.data(), 1, MPI_INT, 0, comm);
  }
  std::vector<int> new_part(size);
  if (rank == 0)
  {
    // Number the parts so that the ranks of a node hold neighboring subdomains. Start from
    // three orders of the parts: along a Hilbert curve through their centroids, with
    // per-axis and with isotropic scaling, and by the mean current global index of their
    // elements.
    std::array<double, 3> lo, hi;
    lo.fill(mfem::infinity());
    hi.fill(-mfem::infinity());
    for (int p = 0; p < size; p++)
    {
      for (int d = 0; d < dim; d++)
      {
        xc[p * dim + d] /= cnt[p];
        lo[d] = std::min(lo[d], xc[p * dim + d]);
        hi[d] = std::max(hi[d], xc[p * dim + d]);
      }
      gidx[p] /= cnt[p];
    }
    double len_max = 0.0;
    for (int d = 0; d < dim; d++)
    {
      len_max = std::max(len_max, hi[d] - lo[d]);
    }
    const int b = (dim == 3) ? 21 : 31;
    const double nb = std::ldexp(1.0, b);
    std::vector<std::vector<int>> orders(3, std::vector<int>(size));
    std::vector<std::uint64_t> key(size, 0);
    for (int s = 0; s < 2; s++)
    {
      for (int p = 0; p < size; p++)
      {
        std::array<std::uint32_t, 3> X = {0, 0, 0};
        for (int d = 0; d < dim; d++)
        {
          const double len = (s == 0) ? hi[d] - lo[d] : len_max;
          const double t = (len > 0.0) ? (xc[p * dim + d] - lo[d]) / len : 0.0;
          X[d] = static_cast<std::uint32_t>(std::clamp(std::floor(t * nb), 0.0, nb - 1.0));
        }
        key[p] = HilbertIndex(X, dim, b);
      }
      std::iota(orders[s].begin(), orders[s].end(), 0);
      std::ranges::stable_sort(orders[s], {}, [&key](int p) { return key[p]; });
    }
    std::iota(orders[2].begin(), orders[2].end(), 0);
    std::ranges::stable_sort(orders[2], {}, [&gidx](int p) { return gidx[p]; });
    new_part = NumberPartsOnNodes(orders, pairs, rank_node);
  }
  Mpi::Broadcast(size, new_part.data(), 0, comm);  // Same numbering on every process
  partitioning.SetSize(ne);
  for (int e = 0; e < ne; e++)
  {
    partitioning[e] = new_part[part[e]];
  }
  return true;
}

}  // namespace

void CompleteMeshEntityCounts(mfem::ParMesh &mesh, MeshEntityCounts &counts)
{
  if (!mesh.Nonconforming())
  {
    return;
  }

  const int dim = mesh.Dimension();
  HYPRE_BigInt true_vertices;
  {
    mfem::H1_FECollection h1_fec(1, dim);
    mfem::ParFiniteElementSpace h1_fespace(&mesh, &h1_fec);
    true_vertices = h1_fespace.GlobalTrueVSize();
  }
  HYPRE_BigInt true_edges;
  {
    mfem::ND_FECollection nd_fec(1, dim);
    mfem::ParFiniteElementSpace nd_fespace(&mesh, &nd_fec);
    true_edges = nd_fespace.GlobalTrueVSize();
  }

  if (Mpi::Root(mesh.GetComm()))
  {
    MFEM_VERIFY(counts.dim == dim && !counts.cells.empty(),
                "Missing gathered mesh entity counts!");
    counts.true_vertices = static_cast<long long>(true_vertices);
    counts.true_edges = static_cast<long long>(true_edges);
    counts.valid = true;
  }
}

namespace
{

// Background writing of the adapted mesh on the root process. The serial mesh is written by
// one process while all others would wait, so when the node has a CPU not used by any MPI
// process, the root writes it from a thread pinned to that CPU and continues meanwhile.
class BackgroundMeshWriter
{
private:
  std::thread thread;
  std::shared_ptr<void> payload;
  std::string error;
  int cpu = -2;  // Idle CPU for the writer (-1 if none, -2 if not yet determined)

public:
  ~BackgroundMeshWriter()
  {
    if (thread.joinable())
    {
      thread.join();
    }
  }

  // Determine (collectively, once) whether background writing is possible.
  bool Available(MPI_Comm comm)
  {
    if (cpu == -2)
    {
      cpu = -1;
#if defined(__linux__)
      MPI_Comm node_comm;
      MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                          &node_comm);
      cpu_set_t set;
      CPU_ZERO(&set);
      sched_getaffinity(0, sizeof(set), &set);
      std::vector<unsigned char> used(CPU_SETSIZE, 0), all_used(CPU_SETSIZE, 0);
      for (int c = 0; c < CPU_SETSIZE; c++)
      {
        used[c] = CPU_ISSET(c, &set) ? 1 : 0;
      }
      MPI_Allreduce(used.data(), all_used.data(), CPU_SETSIZE, MPI_UNSIGNED_CHAR, MPI_BOR,
                    node_comm);
      MPI_Comm_free(&node_comm);
      const int ncpu = std::min<int>(CPU_SETSIZE, std::thread::hardware_concurrency());
      for (int c = ncpu - 1; c >= 0; c--)
      {
        if (!all_used[c])
        {
          cpu = c;
          break;
        }
      }
#endif
      // Only the root writes: all processes agree on the root's decision.
      Mpi::Broadcast(1, &cpu, 0, comm);
      if (!Mpi::Root(comm) && cpu >= 0)
      {
        cpu = 0;
      }
    }
    return cpu >= 0;
  }

  // Wait for a pending write (root), and report its errors.
  void Wait()
  {
    if (thread.joinable())
    {
      thread.join();
    }
    payload.reset();
    if (!error.empty())
    {
      const auto msg = error;
      error.clear();
      MFEM_ABORT("Error writing the adapted mesh: " << msg);
    }
  }

  template <typename T, typename Func>
  void Start(std::unique_ptr<T> &&mesh, Func &&write)
  {
    Wait();
    std::shared_ptr<T> smesh(std::move(mesh));
    payload = smesh;
    thread = std::thread(
        [this, smesh, write = std::forward<Func>(write)]()
        {
#if defined(__linux__)
          cpu_set_t set;
          CPU_ZERO(&set);
          CPU_SET(cpu, &set);
          pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#endif
          try
          {
            write(*smesh);
          }
          catch (const std::exception &e)
          {
            error = e.what();
          }
          catch (...)
          {
            error = "Unknown filesystem error";
          }
        });
  }
};

BackgroundMeshWriter &GetBackgroundMeshWriter()
{
  static BackgroundMeshWriter writer;
  return writer;
}

// Whether the adapted mesh is written in the foreground even when a CPU is idle (for
// testing).
bool ForegroundMeshOutput()
{
  return GetMeshOutputTestSettings().foreground;
}

// Concurrent formatting and compression of a large file written by the root process (the
// adapted mesh), while the other processes wait without spinning. The text is a sequence of
// chunks, each formatted into its own stream (with the stream settings of the sequential
// writer), which threads format and compress independently. Compressed, each chunk is raw
// deflate data ending on a byte boundary (Z_SYNC_FLUSH, Z_FINISH for the last one), so the
// chunks form a single gzip member with the combined CRC-32 (as pigz writes it) and the
// file decompresses to exactly the bytes of the sequential writer.
using TextChunk = std::function<void(std::ostream &)>;

// Call f(k) for k in [0, n) from nt threads (the calling one included).
void ParallelFor(std::size_t n, int nt, const std::function<void(std::size_t)> &f)
{
  std::atomic<std::size_t> next = 0;
  auto Work = [&]()
  {
    for (std::size_t k; (k = next.fetch_add(1)) < n;)
    {
      f(k);
    }
  };
  std::vector<std::thread> threads;
  const auto m = static_cast<int>(std::min<std::size_t>(std::max(nt, 1), n));
  threads.reserve(std::max(m - 1, 0));
  for (int t = 1; t < m; t++)
  {
    threads.emplace_back(Work);
    UnpinThread(threads.back());
  }
  Work();
  for (auto &t : threads)
  {
    t.join();
  }
}

// Items per chunk of the concurrent text output.
std::size_t ChunkItems()
{
  const std::size_t n = GetMeshOutputTestSettings().chunk_items;
  return (n > 0) ? n : (std::size_t(1) << 17);
}

// Add the chunks of the lines of items [0, n), block items per chunk.
template <typename Func>
void AddChunkRanges(std::vector<TextChunk> &chunks, std::size_t n, std::size_t block,
                    Func &&print_range)
{
  for (std::size_t b = 0; b < n; b += block)
  {
    const std::size_t e = std::min(n, b + block);
    chunks.emplace_back([print_range, b, e](std::ostream &os) { print_range(os, b, e); });
  }
}

// Add the chunks of mfem::Vector::Print(os, width) for the values (shared by the chunks).
void AddVectorPrintChunks(std::vector<TextChunk> &chunks,
                          std::shared_ptr<const std::vector<double>> values, int width)
{
  const std::size_t block = ChunkItems();
  const std::size_t n = values->size();
  AddChunkRanges(chunks, n, block,
                 [values, width, n](std::ostream &os, std::size_t b, std::size_t e)
                 {
                   const auto &v = *values;
                   for (std::size_t i = b; i < e; i++)
                   {
                     os << mfem::ZeroSubnormal(v[i]);
                     os << ((i + 1 == n || (i + 1) % width == 0) ? '\n' : ' ');
                   }
                 });
}

// Root: write the chunks to the file (with compression, as gzip) using nt threads. Throws
// on filesystem errors.
void WriteTextChunks(const fs::path &path, const std::vector<TextChunk> &chunks,
                     bool compress, int nt)
{
#if !defined(MFEM_USE_ZLIB)
  compress = false;
#endif
  struct Output
  {
    std::string data;
    unsigned long crc = 0;
    std::size_t size = 0;
    bool ready = false;
  };
  const std::size_t n = chunks.size();
  std::vector<Output> out(n);
  std::mutex mutex;
  std::condition_variable cv;
  std::size_t written = 0;  // Chunks written to the file (guarded by mutex)
  std::string error;        // First formatting or compression error (guarded by mutex)
  const std::size_t window = 2 * static_cast<std::size_t>(std::max(nt, 1)) + 2;
  auto Produce = [&](std::size_t k)
  {
    {
      // Limit the chunks held in memory ahead of the file.
      std::unique_lock lock(mutex);
      cv.wait(lock, [&]() { return k < written + window || !error.empty(); });
      if (!error.empty())
      {
        return;
      }
    }
    Output o;
    try
    {
      std::ostringstream os;
      os << std::scientific;
      os.precision(MSH_FLT_PRECISION);
      chunks[k](os);
      std::string text = std::move(os).str();
      o.size = text.size();
#if defined(MFEM_USE_ZLIB)
      if (compress)
      {
        o.crc = crc32(0L, reinterpret_cast<const Bytef *>(text.data()),
                      static_cast<uInt>(text.size()));
        z_stream zs = {};
        if (deflateInit2(&zs, Z_BEST_SPEED, Z_DEFLATED, -MAX_WBITS, 8,
                         Z_DEFAULT_STRATEGY) != Z_OK)
        {
          throw std::runtime_error("deflateInit2 failed");
        }
        o.data.resize(deflateBound(&zs, static_cast<uLong>(text.size())) + 64);
        zs.next_in = reinterpret_cast<Bytef *>(text.data());
        zs.avail_in = static_cast<uInt>(text.size());
        const int flush = (k + 1 == n) ? Z_FINISH : Z_SYNC_FLUSH;
        std::size_t used = 0;
        for (;;)
        {
          zs.next_out = reinterpret_cast<Bytef *>(o.data.data() + used);
          zs.avail_out = static_cast<uInt>(o.data.size() - used);
          const int ret = deflate(&zs, flush);
          used = o.data.size() - zs.avail_out;
          if (ret == Z_STREAM_END || (flush == Z_SYNC_FLUSH && ret == Z_OK &&
                                      zs.avail_in == 0 && zs.avail_out > 0))
          {
            break;
          }
          if (ret != Z_OK && ret != Z_BUF_ERROR)
          {
            deflateEnd(&zs);
            throw std::runtime_error("deflate failed");
          }
          o.data.resize(2 * o.data.size());
        }
        deflateEnd(&zs);
        o.data.resize(used);
      }
      else
#endif
      {
        o.data = std::move(text);
      }
    }
    catch (const std::exception &e)
    {
      std::lock_guard lock(mutex);
      if (error.empty())
      {
        error = e.what();
      }
      cv.notify_all();
      return;
    }
    {
      std::lock_guard lock(mutex);
      out[k] = std::move(o);
      out[k].ready = true;
    }
    cv.notify_all();
  };

  // The calling thread writes the chunks in order while the workers produce them.
  fs::remove_all(path);
  std::thread producer([&]() { ParallelFor(n, nt, Produce); });
  UnpinThread(producer);
  std::string write_error;
  {
    std::ofstream file(path, std::ios_base::binary);
    if (!file.is_open())
    {
      write_error = "Could not open output file";
    }
    unsigned long crc = 0;
    std::size_t total = 0;
#if defined(MFEM_USE_ZLIB)
    if (compress && write_error.empty())
    {
      // Header as written by zlib for Z_BEST_SPEED: no file name or time, OS = Unix.
      const char header[10] = {'\x1f', '\x8b', '\x08', 0, 0, 0, 0, 0, '\x04', '\x03'};
      file.write(header, sizeof(header));
    }
#endif
    for (std::size_t k = 0; k < n; k++)
    {
      Output o;
      {
        std::unique_lock lock(mutex);
        cv.wait(lock, [&]() { return out[k].ready || !error.empty(); });
        if (!error.empty())
        {
          break;
        }
        o = std::move(out[k]);
        out[k].data = std::string();
      }
      if (write_error.empty())
      {
        file.write(o.data.data(), static_cast<std::streamsize>(o.data.size()));
        if (file.fail())
        {
          write_error = "Could not write output file";
        }
      }
#if defined(MFEM_USE_ZLIB)
      if (compress)
      {
        crc = crc32_combine(crc, o.crc, static_cast<z_off_t>(o.size));
      }
#endif
      total += o.size;
      {
        std::lock_guard lock(mutex);
        written = k + 1;
      }
      cv.notify_all();
    }
#if defined(MFEM_USE_ZLIB)
    if (compress && write_error.empty() && error.empty())
    {
      unsigned char trailer[8];
      const auto isize = static_cast<std::uint32_t>(total & 0xFFFFFFFFu);
      for (int i = 0; i < 4; i++)
      {
        trailer[i] = static_cast<unsigned char>((crc >> (8 * i)) & 0xFF);
        trailer[4 + i] = static_cast<unsigned char>((isize >> (8 * i)) & 0xFF);
      }
      file.write(reinterpret_cast<const char *>(trailer), sizeof(trailer));
    }
#endif
    file.close();
    if (write_error.empty() && file.fail())
    {
      write_error = "Could not write output file";
    }
  }
  producer.join();
  if (!error.empty())
  {
    throw std::runtime_error("Error formatting or compressing \"" + path.string() +
                             "\": " + error);
  }
  if (!write_error.empty())
  {
    throw fs::filesystem_error(write_error, path,
                               std::make_error_code(std::errc::io_error));
  }
}

// Write a file from the root process with one thread per process on the root's node, from
// the chunks added by add_chunks(chunks, threads) on the root (collective).
template <typename Func>
void WriteRootOutputChunks(const fs::path &path, MPI_Comm comm, Func &&add_chunks,
                           bool compress)
{
  BlockTimer bt(Timer::IO);
  int nt = 1;
  {
    MPI_Comm node_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
    nt = Mpi::Size(node_comm);
    MPI_Comm_free(&node_comm);
  }
  std::string local_error;
  if (Mpi::Root(comm))
  {
    try
    {
      std::vector<TextChunk> chunks;
      add_chunks(chunks, nt);
      WriteTextChunks(path, chunks, compress, nt);
    }
    catch (const std::exception &e)
    {
      local_error = e.what();
    }
  }
  SleepyBarrier(comm);
  internal::AbortOnNodeFilesystemError(
      comm, fmt::format("Writing output file \"{}\"", path.string()), local_error);
}

// The refinement hierarchy of a nonconforming ParMesh gathered on the root process, to
// write the serial mesh without constructing it: the file is the nonconforming mesh
// followed by the mesh nodes. Supported for meshes of simplices with linear nodes, for
// which the nodes are the vertex coordinates, and the refinement vertices are edge
// midpoints with the coordinates computed as by the refinement of the nodes.
class SerialNCMesh : public mfem::ParNCMesh
{
private:
  // Root: coordinates of the top-level vertices (by node id) and the header of the nodes
  // grid function.
  std::vector<double> top_coords;
  std::string nodes_header;
  int sdim = 0;
  bool by_vdim = true, valid = false;

public:
  // Whether the serial mesh of the given mesh can be written this way (collective).
  static bool Supported(const mfem::ParMesh &mesh)
  {
    int ok = mesh.Nonconforming() && !mesh.NURBSext;
    const auto simplex =
        (mesh.Dimension() == 3) ? mfem::Geometry::TETRAHEDRON : mfem::Geometry::TRIANGLE;
    for (int e = 0; ok && e < mesh.GetNE(); e++)
    {
      ok = (mesh.GetElementBaseGeometry(e) == simplex);
    }
    if (ok && mesh.GetNodes())
    {
      const auto *fes = mesh.GetNodes()->FESpace();
      ok = (dynamic_cast<const mfem::H1_FECollection *>(fes->FEColl()) != nullptr &&
            fes->GetMaxElementOrder() == 1 && fes->GetVDim() == mesh.SpaceDimension());
    }
    Mpi::GlobalMin(1, &ok, mesh.GetComm());
    return ok;
  }

  // Copy the hierarchy and gather all elements on the root (collective). Valid only if
  // supported (with standard refinement scaling).
  SerialNCMesh(const mfem::ParMesh &mesh) : mfem::ParNCMesh(*mesh.pncmesh)
  {
    MPI_Comm comm = mesh.GetComm();
    int ok = !using_scaling;
    Mpi::GlobalMin(1, &ok, comm);
    valid = ok;
    if (!valid)
    {
      return;
    }
    sdim = mesh.SpaceDimension();
    if (const auto *mesh_nodes = mesh.GetNodes())
    {
      // Coordinates of the local top-level vertices, from the mesh nodes.
      const auto *fes = mesh_nodes->FESpace();
      by_vdim = (fes->GetOrdering() == mfem::Ordering::byVDIM);
      if (Mpi::Root(comm))
      {
        std::ostringstream os;
        fes->Save(os);
        nodes_header = os.str();
      }
      std::vector<double> local;
      mfem::Array<int> dofs;
      for (int v = 0; v < mesh.GetNV(); v++)
      {
        const int id = vertex_nodeId[v];
        if (nodes[id].p1 != nodes[id].p2)
        {
          continue;
        }
        fes->GetVertexDofs(v, dofs);
        local.push_back(id);
        for (int d = 0; d < sdim; d++)
        {
          local.push_back((*mesh_nodes)(fes->DofToVDof(dofs[0], d)));
        }
      }
      // Gathered as records of the node id and the coordinates.
      MPI_Datatype record;
      MPI_Type_contiguous(1 + sdim, MPI_DOUBLE, &record);
      MPI_Type_commit(&record);
      int count = static_cast<int>(local.size() / (1 + sdim));
      std::vector<int> counts(Mpi::Size(comm)), displs(Mpi::Size(comm) + 1, 0);
      MPI_Gather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
      std::vector<double> all;
      if (Mpi::Root(comm))
      {
        for (int r = 0; r < Mpi::Size(comm); r++)
        {
          displs[r + 1] = displs[r] + counts[r];
        }
        all.resize((1 + sdim) * std::size_t(displs.back()));
      }
      MPI_Gatherv(local.data(), count, record, all.data(), counts.data(), displs.data(),
                  record, 0, comm);
      MPI_Type_free(&record);
      if (Mpi::Root(comm))
      {
        for (std::size_t i = 0; i < all.size(); i += 1 + sdim)
        {
          const auto id = static_cast<std::size_t>(all[i]);
          if (top_coords.size() < sdim * (id + 1))
          {
            top_coords.resize(sdim * (id + 1), 0.0);
          }
          std::copy(all.begin() + i + 1, all.begin() + i + 1 + sdim,
                    top_coords.begin() + sdim * id);
        }
      }
    }

    // All elements to the root, as by Rebalance with this partition but without its final
    // pruning of the elements beyond the new ghost layers: the root owns every element, so
    // it would remove none of its elements while passing over all of them, and the other
    // processes discard this copy.
    mfem::Array<int> new_ranks(leaf_elements.Size());
    new_ranks = -1;
    for (int i = 0; i < NElements; i++)
    {
      new_ranks[i] = 0;
    }
    RedistributeElements(new_ranks, -1, false);
  }

  bool Valid() const { return valid; }

  // Root: per-geometry element counts of the serial mesh.
  void CountElements(std::map<mfem::Geometry::Type, long long> &cells) const
  {
    for (int i = 0; i < NElements; i++)
    {
      cells[elements[leaf_elements[i]].Geom()]++;
    }
  }

  // Root: add the text chunks of mfem::NCMesh::Print(os) of a nonconforming mesh (v1.0 and
  // v1.1 formats, without comments; the coordinates section only with print_coordinates,
  // as mfem::Mesh::Print for a mesh with nodes). The counts of the boundary section, which
  // precede it, are computed with nt threads.
  static void AddPrintChunks(const mfem::NCMesh &m, bool print_coordinates, int nt,
                             std::vector<TextChunk> &chunks)
  {
    // The chunks refer to the mesh data through pointers.
    const auto *elements = &(m.*(&SerialNCMesh::elements));
    const auto *faces = &(m.*(&SerialNCMesh::faces));
    const auto *nodes = &(m.*(&SerialNCMesh::nodes));
    const auto *root_state = &(m.*(&SerialNCMesh::root_state));
    const auto *coordinates = &(m.*(&SerialNCMesh::coordinates));
    const int dim = m.*(&SerialNCMesh::Dim);
    const int space_dim = m.*(&SerialNCMesh::spaceDim);
    const int my_rank = m.*(&SerialNCMesh::MyRank);
    const bool scaling = m.*(&SerialNCMesh::using_scaling);
    const std::size_t ne = elements->Size();

    // Header and elements.
    chunks.emplace_back(
        [dim, my_rank, scaling, ne](std::ostream &os)
        {
          os << (scaling ? "MFEM NC mesh v1.1\n\n" : "MFEM NC mesh v1.0\n\n");
          os << "# NCMesh supported geometry types:\n"
                "# SEGMENT     = 1\n"
                "# TRIANGLE    = 2\n"
                "# SQUARE      = 3\n"
                "# TETRAHEDRON = 4\n"
                "# CUBE        = 5\n"
                "# PRISM       = 6\n"
                "# PYRAMID     = 7\n";
          os << "\ndimension\n" << dim << "\n";
          os << "\nrank\n" << my_rank << "\n";
          os << "\n# rank attr geom ref_type nodes/children";
          os << "\nelements\n" << ne << "\n";
        });
    const std::size_t block = ChunkItems();
    AddChunkRanges(chunks, ne, block,
                   [elements](std::ostream &os, std::size_t b, std::size_t e)
                   {
                     for (std::size_t i = b; i < e; i++)
                     {
                       const auto &el = (*elements)[static_cast<int>(i)];
                       os << el.rank << " " << el.attribute << " ";
                       if (el.parent == -2)
                       {
                         os << "-1\n";
                         continue;
                       }
                       os << int(el.geom) << " " << int(el.ref_type);
                       for (int j = 0; j < mfem::NCMesh::MaxElemNodes && el.node[j] >= 0;
                            j++)
                       {
                         os << " " << el.node[j];
                       }
                       os << "\n";
                     }
                   });

    // Boundary (faces of leaf elements with a boundary attribute, in element order).
    auto PrintBoundary =
        [elements, faces, dim](std::ostream *os, std::size_t b, std::size_t e)
    {
      static const int nfv2geom[5] = {mfem::Geometry::INVALID, mfem::Geometry::POINT,
                                      mfem::Geometry::SEGMENT, mfem::Geometry::TRIANGLE,
                                      mfem::Geometry::SQUARE};
      const int deg = (dim == 2) ? 2 : 1;
      std::size_t count = 0;
      for (std::size_t i = b; i < e; i++)
      {
        const auto &el = (*elements)[static_cast<int>(i)];
        if (!el.IsLeaf())
        {
          continue;
        }
        const auto &gi = GI[el.Geom()];
        for (int k = 0; k < gi.nf; k++)
        {
          const int *fv = gi.faces[k];
          const int nfv = gi.nfv[k];
          const auto *face =
              faces->Find(el.node[fv[0]], el.node[fv[1]], el.node[fv[2]], el.node[fv[3]]);
          MFEM_ASSERT(face != NULL, "face not found");
          if (face->Boundary())
          {
            count++;
            if (os)
            {
              (*os) << face->attribute << " " << nfv2geom[nfv];
              for (int j = 0; j < nfv; j++)
              {
                (*os) << " " << el.node[fv[j * deg]];
              }
              (*os) << "\n";
            }
          }
        }
      }
      return count;
    };
    const std::size_t nblocks = (ne + block - 1) / block;
    std::vector<std::size_t> counts(nblocks, 0);
    ParallelFor(
        nblocks, nt, [&](std::size_t k)
        { counts[k] = PrintBoundary(nullptr, k * block, std::min(ne, (k + 1) * block)); });
    const std::size_t nb = std::accumulate(counts.begin(), counts.end(), std::size_t(0));
    if (nb)
    {
      chunks.emplace_back(
          [nb](std::ostream &os)
          {
            os << "\n# attr geom nodes";
            os << "\nboundary\n" << nb << "\n";
          });
      AddChunkRanges(chunks, ne, block,
                     [PrintBoundary](std::ostream &os, std::size_t b, std::size_t e)
                     { PrintBoundary(&os, b, e); });
    }

    // Vertex parents (vertex nodes with parents, in node id order).
    const std::size_t nids = nodes->NumIds();
    std::size_t nvp = 0;
    bool uniform_scaling = true;
    for (std::size_t id = 0; id < nids; id++)
    {
      if (!nodes->IdExists(static_cast<int>(id)))
      {
        continue;
      }
      const auto &node = (*nodes)[static_cast<int>(id)];
      if (node.HasVertex() && node.p1 != node.p2)
      {
        nvp++;
        uniform_scaling = uniform_scaling && (node.GetScale() == 0.5);
      }
    }
    MFEM_VERIFY(scaling || uniform_scaling, "NCMesh has nonuniform "
                                            "scaling. Call Mesh::SetScaledNCMesh first.");
    if (nvp)
    {
      chunks.emplace_back(
          [nvp](std::ostream &os)
          {
            os << "\n# vert_id p1 p2";
            os << "\nvertex_parents\n" << nvp << "\n";
          });
      AddChunkRanges(chunks, nids, block,
                     [nodes, scaling](std::ostream &os, std::size_t b, std::size_t e)
                     {
                       for (std::size_t id = b; id < e; id++)
                       {
                         if (!nodes->IdExists(static_cast<int>(id)))
                         {
                           continue;
                         }
                         const auto &node = (*nodes)[static_cast<int>(id)];
                         if (node.HasVertex() && node.p1 != node.p2)
                         {
                           os << id << " " << node.p1 << " " << node.p2;
                           if (scaling)
                           {
                             os << " " << node.GetScale();
                           }
                           os << "\n";
                         }
                       }
                     });
    }

    // Root element orientations.
    if (std::any_of(root_state->begin(), root_state->end(), [](int r) { return r != 0; }))
    {
      chunks.emplace_back(
          [root_state](std::ostream &os)
          {
            os << "\n# root element orientation";
            os << "\nroot_state\n" << root_state->Size() << "\n";
            for (int i = 0; i < root_state->Size(); i++)
            {
              os << (*root_state)[i] << "\n";
            }
          });
    }

    // Top-level vertex coordinates.
    if (print_coordinates && coordinates->Size())
    {
      const std::size_t nv = coordinates->Size() / 3;
      chunks.emplace_back(
          [nv, space_dim](std::ostream &os)
          {
            os << "\n# top-level node coordinates";
            os << "\ncoordinates\n";
            os << nv << "\n";
            if (nv)
            {
              os << space_dim << "\n";
            }
          });
      AddChunkRanges(
          chunks, nv, block,
          [coordinates, space_dim](std::ostream &os, std::size_t b, std::size_t e)
          {
            for (std::size_t i = b; i < e; i++)
            {
              os << (*coordinates)[static_cast<int>(3 * i)];
              for (int j = 1; j < space_dim; j++)
              {
                os << " " << (*coordinates)[static_cast<int>(3 * i + j)];
              }
              os << "\n";
            }
          });
    }
  }

  // Root: add the text chunks of PrintSerial (the chunks refer to this object).
  void AddSerialChunks(double scale, int nt, std::vector<TextChunk> &chunks)
  {
    AddPrintChunks(*this, true, nt, chunks);
    if (!nodes_header.empty())
    {
      auto values = std::make_shared<std::vector<double>>(SerialNodes(scale));
      chunks.emplace_back(
          [header = nodes_header](std::ostream &os)
          {
            os << "\n# mesh curvature GridFunction";
            os << "\nnodes\n";
            os << header << '\n';
          });
      AddVectorPrintChunks(chunks, values, by_vdim ? sdim : 1);
    }
    chunks.emplace_back([](std::ostream &os) { os << "\nmfem_mesh_end" << std::endl; });
  }

  // Root: the nodes (vertex coordinates by Mesh vertex index, scaled), with the refinement
  // vertices at the midpoints of their parent vertices (as computed by the refinement of
  // linear nodes).
  std::vector<double> SerialNodes(double scale) const
  {
    std::vector<double> coords(sdim * nodes.NumIds());
    std::vector<char> done(nodes.NumIds(), 0);
    std::vector<int> stack;
    auto Compute = [&](int id0)
    {
      stack.assign(1, id0);
      while (!stack.empty())
      {
        const int id = stack.back();
        const auto &nd = nodes[id];
        if (nd.p1 == nd.p2)
        {
          MFEM_VERIFY(top_coords.size() >= sdim * (std::size_t(id) + 1),
                      "Missing coordinates of a top-level vertex!");
          std::copy_n(top_coords.begin() + sdim * id, sdim, coords.begin() + sdim * id);
        }
        else if (!done[nd.p1] || !done[nd.p2])
        {
          if (!done[nd.p1])
          {
            stack.push_back(nd.p1);
          }
          if (!done[nd.p2])
          {
            stack.push_back(nd.p2);
          }
          continue;
        }
        else
        {
          for (int d = 0; d < sdim; d++)
          {
            coords[sdim * id + d] =
                0.5 * coords[sdim * nd.p1 + d] + 0.5 * coords[sdim * nd.p2 + d];
          }
        }
        done[id] = 1;
        stack.pop_back();
      }
    };
    std::vector<double> values(static_cast<std::size_t>(sdim) * NVertices, 0.0);
    for (const auto &node : nodes)
    {
      if (!node.HasVertex() || node.vert_index < 0 || node.vert_index >= NVertices)
      {
        continue;
      }
      const int id = nodes.FindId(node.p1, node.p2);
      if (!done[id])
      {
        Compute(id);
      }
      for (int d = 0; d < sdim; d++)
      {
        const int i =
            by_vdim ? sdim * node.vert_index + d : d * NVertices + node.vert_index;
        values[i] = coords[sdim * id + d] * scale;
      }
    }
    return values;
  }

  // Root: write the serial mesh (as mfem::Mesh::Print), with the coordinates scaled.
  void PrintSerial(std::ostream &os, double scale)
  {
    Print(os);
    if (!nodes_header.empty())
    {
      auto values = SerialNodes(scale);
      os << "\n# mesh curvature GridFunction";
      os << "\nnodes\n";
      os << nodes_header << '\n';
      mfem::Vector v(values.data(), static_cast<int>(values.size()));
      v.Print(os, by_vdim ? sdim : 1);
    }
    os << "\nmfem_mesh_end" << std::endl;
  }
};

// Root: add the text chunks of mfem::Mesh::Print for a nonconforming serial mesh (the
// chunks refer to the mesh).
void AddSerialMeshChunks(const mfem::Mesh &smesh, int nt, std::vector<TextChunk> &chunks)
{
  MFEM_VERIFY(smesh.Nonconforming() && !smesh.NURBSext,
              "Chunked output requires a nonconforming mesh!");
  const auto *nodes = smesh.GetNodes();
  SerialNCMesh::AddPrintChunks(*smesh.ncmesh, nodes == nullptr, nt, chunks);
  if (nodes)
  {
    // As mfem::GridFunction::Save (the nodes of H1 spaces have no negative dof signs, so
    // mfem::ParGridFunction::Save writes the same values).
    const auto *fes = nodes->FESpace();
    std::ostringstream header;
    fes->Save(header);
    const double *data = nodes->HostRead();
    auto values = std::make_shared<std::vector<double>>(data, data + nodes->Size());
    chunks.emplace_back(
        [header = header.str()](std::ostream &os)
        {
          os << "\n# mesh curvature GridFunction";
          os << "\nnodes\n";
          os << header << '\n';
        });
    AddVectorPrintChunks(chunks, values,
                         (fes->GetOrdering() == mfem::Ordering::byNODES) ? 1
                                                                         : fes->GetVDim());
  }
  chunks.emplace_back([](std::ostream &os) { os << "\nmfem_mesh_end" << std::endl; });
}

}  // namespace

MeshOutputTestSettings &GetMeshOutputTestSettings()
{
  static MeshOutputTestSettings settings;
  return settings;
}

void FinishMeshOutput()
{
  GetBackgroundMeshWriter().Wait();
}

double RebalanceMesh(const IoData &iodata, std::unique_ptr<mfem::ParMesh> &mesh,
                     MeshEntityCounts *out_counts)
{
  BlockTimer bt0(Timer::REBALANCE);
  MPI_Comm comm = mesh->GetComm();
  if (iodata.model.refinement.save_adapt_mesh)
  {
    BlockTimer bt1(Timer::ADAPT_MESH_OUTPUT);
    // Create a separate serial mesh to write to disk. When MFEM is built with zlib the mesh
    // is gzip-compressed and written with a ".meshgz" extension; MFEM detects gzip from the
    // file's magic bytes on read, so a restart reopens it transparently regardless of the
    // name. Without zlib it is written uncompressed as ".mesh".
    auto sfile = fs::path(iodata.problem.output) / fs::path(iodata.model.mesh).stem();
#if defined(MFEM_USE_ZLIB)
    constexpr bool compress_mesh = true;
    sfile += ".meshgz";
#else
    constexpr bool compress_mesh = false;
    sfile += ".mesh";
#endif

    auto PrintSerial = [&](mfem::Mesh &smesh)
    {
      // The gathered serial mesh holds all elements on the root rank, so the counts are
      // exact global topological counts. Compute before writing (the counts are
      // scale-invariant, so dimensionalization inside the write is irrelevant).
      if (out_counts != nullptr && Mpi::Root(comm))
      {
        FillMeshEntityCounts(smesh, *out_counts);
      }
      WriteRootOutputFile(
          sfile, comm,
          [&](std::ostream &stream)
          {
            stream << std::scientific;
            stream.precision(MSH_FLT_PRECISION);
            mesh::DimensionalizeMesh(smesh, iodata.units.GetMeshLengthRelativeScale());
            smesh.Mesh::Print(
                stream);  // Do not need to nondimensionalize the temporary mesh
          },
          compress_mesh);
    };

    std::unique_ptr<SerialNCMesh> ncmesh;
    if (mesh->Nonconforming() && SerialNCMesh::Supported(*mesh))
    {
      ncmesh = std::make_unique<SerialNCMesh>(*mesh);
      if (!ncmesh->Valid())
      {
        ncmesh.reset();
      }
    }
    if (ncmesh)
    {
      // The true face count is the global number of RT(0) true dofs, the other counts are
      // those of the gathered hierarchy.
      long long true_faces = 0;
      if (out_counts != nullptr && mesh->Dimension() == 3)
      {
        mfem::RT_FECollection rt_fec(0, 3);
        mfem::ParFiniteElementSpace rt_fespace(mesh.get(), &rt_fec);
        true_faces = rt_fespace.GlobalTrueVSize();
      }
      if (out_counts != nullptr && Mpi::Root(comm))
      {
        *out_counts = MeshEntityCounts{};
        out_counts->dim = mesh->Dimension();
        ncmesh->CountElements(out_counts->cells);
        if (mesh->Dimension() == 3)
        {
          out_counts->true_faces[mfem::Geometry::TRIANGLE] = true_faces;
        }
        out_counts->domain_attributes =
            std::vector<int>(mesh->attributes.begin(), mesh->attributes.end());
        out_counts->boundary_attributes =
            std::vector<int>(mesh->bdr_attributes.begin(), mesh->bdr_attributes.end());
      }
      const double scale = iodata.units.GetMeshLengthRelativeScale();
      auto &writer = GetBackgroundMeshWriter();
      if (writer.Available(comm) && !ForegroundMeshOutput())
      {
        if (Mpi::Root(comm))
        {
          writer.Start(std::move(ncmesh),
                       [sfile, scale](SerialNCMesh &smesh)
                       {
                         WriteOutputFile(
                             sfile,
                             [&](std::ostream &stream)
                             {
                               stream << std::scientific;
                               stream.precision(MSH_FLT_PRECISION);
                               smesh.PrintSerial(stream, scale);
                             },
                             compress_mesh);
                       });
        }
      }
      else
      {
        WriteRootOutputChunks(
            sfile, comm, [&](std::vector<TextChunk> &chunks, int nt)
            { ncmesh->AddSerialChunks(scale, nt, chunks); }, compress_mesh);
      }
    }
    else if (mesh->Nonconforming())
    {
      auto smesh = std::make_unique<mfem::ParMesh>(*mesh);
      mfem::Array<int> serial_partition(mesh->GetNE());
      serial_partition = 0;
      smesh->Rebalance(serial_partition);
      auto &writer = GetBackgroundMeshWriter();
      if (writer.Available(comm) && !ForegroundMeshOutput())
      {
        // The root writes in the background (no further collective communication for the
        // file, errors are reported when the write is completed).
        if (Mpi::Root(comm))
        {
          if (out_counts != nullptr)
          {
            FillMeshEntityCounts(*smesh, *out_counts);
          }
          const double scale = iodata.units.GetMeshLengthRelativeScale();
          writer.Start(std::move(smesh),
                       [sfile, scale](mfem::ParMesh &smesh)
                       {
                         WriteOutputFile(
                             sfile,
                             [&](std::ostream &stream)
                             {
                               stream << std::scientific;
                               stream.precision(MSH_FLT_PRECISION);
                               mesh::DimensionalizeMesh(smesh, scale);
                               smesh.Mesh::Print(stream);
                             },
                             compress_mesh);
                       });
        }
      }
      else
      {
        if (out_counts != nullptr && Mpi::Root(comm))
        {
          FillMeshEntityCounts(*smesh, *out_counts);
        }
        const double scale = iodata.units.GetMeshLengthRelativeScale();
        WriteRootOutputChunks(
            sfile, comm,
            [&](std::vector<TextChunk> &chunks, int nt)
            {
              mesh::DimensionalizeMesh(*smesh, scale);
              AddSerialMeshChunks(*smesh, nt, chunks);
            },
            compress_mesh);
      }
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
      mfem::Array<int> partitioning;
      if (GetNonconformingRebalancePartitioning(*mesh, partitioning))
      {
        mesh->Rebalance(partitioning);
      }
      else
      {
        Mpi::Warning(comm, "Could not partition the nonconforming mesh with PT-Scotch, "
                           "rebalancing it with MFEM's default partition!\n");
        mesh->Rebalance();
      }
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
                                     const config::BoundaryData &boundaries,
                                     bool defer_nodes, bool *finalize_deferred)
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
    mesh = ReadMeshStream(fi, generate_edges, refine, fix_orientation, finalize_deferred);
  }
  else
  {
    // Otherwise, read a Gmsh mesh with the fast reader where it supports the file (unless
    // the environment variable PALACE_GMSH_FAST is 0), or rely on MFEM to load the mesh.
    // Use named_ifgzstream (as MFEM's own Mesh::LoadFromFile does) so a gzip-compressed
    // mesh, such as a ".meshgz" saved-adapt mesh, is transparently decompressed when MFEM
    // is built with zlib; it detects gzip from the file's magic bytes, so a plain-text mesh
    // reads through unchanged.
    const char *fast_env = std::getenv("PALACE_GMSH_FAST");
    if (!fast_env || std::string(fast_env) != "0")
    {
      mesh = mesh::ReadGmshMeshFast(mesh_file, !finalize_deferred);
    }
    if (mesh)
    {
      mesh = FinishReadMesh(std::move(mesh), refine, fix_orientation, finalize_deferred);
    }
    else
    {
      mfem::named_ifgzstream fi(mesh_file);
      if (!fi.good())
      {
        MFEM_ABORT("Unable to open mesh file \"" << mesh_file << "\"!");
      }
      mesh = ReadMeshStream(fi, generate_edges, refine, fix_orientation, finalize_deferred);
    }
  }
  if (remove_curvature)
  {
    mesh->SetCurvature(-1);
  }
  else if (!defer_nodes || mesh->GetNodes() || !boundaries.periodic.boundary_pairs.empty())
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

// Replace the mesh by a new mesh with the domain and boundary elements which are not
// deleted, where the maps give the index in the new mesh of each element, or -1 if deleted.
void DeleteMeshElements(std::unique_ptr<mfem::Mesh> &orig_mesh,
                        const std::vector<int> &elem_delete_map, int new_ne,
                        const std::vector<int> &bdr_elem_delete_map, int new_nbe)
{
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

// Whether the faces and edges of a mesh have been generated (FinalizeTopology).
bool HasTopology(const mfem::Mesh &mesh)
{
  return mesh.Dimension() < 2 || mesh.GetNEdges() > 0 || mesh.GetNE() == 0;
}

// For each boundary element of a mesh without topology, whether it is a face (an edge in
// 2D) of an element of the mesh, with the face identification of MFEM's topology tables
// (the vertex set of a triangle or an edge, the three smallest vertices of a
// quadrilateral).
std::vector<bool> BdrElementIsFace(const mfem::Mesh &mesh)
{
  const int dim = mesh.Dimension(), nv = mesh.GetNV();
  using Key = std::array<int, 3>;
  auto MakeKey = [dim](const int *v, int n)
  {
    Key k = {v[0], v[1], (dim == 3) ? v[2] : -1};
    if (dim == 3 && n == 4)
    {
      // STable3D::Push4/operator(): drop the largest vertex.
      int imax = 0;
      for (int i = 1; i < 4; i++)
      {
        if (v[i] > v[imax])
        {
          imax = i;
        }
      }
      for (int i = 0, j = 0; i < 4; i++)
      {
        if (i != imax)
        {
          k[j++] = v[i];
        }
      }
    }
    std::sort(k.begin(), k.begin() + ((dim == 3) ? 3 : 2));
    return k;
  };

  // Vertex to element adjacency.
  std::vector<int> offsets(nv + 1, 0), elems;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const mfem::Element *el = mesh.GetElement(e);
    for (int i = 0; i < el->GetNVertices(); i++)
    {
      offsets[el->GetVertices()[i] + 1]++;
    }
  }
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  elems.resize(offsets[nv]);
  {
    std::vector<int> pos(offsets.begin(), offsets.end() - 1);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const mfem::Element *el = mesh.GetElement(e);
      for (int i = 0; i < el->GetNVertices(); i++)
      {
        elems[pos[el->GetVertices()[i]]++] = e;
      }
    }
  }

  std::vector<bool> is_face(mesh.GetNBE(), false);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const mfem::Element *bel = mesh.GetBdrElement(be);
    const Key key = MakeKey(bel->GetVertices(), bel->GetNVertices());
    const int v = key[0];
    for (int k = offsets[v]; k < offsets[v + 1] && !is_face[be]; k++)
    {
      const mfem::Element *el = mesh.GetElement(elems[k]);
      const int *ev = el->GetVertices();
      const int nf = (dim == 3) ? el->GetNFaces() : el->GetNEdges();
      for (int f = 0; f < nf; f++)
      {
        int fv[4];
        int n;
        if (dim == 3)
        {
          n = el->GetNFaceVertices(f);
          for (int i = 0; i < n; i++)
          {
            fv[i] = ev[el->GetFaceVertices(f)[i]];
          }
        }
        else
        {
          n = 2;
          fv[0] = ev[el->GetEdgeVertices(f)[0]];
          fv[1] = ev[el->GetEdgeVertices(f)[1]];
        }
        if (MakeKey(fv, n) == key)
        {
          is_face[be] = true;
          break;
        }
      }
    }
  }
  return is_face;
}

std::unique_ptr<mfem::Mesh> ReadMeshStream(std::istream &input, bool generate_edges,
                                           bool refine, bool fix_orientation,
                                           bool *finalize_deferred)
{
  // Read the mesh without finalizing it (see mfem::Mesh::Load), to first remove the
  // boundary elements which are not a face (an edge in 2D) of any element of the mesh:
  // finalizing the mesh, like most operations on it, accesses memory out of bounds for
  // them. No domain element is attached to these boundary elements, which can come for
  // example from a mesher tagging all of an outer boundary when part of it is also the
  // boundary of a region which is not meshed, like a conductor modeled as a hole in the
  // mesh. In 3D, mfem::Mesh::Loader itself aborts on such boundary elements (MFEM 4.10).
  struct UnfinalizedMesh : public mfem::Mesh
  {
    UnfinalizedMesh(std::istream &input, bool generate_edges)
    {
      Loader(input, generate_edges);
    }
  };
  auto mesh = std::make_unique<mfem::Mesh>(UnfinalizedMesh(input, generate_edges));
  return FinishReadMesh(std::move(mesh), refine, fix_orientation, finalize_deferred);
}

// Remove the boundary elements which are not a face of any element (see ReadMeshStream)
// and finalize a mesh read without finalization.
std::unique_ptr<mfem::Mesh> FinishReadMesh(std::unique_ptr<mfem::Mesh> mesh, bool refine,
                                           bool fix_orientation, bool *finalize_deferred)
{
  // The boundary elements of nonconforming and NURBS meshes are faces of the mesh by
  // construction, and those of 1D meshes are vertices.
  if (mesh->Conforming() && !mesh->NURBSext && mesh->Dimension() > 1)
  {
    // A mesh read without its topology is checked with the element faces directly. In 3D,
    // MFEM aborts on a boundary element which is not a face when generating the topology,
    // which is then done for the same error as when reading the mesh with MFEM.
    std::vector<bool> is_face;
    if (!HasTopology(*mesh))
    {
      is_face = BdrElementIsFace(*mesh);
      if (mesh->Dimension() == 3 && std::ranges::find(is_face, false) != is_face.end())
      {
        mesh->FinalizeTopology();
      }
    }
    auto IsFace = [&](int be)
    { return HasTopology(*mesh) ? mesh->GetBdrElementFaceIndex(be) >= 0 : is_face[be]; };
    std::vector<int> elem_delete_map(mesh->GetNE()),
        bdr_elem_delete_map(mesh->GetNBE(), -1);
    std::iota(elem_delete_map.begin(), elem_delete_map.end(), 0);
    int new_nbe = 0;
    std::set<int> bdr_warn_list;
    std::vector<std::string> warn_vertices;
    for (int be = 0; be < mesh->GetNBE(); be++)
    {
      if (IsFace(be))
      {
        bdr_elem_delete_map[be] = new_nbe++;
        continue;
      }
      if (bdr_warn_list.empty())
      {
        mfem::Array<int> verts;
        mesh->GetBdrElementVertices(be, verts);
        for (auto v : verts)
        {
          const auto *coord = mesh->GetVertex(v);
          warn_vertices.push_back(fmt::format(
              "({:.6g})", fmt::join(coord, coord + mesh->SpaceDimension(), ", ")));
        }
      }
      bdr_warn_list.insert(mesh->GetBdrAttribute(be));
    }
    if (new_nbe < mesh->GetNBE())
    {
      Mpi::Warning(
          "Removed {:d} boundary elements which are not faces of any element of the "
          "mesh, the first one with vertices {}!\nNo domain element is attached to "
          "them, so the boundary conditions on them have no effect.\n",
          mesh->GetNBE() - new_nbe, fmt::join(warn_vertices, ", "));
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
      DeleteMeshElements(mesh, elem_delete_map, mesh->GetNE(), bdr_elem_delete_map,
                         new_nbe);
    }
  }
  // A caller can finalize a mesh without nodes later along with its own finalization
  // (refinement marking), in one topology update instead of two.
  if (finalize_deferred && !mesh->GetNodes())
  {
    *finalize_deferred = true;
    return mesh;
  }
  mesh->Finalize(refine, fix_orientation);
  return mesh;
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
  // the domain. Every boundary element is attached to an element of the mesh (see
  // ReadMeshStream), so this is needed only when elements are deleted.
  if (new_ne == orig_mesh->GetNE())
  {
    return;
  }
  if (!HasTopology(*orig_mesh))
  {
    orig_mesh->FinalizeTopology();
  }
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
  DeleteMeshElements(orig_mesh, elem_delete_map, new_ne, bdr_elem_delete_map, new_nbe);
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

}  // namespace

namespace mesh
{

// Conformally bisect a set of mesh edges by inserting a midpoint on each and splitting
// every tetrahedron and boundary element incident on a split edge into two. Unlike
// mfem::Mesh::GeneralRefinement, this performs a purely local edge-fan split with no
// marked-edge/longest-edge selection and no global conformity-closure propagation. This is
// essential for periodic meshes: MFEM's conforming bisection keys the new midpoint on the
// collapsed vertex pair and cannot keep the two periodic copies of a seam edge consistent,
// which corrupts the seam (aborting later in STable3D). A local edge-fan split touches only
// the ring of elements around each edge, leaving the rest of the mesh (and the periodic
// identification) untouched. The inserted midpoints are later duplicated by the regular
// interior-boundary cracking pass, which decouples the two sides of the crack.
int LocalEdgeSplit(std::unique_ptr<mfem::Mesh> &orig_mesh,
                   const std::vector<std::pair<int, int>> &edges)
{
  // Only pure tetrahedral meshes are supported.
  MFEM_VERIFY(orig_mesh->Dimension() == 3, "LocalEdgeSplit only supports 3D meshes!");
  {
    const auto element_types = mesh::CheckElements(*orig_mesh);
    MFEM_VERIFY(element_types.has_simplices && !element_types.has_hexahedra &&
                    !element_types.has_prisms && !element_types.has_pyramids,
                "LocalEdgeSplit only supports pure tetrahedral meshes!");
  }
  if (edges.empty())
  {
    return 0;
  }

  std::unique_ptr<mfem::Table> vert_to_elem(orig_mesh->GetVertexToElementTable());
  auto edge_ring = [&](int v0, int v1)
  {
    const int *e0 = vert_to_elem->GetRow(v0);
    std::unordered_set<int> at_v0(e0, e0 + vert_to_elem->RowSize(v0));
    std::vector<int> ring;
    const int *e1 = vert_to_elem->GetRow(v1);
    for (int i = 0; i < vert_to_elem->RowSize(v1); i++)
    {
      if (at_v0.find(e1[i]) != at_v0.end())
      {
        ring.push_back(e1[i]);
      }
    }
    return ring;
  };

  // Select a maximal independent set of candidate edges: greedily accept an edge only if
  // none of the tetrahedra in its ring already belong to an accepted edge. Each accepted
  // edge is assigned a unique midpoint vertex id (appended after the original vertices).
  // Map from a sorted vertex pair to the midpoint vertex id. The ids are assigned in a
  // second pass over the (sorted) map so that they increase in the same order the midpoint
  // vertices are subsequently added via AddVertex below; assigning during the (input-order)
  // selection loop would mismatch the id with the AddVertex order whenever the input edge
  // order differs from the sorted-key order, corrupting midpoint coordinates.
  std::map<std::pair<int, int>, int> edge_midpoint;
  std::unordered_set<int> claimed_elem;
  for (const auto &e : edges)
  {
    int v0 = e.first, v1 = e.second;
    auto ring = edge_ring(v0, v1);
    MFEM_VERIFY(!ring.empty(), "Edge to split has no incident elements!");
    bool conflict = false;
    for (int el : ring)
    {
      if (claimed_elem.contains(el))
      {
        conflict = true;
        break;
      }
    }
    if (conflict)
    {
      continue;  // Defer to a subsequent pass
    }
    for (int el : ring)
    {
      claimed_elem.insert(el);
    }
    edge_midpoint[{std::min(v0, v1), std::max(v0, v1)}] = -1;  // id assigned below
  }
  const int num_split = static_cast<int>(edge_midpoint.size());
  MFEM_VERIFY(num_split > 0, "Failed to select any edge to split!");
  {
    int v_mid = orig_mesh->GetNV();
    for (auto &[edge, mid] : edge_midpoint)
    {
      mid = v_mid++;
    }
  }

  // For each element / boundary element, find the (at most one, by independence) split edge
  // it contains, returning the matching midpoint id and the local endpoint vertices, or -1.
  auto find_split_edge = [&edge_midpoint](const int *v, int nv, const int(*edge_vert)[2],
                                          int nedge, int &lv0, int &lv1)
  {
    for (int le = 0; le < nedge; le++)
    {
      int a = v[edge_vert[le][0]], b = v[edge_vert[le][1]];
      auto it = edge_midpoint.find({std::min(a, b), std::max(a, b)});
      if (it != edge_midpoint.end())
      {
        lv0 = edge_vert[le][0];
        lv1 = edge_vert[le][1];
        return it->second;
      }
    }
    return -1;
  };
  static constexpr int tet_edge_vert[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                              {0, 3}, {1, 3}, {2, 3}};
  static constexpr int tri_edge_vert[3][2] = {{0, 1}, {1, 2}, {2, 0}};

  // Count split elements / boundary elements to size the new mesh.
  int n_split_elem = 0, n_split_bdr = 0;
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    int lv0, lv1;
    if (find_split_edge(orig_mesh->GetElement(e)->GetVertices(), 4, tet_edge_vert, 6, lv0,
                        lv1) >= 0)
    {
      n_split_elem++;
    }
  }
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    int lv0, lv1;
    if (find_split_edge(orig_mesh->GetBdrElement(be)->GetVertices(), 3, tri_edge_vert, 3,
                        lv0, lv1) >= 0)
    {
      n_split_bdr++;
    }
  }

  const int new_nv = orig_mesh->GetNV() + num_split;
  const int new_ne = orig_mesh->GetNE() + n_split_elem;
  const int new_nbe = orig_mesh->GetNBE() + n_split_bdr;
  auto new_mesh = std::make_unique<mfem::Mesh>(orig_mesh->Dimension(), new_nv, new_ne,
                                               new_nbe, orig_mesh->SpaceDimension());
  for (int v = 0; v < orig_mesh->GetNV(); v++)
  {
    new_mesh->AddVertex(orig_mesh->GetVertex(v));
  }
  // Midpoint coordinates (arithmetic average; curved geometry is fixed up below for
  // high-order meshes via per-parent node synthesis). The map iteration order (sorted keys)
  // matches the id assignment above, so the vertex AddVertex returns equals the stored id.
  for (const auto &[edge, mid] : edge_midpoint)
  {
    double coord[3] = {0.0, 0.0, 0.0};
    const double *c0 = orig_mesh->GetVertex(edge.first);
    const double *c1 = orig_mesh->GetVertex(edge.second);
    for (int d = 0; d < orig_mesh->SpaceDimension(); d++)
    {
      coord[d] = 0.5 * (c0[d] + c1[d]);
    }
    int added = new_mesh->AddVertex(coord);
    MFEM_VERIFY(added == mid, "Midpoint vertex id mismatch: AddVertex returned "
                                  << added << " but edge midpoint id is " << mid
                                  << "; vertex add order must match id assignment order!");
  }

  // Add elements: a tet containing a split edge bisects into two children (one per endpoint
  // replaced by the midpoint); all others are copied unchanged.
  for (int e = 0; e < orig_mesh->GetNE(); e++)
  {
    const int *v = orig_mesh->GetElement(e)->GetVertices();
    int lv0, lv1;
    int mid = find_split_edge(v, 4, tet_edge_vert, 6, lv0, lv1);
    if (mid < 0)
    {
      new_mesh->AddElement(orig_mesh->GetElement(e)->Duplicate(new_mesh.get()));
      continue;
    }
    int attr = orig_mesh->GetAttribute(e);
    for (int child = 0; child < 2; child++)
    {
      int cv[4];
      for (int i = 0; i < 4; i++)
      {
        cv[i] = (i == ((child == 0) ? lv0 : lv1)) ? mid : v[i];
      }
      new_mesh->AddTet(cv, attr);
    }
  }
  for (int be = 0; be < orig_mesh->GetNBE(); be++)
  {
    const mfem::Element *bel = orig_mesh->GetBdrElement(be);
    const int *bv = bel->GetVertices();
    int lv0, lv1;
    int mid = find_split_edge(bv, 3, tri_edge_vert, 3, lv0, lv1);
    if (mid < 0)
    {
      new_mesh->AddBdrElement(bel->Duplicate(new_mesh.get()));
      continue;
    }
    int attr = bel->GetAttribute();
    for (int child = 0; child < 2; child++)
    {
      int cv[3];
      for (int i = 0; i < 3; i++)
      {
        cv[i] = (i == ((child == 0) ? lv0 : lv1)) ? mid : bv[i];
      }
      new_mesh->AddBdrTriangle(cv, attr);
    }
  }

  constexpr bool generate_bdr = false;
  new_mesh->FinalizeTopology(generate_bdr);

  // Synthesize high-order node coordinates for the split (child) elements, accounting for
  // the change in element count. Unsplit elements transfer their nodes 1:1; each child of a
  // split tet evaluates the parent's nodal field at the child's nodes mapped into the
  // parent reference frame (same idiom as MeshTetToHex). For linear meshes (no Nodes) the
  // geometric midpoint above is exact and nothing more is needed.
  if (orig_mesh->GetNodes())
  {
    new_mesh->EnsureNodes();
    auto *orig_fespace = orig_mesh->GetNodes()->FESpace();
    new_mesh->SetCurvature(orig_fespace->GetMaxElementOrder(), orig_fespace->IsDGSpace(),
                           orig_mesh->SpaceDimension(), orig_fespace->GetOrdering());
    const int sdim = orig_mesh->SpaceDimension();
    const auto *orig_FE = orig_fespace->GetTypicalFE();
    const auto *child_FE = new_mesh->GetNodes()->FESpace()->GetTypicalFE();

    // Reference-space vertex locations of a tetrahedron.
    static constexpr double ref_vert[4][3] = {
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    mfem::Array<int> orig_dofs, child_dofs;
    mfem::Vector dof_vals(orig_FE->GetDof() * sdim);
    mfem::DenseMatrix dof_vals_mat(dof_vals.GetData(), orig_FE->GetDof(), sdim);
    mfem::DenseMatrix shape(orig_FE->GetDof(), child_FE->GetDof());
    mfem::DenseMatrix point_matrix(child_FE->GetDof(), sdim);
    mfem::Vector col;

    int new_e = 0;
    for (int e = 0; e < orig_mesh->GetNE(); e++)
    {
      const int *v = orig_mesh->GetElement(e)->GetVertices();
      int lv0, lv1;
      int mid = find_split_edge(v, 4, tet_edge_vert, 6, lv0, lv1);
      if (mid < 0)
      {
        // Unsplit: copy nodes directly.
        orig_fespace->GetElementVDofs(e, orig_dofs);
        orig_mesh->GetNodes()->GetSubVector(orig_dofs, dof_vals);
        new_mesh->GetNodes()->FESpace()->GetElementVDofs(new_e, child_dofs);
        new_mesh->GetNodes()->SetSubVector(child_dofs, dof_vals);
        new_e++;
        continue;
      }
      // Split tet: the child replaces one endpoint of the split edge with the edge midpoint
      // (reference midpoint of the two endpoint reference coordinates).
      double ref_mid[3];
      for (int d = 0; d < 3; d++)
      {
        ref_mid[d] = 0.5 * (ref_vert[lv0][d] + ref_vert[lv1][d]);
      }
      orig_mesh->GetNodes()->GetElementDofValues(e, dof_vals);
      for (int child = 0; child < 2; child++)
      {
        double child_ref[4][3];
        for (int i = 0; i < 4; i++)
        {
          for (int d = 0; d < 3; d++)
          {
            child_ref[i][d] = ref_vert[i][d];
          }
        }
        int replace = (child == 0) ? lv0 : lv1;
        for (int d = 0; d < 3; d++)
        {
          child_ref[replace][d] = ref_mid[d];
        }
        for (int j = 0; j < child_FE->GetNodes().Size(); j++)
        {
          const auto &cn = child_FE->GetNodes()[j];
          mfem::IntegrationPoint ip;
          double L0 = 1.0 - cn.x - cn.y - cn.z, L1 = cn.x, L2 = cn.y, L3 = cn.z;
          double p[3];
          for (int d = 0; d < 3; d++)
          {
            p[d] = L0 * child_ref[0][d] + L1 * child_ref[1][d] + L2 * child_ref[2][d] +
                   L3 * child_ref[3][d];
          }
          ip.x = p[0];
          ip.y = p[1];
          ip.z = p[2];
          shape.GetColumnReference(j, col);
          orig_FE->CalcShape(ip, col);
        }
        MultAtB(shape, dof_vals_mat, point_matrix);
        new_mesh->GetNodes()->FESpace()->GetElementVDofs(new_e, child_dofs);
        new_mesh->GetNodes()->SetSubVector(child_dofs, point_matrix.GetData());
        new_e++;
      }
    }
  }

  orig_mesh = std::move(new_mesh);
  return num_split;
}

}  // namespace mesh

namespace
{

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
        auto donor = std::ranges::find(da, attr) != da.end();
        auto receiver = std::ranges::find(ra, attr) != ra.end();
        if (donor || receiver)
        {
          mesh.GetFaceElements(f, &e1, &e2);
          MFEM_VERIFY(e1 >= 0 && e2 >= 0,
                      "Mesh is not periodic on attribute " << attr << "!");
        }
      }
    }
    MFEM_VERIFY((e1 >= 0 && e2 >= 0) || !face_to_be.contains(f),
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
    // The rows of v_to_v and the elements are processed independently by threads (the
    // marking of an element changes only that element).
    mfem::Array<mfem::real_t> lengths;
    if (Nodes)
    {
      GetEdgeLengths2(v_to_v, lengths);
    }
    else
    {
      EdgeLengths2(v_to_v, vertices, spaceDim, NumOfEdges, serial_mesh_threads, lengths);
    }
    const auto min_length = 0.01 * lengths.Min();
    constexpr int block = 4096;
    const int nt = serial_mesh_threads;
    const int nrows = v_to_v.NumberOfRows();
    ParallelForThreads((nrows + block - 1) / block, nt,
                       [&](std::size_t b, int)
                       {
                         const int end = std::min(nrows, static_cast<int>(b + 1) * block);
                         for (int i = static_cast<int>(b) * block; i < end; i++)
                         {
                           for (mfem::DSTable::RowIterator it(v_to_v, i); !it; ++it)
                           {
                             int j = it.Column();
                             if (!refinement_edges.contains({i, j}))
                             {
                               // "Zero" the edge lengths which do not connect vertices on
                               // the interface. Avoid zero-length edges just in case.
                               lengths[it.Index()] = min_length;
                             }
                           }
                         }
                       });

    // Finish marking (see mfem::Mesh::MarkTetMeshForRefinement).
    mfem::Array<int> indices(NumOfEdges);
    std::iota(indices.begin(), indices.end(), 0);
    ParallelForThreads((NumOfElements + block - 1) / block, nt,
                       [&](std::size_t b, int)
                       {
                         const int end =
                             std::min(NumOfElements, static_cast<int>(b + 1) * block);
                         for (int i = static_cast<int>(b) * block; i < end; i++)
                         {
                           if (elements[i]->GetType() == mfem::Element::TETRAHEDRON)
                           {
                             MFEM_ASSERT(dynamic_cast<mfem::Tetrahedron *>(elements[i]),
                                         "Unexpected non-Tetrahedron element type!");
                             static_cast<mfem::Tetrahedron *>(elements[i])
                                 ->MarkEdge(v_to_v, lengths, indices);
                           }
                         }
                       });
    ParallelForThreads((NumOfBdrElements + block - 1) / block, nt,
                       [&](std::size_t b, int)
                       {
                         const int end =
                             std::min(NumOfBdrElements, static_cast<int>(b + 1) * block);
                         for (int i = static_cast<int>(b) * block; i < end; i++)
                         {
                           if (boundary[i]->GetType() == mfem::Element::TRIANGLE)
                           {
                             MFEM_ASSERT(dynamic_cast<mfem::Triangle *>(boundary[i]),
                                         "Unexpected non-Triangle element type!");
                             static_cast<mfem::Triangle *>(boundary[i])
                                 ->MarkEdge(v_to_v, lengths, indices);
                           }
                         }
                       });
  }

public:
  EdgeRefinementMesh(mfem::Mesh &&mesh, const T &refinement_edges)
    : mfem::Mesh(std::move(mesh)), refinement_edges(refinement_edges)
  {
  }

  // Mark the elements for refinement as Finalize(true, false) does, for a mesh without
  // nodes: the orientation check and topology update of Finalize are not needed before a
  // local refinement, which updates the topology itself (and the marking preserves the
  // element orientations).
  void Mark()
  {
    MFEM_VERIFY(!GetNodes(), "Marking without finalization requires a mesh without nodes!");
    MarkForRefinement();
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
  // By default, interior boundaries with boundary conditions applying to each side
  // separately are cracked, as well as impedance and conductivity sheets modeled as two
  // independent surfaces ("Crack"). Lumped ports are never cracked, which would give
  // invalid results.
  const auto crack_boundary_attributes =
      iodata.boundaries.GetMeshCrackAttributes(iodata.model.crack_bdr_elements);

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
  // The element sets live in a pool (they are many and small); the containers behave like
  // the std:: ones, including the iteration order which determines the duplicate numbering.
  std::vector<std::unique_ptr<std::pmr::unsynchronized_pool_resource>> component_pools;
  std::unordered_map<int, std::vector<std::pair<int, std::pmr::unordered_set<int>>>>
      crack_vert_duplicates;
  std::unique_ptr<mfem::Table> vert_to_elem;
  if (!crack_boundary_attributes.empty() && iodata.model.crack_bdr_elements.value_or(true))
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
    // Faces whose boundary element is on the crack, for the connectivity search below.
    std::vector<bool> crack_face(orig_mesh->GetNumFaces(), false);
    for (const auto &[f, be] : face_to_be)
    {
      crack_face[f] = crack_bdr_elem.contains(be);
    }
    // Neighbor of each element across each of its faces (the other element of the face, or
    // -1 for a boundary face or a face on the crack), in the order of the element faces:
    // the connectivity search reads these instead of the face table and face information.
    const int *etf_I = elem_to_face.GetI(), *etf_J = elem_to_face.GetJ();
    std::vector<int> face_nbr(etf_I[orig_mesh->GetNE()]);
    for (int e = 0; e < orig_mesh->GetNE(); e++)
    {
      for (int k = etf_I[e]; k < etf_I[e + 1]; k++)
      {
        const int f = etf_J[k];
        if (crack_face[f])
        {
          // Skip element-element connectivities which cross the crack.
          face_nbr[k] = -1;
          continue;
        }
        int e1, e2;
        orig_mesh->GetFaceElements(f, &e1, &e2);
        MFEM_VERIFY(e == e1 || e == e2,
                    "Unexpected face-element connectivity in internal boundary cracking!");
        face_nbr[k] = (e == e1) ? e2 : e1;
      }
    }
    // The crack vertices, each once, in the order of the boundary elements, and the
    // connected components of the elements around each of them: a search on the graph of
    // the elements connected to the vertex, with the adjacencies given by the face
    // connectivity off the crack. The searches are independent and run on threads, with
    // their results kept in the order of the vertices.
    std::vector<int> crack_verts;
    {
      std::unordered_set<int> seen;
      for (auto be : crack_bdr_elem)
      {
        const mfem::Element *bdr_el = orig_mesh->GetBdrElement(be);
        const int *verts = bdr_el->GetVertices();
        for (int i = 0; i < bdr_el->GetNVertices(); i++)
        {
          if (seen.insert(verts[i]).second)
          {
            crack_verts.push_back(verts[i]);
          }
        }
      }
    }
    const int nt = std::max(
        1, std::min(serial_mesh_threads, static_cast<int>(crack_verts.size() / 256)));
    while (component_pools.size() < static_cast<std::size_t>(nt))
    {
      component_pools.push_back(std::make_unique<std::pmr::unsynchronized_pool_resource>());
    }
    std::vector<std::unique_ptr<std::pmr::monotonic_buffer_resource>> bfs_buffers(nt);
    std::vector<std::vector<int>> ques(nt);
    for (auto &buffer : bfs_buffers)
    {
      buffer = std::make_unique<std::pmr::monotonic_buffer_resource>();
    }
    int new_nv_dups = 0;
    constexpr std::size_t batch = std::size_t(1) << 16;
    std::vector<std::vector<std::pmr::unordered_set<int>>> crack_vert_components(batch);
    for (std::size_t b = 0; b < crack_verts.size(); b += batch)
    {
      const std::size_t nb = std::min(batch, crack_verts.size() - b);
      ParallelForThreads(nb, nt,
                         [&](std::size_t k, int t)
                         {
                           auto &bfs_buffer = *bfs_buffers[t];
                           auto &que = ques[t];
                           const int v = crack_verts[b + k];
                           auto &components = crack_vert_components[k];
                           const int *elems = vert_to_elem->GetRow(v);
                           bfs_buffer.release();
                           std::pmr::unordered_set<int> unvisited(
                               elems, elems + vert_to_elem->RowSize(v), 0, std::hash<int>(),
                               std::equal_to<int>(), &bfs_buffer);
                           while (!unvisited.empty())
                           {
                             auto &component =
                                 components.emplace_back(component_pools[t].get());
                             component.reserve(unvisited.size());
                             que.clear();
                             que.push_back(*unvisited.begin());
                             unvisited.erase(unvisited.begin());
                             for (std::size_t head = 0; head < que.size(); head++)
                             {
                               // Process the current node.
                               int e = que[head];
                               component.insert(e);

                               // Add neighbors.
                               for (int k = etf_I[e]; k < etf_I[e + 1]; k++)
                               {
                                 const int nbr = face_nbr[k];
                                 if (nbr >= 0)
                                 {
                                   auto it = unvisited.find(nbr);
                                   if (it != unvisited.end())
                                   {
                                     que.push_back(nbr);
                                     unvisited.erase(it);
                                   }
                                 }
                               }
                             }
                           }
                         });
      for (std::size_t k = 0; k < nb; k++)
      {
        const int v = crack_verts[b + k];
        auto &components = crack_vert_components[k];
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
        components.clear();
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
          MFEM_ASSERT(crack_vert_duplicates.contains(v0) &&
                          crack_vert_duplicates.contains(v1),
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
      // Remove all seam edges which are on the "outside" of the crack (visited only
      // once).
      std::erase_if(coarse_crack_edge_to_be,
                    [](const auto &kv) { return kv.second.size() == 1; });
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
        const bool periodic = !iodata.boundaries.periodic.boundary_pairs.empty();
        if (!periodic)
        {
          // Non-periodic mesh: use MFEM's marked-edge conforming refinement.
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
          }
          if (mesh::CheckElements(*orig_mesh).has_simplices)
          {
            // Mark tetrahedral mesh for refinement before doing local refinement. This is a
            // bit of a strange pattern to override the standard conforming refinement of
            // the mfem::Mesh class. We want to implement our own edge marking of the
            // tetrahedra, so we move the mesh to a constructed derived class object, mark
            // it, and then move assign it to the original base class object before
            // refining. All of these moves should be cheap without any extra memory
            // allocation. Also, we mark the mesh every time to ensure multiple rounds of
            // refinement target the interior boundary (we don't care about preserving the
            // refinement hierarchy).
            constexpr bool refine = true, fix_orientation = false;
            EdgeRefinementMesh ref_mesh(std::move(*orig_mesh), coarse_crack_edge_to_be);
            if (ref_mesh.GetNodes())
            {
              ref_mesh.Finalize(refine, fix_orientation);
            }
            else
            {
              ref_mesh.Mark();
            }
            *orig_mesh = std::move(ref_mesh);
          }
          orig_mesh->GeneralRefinement(refinements, 0);
        }
        else
        {
          // Periodic mesh: MFEM's GeneralRefinement keys new midpoints on the collapsed
          // vertex pairs and cannot keep the two periodic copies of an identified seam edge
          // consistent, corrupting the seam and aborting later in STable3D. Use a manual,
          // purely local edge-fan bisection (LocalEdgeSplit) that touches only the ring of
          // elements around each flagged edge with no global conformity closure, so the
          // periodic seam is never disturbed. The inserted midpoints are decoupled by the
          // regular cracking pass on the subsequent retry.
          std::vector<std::pair<int, int>> split_edges;
          split_edges.reserve(coarse_crack_edge_to_be.size());
          for (const auto &[edge, adjacent_be] : coarse_crack_edge_to_be)
          {
            split_edges.emplace_back(edge.first, edge.second);
          }
          mesh::LocalEdgeSplit(orig_mesh, split_edges);
        }
        new_ne_ref += orig_mesh->GetNE() - ne;
        new_ref_its++;
        face_to_be.clear();
        return 0;  // Mesh was refined, start over
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
    if (!crack_bdr_elem.empty())
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
      if (face_to_be.contains(f))
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
          if (component.contains(e))
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
      if (crack_bdr_elem.contains(be))
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

  // MFEM identifies a face by its three smallest vertex numbers, so the two faces of a
  // cracked quadrilateral with three vertices on the crack front (at a corner of the crack)
  // would merge. Numbering the crack front vertices last puts each face's own vertex among
  // its three.
  if (std::any_of(crack_bdr_elem.begin(), crack_bdr_elem.end(), [&orig_mesh](int be)
                  { return orig_mesh->GetBdrElement(be)->GetNVertices() == 4; }))
  {
    const int nv = new_mesh->GetNV(), sdim = new_mesh->SpaceDimension();
    std::vector<bool> front(nv, false);
    for (const auto &[v, vert_components] : crack_vert_duplicates)
    {
      front[v] = vert_components.empty();
    }
    std::vector<int> perm(nv);
    int next = 0;
    for (bool last : {false, true})
    {
      for (int v = 0; v < nv; v++)
      {
        if (front[v] == last)
        {
          perm[v] = next++;
        }
      }
    }
    std::vector<double> x(static_cast<std::size_t>(nv) * sdim);
    for (int v = 0; v < nv; v++)
    {
      std::copy(new_mesh->GetVertex(v), new_mesh->GetVertex(v) + sdim,
                x.begin() + v * sdim);
    }
    for (int v = 0; v < nv; v++)
    {
      std::copy(x.begin() + v * sdim, x.begin() + (v + 1) * sdim,
                new_mesh->GetVertex(perm[v]));
    }
    auto renumber = [&perm](mfem::Element &el)
    {
      for (int j = 0; j < el.GetNVertices(); j++)
      {
        el.GetVertices()[j] = perm[el.GetVertices()[j]];
      }
    };
    for (int e = 0; e < new_mesh->GetNE(); e++)
    {
      renumber(*new_mesh->GetElement(e));
    }
    for (int be = 0; be < new_mesh->GetNBE(); be++)
    {
      renumber(*new_mesh->GetBdrElement(be));
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
    // The displacement moves the nodes of a mesh with nodes, not its vertices (see Load).
    if (!iodata.model.remove_curvature)
    {
      new_mesh->EnsureNodes();
    }
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

// Let a worker thread of a process bound to a single CPU run on any online CPU of the node
// (best effort).
void UnpinThread(std::thread &t)
{
#if defined(__linux__)
  cpu_set_t set;
  CPU_ZERO(&set);
  const long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
  for (long c = 0; c < ncpu && c < CPU_SETSIZE; c++)
  {
    CPU_SET(c, &set);
  }
  pthread_setaffinity_np(t.native_handle(), sizeof(set), &set);
#endif
}

// Reorder the elements of a mesh part along the Hilbert curve of their centroids through
// the cube [lo, lo + len]^sdim (that of the whole mesh, so that the element order of every
// part is a subsequence of one order), for the locality of the elements, edges and faces
// (and so of the degrees of freedom) of the parallel mesh built from it, which numbers its
// edges and faces in the order of the elements. The vertices keep their numbering (the
// global order), so that the edge orientations, given by their vertex numbers, and the
// shared entities are those without the reordering. Returns the new order (the element id
// of the part before the reordering, for each element), or an empty order for a part
// described by nodes, which is not reordered.
std::vector<int> ReorderMeshPart(mfem::MeshPart &part, const double *lo, double len)
{
  const int ne = part.num_elements, nv = part.num_vertices, sdim = part.space_dimension;
  if (part.nodes || ne == 0 || part.vertex_coordinates.Size() != nv * sdim)
  {
    return {};
  }

  // Elements are stored by geometry (the "by-type" ids), in the order of the geometries.
  const int dim = part.dimension;
  std::array<int, mfem::Geometry::NumGeom + 1> offsets{};
  for (int g = 0; g < mfem::Geometry::NumGeom; g++)
  {
    const bool elem_geom =
        g >= mfem::Geometry::DimStart[dim] && g < mfem::Geometry::DimStart[dim + 1];
    offsets[g + 1] =
        offsets[g] +
        (elem_geom ? part.entity_to_vertex[g].Size() / mfem::Geometry::NumVerts[g] : 0);
  }
  const bool have_map = (part.element_map.Size() == ne);
  std::vector<double> xc(static_cast<std::size_t>(ne) * sdim, 0.0);
  for (int i = 0; i < ne; i++)
  {
    const int bytype = have_map ? part.element_map[i] : i;
    int g = mfem::Geometry::DimStart[dim];
    while (bytype >= offsets[g + 1])
    {
      g++;
    }
    const int num_verts = mfem::Geometry::NumVerts[g];
    const int *v = part.entity_to_vertex[g].GetData() + (bytype - offsets[g]) * num_verts;
    for (int j = 0; j < num_verts; j++)
    {
      for (int d = 0; d < sdim; d++)
      {
        xc[i * sdim + d] += part.vertex_coordinates[v[j] * sdim + d];
      }
    }
    for (int d = 0; d < sdim; d++)
    {
      xc[i * sdim + d] /= num_verts;
    }
  }
  std::vector<int> order = mesh::HilbertOrder(xc.data(), ne, sdim, lo, len);

  // The new element order, as the map from the element ids to the by-type ids, and the
  // attributes in that order.
  mfem::Array<int> element_map(ne), attributes(ne);
  for (int i = 0; i < ne; i++)
  {
    element_map[i] = have_map ? part.element_map[order[i]] : order[i];
    attributes[i] = part.attributes[order[i]];
  }
  part.element_map = std::move(element_map);
  part.attributes = std::move(attributes);
  return order;
}

// Extract the parts of a partitioned serial mesh and write them in the MFEM parallel mesh
// format, with nt threads (each part is extracted independently). With a partitioning, the
// elements of each part are reordered (see ReorderMeshPart) along the curve through the
// bounding cube of the serial mesh, and origins returns the serial index of each element of
// each reordered part (empty for a part which is not reordered).
std::vector<std::string>
ExtractPartsAsText(mfem::Mesh &smesh, const mfem::MeshPartitioner &partitioner, int nparts,
                   int nt, const int *partitioning = nullptr,
                   std::vector<std::vector<std::int64_t>> *origins = nullptr)
{
  std::vector<std::string> so(nparts);
  std::array<double, 3> lo{};
  double len = 0.0;
  std::vector<std::int64_t> part_offsets, part_elements;
  if (partitioning)
  {
    // Bounding cube of the vertices, and the elements of each part in increasing order
    // (the order of the elements of an extracted part).
    const int sdim = smesh.SpaceDimension();
    std::array<double, 3> hi{};
    lo.fill(mfem::infinity());
    hi.fill(-mfem::infinity());
    for (int v = 0; v < smesh.GetNV(); v++)
    {
      const double *x = smesh.GetVertex(v);
      for (int d = 0; d < sdim; d++)
      {
        lo[d] = std::min(lo[d], x[d]);
        hi[d] = std::max(hi[d], x[d]);
      }
    }
    for (int d = 0; d < sdim; d++)
    {
      len = std::max(len, hi[d] - lo[d]);
    }
    part_offsets.assign(nparts + 1, 0);
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      part_offsets[partitioning[e] + 1]++;
    }
    std::partial_sum(part_offsets.begin(), part_offsets.end(), part_offsets.begin());
    part_elements.resize(smesh.GetNE());
    std::vector<std::int64_t> next(part_offsets.begin(), part_offsets.end() - 1);
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      part_elements[next[partitioning[e]]++] = e;
    }
    origins->assign(nparts, {});
  }
  auto Work = [&](int k)
  {
    for (int i = k; i < nparts; i += nt)
    {
      mfem::MeshPart part;
      partitioner.ExtractPart(i, part);
      if (partitioning)
      {
        const std::vector<int> order = ReorderMeshPart(part, lo.data(), len);
        auto &origin = (*origins)[i];
        origin.resize(order.size());
        for (std::size_t j = 0; j < order.size(); j++)
        {
          origin[j] = part_elements[part_offsets[i] + order[j]];
        }
      }
      std::ostringstream fo(std::stringstream::out);
      // fo << std::fixed;
      fo << std::scientific;
      fo.precision(MSH_FLT_PRECISION);
      part.Print(fo);
      so[i] = fo.str();
    }
  };
  if (nt > 1 && smesh.Dimension() >= 2)
  {
    // Created on first use by mfem::MeshPartitioner::ExtractPart: create it before the
    // threads do.
    smesh.GetEdgeVertexTable();
  }
  std::vector<std::thread> threads;
  threads.reserve(nt - 1);
  for (int k = 1; k < nt; k++)
  {
    threads.emplace_back(Work, k);
    UnpinThread(threads.back());
  }
  Work(0);
  for (auto &t : threads)
  {
    t.join();
  }
  return so;
}

std::unique_ptr<mfem::ParMesh>
DistributeMesh(MPI_Comm comm, std::unique_ptr<mfem::Mesh> &smesh, const int *partitioning,
               const std::string &output_dir, bool local_reorder)
{
  // Take a serial mesh and partitioning on the root process and construct the global
  // parallel mesh. For now, prefer the MPI-based version to the file IO one. When
  // constructing the ParMesh, we mark for refinement since refinement flags are not copied
  // from the serial mesh. Beware that mfem::ParMesh constructor argument order is not the
  // same as mfem::Mesh! Each processor's component gets sent as a byte string.
  constexpr bool generate_edges = false, refine = true, fix_orientation = false;
  std::unique_ptr<mfem::ParMesh> pmesh;

  // A serial mesh with first-order nodes, as from mfem::Mesh::EnsureNodes, is distributed
  // with its vertices moved to the nodes and without the nodes, which are added again to
  // the parallel mesh: extracting the parts is much cheaper without nodes, and the parallel
  // mesh is the same (it otherwise gets its vertices from the nodes of its part).
  bool first_order_nodes = false;
  if (Mpi::Root(comm) && smesh->GetNodes())
  {
    const auto *fespace = smesh->GetNodes()->FESpace();
    const mfem::H1_FECollection fec(1, smesh->Dimension());
    first_order_nodes = fespace->GetOrdering() == mfem::Ordering::byVDIM &&
                        fespace->GetVDim() == smesh->SpaceDimension() &&
                        !std::strcmp(fespace->FEColl()->Name(), fec.Name());
  }
  Mpi::Broadcast(1, &first_order_nodes, 0, comm);
  if (Mpi::Root(comm) && first_order_nodes)
  {
    const mfem::GridFunction &nodes = *smesh->GetNodes();
    const int sdim = smesh->SpaceDimension();
    mfem::Array<int> vdofs;
    for (int v = 0; v < smesh->GetNV(); v++)
    {
      nodes.FESpace()->GetVertexVDofs(v, vdofs);
      double *x = smesh->GetVertex(v);
      for (int d = 0; d < sdim; d++)
      {
        x[d] = nodes(vdofs[d]);
      }
    }
    smesh->SetCurvature(-1);
  }

  // Threads for the root: one per process on its node (PALACE_SETUP_THREADS overrides).
  int num_threads = 1;
  {
    MPI_Comm node_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, Mpi::Rank(comm), MPI_INFO_NULL,
                        &node_comm);
    num_threads = Mpi::Size(node_comm);
    MPI_Comm_free(&node_comm);
    if (const char *env = std::getenv("PALACE_SETUP_THREADS"))
    {
      num_threads = std::max(1, std::atoi(env));
    }
  }
  std::vector<std::vector<std::int64_t>> origins;
  if (Mpi::Root(comm))
  {
    mfem::MeshPartitioner partitioner(*smesh, Mpi::Size(comm),
                                      const_cast<int *>(partitioning));
    // The parts are extracted and written by threads on the CPUs of the other processes
    // of the node, which wait without spinning meanwhile (see below). Only parts without
    // nodes are extracted concurrently (the extraction of the nodes builds finite element
    // spaces from shared collections).
    const int nt = smesh->GetNodes() ? 1 : std::min(num_threads, Mpi::Size(comm));
    std::vector<std::string> so =
        ExtractPartsAsText(*smesh, partitioner, Mpi::Size(comm), nt,
                           local_reorder ? partitioning : nullptr, &origins);
    std::vector<MPI_Request> send_requests(2 * (Mpi::Size(comm) - 1), MPI_REQUEST_NULL);
    for (int i = 1; i < Mpi::Size(comm); i++)
    {
      int slen = static_cast<int>(so[i].length());
      MFEM_VERIFY(so[i].length() == static_cast<std::size_t>(slen),
                  "Overflow error distributing parallel mesh!");
      MPI_Isend(so[i].data(), slen, MPI_CHAR, i, i, comm, &send_requests[i - 1]);
      if (local_reorder)
      {
        // The serial index of each element of the part, after the part (same tag).
        MPI_Isend(origins[i].data(), static_cast<int>(origins[i].size()), MPI_INT64_T, i, i,
                  comm, &send_requests[Mpi::Size(comm) - 1 + i - 1]);
      }
    }
    smesh.reset();
    std::istringstream fi(so[0]);
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
    // Wait for the part without spinning, leaving the CPU to the root's threads.
    for (int flag = 0;;)
    {
      MPI_Iprobe(0, Mpi::Rank(comm), comm, &flag, &status);
      if (flag)
      {
        break;
      }
      std::this_thread::sleep_for(std::chrono::microseconds(200));
    }
    MPI_Get_count(&status, MPI_CHAR, &rlen);
    si.resize(rlen);
    MPI_Recv(si.data(), rlen, MPI_CHAR, 0, Mpi::Rank(comm), comm, MPI_STATUS_IGNORE);
    if (local_reorder)
    {
      MPI_Probe(0, Mpi::Rank(comm), comm, &status);
      MPI_Get_count(&status, MPI_INT64_T, &rlen);
      origins.assign(1, std::vector<std::int64_t>(rlen));
      MPI_Recv(origins[0].data(), rlen, MPI_INT64_T, 0, Mpi::Rank(comm), comm,
               MPI_STATUS_IGNORE);
    }
    std::istringstream fi(si);
    pmesh =
        std::make_unique<mfem::ParMesh>(comm, fi, refine, generate_edges, fix_orientation);
  }
  if (first_order_nodes)
  {
    pmesh->EnsureNodes();
  }
  // The serial index of each local element of a reordered part (none for a part
  // described by nodes, which is not reordered).
  std::vector<std::int64_t> origin;
  if (local_reorder && !origins.empty())
  {
    origin = std::move(origins[0]);
  }
  mesh::SetElementOrigin(*pmesh, std::move(origin));
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

void RegionRefine(const config::RefinementData &refinement, mfem::Mesh &mesh)
{
  int max_region_ref_levels = 0;
  for (const auto &box : refinement.GetBoxes())
  {
    max_region_ref_levels = std::max(max_region_ref_levels, box.ref_levels);
  }
  for (const auto &sphere : refinement.GetSpheres())
  {
    max_region_ref_levels = std::max(max_region_ref_levels, sphere.ref_levels);
  }
  if (max_region_ref_levels == 0)
  {
    return;
  }

  const auto element_types = mesh::CheckElements(mesh);
  if ((element_types.has_hexahedra || element_types.has_prisms ||
       element_types.has_pyramids) &&
      !mesh.Nonconforming())
  {
    // Region refinement of tensor meshes hangs hanging nodes; convert to NC. Simplex
    // meshes don't need NC for GeneralRefinement and stay conforming.
    mesh.EnsureNCMesh(true);
  }

  const bool use_nodes = (mesh.GetNodes() != nullptr);
  const int ref = use_nodes ? mesh.GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  const int dim = mesh.SpaceDimension();
  for (int region_ref_level = 0; region_ref_level < max_region_ref_levels;
       region_ref_level++)
  {
    // Mark elements for refinement: any vertex inside any active refinement region.
    mfem::Array<mfem::Refinement> refinements;
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      bool refine_elem = false;
      mfem::DenseMatrix pointmat;
      if (use_nodes)
      {
        mfem::ElementTransformation &T = *mesh.GetElementTransformation(i);
        mfem::RefinedGeometry *RefG =
            mfem::GlobGeometryRefiner.Refine(T.GetGeometryType(), ref);
        T.Transform(RefG->RefPts, pointmat);
      }
      else
      {
        const int *verts = mesh.GetElement(i)->GetVertices();
        const int nv = mesh.GetElement(i)->GetNVertices();
        pointmat.SetSize(dim, nv);
        for (int j = 0; j < nv; j++)
        {
          const double *coord = mesh.GetVertex(verts[j]);
          for (int d = 0; d < dim; d++)
          {
            pointmat(d, j) = coord[d];
          }
        }
      }
      for (const auto &box : refinement.GetBoxes())
      {
        if (region_ref_level >= box.ref_levels)
        {
          continue;
        }
        for (int j = 0; j < pointmat.Width(); j++)
        {
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
            refine_elem = true;
            break;
          }
        }
        if (refine_elem)
        {
          break;
        }
      }
      if (!refine_elem)
      {
        for (const auto &sphere : refinement.GetSpheres())
        {
          if (region_ref_level >= sphere.ref_levels)
          {
            continue;
          }
          for (int j = 0; j < pointmat.Width(); j++)
          {
            double dist = 0.0;
            for (int d = 0; d < pointmat.Height(); d++)
            {
              double s = pointmat(d, j) - sphere.center[d];
              dist += s * s;
            }
            if (dist <= sphere.r * sphere.r)
            {
              refine_elem = true;
              break;
            }
          }
          if (refine_elem)
          {
            break;
          }
        }
      }
      if (refine_elem)
      {
        refinements.Append(mfem::Refinement(i));
      }
    }
    // For tensor element meshes the refinement may make the mesh nonconforming (adds
    // hanging nodes).
    mesh.GeneralRefinement(refinements, -1);
  }
}

}  // namespace

}  // namespace palace
