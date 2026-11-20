// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/transientsolver.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "models/domainpostoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/dorfler.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

namespace palace
{

using json = nlohmann::json;

// Diagnostic (not for upstream): hashes of the mesh which do not depend on its partition,
// to compare builds. Vertex positions come from the nodes when the mesh has them (the
// vertices array is not kept consistent then). Per element: a hash of the bits of its
// vertex positions, sorted (geometry) and in its vertex order (marks), and the same with
// positions rounded to float (geometry32, marks32), summed over the (marked) elements.
std::array<std::uint64_t, 4> ElementHashes(const mfem::ParMesh &mesh, int e)
{
  auto Mix = [](std::uint64_t h, std::uint64_t v)
  {
    h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
  };
  const mfem::GridFunction *nodes = mesh.GetNodes();
  const mfem::FiniteElementSpace *fes = nodes ? nodes->FESpace() : nullptr;
  const mfem::Element *el = mesh.GetElement(e);
  std::vector<std::array<std::uint64_t, 3>> p64, p32;
  mfem::Array<int> vdofs;
  for (int j = 0; j < el->GetNVertices(); j++)
  {
    const int v = el->GetVertices()[j];
    std::array<double, 3> x = {0.0, 0.0, 0.0};
    if (fes && fes->GetNVDofs() > 0)
    {
      fes->GetVertexVDofs(v, vdofs);
      for (int d = 0; d < vdofs.Size(); d++)
      {
        x[d] = (*nodes)(vdofs[d]);
      }
    }
    else
    {
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
        x[d] = mesh.GetVertex(v)[d];
      }
    }
    std::array<std::uint64_t, 3> b64 = {0, 0, 0}, b32 = {0, 0, 0};
    for (int d = 0; d < 3; d++)
    {
      std::memcpy(&b64[d], &x[d], sizeof(double));
      const float f = static_cast<float>(x[d]);
      std::uint32_t u;
      std::memcpy(&u, &f, sizeof(float));
      b32[d] = u;
    }
    p64.push_back(b64);
    p32.push_back(b32);
  }
  auto H = [&Mix](const std::vector<std::array<std::uint64_t, 3>> &pts)
  {
    std::uint64_t h = 0;
    for (const auto &p : pts)
    {
      for (auto v : p)
      {
        h = Mix(h, v);
      }
    }
    return h;
  };
  std::array<std::uint64_t, 4> out;
  out[1] = H(p64);
  out[3] = H(p32);
  std::sort(p64.begin(), p64.end());
  std::sort(p32.begin(), p32.end());
  out[0] = H(p64);
  out[2] = H(p32);
  return out;
}

void PrintHashSums(const mfem::ParMesh &mesh, const char *what,
                   const mfem::Array<int> *subset)
{
  unsigned long long v[5] = {0, 0, 0, 0, 0};
  const int n = subset ? subset->Size() : mesh.GetNE();
  for (int i = 0; i < n; i++)
  {
    const auto h = ElementHashes(mesh, subset ? (*subset)[i] : i);
    for (int k = 0; k < 4; k++)
    {
      v[k] += h[k];
    }
  }
  v[4] = n;
  MPI_Allreduce(MPI_IN_PLACE, v, 5, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mesh.GetComm());
  Mpi::Print(" Mesh hashes {}: elements {:d}, geometry {:016x}, marks {:016x}, geometry32 "
             "{:016x}, marks32 {:016x}\n",
             what, v[4], v[0], v[1], v[2], v[3]);
}

void PrintMeshHashes(const mfem::ParMesh &mesh, const char *when)
{
  PrintHashSums(mesh, when, nullptr);
}

fs::path SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width)
{
  BlockTimer bt(Timer::IO);
  Mpi::Barrier(comm);  // Wait for all processes to write postprocessing files
  auto step_output = output_dir / fmt::format("iteration{:0{}d}", step, width);
  auto rel_step = step_output.filename();
  auto operation =
      fmt::format("archiving adaptive iteration output in \"{}\"", step_output.string());
  ApplyOnEachNodeFilesystem(
      comm, fmt::format("Creating adaptive iteration output \"{}\"", step_output.string()),
      [&]() { fs::create_directories(step_output); });
  ApplyOnceOnEachNodeFilesystem(
      comm, step_output / ".palace-archive-claim", operation,
      [&]()
      {
        // Snapshot the directory before changing it. Mutating a directory while iterating
        // it has filesystem-dependent behavior, especially when directory entries are
        // cached by a network filesystem.
        std::vector<fs::path> entries;
        for (const auto &entry : fs::directory_iterator(output_dir))
        {
          entries.push_back(entry.path());
        }

        for (const auto &path : entries)
        {
          const auto &fname = path.filename().string();
          if (fname.starts_with("iteration"))
          {
            continue;
          }

          // A path can disappear after the snapshot, and existing archive symlinks have
          // already been moved. Neither case leaves anything to archive.
          std::error_code ec;
          auto status = fs::symlink_status(path, ec);
          if (ec == std::errc::no_such_file_or_directory)
          {
            continue;
          }
          if (ec)
          {
            throw fs::filesystem_error("Could not inspect adaptive iteration output", path,
                                       ec);
          }
          if (!fs::exists(status) || fs::is_symlink(status))
          {
            continue;
          }

          auto dest = step_output / path.filename();
          if (fname == "palace.json")
          {
            // Metadata is written only by the global root. Copy it only there to avoid
            // concurrent copies when the output directory is shared.
            if (Mpi::Root(comm))
            {
              fs::copy(path, dest, fs::copy_options::overwrite_existing);
            }
          }
          else if (fname.ends_with("_resolved.json"))
          {
            // The resolved configuration is a global record of the run, not per-iteration
            // output. Leave it in the top-level output folder rather than moving it into an
            // iteration subfolder and symlinking back.
            continue;
          }
          else
          {
            // Move to the iteration subfolder and leave a relative symlink behind
            // so that the output directory always has accessible results. Remove
            // any existing destination first (e.g. from a previous run).
            fs::remove_all(dest);
            fs::rename(path, dest);
            fs::create_symlink(rel_step / path.filename(), path);
          }
        }
      });
  return step_output;
}

namespace
{

json LoadMetadata(const fs::path &post_dir)
{
  std::string path = post_dir / "palace.json";
  std::ifstream fi(path);
  if (!fi.is_open())
  {
    MFEM_ABORT("Unable to open metadata file \"" << path << "\"!");
  }
  return json::parse(fi);
}

void WriteMetadata(const fs::path &post_dir, const json &meta)
{
  // Write to a temporary file and atomically rename it into place so a concurrent reader
  // (e.g. an external monitor polling palace.json mid-run) never observes a partial file.
  std::string path = post_dir / "palace.json";
  std::string tmp_path = post_dir / "palace.json.tmp";
  {
    std::ofstream fo(tmp_path);
    if (!fo.is_open())
    {
      MFEM_ABORT("Unable to open metadata file \"" << tmp_path << "\"!");
    }
    fo << meta.dump(2) << '\n';
  }
  fs::rename(tmp_path, path);
}

// Map an MFEM geometry type to its lowercase JSON key used in the SavedAdaptedMesh block.
const char *GeometryKey(mfem::Geometry::Type geom)
{
  switch (geom)
  {
    case mfem::Geometry::TRIANGLE:
      return "triangle";
    case mfem::Geometry::SQUARE:
      return "quadrilateral";
    case mfem::Geometry::TETRAHEDRON:
      return "tetrahedron";
    case mfem::Geometry::CUBE:
      return "hexahedron";
    case mfem::Geometry::PRISM:
      return "prism";
    case mfem::Geometry::PYRAMID:
      return "pyramid";
    default:
      return "unknown";
  }
}

// Returns an array of indices corresponding to marked elements.
mfem::Array<int> MarkedElements(const Vector &e, double threshold)
{
  mfem::Array<int> ind;
  ind.Reserve(e.Size());
  for (int i = 0; i < e.Size(); i++)
  {
    if (e[i] >= threshold)
    {
      ind.Append(i);
    }
  }
  return ind;
}

}  // namespace

BaseSolver::BaseSolver(const IoData &iodata, bool root, int size, int num_thread,
                       const char *git_tag)
  : iodata(iodata), post_dir(iodata.problem.output), root(root)
{
  // Initialize simulation metadata for this simulation.
  if (root)
  {
    json meta;
    if (git_tag)
    {
      meta["GitTag"] = std::string(git_tag);
    }
    if (size > 0)
    {
      meta["Problem"]["MPISize"] = size;
    }
    if (num_thread > 0)
    {
      meta["Problem"]["OpenMPThreads"] = num_thread;
    }
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                            MPI_Comm comm) const
{
  if (!(iodata.model.Lc > 0.0))
  {
    iodata.model.Lc = mesh::ComputeReferenceLength(smesh, comm);
  }
  iodata.NondimensionalizeInputs(smesh);
}

bool BaseSolver::ContinueAdaptation(int it, long long int ntdof, double err) const
{
  // Stop when out of iterations, out of unknowns if a limit was set, or below the
  // tolerance.
  const auto &refinement = iodata.model.refinement;
  return it < refinement.max_it &&
         !(refinement.max_size > 0 && ntdof > refinement.max_size) && err >= refinement.tol;
}

void BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = [&]()
  {
    if (refinement.max_it > 0 && iodata.problem.type == ProblemType::TRANSIENT)
    {
      Mpi::Warning("AMR is not currently supported for transient simulations!\n");
      return false;
    }
    return (refinement.max_it > 0);
  }();
  if (use_amr && mesh.size() > 1)
  {
    Mpi::Print("\nFlattening mesh sequence:\n AMR will start from the final mesh in "
               "the sequence of a priori refinements\n");
    mesh.erase(mesh.begin(), mesh.end() - 1);
  }
  MPI_Comm comm = mesh.back()->GetComm();

  // Perform initial solve and estimation.
  PrintMeshHashes(*mesh.back(), "at load");
  auto [indicators, ntdof] = Solve(mesh);
  MFEM_VERIFY(solve_converged_,
              "Initial solve did not converge; no converged result to fall back on!");
  double err = indicators.Norml2(comm);

  // True topological entity counts of the final adapted mesh.
  mesh::MeshEntityCounts mesh_counts;

  // Main AMR loop.
  int it = 0;
  while (use_amr && ContinueAdaptation(it, ntdof, err))
  {
    // Print timing summary.
    Mpi::Print(comm, "\nCumulative timing statistics:\n");
    BlockTimer::Finalize(comm);
    BlockTimer::Print(comm);
    auto peak_mem = memory_reporting::GetPeakMemoryStats(comm);
    auto peak_node_mem = memory_reporting::GetPeakNodeMemoryStats(comm);
    memory_reporting::PrintMemoryUsage(comm, peak_mem);
    memory_reporting::PrintMemoryUsage(comm, peak_node_mem);
    SaveMetadata(BlockTimer::GlobalTimer());
    SaveMetadata(peak_mem);
    SaveMetadata(peak_node_mem);

    BlockTimer bt(Timer::ADAPTATION);
    Mpi::Print("\nAdaptive mesh refinement (AMR) iteration {:d}:\n"
               " Indicator norm = {:.3e}, global unknowns = {:d}\n"
               " Max. iterations = {:d}, tol. = {:.3e}{}\n",
               ++it, err, ntdof, refinement.max_it, refinement.tol,
               (refinement.max_size > 0
                    ? ", max. size = " + std::to_string(refinement.max_size)
                    : ""));

    // Optionally save off the previous solution, then record it as completed. The index is
    // written only once the "iterationXX" archive is complete, so a reader of palace.json
    // never sees an iteration whose output is still being moved.
    if (refinement.save_adapt_iterations)
    {
      SaveAdaptationIteration(
          it, SaveIteration(comm, post_dir, it,
                            1 + static_cast<int>(std::log10(refinement.max_it))));
    }
    else
    {
      SaveAdaptationIteration(it);
    }

    // Mark.
    const auto marked_elements = [&comm, &refinement](const auto &indicators)
    {
      const auto [threshold, marked_error] = utils::ComputeDorflerThreshold(
          comm, indicators.Local(), refinement.update_fraction);
      const auto marked_elements = MarkedElements(indicators.Local(), threshold);
      const auto [glob_marked_elements, glob_elements] =
          linalg::GlobalSize2(comm, marked_elements, indicators.Local());
      Mpi::Print(
          " Marked {:d}/{:d} elements for refinement ({:.2f}% of the error, θ = {:.2f})\n",
          glob_marked_elements, glob_elements, 100 * marked_error,
          refinement.update_fraction);
      return marked_elements;
    }(indicators);

    // Refine.
    {
      mfem::ParMesh &fine_mesh = *mesh.back();
      const auto initial_elem_count = fine_mesh.GetGlobalNE();
      PrintHashSums(fine_mesh, "of the marked elements", &marked_elements);
      fine_mesh.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
      if (fine_mesh.Conforming())
      {
        // New marks, the same for any number of processes; fields on the mesh are lost.
        fine_mesh.Finalize(/*refine=*/true);
      }
      UpdateAfterMeshChange(*mesh.back());
      PrintMeshHashes(fine_mesh, "after refinement");
      const auto final_elem_count = fine_mesh.GetGlobalNE();
      Mpi::Print(" {} mesh refinement added {:d} elements (initial = {:d}, final = {:d})\n",
                 fine_mesh.Nonconforming() ? "Nonconforming" : "Conforming",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }

    // Optionally rebalance and write the adapted mesh to file.
    {
      const auto result = mesh::RebalanceMesh(iodata, *mesh.back(), &mesh_counts);
      if (result.ratio_before > refinement.maximum_imbalance)
      {
        Mpi::Print(" {}\n", mesh::RebalanceMessage(result, refinement.maximum_imbalance));
      }
      mesh.back()->Update();
      UpdateAfterMeshChange(*mesh.back());
      PrintMeshHashes(*mesh.back(), "after rebalance");
    }

    // Record the adapted mesh's topology now, paired with the mesh file RebalanceMesh
    // wrote, so it matches even if the following solve does not converge (SavedAdaptedMesh
    // contract).
    if (refinement.save_adapt_mesh)
    {
      mesh::CompleteMeshEntityCounts(*mesh.back(), mesh_counts);
      if (mesh_counts.valid)
      {
        SaveMetadata(mesh_counts);
      }
    }

    // Print statistics (element counts, size h, and shape regularity kappa) for the
    // newly-refined mesh so the evolution of mesh quality under AMR is visible.
    mesh::PrintMeshInfo(*mesh.back(), iodata, /*full=*/false);

    // Solve + estimate.
    Mpi::Print("\nProceeding with solve/estimate iteration {}...\n", it + 1);
    const auto prev_ntdof = ntdof;
    std::tie(indicators, ntdof) = Solve(mesh);
    if (!solve_converged_)
    {
      // Keep the previous converged iteration for the summary/CSVs (restore its DOF count);
      // the refined mesh file and its topology (recorded above) describe the failed solve.
      Mpi::Warning(
          comm,
          "Solve did not converge after refinement iteration {:d}; halting AMR and "
          "keeping the last converged iteration!\n",
          it);
      ntdof = prev_ntdof;
      break;
    }
    err = indicators.Norml2(comm);
  }

  // The last converged solve is never archived: its output stays in post_dir. A failed
  // solve after refinement leaves iteration it as the last converged one, already recorded
  // above.
  if (solve_converged_)
  {
    SaveAdaptationIteration(it + 1);
  }

  Mpi::Print("\nCompleted {:d} iteration{} of adaptive mesh refinement (AMR):\n"
             " Indicator norm = {:.3e}, global unknowns = {:d}\n"
             " Max. iterations = {:d}, tol. = {:.3e}{}\n",
             it, (it == 1 ? "" : "s"), err, ntdof, refinement.max_it, refinement.tol,
             (refinement.max_size > 0
                  ? ", max. size = " + std::to_string(refinement.max_size)
                  : ""));
}

void BaseSolver::SaveMetadata(const FiniteElementSpaceHierarchy &fespaces) const
{
  const auto &fespace = fespaces.GetFinestFESpace();
  HYPRE_BigInt ne = fespace.GetParMesh().GetNE();
  Mpi::GlobalSum(1, &ne, fespace.GetComm());
  std::vector<HYPRE_BigInt> ndofs(fespaces.GetNumLevels());
  for (std::size_t l = 0; l < fespaces.GetNumLevels(); l++)
  {
    ndofs[l] = fespaces.GetFESpaceAtLevel(l).GlobalTrueVSize();
  }
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["Problem"]["MeshElements"] = ne;
    meta["Problem"]["DegreesOfFreedom"] = ndofs.back();
    meta["Problem"]["MultigridDegreesOfFreedom"] = ndofs;
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveAdaptationIteration(int iteration, const fs::path &archive_dir) const
{
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["Problem"]["Iteration"] = iteration;
    if (!archive_dir.empty())
    {
      // Update the archived copy before post_dir, which is the one readers poll.
      WriteMetadata(archive_dir, meta);
    }
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(const mesh::MeshEntityCounts &counts) const
{
  // Unlike the other SaveMetadata overloads (called on every rank, some running MPI
  // collectives before writing), this one is invoked on the root rank only, since
  // MeshEntityCounts::valid is set only on root. It must contain no MPI-collective calls.
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    json &saved = meta["SavedAdaptedMesh"];
    saved["Dimension"] = counts.dim;
    saved["TrueVertices"] = counts.true_vertices;
    saved["TrueEdges"] = counts.true_edges;
    // Per-geometry counts emit only the geometries that are present; a missing key means
    // zero by contract. TrueFaces is omitted entirely in 2D (the map is empty there).
    if (!counts.true_faces.empty())
    {
      json faces = json::object();
      for (const auto &[geom, count] : counts.true_faces)
      {
        faces[GeometryKey(geom)] = count;
      }
      saved["TrueFaces"] = faces;
    }
    json cells = json::object();
    for (const auto &[geom, count] : counts.cells)
    {
      cells[GeometryKey(geom)] = count;
    }
    saved["Cells"] = cells;
    saved["DomainAttributes"] = counts.domain_attributes;
    saved["BoundaryAttributes"] = counts.boundary_attributes;
    WriteMetadata(post_dir, meta);
  }
}

template <typename SolverType>
void BaseSolver::SaveMetadata(const SolverType &ksp) const
{
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["LinearSolver"]["TotalSolves"] = ksp.NumTotalMult();
    meta["LinearSolver"]["TotalIts"] = ksp.NumTotalMultIterations();
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(const Timer &timer) const
{
  if (root)
  {
    constexpr double to_mb = 1.0 / (1024.0 * 1024.0);
    auto red = BlockTimer::GetReductions();

    json meta = LoadMetadata(post_dir);
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      std::string key(Timer::descriptions[i]);
      std::erase_if(key, isspace);
      meta["ElapsedTime"]["Durations"][key] = timer.Data(static_cast<Timer::Index>(i));
      meta["ElapsedTime"]["Counts"][key] = timer.Counts(static_cast<Timer::Index>(i));
      meta["PeakMemoryGrowthMegabytes"]["Min"][key] = red.rank_mem.min[i] * to_mb;
      meta["PeakMemoryGrowthMegabytes"]["Max"][key] = red.rank_mem.max[i] * to_mb;
      meta["PeakMemoryGrowthMegabytes"]["Sum"][key] = red.rank_mem.sum[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Min"][key] = red.node_mem.min[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Max"][key] = red.node_mem.max[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Sum"][key] = red.node_mem.sum[i] * to_mb;
    }
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(const memory_reporting::MemoryStats &peak_memory) const
{
  if (root)
  {
    constexpr double to_mb = 1.0 / (1024.0 * 1024.0);
    json meta = LoadMetadata(post_dir);

    // Determine key name based on whether this is per-node or per-rank stats
    std::string key = (peak_memory.label.find("per-node") != std::string::npos)
                          ? "PeakNodeMemoryMegabytes"
                          : "PeakMemoryMegabytes";

    meta[key]["Min"] = peak_memory.min * to_mb;
    meta[key]["Max"] = peak_memory.max * to_mb;
    meta[key]["Average"] = peak_memory.avg * to_mb;
    meta[key]["Total"] = peak_memory.sum * to_mb;
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(const PortExcitations &excitation_helper) const
{
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    meta["Excitations"] = excitation_helper;
    WriteMetadata(post_dir, meta);
  }
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
