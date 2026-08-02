// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <array>
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

void SaveIteration(MPI_Comm comm, const fs::path &output_dir, int step, int width)
{
  BlockTimer bt(Timer::IO);
  Mpi::Barrier(comm);  // Wait for all processes to write postprocessing files
  auto step_output = output_dir / fmt::format("iteration{:0{}d}", step, width);
  auto rel_step = step_output.filename();
  ApplyOnEachNodeFilesystem(
      comm,
      fmt::format("Archiving adaptive iteration output in \"{}\"", step_output.string()),
      [&]()
      {
        // On a shared filesystem, the global root moves the output before the other node
        // roots inspect it, so they see symlinks and have nothing left to move. On a
        // node-local filesystem, each node root archives that node's rank files.
        fs::create_directories(step_output);
        for (const auto &f : fs::directory_iterator(output_dir))
        {
          const auto &fname = f.path().filename().string();
          if (fname.rfind("iteration") == 0)
          {
            continue;
          }
          auto dest = step_output / f.path().filename();
          if (f.is_symlink())
          {
            // Skip symlinks left from a previous iteration's save. They still point
            // to valid data and will be overwritten by the next solve.
            continue;
          }
          else if (fname == "palace.json")
          {
            // Metadata is written only by the global root. Copy it only there to avoid
            // concurrent copies when the output directory is shared.
            if (Mpi::Root(comm))
            {
              fs::copy(f, dest, fs::copy_options::overwrite_existing);
            }
          }
          else if (fname.size() >= 14 &&
                   fname.compare(fname.size() - 14, 14, "_resolved.json") == 0)
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
            fs::rename(f, dest);
            fs::create_symlink(rel_step / f.path().filename(), f.path());
          }
        }
      });
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
  std::string path = post_dir / "palace.json";
  std::ofstream fo(path);
  if (!fo.is_open())
  {
    MFEM_ABORT("Unable to open metadata file \"" << path << "\"!");
  }
  fo << meta.dump(2) << '\n';
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

void BaseSolver::ProcessPartitionedMesh(const mfem::ParMesh &,
                                        const mesh::PartitionMetadata &) const
{
}

mesh::PartitionMetadata BaseSolver::GetSourceEntityMetadata() const
{
  return {};
}

mfem::Array<int> BaseSolver::GetRefinementProtection(const mfem::ParMesh &,
                                                     bool *conforming,
                                                     mfem::Array<int> *repair) const
{
  if (conforming)
  {
    *conforming = true;
  }
  if (repair)
  {
    repair->SetSize(0);
  }
  return {};
}

void BaseSolver::ProcessRefinedMesh(const mfem::ParMesh &) const {}

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
  auto [indicators, ntdof] = Solve(mesh);
  const auto ProtectIndicators = [this, &mesh](ErrorIndicator &indicators)
  {
    bool closure_conforming = true;
    auto protection = GetRefinementProtection(mesh.back()->Get(), &closure_conforming);
    MFEM_VERIFY(protection.Size() == 0 || protection.Size() == mesh.back()->GetNE(),
                "AMR refinement-protection marker has an invalid size!");
    MFEM_VERIFY(closure_conforming,
                "The singular refinement closure is nonconforming before AMR!");
    if (protection.Size() > 0 && mesh.back()->Get().Nonconforming())
    {
      // Custom singular spaces do not yet have coarse/fine trace constraints. Keep the
      // enriched patch and its one-face buffer fixed in h so all enriched interfaces remain
      // conforming, and adapt using the smooth-remainder estimator outside that closure.
      indicators.ZeroElements(protection);
    }
    else
    {
      // Conforming refinement can pass through enriched elements because it introduces no
      // coarse/fine singular traces.
      protection.SetSize(0);
    }
    return protection;
  };
  auto protection = ProtectIndicators(indicators);
  double err = indicators.Norml2(comm);

  // Collection of all tests that might exhaust resources.
  auto ExhaustedResources = [&refinement](auto it, auto ntdof)
  {
    bool ret = false;
    // Run out of iterations.
    ret |= (it >= refinement.max_it);
    // Run out of DOFs if a limit was set.
    ret |= (refinement.max_size > 0 && ntdof > refinement.max_size);
    return ret;
  };

  // Main AMR loop.
  int it = 0;
  while (use_amr && !ExhaustedResources(it, ntdof) && err >= refinement.tol)
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

    // Optionally save off the previous solution.
    if (refinement.save_adapt_iterations)
    {
      SaveIteration(comm, post_dir, it,
                    1 + static_cast<int>(std::log10(refinement.max_it)));
    }

    // Mark. A nonempty driver marker identifies a fixed-h refinement closure.
    // For singular enrichment this is the enriched patch plus one face-neighbor layer.
    if (protection.Size() > 0)
    {
      int closure_elements = 0;
      for (int element = 0; element < protection.Size(); element++)
      {
        closure_elements += protection[element] != 0;
      }
      Mpi::GlobalSum(1, &closure_elements, comm);
      Mpi::Print(" Singular refinement closure contains {:d} elements\n", closure_elements);
    }

    const auto marked_elements = [&comm, &refinement, &protection](const auto &indicators)
    {
      const auto [threshold, marked_error] = utils::ComputeDorflerThreshold(
          comm, indicators.Local(), refinement.update_fraction);
      auto marked_elements = MarkedElements(indicators.Local(), threshold);
      for (int element : marked_elements)
      {
        MFEM_VERIFY(protection.Size() == 0 || !protection[element],
                    "AMR marked an element in a protected singular refinement closure!");
      }
      const auto [glob_marked_elements, glob_elements] =
          linalg::GlobalSize2(comm, marked_elements, indicators.Local());
      Mpi::Print(
          " Marked {:d}/{:d} elements for refinement ({:.2f}% of the error, θ = {:.2f})\n",
          glob_marked_elements, glob_elements, 100 * marked_error,
          refinement.update_fraction);
      return marked_elements;
    }(indicators);
    int global_marked_elements = marked_elements.Size();
    Mpi::GlobalSum(1, &global_marked_elements, comm);
    if (global_marked_elements == 0)
    {
      Mpi::Warning(
          comm,
          "AMR stopped because no unprotected elements were marked for refinement.\n");
      break;
    }

    // Refine.
    {
      mfem::ParMesh &fine_mesh = *mesh.back();
      const auto initial_elem_count = fine_mesh.GetGlobalNE();
      fine_mesh.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
      const auto final_elem_count = fine_mesh.GetGlobalNE();
      Mpi::Print(" {} mesh refinement added {:d} elements (initial = {:d}, final = {:d})\n",
                 fine_mesh.Nonconforming() ? "Nonconforming" : "Conforming",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }
    ProcessRefinedMesh(mesh.back()->Get());
    bool refined_closure_conforming = true;
    const auto refined_protection =
        GetRefinementProtection(mesh.back()->Get(), &refined_closure_conforming);
    MFEM_VERIFY(refined_protection.Size() == 0 ||
                    refined_protection.Size() == mesh.back()->GetNE(),
                "Refined AMR protection marker has an invalid size!");
    MFEM_VERIFY(refined_closure_conforming,
                "AMR introduced an unsupported coarse/fine singular trace!");

    // Optionally rebalance and write the adapted mesh to file.
    if (RebalanceRefinedMesh())
    {
      mesh::PartitionMetadata source_metadata;
      auto *metadata = [&]() -> mesh::PartitionMetadata *
      {
        if (!RequiresSourceSerialMeshMetadata())
        {
          return nullptr;
        }
        source_metadata = GetSourceEntityMetadata();
        MFEM_VERIFY(source_metadata.source_vertex_ids.size() ==
                            static_cast<std::size_t>(mesh.back()->Get().GetNV()) &&
                        source_metadata.source_element_ids.size() ==
                            static_cast<std::size_t>(mesh.back()->Get().GetNE()),
                    "AMR rebalancing received incomplete source-entity metadata!");
        return &source_metadata;
      }();
      const auto ratio_pre = mesh::RebalanceMesh(iodata, *mesh.back(), metadata);
      if (ratio_pre > refinement.maximum_imbalance)
      {
        if (metadata)
        {
          ProcessPartitionedMesh(mesh.back()->Get(), *metadata);
        }
        int min_elem, max_elem;
        min_elem = max_elem = mesh.back()->GetNE();
        Mpi::GlobalMin(1, &min_elem, comm);
        Mpi::GlobalMax(1, &max_elem, comm);
        const auto ratio_post = double(max_elem) / min_elem;
        Mpi::Print(" Rebalanced mesh: Ratio {:.3f} exceeded max. allowed value {:.3f} "
                   "(new ratio = {:.3f})\n",
                   ratio_pre, refinement.maximum_imbalance, ratio_post);
      }
    }
    else
    {
      Mpi::Print(" Skipping mesh rebalancing to preserve singular-feature identities\n");
    }
    bool rebalanced_closure_conforming = true;
    const auto rebalanced_protection =
        GetRefinementProtection(mesh.back()->Get(), &rebalanced_closure_conforming);
    MFEM_VERIFY(rebalanced_protection.Size() == 0 ||
                    rebalanced_protection.Size() == mesh.back()->GetNE(),
                "Rebalanced AMR protection marker has an invalid size!");
    MFEM_VERIFY(rebalanced_closure_conforming,
                "Mesh rebalancing introduced a nonconforming singular interface!");
    mesh.back()->Update();

    // Print statistics (element counts, size h, and shape regularity kappa) for the
    // newly-refined mesh so the evolution of mesh quality under AMR is visible.
    mesh::PrintMeshInfo(*mesh.back(), iodata, /*full=*/false);

    // Solve + estimate.
    Mpi::Print("\nProceeding with solve/estimate iteration {}...\n", it + 1);
    std::tie(indicators, ntdof) = Solve(mesh);
    protection = ProtectIndicators(indicators);
    err = indicators.Norml2(comm);
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

template <typename SolverType>
void BaseSolver::SaveMetadata(const SolverType &ksp) const
{
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["LinearSolver"]["TotalSolves"] = ksp.NumTotalMult();
    meta["LinearSolver"]["TotalIts"] = ksp.NumTotalMultIterations();
    meta["LinearSolver"]["FailedSolves"] = ksp.NumFailedMult();
    meta["LinearSolver"]["Converged"] = (ksp.NumFailedMult() == 0);
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
      auto key = Timer::descriptions[i];
      key.erase(std::remove_if(key.begin(), key.end(), isspace), key.end());
      meta["ElapsedTime"]["Durations"][key] = timer.Data((Timer::Index)i);
      meta["ElapsedTime"]["Counts"][key] = timer.Counts((Timer::Index)i);
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

void BaseSolver::SaveMetadata(std::string_view section, const nlohmann::json &data) const
{
  if (root)
  {
    MFEM_VERIFY(!section.empty(), "Metadata section name must not be empty!");
    json meta = LoadMetadata(post_dir);
    meta[section] = data;
    WriteMetadata(post_dir, meta);
  }
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
