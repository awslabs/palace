// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <sstream>
#include <vector>
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

void ReportStableSeedIds(MPI_Comm comm, const mfem::ParMesh &mesh,
                         const mfem::Array<int> &marked,
                         const mesh::PartitionMetadata &metadata)
{
  // Element numbers may be reassigned after refinement/repartitioning. A simplex key made
  // from sorted persistent vertex ancestry IDs remains decomposition-independent across all
  // AMR passes.
  int ranks_with_metadata = !metadata.source_vertex_ids.empty();
  Mpi::GlobalSum(1, &ranks_with_metadata, comm);
  if (ranks_with_metadata == 0)
  {
    return;
  }
  MFEM_VERIFY(mesh.GetNV() == 0 || metadata.source_vertex_ids.size() ==
                                       static_cast<std::size_t>(mesh.GetNV()),
              "AMR stable seed diagnostics received an incomplete vertex identity map!");
  constexpr int record_size = 5;  // vertex count followed by at most four simplex vertices
  std::vector<std::int64_t> local(record_size * marked.Size(), -1);
  mfem::Array<int> vertices;
  for (int i = 0; i < marked.Size(); i++)
  {
    mesh.GetElementVertices(marked[i], vertices);
    MFEM_VERIFY(vertices.Size() >= 2 && vertices.Size() <= 4,
                "AMR stable seed diagnostics require simplex elements!");
    local[record_size * i] = vertices.Size();
    std::vector<std::int64_t> key(vertices.Size());
    for (int j = 0; j < vertices.Size(); j++)
    {
      key[j] = metadata.source_vertex_ids[vertices[j]];
    }
    std::sort(key.begin(), key.end());
    std::copy(key.begin(), key.end(), local.begin() + record_size * i + 1);
  }
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);
  const int local_size = static_cast<int>(local.size());
  std::vector<int> counts(size), offsets(size + 1);
  MPI_Gather(&local_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
  std::vector<std::int64_t> gathered;
  if (rank == 0)
  {
    for (int r = 0; r < size; r++)
    {
      offsets[r + 1] = offsets[r] + counts[r];
    }
    gathered.resize(offsets.back());
  }
  MPI_Gatherv(local.data(), local_size, MPI_INT64_T, gathered.data(), counts.data(),
              offsets.data(), MPI_INT64_T, 0, comm);
  if (rank == 0)
  {
    MFEM_VERIFY(gathered.size() % record_size == 0,
                "AMR stable seed diagnostics gathered malformed simplex keys!");
    std::vector<std::array<std::int64_t, record_size>> keys(gathered.size() / record_size);
    for (std::size_t i = 0; i < keys.size(); i++)
    {
      std::copy_n(gathered.begin() + record_size * i, record_size, keys[i].begin());
    }
    std::sort(keys.begin(), keys.end());
    std::ostringstream ids;
    std::uint64_t hash = 1469598103934665603ULL;
    constexpr std::size_t maximum_printed_keys = 64;
    for (std::size_t i = 0; i < keys.size(); i++)
    {
      if (i < maximum_printed_keys)
      {
        ids << (i == 0 ? "" : ",") << "{";
        for (int j = 0; j < keys[i][0]; j++)
        {
          ids << (j == 0 ? "" : ":") << keys[i][j + 1];
        }
        ids << "}";
      }
      for (const auto value : keys[i])
      {
        hash ^= static_cast<std::uint64_t>(value);
        hash *= 1099511628211ULL;
      }
    }
    Mpi::Print(comm, " Stable AMR seed simplices: {:d} keys, hash {:016x}, sample [{}{}]\n",
               keys.size(), hash, ids.str(),
               keys.size() > maximum_printed_keys ? ",..." : "");
  }
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

void BaseSolver::ObserveRefinementAncestry(const mfem::ParMesh &) const {}

mfem::Array<int> BaseSolver::GetEnrichedElements(const mfem::ParMesh &) const
{
  return {};
}

void BaseSolver::ReportTraceComponents(const mfem::ParMesh &,
                                       const mfem::Array<int> &) const
{
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
  if (use_amr && iodata.solver.singular_elements.Enabled() &&
      mesh.back()->Get().Dimension() == 3 && refinement.nonconformal &&
      refinement.max_nc_levels > 0 && !refinement.singular_repair)
  {
    MFEM_ABORT("Three-dimensional nonconforming singular-element AMR requires "
               "Model.Refinement.MaxNCLevels = 0 or SingularRepair = true. "
               "Two-dimensional electrostatic line-tip AMR supports bounded hanging-node "
               "levels through closure expansion and localized repair.");
  }

  // Report the distribution of the raw indicator across the enriched region, its
  // conservative one-face closure, and the exterior, before any protection mask is
  // applied. This is diagnostic only: it does not change marking. The enriched region is
  // where the singular basis is active; the closure is the additional buffer whose
  // conformity the current implementation requires.
  const auto ReportIndicatorRegions =
      [this, &mesh, comm](const ErrorIndicator &indicators,
                          const mfem::Array<int> &protection, double update_fraction)
  {
    if (protection.Size() == 0)
    {
      return;
    }
    const auto &local = indicators.Local();
    MFEM_VERIFY(local.Size() == protection.Size(),
                "Indicator region diagnostics received an inconsistent marker!");
    const auto enriched = GetEnrichedElements(mesh.back()->Get());
    const bool have_enriched = (enriched.Size() == protection.Size());

    // Region 0: enriched, 1: closure buffer (protected, not enriched), 2: exterior.
    constexpr int regions = 3;
    std::array<long long, regions> counts{};
    std::array<double, regions> squared{};
    std::array<double, regions> maximum{};
    maximum.fill(0.0);
    const auto region_of = [&](int element)
    {
      if (have_enriched && enriched[element])
      {
        return 0;
      }
      return protection[element] ? 1 : 2;
    };
    for (int element = 0; element < local.Size(); element++)
    {
      const int region = region_of(element);
      counts[region]++;
      squared[region] += local[element] * local[element];
      maximum[region] = std::max(maximum[region], local[element]);
    }

    // Dörfler selection on the unmasked indicator, so the report shows where the
    // refinement budget would go if nothing were protected.
    const auto [threshold, marked_error] =
        utils::ComputeDorflerThreshold(comm, local, update_fraction);
    std::array<long long, regions> selected{};
    for (int element = 0; element < local.Size(); element++)
    {
      if (local[element] >= threshold)
      {
        selected[region_of(element)]++;
      }
    }

    Mpi::GlobalSum(regions, counts.data(), comm);
    Mpi::GlobalSum(regions, squared.data(), comm);
    Mpi::GlobalMax(regions, maximum.data(), comm);
    Mpi::GlobalSum(regions, selected.data(), comm);
    const double total_squared = squared[0] + squared[1] + squared[2];
    const long long total_selected = selected[0] + selected[1] + selected[2];
    constexpr std::array<const char *, regions> labels{"enriched", "closure buffer",
                                                       "exterior"};
    Mpi::Print(" Raw indicator by region (before protection, θ = {:.2f}):\n",
               update_fraction);
    for (int region = 0; region < regions; region++)
    {
      Mpi::Print("  {:14s}: {:d} elements, {:.2f}% of squared indicator, max. {:.3e}, "
                 "{:d} of {:d} unmasked Dörfler selections\n",
                 labels[region], counts[region],
                 total_squared > 0.0 ? 100.0 * squared[region] / total_squared : 0.0,
                 maximum[region], selected[region], total_selected);
    }
    if (!have_enriched)
    {
      Mpi::Print("  (enriched subset unavailable; regions 0 and 1 are reported together as "
                 "the protected closure)\n");
    }
    Mpi::Print("  Unmasked Dörfler would capture {:.2f}% of the error\n",
               100 * marked_error);
  };

  // Perform initial solve and estimation.
  auto [indicators, ntdof] = Solve(mesh);
  const auto ProtectIndicators =
      [this, &mesh, &ReportIndicatorRegions, &refinement](ErrorIndicator &indicators)
  {
    bool closure_conforming = true;
    auto protection = GetRefinementProtection(mesh.back()->Get(), &closure_conforming);
    ReportIndicatorRegions(indicators, protection, refinement.update_fraction);
    MFEM_VERIFY(protection.Size() == 0 || protection.Size() == mesh.back()->GetNE(),
                "AMR refinement-protection marker has an invalid size!");
    MFEM_VERIFY(closure_conforming,
                "The singular refinement closure is nonconforming before AMR!");
    const bool legacy_2d_nonconforming =
        protection.Size() > 0 && iodata.problem.type == ProblemType::ELECTROSTATIC &&
        mesh.back()->Get().Dimension() == 2 && mesh.back()->Get().Nonconforming() &&
        !refinement.singular_repair;
    if (protection.Size() > 0 && mesh.back()->Get().Nonconforming() &&
        refinement.singular_repair)
    {
      // Localized conforming repair keeps every singular trace conforming after the fact,
      // so the enriched closure participates in refinement like any other element and its
      // indicators must not be masked.
      protection.SetSize(0);
    }
    else if (protection.Size() > 0 && mesh.back()->Get().Nonconforming() &&
             !legacy_2d_nonconforming)
    {
      // Three-dimensional custom singular spaces do not yet have coarse/fine trace
      // constraints. Keep the enriched patch and its one-face buffer fixed in h.
      indicators.ZeroElements(protection);
    }
    else if (legacy_2d_nonconforming)
    {
      // Preserve the established 2D electrostatic line-tip AMR behavior: do not erase the
      // dominant near-conductor indicators. If any selected seed intersects this closure,
      // marking expands to the complete closure below and repairs its coarse side after
      // refinement.
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

    const auto marked_elements =
        [this, &comm, &refinement, &protection, &mesh](const auto &indicators)
    {
      const auto [threshold, marked_error] = utils::ComputeDorflerThreshold(
          comm, indicators.Local(), refinement.update_fraction);
      auto marked_elements = MarkedElements(indicators.Local(), threshold);
      const bool legacy_2d_nonconforming =
          protection.Size() > 0 && iodata.problem.type == ProblemType::ELECTROSTATIC &&
          mesh.back()->Get().Dimension() == 2 && mesh.back()->Get().Nonconforming() &&
          !refinement.singular_repair;
      if (legacy_2d_nonconforming)
      {
        bool intersects_closure = false;
        for (int element : marked_elements)
        {
          intersects_closure = intersects_closure || protection[element] != 0;
        }
        Mpi::GlobalOr(1, &intersects_closure, comm);
        if (intersects_closure)
        {
          mfem::Array<int> expanded;
          expanded.Reserve(marked_elements.Size() + protection.Size());
          mfem::Array<int> present(protection.Size());
          present = 0;
          for (int element : marked_elements)
          {
            expanded.Append(element);
            present[element] = 1;
          }
          for (int element = 0; element < protection.Size(); element++)
          {
            if (protection[element] && !present[element])
            {
              expanded.Append(element);
            }
          }
          marked_elements = std::move(expanded);
          Mpi::Print(" Expanded 2D nonconforming marks to the complete singular "
                     "refinement closure\n");
        }
      }
      else
      {
        for (int element : marked_elements)
        {
          MFEM_VERIFY(protection.Size() == 0 || !protection[element],
                      "AMR marked an element in a protected singular refinement closure!");
        }
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

    ReportStableSeedIds(comm, mesh.back()->Get(), marked_elements,
                        GetSourceEntityMetadata());
    ReportTraceComponents(mesh.back()->Get(), marked_elements);

    // Squared-indicator bookkeeping for the seed set, so the closure's effective marking
    // fraction can be reported after refinement. The Dörfler guarantee applies to the seed;
    // conforming closure refines a superset, and the interesting question is how much
    // additional error that superset actually captures versus how many elements it costs.
    double seed_squared = 0.0, total_squared = 0.0;
    {
      const auto &local = indicators.Local();
      for (int element = 0; element < local.Size(); element++)
      {
        total_squared += local[element] * local[element];
      }
      for (int element : marked_elements)
      {
        seed_squared += local[element] * local[element];
      }
      Mpi::GlobalSum(1, &seed_squared, comm);
      Mpi::GlobalSum(1, &total_squared, comm);
    }

    // Growth preflight. Conforming closure size is not predictable analytically, so measure
    // it on a throwaway copy of the mesh and stop the loop cleanly BEFORE the real mesh is
    // mutated or the oversized spaces are allocated. Only runs when a budget is configured,
    // because the copy costs a mesh duplication.
    if (refinement.max_growth_factor > 0.0)
    {
      const auto before = mesh.back()->Get().GetGlobalNE();
      double growth = 1.0;
      {
        mfem::ParMesh trial(mesh.back()->Get());
        trial.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
        growth = before > 0 ? static_cast<double>(trial.GetGlobalNE()) /
                                  static_cast<double>(before)
                            : 1.0;
      }
      if (growth > refinement.max_growth_factor)
      {
        Mpi::Warning(
            comm,
            "Stopping adaptive refinement: this iteration would grow the mesh by "
            "{:.3f}x, exceeding Model.Refinement.MaxGrowthFactor of {:.3f}. "
            "Conforming closure propagates beyond the marked set; inspect estimator "
            "locality/closure or explicitly raise the budget.\n",
            growth, refinement.max_growth_factor);
        break;
      }
      Mpi::Print(" Growth preflight: {:.3f}x (budget {:.3f}x)\n", growth,
                 refinement.max_growth_factor);
    }

    // Refine. Capture exact refinement ancestry while the coarse mesh is still intact.
    const auto elements_before_primary = mesh.back()->Get().GetGlobalNE();
    long long primary_added = 0;
    {
      mfem::ParMesh &fine_mesh = *mesh.back();
      ObserveRefinementAncestry(fine_mesh);
      const auto initial_elem_count = fine_mesh.GetGlobalNE();

      // Count the coarse parents the refinement actually touched, and the squared indicator
      // they carry. For conforming refinement the closure is a superset of the seed, so
      // this separates "elements added" (a child count) from "parents refined" (the true
      // closure size) and gives the effective marking fraction theta_closure.
      long long refined_parents = 0;
      double closure_squared = 0.0;
      {
        const auto &local = indicators.Local();
        std::vector<char> seeded(local.Size(), 0);
        for (int element : marked_elements)
        {
          if (element >= 0 && element < local.Size())
          {
            seeded[element] = 1;
          }
        }
        fine_mesh.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
        const auto &transforms = fine_mesh.GetRefinementTransforms();
        std::vector<char> touched(local.Size(), 0);
        for (int child = 0; child < transforms.embeddings.Size(); child++)
        {
          const int parent = transforms.embeddings[child].parent;
          if (parent >= 0 && parent < local.Size())
          {
            touched[parent] = 1;
          }
        }
        // A parent is "refined" when it has more than one child, or when it has one child
        // that is not the identity embedding. Counting children per parent is the robust
        // test and avoids relying on point-matrix identity detection.
        std::vector<int> children(local.Size(), 0);
        for (int child = 0; child < transforms.embeddings.Size(); child++)
        {
          const int parent = transforms.embeddings[child].parent;
          if (parent >= 0 && parent < local.Size())
          {
            children[parent]++;
          }
        }
        for (int element = 0; element < local.Size(); element++)
        {
          if (children[element] > 1)
          {
            refined_parents++;
            closure_squared += local[element] * local[element];
          }
        }
        Mpi::GlobalSum(1, &refined_parents, comm);
        Mpi::GlobalSum(1, &closure_squared, comm);
      }
      const auto final_elem_count = fine_mesh.GetGlobalNE();
      const long long seed_count = global_marked_elements;
      Mpi::Print(
          " Closure: {:d} seed -> {:d} refined parents (closure ratio {:.1f}), seed error "
          "{:.2f}% -> closure error {:.2f}% (θ = {:.2f}, θ_closure = {:.3f}, error "
          "amplification {:.2f})\n",
          seed_count, refined_parents,
          seed_count > 0 ? static_cast<double>(refined_parents) / seed_count : 0.0,
          total_squared > 0.0 ? 100.0 * seed_squared / total_squared : 0.0,
          total_squared > 0.0 ? 100.0 * closure_squared / total_squared : 0.0,
          refinement.update_fraction,
          total_squared > 0.0 ? closure_squared / total_squared : 0.0,
          seed_squared > 0.0 ? closure_squared / seed_squared : 0.0);
      primary_added = static_cast<long long>(final_elem_count - initial_elem_count);
      // Realized growth, for the record. The budget was already enforced by the preflight
      // above, so this is reporting only.
      const double growth = initial_elem_count > 0
                                ? static_cast<double>(final_elem_count) /
                                      static_cast<double>(initial_elem_count)
                                : 1.0;
      Mpi::Print(" Mesh growth this iteration: {:.3f}x\n", growth);
      Mpi::Print(" {} mesh refinement added {:d} elements (initial = {:d}, final = {:d})\n",
                 fine_mesh.Nonconforming() ? "Nonconforming" : "Conforming",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }
    ProcessRefinedMesh(mesh.back()->Get());

    // Localized conforming repair of the singular subcomplex. Every nonconforming face
    // carrying a singular trace is repaired by refining the coarse side, then identities
    // and features are rebuilt and detection repeats. Rebalancing deliberately happens
    // only after the subcomplex is conforming.
    long long repair_added = 0;
    {
      bool conforming = true;
      mfem::Array<int> repair;
      GetRefinementProtection(mesh.back()->Get(), &conforming, &repair);
      const bool legacy_2d_repair = iodata.problem.type == ProblemType::ELECTROSTATIC &&
                                    mesh.back()->Get().Dimension() == 2 &&
                                    refinement.nonconformal && !refinement.singular_repair;
      const bool repair_requested = refinement.singular_repair || legacy_2d_repair;
      const bool use_repair = repair_requested && !conforming;
      int pass = 0;
      // Detect a self-sustaining repair cycle: refining a coarse master can create new
      // hanging singular faces, so the violation count can plateau instead of decaying.
      // Stall detection reports that immediately rather than exhausting the pass budget.
      long long best_marks = -1;
      int passes_since_improvement = 0;
      while (repair_requested && !conforming)
      {
        const int maximum_passes =
            legacy_2d_repair ? 16 : refinement.singular_repair_max_passes;
        MFEM_VERIFY(pass < maximum_passes,
                    "Localized singular repair did not converge within "
                        << maximum_passes
                        << " passes; the singular subcomplex is still nonconforming!");
        mfem::Array<int> repair_marks;
        for (int element = 0; element < repair.Size(); element++)
        {
          if (repair[element])
          {
            repair_marks.Append(element);
          }
        }
        long long global_repair_marks = repair_marks.Size();
        Mpi::GlobalSum(1, &global_repair_marks, comm);
        // Fail closed: violations remain but no rank can repair anything.
        MFEM_VERIFY(global_repair_marks > 0,
                    "The singular subcomplex is nonconforming but no repair candidate "
                    "exists on any rank!");
        // Track the violation trajectory so a plateau is visible in the log. Refining a
        // coarse master can create new hanging singular faces, so the count may decay
        // geometrically at first and then asymptote. The pass limit and the amplification
        // and growth budgets below are the go/no-go guards; no progress heuristic is
        // applied here because a slowly-improving plateau is indistinguishable from
        // genuine convergence over a short window.
        if (best_marks < 0 || global_repair_marks < best_marks)
        {
          best_marks = global_repair_marks;
          passes_since_improvement = 0;
        }
        else
        {
          passes_since_improvement++;
        }

        mfem::ParMesh &fine_mesh = *mesh.back();
        ObserveRefinementAncestry(fine_mesh);
        const auto before = fine_mesh.GetGlobalNE();
        fine_mesh.GeneralRefinement(repair_marks, -1, refinement.max_nc_levels);
        const auto after = fine_mesh.GetGlobalNE();
        repair_added += static_cast<long long>(after - before);
        ProcessRefinedMesh(fine_mesh);
        Mpi::Print(" Singular repair pass {:d}: {:d} coarse elements marked (best {:d}, "
                   "{:d} passes without improvement), {:d} elements added (total = {:d})\n",
                   pass + 1, global_repair_marks, best_marks, passes_since_improvement,
                   after - before, after);

        // Respect the DOF/size ceiling before continuing to grow the mesh.
        MFEM_VERIFY(refinement.max_size <= 0 || after <= refinement.max_size,
                    "Localized singular repair exceeded Model.Refinement.MaxSize!");
        GetRefinementProtection(mesh.back()->Get(), &conforming, &repair);
        pass++;
      }
      if (use_repair)
      {
        const double amplification =
            static_cast<double>(repair_added) / std::max(primary_added, 1LL);
        const double growth = elements_before_primary > 0
                                  ? static_cast<double>(repair_added) /
                                        static_cast<double>(elements_before_primary)
                                  : 0.0;
        Mpi::Print(" Singular repair summary: {:d} passes, {:d} elements added, "
                   "amplification {:.2f}, global growth {:.2f}%\n",
                   pass, repair_added, amplification, 100 * growth);
        if (amplification > 3.0)
        {
          Mpi::Warning(comm,
                       "Localized singular repair amplification {:.2f} exceeds the "
                       "advisory limit of 3.\n",
                       amplification);
        }
        if (refinement.singular_repair)
        {
          MFEM_VERIFY(refinement.singular_repair_max_amplification <= 0.0 ||
                          amplification <= refinement.singular_repair_max_amplification,
                      "Localized singular repair amplification "
                          << amplification << " exceeds the configured budget of "
                          << refinement.singular_repair_max_amplification << "!");
          MFEM_VERIFY(refinement.singular_repair_max_growth <= 0.0 ||
                          growth <= refinement.singular_repair_max_growth,
                      "Localized singular repair grew the mesh by "
                          << 100 * growth << "%, exceeding the configured budget of "
                          << 100 * refinement.singular_repair_max_growth << "%!");
        }
      }
    }

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
