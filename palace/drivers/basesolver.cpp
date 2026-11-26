// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <fmt/ranges.h>
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
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "models/surfaceresponsefilmrule.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "models/surfaceresponseresolution.hpp"
#include "utils/communication.hpp"
#include "utils/dorfler.hpp"
#include "utils/edgedistance.hpp"
#include "utils/enum_string.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/metadatafile.hpp"
#include "utils/outputdir.hpp"
#include "utils/prefastprofile.hpp"
#include "utils/timer.hpp"

namespace palace
{

using json = nlohmann::json;

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
          else if (fname.ends_with("_resolved.json") || fname == kSurfaceParticipationFile)
          {
            // CUT-PLAN step 4 (lane v1table): the per-iteration table holds every level
            // and is appended at each one: it stays at the top level too.
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

void WriteMetadata(const fs::path &post_dir, const json &meta,
                   const JsonArraySpill *spill = nullptr)
{
  // Replaced through a temporary file, so a concurrent reader (e.g. an external monitor
  // polling palace.json mid-run) never observes a partial file. Round 2a R1 (decision 1475
  // (1); OI-METADATA-WRITE-UNCHECKED-ENOSPC): streamed (the bytes of dump(2) + '\n'
  // without the whole text as a string) and checked: a short write (a full disk, an I/O
  // error) fails the stream, and WriteOutputFile aborts instead of leaving an incomplete
  // palace.json behind.
  WriteOutputFile(post_dir / "palace.json",
                  [&](std::ostream &fo)
                  {
                    if (!StreamJsonMetadata(fo, meta, spill))
                    {
                      fo.setstate(std::ios::failbit);
                    }
                  });
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

// Lane diagdir (decision 1685): the ConductorConsistency Records written to `path` as the
// JSON document {"SurfaceResponse": {"Diagnostics": {"ConductorConsistency": {"Records":
// [...]}}}} (dump(2) + '\n'; the spilled entries spliced in as in palace.json, so the
// array's text is the one palace.json holds by default), checked like palace.json.
void WriteConductorConsistencyRecords(const fs::path &path, const json &records,
                                      const JsonArraySpill *spill)
{
  // A short write fails the stream, and WriteOutputFile aborts (as for palace.json).
  WriteOutputFile(path,
                  [&](std::ostream &fo)
                  {
                    if (!StreamJsonMetadata(fo, records, spill))
                    {
                      fo.setstate(std::ios::failbit);
                    }
                  });
}

// The ConductorConsistency object of palace.json's SurfaceResponse.Diagnostics when it
// holds a Records member, else null (nothing is created in `meta`).
json *FindConductorConsistencyWithRecords(json &meta)
{
  const auto response = meta.find("SurfaceResponse");
  if (response == meta.end() || !response->is_object())
  {
    return nullptr;
  }
  const auto diagnostics = response->find("Diagnostics");
  if (diagnostics == response->end() || !diagnostics->is_object())
  {
    return nullptr;
  }
  const auto gate = diagnostics->find("ConductorConsistency");
  if (gate == diagnostics->end() || !gate->is_object() || !gate->contains("Records"))
  {
    return nullptr;
  }
  return &*gate;
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

bool HasEdgeRefinement(const config::BoundaryPostData &postpro)
{
  return std::any_of(postpro.dielectric.begin(), postpro.dielectric.end(),
                     [](const auto &entry)
                     { return entry.second.edge_refinement.has_value(); });
}

void RefineEdgeTubes(const IoData &iodata, Mesh &mesh)
{
  if (!HasEdgeRefinement(iodata.boundaries.postpro))
  {
    return;
  }

  auto &par_mesh = mesh.Get();
  MPI_Comm comm = par_mesh.GetComm();
  constexpr int max_passes = 32;
  int pass = 0;
  // Edge-distance trees contain only immutable physical coordinates and remain valid as
  // h-refinement changes the volume mesh. Reuse them across all bounded tube-refinement
  // passes instead of repeating global metal-edge extraction on every pass.
  const auto contexts = BuildEdgeRefinementContexts(par_mesh, iodata.boundaries);
  while (true)
  {
    const auto marked = MarkEdgeRefinementElements(par_mesh, contexts);
    long long int global_marked = marked.Size();
    Mpi::GlobalSum(1, &global_marked, comm);
    if (global_marked == 0)
    {
      break;
    }
    MFEM_VERIFY(pass++ < max_passes,
                "Geometry-driven edge refinement failed to reach its target after "
                    << max_passes << " passes!");

    const auto initial_elements = par_mesh.GetGlobalNE();
    par_mesh.GeneralRefinement(marked, -1, iodata.model.refinement.max_nc_levels);
    mesh.Update();
    const auto final_elements = par_mesh.GetGlobalNE();
    Mpi::Print(" Edge-tube refinement pass {:d}: marked {:d} elements and added {:d} "
               "(initial = {:d}, final = {:d})\n",
               pass, global_marked, final_elements - initial_elements, initial_elements,
               final_elements);
  }
  if (pass > 0)
  {
    mesh::RebalanceMesh(iodata, mesh);
    mesh.Update();
    mesh::PrintMeshInfo(par_mesh, iodata, /*full=*/false, "Edge-refined ");
  }
}

Vector BuildMarkingIndicators(const ErrorIndicator &indicators, const mfem::ParMesh &mesh,
                              const config::BoundaryData &boundaries)
{
  Vector marking(indicators.Local());
  const auto contexts = BuildEdgeRefinementContexts(mesh, boundaries);
  if (!contexts.empty())
  {
    long long int weighted = WeightEdgeCoreIndicators(mesh, contexts, marking);
    Mpi::GlobalSum(1, &weighted, mesh.GetComm());
    if (weighted > 0)
    {
      Mpi::Print(" Downweighted {:d} elements wholly inside corrected edge cores\n",
                 weighted);
    }
  }
  return marking;
}

// The conductor-consistency end banners (idfix R0-a, DESIGN A.1 (c) / B.3 (e); USER
// decision 1000 (b)): the LOW-confidence banner whenever the gate excluded a patch on
// either path (RecordRaw: its raw claim served RAW as the ConductorConsistency share;
// the raw share per Type is the energy table's); a RecordOnly run is named a DIAGNOSTIC,
// never a run of record. Formed from the operator's statistics BEFORE SaveMetadata moves
// the diagnostics record into palace.json on the root (decision 2231, lane ccpx:
// formed after it, they read a moved-from record and never printed).
std::string ConductorConsistencyBanner(const json &statistics)
{
  std::string banner;
  const auto diagnostics = statistics.find("Diagnostics");
  if (diagnostics == statistics.end() || !diagnostics->is_object() ||
      !diagnostics->contains("ConductorConsistency") ||
      !diagnostics->at("ConductorConsistency").is_object())
  {
    return banner;
  }
  const auto &gate = diagnostics->at("ConductorConsistency");
  const int excluded = gate.value("Count", 0);
  if (excluded > 0 && gate.value("Exclusion", "") == "PerExcitation")
  {
    // Lane ccpx (decisions 2226 / 2229): each exclusion acts at its own excitation.
    int verdicts = 0, excitations = 0;
    for (const auto &entry : gate.contains("ExcludedByExcitation")
                                 ? gate.at("ExcludedByExcitation")
                                 : json::array())
    {
      verdicts += entry.value("Count", 0);
      excitations++;
    }
    banner += fmt::format(
        "\n*** LOW CONFIDENCE: the conductor-consistency gate (the device potential on "
        "each coupon's metal must match its conductor) EXCLUDED {:d} (patch, "
        "excitation) verdict(s) on {:d} distinct patch(es) at {:d} excitation(s), each "
        "at its own excitation only (Exclusion PerExcitation; {:d} patch tests; claims "
        "{:.6e}, cells {:.6e} mesh-file units of the distinct patches, served RAW as "
        "the ConductorConsistency share where excluded): see "
        "SurfaceResponse.Diagnostics.ConductorConsistency.ExcludedByExcitation of "
        "palace.json and the raw share per Type of each excitation in "
        "surface-response-conductor-consistency-energy.csv ***\n",
        verdicts, excluded, excitations, gate.value("TestedPatches", 0),
        gate.value("ClaimLength", 0.0), gate.value("CellLength", 0.0));
  }
  else if (excluded > 0)
  {
    banner += fmt::format(
        "\n*** LOW CONFIDENCE: the conductor-consistency gate (the device potential on "
        "each coupon's metal must match its conductor) EXCLUDED {:d} patch(es) of {:d} "
        "tested (claims {:.6e}, cells {:.6e} "
        "mesh-file units served RAW as the ConductorConsistency share): see "
        "SurfaceResponse.Diagnostics.ConductorConsistency of palace.json and the "
        "raw share per Type in surface-response-conductor-consistency-energy.csv "
        "***\n",
        excluded, gate.value("TestedPatches", 0), gate.value("ClaimLength", 0.0),
        gate.value("CellLength", 0.0));
  }
  if (gate.value("RecordOnly", false))
  {
    banner += fmt::format(
        "\n*** DIAGNOSTIC: SurfaceResponseCorrection.ConductorConsistency.RecordOnly "
        "was set: the gate recorded {:d} verdict(s) and EXCLUDED NOTHING; this is "
        "NOT a run of record ***\n",
        gate.contains("RecordOnlyVerdicts")
            ? static_cast<int>(gate.at("RecordOnlyVerdicts").size())
            : 0);
  }
  return banner;
}

// CUT-PLAN step 4 (lane v1table): the per-iteration table of a corrected electrostatic,
// eigenmode or driven run, with the verdict data of its Resolution.BarsFile (the built-in
// data, palace/models/surfaceparticipationverdict.json, without one or without a "Verdict"
// member); null otherwise.
std::unique_ptr<SurfaceParticipationTable>
NewSurfaceParticipationTable(const IoData &iodata, const fs::path &dir)
{
  auto Table = [&](const auto &request, const char *problem)
  {
    const auto &bars_file = request.resolution.bars_file;
    return std::make_unique<SurfaceParticipationTable>(
        (dir / kSurfaceParticipationFile).string(), problem, iodata.solver.order,
        bars_file.empty() ? BuiltInParticipationVerdictBars()
                          : ReadParticipationVerdictBars(bars_file));
  };
  const auto &electrostatic = iodata.solver.electrostatic.response_correction;
  const auto &maxwell = iodata.solver.surface_response_correction;
  switch (iodata.problem.type)
  {
    case ProblemType::ELECTROSTATIC:
      return electrostatic ? Table(*electrostatic, "Electrostatic") : nullptr;
    case ProblemType::EIGENMODE:
      return maxwell ? Table(*maxwell, "Eigenmode") : nullptr;
    case ProblemType::DRIVEN:
      return maxwell ? Table(*maxwell, "Driven") : nullptr;
    default:
      return nullptr;
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

bool BaseSolver::ContinueAdaptation(int it, long long int ntdof, double err) const
{
  // Stop when out of iterations, out of unknowns if a limit was set, or below the
  // tolerance.
  const auto &refinement = iodata.model.refinement;
  return it < refinement.max_it &&
         !(refinement.max_size > 0 && ntdof > refinement.max_size) && err >= refinement.tol;
}

bool BaseSolver::PrepareFirstLevelMesh(std::vector<std::unique_ptr<Mesh>> &mesh) const
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
  if (HasEdgeRefinement(iodata.boundaries.postpro) && mesh.size() > 1)
  {
    Mpi::Print("\nFlattening mesh sequence:\n Edge-tube refinement will start from the "
               "final mesh in the sequence of a priori refinements\n");
    mesh.erase(mesh.begin(), mesh.end() - 1);
  }

  // Resolve every configured correction tube before the first field solve, independent of
  // the relative field energy around each edge.
  RefineEdgeTubes(iodata, *mesh.back());
  return use_amr;
}

void BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = PrepareFirstLevelMesh(mesh);
  MPI_Comm comm = mesh.back()->GetComm();

  // CUT-PLAN step 4 (lane v1table): the per-iteration table of a corrected run.
  surface_participation = NewSurfaceParticipationTable(iodata, post_dir);

  // Perform initial solve and estimation.
  solve_level = 1;
  solve_levels = use_amr ? refinement.max_it + 1 : 1;
  auto [indicators, ntdof] = Solve(mesh);
  MFEM_VERIFY(solve_converged_,
              "Initial solve did not converge; no converged result to fall back on!");
  FinishSurfaceParticipationLevel(ntdof, comm);
  auto marking_indicators =
      BuildMarkingIndicators(indicators, mesh.back()->Get(), iodata.boundaries);
  double err = linalg::Norml2(comm, marking_indicators);

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

    // Mark. With Refinement.MinElementSize the Dörfler threshold is computed over the
    // eligible elements (h >= 2 x floor) only; without a floor the marking is unchanged.
    mfem::Array<int> marked_elements;
    if (!(refinement.min_element_size > 0.0))
    {
      marked_elements = [&comm, &refinement](const auto &marking_indicators)
      {
        const auto [threshold, marked_error] = utils::ComputeDorflerThreshold(
            comm, marking_indicators, refinement.update_fraction);
        const auto marked_elements = MarkedElements(marking_indicators, threshold);
        const auto [glob_marked_elements, glob_elements] =
            linalg::GlobalSize2(comm, marked_elements, marking_indicators);
        Mpi::Print(" Marked {:d}/{:d} elements for refinement ({:.2f}% of the error, θ = "
                   "{:.2f})\n",
                   glob_marked_elements, glob_elements, 100 * marked_error,
                   refinement.update_fraction);
        return marked_elements;
      }(marking_indicators);
    }
    else
    {
      const double min_size = refinement.min_element_size;
      mfem::ParMesh &fine_mesh = *mesh.back();
      const auto eligible = ElementsEligibleForMarking(fine_mesh, min_size);
      auto marking = MarkEligibleElements(comm, marking_indicators, eligible,
                                          refinement.update_fraction);
      const auto [glob_marked_elements, glob_elements] =
          linalg::GlobalSize2(comm, marking.marked_elements, marking_indicators);
      Mpi::Print(
          " Refinement floor MinElementSize = {:.3e} (eligible h >= {:.3e}): {:d}/{:d} "
          "elements ineligible for marking ({:.2f}% of the error)\n",
          min_size, 2.0 * min_size, marking.ineligible_elements, glob_elements,
          100 * marking.ineligible_error);
      if (!marking.markable)
      {
        Mpi::Print(" {}: stopping AMR, the last solved mesh is the final result\n",
                   (marking.ineligible_elements == glob_elements)
                       ? "No element is eligible for refinement under the floor"
                       : "Every eligible element under the floor has a zero error "
                         "indicator");
        --it;
        break;
      }
      Mpi::Print(
          " Marked {:d}/{:d} elements for refinement ({:.2f}% of the eligible error, "
          "θ = {:.2f})\n",
          glob_marked_elements, glob_elements, 100 * marking.marked_error,
          refinement.update_fraction);
      marked_elements = std::move(marking.marked_elements);
    }

    // Refine.
    const auto [initial_elem_count, final_elem_count] =
        [this, &mesh, &refinement, &marked_elements]
    {
      mfem::ParMesh &fine_mesh = *mesh.back();
      const long long int initial_elem_count = fine_mesh.GetGlobalNE();
      fine_mesh.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
      UpdateAfterMeshChange(*mesh.back());
      const long long int final_elem_count = fine_mesh.GetGlobalNE();
      Mpi::Print(" {} mesh refinement added {:d} elements (initial = {:d}, final = {:d})\n",
                 fine_mesh.Nonconforming() ? "Nonconforming" : "Conforming",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
      if (refinement.min_element_size > 0.0)
      {
        // Children of marked elements stay at or above the floor (h >= 2 x floor before
        // the refinement); the conformity closure may still create smaller elements.
        Mpi::Print(" Refinement floor: {:d} elements below MinElementSize = {:.3e} after "
                   "refinement (created by the conformity closure)\n",
                   CountElementsBelowSize(fine_mesh, refinement.min_element_size),
                   refinement.min_element_size);
      }
      return std::make_pair(initial_elem_count, final_elem_count);
    }();

    // Stop before a solve whose predicted size exceeds the maximum size: the last solved
    // mesh is the final result and its postprocessing output stays in place.
    if (ExceedsPredictedMaxSize(refinement, ntdof, initial_elem_count, final_elem_count))
    {
      const auto predicted_ntdof =
          PredictRefinedSize(ntdof, initial_elem_count, final_elem_count);
      Mpi::Print(
          " Predicted global unknowns of the refined mesh = {:d} (last solved {:d} x "
          "element growth {:d}/{:d} = {:.3f}) exceed max. size = {:d}: stopping AMR "
          "before the solve, the last solved mesh is the final result\n",
          predicted_ntdof, ntdof, final_elem_count, initial_elem_count,
          static_cast<double>(final_elem_count) / initial_elem_count, refinement.max_size);
      SaveMaxSizeStopMetadata({{"Iteration", it},
                               {"LastSolvedDegreesOfFreedom", ntdof},
                               {"PredictedDegreesOfFreedom", predicted_ntdof},
                               {"InitialMeshElements", initial_elem_count},
                               {"RefinedMeshElements", final_elem_count},
                               {"MaxSize", refinement.max_size}});
      // The in-memory mesh is left refined and is not used again: Run() only writes
      // metadata after this loop and the last solved results stay in place.
      --it;
      break;
    }

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

    // Optionally rebalance and write the adapted mesh to file.
    {
      const auto ratio_pre = mesh::RebalanceMesh(iodata, *mesh.back(), &mesh_counts);
      if (ratio_pre > refinement.maximum_imbalance)
      {
        int min_elem, max_elem;
        min_elem = max_elem = mesh.back()->GetNE();
        Mpi::GlobalMin(1, &min_elem, comm);
        Mpi::GlobalMax(1, &max_elem, comm);
        const auto ratio_post = double(max_elem) / min_elem;
        Mpi::Print(" Rebalanced mesh: Ratio {:.3f} exceeded max. allowed value {:.3f} "
                   "(new ratio = {:.3f})\n",
                   ratio_pre, refinement.maximum_imbalance, ratio_post);
      }
      mesh.back()->Update();
      UpdateAfterMeshChange(*mesh.back());
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
    solve_level = it + 1;
    std::tie(indicators, ntdof) = Solve(mesh);
    FinishSurfaceParticipationLevel(ntdof, comm);
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
    marking_indicators =
        BuildMarkingIndicators(indicators, mesh.back()->Get(), iodata.boundaries);
    err = linalg::Norml2(comm, marking_indicators);
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
  // Erratum 1784 E-3: the AMR loop has ended (or no AMR ran): the root's coverage block is
  // the final level's. Its Final flag is written by the end-of-run timer write
  // (SaveMetadata(const Timer &), driver.cpp), not by a load of its own (review mn-3).
  coverage_final = coverage_block_written;
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

void BaseSolver::SaveLinearSolverMetadata(MPI_Comm comm, long long int solves,
                                          long long int iterations) const
{
  long long int minimum[2] = {solves, iterations};
  long long int maximum[2] = {solves, iterations};
  Mpi::GlobalMin(2, minimum, comm);
  Mpi::GlobalMax(2, maximum, comm);
  MFEM_VERIFY(minimum[0] == maximum[0] && minimum[1] == maximum[1],
              "Rank-inconsistent linear solver statistics!");
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["LinearSolver"]["TotalSolves"] = minimum[0];
    meta["LinearSolver"]["TotalIts"] = minimum[1];
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
      // The impedance-metal rule's categories (Record only) are written when counted: a
      // run without the rule keeps today's key set.
      if ((i == Timer::CONSTRUCT_IMPEDANCE_RULE || i == Timer::POSTPRO_IMPEDANCE_RULE) &&
          timer.Counts(static_cast<Timer::Index>(i)) == 0)
      {
        continue;
      }
      std::string key(Timer::descriptions[i]);
      std::erase_if(key, isspace);
      meta["ElapsedTime"]["Durations"][key] = timer.Data(static_cast<Timer::Index>(i));
      meta["ElapsedTime"]["MinimumDurations"][key] = red.time.min[i];
      meta["ElapsedTime"]["MaximumDurations"][key] = red.time.max[i];
      meta["ElapsedTime"]["AverageDurations"][key] = red.time.avg[i];
      meta["ElapsedTime"]["Counts"][key] = timer.Counts(static_cast<Timer::Index>(i));
      meta["PeakMemoryGrowthMegabytes"]["Min"][key] = red.rank_mem.min[i] * to_mb;
      meta["PeakMemoryGrowthMegabytes"]["Max"][key] = red.rank_mem.max[i] * to_mb;
      meta["PeakMemoryGrowthMegabytes"]["Sum"][key] = red.rank_mem.sum[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Min"][key] = red.node_mem.min[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Max"][key] = red.node_mem.max[i] * to_mb;
      meta["PeakNodeMemoryGrowthMegabytes"]["Sum"][key] = red.node_mem.sum[i] * to_mb;
    }
    // Erratum 1784 E-3 / review mn-3: after the AMR loop the root's coverage block gets
    // "Final": true in this load (the first end-of-run write; the in-loop timer writes
    // come before the flag is set).
    if (coverage_final)
    {
      meta["SurfaceResponse"]["Coverage"]["Final"] = true;
      // Lane rescensus-product (decision 1817): the census block of the same last level.
      if (meta["SurfaceResponse"].contains("Resolution"))
      {
        meta["SurfaceResponse"]["Resolution"]["Final"] = true;
      }
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

void BaseSolver::SaveImpedanceRuleCandidateMetadata(const nlohmann::json &block) const
{
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    meta["SurfaceResponse"]["ImpedanceRuleCandidate"] = block;
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(SurfaceResponseOperator &response) const
{
  // Statistics reduction is collective and must occur before the root-only file update.
  // The final write (round 2a R1, decision 1475 (1)): the operator's diagnostics record is
  // moved, not copied (ReleaseStatistics; the record is the same).
  auto statistics = response.ReleaseStatistics();
  const std::string conductor_consistency_banner = ConductorConsistencyBanner(statistics);
  std::string coverage_banner, resolution_banner;
  std::vector<std::string> near_edge_warnings;
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    for (auto &[key, value] : statistics.items())
    {
      if (key == "Diagnostics" && value.is_object())
      {
        // Merged: a diagnostics record written before the operator's (the EdgeExclusion
        // coverage of the interface dielectrics) survives.
        for (auto &[name, entry] : value.items())
        {
          meta["SurfaceResponse"][key][name] = std::move(entry);
        }
        continue;
      }
      meta["SurfaceResponse"][key] = std::move(value);
    }
    // The v1 default path (decisions 2377 / 2379): a run with the developer switch records
    // the parked options it selects (the list of its "DEVELOPER OPTIONS" log line).
    for (const auto *config : {iodata.solver.electrostatic.response_correction
                                   ? &*iodata.solver.electrostatic.response_correction
                                   : nullptr,
                               iodata.solver.surface_response_correction
                                   ? &*iodata.solver.surface_response_correction
                                   : nullptr})
    {
      if (config && config->developer)
      {
        meta["SurfaceResponse"]["Developer"] = config->developer_options;
      }
    }
    // Lane L1 (decision 1634 (4); DESIGN-v2 s4.5): the coverage of the run's cells (the
    // writers' record, root only); an empty record (no cell measured, e.g. a raw
    // fallback) writes its block too, with Cells 0 (review m-6: never silent).
    // Ruling 1738: a raw fallback measures no cell: its reason (the run-level geometry
    // fallback classes) is named, so that zero cells are never read as complete coverage.
    if (response.IsRawFallback() && response.GetCoverageRunRecord().entries.empty())
    {
      std::string classes;
      for (const auto &fallback : response.GetGeometryFallbacks())
      {
        if (fallback.what == "run" && classes.find(fallback.cls) == std::string::npos)
        {
          classes += (classes.empty() ? "" : "+") + fallback.cls;
        }
      }
      response.GetCoverageRunRecord().no_cells_reason =
          "RawFallback" + (classes.empty() ? std::string() : ": " + classes) +
          ": no correction applied";
    }
    // Lane rescensus-product (review r1817 mn-1): the census block of a run without a cell
    // names the same reason as the coverage block.
    if (response.GetResolutionRunRecord().entries.empty())
    {
      response.GetResolutionRunRecord().no_cells_reason =
          response.GetCoverageRunRecord().no_cells_reason;
    }
    meta["SurfaceResponse"]["Coverage"] = CoverageRunMetadata(
        response.GetCoverageRunRecord(), response.GetCoverageThresholds(),
        response.GetCoverageValidatedRegime(),
        iodata.units.Dimensionalize<Units::ValueType::ENERGY>(1.0), solve_level,
        solve_levels, response.CoverageMarkNotApplied());
    coverage_block_written = true;
    coverage_banner = CoverageRunBanner(meta["SurfaceResponse"]["Coverage"]);
    // Ruling 2276 R7 (b): the near-edge warnings of the run's cells (root: the record).
    near_edge_warnings = CoverageNearEdgeWarnings(response.GetCoverageRunRecord());
    // Lane rescensus-product (decision 1817; ruling 1819): the edge-resolution census of
    // the level (written beside the coverage block; "Final" set with the coverage block's).
    meta["SurfaceResponse"]["Resolution"] =
        ResolutionRunMetadata(response.GetResolutionRunRecord(), solve_level, solve_levels);
    resolution_banner = ResolutionRunBanner(meta["SurfaceResponse"]["Resolution"]);
    // The impedance-metal rule's "Apply" (decision 2206 (a); DESIGN-APPLY s6): its record
    // beside the coverage block, and one banner line when it served a feature on an
    // impedance sheet (the tables of record keep their headers and status words).
    if (const auto *film = response.GetImpedanceApplyRule())
    {
      auto &record = meta["SurfaceResponse"]["ImpedanceMetalRule"];
      const auto &features = response.GetImpedanceLawFeatures();
      record = ImpedanceApplyRecord(*film, features, iodata.units, false);
      const auto counts = CountImpedanceApplyFeatures(features);
      if (counts.corrected > 0)
      {
        std::vector<std::string> inductances;
        for (const double inductance : film->served_inductances)
        {
          inductances.push_back(fmt::format("{:.6g}", inductance));
        }
        coverage_banner +=
            fmt::format("*** IMPEDANCE-METAL RULE (APPLY): {}; served Ls: {} H/sq ***\n",
                        ImpedanceApplyCountText(counts), fmt::join(inductances, ", "));
      }
    }
    // Round 2a R2 (decision 1475 (1)): the Records spilled by the root are spliced in place
    // of the empty in-memory array (the same bytes), then the scratch file is removed.
    auto *spill = response.GetConductorConsistencyRecordsSpill();
    if (spill)
    {
      MFEM_VERIFY(spill->Flush(),
                  "Writing the conductor-consistency records to the scratch "
                  "file \""
                      << spill->Path() << "\" failed (a short write)!");
      MFEM_VERIFY(spill->Size() == 0 || HoldsSpillPlace(meta, *spill),
                  "The metadata holds no empty SurfaceResponse.Diagnostics."
                  "ConductorConsistency.Records array for the "
                      << spill->Size() << " spilled conductor-consistency record(s)!");
    }
    json *gate = response.ConductorConsistencyRecordsInDirectory()
                     ? FindConductorConsistencyWithRecords(meta)
                     : nullptr;
    if (gate)
    {
      // Decision 1685 (Diagnostics.ConductorConsistencyRecords "Directory"): palace.json
      // without the Records member (everything else as above); the array (the spilled
      // entries or the in-memory ones) at the same key path in the records file.
      json records;
      records["SurfaceResponse"]["Diagnostics"]["ConductorConsistency"]["Records"] =
          std::move((*gate)["Records"]);
      gate->erase("Records");
      WriteMetadata(post_dir, meta);
      WriteConductorConsistencyRecords(
          response.DiagnosticsPath(
              post_dir, SurfaceResponseOperator::kConductorConsistencyRecordsFile,
              SurfaceResponseOperator::kConductorConsistencyRecordsDirectory),
          records, spill);
    }
    else
    {
      WriteMetadata(post_dir, meta, spill);
    }
    if (spill)
    {
      MFEM_VERIFY(spill->Remove(), "Removing the conductor-consistency records scratch "
                                   "file \""
                                       << spill->Path() << "\" failed!");
    }
  }
  Mpi::Print("{}", conductor_consistency_banner);
  // Geometry fallbacks (USER decision 1000 (b) / 1001): the LOW-confidence banner closes
  // the run summary whenever a converted site fired; never silent.
  if (const auto record = response.GetGeometryFallbackRecord(); !record.is_null())
  {
    Mpi::Print(
        "{}",
        GeometryFallbackBanner(
            record, "SurfaceResponse.Diagnostics.GeometryFallbacks of palace.json" +
                        std::string(response.IsRawFallback()
                                        ? " (no correction applied: the thin solve)"
                                        : " and surface-response-geometry-fallbacks.csv")));
  }
  // Lane L1 (DESIGN-v2 s4.5): the coverage banner closes the run summary (a line per
  // cell whose coverage status is not OK or whose regime is not IN; never silent).
  if (!coverage_banner.empty())
  {
    Mpi::Print("{}", coverage_banner);
  }
  for (const auto &line : near_edge_warnings)
  {
    Mpi::Warning("{}\n", line);
  }
  if (!resolution_banner.empty())
  {
    Mpi::Print("{}", resolution_banner);
  }
}

void BaseSolver::SaveMetadata(const SurfacePostOperator &surf_post_op) const
{
  // The perimeter is gathered globally by every rank (GetBoundaryEdgeSegments), so the
  // record is rank-consistent and the root writes it.
  if (!root)
  {
    return;
  }
  nlohmann::json interfaces = nlohmann::json::object();
  const double length_scale = iodata.units.GetMeshLengthRelativeScale();
  for (const auto &[idx, data] : surf_post_op.eps_surfs)
  {
    if (!data.edge_exclusion)
    {
      continue;
    }
    const auto &exclusion = *data.edge_exclusion;
    nlohmann::json entry;
    entry["ExcludeSegmentsCount"] = exclusion.exclude_segments_count;
    entry["ExcludeSegmentsLength"] = exclusion.exclude_segments_length * length_scale;
    entry["PerimeterCount"] = exclusion.perimeter_count;
    entry["PerimeterLength"] = exclusion.perimeter_length * length_scale;
    entry["ExcludedPerimeterCount"] = exclusion.excluded_perimeter_count;
    entry["ExcludedPerimeterLength"] = exclusion.excluded_perimeter_length * length_scale;
    entry["StraddlingPerimeterLength"] =
        exclusion.straddling_perimeter_length * length_scale;
    entry["RemovedRuns"] = exclusion.removed_runs;
    entry["ListedChains"] = exclusion.listed_chains;
    entry["Coverage"] =
        exclusion.exclude_segments_length > 0.0
            ? exclusion.excluded_perimeter_length / exclusion.exclude_segments_length
            : 0.0;
    entry["EndAllowanceCoverage"] = exclusion.exclude_segments_length > 0.0
                                        ? (exclusion.excluded_perimeter_length +
                                           exclusion.straddling_perimeter_length) /
                                              exclusion.exclude_segments_length
                                        : 0.0;
    entry["Tolerance"] = exclusion.tolerance * length_scale;
    entry["MaxSagitta"] = exclusion.max_sagitta * length_scale;
    interfaces[std::to_string(idx)] = std::move(entry);
  }
  if (interfaces.empty())
  {
    return;
  }
  nlohmann::json record;
  record["Interfaces"] = std::move(interfaces);
  record["Rule"] =
      "decision 657 (b): per interface dielectric with EdgeExcludeSegments, the listed "
      "segments' total length (ExcludeSegmentsLength) and the perimeter length the "
      "sagitta-aware rule-B5 coincidence test removed (ExcludedPerimeterLength), both in "
      "mesh length units; Coverage = ExcludedPerimeterLength / ExcludeSegmentsLength (the "
      "perimeter is the single outline of the edge attributes, so a fully excluded list "
      "reads 1); StraddlingPerimeterLength = the length of the KEPT perimeter segments "
      "that contain a listed chain's END point within the Tolerance, at most one kept "
      "segment per listed end (decision 667 (1); a listed end = an endpoint of a listed "
      "segment no other listed segment's endpoint meets within the Tolerance), and "
      "EndAllowanceCoverage = (ExcludedPerimeterLength + StraddlingPerimeterLength) / "
      "ExcludeSegmentsLength (decision 665: a listed end falling inside a perimeter "
      "segment keeps that segment whole, a deficit of up to one segment per listed end); "
      "RemovedRuns = the maximal runs of removed perimeter segments chained by shared "
      "vertices and ListedChains = the chains of listed segments joined end to end within "
      "the Tolerance (information: a correct exclusion reads RemovedRuns = ListedChains); "
      "Tolerance = EdgeExcludeSegmentTolerance, MaxSagitta = the largest "
      "EdgeExcludeSegmentSagittas entry";
  nlohmann::json meta = LoadMetadata(post_dir);
  meta["SurfaceResponse"]["Diagnostics"]["EdgeExclusion"] = std::move(record);
  WriteMetadata(post_dir, meta);
}

void BaseSolver::SaveSurfaceResponseSolverMetadata(
    MPI_Comm comm, const std::string &name, long long int solves, long long int iterations,
    std::optional<long long int> unconverged) const
{
  long long int minimum[3] = {solves, iterations, unconverged.value_or(0)};
  long long int maximum[3] = {solves, iterations, unconverged.value_or(0)};
  Mpi::GlobalMin(3, minimum, comm);
  Mpi::GlobalMax(3, maximum, comm);
  MFEM_VERIFY(minimum[0] == maximum[0] && minimum[1] == maximum[1] &&
                  minimum[2] == maximum[2],
              "Rank-inconsistent corrected surface-response solver statistics!");
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    meta["SurfaceResponse"]["CorrectedSolver"][name]["LinearSolves"] = minimum[0];
    meta["SurfaceResponse"]["CorrectedSolver"][name]["LinearIterations"] = minimum[1];
    if (unconverged)
    {
      meta["SurfaceResponse"]["CorrectedSolver"][name]["NonConvergedLinearSolves"] =
          minimum[2];
    }
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::FinishSurfaceParticipationLevel(long long int ntdof, MPI_Comm comm) const
{
  if (!surface_participation)
  {
    return;
  }
  surface_participation->FinishLevel(solve_level, ntdof, solve_converged_, comm);
  if (surface_participation->Started())
  {
    SaveSurfaceResponseMetadata("ParticipationTable", surface_participation->Metadata());
  }
}

void BaseSolver::SaveSurfaceResponseMetadata(const std::string &key,
                                             const nlohmann::json &value) const
{
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    meta["SurfaceResponse"][key] = value;
    WriteMetadata(post_dir, meta);
  }
}

PreflightCensusOutcome BaseSolver::PreflightResolutionCensus(
    std::vector<std::unique_ptr<Mesh>> &mesh, bool /*fail_closed*/,
    const std::shared_ptr<const SurfaceResponseGeometry> & /*identified*/) const
{
  ReportResolutionCensusNotEvaluated(
      mesh.back()->GetComm(),
      "Resolution census: not evaluated in PreflightOnly for this problem type "
      "(PreflightOnly evaluates it for Electrostatic, Eigenmode and Driven)");
  return {};
}

void BaseSolver::ReportResolutionCensusNotEvaluated(MPI_Comm comm,
                                                    const std::string &line) const
{
  Mpi::Print(comm, "\n{}\n", line);
  SaveSurfaceResponseMetadata("ResolutionCensusNotEvaluated", line);
}

void BaseSolver::PrepareCensusFirstLevelMesh(
    std::vector<std::unique_ptr<Mesh>> &mesh,
    const std::shared_ptr<const SurfaceResponseGeometry> &identified,
    std::shared_ptr<const SurfaceResponseGeometry> &response_geometry) const
{
  const Mesh *identified_mesh = mesh.back().get();
  const long int identified_sequence = identified_mesh->Get().GetSequence();
  PrepareFirstLevelMesh(mesh);
  if (identified && mesh.back().get() == identified_mesh &&
      mesh.back()->Get().GetSequence() == identified_sequence)
  {
    response_geometry = identified;
  }
}

bool BaseSolver::PreflightCensusStopsFailClosed(const SurfaceResponseOperator &response,
                                                bool fail_closed, MPI_Comm comm) const
{
  // The ledger is replicated on every rank (every rank decides alike).
  const auto &fallbacks = response.GetGeometryFallbacks();
  if (!fail_closed || fallbacks.empty())
  {
    return false;
  }
  ReportResolutionCensusNotEvaluated(
      comm,
      fmt::format(
          "Resolution census: not evaluated: GeometryFailurePolicy \"FailClosed\" with "
          "{:d} geometry fallback(s) recorded by the response operator's construction "
          "(a solve-form run stops there); \"RecordRaw\" gives it",
          fallbacks.size()));
  return true;
}

PreflightCensusOutcome BaseSolver::PreflightMaxwellResolutionCensus(
    const std::vector<std::unique_ptr<Mesh>> &mesh,
    std::shared_ptr<const SurfaceResponseGeometry> *response_geometry,
    bool fail_closed) const
{
  BlockTimer bt(Timer::CONSTRUCT);
  prefast::Mark("census-spaceop-B");
  SpaceOperator space_op(iodata, mesh);
  prefast::Mark("census-spaceop-E");
  auto response =
      std::make_unique<SurfaceResponseOperator>(iodata, space_op, response_geometry);
  prefast::Mark("census-response-E");
  response->SetPatchRecordResidency(
      SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY);
  if (PreflightCensusStopsFailClosed(*response, fail_closed, space_op.GetComm()))
  {
    return {};
  }
  // The PostOperator's check of a corrected Maxwell run (before its first solve): a solve
  // form stops there (ruling R1, P7).
  if (!response->IsRawFallback() && !response->HasSurfaceResponse())
  {
    ReportResolutionCensusNotEvaluated(
        space_op.GetComm(),
        "Resolution census: not evaluated: the selected coupon models carry no "
        "fabricated and thin surface response matrices (a solve-form run stops there)");
    return {};
  }
  return WritePreflightResolutionCensus(*response, space_op.GlobalTrueVSize(),
                                        space_op.GetComm());
}

PreflightCensusOutcome
BaseSolver::WritePreflightResolutionCensus(SurfaceResponseOperator &response,
                                           long long int unknowns, MPI_Comm comm) const
{
  if (response.IsRawFallback())
  {
    // The solve form makes no census on a whole-run raw fallback either.
    ReportResolutionCensusNotEvaluated(
        comm, "Resolution census: not evaluated: the correction falls back to the raw run "
              "(RawFallback: no patch applied)");
    return {};
  }
  // The census Types of the solve form (PostOperator::GetInterfaceTypes over the coverage
  // interfaces: each entry's Type, read here from the configuration it copies).
  const auto &dielectric = iodata.boundaries.postpro.dielectric;
  std::map<int, std::string> types;
  std::set<int> regions;
  for (const int interface : response.GetCoverageInterfaces())
  {
    const auto entry = dielectric.find(interface);
    if (entry != dielectric.end())
    {
      types[interface] = ToString(entry->second.type);
      if (entry->second.region)
      {
        regions.insert(interface);
      }
    }
  }
  response.ConfigureResolutionTypes(types);
  const auto &request = iodata.problem.type == ProblemType::ELECTROSTATIC
                            ? iodata.solver.electrostatic.response_correction
                            : iodata.solver.surface_response_correction;
  if (request && request->patch_energy)
  {
    // The solve form's per-patch census table (PatchEnergy only; collective).
    const auto path = response.CollectiveDiagnosticsPath(
        post_dir, "surface-response-resolution-patches.csv");
    response.WriteResolutionPatchTable(path.string());
  }
  auto &record = response.GetResolutionRunRecord();
  record.no_cells_reason = "PreflightOnly: no solve";
  auto block = ResolutionRunMetadata(record, 1, 1);
  block["Final"] = true;
  SaveSurfaceResponseMetadata("Resolution", block);
  // CUT-PLAN step 4b with step 3 (Ls sheets served by default): the impedance-metal rule's
  // record of the solve form (its mode, regime and served-Ls lines), which that run writes
  // beside its first coverage block from this operator alone.
  if (const auto *film = response.GetImpedanceApplyRule())
  {
    SaveSurfaceResponseMetadata("ImpedanceMetalRule",
                                ImpedanceApplyRecord(*film,
                                                     response.GetImpedanceLawFeatures(),
                                                     iodata.units, false));
  }
  surface_participation = NewSurfaceParticipationTable(iodata, post_dir);
  if (surface_participation)
  {
    surface_participation->WritePreflightLevel(
        PreflightParticipationRows(types, regions, record.preflight), unknowns, comm);
    SaveSurfaceResponseMetadata("ParticipationTable", surface_participation->Metadata());
  }
  PreflightCensusOutcome outcome;
  outcome.evaluated = true;
  std::set<int> fallback_features;
  for (const auto &fallback : response.GetGeometryFallbacks())
  {
    fallback_features.insert(fallback.features.begin(), fallback.features.end());
  }
  outcome.fallback_features = fallback_features.size();
  return outcome;
}

void BaseSolver::SaveMaxSizeStopMetadata(const nlohmann::json &record) const
{
  if (root)
  {
    nlohmann::json meta = LoadMetadata(post_dir);
    meta["AdaptiveMeshRefinement"]["MaxSizeStop"] = record;
    WriteMetadata(post_dir, meta);
  }
}

long long int PredictRefinedSize(long long int ntdof, long long int initial_elem_count,
                                 long long int final_elem_count)
{
  MFEM_VERIFY(ntdof >= 0 && initial_elem_count > 0 &&
                  final_elem_count >= initial_elem_count,
              "Invalid sizes for the predicted size of a refined mesh (ntdof = "
                  << ntdof << ", elements " << initial_elem_count << " -> "
                  << final_elem_count << ")!");
  const double predicted = static_cast<double>(ntdof) *
                           static_cast<double>(final_elem_count) /
                           static_cast<double>(initial_elem_count);
  return static_cast<long long int>(std::ceil(predicted));
}

bool ExceedsPredictedMaxSize(const config::RefinementData &refinement, long long int ntdof,
                             long long int initial_elem_count,
                             long long int final_elem_count)
{
  if (!refinement.max_size_predicted || refinement.max_size <= 0)
  {
    return false;
  }
  return PredictRefinedSize(ntdof, initial_elem_count, final_elem_count) >
         refinement.max_size;
}

std::vector<bool> ElementsEligibleForMarking(mfem::Mesh &mesh, double min_element_size)
{
  // GetElementSize(i, 0) is the volume-based size pow(|det J|, 1/dim) at the element center
  // that mfem::ParMesh::PrintInfo reports as h in the mesh statistics.
  std::vector<bool> eligible(mesh.GetNE());
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    eligible[i] = (mesh.GetElementSize(i, 0) >= 2.0 * min_element_size);
  }
  return eligible;
}

long long int CountElementsBelowSize(mfem::ParMesh &mesh, double min_element_size)
{
  long long int count = 0;
  for (int i = 0; i < mesh.GetNE(); i++)
  {
    count += (mesh.GetElementSize(i, 0) < min_element_size);
  }
  Mpi::GlobalSum(1, &count, mesh.GetComm());
  return count;
}

FloorMarking MarkEligibleElements(MPI_Comm comm, const Vector &e,
                                  const std::vector<bool> &eligible, double fraction)
{
  MFEM_VERIFY(static_cast<int>(eligible.size()) == e.Size(),
              "Marking eligibility must have one entry per local element!");
  FloorMarking marking;
  e.HostRead();
  Vector eligible_e(e.Size());
  int num_eligible = 0;
  double ineligible_error = 0.0, eligible_error = 0.0, total_error = 0.0;
  for (int i = 0; i < e.Size(); i++)
  {
    total_error += e[i] * e[i];
    if (eligible[i])
    {
      eligible_e[num_eligible++] = e[i];
      eligible_error += e[i] * e[i];
    }
    else
    {
      ineligible_error += e[i] * e[i];
    }
  }
  eligible_e.SetSize(num_eligible);
  long long int glob_eligible = num_eligible;
  Mpi::GlobalSum(1, &glob_eligible, comm);
  std::array<double, 3> error = {ineligible_error, eligible_error, total_error};
  Mpi::GlobalSum(3, error.data(), comm);
  marking.ineligible_elements = e.Size() - num_eligible;
  Mpi::GlobalSum(1, &marking.ineligible_elements, comm);
  marking.ineligible_error = (error[2] > 0.0) ? error[0] / error[2] : 0.0;
  // An eligible set whose indicators are all zero has no Dörfler threshold (0 / 0): it is
  // treated like the empty eligible set, nothing marked.
  marking.markable = (glob_eligible > 0 && error[1] > 0.0);
  if (!marking.markable)
  {
    return marking;
  }
  const auto [threshold, marked_error] =
      utils::ComputeDorflerThreshold(comm, eligible_e, fraction);
  marking.threshold = threshold;
  marking.marked_error = marked_error;
  marking.marked_elements.Reserve(num_eligible);
  for (int i = 0; i < e.Size(); i++)
  {
    if (eligible[i] && e[i] >= marking.threshold)
    {
      marking.marked_elements.Append(i);
    }
  }
  return marking;
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
