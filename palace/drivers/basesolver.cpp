// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <algorithm>
#include <complex>
#include <numeric>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/domainpostoperator.hpp"
#include "fem/postoperator.hpp"
#include "fem/surfacepostoperator.hpp"
#include "linalg/errorestimator.hpp"
#include "utils/communication.hpp"
#include "utils/errorindicators.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using json = nlohmann::json;

namespace
{

std::string GetPostDir(const std::string &output)
{
  return (output.length() > 0 && output.back() != '/') ? output + '/' : output;
}

json LoadMetadata(const std::string &post_dir)
{
  std::string path = post_dir + "palace.json";
  std::ifstream fi(path);
  if (!fi.is_open())
  {
    MFEM_ABORT("Unable to open metadata file \"" << path << "\"!");
  }
  return json::parse(fi);
}

void WriteMetadata(const std::string &post_dir, const json &meta)
{
  std::string path = post_dir + "palace.json";
  std::ofstream fo(path);
  if (!fo.is_open())
  {
    MFEM_ABORT("Unable to open metadata file \"" << path << "\"!");
  }
  fo << meta.dump(2) << '\n';
}

}  // namespace

BaseSolver::BaseSolver(const IoData &iodata_, bool root_, int size, int num_thread,
                       const char *git_tag)
  : iodata(iodata_), post_dir(GetPostDir(iodata_.problem.output)), root(root_),
    table(8, 9, 6)
{
  // Create directory for output.
  if (root && !std::filesystem::exists(post_dir))
  {
    std::filesystem::create_directories(post_dir);
  }

  // Initialize simulation metadata for this simulation.
  if (root && post_dir.length() > 0)
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

std::string BaseSolver::IterationPostDir(int iter) const
{
  std::string dir = fmt::format("{}Adapt{:0>4d}/", post_dir, iter);
  // Create directory for output.
  if (root && !std::filesystem::exists(dir))
  {
    std::filesystem::create_directories(dir);
  }
  return dir;
}

namespace
{

// Given a vector of estimates local to this rank, compute an error threshold
// such that if elements with greater error are marked for refinement, a Dorfler
// marking is achieved.
double ComputeRefineThreshold(double fraction, const mfem::Vector &e)
{
  // Pre compute the sort and partial sum to make evaluating a candidate
  // partition very fast.
  std::vector<double> estimates(e.begin(), e.end());
  std::sort(estimates.begin(), estimates.end());

  MFEM_ASSERT(estimates.front() >= 0, "Indicators must be non-negative");

  std::vector<double> sum;
  sum.reserve(estimates.size());
  std::partial_sum(estimates.begin(), estimates.end(), std::back_inserter(sum));

  auto pivot = std::distance(
      sum.begin(), std::lower_bound(sum.begin(), sum.end(), (1 - fraction) * sum.back()));

  double error_threshold = estimates[pivot];

  auto marked = [&estimates, &sum](double e) -> std::pair<std::size_t, double>
  {
    const auto lb = std::lower_bound(estimates.begin(), estimates.end(), e);
    const auto elems_marked = std::distance(lb, estimates.end());
    const double error_unmarked =
        lb != estimates.begin() ? sum[sum.size() - elems_marked - 1] : 0;
    const double error_marked = sum.back() - error_unmarked;
    return {elems_marked, error_marked};
  };

  // Each rank has computed a different threshold, if a given rank has lots of
  // low error elements, their value will be lower and if a rank has high error,
  // their value will be higher. Using the value from the low error rank will
  // give too many elements, and using the value from the high error rank will
  // give too few. The correct threshold value will be an intermediate
  // between the min and max over ranks.

  auto comm = Mpi::World();

  // Send initial information to all
  const double total_error = [&]()
  {
    double total_error = sum.back();
    Mpi::GlobalSum(1, &total_error, comm);
    return total_error;
  }();
  const double max_indicator = [&]()
  {
    double max_indicator = estimates.back();
    Mpi::GlobalSum(1, &max_indicator, comm);
    return max_indicator;
  }();

  double min_threshold = error_threshold;
  double max_threshold = error_threshold;
  Mpi::GlobalMin(1, &min_threshold, comm);
  Mpi::GlobalMax(1, &max_threshold, comm);
  const std::size_t total_elem = [&]()
  {
    std::size_t tmp = estimates.size();
    Mpi::GlobalSum(1, &tmp, comm);
    return tmp;
  }();

  MFEM_ASSERT(min_threshold <= max_threshold,
              "min: " << min_threshold << " max " << max_threshold);

  auto [elem_marked, error_marked] = marked(error_threshold);

  // Keep track of the number of elements marked by the threshold bounds. If the
  // top and bottom values are equal, there's no point further bisecting.
  auto [max_elem_marked, min_error_marked] = marked(min_threshold);
  auto [min_elem_marked, max_error_marked] = marked(max_threshold);
  Mpi::GlobalSum(1, &min_elem_marked, comm);
  Mpi::GlobalSum(1, &max_elem_marked, comm);

  constexpr int maxiter = 100;  // Maximum limit to prevent runaway.
  for (int i = 0; i < maxiter; ++i)
  {
    error_threshold = (min_threshold + max_threshold) / 2;

    std::tie(elem_marked, error_marked) = marked(error_threshold);

    // All processors need the values used for the stopping criteria
    Mpi::GlobalSum(1, &elem_marked, comm);
    Mpi::GlobalSum(1, &error_marked, comm);

    MFEM_ASSERT(elem_marked > 0, "Some elements must have been marked");
    MFEM_ASSERT(error_marked > 0, "Some error must have been marked");

    const auto candidate_fraction = error_marked / total_error;

    Mpi::Debug("Threshold: {:e} < {:e} < {:e}\n", min_threshold, error_threshold,
               max_threshold);
    Mpi::Debug("Marked Elems: {} <= {} <= {}\n", min_elem_marked, elem_marked,
               max_elem_marked);

    // set the tolerance based off of the largest local indicator value. These
    // tolerance values are chosen based on testing, opt not to expose them.
    const double frac_tol = 1e-3 * fraction;
    const double error_tol = 1e-6 * max_indicator;
    if (std::abs(max_threshold - min_threshold) < error_tol ||
        std::abs(candidate_fraction - fraction) < frac_tol ||
        max_elem_marked == min_elem_marked)
    {
      // candidate fraction matches to tolerance, or the number of marked
      // elements is no longer changing.
      Mpi::Debug(
          "ΔFraction: {:.3e}, Tol {:.3e}, ΔThreshold: {:.3e}, Tol {:.3e},  ΔElements: {}\n",
          candidate_fraction - fraction, frac_tol, max_threshold - min_threshold, error_tol,
          max_elem_marked - min_elem_marked);
      break;
    }

    // Update in preparation for next iteration. The logic here looks inverted
    // compared to a usual binary search, because a smaller value marks a larger
    // number of elements and thus a greater fraction of error.
    if (candidate_fraction > fraction)
    {
      // This candidate marked too much, raise the lower value.
      min_threshold = error_threshold;
      max_elem_marked = elem_marked;
    }
    else if (candidate_fraction < fraction)
    {
      // This candidate marked too little, lower the upper bound
      max_threshold = error_threshold;
      min_elem_marked = elem_marked;
    }
  }

  Mpi::Print("Indicator threshold {:.3e} marked {} of {} elements and {:.3f}\% error\n",
             error_threshold, elem_marked, total_elem, 100 * error_marked / total_error);

  MFEM_ASSERT(error_threshold > 0, "error_threshold must be positive");
  return error_threshold;
}

mfem::Array<int> MarkedElements(double threshold, const mfem::Vector &e, bool gt = true)
{
  mfem::Array<int> ind;
  ind.Reserve(e.Size());
  for (int i = 0; i < e.Size(); i++)
  {
    if (gt && e[i] >= threshold)  // Used for refinement marking.
    {
      ind.Append(i);
    }

    if (!gt && e[i] <= threshold)  // Used for coarsening marking.
    {
      ind.Append(i);
    }
  }
  return ind;
}

void RebalanceMesh(std::unique_ptr<mfem::ParMesh> &mesh)
{
  auto comm = Mpi::World();
  if (Mpi::Size(comm) > 1)
  {
    if (mesh->Nonconforming())
    {
      mesh->Rebalance();
    }
    else
    {
      // Without access to a refinement tree, partitioning must be done on the
      // root processor and then redistributed.
      mesh::RebalanceConformalMesh(mesh);
    }
  }
}

}  // namespace

ErrorIndicators
BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                                    Timer &timer) const
{
  // Helper to save off postprocess data.
  auto save_postprocess = [&, this](int iter)
  {
    if (Mpi::Root(Mpi::World()))
    {
      namespace fs = std::filesystem;
      // Create a subfolder denoting the results of this adaptation.
      const auto out_dir = fs::path(IterationPostDir(iter));

      const auto options =
          fs::copy_options::recursive | fs::copy_options::overwrite_existing;
      const fs::path root_dir(post_dir);
      for (const auto &f : fs::directory_iterator(root_dir))
      {
        if (f.is_regular_file())
        {
          fs::copy(f, out_dir / fs::path(f).filename(), options);
        }
      }
      if (fs::exists(root_dir / "paraview"))
      {
        fs::copy(root_dir / "paraview", out_dir / "paraview", options);
      }
    }
  };

  const auto &param = iodata.model.refinement.adaptation;
  const bool use_amr = param.max_its > 0;
  const bool use_coarsening = param.coarsening_fraction > 0;

  int iter = 0;
  auto indicators = Solve(mesh, timer);

  if (use_amr)
  {
    Mpi::Print("\nAdaptive Mesh Refinement Parameters:\n");
    Mpi::Print("MinIter: {}, MaxIter: {}, Tolerance: {:.3e}, DOFLimit: {}\n\n",
               param.min_its, param.max_its, param.tolerance, param.dof_limit);
    save_postprocess(iter);  // Save the initial solution
  }

  // Collection of all tests that might exhaust resources.
  auto exhausted_resources = [&]()
  {
    bool ret = false;
    // run out of DOFs, and coarsening isn't allowed.
    ret |= (!use_coarsening && indicators.ndof > param.dof_limit);
    ret |= iter > param.max_its;

    return ret;
  };

  while ((iter < param.min_its || indicators.global_error_indicator > param.tolerance) &&
         !exhausted_resources())
  {
    Mpi::Print("Adaptation iteration {}: Error Indicator: {:.3e}, DOF: {}\n", iter,
               indicators.global_error_indicator, indicators.ndof);
    if (indicators.ndof < param.dof_limit)
    {
      // refinement mark
      const auto threshold =
          ComputeRefineThreshold(param.update_fraction, indicators.local_error_indicators);
      const auto marked_elements =
          MarkedElements(threshold, indicators.local_error_indicators);

      if (param.construct_geometric_multigrid)
      {
        // If building a geometric multigrid, duplicate the final mesh.
        mesh.emplace_back(std::make_unique<mfem::ParMesh>(*mesh.back(), true));
      }

      // refine
      mesh.back()->GeneralRefinement(marked_elements, 1, param.max_nc_levels);
    }
    else if (use_coarsening)
    {
      MFEM_VERIFY(false, "Coarsening not implemented yet.");
      // coarsen mark
      // auto marked_elements = CoarsenMarker(
      //   iodata.model.refinement.adaptive.update_fraction,
      //   indicators.error_indicators);

      // mfem::Array<int> marked_elements;  // PLACEHOLDER

      // TODO: Compute error threshold that will trigger sufficient derefinement
      // const double threshold = ComputeThreshold(indicators.error_indicators);

      // const double threshold = 0;  // PLACEHOLDER

      // TODO: Figure out the method to expose here.
      // mesh.back()->NonconformingDerefinement(indicators.error_indicators,
      // param.max_nc_levels, 1);
    }

    // TODO: Measure this/make optional.
    // TODO: This is buggy, sometimes getting invalid elements.
    RebalanceMesh(mesh.back());

    // Solve + estimate.
    indicators = Solve(mesh, timer);

    ++iter;

    // Optionally save solution off
    if (param.save_step > 0 && iter % param.save_step == 0)
    {
      save_postprocess(iter);
    }
  }

  Mpi::Print("\nError Indicator: {:.3e}, DOF: {}\n", indicators.global_error_indicator,
             indicators.ndof);

  return indicators;
}

void BaseSolver::SaveMetadata(const mfem::ParFiniteElementSpace &fespace) const
{
  if (post_dir.length() == 0)
  {
    return;
  }
  const HYPRE_BigInt ndof = fespace.GlobalTrueVSize();
  HYPRE_BigInt ne = fespace.GetParMesh()->GetNE();
  Mpi::GlobalSum(1, &ne, fespace.GetComm());
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["Problem"]["MeshElements"] = ne;
    meta["Problem"]["DegreesOfFreedom"] = ndof;
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(int ksp_mult, int ksp_it) const
{
  if (post_dir.length() == 0)
  {
    return;
  }
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["LinearSolver"]["TotalSolves"] = ksp_mult;
    meta["LinearSolver"]["TotalIts"] = ksp_it;
    WriteMetadata(post_dir, meta);
  }
}

void BaseSolver::SaveMetadata(const Timer &timer) const
{
  if (post_dir.length() == 0)
  {
    return;
  }
  if (root)
  {
    json meta = LoadMetadata(post_dir);
    meta["ElapsedTime"]["Initialization"] = timer.GetAvgTime(Timer::INIT);
    meta["ElapsedTime"]["OperatorConstruction"] = timer.GetAvgTime(Timer::CONSTRUCT);
    meta["ElapsedTime"]["Solve"] = timer.GetAvgTime(Timer::SOLVE);
    meta["ElapsedTime"]["Postprocessing"] = timer.GetAvgTime(Timer::POSTPRO);
    meta["ElapsedTime"]["DiskIO"] = timer.GetAvgTime(Timer::IO);
    meta["ElapsedTime"]["Total"] = timer.GetAvgTime(Timer::TOTAL);
    WriteMetadata(post_dir, meta);
  }
}

namespace
{

struct EpsData
{
  const int idx;    // Domain or interface index
  const double pl;  // Participation ratio
  const double Ql;  // Quality factor
};

struct CapData
{
  const int idx;     // Surface index
  const double Cij;  // Capacitance (integrated charge)
};

struct IndData
{
  const int idx;     // Surface index
  const double Mij;  // Inductance (integrated flux)
};

struct ProbeData
{
  const int idx;                          // Probe index
  const std::complex<double> Fx, Fy, Fz;  // Field values at probe location
};

}  // namespace

void BaseSolver::PostprocessDomains(const PostOperator &postop, const std::string &name,
                                    int step, double time, double E_elec, double E_mag,
                                    double E_cap, double E_ind) const
{
  // If domains have been specified for postprocessing, compute the corresponding values
  // and write out to disk.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the field and lumped element energies.
  if (root)
  {
    std::string path = post_dir + "domain-E.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      // clang-format off
      output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s},{:>{}s}\n",
                   name, table.w1,
                   "E_elec", table.w,
                   "E_mag", table.w,
                   "E_cap", table.w,
                   "E_ind", table.w);
      // clang-format on
    }
    // clang-format off
    output.print("{:{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}\n",
                 time, table.w1, table.p1,
                 E_elec, table.w, table.p,
                 E_mag, table.w, table.p,
                 E_cap, table.w, table.p,
                 E_ind, table.w, table.p);
    // clang-format on
  }

  // Write the Q-factors due to bulk dielectric loss.
  std::vector<EpsData> eps_data;
  eps_data.reserve(postop.GetDomainPostOp().GetEps().size());
  for (const auto &[idx, data] : postop.GetDomainPostOp().GetEps())
  {
    const double pl = postop.GetBulkParticipation(idx, E_elec + E_cap);
    const double Ql = postop.GetBulkQualityFactor(idx, E_elec + E_cap);
    eps_data.push_back({idx, pl, Ql});
  }
  if (root && !eps_data.empty())
  {
    std::string path = post_dir + "domain-Q.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", name, table.w1);
      for (const auto &data : eps_data)
      {
        // clang-format off
        output.print("{:>{}s},{:>{}s}{}",
                     "p_bulk[" + std::to_string(data.idx) + "]", table.w,
                     "Q_bulk[" + std::to_string(data.idx) + "]", table.w,
                     (data.idx == eps_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", time, table.w1, table.p1);
    for (const auto &data : eps_data)
    {
      // clang-format off
      output.print("{:+{}.{}e},{:+{}.{}e}{}",
                   data.pl, table.w, table.p,
                   data.Ql, table.w, table.p,
                   (data.idx == eps_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void BaseSolver::PostprocessSurfaces(const PostOperator &postop, const std::string &name,
                                     int step, double time, double E_elec, double E_mag,
                                     double Vinc, double Iinc) const
{
  // If surfaces have been specified for postprocessing, compute the corresponding values
  // and write out to disk. This output uses the complex magnitude of the computed charge
  // and flux for frequency domain simulations. For capacitance/inductance, use the
  // excitation voltage or current across all sources and excited ports. The passed in
  // E_elec is the sum of the E-field and lumped capacitor energies, and E_mag is the same
  // for the B-field and lumped inductors.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the Q-factors due to interface dielectric loss.
  std::vector<EpsData> eps_data;
  eps_data.reserve(postop.GetSurfacePostOp().GetEps().size());
  for (const auto &[idx, data] : postop.GetSurfacePostOp().GetEps())
  {
    const double pl = postop.GetInterfaceParticipation(idx, E_elec);
    const double tandelta = postop.GetSurfacePostOp().GetInterfaceLossTangent(idx);
    const double Ql =
        (pl == 0.0 || tandelta == 0.0) ? mfem::infinity() : 1.0 / (tandelta * pl);
    eps_data.push_back({idx, pl, Ql});
  }
  if (root && !eps_data.empty())
  {
    std::string path = post_dir + "surface-Q.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", name, table.w1);
      for (const auto &data : eps_data)
      {
        // clang-format off
        output.print("{:>{}s},{:>{}s}{}",
                     "p_surf[" + std::to_string(data.idx) + "]", table.w,
                     "Q_surf[" + std::to_string(data.idx) + "]", table.w,
                     (data.idx == eps_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", time, table.w1, table.p1);
    for (const auto &data : eps_data)
    {
      // clang-format off
      output.print("{:+{}.{}e},{:+{}.{}e}{}",
                   data.pl, table.w, table.p,
                   data.Ql, table.w, table.p,
                   (data.idx == eps_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }

  // Write the surface capacitance (integrated charge).
  std::vector<CapData> cap_data;
  cap_data.reserve(postop.GetSurfacePostOp().GetCap().size());
  for (const auto &[idx, data] : postop.GetSurfacePostOp().GetCap())
  {
    const double Cij = (std::abs(Vinc) > 0.0) ? postop.GetSurfaceCharge(idx) / Vinc : 0.0;
    cap_data.push_back(
        {idx, iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, Cij)});
  }
  if (root && !cap_data.empty())
  {
    std::string path = post_dir + "surface-C.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", name, table.w1);
      for (const auto &data : cap_data)
      {
        // clang-format off
        output.print("{:>{}s}{}",
                     "C[" + std::to_string(data.idx) + "] (F)", table.w,
                     (data.idx == cap_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", time, table.w1, table.p1);
    for (const auto &data : cap_data)
    {
      // clang-format off
      output.print("{:+{}.{}e}{}",
                   data.Cij, table.w, table.p,
                   (data.idx == cap_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }

  // Write the surface inductance (integrated flux).
  std::vector<IndData> ind_data;
  ind_data.reserve(postop.GetSurfacePostOp().GetInd().size());
  for (const auto &[idx, data] : postop.GetSurfacePostOp().GetInd())
  {
    const double Mij = (std::abs(Iinc) > 0.0) ? postop.GetSurfaceFlux(idx) / Iinc : 0.0;
    ind_data.push_back(
        {idx, iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, Mij)});
  }
  if (root && !ind_data.empty())
  {
    std::string path = post_dir + "surface-M.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", name, table.w1);
      for (const auto &data : ind_data)
      {
        // clang-format off
        output.print("{:>{}s}{}",
                     "M[" + std::to_string(data.idx) + "] (H)", table.w,
                     (data.idx == ind_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", time, table.w1, table.p1);
    for (const auto &data : ind_data)
    {
      // clang-format off
      output.print("{:+{}.{}e}{}",
                   data.Mij, table.w, table.p,
                   (data.idx == ind_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void BaseSolver::PostprocessProbes(const PostOperator &postop, const std::string &name,
                                   int step, double time) const
{
#if defined(MFEM_USE_GSLIB)
  // If probe locations have been specified for postprocessing, compute the corresponding
  // field values and write out to disk.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the computed field values at probe locations.
  if (postop.GetProbes().size() == 0)
  {
    return;
  }
  const bool has_imaginary = postop.HasImaginary();
  for (int f = 0; f < 2; f++)
  {
    // Probe data is ordered as [Fx1, Fy1, Fz1, Fx2, Fy2, Fz2, ...].
    if (f == 0 && !postop.HasE())
    {
      continue;
    }
    if (f == 1 && !postop.HasB())
    {
      continue;
    }
    std::vector<ProbeData> probe_data;
    probe_data.reserve(postop.GetProbes().size());
    const std::vector<std::complex<double>> vF =
        (f == 0) ? postop.ProbeEField() : postop.ProbeBField();
    const int dim = vF.size() / postop.GetProbes().size();
    int i = 0;
    for (const auto &idx : postop.GetProbes())
    {
      probe_data.push_back(
          {idx, vF[i * dim], vF[i * dim + 1], (dim == 3) ? vF[i * dim + 2] : 0.0});
      i++;
    }
    const std::string F = (f == 0) ? "E" : "B";
    if (root && !probe_data.empty())
    {
      std::string path = post_dir + "probe-" + F + ".csv";
      auto output = OutputFile(path, (step > 0));
      if (step == 0)
      {
        output.print("{:>{}s},", name, table.w1);
        if (has_imaginary)
        {
          for (const auto &data : probe_data)
          {
            // clang-format off
            output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s}",
                         "Re{" + F + "_x[" + std::to_string(data.idx) + "]}", table.w,
                         "Im{" + F + "_x[" + std::to_string(data.idx) + "]}", table.w,
                         "Re{" + F + "_y[" + std::to_string(data.idx) + "]}", table.w,
                         "Im{" + F + "_y[" + std::to_string(data.idx) + "]}", table.w);
            // clang-format on
            if (dim == 3)
            {
              // clang-format off
              output.print(",{:>{}s},{:>{}s}{}",
                           "Re{" + F + "_z[" + std::to_string(data.idx) + "]}", table.w,
                           "Im{" + F + "_z[" + std::to_string(data.idx) + "]}", table.w,
                           (data.idx == probe_data.back().idx) ? "" : ",");
              // clang-format on
            }
            else
            {
              // clang-format off
              output.print("{}",
                           (data.idx == probe_data.back().idx) ? "" : ",");
              // clang-format on
            }
          }
        }
        else
        {
          for (const auto &data : probe_data)
          {
            // clang-format off
            output.print("{:>{}s},{:>{}s}",
                         F + "_x[" + std::to_string(data.idx) + "]", table.w,
                         F + "_y[" + std::to_string(data.idx) + "]", table.w);
            // clang-format on
            if (dim == 3)
            {
              // clang-format off
              output.print(",{:>{}s}{}",
                           F + "_z[" + std::to_string(data.idx) + "]", table.w,
                           (data.idx == probe_data.back().idx) ? "" : ",");
              // clang-format on
            }
            else
            {
              // clang-format off
              output.print("{}",
                           (data.idx == probe_data.back().idx) ? "" : ",");
              // clang-format on
            }
          }
        }
        output.print("\n");
      }
      output.print("{:{}.{}e},", time, table.w1, table.p1);
      if (has_imaginary)
      {
        for (const auto &data : probe_data)
        {
          // clang-format off
          output.print("{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}",
                       data.Fx.real(), table.w, table.p,
                       data.Fx.imag(), table.w, table.p,
                       data.Fy.real(), table.w, table.p,
                       data.Fy.imag(), table.w, table.p);
          // clang-format on
          if (dim == 3)
          {
            // clang-format off
            output.print(",{:+{}.{}e},{:+{}.{}e}{}",
                         data.Fz.real(), table.w, table.p,
                         data.Fz.imag(), table.w, table.p,
                         (data.idx == probe_data.back().idx) ? "" : ",");
            // clang-format on
          }
          else
          {
            // clang-format off
            output.print("{}",
                         (data.idx == probe_data.back().idx) ? "" : ",");
            // clang-format on
          }
        }
      }
      else
      {
        for (const auto &data : probe_data)
        {
          // clang-format off
          output.print("{:+{}.{}e},{:+{}.{}e}",
                       data.Fx.real(), table.w, table.p,
                       data.Fy.real(), table.w, table.p);
          // clang-format on
          if (dim == 3)
          {
            // clang-format off
            output.print(",{:+{}.{}e}{}",
                         data.Fz.real(), table.w, table.p,
                         (data.idx == probe_data.back().idx) ? "" : ",");
            // clang-format on
          }
          else
          {
            // clang-format off
            output.print("{}",
                         (data.idx == probe_data.back().idx) ? "" : ",");
            // clang-format on
          }
        }
      }
      output.print("\n");
    }
  }
#endif
}

void BaseSolver::PostprocessFields(const PostOperator &postop, int step, double time) const
{
  // Save the computed fields in parallel in format for viewing with ParaView.
  if (post_dir.length() == 0)
  {
    Mpi::Warning("No file specified under [\"Problem\"][\"Output\"]!\nSkipping saving of "
                 "fields to disk!\n");
    return;
  }
  postop.WriteFields(step, time);
  Mpi::Barrier();
}

}  // namespace palace
