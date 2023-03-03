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
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
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
  : iodata(iodata_), post_dir_(GetPostDir(iodata_.problem.output)), root(root_),
    table(8, 9, 6)
{
  // Create directory for output.
  if (root && !std::filesystem::exists(post_dir_))
  {
    std::filesystem::create_directories(post_dir_);
  }

  // Initialize simulation metadata for this simulation.
  if (root && post_dir_.length() > 0)
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
    WriteMetadata(post_dir_, meta);
  }
}

namespace
{
// Given a vector of estimates local to this rank, compute a threshold that
// results in a Dorfler marking across processor ranks.
double ComputeRefineThreshold(double fraction, std::vector<double> estimates)
{

  std::sort(estimates.begin(), estimates.end());

  const auto pivot = static_cast<std::size_t>(
      fraction * std::distance(estimates.begin(), estimates.end()));

  double threshold = estimates[pivot];

  auto count_elems_in_part = [&estimates](double x)
  {
    const auto upper = std::upper_bound(estimates.begin(), estimates.end(), x);
    return std::distance(estimates.begin(), upper);
  };

  // Each rank has computed a different threshold, if a given rank has lots of
  // low error elements, their value will be lower and if a rank has high error,
  // their value will be higher. Using the value from the low error rank will
  // give too few elements, and using the value from the high error rank will
  // give too many. The correct threshold value will be an intermediate
  // between the min and max over ranks. To compute the correct value, can
  // perform bisection search on rank 0.

  auto comm = Mpi::World();

  double min_threshold;
  double max_threshold;
  long long num_elem = estimates.size();
  const long long max_elem = [&]() {  // IILE for const
    long long max_elem = 0;
    MPI_Reduce(&num_elem, &max_elem, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
    return max_elem;
  }();

  long long candidate_total = 0, old_candidate_total = max_elem;

  // Collect to rank zero the initial threshold bounds, and the total number of elements.
  MPI_Reduce(&threshold, &min_threshold, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&threshold, &max_threshold, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

  bool continue_bisecting = true;
  MPI_Request request = MPI_REQUEST_NULL;

  while (continue_bisecting)
  {
    if (Mpi::Root(comm))
    {
      // Root is going to bisect, and send thresholds to other ranks
      constexpr double tol = 1e-3;
      while ((max_threshold - min_threshold) > tol)
      {
        // Compute a new threshold and send to all ranks.
        threshold = (min_threshold + max_threshold) / 2;

        // Send to all processors
        for (int i = 1; i < Mpi::Size(comm); i++)
        {
          MPI_Isend(&threshold, 1, MPI_DOUBLE, i, 0, comm, &request);
        }
      }
    }
    else
    {
      // Receive the new candidate threshold value from rank 0
      MPI_Irecv(&threshold, 1, MPI_DOUBLE, 0, 0, comm, &request);
    }

    // Given the now agreed upon threshold value, count the number of elements
    // this would mark across all ranks.
    num_elem = count_elems_in_part(threshold);
    MPI_Reduce(&num_elem, &candidate_total, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);

    if (Mpi::Root(comm))
    {

      // Check 1: That the number of candidate elements changed. This is important
      // given working with discrete sets. The threshold value might keep changing
      // but no outward effect.

      if (candidate_total == old_candidate_total)
      {
        // changes in the threshold value are no longer making a difference.
        continue_bisecting = false;
      }
      else
      {
        // Root will check the effective fraction given this threshold
        double candidate_fraction =
            static_cast<double>(static_cast<long double>(candidate_total) / max_elem);

        constexpr double tol = 1e-3;
        if (std::abs(candidate_fraction - fraction) < tol)
        {
          // Close enough to the threshold value
          continue_bisecting = false;
        }
        else if (candidate_fraction < fraction)
        {
          min_threshold = threshold;
        }
        else if (candidate_fraction > fraction)
        {
          max_threshold = threshold;
        }
      }
    }

    // Inform all ranks if the loop will be continuing. All ranks already have
    // the threshold.
    MPI_Bcast(&continue_bisecting, 1, MPI_LOGICAL, 0, comm);
  }

  return threshold;
}
}  // namespace

BaseSolver::ErrorIndicators
BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                                    Timer &timer) const
{
  const auto &param = iodata.model.refinement.adaptation;

  const bool use_coarsening = param.coarsening_fraction > 0;
  const int use_nc = use_coarsening ? 1 : -1;

  int iter = 0;
  auto estimates = Solve(mesh, timer, iter++);

  while (iter < param.min_its ||
         (iter < param.max_its && estimates.global_error_indicator > param.tolerance))
  {
    if (estimates.ndof < param.dof_limit)
    {
      // refinement mark
      // auto marked_elements = RefineMarker(
      //   iodata.model.refinement.adaptive.update_fraction,
      //   estimates.error_estimates);

      mfem::Array<int> marked_elements;  // PLACEHOLDER

      // refine
      if (param.construct_geometric_multigrid)
      {
        mesh.emplace_back(std::make_unique<mfem::ParMesh>(*mesh.back(), true));
      }

      // If coarsening, need NC, otherwise let the mesh decide
      mesh.back()->GeneralRefinement(marked_elements, use_nc, param.max_nc_levels);
    }
    else if (use_coarsening)
    {
      // coarsen mark
      // auto marked_elements = CoarsenMarker(
      //   iodata.model.refinement.adaptive.update_fraction,
      //   estimates.error_estimates);

      mfem::Array<int> marked_elements;  // PLACEHOLDER

      // TODO: Compute error threshold that will trigger sufficient derefinement
      // const double threshold = ComputeThreshold(estimates.error_estimates);

      const double threshold = 0;  // PLACEHOLDER

      // TODO: Figure out the method to expose here.
      // mesh.back()->NonconformingDerefinement(estimates.error_estimates,
      // param.max_nc_levels, 1);
    }

    // solve + estimate
    estimates = Solve(mesh, timer, iter);
    iter++;
  }

  return estimates;
}

void BaseSolver::SaveMetadata(const std::string &post_dir,
                              const mfem::ParFiniteElementSpace &fespace) const
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

void BaseSolver::SaveMetadata(const std::string &post_dir, int ksp_mult, int ksp_it) const
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

void BaseSolver::SaveMetadata(const std::string &post_dir, const Timer &timer) const
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

void BaseSolver::PostprocessDomains(const std::string &post_dir, const PostOperator &postop,
                                    const std::string &name, int step, double time,
                                    double E_elec, double E_mag, double E_cap,
                                    double E_ind) const
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

void BaseSolver::PostprocessSurfaces(const std::string &post_dir,
                                     const PostOperator &postop, const std::string &name,
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

void BaseSolver::PostprocessProbes(const std::string &post_dir, const PostOperator &postop,
                                   const std::string &name, int step, double time) const
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

void BaseSolver::PostprocessFields(const std::string &post_dir, const PostOperator &postop,
                                   int step, double time) const
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

void BaseSolver::ErrorReductionOperator::operator()(BaseSolver::ErrorIndicators &ebar,
                                                    std::vector<double> &&ind) const
{

  // Update the maximum global error.
  ebar.global_error_indicator =
      std::max(ebar.global_error_indicator, std::accumulate(ind.begin(), ind.end(), 0.0));

  // update the average local indicator. Using running average update rather
  // than sum and final division to maintain validity at all times.
  auto running_average = [](const auto &xbar, const auto &x)
  { return (xbar * n + x) / (n + 1); };

  if (n > 0)
  {
    MFEM_VERIFY(ebar.local_error_indicators.size() == ind.size(),
                "Local error indicator vectors mismatch.");
    // Combine these error indicators into the current average.
    std::transform(ebar.local_error_indicators.begin(), ebar.local_error_indicators.end(),
                   ind.begin(), ebar.local_error_indicators.begin(), running_average);
  }
  else
  {
    // This is the first sample, just steal the data.
    ebar.local_error_indicators = std::move(ind);
  }

  // Another sample has been added, increment for the running average lambda.
  ++n;
}

}  // namespace palace
