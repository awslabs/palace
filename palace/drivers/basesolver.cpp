// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <complex>
#include <numeric>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/transientsolver.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "models/domainpostoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/dorfler.hpp"
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

BaseSolver::BaseSolver(const IoData &iodata, bool root, int size, int num_thread,
                       const char *git_tag)
  : iodata(iodata), post_dir(GetPostDir(iodata.problem.output)), root(root), table(8, 9, 6)
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

void RebalanceMesh(std::unique_ptr<mfem::ParMesh> &mesh, double maximum_imbalance)
{
  auto comm = Mpi::World();
  if (Mpi::Size(comm) > 1)
  {
    int min_elem, max_elem;
    min_elem = max_elem = mesh->GetNE();

    Mpi::GlobalMin(1, &min_elem, comm);
    Mpi::GlobalMax(1, &max_elem, comm);

    const double ratio = double(max_elem)/min_elem;
    Mpi::Print("Min Elem per processor: {}, Max Elem per processor: {}, Ratio: {:.3e}\n",
                min_elem, max_elem, double(max_elem)/min_elem);

    if (ratio > maximum_imbalance)
    {
      if (mesh->Nonconforming())
      {
        Mpi::Print("Ratio {:.3e} exceeds maximum allowed value {} -> Rebalancing.", ratio, maximum_imbalance);
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

  mesh->FinalizeTopology();
  mesh->Finalize(true);

  // If the mesh is higher order, synchronize through the nodal grid function.
  // This will in turn call the mesh exchange of face neighbor data.
  if (mesh->GetNodes())
  {
    auto *pgf = dynamic_cast<mfem::ParGridFunction *>(mesh->GetNodes());
    MFEM_ASSERT(pgf, "The grid function must be castable to a ParGridFunction");
    pgf->ExchangeFaceNbrData();
  }
  else
  {
    mesh->ExchangeFaceNbrData();
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
  const bool use_coarsening = mesh.back()->Nonconforming() && param.coarsening_fraction > 0;

  if (param.coarsening_fraction > 0 && mesh.back()->Conforming())
  {
    Mpi::Warning("{}\n{}\n", "Non-zero coarsening fraction is being ignored.",
                 "Coarsening can only occur if a mesh is in nonconforming mode.");
  }

  if (use_amr && mesh.size() > 1)
  {
    Mpi::Warning("{}\n", "Flattening mesh sequence: AMR will solve only use the final "
                         "mesh from the refinement sequence.");
    mesh.erase(mesh.begin(), mesh.end() - 1);
  }

  mesh.back()->ExchangeFaceNbrData();

  int iter = 0;
  auto indicators = Solve(mesh, timer);

  if (use_amr)
  {
    const auto is_transient = dynamic_cast<const TransientSolver *>(this) != nullptr;
    if (is_transient)
    {
      Mpi::Warning("{}\n", "AMR is not currently supported for transient simulations");
    }
    else
    {
      Mpi::Print("\nAdaptive Mesh Refinement Parameters:\n");
      Mpi::Print("MinIter: {}, MaxIter: {}, Tolerance: {:.3e}, DOFLimit: {}\n\n",
                 param.min_its, param.max_its, param.tolerance, param.dof_limit);
      save_postprocess(iter);  // Save the initial solution
    }
  }

  // Collection of all tests that might exhaust resources.
  auto exhausted_resources = [&]()
  {
    bool ret = false;
    // run out of DOFs, and coarsening isn't allowed.
    ret |= (!use_coarsening && indicators.ndof > param.dof_limit && param.dof_limit > 0);
    // run out of iterations
    ret |= iter >= param.max_its;
    // There are no resources for transient amr.
    // TODO: remove this once transient simulations are supported.
    ret |= dynamic_cast<const TransientSolver *>(this) != nullptr;

    return ret;
  };

  while ((iter < param.min_its || indicators.global_error_indicator > param.tolerance) &&
         !exhausted_resources())
  {
    Mpi::Print("Adaptation iteration {}: Error Indicator: {:.3e}, DOF: {}\n", iter,
               indicators.global_error_indicator, indicators.ndof);
    if (indicators.ndof < param.dof_limit)
    {
      const auto threshold = utils::ComputeDorflerThreshold(
          param.update_fraction, indicators.local_error_indicators);
      const auto marked_elements =
          MarkedElements(threshold, indicators.local_error_indicators);

      const auto initial_elem_count = mesh.back()->GetGlobalNE();
      mesh.back()->GeneralRefinement(marked_elements, 1, param.max_nc_levels);
      const auto final_elem_count = mesh.back()->GetGlobalNE();
      Mpi::Print("Mesh refinement added {} elements. Initial: {}, Final: {}\n",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }
    else if (use_coarsening)
    {
      // Perform a Dorfler style marking looking for the largest number of
      // derefinement opportunities to represent a fraction of the derefinable error.
      const auto &derefinement_table = mesh.back()->pncmesh->GetDerefinementTable();

      mfem::Vector coarse_error(derefinement_table.Size());
      for (int i = 0; i < derefinement_table.Size(); ++i)
      {
        mfem::Array<int> row;
        derefinement_table.GetRow(i, row);

        // sum the error for all sub elements that can be combined
        coarse_error[i] =
            std::accumulate(row.begin(), row.end(), 0.0,
                            [&indicators](double s, int i)
                            { return s += indicators.local_error_indicators[i]; });
      }

      // Given the coarse errors, we use the Dorfler marking strategy to
      // identify the smallest set of elements that make up (1 - θ) of the total
      // error, where θ is the coarsening fraction. The complement of this set
      // is then the largest number of elements that make up θ of the total error.

      const double threshold =
          utils::ComputeDorflerThreshold(1 - param.coarsening_fraction, coarse_error);

      const auto initial_elem_count = mesh.back()->GetGlobalNE();
      mesh.back()->DerefineByError(indicators.local_error_indicators, threshold,
                                   param.max_nc_levels);
      const auto final_elem_count = mesh.back()->GetGlobalNE();
      Mpi::Print("Mesh coarsening removed {} elements. Initial: {}, Final: {}\n",
                 initial_elem_count - final_elem_count, initial_elem_count,
                 final_elem_count);
    }

    RebalanceMesh(mesh.back(), iodata.model.refinement.adaptation.maximum_imbalance);

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
void BaseSolver::SaveMetadata(const mfem::ParFiniteElementSpaceHierarchy &fespaces) const
{
  if (post_dir.length() == 0)
  {
    return;
  }
  const mfem::ParFiniteElementSpace &fespace = fespaces.GetFinestFESpace();
  HYPRE_BigInt ne = fespace.GetParMesh()->GetNE();
  Mpi::GlobalSum(1, &ne, fespace.GetComm());
  std::vector<HYPRE_BigInt> ndofs(fespaces.GetNumLevels());
  for (int l = 0; l < fespaces.GetNumLevels(); l++)
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
  if (post_dir.length() == 0)
  {
    return;
  }
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

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
