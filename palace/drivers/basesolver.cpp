// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "basesolver.hpp"

#include <array>
#include <complex>
#include <numeric>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "drivers/transientsolver.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/ksp.hpp"
#include "models/domainpostoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/communication.hpp"
#include "utils/dorfler.hpp"
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

std::string GetIterationPostDir(const std::string &output, int step, int width)
{
  return fmt::format("{}iteration{:0{}d}/", output, step, width);
}

void SaveIteration(MPI_Comm comm, const std::string &output, int step, int width)
{
  namespace fs = std::filesystem;
  BlockTimer bt(Timer::IO);
  Mpi::Barrier(comm);  // Wait for all processes to write postprocessing files
  if (Mpi::Root(comm))
  {
    // Create a subfolder for the results of this adaptation.
    const std::string step_output = GetIterationPostDir(output, step, width);
    if (!fs::exists(step_output))
    {
      fs::create_directories(step_output);
    }
    constexpr auto options =
        fs::copy_options::recursive | fs::copy_options::overwrite_existing;
    for (const auto &f : fs::directory_iterator(output))
    {
      if (f.path().filename().string().rfind("iteration") == 0)
      {
        continue;
      }
      fs::copy(f, step_output + f.path().filename().string(), options);
    }
  }
  Mpi::Barrier(comm);
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
  : iodata(iodata), post_dir(GetPostDir(iodata.problem.output)), root(root), table(8, 9, 9)
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

void BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &refinement = iodata.model.refinement;
  const bool use_amr = [&]()
  {
    if (refinement.max_it > 0 && dynamic_cast<const TransientSolver *>(this) != nullptr)
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
  while (!ExhaustedResources(it, ntdof) && err >= refinement.tol)
  {
    // Print timing summary.
    Mpi::Print("\nCumulative timing statistics:\n");
    BlockTimer::Print(comm);
    SaveMetadata(BlockTimer::GlobalTimer());

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
      fine_mesh.GeneralRefinement(marked_elements, -1, refinement.max_nc_levels);
      const auto final_elem_count = fine_mesh.GetGlobalNE();
      Mpi::Print(" {} mesh refinement added {:d} elements (initial = {:d}, final = {:d})\n",
                 fine_mesh.Nonconforming() ? "Nonconforming" : "Conforming",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }

    // Optionally rebalance and write the adapted mesh to file.
    {
      const auto ratio_pre = mesh::RebalanceMesh(iodata, *mesh.back());
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
    }

    // Solve + estimate.
    Mpi::Print("\nProceeding with solve/estimate iteration {}...\n", it + 1);
    std::tie(indicators, ntdof) = Solve(mesh);
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
  if (post_dir.length() == 0)
  {
    return;
  }
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
    for (int i = Timer::INIT; i < Timer::NUM_TIMINGS; i++)
    {
      auto key = Timer::descriptions[i];
      key.erase(std::remove_if(key.begin(), key.end(), isspace), key.end());
      meta["ElapsedTime"]["Durations"][key] = timer.Data((Timer::Index)i);
      meta["ElapsedTime"]["Counts"][key] = timer.Counts((Timer::Index)i);
    }
    WriteMetadata(post_dir, meta);
  }
}

namespace
{

struct EnergyData
{
  const int idx;        // Domain index
  const double E_elec;  // Electric field energy
  const double E_mag;   // Magnetic field energy
};

struct FluxData
{
  const int idx;                   // Surface index
  const std::complex<double> Phi;  // Integrated flux
  const SurfaceFluxType type;      // Flux type
};

struct EpsData
{
  const int idx;   // Interface index
  const double p;  // Participation ratio
  const double Q;  // Quality factor
};

struct ProbeData
{
  const int idx;                          // Probe index
  const std::complex<double> Fx, Fy, Fz;  // Field values at probe location
};

}  // namespace

void BaseSolver::PostprocessDomains(const PostOperator &post_op, const std::string &name,
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
  std::vector<EnergyData> energy_data;
  energy_data.reserve(post_op.GetDomainPostOp().M_i.size());
  for (const auto &[idx, data] : post_op.GetDomainPostOp().M_i)
  {
    const double E_elec_i = (E_elec > 0.0) ? post_op.GetEFieldEnergy(idx) : 0.0;
    const double E_mag_i = (E_mag > 0.0) ? post_op.GetHFieldEnergy(idx) : 0.0;
    energy_data.push_back({idx, E_elec_i, E_mag_i});
  }
  if (root)
  {
    std::string path = post_dir + "domain-E.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      // clang-format off
      output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s},{:>{}s}{}",
                   name, table.w1,
                   "E_elec (J)", table.w,
                   "E_mag (J)", table.w,
                   "E_cap (J)", table.w,
                   "E_ind (J)", table.w,
                   energy_data.empty() ? "" : ",");
      // clang-format on
      for (const auto &data : energy_data)
      {
        // clang-format off
        output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s}{}",
                     "E_elec[" + std::to_string(data.idx) + "] (J)", table.w,
                     "p_elec[" + std::to_string(data.idx) + "]", table.w,
                     "E_mag[" + std::to_string(data.idx) + "] (J)", table.w,
                     "p_mag[" + std::to_string(data.idx) + "]", table.w,
                     (data.idx == energy_data.back().idx) ? "" : ",");
        // clang-format on
      }
      output.print("\n");
    }
    // clang-format off
    output.print("{:{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}{}",
                 time, table.w1, table.p1,
                 iodata.DimensionalizeValue(IoData::ValueType::ENERGY, E_elec),
                 table.w, table.p,
                 iodata.DimensionalizeValue(IoData::ValueType::ENERGY, E_mag),
                 table.w, table.p,
                 iodata.DimensionalizeValue(IoData::ValueType::ENERGY, E_cap),
                 table.w, table.p,
                 iodata.DimensionalizeValue(IoData::ValueType::ENERGY, E_ind),
                 table.w, table.p,
                 energy_data.empty() ? "" : ",");
    // clang-format on
    for (const auto &data : energy_data)
    {
      // clang-format off
      output.print("{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}{}",
                   iodata.DimensionalizeValue(IoData::ValueType::ENERGY, data.E_elec),
                   table.w, table.p,
                   (std::abs(E_elec) > 0.0) ? data.E_elec / E_elec : 0.0, table.w, table.p,
                   iodata.DimensionalizeValue(IoData::ValueType::ENERGY, data.E_mag),
                   table.w, table.p,
                   (std::abs(E_mag) > 0.0) ? data.E_mag / E_mag : 0.0, table.w, table.p,
                   (data.idx == energy_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void BaseSolver::PostprocessSurfaces(const PostOperator &post_op, const std::string &name,
                                     int step, double time, double E_elec,
                                     double E_mag) const
{
  // If surfaces have been specified for postprocessing, compute the corresponding values
  // and write out to disk. The passed in E_elec is the sum of the E-field and lumped
  // capacitor energies, and E_mag is the same for the B-field and lumped inductors.
  if (post_dir.length() == 0)
  {
    return;
  }

  // Write the integrated surface flux.
  const bool has_imaginary = post_op.HasImag();
  std::vector<FluxData> flux_data;
  flux_data.reserve(post_op.GetSurfacePostOp().flux_surfs.size());
  for (const auto &[idx, data] : post_op.GetSurfacePostOp().flux_surfs)
  {
    const std::complex<double> Phi = post_op.GetSurfaceFlux(idx);
    double scale = 1.0;
    switch (data.type)
    {
      case SurfaceFluxType::ELECTRIC:
        scale = iodata.DimensionalizeValue(IoData::ValueType::CAPACITANCE, 1.0);
        scale *= iodata.DimensionalizeValue(IoData::ValueType::VOLTAGE, 1.0);
        break;
      case SurfaceFluxType::MAGNETIC:
        scale = iodata.DimensionalizeValue(IoData::ValueType::INDUCTANCE, 1.0);
        scale *= iodata.DimensionalizeValue(IoData::ValueType::CURRENT, 1.0);
        break;
      case SurfaceFluxType::POWER:
        scale = iodata.DimensionalizeValue(IoData::ValueType::POWER, 1.0);
        break;
    }
    flux_data.push_back({idx, Phi * scale, data.type});
  }
  if (root && !flux_data.empty())
  {
    std::string path = post_dir + "surface-F.csv";
    auto output = OutputFile(path, (step > 0));
    if (step == 0)
    {
      output.print("{:>{}s},", name, table.w1);
      for (const auto &data : flux_data)
      {
        std::string name, unit;
        switch (data.type)
        {
          case SurfaceFluxType::ELECTRIC:
            name = "Φ_elec";
            unit = "(C)";
            break;
          case SurfaceFluxType::MAGNETIC:
            name = "Φ_mag";
            unit = "(Wb)";
            break;
          case SurfaceFluxType::POWER:
            name = "Φ_pow";
            unit = "(W)";
            break;
        }
        if (has_imaginary && data.type != SurfaceFluxType::POWER)
        {
          // clang-format off
          output.print("{:>{}s},{:>{}s}{}",
                       "Re{" + name + "[" + std::to_string(data.idx) + "]} " + unit, table.w,
                       "Im{" + name + "[" + std::to_string(data.idx) + "]} " + unit, table.w,
                       (data.idx == flux_data.back().idx) ? "" : ",");
          // clang-format on
        }
        else
        {
          // clang-format off
          output.print("{:>{}s}{}",
                       name + "[" + std::to_string(data.idx) + "] " + unit, table.w,
                       (data.idx == flux_data.back().idx) ? "" : ",");
          // clang-format on
        }
      }
      output.print("\n");
    }
    output.print("{:{}.{}e},", time, table.w1, table.p1);
    for (const auto &data : flux_data)
    {
      if (has_imaginary && data.type != SurfaceFluxType::POWER)
      {
        // clang-format off
        output.print("{:+{}.{}e},{:+{}.{}e}{}",
                     data.Phi.real(), table.w, table.p,
                     data.Phi.imag(), table.w, table.p,
                     (data.idx == flux_data.back().idx) ? "" : ",");
        // clang-format on
      }
      else
      {
        // clang-format off
        output.print("{:+{}.{}e}{}",
                     data.Phi.real(), table.w, table.p,
                     (data.idx == flux_data.back().idx) ? "" : ",");
        // clang-format on
      }
    }
    output.print("\n");
  }

  // Write the Q-factors due to interface dielectric loss.
  std::vector<EpsData> eps_data;
  eps_data.reserve(post_op.GetSurfacePostOp().eps_surfs.size());
  for (const auto &[idx, data] : post_op.GetSurfacePostOp().eps_surfs)
  {
    const double p = post_op.GetInterfaceParticipation(idx, E_elec);
    const double tandelta = post_op.GetSurfacePostOp().GetInterfaceLossTangent(idx);
    const double Q =
        (p == 0.0 || tandelta == 0.0) ? mfem::infinity() : 1.0 / (tandelta * p);
    eps_data.push_back({idx, p, Q});
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
                   data.p, table.w, table.p,
                   data.Q, table.w, table.p,
                   (data.idx == eps_data.back().idx) ? "" : ",");
      // clang-format on
    }
    output.print("\n");
  }
}

void BaseSolver::PostprocessProbes(const PostOperator &post_op, const std::string &name,
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
  if (post_op.GetProbes().size() == 0)
  {
    return;
  }
  const bool has_imaginary = post_op.HasImag();
  for (int f = 0; f < 2; f++)
  {
    // Probe data is ordered as [Fx1, Fy1, Fz1, Fx2, Fy2, Fz2, ...].
    if (f == 0 && !post_op.HasE())
    {
      continue;
    }
    if (f == 1 && !post_op.HasB())
    {
      continue;
    }
    const std::string F = (f == 0) ? "E" : "B";
    const std::string unit = (f == 0) ? "(V/m)" : "(Wb/m²)";
    const auto type = (f == 0) ? IoData::ValueType::FIELD_E : IoData::ValueType::FIELD_B;
    const auto vF = (f == 0) ? post_op.ProbeEField() : post_op.ProbeBField();
    const int dim = vF.size() / post_op.GetProbes().size();
    std::vector<ProbeData> probe_data;
    probe_data.reserve(post_op.GetProbes().size());
    int i = 0;
    for (const auto &idx : post_op.GetProbes())
    {
      probe_data.push_back(
          {idx, iodata.DimensionalizeValue(type, vF[i * dim]),
           iodata.DimensionalizeValue(type, vF[i * dim + 1]),
           (dim == 3) ? iodata.DimensionalizeValue(type, vF[i * dim + 2]) : 0.0});
      i++;
    }
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
                         "Re{" + F + "_x[" + std::to_string(data.idx) + "]} " + unit, table.w,
                         "Im{" + F + "_x[" + std::to_string(data.idx) + "]} " + unit, table.w,
                         "Re{" + F + "_y[" + std::to_string(data.idx) + "]} " + unit, table.w,
                         "Im{" + F + "_y[" + std::to_string(data.idx) + "]} " + unit, table.w);
            // clang-format on
            if (dim == 3)
            {
              // clang-format off
              output.print(",{:>{}s},{:>{}s}{}",
                           "Re{" + F + "_z[" + std::to_string(data.idx) + "]} " + unit, table.w,
                           "Im{" + F + "_z[" + std::to_string(data.idx) + "]} " + unit, table.w,
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
                         F + "_x[" + std::to_string(data.idx) + "] " + unit, table.w,
                         F + "_y[" + std::to_string(data.idx) + "] " + unit, table.w);
            // clang-format on
            if (dim == 3)
            {
              // clang-format off
              output.print(",{:>{}s}{}",
                           F + "_z[" + std::to_string(data.idx) + "] " + unit, table.w,
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

void BaseSolver::PostprocessFields(const PostOperator &post_op, int step, double time) const
{
  // Save the computed fields in parallel in format for viewing with ParaView.
  BlockTimer bt(Timer::IO);
  if (post_dir.length() == 0)
  {
    Mpi::Warning(post_op.GetComm(),
                 "No file specified under [\"Problem\"][\"Output\"]!\nSkipping saving of "
                 "fields to disk!\n");
    return;
  }
  post_op.WriteFields(step, time);
  Mpi::Barrier(post_op.GetComm());
}

void BaseSolver::PostprocessErrorIndicator(const PostOperator &post_op,
                                           const ErrorIndicator &indicator,
                                           bool fields) const
{
  // Write the indicator statistics.
  if (post_dir.length() == 0)
  {
    return;
  }
  MPI_Comm comm = post_op.GetComm();
  std::array<double, 4> data = {indicator.Norml2(comm), indicator.Min(comm),
                                indicator.Max(comm), indicator.Mean(comm)};
  if (root)
  {
    std::string path = post_dir + "error-indicators.csv";
    auto output = OutputFile(path, false);
    // clang-format off
    output.print("{:>{}s},{:>{}s},{:>{}s},{:>{}s}\n",
                 "Norm", table.w,
                 "Minimum", table.w,
                 "Maximum", table.w,
                 "Mean", table.w);
    output.print("{:+{}.{}e},{:+{}.{}e},{:+{}.{}e},{:+{}.{}e}\n",
                 data[0], table.w, table.p,
                 data[1], table.w, table.p,
                 data[2], table.w, table.p,
                 data[3], table.w, table.p);
    // clang-format on
  }
  if (fields)
  {
    BlockTimer bt(Timer::IO);
    post_op.WriteFieldsFinal(&indicator);
    Mpi::Barrier(post_op.GetComm());
  }
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
