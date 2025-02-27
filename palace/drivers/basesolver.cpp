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
#include "models/spaceoperator.hpp"
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

std::string GetIterationPostDir(const std::string &output, int step, int width)
{
  return fs::path(output) / fmt::format("iteration{:0{}d}/", step, width);
}

void SaveIteration(MPI_Comm comm, const std::string &output, int step, int width)
{
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
  : iodata(iodata), post_dir(iodata.problem.output), root(root)
{
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

BaseSolver::DomainsPostPrinter::DomainsPostPrinter(bool do_measurement, bool root,
                                                   const std::string &post_dir,
                                                   const DomainPostOperator &dom_post_op,
                                                   const std::string &idx_col_name,
                                                   int n_expected_rows)
  : do_measurement_{do_measurement}, root_{root}
{
  do_measurement_ = do_measurement_            //
                    && post_dir.length() > 0;  // Valid output dir

  if (!do_measurement_ || !root_)
  {
    return;
  }
  using fmt::format;

  domain_E = TableWithCSVFile(post_dir + "domain-E.csv");
  domain_E.table.reserve(n_expected_rows, 4 + dom_post_op.M_i.size());
  domain_E.table.insert_column(Column("idx", idx_col_name, 0, {}, {}, ""));

  domain_E.table.insert_column("Ee", "E_elec (J)");
  domain_E.table.insert_column("Em", "E_mag (J)");
  domain_E.table.insert_column("Ec", "E_cap (J)");
  domain_E.table.insert_column("Ei", "E_ind (J)");

  for (const auto &[idx, data] : dom_post_op.M_i)
  {
    domain_E.table.insert_column(format("Ee{}", idx), format("E_elec[{}] (J)", idx));
    domain_E.table.insert_column(format("pe{}", idx), format("p_elec[{}]", idx));
    domain_E.table.insert_column(format("Em{}", idx), format("E_mag[{}] (J)", idx));
    domain_E.table.insert_column(format("pm{}", idx), format("p_mag[{}]", idx));
  }
  domain_E.AppendHeader();
}

void BaseSolver::DomainsPostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                    const PostOperator &post_op,
                                                    double E_elec, double E_mag,
                                                    double E_cap, double E_ind,
                                                    const IoData &iodata)
{
  if (!do_measurement_)
  {
    return;
  }

  // MPI Gather
  std::vector<EnergyData> energy_data;
  energy_data.reserve(post_op.GetDomainPostOp().M_i.size());
  for (const auto &[idx, data] : post_op.GetDomainPostOp().M_i)
  {
    double E_elec_i = (E_elec > 0.0) ? post_op.GetEFieldEnergy(idx) : 0.0;
    double E_mag_i = (E_mag > 0.0) ? post_op.GetHFieldEnergy(idx) : 0.0;
    energy_data.emplace_back(EnergyData{idx, E_elec_i, E_mag_i});
  }

  if (!root_)
  {
    return;
  }

  using VT = IoData::ValueType;
  using fmt::format;

  domain_E.table["idx"] << idx_value_dimensionful;

  domain_E.table["Ee"] << iodata.DimensionalizeValue(VT::ENERGY, E_elec);
  domain_E.table["Em"] << iodata.DimensionalizeValue(VT::ENERGY, E_mag);
  domain_E.table["Ec"] << iodata.DimensionalizeValue(VT::ENERGY, E_cap);
  domain_E.table["Ei"] << iodata.DimensionalizeValue(VT::ENERGY, E_ind);

  // Write the field and lumped element energies.
  for (const auto &[idx, E_e, E_m] : energy_data)
  {
    domain_E.table[format("Ee{}", idx)] << iodata.DimensionalizeValue(VT::ENERGY, E_e);
    domain_E.table[format("pe{}", idx)]
        << ((std::abs(E_elec) > 0.0) ? (E_e / E_elec) : 0.0);
    domain_E.table[format("Em{}", idx)] << iodata.DimensionalizeValue(VT::ENERGY, E_m);
    domain_E.table[format("pm{}", idx)] << ((std::abs(E_mag) > 0.0) ? E_m / E_mag : 0.0);
  }

  domain_E.WriteFullTableTrunc();
}

BaseSolver::SurfacesPostPrinter::SurfacesPostPrinter(bool do_measurement, bool root,
                                                     const std::string &post_dir,
                                                     const PostOperator &post_op,
                                                     const std::string &idx_col_name,
                                                     int n_expected_rows)
  : root_{root},
    do_measurement_flux_(do_measurement            //
                         && post_dir.length() > 0  // Valid output dir
                         && post_op.GetSurfacePostOp().flux_surfs.size() > 0  // Has flux
                         ),
    do_measurement_eps_(do_measurement            //
                        && post_dir.length() > 0  // Valid output dir
                        && post_op.GetSurfacePostOp().eps_surfs.size() > 0  // Has eps
    )
{
  if (!root_)
  {
    return;
  }
  using fmt::format;

  if (do_measurement_flux_)
  {
    surface_F = TableWithCSVFile(post_dir + "surface-F.csv");
    surface_F.table.reserve(n_expected_rows,
                            2 * post_op.GetSurfacePostOp().flux_surfs.size() + 1);
    surface_F.table.insert_column(Column("idx", idx_col_name, 0, {}, {}, ""));

    bool has_imaginary = post_op.HasImag();
    for (const auto &[idx, data] : post_op.GetSurfacePostOp().flux_surfs)
    {
      switch (data.type)
      {
        case SurfaceFluxType::ELECTRIC:
          if (has_imaginary)
          {
            surface_F.table.insert_column(format("F_{}_re", idx),
                                          format("Re{{Φ_elec[{}]}} (C)", idx));
            surface_F.table.insert_column(format("F_{}_im", idx),
                                          format("Im{{Φ_elec[{}]}} (C)", idx));
          }
          else
          {
            surface_F.table.insert_column(format("F_{}_re", idx),
                                          format("Φ_elec[{}] (C)", idx));
          }
          break;
        case SurfaceFluxType::MAGNETIC:
          if (has_imaginary)
          {
            surface_F.table.insert_column(format("F_{}_re", idx),
                                          format("Re{{Φ_mag[{}]}} (Wb)", idx));
            surface_F.table.insert_column(format("F_{}_im", idx),
                                          format("Im{{Φ_mag[{}]}} (Wb)", idx));
          }
          else
          {
            surface_F.table.insert_column(format("F_{}_re", idx),
                                          format("Φ_mag[{}] (Wb)", idx));
          }
          break;
        case SurfaceFluxType::POWER:
          surface_F.table.insert_column(format("F_{}_re", idx),
                                        format("Φ_pow[{}] (W)", idx));
          break;
      }
    }
    surface_F.AppendHeader();
  }

  if (do_measurement_eps_)
  {
    surface_Q = TableWithCSVFile(post_dir + "surface-Q.csv");
    surface_Q.table.reserve(n_expected_rows,
                            2 * post_op.GetSurfacePostOp().eps_surfs.size() + 1);
    surface_Q.table.insert_column(Column("idx", idx_col_name, 0, {}, {}, ""));

    for (const auto &[idx, data] : post_op.GetSurfacePostOp().eps_surfs)
    {
      surface_Q.table.insert_column(format("p_{}", idx), format("p_surf[{}]", idx));
      surface_Q.table.insert_column(format("Q_{}", idx), format("Q_surf[{}]", idx));
    }
  }
}

void BaseSolver::SurfacesPostPrinter::AddMeasurementFlux(double idx_value_dimensionful,
                                                         const PostOperator &post_op,
                                                         double E_elec, double E_mag,
                                                         const IoData &iodata)
{
  if (!do_measurement_flux_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Gather
  const bool has_imaginary = post_op.HasImag();
  std::vector<FluxData> flux_data;
  flux_data.reserve(post_op.GetSurfacePostOp().flux_surfs.size());
  for (const auto &[idx, data] : post_op.GetSurfacePostOp().flux_surfs)
  {
    auto Phi = post_op.GetSurfaceFlux(idx);
    switch (data.type)
    {
      case SurfaceFluxType::ELECTRIC:
        Phi *= iodata.DimensionalizeValue(VT::CAPACITANCE, 1.0);
        Phi *= iodata.DimensionalizeValue(VT::VOLTAGE, 1.0);
        break;
      case SurfaceFluxType::MAGNETIC:
        Phi *= iodata.DimensionalizeValue(VT::INDUCTANCE, 1.0);
        Phi *= iodata.DimensionalizeValue(VT::CURRENT, 1.0);
        break;
      case SurfaceFluxType::POWER:
        Phi *= iodata.DimensionalizeValue(VT::POWER, 1.0);
        break;
    }
    flux_data.emplace_back(FluxData{idx, Phi, data.type});
  }

  if (!root_)
  {
    return;
  }

  surface_F.table["idx"] << idx_value_dimensionful;

  for (const auto &[idx, Phi, data_type] : flux_data)
  {
    surface_F.table[format("F_{}_re", idx)] << Phi.real();
    if (has_imaginary &&
        (data_type == SurfaceFluxType::ELECTRIC || data_type == SurfaceFluxType::MAGNETIC))
    {
      surface_F.table[format("F_{}_im", idx)] << Phi.imag();
    }
  }
  surface_F.WriteFullTableTrunc();
}

void BaseSolver::SurfacesPostPrinter::AddMeasurementEps(double idx_value_dimensionful,
                                                        const PostOperator &post_op,
                                                        double E_elec, double E_mag,
                                                        const IoData &iodata)
{
  if (!do_measurement_eps_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  // Gather
  std::vector<EpsData> eps_data;
  eps_data.reserve(post_op.GetSurfacePostOp().eps_surfs.size());
  for (const auto &[idx, data] : post_op.GetSurfacePostOp().eps_surfs)
  {
    double p = post_op.GetInterfaceParticipation(idx, E_elec);
    double tandelta = post_op.GetSurfacePostOp().GetInterfaceLossTangent(idx);
    double Q = (p == 0.0 || tandelta == 0.0) ? mfem::infinity() : 1.0 / (tandelta * p);
    eps_data.emplace_back(EpsData{idx, p, Q});
  }

  if (!root_)
  {
    return;
  }

  surface_Q.table["idx"] << idx_value_dimensionful;
  for (const auto &[idx, p, Q] : eps_data)
  {
    surface_Q.table[format("p_{}", idx)] << p;
    surface_Q.table[format("Q_{}", idx)] << Q;
  }
  surface_Q.WriteFullTableTrunc();
}

void BaseSolver::SurfacesPostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                     const PostOperator &post_op,
                                                     double E_elec, double E_mag,
                                                     const IoData &iodata)
{
  // If surfaces have been specified for postprocessing, compute the corresponding values
  // and write out to disk. The passed in E_elec is the sum of the E-field and lumped
  // capacitor energies, and E_mag is the same for the B-field and lumped inductors.
  AddMeasurementFlux(idx_value_dimensionful, post_op, E_elec, E_mag, iodata);
  AddMeasurementEps(idx_value_dimensionful, post_op, E_elec, E_mag, iodata);
}

BaseSolver::ProbePostPrinter::ProbePostPrinter(bool do_measurement, bool root,
                                               const std::string &post_dir,
                                               const PostOperator &post_op,
                                               const std::string &idx_col_name,
                                               int n_expected_rows)
  : root_{root}, do_measurement_E_{do_measurement}, do_measurement_B_{do_measurement},
    has_imag{post_op.HasImag()}, v_dim{post_op.GetInterpolationOpVDim()}
{
#if defined(MFEM_USE_GSLIB)
  do_measurement_E_ = do_measurement_E_                    //
                      && (post_dir.length() > 0)           // Valid output dir
                      && (post_op.GetProbes().size() > 0)  // Has probes defined
                      && post_op.HasE();                   // Has E fields

  do_measurement_B_ = do_measurement_B_                    //
                      && (post_dir.length() > 0)           // Valid output dir
                      && (post_op.GetProbes().size() > 0)  // Has probes defined
                      && post_op.HasB();                   // Has B fields

  if (!root_ || (!do_measurement_E_ && !do_measurement_B_))
  {
    return;
  }
  using fmt::format;
  int scale_col = (has_imag ? 2 : 1) * v_dim;
  auto dim_labeler = [](int i) -> std::string
  {
    switch (i)
    {
      // Note: Zero-based indexing here
      case 0:
        return "x";
      case 1:
        return "y";
      case 2:
        return "z";
      default:
        return format("d{}", i);
    }
  };

  if (do_measurement_E_)
  {
    probe_E = TableWithCSVFile(post_dir + "probe-E.csv");
    probe_E.table.reserve(n_expected_rows, scale_col * post_op.GetProbes().size());
    probe_E.table.insert_column(Column("idx", idx_col_name, 0, {}, {}, ""));

    for (const auto &idx : post_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if (has_imag)
        {
          probe_E.table.insert_column(
              format("E{}_{}_re", idx, i_dim),
              format("Re{{E_{}[{}]}} (V/m)", dim_labeler(i_dim), idx));
          probe_E.table.insert_column(
              format("E{}_{}_im", idx, i_dim),
              format("Im{{E_{}[{}]}} (V/m)", dim_labeler(i_dim), idx));
        }
        else
        {
          probe_E.table.insert_column(format("E{}_{}_re", idx, i_dim),
                                      format("E_{}[{}] (V/m)", dim_labeler(i_dim), idx));
        }
      }
    }
    probe_E.AppendHeader();
  }

  if (do_measurement_B_)
  {
    probe_B = TableWithCSVFile(post_dir + "probe-B.csv");
    probe_B.table.reserve(n_expected_rows, scale_col * post_op.GetProbes().size());
    probe_B.table.insert_column(Column("idx", idx_col_name, 0, {}, {}, ""));

    for (const auto &idx : post_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if (has_imag)
        {
          probe_B.table.insert_column(
              format("B{}_{}_re", idx, i_dim),
              format("Re{{B_{}[{}]}} (Wb/m²)", dim_labeler(i_dim), idx));
          probe_B.table.insert_column(
              format("B{}_{}_im", idx, i_dim),
              format("Im{{B_{}[{}]}} (Wb/m²)", dim_labeler(i_dim), idx));
        }
        else
        {
          probe_B.table.insert_column(format("B{}_{}_re", idx, i_dim),
                                      format("B_{}[{}] (Wb/m²)", dim_labeler(i_dim), idx));
        }
      }
    }
    probe_B.AppendHeader();
  }
#endif
}

void BaseSolver::ProbePostPrinter::AddMeasurementE(double idx_value_dimensionful,
                                                   const PostOperator &post_op,
                                                   const IoData &iodata)
{
  if (!do_measurement_E_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  auto probe_field = post_op.ProbeEField();
  MFEM_VERIFY(probe_field.size() == v_dim * post_op.GetProbes().size(),
              format("Size mismatch: expect vector field to ahve size {} * {} = {}; got {}",
                     v_dim, post_op.GetProbes().size(), v_dim * post_op.GetProbes().size(),
                     probe_field.size()))

  if (!root_)
  {
    return;
  }

  probe_E.table["idx"] << idx_value_dimensionful;
  size_t i = 0;
  for (const auto &idx : post_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = iodata.DimensionalizeValue(VT::FIELD_E, probe_field[i * v_dim + i_dim]);
      probe_E.table[format("E{}_{}_re", idx, i_dim)] << val.real();
      if (has_imag)
      {
        probe_E.table[format("E{}_{}_im", idx, i_dim)] << val.imag();
      }
    }
    i++;
  }
  probe_E.WriteFullTableTrunc();
}

void BaseSolver::ProbePostPrinter::AddMeasurementB(double idx_value_dimensionful,
                                                   const PostOperator &post_op,
                                                   const IoData &iodata)
{
  if (!do_measurement_B_)
  {
    return;
  }
  using VT = IoData::ValueType;
  using fmt::format;

  auto probe_field = post_op.ProbeBField();
  MFEM_VERIFY(probe_field.size() == v_dim * post_op.GetProbes().size(),
              format("Size mismatch: expect vector field to ahve size {} * {} = {}; got {}",
                     v_dim, post_op.GetProbes().size(), v_dim * post_op.GetProbes().size(),
                     probe_field.size()))

  if (!root_)
  {
    return;
  }

  probe_B.table["idx"] << idx_value_dimensionful;
  size_t i = 0;
  for (const auto &idx : post_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = iodata.DimensionalizeValue(VT::FIELD_B, probe_field[i * v_dim + i_dim]);
      probe_B.table[format("B{}_{}_re", idx, i_dim)] << val.real();
      if (has_imag)
      {
        probe_B.table[format("B{}_{}_im", idx, i_dim)] << val.imag();
      }
    }
    i++;
  }
  probe_B.WriteFullTableTrunc();
}

void BaseSolver::ProbePostPrinter::AddMeasurement(double idx_value_dimensionful,
                                                  const PostOperator &post_op,
                                                  const IoData &iodata)
{
#if defined(MFEM_USE_GSLIB)
  AddMeasurementE(idx_value_dimensionful, post_op, iodata);
  AddMeasurementB(idx_value_dimensionful, post_op, iodata);
#endif
}

BaseSolver::ErrorIndicatorPostPrinter::ErrorIndicatorPostPrinter(
    bool do_measurement, bool root, const std::string &post_dir)
  : root_{root}, do_measurement_{
                     do_measurement            //
                     && post_dir.length() > 0  // Valid output dir
                 }
{
  if (!do_measurement_ || !root_)
  {
    return;
  }

  error_indicator = TableWithCSVFile(post_dir + "error-indicators.csv");
  error_indicator.table.reserve(1, 4);

  error_indicator.table.insert_column("norm", "Norm");
  error_indicator.table.insert_column("min", "Minimum");
  error_indicator.table.insert_column("max", "Maximum");
  error_indicator.table.insert_column("mean", "Mean");
}

void BaseSolver::ErrorIndicatorPostPrinter::PrintIndicatorStatistics(
    const PostOperator &post_op, const ErrorIndicator &indicator)
{
  if (!do_measurement_)
  {
    return;
  }

  // MPI Gather
  MPI_Comm comm = post_op.GetComm();
  std::array<double, 4> data = {indicator.Norml2(comm), indicator.Min(comm),
                                indicator.Max(comm), indicator.Mean(comm)};

  if (!root_)
  {
    return;
  }

  error_indicator.table["norm"] << data[0];
  error_indicator.table["min"] << data[1];
  error_indicator.table["max"] << data[2];
  error_indicator.table["mean"] << data[3];

  error_indicator.WriteFullTableTrunc();
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
