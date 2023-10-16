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
#include "linalg/errorestimator.hpp"
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

std::string BaseSolver::IterationPostDir(int iter) const
{
  std::string dir = fmt::format("{}Adapt{:0>3d}/", post_dir, iter);
  // Create directory for output.
  if (root && !std::filesystem::exists(dir))
  {
    std::filesystem::create_directories(dir);
  }
  return dir;
}

namespace
{

// Helper function that returns an array of indices corresponding to marked elements.
mfem::Array<int> MarkedElements(double threshold, const mfem::Vector &e)
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

// Helper function responsible for rebalancing the mesh, and optionally writing meshes from
// the intermediate stages to disk.
void RebalanceMesh(std::unique_ptr<mfem::ParMesh> &mesh, const IoData &iodata,
                   std::string output_dir)
{
  const auto &param = iodata.model.refinement.adaptation;
  auto comm = mesh->GetComm();
  const int width = 1 + static_cast<int>(std::log10(Mpi::Size(comm) - 1));
  if (output_dir.back() != '/')
  {
    output_dir += '/';
  }
  std::string serial_mesh_filename;
  if (param.write_serial_mesh)
  {
    serial_mesh_filename += (output_dir + "serial.mesh");
    if (Mpi::Size(comm) == 1)
    {
      std::ofstream serial(serial_mesh_filename);
      iodata.DimensionalizeMesh(*mesh);
      mesh->Mesh::Print(serial);
      iodata.NondimensionalizeMesh(*mesh);
    }
  }
  if (param.write_pre_balance_mesh)
  {
    std::ofstream pfile(
        mfem::MakeParFilename(output_dir + "prebalance.", Mpi::Rank(comm), ".mesh", width));
    mesh->SetPrintShared(false);  // Do not mark processor boundaries in the save
    BlockTimer bt(Timer::IO);
    iodata.DimensionalizeMesh(*mesh);
    mesh->ParPrint(pfile);
    iodata.NondimensionalizeMesh(*mesh);
  }

  // If there is more than one processor, may perform rebalancing.
  if (Mpi::Size(comm) > 1)
  {
    BlockTimer bt(Timer::REBALANCE);
    mesh->ExchangeFaceNbrData();
    int min_elem, max_elem;
    min_elem = max_elem = mesh->GetNE();
    Mpi::GlobalMin(1, &min_elem, comm);
    Mpi::GlobalMax(1, &max_elem, comm);
    const double ratio = double(max_elem) / min_elem;
    Mpi::Print("Min Elem per processor: {}, Max Elem per processor: {}, Ratio: {:.3e}\n",
               min_elem, max_elem, double(max_elem) / min_elem);
    if (ratio > param.maximum_imbalance)
    {
      if (mesh->Nonconforming() && param.write_serial_mesh)
      {
        // Do not need to duplicate the mesh, as rebalancing will undo this.
        BlockTimer bt_post(Timer::POSTPRO);
        mfem::Array<int> serial_partition(mesh->GetNE());
        serial_partition = 0;
        mesh->Rebalance(serial_partition);
        BlockTimer bt_io(Timer::IO);
        if (Mpi::Root(comm))
        {
          std::ofstream serial(serial_mesh_filename);
          iodata.DimensionalizeMesh(*mesh);
          mesh->Mesh::Print(serial);
          iodata.NondimensionalizeMesh(*mesh);
        }
        MPI_Barrier(comm);
      }

      Mpi::Print("Ratio {:.3f} exceeds maximum allowed value {}: Performing rebalancing.\n",
                 ratio, param.maximum_imbalance);
      if (mesh->Nonconforming())
      {
        mesh->Rebalance();
      }
      else
      {
        // Without access to a refinement tree, partitioning must be done on the
        // root processor and then redistributed.
        mesh::RebalanceConformalMesh(mesh, iodata, serial_mesh_filename);
      }
      if (param.write_post_balance_mesh)
      {
        BlockTimer bt_io(Timer::IO);
        std::ofstream pfile(mfem::MakeParFilename(output_dir + "postbalance.",
                                                  Mpi::Rank(comm), ".mesh", width));
        mesh->SetPrintShared(false);  // Do not mark processor boundaries in the save
        iodata.DimensionalizeMesh(*mesh);
        mesh->ParPrint(pfile);
        iodata.NondimensionalizeMesh(*mesh);
      }
    }
    else if (param.write_serial_mesh)
    {
      if (mesh->Nonconforming())
      {
        // Given no rebalancing will be done, need to handle the serial write more
        // carefully. For NC, this requires creating a duplicate mesh.
        BlockTimer bt_post(Timer::POSTPRO);
        mfem::ParMesh smesh(*mesh);
        mfem::Array<int> serial_partition(mesh->GetNE());
        smesh.FinalizeTopology();
        smesh.Finalize();
        smesh.ExchangeFaceNbrData();
        serial_partition = 0;
        smesh.Rebalance(serial_partition);
        BlockTimer bt_io(Timer::IO);
        if (Mpi::Root(comm))
        {
          std::ofstream serial(serial_mesh_filename);
          iodata.DimensionalizeMesh(smesh);
          smesh.Mesh::Print(serial);
        }
      }
      else
      {
        std::ofstream serial(serial_mesh_filename);
        serial.precision(std::numeric_limits<double>::max_digits10);
        iodata.DimensionalizeMesh(*mesh);
        mesh->PrintAsSerial(serial);
        iodata.NondimensionalizeMesh(*mesh);
      }
    }
    mesh->ExchangeFaceNbrData();
  }
}

}  // namespace

ErrorIndicator
BaseSolver::SolveEstimateMarkRefine(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const
{
  auto comm = mesh.back()->GetComm();
  // Helper to save off postprocess data.
  auto SavePostProcess = [&](int iter)
  {
    Mpi::Barrier(comm);
    BlockTimer bt_io(Timer::IO);
    if (Mpi::Root(comm))
    {
      // Create a subfolder for the results of this adaptation.
      namespace fs = std::filesystem;
      const auto out_dir = fs::path(IterationPostDir(iter));
      constexpr auto options =
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
    Mpi::Barrier(comm);
  };

  const auto &param = iodata.model.refinement.adaptation;
  const bool use_amr = param.max_its > 0;
  const bool use_coarsening = mesh.back()->Nonconforming() && param.use_coarsening;
  if (use_amr && mesh.size() > 1)
  {
    Mpi::Print("{}\n", "Flattening mesh sequence: AMR will start from the final mesh in "
                       "the refinement sequence.");
    mesh.erase(mesh.begin(), mesh.end() - 1);
    constexpr bool refine = true, fix_orientation = true;
    mesh.back()->Finalize(refine, fix_orientation);
  }

  int iter = 0;
  auto indicators_and_ntdof = Solve(mesh);
  const auto &indicators = indicators_and_ntdof.first;
  const auto &ntdof = indicators_and_ntdof.second;
  if (use_amr)
  {
    const auto is_transient = dynamic_cast<const TransientSolver *>(this) != nullptr;
    if (is_transient)
    {
      MFEM_ABORT("AMR is not currently supported for transient simulations!");
    }
    else
    {
      Mpi::Print("\nAdaptive Mesh Refinement Parameters:\n");
      Mpi::Print("MaxIter: {}, Tolerance: {:.3e}\n\n", param.max_its, param.tolerance);
      SavePostProcess(iter);  // Save the initial solution
    }
  }

  // Collection of all tests that might exhaust resources.
  auto exhausted_resources = [&]()
  {
    bool ret = false;
    // Run out of DOFs, and coarsening isn't allowed.
    ret |= (!use_coarsening && ntdof > param.dof_limit);
    // Run out of iterations.
    ret |= iter >= param.max_its;
    return ret;
  };

  while (indicators.Norml2(comm) > param.tolerance && !exhausted_resources())
  {
    BlockTimer bt_adapt(Timer::ADAPT);
    Mpi::Print("Adaptation iteration {}: Initial Error Indicator: {:.3e}, DOF: {}, DOF "
               "Limit: {}\n",
               ++iter, indicators.Norml2(comm), ntdof, param.dof_limit);
    if (ntdof < param.dof_limit)
    {
      // Mark.
      const auto threshold =
          utils::ComputeDorflerThreshold(param.update_fraction, indicators.Local());
      const auto marked_elements = MarkedElements(threshold, indicators.Local());

      // Refine.
      const auto initial_elem_count = mesh.back()->GetGlobalNE();
      mesh.back()->GeneralRefinement(marked_elements, -1, param.max_nc_levels);
      const auto final_elem_count = mesh.back()->GetGlobalNE();
      Mpi::Print("Mesh refinement added {} elements. Initial: {}, Final: {}\n",
                 final_elem_count - initial_elem_count, initial_elem_count,
                 final_elem_count);
    }
    else if (use_coarsening)
    {
      // Perform a Dörfler style marking looking for the largest number of derefinement
      // opportunities to represent a fraction of the derefinable error.
      const auto &derefinement_table = mesh.back()->pncmesh->GetDerefinementTable();
      mfem::Vector coarse_error(derefinement_table.Size());
      mfem::Array<int> row;
      for (int i = 0; i < derefinement_table.Size(); i++)
      {
        // Sum the error for all sub elements that can be combined.
        derefinement_table.GetRow(i, row);
        coarse_error[i] =
            std::sqrt(std::accumulate(row.begin(), row.end(), 0.0,
                                      [&indicators](double s, int i) {
                                        return s += std::pow(indicators.Local()[i], 2.0);
                                      }));
      }

      // Given the coarse errors, we use the Dörfler marking strategy to identify the
      // smallest set of original elements that make up (1 - θ) of the total error. The
      // complement of this set is then the largest number of elements that make up θ of the
      // total error.
      const double threshold =
          utils::ComputeDorflerThreshold(1 - param.update_fraction, coarse_error);

      const auto initial_elem_count = mesh.back()->GetGlobalNE();
      constexpr int aggregate_operation = 3;  // sum of squares
      mesh.back()->DerefineByError(indicators.Local(), threshold, param.max_nc_levels,
                                   aggregate_operation);
      const auto final_elem_count = mesh.back()->GetGlobalNE();
      Mpi::Print("Mesh coarsening removed {} elements. Initial: {}, Final: {}\n",
                 initial_elem_count - final_elem_count, initial_elem_count,
                 final_elem_count);
    }
    RebalanceMesh(mesh.back(), iodata, post_dir);

    // Solve + estimate.
    indicators_and_ntdof = Solve(mesh);

    // Optionally save solution off.
    if (param.save_step > 0 && iter % param.save_step == 0)
    {
      SavePostProcess(iter);
    }
  }
  Mpi::Print("\nFinal Error Indicator: {:.3e}, DOF: {}\n", indicators.Norml2(comm), ntdof);
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
    for (int i = Timer::INIT; i < Timer::NUMTIMINGS; i++)
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

void BaseSolver::PostprocessFields(const PostOperator &postop, int step, double time,
                                   const ErrorIndicator *indicator) const
{
  // Save the computed fields in parallel in format for viewing with ParaView.
  BlockTimer bt(Timer::IO);
  if (post_dir.length() == 0)
  {
    Mpi::Warning(postop.GetComm(),
                 "No file specified under [\"Problem\"][\"Output\"]!\nSkipping saving of "
                 "fields to disk!\n");
    return;
  }
  postop.WriteFields(step, time, indicator);
  Mpi::Barrier(postop.GetComm());
}

void BaseSolver::PostprocessErrorIndicator(const PostOperator &postop,
                                           const ErrorIndicator &indicator) const
{
  // Write the indicator statistics.
  if (post_dir.length() == 0)
  {
    return;
  }
  MPI_Comm comm = postop.GetComm();
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
}

template void BaseSolver::SaveMetadata<KspSolver>(const KspSolver &) const;
template void BaseSolver::SaveMetadata<ComplexKspSolver>(const ComplexKspSolver &) const;

}  // namespace palace
