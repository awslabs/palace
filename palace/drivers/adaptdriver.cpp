// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "adaptdriver.hpp"

#include <cstdlib>
#include <string>
#include <fmt/format.h>
#include <mfem.hpp>
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/iodata.hpp"

namespace palace
{

std::string ExternalAdaptationDriver::Adapt(const std::string &mesh_file,
                                            const std::string &metric_file,
                                            const std::string &geometry_file,
                                            const IoData &iodata,
                                            MPI_Comm comm)
{
  const auto &refinement = iodata.model.refinement;
  const auto &adapter = refinement.external_adapter;
  MFEM_VERIFY(!adapter.empty(),
              "ExternalAdaptationDriver::Adapt called with no external adapter configured!");

  // Derive the output path from the input mesh path by replacing the extension.
  const fs::path input_path(mesh_file);
  const fs::path output_path =
      input_path.parent_path() / (input_path.stem().string() + "_adapted.meshb");
  const std::string output_file = output_path.string();

  // Only rank 0 runs the external subprocess.
  int status = 0;
  if (Mpi::Root(comm))
  {
    std::string cmd;
    const std::string &binary = refinement.adapter_path.empty() ? adapter
                                                                : refinement.adapter_path;

    if (adapter == "metris")
    {
      // metris --in mesh.meshb --met metric.solb [--cad model.egads]
      //        --adapt N --tardeg D --out adapted.meshb
      cmd = fmt::format("{} --in {} --met {}", binary, mesh_file, metric_file);
      if (!geometry_file.empty())
      {
        cmd += fmt::format(" --cad {}", geometry_file);
      }
      cmd += fmt::format(" --adapt {} --tardeg {} --out {}",
                          refinement.adapt_iterations, refinement.target_degree,
                          output_file);
    }
    else if (adapter == "refine")
    {
      // refine mesh.meshb [geometry.egads] -m metric.solb -o adapted.meshb
      cmd = fmt::format("{} {}", binary, mesh_file);
      if (!geometry_file.empty())
      {
        cmd += fmt::format(" {}", geometry_file);
      }
      cmd += fmt::format(" -m {} -o {}", metric_file, output_file);
    }
    else
    {
      MFEM_ABORT("Unknown external adapter \"" << adapter
                                                << "\": expected \"metris\" or \"refine\"!");
    }

    Mpi::Print(" Running external adaptation: {}\n", cmd);
    status = std::system(cmd.c_str());
    if (status != 0)
    {
      Mpi::Warning("External adapter returned non-zero exit status {}!\n", status);
    }

    // Verify the output file was created.
    MFEM_VERIFY(fs::exists(output_path),
                "External adapter did not produce output mesh \"" << output_file << "\"!");
  }

  // Broadcast the exit status so all ranks can abort together if needed.
  Mpi::Broadcast(1, &status, 0, comm);
  MFEM_VERIFY(status == 0,
              "External mesh adaptation failed with exit status " << status << "!");

  // Broadcast the output file path (all ranks need it to load the adapted mesh).
  // The path length is fixed at compile time on rank 0, so broadcast its length first.
  int path_len = static_cast<int>(output_file.size());
  Mpi::Broadcast(1, &path_len, 0, comm);
  std::string result(path_len, '\0');
  if (Mpi::Root(comm))
  {
    result = output_file;
  }
  MPI_Bcast(result.data(), path_len, MPI_CHAR, 0, comm);

  return result;
}

}  // namespace palace
