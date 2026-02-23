// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_ADAPT_DRIVER_HPP
#define PALACE_DRIVERS_ADAPT_DRIVER_HPP

#include <string>
#include <mpi.h>

namespace palace
{

class IoData;

//
// Driver for invoking an external mesh adaptation tool (Metris or refine) as a subprocess.
// The adapted mesh is written to disk and its path returned for subsequent loading.
//
class ExternalAdaptationDriver
{
public:
  // Run an external adaptation tool on a mesh + metric pair.
  //
  // Only MPI rank 0 executes the subprocess; all other ranks wait at a barrier. The
  // adapted mesh file path is broadcast to all ranks before returning.
  //
  // Parameters:
  //   mesh_file     - Path to the input mesh in .meshb format.
  //   metric_file   - Path to the metric tensor field in .solb format.
  //   geometry_file - Path to the CAD geometry file (.egads). May be empty if no geometry
  //                   projection is needed.
  //   iodata        - Configuration data (reads external adapter settings from
  //                   iodata.model.refinement).
  //   comm          - MPI communicator.
  //
  // Returns the path to the adapted mesh file (.meshb).
  static std::string Adapt(const std::string &mesh_file,
                           const std::string &metric_file,
                           const std::string &geometry_file,
                           const IoData &iodata,
                           MPI_Comm comm);
};

}  // namespace palace

#endif  // PALACE_DRIVERS_ADAPT_DRIVER_HPP
