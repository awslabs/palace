// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_session.hpp>
#include "fem/libceed/ceed.hpp"
#include "utils/communication.hpp"

using namespace palace;

// Global test options configurable from command line.
int benchmark_ref_levels = 0;
int benchmark_order = 4;
bool benchmark_assemble_q_data = false;
bool benchmark_no_fa = false;
bool benchmark_no_mfem_pa = false;

int main(int argc, char *argv[])
{
  // Initialize MPI.
  Mpi::Init(argc, argv);

  // See https://github.com/catchorg/Catch2/blob/devel/docs/own-main.md.
  Catch::Session session;

  // Extra command line arguments.
  std::string device_str("cpu");          // MFEM device
  std::string ceed_backend("/cpu/self");  // libCEED backend

  // Build a new parser on top of Catch2's.
  using namespace Catch::Clara;
  auto cli = session.cli() |
             Opt(device_str, "device")["--device"]("MFEM device (default: \"cpu\")") |
             Opt(ceed_backend,
                 "backend")["--backend"]("libCEED backend (default: \"/cpu/self\")") |
             Opt(benchmark_ref_levels, "levels")["--benchmark-ref-levels"](
                 "Levels of uniform mesh refinement for benchmarks (default: 0)") |
             Opt(benchmark_order, "order")["--benchmark-order"](
                 "Element order for benchmarks (default: 4)") |
             Opt(benchmark_assemble_q_data)["--benchmark-assemble-q-data"](
                 "Assemble quadrature data for benchmark operators") |
             Opt(benchmark_no_fa)["--benchmark-skip-full-assembly"](
                 "Skip full assembly tests in benchmarks") |
             Opt(benchmark_no_mfem_pa)["--benchmark-skip-mfem-partial-assembly"](
                 "Skip MFEM partial assembly tests in benchmarks");

  // Now pass the new composite back to Catch2 so it uses that.
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line.
  int result = session.applyCommandLine(argc, argv);
  if (result != 0)
  {
    return result;
  }

  // Run the tests.
  mfem::Device device(device_str.c_str());
  ceed::Initialize(ceed_backend.c_str(), PALACE_LIBCEED_JIT_SOURCE_DIR);
  std::ostringstream resource(std::stringstream::out);
  device.Print(resource);
  resource << "libCEED backend: " << ceed::Print();
  Mpi::Print("{}\n", resource.str());
  result = session.run();
  ceed::Finalize();

  return result;
}
