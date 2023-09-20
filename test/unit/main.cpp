// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <string>
#include <mfem.hpp>
#include <catch2/catch_session.hpp>
#include "fem/libceed/utils.hpp"

using namespace palace;

int main(int argc, char *argv[])
{
  // See https://github.com/catchorg/Catch2/blob/devel/docs/own-main.md.
  Catch::Session session;

  // Extra command line arguments.
  std::string device_str("cpu");          // MFEM device
  std::string ceed_backend("/cpu/self");  // libCEED backend

  // Build a new parser on top of Catch2's.
  using namespace Catch::Clara;
  auto cli =
      session.cli() |
      Opt(device_str, "device")["--device"]("MFEM device (default: \"cpu\")") |
      Opt(ceed_backend, "backend")["--backend"]("libCEED backend (default: \"/cpu/self\")");

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
  result = session.run();
  ceed::Finalize();

  return result;
}
