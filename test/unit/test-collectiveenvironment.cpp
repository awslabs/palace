// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cstdlib>
#include <optional>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "models/surfaceresponseoperator.hpp"
#include "utils/collectiveenvironment.hpp"
#include "utils/communication.hpp"

using namespace palace;

namespace
{

// Set (value != nullptr) or unset the variable on THIS rank only.
void SetLocal(const char *name, const char *value)
{
  if (value)
  {
    setenv(name, value, 1);
  }
  else
  {
    unsetenv(name);
  }
}

}  // namespace

// Lane fix-response-memory round 2 C10 (decision 1471 (2); OI-MPIRUN-ENV-NOT-FORWARDED-TO-
// REMOTE-NODES): a variable exported by a PBS job script reaches only the head node's ranks
// under mpirun without -x. A gate whose path enters a collective must give every rank the
// ROOT's value; here the ranks see DIFFERENT values and every rank must read the root's.
TEST_CASE("Collective environment gates read the root's value on every rank",
          "[collectiveenvironment][memfix][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const bool root = Mpi::Root(comm);

  // The census gate (PBS 104703 deadlocked on a per-rank read): set on the root only, then
  // on the other ranks only.
  SetLocal("PALACE_RESPONSE_MEMORY_CENSUS", root ? "1" : nullptr);
  CHECK(ResponseMemoryCensusEnabled(comm));
  CHECK(CollectiveEnvironmentFlag(comm, "PALACE_RESPONSE_MEMORY_CENSUS"));
  SetLocal("PALACE_RESPONSE_MEMORY_CENSUS", root ? nullptr : "1");
  CHECK(!ResponseMemoryCensusEnabled(comm));
  SetLocal("PALACE_RESPONSE_MEMORY_CENSUS", root ? "0" : "1");
  CHECK(!ResponseMemoryCensusEnabled(comm));
  SetLocal("PALACE_RESPONSE_MEMORY_CENSUS", nullptr);
  CHECK(!ResponseMemoryCensusEnabled(comm));

  // The containment block size: every rank walks the root's blocks.
  SetLocal("PALACE_RESPONSE_CONTAINMENT_BLOCK_POINTS", root ? "65536" : "7");
  CHECK(DomainBoundaryContainmentBlockPoints(comm) == 65536);
  SetLocal("PALACE_RESPONSE_CONTAINMENT_BLOCK_POINTS", root ? nullptr : "7");
  CHECK(DomainBoundaryContainmentBlockPoints(comm) == (std::size_t(1) << 23));
  SetLocal("PALACE_RESPONSE_CONTAINMENT_BLOCK_POINTS", nullptr);

  // The raw read: the root's bytes (incl. an empty value) or nothing, on every rank.
  SetLocal("PALACE_RESPONSE_COLLECTIVE_TEST", root ? "the root's value" : "another value");
  const auto value = CollectiveGetenv(comm, "PALACE_RESPONSE_COLLECTIVE_TEST");
  REQUIRE(value);
  CHECK(*value == "the root's value");
  SetLocal("PALACE_RESPONSE_COLLECTIVE_TEST", root ? "" : "another value");
  const auto empty = CollectiveGetenv(comm, "PALACE_RESPONSE_COLLECTIVE_TEST");
  REQUIRE(empty);
  CHECK(empty->empty());
  SetLocal("PALACE_RESPONSE_COLLECTIVE_TEST", root ? nullptr : "another value");
  CHECK(!CollectiveGetenv(comm, "PALACE_RESPONSE_COLLECTIVE_TEST"));
  SetLocal("PALACE_RESPONSE_COLLECTIVE_TEST", nullptr);
}
