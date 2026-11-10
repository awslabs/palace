// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_COLLECTIVE_ENVIRONMENT_HPP
#define PALACE_UTILS_COLLECTIVE_ENVIRONMENT_HPP

#include <cstdlib>
#include <optional>
#include <string>
#include "utils/communication.hpp"

namespace palace
{

// Lane fix-response-memory round 2 C10 (decision 1471 (2); OI-MPIRUN-ENV-NOT-FORWARDED-TO-
// REMOTE-NODES): the value of the environment variable `name` AS THE ROOT RANK OF `comm`
// READS IT, on every rank of `comm`. Collective (a broadcast of the length, then of the
// bytes): every rank must call it at the same point. A gate whose path enters a collective
// reads its variable through this, so that a per-rank environment mismatch (a variable
// exported by the job script reaches only the head node's ranks under mpirun without -x)
// can never send the ranks down different collective paths.
inline std::optional<std::string> CollectiveGetenv(MPI_Comm comm, const char *name)
{
  int length = -1;
  std::string value;
  if (Mpi::Root(comm))
  {
    if (const char *local = std::getenv(name))
    {
      value = local;
      length = static_cast<int>(value.size());
    }
  }
  Mpi::Broadcast(1, &length, 0, comm);
  if (length < 0)
  {
    return std::nullopt;
  }
  value.resize(static_cast<std::size_t>(length));
  if (length > 0)
  {
    Mpi::Broadcast(length, value.data(), 0, comm);
  }
  return value;
}

// The root's "<name>=1" on every rank (collective).
inline bool CollectiveEnvironmentFlag(MPI_Comm comm, const char *name)
{
  const auto value = CollectiveGetenv(comm, name);
  return value && *value == "1";
}

}  // namespace palace

#endif  // PALACE_UTILS_COLLECTIVE_ENVIRONMENT_HPP
