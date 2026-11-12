// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/communication.hpp"
#include "utils/perfstats.hpp"

namespace palace::perf
{

namespace detail
{

namespace
{

struct CommTallies
{
  int reduce, barrier, broadcast;
};

// Communication tallies, registered together and with all their rows at once: which rows
// communicate can depend on the rank (collectives on sub-communicators). Barriers and
// broadcasts are rare outside of setup: they are not split by row.
const CommTallies &Tallies()
{
  static const CommTallies t = {RegisterPhaseTally("Reduce"), RegisterTally("Barrier"),
                                RegisterTally("Broadcast")};
  return t;
}

}  // namespace

void AddCollective(Collective kind, double seconds)
{
  if (registry.suspended)
  {
    return;
  }
  const auto &t = Tallies();
  switch (kind)
  {
    case Collective::REDUCE:
      AddTally(t.reduce, seconds);
      break;
    case Collective::BARRIER:
      AddTally(t.barrier, seconds);
      break;
    case Collective::BROADCAST:
      AddTally(t.broadcast, seconds);
      break;
  }
}

}  // namespace detail

void RegisterCommunicationTallies()
{
  detail::Tallies();
}

}  // namespace palace::perf
