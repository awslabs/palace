// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
#define PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP

#include "fem/errorindicator.hpp"

namespace palace
{

class IoData;
class SpaceOperator;

// The frequency sweep of a driven problem with substructuring (config["Solver"]
// ["Substructuring"]), uniform or adaptive, offline or online.
ErrorIndicator SweepSubstructured(const IoData &iodata, SpaceOperator &space_op,
                                  bool adaptive);

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
