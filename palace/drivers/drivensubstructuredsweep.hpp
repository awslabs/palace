// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
#define PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP

#include <memory>
#include "fem/errorindicator.hpp"

namespace palace
{

class EnvironmentOperator;
class IoData;
class SpaceOperator;

// The offline exact sweep of a driven problem with substructuring (config["Solver"]
// ["Substructuring"]): at each frequency, condense the environment (saving the record of
// the model with "ModelFile"), solve the region against it and expand the environment, for
// the outputs of the whole device.
ErrorIndicator SweepSubstructured(const IoData &iodata, SpaceOperator &space_op);

// The offline adaptive sampling of the environment's rational model (saved with
// "ModelFile"), as the environment of the region of an adaptive sweep, which then reduces
// the region (see DrivenSolver::SweepAdaptive).
std::unique_ptr<EnvironmentOperator> SampleEnvironmentModel(const IoData &iodata,
                                                            SpaceOperator &space_op);

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
