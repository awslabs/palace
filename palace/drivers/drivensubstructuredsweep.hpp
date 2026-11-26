// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
#define PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP

#include <memory>

namespace palace
{

class EnvironmentOperator;
class IoData;
class SpaceOperator;

// The environment of the region of an offline exact sweep of a driven problem with
// substructuring (config["Solver"]["Substructuring"]): an exact model at the frequencies of
// the sweep, which condenses the environment at each one when the region's sweep reaches it
// (appending the record to the saved model, if any). It holds the factorization of the
// environment, on the operators of the full problem, which must outlive it.
std::unique_ptr<EnvironmentOperator> ExactEnvironmentModel(const IoData &iodata,
                                                           SpaceOperator &space_op);

// The offline adaptive sampling of the environment's rational model (saved with
// "ModelFile"), as the environment of the region of an adaptive sweep, which then reduces
// the region (see DrivenSolver::SweepAdaptive).
std::unique_ptr<EnvironmentOperator> SampleEnvironmentModel(const IoData &iodata,
                                                            SpaceOperator &space_op);

}  // namespace palace

#endif  // PALACE_DRIVERS_DRIVEN_SUBSTRUCTURED_SWEEP_HPP
