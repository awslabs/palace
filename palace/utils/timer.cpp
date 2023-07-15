// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils/timer.hpp"
#include <string>
#include <vector>

namespace palace
{

const std::vector<std::string> Timer::descriptions({"Initialization",
                                                    "Operator Construction", "Solve",
                                                    "Postprocessing", "Disk IO", "Total"});

Timer TimedBlock::timer;
std::vector<int> TimedBlock::stack;

}  // namespace palace