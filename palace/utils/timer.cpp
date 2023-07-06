// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <string>
#include <vector>
#include "utils/timer.hpp"

namespace palace
{

const std::vector<std::string> Timer::descriptions ({
  "Initialization",
  "Operator Construction",
  "Solve",
  "Postprocessing",
  "Disk IO",
  "Total"
});

}