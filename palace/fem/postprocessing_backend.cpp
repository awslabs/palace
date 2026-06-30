// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "fem/postprocessing_backend.hpp"

#include <cstdlib>

namespace palace
{

namespace fem
{

bool LibceedPostprocessingEnabled()
{
  static const bool enabled = !std::getenv("PALACE_LEGACY_SURFACE_POSTPRO");
  return enabled;
}

}  // namespace fem

}  // namespace palace
