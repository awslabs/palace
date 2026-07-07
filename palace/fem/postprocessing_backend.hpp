// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_POSTPROCESSING_BACKEND_HPP
#define PALACE_FEM_POSTPROCESSING_BACKEND_HPP

namespace palace
{

namespace fem
{

// Global development/diagnostic gate for libCEED-backed postprocessing paths. Legacy
// coefficient paths remain available as correctness oracles and selected fallbacks during
// this refactor, but supported production paths should prefer the libCEED backend.
bool LibceedPostprocessingEnabled();

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_POSTPROCESSING_BACKEND_HPP
