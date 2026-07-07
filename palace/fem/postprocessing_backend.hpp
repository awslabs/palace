// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_POSTPROCESSING_BACKEND_HPP
#define PALACE_FEM_POSTPROCESSING_BACKEND_HPP

namespace palace
{

namespace fem
{

// Central query for libCEED-backed postprocessing availability. Supported production
// postprocessing paths use the libCEED backend directly; coefficient implementations are
// retained only as reference semantics in tests and for explicitly unsupported features.
bool LibceedPostprocessingEnabled();

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_POSTPROCESSING_BACKEND_HPP
