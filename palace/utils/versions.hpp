// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_VERSIONS_HPP
#define PALACE_UTILS_VERSIONS_HPP

#include <string>
#include <utility>
#include <vector>

namespace palace
{

// Return (name, version) pairs for the dependencies Palace was built with, as reported by
// each library's runtime version query or, where none exists, its headers.
std::vector<std::pair<std::string, std::string>> GetDependencyVersions();

}  // namespace palace

#endif  // PALACE_UTILS_VERSIONS_HPP
