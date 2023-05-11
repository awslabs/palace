// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_FILESYSTEM_HPP
#define PALACE_UTILS_FILESYSTEM_HPP

#if defined(__cpp_lib_filesystem) || defined(__has_include) && __has_include(<filesystem>)
#include <filesystem>
#elif defined(__cpp_lib_experimental_filesystem) || \
    defined(__has_include) && __has_include(<experimental/filesystem>)
// clang-format off
#include <experimental/filesystem>
namespace std { namespace filesystem = experimental::filesystem; }
// clang-format on
#else
#error "Could not find system header <filesystem> or <experimental/filesystem>!"
#endif

#endif  // PALACE_UTILS_FILESYSTEM_HPP
