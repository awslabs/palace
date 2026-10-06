// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// gslib.h leaks C99 compatibility macros (e.g. an empty "inline") into C++, so it is kept
// out of every other translation unit.
#include <gslib.h>

// Not defined before GSLIB v1.0.8 (same fallback as MFEM)
#if !defined(GSLIB_RELEASE_VERSION)
#define GSLIB_RELEASE_VERSION 10007
#endif

namespace palace
{

int GetGslibReleaseVersion()
{
  return GSLIB_RELEASE_VERSION;
}

}  // namespace palace
