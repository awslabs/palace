// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test-only shared library published as libmmg3d so the fixture adapter has a
// runtime-resolved MMG3D dependency whose path and hash must be bound.
int tiny_mmg3d_fixture_version(void)
{
  return 1;
}
