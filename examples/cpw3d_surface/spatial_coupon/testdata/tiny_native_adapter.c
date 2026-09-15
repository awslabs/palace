// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test-only native stand-in for the reviewed MMG adapter command shape. It is
// linked against the fixture MMG3D library so the wrapper's link/rpath library
// resolution is exercised exactly as for the production adapter.
#include <stdio.h>

int tiny_mmg3d_fixture_version(void);

int main(int argc, char **argv)
{
  if (argc < 9 || tiny_mmg3d_fixture_version() != 1)
  {
    return 2;
  }
  FILE *existing = fopen(argv[4], "rb");
  if (existing != NULL)
  {
    fclose(existing);
    return 2;
  }
  FILE *input = fopen(argv[1], "rb");
  FILE *output = fopen(argv[4], "wb");
  if (input == NULL || output == NULL)
  {
    return 2;
  }
  char buffer[65536];
  size_t count;
  while ((count = fread(buffer, 1, sizeof buffer, input)) > 0)
  {
    if (fwrite(buffer, 1, count, output) != count)
    {
      return 2;
    }
  }
  fclose(input);
  return fclose(output) == 0 ? 0 : 2;
}
