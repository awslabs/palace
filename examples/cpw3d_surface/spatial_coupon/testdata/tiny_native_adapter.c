// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Test-only native stand-in for the reviewed MMG adapter command shape. It is
// linked against the fixture MMG3D library so the wrapper's link/rpath library
// resolution is exercised exactly as for the production adapter. Like the
// adapter it takes the metric stage's required-tetrahedra list as the flag
// --required-tetrahedra FILE (removed from the positional arguments).
#include <stdio.h>
#include <string.h>

int tiny_mmg3d_fixture_version(void);

int main(int argc, char **argv)
{
  char *positional[16];
  int count = 0, required = 0;
  for (int i = 0; i < argc; i++)
  {
    if (i > 0 && strcmp(argv[i], "--required-tetrahedra") == 0)
    {
      if (i + 1 >= argc || required)
      {
        return 2;
      }
      FILE *list = fopen(argv[++i], "rb");
      if (list == NULL)
      {
        return 2;
      }
      fclose(list);
      required = 1;
      continue;
    }
    if (count == 16)
    {
      return 2;
    }
    positional[count++] = argv[i];
  }
  if (count < 9 || !required || tiny_mmg3d_fixture_version() != 1)
  {
    return 2;
  }
  FILE *existing = fopen(positional[4], "rb");
  if (existing != NULL)
  {
    fclose(existing);
    return 2;
  }
  FILE *input = fopen(positional[1], "rb");
  FILE *output = fopen(positional[4], "wb");
  if (input == NULL || output == NULL)
  {
    return 2;
  }
  char buffer[65536];
  size_t bytes;
  while ((bytes = fread(buffer, 1, sizeof buffer, input)) > 0)
  {
    if (fwrite(buffer, 1, bytes, output) != bytes)
    {
      return 2;
    }
  }
  fclose(input);
  return fclose(output) == 0 ? 0 : 2;
}
