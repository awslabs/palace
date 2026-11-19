// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <cstdlib>
#include <optional>
#include <string>
#include <mfem.hpp>

mfem::Mesh SingleTetMesh();

// Path to the unit test mesh with the given file name. Most test meshes are stored in
// PALACE_TEST_DATA_DIR/mesh. "star-tri.mesh" and "fichera-tet.mesh" are not stored: they
// are "star-quad.mesh" and "fichera-hex.mesh" after mfem::Mesh::MakeSimplicial, written on
// first use to a temporary file private to the calling process and removed at exit.
std::string TestMeshPath(const std::string &name);

// Sets an environment variable of the calling process, or unsets it when the value is
// null, and restores its previous state on destruction.
class ScopedEnvironmentVariable
{
private:
  std::string name;
  std::optional<std::string> previous;

  void Set(const char *value)
  {
    if (value)
    {
      setenv(name.c_str(), value, 1);
    }
    else
    {
      unsetenv(name.c_str());
    }
  }

public:
  ScopedEnvironmentVariable(const std::string &name, const char *value) : name(name)
  {
    if (const char *old = std::getenv(name.c_str()))
    {
      previous = old;
    }
    Set(value);
  }
  ScopedEnvironmentVariable(const ScopedEnvironmentVariable &) = delete;
  ScopedEnvironmentVariable &operator=(const ScopedEnvironmentVariable &) = delete;
  ~ScopedEnvironmentVariable() { Set(previous ? previous->c_str() : nullptr); }
};

#endif  // PALACE_TEST_HELPERS_HPP
