// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_WORKSPACE_HPP
#define PALACE_UTILS_WORKSPACE_HPP

#include <mfem/general/workspace.hpp>
#include "linalg/vector.hpp"

namespace palace::workspace
{

//
// Utility class and functions for complex-valued temporary workspace vectors. See
// mfem::Workspace and mfem::WorkspaceVector.
//

class ComplexWorkspaceVector : public ComplexVector
{
private:
  mfem::WorkspaceVector data;

public:
  ComplexWorkspaceVector(int size) : data(mfem::Workspace::NewVector(2 * size))
  {
    MakeRef(data, 0, size);
  }

  // No copy constructor.
  ComplexWorkspaceVector(ComplexWorkspaceVector &other) = delete;

  // Copy assignment: copy contents of vector, not metadata.
  ComplexWorkspaceVector &operator=(ComplexWorkspaceVector &other)
  {
    ComplexVector::operator=(other);
    return *this;
  }

  // Cannot move to an existing ComplexWorkspaceVector.
  ComplexWorkspaceVector &operator=(ComplexWorkspaceVector &&other) = delete;

  // All other operator=, inherit from ComplexVector.
  using ComplexVector::operator=;
};

template <typename T>
auto NewVector(int size);

template <>
inline auto NewVector<Vector>(int size)
{
  return mfem::Workspace::NewVector(size);
}

template <>
inline auto NewVector<ComplexVector>(int size)
{
  return ComplexWorkspaceVector(size);
}

}  // namespace palace::workspace

#endif  // PALACE_UTILS_WORKSPACE_HPP
