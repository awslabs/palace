// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <mfem.hpp>

namespace palace
{

//
// Wrapper for MFEM's ParFiniteElementSpace class, where the finite element space object
// is constructed with a unique ID associated with it. This is useful for defining equality
// operations between spaces (either different spaces on the same mesh, or the same space
// type on different meshes).
//
class FiniteElementSpace : public mfem::ParFiniteElementSpace
{
private:
  mutable std::string id;
  mutable long int prev_sequence;

  // Generate a unique UUID (using system generator).
  static std::string GenerateId();

public:
  using mfem::ParFiniteElementSpace::ParFiniteElementSpace;
  FiniteElementSpace(const mfem::ParFiniteElementSpace &fespace)
    : mfem::ParFiniteElementSpace(fespace)
  {
  }

  // Get the ID associated with the instance of this class. If the underlying sequence has
  // changed (due to a mesh update, for example), regenerate the ID.
  std::string GetId() const
  {
    if (id.empty() || GetSequence() != prev_sequence)
    {
      id = GenerateId();
      prev_sequence = GetSequence();
    }
    return id;
  }
};

}  // namespace palace

#endif  // PALACE_FEM_FESPACE_HPP
