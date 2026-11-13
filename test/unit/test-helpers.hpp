// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <string>
#include <mfem.hpp>
#include "linalg/complexcoarse.hpp"

mfem::Mesh SingleTetMesh();

// Path to the unit test mesh with the given file name. Most test meshes are stored in
// PALACE_TEST_DATA_DIR/mesh. "star-tri.mesh" and "fichera-tet.mesh" are not stored: they
// are "star-quad.mesh" and "fichera-hex.mesh" after mfem::Mesh::MakeSimplicial, written on
// first use to a temporary file private to the calling process and removed at exit.
std::string TestMeshPath(const std::string &name);

namespace palace
{

// Forget the measured choice of the complex-valued coarse solve (see
// MeasuredComplexCoarseSolver) before and after a test which builds a complex-valued sparse
// direct coarse solver with ComplexCoarseSolve omitted, so that no choice leaks between
// tests.
struct ComplexCoarseChoiceReset
{
  ComplexCoarseChoiceReset() { MeasuredComplexCoarseSolver::ResetChoice(); }
  ~ComplexCoarseChoiceReset() { MeasuredComplexCoarseSolver::ResetChoice(); }
};

}  // namespace palace

#endif  // PALACE_TEST_HELPERS_HPP
