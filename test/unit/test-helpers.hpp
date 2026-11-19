// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <string>
#include <mfem.hpp>

mfem::Mesh SingleTetMesh();

// A box [0, 2] x [0, 1]^2 of 4 x 2 x 2 hexahedra: the domains x < 1 (attribute 1) and x > 1
// (2), the faces x = 0 (1), x = 2 (2) and the others (3), and the interior boundaries x = 1
// between the domains (4) and x = 3/2 inside the second (5).
mfem::Mesh MakeInteriorBoundaryBox();

// Path to the unit test mesh with the given file name. Most test meshes are stored in
// PALACE_TEST_DATA_DIR/mesh. "star-tri.mesh" and "fichera-tet.mesh" are not stored: they
// are "star-quad.mesh" and "fichera-hex.mesh" after mfem::Mesh::MakeSimplicial, written on
// first use to a temporary file private to the calling process and removed at exit.
std::string TestMeshPath(const std::string &name);

#endif  // PALACE_TEST_HELPERS_HPP
