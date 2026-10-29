// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_TEST_HELPERS_HPP
#define PALACE_TEST_HELPERS_HPP

#include <string>
#include <mfem.hpp>

mfem::Mesh SingleTetMesh();

// Path to the mesh with the given file name in MFEM's data directory. "star-tri.mesh" and
// "fichera-tet.mesh" are not in MFEM: they are "star.mesh" and "fichera.mesh" after
// mfem::Mesh::MakeSimplicial, written on first use to a temporary file private to the
// calling process and removed at exit.
std::string TestMeshPath(const std::string &name);

#endif  // PALACE_TEST_HELPERS_HPP
