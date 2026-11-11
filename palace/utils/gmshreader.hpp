// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_GMSH_READER_HPP
#define PALACE_UTILS_GMSH_READER_HPP

#include <memory>
#include <string>

namespace mfem
{

class Mesh;

}  // namespace mfem

namespace palace::mesh
{

// Fast reader for Gmsh meshes in format 2.2 (ASCII or binary) with first-order elements and
// no periodicity, which builds the same mesh as mfem::Mesh::Loader (vertex, element, and
// boundary element order, attributes, attribute set names) up to, and including, the
// topology. The file is mapped into memory and parsed without iostreams. Returns nullptr
// for any file this reader does not handle (another format or version, compressed,
// high-order or unsupported elements, periodic meshes, unknown node tags), which is then
// read by MFEM. Without finalize_topology, the topology (faces, edges) of a mesh which
// mfem::Mesh::Finalize updates anyway is not generated (see mesh::HasTopology).
std::unique_ptr<mfem::Mesh> ReadGmshMeshFast(const std::string &filename,
                                             bool finalize_topology = true);

}  // namespace palace::mesh

#endif  // PALACE_UTILS_GMSH_READER_HPP
