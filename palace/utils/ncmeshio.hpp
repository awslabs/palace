// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_NCMESHIO_HPP
#define PALACE_UTILS_NCMESHIO_HPP

#include <memory>
#include <string>
#include <mfem.hpp>

namespace palace::mesh
{

struct MeshEntityCounts;

// Write a nonconforming mesh as one serial "MFEM NC mesh v1.0" file, which every process
// prints its part of, without gathering the mesh on one process, with its coordinates
// scaled by length_scale. Collective. Returns false on every process, before writing
// anything, if the mesh cannot be written this way, with the reason in `reason`. When
// `counts` is non-null on every process, it is filled on the root, as from a gathered mesh.
bool PrintParallelNCMesh(const mfem::ParMesh &mesh, double length_scale,
                         const std::string &path, bool compress,
                         MeshEntityCounts *counts = nullptr, std::string *reason = nullptr);

// Rebalance a nonconforming mesh as MFEM's default partition does (consecutive ranges of
// the leaves, in the order of the processes), but with each root of the refinement tree on
// the process its middle leaf goes to, so that saving the mesh moves no root; unless that
// is less balanced than half the tolerance on the largest load over the mean, or than the
// tolerance on the largest load over the smallest.
void RebalanceRootAligned(mfem::ParMesh &mesh, double tol);

// Scale the top-level vertex coordinates which a nonconforming mesh without nodes keeps in
// its NCMesh: MFEM's mesh transformations only move the vertices.
void ScaleNCMeshCoordinates(mfem::Mesh &mesh, double length_scale);

// The coordinates which ScaleNCMeshCoordinates scales, or null if it scales none.
mfem::Array<mfem::real_t> *GetNCMeshCoordinates(mfem::Mesh &mesh);

// What ReadNCMeshCoarse keeps of a mesh file for PartitionNCMeshFile.
struct NCMeshFile;

// Read in parallel a nonconforming mesh which PrintParallelNCMesh wrote, if `enabled`:
// every process reads a slice of an uncompressed file, and the root decompresses a
// compressed file and sends each process its slice. Returns, on every process, the text of
// its coarse mesh in `coarse` (the roots of the refinement tree, as a nonconforming mesh
// without refinements, with straight sides), and what PartitionNCMeshFile needs: the
// refinement trees of the roots which go to this process, and its slice of the nodes.
// Returns null on every process if the file cannot be read this way. Collective.
std::shared_ptr<NCMeshFile> ReadNCMeshCoarse(MPI_Comm comm, const std::string &path,
                                             bool enabled, std::string &coarse);

// Distribute that coarse mesh (on every process, which frees it once it is distributed)
// along the trees of its roots, and refine it into the mesh of the file, with the nodes of
// the file divided by the length scale; then rebalance it if the largest number of elements
// on a process is more than `tol` times the smallest. Returns null on every process if the
// file does not match (for example, if its numbering is not that of PrintParallelNCMesh).
// Collective.
std::unique_ptr<mfem::ParMesh> PartitionNCMeshFile(MPI_Comm comm, NCMeshFile &file,
                                                   std::unique_ptr<mfem::Mesh> coarse,
                                                   double length_scale, double tol);

}  // namespace palace::mesh

#endif  // PALACE_UTILS_NCMESHIO_HPP
