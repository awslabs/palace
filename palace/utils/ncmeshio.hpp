// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_NCMESHIO_HPP
#define PALACE_UTILS_NCMESHIO_HPP

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

// Rebalance a nonconforming mesh as the partition `base` of its elements does, but with
// each root of the refinement tree on one process, so that saving the mesh moves no root:
// the process which gets the most of its leaves in `base`. A root already split between
// processes keeps the partition of `base`. Without `base`, as MFEM's default partition
// does (consecutive ranges of the leaves, in the order of the processes), with each root
// on the process its middle leaf goes to. Rebalance as `base` does instead if that is less
// balanced than half the tolerance on the largest load over the mean, or than the
// tolerance on the largest load over the smallest.
void RebalanceRootAligned(mfem::ParMesh &mesh, double tol,
                          const mfem::Array<int> *base = nullptr);

// Scale the top-level vertex coordinates which a nonconforming mesh without nodes keeps in
// its NCMesh: MFEM's mesh transformations only move the vertices.
void ScaleNCMeshCoordinates(mfem::Mesh &mesh, double length_scale);

}  // namespace palace::mesh

#endif  // PALACE_UTILS_NCMESHIO_HPP
