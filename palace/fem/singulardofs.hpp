// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARDOFS_HPP
#define PALACE_FEM_SINGULARDOFS_HPP

#include <array>
#include <cstddef>
#include <vector>
#include <mfem.hpp>

#include "fem/singularelements.hpp"
#include "fem/singularfeatures.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

// Canonical mesh entity represented only by sorted stable vertex IDs. Unused
// entries are -1. This avoids dependence on MFEM's traversal-dependent local
// edge, face, and element numbers.
struct EntityKey
{
  std::size_t size = 0;
  std::array<GlobalVertexId, 4> vertices{-1, -1, -1, -1};

  bool operator==(const EntityKey &other) const;
  bool operator<(const EntityKey &other) const;
};

// Global identity of one retained singular basis function. The singular entity
// is a mesh vertex for node families and a mesh edge for edge families. Support
// is the interpolation entity on which the function has a conforming trace:
// edge, face, or element. Component distinguishes retained vector directions
// at interpolation points with more than one basis function.
struct DofKey
{
  HigherOrderBasisFamily family;
  int order;
  EntityKey singular_entity;
  EntityKey support_entity;
  EntityKey component_entity;
  std::array<int, 4> interpolation_weights{0, 0, 0, 0};

  bool operator==(const DofKey &other) const;
  bool operator<(const DofKey &other) const;
};

struct ElementDof
{
  // Rank-local index into DofTopology::{h1,nd}_dofs. A parallel true-DOF
  // number is deliberately not assigned at this stage.
  std::size_t dof;
  HigherOrderBasis basis;
};

struct ElementDofMap
{
  // Scalar singular-potential DOFs. Every entry is a gradient family.
  std::vector<ElementDof> h1;

  // Vector singular-field DOFs. This contains the same gradient entries as h1
  // plus the rotational families.
  std::vector<ElementDof> nd;
};

struct DofTopology
{
  // Sorted unique keys present on this mesh or rank.
  std::vector<DofKey> h1_dofs;
  std::vector<DofKey> nd_dofs;

  // Exact enriched discrete gradient: h1_to_nd[i] is the ND basis which is the
  // analytic gradient of scalar H1 singular basis i. Its coefficient is +1.
  std::vector<std::size_t> h1_to_nd;
  std::vector<ElementDofMap> elements;
};

// Lowest-order triangular basis descriptor. The gradient basis with nodes
// {i,j,k} is grad(phi_ij); the rotational basis is N_i,r with the same
// canonical ordering. Higher triangular singular orders require separate
// Graglia-Lombardi interpolation tuples and are not represented by this type.
struct TriangleBasis
{
  HigherOrderBasisFamily family;
  std::array<int, 3> nodes;
  int order;
  double nu;
};

struct TriangleElementDof
{
  std::size_t dof;
  TriangleBasis basis;
};

struct TriangleElementDofMap
{
  std::vector<TriangleElementDof> h1;
  std::vector<TriangleElementDof> nd;
};

struct TriangleDofTopology
{
  std::vector<DofKey> h1_dofs;
  std::vector<DofKey> nd_dofs;
  std::vector<std::size_t> h1_to_nd;
  std::vector<TriangleElementDofMap> elements;
};

struct TrueDofMap
{
  // Partition of the rank-local enrichment L-vectors. Shared functions occur
  // once on every participating rank in this numbering.
  HYPRE_BigInt global_local_size = 0;
  HYPRE_BigInt local_offset = 0;
  HYPRE_BigInt local_size = 0;

  // Partition of the uniquely owned parallel true DOFs.
  HYPRE_BigInt global_size = 0;
  HYPRE_BigInt owned_offset = 0;
  HYPRE_BigInt owned_size = 0;

  // Ownership and global true-DOF number for every rank-local canonical key.
  std::vector<int> owner;
  std::vector<HYPRE_BigInt> local_to_true;
};

struct ParallelDofNumbering
{
  TrueDofMap h1;
  TrueDofMap nd;

  // Global ND true DOF for the analytic gradient of each rank-local H1
  // singular basis. The coefficient remains exactly +1.
  std::vector<HYPRE_BigInt> h1_to_nd_true;
};

bool IsGradientFamily(HigherOrderBasisFamily family);

// Map every rank-local H1 enrichment key to the straight singular features
// which own its singular entity. Edge-gradient keys belong to exactly one
// straight feature. Node-gradient keys belong to every feature incident on
// the singular vertex, producing the intentional overlap at polygon corners.
std::vector<std::vector<std::size_t>>
BuildH1DofFeatureMembership(const FeatureTopology &features, const DofTopology &topology);

// Map every triangular H1 enrichment key to its singular line-tip patch. Each
// physical tip is one feature because the initial 2D extractor supports only
// straight, unbranched selected PEC line chains.
std::vector<std::vector<std::size_t>>
BuildTriangleH1DofFeatureMembership(const TriangleFeatureTopology &features,
                                    const TriangleDofTopology &topology);

// Recover decomposition-independent serial vertex IDs on a ParMesh constructed
// from serial_mesh with the supplied element partition. The map is exact and
// topological: it uses the preserved serial-to-local element ordering and
// vertex incidence, with coordinates used only as a consistency check.
//
// The result is invalidated by subsequent refinement or rebalancing.
std::vector<GlobalVertexId>
MapPartitionedSerialVertexIds(const mfem::Mesh &serial_mesh,
                              const mfem::ParMesh &parallel_mesh, const int *partitioning);

// Enumerate the retained equations (28), (30), (32), and (33) on every
// enriched tetrahedron and assign deterministic serial global DOFs. The
// feature topology and mesh must come from the same serial mesh. Order one is
// exactly the first-order equations (6), (11), (21), and (23).
DofTopology BuildSerialDofTopology(const mfem::Mesh &mesh, const FeatureTopology &features,
                                   int order);

// Enumerate rank-local retained bases using decomposition-independent serial
// vertex IDs. This constructs local key tables and element maps only; MPI
// ownership and true-DOF numbering are a subsequent phase.
DofTopology BuildLocalDofTopology(
    const mfem::Mesh &mesh, const FeatureTopology &features,
    const std::vector<GlobalVertexId> &serial_vertex_ids, int order,
    const std::vector<std::array<GlobalVertexId, 3>> &excluded_boundary_trace_faces = {});

// Enumerate the two edge-trace gradient functions and one element-interior
// rotational function associated with each singular tip/triangle incidence.
// The initial triangular implementation supports only singular order one.
TriangleDofTopology BuildSerialTriangleDofTopology(const mfem::Mesh &mesh,
                                                   const TriangleFeatureTopology &features,
                                                   int order);

TriangleDofTopology BuildLocalTriangleDofTopology(
    const mfem::Mesh &mesh, const TriangleFeatureTopology &features,
    const std::vector<GlobalVertexId> &serial_vertex_ids, int order);

// Assign MPI ownership and true-DOF numbers to the canonical keys present on
// each rank. The minimum participating rank owns a shared key, and owned keys
// occupy one contiguous block per rank.
//
// Fixed-width integer keys use a deterministic distributed rendezvous. Each
// owner sorts its keys before assigning true-DOF numbers, preserving
// decomposition-independent numbering without replicating the global key
// table.
ParallelDofNumbering BuildParallelDofNumbering(MPI_Comm comm, const DofTopology &topology);

ParallelDofNumbering BuildParallelDofNumbering(MPI_Comm comm,
                                               const TriangleDofTopology &topology);

// Return rank-local owned H1 true DOFs whose scalar enrichment trace is
// nonzero on a selected PEC sheet face. The classification is topological: a
// retained interpolation support entity lies on a sheet face if all of its
// stable source vertices belong to that face.
mfem::Array<int> GetEssentialH1TrueDofs(MPI_Comm comm, const FeatureTopology &features,
                                        const DofTopology &topology,
                                        const ParallelDofNumbering &parallel_numbering);
mfem::Array<int>
GetEssentialH1TrueDofs(MPI_Comm comm, const FeatureTopology &features,
                       const DofTopology &topology,
                       const ParallelDofNumbering &parallel_numbering,
                       const mfem::Array<int> &essential_boundary_attributes);

// Return rank-local owned H1 true DOFs whose scalar enrichment trace is
// nonzero on any supplied boundary face. Faces use sorted stable source-vertex
// IDs and need not belong to the attributes selected for singular enrichment.
mfem::Array<int> GetEssentialH1TrueDofsOnFaces(
    MPI_Comm comm, const std::vector<std::array<GlobalVertexId, 3>> &boundary_faces,
    const DofTopology &topology, const ParallelDofNumbering &parallel_numbering);

// Return rank-local owned ND true DOFs whose vector enrichment has nonzero
// tangential trace on a selected PEC sheet face. As for H1, trace support is
// classified from stable source vertices, so the result is independent of
// element ordering and mesh partitioning.
mfem::Array<int> GetEssentialNDTrueDofs(MPI_Comm comm, const FeatureTopology &features,
                                        const DofTopology &topology,
                                        const ParallelDofNumbering &parallel_numbering);
mfem::Array<int>
GetEssentialNDTrueDofs(MPI_Comm comm, const FeatureTopology &features,
                       const DofTopology &topology,
                       const ParallelDofNumbering &parallel_numbering,
                       const mfem::Array<int> &essential_boundary_attributes);

// Return owned enriched H1 true DOFs whose nonzero trace lies on a selected
// zero-thickness PEC line segment.
mfem::Array<int>
GetEssentialTriangleH1TrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering);
mfem::Array<int>
GetEssentialTriangleH1TrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering,
                               const mfem::Array<int> &essential_boundary_attributes);

// Return rank-local owned H1 true DOFs whose scalar enrichment trace is
// nonzero on any supplied boundary segment. Segments use sorted stable
// source-vertex IDs and need not be selected singular line segments.
mfem::Array<int> GetEssentialTriangleH1TrueDofsOnSegments(
    MPI_Comm comm, const std::vector<std::array<GlobalVertexId, 2>> &boundary_segments,
    const TriangleDofTopology &topology, const ParallelDofNumbering &parallel_numbering);

// Return owned enriched ND true DOFs with nonzero tangential trace on a
// selected zero-thickness PEC line segment. Rotational triangle enrichments
// have element-interior support and are therefore never included.
mfem::Array<int>
GetEssentialTriangleNDTrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering);
mfem::Array<int>
GetEssentialTriangleNDTrueDofs(MPI_Comm comm, const TriangleFeatureTopology &features,
                               const TriangleDofTopology &topology,
                               const ParallelDofNumbering &parallel_numbering,
                               const mfem::Array<int> &essential_boundary_attributes);

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARDOFS_HPP
