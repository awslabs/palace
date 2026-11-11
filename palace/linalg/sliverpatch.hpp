// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_SLIVER_PATCH_HPP
#define PALACE_LINALG_SLIVER_PATCH_HPP

#include <memory>
#include <type_traits>
#include <vector>
#include <mpi.h>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace mfem
{

template <typename T>
class Array;
class ParMesh;
class SparseMatrix;

}  // namespace mfem

namespace palace
{

class FiniteElementSpace;

// Minimum dihedral angle (in degrees) of the straight-sided tetrahedron with the given
// vertices.
double MinDihedralAngle(const double *v0, const double *v1, const double *v2,
                        const double *v3);

// Statistics of the minimum dihedral angles of the tetrahedra of a mesh (global).
struct SliverStatistics
{
  // Number of tetrahedra, of those with a minimum dihedral angle below 1, 2 and 5 degrees,
  // and below the threshold.
  long long num_tets = 0, num_below_1 = 0, num_below_2 = 0, num_below_5 = 0,
            num_slivers = 0;

  // Smallest minimum dihedral angle of the tetrahedra (degrees), 180 without tetrahedra.
  double min_angle = 180.0;
};

// The local tetrahedra of the mesh with a minimum dihedral angle below min_angle (degrees),
// using the vertex coordinates. Collective, for the global statistics.
std::vector<int> FindSliverElements(const mfem::ParMesh &mesh, double min_angle,
                                    SliverStatistics *stats = nullptr);

// Whether the MFEM device supports the sliver patch solves: with device or OpenMP backends,
// a SliverPatchSolver has no patch.
bool SliverPatchesSupported();

//
// Exact local solves on patches of degenerate (sliver) tetrahedra, for the levels of a
// geometric multigrid preconditioner: y += Σₚ Rₚᵀ Aₚₚ⁻¹ Rₚ r. The diagonal scaling of
// Chebyshev smoothers cannot correct the local modes of a sliver, whose basis functions are
// close to linearly dependent, and an exact solve on a patch around it removes them where
// they live.
//
// The patch of a sliver holds the true dofs of its own mesh entities (vertices, edges,
// faces and interior, with layers = 0) or of the entities of the process incident to its
// vertices (with layers = 1). Patches are merged while they overlap and the merged patch
// has at most max_size true dofs (otherwise the true dofs already in a patch are left out),
// so that the patches are disjoint. True dofs owned by other processes belong to patches as
// well: the process of the sliver gathers the patch block Aₚₚ of the assembled operator,
// with the contributions of the element matrices of all processes (see
// ceed::AssembleMarkedEntries), and factors it, and at each application gathers the
// residual on the patch from the owners of its true dofs and returns them their part of
// the correction. Constrained L-dofs (of nonconforming meshes) contribute to the blocks
// through the prolongation, and essential true dofs are left out.
//
template <OperatorType OperType>
class SliverPatchSolver
{
  using VecType =
      std::conditional_t<std::is_same_v<OperType, ComplexOperator>, ComplexVector, Vector>;

private:
  // Communicator of the space (duplicated, for the point-to-point messages of the
  // applications), or MPI_COMM_NULL without patches.
  MPI_Comm comm = MPI_COMM_NULL;

  // Patches owned by this process: patch p has the entries ptr[p], ..., ptr[p + 1] - 1 of
  // loc, the local true dof of each patch true dof owned by this process, or -1 - k for the
  // k-th entry of the values received from (and returned to) other processes.
  std::vector<int> ptr, loc;
  std::vector<long long> gtdofs;

  // Messages of the applications: residual values received from the processes in
  // recv_ranks (recv_ptr offsets into the received values) for the owned patches, with
  // their corrections returned in the same layout, and the residual values of the local
  // true dofs send_idx sent to the owners of patches in send_ranks (send_ptr offsets),
  // whose corrections are received in the same layout.
  std::vector<int> recv_ranks, recv_ptr, send_ranks, send_ptr, send_idx;

  // For the assembly of the patch blocks: the patch entries (global patch, position in the
  // patch, prolongation weight) of each marked L-dof (CSR, ldof_ent_ptr), the marker of
  // these L-dofs, the owner and size of each global patch, and the patch of this process
  // of each global patch (or -1).
  std::vector<int> ldof_ent_ptr, ldof_ent_patch, ldof_ent_pos;
  std::vector<double> ldof_ent_w;
  std::vector<char> marker;
  std::vector<int> patch_owner, patch_size, local_patch;

  // LU factorizations of the patch blocks (real or complex).
  struct Factors;
  std::unique_ptr<Factors> factors;
  bool active = false;

  // Residual updates (see AddMult): the L-dofs of the elements with patch L-dofs are
  // marked as well, and the owner of each patch gathers the entries of the operator
  // between the patch true dofs (columns) and all true dofs coupled to them (rows: first
  // those of the patch, in the block, then the others, in CSR form by row, cp_ptr), and
  // returns the update A[rows, patch] δ of the residual to the owners of the rows: rows
  // up_rows of the owned patches (offsets up_ptr) are local true dofs or slots -1 - k of
  // the update messages to up_send_ranks (up_send_ptr), and the updates of the local true
  // dofs up_recv_idx are received from up_recv_ranks (up_recv_ptr).
  bool update_residual = false;
  std::unique_ptr<mfem::SparseMatrix> Pd, Po;
  std::vector<long long> offd_cols, starts;
  long long my_start = 0;
  std::vector<int> up_ptr, up_rows, cp_ptr, cp_col, up_send_ranks, up_send_ptr,
      up_recv_ranks, up_recv_ptr, up_recv_idx;
  std::vector<double> cp_val[2], block_val[2];
  std::vector<std::size_t> block_off;
  std::vector<char> essential;

  // Statistics: slivers covered by the patches of this process, patch true dofs owned by
  // other processes, and L-dofs of the patch entities left out (constrained, essential).
  int num_slivers = 0;
  long long num_remote = 0, num_excluded[2] = {0, 0};

  // Work buffers for the applications.
  mutable std::vector<double> rbuf, sbuf, ubuf_s, ubuf_r, delta;
  mutable std::vector<MPI_Request> reqs;

  void Apply(const VecType &r, VecType &y, VecType *r_update) const;

public:
  // Construct the patches of the given local sliver elements of the mesh of the space.
  // Collective (the same layers and max_size on all processes).
  SliverPatchSolver(const FiniteElementSpace &fespace, const std::vector<int> &slivers,
                    int layers, const mfem::Array<int> *dbc_tdof_list, int max_size,
                    bool update_residual = false);
  ~SliverPatchSolver();

  // Number of patches owned by this process, of their true dofs, of the slivers they
  // cover, of their true dofs owned by other processes, and the largest patch size.
  int NumPatches() const { return static_cast<int>(ptr.size()) - 1; }
  int NumTrueDofs() const { return static_cast<int>(loc.size()); }
  int NumSlivers() const { return num_slivers; }
  long long NumRemoteTrueDofs() const { return num_remote; }
  int MaxPatchSize() const;

  // Sum of the squared sizes of the patches owned by this process, the dense work of an
  // application on this process (in units of the operations per patch entry).
  double PatchWork() const;

  // L-dofs of the patch entities left out of the patches of this process: constrained and
  // essential ones.
  const long long *NumExcludedLDofs() const { return num_excluded; }

  // The local true dofs of the patches owned by this process (-1 - k for true dofs owned by
  // other processes) and their global true dofs, with the patch offsets, for testing.
  const auto &GetPatchOffsets() const { return ptr; }
  const auto &GetPatchTrueDofs() const { return loc; }
  const auto &GetPatchGlobalTrueDofs() const { return gtdofs; }

  // Assemble and factor the patch blocks of the operator, a ParOperator or
  // ComplexParOperator on the space of the constructor (with its essential true dofs).
  // Collective. Returns false and leaves the solver inactive on all processes if a local
  // operator is not supported, and true if all patch blocks are numerically nonsingular
  // (patches with singular blocks are skipped).
  bool SetOperator(const OperType &op);

  // Whether this process takes part in the applications (owns patches or true dofs of
  // patches of other processes) after a successful SetOperator.
  bool IsActive() const { return active; }

  // y += Σₚ Rₚᵀ Aₚₚ⁻¹ Rₚ r. Collective over the processes which take part (IsActive).
  void AddMult(const VecType &r, VecType &y) const { Apply(r, y, nullptr); }

  // y += δ = Σₚ Rₚᵀ Aₚₚ⁻¹ Rₚ r and r -= A δ, the residual of the corrected y (zero on the
  // patches, up to the couplings between them), for a solver constructed with
  // update_residual. Collective over the processes which take part (IsActive).
  void AddMultUpdate(VecType &r, VecType &y) const;
};

}  // namespace palace

#endif  // PALACE_LINALG_SLIVER_PATCH_HPP
