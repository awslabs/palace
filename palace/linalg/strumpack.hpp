// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_HPP
#define PALACE_LINALG_STRUMPACK_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include "linalg/operator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

// Set the compression options of STRUMPACK, for the real-valued solvers and for the solver
// in complex arithmetic (StrumpackComplexSolver).
template <typename ScalarType>
void ConfigureStrumpackCompression(strumpack::SPOptions<ScalarType> &opts,
                                   SparseCompression compression, double lr_tol,
                                   int butterfly_l, int lossy_prec);

// The column ordering of STRUMPACK for "ColumnOrdering": "Default", in the real-valued
// solver and in the solver in complex arithmetic (the coarse solvers and the solver of the
// PML regions): a METIS nested dissection which is parallel over the separator tree
// (ReorderingStrategy::METIS_TREE). With METIS, the default of STRUMPACK, one process
// orders the graph of the whole matrix while all other processes wait. With this ordering,
// that process computes only the top separators, and the processes then order different
// parts of the graph at the same time, for factors of about the same size. On one process,
// or for a matrix of fewer than 4000 rows, it is the METIS ordering. The ordering is a
// patch in the STRUMPACK builds of Palace (superbuild and Spack package), and STRUMPACK
// defines no macro for it: with a STRUMPACK which lacks it, the default is METIS.
template <typename Strategy = strumpack::ReorderingStrategy>
constexpr Strategy GetStrumpackDefaultReordering()
{
  if constexpr (requires { Strategy::METIS_TREE; })
  {
    return Strategy::METIS_TREE;
  }
  else
  {
    return Strategy::METIS;
  }
}

// Set the column ordering of STRUMPACK, for the real-valued solver and for the solver in
// complex arithmetic (StrumpackComplexSolver). The mixed precision solver always uses
// METIS.
template <typename ScalarType>
void ConfigureStrumpackReordering(strumpack::SPOptions<ScalarType> &opts,
                                  SymbolicFactorization reorder);

// Enable the STRUMPACK option which stores the fill-reducing permutation of the distributed
// solver and its inverse (two integers per matrix row) once per node, in an MPI shared
// memory window, instead of on every process. Palace enables it in all its STRUMPACK
// solvers, in real and in complex arithmetic. STRUMPACK does this after a sequential
// ordering or the default ordering above (not ParMETIS or PT-Scotch) without compression,
// and falls back to a copy per process, with a warning, if the window cannot be allocated
// (with Open MPI, OMPI_MCA_osc=^sm,ucx makes every process keep its copy). The option is a
// patch in the STRUMPACK builds of Palace (superbuild and Spack package), and STRUMPACK
// defines no macro for it: this does nothing with a STRUMPACK which lacks it.
template <typename ScalarType>
inline void EnableStrumpackSharedPermutation(strumpack::SPOptions<ScalarType> &opts)
{
  if constexpr (requires { opts.enable_shared_permutation(); })
  {
    opts.enable_shared_permutation();
  }
}

// Make STRUMPACK map the subtrees of its elimination tree to the processes in proportion to
// the size of their factors instead of the operations of their factorization (its default).
// A triangular solve costs the factor size, not the factorization operations: the default
// can leave one process a subtree with many times the average factor size, and every solve
// then waits for that process. Palace sets this in its coarse solvers (StrumpackSolver and
// StrumpackComplexSolver), where a factorization serves many solves.
template <typename ScalarType>
inline void SetStrumpackFactorSizeMapping(strumpack::SPOptions<ScalarType> &opts)
{
  opts.set_proportional_mapping(strumpack::ProportionalMapping::FACTOR_MEMORY);
}

// Enable the STRUMPACK option which replaces the triangular factors of the fronts
// distributed over several processes by their inverses after each factorization (selective
// inversion). The solves with these fronts then use matrix products, which need a constant
// number of messages, instead of distributed triangular solves: each solve gets faster on
// many processes, and each factorization slower. STRUMPACK inverts the factors T of a front
// only if their condition number || |T⁻¹| |T| ||_∞ is at most a limit, as a solve with the
// inverted factors can lose backward accuracy up to about that factor, and keeps the
// triangular solves in the other fronts. Palace sets the limit to 1e8 instead of the
// STRUMPACK default of 1e3: on the eigenmode problems this was measured on, 1e8 inverts
// every front where 1e3 keeps the triangular solves in the largest ones, where most of the
// solve time is spent, and on a nearly singular low-frequency driven problem the results
// were accurate to 2e-7 with limits of 1e8 and 1e10 and wrong with 1e12. The compressed
// fronts are not changed, and STRUMPACK ignores the option when built with SLATE. As for
// the shared permutation, this does nothing with a STRUMPACK which lacks the option.
template <typename ScalarType>
inline void EnableStrumpackSelectiveInversion(strumpack::SPOptions<ScalarType> &opts)
{
  constexpr double max_cond = 1.0e8;
  if constexpr (requires {
                  opts.enable_selective_inversion();
                  opts.set_selective_inversion_max_cond(max_cond);
                })
  {
    opts.enable_selective_inversion();
    opts.set_selective_inversion_max_cond(max_cond);
  }
}

//
// A wrapper for the STRUMPACK direct solver package.
//
template <typename StrumpackSolverType>
class StrumpackSolverBase : public StrumpackSolverType
{
private:
  MPI_Comm comm;

public:
  StrumpackSolverBase(MPI_Comm comm, SymbolicFactorization reorder,
                      SparseCompression compression, double lr_tol, int butterfly_l,
                      int lossy_prec, bool reorder_reuse, bool selective_inversion,
                      int print, bool factor_size_mapping = true);

  StrumpackSolverBase(const IoData &iodata, MPI_Comm comm, int print)
    : StrumpackSolverBase(comm, iodata.solver.linear.sym_factorization,
                          iodata.solver.linear.strumpack_compression_type,
                          iodata.solver.linear.strumpack_lr_tol,
                          iodata.solver.linear.strumpack_butterfly_l,
                          iodata.solver.linear.strumpack_lossy_precision,
                          iodata.solver.linear.reorder_reuse,
                          iodata.solver.linear.strumpack_selective_inversion > 0, print)
  {
  }

  void SetOperator(const Operator &op) override;

  void SetReorderReuse(bool reorder_reuse);

  // The ordering which STRUMPACK uses: for the mixed precision solver, that of its single
  // precision solver, which computes the ordering (METIS for every value of the option).
  strumpack::ReorderingStrategy GetReorderingStrategy() const;
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
