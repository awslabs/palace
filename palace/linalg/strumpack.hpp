// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_STRUMPACK_HPP
#define PALACE_LINALG_STRUMPACK_HPP

#include <mfem.hpp>

#if defined(MFEM_USE_STRUMPACK)

#include "linalg/operator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

// Set the compression options of STRUMPACK, for the real-valued solvers and for the solver
// in complex arithmetic (StrumpackComplexSolver).
template <typename ScalarType>
void ConfigureStrumpackCompression(strumpack::SPOptions<ScalarType> &opts,
                                   SparseCompression compression, double lr_tol,
                                   int butterfly_l, int lossy_prec);

// Enable the STRUMPACK option which stores the fill-reducing permutation of the distributed
// solver and its inverse (two integers per matrix row) once per node, in an MPI shared
// memory window, instead of on every process. Palace enables it in all its STRUMPACK
// solvers, in real and in complex arithmetic. STRUMPACK does this after a sequential
// ordering (not ParMETIS or PT-Scotch) without compression, and falls back to a copy per
// process, with a warning, if the window cannot be allocated (with Open MPI,
// OMPI_MCA_osc=^sm,ucx makes every process keep its copy). The option is a patch in the
// STRUMPACK builds of Palace (superbuild and Spack package), and STRUMPACK defines no macro
// for it: this does nothing with a STRUMPACK which lacks it.
template <typename ScalarType>
inline void EnableStrumpackSharedPermutation(strumpack::SPOptions<ScalarType> &opts)
{
  if constexpr (requires { opts.enable_shared_permutation(); })
  {
    opts.enable_shared_permutation();
  }
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

// LAB (si3), not for a pull request. A value of "STRUMPACKButterflyLevels" of 1000 or more
// carries lab flags (value - 1000), and the number of butterfly levels is then 1. Flag 16:
// STRUMPACK maps the elimination tree to the processes by factor size instead of by the
// operations of the factorization. The flags act on the solver objects of the coarse solver
// which can get the selective inversion (real-valued and complex arithmetic), not on the
// mixed precision solver nor on the mode solver.
inline int StrumpackLabFlags(int butterfly_l)
{
  return (butterfly_l >= 1000) ? butterfly_l - 1000 : 0;
}

inline int StrumpackLabButterflyLevels(int butterfly_l)
{
  return (butterfly_l >= 1000) ? 1 : butterfly_l;
}

template <typename ScalarType>
inline void ApplyStrumpackLabFlags(strumpack::SPOptions<ScalarType> &opts, int butterfly_l,
                                   MPI_Comm comm)
{
  if (StrumpackLabFlags(butterfly_l) & 16)
  {
    opts.set_proportional_mapping(strumpack::ProportionalMapping::FACTOR_MEMORY);
    Mpi::Print(comm, "LAB si3: proportional mapping by factor size\n");
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
                      int print);

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
};

using StrumpackSolver = StrumpackSolverBase<mfem::STRUMPACKSolver>;

using StrumpackMixedPrecisionSolver =
    StrumpackSolverBase<mfem::STRUMPACKMixedPrecisionSolver>;

}  // namespace palace

#endif

#endif  // PALACE_LINALG_STRUMPACK_HPP
