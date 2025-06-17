// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_LABELS_HPP
#define PALACE_UTILS_LABELS_HPP

namespace palace
{

// Usable coordinate systems.
enum class CoordinateSystem : char
{
  CARTESIAN,
  CYLINDRICAL
};

// Problem types Palace is able to solve.
enum class ProblemType : char
{
  DRIVEN,
  EIGENMODE,
  ELECTROSTATIC,
  MAGNETOSTATIC,
  TRANSIENT
};

// Eigenvalue solver type.
enum class EigenSolverType : char
{
  DEFAULT,
  SLEPC,
  ARPACK
};

// Surface fluxes.
enum class SurfaceFluxType : char
{
  ELECTRIC,
  MAGNETIC,
  POWER
};

// Interface dielectrics for computing electric field energy participation ratios.
enum class InterfaceDielectricType : char
{
  DEFAULT,
  MA,
  MS,
  SA
};

// Frequency sampling schemes.
enum class FrequencySampleType : char
{
  LINEAR,
  LOG,
  POINT,
  DEFAULT = LINEAR
};

// Time integration scheme type.
enum class TransientSolverType : char
{
  GEN_ALPHA,
  RUNGE_KUTTA,
  ARKODE,
  CVODE,
  DEFAULT = GEN_ALPHA
};

// Excitation type for port excitation.
enum class ExcitationType : char
{
  SINUSOIDAL,
  GAUSSIAN,
  DIFF_GAUSSIAN,
  MOD_GAUSSIAN,
  RAMP_STEP,
  SMOOTH_STEP
};

// Type of linear solver.
enum class LinearSolverType : char
{
  DEFAULT,
  AMS,
  BOOMER_AMG,
  MUMPS,
  SUPERLU,
  STRUMPACK,
  STRUMPACK_MP,
  JACOBI
};

// Krylov solver type in linear solver.
enum class KrylovSolver : char
{
  DEFAULT,
  CG,
  MINRES,
  GMRES,
  FGMRES,
  BICGSTAB
};

// Type of coarsening for p-multigrid.
enum class MultigridCoarsening : char
{
  LINEAR,
  LOGARITHMIC
};

// Preconditioning side.
enum class PreconditionerSide : char
{
  DEFAULT,
  RIGHT,
  LEFT
};

// Column ordering method in the symbolic factorization for sparse direct solvers.
enum class SymbolicFactorization : char
{
  DEFAULT,
  METIS,
  PARMETIS,
  SCOTCH,
  PTSCOTCH,
  PORD,
  AMD,
  RCM
};

// Low-rank and butterfly compression scheme for sparse direct solvers which support it
// (mainly STRUMPACK).
enum class SparseCompression : char
{
  NONE,
  BLR,
  HSS,
  HODLR,
  ZFP,
  BLR_HODLR,
  ZFP_BLR_HODLR
};

// Variations of Gram-Schmidt orthogonalization for GMRES/FGMRES iterative solvers and SLEPc
// eigenvalue solver.
enum class Orthogonalization : char
{
  MGS,
  CGS,
  CGS2
};

// Device used to configure MFEM.
enum class Device : char
{
  CPU,
  GPU,
  DEBUG
};

}  // namespace palace

#endif  // PALACE_UTILS_LABELS_HPP
