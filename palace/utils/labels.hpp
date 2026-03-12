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

// The types of problem Palace is able to solve.
enum class ProblemType : char
{
  DRIVEN,
  EIGENMODE,
  ELECTROSTATIC,
  MAGNETOSTATIC,
  TRANSIENT,
  BOUNDARYMODE
};

// Eigenvalue solver type.
enum class EigenSolverBackend : char
{
  DEFAULT,
  SLEPC,
  ARPACK
};

// Nonlinear eigenvalue solver type.
enum class NonlinearEigenSolver : char
{
  HYBRID,
  SLP
};

// Surface fluxes.
enum class SurfaceFlux : char
{
  ELECTRIC,
  MAGNETIC,
  POWER
};

// Interface dielectrics for computing electric field energy participation ratios.
enum class InterfaceDielectric : char
{
  DEFAULT,
  MA,
  MS,
  SA
};

// Frequency sampling schemes.
enum class FrequencySampling : char
{
  LINEAR,
  LOG,
  POINT,
  DEFAULT = LINEAR
};

// Time integration scheme type.
enum class TimeSteppingScheme : char
{
  GEN_ALPHA,
  RUNGE_KUTTA,
  ARKODE,
  CVODE,
  DEFAULT = GEN_ALPHA
};

// Excitation type for port excitation.
enum class Excitation : char
{
  SINUSOIDAL,
  GAUSSIAN,
  DIFF_GAUSSIAN,
  MOD_GAUSSIAN,
  RAMP_STEP,
  SMOOTH_STEP
};

// Possible linear solvers
enum class LinearSolver : char
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

// Krylov solvers to use in the linear solver.
enum class KrylovSolver : char
{
  DEFAULT,
  CG,
  MINRES,
  GMRES,
  FGMRES,
  BICGSTAB
};

// Method of coarsening for p-multigrid.
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

// Domain orthogonalization type for adaptive circuit synthesis.
enum class DomainOrthogonalizationWeight : char
{
  ENERGY,
  FE_BASIS_IDENTITY,
  SPACE_OVERLAP
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
