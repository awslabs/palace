// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_LABELS_HPP
#define PALACE_UTILS_LABELS_HPP

namespace palace
{

enum class CoordinateSystem : uint8_t
{
  CARTESIAN,
  CYLINDRICAL
};

enum class ProblemDataType : uint8_t
{
  DRIVEN,
  EIGENMODE,
  ELECTROSTATIC,
  MAGNETOSTATIC,
  TRANSIENT
};

enum class EigenSolverType : uint8_t
{
  DEFAULT,
  SLEPC,
  ARPACK
};

enum class SurfaceFluxType : uint8_t
{
  ELECTRIC,
  MAGNETIC,
  POWER
};

enum class InterfaceDielectricType : uint8_t
{
  DEFAULT,
  MA,
  MS,
  SA
};

enum class FrequencySampleType : u_int8_t
{
  LINEAR,
  LOG,
  POINT,
  DEFAULT = LINEAR
};

enum class TransientSolverType : uint8_t
{
  GEN_ALPHA,
  RUNGE_KUTTA,
  ARKODE,
  CVODE,
  DEFAULT = GEN_ALPHA
};

enum class ExcitationType : uint8_t
{
  SINUSOIDAL,
  GAUSSIAN,
  DIFF_GAUSSIAN,
  MOD_GAUSSIAN,
  RAMP_STEP,
  SMOOTH_STEP
};

enum class LinearSolverType : uint8_t
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

enum class Device : uint8_t
{
  CPU,
  GPU,
  DEBUG
};

}  // namespace palace

#endif  // PALACE_UTILS_LABELS_HPP
