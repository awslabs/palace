// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP
#define PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP

#include <complex>
#include <memory>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class FiniteElementSpace;
class MaterialOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;
class SurfaceRationalImpedanceOperator;

namespace mode_assembly
{

using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                         std::unique_ptr<mfem::HypreParMatrix>>;

enum class CoefficientType
{
  CONSTANT,
  OMEGA,
  OMEGA_SQUARED,
  SHIFT,
  SURFACE_CONDUCTIVITY,
  RATIONAL_IMPEDANCE
};

// One complete, BC-eliminated component of the mode eigenproblem matrix A. Exact assembly
// combines these sparse matrices, while the wave-port ROM projects the same operators.
struct OperatorComponent
{
  CoefficientType type;
  int index;
  std::unique_ptr<mfem::HypreParMatrix> Ar, Ai;
  std::unique_ptr<ComplexWrapperOperator> op;

  OperatorComponent(CoefficientType type, int index,
                    std::unique_ptr<mfem::HypreParMatrix> &&Ar,
                    std::unique_ptr<mfem::HypreParMatrix> &&Ai);
  OperatorComponent(OperatorComponent &&) = default;
  OperatorComponent &operator=(OperatorComponent &&) = default;
  OperatorComponent(const OperatorComponent &) = delete;
  OperatorComponent &operator=(const OperatorComponent &) = delete;
};

// Frequency-parametric mode eigenproblem matrix. The complete fixed components are the
// single source of truth for both exact sparse assembly and wave-port projection.
class ModeOperatorModel
{
private:
  SurfaceConductivityOperator &surf_sigma_op;
  SurfaceRationalImpedanceOperator &surf_rz_op;
  std::vector<OperatorComponent> components;

public:
  ModeOperatorModel(const FiniteElementSpace &nd_fespace,
                    const FiniteElementSpace &h1_fespace, const MaterialOperator &mat_op,
                    const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
                    FarfieldBoundaryOperator &farfield_op,
                    SurfaceConductivityOperator &surf_sigma_op,
                    SurfaceRationalImpedanceOperator &surf_rz_op,
                    const mfem::HypreParMatrix &Bttr, const mfem::HypreParMatrix *Atnr,
                    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Btnr,
                    const mfem::Array<int> &dbc_tdof_list);

  const std::vector<OperatorComponent> &GetComponents() const { return components; }
  std::complex<double> EvaluateCoefficient(const OperatorComponent &component,
                                           std::complex<double> omega, double sigma) const;
  std::unique_ptr<ComplexOperator> Assemble(std::complex<double> omega, double sigma) const;
};

ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op);
ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op);

ComplexHypreParMatrix BuildSystemMatrixA(
    const FiniteElementSpace &nd_fespace, const FiniteElementSpace &h1_fespace,
    const mfem::Array<int> &dbc_tdof_list, const mfem::HypreParMatrix *Attr,
    const mfem::HypreParMatrix *Atti, const mfem::HypreParMatrix *Atnr,
    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Annr,
    const mfem::HypreParMatrix *Anni, const mfem::HypreParMatrix *shifted_Btnr = nullptr,
    Operator::DiagonalPolicy diag_policy = Operator::DIAG_ONE);

void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en);

}  // namespace mode_assembly
}  // namespace palace

#endif  // PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP
