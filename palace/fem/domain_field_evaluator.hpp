// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_DOMAIN_FIELD_EVALUATOR_HPP
#define PALACE_FEM_DOMAIN_FIELD_EVALUATOR_HPP

#include <deque>
#include <vector>
#include <mfem.hpp>
#include "fem/output_functionals.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class GridFunction;
class MaterialOperator;
class Mesh;

//
// Class to evaluate derived field quantities (energy densities, Poynting vector) at the
// nodal points of an interpolatory output space using libCEED, filling a grid function
// for visualization output without per-point host coefficient evaluation. Follows the
// conventions of EnergyDensityCoefficient and PoyntingVectorCoefficient.
//
class DomainFieldEvaluator
{
public:
  enum class Kind
  {
    ENERGY_E,  // 1/2 (ε E)ᴴ E, scalar output
    ENERGY_M,  // 1/2 (μ⁻¹ B)ᴴ B, scalar output
    POYNTING   // Re{E x (μ⁻¹ B)⥁}, vector output
  };

private:
  Kind kind;

  // Source field finite element spaces (not owned): nd_fespace for H(curl), rt_fespace
  // for H(div)/L2 B; either may be nullptr depending on the kind.
  const mfem::ParFiniteElementSpace *nd_fespace;
  const mfem::ParFiniteElementSpace *rt_fespace;

  // Whether the evaluator could be assembled.
  bool valid = true;

  // Per-geometry assembled libCEED operators, evaluating at the target space nodal
  // points and scattering directly into the output grid function. The element attribute
  // vectors are operator inputs and must outlive the operators.
  std::vector<fem::CeedGroupOperator> groups;
  std::deque<Vector> elem_attrs;

  // Staging vector used to initialize the field input CeedVectors at construction.
  mutable Vector field_staging;

  void Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                const mfem::ParFiniteElementSpace &target_fespace, double scaling);

public:
  // Construct an evaluator filling grid functions on target_fespace (an interpolatory
  // L2 space). The scaling multiplies the output as for the legacy coefficients.
  DomainFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                       const mfem::ParFiniteElementSpace *nd_fespace,
                       const mfem::ParFiniteElementSpace *rt_fespace,
                       const mfem::ParFiniteElementSpace &target_fespace, double scaling);
  ~DomainFieldEvaluator();

  DomainFieldEvaluator(const DomainFieldEvaluator &) = delete;
  DomainFieldEvaluator &operator=(const DomainFieldEvaluator &) = delete;

  // Whether the evaluator was successfully assembled.
  bool IsValid() const { return valid; }

  // Fill the output vector (L-vector of the target space, e.g. a GridFunction) with
  // the pointwise quantity. Real and imaginary field contributions add. Local
  // operation (no MPI communication).
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;
};

}  // namespace palace

#endif  // PALACE_FEM_DOMAIN_FIELD_EVALUATOR_HPP
