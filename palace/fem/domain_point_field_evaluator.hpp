// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_DOMAIN_POINT_FIELD_EVALUATOR_HPP
#define PALACE_FEM_DOMAIN_POINT_FIELD_EVALUATOR_HPP

#include <deque>
#include <vector>
#include <mfem.hpp>
#include "fem/ceed_group_operator.hpp"
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
class DomainPointFieldEvaluator
{
public:
  enum class Kind
  {
    FIELD_E,   // E field value, vector output
    FIELD_B,   // B field value, vector output
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
  // points and scattering directly into the output grid function or into a point-major
  // ParaView point buffer. The element attribute vectors are operator inputs and must
  // outlive the operators.
  std::vector<fem::CeedGroupOperator> groups, buffer_groups;
  std::deque<Vector> elem_attrs;

  // Point-major buffer layout for domain visualization fields. Buffer bases are point
  // offsets in the same element/refined-point order used by MFEM's VTU writer.
  int buffer_size = 0, buffer_num_comp = 0;
  std::vector<int> buffer_bases;

  // Staging vector used to initialize the field input CeedVectors at construction.
  mutable Vector field_staging;

  void Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                const mfem::ParFiniteElementSpace &target_fespace, double scaling);

public:
  // Construct an evaluator filling grid functions on target_fespace (an interpolatory
  // L2 space). The scaling multiplies the output as for the legacy coefficients.
  DomainPointFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                       const mfem::ParFiniteElementSpace *nd_fespace,
                       const mfem::ParFiniteElementSpace *rt_fespace,
                       const mfem::ParFiniteElementSpace &target_fespace, double scaling);
  ~DomainPointFieldEvaluator();

  DomainPointFieldEvaluator(const DomainPointFieldEvaluator &) = delete;
  DomainPointFieldEvaluator &operator=(const DomainPointFieldEvaluator &) = delete;

  // Whether the evaluator was successfully assembled.
  bool IsValid() const { return valid; }

  // Total buffer size (all elements, lattice points, components), number of components,
  // and per-element point-base offsets for the domain visualization field. Vector
  // buffers are point-major: xyzxyz... in VTK tuple order.
  int BufferSize() const { return buffer_size; }
  int BufferNumComp() const { return buffer_num_comp; }
  const std::vector<int> &BufferBases() const { return buffer_bases; }

  // Fill the output vector (L-vector of the target space, e.g. a GridFunction) with
  // the pointwise quantity. Real and imaginary field contributions add. Local
  // operation (no MPI communication).
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;

  // Fill the point-major domain visualization point buffer for a single linear field.
  // Local operation (no MPI communication).
  void EvalBuffer(const Vector &u, Vector &buffer) const;

  // Fill the point-major domain visualization point buffer. Real and imaginary field
  // contributions add. Local operation (no MPI communication).
  void EvalBuffer(const GridFunction *E, const GridFunction *B, Vector &buffer) const;
};

}  // namespace palace

#endif  // PALACE_FEM_DOMAIN_POINT_FIELD_EVALUATOR_HPP
