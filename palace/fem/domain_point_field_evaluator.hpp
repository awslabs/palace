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
// Class to evaluate primary and derived fields using libCEED, either into an
// interpolatory output grid function or directly into MFEM's refined VTU point order.
//
class DomainPointFieldEvaluator
{
public:
  enum class Kind
  {
    FIELD_E,   // H(curl) vector field value
    FIELD_B,   // Magnetic flux value, scalar in 2D and vector in 3D
    FIELD_H1,  // H1 scalar field value
    ENERGY_E,  // 1/2 (ε E)ᴴ E, scalar output
    ENERGY_M,  // 1/2 (μ⁻¹ B)ᴴ B, scalar output
    POYNTING,  // Re{E x (μ⁻¹ B)⥁}, vector output
    MODE_SN    // Boundary-mode Re{Eₜ x (μ⁻¹Bₜ)⋆} ⋅ z, scalar output
  };

private:
  Kind kind;

  // Source field finite element spaces (not owned): nd_fespace for H(curl), rt_fespace
  // for H(div)/L2 B; either may be nullptr depending on the kind.
  const mfem::ParFiniteElementSpace *nd_fespace;
  const mfem::ParFiniteElementSpace *rt_fespace;

  // Per-geometry assembled libCEED operators for grid-function and direct point-buffer
  // output. The element attribute vectors are operator inputs and must outlive them.
  std::vector<fem::CeedGroupOperator> groups, buffer_groups;
  std::deque<Vector> elem_attrs;

  // Point-major VTU buffer layout: element order, refined point order, component.
  int buffer_size = 0, buffer_num_comp = 0;
  std::vector<int> buffer_bases;

  // Staging vector used to initialize the field input CeedVectors at construction.
  mutable Vector field_staging;

  void Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                const mfem::ParFiniteElementSpace &target_fespace, double scaling,
                bool build_gridfunction, bool build_buffer);

public:
  // Construct an evaluator filling grid functions on target_fespace (an interpolatory
  // L2 space). The scaling multiplies the output as for the legacy coefficients.
  DomainPointFieldEvaluator(Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
                            const mfem::ParFiniteElementSpace *nd_fespace,
                            const mfem::ParFiniteElementSpace *rt_fespace,
                            const mfem::ParFiniteElementSpace &target_fespace,
                            double scaling, bool build_gridfunction = true,
                            bool build_buffer = false);
  ~DomainPointFieldEvaluator();

  DomainPointFieldEvaluator(const DomainPointFieldEvaluator &) = delete;
  DomainPointFieldEvaluator &operator=(const DomainPointFieldEvaluator &) = delete;

  // Fill the output vector (L-vector of the target space, e.g. a GridFunction) with
  // the pointwise quantity. Real and imaginary field contributions add. Local
  // operation (no MPI communication).
  void Eval(const GridFunction *E, const GridFunction *B, Vector &out) const;

  int BufferSize() const { return buffer_size; }
  int BufferNumComp() const { return buffer_num_comp; }

  // Fill a point-major domain visualization buffer for a primary field.
  void EvalBuffer(const Vector &u, Vector &buffer) const;

  // Fill a point-major domain visualization buffer for a derived field.
  void EvalBuffer(const GridFunction *E, const GridFunction *B, Vector &buffer) const;
};

}  // namespace palace

#endif  // PALACE_FEM_DOMAIN_POINT_FIELD_EVALUATOR_HPP
