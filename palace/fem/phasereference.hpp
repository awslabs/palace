// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_PHASE_REFERENCE_HPP
#define PALACE_FEM_PHASE_REFERENCE_HPP

#include <array>
#include <complex>
#include <mfem.hpp>
#include "linalg/vector.hpp"
#include "utils/communication.hpp"

namespace palace
{

class FiniteElementSpace;

//
// Global phase of the eigenvectors of a Nédélec space, independent of the basis of the
// space (the orientations of its edges and faces), the order of the elements and the
// partition: an eigenvector E is rotated so that ∫ w ⋅ E dΩ is real and positive for a
// fixed smooth reference field w, which makes the eigenvector of a lossless mode real, with
// a definite sign. The reference fields are defined in the coordinates ξ ∈ [-1, 1]ᵈ of the
// bounding box of the mesh: w₁ has in each component a polynomial term of every parity
// under the reflections ξₖ → -ξₖ (1, ξₖ, ξₖξₗ, ...), so that symmetric and antisymmetric
// modes alike have a projection on it, and w₂ is a product of sines and cosines. w₂ is used
// when the projection on w₁ is negligible, |∫ w₁ ⋅ E dΩ| ≤ δ ‖b₁‖ ‖E‖ (δ = 10⁻⁶, the norms
// of the degree of freedom vectors with b₁ the linear form of w₁), and the mean of the
// degrees of freedom (the basis dependent convention used before) when both are.
//
class PhaseReference
{
public:
  static constexpr int num_fields = 2;

  // Threshold of the projection on a reference field, relative to the Euclidean norms of
  // the vectors of degrees of freedom: these depend on the mesh and, from the second order,
  // slightly on the orientations of the faces. A mass-weighted (L²) criterion would be
  // independent of the basis, but needs a mass operator, which the eigenmode solver has
  // released before the postprocessing.
  static constexpr double tol = 1.0e-6;

private:
  FiniteElementSpace &fespace;

  // Center and half sizes of the bounding box of the mesh.
  mfem::Vector center, half;

  // Linear forms of the reference fields on the true degrees of freedom, (∫ wₖ ⋅ φᵢ dΩ),
  // and their (global) Euclidean norms, assembled when first needed.
  mutable std::array<Vector, num_fields> b;
  mutable std::array<double, num_fields> b_norm = {};
  mutable std::array<bool, num_fields> assembled = {};

public:
  PhaseReference(FiniteElementSpace &nd_fespace);

  // Values of the reference fields at the point with bounding box coordinates ξ.
  static void EvalField(int k, const double *xi, int dim, double *w);

  // Whether the reference integral f of a vector of norm x_norm, for a reference field with
  // linear form of norm b_norm, sets the phase: |f| > δ ‖b‖ ‖x‖.
  static bool Accept(std::complex<double> f, double b_norm, double x_norm)
  {
    return std::abs(f) > tol * b_norm * x_norm;
  }

  // The linear form of the reference field k on the true degrees of freedom. Collective.
  const Vector &GetLinearForm(int k) const;

  // Rotate x so that the first accepted reference integral, or else the sum of the degrees
  // of freedom, is real and positive. One global reduction when the first reference field
  // is accepted. Returns the reference used, num_fields for the sum of the degrees of
  // freedom, or -1 for a zero vector. Collective.
  int Normalize(ComplexVector &x) const;
};

}  // namespace palace

#endif  // PALACE_FEM_PHASE_REFERENCE_HPP
