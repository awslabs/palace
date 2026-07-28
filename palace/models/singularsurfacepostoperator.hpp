// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SINGULAR_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SINGULAR_SURFACE_POST_OPERATOR_HPP

#include <map>
#include <vector>
#include <mfem.hpp>

#include "fem/singularassembly.hpp"
#include "utils/labels.hpp"

namespace palace
{

class MaterialOperator;

namespace config
{

struct BoundaryPostData;

}  // namespace config

namespace fem::singular
{

class EnrichedH1FieldEvaluator;
class EnrichedNDFieldEvaluator;
class TriangleEnrichedNDFieldEvaluator;
class TriangleEnrichedH1FieldEvaluator;

}  // namespace fem::singular

// Surface dielectric participation for a combined standard-plus-singular
// triangular H(curl) field. This operator evaluates traces from adjacent
// volume elements and uses endpoint-weighted quadrature for integrable
// Meixner corner singularities.
class TriangleSingularSurfacePostOperator
{
public:
  struct Measurement
  {
    int index;
    double energy;
    double loss_tangent;
    double participation;
    double quality_factor;
  };

private:
  struct InterfaceData
  {
    InterfaceDielectric type;
    double thickness;
    double permittivity;
    double loss_tangent;
    double edge_cutoff;
    mfem::Array<int> attribute_marker;
  };

  const MaterialOperator &material;
  mfem::ParFiniteElementSpace &fespace;
  std::map<int, InterfaceData> interfaces;

  double IntegrateInterface(
      const InterfaceData &interface,
      fem::singular::TriangleEnrichedNDFieldEvaluator *real_evaluator,
      fem::singular::TriangleEnrichedNDFieldEvaluator *imaginary_evaluator,
      fem::singular::TriangleEnrichedH1FieldEvaluator *real_gradient_evaluator,
      fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_gradient_evaluator,
      fem::singular::TriangleEnrichedH1FieldEvaluator *real_longitudinal_evaluator,
      fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_longitudinal_evaluator,
      const fem::singular::AdaptiveAssemblyOptions &options) const;

public:
  TriangleSingularSurfacePostOperator(const config::BoundaryPostData &postpro,
                                      const MaterialOperator &material,
                                      mfem::ParFiniteElementSpace &fespace);

  bool Empty() const { return interfaces.empty(); }

  std::vector<Measurement> Measure(
      fem::singular::TriangleEnrichedNDFieldEvaluator &real_evaluator,
      fem::singular::TriangleEnrichedNDFieldEvaluator &imaginary_evaluator,
      double total_electric_energy, const fem::singular::AdaptiveAssemblyOptions &options,
      fem::singular::TriangleEnrichedH1FieldEvaluator *real_longitudinal_evaluator =
          nullptr,
      fem::singular::TriangleEnrichedH1FieldEvaluator *imaginary_longitudinal_evaluator =
          nullptr) const;

  std::vector<Measurement>
  MeasureElectrostatic(fem::singular::TriangleEnrichedH1FieldEvaluator &real_evaluator,
                       fem::singular::TriangleEnrichedH1FieldEvaluator &imaginary_evaluator,
                       double total_electric_energy,
                       const fem::singular::AdaptiveAssemblyOptions &options) const;
};

// Three-dimensional electrostatic surface participation evaluated from the
// combined standard-plus-singular H1 gradient. Singular edge charts use
// Gauss-Jacobi quadrature when the trace is integrable and a logarithmic radial
// map when an explicit physical cutoff is present.
class TetrahedronSingularSurfacePostOperator
{
public:
  using Measurement = TriangleSingularSurfacePostOperator::Measurement;

private:
  struct InterfaceData
  {
    InterfaceDielectric type;
    double thickness;
    double permittivity;
    double loss_tangent;
    double edge_cutoff;
    mfem::Array<int> attribute_marker;
  };

  const MaterialOperator &material;
  mfem::ParFiniteElementSpace &fespace;
  std::map<int, InterfaceData> interfaces;

  double
  IntegrateInterface(const InterfaceData &interface,
                     fem::singular::EnrichedNDFieldEvaluator *real_evaluator,
                     fem::singular::EnrichedNDFieldEvaluator *imaginary_evaluator,
                     fem::singular::EnrichedH1FieldEvaluator *real_gradient_evaluator,
                     fem::singular::EnrichedH1FieldEvaluator *imaginary_gradient_evaluator,
                     const fem::singular::AdaptiveAssemblyOptions &options) const;

public:
  TetrahedronSingularSurfacePostOperator(const config::BoundaryPostData &postpro,
                                         const MaterialOperator &material,
                                         mfem::ParFiniteElementSpace &fespace);

  bool Empty() const { return interfaces.empty(); }

  std::vector<Measurement>
  Measure(fem::singular::EnrichedNDFieldEvaluator &real_evaluator,
          fem::singular::EnrichedNDFieldEvaluator &imaginary_evaluator,
          double total_electric_energy,
          const fem::singular::AdaptiveAssemblyOptions &options) const;

  std::vector<Measurement>
  MeasureElectrostatic(fem::singular::EnrichedH1FieldEvaluator &real_evaluator,
                       fem::singular::EnrichedH1FieldEvaluator &imaginary_evaluator,
                       double total_electric_energy,
                       const fem::singular::AdaptiveAssemblyOptions &options) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SINGULAR_SURFACE_POST_OPERATOR_HPP
