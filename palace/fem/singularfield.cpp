// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularfield.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

void ValidateInteriorPoint(const mfem::IntegrationPoint &point)
{
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  if (!std::all_of(lambda.begin(), lambda.end(),
                   [](double value) { return std::isfinite(value) && value > 0.0; }))
  {
    throw std::invalid_argument(
        "Singular H1 field evaluation requires a strictly interior tetrahedron point!");
  }
}

double Dot(const Vector3 &left, const Vector3 &right)
{
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

double DistanceToLine(const mfem::Vector &point, const double *first, const double *second)
{
  Vector3 direction, offset;
  for (int d = 0; d < 3; d++)
  {
    direction[d] = second[d] - first[d];
    offset[d] = point[d] - first[d];
  }
  const Vector3 cross{offset[1] * direction[2] - offset[2] * direction[1],
                      offset[2] * direction[0] - offset[0] * direction[2],
                      offset[0] * direction[1] - offset[1] * direction[0]};
  const double direction_norm = std::sqrt(Dot(direction, direction));
  if (!std::isfinite(direction_norm) || !(direction_norm > 0.0))
  {
    throw std::runtime_error("Singular edge diagnostic found a degenerate mesh edge!");
  }
  return std::sqrt(Dot(cross, cross)) / direction_norm;
}

class PotentialCoefficient final : public mfem::Coefficient
{
private:
  EnrichedH1FieldEvaluator &evaluator;

public:
  explicit PotentialCoefficient(EnrichedH1FieldEvaluator &evaluator) : evaluator(evaluator)
  {
  }

  double Eval(mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &point) override
  {
    return evaluator.Evaluate(transformation.ElementNo, point).potential;
  }
};

class ElectricFieldCoefficient final : public mfem::VectorCoefficient
{
private:
  EnrichedH1FieldEvaluator &evaluator;

public:
  explicit ElectricFieldCoefficient(EnrichedH1FieldEvaluator &evaluator)
    : mfem::VectorCoefficient(3), evaluator(evaluator)
  {
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    const auto field = evaluator.Evaluate(transformation.ElementNo, point);
    value.SetSize(3);
    for (int d = 0; d < 3; d++)
    {
      value[d] = -field.gradient[d];
    }
  }
};

}  // namespace

H1FieldValue EvaluateElementH1Enrichment(const ElementDofMap &element_dofs,
                                         const mfem::Vector &local_coefficients,
                                         const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda)
{
  H1FieldValue result;
  for (const auto &element_dof : element_dofs.h1)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()) ||
        !IsGradientFamily(element_dof.basis.family))
    {
      throw std::invalid_argument(
          "Singular H1 field evaluation received an invalid element DOF map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Singular H1 field evaluation received a nonfinite coefficient!");
    }
    const auto basis = EvaluateHigherOrderBasis(lambda, grad_lambda, element_dof.basis);
    result.potential +=
        coefficient * EvaluateHigherOrderGradientPotential(lambda, element_dof.basis);
    for (int d = 0; d < 3; d++)
    {
      result.gradient[d] += coefficient * basis.value[d];
    }
  }
  return result;
}

EnrichedH1FieldEvaluator::EnrichedH1FieldEvaluator(const DofTopology &topology,
                                                   const ParallelDofNumbering &numbering,
                                                   mfem::ParFiniteElementSpace &fespace)
  : topology(topology), numbering(numbering), fespace(fespace), standard_field(&fespace),
    enrichment_prolongation(
        BuildParallelEnrichmentProlongation(fespace.GetComm(), numbering.h1)),
    h1_exponents(topology.h1_dofs.size(), std::numeric_limits<double>::quiet_NaN()),
    initialized(false)
{
  bool valid =
      fespace.GetParMesh() &&
      topology.elements.size() == static_cast<std::size_t>(fespace.GetNE()) &&
      topology.h1_dofs.size() == numbering.h1.local_to_true.size() &&
      numbering.h1.local_size == static_cast<HYPRE_BigInt>(topology.h1_dofs.size()) &&
      numbering.h1.owned_size <=
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      enrichment_prolongation &&
      enrichment_prolongation->Height() == static_cast<int>(topology.h1_dofs.size()) &&
      enrichment_prolongation->Width() == static_cast<int>(numbering.h1.owned_size);
  if (valid)
  {
    for (const auto &element : topology.elements)
    {
      for (const auto &element_dof : element.h1)
      {
        if (element_dof.dof >= h1_exponents.size() ||
            !IsGradientFamily(element_dof.basis.family) ||
            topology.h1_dofs[element_dof.dof].family != element_dof.basis.family ||
            topology.h1_dofs[element_dof.dof].order != element_dof.basis.order ||
            !std::isfinite(element_dof.basis.nu) || !(element_dof.basis.nu > 0.0) ||
            !(element_dof.basis.nu < 1.0))
        {
          valid = false;
          break;
        }
        double &exponent = h1_exponents[element_dof.dof];
        if (std::isnan(exponent))
        {
          exponent = element_dof.basis.nu;
        }
        else if (exponent != element_dof.basis.nu)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
    valid = valid && std::all_of(h1_exponents.begin(), h1_exponents.end(),
                                 [](double exponent) { return std::isfinite(exponent); });
  }
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Singular H1 field evaluator received inconsistent space and DOF topology!");
  }
  local_enrichment.SetSize(static_cast<int>(topology.h1_dofs.size()));
}

void EnrichedH1FieldEvaluator::SetFromTrueDofs(const mfem::Vector &combined_true_dofs)
{
  const int standard_size = fespace.GetTrueVSize();
  const int enrichment_size = static_cast<int>(numbering.h1.owned_size);
  bool valid = combined_true_dofs.Size() == standard_size + enrichment_size;
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Combined singular H1 vector has inconsistent process-local dimensions!");
  }

  mfem::Vector standard_true_dofs(standard_size);
  mfem::Vector enrichment_true_dofs(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard_true_dofs[i] = combined_true_dofs[i];
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_true_dofs[i] = combined_true_dofs[standard_size + i];
  }
  standard_field.SetFromTrueDofs(standard_true_dofs);
  enrichment_prolongation->Mult(enrichment_true_dofs, local_enrichment);
  initialized = true;
}

H1FieldValue EnrichedH1FieldEvaluator::Evaluate(int element,
                                                const mfem::IntegrationPoint &point)
{
  if (!initialized)
  {
    throw std::logic_error("Singular H1 field evaluator has no combined true-DOF vector!");
  }
  if (element < 0 || element >= fespace.GetNE())
  {
    throw std::out_of_range("Singular H1 field evaluation element is out of range!");
  }
  ValidateInteriorPoint(point);

  auto *transformation = fespace.GetElementTransformation(element);
  if (!transformation)
  {
    throw std::runtime_error("Singular H1 field element has no transformation!");
  }
  double jacobian_determinant;
  const auto grad_lambda =
      GetAffineBarycentricGradients(*transformation, jacobian_determinant);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};

  transformation->SetIntPoint(&point);
  H1FieldValue result;
  result.potential = standard_field.GetValue(*transformation, point);
  mfem::Vector standard_gradient(3);
  standard_field.GetGradient(*transformation, standard_gradient);
  for (int d = 0; d < 3; d++)
  {
    result.gradient[d] = standard_gradient[d];
  }

  const auto enrichment = EvaluateElementH1Enrichment(
      topology.elements[element], local_enrichment, lambda, grad_lambda);
  result.potential += enrichment.potential;
  for (int d = 0; d < 3; d++)
  {
    result.gradient[d] += enrichment.gradient[d];
  }
  return result;
}

AdaptiveQuadratureResult EnrichedH1FieldEvaluator::IntegrateElementGradientEnergy(
    int element, double electric_coefficient, const AdaptiveAssemblyOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error("Singular H1 energy evaluator has no combined true-DOF vector!");
  }
  if (!std::isfinite(electric_coefficient) || !(electric_coefficient > 0.0))
  {
    throw std::invalid_argument(
        "Singular H1 energy integration requires a positive electric coefficient!");
  }
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Singular H1 energy integration has invalid adaptive quadrature options!");
  }
  if (element < 0 || element >= fespace.GetNE())
  {
    throw std::out_of_range("Singular H1 energy integration element is out of range!");
  }
  auto *transformation = fespace.GetElementTransformation(element);
  if (!transformation)
  {
    throw std::runtime_error("Singular H1 energy element has no transformation!");
  }
  double jacobian_determinant;
  GetAffineBarycentricGradients(*transformation, jacobian_determinant);

  if (topology.elements[element].h1.empty())
  {
    const auto *finite_element = fespace.GetFE(element);
    if (!finite_element || finite_element->GetGeomType() != mfem::Geometry::TETRAHEDRON ||
        finite_element->GetMapType() != mfem::FiniteElement::VALUE)
    {
      throw std::runtime_error(
          "Singular H1 energy integration requires a tetrahedral value element!");
    }
    const int integration_order = 2 * std::max(0, finite_element->GetOrder() - 1);
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, integration_order);
    mfem::Vector gradient(3);
    double value = 0.0;
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &point = rule.IntPoint(q);
      transformation->SetIntPoint(&point);
      standard_field.GetGradient(*transformation, gradient);
      value += point.weight * electric_coefficient * jacobian_determinant *
               (gradient * gradient);
    }
    return {value, 0.0, 1, 0, true};
  }

  mfem::IntegrationPoint point;
  return IntegrateReferenceTetrahedronAdaptive(
      options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
      options.maximum_subdivisions,
      [&](const BarycentricPoint &lambda)
      {
        point.Set3(lambda[1], lambda[2], lambda[3]);
        const auto value = Evaluate(element, point);
        return electric_coefficient * jacobian_determinant *
               Dot(value.gradient, value.gradient);
      });
}

void EnrichedH1FieldEvaluator::ProjectToDiscontinuousGridFunctions(
    mfem::ParGridFunction &potential, mfem::ParGridFunction &electric_field)
{
  const auto *potential_space = potential.ParFESpace();
  const auto *electric_space = electric_field.ParFESpace();
  if (!initialized || !potential_space || !electric_space ||
      potential_space->GetParMesh() != fespace.GetParMesh() ||
      electric_space->GetParMesh() != fespace.GetParMesh() ||
      potential_space->GetVDim() != 1 || electric_space->GetVDim() != 3 ||
      potential_space->FEColl()->GetMapType(3) != mfem::FiniteElement::VALUE ||
      electric_space->FEColl()->GetMapType(3) != mfem::FiniteElement::VALUE)
  {
    throw std::invalid_argument(
        "Singular field sampling requires initialized scalar and vector value spaces "
        "on the solve mesh!");
  }

  PotentialCoefficient potential_coefficient(*this);
  ElectricFieldCoefficient electric_coefficient(*this);
  potential.ProjectCoefficient(potential_coefficient);
  electric_field.ProjectCoefficient(electric_coefficient);
}

std::vector<H1CoefficientDiagnostic>
EnrichedH1FieldEvaluator::GetOwnedCoefficientDiagnostics() const
{
  if (!initialized)
  {
    throw std::logic_error(
        "Singular H1 coefficient diagnostics require a combined true-DOF vector!");
  }

  std::vector<H1CoefficientDiagnostic> diagnostics;
  diagnostics.reserve(static_cast<std::size_t>(numbering.h1.owned_size));
  const int rank = Mpi::Rank(fespace.GetComm());
  for (std::size_t local = 0; local < topology.h1_dofs.size(); local++)
  {
    if (numbering.h1.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.h1.local_to_true[local];
    if (true_dof < numbering.h1.owned_offset ||
        true_dof >= numbering.h1.owned_offset + numbering.h1.owned_size ||
        !std::isfinite(local_enrichment[static_cast<int>(local)]) ||
        !std::isfinite(h1_exponents[local]))
    {
      throw std::runtime_error(
          "Singular H1 coefficient diagnostics found inconsistent owned data!");
    }
    diagnostics.push_back({true_dof, topology.h1_dofs[local], h1_exponents[local],
                           local_enrichment[static_cast<int>(local)]});
  }
  std::sort(diagnostics.begin(), diagnostics.end(), [](const auto &left, const auto &right)
            { return left.true_dof < right.true_dof; });
  if (diagnostics.size() != static_cast<std::size_t>(numbering.h1.owned_size))
  {
    throw std::runtime_error(
        "Singular H1 coefficient diagnostics do not cover every owned true DOF!");
  }
  return diagnostics;
}

std::vector<H1EdgeSlopeDiagnostic> EnrichedH1FieldEvaluator::FitEdgeSlopes(
    const FeatureTopology &features, const std::vector<GlobalVertexId> &source_vertex_ids,
    const std::vector<GlobalVertexId> &source_element_ids, const EdgeSlopeOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Singular edge slope diagnostics require a combined true-DOF vector!");
  }
  if (features.elements.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      source_vertex_ids.size() != static_cast<std::size_t>(fespace.GetParMesh()->GetNV()) ||
      source_element_ids.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      options.sample_count < 3 || !std::isfinite(options.minimum_barycentric_radius) ||
      !std::isfinite(options.maximum_barycentric_radius) ||
      !(options.minimum_barycentric_radius > 0.0) ||
      !(options.minimum_barycentric_radius < options.maximum_barycentric_radius) ||
      !(options.maximum_barycentric_radius < 1.0))
  {
    throw std::invalid_argument(
        "Singular edge slope diagnostics received inconsistent topology or options!");
  }

  const double log_minimum_radius = std::log(options.minimum_barycentric_radius);
  const double log_maximum_radius = std::log(options.maximum_barycentric_radius);
  std::vector<H1EdgeSlopeDiagnostic> diagnostics;
  mfem::IntegrationPoint point;
  mfem::Vector physical_point(3);
  for (int element = 0; element < fespace.GetNE(); element++)
  {
    const auto *tetrahedron = fespace.GetParMesh()->GetElement(element);
    auto *transformation = fespace.GetElementTransformation(element);
    if (!tetrahedron || tetrahedron->GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
        !transformation)
    {
      throw std::runtime_error(
          "Singular edge slope diagnostics require affine tetrahedra!");
    }
    for (const auto &edge : features.elements[element].edges)
    {
      if (edge.feature >= features.features.size() ||
          edge.segment >= features.segments.size() ||
          features.segments[edge.segment].feature != edge.feature)
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostics found invalid feature incidence!");
      }
      std::array<bool, 4> seen{false, false, false, false};
      std::array<GlobalVertexId, 4> canonical_vertices;
      for (int local = 0; local < 4; local++)
      {
        const int node = edge.canonical_nodes[local];
        if (node < 0 || node >= 4 || seen[node])
        {
          throw std::invalid_argument(
              "Singular edge slope diagnostics found invalid canonical nodes!");
        }
        seen[node] = true;
        const int mesh_vertex = tetrahedron->GetVertices()[node];
        if (mesh_vertex < 0 || mesh_vertex >= static_cast<int>(source_vertex_ids.size()))
        {
          throw std::invalid_argument(
              "Singular edge slope diagnostics found invalid mesh vertices!");
        }
        canonical_vertices[local] = source_vertex_ids[mesh_vertex];
      }
      std::array<GlobalVertexId, 2> incidence_edge{canonical_vertices[0],
                                                   canonical_vertices[1]};
      std::array<GlobalVertexId, 2> source_edge{
          features.segments[edge.segment].mesh_vertices[0],
          features.segments[edge.segment].mesh_vertices[1]};
      std::sort(incidence_edge.begin(), incidence_edge.end());
      std::sort(source_edge.begin(), source_edge.end());
      if (incidence_edge != source_edge)
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostic incidence does not match its source segment!");
      }

      const int first_mesh_vertex = tetrahedron->GetVertices()[edge.canonical_nodes[0]];
      const int second_mesh_vertex = tetrahedron->GetVertices()[edge.canonical_nodes[1]];
      const double *first = fespace.GetParMesh()->GetVertex(first_mesh_vertex);
      const double *second = fespace.GetParMesh()->GetVertex(second_mesh_vertex);
      std::vector<double> log_distance(options.sample_count);
      std::vector<double> log_field_norm(options.sample_count);
      std::vector<double> distance(options.sample_count);
      std::vector<double> field_norm(options.sample_count);
      bool valid = true;
      for (int sample = 0; sample < options.sample_count; sample++)
      {
        const double fraction =
            static_cast<double>(sample) / static_cast<double>(options.sample_count - 1);
        const double radius = std::exp(
            log_minimum_radius + fraction * (log_maximum_radius - log_minimum_radius));
        BarycentricPoint lambda;
        lambda[edge.canonical_nodes[0]] = 0.5 * (1.0 - radius);
        lambda[edge.canonical_nodes[1]] = 0.5 * (1.0 - radius);
        lambda[edge.canonical_nodes[2]] = 0.5 * radius;
        lambda[edge.canonical_nodes[3]] = 0.5 * radius;
        point.Set3(lambda[1], lambda[2], lambda[3]);
        const auto field = Evaluate(element, point);
        transformation->Transform(point, physical_point);
        distance[sample] = DistanceToLine(physical_point, first, second);
        field_norm[sample] = std::sqrt(Dot(field.gradient, field.gradient));
        if (!std::isfinite(distance[sample]) || !(distance[sample] > 0.0) ||
            !std::isfinite(field_norm[sample]) || !(field_norm[sample] > 0.0))
        {
          valid = false;
          break;
        }
        log_distance[sample] = std::log(distance[sample]);
        log_field_norm[sample] = std::log(field_norm[sample]);
      }

      double fitted_slope = 0.0;
      double r_squared = 0.0;
      if (valid)
      {
        long double mean_distance = 0.0L;
        long double mean_field = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          mean_distance += log_distance[sample];
          mean_field += log_field_norm[sample];
        }
        mean_distance /= options.sample_count;
        mean_field /= options.sample_count;
        long double distance_variance = 0.0L;
        long double covariance = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          const long double centered_distance = log_distance[sample] - mean_distance;
          distance_variance += centered_distance * centered_distance;
          covariance += centered_distance * (log_field_norm[sample] - mean_field);
        }
        if (!(distance_variance > 0.0L))
        {
          valid = false;
        }
        else
        {
          fitted_slope = static_cast<double>(covariance / distance_variance);
          const long double intercept = mean_field - fitted_slope * mean_distance;
          long double residual_sum = 0.0L;
          long double total_sum = 0.0L;
          for (int sample = 0; sample < options.sample_count; sample++)
          {
            const long double residual =
                log_field_norm[sample] - (intercept + fitted_slope * log_distance[sample]);
            const long double centered = log_field_norm[sample] - mean_field;
            residual_sum += residual * residual;
            total_sum += centered * centered;
          }
          r_squared = (total_sum > 0.0L)
                          ? static_cast<double>(1.0L - residual_sum / total_sum)
                          : ((residual_sum == 0.0L) ? 1.0 : 0.0);
          valid = std::isfinite(fitted_slope) && std::isfinite(r_squared);
        }
      }

      const double exponent = features.features[edge.feature].nu;
      if (!std::isfinite(exponent) || !(exponent > 0.0) || !(exponent < 1.0))
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostics found an invalid feature exponent!");
      }
      diagnostics.push_back(
          {source_element_ids[element], edge.feature, edge.segment, canonical_vertices,
           options.sample_count, exponent, exponent - 1.0, fitted_slope, r_squared,
           valid ? distance.front() : 0.0, valid ? distance.back() : 0.0,
           valid ? field_norm.front() : 0.0, valid ? field_norm.back() : 0.0, valid});
    }
  }
  return diagnostics;
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
