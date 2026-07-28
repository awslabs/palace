// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SINGULARGEOMETRY_HPP
#define PALACE_FEM_SINGULARGEOMETRY_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <mfem.hpp>

#include "fem/singularelements.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

inline double
GetPositiveSignedJacobianDeterminant(mfem::ElementTransformation &transformation)
{
  const auto &jacobian = transformation.Jacobian();
  if (jacobian.Height() != jacobian.Width() || jacobian.Height() < 1)
  {
    throw std::invalid_argument(
        "Singular element geometry requires a nonempty square Jacobian!");
  }

  const double signed_determinant = jacobian.Det();
  const double weight = transformation.Weight();
  const double scale = std::max({1.0, std::abs(signed_determinant), std::abs(weight)});
  if (!std::isfinite(signed_determinant) || !std::isfinite(weight) ||
      !(signed_determinant > 0.0) || !(weight > 0.0) ||
      std::abs(signed_determinant - weight) >
          256.0 * std::numeric_limits<double>::epsilon() * scale)
  {
    throw std::invalid_argument(
        "Singular element geometry requires a positive, orientation-preserving "
        "Jacobian!");
  }
  return signed_determinant;
}

// High-order mesh generators commonly represent an exactly affine simplex with
// a high-order nodal map. Detect that case from the Jacobian polynomial rather
// than from the declared geometry order so it can retain the reference-table
// assembly path.
inline double GetElementTransformationRelativeJacobianVariation(
    mfem::ElementTransformation &transformation)
{
  if (transformation.OrderJ() == 0)
  {
    const auto &rule = mfem::IntRules.Get(transformation.GetGeometryType(), 1);
    if (rule.GetNPoints() < 1)
    {
      throw std::runtime_error(
          "Singular element affine-geometry detection found an empty integration rule!");
    }
    transformation.SetIntPoint(&rule.IntPoint(0));
    const auto &jacobian = transformation.Jacobian();
    if (jacobian.Height() == jacobian.Width())
    {
      (void)GetPositiveSignedJacobianDeterminant(transformation);
    }
    return 0.0;
  }

  const int order =
      std::max(2, 2 * std::max(transformation.Order(), transformation.OrderJ()) + 2);
  const auto &rule = mfem::IntRules.Get(transformation.GetGeometryType(), order);
  if (rule.GetNPoints() < 1)
  {
    throw std::runtime_error(
        "Singular element affine-geometry detection found an empty integration rule!");
  }

  mfem::DenseMatrix reference_jacobian;
  double jacobian_scale = 0.0;
  double maximum_difference = 0.0;
  for (int q = 0; q < rule.GetNPoints(); q++)
  {
    const auto &point = rule.IntPoint(q);
    transformation.SetIntPoint(&point);
    const auto &jacobian = transformation.Jacobian();
    if (jacobian.Height() == jacobian.Width())
    {
      (void)GetPositiveSignedJacobianDeterminant(transformation);
    }
    if (q == 0)
    {
      reference_jacobian = jacobian;
    }
    if (jacobian.Height() != reference_jacobian.Height() ||
        jacobian.Width() != reference_jacobian.Width())
    {
      throw std::runtime_error(
          "Singular element transformation changed Jacobian dimensions!");
    }
    for (int i = 0; i < jacobian.Height(); i++)
    {
      for (int j = 0; j < jacobian.Width(); j++)
      {
        if (!std::isfinite(jacobian(i, j)))
        {
          return std::numeric_limits<double>::infinity();
        }
        jacobian_scale = std::max(
            {jacobian_scale, std::abs(reference_jacobian(i, j)), std::abs(jacobian(i, j))});
        maximum_difference = std::max(maximum_difference,
                                      std::abs(jacobian(i, j) - reference_jacobian(i, j)));
      }
    }
  }
  if (!(jacobian_scale > std::numeric_limits<double>::min()))
  {
    return std::numeric_limits<double>::infinity();
  }
  return maximum_difference / jacobian_scale;
}

inline bool IsAffineElementTransformation(mfem::ElementTransformation &transformation,
                                          double relative_tolerance = 1.0e-9)
{
  if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0)
  {
    throw std::invalid_argument(
        "Singular element affine-geometry detection requires a finite nonnegative "
        "tolerance!");
  }
  return GetElementTransformationRelativeJacobianVariation(transformation) <=
         relative_tolerance;
}

namespace detail
{

inline void SplitBernsteinCoefficients(const mfem::Vector &coefficients, mfem::Vector &left,
                                       mfem::Vector &right)
{
  const int order = coefficients.Size() - 1;
  left.SetSize(order + 1);
  right.SetSize(order + 1);
  mfem::Vector work(coefficients);
  left[0] = work[0];
  right[order] = work[order];
  for (int level = 1; level <= order; level++)
  {
    for (int i = 0; i <= order - level; i++)
    {
      work[i] = 0.5 * (work[i] + work[i + 1]);
    }
    left[level] = work[0];
    right[order - level] = work[order - level];
  }
}

// Bernstein derivative coefficients bound a polynomial derivative over the
// complete parameter interval. Recursive de Casteljau subdivision makes this
// positivity test complete for every derivative bounded strictly away from
// zero; a zero or sign change remains uncertified and is rejected.
inline bool HasStrictlyPositiveBernsteinDerivative(const mfem::Vector &coefficients,
                                                   double tolerance, int depth = 0)
{
  constexpr int maximum_depth = 30;
  const int order = coefficients.Size() - 1;
  if (order < 1)
  {
    return false;
  }

  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < order; i++)
  {
    // Subdivision changes the local parameter interval. Rescale back to the
    // derivative with respect to the original segment coordinate.
    const double derivative =
        std::ldexp(order * (coefficients[i + 1] - coefficients[i]), depth);
    if (!std::isfinite(derivative))
    {
      return false;
    }
    minimum = std::min(minimum, derivative);
    maximum = std::max(maximum, derivative);
  }
  if (minimum > tolerance)
  {
    return true;
  }
  if (maximum <= tolerance || depth == maximum_depth)
  {
    return false;
  }

  mfem::Vector left, right;
  SplitBernsteinCoefficients(coefficients, left, right);
  return HasStrictlyPositiveBernsteinDerivative(left, tolerance, depth + 1) &&
         HasStrictlyPositiveBernsteinDerivative(right, tolerance, depth + 1);
}

}  // namespace detail

// A high-order segment can have a nonconstant parametrization while its physical
// image is still a straight line. This is admissible for a straight Meixner
// feature, provided the map remains regular and consistently oriented. Convert
// the polynomial map to Bernstein control points: their convex-hull property
// bounds the complete physical image, and recursive derivative bounds certify
// regularity over the complete parameter interval rather than only at samples.
inline bool
IsGeometricallyStraightSegmentTransformation(mfem::ElementTransformation &transformation,
                                             double relative_tolerance = 1.0e-10)
{
  if (transformation.GetGeometryType() != mfem::Geometry::SEGMENT ||
      transformation.GetDimension() != 1 || transformation.GetSpaceDim() < 1)
  {
    throw std::invalid_argument(
        "Singular feature geometry requires a physical segment transformation!");
  }
  if (!std::isfinite(relative_tolerance) || relative_tolerance < 0.0)
  {
    throw std::invalid_argument(
        "Singular feature straightness detection requires a finite nonnegative "
        "tolerance!");
  }

  const int space_dimension = transformation.GetSpaceDim();
  mfem::IntegrationPoint point;
  mfem::Vector start(space_dimension), end(space_dimension);
  point.x = 0.0;
  transformation.Transform(point, start);
  point.x = 1.0;
  transformation.Transform(point, end);

  mfem::Vector chord(end);
  chord -= start;
  const double chord_length = chord.Norml2();
  double coordinate_scale = chord_length;
  for (int d = 0; d < space_dimension; d++)
  {
    if (!std::isfinite(start[d]) || !std::isfinite(end[d]))
    {
      return false;
    }
    coordinate_scale = std::max({coordinate_scale, std::abs(start[d]), std::abs(end[d])});
  }
  if (!std::isfinite(chord_length) ||
      !(chord_length > 256.0 * std::numeric_limits<double>::epsilon() *
                           std::max(coordinate_scale, std::numeric_limits<double>::min())))
  {
    return false;
  }

  const int map_order = std::max(1, transformation.Order());
  const int coefficient_count = map_order + 1;
  mfem::DenseMatrix bernstein_values(coefficient_count);
  mfem::DenseMatrix physical_values(coefficient_count, space_dimension);
  mfem::Vector bernstein_shape(coefficient_count);
  mfem::Vector physical_point(space_dimension);
  for (int q = 0; q < coefficient_count; q++)
  {
    point.x = static_cast<double>(q) / map_order;
    mfem::Poly_1D::CalcBernstein(map_order, point.x, bernstein_shape);
    bernstein_values.SetRow(q, bernstein_shape);
    transformation.Transform(point, physical_point);
    for (int d = 0; d < space_dimension; d++)
    {
      if (!std::isfinite(physical_point[d]))
      {
        return false;
      }
      physical_values(q, d) = physical_point[d];
    }
  }
  mfem::DenseMatrix control_points;
  mfem::DenseMatrixInverse(bernstein_values).Mult(physical_values, control_points);
  if (control_points.Height() != coefficient_count ||
      control_points.Width() != space_dimension)
  {
    throw std::runtime_error(
        "Singular feature Bernstein conversion produced invalid dimensions!");
  }

  const double chord_length_squared = chord_length * chord_length;
  mfem::Vector longitudinal_control(coefficient_count);
  mfem::Vector perpendicular(space_dimension);
  for (int i = 0; i < coefficient_count; i++)
  {
    double longitudinal = 0.0;
    for (int d = 0; d < space_dimension; d++)
    {
      if (!std::isfinite(control_points(i, d)))
      {
        return false;
      }
      coordinate_scale = std::max(coordinate_scale, std::abs(control_points(i, d)));
      longitudinal += (control_points(i, d) - start[d]) * chord[d];
    }
    longitudinal /= chord_length_squared;
    longitudinal_control[i] = longitudinal;
    for (int d = 0; d < space_dimension; d++)
    {
      perpendicular[d] = control_points(i, d) - start[d] - longitudinal * chord[d];
    }
    const double distance_tolerance =
        relative_tolerance * chord_length +
        512.0 * std::numeric_limits<double>::epsilon() *
            std::max(coordinate_scale, std::numeric_limits<double>::min());
    if (!std::isfinite(longitudinal) || perpendicular.Norml2() > distance_tolerance)
    {
      return false;
    }
  }

  const double endpoint_tolerance =
      relative_tolerance + 512.0 * std::numeric_limits<double>::epsilon() *
                               std::max({1.0, std::abs(longitudinal_control[0]),
                                         std::abs(longitudinal_control[map_order])});
  if (std::abs(longitudinal_control[0]) > endpoint_tolerance ||
      std::abs(longitudinal_control[map_order] - 1.0) > endpoint_tolerance)
  {
    return false;
  }
  longitudinal_control[0] = 0.0;
  longitudinal_control[map_order] = 1.0;

  double derivative_scale = 1.0;
  for (int i = 0; i < map_order; i++)
  {
    derivative_scale = std::max(
        derivative_scale,
        std::abs(map_order * (longitudinal_control[i + 1] - longitudinal_control[i])));
  }
  const double derivative_tolerance =
      1024.0 * std::numeric_limits<double>::epsilon() * derivative_scale;
  if (!detail::HasStrictlyPositiveBernsteinDerivative(longitudinal_control,
                                                      derivative_tolerance))
  {
    return false;
  }

  // Verify the computed Bernstein representation independently at a rule which
  // is more than sufficient for the declared map order. This turns an
  // ill-conditioned basis conversion into a clean rejection rather than a
  // false geometry classification.
  const auto &rule = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * map_order + 2);
  mfem::Vector reconstructed(space_dimension);
  for (int q = -2; q < rule.GetNPoints(); q++)
  {
    if (q < 0)
    {
      point.x = q == -2 ? 0.0 : 1.0;
    }
    else
    {
      point = rule.IntPoint(q);
    }
    transformation.Transform(point, physical_point);
    mfem::Poly_1D::CalcBernstein(map_order, point.x, bernstein_shape);
    reconstructed = 0.0;
    for (int i = 0; i < coefficient_count; i++)
    {
      for (int d = 0; d < space_dimension; d++)
      {
        reconstructed[d] += bernstein_shape[i] * control_points(i, d);
      }
    }
    reconstructed -= physical_point;
    const double representation_tolerance =
        relative_tolerance * chord_length +
        1024.0 * std::numeric_limits<double>::epsilon() *
            std::max(coordinate_scale, std::numeric_limits<double>::min());
    if (reconstructed.Norml2() > representation_tolerance)
    {
      return false;
    }
    transformation.SetIntPoint(&point);
    const auto &jacobian = transformation.Jacobian();
    if (jacobian.Height() != space_dimension || jacobian.Width() != 1)
    {
      throw std::runtime_error(
          "Singular feature segment transformation changed Jacobian dimensions!");
    }
    for (int d = 0; d < space_dimension; d++)
    {
      if (!std::isfinite(jacobian(d, 0)))
      {
        return false;
      }
    }
  }
  return true;
}

inline BarycentricGradients
GetBarycentricGradients(mfem::ElementTransformation &transformation,
                        const mfem::IntegrationPoint &point, double &jacobian_determinant)
{
  if (transformation.GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
      transformation.GetDimension() != 3 || transformation.GetSpaceDim() != 3)
  {
    throw std::invalid_argument(
        "Singular element geometry requires a three-dimensional tetrahedron!");
  }

  transformation.SetIntPoint(&point);
  jacobian_determinant = GetPositiveSignedJacobianDeterminant(transformation);
  const auto &inverse_jacobian = transformation.InverseJacobian();
  if (inverse_jacobian.Height() != 3 || inverse_jacobian.Width() != 3 ||
      !std::isfinite(jacobian_determinant))
  {
    throw std::invalid_argument(
        "Singular element geometry requires a positive square Jacobian!");
  }

  BarycentricGradients grad_lambda{};
  for (int i = 1; i < 4; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      grad_lambda[i][d] = inverse_jacobian(i - 1, d);
      grad_lambda[0][d] -= grad_lambda[i][d];
    }
  }
  return grad_lambda;
}

inline TriangleBarycentricGradients
GetTriangleBarycentricGradients(mfem::ElementTransformation &transformation,
                                const mfem::IntegrationPoint &point,
                                double &jacobian_determinant)
{
  if (transformation.GetGeometryType() != mfem::Geometry::TRIANGLE ||
      transformation.GetDimension() != 2 || transformation.GetSpaceDim() != 2)
  {
    throw std::invalid_argument(
        "Triangular singular element geometry requires a two-dimensional triangle!");
  }

  transformation.SetIntPoint(&point);
  jacobian_determinant = GetPositiveSignedJacobianDeterminant(transformation);
  const auto &inverse_jacobian = transformation.InverseJacobian();
  if (inverse_jacobian.Height() != 2 || inverse_jacobian.Width() != 2 ||
      !std::isfinite(jacobian_determinant))
  {
    throw std::invalid_argument(
        "Triangular singular element geometry requires a positive square Jacobian!");
  }

  TriangleBarycentricGradients grad_lambda{};
  for (int i = 1; i < 3; i++)
  {
    for (int d = 0; d < 2; d++)
    {
      grad_lambda[i][d] = inverse_jacobian(i - 1, d);
      grad_lambda[0][d] -= grad_lambda[i][d];
    }
  }
  return grad_lambda;
}

}  // namespace singular

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_SINGULARGEOMETRY_HPP
