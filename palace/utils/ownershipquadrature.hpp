// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_OWNERSHIPQUADRATURE_HPP
#define PALACE_UTILS_OWNERSHIPQUADRATURE_HPP

#include <mfem.hpp>

namespace palace
{

// Nonnegative surface rules for partitioned Gram matrices. MFEM's very-high-order
// simplex rules can have negative weights, which are unsuitable for the square-root
// weighted Gram construction. Symmetrize a Duffy-mapped Gauss product rule over
// three barycentric rotations; the inner Gauss rule supplies reflection symmetry.
class OwnershipQuadrature
{
  mfem::IntegrationRule triangle, square;

public:
  explicit OwnershipQuadrature(int order)
  {
    MFEM_VERIFY(order > 0 && order <= 100, "Invalid ownership quadrature order!");
    int n = (order + 3) / 2;
    // Avoid placing an entire Gauss row on common symmetry/bisector planes.
    n += n % 2;
    const auto &line = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 2 * n - 1);
    const int count = line.GetNPoints();
    triangle.SetSize(3 * count * count);
    square.SetSize(count * count);
    int t = 0, s = 0;
    for (int i = 0; i < count; i++)
    {
      const auto &u = line.IntPoint(i);
      for (int j = 0; j < count; j++)
      {
        const auto &v = line.IntPoint(j);
        auto &quad = square.IntPoint(s++);
        quad.x = u.x;
        quad.y = v.x;
        quad.z = 0.0;
        quad.weight = u.weight * v.weight;
        const double bary[3] = {u.x, (1.0 - u.x) * v.x, (1.0 - u.x) * (1.0 - v.x)};
        for (int rotation = 0; rotation < 3; rotation++)
        {
          auto &point = triangle.IntPoint(t++);
          point.x = bary[rotation];
          point.y = bary[(rotation + 1) % 3];
          point.z = 0.0;
          point.weight = u.weight * v.weight * (1.0 - u.x) / 3.0;
          MFEM_VERIFY(point.weight > 0.0, "Nonpositive ownership quadrature weight!");
        }
      }
    }
    triangle.SetPointIndices();
    square.SetPointIndices();
    triangle.SetOrder(order);
    square.SetOrder(order);
  }

  const mfem::IntegrationRule &Get(mfem::Geometry::Type geometry) const
  {
    MFEM_VERIFY(geometry == mfem::Geometry::TRIANGLE || geometry == mfem::Geometry::SQUARE,
                "Ownership integration requires triangle or quadrilateral faces!");
    return geometry == mfem::Geometry::TRIANGLE ? triangle : square;
  }
};

}  // namespace palace

#endif
