// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_LINE_INTEGRATOR_HPP
#define PALACE_FEM_LINE_INTEGRATOR_HPP

#include <array>
#include <complex>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class GridFunction;

namespace fem
{

//
// Element-aware quadrature for the line integral ∫ F · dl of a finite element vector field
// along a piecewise linear path on a planar mesh: either a 2D mesh, or a planar 2D surface
// mesh embedded in 3D (for example a wave port submesh).
//
// Each path segment is split at every element boundary it crosses, so each sub-interval
// lies inside a single element and is integrated with its own Gauss-Legendre rule. For
// straight-sided elements the integrand is a polynomial on every sub-interval and the
// discrete line integral is computed exactly (to round-off) by a rule of sufficient order,
// independent of how finely the mesh is refined along the path. This is in contrast to a
// single Gauss-Legendre rule over the whole segment, whose error grows as elements become
// small compared to the sample spacing, as happens under adaptive refinement near metal
// edges where the field is singular, and which makes the result change erratically from
// one mesh to the next. For curved (high-order) elements, the crossings with the element
// boundaries are located on the exact curved edges, and the (no longer polynomial)
// integrand on each piece is integrated with the given rule.
//
// Sub-intervals which lie on an edge shared by two elements (a path running along element
// edges) are owned by exactly one element on one rank, so nothing is counted twice. The
// tangential component of an H(curl) field is continuous across such an edge, so the
// choice of owner does not matter. Portions of the path which are not covered by any
// element (for example the interior of an unmeshed conductor) contribute nothing.
//
class LineIntegrator
{
public:
  struct Sample
  {
    int element;                             // Local element index
    mfem::IntegrationPoint point;            // Reference coordinates in the element
    std::array<double, 3> weighted_tangent;  // Quadrature weight times dl/dt
  };

private:
  MPI_Comm comm;
  int sdim;
  const mfem::Mesh *mesh_ptr;
  std::vector<Sample> samples;

  // Total number of sub-intervals (over all segments and ranks), and the length of the path
  // and the portion of it not covered by any mesh element.
  int num_subintervals = 0;
  double path_length = 0.0, uncovered_length = 0.0;

public:
  // Build the quadrature for the path given by consecutive points (each of size equal to
  // the mesh space dimension). When closed is true, a segment from the last point back to
  // the first is added. Each sub-interval is integrated with a Gauss-Legendre rule of the
  // given polynomial order. When min_points_per_segment > 0, every segment is additionally
  // split uniformly so that it receives at least this many quadrature points, which can
  // help for curved (high-order) elements where the integrand is not polynomial.
  LineIntegrator(const mfem::ParMesh &mesh, const std::vector<mfem::Vector> &path,
                 bool closed, int order, int min_points_per_segment = 0);

  // Integrate the given vector field on the mesh along the path. The result is summed over
  // all processes of the mesh communicator (collective).
  double Integrate(const mfem::ParGridFunction &F) const;
  std::complex<double> Integrate(const GridFunction &F) const;

  // Local contribution only (not summed over processes).
  double LocalIntegrate(const mfem::ParGridFunction &F) const;

  const auto &GetSamples() const { return samples; }
  int GetNumSubIntervals() const { return num_subintervals; }
  double GetPathLength() const { return path_length; }
  double GetUncoveredLength() const { return uncovered_length; }
};

}  // namespace fem

}  // namespace palace

#endif  // PALACE_FEM_LINE_INTEGRATOR_HPP
