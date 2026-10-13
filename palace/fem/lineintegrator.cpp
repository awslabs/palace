// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lineintegrator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include "fem/gridfunction.hpp"
#include "utils/communication.hpp"

namespace palace::fem
{

namespace
{

// Tolerances, relative to the parameter t ∈ [0, 1] along a segment: breakpoints closer than
// this are merged, and sub-intervals shorter than this are ignored.
constexpr double T_TOL = 1.0e-12;

// Relative tolerances for geometric predicates, relative to the local element size.
constexpr double PARALLEL_TOL = 1.0e-10;  // Segment parallel to an element edge
constexpr double ON_EDGE_TOL = 1.0e-8;    // Segment lying on an element edge
constexpr double PLANE_TOL = 1.0e-4;      // Segment lying in the plane of an element

using Vec3 = std::array<double, 3>;

inline double Dot(const Vec3 &a, const Vec3 &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Vec3 Sub(const Vec3 &a, const Vec3 &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Vec3 Cross(const Vec3 &a, const Vec3 &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

inline double Norm(const Vec3 &a)
{
  return std::sqrt(Dot(a, a));
}

// Geometry of a straight-sided (planar, convex) polygonal element used for clipping. For
// curved elements the vertices are exact and the edges are approximated as straight.
struct ElementPolygon
{
  int nv = 0;
  std::array<Vec3, 4> x;  // Vertices
  std::array<Vec3, 4> n;  // Unit in-plane outward normals of edges (x[i], x[i + 1])
  Vec3 normal;            // Unit element normal (z for 2D meshes)
  Vec3 bb_min, bb_max;    // Bounding box
  double h = 0.0;         // Maximum edge length
};

ElementPolygon GetElementPolygon(mfem::ParMesh &mesh, int e)
{
  const int sdim = mesh.SpaceDimension();
  const auto geom = mesh.GetElementBaseGeometry(e);
  MFEM_VERIFY(geom == mfem::Geometry::TRIANGLE || geom == mfem::Geometry::SQUARE,
              "LineIntegrator only supports triangle and quadrilateral elements!");
  ElementPolygon poly;
  poly.nv = mfem::Geometry::NumVerts[geom];
  const auto *ref_verts = mfem::Geometries.GetVertices(geom);
  auto &T = *mesh.GetElementTransformation(e);
  mfem::Vector xv(sdim);
  for (int i = 0; i < poly.nv; i++)
  {
    T.Transform(ref_verts->IntPoint(i), xv);
    poly.x[i] = {xv(0), xv(1), (sdim == 3) ? xv(2) : 0.0};
  }
  for (int d = 0; d < 3; d++)
  {
    poly.bb_min[d] = poly.bb_max[d] = poly.x[0][d];
    for (int i = 1; i < poly.nv; i++)
    {
      poly.bb_min[d] = std::min(poly.bb_min[d], poly.x[i][d]);
      poly.bb_max[d] = std::max(poly.bb_max[d], poly.x[i][d]);
    }
  }

  // Element normal from the diagonals (quadrilateral) or two edges (triangle).
  if (sdim == 2)
  {
    poly.normal = {0.0, 0.0, 1.0};
  }
  else
  {
    const Vec3 a = (poly.nv == 4) ? Sub(poly.x[2], poly.x[0]) : Sub(poly.x[1], poly.x[0]);
    const Vec3 b = (poly.nv == 4) ? Sub(poly.x[3], poly.x[1]) : Sub(poly.x[2], poly.x[0]);
    poly.normal = Cross(a, b);
    const double s = Norm(poly.normal);
    MFEM_VERIFY(s > 0.0, "Degenerate element in LineIntegrator!");
    for (auto &c : poly.normal)
    {
      c /= s;
    }
  }

  // In-plane outward edge normals.
  Vec3 c = {0.0, 0.0, 0.0};
  for (int i = 0; i < poly.nv; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      c[d] += poly.x[i][d] / poly.nv;
    }
  }
  for (int i = 0; i < poly.nv; i++)
  {
    const Vec3 t = Sub(poly.x[(i + 1) % poly.nv], poly.x[i]);
    poly.h = std::max(poly.h, Norm(t));
    Vec3 n = Cross(t, poly.normal);
    const double s = Norm(n);
    MFEM_VERIFY(s > 0.0, "Degenerate element edge in LineIntegrator!");
    for (auto &v : n)
    {
      v /= s;
    }
    if (Dot(n, Sub(c, poly.x[i])) > 0.0)
    {
      for (auto &v : n)
      {
        v = -v;
      }
    }
    poly.n[i] = n;
  }
  return poly;
}

// Clip the segment x(t) = p0 + t d, t ∈ [0, 1], against a convex element polygon. Returns
// false if the intersection is empty or degenerate, otherwise the parameter interval
// [t_lo, t_hi]. Segments lying on an element edge (within a small tolerance) are considered
// to intersect the element, so both neighbors of the edge report the interval and the
// caller resolves ownership.
bool ClipSegment(const ElementPolygon &poly, const Vec3 &p0, const Vec3 &d, double &t_lo,
                 double &t_hi)
{
  const double len = Norm(d);
  if (len == 0.0)
  {
    return false;
  }

  // The segment must lie in the plane of the element (only relevant for surface meshes).
  {
    const double d0 = Dot(poly.normal, Sub(p0, poly.x[0]));
    const double d1 = d0 + Dot(poly.normal, d);
    if (std::abs(d0) > PLANE_TOL * poly.h || std::abs(d1) > PLANE_TOL * poly.h)
    {
      return false;
    }
  }

  t_lo = 0.0;
  t_hi = 1.0;
  for (int i = 0; i < poly.nv; i++)
  {
    // Inside: n · (p0 + t d - x_i) ≤ 0  ⇔  t (n · d) ≤ n · (x_i - p0).
    const double num = Dot(poly.n[i], Sub(poly.x[i], p0));
    const double den = Dot(poly.n[i], d);
    if (std::abs(den) <= PARALLEL_TOL * len)
    {
      // Parallel to this edge: either entirely inside its half-plane (or on the edge
      // itself), or entirely outside.
      if (num < -ON_EDGE_TOL * poly.h)
      {
        return false;
      }
      continue;
    }
    const double t = num / den;
    if (den > 0.0)
    {
      t_hi = std::min(t_hi, t);
    }
    else
    {
      t_lo = std::max(t_lo, t);
    }
    if (t_hi - t_lo <= T_TOL)
    {
      return false;
    }
  }
  return (t_hi - t_lo > T_TOL);
}

}  // namespace

LineIntegrator::LineIntegrator(const mfem::ParMesh &mesh_,
                               const std::vector<mfem::Vector> &path, bool closed,
                               int order, int min_points_per_segment)
  : comm(mesh_.GetComm()), sdim(mesh_.SpaceDimension()), mesh_ptr(&mesh_)
{
  // Non-const access is required for element transformations only.
  auto &mesh = const_cast<mfem::ParMesh &>(mesh_);
  MFEM_VERIFY(mesh.Dimension() == 2 && (sdim == 2 || sdim == 3),
              "LineIntegrator requires a 2D mesh or a 2D surface mesh in 3D!");
  MFEM_VERIFY(path.size() >= 2, "LineIntegrator requires a path with at least 2 points!");
  for (const auto &p : path)
  {
    MFEM_VERIFY(p.Size() == sdim,
                "LineIntegrator path points must match the mesh space dimension!");
  }
  const auto &ir = mfem::IntRules.Get(mfem::Geometry::SEGMENT, std::max(order, 1));
  const int nq = ir.GetNPoints();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);

  // Precompute the element polygons once for all segments.
  const int ne = mesh.GetNE();
  std::vector<ElementPolygon> polys;
  polys.reserve(ne);
  for (int e = 0; e < ne; e++)
  {
    polys.push_back(GetElementPolygon(mesh, e));
  }

  mfem::InverseElementTransformation inv;
  inv.SetInitialGuessType(mfem::InverseElementTransformation::ClosestPhysNode);
  inv.SetSolverType(mfem::InverseElementTransformation::NewtonElementProject);
  mfem::Vector x(sdim);

  const int num_segments =
      static_cast<int>(path.size()) - ((closed && path.size() > 2) ? 0 : 1);
  for (int s = 0; s < num_segments; s++)
  {
    const auto &q0 = path[s];
    const auto &q1 = path[(s + 1) % path.size()];
    const Vec3 p0 = {q0(0), q0(1), (sdim == 3) ? q0(2) : 0.0};
    const Vec3 p1 = {q1(0), q1(1), (sdim == 3) ? q1(2) : 0.0};
    const Vec3 d = Sub(p1, p0);
    const double len = Norm(d);
    MFEM_VERIFY(len > 0.0, "LineIntegrator path segments must have nonzero length!");
    path_length += len;
    Vec3 bb_min, bb_max;
    for (int k = 0; k < 3; k++)
    {
      bb_min[k] = std::min(p0[k], p1[k]);
      bb_max[k] = std::max(p0[k], p1[k]);
    }

    // Find the parameter interval of the segment within each local element.
    struct Interval
    {
      double a, b;
      int e;
    };
    std::vector<Interval> intervals;
    std::vector<double> breakpoints = {0.0, 1.0};
    for (int e = 0; e < ne; e++)
    {
      const auto &poly = polys[e];
      const double tol = ON_EDGE_TOL * poly.h + PLANE_TOL * poly.h;
      bool overlap = true;
      for (int k = 0; k < 3; k++)
      {
        overlap = overlap && (poly.bb_min[k] <= bb_max[k] + tol) &&
                  (poly.bb_max[k] >= bb_min[k] - tol);
      }
      double a, b;
      if (overlap && ClipSegment(poly, p0, d, a, b))
      {
        intervals.push_back({a, b, e});
        breakpoints.push_back(a);
        breakpoints.push_back(b);
      }
    }
    if (min_points_per_segment > 0)
    {
      const int n_split = (min_points_per_segment + nq - 1) / nq;
      for (int k = 1; k < n_split; k++)
      {
        breakpoints.push_back(static_cast<double>(k) / n_split);
      }
    }

    // Gather the breakpoints from all processes, so every process sees the same global
    // partition of the segment into sub-intervals.
    std::vector<double> all_breakpoints;
    {
      int count = static_cast<int>(breakpoints.size());
      std::vector<int> counts(size), displs(size);
      Mpi::Allgather(1, &count, counts.data(), comm);
      int total = 0;
      for (int r = 0; r < size; r++)
      {
        displs[r] = total;
        total += counts[r];
      }
      all_breakpoints.resize(total);
      Mpi::Allgatherv(count, breakpoints.data(), all_breakpoints.data(), counts.data(),
                      displs.data(), comm);
    }
    std::sort(all_breakpoints.begin(), all_breakpoints.end());
    std::vector<double> t;
    t.reserve(all_breakpoints.size());
    for (auto v : all_breakpoints)
    {
      v = std::clamp(v, 0.0, 1.0);
      if (t.empty() || v - t.back() > T_TOL)
      {
        t.push_back(v);
      }
    }
    t.front() = 0.0;
    t.back() = 1.0;
    const int nsub = static_cast<int>(t.size()) - 1;
    num_subintervals += nsub;

    // Assign each sub-interval to the first local element containing it, then to the
    // lowest rank with such an element.
    std::vector<int> local_elem(nsub, -1), owner(nsub, size);
    for (const auto &[a, b, e] : intervals)
    {
      auto it = std::lower_bound(t.begin(), t.end(), a - 10.0 * T_TOL);
      for (auto j = static_cast<int>(it - t.begin());
           j < nsub && t[j + 1] <= b + 10.0 * T_TOL; j++)
      {
        if (local_elem[j] < 0)
        {
          local_elem[j] = e;
          owner[j] = rank;
        }
      }
    }
    Mpi::GlobalMin(nsub, owner.data(), comm);

    // Build the quadrature points for the sub-intervals owned by this process.
    for (int j = 0; j < nsub; j++)
    {
      const double dt = t[j + 1] - t[j];
      if (owner[j] == size)
      {
        uncovered_length += dt * len;
        continue;
      }
      if (owner[j] != rank)
      {
        continue;
      }
      const int e = local_elem[j];
      inv.SetTransformation(*mesh.GetElementTransformation(e));
      for (int i = 0; i < nq; i++)
      {
        const auto &ip = ir.IntPoint(i);
        const double tq = t[j] + ip.x * dt;
        for (int k = 0; k < sdim; k++)
        {
          x(k) = p0[k] + tq * d[k];
        }
        Sample sample;
        sample.element = e;
        const int res = inv.Transform(x, sample.point);
        MFEM_VERIFY(res != mfem::InverseElementTransformation::Unknown,
                    "LineIntegrator failed to locate a quadrature point in its element!");
        for (int k = 0; k < 3; k++)
        {
          sample.weighted_tangent[k] = ip.weight * dt * d[k];
        }
        samples.push_back(sample);
      }
    }
  }
}

double LineIntegrator::LocalIntegrate(const mfem::ParGridFunction &F) const
{
  MFEM_VERIFY(F.FESpace()->GetMesh() == mesh_ptr,
              "LineIntegrator field does not live on the integration mesh!");
  F.HostRead();
  const int ndot = std::min(F.VectorDim(), sdim);
  mfem::Vector val(F.VectorDim());
  double result = 0.0;
  for (const auto &sample : samples)
  {
    F.GetVectorValue(sample.element, sample.point, val);
    for (int k = 0; k < ndot; k++)
    {
      result += sample.weighted_tangent[k] * val(k);
    }
  }
  return result;
}

double LineIntegrator::Integrate(const mfem::ParGridFunction &F) const
{
  double result = LocalIntegrate(F);
  Mpi::GlobalSum(1, &result, comm);
  return result;
}

std::complex<double> LineIntegrator::Integrate(const GridFunction &F) const
{
  double result[2] = {LocalIntegrate(F.Real()),
                      F.HasImag() ? LocalIntegrate(F.Imag()) : 0.0};
  Mpi::GlobalSum(2, result, comm);
  return {result[0], result[1]};
}

}  // namespace palace::fem
