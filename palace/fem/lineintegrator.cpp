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

// Boundary of a (possibly curved) planar element, as a closed polyline of points on the
// true element boundary. Straight edges are a single piece, curved edges are split into
// several pieces, and the crossings of a path segment with curved pieces are refined with
// Newton's method on the exact edge parametrization.
struct ElementBoundary
{
  std::vector<Vec3> x;                       // Polyline points, x.back() == x.front()
  std::vector<int> edge;                     // Element edge of each piece
  std::vector<double> s0, s1;                // Edge parameter range of each piece
  std::array<mfem::IntegrationPoint, 4> rv;  // Reference vertices
  int nv = 0;
  bool curved = false;
  Vec3 normal;          // Unit element normal (z for 2D meshes)
  Vec3 e1, e2;          // Unit in-plane basis
  Vec3 bb_min, bb_max;  // Bounding box of the polyline points
  double h = 0.0;       // Maximum polyline extent
};

ElementBoundary GetElementBoundary(mfem::ParMesh &mesh, int e, int pieces_per_edge)
{
  const int sdim = mesh.SpaceDimension();
  const auto geom = mesh.GetElementBaseGeometry(e);
  MFEM_VERIFY(geom == mfem::Geometry::TRIANGLE || geom == mfem::Geometry::SQUARE,
              "LineIntegrator only supports triangle and quadrilateral elements!");
  ElementBoundary bdr;
  bdr.nv = mfem::Geometry::NumVerts[geom];
  bdr.curved = (pieces_per_edge > 1);
  const auto *ref_verts = mfem::Geometries.GetVertices(geom);
  for (int i = 0; i < bdr.nv; i++)
  {
    bdr.rv[i] = ref_verts->IntPoint(i);
  }
  auto &T = *mesh.GetElementTransformation(e);
  mfem::Vector xv(sdim);
  mfem::IntegrationPoint ip;
  for (int k = 0; k < bdr.nv; k++)
  {
    const auto &ra = bdr.rv[k], &rb = bdr.rv[(k + 1) % bdr.nv];
    for (int j = 0; j < pieces_per_edge; j++)
    {
      const double sv = static_cast<double>(j) / pieces_per_edge;
      ip.Set2(ra.x + sv * (rb.x - ra.x), ra.y + sv * (rb.y - ra.y));
      T.Transform(ip, xv);
      bdr.x.push_back({xv(0), xv(1), (sdim == 3) ? xv(2) : 0.0});
      bdr.edge.push_back(k);
      bdr.s0.push_back(sv);
      bdr.s1.push_back(static_cast<double>(j + 1) / pieces_per_edge);
    }
  }
  bdr.x.push_back(bdr.x.front());
  for (int d = 0; d < 3; d++)
  {
    bdr.bb_min[d] = bdr.bb_max[d] = bdr.x[0][d];
    for (const auto &x : bdr.x)
    {
      bdr.bb_min[d] = std::min(bdr.bb_min[d], x[d]);
      bdr.bb_max[d] = std::max(bdr.bb_max[d], x[d]);
    }
  }
  for (int d = 0; d < 3; d++)
  {
    bdr.h = std::max(bdr.h, bdr.bb_max[d] - bdr.bb_min[d]);
  }
  MFEM_VERIFY(bdr.h > 0.0, "Degenerate element in LineIntegrator!");

  // Element normal and in-plane basis, from the vertices.
  const auto &x0 = bdr.x[0];
  const auto &x1 = bdr.x[pieces_per_edge];
  const auto &x2 = bdr.x[2 * pieces_per_edge];
  if (sdim == 2)
  {
    bdr.normal = {0.0, 0.0, 1.0};
    bdr.e1 = {1.0, 0.0, 0.0};
    bdr.e2 = {0.0, 1.0, 0.0};
  }
  else
  {
    const Vec3 a = (bdr.nv == 4) ? Sub(x2, x0) : Sub(x1, x0);
    const Vec3 b = (bdr.nv == 4) ? Sub(bdr.x[3 * pieces_per_edge], x1) : Sub(x2, x0);
    bdr.normal = Cross(a, b);
    const double n = Norm(bdr.normal);
    MFEM_VERIFY(n > 0.0, "Degenerate element in LineIntegrator!");
    for (auto &c : bdr.normal)
    {
      c /= n;
    }
    bdr.e1 = Sub(x1, x0);
    const double n1 = Norm(bdr.e1);
    for (auto &c : bdr.e1)
    {
      c /= n1;
    }
    bdr.e2 = Cross(bdr.normal, bdr.e1);
  }
  return bdr;
}

// In-plane cross product.
inline double Cross2(const ElementBoundary &bdr, const Vec3 &a, const Vec3 &b)
{
  return Dot(bdr.normal, Cross(a, b));
}

// Refine the crossing of the segment x(t) = p0 + t d with a curved element edge using
// Newton's method on the exact edge parametrization X(s), starting from the crossing with
// the polyline piece. Returns false if Newton does not converge near the initial guess.
bool RefineCrossing(mfem::ElementTransformation &T, const ElementBoundary &bdr, int piece,
                    const Vec3 &p0, const Vec3 &d, double &s, double &t)
{
  const int sdim = T.GetSpaceDim();
  const auto &ra = bdr.rv[bdr.edge[piece]], &rb = bdr.rv[(bdr.edge[piece] + 1) % bdr.nv];
  const double dr[2] = {rb.x - ra.x, rb.y - ra.y};
  const double ds_piece = bdr.s1[piece] - bdr.s0[piece];
  mfem::IntegrationPoint ip;
  mfem::Vector xv(sdim);
  for (int it = 0; it < 20; it++)
  {
    ip.Set2(ra.x + s * dr[0], ra.y + s * dr[1]);
    T.SetIntPoint(&ip);
    T.Transform(ip, xv);
    const auto &J = T.Jacobian();
    Vec3 F = {0.0, 0.0, 0.0}, dX = {0.0, 0.0, 0.0};
    for (int k = 0; k < sdim; k++)
    {
      F[k] = xv(k) - p0[k] - t * d[k];
      dX[k] = J(k, 0) * dr[0] + J(k, 1) * dr[1];
    }
    // Solve [dX, -d] [δs, δt]^T = -F in the element plane.
    const double a11 = Dot(dX, bdr.e1), a12 = -Dot(d, bdr.e1);
    const double a21 = Dot(dX, bdr.e2), a22 = -Dot(d, bdr.e2);
    const double det = a11 * a22 - a12 * a21;
    if (std::abs(det) <= std::numeric_limits<double>::epsilon() * Norm(dX) * Norm(d))
    {
      return false;
    }
    const double b1 = -Dot(F, bdr.e1), b2 = -Dot(F, bdr.e2);
    const double delta_s = (b1 * a22 - a12 * b2) / det;
    const double delta_t = (a11 * b2 - a21 * b1) / det;
    s += delta_s;
    t += delta_t;
    if (s < bdr.s0[piece] - ds_piece || s > bdr.s1[piece] + ds_piece)
    {
      return false;
    }
    if (std::abs(delta_s) < 1.0e-14 && std::abs(delta_t) < 1.0e-14)
    {
      return true;
    }
  }
  return false;
}

// Whether the physical point x lies in the element (on its boundary within a tolerance),
// using the exact (possibly curved) element geometry.
bool InElement(mfem::InverseElementTransformation &inv, mfem::ElementTransformation &T,
               const ElementBoundary &bdr, const mfem::Vector &x)
{
  mfem::IntegrationPoint ip;
  const int res = inv.Transform(x, ip);
  if (res == mfem::InverseElementTransformation::Inside)
  {
    return true;
  }
  if (res == mfem::InverseElementTransformation::Unknown)
  {
    return false;
  }
  mfem::Vector y(x.Size());
  T.Transform(ip, y);
  return y.DistanceTo(x) <= ON_EDGE_TOL * bdr.h;
}

// Compute the parameter intervals of the segment x(t) = p0 + t d, t ∈ [0, 1], lying in the
// element (curved elements are not necessarily convex, so there may be several). Segments
// lying on an element edge (within a small tolerance) are considered to be in the element,
// so both neighbors of the edge report the interval and the caller resolves ownership.
void IntersectSegment(mfem::ElementTransformation &T,
                      mfem::InverseElementTransformation &inv, const ElementBoundary &bdr,
                      const Vec3 &p0, const Vec3 &d,
                      std::vector<std::pair<double, double>> &out)
{
  out.clear();
  const int sdim = T.GetSpaceDim();
  const double len = Norm(d);

  // The segment must lie in the plane of the element (only relevant for surface meshes).
  {
    const double d0 = Dot(bdr.normal, Sub(p0, bdr.x[0]));
    const double d1 = d0 + Dot(bdr.normal, d);
    if (std::abs(d0) > PLANE_TOL * bdr.h || std::abs(d1) > PLANE_TOL * bdr.h)
    {
      return;
    }
  }

  // Quick rejection: the element lies entirely on one side of the segment's line. Curved
  // edges may bulge slightly beyond the polyline points between them.
  {
    const double margin = (bdr.curved ? 0.05 : ON_EDGE_TOL) * bdr.h * len;
    bool left = true, right = true;
    for (const auto &x : bdr.x)
    {
      const double c = Cross2(bdr, d, Sub(x, p0));
      left = left && (c > margin);
      right = right && (c < -margin);
    }
    if (left || right)
    {
      return;
    }
  }

  // Crossings of the segment with the element boundary.
  std::vector<double> t = {0.0, 1.0};
  const int npieces = static_cast<int>(bdr.x.size()) - 1;
  for (int j = 0; j < npieces; j++)
  {
    const Vec3 &A = bdr.x[j];
    const Vec3 w = Sub(bdr.x[j + 1], A);
    const double wlen = Norm(w);
    const Vec3 r = Sub(A, p0);
    const double den = Cross2(bdr, d, w);
    if (std::abs(den) <= PARALLEL_TOL * len * wlen)
    {
      // Parallel: if collinear, the overlap endpoints are crossings.
      if (std::abs(Cross2(bdr, r, d)) <= ON_EDGE_TOL * bdr.h * len)
      {
        for (const auto &P : {A, bdr.x[j + 1]})
        {
          const double tp = Dot(Sub(P, p0), d) / (len * len);
          if (tp > 0.0 && tp < 1.0)
          {
            t.push_back(tp);
          }
        }
      }
      continue;
    }
    double tc = Cross2(bdr, r, w) / den;
    const double u = Cross2(bdr, r, d) / den;
    if (bdr.curved)
    {
      // Refine on the exact edge, accepting crossings slightly outside the piece which
      // the polyline may miss.
      if (u < -0.5 || u > 1.5)
      {
        continue;
      }
      double sc = bdr.s0[j] + u * (bdr.s1[j] - bdr.s0[j]);
      double tn = tc;
      if (RefineCrossing(T, bdr, j, p0, d, sc, tn) && sc >= 0.0 && sc <= 1.0)
      {
        tc = tn;
      }
      else if (u < 0.0 || u > 1.0)
      {
        continue;
      }
    }
    else if (u < -ON_EDGE_TOL || u > 1.0 + ON_EDGE_TOL)
    {
      continue;
    }
    if (tc > 0.0 && tc < 1.0)
    {
      t.push_back(tc);
    }
  }
  std::sort(t.begin(), t.end());

  // Classify the pieces between consecutive crossings by testing their midpoints with the
  // exact geometry, and merge consecutive pieces inside the element.
  mfem::Vector xm(sdim);
  double a = -1.0, b = -1.0;
  for (std::size_t i = 0; i + 1 < t.size(); i++)
  {
    if (t[i + 1] - t[i] <= T_TOL)
    {
      continue;
    }
    const double tm = 0.5 * (t[i] + t[i + 1]);
    bool in_bb = true;
    for (int k = 0; k < sdim; k++)
    {
      xm(k) = p0[k] + tm * d[k];
      in_bb = in_bb && (xm(k) >= bdr.bb_min[k] - 0.05 * bdr.h) &&
              (xm(k) <= bdr.bb_max[k] + 0.05 * bdr.h);
    }
    if (in_bb && InElement(inv, T, bdr, xm))
    {
      if (a < 0.0)
      {
        a = t[i];
      }
      b = t[i + 1];
    }
    else if (a >= 0.0)
    {
      out.emplace_back(a, b);
      a = b = -1.0;
    }
  }
  if (a >= 0.0)
  {
    out.emplace_back(a, b);
  }
}

}  // namespace

LineIntegrator::LineIntegrator(const mfem::ParMesh &mesh_,
                               const std::vector<mfem::Vector> &path, bool closed,
                               int field_order)
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
  const int mesh_order =
      mesh.GetNodes() ? mesh.GetNodes()->FESpace()->GetMaxElementOrder() : 1;
  const auto &ir = mfem::IntRules.Get(mfem::Geometry::SEGMENT,
                                      2 * (std::max(field_order, 0) + mesh_order));
  const int nq = ir.GetNPoints();
  const int rank = Mpi::Rank(comm), size = Mpi::Size(comm);

  // Precompute the element boundaries once for all segments. Curved (high-order) mesh
  // elements use several polyline pieces per edge.
  const int ne = mesh.GetNE();
  const int pieces_per_edge = (mesh_order > 1) ? 2 * mesh_order : 1;
  std::vector<ElementBoundary> bdrs;
  bdrs.reserve(ne);
  for (int e = 0; e < ne; e++)
  {
    bdrs.push_back(GetElementBoundary(mesh, e, pieces_per_edge));
  }

  mfem::InverseElementTransformation inv;
  inv.SetInitialGuessType(mfem::InverseElementTransformation::ClosestPhysNode);
  inv.SetSolverType(mfem::InverseElementTransformation::NewtonElementProject);
  mfem::Vector x(sdim);
  std::vector<std::pair<double, double>> element_intervals;

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
      const auto &bdr = bdrs[e];
      const double tol = 0.05 * bdr.h;
      bool overlap = true;
      for (int k = 0; k < 3; k++)
      {
        overlap = overlap && (bdr.bb_min[k] <= bb_max[k] + tol) &&
                  (bdr.bb_max[k] >= bb_min[k] - tol);
      }
      if (!overlap)
      {
        continue;
      }
      auto &T = *mesh.GetElementTransformation(e);
      inv.SetTransformation(T);
      IntersectSegment(T, inv, bdr, p0, d, element_intervals);
      for (const auto &[a, b] : element_intervals)
      {
        intervals.push_back({a, b, e});
        breakpoints.push_back(a);
        breakpoints.push_back(b);
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
