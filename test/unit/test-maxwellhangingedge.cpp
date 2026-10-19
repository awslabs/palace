// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

using Catch::Matchers::WithinAbs;

namespace
{

using Point = std::array<double, 3>;
using Interval = std::pair<double, double>;

Point Combine(const Point &origin, double a, const Point &u, double b, const Point &n,
              double c, const Point &w)
{
  Point p;
  for (int d = 0; d < 3; d++)
  {
    p[d] = origin[d] + a * u[d] + b * n[d] + c * w[d];
  }
  return p;
}

// The fan of tetrahedra around one sheet-plane mesh edge E = [p1, q] of the C4 iteration-11
// mesh of record (decision 862, OI-LOC-HANG), in the device's own coordinates (the trace
// corner p1 of contour line 87, the edge direction -u, the in-plane normal n towards the
// gap, the sheet normal w; mesh units, Lc = 1.055 mm, so ~6.7 against edges of 1e-4 and
// elements of 2e-5 .. 1e-4). On the trace side of E (-n) the cracked PEC sheet: the
// coarse tets T_up / T_lo with the sheet faces (p1, q, s) at z0 +- delta (the two displaced
// copies). On the gap side (+n) the UNCRACKED gap plane z0, refined one level more: the
// fine tets through the hanging vertex m = (p1 + q) / 2 (a nonconforming fan: m is a slave
// vertex of the coarse E-faces (p1, q, a) of T_up / T_lo, whose fine neighbours W1 / W2
// tile them). Like every hanging vertex, m lies on E only to its rounding: here it is
// placed epsilon_m = 3e-15 (2 ULP of the coordinates) towards the gap, so the fine
// E-faces (p1, m, .) / (m, q, .) are tilted by epsilon_m / (l / 2) = 6e-11 rad against the
// coarse ones through E. The contour line is 1e-2 long, ends at p1 offset 1e-16 towards the
// trace side and deviates from E by theta = 5e-11 rad towards the gap (the identification's
// rounding class: line 87 measured 7e-11 rad, 2-8e-15 off the chain vertices), so over the
// edge it runs through the sliver between the coarse acceptance planes (through E) and the
// fine ones (through (p1, m)): no element's exact clip holds it there.
struct HangingFan
{
  static constexpr double kEdgeLength = 1.0756e-4;   // l: the coarse edge [p1, q]
  static constexpr double kLineLength = 1.0e-2;      // L: the contour line (L / l = 93)
  static constexpr double kLineAngle = 5.0e-11;      // theta: line vs E (towards the gap)
  static constexpr double kEndOffset = 1.0e-16;      // the line's end off p1 (trace side)
  static constexpr double kHangingOffset = 3.0e-15;  // epsilon_m: m off E (gap side)
  static constexpr double kCrackDisplacement = 3.2e-17;           // delta: the sheet copies
  static constexpr double kCoordinateScale = 6.7500948047023135;  // the C4 run's

  Point p1 = {-6.6832458113646211, 3.5962729328964325, 4.5506257110352654e-03};
  Point u{}, n{}, w = {0.0, 0.0, 1.0};
  std::vector<Point> vertices;
  std::vector<std::array<int, 4>> tets;
  std::vector<int> upper;  // 1 for a tet above the plane, 0 below
  Point line_begin{}, line_end{};

  HangingFan()
  {
    const double ux = -6.6812569713163139 - p1[0], uy = 3.5949149113331251 - p1[1];
    const double norm = std::sqrt(ux * ux + uy * uy);
    u = {ux / norm, uy / norm, 0.0};  // from p1 towards the far end of line 87
    n = {u[1], -u[0], 0.0};           // the in-plane normal: + = the gap side
    const double l = kEdgeLength;
    auto Add = [&](const Point &p)
    {
      vertices.push_back(p);
      return static_cast<int>(vertices.size()) - 1;
    };
    const int P1 = Add(p1);
    const int Q = Add(Combine(p1, l, u, 0.0, n, 0.0, w));
    const int M = Add(Combine(p1, 0.5 * l, u, kHangingOffset, n, 0.0, w));
    // The trace side: the sheet copies and the apexes of the coarse tets.
    const int S_UP = Add(Combine(p1, 0.5 * l, u, -0.8 * l, n, kCrackDisplacement, w));
    const int S_LO = Add(Combine(p1, 0.5 * l, u, -0.8 * l, n, -kCrackDisplacement, w));
    const int A_T = Add(Combine(p1, 0.5 * l, u, -0.3 * l, n, 1.0 * l, w));
    const int B_T = Add(Combine(p1, 0.5 * l, u, -0.3 * l, n, -1.0 * l, w));
    // The gap side, one level finer: the gap-plane vertices and the apexes.
    const int G1 = Add(Combine(p1, 0.25 * l, u, 0.4 * l, n, 0.0, w));
    const int G2 = Add(Combine(p1, 0.75 * l, u, 0.4 * l, n, 0.0, w));
    const int C1_UP = Add(Combine(p1, 0.25 * l, u, 0.15 * l, n, 0.5 * l, w));
    const int C2_UP = Add(Combine(p1, 0.75 * l, u, 0.15 * l, n, 0.5 * l, w));
    const int C1_LO = Add(Combine(p1, 0.25 * l, u, 0.15 * l, n, -0.5 * l, w));
    const int C2_LO = Add(Combine(p1, 0.75 * l, u, 0.15 * l, n, -0.5 * l, w));
    auto Tet = [&](int a, int b, int c, int d, bool above)
    {
      tets.push_back({a, b, c, d});
      upper.push_back(above ? 1 : 0);
    };
    // Above the plane: the coarse trace tet, the two fine wedges on its E-face, the two
    // fine gap tets.
    Tet(P1, Q, S_UP, A_T, true);
    Tet(P1, M, A_T, C1_UP, true);
    Tet(M, Q, A_T, C2_UP, true);
    Tet(P1, M, C1_UP, G1, true);
    Tet(M, Q, C2_UP, G2, true);
    // Below the plane: the same fan with the lower sheet copy.
    Tet(P1, Q, S_LO, B_T, false);
    Tet(P1, M, B_T, C1_LO, false);
    Tet(M, Q, B_T, C2_LO, false);
    Tet(P1, M, C1_LO, G1, false);
    Tet(M, Q, C2_LO, G2, false);
    // The contour line: from the far end (on the gap side by theta L) to p1 (on the trace
    // side by the end offset), in the sheet plane.
    line_end = Combine(p1, 0.0, u, -kEndOffset, n, 0.0, w);
    line_begin = Combine(p1, kLineLength, u, kLineAngle * kLineLength, n, 0.0, w);
  }

  // The parameter at which the line passes the far end of the edge (q).
  double EdgeBegin() const { return 1.0 - kEdgeLength / kLineLength; }

  mfem::Mesh Serial() const
  {
    mfem::Mesh mesh(3, static_cast<int>(vertices.size()), static_cast<int>(tets.size()), 0,
                    3);
    for (const auto &p : vertices)
    {
      mesh.AddVertex(p.data());
    }
    for (std::size_t t = 0; t < tets.size(); t++)
    {
      mesh.AddTet(tets[t].data(), 1 + upper[t]);
    }
    mesh.FinalizeTetMesh(1, 0, true);
    return mesh;
  }

  // The fan distributed with the plane as the partition: the upper tets on rank 0, the
  // lower ones on the last rank (one rank: everything on it).
  std::unique_ptr<mfem::ParMesh> Parallel(MPI_Comm comm) const
  {
    auto serial = Serial();
    mfem::Array<int> partitioning(serial.GetNE());
    for (int e = 0; e < serial.GetNE(); e++)
    {
      partitioning[e] = serial.GetAttribute(e) == 2 ? 0 : Mpi::Size(comm) - 1;
    }
    return std::make_unique<mfem::ParMesh>(comm, serial, partitioning.GetData());
  }
};

// The intervals of every rank, gathered.
std::vector<Interval> Gather(const std::vector<Interval> &local, MPI_Comm comm)
{
  const int size = Mpi::Size(comm);
  std::vector<int> counts(size), offsets(size);
  const int count = 2 * static_cast<int>(local.size());
  MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
  int total = 0;
  for (int rank = 0; rank < size; rank++)
  {
    offsets[rank] = total;
    total += counts[rank];
  }
  std::vector<double> send(count), receive(total);
  for (std::size_t i = 0; i < local.size(); i++)
  {
    send[2 * i] = local[i].first;
    send[2 * i + 1] = local[i].second;
  }
  MPI_Allgatherv(send.data(), count, MPI_DOUBLE, receive.data(), counts.data(),
                 offsets.data(), MPI_DOUBLE, comm);
  std::vector<Interval> gathered(total / 2);
  for (std::size_t i = 0; i < gathered.size(); i++)
  {
    gathered[i] = {receive[2 * i], receive[2 * i + 1]};
  }
  return gathered;
}

}  // namespace

TEST_CASE("SurfaceResponseOperatorMaxwellLineAlongNonconformingEdge",
          "[surfaceresponseoperator][maxwell][Serial][Parallel]")
{
  // Decision 862 (OI-LOC-HANG): a contour line following a sheet-plane mesh edge of a
  // nonconforming fan must be covered over the whole edge by the in-face intervals. The
  // coarse tets' E-faces (through E) and the fine tets' E-faces (through the hanging
  // vertex m, 2 ULP off E) are grazing faces of the line: their exact clips hold the line
  // on disjoint sides (the sliver between the two planes holds it for neither), and a
  // coplanarity test at the line's far ends (|lambda| = theta L / h = 5e-9 > 1e-9 there)
  // does not skip them, so their ill-conditioned cuts leave the edge uncovered (the C4
  // iteration-11 abort). Tested over the element's own extent (|lambda| <= theta l / h =
  // 5e-11) the grazing faces are coplanar and the in-face interval of every tet of the fan
  // is the clip by its transversal faces: the edge's extent.
  const auto comm = Mpi::World();
  HangingFan fan;
  auto mesh = fan.Parallel(comm);
  REQUIRE(mesh->GetGlobalNE() == 10);
  const auto local = FindMaxwellSegmentIntersections(*mesh, 3, fan.line_begin, fan.line_end,
                                                     HangingFan::kCoordinateScale);
  const auto exact = Gather(local.exact, comm);
  const auto in_face = Gather(local.in_face, comm);
  const double edge_begin = fan.EdgeBegin();
  const double line_length = HangingFan::kLineLength;

  // The far 99 % of the line crosses no element of the fan: taken as covered, so the
  // resolution reads the fan's edge alone.
  auto with_far_part = [&](std::vector<Interval> intervals)
  {
    intervals.emplace_back(0.0, edge_begin);
    return intervals;
  };
  const auto coverage = ResolveMaxwellLineCoverage(
      with_far_part(exact), in_face, line_length, HangingFan::kCoordinateScale, 87,
      fan.line_begin, fan.line_end);
  // The edge is covered: the union of the resolved intervals reaches from the edge's far
  // end to p1 with no piece missing (every breakpoint is a transversal face's).
  std::vector<Interval> resolved(coverage.intervals);
  std::sort(resolved.begin(), resolved.end());
  double reach = 0.0;
  for (const auto &[begin, end] : resolved)
  {
    REQUIRE(begin <= reach + 1.0e-10);
    reach = std::max(reach, end);
  }
  CHECK_THAT(reach, WithinAbs(1.0, 1.0e-12));
  // The part of the edge the exact clips reject (the crack on the trace side, the sliver
  // between the coarse and the fine acceptance planes) is integrated from the in-face
  // intervals: one stretch, shorter than the edge (the fine gap tets hold the line exactly
  // where it is on their side of the sliver).
  CHECK(coverage.in_face_count == 1);
  CHECK(coverage.in_face_length > 0.1 * HangingFan::kEdgeLength);
  CHECK(coverage.in_face_length < HangingFan::kEdgeLength);
  CHECK(coverage.closed_crack_count == 0);
  // Every tet of the fan contributes its whole footprint along the edge as an in-face
  // interval: the coarse ones [edge_begin, 1], the fine ones the half edges.
  int whole_edge = 0, first_half = 0, second_half = 0;
  for (const auto &[begin, end] : in_face)
  {
    const bool at_p1 = std::abs(end - 1.0) <= 1.0e-10;
    const bool from_q = std::abs(begin - edge_begin) <= 1.0e-10;
    const bool from_m = std::abs(begin - 0.5 * (edge_begin + 1.0)) <= 1.0e-10;
    const bool to_m = std::abs(end - 0.5 * (edge_begin + 1.0)) <= 1.0e-10;
    whole_edge += (from_q && at_p1) ? 1 : 0;
    first_half += (from_m && at_p1) ? 1 : 0;
    second_half += (from_q && to_m) ? 1 : 0;
  }
  CHECK(whole_edge == 2);   // T_up, T_lo
  CHECK(first_half == 4);   // W1 / G1 above and below
  CHECK(second_half == 4);  // W2 / G2 above and below
  CHECK(in_face.size() == 10);
}

}  // namespace palace
