
#include "geodata.hpp"

#include "geodata_impl.hpp"

#include <array>
#include <limits>
#include <random>
#include <unordered_set>
#include <Eigen/Dense>
#include "utils/communication.hpp"
#include "utils/omp.hpp"

namespace palace::mesh
{

using Vector3dMap = Eigen::Map<Eigen::Vector3d>;
using CVector3dMap = Eigen::Map<const Eigen::Vector3d>;

// Compute a lexicographic comparison of Eigen Vector3d.
bool EigenLE(const Eigen::Vector3d &x, const Eigen::Vector3d &y)
{
  return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

// Helper for collecting a point cloud from a mesh, used in calculating bounding boxes and
// bounding balls. Returns the dominant rank, for which the vertices argument will be
// filled, while all other ranks will have an empty vector. Vertices are de-duplicated to a
// certain floating point precision.
int CollectPointCloudOnRoot(const mfem::ParMesh &mesh, const mfem::Array<int> &marker,
                            bool bdr, std::vector<Eigen::Vector3d> &vertices)
{
  if (!mesh.GetNodes())
  {
    // Linear mesh, work with element vertices directly.
    PalacePragmaOmp(parallel)
    {
      std::unordered_set<int> vertex_indices;
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetBdrElement(i)->GetVertices();
          vertex_indices.insert(verts, verts + mesh.GetBdrElement(i)->GetNVertices());
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          const int *verts = mesh.GetElement(i)->GetVertices();
          vertex_indices.insert(verts, verts + mesh.GetElement(i)->GetNVertices());
        }
      }
      PalacePragmaOmp(critical(PointCloud))
      {
        for (auto i : vertex_indices)
        {
          const auto &vx = mesh.GetVertex(i);
          vertices.emplace_back(vx[0], vx[1], vx[2]);
        }
      }
    }
  }
  else
  {
    // Curved mesh, need to process point matrices.
    const int ref = mesh.GetNodes()->FESpace()->GetMaxElementOrder();
    auto AddPoints = [&](mfem::GeometryRefiner &refiner, mfem::ElementTransformation &T,
                         mfem::DenseMatrix &pointmat,
                         std::vector<Eigen::Vector3d> &loc_vertices)
    {
      mfem::RefinedGeometry *RefG = refiner.Refine(T.GetGeometryType(), ref);
      T.Transform(RefG->RefPts, pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        loc_vertices.emplace_back(pointmat(0, j), pointmat(1, j), pointmat(2, j));
      }
    };
    PalacePragmaOmp(parallel)
    {
      mfem::GeometryRefiner refiner;
      mfem::IsoparametricTransformation T;
      mfem::DenseMatrix pointmat;  // 3 x N
      std::vector<Eigen::Vector3d> loc_vertices;
      if (bdr)
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNBE(); i++)
        {
          if (!marker[mesh.GetBdrAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetBdrElementTransformation(i, &T);
          AddPoints(refiner, T, pointmat, loc_vertices);
        }
      }
      else
      {
        PalacePragmaOmp(for schedule(static))
        for (int i = 0; i < mesh.GetNE(); i++)
        {
          if (!marker[mesh.GetAttribute(i) - 1])
          {
            continue;
          }
          mesh.GetElementTransformation(i, &T);
          AddPoints(refiner, T, pointmat, loc_vertices);
        }
      }
      PalacePragmaOmp(critical(PointCloud))
      {
        for (const auto &v : loc_vertices)
        {
          vertices.push_back(v);
        }
      }
    }
  }

  // dominant_rank will perform the calculation.
  MPI_Comm comm = mesh.GetComm();
  const auto num_vertices = int(vertices.size());
  const int dominant_rank = [&]()
  {
    int vert = num_vertices, rank = Mpi::Rank(comm);
    Mpi::GlobalMaxLoc(1, &vert, &rank, comm);
    return rank;
  }();
  std::vector<int> recv_counts(Mpi::Size(comm)), displacements;
  std::vector<Eigen::Vector3d> collected_vertices;
  MPI_Gather(&num_vertices, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, dominant_rank,
             comm);
  if (dominant_rank == Mpi::Rank(comm))
  {
    // First displacement is zero, then after is the partial sum of recv_counts.
    displacements.resize(Mpi::Size(comm));
    displacements[0] = 0;
    std::partial_sum(recv_counts.begin(), recv_counts.end() - 1, displacements.begin() + 1);

    // Add on slots at the end of vertices for the incoming data.
    collected_vertices.resize(std::accumulate(recv_counts.begin(), recv_counts.end(), 0));

    // MPI transfer will be done with MPI_DOUBLE, so duplicate all these values.
    for (auto &x : displacements)
    {
      x *= 3;
    }
    for (auto &x : recv_counts)
    {
      x *= 3;
    }
  }

  // Gather the data to the dominant rank.
  static_assert(sizeof(Eigen::Vector3d) == 3 * sizeof(double));
  MPI_Gatherv(vertices.data(), 3 * num_vertices, MPI_DOUBLE, collected_vertices.data(),
              recv_counts.data(), displacements.data(), MPI_DOUBLE, dominant_rank, comm);

  // Deduplicate vertices. Given floating point precision, need a tolerance.
  if (dominant_rank == Mpi::Rank(comm))
  {
    auto vertex_equality = [](const auto &x, const auto &y)
    {
      constexpr double tolerance = 10.0 * std::numeric_limits<double>::epsilon();
      return std::abs(x[0] - y[0]) < tolerance && std::abs(x[1] - y[1]) < tolerance &&
             std::abs(x[2] - y[2]) < tolerance;
    };
    vertices = std::move(collected_vertices);
    std::sort(vertices.begin(), vertices.end(), EigenLE);
    vertices.erase(std::unique(vertices.begin(), vertices.end(), vertex_equality),
                   vertices.end());
  }
  else
  {
    vertices.clear();
  }

  return dominant_rank;
}

// Compute the distance from a point orthogonal to the list of normal axes, relative to
// the given origin.
auto PerpendicularDistance(const std::initializer_list<Eigen::Vector3d> &normals,
                           const Eigen::Vector3d &origin, const Eigen::Vector3d &v)
{
  Eigen::Vector3d v0 = v - origin;
  for (const auto &n : normals)
  {
    v0 -= n.dot(v0) * n;
  }
  return v0.norm();
}

// Calculates a bounding box from a point cloud, result is broadcast across all processes.
BoundingBox BoundingBoxFromPointCloud(MPI_Comm comm,
                                      const std::vector<Eigen::Vector3d> &vertices,
                                      int dominant_rank)
{
  BoundingBox box;
  if (dominant_rank == Mpi::Rank(comm))
  {
    // Pick a candidate 000 vertex using lexicographic sort. This can be vulnerable to
    // floating point precision if the box is axis aligned, but not floating point exact.
    // Pick candidate 111 as the furthest from this candidate, then reassign 000 as the
    // furthest from 111. Such a pair has to form the diagonal for a point cloud defining a
    // box. Verify that p_111 is also the maximum distance from p_000 -> a diagonal is
    // found.
    MFEM_VERIFY(vertices.size() >= 4,
                "A bounding box requires a minimum of four vertices for this algorithm!");
    auto p_000 = std::min_element(vertices.begin(), vertices.end(), EigenLE);
    auto p_111 =
        std::max_element(vertices.begin(), vertices.end(),
                         [p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                         { return (x - *p_000).norm() < (y - *p_000).norm(); });
    p_000 = std::max_element(vertices.begin(), vertices.end(),
                             [p_111](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                             { return (x - *p_111).norm() < (y - *p_111).norm(); });
    MFEM_ASSERT(std::max_element(vertices.begin(), vertices.end(),
                                 [p_000](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                 { return (x - *p_000).norm() < (y - *p_000).norm(); }) ==
                    p_111,
                "p_000 and p_111 must be mutually opposing points!");

    // Define a diagonal of the ASSUMED cuboid bounding box.
    const auto &v_000 = *p_000;
    const auto &v_111 = *p_111;
    MFEM_VERIFY(&v_000 != &v_111, "Minimum and maximum extents cannot be identical!");
    const Eigen::Vector3d origin = v_000;
    const Eigen::Vector3d n_1 = (v_111 - v_000).normalized();

    // Find the vertex furthest from the diagonal axis. We cannot know yet if this defines
    // (001) or (011).
    const auto &t_0 = *std::max_element(vertices.begin(), vertices.end(),
                                        [&](const auto &x, const auto &y)
                                        {
                                          return PerpendicularDistance({n_1}, origin, x) <
                                                 PerpendicularDistance({n_1}, origin, y);
                                        });
    MFEM_VERIFY(&t_0 != &v_000 && &t_0 != &v_111, "Vertices are degenerate!");

    // Use the discovered vertex to define a second direction and thus a plane. n_1 and n_2
    // now define a planar coordinate system intersecting the main diagonal, and two
    // opposite edges of the cuboid.
    const Eigen::Vector3d n_2 =
        ((t_0 - origin) - ((t_0 - origin).dot(n_1) * n_1)).normalized();

    // Collect the furthest point from the plane to determine if the box is planar. Look for
    // a component that maximizes distance from the planar system: complete the axes with a
    // cross, then use a dot product to pick the greatest deviation.
    constexpr double rel_tol = 1.0e-6;
    auto max_distance = PerpendicularDistance(
        {n_1, n_2}, origin,
        *std::max_element(vertices.begin(), vertices.end(),
                          [&](const auto &x, const auto &y)
                          {
                            return PerpendicularDistance({n_1, n_2}, origin, x) <
                                   PerpendicularDistance({n_1, n_2}, origin, y);
                          }));
    box.planar = (max_distance < rel_tol * (v_111 - v_000).norm());

    // For the non-planar case, collect points furthest from the plane and choose the one
    // closest to the origin as the next vertex which might be (001) or (011).
    const auto &t_1 = [&]()
    {
      if (box.planar)
      {
        return t_0;
      }
      std::vector<Eigen::Vector3d> vertices_out_of_plane;
      std::copy_if(vertices.begin(), vertices.end(),
                   std::back_inserter(vertices_out_of_plane),
                   [&](const auto &v)
                   {
                     return std::abs(PerpendicularDistance({n_1, n_2}, origin, v) -
                                     max_distance) < rel_tol * max_distance;
                   });
      return *std::min_element(vertices_out_of_plane.begin(), vertices_out_of_plane.end(),
                               [&](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                               { return (x - origin).norm() < (y - origin).norm(); });
    }();

    // Given candidates t_0 and t_1, the closer to origin defines v_001.
    const bool t_0_gt_t_1 = (t_0 - origin).norm() > (t_1 - origin).norm();
    const auto &v_001 = t_0_gt_t_1 ? t_1 : t_0;
    const auto &v_011 = box.planar ? v_111 : (t_0_gt_t_1 ? t_0 : t_1);

    // Compute the center as halfway along the main diagonal.
    Vector3dMap(box.center.data()) = 0.5 * (v_000 + v_111);

    if constexpr (false)
    {
      fmt::print("box.center {}!\n", box.center);
      fmt::print("v_000 {}!\n", v_000);
      fmt::print("v_001 {}!\n", v_001);
      fmt::print("v_011 {}!\n", v_011);
      fmt::print("v_111 {}!\n", v_111);
    }

    // Compute the box axes. Using the 4 extremal points, we find the first two axes as the
    // edges which are closest to perpendicular. For a perfect rectangular prism point
    // cloud, we could instead compute the axes and length in each direction using the
    // found edges of the cuboid, but this does not work for non-rectangular prism
    // cross-sections or pyramid shapes.
    {
      const auto [e_0, e_1] = [&v_000, &v_001, &v_011, &v_111]()
      {
        std::array<const Eigen::Vector3d *, 4> verts = {&v_000, &v_001, &v_011, &v_111};
        Eigen::Vector3d e_0 = Eigen::Vector3d::Zero(), e_1 = Eigen::Vector3d::Zero();
        double dot_min = mfem::infinity();
        for (int i_0 = 0; i_0 < 4; i_0++)
        {
          for (int j_0 = i_0 + 1; j_0 < 4; j_0++)
          {
            for (int i_1 = 0; i_1 < 4; i_1++)
            {
              for (int j_1 = i_1 + 1; j_1 < 4; j_1++)
              {
                if ((i_1 == i_0 && j_1 == j_0) || verts[i_0] == verts[j_0] ||
                    verts[i_1] == verts[j_1])
                {
                  continue;
                }
                const auto e_ij_0 = (*verts[j_0] - *verts[i_0]).normalized();
                const auto e_ij_1 = (*verts[j_1] - *verts[i_1]).normalized();
                const auto dot = std::abs(e_ij_0.dot(e_ij_1));
                if (dot < dot_min)
                {
                  if constexpr (false)
                  {
                    fmt::print("i_0 {} i_1 {} j_0 {} j_1 {}\n", i_0, i_1, j_0, j_1);
                    fmt::print("e_ij_0 {}, e_ij_1 {}!\n", e_ij_0, e_ij_1);
                  }
                  dot_min = dot;
                  e_0 = e_ij_0;
                  e_1 = e_ij_1;
                  if (dot_min < rel_tol)
                  {
                    return std::make_pair(e_0, e_1);
                  }
                }
              }
            }
          }
        }
        return std::make_pair(e_0, e_1);
      }();

      if constexpr (false)
      {
        fmt::print("e_0 {}, e_1 {}!\n", e_0, e_1);
      }

      Vector3dMap(box.axes[0].data()) = e_0;
      Vector3dMap(box.axes[1].data()) = e_1;
      Vector3dMap(box.axes[2].data()) =
          box.planar ? Eigen::Vector3d::Zero() : e_0.cross(e_1);
    }

    // Scale axes by length of the box in each direction.
    std::array<double, 3> l = {0.0};
    for (const auto &v : {v_000, v_001, v_011, v_111})
    {
      const auto v_0 = v - Vector3dMap(box.center.data());
      l[0] = std::max(l[0], std::abs(v_0.dot(Vector3dMap(box.axes[0].data()))));
      l[1] = std::max(l[1], std::abs(v_0.dot(Vector3dMap(box.axes[1].data()))));
      l[2] = std::max(l[2], std::abs(v_0.dot(Vector3dMap(box.axes[2].data()))));
    }
    Vector3dMap(box.axes[0].data()) *= l[0];
    Vector3dMap(box.axes[1].data()) *= l[1];
    Vector3dMap(box.axes[2].data()) *= l[2];

    // Make sure the longest dimension comes first.
    std::sort(box.axes.begin(), box.axes.end(), [](const auto &x, const auto &y)
              { return CVector3dMap(x.data()).norm() > CVector3dMap(y.data()).norm(); });
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, box.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, box.axes.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &box.planar, dominant_rank, comm);

  return box;
}

// Use 4 points to define a sphere in 3D. If the points are coplanar, 3 of them are used to
// define a circle which is interpreted as the equator of the sphere. We assume the points
// are unique and not collinear.
BoundingBall SphereFromPoints(const std::vector<std::size_t> &indices,
                              const std::vector<Eigen::Vector3d> &vertices)
{
  // Given 0 or 1 points, just return a radius of 0.
  MFEM_VERIFY(
      indices.size() <= 4,
      "Determining a sphere in 3D requires 4 points (and a circle requires 3 points)!");
  BoundingBall ball;
  ball.planar = (indices.size() < 4);
  if (indices.size() < 2)
  {
    ball.origin = Eigen::Vector3d::Zero();
    ball.radius = 0.0;
    return ball;
  }

  // For two points, construct a circle with the segment as its diameter. This could also
  // handle the collinear case for more than 2 points.
  if (indices.size() == 2)
  {
    ball.origin = 0.5 * (vertices[indices[0]] + vertices[indices[1]]);
    ball.radius = (vertices[indices[0]] - ball.origin).norm();
    return ball;
  }

  // Check for coplanarity.
  constexpr double rel_tol = 1.0e-6;
  const Eigen::Vector3d AB = vertices[indices[1]] - vertices[indices[0]];
  const Eigen::Vector3d AC = vertices[indices[2]] - vertices[indices[0]];
  const Eigen::Vector3d ABAC = AB.cross(AC);
  Eigen::Vector3d AD = Eigen::Vector3d::Zero();
  if (!ball.planar)
  {
    AD = vertices[indices[3]] - vertices[indices[0]];
    ball.planar = (std::abs(AD.dot(ABAC)) < rel_tol * AD.norm() * ABAC.norm());
  }

  // Construct a circle passing through 3 points.
  // See: https://en.wikipedia.org/wiki/Circumcircle#Higher_dimensions.
  if (ball.planar)
  {
    ball.origin = (0.5 / ABAC.squaredNorm()) *
                  ((AB.squaredNorm() * AC) - (AC.squaredNorm() * AB)).cross(ABAC);
    ball.radius = ball.origin.norm();
    ball.origin += vertices[indices[0]];
#if defined(MFEM_DEBUG)
    const auto r1 = (vertices[indices[1]] - ball.origin).norm();
    const auto r2 = (vertices[indices[2]] - ball.origin).norm();
    MFEM_VERIFY((1.0 - rel_tol) * ball.radius < r1 && r1 < (1.0 + rel_tol) * ball.radius &&
                    (1.0 - rel_tol) * ball.radius < r2 &&
                    r2 < (1.0 + rel_tol) * ball.radius,
                "Invalid circle calculated from 3 points!");
#endif
    return ball;
  }

  // Construct a sphere passing through 4 points.
  // See: https://steve.hollasch.net/cgindex/geometry/sphere4pts.html.
  Eigen::Matrix3d C;
  Eigen::Vector3d d;
  const auto s = vertices[indices[0]].squaredNorm();
  C.row(0) = AB.transpose();
  C.row(1) = AC.transpose();
  C.row(2) = AD.transpose();
  d(0) = 0.5 * (vertices[indices[1]].squaredNorm() - s);
  d(1) = 0.5 * (vertices[indices[2]].squaredNorm() - s);
  d(2) = 0.5 * (vertices[indices[3]].squaredNorm() - s);
  ball.origin = C.inverse() * d;  // 3x3 matrix inverse might be faster than general LU
                                  // if Eigen uses the explicit closed-form solution
  ball.radius = (vertices[indices[0]] - ball.origin).norm();
#if defined(MFEM_DEBUG)
  const auto r1 = (vertices[indices[1]] - ball.origin).norm();
  const auto r2 = (vertices[indices[2]] - ball.origin).norm();
  const auto r3 = (vertices[indices[3]] - ball.origin).norm();
  MFEM_VERIFY((1.0 - rel_tol) * ball.radius < r1 && r1 < (1.0 + rel_tol) * ball.radius &&
                  (1.0 - rel_tol) * ball.radius < r2 &&
                  r2 < (1.0 + rel_tol) * ball.radius &&
                  (1.0 - rel_tol) * ball.radius < r3 && r3 < (1.0 + rel_tol) * ball.radius,
              "Invalid sphere calculated from 3 points!");
#endif
  return ball;
}

BoundingBall Welzl(std::vector<std::size_t> P, std::vector<std::size_t> R,
                   const std::vector<Eigen::Vector3d> &vertices)
{
  // Base case.
  if (R.size() == 4 || P.empty())
  {
    return SphereFromPoints(R, vertices);
  }

  // Choose a p ∈ P randomly, and recurse for (P \ {p}, R). The set P has already been
  // randomized on input.
  const std::size_t p = P.back();
  P.pop_back();
  BoundingBall D = Welzl(P, R, vertices);

  // If p is outside the sphere, recurse for (P \ {p}, R U {p}).
  constexpr double rel_tol = 1.0e-6;
  if ((vertices[p] - D.origin).norm() >= (1.0 + rel_tol) * D.radius)
  {
    R.push_back(p);
    D = Welzl(P, R, vertices);
  }

  return D;
}

// Calculates a bounding ball from a point cloud using Welzl's algorithm, result is
// broadcast across all processes. We don't operate on the convex hull, since the number of
// points should be small enough that operating on the full set should be OK. If only three
// points are provided, the bounding circle is computed (likewise for if the points are
// coplanar).
BoundingBox BoundingBallFromPointCloud(MPI_Comm comm,
                                       const std::vector<Eigen::Vector3d> &vertices,
                                       int dominant_rank)
{
  BoundingBox ball;
  if (dominant_rank == Mpi::Rank(comm))
  {
    MFEM_VERIFY(vertices.size() >= 3,
                "A bounding ball requires a minimum of three vertices for this algorithm!");
    std::vector<std::size_t> indices(vertices.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Acceleration from https://informatica.vu.lt/journal/INFORMATICA/article/1251. Allow
    // for duplicate points and just add the 4 points to the end of the indicies list to be
    // considered first. The two points are not necessarily the maximizer of the distance
    // between all pairs, but they should be a good estimate.
    {
      auto p_1 = std::min_element(vertices.begin(), vertices.end(), EigenLE);
      auto p_2 = std::max_element(vertices.begin(), vertices.end(),
                                  [p_1](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                  { return (x - *p_1).norm() < (y - *p_1).norm(); });
      p_1 = std::max_element(vertices.begin(), vertices.end(),
                             [p_2](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                             { return (x - *p_2).norm() < (y - *p_2).norm(); });

      // Find the next point as the vertex furthest from the initial axis.
      const Eigen::Vector3d n_1 = (*p_2 - *p_1).normalized();
      auto p_3 = std::max_element(vertices.begin(), vertices.end(),
                                  [&](const auto &x, const auto &y)
                                  {
                                    return PerpendicularDistance({n_1}, *p_1, x) <
                                           PerpendicularDistance({n_1}, *p_1, y);
                                  });
      auto p_4 = std::max_element(vertices.begin(), vertices.end(),
                                  [p_3](const Eigen::Vector3d &x, const Eigen::Vector3d &y)
                                  { return (x - *p_3).norm() < (y - *p_3).norm(); });
      MFEM_VERIFY(p_3 != p_1 && p_3 != p_2 && p_4 != p_1 && p_4 != p_2,
                  "Vertices are degenerate!");

      // Start search with these points, which should be roughly extremal. With the search
      // for p_3 done in an orthogonal direction, p_1, p_2, p_3, and p_4 should all be
      // unique.
      std::swap(indices[indices.size() - 1], indices[p_1 - vertices.begin()]);
      std::swap(indices[indices.size() - 2], indices[p_2 - vertices.begin()]);
      std::swap(indices[indices.size() - 3], indices[p_3 - vertices.begin()]);
      std::swap(indices[indices.size() - 4], indices[p_4 - vertices.begin()]);
    }

    // Randomly permute the point set.
    {
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(indices.begin(), indices.end() - 4, g);
    }

    // Compute the bounding ball.
    BoundingBall min_ball = Welzl(indices, {}, vertices);
    Vector3dMap(ball.center.data()) = min_ball.origin;
    Vector3dMap(ball.axes[0].data()) = Eigen::Vector3d(min_ball.radius, 0.0, 0.0);
    Vector3dMap(ball.axes[1].data()) = Eigen::Vector3d(0.0, min_ball.radius, 0.0);
    Vector3dMap(ball.axes[2].data()) = Eigen::Vector3d(0.0, 0.0, min_ball.radius);
    ball.planar = min_ball.planar;
  }

  // Broadcast result to all processors.
  Mpi::Broadcast(3, ball.center.data(), dominant_rank, comm);
  Mpi::Broadcast(3 * 3, ball.axes.data()->data(), dominant_rank, comm);
  Mpi::Broadcast(1, &ball.planar, dominant_rank, comm);

  return ball;
}

// Compute a normal vector from an element transformation, optionally ensure aligned
// (| normal ⋅ align | > 0)
void Normal(mfem::ElementTransformation &T, mfem::Vector &normal,
            const mfem::Vector *const align)
{
  const mfem::IntegrationPoint &ip = mfem::Geometries.GetCenter(T.GetGeometryType());
  T.SetIntPoint(&ip);
  mfem::CalcOrtho(T.Jacobian(), normal);
  normal /= normal.Norml2();
  if (align && (normal * (*align) < 0))
  {
    normal *= -1;
  }
}

// Compute the centroid of a container of vertices.
template <typename T>
mfem::Vector ComputeCentroid(const std::unique_ptr<mfem::Mesh> &mesh, const T &vertidxs)
{
  mfem::Vector centroid(mesh->SpaceDimension());
  centroid = 0.0;
  int c = 0;
  for (const int v : vertidxs)
  {
    centroid += mfem::Vector(mesh->GetVertex(v), 3);
    c++;
  }
  return centroid /= c;
}

// Compute the normal vector for a set of elements. If "inside" is true, normal will
// point inside the mesh, otherwise it will point outside the mesh.
template <typename T>
mfem::Vector ComputeNormal(const std::unique_ptr<mfem::Mesh> &mesh, const T &elem_set,
                           bool inside, bool check_planar = true)
{
  const int sdim = mesh->SpaceDimension();
  mfem::IsoparametricTransformation trans;
  mfem::Vector normal(sdim), last_normal(sdim), align(sdim);
  normal = 0.0;
  last_normal = 0.0;
  mfem::Array<int> vert_bdr;

  // Ensure that the computed normal points "inside" or "outside".
  auto Alignment = [&](int el, auto &align)
  {
    int eladj, info;
    mesh->GetBdrElementAdjacentElement(el, eladj, info);
    mesh->GetElementCenter(eladj, align);
    mesh->GetBdrElementVertices(el, vert_bdr);
    align -= ComputeCentroid(mesh, vert_bdr);
    if (!inside)  // align points inwards
    {
      align *= -1;
    }
  };

  for (auto elem : elem_set)
  {
    Alignment(elem, align);
    mesh->GetBdrElementTransformation(elem, &trans);
    mesh::Normal(trans, normal, &align);
    if (!check_planar)
    {
      break;  // If not checking planar, use the first.
    }
    MFEM_VERIFY((last_normal * last_normal == 0.0) || ((last_normal * normal - 1) < 1e-8),
                "Periodic boundary mapping is only supported for planar boundaries!");
    last_normal = normal;
  }
  return normal;
}

struct Frame
{
  Frame(const mfem::Vector &o) : origin(o)
  {
    for (auto &x : basis)
    {
      x.SetSize(3);
      x = 0.0;
    }
  }
  mfem::Vector origin;
  std::array<mfem::Vector, 3> basis;
};

Frame Find3DFrame(std::unique_ptr<mfem::Mesh> &mesh,
                  const std::unordered_set<int> &vertidxs, const mfem::Vector &centroid,
                  const mfem::Vector &normal, double mesh_dim)
{
  Frame frame(centroid);
  frame.basis[0] = normal;

  // For each point, compute its distance to the centroid.
  std::map<int, std::vector<int>, std::greater<int>> dist2points;
  for (const int v : vertidxs)
  {
    auto dist = centroid.DistanceTo(mesh->GetVertex(v));
    // Convert dist to integer to avoid floating point differences.
    dist2points[std::round(dist / mesh_dim * 1e8)].push_back(v);
  }

  for (const auto &[dist, verts] : dist2points)
  {
    if (verts.size() > 1 || dist == 0)
    {
      continue;
    }
    frame.basis[1] = mesh->GetVertex(verts.front());
    frame.basis[1] -= centroid;
    frame.basis[1] /= frame.basis[1].Norml2();
    break;
  }

  // Define final point by computing the cross product.
  frame.basis[0].cross3D(frame.basis[1], frame.basis[2]);

  return frame;
}

// Calculate the rotation matrix between two vectors.
void ComputeRotation(const mfem::Vector &normal1, const mfem::Vector &normal2,
                     mfem::DenseMatrix &transformation)
{
  mfem::DenseMatrix R(3), vx(3), vx2(3);

  mfem::Vector v(normal1.Size());
  normal1.cross3D(normal2, v);
  double c = normal1 * normal2;

  vx(0, 1) = -v[2];
  vx(0, 2) = v[1];
  vx(1, 0) = v[2];
  vx(1, 2) = -v[0];
  vx(2, 0) = -v[1];
  vx(2, 1) = v[0];

  R(0, 0) = R(1, 1) = R(2, 2) = 1.0;
  R += vx;
  Mult(vx, vx, vx2);
  if (std::abs(1.0 + c) > 1e-8)
  {
    vx2.Set(1.0 / (1.0 + c), vx2);
  }
  R += vx2;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      transformation(i, j) = R(i, j);
    }
  }
}

mfem::DenseMatrix ComputeAffineTransformationMatrix(const Frame &donor,
                                                    const Frame &receiver)
{

  mfem::DenseMatrix A(4, 4);
  A = 0.0;
  if (donor.basis[1].Norml2() > 0.0 && receiver.basis[1].Norml2() > 0.0)
  {
    // Stably compute the rotation matrix from unit vectors.
    Eigen::Matrix3d source, target;
    Eigen::Vector3d source_centroid, target_centroid;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        source(i, j) = donor.basis[j](i);
        target(i, j) = receiver.basis[j](i);
      }
      source_centroid(i) = donor.origin(i);
      target_centroid(i) = receiver.origin(i);
    }
    Eigen::Matrix3d R = source * target.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();

    // Account for possible reflection in R (det(R) = -1).
    Eigen::DiagonalMatrix<double, 3> Diag(1.0, 1.0, R.determinant());
    R = svd.matrixV() * Diag * svd.matrixU().transpose();

    // Compute translation and form transformation matrix.
    const Eigen::Vector3d translation = target_centroid - R * source_centroid;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        A(i, j) = R(i, j);
      }
      A(i, 3) = translation(i);
    }
    A(3, 3) = 1.0;
  }
  else
  {
    // If the donor or receiver basis is ambiguous, we assume no rotation around the
    // normals, and that the rotation comes only from realigning normal vectors.
    ComputeRotation(donor.basis[0], receiver.basis[0], A);
    for (int i = 0; i < 3; i++)
    {
      A(i, 3) = receiver.origin(i) - donor.origin(i);
    }
    A(3, 3) = 1.0;
  }

  return A;
}

// Create the vertex mapping between sets of donor and receiver pts related
// by an affine transformation matrix.
std::vector<int> CreatePeriodicVertexMapping(std::unique_ptr<mfem::Mesh> &mesh,
                                             const std::unordered_set<int> &donor_v,
                                             const std::unordered_set<int> &receiver_v,
                                             const mfem::DenseMatrix &transform,
                                             double tol = 1e-6)
{
  // Similar to MFEM's CreatePeriodicVertexMapping, maps from replica to primary vertex.
  std::unordered_map<int, int> replica2primary;

  // KD-tree containing all the receiver points.
  mfem::KDTree3D kdtree;
  for (const int v : receiver_v)
  {
    kdtree.AddPoint(mesh->GetVertex(v), v);
  }
  kdtree.Sort();

  // Loop over donor points and find the corresponding receiver point.
  mfem::Vector from(4), to(4);
  for (int vi : donor_v)
  {
    // TODO: mfem patch to allow SetVector direct from pointer.
    std::copy(mesh->GetVertex(vi), mesh->GetVertex(vi) + 3, from.begin());
    from[3] = 1.0;             // reset
    transform.Mult(from, to);  // receiver = transform * donor

    const int vj = kdtree.FindClosestPoint(to.GetData());
    std::copy(mesh->GetVertex(vj), mesh->GetVertex(vj) + 3, from.begin());
    from -= to;  // Check that the loaded vertex is identical to the transformed
    MFEM_VERIFY(from.Norml2() < tol,
                "Could not match points on periodic boundaries, transformed donor point "
                "does not correspond to a receiver point!");
    MFEM_VERIFY(replica2primary.find(vj) == replica2primary.end(),
                "Could not match points on periodic boundaries, multiple donor points map "
                "to the same receiver point!")
    replica2primary[vj] = vi;
  }

  std::vector<int> v2v(mesh->GetNV());
  std::iota(v2v.begin(), v2v.end(), 0);
  for (const auto &[r, p] : replica2primary)
  {
    v2v[r] = p;
  }
  return v2v;
}

// Determine the vertex mapping between donor and receiver boundary attributes.
// Uses the translation vector or affine transformation matrix specified in the
// configuration file. If not provided, attempts to automatically detect the
// affine transformation between donor and receiver boundary vertices.
std::vector<int>
DeterminePeriodicVertexMapping(std::unique_ptr<mfem::Mesh> &mesh,
                               const struct palace::config::PeriodicData &data,
                               const double tol)
{
  // Get mesh dimensions, will be used to define a reasonable tolerance in mesh units.
  mfem::Vector bbmin, bbmax;
  mesh->GetBoundingBox(bbmin, bbmax);
  bbmax -= bbmin;
  const double mesh_dim = bbmax.Norml2();
  const double mesh_tol = tol * mesh_dim;

  // Identify donor and receiver vertices and elements.
  const auto &da = data.donor_attributes, &ra = data.receiver_attributes;
  std::unordered_set<int> bdr_v_donor, bdr_v_receiver;
  std::unordered_set<int> bdr_e_donor, bdr_e_receiver;
  for (int be = 0; be < mesh->GetNBE(); be++)
  {
    int attr = mesh->GetBdrAttribute(be);
    auto donor = std::find(da.begin(), da.end(), attr) != da.end();
    auto receiver = std::find(ra.begin(), ra.end(), attr) != ra.end();
    if (donor || receiver)
    {
      int el, info;
      mesh->GetBdrElementAdjacentElement(be, el, info);
      mfem::Array<int> vertidxs;
      mesh->GetBdrElementVertices(be, vertidxs);
      (donor ? bdr_e_donor : bdr_e_receiver).insert(be);
      (donor ? bdr_v_donor : bdr_v_receiver).insert(vertidxs.begin(), vertidxs.end());
    }
  }

  MFEM_VERIFY(bdr_v_donor.size() == bdr_v_receiver.size(),
              "Different number of "
              "vertices on donor and receiver boundaries. Cannot create periodic mesh.");

  // Check if mesh has enough elements in periodic direction. MFEM's periodicity
  // fails for meshes with <=2 elements in the period direction.
  // Compare the number of mesh elements to the number of periodic boundary
  // elements.
  const int num_periodic_bc_elems = bdr_e_donor.size() + bdr_e_receiver.size();
  mfem::Array<mfem::Geometry::Type> geoms;
  mesh->GetGeometries(3, geoms);
  if (geoms.Size() == 1 && geoms[0] == mfem::Geometry::TETRAHEDRON)
  {
    // Pure tet mesh.
    MFEM_VERIFY(mesh->GetNE() > 3 * num_periodic_bc_elems,
                "Not enough mesh elements in periodic direction!");
  }
  else
  {
    // No tets.
    MFEM_VERIFY(mesh->GetNE() > num_periodic_bc_elems,
                "Not enough mesh elements in periodic direction!");
  }

  // Determine the affine transformation between donor and receiver points.
  // Use the translation vector or affine transformation matrix if provided
  // in the config file, otherwise automatically detect the transformation.
  mfem::DenseMatrix transformation(4);
  if (std::any_of(data.affine_transform.begin(), data.affine_transform.end(),
                  [](auto y) { return std::abs(y) > 0.0; }))
  {
    // Use user-provided affine transformation matrix.
    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        transformation(i, j) = data.affine_transform[i * 4 + j];  // row major conversion
      }
    }
  }
  else
  {
    // Automatically detect transformation.
    // Compute the centroid for each boundary.
    auto donor_centroid = ComputeCentroid(mesh, bdr_v_donor);
    auto receiver_centroid = ComputeCentroid(mesh, bdr_v_receiver);

    // Compute the normal vector for each boundary.
    auto donor_normal = ComputeNormal(mesh, bdr_e_donor, true);
    auto receiver_normal = ComputeNormal(mesh, bdr_e_receiver, false);

    // Return empty mapping if centroids and normal vectors are the same (up to a sign).
    mfem::Vector diff = donor_centroid;
    diff -= receiver_centroid;
    double dot = donor_normal * receiver_normal;
    if (diff.Norml2() < mesh_tol && std::abs(std::abs(dot) - 1.0) < mesh_tol)
    {
      return {};
    }

    // Compute a frame (origin, normal, and two in plane points) for each boundary.
    auto donor_frame =
        Find3DFrame(mesh, bdr_v_donor, donor_centroid, donor_normal, mesh_dim);
    auto receiver_frame =
        Find3DFrame(mesh, bdr_v_receiver, receiver_centroid, receiver_normal, mesh_dim);

    // Compute the affine transformation matrix.
    transformation = ComputeAffineTransformationMatrix(donor_frame, receiver_frame);
  }
  return CreatePeriodicVertexMapping(mesh, bdr_v_donor, bdr_v_receiver, transformation,
                                     mesh_tol);
}

}  // namespace palace::mesh