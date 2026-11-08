// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class FiniteElementSpace;

// Abort unless every domain attribute of the mesh is in the region or the environment, and
// every region attribute is one of the mesh (the mesh of the region written for online runs
// holds only part of the environment).
void CheckSubstructuringAttributes(const mfem::ParMesh &mesh,
                                   const std::vector<int> &region_attributes,
                                   const std::vector<int> &environment_attributes);

// Parallel-capable interface identification in true-DOF space. Given the parent finite
// element space and the region/environment domain attribute sets, mark each parent true DOF
// touched by region element support, by environment support, and by both (the shared
// interface). Works on any MPI partition: L-vector attribute marks are reduced to true DOFs
// across ranks via the parallel prolongation. Returns the global interface true-DOF count.
int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker,
                          mfem::Array<int> &interface_marker);

// The same on a Palace space, which may be broken across interior boundaries: a version of
// a DOF is marked by the elements reading it.
int MarkInterfaceTrueDofs(const FiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker,
                          mfem::Array<int> &interface_marker);

// Geometric signatures of a set of true DOFs, replicated on all ranks (n x width,
// row-major, width = SignatureWidth(fespace)), identifying them independently of the DOF
// numbering: the DOF coordinates for H1 (width 3), the moments dof(e_b) and dof(x_a e_b)
// for H(curl) (width 12), which identify an edge DOF up to an orientation flip (a flip
// negates all of them) and give the point of a point-tangent DOF. index[i] is the position
// of local true DOF i in the set (-1 if not in it), n the set size. Collective.
int SignatureWidth(const mfem::ParFiniteElementSpace &fespace);
std::vector<double> TrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<int> &index, int n);

// Change of basis between the saved and the current DOF sets: current DOF values (of a
// field) u_cur = M u_saved, with M a signed permutation (signed one-to-one matches) except
// on groups of DOFs at a common point whose basis depends on the partition (the face DOFs
// of second-order Nédélec elements on tetrahedra, oriented by their face), where M has a
// dense block. Dual quantities (right-hand sides, Schur complements) transform with M^-T:
// rows[g] holds row g of M^-T as (saved DOF, coefficient).
struct SignatureMap
{
  std::vector<std::vector<std::pair<int, double>>> rows;

  // X_cur = M^-T X_saved for dual vectors X (n x k, row-major), and
  // S_cur = M^-T S_saved M^-1 for a symmetric dual matrix S (n x n) from its saved rows in
  // order, read_row(a, x) filling row a: one saved row is held at a time.
  template <typename T>
  std::vector<T> DualRows(const std::vector<T> &X, int k) const;
  template <typename T, typename ReadRow>
  std::vector<T> DualMatrix(ReadRow &&read_row) const;
};

// Which of a set of DOFs (H1 or H(curl) signatures of the given width) lie at the point of
// one of a saved set of DOFs (of saved_width, 13 with the side moments of split DOFs),
// within the tolerance of MatchSignatureBasis.
std::vector<char> MatchSignaturePoints(const std::vector<double> &candidates, int width,
                                       const std::vector<double> &saved, int saved_width);

// The map from the saved to the current signatures, aborting unless the current signatures
// equal M times the saved ones within a relative tolerance (H1: a permutation). A 13th
// column (width 13) tells the versions of a DOF of a broken space apart, which must match
// one to one: a side moment, zero for the other DOFs, which flips with the orientation as
// the others.
SignatureMap MatchSignatureBasis(const std::vector<double> &current,
                                 const std::vector<double> &saved, int width);

template <typename T>
std::vector<T> SignatureMap::DualRows(const std::vector<T> &X, int k) const
{
  std::vector<T> Y(X.size(), T(0.0));
  for (std::size_t g = 0; g < rows.size(); g++)
  {
    for (const auto &[a, c] : rows[g])
    {
      for (int q = 0; q < k; q++)
      {
        Y[g * k + q] += c * X[static_cast<std::size_t>(a) * k + q];
      }
    }
  }
  return Y;
}

template <typename T, typename ReadRow>
std::vector<T> SignatureMap::DualMatrix(ReadRow &&read_row) const
{
  // (M^-T S M^-1)_gh = sum_a (M^-T)_ga z_h, z_h = sum_b (M^-T)_hb S_ab: each saved row a,
  // mapped on its columns (z), adds to the current rows g with (M^-T)_ga != 0.
  const std::size_t n = rows.size();
  std::vector<std::vector<std::pair<std::size_t, double>>> to(n);
  for (std::size_t g = 0; g < n; g++)
  {
    for (const auto &[a, c] : rows[g])
    {
      to[a].emplace_back(g, c);
    }
  }
  std::vector<T> Y(n * n, T(0.0)), x(n), z(n);
  for (std::size_t a = 0; a < n; a++)
  {
    read_row(a, x.data());
    for (std::size_t h = 0; h < n; h++)
    {
      T v(0.0);
      for (const auto &[b, c] : rows[h])
      {
        v += c * x[b];
      }
      z[h] = v;
    }
    for (const auto &[g, c] : to[a])
    {
      T *y = &Y[g * n];
      for (std::size_t h = 0; h < n; h++)
      {
        y[h] += c * z[h];
      }
    }
  }
  return Y;
}

// Local true DOFs touched by the boundary elements with the given attributes (a true DOF
// shared between ranks is marked on its owner).
std::vector<char> BoundaryTrueDofs(mfem::ParFiniteElementSpace &fespace,
                                   const std::vector<int> &bdr_attributes);

// A with the rows and columns of the given local true DOFs replaced by the identity. The
// identity is added rather than set by the elimination: a matrix assembled on the device
// stores no diagonal entry in a row without element contributions.
std::unique_ptr<mfem::HypreParMatrix> EliminateToIdentity(const mfem::HypreParMatrix &A,
                                                          const mfem::Array<int> &dofs);

// The boundary attributes of the mesh whose elements all border elements of the region
// only, so that a field known in the region determines quantities on them. Collective.
std::vector<int> BoundaryAttributesInRegion(const mfem::ParMesh &mesh,
                                            const std::vector<int> &region_attributes);

// Whether a point (mesh coordinates) lies in an element of the region. Collective.
bool PointInRegion(const mfem::ParMesh &mesh, const std::array<double, 3> &x,
                   const std::vector<int> &region_attributes);

// The mesh of the region written next to a saved substructuring model: <model>_region.mesh.
std::string RegionMeshPath(const std::string &model_path);

// Write a mesh of the region for online substructuring runs, from rank 0 in mesh length
// units (coordinates times length_scale): the region, and the environment elements touching
// the marked local true DOFs, with their boundary elements. The local boundary elements
// marked in plain_bdr (if not empty) are written with the attribute of the boundary faces
// the submesh generates. Returns the local environment elements written. Collective.
std::vector<char> WriteRegionMesh(const std::string &path,
                                  mfem::ParFiniteElementSpace &fespace,
                                  const std::vector<int> &region_attributes,
                                  const std::vector<int> &environment_attributes,
                                  const std::vector<char> &keep, double length_scale,
                                  const std::vector<char> &plain_bdr = {});

// The number, volume, and volume times the centroid distance from the origin and the
// centroid of the elements with the given attributes (and mark[e], if not empty), summed
// over the processes: a description of a set of elements independent of the partition and
// of the element numbering. Collective.
std::array<double, 6> ElementMoments(mfem::ParMesh &mesh,
                                     const std::vector<int> &attributes,
                                     const std::vector<char> &mark = {});

// Whether two ElementMoments describe the same elements within the relative tolerance (the
// centroid moments relative to the distance moment, as they may vanish).
bool MomentsMatch(const std::array<double, 6> &a, const std::array<double, 6> &b,
                  double tol);

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
