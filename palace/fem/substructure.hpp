// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <utility>
#include <vector>
#include <mfem.hpp>

namespace palace
{

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
  // S_cur = M^-T S_saved M^-1 for a symmetric dual matrix S (n x n).
  template <typename T>
  std::vector<T> DualRows(const std::vector<T> &X, int k) const;
  template <typename T>
  std::vector<T> DualMatrix(const std::vector<T> &S) const;
};

// The map from the saved to the current H(curl) signatures (width 12), aborting unless the
// current signatures equal M times the saved ones within a relative tolerance.
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

template <typename T>
std::vector<T> SignatureMap::DualMatrix(const std::vector<T> &S) const
{
  // (M^-T S M^-1)_ij = sum_ab (M^-T)_ia S_ab (M^-T)_jb, in one pass (S symmetric).
  const std::size_t n = rows.size();
  std::vector<T> Y(n * n, T(0.0));
  for (std::size_t j = 0; j < n; j++)
  {
    for (std::size_t i = 0; i < n; i++)
    {
      T v(0.0);
      for (const auto &[a, ca] : rows[i])
      {
        for (const auto &[b, cb] : rows[j])
        {
          v += ca * cb * S[static_cast<std::size_t>(b) * n + a];
        }
      }
      Y[j * n + i] = v;
    }
  }
  return Y;
}

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
