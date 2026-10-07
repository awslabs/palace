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
// numbering and the partition: the DOF coordinates for H1 (width 3), the moments dof(e_b)
// and dof(x_a e_b) for H(curl), which identify an edge or face DOF up to an orientation
// flip (a flip negates all of them; width 12). index[i] is the position of local true DOF i
// in the set (-1 if not in it), n the set size. Collective.
int SignatureWidth(const mfem::ParFiniteElementSpace &fespace);
std::vector<double> TrueDofSignatures(const mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<int> &index, int n);

// The saved DOFs of the current ones by signature: current DOF g is saved DOF perm[g], with
// orientation sgn[g] (+/-1 for H(curl), +1 for H1), so that a saved vector v maps to
// v_cur[g] = sgn[g] v[perm[g]]. Aborts unless the match is a bijection within a relative
// tolerance of the signature scale.
void MatchSignatures(const std::vector<double> &current, const std::vector<double> &saved,
                     int width, bool signed_match, std::vector<int> &perm,
                     std::vector<double> &sgn);

// Change of basis between the saved and the current DOF sets: current DOF values (of a
// field) u_cur = M u_saved, with M a signed permutation (the matches of MatchSignatures)
// except on groups of DOFs at a common point whose basis depends on the partition (the face
// DOFs of second-order Nédélec elements on tetrahedra, oriented by their face), where M has
// a dense block. Dual quantities (right-hand sides, Schur complements) transform with
// M^-T: rows[g] holds row g of M^-T as (saved DOF, coefficient). Aborts unless the current
// signatures equal M times the saved ones within a relative tolerance.
struct SignatureMap
{
  std::vector<std::vector<std::pair<int, double>>> rows;

  // x_cur = M^-T x_saved for a dual vector, and S_cur = M^-T S_saved M^-1 for a symmetric
  // dual matrix (column-major), of size n.
  template <typename T>
  std::vector<T> Dual(const T *x) const;
  template <typename T>
  std::vector<T> DualMatrix(const std::vector<T> &S) const;
};
SignatureMap MatchSignatureBasis(const std::vector<double> &current,
                                 const std::vector<double> &saved, int width);

template <typename T>
std::vector<T> SignatureMap::Dual(const T *x) const
{
  std::vector<T> y(rows.size(), T(0.0));
  for (std::size_t g = 0; g < rows.size(); g++)
  {
    for (const auto &[a, c] : rows[g])
    {
      y[g] += c * x[a];
    }
  }
  return y;
}

template <typename T>
std::vector<T> SignatureMap::DualMatrix(const std::vector<T> &S) const
{
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
