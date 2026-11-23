// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_SUBSTRUCTURE_HPP
#define PALACE_FEM_SUBSTRUCTURE_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json_fwd.hpp>

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
// touched by region element support and by environment support (both on the shared
// interface). Works on any MPI partition: L-vector attribute marks are reduced to true DOFs
// across ranks via the parallel prolongation. Returns the global interface true-DOF count.
int MarkInterfaceTrueDofs(mfem::ParFiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker);

// The same on a Palace space, which may be broken across interior boundaries: a version of
// a DOF is marked by the elements reading it.
int MarkInterfaceTrueDofs(const FiniteElementSpace &parent_fespace,
                          const mfem::Array<int> &region_attrs,
                          const mfem::Array<int> &environment_attrs,
                          mfem::Array<int> &region_marker, mfem::Array<int> &env_marker);

// The classes of the local true DOFs of a substructured space, from the markers of
// MarkInterfaceTrueDofs, the Dirichlet DOFs and, on the mesh of the region, the region's
// DOFs at the saved interface DOFs (cut, see RegionMeshInterface; empty otherwise): Γ (the
// DOFs of both sides, or cut), and the environment and region interiors without Γ, none of
// them Dirichlet.
struct InterfaceClasses
{
  std::vector<char> gamma, env_int, region_int;
};
InterfaceClasses ClassifyInterfaceTrueDofs(const mfem::Array<int> &region_marker,
                                           const mfem::Array<int> &env_marker,
                                           const std::vector<char> &dirichlet,
                                           const std::vector<char> &cut);

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

  // X_cur = M^-T X_saved for dual vectors X (n x k, row-major), and the rows [g0, g0 + m)
  // of S_cur = M^-T S_saved M^-1 for a symmetric dual matrix S (n x n, m x n row-major)
  // from the saved rows they combine, SavedRows(g0, m), read_row(a, x) filling each saved
  // row a in their order: one saved row is held at a time.
  template <typename T>
  std::vector<T> DualRows(const std::vector<T> &X, int k) const;
  std::vector<int> SavedRows(int g0, int m) const;
  template <typename T, typename ReadRow>
  std::vector<T> DualMatrix(int g0, int m, ReadRow &&read_row) const;
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
std::vector<T> SignatureMap::DualMatrix(int g0, int m, ReadRow &&read_row) const
{
  // (M^-T S M^-1)_gh = sum_a (M^-T)_ga z_h, z_h = sum_b (M^-T)_hb S_ab: each saved row a,
  // mapped on its columns (z), adds to the rows g with (M^-T)_ga != 0.
  const std::size_t n = rows.size();
  const std::vector<int> saved = SavedRows(g0, m);
  std::vector<std::vector<std::pair<std::size_t, double>>> to(saved.size());
  for (int g = g0; g < g0 + m; g++)
  {
    for (const auto &[a, c] : rows[g])
    {
      const auto i = std::lower_bound(saved.begin(), saved.end(), a) - saved.begin();
      to[i].emplace_back(g - g0, c);
    }
  }
  std::vector<T> Y(static_cast<std::size_t>(m) * n, T(0.0)), x(n), z(n);
  for (std::size_t a = 0; a < saved.size(); a++)
  {
    read_row(saved[a], x.data());
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

// The mesh of the region written next to a saved substructuring model: <model>_region.mesh.
std::string RegionMeshPath(const std::string &model_path);

// Write a mesh of the region for online substructuring runs from a conforming mesh, from
// rank 0 in mesh length units (coordinates times length_scale): the region, and the
// environment elements touching
// the marked local true DOFs, with their boundary elements. The local boundary elements
// marked in plain_bdr (if not empty) are written with the attribute of the boundary faces
// the submesh generates. Returns the local environment elements written. Collective.
std::vector<char> WriteRegionMesh(const std::string &path,
                                  mfem::ParFiniteElementSpace &fespace,
                                  const std::vector<int> &region_attributes,
                                  const std::vector<int> &environment_attributes,
                                  const std::vector<char> &keep, double length_scale,
                                  const std::vector<char> &plain_bdr = {});

// Whether two moments of AttributeMoments describe the same elements within the relative
// tolerance (the centroid moments relative to the distance moment, as they may vanish).
bool MomentsMatch(const std::array<double, 6> &a, const std::array<double, 6> &b,
                  double tol);

// For the elements with each of the given attributes (and mark[e], if not empty), in the
// order of the attributes, their number, volume, and volume times the distance of their
// center from the origin and times their center, summed over the processes: a description
// of the elements independent of the partition and of their numbering. Collective.
std::vector<std::array<double, 6>> AttributeMoments(mfem::ParMesh &mesh,
                                                    const std::vector<int> &attributes,
                                                    const std::vector<char> &mark = {});

// Add an element of volume (or area) v and center c to moments of AttributeMoments.
void AddMoments(std::array<double, 6> &m, double v, const mfem::Vector &c);

// Whether two moments of AttributeMoments agree within the relative tolerance in what a
// refinement of their elements keeps: their volumes and first moments.
bool RefinedMomentsMatch(const std::array<double, 6> &a, const std::array<double, 6> &b,
                         double tol);

// The boundary attributes of the boundary elements next to an element with one of the given
// domain attributes (on either side, for an interior boundary), sorted. Collective.
std::vector<int> BoundaryAttributesNextTo(mfem::ParMesh &mesh,
                                          const std::vector<int> &attributes);

// The moments of AttributeMoments of the boundary elements with each of the given boundary
// attributes whose neighboring elements all have one of the given domain attributes (with
// their areas for volumes; a refinement of the other domains does not change them), in the
// order of the boundary attributes. Collective.
std::vector<std::array<double, 6>>
BoundaryAttributeMoments(mfem::ParMesh &mesh, const std::vector<int> &bdr_attributes,
                         const std::vector<int> &attributes);

// A hash of the configuration of a substructured problem on its environment: the problem
// type, the length unit, the order and quadrature of the solver, the cracking of the mesh,
// and the entries of the materials on the given domain attributes and of the boundaries
// (but postprocessing) on the given boundary attributes, which are the entries listing one
// of them (in a list whose key ends in "Attributes", also in their "Elements", and in
// groups of entries), with these lists cut to them. The entries are taken as parsed where
// they can be written in several ways (the tensors of the materials, the numbers, and the
// aliases of "PEC" and "PMC"), and lists and entries are sorted, so that their order in the
// configuration does not matter. Without the excitations, which are matched to the modes of
// a model by their own fingerprints, the hash leaves out the indices of the entries and the
// sources (surface currents and flux loops): it is the configuration the operator of the
// environment depends on.
std::uint64_t EnvironmentConfigHash(const nlohmann::json &config,
                                    const std::vector<int> &attributes,
                                    const std::vector<int> &bdr_attributes,
                                    bool excitations);

// The messages of an environment that differs from that of a saved substructuring model:
// as its operator shows, by its mesh, and by its configuration.
inline constexpr const char *kEnvironmentDiffers =
    "The environment differs from the one the saved substructuring model was condensed "
    "from (environment mesh, materials, boundary conditions, order or problem type "
    "changed): rerun offline (without \"ModelFile\") to condense it again!";
inline constexpr const char *kEnvironmentMeshDiffers =
    "The environment differs from the one the saved substructuring model was condensed "
    "from: its mesh changed. Rerun offline (without \"ModelFile\") to condense it again!";
inline constexpr const char *kEnvironmentConfigDiffers =
    "The environment differs from the one the saved substructuring model was condensed "
    "from: its attributes, materials, boundary conditions or discretization changed. Rerun "
    "offline (without \"ModelFile\") to condense it again!";

// The message of excitations of the environment which differ from those of a saved
// substructuring model, on the mesh of the region (the full mesh computes their modes).
inline constexpr const char *kEnvironmentExcitationsDiffer =
    "The excitations of the environment differ from those of the saved substructuring "
    "model: run on the full mesh!";

// A description of the environment of a substructured problem which does not depend on the
// partition or on the numbering, for an online run to check that its environment is the one
// a model was condensed from.
struct EnvironmentRecord
{
  // The environment attributes, and the boundary attributes next to them.
  std::vector<int> attributes, bdr_attributes;

  // The AttributeMoments of the environment elements and of those written in the mesh of
  // the region, and the BoundaryAttributeMoments of the boundary attributes.
  std::vector<std::array<double, 6>> moments, region_moments, bdr_moments;

  // The EnvironmentConfigHash on these attributes, without and with the excitations.
  std::uint64_t config_hash = 0, excitations_hash = 0;

  // The record of the environment with the given attributes, with the configuration of the
  // run, and written[e] marking the elements written in the mesh of the region. Collective.
  static EnvironmentRecord Make(mfem::ParMesh &mesh, const std::vector<int> &attributes,
                                const nlohmann::json &config,
                                const std::vector<char> &written);

  // The record as bytes, and back (nullopt if the bytes are not a record).
  std::vector<char> Serialize() const;
  static std::optional<EnvironmentRecord> Deserialize(const std::vector<char> &bytes);

  // Online, check that the mesh holds the environment within the relative tolerance of the
  // moments: all its elements and their boundaries (those the configuration names by their
  // attributes, the others in the order of theirs), or the elements written in the mesh of
  // the region (that mesh, or a re-meshing of its region, which is not refined again).
  // Returns whether it is the mesh of the region. Collective.
  bool CheckMesh(mfem::ParMesh &mesh, const nlohmann::json &config, double tol) const;

  // Online, check that the environment attributes of the run are those of the record (on
  // the mesh of the region, some of them), and that the configuration on it is unchanged.
  // Returns whether the excitations of the environment are those of the record.
  bool CheckConfig(const nlohmann::json &config,
                   const std::vector<int> &environment_attributes, bool region_mesh) const;
};

// The message of an a priori refinement of the mesh of the region, written refined.
inline constexpr const char *kRegionMeshRefined =
    "The mesh of the region is written refined: an online run on it takes no "
    "\"UniformLevels\", \"SerialUniformLevels\", \"Boxes\" or \"Spheres\" refinement!";

// On a mesh of the region, the interface: the DOFs t of the region side (candidate[t])
// whose parent true DOF (parent[t], or t if parent is empty, for a space broken across
// interior boundaries) lies on the boundary of the mesh at the point of one of the saved
// interface DOFs (signatures of saved_width), where the environment was cut off.
// Collective.
std::vector<char> RegionMeshInterface(mfem::ParFiniteElementSpace &parent_fespace,
                                      const std::vector<char> &candidate,
                                      const std::vector<int> &parent,
                                      const std::vector<double> &saved_signatures,
                                      int saved_width);

// Centers for the three fields of FingerprintField, at fixed fractions of the bounding box
// of the centers of the elements with the given attributes (the origin without any).
// Collective.
std::array<std::array<double, 3>, 3> FingerprintCenters(mfem::ParMesh &mesh,
                                                        const std::vector<int> &attributes);

// The k-th (0, 1 or 2) of three fixed polynomial fields of x - center, vector fields for
// H(curl) and scalar ones for H1, as true DOFs of the space: to fingerprint operators and
// sources independently of the partition and of the DOF numbering. A center of its own for
// each field keeps their energies from being even about one common point.
mfem::Vector FingerprintField(const mfem::ParFiniteElementSpace &fespace, int k,
                              const std::array<double, 3> &center);

// Whether two fingerprints agree: their first n_counts entries exactly, the others (real,
// or complex as pairs (re, im)) each relative to the larger of the two.
bool FingerprintsMatch(std::span<const double> a, std::span<const double> b,
                       std::size_t n_counts, double tol, bool complex = false);

// Whether two sets of pairings of a source with fixed fields agree, relative to the largest
// of them, since a pairing can vanish up to a rounding that depends on the partition: real
// values, or complex as pairs (re, im).
bool PairingsMatch(std::span<const double> a, std::span<const double> b, double tol,
                   bool complex = false);

}  // namespace palace

#endif  // PALACE_FEM_SUBSTRUCTURE_HPP
