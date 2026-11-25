// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
#define PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP

#include <map>
#include <set>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class Mesh;
class MaterialOperator;
class Units;

namespace config
{

struct BoundaryData;
struct PeriodicBoundaryData;
struct PlaneWaveData;

}  // namespace config

//
// Incident plane wave excitation with the total-field/scattered-field (TF/SF) formulation.
// The unknown is the total field E on the true dofs of the closure of the total-field
// region T (dofs on its boundary S included) and the scattered field E - E_inc on all other
// dofs O. With the system operator A of the total field (A E = 0 without other sources),
// substituting E = x + E_inc on O in the rows of T and x = E - E_inc on T in the rows of O
// gives
//
//   A x = b,  b_T = -(A g_O)_T,  b_O = (A g_T)_O,
//
// where g is the interpolant of the incident field and g_T, g_O its restrictions to T and
// O. Only the elements outside of T which touch S couple T and O dofs, so b is the
// discrete counterpart of the equivalent (Love) currents on S, and the incident field is
// only evaluated in that layer of elements. On nonconforming meshes the coupling follows
// the prolongation P: the layer takes the elements whose local dofs depend on both T and
// O true dofs, through the hanging-entity constraints. These must be one homogeneous,
// isotropic, lossless and non-dispersive material (the background), whose speed of light
// sets the wavenumber of the incident field. The spaces, the mesh, and the operator are
// unchanged.
//
class PlaneWaveOperator
{
public:
  struct PlaneWaveSource
  {
    mfem::Vector direction, polarization, origin;  // Nondimensional
    double amplitude;                              // Nondimensional
    double c = 1.0;  // Speed of light of the background (nondimensional), vacuum by default
    // Floquet wave vector of the solution (nondimensional), empty without one. When
    // kf_scaled, it is k_F / omega and scales with the frequency.
    mfem::Vector kf;
    bool kf_scaled = false;
    mfem::Array<int> attr_marker;
  };

private:
  std::map<int, PlaneWaveSource> sources;

  // Nondimensional length of one mesh length unit, for the messages.
  double mesh_scale = 1.0;

  // Frequency in GHz of a unit nondimensional angular frequency, for the messages.
  double unit_GHz = 1.0;

  // Whether the setup is printed.
  bool print = true;

  // Whether the domain has periodic boundaries (see CheckPeriodicity).
  bool periodic_cell = false;

  // Whether the quasi-periodicity of the incident field is checked on the interpolated
  // field at every frequency, when some periodic boundary pair has no explicit translation.
  bool check_interpolated = false;

  // The plane waves already warned of the jump of their interpolated field.
  mutable std::set<int> warned_interpolated;

  // Cached true-dof marker (1 on the dofs of the closure of the total-field region, 0
  // elsewhere) for each source.
  mutable std::map<int, Vector> tf_markers;

  // The incident field of the last plane wave and frequency, shared by the right-hand side
  // and the reconstruction of the total and scattered fields of the same solve.
  mutable struct
  {
    int idx = -1;
    double omega = 0.0;
    ComplexVector g;
  } last_incident;
  const ComplexVector &GetLastIncidentField(int idx, double omega,
                                            const FiniteElementSpace &nd_fespace) const;

  // Marker of the local dofs of the elements which couple the closure of the total-field
  // region and the rest of the domain (see GetLayer), and optionally these elements outside
  // of the region.
  void GetCoupling(int idx, const FiniteElementSpace &nd_fespace, mfem::Array<int> &l_c,
                   mfem::Array<int> *coupling_elements = nullptr) const;

  // Marker of the local (vector) dofs of the closure of the total-field region, combined
  // across the processes sharing each dof.
  void GetLocalTotalFieldMarker(int idx, const mfem::ParFiniteElementSpace &fes,
                                mfem::Array<int> &lmark) const;

public:
  // Interpolate the incident field of the plane wave with the given excitation index at the
  // angular frequency omega into the true dofs of the given space, as the periodic envelope
  // E exp(+i k_F . x) when the solution has a Floquet wave vector k_F. Aborts if the
  // interpolated field is not continuous across the periodic boundaries, which happens
  // when it is not quasi-periodic with k_F. Returns false if there is no such plane wave.
  bool GetIncidentField(int idx, double omega, const FiniteElementSpace &nd_fespace,
                        ComplexVector &g) const;

  // The same, only on the true dofs of the given local elements (zero elsewhere), without
  // the check of the quasi-periodicity. Every true dof of these elements is set on its
  // owning process if one of its local elements is in the list.
  void GetIncidentField(int idx, double omega, const FiniteElementSpace &nd_fespace,
                        const mfem::Array<int> &elements, ComplexVector &g) const;

  // The local true dofs which enter the right-hand side of the plane wave with the given
  // index, which the system matrix couples between the closure of the total-field region
  // and the rest of the domain (on a nonconforming mesh, also through the constraints of
  // hanging entities), and the local elements on which the incident field is interpolated
  // for them. Optionally, the local elements outside of the total-field region which carry
  // this coupling, the layer of elements around the region.
  void GetLayer(int idx, const FiniteElementSpace &nd_fespace, mfem::Array<int> &tdofs,
                mfem::Array<int> &elements,
                mfem::Array<int> *coupling_elements = nullptr) const;

  // The local elements outside of the total-field region of the plane wave with the given
  // index, and the local boundary elements, which have a vertex in the closure of the
  // region (on a nonconforming mesh, also a hanging vertex constrained by it), across the
  // processes. With inside = true, the local elements of the region which have a vertex in
  // the closure of the rest of the domain instead (and the boundary elements likewise).
  void GetVertexNeighbors(int idx, mfem::ParMesh &mesh, mfem::Array<int> &elements,
                          mfem::Array<int> &bdr_elements, bool inside = false) const;

  // For a nonconforming refinement of the given mesh, add to the marked local elements the
  // whole layer of elements on both sides of the boundary of the total-field region of
  // each plane wave (the elements with a vertex on it, and the coupling elements) when part
  // of it is marked, so that the boundary stays conforming.
  void AddLayerToMarkedElements(Mesh &mesh, mfem::Array<int> &marked_elements) const;

  PlaneWaveOperator(const std::map<int, config::PlaneWaveData> &plane_wave,
                    const Units &units, const mfem::ParMesh &mesh, bool print = true);

  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto Size() const { return sources.size(); }
  bool Empty() const { return sources.empty(); }

  // Add the TF/SF right-hand side of the plane wave with the given excitation index, if
  // any, for the system operator A at the angular frequency omega (nondimensional). Returns
  // true if a plane wave contributed.
  bool AddExcitationVector(int idx, double omega, const ComplexOperator &A,
                           const FiniteElementSpace &nd_fespace, ComplexVector &RHS) const;

  // The total field E_t = x + (1 - M) g and the scattered field E_s = x - M g everywhere,
  // from the solution x of the plane wave with the given excitation index at the angular
  // frequency omega (total field on the dofs of the closure of the total-field region,
  // M = 1, and scattered field elsewhere), with the incident field g as for the
  // right-hand side. Returns false, and leaves E_t and E_s unchanged, if there is no such
  // plane wave.
  bool GetTotalAndScatteredField(int idx, double omega,
                                 const FiniteElementSpace &nd_fespace,
                                 const ComplexVector &x, ComplexVector &E_t,
                                 ComplexVector &E_s) const;

  // Set the background material of each plane wave from the elements outside of its
  // total-field region which touch its boundary, where the incident field is evaluated,
  // and the Floquet wave vector of the solution, if any. Aborts if these elements are not
  // one homogeneous, isotropic, lossless and non-dispersive material, or if a PML region
  // overlaps the total-field region or touches its closure.
  void SetUpBackground(const MaterialOperator &mat_op,
                       const FiniteElementSpace &nd_fespace);

  // Check that the incident field is quasi-periodic with the Floquet wave vector of the
  // solution (zero without one) across the periodic boundary pairs, at all the given
  // angular frequencies (nondimensional): (k khat - k_F) . a must be a multiple of 2 pi for
  // the translation a of each pair, and zero when k_F scales with the frequency. Aborts
  // on a phase mismatch above 1e-2 and warns above 1e-6. Pairs with a transformation
  // which is not a translation are rejected, and without an explicit translation the
  // check is done on the interpolated field at each frequency instead.
  void CheckPeriodicity(const config::PeriodicBoundaryData &periodic,
                        const std::vector<double> &omega, mfem::ParMesh &mesh);

  // Check that no boundary with a condition touches the layer of elements around the
  // total-field region of any plane wave, where the incident field does not satisfy it,
  // except PEC or PMC symmetry planes along the direction of propagation which it
  // satisfies, and that no absorbing boundary or Floquet port, which truncates the
  // scattered field, touches the total-field region anywhere.
  void CheckBoundaries(const config::BoundaryData &boundaries,
                       const FiniteElementSpace &nd_fespace) const;

  // Whether the domain is a periodic cell, with periodic boundaries.
  bool IsPeriodicCell() const { return periodic_cell; }

  const PlaneWaveSource &GetSource(int idx) const { return sources.at(idx); }

  // True-dof marker of the plane wave with the given index: 1 on the dofs of the closure
  // of the total-field region, 0 elsewhere.
  const Vector &GetTotalFieldMarker(int idx, const FiniteElementSpace &nd_fespace) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
