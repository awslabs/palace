// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
#define PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP

#include <map>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class MaterialOperator;
class Units;

namespace config
{

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
// only evaluated in that layer of elements. These must be one homogeneous, isotropic and
// lossless material (the background), whose speed of light sets the wavenumber of the
// incident field. The spaces, the mesh, and the operator are unchanged.
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

  // Whether the quasi-periodicity of the incident field is checked on the interpolated
  // field at every frequency, when some periodic boundary pair has no explicit translation.
  bool check_interpolated = false;

  // Cached true-dof marker (1 on the dofs of the closure of the total-field region, 0
  // elsewhere) for each source.
  mutable std::map<int, Vector> tf_markers;

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

  PlaneWaveOperator(const std::map<int, config::PlaneWaveData> &plane_wave,
                    const Units &units, const mfem::ParMesh &mesh);

  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto Size() const { return sources.size(); }
  bool Empty() const { return sources.empty(); }

  // Add the TF/SF right-hand side of the plane wave with the given excitation index, if
  // any, for the system operator A at the angular frequency omega (nondimensional). Returns
  // true if a plane wave contributed.
  bool AddExcitationVector(int idx, double omega, const ComplexOperator &A,
                           const FiniteElementSpace &nd_fespace, ComplexVector &RHS) const;

  // Whether the boundary elements with the given attributes (marker) carry any dof of the
  // closure of the total-field region of the plane wave with the given index, on any
  // process. There the solution is the total field, not the scattered field.
  bool TouchesTotalField(int idx, const FiniteElementSpace &nd_fespace,
                         const mfem::Array<int> &bdr_attr_marker) const;

  // Set the background material of each plane wave from the elements outside of its
  // total-field region which touch its boundary, where the incident field is evaluated,
  // and the Floquet wave vector of the solution, if any. Aborts if these elements are not
  // one homogeneous, isotropic and lossless material.
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
                        const std::vector<double> &omega);

  const PlaneWaveSource &GetSource(int idx) const { return sources.at(idx); }

  // True-dof marker of the plane wave with the given index: 1 on the dofs of the closure
  // of the total-field region, 0 elsewhere.
  const Vector &GetTotalFieldMarker(int idx, const FiniteElementSpace &nd_fespace) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
