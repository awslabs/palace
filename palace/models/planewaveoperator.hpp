// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
#define PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP

#include <map>
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
    mfem::Array<int> attr_marker;
  };

private:
  std::map<int, PlaneWaveSource> sources;

  // Cached true-dof marker (1 on the dofs of the closure of the total-field region, 0
  // elsewhere) for each source.
  mutable std::map<int, Vector> tf_markers;

  // Marker of the local (vector) dofs of the closure of the total-field region, combined
  // across the processes sharing each dof.
  void GetLocalTotalFieldMarker(int idx, const mfem::ParFiniteElementSpace &fes,
                                mfem::Array<int> &lmark) const;

public:
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
  // total-field region which touch its boundary, where the incident field is evaluated.
  // Aborts if these are not one homogeneous, isotropic and lossless material.
  void SetUpBackground(const MaterialOperator &mat_op,
                       const FiniteElementSpace &nd_fespace);

  const PlaneWaveSource &GetSource(int idx) const { return sources.at(idx); }

  // True-dof marker of the plane wave with the given index: 1 on the dofs of the closure
  // of the total-field region, 0 elsewhere.
  const Vector &GetTotalFieldMarker(int idx, const FiniteElementSpace &nd_fespace) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
