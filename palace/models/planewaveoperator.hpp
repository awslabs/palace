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
// only evaluated in that layer of elements, which must be vacuum. The spaces, the mesh, and
// the operator are unchanged.
//
class PlaneWaveOperator
{
public:
  struct PlaneWaveSource
  {
    mfem::Vector direction, polarization, origin;  // Nondimensional
    double amplitude;                              // Nondimensional
    mfem::Array<int> attr_marker;
  };

private:
  std::map<int, PlaneWaveSource> sources;

  // Cached true-dof marker (1 on the dofs of the closure of the total-field region, 0
  // elsewhere) for each source.
  mutable std::map<int, Vector> tf_markers;

  const Vector &GetTotalFieldMarker(int idx, const FiniteElementSpace &nd_fespace) const;

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
};

}  // namespace palace

#endif  // PALACE_MODELS_PLANE_WAVE_OPERATOR_HPP
