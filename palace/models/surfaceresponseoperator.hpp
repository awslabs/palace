// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class GridFunction;
class IoData;
class LaplaceOperator;
class SpaceOperator;

// Mesh-independent automatic coupon layout. A solver retains this across AMR iterations;
// finite-element point interpolation is still rebuilt for every refined mesh.
class SurfaceResponseGeometry
{
private:
  struct Impl;
  std::shared_ptr<const Impl> impl;

  explicit SurfaceResponseGeometry(std::shared_ptr<const Impl> impl_)
    : impl(std::move(impl_))
  {
  }

  friend class SurfaceResponseOperator;
};

//
// Low-rank Schur-complement correction which replaces the local response of an ideal
// thin-metal edge by a fabrication-resolved coupon response. The correction has the form
//
//                         C = Pᵀ (S_fabricated - S_thin) P,
//
// where P evaluates the global potential relative to the metal potential at coupon
// contour knots. The assembled thin-metal Laplace operator remains the preconditioner.
//
class SurfaceResponseOperator : public Operator
{
private:
  struct MaxwellLine
  {
    int point_offset = 0;
    int point_count = 0;
  };

  struct MaxwellLineGeometry
  {
    std::array<double, 3> begin{};
    std::array<double, 3> end{};
  };

  struct MaxwellContourPath
  {
    int anchor_line = -1;
    int contour_line_offset = 0;
    int contour_line_count = 0;
    int end_line = -1;
    int conductor_line = -1;
    int closure_sign = 0;
    bool closed = true;
    std::vector<int> trace_indices;
  };

  struct MaxwellConductorPath
  {
    int line = -1;
    int trace_offset = 0;
  };

  struct OpenContourPath
  {
    std::vector<int> indices;
    int start_conductor = 0;
    int end_conductor = 0;
  };

  struct ResponseModel
  {
    int idx = 0;
    int contour_size = 0;
    int basis_size = 0;
    int conductor_state_count = 0;
    mfem::DenseMatrix fabricated_domain;
    mfem::DenseMatrix thin_domain;
    mfem::DenseMatrix domain_defect;
    mfem::DenseMatrix fixed_flux_transform;
    std::map<int, mfem::DenseMatrix> fabricated_surfaces;
    std::map<int, mfem::DenseMatrix> surface_defects;
    bool spatial_basis = false;
    std::vector<int> contour_groups;
    std::vector<int> zero_trace_indices;
    std::vector<OpenContourPath> open_contour_paths;
  };

  struct Patch
  {
    int model = -1;
    int point_offset = 0;
    int trace_offset = 0;
    double weight = 1.0;
  };

  const FiniteElementSpace &fespace;
  mfem::Array<int> dbc_tdof_list;
  std::vector<ResponseModel> models;
  std::vector<Patch> patches;
  int basis_size;
  int global_basis_size = 0;
  int global_patch_count = 0;
  int point_query_count = 0;
  std::vector<int> point_send_counts;
  std::vector<int> point_send_offsets;
  std::vector<int> point_receive_counts;
  std::vector<int> point_receive_offsets;
  std::vector<int> point_send_indices;
  std::vector<int> point_dof_offsets;
  std::vector<int> point_dofs;
  std::vector<double> point_weights;
  bool maxwell = false;

  // Maxwell postprocessing uses local coupon contour line integrals instead of H1 point
  // values. The line quadrature representation supplies both the trace action and its
  // transpose, so the same map is used for postprocessing and self-consistent correction.
  std::vector<std::vector<mfem::Vector>> maxwell_contours;
  std::vector<mfem::Vector> maxwell_anchors;
  std::vector<mfem::Vector> maxwell_secondary_anchors;
  std::vector<MaxwellLine> maxwell_lines;
  std::vector<MaxwellContourPath> maxwell_paths;
  std::vector<MaxwellConductorPath> maxwell_conductor_paths;
  std::vector<std::pair<int, int>> maxwell_patch_paths;
  int maxwell_quadrature_order = 0;

  double matching_radius = 0.0;
  double minimum_wave_speed = mfem::infinity();
  double matched_length_fraction = 1.0;
  double corner_neighborhood_fraction = 0.0;
  std::map<int, double> matched_length_fraction_by_interface;
  std::map<int, double> corner_neighborhood_fraction_by_interface;
  double maximum_curvature_ratio = 0.0;
  double maximum_library_distance = 0.0;

  mutable Vector x_free, local_x, local_y, trace, response, correction;

  void ConfigurePointCommunication(
      const mfem::Vector &xyz, int dimension,
      const std::vector<std::array<double, 3>> *weighted_tangents = nullptr);
  void ConfigureMaxwellLines(const std::vector<MaxwellLineGeometry> &line_geometry);
  void EvaluatePointValues(const Vector &x, Vector &values) const;
  void AddPointValuesTranspose(const Vector &values, Vector &y) const;
  void EvaluatePoints(const Vector &x, Vector &values) const;
  void EvaluateMaxwellLines(const Vector &x, Vector &values) const;
  void AddMaxwellLinesTranspose(const Vector &values, Vector &y) const;
  void BuildMaxwellTrace(const Vector &line_values, Vector &values) const;
  void BuildMaxwellTraceTranspose(const Vector &values, Vector &line_values) const;
  void ApplyTrace(const Vector &x, Vector &values) const;
  void ApplyTraceTranspose(const Vector &values, Vector &y) const;
  void ApplyUneliminated(const Vector &x, Vector &y) const;

public:
  struct EnergyCorrection
  {
    double domain = 0.0;
    std::map<int, double> interfaces;
  };

  struct ElectrostaticResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;
    std::map<int, double> trace_closure_spread;
    double maximum_trace_closure_spread = 0.0;
  };

  struct MaxwellResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;

    double kR = 0.0;
    double loop_residual = 0.0;
    double response_weighted_loop_residual = 0.0;
    double loop_response_failure_fraction = 0.0;
    double matched_length_fraction = 1.0;
    double corner_neighborhood_fraction = 0.0;
    std::map<int, double> matched_length_fraction_by_interface;
    std::map<int, double> corner_neighborhood_fraction_by_interface;
    double maximum_curvature_ratio = 0.0;
    double maximum_library_distance = 0.0;
    double maximum_trace_closure_spread = 0.0;
    bool confident = true;
  };

  SurfaceResponseOperator(
      const IoData &iodata, const LaplaceOperator &laplace_op,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry = nullptr);
  SurfaceResponseOperator(
      const IoData &iodata, const SpaceOperator &space_op,
      std::shared_ptr<const SurfaceResponseGeometry> *automatic_geometry = nullptr);

  void Mult(const Vector &x, Vector &y) const override;
  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }

  // Add the contribution from prescribed essential values to an already assembled
  // thin-metal right-hand side.
  void EliminateRHS(const Vector &x, Vector &rhs) const;

  // Evaluate the nondimensional domain- and surface-energy defects for a global field.
  EnergyCorrection GetEnergyCorrection(const Vector &x) const;

  // Evaluate the complete fabricated-coupon surface energy for every mapped target
  // interface. Corrected participation replaces the measured global core with this data.
  std::map<int, double> GetFabricatedSurfaceEnergy(const Vector &x) const;

  // Evaluate fixed-trace and fixed-flux coupon responses on an unchanged electrostatic
  // thin-metal potential. This is the electrostatic analogue of the postprocessing-only
  // Maxwell correction.
  ElectrostaticResponse GetElectrostaticResponse(const Vector &x,
                                                 bool include_fixed_flux = true) const;

  // Evaluate a postprocessing-only response for a complex Nedelec Maxwell field. Coupon
  // voltages are reconstructed from transverse contour integrals and applied through
  // Hermitian response-matrix quadratic forms.
  MaxwellResponse GetMaxwellResponse(const GridFunction &E,
                                     std::complex<double> omega) const;

  bool HasSurfaceResponse() const;
  std::set<int> GetTargetInterfaces() const;

  int GetBasisSize() const { return global_basis_size; }
  int GetPatchCount() const { return global_patch_count; }
  int GetEdgeCount() const { return GetPatchCount(); }
  double GetPatchWeight() const;
  double GetMatchingRadius() const { return matching_radius; }
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
