// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <set>
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
  struct PointEvaluation
  {
    int element = -1;
    mfem::IntegrationPoint point;
    bool local = false;
  };

  struct ResponseModel
  {
    int idx = 0;
    int basis_size = 0;
    mfem::DenseMatrix fabricated_domain;
    mfem::DenseMatrix thin_domain;
    mfem::DenseMatrix domain_defect;
    mfem::DenseMatrix fixed_flux_transform;
    std::map<int, mfem::DenseMatrix> fabricated_surfaces;
    std::map<int, mfem::DenseMatrix> surface_defects;
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
  std::vector<PointEvaluation> points;

  // Maxwell postprocessing uses local coupon contour line integrals instead of H1 point
  // values. Contours are stored in the same order as patches.
  std::vector<std::vector<mfem::Vector>> maxwell_contours;
  std::vector<mfem::Vector> maxwell_anchors;
  mutable std::unique_ptr<mfem::FindPointsGSLIB> maxwell_finder;
  int maxwell_quadrature_order = 0;

  double matching_radius = 0.0;
  double minimum_wave_speed = mfem::infinity();
  double matched_length_fraction = 1.0;
  double corner_neighborhood_fraction = 0.0;
  double maximum_curvature_ratio = 0.0;
  double maximum_library_distance = 0.0;

  mutable Vector x_free, local_x, local_y, trace, response, correction;
  mutable mfem::Array<int> element_dofs;
  mutable Vector shape, element_values;

  void EvaluatePoints(const Vector &x, Vector &values) const;
  void AddPointTranspose(int point, double value, Vector &y) const;
  void ApplyTrace(const Vector &x, Vector &values) const;
  void ApplyTraceTranspose(const Vector &values, Vector &y) const;
  void ApplyUneliminated(const Vector &x, Vector &y) const;

public:
  struct EnergyCorrection
  {
    double domain = 0.0;
    std::map<int, double> interfaces;
  };

  struct MaxwellResponse
  {
    double domain_correction = 0.0;
    double domain_correction_fixed_flux = 0.0;
    std::map<int, double> fabricated_surface_energy;
    std::map<int, double> fabricated_surface_energy_fixed_flux;

    double kR = 0.0;
    double loop_residual = 0.0;
    double matched_length_fraction = 1.0;
    double corner_neighborhood_fraction = 0.0;
    double maximum_curvature_ratio = 0.0;
    double maximum_library_distance = 0.0;
    double maximum_trace_closure_spread = 0.0;
    bool confident = true;
  };

  SurfaceResponseOperator(const IoData &iodata, const LaplaceOperator &laplace_op);
  SurfaceResponseOperator(const IoData &iodata, const SpaceOperator &space_op);

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

  // Evaluate a postprocessing-only response for a complex Nedelec Maxwell field. Coupon
  // voltages are reconstructed from transverse contour integrals and applied through
  // Hermitian response-matrix quadratic forms.
  MaxwellResponse GetMaxwellResponse(const GridFunction &E,
                                     std::complex<double> omega) const;

  bool HasSurfaceResponse() const;
  std::set<int> GetTargetInterfaces() const;

  int GetBasisSize() const { return basis_size; }
  int GetPatchCount() const { return static_cast<int>(patches.size()); }
  int GetEdgeCount() const { return GetPatchCount(); }
  double GetPatchWeight() const;
  double GetMatchingRadius() const { return matching_radius; }
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_RESPONSE_OPERATOR_HPP
