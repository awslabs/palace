// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/output_functionals.hpp"

namespace palace
{

class GridFunction;
class IoData;
class MaterialOperator;

namespace config
{

struct BoundaryPostData;
struct SurfaceFluxData;
struct InterfaceDielectricData;
struct FarFieldPostData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling boundary surface postprocessing.
//
class SurfacePostOperator
{
private:
  // Mapping from surface index to data structure containing surface postprocessing
  // information for surface flux or interface dielectric participation.
  struct SurfaceData
  {
    mfem::Array<int> attr_list;

    virtual ~SurfaceData() = default;
  };
  struct SurfaceFluxData : public SurfaceData
  {
    SurfaceFlux type;
    bool two_sided;
    mfem::Vector center;

    SurfaceFluxData(const config::SurfaceFluxData &data, const mfem::ParMesh &mesh,
                    const mfem::Array<int> &bdr_attr_marker);

    std::unique_ptr<mfem::Coefficient> GetCoefficient(const mfem::ParGridFunction *E,
                                                      const mfem::ParGridFunction *B,
                                                      const MaterialOperator &mat_op,
                                                      bool imag_permittivity = false) const;
  };
  struct InterfaceDielectricData : public SurfaceData
  {
    InterfaceDielectric type;
    double t, epsilon, tandelta;

    // Marker of boundary attributes of interior boundaries separating the fields on their
    // two sides (thin metal sheets), for which the energies of both sides are summed.
    mfem::Array<int> sum_sides_marker;

    InterfaceDielectricData(const config::InterfaceDielectricData &data,
                            const mfem::ParMesh &mesh,
                            const mfem::Array<int> &bdr_attr_marker);

    std::unique_ptr<mfem::Coefficient> GetCoefficient(const GridFunction &E,
                                                      const MaterialOperator &mat_op) const;
  };
  struct FarFieldData : public SurfaceData
  {
    std::vector<std::pair<double, double>> thetaphis;

    FarFieldData() = default;
    FarFieldData(const config::FarFieldPostData &data, const mfem::ParMesh &mesh,
                 const mfem::Array<int> &bdr_attr_marker);

    std::size_t size() const { return thetaphis.size(); }
  };

  // Reference to material property operator (not owned).
  const MaterialOperator &mat_op;

  // Reference to scalar finite element space used for computing surface integrals (not
  // owned).
  mfem::ParFiniteElementSpace &h1_fespace;

  // Reference to vector finite element space used for computing far-field integrals (not
  // owned).
  mfem::ParFiniteElementSpace &nd_fespace;

  double GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                 const mfem::Array<int> &attr_marker) const;

  // Linear map from the dofs of one element attached to a boundary element of an
  // interface to the interface energy amplitudes at the quadrature points of that boundary
  // element: a[offset, offset + A.Height()) += A x_e, with x_e the element dofs of the
  // field (on a face neighbor element when ghost) and the element dof transformation folded
  // into A. See GetInterfaceElectricFieldEnergyMatrices.
  struct InterfaceAmplitudeMap
  {
    bool ghost;
    int offset;
    mfem::Array<int> vdofs;
    mfem::DenseMatrix A;
  };
  std::vector<InterfaceAmplitudeMap>
  GetInterfaceAmplitudeMaps(const InterfaceDielectricData &data,
                            const std::vector<int> &bdr_elems, int &num_amplitudes) const;

  // libCEED surface functionals for flux and interface dielectric postprocessing,
  // constructed lazily per surface index to replace per-call coefficient evaluation and
  // boundary LinearForm assembly in the legacy paths when supported.
  mutable std::map<int, std::unique_ptr<SurfaceFunctional>> flux_funcs, eps_funcs;
  // Im{ε}-weighted electric flux functionals for complex fields in lossy dielectrics.
  mutable std::map<int, std::unique_ptr<SurfaceFunctional>> flux_funcs_imag;
  mutable std::unique_ptr<SurfaceFunctional> farfield_func;
  mutable std::vector<std::array<double, 3>> farfield_func_dirs;

public:
  // Data structures for postprocessing the surface with the given type.
  std::map<int, SurfaceFluxData> flux_surfs;
  std::map<int, InterfaceDielectricData> eps_surfs;
  FarFieldData farfield;

  // Interface dielectric energies on interior boundaries with attributes in
  // sheet_attributes, which separate the fields on their two sides, sum the energies of
  // both sides.
  SurfacePostOperator(const config::BoundaryPostData &postpro, ProblemType problem_type,
                      const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &h1_fespace,
                      mfem::ParFiniteElementSpace &nd_fespace,
                      const std::vector<int> &sheet_attributes = {});
  SurfacePostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &h1_fespace,
                      mfem::ParFiniteElementSpace &nd_fespace);

  // Get surface integrals computing electric or magnetic field flux, or the stationary real
  // power flux, through a boundary. For complex fields the electric flux uses the complex
  // permittivity D = (Re{ε} + i Im{ε}) E of a lossy dielectric.
  std::complex<double> GetSurfaceFlux(int idx, const GridFunction *E,
                                      const GridFunction *B) const;

  // Batch version for multiple theta/phi pairs
  std::vector<std::array<std::complex<double>, 3>>
  GetFarFieldrE(const std::vector<std::pair<double, double>> &theta_phi_pairs,
                const GridFunction &E, const GridFunction &B,
                std::complex<double> omega) const;

  // Get surface integrals computing interface dielectric energy.
  double GetInterfaceLossTangent(int idx) const;
  double GetInterfaceElectricFieldEnergy(int idx, const GridFunction &E) const;

  // Interface dielectric energy forms over real fields e_0, ..., e_{n-1}: for each
  // interface index, the n × n symmetric positive semidefinite matrix B with B_kl = U(e_k,
  // e_l), where U is the symmetric bilinear form of the interface energy, U(e, e) =
  // GetInterfaceElectricFieldEnergy(idx, e) with the same quadrature points, evaluation
  // sides, and coefficients. For complex coefficients z, zᴴ B z is the interface energy of
  // the complex field Σ_k z_k e_k (the energies of its real and imaginary parts add). The
  // callback field(k) returns the k-th field as a grid function on the H(curl) space of
  // this operator, whose face neighbor data is exchanged here. Each field is evaluated into
  // its amplitudes a(e) at the quadrature points of the interfaces, with U(e, e) = |a(e)|²,
  // so that B = Aᵀ A with A = [a(e_0), ...]. The amplitudes of all fields are stored at
  // once, within memory_budget bytes per process: larger interfaces are evaluated in
  // several passes, each of which calls field(k) for k = 0, ..., n - 1 in order.
  // Collective.
  std::map<int, Eigen::MatrixXd> GetInterfaceElectricFieldEnergyMatrices(
      int num_fields, const std::function<mfem::ParGridFunction &(int)> &field,
      std::size_t memory_budget = std::size_t(1) << 27) const;

  int GetVDim() const { return mat_op.SpaceDimension(); };
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
