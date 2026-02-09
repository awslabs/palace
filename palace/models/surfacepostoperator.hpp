// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"

namespace palace
{

class GridFunction;
class IoData;
class MaterialOperator;
class Mesh;

namespace config
{

struct SurfaceFluxData;
struct InterfaceDielectricData;
struct FarFieldPostData;

}  // namespace config

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

    SurfaceFluxData(const config::SurfaceFluxData &data, const Mesh &mesh,
                    const mfem::Array<int> &bdr_attr_marker);

    std::unique_ptr<mfem::Coefficient> GetCoefficient(const Mesh &mesh,
                                                      const mfem::ParGridFunction *E,
                                                      const mfem::ParGridFunction *B,
                                                      const MaterialOperator &mat_op) const;
  };
  struct InterfaceDielectricData : public SurfaceData
  {
    InterfaceDielectric type;
    double t, epsilon, tandelta;

    InterfaceDielectricData(const config::InterfaceDielectricData &data, const Mesh &mesh,
                            const mfem::Array<int> &bdr_attr_marker);

    std::unique_ptr<mfem::Coefficient> GetCoefficient(const Mesh &mesh,
                                                      const GridFunction &E,
                                                      const MaterialOperator &mat_op) const;
  };
  struct FarFieldData : public SurfaceData
  {
    std::vector<std::pair<double, double>> thetaphis;

    FarFieldData() = default;
    FarFieldData(const config::FarFieldPostData &data, const Mesh &mesh,
                 const mfem::Array<int> &bdr_attr_marker);

    size_t size() const { return thetaphis.size(); }
  };

  // Reference to material property operator (not owned).
  const MaterialOperator &mat_op;

  // Reference to mesh (not owned).
  const Mesh &mesh;

  // Reference to vector finite element space used for computing far-field integrals (not
  // owned).
  FiniteElementSpace &nd_fespace;

public:
  // Data structures for postprocessing the surface with the given type.
  std::map<int, SurfaceFluxData> flux_surfs;
  std::map<int, InterfaceDielectricData> eps_surfs;
  FarFieldData farfield;

  SurfacePostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      const Mesh &mesh, FiniteElementSpace &nd_fespace);

  // Get surface integrals computing electric or magnetic field flux through a boundary.
  std::complex<double> GetSurfaceFlux(int idx, const GridFunction *E,
                                      const GridFunction *B) const;

  // Batch version for multiple theta/phi pairs
  std::vector<std::array<std::complex<double>, 3>>
  GetFarFieldrE(const std::vector<std::pair<double, double>> &theta_phi_pairs,
                const GridFunction &E, const GridFunction &B, double omega_re,
                double omega_im) const;

  // Get surface integrals computing interface dielectric energy.
  double GetInterfaceLossTangent(int idx) const;
  double GetInterfaceElectricFieldEnergy(int idx, const GridFunction &E) const;

  int GetVDim() const { return mat_op.SpaceDimension(); };
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
