// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"

// XX TODO: Rename BoundaryPostOperator for config file consistency?

namespace palace
{

class IoData;
class MaterialOperator;

namespace config
{

struct InterfaceDielectricData;
struct CapacitanceData;
struct InductanceData;

}  // namespace config

//
// A class handling boundary surface postprocessing.
//
class SurfacePostOperator
{
private:
  // Mapping from surface index to data structure containing surface postprocessing
  // information for surface loss, charge, or magnetic flux.
  struct SurfaceData
  {
    std::vector<mfem::Array<int>> attr_lists;

    virtual ~SurfaceData() = default;

    virtual std::unique_ptr<mfem::Coefficient>
    GetCoefficient(std::size_t i, const mfem::ParGridFunction &U,
                   const MaterialOperator &mat_op) const = 0;
  };
  struct InterfaceDielectricData : public SurfaceData
  {
    DielectricInterfaceType type;
    double epsilon, ts, tandelta;
    std::vector<mfem::Vector> sides;

    InterfaceDielectricData(const config::InterfaceDielectricData &data,
                            const mfem::ParMesh &mesh);

    std::unique_ptr<mfem::Coefficient>
    GetCoefficient(std::size_t i, const mfem::ParGridFunction &U,
                   const MaterialOperator &mat_op) const override;
  };
  struct SurfaceChargeData : public SurfaceData
  {
    SurfaceChargeData(const config::CapacitanceData &data, const mfem::ParMesh &mesh);

    std::unique_ptr<mfem::Coefficient>
    GetCoefficient(std::size_t i, const mfem::ParGridFunction &U,
                   const MaterialOperator &mat_op) const override;
  };
  struct SurfaceFluxData : public SurfaceData
  {
    mfem::Vector direction;

    SurfaceFluxData(const config::InductanceData &data, const mfem::ParMesh &mesh);

    std::unique_ptr<mfem::Coefficient>
    GetCoefficient(std::size_t i, const mfem::ParGridFunction &U,
                   const MaterialOperator &mat_op) const override;
  };
  std::map<int, InterfaceDielectricData> eps_surfs;
  std::map<int, SurfaceChargeData> charge_surfs;
  std::map<int, SurfaceFluxData> flux_surfs;

  // Reference to material property operator (not owned).
  const MaterialOperator &mat_op;

  // Unit function used for computing surface integrals.
  mutable mfem::GridFunction ones;

  double GetLocalSurfaceIntegral(const SurfaceData &data,
                                 const mfem::ParGridFunction &U) const;

public:
  SurfacePostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &h1_fespace);

  // Access data structures for postprocessing the surface with the given type.
  const auto &GetEps() const { return eps_surfs; }
  const auto &GetCap() const { return charge_surfs; }
  const auto &GetInd() const { return flux_surfs; }

  // Get surface integrals computing dielectric interface energy, surface charge, or
  // surface magnetic flux.
  double GetInterfaceLossTangent(int idx) const;
  double GetInterfaceElectricFieldEnergy(int idx,
                                         const mfem::ParComplexGridFunction &E) const;
  double GetInterfaceElectricFieldEnergy(int idx, const mfem::ParGridFunction &E) const;
  double GetSurfaceElectricCharge(int idx, const mfem::ParComplexGridFunction &E) const;
  double GetSurfaceElectricCharge(int idx, const mfem::ParGridFunction &E) const;
  double GetSurfaceMagneticFlux(int idx, const mfem::ParComplexGridFunction &B) const;
  double GetSurfaceMagneticFlux(int idx, const mfem::ParGridFunction &B) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
