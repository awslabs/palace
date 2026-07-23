// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <map>
#include <memory>
#include <unordered_set>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"

namespace palace
{

class GridFunction;
class IoData;
class MaterialOperator;
class EdgeDistanceTree;

namespace config
{

struct BoundaryPostData;
struct BoundaryData;
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
                                                      const MaterialOperator &mat_op) const;
  };
  struct InterfaceDielectricData : public SurfaceData
  {
    InterfaceDielectric type;
    double t, epsilon, tandelta;
    bool flux_recovery;
    std::vector<double> edge_distances;
    double edge_distance_smoothing;
    bool localize_edge_energy;
    std::array<double, 3> edge_frame_normal;
    std::shared_ptr<const EdgeDistanceTree> edge_distance_tree;

    InterfaceDielectricData(const config::InterfaceDielectricData &data,
                            const mfem::Array<int> &bdr_attr_marker);

    std::unique_ptr<mfem::Coefficient> GetCoefficient(
        const GridFunction &E, const GridFunction *D, const MaterialOperator &mat_op,
        InterfaceDielectricComponent component = InterfaceDielectricComponent::TOTAL) const;
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

  struct LocalVolumeEdgeEnergyCache
  {
    const EdgeDistanceTree *edge_distance_tree;
    std::vector<double> edge_distances;
    double edge_distance_smoothing;
    std::array<double, 3> edge_frame_normal;
    std::vector<double> energies;
    std::vector<double> vertex_energies;
    std::array<std::vector<double>, 6> polarized_energies;
  };
  mutable std::vector<LocalVolumeEdgeEnergyCache> local_volume_edge_energy_cache;

  double GetLocalSurfaceIntegral(mfem::Coefficient &f,
                                 const mfem::Array<int> &attr_marker) const;
  const LocalVolumeEdgeEnergyCache &
  GetLocalVolumeEdgeElectricFieldEnergies(const InterfaceDielectricData &data,
                                          const GridFunction &E) const;

public:
  struct InterfaceEdgeEnergy
  {
    double distance;
    double energy_outside;
    double energy_annulus;
  };
  struct InterfaceLocalEdgeEnergy
  {
    int edge;
    bool automatic;
    int component;
    int chain;
    std::array<int, 2> vertex_types;
    std::array<double, 3> process_normal;
    std::array<double, 3> p0;
    std::array<double, 3> p1;
    double length;
    double distance;
    double energy_total;
    double energy_inside;
    double energy_annulus;
    double energy_vertex_inside;
    std::array<double, 2> energy_total_polarized;
    std::array<double, 2> energy_inside_polarized;
    std::array<double, 2> energy_annulus_polarized;
    double energy_volume_annulus;
    double energy_volume_vertex_annulus;
    std::array<double, 6> energy_volume_annulus_polarized;
  };

  // Data structures for postprocessing the surface with the given type.
  std::map<int, SurfaceFluxData> flux_surfs;
  std::map<int, InterfaceDielectricData> eps_surfs;
  FarFieldData farfield;

  SurfacePostOperator(const config::BoundaryPostData &postpro, ProblemType problem_type,
                      const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &h1_fespace,
                      mfem::ParFiniteElementSpace &nd_fespace,
                      const std::unordered_set<int> *cracked_attributes = nullptr,
                      const config::BoundaryData *boundaries = nullptr);
  SurfacePostOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &h1_fespace,
                      mfem::ParFiniteElementSpace &nd_fespace);

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
  double GetInterfaceElectricFieldEnergy(int idx, const GridFunction &E,
                                         const GridFunction *D = nullptr) const;
  std::vector<InterfaceEdgeEnergy>
  GetInterfaceEdgeElectricFieldEnergies(int idx, const GridFunction &E,
                                        const GridFunction *D = nullptr) const;
  InterfaceEdgeEnergy
  GetInterfaceOuterElectricFieldEnergy(int idx, const GridFunction &E,
                                       const GridFunction *D = nullptr) const;
  std::vector<InterfaceLocalEdgeEnergy>
  GetInterfaceLocalEdgeElectricFieldEnergies(int idx, const GridFunction &E,
                                             const GridFunction *D = nullptr,
                                             bool include_volume = true) const;

  std::size_t GetNInterfaceEdgeEntries() const;
  std::size_t GetNInterfaceLocalEdgeEntries() const;
  bool NeedsFluxRecovery() const;
  void ResetInterfaceLocalEdgeEnergyCache() const;

  int GetVDim() const { return mat_op.SpaceDimension(); };
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
