// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <array>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "utils/edgedistance.hpp"

namespace palace
{

class GridFunction;
class IoData;
class MaterialOperator;
class InterfaceOwnershipPartition;
class InterfaceRegion;
class OwnershipQuadrature;

namespace config
{

struct BoundaryPostData;
struct BoundaryData;
struct SurfaceFluxData;
struct InterfaceDielectricData;
struct FarFieldPostData;
struct SolverData;

}  // namespace config

enum class ProblemType : char;

// Whether the run's surface-response correction (electrostatic ResponseCorrection or
// Maxwell SurfaceResponseCorrection) requests the exclusion-energy export
// (Diagnostics.ExclusionEnergy; idfix DESIGN section 8): the SurfacePostOperator builds the
// CLASSIFIED edge trees the export partitions against only then.
bool ExclusionEnergyRequested(const config::SolverData &solver);

// Mesh-independent automatic edge-distance geometry retained across AMR iterations. Only
// aggregate edge-distance trees are cached; mesh-segment-localized diagnostics are rebuilt.
class SurfacePostGeometry
{
private:
  std::map<int, std::shared_ptr<const EdgeDistanceTree>> automatic_edge_distance_trees;
  // The exclusion-energy export's CLASSIFIED trees (every metal edge segment of the
  // interface, physical or not; idfix DESIGN section 8, decision 1076), beside.
  std::map<int, std::shared_ptr<const EdgeDistanceTree>> classified_edge_distance_trees;

  friend class SurfacePostOperator;
};

//
// A class handling boundary surface postprocessing.
//
class SurfacePostOperator
{
public:
  struct InterfaceQuadratureRule
  {
    int interface_index, geometry, requested_order, actual_order, points;
    long long local_faces;
    double minimum_weight;
    bool ownership_rule;
  };

  struct InterfaceResponseMatrix
  {
    double distance;
    mfem::DenseMatrix energy_total;
    mfem::DenseMatrix energy_inside;
    mfem::DenseMatrix energy_total_normal;
    mfem::DenseMatrix energy_inside_normal;
    mfem::DenseMatrix energy_total_tangential;
    mfem::DenseMatrix energy_inside_tangential;
  };

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
    bool save_local_edge_energy;
    std::array<double, 3> edge_frame_normal;
    std::shared_ptr<const InterfaceOwnershipPartition> ownership;
    std::shared_ptr<const OwnershipQuadrature> ownership_rule;
    int ownership_slot;
    int ownership_quadrature_order;
    // The optional quadrature-level spatial filter (decision 352 follow-up (2)).
    std::shared_ptr<const InterfaceRegion> region;
    std::shared_ptr<const EdgeDistanceTree> edge_distance_tree;
    // The exclusion-energy export's tree (decision 1076): EVERY metal edge segment of the
    // interface the identification classifies (physical, fold, non-manifold, truncation,
    // port), so that the field near a cut or fold is attributed to that edge and not, as
    // in `edge_distance_tree` and every energy table reading it, to the nearest PHYSICAL
    // edge. Automatic-edge interfaces only (null otherwise).
    std::shared_ptr<const EdgeDistanceTree> classified_edge_distance_tree;
    // The EdgeExcludeSegments exclusion's coverage record of this entry's perimeter
    // (decision 657 (b): the listed length and the perimeter length the coincidence test
    // removed); nullopt for an entry without EdgeExcludeSegments.
    std::optional<EdgeExclusionRecord> edge_exclusion;

    InterfaceDielectricData(const config::InterfaceDielectricData &data,
                            const mfem::Array<int> &bdr_attr_marker);

    // The entry's energy density coefficient on its attributes, zero at the quadrature
    // points outside its ownership slot and its region.
    std::unique_ptr<mfem::Coefficient> GetCoefficient(
        const GridFunction &E, const GridFunction *D, const MaterialOperator &mat_op,
        InterfaceDielectricComponent component = InterfaceDielectricComponent::TOTAL) const;

    // Whether a quadrature point enters this entry's integrals (the ownership slot and the
    // region); the explicit sample loops apply the same filter as GetCoefficient.
    bool Includes(const mfem::Vector &point) const;
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

  // Axisymmetric (r, z) mesh: every MFEM-side surface and volume quadrature weight carries
  // the revolution measure 2 pi x (the libCEED operators carry it in their geometry data).
  const bool axisymmetric;

  // The exclusion-energy export's CLASSIFIED edge tree (idfix DESIGN section 8; decision
  // 1076) is built per automatic-edge interface only when Diagnostics.ExclusionEnergy is
  // requested (the d1073 review's MINOR-5: an unrequested tree is memory with no reader).
  const bool classified_edge_trees;

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

  double GetLocalSurfaceIntegral(mfem::Coefficient &f, const mfem::Array<int> &attr_marker,
                                 int quadrature_order = 0) const;
  const LocalVolumeEdgeEnergyCache &
  GetLocalVolumeEdgeElectricFieldEnergies(const InterfaceDielectricData &data,
                                          const GridFunction &E) const;

  template <InterfaceDielectric Type>
  std::vector<InterfaceResponseMatrix> GetInterfaceElectricFieldEnergyMatricesImpl(
      const InterfaceDielectricData &data, const std::vector<const GridFunction *> &E,
      const std::vector<const GridFunction *> &D, int interface_index, int quadrature_extra,
      std::vector<InterfaceQuadratureRule> *quadrature_rules) const;

public:
  // The rank-local quadrature samples of one interface response matrix (the streaming
  // one-pass Gram of the archive reduction): the same boundary elements, quadrature
  // rule, points, physical weights, edge distances and ownership selection the batched
  // GetInterfaceElectricFieldEnergyMatrices traversal visits, cached once so that every
  // basis field is evaluated once into its amplitude row (Count() samples x components
  // doubles: the normal amplitude, then the tangential ones, per sample).
  struct InterfaceResponseSamples
  {
    int interface_index;
    const InterfaceDielectricData *data;
    int components;
    std::vector<int> elements;
    std::vector<const mfem::IntegrationRule *> rules;
    std::vector<int> points;
    std::vector<double> weights;
    std::vector<double> distances;

    std::size_t Count() const { return weights.size(); }
    std::size_t RowSize() const { return Count() * static_cast<std::size_t>(components); }
  };

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

  SurfacePostOperator(
      const config::BoundaryPostData &postpro, ProblemType problem_type,
      const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &h1_fespace,
      mfem::ParFiniteElementSpace &nd_fespace,
      const std::unordered_set<int> *cracked_attributes = nullptr,
      const config::BoundaryData *boundaries = nullptr,
      std::shared_ptr<const SurfacePostGeometry> *automatic_geometry = nullptr,
      bool classified_edge_trees = false);
  SurfacePostOperator(
      const IoData &iodata, const MaterialOperator &mat_op,
      mfem::ParFiniteElementSpace &h1_fespace, mfem::ParFiniteElementSpace &nd_fespace,
      std::shared_ptr<const SurfacePostGeometry> *automatic_geometry = nullptr);

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
  std::map<int, InterfaceEdgeEnergy>
  GetInterfaceOuterElectricFieldEnergies(const std::set<int> &indices,
                                         const GridFunction &E,
                                         const GridFunction *D = nullptr) const;

  // Uncovered requirements (decision 394 F2): a sub-segment of the metal perimeter (mesh
  // coordinates) claimed by a feature the response library has no model for, with the
  // feature's type and id.
  struct UncoveredPerimeterPortion
  {
    std::array<double, 3> p0{}, p1{};
    std::string type;
    int feature = -1;
  };
  // The interface energy within the largest configured edge distance R whose nearest
  // perimeter point lies on an uncovered portion (the same edge-distance window weight as
  // the outside energy, so that outside + uncovered + the modelled within-R partition the
  // interface), in total and per feature type (the total is the sum of the types).
  struct UncoveredEdgeEnergy
  {
    double energy = 0.0;
    std::map<std::string, double> by_type;
  };
  // Per requested interface (each with EdgeDistances): the quadrature point's nearest
  // perimeter segment comes from the interface's own edge-distance tree and its foot
  // parameter is tested against the closed intervals of the uncovered portions found on
  // that segment (a portion within kSignatureParameterToleranceOverRadius x R of a segment
  // of this perimeter; a portion of another plane finds no segment and is skipped), an
  // exact partition without a distance tolerance: a foot at a shared vertex belongs to
  // both incident portions of one feature (a corner's two arms), an interior boundary has
  // measure zero.
  std::map<int, UncoveredEdgeEnergy> GetInterfaceUncoveredEdgeEnergies(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<UncoveredPerimeterPortion> &portions) const;
  // PatchEnergy (decision 922 B-1; the 851 (C) export): the within-R energy of each
  // requested interface partitioned EXACTLY over the portions' labels (`type` = a
  // '+'-joined "m<idx>" model-set label): per label (the union label where portions of
  // different labels overlap on a tree segment) the energy whose nearest perimeter foot
  // lies on a piece of that label, plus the key "unclaimed" for the within-R energy whose
  // foot lies on no portion. The labels and "unclaimed" sum to the within-R energy (the
  // interface energy minus the outside energy) exactly, up to the quadrature of the two
  // integrals.
  // The exclusion-energy export (idfix DESIGN section 8; decisions 1074 / 1076): the same
  // exact partition read against the CLASSIFIED tree (every metal edge segment of the
  // interface), the portions' `type` = the exclusion Class; the remainder key is "Claimed"
  // (the within-R energy whose nearest classified edge lies on no exclusion portion: the
  // feature-claimed perimeter, modelled or uncovered). The classes and "Claimed" sum to the
  // within-R energy of the interface under this tree's distance.
  std::map<int, std::map<std::string, double>> GetInterfaceClassifiedEdgeEnergies(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<UncoveredPerimeterPortion> &portions) const;
  std::map<int, std::map<std::string, double>> GetInterfacePartitionedEdgeEnergies(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<UncoveredPerimeterPortion> &portions) const;
  // OI-MX-LOOPGATE (decision 955): the within-R energy of each requested interface on the
  // portions, WEIGHTED per portion (`weights`, parallel to `portions`: the patch's rho_p)
  // and grouped by the portions' `type` label (the model row "m<idx>"): per label the
  // integral of rho x the within-R energy density over the pieces whose nearest perimeter
  // foot lies on a portion of that label. A piece covered by several portions takes the
  // LARGEST weight among them and is attributed to that portion's label (the design's
  // co-located rule). One surface integral per (interface, label).
  std::map<int, std::map<std::string, double>> GetInterfaceWeightedPortionEdgeEnergies(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<UncoveredPerimeterPortion> &portions,
      const std::vector<double> &weights) const;
  std::vector<InterfaceLocalEdgeEnergy>
  GetInterfaceLocalEdgeElectricFieldEnergies(int idx, const GridFunction &E,
                                             const GridFunction *D = nullptr,
                                             bool include_volume = true) const;

  // Assemble aggregate interface response matrices for many real basis fields in one
  // quadrature pass per interface. This avoids reevaluating the surface postprocessor for
  // every pair of basis fields. Diagnostic extra order affects only this local surface
  // rule selection (never global assembly defaults). Optional rule usage is rank-local;
  // diagnostic rules must have strictly positive weights. Ownership overrides are rejected.
  std::map<int, std::vector<InterfaceResponseMatrix>>
  GetInterfaceElectricFieldEnergyMatrices(
      const std::vector<const GridFunction *> &E,
      const std::vector<const GridFunction *> &D = {}, int quadrature_extra = 0,
      std::vector<InterfaceQuadratureRule> *quadrature_rules = nullptr) const;

  // Streaming one-pass assembly of the same matrices (the archive reduction): cache the
  // samples of every localized interface once (the default quadrature rule of the batched
  // traversal, no extra order), evaluate each basis field once into its amplitude row
  // (already scaled by sqrt(0.5 x physical weight), so that the Gram of the rows is the
  // energy matrix), then assemble S = F W F^T per interface with a cache-blocked symmetric
  // rank-k update - the total matrices and the inside-radius ones (W = the edge-distance
  // window per sample) - and reduce over `comm`. The quadrature rule, order, weights,
  // sample set and ownership selection are those of the batched traversal.
  std::vector<InterfaceResponseSamples> CacheInterfaceResponseSamples() const;
  void EvaluateInterfaceResponseRow(const InterfaceResponseSamples &samples,
                                    const GridFunction &E, const GridFunction *D,
                                    double *row) const;

private:
  template <InterfaceDielectric Type>
  void EvaluateInterfaceResponseRowImpl(const InterfaceResponseSamples &samples,
                                        const GridFunction &E, const GridFunction *D,
                                        double *row) const;
  // The exact partition of the within-R energy over the portions' labels against the tree
  // `classified` selects (the physical tree or the classified one), the remainder keyed
  // `remainder_label`: GetInterfacePartitionedEdgeEnergies and
  // GetInterfaceClassifiedEdgeEnergies.
  std::map<int, std::map<std::string, double>>
  PartitionInterfaceEdgeEnergies(const std::set<int> &indices, const GridFunction &E,
                                 const GridFunction *D,
                                 const std::vector<UncoveredPerimeterPortion> &portions,
                                 bool classified, const std::string &remainder_label) const;

public:
  std::vector<InterfaceResponseMatrix>
  AssembleInterfaceResponseMatrices(const InterfaceResponseSamples &samples,
                                    const double *rows, int basis_size,
                                    MPI_Comm comm) const;

  std::size_t GetNInterfaceEdgeEntries() const;
  std::size_t GetNInterfaceLocalEdgeEntries() const;
  bool NeedsFluxRecovery() const;
  void ResetInterfaceLocalEdgeEnergyCache() const;

  int GetVDim() const { return mat_op.SpaceDimension(); };
};

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
