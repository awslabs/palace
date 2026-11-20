// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_POST_OPERATOR_HPP
#define PALACE_MODELS_SURFACE_POST_OPERATOR_HPP

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/output_functionals.hpp"
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

// The USER's warning (decisions 2197 (4) / 2204): the attributes of an interface
// dielectric entry of Type `type` which lie on an interior (unseparated) face of a metal
// sheet (`metal_sheet_attributes`: config::BoundaryData::GetMetalSheetAttributes) with
// vacuum on both sides, by the vacuum test of the SA coefficient (both neighbouring
// elements' maximum light speed >= 1 - 1e-6): there the SA coefficient averages the fields
// of the two sides, which read near-zero raw energy for a charged sheet. Empty unless
// `type` is SA. The faces which a crack separates are exterior and not listed, but a crack
// can leave faces of a cracked sheet interior (for example a face between the crack tip
// and a junction with an uncracked sheet), and those are listed. Collective (the same
// inputs on every rank).
std::vector<int> GetVacuumMetalSheetSubstrateAirAttributes(
    mfem::ParMesh &mesh, const MaterialOperator &mat_op, InterfaceDielectric type,
    const std::vector<int> &entry_attributes,
    const std::vector<int> &metal_sheet_attributes);

// The one-line text of that warning for the interface dielectric entry `index`. It names
// the "Crack" setting only for those of `attributes` which belong to an Impedance,
// RationalImpedance or Conductivity boundary of `boundaries` (the boundaries which have
// that setting; perfect conductors are not cracked) and which this run did not crack
// (`boundaries.cracked_attributes`, filled by the mesh read); for the cracked ones it says
// that the crack leaves these faces joined.
std::string VacuumMetalSheetSubstrateAirWarning(int index,
                                                const std::vector<int> &attributes,
                                                const config::BoundaryData &boundaries);

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
// Ruling 1715 (lane L1): a served claim's end is snapped to a tree-segment end only when
// the two POINTS differ by arithmetic noise alone (16 units in the last place of the larger
// of the coordinate and R, per coordinate: the round-off of a point formed as key0 + s / L
// (key1 - key0)); never a geometric tolerance.
inline constexpr double kClaimEndRoundoffUlps = 16.0;
// Ruling 2277: an unclaimed point counts "at a claim end" when its clamped foot lies within
// this fraction of R of an end of a piece on its segment (the residual lane's reader's
// rule, so that the kinds compare with its dump; a reporting threshold, never a
// correction).
inline constexpr double kCoverageClaimEndOverRadius = 0.05;

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
                                                      const MaterialOperator &mat_op,
                                                      bool imag_permittivity = false) const;
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

    // Marker of boundary attributes of interior boundaries separating the fields on their
    // two sides (thin metal sheets), for which the energies of both sides are summed.
    mfem::Array<int> sum_sides_marker;

    InterfaceDielectricData(const config::InterfaceDielectricData &data,
                            const mfem::Array<int> &bdr_attr_marker);

    // The entry's energy density coefficient on its attributes, zero at the quadrature
    // points outside its ownership slot and (unless apply_region is false) its region.
    std::unique_ptr<mfem::Coefficient> GetCoefficient(
        const GridFunction &E, const GridFunction *D, const MaterialOperator &mat_op,
        InterfaceDielectricComponent component = InterfaceDielectricComponent::TOTAL,
        bool apply_region = true) const;

    // Whether a quadrature point enters this entry's integrals (the ownership slot and the
    // region); the explicit sample loops apply the same filter as GetCoefficient.
    bool Includes(const mfem::Vector &point) const;
  };
  struct FarFieldData : public SurfaceData
  {
    std::vector<std::pair<double, double>> thetaphis;

    // Interior surface (not only external boundaries), which requires the libCEED path.
    bool interior = false;

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

  // libCEED surface functionals for flux and interface dielectric postprocessing,
  // constructed lazily per surface index to replace per-call coefficient evaluation and
  // boundary LinearForm assembly in the legacy paths when supported.
  mutable std::map<int, std::unique_ptr<SurfaceFunctional>> flux_funcs, eps_funcs;
  // Im{ε}-weighted electric flux functionals for complex fields in lossy dielectrics.
  mutable std::map<int, std::unique_ptr<SurfaceFunctional>> flux_funcs_imag;
  mutable std::unique_ptr<SurfaceFunctional> farfield_func;
  mutable std::vector<std::array<double, 3>> farfield_func_dirs;

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

  // Interface dielectric energies on interior boundaries with attributes in
  // sheet_attributes, which separate the fields on their two sides, sum the energies of
  // both sides. The interface dielectric "FluxRecovery" is refused on the mesh crack set of
  // the optional boundaries (the cracked internal boundaries) and where an entry's sheet
  // attributes have an interior face (cracked or not, a zero-thickness sheet has no
  // controlled normal trace of the volume projection); an exterior sheet boundary keeps it.
  SurfacePostOperator(
      const config::BoundaryPostData &postpro, ProblemType problem_type,
      const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &h1_fespace,
      mfem::ParFiniteElementSpace &nd_fespace,
      const std::vector<int> &sheet_attributes = {},
      const config::BoundaryData *boundaries = nullptr,
      std::shared_ptr<const SurfacePostGeometry> *automatic_geometry = nullptr,
      bool classified_edge_trees = false);
  SurfacePostOperator(
      const IoData &iodata, const MaterialOperator &mat_op,
      mfem::ParFiniteElementSpace &h1_fespace, mfem::ParFiniteElementSpace &nd_fespace,
      std::shared_ptr<const SurfacePostGeometry> *automatic_geometry = nullptr);

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
  // Lane L1 (ruling 1682; decisions 1675 / 1439 / 1541 on replicated memory): a
  // perimeter-portion list in a compact EXACT form, for the applied patches' raw claims
  // (the served bins of the coverage partition, replicated on every rank): each distinct
  // end point once, found by equality of the BIT PATTERNS of its three coordinates (never a
  // tolerance: the stored doubles are the appended ones, so -0.0 and +0.0 are two points,
  // each kept as given; a NaN end is refused), each distinct label once, and per portion
  // the two end-point indices, the label index and the feature. Reading portion i
  // (P0 / P1 / Label / Feature) gives back the appended values bit for bit, in the
  // appended order.
  class CompactPerimeterPortions
  {
  public:
    void Append(const std::array<double, 3> &p0, const std::array<double, 3> &p1,
                const std::string &label, int feature);
    // Releases the append lookups and trims every buffer to its size (a later Append
    // rebuilds the lookups).
    void FinishAppending();
    std::size_t Size() const { return portions.size(); }
    bool Empty() const { return portions.empty(); }
    const std::array<double, 3> &P0(std::size_t i) const { return points[portions[i].p0]; }
    const std::array<double, 3> &P1(std::size_t i) const { return points[portions[i].p1]; }
    const std::string &Label(std::size_t i) const { return labels[portions[i].label]; }
    int Feature(std::size_t i) const { return portions[i].feature; }
    std::size_t PointCount() const { return points.size(); }
    std::size_t LabelCount() const { return labels.size(); }
    // The heap bytes held: every buffer's capacity, the label strings' own storage and,
    // while appending, the lookups.
    std::size_t HeapBytes() const;
    // The portions in the plain form (the appended values, in order).
    std::vector<UncoveredPerimeterPortion> Expand() const;

  private:
    struct Portion
    {
      std::uint32_t p0 = 0, p1 = 0;
      std::int32_t label = -1;
      std::int32_t feature = -1;
    };
    std::vector<std::array<double, 3>> points;
    std::vector<std::string> labels;
    std::vector<Portion> portions;
    // Open addressing over `points` by bit pattern (kEmptySlot = free) and the label index.
    static constexpr std::uint32_t kEmptySlot = 0xFFFFFFFFu;
    std::vector<std::uint32_t> point_slots;
    std::unordered_map<std::string, std::int32_t> label_lookup;
    std::uint32_t FindOrAddPoint(const std::array<double, 3> &point);
  };
  // A read-only view of a portion list in either form; the portion readers below take it,
  // so a plain list and the compact claim list pass alike.
  class PerimeterPortionView
  {
  public:
    PerimeterPortionView(const std::vector<UncoveredPerimeterPortion> &portions)
      : plain(&portions)
    {
    }
    PerimeterPortionView(const CompactPerimeterPortions &portions) : compact(&portions) {}
    std::size_t Size() const { return plain ? plain->size() : compact->Size(); }
    bool Empty() const { return Size() == 0; }
    const std::array<double, 3> &P0(std::size_t i) const
    {
      return plain ? (*plain)[i].p0 : compact->P0(i);
    }
    const std::array<double, 3> &P1(std::size_t i) const
    {
      return plain ? (*plain)[i].p1 : compact->P1(i);
    }
    const std::string &Label(std::size_t i) const
    {
      return plain ? (*plain)[i].type : compact->Label(i);
    }
    int Feature(std::size_t i) const
    {
      return plain ? (*plain)[i].feature : compact->Feature(i);
    }

  private:
    const std::vector<UncoveredPerimeterPortion> *plain = nullptr;
    const CompactPerimeterPortions *compact = nullptr;
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
  std::map<int, std::map<std::string, double>>
  GetInterfacePartitionedEdgeEnergies(const std::set<int> &indices, const GridFunction &E,
                                      const GridFunction *D,
                                      PerimeterPortionView portions) const;
  // Diagnostics.UncoveredFeatureEnergy (lane ctx002-miss, decision 1567 (b); diagnostic
  // only): the uncovered within-R energy of GetInterfaceUncoveredEdgeEnergies split per
  // TYPE by FEATURE, keyed (type, label). Per type, the exact piece partition of that
  // type's portions (PartitionedEdgePieces, the PatchEnergy pieces) labels every piece by
  // the ascending '+'-joined ids of the features whose portions cover it; a quadrature
  // point counts for the label of the piece holding its nearest perimeter foot. ONE
  // quadrature pass per interface with GetLocalSurfaceIntegral's rule and weight. Per
  // (interface, type) the labels sum to that type's by_type entry of
  // GetInterfaceUncoveredEdgeEnergies up to the quadrature round-off (and the points whose
  // foot lies exactly on an interior piece end: a set of measure zero).
  std::map<int, std::map<std::pair<std::string, std::string>, double>>
  GetInterfaceUncoveredFeatureEdgeEnergies(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<UncoveredPerimeterPortion> &portions) const;
  // Lane L1 (decision 1634 (1), ruling 1641): THE CLOSED COVERAGE PARTITION of one field.
  // ONE quadrature pass per interface (the rule and weight of
  // GetInterfaceUncoveredFeatureEdgeEnergies); every within-R point counts in exactly ONE
  // bin, by the correction's own ownership: first the KEPT-RAW classes in the given order
  // (the corrected total keeps their raw energy whatever claim lies under them; within a
  // class the exact piece partition of PartitionedEdgePieces labels a piece by the
  // ascending '+'-joined labels of the class portions covering it: "m<feature>" tokens for
  // feature labels), then the SERVED claims (the applied patches' claims by model-set
  // label), then the `remainder_classes` (labels that only subdivide the points no kept
  // class and no served claim holds: ruling 1649 Q3), else the unclaimed remainder.
  // `within` is the total of the same pass: per interface the bins sum to it up to the
  // summation round-off.
  struct CoveragePartitionClass
  {
    std::string name;
    const std::vector<UncoveredPerimeterPortion> *portions = nullptr;
  };
  struct InterfaceCoveragePartition
  {
    double within = 0.0;
    std::map<std::pair<std::string, std::string>, double> classes;  // (class, label)
    std::map<std::string, double> served;                           // model-set label
    double unclaimed = 0.0;
    // The loop remainder's energy no served label held (MoveCoverageLoopRemainder): counted
    // in the class "LoopRemainder" and in another bin both (the bins then sum to `within`
    // plus this amount).
    double loop_remainder_unserved = 0.0;
    // Ruling 2276 (OI-COVERAGE-INSIDE-BOX-DROP-COUNTED-UNCOVERED; with dropped windows
    // only): the points a matched cluster box holds AS CONTEXT, by the label of the holding
    // boxes' model tokens ('+'-joined, ascending model index): `context` held by one box,
    // `context_overlap` by two or more boxes (a box overlap; the label joins every holder);
    // `context_beyond` the points of an inside piece that no holding box holds;
    // `inside_excluded_box` (ruling 2295 (A)) those that no serving holder holds but a
    // holder whose patch the conductor-consistency test excludes at this measurement
    // contains (such a box holds nothing). These four are bins of their own (out of
    // `unclaimed`). `kept_inside_box` = the points of the
    // dropped windows' kept class (CoverageDroppedWindows::kept_class) that lie inside a
    // matched box: a PART of that class's bins, not a bin.
    std::map<std::string, double> context, context_overlap;
    double context_beyond = 0.0, inside_excluded_box = 0.0, kept_inside_box = 0.0;
    // Ruling 2277: `unclaimed` by kind (their sum is `unclaimed` up to round-off): the
    // clamped foot within kCoverageClaimEndOverRadius R of an end of a piece of a kept
    // class, a served claim or a remainder class on its segment; no such piece on its
    // segment; elsewhere on a segment with such pieces.
    double unclaimed_claim_end = 0.0, unclaimed_no_claims = 0.0, unclaimed_elsewhere = 0.0;
    // Rulings 1714 / 1715 / 1731 (lane L1): the census of portion ends at a tree-segment
    // end, for the SERVED claims and for the KEPT classes, by DIRECTION. "Round-off short":
    // the end POINT agrees with the segment's end point within kClaimEndRoundoffUlps ulps
    // per coordinate but not exactly AND its clamped foot falls short of the segment end
    // (inside the segment; a served end is snapped to the segment end; a kept end is
    // counted and never changed: OI-FOOT-TEST-VERTEX-ROUNDOFF); its energy is EXACT: the
    // within-R raw energy whose clamped foot lies at those segment ends that the snap gives
    // back to the served claims (no unsnapped served claim reaches the end) or that no kept
    // class holds. "Round-off beyond": within that rule but its clamped foot reaches the
    // segment end exactly (beyond the end, or short by less than the foot's own round-off):
    // the foot test already holds the vertex, nothing is snapped, nothing is at stake
    // (count and largest distance only). "Near": beyond the round-off rule but within
    // kSupportEndCoincidenceOverRadius R (never snapped; its energy, the within-R raw
    // energy at those segment ends, is an upper bound). Distances in mesh units.
    struct EndCensus
    {
      int roundoff_count = 0, near_count = 0;
      double roundoff_distance = 0.0, near_distance = 0.0;
      double roundoff_energy = 0.0, near_energy = 0.0;
      int roundoff_beyond_count = 0;
      double roundoff_beyond_distance = 0.0;
    };
    EndCensus served_ends, kept_ends;
  };
  // Rulings 1714 / 1715: a portion end against a tree-segment end point. EXACT: every
  // coordinate equal; ROUNDOFF: |p_d - v_d| <= kClaimEndRoundoffUlps eps max(|v_d|, R) for
  // every d; NEAR: otherwise within kSupportEndCoincidenceOverRadius R; BEYOND: beyond
  // (not FAR: zlib's zconf.h defines FAR as an empty macro).
  enum class ClaimEndMatch : char
  {
    EXACT,
    ROUNDOFF,
    NEAR,
    BEYOND
  };
  static std::pair<ClaimEndMatch, double> ClassifyClaimEnd(const std::array<double, 3> &p,
                                                           const std::array<double, 3> &v,
                                                           double radius);
  // Ruling 2276 (OI-COVERAGE-INSIDE-BOX-DROP-COUNTED-UNCOVERED): the windows of the coupons
  // the SpatialClusterPriority rule dropped, as the coverage partition reads them. `inside`
  // = the pieces of the windows that the clip by the matched cluster boxes removed (`type`
  // = the window's label "p<patch>"); `holders` per window label the indices into `boxes`
  // of the boxes whose clip removed a piece of it; `boxes` = every matched cluster box (its
  // patch-frame box, the strict-interior test of the clip; `min` / `max` its mesh AABB, a
  // prefilter only; `token` = "m<model index>"; `patch` = its support's global patch
  // index). `kept_class` = the name of the kept class that holds the windows' parts
  // outside every box. Mesh units.
  struct CoverageContextBox
  {
    std::string token;
    int model = -1, patch = -1;
    std::array<std::array<double, 3>, 3> axes{};
    std::array<double, 3> frame_min{}, frame_max{}, min{}, max{};
  };
  struct CoverageDroppedWindows
  {
    std::vector<UncoveredPerimeterPortion> inside;
    std::map<std::string, std::vector<std::size_t>> holders;
    std::vector<CoverageContextBox> boxes;
    std::string kept_class = "SpatialClusterPriority";
    bool Empty() const { return inside.empty() && boxes.empty(); }
  };
  // A point strictly inside a box's patch-frame box (LineInsideBox's arithmetic at a
  // point).
  static bool PointInsideContextBox(const mfem::Vector &point,
                                    const CoverageContextBox &box);
  // Ruling 2295 (A): `excluded_patches` = the patches the conductor-consistency test
  // excludes at this measurement (the rule in force); a box whose patch is one of them
  // holds nothing.
  std::map<int, InterfaceCoveragePartition> GetInterfaceCoveragePartition(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<CoveragePartitionClass> &kept_classes, PerimeterPortionView served,
      const std::vector<CoveragePartitionClass> &remainder_classes = {},
      const CoverageDroppedWindows *dropped_windows = nullptr,
      const std::set<int> &excluded_patches = {}) const;
  // Ruling 2251 (OI-REGION-TARGET-CORRECTED-NOT-FILTERED): the footprint of an applied
  // patch on an entry with a Region = the entry's within-R points of the coverage partition
  // (the same rule, nearest tree segment, clamped foot and within weight) whose foot lies
  // on the patch's own claim; a portion held by several patches counts for each. The region
  // test is the entry's own 0 / 1 test at the same point (RegionCoefficient), so an `_in`
  // sum holds the raw filter's values of its points. Reduced (sum) on every rank.
  struct RegionPatchFootprint
  {
    double open_energy = 0.0, open_energy_in = 0.0;  // the points no kept class holds
    double energy = 0.0, energy_in = 0.0;            // every point of the claim
    double area = 0.0, area_in = 0.0;                // quadrature weight x within
    bool reference_in = false;                       // the patch's reference point
  };
  // Per entry of `indices` (each with a Region), per patch of the claim record (claims[i]
  // held by claim_patches[i]; reference_points by patch index): its footprint. The kept
  // classes take precedence as in GetInterfaceCoveragePartition. Collective.
  std::map<int, std::vector<RegionPatchFootprint>> GetInterfaceRegionPatchFootprints(
      const std::set<int> &indices, const GridFunction &E, const GridFunction *D,
      const std::vector<CoveragePartitionClass> &kept_classes,
      const std::vector<UncoveredPerimeterPortion> &claims,
      const std::vector<std::vector<std::size_t>> &claim_patches,
      const std::vector<std::array<double, 3>> &reference_points) const;
  // The entries of `indices` that carry a Region.
  std::set<int> GetRegionInterfaces(const std::set<int> &indices) const;
  // Ruling 1688 (lane L1, the stage-0 Missing LENGTH line): per interface the length of its
  // edge-distance tree and, per portion list, the length of the union of the list's
  // intervals on that tree (the energies' foot assignment: a portion on its midpoint's
  // segment; each segment once, overlapping intervals counted once). Mesh units; local
  // (the trees and the lists are replicated): no reduction.
  struct InterfacePerimeterLengths
  {
    double total = 0.0;
    std::vector<double> lists;
  };
  std::map<int, InterfacePerimeterLengths>
  GetInterfacePerimeterLengths(const std::set<int> &indices,
                               const std::vector<PerimeterPortionView> &lists) const;
  // The uncovered portions relabelled by feature ("m<feature>"): the Missing class's bins
  // of the coverage partition, one per feature (joined by '+' where portions overlap).
  static std::vector<UncoveredPerimeterPortion>
  CoverageMissingByFeature(const std::vector<UncoveredPerimeterPortion> &uncovered);
  // Moves a served label that holds no model claim (the PatchEnergy list's `unpatched`
  // label alone, fix round P2 2B) to the class `name` of the same partitions.
  static void MoveCoverageServedLabel(std::map<int, InterfaceCoveragePartition> &partitions,
                                      const std::string &label, const std::string &name);
  // Maxwell (OI-MX-LOOPGATE, decision 955): the loop remainders' raw energy (per interface,
  // by model row "m<idx>") is kept raw by the corrected total although it lies under the
  // served claims; it moves from the served bins to the class "LoopRemainder" by label.
  // A remainder row's energy is taken from the served labels holding the row's model token,
  // in proportion to their energy, at most their total; the rest no served label held
  // (`loop_remainder_unserved`) lies on kept-raw or remainder points already binned.
  static void
  MoveCoverageLoopRemainder(std::map<int, InterfaceCoveragePartition> &partitions,
                            const std::map<int, std::map<std::string, double>> &remainder);
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

  // P3-16: the batched amplitudes of the two forms above cannot represent the summed
  // energies of both sides of an interior sheet (SD-401 (d)), so the runs that use them
  // refuse such a sample before the solve: fails closed (collectively) when any sample of
  // the batched traversal (CacheInterfaceResponseSamples) is a summed sample
  // (InterfaceDielectricCoefficient::IsSummedSample). E is any field on the ND space (it
  // carries the mesh; its values are not read).
  void VerifyBatchedInterfaceResponse(const GridFunction &E) const;

private:
  template <InterfaceDielectric Type>
  int CountSummedInterfaceResponseSamples(const InterfaceResponseSamples &samples,
                                          const GridFunction &E) const;
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
                                 const GridFunction *D, PerimeterPortionView portions,
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
