// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_POST_OPERATOR_HPP
#define PALACE_MODELS_POST_OPERATOR_HPP

#include <algorithm>
#include <complex>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/boundary_derived_field_bundle.hpp"
#include "fem/boundary_physical_trace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "fem/lineintegrator.hpp"
#include "fem/point_field_evaluator.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/domainpostoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/surfacepostoperator.hpp"
#include "models/surfaceresponsefilmrule.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/ceedparaviewdatacollection.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/units.hpp"

namespace palace
{

class CurlCurlOperator;
class RomOperator;
class ErrorIndicator;
class FaceSamplingPlan;
class IoData;
class LaplaceOperator;
class MaterialOperator;
class SpaceOperator;
class SurfaceCurrentOperator;
class WavePortOperator;

namespace config
{

struct BoundaryData;
struct DomainData;
struct ProblemData;
struct SolverData;

}  // namespace config

// D6 (USER 1221; decision 1229): the first-order fixed-trace frequency shift of the raw
// mode from a postprocessed Maxwell response, d omega / omega = (U_m,defect - U_e,defect) /
// (2 U_raw) with U_e,defect = response.domain_correction (the electric fixed-trace domain
// correction, lowering the frequency when positive), U_m,defect =
// response.magnetic_domain_correction (the magnetic one, raising it) and U_raw the raw
// normalization energy (1/2 E^H M E plus the lumped-port capacitor energy); NaN unless the
// response carries the magnetic term and passes magnetic_confident (the magnetic term's
// own loop gate). Decision 1250 (mx-sc-magnetic Phase 2): the shift follows the SAME
// coverage rule as the electric fixed-trace correction it extends, which the binary prints
// on partial coverage with the uncovered requirements kept raw (decision 394 F2); it is NOT
// gated by closure_independent_confident (whose coverage criterion >= 1 - 1e-10 withheld
// it on every device with ends) nor by the trace-closure spread: those flags are printed
// beside it, and `fixed-trace shift coverage` (FixedTraceShiftCoverage) names the matched
// energy share the shift rests on. One function, so that the column of record and its
// unit test read the same formula (review d1227 MAJOR-2).
double FixedTraceFrequencyShift(const SurfaceResponseOperator::MaxwellResponse &response,
                                double raw_normalization_energy);

// Decision 1250: the matched energy share used by the fixed-trace shift = 1 - (the raw
// within-R interface energy the electric correction keeps RAW: the uncovered requirements,
// the DomainBoundary claims, the spatial-cluster priority and conductor-consistency
// exclusions and the loop remainder) / (the raw within-R energy of the target interfaces);
// 1 at full coverage, NaN when there is no within-R energy to cover.
double FixedTraceShiftCoverage(double within_radius_raw_energy, double kept_raw_energy);

// The recovered electric flux D of the interface "FluxRecovery" postprocessing lives on the
// flux recovery space of the error estimators (BuildBrokenFluxSpace: broken across the
// interior faces of the sheet attributes whenever the mesh has any, SD-403); a vector of
// another size is refused (GetRecoveredElectricFluxTrueVSize).
inline constexpr const char *kRecoveredFluxSpaceMessage =
    "The recovered electric flux does not match the flux recovery space of the error "
    "estimators (broken across the interior boundary sheet faces when the mesh has any)!";

// Statically specify if solver uses real or complex fields.

template <ProblemType solver_t>
consteval bool HasComplexGridFunction()
{
  return solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
         solver_t == ProblemType::BOUNDARYMODE;
}

// Statically specify what fields a solver uses

template <ProblemType solver_t>
consteval bool HasVGridFunction()
{
  return solver_t == ProblemType::ELECTROSTATIC;
}

template <ProblemType solver_t>
consteval bool HasAGridFunction()
{
  return solver_t == ProblemType::MAGNETOSTATIC;
}

template <ProblemType solver_t>
consteval bool HasEGridFunction()
{
  return solver_t != ProblemType::MAGNETOSTATIC;
}

template <ProblemType solver_t>
consteval bool HasBGridFunction()
{
  return solver_t != ProblemType::ELECTROSTATIC;
}

// Scale a grid function after redimensionalizing its mesh by L, undoing the effect of the
// mesh scaling on the field values mapped from the reference element through the element
// Jacobians. The scaling is L^p, where p depends on the map type of the finite element
// space: 0 for H1, 1 for H(curl), dim - 1 for H(div), and dim for L2 with INTEGRAL map
// type, with dim the dimension of the mesh. Boundary output evaluates fields in the
// neighboring domain elements, so the same scaling applies to the domain and boundary data
// collections. Wave port fields need no rescaling, since their submesh remains
// nondimensionalized.
void ScaleGridFunction(double L, GridFunction &u);

// Apply ScaleGridFunction to each of the given grid functions which is not null.
void ScaleGridFunctions(double L, std::initializer_list<GridFunction *> fields);

// Scale gridfunctions from non-dimensional to SI units.
void DimensionalizeGridFunctions(Units &units, std::unique_ptr<GridFunction> &E,
                                 std::unique_ptr<GridFunction> &B,
                                 std::unique_ptr<GridFunction> &V,
                                 std::unique_ptr<GridFunction> &A);

// Scale gridfunctions from SI units to non-dimensional.
void NondimensionalizeGridFunctions(Units &units, std::unique_ptr<GridFunction> &E,
                                    std::unique_ptr<GridFunction> &B,
                                    std::unique_ptr<GridFunction> &V,
                                    std::unique_ptr<GridFunction> &A);

//
// A class to handle solution postprocessing for all solvers.
//
template <ProblemType solver_t>
class PostOperator
{
protected:
  // Pointer to operator handling discretization and FEM space appropriate to solver. It
  // also contains the reference to all domains, boundary conditions, etc. needed for
  // measurement and printing.
  fem_op_t<solver_t> *fem_op;

  // Unit converter from IOData to scale mesh and measurements. Lightweight class so it is
  // cheap to copy, rather than keep another reference to IOData.
  Units units;

  // Base post-op output directory.
  fs::path post_dir;

  // P3-14 (SD-3202): with FluxRecovery on a mesh with interior sheet faces the recovered D
  // lives on the estimators' flux space broken across them; flux_space is that view, built
  // by the same rule from the same RT space and sheet attributes (null when unbroken: D is
  // then held on the RT space, as before). D is held for the energies per element (each
  // element reads the version of its own side) on a discontinuous RT space of the same
  // elements, flux_sides_space (ghost sides by the face-neighbour exchange); flux_first
  // maps a broken D to its FIRST version on the RT space (the output fields). Declared
  // before the grid functions that live on these spaces (destroyed after them).
  std::unique_ptr<FiniteElementSpace> flux_space;
  std::unique_ptr<mfem::FiniteElementCollection> flux_sides_fec;
  std::unique_ptr<FiniteElementSpace> flux_sides_space;
  std::unique_ptr<mfem::HypreParMatrix> flux_first;

  // A grid function holding a recovered flux for the energies: on the RT space, or per side
  // on flux_sides_space when the flux space is broken; and its setter from the true DOFs of
  // the flux space (MFEM_VERIFY on the size).
  std::unique_ptr<GridFunction> MakeRecoveredFluxGridFunction(bool complex);
  void SetRecoveredFluxGridFunction(const Vector &d, mfem::ParGridFunction &flux) const;

  // Refuse a recovered flux whose true-dof size is not the flux recovery space's, on EVERY
  // process (review r2152 A S4-5): the sizes are rank-local, so the mismatch is reduced
  // over the communicator first (a process whose sizes agree, e.g. one without a sheet DOF,
  // would otherwise go on alone into the face-neighbour exchange). Collective.
  void VerifyRecoveredFluxTrueVSize(int size) const;

  // Fields: Electric, recovered electric displacement, Magnetic, Scalar Potential, Vector
  // Potential.
  std::unique_ptr<GridFunction> E, D_recovered, B, V, A;

  // Whether SetRecoveredElectricFlux ran since the last transient measurement: the
  // transient MeasureAndPrintAll takes no D, so it fails closed without a D for its time
  // step (SD-504).
  bool recovered_flux_set = false;

  // Mode analysis: normal (out-of-plane) E component on H1 space, and in-plane B field
  // on ND space for visualization. The in-plane B is the dominant component:
  //   Bt = -(kn/omega)(z_hat x Et) + (1/(i*omega))(grad_t(Ez) x z_hat)
  std::unique_ptr<GridFunction> En, Bt_inplane;

  // Field output format control flags.
  bool enable_paraview_output = false;
  bool enable_gridfunction_output = false;
  bool field_coefficients_initialized = false;

  // How many / which fields to output.
  int output_delta_post = 0;                          // printing rate (TRANSIENT)
  int output_n_post = 0;                              // max printing (OTHER SOLVERS)
  std::vector<std::size_t> output_save_indices = {};  // explicit saves

  // Whether any output formats were specified.
  bool AnyOutputFormats() const
  {
    return enable_paraview_output || enable_gridfunction_output;
  }
  bool AnythingToSave() const
  {
    return (output_delta_post > 0) || (output_n_post > 0) || !output_save_indices.empty();
  }

  // Whether any fields should be written at all.
  bool ShouldWriteFields() const { return AnyOutputFormats() && AnythingToSave(); }

  // Whether any fields should be written for this step.
  bool ShouldWriteFields(int step) const
  {
    return AnyOutputFormats() &&
           ((output_delta_post > 0 && step % output_delta_post == 0) ||
            (output_n_post > 0 && step < output_n_post) ||
            std::binary_search(output_save_indices.cbegin(), output_save_indices.cend(),
                               step));
  }

  // Whether fields should be written for a particular output format (at a given step).
  bool ShouldWriteParaviewFields() const
  {
    return enable_paraview_output && AnythingToSave();
  }
  bool ShouldWriteParaviewFields(int step) const
  {
    return enable_paraview_output && ShouldWriteFields(step);
  }
  bool ShouldWriteGridFunctionFields() const
  {
    return enable_gridfunction_output && AnythingToSave();
  }
  bool ShouldWriteGridFunctionFields(int step) const
  {
    return enable_gridfunction_output && ShouldWriteFields(step);
  }

  // ParaView data collection: writing fields to disk for visualization.
  // This is an optional, since ParaViewDataCollection has no default (empty) ctor,
  // and we only want initialize it if ShouldWriteParaviewFields() returns true.
  std::optional<CeedParaViewDataCollection> paraview;
  std::optional<CeedParaViewDataCollection> paraview_bdr;

  // Largest cycle written, so the final save never overwrites a step written out of order.
  int paraview_max_cycle = -1;

  // MFEM grid function output details.
  std::string gridfunction_output_dir;
  static constexpr std::size_t pad_digits_default = 6;

  // Measurements of field solution for ParaView files (full domain or surfaces).

  // Poynting Coefficient (vector for 3D), Electric Boundary Field (re+im), Magnetic
  // Boundary Field (re+im), Surface Current (re+im).
  std::unique_ptr<mfem::VectorCoefficient> S, E_sr, E_si, B_sr, B_si, J_sr, J_si;

  // Electric Energy Density, Magnetic Energy Density, Surface Charge (re+im).
  std::unique_ptr<mfem::Coefficient> U_e, U_m, Q_sr, Q_si;

  // libCEED evaluators and optional output grid functions for the domain coefficient
  // fields (U_e, U_m, S), replacing per-point host coefficient evaluation at ParaView
  // save time when supported. ParaView output uses lazy component-major point buffers;
  // gridfunction output still uses the GridFunction path when requested.
  std::unique_ptr<mfem::L2_FECollection> viz_fec;
  std::unique_ptr<mfem::ParFiniteElementSpace> viz_scalar_fespace, viz_vector_fespace;
  std::unique_ptr<mfem::ParGridFunction> U_e_gf, U_m_gf, S_gf, Sn_gf;
  std::unique_ptr<PointFieldEvaluator> E_domain_eval, B_domain_eval, En_domain_eval,
      Bt_domain_eval, V_domain_eval, A_domain_eval, U_e_eval, U_m_eval, S_eval, Sn_eval;

  // libCEED evaluators for boundary collection fields (E_s, B_s, Q_s, J_s, U_e, U_m,
  // S). CeedParaViewDataCollection evaluates these lazily into one temporary buffer per
  // field in the same integer boundary-element/refined-point order used for writing,
  // avoiding coefficient adapters or floating-point point lookup at save time.
  std::unique_ptr<PointFieldEvaluator> E_bdr_eval, B_bdr_eval, V_bdr_eval, A_bdr_eval,
      Q_bdr_eval, J_bdr_eval, Ue_bdr_eval, Um_bdr_eval, S_bdr_eval;
  // Im{ε}-weighted surface charge evaluator for complex fields in lossy dielectrics
  // (Q_s = (Re{ε} + i Im{ε}) E ⋅ n), used when the fused boundary bundle is unavailable.
  std::unique_ptr<PointFieldEvaluator> Q_bdr_eval_imag;
  // One immutable boundary lattice/trace plan is shared by every compatible ParaView
  // evaluator built for this PostOperator setup.
  std::shared_ptr<const FaceSamplingPlan> bdr_sampling_plan;
  // Save-generation-scoped physical ND/RT/scalar traces shared by every compatible
  // boundary point evaluator. The cache keeps canonical side traces device-resident;
  // individual fields retain their legacy output routing.
  std::shared_ptr<BoundaryPhysicalTraceCache> bdr_trace_cache;
  // 3D E/B boundary family bundle. Compatible direct-VTU fields register non-owning
  // packed slices from this result; unrelated fields retain individual evaluators.
  std::shared_ptr<BoundaryDerivedFieldBundle> bdr_derived_bundle;

  // Wave port boundary mode field postprocessing.
  struct WavePortFieldData
  {
    std::unique_ptr<mfem::VectorCoefficient> E0r, E0i;
  };
  std::map<int, WavePortFieldData> port_E0;

  // Setup coefficients/evaluators for field output after the first solution is ready.
  // Construction-only FE staging is detached and freed immediately; assembled libCEED
  // operators remain cached for later writes in a frequency or mode sweep.
  void SetupFieldCoefficients();
  void EnsureFieldCoefficientsSetup();

  // Initialize Paraview, register all fields to write.
  void InitializeParaviewDataCollection(const fs::path &sub_folder_name = "");
  void EnsureParaviewDataCollection();

  // Subfolder for the field output of a driven excitation, the same for ParaView and grid
  // function output: none for a single excitation.
  fs::path ExcitationSubFolder(int ex_idx) const
    requires(solver_t == ProblemType::DRIVEN);

public:
  // Public overload for the driven solver only, that takes in an excitation index and
  // sets the correct sub_folder_name path for the primary function above.
  void InitializeParaviewDataCollection(int ex_idx)
    requires(solver_t == ProblemType::DRIVEN);

  // Whether this solve requests any field output. Drivers use this to release solver
  // memory before constructing the field-output evaluators only when necessary, and the
  // adaptive online sweep uses it to choose an ordering that keeps frequency-dependent port
  // state hot across excitations.
  bool WillWriteFields() const { return ShouldWriteFields(); }
  bool WillWriteFields(int step) const { return ShouldWriteFields(step); }

protected:
  // Write to disk the E- and B-fields extracted from the solution vectors. Note that
  // fields are not redimensionalized, to do so one needs to compute: B <= B * (μ₀ H₀), E
  // <= E * (Z₀ H₀), V <= V * (Z₀ H₀ L₀), etc.
  void WriteParaviewFields(double time, int step);
  void WriteParaviewFieldsFinal(const ErrorIndicator *indicator = nullptr);
  void WriteMFEMGridFunctions(double time, int step, const fs::path &sub_folder_name = "");
  void WriteMFEMGridFunctionsFinal(const ErrorIndicator *indicator = nullptr);

  // CSV Measure & Print.

  // PostOperatorCSV<solver_t> is a class that contains csv tables and printers of
  // measurements. Conceptually, its members could be a part of this class, like the
  // ParaView fields and functions above. It has been separated out for code readability. To
  // achieve this, it is has a pointer back to its "parent" PostOperator class and is a
  // friend class so it can access the private measurement_cache and references of the
  // system from fem_op.
  friend PostOperatorCSV<solver_t>;

  PostOperatorCSV<solver_t> post_op_csv;

  // Helper classes that actually do some measurements that will be saved to csv files.

  DomainPostOperator dom_post_op;  // Energy in bulk

  // L-vector of a magnetostatic solution on a broken space (empty otherwise), for the
  // domain energies.
  Vector a_l;
  SurfacePostOperator surf_post_op;  // Dielectric Interface Energy, Flux, and FarField
  mutable InterpolationOperator interp_op;  // E & B fields: mutates during measure

  mutable Measurement measurement_cache;

  // Optional Maxwell surface response correction. The ordinary measurement cache,
  // surface-Q.csv, and primary E/D fields remain the raw thin-metal result.
  std::unique_ptr<SurfaceResponseOperator> surface_response_op;
  // The impedance-metal rule (Solver.SurfaceResponseCorrection.ImpedanceMetalRule "Record";
  // decisions 2006 / 2012 / 2024; EXPERIMENTAL): the film resolution; the CANDIDATE
  // operator ALWAYS built beside the operator of record in Record (a twin of it when the
  // rule acts on no feature; the law reading runs on it only); its kept-raw portion lists;
  // the candidate tables and the palace.json block. NEVER OF RECORD: no value or file of
  // record reads them. All empty / null under "Off".
  std::unique_ptr<FilmRuleResolution> impedance_film_rule;
  std::unique_ptr<SurfaceResponseOperator> impedance_candidate_op;
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> impedance_candidate_uncovered,
      impedance_candidate_domain_boundary, impedance_candidate_spatial_cluster_priority;
  mutable std::unique_ptr<TableWithCSVFile> impedance_candidate_table, impedance_law_table;
  mutable nlohmann::json impedance_candidate_record;
  // The coverage status of record per interface of the last printed measurement (root;
  // Record only; read by the candidate rows, ruling 2065 (2)).
  mutable std::map<int, std::string> impedance_record_coverage_status;
  mutable double impedance_candidate_measure_seconds = 0.0;
  // Whether the rule acts on any candidate feature (false = a twin: the self-check).
  bool impedance_candidate_acts = false;
  std::shared_ptr<const SurfaceResponseGeometry> *impedance_candidate_geometry = nullptr;
  // RECORD: the film rule, the printed lines and the candidate operator (always built),
  // with the setup cost of both operators measured.
  void SetUpImpedanceRuleCandidate(const IoData &iodata, const SpaceOperator &space_op,
                                   double record_setup_seconds, long record_setup_bytes);
  // RECORD, after the measurement of record is printed: the candidate's gate, response and
  // kept-raw energies, the law reading, the self-check (a twin), the candidate rows.
  void MeasureAndPrintImpedanceRuleCandidate(double output_index, int excitation,
                                             int source, const ComplexVector &e) const;
  // APPLY (decision 2206 (a)): the per-feature law service of the operator of record
  // (surface-response-impedance-rule-features.csv under the diagnostics directory, or its
  // "impedance-rule" sub-directory when unset: never at top level) and its count line.
  void WriteImpedanceApplyFeatures() const;
  // The uncovered requirements (decision 394 F2) and the DomainBoundary raw claims (F-DB-a)
  // of the response operator as the surface postoperator's perimeter portions: their raw
  // within-R interface energy is kept in the corrected interface energies as the
  // electrostatic driver keeps it (decision 751; the DomainBoundary portions are empty on
  // the Maxwell path today, which runs no containment pass).
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion> surface_response_uncovered;
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
      surface_response_domain_boundary;
  // The SpatialClusterPriority raw remainders (OI-SCP-HOLE, decision 922): the dropped
  // coupons' windows outside every cluster box, kept raw like the two lists above.
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
      surface_response_spatial_cluster_priority;
  // PatchEnergy (decision 922 B-1; the 851 (C) export): the applied patches' raw claims by
  // model set (empty without PatchEnergy), whose within-R raw energy is partitioned per
  // model set (surface-response-model-raw-energy.csv).
  SurfacePostOperator::CompactPerimeterPortions surface_response_model_raw;
  // OI-COVERAGE-EXCLUSION-ENERGY (idfix DESIGN section 8; decisions 1064 (2) / 1074): the
  // excluded / unclaimed perimeter portions by class (empty without
  // Diagnostics.ExclusionEnergy), whose within-R raw energy is partitioned per class
  // (surface-response-exclusion-energy.csv). Read fresh at every measurement: the
  // MaxwellPatchExcluded part is set at the construction, the rest at the placement.
  std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
      surface_response_exclusion_energy;
  bool surface_response_exclusion_energy_enabled = false;
  bool surface_response_within_r_volume_energy_enabled = false;  // COL-3 (R-MAG-15)
  // Diagnostics.UncoveredFeatureEnergy (decision 1567 (b)): off by default (every output
  // then byte-identical to a run without the option).
  bool surface_response_uncovered_feature_energy_enabled = false;
  bool surface_response_postprocess_enabled = false;
  bool surface_response_self_consistent_enabled = false;
  // PatchEnergy on the Maxwell path (decisions 913 / 920 / 922 B): the per-patch export
  // (surface-response-patches.csv at construction; surface-response-patch-energy.csv,
  // -traces.csv and -patch-loops.csv per print) and the per-model raw partition. Off by
  // default: every output is then byte-identical to a run without the option.
  bool surface_response_patch_energy_enabled = false;
  std::unique_ptr<GridFunction> surface_response_corrected_field;
  std::unique_ptr<GridFunction> surface_response_corrected_flux;
  // P3-14 (SD-3202): on a broken flux space the corrected flux is held twice: the output
  // field above (the FIRST version on the RT space, upstream's rule for a broken field) and
  // the per-side field for the energies (each element reads its own side).
  std::unique_ptr<GridFunction> surface_response_corrected_flux_sides;
  bool has_surface_response_corrected_field = false;
  bool has_surface_response_corrected_flux = false;
  bool surface_response_field_paraview_registered = false;
  bool surface_response_flux_paraview_registered = false;
  std::optional<std::complex<double>> surface_response_corrected_frequency;
  std::optional<double> surface_response_corrected_mode_overlap;
  struct MaxwellSurfaceResponseMeasurement
  {
    struct Interface
    {
      double raw_energy = 0.0;
      double corrected_energy = 0.0;
      double corrected_energy_fixed_flux = 0.0;
      double self_consistent_energy = 0.0;
      double trace_closure_spread = 0.0;
      double loss_tangent = 0.0;
      // OI-ZERO-IFACE-NAN (decisions 891 / 922 C): a postprocessing interface that is not a
      // response target carries its RAW energy in the corrected / fixed-flux entries (the
      // electrostatic driver's convention) with spread 0, and reads 0 / 0 for the matched /
      // corner fractions; the surface-Q-corrected.csv "response target" column (written
      // only when such an interface exists) marks it.
      bool response_target = true;
    };
    // Per target interface, the uncovered requirements' and the DomainBoundary claims' raw
    // within-R energy of the raw field (fixed trace / fixed flux) and of the corrected
    // field (self-consistent), in total and per feature type: the rows of
    // surface-response-uncovered-energy.csv and -domain-boundary-energy.csv.
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> uncovered;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> self_consistent_uncovered;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> domain_boundary;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> self_consistent_domain_boundary;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> spatial_cluster_priority;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
        self_consistent_spatial_cluster_priority;
    // R0-a on the Maxwell path (idfix DESIGN B.3 (e); S7diag R1, decision 758): the raw
    // claims of the patches the conductor-consistency gate has excluded so far, kept raw
    // like the three lists above (surface-response-conductor-consistency-energy.csv).
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> conductor_consistency;
    std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>
        self_consistent_conductor_consistency;
    // PatchEnergy (decision 922 B-1): per target interface, the raw field's within-R energy
    // per model-set label (+ "unclaimed").
    std::map<int, std::map<std::string, double>> model_raw_partition;
    // Section 8: per target interface, the raw field's within-R energy per exclusion class
    // (surface-response-exclusion-energy.csv); empty without Diagnostics.ExclusionEnergy.
    std::map<int, std::map<std::string, double>> exclusion_energy_partition;
    // Diagnostics.UncoveredFeatureEnergy (decision 1567 (b)): per target interface, the raw
    // field's uncovered within-R energy per (type, feature label); empty when off.
    std::map<int, std::map<std::pair<std::string, std::string>, double>> uncovered_feature;
    // Lane L1 (decision 1634 (1); rulings 1641 / 1649 Q3): THE CLOSED COVERAGE PARTITION of
    // the raw field, per target interface, ALWAYS measured with the postprocessed
    // correction: the kept-raw classes "Missing" (by "m<feature>" label, the key records of
    // SurfaceResponseOperator::GetUncoveredKeys), "DomainBoundary",
    // "SpatialClusterPriority" and "ConductorConsistency" (by portion type), then the
    // served claims (by model-set label), then "UnpatchedPairs" (the served list's
    // `unpatched` label alone) and "Excluded" (by exclusion class;
    // Diagnostics.ExclusionEnergy only), then the unclaimed remainder. The loop remainder
    // (`loop_remainder`) lies inside the served bins: the coverage shares move it out of
    // them.
    std::map<int, SurfacePostOperator::InterfaceCoveragePartition> coverage;
    // Lane L1 (decision 1634 (2) / (4)): the coverage of this measurement per Type (W, U
    // and its shares, the interim coverage status, the covered correction factor and the
    // regime status) and the key rows, from `coverage` and the corrected bookkeeping's W
    // and fabricated energy (ComputeCoverageMeasurement). ROOT ONLY (the writers' rank;
    // the library's qualification records are read there): empty elsewhere. The writers of
    // surface-Q-corrected.csv, surface-response-confidence.csv and the coverage tables read
    // the per-cell status from here (a coverage policy's hook, lane L4).
    CoverageMeasurement coverage_measurement;
    // Lane rescensus-product (decision 1817): the census cells of this measurement (root).
    std::vector<ResolutionCell> resolution_cells;
    // COL-3 (mx-sc-magnetic Phase 3 (b), R-MAG-15; Diagnostics.WithinRVolumeEnergy): per
    // placed patch the raw field's within-box volume energies (empty when off).
    std::vector<SurfaceResponseOperator::PatchVolumeEnergy> within_r_volume_energies;
    // OI-MX-LOOPGATE (decision 955): per target interface, the circulating remainders' raw
    // within-R energy per model row ("m<idx>": rho_p x the raw energy on the patch's
    // claims; the Total is their sum), of the raw field (both closures) and of the
    // corrected field (self-consistent). The rows' patch counts and rho-weighted claim
    // lengths are recomputed by the printer (PrintSurfaceResponseLoopFit).
    std::map<int, std::map<std::string, double>> loop_remainder;
    std::map<int, std::map<std::string, double>> self_consistent_loop_remainder;
    // Lane L1 (ruling 1694 (i); OI-CORRECTED-TOTAL-CC-LOOPREMAINDER-OVERLAP): the part of
    // `loop_remainder` on the conductor-consistency-EXCLUDED patches, by its own integral
    // (the same weighted portion integral restricted to those patches' remainders): the
    // corrected total keeps it raw twice. A coverage output only; no corrected value reads
    // it.
    std::map<int, std::map<std::string, double>> loop_remainder_excluded;

    double raw_normalization_energy = 0.0;
    double corrected_normalization_energy = 0.0;
    double corrected_normalization_energy_fixed_flux = 0.0;
    double self_consistent_normalization_energy = 0.0;
    std::complex<double> self_consistent_frequency = 0.0;
    double self_consistent_mode_overlap = 0.0;
    bool has_postprocessed = false;
    bool has_self_consistent = false;
    // The magnetic (inductive) term (mx-sc-magnetic DESIGN-v2; USER 1221 D6; decision
    // 1229): the first-order fixed-trace frequency shift of the raw mode / frequency,
    // d omega / omega = (U_m,defect - U_e,defect) / (2 U_raw) with U_e,defect the electric
    // domain correction, U_m,defect the magnetic one and U_raw the raw normalization
    // energy; gated by the magnetic term's own loop gate only (decision 1250: NOT by
    // closure_independent_confident nor the trace-closure spread, which are printed beside
    // it; the shift follows the electric correction's coverage rule): NaN when withheld
    // (FixedTraceFrequencyShift). Written only when a magnetic leaf was read (the columns
    // exist then; NaN when no model carries the term). fixed_trace_shift_coverage = the
    // matched energy share the shift rests on (FixedTraceShiftCoverage; decision 1250).
    double fixed_trace_frequency_shift = std::numeric_limits<double>::quiet_NaN();
    double fixed_trace_shift_coverage = std::numeric_limits<double>::quiet_NaN();
    double magnetic_domain_correction = 0.0;
    std::map<int, Interface> interfaces;
    SurfaceResponseOperator::MaxwellResponse confidence;
    SurfaceResponseOperator::MaxwellResponse self_consistent_confidence;
  };
  mutable std::optional<MaxwellSurfaceResponseMeasurement> surface_response_measurement;
  mutable std::unique_ptr<TableWithCSVFile> surface_Q_corrected;
  // surface-Q-corrected.csv carries the "response target" column (decision 922 C-3 (ii)):
  // decided at the first print from the measured interfaces.
  mutable bool surface_response_untargeted_interface = false;
  mutable std::unique_ptr<TableWithCSVFile> surface_response_confidence;
  // Lane L4 (decision 1829; DESIGN-L4-v2 s2.3 M5): whether the last printed measurement's
  // whole-measurement values were marked (Coverage.Policy "MarkNotApplied"; root only;
  // unset when no measurement was printed since the last
  // TakeSurfaceResponseMeasurementMarked).
  mutable std::optional<bool> surface_response_measurement_marked;
  // surface-response-model-energy.csv on the Maxwell path (decisions 845 O-C / 846): the
  // electrostatic driver's column set EXACTLY and in order as a strict prefix (`source` =
  // the excitation index on Driven, the mode index on Eigenmode / BoundaryMode), then the
  // Maxwell key columns "f (GHz)" / "m" and "exc" so that every row of a multi-frequency
  // sweep is unambiguous; rows accumulate over the prints like surface-Q-corrected.csv.
  mutable std::unique_ptr<TableWithCSVFile> surface_response_model_energy;
  // PatchEnergy (decision 922 B): the electrostatic driver's per-patch columns EXACTLY as a
  // strict prefix (surface-response-patch-energy.csv: source, evaluation, patch, model,
  // feature, weight, cell begin / end (m), origin x / y / z (m), domain correction (J),
  // fabricated surface energy[k] (J) per postprocessing interface), then "f (GHz)" / "m"
  // and "exc", then (decision 920) the patch's "max loop residual" (its contour loops' max;
  // 0 for a patch without a path record); surface-response-traces.csv: the electrostatic
  // spatial-trace columns (source, patch, model, coefficient, conductor state, value (V) =
  // the real part) then "value imag (V)", "f (GHz)" / "m", "exc";
  // surface-response-patch-loops.csv (920): source, patch, model, path, kind (0 contour
  // loop / 1 cap-interior hat), loop residual, loop scale (V), "f (GHz)" / "m", "exc".
  mutable std::unique_ptr<TableWithCSVFile> surface_response_patch_energy;
  mutable std::unique_ptr<TableWithCSVFile> surface_response_traces;
  mutable std::unique_ptr<TableWithCSVFile> surface_response_patch_loops;
  // surface-response-magnetic.csv (DESIGN-v2 section 3; decision 1229 Q4): only when a
  // magnetic leaf was read (every verdict, also when every leaf is refused; decision 1235);
  // rewritten at every print (root only) with the per-model leaf status rows and the dual
  // R0-a gate's per-patch rows so far.
  mutable std::unique_ptr<TableWithCSVFile> surface_response_magnetic;
  // COL-1 (mx-sc-magnetic Phase 3 (a), R-MAG-14): the per-model "energy" rows of
  // surface-response-magnetic.csv accumulated over the prints (one per model carrying the
  // term, per source and frequency / mode): the model's applied magnetic domain correction
  // and its within-R magnetic thin energy.
  struct MagneticEnergyRow
  {
    int model = 0;
    int source = 0;
    double frequency_GHz = 0.0;
    double domain_correction = 0.0;
    double thin_energy = 0.0;
  };
  mutable std::vector<MagneticEnergyRow> surface_response_magnetic_energy_rows;
  void PrintSurfaceResponseMagnetic() const;
  void PrintSurfaceResponsePatchEnergies(double output_index, int excitation) const;
  mutable bool surface_response_warning_printed = false;
  mutable bool surface_response_portion_tables_started = false;
  // OI-MAXWELL-CC-CSV-HEADERLESS-AFTER-FIRST-PRINT (lane L1, decision 1634 (6)): the paths
  // of the plain appended tables (the per-type portion tables, the loop-fit tables) whose
  // header is written; a table writes its header with its FIRST rows, whichever print
  // that is.
  mutable std::set<std::string> surface_response_tables_with_header;
  // The per-type portion tables of the Maxwell correction (decision 751): the uncovered
  // requirements' and the DomainBoundary claims' raw within-R energies and shares per
  // evaluation, written like the electrostatic driver's tables of the same names.
  void PrintSurfaceResponsePortionEnergies(double output_index, int excitation) const;
  // The per-model rows of the Maxwell correction (M-5): per model (ev 0 fixed trace, ev 1
  // fixed flux on the raw field; ev 2 self-consistent on the corrected field when it
  // exists) as the electrostatic driver writes them, from the response's ModelContribution
  // entries.
  void PrintSurfaceResponseModelEnergies(double output_index, int excitation) const;
  // PatchEnergy (decision 922 B-1): surface-response-model-raw-energy.csv (per model-set
  // label the raw within-R energy per postprocessing interface) and
  // surface-response-model-raw-groups.csv (per label and model: patches, weight x length),
  // plain text tables appended per print, root only.
  void PrintSurfaceResponseModelRawEnergies(double output_index, int excitation) const;
  mutable bool surface_response_model_raw_started = false;
  // Section 8: surface-response-exclusion-energy.csv (per exclusion class x interface type
  // the raw within-R energy, the type's interface raw energy and the share, with the
  // Maxwell keys), a plain text table appended per print, root only.
  void PrintSurfaceResponseExclusionEnergies(double output_index, int excitation) const;
  mutable bool surface_response_exclusion_energy_started = false;
  // Diagnostics.UncoveredFeatureEnergy (decision 1567 (b)):
  // surface-response-uncovered-feature-energy.csv (per type and feature label of the
  // uncovered requirements: portions, length, bounding box and the raw field's within-R
  // energy per target interface), a plain text table appended per print, root only.
  void PrintSurfaceResponseUncoveredFeatureEnergies(double output_index,
                                                    int excitation) const;
  mutable bool surface_response_uncovered_feature_started = false;
  // OI-MX-LOOPGATE (decision 955): surface-response-loop-fit.csv (per evaluation the
  // maximum and response-weighted circulating share, the remainder / fallback patch counts
  // and the rho-weighted claim length) and surface-response-loop-remainder-energy.csv (per
  // model row with a remainder + Total: patches, weighted claim length, the remainder's raw
  // energy per target interface), plain text tables appended per print, root only; both
  // removed at the first print before the writes (decision 473 MAJOR-3).
  void PrintSurfaceResponseLoopFit(double output_index, int excitation) const;
  mutable bool surface_response_loop_fit_started = false;
  // Lane L1 (decision 1634 (2); DESIGN-v2 s4.2): surface-response-coverage.csv (one row
  // per cell x Type: W, the uncovered / missing / method-gap / pending / no-status shares,
  // the covered correction factor, the validated range, the coverage and regime statuses)
  // and surface-response-uncovered-energy-by-key.csv (per Missing-class key: kind,
  // reason, features, portions, length and the raw within-R energy and share of W per
  // target interface), plain text tables appended per print (each measurement as it
  // completes), root only; both removed at the first print. Records the cells for the
  // end-of-run banner.
  void PrintSurfaceResponseCoverage(double output_index, int excitation) const;
  mutable bool surface_response_coverage_started = false;

  // Reduced-coordinate outputs for adaptive online postprocessing: the electric energy
  // form, the factor R of the curls of the basis vectors (the time-averaged magnetic energy
  // is 1/4 |R y|² / ω²), the domain energy forms (electric in the basis, magnetic in the
  // coordinates R y), and the lumped-port voltages of the basis vectors (one row per port).
  bool reduced_postprocessing_ready = false;
  std::set<int> reduced_postprocessing_checked;  // Excitations checked
  double reduced_postprocessing_tol = 0.0;
  Eigen::MatrixXd reduced_energy_E, reduced_curl_R, reduced_port_V;
  std::map<int, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> reduced_domain_energy;
  // Projections of the unit mass forms of the frequency-dependent materials, weighted at
  // each frequency by their stored-energy susceptibility.
  using ReducedFrequencyDependentForms =
      std::vector<std::pair<std::size_t, Eigen::MatrixXd>>;
  ReducedFrequencyDependentForms reduced_energy_E_fd;
  std::map<int, ReducedFrequencyDependentForms> reduced_domain_energy_E_fd;

  // Per-entry impedance postprocessing configuration (keyed by config index).
  struct ImpedancePostproConfig
  {
    mfem::Array<int> voltage_marker;
    mfem::Array<int> current_marker;
    std::vector<mfem::Vector> voltage_path;
    std::vector<mfem::Vector> current_path;
    bool has_voltage_coordinates = false;
    bool has_current_path = false;
    bool has_current = false;
    std::shared_ptr<fem::LineIntegrator> voltage_integrator, current_integrator;
  };
  std::map<int, ImpedancePostproConfig> impedance_postpro;

  // Per-entry voltage-only postprocessing configuration (keyed by config index).
  struct VoltagePostproConfig
  {
    mfem::Array<int> voltage_marker;
    std::vector<mfem::Vector> voltage_path;
    bool has_coordinates = false;
    std::shared_ptr<fem::LineIntegrator> voltage_integrator;
  };
  std::map<int, VoltagePostproConfig> voltage_postpro;

  // Individual measurements to fill the cache/workspace. Measurements functions are not
  // constrained by solver type in the signature since they are private member functions.
  // They dispatch on solver type within the function itself using `if constexpr`, and do
  // nothing if the measurement is not solver appropriate.
  void MeasureDomainFieldEnergy() const;
  void MeasureLumpedPorts() const;
  void MeasureLumpedPorts(const std::map<int, std::complex<double>> &port_voltages) const;
  void MeasureWavePorts() const;
  void MeasureFloquetPorts() const;
  void MeasureLumpedPortsEig() const;  // Depends: DomainFieldEnergy, LumpedPorts
  void MeasureSParameter() const;      // Depends: LumpedPorts, WavePorts, FloquetPorts
  void MeasureSurfaceFlux() const;
  void MeasureFarField() const;
  void MeasureInterfaceEFieldEnergy() const;      // Depends: LumpedPorts
  void MeasureSurfaceResponseCorrection() const;  // Depends: Domain, ports, interfaces
  // Lane L1: the excitation class of this solver's coverage cells (the validated-regime
  // ranges' key).
  static constexpr CoverageExcitationClass SurfaceResponseCoverageExcitation()
  {
    return solver_t == ProblemType::EIGENMODE ? CoverageExcitationClass::EIGENMODE
           : solver_t == ProblemType::DRIVEN  ? CoverageExcitationClass::DRIVEN
                                              : CoverageExcitationClass::NONE;
  }
  void MeasureProbes() const;
  void PrintSurfaceResponseCorrection(double output_index, int excitation) const;
  void RegisterSurfaceResponseParaviewFields();
  void ClearSurfaceResponseCorrectedField();

  // Helper function called by all solvers. Has to ensure correct call order to deal with
  // dependent measurements.
  void MeasureAllImpl() const
  {
    MeasureDomainFieldEnergy();
    MeasureLumpedPorts();
    MeasureWavePorts();
    MeasureFloquetPorts();
    MeasureLumpedPortsEig();
    MeasureSParameter();
    MeasureSurfaceFlux();
    MeasureInterfaceEFieldEnergy();
    MeasureSurfaceResponseCorrection();
    MeasureProbes();
    MeasureFarField();
  }

  // Setting grid functions.
  //
  // Populate the grid function solutions for the E- and B-field using the solution vectors
  // on the true dofs. For the real-valued overload, the electric scalar potential can be
  // specified too for electrostatic simulations. The output mesh and fields are
  // non-dimensionalized consistently (B ~ E (L₀ ω₀ E₀⁻¹)).
  //
  // These functions are private helper functions. We want to enforce that a caller passes
  // the appropriate ones as part of the MeasureAndPrintAll interface, rather than do a
  // runtime check to see that they have been set.

  void SetEGridFunction(const ComplexVector &e, bool exchange_face_nbr_data = true)
    requires(HasEGridFunction<solver_t>() && HasComplexGridFunction<solver_t>())
  {
    E->Real().SetFromTrueDofs(e.Real());  // Parallel distribute
    E->Imag().SetFromTrueDofs(e.Imag());
    if (exchange_face_nbr_data)
    {
      E->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
      E->Imag().ExchangeFaceNbrData();
    }
  }

  void SetEGridFunction(const Vector &e, bool exchange_face_nbr_data = true)
    requires(HasEGridFunction<solver_t>() && !HasComplexGridFunction<solver_t>())
  {
    E->Real().SetFromTrueDofs(e);
    if (exchange_face_nbr_data)
    {
      E->Real().ExchangeFaceNbrData();
    }
  }

  void SetBGridFunction(const ComplexVector &b, bool exchange_face_nbr_data = true)
    requires(HasBGridFunction<solver_t>() && HasComplexGridFunction<solver_t>())
  {
    B->Real().SetFromTrueDofs(b.Real());  // Parallel distribute
    B->Imag().SetFromTrueDofs(b.Imag());
    if (exchange_face_nbr_data)
    {
      B->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
      B->Imag().ExchangeFaceNbrData();
    }
  }

  void SetBGridFunction(const Vector &b, bool exchange_face_nbr_data = true)
    requires(HasBGridFunction<solver_t>() && !HasComplexGridFunction<solver_t>())
  {
    B->Real().SetFromTrueDofs(b);
    if (exchange_face_nbr_data)
    {
      B->Real().ExchangeFaceNbrData();
    }
  }

  void SetVGridFunction(const Vector &v, bool exchange_face_nbr_data = true)
    requires(HasVGridFunction<solver_t>() && !HasComplexGridFunction<solver_t>())
  {
    V->Real().SetFromTrueDofs(v);
    if (exchange_face_nbr_data)
    {
      V->Real().ExchangeFaceNbrData();
    }
  }

  void SetAGridFunction(const Vector &a, bool exchange_face_nbr_data = true)
    requires(HasAGridFunction<solver_t>() && !HasComplexGridFunction<solver_t>())
  {
    A->Real().SetFromTrueDofs(a);
    if (exchange_face_nbr_data)
    {
      A->Real().ExchangeFaceNbrData();
    }
  }

public:
  PostOperator(const config::ProblemData &problem, const config::SolverData &solver,
               const config::DomainData &domains, const config::BoundaryData &boundaries,
               const Units &units, fem_op_t<solver_t> &fem_op,
               std::shared_ptr<const SurfacePostGeometry> *surface_post_geometry = nullptr);
  // `impedance_candidate_geometry`: the impedance-metal rule candidate's OWN geometry cache
  // (RECORD; held by the driver across AMR iterations like `response_geometry`, so that
  // the candidate keeps the initial mesh's patch layout as the operator of record does).
  explicit PostOperator(
      const IoData &iodata, fem_op_t<solver_t> &fem_op,
      std::shared_ptr<const SurfaceResponseGeometry> *response_geometry = nullptr,
      std::shared_ptr<const SurfacePostGeometry> *surface_post_geometry = nullptr,
      std::shared_ptr<const SurfaceResponseGeometry> *impedance_candidate_geometry =
          nullptr);

  bool NeedsRecoveredElectricFlux() const { return D_recovered != nullptr; }

  // The true-dof size of a recovered electric flux D the postprocessing accepts (the
  // drivers recover it on the error estimator's flux space, SD-403; the archive reduction
  // reads it).
  int GetRecoveredElectricFluxTrueVSize() const
  {
    return flux_space ? flux_space->GetTrueVSize() : fem_op->GetRTSpace().GetTrueVSize();
  }

  // Lane L1 (decision 1634 (2)): the configured Type of every postprocessing interface
  // (only those in `targets` when given): the coverage cells' Types.
  std::map<int, std::string>
  GetInterfaceTypes(const std::set<int> *targets = nullptr) const;

  // PatchEnergy: the applied patches' raw claims as the partition reads them (taken from
  // the operator at the construction, round 2b B1).
  const SurfacePostOperator::CompactPerimeterPortions &
  GetSurfaceResponseModelRawPortions() const
  {
    return surface_response_model_raw;
  }
  SurfaceResponseOperator *GetSurfaceResponseOperator()
  {
    return surface_response_op.get();
  }
  // The impedance-metal rule's palace.json block (SurfaceResponse.ImpedanceRuleCandidate;
  // RECORD only, root's; null under "Off").
  const nlohmann::json *GetImpedanceRuleCandidateRecord() const
  {
    return impedance_film_rule ? &impedance_candidate_record : nullptr;
  }
  // RECORD, at the end of the solve (root; after the candidate block): the candidate's
  // spilled conductor-consistency Records written to <candidate diagnostics directory>/
  // surface-response-impedance-conductor-consistency-records.json, the scratch removed.
  void WriteImpedanceRuleCandidateRecords() const;
  // Lane L4 (M5): the mark of the last printed measurement (root; nullopt when none was
  // printed since the last call), then cleared.
  std::optional<bool> TakeSurfaceResponseMeasurementMarked() const
  {
    auto marked = surface_response_measurement_marked;
    surface_response_measurement_marked.reset();
    return marked;
  }

  const SurfacePostOperator &GetSurfacePostOp() const { return surf_post_op; }

  // Supply a response-corrected Maxwell field for participation postprocessing while the
  // ordinary MeasureAndPrintAll field and CSV output remain the raw thin-metal result.
  template <ProblemType U = solver_t>
  auto SetSurfaceResponseCorrectedField(
      const ComplexVector &e, const ComplexVector *d = nullptr,
      std::optional<std::complex<double>> frequency = std::nullopt,
      std::optional<double> mode_overlap = std::nullopt)
      -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::EIGENMODE, void>
  {
    if (!surface_response_corrected_field)
    {
      surface_response_corrected_field =
          std::make_unique<GridFunction>(fem_op->GetNDSpace(), true);
    }
    surface_response_corrected_field->Real().SetFromTrueDofs(e.Real());
    surface_response_corrected_field->Imag().SetFromTrueDofs(e.Imag());
    surface_response_corrected_field->Real().ExchangeFaceNbrData();
    surface_response_corrected_field->Imag().ExchangeFaceNbrData();
    if (d)
    {
      VerifyRecoveredFluxTrueVSize(d->Size());
      if (!surface_response_corrected_flux)
      {
        surface_response_corrected_flux =
            std::make_unique<GridFunction>(fem_op->GetRTSpace(), true);
      }
      if (flux_space)
      {
        // The output field: the FIRST version of the broken D on the RT space.
        if (!flux_first)
        {
          flux_first = flux_space->BuildVersionMap(FiniteElementSpace::VersionMap::FIRST);
        }
        Vector first(fem_op->GetRTSpace().GetTrueVSize());
        flux_first->Mult(d->Real(), first);
        surface_response_corrected_flux->Real().SetFromTrueDofs(first);
        flux_first->Mult(d->Imag(), first);
        surface_response_corrected_flux->Imag().SetFromTrueDofs(first);
        if (!surface_response_corrected_flux_sides)
        {
          surface_response_corrected_flux_sides = MakeRecoveredFluxGridFunction(true);
        }
        SetRecoveredFluxGridFunction(d->Real(),
                                     surface_response_corrected_flux_sides->Real());
        SetRecoveredFluxGridFunction(d->Imag(),
                                     surface_response_corrected_flux_sides->Imag());
      }
      else
      {
        surface_response_corrected_flux->Real().SetFromTrueDofs(d->Real());
        surface_response_corrected_flux->Imag().SetFromTrueDofs(d->Imag());
      }
      surface_response_corrected_flux->Real().ExchangeFaceNbrData();
      surface_response_corrected_flux->Imag().ExchangeFaceNbrData();
      has_surface_response_corrected_flux = true;
    }
    else
    {
      MFEM_VERIFY(!D_recovered,
                  "Self-consistent Maxwell surface-response postprocessing requires "
                  "recovered electric flux for the corrected field!");
      has_surface_response_corrected_flux = false;
    }
    has_surface_response_corrected_field = true;
    surface_response_corrected_frequency = frequency;
    surface_response_corrected_mode_overlap = mode_overlap;
    RegisterSurfaceResponseParaviewFields();
  }

  template <ProblemType U = solver_t>
  auto SetRecoveredElectricFlux(const ComplexVector &d)
      -> std::enable_if_t<HasEGridFunction<U>() && HasComplexGridFunction<U>(), void>
  {
    MFEM_VERIFY(D_recovered, "Recovered electric flux was not requested!");
    SetRecoveredFluxGridFunction(d.Real(), D_recovered->Real());
    SetRecoveredFluxGridFunction(d.Imag(), D_recovered->Imag());
    recovered_flux_set = true;
  }

  template <ProblemType U = solver_t>
  auto SetRecoveredElectricFlux(const Vector &d)
      -> std::enable_if_t<HasEGridFunction<U>() && !HasComplexGridFunction<U>(), void>
  {
    MFEM_VERIFY(D_recovered, "Recovered electric flux was not requested!");
    SetRecoveredFluxGridFunction(d, D_recovered->Real());
    recovered_flux_set = true;
  }

  // MeasureAndPrintAll is the primary public interface of this class. It is specialized by
  // solver type, since each solver has different fields and extra data required. These
  // functions all:
  // 1) Set the GridFunctions which have to be passed as part of the call.
  // 2) Perform all measurements and populate measurement_cache with temporary results. This
  //    cache structure exists since measurements have dependencies; we may use some
  //    measurement results in later measurements.
  // 3) Pass the measurement cache to the csv printer which will add the appropriate
  //    rows/cols to the csv tables and print to file.
  // 4) Trigger ParaView field computation and save.
  //
  // The functions return the total domain energy which is the only thing needed in the
  // solver to normalize the error indicator, empty for the frequency-domain solvers without
  // a stored electric energy (Measurement::TotalFieldEnergy). If more measurements were
  // needed by the solver loop, we could imagine passing a small struct (like Measurement
  // above or some sub-set therefore).
  //
  // The measure functions will also do logging of (some) measurements to stdout.
  std::optional<double> MeasureAndPrintAll(int ex_idx, int step, const ComplexVector &e,
                                           const ComplexVector &b,
                                           std::complex<double> omega)
    requires(solver_t == ProblemType::DRIVEN);

  // Configure and evaluate the exact reduced-coordinate default output path (domain
  // energies plus port S-parameters). Unsupported configured measurements leave the path
  // disabled and use full postprocessing. The wave port S-parameters use the mode pairings
  // of the last RomOperator::SolvePROM, which must be at omega. The expanded solution e is
  // needed only when ReducedPostprocessingNeedsField is true, otherwise e may be nullptr.
  // tol is the PROM tolerance.
  void ConfigureReducedPostprocessing(const RomOperator &rom_op, double tol)
    requires(solver_t == ProblemType::DRIVEN);
  bool ReducedPostprocessingNeedsField(int ex_idx) const
    requires(solver_t == ProblemType::DRIVEN);
  void MeasureAndPrintReduced(int ex_idx, int step, const ComplexVector *e,
                              std::complex<double> omega, const Eigen::VectorXcd &y,
                              const RomOperator &rom_op)
    requires(solver_t == ProblemType::DRIVEN);
  bool HasReducedPostprocessing() const { return reduced_postprocessing_ready; }

  std::optional<double> MeasureAndPrintAll(int step, const ComplexVector &e,
                                           const ComplexVector &b,
                                           std::complex<double> omega, double error_abs,
                                           double error_bkwd, int num_conv)
    requires(solver_t == ProblemType::EIGENMODE);

  double MeasureAndPrintAll(int step, const Vector &v, const Vector &e, int idx)
    requires(solver_t == ProblemType::ELECTROSTATIC);

  double MeasureAndPrintAll(int step, const Vector &a, const Vector &b, int idx)
    requires(solver_t == ProblemType::MAGNETOSTATIC);

  double MeasureAndPrintAll(int step, const Vector &e, const Vector &b, double time,
                            double J_coef)
    requires(solver_t == ProblemType::TRANSIENT);

  // Mode analysis: complex E-field tangential (ND) and normal (H1) eigenvectors,
  // propagation constant kn, and eigensolver error estimates. The Poynting power is
  // computed internally via fem_op->ComputePoyntingPower.
  std::optional<double> MeasureAndPrintAll(int step, const ComplexVector &et,
                                           const ComplexVector &en, std::complex<double> kn,
                                           double omega, double error_abs,
                                           double error_bkwd, int num_conv)
    requires(solver_t == ProblemType::BOUNDARYMODE);

  // Write error indicator into ParaView file and print summary statistics to csv. Should be
  // called once at the end of the solver loop.
  void MeasureFinalize(const ErrorIndicator &indicator);

  // Measurement of the domain energy without printing. This is needed during the driven
  // simulation with PROM. There samples are taken and we need the total domain energy for
  // the error indicator, but no other measurement / printing should be done.

  std::optional<double> MeasureDomainFieldEnergyOnly(const ComplexVector &e,
                                                     const ComplexVector &b, double omega)
    requires(solver_t == ProblemType::DRIVEN);

  // Access grid functions for field solutions. Note that these are NOT const functions. The
  // electrostatics / magnetostatics solver do measurements of the capacitance/ inductance
  // matrix globally at the end of all solves. This is done in the solver class, but uses
  // the GridFunctions in this (PostOp) class as already allocated scratch workspace.
  //
  // Future: Consider moving those cap/ind measurements into this class and MeasureFinalize?
  // Would need to store vector of V,A.

  GridFunction &GetEGridFunction()
    requires(HasEGridFunction<solver_t>())
  {
    return *E;
  }

  GridFunction &GetBGridFunction()
    requires(HasBGridFunction<solver_t>())
  {
    return *B;
  }

  GridFunction &GetVGridFunction()
    requires(HasVGridFunction<solver_t>())
  {
    return *V;
  }

  GridFunction &GetAGridFunction()
    requires(HasAGridFunction<solver_t>())
  {
    return *A;
  }

  struct ElectrostaticEnergyData
  {
    struct Interface
    {
      double energy;
      double loss_tangent;
      std::vector<SurfacePostOperator::InterfaceEdgeEnergy> edge_energies;
    };

    double domain = 0.0;
    std::map<int, Interface> interfaces;
  };

  // Reuse the raw electrostatic energies already measured by MeasureAndPrintAll.
  template <ProblemType U = solver_t>
  auto GetCachedElectrostaticEnergies() const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, ElectrostaticEnergyData>
  {
    ElectrostaticEnergyData energies;
    MFEM_VERIFY(measurement_cache.domain_E_field_energy_all,
                "Missing cached electrostatic domain energy!");
    energies.domain = *measurement_cache.domain_E_field_energy_all;
    for (const auto &interface : measurement_cache.interface_eps_i)
    {
      energies.interfaces.emplace(interface.idx,
                                  typename ElectrostaticEnergyData::Interface{
                                      interface.energy, interface.tandelta, {}});
    }
    for (const auto &edge : measurement_cache.interface_edge_i)
    {
      auto interface = energies.interfaces.find(edge.idx);
      MFEM_VERIFY(interface != energies.interfaces.end(),
                  "Missing cached electrostatic interface energy!");
      interface->second.edge_energies.push_back(SurfacePostOperator::InterfaceEdgeEnergy{
          edge.distance, edge.energy_outside, edge.energy_annulus});
    }
    return energies;
  }

  // Evaluate electrostatic domain and interface energies without producing ordinary
  // measurement output. This is used to assemble fabrication-corrected observables while
  // preserving the historical raw thin-metal CSV output.
  template <ProblemType U = solver_t>
  auto GetElectrostaticEnergies(const Vector &v, const Vector &e, const Vector *d = nullptr,
                                const std::set<int> *response_targets = nullptr)
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, ElectrostaticEnergyData>
  {
    SetVGridFunction(v);
    SetEGridFunction(e);
    if (NeedsRecoveredElectricFlux())
    {
      MFEM_VERIFY(d, "Recovered electric flux is required for interface postprocessing!");
      SetRecoveredElectricFlux(*d);
    }

    ElectrostaticEnergyData energies;
    energies.domain = dom_post_op.GetElectricFieldEnergy(*V, 0.0);
    surf_post_op.ResetInterfaceLocalEdgeEnergyCache();
    const auto outside_energies =
        response_targets ? surf_post_op.GetInterfaceOuterElectricFieldEnergies(
                               *response_targets, *E, D_recovered.get())
                         : std::map<int, SurfacePostOperator::InterfaceEdgeEnergy>{};
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      if (response_targets && response_targets->find(idx) != response_targets->end())
      {
        auto outside = outside_energies.at(idx);
        energies.interfaces.emplace(
            idx,
            typename ElectrostaticEnergyData::Interface{
                outside.energy_outside, surf_post_op.GetInterfaceLossTangent(idx),
                std::vector<SurfacePostOperator::InterfaceEdgeEnergy>{std::move(outside)}});
      }
      else
      {
        energies.interfaces.emplace(
            idx, typename ElectrostaticEnergyData::Interface{
                     surf_post_op.GetInterfaceElectricFieldEnergy(idx, *E, 0.0,
                                                                  D_recovered.get()),
                     surf_post_op.GetInterfaceLossTangent(idx),
                     surf_post_op.GetInterfaceEdgeElectricFieldEnergies(
                         idx, *E, D_recovered.get())});
      }
    }
    return energies;
  }

  // Uncovered-requirement energies (decision 394 F2) of the field the grid functions hold
  // (the raw field after MeasureAndPrintAll, the corrected field after
  // GetElectrostaticEnergies): per response target interface, the raw within-R interface
  // energy whose nearest perimeter point lies on an uncovered portion, in total and per
  // feature type.
  template <ProblemType U = solver_t>
  auto GetElectrostaticUncoveredEdgeEnergies(
      const std::set<int> &response_targets,
      const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &portions) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::map<int, SurfacePostOperator::UncoveredEdgeEnergy>>
  {
    return surf_post_op.GetInterfaceUncoveredEdgeEnergies(response_targets, *E,
                                                          D_recovered.get(), portions);
  }
  // The exclusion-energy export (idfix DESIGN section 8; decision 1076): the same field's
  // within-R energy partitioned per exclusion class against the CLASSIFIED edge tree
  // (+ "Claimed").
  template <ProblemType U = solver_t>
  auto GetElectrostaticClassifiedEdgeEnergies(
      const std::set<int> &response_targets,
      const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &portions) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::map<int, std::map<std::string, double>>>
  {
    return surf_post_op.GetInterfaceClassifiedEdgeEnergies(response_targets, *E,
                                                           D_recovered.get(), portions);
  }
  // PatchEnergy (decision 922 B-1): the same field's within-R energy partitioned per
  // model-set label over the applied patches' raw claims.
  template <ProblemType U = solver_t>
  auto GetElectrostaticPartitionedEdgeEnergies(
      const std::set<int> &response_targets,
      SurfacePostOperator::PerimeterPortionView portions) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::map<int, std::map<std::string, double>>>
  {
    return surf_post_op.GetInterfacePartitionedEdgeEnergies(response_targets, *E,
                                                            D_recovered.get(), portions);
  }
  // Lane L1 (decision 1634 (1); rulings 1641 / 1649 Q3): the closed coverage partition of
  // the same field (SurfacePostOperator::GetInterfaceCoveragePartition).
  template <ProblemType U = solver_t>
  auto GetElectrostaticCoveragePartition(
      const std::set<int> &response_targets,
      const std::vector<SurfacePostOperator::CoveragePartitionClass> &kept_classes,
      SurfacePostOperator::PerimeterPortionView served,
      const std::vector<SurfacePostOperator::CoveragePartitionClass> &remainder_classes,
      const SurfacePostOperator::CoverageDroppedWindows *dropped_windows = nullptr,
      const std::set<int> &excluded_patches = {}) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::map<int, SurfacePostOperator::InterfaceCoveragePartition>>
  {
    return surf_post_op.GetInterfaceCoveragePartition(
        response_targets, *E, D_recovered.get(), kept_classes, served, remainder_classes,
        dropped_windows, excluded_patches);
  }
  // Ruling 2251: the applied patches' footprints on the targets with a Region, of the same
  // field (SurfacePostOperator::GetInterfaceRegionPatchFootprints).
  template <ProblemType U = solver_t>
  auto GetElectrostaticRegionPatchFootprints(
      const std::set<int> &region_targets,
      const std::vector<SurfacePostOperator::CoveragePartitionClass> &kept_classes,
      const std::vector<SurfacePostOperator::UncoveredPerimeterPortion> &claims,
      const std::vector<std::vector<std::size_t>> &claim_patches,
      const std::vector<std::array<double, 3>> &reference_points) const
      -> std::enable_if_t<
          U == ProblemType::ELECTROSTATIC,
          std::map<int, std::vector<SurfacePostOperator::RegionPatchFootprint>>>
  {
    return surf_post_op.GetInterfaceRegionPatchFootprints(
        region_targets, *E, D_recovered.get(), kept_classes, claims, claim_patches,
        reference_points);
  }
  // Ruling 2251: the entries of `indices` that carry a Region.
  std::set<int> GetRegionInterfaces(const std::set<int> &indices) const
  {
    return surf_post_op.GetRegionInterfaces(indices);
  }

  // Evaluate localized interface energies for an electrostatic field without producing
  // ordinary measurement output. Volume-annulus diagnostics are omitted because response
  // matrix calibration uses only the surface energies.
  template <ProblemType U = solver_t>
  auto GetInterfaceLocalEdgeElectricFieldEnergies(const Vector &e,
                                                  const Vector *d = nullptr)
      -> std::enable_if_t<
          U == ProblemType::ELECTROSTATIC,
          std::map<int, std::vector<SurfacePostOperator::InterfaceLocalEdgeEnergy>>>
  {
    SetEGridFunction(e);
    if (NeedsRecoveredElectricFlux())
    {
      MFEM_VERIFY(d, "Recovered electric flux is required for interface postprocessing!");
      SetRecoveredElectricFlux(*d);
    }

    std::map<int, std::vector<SurfacePostOperator::InterfaceLocalEdgeEnergy>> energies;
    surf_post_op.ResetInterfaceLocalEdgeEnergyCache();
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      auto interface_energy = surf_post_op.GetInterfaceLocalEdgeElectricFieldEnergies(
          idx, *E, D_recovered.get(), false);
      if (!interface_energy.empty())
      {
        energies.emplace(idx, std::move(interface_energy));
      }
    }
    return energies;
  }

  // Batch aggregate interface response matrices for electrostatic basis fields. Grid
  // functions are materialized once and reused for one quadrature traversal per interface.
  template <ProblemType U = solver_t>
  auto GetInterfaceElectricFieldEnergyMatrices(
      const std::vector<Vector> &e, const std::vector<Vector> *d = nullptr,
      int quadrature_extra = 0,
      std::vector<SurfacePostOperator::InterfaceQuadratureRule> *quadrature_rules = nullptr)
      -> std::enable_if_t<
          U == ProblemType::ELECTROSTATIC,
          std::map<int, std::vector<SurfacePostOperator::InterfaceResponseMatrix>>>
  {
    MFEM_VERIFY(!e.empty() && (!d || d->size() == e.size()),
                "Invalid batched electrostatic response fields!");
    MFEM_VERIFY(!NeedsRecoveredElectricFlux() || d,
                "Recovered electric flux is required for interface postprocessing!");
    MFEM_VERIFY(!d || D_recovered,
                "Recovered electric flux fields were supplied but not configured!");

    std::vector<std::unique_ptr<GridFunction>> e_grid, d_grid;
    std::vector<const GridFunction *> e_ptr, d_ptr;
    e_grid.reserve(e.size());
    e_ptr.reserve(e.size());
    if (d)
    {
      d_grid.reserve(d->size());
      d_ptr.reserve(d->size());
    }
    for (std::size_t i = 0; i < e.size(); i++)
    {
      auto field = std::make_unique<GridFunction>(*E->ParFESpace());
      field->Real().SetFromTrueDofs(e[i]);
      field->Real().ExchangeFaceNbrData();
      e_ptr.push_back(field.get());
      e_grid.push_back(std::move(field));
      if (d)
      {
        auto flux = MakeRecoveredFluxGridFunction(false);
        SetRecoveredFluxGridFunction((*d)[i], flux->Real());
        d_ptr.push_back(flux.get());
        d_grid.push_back(std::move(flux));
      }
    }
    return surf_post_op.GetInterfaceElectricFieldEnergyMatrices(
        e_ptr, d_ptr, quadrature_extra, quadrature_rules);
  }

  // P3-16: a batched interface response (the two forms above and below) cannot represent
  // the summed energies of both sides of an interior sheet, so the runs that use it refuse
  // such a sample before the solve (SurfacePostOperator::VerifyBatchedInterfaceResponse):
  // the constructor for the aggregated response matrix (also the response archive's
  // writer and reducer), the archive estimate itself.
  template <ProblemType U = solver_t>
  auto VerifyBatchedInterfaceResponse() const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, void>
  {
    surf_post_op.VerifyBatchedInterfaceResponse(*E);
  }

  // Streaming one-pass form of the above (the archive reduction, decision 62(4)): the
  // interface samples are cached once, each basis field is set as the E grid function
  // (and the recovered flux, when archived) and evaluated once into its amplitude rows,
  // then the Gram matrices are assembled per interface from the rows.
  template <ProblemType U = solver_t>
  auto CacheInterfaceResponseSamples() const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::vector<SurfacePostOperator::InterfaceResponseSamples>>
  {
    return surf_post_op.CacheInterfaceResponseSamples();
  }

  template <ProblemType U = solver_t>
  auto SetInterfaceResponseField(const Vector &e, const Vector *d)
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, void>
  {
    MFEM_VERIFY(!NeedsRecoveredElectricFlux() || d,
                "Recovered electric flux is required for interface postprocessing!");
    MFEM_VERIFY(!d || D_recovered,
                "Recovered electric flux fields were supplied but not configured!");
    E->Real().SetFromTrueDofs(e);
    E->Real().ExchangeFaceNbrData();
    if (d)
    {
      SetRecoveredElectricFlux(*d);
    }
  }

  template <ProblemType U = solver_t>
  auto
  EvaluateInterfaceResponseRow(const SurfacePostOperator::InterfaceResponseSamples &samples,
                               double *row) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, void>
  {
    surf_post_op.EvaluateInterfaceResponseRow(samples, *E, D_recovered.get(), row);
  }

  template <ProblemType U = solver_t>
  auto AssembleInterfaceResponseMatrices(
      const SurfacePostOperator::InterfaceResponseSamples &samples, const double *rows,
      int basis_size) const
      -> std::enable_if_t<U == ProblemType::ELECTROSTATIC,
                          std::vector<SurfacePostOperator::InterfaceResponseMatrix>>
  {
    return surf_post_op.AssembleInterfaceResponseMatrices(samples, rows, basis_size,
                                                          GetComm());
  }

  // Whether impedance/voltage postprocessing is configured (mode analysis).
  bool HasImpedancePostprocessing() const { return !impedance_postpro.empty(); }
  bool HasCurrent() const
  {
    return std::any_of(impedance_postpro.begin(), impedance_postpro.end(),
                       [](const auto &p) { return p.second.has_current; });
  }
  bool HasVoltagePostprocessing() const { return !voltage_postpro.empty(); }

  // Access to number of padding digits.
  static constexpr auto GetPadDigitsDefault() { return pad_digits_default; }

  // Access to domain postprocessing objects. Use in electrostatic & magnetostatic matrix
  // measurement (see above).
  const auto &GetDomainPostOp() const { return dom_post_op; }

  // The electric and magnetic field energies of the last measurement, nondimensional: the
  // electric field energy with the high-frequency permittivity ε∞ of the
  // frequency-dependent materials
  // (DomainPostOperator::GetHighFrequencyElectricFieldEnergy), which, unlike the stored and
  // the material electric energies, is available for every material and does not vanish at
  // a zero of the permittivity, and the magnetic energy.
  std::pair<double, double> GetDomainFieldEnergies() const
    requires(HasEGridFunction<solver_t>())
  {
    return {dom_post_op.GetHighFrequencyElectricFieldEnergy(*E),
            measurement_cache.domain_H_field_energy_all};
  }

  // The L-vector of a magnetostatic solution, on the space of the domain energy operators
  // (broken across two-sided superconductor films, if any).
  void GetSolutionLVector(const Vector &a, Vector &a_l_out) const
    requires(solver_t == ProblemType::MAGNETOSTATIC)
  {
    const auto &nd_fespace = fem_op->GetSolutionNDSpace();
    a_l_out.SetSize(nd_fespace.GetVSize());
    a_l_out.UseDevice(true);
    nd_fespace.GetProlongationMatrix()->Mult(a, a_l_out);
  }

  // Expose MPI communicator from fem_op for electrostatic & magnetostatic matrix processing
  // (see above).
  auto GetComm() const { return fem_op->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_HPP
