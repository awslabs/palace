// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "postoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <complex>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "drivers/boundarymodesolver.hpp"
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/face_sampling_plan.hpp"
#include "fem/interpolator.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/operator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/floquetportoperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/constants.hpp"
#include "utils/enum_string.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

#include <mfem/config/config.hpp>

namespace palace
{

using namespace std::complex_literals;

namespace
{

constexpr double minimum_eigenmode_response_overlap = 0.8;

// The RT elements of a space (the same finite elements, so the same basis and Piola map)
// with every DOF interior to its element: the space of a field which is discontinuous
// across every face, in the elements' native basis (no DOF transformation). It holds the
// restriction to each element of a field of the flux space broken across the interior
// sheet faces (P3-14, SD-3202), so that the coefficients that read D per element see the
// version of the element's own side.
class DiscontinuousRTFECollection : public mfem::FiniteElementCollection
{
  const mfem::FiniteElementCollection &rt_fec;
  const int dim;
  std::string name;

public:
  DiscontinuousRTFECollection(const mfem::FiniteElementCollection &rt_fec, int dim)
    : mfem::FiniteElementCollection(rt_fec.GetOrder()), rt_fec(rt_fec), dim(dim),
      name(std::string("Discontinuous_") + rt_fec.Name())
  {
  }

  const mfem::FiniteElement *
  FiniteElementForGeometry(mfem::Geometry::Type geom) const override
  {
    return (mfem::Geometry::Dimension[geom] == dim) ? rt_fec.FiniteElementForGeometry(geom)
                                                    : nullptr;
  }

  int DofForGeometry(mfem::Geometry::Type geom) const override
  {
    const auto *fe = FiniteElementForGeometry(geom);
    return fe ? fe->GetDof() : 0;
  }

  const int *DofOrderForOrientation(mfem::Geometry::Type, int) const override
  {
    return nullptr;
  }

  const char *Name() const override { return name.c_str(); }

  int GetContType() const override { return DISCONTINUOUS; }
};

bool IsSupportedDomainOutputDimension(const mfem::ParMesh &mesh)
{
  return (mesh.Dimension() == 2 && mesh.SpaceDimension() == 2) ||
         (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3);
}

bool IsSupportedBoundaryOutputDimension(const mfem::ParMesh &mesh)
{
  return IsSupportedDomainOutputDimension(mesh);
}

bool UseCeedDomainParaviewPointFields(const mfem::ParMesh &mesh)
{
  return IsSupportedDomainOutputDimension(mesh);
}

bool UseCeedBoundaryParaviewPointFields(const mfem::ParMesh &mesh)
{
  return IsSupportedBoundaryOutputDimension(mesh);
}

void RequireCeedPointFieldEvaluator(const PointFieldEvaluator *eval,
                                    const std::string &what)
{
  MFEM_VERIFY(eval && eval->IsValid(), "libCEED postprocessing was expected for " + what +
                                           ", but PointFieldEvaluator could not assemble!");
}

std::string OutputFolderName(const ProblemType solver_t)
{
  switch (solver_t)
  {
    case ProblemType::DRIVEN:
      return "driven";
    case ProblemType::EIGENMODE:
      return "eigenmode";
    case ProblemType::ELECTROSTATIC:
      return "electrostatic";
    case ProblemType::MAGNETOSTATIC:
      return "magnetostatic";
    case ProblemType::TRANSIENT:
      return "transient";
    case ProblemType::BOUNDARYMODE:
      return "boundarymode";
    default:
      return "unknown";
  }
}

void ScaleSurfaceResponseGridFunctions(double L, int dim, GridFunction *E, GridFunction *D)
{
  if (E)
  {
    E->Real() *= L;
    E->Real().FaceNbrData() *= L;
    E->Imag() *= L;
    E->Imag().FaceNbrData() *= L;
  }
  if (D)
  {
    const auto Ld = std::pow(L, dim - 1);
    D->Real() *= Ld;
    D->Real().FaceNbrData() *= Ld;
    D->Imag() *= Ld;
    D->Imag().FaceNbrData() *= Ld;
  }
}

void DimensionalizeSurfaceResponseGridFunctions(Units &units, GridFunction *E,
                                                GridFunction *D)
{
  if (E)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_E>(*E);
  }
  if (D)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_D>(*D);
  }
}

void NondimensionalizeSurfaceResponseGridFunctions(Units &units, GridFunction *E,
                                                   GridFunction *D)
{
  if (E)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(*E);
  }
  if (D)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_D>(*D);
  }
}

}  // namespace

double FixedTraceFrequencyShift(const SurfaceResponseOperator::MaxwellResponse &response,
                                double raw_normalization_energy)
{
  // Decision 1250: not gated by closure_independent_confident (the electric correction's
  // coverage rule applies: printed on partial coverage, the uncovered part kept raw).
  if (!response.has_magnetic || !response.magnetic_confident)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  MFEM_VERIFY(raw_normalization_energy > 0.0,
              "The fixed-trace frequency shift requires a positive raw normalization "
              "energy!");
  return (response.magnetic_domain_correction - response.domain_correction) /
         (2.0 * raw_normalization_energy);
}

double FixedTraceShiftCoverage(double within_radius_raw_energy, double kept_raw_energy)
{
  if (!(within_radius_raw_energy > 0.0))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 1.0 - kept_raw_energy / within_radius_raw_energy;
}

template <ProblemType solver_t>
PostOperator<solver_t>::PostOperator(
    const config::ProblemData &problem, const config::SolverData &solver,
    const config::DomainData &domains, const config::BoundaryData &boundaries,
    const Units &units_, fem_op_t<solver_t> &fem_op_,
    std::shared_ptr<const SurfacePostGeometry> *surface_post_geometry)
  : fem_op(&fem_op_), units(units_), post_dir(problem.output),
    post_op_csv(problem, solver, boundaries, units_, fem_op_),
    // dom_post_op does not have a default ctor so specialize via immediate lambda.
    dom_post_op(std::move(
        [&domains, &fem_op_]()
        {
          if constexpr (solver_t == ProblemType::ELECTROSTATIC)
          {
            return DomainPostOperator(domains.postpro, fem_op_.GetMaterialOp(),
                                      fem_op_.GetH1Space());
          }
          else if constexpr (solver_t == ProblemType::MAGNETOSTATIC)
          {
            // The space of the solution: broken across two-sided superconductor films.
            return DomainPostOperator(domains.postpro, fem_op_.GetMaterialOp(),
                                      fem_op_.GetSolutionNDSpace());
          }
          else if constexpr (solver_t == ProblemType::BOUNDARYMODE)
          {
            // Mode analysis: E on ND space, B (Hz) on L2 curl space.
            return DomainPostOperator(domains.postpro, fem_op_.GetMaterialOp(),
                                      fem_op_.GetNDSpace(), fem_op_.GetCurlSpace());
          }
          else
          {
            return DomainPostOperator(domains.postpro, fem_op_.GetMaterialOp(),
                                      fem_op_.GetNDSpace(), fem_op_.GetCurlSpace());
          }
        }())),
    surf_post_op(boundaries.postpro, solver_t, fem_op->GetMaterialOp(),
                 fem_op->GetH1Space(), fem_op->GetNDSpace(),
                 boundaries.GetSheetAttributes(), &boundaries, surface_post_geometry,
                 ExclusionEnergyRequested(solver)),
    interp_op(domains.postpro.probe, units_, fem_op->GetNDSpace())
{
  // Define primary grid-functions.
  if constexpr (HasVGridFunction<solver_t>())
  {
    V = std::make_unique<GridFunction>(fem_op->GetH1Space());
  }
  if constexpr (HasAGridFunction<solver_t>())
  {
    A = std::make_unique<GridFunction>(fem_op->GetNDSpace());
  }
  if constexpr (HasEGridFunction<solver_t>())
  {
    E = std::make_unique<GridFunction>(fem_op->GetNDSpace(),
                                       HasComplexGridFunction<solver_t>());
    if (surf_post_op.NeedsFluxRecovery())
    {
      // The drivers recover D with the error estimator, on its flux space: broken across
      // the interior faces of the sheet attributes whenever the mesh has any (SD-403). The
      // same rule on the same RT space and attributes gives the same view here (P3-14,
      // SD-3202; OI-UPSYNC-FLUXRECOVERY-REFUSED-ON-BROKEN-RECOVERY-SPACE of decision 1578
      // lifted): D is then held per side for the energies.
      flux_space =
          BuildBrokenFluxSpace(fem_op->GetRTSpace(), boundaries.GetSheetAttributes());
      D_recovered = MakeRecoveredFluxGridFunction(HasComplexGridFunction<solver_t>());
    }
    if constexpr (solver_t == ProblemType::ELECTROSTATIC)
    {
      if (solver.electrostatic.response_matrix &&
          solver.electrostatic.aggregate_response_matrix)
      {
        VerifyBatchedInterfaceResponse();
      }
    }
  }
  if constexpr (HasBGridFunction<solver_t>())
  {
    B = std::make_unique<GridFunction>(fem_op_.GetCurlSpace(),
                                       HasComplexGridFunction<solver_t>());
  }

  // Mode analysis: create grid functions for out-of-plane E (H1) and in-plane B (ND).
  if constexpr (solver_t == ProblemType::BOUNDARYMODE)
  {
    En = std::make_unique<GridFunction>(fem_op->GetH1Space(), true);
    Bt_inplane = std::make_unique<GridFunction>(fem_op->GetNDSpace(), true);
  }

  // Add wave port boundary mode postprocessing, if available.
  if constexpr (std::is_same_v<fem_op_t<solver_t>, SpaceOperator>)
  {
    // Add scaling factor to output dimensional wave port electric fields.
    const double scaling = units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      auto ret = port_E0.emplace(idx, WavePortFieldData());
      ret.first->second.E0r = data.GetModeFieldCoefficientReal(scaling);
      ret.first->second.E0i = data.GetModeFieldCoefficientImag(scaling);
    }
  }

  // Prepare for saving fields
  enable_paraview_output = problem.output_formats.paraview;
  enable_gridfunction_output = problem.output_formats.gridfunction;
  if (solver_t == ProblemType::DRIVEN)
  {
    output_save_indices = solver.driven.save_indices;
  }
  else if (solver_t == ProblemType::EIGENMODE)
  {
    output_n_post = solver.eigenmode.n_post;
  }
  else if (solver_t == ProblemType::ELECTROSTATIC)
  {
    output_n_post = solver.electrostatic.n_post;
  }
  else if (solver_t == ProblemType::MAGNETOSTATIC)
  {
    output_n_post = solver.magnetostatic.n_post;
  }
  else if (solver_t == ProblemType::TRANSIENT)
  {
    output_delta_post = solver.transient.delta_post;
  }
  else if (solver_t == ProblemType::BOUNDARYMODE)
  {
    output_n_post = solver.boundary_mode.n_post;
  }

  // Mode analysis impedance postprocessing setup.
  if constexpr (solver_t == ProblemType::BOUNDARYMODE)
  {
    // Helper to convert coordinate path from config to mfem::Vector list.
    auto ParsePath = [](const std::vector<std::vector<double>> &config_path)
    {
      std::vector<mfem::Vector> path;
      for (const auto &pt : config_path)
      {
        mfem::Vector p(pt.size());
        for (int d = 0; d < static_cast<int>(pt.size()); d++)
        {
          p(d) = pt[d];  // Already nondimensionalized.
        }
        path.push_back(std::move(p));
      }
      return path;
    };

    const auto &pmesh = fem_op->GetNDSpace().GetParMesh();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;

    // Coordinate-path line integrals are split at element boundaries and each piece is
    // integrated with a Gauss rule exact for the field on affine elements, see
    // fem::LineIntegrator.
    const int field_order = fem_op->GetNDSpace().Get().GetMaxElementOrder();

    // A path which misses the mesh would give zero, so reject it as for wave port voltage
    // paths, and warn if only part of it lies on the mesh.
    auto MakeLineIntegrator =
        [&](const std::vector<mfem::Vector> &path, bool closed, const char *name, int idx)
    {
      auto integrator =
          std::make_shared<fem::LineIntegrator>(pmesh, path, closed, field_order);
      const double rel_uncovered =
          integrator->GetUncoveredLength() / integrator->GetPathLength();
      MFEM_VERIFY(rel_uncovered < 1.0 - 1.0e-12, name << " of postprocessing index " << idx
                                                      << " does not intersect the mesh!");
      if (rel_uncovered > 1.0e-9)
      {
        Mpi::Warning(pmesh.GetComm(),
                     "{:.3e}% of the {} of postprocessing index {:d} lies outside the "
                     "mesh and does not contribute!\n",
                     100.0 * rel_uncovered, name, idx);
      }
      return integrator;
    };

    // Impedance postprocessing setup (all configured entries).
    for (const auto &[idx, imp] : boundaries.postpro.impedance)
    {
      auto &cfg = impedance_postpro[idx];
      if (!imp.voltage_attributes.empty())
      {
        cfg.voltage_marker = mesh::AttrToMarker(bdr_attr_max, imp.voltage_attributes);
      }
      if (!imp.current_attributes.empty())
      {
        cfg.current_marker = mesh::AttrToMarker(bdr_attr_max, imp.current_attributes);
        cfg.has_current = true;
      }
      if (imp.voltage_path.size() >= 2)
      {
        cfg.has_voltage_coordinates = true;
        cfg.voltage_path = ParsePath(imp.voltage_path);
        cfg.voltage_integrator =
            MakeLineIntegrator(cfg.voltage_path, false, "Impedance \"VoltagePath\"", idx);
      }
      if (!imp.current_path.empty())
      {
        cfg.has_current_path = true;
        cfg.has_current = true;
        cfg.current_path = ParsePath(imp.current_path);
        cfg.current_integrator =
            MakeLineIntegrator(cfg.current_path, true, "Impedance \"CurrentPath\"", idx);
      }
    }

    // Voltage-only postprocessing setup (all configured entries).
    for (const auto &[idx, vol] : boundaries.postpro.voltage)
    {
      auto &cfg = voltage_postpro[idx];
      if (!vol.voltage_attributes.empty())
      {
        cfg.voltage_marker = mesh::AttrToMarker(bdr_attr_max, vol.voltage_attributes);
      }
      if (vol.voltage_path.size() >= 2)
      {
        cfg.has_coordinates = true;
        cfg.voltage_path = ParsePath(vol.voltage_path);
        cfg.voltage_integrator =
            MakeLineIntegrator(cfg.voltage_path, false, "Voltage \"VoltagePath\"", idx);
      }
    }
  }

  // Remove previous symlink/directory before writing (node-local filesystem aware, so
  // stale artifacts are cleared on every node rather than only where the global root
  // lives). Removing directory only happens when the output folder is dirty.
  auto gridfunction_root = post_dir / "gridfunction";
  RemovePreviousOutput(gridfunction_root, fem_op->GetComm());
  gridfunction_output_dir = (gridfunction_root / OutputFolderName(solver_t)).string();

  // Clear previous ParaView output once per solve. Switching driven excitations only
  // changes subdirectories and must preserve output from earlier excitations. Field
  // evaluators and collections remain lazy so setup is paid only after a solution exists.
  if (ShouldWriteParaviewFields())
  {
    RemovePreviousOutput(post_dir / "paraview", fem_op->GetComm());
  }

  // Initialize CSV files for measurements.
  post_op_csv.InitializeCSVDataCollection(*this);
}

template <ProblemType solver_t>
PostOperator<solver_t>::PostOperator(
    const IoData &iodata, fem_op_t<solver_t> &fem_op_,
    std::shared_ptr<const SurfaceResponseGeometry> *response_geometry,
    std::shared_ptr<const SurfacePostGeometry> *surface_post_geometry)
  : PostOperator(iodata.problem, iodata.solver, iodata.domains, iodata.boundaries,
                 iodata.units, fem_op_, surface_post_geometry)
{
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::BOUNDARYMODE)
  {
    if (iodata.solver.surface_response_correction)
    {
      const auto &response = *iodata.solver.surface_response_correction;
      surface_response_postprocess_enabled = response.IncludesPostprocessing();
      surface_response_self_consistent_enabled =
          solver_t != ProblemType::BOUNDARYMODE && response.IncludesSelfConsistent();
      surface_response_patch_energy_enabled = response.patch_energy;
      surface_response_uncovered_feature_energy_enabled = response.uncovered_feature_energy;
      surface_response_op =
          std::make_unique<SurfaceResponseOperator>(iodata, fem_op_, response_geometry);
      // Lane fix-response-memory (decisions 1446 / 1447): the gathered per-patch records of
      // every measurement live on the root only (their consumers are the root's prints).
      surface_response_op->SetPatchRecordResidency(
          SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY);
      // Round 2a R2 (decision 1475 (1)): the root's Records spilled to the output
      // directory.
      surface_response_op->SpillConductorConsistencyRecords(iodata.problem.output);
      if (surface_response_op->IsRawFallback())
      {
        // A whole-run geometry fallback (USER decision 1000 (b) / 1001): no model / patch
        // applied; the drivers save the operator's record + banner and the Maxwell outputs
        // are the raw run's (nothing measured or printed for the correction).
        surface_response_postprocess_enabled = false;
        surface_response_self_consistent_enabled = false;
        surface_response_patch_energy_enabled = false;
        surface_response_uncovered_feature_energy_enabled = false;
      }
      else
      {
        MFEM_VERIFY(surface_response_op->HasSurfaceResponse(),
                    "Maxwell surface response correction requires fabricated and thin "
                    "surface response matrices in the selected coupon models!");
      }
      for (const auto &portion : surface_response_op->GetUncoveredPortions())
      {
        surface_response_uncovered.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
      for (const auto &portion : surface_response_op->GetDomainBoundaryPortions())
      {
        surface_response_domain_boundary.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
      for (const auto &portion : surface_response_op->GetSpatialClusterPriorityPortions())
      {
        surface_response_spatial_cluster_priority.push_back(
            {portion.p0, portion.p1, portion.topology, portion.feature});
      }
      // Round 2b B1 (decision 1484): moved out of the operator (the same portions in the
      // same order; the operator reads them nowhere after the construction), held once,
      // in the compact exact form (lane L1, ruling 1682).
      surface_response_model_raw = surface_response_op->ReleaseModelRawClaimPortionList();
      // Ruling 1688 (lane L1): the stage-0 Missing LENGTH line (root; log only).
      if (!surface_response_op->IsRawFallback())
      {
        const auto [missing, method_gap] = SplitCoverageUncoveredByKeyClass(
            surface_response_uncovered, surface_response_op->GetCoverageKeyRecords());
        const auto targets = surface_response_op->GetTargetInterfaces();
        Mpi::Print(
            "\n{}",
            CoverageStageZeroLine(
                GetInterfaceTypes(), targets,
                surf_post_op.GetInterfacePerimeterLengths(targets, {missing, method_gap}),
                surface_response_op->GetIdentificationFeatureLengths(), iodata.units));
        const auto not_served = surface_response_op->GetConfiguredTargetsNotServed();
        if (!not_served.empty())
        {
          Mpi::Warning("{}",
                       CoverageTargetsNotServedWarning(GetInterfaceTypes(&not_served)));
        }
        // Lane rescensus-product (decision 1817; ruling 1819 F1): the census Types (the
        // coverage interfaces') and the PREFLIGHT line (collective; before the first
        // solve).
        const auto census_interfaces = surface_response_op->GetCoverageInterfaces();
        surface_response_op->ConfigureResolutionTypes(
            GetInterfaceTypes(&census_interfaces));
      }
      surface_response_exclusion_energy_enabled =
          surface_response_op->ExclusionEnergyRequested();
      surface_response_within_r_volume_energy_enabled =
          surface_response_patch_energy_enabled &&
          surface_response_op->WithinRVolumeEnergyRequested();
      if (surface_response_exclusion_energy_enabled)
      {
        for (const auto &portion : surface_response_op->GetExclusionEnergyPortions())
        {
          surface_response_exclusion_energy.push_back(
              {portion.p0, portion.p1, portion.topology, portion.feature});
        }
      }
      if (surface_response_patch_energy_enabled && Mpi::Root(fem_op->GetComm()) &&
          !surface_response_op->GetPatchAssignments().empty())
      {
        // PatchEnergy (decision 922 B): the applied patch table at construction, the
        // electrostatic driver's column form (no solve-time gate rewrites it here).
        TableWithCSVFile patches(
            surface_response_op->DiagnosticsPath(post_dir, "surface-response-patches.csv"));
        patches.table.insert(Column::Integer("model", "model"));
        for (const char *coordinate : {"x", "y", "z"})
        {
          patches.table.insert(fmt::format("origin_{}", coordinate),
                               fmt::format("origin {} (m)", coordinate));
        }
        for (const char *axis : {"u", "v", "w"})
        {
          for (const char *coordinate : {"x", "y", "z"})
          {
            patches.table.insert(fmt::format("axis_{}_{}", axis, coordinate),
                                 fmt::format("axis {} {}", axis, coordinate));
          }
        }
        patches.table.insert("weight", "weight");
        for (const auto &patch : surface_response_op->GetPatchAssignments())
        {
          patches.table["model"] << patch.model;
          for (int d = 0; d < 3; d++)
          {
            patches.table[fmt::format("origin_{}", "xyz"[d])]
                << units.Dimensionalize<Units::ValueType::LENGTH>(patch.origin[d]);
            patches.table[fmt::format("axis_u_{}", "xyz"[d])] << patch.axis_u[d];
            patches.table[fmt::format("axis_v_{}", "xyz"[d])] << patch.axis_v[d];
            patches.table[fmt::format("axis_w_{}", "xyz"[d])] << patch.axis_w[d];
          }
          patches.table["weight"] << patch.weight;
        }
        patches.WriteFullTableTrunc();
      }
      if (surface_response_patch_energy_enabled)
      {
        // Lane rescensus-product (decision 1817): the census's per-patch readings
        // (PatchEnergy only; collective).
        surface_response_op->WriteResolutionPatchTable(
            surface_response_op
                ->DiagnosticsPath(post_dir, "surface-response-resolution-patches.csv")
                .string());
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::InitializeParaviewDataCollection(int ex_idx)
  requires(solver_t == ProblemType::DRIVEN)
{
  fs::path sub_folder_name = "";
  auto nr_excitations = fem_op->GetPortExcitations().Size();
  if ((nr_excitations > 1) && (ex_idx > 0))
  {
    int spacing = 1 + int(std::log10(nr_excitations));
    sub_folder_name = fmt::format(FMT_STRING("excitation_{:0>{}}"), ex_idx, spacing);
  }
  InitializeParaviewDataCollection(sub_folder_name);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::SetupFieldCoefficients()
{
  // We currently don't use the dependent grid functions apart from saving fields, so only
  // initialize if needed.
  if (!ShouldWriteFields())
  {
    return;
  }
  MFEM_VERIFY(!field_coefficients_initialized,
              "Field coefficients should only be initialized once!");
  field_coefficients_initialized = true;

  // Initialize the (interpolatory L2) output spaces for the libCEED-evaluated
  // visualization fields. The order matches the ParaView output sampling lattice
  // (SetLevelsOfDetail below), preserving the existing VTU point layout.
  // TODO: Revisit the output order choice: order p smooths the (degree 2p) energy
  // densities slightly relative to the legacy exact-at-lattice sampling (U_e showed
  // ~2e-3 relative on the CPW example, visualization only). Order 2p reproduces the
  // legacy values exactly on affine elements at ~3.5x the field storage, which may be
  // immaterial next to the solver memory high water mark.
  // All direct VTU point fields and the MFEM mesh writer must use one sampling lattice.
  // This matters in 2D magnetostatics, where A is ND order p while scalar B is L2 order
  // p-1. Use the maximum primary-field order rather than whichever field happens to be
  // initialized first.
  int paraview_refine_order = 0;
  for (const auto *field : {E.get(), B.get(), V.get(), A.get()})
  {
    if (field)
    {
      paraview_refine_order =
          std::max(paraview_refine_order, field->ParFESpace()->GetMaxElementOrder());
    }
  }
  MFEM_VERIFY(paraview_refine_order > 0, "Unable to determine field output order!");
  auto InitializeVizSpaces = [&](mfem::ParFiniteElementSpace &src_fespace)
  {
    if (viz_fec)
    {
      return;
    }
    auto *pmesh = src_fespace.GetParMesh();
    viz_fec =
        std::make_unique<mfem::L2_FECollection>(paraview_refine_order, pmesh->Dimension());
    viz_scalar_fespace =
        std::make_unique<mfem::ParFiniteElementSpace>(pmesh, viz_fec.get());
    viz_vector_fespace = std::make_unique<mfem::ParFiniteElementSpace>(
        pmesh, viz_fec.get(), pmesh->SpaceDimension());
  };
  const auto &output_mesh = fem_op->GetNDSpace().GetParMesh();
  const bool use_ceed_domain_paraview_fields =
      ShouldWriteParaviewFields() && UseCeedDomainParaviewPointFields(output_mesh);
  const bool use_ceed_domain_fields =
      ShouldWriteGridFunctionFields() || use_ceed_domain_paraview_fields;
  const bool use_ceed_boundary_fields =
      ShouldWriteParaviewFields() && UseCeedBoundaryParaviewPointFields(output_mesh);

  auto MakeFieldEvaluator = [&](PointFieldEvaluator::Kind kind,
                                mfem::ParFiniteElementSpace *e_fespace,
                                mfem::ParFiniteElementSpace *b_fespace, double scaling,
                                std::unique_ptr<PointFieldEvaluator> &eval,
                                std::unique_ptr<mfem::ParGridFunction> &gf)
  {
    InitializeVizSpaces(e_fespace ? *e_fespace : *b_fespace);
    auto &target = (kind == PointFieldEvaluator::Kind::POYNTING) ? *viz_vector_fespace
                                                                 : *viz_scalar_fespace;
    eval = std::make_unique<PointFieldEvaluator>(
        kind, fem_op->GetMaterialOp().GetMesh(), fem_op->GetMaterialOp(), e_fespace,
        b_fespace, target, scaling, ShouldWriteGridFunctionFields(),
        use_ceed_domain_paraview_fields);
    if (eval->IsValid())
    {
      if (ShouldWriteGridFunctionFields())
      {
        gf = std::make_unique<mfem::ParGridFunction>(&target);
        gf->UseDevice(true);
      }
    }
    else
    {
      RequireCeedPointFieldEvaluator(eval.get(), "domain field visualization");
    }
  };
  auto MakeBaseDomainFieldEvaluator = [&](PointFieldEvaluator::Kind kind,
                                          mfem::ParFiniteElementSpace &fespace,
                                          std::unique_ptr<PointFieldEvaluator> &eval)
  {
    if (!use_ceed_domain_paraview_fields)
    {
      return;
    }
    InitializeVizSpaces(fespace);
    const bool scalar_b_field = kind == PointFieldEvaluator::Kind::FIELD_B &&
                                fespace.GetParMesh()->Dimension() == 2;
    auto &target = (kind == PointFieldEvaluator::Kind::FIELD_H1 || scalar_b_field)
                       ? *viz_scalar_fespace
                       : *viz_vector_fespace;
    eval = std::make_unique<PointFieldEvaluator>(
        kind, fem_op->GetMaterialOp().GetMesh(), fem_op->GetMaterialOp(),
        (kind == PointFieldEvaluator::Kind::FIELD_E ||
         kind == PointFieldEvaluator::Kind::FIELD_H1)
            ? &fespace
            : nullptr,
        kind == PointFieldEvaluator::Kind::FIELD_B ? &fespace : nullptr, target, 1.0, false,
        true);
    if (!eval->IsValid())
    {
      RequireCeedPointFieldEvaluator(eval.get(), "base domain field visualization");
    }
  };
  auto MakeBdrSamplingPlan = [&]()
  {
    if (bdr_sampling_plan)
    {
      return;
    }
    const auto &mesh = fem_op->GetMaterialOp().GetMesh();
    const auto &pmesh = mesh.Get();
    const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    bdr_sampling_plan =
        std::make_shared<FaceSamplingPlan>(mesh, marker, paraview_refine_order);
    bdr_trace_cache =
        std::make_shared<BoundaryPhysicalTraceCache>(mesh, marker, bdr_sampling_plan);
  };
  auto MakeBdrFieldEvaluator = [&](PointFieldEvaluator::Kind kind,
                                   mfem::ParFiniteElementSpace &fespace,
                                   std::unique_ptr<PointFieldEvaluator> &eval)
  {
    const auto &mesh = fem_op->GetMaterialOp().GetMesh();
    const auto &pmesh = mesh.Get();
    const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    MakeBdrSamplingPlan();
    eval = std::make_unique<PointFieldEvaluator>(kind, mesh, marker, fespace,
                                                 paraview_refine_order, bdr_sampling_plan,
                                                 bdr_trace_cache, bdr_derived_bundle);
    if (!eval->IsValid())
    {
      RequireCeedPointFieldEvaluator(eval.get(), "boundary field visualization");
    }
  };
  auto MakeBdrCoeffEvaluator = [&](PointFieldEvaluator::Kind kind,
                                   mfem::ParFiniteElementSpace &fespace, double scaling,
                                   std::unique_ptr<PointFieldEvaluator> &eval)
  {
    const auto &mesh = fem_op->GetMaterialOp().GetMesh();
    const auto &pmesh = mesh.Get();
    const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    MakeBdrSamplingPlan();
    eval = std::make_unique<PointFieldEvaluator>(
        kind, mesh, marker, fespace, fem_op->GetMaterialOp(), paraview_refine_order,
        scaling, bdr_sampling_plan, bdr_trace_cache, bdr_derived_bundle);
    if (!eval->IsValid())
    {
      RequireCeedPointFieldEvaluator(eval.get(), "boundary coefficient visualization");
    }
  };
  auto MakeBdrPoyntingEvaluator =
      [&](mfem::ParFiniteElementSpace &e_fespace, mfem::ParFiniteElementSpace &b_fespace,
          double scaling, std::unique_ptr<PointFieldEvaluator> &eval)
  {
    const auto &mesh = fem_op->GetMaterialOp().GetMesh();
    const auto &pmesh = mesh.Get();
    const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> marker(bdr_attr_max);
    marker = 1;
    MakeBdrSamplingPlan();
    eval = std::make_unique<PointFieldEvaluator>(
        PointFieldEvaluator::Kind::POYNTING, mesh, marker, e_fespace, b_fespace,
        fem_op->GetMaterialOp(), paraview_refine_order, scaling, bdr_sampling_plan,
        bdr_trace_cache, bdr_derived_bundle);
    if (!eval->IsValid())
    {
      RequireCeedPointFieldEvaluator(eval.get(), "boundary Poynting visualization");
    }
  };

  // The fixed boundary lattice is independent of the field/QFunction. Build it once
  // after determining the ParaView refinement order, then hand the same immutable plan
  // to every boundary evaluator below.
  if (use_ceed_boundary_fields)
  {
    MakeBdrSamplingPlan();
    // Only the 3D vector E/B family is fused. The 2D scalar integral-mapped B path
    // deliberately retains the independent trace-backed kernels.
    if constexpr (HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
    {
      if (output_mesh.Dimension() == 3 && output_mesh.SpaceDimension() == 3)
      {
        const double electric_scaling =
            units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) /
            units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
        const double magnetic_scaling =
            units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
            units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);
        const auto &mesh = fem_op->GetMaterialOp().GetMesh();
        const auto &pmesh = mesh.Get();
        const int bdr_attr_max =
            pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
        mfem::Array<int> marker(bdr_attr_max);
        marker = 1;
        bdr_derived_bundle = std::make_shared<BoundaryDerivedFieldBundle>(
            mesh, marker, fem_op->GetMaterialOp(), *E->ParFESpace(), *B->ParFESpace(),
            bdr_sampling_plan, bdr_trace_cache, *E, *B, electric_scaling, magnetic_scaling);
        if (!bdr_derived_bundle->IsValid())
        {
          bdr_derived_bundle.reset();
        }
      }
    }
  }

  // Set-up evaluators for the paraview output / measurement.
  if (En)
  {
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_H1, *En->ParFESpace(),
                                 En_domain_eval);
  }

  if (Bt_inplane)
  {
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_E,
                                 *Bt_inplane->ParFESpace(), Bt_domain_eval);
  }

  if (E && Bt_inplane && use_ceed_domain_fields)
  {
    MakeFieldEvaluator(PointFieldEvaluator::Kind::MODE_SN, E->ParFESpace(),
                       Bt_inplane->ParFESpace(), 1.0, Sn_eval, Sn_gf);
  }

  if constexpr (HasVGridFunction<solver_t>())
  {
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_H1, *V->ParFESpace(),
                                 V_domain_eval);
    if (use_ceed_boundary_fields)
    {
      MakeBdrFieldEvaluator(PointFieldEvaluator::Kind::FIELD_H1, *V->ParFESpace(),
                            V_bdr_eval);
    }
  }

  if constexpr (HasAGridFunction<solver_t>())
  {
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_E, *A->ParFESpace(),
                                 A_domain_eval);
    if (use_ceed_boundary_fields)
    {
      MakeBdrFieldEvaluator(PointFieldEvaluator::Kind::FIELD_E, *A->ParFESpace(),
                            A_bdr_eval);
    }
  }

  if constexpr (HasEGridFunction<solver_t>())
  {
    // If E is dimensionalized when the coefficients are evaluated, the scaling only needs
    // to account for the remaining ε_0 = D / E. This assumes ProjectCoefficient(),
    // ProjectBdrCoefficient(), or paraview->Save() for U_e, Q_sr, and Q_si are always
    // called after the E GridFunction has been dimensionalized. To output nondimensional
    // coefficients, omit the scaling argument and make sure E is nondimensional when the
    // coefficients are evaluated.
    const double scaling = units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) /
                           units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);

    // Electric Energy Density.
    // U_e = 1/2 Dᴴ E = 1/2 ε_0 Eᴴ E.
    U_e = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(
        *E, fem_op->GetMaterialOp(), scaling);
    if (use_ceed_domain_fields)
    {
      MakeFieldEvaluator(PointFieldEvaluator::Kind::ENERGY_E, E->ParFESpace(), nullptr,
                         scaling, U_e_eval, U_e_gf);
    }
    if (use_ceed_boundary_fields)
    {
      MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::ENERGY_E, *E->ParFESpace(), scaling,
                            Ue_bdr_eval);
    }

    // Electric Boundary Field & Surface Charge.
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_E, *E->ParFESpace(),
                                 E_domain_eval);
    if (use_ceed_boundary_fields)
    {
      MakeBdrFieldEvaluator(PointFieldEvaluator::Kind::FIELD_E, *E->ParFESpace(),
                            E_bdr_eval);
    }
    // Q_s = D ⋅ n = ε_0 ε E ⋅ n, with the complex permittivity ε = Re{ε} + i Im{ε} of a
    // lossy dielectric for complex fields (the fused bundle applies it in-kernel; the
    // independent evaluators need the Im{ε}-weighted counterpart, combined at write time).
    if (use_ceed_boundary_fields)
    {
      MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::FLUX_Q, *E->ParFESpace(), scaling,
                            Q_bdr_eval);
      if (HasComplexGridFunction<solver_t>() && fem_op->GetMaterialOp().HasLossTangent() &&
          !bdr_derived_bundle)
      {
        const auto &mesh = fem_op->GetMaterialOp().GetMesh();
        const auto &pmesh = mesh.Get();
        const int bdr_attr_max =
            pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
        mfem::Array<int> marker(bdr_attr_max);
        marker = 1;
        Q_bdr_eval_imag = std::make_unique<PointFieldEvaluator>(
            PointFieldEvaluator::Kind::FLUX_Q, mesh, marker, *E->ParFESpace(),
            fem_op->GetMaterialOp(), paraview_refine_order, scaling, bdr_sampling_plan,
            bdr_trace_cache, nullptr, /*imag_permittivity*/ true);
        RequireCeedPointFieldEvaluator(Q_bdr_eval_imag.get(),
                                       "lossy boundary surface charge visualization");
      }
    }
  }

  if constexpr (HasBGridFunction<solver_t>())
  {
    // If B is dimensionalized when the coefficients are evaluated, the scaling only needs
    // to account for the remaining μ⁻¹ = H / B. This assumes ProjectCoefficient(),
    // ProjectBdrCoefficient(), or paraview->Save() for U_m, J_sr, and J_si are always
    // called after the B GridFunction has been dimensionalized. To output nondimensional
    // coefficients, omit the scaling argument and make sure B is nondimensional when the
    // coefficients are evaluated.
    const double scaling = units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
                           units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);

    // Magnetic Energy Density.
    // U_m = 1/2 Hᴴ B = 1/2 μ⁻¹ Bᴴ B.
    U_m = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(
        *B, fem_op->GetMaterialOp(), scaling);
    if (use_ceed_domain_fields)
    {
      MakeFieldEvaluator(PointFieldEvaluator::Kind::ENERGY_M, nullptr, B->ParFESpace(),
                         scaling, U_m_eval, U_m_gf);
    }
    if (use_ceed_boundary_fields)
    {
      MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::ENERGY_M, *B->ParFESpace(), scaling,
                            Um_bdr_eval);
    }

    // Magnetic Field, Boundary Field & Surface Current. In 2D, B is scalar (L2), so the
    // domain field is scalar and there is no boundary B vector. The surface current
    // J_s = n x (H_z ẑ) remains a valid in-plane tangential vector.
    MakeBaseDomainFieldEvaluator(PointFieldEvaluator::Kind::FIELD_B, *B->ParFESpace(),
                                 B_domain_eval);
    if (B->Real().VectorDim() > 1 && use_ceed_boundary_fields)
    {
      MakeBdrFieldEvaluator(PointFieldEvaluator::Kind::FIELD_B, *B->ParFESpace(),
                            B_bdr_eval);
    }
    if (use_ceed_boundary_fields)
    {
      MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::CURRENT_J, *B->ParFESpace(), scaling,
                            J_bdr_eval);
    }
  }

  if constexpr (HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
  {
    if constexpr (solver_t != ProblemType::BOUNDARYMODE)
    {
      // 3D Poynting Vector: S = Re{E x H⋆} = Re{E x μ⁻¹B⋆}.
      // E and B will be dimensionalized when the coefficient is evaluated, so the scaling
      // only needs to account for the remaining μ⁻¹ = H / B. As mentioned above,
      // ProjectCoefficient(*S.get()) or paraview->Save() should always be called after
      // E and B have been dimensionalized.
      const double scaling = units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
                             units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);
      S = std::make_unique<PoyntingVectorCoefficient>(*E, *B, fem_op->GetMaterialOp(),
                                                      scaling);
      if (use_ceed_domain_fields)
      {
        MakeFieldEvaluator(PointFieldEvaluator::Kind::POYNTING, E->ParFESpace(),
                           B->ParFESpace(), scaling, S_eval, S_gf);
      }
      if (use_ceed_boundary_fields)
      {
        MakeBdrPoyntingEvaluator(*E->ParFESpace(), *B->ParFESpace(), scaling, S_bdr_eval);
      }
    }
    // For boundary mode, Sn = Re{Et · (ẑ × Ht*)} uses the in-plane B field and is
    // registered below through a domain point evaluator.
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::EnsureFieldCoefficientsSetup()
{
  if (!field_coefficients_initialized)
  {
    SetupFieldCoefficients();
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::EnsureParaviewDataCollection()
{
  EnsureFieldCoefficientsSetup();
  if (ShouldWriteParaviewFields() && !paraview)
  {
    InitializeParaviewDataCollection();
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::InitializeParaviewDataCollection(
    const fs::path &sub_folder_name)
{
  EnsureFieldCoefficientsSetup();
  if (!ShouldWriteParaviewFields())
  {
    return;
  }
  // Close any files held by the previous collection before switching output paths.
  paraview.reset();
  paraview_bdr.reset();

  auto paraview_root = post_dir / "paraview";
  fs::path paraview_dir_v = paraview_root / OutputFolderName(solver_t);
  fs::path paraview_dir_b =
      paraview_root / fmt::format("{}_boundary", OutputFolderName(solver_t));
  if (!sub_folder_name.empty())
  {
    paraview_dir_v /= sub_folder_name;
    paraview_dir_b /= sub_folder_name;
  }
  // Set up postprocessing for output to disk.
  paraview = {paraview_dir_v, &fem_op->GetNDSpace().GetParMesh()};
  paraview_bdr = {paraview_dir_b, &fem_op->GetNDSpace().GetParMesh()};
  surface_response_field_paraview_registered = false;
  surface_response_flux_paraview_registered = false;

  const mfem::VTKFormat format = mfem::VTKFormat::BINARY32;
  const bool use_ho = true;
  MFEM_VERIFY(viz_scalar_fespace,
              "ParaView point-field space must be initialized before the collection!");
  const int refine_ho = viz_scalar_fespace->GetMaxElementOrder();
  const bool use_ceed_domain_paraview =
      UseCeedDomainParaviewPointFields(fem_op->GetNDSpace().GetParMesh());
  const bool use_ceed_boundary_paraview =
      UseCeedBoundaryParaviewPointFields(fem_op->GetNDSpace().GetParMesh());

  // Output mesh coordinate units same as input.
  paraview->SetCycle(-1);
  paraview_max_cycle = -1;
  paraview->SetDataFormat(format);
  MFEM_VERIFY(use_ceed_domain_paraview,
              "ParaView domain output requires the libCEED point-field path!");
  // Direct point fields use appended raw VTU arrays so their payload can be written
  // without an extra full-field base64 staging buffer.
  paraview->SetCompressionLevel(0);
  paraview->SetHighOrderOutput(use_ho);
  paraview->SetLevelsOfDetail(refine_ho);

  paraview_bdr->SetBoundaryOutput(true);
  paraview_bdr->SetCycle(-1);
  paraview_bdr->SetDataFormat(format);
  MFEM_VERIFY(use_ceed_boundary_paraview,
              "ParaView boundary output requires the libCEED point-field path!");
  paraview_bdr->SetCompressionLevel(0);
  paraview_bdr->SetHighOrderOutput(use_ho);
  paraview_bdr->SetLevelsOfDetail(refine_ho);

  // Register libCEED domain and boundary visualization fields lazily: ParaView Save()
  // evaluates one field at a time into a reusable device-capable buffer and then writes
  // that buffer to the VTU file. Domain direct-VTU buffers are point-major to match VTK
  // tuple order; boundary buffers keep the existing component-major packing path. This
  // keeps peak memory proportional to one output field instead of the sum of all derived
  // fields, while still doing the sampling/transforms in libCEED. Nonconforming AMR
  // meshes use the same point-field path; the static mesh/header caches above avoid
  // repeating VTU setup across mode saves.
  auto RegisterDomainEvalField =
      [&](const std::string &name, const std::unique_ptr<PointFieldEvaluator> &eval,
          const GridFunction *E_field, const GridFunction *B_field)
  {
    const auto *eval_ptr = eval.get();
    const auto *E_ptr = E_field;
    const auto *B_ptr = B_field;
    paraview->RegisterDomainPointEvaluator(
        name, [eval_ptr, E_ptr, B_ptr](Vector &buffer)
        { eval_ptr->EvalBuffer(E_ptr, B_ptr, buffer); }, eval_ptr->BufferBases(),
        eval_ptr->BufferNumComp(), eval_ptr->BufferSize(), true);
  };
  auto RegisterDomainBaseEvalField = [&](const std::string &name,
                                         const std::unique_ptr<PointFieldEvaluator> &eval,
                                         const mfem::ParGridFunction &field)
  {
    RequireCeedPointFieldEvaluator(eval.get(),
                                   fmt::format("domain {} visualization", name));
    const auto *eval_ptr = eval.get();
    const auto *field_ptr = &field;
    paraview->RegisterDomainPointEvaluator(
        name,
        [eval_ptr, field_ptr](Vector &buffer) { eval_ptr->EvalBuffer(*field_ptr, buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize(), true);
  };
  auto RegisterBdrEvalField = [&](const std::string &name,
                                  const std::unique_ptr<PointFieldEvaluator> &eval,
                                  const auto &field)
  {
    const auto *eval_ptr = eval.get();
    const auto *field_ptr = &field;
    paraview_bdr->RegisterBoundaryPointEvaluator(
        name,
        [eval_ptr, field_ptr](Vector &buffer) { eval_ptr->EvalBuffer(*field_ptr, buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize());
  };
  auto RegisterBdrPoyntingField = [&](const std::string &name)
  {
    const auto *eval_ptr = S_bdr_eval.get();
    const auto *E_ptr = E.get();
    const auto *B_ptr = B.get();
    paraview_bdr->RegisterBoundaryPointEvaluator(
        name, [eval_ptr, E_ptr, B_ptr](Vector &buffer)
        { eval_ptr->EvalBuffer(E_ptr, B_ptr, buffer); }, eval_ptr->BufferBases(),
        eval_ptr->BufferNumComp(), eval_ptr->BufferSize());
  };

  // The registered 3D E/B bundle is the sole provider for compatible boundary output
  // fields. Keep all other boundary fields on their established individual evaluator
  // path; this includes H1/V/A, 2D scalar B, and arbitrary-source fallbacks.
  std::vector<CeedParaViewDataCollection::BoundaryPointFieldDescriptor> bdr_bundle_fields;
  auto RegisterBdrBundleField =
      [&](const std::string &name, PointFieldEvaluator::Kind kind, bool imaginary_phase)
  {
    if (!bdr_derived_bundle)
    {
      return false;
    }
    switch (kind)
    {
      case PointFieldEvaluator::Kind::FIELD_E:
      case PointFieldEvaluator::Kind::FIELD_B:
      case PointFieldEvaluator::Kind::FLUX_Q:
      case PointFieldEvaluator::Kind::CURRENT_J:
      case PointFieldEvaluator::Kind::ENERGY_E:
      case PointFieldEvaluator::Kind::ENERGY_M:
      case PointFieldEvaluator::Kind::POYNTING:
        break;
      default:
        return false;
    }
    const auto slice = bdr_derived_bundle->RegisteredSlice(kind, imaginary_phase);
    bdr_bundle_fields.push_back({name, slice.offset, slice.size, slice.num_comp});
    return true;
  };

  // Output fields @ phase = 0 and π/2 for frequency domain (rather than, for example,
  // peak phasors or magnitude = sqrt(2) * RMS). Boundary fields use the documented
  // trace averages or oriented jumps on internal interfaces.
  if (E)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      RegisterDomainBaseEvalField("E_real", E_domain_eval, E->Real());
      RegisterDomainBaseEvalField("E_imag", E_domain_eval, E->Imag());
      RequireCeedPointFieldEvaluator(E_bdr_eval.get(), "boundary E visualization");
      if (!RegisterBdrBundleField("E_real", PointFieldEvaluator::Kind::FIELD_E, false))
      {
        RegisterBdrEvalField("E_real", E_bdr_eval, E->Real());
      }
      if (!RegisterBdrBundleField("E_imag", PointFieldEvaluator::Kind::FIELD_E, true))
      {
        RegisterBdrEvalField("E_imag", E_bdr_eval, E->Imag());
      }
    }
    else
    {
      RegisterDomainBaseEvalField("E", E_domain_eval, E->Real());
      RequireCeedPointFieldEvaluator(E_bdr_eval.get(), "boundary E visualization");
      if (!RegisterBdrBundleField("E", PointFieldEvaluator::Kind::FIELD_E, false))
      {
        RegisterBdrEvalField("E", E_bdr_eval, E->Real());
      }
    }
  }
  if (En)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      RegisterDomainBaseEvalField("En_real", En_domain_eval, En->Real());
      RegisterDomainBaseEvalField("En_imag", En_domain_eval, En->Imag());
    }
    else
    {
      RegisterDomainBaseEvalField("En", En_domain_eval, En->Real());
    }
  }
  if (B)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      RegisterDomainBaseEvalField("B_real", B_domain_eval, B->Real());
      RegisterDomainBaseEvalField("B_imag", B_domain_eval, B->Imag());
      if (B_bdr_eval)
      {
        if (!RegisterBdrBundleField("B_real", PointFieldEvaluator::Kind::FIELD_B, false))
        {
          RegisterBdrEvalField("B_real", B_bdr_eval, B->Real());
        }
        if (!RegisterBdrBundleField("B_imag", PointFieldEvaluator::Kind::FIELD_B, true))
        {
          RegisterBdrEvalField("B_imag", B_bdr_eval, B->Imag());
        }
      }
    }
    else
    {
      RegisterDomainBaseEvalField("B", B_domain_eval, B->Real());
      if (B_bdr_eval)
      {
        if (!RegisterBdrBundleField("B", PointFieldEvaluator::Kind::FIELD_B, false))
        {
          RegisterBdrEvalField("B", B_bdr_eval, B->Real());
        }
      }
    }
  }
  // In-plane B field for mode analysis (2D vector on ND space).
  if (Bt_inplane)
  {
    RegisterDomainBaseEvalField("Bt_real", Bt_domain_eval, Bt_inplane->Real());
    RegisterDomainBaseEvalField("Bt_imag", Bt_domain_eval, Bt_inplane->Imag());
  }
  if (V)
  {
    RegisterDomainBaseEvalField("V", V_domain_eval, V->Real());
    RequireCeedPointFieldEvaluator(V_bdr_eval.get(), "boundary V visualization");
    RegisterBdrEvalField("V", V_bdr_eval, V->Real());
  }
  if (A)
  {
    RegisterDomainBaseEvalField("A", A_domain_eval, A->Real());
    RequireCeedPointFieldEvaluator(A_bdr_eval.get(), "boundary A visualization");
    RegisterBdrEvalField("A", A_bdr_eval, A->Real());
  }

  // Extract energy density field for electric field energy 1/2 Dᴴ E or magnetic field
  // energy 1/2 Hᴴ B. Also Poynting vector S = E x H⋆.
  if (U_e)
  {
    RequireCeedPointFieldEvaluator(U_e_eval.get(), "domain electric energy visualization");
    RegisterDomainEvalField("U_e", U_e_eval, E.get(), nullptr);
    RequireCeedPointFieldEvaluator(Ue_bdr_eval.get(),
                                   "boundary electric energy visualization");
    if (!RegisterBdrBundleField("U_e", PointFieldEvaluator::Kind::ENERGY_E, false))
    {
      RegisterBdrEvalField("U_e", Ue_bdr_eval, *E);
    }
  }
  if (U_m)
  {
    RequireCeedPointFieldEvaluator(U_m_eval.get(), "domain magnetic energy visualization");
    RegisterDomainEvalField("U_m", U_m_eval, nullptr, B.get());
    RequireCeedPointFieldEvaluator(Um_bdr_eval.get(),
                                   "boundary magnetic energy visualization");
    if (!RegisterBdrBundleField("U_m", PointFieldEvaluator::Kind::ENERGY_M, false))
    {
      RegisterBdrEvalField("U_m", Um_bdr_eval, *B);
    }
  }
  if (S)
  {
    RequireCeedPointFieldEvaluator(S_eval.get(), "domain Poynting visualization");
    RegisterDomainEvalField("S", S_eval, E.get(), B.get());
    RequireCeedPointFieldEvaluator(S_bdr_eval.get(), "boundary Poynting visualization");
    if (!RegisterBdrBundleField("S", PointFieldEvaluator::Kind::POYNTING, false))
    {
      RegisterBdrPoyntingField("S");
    }
  }

  // Extract surface charge from normally discontinuous ND E-field. Also extract surface
  // currents from tangentially discontinuous RT B-field The surface charge and surface
  // currents are single-valued at internal boundaries.
  if (Q_bdr_eval)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      // Independent evaluators: Q_r = Re{ε} E_r ⋅ n - Im{ε} E_i ⋅ n and
      // Q_i = Re{ε} E_i ⋅ n + Im{ε} E_r ⋅ n (the Im{ε} evaluator exists only for lossy
      // materials; the fused bundle computes the same combination in-kernel).
      auto RegisterBdrChargeField = [&](const std::string &name, bool imaginary_phase)
      {
        const auto *eval_ptr = Q_bdr_eval.get();
        const auto *eval_imag_ptr = Q_bdr_eval_imag.get();
        const auto *E_ptr = E.get();
        paraview_bdr->RegisterBoundaryPointEvaluator(
            name,
            [eval_ptr, eval_imag_ptr, E_ptr, imaginary_phase](Vector &buffer)
            {
              eval_ptr->EvalBuffer(imaginary_phase ? E_ptr->Imag() : E_ptr->Real(), buffer);
              if (eval_imag_ptr)
              {
                Vector cross(buffer.Size());
                cross.UseDevice(buffer.UseDevice());
                eval_imag_ptr->EvalBuffer(imaginary_phase ? E_ptr->Real() : E_ptr->Imag(),
                                          cross);
                buffer.Add(imaginary_phase ? 1.0 : -1.0, cross);
              }
            },
            eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize());
      };
      if (!RegisterBdrBundleField("Q_s_real", PointFieldEvaluator::Kind::FLUX_Q, false))
      {
        RegisterBdrChargeField("Q_s_real", false);
      }
      if (!RegisterBdrBundleField("Q_s_imag", PointFieldEvaluator::Kind::FLUX_Q, true))
      {
        RegisterBdrChargeField("Q_s_imag", true);
      }
    }
    else
    {
      if (!RegisterBdrBundleField("Q_s", PointFieldEvaluator::Kind::FLUX_Q, false))
      {
        RegisterBdrEvalField("Q_s", Q_bdr_eval, E->Real());
      }
    }
  }
  else if (Q_sr)
  {
    MFEM_ABORT(
        "Boundary surface-charge ParaView output requires a libCEED point evaluator!");
  }
  if (J_bdr_eval)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      if (!RegisterBdrBundleField("J_s_real", PointFieldEvaluator::Kind::CURRENT_J, false))
      {
        RegisterBdrEvalField("J_s_real", J_bdr_eval, B->Real());
      }
      if (!RegisterBdrBundleField("J_s_imag", PointFieldEvaluator::Kind::CURRENT_J, true))
      {
        RegisterBdrEvalField("J_s_imag", J_bdr_eval, B->Imag());
      }
    }
    else
    {
      if (!RegisterBdrBundleField("J_s", PointFieldEvaluator::Kind::CURRENT_J, false))
      {
        RegisterBdrEvalField("J_s", J_bdr_eval, B->Real());
      }
    }
  }
  else if (J_sr)
  {
    MFEM_ABORT(
        "Boundary surface-current ParaView output requires a libCEED point evaluator!");
  }

  if (Sn_eval)
  {
    RequireCeedPointFieldEvaluator(Sn_eval.get(), "boundary-mode Sn visualization");
    RegisterDomainEvalField("Sn", Sn_eval, E.get(), Bt_inplane.get());
  }

  if (!bdr_bundle_fields.empty())
  {
    // Capture the shared bundle so deregistering/re-registering collection fields cannot
    // leave a writer callback with a dangling packed-vector provider.
    const auto bundle = bdr_derived_bundle;
    paraview_bdr->RegisterBoundaryPointFieldBatch(
        bdr_bundle_fields, bundle->BufferBases(), [bundle]() -> const Vector &
        { return bundle->GetRegisteredPacked(); }, bundle->RegisteredPackedSize());
  }

  // Add wave port boundary mode postprocessing when available.
  for (const auto &[idx, data] : port_E0)
  {
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_real", idx), data.E0r.get());
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_imag", idx), data.E0i.get());
  }
  RegisterSurfaceResponseParaviewFields();
}

template <ProblemType solver_t>
std::unique_ptr<GridFunction>
PostOperator<solver_t>::MakeRecoveredFluxGridFunction(bool complex)
{
  auto &rt_fespace = fem_op->GetRTSpace();
  if (!flux_space)
  {
    return std::make_unique<GridFunction>(rt_fespace, complex);
  }
  if (!flux_sides_space)
  {
    flux_sides_fec = std::make_unique<DiscontinuousRTFECollection>(
        rt_fespace.GetFEColl(), rt_fespace.GetMesh().Dimension());
    flux_sides_space =
        std::make_unique<FiniteElementSpace>(rt_fespace.GetMesh(), flux_sides_fec.get());
  }
  return std::make_unique<GridFunction>(*flux_sides_space, complex);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::SetRecoveredFluxGridFunction(const Vector &d,
                                                          mfem::ParGridFunction &flux) const
{
  MFEM_VERIFY(d.Size() == GetRecoveredElectricFluxTrueVSize(), kRecoveredFluxSpaceMessage);
  if (!flux_space)
  {
    flux.SetFromTrueDofs(d);
    flux.ExchangeFaceNbrData();
    return;
  }
  // Each element takes the values of the broken L-vector at its signed element DOFs (the
  // versions of its side for split DOFs), in its native basis (undoing the DOF
  // transformation of the RT space); the other side of a sheet face on another process
  // arrives with the face-neighbour data.
  Vector l(flux_space->GetVSize());
  flux_space->GetProlongationMatrix()->Mult(d, l);
  const auto *hl = l.HostRead();
  auto *hf = flux.HostWrite();
  mfem::Array<int> dofs, element_dofs;
  mfem::DofTransformation dof_trans;
  mfem::Vector values;
  for (int e = 0; e < flux_space->GetMesh().GetNE(); e++)
  {
    flux_space->GetElementDofs(e, dofs, dof_trans);
    values.SetSize(dofs.Size());
    for (int i = 0; i < dofs.Size(); i++)
    {
      const int dof = dofs[i];
      values(i) = (dof >= 0) ? hl[dof] : -hl[-1 - dof];
    }
    if (dof_trans.GetDofTransformation())
    {
      dof_trans.InvTransformPrimal(values);
    }
    flux.ParFESpace()->GetElementDofs(e, element_dofs);
    MFEM_ASSERT(element_dofs.Size() == values.Size(),
                "Unexpected element DOF count of the per-side recovered flux!");
    for (int i = 0; i < element_dofs.Size(); i++)
    {
      hf[element_dofs[i]] = values(i);
    }
  }
  flux.ExchangeFaceNbrData();
}

template <ProblemType solver_t>
void PostOperator<solver_t>::RegisterSurfaceResponseParaviewFields()
{
  if (!paraview)
  {
    return;
  }
  if (has_surface_response_corrected_field && !surface_response_field_paraview_registered)
  {
    paraview->RegisterField("E_corrected_real", &surface_response_corrected_field->Real());
    paraview->RegisterField("E_corrected_imag", &surface_response_corrected_field->Imag());
    surface_response_field_paraview_registered = true;
  }
  if (has_surface_response_corrected_flux && !surface_response_flux_paraview_registered)
  {
    paraview->RegisterField("D_corrected_real", &surface_response_corrected_flux->Real());
    paraview->RegisterField("D_corrected_imag", &surface_response_corrected_flux->Imag());
    surface_response_flux_paraview_registered = true;
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::ClearSurfaceResponseCorrectedField()
{
  if (paraview && surface_response_field_paraview_registered)
  {
    paraview->DeregisterField("E_corrected_real");
    paraview->DeregisterField("E_corrected_imag");
  }
  if (paraview && surface_response_flux_paraview_registered)
  {
    paraview->DeregisterField("D_corrected_real");
    paraview->DeregisterField("D_corrected_imag");
  }
  has_surface_response_corrected_field = false;
  has_surface_response_corrected_flux = false;
  surface_response_field_paraview_registered = false;
  surface_response_flux_paraview_registered = false;
  surface_response_corrected_frequency.reset();
  surface_response_corrected_mode_overlap.reset();
}

namespace
{

// Exponent p of the factor L^p which undoes the effect of scaling the mesh by L on the
// values of a field mapped from the reference element. The mapping is set by the map type
// of the finite element space: the identity for VALUE (H1, or L2 with VALUE map type),
// J^-T for H(curl), J / |J| for H(div), and 1 / |J| for L2 with INTEGRAL map type, where
// the Jacobian J scales as L and its determinant |J| as L^dim.
int PiolaExponent(const mfem::ParFiniteElementSpace &fespace)
{
  const int dim = fespace.GetParMesh()->Dimension();
  switch (fespace.FEColl()->GetMapType(dim))
  {
    case mfem::FiniteElement::VALUE:
      return 0;
    case mfem::FiniteElement::H_CURL:
      return 1;
    case mfem::FiniteElement::H_DIV:
      return dim - 1;
    case mfem::FiniteElement::INTEGRAL:
      return dim;
    default:
      MFEM_ABORT("Unsupported finite element map type for mesh scaling!");
      return 0;
  }
}

}  // namespace

void ScaleGridFunction(double L, GridFunction &u)
{
  const int p = PiolaExponent(*u.ParFESpace());
  if (p != 0)
  {
    u *= std::pow(L, p);
  }
}

void ScaleGridFunctions(double L, std::initializer_list<GridFunction *> fields)
{
  for (auto *u : fields)
  {
    if (u)
    {
      ScaleGridFunction(L, *u);
    }
  }
}

void DimensionalizeGridFunctions(Units &units, std::unique_ptr<GridFunction> &E,
                                 std::unique_ptr<GridFunction> &B,
                                 std::unique_ptr<GridFunction> &V,
                                 std::unique_ptr<GridFunction> &A)
{
  if (E)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_E>(*E);
  }
  if (B)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_B>(*B);
  }
  if (A)
  {
    units.DimensionalizeInPlace<Units::ValueType::CURRENT>(*A);
  }
  if (V)
  {
    units.DimensionalizeInPlace<Units::ValueType::VOLTAGE>(*V);
  }
}

void NondimensionalizeGridFunctions(Units &units, std::unique_ptr<GridFunction> &E,
                                    std::unique_ptr<GridFunction> &B,
                                    std::unique_ptr<GridFunction> &V,
                                    std::unique_ptr<GridFunction> &A)
{
  if (E)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(*E);
  }
  if (B)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_B>(*B);
  }
  if (A)
  {
    units.NondimensionalizeInPlace<Units::ValueType::CURRENT>(*A);
  }
  if (V)
  {
    units.NondimensionalizeInPlace<Units::ValueType::VOLTAGE>(*V);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteParaviewFields(double time, int step)
{
  BlockTimer bt(Timer::POSTPRO_PARAVIEW);
  EnsureParaviewDataCollection();

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization. Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  DimensionalizeGridFunctions(units, E, B, V, A);
  if (has_surface_response_corrected_field)
  {
    ScaleSurfaceResponseGridFunctions(
        mesh_Lc0, mesh.Dimension(), surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
    DimensionalizeSurfaceResponseGridFunctions(
        units, surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
  }
  if (En)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_E>(*En);
  }
  if (Bt_inplane)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_B>(*Bt_inplane);
  }
  // libCEED-evaluated derived domain fields are registered as lazy point evaluators;
  // Save() evaluates them now, after dimensionalization, without materializing
  // GridFunctions. A solution Vector can keep its address across writes, so explicitly
  // advance the boundary trace save generation before any lazy evaluator asks for it.
  if (bdr_derived_bundle)
  {
    bdr_derived_bundle->Invalidate();
  }
  if (bdr_trace_cache)
  {
    bdr_trace_cache->Invalidate();
  }
  double paraview_time = time;
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    paraview_time =
        units.Dimensionalize<Units::ValueType::FREQUENCY>(time) / (2.0 * std::numbers::pi);
  }
  paraview->SetCycle(step);
  paraview_max_cycle = std::max(paraview_max_cycle, step);
  paraview->SetTime(paraview_time);
  paraview_bdr->SetCycle(step);
  paraview_bdr->SetTime(paraview_time);
  paraview->Save();
  paraview_bdr->Save();
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  NondimensionalizeGridFunctions(units, E, B, V, A);
  if (has_surface_response_corrected_field)
  {
    ScaleSurfaceResponseGridFunctions(
        1.0 / mesh_Lc0, mesh.Dimension(), surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
    NondimensionalizeSurfaceResponseGridFunctions(
        units, surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
  }
  if (En)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(*En);
  }
  if (Bt_inplane)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_B>(*Bt_inplane);
  }
  Mpi::Barrier(fem_op->GetComm());
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteParaviewFieldsFinal(const ErrorIndicator *indicator)
{
  BlockTimer bt(Timer::POSTPRO_PARAVIEW);
  EnsureParaviewDataCollection();

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Write the mesh partitioning and (optionally) error indicators at the final step. No
  // need for these to be parallel objects, since the data is local to each process and
  // there isn't a need to ever access the element neighbors. We set the time to some
  // non-used value to make the step identifiable within the data collection.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  paraview->SetCycle(paraview_max_cycle + 1);
  if (paraview->GetTime() < 1.0)
  {
    paraview->SetTime(99.0);
  }
  else
  {
    // 1 -> 99, 10 -> 999, etc.
    paraview->SetTime(
        std::pow(10.0, 2.0 + static_cast<int>(std::log10(paraview->GetTime()))) - 1.0);
  }
  auto field_map = paraview->GetFieldMap();  // Copy, so can reregister later
  for (const auto &[name, gf] : field_map)
  {
    paraview->DeregisterField(name);
  }
  auto coeff_field_map = paraview->GetCoeffFieldMap();
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview->DeregisterCoeffField(name);
  }
  auto vcoeff_field_map = paraview->GetVCoeffFieldMap();
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview->DeregisterVCoeffField(name);
  }
  Vector rank(mesh.GetNE());
  rank = mesh.GetMyRank() + 1;
  paraview->RegisterDomainCellField("Rank", rank);
  if (indicator)
  {
    MFEM_VERIFY(mesh.GetNE() == indicator->Local().Size(),
                "Size mismatch for provided ErrorIndicator for postprocessing!");
    paraview->RegisterDomainCellField("Indicator", indicator->Local());
  }
  paraview->DeregisterDomainPointFields();
  paraview->Save();
  paraview->DeregisterDomainCellField("Rank");
  if (indicator)
  {
    paraview->DeregisterDomainCellField("Indicator");
  }
  for (const auto &[name, gf] : field_map)
  {
    paraview->RegisterField(name, gf);
  }
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview->RegisterCoeffField(name, gf);
  }
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview->RegisterVCoeffField(name, gf);
  }
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  Mpi::Barrier(fem_op->GetComm());
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteMFEMGridFunctions(double time, int step)
{
  BlockTimer bt(Timer::POSTPRO_GRIDFUNCTION);
  EnsureFieldCoefficientsSetup();

  // Create output directory if it doesn't exist (node-local filesystem aware).
  EnsureDirectory(fs::path(gridfunction_output_dir), fem_op->GetComm());

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization. Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  DimensionalizeGridFunctions(units, E, B, V, A);
  if (has_surface_response_corrected_field)
  {
    ScaleSurfaceResponseGridFunctions(
        mesh_Lc0, mesh.Dimension(), surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
    DimensionalizeSurfaceResponseGridFunctions(
        units, surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
  }
  if (En)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_E>(*En);
  }
  if (Bt_inplane)
  {
    units.DimensionalizeInPlace<Units::ValueType::FIELD_B>(*Bt_inplane);
  }
  // Create grid function for vector coefficients.
  mfem::ParFiniteElementSpace &fespace = E ? *E->ParFESpace() : *B->ParFESpace();
  mfem::ParGridFunction gridfunc_vector(&fespace);

  // Create grid function for scalar coefficients.
  mfem::H1_FECollection h1_fec(fespace.GetMaxElementOrder(), mesh.Dimension());
  mfem::ParFiniteElementSpace h1_fespace(&mesh, &h1_fec);
  mfem::ParGridFunction gridfunc_scalar(&h1_fespace);

  const int local_rank = mesh.GetMyRank();

  auto write_grid_function = [&](const auto &gridfunc, const std::string &name)
  {
    auto path = fs::path(gridfunction_output_dir) /
                fmt::format("{}_{:0{}d}.gf.{:0{}d}", name, step, pad_digits_default,
                            local_rank, pad_digits_default);
    std::ofstream file(path);
    gridfunc.Save(file);
  };

  // Write grid functions using MFEM's built-in Save method.
  // Use 6-digit padding to match MFEM's pad_digits_default.
  if constexpr (HasEGridFunction<solver_t>())
  {
    if (E)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        // Write real and imaginary parts separately.
        write_grid_function(E->Real(), "E_real");
        write_grid_function(E->Imag(), "E_imag");
      }
      else
      {
        // Write real part only.
        write_grid_function(E->Real(), "E");
      }
    }
  }
  if (has_surface_response_corrected_field)
  {
    write_grid_function(surface_response_corrected_field->Real(), "E_corrected_real");
    write_grid_function(surface_response_corrected_field->Imag(), "E_corrected_imag");
    if (has_surface_response_corrected_flux)
    {
      write_grid_function(surface_response_corrected_flux->Real(), "D_corrected_real");
      write_grid_function(surface_response_corrected_flux->Imag(), "D_corrected_imag");
    }
  }

  if (En)
  {
    if constexpr (HasComplexGridFunction<solver_t>())
    {
      write_grid_function(En->Real(), "En_real");
      write_grid_function(En->Imag(), "En_imag");
    }
    else
    {
      write_grid_function(En->Real(), "En");
    }
  }

  if constexpr (HasBGridFunction<solver_t>())
  {
    if (B)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        write_grid_function(B->Real(), "B_real");
        write_grid_function(B->Imag(), "B_imag");
      }
      else
      {
        write_grid_function(B->Real(), "B");
      }
    }
  }

  if (Bt_inplane)
  {
    write_grid_function(Bt_inplane->Real(), "Bt_real");
    write_grid_function(Bt_inplane->Imag(), "Bt_imag");
  }

  if constexpr (HasVGridFunction<solver_t>())
  {
    if (V)
    {
      write_grid_function(V->Real(), "V");
    }
  }

  if constexpr (HasAGridFunction<solver_t>())
  {
    if (A)
    {
      write_grid_function(A->Real(), "A");
    }
  }

  if (U_e)
  {
    if (U_e_eval)
    {
      // Device-evaluated (libCEED) energy density into the interpolatory L2 output
      // space, reusing the ParaView evaluator instead of a host coefficient projection.
      U_e_eval->Eval(E.get(), nullptr, *U_e_gf);
      write_grid_function(*U_e_gf, "U_e");
    }
    else
    {
      gridfunc_scalar = 0.0;
      gridfunc_scalar.ProjectCoefficient(*U_e.get());
      write_grid_function(gridfunc_scalar, "U_e");
    }
  }

  if (U_m)
  {
    if (U_m_eval)
    {
      U_m_eval->Eval(nullptr, B.get(), *U_m_gf);
      write_grid_function(*U_m_gf, "U_m");
    }
    else
    {
      gridfunc_scalar = 0.0;
      gridfunc_scalar.ProjectCoefficient(*U_m.get());
      write_grid_function(gridfunc_scalar, "U_m");
    }
  }

  if (S)
  {
    if (S_eval)
    {
      S_eval->Eval(E.get(), B.get(), *S_gf);
      write_grid_function(*S_gf, "S");
    }
    else
    {
      gridfunc_vector = 0.0;
      gridfunc_vector.ProjectCoefficient(*S.get());
      write_grid_function(gridfunc_vector, "S");
    }
  }
  if (Sn_eval)
  {
    Sn_eval->Eval(E.get(), Bt_inplane.get(), *Sn_gf);
    write_grid_function(*Sn_gf, "Sn");
  }

  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  NondimensionalizeGridFunctions(units, E, B, V, A);
  if (has_surface_response_corrected_field)
  {
    ScaleSurfaceResponseGridFunctions(
        1.0 / mesh_Lc0, mesh.Dimension(), surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
    NondimensionalizeSurfaceResponseGridFunctions(
        units, surface_response_corrected_field.get(),
        has_surface_response_corrected_flux ? surface_response_corrected_flux.get()
                                            : nullptr);
  }
  if (En)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(*En);
  }
  if (Bt_inplane)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_B>(*Bt_inplane);
  }
  Mpi::Barrier(fem_op->GetComm());
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteMFEMGridFunctionsFinal(const ErrorIndicator *indicator)
{
  BlockTimer bt(Timer::POSTPRO_GRIDFUNCTION);
  EnsureFieldCoefficientsSetup();

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Write the mesh partitioning and (optionally) error indicators at the final step.
  // Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);

  // Create output directory if it doesn't exist (node-local filesystem aware).
  EnsureDirectory(fs::path(gridfunction_output_dir), fem_op->GetComm());

  // Create piecewise constant finite element space for rank and error indicator.
  mfem::L2_FECollection pwconst_fec(0, mesh.Dimension());
  mfem::FiniteElementSpace pwconst_fespace(&mesh, &pwconst_fec);

  const int local_rank = mesh.GetMyRank();

  auto write_grid_function = [&](const auto &gridfunc, const std::string &name)
  {
    auto path = fs::path(gridfunction_output_dir) /
                fmt::format("{}.gf.{:0{}d}", name, local_rank, pad_digits_default);
    std::ofstream file(path);
    gridfunc.Save(file);
  };

  // Write mesh partitioning (rank information).
  {
    mfem::GridFunction rank(&pwconst_fespace);
    rank = local_rank + 1;
    write_grid_function(rank, "rank");
  }

  // Write error indicator if provided.
  if (indicator)
  {
    mfem::GridFunction eta(&pwconst_fespace);
    MFEM_VERIFY(eta.Size() == indicator->Local().Size(),
                "Size mismatch for provided ErrorIndicator for postprocessing!");
    eta = indicator->Local();
    write_grid_function(eta, "indicator");
  }

  // Save ParMesh files; necessary to visualize grid functions.
  fs::path mesh_filename = fs::path(gridfunction_output_dir) / "mesh";
  mesh.Save(mesh_filename);

  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  Mpi::Barrier(fem_op->GetComm());
}

// Measurements.

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureDomainFieldEnergy() const
{
  measurement_cache.domain_E_field_energy_i.clear();
  measurement_cache.domain_H_field_energy_i.clear();

  measurement_cache.domain_E_field_energy_i.reserve(dom_post_op.M_i.size());
  measurement_cache.domain_H_field_energy_i.reserve(dom_post_op.M_i.size());

  if constexpr (HasEGridFunction<solver_t>())
  {
    // Use V if it has it rather than E.
    auto &field = V ? *V : *E;
    auto energy = dom_post_op.GetElectricFieldEnergy(field);
    measurement_cache.domain_E_field_energy_all = energy;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      auto energy_i = dom_post_op.GetDomainElectricFieldEnergy(idx, field);
      auto participation_ratio = std::abs(energy_i) > 0.0 ? energy_i / energy : 0.0;
      measurement_cache.domain_E_field_energy_i.emplace_back(
          Measurement::DomainData{idx, energy_i, participation_ratio});
    }
  }
  else
  {
    // Magnetic field only.
    measurement_cache.domain_E_field_energy_all = 0.0;
    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      measurement_cache.domain_E_field_energy_i.emplace_back(
          Measurement::DomainData{idx, 0.0, 0.0});
    }
  }

  if (HasBGridFunction<solver_t>())
  {
    // A broken magnetostatic solution has no grid function, but its L-vector.
    auto &field = A ? *A : *B;
    const bool use_l = (a_l.Size() > 0);
    auto energy = use_l ? dom_post_op.GetMagneticFieldEnergy(a_l, field.GetComm())
                        : dom_post_op.GetMagneticFieldEnergy(field);
    measurement_cache.domain_H_field_energy_all = energy;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      auto energy_i =
          use_l ? dom_post_op.GetDomainMagneticFieldEnergy(idx, a_l, field.GetComm())
                : dom_post_op.GetDomainMagneticFieldEnergy(idx, field);
      auto participation_ratio = std::abs(energy) > 0.0 ? energy_i / energy : 0.0;
      measurement_cache.domain_H_field_energy_i.emplace_back(
          Measurement::DomainData{idx, energy_i, participation_ratio});
    }
  }
  else
  {
    // Electric field only.
    measurement_cache.domain_H_field_energy_all = 0.0;
    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      measurement_cache.domain_H_field_energy_i.emplace_back(
          Measurement::DomainData{idx, 0.0, 0.0});
    }
  }

  // Log Domain Energy.
  const auto domain_E = units.Dimensionalize<Units::ValueType::ENERGY>(
      measurement_cache.domain_E_field_energy_all);
  const auto domain_H = units.Dimensionalize<Units::ValueType::ENERGY>(
      measurement_cache.domain_H_field_energy_all);
  if constexpr (HasEGridFunction<solver_t>() && !HasBGridFunction<solver_t>())
  {
    Mpi::Print(" Field energy E = {:.3e} J\n", domain_E);
  }
  else if constexpr (!HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
  {
    Mpi::Print(" Field energy H = {:.3e} J\n", domain_H);
  }
  else if constexpr (solver_t != ProblemType::EIGENMODE &&
                     solver_t != ProblemType::BOUNDARYMODE)
  {
    Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", domain_E, domain_H,
               domain_E + domain_H);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPorts() const
{
  measurement_cache.lumped_port_vi.clear();
  measurement_cache.lumped_port_inductor_energy = 0.0;
  measurement_cache.lumped_port_capacitor_energy = 0.0;

  if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN ||
                solver_t == ProblemType::TRANSIENT)
  {
    const auto port_powers = fem_op->GetLumpedPortOp().GetPowers(*E, *B);
    const auto port_voltages = fem_op->GetLumpedPortOp().GetVoltages(*E);
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      auto &vi = measurement_cache.lumped_port_vi[idx];
      vi.P = port_powers.at(idx);
      vi.V = port_voltages.at(idx);
      if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN)
      {
        // Compute current from the port impedance, separate contributions for R, L, C
        // branches.
        // Get value and make real: Matches current behaviour (even for eigensolver!).
        MFEM_VERIFY(
            measurement_cache.freq.real() > 0.0,
            "Frequency domain lumped port postprocessing requires nonzero frequency!");
        vi.I_RLC[0] =
            (std::abs(data.R) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::R)
                : 0.0;
        vi.I_RLC[1] =
            (std::abs(data.L) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::L)
                : 0.0;
        vi.I_RLC[2] =
            (std::abs(data.C) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::C)
                : 0.0;
        vi.I = std::accumulate(vi.I_RLC.begin(), vi.I_RLC.end(),
                               std::complex<double>{0.0, 0.0});
        vi.S = vi.V / std::sqrt(data.GetExcitationRefResistance());

        // Add contribution due to all inductive lumped boundaries in the model:
        //                      E_ind = ∑_j 1/2 L_j I_mj².
        if (std::abs(data.L) > 0.0)
        {
          std::complex<double> I_mj = vi.I_RLC[1];
          vi.inductor_energy = 0.5 * std::abs(data.L) * std::real(I_mj * std::conj(I_mj));
          measurement_cache.lumped_port_inductor_energy += vi.inductor_energy;
        }

        // Add contribution due to all capacitive lumped boundaries in the model:
        //                      E_cap = ∑_j 1/2 C_j V_mj².
        if (std::abs(data.C) > 0.0)
        {
          std::complex<double> V_mj = vi.V;
          vi.capacitor_energy = 0.5 * std::abs(data.C) * std::real(V_mj * std::conj(V_mj));
          measurement_cache.lumped_port_capacitor_energy += vi.capacitor_energy;
        }
      }
      else
      {
        // Compute current from P = V I^* since there is no frequency & characteristic
        // impedance of the lumped element.
        vi.I = (std::abs(vi.V) > 0.0) ? std::conj(vi.P / vi.V) : 0.0;
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPortsEig() const
{
  // Depends on MeasureLumpedPorts.
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    auto freq_re = measurement_cache.freq.real();
    auto energy_electric_all = measurement_cache.domain_E_field_energy_all +
                               measurement_cache.lumped_port_capacitor_energy;
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      // Get previously computed data: should never fail as defined by MeasureLumpedPorts.
      auto &vi = measurement_cache.lumped_port_vi.at(idx);

      // Resistive Lumped Ports:
      // Compute participation ratio of external ports (given as any port boundary with
      // nonzero resistance). Currently no reactance of the ports is supported. The
      // magnitude of the κ of the port is the decay rate due to the power P_mj dissipated
      // in the resistive sheet, and its sign is the sign of Re{I_mj}, which is not
      // meaningful when Re{I_mj} ≈ 0, for example when V_mj ≈ 0:
      //                   κ_mj = sign(Re{I_mj}) P_mj / E_m,
      // from which the mode coupling quality factor is computed as:
      //                         Q_mj = Re{ω_m} / |κ_mj|.
      // P_mj = 1/2 ∫ |E_t|² / R_s dS is the boundary term of the damping matrix, so for
      // active ports with R_j > 0, ∑_j |κ_mj| = 2 Im{ω_m} when the ports are the only
      // loss. The voltage estimate 1/2 R_j |I_mj|² = |V_mj|² / (2 R_j) is smaller when the
      // port field does not match the uniform or coaxial port mode, and zero when the mode
      // voltage cancels. Inactive ports are not in the damping matrix and report the loss
      // which the port termination would add.
      if (std::abs(data.R) > 0.0)
      {
        std::complex<double> I_mj = vi.I_RLC[0];
        const double resistor_power = data.GetResistorPower(fem_op->GetNDSpace(), *E);
        vi.mode_port_kappa =
            std::copysign(resistor_power / energy_electric_all, I_mj.real());
        vi.quality_factor = (vi.mode_port_kappa == 0.0)
                                ? mfem::infinity()
                                : freq_re / std::abs(vi.mode_port_kappa);
      }

      // Inductive Lumped Ports:
      // Compute energy-participation ratio of junction given by index idx for the field
      // mode. We first get the port line voltage, and use lumped port circuit impedance to
      // get peak current through the inductor: I_mj = V_mj / Z_mj,  Z_mj = i ω_m L_j. E_m
      // is the total energy in mode m: E_m = E_elec + E_cap = E_mag + E_ind. The signed EPR
      // for a lumped inductive element is computed as:
      //                            p_mj = 1/2 L_j I_mj² / E_m.
      // An element with no assigned inductance will be treated as having zero admittance
      // and thus zero current.
      if (std::abs(data.L) > 0.0)
      {
        std::complex<double> I_mj = vi.I_RLC[1];
        vi.inductive_energy_participation =
            std::copysign(vi.inductor_energy / energy_electric_all, I_mj.real());
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureWavePorts() const
{
  measurement_cache.wave_port_vi.clear();

  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      // Get value and make real: Matches current behaviour.
      auto freq_re = measurement_cache.freq.real();  // TODO: Fix
      MFEM_VERIFY(freq_re > 0.0,
                  "Frequency domain wave port postprocessing requires nonzero frequency!");
      auto &vi = measurement_cache.wave_port_vi[idx];
      vi.S = data.GetSParameter(*E);
      if (data.HasVoltageCoords())
      {
        vi.P = data.GetPower(*E, *B);
        vi.V = data.GetVoltage(*E);
        vi.Z_PV = data.GetCharacteristicImpedance();
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureFloquetPorts() const
{
  measurement_cache.floquet_port_s.clear();

  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    if (fem_op->GetFloquetPortOp().Empty())
    {
      return;
    }
    // Determine if the excitation uses circular polarization (for output basis).
    bool circular_output = false;
    for (const auto &[idx, port] : fem_op->GetFloquetPortOp())
    {
      if (port.excitation == measurement_cache.ex_idx)
      {
        circular_output = port.HasCircularExcitation();
        break;
      }
    }
    measurement_cache.floquet_circular_output = circular_output;
    for (const auto &[idx, data] : fem_op->GetFloquetPortOp())
    {
      bool is_driving = (data.excitation == measurement_cache.ex_idx);
      auto S_all = data.GetAllSParameters(*E, is_driving);

      // For circular excitation, rotate TE/TM output to RHC/LHC basis.
      // Unitary: |S_RHC|² + |S_LHC|² = |S_TE|² + |S_TM|² (energy conserved).
      if (circular_output)
      {
        const double inv_sqrt2 = 1.0 / std::numbers::sqrt2;
        std::map<std::tuple<int, int, bool>, std::complex<double>> circ;
        std::set<std::pair<int, int>> orders;
        for (const auto &[key, S] : S_all)
        {
          orders.insert({std::get<0>(key), std::get<1>(key)});
        }
        for (const auto &[m, n] : orders)
        {
          auto it_te = S_all.find({m, n, true});
          auto it_tm = S_all.find({m, n, false});
          auto s_te = (it_te != S_all.end()) ? it_te->second : 0.0;
          auto s_tm = (it_tm != S_all.end()) ? it_tm->second : 0.0;
          circ[{m, n, true}] = (s_te + 1i * s_tm) * inv_sqrt2;   // RHC
          circ[{m, n, false}] = (s_te - 1i * s_tm) * inv_sqrt2;  // LHC
        }
        S_all = std::move(circ);
      }

      measurement_cache.floquet_port_s[idx] = std::move(S_all);
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureSParameter() const
{
  // Depends on LumpedPorts, WavePorts, FloquetPorts.
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    using fmt::format;
    using std::complex_literals::operator""i;

    // Don't measure S-Matrix unless there is only one excitation per port. Mixed
    // lumped/wave/Floquet S-parameters are supported with a power normalization bridge:
    // lumped/wave ports normalize to unit peak power (∫E×H*·n dS = 1), while Floquet ports
    // normalize to unit time-averaged power (½Re{∫E×H*·n dS} = 1). The √2 bridge factor
    // accounts for this 2× difference in incident power convention.
    if (!fem_op->GetPortExcitations().IsMultipleSimple())
    {
      return;
    }

    // S-parameter computation with per-port reference impedances (Kurokawa power-wave
    // generalized S-parameters).
    //
    // Each port type uses its natural real reference impedance: lumped ports use
    // GetExcitationRefResistance() (the user-specified R, or the unit internal reference
    // for a purely reactive R == 0 port — any port reactance acts through the
    // system-matrix termination, not this normalization); wave ports implicitly use the
    // line's characteristic impedance encoded in the unit-power normalization of the
    // boundary mode (∫|E_mode×H_mode⋆|·n dS = 1).
    //
    // GetSParameter() returns the b-amplitude in Kurokawa convention for both port types:
    //   - Lumped: ∫E·(E_inc/R_ref) dS = V/V_inc, where V_inc encodes R_ref
    //   - Wave:   ∫(E×H_mode)·n dS / ∫(E_mode×H_mode)·n dS, the modal reaction overlap,
    //             which separates the incident and reflected modal amplitudes also for a
    //             lossy port mode, and is the power overlap ∫(E×H_mode⋆)·n dS for a
    //             lossless one
    // Both are directly the Kurokawa b/a ratio, so no impedance scaling is needed.
    //
    // Both measure the response to a unit-power incident wave. When the drive port carries
    // an excitation weight a = A exp(iφ) (lumped ports only in this version), the measured
    // fields respond to a times the unit-power wave, so the projections on every
    // observation port (lumped, wave, and Floquet) are divided by a. Then
    // S = proj / a - δ on the driven port and S = proj / a on the others, which is
    // independent of A and φ.

    // Get information about excited port.
    auto [drive_is_simple, drive_port_type, drive_port_idx] =
        fem_op->GetPortExcitations().excitations.at(measurement_cache.ex_idx).IsSimple();

    if (drive_port_type != PortType::LumpedPort && drive_port_type != PortType::WavePort &&
        drive_port_type != PortType::FloquetPort)
    {
      return;
    }

    // Drive weight a = A exp(iφ) of the excited port (1 for wave and Floquet ports, which
    // do not carry a weight in this version).
    const std::complex<double> drive_weight =
        (drive_port_type == PortType::LumpedPort)
            ? fem_op->GetLumpedPortOp().GetPort(drive_port_idx).GetExcitationWeight()
            : std::complex<double>{1.0, 0.0};

    // Per-side wave-port de-embedding factor: exp(i·kₙ·d_offset) when the port is a
    // wave port with a non-zero offset; 1 otherwise (lumped ports have no offset).
    // Applied independently on the source and observation sides, so cross-type S
    // (lumped↔wave) gets the correct single-sided phase.
    const std::complex<double> src_deembed =
        (drive_port_type == PortType::WavePort)
            ? std::exp(1i * fem_op->GetWavePortOp().GetPort(drive_port_idx).kn0 *
                       fem_op->GetWavePortOp().GetPort(drive_port_idx).d_offset)
            : std::complex<double>{1.0, 0.0};

    // Power normalization bridge factor for mixed Floquet + lumped/wave configurations.
    // Floquet ports normalize incident waves to unit time-averaged power (P_avg = 1),
    // while lumped/wave ports normalize to unit peak power (∫E×H*·n dS = 1, P_avg = ½).
    // Cross-type observations require a √2 correction:
    //   - Floquet drives, lumped/wave observes: lumped |b|² = 2×P_avg, divide by √2
    //   - Lumped/wave drives, Floquet observes: Floquet |S|² = P_avg, multiply by √2
    const double inv_sqrt2 = 1.0 / std::numbers::sqrt2;
    const bool floquet_drives = (drive_port_type == PortType::FloquetPort);
    const bool lumped_or_wave_drives =
        (drive_port_type == PortType::LumpedPort || drive_port_type == PortType::WavePort);

    // Iterate over observation lumped ports.
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      auto &vi = measurement_cache.lumped_port_vi.at(idx);
      vi.S /= drive_weight;
      if (drive_port_type == PortType::LumpedPort && idx == drive_port_idx)
      {
        vi.S.real(vi.S.real() - 1.0);
      }
      // Lumped observation has no d_offset — only the source-side factor applies.
      vi.S *= src_deembed;
      if (floquet_drives)
      {
        vi.S *= inv_sqrt2;
      }

      Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                 format("S[{}][{}]", idx, drive_port_idx), vi.S.real(), vi.S.imag(),
                 Measurement::Magnitude(vi.S), Measurement::Phase(vi.S));
    }

    // Iterate over observation wave ports.
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      auto &vi = measurement_cache.wave_port_vi.at(idx);
      vi.S /= drive_weight;
      if (drive_port_type == PortType::WavePort && idx == drive_port_idx)
      {
        vi.S.real(vi.S.real() - 1.0);
      }
      // Apply both source and observation de-embedding factors.
      vi.S *= src_deembed;
      vi.S *= std::exp(1i * data.kn0 * data.d_offset);
      if (floquet_drives)
      {
        vi.S *= inv_sqrt2;
      }

      Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                 format("S[{}][{}]", idx, drive_port_idx), vi.S.real(), vi.S.imag(),
                 Measurement::Magnitude(vi.S), Measurement::Phase(vi.S));
    }

    // Floquet port S-parameters (already post-processed in MeasureFloquetPorts). Apply
    // the reciprocal bridge factor when lumped/wave ports are driving.
    for (auto &[port_idx, S_map] : measurement_cache.floquet_port_s)
    {
      for (auto &[key, S] : S_map)
      {
        S /= drive_weight;
        if (lumped_or_wave_drives)
        {
          S *= std::numbers::sqrt2;
        }
        auto [m, n, is_te] = key;
        auto pol = measurement_cache.floquet_circular_output ? (is_te ? "RHC" : "LHC")
                                                             : (is_te ? "TE" : "TM");
        Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                   format("S[P{}({},{}){}][{}]", port_idx, m, n, pol, drive_port_idx),
                   S.real(), S.imag(), Measurement::Magnitude(S), Measurement::Phase(S));
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureSurfaceFlux() const
{
  // Compute the flux through a surface as Φ_j = ∫ F ⋅ n_j dS, with F = B, F = ε D, or F =
  // E x H. The special coefficient is used to avoid issues evaluating MFEM GridFunctions
  // which are discontinuous at interior boundary elements.
  measurement_cache.surface_flux_i.clear();
  measurement_cache.surface_flux_i.reserve(surf_post_op.flux_surfs.size());
  for (const auto &[idx, data] : surf_post_op.flux_surfs)
  {
    measurement_cache.surface_flux_i.emplace_back(Measurement::FluxData{
        idx, surf_post_op.GetSurfaceFlux(idx, E.get(), B.get()), data.type});
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureFarField() const
{
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE)
  {
    measurement_cache.farfield.thetaphis = surf_post_op.farfield.thetaphis;

    // NOTE: measurement_cache.freq is omega (it has a factor of 2pi).
    measurement_cache.farfield.E_field = surf_post_op.GetFarFieldrE(
        measurement_cache.farfield.thetaphis, *E, *B, measurement_cache.freq);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureInterfaceEFieldEnergy() const
{
  // Depends on Lumped Port Energy since this is used in normalization of participation
  // ratio.

  // Compute the surface dielectric participation ratio and associated quality factor for
  // the material interface given by index idx. We have:
  //                            1/Q_mj = p_mj tan(δ)_j
  // with:
  //          p_mj = 1/2 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} / (E_elec + E_cap).
  measurement_cache.interface_eps_i.clear();
  measurement_cache.interface_edge_i.clear();
  measurement_cache.interface_local_edge_i.clear();
  if constexpr (HasEGridFunction<solver_t>())
  {
    // Domain and port energies must have been measured first. E_cap returns zero if the
    // solver does not support lumped ports.
    //
    // TODO: Should this not include other types of energy too (surface impedance case)?
    auto energy_electric_all = measurement_cache.domain_E_field_energy_all +
                               measurement_cache.lumped_port_capacitor_energy;

    measurement_cache.interface_eps_i.reserve(surf_post_op.eps_surfs.size());
    measurement_cache.interface_edge_i.reserve(surf_post_op.GetNInterfaceEdgeEntries());
    measurement_cache.interface_local_edge_i.reserve(
        surf_post_op.GetNInterfaceLocalEdgeEntries());
    surf_post_op.ResetInterfaceLocalEdgeEnergyCache();
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      auto energy =
          surf_post_op.GetInterfaceElectricFieldEnergy(idx, *E, D_recovered.get());

      auto energy_participation_p = energy / energy_electric_all;
      auto loss_tangent_delta = surf_post_op.GetInterfaceLossTangent(idx);
      auto quality_factor_Q = (energy_participation_p == 0.0 || loss_tangent_delta == 0.0)
                                  ? mfem::infinity()
                                  : 1.0 / (loss_tangent_delta * energy_participation_p);

      measurement_cache.interface_eps_i.emplace_back(Measurement::InterfaceData{
          idx, energy, loss_tangent_delta, energy_participation_p, quality_factor_Q});

      for (const auto &edge :
           surf_post_op.GetInterfaceEdgeElectricFieldEnergies(idx, *E, D_recovered.get()))
      {
        measurement_cache.interface_edge_i.emplace_back(Measurement::InterfaceEdgeData{
            idx, edge.distance, edge.energy_outside,
            edge.energy_outside / energy_electric_all, edge.energy_annulus,
            edge.energy_annulus / energy_electric_all});
      }
      if (!data.save_local_edge_energy)
      {
        continue;
      }
      for (const auto &edge : surf_post_op.GetInterfaceLocalEdgeElectricFieldEnergies(
               idx, *E, D_recovered.get()))
      {
        measurement_cache.interface_local_edge_i.emplace_back(
            Measurement::InterfaceLocalEdgeData{
                idx,
                edge.edge,
                edge.automatic,
                edge.component,
                edge.chain,
                edge.vertex_types,
                edge.process_normal,
                edge.p0,
                edge.p1,
                edge.length,
                edge.distance,
                edge.energy_total,
                edge.energy_total / energy_electric_all,
                edge.energy_total_polarized,
                {edge.energy_total_polarized[0] / energy_electric_all,
                 edge.energy_total_polarized[1] / energy_electric_all},
                edge.energy_inside,
                edge.energy_inside / energy_electric_all,
                edge.energy_inside_polarized,
                {edge.energy_inside_polarized[0] / energy_electric_all,
                 edge.energy_inside_polarized[1] / energy_electric_all},
                edge.energy_annulus,
                edge.energy_annulus / energy_electric_all,
                edge.energy_annulus_polarized,
                {edge.energy_annulus_polarized[0] / energy_electric_all,
                 edge.energy_annulus_polarized[1] / energy_electric_all},
                edge.energy_vertex_inside,
                edge.energy_vertex_inside / energy_electric_all,
                edge.energy_volume_annulus,
                edge.energy_volume_annulus / energy_electric_all,
                edge.energy_volume_vertex_annulus,
                edge.energy_volume_vertex_annulus / energy_electric_all,
                edge.energy_volume_annulus_polarized,
                {edge.energy_volume_annulus_polarized[0] / energy_electric_all,
                 edge.energy_volume_annulus_polarized[1] / energy_electric_all,
                 edge.energy_volume_annulus_polarized[2] / energy_electric_all,
                 edge.energy_volume_annulus_polarized[3] / energy_electric_all,
                 edge.energy_volume_annulus_polarized[4] / energy_electric_all,
                 edge.energy_volume_annulus_polarized[5] / energy_electric_all}});
      }
    }
  }
}

template <ProblemType solver_t>
std::map<int, std::string>
PostOperator<solver_t>::GetInterfaceTypes(const std::set<int> *targets) const
{
  std::map<int, std::string> types;
  for (const auto &[idx, data] : surf_post_op.eps_surfs)
  {
    if (!targets || targets->count(idx) > 0)
    {
      types[idx] = ToString(data.type);
    }
  }
  return types;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureSurfaceResponseCorrection() const
{
  surface_response_measurement.reset();
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::BOUNDARYMODE)
  {
    if (!surface_response_op || surface_response_op->IsRawFallback())
    {
      return;
    }
    BlockTimer response_timer(Timer::POSTPRO_RESPONSE);

    MaxwellSurfaceResponseMeasurement result;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // R0-a (idfix DESIGN B.3 (e)): the gate's exclusions so far, read fresh at every
    // measurement (the list grows with the sources; empty with 0 exclusions).
    std::vector<SurfacePostOperator::UncoveredPerimeterPortion>
        conductor_consistency_portions;
    for (const auto &portion : surface_response_op->GetConductorConsistencyPortions())
    {
      conductor_consistency_portions.push_back(
          {portion.p0, portion.p1, portion.topology, portion.feature});
    }
    result.raw_normalization_energy = measurement_cache.domain_E_field_energy_all +
                                      measurement_cache.lumped_port_capacitor_energy;
    MFEM_VERIFY(result.raw_normalization_energy > 0.0,
                "Maxwell surface-response normalization energy must be positive!");
    result.corrected_normalization_energy = nan;
    result.corrected_normalization_energy_fixed_flux = nan;
    result.self_consistent_normalization_energy = nan;
    result.self_consistent_frequency = {nan, nan};
    result.self_consistent_mode_overlap = nan;
    // Lane rescensus-product (decision 1817): this measurement's census sums (the
    // postprocessed response of the raw field = evaluation 0).
    ResolutionAccumulator resolution_sums = surface_response_op->NewResolutionAccumulator();
    if (surface_response_postprocess_enabled)
    {
      result.confidence = surface_response_op->GetMaxwellResponse(
          *E, measurement_cache.freq, surface_response_patch_energy_enabled,
          &resolution_sums);
      result.corrected_normalization_energy =
          result.raw_normalization_energy + result.confidence.domain_correction;
      result.corrected_normalization_energy_fixed_flux =
          result.raw_normalization_energy + result.confidence.domain_correction_fixed_flux;
      MFEM_VERIFY(result.corrected_normalization_energy > 0.0 &&
                      result.corrected_normalization_energy_fixed_flux > 0.0,
                  "Maxwell postprocessed surface-response normalization energy must be "
                  "positive!");
      result.has_postprocessed = true;
      if (result.confidence.has_magnetic)
      {
        // D6 (USER 1221; decision 1229): the corrected frequency of record.
        result.magnetic_domain_correction = result.confidence.magnetic_domain_correction;
        result.fixed_trace_frequency_shift =
            FixedTraceFrequencyShift(result.confidence, result.raw_normalization_energy);
      }
    }

    const auto target_interfaces = surface_response_op->GetTargetInterfaces();
    // Ruling 1738 (lane L1): the coverage cells' interfaces = the served targets and the
    // configured ones no loaded model serves (W = raw - outside R, nothing fabricated:
    // uncovered); every corrected value and every pre-existing table stays on the served.
    const auto coverage_interfaces = surface_response_op->GetCoverageInterfaces();
    const auto not_served_interfaces = surface_response_op->GetConfiguredTargetsNotServed();
    // The raw within-R energy of every uncovered requirement (decision 394 F2) and
    // DomainBoundary claim (F-DB-a) of the same field, kept in the corrected interface
    // energies as the electrostatic driver keeps it (decision 751): a feature the library
    // has no model for is neither modelled nor dropped.
    auto PortionEnergy =
        [](const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> &energies,
           int target)
    {
      const auto it = energies.find(target);
      return it != energies.end() ? it->second.energy : 0.0;
    };
    auto RemainderPerimeter =
        [](const SurfaceResponseOperator::LoopRemainderPortions &remainder_portions)
    {
      std::vector<SurfacePostOperator::UncoveredPerimeterPortion> perimeter;
      for (const auto &portion : remainder_portions.portions)
      {
        perimeter.push_back({portion.p0, portion.p1, portion.topology, portion.feature});
      }
      return perimeter;
    };
    // Lane L1 (decision 1634 (1); rulings 1641 / 1649 Q3): the closed coverage partition of
    // the RAW field (one pass per target interface), the served list's `unpatched` label
    // alone (no model claim under the point) moved to the UnpatchedPairs class (fix round
    // P2 2B); the loop remainder is moved by the caller.
    auto RawCoveragePartition = [&]()
    {
      const auto missing_by_feature =
          SurfacePostOperator::CoverageMissingByFeature(surface_response_uncovered);
      std::vector<SurfacePostOperator::CoveragePartitionClass> kept_classes = {
          {"Missing", &missing_by_feature},
          {"DomainBoundary", &surface_response_domain_boundary},
          {"SpatialClusterPriority", &surface_response_spatial_cluster_priority},
          {"ConductorConsistency", &conductor_consistency_portions}};
      std::vector<SurfacePostOperator::CoveragePartitionClass> remainder_classes;
      if (surface_response_exclusion_energy_enabled)
      {
        remainder_classes.push_back({"Excluded", &surface_response_exclusion_energy});
      }
      auto coverage = surf_post_op.GetInterfaceCoveragePartition(
          coverage_interfaces, *E, D_recovered.get(), kept_classes,
          surface_response_model_raw, remainder_classes);
      SurfacePostOperator::MoveCoverageServedLabel(coverage, "unpatched", "UnpatchedPairs");
      return coverage;
    };
    auto RootCoverageMeasurement =
        [&](const std::map<int, SurfacePostOperator::InterfaceCoveragePartition> &coverage,
            const std::map<int, double> &within, const std::map<int, double> &fabricated)
    {
      // Lane L1 (decision 1634 (2) / (4)): the coverage per Type of this measurement, on
      // the root (the writers' rank), before any writer runs.
      std::map<int, double> loop_remainder_excluded;
      for (const auto &[target, rows] : result.loop_remainder_excluded)
      {
        for (const auto &[label, energy] : rows)
        {
          (void)label;
          loop_remainder_excluded[target] += energy;
        }
      }
      auto measurement = ComputeCoverageMeasurement(
          coverage, GetInterfaceTypes(&coverage_interfaces), within, fabricated,
          surface_response_op->GetCoverageKeyRecords(),
          surface_response_op->GetCoverageModelQualifications(),
          surface_response_op->GetCoverageThresholds(),
          surface_response_op->GetCoverageValidatedRegime(),
          SurfaceResponseCoverageExcitation(), loop_remainder_excluded);
      MarkCoverageTargetsNotServed(measurement, GetInterfaceTypes(&not_served_interfaces));
      return measurement;
    };
    // Ruling 1694 (1): the untargeted postprocessing interfaces (their Types and raw
    // interface energies) in every measurement of the root.
    auto AddUntargeted = [&](CoverageMeasurement &measurement)
    {
      std::map<int, std::string> untargeted;
      for (const auto &[interface, type] : GetInterfaceTypes())
      {
        if (coverage_interfaces.count(interface) == 0)
        {
          untargeted[interface] = type;
        }
      }
      std::map<int, double> raw_energy;
      for (const auto &raw : measurement_cache.interface_eps_i)
      {
        raw_energy[raw.idx] = raw.energy;
      }
      AddCoverageUntargetedInterfaces(measurement, untargeted, raw_energy);
    };
    if (surface_response_postprocess_enabled)
    {
      BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
      if (!surface_response_uncovered.empty())
      {
        result.uncovered = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), surface_response_uncovered);
      }
      if (!surface_response_domain_boundary.empty())
      {
        result.domain_boundary = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), surface_response_domain_boundary);
      }
      if (!surface_response_spatial_cluster_priority.empty())
      {
        result.spatial_cluster_priority = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
            target_interfaces, *E, D_recovered.get(),
            surface_response_spatial_cluster_priority);
      }
      if (!conductor_consistency_portions.empty())
      {
        // R0-a (idfix DESIGN B.3 (e)): the conductor-consistency-excluded patches' raw
        // claims (grown by the gate before this measurement, S7diag R1), on the raw field.
        result.conductor_consistency = surf_post_op.GetInterfaceUncoveredEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), conductor_consistency_portions);
      }
      if (surface_response_patch_energy_enabled && !surface_response_model_raw.Empty())
      {
        // PatchEnergy (decision 922 B-1): the raw field's within-R energy per model set
        // (the claims exist on every run since lane L1; the table with PatchEnergy only).
        result.model_raw_partition = surf_post_op.GetInterfacePartitionedEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), surface_response_model_raw);
      }
      if (surface_response_exclusion_energy_enabled)
      {
        // Section 8 (decision 1076): the raw field's within-R energy per exclusion class
        // against the CLASSIFIED edge tree (+ "Claimed"); a diagnostic partition no energy
        // column of any other table reads.
        result.exclusion_energy_partition = surf_post_op.GetInterfaceClassifiedEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), surface_response_exclusion_energy);
        // The "Total" row: the within-R energy under the classified tree by ONE integral
        // (the P == S identity of the table).
        for (const auto &[interface, total] :
             surf_post_op.GetInterfaceClassifiedEdgeEnergies(target_interfaces, *E,
                                                             D_recovered.get(), {}))
        {
          for (const auto &[cls, energy] : total)
          {
            (void)cls;
            result.exclusion_energy_partition[interface]["Total"] += energy;
          }
        }
      }
      if (surface_response_uncovered_feature_energy_enabled &&
          !surface_response_uncovered.empty())
      {
        // Diagnostics.UncoveredFeatureEnergy (decision 1567 (b)): the uncovered energy of
        // the same field split per type by feature; a diagnostic partition no energy column
        // of any other table reads.
        result.uncovered_feature = surf_post_op.GetInterfaceUncoveredFeatureEdgeEnergies(
            target_interfaces, *E, D_recovered.get(), surface_response_uncovered);
      }
      if (surface_response_within_r_volume_energy_enabled)
      {
        // COL-3 (R-MAG-15): the raw field's within-box volume energies per patch (E and
        // B of this measurement; the material operator's eps / mu^-1).
        result.within_r_volume_energies =
            surface_response_op->GetWithinRVolumeEnergies(*E, *B, fem_op->GetMaterialOp());
      }
      {
        // Lane L1 (decision 1634 (1); rulings 1641 / 1649 Q3): the closed coverage
        // partition, every point in ONE bin by the correction's ownership (the kept-raw
        // classes first: the corrected total keeps them raw whatever claim lies under
        // them, e.g. a conductor-consistency-excluded patch's claim, REPORT STEP 10).
        result.coverage = RawCoveragePartition();
      }
      if (!result.confidence.loop_remainders.empty())
      {
        // OI-MX-LOOPGATE (decision 955): the circulating remainders' share of their
        // patches' raw within-R energy, kept raw like the uncovered class.
        const auto remainder_portions = surface_response_op->GetLoopRemainderPortions(
            result.confidence.loop_remainders);
        result.loop_remainder = surf_post_op.GetInterfaceWeightedPortionEdgeEnergies(
            target_interfaces, *E, D_recovered.get(),
            RemainderPerimeter(remainder_portions), remainder_portions.weights);
        // Lane L1 (ruling 1649 Q3): the loop remainder is kept raw by the corrected total:
        // out of the coverage partition's served bins.
        SurfacePostOperator::MoveCoverageLoopRemainder(result.coverage,
                                                       result.loop_remainder);
        // Ruling 1694 (i): the remainders of the conductor-consistency-EXCLUDED patches
        // (global patch indices of the gate's portion record), integrated alone.
        std::set<std::size_t> excluded_patches;
        for (const auto &patches :
             surface_response_op->GetConductorConsistencyPortionRecord().patches)
        {
          excluded_patches.insert(patches.begin(), patches.end());
        }
        std::vector<SurfaceResponseOperator::MaxwellResponse::LoopRemainder>
            excluded_remainders;
        for (const auto &remainder : result.confidence.loop_remainders)
        {
          if (remainder.patch >= 0 &&
              excluded_patches.count(static_cast<std::size_t>(remainder.patch)) > 0)
          {
            excluded_remainders.push_back(remainder);
          }
        }
        if (!excluded_remainders.empty())
        {
          const auto excluded_portions =
              surface_response_op->GetLoopRemainderPortions(excluded_remainders);
          result.loop_remainder_excluded =
              surf_post_op.GetInterfaceWeightedPortionEdgeEnergies(
                  target_interfaces, *E, D_recovered.get(),
                  RemainderPerimeter(excluded_portions), excluded_portions.weights);
        }
      }
    }
    auto RemainderEnergy =
        [](const std::map<int, std::map<std::string, double>> &energies, int target)
    {
      const auto it = energies.find(target);
      double total = 0.0;
      if (it != energies.end())
      {
        for (const auto &[label, energy] : it->second)
        {
          (void)label;
          total += energy;
        }
      }
      return total;
    };
    for (const auto &raw : measurement_cache.interface_eps_i)
    {
      auto &interface = result.interfaces[raw.idx];
      interface.raw_energy = raw.energy;
      interface.corrected_energy = nan;
      interface.corrected_energy_fixed_flux = nan;
      interface.self_consistent_energy = nan;
      interface.trace_closure_spread = nan;
      interface.loss_tangent = raw.tandelta;
      interface.response_target = target_interfaces.count(raw.idx) > 0;
      if (surface_response_postprocess_enabled && !interface.response_target)
      {
        // OI-ZERO-IFACE-NAN (decisions 891 / 922 C): an untargeted interface's corrected
        // energies are its raw energy, as the electrostatic driver carries them (no
        // response applies to it; nothing to close). The confidence gate never reads it:
        // the closure spread is formed over the response's fabricated (target) interfaces
        // and the matched / corner fractions are the identification's.
        interface.corrected_energy = raw.energy;
        interface.corrected_energy_fixed_flux = raw.energy;
        interface.trace_closure_spread = 0.0;
      }
    }
    if (surface_response_postprocess_enabled)
    {
      // Decision 1250: the matched energy share of the fixed-trace shift, accumulated over
      // the target interfaces from the same partition the electric correction uses.
      double shift_coverage_within_radius = 0.0, shift_coverage_kept_raw = 0.0;
      // Lane L1: W and the fabricated energy per target, the coverage cells' inputs.
      std::map<int, double> coverage_within, coverage_fabricated;
      for (const int target : target_interfaces)
      {
        auto interface = result.interfaces.find(target);
        MFEM_VERIFY(interface != result.interfaces.end(),
                    "Maxwell surface response refers to target interface "
                        << target << " which is not configured for postprocessing!");
        auto fabricated = result.confidence.fabricated_surface_energy.find(target);
        MFEM_VERIFY(fabricated != result.confidence.fabricated_surface_energy.end(),
                    "Maxwell surface response produced no fabricated energy for target "
                        << target << "!");
        auto fabricated_fixed_flux =
            result.confidence.fabricated_surface_energy_fixed_flux.find(target);
        MFEM_VERIFY(
            fabricated_fixed_flux !=
                result.confidence.fabricated_surface_energy_fixed_flux.end(),
            "Maxwell surface response produced no fixed-flux fabricated energy for target "
                << target << "!");

        const double radius = surface_response_op->GetMatchingRadius();
        auto outside = std::find_if(measurement_cache.interface_edge_i.begin(),
                                    measurement_cache.interface_edge_i.end(),
                                    [&](const auto &edge)
                                    {
                                      return edge.idx == target &&
                                             std::abs(edge.distance - radius) <=
                                                 1.0e-10 * std::max(edge.distance, radius);
                                    });
        MFEM_VERIFY(
            outside != measurement_cache.interface_edge_i.end(),
            "Maxwell response-corrected target interface "
                << target
                << " requires EdgeDistances containing the coupon matching radius!");
        // EdgeDistances contain the coupon matching radius: the largest configured radius
        // is the coupon's. The corrected energy keeps the raw energy outside R, replaces
        // the modelled within-R energy by the coupon response, and keeps the raw within-R
        // energy of the uncovered requirements and the DomainBoundary claims.
        const double kept_raw = PortionEnergy(result.uncovered, target) +
                                PortionEnergy(result.domain_boundary, target) +
                                PortionEnergy(result.spatial_cluster_priority, target) +
                                PortionEnergy(result.conductor_consistency, target) +
                                RemainderEnergy(result.loop_remainder, target);
        interface->second.corrected_energy =
            outside->energy_outside + fabricated->second + kept_raw;
        interface->second.corrected_energy_fixed_flux =
            outside->energy_outside + fabricated_fixed_flux->second + kept_raw;
        shift_coverage_within_radius +=
            interface->second.raw_energy - outside->energy_outside;
        shift_coverage_kept_raw += kept_raw;
        coverage_within[target] = interface->second.raw_energy - outside->energy_outside;
        coverage_fabricated[target] = fabricated->second;
        const double closure_scale =
            std::max(std::abs(fabricated->second), std::abs(fabricated_fixed_flux->second));
        if (closure_scale > 0.0)
        {
          interface->second.trace_closure_spread =
              std::abs(fabricated->second - fabricated_fixed_flux->second) / closure_scale;
        }
      }
      result.fixed_trace_shift_coverage =
          FixedTraceShiftCoverage(shift_coverage_within_radius, shift_coverage_kept_raw);
      // Ruling 1738: a configured target no loaded model serves: its W as a served target's
      // (raw - outside at R), nothing fabricated; its corrected values stay its raw ones.
      for (const int target : not_served_interfaces)
      {
        const auto interface = result.interfaces.find(target);
        const double radius = surface_response_op->GetMatchingRadius();
        const auto outside = std::find_if(
            measurement_cache.interface_edge_i.begin(),
            measurement_cache.interface_edge_i.end(),
            [&](const auto &edge)
            {
              return edge.idx == target && std::abs(edge.distance - radius) <=
                                               1.0e-10 * std::max(edge.distance, radius);
            });
        if (interface != result.interfaces.end() &&
            outside != measurement_cache.interface_edge_i.end())
        {
          coverage_within[target] = interface->second.raw_energy - outside->energy_outside;
          coverage_fabricated[target] = 0.0;
        }
      }
      if (Mpi::Root(fem_op->GetComm()))
      {
        result.coverage_measurement =
            RootCoverageMeasurement(result.coverage, coverage_within, coverage_fabricated);
        AddUntargeted(result.coverage_measurement);
      }
    }
    else
    {
      // Ruling 1689 M-6 (ruling 1654 (III) Q2 (b)): a SelfConsistent-only solve evaluates
      // no response on the raw field; until a bitwise proof of an extra raw evaluation and
      // its stated cost exist, its cells are NOT-MEASURED (no partition, no shares).
      if (Mpi::Root(fem_op->GetComm()))
      {
        result.coverage_measurement = NotMeasuredCoverageMeasurement(
            GetInterfaceTypes(&coverage_interfaces), SurfaceResponseCoverageExcitation());
        MarkCoverageTargetsNotServed(result.coverage_measurement,
                                     GetInterfaceTypes(&not_served_interfaces));
        AddUntargeted(result.coverage_measurement);
      }
    }
    // Lane rescensus-product (decision 1817): this measurement's census cells (collective;
    // the root's; a SelfConsistent-only solve has no evaluation-0 energy).
    result.resolution_cells = surface_response_op->MeasureResolution(
        surface_response_postprocess_enabled ? &resolution_sums : nullptr,
        units.Dimensionalize<Units::ValueType::ENERGY>(1.0));

    if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE)
    {
      if (surface_response_self_consistent_enabled)
      {
        MFEM_VERIFY(has_surface_response_corrected_field,
                    "Self-consistent surface response mode requires a corrected field!");
        auto &corrected_field = *surface_response_corrected_field;
        // The energies read the per-side corrected flux on a broken flux space (P3-14).
        const GridFunction *corrected_flux =
            !has_surface_response_corrected_flux ? nullptr
            : surface_response_corrected_flux_sides
                ? surface_response_corrected_flux_sides.get()
                : surface_response_corrected_flux.get();

        double capacitor_energy = 0.0;
        for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
        {
          (void)idx;
          if (std::abs(data.C) > 0.0)
          {
            const auto voltage = data.GetVoltage(corrected_field);
            capacitor_energy +=
                0.5 * std::abs(data.C) * std::real(voltage * std::conj(voltage));
          }
        }
        const auto corrected_frequency =
            surface_response_corrected_frequency.value_or(measurement_cache.freq);
        const auto corrected_response = surface_response_op->GetMaxwellResponse(
            corrected_field, corrected_frequency, surface_response_patch_energy_enabled);
        result.self_consistent_confidence = corrected_response;
        result.self_consistent_frequency = corrected_frequency;
        result.self_consistent_mode_overlap =
            surface_response_corrected_mode_overlap.value_or(1.0);
        result.self_consistent_normalization_energy =
            dom_post_op.GetElectricFieldEnergy(corrected_field) + capacitor_energy +
            corrected_response.domain_correction;
        MFEM_VERIFY(result.self_consistent_normalization_energy > 0.0,
                    "Self-consistent Maxwell surface-response normalization energy must "
                    "be positive!");

        for (auto &[index, interface] : result.interfaces)
        {
          interface.self_consistent_energy = surf_post_op.GetInterfaceElectricFieldEnergy(
              index, corrected_field, corrected_flux);
        }
        {
          BlockTimer energy_timer(Timer::POSTPRO_RESPONSE_ENERGY);
          if (!surface_response_uncovered.empty())
          {
            result.self_consistent_uncovered =
                surf_post_op.GetInterfaceUncoveredEdgeEnergies(
                    target_interfaces, corrected_field, corrected_flux,
                    surface_response_uncovered);
          }
          if (!surface_response_domain_boundary.empty())
          {
            result.self_consistent_domain_boundary =
                surf_post_op.GetInterfaceUncoveredEdgeEnergies(
                    target_interfaces, corrected_field, corrected_flux,
                    surface_response_domain_boundary);
          }
          if (!surface_response_spatial_cluster_priority.empty())
          {
            result.self_consistent_spatial_cluster_priority =
                surf_post_op.GetInterfaceUncoveredEdgeEnergies(
                    target_interfaces, corrected_field, corrected_flux,
                    surface_response_spatial_cluster_priority);
          }
          if (!conductor_consistency_portions.empty())
          {
            result.self_consistent_conductor_consistency =
                surf_post_op.GetInterfaceUncoveredEdgeEnergies(
                    target_interfaces, corrected_field, corrected_flux,
                    conductor_consistency_portions);
          }
          if (!corrected_response.loop_remainders.empty())
          {
            const auto remainder_portions = surface_response_op->GetLoopRemainderPortions(
                corrected_response.loop_remainders);
            result.self_consistent_loop_remainder =
                surf_post_op.GetInterfaceWeightedPortionEdgeEnergies(
                    target_interfaces, corrected_field, corrected_flux,
                    RemainderPerimeter(remainder_portions), remainder_portions.weights);
          }
        }
        for (const int target : target_interfaces)
        {
          auto fabricated = corrected_response.fabricated_surface_energy.find(target);
          MFEM_VERIFY(fabricated != corrected_response.fabricated_surface_energy.end(),
                      "Self-consistent Maxwell surface response produced no fabricated "
                      "energy for target "
                          << target << "!");
          const auto edge_energies = surf_post_op.GetInterfaceEdgeElectricFieldEnergies(
              target, corrected_field, corrected_flux);
          const double radius = surface_response_op->GetMatchingRadius();
          auto outside = std::find_if(edge_energies.begin(), edge_energies.end(),
                                      [&](const auto &edge)
                                      {
                                        return std::abs(edge.distance - radius) <=
                                               1.0e-10 * std::max(edge.distance, radius);
                                      });
          MFEM_VERIFY(outside != edge_energies.end(),
                      "Self-consistent Maxwell response-corrected target interface "
                          << target
                          << " requires EdgeDistances containing the coupon matching "
                             "radius!");
          result.interfaces.at(target).self_consistent_energy =
              outside->energy_outside + fabricated->second +
              PortionEnergy(result.self_consistent_uncovered, target) +
              PortionEnergy(result.self_consistent_domain_boundary, target) +
              PortionEnergy(result.self_consistent_spatial_cluster_priority, target) +
              PortionEnergy(result.self_consistent_conductor_consistency, target) +
              RemainderEnergy(result.self_consistent_loop_remainder, target);
        }
        result.has_self_consistent = true;
      }
    }

    bool self_consistent_confident =
        result.has_self_consistent &&
        result.self_consistent_confidence.closure_independent_confident;
    if constexpr (solver_t == ProblemType::EIGENMODE)
    {
      self_consistent_confident =
          self_consistent_confident &&
          result.self_consistent_mode_overlap >= minimum_eigenmode_response_overlap;
    }
    if (result.has_postprocessed && !result.confidence.confident &&
        !surface_response_warning_printed)
    {
      const std::string self_consistent_status =
          !result.has_self_consistent
              ? "No self-consistent corrected field is available for this solve."
          : self_consistent_confident
              ? "The self-consistent corrected result passes the non-closure "
                "confidence gates and is preferred."
              : "The self-consistent corrected result also fails a non-closure "
                "confidence gate or, for an eigenmode, the minimum mode-overlap gate.";
      Mpi::Warning(
          "Maxwell postprocessing-only surface-response confidence limits were exceeded: "
          "kR = {:.3e}, loop residual = {:.3e}, response-weighted loop residual = "
          "{:.3e}, loop-response fraction above limit = {:.3e}, matched fraction = "
          "{:.6f}, "
          "unmodeled corner fraction = {:.3e}, max R/rho = {:.3e}, "
          "library distance = {:.3e}, boundary-law parameters verified = {:d}, "
          "max interface trace-closure spread = {:.3e}, response-weighted local "
          "trace-closure spread = {:.3e}, trace-closure response fraction above limit = "
          "{:.3e}. "
          "Fixed-trace and fixed-flux values are reported but should not be treated as "
          "validated. {}\n",
          result.confidence.kR, result.confidence.loop_residual,
          result.confidence.response_weighted_loop_residual,
          result.confidence.loop_response_failure_fraction,
          result.confidence.matched_length_fraction,
          result.confidence.corner_neighborhood_fraction,
          result.confidence.maximum_curvature_ratio,
          result.confidence.maximum_library_distance,
          result.confidence.boundary_law_verified ? 1 : 0,
          result.confidence.maximum_trace_closure_spread,
          result.confidence.response_weighted_trace_closure_spread,
          result.confidence.trace_closure_response_failure_fraction,
          self_consistent_status);
      surface_response_warning_printed = true;
    }
    if (!result.has_postprocessed && result.has_self_consistent &&
        !self_consistent_confident && !surface_response_warning_printed)
    {
      Mpi::Warning(
          "Self-consistent Maxwell surface-response confidence limits were exceeded: "
          "kR = {:.3e}, response-weighted loop residual = {:.3e}, loop-response "
          "fraction above limit = {:.3e}, matched fraction = {:.6f}, unmodeled corner "
          "fraction = {:.3e}, max R/rho = {:.3e}, library distance = {:.3e}, "
          "boundary-law parameters verified = {:d}, mode overlap = {:.3e}. Corrected "
          "values are reported but should not be treated as validated.\n",
          result.self_consistent_confidence.kR,
          result.self_consistent_confidence.response_weighted_loop_residual,
          result.self_consistent_confidence.loop_response_failure_fraction,
          result.self_consistent_confidence.matched_length_fraction,
          result.self_consistent_confidence.corner_neighborhood_fraction,
          result.self_consistent_confidence.maximum_curvature_ratio,
          result.self_consistent_confidence.maximum_library_distance,
          result.self_consistent_confidence.boundary_law_verified ? 1 : 0,
          result.self_consistent_mode_overlap);
      surface_response_warning_printed = true;
    }
    surface_response_measurement = std::move(result);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseCorrection(double output_index,
                                                            int excitation) const
{
  if (!surface_response_measurement || !Mpi::Root(fem_op->GetComm()))
  {
    return;
  }
  // Lane L4 (decision 1829; DESIGN-L4-v2 s2.3 M1-M4): under Coverage.Policy
  // "MarkNotApplied" the corrected columns below print nan where their cell is marked and
  // the value moves to the appended partial column (the column exists only under that
  // policy); under "Warn" every cell and column is printed as before.
  const bool marking = surface_response_op->CoverageMarkNotApplied();
  // M1 (per interface: its Type's cell) and M2 (per measurement: any marked cell).
  static constexpr std::array<const char *, 3> marked_measurement_columns = {
      "domain_corrected", "domain_corrected_fixed_flux", "domain_self_consistent"};
  static constexpr std::array<const char *, 9> marked_interface_columns = {
      "energy_corrected",
      "participation_corrected",
      "quality_corrected",
      "energy_corrected_fixed_flux",
      "participation_corrected_fixed_flux",
      "quality_corrected_fixed_flux",
      "energy_self_consistent",
      "participation_self_consistent",
      "quality_self_consistent"};
  // M3 (surface-response-confidence.csv; per measurement), the columns this run has.
  auto MarkedConfidenceColumns = [this]()
  {
    std::vector<std::string> names;
    if constexpr (solver_t == ProblemType::EIGENMODE)
    {
      names = {"self_consistent_frequency_re", "self_consistent_frequency_im"};
    }
    if (surface_response_op->HasMagneticLeaves())
    {
      for (const char *name :
           {"fixed_trace_frequency_shift", "fixed_trace_corrected_frequency",
            "magnetic_domain_correction", "self_consistent_frequency_magnetic_re"})
      {
        names.push_back(name);
      }
    }
    return names;
  };
  auto AppendPartialColumns = [](Table &table, const auto &names)
  {
    for (const auto &name : names)
    {
      const std::string header = table[name].header_text;
      table.insert(fmt::format("{}_partial", name), header + kCoveragePartialColumnSuffix);
    }
  };
  // A marked cell prints nan, its value in the partial column (nan when not marked).
  auto PrintMarked =
      [marking](Table &table, const std::string &name, double value, bool marked)
  {
    if (!marking)
    {
      table[name] << value;
      return;
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    table[name] << (marked ? nan : value);
    table[name + "_partial"] << (marked ? value : nan);
  };
  if (!surface_Q_corrected)
  {
    surface_Q_corrected =
        std::make_unique<TableWithCSVFile>(post_dir / "surface-Q-corrected.csv");
    auto &table = surface_Q_corrected->table;
    table.insert("idx", solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m");
    table.insert("excitation", "exc");
    table.insert("interface", "interface");
    table.insert("domain_raw", "E_norm raw (J)");
    table.insert("domain_corrected", "E_norm corrected (J)");
    table.insert("domain_corrected_fixed_flux", "E_norm corrected fixed-flux (J)");
    table.insert("domain_self_consistent", "E_norm self-consistent (J)");
    table.insert("energy_raw", "E_surf raw (J)");
    table.insert("participation_raw", "p_surf raw");
    table.insert("quality_raw", "Q_surf raw");
    table.insert("energy_corrected", "E_surf corrected (J)");
    table.insert("participation_corrected", "p_surf corrected");
    table.insert("quality_corrected", "Q_surf corrected");
    table.insert("energy_corrected_fixed_flux", "E_surf corrected fixed-flux (J)");
    table.insert("participation_corrected_fixed_flux", "p_surf corrected fixed-flux");
    table.insert("quality_corrected_fixed_flux", "Q_surf corrected fixed-flux");
    table.insert("energy_self_consistent", "E_surf self-consistent (J)");
    table.insert("participation_self_consistent", "p_surf self-consistent");
    table.insert("quality_self_consistent", "Q_surf self-consistent");
    table.insert("trace_closure_spread", "trace closure spread");
    table.insert("matched_fraction", "matched edge length fraction");
    table.insert("corner_fraction", "unmodeled corner neighborhood fraction");
    // OI-ZERO-IFACE-NAN (decision 922 C-3 (ii)): the "response target" flag column only
    // when a postprocessing interface is not a target (a run whose interfaces are all
    // targets keeps its column set byte for byte).
    const auto &measured = *surface_response_measurement;
    surface_response_untargeted_interface =
        std::any_of(measured.interfaces.begin(), measured.interfaces.end(),
                    [](const auto &entry) { return !entry.second.response_target; });
    if (surface_response_untargeted_interface)
    {
      table.insert("response_target", "response target");
      table["response_target"].print_as_int = true;
    }
    // Lane L1 (DESIGN-v2 s4.3; decision 1636 (1)): the coverage columns, APPENDED (every
    // column above keeps its position, name and value): this interface's Type, W and
    // uncovered share, then its Type's cell (repeated on each of its interfaces).
    table.insert(Column::Text("coverage_type", "coverage Type"));
    table.insert("coverage_within", "W within-R raw (J)");
    table.insert("coverage_uncovered_share", "uncovered raw share of W");
    table.insert("coverage_type_uncovered_share", "Type uncovered raw share of W");
    table.insert("coverage_type_method_gap_share", "Type method-gap raw share of W");
    table.insert("coverage_type_pending_share", "Type pending-model share of W");
    table.insert("coverage_type_factor", "Type covered correction factor");
    table.insert(Column::Text("coverage_status", "coverage status"));
    table.insert(Column::Text("regime_status", "regime status"));
    // Rulings 1689 M-2 / 1694 (i): this interface's loop-remainder terms (J).
    table.insert("coverage_loop_remainder_unserved", "loop remainder unserved (J)");
    table.insert("coverage_loop_remainder_excluded",
                 "loop remainder of excluded patches kept raw twice (J)");
    if (marking)
    {
      AppendPartialColumns(table, marked_measurement_columns);
      AppendPartialColumns(table, marked_interface_columns);
      // Lane L4 round 2 (decision 1914): the explicit marks, last.
      table.insert("coverage_mark_normalization", kCoverageMarkNormalizationColumn);
      table.insert("coverage_mark_interface", kCoverageMarkInterfaceColumn);
      table["coverage_mark_normalization"].print_as_int = true;
      table["coverage_mark_interface"].print_as_int = true;
    }
    table["excitation"].print_as_int = true;
    table["interface"].print_as_int = true;
    if constexpr (solver_t == ProblemType::EIGENMODE ||
                  solver_t == ProblemType::BOUNDARYMODE)
    {
      table["idx"].print_as_int = true;
    }
  }
  if (!surface_response_confidence)
  {
    surface_response_confidence =
        std::make_unique<TableWithCSVFile>(post_dir / "surface-response-confidence.csv");
    auto &table = surface_response_confidence->table;
    table.insert("idx", solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m");
    table.insert("excitation", "exc");
    table.insert("kR", "kR");
    table.insert("loop_residual", "loop residual");
    table.insert("weighted_loop_residual", "response-weighted loop residual");
    table.insert("loop_failure_fraction", "loop-response fraction above limit");
    table.insert("matched_fraction", "matched edge length fraction");
    table.insert("corner_fraction", "unmodeled corner neighborhood fraction");
    table.insert("curvature", "max R/rho");
    table.insert("library_distance", "max library distance");
    table.insert("boundary_law_verified", "boundary-law parameters verified");
    table.insert("trace_closure_spread", "max trace closure spread");
    table.insert("weighted_trace_closure_spread",
                 "response-weighted local trace closure spread");
    table.insert("trace_closure_failure_fraction",
                 "trace-closure response fraction above limit");
    table.insert("confident", "confidence pass");
    table.insert("self_consistent_weighted_loop_residual",
                 "self-consistent response-weighted loop residual");
    table.insert("self_consistent_loop_failure_fraction",
                 "self-consistent loop-response fraction above limit");
    table.insert("self_consistent_confident", "self-consistent confidence pass");
    table["excitation"].print_as_int = true;
    table["confident"].print_as_int = true;
    table["boundary_law_verified"].print_as_int = true;
    if constexpr (solver_t == ProblemType::EIGENMODE ||
                  solver_t == ProblemType::BOUNDARYMODE)
    {
      table["idx"].print_as_int = true;
    }
    if constexpr (solver_t == ProblemType::EIGENMODE)
    {
      table.insert("self_consistent_frequency_re", "Re{f self-consistent} (GHz)");
      table.insert("self_consistent_frequency_im", "Im{f self-consistent} (GHz)");
      table.insert("self_consistent_mode_overlap", "self-consistent M-overlap");
    }
    // Cap-interior hats (decision 818 R2), the last columns: information, not a gate
    // (decision 814 (1)).
    table.insert("hat_residual", "hat residual");
    table.insert("weighted_hat_residual", "response-weighted hat residual");
    if (surface_response_op->HasMagneticLeaves())
    {
      // The magnetic (inductive) term's columns (DESIGN-v2; decisions 1229 Q1 / Q3, 1235):
      // only when a magnetic leaf was READ (also when every leaf is refused: NaN numerics
      // then), so every other run keeps its column set byte for byte. `fixed-trace
      // frequency shift` is the corrected frequency of record (USER 1221 D6); `magnetic
      // term applied` says whether a non-zero dK was in the corrected operator the
      // self-consistent columns report; the fallback columns are the magnetic image of
      // the loop-fit CSV's withheld correction (decision 955 D6); `fixed-trace shift
      // coverage` is the matched energy share the shift rests on (decision 1250: the shift
      // follows the electric correction's coverage rule, the confidence flags apart).
      table.insert("fixed_trace_frequency_shift", "fixed-trace frequency shift");
      table.insert("fixed_trace_corrected_frequency", "f fixed-trace corrected (GHz)");
      table.insert("fixed_trace_shift_coverage", "fixed-trace shift coverage");
      table.insert("magnetic_domain_correction", "magnetic domain correction (J)");
      table.insert("magnetic_loop_residual", "magnetic loop residual");
      table.insert("magnetic_weighted_loop_residual",
                   "magnetic response-weighted loop residual");
      table.insert("magnetic_loop_failure_fraction",
                   "magnetic loop-response fraction above limit");
      table.insert("magnetic_confident", "magnetic confidence pass");
      table.insert("magnetic_patches", "magnetic patches");
      table.insert("magnetic_excluded_patches", "magnetic patches excluded");
      table.insert("magnetic_fallback_patches", "magnetic fallback patches");
      table.insert("withheld_magnetic_domain_correction",
                   "withheld magnetic domain correction (J)");
      table.insert("magnetic_term_applied", "magnetic term applied");
      table.insert("self_consistent_frequency_magnetic_re",
                   "Re{f self-consistent (magnetic)} (GHz)");
      table["magnetic_confident"].print_as_int = true;
      table["magnetic_patches"].print_as_int = true;
      table["magnetic_excluded_patches"].print_as_int = true;
      table["magnetic_fallback_patches"].print_as_int = true;
      table["magnetic_term_applied"].print_as_int = true;
    }
    // Lane L1 (DESIGN-v2 s4.3): the worst cell over Types, APPENDED last.
    table.insert("coverage_worst_uncovered_share", "worst uncovered raw share of W");
    table.insert(Column::Text("coverage_worst_status", "worst coverage status"));
    table.insert(Column::Text("coverage_worst_regime", "worst regime status"));
    if (marking)
    {
      // Lane L4 (M3 / M4): the partial columns, then the raw within-R energy the marked
      // cells leave uncorrected over E_norm raw (the exposure of the corrected
      // normalization: every unmarked Type's p / Q in this measurement is built on it).
      AppendPartialColumns(table, MarkedConfidenceColumns());
      table.insert("coverage_marked_uncovered_share",
                   "not-applied uncovered raw share of E_norm");
      // Lane L4 round 2 (decision 1914): the explicit mark of the row's measurement cells.
      table.insert("coverage_mark_measurement", kCoverageMarkMeasurementColumn);
      table["coverage_mark_measurement"].print_as_int = true;
    }
  }

  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  const auto &coverage = result.coverage_measurement;
  const auto coverage_types = GetInterfaceTypes();
  const auto &confidence = result.confidence;
  // M2 / M3: the measurement's whole values are marked when any cell is (fail-closed: a
  // measurement without cells is NOT-MEASURED).
  const bool measurement_marked = marking && coverage.MarksWholeMeasurement();
  surface_response_measurement_marked = measurement_marked;
  for (const auto &[interface_index, interface] : result.interfaces)
  {
    const double p_raw = interface.raw_energy / result.raw_normalization_energy;
    const double p_corrected =
        interface.corrected_energy / result.corrected_normalization_energy;
    const double p_corrected_fixed_flux = interface.corrected_energy_fixed_flux /
                                          result.corrected_normalization_energy_fixed_flux;
    const double q_raw = p_raw == 0.0 || interface.loss_tangent == 0.0
                             ? mfem::infinity()
                             : 1.0 / (p_raw * interface.loss_tangent);
    const double q_corrected = p_corrected == 0.0 || interface.loss_tangent == 0.0
                                   ? mfem::infinity()
                                   : 1.0 / (p_corrected * interface.loss_tangent);
    const double q_corrected_fixed_flux =
        p_corrected_fixed_flux == 0.0 || interface.loss_tangent == 0.0
            ? mfem::infinity()
            : 1.0 / (p_corrected_fixed_flux * interface.loss_tangent);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double p_self_consistent =
        result.has_self_consistent
            ? interface.self_consistent_energy / result.self_consistent_normalization_energy
            : nan;
    const double q_self_consistent =
        !result.has_self_consistent
            ? nan
            : (p_self_consistent == 0.0 || interface.loss_tangent == 0.0
                   ? mfem::infinity()
                   : 1.0 / (p_self_consistent * interface.loss_tangent));
    const auto matched_fraction =
        confidence.matched_length_fraction_by_interface.find(interface_index);
    const auto corner_fraction =
        confidence.corner_neighborhood_fraction_by_interface.find(interface_index);
    auto &table = surface_Q_corrected->table;
    // M1: a coverage interface (a target, or a configured target no model serves) of a
    // marked Type (its Type's cell; a target without a cell reads NOT-MEASURED: marked).
    const auto marked_type = coverage_types.find(interface_index);
    const CoverageCell *marked_cell = marked_type != coverage_types.end()
                                          ? coverage.FindCell(marked_type->second)
                                          : nullptr;
    const bool interface_marked =
        marking &&
        (interface.response_target || coverage.not_served.count(interface_index) > 0) &&
        (!marked_cell || coverage.IsTypeMarked(marked_type->second));
    table["idx"] << output_index;
    table["excitation"] << excitation;
    table["interface"] << interface_index;
    table["domain_raw"] << units.Dimensionalize<VT::ENERGY>(
        result.raw_normalization_energy);
    PrintMarked(table, "domain_corrected",
                units.Dimensionalize<VT::ENERGY>(result.corrected_normalization_energy),
                measurement_marked);
    PrintMarked(
        table, "domain_corrected_fixed_flux",
        units.Dimensionalize<VT::ENERGY>(result.corrected_normalization_energy_fixed_flux),
        measurement_marked);
    PrintMarked(
        table, "domain_self_consistent",
        result.has_self_consistent
            ? units.Dimensionalize<VT::ENERGY>(result.self_consistent_normalization_energy)
            : nan,
        measurement_marked);
    table["energy_raw"] << units.Dimensionalize<VT::ENERGY>(interface.raw_energy);
    table["participation_raw"] << p_raw;
    table["quality_raw"] << q_raw;
    PrintMarked(table, "energy_corrected",
                units.Dimensionalize<VT::ENERGY>(interface.corrected_energy),
                interface_marked);
    PrintMarked(table, "participation_corrected", p_corrected, interface_marked);
    PrintMarked(table, "quality_corrected", q_corrected, interface_marked);
    PrintMarked(table, "energy_corrected_fixed_flux",
                units.Dimensionalize<VT::ENERGY>(interface.corrected_energy_fixed_flux),
                interface_marked);
    PrintMarked(table, "participation_corrected_fixed_flux", p_corrected_fixed_flux,
                interface_marked);
    PrintMarked(table, "quality_corrected_fixed_flux", q_corrected_fixed_flux,
                interface_marked);
    PrintMarked(table, "energy_self_consistent",
                result.has_self_consistent
                    ? units.Dimensionalize<VT::ENERGY>(interface.self_consistent_energy)
                    : nan,
                interface_marked);
    PrintMarked(table, "participation_self_consistent", p_self_consistent,
                interface_marked);
    PrintMarked(table, "quality_self_consistent", q_self_consistent, interface_marked);
    table["trace_closure_spread"] << interface.trace_closure_spread;
    // An untargeted interface (922 C-2): nothing is matched or modelled there, 0 / 0.
    const bool untargeted = result.has_postprocessed && !interface.response_target;
    table["matched_fraction"]
        << (matched_fraction != confidence.matched_length_fraction_by_interface.end()
                ? matched_fraction->second
                : (untargeted ? 0.0 : nan));
    table["corner_fraction"]
        << (corner_fraction != confidence.corner_neighborhood_fraction_by_interface.end()
                ? corner_fraction->second
                : (untargeted ? 0.0 : nan));
    if (surface_response_untargeted_interface)
    {
      table["response_target"] << (interface.response_target ? 1.0 : 0.0);
    }
    {
      // Lane L1: an untargeted interface is in no Type (NOT-TARGET); a target of an
      // unmeasured cell reads its cell's NOT-MEASURED row.
      const auto type = coverage_types.find(interface_index);
      table["coverage_type"] << (type != coverage_types.end() ? type->second : "Default");
      const auto summary = coverage.interfaces.find(interface_index);
      const CoverageCell *cell =
          type != coverage_types.end() ? coverage.FindCell(type->second) : nullptr;
      // Ruling 1738: a configured target no loaded model serves is a coverage interface
      // (its W, its uncovered share, its Type's TARGET-NOT-SERVED cell), not NOT-TARGET.
      if (!interface.response_target && coverage.not_served.count(interface_index) == 0)
      {
        // Ruling 1694 (A): no target = no correction radius = an empty within-R region
        // (W = 0); the shares and the factor nan; its raw energy stays in E_surf raw.
        table["coverage_within"] << 0.0;
        table["coverage_uncovered_share"] << nan;
        table["coverage_type_uncovered_share"] << nan;
        table["coverage_type_method_gap_share"] << nan;
        table["coverage_type_pending_share"] << nan;
        table["coverage_type_factor"] << nan;
        table["coverage_status"]
            << std::string(CoverageStatusName(CoverageStatus::NOT_TARGET));
        table["regime_status"] << std::string(RegimeStatusName(RegimeStatus::NO_MODEL));
        table["coverage_loop_remainder_unserved"] << 0.0;
        table["coverage_loop_remainder_excluded"] << 0.0;
      }
      else if (!cell)
      {
        table["coverage_within"] << nan;
        table["coverage_uncovered_share"] << nan;
        table["coverage_type_uncovered_share"] << nan;
        table["coverage_type_method_gap_share"] << nan;
        table["coverage_type_pending_share"] << nan;
        table["coverage_type_factor"] << nan;
        table["coverage_status"]
            << std::string(CoverageStatusName(CoverageStatus::NOT_MEASURED));
        table["regime_status"] << std::string(RegimeStatusName(RegimeStatus::NO_MODEL));
        table["coverage_loop_remainder_unserved"] << nan;
        table["coverage_loop_remainder_excluded"] << nan;
      }
      else
      {
        const bool measured = summary != coverage.interfaces.end();
        const double within = measured ? summary->second.within : nan;
        table["coverage_within"] << units.Dimensionalize<VT::ENERGY>(within);
        table["coverage_uncovered_share"]
            << (measured && within > 0.0 ? (within - summary->second.served) / within
                                         : nan);
        table["coverage_type_uncovered_share"] << cell->uncovered_share;
        table["coverage_type_method_gap_share"] << cell->method_gap_share;
        table["coverage_type_pending_share"] << cell->pending_share;
        table["coverage_type_factor"] << cell->factor;
        table["coverage_status"] << std::string(CoverageStatusName(cell->status));
        table["regime_status"] << std::string(RegimeStatusName(cell->regime));
        table["coverage_loop_remainder_unserved"]
            << (measured ? units.Dimensionalize<VT::ENERGY>(
                               summary->second.loop_remainder_unserved)
                         : nan);
        table["coverage_loop_remainder_excluded"]
            << (measured ? units.Dimensionalize<VT::ENERGY>(
                               summary->second.loop_remainder_excluded)
                         : nan);
      }
    }
    if (marking)
    {
      table["coverage_mark_normalization"] << (measurement_marked ? 1.0 : 0.0);
      table["coverage_mark_interface"] << (interface_marked ? 1.0 : 0.0);
    }
  }
  surface_Q_corrected->WriteFullTableTrunc();

  auto &table = surface_response_confidence->table;
  table["idx"] << output_index;
  table["excitation"] << excitation;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  table["kR"] << (result.has_postprocessed ? confidence.kR : nan);
  table["loop_residual"] << (result.has_postprocessed ? confidence.loop_residual : nan);
  table["weighted_loop_residual"]
      << (result.has_postprocessed ? confidence.response_weighted_loop_residual : nan);
  table["loop_failure_fraction"]
      << (result.has_postprocessed ? confidence.loop_response_failure_fraction : nan);
  table["matched_fraction"] << (result.has_postprocessed
                                    ? confidence.matched_length_fraction
                                    : nan);
  table["corner_fraction"] << (result.has_postprocessed
                                   ? confidence.corner_neighborhood_fraction
                                   : nan);
  table["curvature"] << (result.has_postprocessed ? confidence.maximum_curvature_ratio
                                                  : nan);
  table["library_distance"] << (result.has_postprocessed
                                    ? confidence.maximum_library_distance
                                    : nan);
  table["boundary_law_verified"]
      << (result.has_postprocessed ? (confidence.boundary_law_verified ? 1.0 : 0.0) : nan);
  table["trace_closure_spread"]
      << (result.has_postprocessed ? confidence.maximum_trace_closure_spread : nan);
  table["weighted_trace_closure_spread"]
      << (result.has_postprocessed ? confidence.response_weighted_trace_closure_spread
                                   : nan);
  table["trace_closure_failure_fraction"]
      << (result.has_postprocessed ? confidence.trace_closure_response_failure_fraction
                                   : nan);
  table["confident"] << (result.has_postprocessed ? (confidence.confident ? 1.0 : 0.0)
                                                  : nan);
  table["self_consistent_weighted_loop_residual"]
      << (result.has_self_consistent
              ? result.self_consistent_confidence.response_weighted_loop_residual
              : nan);
  table["self_consistent_loop_failure_fraction"]
      << (result.has_self_consistent
              ? result.self_consistent_confidence.loop_response_failure_fraction
              : nan);
  bool self_consistent_confident =
      result.has_self_consistent &&
      result.self_consistent_confidence.closure_independent_confident;
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    self_consistent_confident =
        self_consistent_confident &&
        result.self_consistent_mode_overlap >= minimum_eigenmode_response_overlap;
    PrintMarked(table, "self_consistent_frequency_re",
                result.has_self_consistent ? units.Dimensionalize<VT::FREQUENCY>(
                                                 result.self_consistent_frequency.real()) /
                                                 (2 * M_PI)
                                           : nan,
                measurement_marked);
    PrintMarked(table, "self_consistent_frequency_im",
                result.has_self_consistent ? units.Dimensionalize<VT::FREQUENCY>(
                                                 result.self_consistent_frequency.imag()) /
                                                 (2 * M_PI)
                                           : nan,
                measurement_marked);
    table["self_consistent_mode_overlap"]
        << (result.has_self_consistent ? result.self_consistent_mode_overlap : nan);
  }
  table["self_consistent_confident"]
      << (result.has_self_consistent ? (self_consistent_confident ? 1.0 : 0.0) : nan);
  table["hat_residual"] << (result.has_postprocessed ? confidence.hat_residual : nan);
  table["weighted_hat_residual"]
      << (result.has_postprocessed ? confidence.response_weighted_hat_residual : nan);
  if (surface_response_op->HasMagneticLeaves())
  {
    const bool magnetic = result.has_postprocessed && confidence.has_magnetic;
    const double shift = result.fixed_trace_frequency_shift;
    // The raw frequency of this row in GHz (the driven sweep's frequency, the mode's).
    const double frequency_GHz =
        units.Dimensionalize<VT::FREQUENCY>(measurement_cache.freq.real()) / (2 * M_PI);
    PrintMarked(table, "fixed_trace_frequency_shift", magnetic ? shift : nan,
                measurement_marked);
    PrintMarked(table, "fixed_trace_corrected_frequency",
                magnetic && std::isfinite(shift) ? frequency_GHz * (1.0 + shift) : nan,
                measurement_marked);
    // Decision 1250: the matched energy share the shift rests on (1 at full coverage).
    table["fixed_trace_shift_coverage"]
        << (magnetic ? result.fixed_trace_shift_coverage : nan);
    PrintMarked(
        table, "magnetic_domain_correction",
        magnetic ? units.Dimensionalize<VT::ENERGY>(confidence.magnetic_domain_correction)
                 : nan,
        measurement_marked);
    table["magnetic_loop_residual"] << (magnetic ? confidence.magnetic_loop_residual : nan);
    table["magnetic_weighted_loop_residual"]
        << (magnetic ? confidence.magnetic_weighted_loop_residual : nan);
    table["magnetic_loop_failure_fraction"]
        << (magnetic ? confidence.magnetic_loop_response_failure_fraction : nan);
    table["magnetic_confident"]
        << (magnetic ? (confidence.magnetic_confident ? 1.0 : 0.0) : nan);
    table["magnetic_patches"] << (magnetic
                                      ? static_cast<double>(confidence.magnetic_patch_count)
                                      : nan);
    table["magnetic_excluded_patches"]
        << (magnetic ? static_cast<double>(confidence.magnetic_excluded_patch_count) : nan);
    table["magnetic_fallback_patches"]
        << (magnetic ? static_cast<double>(confidence.magnetic_fallback_patch_count) : nan);
    table["withheld_magnetic_domain_correction"]
        << (magnetic ? units.Dimensionalize<VT::ENERGY>(
                           confidence.loop_fallback_magnetic_domain_correction)
                     : nan);
    // The self-consistent columns above report the corrected operator AS CONFIGURED
    // (decision 1229 Q1 (ii)): 1 when it carried a NON-ZERO dK, i.e. a patch carries the
    // magnetic term after the dual gate's exclusions (decision 1235; `magnetic patches`
    // counts the carrying patches, the excluded ones are counted apart).
    table["magnetic_term_applied"]
        << (result.has_self_consistent
                ? (magnetic && confidence.magnetic_patch_count > 0 ? 1.0 : 0.0)
                : nan);
    PrintMarked(
        table, "self_consistent_frequency_magnetic_re",
        result.has_self_consistent && magnetic && confidence.magnetic_patch_count > 0
            ? units.Dimensionalize<VT::FREQUENCY>(result.self_consistent_frequency.real()) /
                  (2 * M_PI)
            : nan,
        measurement_marked);
  }
  {
    // Lane L1: the worst cell over Types (the largest finite uncovered share; the worst
    // statuses by CoverageStatusRank / RegimeStatusRank).
    double worst_share = nan;
    for (const auto &cell : coverage.cells)
    {
      if (std::isfinite(cell.uncovered_share) && !(cell.uncovered_share <= worst_share))
      {
        worst_share = cell.uncovered_share;
      }
    }
    const auto *worst_status = coverage.WorstCoverageCell();
    const auto *worst_regime = coverage.WorstRegimeCell();
    table["coverage_worst_uncovered_share"] << worst_share;
    table["coverage_worst_status"]
        << std::string(worst_status ? CoverageStatusName(worst_status->status)
                                    : CoverageStatusName(CoverageStatus::NOT_MEASURED));
    table["coverage_worst_regime"]
        << std::string(worst_regime ? RegimeStatusName(worst_regime->regime)
                                    : RegimeStatusName(RegimeStatus::NO_MODEL));
    if (marking)
    {
      table["coverage_marked_uncovered_share"]
          << coverage.MarkedUncoveredEnergy() / result.raw_normalization_energy;
      table["coverage_mark_measurement"] << (measurement_marked ? 1.0 : 0.0);
    }
  }
  surface_response_confidence->WriteFullTableTrunc();
  if (surface_response_op->HasMagneticLeaves())
  {
    // COL-1 (R-MAG-14): the per-model energy rows of this measurement (the raw field's
    // postprocessing; models in model-index order as the model contributions).
    if (result.has_postprocessed && confidence.has_magnetic)
    {
      const double frequency_GHz =
          units.Dimensionalize<VT::FREQUENCY>(measurement_cache.freq.real()) / (2 * M_PI);
      const int source = solver_t == ProblemType::DRIVEN
                             ? excitation
                             : static_cast<int>(std::lround(output_index));
      for (std::size_t m = 0; m < confidence.model_contributions.size() &&
                              m < confidence.model_magnetic_domain_corrections.size() &&
                              m < confidence.model_magnetic_thin_energies.size();
           m++)
      {
        if (confidence.model_magnetic_thin_energies[m] == 0.0 &&
            confidence.model_magnetic_domain_corrections[m] == 0.0)
        {
          continue;  // a model without the term (no leaf, or no patch carrying it)
        }
        surface_response_magnetic_energy_rows.push_back(
            {confidence.model_contributions[m].model, source, frequency_GHz,
             units.Dimensionalize<VT::ENERGY>(
                 confidence.model_magnetic_domain_corrections[m]),
             units.Dimensionalize<VT::ENERGY>(confidence.model_magnetic_thin_energies[m])});
      }
    }
    PrintSurfaceResponseMagnetic();
  }
  PrintSurfaceResponseModelEnergies(output_index, excitation);
  if (surface_response_patch_energy_enabled)
  {
    PrintSurfaceResponsePatchEnergies(output_index, excitation);
  }
  PrintSurfaceResponsePortionEnergies(output_index, excitation);
  PrintSurfaceResponseModelRawEnergies(output_index, excitation);
  PrintSurfaceResponseExclusionEnergies(output_index, excitation);
  PrintSurfaceResponseUncoveredFeatureEnergies(output_index, excitation);
  PrintSurfaceResponseLoopFit(output_index, excitation);
  PrintSurfaceResponseCoverage(output_index, excitation);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseMagnetic() const
{
  // surface-response-magnetic.csv (DESIGN-v2 section 3; decision 1229 Q4), root only,
  // rewritten whole at every print: one "model" row per loaded model (the leaf's verdict)
  // and one "gate" row per dual R0-a reading so far (global patch order within a source).
  if (!surface_response_op || !surface_response_op->HasMagneticLeaves() ||
      !Mpi::Root(fem_op->GetComm()))
  {
    return;
  }
  const auto path = post_dir / "surface-response-magnetic.csv";
  // OI-MAGNETIC-CSV-WRITE-THROUGH-LINK (fix round P4 re-freeze; decision 1368; the
  // OI-ADAPT-SAVEMESH-OVERWRITE form): SaveIteration leaves a symlink at this path into the
  // archived iteration folder; truncating through it would overwrite the archived table
  // and leave no real top-level file. Remove the link first, then write a real file.
  mesh::RemoveSymlinkBeforeWrite(path);
  std::ofstream output(path, std::ios::trunc);
  MFEM_VERIFY(output, "Unable to open " << path.string() << " for writing!");
  // Phase 3 (a) COL-1 (R-MAG-14): two trailing columns + one "energy" row per model
  // carrying the term, per source and frequency / mode (the model's applied magnetic domain
  // correction and its within-R magnetic thin energy, J); the model and gate rows leave
  // them empty.
  output << "kind,model,name,status,rule,psd ratio,negative eigenvalues,same-mesh "
            "deviation,patches,patch,source,f (GHz),max ratio,excluded,reason,magnetic "
            "domain correction (J),magnetic thin energy (J)\n";
  auto Quote = [](const std::string &text)
  {
    std::string quoted = "\"";
    for (const char c : text)
    {
      quoted += c == '"' ? std::string("\"\"") : std::string(1, c);
    }
    return quoted + "\"";
  };
  for (const auto &record : surface_response_op->GetMagneticModelRecords())
  {
    output << fmt::format("model,{},{},{},{},{:.6e},{},{:.6e},{},,,,,,{},,\n", record.model,
                          Quote(record.name), record.status, record.rule, record.psd_ratio,
                          record.negative_eigenvalue_count, record.same_mesh_deviation,
                          record.patches, Quote(record.reason));
  }
  for (const auto &record : surface_response_op->GetMagneticGateRecords())
  {
    output << fmt::format("gate,{},,,,,,,,{},{},{:.8e},{:.6e},{},,,\n", record.model,
                          record.patch, record.source, record.frequency, record.max_ratio,
                          record.excluded ? 1 : 0);
  }
  for (const auto &row : surface_response_magnetic_energy_rows)
  {
    output << fmt::format("energy,{},,,,,,,,,{},{:.8e},,,,{:.15e},{:.15e}\n", row.model,
                          row.source, row.frequency_GHz, row.domain_correction,
                          row.thin_energy);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseExclusionEnergies(double output_index,
                                                                   int excitation) const
{
  // OI-COVERAGE-EXCLUSION-ENERGY (idfix DESIGN section 8; decisions 1064 (2) / 1074). Root
  // only (the caller returns on the other ranks). A plain text table (the class and type
  // columns are text): rows = the Maxwell key (f (GHz) / m), exc, class, type, the class's
  // portions and summed length (m) over every type, the raw field's within-R energy of the
  // class on the type's target interfaces (J), the type's interface raw energy (J) and the
  // share; rows accumulate over the prints like the other Maxwell tables. The electrostatic
  // driver writes the same columns with the `source` key.
  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  if (result.exclusion_energy_partition.empty())
  {
    return;
  }
  const auto path = post_dir / "surface-response-exclusion-energy.csv";
  if (!surface_response_exclusion_energy_started)
  {
    surface_response_exclusion_energy_started = true;
    // OI-ES-EXCLUSION-CSV-SYMLINK-OVERWRITE (fix round P2): a previous adaptive iteration's
    // SaveIteration left a symlink here into its archive folder; truncating through it
    // would overwrite the archived table and leave no real top-level file (the loop-fit
    // writer's form: remove first, then write a real file).
    fs::remove(path);
    std::ofstream output(path);
    MFEM_VERIFY(output, "Unable to open " << path.string() << " for writing!");
    output
        << (solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m")
        << ",exc,class,type,class portions,class length (m),raw energy (J),interface raw "
           "energy (J),share\n";
  }
  const auto &types = surface_response_op->GetGeometryFallbackInterfaceTypes();
  auto TypeName = [&](int interface)
  {
    const auto type = types.find(interface);
    return type == types.end() ? std::string("Default") : ToString(type->second);
  };
  std::map<std::string, double> interface_raw_by_type;
  std::map<std::string, std::map<std::string, double>> excluded_by_class_type;
  for (const auto &[interface, partition] : result.exclusion_energy_partition)
  {
    const std::string name = TypeName(interface);
    const auto data = result.interfaces.find(interface);
    interface_raw_by_type[name] +=
        data != result.interfaces.end() ? data->second.raw_energy : 0.0;
    for (const auto &[cls, energy] : partition)
    {
      excluded_by_class_type[cls][name] += energy;
    }
  }
  std::map<std::string, std::pair<int, double>> portions_by_class;
  for (const auto &portion : surface_response_exclusion_energy)
  {
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
    }
    for (const std::string &key : {portion.type, std::string("Total")})
    {
      auto &[count, length] = portions_by_class[key];
      count++;
      length += std::sqrt(length2);
    }
  }
  std::ofstream output(path, std::ios::app);
  MFEM_VERIFY(output, "Unable to open " << path.string() << " for appending!");
  for (const auto &[cls, by_type] : excluded_by_class_type)
  {
    const auto [count, length] = portions_by_class[cls];
    for (const auto &[name, energy] : by_type)
    {
      const double interface_raw = interface_raw_by_type[name];
      if constexpr (solver_t == ProblemType::DRIVEN)
      {
        output << fmt::format("{:.8e},{}", output_index, excitation);
      }
      else
      {
        output << fmt::format("{},{}", static_cast<int>(std::lround(output_index)),
                              excitation);
      }
      output << fmt::format(",{},{},{},{},{},{},{}\n", cls, name, count,
                            CsvNumber(units.Dimensionalize<VT::LENGTH>(length)),
                            CsvNumber(units.Dimensionalize<VT::ENERGY>(energy)),
                            CsvNumber(units.Dimensionalize<VT::ENERGY>(interface_raw)),
                            CsvNumber(interface_raw > 0.0 ? energy / interface_raw : 0.0));
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseUncoveredFeatureEnergies(
    double output_index, int excitation) const
{
  // Diagnostics.UncoveredFeatureEnergy (lane ctx002-miss, decision 1567 (b); diagnostic
  // only). Root only (the caller returns on the other ranks). A plain text table (the type
  // and features columns are text): rows = the Maxwell key (f (GHz) / m), exc, the
  // uncovered requirement's type, its feature label (an identification feature id, or the
  // '+'-joined ids where portions of several features cover one perimeter piece), the
  // label's portions, their summed length and bounding box (m), then per target interface
  // the raw field's within-R energy (J) whose nearest perimeter foot lies on the label's
  // pieces. Per (type, interface) the rows sum to the raw-field rows (evaluations 0 and 1)
  // of surface-response-uncovered-energy.csv. The file is removed at the first print
  // whatever the option (a stale table of an earlier run never reads as this run's).
  using VT = Units::ValueType;
  const auto path = post_dir / "surface-response-uncovered-feature-energy.csv";
  const bool first_print = !surface_response_uncovered_feature_started;
  if (first_print)
  {
    fs::remove(path);
    surface_response_uncovered_feature_started = true;
  }
  const auto &result = *surface_response_measurement;
  if (result.uncovered_feature.empty())
  {
    return;
  }
  const auto target_interfaces = surface_response_op->GetTargetInterfaces();
  if (first_print || !fs::exists(path))
  {
    std::ofstream output(path);
    MFEM_VERIFY(output, "Unable to open " << path.string() << " for writing!");
    output << (solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m")
           << ",exc,type,features,portions,length (m),x min (m),y min (m),z min (m),x max "
              "(m),y max (m),z max (m)";
    for (const int interface : target_interfaces)
    {
      output << fmt::format(",uncovered raw energy[{}] (J)", interface);
    }
    output << '\n';
  }
  // The portions per (type, feature): count, length and bounding box.
  struct PortionSummary
  {
    int portions = 0;
    double length = 0.0;
    std::array<double, 3> min{}, max{};
  };
  std::map<std::pair<std::string, int>, PortionSummary> summaries;
  for (const auto &portion : surface_response_uncovered)
  {
    auto &summary = summaries[{portion.type, portion.feature}];
    const bool first = summary.portions == 0;
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
      const double low = std::min(portion.p0[d], portion.p1[d]);
      const double high = std::max(portion.p0[d], portion.p1[d]);
      summary.min[d] = first ? low : std::min(summary.min[d], low);
      summary.max[d] = first ? high : std::max(summary.max[d], high);
    }
    summary.portions++;
    summary.length += std::sqrt(length2);
  }
  std::set<std::pair<std::string, std::string>> keys;
  for (const auto &[interface, energies] : result.uncovered_feature)
  {
    (void)interface;
    for (const auto &[key, energy] : energies)
    {
      (void)energy;
      keys.insert(key);
    }
  }
  std::ofstream output(path, std::ios::app);
  MFEM_VERIFY(output, "Unable to open " << path.string() << " for appending!");
  for (const auto &key : keys)
  {
    const auto &[type, label] = key;
    PortionSummary label_summary;
    std::size_t begin = 0;
    while (begin <= label.size())
    {
      const auto end = label.find('+', begin);
      const int feature = std::stoi(
          label.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
      const auto summary = summaries.find({type, feature});
      MFEM_VERIFY(summary != summaries.end(), "Uncovered feature label \""
                                                  << label << "\" of type " << type
                                                  << " names no uncovered portion!");
      for (int d = 0; d < 3; d++)
      {
        label_summary.min[d] = label_summary.portions == 0
                                   ? summary->second.min[d]
                                   : std::min(label_summary.min[d], summary->second.min[d]);
        label_summary.max[d] = label_summary.portions == 0
                                   ? summary->second.max[d]
                                   : std::max(label_summary.max[d], summary->second.max[d]);
      }
      label_summary.portions += summary->second.portions;
      label_summary.length += summary->second.length;
      if (end == std::string::npos)
      {
        break;
      }
      begin = end + 1;
    }
    if constexpr (solver_t == ProblemType::DRIVEN)
    {
      output << fmt::format("{:.8e},{}", output_index, excitation);
    }
    else
    {
      output << fmt::format("{},{}", static_cast<int>(std::lround(output_index)),
                            excitation);
    }
    output << fmt::format(
        ",{},{},{},{}", type, label, label_summary.portions,
        CsvNumber(units.Dimensionalize<VT::LENGTH>(label_summary.length)));
    for (const auto &bound : {label_summary.min, label_summary.max})
    {
      for (int d = 0; d < 3; d++)
      {
        output << "," << CsvNumber(units.Dimensionalize<VT::LENGTH>(bound[d]));
      }
    }
    for (const int interface : target_interfaces)
    {
      double energy = 0.0;
      const auto energies = result.uncovered_feature.find(interface);
      if (energies != result.uncovered_feature.end())
      {
        const auto entry = energies->second.find(key);
        energy = entry != energies->second.end() ? entry->second : 0.0;
      }
      output << "," << CsvNumber(units.Dimensionalize<VT::ENERGY>(energy));
    }
    output << '\n';
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseModelRawEnergies(double output_index,
                                                                  int excitation) const
{
  // PatchEnergy (decision 922 B-1; the 851 (C) export). Root only (the caller returns on
  // the other ranks). Plain text tables (the label column is text):
  // surface-response-model-raw-energy.csv rows = source, evaluation (0: the raw field),
  // label, then per postprocessing interface the raw within-R energy on the label's pieces
  // (0 for a non-target interface), then the Maxwell keys "f (GHz)" / "m" and "exc";
  // surface-response-model-raw-groups.csv rows = label, model, patches, weight x length
  // (m), written once (the claims are mesh-independent per configuration). Rows accumulate
  // over the prints like the other Maxwell tables.
  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  if (result.model_raw_partition.empty())
  {
    return;
  }
  const auto energy_path = post_dir / "surface-response-model-raw-energy.csv";
  const auto groups_path = post_dir / "surface-response-model-raw-groups.csv";
  const int source = solver_t == ProblemType::DRIVEN
                         ? excitation
                         : static_cast<int>(std::lround(output_index));
  if (!surface_response_model_raw_started)
  {
    surface_response_model_raw_started = true;
    // OI-ES-EXCLUSION-CSV-SYMLINK-OVERWRITE (fix round P2): never write through the symlink
    // SaveIteration left (see PrintSurfaceResponseExclusionEnergies).
    fs::remove(energy_path);
    fs::remove(groups_path);
    std::ofstream energy(energy_path);
    MFEM_VERIFY(energy, "Unable to open " << energy_path.string() << " for writing!");
    energy << "source,evaluation,label";
    for (const auto &[interface, data] : result.interfaces)
    {
      (void)data;
      energy << fmt::format(",raw energy[{}] (J)", interface);
    }
    energy << (solver_t == ProblemType::DRIVEN ? ",f (GHz),exc\n" : ",m,exc\n");
    std::ofstream groups(groups_path);
    MFEM_VERIFY(groups, "Unable to open " << groups_path.string() << " for writing!");
    groups << "label,model,patches,weight length (m)\n";
    for (const auto &member : surface_response_op->GetModelRawClaimPortions().members)
    {
      groups << fmt::format(
          "{},{},{},{}\n", member.label, member.model, member.patches,
          CsvNumber(units.Dimensionalize<VT::LENGTH>(member.weight_length)));
    }
  }
  std::set<std::string> labels;
  for (const auto &[interface, partition] : result.model_raw_partition)
  {
    (void)interface;
    for (const auto &[label, energy] : partition)
    {
      (void)energy;
      labels.insert(label);
    }
  }
  std::ofstream energy(energy_path, std::ios::app);
  MFEM_VERIFY(energy, "Unable to open " << energy_path.string() << " for appending!");
  for (const auto &label : labels)
  {
    energy << fmt::format("{},{},{}", source, 0, label);
    for (const auto &[interface, data] : result.interfaces)
    {
      (void)data;
      const auto partition = result.model_raw_partition.find(interface);
      double value = 0.0;
      if (partition != result.model_raw_partition.end())
      {
        const auto entry = partition->second.find(label);
        value = entry != partition->second.end() ? entry->second : 0.0;
      }
      energy << "," << CsvNumber(units.Dimensionalize<VT::ENERGY>(value));
    }
    if constexpr (solver_t == ProblemType::DRIVEN)
    {
      energy << fmt::format(",{:.8e},{}\n", output_index, excitation);
    }
    else
    {
      energy << fmt::format(",{},{}\n", static_cast<int>(std::lround(output_index)),
                            excitation);
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseLoopFit(double output_index,
                                                         int excitation) const
{
  // OI-MX-LOOPGATE (decision 955). Root only (the caller returns on the other ranks). Plain
  // text tables (the label column is text). surface-response-loop-fit.csv rows = source,
  // evaluation (0 fixed trace / 1 fixed flux on the raw field; 2 self-consistent on the
  // corrected field), the Maxwell key "f (GHz)" / "m" and "exc", max loop share,
  // response-weighted loop share, remainder patches, fallback patches, weighted claim
  // length (m) (sum over the remainder patches of rho_p x the raw-claim length), withheld
  // domain correction (J) (the fallback patches' domain defects of that closure).
  // surface-response-loop-remainder-energy.csv rows = source, evaluation, model row
  // ("m<idx>" / "Total"), patches, weighted claim length (m), the remainder's raw energy
  // per target interface (J), then the keys; a clean window writes the Total row with 0.
  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  const auto fit_path = post_dir / "surface-response-loop-fit.csv";
  const auto remainder_path = post_dir / "surface-response-loop-remainder-energy.csv";
  const int source = solver_t == ProblemType::DRIVEN
                         ? excitation
                         : static_cast<int>(std::lround(output_index));
  const auto target_interfaces = surface_response_op->GetTargetInterfaces();
  const bool first_print = !surface_response_loop_fit_started;
  if (first_print)
  {
    surface_response_loop_fit_started = true;
    fs::remove(fit_path);
    fs::remove(remainder_path);
  }
  if (!result.has_postprocessed && !result.has_self_consistent)
  {
    return;
  }
  auto Keys = [&]()
  {
    if constexpr (solver_t == ProblemType::DRIVEN)
    {
      return fmt::format("{},{}", CsvNumber(output_index), excitation);
    }
    else
    {
      return fmt::format("{},{}", static_cast<int>(std::lround(output_index)), excitation);
    }
  };
  const char *key_header = solver_t == ProblemType::DRIVEN ? "f (GHz),exc" : "m,exc";
  // OI-MAXWELL-CC-CSV-HEADERLESS-AFTER-FIRST-PRINT (lane L1, decision 1634 (6)): the
  // headers go with the first rows (a run whose first print has neither evaluation wrote
  // both tables without them).
  const bool write_header =
      surface_response_tables_with_header.insert(fit_path.string()).second;
  std::ofstream fit(fit_path, write_header ? std::ios::out | std::ios::trunc
                                           : std::ios::out | std::ios::app);
  MFEM_VERIFY(fit, "Unable to open " << fit_path.string() << " for writing!");
  std::ofstream remainder(remainder_path, write_header ? std::ios::out | std::ios::trunc
                                                       : std::ios::out | std::ios::app);
  MFEM_VERIFY(remainder, "Unable to open " << remainder_path.string() << " for writing!");
  if (write_header)
  {
    fit << "source,evaluation," << key_header
        << ",max loop share,response-weighted loop share,remainder patches,fallback "
           "patches,weighted claim length (m),withheld domain correction (J)\n";
    remainder << "source,evaluation,model,patches,weighted claim length (m)";
    for (const int interface : target_interfaces)
    {
      remainder << fmt::format(",remainder energy[{}] (J)", interface);
    }
    remainder << "," << key_header << "\n";
  }
  auto Rows = [&](const SurfaceResponseOperator::MaxwellResponse &response)
  { return surface_response_op->GetLoopRemainderPortions(response.loop_remainders).rows; };
  auto Write = [&](int evaluation, const SurfaceResponseOperator::MaxwellResponse &response,
                   const std::map<int, std::map<std::string, double>> &energies)
  {
    const auto rows = Rows(response);
    SurfaceResponseOperator::LoopRemainderPortions::Row total;
    for (const auto &[label, row] : rows)
    {
      (void)label;
      total.patches += row.patches;
      total.weighted_length += row.weighted_length;
    }
    fit << fmt::format(
        "{},{},{},{},{},{},{},{},{}\n", source, evaluation, Keys(),
        CsvNumber(response.maximum_loop_fit_share),
        CsvNumber(response.response_weighted_loop_fit_share),
        response.loop_remainder_patch_count, response.loop_fallback_patch_count,
        CsvNumber(units.Dimensionalize<VT::LENGTH>(total.weighted_length)),
        CsvNumber(units.Dimensionalize<VT::ENERGY>(
            evaluation == 1 ? response.loop_fallback_domain_correction_fixed_flux
                            : response.loop_fallback_domain_correction)));
    auto Line = [&](const std::string &label,
                    const SurfaceResponseOperator::LoopRemainderPortions::Row &row)
    {
      remainder << fmt::format(
          "{},{},{},{},{}", source, evaluation, label, row.patches,
          CsvNumber(units.Dimensionalize<VT::LENGTH>(row.weighted_length)));
      for (const int interface : target_interfaces)
      {
        double energy = 0.0;
        const auto it = energies.find(interface);
        if (it != energies.end())
        {
          if (label == "Total")
          {
            for (const auto &[key, value] : it->second)
            {
              (void)key;
              energy += value;
            }
          }
          else
          {
            const auto entry = it->second.find(label);
            energy = entry != it->second.end() ? entry->second : 0.0;
          }
        }
        remainder << "," << CsvNumber(units.Dimensionalize<VT::ENERGY>(energy));
      }
      remainder << "," << Keys() << "\n";
    };
    for (const auto &[label, row] : rows)
    {
      Line(label, row);
    }
    Line("Total", total);
  };
  if (result.has_postprocessed)
  {
    Write(0, result.confidence, result.loop_remainder);
    Write(1, result.confidence, result.loop_remainder);
  }
  if (result.has_self_consistent)
  {
    Write(2, result.self_consistent_confidence, result.self_consistent_loop_remainder);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseCoverage(double output_index,
                                                          int excitation) const
{
  // Lane L1 (decision 1634 (2); DESIGN-v2 s4.2; rulings 1641 / 1649 / 1654 (III)). Root
  // only (the caller returns on the other ranks). Plain text tables (the Type, key and
  // status columns are text; WriteCoverageRows / WriteCoverageByKeyRows). Both files are
  // removed at the first print whatever is written (a stale table of an earlier run never
  // reads as this run's); each writes its header with its first rows.
  const auto coverage_path = post_dir / "surface-response-coverage.csv";
  const auto by_key_path = post_dir / "surface-response-uncovered-energy-by-key.csv";
  // Lane rescensus-product (decision 1817): the census table beside them (top level).
  const auto resolution_path = post_dir / "surface-response-resolution.csv";
  if (!surface_response_coverage_started)
  {
    fs::remove(coverage_path);
    fs::remove(by_key_path);
    fs::remove(resolution_path);
    surface_response_coverage_started = true;
  }
  const auto &coverage = surface_response_measurement->coverage_measurement;
  const std::string keys =
      solver_t == ProblemType::DRIVEN
          ? fmt::format("{},{}", CsvNumber(output_index),
                        excitation)  // the other tables' form
          : fmt::format("{},{}", static_cast<int>(std::lround(output_index)), excitation);
  const std::string key_header = solver_t == ProblemType::DRIVEN ? "f (GHz),exc" : "m,exc";
  // Ruling 1738: the coverage cells' interfaces (the served and the configured not served).
  const auto target_interfaces = surface_response_op->GetCoverageInterfaces();
  for (const auto &path : {coverage_path, by_key_path})
  {
    const bool write_header =
        surface_response_tables_with_header.insert(path.string()).second;
    std::ofstream output(path, write_header ? std::ios::out | std::ios::trunc
                                            : std::ios::out | std::ios::app);
    MFEM_VERIFY(output, "Unable to open " << path.string() << " for writing!");
    if (path == coverage_path)
    {
      output << (write_header ? CoverageTableHeader(key_header) : "");
      WriteCoverageRows(output, keys, coverage, units);
    }
    else
    {
      output << (write_header ? CoverageByKeyTableHeader(key_header, target_interfaces)
                              : "");
      WriteCoverageByKeyRows(output, keys, coverage, target_interfaces,
                             CoveragePortionsByFeature(surface_response_uncovered), units);
    }
  }
  std::vector<std::pair<std::string, double>> cell_keys;
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    cell_keys = {{"f (GHz)", output_index}, {"exc", static_cast<double>(excitation)}};
  }
  else
  {
    cell_keys = {{"m", std::round(output_index)}};
  }
  surface_response_op->GetCoverageRunRecord().Add(cell_keys, coverage);
  {
    const auto &cells = surface_response_measurement->resolution_cells;
    const bool write_header =
        surface_response_tables_with_header.insert(resolution_path.string()).second;
    std::ofstream output(resolution_path, write_header ? std::ios::out | std::ios::trunc
                                                       : std::ios::out | std::ios::app);
    MFEM_VERIFY(output, "Unable to open " << resolution_path.string() << " for writing!");
    output << (write_header ? ResolutionTableHeader(key_header) : "");
    WriteResolutionRows(output, keys, cells);
    for (const auto &cell : cells)
    {
      surface_response_op->GetResolutionRunRecord().entries.push_back({cell_keys, cell});
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponseModelEnergies(double output_index,
                                                               int excitation) const
{
  // Root only (the caller returns on the other ranks); the ModelContribution entries are
  // replicated by the response's batched reduction. Decisions 845 O-C / 846: the
  // electrostatic driver's surface-response-model-energy.csv column set EXACTLY and in
  // order as a strict prefix (`source`, `evaluation`, `model`, `patch count`,
  // `patch weight`, `domain correction (J)`, `fabricated surface energy[k] (J)` for every
  // postprocessing interface, 0 for a non-target interface), then the Maxwell keys
  // "f (GHz)" / "m" and "exc". `source` is the excitation index on a driven solve and
  // the mode index on an eigenmode / boundary-mode solve. Rows per print as the
  // electrostatic driver writes them per source: every model at evaluation 0 (fixed trace)
  // and 1 (fixed flux) on the raw field when postprocessing ran; every model at
  // evaluation 2 on the self-consistent corrected field when it exists (its `domain
  // correction` is the fixed-trace defect the self-consistent normalization carries = the
  // electrostatic driver's evaluation-2 row under the default TranslationalDomainCorrection
  // FixedTrace; its energies the fixed-trace ones). Rows accumulate over the prints: a
  // reader of a multi-frequency sweep keys the rows by "f (GHz)" / "m" and "exc".
  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  if (!surface_response_model_energy)
  {
    surface_response_model_energy =
        std::make_unique<TableWithCSVFile>(post_dir / "surface-response-model-energy.csv");
    auto &table = surface_response_model_energy->table;
    table.insert(Column::Integer("source", "source"));
    table.insert(Column::Integer("evaluation", "evaluation"));
    table.insert(Column::Integer("model", "model"));
    table.insert(Column::Integer("patch_count", "patch count"));
    table.insert("patch_weight", "patch weight");
    table.insert("domain_correction", "domain correction (J)");
    for (const auto &[interface, data] : result.interfaces)
    {
      (void)data;
      table.insert(fmt::format("interface_{}", interface),
                   fmt::format("fabricated surface energy[{}] (J)", interface));
    }
    table.insert("idx", solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m");
    table.insert("excitation", "exc");
    table["excitation"].print_as_int = true;
    if constexpr (solver_t == ProblemType::EIGENMODE ||
                  solver_t == ProblemType::BOUNDARYMODE)
    {
      table["idx"].print_as_int = true;
    }
  }
  auto &table = surface_response_model_energy->table;
  const int source = solver_t == ProblemType::DRIVEN
                         ? excitation
                         : static_cast<int>(std::lround(output_index));
  auto Append = [&](int evaluation, const SurfaceResponseOperator::ModelContribution &data,
                    bool fixed_flux)
  {
    table["source"] << source;
    table["evaluation"] << evaluation;
    table["model"] << data.model;
    table["patch_count"] << data.patch_count;
    table["patch_weight"] << data.patch_weight;
    table["domain_correction"] << units.Dimensionalize<VT::ENERGY>(
        fixed_flux ? data.domain_correction_fixed_flux : data.domain_correction);
    const auto &energies = fixed_flux ? data.fabricated_surface_energy_fixed_flux
                                      : data.fabricated_surface_energy;
    for (const auto &[interface, interface_data] : result.interfaces)
    {
      (void)interface_data;
      const auto energy = energies.find(interface);
      table[fmt::format("interface_{}", interface)] << units.Dimensionalize<VT::ENERGY>(
          energy != energies.end() ? energy->second : 0.0);
    }
    table["idx"] << output_index;
    table["excitation"] << excitation;
  };
  if (result.has_postprocessed)
  {
    for (const auto &contribution : result.confidence.model_contributions)
    {
      Append(0, contribution, false);
      Append(1, contribution, true);
    }
  }
  if (result.has_self_consistent)
  {
    for (const auto &contribution : result.self_consistent_confidence.model_contributions)
    {
      Append(2, contribution, false);
    }
  }
  surface_response_model_energy->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponsePatchEnergies(double output_index,
                                                               int excitation) const
{
  // PatchEnergy (decisions 913 / 920 / 922 B). Root only (the caller returns on the other
  // ranks); the records are replicated by the response's gathers. The electrostatic
  // driver's per-patch rows per print as it writes them per source: every applied patch at
  // evaluation 0 (fixed trace) and 1 (fixed flux) on the raw field, at evaluation 2 on the
  // self-consistent corrected field when it exists; `patch` = the 1-based placed patch
  // index of the traces file; the origin from the patch table; then the Maxwell keys and
  // the patch's max contour-loop residual (decision 920; 0 without a path record). The
  // spatial traces (real part in the electrostatic column, the imaginary part after the
  // prefix) and the per-path loop residuals follow the same key form. The rows of each
  // print are appended to the files and released (TableWithCSVFile::AppendRowsAndClear;
  // lane fix-response-memory C4: the file after print m = header + the rows of prints 1..m,
  // byte for byte the former whole-table rewrite per print).
  using VT = Units::ValueType;
  const auto &result = *surface_response_measurement;
  std::map<int, const SurfaceResponseOperator::PatchAssignment *> assignments;
  for (const auto &assignment : surface_response_op->GetPatchAssignments())
  {
    assignments.emplace(assignment.global_index, &assignment);
  }
  const int source = solver_t == ProblemType::DRIVEN
                         ? excitation
                         : static_cast<int>(std::lround(output_index));
  const char *index_label = solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m";
  if (!surface_response_patch_energy)
  {
    surface_response_patch_energy =
        std::make_unique<TableWithCSVFile>(surface_response_op->DiagnosticsPath(
            post_dir, "surface-response-patch-energy.csv"));
    auto &table = surface_response_patch_energy->table;
    table.insert(Column::Integer("source", "source"));
    table.insert(Column::Integer("evaluation", "evaluation"));
    table.insert(Column::Integer("patch", "patch"));
    table.insert(Column::Integer("model", "model"));
    table.insert(Column::Integer("feature", "feature"));
    table.insert("weight", "weight");
    table.insert("cell_begin", "cell begin (m)");
    table.insert("cell_end", "cell end (m)");
    for (const char *coordinate : {"x", "y", "z"})
    {
      table.insert(fmt::format("origin_{}", coordinate),
                   fmt::format("origin {} (m)", coordinate));
    }
    table.insert("domain_correction", "domain correction (J)");
    for (const auto &[interface, data] : result.interfaces)
    {
      (void)data;
      table.insert(fmt::format("interface_{}", interface),
                   fmt::format("fabricated surface energy[{}] (J)", interface));
    }
    if (surface_response_op->HasMagneticLeaves())
    {
      // COL-2 (Phase 3 (a), R-MAG-14): the patch's applied magnetic domain correction and
      // its within-R magnetic thin energy (NaN on a patch without a magnetic leaf).
      table.insert("magnetic_domain_correction", "magnetic domain correction (J)");
      table.insert("magnetic_thin_energy", "magnetic thin energy (J)");
    }
    if (surface_response_within_r_volume_energy_enabled)
    {
      // COL-3 (Phase 3 (b), R-MAG-15): the raw field's energies in the patch's own box by
      // the volume integral at quadrature resolution, and the box volume it resolved (the
      // parent of any per-volume ratio); NaN on a patch without a box.
      table.insert("within_r_magnetic_energy", "within-R magnetic energy (J)");
      table.insert("within_r_electric_energy", "within-R electric energy (J)");
      table.insert("within_r_volume", "within-R volume (m^3)");
    }
    table.insert("idx", index_label);
    table.insert("excitation", "exc");
    table.insert("max_loop_residual", "max loop residual");
    table.insert("loop_share", "loop share");
    table.insert("loop_fallback", "loop fallback");
    table.insert("loop_weight", "loop weight (J)");
    table["excitation"].print_as_int = true;
    table["loop_fallback"].print_as_int = true;
    if constexpr (solver_t == ProblemType::EIGENMODE ||
                  solver_t == ProblemType::BOUNDARYMODE)
    {
      table["idx"].print_as_int = true;
    }
  }
  {
    auto &table = surface_response_patch_energy->table;
    // OI-PATCH-ENERGY-CSV-QUADRATIC-ROOT-WRITER (decision 1544): the per-patch max loop
    // residuals come built once per response (ForEachMaxwellPatchEnergyRow); they were
    // rebuilt over all of the response's loop records for every row (rows x loops).
    auto Append =
        [&](int evaluation, const SurfaceResponseOperator::MaxwellResponse &response,
            const std::map<int, double> &max_loops,
            const SurfaceResponseOperator::PatchContribution &data, bool fixed_flux)
    {
      const auto assignment = assignments.find(data.patch);
      MFEM_VERIFY(assignment != assignments.end(),
                  "Per-patch Maxwell response energy of an unknown patch " << data.patch + 1
                                                                           << "!");
      table["source"] << source;
      table["evaluation"] << evaluation;
      table["patch"] << data.patch + 1;
      table["model"] << data.model;
      table["feature"] << data.feature;
      table["weight"] << data.weight;
      table["cell_begin"] << units.Dimensionalize<VT::LENGTH>(data.cell[0]);
      table["cell_end"] << units.Dimensionalize<VT::LENGTH>(data.cell[1]);
      for (int d = 0; d < 3; d++)
      {
        table[fmt::format("origin_{}", "xyz"[d])]
            << units.Dimensionalize<VT::LENGTH>(assignment->second->origin[d]);
      }
      table["domain_correction"] << units.Dimensionalize<VT::ENERGY>(
          fixed_flux ? data.domain_correction_fixed_flux : data.domain_correction);
      const auto &energies = fixed_flux ? data.fabricated_surface_energy_fixed_flux
                                        : data.fabricated_surface_energy;
      for (const auto &[interface, interface_data] : result.interfaces)
      {
        (void)interface_data;
        const auto energy = energies.find(interface);
        table[fmt::format("interface_{}", interface)] << units.Dimensionalize<VT::ENERGY>(
            energy != energies.end() ? energy->second : 0.0);
      }
      if (surface_response_op->HasMagneticLeaves())
      {
        const auto magnetic = std::lower_bound(response.patch_magnetic_energies.begin(),
                                               response.patch_magnetic_energies.end(),
                                               data.patch, [](const auto &entry, int patch)
                                               { return entry.patch < patch; });
        const bool found = magnetic != response.patch_magnetic_energies.end() &&
                           magnetic->patch == data.patch;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        table["magnetic_domain_correction"]
            << (found ? units.Dimensionalize<VT::ENERGY>(magnetic->domain_correction)
                      : nan);
        table["magnetic_thin_energy"]
            << (found ? units.Dimensionalize<VT::ENERGY>(magnetic->thin_energy) : nan);
      }
      if (surface_response_within_r_volume_energy_enabled)
      {
        const auto &volumes = result.within_r_volume_energies;
        const auto entry = std::lower_bound(volumes.begin(), volumes.end(), data.patch,
                                            [](const auto &record, int patch)
                                            { return record.patch < patch; });
        const bool found =
            entry != volumes.end() && entry->patch == data.patch && entry->defined;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double length_scale = units.Dimensionalize<VT::LENGTH>(1.0);
        table["within_r_magnetic_energy"]
            << (found ? units.Dimensionalize<VT::ENERGY>(entry->magnetic) : nan);
        table["within_r_electric_energy"]
            << (found ? units.Dimensionalize<VT::ENERGY>(entry->electric) : nan);
        table["within_r_volume"]
            << (found ? entry->volume * length_scale * length_scale * length_scale : nan);
      }
      table["idx"] << output_index;
      table["excitation"] << excitation;
      const auto max_loop = max_loops.find(data.patch);
      table["max_loop_residual"] << (max_loop != max_loops.end() ? max_loop->second : 0.0);
      // OI-MX-LOOPGATE (decision 955): the patch's circulating share, the fallback flag and
      // the detector weight (its response energy on the chain trace: the loop-RMS weight).
      const auto share = std::lower_bound(
          response.patch_loop_shares.begin(), response.patch_loop_shares.end(), data.patch,
          [](const auto &entry, int patch) { return entry.patch < patch; });
      MFEM_VERIFY(share != response.patch_loop_shares.end() && share->patch == data.patch,
                  "Per-patch Maxwell loop share of patch " << data.patch + 1
                                                           << " missing!");
      table["loop_share"] << share->share;
      table["loop_fallback"] << share->fallback;
      table["loop_weight"] << units.Dimensionalize<VT::ENERGY>(share->detector_energy);
    };
    ForEachMaxwellPatchEnergyRow(result.confidence, result.has_postprocessed,
                                 result.self_consistent_confidence,
                                 result.has_self_consistent, Append);
    surface_response_patch_energy->AppendRowsAndClear();
  }
  if (result.has_postprocessed && !result.confidence.spatial_traces.empty())
  {
    if (!surface_response_traces)
    {
      surface_response_traces = std::make_unique<TableWithCSVFile>(
          surface_response_op->DiagnosticsPath(post_dir, "surface-response-traces.csv"));
      auto &table = surface_response_traces->table;
      table.insert(Column::Integer("source", "i"));
      table.insert(Column::Integer("patch", "patch"));
      table.insert(Column::Integer("model", "model"));
      table.insert(Column::Integer("coefficient", "coefficient"));
      table.insert(Column::Integer("conductor_state", "conductor state"));
      table.insert("value", "value (V)");
      table.insert("value_imag", "value imag (V)");
      table.insert("idx", index_label);
      table.insert("excitation", "exc");
      table["excitation"].print_as_int = true;
      if constexpr (solver_t == ProblemType::EIGENMODE ||
                    solver_t == ProblemType::BOUNDARYMODE)
      {
        table["idx"].print_as_int = true;
      }
    }
    auto &table = surface_response_traces->table;
    for (const auto &entry : result.confidence.spatial_traces)
    {
      for (std::size_t i = 0; i < entry.coefficients.size(); i++)
      {
        table["source"] << source;
        table["patch"] << entry.patch + 1;
        table["model"] << entry.model;
        table["coefficient"] << static_cast<int>(i) + 1;
        table["conductor_state"] << (static_cast<int>(i) >= entry.contour_size ? 1 : 0);
        table["value"] << units.Dimensionalize<VT::VOLTAGE>(entry.coefficients[i].real());
        table["value_imag"] << units.Dimensionalize<VT::VOLTAGE>(
            entry.coefficients[i].imag());
        table["idx"] << output_index;
        table["excitation"] << excitation;
      }
    }
    surface_response_traces->AppendRowsAndClear();
  }
  if (result.has_postprocessed && !result.confidence.path_loop_residuals.empty())
  {
    if (!surface_response_patch_loops)
    {
      surface_response_patch_loops =
          std::make_unique<TableWithCSVFile>(surface_response_op->DiagnosticsPath(
              post_dir, "surface-response-patch-loops.csv"));
      auto &table = surface_response_patch_loops->table;
      table.insert(Column::Integer("source", "source"));
      table.insert(Column::Integer("patch", "patch"));
      table.insert(Column::Integer("model", "model"));
      table.insert(Column::Integer("path", "path"));
      table.insert(Column::Integer("kind", "kind"));
      table.insert("residual", "loop residual");
      table.insert("scale", "loop scale (V)");
      table.insert("idx", index_label);
      table.insert("excitation", "exc");
      table["excitation"].print_as_int = true;
      if constexpr (solver_t == ProblemType::EIGENMODE ||
                    solver_t == ProblemType::BOUNDARYMODE)
      {
        table["idx"].print_as_int = true;
      }
    }
    auto &table = surface_response_patch_loops->table;
    for (const auto &loop : result.confidence.path_loop_residuals)
    {
      table["source"] << source;
      table["patch"] << loop.patch + 1;
      table["model"] << loop.model;
      table["path"] << loop.path;
      table["kind"] << loop.kind;
      table["residual"] << loop.residual;
      table["scale"] << units.Dimensionalize<VT::VOLTAGE>(loop.scale);
      table["idx"] << output_index;
      table["excitation"] << excitation;
    }
    surface_response_patch_loops->AppendRowsAndClear();
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::PrintSurfaceResponsePortionEnergies(double output_index,
                                                                 int excitation) const
{
  // The electrostatic driver's per-type portion tables (decision 394 F2 / F-DB-a) for the
  // Maxwell correction: `type` with a Total row, `portions`, `length (m)`, the extra
  // column, then per target interface the raw within-R energy and its share of the
  // evaluation's corrected interface energy (0 fixed trace, 1 fixed flux on the raw field,
  // 2 self-consistent on the corrected field). A plain writer (the type column is text).
  // The three tables are removed UNCONDITIONALLY at the first print before the conditional
  // writes (decision 473 MAJOR-3): a stale table of an earlier run in the output directory
  // would otherwise read as a term this solve did not have; later prints append rows.
  using VT = Units::ValueType;
  const fs::path uncovered_path = post_dir / "surface-response-uncovered-energy.csv";
  const fs::path domain_boundary_path =
      post_dir / "surface-response-domain-boundary-energy.csv";
  const fs::path spatial_cluster_priority_path =
      post_dir / "surface-response-spatial-cluster-priority-energy.csv";
  const fs::path conductor_consistency_path =
      post_dir / "surface-response-conductor-consistency-energy.csv";
  const bool first_print = !surface_response_portion_tables_started;
  if (first_print)
  {
    fs::remove(uncovered_path);
    fs::remove(domain_boundary_path);
    fs::remove(spatial_cluster_priority_path);
    fs::remove(conductor_consistency_path);
    surface_response_portion_tables_started = true;
  }
  const auto &result = *surface_response_measurement;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  struct PortionTypeRow
  {
    int portions = 0;
    double length = 0.0;
    double extra = 0.0;
  };
  auto PortionLength = [](const auto &portion)
  {
    double length2 = 0.0;
    for (int d = 0; d < 3; d++)
    {
      length2 += (portion.p1[d] - portion.p0[d]) * (portion.p1[d] - portion.p0[d]);
    }
    return std::sqrt(length2);
  };
  const auto target_interfaces = surface_response_op->GetTargetInterfaces();
  auto WritePortionEnergyTable =
      [&](const fs::path &path,
          const std::map<std::string, PortionTypeRow> &portions_by_type,
          const std::string &extra_header, const std::string &energy_label,
          const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> &raw_energies,
          const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> &corrected_energies)
  {
    PortionTypeRow portions_total;
    for (const auto &[type, entry] : portions_by_type)
    {
      (void)type;
      portions_total.portions += entry.portions;
      portions_total.length += entry.length;
      portions_total.extra += entry.extra;
    }
    // OI-MAXWELL-CC-CSV-HEADERLESS-AFTER-FIRST-PRINT (lane L1, decision 1634 (6)): the
    // header goes with the table's first rows (a table first written at a later print, e.g.
    // the conductor-consistency table of a run whose first exclusion comes at a later mode
    // or frequency, was appended without one). Every table first written at the first print
    // is unchanged.
    const bool write_header =
        surface_response_tables_with_header.insert(path.string()).second;
    std::ofstream output(path, write_header ? std::ios::out | std::ios::trunc
                                            : std::ios::out | std::ios::app);
    MFEM_VERIFY(output, "Unable to open " << path.filename().string() << " for writing!");
    if (write_header)
    {
      output << (solver_t == ProblemType::DRIVEN ? "f (GHz)" : "m")
             << ",exc,evaluation,type,portions,length (m)," << extra_header;
      for (const int interface : target_interfaces)
      {
        output << fmt::format(",{} raw energy[{}] (J),{} share[{}]", energy_label,
                              interface, energy_label, interface);
      }
      output << '\n';
    }
    const bool extra_is_length = extra_header.find("(m)") != std::string::npos;
    auto CorrectedEnergy = [&](int evaluation, int interface)
    {
      const auto it = result.interfaces.find(interface);
      if (it == result.interfaces.end())
      {
        return nan;
      }
      return evaluation == 0   ? it->second.corrected_energy
             : evaluation == 1 ? it->second.corrected_energy_fixed_flux
                               : it->second.self_consistent_energy;
    };
    auto Append =
        [&](int evaluation, const std::string &type, const PortionTypeRow &portions,
            const std::map<int, SurfacePostOperator::UncoveredEdgeEnergy> &energies,
            bool available)
    {
      if constexpr (solver_t == ProblemType::DRIVEN)
      {
        output << CsvNumber(output_index);
      }
      else
      {
        output << static_cast<int>(std::lround(output_index));
      }
      output << fmt::format(",{},{},{},{},{},", excitation, evaluation, type,
                            portions.portions,
                            CsvNumber(units.Dimensionalize<VT::LENGTH>(portions.length)));
      if (extra_is_length)
      {
        output << CsvNumber(units.Dimensionalize<VT::LENGTH>(portions.extra));
      }
      else
      {
        output << fmt::format("{}", static_cast<int>(std::lround(portions.extra)));
      }
      for (const int interface : target_interfaces)
      {
        double energy = nan, share = nan;
        const auto it = energies.find(interface);
        if (available && it != energies.end())
        {
          if (type == "Total")
          {
            energy = it->second.energy;
          }
          else
          {
            const auto by_type = it->second.by_type.find(type);
            energy = by_type != it->second.by_type.end() ? by_type->second : 0.0;
          }
          const double corrected = CorrectedEnergy(evaluation, interface);
          share = corrected > 0.0 ? energy / corrected : nan;
          energy = units.Dimensionalize<VT::ENERGY>(energy);
        }
        output << fmt::format(",{},{}", CsvNumber(energy), CsvNumber(share));
      }
      output << '\n';
    };
    for (const auto &[type, portions] : portions_by_type)
    {
      Append(0, type, portions, raw_energies, result.has_postprocessed);
      Append(1, type, portions, raw_energies, result.has_postprocessed);
      Append(2, type, portions, corrected_energies, result.has_self_consistent);
    }
    Append(0, "Total", portions_total, raw_energies, result.has_postprocessed);
    Append(1, "Total", portions_total, raw_energies, result.has_postprocessed);
    Append(2, "Total", portions_total, corrected_energies, result.has_self_consistent);
  };
  if (!surface_response_uncovered.empty())
  {
    std::map<std::string, PortionTypeRow> portions_by_type;
    for (const auto &portion : surface_response_op->GetUncoveredPortions())
    {
      auto &entry = portions_by_type[portion.topology];
      entry.portions++;
      entry.length += PortionLength(portion);
    }
    for (const auto &clip : surface_response_op->GetUncoveredSpatialSupportClipping().clips)
    {
      portions_by_type[clip.topology].extra += clip.length;
    }
    WritePortionEnergyTable(uncovered_path, portions_by_type,
                            "clipped by spatial support (m)", "uncovered", result.uncovered,
                            result.self_consistent_uncovered);
  }
  if (!surface_response_domain_boundary.empty())
  {
    // F-DB-a: the DomainBoundary share per type (type = the excluded patch's model
    // topology; the geometric cell count = translational own-edge intervals).
    std::map<std::string, PortionTypeRow> portions_by_type;
    const auto &record = surface_response_op->GetDomainBoundaryPortionRecord();
    for (std::size_t i = 0; i < record.portions.size(); i++)
    {
      auto &entry = portions_by_type[record.portions[i].topology];
      entry.portions++;
      entry.length += PortionLength(record.portions[i]);
      entry.extra += record.translational[i] ? 1.0 : 0.0;
    }
    WritePortionEnergyTable(domain_boundary_path, portions_by_type, "geometric cells",
                            "domainboundary", result.domain_boundary,
                            result.self_consistent_domain_boundary);
  }
  if (!surface_response_spatial_cluster_priority.empty())
  {
    // OI-SCP-HOLE (decision 922): the SpatialClusterPriority share per type (the
    // DomainBoundary table's columns).
    std::map<std::string, PortionTypeRow> portions_by_type;
    const auto &record =
        surface_response_op->GetSpatialClusterPriorityExclusions().portions;
    for (std::size_t i = 0; i < record.portions.size(); i++)
    {
      auto &entry = portions_by_type[record.portions[i].topology];
      entry.portions++;
      entry.length += PortionLength(record.portions[i]);
      entry.extra += record.translational[i] ? 1.0 : 0.0;
    }
    WritePortionEnergyTable(spatial_cluster_priority_path, portions_by_type,
                            "geometric cells", "spatialclusterpriority",
                            result.spatial_cluster_priority,
                            result.self_consistent_spatial_cluster_priority);
  }
  if (!result.conductor_consistency.empty())
  {
    // R0-a (idfix DESIGN B.3 (e); S7diag R1, decision 758): the ConductorConsistency share
    // per type (the DomainBoundary table's columns), written once something is excluded.
    std::map<std::string, PortionTypeRow> portions_by_type;
    const auto &record = surface_response_op->GetConductorConsistencyPortionRecord();
    for (std::size_t i = 0; i < record.portions.size(); i++)
    {
      auto &entry = portions_by_type[record.portions[i].topology];
      entry.portions++;
      entry.length += PortionLength(record.portions[i]);
      entry.extra += record.translational[i] ? 1.0 : 0.0;
    }
    WritePortionEnergyTable(conductor_consistency_path, portions_by_type, "geometric cells",
                            "conductorconsistency", result.conductor_consistency,
                            result.self_consistent_conductor_consistency);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureProbes() const
{
  measurement_cache.probe_E_field.clear();
  measurement_cache.probe_En_field.clear();
  measurement_cache.probe_Bt_field.clear();
  measurement_cache.probe_B_field.clear();

#if defined(MFEM_USE_GSLIB)
  if constexpr (HasEGridFunction<solver_t>())
  {
    if (!interp_op.GetProbes().empty() && E)
    {
      measurement_cache.probe_E_field = interp_op.ProbeField(*E);
    }
  }
  if (En && !interp_op.GetProbes().empty())
  {
    measurement_cache.probe_En_field = interp_op.ProbeField(*En);
  }
  if (Bt_inplane && !interp_op.GetProbes().empty())
  {
    measurement_cache.probe_Bt_field = interp_op.ProbeField(*Bt_inplane);
  }
  if constexpr (HasBGridFunction<solver_t>())
  {
    if (!interp_op.GetProbes().empty() && B)
    {
      measurement_cache.probe_B_field = interp_op.ProbeField(*B);
    }
  }
#endif
}

using fmt::format;

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int ex_idx, int step,
                                                  const ComplexVector &e,
                                                  const ComplexVector &b,
                                                  std::complex<double> omega)
  requires(solver_t == ProblemType::DRIVEN)
{
  BlockTimer bt0(Timer::POSTPRO);
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.ex_idx = ex_idx;
  MeasureAllImpl();

  std::complex<double> freq =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * std::numbers::pi);
  post_op_csv.PrintAllCSVData(*this, measurement_cache, freq.real(), step, ex_idx);
  PrintSurfaceResponseCorrection(freq.real(), ex_idx);
  if (ShouldWriteParaviewFields(step))
  {
    Mpi::Print("\n");
    auto ind = 1 + std::distance(output_save_indices.begin(),
                                 std::ranges::lower_bound(output_save_indices, step));
    WriteParaviewFields(omega.real(), ind);
    Mpi::Print(" Wrote fields to disk (Paraview) at step {:d}\n", step + 1);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    Mpi::Print("\n");
    auto ind = 1 + std::distance(output_save_indices.begin(),
                                 std::ranges::lower_bound(output_save_indices, step));
    WriteMFEMGridFunctions(freq.real(), ind);
    Mpi::Print(" Wrote fields to disk (grid function) at step {:d}\n", step + 1);
  }
  ClearSurfaceResponseCorrectedField();
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::ConfigureReducedPostprocessing(const RomOperator &rom_op)
  requires(solver_t == ProblemType::DRIVEN)
{
  reduced_postprocessing_ready = false;
  reduced_postprocessing_checked = false;
  // The reduced magnetic form contains curl(E) only, while a Floquet wave vector adds a
  // frequency-dependent k × E correction to the full magnetic flux density.
  // A surface response correction (the corrected energies, the response tables and the
  // self-consistent field) is measured on the full field only.
  if (fem_op->GetMaterialOp().HasWaveVector() || !fem_op->GetFloquetPortOp().Empty() ||
      fem_op->GetSurfaceCurrentOp().Size() > 0 || !surf_post_op.flux_surfs.empty() ||
      !surf_post_op.eps_surfs.empty() || surf_post_op.farfield.size() > 0 ||
      !interp_op.GetProbes().empty() || !fem_op->GetPortExcitations().IsMultipleSimple() ||
      surface_response_op)
  {
    return;
  }
  for (const auto &[excitation_idx, excitation_spec] : fem_op->GetPortExcitations())
  {
    const auto [is_simple, port_type, port_idx] = excitation_spec.IsSimple();
    if (!is_simple ||
        (port_type != PortType::LumpedPort && port_type != PortType::WavePort))
    {
      return;
    }
  }
  for (const auto &[idx, data] : fem_op->GetWavePortOp())
  {
    if (data.HasVoltageCoords())
    {
      return;
    }
  }

  const auto &basis = rom_op.GetBasis();
  const std::size_t n = basis.size();
  if (n == 0)
  {
    return;
  }
  std::vector<Vector> basis_E, basis_B;
  basis_E.reserve(n);
  basis_B.reserve(n);
  const auto &curl = fem_op->GetCurlMatrix();
  Vector b_true(curl.Height());
  for (const auto &v : basis)
  {
    E->Real().SetFromTrueDofs(v);
    basis_E.emplace_back(E->Real());
    curl.Mult(v, b_true);
    B->Real().SetFromTrueDofs(b_true);
    basis_B.emplace_back(B->Real());
  }

  auto project = [&](const Operator *op,
                     const std::vector<Vector> &local_basis) -> Eigen::MatrixXd
  {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
    if (!op)
    {
      return result;
    }
    Vector z(op->Height());
    for (std::size_t j = 0; j < n; j++)
    {
      op->Mult(local_basis[j], z);
      for (std::size_t i = 0; i < n; i++)
      {
        result(i, j) = linalg::LocalDot(z, local_basis[i]);
      }
    }
    Mpi::GlobalSum(static_cast<int>(n * n), result.data(), fem_op->GetComm());
    return 0.5 * (result + result.transpose()).eval();
  };

  reduced_energy_E = project(dom_post_op.M_elec.get(), basis_E);
  reduced_energy_H = project(dom_post_op.M_mag.get(), basis_B);
  reduced_domain_energy.clear();
  for (const auto &[idx, ops] : dom_post_op.M_i)
  {
    reduced_domain_energy.emplace(idx, std::make_pair(project(ops.first.get(), basis_E),
                                                      project(ops.second.get(), basis_B)));
  }
  reduced_postprocessing_ready = true;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureAndPrintReduced(int ex_idx, int step,
                                                    const ComplexVector &e,
                                                    std::complex<double> omega,
                                                    const Eigen::VectorXcd &y)
  requires(solver_t == ProblemType::DRIVEN)
{
  BlockTimer bt0(Timer::POSTPRO);
  MFEM_VERIFY(reduced_postprocessing_ready && y.size() == reduced_energy_E.rows(),
              "Invalid reduced postprocessing state!");
  SetEGridFunction(e);

  auto quadratic = [&y](const Eigen::MatrixXd &Q)
  { return 0.5 * std::real(y.dot(Q.cast<std::complex<double>>() * y)); };
  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.ex_idx = ex_idx;
  measurement_cache.domain_E_field_energy_all = quadratic(reduced_energy_E);
  measurement_cache.domain_H_field_energy_all =
      quadratic(reduced_energy_H) / (omega.real() * omega.real());
  for (const auto &[idx, forms] : reduced_domain_energy)
  {
    const double energy_E = quadratic(forms.first);
    const double energy_H = quadratic(forms.second) / (omega.real() * omega.real());
    measurement_cache.domain_E_field_energy_i.push_back(
        {idx, energy_E,
         (measurement_cache.domain_E_field_energy_all != 0.0)
             ? energy_E / measurement_cache.domain_E_field_energy_all
             : 0.0});
    measurement_cache.domain_H_field_energy_i.push_back(
        {idx, energy_H,
         (measurement_cache.domain_H_field_energy_all != 0.0)
             ? energy_H / measurement_cache.domain_H_field_energy_all
             : 0.0});
  }
  if (!reduced_postprocessing_checked)
  {
    ComplexVector b(fem_op->GetCurlMatrix().Height());
    b.UseDevice(true);
    fem_op->GetCurlMatrix().Mult(e.Real(), b.Real());
    fem_op->GetCurlMatrix().Mult(e.Imag(), b.Imag());
    b *= -1.0 / (1i * omega);
    SetBGridFunction(b);

    double max_error = 0.0, max_reference = 0.0;
    auto compare = [&](double reduced, double exact)
    {
      max_error = std::max(max_error, std::abs(reduced - exact));
      max_reference = std::max(max_reference, std::abs(exact));
    };
    const double exact_E = dom_post_op.GetElectricFieldEnergy(*E);
    const double exact_H = dom_post_op.GetMagneticFieldEnergy(*B);
    compare(measurement_cache.domain_E_field_energy_all, exact_E);
    compare(measurement_cache.domain_H_field_energy_all, exact_H);
    std::vector<Measurement::DomainData> exact_domain_E, exact_domain_H;
    for (const auto &[idx, forms] : reduced_domain_energy)
    {
      const double energy_E = dom_post_op.GetDomainElectricFieldEnergy(idx, *E);
      const double energy_H = dom_post_op.GetDomainMagneticFieldEnergy(idx, *B);
      compare(measurement_cache.domain_E_field_energy_i[exact_domain_E.size()].energy,
              energy_E);
      compare(measurement_cache.domain_H_field_energy_i[exact_domain_H.size()].energy,
              energy_H);
      exact_domain_E.push_back(
          {idx, energy_E, (exact_E != 0.0) ? energy_E / exact_E : 0.0});
      exact_domain_H.push_back(
          {idx, energy_H, (exact_H != 0.0) ? energy_H / exact_H : 0.0});
    }
    reduced_postprocessing_checked = true;
    if (max_error > 1.0e-9 * std::max(max_reference, 1.0e-300))
    {
      reduced_postprocessing_ready = false;
      measurement_cache.domain_E_field_energy_all = exact_E;
      measurement_cache.domain_H_field_energy_all = exact_H;
      measurement_cache.domain_E_field_energy_i = std::move(exact_domain_E);
      measurement_cache.domain_H_field_energy_i = std::move(exact_domain_H);
      Mpi::Warning("Reduced domain-energy postprocessing disagrees with full evaluation "
                   "(relative error {:.3e}); reverting to full postprocessing.\n",
                   max_error / std::max(max_reference, 1.0e-300));
    }
  }
  const auto port_voltages = fem_op->GetLumpedPortOp().GetVoltages(*E);
  for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
  {
    auto &vi = measurement_cache.lumped_port_vi[idx];
    vi.V = port_voltages.at(idx);
    vi.I_RLC[0] = (std::abs(data.R) > 0.0)
                      ? vi.V / data.GetCharacteristicImpedance(omega.real(),
                                                               LumpedPortData::Branch::R)
                      : 0.0;
    vi.I_RLC[1] = (std::abs(data.L) > 0.0)
                      ? vi.V / data.GetCharacteristicImpedance(omega.real(),
                                                               LumpedPortData::Branch::L)
                      : 0.0;
    vi.I_RLC[2] = (std::abs(data.C) > 0.0)
                      ? vi.V / data.GetCharacteristicImpedance(omega.real(),
                                                               LumpedPortData::Branch::C)
                      : 0.0;
    vi.I =
        std::accumulate(vi.I_RLC.begin(), vi.I_RLC.end(), std::complex<double>{0.0, 0.0});
    vi.S = vi.V / std::sqrt(data.GetExcitationRefResistance());
    if (std::abs(data.L) > 0.0)
    {
      vi.inductor_energy = 0.5 * std::abs(data.L) * std::norm(vi.I_RLC[1]);
      measurement_cache.lumped_port_inductor_energy += vi.inductor_energy;
    }
    if (std::abs(data.C) > 0.0)
    {
      vi.capacitor_energy = 0.5 * std::abs(data.C) * std::norm(vi.V);
      measurement_cache.lumped_port_capacitor_energy += vi.capacitor_energy;
    }
  }
  for (const auto &[idx, data] : fem_op->GetWavePortOp())
  {
    measurement_cache.wave_port_vi[idx].S = data.GetSParameter(*E);
  }
  MeasureSParameter();

  const std::complex<double> freq =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * std::numbers::pi);
  post_op_csv.PrintReducedCSVData(*this, measurement_cache, freq.real(), step, ex_idx);
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const ComplexVector &e,
                                                  const ComplexVector &b,
                                                  std::complex<double> omega,
                                                  double error_abs, double error_bkwd,
                                                  int num_conv)
  requires(solver_t == ProblemType::EIGENMODE)
{
  BlockTimer bt0(Timer::POSTPRO);
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.eigenmode_Q =
      (omega == 0.0) ? mfem::infinity() : 0.5 * std::abs(omega) / std::abs(omega.imag());
  measurement_cache.error_abs = error_abs;
  measurement_cache.error_bkwd = error_bkwd;

  // Mini pretty-print table of eig summaries: always print with header since other
  // measurements may log their results.
  if (Mpi::Root(fem_op->GetComm()))
  {
    Table table;
    int idx_pad = 1 + static_cast<int>(std::log10(num_conv));
    table.col_options = {6, 6};
    table.insert(Column("idx", "m", idx_pad, {}, {}, "") << step + 1);
    table.insert(Column("f_re", "Re{f} (GHz)")
                 << (units.Dimensionalize<Units::ValueType::FREQUENCY>(omega.real())) /
                        (2 * std::numbers::pi));
    table.insert(Column("f_im", "Im{f} (GHz)")
                 << (units.Dimensionalize<Units::ValueType::FREQUENCY>(omega.imag())) /
                        (2 * std::numbers::pi));
    table.insert(Column("q", "Q") << measurement_cache.eigenmode_Q);
    table.insert(Column("err_back", "Error (Bkwd.)") << error_bkwd);
    table.insert(Column("err_abs", "Error (Abs.)") << error_abs);
    table[0].print_as_int = true;
    Mpi::Print("{}", (step == 0) ? table.format_table() : table.format_row(0));
  }
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  PrintSurfaceResponseCorrection(print_idx, 0);
  if (ShouldWriteParaviewFields(step))
  {
    WriteParaviewFields(step, print_idx);
    Mpi::Print(" Wrote mode {:d} to disk (Paraview)\n", print_idx);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    WriteMFEMGridFunctions(step, print_idx);
    Mpi::Print(" Wrote mode {:d} to disk (grid function)\n", print_idx);
  }
  ClearSurfaceResponseCorrectedField();
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &v,
                                                  const Vector &e, int idx)
  requires(solver_t == ProblemType::ELECTROSTATIC)
{
  BlockTimer bt0(Timer::POSTPRO);
  SetVGridFunction(v);
  SetEGridFunction(e);

  measurement_cache = {};
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  if (ShouldWriteParaviewFields(step))
  {
    Mpi::Print("\n");
    WriteParaviewFields(step, idx);
    Mpi::Print(" Wrote fields to disk (Paraview) for source {:d}\n", idx);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    Mpi::Print("\n");
    WriteMFEMGridFunctions(step, idx);
    Mpi::Print(" Wrote fields to disk (grid function) for source {:d}\n", idx);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}
template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &a,
                                                  const Vector &b, int idx)
  requires(solver_t == ProblemType::MAGNETOSTATIC)
{
  BlockTimer bt0(Timer::POSTPRO);
  if (fem_op->IsBroken())
  {
    // The output fields are those of the original DOFs (one side of each two-sided film),
    // the energies are those of the broken solution.
    Vector a0, b0;
    fem_op->GetFirstVersion(fem_op->GetSolutionNDSpace(), a, a0);
    fem_op->GetFirstVersion(fem_op->GetSolutionCurlSpace(), b, b0);
    SetAGridFunction(a0);
    SetBGridFunction(b0);
    GetSolutionLVector(a, a_l);
  }
  else
  {
    SetAGridFunction(a);
    SetBGridFunction(b);
  }

  measurement_cache = {};
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  if (ShouldWriteParaviewFields(step))
  {
    Mpi::Print("\n");
    WriteParaviewFields(step, idx);
    Mpi::Print(" Wrote fields to disk (Paraview) for source {:d}\n", idx);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    Mpi::Print("\n");
    WriteMFEMGridFunctions(step, idx);
    Mpi::Print(" Wrote fields to disk (grid function) for source {:d}\n", idx);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &e,
                                                  const Vector &b, double time,
                                                  double J_coef)
  requires(solver_t == ProblemType::TRANSIENT)
{
  BlockTimer bt0(Timer::POSTPRO);
  // SD-504 (ruling 1574 (3)): with interface "FluxRecovery" the measurement reads the
  // recovered D of this time step, which the driver sets before measuring.
  MFEM_VERIFY(!NeedsRecoveredElectricFlux() || recovered_flux_set,
              "Transient interface \"FluxRecovery\" postprocessing requires the recovered "
              "electric flux of the current time step (SetRecoveredElectricFlux)!");
  recovered_flux_set = false;
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.Jcoeff_excitation = J_coef;
  MeasureAllImpl();

  // Time must be converted before passing into csv due to the shared PrintAllCSVData
  // method.
  time = units.Dimensionalize<Units::ValueType::TIME>(time);
  post_op_csv.PrintAllCSVData(*this, measurement_cache, time, step);
  if (ShouldWriteParaviewFields(step))
  {
    Mpi::Print("\n");
    WriteParaviewFields(time, double(step) / output_delta_post);
    Mpi::Print(" Wrote fields to disk (Paraview) at step {:d}\n", step + 1);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    Mpi::Print("\n");
    WriteMFEMGridFunctions(time, double(step) / output_delta_post);
    Mpi::Print(" Wrote fields to disk (grid function) at step {:d}\n", step + 1);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureFinalize(const ErrorIndicator &indicator)
{
  BlockTimer bt0(Timer::POSTPRO);
  post_op_csv.FinalizeCSVData();
  // Pass nullptr for the indicator if it is empty (no AMR), so the write functions
  // skip the indicator grid function but still save the mesh and rank partition.
  const ErrorIndicator *ind_ptr = (indicator.Local().Size() > 0) ? &indicator : nullptr;
  if (ind_ptr)
  {
    auto indicator_stats = indicator.GetSummaryStatistics(fem_op->GetComm());
    post_op_csv.PrintErrorIndicator(Mpi::Root(fem_op->GetComm()), indicator_stats);
  }
  if (ShouldWriteParaviewFields())
  {
    WriteParaviewFieldsFinal(ind_ptr);
  }
  if (ShouldWriteGridFunctionFields())
  {
    WriteMFEMGridFunctionsFinal(ind_ptr);
  }
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureDomainFieldEnergyOnly(const ComplexVector &e,
                                                            const ComplexVector &b)
  requires(solver_t == ProblemType::DRIVEN)
{
  SetEGridFunction(e);
  SetBGridFunction(b);
  MeasureDomainFieldEnergy();
  Mpi::Barrier(fem_op->GetComm());

  // Return total domain energy for normalizing error indicator.
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const ComplexVector &et,
                                                  const ComplexVector &en,
                                                  std::complex<double> kn, double omega,
                                                  double error_abs, double error_bkwd,
                                                  int num_conv)
  requires(solver_t == ProblemType::BOUNDARYMODE)
{
  BlockTimer bt0(Timer::POSTPRO);

  // Poynting power for impedance postprocessing. After the driver's power normalization
  // |P| ≈ 1, but phase matters for the Z0 = V*conj(V) / (2P) formula.
  const std::complex<double> P = fem_op->ComputePoyntingPower(omega, kn, et, en);

  SetEGridFunction(et);

  // Set the normal (out-of-plane) E component on H1 space.
  En->Real().SetFromTrueDofs(en.Real());
  En->Imag().SetFromTrueDofs(en.Imag());
  En->Real().ExchangeFaceNbrData();
  En->Imag().ExchangeFaceNbrData();

  // Compute Bz = curl_t(Et) / (iω) as a scalar (Hz) on the L2 curl space.
  // B_r = curl(Et_i) / ω,  B_i = -curl(Et_r) / ω  (from 1/(iω) = -i/ω).
  {
    auto &nd_fs = fem_op->GetNDSpace();
    auto &l2_fs = fem_op->GetCurlSpace();
    const auto &Curl = l2_fs.GetDiscreteInterpolator(nd_fs);

    const int nd_sz = nd_fs.GetTrueVSize();
    const int l2_sz = l2_fs.GetTrueVSize();
    ComplexVector b(l2_sz);

    // Get Et true DOFs.
    Vector etr(nd_sz), eti(nd_sz);
    E->Real().GetTrueDofs(etr);
    E->Imag().GetTrueDofs(eti);

    // Apply curl: curl(Et_r) and curl(Et_i).
    Vector curl_etr(l2_sz), curl_eti(l2_sz);
    Curl.Mult(etr, curl_etr);
    Curl.Mult(eti, curl_eti);

    // B = curl(Et) / (iω) = curl(Et) × (-i/ω)
    // B_r = curl(Et_i) / ω,  B_i = -curl(Et_r) / ω
    b.Real() = curl_eti;
    b.Real() *= 1.0 / omega;
    b.Imag() = curl_etr;
    b.Imag() *= -1.0 / omega;

    SetBGridFunction(b);
  }

  // Compute full in-plane B: Bt = -(kn/omega)(z_hat x Et) + (1/(i*omega))(grad_t(Ez) x
  // z_hat) This is the dominant B component for a propagating mode. Project onto ND space.
  {
    ModeInPlaneBCoefficient bt_real_coeff(E->Real(), E->Imag(), En->Real(), En->Imag(),
                                          kn.real(), kn.imag(), omega, true);
    ModeInPlaneBCoefficient bt_imag_coeff(E->Real(), E->Imag(), En->Real(), En->Imag(),
                                          kn.real(), kn.imag(), omega, false);
    Bt_inplane->Real().ProjectCoefficient(bt_real_coeff);
    Bt_inplane->Imag().ProjectCoefficient(bt_imag_coeff);
    Bt_inplane->Real().ExchangeFaceNbrData();
    Bt_inplane->Imag().ExchangeFaceNbrData();
  }

  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.error_abs = error_abs;
  measurement_cache.error_bkwd = error_bkwd;

  // Store nondimensional mode data (dimensionalized in Measurement::Dimensionalize).
  std::complex<double> n_eff = kn / omega;
  measurement_cache.mode_data.kn = kn;
  measurement_cache.mode_data.n_eff = n_eff;

  struct LineIntegralKey
  {
    bool has_coords = false;
    std::vector<double> path_data;
    std::vector<int> marker_data;

    auto operator<=>(const LineIntegralKey &) const = default;
  };
  auto MakeLineIntegralKey = [](const std::vector<mfem::Vector> &path, bool has_coords,
                                const mfem::Array<int> &marker)
  {
    LineIntegralKey key;
    key.has_coords = has_coords;
    if (has_coords)
    {
      for (const auto &pt : path)
      {
        for (int d = 0; d < pt.Size(); d++)
        {
          key.path_data.push_back(pt(d));
        }
      }
    }
    else
    {
      key.marker_data.reserve(marker.Size());
      for (int i = 0; i < marker.Size(); i++)
      {
        key.marker_data.push_back(marker[i]);
      }
    }
    return key;
  };
  std::map<LineIntegralKey, std::complex<double>> voltage_cache, current_cache;

  // Helper: compute voltage V = ∫ E · dl via coordinate path or boundary attributes.
  auto ComputeVoltage = [this, &MakeLineIntegralKey, &voltage_cache](
                            const std::vector<mfem::Vector> &path, bool has_coords,
                            mfem::Array<int> marker, const fem::LineIntegrator *integrator,
                            const GridFunction &field) -> std::complex<double>
  {
    auto key = MakeLineIntegralKey(path, has_coords, marker);
    if (auto it = voltage_cache.find(key); it != voltage_cache.end())
    {
      return it->second;
    }

    std::complex<double> V(0.0, 0.0);
    if (has_coords)
    {
      MFEM_VERIFY(integrator, "Missing line integrator for coordinate voltage path!");
      V = integrator->Integrate(field);
    }
    else if (marker.Size() > 0)
    {
      auto &nd_fespace = fem_op->GetNDSpace();
      const int nd_size = nd_fespace.GetTrueVSize();
      mfem::ConstantCoefficient one(1.0);
      mfem::ParLinearForm v_lf(&nd_fespace.Get());
      v_lf.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(one), marker);
      v_lf.UseFastAssembly(false);
      v_lf.Assemble();
      Vector v_lf_tdof(nd_size);
      v_lf.ParallelAssemble(v_lf_tdof);
      Vector fr(nd_size), fi(nd_size);
      field.Real().GetTrueDofs(fr);
      field.Imag().GetTrueDofs(fi);
      V.real(linalg::Dot(nd_fespace.GetComm(), v_lf_tdof, fr));
      V.imag(linalg::Dot(nd_fespace.GetComm(), v_lf_tdof, fi));
    }
    voltage_cache.emplace(std::move(key), V);
    return V;
  };

  // Helper: compute current I = ∮ H · t dl via coordinate path or boundary attributes.
  auto ComputeCurrent = [this, &MakeLineIntegralKey, &current_cache](
                            const std::vector<mfem::Vector> &path, bool has_path,
                            mfem::Array<int> marker, const fem::LineIntegrator *integrator,
                            const GridFunction &field) -> std::complex<double>
  {
    auto key = MakeLineIntegralKey(path, has_path, marker);
    if (auto it = current_cache.find(key); it != current_cache.end())
    {
      return it->second;
    }

    std::complex<double> I(0.0, 0.0);
    if (has_path)
    {
      MFEM_VERIFY(integrator, "Missing line integrator for coordinate current path!");
      I = integrator->Integrate(field);
    }
    else if (marker.Size() > 0)
    {
      auto &nd_fespace = fem_op->GetNDSpace();
      const int nd_size = nd_fespace.GetTrueVSize();
      mfem::ConstantCoefficient one(1.0);
      mfem::ParLinearForm i_lf(&nd_fespace.Get());
      i_lf.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(one), marker);
      i_lf.UseFastAssembly(false);
      i_lf.Assemble();
      Vector i_lf_tdof(nd_size);
      i_lf.ParallelAssemble(i_lf_tdof);
      Vector fr(nd_size), fi(nd_size);
      field.Real().GetTrueDofs(fr);
      field.Imag().GetTrueDofs(fi);
      I.real(linalg::Dot(nd_fespace.GetComm(), i_lf_tdof, fr));
      I.imag(linalg::Dot(nd_fespace.GetComm(), i_lf_tdof, fi));
    }
    current_cache.emplace(std::move(key), I);
    return I;
  };

  auto SameMarker = [](const mfem::Array<int> &a, const mfem::Array<int> &b)
  {
    if (a.Size() != b.Size())
    {
      return false;
    }
    for (int i = 0; i < a.Size(); i++)
    {
      if (a[i] != b[i])
      {
        return false;
      }
    }
    return true;
  };
  auto SamePath = [](const std::vector<mfem::Vector> &a, const std::vector<mfem::Vector> &b)
  {
    if (a.size() != b.size())
    {
      return false;
    }
    for (std::size_t i = 0; i < a.size(); i++)
    {
      if (a[i].Size() != b[i].Size())
      {
        return false;
      }
      for (int d = 0; d < a[i].Size(); d++)
      {
        if (a[i](d) != b[i](d))
        {
          return false;
        }
      }
    }
    return true;
  };
  auto SameVoltageConfig =
      [&](const ImpedancePostproConfig &imp, const VoltagePostproConfig &vol)
  {
    return imp.has_voltage_coordinates == vol.has_coordinates &&
           (imp.has_voltage_coordinates
                ? SamePath(imp.voltage_path, vol.voltage_path)
                : SameMarker(imp.voltage_marker, vol.voltage_marker));
  };
  std::set<int> filled_voltage_postpro;

  // Compute impedance for each configured entry. P is the Poynting power of the
  // power-normalized mode (|P| ≈ 1), computed by the caller on the mode eigensolver.
  for (const auto &[idx, cfg] : impedance_postpro)
  {
    auto V = ComputeVoltage(cfg.voltage_path, cfg.has_voltage_coordinates,
                            cfg.voltage_marker, cfg.voltage_integrator.get(), *E);
    auto &result = measurement_cache.mode_data.impedance[idx];
    if (auto vol_it = voltage_postpro.find(idx);
        vol_it != voltage_postpro.end() && SameVoltageConfig(cfg, vol_it->second))
    {
      measurement_cache.mode_data.voltage[idx] = {V};
      filled_voltage_postpro.insert(idx);
    }

    // Power-voltage impedance Z_PV = |V|² / 2P* with the complex (time-averaged) power P of
    // the mode: for a TEM mode, P = V I* / 2 and Z_PV = V / I is the complex characteristic
    // impedance of the line, real for a lossless mode.
    if (std::abs(P) > 1e-30)
    {
      result.Z0 = (V * std::conj(V)) / (2.0 * std::conj(P));
      result.has_impedance = true;
    }

    // Voltage-current impedance Z_VI = V / I, oriented with the mode power (Re{V I* / P} >
    // 0, as P = V I* / 2 for a TEM mode), so that it does not depend on the directions of
    // the voltage path and of the current loop. For a TEM mode, Z_VI = Z_PV.
    if (cfg.has_current)
    {
      auto I = ComputeCurrent(cfg.current_path, cfg.has_current_path, cfg.current_marker,
                              cfg.current_integrator.get(), *Bt_inplane);
      if (std::abs(I) > 1e-30)
      {
        result.Z_VI = V / I;
        if ((V * std::conj(I) * std::conj(P)).real() < 0.0)
        {
          result.Z_VI = -result.Z_VI;
        }
        result.has_vi_impedance = true;
      }
    }
  }

  // Compute voltage for each configured voltage-only entry.
  for (const auto &[idx, cfg] : voltage_postpro)
  {
    if (filled_voltage_postpro.contains(idx))
    {
      continue;
    }
    auto V = ComputeVoltage(cfg.voltage_path, cfg.has_coordinates, cfg.voltage_marker,
                            cfg.voltage_integrator.get(), *E);
    measurement_cache.mode_data.voltage[idx] = {V};
  }

  // Compute electric field energy via the domain post operator.
  MeasureAllImpl();

  // Pretty-print mode summary using Table (matches eigenmode solver style).
  if (Mpi::Root(fem_op->GetComm()))
  {
    Table table;
    int idx_pad = 1 + static_cast<int>(std::log10(std::max(num_conv, 1)));
    table.col_options = {6, 6};
    table.insert(Column("m", "m", idx_pad, {}, {}, "") << (step + 1));
    const double kc = 1.0 / units.Dimensionalize<Units::ValueType::LENGTH>(1.0);
    table.insert(Column("kn_re", "Re{kn} (1/m)") << kn.real() * kc);
    table.insert(Column("kn_im", "Im{kn} (1/m)") << kn.imag() * kc);
    table.insert(Column("neff_re", "Re{n_eff}") << n_eff.real());
    table.insert(Column("neff_im", "Im{n_eff}") << n_eff.imag());
    table.insert(Column("err_back", "Error (Bkwd.)") << error_bkwd);
    table.insert(Column("err_abs", "Error (Abs.)") << error_abs);
    table[0].print_as_int = true;
    Mpi::Print("{}", (step == 0) ? table.format_table() : table.format_row(0));
  }

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  PrintSurfaceResponseCorrection(print_idx, 0);
  if (ShouldWriteParaviewFields(step))
  {
    WriteParaviewFields(step, print_idx);
    Mpi::Print(" Wrote mode {:d} to disk (Paraview)\n", print_idx);
  }
  if (ShouldWriteGridFunctionFields(step))
  {
    WriteMFEMGridFunctions(step, print_idx);
    Mpi::Print(" Wrote mode {:d} to disk (grid function)\n", print_idx);
  }
  return measurement_cache.domain_E_field_energy_all;
}

// Explicit template instantiation.

template class PostOperator<ProblemType::DRIVEN>;
template class PostOperator<ProblemType::EIGENMODE>;
template class PostOperator<ProblemType::ELECTROSTATIC>;
template class PostOperator<ProblemType::MAGNETOSTATIC>;
template class PostOperator<ProblemType::TRANSIENT>;
template class PostOperator<ProblemType::BOUNDARYMODE>;
}  // namespace palace
