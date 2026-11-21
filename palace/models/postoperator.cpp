// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "postoperator.hpp"

#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <tuple>
#include <vector>
#include "drivers/boundarymodesolver.hpp"
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "fem/face_sampling_plan.hpp"
#include "fem/interpolator.hpp"
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp"
#include "linalg/rap.hpp"
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

// The fraction of a total energy in a part of the domain, zero for a zero total.
double ParticipationRatio(double energy, double total)
{
  return (total != 0.0) ? energy / total : 0.0;
}

}  // namespace

template <ProblemType solver_t>
PostOperator<solver_t>::PostOperator(const config::ProblemData &problem,
                                     const config::SolverData &solver,
                                     const config::DomainData &domains,
                                     const config::BoundaryData &boundaries,
                                     const Units &units_, fem_op_t<solver_t> &fem_op_)
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
            // Share the unit masses B_m of SpaceOperator, rather than assemble them again.
            const auto &mat_op = fem_op_.GetMaterialOp();
            std::vector<const Operator *> unit_mass(
                mat_op.NumFrequencyDependentPermittivityMaterials());
            for (std::size_t m = 0; m < unit_mass.size(); m++)
            {
              unit_mass[m] = fem_op_.GetFrequencyDependentPermittivityLocalMassMatrix(m);
            }
            return DomainPostOperator(domains.postpro, mat_op, fem_op_.GetNDSpace(),
                                      fem_op_.GetCurlSpace(), unit_mass);
          }
        }())),
    surf_post_op(boundaries.postpro, solver_t, fem_op->GetMaterialOp(),
                 fem_op->GetH1Space(), fem_op->GetNDSpace(),
                 boundaries.GetSheetAttributes()),
    interp_op(domains.postpro.probe, units_, fem_op->GetNDSpace(),
              solver.substructuring ? solver.substructuring->region_attributes
                                    : std::vector<int>{})
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
    // Add scaling factor to output dimensional wave port electric fields: the field of the
    // incident wave of 1 W time-averaged power, sqrt(2) times the unit-overlap mode.
    const double scaling = electromagnetics::UnitPowerAmplitude(true) *
                           units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
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

  // The electric energies a frequency-dependent material without a stored electric energy
  // would enter are omitted, with the outputs normalized by them. Decided once for the run,
  // the same on every process, for the materials on the mesh only.
  {
    const auto &mat_op = fem_op->GetMaterialOp();
    std::vector<std::size_t> materials;
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (mat_op.HasFrequencyDependentPermittivityModel(m) &&
          !mat_op.HasStoredElectricEnergy(m))
      {
        materials.push_back(m + 1);
      }
    }
    if (!materials.empty())
    {
      std::vector<int> regions;
      for (const auto &[idx, data] : dom_post_op.M_i)
      {
        if (!dom_post_op.HasStoredDomainElectricFieldEnergy(idx))
        {
          regions.push_back(idx);
        }
      }
      Mpi::Warning(
          fem_op->GetComm(),
          "Frequency-dependent material(s) {} have no stored electric energy (their "
          "PoleResidue terms specify none)!\nOmitted outputs: the total electric "
          "energy E_elec and the electric participations p_elec, the electric energy "
          "density U_e, the "
          "outputs normalized by the electric energy (surface-Q.csv, port-Q.csv, "
          "port-EPR.csv){}. The field solution, S-parameters, port voltages, "
          "currents, and powers, fluxes, eigenvalues and their Q, and magnetic "
          "energies are reported.\n",
          materials,
          regions.empty()
              ? std::string()
              : fmt::format(", and E_elec of the postprocessing domains {}", regions));
    }
  }

  // Initialize CSV files for measurements.
  post_op_csv.InitializeCSVDataCollection(*this);
}

template <ProblemType solver_t>
PostOperator<solver_t>::PostOperator(const IoData &iodata, fem_op_t<solver_t> &fem_op_)
  : PostOperator(iodata.problem, iodata.solver, iodata.domains, iodata.boundaries,
                 iodata.units, fem_op_)
{
}

template <ProblemType solver_t>
fs::path PostOperator<solver_t>::ExcitationSubFolder(int ex_idx) const
  requires(solver_t == ProblemType::DRIVEN)
{
  fs::path sub_folder_name = "";
  auto nr_excitations = fem_op->GetPortExcitations().Size();
  if ((nr_excitations > 1) && (ex_idx > 0))
  {
    int spacing = 1 + int(std::log10(nr_excitations));
    sub_folder_name = fmt::format(FMT_STRING("excitation_{:0>{}}"), ex_idx, spacing);
  }
  return sub_folder_name;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::InitializeParaviewDataCollection(int ex_idx, bool restart)
  requires(solver_t == ProblemType::DRIVEN)
{
  InitializeParaviewDataCollection(ExcitationSubFolder(ex_idx), restart);
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

  // Quadratic field outputs (energy densities, Poynting vectors) are physical time averages
  // for the complex peak phasors of frequency domain solvers and instantaneous values for
  // real-valued fields. The coefficient and libCEED kernels evaluate the raw same-phase
  // sums (1/2 DᴴE, E x H⋆, ...) and the owner applies the time-average weight through their
  // output scaling, consistently across the CPU coefficient, libCEED domain, boundary, and
  // fused boundary bundle paths. Linear field outputs (E, B, D ⋅ n, n x H) are unaffected.
  constexpr double quadratic_weight =
      electromagnetics::TimeAverageWeight(HasComplexGridFunction<solver_t>());

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
            bdr_sampling_plan, bdr_trace_cache, *E, *B, electric_scaling, magnetic_scaling,
            quadratic_weight);
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
    // Boundary mode longitudinal power density Sn = 1/2 Re{Et · (ẑ × Ht⋆)} for
    // complex peak phasors; also convert the dimensional B input to H.
    MakeFieldEvaluator(
        PointFieldEvaluator::Kind::MODE_SN, E->ParFESpace(), Bt_inplane->ParFESpace(),
        quadratic_weight * units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
            units.Dimensionalize<Units::ValueType::FIELD_B>(1.0),
        Sn_eval, Sn_gf);
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
    // U_e = 1/2 Dᴴ E = 1/2 ε_0 Eᴴ E (time-averaged, 1/4 ε_0 Eᴴ E, for complex peak
    // phasors). Omitted without the stored electric energy of every frequency-dependent
    // material, the same on every process.
    if (fem_op->GetMaterialOp().HasStoredElectricEnergy())
    {
      U_e = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(
          *E, fem_op->GetMaterialOp(), quadratic_weight * scaling);
      if (use_ceed_domain_fields)
      {
        MakeFieldEvaluator(PointFieldEvaluator::Kind::ENERGY_E, E->ParFESpace(), nullptr,
                           quadratic_weight * scaling, U_e_eval, U_e_gf);
      }
      if (use_ceed_boundary_fields)
      {
        MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::ENERGY_E, *E->ParFESpace(),
                              quadratic_weight * scaling, Ue_bdr_eval);
      }
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
      if (HasComplexGridFunction<solver_t>() &&
          (fem_op->GetMaterialOp().HasLossTangent() ||
           fem_op->GetMaterialOp().HasFrequencyDependentPermittivity()) &&
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
    // U_m = 1/2 Hᴴ B = 1/2 μ⁻¹ Bᴴ B (time-averaged, 1/4 μ⁻¹ Bᴴ B, for complex peak
    // phasors).
    U_m = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(
        *B, fem_op->GetMaterialOp(), quadratic_weight * scaling);
    if (use_ceed_domain_fields)
    {
      MakeFieldEvaluator(PointFieldEvaluator::Kind::ENERGY_M, nullptr, B->ParFESpace(),
                         quadratic_weight * scaling, U_m_eval, U_m_gf);
    }
    if (use_ceed_boundary_fields)
    {
      MakeBdrCoeffEvaluator(PointFieldEvaluator::Kind::ENERGY_M, *B->ParFESpace(),
                            quadratic_weight * scaling, Um_bdr_eval);
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
      // 3D Poynting Vector: S = Re{E x H⋆} = Re{E x μ⁻¹B⋆} (time-averaged, 1/2 Re{E x H⋆},
      // for complex peak phasors).
      // E and B will be dimensionalized when the coefficient is evaluated, so the scaling
      // only needs to account for the remaining μ⁻¹ = H / B. As mentioned above,
      // ProjectCoefficient(*S.get()) or paraview->Save() should always be called after
      // E and B have been dimensionalized.
      const double scaling = quadratic_weight *
                             units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
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
    // For boundary mode, Sn = 1/2 Re{Et · (ẑ × Ht*)} uses the in-plane B field and is
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
    const fs::path &sub_folder_name, bool restart)
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
  paraview->UseRestartMode(restart);
  paraview_bdr->UseRestartMode(restart);

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
  // repeating VTU setup across mode saves. The fields are evaluated at the real part of
  // the measured frequency (of an eigenmode), for the surface charge and electric energy
  // density of a frequency-dependent material.
  const auto *freq = &measurement_cache.freq;
  auto RegisterDomainEvalField =
      [&](const std::string &name, const std::unique_ptr<PointFieldEvaluator> &eval,
          const GridFunction *E_field, const GridFunction *B_field)
  {
    const auto *eval_ptr = eval.get();
    const auto *E_ptr = E_field;
    const auto *B_ptr = B_field;
    paraview->RegisterDomainPointEvaluator(
        name, [eval_ptr, E_ptr, B_ptr, freq](Vector &buffer)
        { eval_ptr->EvalBuffer(E_ptr, B_ptr, freq->real(), buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize(), true);
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
        name, [eval_ptr, field_ptr, freq](Vector &buffer)
        { eval_ptr->EvalBuffer(*field_ptr, freq->real(), buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize(), true);
  };
  auto RegisterBdrEvalField = [&](const std::string &name,
                                  const std::unique_ptr<PointFieldEvaluator> &eval,
                                  const auto &field)
  {
    const auto *eval_ptr = eval.get();
    const auto *field_ptr = &field;
    paraview_bdr->RegisterBoundaryPointEvaluator(
        name, [eval_ptr, field_ptr, freq](Vector &buffer)
        { eval_ptr->EvalBuffer(*field_ptr, freq->real(), buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize());
  };
  auto RegisterBdrPoyntingField = [&](const std::string &name)
  {
    const auto *eval_ptr = S_bdr_eval.get();
    const auto *E_ptr = E.get();
    const auto *B_ptr = B.get();
    paraview_bdr->RegisterBoundaryPointEvaluator(
        name, [eval_ptr, E_ptr, B_ptr, freq](Vector &buffer)
        { eval_ptr->EvalBuffer(E_ptr, B_ptr, freq->real(), buffer); },
        eval_ptr->BufferBases(), eval_ptr->BufferNumComp(), eval_ptr->BufferSize());
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
  // energy 1/2 Hᴴ B, and Poynting vector S = E x H⋆, as instantaneous values for real
  // fields and time averages (1/4 Dᴴ E, 1/4 Hᴴ B, 1/2 Re{E x H⋆}) for complex peak phasors
  // (see SetupFieldCoefficients).
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
            [eval_ptr, eval_imag_ptr, E_ptr, imaginary_phase, freq](Vector &buffer)
            {
              eval_ptr->EvalBuffer(imaginary_phase ? E_ptr->Imag() : E_ptr->Real(),
                                   freq->real(), buffer);
              if (eval_imag_ptr)
              {
                Vector cross(buffer.Size());
                cross.UseDevice(buffer.UseDevice());
                eval_imag_ptr->EvalBuffer(imaginary_phase ? E_ptr->Real() : E_ptr->Imag(),
                                          freq->real(), cross);
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
        bdr_bundle_fields, bundle->BufferBases(), [bundle, freq]() -> const Vector &
        { return bundle->GetRegisteredPacked(freq->real()); },
        bundle->RegisteredPackedSize());
  }

  // Add wave port boundary mode postprocessing when available.
  for (const auto &[idx, data] : port_E0)
  {
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_real", idx), data.E0r.get());
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_imag", idx), data.E0i.get());
  }
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
void PostOperator<solver_t>::WriteMFEMGridFunctions(double time, int step,
                                                    const fs::path &sub_folder_name)
{
  BlockTimer bt(Timer::POSTPRO_GRIDFUNCTION);
  EnsureFieldCoefficientsSetup();

  // Create output directory if it doesn't exist (node-local filesystem aware).
  fs::path output_dir = gridfunction_output_dir;
  if (!sub_folder_name.empty())
  {
    output_dir /= sub_folder_name;
  }
  EnsureDirectory(output_dir, fem_op->GetComm());

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization. Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  DimensionalizeGridFunctions(units, E, B, V, A);
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

  std::string error;
  auto write_grid_function = [&](const auto &gridfunc, const std::string &name)
  {
    auto path =
        output_dir / fmt::format("{}_{:0{}d}.gf.{:0{}d}", name, step, pad_digits_default,
                                 local_rank, pad_digits_default);
    WriteRankOutputFile(path, [&](std::ostream &os) { gridfunc.Save(os); }, error);
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

  // The fields are evaluated at the real part of the measured frequency (of an eigenmode).
  const double omega = measurement_cache.freq.real();
  if (U_e)
  {
    if (U_e_eval)
    {
      // Device-evaluated (libCEED) energy density into the interpolatory L2 output
      // space, reusing the ParaView evaluator instead of a host coefficient projection.
      U_e_eval->Eval(E.get(), nullptr, omega, *U_e_gf);
      write_grid_function(*U_e_gf, "U_e");
    }
    else
    {
      const auto &mat_op = fem_op->GetMaterialOp();
      if (mat_op.HasFrequencyDependentPermittivity())
      {
        static_cast<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> &>(*U_e)
            .SetStoredEnergyPermittivity(mat_op.EvaluateStoredEnergyPermittivity(omega));
      }
      gridfunc_scalar = 0.0;
      gridfunc_scalar.ProjectCoefficient(*U_e.get());
      write_grid_function(gridfunc_scalar, "U_e");
    }
  }

  if (U_m)
  {
    if (U_m_eval)
    {
      U_m_eval->Eval(nullptr, B.get(), omega, *U_m_gf);
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
      S_eval->Eval(E.get(), B.get(), omega, *S_gf);
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
    Sn_eval->Eval(E.get(), Bt_inplane.get(), omega, *Sn_gf);
    write_grid_function(*Sn_gf, "Sn");
  }

  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0,
                     {E.get(), B.get(), V.get(), A.get(), En.get(), Bt_inplane.get()});
  NondimensionalizeGridFunctions(units, E, B, V, A);
  if (En)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(*En);
  }
  if (Bt_inplane)
  {
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_B>(*Bt_inplane);
  }
  AbortOnOutputFileErrors(fem_op->GetComm(), error);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteMFEMGridFunctionsFinal(const ErrorIndicator *indicator)
{
  BlockTimer bt(Timer::POSTPRO_GRIDFUNCTION);
  EnsureFieldCoefficientsSetup();

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Write the mesh partitioning and (optionally) error indicators at the final step.
  // Write the mesh coordinates in the same units as originally input. These are written
  // once, outside of any excitation subfolders, since all the excitations of a driven
  // simulation share the mesh, and the error indicator combines them.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);

  // Create output directory if it doesn't exist (node-local filesystem aware).
  EnsureDirectory(fs::path(gridfunction_output_dir), fem_op->GetComm());

  // Create piecewise constant finite element space for rank and error indicator.
  mfem::L2_FECollection pwconst_fec(0, mesh.Dimension());
  mfem::FiniteElementSpace pwconst_fespace(&mesh, &pwconst_fec);

  const int local_rank = mesh.GetMyRank();

  std::string error;
  auto write_grid_function = [&](const auto &gridfunc, const std::string &name)
  {
    auto path = fs::path(gridfunction_output_dir) /
                fmt::format("{}.gf.{:0{}d}", name, local_rank, pad_digits_default);
    WriteRankOutputFile(path, [&](std::ostream &os) { gridfunc.Save(os); }, error);
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

  // Save ParMesh files; necessary to visualize grid functions. Written as
  // mfem::ParMesh::Save writes them, which does not check the write.
  WriteRankOutputFile(
      fs::path(gridfunction_output_dir) /
          fmt::format("mesh.{:0{}d}", local_rank, pad_digits_default),
      [&](std::ostream &os)
      {
        os.precision(16);
        mesh.Print(os);
      },
      error);

  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  AbortOnOutputFileErrors(fem_op->GetComm(), error);
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
    // The stored energy of a frequency-dependent material is evaluated at the real part of
    // the frequency (of an eigenmode).
    // Without a stored electric energy of a frequency-dependent material (the same on
    // every process), the energies with it, and the participations in the total, are
    // left empty. A given total (substructuring, whose materials are constant) is both.
    const double omega = measurement_cache.freq.real();
    const auto [energy, material] =
        total_field_energy ? DomainPostOperator::ElectricFieldEnergies{*total_field_energy,
                                                                       *total_field_energy}
                           : dom_post_op.GetElectricFieldEnergies(field, omega);
    measurement_cache.domain_E_field_energy_all = energy;
    measurement_cache.domain_E_field_energy_material = material;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      Measurement::DomainData &d = measurement_cache.domain_E_field_energy_i.emplace_back(
          Measurement::DomainData{idx, std::nullopt, std::nullopt});
      if (dom_post_op.HasStoredDomainElectricFieldEnergy(idx))
      {
        d.energy = dom_post_op.GetDomainElectricFieldEnergy(idx, field, omega);
        if (energy)
        {
          d.participation_ratio = ParticipationRatio(*d.energy, *energy);
        }
      }
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
    auto energy = total_field_energy ? *total_field_energy
                  : use_l ? dom_post_op.GetMagneticFieldEnergy(a_l, field.GetComm())
                          : dom_post_op.GetMagneticFieldEnergy(field);
    measurement_cache.domain_H_field_energy_all = energy;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      auto energy_i =
          use_l ? dom_post_op.GetDomainMagneticFieldEnergy(idx, a_l, field.GetComm())
                : dom_post_op.GetDomainMagneticFieldEnergy(idx, field);
      measurement_cache.domain_H_field_energy_i.emplace_back(
          Measurement::DomainData{idx, energy_i, ParticipationRatio(energy_i, energy)});
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
  if (dom_post_op.partial)
  {
    return;
  }
  const auto &E_all = measurement_cache.domain_E_field_energy_all;
  const auto domain_H = units.Dimensionalize<Units::ValueType::ENERGY>(
      measurement_cache.domain_H_field_energy_all);
  if constexpr (HasEGridFunction<solver_t>() && !HasBGridFunction<solver_t>())
  {
    if (E_all)
    {
      Mpi::Print(" Field energy E = {:.3e} J\n",
                 units.Dimensionalize<Units::ValueType::ENERGY>(*E_all));
    }
  }
  else if constexpr (!HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
  {
    Mpi::Print(" Field energy H = {:.3e} J\n", domain_H);
  }
  else if constexpr (solver_t != ProblemType::EIGENMODE &&
                     solver_t != ProblemType::BOUNDARYMODE)
  {
    if (E_all)
    {
      const auto domain_E = units.Dimensionalize<Units::ValueType::ENERGY>(*E_all);
      Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", domain_E,
                 domain_H, domain_E + domain_H);
    }
    else
    {
      Mpi::Print(" Field energy H = {:.3e} J (no stored electric energy)\n", domain_H);
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPorts() const
{
  if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN ||
                solver_t == ProblemType::TRANSIENT)
  {
    // The voltages given by SetLumpedPortVoltages replace those of the field.
    auto voltages = fem_op->GetLumpedPortOp().GetVoltages(*E);
    for (const auto &[idx, given] : lumped_port_voltages)
    {
      voltages[idx] = given;
    }
    MeasureLumpedPorts(voltages);
  }
  else
  {
    measurement_cache.lumped_port_vi.clear();
    measurement_cache.lumped_port_inductor_energy = 0.0;
    measurement_cache.lumped_port_capacitor_energy = 0.0;
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPorts(
    const std::map<int, std::complex<double>> &port_voltages) const
{
  measurement_cache.lumped_port_vi.clear();
  measurement_cache.lumped_port_inductor_energy = 0.0;
  measurement_cache.lumped_port_capacitor_energy = 0.0;

  if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN)
  {
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      auto &vi = measurement_cache.lumped_port_vi[idx];
      vi.V = port_voltages.at(idx);

      // Compute current from the port impedance, separate contributions for R, L, C
      // branches.
      // Get value and make real: Matches current behaviour (even for eigensolver!).
      const double freq_re = measurement_cache.freq.real();
      MFEM_VERIFY(
          freq_re > 0.0,
          "Frequency domain lumped port postprocessing requires nonzero frequency!");
      using Branch = LumpedPortData::Branch;
      const std::array<std::pair<double, Branch>, 3> branches = {
          {{data.R, Branch::R}, {data.L, Branch::L}, {data.C, Branch::C}}};
      for (std::size_t k = 0; k < branches.size(); k++)
      {
        const auto [value, branch] = branches[k];
        vi.I_RLC[k] = (std::abs(value) > 0.0)
                          ? vi.V / data.GetCharacteristicImpedance(freq_re, branch)
                          : 0.0;
      }
      vi.I =
          std::accumulate(vi.I_RLC.begin(), vi.I_RLC.end(), std::complex<double>{0.0, 0.0});
      vi.S = vi.V / data.GetReferenceVoltage();

      // Add contribution due to all inductive lumped boundaries in the model, as the
      // time-averaged stored energy for the peak phasor branch current:
      //                      E_ind = ∑_j 1/2 L_j <i_j(t)²> = ∑_j 1/4 L_j |I_mj|².
      if (std::abs(data.L) > 0.0)
      {
        vi.inductor_energy = 0.5 * electromagnetics::TimeAverageWeight(true) *
                             std::abs(data.L) * std::norm(vi.I_RLC[1]);
        measurement_cache.lumped_port_inductor_energy += vi.inductor_energy;
      }

      // Add contribution due to all capacitive lumped boundaries in the model, as the
      // time-averaged stored energy for the peak phasor port voltage:
      //                      E_cap = ∑_j 1/2 C_j <v_j(t)²> = ∑_j 1/4 C_j |V_mj|².
      if (std::abs(data.C) > 0.0)
      {
        vi.capacitor_energy = 0.5 * electromagnetics::TimeAverageWeight(true) *
                              std::abs(data.C) * std::norm(vi.V);
        measurement_cache.lumped_port_capacitor_energy += vi.capacitor_energy;
      }
    }
  }
  else if constexpr (solver_t == ProblemType::TRANSIENT)
  {
    // Compute current from P = V I^* since there is no frequency & characteristic
    // impedance of the lumped element.
    const auto port_powers = fem_op->GetLumpedPortOp().GetPowers(*E, *B);
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      auto &vi = measurement_cache.lumped_port_vi[idx];
      vi.P = port_powers.at(idx);
      vi.V = port_voltages.at(idx);
      vi.I = (std::abs(vi.V) > 0.0) ? std::conj(vi.P / vi.V) : 0.0;
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPortsEig() const
{
  // Depends on MeasureLumpedPorts.
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    // The coupling rates and participations are normalized by the time-averaged
    // electric-side energy U (Measurement::ElectricEnergyNormalization), which is
    // E_elec + E_cap for constant materials, and are left empty without one. The total
    // stored energy uses the weak-loss equipartition 2 U (E_elec + E_cap = E_mag + E_ind),
    // which leaves out distributed sheet storage.
    const auto energy_electric_all = measurement_cache.ElectricEnergyNormalization();
    if (!energy_electric_all)
    {
      return;
    }
    const double energy_total = 2.0 * *energy_electric_all;
    auto freq_re = measurement_cache.freq.real();
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      // Get previously computed data: should never fail as defined by MeasureLumpedPorts.
      auto &vi = measurement_cache.lumped_port_vi.at(idx);

      // Resistive Lumped Ports:
      // Compute participation ratio of external ports (given as any port boundary with
      // nonzero resistance). Currently no reactance of the ports is supported. The
      // magnitude of the κ of the port is the decay rate due to the time-averaged power
      // P_mj dissipated in the resistive sheet, relative to the total stored energy E_m of
      // the mode, and its sign is the sign of Re{I_mj}, which is not meaningful when
      // Re{I_mj} ≈ 0, for example when V_mj ≈ 0:
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
        const double kappa = std::copysign(resistor_power / energy_total, I_mj.real());
        vi.mode_port_kappa = kappa;
        vi.quality_factor = (kappa == 0.0) ? mfem::infinity() : freq_re / std::abs(kappa);
      }

      // Inductive Lumped Ports:
      // Compute energy-participation ratio of junction given by index idx for the field
      // mode. We first get the port line voltage, and use lumped port circuit impedance to
      // get peak current through the inductor: I_mj = V_mj / Z_mj,  Z_mj = i ω_m L_j. The
      // signed EPR for a lumped inductive element is the ratio of the time-averaged
      // inductor energy to the time-averaged electric-side energy U of the mode
      // (Measurement::ElectricEnergyNormalization), E_elec + E_cap for constant materials:
      //                            p_mj = 1/4 L_j |I_mj|² / U.
      // An element with no assigned inductance will be treated as having zero admittance
      // and thus zero current.
      if (std::abs(data.L) > 0.0)
      {
        std::complex<double> I_mj = vi.I_RLC[1];
        vi.inductive_energy_participation =
            std::copysign(vi.inductor_energy / *energy_electric_all, I_mj.real());
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
      const auto given = wave_port_overlaps.find(idx);
      vi.S = (given != wave_port_overlaps.end()) ? given->second : data.GetSParameter(*E);
      if (data.HasVoltageCoords())
      {
        // The power and the voltage need the field at the port, which a port with a given
        // overlap (of the environment of a substructuring model) has not: NaN, not 0.
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        vi.P = (given != wave_port_overlaps.end()) ? std::complex<double>(nan, nan)
                                                   : data.GetPower(*E, *B);
        vi.V = (given != wave_port_overlaps.end()) ? std::complex<double>(nan, nan)
                                                   : data.GetVoltage(*E);
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

    // Don't measure S-Matrix unless there is only one excitation per port. All port
    // families (lumped, wave, Floquet) normalize their incident waves to unit time-averaged
    // power, so mixed lumped/wave/Floquet S-parameters need no cross-family scaling.
    if (!fem_op->GetPortExcitations().HasSParameters())
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
    // line's characteristic impedance encoded in the stored unit-overlap mode:
    // |∫ E_mode×H_mode⋆ ⋅ n dS| = 1. The incident wave is sqrt(2) times this mode,
    // carrying 1 W of time-averaged power for a lossless propagating mode.
    //
    // GetSParameter() returns the b-amplitude in Kurokawa convention for both port types:
    //   - Lumped: V / V_ref with V_ref = sqrt(2 R_ref) the peak incident voltage of the
    //     unit-power wave (see LumpedPortData::GetReferenceVoltage)
    //   - Wave:   ∫(E×H_mode)·n dS / (sqrt(2) ∫(E_mode×H_mode)·n dS), the modal reaction
    //             overlap divided by the amplitude of the incident wave relative to the
    //             mode, so that the incident wave itself projects to 1. The reaction
    //             separates the incident and reflected modal amplitudes also for a lossy
    //             port mode, and is the power overlap ∫(E×H_mode⋆)·n dS / sqrt(2) for a
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

      Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                 format("S[{}][{}]", idx, drive_port_idx), vi.S.real(), vi.S.imag(),
                 Measurement::Magnitude(vi.S), Measurement::Phase(vi.S));
    }

    // Floquet port S-parameters (already post-processed in MeasureFloquetPorts).
    for (auto &[port_idx, S_map] : measurement_cache.floquet_port_s)
    {
      for (auto &[key, S] : S_map)
      {
        S /= drive_weight;
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
        idx,
        surf_post_op.GetSurfaceFlux(idx, E.get(), B.get(), measurement_cache.freq.real()),
        data.type});
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
  // with the ratio of the time-averaged interface energy to the time-averaged
  // electric-side energy U of the mode (Measurement::ElectricEnergyNormalization):
  //          p_mj = 1/4 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} / U,
  // with U = E_elec + E_cap for constant materials.
  measurement_cache.interface_eps_i.clear();
  if constexpr (HasEGridFunction<solver_t>())
  {
    // Domain and port energies must have been measured first. E_cap returns zero if the
    // solver does not support lumped ports.
    //
    // TODO: Should this not include other types of energy too (surface impedance case)?
    // Without U, the participations are not measured.
    const auto energy_electric_all = measurement_cache.ElectricEnergyNormalization();
    if (!energy_electric_all)
    {
      return;
    }

    measurement_cache.interface_eps_i.reserve(surf_post_op.eps_surfs.size());
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      auto energy = surf_post_op.GetInterfaceElectricFieldEnergy(
          idx, *E, measurement_cache.freq.real());

      auto energy_participation_p = energy / *energy_electric_all;
      auto loss_tangent_delta = surf_post_op.GetInterfaceLossTangent(idx);
      auto quality_factor_Q = (energy_participation_p == 0.0 || loss_tangent_delta == 0.0)
                                  ? mfem::infinity()
                                  : 1.0 / (loss_tangent_delta * energy_participation_p);

      measurement_cache.interface_eps_i.emplace_back(Measurement::InterfaceData{
          idx, energy, loss_tangent_delta, energy_participation_p, quality_factor_Q});
    }
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
std::optional<double> PostOperator<solver_t>::MeasureAndPrintAll(int ex_idx, int step,
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
    WriteMFEMGridFunctions(freq.real(), ind, ExcitationSubFolder(ex_idx));
    Mpi::Print(" Wrote fields to disk (grid function) at step {:d}\n", step + 1);
  }
  return measurement_cache.TotalFieldEnergy();
}

template <ProblemType solver_t>
void PostOperator<solver_t>::ConfigureReducedPostprocessing(const RomOperator &rom_op,
                                                            double tol)
  requires(solver_t == ProblemType::DRIVEN)
{
  BlockTimer bt0(Timer::POSTPRO);
  reduced_postprocessing_ready = false;
  reduced_postprocessing_checked.clear();
  // A tenth of the PROM tolerance, but not below the round-off of the reduced outputs.
  reduced_postprocessing_tol =
      std::max(0.1 * tol, 1.0e7 * std::numeric_limits<double>::epsilon());
  // The reduced magnetic form contains curl(E) only, while a Floquet wave vector adds a
  // frequency-dependent k × E correction to the full magnetic flux density.
  if (fem_op->GetMaterialOp().HasWaveVector() || !fem_op->GetFloquetPortOp().Empty() ||
      fem_op->GetSurfaceCurrentOp().Size() > 0 || !surf_post_op.flux_surfs.empty() ||
      !surf_post_op.eps_surfs.empty() || surf_post_op.farfield.size() > 0 ||
      !interp_op.GetProbes().empty() || !fem_op->GetPortExcitations().HasSParameters())
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
  // Lumped-port voltages are linear in E and the basis is real, so the voltages at any
  // frequency are the reduced solution applied to the voltages of the basis vectors.
  const auto &lumped_port_op = fem_op->GetLumpedPortOp();
  reduced_port_V.resize(lumped_port_op.Size(), n);
  for (std::size_t j = 0; j < n && lumped_port_op.Size() > 0; j++)
  {
    E->Real().SetFromTrueDofs(basis[j]);
    E->Real().ExchangeFaceNbrData();
    Eigen::Index i = 0;
    for (const auto &[idx, voltage] : lumped_port_op.GetVoltages(*E))
    {
      reduced_port_V(i++, j) = voltage.real();
    }
  }

  // The energy forms on the true dofs, the whole domain first. The electric forms end with
  // the unit mass forms of the frequency-dependent materials, over the whole domain and
  // then over each postprocessing domain, weighted by their susceptibility at each
  // frequency.
  std::vector<std::unique_ptr<ParOperator>> ops_E, ops_H;
  const auto &nd_fespace = fem_op->GetNDSpace();
  const auto &curl_fespace = fem_op->GetCurlSpace();
  ops_E.push_back(std::make_unique<ParOperator>(*dom_post_op.M_elec, nd_fespace));
  ops_H.push_back(std::make_unique<ParOperator>(*dom_post_op.M_mag, curl_fespace));
  for (const auto &[idx, ops] : dom_post_op.M_i)
  {
    ops_E.push_back(std::make_unique<ParOperator>(*ops.first, nd_fespace));
    ops_H.push_back(std::make_unique<ParOperator>(*ops.second, curl_fespace));
  }
  for (const auto &[m, op] : dom_post_op.M_elec_fd)
  {
    ops_E.push_back(std::make_unique<ParOperator>(*op, nd_fespace));
  }
  for (const auto &[idx, forms] : dom_post_op.M_elec_fd_i)
  {
    for (const auto &[m, op] : forms)
    {
      ops_E.push_back(std::make_unique<ParOperator>(*op, nd_fespace));
    }
  }
  auto pointers = [](const auto &ops, std::size_t first)
  {
    std::vector<const Operator *> p;
    for (std::size_t i = first; i < ops.size(); i++)
    {
      p.push_back(ops[i].get());
    }
    return p;
  };
  const MPI_Comm comm = fem_op->GetComm();
  auto forms_E = ProjectSymmetricOperators(comm, basis, pointers(ops_E, 0));

  // The curls of the basis vectors nearly cancel in the magnetic energy at low frequency,
  // so the energy is not evaluated from Vᵀ Cᵀ M C V but as |R y|² from the factorization
  // C V = Q R by CGS2, with Q orthonormal in the inner product of the magnetic energy and
  // without the columns numerically in the span of the previous ones. The domain energies
  // are (R y)ᴴ Qᵀ Mᵢ Q (R y).
  const auto &curl = fem_op->GetCurlMatrix();
  const auto &M_mag = *ops_H[0];
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(n, n);
  std::vector<Vector> Q;
  ComplexVector work;
  for (std::size_t j = 0; j < n; j++)
  {
    Q.emplace_back(curl.Height());
    Q.back().UseDevice(true);
    curl.Mult(basis[j], Q.back());
    const auto [norm0, norm] = linalg::OrthogonalizeColumnCGS(
        comm, std::span(Q), R.col(j).data(), true, &M_mag, &work);
    if (!(norm <= 100.0 * std::numeric_limits<double>::epsilon() * norm0))
    {
      R(Q.size() - 1, j) = norm;
      Q.back() *= 1.0 / norm;
    }
    else
    {
      Q.pop_back();
    }
  }
  const auto r = Q.size();
  reduced_curl_R = R.topRows(r);
  auto forms_H = (r > 0) ? ProjectSymmetricOperators(comm, Q, pointers(ops_H, 1))
                         : std::vector<Eigen::MatrixXd>(ops_H.size() - 1);
  reduced_energy_E = std::move(forms_E[0]);
  reduced_domain_energy.clear();
  for (std::size_t i = 0; const auto &[idx, ops] : dom_post_op.M_i)
  {
    reduced_domain_energy.emplace(
        idx, std::make_pair(std::move(forms_E[i + 1]), std::move(forms_H[i])));
    i++;
  }
  auto form_E = forms_E.begin() + 1 + dom_post_op.M_i.size();
  auto unpack_fd = [&form_E](const DomainPostOperator::FrequencyDependentForms &forms)
  {
    ReducedFrequencyDependentForms reduced;
    for (const auto &[m, op] : forms)
    {
      reduced.emplace_back(m, std::move(*form_E++));
    }
    return reduced;
  };
  reduced_energy_E_fd = unpack_fd(dom_post_op.M_elec_fd);
  reduced_domain_energy_E_fd.clear();
  for (const auto &[idx, forms] : dom_post_op.M_elec_fd_i)
  {
    reduced_domain_energy_E_fd.emplace(idx, unpack_fd(forms));
  }
  reduced_postprocessing_ready = true;
}

template <ProblemType solver_t>
bool PostOperator<solver_t>::ReducedPostprocessingNeedsField(int ex_idx) const
  requires(solver_t == ProblemType::DRIVEN)
{
  // The full solution is needed only for the check against the full evaluation at the
  // first point of each excitation.
  return !reduced_postprocessing_checked.contains(ex_idx);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureAndPrintReduced(int ex_idx, int step,
                                                    const ComplexVector *e,
                                                    std::complex<double> omega,
                                                    const Eigen::VectorXcd &y,
                                                    const RomOperator &rom_op)
  requires(solver_t == ProblemType::DRIVEN)
{
  BlockTimer bt0(Timer::POSTPRO);
  MFEM_VERIFY(reduced_postprocessing_ready && y.size() == reduced_energy_E.rows(),
              "Invalid reduced postprocessing state!");
  MFEM_VERIFY(e || !ReducedPostprocessingNeedsField(ex_idx),
              "Reduced postprocessing requires the expanded solution!");
  if (e)
  {
    SetEGridFunction(*e);
  }

  // The quadratic forms in real arithmetic, on the columns [Re{y}, Im{y}], with the
  // field-energy weight of the full evaluation.
  Eigen::MatrixX2d Y(y.size(), 2);
  Y << y.real(), y.imag();
  const Eigen::MatrixX2d RY = reduced_curl_R * Y;
  auto quadratic = [](const Eigen::MatrixXd &A, const Eigen::MatrixX2d &X)
  {
    return DomainPostOperator::EnergyWeight(true) *
           (X.array() * (A.selfadjointView<Eigen::Lower>() * X).array()).sum();
  };
  auto quadratic_E = [&](const Eigen::MatrixXd &A, const ReducedFrequencyDependentForms &fd,
                         DomainPostOperator::ElectricEnergy type)
  {
    double energy = quadratic(A, Y);
    for (const auto &[m, A_m] : fd)
    {
      energy += dom_post_op.GetFrequencyDependentEnergyWeight(m, omega.real(), type) *
                quadratic(A_m, Y);
    }
    return energy;
  };
  const double scale_H = 1.0 / (omega.real() * omega.real());
  // The same omissions as the full measurement (MeasureDomainFieldEnergy): without the
  // stored electric energy of a material, the electric energies with it and the electric
  // participations are left empty.
  const bool stored_E = dom_post_op.HasStoredElectricFieldEnergy();
  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.ex_idx = ex_idx;
  if (stored_E)
  {
    measurement_cache.domain_E_field_energy_all = quadratic_E(
        reduced_energy_E, reduced_energy_E_fd, DomainPostOperator::ElectricEnergy::STORED);
  }
  measurement_cache.domain_E_field_energy_material = quadratic_E(
      reduced_energy_E, reduced_energy_E_fd, DomainPostOperator::ElectricEnergy::MATERIAL);
  measurement_cache.domain_H_field_energy_all =
      DomainPostOperator::EnergyWeight(true) * RY.squaredNorm() * scale_H;
  for (const auto &[idx, forms] : reduced_domain_energy)
  {
    auto &d = measurement_cache.domain_E_field_energy_i.emplace_back(
        Measurement::DomainData{idx, std::nullopt, std::nullopt});
    if (dom_post_op.HasStoredDomainElectricFieldEnergy(idx))
    {
      d.energy = quadratic_E(forms.first, reduced_domain_energy_E_fd.at(idx),
                             DomainPostOperator::ElectricEnergy::STORED);
      if (stored_E)
      {
        d.participation_ratio =
            ParticipationRatio(*d.energy, *measurement_cache.domain_E_field_energy_all);
      }
    }
    const double energy_H = quadratic(forms.second, RY) * scale_H;
    measurement_cache.domain_H_field_energy_i.push_back(
        {idx, energy_H,
         ParticipationRatio(energy_H, measurement_cache.domain_H_field_energy_all)});
  }

  // Lumped-port voltages and wave-port S-parameters are linear in the solution. The wave
  // port S-parameter is sᵀe / (sqrt(2) R) for the port mode s and its reaction R, relative
  // to the incident wave of sqrt(2) times the mode (WavePortData::GetSParameter), and every
  // PROM solve caches the pairings sᵀv of the basis vectors with the modes of all wave
  // ports for its wave-port terms.
  // The voltages and overlaps given by SetLumpedPortVoltages and SetWavePortOverlaps (of
  // the environment of a substructured region) replace those of the solution, as in the
  // full measurement.
  const Eigen::MatrixX2d port_V = reduced_port_V * Y;
  std::map<int, std::complex<double>> port_voltages;
  for (Eigen::Index i = 0; const auto &[idx, data] : fem_op->GetLumpedPortOp())
  {
    port_voltages[idx] = {port_V(i, 0), port_V(i, 1)};
    i++;
  }
  for (const auto &[idx, given] : lumped_port_voltages)
  {
    port_voltages[idx] = given;
  }
  MeasureLumpedPorts(port_voltages);
  for (const auto &[idx, data] : fem_op->GetWavePortOp())
  {
    if (const auto given = wave_port_overlaps.find(idx); given != wave_port_overlaps.end())
    {
      measurement_cache.wave_port_vi[idx].S = given->second;
      continue;
    }
    const auto *sV = rom_op.GetWavePortModePairing(idx);
    MFEM_VERIFY(sV && sV->size() == y.size(),
                "Missing or stale PROM mode pairing for wave port " << idx << "!");
    measurement_cache.wave_port_vi[idx].S =
        sV->cwiseProduct(y).sum() /
        (electromagnetics::UnitPowerAmplitude(true) * data.modal_reaction);
  }
  MeasureSParameter();

  const std::complex<double> freq =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * std::numbers::pi);
  if (reduced_postprocessing_checked.contains(ex_idx))
  {
    post_op_csv.PrintAllCSVData(*this, measurement_cache, freq.real(), step, ex_idx);
    return;
  }

  // At the first point of each excitation, output the full evaluation and compare it with
  // the reduced one. Each kind of output is compared relative to its scale: the total
  // electric or magnetic energy for the energies, the largest value for the port values.
  // The reduced outputs need to match to a tenth of the PROM tolerance only (or to
  // round-off, about 1e-9, for tighter tolerances), while implementation errors give O(1)
  // differences. A NaN error fails the check.
  reduced_postprocessing_checked.insert(ex_idx);
  const Measurement reduced = std::move(measurement_cache);
  ComplexVector b(fem_op->GetCurlMatrix().Height());
  b.UseDevice(true);
  fem_op->GetCurlMatrix().Mult(e->Real(), b.Real());
  fem_op->GetCurlMatrix().Mult(e->Imag(), b.Imag());
  b *= -1.0 / (1i * omega);
  MeasureAndPrintAll(ex_idx, step, *e, b, omega);
  const auto &full = measurement_cache;
  auto worst = [](double a, double b) { return (std::isnan(a) || a >= b) ? a : b; };
  double rel_error = 0.0;
  auto compare = [&](const std::optional<double> &value, const std::optional<double> &exact,
                     double scale)
  {
    MFEM_VERIFY(value.has_value() == exact.has_value(),
                "Reduced and full postprocessing report different energies!");
    if (exact)
    {
      rel_error =
          worst(rel_error, std::abs(*value - *exact) / std::max(std::abs(scale), 1.0e-300));
    }
  };
  // The electric energies relative to the stored electric energy, or without one (a
  // frequency-dependent material without a stored energy) to the electric energy at the
  // permittivity of the materials.
  const double E_all =
      full.domain_E_field_energy_all.value_or(full.domain_E_field_energy_material);
  const double H_all = full.domain_H_field_energy_all;
  compare(reduced.domain_E_field_energy_all, full.domain_E_field_energy_all, E_all);
  compare(reduced.domain_E_field_energy_material, full.domain_E_field_energy_material,
          E_all);
  compare(reduced.domain_H_field_energy_all, H_all, H_all);
  for (std::size_t i = 0; i < full.domain_E_field_energy_i.size(); i++)
  {
    compare(reduced.domain_E_field_energy_i[i].energy,
            full.domain_E_field_energy_i[i].energy, E_all);
    compare(reduced.domain_H_field_energy_i[i].energy,
            full.domain_H_field_energy_i[i].energy, H_all);
  }
  auto compare_ports = [&](const auto &values, const auto &exact, auto get)
  {
    double error = 0.0, scale = 0.0;
    for (const auto &[idx, vi] : exact)
    {
      error = worst(error, std::abs(get(values.at(idx)) - get(vi)));
      scale = worst(scale, std::abs(get(vi)));
    }
    compare(error, 0.0, scale);
  };
  compare_ports(reduced.lumped_port_vi, full.lumped_port_vi,
                [](const auto &vi) { return vi.V; });
  compare_ports(reduced.wave_port_vi, full.wave_port_vi,
                [](const auto &vi) { return vi.S; });
  // Take the same decision on every process, since only the full path communicates.
  rel_error = std::isnan(rel_error) ? mfem::infinity() : rel_error;
  Mpi::GlobalMax(1, &rel_error, fem_op->GetComm());
  if (rel_error > reduced_postprocessing_tol)
  {
    reduced_postprocessing_ready = false;
    Mpi::Warning("Reduced postprocessing disagrees with full evaluation (relative error "
                 "{:.3e}); reverting to full postprocessing.\n",
                 rel_error);
  }
}

template <ProblemType solver_t>
std::optional<double> PostOperator<solver_t>::MeasureAndPrintAll(
    int step, const ComplexVector &e, const ComplexVector &b, std::complex<double> omega,
    double error_abs, double error_bkwd, int num_conv)
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
  return measurement_cache.TotalFieldEnergy();
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
  return measurement_cache.TotalFieldEnergy().value();
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
  return measurement_cache.TotalFieldEnergy().value();
}

template <ProblemType solver_t>
double PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &e,
                                                  const Vector &b, double time,
                                                  double J_coef)
  requires(solver_t == ProblemType::TRANSIENT)
{
  BlockTimer bt0(Timer::POSTPRO);
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
  return measurement_cache.TotalFieldEnergy().value();
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
std::optional<double>
PostOperator<solver_t>::MeasureDomainFieldEnergyOnly(const ComplexVector &e,
                                                     const ComplexVector &b, double omega)
  requires(solver_t == ProblemType::DRIVEN)
{
  SetEGridFunction(e);
  SetBGridFunction(b);
  measurement_cache.freq = omega;
  MeasureDomainFieldEnergy();
  Mpi::Barrier(fem_op->GetComm());

  // Return total domain energy for normalizing error indicator.
  return measurement_cache.TotalFieldEnergy();
}

template <ProblemType solver_t>
std::optional<double> PostOperator<solver_t>::MeasureAndPrintAll(
    int step, const ComplexVector &et, const ComplexVector &en, std::complex<double> kn,
    double omega, double error_abs, double error_bkwd, int num_conv)
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
