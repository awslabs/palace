// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "drivensolver.hpp"

#include <complex>
#include <cstddef>
#include <iostream>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "fem/singularfield.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/floquetportoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/romoperator.hpp"
#include "models/singularsurfacepostoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

void DrivenSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                              MPI_Comm comm) const
{
  BaseSolver::Preprocess(iodata, smesh, comm);
  singular_features.Preprocess(iodata, smesh, comm);
}

bool DrivenSolver::RequiresSourceSerialMeshMetadata() const
{
  return iodata.solver.singular_elements.Enabled();
}

void DrivenSolver::ProcessPartitionedMesh(const mfem::ParMesh &parallel_mesh,
                                          const mesh::PartitionMetadata &metadata) const
{
  singular_features.ProcessPartitionedMesh(iodata, parallel_mesh, metadata);
}

std::pair<ErrorIndicator, long long int>
DrivenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Set up the spatial discretization and frequency sweep.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh, singular_features.GetSheetFeatures(),
                         singular_features.GetLineFeatures(),
                         singular_features.GetSourceVertexIds());
  const auto &port_excitations = space_op.GetPortExcitations();
  SaveMetadata(port_excitations);

  const auto &omega_sample = iodata.solver.driven.sample_f;

  bool adaptive = (iodata.solver.driven.adaptive_tol > 0.0);
  if (adaptive && omega_sample.size() <= iodata.solver.driven.prom_indices.size() &&
      !iodata.solver.driven.adaptive_circuit_synthesis)
  {
    Mpi::Warning("Adaptive frequency sweep requires > {} total frequency samples!\n"
                 "Reverting to uniform sweep!\n",
                 iodata.solver.driven.prom_indices.size());
    adaptive = false;
  }
  SaveMetadata(space_op.GetNDSpaces());
  Mpi::Print("\nComputing {}frequency response for:\n{}", adaptive ? "adaptive fast " : "",
             port_excitations.FmtLog());

  std::size_t restart = iodata.solver.driven.restart;
  if (restart != 1)
  {
    std::size_t max_iter = omega_sample.size() * space_op.GetPortExcitations().Size();
    MFEM_VERIFY(
        restart - 1 < max_iter,
        fmt::format("\"Restart\" ({}) is greater than the number of total samples ({})!",
                    restart, max_iter));

    Mpi::Print("\nRestarting from solve {}", iodata.solver.driven.restart);
  }

  // Main frequency sweep loop.
  return {adaptive ? SweepAdaptive(space_op) : SweepUniform(space_op),
          space_op.GlobalTrueVSize()};
}

ErrorIndicator DrivenSolver::SweepUniform(SpaceOperator &space_op) const
{
  if (space_op.HasSingularEnrichment())
  {
    return SweepUniformSingular(space_op);
  }

  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogeneous, no special elimination is required on
  // the RHS. Assemble the linear system for the initial frequency (so we can call
  // KspSolver::SetOperators). Compute everything at the first frequency step.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto &Curl = space_op.GetCurlMatrix();

  // Set up the linear solver.
  // The operators are constructed for each frequency step and used to initialize the ksp.
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());

  // Set up RHS vector for the incident field at port boundaries, and the vector for the
  // first frequency step.
  ComplexVector RHS(Curl.Width()), E(Curl.Width()), B(Curl.Height());
  RHS.UseDevice(true);
  E.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    else
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
  };
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Main excitation and frequency loop.
  auto t0 = Timer::Now();
  std::size_t excitation_counter = 0;
  const std::size_t excitation_restart_counter =
      ((iodata.solver.driven.restart - 1) / omega_sample.size()) + 1;
  const std::size_t freq_restart_idx =
      (iodata.solver.driven.restart - 1) % omega_sample.size();
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (++excitation_counter < excitation_restart_counter)
    {
      continue;
    }
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 excitation_counter, port_excitations.Size());
    }
    // Switch paraview subfolders: one for each excitation, if nr_excitations > 1.
    post_op.InitializeParaviewDataCollection(excitation_idx);

    // Frequency loop.
    for (std::size_t omega_i =
             ((excitation_counter == excitation_restart_counter) ? freq_restart_idx : 0);
         omega_i < omega_sample.size(); omega_i++)
    {
      auto omega = omega_sample[omega_i];
      // Assemble frequency dependent matrices and initialize operators in linear
      // solver.
      auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                        K.get(), C.get(), M.get(), A2.get());
      auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
          1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, omega);
      ksp.SetOperators(*A, *P);

      Mpi::Print(
          "\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s{})\n",
          omega_i + 1, omega_sample.size(),
          iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI),
          Timer::Duration(Timer::Now() - t0).count(),
          (port_excitations.Size() > 1)
              ? fmt::format(", solve {:d}/{:d}",
                            1 + omega_i + (excitation_counter - 1) * omega_sample.size(),
                            omega_sample.size() * port_excitations.Size())
              : "");

      // Solve linear system.
      space_op.GetExcitationVector(excitation_idx, omega, RHS);

      Mpi::Print("\n");
      ksp.Mult(RHS, E);

      // Start Post-processing.
      BlockTimer bt0(Timer::POSTPRO);
      Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e})\n",
                 linalg::Norml2(space_op.GetComm(), E),
                 linalg::Norml2(space_op.GetComm(), RHS));

      // Compute B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E.Real(), B.Real());
      Curl.Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (space_op.GetMaterialOp().HasWaveVector())
      {
        // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
        // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
        floquet_corr->AddMult(
            E, B,
            space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
      }

      auto total_domain_energy =
          post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);

      // Calculate and record the error indicators.
      Mpi::Print(" Updating solution error estimates\n");
      AddEstimate(E, B, total_domain_energy, indicator);
    }

    // Final postprocessing & printing.
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(ksp);
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

ErrorIndicator DrivenSolver::SweepUniformSingular(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto K_energy = space_op.GetStiffnessMatrix<Operator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M_energy = space_op.GetMassMatrix<Operator>(Operator::DIAG_ZERO);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K->Real() && !K->Imag() && (!C || (C->Real() && !C->Imag())) && M->Real(),
              "Driven singular simulations require real stiffness and damping with a "
              "complex electric mass operator!");

  auto ksp = MakeSingularComplexKspSolver(iodata, space_op.GetComm());
  ComplexVector rhs(space_op.GetNDTrueVSize()), electric_field(space_op.GetNDTrueVSize()),
      residual(space_op.GetNDTrueVSize());
  rhs.UseDevice(true);
  electric_field.UseDevice(true);
  residual.UseDevice(true);
  electric_field = 0.0;

  const fem::singular::AdaptiveAssemblyOptions surface_options{
      iodata.solver.singular_elements.quadrature_order,
      iodata.solver.singular_elements.abs_tol, iodata.solver.singular_elements.rel_tol,
      iodata.solver.singular_elements.max_subdivisions};
  std::unique_ptr<fem::singular::EnrichedNDFieldEvaluator> tetrahedral_real_evaluator,
      tetrahedral_imaginary_evaluator;
  std::unique_ptr<fem::singular::TriangleEnrichedNDFieldEvaluator>
      triangular_real_evaluator, triangular_imaginary_evaluator;
  std::unique_ptr<TetrahedronSingularSurfacePostOperator> tetrahedral_surface_postoperator;
  std::unique_ptr<TriangleSingularSurfacePostOperator> triangular_surface_postoperator;
  if (space_op.GetMesh().Dimension() == 3)
  {
    const auto *topology = space_op.GetSingularDofTopology();
    MFEM_VERIFY(topology,
                "Three-dimensional driven singular surface postprocessing requires "
                "tetrahedral singular DOF topology!");
    tetrahedral_real_evaluator = std::make_unique<fem::singular::EnrichedNDFieldEvaluator>(
        *topology, space_op.GetSingularParallelNumbering(), space_op.GetNDSpace().Get());
    tetrahedral_imaginary_evaluator =
        std::make_unique<fem::singular::EnrichedNDFieldEvaluator>(
            *topology, space_op.GetSingularParallelNumbering(),
            space_op.GetNDSpace().Get());
    tetrahedral_surface_postoperator =
        std::make_unique<TetrahedronSingularSurfacePostOperator>(
            iodata.boundaries.postpro, space_op.GetMaterialOp(),
            space_op.GetNDSpace().Get());
  }
  else
  {
    const auto *topology = space_op.GetTriangleSingularDofTopology();
    MFEM_VERIFY(topology, "Two-dimensional driven singular surface postprocessing requires "
                          "triangular singular DOF topology!");
    triangular_real_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            *topology, space_op.GetSingularParallelNumbering(),
            space_op.GetNDSpace().Get());
    triangular_imaginary_evaluator =
        std::make_unique<fem::singular::TriangleEnrichedNDFieldEvaluator>(
            *topology, space_op.GetSingularParallelNumbering(),
            space_op.GetNDSpace().Get());
    triangular_surface_postoperator = std::make_unique<TriangleSingularSurfacePostOperator>(
        iodata.boundaries.postpro, space_op.GetMaterialOp(), space_op.GetNDSpace().Get());
  }

  TableWithCSVFile output, surface_output;
  if (root)
  {
    output = TableWithCSVFile(post_dir / "singular-driven.csv");
    output.table.insert("excitation", "Excitation", -1, 0, 0, "");
    output.table.insert("index", "Frequency index", -1, 0, 0, "");
    output.table.insert("frequency", "Frequency (GHz)");
    output.table.insert("electric_energy", "Electric field energy (J)");
    output.table.insert("magnetic_energy", "Magnetic field energy (J)");
    output.table.insert("relative_residual", "Relative residual");
    for (const auto &[port_index, port] : space_op.GetLumpedPortOp())
    {
      output.table.insert(fmt::format("voltage_real_{}", port_index),
                          fmt::format("Re{{V[{}]}} (V)", port_index));
      output.table.insert(fmt::format("voltage_imag_{}", port_index),
                          fmt::format("Im{{V[{}]}} (V)", port_index));
      output.table.insert(fmt::format("current_real_{}", port_index),
                          fmt::format("Re{{I[{}]}} (A)", port_index));
      output.table.insert(fmt::format("current_imag_{}", port_index),
                          fmt::format("Im{{I[{}]}} (A)", port_index));
      output.table.insert(fmt::format("s_real_{}", port_index),
                          fmt::format("Re{{S[{}]}}", port_index));
      output.table.insert(fmt::format("s_imag_{}", port_index),
                          fmt::format("Im{{S[{}]}}", port_index));
    }
    output.table[0].print_as_int = true;
    output.table[1].print_as_int = true;
    output.WriteFullTableTrunc();
    const bool has_surface = tetrahedral_surface_postoperator
                                 ? !tetrahedral_surface_postoperator->Empty()
                                 : !triangular_surface_postoperator->Empty();
    if (has_surface)
    {
      surface_output = TableWithCSVFile(post_dir / "surface-Q.csv");
      surface_output.table.insert("excitation", "Excitation", -1, 0, 0, "");
      surface_output.table.insert("frequency", "f (GHz)");
      surface_output.table[0].print_as_int = true;
      for (const auto &[index, data] : iodata.boundaries.postpro.dielectric)
      {
        surface_output.table.insert(fmt::format("p_{}", index),
                                    fmt::format("p_surf[{}]", index));
        surface_output.table.insert(fmt::format("Q_{}", index),
                                    fmt::format("Q_surf[{}]", index));
      }
      surface_output.WriteFullTableTrunc();
    }
  }

  auto t0 = Timer::Now();
  std::size_t excitation_counter = 0;
  const std::size_t excitation_restart_counter =
      ((iodata.solver.driven.restart - 1) / omega_sample.size()) + 1;
  const std::size_t frequency_restart_index =
      (iodata.solver.driven.restart - 1) % omega_sample.size();
  for (const auto &[excitation_index, excitation_spec] : port_excitations)
  {
    if (++excitation_counter < excitation_restart_counter)
    {
      continue;
    }
    for (std::size_t frequency_index =
             excitation_counter == excitation_restart_counter ? frequency_restart_index : 0;
         frequency_index < omega_sample.size(); frequency_index++)
    {
      const double omega = omega_sample[frequency_index];
      MFEM_VERIFY(omega > 0.0 && std::isfinite(omega),
                  "Driven singular simulations require positive finite frequencies!");
      auto extra = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      MFEM_VERIFY(!extra,
                  "Driven singular simulations do not support extra boundary operators!");
      auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i,
                                        K.get(), C.get(), M.get());
      ksp->SetOperators(*A, *A);

      Mpi::Print("\nIt {:d}/{:d}: omega/2pi = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
                 frequency_index + 1, omega_sample.size(),
                 iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                     (2 * M_PI),
                 Timer::Duration(Timer::Now() - t0).count());
      space_op.GetExcitationVector(excitation_index, omega, rhs);
      ksp->Mult(rhs, electric_field);

      A->Mult(electric_field, residual);
      linalg::AXPY(-1.0, rhs, residual);
      const double relative_residual = linalg::Norml2(space_op.GetComm(), residual) /
                                       std::max(linalg::Norml2(space_op.GetComm(), rhs),
                                                std::numeric_limits<double>::min());
      const auto energy = MeasureSingularFullWaveEnergy(space_op.GetComm(), *M_energy,
                                                        *K_energy, electric_field, omega);
      std::vector<TriangleSingularSurfacePostOperator::Measurement> surface_measurements;
      if (tetrahedral_surface_postoperator && !tetrahedral_surface_postoperator->Empty())
      {
        tetrahedral_real_evaluator->SetFromTrueDofs(electric_field.Real());
        tetrahedral_imaginary_evaluator->SetFromTrueDofs(electric_field.Imag());
        surface_measurements = tetrahedral_surface_postoperator->Measure(
            *tetrahedral_real_evaluator, *tetrahedral_imaginary_evaluator, energy.electric,
            surface_options);
      }
      else if (triangular_surface_postoperator && !triangular_surface_postoperator->Empty())
      {
        triangular_real_evaluator->SetFromTrueDofs(electric_field.Real());
        triangular_imaginary_evaluator->SetFromTrueDofs(electric_field.Imag());
        surface_measurements = triangular_surface_postoperator->Measure(
            *triangular_real_evaluator, *triangular_imaginary_evaluator, energy.electric,
            surface_options);
      }
      const double electric_energy =
          iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy.electric);
      const double magnetic_energy =
          iodata.units.Dimensionalize<Units::ValueType::ENERGY>(energy.magnetic);
      struct LumpedPortMeasurement
      {
        int index;
        std::complex<double> voltage;
        std::complex<double> current;
        std::complex<double> scattering;
      };
      std::vector<LumpedPortMeasurement> lumped_port_measurements;
      const auto [drive_is_simple, drive_port_type, drive_port_index] =
          excitation_spec.IsSimple();
      for (const auto &[port_index, port] : space_op.GetLumpedPortOp())
      {
        const auto voltage =
            space_op.GetSingularLumpedPortVoltage(port_index, electric_field);
        const auto current = voltage / port.GetCharacteristicImpedance(omega);
        auto scattering =
            space_op.GetSingularLumpedPortSParameter(port_index, electric_field);
        if (drive_is_simple && drive_port_type == PortType::LumpedPort &&
            drive_port_index == port_index)
        {
          scattering.real(scattering.real() - 1.0);
        }
        lumped_port_measurements.push_back({port_index, voltage, current, scattering});
        Mpi::Print(" Port {:d}: V = {:.6e}{:+.6e}i, I = {:.6e}{:+.6e}i, "
                   "S = {:.6e}{:+.6e}i\n",
                   port_index, voltage.real(), voltage.imag(), current.real(),
                   current.imag(), scattering.real(), scattering.imag());
      }
      Mpi::Print(" Sol. ||E|| = {:.6e} (||RHS|| = {:.6e}, rel. residual = {:.3e})\n"
                 " Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n",
                 linalg::Norml2(space_op.GetComm(), electric_field),
                 linalg::Norml2(space_op.GetComm(), rhs), relative_residual,
                 electric_energy, magnetic_energy, electric_energy + magnetic_energy);

      if (root)
      {
        output.table["excitation"] << excitation_index;
        output.table["index"] << frequency_index + 1;
        output.table["frequency"]
            << iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) / (2 * M_PI);
        output.table["electric_energy"] << electric_energy;
        output.table["magnetic_energy"] << magnetic_energy;
        output.table["relative_residual"] << relative_residual;
        for (const auto &measurement : lumped_port_measurements)
        {
          const auto voltage =
              iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(measurement.voltage);
          const auto current =
              iodata.units.Dimensionalize<Units::ValueType::CURRENT>(measurement.current);
          output.table[fmt::format("voltage_real_{}", measurement.index)] << voltage.real();
          output.table[fmt::format("voltage_imag_{}", measurement.index)] << voltage.imag();
          output.table[fmt::format("current_real_{}", measurement.index)] << current.real();
          output.table[fmt::format("current_imag_{}", measurement.index)] << current.imag();
          output.table[fmt::format("s_real_{}", measurement.index)]
              << measurement.scattering.real();
          output.table[fmt::format("s_imag_{}", measurement.index)]
              << measurement.scattering.imag();
        }
        output.WriteFullTableTrunc();
        if (!surface_measurements.empty())
        {
          surface_output.table["excitation"] << excitation_index;
          surface_output.table["frequency"]
              << iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                     (2 * M_PI);
          for (const auto &measurement : surface_measurements)
          {
            surface_output.table[fmt::format("p_{}", measurement.index)]
                << measurement.participation;
            surface_output.table[fmt::format("Q_{}", measurement.index)]
                << measurement.quality_factor;
          }
          surface_output.WriteFullTableTrunc();
        }
      }
    }
  }

  SaveMetadata(*ksp);
  SaveMetadata(
      "SingularFullWave",
      nlohmann::json{
          {"Enabled", true},
          {"Dimension", space_op.GetMesh().Dimension()},
          {"Output", "singular-driven.csv"},
          {"ElectricEnergy", "0.5 E^H M_epsilon E"},
          {"MagneticEnergy", "0.5 E^H K_mu^-1 E / |omega|^2"},
          {"BulkDielectricLoss", space_op.GetMaterialOp().HasLossTangent()},
          {"LumpedPorts", space_op.GetLumpedPortOp().Size()},
          {"FieldGridOutput", false},
          {"ErrorEstimator", false},
          {"SurfaceIntegrability", space_op.GetMesh().Dimension() == 3
                                       ? GetSingularSurfaceIntegrabilityMetadata(
                                             *singular_features.GetSheetFeatures())
                                       : GetSingularSurfaceIntegrabilityMetadata(
                                             *singular_features.GetLineFeatures())},
          {"SurfaceParticipation", GetSingularSurfaceParticipationMetadata(iodata)}});
  return {};
}

ErrorIndicator DrivenSolver::SweepAdaptive(SpaceOperator &space_op) const
{
  const auto &port_excitations = space_op.GetPortExcitations();
  const auto &omega_sample = iodata.solver.driven.sample_f;

  // Initialize postprocessing for measurement and printers.
  // Initialize write directory with default path; will be changed for multi-excitations.
  PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

  // Configure PROM parameters if not specified.
  double offline_tol = iodata.solver.driven.adaptive_tol;
  std::size_t convergence_memory = iodata.solver.driven.adaptive_memory;
  std::size_t max_size_per_excitation = iodata.solver.driven.adaptive_max_size;
  std::size_t nprom_indices = iodata.solver.driven.prom_indices.size();
  MFEM_VERIFY(max_size_per_excitation <= 0 || max_size_per_excitation >= nprom_indices,
              "Adaptive frequency sweep must sample at least " << nprom_indices
                                                               << " frequency points!");

  // Allocate negative curl matrix for postprocessing the B-field and vectors for the
  // high-dimensional field solution.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector E(Curl.Width()), Eh(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  Eh.UseDevice(true);
  B.UseDevice(true);
  E = 0.0;
  Eh = 0.0;
  B = 0.0;

  // Initialize structures for storing and reducing the results of error estimation.
  const bool is_2d = (space_op.GetNDSpace().Dimension() < 3);
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator_3d;
  std::unique_ptr<BoundaryModeFluxErrorEstimator<ComplexVector>> estimator_2d;
  if (is_2d)
  {
    estimator_2d = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        space_op.GetCurlSpace(), space_op.GetH1Spaces(), iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg);
  }
  else
  {
    estimator_3d = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
        iodata.solver.linear.estimator_mg);
  }
  auto AddEstimate =
      [&](const ComplexVector &E, const ComplexVector &B, double Et, ErrorIndicator &ind)
  {
    if (is_2d)
      estimator_2d->AddErrorIndicator(E, B, Et, ind);
    else
      estimator_3d->AddErrorIndicator(E, B, Et, ind);
  };
  ErrorIndicator indicator;

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Configure the PROM operator which performs the parameter space sampling and basis
  // construction during the offline phase as well as the PROM solution during the online
  // phase.
  auto t0 = Timer::Now();
  const double unit_GHz =
      iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) / (2 * M_PI);
  Mpi::Print("\nBeginning PROM construction offline phase:\n"
             " {:d} points for frequency sweep over [{:.3e}, {:.3e}] GHz\n",
             omega_sample.size(), omega_sample.front() * unit_GHz,
             omega_sample.back() * unit_GHz);
  RomOperator prom_op(iodata, space_op, max_size_per_excitation);
  space_op.GetWavePortOp().SetSuppressOutput(true);

  // Add ports to PROM if we do synthesis.
  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.AddLumpedPortModesForSynthesis();
    if (space_op.GetWavePortOp().Size() > 0)
    {
      // Use the band center as the reference frequency for seeding wave-port modes.
      // The choice rescales the basis vector but does not change correctness.
      const double omega_ref = 0.5 * (omega_sample.front() + omega_sample.back());
      prom_op.AddWavePortModesForSynthesis(omega_ref);
    }
  }

  // Initialize the basis with samples from the top and bottom of the frequency
  // range of interest. Each call for an HDM solution adds the frequency sample to P_S and
  // removes it from P \ P_S. Timing for the HDM construction and solve is handled inside
  // of the RomOperator.
  auto UpdatePROM = [&](int excitation_idx, double omega, std::size_t sample_idx)
  {
    // Add the HDM solution to the PROM reduced basis.
    prom_op.UpdatePROM(E, fmt::format("sample_e{:d}_s{:d}", excitation_idx, sample_idx));
    prom_op.UpdateMRI(excitation_idx, omega, E);

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for energy postprocessing and error estimation.
    BlockTimer bt0(Timer::POSTPRO);
    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
      // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k₀, so scale = 1.
      floquet_corr->AddMult(
          E, B, space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
    }

    // Measure domain energies for the error indicator only. Don't exchange face_nbr_data,
    // unless printing paraview fields.
    auto total_domain_energy = post_op.MeasureDomainFieldEnergyOnly(E, B);
    AddEstimate(E, B, total_domain_energy, indicator);
  };

  // Loop excitations to add to PROM.
  //
  // Restart should not really be used for adaptive sweeps, but must work. Construct PROM in
  // the same way same regardless of restart for consistency. Don't shift excitation start.
  int excitation_counter = 0;
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nAdding excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++excitation_counter, port_excitations.Size());
    }
    prom_op.SetExcitationIndex(excitation_idx);  // Pre-compute RHS1

    // Initialize PROM with explicit HDM samples, record the estimate but do not act on it.
    std::vector<double> max_errors;
    std::size_t counter_rom_sample = 0;
    for (auto i : iodata.solver.driven.prom_indices)
    {
      auto omega = omega_sample[i];
      prom_op.SolveHDM(excitation_idx, omega, E);
      prom_op.SolvePROM(excitation_idx, omega, Eh);
      linalg::AXPY(-1.0, E, Eh);
      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      UpdatePROM(excitation_idx, omega, counter_rom_sample);
      counter_rom_sample++;
    }
    // The estimates associated to the end points are assumed inaccurate.
    max_errors[0] = std::numeric_limits<double>::infinity();
    max_errors[1] = std::numeric_limits<double>::infinity();
    auto memory = std::distance(max_errors.rbegin(),
                                std::find_if(max_errors.rbegin(), max_errors.rend(),
                                             [=](auto x) { return x > offline_tol; }));
    memory = std::max(0L, memory);  // Ensure memory >= 0 as it should be.

    // Greedy procedure for basis construction (offline phase). Basis is initialized with
    // solutions at frequency sweep endpoints and explicit sample frequencies.
    std::size_t it = max_errors.size();
    for (std::size_t it0 = it; it < max_size_per_excitation && memory < convergence_memory;
         it++)
    {
      // Compute the location of the maximum error in parameter domain (bounded by the
      // previous samples).
      double omega_star = prom_op.FindMaxError(excitation_idx)[0];

      // Sample HDM and add solution to basis.
      prom_op.SolveHDM(excitation_idx, omega_star, E);
      prom_op.SolvePROM(excitation_idx, omega_star, Eh);
      linalg::AXPY(-1.0, E, Eh);

      max_errors.push_back(linalg::Norml2(space_op.GetComm(), Eh) /
                           linalg::Norml2(space_op.GetComm(), E));
      memory = max_errors.back() < offline_tol ? memory + 1 : 0;

      Mpi::Print("\nGreedy iteration {:d} (n = {:d}): ω* = {:.3e} GHz ({:.3e}), error = "
                 "{:.3e}, memory = {:d}/{:d}\n",
                 it - it0 + 1, prom_op.GetReducedDimension(), omega_star * unit_GHz,
                 omega_star, max_errors.back(), memory, convergence_memory);
      UpdatePROM(excitation_idx, omega_star, counter_rom_sample);
      counter_rom_sample++;
    }
    Mpi::Print("\nAdaptive sampling{} {:d} frequency samples:\n"
               " n = {:d}, error = {:.3e}, tol = {:.3e}, memory = {:d}/{:d}\n",
               (it == max_size_per_excitation) ? " reached maximum" : " converged with", it,
               prom_op.GetReducedDimension(), max_errors.back(), offline_tol, memory,
               convergence_memory);
    utils::PrettyPrint(prom_op.GetSamplePoints(excitation_idx), unit_GHz,
                       " Sampled frequencies (GHz):");
    utils::PrettyPrint(max_errors, 1.0, " Sample errors:");
  }

  Mpi::Print(" Total offline phase elapsed time: {:.2e} s\n",
             Timer::Duration(Timer::Now() - t0).count());  // Timing on root

  if (iodata.solver.driven.adaptive_circuit_synthesis)
  {
    prom_op.PrintPROMMatrices(iodata.units, iodata.problem.output);
  }

  // Main fast frequency sweep loop (online phase).
  Mpi::Print("\nBeginning fast frequency sweep online phase\n");
  space_op.GetWavePortOp().SetSuppressOutput(false);  // Disable output suppression
  excitation_counter = 0;
  for (const auto &[excitation_idx, excitation_spec] : port_excitations)
  {
    if (port_excitations.Size() > 1)
    {
      Mpi::Print("\nSweeping excitation index {:d} ({:d}/{:d}):\n", excitation_idx,
                 ++excitation_counter, port_excitations.Size());
    }
    // Switch paraview subfolders: one for each excitation, if nr_excitations > 1.
    post_op.InitializeParaviewDataCollection(excitation_idx);

    // Frequency loop.
    for (std::size_t omega_i = 0; omega_i < omega_sample.size(); omega_i++)
    {
      auto omega = omega_sample[omega_i];
      Mpi::Print("\nIt {:d}/{:d}: ω/2π = {:.3e} GHz (total elapsed time = {:.2e} s)\n",
                 omega_i + 1, omega_sample.size(),
                 iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                     (2 * M_PI),
                 Timer::Duration(Timer::Now() - t0).count());

      // Assemble and solve the PROM linear system.
      prom_op.SolvePROM(excitation_idx, omega, E);
      Mpi::Print("\n");

      // Start Post-processing.
      BlockTimer bt0(Timer::POSTPRO);
      Mpi::Print(" Sol. ||E|| = {:.6e}\n", linalg::Norml2(space_op.GetComm(), E));

      // Compute B = -1/(iω) ∇ x E on the true dofs.
      Curl.Mult(E.Real(), B.Real());
      Curl.Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      if (space_op.GetMaterialOp().HasWaveVector())
      {
        // Calculate B field correction for Floquet BCs: B += k_F(ω)/ω × E.
        // With k₀ = k_F_ref/ω_ref stored, k_F(ω)/ω = k_F_ref/ω_ref = k₀, so scale = 1.
        floquet_corr->AddMult(
            E, B,
            space_op.GetMaterialOp().HasFloquetFrequencyScaling() ? 1.0 : 1.0 / omega);
      }
      post_op.MeasureAndPrintAll(excitation_idx, int(omega_i), E, B, omega);
    }

    // Final postprocessing & printing: no change to indicator since these are in PROM.
    BlockTimer bt0(Timer::POSTPRO);
    SaveMetadata(prom_op.GetLinearSolver());
  }
  post_op.MeasureFinalize(indicator);
  return indicator;
}

}  // namespace palace
