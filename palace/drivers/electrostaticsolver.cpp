// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "electrostaticsolver.hpp"

#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfaceresponseoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
ElectrostaticSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct the system matrix defining the linear operator. Dirichlet boundaries are
  // handled eliminating the rows and columns of the system matrix for the corresponding
  // dofs. The eliminated matrix is stored in order to construct the RHS vector for nonzero
  // prescribed BC values.
  BlockTimer bt0(Timer::CONSTRUCT);
  LaplaceOperator laplace_op(iodata, mesh);
  auto K = laplace_op.GetStiffnessMatrix();
  std::unique_ptr<SurfaceResponseOperator> response_correction;
  std::unique_ptr<SumOperator> corrected_K;
  const Operator *system_K = K.get();
  if (iodata.solver.electrostatic.response_correction && final_postprocessing_pass)
  {
    response_correction = std::make_unique<SurfaceResponseOperator>(iodata, laplace_op);
    corrected_K = std::make_unique<SumOperator>(*K, *response_correction);
    system_K = corrected_K.get();
  }
  const auto &Grad = laplace_op.GetGradMatrix();
  SaveMetadata(laplace_op.GetH1Spaces());

  // Preserve the historical thin-metal solve and outputs even when response correction is
  // enabled. The corrected solve uses the same assembled thin operator as its
  // preconditioner.
  KspSolver ksp(iodata, laplace_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);

  // Source indices are either equipotential terminals or prescribed potential traces.
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
  int n_step = static_cast<int>(laplace_op.GetSources().size());
  MFEM_VERIFY(n_step > 0,
              "No terminal or prescribed potential boundaries specified for electrostatic "
              "simulation!");

  // Right-hand side term and solution vector storage.
  Vector RHS(Grad.Width()), E(Grad.Height()), D(laplace_op.GetRTSpace().GetTrueVSize());
  D.UseDevice(true);
  std::vector<Vector> V(n_step);
  std::vector<Vector> V_corrected(response_correction ? n_step : 0);
  std::vector<Vector> D_basis(post_op.NeedsRecoveredElectricFlux() &&
                                      iodata.solver.electrostatic.response_matrix
                                  ? n_step
                                  : 0);
  using EnergyData = PostOperator<ProblemType::ELECTROSTATIC>::ElectrostaticEnergyData;
  struct CorrectedResult
  {
    int source;
    EnergyData raw;
    EnergyData postprocessed_fixed_trace;
    EnergyData postprocessed_fixed_flux;
    EnergyData corrected;
    std::map<int, double> trace_closure_spread;
  };
  std::vector<CorrectedResult> corrected_results;
  corrected_results.reserve(response_correction ? n_step : 0);

  // Initialize structures for storing and reducing the results of error estimation.
  GradFluxErrorEstimator estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Main loop over terminal boundaries.
  Mpi::Print("\nComputing electrostatic fields for {:d} {}\n", n_step,
             (n_step > 1) ? "sources" : "source");
  int step = 0;
  auto t0 = Timer::Now();
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    Mpi::Print("\nIt {:d}/{:d}: Index = {:d} (elapsed time = {:.2e} s)\n", step + 1, n_step,
               idx, Timer::Duration(Timer::Now() - t0).count());

    // Form and solve the linear system for a prescribed nonzero voltage on the specified
    // terminal.
    Mpi::Print("\n");
    laplace_op.GetExcitationVector(idx, *K, V[step], RHS);
    Vector corrected_rhs;
    if (response_correction)
    {
      V_corrected[step] = V[step];
      corrected_rhs = RHS;
      response_correction->EliminateRHS(V_corrected[step], corrected_rhs);
    }
    ksp.Mult(RHS, V[step]);

    // Start Post-processing.
    BlockTimer bt2(Timer::POSTPRO);
    Mpi::Print(" Sol. ||V|| = {:.6e} (||RHS|| = {:.6e})\n",
               linalg::Norml2(laplace_op.GetComm(), V[step]),
               linalg::Norml2(laplace_op.GetComm(), RHS));

    // Compute E = -∇V on the true dofs.
    E = 0.0;
    Grad.AddMult(V[step], E, -1.0);

    if (post_op.NeedsRecoveredElectricFlux())
    {
      Mpi::Print(" Recovering electric flux for interface postprocessing\n");
      estimator.RecoverFlux(E, D);
      post_op.SetRecoveredElectricFlux(D);
      if (!D_basis.empty())
      {
        D_basis[step] = D;
      }
    }

    // Measurement and printing.
    auto total_domain_energy = post_op.MeasureAndPrintAll(step, V[step], E, idx);

    EnergyData raw_energies;
    if (response_correction)
    {
      raw_energies = post_op.GetElectrostaticEnergies(
          V[step], E, post_op.NeedsRecoveredElectricFlux() ? &D : nullptr);
    }

    if (response_correction)
    {
      const auto response = response_correction->GetElectrostaticResponse(V[step]);
      auto ApplyResponse = [&](EnergyData energies, double domain_correction,
                               const std::map<int, double> &fabricated_surface)
      {
        energies.domain += domain_correction;
        for (const auto &[interface, energy] : fabricated_surface)
        {
          auto it = energies.interfaces.find(interface);
          MFEM_VERIFY(it != energies.interfaces.end(),
                      "Response correction refers to target interface "
                          << interface << " which is not configured for postprocessing!");
          MFEM_VERIFY(
              !it->second.edge_energies.empty(),
              "Response-corrected target interface "
                  << interface << " requires EdgeDistances and EdgeAttributes or AutomaticEdges!");
          // EdgeDistances are sorted. The largest configured radius is the matching
          // distance of the coupon response model.
          it->second.energy = it->second.edge_energies.back().energy_outside + energy;
        }
        MFEM_VERIFY(energies.domain > 0.0,
                    "Response-corrected electrostatic energy is not positive!");
        return energies;
      };
      auto postprocessed_fixed_trace = ApplyResponse(
          raw_energies, response.domain_correction, response.fabricated_surface_energy);
      auto postprocessed_fixed_flux =
          ApplyResponse(raw_energies, response.domain_correction_fixed_flux,
                        response.fabricated_surface_energy_fixed_flux);
      constexpr double maximum_trace_closure_spread = 0.05;
      if (response.maximum_trace_closure_spread > maximum_trace_closure_spread)
      {
        Mpi::Warning(
            "Electrostatic postprocessing-only surface-response trace-closure spread "
            "exceeds 5% ({:.3e}). Corrected values are reported, but the raw thin-metal "
            "field does not determine a closure-independent local response.\n",
            response.maximum_trace_closure_spread);
      }

      Mpi::Print(" Solving fabrication-response corrected field\n");
      V_corrected[step] = V[step];
      ksp.SetOperator(*system_K);
      ksp.SetInitialGuess(true);
      ksp.Mult(corrected_rhs, V_corrected[step]);
      ksp.SetInitialGuess(iodata.solver.linear.initial_guess);
      ksp.SetOperator(*K);
      Vector E_corrected(Grad.Height()), D_corrected;
      E_corrected = 0.0;
      Grad.AddMult(V_corrected[step], E_corrected, -1.0);
      const Vector *D_corrected_ptr = nullptr;
      if (post_op.NeedsRecoveredElectricFlux())
      {
        D_corrected.SetSize(laplace_op.GetRTSpace().GetTrueVSize());
        D_corrected.UseDevice(true);
        estimator.RecoverFlux(E_corrected, D_corrected);
        D_corrected_ptr = &D_corrected;
      }

      auto corrected_energies =
          post_op.GetElectrostaticEnergies(V_corrected[step], E_corrected, D_corrected_ptr);
      const auto corrected_response =
          response_correction->GetElectrostaticResponse(V_corrected[step]);
      corrected_energies =
          ApplyResponse(std::move(corrected_energies), corrected_response.domain_correction,
                        corrected_response.fabricated_surface_energy);

      if (response_correction->HasSurfaceResponse())
      {
        corrected_results.push_back(CorrectedResult{
            idx, std::move(raw_energies), std::move(postprocessed_fixed_trace),
            std::move(postprocessed_fixed_flux), std::move(corrected_energies),
            response.trace_closure_spread});
      }
    }

    // Keep AMR driven by the historical raw thin-metal solution. This ensures enabling
    // response correction does not alter the mesh sequence or any ordinary output.
    Mpi::Print(" Updating solution error estimates\n");
    if (post_op.NeedsRecoveredElectricFlux())
    {
      estimator.AddErrorIndicator(E, D, total_domain_energy, indicator);
    }
    else
    {
      estimator.AddErrorIndicator(E, total_domain_energy, indicator);
    }

    // Next terminal.
    step++;
  }

  // Postprocess the capacitance matrix only for equipotential terminal solutions.
  BlockTimer bt1(Timer::POSTPRO);
  SaveMetadata(ksp);
  if (iodata.boundaries.prescribed_potential.empty())
  {
    PostprocessTerminals(post_op, laplace_op.GetSources(), V);
  }
  else if (iodata.solver.electrostatic.response_matrix)
  {
    PostprocessResponseMatrix(post_op, laplace_op, Grad, V, D_basis);
  }

  if (root && !corrected_results.empty())
  {
    using VT = Units::ValueType;
    TableWithCSVFile output(post_dir / "surface-Q-corrected.csv");
    output.table.insert(Column("source", "i", 0, 0, 2, ""));
    output.table.insert("domain_raw", "E_elec raw (J)");
    output.table.insert("domain_postprocessed_fixed_trace",
                        "E_elec postprocessed fixed-trace (J)");
    output.table.insert("domain_postprocessed_fixed_flux",
                        "E_elec postprocessed fixed-flux (J)");
    output.table.insert("domain_corrected", "E_elec corrected (J)");
    const auto &interfaces = corrected_results.front().raw.interfaces;
    for (const auto &[interface, data] : interfaces)
    {
      output.table.insert(fmt::format("energy_raw_{}", interface),
                          fmt::format("E_surf raw[{}] (J)", interface));
      output.table.insert(fmt::format("participation_raw_{}", interface),
                          fmt::format("p_surf raw[{}]", interface));
      output.table.insert(fmt::format("quality_raw_{}", interface),
                          fmt::format("Q_surf raw[{}]", interface));
      output.table.insert(
          fmt::format("energy_postprocessed_fixed_trace_{}", interface),
          fmt::format("E_surf postprocessed fixed-trace[{}] (J)", interface));
      output.table.insert(
          fmt::format("participation_postprocessed_fixed_trace_{}", interface),
          fmt::format("p_surf postprocessed fixed-trace[{}]", interface));
      output.table.insert(fmt::format("quality_postprocessed_fixed_trace_{}", interface),
                          fmt::format("Q_surf postprocessed fixed-trace[{}]", interface));
      output.table.insert(
          fmt::format("energy_postprocessed_fixed_flux_{}", interface),
          fmt::format("E_surf postprocessed fixed-flux[{}] (J)", interface));
      output.table.insert(
          fmt::format("participation_postprocessed_fixed_flux_{}", interface),
          fmt::format("p_surf postprocessed fixed-flux[{}]", interface));
      output.table.insert(fmt::format("quality_postprocessed_fixed_flux_{}", interface),
                          fmt::format("Q_surf postprocessed fixed-flux[{}]", interface));
      output.table.insert(fmt::format("energy_corrected_{}", interface),
                          fmt::format("E_surf corrected[{}] (J)", interface));
      output.table.insert(fmt::format("participation_corrected_{}", interface),
                          fmt::format("p_surf corrected[{}]", interface));
      output.table.insert(fmt::format("quality_corrected_{}", interface),
                          fmt::format("Q_surf corrected[{}]", interface));
      output.table.insert(fmt::format("trace_closure_spread_{}", interface),
                          fmt::format("trace closure spread[{}]", interface));
    }
    for (const auto &result : corrected_results)
    {
      output.table["source"] << result.source;
      output.table["domain_raw"]
          << iodata.units.Dimensionalize<VT::ENERGY>(result.raw.domain);
      output.table["domain_postprocessed_fixed_trace"]
          << iodata.units.Dimensionalize<VT::ENERGY>(
                 result.postprocessed_fixed_trace.domain);
      output.table["domain_postprocessed_fixed_flux"]
          << iodata.units.Dimensionalize<VT::ENERGY>(
                 result.postprocessed_fixed_flux.domain);
      output.table["domain_corrected"]
          << iodata.units.Dimensionalize<VT::ENERGY>(result.corrected.domain);
      MFEM_VERIFY(
          result.raw.interfaces.size() == interfaces.size() &&
              result.postprocessed_fixed_trace.interfaces.size() == interfaces.size() &&
              result.postprocessed_fixed_flux.interfaces.size() == interfaces.size() &&
              result.corrected.interfaces.size() == interfaces.size(),
          "Inconsistent corrected surface response entries!");
      for (const auto &[interface, raw] : result.raw.interfaces)
      {
        const auto fixed_trace = result.postprocessed_fixed_trace.interfaces.at(interface);
        const auto fixed_flux = result.postprocessed_fixed_flux.interfaces.at(interface);
        const auto corrected = result.corrected.interfaces.at(interface);
        const double p_raw = raw.energy / result.raw.domain;
        const double p_fixed_trace =
            fixed_trace.energy / result.postprocessed_fixed_trace.domain;
        const double p_fixed_flux =
            fixed_flux.energy / result.postprocessed_fixed_flux.domain;
        const double p_corrected = corrected.energy / result.corrected.domain;
        auto Quality = [](double participation, double loss_tangent)
        {
          return participation == 0.0 || loss_tangent == 0.0
                     ? mfem::infinity()
                     : 1.0 / (participation * loss_tangent);
        };
        output.table[fmt::format("energy_raw_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(raw.energy);
        output.table[fmt::format("participation_raw_{}", interface)] << p_raw;
        output.table[fmt::format("quality_raw_{}", interface)]
            << Quality(p_raw, raw.loss_tangent);
        output.table[fmt::format("energy_postprocessed_fixed_trace_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(fixed_trace.energy);
        output.table[fmt::format("participation_postprocessed_fixed_trace_{}", interface)]
            << p_fixed_trace;
        output.table[fmt::format("quality_postprocessed_fixed_trace_{}", interface)]
            << Quality(p_fixed_trace, fixed_trace.loss_tangent);
        output.table[fmt::format("energy_postprocessed_fixed_flux_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(fixed_flux.energy);
        output.table[fmt::format("participation_postprocessed_fixed_flux_{}", interface)]
            << p_fixed_flux;
        output.table[fmt::format("quality_postprocessed_fixed_flux_{}", interface)]
            << Quality(p_fixed_flux, fixed_flux.loss_tangent);
        output.table[fmt::format("energy_corrected_{}", interface)]
            << iodata.units.Dimensionalize<VT::ENERGY>(corrected.energy);
        output.table[fmt::format("participation_corrected_{}", interface)] << p_corrected;
        output.table[fmt::format("quality_corrected_{}", interface)]
            << Quality(p_corrected, corrected.loss_tangent);
        const auto closure_spread = result.trace_closure_spread.find(interface);
        output.table[fmt::format("trace_closure_spread_{}", interface)]
            << (closure_spread == result.trace_closure_spread.end()
                    ? 0.0
                    : closure_spread->second);
      }
    }
    output.WriteFullTableTrunc();
  }
  post_op.MeasureFinalize(indicator);
  return {indicator, laplace_op.GlobalTrueVSize()};
}

void ElectrostaticSolver::PostprocessResponseMatrix(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op, const LaplaceOperator &laplace_op,
    const Operator &Grad, const std::vector<Vector> &V, const std::vector<Vector> &D) const
{
  using LocalEnergy = SurfacePostOperator::InterfaceLocalEdgeEnergy;
  using EnergyMap = std::map<int, std::vector<LocalEnergy>>;

  const auto &sources = laplace_op.GetSources();
  MFEM_VERIFY(V.size() == sources.size(),
              "Unexpected prescribed-potential basis field count!");
  MFEM_VERIFY(D.empty() || D.size() == V.size(),
              "Unexpected recovered electric flux basis field count!");

  std::vector<int> basis_indices;
  basis_indices.reserve(sources.size());
  for (const auto &[idx, data] : sources)
  {
    basis_indices.push_back(idx);
  }

  Vector E_sum(Grad.Height()), D_sum;
  if (!D.empty())
  {
    D_sum.SetSize(D.front().Size());
    D_sum.UseDevice(true);
  }

  auto Evaluate = [&](std::size_t i, std::optional<std::size_t> j = std::nullopt)
  {
    E_sum = 0.0;
    Grad.AddMult(V[i], E_sum, -1.0);
    if (j)
    {
      Grad.AddMult(V[*j], E_sum, -1.0);
    }

    const Vector *D_ptr = nullptr;
    if (!D.empty())
    {
      D_sum = D[i];
      if (j)
      {
        D_sum.Add(1.0, D[*j]);
      }
      D_ptr = &D_sum;
    }
    return post_op.GetInterfaceLocalEdgeElectricFieldEnergies(E_sum, D_ptr);
  };

  Mpi::Print("\nAssembling localized interface response matrix for {:d} basis fields\n",
             static_cast<int>(V.size()));
  std::vector<EnergyMap> diagonal;
  diagonal.reserve(V.size());
  for (std::size_t i = 0; i < V.size(); i++)
  {
    diagonal.push_back(Evaluate(i));
  }

  TableWithCSVFile output;
  if (root)
  {
    output = TableWithCSVFile(post_dir / "surface-response-matrix.csv");
    output.table.insert(Column("interface", "interface", 0, 0, 2, ""));
    output.table.insert(Column("edge", "edge", 0, 0, 2, ""));
    output.table.insert("distance", "R (m)");
    output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    output.table.insert("Q", "Q_ij (J)");
    output.table.insert("Q_normal", "Q_ij normal (J)");
    output.table.insert("Q_tangential", "Q_ij tangential (J)");
    output.table.insert("Q_total", "Q_total_ij (J)");
    output.table.insert("Q_total_normal", "Q_total_ij normal (J)");
    output.table.insert("Q_total_tangential", "Q_total_ij tangential (J)");
    output.table.reserve(V.size() * (V.size() + 1) / 2, 11);
  }

  using VT = Units::ValueType;
  auto Append = [&](int interface, const LocalEnergy &energy, std::size_t i, std::size_t j,
                    double q, double q_normal, double q_tangential, double q_total,
                    double q_total_normal, double q_total_tangential)
  {
    if (!root)
    {
      return;
    }
    output.table["interface"] << interface;
    output.table["edge"] << energy.edge;
    output.table["distance"] << iodata.units.Dimensionalize<VT::LENGTH>(energy.distance);
    output.table["basis_i"] << basis_indices[i];
    output.table["basis_j"] << basis_indices[j];
    output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(q);
    output.table["Q_normal"] << iodata.units.Dimensionalize<VT::ENERGY>(q_normal);
    output.table["Q_tangential"] << iodata.units.Dimensionalize<VT::ENERGY>(q_tangential);
    output.table["Q_total"] << iodata.units.Dimensionalize<VT::ENERGY>(q_total);
    output.table["Q_total_normal"]
        << iodata.units.Dimensionalize<VT::ENERGY>(q_total_normal);
    output.table["Q_total_tangential"]
        << iodata.units.Dimensionalize<VT::ENERGY>(q_total_tangential);
  };

  for (std::size_t i = 0; i < V.size(); i++)
  {
    for (std::size_t j = i; j < V.size(); j++)
    {
      std::optional<EnergyMap> combined;
      if (i != j)
      {
        combined = Evaluate(i, j);
      }
      MFEM_VERIFY(diagonal[i].size() == diagonal[j].size() &&
                      (!combined || combined->size() == diagonal[i].size()),
                  "Inconsistent localized interface response data!");

      for (const auto &[interface, energy_i] : diagonal[i])
      {
        const auto it_j = diagonal[j].find(interface);
        const std::vector<LocalEnergy> *energy_sum = nullptr;
        if (combined)
        {
          const auto it_sum = combined->find(interface);
          MFEM_VERIFY(it_sum != combined->end(),
                      "Missing combined-field localized interface response entries!");
          energy_sum = &it_sum->second;
        }
        MFEM_VERIFY(it_j != diagonal[j].end() && energy_i.size() == it_j->second.size() &&
                        (!energy_sum || energy_i.size() == energy_sum->size()),
                    "Inconsistent localized interface response entries!");
        for (std::size_t entry = 0; entry < energy_i.size(); entry++)
        {
          const auto &ei = energy_i[entry];
          const auto &ej = it_j->second[entry];
          MFEM_VERIFY(ei.edge == ej.edge && ei.distance == ej.distance,
                      "Inconsistent localized interface response metadata!");
          if (i == j)
          {
            Append(interface, ei, i, j, ei.energy_inside, ei.energy_inside_polarized[0],
                   ei.energy_inside_polarized[1], ei.energy_total,
                   ei.energy_total_polarized[0], ei.energy_total_polarized[1]);
          }
          else
          {
            const auto &es = (*energy_sum)[entry];
            MFEM_VERIFY(ei.edge == es.edge && ei.distance == es.distance,
                        "Inconsistent localized interface response metadata!");
            Append(interface, ei, i, j,
                   0.5 * (es.energy_inside - ei.energy_inside - ej.energy_inside),
                   0.5 * (es.energy_inside_polarized[0] - ei.energy_inside_polarized[0] -
                          ej.energy_inside_polarized[0]),
                   0.5 * (es.energy_inside_polarized[1] - ei.energy_inside_polarized[1] -
                          ej.energy_inside_polarized[1]),
                   0.5 * (es.energy_total - ei.energy_total - ej.energy_total),
                   0.5 * (es.energy_total_polarized[0] - ei.energy_total_polarized[0] -
                          ej.energy_total_polarized[0]),
                   0.5 * (es.energy_total_polarized[1] - ei.energy_total_polarized[1] -
                          ej.energy_total_polarized[1]));
          }
        }
      }
    }
  }

  if (root)
  {
    output.WriteFullTableTrunc();
  }

  Mpi::Print("Assembling domain-energy response matrix\n");
  TableWithCSVFile domain_output;
  if (root)
  {
    domain_output = TableWithCSVFile(post_dir / "domain-response-matrix.csv");
    domain_output.table.insert(Column("basis_i", "basis_i", 0, 0, 2, ""));
    domain_output.table.insert(Column("basis_j", "basis_j", 0, 0, 2, ""));
    domain_output.table.insert("Q", "Q_ij (J)");
    domain_output.table.reserve(V.size() * (V.size() + 1) / 2, 3);
  }

  auto &V_gf = post_op.GetVGridFunction().Real();
  auto &D_gf = post_op.GetDomainPostOp().D;
  for (std::size_t i = 0; i < V.size(); i++)
  {
    V_gf.SetFromTrueDofs(V[i]);
    post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    for (std::size_t j = i; j < V.size(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      const double q = 0.5 * linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);
      if (root)
      {
        domain_output.table["basis_i"] << basis_indices[i];
        domain_output.table["basis_j"] << basis_indices[j];
        domain_output.table["Q"] << iodata.units.Dimensionalize<VT::ENERGY>(q);
      }
    }
  }
  if (root)
  {
    domain_output.WriteFullTableTrunc();
  }
}

void ElectrostaticSolver::PostprocessTerminals(
    PostOperator<ProblemType::ELECTROSTATIC> &post_op,
    const std::map<int, mfem::Array<int>> &terminal_sources,
    const std::vector<Vector> &V) const
{
  // Postprocess the Maxwell capacitance matrix. See p. 97 of the COMSOL AC/DC Module manual
  // for the associated formulas based on the electric field energy based on a unit voltage
  // excitation for each terminal. Alternatively, we could compute the resulting terminal
  // charges from the prescribed voltage to get C directly as:
  //         Q_i = ∫ ρ dV = ∫ ∇ ⋅ (ε E) dV = ∫ (ε E) ⋅ n dS
  // and C_ij = Q_i/V_j. The energy formulation avoids having to locally integrate E = -∇V.
  mfem::DenseMatrix C(V.size()), Cm(V.size());
  for (int i = 0; i < C.Height(); i++)
  {
    // Diagonal: Cᵢᵢ = 2 Uₑ(Vᵢ) / Vᵢ² = (Vᵢᵀ K Vᵢ) / Vᵢ² (with ∀i, Vᵢ = 1)
    auto &V_gf = post_op.GetVGridFunction().Real();
    auto &D_gf = post_op.GetDomainPostOp().D;
    V_gf.SetFromTrueDofs(V[i]);
    post_op.GetDomainPostOp().M_elec->Mult(V_gf, D_gf);
    C(i, i) = Cm(i, i) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);

    // Off-diagonals: Cᵢⱼ = Uₑ(Vᵢ + Vⱼ) / (Vᵢ Vⱼ) - 1/2 (Vᵢ/Vⱼ Cᵢᵢ + Vⱼ/Vᵢ Cⱼⱼ)
    //                    = (Vⱼᵀ K Vᵢ) / (Vᵢ Vⱼ)
    for (int j = i + 1; j < C.Width(); j++)
    {
      V_gf.SetFromTrueDofs(V[j]);
      C(i, j) = linalg::Dot<Vector>(post_op.GetComm(), V_gf, D_gf);
      Cm(i, j) = -C(i, j);
      Cm(i, i) -= Cm(i, j);
    }

    // Copy lower triangle from already computed upper triangle.
    for (int j = 0; j < i; j++)
    {
      C(i, j) = C(j, i);
      Cm(i, j) = Cm(j, i);
      Cm(i, i) -= Cm(i, j);
    }
  }
  mfem::DenseMatrix Cinv(C);
  Cinv.Invert();  // In-place, uses LAPACK (when available) and should be cheap

  // Only root writes to disk (every process has full matrices).
  if (!root)
  {
    return;
  }
  using VT = Units::ValueType;

  // Write capacitance matrix data.
  auto PrintMatrix = [&terminal_sources, this](const std::string &file,
                                               const std::string &name,
                                               const std::string &unit,
                                               const mfem::DenseMatrix &mat, double scale)
  {
    TableWithCSVFile output(post_dir / file);
    output.table.insert(Column("i", "i", 0, 0, 2, ""));
    int j = 0;
    for (const auto &[idx2, data2] : terminal_sources)
    {
      output.table.insert(fmt::format("i2{}", idx2),
                          fmt::format("{}[i][{}] {}", name, idx2, unit));
      // Use the fact that iterator over i and j is the same span.
      output.table["i"] << idx2;

      auto &col = output.table[fmt::format("i2{}", idx2)];
      for (std::size_t i = 0; i < terminal_sources.size(); i++)
      {
        col << mat(i, j) * scale;
      }
      j++;
    }
    output.WriteFullTableTrunc();
  };
  const double F = iodata.units.Dimensionalize<VT::CAPACITANCE>(1.0);
  PrintMatrix("terminal-C.csv", "C", "(F)", C, F);
  PrintMatrix("terminal-Cinv.csv", "C⁻¹", "(1/F)", Cinv, 1.0 / F);
  PrintMatrix("terminal-Cm.csv", "C_m", "(F)", Cm, F);

  // Also write out a file with terminal voltage excitations.
  {
    TableWithCSVFile terminal_V(post_dir / "terminal-V.csv");
    terminal_V.table.insert(Column("i", "i", 0, 0, 2, ""));
    terminal_V.table.insert("Vinc", "V_inc[i] (V)");
    for (const auto &[idx, data] : terminal_sources)
    {
      terminal_V.table["i"] << double(idx);
      terminal_V.table["Vinc"] << iodata.units.Dimensionalize<VT::VOLTAGE>(1.0);
    }
    terminal_V.WriteFullTableTrunc();
  }
}

}  // namespace palace
