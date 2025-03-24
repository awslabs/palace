#include "postoperatorcsv.hpp"

#include <mfem.hpp>

#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"

namespace palace
{

namespace
{

// TODO(C++20): Do constexpr with string
std::string DimLabel(int i)
{
  switch (i)
  {
    // Note: Zero-based indexing here
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      return fmt::format("d{}", i);
  }
}

// TODO(C++20): Do constexpr with string
std::string LabelIndexCol(const config::ProblemData::Type solver_t)
{
  switch (solver_t)
  {
    case config::ProblemData::Type::DRIVEN:
      return "f (GHz)";
    case config::ProblemData::Type::EIGENMODE:
      return "m";
    case config::ProblemData::Type::ELECTROSTATIC:
    case config::ProblemData::Type::MAGNETOSTATIC:
      return "i";
    case config::ProblemData::Type::TRANSIENT:
      return "t (ns)";
    default:
      return "unkown";
  }
}

}  // namespace

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeDomainE()
{
  using fmt::format;
  domain_E = TableWithCSVFile(post_op->post_dir / "domain-E.csv");
  domain_E->table.reserve(nr_expected_measurement_rows,
                          4 * (1 + post_op->dom_post_op.M_i.size()));
  domain_E->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));

  domain_E->table.insert("Ee", "E_elec (J)");
  domain_E->table.insert("Em", "E_mag (J)");
  domain_E->table.insert("Ec", "E_cap (J)");
  domain_E->table.insert("Ei", "E_ind (J)");

  for (const auto &[idx, data] : post_op->dom_post_op.M_i)
  {
    domain_E->table.insert(format("Ee{}", idx), format("E_elec[{}] (J)", idx));
    domain_E->table.insert(format("pe{}", idx), format("p_elec[{}]", idx));
    domain_E->table.insert(format("Em{}", idx), format("E_mag[{}] (J)", idx));
    domain_E->table.insert(format("pm{}", idx), format("p_mag[{}]", idx));
  }
  domain_E->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintDomainE()
{
  if (!domain_E)  // trivial check: always written and we are always on root
  {
    return;
  }
  using fmt::format;
  domain_E->table["idx"] << current_idx_value_dimensionful;
  domain_E->table["Ee"] << post_op->measurement_cache.domain_E_field_energy_all;
  domain_E->table["Em"] << post_op->measurement_cache.domain_H_field_energy_all;
  domain_E->table["Ec"] << post_op->measurement_cache.lumped_port_capacitor_energy;
  domain_E->table["Ei"] << post_op->measurement_cache.lumped_port_inductor_energy;
  for (const auto &data : post_op->measurement_cache.domain_E_field_energy_i)
  {
    domain_E->table[format("Ee{}", data.idx)] << data.energy;
    domain_E->table[format("pe{}", data.idx)] << data.participation_ratio;
  }
  for (const auto &data : post_op->measurement_cache.domain_H_field_energy_i)
  {
    domain_E->table[format("Em{}", data.idx)] << data.energy;
    domain_E->table[format("pm{}", data.idx)] << data.participation_ratio;
  }
  domain_E->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceF()
{
  if (!(post_op->surf_post_op.flux_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_F = TableWithCSVFile(post_op->post_dir / "surface-F.csv");
  surface_F->table.reserve(nr_expected_measurement_rows,
                           1 + (HasComplexGridFunction<solver_t>() ? 2 : 1) *
                                   post_op->surf_post_op.flux_surfs.size());
  surface_F->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));
  for (const auto &[idx, data] : post_op->surf_post_op.flux_surfs)
  {
    switch (data.type)
    {
      case SurfaceFluxType::ELECTRIC:
        if (HasComplexGridFunction<solver_t>())
        {
          surface_F->table.insert(format("F_{}_re", idx),
                                  format("Re{{Φ_elec[{}]}} (C)", idx));
          surface_F->table.insert(format("F_{}_im", idx),
                                  format("Im{{Φ_elec[{}]}} (C)", idx));
        }
        else
        {
          surface_F->table.insert(format("F_{}_re", idx), format("Φ_elec[{}] (C)", idx));
        }
        break;
      case SurfaceFluxType::MAGNETIC:
        if (HasComplexGridFunction<solver_t>())
        {
          surface_F->table.insert(format("F_{}_re", idx),
                                  format("Re{{Φ_mag[{}]}} (Wb)", idx));
          surface_F->table.insert(format("F_{}_im", idx),
                                  format("Im{{Φ_mag[{}]}} (Wb)", idx));
        }
        else
        {
          surface_F->table.insert(format("F_{}_re", idx), format("Φ_mag[{}] (Wb)", idx));
        }
        break;
      case SurfaceFluxType::POWER:
        surface_F->table.insert(format("F_{}_re", idx), format("Φ_pow[{}] (W)", idx));
        break;
    }
  }
  surface_F->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceF()
{
  if (!surface_F)
  {
    return;
  }
  using fmt::format;
  surface_F->table["idx"] << current_idx_value_dimensionful;
  for (const auto &flux_data : post_op->measurement_cache.surface_flux_i)
  {
    surface_F->table[format("F_{}_re", flux_data.idx)] << flux_data.Phi.real();
    if (HasComplexGridFunction<solver_t>() &&
        (flux_data.type == SurfaceFluxType::ELECTRIC ||
         flux_data.type == SurfaceFluxType::MAGNETIC))
    {
      surface_F->table[format("F_{}_im", flux_data.idx)] << flux_data.Phi.imag();
    }
  }
  surface_F->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceQ()
{
  if (!(post_op->surf_post_op.eps_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_Q = TableWithCSVFile(post_op->post_dir / "surface-Q.csv");
  surface_Q->table.reserve(nr_expected_measurement_rows,
                           1 + 2 * post_op->surf_post_op.eps_surfs.size());
  surface_Q->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));

  for (const auto &[idx, data] : post_op->surf_post_op.eps_surfs)
  {
    surface_Q->table.insert(format("p_{}", idx), format("p_surf[{}]", idx));
    surface_Q->table.insert(format("Q_{}", idx), format("Q_surf[{}]", idx));
  }
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceQ()
{
  if (!surface_Q)
  {
    return;
  }
  using fmt::format;
  surface_Q->table["idx"] << current_idx_value_dimensionful;
  for (const auto &eps_data : post_op->measurement_cache.interface_eps_i)
  {
    surface_Q->table[format("p_{}", eps_data.idx)] << eps_data.energy_participation;
    surface_Q->table[format("Q_{}", eps_data.idx)] << eps_data.quality_factor;
  }
  surface_Q->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeProbeE()
{
  if (!(post_op->interp_op.GetProbes().size() > 0) || !HasEGridFunction<solver_t>())
  {
    return;
  }
  using fmt::format;
  probe_E = TableWithCSVFile(post_op->post_dir / "probe-E.csv");
  auto v_dim = post_op->interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;

  probe_E->table.reserve(nr_expected_measurement_rows,
                         scale_col * post_op->interp_op.GetProbes().size());
  probe_E->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));

  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        probe_E->table.insert(format("E{}_{}_re", idx, i_dim),
                              format("Re{{E_{}[{}]}} (V/m)", DimLabel(i_dim), idx));
        probe_E->table.insert(format("E{}_{}_im", idx, i_dim),
                              format("Im{{E_{}[{}]}} (V/m)", DimLabel(i_dim), idx));
      }
      else
      {
        probe_E->table.insert(format("E{}_{}_re", idx, i_dim),
                              format("E_{}[{}] (V/m)", DimLabel(i_dim), idx));
      }
    }
  }
  probe_E->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintProbeE()
{
  if (!probe_E)
  {
    return;
  }
  using fmt::format;
  auto v_dim = post_op->interp_op.GetVDim();
  auto probe_field = post_op->measurement_cache.probe_E_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  probe_E->table["idx"] << current_idx_value_dimensionful;
  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_E->table[format("E{}_{}_re", idx, i_dim)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_E->table[format("E{}_{}_im", idx, i_dim)] << val.imag();
      }
    }
    i++;
  }
  probe_E->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeProbeB()
{
  if (!(post_op->interp_op.GetProbes().size() > 0) || !HasBGridFunction<solver_t>())
  {
    return;
  }
  using fmt::format;
  probe_B = TableWithCSVFile(post_op->post_dir / "probe-B.csv");
  auto v_dim = post_op->interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;

  probe_B->table.reserve(nr_expected_measurement_rows,
                         scale_col * post_op->interp_op.GetProbes().size());
  probe_B->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));

  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      if (HasComplexGridFunction<solver_t>())
      {
        probe_B->table.insert(format("B{}_{}_re", idx, i_dim),
                              format("Re{{B_{}[{}]}} (Wb/m²)", DimLabel(i_dim), idx));
        probe_B->table.insert(format("B{}_{}_im", idx, i_dim),
                              format("Im{{B_{}[{}]}} (Wb/m²)", DimLabel(i_dim), idx));
      }
      else
      {
        probe_B->table.insert(format("B{}_{}_re", idx, i_dim),
                              format("B_{}[{}] (Wb/m²)", DimLabel(i_dim), idx));
      }
    }
  }
  probe_B->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintProbeB()
{
  if (!probe_B)
  {
    return;
  }
  using fmt::format;

  auto v_dim = post_op->interp_op.GetVDim();
  auto probe_field = post_op->measurement_cache.probe_B_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  probe_B->table["idx"] << current_idx_value_dimensionful;
  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_B->table[format("B{}_{}_re", idx, i_dim)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_B->table[format("B{}_{}_im", idx, i_dim)] << val.imag();
      }
    }
    i++;
  }
  probe_B->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeSurfaceI()
    -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN ||
                            U == config::ProblemData::Type::TRANSIENT,
                        void>
{
  if (!(post_op->fem_op->GetSurfaceCurrentOp().Size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_I = TableWithCSVFile(post_op->post_dir / "surface-I.csv");
  const auto &surf_j_op = post_op->fem_op->GetSurfaceCurrentOp();
  surface_I->table.reserve(nr_expected_measurement_rows, surf_j_op.Size());
  surface_I->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));
  for (const auto &[idx, data] : surf_j_op)
  {
    surface_I->table.insert(format("I_{}", idx), format("I_inc[{}] (A)", idx));
  }
  surface_I->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintSurfaceI()
    -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN ||
                            U == config::ProblemData::Type::TRANSIENT,
                        void>
{
  if (!surface_I)
  {
    return;
  }
  using fmt::format;
  surface_I->table["idx"] << current_idx_value_dimensionful;
  for (const auto &[idx, data] : post_op->fem_op->GetSurfaceCurrentOp())
  {
    auto I_inc_raw =
        data.GetExcitationCurrent() * post_op->measurement_cache.Jcoeff_excitation;
    auto I_inc =
        post_op->units.template Dimensionalize<Units::ValueType::CURRENT>(I_inc_raw);
    surface_I->table[format("I_{}", idx)] << I_inc;
  }
  surface_I->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializePortVI()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE ||
                            U == config::ProblemData::Type::DRIVEN ||
                            U == config::ProblemData::Type::TRANSIENT,
                        void>
{
  if (!(post_op->fem_op->GetLumpedPortOp().Size() > 0))
  {
    return;
  }
  using fmt::format;
  // Currently only works for lumped ports
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  port_V = TableWithCSVFile(post_op->post_dir / "port-V.csv");
  port_I = TableWithCSVFile(post_op->post_dir / "port-I.csv");

  port_V->table.reserve(nr_expected_measurement_rows, 1 + lumped_port_op.Size());
  port_I->table.reserve(nr_expected_measurement_rows, 1 + lumped_port_op.Size());

  port_V->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));
  port_I->table.insert(Column("idx", LabelIndexCol(solver_t), 0, {}, {}, ""));

  for (const auto &[idx, data] : lumped_port_op)
  {
    if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                  solver_t == config::ProblemData::Type::TRANSIENT)
    {
      if (data.excitation)
      {
        // Incident voltage is currently always real
        port_V->table.insert(format("inc{}", idx), format("V_inc[{}] (V)", idx));
        port_I->table.insert(format("inc{}", idx), format("I_inc[{}] (A)", idx));
      }
    }
    if constexpr (HasComplexGridFunction<solver_t>())
    {
      port_V->table.insert(format("re{}", idx), format("Re{{V[{}]}} (V)", idx));
      port_V->table.insert(format("im{}", idx), format("Im{{V[{}]}} (V)", idx));

      port_I->table.insert(format("re{}", idx), format("Re{{I[{}]}} (A)", idx));
      port_I->table.insert(format("im{}", idx), format("Im{{I[{}]}} (A)", idx));
    }
    else
    {
      port_V->table.insert(format("re{}", idx), format("V[{}] (V)", idx));
      port_I->table.insert(format("re{}", idx), format("I[{}] (A)", idx));
    }
  }
  port_V->WriteFullTableTrunc();
  port_I->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintPortVI()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE ||
                            U == config::ProblemData::Type::DRIVEN ||
                            U == config::ProblemData::Type::TRANSIENT,
                        void>
{
  if (!port_V)  // no need to recheck port_I
  {
    return;
  }
  using fmt::format;
  // Currently only works for lumped ports
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).

  port_V->table["idx"] << current_idx_value_dimensionful;
  port_I->table["idx"] << current_idx_value_dimensionful;

  auto unit_V = post_op->units.template GetScaleFactor<Units::ValueType::VOLTAGE>();
  auto unit_A = post_op->units.template GetScaleFactor<Units::ValueType::CURRENT>();

  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    for (const auto &[idx, data] : lumped_port_op)
    {
      if (data.excitation)
      {
        auto Jcoeff = post_op->measurement_cache.Jcoeff_excitation;
        double V_inc = data.GetExcitationVoltage() * Jcoeff;
        double I_inc = (std::abs(V_inc) > 0.0)
                           ? data.GetExcitationPower() * Jcoeff * Jcoeff / V_inc
                           : 0.0;

        port_V->table[format("inc{}", idx)] << V_inc * unit_V;
        port_I->table[format("inc{}", idx)] << I_inc * unit_A;
      }
    }
  }

  for (const auto &[idx, data] : post_op->measurement_cache.lumped_port_vi)
  {
    port_V->table[fmt::format("re{}", idx)] << data.V.real() * unit_V;
    port_I->table[fmt::format("re{}", idx)] << data.I.real() * unit_A;

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      port_V->table[fmt::format("im{}", idx)] << data.V.imag() * unit_V;
      port_I->table[fmt::format("im{}", idx)] << data.I.imag() * unit_A;
    }
  }
  port_V->WriteFullTableTrunc();
  port_I->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializePortS()
    -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, void>
{
  if (!((post_op->fem_op->GetLumpedPortOp().Size() > 0) xor
        (post_op->fem_op->GetWavePortOp().Size() > 0)))
  {
    return;
  }
  // Get source index of single excitation from space_op
  driven_source_index = -1;
  // Get excitation index as is currently done: if -1 then no excitation
  // Already ensured that one of lumped or wave ports are empty
  for (const auto &[idx, data] : post_op->fem_op->GetLumpedPortOp())
  {
    if (data.excitation)
    {
      driven_source_index = idx;
    }
  }
  for (const auto &[idx, data] : post_op->fem_op->GetWavePortOp())
  {
    if (data.excitation)
    {
      driven_source_index = idx;
    }
  }
  if (!(driven_source_index > 0))
  {
    return;
  }
  using fmt::format;
  port_S = TableWithCSVFile(post_op->post_dir / "port-S.csv");

  auto nr_ports = std::max(post_op->fem_op->GetLumpedPortOp().Size(),
                           post_op->fem_op->GetWavePortOp().Size());
  port_S->table.reserve(nr_expected_measurement_rows, nr_ports);
  port_S->table.insert(Column("idx", "f (GHz)", 0, {}, {}, ""));

  // Already ensured that one of lumped or wave ports are empty
  for (const auto &[o_idx, data] : post_op->fem_op->GetLumpedPortOp())
  {
    port_S->table.insert(format("abs_{}_{}", o_idx, driven_source_index),
                         format("|S[{}][{}]| (dB)", o_idx, driven_source_index));
    port_S->table.insert(format("arg_{}_{}", o_idx, driven_source_index),
                         format("arg(S[{}][{}]) (deg.)", o_idx, driven_source_index));
  }
  for (const auto &[o_idx, data] : post_op->fem_op->GetWavePortOp())
  {
    port_S->table.insert(format("abs_{}_{}", o_idx, driven_source_index),
                         format("|S[{}][{}]| (dB)", o_idx, driven_source_index));
    port_S->table.insert(format("arg_{}_{}", o_idx, driven_source_index),
                         format("arg(S[{}][{}]) (deg.)", o_idx, driven_source_index));
  }
  port_S->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintPortS()
    -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, void>
{
  if (!port_S)
  {
    return;
  }
  using fmt::format;
  port_S->table["idx"] << current_idx_value_dimensionful;
  for (const auto &[idx, data] : post_op->measurement_cache.lumped_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, driven_source_index)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, driven_source_index)] << data.arg_S_ij;
  }
  for (const auto &[idx, data] : post_op->measurement_cache.wave_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, driven_source_index)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, driven_source_index)] << data.arg_S_ij;
  }
  port_S->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeEig()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  using fmt::format;
  eig = TableWithCSVFile(post_op->post_dir / "eig.csv");
  eig->table.reserve(nr_expected_measurement_rows, 6);
  eig->table.insert(Column("idx", "m", 0, {}, {}, ""));
  eig->table.insert("f_re", "Re{f} (GHz)");
  eig->table.insert("f_im", "Im{f} (GHz)");
  eig->table.insert("q", "Q");
  eig->table.insert("err_back", "Error (Bkwd.)");
  eig->table.insert("err_abs", "Error (Abs.)");
  eig->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintEig()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  if (!eig)  // trivial check
  {
    return;
  }
  eig->table["idx"] << current_idx_value_dimensionful;
  eig->table["f_re"] << post_op->measurement_cache.freq.real();
  eig->table["f_im"] << post_op->measurement_cache.freq.imag();
  eig->table["q"] << post_op->measurement_cache.eigenmode_Q;
  eig->table["err_back"] << post_op->measurement_cache.error_bkwd;
  eig->table["err_abs"] << post_op->measurement_cache.error_abs;
  eig->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeEigPortEPR()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  // TODO(C++20): Make this a filterd iterator in LumpedPortOp
  for (const auto &[idx, data] : post_op->fem_op->GetLumpedPortOp())
  {
    if (std::abs(data.L) > 0.0)
    {
      ports_with_L.push_back(idx);
    }
  }
  if (ports_with_L.empty())
  {
    return;
  }
  using fmt::format;
  port_EPR = TableWithCSVFile(post_op->post_dir / "port-EPR.csv");
  port_EPR->table.reserve(nr_expected_measurement_rows, 1 + ports_with_L.size());
  port_EPR->table.insert(Column("idx", "m", 0, {}, {}, ""));
  for (const auto idx : ports_with_L)
  {
    port_EPR->table.insert(format("p_{}", idx), format("p[{}]", idx));
  }
  port_EPR->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintEigPortEPR()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  if (!port_EPR)
  {
    return;
  }
  using fmt::format;
  port_EPR->table["idx"] << current_idx_value_dimensionful;
  for (const auto idx : ports_with_L)
  {
    auto vi = post_op->measurement_cache.lumped_port_vi.at(idx);
    port_EPR->table[format("p_{}", idx)] << vi.inductive_energy_participation;
  }
  port_EPR->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeEigPortQ()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  // TODO(C++20): Make this a filtered iterator in LumpedPortOp
  for (const auto &[idx, data] : post_op->fem_op->GetLumpedPortOp())
  {
    if (std::abs(data.R) > 0.0)
    {
      ports_with_R.push_back(idx);
    }
  }
  if (ports_with_R.empty())
  {
    return;
  }
  using fmt::format;
  port_Q = TableWithCSVFile(post_op->post_dir / "port-Q.csv");
  port_Q->table.reserve(nr_expected_measurement_rows, 1 + ports_with_R.size());
  port_Q->table.insert(Column("idx", "m", 0, {}, {}, ""));
  for (const auto idx : ports_with_R)
  {
    port_Q->table.insert(format("Ql_{}", idx), format("Q_ext[{}]", idx));
    port_Q->table.insert(format("Kl_{}", idx), format("κ_ext[{}] (GHz)", idx));
  }
  port_Q->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::PrintEigPortQ()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
{
  if (!port_Q)
  {
    return;
  }
  using fmt::format;
  port_Q->table["idx"] << current_idx_value_dimensionful;
  for (const auto idx : ports_with_R)
  {
    auto vi = post_op->measurement_cache.lumped_port_vi.at(idx);
    port_Q->table[format("Ql_{}", idx)] << vi.quality_factor;
    port_Q->table[format("Kl_{}", idx)] << vi.mode_port_kappa;
  }
  port_Q->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintErrorIndicator(
    const ErrorIndicator::SummaryStatistics &indicator_stats)
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }

  TableWithCSVFile error_indicator(post_op->post_dir / "error-indicators.csv");
  error_indicator.table.reserve(1, 4);

  error_indicator.table.insert(Column("norm", "Norm") << indicator_stats.norm);
  error_indicator.table.insert(Column("min", "Minimum") << indicator_stats.min);
  error_indicator.table.insert(Column("max", "Maximum") << indicator_stats.max);
  error_indicator.table.insert(Column("mean", "Mean") << indicator_stats.mean);

  error_indicator.WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeCSVDataCollection()
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }
  InitializeDomainE();
  InitializeSurfaceF();
  InitializeSurfaceQ();
#if defined(MFEM_USE_GSLIB)
  InitializeProbeE();
  InitializeProbeB();
#endif
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    InitializeSurfaceI();
  }
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::EIGENMODE ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    InitializePortVI();
  }
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN)
  {
    InitializePortS();
  }
  if constexpr (solver_t == config::ProblemData::Type::EIGENMODE)
  {
    InitializeEig();
    InitializeEigPortEPR();
    InitializeEigPortQ();
  }
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintAllCSVData(double idx_value_dimensionful, int step,
                                                int column_block)
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }
  current_idx_value_dimensionful = idx_value_dimensionful;
  current_idx_row = step;
  current_column_block = column_block;

  PrintDomainE();
  PrintSurfaceF();
  PrintSurfaceQ();
#if defined(MFEM_USE_GSLIB)
  PrintProbeE();
  PrintProbeB();
#endif
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    PrintSurfaceI();
  }

  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::EIGENMODE ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    PrintPortVI();
  }
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN)
  {
    PrintPortS();
  }
  if constexpr (solver_t == config::ProblemData::Type::EIGENMODE)
  {
    PrintEig();
    PrintEigPortEPR();
    PrintEigPortQ();
  }
}

// Explict template instantiation

template class PostOperatorCSV<config::ProblemData::Type::DRIVEN>;
template class PostOperatorCSV<config::ProblemData::Type::EIGENMODE>;
template class PostOperatorCSV<config::ProblemData::Type::ELECTROSTATIC>;
template class PostOperatorCSV<config::ProblemData::Type::MAGNETOSTATIC>;
template class PostOperatorCSV<config::ProblemData::Type::TRANSIENT>;

}  // namespace palace
