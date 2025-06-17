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

// TODO(C++20): Do constexpr with string.
std::string DimLabel(int i)
{
  switch (i)
  {
    // Note: Zero-based indexing here.
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

// TODO(C++20): Do constexpr with string.
std::string LabelIndexCol(const ProblemType solver_t)
{
  switch (solver_t)
  {
    case ProblemType::DRIVEN:
      return "f (GHz)";
    case ProblemType::EIGENMODE:
      return "m";
    case ProblemType::ELECTROSTATIC:
    case ProblemType::MAGNETOSTATIC:
      return "i";
    case ProblemType::TRANSIENT:
      return "t (ns)";
    default:
      return "unkown";
  }
}
int PrecIndexCol(const ProblemType solver_t)
{
  switch (solver_t)
  {
    case ProblemType::DRIVEN:
    case ProblemType::TRANSIENT:
      return 8;
    case ProblemType::EIGENMODE:
    case ProblemType::ELECTROSTATIC:
    case ProblemType::MAGNETOSTATIC:
      return 2;
    default:
      return 8;
  }
}

// Index checking when adding to a new excitation block: When adding data to data_col with
// index idx, checks that idx matches what is already written in the corresponding row of
// idx_col. Adds a new idx row to idx_col if needed.
void CheckAppendIndex(Column &idx_col, double idx_value, size_t m_idx_row)
{
  if (m_idx_row == idx_col.n_rows())
  {
    idx_col << idx_value;
  }
  else
  {
    auto current_idx = idx_col.data.at(m_idx_row);
    MFEM_VERIFY(idx_value == current_idx,
                fmt::format("Writing data table at incorrect index. Data has index {} "
                            "while table is at {}",
                            idx_value, current_idx));
  }
}

}  // namespace

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeDomainE()
{
  using fmt::format;
  domain_E = TableWithCSVFile(post_op->post_dir / "domain-E.csv");
  auto nr_expected_measurement_cols =
      1 + excitation_idx_all.size() * 4 * (1 + post_op->dom_post_op.M_i.size());
  domain_E->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  domain_E->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));

  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);

    domain_E->table.insert(format("Ee_{}", ex_idx), format("E_elec{} (J)", ex_label));
    domain_E->table.insert(format("Em_{}", ex_idx), format("E_mag{} (J)", ex_label));
    domain_E->table.insert(format("Ec_{}", ex_idx), format("E_cap{} (J)", ex_label));
    domain_E->table.insert(format("Ei_{}", ex_idx), format("E_ind{} (J)", ex_label));

    for (const auto &[idx, data] : post_op->dom_post_op.M_i)
    {
      domain_E->table.insert(format("Ee_{}_{}", idx, ex_idx),
                             format("E_elec[{}]{} (J)", idx, ex_label));
      domain_E->table.insert(format("pe_{}_{}", idx, ex_idx),
                             format("p_elec[{}]{}", idx, ex_label));
      domain_E->table.insert(format("Em_{}_{}", idx, ex_idx),
                             format("E_mag[{}]{} (J)", idx, ex_label));
      domain_E->table.insert(format("pm_{}_{}", idx, ex_idx),
                             format("p_mag[{}]{}", idx, ex_label));
    }
  }
  domain_E->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintDomainE()
{
  if (!domain_E)  // trivial check: always written and we are always on root
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(domain_E->table["idx"], m_idx_value, m_idx_row);
  domain_E->table[format("Ee_{}", m_ex_idx)] << MCache().domain_E_field_energy_all;
  domain_E->table[format("Em_{}", m_ex_idx)] << MCache().domain_H_field_energy_all;
  domain_E->table[format("Ec_{}", m_ex_idx)] << MCache().lumped_port_capacitor_energy;
  domain_E->table[format("Ei_{}", m_ex_idx)] << MCache().lumped_port_inductor_energy;
  for (const auto &data : MCache().domain_E_field_energy_i)
  {
    domain_E->table[format("Ee_{}_{}", data.idx, m_ex_idx)] << data.energy;
    domain_E->table[format("pe_{}_{}", data.idx, m_ex_idx)] << data.participation_ratio;
  }
  for (const auto &data : MCache().domain_H_field_energy_i)
  {
    domain_E->table[format("Em_{}_{}", data.idx, m_ex_idx)] << data.energy;
    domain_E->table[format("pm_{}_{}", data.idx, m_ex_idx)] << data.participation_ratio;
  }
  domain_E->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceF()
{
  if (!(post_op->surf_post_op.flux_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_F = TableWithCSVFile(post_op->post_dir / "surface-F.csv");
  auto nr_expected_measurement_cols = 1 + excitation_idx_all.size() *
                                              (HasComplexGridFunction<solver_t>() ? 2 : 1) *
                                              post_op->surf_post_op.flux_surfs.size();
  surface_F->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  surface_F->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : post_op->surf_post_op.flux_surfs)
    {
      switch (data.type)
      {
        case SurfaceFluxType::ELECTRIC:
          if (HasComplexGridFunction<solver_t>())
          {
            surface_F->table.insert(format("F_{}_{}_re", idx, ex_idx),
                                    format("Re{{Φ_elec[{}]{}}} (C)", idx, ex_label));
            surface_F->table.insert(format("F_{}_{}_im", idx, ex_idx),
                                    format("Im{{Φ_elec[{}]{}}} (C)", idx, ex_label));
          }
          else
          {
            surface_F->table.insert(format("F_{}_{}_re", idx, ex_idx),
                                    format("Φ_elec[{}]{} (C)", idx, ex_label));
          }
          break;
        case SurfaceFluxType::MAGNETIC:
          if (HasComplexGridFunction<solver_t>())
          {
            surface_F->table.insert(format("F_{}_{}_re", idx, ex_idx),
                                    format("Re{{Φ_mag[{}]{}}} (Wb)", idx, ex_label));
            surface_F->table.insert(format("F_{}_{}_im", idx, ex_idx),
                                    format("Im{{Φ_mag[{}]{}}} (Wb)", idx, ex_label));
          }
          else
          {
            surface_F->table.insert(format("F_{}_{}_re", idx, ex_idx),
                                    format("Φ_mag[{}]{} (Wb)", idx, ex_label));
          }
          break;
        case SurfaceFluxType::POWER:
          surface_F->table.insert(format("F_{}_{}_re", idx, ex_idx),
                                  format("Φ_pow[{}]{} (W)", idx, ex_label));
          break;
      }
    }
  }
  surface_F->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceF()
{
  if (!surface_F)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_F->table["idx"], m_idx_value, m_idx_row);
  for (const auto &data : MCache().surface_flux_i)
  {
    surface_F->table[format("F_{}_{}_re", data.idx, m_ex_idx)] << data.Phi.real();
    if (HasComplexGridFunction<solver_t>() &&
        (data.type == SurfaceFluxType::ELECTRIC || data.type == SurfaceFluxType::MAGNETIC))
    {
      surface_F->table[format("F_{}_{}_im", data.idx, m_ex_idx)] << data.Phi.imag();
    }
  }
  surface_F->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceQ()
{
  if (!(post_op->surf_post_op.eps_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_Q = TableWithCSVFile(post_op->post_dir / "surface-Q.csv");
  auto nr_expected_measurement_cols =
      1 + excitation_idx_all.size() * (2 * post_op->surf_post_op.eps_surfs.size());
  surface_Q->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  surface_Q->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : post_op->surf_post_op.eps_surfs)
    {
      surface_Q->table.insert(format("p_{}_{}", idx, ex_idx),
                              format("p_surf[{}]{}", idx, ex_label));
      surface_Q->table.insert(format("Q_{}_{}", idx, ex_idx),
                              format("Q_surf[{}]{}", idx, ex_label));
    }
  }
  surface_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceQ()
{
  if (!surface_Q)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_Q->table["idx"], m_idx_value, m_idx_row);

  for (const auto &data : MCache().interface_eps_i)
  {
    surface_Q->table[format("p_{}_{}", data.idx, m_ex_idx)] << data.energy_participation;
    surface_Q->table[format("Q_{}_{}", data.idx, m_ex_idx)] << data.quality_factor;
  }
  surface_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
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
  auto nr_expected_measurement_cols =
      1 + excitation_idx_all.size() * scale_col * post_op->interp_op.GetProbes().size();
  probe_E->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  probe_E->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : post_op->interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if constexpr (HasComplexGridFunction<solver_t>())
        {
          probe_E->table.insert(
              format("E{}_{}_{}_re", i_dim, idx, ex_idx),
              format("Re{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label));
          probe_E->table.insert(
              format("E{}_{}_{}_im", i_dim, idx, ex_idx),
              format("Im{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label));
        }
        else
        {
          probe_E->table.insert(format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                                format("E_{}[{}]{} (V/m)", DimLabel(i_dim), idx, ex_label));
        }
      }
    }
  }
  probe_E->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintProbeE()
{
  if (!probe_E)
  {
    return;
  }
  using fmt::format;
  auto v_dim = post_op->interp_op.GetVDim();
  auto probe_field = MCache().probe_E_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_E->table["idx"], m_idx_value, m_idx_row);

  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_E->table[format("E{}_{}_{}_re", i_dim, idx, m_ex_idx)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_E->table[format("E{}_{}_{}_im", i_dim, idx, m_ex_idx)] << val.imag();
      }
    }
    i++;
  }
  probe_E->WriteFullTableTrunc();
}

template <ProblemType solver_t>
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
  auto nr_expected_measurement_cols =
      1 + excitation_idx_all.size() * scale_col * post_op->interp_op.GetProbes().size();
  probe_B->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  probe_B->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : post_op->interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if (HasComplexGridFunction<solver_t>())
        {
          probe_B->table.insert(
              format("B{}_{}_{}_re", i_dim, idx, ex_idx),
              format("Re{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
          probe_B->table.insert(
              format("B{}_{}_{}_im", i_dim, idx, ex_idx),
              format("Im{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
        }
        else
        {
          probe_B->table.insert(
              format("B{}_{}_{}_re", i_dim, idx, ex_idx),
              format("B_{}[{}]{} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
        }
      }
    }
  }
  probe_B->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintProbeB()
{
  if (!probe_B)
  {
    return;
  }
  using fmt::format;

  auto v_dim = post_op->interp_op.GetVDim();
  auto probe_field = MCache().probe_B_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_B->table["idx"], m_idx_value, m_idx_row);

  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_B->table[format("B{}_{}_{}_re", i_dim, idx, m_ex_idx)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_B->table[format("B{}_{}_{}_im", i_dim, idx, m_ex_idx)] << val.imag();
      }
    }
    i++;
  }
  probe_B->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeSurfaceI()
    -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>
{
  if (!(post_op->fem_op->GetSurfaceCurrentOp().Size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_I = TableWithCSVFile(post_op->post_dir / "surface-I.csv");
  const auto &surf_j_op = post_op->fem_op->GetSurfaceCurrentOp();
  auto nr_expected_measurement_cols = 1 + excitation_idx_all.size() * surf_j_op.Size();
  surface_I->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  surface_I->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : surf_j_op)
    {
      surface_I->table.insert(format("I_{}_{}", idx, ex_idx),
                              format("I_inc[{}]{} (A)", idx, ex_label));
    }
  }
  surface_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintSurfaceI()
    -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>
{
  if (!surface_I)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_I->table["idx"], m_idx_value, m_idx_row);
  for (const auto &[idx, data] : post_op->fem_op->GetSurfaceCurrentOp())
  {
    auto I_inc_raw = data.GetExcitationCurrent() * MCache().Jcoeff_excitation;
    auto I_inc =
        post_op->units.template Dimensionalize<Units::ValueType::CURRENT>(I_inc_raw);
    surface_I->table[format("I_{}_{}", idx, m_ex_idx)] << I_inc;
  }
  surface_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializePortVI()
    -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                            U == ProblemType::TRANSIENT,
                        void>
{
  if (!(post_op->fem_op->GetLumpedPortOp().Size() > 0))
  {
    return;
  }
  using fmt::format;
  // Currently only works for lumped ports.
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  port_V = TableWithCSVFile(post_op->post_dir / "port-V.csv");
  port_I = TableWithCSVFile(post_op->post_dir / "port-I.csv");

  auto nr_expected_measurement_cols = 1 + excitation_idx_all.size() * lumped_port_op.Size();
  port_V->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  port_I->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);

  port_V->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  port_I->table.insert(
      Column("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto ex_idx : excitation_idx_all)
  {
    std::string ex_label = SingleColBlock() ? "" : format("[{}]", ex_idx);

    // Print incident signal, if solver supports excitation on ports.
    if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
    {
      auto ex_spec = post_op->fem_op->GetPortExcitations().excitations.at(ex_idx);
      for (const auto &idx : ex_spec.lumped_port)
      {
        port_V->table.insert(format("inc{}_{}", idx, ex_idx),
                             format("V_inc[{}]{} (V)", idx, ex_label));
        port_I->table.insert(format("inc{}_{}", idx, ex_idx),
                             format("I_inc[{}]{} (A)", idx, ex_label));
      }
    }
    for (const auto &[idx, data] : lumped_port_op)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        port_V->table.insert(format("re{}_{}", idx, ex_idx),
                             format("Re{{V[{}]{}}} (V)", idx, ex_label));
        port_V->table.insert(format("im{}_{}", idx, ex_idx),
                             format("Im{{V[{}]{}}} (V)", idx, ex_label));
        port_I->table.insert(format("re{}_{}", idx, ex_idx),
                             format("Re{{I[{}]{}}} (A)", idx, ex_label));
        port_I->table.insert(format("im{}_{}", idx, ex_idx),
                             format("Im{{I[{}]{}}} (A)", idx, ex_label));
      }
      else
      {
        port_V->table.insert(format("re{}_{}", idx, ex_idx),
                             format("V[{}]{} (V)", idx, ex_label));
        port_I->table.insert(format("re{}_{}", idx, ex_idx),
                             format("I[{}]{} (A)", idx, ex_label));
      }
    }
  }
  port_V->WriteFullTableTrunc();
  port_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintPortVI()
    -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                            U == ProblemType::TRANSIENT,
                        void>
{
  if (!port_V)  // no need to recheck port_I
  {
    return;
  }
  using fmt::format;
  // Currently only works for lumped ports.
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).

  CheckAppendIndex(port_V->table["idx"], m_idx_value, m_idx_row);
  CheckAppendIndex(port_I->table["idx"], m_idx_value, m_idx_row);

  auto unit_V = post_op->units.template GetScaleFactor<Units::ValueType::VOLTAGE>();
  auto unit_A = post_op->units.template GetScaleFactor<Units::ValueType::CURRENT>();

  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    for (const auto &[idx, data] : lumped_port_op)
    {
      if (data.excitation == m_ex_idx)
      {
        auto Jcoeff = MCache().Jcoeff_excitation;
        double V_inc = data.GetExcitationVoltage() * Jcoeff;
        double I_inc = (std::abs(V_inc) > 0.0)
                           ? data.GetExcitationPower() * Jcoeff * Jcoeff / V_inc
                           : 0.0;

        port_V->table[format("inc{}_{}", idx, m_ex_idx)] << V_inc * unit_V;
        port_I->table[format("inc{}_{}", idx, m_ex_idx)] << I_inc * unit_A;
      }
    }
  }

  for (const auto &[idx, data] : MCache().lumped_port_vi)
  {
    port_V->table[fmt::format("re{}_{}", idx, m_ex_idx)] << data.V.real() * unit_V;
    port_I->table[fmt::format("re{}_{}", idx, m_ex_idx)] << data.I.real() * unit_A;

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      port_V->table[fmt::format("im{}_{}", idx, m_ex_idx)] << data.V.imag() * unit_V;
      port_I->table[fmt::format("im{}_{}", idx, m_ex_idx)] << data.I.imag() * unit_A;
    }
  }
  port_V->WriteFullTableTrunc();
  port_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializePortS()
    -> std::enable_if_t<U == ProblemType::DRIVEN, void>
{
  if (!post_op->fem_op->GetPortExcitations().IsMultipleSimple() ||
      !((post_op->fem_op->GetLumpedPortOp().Size() > 0) xor
        (post_op->fem_op->GetWavePortOp().Size() > 0)))
  {
    return;
  }
  using fmt::format;
  port_S = TableWithCSVFile(post_op->post_dir / "port-S.csv");
  auto nr_ports =
      post_op->fem_op->GetLumpedPortOp().Size() + post_op->fem_op->GetWavePortOp().Size();

  auto nr_expected_measurement_cols = 1 + excitation_idx_all.size() * nr_ports;
  port_S->table.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  port_S->table.insert(Column("idx", "f (GHz)", 0, PrecIndexCol(solver_t), {}, ""));

  for (const auto ex_idx : excitation_idx_all)
  {
    // TODO(C++20): Combine identical loops with ranges + projection.
    for (const auto &[o_idx, data] : post_op->fem_op->GetLumpedPortOp())
    {
      port_S->table.insert(format("abs_{}_{}", o_idx, ex_idx),
                           format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      port_S->table.insert(format("arg_{}_{}", o_idx, ex_idx),
                           format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
    for (const auto &[o_idx, data] : post_op->fem_op->GetWavePortOp())
    {
      port_S->table.insert(format("abs_{}_{}", o_idx, ex_idx),
                           format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      port_S->table.insert(format("arg_{}_{}", o_idx, ex_idx),
                           format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
  }
  port_S->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintPortS()
    -> std::enable_if_t<U == ProblemType::DRIVEN, void>
{
  if (!port_S)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(port_S->table["idx"], m_idx_value, m_idx_row);
  for (const auto &[idx, data] : MCache().lumped_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, m_ex_idx)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, m_ex_idx)] << data.arg_S_ij;
  }
  for (const auto &[idx, data] : MCache().wave_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, m_ex_idx)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, m_ex_idx)] << data.arg_S_ij;
  }
  port_S->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEig()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  using fmt::format;
  eig = TableWithCSVFile(post_op->post_dir / "eig.csv");
  eig->table.reserve(nr_expected_measurement_rows, 6);
  eig->table.insert(Column("idx", "m", 0, PrecIndexCol(solver_t), {}, ""));
  eig->table.insert("f_re", "Re{f} (GHz)");
  eig->table.insert("f_im", "Im{f} (GHz)");
  eig->table.insert("q", "Q");
  eig->table.insert("err_back", "Error (Bkwd.)");
  eig->table.insert("err_abs", "Error (Abs.)");
  eig->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintEig()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  if (!eig)  // trivial check
  {
    return;
  }
  eig->table["idx"] << m_idx_value;
  eig->table["f_re"] << MCache().freq.real();
  eig->table["f_im"] << MCache().freq.imag();
  eig->table["q"] << MCache().eigenmode_Q;
  eig->table["err_back"] << MCache().error_bkwd;
  eig->table["err_abs"] << MCache().error_abs;
  eig->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEigPortEPR()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  // TODO(C++20): Make this a filtered iterator in LumpedPortOp.
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
  port_EPR->table.insert(Column("idx", "m", 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto idx : ports_with_L)
  {
    port_EPR->table.insert(format("p_{}", idx), format("p[{}]", idx));
  }
  port_EPR->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintEigPortEPR()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  if (!port_EPR)
  {
    return;
  }
  using fmt::format;
  port_EPR->table["idx"] << m_idx_value;
  for (const auto idx : ports_with_L)
  {
    auto vi = MCache().lumped_port_vi.at(idx);
    port_EPR->table[format("p_{}", idx)] << vi.inductive_energy_participation;
  }
  port_EPR->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEigPortQ()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  // TODO(C++20): Make this a filtered iterator in LumpedPortOp.
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
  port_Q->table.insert(Column("idx", "m", 0, PrecIndexCol(solver_t), {}, ""));
  for (const auto idx : ports_with_R)
  {
    port_Q->table.insert(format("Ql_{}", idx), format("Q_ext[{}]", idx));
    port_Q->table.insert(format("Kl_{}", idx), format("κ_ext[{}] (GHz)", idx));
  }
  port_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintEigPortQ()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  if (!port_Q)
  {
    return;
  }
  using fmt::format;
  port_Q->table["idx"] << m_idx_value;
  for (const auto idx : ports_with_R)
  {
    auto vi = MCache().lumped_port_vi.at(idx);
    port_Q->table[format("Ql_{}", idx)] << vi.quality_factor;
    port_Q->table[format("Kl_{}", idx)] << vi.mode_port_kappa;
  }
  port_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
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

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeCSVDataCollection()
{
  // Initialize multi-excitation column block index. Only driven or transient support
  // excitations; for other solvers this is default to a single idx=0.
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    auto excitation_helper = post_op->fem_op->GetPortExcitations();
    excitation_idx_all.clear();
    excitation_idx_all.reserve(excitation_helper.Size());
    std::transform(excitation_helper.begin(), excitation_helper.end(),
                   std::back_inserter(excitation_idx_all),
                   [](const auto &pair) { return pair.first; });
    // Default to the first excitation.
    m_ex_idx = excitation_idx_all.front();
  }

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
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    InitializeSurfaceI();
  }
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::TRANSIENT)
  {
    InitializePortVI();
  }
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    InitializePortS();
  }
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    InitializeEig();
    InitializeEigPortEPR();
    InitializeEigPortQ();
  }
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintAllCSVData(double idx_value_dimensionful, int step)
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }
  m_idx_value = idx_value_dimensionful;
  m_idx_row = step;

  PrintDomainE();
  PrintSurfaceF();
  PrintSurfaceQ();
#if defined(MFEM_USE_GSLIB)
  PrintProbeE();
  PrintProbeB();
#endif
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    PrintSurfaceI();
  }

  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::TRANSIENT)
  {
    PrintPortVI();
  }
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    PrintPortS();
  }
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    PrintEig();
    PrintEigPortEPR();
    PrintEigPortQ();
  }
}

// Explicit template instantiation.
template class PostOperatorCSV<ProblemType::DRIVEN>;
template class PostOperatorCSV<ProblemType::EIGENMODE>;
template class PostOperatorCSV<ProblemType::ELECTROSTATIC>;
template class PostOperatorCSV<ProblemType::MAGNETOSTATIC>;
template class PostOperatorCSV<ProblemType::TRANSIENT>;

}  // namespace palace
