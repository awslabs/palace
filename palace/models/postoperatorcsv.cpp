#include "postoperatorcsv.hpp"

#include <mfem.hpp>

#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/iodata.hpp"

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

int PrecIndexCol(const config::ProblemData::Type solver_t)
{
  switch (solver_t)
  {
    case config::ProblemData::Type::DRIVEN:
    case config::ProblemData::Type::TRANSIENT:
      return 8;
    case config::ProblemData::Type::EIGENMODE:
    case config::ProblemData::Type::ELECTROSTATIC:
    case config::ProblemData::Type::MAGNETOSTATIC:
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

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::MoveTableValidateReload(TableWithCSVFile &t_csv_base,
                                                        Table &&t_ref)
{
  if (!may_reload_table())
  {
    t_csv_base.table = std::move(t_ref);
    return;
  }
  // t_base is empty. This happens if there was no file to reload, the file was empty, or
  // the read from disk was invalid. Just check the simulation does not expect a restart.
  auto file = t_csv_base.get_csv_filepath();
  Table &t_base = t_csv_base.table;
  if (t_base.empty())
  {
    if ((ex_idx_i != 0) || (row_i != 0))  // Non-trivial restart
    {
      MFEM_ABORT(fmt::format("The data table loaded from path {} was empty, but the "
                             "simulation expects a non-trivial restart!",
                             file))
    }
    t_base = std::move(t_ref);  // Initializing a new run, resused expected table.
    return;
  }

  // t_base has data in it. We need to verify that (a) the structure of the table is valid,
  // (b) the t_base cursor location matches the expected restart location.
  auto err_msg = fmt::format("The results table loaded from path {} contains pre-existing "
                             "data, but it doest not match the "
                             "expected table structure.",
                             file);
  if (t_base.n_cols() != t_ref.n_cols())
  {
    MFEM_ABORT(fmt::format("{} [Mismatched number of columns: expected {}, got {}.]",
                           err_msg, t_base.n_cols(), t_ref.n_cols()))
  }
  // Large number using signed to unsized conversion.
  std::vector<std::size_t> base_ex_idx_nrows;
  long current_ex_idx_v = std::numeric_limits<long>::max();
  for (std::size_t i = 0; i < t_base.n_cols(); i++)
  {
    auto &t_base_i = t_base[i];
    auto &t_ref_i = t_ref[i];

    if (t_base_i.header_text != t_ref_i.header_text)
    {
      MFEM_ABORT(fmt::format("{} [Mismatched column header: expected {}, got {}.]", err_msg,
                             t_base_i.header_text, t_ref_i.header_text))
    }
    // Since we cannot parse the column name (internal label) from the printed csv file or
    // other options from csv file, we will over-write them from t_ref.
    t_base_i.name = t_ref_i.name;
    t_base_i.column_group_idx = t_ref_i.column_group_idx;
    t_base_i.min_left_padding = t_ref_i.min_left_padding;
    t_base_i.float_precision = t_ref_i.float_precision;
    t_base_i.fmt_sign = t_ref_i.fmt_sign;
    t_base_i.print_as_int = t_ref_i.print_as_int;

    // Check that columns in same group have the same row number. Assumes that column groups
    // are contiguous.
    if (t_base_i.column_group_idx != current_ex_idx_v)
    {
      current_ex_idx_v = t_base_i.column_group_idx;
      base_ex_idx_nrows.push_back(t_base_i.n_rows());
    }
    else
    {
      if (t_base_i.n_rows() != base_ex_idx_nrows.back())
      {
        MFEM_ABORT(fmt::format("{} [Mismatched rows in excitation {}.]", err_msg,
                               current_ex_idx_v))
      }
    }
  }
  // Match expected column group pattern.
  std::vector<std::size_t> expected_ex_idx_nrows;
  expected_ex_idx_nrows.reserve(ex_idx_v_all.size());
  for (std::size_t i = 0; i < ex_idx_v_all.size(); i++)
  {
    if (i < ex_idx_i)
    {
      expected_ex_idx_nrows.emplace_back(nr_expected_measurement_rows);
    }
    else if (i == ex_idx_i)
    {
      expected_ex_idx_nrows.emplace_back(row_i);
    }
    else
    {
      expected_ex_idx_nrows.emplace_back(0);
    }
  }
  MFEM_VERIFY(base_ex_idx_nrows == expected_ex_idx_nrows,
              fmt::format("{} [Specified restart position is incompatible with reloaded "
                          "file. Row filling by excitation expected {}, got {}]",
                          err_msg, expected_ex_idx_nrows, base_ex_idx_nrows))

  // Don't check index column (frequency) values or size. Size should match with sizing from
  // cursor from printer below. Values will be checked as new frequencies are written.
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::InitializeDomainE()
{
  using fmt::format;
  domain_E = TableWithCSVFile(post_op->post_dir / "domain-E.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * 4 * (1 + post_op->dom_post_op.M_i.size());
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);

    t.insert(format("Ee_{}", ex_idx), format("E_elec{} (J)", ex_label));
    t.insert(format("Em_{}", ex_idx), format("E_mag{} (J)", ex_label));
    t.insert(format("Ec_{}", ex_idx), format("E_cap{} (J)", ex_label));
    t.insert(format("Ei_{}", ex_idx), format("E_ind{} (J)", ex_label));

    for (const auto &[idx, data] : post_op->dom_post_op.M_i)
    {
      t.insert(format("Ee_{}_{}", idx, ex_idx), format("E_elec[{}]{} (J)", idx, ex_label));
      t.insert(format("pe_{}_{}", idx, ex_idx), format("p_elec[{}]{}", idx, ex_label));
      t.insert(format("Em_{}_{}", idx, ex_idx), format("E_mag[{}]{} (J)", idx, ex_label));
      t.insert(format("pm_{}_{}", idx, ex_idx), format("p_mag[{}]{}", idx, ex_label));
    }
  }
  MoveTableValidateReload(*domain_E, std::move(t));
  // TODO: domain_E->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintDomainE()
{
  if (!domain_E)  // trivial check: always written and we are always on root
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(domain_E->table["idx"], row_idx_v, row_i);
  domain_E->table[format("Ee_{}", ex_idx_v)] << MCache().domain_E_field_energy_all;
  domain_E->table[format("Em_{}", ex_idx_v)] << MCache().domain_H_field_energy_all;
  domain_E->table[format("Ec_{}", ex_idx_v)] << MCache().lumped_port_capacitor_energy;
  domain_E->table[format("Ei_{}", ex_idx_v)] << MCache().lumped_port_inductor_energy;
  for (const auto &data : MCache().domain_E_field_energy_i)
  {
    domain_E->table[format("Ee_{}_{}", data.idx, ex_idx_v)] << data.energy;
    domain_E->table[format("pe_{}_{}", data.idx, ex_idx_v)] << data.participation_ratio;
  }
  for (const auto &data : MCache().domain_H_field_energy_i)
  {
    domain_E->table[format("Em_{}_{}", data.idx, ex_idx_v)] << data.energy;
    domain_E->table[format("pm_{}_{}", data.idx, ex_idx_v)] << data.participation_ratio;
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
  surface_F = TableWithCSVFile(post_op->post_dir / "surface-F.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() *
                                              (HasComplexGridFunction<solver_t>() ? 2 : 1) *
                                              post_op->surf_post_op.flux_surfs.size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : post_op->surf_post_op.flux_surfs)
    {
      switch (data.type)
      {
        case SurfaceFluxType::ELECTRIC:
          if (HasComplexGridFunction<solver_t>())
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Re{{Φ_elec[{}]{}}} (C)", idx, ex_label));
            t.insert(format("F_{}_{}_im", idx, ex_idx),
                     format("Im{{Φ_elec[{}]{}}} (C)", idx, ex_label));
          }
          else
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Φ_elec[{}]{} (C)", idx, ex_label));
          }
          break;
        case SurfaceFluxType::MAGNETIC:
          if (HasComplexGridFunction<solver_t>())
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Re{{Φ_mag[{}]{}}} (Wb)", idx, ex_label));
            t.insert(format("F_{}_{}_im", idx, ex_idx),
                     format("Im{{Φ_mag[{}]{}}} (Wb)", idx, ex_label));
          }
          else
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Φ_mag[{}]{} (Wb)", idx, ex_label));
          }
          break;
        case SurfaceFluxType::POWER:
          t.insert(format("F_{}_{}_re", idx, ex_idx),
                   format("Φ_pow[{}]{} (W)", idx, ex_label));
          break;
      }
    }
  }
  MoveTableValidateReload(*surface_F, std::move(t));
  // surface_F->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceF()
{
  if (!surface_F)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_F->table["idx"], row_idx_v, row_i);
  for (const auto &data : MCache().surface_flux_i)
  {
    surface_F->table[format("F_{}_{}_re", data.idx, ex_idx_v)] << data.Phi.real();
    if (HasComplexGridFunction<solver_t>() &&
        (data.type == SurfaceFluxType::ELECTRIC || data.type == SurfaceFluxType::MAGNETIC))
    {
      surface_F->table[format("F_{}_{}_im", data.idx, ex_idx_v)] << data.Phi.imag();
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
  surface_Q = TableWithCSVFile(post_op->post_dir / "surface-Q.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * (2 * post_op->surf_post_op.eps_surfs.size());
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : post_op->surf_post_op.eps_surfs)
    {
      t.insert(format("p_{}_{}", idx, ex_idx), format("p_surf[{}]{}", idx, ex_label));
      t.insert(format("Q_{}_{}", idx, ex_idx), format("Q_surf[{}]{}", idx, ex_label));
    }
  }
  MoveTableValidateReload(*surface_Q, std::move(t));
  // surface_Q->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceQ()
{
  if (!surface_Q)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_Q->table["idx"], row_idx_v, row_i);

  for (const auto &data : MCache().interface_eps_i)
  {
    surface_Q->table[format("p_{}_{}", data.idx, ex_idx_v)] << data.energy_participation;
    surface_Q->table[format("Q_{}_{}", data.idx, ex_idx_v)] << data.quality_factor;
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
  probe_E = TableWithCSVFile(post_op->post_dir / "probe-E.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.
  auto v_dim = post_op->interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * scale_col * post_op->interp_op.GetProbes().size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : post_op->interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if constexpr (HasComplexGridFunction<solver_t>())
        {
          t.insert(format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("Re{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label));
          t.insert(format("E{}_{}_{}_im", i_dim, idx, ex_idx),
                   format("Im{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label));
        }
        else
        {
          t.insert(format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("E_{}[{}]{} (V/m)", DimLabel(i_dim), idx, ex_label));
        }
      }
    }
  }
  MoveTableValidateReload(*probe_E, std::move(t));
  // probe_E->WriteFullTableTrunc();
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
  auto probe_field = MCache().probe_E_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_E->table["idx"], row_idx_v, row_i);

  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_E->table[format("E{}_{}_{}_re", i_dim, idx, ex_idx_v)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_E->table[format("E{}_{}_{}_im", i_dim, idx, ex_idx_v)] << val.imag();
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
  probe_B = TableWithCSVFile(post_op->post_dir / "probe-B.csv", may_reload_table());
  Table t;  // Define table locally first due to potential reload.
  auto v_dim = post_op->interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * scale_col * post_op->interp_op.GetProbes().size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : post_op->interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if (HasComplexGridFunction<solver_t>())
        {
          t.insert(format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("Re{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
          t.insert(format("B{}_{}_{}_im", i_dim, idx, ex_idx),
                   format("Im{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
        }
        else
        {
          t.insert(format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("B_{}[{}]{} (Wb/m²)", DimLabel(i_dim), idx, ex_label));
        }
      }
    }
  }
  MoveTableValidateReload(*probe_B, std::move(t));
  // probe_B->WriteFullTableTrunc();
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
  auto probe_field = MCache().probe_B_field;
  MFEM_VERIFY(probe_field.size() == v_dim * post_op->interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, post_op->interp_op.GetProbes().size(),
                     v_dim * post_op->interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_B->table["idx"], row_idx_v, row_i);

  size_t i = 0;
  for (const auto &idx : post_op->interp_op.GetProbes())
  {
    for (int i_dim = 0; i_dim < v_dim; i_dim++)
    {
      auto val = probe_field[i * v_dim + i_dim];
      probe_B->table[format("B{}_{}_{}_re", i_dim, idx, ex_idx_v)] << val.real();
      if (HasComplexGridFunction<solver_t>())
      {
        probe_B->table[format("B{}_{}_{}_im", i_dim, idx, ex_idx_v)] << val.imag();
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
  surface_I = TableWithCSVFile(post_op->post_dir / "surface-I.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.
  const auto &surf_j_op = post_op->fem_op->GetSurfaceCurrentOp();
  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * surf_j_op.Size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : surf_j_op)
    {
      t.insert(format("I_{}_{}", idx, ex_idx), format("I_inc[{}]{} (A)", idx, ex_label));
    }
  }
  MoveTableValidateReload(*surface_I, std::move(t));
  // surface_I->WriteFullTableTrunc();
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
  CheckAppendIndex(surface_I->table["idx"], row_idx_v, row_i);
  for (const auto &[idx, data] : post_op->fem_op->GetSurfaceCurrentOp())
  {
    auto I_inc_raw = data.GetExcitationCurrent() * MCache().Jcoeff_excitation;
    auto I_inc =
        post_op->units.template Dimensionalize<Units::ValueType::CURRENT>(I_inc_raw);
    surface_I->table[format("I_{}_{}", idx, ex_idx_v)] << I_inc;
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
  // Currently only works for lumped ports.
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  port_V = TableWithCSVFile(post_op->post_dir / "port-V.csv", may_reload_table());
  port_I = TableWithCSVFile(post_op->post_dir / "port-I.csv", may_reload_table());

  Table tV;  // Define table locally first due to potential reload.
  Table tI;

  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * lumped_port_op.Size();
  tV.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  tI.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);

  tV.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  tI.insert("idx", LabelIndexCol(solver_t), 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = SingleExIdx() ? "" : format("[{}]", ex_idx);

    // Print incident signal, if solver supports excitation on ports.
    if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                  solver_t == config::ProblemData::Type::TRANSIENT)
    {
      auto ex_spec = post_op->fem_op->GetPortExcitations().excitations.at(ex_idx);
      for (const auto &idx : ex_spec.lumped_port)
      {
        tV.insert(format("inc{}_{}", idx, ex_idx),
                  format("V_inc[{}]{} (V)", idx, ex_label));
        tI.insert(format("inc{}_{}", idx, ex_idx),
                  format("I_inc[{}]{} (A)", idx, ex_label));
      }
    }
    for (const auto &[idx, data] : lumped_port_op)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        tV.insert(format("re{}_{}", idx, ex_idx),
                  format("Re{{V[{}]{}}} (V)", idx, ex_label));
        tV.insert(format("im{}_{}", idx, ex_idx),
                  format("Im{{V[{}]{}}} (V)", idx, ex_label));
        tI.insert(format("re{}_{}", idx, ex_idx),
                  format("Re{{I[{}]{}}} (A)", idx, ex_label));
        tI.insert(format("im{}_{}", idx, ex_idx),
                  format("Im{{I[{}]{}}} (A)", idx, ex_label));
      }
      else
      {
        tV.insert(format("re{}_{}", idx, ex_idx), format("V[{}]{} (V)", idx, ex_label));
        tI.insert(format("re{}_{}", idx, ex_idx), format("I[{}]{} (A)", idx, ex_label));
      }
    }
  }
  MoveTableValidateReload(*port_V, std::move(tV));
  MoveTableValidateReload(*port_I, std::move(tI));
  // port_V->WriteFullTableTrunc();
  // port_I->WriteFullTableTrunc();
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
  // Currently only works for lumped ports.
  const auto &lumped_port_op = post_op->fem_op->GetLumpedPortOp();
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).

  CheckAppendIndex(port_V->table["idx"], row_idx_v, row_i);
  CheckAppendIndex(port_I->table["idx"], row_idx_v, row_i);

  auto unit_V = post_op->units.template GetScaleFactor<Units::ValueType::VOLTAGE>();
  auto unit_A = post_op->units.template GetScaleFactor<Units::ValueType::CURRENT>();

  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    for (const auto &[idx, data] : lumped_port_op)
    {
      if (data.excitation == ex_idx_v)
      {
        auto Jcoeff = MCache().Jcoeff_excitation;
        double V_inc = data.GetExcitationVoltage() * Jcoeff;
        double I_inc = (std::abs(V_inc) > 0.0)
                           ? data.GetExcitationPower() * Jcoeff * Jcoeff / V_inc
                           : 0.0;

        port_V->table[format("inc{}_{}", idx, ex_idx_v)] << V_inc * unit_V;
        port_I->table[format("inc{}_{}", idx, ex_idx_v)] << I_inc * unit_A;
      }
    }
  }

  for (const auto &[idx, data] : MCache().lumped_port_vi)
  {
    port_V->table[fmt::format("re{}_{}", idx, ex_idx_v)] << data.V.real() * unit_V;
    port_I->table[fmt::format("re{}_{}", idx, ex_idx_v)] << data.I.real() * unit_A;

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      port_V->table[fmt::format("im{}_{}", idx, ex_idx_v)] << data.V.imag() * unit_V;
      port_I->table[fmt::format("im{}_{}", idx, ex_idx_v)] << data.I.imag() * unit_A;
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
  if (!post_op->fem_op->GetPortExcitations().IsMultipleSimple() ||
      !((post_op->fem_op->GetLumpedPortOp().Size() > 0) xor
        (post_op->fem_op->GetWavePortOp().Size() > 0)))
  {
    return;
  }
  using fmt::format;
  port_S = TableWithCSVFile(post_op->post_dir / "port-S.csv", may_reload_table());

  Table t;  // Define table locally first due to potential reload.

  auto nr_ports =
      post_op->fem_op->GetLumpedPortOp().Size() + post_op->fem_op->GetWavePortOp().Size();

  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * nr_ports;
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", "f (GHz)", 0, PrecIndexCol(solver_t), "");

  for (const auto ex_idx : ex_idx_v_all)
  {
    // TODO(C++20): Combine identical loops with ranges + projection.
    for (const auto &[o_idx, data] : post_op->fem_op->GetLumpedPortOp())
    {
      t.insert(format("abs_{}_{}", o_idx, ex_idx),
               format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      t.insert(format("arg_{}_{}", o_idx, ex_idx),
               format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
    for (const auto &[o_idx, data] : post_op->fem_op->GetWavePortOp())
    {
      t.insert(format("abs_{}_{}", o_idx, ex_idx),
               format("|S[{}][{}]| (dB)", o_idx, ex_idx));
      t.insert(format("arg_{}_{}", o_idx, ex_idx),
               format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx));
    }
  }
  MoveTableValidateReload(*port_S, std::move(t));
  // port_S->WriteFullTableTrunc();
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
  CheckAppendIndex(port_S->table["idx"], row_idx_v, row_i);
  for (const auto &[idx, data] : MCache().lumped_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, ex_idx_v)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, ex_idx_v)] << data.arg_S_ij;
  }
  for (const auto &[idx, data] : MCache().wave_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, ex_idx_v)] << data.abs_S_ij;
    port_S->table[format("arg_{}_{}", idx, ex_idx_v)] << data.arg_S_ij;
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
  eig->table.insert("idx", "m", 0, PrecIndexCol(solver_t), "");
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
  eig->table["idx"] << row_idx_v;
  eig->table["f_re"] << MCache().freq.real();
  eig->table["f_im"] << MCache().freq.imag();
  eig->table["q"] << MCache().eigenmode_Q;
  eig->table["err_back"] << MCache().error_bkwd;
  eig->table["err_abs"] << MCache().error_abs;
  eig->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeEigPortEPR()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
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
  port_EPR->table.insert("idx", "m", 0, PrecIndexCol(solver_t), "");
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
  port_EPR->table["idx"] << row_idx_v;
  for (const auto idx : ports_with_L)
  {
    auto vi = MCache().lumped_port_vi.at(idx);
    port_EPR->table[format("p_{}", idx)] << vi.inductive_energy_participation;
  }
  port_EPR->WriteFullTableTrunc();
}

template <config::ProblemData::Type solver_t>
template <config::ProblemData::Type U>
auto PostOperatorCSV<solver_t>::InitializeEigPortQ()
    -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, void>
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
  port_Q->table.insert("idx", "m", 0, PrecIndexCol(solver_t), "");
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
  port_Q->table["idx"] << row_idx_v;
  for (const auto idx : ports_with_R)
  {
    auto vi = MCache().lumped_port_vi.at(idx);
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
void PostOperatorCSV<solver_t>::PrintAllCSVData(double idx_value_dimensionful,
                                                int print_row)
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }
  row_idx_v = idx_value_dimensionful;
  row_i = print_row;

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

template <config::ProblemData::Type solver_t>
void PostOperatorCSV<solver_t>::SetUpAndInitialize(const IoData &iodata)
{
  if (!Mpi::Root(post_op->fem_op->GetComm()))
  {
    return;
  }
  // Initialize multi-excitation column group index. Only driven or transient support
  // excitations; for other solvers this is default to a single idx=0.
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN ||
                solver_t == config::ProblemData::Type::TRANSIENT)
  {
    auto excitation_helper = post_op->fem_op->GetPortExcitations();
    ex_idx_v_all.clear();
    ex_idx_v_all.reserve(excitation_helper.Size());
    std::transform(excitation_helper.begin(), excitation_helper.end(),
                   std::back_inserter(ex_idx_v_all),
                   [](const auto &pair) { return pair.first; });
    // Default to the first excitation.
    ex_idx_i = 0;
    ex_idx_v = ex_idx_v_all.front();
  }

  // Driven solver: can have non-trivial restart.
  if constexpr (solver_t == config::ProblemData::Type::DRIVEN)
  {
    nr_expected_measurement_rows = iodata.solver.driven.sample_f.size();

    row_i = std::size_t(iodata.solver.driven.restart - 1) % nr_expected_measurement_rows;
    ex_idx_i = std::size_t(iodata.solver.driven.restart - 1) / nr_expected_measurement_rows;
    ex_idx_v = ex_idx_v_all.at(ex_idx_i);
  }

  // Non-driven solver: get nr_expected_measurement_rows to reserve table space.
  if (solver_t == config::ProblemData::Type::EIGENMODE)
  {
    nr_expected_measurement_rows = iodata.solver.eigenmode.n;
  }
  else if (solver_t == config::ProblemData::Type::ELECTROSTATIC)
  {
    nr_expected_measurement_rows = iodata.solver.electrostatic.n_post;
  }
  else if (solver_t == config::ProblemData::Type::MAGNETOSTATIC)
  {
    nr_expected_measurement_rows = iodata.solver.magnetostatic.n_post;
  }
  else if (solver_t == config::ProblemData::Type::TRANSIENT)
  {
    // Estimate number for fixed (linear) stepping.
    nr_expected_measurement_rows =
        std::size_t(iodata.solver.transient.max_t / iodata.solver.transient.delta_t) + 1;
  }

  InitializeCSVDataCollection();
}

// Explicit template instantiation.
template class PostOperatorCSV<config::ProblemData::Type::DRIVEN>;
template class PostOperatorCSV<config::ProblemData::Type::EIGENMODE>;
template class PostOperatorCSV<config::ProblemData::Type::ELECTROSTATIC>;
template class PostOperatorCSV<config::ProblemData::Type::MAGNETOSTATIC>;
template class PostOperatorCSV<config::ProblemData::Type::TRANSIENT>;

}  // namespace palace
