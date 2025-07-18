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

// static
Measurement Measurement::Dimensionalize(const Units &units,
                                        const Measurement &nondim_measurement_cache)
{
  Measurement measurement_cache;
  measurement_cache.freq =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(nondim_measurement_cache.freq);
  measurement_cache.ex_idx = nondim_measurement_cache.ex_idx;                        // NONE
  measurement_cache.Jcoeff_excitation = nondim_measurement_cache.Jcoeff_excitation;  // NONE
  measurement_cache.eigenmode_Q = nondim_measurement_cache.eigenmode_Q;              // NONE
  measurement_cache.error_abs = nondim_measurement_cache.error_abs;                  // NONE
  measurement_cache.error_bkwd = nondim_measurement_cache.error_bkwd;                // NONE

  measurement_cache.domain_E_field_energy_all =
      units.Dimensionalize<Units::ValueType::ENERGY>(
          nondim_measurement_cache.domain_E_field_energy_all);
  measurement_cache.domain_H_field_energy_all =
      units.Dimensionalize<Units::ValueType::ENERGY>(
          nondim_measurement_cache.domain_H_field_energy_all);
  for (const auto &e : nondim_measurement_cache.domain_E_field_energy_i)
  {
    measurement_cache.domain_E_field_energy_i.emplace_back(Measurement::DomainData{
        e.idx, units.Dimensionalize<Units::ValueType::ENERGY>(e.energy),
        e.participation_ratio});
  }
  for (const auto &e : nondim_measurement_cache.domain_H_field_energy_i)
  {
    measurement_cache.domain_H_field_energy_i.emplace_back(Measurement::DomainData{
        e.idx, units.Dimensionalize<Units::ValueType::ENERGY>(e.energy),
        e.participation_ratio});
  }
  measurement_cache.lumped_port_capacitor_energy =
      units.Dimensionalize<Units::ValueType::ENERGY>(
          nondim_measurement_cache.lumped_port_capacitor_energy);
  measurement_cache.lumped_port_inductor_energy =
      units.Dimensionalize<Units::ValueType::ENERGY>(
          nondim_measurement_cache.lumped_port_inductor_energy);

  auto dimensionalize_port_post_data =
      [&units](const std::map<int, Measurement::PortPostData> &nondim)
  {
    std::map<int, Measurement::PortPostData> dim;
    for (const auto &[k, data] : nondim)
    {
      dim[k] = Measurement::PortPostData();
      dim[k].P = units.Dimensionalize<Units::ValueType::POWER>(data.P);
      dim[k].V = units.Dimensionalize<Units::ValueType::VOLTAGE>(data.V),
      dim[k].I = units.Dimensionalize<Units::ValueType::CURRENT>(data.I),
      dim[k].I_RLC = {units.Dimensionalize<Units::ValueType::CURRENT>(data.I_RLC[0]),
                      units.Dimensionalize<Units::ValueType::CURRENT>(data.I_RLC[1]),
                      units.Dimensionalize<Units::ValueType::CURRENT>(data.I_RLC[2])};
      dim[k].S = data.S;  // NONE

      dim[k].inductor_energy =
          units.Dimensionalize<Units::ValueType::ENERGY>(data.inductor_energy);
      dim[k].capacitor_energy =
          units.Dimensionalize<Units::ValueType::ENERGY>(data.capacitor_energy);

      dim[k].mode_port_kappa =
          units.Dimensionalize<Units::ValueType::FREQUENCY>(data.mode_port_kappa);
      dim[k].quality_factor = data.quality_factor;                                  // NONE
      dim[k].inductive_energy_participation = data.inductive_energy_participation;  // NONE
    }
    return dim;
  };
  measurement_cache.lumped_port_vi =
      dimensionalize_port_post_data(nondim_measurement_cache.lumped_port_vi);
  measurement_cache.wave_port_vi =
      dimensionalize_port_post_data(nondim_measurement_cache.wave_port_vi);

  measurement_cache.probe_E_field = units.Dimensionalize<Units::ValueType::FIELD_E>(
      nondim_measurement_cache.probe_E_field);
  measurement_cache.probe_B_field = units.Dimensionalize<Units::ValueType::FIELD_B>(
      nondim_measurement_cache.probe_B_field);

  for (const auto &data : nondim_measurement_cache.surface_flux_i)
  {
    auto &flux = measurement_cache.surface_flux_i.emplace_back(data);
    if (data.type == SurfaceFlux::ELECTRIC)
    {
      flux.Phi *= units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
      flux.Phi *= units.GetScaleFactor<Units::ValueType::VOLTAGE>();
    }
    else if (data.type == SurfaceFlux::MAGNETIC)
    {
      flux.Phi *= units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
      flux.Phi *= units.GetScaleFactor<Units::ValueType::CURRENT>();
    }
    else if (data.type == SurfaceFlux::POWER)
    {
      flux.Phi *= units.GetScaleFactor<Units::ValueType::POWER>();
    }
  }

  for (const auto &data : nondim_measurement_cache.interface_eps_i)
  {
    auto &eps = measurement_cache.interface_eps_i.emplace_back(data);
    eps.energy = units.Dimensionalize<Units::ValueType::ENERGY>(data.energy);
  }
  return measurement_cache;
}

// static
Measurement Measurement::Nondimensionalize(const Units &units,
                                           const Measurement &dim_measurement_cache)
{
  Measurement measurement_cache;
  measurement_cache.freq =
      units.Nondimensionalize<Units::ValueType::FREQUENCY>(dim_measurement_cache.freq);
  measurement_cache.ex_idx = dim_measurement_cache.ex_idx;                        // NONE
  measurement_cache.Jcoeff_excitation = dim_measurement_cache.Jcoeff_excitation;  // NONE
  measurement_cache.eigenmode_Q = dim_measurement_cache.eigenmode_Q;              // NONE
  measurement_cache.error_abs = dim_measurement_cache.error_abs;                  // NONE
  measurement_cache.error_bkwd = dim_measurement_cache.error_bkwd;                // NONE

  measurement_cache.domain_E_field_energy_all =
      units.Nondimensionalize<Units::ValueType::ENERGY>(
          dim_measurement_cache.domain_E_field_energy_all);
  measurement_cache.domain_H_field_energy_all =
      units.Nondimensionalize<Units::ValueType::ENERGY>(
          dim_measurement_cache.domain_H_field_energy_all);
  for (const auto &e : dim_measurement_cache.domain_E_field_energy_i)
  {
    measurement_cache.domain_E_field_energy_i.emplace_back(Measurement::DomainData{
        e.idx, units.Nondimensionalize<Units::ValueType::ENERGY>(e.energy),
        e.participation_ratio});
  }
  for (const auto &e : dim_measurement_cache.domain_H_field_energy_i)
  {
    measurement_cache.domain_H_field_energy_i.emplace_back(Measurement::DomainData{
        e.idx, units.Nondimensionalize<Units::ValueType::ENERGY>(e.energy),
        e.participation_ratio});
  }
  measurement_cache.lumped_port_capacitor_energy =
      units.Nondimensionalize<Units::ValueType::ENERGY>(
          dim_measurement_cache.lumped_port_capacitor_energy);
  measurement_cache.lumped_port_inductor_energy =
      units.Nondimensionalize<Units::ValueType::ENERGY>(
          dim_measurement_cache.lumped_port_inductor_energy);

  auto dimensionalize_port_post_data =
      [&units](const std::map<int, Measurement::PortPostData> &nondim)
  {
    std::map<int, Measurement::PortPostData> dim;
    for (const auto &[k, data] : nondim)
    {
      dim[k] = Measurement::PortPostData();
      dim[k].P = units.Nondimensionalize<Units::ValueType::POWER>(data.P);
      dim[k].V = units.Nondimensionalize<Units::ValueType::VOLTAGE>(data.V),
      dim[k].I = units.Nondimensionalize<Units::ValueType::CURRENT>(data.I),
      dim[k].I_RLC = {units.Nondimensionalize<Units::ValueType::CURRENT>(data.I_RLC[0]),
                      units.Nondimensionalize<Units::ValueType::CURRENT>(data.I_RLC[1]),
                      units.Nondimensionalize<Units::ValueType::CURRENT>(data.I_RLC[2])};
      dim[k].S = data.S;  // NONE

      dim[k].inductor_energy =
          units.Nondimensionalize<Units::ValueType::ENERGY>(data.inductor_energy);
      dim[k].capacitor_energy =
          units.Nondimensionalize<Units::ValueType::ENERGY>(data.capacitor_energy);

      dim[k].mode_port_kappa =
          units.Nondimensionalize<Units::ValueType::FREQUENCY>(data.mode_port_kappa);
      dim[k].quality_factor = data.quality_factor;                                  // NONE
      dim[k].inductive_energy_participation = data.inductive_energy_participation;  // NONE
    }
    return dim;
  };
  measurement_cache.lumped_port_vi =
      dimensionalize_port_post_data(dim_measurement_cache.lumped_port_vi);
  measurement_cache.wave_port_vi =
      dimensionalize_port_post_data(dim_measurement_cache.wave_port_vi);

  measurement_cache.probe_E_field = units.Nondimensionalize<Units::ValueType::FIELD_E>(
      dim_measurement_cache.probe_E_field);
  measurement_cache.probe_B_field = units.Nondimensionalize<Units::ValueType::FIELD_B>(
      dim_measurement_cache.probe_B_field);

  for (const auto &data : dim_measurement_cache.surface_flux_i)
  {
    auto &flux = measurement_cache.surface_flux_i.emplace_back(data);
    if (data.type == SurfaceFlux::ELECTRIC)
    {
      flux.Phi /= units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
      flux.Phi /= units.GetScaleFactor<Units::ValueType::VOLTAGE>();
    }
    else if (data.type == SurfaceFlux::MAGNETIC)
    {
      flux.Phi /= units.GetScaleFactor<Units::ValueType::INDUCTANCE>();
      flux.Phi /= units.GetScaleFactor<Units::ValueType::CURRENT>();
    }
    else if (data.type == SurfaceFlux::POWER)
    {
      flux.Phi /= units.GetScaleFactor<Units::ValueType::POWER>();
    }
  }

  for (const auto &data : dim_measurement_cache.interface_eps_i)
  {
    auto &eps = measurement_cache.interface_eps_i.emplace_back(data);
    eps.energy = units.Nondimensionalize<Units::ValueType::ENERGY>(data.energy);
  }
  return measurement_cache;
}

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
      return "unknown";
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

std::vector<std::size_t> _impl::table_expected_filling(std::size_t m_idx_row,
                                                       std::size_t ex_idx_i,
                                                       std::size_t nr_rows,
                                                       std::size_t nr_col_blocks)
{
  // Expected column group filling pattern. Include leading index (freq, ...)
  std::vector<std::size_t> filling_pattern(nr_col_blocks + 1, 0);
  filling_pattern.at(0) = (ex_idx_i == 0) ? m_idx_row : nr_rows;  // index column
  for (std::size_t i = 1; i < ex_idx_i + 1; i++)
  {
    filling_pattern.at(i) = nr_rows;
  }
  filling_pattern.at(ex_idx_i + 1) = m_idx_row;
  return filling_pattern;
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::MoveTableValidateReload(TableWithCSVFile &t_csv_base,
                                                        Table &&t_ref)
{
  // For non-driven solvers or driven with default restart, no table was loaded.
  if (!reload_table)
  {
    t_csv_base.table = std::move(t_ref);
    return;
  }

  // At this point we have a non-default restart. We need to verify that (a) the structure
  // of the table is valid, (b) the cursor location matches the expected restart location.

  auto file = t_csv_base.get_csv_filepath();
  Table &t_base = t_csv_base.table;
  MFEM_VERIFY(!t_base.empty(),
              fmt::format("The data table loaded from path {} was empty, but the "
                          "simulation expected a restart with existing data!",
                          file));

  auto err_msg = fmt::format("The results table loaded from path {} contains pre-existing "
                             "data, but it doest not match the "
                             "expected table structure.",
                             file);
  if (t_base.n_cols() != t_ref.n_cols())
  {
    MFEM_ABORT(fmt::format("{} [Mismatched number of columns: expected {}, got {}.]",
                           err_msg, t_base.n_cols(), t_ref.n_cols()))
  }
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
    // are contiguous. If no error, save row number to compare to expected pattern.
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
  auto expected_ex_idx_nrows = _impl::table_expected_filling(
      row_i, ex_idx_i, nr_expected_measurement_rows, ex_idx_v_all.size());

  // Copy over other options from reference table, since we dont't recover them on load.
  t_csv_base.table.col_options = t_ref.col_options;

  MFEM_VERIFY(base_ex_idx_nrows == expected_ex_idx_nrows,
              fmt::format("{} [Specified restart position is incompatible with reloaded "
                          "file. Row filling by excitation expected {}, got {}]",
                          err_msg, expected_ex_idx_nrows, base_ex_idx_nrows))

  // Don't check index column (frequency) values or size. Size should match with sizing from
  // cursor from printer below. Values will be checked as new frequencies are written.
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeDomainE(const DomainPostOperator &dom_post_op)
{
  using fmt::format;
  domain_E = TableWithCSVFile(post_dir / "domain-E.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * 4 * (1 + dom_post_op.M_i.size());
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);

    t.insert(format("Ee_{}", ex_idx), format("E_elec{} (J)", ex_label), ex_idx);
    t.insert(format("Em_{}", ex_idx), format("E_mag{} (J)", ex_label), ex_idx);
    t.insert(format("Ec_{}", ex_idx), format("E_cap{} (J)", ex_label), ex_idx);
    t.insert(format("Ei_{}", ex_idx), format("E_ind{} (J)", ex_label), ex_idx);

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      t.insert(format("Ee_{}_{}", idx, ex_idx), format("E_elec[{}]{} (J)", idx, ex_label),
               ex_idx);
      t.insert(format("pe_{}_{}", idx, ex_idx), format("p_elec[{}]{}", idx, ex_label),
               ex_idx);
      t.insert(format("Em_{}_{}", idx, ex_idx), format("E_mag[{}]{} (J)", idx, ex_label),
               ex_idx);
      t.insert(format("pm_{}_{}", idx, ex_idx), format("p_mag[{}]{}", idx, ex_label),
               ex_idx);
    }
  }
  MoveTableValidateReload(*domain_E, std::move(t));
  // No longer WriteFullTableTrunc here. We want to potentially reload and check all
  // measurement tables, which may have existing values, before overwriting anything. Just
  // write on first measurement.
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintDomainE()
{
  if (!domain_E)  // trivial check: always written and we are always on root
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(domain_E->table["idx"], row_idx_v, row_i);
  domain_E->table[format("Ee_{}", m_ex_idx)] << measurement_cache.domain_E_field_energy_all;
  domain_E->table[format("Em_{}", m_ex_idx)] << measurement_cache.domain_H_field_energy_all;
  domain_E->table[format("Ec_{}", m_ex_idx)]
      << measurement_cache.lumped_port_capacitor_energy;
  domain_E->table[format("Ei_{}", m_ex_idx)]
      << measurement_cache.lumped_port_inductor_energy;
  for (const auto &data : measurement_cache.domain_E_field_energy_i)
  {
    domain_E->table[format("Ee_{}_{}", data.idx, m_ex_idx)] << data.energy;
    domain_E->table[format("pe_{}_{}", data.idx, m_ex_idx)] << data.participation_ratio;
  }
  for (const auto &data : measurement_cache.domain_H_field_energy_i)
  {
    domain_E->table[format("Em_{}_{}", data.idx, m_ex_idx)] << data.energy;
    domain_E->table[format("pm_{}_{}", data.idx, m_ex_idx)] << data.participation_ratio;
  }
  domain_E->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceF(const SurfacePostOperator &surf_post_op)
{
  if (!(surf_post_op.flux_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_F = TableWithCSVFile(post_dir / "surface-F.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() *
                                              (HasComplexGridFunction<solver_t>() ? 2 : 1) *
                                              surf_post_op.flux_surfs.size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : surf_post_op.flux_surfs)
    {
      switch (data.type)
      {
        case SurfaceFlux::ELECTRIC:
          if (HasComplexGridFunction<solver_t>())
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Re{{Φ_elec[{}]{}}} (C)", idx, ex_label), ex_idx);
            t.insert(format("F_{}_{}_im", idx, ex_idx),
                     format("Im{{Φ_elec[{}]{}}} (C)", idx, ex_label), ex_idx);
          }
          else
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Φ_elec[{}]{} (C)", idx, ex_label), ex_idx);
          }
          break;
        case SurfaceFlux::MAGNETIC:
          if (HasComplexGridFunction<solver_t>())
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Re{{Φ_mag[{}]{}}} (Wb)", idx, ex_label), ex_idx);
            t.insert(format("F_{}_{}_im", idx, ex_idx),
                     format("Im{{Φ_mag[{}]{}}} (Wb)", idx, ex_label), ex_idx);
          }
          else
          {
            t.insert(format("F_{}_{}_re", idx, ex_idx),
                     format("Φ_mag[{}]{} (Wb)", idx, ex_label), ex_idx);
          }
          break;
        case SurfaceFlux::POWER:
          t.insert(format("F_{}_{}_re", idx, ex_idx),
                   format("Φ_pow[{}]{} (W)", idx, ex_label), ex_idx);
          break;
      }
    }
  }
  MoveTableValidateReload(*surface_F, std::move(t));
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceF()
{
  if (!surface_F)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_F->table["idx"], row_idx_v, row_i);
  for (const auto &data : measurement_cache.surface_flux_i)
  {
    surface_F->table[format("F_{}_{}_re", data.idx, m_ex_idx)] << data.Phi.real();
    if (HasComplexGridFunction<solver_t>() &&
        (data.type == SurfaceFlux::ELECTRIC || data.type == SurfaceFlux::MAGNETIC))
    {
      surface_F->table[format("F_{}_{}_im", data.idx, m_ex_idx)] << data.Phi.imag();
    }
  }
  surface_F->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeSurfaceQ(const SurfacePostOperator &surf_post_op)
{
  if (!(surf_post_op.eps_surfs.size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_Q = TableWithCSVFile(post_dir / "surface-Q.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * (2 * surf_post_op.eps_surfs.size());
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      t.insert(format("p_{}_{}", idx, ex_idx), format("p_surf[{}]{}", idx, ex_label),
               ex_idx);
      t.insert(format("Q_{}_{}", idx, ex_idx), format("Q_surf[{}]{}", idx, ex_label),
               ex_idx);
    }
  }
  MoveTableValidateReload(*surface_Q, std::move(t));
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintSurfaceQ()
{
  if (!surface_Q)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_Q->table["idx"], row_idx_v, row_i);

  for (const auto &data : measurement_cache.interface_eps_i)
  {
    surface_Q->table[format("p_{}_{}", data.idx, m_ex_idx)] << data.energy_participation;
    surface_Q->table[format("Q_{}_{}", data.idx, m_ex_idx)] << data.quality_factor;
  }
  surface_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeProbeE(const InterpolationOperator &interp_op)
{
  if (!(interp_op.GetProbes().size() > 0) || !HasEGridFunction<solver_t>())
  {
    return;
  }
  using fmt::format;
  probe_E = TableWithCSVFile(post_dir / "probe-E.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.
  auto v_dim = interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * scale_col * interp_op.GetProbes().size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if constexpr (HasComplexGridFunction<solver_t>())
        {
          t.insert(format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("Re{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label),
                   ex_idx);
          t.insert(format("E{}_{}_{}_im", i_dim, idx, ex_idx),
                   format("Im{{E_{}[{}]{}}} (V/m)", DimLabel(i_dim), idx, ex_label),
                   ex_idx);
        }
        else
        {
          t.insert(format("E{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("E_{}[{}]{} (V/m)", DimLabel(i_dim), idx, ex_label), ex_idx);
        }
      }
    }
  }
  MoveTableValidateReload(*probe_E, std::move(t));
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintProbeE(const InterpolationOperator &interp_op)
{
  if (!probe_E)
  {
    return;
  }
  using fmt::format;
  auto v_dim = interp_op.GetVDim();
  auto probe_field = measurement_cache.probe_E_field;
  MFEM_VERIFY(probe_field.size() == v_dim * interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, interp_op.GetProbes().size(),
                     v_dim * interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_E->table["idx"], row_idx_v, row_i);

  size_t i = 0;
  for (const auto &idx : interp_op.GetProbes())
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
void PostOperatorCSV<solver_t>::InitializeProbeB(const InterpolationOperator &interp_op)
{
  if (!(interp_op.GetProbes().size() > 0) || !HasBGridFunction<solver_t>())
  {
    return;
  }
  using fmt::format;
  probe_B = TableWithCSVFile(post_dir / "probe-B.csv", reload_table);
  Table t;  // Define table locally first due to potential reload.
  auto v_dim = interp_op.GetVDim();
  int scale_col = (HasComplexGridFunction<solver_t>() ? 2 : 1) * v_dim;
  auto nr_expected_measurement_cols =
      1 + ex_idx_v_all.size() * scale_col * interp_op.GetProbes().size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &idx : interp_op.GetProbes())
    {
      for (int i_dim = 0; i_dim < v_dim; i_dim++)
      {
        if (HasComplexGridFunction<solver_t>())
        {
          t.insert(format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("Re{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label),
                   ex_idx);
          t.insert(format("B{}_{}_{}_im", i_dim, idx, ex_idx),
                   format("Im{{B_{}[{}]{}}} (Wb/m²)", DimLabel(i_dim), idx, ex_label),
                   ex_idx);
        }
        else
        {
          t.insert(format("B{}_{}_{}_re", i_dim, idx, ex_idx),
                   format("B_{}[{}]{} (Wb/m²)", DimLabel(i_dim), idx, ex_label), ex_idx);
        }
      }
    }
  }
  MoveTableValidateReload(*probe_B, std::move(t));
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintProbeB(const InterpolationOperator &interp_op)
{
  if (!probe_B)
  {
    return;
  }
  using fmt::format;

  auto v_dim = interp_op.GetVDim();
  auto probe_field = measurement_cache.probe_B_field;
  MFEM_VERIFY(probe_field.size() == v_dim * interp_op.GetProbes().size(),
              format("Size mismatch: expect vector field to have size {} * {} = {}; got {}",
                     v_dim, interp_op.GetProbes().size(),
                     v_dim * interp_op.GetProbes().size(), probe_field.size()))

  CheckAppendIndex(probe_B->table["idx"], row_idx_v, row_i);

  size_t i = 0;
  for (const auto &idx : interp_op.GetProbes())
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
auto PostOperatorCSV<solver_t>::InitializeSurfaceI(const SurfaceCurrentOperator &surf_j_op)
    -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>
{
  if (!(surf_j_op.Size() > 0))
  {
    return;
  }
  using fmt::format;
  surface_I = TableWithCSVFile(post_dir / "surface-I.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.
  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * surf_j_op.Size();
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);
    for (const auto &[idx, data] : surf_j_op)
    {
      t.insert(format("I_{}_{}", idx, ex_idx), format("I_inc[{}]{} (A)", idx, ex_label),
               ex_idx);
    }
  }
  MoveTableValidateReload(*surface_I, std::move(t));
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintSurfaceI(const SurfaceCurrentOperator &surf_j_op,
                                              const Units &units)
    -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>
{
  if (!surface_I)
  {
    return;
  }
  using fmt::format;
  CheckAppendIndex(surface_I->table["idx"], row_idx_v, row_i);
  for (const auto &[idx, data] : surf_j_op)
  {
    auto I_inc_raw = data.GetExcitationCurrent() * measurement_cache.Jcoeff_excitation;
    auto I_inc = units.Dimensionalize<Units::ValueType::CURRENT>(I_inc_raw);
    surface_I->table[format("I_{}_{}", idx, m_ex_idx)] << I_inc;
  }
  surface_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializePortVI(const SpaceOperator &fem_op)
    -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                            U == ProblemType::TRANSIENT,
                        void>
{
  if (!(fem_op.GetLumpedPortOp().Size() > 0))
  {
    return;
  }
  using fmt::format;
  // Currently only works for lumped ports.
  const auto &lumped_port_op = fem_op.GetLumpedPortOp();
  port_V = TableWithCSVFile(post_dir / "port-V.csv", reload_table);
  port_I = TableWithCSVFile(post_dir / "port-I.csv", reload_table);

  Table tV;  // Define table locally first due to potential reload.
  Table tI;

  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * lumped_port_op.Size();
  tV.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  tI.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);

  tV.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  tI.insert("idx", LabelIndexCol(solver_t), -1, 0, PrecIndexCol(solver_t), "");
  for (const auto ex_idx : ex_idx_v_all)
  {
    std::string ex_label = HasSingleExIdx() ? "" : format("[{}]", ex_idx);

    // Print incident signal, if solver supports excitation on ports.
    if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
    {
      auto ex_spec = fem_op.GetPortExcitations().excitations.at(ex_idx);
      for (const auto &idx : ex_spec.lumped_port)
      {
        tV.insert(format("inc{}_{}", idx, ex_idx), format("V_inc[{}]{} (V)", idx, ex_label),
                  ex_idx);
        tI.insert(format("inc{}_{}", idx, ex_idx), format("I_inc[{}]{} (A)", idx, ex_label),
                  ex_idx);
      }
    }
    for (const auto &[idx, data] : lumped_port_op)
    {
      if constexpr (HasComplexGridFunction<solver_t>())
      {
        tV.insert(format("re{}_{}", idx, ex_idx),
                  format("Re{{V[{}]{}}} (V)", idx, ex_label), ex_idx);
        tV.insert(format("im{}_{}", idx, ex_idx),
                  format("Im{{V[{}]{}}} (V)", idx, ex_label), ex_idx);
        tI.insert(format("re{}_{}", idx, ex_idx),
                  format("Re{{I[{}]{}}} (A)", idx, ex_label), ex_idx);
        tI.insert(format("im{}_{}", idx, ex_idx),
                  format("Im{{I[{}]{}}} (A)", idx, ex_label), ex_idx);
      }
      else
      {
        tV.insert(format("re{}_{}", idx, ex_idx), format("V[{}]{} (V)", idx, ex_label),
                  ex_idx);
        tI.insert(format("re{}_{}", idx, ex_idx), format("I[{}]{} (A)", idx, ex_label),
                  ex_idx);
      }
    }
  }
  MoveTableValidateReload(*port_V, std::move(tV));
  MoveTableValidateReload(*port_I, std::move(tI));
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::PrintPortVI(const LumpedPortOperator &lumped_port_op,
                                            const Units &units)
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
  // Postprocess the frequency domain lumped port voltages and currents (complex magnitude
  // = sqrt(2) * RMS).

  CheckAppendIndex(port_V->table["idx"], row_idx_v, row_i);
  CheckAppendIndex(port_I->table["idx"], row_idx_v, row_i);

  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    for (const auto &[idx, data] : lumped_port_op)
    {
      if (data.excitation == m_ex_idx)
      {
        auto Jcoeff = measurement_cache.Jcoeff_excitation;
        double V_inc = data.GetExcitationVoltage() * Jcoeff;
        double I_inc = (std::abs(V_inc) > 0.0)
                           ? data.GetExcitationPower() * Jcoeff * Jcoeff / V_inc
                           : 0.0;

        port_V->table[format("inc{}_{}", idx, m_ex_idx)]
            << units.Dimensionalize<Units::ValueType::VOLTAGE>(V_inc);
        port_I->table[format("inc{}_{}", idx, m_ex_idx)]
            << units.Dimensionalize<Units::ValueType::CURRENT>(I_inc);
      }
    }
  }

  for (const auto &[idx, data] : measurement_cache.lumped_port_vi)
  {
    port_V->table[fmt::format("re{}_{}", idx, m_ex_idx)] << data.V.real();
    port_I->table[fmt::format("re{}_{}", idx, m_ex_idx)] << data.I.real();

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      port_V->table[fmt::format("im{}_{}", idx, m_ex_idx)] << data.V.imag();
      port_I->table[fmt::format("im{}_{}", idx, m_ex_idx)] << data.I.imag();
    }
  }
  port_V->WriteFullTableTrunc();
  port_I->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializePortS(const SpaceOperator &fem_op)
    -> std::enable_if_t<U == ProblemType::DRIVEN, void>
{
  if (!fem_op.GetPortExcitations().IsMultipleSimple() ||
      !((fem_op.GetLumpedPortOp().Size() > 0) xor (fem_op.GetWavePortOp().Size() > 0)))
  {
    return;
  }
  using fmt::format;
  port_S = TableWithCSVFile(post_dir / "port-S.csv", reload_table);

  Table t;  // Define table locally first due to potential reload.

  auto nr_ports = fem_op.GetLumpedPortOp().Size() + fem_op.GetWavePortOp().Size();

  auto nr_expected_measurement_cols = 1 + ex_idx_v_all.size() * nr_ports;
  t.reserve(nr_expected_measurement_rows, nr_expected_measurement_cols);
  t.insert("idx", "f (GHz)", -1, 0, PrecIndexCol(solver_t), "");

  for (const auto ex_idx : ex_idx_v_all)
  {
    // TODO(C++20): Combine identical loops with ranges + projection.
    for (const auto &[o_idx, data] : fem_op.GetLumpedPortOp())
    {
      t.insert(format("abs_{}_{}", o_idx, ex_idx),
               format("|S[{}][{}]| (dB)", o_idx, ex_idx), ex_idx);
      t.insert(format("arg_{}_{}", o_idx, ex_idx),
               format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx), ex_idx);
    }
    for (const auto &[o_idx, data] : fem_op.GetWavePortOp())
    {
      t.insert(format("abs_{}_{}", o_idx, ex_idx),
               format("|S[{}][{}]| (dB)", o_idx, ex_idx), ex_idx);
      t.insert(format("arg_{}_{}", o_idx, ex_idx),
               format("arg(S[{}][{}]) (deg.)", o_idx, ex_idx), ex_idx);
    }
  }
  MoveTableValidateReload(*port_S, std::move(t));
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
  CheckAppendIndex(port_S->table["idx"], row_idx_v, row_i);
  for (const auto &[idx, data] : measurement_cache.lumped_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, m_ex_idx)] << Measurement::Magnitude(data.S);
    port_S->table[format("arg_{}_{}", idx, m_ex_idx)] << Measurement::Phase(data.S);
  }
  for (const auto &[idx, data] : measurement_cache.wave_port_vi)
  {
    port_S->table[format("abs_{}_{}", idx, m_ex_idx)] << Measurement::Magnitude(data.S);
    port_S->table[format("arg_{}_{}", idx, m_ex_idx)] << Measurement::Phase(data.S);
  }
  port_S->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEig()
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  using fmt::format;
  eig = TableWithCSVFile(post_dir / "eig.csv");
  eig->table.reserve(nr_expected_measurement_rows, 6);
  eig->table.insert("idx", "m", -1, 0, PrecIndexCol(solver_t), "");
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
  eig->table["idx"] << row_idx_v;
  eig->table["f_re"] << measurement_cache.freq.real();
  eig->table["f_im"] << measurement_cache.freq.imag();
  eig->table["q"] << measurement_cache.eigenmode_Q;
  eig->table["err_back"] << measurement_cache.error_bkwd;
  eig->table["err_abs"] << measurement_cache.error_abs;
  eig->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEigPortEPR(
    const LumpedPortOperator &lumped_port_op)
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  // TODO(C++20): Make this a filtered iterator in LumpedPortOp.
  for (const auto &[idx, data] : lumped_port_op)
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
  port_EPR = TableWithCSVFile(post_dir / "port-EPR.csv");
  port_EPR->table.reserve(nr_expected_measurement_rows, 1 + ports_with_L.size());
  port_EPR->table.insert("idx", "m", -1, 0, PrecIndexCol(solver_t), "");
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
  port_EPR->table["idx"] << row_idx_v;
  for (const auto idx : ports_with_L)
  {
    auto vi = measurement_cache.lumped_port_vi.at(idx);
    port_EPR->table[format("p_{}", idx)] << vi.inductive_energy_participation;
  }
  port_EPR->WriteFullTableTrunc();
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperatorCSV<solver_t>::InitializeEigPortQ(const LumpedPortOperator &lumped_port_op)
    -> std::enable_if_t<U == ProblemType::EIGENMODE, void>
{
  // TODO(C++20): Make this a filtered iterator in LumpedPortOp.
  for (const auto &[idx, data] : lumped_port_op)
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
  port_Q = TableWithCSVFile(post_dir / "port-Q.csv");
  port_Q->table.reserve(nr_expected_measurement_rows, 1 + ports_with_R.size());
  port_Q->table.insert("idx", "m", -1, 0, PrecIndexCol(solver_t), "");
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
  port_Q->table["idx"] << row_idx_v;
  for (const auto idx : ports_with_R)
  {
    auto vi = measurement_cache.lumped_port_vi.at(idx);
    port_Q->table[format("Ql_{}", idx)] << vi.quality_factor;
    port_Q->table[format("Kl_{}", idx)] << vi.mode_port_kappa;
  }
  port_Q->WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintErrorIndicator(
    bool is_root, const ErrorIndicator::SummaryStatistics &indicator_stats)
{
  if (!is_root)
  {
    return;
  }

  TableWithCSVFile error_indicator(post_dir / "error-indicators.csv");
  error_indicator.table.reserve(1, 4);

  error_indicator.table.insert(Column("norm", "Norm") << indicator_stats.norm);
  error_indicator.table.insert(Column("min", "Minimum") << indicator_stats.min);
  error_indicator.table.insert(Column("max", "Maximum") << indicator_stats.max);
  error_indicator.table.insert(Column("mean", "Mean") << indicator_stats.mean);

  error_indicator.WriteFullTableTrunc();
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::InitializeCSVDataCollection(
    const PostOperator<solver_t> &post_op)
{
  if (!Mpi::Root(post_op.fem_op->GetComm()))
  {
    return;
  }
  InitializeDomainE(post_op.dom_post_op);
  InitializeSurfaceF(post_op.surf_post_op);
  InitializeSurfaceQ(post_op.surf_post_op);
#if defined(MFEM_USE_GSLIB)
  InitializeProbeE(post_op.interp_op);
  InitializeProbeB(post_op.interp_op);
#endif
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    InitializeSurfaceI(post_op.fem_op->GetSurfaceCurrentOp());
  }
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::TRANSIENT)
  {
    InitializePortVI(*post_op.fem_op);
  }
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    InitializePortS(*post_op.fem_op);
  }
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    InitializeEig();
    InitializeEigPortEPR(post_op.fem_op->GetLumpedPortOp());
    InitializeEigPortQ(post_op.fem_op->GetLumpedPortOp());
  }
}

template <ProblemType solver_t>
void PostOperatorCSV<solver_t>::PrintAllCSVData(
    const PostOperator<solver_t> &post_op, const Measurement &non_dim_measurement_cache,
    double idx_value_dimensionful, int step)
{
  if (!Mpi::Root(post_op.fem_op->GetComm()))
  {
    return;
  }
  row_idx_v = idx_value_dimensionful;
  row_i = step;

  // PostOperator acts on a nondimensional measurement cache, we write a dimensional cache.
  measurement_cache = Measurement::Dimensionalize(post_op.units, non_dim_measurement_cache);

  PrintDomainE();
  PrintSurfaceF();
  PrintSurfaceQ();
#if defined(MFEM_USE_GSLIB)
  PrintProbeE(post_op.interp_op);
  PrintProbeB(post_op.interp_op);
#endif
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    PrintSurfaceI(post_op.fem_op->GetSurfaceCurrentOp(), post_op.units);
  }

  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::EIGENMODE ||
                solver_t == ProblemType::TRANSIENT)
  {
    PrintPortVI(post_op.fem_op->GetLumpedPortOp(), post_op.units);
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

template <ProblemType solver_t>
PostOperatorCSV<solver_t>::PostOperatorCSV(const IoData &iodata,
                                           const fem_op_t<solver_t> &fem_op)
{
  if (!Mpi::Root(fem_op.GetComm()))
  {
    return;
  }

  post_dir = iodata.problem.output;

  // Initialize multi-excitation column group index. Only driven or transient support
  // excitations; for other solvers this is default to a single idx=0.
  if constexpr (solver_t == ProblemType::DRIVEN || solver_t == ProblemType::TRANSIENT)
  {
    auto excitation_helper = fem_op.GetPortExcitations();
    ex_idx_v_all.clear();
    ex_idx_v_all.reserve(excitation_helper.Size());
    std::transform(excitation_helper.begin(), excitation_helper.end(),
                   std::back_inserter(ex_idx_v_all),
                   [](const auto &pair) { return pair.first; });
    // Default to the first excitation.
    ex_idx_i = 0;
    m_ex_idx = ex_idx_v_all.front();
  }

  // Driven solver: can have non-trivial restart.
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    nr_expected_measurement_rows = iodata.solver.driven.sample_f.size();
    reload_table = (iodata.solver.driven.restart != 1);

    row_i = std::size_t(iodata.solver.driven.restart - 1) % nr_expected_measurement_rows;
    ex_idx_i = std::size_t(iodata.solver.driven.restart - 1) / nr_expected_measurement_rows;
    m_ex_idx = ex_idx_v_all.at(ex_idx_i);
  }

  // Non-driven solver: get nr_expected_measurement_rows to reserve table space.
  if (solver_t == ProblemType::EIGENMODE)
  {
    nr_expected_measurement_rows = iodata.solver.eigenmode.n;
  }
  else if (solver_t == ProblemType::ELECTROSTATIC)
  {
    nr_expected_measurement_rows = iodata.solver.electrostatic.n_post;
  }
  else if (solver_t == ProblemType::MAGNETOSTATIC)
  {
    nr_expected_measurement_rows = iodata.solver.magnetostatic.n_post;
  }
  else if (solver_t == ProblemType::TRANSIENT)
  {
    // Estimate number for fixed (linear) stepping.
    nr_expected_measurement_rows =
        std::size_t(iodata.solver.transient.max_t / iodata.solver.transient.delta_t) + 1;
  }
}

// Explicit template instantiation.
template class PostOperatorCSV<ProblemType::DRIVEN>;
template class PostOperatorCSV<ProblemType::EIGENMODE>;
template class PostOperatorCSV<ProblemType::ELECTROSTATIC>;
template class PostOperatorCSV<ProblemType::MAGNETOSTATIC>;
template class PostOperatorCSV<ProblemType::TRANSIENT>;

// Function explict needed testing since everywhere it's through PostOperator.
// TODO(C++20): with requires, we won't need a second template.

template auto PostOperatorCSV<ProblemType::DRIVEN>::InitializePortVI<ProblemType::DRIVEN>(
    const SpaceOperator &fem_op) -> void;

}  // namespace palace
