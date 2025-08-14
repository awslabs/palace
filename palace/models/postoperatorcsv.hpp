#ifndef PALACE_MODELS_POST_OPERATOR_CSV_HPP
#define PALACE_MODELS_POST_OPERATOR_CSV_HPP

#include <memory>
#include <optional>
#include "fem/errorindicator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/tablecsv.hpp"
#include "utils/units.hpp"

namespace palace
{
class IoData;
class DomainPostOperator;
class SurfacePostOperator;
class InterpolationOperator;
class SurfaceCurrentOperator;
class LumpedPortOperator;

// Advance declaration.
template <ProblemType solver_t>
class PostOperator;

// Statically map solver (ProblemType) to finite element operator.

template <ProblemType solver_t>
struct fem_op_map_type
{
  using type = SpaceOperator;
};
template <>
struct fem_op_map_type<ProblemType::ELECTROSTATIC>
{
  using type = LaplaceOperator;
};
template <>
struct fem_op_map_type<ProblemType::MAGNETOSTATIC>
{
  using type = CurlCurlOperator;
};

template <ProblemType solver_t>
using fem_op_t = typename fem_op_map_type<solver_t>::type;

// Results of measurements on fields. Not all measurements are sensible to define for all
// solvers.
struct Measurement
{
  // Mini storage structs for data measurements.
  struct DomainData
  {
    int idx;
    double energy;
    double participation_ratio;
  };

  struct FluxData
  {
    int idx;                   // Surface index
    std::complex<double> Phi;  // Integrated flux
    SurfaceFlux type;
  };

  struct InterfaceData
  {
    int idx;                      // Interface index
    double energy;                // Surface Electric Field Energy
    double tandelta;              // Dissipation tangent tan(δ)
    double energy_participation;  // ratio of interface energy / total_energy
    double quality_factor;        // 1 / (energy_participation * tan δ)
  };

  // Data for both lumped and wave port.
  struct PortPostData
  {
    std::complex<double> P = 0.0;
    std::complex<double> V = 0.0;
    std::complex<double> I = 0.0;
    // Separate R, L, and C branches for current via Z.
    std::array<std::complex<double>, 3> I_RLC = {0.0, 0.0, 0.0};

    // S-Parameter.
    std::complex<double> S = 0.0;

    // Energies (currently only for lumped port).
    double inductor_energy = 0.0;   // E_ind = ∑_j 1/2 L_j I_mj².
    double capacitor_energy = 0.0;  // E_cap = ∑_j 1/2 C_j V_mj².

    // Resistive lumped port (only eigenmode).
    double mode_port_kappa = 0.0;
    double quality_factor = mfem::infinity();

    // Inductive lumped port (only eigenmode).
    double inductive_energy_participation = 0.0;
  };

  // "Pseudo-measurements": input required during measurement or data which is stored here
  // in order to pass it along to the printers.

  int ex_idx = 0;  // driven

  std::complex<double> freq = {0.0, 0.0};  // driven || eigenvalue.

  // Modulation factor for input excitation:
  // - I_inc(t) = J(t) I_in, for transient
  // - I_inc(omega) = I_in, for driven so that Jcoeff_excitation = 1.0
  double Jcoeff_excitation = 1.0;  // transient || driven

  // Eigenmode data including error from solver.
  double eigenmode_Q = 0.0;
  double error_bkwd = 0.0;
  double error_abs = 0.0;

  // "Actual measurements".

  double domain_E_field_energy_all = 0.0;
  double domain_H_field_energy_all = 0.0;

  std::vector<DomainData> domain_E_field_energy_i;
  std::vector<DomainData> domain_H_field_energy_i;

  double lumped_port_capacitor_energy = 0.0;
  double lumped_port_inductor_energy = 0.0;

  std::map<int, PortPostData> lumped_port_vi;
  std::map<int, PortPostData> wave_port_vi;

  // Probe data is ordered as [Fx1, Fy1, Fz1, Fx2, Fy2, Fz2, ...].
  // TODO: Replace with proper matrix: mdspan (C++23) / Eigen.
  std::vector<std::complex<double>> probe_E_field;
  std::vector<std::complex<double>> probe_B_field;

  std::vector<FluxData> surface_flux_i;
  std::vector<InterfaceData> interface_eps_i;

  // Dimensionalize and nondimensionalize a set of measurements
  static Measurement Dimensionalize(const Units &units,
                                    const Measurement &nondim_measurement_cache);
  static Measurement Nondimensionalize(const Units &units,
                                       const Measurement &dim_measurement_cache);
  // Helpers for converting complex variable to magnitude in dB and phase.
  static double Magnitude(std::complex<double> x) { return 20.0 * std::log10(std::abs(x)); }
  static double Phase(std::complex<double> x) { return std::arg(x) * 180.0 / M_PI; }
};

namespace _impl
{
// Filling pattern of rows of column groups — needed to validate reload position of
// previous data. Make it public in an _impl namespace for testing.
std::vector<std::size_t> table_expected_filling(std::size_t m_idx_row, std::size_t ex_idx_i,
                                                std::size_t nr_rows,
                                                std::size_t nr_col_blocks);

}  // namespace _impl

// Helper class to PostOperator to collect csv tables and printers for measurement that will
// be saved to file. This class contains a pointer to the corresponding PostOperator class
// and is a friend to a PostOperator class; this is equivalent to having these members
// and methods in PostOperator. It exists for code clarity.
template <ProblemType solver_t>
class PostOperatorCSV
{
protected:
  // Copy savepath from PostOperator for simpler dependencies.
  fs::path post_dir;
  bool reload_table = false;  // True only for driven simulation with non-default restart

  // Dimensionalized measurement cache. Converted from the PostOperator member variable.
  Measurement measurement_cache;

  // Cursor location & cursor value.

  std::size_t row_i = 0;     // Plain count of current row  (measurement index)
  std::size_t ex_idx_i = 0;  // Plain count of current column group (excitation)

  double row_idx_v;          // Value of row index (time, freq..); must be dimensionful
  std::size_t m_ex_idx = 0;  // ex_idx_v: Excitation index value (= ex_idx_v_all[ex_idx_i])

  // Required in validation of re-loaded table (driven), otherwise just to reserve space.
  // Transient (adaptive time-stepping) or eigenvalue (converged eigenvalues) solver output
  // may differ from expectation.
  std::size_t nr_expected_measurement_rows = 1;

  // Stored column groups (excitations). Default single "0" for solvers without excitations.
  std::vector<std::size_t> ex_idx_v_all = {std::size_t(0)};
  bool HasSingleExIdx() const { return ex_idx_v_all.size() == 1; }

  void MoveTableValidateReload(TableWithCSVFile &t_csv_base, Table &&t_ref);

  // Data tables.
  //
  // These are all std::optional since: (a) should only be instantiated on the root mpi
  // process, (b) they should only be written if the data is non-empty.

  // Initialize and print methods for various output quantities: The initialize methods
  // prepare the tables for data insertion, whilst the print methods insert data
  // appropriately. Methods are only enabled when valid given the problem type.

  // Base (all solvers)
  std::optional<TableWithCSVFile> domain_E;
  void InitializeDomainE(const DomainPostOperator &dom_post_op);
  void PrintDomainE();

  std::optional<TableWithCSVFile> surface_F;
  void InitializeSurfaceF(const SurfacePostOperator &surf_post_op);
  void PrintSurfaceF();

  std::optional<TableWithCSVFile> surface_Q;
  void InitializeSurfaceQ(const SurfacePostOperator &surf_post_op);
  void PrintSurfaceQ();

  std::optional<TableWithCSVFile> probe_E;
  void InitializeProbeE(const InterpolationOperator &interp_op);
  void PrintProbeE(const InterpolationOperator &interp_op);

  std::optional<TableWithCSVFile> probe_B;
  void InitializeProbeB(const InterpolationOperator &interp_op);
  void PrintProbeB(const InterpolationOperator &interp_op);

  // TODO(C++20): Upgrade SFINAE to C++20 concepts to simplify static selection since we can
  // just use `void Function(...) requires (solver_t == Type::A);`.

  // Driven + Transient
  std::optional<TableWithCSVFile> surface_I;
  template <ProblemType U = solver_t>
  auto InitializeSurfaceI(const SurfaceCurrentOperator &surf_j_op)
      -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>;
  template <ProblemType U = solver_t>
  auto PrintSurfaceI(const SurfaceCurrentOperator &surf_j_op, const Units &units)
      -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>;

  // Eigenmode + Driven + Transient
  std::optional<TableWithCSVFile> port_V;
  std::optional<TableWithCSVFile> port_I;
  template <ProblemType U = solver_t>
  auto InitializePortVI(const SpaceOperator &fem_op)
      -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                              U == ProblemType::TRANSIENT,
                          void>;
  template <ProblemType U = solver_t>
  auto PrintPortVI(const LumpedPortOperator &lumped_port_op, const Units &units)
      -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                              U == ProblemType::TRANSIENT,
                          void>;

  // Driven
  std::optional<TableWithCSVFile> port_S;
  template <ProblemType U = solver_t>
  auto InitializePortS(const SpaceOperator &fem_op)
      -> std::enable_if_t<U == ProblemType::DRIVEN, void>;
  template <ProblemType U = solver_t>
  auto PrintPortS() -> std::enable_if_t<U == ProblemType::DRIVEN, void>;

  // Eigenmode
  std::optional<TableWithCSVFile> eig;
  template <ProblemType U = solver_t>
  auto InitializeEig() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;
  template <ProblemType U = solver_t>
  auto PrintEig() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;

  std::vector<int> ports_with_L;
  std::vector<int> ports_with_R;
  std::optional<TableWithCSVFile> port_EPR;
  template <ProblemType U = solver_t>
  auto InitializeEigPortEPR(const LumpedPortOperator &lumped_port_op)
      -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;
  template <ProblemType U = solver_t>
  auto PrintEigPortEPR() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;

  std::optional<TableWithCSVFile> port_Q;
  template <ProblemType U = solver_t>
  auto InitializeEigPortQ(const LumpedPortOperator &lumped_port_op)
      -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;
  template <ProblemType U = solver_t>
  auto PrintEigPortQ() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;

public:
  // Print all data from nondim_measurement_cache.
  void PrintAllCSVData(const PostOperator<solver_t> &post_op,
                       const Measurement &nondim_measurement_cache,
                       double idx_value_dimensionful, int step);

  // Driven specific overload for specifying excitation index
  template <ProblemType U = solver_t>
  auto PrintAllCSVData(const PostOperator<solver_t> &post_op,
                       const Measurement &nondim_measurement_cache,
                       double idx_value_dimensionful, int step, int ex_idx)
      -> std::enable_if_t<U == ProblemType::DRIVEN, void>
  {
    m_ex_idx = ex_idx;
    PrintAllCSVData(post_op, nondim_measurement_cache, idx_value_dimensionful, step);
  }

  // Special case of global indicator — init and print all at once.
  void PrintErrorIndicator(bool is_root,
                           const ErrorIndicator::SummaryStatistics &indicator_stats);

  // "Delayed ctor" so that PostOperator can call it once it is fully constructed.
  // Set-up all files to be called from post_op.
  void InitializeCSVDataCollection(const PostOperator<solver_t> &post_op);

  explicit PostOperatorCSV(const IoData &iodata, const fem_op_t<solver_t> &fem_op);
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_CSV_HPP
