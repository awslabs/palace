#ifndef PALACE_MODELS_POST_OPERATOR_CSV_HPP
#define PALACE_MODELS_POST_OPERATOR_CSV_HPP

#include <memory>
#include <optional>
#include "fem/errorindicator.hpp"
#include "utils/configfile.hpp"
#include "utils/tablecsv.hpp"
#include "utils/units.hpp"

namespace palace
{

// Advance declaration.
template <ProblemType solver_t>
class PostOperator;

// Results of measurements on fields. All values here should be in SI units (i.e. we call
// Dimensionalize<units>(value) on the results before storage here). Not all measurements
// are sensible to define for all solvers.
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
    double abs_S_ij = 0.0;
    double arg_S_ij = 0.0;

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
};

// Helper class to PostOperator to collect csv tables and printers for measurement that will
// be saved to file. This class contains a pointer to the corresponding PostOperator class
// and is a friend to a PostOperator class; this is equivalent to having these members
// and methods in PostOperator. It exists for code clarity.
template <ProblemType solver_t>
class PostOperatorCSV
{
  PostOperator<solver_t> *post_op = nullptr;
  int nr_expected_measurement_rows = 1;

  // Dimensionalized measurement cache. Converted from the PostOperator member variable.
  Measurement measurement_cache;

  // Current measurement step index being written.
  int m_idx_row;

  // Current measurement index values for printing of primary "index" column; assumed
  // dimensionful like all measurements.
  double m_idx_value;

  // List of all excitations (column blocks). Single "0" default for solvers that don't
  // support excitations.
  std::vector<int> excitation_idx_all = {int(0)};
  bool SingleColBlock() const { return excitation_idx_all.size() == 1; }

  // Current measurement excitation index.
  int m_ex_idx = 0;

  // These are all std::optional since: (a) should only be instantiated on the root mpi
  // process, (b) they should only be written if the data is non-empty.

  // Base (all solvers)
  std::optional<TableWithCSVFile> domain_E;
  void InitializeDomainE();
  void PrintDomainE();

  std::optional<TableWithCSVFile> surface_F;
  void InitializeSurfaceF();
  void PrintSurfaceF();

  std::optional<TableWithCSVFile> surface_Q;
  void InitializeSurfaceQ();
  void PrintSurfaceQ();

  std::optional<TableWithCSVFile> probe_E;
  void InitializeProbeE();
  void PrintProbeE();

  std::optional<TableWithCSVFile> probe_B;
  void InitializeProbeB();
  void PrintProbeB();

  // TODO: Upgrade SFINAE to C++20 concepts to simplify static selection since we can just
  // use `void Function(...) requires (solver_t == Type::A);`.

  // Eigenmode + Driven + Transient
  std::optional<TableWithCSVFile> port_V;
  std::optional<TableWithCSVFile> port_I;

  // Initialize and print methods for various output quantities. The initialize methods
  // prepare the tables for data insertion, whilst the print methods insert data
  // appropriately. Methods are only enabled when valid given the problem type.

  template <ProblemType U = solver_t>
  auto InitializePortVI()
      -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                              U == ProblemType::TRANSIENT,
                          void>;

  template <ProblemType U = solver_t>
  auto PrintPortVI()
      -> std::enable_if_t<U == ProblemType::EIGENMODE || U == ProblemType::DRIVEN ||
                              U == ProblemType::TRANSIENT,
                          void>;

  // Driven + Transient
  std::optional<TableWithCSVFile> surface_I;

  template <ProblemType U = solver_t>
  auto InitializeSurfaceI()
      -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>;

  template <ProblemType U = solver_t>
  auto PrintSurfaceI()
      -> std::enable_if_t<U == ProblemType::DRIVEN || U == ProblemType::TRANSIENT, void>;

  // Driven
  int driven_source_index = -1;
  std::optional<TableWithCSVFile> port_S;

  template <ProblemType U = solver_t>
  auto InitializePortS() -> std::enable_if_t<U == ProblemType::DRIVEN, void>;
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
  auto InitializeEigPortEPR() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;
  template <ProblemType U = solver_t>
  auto PrintEigPortEPR() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;

  std::optional<TableWithCSVFile> port_Q;

  template <ProblemType U = solver_t>
  auto InitializeEigPortQ() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;
  template <ProblemType U = solver_t>
  auto PrintEigPortQ() -> std::enable_if_t<U == ProblemType::EIGENMODE, void>;

public:
  // Set-up all files to be called from post_op.
  void InitializeCSVDataCollection();

  // Print all data from nondim_measurement_cache.
  void PrintAllCSVData(const Measurement &nondim_measurement_cache,
                       double idx_value_dimensionful, int step);

  // Driven specific overload for specifying excitation index
  template <ProblemType U = solver_t>
  auto PrintAllCSVData(const Measurement &nondim_measurement_cache,
                       double idx_value_dimensionful, int step, int ex_idx)
      -> std::enable_if_t<U == ProblemType::DRIVEN, void>
  {
    m_ex_idx = ex_idx;
    PrintAllCSVData(nondim_measurement_cache, idx_value_dimensionful, step);
  }

  // Special case of global indicator — init and print all at once.
  void PrintErrorIndicator(const ErrorIndicator::SummaryStatistics &indicator_stats);

  PostOperatorCSV() = delete;
  PostOperatorCSV(PostOperator<solver_t> *post_op_, int nr_expected_measurement_rows_ = 1)
    : post_op(post_op_), nr_expected_measurement_rows(nr_expected_measurement_rows_) {};
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_CSV_HPP
