#ifndef PALACE_MODELS_POST_OPERATOR_CSV_HPP
#define PALACE_MODELS_POST_OPERATOR_CSV_HPP

#include <optional>
#include "fem/errorindicator.hpp"
#include "utils/configfile.hpp"
#include "utils/tablecsv.hpp"

namespace palace
{

// Advance declaration.
template <config::ProblemData::Type solver_t>
class PostOperator;

// Helper class to PostOperator to collect csv tables and printers for measurement that will
// be saved to file. This class contains a pointer to the corresponding PostOperator class
// and is a friend to a PostOperator class; this is equivalent to having these members
// and methods in PostOperator. It exists for code clarity.
template <config::ProblemData::Type solver_t>
class PostOperatorCSV
{
  PostOperator<solver_t> *post_op = nullptr;
  int nr_expected_measurement_rows = 1;

  // Alias for code clarity
  auto MCache() const { return post_op->measurement_cache; }

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

  // Eigenmode + Driven + Transient
  std::optional<TableWithCSVFile> port_V;
  std::optional<TableWithCSVFile> port_I;

  // Initialize and print methods for various output quantities. The initialize methods
  // prepare the tables for data insertion, whilst the print methods insert data
  // appropriately. Methods are only enabled when valid given the problem type.

  void InitializePortVI()
    requires(solver_t == config::ProblemData::Type::EIGENMODE ||
             solver_t == config::ProblemData::Type::DRIVEN ||
             solver_t == config::ProblemData::Type::TRANSIENT);

  void PrintPortVI()
    requires(solver_t == config::ProblemData::Type::EIGENMODE ||
             solver_t == config::ProblemData::Type::DRIVEN ||
             solver_t == config::ProblemData::Type::TRANSIENT);

  // Driven + Transient
  std::optional<TableWithCSVFile> surface_I;

  void InitializeSurfaceI()
    requires(solver_t == config::ProblemData::Type::DRIVEN ||
             solver_t == config::ProblemData::Type::TRANSIENT);

  void PrintSurfaceI()
    requires(solver_t == config::ProblemData::Type::DRIVEN ||
             solver_t == config::ProblemData::Type::TRANSIENT);

  // Driven
  int driven_source_index = -1;
  std::optional<TableWithCSVFile> port_S;

  void InitializePortS()
    requires(solver_t == config::ProblemData::Type::DRIVEN);
  void PrintPortS()
    requires(solver_t == config::ProblemData::Type::DRIVEN);

  // Eigenmode
  std::optional<TableWithCSVFile> eig;

  void InitializeEig()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);
  void PrintEig()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);

  std::vector<int> ports_with_L;
  std::vector<int> ports_with_R;
  std::optional<TableWithCSVFile> port_EPR;

  void InitializeEigPortEPR()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);
  void PrintEigPortEPR()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);

  std::optional<TableWithCSVFile> port_Q;

  void InitializeEigPortQ()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);
  void PrintEigPortQ()
    requires(solver_t == config::ProblemData::Type::EIGENMODE);

public:
  // Set-up all files to be called from post_op.
  void InitializeCSVDataCollection();

  // Print all data from post_op->measurement_cache.
  void PrintAllCSVData(double idx_value_dimensionful, int step);

  // Driven specific overload for specifying excitation index
  template <config::ProblemData::Type U = solver_t>
  auto PrintAllCSVData(double idx_value_dimensionful, int step, int ex_idx)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, void>
  {
    m_ex_idx = ex_idx;
    PrintAllCSVData(idx_value_dimensionful, step);
  }

  // Special case of global indicator — init and print all at once.
  void PrintErrorIndicator(const ErrorIndicator::SummaryStatistics &indicator_stats);

  PostOperatorCSV() = delete;
  PostOperatorCSV(PostOperator<solver_t> *post_op_, int nr_expected_measurement_rows_)
    : post_op(post_op_), nr_expected_measurement_rows(nr_expected_measurement_rows_) {};
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_CSV_HPP
