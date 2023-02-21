// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_BASE_SOLVER_HPP
#define PALACE_BASE_SOLVER_HPP

#include <memory>
#include <string>
#include <vector>
#include <fmt/os.h>

#include <mfem.hpp>

namespace mfem
{

class ParFiniteElementSpace;
class ParMesh;

}  // namespace mfem

namespace palace
{

class IoData;
class PostOperator;
class Timer;

//
// Base driver class for all simulation types.
//
class BaseSolver
{
protected:
  // Reference to configuration file data (not owned).
  const IoData &iodata;

  // Parameters for writing postprocessing outputs.
  const std::string post_dir;
  const bool root;

  // Table formatting for output files.
  struct Table
  {
    int w;   // Total column width = precision + spaces + 7 extra (signs/exponent)
    int sp;  // Table column spaces
    int p;   // Floating point precision for data
    int w1;  // First column width = precision + 7 extra
    int p1;  // Floating point precision for first column
    Table(int sp_, int p_, int p1_) : w(sp_ + p_ + 7), sp(sp_), p(p_), w1(p1_ + 7), p1(p1_)
    {
    }
  };
  const Table table;

  // Helper method for creating/appending to output files.
  fmt::ostream OutputFile(const std::string &path, bool append) const
  {
    return append ? fmt::output_file(path, fmt::file::WRONLY | fmt::file::APPEND)
                  : fmt::output_file(path, fmt::file::WRONLY | fmt::file::CREATE |
                                               fmt::file::TRUNC);
  }

  // Common postprocessing functions for all simulation types.
  void PostprocessDomains(const PostOperator &postop, const std::string &name, int step,
                          double time, double E_elec, double E_mag, double E_cap,
                          double E_ind) const;
  void PostprocessSurfaces(const PostOperator &postop, const std::string &name, int step,
                           double time, double E_elec, double E_mag, double Vinc,
                           double Iinc) const;
  void PostprocessProbes(const PostOperator &postop, const std::string &name, int step,
                         double time) const;
  void PostprocessFields(const PostOperator &postop, int step, double time) const;

public:
  BaseSolver(const IoData &iodata_, bool root_, int size = 0, int num_thread = 0,
             const char *git_tag = nullptr);
  virtual ~BaseSolver() = default;

  /// Storage for error estimation results from the solve. Required in the
  /// AMR loop. An error indicator is non-negative, whilst an error estimate is
  /// signed.
  struct SolveOutput
  {
    /// Elemental localized error indicators. Used for marking elements for
    /// refinement and coarsening.
    mfem::Array<double> local_error_indicator;
    /// Global error indicator. Used for driving AMR and diagnostics.
    double global_error_indicator;
  };

  virtual SolveOutput Solve(std::vector<std::unique_ptr<mfem::ParMesh>> &mesh,
                            Timer &timer) const = 0;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const mfem::ParFiniteElementSpace &fespace) const;
  void SaveMetadata(int ksp_mult, int ksp_it) const;
  void SaveMetadata(const Timer &timer) const;
};

}  // namespace palace

#endif  // PALACE_BASE_SOLVER_HPP
