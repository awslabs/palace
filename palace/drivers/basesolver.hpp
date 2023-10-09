// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_BASE_SOLVER_HPP
#define PALACE_DRIVERS_BASE_SOLVER_HPP

#include <memory>
#include <string>
#include <vector>
#include <fmt/os.h>

namespace mfem
{

class ParFiniteElementSpaceHierarchy;
class ParMesh;

}  // namespace mfem

namespace palace
{

class ErrorIndicator;
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
    Table(int sp, int p, int p1) : w(sp + p + 7), sp(sp), p(p), w1(p1 + 7), p1(p1) {}
  };
  const Table table;

  // Helper method for creating/appending to output files.
  fmt::ostream OutputFile(const std::string &path, bool append) const
  {
    return append ? fmt::output_file(path, fmt::file::WRONLY | fmt::file::APPEND)
                  : fmt::output_file(path, fmt::file::WRONLY | fmt::file::CREATE |
                                               fmt::file::TRUNC);
  }

  // Common domain postprocessing for all simulation types.
  void PostprocessDomains(const PostOperator &postop, const std::string &name, int step,
                          double time, double E_elec, double E_mag, double E_cap,
                          double E_ind) const;

  // Common surface postprocessing for all simulation types.
  void PostprocessSurfaces(const PostOperator &postop, const std::string &name, int step,
                           double time, double E_elec, double E_mag, double Vinc,
                           double Iinc) const;

  // Common probe postprocessing for all simulation types.
  void PostprocessProbes(const PostOperator &postop, const std::string &name, int step,
                         double time) const;

  // Common field visualization postprocessing for all simulation types.
  void PostprocessFields(const PostOperator &postop, int step, double time,
                         const ErrorIndicator *indicator = nullptr) const;

  // Common error indicator postprocessing for all simulation types.
  void PostprocessErrorIndicator(const PostOperator &postop,
                                 const ErrorIndicator &indicator) const;

public:
  BaseSolver(const IoData &iodata, bool root, int size = 0, int num_thread = 0,
             const char *git_tag = nullptr);
  virtual ~BaseSolver() = default;

  // Performs a solve using the mesh sequence, then reports error indicators.
  virtual ErrorIndicator
  Solve(const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh) const = 0;

  // These methods write different simulation metadata to a JSON file in post_dir.
  void SaveMetadata(const mfem::ParFiniteElementSpaceHierarchy &fespaces) const;
  template <typename SolverType>
  void SaveMetadata(const SolverType &ksp) const;
  void SaveMetadata(const Timer &timer) const;
};

}  // namespace palace

#endif  // PALACE_DRIVERS_BASE_SOLVER_HPP
