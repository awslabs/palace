// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LAPLACE_OPERATOR_HPP
#define PALACE_MODELS_LAPLACE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

class IoData;
class Mesh;

namespace config
{

struct BoundaryData;
struct MaterialData;
struct PecBoundaryData;
struct SolverData;
struct TerminalData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling discretization of Laplace problems for electrostatics.
//
class LaplaceOperator
{
private:
  // Helper variable for log file printing.
  bool print_hdr;

  // Essential boundary condition markers.
  mfem::Array<int> dbc_attr;
  std::vector<mfem::Array<int>> dbc_tdof_lists;

  // Objects defining the finite element spaces for the electrostatic potential (H1) and
  // electric field (Nedelec) on the given mesh. The RT spaces are used for error
  // estimation.
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  std::unique_ptr<mfem::ND_FECollection> nd_fec;
  std::vector<std::unique_ptr<mfem::RT_FECollection>> rt_fecs;
  FiniteElementSpaceHierarchy h1_fespaces;
  FiniteElementSpace nd_fespace;
  FiniteElementSpaceHierarchy rt_fespaces;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Boundary attributes for each terminal index.
  std::map<int, mfem::Array<int>> source_attr_lists;

  // Optional volumetric source ρ for -div(ε ∇V) = ρ (not owned). Null except in
  // manufactured-solution verification tests; there is no config path that sets it.
  mfem::Coefficient *rhs_source = nullptr;

  // Optional prescribed Dirichlet values on the essential boundary (not owned). Null means
  // homogeneous (V = 0). Set only by manufactured-solution verification tests.
  mfem::Coefficient *dbc_source = nullptr;

  // Optional prescribed Neumann flux g = ε ∂V/∂n on non-essential boundaries (not owned),
  // applied over neumann_attr. Null means the natural homogeneous-Neumann BC. Set only by
  // manufactured-solution verification tests.
  mfem::Coefficient *neumann_source = nullptr;
  mfem::Array<int> neumann_attr;

  mfem::Array<int>
  SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                          const std::map<int, config::TerminalData> &terminal,
                          const mfem::ParMesh &mesh);
  std::map<int, mfem::Array<int>>
  ConstructSources(const std::map<int, config::TerminalData> &terminal);

public:
  LaplaceOperator(const config::BoundaryData &boundaries, const config::SolverData &solver,
                  const std::vector<config::MaterialData> &materials,
                  ProblemType problem_type, const std::vector<std::unique_ptr<Mesh>> &mesh);
  LaplaceOperator(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access source attribute lists.
  const auto &GetSources() const { return source_attr_lists; }

  // Return the parallel finite element space objects.
  auto &GetH1Spaces() { return h1_fespaces; }
  const auto &GetH1Spaces() const { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  const auto &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  auto &GetNDSpace() { return nd_fespace; }
  const auto &GetNDSpace() const { return nd_fespace; }
  auto &GetRTSpaces() { return rt_fespaces; }
  const auto &GetRTSpaces() const { return rt_fespaces; }
  auto &GetRTSpace() { return rt_fespaces.GetFinestFESpace(); }
  const auto &GetRTSpace() const { return rt_fespaces.GetFinestFESpace(); }

  // Access the underlying mesh object.
  const auto &GetMesh() const { return GetH1Space().GetMesh(); }

  // Return the number of true (conforming) dofs on the finest H1 space.
  auto GlobalTrueVSize() const { return GetH1Space().GlobalTrueVSize(); }

  // Construct and return system matrix representing discretized Laplace operator for
  // Gauss's law.
  std::unique_ptr<Operator> GetStiffnessMatrix();

  // Construct and return the discrete gradient matrix.
  const Operator &GetGradMatrix() const
  {
    return GetNDSpace().GetDiscreteInterpolator(GetH1Space());
  }

  // Assemble the solution boundary conditions and right-hand side vector for a nonzero
  // prescribed voltage on the specified surface index.
  void GetExcitationVector(int idx, const Operator &K, Vector &X, Vector &RHS);

  // Set the volumetric source ρ (see rhs_source; caller retains ownership).
  void SetRhsSource(mfem::Coefficient &source) { rhs_source = &source; }
  bool HasRhsSource() const { return rhs_source != nullptr; }

  // Set prescribed Dirichlet values on the essential boundary (see dbc_source; caller
  // retains ownership).
  void SetDbcCoefficient(mfem::Coefficient &dbc) { dbc_source = &dbc; }

  // Set a prescribed Neumann flux g = ε ∂V/∂n over the given boundary attributes (see
  // neumann_source; caller retains ownership).
  void SetNeumannCoefficient(mfem::Coefficient &g, const std::vector<int> &attributes)
  {
    neumann_source = &g;
    neumann_attr.SetSize(0);
    for (int attr : attributes)
    {
      neumann_attr.Append(attr);
    }
  }

  // Assemble the RHS from the volumetric source, with Dirichlet BCs on all essential
  // boundaries (homogeneous, or the values from SetDbcCoefficient). Requires SetRhsSource.
  void GetSourceExcitationVector(const Operator &K, Vector &X, Vector &RHS);

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return GetH1Space().GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_LAPLACE_OPERATOR_HPP
