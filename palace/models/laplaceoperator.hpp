// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LAPLACE_OPERATOR_HPP
#define PALACE_MODELS_LAPLACE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/singularassembly.hpp"
#include "fem/singularfeatures.hpp"
#include "fem/singularfield.hpp"
#include "fem/singularsystem.hpp"
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

struct SingularOperatorDiagnostics
{
  int convention_version;
  int standard_order;
  int singular_order;
  int quadrature_order;
  double quadrature_absolute_tolerance;
  double quadrature_relative_tolerance;
  int quadrature_maximum_subdivisions;
  HYPRE_BigInt quadrature_leaf_count;
  int quadrature_maximum_depth;
  HYPRE_BigInt duffy_reference_table_maximum_entries;
  HYPRE_BigInt duffy_reference_cache_hits;
  HYPRE_BigInt h1_enrichment_dofs;
  HYPRE_BigInt nd_enrichment_dofs;
  double standard_diagonal_minimum;
  double standard_diagonal_maximum;
  double standard_diagonal_spread;
  double enrichment_diagonal_minimum;
  double enrichment_diagonal_maximum;
  double enrichment_diagonal_spread;
  double combined_diagonal_spread;
  HYPRE_BigInt feature_patch_count = 0;
  HYPRE_BigInt feature_patch_sum_standard_dofs = 0;
  HYPRE_BigInt feature_patch_sum_enrichment_dofs = 0;
  HYPRE_BigInt feature_patch_maximum_standard_dofs = 0;
  HYPRE_BigInt feature_patch_maximum_enrichment_dofs = 0;
  int feature_patch_minimum_enrichment_multiplicity = 0;
  int feature_patch_maximum_enrichment_multiplicity = 0;
};

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

  // Optional additive singular space and its assembled true-DOF data. Combined vectors
  // use rank-local ordering [standard true DOFs, owned enrichment true DOFs].
  const fem::singular::FeatureTopology *singular_features;
  const fem::singular::TriangleFeatureTopology *triangle_singular_features;
  const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids;
  int standard_order;
  int singular_order;
  fem::singular::AdaptiveAssemblyOptions singular_assembly_options;
  std::unique_ptr<fem::singular::DofTopology> singular_dofs;
  std::unique_ptr<fem::singular::TriangleDofTopology> triangle_singular_dofs;
  std::unique_ptr<fem::singular::ParallelDofNumbering> singular_numbering;
  std::unique_ptr<mfem::HypreParMatrix> singular_gradient;
  std::unique_ptr<mfem::HypreParMatrix> singular_unconstrained_stiffness;
  std::unique_ptr<mfem::HypreParMatrix> singular_eliminated_stiffness;
  std::unique_ptr<mfem::HypreParMatrix> singular_constrained_standard_stiffness;
  fem::singular::ParallelFeaturePatches singular_feature_patches;
  std::unique_ptr<mfem::HypreParMatrix> singular_stiffness_enrichment_error;
  std::unique_ptr<mfem::HypreParMatrix> singular_stiffness_coupling_error;
  mfem::Array<int> singular_essential_true_dofs;
  std::vector<std::unique_ptr<mfem::HypreParMatrix>> singular_h1_prolongations;
  std::vector<mfem::Array<int>> combined_h1_dbc_tdof_lists;
  std::unique_ptr<SingularOperatorDiagnostics> singular_diagnostics;

  mfem::Array<int>
  SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                          const std::map<int, config::TerminalData> &terminal,
                          const mfem::ParMesh &mesh);
  std::map<int, mfem::Array<int>>
  ConstructSources(const std::map<int, config::TerminalData> &terminal);

public:
  LaplaceOperator(
      const config::BoundaryData &boundaries, const config::SolverData &solver,
      const std::vector<config::MaterialData> &materials, ProblemType problem_type,
      const std::vector<std::unique_ptr<Mesh>> &mesh,
      const fem::singular::FeatureTopology *singular_features = nullptr,
      const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids = nullptr,
      const fem::singular::TriangleFeatureTopology *triangle_singular_features = nullptr);
  LaplaceOperator(
      const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh,
      const fem::singular::FeatureTopology *singular_features = nullptr,
      const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids = nullptr,
      const fem::singular::TriangleFeatureTopology *triangle_singular_features = nullptr);

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
  auto &GetMesh() { return GetH1Space().GetMesh(); }
  const auto &GetMesh() const { return GetH1Space().GetMesh(); }

  // Return the number of true (conforming) dofs on the finest H1 space.
  auto GlobalTrueVSize() const
  {
    return GetH1Space().GlobalTrueVSize() +
           (singular_numbering ? singular_numbering->h1.global_size : 0);
  }

  bool HasSingularEnrichment() const
  {
    return singular_features != nullptr || triangle_singular_features != nullptr;
  }
  bool HasTriangleSingularEnrichment() const
  {
    return triangle_singular_features != nullptr;
  }
  const mfem::HypreParMatrix &GetUnconstrainedStiffnessMatrix() const;
  const mfem::HypreParMatrix &GetSingularStandardStiffnessMatrix() const;
  const fem::singular::ParallelFeaturePatches &GetSingularFeaturePatches() const;
  const SingularOperatorDiagnostics &GetSingularDiagnostics() const;
  std::vector<const Operator *> GetCombinedH1ProlongationOperators() const;
  const std::vector<mfem::Array<int>> &GetCombinedH1DbcTDofLists() const
  {
    return combined_h1_dbc_tdof_lists;
  }
  double GetSingularStiffnessEnergyErrorBound(const mfem::Vector &combined_true_dofs) const;
  std::unique_ptr<fem::singular::EnrichedH1FieldEvaluator> GetSingularFieldEvaluator();
  std::unique_ptr<fem::singular::TriangleEnrichedH1FieldEvaluator>
  GetTriangleSingularFieldEvaluator();

  // Construct and return system matrix representing discretized Laplace operator for
  // Gauss's law.
  std::unique_ptr<Operator> GetStiffnessMatrix();

  // Construct and return the discrete gradient matrix.
  const Operator &GetGradMatrix() const;

  // Assemble the solution boundary conditions and right-hand side vector for a nonzero
  // prescribed voltage on the specified surface index.
  void GetExcitationVector(int idx, const Operator &K, Vector &X, Vector &RHS);

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return GetH1Space().GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_LAPLACE_OPERATOR_HPP
