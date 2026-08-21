// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CURL_CURL_OPERATOR_HPP
#define PALACE_MODELS_CURL_CURL_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/surfacefluxoperator.hpp"

namespace palace
{

class IoData;
class Mesh;
template <ProblemType T>
class PostOperator;

namespace config
{

struct BoundaryData;
struct MaterialData;
struct PecBoundaryData;
struct SolverData;

}  // namespace config

enum class ProblemType : char;

//
// A class handling discretization of curl-curl problems for magnetostatics.
//
class CurlCurlOperator
{
private:
  // Helper variable for log file printing.
  bool print_hdr;

  // Essential boundary condition attributes.
  mfem::Array<int> dbc_attr;
  std::vector<mfem::Array<int>> dbc_tdof_lists;

  // Cache of screened stiffness operators keyed by their shorted-attribute set, so that
  // excitation steps that short the same set of inactive ports (e.g. every Open-port
  // excitation shorts the identical set of Short ports) reuse one assembled operator
  // instead of reassembling it per step. Each entry owns the per-level essential true DOF
  // lists as well: ParOperator::SetEssentialTrueDofs stores a shallow reference, so the
  // lists must outlive the operator. std::map nodes are address-stable across later
  // insertions, which keeps those references valid for the cache's lifetime.
  struct ScreenedStiffnessCacheEntry
  {
    std::vector<mfem::Array<int>> dbc_tdof_lists;
    std::unique_ptr<Operator> K;
  };
  std::map<std::vector<int>, ScreenedStiffnessCacheEntry> screened_stiffness_cache;

  // Objects defining the finite element spaces for the magnetic vector potential
  // (Nedelec) and magnetic flux density (Raviart-Thomas) on the given mesh. The H1 spaces
  // are used for various purposes throughout the code including postprocessing.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  std::unique_ptr<mfem::RT_FECollection> rt_fec;
  std::unique_ptr<mfem::L2_FECollection> l2_curl_fec;
  FiniteElementSpaceHierarchy nd_fespaces, h1_fespaces;
  FiniteElementSpace rt_fespace;
  std::unique_ptr<FiniteElementSpace> l2_curl_fespace;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Operator for source current excitation.
  SurfaceCurrentOperator surf_j_op;

  // Operator for flux loop excitation.
  SurfaceFluxOperator surf_flux_op;

  // Cached original matrix for flux loop boundary-interior coupling
  mutable std::unique_ptr<ParOperator> K_orig_;

  mfem::Array<int>
  SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                          const std::map<int, config::FluxLoopData> &fluxloop,
                          const mfem::ParMesh &mesh);
  void CheckBoundaryProperties();

public:
  CurlCurlOperator(const config::BoundaryData &boundaries, const config::SolverData &solver,
                   const std::vector<config::MaterialData> &materials,
                   ProblemType problem_type,
                   const std::vector<std::unique_ptr<Mesh>> &mesh);
  CurlCurlOperator(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access to underlying BC operator objects for postprocessing.
  const auto &GetSurfaceCurrentOp() const { return surf_j_op; }
  const auto &GetSurfaceFluxOp() const { return surf_flux_op; }

  // Return the parallel finite element space objects.
  auto &GetNDSpaces() { return nd_fespaces; }
  const auto &GetNDSpaces() const { return nd_fespaces; }
  auto &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  const auto &GetNDSpace() const { return nd_fespaces.GetFinestFESpace(); }
  auto &GetH1Spaces() { return h1_fespaces; }
  const auto &GetH1Spaces() const { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  const auto &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  auto &GetRTSpace() { return rt_fespace; }
  const auto &GetRTSpace() const { return rt_fespace; }
  // In 2D, curl maps H(curl) → L2 (scalar), so use L2 space for B = curl A.
  auto &GetCurlSpace() { return l2_curl_fespace ? *l2_curl_fespace : rt_fespace; }
  const auto &GetCurlSpace() const
  {
    return l2_curl_fespace ? *l2_curl_fespace : rt_fespace;
  }

  // Access the underlying mesh object.
  const auto &GetMesh() const { return GetNDSpace().GetMesh(); }

  // Return the number of true (conforming) dofs on the finest ND space.
  auto GlobalTrueVSize() const { return GetNDSpace().GlobalTrueVSize(); }

  // Construct and return system matrix representing discretized curl-curl operator for
  // Ampere's law.
  std::unique_ptr<Operator> GetStiffnessMatrix();

  // Return the stiffness matrix with extra essential (PEC) attributes beyond those set at
  // construction, without mutating base boundary state. Used in Short mode to treat
  // inactive surface current ports as PEC for a single excitation step. The operator is
  // cached and owned by this object, keyed on the shorted-attribute set, so repeated calls
  // with the same set return a stable reference to one assembled operator rather than
  // reassembling it.
  const Operator &GetScreenedStiffnessMatrix(const mfem::Array<int> &extra_dbc_attr);

  // Zero v on the merged essential set (base Dirichlet plus extra_dbc_attr), clearing the
  // excitation on shorted inactive ports so DIAG_ONE elimination injects no spurious
  // values.
  void ZeroEssentialTrueDofs(const mfem::Array<int> &extra_dbc_attr, Vector &v) const;

  // Construct and return the discrete curl matrix.
  const Operator &GetCurlMatrix() const
  {
    return GetCurlSpace().GetDiscreteInterpolator(GetNDSpace());
  }

  // Assemble the right-hand side source term vector for a current source applied on
  // specified excited boundaries.
  void GetCurrentExcitationVector(int idx, Vector &RHS);

  // Assemble flux loop excitation vector for specified flux loop index
  template <ProblemType T>
  void GetFluxExcitationVector(int idx, Vector &RHS, PostOperator<T> &post_op);

  template <ProblemType T>
  void GetFluxExcitationVector(int idx, Vector &RHS, PostOperator<T> &post_op,
                               Vector *boundary_values);

  // Solve 2D surface curl problem for flux loop boundary conditions
  template <ProblemType T>
  Vector SolveSurfaceCurlProblem(int flux_loop_idx, PostOperator<T> &post_op) const;

  template <ProblemType T>
  void SolveSurfaceCurlProblem(int flux_loop_idx, PostOperator<T> &post_op,
                               Vector &result) const;

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return GetNDSpace().GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_CURL_CURL_OPERATOR_HPP
