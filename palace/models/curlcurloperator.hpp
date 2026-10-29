// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CURL_CURL_OPERATOR_HPP
#define PALACE_MODELS_CURL_CURL_OPERATOR_HPP

#include <deque>
#include <map>
#include <memory>
#include <set>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/superconductorsheetoperator.hpp"
#include "models/surfacecurlsolver.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/surfacefluxoperator.hpp"

namespace palace
{

class IoData;
class Mesh;
class Units;
template <ProblemType T>
class PostOperator;

namespace config
{

struct BoundaryData;
struct MaterialData;
struct PecBoundaryData;
struct SolverData;
struct SuperconductorData;

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

  // Film attributes (from FilmAttributes) that are also Superconductor sheets ("London flux
  // films"). Their whole surface is a free unknown; the hole fluxoid is realized through
  // the sheet penalty, not by pinning DOFs. Declared before dbc_attr because
  // SetUpBoundaryProperties (which initializes dbc_attr) populates this set, so it must be
  // constructed first.
  std::set<int> london_flux_film_attr_;

  // Essential boundary condition attributes.
  mfem::Array<int> dbc_attr;
  std::vector<mfem::Array<int>> dbc_tdof_lists;

  // Rolling two-entry cache of screened stiffness operators, tagged by shorted-attribute
  // set. Steps that short the same ports are solved back-to-back, so equal keys are
  // adjacent and matching the newest entry reuses its operator. Two entries bound memory to
  // O(1) hierarchies while keeping the previously bound operator alive. std::deque never
  // relocates elements, so the shallow tdof references ParOperator::SetEssentialTrueDofs
  // holds stay valid.
  struct ScreenedStiffnessCacheEntry
  {
    std::vector<int> key;
    std::vector<mfem::Array<int>> dbc_tdof_lists;
    std::unique_ptr<Operator> K;
    // Shifted preconditioner for this screened operator, built lazily and kept alive with
    // it.
    std::unique_ptr<Operator> P;
    std::vector<mfem::Array<int>> h1_dbc_tdof_lists;
  };
  std::deque<ScreenedStiffnessCacheEntry> screened_stiffness_cache;

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

  // With two-sided films, the solution spaces on every level are the ND, H1, and RT spaces
  // broken across the films (see fem/brokenspace.hpp); the unbroken ones remain for
  // postprocessing and error estimation. dbc_tdof_lists then holds the broken essential
  // DOFs, and the sheet terms of each level are assembled matrices built on first use.
  std::vector<int> two_sided_attr;
  FiniteElementSpaceHierarchy nd_broken_fespaces, h1_broken_fespaces;
  std::unique_ptr<FiniteElementSpace> rt_broken_fespace;
  mutable std::vector<std::unique_ptr<mfem::HypreParMatrix>> broken_sheet_mats;

  void SetUpBrokenSpaces();

  // The sheet terms on level l of the broken ND hierarchy.
  const mfem::HypreParMatrix &GetBrokenSheetMatrix(std::size_t l) const;

  // Add the sheet terms of level l to an operator of the broken hierarchy, eliminated on
  // its essential true DOFs (nothing to do on unbroken spaces, where the sheet terms are
  // boundary integrators).
  void AddBrokenSheetMatrix(std::size_t l, ParOperator &K_l) const;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Operator for source current excitation.
  SurfaceCurrentOperator surf_j_op;

  // Operator for flux loop excitation.
  SurfaceFluxOperator surf_flux_op;

  // Operator for thin-film superconductor sheet (kinetic inductance) boundaries.
  SuperconductorSheetOperator sc_sheet_op;

  // Flux-loop indices whose film is (partly) a Superconductor sheet — i.e.
  // London flux films. Populated at construction after surf_flux_op is set.
  std::set<int> london_flux_loops_;

  // Preconditioner-only gauge shift (config LondonPCShift); see GetPreconditionerMatrix.
  double london_pc_shift_;

  // Whether the shifted preconditioner carries auxiliary (gradient) operators for the
  // multigrid smoother, and their H1 essential DOFs.
  bool mg_smooth_aux_;
  std::vector<mfem::Array<int>> h1_dbc_tdof_lists_;

  // Cached sheet-only mass operator M_sheet (the boundary term ∫_Σ (1/L_ksq) A_t·v_t). The
  // London flux excitation RHS = M_sheet·a_h is the Euler-Lagrange source of the shifted
  // penalty ½∫_Σ (1/L_ksq)|A_t − a_h|². Lazily assembled in GetFluxExcitationVector.
  mutable std::unique_ptr<Operator> M_sheet_;

  // Cached London domain mass M_L = ∫ (μλ²)⁻¹ A·v, with AᵀM_L·A twice the supercurrent
  // kinetic energy. Lazily assembled in ApplyLondonMass.
  mutable std::unique_ptr<Operator> M_london_;

  // Fluxoid circulation functional c = Curlᵀ·f_hole per London flux loop (f_hole the RT
  // hole-cap flux functional), so cᵀA = ∮_∂hole A·dl by Stokes. Used to normalize the
  // cohomology generator to cᵀa_h = Φ (the imposed fluxoid) and to measure cᵀA. Zeroed on
  // the PEC essential set.
  std::map<int, Vector> london_flux_constraint_;

  void SetUpLondonFluxConstraints();

  // Assemble the London sheet mass M_sheet = ∫_Σ (1/L_ksq) A_t·v_t on first use.
  void EnsureSheetMass();

  // Add the stiffness coefficients: curl-curl (μ⁻¹), London domains ((μλ²)⁻¹), and the
  // kinetic sheet term (1/L_ksq, unbroken spaces only).
  void AddStiffnessCoefficients(MaterialPropertyCoefficient &df,
                                MaterialPropertyCoefficient &f,
                                MaterialPropertyCoefficient &fb) const;

  // True if a superconductor (film or London domain) is present and the shift is positive.
  bool HasShiftedPreconditioner() const;

  // Assemble the London preconditioner (stiffness + shift·(1/µ)∫|A|²) constraining the
  // given per-level essential DOFs (from attributes pc_dbc_attr); preconditioner-only, so
  // the extracted inductance is unchanged. The auxiliary H1 essential DOFs are written to
  // h1_dbc_tdof_lists, which must outlive the result.
  std::unique_ptr<Operator>
  AssembleShiftedPreconditioner(const std::vector<mfem::Array<int>> &pc_dbc_tdof_lists,
                                const mfem::Array<int> &pc_dbc_attr,
                                std::vector<mfem::Array<int>> &h1_dbc_tdof_lists);

  mfem::Array<int>
  SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                          const std::map<int, config::FluxLoopData> &fluxloop,
                          const std::vector<config::SuperconductorData> &superconductor,
                          const mfem::ParMesh &mesh);
  void CheckBoundaryProperties();

public:
  CurlCurlOperator(const config::BoundaryData &boundaries, const config::SolverData &solver,
                   const std::vector<config::MaterialData> &materials,
                   ProblemType problem_type, const Units &units,
                   const std::vector<std::unique_ptr<Mesh>> &mesh);
  CurlCurlOperator(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access to underlying BC operator objects for postprocessing.
  const auto &GetSurfaceCurrentOp() const { return surf_j_op; }
  const auto &GetSurfaceFluxOp() const { return surf_flux_op; }
  auto &GetSuperconductorOp() { return sc_sheet_op; }
  const auto &GetSuperconductorOp() const { return sc_sheet_op; }

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

  // Return the solution spaces: the ND, H1, and curl spaces broken across two-sided
  // superconductor films, if there are any, and otherwise the same as above. The solution
  // vectors (A, the right-hand sides, and B = curl A) are on their true DOFs.
  bool IsBroken() const { return nd_broken_fespaces.GetNumLevels() > 0; }
  auto &GetSolutionNDSpaces() { return IsBroken() ? nd_broken_fespaces : nd_fespaces; }
  const auto &GetSolutionNDSpaces() const
  {
    return IsBroken() ? nd_broken_fespaces : nd_fespaces;
  }
  auto &GetSolutionNDSpace() { return GetSolutionNDSpaces().GetFinestFESpace(); }
  const auto &GetSolutionNDSpace() const
  {
    return GetSolutionNDSpaces().GetFinestFESpace();
  }
  auto &GetSolutionH1Spaces() { return IsBroken() ? h1_broken_fespaces : h1_fespaces; }
  const auto &GetSolutionH1Spaces() const
  {
    return IsBroken() ? h1_broken_fespaces : h1_fespaces;
  }
  auto &GetSolutionCurlSpace() { return IsBroken() ? *rt_broken_fespace : GetCurlSpace(); }
  const auto &GetSolutionCurlSpace() const
  {
    return IsBroken() ? *rt_broken_fespace : GetCurlSpace();
  }

  // Copy the version 0 of a vector on the true DOFs of a broken solution space (the
  // original DOFs, exact away from the two-sided films and the field of one of their sides
  // on them) to a vector on the true DOFs of the unbroken space, for output fields.
  void GetFirstVersion(const FiniteElementSpace &broken_fespace, const Vector &x,
                       Vector &y) const;

  // Access the underlying mesh object.
  const auto &GetMesh() const { return GetNDSpace().GetMesh(); }

  // Return the number of true (conforming) dofs on the finest ND solution space.
  auto GlobalTrueVSize() const { return GetSolutionNDSpace().GlobalTrueVSize(); }

  // Construct and return system matrix representing discretized curl-curl operator for
  // Ampere's law.
  std::unique_ptr<Operator> GetStiffnessMatrix();

  // True if any flux loop drives a London (finite-λ) superconductor film, whose free
  // interior gives the shifted-penalty operator a residual 1-D gradient null space.
  bool HasLondonFluxLoops() const { return !london_flux_loops_.empty(); }

  // True if any domain has a London depth, adding (μλ²)⁻¹ A·v to the stiffness operator.
  bool HasLondonDomains() const { return mat_op.HasLondonDepth(); }

  // Construct the preconditioner matrix for the London magnetostatic path: the stiffness
  // matrix plus a small reluctance-weighted volume mass (london_pc_shift_ · (1/µ) ∫|A|²)
  // that lifts the residual gradient null space so AMS is SPD-solvable in parallel. This
  // shift is applied ONLY to the preconditioner, so the operator, solution, and extracted
  // inductance are unchanged. Returns nullptr when there is no superconductor or the shift
  // is non-positive, signaling that the stiffness matrix itself should be used as
  // preconditioner.
  std::unique_ptr<Operator> GetPreconditionerMatrix();

  // Return the stiffness matrix with extra essential (PEC) attributes beyond those set at
  // construction, without mutating base boundary state. Used in Short mode to treat
  // inactive surface current ports as PEC for a single excitation step. The operator is
  // cached and owned by this object, keyed on the shorted-attribute set, so repeated calls
  // with the same set return a stable reference to one assembled operator rather than
  // reassembling it.
  const Operator &GetScreenedStiffnessMatrix(const mfem::Array<int> &extra_dbc_attr);

  // Shifted preconditioner for a screened step, matching the essential set of the screened
  // operator and owned by the same cache entry (so both survive back-to-back same-key
  // steps). Must follow GetScreenedStiffnessMatrix with the same set. Null when there is no
  // shift.
  const Operator *GetScreenedPreconditionerMatrix(const mfem::Array<int> &extra_dbc_attr);

  // Zero v on the merged essential set (base Dirichlet plus extra_dbc_attr), clearing the
  // excitation on shorted inactive ports so DIAG_ONE elimination injects no spurious
  // values.
  void ZeroEssentialTrueDofs(const mfem::Array<int> &extra_dbc_attr, Vector &v) const;

  // Construct and return the discrete curl matrix (between the solution spaces).
  const Operator &GetCurlMatrix() const
  {
    return GetSolutionCurlSpace().GetDiscreteInterpolator(GetSolutionNDSpace());
  }

  // Apply the sheet mass M_sheet (∫_Σ (1/L_ksq) A_t·v_t), assembling it on first use. Lets
  // the driver form the kinetic energy (A−a_h)ᵀ M_sheet (A−a_h) directly, avoiding
  // cancellation.
  void ApplySheetMass(const Vector &x, Vector &y)
  {
    EnsureSheetMass();
    M_sheet_->Mult(x, y);
  }

  // Apply the London domain mass M_L, for the kinetic energy AᵀM_L·A of the domains.
  void ApplyLondonMass(const Vector &x, Vector &y);

  // Assemble the right-hand side source term vector for a current source applied on
  // specified excited boundaries.
  void GetCurrentExcitationVector(int idx, Vector &RHS);

  // Assemble flux loop excitation vector for specified flux loop index
  template <ProblemType T>
  void GetFluxExcitationVector(int idx, Vector &RHS, PostOperator<T> &post_op);

  template <ProblemType T>
  void GetFluxExcitationVector(int idx, Vector &RHS, PostOperator<T> &post_op,
                               Vector *boundary_values);

  // True if flux loop idx has a London (finite-λ) film. Its flux excitation uses the
  // shifted London sheet penalty ½∫(1/L_ksq)|A_t − a_h|² (RHS = M_sheet·a_h), and the
  // extracted energy is corrected to |A_t − a_h|².
  bool IsLondonFluxLoop(int idx) const { return london_flux_loops_.contains(idx); }

  // Measure the hole flux ∮A·dl = cᵀA of a candidate London field via the stored fluxoid
  // functional (exact by Stokes; cheaper than integrating B·n over the hole).
  double MeasureLondonHoleFlux(int idx, const Vector &A) const
  {
    return linalg::Dot(GetComm(), london_flux_constraint_.at(idx), A);
  }

  // Build the London cut cohomology generator a_h for the given flux loop.
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
