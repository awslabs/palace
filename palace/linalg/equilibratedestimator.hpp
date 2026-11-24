// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_EQUILIBRATED_ESTIMATOR_HPP
#define PALACE_LINALG_EQUILIBRATED_ESTIMATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "utils/configfile.hpp"

namespace palace
{

class ErrorIndicator;
class MaterialOperator;

namespace equilibrated
{

// Reference tables of the basis functions of the flux and primal spaces and of the
// coordinate map at the quadrature points of one element geometry, for cells with an
// affine coordinate map (the Jacobian evaluated once) or not.
struct QuadratureTables
{
  const mfem::IntegrationRule *rule = nullptr;
  std::vector<mfem::DenseMatrix> flux, map, grad;
  std::vector<mfem::Vector> div, shape;
};

// The cell loops of the estimator: the quadrature rule of order 2p + 2 (+ 2 (q - 1) on
// cells with a curved coordinate map of order q), exact on affine cells for every integrand
// of the bound, or of a given order instead of 2p + 2, and the reference tables of the flux
// and primal spaces.
class CellQuadrature
{
private:
  const FiniteElementSpace &flux, &primal;
  int degree;
  std::vector<bool> affine;
  mutable std::map<std::pair<int, bool>, QuadratureTables> tables;

public:
  CellQuadrature(const FiniteElementSpace &flux, const FiniteElementSpace &primal,
                 int degree = -1);

  int Order(bool curved = true) const;
  bool IsAffine(int e) const { return affine[e]; }
  const QuadratureTables &GetTables(int e, mfem::IsoparametricTransformation &T) const;
};

}  // namespace equilibrated

// The global curl correction of an equilibrated electrostatic flux. Two fluxes in the
// Raviart-Thomas space RT(k) with the same divergence and the same normal trace on the
// Neumann faces differ by a discretely divergence-free field, the curl of a field of the
// Nédélec space ND(k + 1) (the exact sequence H1(k + 1) → ND(k + 1) → RT(k) → L2(k), MFEM
// orders), broken across the interior Dirichlet boundaries as the flux space, whose
// tangential trace vanishes on the Neumann faces (and on interior faces with a Neumann
// attribute), up to finitely many harmonic fields (the capacitive fields, see
// EquilibratedFluxErrorEstimator). For a flux D₀ and the flux D_v = -ε ∇v of a potential
// v, D = D₀ + ∇ × A minimizes
//   η²(A) = ‖ε^{-1/2} (D₀ + ∇ × A - D_v)‖² = η₀² - 2 bᵀ A + Aᵀ K A,
// K = (ε⁻¹ ∇ × ., ∇ × .), b = -(ε⁻¹ (D₀ - D_v), ∇ × .): the curl-curl problem of the
// magnetostatic solver with ε⁻¹ in place of μ⁻¹, singular and consistent (b vanishes on the
// gradients, since v is constant on each Dirichlet boundary), solved by CG from zero with
// its preconditioner, built from the configuration of a magnetostatic solve (geometric
// multigrid with Chebyshev smoothing on ND(k + 1) and ND(1) and AMS for singular operators
// on ND(1), with direct interpolation in its auxiliary space hierarchies). CG minimizes η²
// over its Krylov space, so every iterate is admissible (the curl changes neither the
// divergence nor the normal trace on the Neumann faces) and η decreases monotonically. It
// stops at the relative residual "EstimatorTol", by the least-squares test of CgSolver once
// η is estimated to be within "EstimatorErrorTol" (relative) of the bound of the converged
// correction, or once the decrease of η² is below its rounding error. The bound never
// depends on how far CG is run: Correct keeps the uncorrected flux if the corrected one,
// with the divergence of its curl, has the larger bound (see Correct): the gradient part of
// the iterate, which rounding grows, has a curl exact only to rounding relative to it, so
// that the flux could lose some divergence.
class CurlCorrectionSolver
{
private:
  const MaterialOperator &mat_op;
  const FiniteElementSpace &primal, &flux;
  const equilibrated::CellQuadrature &quadrature;
  const std::vector<int> dirichlet_attr;

  // The spaces ND(1), ND(k + 1) (the finest the caller's if it has that order), H1(1) for
  // AMS, broken across the interior Dirichlet boundaries when the flux space is.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  std::unique_ptr<mfem::H1_FECollection> h1_fec;
  FiniteElementSpaceHierarchy nd_spaces, h1_spaces, nd_broken, h1_broken;
  const FiniteElementSpace *shared = nullptr;  // The finest level, if shared_nd.
  bool broken = false;
  std::vector<mfem::Array<int>> essential;  // Per level.

  // K on every level (the finest is that of CG).
  std::unique_ptr<Operator> K;
  std::unique_ptr<KspSolver> ksp;
  double ls_tol = 0.0;  // Of the least-squares test, from EstimatorErrorTol.

  // The reference curl matrix of each element geometry: ∇ × Nᵢ = Σⱼ Cⱼᵢ φⱼ exactly.
  mutable std::map<int, mfem::DenseMatrix> curl_matrix;

  FiniteElementSpaceHierarchy &NDSpaces() { return broken ? nd_broken : nd_spaces; }
  const FiniteElementSpace &ND() const
  {
    return broken ? nd_broken.GetFinestFESpace() : nd_spaces.GetFinestFESpace();
  }
  const mfem::DenseMatrix &GetCurlMatrix(int e) const;

  // The dual ND L-vector -(ε⁻¹ (D₀ - D_v), ∇ × Nᵢ) of the flux and companion L-vectors;
  // returns η₀² (local).
  double AssembleRHS(const Vector &flux_l, const Vector &companion_l, Vector &b_l) const;

  // The RT L-vector of the curl of the ND L-vector a_l.
  Vector Curl(const Vector &a_l) const;

public:
  CurlCorrectionSolver(const MaterialOperator &mat_op, FiniteElementSpace &primal,
                       const FiniteElementSpace &flux,
                       const std::vector<int> &dirichlet_attr,
                       FiniteElementSpace *shared_nd,
                       const config::LinearSolverData &linear,
                       const equilibrated::CellQuadrature &quadrature);

  // Releases the discrete gradient of the coarsest level cached in a shared finest level
  // (one level at order one), whose auxiliary space is that of this solver. Otherwise a
  // later solver whose auxiliary space reuses the address on some processes but not on
  // others would find the gradient of the freed space there and build it, collectively,
  // only on the others.
  ~CurlCorrectionSolver();

  struct Result
  {
    // η² before the correction and after each CG iteration (from the decrease of the
    // objective), the number of iterations, and the Euclidean norm of the final iterate A
    // (true DOFs), whose gradient part the curl does not see.
    std::vector<double> eta2;
    int iterations = 0;
    double norm = 0.0;
    // Whether the curl was added (see Correct).
    bool accepted = true;
  };

  // Correct the flux (an L-vector of the flux space) in place, for the companion potential
  // (an L-vector of the primal space) with Dirichlet range range[0] ≤ v ≤ range[1] on the
  // Dirichlet boundary, with the least-squares tolerance of the constructor or the given
  // one (if not negative). The curl of the final iterate is added as by AddIfNotWorse, so
  // that the result never depends on how CG behaves. Collective.
  Result Correct(Vector &flux_l, const Vector &companion_l, const double range[2],
                 double tol = -1.0) const;

  // The bound L = (η + (η² + 4 b)^{1/2}) / 2 of a flux (an L-vector of the flux space) for
  // the companion potential with Dirichlet range range[0] ≤ v ≤ range[1], evaluated as by
  // EquilibratedFluxErrorEstimator (η and the balance term b of its cells and Neumann
  // faces). Collective.
  double Bound(const Vector &flux_l, const Vector &companion_l,
               const double range[2]) const;

  // Add the curl L-vector to the flux only if the assembled corrected flux has no larger
  // bound (but for rounding, 1e-12 relative) than the flux; returns whether it was added.
  // Collective.
  bool AddIfNotWorse(Vector &flux_l, const Vector &curl_l, const Vector &companion_l,
                     const double range[2]) const;

  // The tolerance of the least-squares test from EstimatorErrorTol (0 if it is 0).
  double GetLeastSquaresTol() const { return ls_tol; }

  HYPRE_BigInt GlobalTrueVSize() const { return ND().GlobalTrueVSize(); }
};

// Guaranteed upper bound of the energy error of an electrostatic solution: the flux D of
// the Raviart-Thomas space of L2ZZ's recovered flux, RT(p - 1) for solution order p (or
// RT(p), see EstimatorFluxOrder), with the divergence constraint of the exact flux added:
// div D = 0, D ⋅ n = 0 on the Neumann faces. Then for the computed potential u_h, by
// Prager and Synge, |||u - u_h||| ≤ ‖ε^{-1/2} (D + ε ∇u_h)‖, |||v|||² = (ε ∇v, ∇v). D is
// built in three steps:
//  1. An exactly equilibrated RT(0) flux from lowest-order vertex patches (see
//     fem::VertexPatchFlux), embedded exactly into RT(k).
//  2. The global ND(k + 1) curl correction (see CurlCorrectionSolver), for which every CG
//     iterate is a guaranteed bound, decreasing monotonically.
//  3. The capacitive step: the net fluxes into the connected components of the Dirichlet
//     boundary are not changed by a curl, and the divergence-free fields with a zero
//     normal trace on the Neumann faces are the curls and one harmonic field per
//     component, less one. The corrected fluxes D⁽ʲ⁾ of the N terminal solves (unit
//     voltage on terminal j), and the fluxes of the order-one potentials which complete
//     their Dirichlet data to a basis of the potentials constant on each component (a
//     ground, or a terminal, of several parts), corrected as far as their weight in the
//     combination requires (see AddErrorIndicators), are all admissible, so for solve i the
//     best combination D = Σⱼ αⱼ D⁽ʲ⁾ is admissible: with Gⱼₗ = (ε⁻¹ D⁽ʲ⁾, D⁽ˡ⁾) and Qᵢⱼ =
//     (ε⁻¹ D⁽ʲ⁾, D_v⁽ⁱ⁾), α = G⁻¹ Qᵢ.ᵀ and η² = (C_h)ᵢᵢ - Qᵢ. G⁻¹ Qᵢ.ᵀ, with C_h = (ε
//     ∇v⁽ⁱ⁾, ∇v⁽ʲ⁾) the capacitance matrix of the companions.
// When converged, the flux is the dual mixed RT(k) x L2(k) solution. The bound is
// evaluated as for a locally equilibrated flux, accounting for what remains of the
// divergence (the algebraic error of the solve) and for the Dirichlet data. For the
// companion v, which takes the intended Dirichlet values at the essential DOFs, with
// lifting β = |||v - u_h|||, D_v = -ε ∇v and e = u - v:
//   |||e|||² ≤ η |||e||| + b,  η² = Σ_K ‖ε^{-1/2} (D - D_v)‖²_K,
//   b = Σ_K ‖div D‖_K (R |K|^{1/2} + ‖v - c‖_K) + Σ_{F ∈ Γ_N} ‖D ⋅ n‖_F (R |F|^{1/2} +
//       ‖v - c‖_F),
// by the maximum principle (the Dirichlet range [c - R, c + R] bounds u), so |||u - u_h|||
// ≤ L + β, L = (η + (η² + 4 b)^{1/2}) / 2. The element indicators q_K² = (L + β) (m_K / L +
// β_K² / β), m_K = (L / η) η_K² + b_K, sum to (L + β)², and are normalized by a guaranteed
// lower bound of the exact energy norm, |||u|||² ≥ ((2 U_h)^{1/2} - β)₊² - L², so that the
// norm of the indicators of a solve is a guaranteed upper bound of its relative error
// |||u - u_h||| / |||u||| (see RelativeBound), and an adaptive refinement with
// Refinement.Tol stops only when the root mean square over the solves of the true relative
// errors is below it. The capacitance matrix C of the problem as meshed is C_h - E, with E
// the Gram matrix of the errors of the companions, bounded entry by entry from the bounds
// of the solves and of combinations of them (see CapacitanceBrackets), at any tolerance of
// the terminal solves. Integrals use quadrature, exact on affine cells.
class EquilibratedFluxErrorEstimator
{
public:
  struct Summary
  {
    // The bound L + β of |||u - u_h||| and its parts, and the bound of the relative error.
    double bound = 0.0, eta = 0.0, balance = 0.0, lifting = 0.0, relative = 0.0;
    // ‖ε^{-1/2} (D₀ - D_v)‖ of the vertex patch flux, before the curl correction, and the
    // number of CG iterations of the correction.
    double initial = 0.0;
    int iterations = 0;
    // ‖div D‖ and ‖D ⋅ n‖ on the Neumann faces, which vanish for an exactly equilibrated
    // flux (they enter the balance term).
    double divergence = 0.0, neumann = 0.0;
  };

  // Guaranteed upper bound of |||u - u_h||| / |||u||| from the bound A = L + β of
  // |||u - u_h|||, its part L ≥ |||u - v||| and lifting β, and the computed energy U_h: A /
  // (((2 U_h)^{1/2} - β)₊² - L²)^{1/2}; 0 if A = 0, and +∞ if that lower bound of |||u|||²
  // is not positive (the bound certifies no relative accuracy).
  static double RelativeBound(double L, double lifting, double energy);

  // Value of an infinite relative bound in the indicators, so that they stay finite for
  // marking and an adaptive refinement does not stop for any Refinement.Tol below it.
  static constexpr double uncertified = 1.0e8;

private:
  const MaterialOperator &mat_op;
  FiniteElementSpace &h1_fespace;
  FiniteElementSpace *nd_fespace;
  config::LinearSolverData linear;
  std::vector<int> dirichlet_attr;
  mfem::Array<int> essential;

  // The RT(k) flux space (the caller's if it has that order) and its view broken across the
  // interior Dirichlet boundaries, if there are any (the spaces of the patches, the
  // order-one hats and RT(0), are built for them only).
  mfem::H1_FECollection hat_fec;
  mfem::RT_FECollection rt0_fec;
  std::unique_ptr<mfem::RT_FECollection> flux_fec;
  std::unique_ptr<FiniteElementSpace> flux_owned, flux_broken;
  FiniteElementSpace *flux_unbroken;
  bool split = false;
  std::unique_ptr<equilibrated::CellQuadrature> quadrature, low_quadrature;

  // Reference RT(0) → RT(k) embedding of each element geometry.
  mutable std::map<int, mfem::DenseMatrix> embedding;

  // The capacitance matrix of the companions and the entrywise brackets of the exact
  // capacitance and mutual capacitance matrices (see CapacitanceBrackets).
  mfem::DenseMatrix C_h, C_lower, C_upper, Cm_lower, Cm_upper;

  // The connected components of the Dirichlet boundary, the number of fields of the
  // capacitive step besides those of the solves and of those corrected, in the last
  // AddErrorIndicators.
  int num_components = 0, num_harmonic = 0, num_harmonic_corrected = 0;

  const FiniteElementSpace &Flux() const
  {
    return flux_broken ? *flux_broken : *flux_unbroken;
  }

  // The flux L-vector of Flux() of the true-DOF vector x of the RT(0) space rt0.
  Vector Embed(const FiniteElementSpace &rt0, const Vector &x) const;

  // The Gram matrices of the fluxes and companions (L-vectors): G of the fluxes, Q of the
  // companions and the fluxes, and C_h of the companions (global).
  void Gram(const std::vector<Vector> &fluxes, const std::vector<Vector> &companions,
            mfem::DenseMatrix &G, mfem::DenseMatrix &Q, mfem::DenseMatrix &Ch) const;

  // The rows of the Gram matrices of the fields fₛ of the order-one potentials wₛ (RT(0)
  // true DOFs of the patch space rt0, true DOFs of hats), appended to G and Q of the fluxes
  // and companions (L-vectors) of the solves (n x n on input, m x m and n x m on output, m
  // = n + the number of fields), and the mismatch ‖ε^{-1/2} (fₛ + ε ∇wₛ)‖² of each field,
  // with the lower-order quadrature (global).
  void HarmonicGram(const std::vector<Vector> &fluxes,
                    const std::vector<Vector> &companions, const FiniteElementSpace &rt0,
                    const FiniteElementSpace &hats, const std::vector<Vector> &f,
                    const std::vector<Vector> &w, mfem::DenseMatrix &G,
                    mfem::DenseMatrix &Q, std::vector<double> &mismatch) const;

  // The entrywise brackets of C and C_m from the final fluxes and companions (L-vectors),
  // the intended values at the essential DOFs, the summaries of the solves and the Gram
  // matrix H of the mismatches of the fluxes. Collective.
  void CapacitanceBrackets(const std::vector<Vector> &fluxes,
                           const std::vector<Vector> &companions,
                           const std::vector<Vector> &G,
                           const std::vector<Summary> &summaries,
                           const mfem::DenseMatrix &H);

  // The bounds and element indicators of the fluxes (L-vectors) for the companions
  // (L-vectors and true DOFs), liftings (L-vectors), Dirichlet ranges and energies.
  std::vector<Summary>
  Evaluate(const std::vector<Vector> &fluxes, const std::vector<Vector> &companions,
           const std::vector<Vector> &liftings, const std::vector<double> &range,
           const std::vector<double> &energies, ErrorIndicator &indicator) const;

public:
  // The flux space is RT(p - 1) (flux_order_offset -1, the recovery space of L2ZZ) or RT(p)
  // (0), for solution order p; rt_fespace and nd_fespace, RT(p - 1) and ND(p) spaces of the
  // caller on the mesh, are used for the flux and for the finest level of the correction
  // when they are of these orders. They must outlive the estimator. The linear solver data
  // are the resolved ones of the run: the correction uses its multigrid settings, and
  // EstimatorTol, EstimatorErrorTol and EstimatorMaxIts.
  EquilibratedFluxErrorEstimator(const MaterialOperator &mat_op,
                                 FiniteElementSpace &h1_fespace,
                                 FiniteElementSpace &rt_fespace,
                                 FiniteElementSpace &nd_fespace,
                                 const mfem::Array<int> &dirichlet_attr,
                                 int flux_order_offset,
                                 const config::LinearSolverData &linear);
  ~EquilibratedFluxErrorEstimator();

  // For each solve t: V[t] is the computed potential (true DOFs), G[t] holds the intended
  // values at the essential true DOFs of the Dirichlet attributes (other entries are
  // ignored) and energies[t] is the field energy U_h. With one unit voltage per solve, the
  // solves are those of the terminals, which the capacitive step combines. The indicators
  // of the solves are added in order, each normalized so that its norm is the relative
  // bound (at most uncertified). Returns the summaries in the same order. Collective.
  std::vector<Summary> AddErrorIndicators(const std::vector<Vector> &V,
                                          const std::vector<Vector> &G,
                                          const std::vector<double> &energies,
                                          ErrorIndicator &indicator);

  // The capacitance matrix C_h of the companions and the entrywise brackets lower ≤ C ≤
  // upper of the exact capacitance matrix and of the mutual capacitance matrix C_m of the
  // last AddErrorIndicators (in its solve order).
  const auto &GetCapacitance() const { return C_h; }
  const auto &GetCapacitanceLower() const { return C_lower; }
  const auto &GetCapacitanceUpper() const { return C_upper; }
  const auto &GetMutualCapacitanceLower() const { return Cm_lower; }
  const auto &GetMutualCapacitanceUpper() const { return Cm_upper; }

  // The number of connected components of the Dirichlet boundary and of the fields which
  // the capacitive step of the last AddErrorIndicators added to those of the solves.
  int GetNumDirichletComponents() const { return num_components; }
  int GetNumHarmonicFields() const { return num_harmonic; }
  int GetNumCorrectedHarmonicFields() const { return num_harmonic_corrected; }
};

}  // namespace palace

#endif  // PALACE_LINALG_EQUILIBRATED_ESTIMATOR_HPP
