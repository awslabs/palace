// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP
#define PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP

#include <complex>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class MaterialOperator;
class SurfaceTerms;

namespace mode_assembly
{

using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                         std::unique_ptr<mfem::HypreParMatrix>>;

enum class CoefficientType
{
  CONSTANT,
  OMEGA,
  OMEGA_SQUARED,
  SHIFT,
  // Scale of the frequency-dependent boundary mass term index of SurfaceTerms.
  SURFACE_TERM,
  // Nonlinear system coefficient g(iω) = s²χ(s) of the frequency-dependent material index.
  FREQUENCY_DEPENDENT_PERMITTIVITY,
  // Scale β(ω) of the second-order absorbing boundary term index of SurfaceTerms, for its
  // terms which are not in B, and σ β(ω) for those which are (see ModeOperatorModel).
  SECOND_ORDER_FARFIELD,
  SECOND_ORDER_FARFIELD_SHIFT
};

// One complete, BC-eliminated component of the mode eigenproblem matrix A. Exact assembly
// combines these sparse matrices, while the wave-port ROM projects the same operators.
struct OperatorComponent
{
  CoefficientType type;
  int index;
  std::unique_ptr<mfem::HypreParMatrix> Ar, Ai;
  std::unique_ptr<ComplexWrapperOperator> op;

  OperatorComponent(CoefficientType type, int index,
                    std::unique_ptr<mfem::HypreParMatrix> &&Ar,
                    std::unique_ptr<mfem::HypreParMatrix> &&Ai);
  OperatorComponent(OperatorComponent &&) = default;
  OperatorComponent &operator=(OperatorComponent &&) = default;
  OperatorComponent(const OperatorComponent &) = delete;
  OperatorComponent &operator=(const OperatorComponent &) = delete;
};

// Frequency-parametric mode eigenproblem matrix. The complete fixed components are the
// single source of truth for both exact sparse assembly and wave-port projection.
//
// The eigenproblem (A₀(ω) + kₙ² B(ω)) e = 0 is solved for the shifted matrix A(ω, σ) =
// A₀(ω) - σ B(ω), the sum of the components scaled by their coefficients. Without
// second-order absorbing boundaries, B = [Btt, 0; Btn, 0] does not depend on frequency.
// On the edges of the cross-section, the boundary curl-curl term β (n·∇×E, n·∇×F) of the
// second-order absorbing boundary condition of the 3D solve, for E = (Eₜ + ẑEₙ) e^{-ikₙz}
// and F = (Fₜ + ẑFₙ) e^{ikₙz}, has the surface curl n·∇×E = ∂_τEₙ + ikₙE_τ along the
// tangent τ = ẑ × n of the edge (and n·∇×F = ∂_τFₙ - ikₙF_τ). With ẽₙ = ikₙEₙ, and the
// equation of the out-of-plane component scaled by -ikₙ as for the domain terms, it adds
//   β (kₙ² [M, 0; Gᵀ, 0] - [0, G; 0, D])
// with the boundary forms M = (E_τ, F_τ), G = (∂_τẽₙ, F_τ), and D = (∂_τẽₙ, ∂_τFₙ), with
// the coefficient μ⁻¹c₀. So it adds the component -[0, G; 0, D], with coefficient β(ω), and
// the component -[M, 0; Gᵀ, 0], with coefficient σ β(ω) in A and -β(ω) in B (see
// EvaluateBCoefficient), like the domain terms of Btt, Atn, Ann and Btn.
class ModeOperatorModel
{
private:
  const MaterialOperator &mat_op;
  const SurfaceTerms &surface_terms;
  std::vector<OperatorComponent> components;

  // The lower-left (H1 x ND) blocks -Gᵀ of the components -[M, 0; Gᵀ, 0] of second-order
  // absorbing boundaries, without elimination of essential DOFs, and their term index.
  std::vector<std::pair<int, std::unique_ptr<mfem::HypreParMatrix>>> farfield_lower_blocks;

  // Largest light speed of the materials next to second-order absorbing boundaries.
  double farfield_light_speed_max = 0.0;

public:
  ModeOperatorModel(const FiniteElementSpace &nd_fespace,
                    const FiniteElementSpace &h1_fespace, const MaterialOperator &mat_op,
                    const mfem::Vector *normal, SurfaceTerms &surface_terms,
                    const mfem::HypreParMatrix &Bttr, const mfem::HypreParMatrix *Atnr,
                    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Btnr,
                    const mfem::Array<int> &dbc_tdof_list);

  const std::vector<OperatorComponent> &GetComponents() const { return components; }
  std::complex<double> EvaluateCoefficient(const OperatorComponent &component,
                                           std::complex<double> omega, double sigma) const;
  std::unique_ptr<ComplexOperator> Assemble(std::complex<double> omega, double sigma) const;

  // Coefficient of the component in B(ω) = B₀ + Σ bᵢ(ω) Cᵢ, with the frequency-independent
  // B₀ = [Btt, 0; Btn, 0] (zero for the components which are not in B), whether B depends
  // on frequency, and the assembled B(ω) (MPI collective on the FE space communicator).
  std::complex<double> EvaluateBCoefficient(const OperatorComponent &component,
                                            std::complex<double> omega) const;
  bool HasFrequencyDependentB() const;
  std::unique_ptr<ComplexOperator> AssembleB(const mfem::HypreParMatrix &B0,
                                             std::complex<double> omega) const;

  // Add the lower-left (H1 x ND) block of A(ω, σ) from second-order absorbing boundaries,
  // σ β(ω) (-Gᵀ), without elimination of essential DOFs, to the real and imaginary parts of
  // a matrix (MPI collective on the FE space communicator).
  void AddSecondOrderFarfieldLowerBlock(std::complex<double> omega, double sigma,
                                        std::unique_ptr<mfem::HypreParMatrix> &real,
                                        std::unique_ptr<mfem::HypreParMatrix> &imag) const;

  // In the real-valued approximation Ar + Ai of A(ω, σ), replace the terms of second-order
  // absorbing boundaries, β(ω) (kₙ² [M, 0; Gᵀ, 0] - [0, G; 0, D]) with kₙ² = -σ, by the
  // magnitude |β(ω)| of their coefficient, when the shift is above the passivity limit
  // of the boundaries (MPI collective on the FE space communicator). Their blocks have the
  // signs of the stiffness of the blocks of A, the curl-curl and shift terms of the ND
  // block, the coupling terms, and the negative diffusion of the H1 block, so that with
  // |β| they add to it, while with the imaginary part of β = 1/(2iω), which is negative,
  // they cancel part of it. In the ND block, the first-order condition adds ω/c to the
  // -kₙ²c/(2ω) of the second-order one on the edges, times the same edge mass, for the
  // light speed c of the material there: below the passivity limit kₙ = √2 ω/c, their sum
  // is positive, and Ar + Ai keeps its magnitude. Above it, for the slow modes of a
  // cross-section whose edges are in a material of lower index, the sum is negative, as in
  // the real-valued coarse solve of the 3D problem (SpaceOperator::AssemblePreconditioner),
  // and the solves need many more iterations or stall.
  void
  AddSecondOrderFarfieldRealPreconditioner(std::complex<double> omega, double sigma,
                                           std::unique_ptr<mfem::HypreParMatrix> &P) const;
};

// The inverse permeability of the terms of the mode eigenproblem which come from the
// in-plane part of the curl of the field, (∇×E)ₜ = -n×(∇ₜEₙ + ikₙEₜ) for the direction of
// propagation n: ν = [n×]ᵀ μ⁻¹ [n×], μ⁻¹ rotated by 90 degrees about n, which is μ⁻¹ for a
// permeability isotropic in the plane of the cross-section, but in general has the
// in-plane components of μ⁻¹ across those of μ⁻¹ (for a diagonal μ, the x component of ν
// is 1/μ_y). These are the terms of Btt, Atn, Btn and Ann, and ν is the admittance of the
// in-plane magnetic field: n×Hₜ is ν (kₙEₜ - i∇ₜEₙ)/ω, up to the sign of the convention.
// Without a normal, the inverse permeability is the in-plane 2 x 2 μ⁻¹ of a 2D
// cross-section, and ν is Rᵀμ⁻¹R for the rotation R by 90 degrees, the adjugate of μ⁻¹.
// With the normal of a cross-section in 3D (not necessarily of unit length), it is the
// 3 x 3 μ⁻¹.
mfem::DenseMatrix GetTransverseInvPermeability(const mfem::DenseMatrix &muinv,
                                               const mfem::Vector *normal);
mfem::DenseTensor GetTransverseInvPermeability(const mfem::DenseTensor &muinv,
                                               const mfem::Vector *normal);

// Frequency-independent block matrices Atn = -(ν ∇ₜu, v) and Btt = (ν u, v), with ν of
// GetTransverseInvPermeability for the normal of the cross-section (null for a 2D mesh).
ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal);
ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal);

ComplexHypreParMatrix BuildSystemMatrixA(
    const FiniteElementSpace &nd_fespace, const FiniteElementSpace &h1_fespace,
    const mfem::Array<int> &dbc_tdof_list, const mfem::HypreParMatrix *Attr,
    const mfem::HypreParMatrix *Atti, const mfem::HypreParMatrix *Atnr,
    const mfem::HypreParMatrix *Atni, const mfem::HypreParMatrix *Annr,
    const mfem::HypreParMatrix *Anni, const mfem::HypreParMatrix *shifted_Btnr = nullptr,
    Operator::DiagonalPolicy diag_policy = Operator::DIAG_ONE);

void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en);

// Propagation constant kₙ of a boundary mode from kₙ² = −σ − 1/λ at the frequency ω. At a
// real frequency, kₙ is the principal root, with Re{kₙ} ≥ 0, when Re{kₙ²} ≥ 0, and
// otherwise (below cutoff) the root with Im{kₙ} < 0, for which the mode e^{−ikₙxₙ} decays
// in the +n̂ direction: below cutoff, the principal root of a lossless mode would follow the
// sign of the roundoff in Im{kₙ²}. For a passive cross-section, with Im{kₙ²} ≤ 0, both are
// the same root, with Re{kₙ} ≥ 0 and Im{kₙ} ≤ 0, continuous through cutoff. At a complex
// frequency, kₙ is the analytic continuation of the root at the real frequency Re{ω} along
// Re{ω} + it: for kₙ² = n²ω² − k_c² with real n² and k_c², the mode is below cutoff at
// Re{ω} when Re{kₙ²} + Im{kₙ²} Im{ω} / (2 Re{ω}) < 0. For example, kₙ = nω for a TEM mode
// at any complex ω with Re{ω} > 0, while Re{kₙ²} < 0 for |Im{ω}| > Re{ω}.
std::complex<double> PropagationConstant(std::complex<double> kn2,
                                         std::complex<double> omega = 0.0);

// Ranking key of a boundary mode with propagation constant kₙ relative to the real shift
// target kₙ,target = √(−σ): the complex distance |kₙ − kₙ,target|. The exact and reduced
// mode solves both rank modes by ascending distance, so they select the same mode. Unlike
// |Re{kₙ} − kₙ,target|, the distance accounts for attenuation: it is at least |Im{kₙ}|, so
// a strongly evanescent mode of a lossy cross-section ranks behind the low-loss propagating
// modes below the target (distance < kₙ,target) even when its real part lies close to the
// target. For lossless propagating modes below the target, ascending distance is descending
// Re{kₙ}. The distance is invariant under kₙ → conj(kₙ), so it is the same for the two
// roots ±iα of a lossless mode below cutoff. NaN maps to +∞ so the key always defines a
// strict weak ordering.
double TargetDistance(std::complex<double> kn, double kn_target);

}  // namespace mode_assembly
}  // namespace palace

#endif  // PALACE_MODELS_MODE_OPERATOR_ASSEMBLY_HPP
