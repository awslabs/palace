// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <numbers>
#include <string>
#include <vector>
#include <mfem.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/permittivitymodel.hpp"
#include "utils/rationalfit.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const bool has_C = (C != nullptr);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  // Check if there are nonlinear terms and, if so, setup interpolation operator.
  auto funcA2 = [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
  };
  // As funcA2, plus the matrix-free wave-port modal correction W. Used by the SLP and
  // Quasi-Newton solves (the polynomial seed keeps funcA2, since W cannot enter
  // BuildParSumOperator).
  auto funcA2_full =
      [&space_op](std::complex<double> lambda) -> std::unique_ptr<ComplexOperator>
  {
    const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);  // ω = λ/i
    return space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
  };
  auto funcP = [&space_op](std::complex<double> a0, std::complex<double> a1,
                           std::complex<double> a2,
                           std::complex<double> a3) -> std::unique_ptr<ComplexOperator>
  { return space_op.GetPreconditionerMatrix<ComplexOperator>(a0, a1, a2, a3); };
  const double target = iodata.solver.eigenmode.target;
  auto A2 = funcA2(1i * target);
  // Known poles of A2(λ), which the SLP nonlinear eigensolver clears from its local model.
  const auto a2_poles = space_op.GetMaterialOp().GetFrequencyDependentPermittivityPoles();
  // A2 is non-null whenever a frequency-dependent material is present, even if its
  // contribution cancels at the target, so it alone selects the nonlinear solver path.
  bool has_A2 = (A2 != nullptr);

  // Extend K, C, M operators with interpolated A2 operator.
  // K' = K + A2_0, C' = C + A2_1, M' = M + A2_2
  std::unique_ptr<ComplexOperator> Kp, Cp, Mp;
  std::unique_ptr<Interpolation> interp_op;
  std::unique_ptr<ComplexOperator> A2_0, A2_1, A2_2;
  NonlinearEigenSolver nonlinear_type = iodata.solver.eigenmode.nonlinear_type;
  // SLP is only realized through SLEPc's NEP path. Resolve it against the selected backend
  // here, before the HYBRID interpolation below: an ARPACK backend (or a build without
  // SLEPc) has no SLP solver, and switching later would leave Kp/Cp/Mp unbuilt and crash.
  bool slp_available = true;
#if !defined(PALACE_WITH_SLEPC)
  slp_available = false;
#endif
  if (iodata.solver.eigenmode.type == EigenSolverBackend::ARPACK)
  {
    slp_available = false;
  }
  if (nonlinear_type == NonlinearEigenSolver::SLP && !slp_available)
  {
    Mpi::Warning("SLP nonlinear eigensolver requires the SLEPc backend, using Hybrid!\n");
    nonlinear_type = NonlinearEigenSolver::HYBRID;
  }
  if (nonlinear_type == NonlinearEigenSolver::SLP && !has_A2)
  {
    // The SLP nonlinear eigensolver requires a frequency-dependent term A2(λ) (wave ports
    // or absorbing boundaries). Without one the NEP function/Jacobian shells dereference a
    // null A2; fall back to the standard linear/quadratic eigensolver instead.
    Mpi::Warning(
        "SLP nonlinear eigensolver requires a nonlinear system term (wave ports or "
        "absorbing boundaries), none found; using a linear eigensolver!\n");
    nonlinear_type = NonlinearEigenSolver::HYBRID;
  }

  // Every term of A2(λ) is a scalar function of λ times a fixed matrix, f(λ) B. When f is
  // rational, the linearized eigensolvers represent the term exactly, with one auxiliary
  // block per pole (EigenvalueSolver::AddRationalTerm). Problems whose terms are all
  // rational are then linear eigenvalue problems: there is no polynomial seed and no
  // nonlinear refinement. The rational terms are:
  //  - frequency-dependent permittivity pole terms gₘ(λ) Bₘ, gₘ(λ) = Σₖ rₖ λ² / (λ - pₖ),
  //  - second-order absorbing boundaries next to constant materials, -0.5 / λ M_ff,
  //  - rational impedance boundaries g(λ) M_b, g(λ) = λ D(λ) / N(λ).
  // The absorbing boundaries next to frequency-dependent materials (with √(ε(λ) / ε∞)) and
  // the finite conductivity boundaries of zero thickness (with √λ) are fitted with
  // rational functions in a window of the search (utils::FitScalar), valid outside a strip
  // around their branch cuts. Their reported errors are those of the complete operator.
  // Djordjevic-Sarkar terms and wave ports keep the Hybrid seed and the Quasi-Newton
  // refinement. RefineNonlinear has no effect here: the terms are exact, or fitted to the
  // tolerance of the eigensolver, and there is nothing to refine. With PEPLinear = false,
  // SLEPc's polynomial eigensolver linearizes the problem itself, without the auxiliary
  // blocks.
  struct RationalA2Term
  {
    std::unique_ptr<ComplexOperator> B;
    std::array<std::complex<double>, 3> poly = {0.0, 0.0, 0.0};
    std::vector<std::complex<double>> poles, residues;
    std::complex<double> Eval(std::complex<double> l) const
    {
      std::complex<double> f = poly[0] + l * (poly[1] + l * poly[2]);
      for (std::size_t j = 0; j < poles.size(); j++)
      {
        f += residues[j] / (l - poles[j]);
      }
      return f;
    }
    // For a fitted term, its complete scalar function, a name, and its fit.
    std::function<std::complex<double>(std::complex<double>)> f;
    std::string name;
    utils::ScalarFit fit;
  };
  std::vector<RationalA2Term> a2_terms;
  const utils::FitWindow window = {target, iodata.solver.eigenmode.target_upper - target};
  bool fitted = false;
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID &&
      (iodata.solver.eigenmode.type == EigenSolverBackend::ARPACK ||
       iodata.solver.eigenmode.pep_linear))
  {
    const auto &mat_op = space_op.GetMaterialOp();
    const auto &farfield_op = space_op.GetFarfieldOp();
    const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
    const auto &surf_sigma_op = space_op.GetSurfaceConductivityOp();
    bool rational = space_op.GetWavePortOp().Size() == 0;
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (mat_op.HasFrequencyDependentPermittivitySupport(m) &&
          mat_op.HasDjordjevicSarkarPermittivity(m))
      {
        rational = false;
      }
    }
    for (std::size_t m = 0;
         rational && m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (!mat_op.HasFrequencyDependentPermittivitySupport(m))
      {
        continue;  // No contribution (SpaceOperator::GetExtraSystemMatrix)
      }
      auto B = space_op.GetFrequencyDependentPermittivityMassMatrix(m, Operator::DIAG_ZERO);
      if (!B)
      {
        continue;  // No nonzero pole (exact zero poles are in the conductivity of C)
      }
      // A2 holds gₘ(λ)·Bₘ, and the operator is on the imaginary slot: it applies i·Bₘ.
      // With c = -i r: c λ² / (λ - p) = c p + c λ + c p² / (λ - p).
      auto &term = a2_terms.emplace_back();
      for (const auto &pole : mat_op.GetPermittivityPoleTerms(m))
      {
        const std::complex<double> c = -1i * pole.residue;
        term.poly[0] += c * pole.pole;
        term.poly[1] += c;
        term.poles.push_back(pole.pole);
        term.residues.push_back(c * pole.pole * pole.pole);
      }
      term.B = std::move(B);
    }
    if (rational && farfield_op.GetOrder() >= 2)
    {
      if (auto M_ff = space_op.GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(
              Operator::DIAG_ZERO))
      {
        auto &term = a2_terms.emplace_back();
        term.B = std::move(M_ff);
        term.poles = {0.0};
        term.residues = {-0.5};
      }
    }
    for (int idx = 0; rational && idx < surf_rz_op.GetNumBoundaries(); idx++)
    {
      RationalA2Term term;
      if (!surf_rz_op.GetRobinPoleResidue(idx, term.poly, term.poles, term.residues))
      {
        rational = false;
        break;
      }
      term.B = space_op.GetRationalImpedanceBoundaryMassMatrix<ComplexOperator>(
          idx, Operator::DIAG_ZERO);
      if (term.B)
      {
        a2_terms.push_back(std::move(term));
      }
    }

    // The fitted terms f(λ) = α λᵐ g(λ), with g fitted, and the points of the branch cuts
    // of g as ω = λ / i: of √(ε / ε∞) for the absorbing boundaries next to material, of √λ
    // (λ = 0 and the negative real axis of ω) for the conductivity. The fit tolerance is
    // that of the eigenvalue solve: a relative error ε of a term changes the backward
    // error of the complete operator by at most about ε.
    std::string unfitted;
    const double fit_tol = std::clamp(iodata.solver.eigenmode.tol, 1.0e-12, 1.0e-6);
    auto AddFitted = [&](std::unique_ptr<ComplexOperator> &&B, std::string name,
                         std::function<std::complex<double>(std::complex<double>)> g,
                         std::complex<double> alpha, int m,
                         const std::vector<std::complex<double>> &cut)
    {
      if (!B || !rational)
      {
        return;
      }
      auto &term = a2_terms.emplace_back();
      term.fit = utils::FitScalar(g, alpha, m, cut, window, fit_tol);
      term.poly = term.fit.poly;
      term.poles = term.fit.poles;
      term.residues = term.fit.residues;
      term.f = [g, alpha, m](std::complex<double> l)
      { return alpha * std::pow(l, m) * g(l); };
      term.name = std::move(name);
      term.B = std::move(B);
      rational = term.fit.success;
      unfitted = rational ? unfitted : term.name;
      fitted = true;
    };
    for (std::size_t k = 0; rational && k < farfield_op.NumDispersiveMaterials(); k++)
    {
      std::vector<std::complex<double>> p, r;
      const auto m = farfield_op.GetDispersiveMaterial(k);
      if (!mat_op.GetFrequencyDependentWavePermittivityPoles(m, p, r))
      {
        rational = false;  // Djordjevic-Sarkar
        break;
      }
      {
        // The pole-residue form is the wave permittivity of the material.
        const std::complex<double> l(-0.1 * target, 1.3 * target);
        std::complex<double> e = 1.0;
        for (std::size_t j = 0; j < p.size(); j++)
        {
          e += r[j] / (l - p[j]);
        }
        const auto e_ref = mat_op.EvaluateFrequencyDependentWavePermittivity(m, l) /
                           mat_op.GetFrequencyDependentPermittivityHighFrequency(m);
        MFEM_VERIFY(std::abs(e - e_ref) <= 1.0e-10 * std::abs(e_ref),
                    "The poles of the wave permittivity of material "
                        << m << " do not reproduce it!");
      }
      const auto cut = utils::SquareRootCut(p, r, window.radius / 4096.0);
      auto scale = [&farfield_op, k](std::complex<double> l)
      { return farfield_op.GetDispersiveScale(k, l / 1i); };
      auto inv_scale = [&farfield_op, k](std::complex<double> l)
      { return 1.0 / farfield_op.GetDispersiveScale(k, l / 1i); };
      const auto name = fmt::format("absorbing boundary next to material {:d}", m + 1);
      // The scalars ω r, 0.5 / (ω r), ω / r, and 1 / r of the damping, curl-curl, and
      // frequency-scaled Floquet matrices, with ω = -i λ.
      AddFitted(space_op.GetFarfieldDispersiveDampingMatrix(k, Operator::DIAG_ZERO), name,
                scale, -1i, 1, cut);
      AddFitted(space_op.GetFarfieldDispersiveCurlCurlMatrix(k, Operator::DIAG_ZERO), name,
                inv_scale, 0.5i, -1, cut);
      AddFitted(space_op.GetFarfieldDispersiveDampingFloquetMatrix(k, Operator::DIAG_ZERO),
                name, inv_scale, -1i, 1, cut);
      AddFitted(
          space_op.GetFarfieldDispersiveStiffnessFloquetMatrix(k, Operator::DIAG_ZERO),
          name, inv_scale, 1.0, 0, cut);
    }
    for (std::size_t g = 0; rational && g < surf_sigma_op.Size(); g++)
    {
      if (!surf_sigma_op.IsActive(g))
      {
        continue;
      }
      if (surf_sigma_op.GetThickness(g) > 0.0)
      {
        rational = false;  // The thickness correction has other singularities
        break;
      }
      std::vector<std::complex<double>> cut;
      for (int k = 0; k <= 64; k++)
      {
        cut.push_back(-window.radius * k / 16.0);
      }
      // The matrix applies i A_σ, so the scalar of i ω / Z(ω) A_σ is fitted as λ g(λ).
      AddFitted(
          space_op.GetSurfaceConductivityBoundaryMatrix<ComplexOperator>(
              static_cast<int>(g), Operator::DIAG_ZERO),
          fmt::format("finite conductivity boundary {:d}", g + 1),
          [&surf_sigma_op, g](std::complex<double> l)
          { return surf_sigma_op.EvaluateScalar(g, l / 1i) / (1i * l); }, 1.0, 1, cut);
    }
    if (!rational)
    {
      if (!unfitted.empty())
      {
        Mpi::Print(" The {} has no rational fit with at most 24 poles in the window; using "
                   "the Hybrid nonlinear eigensolver\n",
                   unfitted);
      }
      a2_terms.clear();
      fitted = false;
    }
    if (!a2_terms.empty())
    {
      // The terms must be the complete A2(λ): any other contribution is an error here, not
      // an approximation. Checked at one complex λ, with the complete scalar functions of
      // the fitted terms.
      const std::complex<double> l(-0.01 * target, 1.05 * target);
      ComplexVector v, r, r_terms;
      space_op.GetRandomInitialVector(v, 2);
      r.SetSize(v.Size());
      r.UseDevice(true);
      r_terms.SetSize(v.Size());
      r_terms.UseDevice(true);
      funcA2(l)->Mult(v, r);
      r_terms = 0.0;
      for (const auto &term : a2_terms)
      {
        term.B->AddMult(v, r_terms, term.f ? term.f(l) : term.Eval(l));
      }
      const double norm = linalg::Norml2(space_op.GetComm(), r);
      r_terms.AXPY(-1.0, r);
      const double err = linalg::Norml2(space_op.GetComm(), r_terms);
      MFEM_VERIFY(err <= 1.0e-10 * norm,
                  "The rational terms do not reproduce the frequency-dependent operator "
                  "(relative difference "
                      << err / norm << ")!");
    }
  }
  // The fits, and the rule where they are not valid: each strip around a branch cut in the
  // window is named in a warning, and so is each mode found where a fit is not valid (see
  // below). Modes in a strip may be missing, or wrong; the errors reported for them are
  // those of the complete operator.
  auto ToGHz = [this](std::complex<double> omega)
  {
    return iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
           (2.0 * std::numbers::pi);
  };
  for (std::size_t t = 0; t < a2_terms.size(); t++)
  {
    const auto &term = a2_terms[t];
    if (!term.f)
    {
      continue;
    }
    const auto &fit = term.fit;
    Mpi::Print(" Rational fit of the {}: {:d} poles from {:d} samples, error {:.2e} on "
               "{:d} check points\n",
               term.name, fit.poles.size(), fit.num_samples, fit.error, fit.num_checks);
    if (fit.strip == 0.0 ||
        std::ranges::any_of(a2_terms.begin() + t + 1, a2_terms.end(),
                            [&](const auto &other) { return other.name == term.name; }))
    {
      continue;
    }
    // One warning for each boundary, with the widest strip of its fits.
    double strip = 0.0;
    for (const auto &other : a2_terms)
    {
      strip = (other.name == term.name) ? std::max(strip, other.fit.strip) : strip;
    }
    double re[2] = {mfem::infinity(), -mfem::infinity()}, im[2] = {re[0], re[1]};
    for (const auto w : fit.singular)
    {
      if (window.Contains(w) ||
          (window.Contains(w + strip) || window.Contains(w - strip) ||
           window.Contains(w + 1i * strip) || window.Contains(w - 1i * strip)))
      {
        const auto f = ToGHz(w);
        re[0] = std::min(re[0], f.real());
        re[1] = std::max(re[1], f.real());
        im[0] = std::min(im[0], f.imag());
        im[1] = std::max(im[1], f.imag());
      }
    }
    if (re[0] > re[1])
    {
      continue;  // The strip is outside the window
    }
    Mpi::Warning("The {} has a branch cut in the window of the eigenvalue search (Target "
                 "to TargetUpper), at Re{{f}} = {:.4f} to {:.4f} GHz, Im{{f}} = {:.4f} to "
                 "{:.4f} GHz. Its rational fit is not valid within {:.3e} GHz of the cut: "
                 "modes there may be missing or wrong, and their errors are those of the "
                 "complete operator.\n",
                 term.name, re[0], re[1], im[0], im[1], ToGHz(strip).real());
  }
  const bool exact_poles = !a2_terms.empty();
  const bool hybrid_seed =
      has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID && !exact_poles;
  if (hybrid_seed || nonlinear_type == NonlinearEigenSolver::SLP)
  {
    // The nonlinear eigensolver can miss modes or stop near a material resonance, where
    // A2(λ) changes on the scale of its half-width (see the guide). The search window is
    // not changed: that would change the results.
    const auto &mat_op = space_op.GetMaterialOp();
    const double target_max = iodata.solver.eigenmode.target_upper;
    const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                            (2.0 * std::numbers::pi);
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      auto f = mat_op.GetPermittivityResonances(m, target, target_max, 3.0);
      std::ranges::for_each(f, [unit_GHz](double &f_m) { f_m *= unit_GHz; });
      if (mat_op.HasFrequencyDependentPermittivitySupport(m) && !f.empty())
      {
        Mpi::Warning("The frequency-dependent Permittivity of material {:d} has a "
                     "resonance at {:.6g} GHz, within 3 half-widths of the search window "
                     "from Target to TargetUpper, {:.6g} to {:.6g} GHz!\nThe nonlinear "
                     "eigensolver can miss modes or stop near it.\n",
                     m + 1, fmt::join(f, ", "), target * unit_GHz, target_max * unit_GHz);
      }
    }
  }
  if (hybrid_seed)
  {
    const double target_max = iodata.solver.eigenmode.target_upper;
    auto interp = std::make_unique<NewtonInterpolationOperator>(funcA2, A2->Width());
    interp->Interpolate(1i * target, 1i * target_max);
    // Frozen-ABC seed for the NLEPS HYBRID polynomial seed pencil. The 2nd-order farfield
    // ABC contributes pole terms g(λ)·M to A2(λ), with g(λ) = -0.5/λ next to constant
    // materials (M includes the negative sign of β) and 0.5i/(λ r_k(-iλ)) next to
    // frequency-dependent material k, that a polynomial (K' + λC' + λ²M') seed cannot fit
    // accurately. AddFrozenPole removes each term's interpolated contribution and re-adds
    // it frozen at the target into the K-block: next to a frequency-dependent material, the
    // material factor 1/r_k is frozen with the pole, as the seed has no separate
    // representation of it. The freeze is deliberately unconditional (no DetermineFrozen):
    // the fit's pointwise window error is small (the pole at λ = 0 sits outside the
    // window), but its curvature term injects a fictitious rank-deficient λ² contribution
    // into the seed's M-block that displaces spurious roots into the unphysical half-plane
    // — the failure mode is structural, not a pointwise approximation error, so the
    // pointwise metric must not be allowed to choose the fit here.
    for (auto &[M_ff, g] : space_op.GetFarfieldPoleTerms(Operator::DIAG_ZERO))
    {
      interp->AddFrozenPole(std::move(M_ff), g, 1i * target);
    }
    // Fit-or-freeze seed for rational impedance boundaries. Each boundary contributes
    // g(λ)·M_b to A2(λ) with g(λ) = λ·D(λ)/N(λ) = P(λ) + R(λ)/N(λ). The polynomial part
    // P is exactly representable in the seed pencil when deg(P) <= 2, so only the strictly
    // proper pole part R/N is a candidate for freezing. We freeze it at the target when
    // that approximates g over the interpolation window better than the polynomial
    // interpolant does; otherwise keep the fit.
    const auto &surf_rz_op = space_op.GetRationalImpedanceOp();
    for (int idx = 0; idx < surf_rz_op.GetNumBoundaries(); idx++)
    {
      const bool proper_split = (surf_rz_op.GetRobinQuotientDegree(idx) <= 2);
      auto f_full = [&surf_rz_op, idx](std::complex<double> lambda)
      { return surf_rz_op.EvalRobinCoefficient(idx, lambda); };
      auto f_frozen = [&surf_rz_op, idx, proper_split](std::complex<double> lambda)
      {
        return proper_split ? surf_rz_op.EvalRobinRemainder(idx, lambda)
                            : surf_rz_op.EvalRobinCoefficient(idx, lambda);
      };
      double fit_err, freeze_err;
      if (interp->DetermineFrozen(f_full, f_frozen, 1i * target, fit_err, freeze_err))
      {
        Mpi::Print(" Freezing rational impedance boundary (attribute {}) pole part in the "
                   "NLEPS seed (fit error {:.2e} > freeze error {:.2e})\n",
                   fmt::join(surf_rz_op.GetAttrList(idx), ", "), fit_err, freeze_err);
        auto M_b = space_op.GetRationalImpedanceBoundaryMassMatrix<ComplexOperator>(
            idx, Operator::DIAG_ZERO);
        interp->AddFrozenPole(std::move(M_b), f_frozen, 1i * target);
      }
    }
    A2_0 = interp->GetInterpolationOperator(0);
    A2_1 = interp->GetInterpolationOperator(1);
    A2_2 = interp->GetInterpolationOperator(2);
    interp_op = std::move(interp);  // retain: A2_0/A2_1/A2_2 reference its operator DAG
    Kp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {K.get(), A2_0.get()});
    Cp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {C.get(), A2_1.get()});
    Mp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i}, {M.get(), A2_2.get()});
  }

  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);
  ComplexVector B(Curl.Height());
  B.UseDevice(true);

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  const EigenSolverBackend type = iodata.solver.eigenmode.type;
#if !defined(PALACE_WITH_ARPACK) && !defined(PALACE_WITH_SLEPC)
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C || has_A2)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
#endif
  }
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      slepc = std::make_unique<slepc::SlepcNEPSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
    }
    else
    {
      if (C || has_A2)
      {
        if (!iodata.solver.eigenmode.pep_linear)
        {
          slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                          iodata.problem.verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR);
        }
        else
        {
          slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                                iodata.problem.verbose);
          slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
        }
      }
      else
      {
        slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    }
    slepc->SetOrthogonalization(iodata.solver.linear.gs_orthog == Orthogonalization::MGS,
                                iodata.solver.linear.gs_orthog == Orthogonalization::CGS2);
    eigen = std::move(slepc);
#endif
  }
  if (exact_poles)
  {
    std::size_t num_poles = 0;
    for (const auto &term : a2_terms)
    {
      eigen->AddRationalTerm(*term.B, term.poly, term.poles, term.residues);
      num_poles += term.poles.size();
    }
    Mpi::Print(" Frequency-dependent terms linearized {}({:d} term{}, {:d} pole{}), "
               "without nonlinear refinement\n",
               fitted ? "with rational fits " : "exactly ", a2_terms.size(),
               (a2_terms.size() > 1) ? "s" : "", num_poles, (num_poles > 1) ? "s" : "");
    if (fitted)
    {
      // Report the errors of the complete operator, not of its rational model.
      eigen->SetExtraSystemMatrix(funcA2_full);
    }
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (nonlinear_type == NonlinearEigenSolver::SLP)
  {
    if (C)
    {
      eigen->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      eigen->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    eigen->SetExtraSystemMatrix(funcA2_full);
    eigen->SetExtraSystemPoles(a2_poles);
    eigen->SetPreconditionerUpdate(funcP);
  }
  else
  {
    if (hybrid_seed)
    {
      eigen->SetOperators(*Kp, *Cp, *Mp, scale);
    }
    else if (C)
    {
      eigen->SetOperators(*K, *C, *M, scale);
    }
    else
    {
      eigen->SetOperators(*K, *M, scale);
    }
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  const double tol =
      hybrid_seed ? iodata.solver.eigenmode.linear_tol : iodata.solver.eigenmode.tol;
  eigen->SetTol(tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  std::unique_ptr<Operator> KM;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  if (iodata.solver.linear.divfree_max_it > 0 &&
      !space_op.GetMaterialOp().HasWaveVector() &&
      !space_op.GetMaterialOp().HasLondonDepth() &&
      !space_op.GetMaterialOp().HasFrequencyDependentPermittivity())
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }
  else if (iodata.solver.linear.divfree_max_it > 0 &&
           space_op.GetMaterialOp().HasFrequencyDependentPermittivity())
  {
    Mpi::Print(" Skipping divergence-free projection for frequency-dependent permittivity: "
               "the fixed ε∞ projector does not enforce ∇·(ε(ω)E) = 0 at material "
               "interfaces.\n");
  }

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      space_op.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  {
    const double f_target =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) /
        (2 * std::numbers::pi);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (C || has_A2 || nonlinear_type == NonlinearEigenSolver::SLP)
  {
    // Search for eigenvalues closest to λ = iσ.
    eigen->SetShiftInvert(1i * target);
    if (nonlinear_type == NonlinearEigenSolver::SLP)
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    eigen->SetShiftInvert(target * target);
    if (type == EigenSolverBackend::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  std::unique_ptr<ComplexOperator> A;
  std::unique_ptr<ComplexOperator> A2_shift;
  if (hybrid_seed)
  {
    // Invert the same W-aware polynomial used by the seed eigensolver.
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, Kp.get(),
                                 Cp.get(), Mp.get());
  }
  else
  {
    // SLP applies the exact nonlinear operator at the target; with exactly linearized
    // permittivity poles, A2 at the target completes T(σ); linear problems have no A2.
    A2_shift = (nonlinear_type == NonlinearEigenSolver::SLP) ? funcA2_full(1i * target)
                                                             : std::move(A2);
    if (fitted)
    {
      // The linearization requires its own model at the shift, which differs from A2 by
      // the fit error.
      std::vector<std::complex<double>> coeff;
      std::vector<const ComplexParOperator *> ops;
      for (const auto &term : a2_terms)
      {
        coeff.push_back(term.Eval(1i * target));
        ops.push_back(dynamic_cast<const ComplexParOperator *>(term.B.get()));
        MFEM_VERIFY(ops.back(), "Unexpected operator type of a rational term!");
      }
      A2_shift = BuildParSumOperator(coeff, ops);
    }
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, K.get(),
                                 C.get(), M.get(), A2_shift.get());
  }
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(
      1.0 + 0.0i, 1i * target, -target * target + 0.0i, target + 0.0i);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);

  // Initialize structures for storing and reducing the results of error estimation.
  const FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);
  ErrorIndicator indicator;

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  if (hybrid_seed)
  {
    Mpi::Print("\n Refining eigenvalues with Quasi-Newton solver\n");
    auto qn = std::make_unique<QuasiNewtonSolver>(space_op.GetComm(), std::move(eigen),
                                                  num_conv, iodata.problem.verbose,
                                                  iodata.solver.eigenmode.refine_nonlinear);
    qn->SetTol(iodata.solver.eigenmode.tol);
    qn->SetMaxIter(iodata.solver.eigenmode.nonlinear_max_it);
    if (C)
    {
      qn->SetOperators(*K, *C, *M, EigenvalueSolver::ScaleType::NONE);
    }
    else
    {
      qn->SetOperators(*K, *M, EigenvalueSolver::ScaleType::NONE);
    }
    qn->SetExtraSystemMatrix(funcA2_full);
    qn->SetPreconditionerUpdate(funcP);
    qn->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
    qn->SetPreconditionerLag(iodata.solver.eigenmode.preconditioner_lag,
                             iodata.solver.eigenmode.preconditioner_lag_tol);
    qn->SetMaxRestart(iodata.solver.eigenmode.max_restart);
    qn->SetLinearSolver(*ksp);
    qn->SetShiftInvert(1i * target);
    eigen = std::move(qn);

    // Suppress wave port output during nonlinear eigensolver iterations.
    space_op.GetWavePortOp().SetSuppressOutput(true);
    num_conv = eigen->Solve();
    space_op.GetWavePortOp().SetSuppressOutput(false);
  }

  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  const int num_requested = iodata.solver.eigenmode.n;
  const int num_report = num_conv;
  MFEM_VERIFY(num_report >= num_requested, "Eigenmode solve only found "
                                               << num_report << " modes when "
                                               << num_requested << " were requested!");
  if (!KM)
  {
    // Normalize every converged eigenvector with respect to the mass matrix (unit electric
    // field energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_report);
  }
  std::vector<std::complex<double>> eigenvalues(num_report);
  std::vector<double> errors_bkwd(num_report), errors_abs(num_report);
  const bool stage_eigenvectors = post_op.WillWriteFields();
  // Field output can construct large libCEED visualization operators. Only in that case,
  // stage the converged eigenvectors in host memory so the eigensolver and matrices can
  // release their device allocations before postprocessing.
  std::vector<std::vector<std::complex<double>>> eigenvectors;
  if (stage_eigenvectors)
  {
    eigenvectors.assign(num_report, std::vector<std::complex<double>>(Curl.Width()));
  }
  ComplexVector E(Curl.Width());
  E.UseDevice(true);
  for (int i = 0; i < num_report; i++)
  {
    eigenvalues[i] = eigen->GetEigenvalue(i);
    errors_bkwd[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    errors_abs[i] = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    if (stage_eigenvectors)
    {
      eigen->GetEigenvector(i, E);
      E.Get(eigenvectors[i].data(), E.Size(), false);
    }
  }

  // The target window is checked with the configuration; a mode may converge outside it.
  for (std::size_t m = 0; m < iodata.domains.materials.size(); m++)
  {
    const auto &material = iodata.domains.materials[m];
    if (permittivity::IsPassiveByConstruction(material.permittivity_terms))
    {
      continue;
    }
    const auto defects =
        permittivity::FindDefects(material.epsilon_r.s[0], material.permittivity_terms);
    for (int i = 0; i < num_report; i++)
    {
      const double omega = eigenvalues[i].imag();
      if (permittivity::Contains(defects.gain, omega))
      {
        Mpi::Warning("Eigenmode {:d} at {:.6g} GHz lies where the frequency-dependent "
                     "Permittivity of material {:d} has gain, a negative dissipated "
                     "power!\nThe mode solves the model as given; a fit may only be meant "
                     "for some frequencies.\n",
                     i + 1,
                     iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(omega) /
                         (2.0 * std::numbers::pi),
                     m + 1);
      }
    }
  }

  for (int i = 0; i < num_report; i++)
  {
    const std::complex<double> omega = eigenvalues[i] / 1i;
    for (std::size_t t = 0; t < a2_terms.size(); t++)
    {
      const auto &term = a2_terms[t];
      if (!term.f || term.fit.IsValid(window, omega) ||
          std::ranges::any_of(a2_terms.begin(), a2_terms.begin() + t, [&](const auto &other)
                              { return other.f && !other.fit.IsValid(window, omega); }))
      {
        continue;
      }
      const auto f = ToGHz(omega);
      Mpi::Warning(
          "Eigenmode {:d} at {:.6e}{:+.6e}i GHz lies where the rational fit of the "
          "{} is not valid ({}): it solves the fitted problem, and its error is "
          "that of the complete operator.\n",
          i + 1, f.real(), f.imag(), term.name,
          window.Contains(omega) ? "next to its branch cut"
                                 : "outside the window of the fit");
    }
  }

  if (stage_eigenvectors)
  {
    eigen.reset();
    KM.reset();
    ksp.reset();
    P.reset();
    A.reset();
    divfree.reset();
    Kp.reset();
    Cp.reset();
    Mp.reset();
    A2_0.reset();
    A2_1.reset();
    A2_2.reset();
    A2.reset();
    interp_op.reset();
    K.reset();
    C.reset();
    M.reset();
  }
  Mpi::Print("\n");

  // A mode whose magnetic energy is negligible against its electric energy is a
  // longitudinal field, such as the gradient fields which solve the equations at a zero of
  // a frequency-dependent permittivity (see the guide). Their curl, and so their magnetic
  // energy, is the error of the solve: 3e-16 to 1.2e-15 of the electric energy for
  // cavity_lorentz.json at Tol 1e-8 from a Target of 2.905 GHz. The modes of that example
  // which accumulate at its resonance have 1.4e-2. The modes are kept and listed.
  constexpr double longitudinal_tol = 1.0e-8;
  std::vector<int> longitudinal;
  for (int i = 0; i < num_report; i++)
  {
    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigenvalues[i];
    double error_bkwd = errors_bkwd[i];
    double error_abs = errors_abs[i];
    if (!has_C && !has_A2)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    if (stage_eigenvectors)
    {
      E.Set(eigenvectors[i].data(), E.Size(), false);
    }
    else
    {
      eigen->GetEigenvector(i, E);
    }
    linalg::NormalizePhase(space_op.GetComm(), E);

    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E.
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_report);
    if (const auto [energy_elec, energy_mag] = post_op.GetDomainFieldEnergies();
        energy_elec && energy_mag < longitudinal_tol * *energy_elec)
    {
      longitudinal.push_back(i + 1);
    }
    if (i < num_requested)
    {
      // Frequency-dependent materials use Re{ω} of the eigenvalue, as in postprocessing.
      estimator.AddErrorIndicator(E, B, total_domain_energy, omega.real(), indicator);
    }
  }
  if (!longitudinal.empty())
  {
    Mpi::Warning(
        "Eigenmode{} {} {} a magnetic energy below {:.0e} of the electric energy!\n"
        "These are longitudinal fields, which solve the equations without a "
        "magnetic field, such as the gradient fields at a zero of a "
        "frequency-dependent permittivity.\n",
        (longitudinal.size() > 1) ? "s" : "", fmt::join(longitudinal, ", "),
        (longitudinal.size() > 1) ? "have" : "has", longitudinal_tol);
  }
  post_op.MeasureFinalize(indicator);
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
