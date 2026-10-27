// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <numbers>
#include <vector>
#include <mfem.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/fieldtransfer.hpp"
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
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Add the eigenpair (λ, x) to the combinations v = Σᵢ xᵢ / ‖xᵢ‖ and, if v_lambda is given,
// v_lambda = Σᵢ λᵢ xᵢ / ‖xᵢ‖, from which the next adaptation iteration starts. Returns
// whether x was added.
bool AddWarmStartEigenvector(MPI_Comm comm, const ComplexVector &x,
                             std::complex<double> lambda, ComplexVector &v,
                             ComplexVector *v_lambda)
{
  const double norm = linalg::Norml2(comm, x);
  if (!(norm > 0.0))
  {
    return false;
  }
  v.AXPY(1.0 / norm, x);
  if (v_lambda)
  {
    v_lambda->AXPY(lambda / norm, x);
  }
  return true;
}

// Starting vector v of a warm-started solve, and the second block v_lambda if given, from
// the combinations of the previous adaptation iteration (fields 0, 1 and 2, 3 of
// prev_vectors, as real and imaginary parts), plus eta ‖v‖ / ‖r‖ times the default
// starting vector r for modes not found before.
void GetWarmStartVector(const VectorFieldTransfer &prev_vectors, SpaceOperator &space_op,
                        const ComplexVector &r, double eta, ComplexVector &v,
                        ComplexVector *v_lambda)
{
  auto &nd_fespace = space_op.GetNDSpace().Get();
  const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists().back();
  const int n = nd_fespace.GetTrueVSize();
  for (int b = 0; b < (v_lambda ? 2 : 1); b++)
  {
    ComplexVector &u = (b == 0) ? v : *v_lambda;
    u.SetSize(n);
    u.UseDevice(true);
    prev_vectors.GetField(2 * b, nd_fespace, u.Real());
    prev_vectors.GetField(2 * b + 1, nd_fespace, u.Imag());
    linalg::SetSubVector(u, dbc_tdof_list, 0.0);
  }
  const double norm_r = linalg::Norml2(space_op.GetComm(), r);
  if (eta > 0.0 && norm_r > 0.0)
  {
    v.AXPY(eta * linalg::Norml2(space_op.GetComm(), v) / norm_r, r);
  }
}

}  // namespace

EigenSolver::EigenSolver(const IoData &iodata, bool root, int size, int num_thread,
                         const char *git_tag)
  : BaseSolver(iodata, root, size, num_thread, git_tag)
{
}

EigenSolver::~EigenSolver() = default;

void EigenSolver::UpdateAfterMeshChange(Mesh &mesh) const
{
  if (!prev_eigenvectors)
  {
    return;
  }
  if (prev_eigenvectors->GetParMesh() != &mesh.Get())
  {
    // The mesh object was replaced.
    prev_eigenvectors.reset();
    return;
  }
  prev_eigenvectors->Update();
}

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
  std::array<std::unique_ptr<ComplexOperator>, 3> A2_pml;
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
    // The SLP nonlinear eigensolver requires a frequency-dependent term A2(λ). Without one
    // the NEP function/Jacobian shells dereference a null A2; fall back to the standard
    // linear/quadratic eigensolver instead.
    Mpi::Warning("SLP nonlinear eigensolver requires a frequency-dependent system term, "
                 "none found; using a linear eigensolver!\n");
    nonlinear_type = NonlinearEigenSolver::HYBRID;
  }

  // Frequency-dependent permittivity pole terms: their part of A2(λ) is Σₘ gₘ(λ) Bₘ with
  // gₘ(λ) = Σₖ rₖ λ² / (λ - pₖ), rational with known poles, which the linearized
  // eigensolvers represent exactly, with one auxiliary block per pole
  // (EigenvalueSolver::AddPoleTerms). When these are the only frequency-dependent terms,
  // the problem is then a linear eigenvalue problem: there is no polynomial seed and no
  // nonlinear refinement. RefineNonlinear has no effect here: the solution is exact and
  // there is nothing to refine. Djordjevic-Sarkar terms and the frequency-dependent
  // boundary terms are not of this form, and keep the Hybrid seed and the Quasi-Newton
  // refinement: the seed then takes the pole terms exactly as well, and interpolates only
  // the rest of A2, which does not change on the scale of a material resonance (unless a
  // wave port or an absorbing boundary is next to the material). Interpolated over the
  // search window, a pole term has spurious roots and misses the modes which accumulate at
  // the resonance. With PEPLinear = false, SLEPc's polynomial eigensolver linearizes the
  // problem itself, without the auxiliary blocks.
  std::vector<std::unique_ptr<ComplexOperator>> pole_ops;
  std::vector<std::vector<std::complex<double>>> pole_p, pole_c;
  // The pole terms at λ, Σⱼ cₜⱼ λ² / (λ - pₜⱼ) for the operator of each material.
  auto PoleTermCoefficient = [&pole_p, &pole_c](std::size_t t, std::complex<double> lambda)
  {
    std::complex<double> s = 0.0;
    for (std::size_t j = 0; j < pole_p[t].size(); j++)
    {
      s += pole_c[t][j] * lambda * lambda / (lambda - pole_p[t][j]);
    }
    return s;
  };
  bool poles_only = false;
  if (has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID &&
      (iodata.solver.eigenmode.type == EigenSolverBackend::ARPACK ||
       iodata.solver.eigenmode.pep_linear))
  {
    // The frequency-dependent boundary mass terms of SurfaceTerms are those of the finite
    // conductivity and rational impedance boundaries and of the absorbing boundaries next
    // to frequency-dependent materials.
    const auto &mat_op = space_op.GetMaterialOp();
    const auto &surface_terms = space_op.GetSurfaceTerms();
    const auto &farfield_op = surface_terms.GetFarfieldOp();
    poles_only = space_op.GetWavePortOp().Size() == 0 &&
                 surface_terms.GetNumExtraSystemBdrMassTerms() == 0 &&
                 (farfield_op.GetAttrList().Size() == 0 || farfield_op.GetOrder() < 2) &&
                 !mat_op.HasFrequencyDependentPML();
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (mat_op.HasFrequencyDependentPermittivitySupport(m) &&
          mat_op.HasDjordjevicSarkarPermittivity(m))
      {
        poles_only = false;
      }
    }
    for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
    {
      if (!mat_op.HasFrequencyDependentPermittivitySupport(m) ||
          mat_op.GetPermittivityPoleTerms(m).empty())
      {
        continue;  // No contribution (SpaceOperator::GetExtraSystemMatrix), or no pole
      }
      auto B = space_op.GetFrequencyDependentPermittivityMassMatrix(m, Operator::DIAG_ZERO);
      if (!B)
      {
        continue;  // No nonzero pole (exact zero poles are in the conductivity of C)
      }
      // A2 holds gₘ(λ)·Bₘ, and the operator is on the imaginary slot: it applies i·Bₘ.
      auto &p = pole_p.emplace_back();
      auto &c = pole_c.emplace_back();
      for (const auto &term : mat_op.GetPermittivityPoleTerms(m))
      {
        p.push_back(term.pole);
        c.push_back(-1i * term.residue);
      }
      pole_ops.push_back(std::move(B));
    }
    if (poles_only && !pole_ops.empty())
    {
      // The terms must be the complete A2(λ): any other contribution is an error here, not
      // an approximation. Checked at one complex λ, off the interpolation points.
      const std::complex<double> l(-0.01 * target, 1.05 * target);
      ComplexVector v, r, r_poles;
      space_op.GetRandomInitialVector(v, 2);
      r.SetSize(v.Size());
      r.UseDevice(true);
      r_poles.SetSize(v.Size());
      r_poles.UseDevice(true);
      funcA2(l)->Mult(v, r);
      r_poles = 0.0;
      for (std::size_t t = 0; t < pole_ops.size(); t++)
      {
        pole_ops[t]->AddMult(v, r_poles, PoleTermCoefficient(t, l));
      }
      const double norm = linalg::Norml2(space_op.GetComm(), r);
      r_poles.AXPY(-1.0, r);
      const double err = linalg::Norml2(space_op.GetComm(), r_poles);
      MFEM_VERIFY(err <= 1.0e-10 * norm,
                  "The permittivity pole terms do not reproduce the frequency-dependent "
                  "operator (relative difference "
                      << err / norm << ")!");
    }
  }
  const bool exact_poles = poles_only && !pole_ops.empty();
  const bool hybrid_seed =
      has_A2 && nonlinear_type == NonlinearEigenSolver::HYBRID && !exact_poles;
  const bool seed_poles = hybrid_seed && !pole_ops.empty();
  if (hybrid_seed || nonlinear_type == NonlinearEigenSolver::SLP)
  {
    // The nonlinear eigensolver can miss modes or stop near a material resonance, where
    // A2(λ) changes on the scale of its half-width (see the guide). The search window is
    // not changed: that would change the results.
    const auto &mat_op = space_op.GetMaterialOp();
    const double target_max = iodata.solver.eigenmode.target_upper;
    const double unit_GHz = iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
                            (2.0 * std::numbers::pi);
    // A wave port keeps the problem on this route whatever the materials. Its mode is
    // interpolated in the seed, and changes on the scale of a resonance where the port is
    // in a frequency-dependent material: warn with any such material on the mesh, also
    // without a resonance near the window (a Drude, Debye, or Djordjevic-Sarkar material).
    // The terms of an absorbing boundary (of either order) next to a frequency-dependent
    // material are interpolated in the seed in the same way, and also skip modes: one
    // warning names the wave ports, the absorbing boundary, or both. An absorbing boundary
    // next to constant materials only does not change on that scale.
    const bool wave_port =
        space_op.GetWavePortOp().Size() > 0 && mat_op.HasFrequencyDependentPermittivity();
    const bool absorbing =
        space_op.GetSurfaceTerms().GetFarfieldOp().NumDispersiveMaterials() > 0;
    if (wave_port || absorbing)
    {
      Mpi::Warning(
          "The eigenmode problem has {} a frequency-dependent Permittivity, and is "
          "solved with the nonlinear eigensolver!\nIt can miss modes nearest the "
          "target or stop with fewer modes than requested, most of all near a "
          "resonance of a material at {}.\n",
          !absorbing  ? "a wave port and"
          : wave_port ? "a wave port and an absorbing boundary next to"
                      : "an absorbing boundary next to",
          !absorbing  ? "a wave port"
          : wave_port ? "a wave port or the absorbing boundary"
                      : "the absorbing boundary");
    }
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
    // Frozen-stretch seed for frequency-dependent PML regions. Their A2(λ) terms are
    // rational in λ through the stretch factors s(λ) = κ + σ / (α + λ), with a continuum of
    // poles on the negative real axis (σ is graded across the layer). A quadratic fit over
    // the interpolation window, up to the upper target, cannot represent them and seeds
    // spurious, strongly damped modes. Instead, the PML stretch is frozen at the target:
    // the PML terms of the seed pencil are those of a static PML with reference frequency
    // equal to the target, which is exact at the target and a consistent PML (without
    // spurious modes) elsewhere.
    const bool freeze_pml = space_op.GetMaterialOp().HasFrequencyDependentPML();
    // With the pole terms or the PML terms in the seed exactly, only the rest of A2 is
    // interpolated.
    std::function<std::unique_ptr<ComplexOperator>(std::complex<double>)> funcA2_seed =
        funcA2;
    if (seed_poles || freeze_pml)
    {
      funcA2_seed = [&space_op, seed_poles, freeze_pml](std::complex<double> lambda)
      {
        const std::complex<double> omega = lambda / std::complex<double>(0.0, 1.0);
        return space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO,
                                             /*include_permittivity_poles=*/!seed_poles,
                                             /*include_pml=*/!freeze_pml);
      };
    }
    if (!freeze_pml || space_op.HasExtraSystemBdrTerms() ||
        space_op.GetMaterialOp().HasFrequencyDependentPermittivityA2())
    {
      const double target_max = iodata.solver.eigenmode.target_upper;
      auto interp = std::make_unique<NewtonInterpolationOperator>(funcA2_seed, A2->Width());
      interp->Interpolate(1i * target, 1i * target_max);
      // Frozen-ABC seed for the NLEPS HYBRID polynomial seed pencil. The 2nd-order farfield
      // ABC contributes pole terms g(λ)·M to A2(λ), with g(λ) = -0.5/λ next to constant
      // materials (M includes the negative sign of β) and 0.5i/(λ r_k(-iλ)) next to
      // frequency-dependent material k, that a polynomial (K' + λC' + λ²M') seed cannot fit
      // accurately. AddFrozenPole removes each term's interpolated contribution and re-adds
      // it frozen at the target into the K-block: next to a frequency-dependent material,
      // the material factor 1/r_k is frozen with the pole, as the seed has no separate
      // representation of it. The freeze is deliberately unconditional (no
      // DetermineFrozen): the fit's pointwise window error is small (the pole at λ = 0 sits
      // outside the window), but its curvature term injects a fictitious rank-deficient λ²
      // contribution into the seed's M-block that displaces spurious roots into the
      // unphysical half-plane — the failure mode is structural, not a pointwise
      // approximation error, so the pointwise metric must not be allowed to choose the fit
      // here.
      for (auto &[M_ff, g] : space_op.GetFarfieldPoleTerms(Operator::DIAG_ZERO))
      {
        interp->AddFrozenPole(std::move(M_ff), g, 1i * target);
      }
      // Fit-or-freeze seed for the frequency-dependent boundary mass terms with an exact
      // split of their scale (rational impedance boundaries). Each such term contributes
      // c(λ)·M_i to A2(λ) with c(λ) = P(λ) + R(λ), for a rational impedance boundary the
      // Robin coefficient λ·D(λ)/N(λ) with the pole part R(λ)/N(λ). The polynomial part P
      // is exactly representable in the seed pencil when deg(P) <= 2, so only the strictly
      // proper pole part R is a candidate for freezing (otherwise the whole coefficient).
      // We freeze it at the target when that approximates c over the interpolation window
      // better than the polynomial interpolant does; otherwise keep the fit. Terms without
      // a split (finite conductivity) are left to the interpolant.
      const auto &surface_terms = space_op.GetSurfaceTerms();
      for (int i = 0; i < surface_terms.GetNumExtraSystemBdrMassTerms(); i++)
      {
        const auto split = surface_terms.GetExtraSystemBdrMassTermSplit(i);
        if (!split)
        {
          continue;
        }
        const bool proper_split = (split->degree <= 2);
        const auto &f_full = split->coefficient;
        const auto &f_frozen = proper_split ? split->remainder : split->coefficient;
        double fit_err, freeze_err;
        if (interp->DetermineFrozen(f_full, f_frozen, 1i * target, fit_err, freeze_err))
        {
          Mpi::Print(" Freezing {} pole part in the NLEPS seed (fit error {:.2e} > freeze "
                     "error {:.2e})\n",
                     surface_terms.GetExtraSystemBdrMassTermName(i), fit_err, freeze_err);
          auto M_i = space_op.GetExtraSystemBdrMassTermMatrix<ComplexOperator>(
              i, Operator::DIAG_ZERO);
          interp->AddFrozenPole(std::move(M_i), f_frozen, 1i * target);
        }
      }
      A2_0 = interp->GetInterpolationOperator(0);
      A2_1 = interp->GetInterpolationOperator(1);
      A2_2 = interp->GetInterpolationOperator(2);
      interp_op = std::move(interp);  // retain: A2_0/A2_1/A2_2 reference its operator DAG
    }
    if (freeze_pml)
    {
      Mpi::Print(" Freezing the frequency-dependent PML stretch at the target in the NLEPS "
                 "seed\n");
      A2_pml = space_op.GetFrozenPMLMatrices(target, Operator::DIAG_ZERO);
    }
    Kp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i, 1.0 + 0i},
                             {K.get(), A2_0.get(), A2_pml[0].get()});
    if (C || A2_1 || A2_pml[1])
    {
      Cp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i, 1.0 + 0i},
                               {C.get(), A2_1.get(), A2_pml[1].get()});
    }
    Mp = BuildParSumOperator({1.0 + 0i, 1.0 + 0i, 1.0 + 0i},
                             {M.get(), A2_2.get(), A2_pml[2].get()});
  }

  // The (seed) eigenvalue problem is quadratic in λ if it has a damping term, otherwise it
  // is a linear generalized eigenvalue problem for μ = -λ² = ω².
  const bool quadratic = hybrid_seed ? (Cp != nullptr) : (C || has_A2);

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
    if (quadratic)
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
      if (quadratic)
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
  if (exact_poles || seed_poles)
  {
    std::size_t num_poles = 0;
    for (std::size_t t = 0; t < pole_ops.size(); t++)
    {
      eigen->AddPoleTerms(*pole_ops[t], pole_p[t], pole_c[t]);
      num_poles += pole_p[t].size();
    }
    if (exact_poles)
    {
      Mpi::Print(" Frequency-dependent permittivity linearized exactly ({:d} pole{}), "
                 "without nonlinear refinement\n",
                 num_poles, (num_poles > 1) ? "s" : "");
    }
    else
    {
      Mpi::Print(" Frequency-dependent permittivity poles taken exactly in the initial "
                 "guesses of the nonlinear eigensolver ({:d} pole{})\n",
                 num_poles, (num_poles > 1) ? "s" : "");
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
      if (Cp)
      {
        eigen->SetOperators(*Kp, *Cp, *Mp, scale);
      }
      else
      {
        eigen->SetOperators(*Kp, *Mp, scale);
      }
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
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get(), target);
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
      !space_op.GetMaterialOp().HasFrequencyDependentPermittivity() &&
      !space_op.GetMaterialOp().HasPML())
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

  // With adaptive mesh refinement, the SLEPc Krylov-Schur solver can start from the
  // eigenvectors of the previous solve (warm start). Its more accurate first Krylov cycle
  // can cost more than a warm start saves on solves of two cycles: the eigenvectors are
  // kept only if the last solve without a warm start took three cycles or more, and only if
  // the adaptation goes on.
  constexpr int min_cold_its = 3;
  const auto &eigenmode = iodata.solver.eigenmode;
  const auto &refinement = iodata.model.refinement;
  const int adapt_it = num_solves++;
  const bool blocks = quadratic;  // The linearization also takes the block λ x
  const bool can_warm_start = eigenmode.init_v0 && type == EigenSolverBackend::SLEPC &&
                              nonlinear_type != NonlinearEigenSolver::SLP &&
                              (!blocks || eigenmode.pep_linear) &&
                              space_op.GetNDSpace().GetParMesh().Nonconforming();

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  int num_warm = 0;
  if (eigenmode.init_v0)
  {
    ComplexVector v0;
    if (eigenmode.init_v0_const)
    {
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      space_op.GetRandomInitialVector(v0);
    }
    if (prev_eigenvectors &&
        prev_eigenvectors->IsCurrent(space_op.GetNDSpace().GetParMesh()))
    {
      num_warm = prev_num_eigenvectors;
      Mpi::Print(" Using starting vector from {:d} eigenvectors of the previous adaptation "
                 "iteration\n",
                 num_warm);
      constexpr double eta = 1.0e-3;
      ComplexVector v, v_lambda;
      GetWarmStartVector(*prev_eigenvectors, space_op, v0, eta, v,
                         blocks ? &v_lambda : nullptr);
      if (divfree)
      {
        divfree->Mult(v);
        if (blocks)
        {
          divfree->Mult(v_lambda);
        }
      }
      if (blocks)
      {
        eigen->SetInitialSpaceBlocks(v, v_lambda);  // Copies the vectors
      }
      else
      {
        eigen->SetInitialSpace(v);  // Copies the vector
      }

      // A warm-started solve converges after few linear solves, and their error, kept in
      // the Krylov vectors through the restarts, is not seen by the error estimate of the
      // eigenvalue solver: tighten the first cycle by 100, but not below 1e-3 times the
      // tolerance of the eigenvalue solve (of the initial guess for the hybrid nonlinear
      // solver) or 1e-12.
      eigen->SetFirstCycleLinearSolverTol(
          std::max({1.0e-2 * iodata.solver.linear.tol, 1.0e-3 * tol, 1.0e-12}));
    }
    else
    {
      if (adapt_it > 0 && can_warm_start && last_cold_its < min_cold_its)
      {
        Mpi::Print(" No warm start: the last solve without one took {:d} Krylov cycle{}\n",
                   last_cold_its, (last_cold_its > 1) ? "s" : "");
      }
      Mpi::Print(" Using {} starting vector\n",
                 eigenmode.init_v0_const ? "constant" : "random");
      if (divfree)
      {
        divfree->Mult(v0);
      }
      eigen->SetInitialSpace(v0);  // Copies the vector
    }

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }
  prev_eigenvectors.reset();

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  {
    const double f_target =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target) /
        (2 * std::numbers::pi);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  if (quadratic || nonlinear_type == NonlinearEigenSolver::SLP)
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
    // Invert the same W-aware polynomial used by the seed eigensolver, with its pole terms
    // at the target.
    if (seed_poles)
    {
      std::vector<std::complex<double>> coeff;
      std::vector<const ComplexOperator *> ops;
      for (std::size_t t = 0; t < pole_ops.size(); t++)
      {
        coeff.push_back(PoleTermCoefficient(t, 1i * target));
        ops.push_back(pole_ops[t].get());
      }
      A2_shift = BuildParSumOperator(coeff, ops);
    }
    A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * target, -target * target + 0.0i, Kp.get(),
                                 Cp.get(), Mp.get(), A2_shift.get());
  }
  else
  {
    // SLP applies the exact nonlinear operator at the target; with exactly linearized
    // permittivity poles, A2 at the target completes T(σ); linear problems have no A2.
    A2_shift = (nonlinear_type == NonlinearEigenSolver::SLP) ? funcA2_full(1i * target)
                                                             : std::move(A2);
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
  if (num_warm == 0)
  {
    last_cold_its = eigen->GetNumIterations();
  }
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
    qn->SetLinearEigenvalueSquared(!quadratic);
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
  MFEM_VERIFY(num_report > 0, "Eigenmode solve only found 0 converged modes when "
                                  << num_requested << " were requested!");
  if (num_report < num_requested)
  {
    // Postprocess the converged modes rather than abort, for example when clusters of modes
    // (such as the modes of a PML) keep the eigenvalue solver from converging all of the
    // requested modes within the maximum number of iterations.
    Mpi::Warning(space_op.GetComm(),
                 "Eigenmode solve only found {:d} converged mode{} when {:d} were "
                 "requested, postprocessing the converged modes (increase "
                 "config[\"Solver\"][\"Eigenmode\"][\"MaxIts\"] to converge more modes)!\n",
                 num_report, (num_report > 1) ? "s" : "", num_requested);
  }
  if (!KM)
  {
    // Normalize every converged eigenvector with respect to the mass matrix (unit electric
    // field energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get(), target);
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_report);
  }
  std::vector<std::complex<double>> eigenvalues(num_report);
  std::vector<double> errors_bkwd(num_report), errors_abs(num_report);

  // Combine the converged eigenvectors to start the next adaptation iteration: only the
  // combinations follow the refinement and rebalancing of the mesh.
  const bool combine_eigenvectors =
      can_warm_start && last_cold_its >= min_cold_its && adapt_it < refinement.max_it;
  ComplexVector v_comb, v_comb_lambda;
  int num_combined = 0;
  if (combine_eigenvectors)
  {
    for (auto *u : {&v_comb, &v_comb_lambda})
    {
      if (u == &v_comb || blocks)
      {
        u->SetSize(Curl.Width());
        u->UseDevice(true);
        *u = 0.0;
      }
    }
  }
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
    for (auto &A2_pml_k : A2_pml)
    {
      A2_pml_k.reset();
    }
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
  // which accumulate at its resonance have 1.4e-2. The modes are kept and listed. The
  // electric energy is the field energy with ε∞ (PostOperator::GetDomainFieldEnergies):
  // the stored energy is not specified for a PoleResidue term, and the material energy,
  // with Re ε(ω), vanishes at the zero of ε.
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
    if (combine_eigenvectors)
    {
      num_combined += AddWarmStartEigenvector(space_op.GetComm(), E, eigenvalues[i], v_comb,
                                              blocks ? &v_comb_lambda : nullptr);
    }

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
        energy_mag < longitudinal_tol * energy_elec)
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
  if (num_combined > 0 && ContinueAdaptation(adapt_it, space_op.GlobalTrueVSize(),
                                             indicator.Norml2(space_op.GetComm())))
  {
    prev_eigenvectors = std::make_unique<VectorFieldTransfer>(
        space_op.GetNDSpace().GetParMesh(), iodata.solver.order);
    for (auto *u : {&v_comb, &v_comb_lambda})
    {
      if (u == &v_comb || blocks)
      {
        prev_eigenvectors->AddField(space_op.GetNDSpace().Get(), u->Real());
        prev_eigenvectors->AddField(space_op.GetNDSpace().Get(), u->Imag());
      }
    }
    prev_num_eigenvectors = num_combined;
  }
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
