// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <chrono>
#include <mfem.hpp>
#include "linalg/orthog.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

inline void ProjectMatInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexOperator &A, Eigen::MatrixXcd &Ar,
                               ComplexVector &r, int n0)
{
  // Update Ar = Vᴴ A V for the new basis dimension n0 -> n. V is real and thus the result
  // is complex symmetric if A is symmetric (which we assume is the case). Ar is replicated
  // across all processes as a sequential n x n matrix.
  const auto n = Ar.rows();
  MFEM_VERIFY(n0 < n, "Unexpected dimensions in PROM matrix projection!");
  for (int j = n0; j < n; j++)
  {
    // Fill block of Vᴴ A V = [  | Vᴴ A vj ] . We can optimize the matrix-vector product
    // since the columns of V are real.
    MFEM_VERIFY(A.HasReal() || A.HasImag(),
                "Invalid zero ComplexOperator for PROM matrix projection!");
    if (A.HasReal())
    {
      A.Real()->Mult(V[j], r.Real());
    }
    if (A.HasImag())
    {
      A.Imag()->Mult(V[j], r.Imag());
    }
    for (int i = 0; i < n; i++)
    {
      Ar(i, j).real(A.HasReal() ? V[i] * r.Real() : 0.0);  // Local inner product
      Ar(i, j).imag(A.HasImag() ? V[i] * r.Imag() : 0.0);
    }
  }
  Mpi::GlobalSum((n - n0) * n, Ar.data() + n0 * n, comm);

  // Fill lower block of Vᴴ A V = [ ____________  |  ]
  //                              [ vjᴴ A V[1:n0] |  ] .
  for (int j = 0; j < n0; j++)
  {
    for (int i = n0; i < n; i++)
    {
      Ar(i, j) = Ar(j, i);
    }
  }
}

inline void ProjectVecInternal(MPI_Comm comm, const std::vector<Vector> &V,
                               const ComplexVector &b, Eigen::VectorXcd &br, int n0)
{
  // Update br = Vᴴ b for the new basis dimension n0 -> n. br is replicated across all
  // processes as a sequential n-dimensional vector.
  const auto n = br.size();
  MFEM_VERIFY(n0 < n, "Unexpected dimensions in PROM vector projection!");
  for (int i = n0; i < n; i++)
  {
    br(i).real(V[i] * b.Real());  // Local inner product
    br(i).imag(V[i] * b.Imag());
  }
  Mpi::GlobalSum(n - n0, br.data() + n0, comm);
}

}  // namespace

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &spaceop) : spaceop(spaceop)
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogenous, no special elimination is required on
  // the RHS. The damping matrix may be nullptr.
  K = spaceop.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  C = spaceop.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  M = spaceop.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices when constructing PROM!");

  // Set up RHS vector (linear in frequency part) for the incident field at port boundaries,
  // and the vector for the solution, which satisfies the Dirichlet (PEC) BC.
  if (!spaceop.GetExcitationVector1(RHS1))
  {
    RHS1.SetSize(0);
  }
  has_A2 = has_RHS2 = true;

  // Initialize temporary vector storage.
  r.SetSize(K->Height());
  w.SetSize(K->Height());

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  ksp = std::make_unique<ComplexKspSolver>(iodata, spaceop.GetNDSpaces(),
                                           &spaceop.GetH1Spaces());

  // Initialize solver for inner product solves. The system matrix for the inner product is
  // real and SPD. This uses the dual norm from https://ieeexplore.ieee.org/document/5313818
  // in the error estimate.
  if (iodata.solver.driven.adaptive_metric_aposteriori)
  {
    constexpr int curlcurl_verbose = 0;
    kspKM = std::make_unique<WeightedHCurlNormSolver>(
        spaceop.GetMaterialOp(), spaceop.GetNDSpaces(), spaceop.GetH1Spaces(),
        spaceop.GetNDDbcTDofLists(), spaceop.GetH1DbcTDofLists(), iodata.solver.linear.tol,
        iodata.solver.linear.max_it, curlcurl_verbose, iodata.solver.pa_order_threshold,
        iodata.solver.pa_discrete_interp);
  }

  // The initial PROM basis is empty. Orthogonalization uses MGS by default, else CGS2.
  dim_V = 0;
  orthog_mgs =
      (iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::MGS);

  // Seed the random number generator for parameter space sampling.
  engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
}

void RomOperator::Initialize(double start, double delta, int num_steps, int max_dim)
{
  // Initialize P = {ω_L, ω_L+δ, ..., ω_R}. Always insert in ascending order.
  MFEM_VERIFY(PS.empty() && P_m_PS.empty(),
              "RomOperator::Initialize should only be called once!");
  MFEM_VERIFY(
      num_steps > 2,
      "RomOperator adaptive frequency sweep should have more than two frequency steps!");
  if (delta < 0.0)
  {
    start = start + (num_steps - 1) * delta;
    delta = -delta;
  }
  auto it = P_m_PS.begin();
  for (int step = 0; step < num_steps; step++)
  {
    it = P_m_PS.emplace_hint(it, start + step * delta);
  }

  // PROM operators Ar = Vᴴ A V when assembled is complex symmetric for real V. The provided
  // max_dim is the number of sample points (2 basis vectors per point).
  MFEM_VERIFY(max_dim > 0, "Reduced order basis storage must have > 0 columns!");
  V.resize(2 * max_dim, Vector());
}

void RomOperator::SolveHDM(double omega, ComplexVector &e)
{
  // Compute HDM solution at the given frequency. The system matrix, A = K + iω C - ω² M +
  // A2(ω) is built by summing the underlying operator contributions.
  BlockTimer bt0(Timer::CONSTRUCT);
  A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
  has_A2 = (A2 != nullptr);
  auto A = spaceop.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                   std::complex<double>(-omega * omega, 0.0), K.get(),
                                   C.get(), M.get(), A2.get());
  auto P =
      spaceop.GetPreconditionerMatrix<ComplexOperator>(1.0, omega, -omega * omega, omega);
  ksp->SetOperators(*A, *P);

  // The HDM excitation vector is computed as RHS = iω RHS1 + RHS2(ω).
  Mpi::Print("\n");
  if (has_RHS2)
  {
    has_RHS2 = spaceop.GetExcitationVector2(omega, r);
  }
  else
  {
    r = 0.0;
  }
  if (RHS1.Size())
  {
    r.Add(1i * omega, RHS1);
  }

  // Solve the linear system.
  BlockTimer bt1(Timer::SOLVE);
  ksp->Mult(r, e);
}

void RomOperator::AddHDMSample(double omega, ComplexVector &e)
{
  // Use the given HDM solution at the given frequency to update the reduced-order basis
  // updating the PROM operators.
  auto it = P_m_PS.lower_bound(omega);
  MFEM_VERIFY(it != P_m_PS.end(),
              "Sample frequency " << omega << " not found in parameter set!");
  P_m_PS.erase(it);
  auto ret = PS.insert(omega);
  MFEM_VERIFY(ret.second, "Sample frequency "
                              << omega << " already exists in the sampled parameter set!");

  // Update V. The basis is always real (each complex solution adds two basis vectors if it
  // has a nonzero real and imaginary parts).
  const double normr = linalg::Norml2(spaceop.GetComm(), e.Real());
  const double normi = linalg::Norml2(spaceop.GetComm(), e.Imag());
  const bool has_real = (normr > 1.0e-12 * std::sqrt(normr * normr + normi * normi));
  const bool has_imag = (normi > 1.0e-12 * std::sqrt(normr * normr + normi * normi));
  MFEM_VERIFY(dim_V + has_real + has_imag <= static_cast<int>(V.size()),
              "Unable to increase basis storage size, increase maximum number of vectors!");
  const int dim_V0 = dim_V;
  std::vector<double> H(dim_V + 1);
  if (has_real)
  {
    V[dim_V] = e.Real();
    if (orthog_mgs)
    {
      linalg::OrthogonalizeColumnMGS(spaceop.GetComm(), V, V[dim_V], H.data(), dim_V);
    }
    else
    {
      linalg::OrthogonalizeColumnCGS(spaceop.GetComm(), V, V[dim_V], H.data(), dim_V, true);
    }
    V[dim_V] *= 1.0 / linalg::Norml2(spaceop.GetComm(), V[dim_V]);
    dim_V++;
  }
  if (has_imag)
  {
    V[dim_V] = e.Imag();
    if (orthog_mgs)
    {
      linalg::OrthogonalizeColumnMGS(spaceop.GetComm(), V, V[dim_V], H.data(), dim_V);
    }
    else
    {
      linalg::OrthogonalizeColumnCGS(spaceop.GetComm(), V, V[dim_V], H.data(), dim_V, true);
    }
    V[dim_V] *= 1.0 / linalg::Norml2(spaceop.GetComm(), V[dim_V]);
    dim_V++;
  }

  // Update reduced-order operators. Resize preserves the upper dim0 x dim0 block of each
  // matrix and first dim0 entries of each vector and the projection uses the values
  // computed for the unchanged basis vectors.
  Kr.conservativeResize(dim_V, dim_V);
  ProjectMatInternal(spaceop.GetComm(), V, *K, Kr, r, dim_V0);
  if (C)
  {
    Cr.conservativeResize(dim_V, dim_V);
    ProjectMatInternal(spaceop.GetComm(), V, *C, Cr, r, dim_V0);
  }
  Mr.conservativeResize(dim_V, dim_V);
  ProjectMatInternal(spaceop.GetComm(), V, *M, Mr, r, dim_V0);
  Ar.resize(dim_V, dim_V);
  if (RHS1.Size())
  {
    RHS1r.conservativeResize(dim_V);
    ProjectVecInternal(spaceop.GetComm(), V, RHS1, RHS1r, dim_V0);
  }
  RHSr.resize(dim_V);
}

void RomOperator::AssemblePROM(double omega)
{
  // Assemble the PROM linear system at the given frequency. The PROM system is defined by
  // the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ - ω² Mᵣ + Vᴴ A2 V(ω) and source vector RHSᵣ(ω) =
  // iω RHS1ᵣ + Vᴴ RHS2(ω). A2(ω) and RHS2(ω) are constructed only if required and are
  // only nonzero on boundaries, will be empty if not needed.
  if (has_A2)
  {
    A2 = spaceop.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    ProjectMatInternal(spaceop.GetComm(), V, *A2, Ar, r, 0);
  }
  else
  {
    Ar.setZero();
  }
  Ar += Kr;
  if (C)
  {
    Ar += (1i * omega) * Cr;
  }
  Ar += (-omega * omega) * Mr;

  if (has_RHS2)
  {
    spaceop.GetExcitationVector2(omega, RHS2);
    ProjectVecInternal(spaceop.GetComm(), V, RHS2, RHSr, 0);
  }
  else
  {
    RHSr.setZero();
  }
  if (RHS1.Size())
  {
    RHSr += (1i * omega) * RHS1r;
  }
}

void RomOperator::SolvePROM(ComplexVector &e)
{
  // Compute PROM solution at the given frequency and expand into high-dimensional space.
  // The PROM is solved on every process so the matrix-vector product for vector expansion
  // does not require communication.
  RHSr = Ar.partialPivLu().solve(RHSr);
  // RHSr = Ar.ldlt().solve(RHSr);
  // RHSr = Ar.selfadjointView<Eigen::Lower>().ldlt().solve(RHSr);

  e = 0.0;
  for (int j = 0; j < dim_V; j++)
  {
    e.Real().Add(RHSr(j).real(), V[j]);
    e.Imag().Add(RHSr(j).imag(), V[j]);
  }
}

double RomOperator::ComputeError(double omega)
{
  // Compute the error metric associated with the approximate PROM solution at the given
  // frequency. The HDM residual -r = [K + iω C - ω² M + A2(ω)] x - [iω RHS1 + RHS2(ω)] is
  // computed using the most recently computed A2(ω) and RHS2(ω).
  AssemblePROM(omega);
  SolvePROM(w);

  // Residual error.
  r = 0.0;
  if (RHS1.Size())
  {
    r.Add(-1i * omega, RHS1);
  }
  if (has_RHS2)
  {
    r.Add(-1.0, RHS2);
  }
  double den = !kspKM ? linalg::Norml2(spaceop.GetComm(), r) : 0.0;

  K->AddMult(w, r, 1.0);
  if (C)
  {
    C->AddMult(w, r, 1i * omega);
  }
  M->AddMult(w, r, -omega * omega);
  if (has_A2)
  {
    A2->AddMult(w, r, 1.0);
  }

  double num;
  if (!kspKM)
  {
    num = linalg::Norml2(spaceop.GetComm(), r);
  }
  else
  {
    z.SetSize(r.Size());
    kspKM->Mult(r, z);
    auto dot = linalg::Dot(spaceop.GetComm(), z, r);
    MFEM_ASSERT(dot.real() > 0.0 && std::abs(dot.imag()) < 1.0e-9 * dot.real(),
                "Non-positive vector norm in normalization (dot = " << dot << ")!");
    num = std::sqrt(dot.real());
    den = linalg::Norml2(spaceop.GetComm(), w, kspKM->GetOperator(), z);
  }
  MFEM_VERIFY(den > 0.0, "Unexpected zero denominator in HDM residual!");
  return num / den;
}

double RomOperator::ComputeMaxError(int num_cand, double &omega_star)
{
  // Greedy iteration: Find argmax_{ω ∈ P_C} η(e; ω). We sample num_cand candidates from
  // P \ P_S.
  num_cand = std::min(num_cand, static_cast<int>(P_m_PS.size()));
  std::vector<double> PC;
  if (Mpi::Root(spaceop.GetComm()))
  {
    if constexpr (false)
    {
      // Sample with weighted probability by distance from the set of already sampled
      // points.
      std::vector<double> weights(P_m_PS.size());
      PC.reserve(num_cand);
      for (auto sample : PS)
      {
        int i = std::distance(P_m_PS.begin(), P_m_PS.lower_bound(sample));
        int il = i - 1;
        while (il >= 0)
        {
          weights[il] = std::min(weights[il], static_cast<double>(i - il));
          il--;
        }
        int iu = i;
        while (iu < weights.size())
        {
          weights[iu] = std::min(weights[iu], static_cast<double>(1 + iu - i));
          iu++;
        }
      }
      for (int i = 0; i < num_cand; i++)
      {
        std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
        auto res = dist(engine);
        auto it = P_m_PS.begin();
        std::advance(it, res);
        PC.push_back(*it);
        weights[res] = 0.0;  // No replacement
      }
    }
    else
    {
      // Sample with uniform probability.
      PC.reserve(num_cand);
      std::sample(P_m_PS.begin(), P_m_PS.end(), std::back_inserter(PC), num_cand, engine);
    }
  }
  else
  {
    PC.resize(num_cand);
  }
  Mpi::Broadcast(num_cand, PC.data(), 0, spaceop.GetComm());

  // Debug
  // Mpi::Print("Candidate sampling:\n");
  // Mpi::Print(" P_S: {}", PS);
  // Mpi::Print(" P\\P_S: {}\n", P_m_PS);
  // Mpi::Print(" P_C: {}\n", PC);
  // Mpi::Print("\n");

  // For each candidate, compute the PROM solution and associated error metric.
  double err_max = 0.0;
  for (auto omega : PC)
  {
    double err = ComputeError(omega);

    // Debug
    // Mpi::Print("ω = {:.3e}, error = {:.3e}\n", omega, err);

    if (err > err_max)
    {
      err_max = err;
      omega_star = omega;
    }
  }
  return err_max;
}

}  // namespace palace
