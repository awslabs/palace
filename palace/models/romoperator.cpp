// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "romoperator.hpp"

#include <algorithm>
#include <chrono>
#include "fem/freqdomain.hpp"
#include "fem/operator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

RomOperator::RomOperator(const IoData &iodata, SpaceOperator &sp, int nmax)
  : spaceop(sp),
    engine((unsigned)std::chrono::system_clock::now().time_since_epoch().count())
{
  // Construct the system matrices defining the linear operator. PEC boundaries are handled
  // simply by setting diagonal entries of the system matrix for the corresponding dofs.
  // Because the Dirichlet BC is always homogenous, no special elimination is required on
  // the RHS. The damping matrix may be nullptr.
  K = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::STIFFNESS,
                                   mfem::Operator::DIAG_ONE);
  M = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::MASS,
                                   mfem::Operator::DIAG_ZERO);
  C = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::DAMPING,
                                   mfem::Operator::DIAG_ZERO);

  // Set up the linear solver and set operators but don't set the operators yet (this will
  // be done during an HDM solve at a given parameter point). The preconditioner for the
  // complex linear system is constructed from a real approximation to the complex system
  // matrix.
  pc0 = std::make_unique<KspPreconditioner>(iodata, spaceop.GetDbcMarker(),
                                            spaceop.GetNDSpaces(), &spaceop.GetH1Spaces());
  ksp0 = std::make_unique<KspSolver>(K->GetComm(), iodata, "ksp_");
  ksp0->SetPreconditioner(*pc0);

  // Set up RHS vector (linear in frequency part) for the incident field at port boundaries,
  // and the vector for the solution, which satisfies the Dirichlet (PEC) BC.
  RHS1 = std::make_unique<petsc::PetscParVector>(*K);
  if (!spaceop.GetFreqDomainExcitationVector1(*RHS1))
  {
    RHS1.reset();
  }
  init2 = true;
  hasA2 = hasRHS2 = false;

  // Initialize other data structure and storage.
  E0 = std::make_unique<petsc::PetscParVector>(*K);
  R0 = std::make_unique<petsc::PetscParVector>(*K);
  T0 = std::make_unique<petsc::PetscParVector>(*K);

  // Initialize solver for inner product solves. The system matrix for the inner product is
  // real and SPD. This uses the dual norm from https://ieeexplore.ieee.org/document/5313818
  // in the error estimate.
  if (iodata.solver.driven.adaptive_metric_aposteriori)
  {
    constexpr int curlcurl_verbose = 0;
    kspKM = std::make_unique<CurlCurlSolver>(
        spaceop.GetMaterialOp(), spaceop.GetDbcMarker(), spaceop.GetNDSpaces(),
        spaceop.GetH1Spaces(), iodata.solver.linear.tol, iodata.solver.linear.max_it,
        curlcurl_verbose);

    auto KM = std::make_unique<SumOperator>(K->GetNumRows(), K->GetNumCols());
    KM->AddOperator(*K->GetOperator(petsc::PetscParMatrix::ExtractStructure::REAL));
    KM->AddOperator(*M->GetOperator(petsc::PetscParMatrix::ExtractStructure::REAL));
    opKM = std::make_unique<petsc::PetscShellMatrix>(K->GetComm(), std::move(KM));
    opKM->SetRealSymmetric();
  }

  // Construct initial (empty) basis and ROM operators. Ar = Vᴴ A V when assembled is
  // complex symmetric for real V. The provided nmax is the number of sample points(2 basis
  // vectors per point).
  MFEM_VERIFY(K && M, "Invalid empty HDM matrices constructing PROM operators!");
  MFEM_VERIFY(nmax > 0, "Reduced order basis storage must have > 0 columns!");
  dim = 0;
  omega_min = delta_omega = 0.0;
  V = std::make_unique<petsc::PetscDenseMatrix>(K->GetComm(), K->Height(), PETSC_DECIDE,
                                                PETSC_DECIDE, 2 * nmax, nullptr);

  Kr = std::make_unique<petsc::PetscDenseMatrix>(dim, dim, nullptr);
  Kr->CopySymmetry(*K);
  Mr = std::make_unique<petsc::PetscDenseMatrix>(dim, dim, nullptr);
  Mr->CopySymmetry(*M);
  if (C)
  {
    Cr = std::make_unique<petsc::PetscDenseMatrix>(dim, dim, nullptr);
    Cr->CopySymmetry(*C);
  }
  else
  {
    Cr = nullptr;
  }
  Ar = std::make_unique<petsc::PetscDenseMatrix>(dim, dim, nullptr);
  Ar->SetSymmetric(K->GetSymmetric() && M->GetSymmetric() && (!C || C->GetSymmetric()));

  RHS1r = (RHS1) ? std::make_unique<petsc::PetscParVector>(*Ar) : nullptr;
  RHSr = std::make_unique<petsc::PetscParVector>(*Ar);
  Er = std::make_unique<petsc::PetscParVector>(*Ar);

  // Set up the linear solver (dense sequential on all processors). An indefinite LDLᵀ
  // factorization is used when Ar has its symmetry flag set. The default sequential dense
  // matrix uses LAPACK for the factorization.
  int print = 0;
  ksp = std::make_unique<KspSolver>(Ar->GetComm(), print, "rom_");
  ksp->SetType(KspSolver::Type::CHOLESKY);  // Symmetric indefinite factorization
}

void RomOperator::Initialize(int steps, double start, double delta)
{
  // Initialize P = {ω_L, ω_L+δ, ..., ω_R}. Always insert in ascending order.
  MFEM_VERIFY(Ps.empty(), "RomOperator::Initialize should only be called once!");
  MFEM_VERIFY(steps > 2, "RomOperator adaptive frequency sweep should have more than two "
                         "frequency steps!");
  Ps.reserve(steps);
  PmPs.resize(steps);
  if (delta < 0.0)
  {
    start = start + (steps - 1) * delta;
    delta = -delta;
  }
  for (int step = 0; step < steps; step++)
  {
    PmPs[step] = start + step * delta;
  }
  omega_min = start;
  delta_omega = delta;
  A2.resize(steps);
  RHS2.resize(steps);
}

void RomOperator::SolveHDM(double omega, petsc::PetscParVector &E, bool print)
{
  // Compute HDM solution at the given frequency and add solution to the reduced-order
  // basis, updating the PROM operators. Update P_S and P\P_S sets.
  auto it = std::lower_bound(PmPs.begin(), PmPs.end(), omega);
  MFEM_VERIFY(it != PmPs.end(),
              "Sample frequency " << omega << " not found in parameter set!");
  PmPs.erase(it);
  Ps.push_back(omega);

  // Set up HDM system and solve. The system matrix A = K + iω C - ω² M + A2(ω) is built
  // by summing the underlying operator contributions (to save memory).
  {
    const auto step = std::lround((omega - omega_min) / delta_omega);
    MFEM_VERIFY(step >= 0 && static_cast<std::size_t>(step) < A2.size(),
                "Invalid out-of-range frequency for PROM solution!");
    std::vector<std::unique_ptr<mfem::Operator>> P, AuxP;
    A2[step] = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::EXTRA, omega,
                                            mfem::Operator::DIAG_ZERO, print);
    auto A = utils::GetSystemMatrixShell(omega, *K, *M, C.get(), A2[step].get());
    spaceop.GetPreconditionerMatrix(omega, P, AuxP, print);
    pc0->SetOperator(P, &AuxP);
    ksp0->SetOperator(*A);

    Mpi::Print("\n");
    spaceop.GetFreqDomainExcitationVector(omega, *R0);
    E.SetZero();
    ksp0->Mult(*R0, E);
  }

  double norm = E.Normlinf(), ntol = 1.0e-12;
  mfem::Vector Er_(E.GetSize()), Ei_(E.GetSize());
  E.GetToVectors(Er_, Ei_);
  bool has_real = (std::sqrt(mfem::InnerProduct(E.GetComm(), Er_, Er_)) > ntol * norm);
  bool has_imag = (std::sqrt(mfem::InnerProduct(E.GetComm(), Ei_, Ei_)) > ntol * norm);

  // Update V. The basis is always real (each complex solution adds two basis vectors if it
  // has a nonzero real and imaginary parts).
  PetscInt nmax = V->GetGlobalNumCols(), dim0 = dim;
  dim = (has_real) + (has_imag) + static_cast<int>(dim0);
  MFEM_VERIFY(dim <= nmax, "Unable to increase basis storage size, increase maximum number "
                           "of vectors!");
  bool mgs = false, cgs2 = true;
  if (has_real && has_imag)
  {
    {
      petsc::PetscParVector v = V->GetColumn(dim - 2);
      v.SetFromVector(Er_);
      V->RestoreColumn(dim - 2, v);
      if (opKM)
      {
        V->OrthonormalizeColumn(dim - 2, mgs, cgs2, *opKM, *T0);
      }
      else
      {
        V->OrthonormalizeColumn(dim - 2, mgs, cgs2);
      }
    }
    {
      petsc::PetscParVector v = V->GetColumn(dim - 1);
      v.SetFromVector(Ei_);
      V->RestoreColumn(dim - 1, v);
      if (opKM)
      {
        V->OrthonormalizeColumn(dim - 1, mgs, cgs2, *opKM, *T0);
      }
      else
      {
        V->OrthonormalizeColumn(dim - 1, mgs, cgs2);
      }
    }
  }
  else
  {
    {
      petsc::PetscParVector v = V->GetColumn(dim - 1);
      v.Copy(E);
      V->RestoreColumn(dim - 1, v);
      if (opKM)
      {
        V->OrthonormalizeColumn(dim - 1, mgs, cgs2, *opKM, *T0);
      }
      else
      {
        V->OrthonormalizeColumn(dim - 1, mgs, cgs2);
      }
    }
  }

  // Update reduced-order operators. Resize preserves the upper dim0 x dim0 block of each
  // matrix and first dim0 entries of each vector and the projection uses the values
  // computed for the unchanged basis vectors.
  bool init = (dim0 > 0);
  Kr->Resize(dim, dim, init);
  Mr->Resize(dim, dim, init);
  BVMatProjectInternal(*V, *K, *Kr, *R0, dim0, dim);
  BVMatProjectInternal(*V, *M, *Mr, *R0, dim0, dim);
  if (C)
  {
    Cr->Resize(dim, dim, init);
    BVMatProjectInternal(*V, *C, *Cr, *R0, dim0, dim);
  }
  if (RHS1)
  {
    RHS1r->Resize(dim, init);
    BVDotVecInternal(*V, *RHS1, *RHS1r, dim0, dim);
  }
  Ar->Resize(dim, dim);
  RHSr->Resize(dim);
  Er->Resize(dim);
  if (init)
  {
    ksp->Reset();  // Operator size change
  }
  ksp->SetOperator(*Ar);
}

void RomOperator::AssemblePROM(double omega)
{
  // Assemble the PROM linear system at the given frequency. Do some additional set up at
  // the first solve call. The PROM system is defined by the matrix Aᵣ(ω) = Kᵣ + iω Cᵣ
  // - ω² Mᵣ + Vᴴ A2ᵣ V(ω) and source vector RHSᵣ(ω) = iω RHS1ᵣ + Vᴴ RHS2ᵣ(ω) V.
  const auto step = std::lround((omega - omega_min) / delta_omega);
  MFEM_VERIFY(step >= 0 && static_cast<std::size_t>(step) < A2.size(),
              "Invalid out-of-range frequency for PROM solution!");

  // Construct A2(ω) and RHS2(ω) if required (only nonzero on boundaries, will be empty
  // if not needed).
  if (init2)
  {
    auto tA2 = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::EXTRA, omega,
                                            mfem::Operator::DIAG_ZERO, false);
    if (tA2)
    {
      hasA2 = true;
      A2[step] = std::move(tA2);
    }
    auto tRHS2 = std::make_unique<petsc::PetscParVector>(*K);
    if (spaceop.GetFreqDomainExcitationVector2(omega, *tRHS2))
    {
      hasRHS2 = true;
      RHS2[step] = std::move(tRHS2);
    }
    init2 = false;
  }

  // Set up PROM linear system.
  Ar->Scale(0.0);
  if (hasA2)
  {
    if (!A2[step])
    {
      // Debug
      // Mpi::Print("Inserting cache value for omega = {:e}\n", omega);
      A2[step] = spaceop.GetSystemMatrixPetsc(SpaceOperator::OperatorType::EXTRA, omega,
                                              mfem::Operator::DIAG_ZERO, false);
    }
    else
    {
      // Debug
      // Mpi::Print("Found cache value for omega = {:e} (step = {:d})\n", omega, step);
    }
    BVMatProjectInternal(*V, *A2[step], *Ar, *R0, 0, dim);
  }
  Ar->AXPY(1.0, *Kr, petsc::PetscParMatrix::NNZStructure::SAME);
  Ar->AXPY(-omega * omega, *Mr, petsc::PetscParMatrix::NNZStructure::SAME);
  if (C)
  {
    Ar->AXPY(1i * omega, *Cr, petsc::PetscParMatrix::NNZStructure::SAME);
  }

  RHSr->SetZero();
  if (hasRHS2)
  {
    if (!RHS2[step])
    {
      RHS2[step] = std::make_unique<petsc::PetscParVector>(*K);
      spaceop.GetFreqDomainExcitationVector2(omega, *RHS2[step]);
    }
    BVDotVecInternal(*V, *RHS2[step], *RHSr, 0, dim);
  }
  if (RHS1)
  {
    RHSr->AXPY(1i * omega, *RHS1r);
  }
}

void RomOperator::SolvePROM(petsc::PetscParVector &E)
{
  // Compute PROM solution at the given frequency and expand into high- dimensional space.
  // The PROM is solved on every process so the matrix- vector product for vector expansion
  // is sequential.
  ksp->Mult(*RHSr, *Er);
  {
    PetscScalar *pV = V->GetArray(), *pE = E.GetArray();
    petsc::PetscDenseMatrix locV(V->Height(), dim, pV);
    petsc::PetscParVector locE(V->Height(), pE);
    locV.Mult(*Er, locE);
    V->RestoreArray(pV);
    E.RestoreArray(pE);
  }
}

double RomOperator::ComputeMaxError(int Nc, double &omega_star)
{
  // Greedy iteration: Find argmax_{ω ∈ P_C} η(E; ω). We sample Nc candidates from P \ P_S.
  MPI_Comm comm = K->GetComm();
  Nc = std::min(Nc, static_cast<int>(PmPs.size()));
  std::vector<double> Pc;
  if (Mpi::Root(comm))
  {
    // Sample with uniform probability.
    Pc.reserve(Nc);
    std::sample(PmPs.begin(), PmPs.end(), std::back_inserter(Pc), Nc, engine);

#if 0
    // Sample with weighted probability by distance from the set of already sampled
    // points.
    std::vector<double> weights(PmPs.size());
    weights = static_cast<double>(weights.Size());
    Pc.reserve(Nc);
    for (auto sample : Ps)
    {
      int i = std::distance(PmPs.begin(),
                            std::lower_bound(PmPs.begin(), PmPs.end(), sample));
      int il = i-1;
      while (il >= 0)
      {
        weights[il] = std::min(weights[il], static_cast<double>(i-il));
        il--;
      }
      int iu = i;
      while (iu < weights.size())
      {
        weights[iu] = std::min(weights[iu], static_cast<double>(1+iu-i));
        iu++;
      }
    }
    for (int i = 0; i < Nc; i++)
    {
      std::discrete_distribution<int> dist(weights.begin(), weights.end());
      int res = dist(engine);
      Pc.push_back(PmPs[res]);
      weights[res] = 0.0;  // No replacement
    }
#endif
  }
  else
  {
    Pc.resize(Nc);
  }
  Mpi::Broadcast(Nc, Pc.data(), 0, comm);

  // Debug
  // Mpi::Print("Candidate sampling:\n");
  // Mpi::Print(" P_S: {}", Ps);
  // Mpi::Print(" P\\P_S: {}\n", PmPs);
  // Mpi::Print(" P_C: {}\n", Pc);

  // For each candidate, compute the PROM solution and associated error metric.
  double err_max = 0.0;
  for (auto omega : Pc)
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

double RomOperator::ComputeError(double omega)
{
  // Compute the error metric associated with the approximate PROM solution at the given
  // frequency. The HDM residual R = [K + iω C - ω² M + A2(ω)] x - [iω RHS1 + RHS2(ω)] is
  // computed using the most recently computed A2(ω) and RHS2(ω).
  AssemblePROM(omega);
  SolvePROM(*E0);

  // Residual error.
  const auto step = std::lround((omega - omega_min) / delta_omega);
  MFEM_VERIFY(step >= 0 && static_cast<std::size_t>(step) < A2.size(),
              "Invalid out-of-range frequency for PROM solution!");
  double num, den = 1.0;
  R0->SetZero();
  if (RHS1)
  {
    R0->AXPY(-1i * omega, *RHS1);
  }
  if (hasRHS2)
  {
    MFEM_VERIFY(RHS2[step], "Unexpected uncached frequency for RHS2 vector in PROM!");
    R0->AXPY(-1.0, *RHS2[step]);
  }
  if (!kspKM)
  {
    den = R0->Norml2();
  }

  K->MultAdd(*E0, *R0);
  M->Mult(*E0, *T0);
  R0->AXPY(-omega * omega, *T0);
  if (C)
  {
    C->Mult(*E0, *T0);
    R0->AXPY(1i * omega, *T0);
  }
  if (hasA2)
  {
    MFEM_VERIFY(A2[step], "Unexpected uncached frequency for A2 matrix in PROM!");
    A2[step]->MultAdd(*E0, *R0);
  }
  if (!kspKM)
  {
    num = R0->Norml2();
  }
  else
  {
    kspKM->Mult(*R0, *T0);
    num = std::sqrt(std::real(R0->Dot(*T0)));
    opKM->Mult(*E0, *T0);
    den = std::sqrt(std::real(E0->Dot(*T0)));
  }
  MFEM_VERIFY(den > 0.0, "Unexpected zero denominator in HDM residual!");
  return num / den;
}

// static
void RomOperator::BVMatProjectInternal(petsc::PetscDenseMatrix &V, petsc::PetscParMatrix &A,
                                       petsc::PetscDenseMatrix &Ar,
                                       petsc::PetscParVector &r, int n0, int n)
{
  // Update Ar = Vᴴ A V for the new basis dimension n0 => n. We assume V is real and thus
  // the result is complex symmetric if A is symmetric. Ar is replicated across all
  // processes (sequential n x n matrix).
  MFEM_VERIFY(n0 < n, "Unexpected dimensions in BVMatProjectInternal!");
  MFEM_VERIFY(A.GetSymmetric() && Ar.GetSymmetric(),
              "BVMatProjectInternal is specialized for symmetric matrices!");
  mfem::Vector vr(V.Height());
  for (int j = n0; j < n; j++)
  {
    // Fill block of Vᴴ A V = [  | Vᴴ A vj ] . We optimize matrix-vector product since we
    // know columns of V are real.
    {
      petsc::PetscParVector v = V.GetColumn(j);
      v.GetToVector(vr);
      A.Mult(vr, r);
      // A.Mult(v, r);
      V.RestoreColumn(j, v);
    }
    {
      PetscScalar *pV = V.GetArray(), *pr = r.GetArray(), *pAr = Ar.GetArray();
      petsc::PetscDenseMatrix locV(V.Height(), n, pV);
      petsc::PetscParVector locr(V.Height(), pr), arn(n, pAr + j * n);
      locV.MultTranspose(locr, arn);  // Vᴴ = Vᵀ
      V.RestoreArray(pV);
      r.RestoreArray(pr);
      Ar.RestoreArray(pAr);
    }
  }
  // Fill lower block of Vᴴ A V = [ ____________  |  ]
  //                              [ vjᴴ A V[1:n0] |  ] .
  {
    PetscScalar *pAr = Ar.GetArray();
    Mpi::GlobalSum((n - n0) * n, pAr + n0 * n, V.GetComm());
    for (int j = 0; j < n0; j++)
    {
      for (int i = n0; i < n; i++)
      {
        pAr[i + j * n] = pAr[j + i * n];
      }
    }
    Ar.RestoreArray(pAr);
  }
}

// static
void RomOperator::BVDotVecInternal(petsc::PetscDenseMatrix &V, petsc::PetscParVector &b,
                                   petsc::PetscParVector &br, int n0, int n)
{
  // Update br = Vᴴ b for the new basis dimension n0 => n. br is replicated across all
  // processes (sequential n-dimensional vector).
  MFEM_VERIFY(n0 < n, "Unexpected dimensions in BVDotVecInternal!");
  PetscScalar *pV = V.GetArray(), *pb = b.GetArray(), *pbr = br.GetArray();
  petsc::PetscDenseMatrix locV(V.Height(), n - n0, pV + n0 * V.Height());
  petsc::PetscParVector locb(V.Height(), pb), brn(n - n0, pbr + n0);
  locV.MultTranspose(locb, brn);  // Vᴴ = Vᵀ
  V.RestoreArray(pV);
  b.RestoreArray(pb);
  Mpi::GlobalSum(n - n0, pbr + n0, V.GetComm());
  br.RestoreArray(pbr);
}

}  // namespace palace
