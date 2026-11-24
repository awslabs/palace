// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "blockjacobi.hpp"

#include <cmath>
#include <type_traits>
#include <Eigen/Dense>
#include <fmt/format.h>
#include "fem/fespace.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Spectral radius of M⁻¹A for a relaxation M⁻¹ which is not definite (the blocks carry
// signs, and the operator can be complex), by the Arnoldi method: the largest modulus of
// the eigenvalues of the Hessenberg matrix of M⁻¹A on a Krylov space of at most MAX_BASIS
// vectors, restarted from the Ritz vector of that eigenvalue until the estimate changes
// by less than the relative tolerance between two checks. For a definite M and A this is
// the Lanczos estimate of the point Jacobi method. The relaxed operator is not Hermitian
// in general, so its radius is not a spectral norm.
//
// The starting vector is the one of the spectral norm estimate of the point Jacobi
// method (linalg::SetSpectralStartVector): random, or x0 with a small random component if
// x0 has the size of the operator, so that the Ritz vector of a previous operator, which
// can be an exact eigenvector of this one without being the dominant one, does not
// confine the Krylov space. On return x0 is the Ritz vector found.
template <typename OperType, typename VecType, typename RelaxType>
double SpectralRadius(MPI_Comm comm, const OperType &A, const RelaxType &relax, VecType &x0)
{
  constexpr int MAX_BASIS = 20, MAX_APPLICATIONS = 100, CHECK = 5;
  constexpr double TOL = 1.0e-3, BREAKDOWN_TOL = 1.0e-12;
  const int n = A.Height();
  std::vector<VecType> V;
  VecType w(n), t(n);
  w.UseDevice(true);
  t.UseDevice(true);
  V.emplace_back(n);
  V.back().UseDevice(true);
  linalg::SetSpectralStartVector(comm, &x0, V.back());
  x0.SetSize(n);
  x0.UseDevice(true);
  double rho = 0.0, rho_check = 0.0;
  int applications = 0;
  for (bool done = false; !done;)
  {
    // One cycle from the vector in V. The eigenvalues of the Hessenberg matrix are needed
    // at the checks, when the basis is full, and at a breakdown (an invariant subspace).
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(MAX_BASIS + 1, MAX_BASIS);
    Eigen::VectorXcd ritz;
    for (int m = 0; m < MAX_BASIS && !done; m++)
    {
      A.Mult(V[m], t);
      relax(t, w);
      applications++;
      for (int i = 0; i <= m; i++)
      {
        const auto h = linalg::Dot(comm, w, V[i]);
        H(i, m) = h;
        linalg::AXPY(-h, V[i], w);
      }
      const double beta = linalg::Norml2(comm, w);
      H(m + 1, m) = beta;
      const bool breakdown = !(beta > BREAKDOWN_TOL * H.col(m).head(m + 2).norm());
      const bool check = (applications % CHECK == 0);
      const bool last = (m + 1 == MAX_BASIS || applications >= MAX_APPLICATIONS);
      if (breakdown || check || last)
      {
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H.topLeftCorner(m + 1, m + 1));
        int dominant = 0;
        rho = es.eigenvalues().cwiseAbs().maxCoeff(&dominant);
        ritz = es.eigenvectors().col(dominant);
        done = breakdown || !std::isfinite(rho) || applications >= MAX_APPLICATIONS;
        if (check && !done)
        {
          done = (applications >= 2 * CHECK && std::abs(rho - rho_check) <= TOL * rho);
          rho_check = rho;
        }
      }
      if (!done && m + 1 < MAX_BASIS)
      {
        V.emplace_back(w);
        V.back() *= 1.0 / beta;
      }
    }
    // The Ritz vector of the dominant eigenvalue (its real part for a real operator), from
    // which the next cycle starts.
    x0 = 0.0;
    for (int i = 0; i < static_cast<int>(ritz.size()); i++)
    {
      if constexpr (std::is_same_v<VecType, ComplexVector>)
      {
        linalg::AXPY(ritz(i), V[i], x0);
      }
      else
      {
        linalg::AXPY(ritz(i).real(), V[i], x0);
      }
    }
    const double norm = linalg::Norml2(comm, x0);
    if (!done && norm > 0.0 && std::isfinite(norm))
    {
      V.resize(1);
      V[0] = x0;
      V[0] *= 1.0 / norm;
    }
    else
    {
      done = true;
    }
  }
  return rho;
}

// What became of the blocks, for the output.
std::string Describe(const AggregateBlocks::Census &census)
{
  std::string text;
  HYPRE_BigInt num_blocks = 0;
  for (const auto num : census.blocks)
  {
    num_blocks += num;
  }
  if (num_blocks > 0)
  {
    text = fmt::format("{:d} blocks of 1/2/4/8/16/32/64 unknowns: "
                       "{:d}/{:d}/{:d}/{:d}/{:d}/{:d}/{:d} ({:d} unknowns, {:.2f} MB, at "
                       "most {:.2f} MB on a process",
                       num_blocks, census.blocks[0], census.blocks[1], census.blocks[2],
                       census.blocks[3], census.blocks[4], census.blocks[5],
                       census.blocks[6], census.block_dofs, 8.0e-6 * census.entries,
                       8.0e-6 * census.max_entries);
    if (census.split > 0)
    {
      text +=
          fmt::format(", {:d} group{} split", census.split, (census.split == 1) ? "" : "s");
    }
    if (census.gathered > 0)
    {
      text += fmt::format(", {:d} gathered from several processes", census.gathered);
    }
    if (census.dropped > 0)
    {
      text += fmt::format(", {:d} singular directions", census.dropped);
    }
    if (census.negative > 0)
    {
      text += fmt::format(", {:d} negative directions", census.negative);
    }
    text += ")";
  }
  if (census.sheets > 0)
  {
    text += fmt::format("{}{:d} large group{} solved approximately ({:d} unknowns, degree "
                        "{:d}{}, {:.2f} MB)",
                        text.empty() ? "" : ", ", census.sheets,
                        (census.sheets == 1) ? "" : "s", census.sheet_dofs,
                        census.sheet_degree, census.sheet_limited ? " at its limit" : "",
                        12.0e-6 * census.sheet_entries);
  }
  if (census.skipped > 0)
  {
    text += fmt::format("{}{:d} large group{} left to point Jacobi ({:d} unknowns)",
                        text.empty() ? "" : ", ", census.skipped,
                        (census.skipped == 1) ? "" : "s", census.skipped_dofs);
  }
  return text.empty() ? "no blocks (no flat elements), point Jacobi" : text;
}

}  // namespace

template <OperatorType OperType>
BlockJacobiSmoother<OperType>::BlockJacobiSmoother(MPI_Comm comm, int block_size,
                                                   int cheby_order, bool print)
  : JacobiSmoother<OperType>(comm), block_size(block_size),
    sheet_reduction(1.0 / (2.0 * cheby_order + 1.0)), print(print)
{
}

template <OperatorType OperType>
BlockJacobiSmoother<OperType>::~BlockJacobiSmoother() = default;

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::SetOperator(const OperType &op)
{
  // Point Jacobi scaling for the true DOFs in no block.
  JacobiSmoother<OperType>::SetOperator(op);
  blocks.reset();

  // The operator provides what the blocks are constructed from: its real part, and its
  // definite counterpart if it has one.
  const ParOperator *PtAP = nullptr, *PtAP_def = nullptr;
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    const auto *PtAP_complex = dynamic_cast<const ComplexParOperator *>(&op);
    PtAP = dynamic_cast<const ParOperator *>(op.Real());
    MFEM_VERIFY(PtAP_complex && PtAP,
                "BlockJacobiSmoother requires a ComplexParOperator operator with a real "
                "part!");
    PtAP_def = PtAP_complex->GetDefiniteCounterpart();
  }
  else
  {
    PtAP = dynamic_cast<const ParOperator *>(&op);
    MFEM_VERIFY(PtAP, "BlockJacobiSmoother requires a ParOperator operator!");
    PtAP_def = PtAP->GetDefiniteCounterpart();
  }
  // The blocks are constructed for 3D meshes: point Jacobi otherwise.
  std::string result = "no blocks (not a 3D mesh), point Jacobi";
  if (AggregateBlocks::IsSupported(*PtAP))
  {
    blocks =
        std::make_unique<AggregateBlocks>(*PtAP, PtAP_def, block_size, sheet_reduction);
    const auto &census = blocks->GetCensus();
    result = Describe(census);
    if (census.block_dofs == 0 && census.sheets == 0)
    {
      blocks.reset();
    }
  }
  if (blocks)
  {
    // The phases and the sign of the diagonal on the DOFs of the blocks and sheets. On the
    // DOFs of the blocks the point scaling is then set to zero, so that the block product
    // replaces it instead of being added to it and its negative: the reciprocal of a small
    // diagonal entry (which a block of an indefinite operator can have, without being
    // close to singular) would cancel the block product in that sum. A point row requires
    // a finite nonzero diagonal entry, as in the point Jacobi method; a block row does not:
    // if its entry is zero or not finite it has no phase, and the real block is applied as
    // it is.
    const int n = op.Height();
    for (Vector *v : {&wr, &wi, &zr, &zi})
    {
      v->SetSize(n);
      v->UseDevice(false);
      *v = 0.0;
    }
    psi.clear();
    psi_sheet.clear();
    sign_sheet = 1.0;
    if constexpr (std::is_same_v<OperType, ComplexOperator>)
    {
      auto *DR = dinv.Real().HostReadWrite();
      auto *DI = dinv.Imag().HostReadWrite();
      for (int i : blocks->GetSheetDofs())
      {
        const std::complex<double> d(DR[i], DI[i]);
        psi_sheet.push_back(std::sqrt(d / std::abs(d)));
      }
      for (int i : blocks->GetBlockDofs())
      {
        // The sign of the real part is in the block inverse.
        const std::complex<double> d = (DR[i] < 0.0) ? -std::complex<double>(DR[i], DI[i])
                                                     : std::complex<double>(DR[i], DI[i]);
        const double abs_d = std::abs(d);
        psi.push_back((abs_d > 0.0 && std::isfinite(abs_d)) ? std::sqrt(d / abs_d) : 1.0);
        DR[i] = DI[i] = 0.0;
      }
    }
    else
    {
      // The sheets are those of a definite matrix, with the sign of the operator on them.
      auto *D = dinv.HostReadWrite();
      double sum = 0.0;
      for (int i : blocks->GetSheetDofs())
      {
        sum += (D[i] < 0.0) ? -1.0 : 1.0;
      }
      Mpi::GlobalSum(1, &sum, comm);
      sign_sheet = (sum < 0.0) ? -1.0 : 1.0;
      for (int i : blocks->GetBlockDofs())
      {
        D[i] = 0.0;
      }
    }
  }
  if (print && result != summary)
  {
    const auto &fespace = PtAP->TrialFiniteElementSpace();
    Mpi::Print(comm, " Block Jacobi smoother ({}, {:d} unknowns): {}\n",
               fespace.GetFEColl().Name(), fespace.GlobalTrueVSize(), result);
    summary = result;
  }
}

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::AddBlockMult(double s, const VecType &r,
                                                 VecType &d) const
{
  const auto &dofs = blocks->GetBlockDofs();
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // d += s Ψ M_r⁻¹ Ψ r on the DOFs of the blocks.
    const auto *RR = r.Real().HostRead(), *RI = r.Imag().HostRead();
    {
      auto *WR = wr.HostReadWrite(), *WI = wi.HostReadWrite();
      for (std::size_t k = 0; k < dofs.size(); k++)
      {
        const int i = dofs[k];
        const std::complex<double> w = psi[k] * std::complex<double>(RR[i], RI[i]);
        WR[i] = w.real();
        WI[i] = w.imag();
      }
    }
    blocks->MultInverse(wr, zr);
    blocks->MultInverse(wi, zi);
    const auto *ZR = zr.HostRead(), *ZI = zi.HostRead();
    auto *YR = d.Real().HostReadWrite(), *YI = d.Imag().HostReadWrite();
    for (std::size_t k = 0; k < dofs.size(); k++)
    {
      const int i = dofs[k];
      const std::complex<double> z = psi[k] * std::complex<double>(ZR[i], ZI[i]);
      YR[i] += s * z.real();
      YI[i] += s * z.imag();
    }
  }
  else
  {
    blocks->MultInverse(r, zr);
    const auto *Z = zr.HostRead();
    auto *Y = d.HostReadWrite();
    for (int i : dofs)
    {
      Y[i] += s * Z[i];
    }
  }
}

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::Mult(const VecType &x, VecType &y) const
{
  JacobiSmoother<OperType>::Mult(x, y);
  if (blocks)
  {
    AddBlockMult(1.0, x, y);
  }
}

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::ScaledMult(double s, const VecType &r, VecType &d) const
{
  JacobiSmoother<OperType>::ScaledMult(s, r, d);
  if (blocks)
  {
    AddBlockMult(s, r, d);
  }
}

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::ScaledAddMult(double sd, double sr, const VecType &r,
                                                  VecType &d) const
{
  JacobiSmoother<OperType>::ScaledAddMult(sd, sr, r, d);
  if (blocks)
  {
    AddBlockMult(sr, r, d);
  }
}

template <OperatorType OperType>
double BlockJacobiSmoother<OperType>::EstimateMaxEigenvalue(const OperType &A,
                                                            VecType &x0) const
{
  if (!blocks)
  {
    return JacobiSmoother<OperType>::EstimateMaxEigenvalue(A, x0);
  }
  auto relax = [this](const VecType &r, VecType &z) { this->ScaledMult(1.0, r, z); };
  return SpectralRadius(comm, A, relax, x0);
}

template <OperatorType OperType>
bool BlockJacobiSmoother<OperType>::HasSweepCorrection() const
{
  return blocks && blocks->HasSheets();
}

template <OperatorType OperType>
void BlockJacobiSmoother<OperType>::AddSweepCorrection(const VecType &r, VecType &y) const
{
  const auto &dofs = blocks->GetSheetDofs();
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    // y += Ψ S Ψ r on the DOFs of the sheets.
    const auto *RR = r.Real().HostRead(), *RI = r.Imag().HostRead();
    {
      auto *WR = wr.HostReadWrite(), *WI = wi.HostReadWrite();
      auto *ZR = zr.HostReadWrite(), *ZI = zi.HostReadWrite();
      for (std::size_t k = 0; k < dofs.size(); k++)
      {
        const int i = dofs[k];
        const std::complex<double> w = psi_sheet[k] * std::complex<double>(RR[i], RI[i]);
        WR[i] = w.real();
        WI[i] = w.imag();
        ZR[i] = ZI[i] = 0.0;
      }
    }
    blocks->AddSheetSolve(wr, zr);
    blocks->AddSheetSolve(wi, zi);
    const auto *ZR = zr.HostRead(), *ZI = zi.HostRead();
    auto *YR = y.Real().HostReadWrite(), *YI = y.Imag().HostReadWrite();
    for (std::size_t k = 0; k < dofs.size(); k++)
    {
      const int i = dofs[k];
      const std::complex<double> z = psi_sheet[k] * std::complex<double>(ZR[i], ZI[i]);
      YR[i] += z.real();
      YI[i] += z.imag();
    }
  }
  else if (sign_sheet > 0.0)
  {
    blocks->AddSheetSolve(r, y);
  }
  else
  {
    auto *Z = zr.HostReadWrite();
    for (int i : dofs)
    {
      Z[i] = 0.0;
    }
    blocks->AddSheetSolve(r, zr);
    auto *Y = y.HostReadWrite();
    for (int i : dofs)
    {
      Y[i] -= Z[i];
    }
  }
}

template class BlockJacobiSmoother<Operator>;
template class BlockJacobiSmoother<ComplexOperator>;

}  // namespace palace
