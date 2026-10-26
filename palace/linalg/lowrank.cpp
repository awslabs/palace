// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowrank.hpp"

#include <algorithm>
#include <mfem/general/forall.hpp>

namespace palace
{

template <bool Conj>
LowRankComplexUpdate<Conj>::LowRankComplexUpdate(MPI_Comm comm,
                                                 std::span<const ComplexVector *const> v,
                                                 std::span<const std::complex<double>> g)
  : comm(comm), K(static_cast<int>(v.size())), g(g.begin(), g.end())
{
  MFEM_VERIFY(v.size() == g.size(), "Mismatched number of low-rank vectors and weights!");
  const int n = v.empty() ? 0 : v.front()->Size();

  // Union of the exactly nonzero entries over all terms, so dropping the rest is exact.
  std::vector<bool> nonzero(n, false);
  for (const auto *vk : v)
  {
    MFEM_VERIFY(vk->Size() == n, "Mismatched low-rank vector sizes!");
    const auto *VR = vk->Real().HostRead();
    const auto *VI = vk->Imag().HostRead();
    for (int i = 0; i < n; i++)
    {
      nonzero[i] = nonzero[i] || VR[i] != 0.0 || VI[i] != 0.0;
    }
  }
  const int nnz = static_cast<int>(std::count(nonzero.begin(), nonzero.end(), true));
  idx.SetSize(nnz);
  {
    auto *I = idx.HostWrite();
    for (int i = 0, j = 0; i < n; i++)
    {
      if (nonzero[i])
      {
        I[j++] = i;
      }
    }
  }

  vr.SetSize(K * nnz);
  vi.SetSize(K * nnz);
  vr.UseDevice(true);
  vi.UseDevice(true);
  {
    const auto *I = idx.HostRead();
    auto *PR = vr.HostWrite();
    auto *PI = vi.HostWrite();
    for (int k = 0; k < K; k++)
    {
      const auto *VR = v[k]->Real().HostRead();
      const auto *VI = v[k]->Imag().HostRead();
      for (int j = 0; j < nnz; j++)
      {
        PR[k * nnz + j] = VR[I[j]];
        PI[k * nnz + j] = VI[I[j]];
      }
    }
  }

  dots.SetSize(4 * K);
  coeff.SetSize(2 * K);
  dots.UseDevice(true);
  coeff.UseDevice(true);
}

template <bool Conj>
void LowRankComplexUpdate<Conj>::AddMult(const ComplexVector &x, ComplexVector &y,
                                         const std::complex<double> a) const
{
  // K is the same on every rank, so this early exit never skips the collective below.
  if (K == 0)
  {
    return;
  }
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const int nnz = idx.Size();
  const auto *I = idx.Read(use_dev);
  const auto *VR = vr.Read(use_dev);
  const auto *VI = vi.Read(use_dev);
  {
    // Raw local products vr·xr, vi·xi, vr·xi, vi·xr per term, summed in a fixed order.
    const auto *XR = x.Real().Read(use_dev);
    const auto *XI = x.Imag().Read(use_dev);
    auto *D = dots.Write(use_dev);
    mfem::forall_switch(use_dev, K,
                        [=] MFEM_HOST_DEVICE(int k)
                        {
                          double rr = 0.0, ii = 0.0, ri = 0.0, ir = 0.0;
                          for (int j = 0; j < nnz; j++)
                          {
                            const double vr_j = VR[k * nnz + j], vi_j = VI[k * nnz + j];
                            const double xr_j = XR[I[j]], xi_j = XI[I[j]];
                            rr += vr_j * xr_j;
                            ii += vi_j * xi_j;
                            ri += vr_j * xi_j;
                            ir += vi_j * xr_j;
                          }
                          D[4 * k + 0] = rr;
                          D[4 * k + 1] = ii;
                          D[4 * k + 2] = ri;
                          D[4 * k + 3] = ir;
                        });
  }

  // One reduction for all terms, then s_k = a g_k (v_kᴴ x) or a g_k (v_kᵀ x).
  auto *D = dots.HostReadWrite();
  Mpi::GlobalSum(4 * K, D, comm);
  {
    auto *C = coeff.HostWrite();
    for (int k = 0; k < K; k++)
    {
      const double rr = D[4 * k + 0], ii = D[4 * k + 1], ri = D[4 * k + 2],
                   ir = D[4 * k + 3];
      const auto dot = Conj ? std::complex<double>(rr + ii, ri - ir)
                            : std::complex<double>(rr - ii, ri + ir);
      const auto s = a * g[k] * dot;
      C[2 * k + 0] = s.real();
      C[2 * k + 1] = s.imag();
    }
  }

  // y += Σ_k s_k v_k, with each entry accumulating the terms in order.
  const auto *C = coeff.Read(use_dev);
  auto *YR = y.Real().ReadWrite(use_dev);
  auto *YI = y.Imag().ReadWrite(use_dev);
  const int num_terms = K;
  mfem::forall_switch(use_dev, nnz,
                      [=] MFEM_HOST_DEVICE(int j)
                      {
                        double yr = YR[I[j]], yi = YI[I[j]];
                        for (int k = 0; k < num_terms; k++)
                        {
                          const double sr = C[2 * k + 0], si = C[2 * k + 1];
                          const double vr_j = VR[k * nnz + j], vi_j = VI[k * nnz + j];
                          yr += sr * vr_j;
                          yr += -si * vi_j;
                          yi += si * vr_j;
                          yi += sr * vi_j;
                        }
                        YR[I[j]] = yr;
                        YI[I[j]] = yi;
                      });
}

template class LowRankComplexUpdate<true>;
template class LowRankComplexUpdate<false>;

}  // namespace palace
