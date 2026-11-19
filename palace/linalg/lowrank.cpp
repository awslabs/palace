// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowrank.hpp"

#include <algorithm>
#include <climits>
#include <map>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

// Threads per block of the device dot products, and the cap on the number of blocks of
// threads of a block of terms (which bounds the partials copied to the host) unless it has
// more terms.
constexpr int block_size = 256;
constexpr int max_blocks = 1024;

// re + i im += conj(v) x (Conj = true) or v x (Conj = false).
template <bool Conj>
MFEM_HOST_DEVICE inline void DotAdd(double vr, double vi, double xr, double xi, double &re,
                                    double &im)
{
  if constexpr (Conj)
  {
    re += vr * xr;
    re += vi * xi;
    im += vr * xi;
    im -= vi * xr;
  }
  else
  {
    re += vr * xr;
    re -= vi * xi;
    im += vr * xi;
    im += vi * xr;
  }
}

// y += s v.
MFEM_HOST_DEVICE inline void UpdateAdd(double sr, double si, double vr, double vi,
                                       double &yr, double &yi)
{
  yr += sr * vr;
  yr -= si * vi;
  yi += si * vr;
  yi += sr * vi;
}

}  // namespace

template <bool Conj>
LowRankComplexUpdate<Conj>::LowRankComplexUpdate(MPI_Comm comm, int n,
                                                 std::span<const ComplexVector *const> v,
                                                 std::span<const std::complex<double>> g,
                                                 std::span<const int> block)
  : comm(comm), n(n), K(static_cast<int>(v.size())), g(g.begin(), g.end())
{
  MFEM_VERIFY(v.size() == g.size() && v.size() == block.size(),
              "Mismatched numbers of low-rank vectors, weights and block ids!");

  // The terms of each block in order, and the local entries where any of them is nonzero.
  // Blocks without local entries are dropped: their terms contribute zero here.
  std::map<int, std::vector<int>> block_terms;
  for (int k = 0; k < K; k++)
  {
    MFEM_VERIFY(v[k]->Size() == n, "Mismatched low-rank vector sizes!");
    block_terms[block[k]].push_back(k);
  }
  std::vector<std::vector<int>> terms, ents;
  std::vector<char> on(n);
  for (auto &[id, bt] : block_terms)
  {
    std::fill(on.begin(), on.end(), 0);
    for (int k : bt)
    {
      const auto *VR = v[k]->Real().HostRead();
      const auto *VI = v[k]->Imag().HostRead();
      for (int i = 0; i < n; i++)
      {
        on[i] = on[i] || VR[i] != 0.0 || VI[i] != 0.0;
      }
    }
    std::vector<int> e;
    for (int i = 0; i < n; i++)
    {
      if (on[i])
      {
        e.push_back(i);
      }
    }
    if (!e.empty())
    {
      MFEM_VERIFY(static_cast<long long>(bt.size()) * static_cast<long long>(e.size()) <=
                      INT_MAX,
                  "Too many packed low-rank values in a block!");
      terms.push_back(std::move(bt));
      ents.push_back(std::move(e));
    }
  }

  // Pack the values of the terms of each block on its entries, term-major. The device dot
  // products use enough blocks of threads per term for the entries of its block, capped so
  // that the partials stay small. They depend only on the sizes of the block, so the
  // summation order is the same on every device.
  blocks.resize(terms.size());
  int num_partials = 0, num_local = 0;
  for (std::size_t b = 0; b < blocks.size(); b++)
  {
    auto &blk = blocks[b];
    const int nt = static_cast<int>(terms[b].size()),
              nnz = static_cast<int>(ents[b].size());
    blk.terms = std::move(terms[b]);
    blk.idx.SetSize(nnz);
    std::copy(ents[b].begin(), ents[b].end(), blk.idx.HostWrite());
    blk.vr.SetSize(nt * nnz);
    blk.vi.SetSize(nt * nnz);
    blk.vr.UseDevice(true);
    blk.vi.UseDevice(true);
    auto *PR = blk.vr.HostWrite();
    auto *PI = blk.vi.HostWrite();
    for (int t = 0; t < nt; t++)
    {
      const auto *VR = v[blk.terms[t]]->Real().HostRead();
      const auto *VI = v[blk.terms[t]]->Imag().HostRead();
      for (int j = 0; j < nnz; j++)
      {
        PR[t * nnz + j] = VR[ents[b][j]];
        PI[t * nnz + j] = VI[ents[b][j]];
      }
    }
    blk.nb = std::min((nnz + block_size - 1) / block_size, std::max(1, max_blocks / nt));
    blk.offset = num_partials;
    blk.first = num_local;
    num_partials += nt * blk.nb;
    num_local += nt;
  }
  partial.SetSize(2 * num_partials);
  partial.UseDevice(true);
  dots.SetSize(2 * K);
  coeff.SetSize(2 * num_local);
  coeff.UseDevice(true);
}

template <bool Conj>
int LowRankComplexUpdate<Conj>::NumPackedEntries() const
{
  int nnz = 0;
  for (const auto &blk : blocks)
  {
    nnz += blk.idx.Size();
  }
  return nnz;
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
  MFEM_ASSERT(x.Size() == n && y.Size() == n, "Mismatched low-rank vector sizes!");
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const bool on_dev = use_dev && mfem::Device::Allows(mfem::Backend::DEVICE_MASK);

  // Local dot products v_kᴴ x or v_kᵀ x of the terms of the local blocks, summed in a fixed
  // order. The other terms contribute zero.
  auto *D = dots.HostWrite();
  std::fill(D, D + 2 * K, 0.0);
  if (on_dev)
  {
    // Each block of threads reduces a strided slice of the entries of one term in shared
    // memory, and the partials of each term are summed on the host in order.
    const auto *XR = x.Real().Read(), *XI = x.Imag().Read();
    auto *P = partial.Write();
    for (const auto &blk : blocks)
    {
      const int nnz = blk.idx.Size(), NB = blk.nb;
      const auto *I = blk.idx.Read();
      const auto *VR = blk.vr.Read(), *VI = blk.vi.Read();
      auto *PB = P + 2 * blk.offset;
      mfem::forall_2D(static_cast<int>(blk.terms.size()) * NB, block_size, 1,
                      [=] MFEM_HOST_DEVICE(int b)
                      {
                        const int v0 = (b / NB) * nnz;
                        MFEM_SHARED double sm[2][block_size];
                        MFEM_FOREACH_THREAD(tid, x, block_size)
                        {
                          double re = 0.0, im = 0.0;
                          for (int j = (b % NB) * block_size + tid; j < nnz;
                               j += NB * block_size)
                          {
                            const int i = I[j];
                            DotAdd<Conj>(VR[v0 + j], VI[v0 + j], XR[i], XI[i], re, im);
                          }
                          sm[0][tid] = re;
                          sm[1][tid] = im;
                        }
                        for (int w = block_size / 2; w > 0; w /= 2)
                        {
                          MFEM_SYNC_THREAD;
                          MFEM_FOREACH_THREAD(tid, x, w)
                          {
                            sm[0][tid] += sm[0][tid + w];
                            sm[1][tid] += sm[1][tid + w];
                          }
                        }
                        MFEM_FOREACH_THREAD(tid, x, 1)
                        {
                          PB[2 * b + 0] = sm[0][0];
                          PB[2 * b + 1] = sm[1][0];
                        }
                      });
    }
    const auto *PH = partial.HostRead();
    for (const auto &blk : blocks)
    {
      for (std::size_t t = 0; t < blk.terms.size(); t++)
      {
        for (int b = 0; b < blk.nb; b++)
        {
          const int p = blk.offset + static_cast<int>(t) * blk.nb + b;
          D[2 * blk.terms[t] + 0] += PH[2 * p + 0];
          D[2 * blk.terms[t] + 1] += PH[2 * p + 1];
        }
      }
    }
  }
  else
  {
    // On the host, each term sums the entries of its block in four interleaved lanes (entry
    // j in lane j % 4), each with separate sums of the four real products, combined in a
    // fixed order at the end. The loop then carries no dependency from one entry to the
    // next and needs no reduction across vector lanes, which GCC compiles into strictly
    // ordered (slow) instructions for SVE.
    const auto *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
    for (const auto &blk : blocks)
    {
      const int nnz = blk.idx.Size();
      const int *T = blk.terms.data();
      const auto *I = blk.idx.HostRead();
      const auto *VR = blk.vr.HostRead(), *VI = blk.vi.HostRead();
      mfem::forall_switch(use_dev, static_cast<int>(blk.terms.size()),
                          [=] MFEM_HOST_DEVICE(int t)
                          {
                            constexpr int L = 4;
                            const int v0 = t * nnz;
                            double rr[L] = {}, ii[L] = {}, ri[L] = {}, ir[L] = {};
                            int j = 0;
                            for (; j + L <= nnz; j += L)
                            {
                              for (int l = 0; l < L; l++)
                              {
                                const int i = I[j + l];
                                const double vr_j = VR[v0 + j + l], vi_j = VI[v0 + j + l];
                                rr[l] += vr_j * XR[i];
                                ii[l] += vi_j * XI[i];
                                ri[l] += vr_j * XI[i];
                                ir[l] += vi_j * XR[i];
                              }
                            }
                            for (int l = 0; j < nnz; j++, l++)
                            {
                              const int i = I[j];
                              const double vr_j = VR[v0 + j], vi_j = VI[v0 + j];
                              rr[l] += vr_j * XR[i];
                              ii[l] += vi_j * XI[i];
                              ri[l] += vr_j * XI[i];
                              ir[l] += vi_j * XR[i];
                            }
                            for (int l = 1; l < L; l++)
                            {
                              rr[0] += rr[l];
                              ii[0] += ii[l];
                              ri[0] += ri[l];
                              ir[0] += ir[l];
                            }
                            D[2 * T[t] + 0] = Conj ? rr[0] + ii[0] : rr[0] - ii[0];
                            D[2 * T[t] + 1] = Conj ? ri[0] - ir[0] : ri[0] + ir[0];
                          });
    }
  }

  // One reduction for all terms, then the coefficients s_k = a g_k (v_kᴴ x) or
  // a g_k (v_kᵀ x) of the terms of the local blocks.
  Mpi::GlobalSum(2 * K, D, comm);
  auto *CH = coeff.HostWrite();
  for (const auto &blk : blocks)
  {
    for (std::size_t t = 0; t < blk.terms.size(); t++)
    {
      const int k = blk.terms[t];
      const auto s = a * g[k] * std::complex<double>(D[2 * k], D[2 * k + 1]);
      CH[2 * (blk.first + t) + 0] = s.real();
      CH[2 * (blk.first + t) + 1] = s.imag();
    }
  }

  // y += Σ_k s_k v_k, one block after the other (so blocks may share entries), with each
  // entry accumulating the terms of its block in order.
  const auto *DC = coeff.Read(use_dev);
  auto *YR = y.Real().ReadWrite(use_dev), *YI = y.Imag().ReadWrite(use_dev);
  for (const auto &blk : blocks)
  {
    const int nt = static_cast<int>(blk.terms.size()), nnz = blk.idx.Size();
    const auto *I = blk.idx.Read(use_dev);
    const auto *VR = blk.vr.Read(use_dev), *VI = blk.vi.Read(use_dev);
    const auto *C = DC + 2 * blk.first;
    mfem::forall_switch(use_dev, nnz,
                        [=] MFEM_HOST_DEVICE(int j)
                        {
                          const int i = I[j];
                          double yr = YR[i], yi = YI[i];
                          for (int t = 0; t < nt; t++)
                          {
                            UpdateAdd(C[2 * t], C[2 * t + 1], VR[t * nnz + j],
                                      VI[t * nnz + j], yr, yi);
                          }
                          YR[i] = yr;
                          YI[i] = yi;
                        });
  }
}

template class LowRankComplexUpdate<true>;
template class LowRankComplexUpdate<false>;

}  // namespace palace
