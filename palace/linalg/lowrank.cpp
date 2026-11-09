// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lowrank.hpp"

#include <algorithm>
#include <climits>
#include <numeric>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

// Threads per block of the device dot products, and the cap on the total number of blocks
// (which bounds the block partials copied to the host) unless there are more local terms.
constexpr int block_size = 256;
constexpr int max_blocks = 1024;

int FindRoot(std::vector<int> &parent, int k)
{
  while (parent[k] != k)
  {
    parent[k] = parent[parent[k]];
    k = parent[k];
  }
  return k;
}

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
                                                 std::span<const std::complex<double>> g)
  : comm(comm), n(n), K(static_cast<int>(v.size())), g(g.begin(), g.end())
{
  MFEM_VERIFY(v.size() == g.size(), "Mismatched number of low-rank vectors and weights!");

  // Group the terms by connected support: two terms are in the same group if both are
  // nonzero at some local entry, so the groups have disjoint entries (typically one group
  // per port). Each group drops the entries where all of its terms are exactly zero, which
  // is exact. The root of a group is its first term.
  std::vector<int> owner(n, -1), parent(K);
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<bool> supported(K, false);
  for (int k = 0; k < K; k++)
  {
    MFEM_VERIFY(v[k]->Size() == n, "Mismatched low-rank vector sizes!");
    const auto *VR = v[k]->Real().HostRead();
    const auto *VI = v[k]->Imag().HostRead();
    for (int i = 0; i < n; i++)
    {
      if (VR[i] != 0.0 || VI[i] != 0.0)
      {
        supported[k] = true;
        if (owner[i] < 0)
        {
          owner[i] = k;
        }
        else
        {
          const int r1 = FindRoot(parent, owner[i]), r2 = FindRoot(parent, k);
          parent[std::max(r1, r2)] = std::min(r1, r2);
        }
      }
    }
  }
  std::vector<int> group(K, -1), terms, ents;
  int num_groups = 0;
  for (int k = 0; k < K; k++)
  {
    if (supported[k])
    {
      const int r = FindRoot(parent, k);
      if (group[r] < 0)
      {
        group[r] = num_groups++;
      }
      group[k] = group[r];
      terms.push_back(k);
    }
  }
  for (int i = 0; i < n; i++)
  {
    if (owner[i] >= 0)
    {
      ents.push_back(i);
    }
  }

  // Order the local terms and the packed entries by group, keeping the order of the terms
  // and of the local indices within each group.
  std::stable_sort(terms.begin(), terms.end(),
                   [&](int k1, int k2) { return group[k1] < group[k2]; });
  std::stable_sort(ents.begin(), ents.end(),
                   [&](int i1, int i2) { return group[owner[i1]] < group[owner[i2]]; });
  const int num_local = static_cast<int>(terms.size());
  const int num_entries = static_cast<int>(ents.size());
  perm.SetSize(num_local);
  term_group.SetSize(num_local);
  idx.SetSize(num_entries);
  ent_group.SetSize(num_entries);
  group_terms.SetSize(num_groups + 1);
  group_ents.SetSize(num_groups + 1);
  group_vals.SetSize(num_groups + 1);
  auto *PM = perm.HostWrite();
  auto *TG = term_group.HostWrite();
  auto *I = idx.HostWrite();
  auto *EG = ent_group.HostWrite();
  auto *GT = group_terms.HostWrite();
  auto *GE = group_ents.HostWrite();
  auto *GV = group_vals.HostWrite();
  std::fill(GT, GT + num_groups + 1, 0);
  std::fill(GE, GE + num_groups + 1, 0);
  for (int t = 0; t < num_local; t++)
  {
    PM[t] = terms[t];
    TG[t] = group[terms[t]];
    GT[TG[t] + 1]++;
  }
  for (int j = 0; j < num_entries; j++)
  {
    I[j] = ents[j];
    EG[j] = group[owner[ents[j]]];
    GE[EG[j] + 1]++;
  }
  int max_nnz = 0;
  long long num_vals = 0;
  for (int gi = 0; gi < num_groups; gi++)
  {
    num_vals += static_cast<long long>(GT[gi + 1]) * GE[gi + 1];
  }
  // Collective, so that every process fails rather than only the ones over the limit.
  bool too_many = num_vals > INT_MAX;
  Mpi::GlobalOr(1, &too_many, comm);
  MFEM_VERIFY(!too_many, "Too many packed low-rank values on a process!");
  GV[0] = 0;
  for (int gi = 0; gi < num_groups; gi++)
  {
    max_nnz = std::max(max_nnz, GE[gi + 1]);
    GV[gi + 1] = GV[gi] + GT[gi + 1] * GE[gi + 1];
    GT[gi + 1] += GT[gi];
    GE[gi + 1] += GE[gi];
  }

  // Pack the values of the terms of each group on its entries, term-major.
  vr.SetSize(GV[num_groups]);
  vi.SetSize(GV[num_groups]);
  vr.UseDevice(true);
  vi.UseDevice(true);
  auto *PR = vr.HostWrite();
  auto *PI = vi.HostWrite();
  for (int t = 0; t < num_local; t++)
  {
    const int gi = TG[t], e0 = GE[gi], nnz = GE[gi + 1] - e0;
    const int v0 = GV[gi] + (t - GT[gi]) * nnz;
    const auto *VR = v[PM[t]]->Real().HostRead();
    const auto *VI = v[PM[t]]->Imag().HostRead();
    for (int j = 0; j < nnz; j++)
    {
      PR[v0 + j] = VR[I[e0 + j]];
      PI[v0 + j] = VI[I[e0 + j]];
    }
  }

  // Blocks per term of the device dot products: enough for the largest group, capped so
  // that the block partials stay small. They depend only on the sizes, so the summation
  // order is the same on every device.
  nb = (num_local == 0) ? 0
                        : std::min((max_nnz + block_size - 1) / block_size,
                                   std::max(1, max_blocks / num_local));
  partial.SetSize(2 * num_local * nb);
  partial.UseDevice(true);
  dots.SetSize(2 * K);
  coeff.SetSize(2 * num_local);
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
  MFEM_ASSERT(x.Size() == n && y.Size() == n, "Mismatched low-rank vector sizes!");
  const bool use_dev = x.UseDevice() || y.UseDevice();
  const bool on_dev = use_dev && mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  const int num_local = perm.Size(), num_entries = idx.Size();
  const auto *PM = perm.HostRead();

  // Local dot products v_kᴴ x or v_kᵀ x of the local terms, summed in a fixed order. The
  // terms without local support contribute zero.
  auto *D = dots.HostWrite();
  std::fill(D, D + 2 * K, 0.0);
  if (on_dev)
  {
    // Each block reduces a strided slice of the entries of one term in shared memory, and
    // the block partials of each term are summed on the host in order.
    const int NB = nb;
    const auto *TG = term_group.Read();
    const auto *GT = group_terms.Read(), *GE = group_ents.Read(), *GV = group_vals.Read();
    const auto *I = idx.Read();
    const auto *VR = vr.Read(), *VI = vi.Read();
    const auto *XR = x.Real().Read(true), *XI = x.Imag().Read(true);
    auto *P = partial.Write();
    mfem::forall_2D(num_local * NB, block_size, 1,
                    [=] MFEM_HOST_DEVICE(int blk)
                    {
                      const int t = blk / NB, b = blk % NB, gi = TG[t];
                      const int e0 = GE[gi], nnz = GE[gi + 1] - e0;
                      const int v0 = GV[gi] + (t - GT[gi]) * nnz;
                      MFEM_SHARED double sm[2][block_size];
                      MFEM_FOREACH_THREAD(tid, x, block_size)
                      {
                        double re = 0.0, im = 0.0;
                        for (int j = b * block_size + tid; j < nnz; j += NB * block_size)
                        {
                          const int i = I[e0 + j];
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
                        P[2 * blk + 0] = sm[0][0];
                        P[2 * blk + 1] = sm[1][0];
                      }
                    });
    const auto *PH = partial.HostRead();
    for (int t = 0; t < num_local; t++)
    {
      for (int b = 0; b < NB; b++)
      {
        D[2 * PM[t] + 0] += PH[2 * (t * NB + b) + 0];
        D[2 * PM[t] + 1] += PH[2 * (t * NB + b) + 1];
      }
    }
  }
  else
  {
    // On the host, each term accumulates the entries of its group in order.
    const auto *TG = term_group.HostRead();
    const auto *GT = group_terms.HostRead(), *GE = group_ents.HostRead();
    const auto *GV = group_vals.HostRead();
    const auto *I = idx.HostRead();
    const auto *VR = vr.HostRead(), *VI = vi.HostRead();
    const auto *XR = x.Real().HostRead(), *XI = x.Imag().HostRead();
    mfem::forall_switch(use_dev, num_local,
                        [=] MFEM_HOST_DEVICE(int t)
                        {
                          const int gi = TG[t], e0 = GE[gi], nnz = GE[gi + 1] - e0;
                          const int v0 = GV[gi] + (t - GT[gi]) * nnz;
                          double re = 0.0, im = 0.0;
                          for (int j = 0; j < nnz; j++)
                          {
                            const int i = I[e0 + j];
                            DotAdd<Conj>(VR[v0 + j], VI[v0 + j], XR[i], XI[i], re, im);
                          }
                          D[2 * PM[t] + 0] = re;
                          D[2 * PM[t] + 1] = im;
                        });
  }

  // One reduction for all terms, then the coefficients s_k = a g_k (v_kᴴ x) or
  // a g_k (v_kᵀ x) of the local terms.
  Mpi::GlobalSum(2 * K, D, comm);
  auto *C = coeff.HostWrite();
  for (int t = 0; t < num_local; t++)
  {
    const int k = PM[t];
    const auto s = a * g[k] * std::complex<double>(D[2 * k], D[2 * k + 1]);
    C[2 * t + 0] = s.real();
    C[2 * t + 1] = s.imag();
  }

  // y += Σ_k s_k v_k, with each entry accumulating the terms of its group in order.
  const auto *EG = ent_group.Read(use_dev);
  const auto *GT = group_terms.Read(use_dev), *GE = group_ents.Read(use_dev);
  const auto *GV = group_vals.Read(use_dev);
  const auto *I = idx.Read(use_dev);
  const auto *VR = vr.Read(use_dev), *VI = vi.Read(use_dev);
  const auto *DC = coeff.Read(use_dev);
  auto *YR = y.Real().ReadWrite(use_dev), *YI = y.Imag().ReadWrite(use_dev);
  mfem::forall_switch(use_dev, num_entries,
                      [=] MFEM_HOST_DEVICE(int e)
                      {
                        const int gi = EG[e], e0 = GE[gi], nnz = GE[gi + 1] - e0;
                        const int t0 = GT[gi], t1 = GT[gi + 1], v0 = GV[gi] + e - e0;
                        const int i = I[e];
                        double yr = YR[i], yi = YI[i];
                        for (int t = t0; t < t1; t++)
                        {
                          const int v = v0 + (t - t0) * nnz;
                          UpdateAdd(DC[2 * t], DC[2 * t + 1], VR[v], VI[v], yr, yi);
                        }
                        YR[i] = yr;
                        YI[i] = yi;
                      });
}

template class LowRankComplexUpdate<true>;
template class LowRankComplexUpdate<false>;

}  // namespace palace
