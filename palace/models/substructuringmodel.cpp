// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <cmath>
#include "utils/timer.hpp"

namespace palace
{

void SubstructuringSolver::Impl::WriteModel(
    const std::string &path, bool append,
    const std::function<void(std::ofstream &)> &write) const
{
  int ok = 1;
  if (Mpi::Root(parent_fes.GetComm()))
  {
    std::ofstream f(path, std::ios::binary | (append ? std::ios::app : std::ios::trunc));
    if (f)
    {
      write(f);
      f.close();
    }
    ok = !f.fail();
  }
  Mpi::Broadcast(1, &ok, 0, parent_fes.GetComm());
  MFEM_VERIFY(ok, "Cannot write the substructuring model \"" << path << "\"!");
}

void SubstructuringSolver::Impl::AppendSourceModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(src_ids.size());
  const std::vector<double> G_full = GatherRows(SG_rows, K);
  const std::vector<double> H_full = GatherRows(SH_rows, K);
  WriteModel(
      path, true,
      [&](std::ofstream &f)
      {
        f.write(reinterpret_cast<const char *>(&kSourceModesMagic), sizeof(int));
        f.write(reinterpret_cast<const char *>(&K), sizeof(int));
        f.write(reinterpret_cast<const char *>(src_ids.data()), sizeof(int) * K);
        f.write(reinterpret_cast<const char *>(src_fp.data()), sizeof(double) * 2 * K);
        f.write(reinterpret_cast<const char *>(G_full.data()),
                sizeof(double) * G_full.size());
        f.write(reinterpret_cast<const char *>(H_full.data()),
                sizeof(double) * H_full.size());
        f.write(reinterpret_cast<const char *>(SC.data()), sizeof(double) * SC.size());
        f.write(reinterpret_cast<const char *>(SB.data()), sizeof(double) * SB.size());
      });
  src_saved = true;
}

int SubstructuringSolver::Impl::SignatureFileWidth(int sig_type)
{
  return (sig_type == 3) ? 13 : (sig_type == 2) ? 12 : (sig_type == 1) ? 3 : 0;
}

SubstructuringSolver::Impl::SavedHead
SubstructuringSolver::Impl::ReadSavedHead(const std::string &path, MPI_Comm comm)
{
  SavedHead h;
  // Version, order, nG, signature type, width, with the fingerprint, with the moments.
  int head[7] = {-1, 0, 0, 0, 0, 0, 0};
  auto offsets = [&h](std::streamoff nG, std::streamoff width)
  {
    h.S_pos = kHeadInts * static_cast<std::streamoff>(sizeof(int)) +
              static_cast<std::streamoff>(sizeof(double)) * nG * width;
    h.sections_pos = h.S_pos + static_cast<std::streamoff>(sizeof(double)) * nG * nG;
  };
  if (Mpi::Root(comm))
  {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    const std::streamoff size = f ? static_cast<std::streamoff>(f.tellg()) : -1;
    f.seekg(0);
    int v[kHeadInts] = {0, 0, 0, 0, 0};
    head[0] = (size < 0) ? -1 : 0;
    if (f.read(reinterpret_cast<char *>(v), sizeof(v)) && v[0] == kModelMagic)
    {
      head[0] = v[1];
      head[1] = v[2];
      head[2] = v[3];
      head[3] = v[4];
      // The file holds the signatures and S_E (in floating point: any nG read).
      const double n = v[3], w = SignatureFileWidth(v[4]);
      const double end = kHeadInts * sizeof(int) + sizeof(double) * n * (w + n);
      head[4] =
          (v[1] == kModelVersion && n >= 0.0 && w > 0.0 && end <= static_cast<double>(size))
              ? static_cast<int>(w)
              : 0;
    }
    if (head[4] > 0)
    {
      offsets(head[2], head[4]);
      h.sig.resize(static_cast<std::size_t>(head[2]) * head[4]);
      f.read(reinterpret_cast<char *>(h.sig.data()), sizeof(double) * h.sig.size());
      f.seekg(h.sections_pos);
      int magic = 0;
      auto read = [&f](double *x, int n)
      {
        return static_cast<bool>(f.read(reinterpret_cast<char *>(x), sizeof(double) * n));
      };
      while (f.read(reinterpret_cast<char *>(&magic), sizeof(int)))
      {
        if (magic == kLengthMagic)
        {
          read(&h.length, 1);
        }
        else if (magic == kEpsMagic)
        {
          read(&h.regularization, 1);
        }
        else if (magic == kEnvFpMagic)
        {
          h.env_fp.resize(5);
          head[5] = read(h.env_fp.data(), 5);
        }
        else if (magic == kEnvMeshMagic)
        {
          head[6] = read(h.env_moments.data(), 6) && read(h.region_moments.data(), 6);
        }
        else
        {
          break;  // the sections loaded with S_E
        }
      }
    }
  }
  Mpi::Broadcast(7, head, 0, comm);
  h.version = head[0];
  h.order = head[1];
  h.nG = head[2];
  h.sig_type = head[3];
  h.width = (head[5] && head[6]) ? head[4] : 0;
  if (h.width == 0)
  {
    return h;
  }
  Mpi::Broadcast(1, &h.length, 0, comm);
  Mpi::Broadcast(1, &h.regularization, 0, comm);
  offsets(h.nG, h.width);
  h.sig.resize(static_cast<std::size_t>(h.nG) * h.width);
  Mpi::Broadcast(static_cast<int>(h.sig.size()), h.sig.data(), 0, comm);
  h.env_fp.resize(5);
  Mpi::Broadcast(5, h.env_fp.data(), 0, comm);
  Mpi::Broadcast(6, h.env_moments.data(), 0, comm);
  Mpi::Broadcast(6, h.region_moments.data(), 0, comm);
  return h;
}

mfem::Vector SubstructuringSolver::Impl::FingerprintField(int k) const
{
  mfem::Vector r = palace::FingerprintField(parent_fes, k);
  if (!broken)
  {
    return r;
  }
  mfem::Vector rb(nt);
  all_versions->MultTranspose(r, rb);
  return rb;
}

std::array<double, 5> SubstructuringSolver::Impl::EnvironmentFingerprint() const
{
  // The environment does not change within a run: computed once.
  if (env_fp)
  {
    return *env_fp;
  }
  double loc[2] = {0.0, static_cast<double>(env_dbc_local)};
  for (int i = 0; i < nt; i++)
  {
    loc[0] += (is_env_int[i] || is_gamma[i]) ? 1.0 : 0.0;
  }
  const mfem::Vector r2 = FingerprintField(2);
  const std::vector<const mfem::Vector *> r = {&CachedFingerprintField(0),
                                               &CachedFingerprintField(1), &r2};
  const std::vector<mfem::Vector> Ar = EnvApply(r);
  double glob[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
  for (int k = 0; k < 3; k++)
  {
    glob[2 + k] = mfem::InnerProduct(parent_fes.GetComm(), *r[k], Ar[k]);
  }
  env_fp = {glob[0], glob[1], glob[2], glob[3], glob[4]};
  return *env_fp;
}

const mfem::Vector &SubstructuringSolver::Impl::CachedFingerprintField(int k) const
{
  if (fp_fields[k].Size() != nt)
  {
    fp_fields[k] = FingerprintField(k);
  }
  return fp_fields[k];
}

void SubstructuringSolver::Impl::AppendEnvMoments(const std::string &path,
                                                  const std::vector<char> &written)
{
  BlockTimer bt(Timer::IO);
  const auto &env = iodata.solver.substructuring->environment_attributes;
  const std::array<double, 6> m[2] = {ElementMoments(parent, env),
                                      ElementMoments(parent, env, written)};  // collective
  WriteModel(path, true,
             [&](std::ofstream &f)
             {
               f.write(reinterpret_cast<const char *>(&kEnvMeshMagic), sizeof(int));
               f.write(reinterpret_cast<const char *>(m[0].data()), sizeof(double) * 6);
               f.write(reinterpret_cast<const char *>(m[1].data()), sizeof(double) * 6);
             });
}

void SubstructuringSolver::Impl::AppendEnvFingerprint(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const std::array<double, 5> fp = EnvironmentFingerprint();  // collective
  WriteModel(path, true,
             [&](std::ofstream &f)
             {
               f.write(reinterpret_cast<const char *>(&kEnvFpMagic), sizeof(int));
               f.write(reinterpret_cast<const char *>(fp.data()),
                       sizeof(double) * fp.size());
               if (magnetostatic)
               {
                 f.write(reinterpret_cast<const char *>(&kEpsMagic), sizeof(int));
                 f.write(reinterpret_cast<const char *>(&mag_eps), sizeof(double));
               }
             });
}

std::vector<char>
SubstructuringSolver::Impl::WrittenEnvTrueDofs(const std::vector<char> &keep) const
{
  const auto &P = *parent_fes.Dof_TrueDof_Matrix();
  mfem::Vector t(P.Width()), l(P.Height()), lw(P.Height());
  for (int i = 0; i < t.Size(); i++)
  {
    t(i) = keep[i] ? 1.0 : 0.0;
  }
  P.Mult(t, l);
  lw = 0.0;
  mfem::Array<int> dofs;
  auto index = [](int d) { return (d >= 0) ? d : -1 - d; };
  for (int e = 0; e < parent.GetNE(); e++)
  {
    if (ea_arr.Find(parent.GetAttribute(e)) < 0)
    {
      continue;
    }
    parent_fes.GetElementDofs(e, dofs);
    if (std::any_of(dofs.begin(), dofs.end(), [&](int d) { return l(index(d)) != 0.0; }))
    {
      for (int d : dofs)
      {
        lw(index(d)) = 1.0;
      }
    }
  }
  P.AbsMultTranspose(1.0, lw, 0.0, t);
  std::vector<char> written(t.Size());
  for (int i = 0; i < t.Size(); i++)
  {
    written[i] = (t(i) > 0.0);
  }
  return written;
}

std::vector<char>
SubstructuringSolver::Impl::RegionFilmClosure(std::vector<char> &keep,
                                              const std::set<int> &whole) const
{
  // The film faces of the environment (of the films not given whole), and the DOFs of the
  // film faces of the region.
  const int nt0 = parent_fes.GetTrueVSize();
  const auto &P = *parent_fes.Dof_TrueDof_Matrix();
  auto index = [](int d) { return (d >= 0) ? d : -1 - d; };
  std::set<int> films;
  for (const auto &sc : iodata.boundaries.superconductor)
  {
    if (sc.two_sided)
    {
      films.insert(sc.attributes.begin(), sc.attributes.end());
    }
  }
  std::vector<int> faces;
  mfem::Vector lr(P.Height()), tr(nt0);
  lr = 0.0;
  mfem::Array<int> dofs;
  for (int be = 0; be < parent.GetNBE(); be++)
  {
    const int attr = parent.GetBdrAttribute(be);
    if (!films.contains(attr))
    {
      continue;
    }
    int e, info;
    parent.GetBdrElementAdjacentElement(be, e, info);
    if (ea_arr.Find(parent.GetAttribute(e)) >= 0)
    {
      if (!whole.contains(attr))
      {
        faces.push_back(be);
      }
      continue;
    }
    parent_fes.GetBdrElementDofs(be, dofs);
    for (int d : dofs)
    {
      lr(index(d)) = 1.0;
    }
  }
  P.AbsMultTranspose(1.0, lr, 0.0, tr);
  // The interface DOFs on a film crossing Gamma, and those at the edge of a film of the
  // region (kept).
  std::vector<char> split(nt0, 0), gamma(nt0, 0);
  for (int t = 0; t < nt; t++)
  {
    gamma[ParentTrueDof(t)] = gamma[ParentTrueDof(t)] || is_gamma[t];
  }
  const double *r = tr.HostRead();
  for (int t0 = 0; t0 < nt0; t0++)
  {
    split[t0] = gamma[t0] && broken->GetNumVersions(t0) > 1;
    keep[t0] = keep[t0] || (gamma[t0] && !split[t0] && r[t0] > 0.0);
  }
  const auto written = WrittenEnvTrueDofs(keep);
  int touched = 0;
  for (int t0 = 0; t0 < nt0 && !touched; t0++)
  {
    touched = split[t0] && written[t0];
  }
  Mpi::GlobalMax(1, &touched, parent.GetComm());
  std::vector<char> plain(parent.GetNBE(), 0);
  for (int be : faces)
  {
    plain[be] = 1;
  }
  if (!touched)
  {
    return plain;
  }
  mfem::Vector tk(nt0), lk(P.Height()), lm(P.Height()), tm(nt0);
  for (int t0 = 0; t0 < nt0; t0++)
  {
    tk(t0) = split[t0] ? 1.0 : 0.0;
  }
  P.Mult(tk, lk);
  const double *l = lk.HostRead();
  lm = 0.0;
  for (int be : faces)
  {
    parent_fes.GetBdrElementDofs(be, dofs);
    if (std::any_of(dofs.begin(), dofs.end(), [&](int d) { return l[index(d)] != 0.0; }))
    {
      plain[be] = 0;
      for (int d : dofs)
      {
        lm(index(d)) = 1.0;
      }
    }
  }
  P.AbsMultTranspose(1.0, lm, 0.0, tm);
  const double *m = tm.HostRead();
  for (int t0 = 0; t0 < nt0; t0++)
  {
    keep[t0] = keep[t0] || (m[t0] > 0.0);
  }
  return plain;
}

std::vector<char> SubstructuringSolver::Impl::WriteRegionMesh(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  std::vector<char> keep(parent_fes.GetTrueVSize(), 0);
  for (int t = 0; t < nt; t++)
  {
    keep[ParentTrueDof(t)] = keep[ParentTrueDof(t)] || region_dbc[t];
  }
  std::set<int> loop_films;
  if (magnetostatic)
  {
    std::vector<int> attrs;
    for (const auto &[idx, data] : iodata.boundaries.current)
    {
      for (const auto &elem : data.elements)
      {
        attrs.insert(attrs.end(), elem.attributes.begin(), elem.attributes.end());
        if (elem.aperture)
        {
          attrs.insert(attrs.end(), elem.aperture->attributes.begin(),
                       elem.aperture->attributes.end());
        }
      }
    }
    for (const auto &[idx, data] : iodata.boundaries.fluxloop)
    {
      attrs.insert(attrs.end(), data.film_attributes.begin(), data.film_attributes.end());
      attrs.insert(attrs.end(), data.hole_attributes.begin(), data.hole_attributes.end());
      loop_films.insert(data.film_attributes.begin(), data.film_attributes.end());
    }
    const auto mark = BoundaryTrueDofs(parent_fes, attrs);
    for (std::size_t i = 0; i < keep.size(); i++)
    {
      keep[i] = keep[i] || mark[i];
    }
  }
  const std::vector<char> plain =
      broken ? RegionFilmClosure(keep, loop_films) : std::vector<char>();
  const auto &sub = *iodata.solver.substructuring;
  std::vector<char> written = palace::WriteRegionMesh(
      path, parent_fes, sub.region_attributes, sub.environment_attributes, keep,
      iodata.units.GetMeshLengthRelativeScale(), plain);
  Mpi::Print(" Writing the mesh of the region for online runs to {}\n", path);
  return written;
}

void SubstructuringSolver::Impl::CheckEnvFingerprint()
{
  const std::array<double, 5> fp = EnvironmentFingerprint();  // collective
  MFEM_VERIFY(have_env_fp, "The saved substructuring model has no environment "
                           "fingerprint: rerun in \"Offline\" mode to condense it again!");
  bool ok = (fp[0] == saved_env_fp[0] && fp[1] == saved_env_fp[1]);
  for (int q = 2; q < 5; q++)
  {
    ok = ok && std::abs(fp[q] - saved_env_fp[q]) <=
                   kFingerprintTol * std::max(std::abs(fp[q]), std::abs(saved_env_fp[q]));
  }
  MFEM_VERIFY(ok, kEnvironmentDiffers);
}

void SubstructuringSolver::Impl::AppendSK(const std::string &path) const
{
  BlockTimer bt(Timer::IO);
  const std::vector<double> full = GatherRows(SK_rows, nG_global);
  WriteModel(path, true,
             [&](std::ofstream &f)
             {
               f.write(reinterpret_cast<const char *>(&kSKenMagic), sizeof(int));
               f.write(reinterpret_cast<const char *>(full.data()),
                       sizeof(double) * full.size());
             });
}

void SubstructuringSolver::Impl::AppendModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(mode_ids.size());
  const int has_gk = have_gk ? 1 : 0;
  const std::vector<double> G_full = GatherRows(G_rows, K);
  const std::vector<double> GK_full =
      has_gk ? GatherRows(GK_rows, K) : std::vector<double>();
  WriteModel(path, true,
             [&](std::ofstream &f)
             {
               f.write(reinterpret_cast<const char *>(&kModesMagic), sizeof(int));
               f.write(reinterpret_cast<const char *>(&K), sizeof(int));
               f.write(reinterpret_cast<const char *>(mode_ids.data()), sizeof(int) * K);
               f.write(reinterpret_cast<const char *>(mode_fp.data()),
                       sizeof(double) * 2 * K);
               f.write(reinterpret_cast<const char *>(&has_gk), sizeof(int));
               f.write(reinterpret_cast<const char *>(G_full.data()),
                       sizeof(double) * G_full.size());
               f.write(reinterpret_cast<const char *>(GK_full.data()),
                       sizeof(double) * GK_full.size());
               f.write(reinterpret_cast<const char *>(Cmode.data()),
                       sizeof(double) * Cmode.size());
             });
  modes_saved = true;
}

void SubstructuringSolver::Impl::LoadSections(const std::string &path, std::streamoff pos,
                                              const SignatureMap &map)
{
  MPI_Comm comm = parent_fes.GetComm();
  const bool root = (Mpi::Rank(comm) == 0);
  const int nG = nG_global;
  std::ifstream f;
  std::streamoff size = 0;
  if (root)
  {
    f.open(path, std::ios::binary | std::ios::ate);
    size = f.tellg();
    f.seekg(pos);
  }
  // Abort unless the rest of the file holds the next bytes of a section (counted on rank
  // 0 from the sizes it read, in floating point: any size read). Collective.
  auto holds = [&](double bytes)
  {
    int ok = 1;
    if (root)
    {
      ok = f && bytes >= 0.0 &&
           static_cast<double>(size - static_cast<std::streamoff>(f.tellg())) >= bytes;
    }
    Mpi::Broadcast(1, &ok, 0, comm);
    MFEM_VERIFY(ok, "The saved substructuring model \""
                        << path << "\" is incomplete: rerun in \"Offline\" mode!");
  };
  constexpr double I = sizeof(int), D = sizeof(double);
  // Saved interface rows of width K, re-ordered onto the online numbering (rank 0).
  auto read_rows = [&](int K, std::vector<double> &full)
  {
    if (root)
    {
      std::vector<double> off(static_cast<std::size_t>(nG) * K);
      f.read(reinterpret_cast<char *>(off.data()), sizeof(double) * off.size());
      full = map.DualRows(off, K);
    }
  };
  while (true)
  {
    int magic = 0;
    if (root && !f.read(reinterpret_cast<char *>(&magic), sizeof(int)))
    {
      magic = 0;
    }
    MPI_Bcast(&magic, 1, MPI_INT, 0, comm);
    if (magic == kLengthMagic || magic == kEpsMagic || magic == kEnvMeshMagic)
    {
      const int n = (magic == kEnvMeshMagic) ? 12 : 1;
      holds(D * n);
      if (root)
      {
        f.seekg(sizeof(double) * n, std::ios::cur);  // read by ReadSavedHead
      }
    }
    else if (magic == kEnvFpMagic)
    {
      holds(D * 5);
      if (root)
      {
        f.read(reinterpret_cast<char *>(saved_env_fp.data()),
               sizeof(double) * saved_env_fp.size());
      }
      MPI_Bcast(saved_env_fp.data(), 5, MPI_DOUBLE, 0, comm);
      have_env_fp = true;
    }
    else if (magic == kSKenMagic)
    {
      holds(D * nG * static_cast<double>(nG));
      std::vector<double> full;
      if (root)
      {
        full = map.DualMatrix<double>(
            [&](std::size_t, double *x)
            { f.read(reinterpret_cast<char *>(x), sizeof(double) * nG); });
      }
      ScatterRows(full, nG, SK_rows);
      have_sk = true;
    }
    else if (magic == kSourceModesMagic)
    {
      holds(I);
      int K = 0;
      if (root)
      {
        f.read(reinterpret_cast<char *>(&K), sizeof(int));
      }
      // ids, fingerprints, g and h rows, c and B.
      const double k = K;
      holds((k < 0.0) ? -1.0 : I * k + D * (2.0 * k + 2.0 * k * nG + 2.0 * k * k));
      MPI_Bcast(&K, 1, MPI_INT, 0, comm);
      std::vector<int> ids(K);
      std::vector<double> fp(2 * K), C(static_cast<std::size_t>(K) * K), B(C.size()),
          G_full, H_full;
      if (root)
      {
        f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * K);
        f.read(reinterpret_cast<char *>(fp.data()), sizeof(double) * 2 * K);
      }
      read_rows(K, G_full);
      read_rows(K, H_full);
      if (root)
      {
        f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
        f.read(reinterpret_cast<char *>(B.data()), sizeof(double) * B.size());
      }
      MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
      MPI_Bcast(fp.data(), 2 * K, MPI_DOUBLE, 0, comm);
      MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
      MPI_Bcast(B.data(), K * K, MPI_DOUBLE, 0, comm);
      ScatterRows(G_full, K, SG_rows);
      ScatterRows(H_full, K, SH_rows);
      src_ids = ids;
      src_fp = fp;
      SC = C;
      SB = B;
      src_ready = true;
      src_saved = true;  // already in the file
    }
    else if (magic == kModesMagic)
    {
      holds(I);
      int hdr[2] = {0, 0};  // K, has_gk
      if (root)
      {
        f.read(reinterpret_cast<char *>(&hdr[0]), sizeof(int));
      }
      // ids, fingerprints and has_gk, then the G (and G^K) rows and Cmode.
      holds((hdr[0] < 0) ? -1.0 : I * hdr[0] + D * 2.0 * hdr[0] + I);
      std::vector<int> ids;
      std::vector<double> fp;
      if (root)
      {
        ids.resize(hdr[0]);
        fp.resize(2 * hdr[0]);
        f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * hdr[0]);
        f.read(reinterpret_cast<char *>(fp.data()), sizeof(double) * 2 * hdr[0]);
        f.read(reinterpret_cast<char *>(&hdr[1]), sizeof(int));
      }
      const double k = hdr[0];
      holds(D * (k * nG * (hdr[1] ? 2.0 : 1.0) + k * k));
      MPI_Bcast(hdr, 2, MPI_INT, 0, comm);
      const int K = hdr[0];
      ids.resize(K);
      fp.resize(2 * K);
      std::vector<double> C(static_cast<std::size_t>(K) * K), G_full, GK_full;
      read_rows(K, G_full);
      if (hdr[1])
      {
        read_rows(K, GK_full);
      }
      if (root)
      {
        f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
      }
      MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
      MPI_Bcast(fp.data(), 2 * K, MPI_DOUBLE, 0, comm);
      MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
      ScatterRows(G_full, K, G_rows);
      GK_rows.clear();
      have_gk = (hdr[1] != 0);
      if (have_gk)
      {
        ScatterRows(GK_full, K, GK_rows);
      }
      mode_ids = ids;
      mode_fp = fp;
      Cmode = C;
      modes_ready = true;
      modes_saved = true;  // already in the file
    }
    else
    {
      break;
    }
  }
}

void SubstructuringSolver::Impl::LoadModel(const std::string &model_path, int sig_w)
{
  const int nG = nG_global;
  const int rank = Mpi::Rank(parent_fes.GetComm());
  BlockTimer bt_io(Timer::IO);
  MFEM_VERIFY(saved.nG == nG,
              "Saved substructuring model interface size ("
                  << saved.nG << ") does not match this run (" << nG
                  << "); the interface (Gamma) must be identical between the offline and "
                     "online runs.");
  MFEM_VERIFY(saved.width == sig_w,
              "The interface of the saved substructuring model "
                  << ((saved.width == 13) ? "is" : "is not")
                  << " split by \"TwoSided\" superconductor films, unlike in this run; "
                     "the interface (Gamma) must be identical between the offline and "
                     "online runs.");
  // Signature match: the saved interface DOFs on the current ones (dual quantities, S_E
  // and the sections' interface rows, transform as x_cur = M^-T x_saved). H1: a
  // permutation. H(curl): signed matches, and blocks for the DOFs whose basis depends on
  // the face orientation (second-order face DOFs on tetrahedra).
  const SignatureMap map =
      MatchSignatureBasis(InterfaceSignatures(sig_w), saved.sig, saved.width);
  {
    // Rank 0 reads the saved S_E onto the current numbering and scatters its rows.
    std::vector<double> S_full;
    if (rank == 0)
    {
      std::ifstream f(model_path, std::ios::binary);
      f.seekg(saved.S_pos);
      S_full = map.DualMatrix<double>(
          [&](std::size_t, double *x)
          { f.read(reinterpret_cast<char *>(x), sizeof(double) * nG); });
    }
    ScatterRows(S_full, nG, S_rows);
  }
  // Sections after S_E (fingerprint, S^K, modes).
  LoadSections(model_path, saved.sections_pos, map);
  if (region_only)
  {
    Mpi::Print(" The mesh holds the region (the environment is that of the model)\n");
  }
  else
  {
    CheckEnvFingerprint();
  }
}

void SubstructuringSolver::Impl::SaveModel(const std::string &model_path, int sig_w)
{
  const int nG = nG_global;
  BlockTimer bt_io(Timer::IO);
  const int sig_type = magnetostatic ? ((sig_w == 13) ? 3 : 2) : 1;
  const std::vector<double> sig =
      InterfaceSignatures(sig_w);  // collective: all ranks participate
  const std::vector<double> S_full = GatherRows(S_rows, nG);
  WriteModel(
      model_path, false,
      [&](std::ofstream &f)
      {
        const int head[kHeadInts] = {kModelMagic, kModelVersion, iodata.solver.order, nG,
                                     sig_type};
        f.write(reinterpret_cast<const char *>(head), sizeof(head));
        f.write(reinterpret_cast<const char *>(sig.data()), sizeof(double) * nG * sig_w);
        f.write(reinterpret_cast<const char *>(S_full.data()), sizeof(double) * nG * nG);
        f.write(reinterpret_cast<const char *>(&kLengthMagic), sizeof(int));
        f.write(reinterpret_cast<const char *>(&iodata.model.Lc), sizeof(double));
      });
  AppendEnvFingerprint(model_path);  // collective
  AppendEnvMoments(model_path, WriteRegionMesh(RegionMeshPath(model_path)));
  if (have_sk)
  {
    AppendSK(model_path);  // collective
  }
  if (modes_ready)
  {
    AppendModes(model_path);  // collective
  }
}

void SubstructuringSolver::UseModelLength(IoData &iodata, MPI_Comm comm)
{
  const auto &sub = iodata.solver.substructuring;
  if (!sub || sub->mode != SubstructuringMode::ONLINE || sub->save_model.empty())
  {
    return;
  }
  const double L = Impl::ReadSavedHead(sub->save_model, comm).length;
  if (L > 0.0)
  {
    MFEM_VERIFY(!(iodata.model.Lc > 0.0) || std::abs(iodata.model.Lc - L) <= 1.0e-12 * L,
                "\"Model\": \"Lc\" (" << iodata.model.Lc
                                      << ") differs from that of the saved substructuring "
                                         "model ("
                                      << L << ")!");
    iodata.model.Lc = L;
  }
}

}  // namespace palace
