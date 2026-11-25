// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <utility>
#include "linalg/lowertriangle.hpp"
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The doubles of a chunk of rows passed through rank 0 (8 MB).
constexpr std::size_t kIoChunk = std::size_t(1) << 20;

// The payload of a record of modes: K, then its ids and fingerprints, and its interface
// rows of width K and K x K blocks in turn, read with Next.
struct Modes
{
  std::vector<char> p;
  std::size_t at = 0;
  int K = 0;
  template <typename T>
  std::vector<T> Next(std::size_t n)
  {
    std::vector<T> x(n);
    std::memcpy(x.data(), p.data() + at, sizeof(T) * n);
    at += sizeof(T) * n;
    return x;
  }
};

}  // namespace

void SubstructuringSolver::Impl::AppendRecord(
    const std::string &path, std::uint32_t tag, long long size,
    const std::function<void(std::ostream &)> &write) const
{
  ModelFile::Append(path, parent_fes.GetComm(), tag, size, write);
}

void SubstructuringSolver::Impl::AppendTriangle(const std::string &path, std::uint32_t tag,
                                                const std::vector<double> &local) const
{
  MPI_Comm comm = parent_fes.GetComm();
  const int nG = nG_global;
  const std::vector<double> share = ShareFromGrid(*gamma_grid, local);
  const std::vector<long long> size = TriangleShareSizes(nG, comm);
  ModelFile::Append(
      path, comm, tag, static_cast<long long>(sizeof(double)) * TriangleSize(nG),
      [&](std::ostream &f)
      {
        ToRoot<double>(comm, share, size, kIoChunk,
                       [&](const double *x, std::size_t n)
                       {
                         f.write(reinterpret_cast<const char *>(x),
                                 static_cast<std::streamsize>(sizeof(double) * n));
                       });
      },
      [&]() { ToRoot<double>(comm, share, size, kIoChunk, {}); });
}

std::vector<double> SubstructuringSolver::Impl::ReadTriangle(const std::string &path,
                                                             std::streamoff pos,
                                                             const SignatureMap &map) const
{
  MPI_Comm comm = parent_fes.GetComm();
  const int nG = nG_global;
  std::vector<std::vector<std::pair<int, double>>> to(nG);
  for (int g = 0; g < nG; g++)
  {
    for (const auto &[a, c] : map.rows[g])
    {
      to[a].emplace_back(g, c);
    }
  }
  const BlockCyclicGrid &grid = *gamma_grid;
  std::vector<double> local(grid.LocalSize(), 0.0);
  std::ifstream f(path, std::ios::binary);
  bool ok = f.good();
  RedistributeTriangle<double>(
      comm, nG,
      [&](long long e, long long m, double *x)
      {
        f.seekg(pos + static_cast<std::streamoff>(sizeof(double)) * e);
        f.read(reinterpret_cast<char *>(x),
               static_cast<std::streamsize>(sizeof(double) * m));
        ok = ok && f.good();
      },
      to, [&](int g, int h) { return grid.Owner(g, h); },
      [&](int g, int h) { return grid.LocalIndex(g, h); }, local);
  int all = ok;
  Mpi::GlobalMin(1, &all, comm);
  MFEM_VERIFY(all, "Reading the substructuring model \"" << path << "\" failed!");
  return local;
}

void SubstructuringSolver::Impl::AppendSourceModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(src_ids.size());
  const std::vector<double> G_full = GatherRows(SG_rows, K);
  const std::vector<double> H_full = GatherRows(SH_rows, K);
  const std::vector<double> Hkin_full =
      have_skin ? GatherRows(SHkin_rows, K) : std::vector<double>();
  const long long k = K, n = nG_global;
  const long long size =
      static_cast<long long>(sizeof(int)) * (1 + k) +
      static_cast<long long>(sizeof(double)) *
          (2 * k + 2 * n * k + 2 * k * k + (have_skin ? n * k + k * k : 0));
  AppendRecord(
      path, kSourceModesTag, size,
      [&](std::ostream &f)
      {
        f.write(reinterpret_cast<const char *>(&K), sizeof(int));
        f.write(reinterpret_cast<const char *>(src_ids.data()), sizeof(int) * K);
        f.write(reinterpret_cast<const char *>(src_fp.data()), sizeof(double) * 2 * K);
        f.write(reinterpret_cast<const char *>(G_full.data()),
                sizeof(double) * G_full.size());
        f.write(reinterpret_cast<const char *>(H_full.data()),
                sizeof(double) * H_full.size());
        f.write(reinterpret_cast<const char *>(SC.data()), sizeof(double) * SC.size());
        f.write(reinterpret_cast<const char *>(SB.data()), sizeof(double) * SB.size());
        if (have_skin)
        {
          f.write(reinterpret_cast<const char *>(Hkin_full.data()),
                  sizeof(double) * Hkin_full.size());
          f.write(reinterpret_cast<const char *>(SCkin.data()),
                  sizeof(double) * SCkin.size());
        }
      });
}

int SubstructuringSolver::Impl::SignatureFileWidth(int sig_type)
{
  return (sig_type == 3) ? 13 : (sig_type == 2) ? 12 : (sig_type == 1) ? 3 : 0;
}

SubstructuringSolver::Impl::SavedHead
SubstructuringSolver::Impl::ReadSavedHead(const std::string &path, MPI_Comm comm)
{
  SavedHead h;
  h.file = ModelFile::Read(path, comm);
  h.version = h.file.version;
  if (h.version <= 0)
  {
    return h;
  }
  if (h.version != static_cast<int>(ModelFile::kVersion))
  {
    h.version = kModelVersion + 1;
    return h;
  }
  // The payload of the record tag, if it has n values of type T.
  auto values = [&]<typename T>(std::uint32_t tag, std::size_t n, T *x)
  {
    const auto *r = h.file.Find(tag);
    if (!r || r->size != static_cast<long long>(sizeof(T) * n))
    {
      return false;
    }
    const std::vector<char> p = ModelFile::ReadPayload(path, *r, comm);
    if (p.size() != sizeof(T) * n)
    {
      return false;
    }
    std::memcpy(x, p.data(), p.size());
    return true;
  };
  int head[4] = {0, 0, 0, 0};
  h.version = kModelVersion;  // a model without its head is incomplete
  if (!values(kHeadTag, 4, head))
  {
    return h;
  }
  h.version = head[0];
  h.order = head[1];
  h.nG = head[2];
  h.sig_type = head[3];
  if (h.version != kModelVersion)
  {
    return h;
  }
  // A complete model holds its signatures, S_E, environment fingerprint and record.
  const int w = SignatureFileWidth(h.sig_type);
  const long long nG = h.nG;
  const auto *S = h.file.Find(kSETag);
  const auto *fp = h.file.Find(kEnvFpTag);
  const auto *env = h.file.Find(kEnvRecordTag);
  const auto environment =
      env ? EnvironmentRecord::Deserialize(ModelFile::ReadPayload(path, *env, comm))
          : std::nullopt;
  substructuring::InterfaceRecords iface;
  if (h.file.incomplete || w == 0 || nG < 0 || !S ||
      S->size != static_cast<long long>(sizeof(double)) * TriangleSize(nG) || !fp ||
      fp->size != static_cast<long long>(sizeof(double)) * kEnvFpSize || !environment ||
      !iface.Read(path, h.file, comm) || iface.n != nG || iface.width != w)
  {
    return h;
  }
  h.sig = std::move(iface.signatures);
  h.width = w;
  h.S_pos = S->pos;
  h.environment = *environment;
  values(kEpsTag, 1, &h.regularization);
  h.sheet_sides.resize(6);
  if (!values(kSheetSidesTag, 6, h.sheet_sides.data()))
  {
    h.sheet_sides.clear();
  }
  return h;
}

const std::array<std::array<double, 3>, 3> &
SubstructuringSolver::Impl::FingerprintCenters() const
{
  if (!fp_centers)
  {
    fp_centers = palace::FingerprintCenters(parent, env_attributes);
  }
  return *fp_centers;
}

mfem::Vector SubstructuringSolver::Impl::FingerprintField(int k) const
{
  mfem::Vector r = palace::FingerprintField(parent_fes, k, FingerprintCenters()[k]);
  if (!broken)
  {
    return r;
  }
  mfem::Vector rb(nt);
  all_versions->MultTranspose(r, rb);
  return rb;
}

std::array<double, SubstructuringSolver::Impl::kEnvFpSize>
SubstructuringSolver::Impl::EnvironmentFingerprint() const
{
  // The environment does not change within a run: computed once.
  if (env_fp)
  {
    return *env_fp;
  }
  double loc[2] = {0.0, 0.0};
  for (int i = 0; i < nt; i++)
  {
    loc[0] += (is_env_int[i] || is_gamma[i]) ? 1.0 : 0.0;
    loc[1] += env_dbc[i] ? 1.0 : 0.0;
  }
  const mfem::Vector r2 = FingerprintField(2);
  const std::vector<const mfem::Vector *> r = {&CachedFingerprintField(0),
                                               &CachedFingerprintField(1), &r2};
  const std::vector<mfem::Vector> Ar = EnvApply(r);
  std::array<double, kEnvFpSize> fp = {};
  MPI_Allreduce(loc, fp.data(), 2, MPI_DOUBLE, MPI_SUM, parent_fes.GetComm());
  for (int k = 0; k < 3; k++)
  {
    fp[2 + k] = mfem::InnerProduct(parent_fes.GetComm(), *r[k], Ar[k]);
  }
  env_fp = fp;
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

void SubstructuringSolver::Impl::AppendEnvRecord(const std::string &path,
                                                 const std::vector<char> &written)
{
  BlockTimer bt(Timer::IO);
  const std::vector<char> record =
      EnvironmentRecord::Make(parent, env_attributes, *iodata.resolved_config, written)
          .Serialize();  // collective
  AppendRecord(path, kEnvRecordTag, static_cast<long long>(record.size()),
               [&](std::ostream &f) { f.write(record.data(), record.size()); });
}

void SubstructuringSolver::Impl::AppendEnvFingerprint(const std::string &path)
{
  const std::array<double, kEnvFpSize> fp = EnvironmentFingerprint();  // collective
  BlockTimer bt(Timer::IO);
  AppendRecord(
      path, kEnvFpTag, sizeof(double) * fp.size(), [&](std::ostream &f)
      { f.write(reinterpret_cast<const char *>(fp.data()), sizeof(double) * fp.size()); });
  if (magnetostatic)
  {
    AppendRecord(path, kEpsTag, sizeof(double), [&](std::ostream &f)
                 { f.write(reinterpret_cast<const char *>(&mag_eps), sizeof(double)); });
    AppendRecord(path, kSheetSidesTag, sizeof(double) * gamma_sheet_sides.size(),
                 [&](std::ostream &f)
                 {
                   f.write(reinterpret_cast<const char *>(gamma_sheet_sides.data()),
                           sizeof(double) * gamma_sheet_sides.size());
                 });
  }
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
  // On a GPU, P lives on the device and its products are only valid there: read the
  // results with HostRead (an element access l(i) reads the host copy as it is).
  const double *ll = l.HostRead();
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
    if (std::any_of(dofs.begin(), dofs.end(), [&](int d) { return ll[index(d)] != 0.0; }))
    {
      for (int d : dofs)
      {
        lw(index(d)) = 1.0;
      }
    }
  }
  P.AbsMultTranspose(1.0, lw, 0.0, t);
  const double *tw = t.HostRead();
  std::vector<char> written(t.Size());
  for (int i = 0; i < t.Size(); i++)
  {
    written[i] = (tw[i] > 0.0);
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
  if (parent.Nonconforming())
  {
    // The submesh of a nonconforming mesh also selects by the attributes of the parents of
    // refined elements, which the written environment elements do not shift.
    Mpi::Warning(parent.GetComm(),
                 "No mesh of the region is written for a nonconforming mesh: online runs "
                 "need the full mesh!\n");
    return std::vector<char>(parent.GetNE(), 0);
  }
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
    // A flux loop of the environment (env_loops) has no DOF of its film and hole on a
    // region element (none on Gamma either), and no surface current or aperture on them:
    // its excitation, zero on the region, is in the model.
    const std::set<int> current_attrs(attrs.begin(), attrs.end());
    std::vector<char> region_dofs(parent_fes.GetTrueVSize(), 0);
    for (int t : rc_idx)
    {
      region_dofs[ParentTrueDof(t)] = 1;
    }
    env_loops.clear();
    for (const auto &[idx, data] : iodata.boundaries.fluxloop)
    {
      std::vector<int> loop_attrs(data.film_attributes);
      loop_attrs.insert(loop_attrs.end(), data.hole_attributes.begin(),
                        data.hole_attributes.end());
      int kept =
          std::ranges::any_of(loop_attrs, [&](int a) { return current_attrs.contains(a); });
      if (!kept)
      {
        const auto loop_mark = BoundaryTrueDofs(parent_fes, loop_attrs);
        for (std::size_t i = 0; i < loop_mark.size() && !kept; i++)
        {
          kept = loop_mark[i] && region_dofs[i];
        }
        Mpi::GlobalMax(1, &kept, parent.GetComm());
      }
      if (!kept)
      {
        env_loops.push_back(idx);
        continue;
      }
      attrs.insert(attrs.end(), loop_attrs.begin(), loop_attrs.end());
      loop_films.insert(data.film_attributes.begin(), data.film_attributes.end());
    }
    const auto mark = BoundaryTrueDofs(parent_fes, attrs);
    for (std::size_t i = 0; i < keep.size(); i++)
    {
      keep[i] = keep[i] || mark[i];
    }
    // The environment neighbors of the sheet faces on Gamma, so that each of these faces
    // has the same neighbors, and so belongs to the same side, as on the full mesh (see
    // CheckSheets).
    mfem::Vector lmark(parent_fes.GetVSize()), tmark(parent_fes.GetTrueVSize());
    lmark = 0.0;
    mfem::Array<int> dofs;
    for (int be : gamma_sheet_faces)
    {
      parent_fes.GetBdrElementDofs(be, dofs);
      for (int d : dofs)
      {
        lmark((d >= 0) ? d : -1 - d) = 1.0;
      }
    }
    parent_fes.Dof_TrueDof_Matrix()->AbsMultTranspose(1.0, lmark, 0.0, tmark);
    const double *t = tmark.HostRead();
    for (std::size_t i = 0; i < keep.size(); i++)
    {
      keep[i] = keep[i] || (t[i] > 0.0);
    }
  }
  const std::vector<char> plain =
      broken ? RegionFilmClosure(keep, loop_films) : std::vector<char>();
  const auto &sub = *iodata.solver.substructuring;
  std::vector<char> written =
      palace::WriteRegionMesh(path, parent_fes, sub.region_attributes, env_attributes, keep,
                              iodata.units.GetMeshLengthRelativeScale(), plain);
  Mpi::Print(" Writing the mesh of the region for online runs to {}\n", path);
  return written;
}

std::vector<int> SubstructuringSolver::Impl::ReadEnvLoops(const std::string &path,
                                                          MPI_Comm comm)
{
  const ModelFile file = ModelFile::Read(path, comm);
  const auto *r = (file.version == static_cast<int>(ModelFile::kVersion))
                      ? file.Find(kEnvLoopsTag)
                      : nullptr;
  if (!r)
  {
    return {};
  }
  const std::vector<char> p = ModelFile::ReadPayload(path, *r, comm);
  std::vector<int> loops(p.size() / sizeof(int));
  std::memcpy(loops.data(), p.data(), sizeof(int) * loops.size());
  return loops;
}

void SubstructuringSolver::Impl::AppendFluxoids(const std::string &path,
                                                const std::vector<int> &ids,
                                                const std::vector<int> &loop_ids,
                                                const mfem::DenseMatrix &fluxoid) const
{
  const int KL[2] = {static_cast<int>(ids.size()), static_cast<int>(loop_ids.size())};
  std::vector<double> values;
  values.reserve(static_cast<std::size_t>(KL[0]) * KL[1]);
  for (int k = 0; k < KL[0]; k++)
  {
    for (int j = 0; j < KL[1]; j++)
    {
      values.push_back(fluxoid(k, j));
    }
  }
  const long long size = static_cast<long long>(sizeof(int)) * (2 + KL[0] + KL[1]) +
                         static_cast<long long>(sizeof(double)) * values.size();
  AppendRecord(path, kFluxoidsTag, size,
               [&](std::ostream &f)
               {
                 f.write(reinterpret_cast<const char *>(KL), sizeof(KL));
                 f.write(reinterpret_cast<const char *>(ids.data()), sizeof(int) * KL[0]);
                 f.write(reinterpret_cast<const char *>(loop_ids.data()),
                         sizeof(int) * KL[1]);
                 f.write(reinterpret_cast<const char *>(values.data()),
                         sizeof(double) * values.size());
               });
}

void SubstructuringSolver::Impl::CheckEnvFingerprint()
{
  const std::array<double, kEnvFpSize> fp = EnvironmentFingerprint();  // collective
  MFEM_VERIFY(have_env_fp, "The saved substructuring model has no environment fingerprint: "
                           "rerun offline (without \"ModelFile\") to condense it again!");
  // The counts exactly, the energies relative to themselves.
  MFEM_VERIFY(FingerprintsMatch(fp, saved_env_fp, 2, kFingerprintTol), kEnvironmentDiffers);
}

void SubstructuringSolver::Impl::AppendModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(mode_ids.size());
  const std::vector<double> G_full = GatherRows(G_rows, K);
  const long long k = K, n = nG_global;
  AppendRecord(path, kModesTag,
               static_cast<long long>(sizeof(int)) * (1 + k) +
                   static_cast<long long>(sizeof(double)) * (2 * k + n * k + k * k),
               [&](std::ostream &f)
               {
                 f.write(reinterpret_cast<const char *>(&K), sizeof(int));
                 f.write(reinterpret_cast<const char *>(mode_ids.data()), sizeof(int) * K);
                 f.write(reinterpret_cast<const char *>(mode_fp.data()),
                         sizeof(double) * 2 * K);
                 f.write(reinterpret_cast<const char *>(G_full.data()),
                         sizeof(double) * G_full.size());
                 f.write(reinterpret_cast<const char *>(Cmode.data()),
                         sizeof(double) * Cmode.size());
               });
}

void SubstructuringSolver::Impl::LoadSections(const std::string &path,
                                              const SignatureMap &map)
{
  MPI_Comm comm = parent_fes.GetComm();
  const bool root = Mpi::Root(comm);
  const long long nG = nG_global;
  const auto &file = saved.file;
  auto complete = [&](bool ok)
  {
    MFEM_VERIFY(ok, "The saved substructuring model \""
                        << path
                        << "\" is incomplete: rerun offline (without \"ModelFile\")!");
  };
  {
    const std::vector<char> p = ModelFile::ReadPayload(path, *file.Find(kEnvFpTag), comm);
    complete(p.size() == sizeof(double) * kEnvFpSize);
    std::memcpy(saved_env_fp.data(), p.data(), p.size());
    have_env_fp = true;
  }
  for (const auto tag : {kSKenTag, kSKinTag})
  {
    if (const auto *r = file.Find(tag))
    {
      complete(r->size == static_cast<long long>(sizeof(double)) * TriangleSize(nG));
      (tag == kSKenTag ? SK_grid : SKin_grid) = ReadTriangle(path, r->pos, map);
      (tag == kSKenTag ? have_sk : have_skin) = true;
    }
  }

  // The modes of a record, with the interface rows re-ordered onto the online numbering
  // on rank 0.
  auto read_modes = [&](std::uint32_t tag, long long rows, long long blocks)
  {
    Modes m;
    if (const auto *r = file.Find(tag))
    {
      m.p = ModelFile::ReadPayload(path, *r, comm);
      complete(m.p.size() >= sizeof(int));
      std::memcpy(&m.K, m.p.data(), sizeof(int));
      m.at = sizeof(int);
      const long long k = m.K;
      complete(k >= 0 && static_cast<long long>(m.p.size()) ==
                             static_cast<long long>(sizeof(int)) * (1 + k) +
                                 static_cast<long long>(sizeof(double)) *
                                     (2 * k + rows * nG * k + blocks * k * k));
    }
    return m;
  };
  auto interface_rows = [&](Modes &m, std::vector<double> &rows)
  {
    const std::vector<double> saved_rows = m.Next<double>(nG * m.K);
    ScatterRows(root ? map.DualRows(saved_rows, m.K) : std::vector<double>(), m.K, rows);
  };
  if (Modes m = read_modes(kModesTag, 1, 1); !m.p.empty())
  {
    mode_ids = m.Next<int>(m.K);
    mode_fp = m.Next<double>(2 * m.K);
    interface_rows(m, G_rows);
    Cmode = m.Next<double>(static_cast<std::size_t>(m.K) * m.K);
    modes_ready = true;
  }
  if (Modes m = read_modes(kSourceModesTag, have_skin ? 3 : 2, have_skin ? 3 : 2);
      !m.p.empty())
  {
    src_ids = m.Next<int>(m.K);
    src_fp = m.Next<double>(2 * m.K);
    interface_rows(m, SG_rows);
    interface_rows(m, SH_rows);
    SC = m.Next<double>(static_cast<std::size_t>(m.K) * m.K);
    SB = m.Next<double>(static_cast<std::size_t>(m.K) * m.K);
    if (have_skin)
    {
      interface_rows(m, SHkin_rows);
      SCkin = m.Next<double>(static_cast<std::size_t>(m.K) * m.K);
    }
    src_ready = true;
  }
}

void SubstructuringSolver::Impl::LoadModel(const std::string &model_path, int sig_w)
{
  const int nG = nG_global;
  MFEM_VERIFY(saved.nG == nG,
              "Saved substructuring model interface size ("
                  << saved.nG << ") does not match this run (" << nG
                  << "); the interface (Gamma) must be identical between the offline and "
                     "online runs!");
  MFEM_VERIFY(saved.width == sig_w,
              "The interface of the saved substructuring model "
                  << ((saved.width == 13) ? "is" : "is not")
                  << " split by \"TwoSided\" superconductor films, unlike in this run; "
                     "the interface (Gamma) must be identical between the offline and "
                     "online runs!");
  // Signature match: the saved interface DOFs on the current ones (dual quantities, S_E
  // and the sections' interface rows, transform as x_cur = M^-T x_saved). H1: a
  // permutation. H(curl): signed matches, and blocks for the DOFs whose basis depends on
  // the face orientation (second-order face DOFs on tetrahedra).
  const SignatureMap map =
      MatchSignatureBasis(InterfaceSignatures(sig_w), saved.sig, saved.width);
  {
    BlockTimer bt_io(Timer::IO);
    S_grid = ReadTriangle(model_path, saved.S_pos, map);
  }
  {
    // The records after S_E (fingerprint, S^K, modes).
    BlockTimer bt_io(Timer::IO);
    LoadSections(model_path, map);
  }
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
  const int sig_type = magnetostatic ? ((sig_w == 13) ? 3 : 2) : 1;
  const std::vector<double> sig =
      InterfaceSignatures(sig_w);  // collective: all ranks participate
  {
    BlockTimer bt_io(Timer::IO);
    ModelFile::Create(model_path, parent_fes.GetComm());
    const int head[4] = {kModelVersion, iodata.solver.order, nG, sig_type};
    AppendRecord(model_path, kHeadTag, sizeof(head), [&](std::ostream &f)
                 { f.write(reinterpret_cast<const char *>(head), sizeof(head)); });
    substructuring::InterfaceRecords{nG, sig_w, sig, iodata.model.Lc}.Append(
        model_path, parent_fes.GetComm());
    AppendTriangle(model_path, kSETag, S_grid);  // collective
  }
  AppendEnvFingerprint(model_path);  // collective
  AppendEnvRecord(model_path, WriteRegionMesh(RegionMeshPath(model_path)));
  AppendRecord(model_path, kEnvLoopsTag,
               static_cast<long long>(sizeof(int) * env_loops.size()),
               [&](std::ostream &f)
               {
                 f.write(reinterpret_cast<const char *>(env_loops.data()),
                         sizeof(int) * env_loops.size());
               });
  if (have_sk)
  {
    BlockTimer bt_io(Timer::IO);
    AppendTriangle(model_path, kSKenTag, SK_grid);  // collective
    if (have_skin)
    {
      AppendTriangle(model_path, kSKinTag, SKin_grid);
    }
  }
  if (modes_ready)
  {
    AppendModes(model_path);  // collective
  }
}

void SubstructuringSolver::EnvironmentFluxoids(const std::vector<int> &ids,
                                               const std::vector<int> &loop_ids,
                                               mfem::DenseMatrix &fluxoid) const
{
  MFEM_VERIFY(fluxoid.Height() == static_cast<int>(ids.size()) &&
                  fluxoid.Width() == static_cast<int>(loop_ids.size()),
              "EnvironmentFluxoids: one row per excitation and one column per loop!");
  const auto &sub = *impl->iodata.solver.substructuring;
  if (!sub.Online())
  {
    if (const auto path = substructuring::SaveModelPath(impl->iodata); !path.empty())
    {
      impl->AppendFluxoids(path, ids, loop_ids, fluxoid);  // collective
    }
    return;
  }
  if (!impl->region_only || impl->env_loops.empty())
  {
    return;
  }
  using ModelFile = substructuring::ModelFile;
  MPI_Comm comm = impl->parent_fes.GetComm();
  const ModelFile file = ModelFile::Read(sub.model_file, comm);
  const auto *r = file.Find(Impl::kFluxoidsTag);
  const std::vector<char> p =
      r ? ModelFile::ReadPayload(sub.model_file, *r, comm) : std::vector<char>();
  // The record [K][L][excitation ids K][loop ids L][fluxoids K x L].
  std::size_t at = 0;
  auto read = [&]<typename T>(T *x, std::size_t n)
  {
    if (p.size() - at < sizeof(T) * n)
    {
      return false;
    }
    std::memcpy(x, p.data() + at, sizeof(T) * n);
    at += sizeof(T) * n;
    return true;
  };
  int KL[2] = {0, 0};
  std::vector<int> saved_ids, saved_loops;
  std::vector<double> values;
  bool ok = read(KL, 2) && KL[0] >= 0 && KL[1] >= 0;
  if (ok)
  {
    saved_ids.resize(KL[0]);
    saved_loops.resize(KL[1]);
    values.resize(static_cast<std::size_t>(KL[0]) * KL[1]);
    ok = read(saved_ids.data(), saved_ids.size()) &&
         read(saved_loops.data(), saved_loops.size()) &&
         read(values.data(), values.size()) && at == p.size();
  }
  MFEM_VERIFY(ok, "The saved substructuring model \""
                      << sub.model_file
                      << "\" holds no fluxoids of the flux loops of its environment: rerun "
                         "offline (without \"ModelFile\")!");
  auto index = [](const std::vector<int> &v, int x)
  { return static_cast<int>(std::find(v.begin(), v.end(), x) - v.begin()); };
  for (int loop : impl->env_loops)
  {
    const int j = index(loop_ids, loop), q = index(saved_loops, loop);
    MFEM_VERIFY(j < static_cast<int>(loop_ids.size()) && q < KL[1],
                "Flux loop " << loop << " of the environment of the saved substructuring "
                             << "model is not in its fluxoids: rerun offline (without "
                                "\"ModelFile\")!");
    for (std::size_t k = 0; k < ids.size(); k++)
    {
      const int row = index(saved_ids, ids[k]);
      MFEM_VERIFY(row < KL[0], "Flux loop excitation "
                                   << ids[k] << " is not in the fluxoids of the saved "
                                   << "substructuring model: rerun offline (without "
                                      "\"ModelFile\")!");
      fluxoid(static_cast<int>(k), j) = values[static_cast<std::size_t>(row) * KL[1] + q];
    }
  }
}

void SubstructuringSolver::TakeModelExcitations(
    const std::map<int, std::vector<int>> &excitations)
{
  MFEM_VERIFY(impl->mat_dtn,
              "CondenseEnvironment must be called before TakeModelExcitations!");
  impl->model_sources.clear();
  const bool single =
      std::ranges::all_of(excitations, [](const auto &x) { return x.second.size() == 1; });
  if (single && (impl->region_only || impl->env_excitations_unchanged))
  {
    auto has = [](const std::vector<int> &v, int x)
    { return std::ranges::find(v, x) != v.end(); };
    for (const auto &[idx, loops] : excitations)
    {
      if (has(impl->env_loops, loops.front()) && impl->src_ready && has(impl->src_ids, idx))
      {
        impl->model_sources.insert(idx);
      }
    }
  }
}

}  // namespace palace
