// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "substructuringsolverimpl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <system_error>
#include <utility>
#include "utils/outputdir.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The doubles of a chunk of rows passed through rank 0 (8 MB).
constexpr std::size_t kIoChunk = std::size_t(1) << 20;

}  // namespace

void SubstructuringSolver::Impl::WriteModel(
    const std::string &path, bool append,
    const std::function<void(std::ofstream &)> &write) const
{
  ApplyOnRootFilesystem(
      parent_fes.GetComm(), fmt::format("Writing the substructuring model \"{}\"", path),
      [&]()
      {
        std::ofstream f(path,
                        std::ios::binary | (append ? std::ios::app : std::ios::trunc));
        if (!f)
        {
          throw fs::filesystem_error("Could not open the model file", path,
                                     std::make_error_code(std::errc::io_error));
        }
        write(f);
        f.close();
        if (f.fail())
        {
          throw fs::filesystem_error("Could not write the model file", path,
                                     std::make_error_code(std::errc::io_error));
        }
      });
}

namespace substructuring
{

void RowsToRoot(MPI_Comm comm, const std::vector<double> &rows, int w,
                const std::vector<int> &cnt, std::size_t chunk,
                const std::function<void(const double *, int)> &write)
{
  const int per = static_cast<int>(std::max<std::size_t>(1, chunk / w));
  const RowType row(w);
  if (!Mpi::Root(comm))
  {
    const int m = static_cast<int>(rows.size() / w);
    for (int r0 = 0; r0 < m; r0 += per)
    {
      MPI_Send(rows.data() + static_cast<std::size_t>(r0) * w, std::min(per, m - r0),
               row.type, 0, 0, comm);
    }
    return;
  }
  // The chunks of another rank as it sends them, until it has sent its rows.
  std::vector<double> buf;
  for (std::size_t q = 0; q < cnt.size(); q++)
  {
    for (int r0 = 0; r0 < cnt[q];)
    {
      int n = std::min(per, cnt[q] - r0);
      const double *x = rows.data() + static_cast<std::size_t>(r0) * w;
      if (q > 0)
      {
        MPI_Status status;
        MPI_Probe(static_cast<int>(q), 0, comm, &status);
        MPI_Get_count(&status, row.type, &n);
        buf.resize(static_cast<std::size_t>(n) * w);
        MPI_Recv(buf.data(), n, row.type, static_cast<int>(q), 0, comm, MPI_STATUS_IGNORE);
        x = buf.data();
      }
      write(x, n);
      r0 += n;
    }
  }
}

RowsFromRoot::RowsFromRoot(MPI_Comm comm, int w, std::size_t chunk)
  : comm(comm), w(w), per(static_cast<int>(std::max<std::size_t>(1, chunk / w))), row(w)
{
}

void RowsFromRoot::Serve(const std::function<std::vector<int>(int)> &need,
                         const std::function<void(int, double *)> &read)
{
  for (int q = 1; q < Mpi::Size(comm); q++)
  {
    const std::vector<int> a = need(q);
    for (std::size_t i0 = 0; i0 < a.size(); i0 += per)
    {
      const int m = static_cast<int>(std::min<std::size_t>(per, a.size() - i0));
      buf.resize(static_cast<std::size_t>(m) * w);
      for (int k = 0; k < m; k++)
      {
        read(a[i0 + k], buf.data() + static_cast<std::size_t>(k) * w);
      }
      MPI_Send(buf.data(), m, row.type, q, 0, comm);
    }
  }
}

void RowsFromRoot::Next(double *x)
{
  if (i == n)
  {
    MPI_Status status;
    MPI_Probe(0, 0, comm, &status);
    MPI_Get_count(&status, row.type, &n);
    buf.resize(static_cast<std::size_t>(n) * w);
    MPI_Recv(buf.data(), n, row.type, 0, 0, comm, MPI_STATUS_IGNORE);
    i = 0;
  }
  std::copy_n(buf.data() + static_cast<std::size_t>(i++) * w, w, x);
}

}  // namespace substructuring

void SubstructuringSolver::Impl::AppendRows(const std::string &path,
                                            const std::vector<double> &rows) const
{
  MPI_Comm comm = parent_fes.GetComm();
  const int nG = nG_global;
  if (nG == 0)
  {
    return;
  }
  std::vector<int> cnt, disp;
  RowLayout(cnt, disp);
  if (!Mpi::Root(comm))
  {
    substructuring::RowsToRoot(comm, rows, nG, cnt, kIoChunk, {});
  }
  ApplyOnRootFilesystem(
      comm, fmt::format("Writing the substructuring model \"{}\"", path),
      [&]()
      {
        // A failed stream writes nothing.
        std::ofstream f(path, std::ios::binary | std::ios::app);
        substructuring::RowsToRoot(
            comm, rows, nG, cnt, kIoChunk, [&](const double *x, int n)
            { f.write(reinterpret_cast<const char *>(x), sizeof(double) * n * nG); });
        f.close();
        if (f.fail())
        {
          throw fs::filesystem_error("Could not write the model file", path,
                                     std::make_error_code(std::errc::io_error));
        }
      });
}

std::vector<double> SubstructuringSolver::Impl::ReadRows(const std::string &path,
                                                         std::streamoff pos,
                                                         const SignatureMap &map) const
{
  // Rank 0 sends every row even if the file fails, then tells whether it did.
  MPI_Comm comm = parent_fes.GetComm();
  const int nG = nG_global;
  std::vector<double> rows;
  if (nG == 0)
  {
    return rows;
  }
  std::vector<int> cnt, disp;
  RowLayout(cnt, disp);
  substructuring::RowsFromRoot from(comm, nG, kIoChunk);
  int ok = 1;
  if (Mpi::Root(comm))
  {
    std::ifstream f(path, std::ios::binary);
    auto read = [&](int a, double *x)
    {
      f.seekg(pos + static_cast<std::streamoff>(sizeof(double)) * a * nG);
      f.read(reinterpret_cast<char *>(x), sizeof(double) * nG);
    };
    from.Serve([&](int q) { return map.SavedRows(disp[q], cnt[q]); }, read);
    rows = map.DualMatrix<double>(gamma_off, gamma_nloc, read);
    ok = static_cast<bool>(f);
  }
  else
  {
    rows = map.DualMatrix<double>(gamma_off, gamma_nloc,
                                  [&](int, double *x) { from.Next(x); });
  }
  Mpi::Broadcast(1, &ok, 0, comm);
  MFEM_VERIFY(ok, "Reading the substructuring model \"" << path << "\" failed!");
  return rows;
}

void SubstructuringSolver::Impl::AppendSourceModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(src_ids.size());
  const std::vector<double> G_full = GatherRows(SG_rows, K);
  const std::vector<double> H_full = GatherRows(SH_rows, K);
  const std::vector<double> Hkin_full =
      have_skin ? GatherRows(SHkin_rows, K) : std::vector<double>();
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
          f.seekg(sizeof(double) * kEnvFpSize, std::ios::cur);  // read by LoadSections
          head[5] = f && f.tellg() <= size;
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
    // Only the broadcast fields: rank 0 may have read more of an incomplete model.
    SavedHead incomplete;
    incomplete.version = h.version;
    incomplete.order = h.order;
    incomplete.nG = h.nG;
    incomplete.sig_type = h.sig_type;
    return incomplete;
  }
  Mpi::Broadcast(1, &h.length, 0, comm);
  Mpi::Broadcast(1, &h.regularization, 0, comm);
  offsets(h.nG, h.width);
  h.sig.resize(static_cast<std::size_t>(h.nG) * h.width);
  Mpi::Broadcast(static_cast<int>(h.sig.size()), h.sig.data(), 0, comm);
  Mpi::Broadcast(6, h.env_moments.data(), 0, comm);
  Mpi::Broadcast(6, h.region_moments.data(), 0, comm);
  return h;
}

const std::array<std::array<double, 3>, 3> &
SubstructuringSolver::Impl::FingerprintCenters() const
{
  if (!fp_centers)
  {
    double box[6] = {
        std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    mfem::Vector c;
    for (int e = 0; e < parent.GetNE(); e++)
    {
      if (ea_arr.Find(parent.GetAttribute(e)) >= 0)
      {
        parent.GetElementCenter(e, c);
        for (int d = 0; d < c.Size(); d++)
        {
          box[d] = std::min(box[d], c(d));
          box[3 + d] = std::min(box[3 + d], -c(d));
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, box, 6, MPI_DOUBLE, MPI_MIN, parent_fes.GetComm());
    constexpr double f[3][3] = {{0.21, 0.37, 0.16}, {0.62, 0.28, 0.41}, {0.35, 0.73, 0.57}};
    fp_centers.emplace();
    for (int k = 0; k < 3; k++)
    {
      for (int d = 0; d < 3; d++)
      {
        const double lo = box[d], hi = -box[3 + d];
        (*fp_centers)[k][d] = (lo <= hi) ? lo + f[k][d] * (hi - lo) : 0.0;
      }
    }
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
  // The Dirichlet DOFs of the parent space (the same on both sides of the two-sided films).
  std::vector<char> dbc(parent_fes.GetTrueVSize(), 0);
  for (int i = 0; i < nt; i++)
  {
    dbc[ParentTrueDof(i)] |= env_dbc[i];
  }
  const std::array<double, 6> sums =
      TrueDofPointSums(parent_fes, dbc, FingerprintCenters());
  std::copy(sums.begin(), sums.end(), fp.begin() + 5);
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
  const std::array<double, kEnvFpSize> fp = EnvironmentFingerprint();  // collective
  BlockTimer bt(Timer::IO);
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
  std::vector<char> written = palace::WriteRegionMesh(
      path, parent_fes, sub.region_attributes, sub.environment_attributes, keep,
      iodata.units.GetMeshLengthRelativeScale(), plain);
  Mpi::Print(" Writing the mesh of the region for online runs to {}\n", path);
  return written;
}

void SubstructuringSolver::Impl::CheckEnvFingerprint()
{
  const std::array<double, kEnvFpSize> fp = EnvironmentFingerprint();  // collective
  MFEM_VERIFY(have_env_fp, "The saved substructuring model has no environment "
                           "fingerprint: rerun in \"Offline\" mode to condense it again!");
  // The counts exactly, the energies and the absolute sums relative to themselves, and the
  // sums relative to the absolute sums (they may cancel).
  bool ok = (fp[0] == saved_env_fp[0] && fp[1] == saved_env_fp[1]);
  for (int q = 2; q < kEnvFpSize; q++)
  {
    const int scale = (q >= 5 && q < 8) ? q + 3 : q;
    ok = ok &&
         std::abs(fp[q] - saved_env_fp[q]) <=
             kFingerprintTol * std::max(std::abs(fp[scale]), std::abs(saved_env_fp[scale]));
  }
  MFEM_VERIFY(ok, kEnvironmentDiffers);
}

void SubstructuringSolver::Impl::AppendOperator(const std::string &path, int magic,
                                                const std::vector<double> &rows) const
{
  BlockTimer bt(Timer::IO);
  WriteModel(path, true, [&](std::ofstream &f)
             { f.write(reinterpret_cast<const char *>(&magic), sizeof(int)); });
  AppendRows(path, rows);
}

void SubstructuringSolver::Impl::AppendModes(const std::string &path)
{
  BlockTimer bt(Timer::IO);
  const int K = static_cast<int>(mode_ids.size());
  const std::vector<double> G_full = GatherRows(G_rows, K);
  WriteModel(path, true,
             [&](std::ofstream &f)
             {
               f.write(reinterpret_cast<const char *>(&kModesMagic), sizeof(int));
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
      holds(D * kEnvFpSize);
      if (root)
      {
        f.read(reinterpret_cast<char *>(saved_env_fp.data()),
               sizeof(double) * saved_env_fp.size());
      }
      MPI_Bcast(saved_env_fp.data(), kEnvFpSize, MPI_DOUBLE, 0, comm);
      have_env_fp = true;
    }
    else if (magic == kSKenMagic || magic == kSKinMagic)
    {
      holds(D * nG * static_cast<double>(nG));
      std::streamoff at = 0;
      if (root)
      {
        at = f.tellg();
        f.seekg(static_cast<std::streamoff>(sizeof(double)) * nG * nG, std::ios::cur);
      }
      Mpi::Broadcast(1, &at, 0, comm);
      if (magic == kSKenMagic)
      {
        SK_rows = ReadRows(path, at, map);
        have_sk = true;
      }
      else
      {
        SKin_rows = ReadRows(path, at, map);
        have_skin = true;
      }
    }
    else if (magic == kSourceModesMagic)
    {
      holds(I);
      int K = 0;
      if (root)
      {
        f.read(reinterpret_cast<char *>(&K), sizeof(int));
      }
      // ids, fingerprints, g and h rows, c and B, and with S^Kin the h^kin rows and c^kin.
      const double k = K;
      holds((k < 0.0) ? -1.0
                      : I * k + D * (2.0 * k + 2.0 * k * nG + 2.0 * k * k +
                                     (have_skin ? k * nG + k * k : 0.0)));
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
      if (have_skin)
      {
        std::vector<double> Hkin_full, Ckin(C.size());
        read_rows(K, Hkin_full);
        if (root)
        {
          f.read(reinterpret_cast<char *>(Ckin.data()), sizeof(double) * Ckin.size());
        }
        MPI_Bcast(Ckin.data(), K * K, MPI_DOUBLE, 0, comm);
        ScatterRows(Hkin_full, K, SHkin_rows);
        SCkin = std::move(Ckin);
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
    }
    else if (magic == kModesMagic)
    {
      holds(I);
      int K = 0;
      if (root)
      {
        f.read(reinterpret_cast<char *>(&K), sizeof(int));
      }
      // ids, fingerprints, the G rows and Cmode.
      const double k = K;
      holds((k < 0.0) ? -1.0 : I * k + D * (2.0 * k + k * nG + k * k));
      MPI_Bcast(&K, 1, MPI_INT, 0, comm);
      std::vector<int> ids(K);
      std::vector<double> fp(2 * K), C(static_cast<std::size_t>(K) * K), G_full;
      if (root)
      {
        f.read(reinterpret_cast<char *>(ids.data()), sizeof(int) * K);
        f.read(reinterpret_cast<char *>(fp.data()), sizeof(double) * 2 * K);
      }
      read_rows(K, G_full);
      if (root)
      {
        f.read(reinterpret_cast<char *>(C.data()), sizeof(double) * C.size());
      }
      MPI_Bcast(ids.data(), K, MPI_INT, 0, comm);
      MPI_Bcast(fp.data(), 2 * K, MPI_DOUBLE, 0, comm);
      MPI_Bcast(C.data(), K * K, MPI_DOUBLE, 0, comm);
      ScatterRows(G_full, K, G_rows);
      mode_ids = ids;
      mode_fp = fp;
      Cmode = C;
      modes_ready = true;
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
    S_rows = ReadRows(model_path, saved.S_pos, map);
  }
  {
    // Sections after S_E (fingerprint, S^K, modes).
    BlockTimer bt_io(Timer::IO);
    LoadSections(model_path, saved.sections_pos, map);
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
    WriteModel(model_path, false,
               [&](std::ofstream &f)
               {
                 const int head[kHeadInts] = {kModelMagic, kModelVersion,
                                              iodata.solver.order, nG, sig_type};
                 f.write(reinterpret_cast<const char *>(head), sizeof(head));
                 f.write(reinterpret_cast<const char *>(sig.data()),
                         sizeof(double) * nG * sig_w);
               });
    AppendRows(model_path, S_rows);  // collective
    WriteModel(model_path, true,
               [&](std::ofstream &f)
               {
                 f.write(reinterpret_cast<const char *>(&kLengthMagic), sizeof(int));
                 f.write(reinterpret_cast<const char *>(&iodata.model.Lc), sizeof(double));
               });
  }
  AppendEnvFingerprint(model_path);  // collective
  AppendEnvMoments(model_path, WriteRegionMesh(RegionMeshPath(model_path)));
  if (have_sk)
  {
    AppendOperator(model_path, kSKenMagic, SK_rows);  // collective
    if (have_skin)
    {
      AppendOperator(model_path, kSKinMagic, SKin_rows);
    }
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
