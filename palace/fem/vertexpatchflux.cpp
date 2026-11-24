// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "vertexpatchflux.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <Eigen/Dense>
#include <mfem.hpp>
#include "fem/brokenspace.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"

namespace palace::fem
{

namespace
{

// The nonzero entries of the local rows of a parallel matrix, with global column indices.
struct GlobalRows
{
  std::vector<int> offsets;
  std::vector<HYPRE_BigInt> cols;
  std::vector<double> vals;
};

GlobalRows GetGlobalRows(const mfem::HypreParMatrix &A)
{
  // The loops read the arrays of A on the host: A is returned to hypre's memory space (the
  // device in GPU builds) after them.
  A.HostRead();
  hypre_ParCSRMatrix *hA = const_cast<mfem::HypreParMatrix &>(A);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(hA),
                        *offd = hypre_ParCSRMatrixOffd(hA);
  const HYPRE_Int *di = hypre_CSRMatrixI(diag), *dj = hypre_CSRMatrixJ(diag),
                  *oi = hypre_CSRMatrixI(offd), *oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Real *dv = hypre_CSRMatrixData(diag), *ov = hypre_CSRMatrixData(offd);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(hA);
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(hA);
  GlobalRows rows;
  rows.offsets.assign(1, 0);
  for (int i = 0; i < hypre_CSRMatrixNumRows(diag); i++)
  {
    for (HYPRE_Int k = di[i]; k < di[i + 1]; k++)
    {
      if (dv[k] != 0.0)
      {
        rows.cols.push_back(col0 + dj[k]);
        rows.vals.push_back(dv[k]);
      }
    }
    for (HYPRE_Int k = oi[i]; k < oi[i + 1]; k++)
    {
      if (ov[k] != 0.0)
      {
        rows.cols.push_back(cmap[oj[k]]);
        rows.vals.push_back(ov[k]);
      }
    }
    rows.offsets.push_back(static_cast<int>(rows.cols.size()));
  }
  A.HypreRead();
  return rows;
}

// The first global column of each process, and the global number of columns last.
std::vector<HYPRE_BigInt> ColumnStarts(const mfem::HypreParMatrix &A, MPI_Comm comm)
{
  hypre_ParCSRMatrix *hA = const_cast<mfem::HypreParMatrix &>(A);
  HYPRE_BigInt first = hypre_ParCSRMatrixFirstColDiag(hA);
  std::vector<HYPRE_BigInt> starts(Mpi::Size(comm) + 1);
  MPI_Allgather(&first, 1, HYPRE_MPI_BIG_INT, starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  starts.back() = A.GetGlobalNumCols();
  return starts;
}

int Owner(const std::vector<HYPRE_BigInt> &starts, HYPRE_BigInt t)
{
  return static_cast<int>(std::upper_bound(starts.begin(), starts.end() - 1, t) -
                          starts.begin()) -
         1;
}

// A cell of a patch, read from a record of the cell data: for each face, the true DOF of
// rt0 it reads (-1 for a Neumann face) and the coefficient of that true DOF in the cell's
// basis function of the face; the ε⁻¹-weighted Raviart-Thomas mass matrix of the cell and
// the integrals of the divergences of its basis functions; then, for the hat, the faces on
// which the hat does not vanish (bits), and per potential the divergence datum
// ∫_K ∇ψₐ ⋅ D_v and the loads (ε⁻¹ ψₐ D_v, φᵢ)_K.
struct PatchCell
{
  int nf;
  const double *faces, *M, *d, *hat;
};

}  // namespace

std::vector<Vector> VertexPatchFlux(const FiniteElementSpace &primal,
                                    const FiniteElementSpace &hats,
                                    const FiniteElementSpace &rt0,
                                    const MaterialOperator &mat_op,
                                    const std::vector<int> &dirichlet_attr,
                                    const std::vector<Vector> &potentials)
{
  const auto &mesh = primal.GetParMesh();
  const MPI_Comm comm = primal.GetComm();
  const int nt = static_cast<int>(potentials.size()), np = Mpi::Size(comm),
            rank = Mpi::Rank(comm);
  MFEM_VERIFY(mesh.Dimension() == 3 && hats.GetMaxElementOrder() == 1 &&
                  rt0.GetMaxElementOrder() == 1,
              "Vertex patch fluxes require a 3D mesh, the order-one H1 space and the "
              "lowest-order Raviart-Thomas space!");

  // The rows of the prolongations, with the global true DOFs of their entries. Each row of
  // the prolongation of rt0 has a single entry (one DOF per face, constrained by at most
  // one master face).
  const auto &P_h = *hats.Dof_TrueDof_Matrix(), &P_f = *rt0.Dof_TrueDof_Matrix();
  const auto h_rows = GetGlobalRows(P_h), f_rows = GetGlobalRows(P_f);
  const auto h_starts = ColumnStarts(P_h, comm), f_starts = ColumnStarts(P_f, comm);
  for (int l = 0; l < rt0.GetVSize(); l++)
  {
    MFEM_VERIFY(f_rows.offsets[l + 1] - f_rows.offsets[l] <= 1,
                "Unexpected prolongation of the lowest-order Raviart-Thomas space!");
  }
  auto FluxEntry = [&f_rows](int l) -> std::pair<HYPRE_BigInt, double>
  {
    const int k = f_rows.offsets[l];
    return (f_rows.offsets[l + 1] > k) ? std::make_pair(f_rows.cols[k], f_rows.vals[k])
                                       : std::make_pair(HYPRE_BigInt(-1), 0.0);
  };

  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  const int ne = mesh.GetNE();

  // The L-vectors of the potentials.
  std::vector<Vector> v(nt);
  for (int t = 0; t < nt; t++)
  {
    v[t].SetSize(primal.GetVSize());
    v[t].UseDevice(false);
    primal.GetProlongationMatrix()->Mult(potentials[t], v[t]);
  }

  // Neumann faces: exterior faces without a Dirichlet attribute.
  const std::set<int> dirichlet(dirichlet_attr.begin(), dirichlet_attr.end());
  const auto face_to_bdr = mesh.GetFaceToBdrElMap();
  auto Neumann = [&](int f)
  {
    if (!mesh.GetFaceInformation(f).IsBoundary())
    {
      return false;
    }
    const int be = face_to_bdr[f];
    return be < 0 || !dirichlet.count(mesh.GetBdrAttribute(be));
  };

  // The records of the cell data, for the processes owning the hats of each cell: per
  // process and cell [nf, nf x (true DOF, coefficient), M (nf x nf), d (nf), number of
  // hats, then per hat (hat, face bits, nt x (datum, nf loads))].
  std::vector<std::vector<double>> send(np);
  {
    mfem::Array<int> hdofs, pdofs, ents, faces, orients;
    mfem::IsoparametricTransformation T;
    mfem::DenseMatrix vshape, hdshape, pdshape, M, eps_inv(3);
    mfem::Vector divshape, hshape, u, gu(3), Dv(3), es, d;
    std::vector<double> g, r;  // Per vertex and potential: datum, loads.
    std::vector<int> face_of_dof;
    for (int e = 0; e < ne; e++)
    {
      const auto &fe_f = *rt0.Get().GetFE(e), &fe_h = *hats.Get().GetFE(e),
                 &fe_p = *primal.Get().GetFE(e);
      const int nf = fe_f.GetDof(), nv = fe_h.GetDof(), npd = fe_p.GetDof();
      const auto *elem = mesh.GetElement(e);

      // Faces: the true DOF read and its coefficient (with the orientation of the element).
      rt0.GetElementDofs(e, dofs, dof_trans);
      GetElementDofEntities(rt0.Get(), e, ents);
      mesh.GetElementFaces(e, faces, orients);
      const int nve = elem->GetNVertices(), nee = elem->GetNEdges();
      face_of_dof.assign(nf, -1);
      std::vector<double> face_data(2 * nf);
      for (int i = 0; i < nf; i++)
      {
        const int l = (dofs[i] >= 0) ? dofs[i] : -1 - dofs[i];
        const auto [t, c] = FluxEntry(l);
        face_of_dof[i] = ents[i] - nve - nee;
        MFEM_VERIFY(face_of_dof[i] >= 0 && face_of_dof[i] < faces.Size(),
                    "Unexpected lowest-order Raviart-Thomas element DOFs!");
        const bool neumann = Neumann(faces[face_of_dof[i]]);
        face_data[2 * i] = (neumann || t < 0) ? -1.0 : static_cast<double>(t);
        face_data[2 * i + 1] = ((dofs[i] >= 0) ? 1.0 : -1.0) * c;
      }

      // The hats on the cell: the entries of the rows of its vertex L-DOFs.
      hats.Get().GetElementDofs(e, hdofs);
      GetElementDofEntities(hats.Get(), e, ents);
      std::map<HYPRE_BigInt, std::vector<double>> cell_hats;  // Hat -> vertex weights.
      for (int j = 0; j < nv; j++)
      {
        MFEM_VERIFY(ents[j] == j && hdofs[j] >= 0, "Unexpected order-one H1 element DOFs!");
        for (int k = h_rows.offsets[hdofs[j]]; k < h_rows.offsets[hdofs[j] + 1]; k++)
        {
          auto &w = cell_hats[h_rows.cols[k]];
          w.resize(nv, 0.0);
          w[j] += h_rows.vals[k];
        }
      }

      // Cell integrals.
      primal.Get().GetElementDofs(e, pdofs);
      std::vector<mfem::Vector> ue(nt);
      for (int t = 0; t < nt; t++)
      {
        v[t].GetSubVector(pdofs, ue[t]);
      }
      const auto eps = mat_op.GetPermittivityReal(mesh.GetAttribute(e));  // A view.
      eps_inv.Set(1.0, eps);
      eps_inv.Invert();
      M.SetSize(nf);
      M = 0.0;
      d.SetSize(nf);
      d = 0.0;
      g.assign(nv * nt, 0.0);
      r.assign(nv * nt * nf, 0.0);
      vshape.SetSize(nf, 3);
      divshape.SetSize(nf);
      es.SetSize(nf);
      hshape.SetSize(nv);
      hdshape.SetSize(nv, 3);
      pdshape.SetSize(npd, 3);
      mfem::DenseMatrix ev(nf, 3);
      mesh.GetElementTransformation(e, &T);
      const auto &ir =
          mfem::IntRules.Get(fe_f.GetGeomType(), DefaultIntegrationOrder::Get(T));
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        T.SetIntPoint(&ip);
        const double w = ip.weight * T.Weight();
        fe_f.CalcVShape(T, vshape);
        fe_f.CalcPhysDivShape(T, divshape);
        fe_h.CalcShape(ip, hshape);
        fe_h.CalcPhysDShape(T, hdshape);
        fe_p.CalcPhysDShape(T, pdshape);
        mfem::MultABt(vshape, eps_inv, ev);
        mfem::AddMult_a_ABt(w, ev, vshape, M);
        d.Add(w, divshape);
        for (int t = 0; t < nt; t++)
        {
          pdshape.MultTranspose(ue[t], gu);
          eps.Mult(gu, Dv);
          Dv.Neg();  // D_v = -ε ∇v, ε⁻¹ D_v = -∇v.
          vshape.Mult(gu, es);
          for (int j = 0; j < nv; j++)
          {
            double gj = 0.0;
            for (int k = 0; k < 3; k++)
            {
              gj += hdshape(j, k) * Dv(k);
            }
            g[j * nt + t] += w * gj;
            for (int i = 0; i < nf; i++)
            {
              r[(j * nt + t) * nf + i] -= w * hshape(j) * es(i);
            }
          }
        }
      }

      // One record per process owning hats of the cell.
      std::map<int, std::vector<const std::pair<const HYPRE_BigInt, std::vector<double>> *>>
          by_owner;
      for (const auto &entry : cell_hats)
      {
        by_owner[Owner(h_starts, entry.first)].push_back(&entry);
      }
      for (const auto &[p, list] : by_owner)
      {
        auto &buf = send[p];
        buf.push_back(nf);
        buf.insert(buf.end(), face_data.begin(), face_data.end());
        buf.insert(buf.end(), M.Data(), M.Data() + nf * nf);
        buf.insert(buf.end(), d.GetData(), d.GetData() + nf);
        buf.push_back(static_cast<double>(list.size()));
        for (const auto *entry : list)
        {
          const auto &w = entry->second;
          // The faces on which the hat does not vanish (a nonzero weight at a vertex).
          int bits = 0;
          for (int i = 0; i < nf; i++)
          {
            const int fi = face_of_dof[i];
            const int *fv = elem->GetFaceVertices(fi);
            for (int k = 0; k < elem->GetNFaceVertices(fi); k++)
            {
              if (w[fv[k]] != 0.0)
              {
                bits |= (1 << i);
              }
            }
          }
          buf.push_back(static_cast<double>(entry->first));
          buf.push_back(bits);
          for (int t = 0; t < nt; t++)
          {
            double gt = 0.0;
            for (int j = 0; j < nv; j++)
            {
              gt += w[j] * g[j * nt + t];
            }
            buf.push_back(gt);
            for (int i = 0; i < nf; i++)
            {
              double rt = 0.0;
              for (int j = 0; j < nv; j++)
              {
                rt += w[j] * r[(j * nt + t) * nf + i];
              }
              buf.push_back(rt);
            }
          }
        }
      }
    }
  }
  std::vector<int> recv_counts;
  const auto records = mpi::ExchangeAll(comm, send, recv_counts);
  send.assign(np, {});

  // The cells of each owned hat.
  const HYPRE_BigInt h0 = h_starts[rank];
  std::vector<std::vector<PatchCell>> patches(h_starts[rank + 1] - h0);
  for (std::size_t k = 0; k < records.size();)
  {
    PatchCell cell;
    cell.nf = static_cast<int>(records[k]);
    cell.faces = &records[k + 1];
    cell.M = cell.faces + 2 * cell.nf;
    cell.d = cell.M + cell.nf * cell.nf;
    const int nh = static_cast<int>(cell.d[cell.nf]);
    k += 1 + 2 * cell.nf + cell.nf * cell.nf + cell.nf + 1;
    for (int h = 0; h < nh; h++)
    {
      cell.hat = &records[k];
      const auto a = static_cast<HYPRE_BigInt>(records[k]);
      MFEM_ASSERT(a >= h0 && a < h_starts[rank + 1],
                  "Hat record sent to the wrong process!");
      patches[a - h0].push_back(cell);
      k += 2 + nt * (1 + cell.nf);
    }
  }

  // Solve the patch problems. The patch unknowns are the true DOFs of rt0 of the faces of
  // its cells on which the hat does not vanish: their basis functions are supported in the
  // support of the hat (a function of the order-one space which does not vanish on a face
  // does not vanish on the cells next to it, nor on any part of a face which constrains
  // it). The fluxes are summed on the owners of the true DOFs.
  std::vector<Vector> fluxes(nt);
  for (auto &x : fluxes)
  {
    x.SetSize(rt0.GetTrueVSize());
    x.UseDevice(false);
    x = 0.0;
  }
  const HYPRE_BigInt f0 = f_starts[rank];
  for (const auto &cells : patches)
  {
    const int m = static_cast<int>(cells.size());
    if (m == 0)
    {
      continue;
    }
    std::unordered_map<HYPRE_BigInt, int> column;
    std::vector<HYPRE_BigInt> columns;
    for (const auto &cell : cells)
    {
      const int bits = static_cast<int>(cell.hat[1]);
      for (int i = 0; i < cell.nf; i++)
      {
        const auto t = static_cast<HYPRE_BigInt>(cell.faces[2 * i]);
        if (t >= 0 && (bits & (1 << i)) && column.emplace(t, columns.size()).second)
        {
          columns.push_back(t);
        }
      }
    }
    const int n = static_cast<int>(columns.size());
    if (n == 0)
    {
      continue;
    }
    Eigen::MatrixXd Mp = Eigen::MatrixXd::Zero(n, n), B = Eigen::MatrixXd::Zero(m, n),
                    R = Eigen::MatrixXd::Zero(n, nt), G(m, nt);
    std::vector<int> col(8);
    for (int c = 0; c < m; c++)
    {
      const auto &cell = cells[c];
      const int nf = cell.nf;
      col.resize(nf);
      const int bits = static_cast<int>(cell.hat[1]);
      for (int i = 0; i < nf; i++)
      {
        const auto t = static_cast<HYPRE_BigInt>(cell.faces[2 * i]);
        col[i] = (t >= 0 && (bits & (1 << i))) ? column.at(t) : -1;
      }
      for (int i = 0; i < nf; i++)
      {
        if (col[i] < 0)
        {
          continue;
        }
        const double ai = cell.faces[2 * i + 1];
        for (int j = 0; j < nf; j++)
        {
          if (col[j] >= 0)
          {
            Mp(col[i], col[j]) += ai * cell.faces[2 * j + 1] * cell.M[i + nf * j];
          }
        }
        B(c, col[i]) += ai * cell.d[i];
        for (int t = 0; t < nt; t++)
        {
          R(col[i], t) += ai * cell.hat[2 + t * (1 + nf) + 1 + i];
        }
      }
      for (int t = 0; t < nt; t++)
      {
        G(c, t) = cell.hat[2 + t * (1 + nf)];
      }
    }

    // min ½ σᵀ M σ - σᵀ R subject to B σ = G, by the Schur complement S = B M⁻¹ Bᵀ on the
    // cell constraints. Where no unknown carries flux out of the patch (a free hat), the
    // divergences of each basis function sum to zero and so must the data, which they do
    // only up to the algebraic and quadrature error of v: the least-squares solution of the
    // semidefinite system satisfies B σ = G projected onto the range of B, that is with the
    // mean of the data subtracted on each part of the patch joined by unknowns.
    Eigen::LLT<Eigen::MatrixXd> llt(Mp);
    MFEM_VERIFY(llt.info() == Eigen::Success, "Singular vertex patch mass matrix!");
    const Eigen::MatrixXd MiBt = llt.solve(B.transpose()), MiR = llt.solve(R);
    const Eigen::MatrixXd S = B * MiBt;
    const Eigen::MatrixXd lambda = S.completeOrthogonalDecomposition().solve(B * MiR - G);
    const Eigen::MatrixXd sigma = MiR - MiBt * lambda;
    for (int j = 0; j < n; j++)
    {
      const HYPRE_BigInt t = columns[j];
      const int p = Owner(f_starts, t);
      if (p == rank)
      {
        for (int s = 0; s < nt; s++)
        {
          fluxes[s](t - f0) += sigma(j, s);
        }
      }
      else
      {
        auto &buf = send[p];
        buf.push_back(static_cast<double>(t));
        for (int s = 0; s < nt; s++)
        {
          buf.push_back(sigma(j, s));
        }
      }
    }
  }
  const auto sums = mpi::ExchangeAll(comm, send, recv_counts);
  for (std::size_t k = 0; k < sums.size(); k += 1 + nt)
  {
    const auto t = static_cast<HYPRE_BigInt>(sums[k]);
    for (int s = 0; s < nt; s++)
    {
      fluxes[s](t - f0) += sums[k + 1 + s];
    }
  }
  return fluxes;
}

}  // namespace palace::fem
