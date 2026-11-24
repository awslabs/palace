// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_DRIVERS_AMSR_STUDY_HPP
#define PALACE_DRIVERS_AMSR_STUDY_HPP

// STUDY ONLY (ams-robust): instrumentation of the singular magnetostatic solve (gauge /
// gradient component of preconditioner outputs and CG iterates, consistent test right-hand
// side, gauge-handling prototypes). Environment variables:
//   AMSR_RHS=smooth    replace the current excitation by b = K w, w = ND interpolant of a
//                      smooth field (consistent; exact B = Curl w is known)
//   AMSR_PROBE=1       per-component gradient injection of the preconditioner pieces
//   AMSR_MON=k         CG monitor every k iterations (|x|, gradient part, true residual)
//   AMSR_FIX=proj|projamg|final   gauge handling prototypes (see ProjectedPC)
//   AMSR_PROJ_IT=k     AMG V-cycles of the approximate projection (projamg)
// Not for delivery.

#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/ams.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "models/curlcurloperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace::amsr
{

inline int EnvInt(const char *name, int def)
{
  const char *s = std::getenv(name);
  return (s && *s) ? std::atoi(s) : def;
}

inline std::string EnvStr(const char *name)
{
  const char *s = std::getenv(name);
  return s ? s : "";
}

// P_G y = G phi, (G^T M G) phi = G^T M y: the M-orthogonal projection onto the gradients of
// the H1 space with the H1 essential dofs of the ND boundary marker fixed (unit-coefficient
// ND mass M). amg_its > 0: phi from amg_its BoomerAMG V-cycles (a fixed linear operator),
// else PCG + BoomerAMG to rel. tol 1e-12.
class GaugeProjector
{
  MPI_Comm comm;
  std::unique_ptr<mfem::HypreParMatrix> M, L;
  const mfem::HypreParMatrix *G;
  mfem::Array<int> h1_ess;
  std::unique_ptr<mfem::HypreBoomerAMG> amg;
  std::unique_ptr<mfem::HyprePCG> pcg;
  mutable mfem::Vector My, rhs, phi, g, Mg;

public:
  GaugeProjector(FiniteElementSpace &nd, FiniteElementSpace &h1,
                 const mfem::Array<int> &dbc_attr, int amg_its = 0)
    : comm(nd.GetComm())
  {
    {
      mfem::ParBilinearForm m(&nd.Get());
      m.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
      m.Assemble();
      m.Finalize();
      M.reset(m.ParallelAssemble());
    }
    G = &dynamic_cast<const ParOperator &>(nd.GetDiscreteInterpolator(h1)).ParallelAssemble();
    L.reset(mfem::RAP(M.get(), G));
    const auto &pmesh = h1.GetParMesh();
    auto marker = mesh::BdrAttrToMarker(pmesh, dbc_attr, true);
    h1.GetEssentialTrueDofs(marker, h1_ess);
    HYPRE_BigInt ness = h1_ess.Size();
    Mpi::GlobalSum(1, &ness, comm);
    if (ness == 0 && Mpi::Root(comm) && h1.GetTrueVSize() > 0)
    {
      h1_ess.Append(0);  // Natural BCs everywhere: pin one dof (G 1 = 0).
    }
    L->EliminateBC(h1_ess, mfem::Operator::DiagonalPolicy::DIAG_ONE);
    amg = std::make_unique<mfem::HypreBoomerAMG>(*L);
    amg->SetPrintLevel(0);
    if (amg_its > 0)
    {
      amg->SetMaxIter(amg_its);
      amg->SetTol(0.0);
    }
    else
    {
      pcg = std::make_unique<mfem::HyprePCG>(*L);
      pcg->SetTol(1e-12);
      pcg->SetMaxIter(2000);
      pcg->SetPrintLevel(0);
      pcg->SetPreconditioner(*amg);
    }
    My.SetSize(M->Height());
    Mg.SetSize(M->Height());
    g.SetSize(M->Height());
    rhs.SetSize(L->Height());
    phi.SetSize(L->Height());
  }

  const mfem::HypreParMatrix &GetM() const { return *M; }
  const mfem::HypreParMatrix &GetG() const { return *G; }
  const mfem::Array<int> &GetH1Ess() const { return h1_ess; }

  double MNorm(const mfem::Vector &y) const
  {
    M->Mult(y, My);
    return std::sqrt(std::abs(mfem::InnerProduct(comm, y, My)));
  }

  // b <- b - M G L^-1 G^T b: removes the component of a dual vector along M range(G)
  // (G^T b = 0 to the solve tolerance afterwards). Returns |G^T b| before.
  double ProjectDual(mfem::Vector &b) const
  {
    G->MultTranspose(b, rhs);
    for (int i : h1_ess)
    {
      rhs(i) = 0.0;
    }
    const double r0 = std::sqrt(mfem::InnerProduct(comm, rhs, rhs));
    phi = 0.0;
    if (pcg)
    {
      pcg->Mult(rhs, phi);
    }
    else
    {
      amg->Mult(rhs, phi);
    }
    G->Mult(phi, g);
    M->Mult(g, Mg);
    b -= Mg;
    return r0;
  }

  // ny = ||y||_M, ng = ||P_G y||_M; optionally y <- y - P_G y.
  void Measure(mfem::Vector &y, double &ny, double &ng, bool subtract = false) const
  {
    M->Mult(y, My);
    ny = std::sqrt(std::abs(mfem::InnerProduct(comm, y, My)));
    G->MultTranspose(My, rhs);
    for (int i : h1_ess)
    {
      rhs(i) = 0.0;
    }
    phi = 0.0;
    if (pcg)
    {
      pcg->Mult(rhs, phi);
    }
    else
    {
      amg->Mult(rhs, phi);
    }
    G->Mult(phi, g);
    M->Mult(g, Mg);
    ng = std::sqrt(std::abs(mfem::InnerProduct(comm, g, Mg)));
    if (subtract)
    {
      y -= g;
    }
  }
};

// Preconditioner wrapper: z = (I - P_G) B r (the M-orthogonal gradient part of the output
// removed; for consistent residuals (G^T r = 0) this is the symmetric (I - P_G) B (I - P_G)^T).
class ProjectedPC : public Solver<Operator>
{
  const Solver<Operator> &B;
  const GaugeProjector &P;

public:
  ProjectedPC(const Solver<Operator> &B, const GaugeProjector &P) : B(B), P(P) {}
  void SetOperator(const Operator &op) override
  {
    height = op.Height();
    width = op.Width();
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    B.Mult(x, y);
    double ny, ng;
    P.Measure(y, ny, ng, true);
  }
};

// Smooth ND test field w and b = K w (consistent by construction: K G = 0).
inline void SmoothRHS(CurlCurlOperator &op, const Operator &K, Vector &w, Vector &b)
{
  auto &nd = op.GetNDSpace();
  auto &pmesh = const_cast<mfem::ParMesh &>(op.GetMesh().Get());
  mfem::Vector bmin, bmax;
  pmesh.GetBoundingBox(bmin, bmax);
  const int sdim = pmesh.SpaceDimension();
  mfem::Vector c(sdim);
  double diam = 0.0;
  for (int d = 0; d < sdim; d++)
  {
    c(d) = 0.5 * (bmin(d) + bmax(d));
    diam = std::max(diam, bmax(d) - bmin(d));
  }
  const double s = 0.1 * diam;
  auto f = [c, s, sdim](const mfem::Vector &x, mfem::Vector &v)
  {
    mfem::Vector r(x);
    r -= c;
    const double e = std::exp(-(r * r) / (s * s));
    v.SetSize(sdim);
    // A rotation about z plus an x-directed shear: curl is nonzero and smooth.
    v(0) = (-r(1) + 0.3 * r(2)) * e;
    v(1) = (r(0)) * e;
    if (sdim > 2)
    {
      v(2) = (0.5 * r(1) * r(0) / s) * e;
    }
  };
  mfem::VectorFunctionCoefficient fc(sdim, f);
  mfem::ParGridFunction wg(&nd.Get());
  wg.ProjectCoefficient(fc);
  w.SetSize(nd.GetTrueVSize());
  wg.GetTrueDofs(w);
  const auto *Kp = dynamic_cast<const ParOperator *>(
      &dynamic_cast<const MultigridOperator &>(K).GetFinestOperator());
  if (Kp && Kp->GetEssentialTrueDofs())
  {
    for (int i : *Kp->GetEssentialTrueDofs())
    {
      w(i) = 0.0;
    }
  }
  b.SetSize(w.Size());
  K.Mult(w, b);
}

// y = A x for a local CSR matrix (np = 1) with long double accumulation.
class LongDoubleOperator : public Operator
{
  mfem::SparseMatrix A;

public:
  LongDoubleOperator(const mfem::HypreParMatrix &hA) : Operator(hA.Height(), hA.Width())
  {
    MFEM_VERIFY(Mpi::Size(hA.GetComm()) == 1, "LongDoubleOperator is serial only!");
    hA.GetDiag(A);
  }
  void Mult(const Vector &x, Vector &y) const override
  {
    const int *I = A.GetI(), *J = A.GetJ();
    const double *V = A.GetData();
    const double *X = x.HostRead();
    double *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      long double t = 0.0L;
      for (int k = I[i]; k < I[i + 1]; k++)
      {
        t += static_cast<long double>(V[k]) * static_cast<long double>(X[J[k]]);
      }
      Y[i] = static_cast<double>(t);
    }
  }
};

// K through the exact sequence: K = C^T M_RT C (unit coefficient), C the discrete curl
// ND -> RT. Returns the assembled triple product and (factored) the operator applying
// C^T (M (C x)).
struct ExactSequenceK
{
  const mfem::HypreParMatrix *C;
  std::unique_ptr<mfem::HypreParMatrix> M, K;
  class Factored : public Operator
  {
    const mfem::HypreParMatrix &C, &M;
    mutable mfem::Vector t1, t2;

  public:
    Factored(const mfem::HypreParMatrix &C, const mfem::HypreParMatrix &M)
      : Operator(C.Width()), C(C), M(M), t1(C.Height()), t2(C.Height())
    {
    }
    void Mult(const Vector &x, Vector &y) const override
    {
      C.Mult(x, t1);
      M.Mult(t1, t2);
      C.MultTranspose(t2, y);
    }
  };
  std::unique_ptr<Factored> F;
  ExactSequenceK(CurlCurlOperator &op)
  {
    C = &dynamic_cast<const ParOperator &>(op.GetCurlMatrix()).ParallelAssemble();
    mfem::ParBilinearForm m(&op.GetRTSpace().Get());
    m.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
    m.Assemble();
    m.Finalize();
    M.reset(m.ParallelAssemble());
    K.reset(mfem::RAP(M.get(), C));
    F = std::make_unique<Factored>(*C, *M);
  }
};

// ||A G q||_2 / (||A||_2 ||G q||_2) for random q (||A||_2 by 30 power iterations).
inline double Leakage(const Operator &A, const mfem::HypreParMatrix &G,
                      const mfem::Array<int> &h1_ess, MPI_Comm comm, double *nA = nullptr)
{
  Vector x(A.Width()), y(A.Height());
  linalg::SetRandom(comm, x, 5);
  double lam = 0.0;
  for (int it = 0; it < 30; it++)
  {
    x *= 1.0 / linalg::Norml2(comm, x);
    A.Mult(x, y);
    lam = linalg::Norml2(comm, y);
    x = y;
  }
  Vector q(G.Width());
  linalg::SetRandom(comm, q, 7);
  for (int i : h1_ess)
  {
    q(i) = 0.0;
  }
  G.Mult(q, x);
  A.Mult(x, y);
  if (nA)
  {
    *nA = lam;
  }
  return linalg::Norml2(comm, y) / (lam * linalg::Norml2(comm, x));
}

// Per-component gradient injection: for consistent random (K w), consistent smooth
// (restricted RHS) and inconsistent (M G q) inputs of unit l2 norm at each level, the
// M-norm of each component's output and the M-fraction of its gradient part.
inline void Probe(const IoData &iodata, CurlCurlOperator &op, const Operator &K,
                  KspSolver &ksp, const Vector &b)
{
  const auto &linear = iodata.solver.linear;
  const auto &Kmg = dynamic_cast<const MultigridOperator &>(K);
  auto &nd_h = op.GetNDSpaces();
  auto &h1_h = op.GetH1Spaces();
  const int nl = static_cast<int>(nd_h.GetNumLevels());
  MPI_Comm comm = op.GetComm();
  // Restrictions of b to all levels.
  std::vector<Vector> bl(nl);
  bl[nl - 1] = b;
  for (int l = nl - 2; l >= 0; l--)
  {
    bl[l].SetSize(nd_h.GetFESpaceAtLevel(l).GetTrueVSize());
    nd_h.GetProlongationAtLevel(l).MultTranspose(bl[l + 1], bl[l]);
    const auto *Kl = dynamic_cast<const ParOperator *>(&Kmg.GetOperatorAtLevel(l));
    if (Kl->GetEssentialTrueDofs())
    {
      for (int i : *Kl->GetEssentialTrueDofs())
      {
        bl[l](i) = 0.0;
      }
    }
  }
  Mpi::Print(comm, "\nAMSR probe: output M-norm and gradient fraction |P_G y|_M / |y|_M of "
                   "each component, inputs of unit l2 norm\n");
  for (int l = 0; l < nl; l++)
  {
    auto &nd = nd_h.GetFESpaceAtLevel(l);
    auto &h1 = h1_h.GetFESpaceAtLevel(l);
    const auto &Kl = dynamic_cast<const ParOperator &>(Kmg.GetOperatorAtLevel(l));
    GaugeProjector P(nd, h1, op.GetDbcAttr());
    const int n = nd.GetTrueVSize();
    // Inputs.
    Vector w(n), rc(n), rs(bl[l]), ri(n), q(h1.GetTrueVSize()), t(n);
    linalg::SetRandom(comm, w, 11 + l);
    if (Kl.GetEssentialTrueDofs())
    {
      for (int i : *Kl.GetEssentialTrueDofs())
      {
        w(i) = 0.0;
      }
    }
    Kl.Mult(w, rc);
    linalg::SetRandom(comm, q, 23 + l);
    for (int i : P.GetH1Ess())
    {
      q(i) = 0.0;
    }
    P.GetG().Mult(q, t);
    P.GetM().Mult(t, ri);
    std::vector<std::pair<std::string, Vector *>> inputs = {
        {"cons-rand", &rc}, {"cons-rhs", &rs}, {"incons", &ri}};
    for (auto &[name, v] : inputs)
    {
      const double nv = linalg::Norml2(comm, *v);
      if (nv > 0.0)
      {
        *v *= 1.0 / nv;
      }
    }
    // Gradient content of the inputs themselves: |G_0^T r|_2 / |r|_2 (incons input ~ O(1)).
    {
      std::string line;
      Vector gr(h1.GetTrueVSize());
      for (auto &[name, v] : inputs)
      {
        P.GetG().MultTranspose(*v, gr);
        for (int i : P.GetH1Ess())
        {
          gr(i) = 0.0;
        }
        line += fmt::format("  {}: |G^T r|/|r| = {:.3e}", name, linalg::Norml2(comm, gr));
      }
      Mpi::Print(comm, " L{} (n = {}) {:<22s}{}\n", l, nd.GlobalTrueVSize(), "inputs", line);
    }
    {
      double nA;
      const double leak = Leakage(Kl, P.GetG(), P.GetH1Ess(), comm, &nA);
      Mpi::Print(comm, " L{} (n = {}) {:<22s}  ||K G q||/(||K||_2 ||G q||) = {:.3e} (||K||_2 = {:.3e})\n",
                 l, nd.GlobalTrueVSize(), "kernel leakage", leak, nA);
    }
    auto report = [&](const std::string &comp, const Solver<Operator> &C)
    {
      std::string line;
      for (auto &[name, v] : inputs)
      {
        Vector y(n);
        y = 0.0;
        C.Mult(*v, y);
        double ny, ng;
        P.Measure(y, ny, ng);
        line += fmt::format("  {}: |y|_M = {:.3e}, grad = {:.3e}", name, ny,
                            ny > 0.0 ? ng / ny : 0.0);
      }
      Mpi::Print(comm, " L{} (n = {}) {:<22s}{}\n", l, nd.GlobalTrueVSize(), comp, line);
    };
    if (l > 0)
    {
      ChebyshevSmoother<Operator> S(comm, linear.mg_smooth_it, linear.mg_smooth_order,
                                    linear.mg_smooth_sf_max);
      S.SetOperator(Kl);
      report(fmt::format("Chebyshev4(k={})", linear.mg_smooth_order), S);
      ChebyshevSmoother1stKind<Operator> S1(comm, linear.mg_smooth_it,
                                            linear.mg_smooth_order, linear.mg_smooth_sf_max,
                                            linear.mg_smooth_sf_min);
      S1.SetOperator(Kl);
      report(fmt::format("Chebyshev1(k={})", linear.mg_smooth_order), S1);
    }
    else
    {
      for (int singular : {1, 0})
      {
        MfemWrapperSolver<Operator> ams(
            std::make_unique<HypreAmsSolver>(nd, h1, linear.ams_max_it, linear.mg_smooth_it,
                                             linear.ams_vector_interp, singular,
                                             linear.amg_agg_coarsen, 0),
            true, false, false, false);
        ams.SetOperator(Kl);
        report(singular ? "AMS singular" : "AMS full (beta)", ams);
      }
    }
    if (l == nl - 1 && ksp.GetPreconditioner())
    {
      report("full preconditioner", *ksp.GetPreconditioner());
    }
  }
}


// CG monitor: every k iterations, |x|_2, |x|_M, gradient fraction, true and recursive
// residuals, and the B-field error if the exact field is known.
inline void SetMonitor(int every, const Operator &K, const Vector &b, const GaugeProjector &P,
                       const Operator *Curl, const Vector *w, MPI_Comm comm,
                       const Solver<Operator> *B = nullptr)
{
  auto Bw = std::make_shared<Vector>();
  if (Curl && w)
  {
    Bw->SetSize(Curl->Height());
    Curl->Mult(*w, *Bw);
  }
  const double nb = linalg::Norml2(comm, b);
  CgSolver<Operator>::study_monitor =
      [every, &K, &b, &P, Curl, w, comm, Bw, nb, B](int it, const Vector &x, const Vector &p,
                                                    double alpha, const Vector &r)
  {
    if (it % every)
    {
      return;
    }
    Vector xc(x), t(x.Size());
    xc.Add(alpha, p);
    K.Mult(xc, t);
    t -= b;
    const double tr = linalg::Norml2(comm, t) / nb, rr = linalg::Norml2(comm, r) / nb;
    const double nx2 = linalg::Norml2(comm, xc);
    double eB = -1.0;
    if (Curl && w)
    {
      Vector Bx(Curl->Height());
      Curl->Mult(xc, Bx);
      Bx -= *Bw;
      eB = linalg::Norml2(comm, Bx) / linalg::Norml2(comm, *Bw);
    }
    double nx, ng;
    P.Measure(xc, nx, ng);
    std::string zs;
    if (B)
    {
      Vector z(r.Size());
      z = 0.0;
      B->Mult(r, z);
      double nz, nzg;
      P.Measure(z, nz, nzg);
      const double nzr = std::sqrt(std::max(nz * nz - nzg * nzg, 0.0));
      zs = fmt::format(", Br: |z|_M = {:.3e}, grad/range = {:.3e}", nz, nzr > 0.0 ? nzg / nzr : 0.0);
    }
    Mpi::Print(comm,
               "  AMSR mon it {:5d}: |x|_2 = {:.3e}, |x|_M = {:.3e}, grad = {:.3e}, "
               "true res = {:.3e}, rec res = {:.3e}, eB = {:.3e}{}\n",
               it, nx2, nx, nx > 0.0 ? ng / nx : 0.0, tr, rr, eB, zs);
  };
}

}  // namespace palace::amsr

#endif  // PALACE_DRIVERS_AMSR_STUDY_HPP
