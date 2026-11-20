// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>
#include <Eigen/SparseLU>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/vertexpatchflux.hpp"
#include "linalg/equilibratedestimator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"

namespace palace
{

namespace
{

// div(ε ∇u) = 0 in a box of two materials, ε = diag(2, 3, 5) (attribute 1) and 7 times that
// (attribute 2, the cells above the middle of the second direction), u = 1 on a terminal
// (boundary attribute 7, or 9 for a second terminal), u = 0 on the ground (8), natural
// conditions elsewhere (1).
constexpr int terminal = 7, ground = 8, second = 9;
const double eps_diag[3] = {2.0, 3.0, 5.0};

// The quadrature of a run of order p with Palace's default settings
// (QuadratureOrderJacobian false, no extra order), as IoData::CheckConfiguration sets it,
// restored at the end of the scope: the globals keep the values of earlier tests otherwise.
class RunQuadrature
{
private:
  const int p_trial = fem::DefaultIntegrationOrder::p_trial,
            extra_pk = fem::DefaultIntegrationOrder::q_order_extra_pk,
            extra_qk = fem::DefaultIntegrationOrder::q_order_extra_qk;
  const bool jac = fem::DefaultIntegrationOrder::q_order_jac;

public:
  explicit RunQuadrature(int p)
  {
    fem::DefaultIntegrationOrder::p_trial = p;
    fem::DefaultIntegrationOrder::q_order_jac = false;
    fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
    fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  }
  ~RunQuadrature()
  {
    fem::DefaultIntegrationOrder::p_trial = p_trial;
    fem::DefaultIntegrationOrder::q_order_jac = jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = extra_pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = extra_qk;
  }
};

double Eps(int attr, int d)
{
  return eps_diag[d] * ((attr == 2) ? 7.0 : 1.0);
}

std::vector<config::MaterialData> Materials()
{
  std::vector<config::MaterialData> materials(2);
  for (int attr : {1, 2})
  {
    auto &data = materials[attr - 1];
    data.attributes = {attr};
    data.epsilon_r.s = {Eps(attr, 0), Eps(attr, 1), Eps(attr, 2)};
    data.epsilon_r.v = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  }
  return materials;
}

// The coefficient ε of MFEM's diffusion integrator.
struct EpsCoefficient
{
  mfem::VectorConstantCoefficient eps1{mfem::Vector({Eps(1, 0), Eps(1, 1), Eps(1, 2)})},
      eps2{mfem::Vector({Eps(2, 0), Eps(2, 1), Eps(2, 2)})};
  mfem::Array<int> attrs{1, 2};
  mfem::Array<mfem::VectorCoefficient *> coeffs{&eps1, &eps2};
  mfem::PWVectorCoefficient pw{3, attrs, coeffs};
};

using Classifier = std::function<int(const std::array<double, 3> &, bool)>;

// The elements of base, with a boundary element of attribute attribute(center, exterior) >
// 0 on each such face.
mfem::Mesh WithBoundary(const mfem::Mesh &base, const Classifier &attribute)
{
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(), 0, 3);
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    mfem::Array<int> vertices;
    base.GetFaceVertices(f, vertices);
    std::array<double, 3> center{};
    for (int v : vertices)
    {
      for (int d = 0; d < 3; d++)
      {
        center[d] += base.GetVertex(v)[d] / vertices.Size();
      }
    }
    const int attr = attribute(center, base.GetFaceInformation(f).IsBoundary());
    if (attr > 0)
    {
      auto *face = base.GetFace(f)->Duplicate(&mesh);
      face->SetAttribute(attr);
      mesh.AddBdrElement(face);
    }
  }
  mesh.FinalizeTopology();
  mesh.Finalize(false, true);
  mesh.SetAttributes();
  return mesh;
}

// The test geometries: a box with a terminal on part of one side and the ground on the
// opposite side (BOX), the same with an interior ground sheet across the middle instead
// (SHEET) or with a second, disconnected part of the ground on a strip of a third side
// (SPLIT), a thin slab with a terminal strip below and the ground above (SLAB), and a box
// with two terminals (TWO).
struct Geometry
{
  mfem::Element::Type type;
  std::array<int, 3> n;
  std::array<double, 3> s;
  enum Kind
  {
    BOX,
    SHEET,
    SPLIT,
    SLAB,
    TWO
  } kind;

  int Attribute(const std::array<double, 3> &x, bool exterior) const
  {
    const double tol = 1.0e-9 * std::min(s[0], s[2]);
    switch (kind)
    {
      case BOX:
      case SHEET:
        if (!exterior)
        {
          return (kind == SHEET && std::abs(x[0] - 0.5 * s[0]) < tol) ? ground : 0;
        }
        if (x[0] < tol && x[1] < 0.5 * s[1])
        {
          return terminal;
        }
        return (kind == BOX && x[0] > s[0] - tol) ? ground : 1;
      case SPLIT:
        if (!exterior)
        {
          return 0;
        }
        if (x[0] < tol && x[1] < 0.5 * s[1])
        {
          return terminal;
        }
        if (x[0] > s[0] - tol || (x[1] < tol && x[0] > 0.25 * s[0] && x[0] < 0.75 * s[0]))
        {
          return ground;
        }
        return 1;
      case TWO:
        if (!exterior)
        {
          return 0;
        }
        if (x[0] < tol && x[1] < 0.5 * s[1])
        {
          return terminal;
        }
        if (x[0] < tol && x[1] > 0.75 * s[1])
        {
          return second;
        }
        return (x[0] > s[0] - tol) ? ground : 1;
      case SLAB:
        if (!exterior)
        {
          return 0;
        }
        if (x[2] < tol)
        {
          return (x[0] < 0.5 * s[0]) ? terminal : 1;
        }
        return (x[2] > s[2] - tol) ? ground : 1;
    }
    return 0;
  }

  mfem::Mesh Make() const
  {
    auto base = mfem::Mesh::MakeCartesian3D(n[0], n[1], n[2], type, s[0], s[1], s[2]);
    mfem::Vector c(3);
    for (int e = 0; e < base.GetNE(); e++)
    {
      base.GetElementCenter(e, c);
      base.SetAttribute(e, (c(1) > 0.5 * s[1]) ? 2 : 1);
    }
    base.SetAttributes();
    return WithBoundary(base, [this](const std::array<double, 3> &x, bool exterior)
                        { return Attribute(x, exterior); });
  }

  std::vector<int> Terminals() const
  {
    return (kind == TWO) ? std::vector<int>{terminal, second} : std::vector<int>{terminal};
  }

  std::vector<int> Dirichlet() const
  {
    return (kind == TWO) ? std::vector<int>{terminal, second, ground}
                         : std::vector<int>{terminal, ground};
  }
};

// Nonconforming refinement of the cells centered in x0 < x / s < x0 + 1/4 (in a thin box
// only those of the lower layer). The same marking on the serial and on the distributed
// mesh gives the same refined mesh.
template <typename MeshType>
void Refine(MeshType &mesh, double x0, const std::array<double, 3> &hi)
{
  mfem::Vector c(3);
  const bool thin = hi[2] < 0.5 * hi[0];
  mfem::Array<mfem::Refinement> marked;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    if (c(0) > x0 * hi[0] && c(0) < (x0 + 0.25) * hi[0] && (!thin || c(2) < 0.5 * hi[2]))
    {
      marked.Append(mfem::Refinement(e));
    }
  }
  mesh.GeneralRefinement(marked, 1, 0);
}

// A sparse direct solve (Eigen SparseLU) of an MFEM serial matrix.
mfem::Vector DirectSolve(const mfem::SparseMatrix &A, const mfem::Vector &b)
{
  std::vector<Eigen::Triplet<double>> triplets;
  for (int i = 0; i < A.Height(); i++)
  {
    for (int k = A.GetI()[i]; k < A.GetI()[i + 1]; k++)
    {
      triplets.emplace_back(i, A.GetJ()[k], A.GetData()[k]);
    }
  }
  Eigen::SparseMatrix<double> K(A.Height(), A.Width());
  K.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(K);
  REQUIRE(lu.info() == Eigen::Success);
  Eigen::VectorXd r(b.Size());
  for (int i = 0; i < b.Size(); i++)
  {
    r(i) = b(i);
  }
  const Eigen::VectorXd x = lu.solve(r);
  mfem::Vector y(b.Size());
  for (int i = 0; i < b.Size(); i++)
  {
    y(i) = x(i);
  }
  return y;
}

// The potential of the terminal on a serial mesh, at order p (MFEM assembly, direct solve),
// and its energy (ε ∇v, ∇v) = C_h.
struct SerialSolve
{
  std::unique_ptr<mfem::H1_FECollection> fec;
  std::unique_ptr<mfem::FiniteElementSpace> fes;
  std::unique_ptr<mfem::GridFunction> v;
  double energy2 = 0.0;
};

SerialSolve SolveSerial(mfem::Mesh &mesh, const std::vector<int> &dirichlet, int p,
                        int excited)
{
  SerialSolve s;
  s.fec = std::make_unique<mfem::H1_FECollection>(p, 3);
  s.fes = std::make_unique<mfem::FiniteElementSpace>(&mesh, s.fec.get());
  EpsCoefficient coefficient;
  mfem::BilinearForm a(s.fes.get());
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(coefficient.pw));
  a.Assemble();
  a.Finalize();
  const int max = mesh.bdr_attributes.Max();
  mfem::Array<int> essential;
  s.fes->GetEssentialTrueDofs(mesh::AttrToMarker(max, dirichlet), essential);
  s.v = std::make_unique<mfem::GridFunction>(s.fes.get());
  *s.v = 0.0;
  mfem::ConstantCoefficient one(1.0);
  s.v->ProjectBdrCoefficient(one, mesh::AttrToMarker(max, std::vector<int>{excited}));
  mfem::LinearForm b(s.fes.get());
  b = 0.0;
  mfem::SparseMatrix A;
  mfem::Vector X, B;
  a.FormLinearSystem(essential, *s.v, b, A, X, B);
  const mfem::Vector solution = DirectSolve(A, B);
  X = solution;  // A copy: X may reference the data of v.
  a.RecoverFEMSolution(X, b, *s.v);
  mfem::Vector grad(3);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    auto &T = *mesh.GetElementTransformation(e);
    const auto &ir = mfem::IntRules.Get(mesh.GetElementGeometry(e), 2 * p + 2);
    for (int k = 0; k < ir.GetNPoints(); k++)
    {
      T.SetIntPoint(&ir.IntPoint(k));
      s.v->GetGradient(T, grad);
      for (int d = 0; d < 3; d++)
      {
        s.energy2 += ir.IntPoint(k).weight * T.Weight() * Eps(mesh.GetAttribute(e), d) *
                     grad(d) * grad(d);
      }
    }
  }
  return s;
}

// The reference: the dual mixed problem over the RT(q) fluxes of the cells of the serial
// mesh (normal continuous across every interior face except the Dirichlet sheets, on a
// nonconforming face the trace of the coarse cell on each fine face, imposed against the
// face polynomials, which is exact), divergence free (tested against L2(q) of each cell),
// with zero normal trace on the Neumann faces and free on the Dirichlet faces and on both
// sides of a sheet, for the potential v of order p:
//   η² = min ‖ε^{-1/2} (D + ε ∇v)‖²,
// one KKT system, solved directly. The best equilibrated flux of RT(q), independent of the
// patches, the broken spaces and the curl correction.
double MixedReference(mfem::Mesh &mesh, const Geometry &g, int p, int excited, int q)
{
  const auto s = SolveSerial(mesh, g.Dirichlet(), p, excited);
  const auto dirichlet = g.Dirichlet();
  mfem::RT_FECollection rt_fec(q, 3);
  mfem::L2_FECollection l2_fec(q, 3), trace_fec(q, 2);
  const int ne = mesh.GetNE();
  std::vector<int> offset(ne + 1, 0);
  for (int e = 0; e < ne; e++)
  {
    offset[e + 1] =
        offset[e] + rt_fec.FiniteElementForGeometry(mesh.GetElementGeometry(e))->GetDof();
  }
  const int n = offset[ne];
  std::vector<Eigen::Triplet<double>> K;  // KKT, the constraint rows from n.
  Eigen::VectorXd f = Eigen::VectorXd::Zero(n);
  int rows = n;
  mfem::DenseMatrix shape, eps_inv(3);
  mfem::Vector grad(3), div, l2;
  for (int e = 0; e < ne; e++)
  {
    eps_inv = 0.0;
    for (int d = 0; d < 3; d++)
    {
      eps_inv(d, d) = 1.0 / Eps(mesh.GetAttribute(e), d);
    }
    const auto &fe = *rt_fec.FiniteElementForGeometry(mesh.GetElementGeometry(e));
    const auto &le = *l2_fec.FiniteElementForGeometry(mesh.GetElementGeometry(e));
    const int nd = fe.GetDof();
    auto &T = *mesh.GetElementTransformation(e);
    const auto &ir = mfem::IntRules.Get(fe.GetGeomType(), 2 * std::max(p, q) + 4);
    mfem::DenseMatrix M(nd), es_shape(nd, 3);
    M = 0.0;
    shape.SetSize(nd, 3);
    for (int k = 0; k < ir.GetNPoints(); k++)
    {
      const auto &ip = ir.IntPoint(k);
      T.SetIntPoint(&ip);
      const double w = ip.weight * T.Weight();
      fe.CalcVShape(T, shape);
      mfem::MultABt(shape, eps_inv, es_shape);
      mfem::AddMult_a_ABt(w, es_shape, shape, M);
      s.v->GetGradient(T, grad);
      for (int j = 0; j < nd; j++)
      {
        double dot = 0.0;
        for (int d = 0; d < 3; d++)
        {
          dot += grad(d) * shape(j, d);
        }
        f(offset[e] + j) -= w * dot;
      }
    }
    for (int i = 0; i < nd; i++)
    {
      for (int j = 0; j < nd; j++)
      {
        K.emplace_back(offset[e] + i, offset[e] + j, M(i, j));
      }
    }
    // Divergence rows: ∫ l div σ dx = sign(det J) ∫ l div̂ σ̂ dx̂.
    const int nl = le.GetDof();
    div.SetSize(nd);
    l2.SetSize(nl);
    mfem::DenseMatrix D(nl, nd);
    D = 0.0;
    for (int k = 0; k < ir.GetNPoints(); k++)
    {
      const auto &ip = ir.IntPoint(k);
      fe.CalcDivShape(ip, div);
      le.CalcShape(ip, l2);
      mfem::AddMult_a_VWt(ip.weight, l2, div, D);
    }
    for (int i = 0; i < nl; i++, rows++)
    {
      for (int j = 0; j < nd; j++)
      {
        K.emplace_back(rows, offset[e] + j, D(i, j));
        K.emplace_back(offset[e] + j, rows, D(i, j));
      }
    }
  }
  mfem::Vector normal(3), mu;
  mfem::DenseMatrix shape1, shape2;
  for (int fi = 0; fi < mesh.GetNumFaces(); fi++)
  {
    const auto info = mesh.GetFaceInformation(fi);
    if (info.IsNonconformingCoarse())
    {
      continue;  // Its fine faces carry the constraints.
    }
    auto *FT = mesh.GetFaceElementTransformations(fi);
    const auto &tfe = *trace_fec.FiniteElementForGeometry(FT->GetGeometryType());
    mfem::Vector c(3);
    {
      mfem::IntegrationPoint center;
      center.Set3(mfem::Geometries.GetCenter(FT->GetGeometryType()).x,
                  mfem::Geometries.GetCenter(FT->GetGeometryType()).y, 0.0);
      FT->Transform(center, c);
    }
    const bool exterior = info.IsBoundary();
    const int attr = g.Attribute({c(0), c(1), c(2)}, exterior);
    if (std::find(dirichlet.begin(), dirichlet.end(), attr) != dirichlet.end())
    {
      continue;  // Free normal trace (exterior), or no continuity (sheet).
    }
    const int e1 = FT->Elem1No, e2 = exterior ? -1 : FT->Elem2No;
    const auto &fe1 = *rt_fec.FiniteElementForGeometry(mesh.GetElementGeometry(e1));
    shape1.SetSize(fe1.GetDof(), 3);
    const mfem::FiniteElement *fe2 =
        (e2 >= 0) ? rt_fec.FiniteElementForGeometry(mesh.GetElementGeometry(e2)) : nullptr;
    if (fe2)
    {
      shape2.SetSize(fe2->GetDof(), 3);
    }
    const int nt = tfe.GetDof();
    mu.SetSize(nt);
    std::vector<double> row1(nt * fe1.GetDof(), 0.0),
        row2(fe2 ? nt * fe2->GetDof() : 0, 0.0);
    const auto &ir = mfem::IntRules.Get(FT->GetGeometryType(), 2 * std::max(p, q) + 4);
    for (int k = 0; k < ir.GetNPoints(); k++)
    {
      const auto &ip = ir.IntPoint(k);
      FT->SetAllIntPoints(&ip);
      mfem::CalcOrtho(FT->Jacobian(), normal);  // Out of e1, scaled by the area element.
      tfe.CalcShape(ip, mu);
      fe1.CalcVShape(*FT->Elem1, shape1);
      for (int j = 0; j < fe1.GetDof(); j++)
      {
        double sn = 0.0;
        for (int d = 0; d < 3; d++)
        {
          sn += shape1(j, d) * normal(d);
        }
        for (int i = 0; i < nt; i++)
        {
          row1[i * fe1.GetDof() + j] += ip.weight * mu(i) * sn;
        }
      }
      if (fe2)
      {
        fe2->CalcVShape(*FT->Elem2, shape2);
        for (int j = 0; j < fe2->GetDof(); j++)
        {
          double sn = 0.0;
          for (int d = 0; d < 3; d++)
          {
            sn += shape2(j, d) * normal(d);
          }
          for (int i = 0; i < nt; i++)
          {
            row2[i * fe2->GetDof() + j] -= ip.weight * mu(i) * sn;
          }
        }
      }
    }
    for (int i = 0; i < nt; i++, rows++)
    {
      for (int j = 0; j < fe1.GetDof(); j++)
      {
        K.emplace_back(rows, offset[e1] + j, row1[i * fe1.GetDof() + j]);
        K.emplace_back(offset[e1] + j, rows, row1[i * fe1.GetDof() + j]);
      }
      for (int j = 0; fe2 && j < fe2->GetDof(); j++)
      {
        K.emplace_back(rows, offset[e2] + j, row2[i * fe2->GetDof() + j]);
        K.emplace_back(offset[e2] + j, rows, row2[i * fe2->GetDof() + j]);
      }
    }
  }
  Eigen::SparseMatrix<double> KKT(rows, rows);
  KKT.setFromTriplets(K.begin(), K.end());
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(rows);
  rhs.head(n) = f;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(KKT);
  REQUIRE(lu.info() == Eigen::Success);
  const Eigen::VectorXd x = lu.solve(rhs);
  const Eigen::VectorXd sigma = x.head(n);
  const Eigen::SparseMatrix<double> M = KKT.topLeftCorner(n, n);
  CHECK((KKT * x - rhs).norm() < 1.0e-10 * rhs.norm());
  return sigma.dot(M * sigma) - 2.0 * f.dot(sigma) + s.energy2;
}

// The distributed solve of a terminal (CG 1e-14 with AMG), its intended essential values
// and its energy U_h.
struct Solve
{
  Vector V, G;
  double energy = 0.0;
};

Solve SolveTerminal(FiniteElementSpace &primal, const std::vector<int> &dirichlet,
                    int excited)
{
  auto &fes = primal.Get();
  EpsCoefficient coefficient;
  mfem::ParBilinearForm a(&fes);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(coefficient.pw));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> K(a.ParallelAssemble());
  const int max = primal.GetParMesh().bdr_attributes.Max();
  mfem::Array<int> essential;
  fes.GetEssentialTrueDofs(mesh::AttrToMarker(max, dirichlet), essential);
  mfem::ParGridFunction x(&fes);
  x = 0.0;
  mfem::ConstantCoefficient one(1.0);
  x.ProjectBdrCoefficient(one, mesh::AttrToMarker(max, std::vector<int>{excited}));
  Solve s;
  s.G.SetSize(fes.GetTrueVSize());
  x.GetTrueDofs(s.G);
  mfem::ParLinearForm b(&fes);
  b = 0.0;
  mfem::HypreParMatrix A;
  mfem::Vector X, B;
  a.FormLinearSystem(essential, x, b, A, X, B);
  mfem::HypreBoomerAMG amg(A);
  amg.SetPrintLevel(0);
  mfem::CGSolver cg(fes.GetComm());
  cg.SetRelTol(1.0e-14);
  cg.SetMaxIter(4000);
  cg.SetOperator(A);
  cg.SetPreconditioner(amg);
  cg.Mult(B, X);
  a.RecoverFEMSolution(X, b, x);
  s.V.SetSize(fes.GetTrueVSize());
  x.GetTrueDofs(s.V);
  Vector KV(s.V.Size());
  K->Mult(s.V, KV);
  s.energy = 0.5 * mfem::InnerProduct(fes.GetComm(), KV, s.V);
  return s;
}

// Resolved linear solver settings, with the curl correction run to convergence.
config::LinearSolverData Linear()
{
  config::LinearSolverData linear;
  linear.type = LinearSolver::BOOMER_AMG;  // As resolved for electrostatics.
  linear.krylov_solver = KrylovSolver::CG;
  linear.mg_max_levels = 100;
  linear.mg_cycle_it = 1;
  linear.mg_smooth_order = 4;
  linear.amg_agg_coarsen = 1;
  linear.ams_max_it = 2;
  linear.estimator_tol = 1.0e-14;
  linear.estimator_error_tol = 1.0e-10;
  linear.estimator_max_it = 4000;
  return linear;
}

struct Case
{
  Geometry g;
  bool refine;
  double x0;
  int p;
  // A sheet refined on one side keeps its constraint to the coarse face in the broken
  // spaces, which the reference does not have: the bound is admissible, not optimal.
  bool exact = true;
};

std::vector<Case> Cases()
{
  using G = Geometry;
  return {
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, false, 0.0, 1},
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, false, 0.0, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, true, 0.0, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, false, 0.0, 3},
      {G{mfem::Element::HEXAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, true, 0.0, 1},
      {G{mfem::Element::HEXAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::BOX}, false, 0.0, 2},
      {G{mfem::Element::WEDGE, {4, 2, 2}, {1.0, 1.0, 1.0}, G::SHEET}, false, 0.0, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::SHEET}, true, 0.0, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 4, 2}, {1.0, 1.0, 0.025}, G::SLAB}, true, 0.25, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 4, 2}, {1.0, 1.0, 0.025}, G::SLAB}, false, 0.0, 1},
      {G{mfem::Element::TETRAHEDRON, {4, 4, 2}, {1.0, 1.0, 1.0}, G::TWO}, false, 0.0, 2},
      {G{mfem::Element::HEXAHEDRON, {4, 4, 2}, {1.0, 1.0, 1.0}, G::TWO}, true, 0.0, 2},
      {G{mfem::Element::TETRAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::SPLIT}, false, 0.0, 2},
      {G{mfem::Element::HEXAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::SPLIT}, true, 0.0, 1},
      {G{mfem::Element::HEXAHEDRON, {4, 2, 2}, {1.0, 1.0, 1.0}, G::SHEET}, true, 0.25, 2,
       false}};
}

// The serial and distributed meshes of a case, refined alike.
std::pair<mfem::Mesh, std::unique_ptr<Mesh>> MakeMeshes(const Case &c)
{
  auto serial = c.g.Make(), reference = c.g.Make();
  serial.EnsureNCMesh(true);
  reference.EnsureNCMesh(true);
  if (c.refine)
  {
    Refine(reference, c.x0, c.g.s);
  }
  auto mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(Mpi::World(), serial));
  if (c.refine)
  {
    Refine(mesh->Get(), c.x0, c.g.s);
    mesh->Update();
  }
  return {std::move(reference), std::move(mesh)};
}

// The checks of the dual mixed test below for a case with the flux space RT(p +
// flux_order_offset). C_fine holds the capacitances of the fine reference of the terminals,
// which do not depend on the flux space: solved for the first flux order of the case,
// reused for the second.
void CheckDualMixed(const Case &c, int flux_order_offset, std::vector<double> &C_fine)
{
  auto [reference, mesh] = MakeMeshes(c);
  const RunQuadrature run(c.p);
  mfem::H1_FECollection h1_fec(c.p, 3);
  mfem::ND_FECollection nd_fec(c.p, 3);
  mfem::RT_FECollection rt_fec(c.p - 1, 3);
  FiniteElementSpace primal(*mesh, &h1_fec), nd(*mesh, &nd_fec), rt(*mesh, &rt_fec);
  const auto dirichlet = c.g.Dirichlet(), terminals = c.g.Terminals();
  const int n = static_cast<int>(terminals.size());
  std::vector<Vector> V, G;
  std::vector<double> energies;
  for (int t : terminals)
  {
    auto s = SolveTerminal(primal, dirichlet, t);
    V.push_back(s.V);
    G.push_back(s.G);
    energies.push_back(s.energy);
  }
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(Materials(), periodic, ProblemType::ELECTROSTATIC, *mesh);
  mfem::Array<int> dirichlet_attr;
  for (int attr : dirichlet)
  {
    dirichlet_attr.Append(attr);
  }
  EquilibratedFluxErrorEstimator estimator(mat_op, primal, rt, nd, dirichlet_attr,
                                           flux_order_offset, Linear());
  ErrorIndicator indicator;
  const auto summaries = estimator.AddErrorIndicators(V, G, energies, indicator);
  const auto &C_h = estimator.GetCapacitance(),
             &C_low = estimator.GetCapacitanceLowerBound();

  // A fine reference (order p + 1 on the uniformly refined mesh) for the true error:
  // its capacitance is above the exact one, so C_h - C_fine is below the true error².
  if (C_fine.empty())
  {
    mfem::Mesh fine(reference);
    fine.UniformRefinement();
    for (int t = 0; t < n; t++)
    {
      C_fine.push_back(SolveSerial(fine, dirichlet, c.p + 1, terminals[t]).energy2);
    }
  }
  double rms2 = 0.0;
  std::vector<double> refs(n);
  for (int t = 0; t < n; t++)
  {
    const auto &s = summaries[t];
    const double ref = refs[t] = std::sqrt(
        MixedReference(reference, c.g, c.p, terminals[t], c.p + flux_order_offset));
    INFO("kind " << c.g.kind << " type " << c.g.type << " refine " << c.refine << " p "
                 << c.p << " flux RT(p" << flux_order_offset << ") terminal " << t
                 << ": mixed " << ref << ", vertex patches " << s.initial << " -> " << s.eta
                 << " (" << s.iterations << " its)");
    CHECK(s.divergence <= 1.0e-10 * std::sqrt(C_h(t, t)));
    CHECK(s.neumann <= 1.0e-10 * std::sqrt(C_h(t, t)));
    CHECK(s.eta <= s.initial * (1.0 + 1.0e-12));
    CHECK(std::abs(C_h(t, t) - 2.0 * energies[t]) <= 1.0e-10 * C_h(t, t));
    CHECK(C_low(t, t) <= C_h(t, t));
    CHECK(s.bound * s.bound >= (C_h(t, t) - C_fine[t]) * (1.0 - 1.0e-8));
    CHECK(C_low(t, t) <= C_fine[t] * (1.0 + 1.0e-10));
    CHECK(s.relative >= s.bound / std::sqrt(C_h(t, t)));
    rms2 += s.relative * s.relative / n;
    if (!c.exact)
    {
      CHECK(s.eta >= ref * (1.0 - 1.0e-7));
      continue;
    }
    CHECK(std::abs(s.eta - ref) <= 1.0e-7 * ref);
    const double C_star = C_h(t, t) - ref * ref;
    CHECK(std::abs(C_low(t, t) - C_star) <= 1.0e-7 * C_star);
  }
  // The indicator norm is the root mean square of the relative bounds of the solves.
  CHECK(std::abs(indicator.Norml2(Mpi::World()) - std::sqrt(rms2)) <=
        1.0e-10 * std::sqrt(rms2));

  // With the default EstimatorErrorTol the fields of the potentials are corrected only as
  // far as their weight requires: each bound is above the mixed one and within the
  // tolerances of the corrections (a flux error of at most tol' η from the solve's own
  // correction and tol' η from the fields, tol' = 0.14 for tol = 0.01).
  if (c.exact && estimator.GetNumHarmonicFields() > 0)
  {
    auto linear = Linear();
    linear.estimator_error_tol = 1.0e-2;
    EquilibratedFluxErrorEstimator loose(mat_op, primal, rt, nd, dirichlet_attr,
                                         flux_order_offset, linear);
    ErrorIndicator loose_indicator;
    const auto loose_summaries = loose.AddErrorIndicators(V, G, energies, loose_indicator);
    const auto &loose_C_low = loose.GetCapacitanceLowerBound();
    for (int t = 0; t < n; t++)
    {
      const double C_star = C_h(t, t) - refs[t] * refs[t];
      INFO("kind " << c.g.kind << " type " << c.g.type << " p " << c.p << " terminal " << t
                   << ": " << loose.GetNumCorrectedHarmonicFields() << " of "
                   << loose.GetNumHarmonicFields() << " fields corrected, eta "
                   << loose_summaries[t].eta << " against the mixed " << refs[t]);
      CHECK(loose_summaries[t].eta >= refs[t] * (1.0 - 1.0e-7));
      CHECK(loose_summaries[t].eta <= refs[t] * 1.05);
      CHECK(loose_C_low(t, t) <= C_star * (1.0 + 1.0e-7));
    }
  }
}

}  // namespace

// The converged flux, after the capacitive step, is the dual mixed solution of its space:
// the bound equals the one of the direct mixed RT(k) x L2(k) solve, k = p - 1 (the space of
// L2ZZ) or p, and C_low equals the mixed capacitance C* = C_h - η*². Also: the bound is
// above the true error (C_low below the capacitance of a much finer solve), the flux is
// equilibrated, and the indicators sum to the relative bound, which is above the bound
// relative to the computed energy. Tets, hexes and wedges, p 1 to 3, conforming and
// nonconforming, a thin slab, an interior Dirichlet sheet (broken spaces), two terminals.
// A representative subset of the cases: every element type, geometry, polynomial order and
// flux space, the nonconforming refinements of a box (tets and hexes) and of an interior
// boundary, the case where the bound is admissible but not optimal, and hexes with RT(p)
// (whose curl-curl terms of ND(p + 1) need more points than the run's quadrature), for a
// third of the cost of all of them, which the [Regression] case below checks.
TEST_CASE("Equilibrated estimator reaches the dual mixed solution",
          "[equilibratedestimator][Serial][Parallel]")
{
  // Indices into Cases(), with their flux order offsets.
  const std::vector<std::pair<std::size_t, std::vector<int>>> subset = {
      {0, {-1, 0}},  // Tets, box, p 1: RT(0) and RT(1).
      {3, {-1}},     // Tets, box, p 3.
      {4, {0}},      // Nonconforming hexes, box, p 1.
      {5, {0}},      // Hexes, box, p 2, RT(2) ("Primal"): ND(3) on hexes needs more points.
      {6, {-1, 0}},  // Wedges, interior boundary, p 2: RT(1) and RT(2).
      {7, {-1}},     // Nonconforming tets, interior boundary, p 2.
      {9, {-1}},     // Thin slab.
      {10, {-1}},    // Two terminals.
      {13, {-1}},    // Nonconforming hexes, split ground.
      {14, {-1}}};   // Hexes refined on one side of the interior boundary.
  const auto cases = Cases();
  for (const auto &[i, offsets] : subset)
  {
    std::vector<double> C_fine;
    for (const int flux_order_offset : offsets)
    {
      CheckDualMixed(cases[i], flux_order_offset, C_fine);
    }
  }
}

// All the cases, with RT(p - 1) and, but at p 3, RT(p).
TEST_CASE("Equilibrated estimator reaches the dual mixed solution (all cases)",
          "[equilibratedestimator][Serial][Parallel][Regression]")
{
  for (const auto &c : Cases())
  {
    std::vector<double> C_fine;
    for (const int flux_order_offset : {-1, 0})
    {
      if (flux_order_offset == 0 && c.p > 2)
      {
        continue;
      }
      CheckDualMixed(c, flux_order_offset, C_fine);
    }
  }
}

// The vertex patch flux is exactly equilibrated (divergence free, zero normal trace on the
// Neumann faces, at rounding) and, where no patch touches two Dirichlet boundaries (the
// boxes), carries the fluxes of the computed solution into each of them: q_k(D₀) =
// (D₀, ∇w_k) = (D_v, ∇w_k) for the order-one function w_k which is one on the vertices of
// boundary k.
TEST_CASE("Vertex patch flux is equilibrated", "[equilibratedestimator][Serial][Parallel]")
{
  for (const auto &c : Cases())
  {
    auto [reference, mesh] = MakeMeshes(c);
    const RunQuadrature run(c.p);
    mfem::H1_FECollection h1_fec(c.p, 3), hat_fec(1, 3);
    mfem::RT_FECollection rt0_fec(0, 3);
    FiniteElementSpace primal(*mesh, &h1_fec), hats(*mesh, &hat_fec),
        rt0_u(*mesh, &rt0_fec);
    const auto dirichlet = c.g.Dirichlet(), terminals = c.g.Terminals();
    const auto &sides = mesh->GetCrackSides(dirichlet);
    bool split = sides.Any();
    Mpi::GlobalOr(1, &split, Mpi::World());
    std::unique_ptr<FiniteElementSpace> rt0_b;
    if (split)
    {
      rt0_b = std::make_unique<FiniteElementSpace>(rt0_u, sides);
    }
    const auto &rt0 = rt0_b ? *rt0_b : rt0_u;
    std::vector<Vector> V;
    for (int t : terminals)
    {
      V.push_back(SolveTerminal(primal, dirichlet, t).V);
    }
    config::PeriodicBoundaryData periodic;
    MaterialOperator mat_op(Materials(), periodic, ProblemType::ELECTROSTATIC, *mesh);
    const auto D0 = fem::VertexPatchFlux(primal, hats, rt0, mat_op, dirichlet, V);

    // The indicators w_k of the Dirichlet boundaries.
    const auto &pmesh = mesh->Get();
    const int max = pmesh.bdr_attributes.Max();
    std::vector<Vector> w;
    for (int attr : dirichlet)
    {
      mfem::Array<int> list;
      hats.Get().GetEssentialTrueDofs(mesh::AttrToMarker(max, std::vector<int>{attr}),
                                      list);
      Vector t(hats.GetTrueVSize()), l(hats.GetVSize());
      t = 0.0;
      t.SetSubVector(list, 1.0);
      hats.GetProlongationMatrix()->Mult(t, l);
      w.push_back(l);
    }
    const std::set<int> dset(dirichlet.begin(), dirichlet.end());
    const auto face_to_bdr = pmesh.GetFaceToBdrElMap();
    for (std::size_t t = 0; t < V.size(); t++)
    {
      Vector d_l(rt0.GetVSize()), v_l(primal.GetVSize());
      rt0.GetProlongationMatrix()->Mult(D0[t], d_l);
      primal.GetProlongationMatrix()->Mult(V[t], v_l);
      // ‖D₀‖², ‖div D₀‖², the Neumann normal trace, (D₀, ∇w_k) and (D_v, ∇w_k).
      std::vector<double> sums(3 + 2 * w.size(), 0.0);
      mfem::Array<int> dofs, pdofs, hdofs;
      mfem::DofTransformation dt;
      mfem::DenseMatrix vshape, pd, hd;
      mfem::Vector divshape, D(3), gu(3), gw(3);
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
        const auto &fe = *rt0.Get().GetFE(e), &pe = *primal.Get().GetFE(e),
                   &he = *hats.Get().GetFE(e);
        rt0.GetElementDofs(e, dofs, dt);
        mfem::Vector c_e, u_e;
        d_l.GetSubVector(dofs, c_e);
        primal.Get().GetElementDofs(e, pdofs);
        v_l.GetSubVector(pdofs, u_e);
        hats.Get().GetElementDofs(e, hdofs);
        vshape.SetSize(fe.GetDof(), 3);
        divshape.SetSize(fe.GetDof());
        pd.SetSize(pe.GetDof(), 3);
        hd.SetSize(he.GetDof(), 3);
        auto &T = *const_cast<mfem::ParMesh &>(pmesh).GetElementTransformation(e);
        const auto &ir = mfem::IntRules.Get(fe.GetGeomType(), 2 * c.p + 2);
        for (int q = 0; q < ir.GetNPoints(); q++)
        {
          const auto &ip = ir.IntPoint(q);
          T.SetIntPoint(&ip);
          const double wq = ip.weight * T.Weight();
          fe.CalcVShape(T, vshape);
          fe.CalcPhysDivShape(T, divshape);
          pe.CalcPhysDShape(T, pd);
          he.CalcPhysDShape(T, hd);
          vshape.MultTranspose(c_e, D);
          pd.MultTranspose(u_e, gu);
          const double div = divshape * c_e;
          sums[0] += wq * (D * D);
          sums[1] += wq * div * div;
          for (std::size_t k = 0; k < w.size(); k++)
          {
            mfem::Vector w_e;
            w[k].GetSubVector(hdofs, w_e);
            hd.MultTranspose(w_e, gw);
            double Dv_gw = 0.0;
            for (int d = 0; d < 3; d++)
            {
              Dv_gw -= Eps(pmesh.GetAttribute(e), d) * gu(d) * gw(d);
            }
            sums[3 + 2 * k] += wq * (D * gw);
            sums[4 + 2 * k] += wq * Dv_gw;
          }
        }
      }
      mfem::FaceElementTransformations FT;
      mfem::IsoparametricTransformation T1, T2;
      mfem::Vector normal(3);
      for (int f = 0; f < pmesh.GetNumFaces(); f++)
      {
        if (!pmesh.GetFaceInformation(f).IsBoundary() ||
            (face_to_bdr[f] >= 0 && dset.count(pmesh.GetBdrAttribute(face_to_bdr[f]))))
        {
          continue;
        }
        pmesh.GetFaceElementTransformations(f, FT, T1, T2);
        const auto &fe = *rt0.Get().GetFE(FT.Elem1No);
        rt0.GetElementDofs(FT.Elem1No, dofs, dt);
        mfem::Vector c_e;
        d_l.GetSubVector(dofs, c_e);
        vshape.SetSize(fe.GetDof(), 3);
        const auto &ir = mfem::IntRules.Get(FT.GetGeometryType(), 2);
        for (int q = 0; q < ir.GetNPoints(); q++)
        {
          FT.SetAllIntPoints(&ir.IntPoint(q));
          fe.CalcVShape(*FT.Elem1, vshape);
          vshape.MultTranspose(c_e, D);
          mfem::CalcOrtho(FT.Jacobian(), normal);
          const double Dn = (D * normal) / normal.Norml2();
          sums[2] += ir.IntPoint(q).weight * FT.Weight() * Dn * Dn;
        }
      }
      Mpi::GlobalSum(static_cast<int>(sums.size()), sums.data(), Mpi::World());
      INFO("kind " << c.g.kind << " type " << c.g.type << " refine " << c.refine << " p "
                   << c.p << " terminal " << t);
      const double norm = std::sqrt(sums[0]);
      CHECK(std::sqrt(sums[1]) <= 1.0e-10 * norm);
      CHECK(std::sqrt(sums[2]) <= 1.0e-10 * norm);
      if (c.g.kind == Geometry::BOX || c.g.kind == Geometry::TWO)
      {
        double qmax = 0.0;
        for (std::size_t k = 0; k < w.size(); k++)
        {
          qmax = std::max(qmax, std::abs(sums[4 + 2 * k]));
        }
        for (std::size_t k = 0; k < w.size(); k++)
        {
          CHECK(std::abs(sums[3 + 2 * k] - sums[4 + 2 * k]) <= 1.0e-10 * qmax);
        }
      }
    }
  }
}

// Every iterate of the curl correction is admissible and the bound decreases monotonically:
// the η² tracked from the decrease of the CG objective never increases, and after the last
// iteration it equals η² evaluated for the corrected flux.
TEST_CASE("Curl correction decreases the bound monotonically",
          "[equilibratedestimator][Serial][Parallel]")
{
  for (const auto &c : Cases())
  {
    if (c.p != 1)
    {
      continue;  // The flux space RT(0) of the vertex patches, no embedding.
    }
    auto [reference, mesh] = MakeMeshes(c);
    const RunQuadrature run(c.p);
    mfem::H1_FECollection h1_fec(c.p, 3), hat_fec(1, 3);
    mfem::RT_FECollection rt0_fec(0, 3);
    FiniteElementSpace primal(*mesh, &h1_fec), hats(*mesh, &hat_fec),
        rt0_u(*mesh, &rt0_fec);
    const auto dirichlet = c.g.Dirichlet();
    const auto &sides = mesh->GetCrackSides(dirichlet);
    bool split = sides.Any();
    Mpi::GlobalOr(1, &split, Mpi::World());
    std::unique_ptr<FiniteElementSpace> rt0_b;
    if (split)
    {
      rt0_b = std::make_unique<FiniteElementSpace>(rt0_u, sides);
    }
    const auto &rt0 = rt0_b ? *rt0_b : rt0_u;
    const auto s = SolveTerminal(primal, dirichlet, c.g.Terminals()[0]);
    config::PeriodicBoundaryData periodic;
    MaterialOperator mat_op(Materials(), periodic, ProblemType::ELECTROSTATIC, *mesh);
    const auto D0 = fem::VertexPatchFlux(primal, hats, rt0, mat_op, dirichlet, {s.V});
    Vector flux_l(rt0.GetVSize()), v_l(primal.GetVSize());
    rt0.GetProlongationMatrix()->Mult(D0[0], flux_l);
    primal.GetProlongationMatrix()->Mult(s.V, v_l);
    const Vector flux0(flux_l);
    equilibrated::CellQuadrature quadrature(rt0, primal);
    auto linear = Linear();
    linear.estimator_error_tol = 1.0e-3;  // Stopped early: the history is what is checked.
    CurlCorrectionSolver curl(mat_op, primal, rt0, dirichlet, nullptr, linear, quadrature);
    const double unit[2] = {0.0, 1.0};  // The Dirichlet range of the terminal solve.
    const auto result = curl.Correct(flux_l, v_l, unit);

    // η² and ‖div D‖ of a corrected flux, and ‖D‖.
    const auto &pmesh = mesh->Get();
    mfem::Array<int> dofs, pdofs;
    mfem::DofTransformation dt;
    mfem::DenseMatrix vshape, pd;
    mfem::Vector D(3), gu(3), divshape;
    auto Evaluate = [&](const Vector &f)
    {
      std::array<double, 3> sums = {0.0, 0.0, 0.0};
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
        const auto &fe = *rt0.Get().GetFE(e), &pe = *primal.Get().GetFE(e);
        rt0.GetElementDofs(e, dofs, dt);
        mfem::Vector c_e, u_e;
        f.GetSubVector(dofs, c_e);
        primal.Get().GetElementDofs(e, pdofs);
        v_l.GetSubVector(pdofs, u_e);
        vshape.SetSize(fe.GetDof(), 3);
        pd.SetSize(pe.GetDof(), 3);
        divshape.SetSize(fe.GetDof());
        auto &T = *const_cast<mfem::ParMesh &>(pmesh).GetElementTransformation(e);
        const auto &ir = mfem::IntRules.Get(fe.GetGeomType(), 4);
        for (int q = 0; q < ir.GetNPoints(); q++)
        {
          T.SetIntPoint(&ir.IntPoint(q));
          fe.CalcVShape(T, vshape);
          fe.CalcDivShape(ir.IntPoint(q), divshape);
          pe.CalcPhysDShape(T, pd);
          vshape.MultTranspose(c_e, D);
          pd.MultTranspose(u_e, gu);
          double r2 = 0.0;
          for (int d = 0; d < 3; d++)
          {
            const double eps = Eps(pmesh.GetAttribute(e), d);
            const double r = D(d) + eps * gu(d);  // D - D_v.
            r2 += r * r / eps;
          }
          const double w = ir.IntPoint(q).weight * T.Weight(),
                       div = (divshape * c_e) / T.Weight();
          sums[0] += w * r2;
          sums[1] += w * div * div;
          sums[2] += w * (D * D);
        }
      }
      Mpi::GlobalSum(3, sums.data(), Mpi::World());
      return sums;
    };
    const double eta2 = Evaluate(flux_l)[0];
    INFO("kind " << c.g.kind << " type " << c.g.type << " refine " << c.refine << ": "
                 << result.iterations << " iterations, η² " << result.eta2.front() << " -> "
                 << eta2);
    CHECK(result.iterations > 0);
    CHECK(result.accepted);
    REQUIRE(static_cast<int>(result.eta2.size()) == result.iterations + 1);
    bool monotone = true;
    for (std::size_t k = 1; k < result.eta2.size(); k++)
    {
      monotone = monotone && result.eta2[k] <= result.eta2[k - 1] * (1.0 + 1.0e-14);
    }
    CHECK(monotone);
    CHECK(std::abs(result.eta2.back() - eta2) <= 1.0e-8 * result.eta2.front());

    // The decision on the assembled flux: the corrected flux has no larger bound, and three
    // times its curl, admissible but with η² = η₀² + 3 (A, K A) for a CG iterate A ((b, A)
    // = (A, K A)), is rejected and leaves the flux unchanged.
    {
      Vector curl_l(flux_l);
      curl_l -= flux0;
      CHECK(curl.Bound(flux_l, v_l, unit) <= curl.Bound(flux0, v_l, unit));
      Vector rejected(flux0), three(curl_l);
      three *= 3.0;
      CHECK_FALSE(curl.AddIfNotWorse(rejected, three, v_l, unit));
      rejected -= flux0;
      CHECK(rejected.Normlinf() == 0.0);
      Vector kept(flux0);
      CHECK(curl.AddIfNotWorse(kept, curl_l, v_l, unit));
      kept -= flux_l;
      CHECK(kept.Normlinf() == 0.0);
    }

    // Converged, and without a tolerance (the data scaled so that the absolute tolerance of
    // CG is not reached, at most 2000 iterations, stopped once η² no longer decreases): the
    // gradient part of the iterate, which the curl does not see, may grow, but here the
    // corrected flux stays divergence free to rounding and the bound does not grow (Correct
    // would keep the uncorrected flux otherwise).
    {
      constexpr double scale = 1.0e8;
      auto Run = [&](double tol)
      {
        auto long_linear = Linear();
        long_linear.estimator_tol = tol;
        long_linear.estimator_error_tol = 0.0;
        long_linear.estimator_max_it = 2000;
        CurlCorrectionSolver long_curl(mat_op, primal, rt0, dirichlet, nullptr, long_linear,
                                       quadrature);
        Vector long_flux(rt0.GetVSize()), long_v(v_l);
        rt0.GetProlongationMatrix()->Mult(D0[0], long_flux);
        long_flux *= scale;
        long_v *= scale;
        const double scaled[2] = {0.0, scale};
        const auto long_result = long_curl.Correct(long_flux, long_v, scaled);
        long_flux *= 1.0 / scale;
        return std::make_pair(long_result, Evaluate(long_flux));
      };
      const auto [converged, converged_sums] = Run(1.0e-12);
      const auto [untol, sums] = Run(0.0);
      INFO("converged: " << converged.iterations << " iterations, ‖A‖ " << converged.norm
                         << "; without tolerance: " << untol.iterations
                         << " iterations, ‖A‖ " << untol.norm << ", ‖div D‖ "
                         << std::sqrt(sums[1]) << ", ‖D‖ " << std::sqrt(sums[2]) << ", η² "
                         << sums[0]);
      CHECK(converged.accepted);
      CHECK(untol.accepted);
      CHECK(untol.iterations >= converged.iterations);
      CHECK(untol.iterations < 2000);
      CHECK(std::sqrt(sums[1]) <= 1.0e-10 * std::sqrt(sums[2]));
      CHECK(sums[0] <= eta2 * (1.0 + 1.0e-10));
      CHECK(sums[0] <= converged_sums[0] * (1.0 + 1.0e-10));
    }
  }
}

// The relative bound and its uncertified limit.
TEST_CASE("Equilibrated relative bound", "[equilibratedestimator][Serial]")
{
  using E = EquilibratedFluxErrorEstimator;
  CHECK(E::RelativeBound(0.0, 0.0, 1.0) == 0.0);
  // A / ((2 U_h)^{1/2} - β)² - L²)^{1/2}: 0.6 / (1 - 0.36)^{1/2} = 0.75.
  CHECK(std::abs(E::RelativeBound(0.6, 0.0, 0.5) - 0.75) <= 1.0e-14);
  CHECK(std::abs(E::RelativeBound(0.3, 0.3, 0.5) - 0.6 / std::sqrt(0.49 - 0.09)) <=
        1.0e-14);
  CHECK(std::isinf(E::RelativeBound(1.0, 0.0, 0.5)));
  CHECK(std::isinf(E::RelativeBound(0.5, 1.0, 0.5)));
}

}  // namespace palace
