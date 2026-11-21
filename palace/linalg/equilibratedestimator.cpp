// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "equilibratedestimator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <vector>
#include "fem/bilinearform.hpp"
#include "fem/errorindicator.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/vertexpatchflux.hpp"
#include "linalg/iterative.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// The values of the element DOFs of an L-vector, in the reference orientation.
mfem::Vector ElementValues(const FiniteElementSpace &space, int e, const Vector &local)
{
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  space.GetElementDofs(e, dofs, dof_trans);
  mfem::Vector values(dofs.Size());
  local.GetSubVector(dofs, values);
  dof_trans.InvTransformPrimal(values);
  return values;
}

// Set the element DOFs of an L-vector from values in the reference orientation.
void SetElementValues(const FiniteElementSpace &space, int e, mfem::Vector &values,
                      Vector &local)
{
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  space.GetElementDofs(e, dofs, dof_trans);
  dof_trans.TransformPrimal(values);
  local.SetSubVector(dofs, values);
}

// Whether the coordinate map of the element is affine (its nodes are those of the affine
// map of its Jacobian at the first vertex).
bool IsAffine(mfem::IsoparametricTransformation &T)
{
  const auto &fe = *T.GetFE();
  const auto &v0 = mfem::Geometries.GetVertices(T.GetGeometryType())->IntPoint(0);
  T.SetIntPoint(&v0);
  const mfem::DenseMatrix J = T.Jacobian();
  mfem::Vector x0(3), x(3);
  T.Transform(v0, x0);
  double extent = 0.0, deviation = 0.0;
  for (int i = 0; i < fe.GetNodes().GetNPoints(); i++)
  {
    const auto &ip = fe.GetNodes().IntPoint(i);
    T.Transform(ip, x);
    const double xi[3] = {ip.x - v0.x, ip.y - v0.y, ip.z - v0.z};
    for (int k = 0; k < 3; k++)
    {
      extent = std::max(extent, std::abs(x(k) - x0(k)));
      deviation = std::max(deviation, std::abs(x(k) - x0(k) - J(k, 0) * xi[0] -
                                               J(k, 1) * xi[1] - J(k, 2) * xi[2]));
    }
  }
  return deviation <= 1.0e-12 * extent;
}

// The exact reference matrix X of the basis functions of a space in those of a space
// containing it, on the reference element: lo_i = Σ_l X_li hi_l, from the mass matrix of
// hi. Both use the same Piola map, so X holds on every cell of the geometry (with the DOF
// transformations of each space). For the curl, lo is the curl of a Nédélec element.
mfem::DenseMatrix ReferenceEmbedding(const mfem::FiniteElement &lo,
                                     const mfem::FiniteElement &hi, bool curl)
{
  const int nl = lo.GetDof(), nh = hi.GetDof();
  const auto &ir = mfem::IntRules.Get(hi.GetGeomType(), 2 * hi.GetOrder() + 2);
  mfem::DenseMatrix M(nh), G(nh, nl), sh(nh, 3), sl(nl, 3), X(nh, nl);
  M = 0.0;
  G = 0.0;
  auto Shapes = [&](const mfem::IntegrationPoint &ip)
  {
    hi.CalcVShape(ip, sh);
    curl ? lo.CalcCurlShape(ip, sl) : lo.CalcVShape(ip, sl);
  };
  for (int q = 0; q < ir.GetNPoints(); q++)
  {
    Shapes(ir.IntPoint(q));
    mfem::AddMult_a_ABt(ir.IntPoint(q).weight, sh, sh, M);
    mfem::AddMult_a_ABt(ir.IntPoint(q).weight, sh, sl, G);
  }
  mfem::DenseMatrixInverse inverse(M);
  inverse.Mult(G, X);
  double defect = 0.0, scale = 0.0;
  mfem::DenseMatrix fit(nl, 3);
  for (int q = 0; q < ir.GetNPoints(); q++)
  {
    Shapes(ir.IntPoint(q));
    mfem::MultAtB(X, sh, fit);
    for (int i = 0; i < nl; i++)
    {
      for (int d = 0; d < 3; d++)
      {
        defect = std::max(defect, std::abs(fit(i, d) - sl(i, d)));
        scale = std::max(scale, std::abs(sl(i, d)));
      }
    }
  }
  MFEM_VERIFY(defect <= 1.0e-10 * scale, "The reference basis is not contained in the flux "
                                         "space (defect "
                                             << defect << " of " << scale << ")!");
  return X;
}

// The Jacobian J, its determinant and adjugate at quadrature point q of a cell (once on
// affine cells).
struct CellMap
{
  mfem::DenseMatrix J{3}, adjJ{3};
  double det = 0.0;
  void Update(const mfem::DenseMatrix &points, const mfem::DenseMatrix &map_dshape)
  {
    mfem::Mult(points, map_dshape, J);
    det = J.Det();
    mfem::CalcAdjugate(J, adjJ);
  }
};

// ε and ε⁻¹ of the cell (copies: the material operator returns views of its data, which an
// assignment would alias).
struct CellMaterial
{
  mfem::DenseMatrix eps{3}, eps_inv{3};
  void Update(const MaterialOperator &mat_op, int attr)
  {
    eps.Set(1.0, mat_op.GetPermittivityReal(attr));
    eps_inv.Set(1.0, eps);
    eps_inv.Invert();
  }
};

// The connected components of the Dirichlet boundary: the vertices of the boundary
// elements with a Dirichlet attribute, connected by these elements and, on a nonconforming
// mesh, by the constraints of their hanging vertices (the true DOFs of the order-one H1
// space hats in the rows of its prolongation P for these vertices). Each true DOF starts
// with its global index as its label; in each round the owners send the labels to the
// processes referencing the true DOFs, which propose the smallest label of each component
// of their local graph to the owners, until no label changes. Returns the label of each
// owned true DOF of hats (-1 off the Dirichlet boundary), and for each vertex of the mesh
// that of its component (-1 off the Dirichlet boundary). Collective.
std::vector<std::int64_t> DirichletComponents(const FiniteElementSpace &hats,
                                              const std::vector<int> &dirichlet_attr,
                                              std::vector<std::int64_t> &vertex_label)
{
  const auto &mesh = hats.GetParMesh();
  const MPI_Comm comm = hats.GetComm();
  const int nprocs = Mpi::Size(comm);
  const std::set<int> dirichlet(dirichlet_attr.begin(), dirichlet_attr.end());
  const auto &P = *hats.Dof_TrueDof_Matrix();
  P.HostRead();
  hypre_ParCSRMatrix *hP = const_cast<mfem::HypreParMatrix &>(P);
  const hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(hP),
                        *offd = hypre_ParCSRMatrixOffd(hP);
  const HYPRE_Int *di = hypre_CSRMatrixI(diag), *dj = hypre_CSRMatrixJ(diag),
                  *oi = hypre_CSRMatrixI(offd), *oj = hypre_CSRMatrixJ(offd);
  const HYPRE_Real *dv = hypre_CSRMatrixData(diag), *ov = hypre_CSRMatrixData(offd);
  const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(hP);
  const HYPRE_BigInt col0 = hypre_ParCSRMatrixFirstColDiag(hP);
  const int tsize = hats.GetTrueVSize();

  // The true DOFs of the vertices of each Dirichlet boundary element (the L-DOF of a
  // vertex of the order-one space is its index).
  std::vector<int> group_start(1, 0);
  std::vector<std::int64_t> cols;
  std::vector<int> vertex_entry(mesh.GetNV(), -1);
  mfem::Array<int> vertices;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (!dirichlet.count(mesh.GetBdrAttribute(be)))
    {
      continue;
    }
    mesh.GetBdrElementVertices(be, vertices);
    for (int v : vertices)
    {
      const std::size_t first = cols.size();
      for (HYPRE_Int k = di[v]; k < di[v + 1]; k++)
      {
        if (dv[k] != 0.0)
        {
          cols.push_back(col0 + dj[k]);
        }
      }
      for (HYPRE_Int k = (oi ? oi[v] : 0); k < (oi ? oi[v + 1] : 0); k++)
      {
        if (ov[k] != 0.0)
        {
          cols.push_back(cmap[oj[k]]);
        }
      }
      if (cols.size() > first)
      {
        vertex_entry[v] = static_cast<int>(first);
      }
    }
    group_start.push_back(static_cast<int>(cols.size()));
  }
  P.HypreRead();
  std::vector<HYPRE_BigInt> col_starts(nprocs + 1);
  MPI_Allgather(&col0, 1, HYPRE_MPI_BIG_INT, col_starts.data(), 1, HYPRE_MPI_BIG_INT, comm);
  col_starts[nprocs] = hypre_ParCSRMatrixGlobalNumCols(hP);
  auto Owner = [&](std::int64_t g)
  {
    return static_cast<int>(std::upper_bound(col_starts.begin(), col_starts.end(),
                                             static_cast<HYPRE_BigInt>(g)) -
                            col_starts.begin()) -
           1;
  };
  std::vector<std::int64_t> ref(cols);
  std::sort(ref.begin(), ref.end());
  ref.erase(std::unique(ref.begin(), ref.end()), ref.end());
  const int nref = static_cast<int>(ref.size());
  auto Index = [&](std::int64_t g)
  { return static_cast<int>(std::lower_bound(ref.begin(), ref.end(), g) - ref.begin()); };
  std::vector<int> comp(nref);
  {
    std::vector<int> parent(nref);
    std::iota(parent.begin(), parent.end(), 0);
    auto Find = [&](int k)
    {
      while (parent[k] != k)
      {
        parent[k] = parent[parent[k]];
        k = parent[k];
      }
      return k;
    };
    for (std::size_t g = 0; g + 1 < group_start.size(); g++)
    {
      if (group_start[g + 1] == group_start[g])
      {
        continue;
      }
      const int root = Find(Index(cols[group_start[g]]));
      for (int e = group_start[g] + 1; e < group_start[g + 1]; e++)
      {
        const int other = Find(Index(cols[e]));
        if (other != root)
        {
          parent[other] = root;
        }
      }
    }
    for (int k = 0; k < nref; k++)
    {
      comp[k] = Find(k);
    }
  }
  // The referenced true DOFs by owner (ref is sorted, so the replies of the owners, in
  // their order, are in the order of ref).
  std::vector<std::vector<std::int64_t>> request(nprocs);
  for (int k = 0; k < nref; k++)
  {
    request[Owner(ref[k])].push_back(ref[k]);
  }
  std::vector<int> asked_counts, recv_counts;
  const auto asked = mpi::ExchangeAll(comm, request, asked_counts);
  request.clear();
  std::vector<std::int64_t> label(tsize, -1), ghost, low(nref);
  for (const auto t : asked)
  {
    label[t - col0] = t;
  }
  while (true)
  {
    std::vector<std::vector<std::int64_t>> reply(nprocs);
    for (int p = 0, pos = 0; p < nprocs; p++)
    {
      for (int i = 0; i < asked_counts[p]; i++, pos++)
      {
        reply[p].push_back(label[asked[pos] - col0]);
      }
    }
    ghost = mpi::ExchangeAll(comm, reply, recv_counts);
    std::fill(low.begin(), low.end(), std::numeric_limits<std::int64_t>::max());
    for (int k = 0; k < nref; k++)
    {
      low[comp[k]] = std::min(low[comp[k]], ghost[k]);
    }
    std::vector<std::vector<std::int64_t>> propose(nprocs);
    int changed = 0;
    for (int k = 0; k < nref; k++)
    {
      if (low[comp[k]] < ghost[k])
      {
        propose[Owner(ref[k])].insert(propose[Owner(ref[k])].end(), {ref[k], low[comp[k]]});
        changed = 1;
      }
    }
    Mpi::GlobalMax(1, &changed, comm);
    if (!changed)
    {
      break;
    }
    const auto proposed = mpi::ExchangeAll(comm, propose, recv_counts);
    for (std::size_t i = 0; i < proposed.size(); i += 2)
    {
      auto &l = label[proposed[i] - col0];
      l = std::min(l, proposed[i + 1]);
    }
  }
  vertex_label.assign(mesh.GetNV(), -1);
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    if (vertex_entry[v] >= 0)
    {
      vertex_label[v] = low[comp[Index(cols[vertex_entry[v]])]];
    }
  }
  return label;
}

// The order-one potentials which, with the potentials of the solves, span the potentials
// constant on each connected component of the Dirichlet boundary (up to a constant): for
// each component (in order) whose indicator is not in the span of the Dirichlet data of the
// solves (the companions, L-vectors of a space whose L-DOF of a vertex is its index), the
// constant and the indicators of the components chosen before, the solution of div(ε ∇w) =
// 0 with w = 1 on the component and 0 on the rest of the Dirichlet boundary, in the
// order-one space hats with the linear solver of the run. Components on which the data of
// a solve is not constant are left out. Returns the true-DOF vectors of hats, and the
// number of components. Collective.
std::vector<Vector> HarmonicPotentials(FiniteElementSpace &hats,
                                       const MaterialOperator &mat_op,
                                       const std::vector<int> &dirichlet_attr,
                                       const std::vector<Vector> &companions,
                                       const config::LinearSolverData &linear,
                                       int &num_components)
{
  const MPI_Comm comm = hats.GetComm();
  const auto &mesh = hats.GetParMesh();
  std::vector<std::int64_t> vertex_label;
  const auto label = DirichletComponents(hats, dirichlet_attr, vertex_label);

  // The components, by their labels (the smallest true DOF of each), on all processes.
  std::vector<std::int64_t> components;
  {
    const HYPRE_BigInt col0 = hats.Dof_TrueDof_Matrix()->ColPart()[0];
    std::vector<std::int64_t> roots;
    for (std::size_t t = 0; t < label.size(); t++)
    {
      if (label[t] == col0 + static_cast<std::int64_t>(t))
      {
        roots.push_back(label[t]);
      }
    }
    const int nprocs = Mpi::Size(comm);
    std::vector<int> counts(nprocs), displs(nprocs, 0);
    const int nroots = static_cast<int>(roots.size());
    Mpi::Allgather(1, &nroots, counts.data(), comm);
    for (int p = 1; p < nprocs; p++)
    {
      displs[p] = displs[p - 1] + counts[p - 1];
    }
    components.resize(displs[nprocs - 1] + counts[nprocs - 1]);
    MPI_Allgatherv(roots.data(), nroots, mpi::DataType<std::int64_t>(), components.data(),
                   counts.data(), displs.data(), mpi::DataType<std::int64_t>(), comm);
    std::sort(components.begin(), components.end());
  }
  const int nc = static_cast<int>(components.size()),
            n = static_cast<int>(companions.size());
  num_components = nc;
  auto Component = [&](std::int64_t l)
  {
    return static_cast<int>(std::lower_bound(components.begin(), components.end(), l) -
                            components.begin());
  };

  // The Dirichlet data of each solve on each component (minimum and maximum).
  std::vector<double> lo(nc * n, std::numeric_limits<double>::infinity()),
      hi(nc * n, std::numeric_limits<double>::infinity());
  for (int v = 0; v < mesh.GetNV(); v++)
  {
    if (vertex_label[v] < 0)
    {
      continue;
    }
    const int c = Component(vertex_label[v]);
    for (int i = 0; i < n; i++)
    {
      const double value = companions[i].HostRead()[v];
      lo[c * n + i] = std::min(lo[c * n + i], value);
      hi[c * n + i] = std::min(hi[c * n + i], -value);
    }
  }
  Mpi::GlobalMin(nc * n, lo.data(), comm);
  Mpi::GlobalMin(nc * n, hi.data(), comm);

  // Greedy completion of the span of the data of the solves and of the constant (Gram-
  // Schmidt in the space of the values on the components).
  std::vector<int> chosen;
  {
    std::vector<mfem::Vector> basis;
    auto Add = [&](mfem::Vector u)
    {
      for (const auto &b : basis)
      {
        u.Add(-(u * b), b);
      }
      const double norm = u.Norml2();
      if (norm > 1.0e-6)
      {
        u /= norm;
        basis.push_back(u);
        return true;
      }
      return false;
    };
    std::vector<bool> constant(nc, true);
    for (int i = 0; i < n; i++)
    {
      mfem::Vector u(nc);
      double scale = 0.0;
      for (int c = 0; c < nc; c++)
      {
        u(c) = lo[c * n + i];
        scale = std::max(scale, std::abs(u(c)));
      }
      for (int c = 0; c < nc; c++)
      {
        constant[c] = constant[c] &&
                      (-hi[c * n + i] - lo[c * n + i] <= 1.0e-12 * std::max(scale, 1.0));
      }
      if (scale > 0.0)
      {
        u /= scale;
        Add(u);
      }
    }
    mfem::Vector one(nc);
    one = 1.0;
    Add(one);
    for (int c = 0; c < nc && static_cast<int>(basis.size()) < nc; c++)
    {
      mfem::Vector e(nc);
      e = 0.0;
      e(c) = 1.0;
      if (constant[c] && Add(e))
      {
        chosen.push_back(c);
      }
    }
  }
  if (chosen.empty())
  {
    return {};
  }

  // The potentials, with the linear solver of the run for the order-one space.
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> essential;
  hats.Get().GetEssentialTrueDofs(mesh::AttrToMarker(bdr_attr_max, dirichlet_attr),
                                  essential);
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  BilinearForm a(hats);
  a.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
  ParOperator A(a.Assemble(false), hats);
  A.SetEssentialTrueDofs(essential, Operator::DiagonalPolicy::DIAG_ONE);
  FiniteElementSpaceHierarchy spaces;
  spaces.AddLevel(hats);
  config::LinearSolverData solver = linear;
  solver.krylov_solver = KrylovSolver::CG;
  solver.initial_guess = 0;
  solver.tol = 1.0e-12;  // Its residual remains as the divergence of the flux.
  solver.abs_tol = std::numeric_limits<double>::epsilon();
  KspSolver ksp(solver, MatrixSymmetry::SPD, 0, spaces, nullptr, false);
  ksp.SetOperators(A, A);
  std::vector<Vector> potentials;
  for (int c : chosen)
  {
    Vector x(hats.GetTrueVSize()), b(hats.GetTrueVSize()), w(hats.GetTrueVSize());
    x.UseDevice(true);
    b.UseDevice(true);
    w.UseDevice(true);
    x = 0.0;
    {
      auto *h_x = x.HostReadWrite();
      for (int t : essential)
      {
        h_x[t] = (label[t] == components[c]) ? 1.0 : 0.0;
      }
    }
    b = 0.0;
    A.EliminateRHS(x, b);
    w = 0.0;
    ksp.Mult(b, w);
    potentials.push_back(std::move(w));
  }
  return potentials;
}

// The terms of the bound of a flux D (an L-vector of the flux space) for the companion v
// (an L-vector of the primal space) with Dirichlet range c ± R (see
// EquilibratedFluxErrorEstimator): per element, the mismatch η_K², the balance term b_K
// (from the divergence and from the Neumann faces of the element) and, with a lifting
// L-vector, the lifting β_K²; and the totals of ‖div D‖² and of ‖D ⋅ n‖² on the Neumann
// faces (local).
struct BoundTerms
{
  std::vector<double> mismatch, balance, lift;
  double div2 = 0.0, neumann2 = 0.0;
};

BoundTerms EvaluateBoundTerms(const MaterialOperator &mat_op,
                              const FiniteElementSpace &flux,
                              const FiniteElementSpace &primal,
                              const equilibrated::CellQuadrature &quadrature,
                              const std::vector<int> &dirichlet_attr, const Vector &flux_l,
                              const Vector &companion_l, double center, double radius,
                              const Vector *lifting_l = nullptr)
{
  const auto &mesh = primal.GetParMesh();
  const int ne = mesh.GetNE();
  const std::set<int> dirichlet(dirichlet_attr.begin(), dirichlet_attr.end());
  BoundTerms terms;
  terms.mismatch.assign(ne, 0.0);
  terms.balance.assign(ne, 0.0);
  terms.lift.assign(ne, 0.0);

  // Cell integrals.
  mfem::IsoparametricTransformation T;
  CellMap map;
  CellMaterial material;
  mfem::Vector ref(3), D(3), grad_ref(3), grad(3), Dv(3), r(3), Er(3), du;
  for (int e = 0; e < ne; e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &tab = quadrature.GetTables(e, T);
    const auto c = ElementValues(flux, e, flux_l);
    const auto u = ElementValues(primal, e, companion_l);
    if (lifting_l)
    {
      du = ElementValues(primal, e, *lifting_l);
    }
    material.Update(mat_op, mesh.GetAttribute(e));
    double eta2 = 0.0, div2 = 0.0, volume = 0.0, offset2 = 0.0, beta2 = 0.0;
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      if (q == 0 || !quadrature.IsAffine(e))
      {
        map.Update(T.GetPointMat(), tab.map[q]);
      }
      const double w = tab.rule->IntPoint(q).weight * std::abs(map.det);
      tab.flux[q].MultTranspose(c, ref);
      map.J.Mult(ref, D);
      D /= map.det;
      const double div = (tab.div[q] * c) / map.det;
      tab.grad[q].MultTranspose(u, grad_ref);
      map.adjJ.MultTranspose(grad_ref, grad);
      grad /= map.det;
      material.eps.Mult(grad, Dv);
      Dv.Neg();
      subtract(D, Dv, r);
      material.eps_inv.Mult(r, Er);
      eta2 += w * (r * Er);
      div2 += w * div * div;
      volume += w;
      const double offset = (tab.shape[q] * u) - center;
      offset2 += w * offset * offset;
      if (lifting_l)
      {
        tab.grad[q].MultTranspose(du, grad_ref);
        map.adjJ.MultTranspose(grad_ref, grad);
        grad /= map.det;
        material.eps.Mult(grad, Dv);
        beta2 += w * (grad * Dv);
      }
    }
    terms.mismatch[e] = eta2;
    terms.div2 += div2;
    terms.balance[e] = std::sqrt(div2) * (radius * std::sqrt(volume) + std::sqrt(offset2));
    terms.lift[e] = beta2;
  }

  // Neumann faces (exterior faces without a Dirichlet attribute): the normal trace of the
  // flux is zero there up to rounding, its actual value enters the balance term.
  const auto face_to_bdr = mesh.GetFaceToBdrElMap();
  mfem::Vector normal(3), pshape;
  mfem::DenseMatrix shape;
  mfem::FaceElementTransformations FT;
  mfem::IsoparametricTransformation T1, T2;
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    if (!mesh.GetFaceInformation(f).IsBoundary() ||
        (face_to_bdr[f] >= 0 && dirichlet.count(mesh.GetBdrAttribute(face_to_bdr[f]))))
    {
      continue;
    }
    mesh.GetFaceElementTransformations(f, FT, T1, T2);
    const int e = FT.Elem1No;
    const auto &fe = *flux.Get().GetFE(e), &pe = *primal.Get().GetFE(e);
    const auto c = ElementValues(flux, e, flux_l);
    const auto u = ElementValues(primal, e, companion_l);
    shape.SetSize(fe.GetDof(), 3);
    pshape.SetSize(pe.GetDof());
    double flux2 = 0.0, area = 0.0, offset2 = 0.0;
    const auto &rule = mfem::IntRules.Get(FT.GetGeometryType(), quadrature.Order());
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &ip = rule.IntPoint(q);
      FT.SetAllIntPoints(&ip);
      const double w = ip.weight * FT.Weight();
      fe.CalcVShape(*FT.Elem1, shape);
      shape.MultTranspose(c, D);
      mfem::CalcOrtho(FT.Jacobian(), normal);
      const double Dn = (D * normal) / normal.Norml2();
      flux2 += w * Dn * Dn;
      area += w;
      pe.CalcShape(FT.GetElement1IntPoint(), pshape);
      const double offset = (pshape * u) - center;
      offset2 += w * offset * offset;
    }
    terms.neumann2 += flux2;
    terms.balance[e] += std::sqrt(flux2) * (radius * std::sqrt(area) + std::sqrt(offset2));
  }
  return terms;
}

// L = (η + (η² + 4 b)^{1/2}) / 2, the bound of |||u - v||| without the lifting.
double BoundL(double eta2, double b)
{
  const double eta = std::sqrt(std::max(eta2, 0.0));
  return 0.5 * (eta + std::sqrt(eta * eta + 4.0 * std::max(b, 0.0)));
}

}  // namespace

namespace equilibrated
{

CellQuadrature::CellQuadrature(const FiniteElementSpace &flux,
                               const FiniteElementSpace &primal, int degree)
  : flux(flux), primal(primal),
    degree((degree < 0) ? 2 * primal.GetMaxElementOrder() + 2 : degree)
{
  const auto &mesh = primal.GetParMesh();
  mfem::IsoparametricTransformation T;
  affine.resize(mesh.GetNE());
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    affine[e] = palace::IsAffine(T);
  }
}

int CellQuadrature::Order(bool curved) const
{
  // On affine cells the integrands have reference degree 2p + 2 at most (the RT metric and
  // the mismatch |D - D_v|², with RT(p)); a coordinate map of order q adds 2 (q - 1).
  const auto *nodes = primal.GetParMesh().GetNodes();
  const int q = (nodes && curved) ? nodes->FESpace()->GetMaxElementOrder() : 1;
  return degree + 2 * (q - 1);
}

const QuadratureTables &
CellQuadrature::GetTables(int e, mfem::IsoparametricTransformation &T) const
{
  auto &tab = tables[{T.GetGeometryType(), affine[e]}];
  if (!tab.rule)
  {
    tab.rule = &mfem::IntRules.Get(T.GetGeometryType(), Order(!affine[e]));
    const auto &fe = *flux.Get().GetFE(e), &pe = *primal.Get().GetFE(e);
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      const auto &ip = tab.rule->IntPoint(q);
      tab.flux.emplace_back(fe.GetDof(), 3);
      fe.CalcVShape(ip, tab.flux.back());
      tab.div.emplace_back(fe.GetDof());
      fe.CalcDivShape(ip, tab.div.back());
      tab.map.emplace_back(T.GetFE()->GetDof(), 3);
      T.GetFE()->CalcDShape(ip, tab.map.back());
      tab.grad.emplace_back(pe.GetDof(), 3);
      pe.CalcDShape(ip, tab.grad.back());
      tab.shape.emplace_back(pe.GetDof());
      pe.CalcShape(ip, tab.shape.back());
    }
  }
  return tab;
}

}  // namespace equilibrated

CurlCorrectionSolver::CurlCorrectionSolver(const MaterialOperator &mat_op,
                                           FiniteElementSpace &primal,
                                           const FiniteElementSpace &flux,
                                           const std::vector<int> &dirichlet_attr,
                                           FiniteElementSpace *shared_nd,
                                           const config::LinearSolverData &linear,
                                           const equilibrated::CellQuadrature &quadrature)
  : mat_op(mat_op), primal(primal), flux(flux), quadrature(quadrature),
    dirichlet_attr(dirichlet_attr)
{
  auto &mesh = primal.GetMesh();
  const auto &pmesh = mesh.Get();
  const MPI_Comm comm = primal.GetComm();

  // The de Rham partner of the flux space RT(k) (MFEM order k + 1) is ND(k + 1).
  const int p_nd = flux.GetFEColl().GetOrder();
  ParOperator *K_coarse = nullptr;

  // Neumann attributes: those of exterior faces which are not Dirichlet attributes. The
  // tangential trace of the correction vanishes on them (an interior boundary with such an
  // attribute keeps the normal trace of the flux too, which is admissible).
  const std::set<int> dirichlet(dirichlet_attr.begin(), dirichlet_attr.end());
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> neumann_marker(bdr_attr_max);
  neumann_marker = 0;
  {
    const auto face_to_bdr = pmesh.GetFaceToBdrElMap();
    int missing = 0;
    for (int f = 0; f < pmesh.GetNumFaces(); f++)
    {
      if (!pmesh.GetFaceInformation(f).IsBoundary())
      {
        continue;
      }
      if (face_to_bdr[f] < 0)
      {
        missing++;
        continue;
      }
      const int attr = pmesh.GetBdrAttribute(face_to_bdr[f]);
      if (!dirichlet.count(attr))
      {
        neumann_marker[attr - 1] = 1;
      }
    }
    Mpi::GlobalSum(1, &missing, comm);
    MFEM_VERIFY(missing == 0, "The equilibrated estimator requires a boundary element on "
                              "every exterior face!");
    Mpi::GlobalMax(bdr_attr_max, neumann_marker.GetData(), comm);
  }

  // Two p-levels for the correction space, ND(1) for AMS and ND(p) (one if the multigrid
  // of the run has one level), with the order-one H1 space for AMS on the coarsest level,
  // broken across the interior Dirichlet boundaries when the flux space is, and the
  // smoother of the run on ND(p). The intermediate levels of the run's coarsening would add
  // their spaces, broken views, operators and work vectors: with them, on a refined mesh of
  // 4.2M unknowns with an interior boundary at 144 processes, the peak of the estimator
  // per process exceeds the field solve's (1.01 times, against 0.96), for 8-15% less time
  // of the correction on refined meshes and 5-37% more on first meshes.
  nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(p_nd, 3, 1,
                                                               linear.mg_coarsening, false);
  if (linear.mg_max_levels > 1 && p_nd > 1)
  {
    auto coarse = fem::ConstructFECollections<mfem::ND_FECollection>(
        1, 3, 1, linear.mg_coarsening, false);
    nd_fecs.insert(nd_fecs.begin(), std::move(coarse.front()));
  }
  h1_fec = std::make_unique<mfem::H1_FECollection>(nd_fecs.front()->GetOrder(), 3);
  const bool share = shared_nd && &shared_nd->GetMesh() == &primal.GetMesh() &&
                     !shared_nd->IsBroken() && shared_nd->GetVDim() == 1 &&
                     dynamic_cast<const mfem::ND_FECollection *>(&shared_nd->GetFEColl()) &&
                     shared_nd->GetFEColl().GetOrder() == p_nd;
  for (std::size_t l = 0; l < nd_fecs.size(); l++)
  {
    if (share && l + 1 == nd_fecs.size())
    {
      shared = shared_nd;
      nd_spaces.AddLevel(*shared_nd);
    }
    else
    {
      nd_spaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, nd_fecs[l].get()));
    }
  }
  h1_spaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, h1_fec.get()));
  broken = flux.IsBroken();
  if (broken)
  {
    const auto &sides = mesh.GetCrackSides(dirichlet_attr);
    for (std::size_t l = 0; l < nd_fecs.size(); l++)
    {
      nd_broken.AddLevel(
          std::make_unique<FiniteElementSpace>(nd_spaces.GetFESpaceAtLevel(l), sides));
      // The prolongations only enter the preconditioner: partially assembled, a copy which
      // no row of the broken prolongation determines is left at zero (a smaller coarse
      // space) rather than assembling the level.
      nd_broken.GetFESpaceAtLevel(l).SetUndeterminedRestriction(true);
    }
    h1_broken.AddLevel(
        std::make_unique<FiniteElementSpace>(h1_spaces.GetFESpaceAtLevel(0), sides));
  }
  auto &nd_hierarchy = NDSpaces();
  essential.resize(nd_hierarchy.GetNumLevels());
  for (std::size_t l = 0; l < nd_hierarchy.GetNumLevels(); l++)
  {
    nd_hierarchy.GetFESpaceAtLevel(l).GetEssentialTrueDofs(neumann_marker, essential[l]);
  }

  // K = (ε⁻¹ ∇ × ., ∇ × .) on every level with the essential DOFs eliminated, with the
  // geometry factors computed on the fly. The quadrature is that of the run, or the order
  // the curl-curl terms of ND(k + 1) need if higher: exact on affine cells (reference
  // degree 2 (k + 1) - 2 on tetrahedra, 2 (k + 1) on hexahedra and prisms), with 2 (q - 1)
  // more for a coordinate map of order q, as the cell integrals of the estimator (see
  // CellQuadrature). CG minimizes this approximation of η², the estimator evaluates the
  // final flux; without it ND(p + 1) of the RT(p) flux would be under-integrated on
  // hexahedra and prisms with the default quadrature of the run (QuadratureOrderJacobian
  // false).
  {
    const auto &eps = mat_op.GetPermittivityReal();
    mfem::DenseTensor eps_inv(eps.SizeI(), eps.SizeJ(), eps.SizeK());
    for (int k = 0; k < eps.SizeK(); k++)
    {
      mfem::DenseMatrix m(eps(k));
      m.Invert();
      eps_inv(k) = m;
    }
    MaterialPropertyCoefficient epsinv_func(mat_op.GetAttributeToMaterial(), eps_inv);
    BilinearForm k(nd_hierarchy.GetFinestFESpace());
    k.AddDomainIntegrator<CurlCurlIntegrator>(epsinv_func)
        .ComputeGeometryOnTheFly(
            [p_nd](const mfem::IsoparametricTransformation &T)
            {
              const bool simplex = (T.GetGeometryType() == mfem::Geometry::TETRAHEDRON);
              const int q = T.GetFE()->GetOrder();
              return std::max(fem::DefaultIntegrationOrder::Get(T),
                              2 * p_nd - (simplex ? 2 : 0) + 2 * (q - 1));
            });
    auto k_vec = k.Assemble(nd_hierarchy, false);
    auto K_mg = std::make_unique<MultigridOperator>(nd_hierarchy.GetNumLevels());
    for (std::size_t l = 0; l < nd_hierarchy.GetNumLevels(); l++)
    {
      auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]),
                                               nd_hierarchy.GetFESpaceAtLevel(l));
      K_l->SetEssentialTrueDofs(essential[l], Operator::DiagonalPolicy::DIAG_ONE);
      if (l == 0)
      {
        K_coarse = K_l.get();
      }
      K_mg->AddOperator(std::move(K_l));
    }
    K = std::move(K_mg);
  }

  // The solver of the magnetostatic curl-curl problem, configured as for a magnetostatic
  // simulation (see BaseKspSolver): CG, geometric multigrid over the p-levels with
  // Chebyshev smoothing and AMS for the singular curl-curl operator on the coarsest level
  // (with the smoothing sweeps of a singular coarse solve), without auxiliary space
  // smoothing (the gradients are in the kernel of K), with the multigrid settings of the
  // run. CG from zero with the stopping tests of the estimator, and its time counted as the
  // estimator's.
  auto &h1_hierarchy = broken ? h1_broken : h1_spaces;
  config::LinearSolverData solver = linear;
  solver.type = LinearSolver::AMS;
  solver.krylov_solver = KrylovSolver::CG;
  solver.ams_singular_op = 1;
  solver.mg_smooth_aux = 0;
  // Distance-one interpolation in the three auxiliary space hierarchies of AMS (cycle 14),
  // whose setup with extended+i interpolation exchanges with most processes on their coarse
  // levels: on the mesh of 4.2M unknowns of a model with an interior boundary at 144
  // processes, 26 MB per process of MPI buffers against 6 MB, the same CG iterations.
  // Direct interpolation divides by the diagonal and by sums of entries of one sign, which
  // do not cancel, unlike classical interpolation (which diverged on a refined mesh).
  solver.ams_aux_direct_interp = true;
  solver.initial_guess = 0;
  solver.tol = linear.estimator_tol;
  solver.abs_tol = std::numeric_limits<double>::epsilon();
  solver.max_it = linear.estimator_max_it;
  ksp = std::make_unique<KspSolver>(solver, MatrixSymmetry::SPD, 0, nd_hierarchy,
                                    &h1_hierarchy, false);
  // Every iterate gives a guaranteed bound, so the tolerance controls the bound itself: the
  // bound η of an iterate exceeds the converged one η* by η² - η*² = ‖D - D*‖², so η ≤ (1 +
  // tol) η* once the flux error ‖D - D*‖ is below tol' η, tol' = (2 tol + tol²)^{1/2} / (1
  // + tol), the least-squares test of CG with tol'. The test also stops CG once the
  // decrease of η² is below its rounding error (with tol = 0 too): further iterations
  // cannot lower the bound.
  const double tol = linear.estimator_error_tol;
  ls_tol = (tol > 0.0) ? std::sqrt(2.0 * tol + tol * tol) / (1.0 + tol) : 0.0;
  ksp->SetLeastSquaresTol(ls_tol);

  // The coarsest level is only applied through its assembled matrix (for AMS): release its
  // partial assembly once assembled, and return the memory before the smoothers are set up.
  if (nd_hierarchy.GetNumLevels() > 1)
  {
    K_coarse->ParallelAssemble();
    K_coarse->ReleaseLocalOperator();
    memory_reporting::ReleaseFreedMemory();
  }
  ksp->SetOperators(*K, *K);
  memory_reporting::ReleaseFreedMemory();
}

const mfem::DenseMatrix &CurlCorrectionSolver::GetCurlMatrix(int e) const
{
  const auto &rt = *flux.Get().GetFE(e);
  auto it = curl_matrix.find(rt.GetGeomType());
  if (it == curl_matrix.end())
  {
    it = curl_matrix
             .emplace(rt.GetGeomType(), ReferenceEmbedding(*ND().Get().GetFE(e), rt, true))
             .first;
  }
  return it->second;
}

CurlCorrectionSolver::~CurlCorrectionSolver()
{
  ksp.reset();
  if (shared)
  {
    shared->ReleaseDiscreteInterpolator();
  }
}

double CurlCorrectionSolver::AssembleRHS(const Vector &flux_l, const Vector &companion_l,
                                         Vector &b_l) const
{
  const auto &mesh = primal.GetParMesh();
  b_l.SetSize(ND().GetVSize());
  b_l.UseDevice(false);
  b_l = 0.0;
  double eta2 = 0.0;
  mfem::IsoparametricTransformation T;
  CellMap map;
  CellMaterial material;
  mfem::Vector ref(3), D(3), grad_ref(3), grad(3), R(3), eR(3), JtR(3), r, l;
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_trans;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &tab = quadrature.GetTables(e, T);
    const auto &C = GetCurlMatrix(e);
    const auto c = ElementValues(flux, e, flux_l),
               u = ElementValues(primal, e, companion_l);
    material.Update(mat_op, mesh.GetAttribute(e));
    r.SetSize(C.Height());
    r = 0.0;
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      if (q == 0 || !quadrature.IsAffine(e))
      {
        map.Update(T.GetPointMat(), tab.map[q]);
      }
      // R = ε⁻¹ (D₀ - D_v) = ε⁻¹ D₀ + ∇v, with D₀ = J D̂ / det J and ∇v = adj(J)ᵀ ∇̂v / det
      // J.
      const double weight = tab.rule->IntPoint(q).weight;
      tab.flux[q].MultTranspose(c, ref);
      map.J.Mult(ref, D);
      D /= map.det;
      tab.grad[q].MultTranspose(u, grad_ref);
      map.adjJ.MultTranspose(grad_ref, grad);
      grad /= map.det;
      material.eps_inv.Mult(D, R);
      R += grad;
      material.eps.Mult(R, eR);
      eta2 += weight * std::abs(map.det) * (R * eR);
      // (R, φⱼ) = sign(det J) ∫ (Jᵀ R) ⋅ φ̂ⱼ dx̂.
      map.J.MultTranspose(R, JtR);
      tab.flux[q].AddMult(JtR, r, (map.det > 0.0 ? 1.0 : -1.0) * weight);
    }
    // b = -Cᵀ r on the reference ND basis, to the L-vector with the DOF transformation.
    l.SetSize(C.Width());
    C.MultTranspose(r, l);
    l.Neg();
    ND().GetElementDofs(e, dofs, dof_trans);
    dof_trans.TransformDual(l);
    b_l.AddElementVector(dofs, l);
  }
  return eta2;
}

Vector CurlCorrectionSolver::Curl(const Vector &a_l) const
{
  // The curl is exact element by element: ∇ × A has the reference RT coefficients C a.
  const auto &mesh = primal.GetParMesh();
  Vector curl_l(flux.GetVSize());
  curl_l.UseDevice(false);
  curl_l = 0.0;
  mfem::Vector d;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const auto &C = GetCurlMatrix(e);
    const auto a = ElementValues(ND(), e, a_l);
    d.SetSize(C.Height());
    C.Mult(a, d);
    SetElementValues(flux, e, d, curl_l);
  }
  return curl_l;
}

double CurlCorrectionSolver::Bound(const Vector &flux_l, const Vector &companion_l,
                                   const double range[2]) const
{
  const auto terms = EvaluateBoundTerms(mat_op, flux, primal, quadrature, dirichlet_attr,
                                        flux_l, companion_l, 0.5 * (range[0] + range[1]),
                                        0.5 * (range[1] - range[0]));
  double sums[2] = {0.0, 0.0};  // η², b.
  for (std::size_t e = 0; e < terms.mismatch.size(); e++)
  {
    sums[0] += terms.mismatch[e];
    sums[1] += terms.balance[e];
  }
  Mpi::GlobalSum(2, sums, primal.GetComm());
  return BoundL(sums[0], sums[1]);
}

bool CurlCorrectionSolver::AddIfNotWorse(Vector &flux_l, const Vector &curl_l,
                                         const Vector &companion_l,
                                         const double range[2]) const
{
  Vector corrected(flux_l);
  corrected.UseDevice(false);
  corrected += curl_l;
  const double L0 = Bound(flux_l, companion_l, range),
               L1 = Bound(corrected, companion_l, range);
  if (L1 <= L0 * (1.0 + 1.0e-12))
  {
    flux_l = corrected;
    return true;
  }
  Mpi::Warning(primal.GetComm(),
               "The corrected flux has a larger bound ({:.3e}) than the uncorrected one "
               "({:.3e}): the flux is kept uncorrected!\n",
               L1, L0);
  return false;
}

CurlCorrectionSolver::Result CurlCorrectionSolver::Correct(Vector &flux_l,
                                                           const Vector &companion_l,
                                                           const double range[2],
                                                           double tol) const
{
  BlockTimer bt(Timer::SOLVE_ESTIMATOR);
  const MPI_Comm comm = primal.GetComm();
  const auto &P = *ND().GetProlongationMatrix();
  Vector b_l;
  double eta2 = AssembleRHS(flux_l, companion_l, b_l);
  Mpi::GlobalSum(1, &eta2, comm);
  Vector b(ND().GetTrueVSize()), x(ND().GetTrueVSize());
  b.UseDevice(true);
  x.UseDevice(true);
  P.MultTranspose(b_l, b);
  b_l.Destroy();
  b.SetSubVector(essential.back(), 0.0);
  x = 0.0;
  Result result;
  result.eta2.push_back(eta2);
  ksp->SetLeastSquaresNorm(eta2);
  if (tol >= 0.0)
  {
    ksp->SetLeastSquaresTol(tol);
  }
  ksp->Mult(b, x);
  if (tol >= 0.0)
  {
    ksp->SetLeastSquaresTol(ls_tol);
  }
  result.iterations = ksp->GetNumIterations();
  const auto decrease = ksp->GetLeastSquaresDecrease();
  for (std::size_t k = 1; k < decrease.size(); k++)
  {
    result.eta2.push_back(eta2 - decrease[k]);
  }
  result.norm = linalg::Norml2(comm, x);
  Vector a_l(ND().GetVSize());
  a_l.UseDevice(true);
  P.Mult(x, a_l);
  // The curl of the final iterate is added only if the assembled corrected flux has no
  // larger bound (evaluated as by the estimator), whatever CG did.
  result.accepted = AddIfNotWorse(flux_l, Curl(a_l), companion_l, range);
  return result;
}

EquilibratedFluxErrorEstimator::EquilibratedFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpace &h1_fespace,
    FiniteElementSpace &rt_fespace, FiniteElementSpace &nd_fespace,
    const mfem::Array<int> &dirichlet_attr, int flux_order_offset,
    const config::LinearSolverData &linear)
  : mat_op(mat_op), h1_fespace(h1_fespace), nd_fespace(&nd_fespace), linear(linear),
    dirichlet_attr(dirichlet_attr.begin(), dirichlet_attr.end()),
    hat_fec(1, h1_fespace.Dimension()), rt0_fec(0, h1_fespace.Dimension())
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  const MPI_Comm comm = h1_fespace.GetComm();
  auto &mesh = h1_fespace.GetMesh();
  MFEM_VERIFY(h1_fespace.Dimension() == 3,
              "The equilibrated error estimator requires a 3D mesh!");
  bool pyramids = mesh.Get().HasGeometry(mfem::Geometry::PYRAMID);
  Mpi::GlobalOr(1, &pyramids, comm);
  MFEM_VERIFY(!pyramids, "The equilibrated error estimator does not support pyramids!");
  MFEM_VERIFY(flux_order_offset == -1 || flux_order_offset == 0,
              "Unexpected flux order of the equilibrated error estimator!");

  // The flux space RT(k), k = p - 1 or p (the caller's RT(p - 1) if it is one), and the
  // views broken across the interior Dirichlet boundaries.
  const int k = h1_fespace.GetMaxElementOrder() + flux_order_offset;
  if (&rt_fespace.GetMesh() == &mesh && !rt_fespace.IsBroken() &&
      rt_fespace.GetVDim() == 1 &&
      dynamic_cast<const mfem::RT_FECollection *>(&rt_fespace.GetFEColl()) &&
      rt_fespace.GetFEColl().GetOrder() == k + 1)
  {
    flux_unbroken = &rt_fespace;
  }
  else
  {
    flux_fec = std::make_unique<mfem::RT_FECollection>(k, h1_fespace.Dimension());
    flux_owned = std::make_unique<FiniteElementSpace>(mesh, flux_fec.get());
    flux_unbroken = flux_owned.get();
  }
  const auto &sides = mesh.GetCrackSides(this->dirichlet_attr);
  split = sides.Any();
  Mpi::GlobalOr(1, &split, comm);
  if (split)
  {
    flux_broken = std::make_unique<FiniteElementSpace>(*flux_unbroken, sides);
  }
  quadrature = std::make_unique<equilibrated::CellQuadrature>(Flux(), h1_fespace);
  // The products with the order-one fields of the capacitive step (RT(0) and gradients of
  // order one, see AddErrorIndicators) have degree k + 2 with the flux space RT(k) and p
  // with the solution order p on affine cells.
  low_quadrature = std::make_unique<equilibrated::CellQuadrature>(
      Flux(), h1_fespace,
      std::max(Flux().GetFEColl().GetOrder() + 1, h1_fespace.GetMaxElementOrder()));

  const int bdr_attr_max =
      mesh.Get().bdr_attributes.Size() ? mesh.Get().bdr_attributes.Max() : 0;
  h1_fespace.Get().GetEssentialTrueDofs(mesh::AttrToMarker(bdr_attr_max, dirichlet_attr),
                                        essential);
}

EquilibratedFluxErrorEstimator::~EquilibratedFluxErrorEstimator() = default;

Vector EquilibratedFluxErrorEstimator::Embed(const FiniteElementSpace &lo,
                                             const Vector &x) const
{
  const auto &hi = Flux();
  Vector x_l(lo.GetVSize()), y(hi.GetVSize());
  x_l.UseDevice(false);
  y.UseDevice(false);
  lo.GetProlongationMatrix()->Mult(x, x_l);
  if (&lo == &hi)
  {
    return x_l;
  }
  y = 0.0;
  const auto &mesh = h1_fespace.GetParMesh();
  mfem::Vector b;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const auto &fl = *lo.Get().GetFE(e), &fh = *hi.Get().GetFE(e);
    auto it = embedding.find(fl.GetGeomType());
    if (it == embedding.end())
    {
      it = embedding.emplace(fl.GetGeomType(), ReferenceEmbedding(fl, fh, false)).first;
    }
    const auto a = ElementValues(lo, e, x_l);
    b.SetSize(it->second.Height());
    it->second.Mult(a, b);
    SetElementValues(hi, e, b, y);
  }
  return y;
}

void EquilibratedFluxErrorEstimator::HarmonicGram(
    const std::vector<Vector> &fluxes, const std::vector<Vector> &companions,
    const FiniteElementSpace &rt0, const FiniteElementSpace &hats,
    const std::vector<Vector> &f, const std::vector<Vector> &w, mfem::DenseMatrix &G,
    mfem::DenseMatrix &Q, std::vector<double> &mismatch) const
{
  // The products have degree k + 2 (RT(k) and RT(0)), p (gradients of order p - 1 and
  // RT(0)) and 2 at most on affine cells: the lower-order quadrature integrates them.
  const auto &mesh = h1_fespace.GetParMesh();
  const auto &cq = *low_quadrature;
  const int n = static_cast<int>(fluxes.size()), mh = static_cast<int>(f.size()),
            m = n + mh;
  std::vector<Vector> f_l(mh), w_l(mh);
  for (int s = 0; s < mh; s++)
  {
    f_l[s].SetSize(rt0.GetVSize());
    f_l[s].UseDevice(false);
    rt0.GetProlongationMatrix()->Mult(f[s], f_l[s]);
    w_l[s].SetSize(hats.GetVSize());
    w_l[s].UseDevice(false);
    hats.GetProlongationMatrix()->Mult(w[s], w_l[s]);
  }
  // The reference RT(0) shapes and hat gradients at the points of each rule.
  std::map<const mfem::IntegrationRule *,
           std::pair<std::vector<mfem::DenseMatrix>, std::vector<mfem::DenseMatrix>>>
      low;
  std::vector<double> sums(n * mh + mh * mh + n * mh + mh, 0.0);
  double *sGt = sums.data(), *sGh = sGt + n * mh, *sQ = sGh + mh * mh, *sM = sQ + n * mh;
  mfem::IsoparametricTransformation T;
  CellMap map;
  CellMaterial material;
  mfem::Vector ref(3), grad_ref(3), grad(3), r(3), Er(3);
  std::vector<mfem::Vector> c(n), u(n), h(mh), g(mh), D(n, mfem::Vector(3)),
      Dv(n, mfem::Vector(3)), F(mh, mfem::Vector(3)), EF(mh, mfem::Vector(3)),
      gw(mh, mfem::Vector(3));
  mfem::Array<int> fdofs, pdofs, rdofs, hdofs;
  mfem::DofTransformation fdt, pdt, rdt, hdt;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &tab = cq.GetTables(e, T);
    auto &[rt0_shape, hat_grad] = low[tab.rule];
    if (rt0_shape.empty())
    {
      const auto &fe = *rt0.Get().GetFE(e), &he = *hats.Get().GetFE(e);
      for (int q = 0; q < tab.rule->GetNPoints(); q++)
      {
        rt0_shape.emplace_back(fe.GetDof(), 3);
        fe.CalcVShape(tab.rule->IntPoint(q), rt0_shape.back());
        hat_grad.emplace_back(he.GetDof(), 3);
        he.CalcDShape(tab.rule->IntPoint(q), hat_grad.back());
      }
    }
    Flux().GetElementDofs(e, fdofs, fdt);
    h1_fespace.GetElementDofs(e, pdofs, pdt);
    for (int t = 0; t < n; t++)
    {
      c[t].SetSize(fdofs.Size());
      fluxes[t].GetSubVector(fdofs, c[t]);
      fdt.InvTransformPrimal(c[t]);
      u[t].SetSize(pdofs.Size());
      companions[t].GetSubVector(pdofs, u[t]);
      pdt.InvTransformPrimal(u[t]);
    }
    rt0.GetElementDofs(e, rdofs, rdt);
    hats.GetElementDofs(e, hdofs, hdt);
    for (int s = 0; s < mh; s++)
    {
      h[s].SetSize(rdofs.Size());
      f_l[s].GetSubVector(rdofs, h[s]);
      rdt.InvTransformPrimal(h[s]);
      g[s].SetSize(hdofs.Size());
      w_l[s].GetSubVector(hdofs, g[s]);
      hdt.InvTransformPrimal(g[s]);
    }
    material.Update(mat_op, mesh.GetAttribute(e));
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      if (q == 0 || !cq.IsAffine(e))
      {
        map.Update(T.GetPointMat(), tab.map[q]);
      }
      const double wq = tab.rule->IntPoint(q).weight * std::abs(map.det);
      for (int t = 0; t < n; t++)
      {
        tab.flux[q].MultTranspose(c[t], ref);
        map.J.Mult(ref, D[t]);
        D[t] /= map.det;
        tab.grad[q].MultTranspose(u[t], grad_ref);
        map.adjJ.MultTranspose(grad_ref, grad);
        grad /= map.det;
        material.eps.Mult(grad, Dv[t]);
        Dv[t].Neg();  // D_v = -ε ∇v.
      }
      for (int s = 0; s < mh; s++)
      {
        rt0_shape[q].MultTranspose(h[s], ref);
        map.J.Mult(ref, F[s]);
        F[s] /= map.det;
        material.eps_inv.Mult(F[s], EF[s]);
        hat_grad[q].MultTranspose(g[s], grad_ref);
        map.adjJ.MultTranspose(grad_ref, gw[s]);
        gw[s] /= map.det;
        material.eps.Mult(gw[s], r);
        r += F[s];  // f + ε ∇w.
        material.eps_inv.Mult(r, Er);
        sM[s] += wq * (r * Er);
      }
      for (int s = 0; s < mh; s++)
      {
        for (int t = 0; t < n; t++)
        {
          sGt[t * mh + s] += wq * (EF[s] * D[t]);
          sQ[t * mh + s] += wq * (EF[s] * Dv[t]);
        }
        for (int l = 0; l < mh; l++)
        {
          sGh[s * mh + l] += wq * (EF[s] * F[l]);
        }
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(sums.size()), sums.data(), h1_fespace.GetComm());
  mfem::DenseMatrix G_n(G), Q_n(Q);
  G.SetSize(m);
  Q.SetSize(n, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      G(i, j) = (i < n && j < n) ? G_n(i, j)
                : (i < n)        ? sGt[i * mh + (j - n)]
                : (j < n)        ? sGt[j * mh + (i - n)]
                                 : sGh[(i - n) * mh + (j - n)];
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      Q(i, j) = (j < n) ? Q_n(i, j) : sQ[i * mh + (j - n)];
    }
  }
  mismatch.assign(sM, sM + mh);
}

void EquilibratedFluxErrorEstimator::Gram(const std::vector<Vector> &fluxes,
                                          const std::vector<Vector> &companions,
                                          mfem::DenseMatrix &G, mfem::DenseMatrix &Q,
                                          mfem::DenseMatrix &Ch) const
{
  const auto &cq = *quadrature;
  const auto &mesh = h1_fespace.GetParMesh();
  const int m = static_cast<int>(fluxes.size()), n = static_cast<int>(companions.size());
  std::vector<double> sums(m * m + n * m + n * n, 0.0);
  double *sG = sums.data(), *sQ = sG + m * m, *sC = sQ + n * m;
  mfem::IsoparametricTransformation T;
  CellMap map;
  CellMaterial material;
  mfem::Vector ref(3), grad_ref(3), grad(3), tmp(3);
  std::vector<mfem::Vector> c(m), u(n), D(m, mfem::Vector(3)), Dv(n, mfem::Vector(3)),
      ED(m, mfem::Vector(3));
  mfem::Array<int> fdofs, pdofs;
  mfem::DofTransformation fdt, pdt;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &tab = cq.GetTables(e, T);
    Flux().GetElementDofs(e, fdofs, fdt);
    h1_fespace.GetElementDofs(e, pdofs, pdt);
    for (int t = 0; t < m; t++)
    {
      c[t].SetSize(fdofs.Size());
      fluxes[t].GetSubVector(fdofs, c[t]);
      fdt.InvTransformPrimal(c[t]);
    }
    for (int t = 0; t < n; t++)
    {
      u[t].SetSize(pdofs.Size());
      companions[t].GetSubVector(pdofs, u[t]);
      pdt.InvTransformPrimal(u[t]);
    }

    material.Update(mat_op, mesh.GetAttribute(e));
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      if (q == 0 || !cq.IsAffine(e))
      {
        map.Update(T.GetPointMat(), tab.map[q]);
      }
      const double w = tab.rule->IntPoint(q).weight * std::abs(map.det);
      for (int t = 0; t < m; t++)
      {
        tab.flux[q].MultTranspose(c[t], ref);
        map.J.Mult(ref, D[t]);
        D[t] /= map.det;
        material.eps_inv.Mult(D[t], ED[t]);
      }
      for (int t = 0; t < n; t++)
      {
        tab.grad[q].MultTranspose(u[t], grad_ref);
        map.adjJ.MultTranspose(grad_ref, grad);
        grad /= map.det;
        material.eps.Mult(grad, Dv[t]);
        Dv[t].Neg();  // D_v = -ε ∇v.
      }
      for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < m; j++)
        {
          sG[i * m + j] += w * (ED[i] * D[j]);
        }
      }

      for (int i = 0; i < n; i++)
      {
        material.eps_inv.Mult(Dv[i], tmp);
        for (int j = 0; j < m; j++)
        {
          sQ[i * m + j] += w * (ED[j] * Dv[i]);
        }
        for (int j = 0; j < n; j++)
        {
          sC[i * n + j] += w * (tmp * Dv[j]);
        }
      }
    }
  }
  Mpi::GlobalSum(static_cast<int>(sums.size()), sums.data(), h1_fespace.GetComm());
  G.SetSize(m);
  Q.SetSize(n, m);
  Ch.SetSize(n);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      G(i, j) = sG[i * m + j];
    }
  }
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      Q(i, j) = sQ[i * m + j];
    }
    for (int j = 0; j < n; j++)
    {
      Ch(i, j) = sC[i * n + j];
    }
  }
}

std::vector<EquilibratedFluxErrorEstimator::Summary>
EquilibratedFluxErrorEstimator::AddErrorIndicators(const std::vector<Vector> &V,
                                                   const std::vector<Vector> &G,
                                                   const std::vector<double> &energies,
                                                   ErrorIndicator &indicator)
{
  BlockTimer bt(Timer::ESTIMATION);
  const MPI_Comm comm = h1_fespace.GetComm();
  const int n = static_cast<int>(V.size());
  std::vector<Summary> summaries(n);

  // The companion potentials, which take the intended Dirichlet values at the essential
  // DOFs, and the Dirichlet ranges (min and -max per solve).
  std::vector<Vector> v(n), companions(n), liftings(n);
  std::vector<double> range(2 * n, std::numeric_limits<double>::infinity());
  for (int t = 0; t < n; t++)
  {
    v[t] = V[t];
    const auto *g = G[t].HostRead();
    auto *x = v[t].HostReadWrite();
    for (int i : essential)
    {
      x[i] = g[i];
      range[2 * t] = std::min(range[2 * t], g[i]);
      range[2 * t + 1] = std::min(range[2 * t + 1], -g[i]);
    }
    const auto &P = *h1_fespace.GetProlongationMatrix();
    companions[t].SetSize(h1_fespace.GetVSize());
    companions[t].UseDevice(false);
    P.Mult(v[t], companions[t]);
    Vector dv(v[t]);
    dv -= V[t];
    liftings[t].SetSize(h1_fespace.GetVSize());
    liftings[t].UseDevice(false);
    P.Mult(dv, liftings[t]);
  }
  Mpi::GlobalMin(2 * n, range.data(), comm);

  // 1. The RT(0) vertex patch flux, embedded into the flux space (RT(0) itself for k = 0),
  // of the solves and of the order-one potentials which complete their Dirichlet data to
  // the potentials constant on each connected component of the Dirichlet boundary (see
  // step 3).
  std::vector<Vector> fluxes(n);
  // The spaces of the patches, the order-one hats and RT(0) (the flux space for k = 0),
  // kept for the fields of the potentials (see steps 2 and 3).
  auto &mesh = h1_fespace.GetMesh();
  FiniteElementSpace hats(mesh, &hat_fec);
  std::unique_ptr<FiniteElementSpace> rt0, rt0_broken;
  const FiniteElementSpace *patch = &Flux();
  if (Flux().GetFEColl().GetOrder() > 1)
  {
    rt0 = std::make_unique<FiniteElementSpace>(mesh, &rt0_fec);
    if (split)
    {
      rt0_broken =
          std::make_unique<FiniteElementSpace>(*rt0, mesh.GetCrackSides(dirichlet_attr));
    }
    patch = rt0_broken ? rt0_broken.get() : rt0.get();
  }
  // The potentials (true DOFs of hats) and their RT(0) fields (true DOFs of the patch
  // space).
  std::vector<Vector> harmonic_w, harmonic_f;
  {
    // The order-one potentials, interpolated exactly into the primal space, have their
    // patch fluxes computed with those of the solves.
    std::vector<Vector> potentials(v);
    harmonic_w = HarmonicPotentials(hats, mat_op, dirichlet_attr, companions, linear,
                                    num_components);
    mfem::ParGridFunction w1(&hats.Get()), wp(&h1_fespace.Get());
    for (const auto &wk : harmonic_w)
    {
      w1.SetFromTrueDofs(wk);
      wp.ProjectGridFunction(w1);
      wp.GetTrueDofs(potentials.emplace_back(h1_fespace.GetTrueVSize()));
    }
    auto D0 =
        fem::VertexPatchFlux(h1_fespace, hats, *patch, mat_op, dirichlet_attr, potentials);
    for (int t = 0; t < n; t++)
    {
      fluxes[t] = Embed(*patch, D0[t]);
    }
    harmonic_f.assign(std::make_move_iterator(D0.begin() + n),
                      std::make_move_iterator(D0.end()));
  }
  num_harmonic = static_cast<int>(harmonic_f.size());
  memory_reporting::ReleaseFreedMemory();

  // The Gram matrices of the capacitive step over the fluxes of the solves and the fields
  // (see step 3), and whether the fields are corrected (then embedded into the flux space
  // and appended to the fluxes).
  const int m = n + num_harmonic;
  mfem::DenseMatrix Gm, Qm;
  num_harmonic_corrected = 0;

  // 2. The global curl correction of the fluxes of the solves, then of the fields of the
  // potentials as far as the fluxes of the solves need it. In the capacitive step the field
  // Fₛ enters the flux of solve i as αₛᵢ Fₛ, so its correction error eₛ (a curl) adds at
  // most |αₛᵢ| ‖eₛ‖ to the flux error of solve i, and ‖eₛ‖ ≤ ηₛ⁰ = ‖ε^{-1/2} (Fₛ + ε ∇wₛ)‖
  // for its potential wₛ (the correction does not depend on a potential constant on each
  // component). With the α of the capacitive step over the corrected fluxes of the solves
  // and the uncorrected fields, the fields are corrected (all of them, each until max_i
  // |αₛᵢ| ‖eₛ‖ / ηᵢ ≤ tol' / m, m fields, tol' the least-squares tolerance of
  // EstimatorErrorTol) only if one of them does not already meet that bound, so that
  // together they add at most tol' ηᵢ, the flux error allowed to the correction of the
  // solve itself; when EstimatorErrorTol is 0 they are corrected to convergence.
  {
    std::unique_ptr<CurlCorrectionSolver> curl;
    {
      BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
      curl = std::make_unique<CurlCorrectionSolver>(
          mat_op, h1_fespace, Flux(), dirichlet_attr, nd_fespace, linear, *quadrature);
    }
    for (int t = 0; t < n; t++)
    {
      const double r[2] = {range[2 * t], -range[2 * t + 1]};
      const auto result = curl->Correct(fluxes[t], companions[t], r);
      summaries[t].initial = std::sqrt(std::max(result.eta2.front(), 0.0));
      summaries[t].iterations = result.iterations;
    }
    Gram(fluxes, companions, Gm, Qm, C_h);
    if (num_harmonic > 0)
    {
      std::vector<double> mismatch;
      HarmonicGram(fluxes, companions, *patch, hats, harmonic_f, harmonic_w, Gm, Qm,
                   mismatch);
      mfem::DenseMatrix Qt, GiQt(m, n);
      Qt.Transpose(Qm);
      mfem::DenseMatrixInverse inverse(Gm);
      inverse.Mult(Qt, GiQt);  // Column i: α of solve i.
      std::vector<double> eta(n);
      for (int i = 0; i < n; i++)
      {
        double q = 0.0;
        for (int s = 0; s < m; s++)
        {
          q += Qm(i, s) * GiQt(s, i);
        }
        eta[i] = std::sqrt(std::max(C_h(i, i) - q, 0.0));
      }
      const double ls = curl->GetLeastSquaresTol();
      std::vector<double> tol(num_harmonic, 0.0);
      bool correct = false;
      for (int s = 0; s < num_harmonic; s++)
      {
        double weight = 0.0;  // max_i |αₛᵢ| ηₛ⁰ / ηᵢ.
        for (int i = 0; i < n; i++)
        {
          const double a = std::abs(GiQt(n + s, i)) * std::sqrt(std::max(mismatch[s], 0.0));
          weight =
              std::max(weight, (eta[i] > 0.0) ? a / eta[i]
                                              : std::numeric_limits<double>::infinity());
        }
        tol[s] = (ls > 0.0) ? ls / (num_harmonic * weight) : 0.0;
        correct = correct || (tol[s] < 1.0);
      }
      if (correct)
      {
        const double unit[2] = {0.0, 1.0};
        mfem::ParGridFunction w1(&hats.Get()), wp(&h1_fespace.Get());
        Vector w_t(h1_fespace.GetTrueVSize()), w_l(h1_fespace.GetVSize());
        w_l.UseDevice(false);
        for (int s = 0; s < num_harmonic; s++)
        {
          fluxes.push_back(Embed(*patch, harmonic_f[s]));
          if (tol[s] < 1.0)
          {
            w1.SetFromTrueDofs(harmonic_w[s]);
            wp.ProjectGridFunction(w1);
            wp.GetTrueDofs(w_t);
            h1_fespace.GetProlongationMatrix()->Mult(w_t, w_l);
            curl->Correct(fluxes.back(), w_l, unit, tol[s]);
            num_harmonic_corrected++;
          }
        }
        Gram(fluxes, companions, Gm, Qm, C_h);
      }
    }
  }
  memory_reporting::ReleaseFreedMemory();

  // 3. The capacitive step over the corrected fluxes of the solves and the fields of the
  // potentials, for the fluxes of the solves, in place (the fields not appended to the
  // fluxes are embedded one at a time).
  {
    mfem::DenseMatrix Qt(Qm), GiQt(m, n);
    Qt.Transpose();
    mfem::DenseMatrixInverse inverse(Gm);
    inverse.Mult(Qt, GiQt);  // Column i: α of solve i.
    C_low.SetSize(n);
    mfem::Mult(Qm, GiQt, C_low);
    const int mf = static_cast<int>(fluxes.size());
    std::vector<double> h(mf);
    std::vector<double *> f(mf);
    for (int t = 0; t < mf; t++)
    {
      f[t] = fluxes[t].HostReadWrite();
    }
    for (int i = 0; i < fluxes[0].Size(); i++)
    {
      for (int s = 0; s < mf; s++)
      {
        h[s] = f[s][i];
      }
      for (int t = 0; t < n; t++)
      {
        double value = 0.0;
        for (int s = 0; s < mf; s++)
        {
          value += GiQt(s, t) * h[s];
        }
        f[t][i] = value;
      }
    }
    fluxes.resize(n);
    for (int s = mf - n; s < num_harmonic; s++)
    {
      const auto field = Embed(*patch, harmonic_f[s]);
      for (int t = 0; t < n; t++)
      {
        fluxes[t].Add(GiQt(n + s, t), field);
      }
    }
  }

  // The bounds and indicators of the final fluxes.
  auto evaluated = Evaluate(fluxes, companions, liftings, range, energies, indicator);
  for (int t = 0; t < n; t++)
  {
    evaluated[t].initial = summaries[t].initial;
    evaluated[t].iterations = summaries[t].iterations;
  }
  return evaluated;
}

std::vector<EquilibratedFluxErrorEstimator::Summary>
EquilibratedFluxErrorEstimator::Evaluate(const std::vector<Vector> &fluxes,
                                         const std::vector<Vector> &companions,
                                         const std::vector<Vector> &liftings,
                                         const std::vector<double> &range,
                                         const std::vector<double> &energies,
                                         ErrorIndicator &indicator) const
{
  const MPI_Comm comm = h1_fespace.GetComm();
  const auto &mesh = h1_fespace.GetParMesh();
  const int n = static_cast<int>(fluxes.size()), ne = mesh.GetNE();
  std::vector<Summary> summaries(n);
  for (int t = 0; t < n; t++)
  {
    MFEM_VERIFY(std::isfinite(range[2 * t]) && std::isfinite(range[2 * t + 1]),
                "The equilibrated error estimator requires Dirichlet data!");
    const double center = 0.5 * (range[2 * t] - range[2 * t + 1]),
                 radius = -0.5 * (range[2 * t] + range[2 * t + 1]);

    const auto terms =
        EvaluateBoundTerms(mat_op, Flux(), h1_fespace, *quadrature, dirichlet_attr,
                           fluxes[t], companions[t], center, radius, &liftings[t]);
    const auto &mismatch = terms.mismatch, &balance = terms.balance, &lift = terms.lift;
    double sums[5] = {0.0, 0.0, 0.0, terms.div2, terms.neumann2};
    for (int e = 0; e < ne; e++)
    {
      sums[0] += mismatch[e];
      sums[1] += balance[e];
      sums[2] += lift[e];
    }
    Mpi::GlobalSum(5, sums, comm);
    auto &summary = summaries[t];
    summary.divergence = std::sqrt(sums[3]);
    summary.neumann = std::sqrt(sums[4]);
    summary.eta = std::sqrt(sums[0]);
    summary.balance = sums[1];
    summary.lifting = std::sqrt(sums[2]);
    const double L = BoundL(sums[0], sums[1]);
    const double A = L + summary.lifting;
    MFEM_VERIFY(std::isfinite(A) && std::isfinite(energies[t]) && energies[t] >= 0.0,
                "Nonfinite equilibrated error estimate or field energy!");
    summary.bound = A;
    summary.relative = RelativeBound(L, summary.lifting, energies[t]);

    // The indicators, q_K² = A (m_K / L + β_K² / β) with m_K = (L / η) η_K² + b_K, sum to
    // A², scaled to sum to the relative bound squared.
    const double weight = (summary.eta > 0.0) ? L / summary.eta : 0.0;
    const double scale = (A > 0.0) ? std::min(summary.relative, uncertified) / A : 0.0;
    Vector estimates(ne);
    estimates.UseDevice(false);
    for (int e = 0; e < ne; e++)
    {
      double q2 = 0.0;
      if (L > 0.0)
      {
        q2 += A * (weight * mismatch[e] + balance[e]) / L;
      }
      if (summary.lifting > 0.0)
      {
        q2 += A * lift[e] / summary.lifting;
      }
      estimates[e] = std::sqrt(q2) * scale;
    }
    indicator.AddIndicator(estimates);
  }
  return summaries;
}

double EquilibratedFluxErrorEstimator::RelativeBound(double L, double lifting,
                                                     double energy)
{
  const double A = L + lifting;
  if (A <= 0.0)
  {
    return 0.0;
  }
  const double norm = std::max(std::sqrt(2.0 * energy) - lifting, 0.0);
  const double lower2 = norm * norm - L * L;  // ≤ |||u|||².
  return (lower2 > 0.0) ? A / std::sqrt(lower2) : std::numeric_limits<double>::infinity();
}

}  // namespace palace
