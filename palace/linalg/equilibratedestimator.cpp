// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "equilibratedestimator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include "fem/bilinearform.hpp"
#include "fem/errorindicator.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/vertexpatchflux.hpp"
#include "linalg/ams.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
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

}  // namespace

namespace equilibrated
{

CellQuadrature::CellQuadrature(const FiniteElementSpace &flux,
                               const FiniteElementSpace &primal)
  : flux(flux), primal(primal)
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
  return 2 * primal.GetMaxElementOrder() + 2 + 2 * (q - 1);
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
  : mat_op(mat_op), primal(primal), flux(flux), quadrature(quadrature)
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

  // The p-levels of the correction space as for the field solve, with the order-one H1
  // space for AMS on the coarsest level, broken across the interior Dirichlet boundaries
  // when the flux space is.
  nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      p_nd, 3, linear.mg_max_levels, linear.mg_coarsening, false);
  h1_fec = std::make_unique<mfem::H1_FECollection>(nd_fecs.front()->GetOrder(), 3);
  const bool share = shared_nd && &shared_nd->GetMesh() == &primal.GetMesh() &&
                     !shared_nd->IsBroken() && shared_nd->GetVDim() == 1 &&
                     dynamic_cast<const mfem::ND_FECollection *>(&shared_nd->GetFEColl()) &&
                     shared_nd->GetFEColl().GetOrder() == p_nd;
  for (std::size_t l = 0; l < nd_fecs.size(); l++)
  {
    if (share && l + 1 == nd_fecs.size())
    {
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
  // geometry factors computed on the fly. The quadrature is that of the run (exact on
  // affine tetrahedra); CG minimizes its approximation of η², the estimator evaluates the
  // final flux.
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
    k.AddDomainIntegrator<CurlCurlIntegrator>(epsinv_func).ComputeGeometryOnTheFly();
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

  // The preconditioner of the magnetostatic solve with the multigrid settings of the run:
  // geometric multigrid over the p-levels with Chebyshev smoothing and AMS for the singular
  // curl-curl operator on the coarsest level, without auxiliary space smoothing (the
  // gradients are in the kernel of K, where the right-hand side vanishes). CG from zero
  // with the stopping tests of the estimator.
  auto &h1_hierarchy = broken ? h1_broken : h1_spaces;
  const bool multigrid = (nd_hierarchy.GetNumLevels() > 1);
  auto ams = std::make_unique<MfemWrapperSolver<Operator>>(std::make_unique<HypreAmsSolver>(
      nd_hierarchy.GetFESpaceAtLevel(0), h1_hierarchy.GetFESpaceAtLevel(0),
      multigrid ? linear.ams_max_it : linear.mg_cycle_it, linear.mg_smooth_it,
      linear.ams_vector_interp, true, linear.amg_agg_coarsen,
      linear.ams_agglomeration_target, 0));
  std::unique_ptr<Solver<Operator>> pc;
  if (multigrid)
  {
    pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        comm, std::move(ams), nd_hierarchy.GetProlongationOperators(), nullptr,
        linear.mg_cycle_it, linear.mg_smooth_it, linear.mg_smooth_order,
        linear.mg_smooth_sf_max, linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  }
  else
  {
    pc = std::move(ams);
  }
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(linear.estimator_tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(linear.estimator_max_it);
  // Every iterate gives a guaranteed bound, so the tolerance controls the bound itself: the
  // bound η of an iterate exceeds the converged one η* by η² - η*² = ‖D - D*‖², so η ≤ (1 +
  // tol) η* once the flux error ‖D - D*‖ is below tol' η, tol' = (2 tol + tol²)^{1/2} / (1
  // + tol), the least-squares test of CG with tol'.
  const double tol = linear.estimator_error_tol;
  least_squares = (tol > 0.0);
  pcg->SetLeastSquaresTol(least_squares ? std::sqrt(2.0 * tol + tol * tol) / (1.0 + tol)
                                        : -1.0);
  ksp = std::make_unique<KspSolver>(std::move(pcg), std::move(pc));

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

void CurlCorrectionSolver::AddCurl(const Vector &a_l, Vector &flux_l) const
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
  flux_l += curl_l;
}

CurlCorrectionSolver::Result CurlCorrectionSolver::Correct(Vector &flux_l,
                                                           const Vector &companion_l) const
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
  if (least_squares)
  {
    ksp->SetLeastSquaresNorm(eta2);
  }
  ksp->Mult(b, x);
  result.iterations = ksp->GetNumIterations();
  if (least_squares)
  {
    const auto decrease = ksp->GetLeastSquaresDecrease();
    for (std::size_t k = 1; k < decrease.size(); k++)
    {
      result.eta2.push_back(eta2 - decrease[k]);
    }
  }
  Vector a_l(ND().GetVSize());
  a_l.UseDevice(true);
  P.Mult(x, a_l);
  AddCurl(a_l, flux_l);
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

void EquilibratedFluxErrorEstimator::Gram(const std::vector<Vector> &fluxes,
                                          const std::vector<Vector> &companions,
                                          mfem::DenseMatrix &G, mfem::DenseMatrix &Q,
                                          mfem::DenseMatrix &Ch) const
{
  const auto &mesh = h1_fespace.GetParMesh();
  const int n = static_cast<int>(fluxes.size());
  std::vector<double> sums(3 * n * n, 0.0);
  mfem::IsoparametricTransformation T;
  CellMap map;
  CellMaterial material;
  mfem::Vector ref(3), grad_ref(3), grad(3), tmp(3);
  std::vector<mfem::Vector> c(n), u(n), D(n, mfem::Vector(3)), Dv(n, mfem::Vector(3)),
      ED(n, mfem::Vector(3));
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementTransformation(e, &T);
    const auto &tab = quadrature->GetTables(e, T);
    for (int t = 0; t < n; t++)
    {
      c[t] = ElementValues(Flux(), e, fluxes[t]);
      u[t] = ElementValues(h1_fespace, e, companions[t]);
    }
    material.Update(mat_op, mesh.GetAttribute(e));
    for (int q = 0; q < tab.rule->GetNPoints(); q++)
    {
      if (q == 0 || !quadrature->IsAffine(e))
      {
        map.Update(T.GetPointMat(), tab.map[q]);
      }
      const double w = tab.rule->IntPoint(q).weight * std::abs(map.det);
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
        material.eps_inv.Mult(D[t], ED[t]);
      }
      for (int i = 0; i < n; i++)
      {
        material.eps_inv.Mult(Dv[i], tmp);
        for (int j = 0; j < n; j++)
        {
          sums[i * n + j] += w * (ED[i] * D[j]);
          sums[n * n + i * n + j] += w * (ED[j] * Dv[i]);
          sums[2 * n * n + i * n + j] += w * (tmp * Dv[j]);
        }
      }
    }
  }
  Mpi::GlobalSum(3 * n * n, sums.data(), h1_fespace.GetComm());
  G.SetSize(n);
  Q.SetSize(n);
  Ch.SetSize(n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      G(i, j) = sums[i * n + j];
      Q(i, j) = sums[n * n + i * n + j];
      Ch(i, j) = sums[2 * n * n + i * n + j];
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

  // 1. The RT(0) vertex patch flux, embedded into the flux space (RT(0) itself for k = 0).
  std::vector<Vector> fluxes(n);
  {
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
    const auto D0 =
        fem::VertexPatchFlux(h1_fespace, hats, *patch, mat_op, dirichlet_attr, v);
    for (int t = 0; t < n; t++)
    {
      fluxes[t] = Embed(*patch, D0[t]);
    }
  }
  memory_reporting::ReleaseFreedMemory();

  // 2. The global curl correction of each flux.
  {
    std::unique_ptr<CurlCorrectionSolver> curl;
    {
      BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
      curl = std::make_unique<CurlCorrectionSolver>(
          mat_op, h1_fespace, Flux(), dirichlet_attr, nd_fespace, linear, *quadrature);
    }
    for (int t = 0; t < n; t++)
    {
      const auto result = curl->Correct(fluxes[t], companions[t]);
      summaries[t].initial = std::sqrt(std::max(result.eta2.front(), 0.0));
      summaries[t].iterations = result.iterations;
    }
  }
  memory_reporting::ReleaseFreedMemory();

  // 3. The capacitive step over the corrected fluxes of the terminal solves, in place.
  {
    mfem::DenseMatrix Gm, Qm;
    Gram(fluxes, companions, Gm, Qm, C_h);
    mfem::DenseMatrix Qt(Qm), GiQt(n), alpha(n);
    Qt.Transpose();
    mfem::DenseMatrixInverse inverse(Gm);
    inverse.Mult(Qt, GiQt);  // Column i: α of solve i.
    C_low.SetSize(n);
    mfem::Mult(Qm, GiQt, C_low);
    std::vector<double> h(n);
    std::vector<double *> f(n);
    for (int t = 0; t < n; t++)
    {
      f[t] = fluxes[t].HostReadWrite();
    }
    for (int i = 0; i < fluxes[0].Size(); i++)
    {
      for (int s = 0; s < n; s++)
      {
        h[s] = f[s][i];
      }
      for (int t = 0; t < n; t++)
      {
        double value = 0.0;
        for (int s = 0; s < n; s++)
        {
          value += GiQt(s, t) * h[s];
        }
        f[t][i] = value;
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
  const std::set<int> dirichlet(dirichlet_attr.begin(), dirichlet_attr.end());
  const auto face_to_bdr = mesh.GetFaceToBdrElMap();
  std::vector<Summary> summaries(n);
  for (int t = 0; t < n; t++)
  {
    MFEM_VERIFY(std::isfinite(range[2 * t]) && std::isfinite(range[2 * t + 1]),
                "The equilibrated error estimator requires Dirichlet data!");
    const double center = 0.5 * (range[2 * t] - range[2 * t + 1]),
                 radius = -0.5 * (range[2 * t] + range[2 * t + 1]);

    // Cell integrals: the mismatch η_K², the balance term b_K (from the divergence, the
    // Neumann faces added below) and the lifting β_K².
    std::vector<double> mismatch(ne, 0.0), balance(ne, 0.0), lift(ne, 0.0);
    double div_total = 0.0, neumann_total = 0.0;
    mfem::IsoparametricTransformation T;
    CellMap map;
    CellMaterial material;
    mfem::Vector ref(3), D(3), grad_ref(3), grad(3), Dv(3), r(3), Er(3);
    for (int e = 0; e < ne; e++)
    {
      mesh.GetElementTransformation(e, &T);
      const auto &tab = quadrature->GetTables(e, T);
      const auto c = ElementValues(Flux(), e, fluxes[t]);
      const auto u = ElementValues(h1_fespace, e, companions[t]);
      const auto du = ElementValues(h1_fespace, e, liftings[t]);
      material.Update(mat_op, mesh.GetAttribute(e));
      double eta2 = 0.0, div2 = 0.0, volume = 0.0, offset2 = 0.0, beta2 = 0.0;
      for (int q = 0; q < tab.rule->GetNPoints(); q++)
      {
        if (q == 0 || !quadrature->IsAffine(e))
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
        tab.grad[q].MultTranspose(du, grad_ref);
        map.adjJ.MultTranspose(grad_ref, grad);
        grad /= map.det;
        material.eps.Mult(grad, Dv);
        beta2 += w * (grad * Dv);
      }
      mismatch[e] = eta2;
      div_total += div2;
      balance[e] = std::sqrt(div2) * (radius * std::sqrt(volume) + std::sqrt(offset2));
      lift[e] = beta2;
    }

    // Neumann faces (exterior faces without a Dirichlet attribute): the normal trace of the
    // flux is zero there up to rounding, its actual value enters the balance term.
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
      const auto &fe = *Flux().Get().GetFE(e), &pe = *h1_fespace.Get().GetFE(e);
      const auto c = ElementValues(Flux(), e, fluxes[t]);
      const auto u = ElementValues(h1_fespace, e, companions[t]);
      shape.SetSize(fe.GetDof(), 3);
      pshape.SetSize(pe.GetDof());
      double flux2 = 0.0, area = 0.0, offset2 = 0.0;
      const auto &rule = mfem::IntRules.Get(FT.GetGeometryType(), quadrature->Order());
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
      neumann_total += flux2;
      balance[e] += std::sqrt(flux2) * (radius * std::sqrt(area) + std::sqrt(offset2));
    }

    double sums[5] = {0.0, 0.0, 0.0, div_total, neumann_total};
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
    const double L = 0.5 * (summary.eta + std::sqrt(sums[0] + 4.0 * sums[1]));
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
