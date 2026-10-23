// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurlsolver.hpp"
#include "surfacefluxoperator.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_set>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"

namespace palace
{

namespace
{

// Build the London drive a_h as a 3D curl-free cohomology generator carrying the hole
// fluxoid ∮_∂hole a_h·dl = Φ, about the flux line L along the loop Direction n through the
// hole centroid. L threads the hole opening off the film Σ, so a_h|Σ is curl-free with
// circulation Φ; its gradient part is absorbable by A → A + ∇χ, so the extracted inductance
// depends only on the cohomology class (Φ) and is gauge-invariant.
//
// The cut is the half-plane bounded by L in the plane ⊥ n, located via a right-handed
// in-plane basis (e1, e2, n), so the construction handles any planar hole orientation. The
// start is the lowest-order (Whitney) cut cochain c0: ±Φ on each edge crossing the cut, 0
// elsewhere, written on the true (unconstrained) edge DOFs. On a conforming mesh c0 is
// already curl-free on Σ. With hanging edges the slave edges are interpolated from the
// masters and the film picks up a residual curl (0.07–0.10 Φ, biasing L), so the curl is
// removed in the constrained space: with C the discrete curl from the ND_1 true DOFs to the
// film's RT_0 face true DOFs, a_h = c0 − Cᵀψ with (C Cᵀ) ψ = C c0 has C a_h = 0 exactly on
// any mesh (conforming, hanging edges, after rebalance), and the same cohomology class as
// c0. a_h is then projected into the order-p ND space (ND_1 ⊂ ND_p). The cut shape and
// overall sign are irrelevant to L: the downstream normalization rescales a_h so cᵀa_h = Φ.
Vector BuildCutCohomologyGenerator(const SurfaceFluxData &flux_data,
                                   const mfem::ParFiniteElementSpace &ndp_fespace,
                                   const Mesh &mesh, CohomologyGeneratorCache *, int)
{
  auto &pmesh = const_cast<mfem::ParMesh &>(mesh.Get());
  MPI_Comm comm = pmesh.GetComm();
  const int sdim = pmesh.SpaceDimension();
  MFEM_VERIFY(sdim == 3, "London cut cohomology generator requires a 3D mesh!");

  // Hole axis n = normalized Direction. Build a right-handed in-plane basis (e1, e2, n), e1
  // seeded from the Cartesian axis least aligned with n (avoids a degenerate cross
  // product); for n = ±ẑ this gives e1 = x̂, e2 = ±ŷ, reproducing the xy construction.
  MFEM_VERIFY(flux_data.direction.size() >= 3, "London flux loop needs a 3D Direction!");
  double n[3] = {flux_data.direction[0], flux_data.direction[1], flux_data.direction[2]};
  const double nrm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  MFEM_VERIFY(nrm > 1.0e-12, "London flux loop Direction must be nonzero!");
  n[0] /= nrm;
  n[1] /= nrm;
  n[2] /= nrm;
  double e1[3], e2[3];
  {
    int imin = 0;
    if (std::abs(n[1]) < std::abs(n[imin]))
    {
      imin = 1;
    }
    if (std::abs(n[2]) < std::abs(n[imin]))
    {
      imin = 2;
    }
    double seed[3] = {0.0, 0.0, 0.0};
    seed[imin] = 1.0;
    const double sn = seed[0] * n[0] + seed[1] * n[1] + seed[2] * n[2];
    for (int i = 0; i < 3; i++)
    {
      e1[i] = seed[i] - sn * n[i];
    }
    const double e1n = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (int i = 0; i < 3; i++)
    {
      e1[i] /= e1n;
    }
    e2[0] = n[1] * e1[2] - n[2] * e1[1];
    e2[1] = n[2] * e1[0] - n[0] * e1[2];
    e2[2] = n[0] * e1[1] - n[1] * e1[0];
  }

  // Total imposed fluxoid Φ.
  double phi = 0.0;
  for (double f : flux_data.flux_amounts)
  {
    phi += f;
  }

  // Hole centroid c: average of hole-boundary vertex coordinates. The flux line L = c + s·n
  // passes through the hole opening (empty region), never the film.
  std::unordered_set<int> hole_attrs(flux_data.hole_attributes.begin(),
                                     flux_data.hole_attributes.end());
  double csum[3] = {0.0, 0.0, 0.0};
  double cnt = 0.0;
  mfem::Array<int> bverts;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    if (!hole_attrs.count(pmesh.GetBdrAttribute(be)))
    {
      continue;
    }
    pmesh.GetBdrElementVertices(be, bverts);
    for (int v : bverts)
    {
      const double *x = pmesh.GetVertex(v);
      csum[0] += x[0];
      csum[1] += x[1];
      csum[2] += x[2];
      cnt += 1.0;
    }
  }
  Mpi::GlobalSum(3, csum, comm);
  Mpi::GlobalSum(1, &cnt, comm);
  MFEM_VERIFY(cnt > 0.0, "No hole boundary elements found for London cut generator!");
  double c[3] = {csum[0] / cnt, csum[1] / cnt, csum[2] / cnt};

  // In-plane coordinates (s, t) of a point about L: components along (e1, e2).
  auto inplane = [&](const double *x, double &s, double &t)
  {
    const double dx = x[0] - c[0], dy = x[1] - c[1], dz = x[2] - c[2];
    s = dx * e1[0] + dy * e1[1] + dz * e1[2];
    t = dx * e2[0] + dy * e2[1] + dz * e2[2];
  };

  {
    // θ is undefined on L. Measure the hole's radial extent and the closest approach of a
    // hole vertex to L, then shift L off any vertex sitting on it — any interior point of
    // the hole represents the same cohomology class.
    double r_min = mfem::infinity(), r_max = 0.0;
    for (int be = 0; be < pmesh.GetNBE(); be++)
    {
      if (!hole_attrs.count(pmesh.GetBdrAttribute(be)))
      {
        continue;
      }
      pmesh.GetBdrElementVertices(be, bverts);
      for (int v : bverts)
      {
        const double *x = pmesh.GetVertex(v);
        double s, t;
        inplane(x, s, t);
        const double r = std::hypot(s, t);
        r_min = std::min(r_min, r);
        r_max = std::max(r_max, r);
      }
    }
    Mpi::GlobalMin(1, &r_min, comm);
    Mpi::GlobalMax(1, &r_max, comm);
    if (r_min < 1.0e-6 * r_max)
    {
      for (int i = 0; i < 3; i++)
      {
        c[i] += r_max * (3.7e-3 * e1[i] + 2.3e-3 * e2[i]);
      }
      Mpi::Print(" London a_h: flux line met a hole vertex, offset to ({:.6e}, {:.6e}, "
                 "{:.6e})\n",
                 c[0], c[1], c[2]);
    }
  }

  // Film marker Σ (both faces of a two-sided film share the attribute).
  const int bmax = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> film_marker(bmax);
  film_marker = 0;
  for (int a : flux_data.film_attributes)
  {
    if (a >= 1 && a <= bmax)
    {
      film_marker[a - 1] = 1;
    }
  }

  // Lowest-order (Whitney) cut cochain on an ND_1 space over the full 3D mesh. Each rank
  // sets its local edge DOFs from vertex COORDINATES (identical on all ranks for a shared
  // edge) in the canonical GetEdgeVertices orientation (ev0→ev1), so the owner's true-DOF
  // value is correct and GetTrueDofs needs no cross-rank sign reconciliation. An edge
  // crosses the plane s = 0 when its vertices lie on different sides (s = 0 counts as the
  // positive side), at t = t0 + (t1 − t0)·s0/(s0 − s1); the cut is the half t > 0. Around
  // any face not pierced by L the signed crossings cancel, so c0 is curl-free on Σ wherever
  // the mesh is conforming.
  mfem::ND_FECollection nd1_fec(1, sdim);
  mfem::ParFiniteElementSpace nd1_fespace(&pmesh, &nd1_fec);
  mfem::ParGridFunction ah1(&nd1_fespace);
  ah1.UseDevice(false);
  ah1 = 0.0;
  {
    mfem::Array<int> ev, edofs;
    for (int e = 0; e < pmesh.GetNEdges(); e++)
    {
      pmesh.GetEdgeVertices(e, ev);
      double s0, t0, s1, t1;
      inplane(pmesh.GetVertex(ev[0]), s0, t0);
      inplane(pmesh.GetVertex(ev[1]), s1, t1);
      const bool side0 = (s0 >= 0.0), side1 = (s1 >= 0.0);
      if (side0 == side1 || t0 + (t1 - t0) * s0 / (s0 - s1) <= 0.0)
      {
        continue;
      }
      nd1_fespace.GetEdgeDofs(e, edofs);
      ah1(edofs[0]) = side1 ? phi : -phi;
    }
  }
  Vector c0(nd1_fespace.GetTrueVSize());
  c0.UseDevice(false);
  ah1.GetTrueDofs(c0);

  // Remove the curl of c0 on Σ in the constrained space. C is the lowest-order discrete
  // curl (ND_1 → RT_0) on true DOFs, so hanging edges and faces are constrained, with rows
  // restricted to the film faces (for a hanging film face its master, whose RT DOF carries
  // the curl of all its slaves). K = C Cᵀ has one unknown per film face and is symmetric
  // positive semidefinite: definite on an open film, and on a closed one (a two-sided film
  // seamed at both rims) C c0 ∈ range(K) still, so Jacobi-PCG converges and a_h = c0 − Cᵀψ
  // is unique either way. K is tiny next to the 3D system, so it does not warrant an AMG
  // setup.
  mfem::RT_FECollection rt0_fec(0, sdim);
  mfem::ParFiniteElementSpace rt0_fespace(&pmesh, &rt0_fec);
  std::unique_ptr<mfem::HypreParMatrix> C, Ct, K;
  {
    mfem::ParDiscreteLinearOperator curl(&nd1_fespace, &rt0_fespace);
    curl.AddDomainInterpolator(new mfem::CurlInterpolator);
    curl.Assemble();
    curl.Finalize();
    // The lowest-order curl is the face-edge incidence (entries 0, ±1), but MFEM evaluates
    // it through a numerically inverted reference basis; round it back so that C c0
    // vanishes exactly on conforming faces.
    auto *data = curl.SpMat().HostReadWriteData();
    for (int i = 0; i < curl.SpMat().NumNonZeroElems(); i++)
    {
      data[i] = std::round(data[i]);
    }
    C.reset(curl.ParallelAssemble());
  }
  mfem::Array<int> film_faces;
  rt0_fespace.GetEssentialTrueDofs(film_marker, film_faces);
  {
    HYPRE_BigInt n_film = film_faces.Size();
    Mpi::GlobalSum(1, &n_film, comm);
    MFEM_VERIFY(n_film > 0, "No film faces found for London cut generator!");
    Ct.reset(C->Transpose());
    std::unique_ptr<mfem::HypreParMatrix> CCt(mfem::ParMult(C.get(), Ct.get(), true));
    K.reset(CCt->ExtractSubmatrix(film_faces));
    // GetDiag (the Jacobi diagonal) reads the first entry of each row.
    hypre_CSRMatrixReorder(hypre_ParCSRMatrixDiag(static_cast<hypre_ParCSRMatrix *>(*K)));
  }
  Vector curl_c0(rt0_fespace.GetTrueVSize()), psi(rt0_fespace.GetTrueVSize());
  Vector r(film_faces.Size()), psi_film(film_faces.Size());
  curl_c0.UseDevice(false);
  psi.UseDevice(false);
  r.UseDevice(false);
  psi_film.UseDevice(false);
  C->Mult(c0, curl_c0);
  curl_c0.GetSubVector(film_faces, r);
  const double curl_before = linalg::Norml2(comm, r);
  int its = 0;
  if (curl_before > 0.0)
  {
    JacobiSmoother<Operator> jacobi(comm);
    jacobi.SetOperator(*K);
    CgSolver<Operator> pcg(comm, 0);
    pcg.SetOperator(*K);
    pcg.SetPreconditioner(jacobi);
    pcg.SetInitialGuess(false);
    pcg.SetRelTol(1.0e-12);
    pcg.SetAbsTol(0.0);
    pcg.SetMaxIter(10000);
    pcg.Mult(r, psi_film);
    MFEM_VERIFY(
        pcg.GetConverged(),
        "London cut generator curl projection did not converge, norm(Ax-b)/norm(b) = "
            << pcg.GetFinalRes() / pcg.GetInitialRes() << "!");
    its = pcg.GetNumIterations();
    psi = 0.0;
    psi.SetSubVector(film_faces, psi_film);
    Ct->AddMult(psi, c0, -1.0);
    C->Mult(c0, curl_c0);
    curl_c0.GetSubVector(film_faces, r);
  }
  Mpi::Print(" London a_h: film curl ‖C a_h‖ = {:.3e} → {:.3e} ({:d} PCG iterations)\n",
             curl_before, linalg::Norml2(comm, r), its);
  ah1.SetFromTrueDofs(c0);

  // Project the ND_1 cut field exactly into the order-p ND space. The per-element field is
  // a degree-1 polynomial, integrated exactly by the ND_p DOF functionals, so ND_1 ⊂ ND_p
  // is reproduced and tangential continuity (hence curl-free-on-Σ, circulation Φ) is
  // preserved.
  mfem::VectorGridFunctionCoefficient ah1_coeff(&ah1);
  mfem::ParGridFunction ahp(const_cast<mfem::ParFiniteElementSpace *>(&ndp_fespace));
  ahp.UseDevice(false);
  ahp = 0.0;
  ahp.ProjectCoefficient(ah1_coeff);

  Vector result(ndp_fespace.GetTrueVSize());
  result.UseDevice(true);
  ahp.GetTrueDofs(result);

  // Restrict a_h to this loop's film DOFs. It is consumed only through M_sheet and c (both
  // on the film), so this is a no-op for a single film but stops the cut line from driving
  // a spurious vortex in another Superconductor sheet the line threads.
  {
    mfem::Array<int> film_tdofs;
    const_cast<mfem::ParFiniteElementSpace &>(ndp_fespace)
        .GetEssentialTrueDofs(film_marker, film_tdofs);
    Vector keep(result.Size());
    keep.UseDevice(true);
    keep = 0.0;
    auto *k = keep.HostReadWrite();
    const auto *r = result.HostRead();
    for (int i = 0; i < film_tdofs.Size(); i++)
    {
      k[film_tdofs[i]] = r[film_tdofs[i]];
    }
    result = keep;
  }
  return result;
}

}  // namespace

CohomologyGeneratorCache::CohomologyGeneratorCache() = default;

CohomologyGeneratorCache::~CohomologyGeneratorCache() = default;

const mfem::ParGridFunction *CohomologyGeneratorCache::Find(int key,
                                                            const mfem::ParMesh &mesh) const
{
  auto it = entries.find(key);
  return (it != entries.end() && it->second.mesh == &mesh &&
          it->second.sequence == mesh.GetSequence())
             ? it->second.a.get()
             : nullptr;
}

void CohomologyGeneratorCache::Store(int key, const mfem::ParGridFunction &a)
{
  auto &mesh = *a.ParFESpace()->GetParMesh();
  const int sdim = mesh.SpaceDimension();
  auto &entry = entries[key];
  entry.a.reset();
  entry.fespace.reset();
  entry.mesh = &mesh;
  entry.sequence = mesh.GetSequence();
  entry.fec = std::make_unique<mfem::L2_FECollection>(1, mesh.Dimension());
  entry.fespace =
      std::make_unique<mfem::ParFiniteElementSpace>(&mesh, entry.fec.get(), sdim);
  entry.a = std::make_unique<mfem::ParGridFunction>(entry.fespace.get());
  mfem::VectorGridFunctionCoefficient a_coeff(&a);
  entry.a->ProjectCoefficient(a_coeff);
}

void CohomologyGeneratorCache::Update(mfem::ParMesh &mesh)
{
  // MFEM transfers a space across exactly one mesh modification; conformal meshes rebuild
  // the exact generator and do not need it.
  std::erase_if(entries,
                [&mesh](const auto &kv)
                {
                  const auto &entry = kv.second;
                  return entry.mesh != &mesh || !mesh.Nonconforming() ||
                         mesh.GetSequence() > entry.sequence + 1;
                });
  for (auto &[key, entry] : entries)
  {
    if (mesh.GetSequence() == entry.sequence + 1)
    {
      entry.fespace->Update();
      entry.a->Update();
      entry.sequence = mesh.GetSequence();
    }
  }
}

Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                               CohomologyGeneratorCache *cache, int key)
{
  Vector result;
  SolveSurfaceCurlProblem(flux_data, mesh, nd_fespace, post_op, result, cache, key);
  return result;
}

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                             const FiniteElementSpace &nd_fespace,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result, CohomologyGeneratorCache *cache, int key)
{
  // London drive: build a_h directly on the 3D ND space as a curl-free cohomology generator
  // via a topological cut cochain (gauge-invariant; see BuildCutCohomologyGenerator).
  result = BuildCutCohomologyGenerator(flux_data, nd_fespace.Get(), mesh, cache, key);
  // Populate the post_op A buffer (used only as scratch downstream).
  auto &A_3d = post_op.GetAGridFunction().Real();
  A_3d.SetFromTrueDofs(result);
}

double ComputeFluxThroughSurface(const mfem::ParGridFunction &B_gf,
                                 const std::vector<int> &attributes, const Mesh &mesh,
                                 const MaterialOperator &mat_op,
                                 const mfem::Vector &flux_direction, MPI_Comm comm)
{
  if (attributes.empty())
  {
    return 0.0;
  }

  // Ensure face neighbor data is available for shared boundary evaluation.
  const_cast<mfem::ParGridFunction &>(B_gf).ExchangeFaceNbrData();

  // Create magnetic flux coefficient with direction-based orientation
  BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC> flux_coeff(
      nullptr, &B_gf, mat_op, false, flux_direction,
      BdrSurfaceFluxCoefficient<SurfaceFlux::MAGNETIC>::OrientationMode::DIRECTION_BASED);

  mfem::ParMesh *pmesh = const_cast<mfem::ParMesh *>(&mesh.Get());
  mfem::ParFiniteElementSpace *fes =
      const_cast<mfem::ParFiniteElementSpace *>(B_gf.ParFESpace());

  // Integrate flux_coeff over the selected boundaries
  double local_flux = 0.0;
  int nbdr = pmesh->GetNBE();
  // Precompute integration rules for each geometry/order combination
  std::map<std::pair<mfem::Geometry::Type, int>, const mfem::IntegrationRule *> ir_map;
  for (int be = 0; be < nbdr; ++be)
  {
    int attr = pmesh->GetBdrAttribute(be);
    if (std::ranges::find(attributes, attr) == attributes.end())
    {
      continue;
    }
    const mfem::FiniteElement *fe = fes->GetBE(be);
    mfem::ElementTransformation *Tr = pmesh->GetBdrElementTransformation(be);
    const int int_order = fem::DefaultIntegrationOrder::Get(*Tr);
    std::pair geom_order{fe->GetGeomType(), int_order};
    const mfem::IntegrationRule *ir;
    auto it = ir_map.find(geom_order);
    if (it == ir_map.end())
    {
      ir = &mfem::IntRules.Get(fe->GetGeomType(), int_order);
      ir_map[geom_order] = ir;
    }
    else
    {
      ir = it->second;
    }
    for (int q = 0; q < ir->GetNPoints(); ++q)
    {
      const mfem::IntegrationPoint &ip = ir->IntPoint(q);
      Tr->SetIntPoint(&ip);
      double val = flux_coeff.Eval(*Tr, ip);
      local_flux += val * ip.weight * Tr->Weight();
    }
  }
  double computed_flux = local_flux;
  Mpi::GlobalSum(1, &computed_flux, comm);
  return computed_flux;
}

}  // namespace palace
