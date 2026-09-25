// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacecurlsolver.hpp"
#include "surfacefluxoperator.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

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
// The angle θ is measured in the plane ⊥ n via a right-handed in-plane basis (e1, e2, n),
// so the construction handles any planar hole orientation and reduces to the xy form when
// n = ±ẑ. a_h is the cut cochain a_h = Grad ψ − a_angle: a lowest-order (Whitney) edge DOF
// of ±Φ across the cut half-plane bounded by L, projected into the order-p ND space
// (ND_1 ⊂ ND_p). It is formed as the discrete gradient of a conformed nodal branch
// potential ψ minus the smooth angle interpolant a_angle so it is non-conformal-safe (see
// below). The cut shape and overall sign are irrelevant to L: the downstream normalization
// rescales a_h so cᵀa_h = Φ.
Vector BuildCutCohomologyGenerator(const SurfaceFluxData &flux_data,
                                   const mfem::ParFiniteElementSpace &ndp_fespace,
                                   const Mesh &mesh)
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

  // Lowest-order (Whitney) cochain on an ND_1 space over the full 3D mesh. Each rank sets
  // its local edge DOFs from vertex COORDINATES (identical on all ranks for a shared edge)
  // in the canonical GetEdgeVertices orientation (ev0→ev1), so the owner's true-DOF value
  // is correct and GetTrueDofs needs no cross-rank sign reconciliation.
  mfem::ND_FECollection nd1_fec(1, sdim);
  mfem::ParFiniteElementSpace nd1_fespace(&pmesh, &nd1_fec);
  mfem::ParGridFunction ah1(&nd1_fespace);
  ah1.UseDevice(false);
  ah1 = 0.0;

  // Nodal branch-cut potential ψ = (Φ/2π)·atan2(d, -s), whose branch discontinuity is
  // exactly the cut half-plane, so cut = Grad ψ - a_angle. Carrying the O(1) step as a
  // discrete gradient is what makes this non-conformal-safe: Grad of a *conformed* H1
  // function is exactly conforming and exactly curl-free, whereas a step written on edges
  // is neither, and cᵀGrad = 0 identically, so ψ cannot perturb the circulation.
  mfem::H1_FECollection h1_fec(1, sdim);
  mfem::ParFiniteElementSpace h1_fespace(&pmesh, &h1_fec);
  mfem::ParGridFunction psi(&h1_fespace);
  psi.UseDevice(false);
  psi = 0.0;
  {
    mfem::Array<int> vdofs;
    for (int v = 0; v < pmesh.GetNV(); v++)
    {
      const double *x = pmesh.GetVertex(v);
      h1_fespace.GetVertexDofs(v, vdofs);
      double s, t;
      inplane(x, s, t);
      psi(vdofs[0]) = phi * std::atan2(s, -t) / (2.0 * M_PI);
    }
    mfem::Vector t(h1_fespace.GetTrueVSize());
    t.UseDevice(false);
    psi.GetTrueDofs(t);
    psi.SetFromTrueDofs(t);
  }

  mfem::Array<int> ev, edofs, vd0, vd1;
  for (int e = 0; e < pmesh.GetNEdges(); e++)
  {
    pmesh.GetEdgeVertices(e, ev);
    const double *x0 = pmesh.GetVertex(ev[0]);
    const double *x1 = pmesh.GetVertex(ev[1]);
    // DOF ∫_e ∇θ·t = angle subtended at L by ev0→ev1. atan2 returns the rotation of
    // magnitude ≤ π, i.e. the branch consistent with the straight edge. Only the in-plane
    // components enter; edges parallel to n subtend zero angle.
    double s0, t0, s1, t1;
    inplane(x0, s0, t0);
    inplane(x1, s1, t1);
    const double a_angle =
        phi * std::atan2(s0 * t1 - t0 * s1, s0 * s1 + t0 * t1) / (2.0 * M_PI);
    nd1_fespace.GetEdgeDofs(e, edofs);
    // cut = Grad ψ - a_angle. ψ is continuous except across the cut half-plane, where it
    // jumps by Φ, so this reproduces the ±Φ step cochain exactly on a conformal mesh.
    h1_fespace.GetVertexDofs(ev[0], vd0);
    h1_fespace.GetVertexDofs(ev[1], vd1);
    ah1(edofs[0]) = psi(vd1[0]) - psi(vd0[0]) - a_angle;
  }

  // Conform the cochain before projecting. Grad ψ is a discrete gradient of a conformed H1
  // function, so it survives the round trip exactly (circulation and curl-free-on-Σ
  // intact); only the smooth a_angle term is re-interpolated onto slaves. On a
  // non-conformal mesh that leaves a small residual curl on Σ near graded refinement, an
  // O(h) discretization error that vanishes under refinement while the downstream cᵀa_h = Φ
  // normalization keeps the fluxoid exact; on a conformal mesh a_h is the exact integer
  // step, curl-free.
  {
    mfem::Vector t(nd1_fespace.GetTrueVSize());
    t.UseDevice(false);
    ah1.GetTrueDofs(t);
    ah1.SetFromTrueDofs(t);
  }

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
    const int bmax = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> marker(bmax);
    marker = 0;
    for (int a : flux_data.film_attributes)
    {
      if (a >= 1 && a <= bmax)
      {
        marker[a - 1] = 1;
      }
    }
    mfem::Array<int> film_tdofs;
    const_cast<mfem::ParFiniteElementSpace &>(ndp_fespace)
        .GetEssentialTrueDofs(marker, film_tdofs);
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

Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op)
{
  Vector result;
  SolveSurfaceCurlProblem(flux_data, mesh, nd_fespace, post_op, result);
  return result;
}

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                             const FiniteElementSpace &nd_fespace,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result)
{
  // London drive: build a_h directly on the 3D ND space as a curl-free cohomology generator
  // via a topological cut cochain (gauge-invariant; see BuildCutCohomologyGenerator).
  result = BuildCutCohomologyGenerator(flux_data, nd_fespace.Get(), mesh);
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
    if (std::find(attributes.begin(), attributes.end(), attr) == attributes.end())
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

void VerifyFluxThroughHoles(const mfem::ParGridFunction &B_gf,
                            const std::vector<int> &hole_attributes,
                            const std::vector<double> &target_fluxes, const Mesh &mesh,
                            const MaterialOperator &mat_op,
                            const mfem::Vector &flux_direction, MPI_Comm comm)
{
  for (std::size_t h = 0; h < hole_attributes.size(); h++)
  {
    int hole_attr = hole_attributes[h];
    double target_flux = target_fluxes[h];

    double computed_flux =
        ComputeFluxThroughSurface(B_gf, {hole_attr}, mesh, mat_op, flux_direction, comm);

    if (Mpi::Root(comm))
    {
      Mpi::Print("Hole attribute {:d}: Target flux = {:.6e}, Computed flux = {:.6e}, Error "
                 "= {:.6e}\n",
                 hole_attr, target_flux, computed_flux,
                 std::abs(computed_flux - target_flux));
    }
  }
}

void VerifyFluxThroughAllHoles(const mfem::ParGridFunction &B_gf, const IoData &iodata,
                               int current_flux_loop_idx, const Mesh &mesh,
                               const MaterialOperator &mat_op, MPI_Comm comm)
{
  if (Mpi::Root(comm))
  {
    Mpi::Print("FluxLoop {:d} excitation - Flux through all holes:\n",
               current_flux_loop_idx);
  }

  // Compute flux through all holes in all flux loops
  for (const auto &[loop_idx, flux_data] : iodata.boundaries.fluxloop)
  {
    for (std::size_t h = 0; h < flux_data.hole_attributes.size(); h++)
    {
      int hole_attr = flux_data.hole_attributes[h];
      double target_flux =
          (loop_idx == current_flux_loop_idx) ? flux_data.flux_amounts[h] : 0.0;

      mfem::Vector flux_direction(const_cast<double *>(flux_data.direction.data()), 3);
      double computed_flux =
          ComputeFluxThroughSurface(B_gf, {hole_attr}, mesh, mat_op, flux_direction, comm);

      Mpi::Print(
          "  Loop {:d} Hole {:d}: Target = {:.6e}, Computed = {:.6e}, Error = {:.6e}\n",
          loop_idx, hole_attr, target_flux, computed_flux,
          std::abs(computed_flux - target_flux));
    }
  }
}

}  // namespace palace
