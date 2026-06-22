// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/surfacefunctional.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/labels.hpp"

namespace palace
{

namespace
{

std::unique_ptr<Mesh> MakeSmallTetInterfaceMesh(MPI_Comm comm)
{
  // 2 x 2 x 2 cubes split into tetrahedra: small enough for a unit test, but it has
  // exterior boundaries, an interior material interface, and enough face orientations to
  // exercise the boundary-normal conventions.
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }

  // Add boundary elements on the interior material interface z = 0.5. This makes one
  // SurfaceFunctional instance cover both one-sided and two-sided boundary-viz kernels.
  for (int f = 0; f < smesh.GetNumFaces(); f++)
  {
    int e1, e2;
    smesh.GetFaceElements(f, &e1, &e2);
    if (e1 >= 0 && e2 >= 0 && smesh.GetAttribute(e1) != smesh.GetAttribute(e2))
    {
      auto *face_elem = smesh.GetFace(f)->Duplicate(&smesh);
      face_elem->SetAttribute(7);
      smesh.AddBdrElement(face_elem);
    }
  }
  smesh.FinalizeTopology();
  smesh.Finalize();
  smesh.SetAttributes();
  smesh.EnsureNodes();

  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

std::vector<config::MaterialData> MakeTwoMaterials()
{
  config::MaterialData lower, upper;
  lower.attributes = {1};
  upper.attributes = {2};

  // Use diagonal anisotropy and different values on each side so material lookup and
  // component ordering are tested, not just the identity/isotropic path.
  lower.epsilon_r.s = {1.1, 1.2, 1.3};
  lower.mu_r.s = {0.9, 1.0, 1.1};
  upper.epsilon_r.s = {11.7, 3.1, 2.4};
  upper.mu_r.s = {1.4, 1.8, 2.2};

  return {lower, upper};
}

void ProjectTestFields(GridFunction &E, GridFunction &B)
{
  mfem::VectorFunctionCoefficient er(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(1.3 * x(1)) + x(2) * x(2) + 0.2;
                                       v(1) = std::cos(0.7 * x(2)) + 0.4 * x(0);
                                       v(2) = x(0) * x(1) - 0.3 * x(2) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient ei(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) * x(2) - 0.5;
                                       v(1) = std::sin(x(0)) - 0.25 * x(2);
                                       v(2) = std::cos(1.1 * x(1)) + x(0) * x(0);
                                     });
  mfem::VectorFunctionCoefficient br(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) - 0.3 * x(2) + 0.4;
                                       v(1) = std::sin(0.9 * x(2)) + 0.5;
                                       v(2) = std::cos(x(0)) - x(1) * x(2);
                                     });
  mfem::VectorFunctionCoefficient bi(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::cos(x(2)) - 0.2;
                                       v(1) = x(0) * x(2) + 0.1;
                                       v(2) = std::sin(x(1)) - 0.4 * x(0);
                                     });

  E.Real().ProjectCoefficient(er);
  E.Imag().ProjectCoefficient(ei);
  B.Real().ProjectCoefficient(br);
  B.Imag().ProjectCoefficient(bi);

  // Legacy BdrGridFunctionCoefficient evaluation reads face-neighbor data in parallel.
  E.Real().ExchangeFaceNbrData();
  E.Imag().ExchangeFaceNbrData();
  B.Real().ExchangeFaceNbrData();
  B.Imag().ExchangeFaceNbrData();
}

struct ErrorStats
{
  long long count = 0;
  double max_abs = 0.0;
  double max_scaled = 0.0;
  double sum_sq = 0.0;

  void Add(double val, double ref, double rtol, double atol)
  {
    const double err = std::abs(val - ref);
    const double denom = atol + rtol * std::max({1.0, std::abs(val), std::abs(ref)});
    max_abs = std::max(max_abs, err);
    max_scaled = std::max(max_scaled, err / denom);
    sum_sq += err * err;
    count++;
  }

  void Reduce(MPI_Comm comm)
  {
    Mpi::GlobalSum(1, &count, comm);
    Mpi::GlobalMax(1, &max_abs, comm);
    Mpi::GlobalMax(1, &max_scaled, comm);
    Mpi::GlobalSum(1, &sum_sq, comm);
  }

  double Rms() const { return count ? std::sqrt(sum_sq / static_cast<double>(count)) : 0.0; }
};

void RequireGloballyConsistentValidity(const SurfaceFunctional &viz, MPI_Comm comm)
{
  bool valid = viz.IsValid();
  bool valid_and = valid, valid_or = valid;
  Mpi::GlobalAnd(1, &valid_and, comm);
  Mpi::GlobalOr(1, &valid_or, comm);
  REQUIRE(valid_and == valid_or);
}

void CheckStats(const std::string &name, ErrorStats stats)
{
  INFO("field = " << name);
  INFO("count = " << stats.count);
  INFO("max_abs = " << stats.max_abs);
  INFO("max_scaled = " << stats.max_scaled);
  INFO("rms = " << stats.Rms());
  CHECK(stats.count > 0);
  CHECK(stats.max_scaled <= 1.0);
}

ErrorStats CompareScalarBuffer(mfem::ParMesh &pmesh, const mfem::Array<int> &marker,
                               int lod, const SurfaceFunctional &viz,
                               const Vector &buffer, mfem::Coefficient &legacy,
                               double rtol, double atol)
{
  ErrorStats stats;
  const double *buf = buffer.HostRead();
  const auto &bases = viz.BufferBases();
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    MFEM_VERIFY(i < static_cast<int>(bases.size()) && bases[i] >= 0,
                "Missing boundary visualization buffer slot for marked element!");
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
    auto *T = pmesh.GetBdrElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      stats.Add(buf[bases[i] + j], legacy.Eval(*T, ip), rtol, atol);
    }
  }
  stats.Reduce(pmesh.GetComm());
  return stats;
}

ErrorStats CompareVectorBuffer(mfem::ParMesh &pmesh, const mfem::Array<int> &marker,
                               int lod, const SurfaceFunctional &viz,
                               const Vector &buffer, mfem::VectorCoefficient &legacy,
                               double rtol, double atol)
{
  ErrorStats stats;
  const double *buf = buffer.HostRead();
  const auto &bases = viz.BufferBases();
  mfem::Vector ref(3);
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    if (!marker[pmesh.GetBdrAttribute(i) - 1])
    {
      continue;
    }
    MFEM_VERIFY(i < static_cast<int>(bases.size()) && bases[i] >= 0,
                "Missing boundary visualization buffer slot for marked element!");
    const auto &RefG =
        *mfem::GlobGeometryRefiner.Refine(pmesh.GetBdrElementGeometry(i), lod, 1);
    auto *T = pmesh.GetBdrElementTransformation(i);
    for (int j = 0; j < RefG.RefPts.GetNPoints(); j++)
    {
      const auto &ip = RefG.RefPts.IntPoint(j);
      T->SetIntPoint(&ip);
      legacy.Eval(ref, *T, ip);
      for (int c = 0; c < 3; c++)
      {
        stats.Add(buf[bases[i] + 3 * j + c], ref(c), rtol, atol);
      }
    }
  }
  stats.Reduce(pmesh.GetComm());
  return stats;
}

template <typename Kind, typename = void>
struct HasBoundaryPoyntingKind : std::false_type
{
};

template <typename Kind>
struct HasBoundaryPoyntingKind<Kind, std::void_t<decltype(Kind::BDR_POYNTING)>>
  : std::true_type
{
};

template <typename Kind = SurfaceFunctional::Kind,
          typename std::enable_if<HasBoundaryPoyntingKind<Kind>::value, int>::type = 0>
void CheckBoundaryPoyntingIfAvailable(Mesh &mesh, mfem::ParMesh &pmesh,
                                      const mfem::Array<int> &marker,
                                      const mfem::ParFiniteElementSpace &nd_fespace,
                                      const mfem::ParFiniteElementSpace &rt_fespace,
                                      const MaterialOperator &mat_op, const GridFunction &E,
                                      const GridFunction &B, int lod, double scaling,
                                      double rtol, double atol)
{
  SurfaceFunctional viz(Kind::BDR_POYNTING, mesh, marker, nd_fespace, rt_fespace, mat_op,
                        lod, scaling);
  RequireGloballyConsistentValidity(viz, pmesh.GetComm());
  if (!viz.IsValid())
  {
    return;
  }
  Vector buffer(viz.BufferSize());
  buffer.UseDevice(true);
  viz.EvalBuffer(E, B, buffer);
  PoyntingVectorCoefficient legacy(E, B, mat_op, scaling);
  CheckStats("S", CompareVectorBuffer(pmesh, marker, lod, viz, buffer, legacy, rtol, atol));
}

template <typename Kind = SurfaceFunctional::Kind,
          typename std::enable_if<!HasBoundaryPoyntingKind<Kind>::value, int>::type = 0>
void CheckBoundaryPoyntingIfAvailable(Mesh &, mfem::ParMesh &, const mfem::Array<int> &,
                                      const mfem::ParFiniteElementSpace &,
                                      const mfem::ParFiniteElementSpace &,
                                      const MaterialOperator &, const GridFunction &,
                                      const GridFunction &, int, double, double, double)
{
  // Older branches do not have a boundary libCEED Poynting buffer yet. The same test
  // file still validates the boundary kernels that exist there; branches with
  // BDR_POYNTING automatically instantiate the overload above.
}

}  // namespace

TEST_CASE("SurfaceFunctional boundary visualization buffers match legacy coefficients",
          "[surfacefunctional][boundary-viz][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  constexpr int order = 2;
  constexpr double rtol = 1.0e-8;
  constexpr double atol = 1.0e-10;

  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeSmallTetInterfaceMesh(comm);
  auto &pmesh = mesh->Get();
  auto materials = MakeTwoMaterials();
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op(materials, periodic, ProblemType::DRIVEN, *mesh);

  mfem::ND_FECollection nd_fec(order, 3);
  mfem::RT_FECollection rt_fec(order - 1, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  GridFunction E(nd_fespace, true), B(rt_fespace, true);
  ProjectTestFields(E, B);

  const int lod = order;
  const int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
  mfem::Array<int> marker(bdr_attr_max);
  marker = 1;

  auto CheckField = [&](const std::string &name, SurfaceFunctional::Kind kind,
                        const mfem::ParFiniteElementSpace &fes,
                        const mfem::ParGridFunction &u,
                        mfem::VectorCoefficient &legacy)
  {
    SurfaceFunctional viz(kind, *mesh, marker, fes, lod);
    RequireGloballyConsistentValidity(viz, comm);
    if (!viz.IsValid())
    {
      return;
    }
    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    viz.EvalBuffer(u, buffer);
    CheckStats(name, CompareVectorBuffer(pmesh, marker, lod, viz, buffer, legacy, rtol,
                                         atol));
  };

  auto CheckCoeff = [&](const std::string &name, SurfaceFunctional::Kind kind,
                        const mfem::ParFiniteElementSpace &fes, const Vector &u,
                        mfem::Coefficient *legacy_s, mfem::VectorCoefficient *legacy_v,
                        double scaling)
  {
    SurfaceFunctional viz(kind, *mesh, marker, fes, mat_op, lod, scaling);
    RequireGloballyConsistentValidity(viz, comm);
    if (!viz.IsValid())
    {
      return;
    }
    Vector buffer(viz.BufferSize());
    buffer.UseDevice(true);
    viz.EvalBuffer(u, buffer);
    if (legacy_s)
    {
      CheckStats(name, CompareScalarBuffer(pmesh, marker, lod, viz, buffer, *legacy_s,
                                           rtol, atol));
    }
    else
    {
      REQUIRE(legacy_v);
      CheckStats(name, CompareVectorBuffer(pmesh, marker, lod, viz, buffer, *legacy_v,
                                           rtol, atol));
    }
  };

  // Direct boundary field buffers: legacy BdrFieldVectorCoefficient averages interior
  // boundary values and evaluates exterior values from the attached volume element.
  BdrFieldVectorCoefficient E_real_legacy(E.Real()), E_imag_legacy(E.Imag());
  BdrFieldVectorCoefficient B_real_legacy(B.Real()), B_imag_legacy(B.Imag());
  CheckField("E_real", SurfaceFunctional::Kind::BDR_FIELD_E, nd_fespace.Get(), E.Real(),
             E_real_legacy);
  CheckField("E_imag", SurfaceFunctional::Kind::BDR_FIELD_E, nd_fespace.Get(), E.Imag(),
             E_imag_legacy);
  CheckField("B_real", SurfaceFunctional::Kind::BDR_FIELD_B, rt_fespace.Get(), B.Real(),
             B_real_legacy);
  CheckField("B_imag", SurfaceFunctional::Kind::BDR_FIELD_B, rt_fespace.Get(), B.Imag(),
             B_imag_legacy);

  // Material/normal-dependent buffers. Use non-unit scales to exercise the QFunction
  // context slots used by production dimensionalization.
  constexpr double eps_scaling = 1.3;
  constexpr double invmu_scaling = 0.7;
  mfem::Vector unused_x0;
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_real_legacy(
      &E.Real(), nullptr, mat_op, true, unused_x0, eps_scaling);
  BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC> Q_imag_legacy(
      &E.Imag(), nullptr, mat_op, true, unused_x0, eps_scaling);
  BdrSurfaceCurrentVectorCoefficient J_real_legacy(B.Real(), mat_op, invmu_scaling);
  BdrSurfaceCurrentVectorCoefficient J_imag_legacy(B.Imag(), mat_op, invmu_scaling);
  CheckCoeff("Q_s_real", SurfaceFunctional::Kind::BDR_FLUX_Q, nd_fespace.Get(), E.Real(),
             &Q_real_legacy, nullptr, eps_scaling);
  CheckCoeff("Q_s_imag", SurfaceFunctional::Kind::BDR_FLUX_Q, nd_fespace.Get(), E.Imag(),
             &Q_imag_legacy, nullptr, eps_scaling);
  CheckCoeff("J_s_real", SurfaceFunctional::Kind::BDR_CURRENT_J, rt_fespace.Get(),
             B.Real(), nullptr, &J_real_legacy, invmu_scaling);
  CheckCoeff("J_s_imag", SurfaceFunctional::Kind::BDR_CURRENT_J, rt_fespace.Get(),
             B.Imag(), nullptr, &J_imag_legacy, invmu_scaling);

  // Quadratic buffers use the GridFunction overload so real and imaginary contributions
  // accumulate into a single output buffer, matching production output semantics.
  {
    SurfaceFunctional viz(SurfaceFunctional::Kind::BDR_ENERGY_E, *mesh, marker,
                          nd_fespace.Get(), mat_op, lod, eps_scaling);
    RequireGloballyConsistentValidity(viz, comm);
    if (viz.IsValid())
    {
      Vector buffer(viz.BufferSize());
      buffer.UseDevice(true);
      viz.EvalBuffer(E, buffer);
      EnergyDensityCoefficient<EnergyDensityType::ELECTRIC> legacy(E, mat_op, eps_scaling);
      CheckStats("U_e", CompareScalarBuffer(pmesh, marker, lod, viz, buffer, legacy, rtol,
                                            atol));
    }
  }
  {
    SurfaceFunctional viz(SurfaceFunctional::Kind::BDR_ENERGY_M, *mesh, marker,
                          rt_fespace.Get(), mat_op, lod, invmu_scaling);
    RequireGloballyConsistentValidity(viz, comm);
    if (viz.IsValid())
    {
      Vector buffer(viz.BufferSize());
      buffer.UseDevice(true);
      viz.EvalBuffer(B, buffer);
      EnergyDensityCoefficient<EnergyDensityType::MAGNETIC> legacy(B, mat_op,
                                                                    invmu_scaling);
      CheckStats("U_m", CompareScalarBuffer(pmesh, marker, lod, viz, buffer, legacy, rtol,
                                            atol));
    }
  }

  CheckBoundaryPoyntingIfAvailable(*mesh, pmesh, marker, nd_fespace.Get(),
                                   rt_fespace.Get(), mat_op, E, B, lod, invmu_scaling,
                                   rtol, atol);
}

}  // namespace palace
