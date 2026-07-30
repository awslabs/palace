// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/facenbrexchange.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/interpolator.hpp"
#include "fem/mesh.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// Build a simple Cartesian probe mesh.
std::unique_ptr<Mesh> MakeProbeMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

std::unique_ptr<Mesh> MakeFaceExchangeMesh2D(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian2D(2, 2, elem_type);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

std::unique_ptr<Mesh> MakeFaceExchangeMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  smesh.EnsureNodes();
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

// Nonconformally refine one half of a Cartesian mesh to exercise ghost
// master/slave elements without introducing boundary-element fixtures.
std::unique_ptr<Mesh> MakeNCFaceExchangeMesh(MPI_Comm comm, mfem::Element::Type elem_type)
{
  mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, elem_type);
  smesh.EnsureNodes();
  smesh.EnsureNCMesh(true);
  mfem::Array<int> refs;
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center(3);
    smesh.GetElementCenter(e, center);
    if (center(0) < 0.5)
    {
      refs.Append(e);
    }
  }
  smesh.GeneralRefinement(refs, 1);
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  return std::make_unique<Mesh>(std::move(pmesh));
}

}  // namespace

#if defined(MFEM_USE_GSLIB)
TEST_CASE("InterpolationOperator Ceed Probes", "[interpolator][Serial][Parallel][GPU]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  CAPTURE(elem_type, order);
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  auto mesh = MakeProbeMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  GridFunction E(nd_fespace, true);
  mfem::VectorFunctionCoefficient fer(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = std::sin(x(1)) + x(2) * x(2);
                                        v(1) = std::cos(x(2)) + x(0);
                                        v(2) = x(0) * x(1) + 1.0;
                                      });
  mfem::VectorFunctionCoefficient fei(3,
                                      [](const mfem::Vector &x, mfem::Vector &v)
                                      {
                                        v(0) = x(1) * x(2) - 0.5;
                                        v(1) = std::sin(x(0)) - x(2);
                                        v(2) = std::cos(x(1)) + x(0) * x(0);
                                      });
  E.Real().ProjectCoefficient(fer);
  E.Imag().ProjectCoefficient(fei);

  // Probe points are in element interiors and both material regions. Use enough points
  // to select the libCEED path, and compare against direct GSLIB interpolation.
  Units units(1.0, 1.0);
  auto CheckProbes = [&](const std::array<std::array<double, 3>, 8> &pts)
  {
    std::map<int, config::ProbeData> probes;
    for (std::size_t i = 0; i < pts.size(); i++)
    {
      config::ProbeData data;
      data.center = pts[i];
      probes.emplace(static_cast<int>(i + 1), data);
    }
    InterpolationOperator interp(probes, units, nd_fespace);

    const int npts = static_cast<int>(pts.size());
    mfem::Vector xyz(npts * 3), vals_r(npts * 3), vals_i(npts * 3);
    for (int i = 0; i < npts; i++)
    {
      for (int d = 0; d < 3; d++)
      {
        xyz(d * npts + i) = pts[i][d];
      }
    }
    fem::InterpolateFunction(xyz, E.Real(), vals_r, mfem::Ordering::byNODES);
    fem::InterpolateFunction(xyz, E.Imag(), vals_i, mfem::Ordering::byNODES);

    auto vals = interp.ProbeField(E);
    REQUIRE(static_cast<int>(vals.size()) == npts * 3);
    for (int i = 0; i < npts; i++)
    {
      for (int d = 0; d < 3; d++)
      {
        // ProbeField returns byVDIM; the reference uses the requested byNODES ordering.
        const auto val = vals[3 * i + d];
        const double ref_r = vals_r(d * npts + i), ref_i = vals_i(d * npts + i);
        CAPTURE(i, d, ref_r, ref_i);
        CHECK(val.real() == Catch::Approx(ref_r).epsilon(1.0e-9).margin(1.0e-12));
        CHECK(val.imag() == Catch::Approx(ref_i).epsilon(1.0e-9).margin(1.0e-12));
      }
    }
  };

  const std::array<std::array<double, 3>, 8> pts_1 = {{{0.21, 0.37, 0.23},
                                                       {0.74, 0.52, 0.81},
                                                       {0.53, 0.48, 0.52},
                                                       {0.13, 0.68, 0.34},
                                                       {0.32, 0.82, 0.67},
                                                       {0.87, 0.16, 0.42},
                                                       {0.61, 0.29, 0.73},
                                                       {0.43, 0.57, 0.89}}};
  const std::array<std::array<double, 3>, 8> pts_2 = {{{0.17, 0.31, 0.27},
                                                       {0.78, 0.56, 0.84},
                                                       {0.58, 0.44, 0.54},
                                                       {0.11, 0.72, 0.38},
                                                       {0.36, 0.86, 0.63},
                                                       {0.83, 0.19, 0.46},
                                                       {0.66, 0.24, 0.77},
                                                       {0.47, 0.61, 0.92}}};
  CheckProbes(pts_1);
  // CheckProbes destroys its InterpolationOperator. Reconstruct one with different
  // coordinates while nd_fespace and its MFEM DofToQuad caches remain alive.
  CheckProbes(pts_2);
}
#endif

TEST_CASE("FaceNbrFieldExchange 2D", "[facenbrexchange][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto mesh = MakeFaceExchangeMesh2D(comm, mfem::Element::TRIANGLE);
  auto &pmesh = mesh->Get();

  constexpr int order = 2;
  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec);
  mfem::ParGridFunction E(&nd_fespace.Get());
  mfem::VectorFunctionCoefficient field(2,
                                        [](const mfem::Vector &x, mfem::Vector &v)
                                        {
                                          v(0) = std::sin(x(1)) + 0.25 * x(0);
                                          v(1) = std::cos(x(0)) + x(0) * x(1);
                                        });
  E.ProjectCoefficient(field);

  const int num_ghost = pmesh.GetNFaceNeighborElements();
  std::vector<FaceNbrFieldExchange::Request> requests;
  for (int fn = 0; fn < num_ghost; fn++)
  {
    auto &req = requests.emplace_back();
    req.face_nbr_elem = fn;
    req.source_mask = 0b01u;
    req.point_key = {mfem::Element::TRIANGLE, order, 3};
    req.pts.resize(3);
    req.pts[0].Set2(0.1, 0.2);
    req.pts[1].Set2(0.2, 0.1);
    req.pts[2].Set2(0.25, 0.25);
  }
  FaceNbrFieldExchange exchange(*mesh, {&nd_fespace.Get(), nullptr, nullptr, nullptr},
                                requests);
  exchange.Exchange({&E, nullptr, nullptr, nullptr});

  E.ExchangeFaceNbrData();
  const double *values = exchange.Imported().HostRead();
  mfem::Vector ref(2);
  int num_checked = 0;
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    const int offset = exchange.ImportOffset(static_cast<int>(r), 0);
    REQUIRE(offset >= 0);
    for (std::size_t j = 0; j < req.pts.size(); j++)
    {
      E.GetVectorValue(pmesh.GetNE() + req.face_nbr_elem, req.pts[j], ref);
      for (int c = 0; c < 2; c++)
      {
        CAPTURE(r, j, c);
        CHECK(values[offset + 2 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
        num_checked++;
      }
    }
  }
  int num_global = num_checked;
  Mpi::GlobalSum(1, &num_global, comm);
  CHECK((Mpi::Size(comm) == 1 || num_global > 0));
}

TEST_CASE("FaceNbrFieldExchange", "[facenbrexchange][Serial][Parallel]")
{
  MPI_Comm comm = MPI_COMM_WORLD;
  auto elem_type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  auto order = GENERATE(1, 2);
  auto nonconformal = GENERATE(false, true);
  CAPTURE(elem_type, order, nonconformal);

  auto mesh = nonconformal ? MakeNCFaceExchangeMesh(comm, elem_type)
                           : MakeFaceExchangeMesh(comm, elem_type);
  auto &pmesh = mesh->Get();

  mfem::ND_FECollection nd_fec(order, pmesh.Dimension());
  mfem::RT_FECollection rt_fec(order - 1, pmesh.Dimension());
  FiniteElementSpace nd_fespace(*mesh, &nd_fec), rt_fespace(*mesh, &rt_fec);

  // Project non-trivial smooth fields. The reference values are computed from the same
  // projected grid functions through the legacy mfem face neighbor dof exchange, so the
  // comparison is exact up to roundoff (not projection error).
  mfem::ParGridFunction E(&nd_fespace.Get()), B(&rt_fespace.Get());
  mfem::VectorFunctionCoefficient fe(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = std::sin(x(1)) + x(2) * x(2);
                                       v(1) = std::cos(x(2)) + x(0);
                                       v(2) = x(0) * x(1) + 1.0;
                                     });
  mfem::VectorFunctionCoefficient fb(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(1) * x(2) - 0.5;
                                       v(1) = std::sin(x(0)) - x(2);
                                       v(2) = std::cos(x(1)) + x(0) * x(0);
                                     });
  E.ProjectCoefficient(fe);
  B.ProjectCoefficient(fb);

  // Request E (slot 0) at a few reference points of every ghost element, and both E
  // and B (slot 1) for every other ghost element (exercising the value layouts). The
  // points are valid reference coordinates for both tetrahedra and hexahedra.
  const int num_ghost = pmesh.GetNFaceNeighborElements();
  std::vector<FaceNbrFieldExchange::Request> requests;
  for (int fn = 0; fn < num_ghost; fn++)
  {
    auto &req = requests.emplace_back();
    req.face_nbr_elem = fn;
    req.source_mask = (fn % 2 == 0) ? 0b01u : 0b11u;
    req.point_key = {static_cast<long long>(elem_type), static_cast<long long>(order), 4};
    req.pts.resize(4);
    req.pts[0].Set3(0.1, 0.2, 0.3);
    req.pts[1].Set3(0.25, 0.25, 0.25);
    req.pts[2].Set3(0.05, 0.1, 0.7);
    req.pts[3].Set3(0.3, 0.05, 0.05);
  }
  FaceNbrFieldExchange exchange(
      *mesh, {&nd_fespace.Get(), &rt_fespace.Get(), nullptr, nullptr}, requests);
  exchange.Exchange({&E, &B, nullptr, nullptr});

  // Reference: evaluate the ghost elements through the legacy dof exchange.
  E.ExchangeFaceNbrData();
  B.ExchangeFaceNbrData();
  std::vector<double> vals(exchange.Imported().HostRead(),
                           exchange.Imported().HostRead() + exchange.ImportSize());
  mfem::Vector ref(3);
  int num_checked = 0;
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    for (int s = 0; s < 2; s++)
    {
      const int offset = exchange.ImportOffset(static_cast<int>(r), s);
      if (!(req.source_mask & (1u << s)))
      {
        CHECK(offset < 0);
        continue;
      }
      const auto &U = (s == 0) ? E : B;
      for (std::size_t j = 0; j < req.pts.size(); j++)
      {
        U.GetVectorValue(pmesh.GetNE() + req.face_nbr_elem, req.pts[j], ref);
        for (int c = 0; c < 3; c++)
        {
          CAPTURE(r, s, j, c);
          CHECK(vals[offset + 3 * j + c] == Catch::Approx(ref(c)).margin(1.0e-12));
          num_checked++;
        }
      }
    }
  }
  // With more than one process, the partition must produce at least one ghost element
  // somewhere (the exchange is the point of the test).
  int num_global = num_checked;
  Mpi::GlobalSum(1, &num_global, comm);
  CHECK((Mpi::Size(comm) == 1 || num_global > 0));

  // The field inputs are re-pointed at the sources on each call: scaling the field
  // scales the exchanged values.
  E *= 2.0;
  exchange.Exchange({&E, &B, nullptr, nullptr});
  const double *vals2 = exchange.Imported().HostRead();
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const int offset = exchange.ImportOffset(static_cast<int>(r), 0);
    for (std::size_t j = 0; j < 3 * requests[r].pts.size(); j++)
    {
      CHECK(vals2[offset + j] == Catch::Approx(2.0 * vals[offset + j]).margin(1.0e-12));
    }
  }
}

}  // namespace palace
