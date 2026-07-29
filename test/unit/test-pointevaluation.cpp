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

}  // namespace palace
