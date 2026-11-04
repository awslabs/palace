// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <numbers>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/planewaveoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// Unit cube of n x n x n hexahedra split into tetrahedra (n a multiple of 3), with the
// volume attributes 1 for the inner cube [1/3, 2/3]^3, 2 for the slab x < 1/3, and 3 for
// the rest.
std::unique_ptr<mfem::ParMesh> MakeMesh(int n)
{
  auto mesh =
      mfem::Mesh::MakeCartesian3D(n, n, n, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
  mfem::Vector c;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    mesh.GetElementCenter(e, c);
    bool inner = true;
    for (int d = 0; d < 3; d++)
    {
      inner = inner && c(d) > 1.0 / 3.0 && c(d) < 2.0 / 3.0;
    }
    mesh.SetAttribute(e, inner ? 1 : (c(0) < 1.0 / 3.0 ? 2 : 3));
  }
  mesh.SetAttributes();
  return std::make_unique<mfem::ParMesh>(Mpi::World(), mesh);
}

config::PlaneWaveData MakePlaneWave(const std::vector<int> &attributes,
                                    const std::array<double, 3> &direction,
                                    const std::array<double, 3> &polarization,
                                    double amplitude = 1.0,
                                    const std::array<double, 3> &origin = {0.0, 0.0, 0.0})
{
  config::PlaneWaveData data;
  data.attributes = attributes;
  data.direction = direction;
  data.polarization = polarization;
  data.amplitude = amplitude;
  data.origin = origin;
  return data;
}

}  // namespace

TEST_CASE("PlaneWaveOperator source", "[planewaveoperator][Serial]")
{
  auto mesh = MakeMesh(3);
  Units units(0.5, 2.0);

  // The direction is normalized and the polarization is made orthogonal to it.
  std::map<int, config::PlaneWaveData> plane_wave = {
      {1, MakePlaneWave({1, 3}, {0.0, 0.0, 2.0}, {1.0, 0.0, 1.0}, 3.0)}};
  PlaneWaveOperator op(plane_wave, units, *mesh);
  REQUIRE(op.Size() == 1);
  const auto &src = op.GetSource(1);
  const std::array<double, 3> k_ref = {0.0, 0.0, 1.0}, p_ref = {1.0, 0.0, 0.0};
  for (int d = 0; d < 3; d++)
  {
    CHECK(src.direction(d) == Catch::Approx(k_ref[d]).margin(1.0e-14));
    CHECK(src.polarization(d) == Catch::Approx(p_ref[d]).margin(1.0e-14));
  }
  CHECK(src.amplitude ==
        Catch::Approx(units.Nondimensionalize<Units::ValueType::FIELD_E>(3.0)));
  REQUIRE(src.attr_marker.Size() == 3);
  CHECK(src.attr_marker[0] == 1);
  CHECK(src.attr_marker[1] == 0);
  CHECK(src.attr_marker[2] == 1);

  // Invalid sources.
  auto Make = [&](const config::PlaneWaveData &data)
  { return PlaneWaveOperator({{1, data}}, units, *mesh); };
  CHECK_THROWS(Make(MakePlaneWave({1}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0})));
  CHECK_THROWS(Make(MakePlaneWave({1}, {0.0, 0.0, 1.0}, {0.0, 0.0, -2.0})));
  CHECK_THROWS(Make(MakePlaneWave({4}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})));
  CHECK_THROWS(Make(MakePlaneWave({0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})));
}

TEST_CASE("PlaneWaveOperator total-field marker", "[planewaveoperator][Serial][Parallel]")
{
  // The marker selects the true dofs of the closure of the total-field region: with
  // Nédélec spaces of order 2 and above on tetrahedra, face dofs on the triangles shared
  // between processes need the marker to be combined as for essential dofs. The region is
  // either the inner cube or most of the domain, so that in parallel the interfaces
  // between processes cross it.
  const int order = GENERATE(1, 2, 3);
  const auto tf_attributes = GENERATE(std::vector<int>{1}, std::vector<int>{1, 3});
  Mesh mesh(MakeMesh(6));
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  Units units(1.0, 1.0);
  std::map<int, config::PlaneWaveData> plane_wave = {
      {1, MakePlaneWave(tf_attributes, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})}};
  PlaneWaveOperator op(plane_wave, units, mesh.Get());

  const Vector &marker = op.GetTotalFieldMarker(1, nd_fespace);
  REQUIRE(marker.Size() == nd_fespace.GetTrueVSize());
  int counts[2] = {0, 0};  // Marked dofs, entries which are neither 0 nor 1
  for (int i = 0; i < marker.Size(); i++)
  {
    counts[0] += (marker(i) == 1.0);
    counts[1] += (marker(i) != 0.0 && marker(i) != 1.0);
  }
  Mpi::GlobalSum(2, counts, nd_fespace.GetComm());
  CHECK(counts[1] == 0);

  // The true dofs of the closure of the region are the true dofs of the same space on the
  // region alone, independently of the number of processes.
  mfem::Array<int> tf_attr_array;
  for (auto attr : tf_attributes)
  {
    tf_attr_array.Append(attr);
  }
  auto submesh = mfem::ParSubMesh::CreateFromDomain(mesh.Get(), tf_attr_array);
  mfem::ND_FECollection sub_fec(order, 3);
  mfem::ParFiniteElementSpace sub_fespace(&submesh, &sub_fec);
  CHECK(counts[0] == sub_fespace.GlobalTrueVSize());
}

TEST_CASE("PlaneWaveOperator right-hand side", "[planewaveoperator][Serial][Parallel]")
{
  // With x = M g, the incident field inside of S (the marked dofs) and no scattered field
  // outside, the residual of the TF/SF system is A x - b = M A g: it vanishes on all the
  // scattered-field dofs and it is the residual of the discrete incident field on the
  // total-field dofs. The interpolant g of the incident field is computed here
  // independently, so this also checks the phase, amplitude, and polarization conventions.
  const int order = GENERATE(1, 2, 3);
  Mesh mesh(MakeMesh(6));
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  Units units(1.0, 1.0);
  const double omega = 2.0 * std::numbers::pi, E0 = 1.5;
  const std::array<double, 3> x0 = {0.2, 0.1, -0.3};
  std::map<int, config::PlaneWaveData> plane_wave = {
      {1, MakePlaneWave({1}, {1.0, 2.0, 2.0}, {0.0, 1.0, 0.0}, E0, x0)}};
  PlaneWaveOperator op(plane_wave, units, mesh.Get());

  // System operator of the total field in vacuum, without boundary conditions.
  mfem::ParBilinearForm a(&nd_fespace.Get());
  mfem::ConstantCoefficient neg_omega2(-omega * omega);
  a.AddDomainIntegrator(new mfem::CurlCurlIntegrator());
  a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(neg_omega2));
  a.Assemble();
  a.Finalize();
  std::unique_ptr<mfem::HypreParMatrix> K(a.ParallelAssemble());
  ComplexWrapperOperator A(K.get(), nullptr);

  const int n = nd_fespace.GetTrueVSize();
  ComplexVector b(n);
  b = 0.0;
  REQUIRE(op.AddExcitationVector(1, omega, A, nd_fespace, b));
  CHECK_FALSE(op.AddExcitationVector(2, omega, A, nd_fespace, b));

  // Incident field E0 p exp(-i omega khat . (x - x0)), with khat = (1, 2, 2) / 3 and p the
  // normalized part of +Y orthogonal to khat, in nondimensional units.
  const double E0_nd = units.Nondimensionalize<Units::ValueType::FIELD_E>(E0);
  const std::array<double, 3> khat = {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
  const double p_norm = std::sqrt(45.0);
  const std::array<double, 3> p = {-2.0 / p_norm, 5.0 / p_norm, -4.0 / p_norm};
  ComplexVector g(n);
  for (int part = 0; part < 2; part++)
  {
    mfem::VectorFunctionCoefficient coef(
        3,
        [&](const mfem::Vector &x, mfem::Vector &E)
        {
          double phase = 0.0;
          for (int d = 0; d < 3; d++)
          {
            phase += khat[d] * (x(d) - x0[d]);
          }
          const auto f = E0_nd * std::exp(std::complex<double>(0.0, -omega * phase));
          for (int d = 0; d < 3; d++)
          {
            E(d) = ((part == 0) ? f.real() : f.imag()) * p[d];
          }
        });
    mfem::ParGridFunction gf(&nd_fespace.Get());
    gf.ProjectCoefficient(coef);
    gf.ParallelProject((part == 0) ? g.Real() : g.Imag());
  }

  const Vector &marker = op.GetTotalFieldMarker(1, nd_fespace);
  ComplexVector x(n), Ax(n), Ag(n);
  for (int i = 0; i < n; i++)
  {
    x.Real()(i) = marker(i) * g.Real()(i);
    x.Imag()(i) = marker(i) * g.Imag()(i);
  }
  A.Mult(x, Ax);
  A.Mult(g, Ag);
  double err[3] = {0.0, 0.0, 0.0};  // Scattered-field rows, total-field rows, scale
  for (int i = 0; i < n; i++)
  {
    const std::complex<double> r(Ax.Real()(i) - b.Real()(i), Ax.Imag()(i) - b.Imag()(i));
    const std::complex<double> Agi(Ag.Real()(i), Ag.Imag()(i));
    if (marker(i) == 0.0)
    {
      err[0] = std::max(err[0], std::abs(r));
    }
    else
    {
      err[1] = std::max(err[1], std::abs(r - Agi));
    }
    err[2] = std::max(err[2], std::abs(std::complex<double>(b.Real()(i), b.Imag()(i))));
  }
  Mpi::GlobalMax(3, err, nd_fespace.GetComm());
  REQUIRE(err[2] > 0.0);
  CHECK(err[0] <= 1.0e-12 * err[2]);
  CHECK(err[1] <= 1.0e-12 * err[2]);
}

TEST_CASE("PlaneWaveOperator boundary touching the total field",
          "[planewaveoperator][Serial][Parallel]")
{
  const int order = GENERATE(1, 2);
  Mesh mesh(MakeMesh(3));
  mfem::ND_FECollection nd_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  Units units(1.0, 1.0);
  std::map<int, config::PlaneWaveData> plane_wave = {
      {1, MakePlaneWave({1}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})},
      {2, MakePlaneWave({2}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})},
      {3, MakePlaneWave({1, 2, 3}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0})}};
  PlaneWaveOperator op(plane_wave, units, mesh.Get());

  // Boundary attributes of the faces x = 0 and x = 1 of the cube.
  const auto &pmesh = mesh.Get();
  int attr_x[2] = {0, 0};
  mfem::Array<int> v;
  for (int i = 0; i < pmesh.GetNBE(); i++)
  {
    pmesh.GetBdrElementVertices(i, v);
    double xc = 0.0;
    for (int j : v)
    {
      xc += pmesh.GetVertex(j)[0] / v.Size();
    }
    if (std::abs(xc) < 1.0e-12)
    {
      attr_x[0] = pmesh.GetBdrAttribute(i);
    }
    else if (std::abs(xc - 1.0) < 1.0e-12)
    {
      attr_x[1] = pmesh.GetBdrAttribute(i);
    }
  }
  Mpi::GlobalMax(2, attr_x, nd_fespace.GetComm());
  REQUIRE(attr_x[0] > 0);
  REQUIRE(attr_x[1] > 0);
  auto Marker = [&](const std::vector<int> &attributes)
  {
    mfem::Array<int> marker(pmesh.bdr_attributes.Max());
    marker = 0;
    for (auto attr : attributes)
    {
      marker[attr - 1] = 1;
    }
    return marker;
  };
  std::vector<int> all_attributes;
  for (int attr : pmesh.bdr_attributes)
  {
    all_attributes.push_back(attr);
  }

  // The inner cube is separated from the boundary by a layer of elements.
  CHECK_FALSE(op.TouchesTotalField(1, nd_fespace, Marker(all_attributes)));

  // The slab x < 1/3 touches the face x = 0, not the face x = 1.
  CHECK(op.TouchesTotalField(2, nd_fespace, Marker({attr_x[0]})));
  CHECK_FALSE(op.TouchesTotalField(2, nd_fespace, Marker({attr_x[1]})));

  // The whole cube touches every face.
  CHECK(op.TouchesTotalField(3, nd_fespace, Marker({attr_x[1]})));
}

}  // namespace palace
