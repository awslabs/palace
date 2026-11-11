// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/phasereference.hpp"
#include "fixtures.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using json = nlohmann::json;
namespace fs = std::filesystem;
using namespace std::complex_literals;

namespace
{

// A tetrahedral mesh of a box with its elements in a random order, written to a file in the
// given directory.
fs::path WriteShuffledBox(const fs::path &dir)
{
  auto mesh =
      mfem::Mesh::MakeCartesian3D(4, 5, 3, mfem::Element::TETRAHEDRON, 1.0, 1.3, 0.7);
  mfem::Array<int> ordering(mesh.GetNE());
  std::iota(ordering.begin(), ordering.end(), 0);
  std::mt19937 gen(4321);
  std::shuffle(ordering.begin(), ordering.end(), gen);
  mesh.ReorderElements(ordering);
  const auto path = dir / "palace-test-phasereference.mesh";
  if (Mpi::Root(Mpi::World()))
  {
    std::ofstream fo(path);
    fo.precision(17);
    mesh.Print(fo);
  }
  Mpi::Barrier(Mpi::World());
  return path;
}

IoData MakeIoData(const fs::path &mesh_path, bool local_reorder)
{
  json config = {{"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
                 {"Model", {{"Mesh", mesh_path.string()}, {"LocalReorder", local_reorder}}},
                 {"Domains", {{"Materials", {{{"Attributes", {1}}}}}}},
                 {"Boundaries", {{"PEC", {{"Attributes", {1, 2, 3, 4, 5, 6}}}}}},
                 {"Solver", {{"Eigenmode", {{"Target", 1.0}}}}}};
  return IoData(config, false);
}

// Real and imaginary parts of a smooth field of no particular symmetry, and test fields
// for functionals of the normalized field.
void FieldRe(const mfem::Vector &x, mfem::Vector &v)
{
  v(0) = std::sin(2.0 * x(1)) + 0.3 * x(2);
  v(1) = std::cos(1.5 * x(0) * x(2)) - 0.2;
  v(2) = x(0) * x(1) - 0.5 * x(2) * x(2);
}

void FieldIm(const mfem::Vector &x, mfem::Vector &v)
{
  v(0) = 0.4 * x(0) - x(1) * x(2);
  v(1) = std::sin(x(0) + 2.0 * x(2));
  v(2) = 0.1 + std::cos(3.0 * x(1));
}

void TestField(int j, const mfem::Vector &x, mfem::Vector &v)
{
  v(0) = std::cos(1.0 + j + x(0) * x(1));
  v(1) = x(2) - 0.5 * j;
  v(2) = std::sin(2.0 * x(0) - j * x(2));
}

struct Normalized
{
  int ref;                              // reference used
  std::vector<std::complex<double>> g;  // ∫ gⱼ ⋅ E dΩ of the normalized field
  double imag;                          // ‖Im E‖ / ‖E‖ of the normalized field
};

// Interpolate the field (or its real part) multiplied by exp(iθ) on a second-order Nédélec
// space, normalize its phase and measure it.
Normalized NormalizeField(Mesh &mesh, bool real_field, double theta, bool legacy = false)
{
  mfem::ND_FECollection fec(2, 3);
  FiniteElementSpace fespace(mesh, &fec);
  const MPI_Comm comm = fespace.GetComm();
  mfem::VectorFunctionCoefficient er(3, FieldRe), ei(3, FieldIm);
  mfem::ParGridFunction gr(&fespace.Get()), gi(&fespace.Get());
  gr.ProjectCoefficient(er);
  gi.ProjectCoefficient(ei);
  ComplexVector x(fespace.GetTrueVSize());
  gr.GetTrueDofs(x.Real());
  gi.GetTrueDofs(x.Imag());
  if (real_field)
  {
    x.Imag() = 0.0;
  }
  x *= std::exp(1i * theta);

  Normalized res;
  if (legacy)
  {
    linalg::NormalizePhase(comm, x);
    res.ref = -2;
  }
  else
  {
    const PhaseReference phase_ref(fespace);
    res.ref = phase_ref.Normalize(x);
  }
  for (int j = 0; j < 3; j++)
  {
    mfem::VectorFunctionCoefficient gj(3, [j](const mfem::Vector &p, mfem::Vector &v)
                                       { TestField(j, p, v); });
    mfem::ParLinearForm lf(&fespace.Get());
    lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(gj));
    lf.Assemble();
    std::unique_ptr<mfem::HypreParVector> lt(lf.ParallelAssemble());
    const Vector &l = *lt;
    res.g.emplace_back(linalg::Dot(comm, l, x.Real()), linalg::Dot(comm, l, x.Imag()));
  }
  res.imag = linalg::Norml2(comm, x.Imag()) / linalg::Norml2(comm, x);
  return res;
}

double MaxDiff(const Normalized &a, const Normalized &b)
{
  double diff = 0.0, scale = 0.0;
  for (std::size_t j = 0; j < a.g.size(); j++)
  {
    diff = std::max(diff, std::abs(a.g[j] - b.g[j]));
    scale = std::max(scale, std::abs(a.g[j]));
  }
  return diff / scale;
}

}  // namespace

TEST_CASE("Eigenvector phase reference fields", "[phase][Serial]")
{
  // Every component of the first reference field has a projection on every product of
  // coordinates over the cube [-1, 1]³ (one term of each parity under reflections).
  mfem::IntegrationRules rules(0, mfem::Quadrature1D::GaussLegendre);
  const auto &ir = rules.Get(mfem::Geometry::CUBE, 8);
  for (int a = 0; a < 3; a++)
  {
    for (int t = 0; t < 8; t++)
    {
      double s = 0.0;
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        const double xi[3] = {2.0 * ip.x - 1.0, 2.0 * ip.y - 1.0, 2.0 * ip.z - 1.0};
        double w[3], m = 1.0;
        PhaseReference::EvalField(0, xi, 3, w);
        for (int d = 0; d < 3; d++)
        {
          if (t & (1 << d))
          {
            m *= xi[d];
          }
        }
        s += ip.weight * w[a] * m;
      }
      CHECK(s > 0.01);  // At least 0.5 / 3^|t| (mean over the cube)
    }
  }

  // A reference integral sets the phase when above δ ‖b‖ ‖x‖.
  CHECK(PhaseReference::Accept({1.0, 0.5}, 2.0, 1.0));
  CHECK(PhaseReference::Accept({0.0, 3.0e-6}, 2.0, 1.0));
  CHECK(!PhaseReference::Accept({1.0e-6, 1.0e-6}, 2.0, 1.0));
  CHECK(!PhaseReference::Accept({0.0, 0.0}, 2.0, 1.0));
}

TEST_CASE("Eigenvector phase reference fallbacks", "[phase][Serial][Parallel]")
{
  // A vector orthogonal to the linear form of the first reference field takes its phase
  // from the second, one orthogonal to both from the sum of the degrees of freedom.
  auto smesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::TETRAHEDRON);
  Mesh mesh(Mpi::World(), smesh);
  mfem::ND_FECollection fec(2, 3);
  FiniteElementSpace fespace(mesh, &fec);
  const MPI_Comm comm = fespace.GetComm();
  const PhaseReference phase_ref(fespace);
  mfem::VectorFunctionCoefficient er(3, FieldRe), ei(3, FieldIm);
  mfem::ParGridFunction gr(&fespace.Get()), gi(&fespace.Get());
  gr.ProjectCoefficient(er);
  gi.ProjectCoefficient(ei);
  ComplexVector x(fespace.GetTrueVSize());
  gr.GetTrueDofs(x.Real());
  gi.GetTrueDofs(x.Imag());
  auto Orthogonalize = [&](Vector &v, const Vector &q)
  { v.Add(-linalg::Dot(comm, q, v) / linalg::Dot(comm, q, q), q); };
  ComplexVector y(x);
  CHECK(phase_ref.Normalize(y) == 0);
  const Vector &b0 = phase_ref.GetLinearForm(0);
  Vector q1(phase_ref.GetLinearForm(1));
  Orthogonalize(q1, b0);
  for (Vector *v : {&x.Real(), &x.Imag()})
  {
    Orthogonalize(*v, b0);
  }
  y = x;
  CHECK(phase_ref.Normalize(y) == 1);
  for (Vector *v : {&x.Real(), &x.Imag()})
  {
    Orthogonalize(*v, q1);
  }
  y = x;
  CHECK(phase_ref.Normalize(y) == PhaseReference::num_fields);
  const std::complex<double> sum = linalg::Sum(comm, y);
  CHECK(std::abs(sum.imag()) <= 1.0e-12 * std::abs(sum));
  CHECK(sum.real() > 0.0);
  y = 0.0;
  CHECK(phase_ref.Normalize(y) == -1);
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Eigenvector phase is independent of basis, element order and partition",
                 "[phase][Serial][Parallel]")
{
  // The same field, multiplied by different phases, on the shuffled mesh of one process
  // (the reference), and distributed with and without the local reordering of the elements
  // (which changes the orientations of faces and the order of the degrees of freedom).
  const auto mesh_path = WriteShuffledBox(temp_dir);
  for (const bool real_field : {false, true})
  {
    mfem::Mesh serial(mesh_path.string());
    Mesh self_mesh(MPI_COMM_SELF, serial);
    const auto expected = NormalizeField(self_mesh, real_field, 2.5);
    CHECK(expected.ref == 0);
    if (real_field)
    {
      // A real field stays real, with a definite sign.
      CHECK(expected.imag < 1.0e-12);
    }
    for (const bool local_reorder : {false, true})
    {
      auto iodata = MakeIoData(mesh_path, local_reorder);
      Mesh mesh(mesh::ReadMesh(iodata, Mpi::World()));
      const auto res = NormalizeField(mesh, real_field, local_reorder ? -1.7 : 0.3);
      CHECK(res.ref == expected.ref);
      CHECK(MaxDiff(res, expected) < 1.0e-11);
      if (real_field)
      {
        CHECK(res.imag < 1.0e-12);
      }
      if (local_reorder && !real_field)
      {
        // The reordering changes the basis: the mean of the degrees of freedom, the
        // previous convention, gives another phase (for a real field, at most a sign).
        const auto legacy = NormalizeField(mesh, real_field, -1.7, true);
        const auto legacy_expected = NormalizeField(self_mesh, real_field, 2.5, true);
        CHECK(MaxDiff(legacy, legacy_expected) > 1.0e-8);
      }
    }
  }
}

}  // namespace palace
