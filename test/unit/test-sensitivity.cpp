// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <fmt/format.h>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "driver.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/shapederivativeform.hpp"
#include "fixtures.hpp"
#include "linalg/divfree.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/rap.hpp"
#include "linalg/vector.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/sensitivity.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"
#include "utils/timer.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// Refine about half of the elements, chosen by their centers rather than at random, so
// that the refinement depends neither on the partition nor on earlier tests.
void RefineHalf(mfem::ParMesh &pmesh)
{
  constexpr double w[3] = {12.9898, 78.233, 37.719};
  mfem::Array<int> marked;
  mfem::Vector center;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    pmesh.GetElementCenter(e, center);
    double s = 0.0;
    for (int d = 0; d < center.Size(); d++)
    {
      s += w[d] * center(d);
    }
    const double r = std::sin(s) * 43758.5453;
    if (r - std::floor(r) < 0.5)
    {
      marked.Append(e);
    }
  }
  pmesh.GeneralRefinement(marked);
}

// A unit cube of two halves, with attribute 1 below z = 0.5 and 2 above, and the boundary
// attributes of mfem::Mesh::MakeCartesian3D. The nodes are moved by a smooth map which
// keeps the faces planar, so the elements are curved. With amr, the mesh is refined
// nonconformally, with hanging nodes.
std::vector<std::unique_ptr<Mesh>> CubeMesh(IoData &iodata, mfem::Element::Type type,
                                            bool amr)
{
  auto smesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(4, 4, 4, type, 1.0, 1.0, 1.0));
  for (int e = 0; e < smesh->GetNE(); e++)
  {
    mfem::Vector center;
    smesh->GetElementCenter(e, center);
    smesh->SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  smesh->SetAttributes();
  smesh->SetCurvature(2);
  const auto &nodal_fespace = *smesh->GetNodes()->FESpace();
  auto *nodes = smesh->GetNodes()->HostReadWrite();
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    double x[3];
    for (int d = 0; d < 3; d++)
    {
      x[d] = nodes[nodal_fespace.DofToVDof(i, d)];
    }
    for (int d = 0; d < 3; d++)
    {
      nodes[nodal_fespace.DofToVDof(i, d)] +=
          0.12 * x[d] * (1.0 - x[d]) * std::sin(2.0 * x[(d + 1) % 3] + x[(d + 2) % 3] + d);
    }
  }
  iodata.NondimensionalizeInputs(smesh);
  iodata.CheckConfiguration();
  if (amr)
  {
    smesh->EnsureNCMesh(true);
  }
  if (Mpi::Size(Mpi::World()) > smesh->GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  auto pmesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *smesh);
  if (amr)
  {
    RefineHalf(*pmesh);
  }
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(pmesh)));
  return mesh;
}

// The cube with two dielectrics, the upper one anisotropic. Terminals 1 and 2 are the faces
// z = 0 and z = 1, the faces y = 0 and y = 1 are grounded, and the others have zero charge.
struct Problem
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;

  Problem(mfem::Element::Type type, int order, bool amr) : iodata(Units(1.0, 1.0))
  {
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::ELECTROSTATIC;
    auto &lower = iodata.domains.materials.emplace_back();
    lower.attributes = {1};
    auto &upper = iodata.domains.materials.emplace_back();
    upper.attributes = {2};
    upper.epsilon_r.s = {2.0, 3.0, 4.0};
    iodata.boundaries.pec.attributes = {2, 4};
    iodata.boundaries.terminal[1].attributes = {1};
    iodata.boundaries.terminal[2].attributes = {6};
    iodata.solver.order = order;
    iodata.solver.linear.tol = 1.0e-13;
    iodata.solver.linear.max_it = 1000;
    mesh = CubeMesh(iodata, type, amr);
  }
};

// The cube as an eigenmode problem with every kind of term of T(λ) = K + λ C + λ² M +
// A2(λ): an anisotropic lossy dielectric below, an anisotropic magnetic conductor with a
// London depth above, and PEC, impedance, lumped port, second-order farfield and finite
// conductivity boundaries.
struct EigenProblem
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;

  EigenProblem(mfem::Element::Type type, int order, bool amr,
               const std::function<void(IoData &)> &modify = {})
    : iodata(Units(1.0, 1.0))
  {
    iodata.model.L0 = 1.0;  // m, so that every term is of order one
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::EIGENMODE;
    auto &lower = iodata.domains.materials.emplace_back();
    lower.attributes = {1};
    lower.epsilon_r.s = {2.0, 3.0, 4.0};
    lower.tandelta.s = {0.01, 0.02, 0.03};
    auto &upper = iodata.domains.materials.emplace_back();
    upper.attributes = {2};
    upper.mu_r.s = {1.5, 2.0, 2.5};
    upper.sigma.s = {1.0e-3, 2.0e-3, 3.0e-3};
    upper.lambda_L = 0.3;
    iodata.boundaries.pec.attributes = {2};
    auto &impedance = iodata.boundaries.impedance.emplace_back();
    impedance.attributes = {4};
    impedance.Rs = 100.0;
    impedance.Ls = 1.0e-6;
    impedance.Cs = 1.0e-11;
    auto &port = iodata.boundaries.lumpedport[1];
    port.R = 50.0;
    port.L = 1.0e-6;
    port.C = 1.0e-11;
    auto &element = port.elements.emplace_back();
    element.direction = {0.0, 1.0, 0.0};
    element.attributes = {3};
    iodata.boundaries.farfield.order = 2;
    iodata.boundaries.farfield.attributes = {6};
    auto &conductivity = iodata.boundaries.conductivity.emplace_back();
    conductivity.attributes = {5};
    conductivity.sigma = 1.0e3;
    iodata.solver.order = order;
    if (modify)
    {
      modify(iodata);
    }
    mesh = CubeMesh(iodata, type, amr);
  }
};

// The cube as a magnetostatic problem with an anisotropic magnetic material above, with a
// London depth if given. The faces x = 0, z = 0 and x = 1 are a connected PEC, and the
// faces z = 1 and y = 0 are surface current sources 1 and 2 along x, which close current
// loops through it.
struct MagnetostaticProblem
{
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;

  MagnetostaticProblem(mfem::Element::Type type, int order, bool amr, double lambda_L = 0.0)
    : iodata(Units(1.0, 1.0))
  {
    iodata.model.Lc = 1.0;
    iodata.problem.type = ProblemType::MAGNETOSTATIC;
    auto &lower = iodata.domains.materials.emplace_back();
    lower.attributes = {1};
    auto &upper = iodata.domains.materials.emplace_back();
    upper.attributes = {2};
    upper.mu_r.s = {1.5, 2.0, 2.5};
    upper.lambda_L = lambda_L;
    iodata.boundaries.pec.attributes = {1, 3, 5};
    for (const auto &[idx, attr] : {std::pair{1, 6}, std::pair{2, 2}})
    {
      auto &element = iodata.boundaries.current[idx].elements.emplace_back();
      element.direction = {1.0, 0.0, 0.0};
      element.attributes = {attr};
    }
    iodata.solver.order = order;
    iodata.solver.linear.tol = 1.0e-13;
    iodata.solver.linear.max_it = 1000;
    mesh = CubeMesh(iodata, type, amr);
  }
};

// The terminal potentials, as in ElectrostaticSolver::Solve.
std::vector<Vector> SolveTerminals(const IoData &iodata, LaplaceOperator &laplace_op)
{
  auto K = laplace_op.GetStiffnessMatrix();
  KspSolver ksp(iodata, laplace_op.GetH1Spaces());
  ksp.SetOperators(*K, *K);
  std::vector<Vector> V;
  Vector RHS;
  for (const auto &[idx, data] : laplace_op.GetSources())
  {
    auto &v = V.emplace_back();
    laplace_op.GetExcitationVector(idx, *K, v, RHS);
    ksp.Mult(RHS, v);
  }
  return V;
}

// The vector potentials of the surface current sources for unit currents, and the
// inductances M_ij = A_jᵀ K A_i for (i, j) = (0, 0), (0, 1) and (1, 1), as in
// MagnetostaticSolver.
std::pair<std::vector<Vector>, std::vector<double>>
SolveInductances(const IoData &iodata, CurlCurlOperator &curlcurl_op)
{
  auto K = curlcurl_op.GetStiffnessMatrix();
  auto P = curlcurl_op.GetPreconditionerMatrix();
  KspSolver ksp(iodata, curlcurl_op.GetNDSpaces(), &curlcurl_op.GetH1Spaces());
  ksp.SetOperators(*K, P ? *P : *K);
  std::vector<Vector> A;
  Vector RHS;
  for (const auto &[idx, data] : curlcurl_op.GetSurfaceCurrentOp())
  {
    REQUIRE(data.GetExcitationCurrent() == 1.0);
    curlcurl_op.GetCurrentExcitationVector(idx, RHS);
    auto &a = A.emplace_back(RHS.Size());
    a.UseDevice(true);
    a = 0.0;
    ksp.Mult(RHS, a);
  }
  REQUIRE(A.size() == 2);
  std::vector<double> M;
  Vector y(RHS.Size());
  y.UseDevice(true);
  for (const auto &[i, j] : {std::pair{0, 0}, std::pair{0, 1}, std::pair{1, 1}})
  {
    K->Mult(A[i], y);
    M.push_back(linalg::Dot<Vector>(Mpi::World(), A[j], y));
  }
  return {A, M};
}

// The capacitance C_ij = V_jᵀ K V_i, as in ElectrostaticSolver::PostprocessTerminals.
double Capacitance(const LaplaceOperator &laplace_op, const std::vector<int> &attributes,
                   const Vector &Vi, const Vector &Vj)
{
  const auto &h1_fespace = laplace_op.GetH1Space();
  Vector u(h1_fespace.GetVSize()), w(h1_fespace.GetVSize()), y(h1_fespace.GetVSize());
  h1_fespace.GetProlongationMatrix()->Mult(Vi, u);
  h1_fespace.GetProlongationMatrix()->Mult(Vj, w);
  laplace_op.GetStiffness(attributes)->Mult(u, y);
  return linalg::Dot<Vector>(laplace_op.GetComm(), w, y);
}

// A random complex vector on the true DOFs of the H(curl) space, zero on its essential true
// DOFs.
ComplexVector RandomComplexVector(const SpaceOperator &space_op, std::mt19937 &gen)
{
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  ComplexVector x(space_op.GetNDSpace().GetTrueVSize());
  for (int i = 0; i < x.Size(); i++)
  {
    x.Real()(i) = dist(gen);
    x.Imag()(i) = dist(gen);
  }
  for (const int i : space_op.GetNDDbcTDofLists().back())
  {
    x.Real()(i) = x.Imag()(i) = 0.0;
  }
  return x;
}

// The unconjugated product xᵀ A x.
std::complex<double> TransposeProduct(const ComplexOperator &A, const ComplexVector &x)
{
  ComplexVector y(x.Size());
  A.Mult(x, y);
  std::complex<double> dot = {
      linalg::LocalDot(x.Real(), y.Real()) - linalg::LocalDot(x.Imag(), y.Imag()),
      linalg::LocalDot(x.Real(), y.Imag()) + linalg::LocalDot(x.Imag(), y.Real())};
  Mpi::GlobalSum(1, &dot, Mpi::World());
  return dot;
}

// Move the mesh nodes by h W, for W on the true DOFs of the mesh nodal space.
void MoveNodes(Mesh &mesh, const mfem::ParFiniteElementSpace &nodal_fespace,
               const Vector &W, double h)
{
  mfem::Vector W_l(nodal_fespace.GetVSize());
  nodal_fespace.GetProlongationMatrix()->Mult(W, W_l);
  mesh.Get().GetNodes()->Add(h, W_l);
  mesh.Get().NodesUpdated();
  mesh.ResetCeedObjects();
}

// A velocity field, at a point in mesh length units.
using VelocityField = std::function<std::array<double, 3>(const std::array<double, 3> &)>;

// A velocity field at the nodes of a serial mesh with nodes.
mfem::Vector NodalVelocity(const mfem::Mesh &mesh, const VelocityField &V)
{
  const auto &fespace = *mesh.GetNodes()->FESpace();
  const auto *nodes = mesh.GetNodes()->HostRead();
  mfem::Vector v(mesh.GetNodes()->Size());
  for (int i = 0; i < fespace.GetNDofs(); i++)
  {
    std::array<double, 3> x;
    for (int d = 0; d < 3; d++)
    {
      x[d] = nodes[fespace.DofToVDof(i, d)];
    }
    const auto vi = V(x);
    for (int d = 0; d < 3; d++)
    {
      v(fespace.DofToVDof(i, d)) = vi[d];
    }
  }
  return v;
}

// Write a serial mesh with nodes moved by t v, and the table of the velocity v at its
// nodes, from the root process.
void WriteMesh(mfem::Mesh &mesh, const mfem::Vector &v, double t,
               const std::filesystem::path &path)
{
  if (Mpi::Root(Mpi::World()))
  {
    auto &nodes = *mesh.GetNodes();
    const mfem::Vector X(nodes);
    nodes.Add(t, v);
    std::ofstream fo(path);
    fo.precision(17);
    mesh.Print(fo);
    nodes = X;
  }
  Mpi::Barrier(Mpi::World());
}

void WriteVelocityTable(const mfem::Mesh &mesh, const std::string &name,
                        const mfem::Vector &v, const std::filesystem::path &path)
{
  if (Mpi::Root(Mpi::World()))
  {
    const auto &fespace = *mesh.GetNodes()->FESpace();
    const auto *nodes = mesh.GetNodes()->HostRead();
    std::ofstream fo(path);
    fo.precision(17);
    fo << fmt::format("x, y, z, {0}.x, {0}.y, {0}.z\n", name);
    for (int i = 0; i < fespace.GetNDofs(); i++)
    {
      for (int d = 0; d < 3; d++)
      {
        fo << nodes[fespace.DofToVDof(i, d)] << ", ";
      }
      for (int d = 0; d < 3; d++)
      {
        fo << v(fespace.DofToVDof(i, d)) << ((d < 2) ? ", " : "\n");
      }
    }
  }
  Mpi::Barrier(Mpi::World());
}

// Run Palace, returning its output folder.
std::filesystem::path RunPalace(const nlohmann::json &config)
{
  IoData iodata(config, false);
  MPI_Comm comm = Mpi::World();
  MakeOutputFolder(iodata, comm);
  BlockTimer::Reset();
  palace::Run(iodata, comm, 1);
  return iodata.problem.output;
}

// A column of an output CSV file.
std::vector<double> ReadColumn(const std::filesystem::path &file, const std::string &column)
{
  std::ifstream fi(file);
  REQUIRE(fi.is_open());
  std::stringstream buffer;
  buffer << fi.rdbuf();
  Table table(buffer.str());
  for (std::size_t j = 0; j < table.n_cols(); j++)
  {
    if (table[j].header_text == column)
    {
      return table[j].data;
    }
  }
  FAIL("No column \"" << column << "\" in " << file);
  return {};
}

// The sensitivities of a Palace run.
nlohmann::json ReadSensitivities(const std::filesystem::path &output)
{
  std::ifstream fi(output / "sensitivity.json");
  REQUIRE(fi.is_open());
  return nlohmann::json::parse(fi);
}

}  // namespace

TEST_CASE("Electrostatic capacitance sensitivities", "[sensitivity][Serial][Parallel]")
{
  // The derivatives of the capacitances need no adjoint solve. Check them against finite
  // differences of the solved capacitances on moved meshes.
  auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  auto order = GENERATE(1, 2);
  auto amr = GENERATE(false, true);
  INFO("Element type: " << type << ", order: " << order << ", AMR: " << amr);
  Problem problem(type, order, amr);
  auto &mesh = *problem.mesh[0];

  LaplaceOperator laplace_op(problem.iodata, problem.mesh);
  const auto V = SolveTerminals(problem.iodata, laplace_op);
  REQUIRE(V.size() == 2);
  const auto dK = laplace_op.GetStiffnessShapeDerivative();
  const auto &nodal_fespace = dK->GetNodalSpace();
  std::mt19937 gen(1234 + Mpi::Rank(Mpi::World()));
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Vector W(nodal_fespace.GetTrueVSize());
  for (int i = 0; i < W.Size(); i++)
  {
    W(i) = dist(gen);
  }

  const std::vector<std::pair<int, int>> entries = {{0, 0}, {0, 1}, {1, 1}};
  std::vector<double> adjoint;
  Vector g;
  for (const auto &[i, j] : entries)
  {
    dK->Mult(V[i], V[j], g);
    adjoint.push_back(linalg::Dot<Vector>(Mpi::World(), g, W));
  }

  auto CapacitancesAt = [&](double h)
  {
    MoveNodes(mesh, nodal_fespace, W, h);
    LaplaceOperator moved_op(problem.iodata, problem.mesh);
    const auto Vh = SolveTerminals(problem.iodata, moved_op);
    std::vector<double> C;
    for (const auto &[i, j] : entries)
    {
      C.push_back(Capacitance(moved_op, {1, 2}, Vh[i], Vh[j]));
    }
    MoveNodes(mesh, nodal_fespace, W, -h);
    return C;
  };
  auto CentralDifference = [&](double h)
  {
    const auto Cp = CapacitancesAt(h), Cm = CapacitancesAt(-h);
    std::vector<double> d;
    for (std::size_t k = 0; k < Cp.size(); k++)
    {
      d.push_back((Cp[k] - Cm[k]) / (2.0 * h));
    }
    return d;
  };
  const double h = 2.0e-4;
  const auto d1 = CentralDifference(h), d2 = CentralDifference(0.5 * h);
  for (std::size_t k = 0; k < entries.size(); k++)
  {
    const double fd = (4.0 * d2[k] - d1[k]) / 3.0;
    INFO("Entry (" << entries[k].first << ", " << entries[k].second
                   << "): g ⋅ W = " << adjoint[k] << ", finite difference = " << fd);
    CHECK(std::abs(adjoint[k] - fd) <= 1.0e-7 * std::abs(fd));
  }
}

TEST_CASE("Eigenmode system shape derivatives", "[sensitivity][Serial][Parallel]")
{
  // The gradient of xᵀ T(λ) x for a complex λ and vector x, against finite differences of
  // the assembled form on moved meshes.
  auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  auto order = GENERATE(1, 2);
  auto amr = GENERATE(false, true);
  INFO("Element type: " << type << ", order: " << order << ", AMR: " << amr);
  EigenProblem problem(type, order, amr);
  auto &mesh = *problem.mesh[0];
  SpaceOperator space_op(problem.iodata, problem.mesh);
  const std::complex<double> lambda(0.05, 1.3),
      omega = lambda / std::complex<double>(0.0, 1.0);

  std::mt19937 gen(1234 + Mpi::Rank(Mpi::World()));
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  const auto x = RandomComplexVector(space_op, gen);
  auto Form = [&]()
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
    REQUIRE((K && C && M && A2));
    return TransposeProduct(*space_op.GetSystemMatrix(std::complex<double>(1.0), lambda,
                                                      lambda * lambda, K.get(), C.get(),
                                                      M.get(), A2.get()),
                            x);
  };

  const auto dT = space_op.GetSystemShapeDerivative(1.0, lambda, lambda * lambda, omega);
  Vector gr, gi;
  GetComplexShapeGradient(*dT[0], *dT[1], x, gr, gi);
  const auto &nodal_fespace = dT[0]->GetNodalSpace();
  auto Contract = [&](const Vector &W)
  {
    return std::complex<double>(linalg::Dot<Vector>(Mpi::World(), gr, W),
                                linalg::Dot<Vector>(Mpi::World(), gi, W));
  };
  Vector W(nodal_fespace.GetTrueVSize());
  for (int i = 0; i < W.Size(); i++)
  {
    W(i) = dist(gen);
  }
  auto CentralDifference = [&](double h)
  {
    MoveNodes(mesh, nodal_fespace, W, h);
    const auto Fp = Form();
    MoveNodes(mesh, nodal_fespace, W, -2.0 * h);
    const auto Fm = Form();
    MoveNodes(mesh, nodal_fespace, W, h);
    return (Fp - Fm) / (2.0 * h);
  };
  const double h = 2.0e-4;
  const auto fd = (4.0 * CentralDifference(0.5 * h) - CentralDifference(h)) / 3.0;
  const auto adjoint = Contract(W);
  INFO("xᵀ T x = " << Form() << ", g ⋅ W = " << adjoint << ", finite difference = " << fd);
  CHECK(std::abs(adjoint - fd) <= 1.0e-7 * std::abs(fd));
}

TEST_CASE("Eigenmode system material and lumped port derivatives",
          "[sensitivity][Serial][Parallel]")
{
  // The derivatives of xᵀ T(λ) x with respect to relative changes of the material
  // properties and of the lumped port elements, against finite differences of the form with
  // the property scaled in the configuration.
  auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  INFO("Element type: " << type);
  const std::complex<double> lambda(0.05, 1.3),
      omega = lambda / std::complex<double>(0.0, 1.0);
  EigenProblem problem(type, 1, false);
  SpaceOperator space_op(problem.iodata, problem.mesh);
  std::mt19937 gen(1234 + Mpi::Rank(Mpi::World()));
  const auto x = RandomComplexVector(space_op, gen);
  auto Product = [&x](const ComplexOperator &A) { return TransposeProduct(A, x); };
  auto CentralDifference = [&](const std::function<void(IoData &, double)> &scale)
  {
    auto Form = [&](double s)
    {
      EigenProblem scaled(type, 1, false, [&](IoData &iodata) { scale(iodata, s); });
      SpaceOperator op(scaled.iodata, scaled.mesh);
      auto K = op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      auto C = op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      auto M = op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
      auto A2 = op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
      return Product(*op.GetSystemMatrix(std::complex<double>(1.0), lambda, lambda * lambda,
                                         K.get(), C.get(), M.get(), A2.get()));
    };
    auto Difference = [&](double h) { return (Form(1.0 + h) - Form(1.0 - h)) / (2.0 * h); };
    const double h = 1.0e-2;
    return std::pair{(4.0 * Difference(0.5 * h) - Difference(h)) / 3.0,
                     std::abs(Form(1.0))};
  };
  auto Check =
      [](std::complex<double> derivative, const std::pair<std::complex<double>, double> &fd)
  {
    // The finite differences of the form F have a rounding error of about ε |F| / h.
    INFO("derivative = " << derivative << ", finite difference = " << fd.first
                         << ", |xᵀ T x| = " << fd.second);
    CHECK(std::abs(derivative - fd.first) <=
          1.0e-8 * std::abs(fd.first) + 1.0e-12 * fd.second);
  };

  // Each material property, as a scaling of its configuration.
  auto Scale = [](auto &tensor, double s)
  {
    for (auto &v : tensor.s)
    {
      v *= s;
    }
  };
  const std::vector<std::pair<MaterialProperty, std::vector<int>>> properties = {
      {MaterialProperty::PERMITTIVITY, {0, 1}},
      {MaterialProperty::LOSS_TANGENT, {0}},
      {MaterialProperty::CONDUCTIVITY, {1}},
      {MaterialProperty::PERMEABILITY, {0, 1}},
      {MaterialProperty::LONDON_DEPTH, {1}}};
  for (const auto &entry : properties)
  {
    // The lambda below captures the property, which may not be a structured binding with
    // OpenMP.
    const MaterialProperty property = entry.first;
    for (const int k : entry.second)
    {
      INFO("Property " << static_cast<int>(property) << " of material " << k);
      const auto dT = space_op.GetSystemMaterialDerivative(
          property, problem.iodata.domains.materials[k].attributes, 1.0, lambda,
          lambda * lambda, omega);
      REQUIRE(dT);
      Check(Product(*dT), CentralDifference(
                              [&, property, k](IoData &iodata, double s)
                              {
                                auto &material = iodata.domains.materials[k];
                                switch (property)
                                {
                                  case MaterialProperty::PERMITTIVITY:
                                    Scale(material.epsilon_r, s);
                                    break;
                                  case MaterialProperty::LOSS_TANGENT:
                                    Scale(material.tandelta, s);
                                    break;
                                  case MaterialProperty::CONDUCTIVITY:
                                    Scale(material.sigma, s);
                                    break;
                                  case MaterialProperty::PERMEABILITY:
                                    Scale(material.mu_r, s);
                                    break;
                                  case MaterialProperty::LONDON_DEPTH:
                                    material.lambda_L *= s;
                                    break;
                                }
                              }));
    }
  }

  // The resistance, inductance and capacitance of the lumped port.
  const auto dT = space_op.GetLumpedPortDerivatives(1, -1, 1.0, lambda, lambda * lambda);
  for (int e = 0; e < 3; e++)
  {
    INFO("Lumped port element " << e);
    REQUIRE(dT[e]);
    Check(Product(*dT[e]), CentralDifference(
                               [e](IoData &iodata, double s)
                               {
                                 auto &port = iodata.boundaries.lumpedport.at(1);
                                 (e == 0 ? port.R : (e == 1 ? port.L : port.C)) *= s;
                               }));
  }
}

TEST_CASE("Divergence-free projection with losses", "[sensitivity][Serial][Parallel]")
{
  // The projection of the eigenvalue solve is with respect to the permittivity. With
  // losses, the correction for the sensitivities makes a vector divergence-free with
  // respect to B = M + C / λ instead, changing it by a gradient.
  EigenProblem problem(mfem::Element::TETRAHEDRON, 2, false);
  SpaceOperator space_op(problem.iodata, problem.mesh);
  DivFreeSolver<ComplexVector> divfree(space_op.GetMaterialOp(), space_op.GetNDSpace(),
                                       space_op.GetH1Spaces(),
                                       space_op.GetAuxBdrTDofLists(), 1.0e-14, 1000, 0);
  const std::complex<double> lambda(0.05, 1.3);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  const auto B =
      BuildParSumOperator(std::vector<std::complex<double>>{1.0, 1.0 / lambda},
                          std::vector<const ComplexOperator *>{M.get(), C.get()});

  std::mt19937 gen(1234 + Mpi::Rank(Mpi::World()));
  auto x = RandomComplexVector(space_op, gen);
  const int size = x.Size();
  ComplexVector y(size), By(size);
  divfree.Mult(x);
  y = x;
  divfree.Mult(*B, y);

  // The divergence Gᵀ B y vanishes away from the essential boundaries of the potential.
  const auto &Grad = space_op.GetGradMatrix();
  auto Divergence = [&](const ComplexVector &v)
  {
    B->Mult(v, By);
    ComplexVector div(Grad.Width());
    Grad.MultTranspose(By.Real(), div.Real());
    Grad.MultTranspose(By.Imag(), div.Imag());
    linalg::SetSubVector(div, space_op.GetAuxBdrTDofLists().back(), 0.0);
    return linalg::Norml2(Mpi::World(), div);
  };
  const double div_x = Divergence(x);
  CHECK(div_x > 0.0);
  CHECK(Divergence(y) <= 1.0e-10 * div_x);

  // The correction is a gradient, so the curl is unchanged.
  const auto &Curl = space_op.GetCurlMatrix();
  ComplexVector d(size), curl(Curl.Height()), curl_x(Curl.Height());
  d = y;
  d -= x;
  Curl.Mult(d.Real(), curl.Real());
  Curl.Mult(d.Imag(), curl.Imag());
  Curl.Mult(x.Real(), curl_x.Real());
  Curl.Mult(x.Imag(), curl_x.Imag());
  CHECK(linalg::Norml2(Mpi::World(), d) > 0.0);
  CHECK(linalg::Norml2(Mpi::World(), curl) <=
        1.0e-12 * linalg::Norml2(Mpi::World(), curl_x));
}

TEST_CASE("Magnetostatic inductance sensitivities", "[sensitivity][Serial][Parallel]")
{
  // The derivatives of the inductances need no adjoint solve. Check them against finite
  // differences of the solved inductances on moved meshes, with the current densities of
  // the sources fixed.
  // Each element type, order, refinement and London depth, in two of four combinations.
  const auto params = GENERATE(table<mfem::Element::Type, int, bool, double>(
      {{mfem::Element::HEXAHEDRON, 1, false, 0.0},
       {mfem::Element::TETRAHEDRON, 2, false, 0.3},
       {mfem::Element::HEXAHEDRON, 2, true, 0.3},
       {mfem::Element::TETRAHEDRON, 1, true, 0.0}}));
  const auto type = std::get<0>(params);
  const int order = std::get<1>(params);
  const bool amr = std::get<2>(params);
  const double lambda_L = std::get<3>(params);
  INFO("Element type: " << type << ", order: " << order << ", AMR: " << amr
                        << ", London depth: " << lambda_L);
  MagnetostaticProblem problem(type, order, amr, lambda_L);
  auto &mesh = *problem.mesh[0];
  CurlCurlOperator curlcurl_op(problem.iodata, problem.mesh);

  auto Solve = [&]() { return SolveInductances(problem.iodata, curlcurl_op); };
  const auto [A, M] = Solve();

  const auto dK = curlcurl_op.GetStiffnessShapeDerivative();
  std::vector<std::unique_ptr<ShapeDerivativeOperator>> db;
  db.push_back(curlcurl_op.GetCurrentExcitationShapeDerivative(1));
  db.push_back(curlcurl_op.GetCurrentExcitationShapeDerivative(2));
  // A smooth displacement which keeps the source faces planar and their edges off the PEC
  // along x, as the sources must be for their excitations to be compatible.
  const auto &nodal_fespace = dK->GetNodalSpace();
  const auto *nodes = mesh.Get().GetNodes()->HostRead();
  mfem::Vector W_l(nodal_fespace.GetVSize());
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    double x[3];
    for (int d = 0; d < 3; d++)
    {
      x[d] = nodes[nodal_fespace.DofToVDof(i, d)];
    }
    const double w[3] = {std::sin(2.0 * x[1] + x[2]),
                         x[1] * (1.0 - x[1]) * std::cos(x[0] + 2.0 * x[2]),
                         (1.0 - x[2]) * std::sin(x[0] - x[1])};
    for (int d = 0; d < 3; d++)
    {
      W_l(nodal_fespace.DofToVDof(i, d)) = w[d];
    }
  }
  Vector W(nodal_fespace.GetTrueVSize());
  nodal_fespace.GetRestrictionMatrix()->Mult(W_l, W);
  std::vector<double> adjoint;
  Vector g, t;
  for (const auto &[i, j] : {std::pair{0, 0}, std::pair{0, 1}, std::pair{1, 1}})
  {
    // dM_ij = b_i'(A_j) + b_j'(A_i) - K'(A_i, A_j).
    dK->Mult(A[i], A[j], g);
    g *= -1.0;
    db[i]->Mult(A[j], A[j], t);
    g += t;
    db[j]->Mult(A[i], A[i], t);
    g += t;
    adjoint.push_back(linalg::Dot<Vector>(Mpi::World(), g, W));
  }

  auto CentralDifference = [&](double h)
  {
    MoveNodes(mesh, nodal_fespace, W, h);
    const auto Mp = Solve().second;
    MoveNodes(mesh, nodal_fespace, W, -2.0 * h);
    const auto Mm = Solve().second;
    MoveNodes(mesh, nodal_fespace, W, h);
    std::vector<double> d;
    for (std::size_t k = 0; k < Mp.size(); k++)
    {
      d.push_back((Mp[k] - Mm[k]) / (2.0 * h));
    }
    return d;
  };
  const double h = 2.0e-4;
  const auto d1 = CentralDifference(h), d2 = CentralDifference(0.5 * h);
  for (std::size_t k = 0; k < adjoint.size(); k++)
  {
    const double fd = (4.0 * d2[k] - d1[k]) / 3.0;
    INFO("Entry " << k << ": M = " << M[k] << ", g ⋅ W = " << adjoint[k]
                  << ", finite difference = " << fd);
    CHECK(std::abs(adjoint[k] - fd) <= 1.0e-7 * std::abs(fd));
  }
}

TEST_CASE("Magnetostatic material derivatives", "[sensitivity][Serial][Parallel]")
{
  // The derivatives of the inductances with respect to relative changes of the permeability
  // and the London depth of the upper material, as in MagnetostaticSolver, against finite
  // differences.
  auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  INFO("Element type: " << type);
  constexpr int order = 2;
  constexpr double lambda_L = 0.3;
  auto Inductances = [&](double mu_scale, double lambda_scale)
  {
    MagnetostaticProblem problem(type, order, false, lambda_L);
    auto &upper = problem.iodata.domains.materials[1];
    for (auto &mu : upper.mu_r.s)
    {
      mu *= mu_scale;
    }
    upper.lambda_L *= lambda_scale;
    CurlCurlOperator curlcurl_op(problem.iodata, problem.mesh);
    return SolveInductances(problem.iodata, curlcurl_op).second;
  };

  MagnetostaticProblem problem(type, order, false, lambda_L);
  CurlCurlOperator curlcurl_op(problem.iodata, problem.mesh);
  const auto [A, M] = SolveInductances(problem.iodata, curlcurl_op);
  const auto &attributes = problem.iodata.domains.materials[1].attributes;
  const auto K_m = curlcurl_op.GetStiffness(attributes);
  const auto L_m = curlcurl_op.GetLondonMass(attributes);
  const auto &nd_fespace = curlcurl_op.GetNDSpace();
  Vector u(nd_fespace.GetVSize()), w(nd_fespace.GetVSize()), y(nd_fespace.GetVSize());
  u.UseDevice(true);
  w.UseDevice(true);
  y.UseDevice(true);
  std::vector<double> d_mu, d_lambda;
  for (const auto &[i, j] : {std::pair{0, 0}, std::pair{0, 1}, std::pair{1, 1}})
  {
    nd_fespace.GetProlongationMatrix()->Mult(A[i], u);
    nd_fespace.GetProlongationMatrix()->Mult(A[j], w);
    K_m->Mult(u, y);
    d_mu.push_back(linalg::Dot<Vector>(Mpi::World(), w, y));
    L_m->Mult(u, y);
    d_lambda.push_back(2.0 * linalg::Dot<Vector>(Mpi::World(), w, y));
  }

  const double h = 1.0e-3;
  // Plain variables rather than structured bindings, which lambdas can't capture with all
  // compilers.
  for (const auto &property : {std::tuple{"Permeability", &d_mu, 1.0, 0.0},
                               std::tuple{"LondonDepth", &d_lambda, 0.0, 1.0}})
  {
    const char *name = std::get<0>(property);
    const auto *derivative = std::get<1>(property);
    const double mu_step = std::get<2>(property), lambda_step = std::get<3>(property);
    auto CentralDifference = [&](double h)
    {
      const auto Mp = Inductances(1.0 + mu_step * h, 1.0 + lambda_step * h);
      const auto Mm = Inductances(1.0 - mu_step * h, 1.0 - lambda_step * h);
      std::vector<double> d;
      for (std::size_t k = 0; k < Mp.size(); k++)
      {
        d.push_back((Mp[k] - Mm[k]) / (2.0 * h));
      }
      return d;
    };
    const auto d1 = CentralDifference(h), d2 = CentralDifference(0.5 * h);
    for (std::size_t k = 0; k < M.size(); k++)
    {
      const double fd = (4.0 * d2[k] - d1[k]) / 3.0;
      INFO(name << ", entry " << k << ": M = " << M[k]
                << ", derivative = " << (*derivative)[k] << ", finite difference = " << fd);
      CHECK(std::abs((*derivative)[k] - fd) <= 1.0e-7 * std::abs(M[k]));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Eigenmode sensitivities of a lossy cavity",
                 "[sensitivity][Serial][Parallel]")
{
  // f and Q of a cavity with two dielectrics of different loss tangents and a lumped port
  // given by R, L and C, against finite differences of Palace on morphed meshes. The
  // velocity stretches the port and moves the interface of the dielectrics, which changes Q
  // through the irrotational part of the eigenvectors.
#if !defined(MFEM_USE_GSLIB)
  SKIP("Design velocities require GSLIB");
#endif
  auto smesh =
      mfem::Mesh::MakeCartesian3D(4, 3, 4, mfem::Element::HEXAHEDRON, 1.0, 0.8, 0.6);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center;
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.3) ? 1 : 2);
  }
  smesh.SetAttributes();
  smesh.EnsureNodes();
  // The port is the face x = 1 along y, which stays planar.
  const auto v = NodalVelocity(
      smesh,
      [](const std::array<double, 3> &x) -> std::array<double, 3>
      {
        return {0.3 * x[0] * (1.0 - x[0]) * std::sin(2.0 * x[1] + x[2]),
                0.2 * x[1] * (1.0 + x[2]),
                0.15 * x[2] * (1.0 - x[2] / 0.6) * std::cos(x[0] + x[1]) + 0.1 * x[2]};
      });
  WriteVelocityTable(smesh, "p", v, temp_dir / "velocities.csv");

  int count = 0;
  auto Config = [&](double t, bool sensitivity)
  {
    const auto path = temp_dir / fmt::format("cavity-{:d}.mesh", count);
    WriteMesh(smesh, v, t, path);
    nlohmann::json config = {
        {"Problem",
         {{"Type", "Eigenmode"},
          {"Output", (temp_dir / fmt::format("out-{:d}", count++)).string()}}},
        {"Model", {{"Mesh", path.string()}, {"L0", 1.0e-2}, {"Lc", 1.0}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"LossTan", 0.01}},
            {{"Attributes", {2}}, {"Permittivity", 4.0}, {"LossTan", 0.002}}}}}},
        {"Boundaries",
         {{"PEC", {{"Attributes", {1, 2, 4, 5, 6}}}},
          {"LumpedPort",
           {{{"Index", 1},
             {"R", 2000.0},
             {"L", 1.0e-9},
             {"C", 1.0e-13},
             {"Attributes", {3}},
             {"Direction", "+Y"}}}}}},
        {"Solver",
         {{"Order", 2},
          {"Eigenmode", {{"N", 2}, {"Target", 10.0}, {"Tol", 1.0e-11}}},
          {"Linear", {{"Tol", 1.0e-12}, {"MaxIts", 1000}}}}}};
    if (sensitivity)
    {
      config["Sensitivity"] = {{"Velocities", (temp_dir / "velocities.csv").string()}};
    }
    return config;
  };
  const auto output = RunPalace(Config(0.0, true));
  const std::array<const char *, 2> quantities = {"Frequency", "QualityFactor"},
                                    columns = {"Re{f} (GHz)", "Q"};
  const auto sensitivities = ReadSensitivities(output);

  // The derivatives of f and Q of each mode, by Richardson extrapolation of central
  // differences.
  auto CentralDifference = [&](double h)
  {
    const auto output_p = RunPalace(Config(h, false)),
               output_m = RunPalace(Config(-h, false));
    std::array<std::vector<double>, 2> d;
    for (int k = 0; k < 2; k++)
    {
      const auto Jp = ReadColumn(output_p / "eig.csv", columns[k]),
                 Jm = ReadColumn(output_m / "eig.csv", columns[k]);
      for (int m = 0; m < 2; m++)
      {
        d[k].push_back((Jp[m] - Jm[m]) / (2.0 * h));
      }
    }
    return d;
  };
  const double h = 1.0e-3;
  const auto d1 = CentralDifference(h), d2 = CentralDifference(0.5 * h);
  for (int k = 0; k < 2; k++)
  {
    const auto J = ReadColumn(output / "eig.csv", columns[k]);
    for (int m = 0; m < 2; m++)
    {
      const double fd = (4.0 * d2[k][m] - d1[k][m]) / 3.0;
      const double d = sensitivities.at(quantities[k]).at("Derivatives").at("p")[m];
      INFO(quantities[k] << " of mode " << m + 1 << ": " << J[m] << ", derivative " << d
                         << ", finite difference " << fd);
      CHECK(std::abs(d - fd) <= 1.0e-6 * std::abs(J[m]));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "Magnetostatic sensitivities at fixed currents",
                 "[sensitivity][Serial][Parallel]")
{
  // The inductances at fixed currents against finite differences of Palace on morphed
  // meshes, for a velocity which widens both sources, the faces z = 1 and y = 0 along x
  // between the PEC faces x = 0 and x = 1, and keeps them planar.
#if !defined(MFEM_USE_GSLIB)
  SKIP("Design velocities require GSLIB");
#endif
  auto smesh =
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center;
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  smesh.SetAttributes();
  smesh.EnsureNodes();
  const auto v =
      NodalVelocity(smesh,
                    [](const std::array<double, 3> &x) -> std::array<double, 3>
                    {
                      const double bubble =
                          x[0] * (1.0 - x[0]) * x[1] * (1.0 - x[1]) * x[2] * (1.0 - x[2]);
                      return {0.3 * x[0] * (1.0 - x[0]) * std::sin(2.0 * x[1] + x[2]),
                              0.3 * x[1] + 0.5 * bubble,
                              0.2 * x[2] + 0.5 * bubble * std::cos(x[0] + x[1])};
                    });
  WriteVelocityTable(smesh, "p", v, temp_dir / "velocities.csv");

  int count = 0;
  auto Config = [&](double t, bool sensitivity)
  {
    const auto path = temp_dir / fmt::format("cube-{:d}.mesh", count);
    WriteMesh(smesh, v, t, path);
    nlohmann::json config = {
        {"Problem",
         {{"Type", "Magnetostatic"},
          {"Output", (temp_dir / fmt::format("out-{:d}", count++)).string()}}},
        {"Model", {{"Mesh", path.string()}, {"L0", 1.0e-3}, {"Lc", 1.0}}},
        {"Domains",
         {{"Materials",
           {{{"Attributes", {1}}}, {{"Attributes", {2}}, {"Permeability", 2.0}}}}}},
        {"Boundaries",
         {{"PEC", {{"Attributes", {1, 3, 5}}}},
          {"SurfaceCurrent",
           {{{"Index", 1}, {"Attributes", {6}}, {"Direction", "+X"}},
            {{"Index", 2}, {"Attributes", {2}}, {"Direction", "+X"}}}}}},
        {"Solver", {{"Order", 2}, {"Linear", {{"Tol", 1.0e-12}, {"MaxIts", 1000}}}}}};
    if (sensitivity)
    {
      config["Sensitivity"] = {{"Velocities", (temp_dir / "velocities.csv").string()}};
    }
    return config;
  };
  const auto output = RunPalace(Config(0.0, true));
  const auto dM = ReadSensitivities(output).at("Inductance").at("Derivatives").at("p");

  // The derivatives of the columns M[i][j], by Richardson extrapolation of central
  // differences.
  auto Column = [](int j) { return fmt::format("M[i][{:d}] (H)", j + 1); };
  auto CentralDifference = [&](double h)
  {
    const auto output_p = RunPalace(Config(h, false)),
               output_m = RunPalace(Config(-h, false));
    std::array<std::vector<double>, 2> d;
    for (int j = 0; j < 2; j++)
    {
      const auto Mp = ReadColumn(output_p / "terminal-M.csv", Column(j)),
                 Mm = ReadColumn(output_m / "terminal-M.csv", Column(j));
      d[j] = {(Mp[0] - Mm[0]) / (2.0 * h), (Mp[1] - Mm[1]) / (2.0 * h)};
    }
    return d;
  };
  const double h = 1.0e-3;
  const auto d1 = CentralDifference(h), d2 = CentralDifference(0.5 * h);
  for (int j = 0; j < 2; j++)
  {
    const auto M = ReadColumn(output / "terminal-M.csv", Column(j));
    for (int i = 0; i < 2; i++)
    {
      const double fd = (4.0 * d2[j][i] - d1[j][i]) / 3.0;
      const double d = dM[i][j];
      INFO("M[" << i + 1 << "][" << j + 1 << "] = " << M[i] << " H, derivative " << d
                << " H, finite difference " << fd << " H");
      CHECK(std::abs(d - fd) <= 1.0e-6 * std::abs(M[i]));
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Input mesh nodes of a Gmsh file",
                 "[sensitivity][Serial]")
{
  // A third-order tetrahedron, whose nodes MFEM's reader moves from the uniformly spaced
  // points of the file to the Gauss-Lobatto points of its basis, and ReadInputMesh back.
  using Point = std::array<double, 3>;
  const std::array<Point, 4> v = {
      {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
  std::vector<Point> x(v.begin(), v.end());
  for (const auto &[a, b] : {std::pair{0, 1}, std::pair{1, 2}, std::pair{2, 0},
                             std::pair{3, 0}, std::pair{3, 2}, std::pair{3, 1}})
  {
    for (const double t : {1.0 / 3.0, 2.0 / 3.0})
    {
      x.push_back({v[a][0] + t * (v[b][0] - v[a][0]), v[a][1] + t * (v[b][1] - v[a][1]),
                   v[a][2] + t * (v[b][2] - v[a][2])});
    }
  }
  for (const auto &f :
       {std::array{0, 2, 1}, std::array{0, 1, 3}, std::array{0, 3, 2}, std::array{3, 1, 2}})
  {
    Point c = {0.0, 0.0, 0.0};
    for (const int k : f)
    {
      for (int d = 0; d < 3; d++)
      {
        c[d] += v[k][d] / 3.0;
      }
    }
    x.push_back(c);
  }
  const auto path = temp_dir / "tet20.msh";
  {
    std::ofstream fo(path);
    fo.precision(17);
    fo << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << x.size() << "\n";
    for (std::size_t k = 0; k < x.size(); k++)
    {
      fo << k + 1 << " " << x[k][0] << " " << x[k][1] << " " << x[k][2] << "\n";
    }
    fo << "$EndNodes\n$Elements\n1\n1 29 2 1 1";
    for (std::size_t k = 0; k < x.size(); k++)
    {
      fo << " " << k + 1;
    }
    fo << "\n$EndElements\n";
  }
  auto NodePoints = [](const mfem::Mesh &mesh)
  {
    const auto &fespace = *mesh.GetNodes()->FESpace();
    const auto *nodes = mesh.GetNodes()->HostRead();
    std::vector<Point> points;
    for (int i = 0; i < fespace.GetNDofs(); i++)
    {
      points.push_back({nodes[fespace.DofToVDof(i, 0)], nodes[fespace.DofToVDof(i, 1)],
                        nodes[fespace.DofToVDof(i, 2)]});
    }
    return points;
  };
  auto Same = [&](std::vector<Point> points)
  {
    return std::ranges::all_of(points,
                               [&](const Point &p)
                               {
                                 return std::ranges::any_of(
                                     x,
                                     [&](const Point &q)
                                     {
                                       return std::abs(p[0] - q[0]) +
                                                  std::abs(p[1] - q[1]) +
                                                  std::abs(p[2] - q[2]) <
                                              1.0e-14;
                                     });
                               });
  };
  IoData iodata(Units(1.0, 1.0));
  iodata.model.mesh = path.string();
  const auto mesh = mesh::ReadInputMesh(iodata);
  REQUIRE(mesh->GetNodes()->FESpace()->GetNDofs() == static_cast<int>(x.size()));
  CHECK(Same(NodePoints(*mesh)));
  CHECK(!Same(NodePoints(mfem::Mesh(path.string()))));
}

TEST_CASE_METHOD(palace::test::SharedTempDir, "Design velocities",
                 "[sensitivity][Serial][Parallel]")
{
  // The velocities at the mesh nodes are the derivatives of the node positions for input
  // meshes morphed by the table, through the mesh pipeline with uniform and nonconformal
  // refinement of curved meshes, or without their curvature. So are the derivatives of the
  // geometry.
#if !defined(MFEM_USE_GSLIB)
  SKIP("Design velocities require GSLIB");
#endif
  auto type = GENERATE(mfem::Element::HEXAHEDRON, mfem::Element::TETRAHEDRON);
  auto amr = GENERATE(false, true);
  auto remove_curvature = GENERATE(false, true);
  INFO("Element type: " << type << ", AMR: " << amr
                        << ", without curvature: " << remove_curvature);
  auto smesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, type, 1.0, 1.0, 1.0);
  for (int e = 0; e < smesh.GetNE(); e++)
  {
    mfem::Vector center;
    smesh.GetElementCenter(e, center);
    smesh.SetAttribute(e, (center(2) < 0.5) ? 1 : 2);
  }
  smesh.SetAttributes();
  smesh.SetCurvature(2);
  auto &nodes = *smesh.GetNodes();
  const auto &fespace = *nodes.FESpace();
  auto Position = [&](int i)
  {
    const auto *h_nodes = nodes.HostRead();
    std::array<double, 3> x;
    for (int d = 0; d < 3; d++)
    {
      x[d] = h_nodes[fespace.DofToVDof(i, d)];
    }
    return x;
  };
  for (int i = 0; i < fespace.GetNDofs(); i++)
  {
    const auto x = Position(i);
    for (int d = 0; d < 3; d++)
    {
      nodes.HostReadWrite()[fespace.DofToVDof(i, d)] +=
          0.12 * x[d] * (1.0 - x[d]) * std::sin(2.0 * x[(d + 1) % 3] + x[(d + 2) % 3] + d);
    }
  }

  // A smooth velocity and a translation along x at the nodes with x < 0.6, the others
  // fixed.
  mfem::Vector W(nodes.Size()), T(nodes.Size());
  W = 0.0;
  T = 0.0;
  std::ofstream table(temp_dir / "velocities.csv");
  table.precision(17);
  table << "x, y, z, w.x, w.y, w.z, t.x, t.y, t.z\n";
  for (int i = 0; i < fespace.GetNDofs(); i++)
  {
    const auto x = Position(i);
    if (x[0] >= 0.6)
    {
      continue;
    }
    T(fespace.DofToVDof(i, 0)) = 1.0;
    const std::array<double, 3> w = {std::sin(2.0 * x[1] + x[2]),
                                     x[1] * (1.0 - x[1]) * std::cos(x[0] + 2.0 * x[2]),
                                     (1.0 - x[2]) * std::sin(x[0] - x[1])};
    for (int d = 0; d < 3; d++)
    {
      W(fespace.DofToVDof(i, d)) = w[d];
    }
    table << x[0] << ", " << x[1] << ", " << x[2] << ", " << w[0] << ", " << w[1] << ", "
          << w[2] << ", 1, 0, 0\n";
  }
  table.close();

  // The input mesh moved by t times a velocity, through the mesh pipeline.
  const mfem::Vector X(nodes);
  int count = 0;
  auto Read = [&](const mfem::Vector *V, double t)
  {
    nodes = X;
    if (V)
    {
      nodes.Add(t, *V);
    }
    const auto path = temp_dir / fmt::format("cube-{:d}.mesh", count++);
    if (Mpi::Root(Mpi::World()))
    {
      std::ofstream fo(path);
      fo.precision(17);
      smesh.Print(fo);
    }
    Mpi::Barrier(Mpi::World());
    nlohmann::json config = {
        {"Problem", {{"Type", "Electrostatic"}, {"Output", (temp_dir / "out").string()}}},
        {"Model",
         {{"Mesh", path.string()},
          {"L0", 1.0e-3},
          {"Lc", 2.0},
          {"Refinement", {{"UniformLevels", 1}}}}},
        {"Domains", {{"Materials", {{{"Attributes", {1}}}, {{"Attributes", {2}}}}}}},
        {"Boundaries",
         {{"Ground", {{"Attributes", {2, 4}}}},
          {"Terminal", {{{"Index", 1}, {"Attributes", {1}}}}}}},
        {"Solver", {{"Order", 1}}},
        {"Sensitivity", {{"Velocities", (temp_dir / "velocities.csv").string()}}}};
    if (amr)
    {
      config["Model"]["Refinement"]["Nonconformal"] = true;
      config["Model"]["Refinement"]["MaxIts"] = 1;
    }
    if (remove_curvature)
    {
      config["Model"]["RemoveCurvature"] = true;
    }
    auto iodata = std::make_unique<IoData>(config, false);
    auto smesh = mesh::Load(*iodata, Mpi::World());
    iodata->NondimensionalizeInputs(smesh);
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::Partition(*iodata, std::move(smesh), Mpi::World()));
    mesh::RefineMesh(*iodata, mfem_mesh);
    if (amr)
    {
      auto &pmesh = *mfem_mesh.back();
      mfem::Array<int> marked;
      for (int e = 0; e < pmesh.GetNE(); e += 3)
      {
        marked.Append(e);
      }
      pmesh.GeneralRefinement(marked, 1, 0);
    }
    return std::pair{std::move(iodata),
                     std::make_unique<Mesh>(std::move(mfem_mesh.back()))};
  };
  auto Nodes = [](const Mesh &mesh)
  {
    const auto &nodes = *mesh.Get().GetNodes();
    const auto &nodal_fespace =
        dynamic_cast<const mfem::ParFiniteElementSpace &>(*nodes.FESpace());
    Vector X(nodal_fespace.GetTrueVSize());
    nodal_fespace.GetRestrictionMatrix()->Mult(nodes, X);
    return X;
  };
  auto Volume = [](const mfem::ParMesh &mesh)
  {
    double volume = 0.0;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      volume += const_cast<mfem::ParMesh &>(mesh).GetElementVolume(e);
    }
    Mpi::GlobalSum(1, &volume, mesh.GetComm());
    return std::vector<double>{volume};
  };

  const auto [iodata, mesh] = Read(nullptr, 0.0);
  const InputVelocities input(*iodata, Mpi::World());
  const DesignVelocities velocities(&input, *iodata, *mesh);
  REQUIRE(velocities.GetNames() == std::vector<std::string>{"w", "t"});
  const auto dvolume = velocities.Differentiate(Volume);
  const double h = 1.0e-3;
  for (const auto &[p, V] : {std::pair{0, &W}, std::pair{1, &T}})
  {
    const auto [iodata_p, mesh_p] = Read(V, h);
    const auto [iodata_m, mesh_m] = Read(V, -h);
    auto dX = Nodes(*mesh_p);
    dX -= Nodes(*mesh_m);
    dX *= 1.0 / (2.0 * h);
    const auto &v = velocities.GetVelocity(p);
    REQUIRE(v.Size() == dX.Size());
    Vector error(v);
    error -= dX;
    INFO("Parameter " << velocities.GetNames()[p]);
    CHECK(linalg::Norml2(Mpi::World(), error) <= 1.0e-9 * linalg::Norml2(Mpi::World(), dX));

    // The volume is cubic in the positions of the nodes: Richardson extrapolation of the
    // central differences is exact.
    const auto [iodata_p2, mesh_p2] = Read(V, 0.5 * h);
    const auto [iodata_m2, mesh_m2] = Read(V, -0.5 * h);
    const double d1 = (Volume(mesh_p->Get())[0] - Volume(mesh_m->Get())[0]) / (2.0 * h),
                 d2 = (Volume(mesh_p2->Get())[0] - Volume(mesh_m2->Get())[0]) / h;
    const double fd = (4.0 * d2 - d1) / 3.0;
    INFO("Volume derivative " << dvolume[p][0] << ", finite difference " << fd);
    CHECK(std::abs(dvolume[p][0] - fd) <= 1.0e-8 * Volume(mesh->Get())[0]);
  }
}

}  // namespace palace
