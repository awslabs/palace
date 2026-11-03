// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <random>
#include <set>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/shapederivativeform.hpp"
#include "fixtures.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/laplaceoperator.hpp"
#include "models/sensitivityexport.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

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
  auto &nodes = *smesh->GetNodes();
  const auto &nodal_fespace = *nodes.FESpace();
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    double x[3];
    for (int d = 0; d < 3; d++)
    {
      x[d] = nodes(nodal_fespace.DofToVDof(i, d));
    }
    for (int d = 0; d < 3; d++)
    {
      nodes(nodal_fespace.DofToVDof(i, d)) +=
          0.12 * x[d] * (1.0 - x[d]) * std::sin(2.0 * x[(d + 1) % 3] + x[(d + 2) % 3] + d);
    }
  }
  iodata.NondimensionalizeInputs(smesh);
  iodata.CheckConfiguration();
  if (amr)
  {
    smesh->EnsureNCMesh(true);
  }
  auto pmesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *smesh);
  if (amr)
  {
    pmesh->RandomRefinement(0.5);
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

  EigenProblem(mfem::Element::Type type, int order, bool amr) : iodata(Units(1.0, 1.0))
  {
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
  const int size = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector x(size), y(size);
  for (int i = 0; i < size; i++)
  {
    x.Real()(i) = dist(gen);
    x.Imag()(i) = dist(gen);
  }
  for (const int i : space_op.GetNDDbcTDofLists().back())
  {
    x.Real()(i) = x.Imag()(i) = 0.0;
  }
  auto Form = [&]()
  {
    auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A2 = space_op.GetExtraSystemMatrix(omega, Operator::DIAG_ZERO);
    REQUIRE((K && C && M && A2));
    auto T = space_op.GetSystemMatrix(std::complex<double>(1.0), lambda, lambda * lambda,
                                      K.get(), C.get(), M.get(), A2.get());
    T->Mult(x, y);
    std::complex<double> dot = {
        linalg::LocalDot(x.Real(), y.Real()) - linalg::LocalDot(x.Imag(), y.Imag()),
        linalg::LocalDot(x.Real(), y.Imag()) + linalg::LocalDot(x.Imag(), y.Real())};
    Mpi::GlobalSum(1, &dot, Mpi::World());
    return dot;
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

TEST_CASE_METHOD(palace::test::SharedTempDir, "Sensitivity export",
                 "[sensitivity][Serial][Parallel]")
{
  // The exported gradients of the capacitances satisfy the exact discrete identities: zero
  // for translations, and C for a dilation, since C scales with length. So does the sum of
  // the derivatives with respect to the permittivities.
  auto amr = GENERATE(false, true);
  INFO("AMR: " << amr);
  Problem problem(mfem::Element::HEXAHEDRON, 2, amr);
  problem.iodata.sensitivity.shape.emplace().boundary_tags = true;
  LaplaceOperator laplace_op(problem.iodata, problem.mesh);
  const auto V = SolveTerminals(problem.iodata, laplace_op);
  const auto dK = laplace_op.GetStiffnessShapeDerivative();

  SensitivityExport output(problem.iodata, *problem.mesh[0]);
  REQUIRE(output.HasShape());
  std::vector<double> C;
  Vector g;
  for (int i = 0; i < 2; i++)
  {
    C.push_back(Capacitance(laplace_op, {1, 2}, V[i], V[1]));
    dK->Mult(V[i], V[1], g);
    output.AddObjective({{"Index", i + 1},
                         {"Type", "Capacitance"},
                         {"Terminals", {i + 1, 2}},
                         {"Value", C.back()},
                         {"Unit", "F"}},
                        &g, 1.0);
  }
  const auto dir = temp_dir / "sensitivity";
  output.Write(dir);
  const auto num_nodes = dK->GetNodalSpace().GlobalTrueVSize() / 3;  // Collective
  Mpi::Barrier(Mpi::World());
  if (!Mpi::Root(Mpi::World()))
  {
    return;
  }

  std::ifstream fj(dir / "sensitivity.json");
  const auto meta = nlohmann::json::parse(fj);
  const std::size_t n = meta.at("NumPoints");
  CHECK(n == static_cast<std::size_t>(num_nodes));
  auto ReadBinary = [&](const std::string &file, auto value)
  {
    std::ifstream fi(dir / file, std::ios::binary);
    std::vector<decltype(value)> data(std::filesystem::file_size(dir / file) /
                                      sizeof(value));
    fi.read(reinterpret_cast<char *>(data.data()), data.size() * sizeof(value));
    return data;
  };
  const auto points = ReadBinary(meta.at("Points"), 0.0);
  REQUIRE(points.size() == 3 * n);
  for (std::size_t k = 0; k < 2; k++)
  {
    const auto &objective = meta.at("Objectives")[k];
    const auto gradient = ReadBinary(objective.at("Gradient"), 0.0);
    REQUIRE(gradient.size() == 3 * n);
    double dilation = 0.0, scale = 0.0;
    std::array<double, 3> translation = {0.0, 0.0, 0.0};
    for (std::size_t p = 0; p < 3 * n; p++)
    {
      dilation += gradient[p] * points[p];
      translation[p % 3] += gradient[p];
      scale += std::abs(gradient[p]);
    }
    CHECK(std::abs(dilation - C[k]) <= 1.0e-10 * std::abs(C[k]));
    for (auto t : translation)
    {
      CHECK(std::abs(t) <= 1.0e-12 * scale);
    }
    CHECK(objective.at("Value") == C[k]);
  }

  // The nodes of every face of the cube are tagged with its boundary attribute.
  std::ifstream fb(dir / meta.at("Boundary").get<std::string>(), std::ios::binary);
  std::vector<std::int64_t> offsets(n + 1);
  fb.read(reinterpret_cast<char *>(offsets.data()), offsets.size() * sizeof(std::int64_t));
  REQUIRE(fb.good());
  CHECK(offsets.front() == 0);
  CHECK(std::is_sorted(offsets.begin(), offsets.end()));
  std::vector<std::int32_t> attributes(offsets.back());
  fb.read(reinterpret_cast<char *>(attributes.data()),
          attributes.size() * sizeof(std::int32_t));
  CHECK(fb.good());
  CHECK(fb.peek() == std::char_traits<char>::eof());
  CHECK(std::set<int>(attributes.begin(), attributes.end()) ==
        std::set<int>{1, 2, 3, 4, 5, 6});
}

}  // namespace palace
