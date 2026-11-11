// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/operator.hpp"
#include "fem/mesh.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hypre.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "linalg/sliverpatch.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/memoryreporting.hpp"
#include "utils/perfstats.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace Catch::Matchers;
using json = nlohmann::json;

namespace
{

// An octahedron around a sliver tetrahedron ABCD, whose vertices are the corners of the
// unit square with two opposite corners raised by h, and four tetrahedra between its faces
// and the apices E above and F below, refined uniformly. The boundary attribute is 1 below
// and 2 above.
mfem::Mesh MakeSliverMesh(double h, int refinements)
{
  mfem::Mesh mesh(3, 6, 5, 8);
  const double v[6][3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, h},   {1.0, 1.0, 0.0},
                          {0.0, 1.0, h},   {0.5, 0.5, 1.0}, {0.5, 0.5, -1.0}};
  for (const auto &x : v)
  {
    mesh.AddVertex(x);
  }
  enum
  {
    A,
    B,
    C,
    D,
    E,
    F
  };
  const int tets[5][4] = {
      {A, B, C, D}, {A, B, D, E}, {C, D, B, E}, {A, C, B, F}, {A, D, C, F}};
  for (const auto &t : tets)
  {
    mesh.AddTet(t, 1);
  }
  const int bdr[8][3] = {{E, A, B}, {E, B, C}, {E, C, D}, {E, D, A},
                         {F, B, A}, {F, C, B}, {F, D, C}, {F, A, D}};
  for (int i = 0; i < 8; i++)
  {
    mesh.AddBdrTriangle(bdr[i], (i < 4) ? 2 : 1);
  }
  mesh.FinalizeTetMesh(1, 1, true);
  for (int l = 0; l < refinements; l++)
  {
    mesh.UniformRefinement();
  }
  return mesh;
}

void FillRandom(Vector &v, unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  auto *h = v.HostWrite();
  for (int i = 0; i < v.Size(); i++)
  {
    h[i] = dist(gen);
  }
}

constexpr double Degrees(double radians)
{
  return radians * 180.0 / std::numbers::pi;
}

}  // namespace

TEST_CASE("Minimum dihedral angles of tetrahedra", "[sliverpatch][Serial]")
{
  // The corner of a cube: three right angles and three of acos(1/√3).
  const double o[3] = {0.0, 0.0, 0.0}, x[3] = {1.0, 0.0, 0.0}, y[3] = {0.0, 1.0, 0.0},
               z[3] = {0.0, 0.0, 1.0};
  CHECK_THAT(MinDihedralAngle(o, x, y, z),
             WithinRel(Degrees(std::acos(1.0 / std::sqrt(3.0))), 1.0e-12));

  // A regular tetrahedron, for both orientations.
  const double a[3] = {1.0, 1.0, 1.0}, b[3] = {1.0, -1.0, -1.0}, c[3] = {-1.0, 1.0, -1.0},
               d[3] = {-1.0, -1.0, 1.0};
  CHECK_THAT(MinDihedralAngle(a, b, c, d),
             WithinRel(Degrees(std::acos(1.0 / 3.0)), 1.0e-12));
  CHECK_THAT(MinDihedralAngle(b, a, c, d),
             WithinRel(Degrees(std::acos(1.0 / 3.0)), 1.0e-12));

  // A sliver: the corners of the unit square with two opposite ones raised by h, whose
  // dihedral angles at the sides of the square are about 2 h.
  for (double h : {1.0e-1, 1.0e-2, 1.0e-4})
  {
    const double s0[3] = {0.0, 0.0, 0.0}, s1[3] = {1.0, 0.0, h}, s2[3] = {1.0, 1.0, 0.0},
                 s3[3] = {0.0, 1.0, h};
    CHECK_THAT(MinDihedralAngle(s0, s1, s2, s3), WithinRel(Degrees(2.0 * h), 0.02));
  }
}

TEST_CASE("Sliver detection", "[sliverpatch][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  {
    auto smesh = MakeSliverMesh(0.01, 1);
    Mesh mesh(comm, smesh);
    SliverStatistics stats;
    const auto slivers = FindSliverElements(mesh.Get(), 2.0, &stats);
    long long num_local = static_cast<long long>(slivers.size());
    Mpi::GlobalSum(1, &num_local, comm);
    CHECK(stats.num_tets == 40);
    CHECK(stats.num_slivers > 0);
    CHECK(stats.num_slivers == num_local);
    CHECK(stats.num_below_2 == stats.num_slivers);
    CHECK(stats.num_below_1 <= stats.num_below_2);
    CHECK(stats.num_below_2 <= stats.num_below_5);
    CHECK_THAT(stats.min_angle, WithinRel(Degrees(0.02), 0.05));
  }
  {
    // The tetrahedra of a Cartesian mesh have dihedral angles of 45 degrees or more.
    auto smesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
    Mesh mesh(comm, smesh);
    SliverStatistics stats;
    const auto slivers = FindSliverElements(mesh.Get(), 2.0, &stats);
    CHECK(slivers.empty());
    CHECK(stats.num_tets == 48);
    CHECK(stats.num_below_5 == 0);
    CHECK(stats.min_angle >= 45.0 - 1.0e-9);
  }
}

// The entries between marked L-dofs, summed over the elements, are those of the assembled
// local operator (all elements with entries between marked L-dofs have one marked 2).
TEST_CASE("Marked entries of libCEED operators", "[sliverpatch][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const int order = GENERATE(1, 2, 3);
  CAPTURE(order);
  auto smesh = MakeSliverMesh(0.05, 1);
  Mesh mesh(comm, smesh);
  const int p_trial = fem::DefaultIntegrationOrder::p_trial;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection fec(order, 3);
  FiniteElementSpace fespace(mesh, &fec);
  MaterialPropertyCoefficient muinv(mesh.MaxCeedAttribute()), eps(mesh.MaxCeedAttribute());
  muinv.AddMaterialProperty(mesh.GetCeedAttributes(1), 1.5);
  eps.AddMaterialProperty(mesh.GetCeedAttributes(1), -2.0);
  BilinearForm a(fespace);
  a.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv, eps);
  a.AddBoundaryIntegrator<VectorFEMassIntegrator>();
  const auto op = a.PartialAssemble();
  const auto csr = BilinearForm::FullAssemble(*op, false);

  const int n = fespace.GetVSize();
  std::vector<char> marker(n, 0);
  std::mt19937 gen(11 + Mpi::Rank(comm));
  std::uniform_int_distribution<int> dist(0, 4);
  for (int l = 0; l < n; l++)
  {
    marker[l] = (dist(gen) == 0) ? 2 : 0;
  }
  std::map<std::pair<int, int>, double> sum;
  REQUIRE(ceed::AssembleMarkedEntries(*op, marker,
                                      [&](int i, int j, double v) { sum[{i, j}] += v; }));
  hypre_CSRMatrix *mat = *csr;
  const HYPRE_Int *I = hypre_CSRMatrixI(mat), *J = hypre_CSRMatrixJ(mat);
  const double *data = hypre_CSRMatrixData(mat);
  double max_entry = 0.0, err = 0.0;
  std::size_t num_marked = 0;
  for (int i = 0; i < n; i++)
  {
    for (HYPRE_Int k = I[i]; k < I[i + 1]; k++)
    {
      max_entry = std::max(max_entry, std::abs(data[k]));
      if (marker[i] && marker[J[k]])
      {
        auto it = sum.find({i, static_cast<int>(J[k])});
        err = std::max(err, std::abs(data[k] - ((it != sum.end()) ? it->second : 0.0)));
        num_marked += (it != sum.end());
      }
    }
  }
  CHECK(num_marked == sum.size());
  CHECK(err <= 1.0e-12 * max_entry);
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

// The patch solve is exact on its block: for y = Σₚ Rₚᵀ Aₚₚ⁻¹ Rₚ r, A y = r on the patch
// true dofs, and y vanishes elsewhere.
TEST_CASE("Sliver patch solves", "[sliverpatch][Serial][Parallel]")
{
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK | mfem::Backend::OMP_MASK))
  {
    SKIP("Sliver patches are disabled with device and OpenMP backends");
  }
  MPI_Comm comm = Mpi::World();
  const int order = GENERATE(2, 3);
  const int layers = GENERATE(0, 1);
  const bool complex = GENERATE(false, true);
  const bool nonconforming = GENERATE(false, true);
  CAPTURE(order, layers, complex, nonconforming);
  auto smesh = MakeSliverMesh(0.01, 2);
  if (nonconforming)
  {
    // Refine some elements around the slivers, so that patches have constrained L-dofs.
    smesh.EnsureNCMesh(true);
    mfem::Array<int> refine;
    for (int e = 0; e < smesh.GetNE(); e += 3)
    {
      refine.Append(e);
    }
    smesh.GeneralRefinement(refine);
  }
  Mesh mesh(comm, smesh);
  const int p_trial = fem::DefaultIntegrationOrder::p_trial;
  fem::DefaultIntegrationOrder::p_trial = order;
  mfem::ND_FECollection fec(order, 3);
  FiniteElementSpace fespace(mesh, &fec);

  // An indefinite operator (a curl-curl term and a negative mass term), with a mass term
  // as imaginary part.
  MaterialPropertyCoefficient muinv(mesh.MaxCeedAttribute()),
      eps_r(mesh.MaxCeedAttribute()), eps_i(mesh.MaxCeedAttribute());
  muinv.AddMaterialProperty(mesh.GetCeedAttributes(1), 1.0);
  eps_r.AddMaterialProperty(mesh.GetCeedAttributes(1), -3.0);
  eps_i.AddMaterialProperty(mesh.GetCeedAttributes(1), 0.5);
  BilinearForm ar(fespace), ai(fespace);
  ar.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv, eps_r);
  ai.AddDomainIntegrator<VectorFEMassIntegrator>(eps_i);
  mfem::Array<int> bdr_marker(mesh.Get().bdr_attributes.Max()), dbc;
  bdr_marker = 0;
  bdr_marker[0] = 1;
  fespace.Get().GetEssentialTrueDofs(bdr_marker, dbc);
  ComplexParOperator A(ar.PartialAssemble(),
                       complex ? std::unique_ptr<Operator>(ai.PartialAssemble()) : nullptr,
                       fespace);
  A.SetEssentialTrueDofs(dbc, Operator::DiagonalPolicy::DIAG_ONE);

  const auto slivers = FindSliverElements(mesh.Get(), 2.0);
  SliverPatchSolver<ComplexOperator> patch(fespace, slivers, layers, &dbc, 2000);
  long long num[3] = {patch.NumPatches(), patch.NumTrueDofs(), patch.NumRemoteTrueDofs()};
  Mpi::GlobalSum(3, num, comm);
  REQUIRE(num[0] > 0);
  REQUIRE(num[1] > 0);
  const bool ok = patch.SetOperator(A);
  CHECK(ok);

  // One patch at a time (patches can be coupled through the operator), with the residual
  // on the true dofs of the patch owned by each process.
  const int n = fespace.GetTrueVSize();
  const long long my_start = fespace.Get().GetMyTDofOffset();
  std::vector<char> ess(n, 0);
  for (int t : dbc)
  {
    ess[t] = 1;
  }
  ComplexVector r(n), y(n), Ay(n), r0(n);
  FillRandom(r0.Real(), 5 + Mpi::Rank(comm));
  FillRandom(r0.Imag(), 17 + Mpi::Rank(comm));
  std::vector<int> num_patches(Mpi::Size(comm));
  const int my_num = patch.NumPatches();
  Mpi::Allgather(1, &my_num, num_patches.data(), comm);
  const auto &ptr = patch.GetPatchOffsets();
  const auto &gtdofs = patch.GetPatchGlobalTrueDofs();
  for (int q = 0; q < Mpi::Size(comm); q++)
  {
    for (int p = 0; p < num_patches[q]; p++)
    {
      std::vector<long long> g;
      int size = 0;
      if (Mpi::Rank(comm) == q)
      {
        g.assign(gtdofs.begin() + ptr[p], gtdofs.begin() + ptr[p + 1]);
        size = static_cast<int>(g.size());
      }
      MPI_Bcast(&size, 1, MPI_INT, q, comm);
      g.resize(size);
      MPI_Bcast(g.data(), size, MPI_LONG_LONG, q, comm);
      std::vector<char> in_p(n, 0);
      r = 0.0;
      for (long long gt : g)
      {
        const long long t = gt - my_start;
        if (t >= 0 && t < n)
        {
          CHECK(!ess[t]);
          in_p[t] = 1;
          r.Real()[t] = r0.Real()[t];
          r.Imag()[t] = r0.Imag()[t];
        }
      }
      y = 0.0;
      patch.AddMult(r, y);
      A.Mult(y, Ay);
      double err = 0.0, nrm = 0.0, outside = 0.0;
      for (int t = 0; t < n; t++)
      {
        const std::complex<double> res(Ay.Real()[t] - r.Real()[t],
                                       Ay.Imag()[t] - r.Imag()[t]);
        if (in_p[t])
        {
          err = std::max(err, std::abs(res));
          nrm = std::max(nrm, std::abs(std::complex<double>(r.Real()[t], r.Imag()[t])));
        }
        else
        {
          outside =
              std::max(outside, std::abs(std::complex<double>(y.Real()[t], y.Imag()[t])));
        }
      }
      Mpi::GlobalMax(1, &err, comm);
      Mpi::GlobalMax(1, &nrm, comm);
      Mpi::GlobalMax(1, &outside, comm);
      CAPTURE(q, p, size);
      CHECK(nrm > 0.0);
      CHECK(err <= 1.0e-8 * nrm);
      CHECK(outside == 0.0);
    }
  }
  r = r0;
  y = 0.0;
  patch.AddMult(r, y);

  // With residual updates, the updated residual is r - A δ on all true dofs.
  {
    SliverPatchSolver<ComplexOperator> patch_u(fespace, slivers, layers, &dbc, 2000, true);
    CHECK(patch_u.SetOperator(A));
    ComplexVector ru(n), yu(n), ref(n);
    ru = r0;
    yu = 0.0;
    patch_u.AddMultUpdate(ru, yu);
    A.Mult(yu, ref);
    linalg::AXPBY(1.0, r0, -1.0, ref);
    linalg::AXPBY(-1.0, y, 1.0, yu);
    linalg::AXPBY(-1.0, ref, 1.0, ru);
    double err_y = linalg::Norml2(comm, yu), err_r = linalg::Norml2(comm, ru);
    const double nrm_r = linalg::Norml2(comm, r0), nrm_y = linalg::Norml2(comm, y);
    CHECK(err_y <= 1.0e-10 * nrm_y);
    CHECK(err_r <= 1.0e-10 * nrm_r);
  }

  // Real operators use the same blocks.
  if (!complex)
  {
    ParOperator Ar(ar.PartialAssemble(), fespace);
    Ar.SetEssentialTrueDofs(dbc, Operator::DiagonalPolicy::DIAG_ONE);
    SliverPatchSolver<Operator> patch_r(fespace, slivers, layers, &dbc, 2000);
    CHECK(patch_r.SetOperator(Ar));
    Vector yr(n);
    yr = 0.0;
    patch_r.AddMult(r.Real(), yr);
    yr -= y.Real();
    double diff = yr.Normlinf();
    Mpi::GlobalMax(1, &diff, comm);
    CHECK(diff <= 1.0e-12 * std::max(1.0, y.Real().Normlinf()));
  }
  fem::DefaultIntegrationOrder::p_trial = p_trial;
}

namespace
{

// Options of SolveShiftedMaxwell: sliver patches, the Chebyshev smoothing order (the
// eigenmode default unless positive), the loss tangent, the frequency, the absorbing
// boundary attributes (the others are PEC), and the coarse solve: the real-valued
// approximation by default, the complex-valued one with complex_coarse true, and chosen
// at setup without complex_coarse, with the sparse direct solver coarse_type if given;
// configure is applied to the multigrid preconditioner before its setup.
struct ShiftedMaxwellOptions
{
  double sliver_angle = 0.0;
  int sliver_layers = 0, smooth_order = -1;
  double loss_tan = 0.01, omega = 2.0;
  std::vector<int> absorbing = {2};
  std::optional<bool> complex_coarse = false;
  const char *coarse_type = nullptr;
  std::function<void(GeometricMultigridSolver<ComplexOperator> &)> configure;
  json patch = json::object();
};

// Krylov iterations of SolveShiftedMaxwell, whether the multigrid solver used and tried the
// complex-valued coarse solve chosen at setup, and its choice.
struct ShiftedMaxwellResult
{
  int its = 0;
  bool used = false, tried = false;
  std::optional<CoarseSolveChoice> choice;
};

// Solve a shifted Maxwell system with the multigrid preconditioner of a driven or
// eigenmode simulation (p-multigrid with distributive relaxation smoothers).
ShiftedMaxwellResult SolveShiftedMaxwell(const mfem::Mesh &smesh, ComplexVector &x,
                                         const ShiftedMaxwellOptions &opt)
{
  MPI_Comm comm = Mpi::World();
  std::vector<int> pec;
  for (int attr = 1; attr <= smesh.bdr_attributes.Max(); attr++)
  {
    if (std::find(opt.absorbing.begin(), opt.absorbing.end(), attr) == opt.absorbing.end())
    {
      pec.push_back(attr);
    }
  }
  json config = {
      {"Problem", {{"Type", "Eigenmode"}, {"Verbose", 0}}},
      {"Model", {{"Mesh", "test.mesh"}, {"L0", 1.0}}},
      {"Domains",
       {{"Materials",
         {{{"Attributes", {1}}, {"Permittivity", 2.0}, {"LossTan", opt.loss_tan}}}}}},
      {"Boundaries",
       {{"PEC", {{"Attributes", pec}}},
        {"Absorbing", {{"Attributes", opt.absorbing}, {"Order", 1}}}}},
      {"Solver",
       {{"Order", 2},
        {"Eigenmode", {{"Target", 0.2}}},
        {"Linear",
         {{"Tol", 1.0e-10},
          {"MaxIts", 400},
          {"MGSliverAngle", opt.sliver_angle},
          {"MGSliverLayers", opt.sliver_layers}}}}}};
  if (opt.smooth_order > 0)
  {
    config["Solver"]["Linear"]["MGSmoothOrder"] = opt.smooth_order;
  }
  if (opt.complex_coarse)
  {
    config["Solver"]["Linear"]["ComplexCoarseSolve"] = *opt.complex_coarse;
  }
  if (opt.coarse_type)
  {
    config["Solver"]["Linear"]["Type"] = opt.coarse_type;
  }
  config.merge_patch(opt.patch);
  IoData iodata(config, false);
  iodata.model.Lc = 1.0;
  auto serial_mesh = std::make_unique<mfem::Mesh>(smesh);
  iodata.NondimensionalizeInputs(serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(comm, *serial_mesh));
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  using namespace std::complex_literals;
  const double omega = opt.omega;
  auto A = space_op.GetSystemMatrix(1.0 + 0.0i, 1i * omega, -omega * omega + 0.0i, K.get(),
                                    C.get(), M.get());
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0 + 0.0i, 1i * omega,
                                                             -omega * omega + 0.0i, omega);
  ComplexKspSolver ksp(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
  auto *gmg =
      dynamic_cast<GeometricMultigridSolver<ComplexOperator> *>(ksp.GetPreconditioner());
  if (opt.configure)
  {
    REQUIRE(gmg);
    opt.configure(*gmg);
  }
  ksp.SetOperators(*A, *P);
  const int n = A->Height();
  ComplexVector b(n);
  FillRandom(b.Real(), 3 + Mpi::Rank(comm));
  FillRandom(b.Imag(), 23 + Mpi::Rank(comm));
  linalg::SetSubVector(b, space_op.GetNDDbcTDofLists().back(), 0.0);
  x.SetSize(n);
  x = 0.0;
  ksp.Mult(b, x);
  ShiftedMaxwellResult result;
  result.its = ksp.NumTotalMultIterations();
  if (gmg)
  {
    result.used = gmg->UsesComplexCoarseSolver();
    result.tried = gmg->TriedComplexCoarseSolver();
    result.choice = gmg->GetCoarseSolveChoice();
  }
  return result;
}

int SolveShiftedMaxwell(const mfem::Mesh &smesh, double sliver_angle, int sliver_layers,
                        ComplexVector &x, int smooth_order = -1,
                        const json &patch = json::object())
{
  ShiftedMaxwellOptions opt;
  opt.sliver_angle = sliver_angle;
  opt.sliver_layers = sliver_layers;
  opt.smooth_order = smooth_order;
  opt.patch = patch;
  return SolveShiftedMaxwell(smesh, x, opt).its;
}

// The solution x1 with patches enabled, which create no patch, is that of the solve x0
// without them, to a tolerance far below that of the solves: the parallel factorizations of
// MUMPS, the coarse solver of builds without STRUMPACK or SuperLU, are not reproducible bit
// for bit.
void CheckSameSolution(MPI_Comm comm, const ComplexVector &x0, ComplexVector &x1,
                       double tol = 1.0e-10)
{
  x1 -= x0;
  CHECK(linalg::Norml2(comm, x1) <= tol * linalg::Norml2(comm, x0));
}

}  // namespace

TEST_CASE("Raised Chebyshev smoothing orders are checked by a GMRES probe",
          "[sliverpatch][Serial][Parallel]")
{
  // The order-2 eigenmode default raises the Chebyshev order from 4 to 8 on the nearly
  // positive semi-definite levels, and the multigrid preconditioner keeps the raised or the
  // base orders after a GMRES probe. Either way, the preconditioner is then the one with
  // that order given explicitly, so the solve takes the same iterations, and it takes no
  // more than the slower of the two. The coarse solve is the real-valued approximation
  // ("ComplexCoarseSolve": false), whose cycle is where order 4 can be kept.
  const auto smesh = MakeSliverMesh(0.002, 2);
  ComplexVector x, x4, x8;
  const int its = SolveShiftedMaxwell(smesh, 0.0, 0, x);
  const int its4 = SolveShiftedMaxwell(smesh, 0.0, 0, x4, 4);
  const int its8 = SolveShiftedMaxwell(smesh, 0.0, 0, x8, 8);
  INFO("iterations: default " << its << ", order 4 " << its4 << ", order 8 " << its8);
  CHECK((its == its4 || its == its8));
  CHECK(its <= std::max(its4, its8));
  CheckSameSolution(Mpi::World(), (its == its8) ? x8 : x4, x);
}

TEST_CASE("Raised Chebyshev smoothing orders fall back to the base orders",
          "[sliverpatch][Serial][Parallel]")
{
  // With Floquet terms and a poor coarse-grid correction (a Jacobi sweep on the coarse
  // level), the order raised on the nearly positive semi-definite level slows the solve
  // down: the GMRES probe falls back to the base order, and the solve is then the one with
  // order 4 given explicitly.
  auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
  for (int i = 0; i < smesh.GetNBE(); i++)
  {
    smesh.SetBdrAttribute(i, (smesh.GetBdrAttribute(i) == 1) ? 1 : 2);
  }
  smesh.SetAttributes();
  const json patch = {
      {"Boundaries",
       {{"Periodic",
         {{"FloquetWaveVector", {2.0, 0.0, 0.0}}, {"BoundaryPairs", json::array()}}}}},
      {"Solver", {{"Linear", {{"Type", "Jacobi"}}}}}};
  MPI_Comm comm = Mpi::World();
  ComplexVector x, x4, x8;
  BlockTimer::Reset();
  const int its = SolveShiftedMaxwell(smesh, 0.0, 0, x, -1, patch);
  BlockTimer::Finalize(comm);
  json meta = json::object();
  perf::WriteJson(meta);
  BlockTimer::Reset();

  // The order was raised (the probe ran), the probe chose the base order, and the levels
  // use it.
  const auto &ls = meta.at("LinearSolver").at("Main");
  REQUIRE(ls.contains("ChebyshevOrderProbe"));
  const auto &probe = ls.at("ChebyshevOrderProbe");
  CHECK(!probe.at("Raised").get<bool>());
  CHECK(probe.at("PredictedIterationRatio").get<double>() > 1.25);
  for (const auto &level : ls.at("Chebyshev"))
  {
    CHECK(level.at("Primary").at("Order").get<int>() == 4);
  }
  const int its4 = SolveShiftedMaxwell(smesh, 0.0, 0, x4, 4, patch);
  const int its8 = SolveShiftedMaxwell(smesh, 0.0, 0, x8, 8, patch);
  INFO("iterations: default " << its << ", order 4 " << its4 << ", order 8 " << its8);
  CHECK(4 * its8 > 5 * its4);
  CHECK(its == its4);
  CheckSameSolution(comm, x4, x);
}

TEST_CASE("Complex-valued coarse solve chosen at setup", "[sliverpatch][Serial][Parallel]")
{
  // Without ComplexCoarseSolve, the multigrid preconditioner of an eigenmode simulation
  // uses the complex-valued coarse solve when GMRES probes find the real-valued
  // approximation Ar + Ai of the coarse operator poor and the cycle faster with the
  // complex-valued one, and otherwise the approximation: the solve is then that with the
  // explicit choice. The available memory is set, so that the check runs where it is
  // unknown, and all sparse direct solvers are allowed to choose (by default only
  // STRUMPACK does, see ChoosesComplexCoarseSolve).
  std::vector<const char *> types;
#if defined(MFEM_USE_STRUMPACK)
  types.push_back("STRUMPACK");
#endif
#if defined(MFEM_USE_SUPERLU)
  types.push_back("SuperLU");
#endif
#if defined(MFEM_USE_MUMPS)
  types.push_back("MUMPS");
#endif
  if (types.empty())
  {
    SKIP("No sparse direct solver");
  }
  const char *type = GENERATE_COPY(from_range(types));
  struct Settings
  {
    Settings()
    {
      memory_reporting::SetAvailableMemoryForTesting(1L << 40);
      SetComplexCoarseSolversForTesting(
          {LinearSolver::STRUMPACK, LinearSolver::SUPERLU, LinearSolver::MUMPS});
    }
    ~Settings()
    {
      memory_reporting::SetAvailableMemoryForTesting(-1);
      SetComplexCoarseSolversForTesting({});
    }
  } settings;
  MPI_Comm comm = Mpi::World();
  ShiftedMaxwellOptions opt;
  opt.coarse_type = type;
  auto Solve = [&](const mfem::Mesh &smesh, ComplexVector &x, std::optional<bool> complex)
  {
    opt.complex_coarse = complex;
    return SolveShiftedMaxwell(smesh, x, opt);
  };
  auto CheckSameSolve = [&](const ShiftedMaxwellResult &r, const ShiftedMaxwellResult &ref,
                            ComplexVector &x, const ComplexVector &x_ref)
  {
    // The parallel factorizations of SuperLU are not reproducible bit for bit.
    INFO(type << ": iterations " << r.its << " (complex coarse solve " << r.used
              << ", tried " << r.tried << "), " << ref.its << " with the explicit choice");
    if (std::string(type) == "SuperLU")
    {
      CHECK(std::abs(r.its - ref.its) <= 2);
    }
    else
    {
      CHECK(r.its == ref.its);
    }
    CheckSameSolution(comm, x_ref, x, 1.0e-8);
  };
  ComplexVector x, x_ref;
  SECTION("Kept when the cycle is limited by its smoothing")
  {
    // On the mesh with a cluster of slivers, the coarse probe finds the approximation
    // poor, but the complex-valued coarse solve does not speed up the cycle.
    const auto smesh = MakeSliverMesh(0.002, 2);
    const auto r = Solve(smesh, x, std::nullopt);
    const auto ref = Solve(smesh, x_ref, false);
    CHECK(r.tried);
    CHECK_FALSE(r.used);
    CheckSameSolve(r, ref, x, x_ref);
  }
  SECTION("Used with absorbing boundaries")
  {
    // With absorbing boundaries on most of the boundary of a regular mesh, the
    // complex-valued coarse solve needs about three times fewer iterations.
    const auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
    opt.omega = 3.0;
    opt.absorbing = {2, 3, 4, 5, 6};
    const auto r = Solve(smesh, x, std::nullopt);
    const auto ref = Solve(smesh, x_ref, true);
    CHECK(r.tried);
    CHECK(r.used);
    CHECK_FALSE(ref.used);  // An explicit choice is used from the start
    CheckSameSolve(r, ref, x, x_ref);
    REQUIRE(r.choice);
    CHECK(r.choice->complex);

    // A later operator of the problem (for example of a repeated eigenvalue solve) uses
    // the choice without the check.
    ComplexVector x2;
    opt.complex_coarse = std::nullopt;
    opt.configure = [&](GeometricMultigridSolver<ComplexOperator> &gmg)
    { gmg.SetCoarseSolveChoice(*r.choice); };
    const auto r2 = SolveShiftedMaxwell(smesh, x2, opt);
    CHECK_FALSE(r2.tried);
    CHECK(r2.used);
    CheckSameSolve(r2, ref, x2, x_ref);
  }
  SECTION("Not tried without enough memory or by default")
  {
    // Too little available memory keeps the approximation, and so do the solvers which do
    // not choose by default.
    const auto smesh = mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON);
    opt.omega = 3.0;
    opt.absorbing = {2, 3, 4, 5, 6};
    const auto ref = Solve(smesh, x_ref, false);
    memory_reporting::SetAvailableMemoryForTesting(0);
    const auto r = Solve(smesh, x, std::nullopt);
    CHECK_FALSE(r.tried);
    CHECK_FALSE(r.used);
    CheckSameSolve(r, ref, x, x_ref);
    memory_reporting::SetAvailableMemoryForTesting(1L << 40);
    SetComplexCoarseSolversForTesting({});
    ComplexVector x2;
    const auto r2 = Solve(smesh, x2, std::nullopt);
    CHECK(r2.tried == (std::string(type) == "STRUMPACK"));
    if (std::string(type) != "STRUMPACK")
    {
      CHECK_FALSE(r2.used);
      CheckSameSolve(r2, ref, x2, x_ref);
    }
  }
}

TEST_CASE("Sliver patch multigrid", "[sliverpatch][Serial][Parallel]")
{
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK | mfem::Backend::OMP_MASK))
  {
    SKIP("Sliver patches are disabled with device and OpenMP backends");
  }
  MPI_Comm comm = Mpi::World();
  SECTION("Fewer iterations on a mesh with slivers")
  {
    // A cluster of slivers (the children of the refined sliver): the patches of its
    // vertex stars, merged, cover it, which the patches of the own entities do not.
    const int layers = GENERATE(0, 1);
    const auto smesh = MakeSliverMesh(0.002, 2);
    ComplexVector x0, x1;
    const int its0 = SolveShiftedMaxwell(smesh, 0.0, layers, x0);
    const int its1 = SolveShiftedMaxwell(smesh, 2.0, layers, x1);
    INFO("layers " << layers << ": iterations without patches " << its0 << ", with patches "
                   << its1);
    CHECK(its1 < ((layers == 1) ? its0 / 2 : its0));
    x1 -= x0;
    const double diff = linalg::Norml2(comm, x1), ref = linalg::Norml2(comm, x0);
    CHECK(diff <= 1.0e-6 * ref);
  }
  SECTION("No change with expensive automatic patches")
  {
    // The dense work of the vertex-star patches of the cluster is large compared with the
    // smoothing of this small level: the automatic extent keeps no patch solves.
    const auto smesh = MakeSliverMesh(0.002, 2);
    ComplexVector x0, x1;
    const int its0 = SolveShiftedMaxwell(smesh, 0.0, -1, x0);
    const int its1 = SolveShiftedMaxwell(smesh, 2.0, -1, x1);
    CHECK(its1 == its0);
    CheckSameSolution(comm, x0, x1);
  }
  SECTION("No change without slivers")
  {
    auto smesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::TETRAHEDRON);
    for (int i = 0; i < smesh.GetNBE(); i++)
    {
      smesh.SetBdrAttribute(i, (smesh.GetBdrAttribute(i) == 1) ? 1 : 2);
    }
    smesh.SetAttributes();
    ComplexVector x0, x1;
    const int its0 = SolveShiftedMaxwell(smesh, 0.0, 1, x0);
    const int its1 = SolveShiftedMaxwell(smesh, 2.0, 1, x1);
    CHECK(its1 == its0);
    CheckSameSolution(comm, x0, x1);
  }
}

}  // namespace palace
