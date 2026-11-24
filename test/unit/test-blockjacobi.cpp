// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/aggregateblocks.hpp"
#include "linalg/amg.hpp"
#include "linalg/blockjacobi.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "linalg/relaxation.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"
#include "linalg/vertexaggregates.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"

namespace palace
{

using namespace Catch::Matchers;

namespace
{

// Box [0, 1]² × [0, height] with n × n × nz cells: for a small height, a stack of nz layers
// of flat elements. With refine, the first element is refined once (nonconforming), or
// with refine_half the elements of the half x < 1/2 of the box.
Mesh MakeLayeredMesh(MPI_Comm comm, int n, int nz, double height, mfem::Element::Type type,
                     bool refine = false, bool refine_half = false)
{
  // The forms of these tests have at most order 2.
  fem::DefaultIntegrationOrder::p_trial = 2;
  auto smesh = mfem::Mesh::MakeCartesian3D(n, n, nz, type, 1.0, 1.0, height);
  if (refine)
  {
    smesh.EnsureNCMesh(true);
    mfem::Array<int> elems;
    for (int e = 0; e < smesh.GetNE(); e++)
    {
      mfem::Vector center(3);
      smesh.GetElementCenter(e, center);
      if (e == 0 || (refine_half && center(0) < 0.5))
      {
        elems.Append(e);
      }
    }
    smesh.GeneralRefinement(elems, 1);
  }
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  return Mesh(std::make_unique<mfem::ParMesh>(comm, smesh));
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

// Dense matrix of an operator (serial).
template <typename Apply>
Eigen::MatrixXd ToDense(int n, Apply &&apply)
{
  Eigen::MatrixXd M(n, n);
  Vector x(n), y(n);
  x.UseDevice(true);
  y.UseDevice(true);
  for (int j = 0; j < n; j++)
  {
    x = 0.0;
    x.HostReadWrite()[j] = 1.0;
    apply(x, y);
    const auto *Y = y.HostRead();
    for (int i = 0; i < n; i++)
    {
      M(i, j) = Y[i];
    }
  }
  return M;
}

// Dense matrix of the assembled operator (serial).
Eigen::MatrixXd AssembleDense(const ParOperator &op)
{
  mfem::SparseMatrix diag;
  op.ParallelAssemble().GetDiag(diag);
  const int n = diag.Height();
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  const int *I = diag.HostReadI(), *J = diag.HostReadJ();
  const double *A = diag.HostReadData();
  for (int i = 0; i < n; i++)
  {
    for (int k = I[i]; k < I[i + 1]; k++)
    {
      M(i, J[k]) += A[k];
    }
  }
  return M;
}

// The sum of the stiffness (diffusion or curl-curl) and mass matrices of a space: a
// symmetric positive definite operator.
enum class Space
{
  H1,
  ND
};

// The same coefficient in all elements of the mesh.
MaterialPropertyCoefficient ConstantCoefficient(const Mesh &mesh, double value)
{
  mfem::Array<int> attr_mat(mesh.MaxCeedAttribute());
  attr_mat = 0;
  mfem::DenseTensor mat_coeff(1, 1, 1);
  mat_coeff(0)(0, 0) = value;
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

// The largest block and the Chebyshev order of the smoothers of these tests.
constexpr int BLOCK_SIZE = 64, CHEBY_ORDER = 4;

// Coefficients (stiffness, mass) of operators which are positive definite, negative
// definite, and indefinite on the layered meshes of these tests: for the nodal DOFs over
// each other in layers of thickness 0.002, the stiffness matrix is larger than 1e5 times
// the mass matrix except on the functions which are constant through the layers.
constexpr std::pair<double, double> POSITIVE{1.0, 1.0}, NEGATIVE{-1.0, -1.0},
    INDEFINITE{1.0, -1.0e5};

std::unique_ptr<mfem::FiniteElementCollection> MakeFECollection(Space space, int order)
{
  if (space == Space::H1)
  {
    return std::make_unique<mfem::H1_FECollection>(order, 3);
  }
  return std::make_unique<mfem::ND_FECollection>(order, 3);
}

void AddIntegrators(BilinearForm &a, Space space)
{
  if (space == Space::H1)
  {
    a.AddDomainIntegrator<MassIntegrator>();
    a.AddDomainIntegrator<DiffusionIntegrator>();
  }
  else
  {
    a.AddDomainIntegrator<VectorFEMassIntegrator>();
    a.AddDomainIntegrator<CurlCurlIntegrator>();
  }
}

void AddIntegrators(BilinearForm &a, Space space, MaterialPropertyCoefficient &stiffness,
                    MaterialPropertyCoefficient &mass)
{
  if (space == Space::H1)
  {
    a.AddDomainIntegrator<MassIntegrator>(mass);
    a.AddDomainIntegrator<DiffusionIntegrator>(stiffness);
  }
  else
  {
    a.AddDomainIntegrator<VectorFEMassIntegrator>(mass);
    a.AddDomainIntegrator<CurlCurlIntegrator>(stiffness);
  }
}

// Every third true DOF of the space as an essential DOF.
mfem::Array<int> SomeEssentialTrueDofs(const FiniteElementSpace &fespace)
{
  mfem::Array<int> list;
  for (int i = 0; i < fespace.GetTrueVSize(); i += 3)
  {
    list.Append(i);
  }
  return list;
}

// True DOF vector of the functional v -> ∫ f · v (or ∫ f v for a nodal space), which does
// not depend on the partitioning of the mesh or on the orientation of the basis functions.
Vector LoadVector(FiniteElementSpace &fespace, Space space)
{
  mfem::ParLinearForm lf(&fespace.Get());
  mfem::FunctionCoefficient f([](const mfem::Vector &x)
                              { return 1.0 + x(0) - 2.0 * x(1) + 30.0 * x(2) * x(0); });
  mfem::VectorFunctionCoefficient vf(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 1.0 + x(1);
                                       v(1) = 2.0 - x(0) * x(1);
                                       v(2) = 0.5 + 20.0 * x(2) + x(0);
                                     });
  if (space == Space::H1)
  {
    lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(f));
  }
  else
  {
    lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(vf));
  }
  lf.UseFastAssembly(false);
  lf.Assemble();
  Vector b(fespace.GetTrueVSize());
  b.UseDevice(true);
  fespace.GetProlongationMatrix()->MultTranspose(lf, b);
  return b;
}

// The quadratic form bᵀ M⁻¹ b of the block Jacobi smoother for a problem on the given
// communicator.
double SmootherQuadraticForm(MPI_Comm comm, Space space, int order,
                             mfem::Element::Type type, bool refine, bool essential,
                             HYPRE_BigInt *gathered, HYPRE_BigInt *split)
{
  auto mesh = MakeLayeredMesh(comm, 3, 2, 0.01, type, refine);
  const auto fec = MakeFECollection(space, order);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, space);
  ParOperator A(a.PartialAssemble(), fespace);
  mfem::Array<int> ess_list;
  if (essential)
  {
    // Essential DOFs which do not depend on the partitioning: those of a boundary.
    mfem::Array<int> marker(mesh.Get().bdr_attributes.Max());
    marker = 0;
    marker[0] = 1;
    fespace.Get().GetEssentialTrueDofs(marker, ess_list);
    A.SetEssentialTrueDofs(ess_list, Operator::DIAG_ONE);
  }
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  REQUIRE(smoother.HasBlocks());
  *gathered = smoother.GetBlocks()->GetCensus().gathered;
  *split = smoother.GetBlocks()->GetCensus().split;
  Vector b = LoadVector(fespace, space), y(b.Size());
  y.UseDevice(true);
  linalg::SetSubVector(b, ess_list, 0.0);
  smoother.Mult(b, y);
  return linalg::Dot(comm, b, y);
}

// Connected components of the sparsity pattern of a symmetric matrix.
std::vector<std::vector<int>> Components(const Eigen::MatrixXd &M, double tol)
{
  const int n = static_cast<int>(M.rows());
  std::vector<int> comp(n, -1);
  std::vector<std::vector<int>> out;
  for (int s = 0; s < n; s++)
  {
    if (comp[s] >= 0)
    {
      continue;
    }
    std::vector<int> stack{s}, members;
    comp[s] = static_cast<int>(out.size());
    while (!stack.empty())
    {
      const int i = stack.back();
      stack.pop_back();
      members.push_back(i);
      for (int j = 0; j < n; j++)
      {
        if (comp[j] < 0 && std::abs(M(i, j)) > tol)
        {
          comp[j] = comp[s];
          stack.push_back(j);
        }
      }
    }
    std::sort(members.begin(), members.end());
    out.push_back(std::move(members));
  }
  return out;
}

// Symmetric tridiagonal test operator with a varying diagonal.
class TridiagonalOperator : public Operator
{
public:
  TridiagonalOperator(int n) : Operator(n) {}

  double Diagonal(int i) const { return 2.0 + 3.0 * i / height; }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(height);
    for (int i = 0; i < height; i++)
    {
      diag(i) = Diagonal(i);
    }
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    const auto *X = x.HostRead();
    auto *Y = y.HostWrite();
    for (int i = 0; i < height; i++)
    {
      Y[i] = Diagonal(i) * X[i] - ((i > 0) ? X[i - 1] : 0.0) -
             ((i < height - 1) ? X[i + 1] : 0.0);
    }
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

// A relaxation with a sweep correction: point Jacobi on all DOFs but the last ones, which
// are corrected once per sweep with a multiple of their residual.
class CorrectedRelaxation : public JacobiSmoother<Operator>
{
private:
  const int num_corrected;
  const double scale;

  void ZeroCorrected(Vector &d) const
  {
    auto *D = d.HostReadWrite();
    for (int i = d.Size() - num_corrected; i < d.Size(); i++)
    {
      D[i] = 0.0;
    }
  }

public:
  CorrectedRelaxation(MPI_Comm comm, int num_corrected, double scale)
    : JacobiSmoother<Operator>(comm), num_corrected(num_corrected), scale(scale)
  {
  }

  void ScaledMult(double s, const Vector &r, Vector &d) const override
  {
    JacobiSmoother<Operator>::ScaledMult(s, r, d);
    ZeroCorrected(d);
  }

  void ScaledAddMult(double sd, double sr, const Vector &r, Vector &d) const override
  {
    JacobiSmoother<Operator>::ScaledAddMult(sd, sr, r, d);
    ZeroCorrected(d);
  }

  bool HasSweepCorrection() const override { return num_corrected > 0; }

  void AddSweepCorrection(const Vector &r, Vector &y) const override
  {
    const auto *R = r.HostRead();
    auto *Y = y.HostReadWrite();
    for (int i = y.Size() - num_corrected; i < y.Size(); i++)
    {
      Y[i] += scale * R[i];
    }
  }
};

}  // namespace

TEST_CASE("Vertex aggregates follow the flat direction of the elements",
          "[blockjacobi][Serial][Parallel]")
{
  const auto comm = Mpi::World();
  const auto type = GENERATE(mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON);
  constexpr int n = 3, nz = 4;

  SECTION("Flat elements")
  {
    // The vertices above each other through the stack of layers form one aggregate, and
    // the aggregates of different lines are different.
    auto mesh = MakeLayeredMesh(comm, n, nz, 0.004, type);
    const auto &aggregates = GetVertexAggregates(mesh);
    const auto &pmesh = mesh.Get();
    REQUIRE(static_cast<int>(aggregates.aggregate.size()) == pmesh.GetNV());
    std::map<std::pair<long, long>, std::uint64_t> line_aggregate;
    std::set<std::uint64_t> ids;
    int mismatch = 0, trivial = 0;
    for (int v = 0; v < pmesh.GetNV(); v++)
    {
      const double *x = pmesh.GetVertex(v);
      const std::pair<long, long> line{std::lround(x[0] * n), std::lround(x[1] * n)};
      const auto [it, isnew] = line_aggregate.emplace(line, aggregates.aggregate[v]);
      mismatch += (it->second != aggregates.aggregate[v]);
      trivial += !aggregates.nontrivial[v];
      ids.insert(aggregates.aggregate[v]);
    }
    CHECK(mismatch == 0);
    CHECK(trivial == 0);
    CHECK(ids.size() == line_aggregate.size());

    // The aggregates are cached with the mesh.
    CHECK(&GetVertexAggregates(mesh) == &aggregates);
  }

  SECTION("Isotropic elements")
  {
    auto mesh = MakeLayeredMesh(comm, n, nz, static_cast<double>(nz) / n, type);
    const auto &aggregates = GetVertexAggregates(mesh);
    int nontrivial = 0;
    for (const auto flag : aggregates.nontrivial)
    {
      nontrivial += flag;
    }
    Mpi::GlobalSum(1, &nontrivial, comm);
    CHECK(nontrivial == 0);
  }
}

TEST_CASE("Chebyshev smoothers apply the sweep correction of a relaxation",
          "[blockjacobi][smoother][Serial]")
{
  // The transposed smoother is the transpose of the smoother, with the correction of the
  // relaxation after the polynomial in one and before it in the other.
  const auto comm = Mpi::World();
  constexpr int n = 12;
  TridiagonalOperator A(n);
  const bool fourth_kind = GENERATE(true, false);
  const int num_corrected = GENERATE(0, 3);
  auto MakeChebyshev = [&](int corrected) -> std::unique_ptr<Solver<Operator>>
  {
    auto relaxation = std::make_unique<CorrectedRelaxation>(comm, corrected, 0.2);
    if (fourth_kind)
    {
      return std::make_unique<ChebyshevSmoother<Operator>>(std::move(relaxation), 2, 3,
                                                           1.0);
    }
    return std::make_unique<ChebyshevSmoother1stKind<Operator>>(std::move(relaxation), 2, 3,
                                                                1.0, 0.0);
  };
  auto smoother = MakeChebyshev(num_corrected);
  smoother->SetOperator(A);
  const Eigen::MatrixXd B =
      ToDense(n, [&](const Vector &x, Vector &y) { smoother->Mult(x, y); });
  const Eigen::MatrixXd Bt =
      ToDense(n, [&](const Vector &x, Vector &y) { smoother->MultTranspose(x, y); });
  CHECK((Bt - B.transpose()).cwiseAbs().maxCoeff() < 1.0e-12 * B.cwiseAbs().maxCoeff());
  if (num_corrected > 0)
  {
    // The correction changes the smoother, and makes it nonsymmetric.
    CHECK((B - B.transpose()).cwiseAbs().maxCoeff() > 1.0e-3 * B.cwiseAbs().maxCoeff());
  }
  else
  {
    CHECK((B - B.transpose()).cwiseAbs().maxCoeff() < 1.0e-12 * B.cwiseAbs().maxCoeff());
  }
}

// The dense matrices of the relaxation M⁻¹ of a block Jacobi smoother and of the assembled
// operator, and the blocks of M⁻¹ (the connected components of its sparsity).
struct DenseBlocks
{
  Eigen::MatrixXd Minv, A;
  std::vector<std::vector<int>> blocks;

  Eigen::MatrixXd Block(const Eigen::MatrixXd &M, const std::vector<int> &block) const
  {
    const int m = static_cast<int>(block.size());
    Eigen::MatrixXd B(m, m);
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
      {
        B(i, j) = M(block[i], block[j]);
      }
    }
    return B;
  }
};

DenseBlocks GetDenseBlocks(const BlockJacobiSmoother<Operator> &smoother,
                           const ParOperator &A_full)
{
  DenseBlocks dense;
  const int n = A_full.Height();
  dense.Minv = ToDense(n, [&](const Vector &x, Vector &y) { smoother.Mult(x, y); });
  dense.A = AssembleDense(A_full);
  dense.blocks = Components(dense.Minv, 1.0e-12 * dense.Minv.cwiseAbs().maxCoeff());
  return dense;
}

TEST_CASE("Block Jacobi blocks are the diagonal blocks along the flat direction",
          "[blockjacobi][Serial]")
{
  // The block Jacobi method is the exact inverse of the diagonal blocks of the assembled
  // operator for a partition of the true DOFs, which is read from the sparsity of M⁻¹.
  // Without a definite counterpart, the operator need not be positive definite: the
  // inverse of a block then has the signs of the block.
  const auto comm = Mpi::World();
  const auto [space, order] = GENERATE(std::pair{Space::H1, 1}, std::pair{Space::H1, 2},
                                       std::pair{Space::ND, 1}, std::pair{Space::ND, 2});
  const bool refine = GENERATE(false, true);
  const bool essential = GENERATE(false, true);
  const auto coeff = GENERATE(POSITIVE, NEGATIVE, INDEFINITE);
  CAPTURE(static_cast<int>(space), order, refine, essential, coeff.first, coeff.second);
  constexpr int nx = 2, nz = 3;

  auto mesh = MakeLayeredMesh(comm, nx, nz, 0.006, mfem::Element::TETRAHEDRON, refine);
  const auto fec = MakeFECollection(space, order);
  FiniteElementSpace fespace(mesh, fec.get());
  const int n = fespace.GetTrueVSize();
  mfem::Array<int> ess_list;
  if (essential)
  {
    ess_list = SomeEssentialTrueDofs(fespace);
  }
  auto stiffness = ConstantCoefficient(mesh, coeff.first);
  auto mass = ConstantCoefficient(mesh, coeff.second);
  BilinearForm a(fespace);
  AddIntegrators(a, space, stiffness, mass);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  A.SetEssentialTrueDofs(ess_list, Operator::DIAG_ONE);
  A_full.SetEssentialTrueDofs(ess_list, Operator::DIAG_ONE);
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  const auto dense = GetDenseBlocks(smoother, A_full);

  int num_lines = 0, num_indefinite = 0;
  double error = 0.0;
  for (const auto &block : dense.blocks)
  {
    const int m = static_cast<int>(block.size());
    const Eigen::MatrixXd Ab = dense.Block(dense.A, block),
                          Mb = dense.Block(dense.Minv, block);
    error =
        std::max(error, (Mb * Ab - Eigen::MatrixXd::Identity(m, m)).cwiseAbs().maxCoeff());
    num_lines += (m == nz + 1);
    const Eigen::VectorXd lambda =
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Ab).eigenvalues();
    num_indefinite += (lambda(0) < 0.0 && lambda(m - 1) > 0.0);
    CHECK(m <= BLOCK_SIZE);
  }
  CHECK(error < 1.0e-7);
  CHECK(dense.blocks.size() < static_cast<std::size_t>(n));
  REQUIRE(smoother.HasBlocks());
  const auto &census = smoother.GetBlocks()->GetCensus();
  CHECK(census.dropped == 0);
  if (coeff == POSITIVE)
  {
    CHECK(num_indefinite == 0);
    CHECK(census.negative == 0);
  }
  else if (coeff == NEGATIVE)
  {
    CHECK(census.negative == census.block_dofs);
  }
  else if (space == Space::H1)
  {
    CHECK(num_indefinite > 0);
  }
  if (space == Space::H1 && order == 1 && !refine && !essential)
  {
    // The blocks of the nodal space of lowest order are the lines of vertices through the
    // stack of layers.
    CHECK(static_cast<int>(dense.blocks.size()) == (nx + 1) * (nx + 1));
    CHECK(num_lines == (nx + 1) * (nx + 1));
  }
}

TEST_CASE("Block Jacobi blocks fill the size classes", "[blockjacobi][Serial][Parallel]")
{
  // The blocks of an edge space have at most the given number of DOFs, for every size
  // class: larger groups are split, and what is stored is one triangle of every block
  // padded into its class. A group of the nodal space which is too large is not split, and
  // is solved approximately as a whole.
  const auto comm = Mpi::World();
  const int block_size = GENERATE(2, 4, 8, 16, 32, 64);
  const auto space = GENERATE(Space::ND, Space::H1);
  CAPTURE(block_size, static_cast<int>(space));
  constexpr int nx = 2, nz = 7;
  auto mesh = MakeLayeredMesh(comm, nx, nz, 0.014, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(space, 2);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, space);
  ParOperator A(a.PartialAssemble(), fespace);
  const AggregateBlocks blocks(A, nullptr, block_size, 1.0 / (2.0 * CHEBY_ORDER + 1.0));
  const auto &census = blocks.GetCensus();
  HYPRE_BigInt num_blocks = 0, entries = 0, capacity = 0;
  for (int q = 0; q < AggregateBlocks::NUM_CLASSES; q++)
  {
    const HYPRE_BigInt c = HYPRE_BigInt(1) << q;
    CAPTURE(q, census.blocks[q]);
    CHECK((c <= block_size || census.blocks[q] == 0));
    num_blocks += census.blocks[q];
    entries += census.blocks[q] * (c * (c + 1) / 2);
    capacity += census.blocks[q] * c;
  }
  HYPRE_BigInt num_dofs = blocks.GetBlockDofs().size();
  Mpi::GlobalSum(1, &num_dofs, comm);
  CHECK(num_dofs == census.block_dofs);
  CHECK(census.entries == entries);
  // The padding: a block has more DOFs than half the size of its class.
  CHECK(capacity >= census.block_dofs);
  CHECK((capacity < 2 * census.block_dofs || census.block_dofs == 0));
  if (space == Space::ND)
  {
    // The largest groups, of the edges and faces between two lines of vertices through the
    // stack of 7 layers, have 58 DOFs.
    CHECK(census.sheets == 0);
    CHECK(census.block_dofs > 0);
    CHECK((census.split > 0) == (block_size < 64));
    if (Mpi::Size(comm) == 1)
    {
      // The largest class is used (in parallel, the pieces of a group which is cut at
      // the boundaries of the processes can all be smaller).
      CHECK(census.blocks[static_cast<int>(std::log2(block_size))] > 0);
    }
  }
  else
  {
    // The lines of vertices and edges through the stack have 15 DOFs.
    CHECK((census.sheets > 0) == (block_size < 16));
    CHECK(census.split == 0);
    CHECK((census.sheet_degree > 0) == (census.sheets > 0));
    CHECK((census.block_dofs > 0) == (block_size >= 16));
  }
}

TEST_CASE("Block Jacobi splits large groups into connected blocks", "[blockjacobi][Serial]")
{
  // The pieces of a group of an edge space which is too large are blocks of the operator
  // themselves: M⁻¹ is the inverse of the diagonal blocks of the assembled operator on
  // them, and every piece is connected in the graph of the operator (the DOFs of the
  // edges along a line of vertices, which are one group and not coupled, are blocks of
  // one DOF each at order 1).
  const auto comm = Mpi::World();
  const int block_size = GENERATE(4, 16);
  const int order = GENERATE(1, 2);
  CAPTURE(block_size, order);
  constexpr int nx = 2, nz = 10;
  auto mesh = MakeLayeredMesh(comm, nx, nz, 0.02, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(Space::ND, order);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, Space::ND);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  BlockJacobiSmoother<Operator> smoother(comm, block_size, CHEBY_ORDER);
  smoother.SetOperator(A);
  REQUIRE(smoother.HasBlocks());
  const auto &census = smoother.GetBlocks()->GetCensus();
  CHECK(census.split > 0);
  const auto dense = GetDenseBlocks(smoother, A_full);
  double error = 0.0;
  int num_in_blocks = 0;
  for (const auto &block : dense.blocks)
  {
    const int m = static_cast<int>(block.size());
    REQUIRE(m <= block_size);
    const Eigen::MatrixXd Ab = dense.Block(dense.A, block),
                          Mb = dense.Block(dense.Minv, block);
    error =
        std::max(error, (Mb * Ab - Eigen::MatrixXd::Identity(m, m)).cwiseAbs().maxCoeff());
    num_in_blocks += (m > 1) ? m : 0;
  }
  CHECK(error < 1.0e-7);
  // Every block is connected (the blocks which are counted are those of M⁻¹, where a
  // block which is not connected would show as several), no DOF is lost, and the split
  // does not cut more than it has to: the blocks of more than one DOF are more than half
  // full on average.
  HYPRE_BigInt num_blocks = 0;
  for (const auto num : census.blocks)
  {
    num_blocks += num;
  }
  CHECK(num_in_blocks + census.blocks[0] == census.block_dofs);
  CHECK(2 * (census.block_dofs - census.blocks[0]) >
        (num_blocks - census.blocks[0]) * block_size);
}

TEST_CASE("Block Jacobi blocks of a singular operator have their pseudo-inverse",
          "[blockjacobi][Serial]")
{
  // The blocks of the curl-curl operator without a mass term at order 3 contain gradients,
  // their kernel. M⁻¹ is the pseudo-inverse of every block: M⁻¹ A M⁻¹ = M⁻¹ on it, with as
  // many zero eigenvalues as were dropped.
  const auto comm = Mpi::World();
  auto mesh = MakeLayeredMesh(comm, 2, 2, 0.004, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(Space::ND, 3);
  FiniteElementSpace fespace(mesh, fec.get());
  auto stiffness = ConstantCoefficient(mesh, 1.0);
  BilinearForm a(fespace);
  a.AddDomainIntegrator<CurlCurlIntegrator>(stiffness);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  REQUIRE(smoother.HasBlocks());
  const auto &census = smoother.GetBlocks()->GetCensus();
  CHECK(census.dropped > 0);
  const auto dense = GetDenseBlocks(smoother, A_full);
  double error = 0.0, scale = 0.0;
  HYPRE_BigInt num_zero = 0;
  for (const auto &block : dense.blocks)
  {
    if (block.size() == 1)
    {
      continue;
    }
    const Eigen::MatrixXd Ab = dense.Block(dense.A, block),
                          Mb = dense.Block(dense.Minv, block);
    error = std::max(error, (Mb * Ab * Mb - Mb).cwiseAbs().maxCoeff());
    scale = std::max(scale, Mb.cwiseAbs().maxCoeff());
    const Eigen::VectorXd lambda =
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Mb).eigenvalues();
    num_zero += (lambda.array().abs() < 1.0e-10 * lambda.cwiseAbs().maxCoeff()).count();
  }
  CHECK(error < 1.0e-8 * scale);
  CHECK(num_zero == census.dropped);
}

TEST_CASE("Block Jacobi blocks of an indefinite operator take its signs",
          "[blockjacobi][Serial]")
{
  // With a definite counterpart, the blocks are those of the counterpart with the signs of
  // the operator: on every block the relaxed operator M⁻¹ A has real eigenvalues in
  // (0, 1], where the inverse of the block of the counterpart alone gives negative ones.
  const auto comm = Mpi::World();
  const auto [space, order] = GENERATE(std::pair{Space::H1, 2}, std::pair{Space::ND, 3});
  CAPTURE(static_cast<int>(space), order);
  auto mesh = MakeLayeredMesh(comm, 2, 2, 0.004, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(space, order);
  FiniteElementSpace fespace(mesh, fec.get());
  auto stiffness = ConstantCoefficient(mesh, 1.0);
  auto mass = ConstantCoefficient(mesh, (space == Space::H1) ? -1.0e5 : -1.0);
  auto abs_mass = ConstantCoefficient(mesh, (space == Space::H1) ? 1.0e5 : 1.0);
  BilinearForm a(fespace), a_def(fespace);
  AddIntegrators(a, space, stiffness, mass);
  AddIntegrators(a_def, space, stiffness, abs_mass);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  A.SetDefiniteCounterpart(std::make_unique<ParOperator>(a_def.PartialAssemble(), fespace));
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  REQUIRE(smoother.HasBlocks());
  const auto &census = smoother.GetBlocks()->GetCensus();
  CHECK(census.negative > 0);
  CHECK(census.negative < census.block_dofs);
  const auto dense = GetDenseBlocks(smoother, A_full);
  double min_eig = 1.0, max_eig = 0.0, max_imag = 0.0;
  HYPRE_BigInt num_negative = 0;
  for (const auto &block : dense.blocks)
  {
    if (block.size() == 1)
    {
      continue;
    }
    const Eigen::MatrixXd Ab = dense.Block(dense.A, block),
                          Mb = dense.Block(dense.Minv, block);
    const Eigen::VectorXcd mu = Eigen::EigenSolver<Eigen::MatrixXd>(Mb * Ab).eigenvalues();
    min_eig = std::min(min_eig, mu.real().minCoeff());
    max_eig = std::max(max_eig, mu.real().maxCoeff());
    max_imag = std::max(max_imag, mu.imag().cwiseAbs().maxCoeff());
    num_negative +=
        (Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Mb).eigenvalues().array() < 0.0)
            .count();
  }
  CHECK(min_eig > 0.0);
  CHECK(min_eig < 0.99);  // The blocks of the counterpart, not of the operator
  CHECK(max_eig < 1.0 + 1.0e-8);
  CHECK(max_imag < 1.0e-8);
  CHECK(num_negative == census.negative);
}

TEST_CASE("Block Jacobi blocks of a complex operator are applied between its phases",
          "[blockjacobi][Serial]")
{
  // For a complex symmetric operator A_r + i A_i with a definite counterpart, the blocks
  // are the real ones of A_r, applied between the phases which the diagonal D of the
  // operator has on top of the sign of its real part: M⁻¹ = Ψ M_r⁻¹ Ψ with
  // Ψ² = sign(Re D) D̅ / |D|. A DOF in no block has 1 / D.
  const auto comm = Mpi::World();
  const auto [space, order] = GENERATE(std::pair{Space::H1, 2}, std::pair{Space::ND, 2});
  CAPTURE(static_cast<int>(space), order);
  auto mesh = MakeLayeredMesh(comm, 2, 2, 0.004, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(space, order);
  FiniteElementSpace fespace(mesh, fec.get());
  const int n = fespace.GetTrueVSize();
  auto stiffness = ConstantCoefficient(mesh, 1.0);
  auto mass = ConstantCoefficient(mesh, -1.0e5);
  auto abs_mass = ConstantCoefficient(mesh, 1.0e5);
  auto damping = ConstantCoefficient(mesh, 1.0e6);
  BilinearForm ar(fespace), ai(fespace), a_def(fespace);
  AddIntegrators(ar, space, stiffness, mass);
  AddIntegrators(a_def, space, stiffness, abs_mass);
  if (space == Space::H1)
  {
    ai.AddDomainIntegrator<MassIntegrator>(damping);
  }
  else
  {
    ai.AddDomainIntegrator<VectorFEMassIntegrator>(damping);
  }
  ComplexParOperator A(ar.PartialAssemble(), ai.PartialAssemble(), fespace);
  A.SetDefiniteCounterpart(std::make_unique<ParOperator>(a_def.PartialAssemble(), fespace));
  ParOperator A_real(ar.PartialAssemble(), fespace);
  A_real.SetDefiniteCounterpart(
      std::make_unique<ParOperator>(a_def.PartialAssemble(), fespace));
  BlockJacobiSmoother<ComplexOperator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  BlockJacobiSmoother<Operator> smoother_real(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother_real.SetOperator(A_real);
  REQUIRE(smoother.HasBlocks());
  CHECK(smoother.GetBlocks()->GetCensus().negative > 0);

  // The dense matrices: M⁻¹ column by column, and the diagonal of the operator.
  Eigen::MatrixXcd Minv(n, n);
  {
    ComplexVector x(n), y(n);
    x.UseDevice(true);
    y.UseDevice(true);
    for (int j = 0; j < n; j++)
    {
      x = 0.0;
      x.Real().HostReadWrite()[j] = 1.0;
      smoother.Mult(x, y);
      const auto *YR = y.Real().HostRead(), *YI = y.Imag().HostRead();
      for (int i = 0; i < n; i++)
      {
        Minv(i, j) = std::complex<double>(YR[i], YI[i]);
      }
    }
  }
  const Eigen::MatrixXd Minv_real =
      ToDense(n, [&](const Vector &x, Vector &y) { smoother_real.Mult(x, y); });
  ComplexVector diag(n);
  diag.UseDevice(true);
  A.AssembleDiagonal(diag);
  const auto *DR = diag.Real().HostRead(), *DI = diag.Imag().HostRead();
  std::vector<char> in_block(n, 0);
  for (int i : smoother.GetBlocks()->GetBlockDofs())
  {
    in_block[i] = 1;
  }
  std::vector<std::complex<double>> psi(n);
  for (int i = 0; i < n; i++)
  {
    const std::complex<double> d(DR[i], DI[i]);
    psi[i] = std::sqrt(((d.real() < 0.0) ? -1.0 : 1.0) * std::conj(d) / std::abs(d));
  }
  double error = 0.0, scale = 0.0, phase = 0.0;
  int num_in_blocks = 0;
  for (int i = 0; i < n; i++)
  {
    num_in_blocks += in_block[i];
    phase = std::max(phase, std::abs(psi[i].imag()));
    for (int j = 0; j < n; j++)
    {
      std::complex<double> expected = 0.0;
      if (in_block[i] && in_block[j])
      {
        expected = psi[i] * Minv_real(i, j) * psi[j];
      }
      else if (i == j)
      {
        expected = 1.0 / std::complex<double>(DR[i], DI[i]);
      }
      error = std::max(error, std::abs(Minv(i, j) - expected));
      scale = std::max(scale, std::abs(expected));
    }
  }
  CHECK(num_in_blocks > 0);
  CHECK(phase > 0.1);
  CHECK(error < 1.0e-10 * scale);
}

// A single layer of flat cells and its nodal space of lowest order, where the blocks are
// the pairs of vertices above each other, and operators which are given by a symmetric
// 2 x 2 matrix for each pair: an assembled local matrix with the rows (bottom, across) and
// (across, top) for the two vertices of a pair. Serial.
struct PairSpace
{
  Mesh mesh;
  mfem::H1_FECollection fec;
  FiniteElementSpace fespace;
  std::vector<int> partner;
  std::vector<char> bottom;

  PairSpace(MPI_Comm comm)
    : mesh(MakeLayeredMesh(comm, 2, 1, 0.001, mfem::Element::TETRAHEDRON)), fec(1, 3),
      fespace(mesh, &fec)
  {
    const auto &pmesh = mesh.Get();
    const int n = pmesh.GetNV();
    partner.assign(n, -1);
    bottom.assign(n, 0);
    for (int i = 0; i < n; i++)
    {
      const double *p = pmesh.GetVertex(i);
      bottom[i] = (p[2] == 0.0);
      for (int j = 0; j < n; j++)
      {
        const double *q = pmesh.GetVertex(j);
        if (j != i && q[0] == p[0] && q[1] == p[1])
        {
          partner[i] = j;
        }
      }
      REQUIRE(partner[i] >= 0);
    }
  }

  int Size() const { return static_cast<int>(partner.size()); }

  // The operator with the block (bottom(pair), across(pair); across(pair), top(pair)) for
  // the pair of vertex i, where pair is the index of its bottom vertex.
  template <typename Bottom, typename Across, typename Top>
  std::unique_ptr<Operator> Local(Bottom &&bottom_value, Across &&across_value,
                                  Top &&top_value) const
  {
    const int n = Size();
    auto a = std::make_unique<hypre::HypreCSRMatrix>(n, n, 2 * n);
    for (int i = 0; i < n; i++)
    {
      const int pair = bottom[i] ? i : partner[i];
      a->GetI()[i] = 2 * i;
      a->GetJ()[2 * i] = i;
      a->GetData()[2 * i] = bottom[i] ? bottom_value(pair) : top_value(pair);
      a->GetJ()[2 * i + 1] = partner[i];
      a->GetData()[2 * i + 1] = across_value(pair);
    }
    a->GetI()[n] = 2 * n;
    return a;
  }
};

TEST_CASE("Block Jacobi smoothers replace the point scaling on the rows of the blocks",
          "[blockjacobi][smoother][Serial]")
{
  // The block product is not added to the point Jacobi product and its negative: for a
  // well-conditioned indefinite block with a small (or zero) diagonal entry, the
  // reciprocal of that entry would cancel it. The applications of the smoother are the
  // product with the block inverse on the rows of the blocks, for a real operator and
  // for a complex one (between the phases of its diagonal).
  const auto comm = Mpi::World();
  const double small = GENERATE(std::ldexp(1.0, -60), 0.0);
  CAPTURE(small);
  PairSpace space(comm);
  const int n = space.Size();
  auto Constant = [](double value) { return [value](int) { return value; }; };
  auto MakeOperator = [&]()
  {
    auto A = std::make_unique<ParOperator>(
        space.Local(Constant(small), Constant(1.0), Constant(1.0)), space.fespace);
    A->SetDefiniteCounterpart(std::make_unique<ParOperator>(
        space.Local(Constant(2.0), Constant(0.5), Constant(4.0)), space.fespace));
    return A;
  };
  const auto A = MakeOperator();
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(*A);
  REQUIRE(smoother.HasBlocks());
  REQUIRE(static_cast<int>(smoother.GetBlocks()->GetBlockDofs().size()) == n);

  // The block inverse, which has entries of order one.
  Vector r(n), z(n), y(n);
  r.UseDevice(true);
  z.UseDevice(true);
  y.UseDevice(true);
  r = 1.0;
  z = 0.0;
  smoother.GetBlocks()->MultInverse(r, z);
  {
    const auto *Z = z.HostRead();
    CHECK(*std::max_element(Z, Z + n) < 10.0);
    CHECK(*std::min_element(Z, Z + n) > 0.05);
  }
  auto MaxError = [&](const Vector &u, double su, const Vector &v, double sv, double shift)
  {
    const auto *U = u.HostRead(), *V = v.HostRead();
    double error = 0.0;
    for (int i = 0; i < n; i++)
    {
      error = std::max(error, std::abs(su * U[i] - sv * V[i] - shift));
    }
    return error;
  };
  y = 0.0;
  smoother.Mult(r, y);
  CHECK(MaxError(y, 1.0, z, 1.0, 0.0) < 1.0e-14);
  y = 0.0;
  smoother.ScaledMult(3.0, r, y);
  CHECK(MaxError(y, 1.0, z, 3.0, 0.0) < 1.0e-14);
  y = 5.0;
  smoother.ScaledAddMult(0.5, 2.0, r, y);
  CHECK(MaxError(y, 1.0, z, 2.0, 2.5) < 1.0e-14);

  // The complex operator with the same real part and an imaginary part with the same
  // diagonal: the phases are those of 1 + i on every row (none for a zero diagonal).
  ComplexParOperator A_complex(space.Local(Constant(small), Constant(1.0), Constant(1.0)),
                               space.Local(Constant(small), Constant(0.0), Constant(1.0)),
                               space.fespace);
  A_complex.SetDefiniteCounterpart(std::make_unique<ParOperator>(
      space.Local(Constant(2.0), Constant(0.5), Constant(4.0)), space.fespace));
  BlockJacobiSmoother<ComplexOperator> smoother_complex(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother_complex.SetOperator(A_complex);
  REQUIRE(smoother_complex.HasBlocks());
  // M⁻¹ r = Ψ M_r⁻¹ Ψ r with Ψ² the conjugate phase of the diagonal: -π/4 for 1 + i.
  std::vector<std::complex<double>> psi(n);
  for (int i = 0; i < n; i++)
  {
    const bool no_phase = (small == 0.0 && space.bottom[i]);
    psi[i] = no_phase ? 1.0 : std::polar(1.0, -M_PI / 8.0);
  }
  Vector wr(n), wi(n), zr(n), zi(n);
  for (Vector *v : {&wr, &wi, &zr, &zi})
  {
    v->UseDevice(true);
    *v = 0.0;
  }
  {
    auto *WR = wr.HostReadWrite(), *WI = wi.HostReadWrite();
    for (int i = 0; i < n; i++)
    {
      WR[i] = psi[i].real();
      WI[i] = psi[i].imag();
    }
  }
  smoother.GetBlocks()->MultInverse(wr, zr);
  smoother.GetBlocks()->MultInverse(wi, zi);
  ComplexVector rc(n), yc(n);
  rc.UseDevice(true);
  yc.UseDevice(true);
  rc.Real() = 1.0;
  rc.Imag() = 0.0;
  yc = 0.0;
  smoother_complex.ScaledMult(1.0, rc, yc);
  const auto *ZR = zr.HostRead(), *ZI = zi.HostRead();
  const auto *YR = yc.Real().HostRead(), *YI = yc.Imag().HostRead();
  double error = 0.0, scale = 0.0;
  for (int i = 0; i < n; i++)
  {
    const std::complex<double> expected = psi[i] * std::complex<double>(ZR[i], ZI[i]);
    error = std::max(error, std::abs(std::complex<double>(YR[i], YI[i]) - expected));
    scale = std::max(scale, std::abs(expected));
  }
  CHECK(scale > 0.05);
  CHECK(error < 1.0e-14);
}

TEST_CASE("Block Jacobi smoothers find an eigenvalue which has become dominant",
          "[blockjacobi][smoother][Serial]")
{
  // The relaxed operator is c I on every block when the block of the operator is c times
  // the block of its definite counterpart. With c = 4 on one block and 1 on the others,
  // the Ritz vector of the estimate is an eigenvector on that block. It stays an exact
  // eigenvector, for the eigenvalue 4, when c becomes 100 on another block: the estimate
  // which starts from it still has to find 100.
  const auto comm = Mpi::World();
  PairSpace space(comm);
  int first = -1, second = -1;
  for (int i = 0; i < space.Size(); i++)
  {
    if (space.bottom[i])
    {
      (first < 0 ? first : second) = i;
    }
  }
  auto MakeOperator = [&](double c_second)
  {
    auto Factor = [&, c_second](int pair)
    {
      return (pair == first) ? 4.0 : (pair == second) ? c_second : 1.0;
    };
    auto A = std::make_unique<ParOperator>(
        space.Local([&](int pair) { return 2.0 * Factor(pair); },
                    [&](int pair) { return 0.5 * Factor(pair); },
                    [&](int pair) { return 4.0 * Factor(pair); }),
        space.fespace);
    A->SetDefiniteCounterpart(
        std::make_unique<ParOperator>(space.Local([](int) { return 2.0; }, [](int)
                                                  { return 0.5; }, [](int) { return 4.0; }),
                                      space.fespace));
    return A;
  };
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  Vector x0;
  {
    const auto A = MakeOperator(1.0);
    smoother.SetOperator(*A);
    REQUIRE(smoother.HasBlocks());
    CHECK_THAT(smoother.EstimateMaxEigenvalue(*A, x0), WithinRel(4.0, 1.0e-6));
    // The Ritz vector is supported on the first block.
    REQUIRE(x0.Size() == space.Size());
    const auto *X0 = x0.HostRead();
    double outside = 0.0, inside = 0.0;
    for (int i = 0; i < space.Size(); i++)
    {
      auto &side = (i != first && i != space.partner[first]) ? outside : inside;
      side = std::max(side, std::abs(X0[i]));
    }
    CHECK(outside < 1.0e-8 * inside);
  }
  {
    const auto A = MakeOperator(100.0);
    smoother.SetOperator(*A);
    Vector cold;
    CHECK_THAT(smoother.EstimateMaxEigenvalue(*A, x0), WithinRel(100.0, 1.0e-6));
    CHECK_THAT(smoother.EstimateMaxEigenvalue(*A, cold), WithinRel(100.0, 1.0e-6));
  }
}

TEST_CASE("Block Jacobi blocks are closed under the constraints of hanging entities",
          "[blockjacobi][Serial]")
{
  // On a nonconforming mesh, the DOFs of a column (the local DOFs over each other through
  // the stack of flat elements) are rows of the prolongation P. The strong part of the
  // operator penalizes their differences, so that the functions which it does not see in
  // a block are those for which all DOFs of a column agree: summed over the true DOFs of
  // any block, the rows of P of a column are the same. Checked for the nodal space of
  // lowest order, where the column of a local DOF is the aggregate of its vertex, for the
  // dense blocks and for the groups which are solved approximately.
  const auto comm = Mpi::World();
  const bool refine_half = GENERATE(false, true);
  CAPTURE(refine_half);
  auto mesh =
      MakeLayeredMesh(comm, 4, 3, 0.006, mfem::Element::TETRAHEDRON, true, refine_half);
  mfem::H1_FECollection fec(1, 3);
  FiniteElementSpace fespace(mesh, &fec);
  const int n = fespace.GetTrueVSize(), nl = fespace.GetVSize();
  REQUIRE(n < nl);
  BilinearForm a(fespace);
  AddIntegrators(a, Space::H1);
  ParOperator A(a.PartialAssemble(), fespace);
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  Eigen::MatrixXd Minv =
      ToDense(n, [&](const Vector &x, Vector &y) { smoother.Mult(x, y); }).cwiseAbs();
  Minv /= Minv.maxCoeff();
  // With one half of the mesh refined, the groups at the hanging entities are too large
  // for a block and solved approximately: they are then those of the approximate solve.
  CHECK(smoother.HasSweepCorrection() == refine_half);
  if (smoother.HasSweepCorrection())
  {
    Eigen::MatrixXd S = ToDense(n,
                                [&](const Vector &x, Vector &y)
                                {
                                  y = 0.0;
                                  smoother.AddSweepCorrection(x, y);
                                })
                            .cwiseAbs();
    Minv += S / S.maxCoeff();
  }
  const auto blocks = Components(Minv, 1.0e-12);
  std::vector<int> block_of(n);
  for (std::size_t b = 0; b < blocks.size(); b++)
  {
    for (int t : blocks[b])
    {
      block_of[t] = static_cast<int>(b);
    }
  }

  // Rows of P summed over the blocks, by column.
  mfem::SparseMatrix P;
  fespace.Dof_TrueDof_Matrix()->GetDiag(P);
  const auto &va = GetVertexAggregates(mesh);
  std::map<std::uint64_t, std::map<int, double>> column_weights;
  int num_constrained = 0;
  double error = 0.0;
  for (int i = 0; i < nl; i++)
  {
    if (!va.nontrivial[i])
    {
      continue;
    }
    std::map<int, double> weights;
    for (int k = P.GetI()[i]; k < P.GetI()[i + 1]; k++)
    {
      weights[block_of[P.GetJ()[k]]] += P.GetData()[k];
    }
    num_constrained += (P.RowSize(i) > 1);
    const auto [it, first] = column_weights.try_emplace(va.aggregate[i], weights);
    if (first)
    {
      continue;
    }
    for (const auto &reference : {it->second, weights})
    {
      for (const auto &[b, w] : reference)
      {
        const auto other = weights.find(b), ref = it->second.find(b);
        const double w_other = (other != weights.end()) ? other->second : 0.0;
        const double w_ref = (ref != it->second.end()) ? ref->second : 0.0;
        error = std::max(error, std::abs(w_other - w_ref));
      }
    }
  }
  CHECK(num_constrained > 0);
  CHECK(error < 1.0e-9);
}

TEST_CASE("Block Jacobi blocks across processes are gathered", "[blockjacobi][Parallel]")
{
  // The quadratic form bᵀ M⁻¹ b for a load vector b is the same for the mesh distributed
  // over all processes and for a copy of the mesh on every process: a block with DOFs on
  // more than one process is gathered on one of them. (A group which is too large for a
  // block is cut at the boundaries of the processes instead: these meshes have none, with
  // the edge space at order 1 where its groups grow at the hanging entities.)
  HYPRE_BigInt total_gathered = 0;
  for (const auto space : {Space::H1, Space::ND})
  {
    for (const auto type : {mfem::Element::TETRAHEDRON, mfem::Element::HEXAHEDRON})
    {
      for (const bool refine : {false, true})
      {
        for (const bool essential : {false, true})
        {
          const int order = (space == Space::ND && refine) ? 1 : 2;
          CAPTURE(static_cast<int>(space), order, static_cast<int>(type), refine,
                  essential);
          HYPRE_BigInt gathered = 0, split = 0, gathered_local = 0, split_local = 0;
          const double distributed = SmootherQuadraticForm(
              Mpi::World(), space, order, type, refine, essential, &gathered, &split);
          const double local =
              SmootherQuadraticForm(MPI_COMM_SELF, space, order, type, refine, essential,
                                    &gathered_local, &split_local);
          CHECK(split + split_local == 0);
          CHECK(gathered_local == 0);
          CHECK_THAT(distributed, WithinRel(local, 1.0e-9));
          total_gathered += gathered;
        }
      }
    }
  }
  CHECK(total_gathered > 0);
}

TEST_CASE("Block Jacobi cuts large groups at the boundaries of the processes",
          "[blockjacobi][Serial][Parallel]")
{
  // A group of an edge space which is too large for a block is cut where the processes
  // meet, and its pieces are local: no block is gathered, and the DOFs of the blocks of a
  // process are exactly those of its own true DOFs whose group has at least two of them
  // on this process (a single one is not a block), for the actual partition.
  const auto comm = Mpi::World();
  constexpr int block_size = 4, nx = 3, nz = 6;
  auto mesh = MakeLayeredMesh(comm, nx, nz, 0.012, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(Space::ND, 1);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, Space::ND);
  ParOperator A(a.PartialAssemble(), fespace);
  const AggregateBlocks blocks(A, nullptr, block_size, 1.0 / (2.0 * CHEBY_ORDER + 1.0));
  const auto &census = blocks.GetCensus();
  CHECK(census.split > 0);
  CHECK(census.gathered == 0);
  CHECK(census.sheets == 0);

  // The group of an edge, independently of the smoother: the lines of vertices through
  // the stack (their positions in the plane) which its two vertices are on. Every group
  // of this mesh is too large for a block (the smallest, the edges along a line of
  // vertices, have nz = 6 DOFs).
  const auto &pmesh = mesh.Get();
  const auto &fes = fespace.Get();
  auto Line = [&](int v)
  {
    const double *x = pmesh.GetVertex(v);
    return static_cast<int>(std::lround(nx * x[0])) * (nx + 1) +
           static_cast<int>(std::lround(nx * x[1]));
  };
  std::map<std::pair<int, int>, std::vector<int>> groups;
  mfem::Array<int> edge_dofs, vert;
  for (int e = 0; e < pmesh.GetNEdges(); e++)
  {
    fes.GetEdgeDofs(e, edge_dofs);
    REQUIRE(edge_dofs.Size() == 1);
    const int ldof = (edge_dofs[0] < 0) ? -1 - edge_dofs[0] : edge_dofs[0];
    const int tdof = fes.GetLocalTDofNumber(ldof);
    if (tdof >= 0)
    {
      pmesh.GetEdgeVertices(e, vert);
      const int l0 = Line(vert[0]), l1 = Line(vert[1]);
      groups[{std::min(l0, l1), std::max(l0, l1)}].push_back(tdof);
    }
  }
  std::set<int> expected;
  HYPRE_BigInt num_single = 0;
  for (const auto &[lines, tdofs] : groups)
  {
    if (tdofs.size() > 1)
    {
      expected.insert(tdofs.begin(), tdofs.end());
    }
    else
    {
      num_single++;
    }
  }

  // Ownership and disjointness: the DOFs of the blocks are distinct true DOFs of this
  // process. Coverage: they are the expected ones.
  const auto &dofs = blocks.GetBlockDofs();
  const std::set<int> in_blocks(dofs.begin(), dofs.end());
  CHECK(in_blocks.size() == dofs.size());
  CHECK(std::all_of(dofs.begin(), dofs.end(),
                    [&](int i) { return i >= 0 && i < fespace.GetTrueVSize(); }));
  CHECK(in_blocks == expected);
  HYPRE_BigInt num_local = dofs.size();
  Mpi::GlobalSum(1, &num_local, comm);
  Mpi::GlobalSum(1, &num_single, comm);
  CHECK(num_local == census.block_dofs);
  CHECK(census.block_dofs + num_single == fespace.GlobalTrueVSize());
}

TEST_CASE("The approximate solve of a large nodal group is a fixed definite operator",
          "[blockjacobi][Serial][Parallel]")
{
  // With blocks of two DOFs, the lines of vertices through a stack of layers are solved
  // approximately as a whole: S is symmetric positive definite, and I - S A_b reduces the
  // error in the A_b-norm by the factor asked for, where A_b is the operator between the
  // DOFs of the lines.
  const auto comm = Mpi::World();
  const double reduction = GENERATE(1.0 / 9.0, 1.0 / 13.0);
  CAPTURE(reduction);
  auto mesh = MakeLayeredMesh(comm, 3, 6, 0.012, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(Space::H1, 1);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, Space::H1);
  ParOperator A(a.PartialAssemble(), fespace);
  const AggregateBlocks blocks(A, nullptr, 2, reduction);
  const auto &census = blocks.GetCensus();
  REQUIRE(blocks.HasSheets());
  CHECK(census.sheets == 16);
  CHECK(census.sheet_dofs == 16 * 7);
  CHECK(census.sheet_degree >= 1);
  CHECK(census.sheet_entries > 0);

  const int n = fespace.GetTrueVSize();
  auto Solve = [&](const Vector &r)
  {
    Vector y(n);
    y.UseDevice(true);
    y = 0.0;
    blocks.AddSheetSolve(r, y);
    return y;
  };
  Vector u(n), v(n);
  u.UseDevice(true);
  v.UseDevice(true);
  linalg::SetRandom(comm, u, 1);
  linalg::SetRandom(comm, v, 2);
  const Vector Su = Solve(u), Sv = Solve(v);
  CHECK_THAT(linalg::Dot(comm, v, Su), WithinRel(linalg::Dot(comm, u, Sv), 1.0e-10));
  CHECK(linalg::Dot(comm, u, Su) > 0.0);
  CHECK(linalg::Dot(comm, v, Sv) > 0.0);
  // The same operator at every application.
  const Vector Su2 = Solve(u);
  CHECK(linalg::Dot(comm, Su2, Su2) == linalg::Dot(comm, Su, Su));

  // A_b e for an error e on the lines: the operator applied to e, with the entries
  // between different lines removed. On this mesh every vertex is on a line, and two
  // vertices of different lines are coupled only through elements which a third,
  // separate application for each line would separate; the lines are instead read from
  // the solve itself, serially.
  if (Mpi::Size(comm) == 1)
  {
    ParOperator A_full(a.FullAssemble(false), fespace);
    const Eigen::MatrixXd A_dense = AssembleDense(A_full);
    const Eigen::MatrixXd S = ToDense(n, [&](const Vector &x, Vector &y) { y = Solve(x); });
    const auto lines = Components(S, 1.0e-12 * S.cwiseAbs().maxCoeff());
    CHECK(lines.size() == 16);
    Eigen::MatrixXd Ab = Eigen::MatrixXd::Zero(n, n);
    for (const auto &line : lines)
    {
      for (int i : line)
      {
        for (int j : line)
        {
          Ab(i, j) = A_dense(i, j);
        }
      }
    }
    // The norm of I - S A_b in the A_b-inner product: the largest |1 - λ| over the
    // eigenvalues λ of the symmetric matrix A_b^(1/2) S A_b^(1/2).
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ab);
    const Eigen::MatrixXd root = es.operatorSqrt();
    const Eigen::VectorXd lambda =
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(root * S * root).eigenvalues();
    CHECK(lambda.minCoeff() > 0.0);
    CHECK((1.0 - lambda.array()).abs().maxCoeff() <= reduction);
  }

  // The degree of the polynomial for an interval, and the limits of the work: for the
  // reduction 1 / (2k + 1) of the smoother of order k = 4, the condition number 81, the
  // degree 13 for it, and 165 steps of the Lanczos run.
  CHECK(!census.sheet_limited);
  const auto limits = AggregateBlocks::GetSheetLimits(1.0 / 9.0);
  CHECK_THAT(limits.condition, WithinRel(81.0, 1.0e-12));
  CHECK(limits.degree == 13);
  CHECK(limits.degree == AggregateBlocks::SheetDegree(1.0, 81.0, 1.0 / 9.0));
  CHECK(limits.lanczos_steps == 165);
  CHECK(AggregateBlocks::GetSheetLimits(1.0 / 13.0).degree == 22);
  CHECK(AggregateBlocks::SheetDegree(1.0, 1.0, 0.1) == 1);
  CHECK(AggregateBlocks::SheetDegree(1.0, 4.0, 1.0 / 9.0) == 3);
  CHECK(AggregateBlocks::SheetDegree(1.0, 100.0, 1.0 / 9.0) == 15);
  CHECK(AggregateBlocks::SheetDegree(1.0, 100.0, 1.0 / 13.0) >
        AggregateBlocks::SheetDegree(1.0, 100.0, 1.0 / 9.0));
}

TEST_CASE("The approximate solve of a large nodal group has a bounded degree",
          "[blockjacobi][Serial]")
{
  // On a mesh with one half refined, the nodal group at the hanging entities is a sheet
  // for which one AMG cycle leaves a condition number above the bound (2k + 1)² = 81: the
  // degree is the largest one, 13, and the solve is still a fixed symmetric positive
  // definite operator S which does not amplify any component: the eigenvalues of S A_b
  // are in (0, 2), where A_b is the operator between the DOFs of the sheet.
  const auto comm = Mpi::World();
  constexpr double reduction = 1.0 / 9.0;
  auto mesh = MakeLayeredMesh(comm, 4, 3, 0.006, mfem::Element::TETRAHEDRON, true, true);
  const auto fec = MakeFECollection(Space::H1, 1);
  FiniteElementSpace fespace(mesh, fec.get());
  BilinearForm a(fespace);
  AddIntegrators(a, Space::H1);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  const AggregateBlocks blocks(A, nullptr, BLOCK_SIZE, reduction);
  const auto &census = blocks.GetCensus();
  REQUIRE(census.sheets == 1);
  CHECK(census.sheet_limited);
  CHECK(census.sheet_degree == AggregateBlocks::GetSheetLimits(reduction).degree);
  CHECK(census.sheet_degree == 13);

  const int n = fespace.GetTrueVSize();
  const auto &dofs = blocks.GetSheetDofs();
  const int m = static_cast<int>(dofs.size());
  REQUIRE(m == census.sheet_dofs);
  auto Solve = [&](const Vector &r, Vector &y)
  {
    y = 0.0;
    blocks.AddSheetSolve(r, y);
  };
  const Eigen::MatrixXd S_full = ToDense(n, Solve), S_again = ToDense(n, Solve);
  CHECK((S_full - S_again).cwiseAbs().maxCoeff() == 0.0);
  const Eigen::MatrixXd A_dense = AssembleDense(A_full);
  Eigen::MatrixXd S(m, m), Ab(m, m), S_outside = S_full;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      S(i, j) = S_full(dofs[i], dofs[j]);
      Ab(i, j) = A_dense(dofs[i], dofs[j]);
      S_outside(dofs[i], dofs[j]) = 0.0;
    }
  }
  CHECK(S_outside.cwiseAbs().maxCoeff() == 0.0);  // Nothing outside of the sheet
  CHECK((S - S.transpose()).cwiseAbs().maxCoeff() < 1.0e-10 * S.cwiseAbs().maxCoeff());
  const Eigen::MatrixXd root =
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Ab).operatorSqrt();
  const Eigen::VectorXd lambda = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(
                                     root * (0.5 * (S + S.transpose())) * root)
                                     .eigenvalues();
  CAPTURE(lambda.minCoeff(), lambda.maxCoeff());
  CHECK(lambda.minCoeff() > 0.0);
  CHECK(lambda.maxCoeff() < 2.0);
  // The error polynomial is within the required reduction on the upper part of the
  // spectrum, for which the polynomial is made, and between 0 and 1 below it: here the
  // limit is in effect, and the reduction of the weakest component is less than required.
  CHECK(lambda.maxCoeff() <= 1.0 + reduction + 1.0e-8);
  CHECK(1.0 - lambda.minCoeff() < 1.0);
  CHECK(1.0 - lambda.minCoeff() > reduction);
}

TEST_CASE("Block Jacobi smoothers estimate the spectral radius of the relaxed operator",
          "[blockjacobi][smoother][Serial]")
{
  // The estimate of the largest eigenvalue of M⁻¹A, in magnitude, against the eigenvalues
  // of the dense matrices, for a definite operator and for an indefinite one with its
  // definite counterpart.
  const auto comm = Mpi::World();
  const auto [space, order] = GENERATE(std::pair{Space::H1, 2}, std::pair{Space::ND, 1});
  const bool indefinite = GENERATE(false, true);
  CAPTURE(static_cast<int>(space), order, indefinite);
  auto mesh = MakeLayeredMesh(comm, 2, 3, 0.006, mfem::Element::TETRAHEDRON);
  const auto fec = MakeFECollection(space, order);
  FiniteElementSpace fespace(mesh, fec.get());
  const int n = fespace.GetTrueVSize();
  auto stiffness = ConstantCoefficient(mesh, 1.0);
  auto mass = ConstantCoefficient(mesh, indefinite ? -1.0e5 : 1.0);
  auto abs_mass = ConstantCoefficient(mesh, indefinite ? 1.0e5 : 1.0);
  BilinearForm a(fespace), a_def(fespace);
  AddIntegrators(a, space, stiffness, mass);
  AddIntegrators(a_def, space, stiffness, abs_mass);
  ParOperator A(a.PartialAssemble(), fespace), A_full(a.FullAssemble(false), fespace);
  if (indefinite)
  {
    A.SetDefiniteCounterpart(
        std::make_unique<ParOperator>(a_def.PartialAssemble(), fespace));
  }
  BlockJacobiSmoother<Operator> smoother(comm, BLOCK_SIZE, CHEBY_ORDER);
  smoother.SetOperator(A);
  REQUIRE(smoother.HasBlocks());
  Vector x0;
  const double estimate = smoother.EstimateMaxEigenvalue(A, x0);
  const auto dense = GetDenseBlocks(smoother, A_full);
  const double radius = Eigen::EigenSolver<Eigen::MatrixXd>(dense.Minv * dense.A)
                            .eigenvalues()
                            .cwiseAbs()
                            .maxCoeff();
  CHECK(estimate <= radius * (1.0 + 1.0e-5));
  CHECK(estimate >= 0.98 * radius);
  CHECK(x0.Size() == n);
}

TEST_CASE("Multigrid smoothers without blocks are those of point Jacobi",
          "[blockjacobi][smoother][Serial][Parallel]")
{
  // One cycle of p-multigrid with an AMG coarse solve, with the smoother resolved for the
  // spaces as the linear solver does: with blocks of one DOF the result is that of the
  // point Jacobi smoothers bit for bit, and so it is with larger blocks on a mesh without
  // flat elements, where the smoother is resolved to point Jacobi without a smoother
  // around the coarse solve.
  const auto comm = Mpi::World();
  const bool flat = GENERATE(false, true);
  const int block_size = GENERATE(1, 64);
  CAPTURE(flat, block_size);
  auto mesh = MakeLayeredMesh(comm, 4, 4, flat ? 0.004 : 1.0, mfem::Element::TETRAHEDRON);
  mfem::H1_FECollection fec1(1, 3), fec2(2, 3);
  FiniteElementSpaceHierarchy fespaces(std::make_unique<FiniteElementSpace>(mesh, &fec1));
  fespaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &fec2));
  BilinearForm k(fespaces.GetFinestFESpace());
  k.AddDomainIntegrator<DiffusionIntegrator>();
  k.AddDomainIntegrator<MassIntegrator>();
  auto k_vec = k.Assemble(fespaces, false);
  MultigridOperator K(2);
  for (int l = 0; l < 2; l++)
  {
    K.AddOperator(
        std::make_unique<ParOperator>(std::move(k_vec[l]), fespaces.GetFESpaceAtLevel(l)));
  }
  config::SmootherData configured;
  configured.block_size = block_size;
  configured.cheby_order = 2;
  configured.coarse = (block_size > 1);
  const auto smoother = ResolveSmoother(configured, fespaces);
  CHECK(smoother.block_size == ((flat) ? block_size : 1));
  CHECK(smoother.coarse == (flat && block_size > 1));
  CHECK(smoother.cheby_order == configured.cheby_order);
  auto Cycle = [&](bool resolved)
  {
    auto coarse = std::make_unique<MfemWrapperSolver<Operator>>(
        std::make_unique<BoomerAmgSolver>(1, 1, 1, 0));
    std::unique_ptr<GeometricMultigridSolver<Operator>> pc;
    if (resolved)
    {
      pc = std::make_unique<GeometricMultigridSolver<Operator>>(
          comm, std::move(coarse), fespaces.GetProlongationOperators(), nullptr, 1,
          smoother);
    }
    else
    {
      pc = std::make_unique<GeometricMultigridSolver<Operator>>(
          comm, std::move(coarse), fespaces.GetProlongationOperators(), nullptr, 1, 1, 2,
          1.0, 0.0, true);
    }
    pc->SetOperator(K);
    Vector b(K.Height()), x(K.Height());
    b.UseDevice(true);
    x.UseDevice(true);
    linalg::SetRandom(comm, b, 1);
    pc->Mult(b, x);
    return x;
  };
  const Vector x_blocks = Cycle(true), x_point = Cycle(false);
  Vector diff(x_blocks);
  diff -= x_point;
  const double difference = linalg::Norml2(comm, diff);
  if (block_size == 1 || !flat)
  {
    CHECK(difference == 0.0);
  }
  else
  {
    CHECK(difference > 0.0);
  }
}

TEST_CASE("Multigrid smoothers are resolved before they are constructed",
          "[blockjacobi][smoother][Serial][Parallel]")
{
  // Whether the levels use blocks is decided from the options and the meshes of the
  // spaces. A single level without a smoother around its coarse solve (a sparse direct
  // solver) has no smoother: point Jacobi, and its mesh is not examined.
  const auto comm = Mpi::World();
  mfem::H1_FECollection fec1(1, 3), fec2(2, 3);
  config::SmootherData configured;
  configured.block_size = 64;
  {
    auto mesh = MakeLayeredMesh(comm, 3, 3, 0.003, mfem::Element::TETRAHEDRON);
    FiniteElementSpaceHierarchy fespaces(std::make_unique<FiniteElementSpace>(mesh, &fec1));
    configured.coarse = false;
    const auto direct = ResolveSmoother(configured, fespaces);
    CHECK(direct.block_size == 1);
    CHECK(!direct.coarse);
    CHECK(mesh.GetVertexAggregates() == nullptr);
    configured.coarse = true;
    const auto amg = ResolveSmoother(configured, fespaces);
    CHECK(amg.block_size == 64);
    CHECK(amg.coarse);
    CHECK(mesh.GetVertexAggregates() != nullptr);
  }
  {
    // Two levels, without and with a smoother around the coarse solve.
    auto mesh = MakeLayeredMesh(comm, 3, 3, 0.003, mfem::Element::TETRAHEDRON);
    FiniteElementSpaceHierarchy fespaces(std::make_unique<FiniteElementSpace>(mesh, &fec1));
    fespaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &fec2));
    for (const bool coarse : {false, true})
    {
      configured.coarse = coarse;
      const auto smoother = ResolveSmoother(configured, fespaces);
      CHECK(smoother.block_size == 64);
      CHECK(smoother.coarse == coarse);
    }
    configured.block_size = 1;
    configured.coarse = false;
    CHECK(ResolveSmoother(configured, fespaces).block_size == 1);
  }
}

TEST_CASE("Multigrid with block Jacobi smoothing is robust for flat elements",
          "[blockjacobi][Serial][Parallel]")
{
  // Conjugate gradient iterations for the Poisson problem on a stack of layers,
  // preconditioned by p-multigrid with an AMG coarse solve: with the block Jacobi smoother
  // the number of iterations does not grow when the layers become thin, with the point
  // Jacobi smoother it does. With a single level, the smoother acts around the coarse
  // solve.
  const auto comm = Mpi::World();
  const int num_levels = GENERATE(1, 2);
  const bool refine = GENERATE(false, true);
  auto Iterations = [&](double height, int block_size)
  {
    auto mesh =
        MakeLayeredMesh(comm, 4, 4, height, mfem::Element::TETRAHEDRON, refine, refine);
    mfem::H1_FECollection fec1(1, 3), fec2(2, 3);
    FiniteElementSpaceHierarchy fespaces(std::make_unique<FiniteElementSpace>(mesh, &fec1));
    if (num_levels > 1)
    {
      fespaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &fec2));
    }
    std::vector<mfem::Array<int>> ess_lists(num_levels);
    mfem::Array<int> marker(mesh.Get().bdr_attributes.Max());
    marker = 0;
    marker[0] = 1;
    BilinearForm k(fespaces.GetFinestFESpace());
    k.AddDomainIntegrator<DiffusionIntegrator>();
    auto k_vec = k.Assemble(fespaces, false);
    auto K = std::make_unique<MultigridOperator>(num_levels);
    for (int l = 0; l < num_levels; l++)
    {
      const auto &fespace_l = fespaces.GetFESpaceAtLevel(l);
      fespace_l.Get().GetEssentialTrueDofs(marker, ess_lists[l]);
      auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), fespace_l);
      K_l->SetEssentialTrueDofs(ess_lists[l], Operator::DIAG_ONE);
      K->AddOperator(std::move(K_l));
    }

    // The smoother as the linear solver resolves it for the spaces: without flat
    // elements, point Jacobi and no smoother around the coarse solve.
    config::SmootherData configured;
    configured.block_size = block_size;
    configured.cheby_order = 2;
    configured.coarse = (block_size > 1);
    const auto smoother = ResolveSmoother(configured, fespaces);
    auto pc = std::make_unique<GeometricMultigridSolver<Operator>>(
        comm,
        std::make_unique<MfemWrapperSolver<Operator>>(
            std::make_unique<BoomerAmgSolver>(1, 1, 1, 0)),
        fespaces.GetProlongationOperators(), nullptr, 1, smoother);
    auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
    pcg->SetInitialGuess(false);
    pcg->SetRelTol(1.0e-8);
    pcg->SetMaxIter(500);
    KspSolver ksp(std::move(pcg), std::move(pc));
    ksp.SetOperators(*K, *K);

    Vector b(K->Height()), x(K->Height());
    b.UseDevice(true);
    x.UseDevice(true);
    linalg::SetRandom(comm, b, 1);
    linalg::SetSubVector(b, ess_lists.back(), 0.0);
    ksp.Mult(b, x);
    return ksp.NumTotalMultIterations();
  };
  const int block_isotropic = Iterations(1.0, BLOCK_SIZE);
  const int block_flat = Iterations(0.004, BLOCK_SIZE);
  const int small_flat = Iterations(0.004, 2);
  const int point_flat = Iterations(0.004, 1);
  CAPTURE(num_levels, refine, block_isotropic, block_flat, small_flat, point_flat);
  CHECK(block_flat <= block_isotropic + 4);
  CHECK(small_flat <= block_isotropic + 6);
  if (num_levels > 1)
  {
    CHECK(point_flat >= 2 * block_flat);
  }
  else
  {
    // The AMG solve alone is already a good preconditioner for this problem.
    CHECK(point_flat > block_flat);
  }
}

}  // namespace palace
