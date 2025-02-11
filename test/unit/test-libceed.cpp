// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <memory>
#include <sstream>
#include <string>
#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/hypre.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"

extern int benchmark_ref_levels;
extern int benchmark_order;
extern bool benchmark_assemble_q_data;
extern bool benchmark_no_fa;
extern bool benchmark_no_mfem_pa;

namespace palace
{

namespace
{

auto Initialize(MPI_Comm comm, const std::string &input, int ref_levels, bool amr)
{
  // Load the mesh.
  mfem::Mesh smesh(input, 1, 1);
  smesh.EnsureNodes();

  // Configure attributes for piecewise coefficients (input mesh is always conformal, so
  // this is OK).
  const int max_attr = (smesh.GetNE() + 1) / 2;
  const int max_bdr_attr = (smesh.GetNBE() + 1) / 2;
  for (int i = 0; i < smesh.GetNE(); i++)
  {
    smesh.SetAttribute(i, 1 + (i % max_attr));
  }
  for (int i = 0; i < smesh.GetNBE(); i++)
  {
    smesh.SetBdrAttribute(i, 1 + (i % max_bdr_attr));
  }
  smesh.SetAttributes();

  // Construct nonconforming mesh for AMR.
  if (amr)
  {
    smesh.EnsureNCMesh(true);
  }

  // Construct the parallel mesh.
  REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
  auto pmesh = std::make_unique<mfem::ParMesh>(comm, smesh);
  for (int l = 0; l < ref_levels; l++)
  {
    pmesh->UniformRefinement();
  }

  // Perform nonconforming AMR (two levels of refinement with no hanging node restritions).
  if (amr)
  {
    pmesh->RandomRefinement(0.5);
    pmesh->RandomRefinement(0.5);
  }

  return Mesh(std::move(pmesh));
}

enum class CoeffType
{
  Const,
  Scalar,
  Matrix
};

auto ToString(CoeffType type)
{
  switch (type)
  {
    case CoeffType::Const:
      return "Constant";
    case CoeffType::Scalar:
      return "Scalar";
    case CoeffType::Matrix:
      return "Matrix";
  }
  return "";
}

class PWCoefficient : public mfem::Coefficient, public mfem::MatrixCoefficient
{
private:
  mfem::DenseTensor C;

public:
  PWCoefficient(const mfem::DenseTensor &C)
    : mfem::Coefficient(), mfem::MatrixCoefficient(C.SizeI(), C.SizeJ()), C(C)
  {
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    return C(0, 0, T.Attribute - 1);
  }

  void Eval(mfem::DenseMatrix &K, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    K = C(T.Attribute - 1);
  }
};

void BuildCoefficientHelper(const mfem::Mesh &mesh, bool bdr_integ, CoeffType coeff_type,
                            mfem::Array<int> &attr_mat, mfem::DenseTensor &mat_coeff)
{
  // Assign material properties to domain or boundary attributes, based on the global
  // attributes.
  const auto &attributes = bdr_integ ? mesh.bdr_attributes : mesh.attributes;
  attr_mat.SetSize(attributes.Size() ? attributes.Max() : 0);
  const int num_mat = std::min(attributes.Size() ? attributes.Max() : 0, 4);
  for (int i = 0; i < attributes.Size(); i++)
  {
    attr_mat[i] = i % num_mat;
  }

  // Generate material properties.
  const auto dim = (coeff_type == CoeffType::Scalar) ? 1 : mesh.Dimension();
  mat_coeff.SetSize(dim, dim, num_mat);
  for (int k = 0; k < num_mat; k++)
  {
    mat_coeff(k) = 0.1;
    for (int d = 0; d < dim; d++)
    {
      mat_coeff(d, d, k) = 10.0 * k + (d + 1.0);
    }
  }
}

auto BuildCoefficient(const Mesh &mesh, bool bdr_integ, CoeffType coeff_type)
{
  if (coeff_type == CoeffType::Const)
  {
    return MaterialPropertyCoefficient(0);
  }
  mfem::Array<int> attr_mat;
  mfem::DenseTensor mat_coeff;
  BuildCoefficientHelper(mesh, bdr_integ, coeff_type, attr_mat, mat_coeff);

  // Convert attribute to material mapping from global MFEM attributes to local libCEED
  // ones.
  mfem::Array<int> loc_attr_mat(bdr_integ ? mesh.MaxCeedBdrAttribute()
                                          : mesh.MaxCeedAttribute());
  loc_attr_mat = -1;
  for (int i = 0; i < attr_mat.Size(); i++)
  {
    for (auto attr :
         (bdr_integ ? mesh.GetCeedBdrAttributes(i + 1) : mesh.GetCeedAttributes(i + 1)))
    {
      loc_attr_mat[attr - 1] = attr_mat[i];
    }
  }
  return MaterialPropertyCoefficient(loc_attr_mat, mat_coeff);
}

auto BuildCoefficientRef(const Mesh &mesh, bool bdr_integ, CoeffType coeff_type)
{
  if (coeff_type == CoeffType::Const)
  {
    return PWCoefficient(mfem::DenseTensor());
  }
  mfem::Array<int> attr_mat;
  mfem::DenseTensor mat_coeff;
  BuildCoefficientHelper(mesh, bdr_integ, coeff_type, attr_mat, mat_coeff);

  mfem::DenseTensor C(mat_coeff.SizeI(), mat_coeff.SizeJ(), attr_mat.Size());
  for (int i = 0; i < attr_mat.Size(); i++)
  {
    C(i) = mat_coeff(attr_mat[i]);
  }
  return PWCoefficient(C);
}

template <typename T1, typename T2, typename U>
void AddIntegrators(bool bdr_integ, BilinearForm &a_test, U &a_ref)
{
  if (bdr_integ)
  {
    a_test.AddBoundaryIntegrator<T1>();
    a_ref.AddBoundaryIntegrator(new T2());
  }
  else
  {
    a_test.AddDomainIntegrator<T1>();
    a_ref.AddDomainIntegrator(new T2());
  }
}

template <typename T1, typename T2, typename U, typename V>
void AddIntegrators(bool bdr_integ, BilinearForm &a_test, U &a_ref,
                    MaterialPropertyCoefficient &Q, V &Q_ref)
{
  if (bdr_integ)
  {
    a_test.AddBoundaryIntegrator<T1>(Q);
    a_ref.AddBoundaryIntegrator(new T2(Q_ref));
  }
  else
  {
    a_test.AddDomainIntegrator<T1>(Q);
    a_ref.AddDomainIntegrator(new T2(Q_ref));
  }
}

void TestCeedOperatorMult(const Operator &op_test, const Operator &op_ref,
                          bool test_transpose, double scaling = 1.0)
{
  Vector x(op_ref.Width()), y_ref(op_ref.Height()), y_test(op_ref.Height());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  y_test.UseDevice(true);
  {
    x.Randomize(1);

    op_ref.Mult(x, y_ref);
    op_test.Mult(x, y_test);

    y_test *= scaling;
    y_test -= y_ref;

    // REQUIRE(y_ref * y_ref > 0.0);
    REQUIRE(y_test * y_test < 1.0e-12 * std::max(y_ref * y_ref, 1.0));
  }
  if (test_transpose)
  {
    Vector x_t(op_ref.Height()), y_t_ref(op_ref.Width()), y_t_test(op_ref.Width());
    x_t.UseDevice(true);
    y_t_ref.UseDevice(true);
    y_t_test.UseDevice(true);

    x_t.Randomize(1);

    op_ref.MultTranspose(x_t, y_t_ref);
    op_test.MultTranspose(x_t, y_t_test);

    y_t_test *= scaling;
    y_t_test -= y_t_ref;

    // REQUIRE(y_t_ref * y_t_ref > 0.0);
    REQUIRE(y_t_test * y_t_test < 1.0e-12 * std::max(y_t_ref * y_t_ref, 1.0));
  }
}

void TestCeedOperatorFullAssemble(mfem::SparseMatrix &mat_test, mfem::SparseMatrix &mat_ref,
                                  double scaling = 1.0)
{
  // Ensure host memory is up to date (mfem::Add is missing the device to host copy).
  mat_test.HostReadI();
  mat_test.HostReadJ();
  mat_test.HostReadData();
  mat_ref.HostReadI();
  mat_ref.HostReadJ();
  mat_ref.HostReadData();

  std::unique_ptr<mfem::SparseMatrix> mat_diff(mfem::Add(scaling, mat_test, -1.0, mat_ref));

  // REQUIRE(mat_ref.MaxNorm() > 0.0);
  REQUIRE(mat_diff->MaxNorm() < 1.0e-12 * std::max(mat_ref.MaxNorm(), 1.0));
}

void TestCeedOperatorFullAssemble(hypre::HypreCSRMatrix &mat_test,
                                  mfem::SparseMatrix &mat_ref, double scaling = 1.0)
{
  // Copy test matrix into MFEM's sparse matrix data type.
  hypre_CSRMatrixMigrate(mat_test, HYPRE_MEMORY_HOST);
  mfem::SparseMatrix mat_test_sp(mat_test.GetI(), mat_test.GetJ(), mat_test.GetData(),
                                 mat_test.Height(), mat_test.Width(), false, false, false);

  // Perform the test.
  TestCeedOperatorFullAssemble(mat_test_sp, mat_ref, scaling);
}

void TestCeedOperatorFullAssemble(hypre::HypreCSRMatrix &mat_test,
                                  hypre::HypreCSRMatrix &mat_ref, double scaling = 1.0)
{
  // Copy test and reference matrix into MFEM's sparse matrix data type.
  hypre_CSRMatrixMigrate(mat_test, HYPRE_MEMORY_HOST);
  hypre_CSRMatrixMigrate(mat_ref, HYPRE_MEMORY_HOST);
  mfem::SparseMatrix mat_test_sp(mat_test.GetI(), mat_test.GetJ(), mat_test.GetData(),
                                 mat_test.Height(), mat_test.Width(), false, false, false);
  mfem::SparseMatrix mat_ref_sp(mat_ref.GetI(), mat_ref.GetJ(), mat_ref.GetData(),
                                mat_ref.Height(), mat_ref.Width(), false, false, false);

  // Perform the test.
  TestCeedOperatorFullAssemble(mat_test_sp, mat_ref_sp, scaling);
}

template <typename T1, typename T2>
void TestCeedOperator(T1 &a_test, T2 &a_ref, bool test_transpose, bool skip_zeros,
                      double scaling = 1.0)
{
  a_ref.Assemble(skip_zeros);
  a_ref.Finalize(skip_zeros);
  auto *mat_ref = &a_ref.SpMat();
  auto *op_ref = mat_ref;

  // Test operator application.
  auto op_test = a_test.PartialAssemble();
  TestCeedOperatorMult(*op_test, *op_ref, test_transpose, scaling);

  // Test full assembly.
  auto mat_test = a_test.FullAssemble(*op_test, skip_zeros);
  TestCeedOperatorFullAssemble(*mat_test, *mat_ref, scaling);

  // Test diagonal assembly if possible.
  if (&a_test.GetTrialSpace() == &a_test.GetTestSpace())
  {
    Vector d_ref(mat_ref->Height()), d_test(mat_ref->Height());
    d_ref.UseDevice(true);
    d_test.UseDevice(true);

    mat_ref->GetDiag(d_ref);
    op_test->AssembleDiagonal(d_test);

    d_test *= scaling;
    d_test -= d_ref;

    // Diagonal assembly for high-order Nedelec spaces is only approximate due to face
    // dofs in 3D.
    double rtol = 1.0e-12;
    const auto &trial_fespace = a_test.GetTrialSpace();
    const auto &test_fespace = a_test.GetTestSpace();
    const auto &trial_fec = trial_fespace.GetFEColl();
    const auto &test_fec = test_fespace.GetFEColl();
    if (trial_fespace.Dimension() == 3 &&
        ((dynamic_cast<const mfem::ND_FECollection *>(&trial_fec) &&
          trial_fec.GetOrder() > 1 && !mfem::UsesTensorBasis(trial_fespace)) ||
         (dynamic_cast<const mfem::ND_FECollection *>(&test_fec) &&
          test_fec.GetOrder() > 1 && !mfem::UsesTensorBasis(test_fespace))))
    {
      rtol = 1.0;
    }

    // REQUIRE(d_ref * d_ref > 0.0);
    REQUIRE(d_test * d_test < rtol * std::max(d_ref * d_ref, 1.0));
  }
}

void TestCeedOperator(BilinearForm &op_test, mfem::BilinearForm &op_ref,
                      double scaling = 1.0)
{
  TestCeedOperator(op_test, op_ref, false, false, scaling);
}

void TestCeedOperator(BilinearForm &op_test, mfem::MixedBilinearForm &op_ref,
                      double scaling = 1.0)
{
  TestCeedOperator(op_test, op_ref, false, false, scaling);
}

void TestCeedOperator(DiscreteLinearOperator &op_test, mfem::DiscreteLinearOperator &op_ref,
                      double scaling = 1.0)
{
  TestCeedOperator(op_test, op_ref, true, true, scaling);
}

template <typename T1, typename T2, typename T3>
void BenchmarkCeedIntegrator(FiniteElementSpace &fespace, T1 AssembleTest,
                             T2 AssembleTestRef, T3 AssembleRef, int q_data_size)
{
  const bool skip_zeros = false;
  Vector x(fespace.GetVSize()), y_ref(fespace.GetVSize()), y_test(fespace.GetVSize());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  y_test.UseDevice(true);
  x.Randomize(1);

  // Check correctness (with boundary integrators).
  std::size_t nnz = 0;
  if (!benchmark_no_fa)
  {
    constexpr bool bdr_integ = true;
    auto op_test = AssembleTest(fespace, bdr_integ);
    auto op_test_ref = AssembleTestRef(fespace, bdr_integ);
    auto mat_test = BilinearForm::FullAssemble(*op_test, skip_zeros);
    auto mat_test_ref = BilinearForm::FullAssemble(*op_test_ref, skip_zeros);
    nnz = mat_test->NNZ();
    TestCeedOperatorFullAssemble(*mat_test, *mat_test_ref);
  }

  // Benchmark MFEM legacy assembly.
  if (!benchmark_no_fa)
  {
    BENCHMARK("Assemble (MFEM Legacy)")
    {
      auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      return op_ref->Height();
    };
    {
      auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Legacy)")
      {
        op_ref->AddMult(x, y_ref);
        return y_ref.Size();
      };
    }
  }

  // Benchmark MFEM PA (tensor-product elements only).
  if (!benchmark_no_mfem_pa && mfem::UsesTensorBasis(fespace))
  {
    BENCHMARK("Assemble (MFEM Partial)")
    {
      auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::PARTIAL, skip_zeros);
      return op_ref->Height();
    };
    {
      auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::PARTIAL, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Partial)")
      {
        // MFEM PA does not implement AddMult from BilinearForm.
        op_ref->Mult(x, y_test);
        y_ref += y_test;
        return y_ref.Size();
      };
    }
  }

  // Benchmark libCEED assembly.
  BENCHMARK("Assemble (libCEED)")
  {
    auto op_test = AssembleTest(fespace);
    return op_test->Height();
  };
  {
    auto op_test = AssembleTest(fespace);
    y_test = 0.0;
    BENCHMARK("AddMult (libCEED)")
    {
      op_test->AddMult(x, y_test);
      return y_test.Size();
    };
  }
  if (!benchmark_no_fa)
  {
    BENCHMARK("Full Assemble (libCEED)")
    {
      auto op_test = AssembleTest(fespace);
      auto mat_test = BilinearForm::FullAssemble(*op_test, skip_zeros);
      return mat_test->NNZ();
    };
  }

  // Memory estimate (only for non-mixed meshes).
  mfem::ParMesh &mesh = fespace.GetParMesh();
  if (mesh.GetNumGeometries(mesh.Dimension()) == 1)
  {
    // Integration rule gives the complete non-tensor number of points.
    const mfem::FiniteElement &fe = *fespace.Get().GetFE(0);
    const mfem::ElementTransformation &T = *mesh.GetElementTransformation(0);
    const int q_order = fem::DefaultIntegrationOrder::Get(T);
    const int Q = mfem::IntRules.Get(mesh.GetElementGeometry(0), q_order).GetNPoints();
    const int P = fe.GetDof();

    // Rough estimate for memory consumption as quadrature data + offsets for element
    // restriction.
    std::size_t mem_ref = nnz * (8 + 4) + (y_ref.Size() + 1) * 4;
    std::size_t mem_test = (Q * q_data_size * 8 + P * 4) * (std::size_t)mesh.GetNE();
    std::stringstream msg;
    msg << "benchmark memory estimate:\n"
        << "  N = " << fespace.GetVSize() << " (NE = " << mesh.GetNE() << ", P = " << P
        << ", Q = " << Q << ")\n";
    if (nnz > 0)
    {
      msg << "  Full Assembly = " << mem_ref / (double)(1024 * 1024) << " MB (" << nnz
          << " NNZ)\n";
    }
    else
    {
      msg << "  Full Assembly = N/A (skipped)\n";
    }
    msg << "  Partial Assembly = " << mem_test / (double)(1024 * 1024) << " MB\n";
    WARN(msg.str());
  }
}

template <typename T1, typename T2>
void BenchmarkCeedInterpolator(FiniteElementSpace &trial_fespace,
                               FiniteElementSpace &test_fespace, T1 AssembleTest,
                               T2 AssembleRef)
{
  const bool skip_zeros_interp = true;
  Vector x(trial_fespace.GetVSize()), y_ref(test_fespace.GetVSize()),
      y_test(test_fespace.GetVSize());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  y_test.UseDevice(true);
  x.Randomize(1);

  // Check correctness.
  std::size_t nnz = 0;
  if (!benchmark_no_fa)
  {
    auto op_test = AssembleTest(trial_fespace, test_fespace);
    auto op_ref = AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY,
                              skip_zeros_interp);
    auto mat_test = DiscreteLinearOperator::FullAssemble(*op_test, skip_zeros_interp);
    auto *mat_ref = &op_ref->SpMat();
    nnz = mat_test->NNZ();
    TestCeedOperatorFullAssemble(*mat_test, *mat_ref);
  }

  // Benchmark MFEM legacy assembly.
  if (!benchmark_no_fa)
  {
    BENCHMARK("Assemble (MFEM Legacy)")
    {
      auto op_ref = AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY,
                                skip_zeros_interp);
      return op_ref->Height();
    };
    {
      auto op_ref = AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY,
                                skip_zeros_interp);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Legacy)")
      {
        op_ref->AddMult(x, y_ref);
        return y_ref.Size();
      };
    }
  }

  // Benchmark MFEM PA (tensor-product elements only).
  if (!benchmark_no_mfem_pa && mfem::UsesTensorBasis(trial_fespace) &&
      mfem::UsesTensorBasis(test_fespace))
  {
    BENCHMARK("Assemble (MFEM Partial)")
    {
      auto op_ref = AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::PARTIAL,
                                skip_zeros_interp);
      return op_ref->Height();
    };
    {
      auto op_ref = AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::PARTIAL,
                                skip_zeros_interp);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Partial)")
      {
        // MFEM PA does not implement AddMult from BilinearForm.
        op_ref->Mult(x, y_test);
        y_ref += y_test;
        return y_ref.Size();
      };
    }
  }

  // Benchmark libCEED assembly.
  BENCHMARK("Assemble (libCEED)")
  {
    auto op_test = AssembleTest(trial_fespace, test_fespace);
    return op_test->Height();
  };
  {
    auto op_test = AssembleTest(trial_fespace, test_fespace);
    y_test = 0.0;
    BENCHMARK("AddMult (libCEED)")
    {
      op_test->AddMult(x, y_test);
      return y_test.Size();
    };
  }
  if (!benchmark_no_fa)
  {
    BENCHMARK("Full Assemble (libCEED)")
    {
      auto op_test = AssembleTest(trial_fespace, test_fespace);
      auto mat_test = DiscreteLinearOperator::FullAssemble(*op_test, skip_zeros_interp);
      return mat_test->NNZ();
    };
  }

  // Memory estimate (only for non-mixed meshes).
  mfem::ParMesh &mesh = trial_fespace.GetParMesh();
  if (mesh.GetNumGeometries(mesh.Dimension()) == 1)
  {
    const mfem::FiniteElement &trial_fe = *trial_fespace.Get().GetFE(0);
    const mfem::FiniteElement &test_fe = *test_fespace.Get().GetFE(0);
    const int trial_P = trial_fe.GetDof();
    const int test_P = test_fe.GetDof();

    // Rough estimate for memory consumption as quadrature data + offsets for element
    // restriction.
    std::size_t mem_ref = nnz * (8 + 4) + (y_ref.Size() + 1) * 4;
    std::size_t mem_test = (trial_P * 4 + test_P * 4) * (std::size_t)mesh.GetNE();
    std::stringstream msg;
    msg << "benchmark memory estimate:\n"
        << "  N = " << trial_fespace.GetVSize() << ", " << test_fespace.GetVSize()
        << " (NE = " << mesh.GetNE() << ", P = " << trial_P << ", " << test_P << ")\n";
    if (nnz > 0)
    {
      msg << "  Full Assembly = " << mem_ref / (double)(1024 * 1024) << " MB (" << nnz
          << " NNZ)\n";
    }
    else
    {
      msg << "  Full Assembly = N/A (skipped)\n";
    }
    msg << "  Partial Assembly = " << mem_test / (double)(1024 * 1024) << " MB\n";
    WARN(msg.str());
  }
}

void RunCeedIntegratorTests(MPI_Comm comm, const std::string &input, int ref_levels,
                            bool amr, int order)
{
  // Load the mesh.
  auto mesh = Initialize(comm, input, ref_levels, amr);
  const int dim = mesh.Dimension();

  // Match MFEM's default integration orders.
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Run the tests.
  auto bdr_integ = GENERATE(false, true);
  auto coeff_type = GENERATE(CoeffType::Const, CoeffType::Scalar, CoeffType::Matrix);
  std::string section =
      "Mesh: " + input + "\n" + "Refinement levels: " + std::to_string(ref_levels) + "\n" +
      "AMR: " + std::to_string(amr) + "\n" + "Order: " + std::to_string(order) + "\n" +
      "Integrator: " + (bdr_integ ? "Boundary" : "Domain") + "\n" +
      "Coefficient: " + ToString(coeff_type) + "\n";
  INFO(section);

  // Initialize coefficients.
  auto Q = BuildCoefficient(mesh, bdr_integ, coeff_type);
  auto Q_ref = BuildCoefficientRef(mesh, bdr_integ, coeff_type);

  // Tests on H1 spaces.
  SECTION("H1 Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), h1d_fespace(mesh, &h1_fec, dim);
    SECTION("H1 Mass Integrator")
    {
      BilinearForm a_test(h1_fespace);
      mfem::BilinearForm a_ref(&h1_fespace.Get());
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref);
          break;
        case CoeffType::Scalar:
          AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                               Q_ref);
          break;
        case CoeffType::Matrix:
          break;  // Good to test empty operators
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Vector H1 Mass Integrator")
    {
      BilinearForm a_test(h1d_fespace);
      mfem::BilinearForm a_ref(&h1d_fespace.Get());
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators<MassIntegrator, mfem::VectorMassIntegrator>(bdr_integ, a_test,
                                                                     a_ref);
          break;
        case CoeffType::Scalar:
          AddIntegrators<MassIntegrator, mfem::VectorMassIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
          break;
        case CoeffType::Matrix:
          AddIntegrators<MassIntegrator, mfem::VectorMassIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("H1 Diffusion Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_jac = false;
      fem::DefaultIntegrationOrder::q_order_extra_pk = -2;
      fem::DefaultIntegrationOrder::q_order_extra_qk = dim - bdr_integ - 1;
      mesh.ResetCeedObjects();
      BilinearForm a_test(h1_fespace);
      mfem::BilinearForm a_ref(&h1_fespace.Get());
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators<DiffusionIntegrator, mfem::DiffusionIntegrator>(bdr_integ, a_test,
                                                                         a_ref);
          break;
        case CoeffType::Scalar:
          AddIntegrators<DiffusionIntegrator, mfem::DiffusionIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
          break;
        case CoeffType::Matrix:
          AddIntegrators<DiffusionIntegrator, mfem::DiffusionIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on H(curl) spaces.
  SECTION("H(curl) Integrators")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace nd_fespace(mesh, &nd_fec);
    SECTION("ND Mass Integrator")
    {
      BilinearForm a_test(nd_fespace);
      mfem::BilinearForm a_ref(&nd_fespace.Get());
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
              bdr_integ, a_test, a_ref);
          break;
        case CoeffType::Scalar:
          AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
          break;
        case CoeffType::Matrix:
          AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
              bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("ND Curl-Curl Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_jac = false;
      fem::DefaultIntegrationOrder::q_order_extra_pk = -2;
      fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
      mesh.ResetCeedObjects();
      BilinearForm a_test(nd_fespace);
      mfem::BilinearForm a_ref(&nd_fespace.Get());
      if (dim == 3 || (dim == 2 && !bdr_integ))  // No 1D ND curl shape
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<CurlCurlIntegrator, mfem::CurlCurlIntegrator>(bdr_integ, a_test,
                                                                         a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<CurlCurlIntegrator, mfem::CurlCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            if (dim == 3 && !bdr_integ)
            {
              AddIntegrators<CurlCurlIntegrator, mfem::CurlCurlIntegrator>(
                  bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            }
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on H(div) spaces.
  SECTION("H(div) Integrators")
  {
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace rt_fespace(mesh, &rt_fec);
    SECTION("RT Mass Integrator")
    {
      BilinearForm a_test(rt_fespace);
      mfem::BilinearForm a_ref(&rt_fespace.Get());
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
                bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<VectorFEMassIntegrator, mfem::VectorFEMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("RT Div-Div Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_jac = false;
      fem::DefaultIntegrationOrder::q_order_extra_pk = -2;
      fem::DefaultIntegrationOrder::q_order_extra_qk = -2;
      mesh.ResetCeedObjects();
      BilinearForm a_test(rt_fespace);
      mfem::BilinearForm a_ref(&rt_fespace.Get());
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<DivDivIntegrator, mfem::DivDivIntegrator>(bdr_integ, a_test,
                                                                     a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<DivDivIntegrator, mfem::DivDivIntegrator>(bdr_integ, a_test,
                                                                     a_ref, Q, Q_ref);
            break;
          case CoeffType::Matrix:
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on mixed H1-H(curl) spaces.
  SECTION("H1-H(curl) Mixed Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);
    SECTION("Mixed Vector Gradient Integrator")
    {
      BilinearForm a_test(h1_fespace, nd_fespace);
      mfem::MixedBilinearForm a_ref(&h1_fespace.Get(), &nd_fespace.Get());
      if (dim == 3 || (dim == 2 && !bdr_integ))  // Only in 2D or 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                                (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Weak Divergence Integrator")
    {
      BilinearForm a_test(nd_fespace, h1_fespace);
      mfem::MixedBilinearForm a_ref(&nd_fespace.Get(), &h1_fespace.Get());
      if (dim == 3 || (dim == 2 && !bdr_integ))  // Only in 2D or 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorWeakDivergenceIntegrator,
                           mfem::MixedVectorWeakDivergenceIntegrator>(bdr_integ, a_test,
                                                                      a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorWeakDivergenceIntegrator,
                           mfem::MixedVectorWeakDivergenceIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorWeakDivergenceIntegrator,
                           mfem::MixedVectorWeakDivergenceIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on mixed H1-H(div) spaces.
  SECTION("H1-H(div) Mixed Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), rt_fespace(mesh, &rt_fec);
    SECTION("Mixed Vector Gradient Integrator")
    {
      BilinearForm a_test(h1_fespace, rt_fespace);
      mfem::MixedBilinearForm a_ref(&h1_fespace.Get(), &rt_fespace.Get());
      if (!bdr_integ)
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                                (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorGradientIntegrator,
                           mfem::MixedVectorGradientIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
        TestCeedOperator(a_test, a_ref);
      }
    }
  }

  // Tests on mixed H(curl)-H(div) spaces.
  SECTION("H(curl)-H(div) Mixed Integrators")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace nd_fespace(mesh, &nd_fec), rt_fespace(mesh, &rt_fec);
    SECTION("Mixed H(curl)-H(div) Mass Integrator")
    {
      BilinearForm a_test(nd_fespace, rt_fespace);
      mfem::MixedBilinearForm a_ref(&nd_fespace.Get(), &rt_fespace.Get());
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed H(div)-H(curl) Mass Integrator")
    {
      BilinearForm a_test(rt_fespace, nd_fespace);
      mfem::MixedBilinearForm a_ref(&rt_fespace.Get(), &nd_fespace.Get());
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<VectorFEMassIntegrator, mfem::MixedVectorMassIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Curl Integrator")
    {
      BilinearForm a_test(nd_fespace, rt_fespace);
      mfem::MixedBilinearForm a_ref(&nd_fespace.Get(), &rt_fespace.Get());
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Weak Curl Integrator")
    {
      BilinearForm a_test(rt_fespace, nd_fespace);
      mfem::MixedBilinearForm a_ref(&rt_fespace.Get(), &nd_fespace.Get());
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                                (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref, -1.0);
    }
    SECTION("Mixed Vector Curl Integrator (H(curl) range)")
    {
      BilinearForm a_test(nd_fespace, nd_fespace);
      mfem::MixedBilinearForm a_ref(&nd_fespace.Get(), &nd_fespace.Get());
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorCurlIntegrator, mfem::MixedVectorCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Weak Curl Integrator (H(curl) domain)")
    {
      BilinearForm a_test(nd_fespace, nd_fespace);
      mfem::MixedBilinearForm a_ref(&nd_fespace.Get(), &nd_fespace.Get());
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(bdr_integ, a_test, a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                                (mfem::Coefficient &)Q_ref);
            break;
          case CoeffType::Matrix:
            AddIntegrators<MixedVectorWeakCurlIntegrator,
                           mfem::MixedVectorWeakCurlIntegrator>(
                bdr_integ, a_test, a_ref, Q, (mfem::MatrixCoefficient &)Q_ref);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref, -1.0);
    }
  }

  // Tests on mixed H1-(H1)ᵈ spaces.
  SECTION("Mixed H1-Vector H1 Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), h1d_fespace(mesh, &h1_fec, dim);
    SECTION("Mixed H1 Gradient Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_jac = true;
      fem::DefaultIntegrationOrder::q_order_extra_pk = -1;
      fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
      mesh.ResetCeedObjects();
      BilinearForm a_test(h1_fespace, h1d_fespace);
      mfem::MixedBilinearForm a_ref(&h1_fespace.Get(), &h1d_fespace.Get());
      if (!bdr_integ)  // MFEM's GradientIntegrator only supports square Jacobians
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators<GradientIntegrator, mfem::GradientIntegrator>(bdr_integ, a_test,
                                                                         a_ref);
            break;
          case CoeffType::Scalar:
            AddIntegrators<GradientIntegrator, mfem::GradientIntegrator>(bdr_integ, a_test,
                                                                         a_ref, Q, Q_ref);
            break;
          case CoeffType::Matrix:
            break;  // No support for non-scalar coefficients in MFEM's GradientIntegrator
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }
}

void RunCeedInterpolatorTests(MPI_Comm comm, const std::string &input, int ref_levels,
                              bool amr, int order)
{
  // Load the mesh.
  auto mesh = Initialize(comm, input, ref_levels, amr);
  const int dim = mesh.Dimension();

  // Match MFEM's default integration orders.
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Run the tests.
  std::string section =
      "Mesh: " + input + "\n" + "Refinement levels: " + std::to_string(ref_levels) + "\n" +
      "AMR: " + std::to_string(amr) + "\n" + "Order: " + std::to_string(order) + "\n";
  INFO(section);

  // Linear interpolators for prolongation.
  SECTION("H1 Prolongation")
  {
    mfem::H1_FECollection coarse_h1_fec(order, dim), fine_h1_fec(order + 1, dim);
    FiniteElementSpace coarse_h1_fespace(mesh, &coarse_h1_fec),
        fine_h1_fespace(mesh, &fine_h1_fec);
    DiscreteLinearOperator id_test(coarse_h1_fespace, fine_h1_fespace);
    id_test.AddDomainInterpolator<IdentityInterpolator>();
    mfem::PRefinementTransferOperator id_ref(coarse_h1_fespace, fine_h1_fespace);
    TestCeedOperatorMult(*id_test.PartialAssemble(), id_ref, true);
  }
  SECTION("H(curl) Prolongation")
  {
    mfem::ND_FECollection coarse_nd_fec(order, dim), fine_nd_fec(order + 1, dim);
    FiniteElementSpace coarse_nd_fespace(mesh, &coarse_nd_fec),
        fine_nd_fespace(mesh, &fine_nd_fec);
    DiscreteLinearOperator id_test(coarse_nd_fespace, fine_nd_fespace);
    id_test.AddDomainInterpolator<IdentityInterpolator>();
    mfem::PRefinementTransferOperator id_ref(coarse_nd_fespace, fine_nd_fespace);
    TestCeedOperatorMult(*id_test.PartialAssemble(), id_ref, true);
  }
  SECTION("H(div) Prolongation")
  {
    mfem::RT_FECollection coarse_rt_fec(order - 1, dim), fine_rt_fec(order, dim);
    FiniteElementSpace coarse_rt_fespace(mesh, &coarse_rt_fec),
        fine_rt_fespace(mesh, &fine_rt_fec);
    DiscreteLinearOperator id_test(coarse_rt_fespace, fine_rt_fespace);
    id_test.AddDomainInterpolator<IdentityInterpolator>();
    mfem::PRefinementTransferOperator id_ref(coarse_rt_fespace, fine_rt_fespace);
    TestCeedOperatorMult(*id_test.PartialAssemble(), id_ref, true);
  }

  // Linear interpolators for differentiation.
  SECTION("H1-H(curl) Discrete Gradient")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);
    DiscreteLinearOperator grad_test(h1_fespace, nd_fespace);
    mfem::DiscreteLinearOperator grad_ref(&h1_fespace.Get(), &nd_fespace.Get());
    grad_test.AddDomainInterpolator<GradientInterpolator>();
    grad_ref.AddDomainInterpolator(new mfem::GradientInterpolator());
    TestCeedOperator(grad_test, grad_ref);
  }
  SECTION("H(curl)-H(div) Discrete Curl")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace nd_fespace(mesh, &nd_fec), rt_fespace(mesh, &rt_fec);
    DiscreteLinearOperator curl_test(nd_fespace, rt_fespace);
    mfem::DiscreteLinearOperator curl_ref(&nd_fespace.Get(), &rt_fespace.Get());
    if (dim == 3)
    {
      curl_test.AddDomainInterpolator<CurlInterpolator>();
      curl_ref.AddDomainInterpolator(new mfem::CurlInterpolator());
    }
    TestCeedOperator(curl_test, curl_ref);
  }
}

void RunCeedBenchmarks(MPI_Comm comm, const std::string &input, int ref_levels, bool amr,
                       int order)
{
  // Load the mesh.
  auto mesh = Initialize(comm, input, ref_levels, amr);
  const int dim = mesh.Dimension();

  // Match MFEM's default integration orders.
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = false;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Run the benchmarks.
  std::string section =
      "Mesh: " + input + "\n" + "Refinement levels: " + std::to_string(ref_levels) + "\n" +
      "AMR: " + std::to_string(amr) + "\n" + "Order: " + std::to_string(order) + "\n";
  INFO(section);
  if (Mpi::Root(comm))
  {
    auto pos = input.find_last_of('/');
    WARN("benchmark input mesh: " << input.substr(pos + 1) << "\n");
  }

  // Initialize coefficients.
  auto Q = BuildCoefficient(mesh, false, CoeffType::Scalar);
  auto MQ = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto Q_ref = BuildCoefficientRef(mesh, false, CoeffType::Scalar);
  auto MQ_ref = BuildCoefficientRef(mesh, false, CoeffType::Matrix);

  // Diffusion + mass benchmark.
  SECTION("Diffusion + Mass Integrator Benchmark")
  {
    auto AssembleTest = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test(fespace);
      a_test.AddDomainIntegrator<DiffusionMassIntegrator>(MQ, Q);
      if (bdr_integ)
      {
        a_test.AddBoundaryIntegrator<MassIntegrator>();
      }
      if (benchmark_assemble_q_data)
      {
        a_test.AssembleQuadratureData();
      }
      return a_test.PartialAssemble();
    };
    auto AssembleTestRef = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test_ref(fespace);
      a_test_ref.AddDomainIntegrator<DiffusionIntegrator>(MQ);
      a_test_ref.AddDomainIntegrator<MassIntegrator>(Q);
      if (bdr_integ)
      {
        a_test_ref.AddBoundaryIntegrator<MassIntegrator>();
      }
      return a_test_ref.PartialAssemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace.Get());
      a_ref->AddDomainIntegrator(
          new mfem::DiffusionIntegrator((mfem::MatrixCoefficient &)MQ_ref));
      a_ref->AddDomainIntegrator(new mfem::MassIntegrator(Q_ref));
      if (bdr_integ)
      {
        a_ref->AddBoundaryIntegrator(new mfem::MassIntegrator());
      }
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec);
    if (Mpi::Root(comm))
    {
      BenchmarkCeedIntegrator(h1_fespace, AssembleTest, AssembleTestRef, AssembleRef,
                              (dim * (dim + 1)) / 2 + 1);
    }
  }

  // Curl-curl + mass benchmark.
  SECTION("Curl-Curl + Mass Integrator Benchmark")
  {
    auto AssembleTest = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test(fespace);
      a_test.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
      if (bdr_integ)
      {
        a_test.AddBoundaryIntegrator<VectorFEMassIntegrator>();
      }
      if (benchmark_assemble_q_data)
      {
        a_test.AssembleQuadratureData();
      }
      return a_test.PartialAssemble();
    };
    auto AssembleTestRef = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test_ref(fespace);
      a_test_ref.AddDomainIntegrator<CurlCurlIntegrator>(MQ);
      a_test_ref.AddDomainIntegrator<VectorFEMassIntegrator>(Q);
      if (bdr_integ)
      {
        a_test_ref.AddBoundaryIntegrator<VectorFEMassIntegrator>();
      }
      return a_test_ref.PartialAssemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace.Get());
      a_ref->AddDomainIntegrator(
          new mfem::CurlCurlIntegrator((mfem::MatrixCoefficient &)MQ_ref));
      a_ref->AddDomainIntegrator(
          new mfem::VectorFEMassIntegrator((mfem::Coefficient &)Q_ref));
      if (bdr_integ)
      {
        a_ref->AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator());
      }
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace nd_fespace(mesh, &nd_fec);
    if (Mpi::Root(comm))
    {
      BenchmarkCeedIntegrator(nd_fespace, AssembleTest, AssembleTestRef, AssembleRef,
                              2 * (dim * (dim + 1)) / 2);
    }
  }

  // Div-div + mass benchmark.
  SECTION("Div-Div + Mass Integrator Benchmark")
  {
    auto AssembleTest = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test(fespace);
      a_test.AddDomainIntegrator<DivDivMassIntegrator>(Q, MQ);
      if (benchmark_assemble_q_data)
      {
        a_test.AssembleQuadratureData();
      }
      return a_test.PartialAssemble();
    };
    auto AssembleTestRef = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test_ref(fespace);
      a_test_ref.AddDomainIntegrator<DivDivIntegrator>(Q);
      a_test_ref.AddDomainIntegrator<VectorFEMassIntegrator>(MQ);
      return a_test_ref.PartialAssemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace.Get());
      a_ref->AddDomainIntegrator(new mfem::DivDivIntegrator(Q_ref));
      a_ref->AddDomainIntegrator(
          new mfem::VectorFEMassIntegrator((mfem::MatrixCoefficient &)MQ_ref));
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace rt_fespace(mesh, &rt_fec);
    if (Mpi::Root(comm))
    {
      BenchmarkCeedIntegrator(rt_fespace, AssembleTest, AssembleTestRef, AssembleRef, 2);
    }
  }

  // Discrete gradient benchmark.
  SECTION("Discrete Gradient Benchmark")
  {
    auto AssembleTest =
        [](const FiniteElementSpace &trial_fespace, const FiniteElementSpace &test_fespace)
    {
      DiscreteLinearOperator a_test(trial_fespace, test_fespace);
      a_test.AddDomainInterpolator<GradientInterpolator>();
      return a_test.PartialAssemble();
    };
    auto AssembleRef = [](FiniteElementSpace &trial_fespace,
                          FiniteElementSpace &test_fespace,
                          mfem::AssemblyLevel assembly_level, bool skip_zeros)
    {
      auto a_ref = std::make_unique<mfem::DiscreteLinearOperator>(&trial_fespace.Get(),
                                                                  &test_fespace.Get());
      a_ref->AddDomainInterpolator(new mfem::GradientInterpolator());
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec), nd_fespace(mesh, &nd_fec);
    if (Mpi::Root(comm))
    {
      BenchmarkCeedInterpolator(h1_fespace, nd_fespace, AssembleTest, AssembleRef);
    }
  }

  // Wait before returning.
  Mpi::Barrier(comm);
}

}  // namespace

TEST_CASE("2D libCEED Operators", "[libCEED]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                         amr, order);
}

TEST_CASE("3D libCEED Operators", "[libCEED]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                         amr, order);
}

TEST_CASE("2D libCEED Interpolators", "[libCEED][Interpolator]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                           amr, order);
}

TEST_CASE("3D libCEED Interpolators", "[libCEED][Interpolator]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                           amr, order);
}

TEST_CASE("3D libCEED Benchmarks", "[libCEED][Benchmark]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh");
  RunCeedBenchmarks(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh,
                    benchmark_ref_levels, false, benchmark_order);
}

}  // namespace palace
