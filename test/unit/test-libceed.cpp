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
#include "utils/communication.hpp"

extern int benchmark_ref_levels;
extern int benchmark_order;
extern bool benchmark_no_fa;
extern bool benchmark_no_mfem_pa;

namespace palace
{

namespace
{

enum class CoeffType
{
  Const,
  Scalar,
  Vector,
  Matrix
};

std::string ToString(CoeffType type)
{
  switch (type)
  {
    case CoeffType::Const:
      return "Constant";
    case CoeffType::Scalar:
      return "Scalar";
    case CoeffType::Vector:
      return "Vector";
    case CoeffType::Matrix:
      return "Matrix";
  }
  return "";
}

// Scalar coefficient
double coeff_function(const Vector &x)
{
  return 1.0 + x[0] * x[0];
}

// Vector coefficient
void vector_coeff_function(const Vector &x, Vector &v)
{
  const int dim = x.Size();
  const double w = coeff_function(x);
  switch (dim)
  {
    case 1:
      v(0) = w;
      break;
    case 2:
      v(0) = w * sqrt(2.0 / 3.0);
      v(1) = w * sqrt(1.0 / 3.0);
      break;
    case 3:
      v(0) = w * sqrt(3.0 / 6.0);
      v(1) = w * sqrt(2.0 / 6.0);
      v(2) = w * sqrt(1.0 / 6.0);
      break;
  }
}

// Matrix coefficient
void matrix_coeff_function(const Vector &x, mfem::DenseMatrix &m)
{
  const int dim = x.Size();
  Vector v(dim);
  vector_coeff_function(x, v);
  m.SetSize(dim);
  m = 0.1;
  for (int i = 0; i < dim; i++)
  {
    m(i, i) = 1.0 + v(i);
  }
}

template <typename T1, typename T2, typename T3>
void AddIntegrators(BilinearForm &a_test, T1 &a_ref, std::unique_ptr<T2> &&blfi_test,
                    T3 *blfi_ref, bool bdr_integ)
{
  if (bdr_integ)
  {
    a_test.AddBoundaryIntegrator(std::move(blfi_test));
    a_ref.AddBoundaryIntegrator(blfi_ref);
  }
  else
  {
    a_test.AddDomainIntegrator(std::move(blfi_test));
    a_ref.AddDomainIntegrator(blfi_ref);
  }
}

void TestCeedOperatorMult(const Operator &op_test, const Operator &op_ref,
                          bool test_transpose)
{
  Vector x(op_ref.Width()), y_ref(op_ref.Height()), y_test(op_ref.Height());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  y_test.UseDevice(true);
  {
    x.Randomize(1);

    op_ref.Mult(x, y_ref);
    op_test.Mult(x, y_test);

    y_test -= y_ref;

    // REQUIRE(y_ref.Norml2() > 0.0);
    REQUIRE(y_test.Norml2() < 1.0e-12 * std::max(y_ref.Norml2(), 1.0));
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

    y_t_test -= y_t_ref;

    // REQUIRE(y_t_ref.Norml2() > 0.0);
    REQUIRE(y_t_test.Norml2() < 1.0e-12 * std::max(y_t_ref.Norml2(), 1.0));
  }
}

void TestCeedOperatorFullAssemble(const mfem::SparseMatrix &mat_test,
                                  const mfem::SparseMatrix &mat_ref)
{
  // Ensure host memory is up to date (mfem::Add is missing the device to host copy).
  mat_test.HostReadI();
  mat_test.HostReadJ();
  mat_test.HostReadData();
  mat_ref.HostReadI();
  mat_ref.HostReadJ();
  mat_ref.HostReadData();

  std::unique_ptr<mfem::SparseMatrix> mat_diff(mfem::Add(1.0, mat_test, -1.0, mat_ref));

  // REQUIRE(mat_ref.MaxNorm() > 0.0);
  REQUIRE(mat_diff->MaxNorm() < 1.0e-12 * std::max(mat_ref.MaxNorm(), 1.0));
}

template <typename T1, typename T2>
void TestCeedOperator(T1 &a_test, T2 &a_ref, bool test_transpose, bool skip_zeros)
{
  a_ref.Assemble(skip_zeros);
  a_ref.Finalize(skip_zeros);
  const mfem::SparseMatrix *mat_ref = &a_ref.SpMat();
  const Operator *op_ref = mat_ref;

  // Test operator application.
  const auto op_test = a_test.Assemble();
  TestCeedOperatorMult(*op_test, *op_ref, test_transpose);

  // Test full assembly.
  const auto mat_test = a_test.FullAssemble(*op_test, skip_zeros);
  TestCeedOperatorFullAssemble(*mat_test, *mat_ref);

  // Test diagonal assembly if possible.
  if (&a_test.GetTrialSpace() == &a_test.GetTestSpace())
  {
    Vector d_ref(mat_ref->Height()), d_test(mat_ref->Height());
    d_ref.UseDevice(true);
    d_test.UseDevice(true);

    mat_ref->GetDiag(d_ref);
    op_test->AssembleDiagonal(d_test);

    d_test -= d_ref;

    // Diagonal assembly for high-order Nedelec spaces is only approximate due to face
    // dofs in 3D.
    double rtol = 1.0e-12;
    const auto &trial_fespace = a_test.GetTrialSpace();
    const auto &test_fespace = a_test.GetTestSpace();
    const auto &trial_fec = *trial_fespace.FEColl();
    const auto &test_fec = *test_fespace.FEColl();
    if (trial_fespace.GetParMesh()->Dimension() == 3 &&
        ((dynamic_cast<const mfem::ND_FECollection *>(&trial_fec) &&
          trial_fec.GetOrder() > 1 && !mfem::UsesTensorBasis(trial_fespace)) ||
         (dynamic_cast<const mfem::ND_FECollection *>(&test_fec) &&
          test_fec.GetOrder() > 1 && !mfem::UsesTensorBasis(test_fespace))))
    {
      rtol = 1.0;
    }

    // REQUIRE(d_ref.Norml2() > 0.0);
    REQUIRE(d_test.Norml2() < rtol * std::max(d_ref.Norml2(), 1.0));
  }
}

void TestCeedOperator(BilinearForm &op_test, mfem::BilinearForm &op_ref)
{
  TestCeedOperator(op_test, op_ref, false, false);
}

void TestCeedOperator(BilinearForm &op_test, mfem::MixedBilinearForm &op_ref)
{
  TestCeedOperator(op_test, op_ref, false, false);
}

void TestCeedOperator(DiscreteLinearOperator &op_test, mfem::DiscreteLinearOperator &op_ref)
{
  TestCeedOperator(op_test, op_ref, true, true);
}

template <typename T1, typename T2, typename T3>
void BenchmarkCeedIntegrator(FiniteElementSpace &fespace, T1 AssembleTest,
                             T2 AssembleTestRef, T3 AssembleRef, int qdata_size)
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
    const auto op_test = AssembleTest(fespace, true);
    const auto op_test_ref = AssembleTestRef(fespace, true);
    const auto mat_test = BilinearForm::FullAssemble(*op_test, skip_zeros);
    const auto mat_test_ref = BilinearForm::FullAssemble(*op_test_ref, skip_zeros);
    nnz = mat_test->NumNonZeroElems();
    TestCeedOperatorFullAssemble(*mat_test, *mat_test_ref);
  }

  // Benchmark MFEM legacy assembly.
  if (!benchmark_no_fa)
  {
    BENCHMARK("Assemble (MFEM Legacy)")
    {
      const auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      return op_ref->Height();
    };
    {
      const auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Legacy)")
      {
        op_ref->AddMult(x, y_ref);
      };
    }
  }

  // Benchmark MFEM PA (tensor-product elements only).
  if (!benchmark_no_mfem_pa && mfem::UsesTensorBasis(fespace))
  {
    BENCHMARK("Assemble (MFEM Partial)")
    {
      const auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::PARTIAL, skip_zeros);
      return op_ref->Height();
    };
    {
      const auto op_ref = AssembleRef(fespace, mfem::AssemblyLevel::PARTIAL, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Partial)")
      {
        // MFEM PA does not implement AddMult from BilinearForm.
        op_ref->Mult(x, y_test);
        y_ref += y_test;
      };
    }
  }

  // Benchmark libCEED assembly.
  BENCHMARK("Assemble (libCEED)")
  {
    const auto op_test = AssembleTest(fespace);
    return op_test->Height();
  };
  {
    const auto op_test = AssembleTest(fespace);
    y_test = 0.0;
    BENCHMARK("AddMult (libCEED)")
    {
      op_test->AddMult(x, y_test);
    };
  }
  if (!benchmark_no_fa)
  {
    BENCHMARK("Full Assemble (libCEED)")
    {
      const auto op_test = AssembleTest(fespace);
      const auto mat_test = BilinearForm::FullAssemble(*op_test, skip_zeros);
      return mat_test->NumNonZeroElems();
    };
  }

  // Memory estimate (only for non-mixed meshes).
  mfem::ParMesh &mesh = *fespace.GetParMesh();
  if (mesh.GetNumGeometries(mesh.Dimension()) == 1)
  {
    // Integration rule gives the complete non-tensor number of points.
    const mfem::FiniteElement &fe = *fespace.GetFE(0);
    const mfem::ElementTransformation &T = *mesh.GetElementTransformation(0);
    const int q_order = fem::GetDefaultIntegrationOrder(fe, fe, T);
    const int Q = mfem::IntRules.Get(mesh.GetElementGeometry(0), q_order).GetNPoints();
    const int P = fe.GetDof();

    // Rough estimate for memory consumption as quadrature data + offsets for element
    // restriction.
    std::size_t mem_ref = nnz * (8 + 4) + (y_ref.Size() + 1) * 4;
    std::size_t mem_test = (Q * qdata_size * 8 + P * 4) * (std::size_t)mesh.GetNE();
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
  const bool skip_zeros = true;
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
    const auto op_test = AssembleTest(trial_fespace, test_fespace);
    const auto op_ref =
        AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
    const auto mat_test = DiscreteLinearOperator::FullAssemble(*op_test, skip_zeros);
    const auto *mat_ref = &op_ref->SpMat();
    nnz = mat_test->NumNonZeroElems();
    TestCeedOperatorFullAssemble(*mat_test, *mat_ref);
  }

  // Benchmark MFEM legacy assembly.
  if (!benchmark_no_fa)
  {
    BENCHMARK("Assemble (MFEM Legacy)")
    {
      const auto op_ref =
          AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      return op_ref->Height();
    };
    {
      const auto op_ref =
          AssembleRef(trial_fespace, test_fespace, mfem::AssemblyLevel::LEGACY, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Legacy)")
      {
        op_ref->AddMult(x, y_ref);
      };
    }
  }

  // Benchmark MFEM PA (tensor-product elements only).
  if (!benchmark_no_mfem_pa && mfem::UsesTensorBasis(trial_fespace) &&
      mfem::UsesTensorBasis(test_fespace))
  {
    BENCHMARK("Assemble (MFEM Partial)")
    {
      const auto op_ref = AssembleRef(trial_fespace, test_fespace,
                                      mfem::AssemblyLevel::PARTIAL, skip_zeros);
      return op_ref->Height();
    };
    {
      const auto op_ref = AssembleRef(trial_fespace, test_fespace,
                                      mfem::AssemblyLevel::PARTIAL, skip_zeros);
      y_ref = 0.0;
      BENCHMARK("AddMult (MFEM Partial)")
      {
        // MFEM PA does not implement AddMult from BilinearForm.
        op_ref->Mult(x, y_test);
        y_ref += y_test;
      };
    }
  }

  // Benchmark libCEED assembly.
  BENCHMARK("Assemble (libCEED)")
  {
    const auto op_test = AssembleTest(trial_fespace, test_fespace);
    return op_test->Height();
  };
  {
    const auto op_test = AssembleTest(trial_fespace, test_fespace);
    y_test = 0.0;
    BENCHMARK("AddMult (libCEED)")
    {
      op_test->AddMult(x, y_test);
    };
  }
  if (!benchmark_no_fa)
  {
    BENCHMARK("Full Assemble (libCEED)")
    {
      const auto op_test = AssembleTest(trial_fespace, test_fespace);
      const auto mat_test = DiscreteLinearOperator::FullAssemble(*op_test, skip_zeros);
      return mat_test->NumNonZeroElems();
    };
  }

  // Memory estimate (only for non-mixed meshes).
  mfem::ParMesh &mesh = *trial_fespace.GetParMesh();
  if (mesh.GetNumGeometries(mesh.Dimension()) == 1)
  {
    const mfem::FiniteElement &trial_fe = *trial_fespace.GetFE(0);
    const mfem::FiniteElement &test_fe = *test_fespace.GetFE(0);
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
                            int order)
{
  // Load the mesh.
  std::unique_ptr<mfem::ParMesh> mesh;
  {
    mfem::Mesh smesh(input, 1, 1);
    smesh.EnsureNodes();
    REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
    mesh = std::make_unique<mfem::ParMesh>(comm, smesh);
    for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }
  }
  const int dim = mesh->Dimension();

  // Initialize coefficients.
  mfem::FunctionCoefficient Q(coeff_function);
  mfem::VectorFunctionCoefficient VQ(dim, vector_coeff_function);
  mfem::MatrixFunctionCoefficient MQ(dim, matrix_coeff_function);

  // Run the tests.
  auto coeff_type =
      GENERATE(CoeffType::Const, CoeffType::Scalar, CoeffType::Vector, CoeffType::Matrix);
  auto bdr_integ = GENERATE(false, true);
  std::string section =
      "Mesh: " + input + "\n" + "Refinement levels: " + std::to_string(ref_levels) + "\n" +
      "Order: " + std::to_string(order) + "\n" + "Coefficient: " + ToString(coeff_type) +
      "\n" + "Integrator: " + (bdr_integ ? "Boundary" : "Domain") + "\n";
  INFO(section);

  // Used in some hacks to match MFEM's default integration order on mixed meshes.
  const int mesh_order = mesh->GetNodalFESpace()->GetMaxElementOrder();
  const int order_w_pk = (mesh_order - 1) * (dim - bdr_integ);
  const int order_w_qk = mesh_order * (dim - bdr_integ) - 1;

  // Tests on H1 spaces.
  SECTION("H1 Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh.get(), &h1_fec),
        vector_h1_fespace(mesh.get(), &h1_fec, dim);
    SECTION("H1 Mass Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(h1_fespace, q_extra);
      mfem::BilinearForm a_ref(&h1_fespace);
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(),
                         new mfem::MassIntegrator(), bdr_integ);
          break;
        case CoeffType::Scalar:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(Q),
                         new mfem::MassIntegrator(Q), bdr_integ);
          break;
        case CoeffType::Vector:
        case CoeffType::Matrix:
          break;  // Good to test empty operators
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Vector H1 Mass Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(vector_h1_fespace, q_extra);
      mfem::BilinearForm a_ref(&vector_h1_fespace);
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(),
                         new mfem::VectorMassIntegrator(), bdr_integ);
          break;
        case CoeffType::Scalar:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(Q),
                         new mfem::VectorMassIntegrator(Q), bdr_integ);
          break;
        case CoeffType::Vector:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(VQ),
                         new mfem::VectorMassIntegrator(VQ), bdr_integ);
          break;
        case CoeffType::Matrix:
          AddIntegrators(a_test, a_ref, std::make_unique<MassIntegrator>(MQ),
                         new mfem::VectorMassIntegrator(MQ), bdr_integ);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("H1 Diffusion Integrator")
    {
      const auto q_extra_pk = -2 - order_w_pk,
                 q_extra_qk = dim - bdr_integ - 1 - order_w_qk;
      BilinearForm a_test(h1_fespace, q_extra_pk, q_extra_qk);
      mfem::BilinearForm a_ref(&h1_fespace);
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators(a_test, a_ref, std::make_unique<DiffusionIntegrator>(),
                         new mfem::DiffusionIntegrator(), bdr_integ);
          break;
        case CoeffType::Scalar:
          AddIntegrators(a_test, a_ref, std::make_unique<DiffusionIntegrator>(Q),
                         new mfem::DiffusionIntegrator(Q), bdr_integ);
          break;
        case CoeffType::Vector:
          AddIntegrators(a_test, a_ref, std::make_unique<DiffusionIntegrator>(VQ),
                         new mfem::DiffusionIntegrator(VQ), bdr_integ);
          break;
        case CoeffType::Matrix:
          AddIntegrators(a_test, a_ref, std::make_unique<DiffusionIntegrator>(MQ),
                         new mfem::DiffusionIntegrator(MQ), bdr_integ);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on H(curl) spaces.
  SECTION("H(curl) Integrators")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace nd_fespace(mesh.get(), &nd_fec);
    SECTION("ND Mass Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(nd_fespace, q_extra);
      mfem::BilinearForm a_ref(&nd_fespace);
      switch (coeff_type)
      {
        case CoeffType::Const:
          AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(),
                         new mfem::VectorFEMassIntegrator(), bdr_integ);
          break;
        case CoeffType::Scalar:
          AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(Q),
                         new mfem::VectorFEMassIntegrator(Q), bdr_integ);
          break;
        case CoeffType::Vector:
          AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(VQ),
                         new mfem::VectorFEMassIntegrator(VQ), bdr_integ);
          break;
        case CoeffType::Matrix:
          AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(MQ),
                         new mfem::VectorFEMassIntegrator(MQ), bdr_integ);
          break;
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("ND Curl-Curl Integrator")
    {
      const auto q_extra_pk = -2 - order_w_pk, q_extra_qk = -order_w_qk;
      BilinearForm a_test(nd_fespace, q_extra_pk, q_extra_qk);
      mfem::BilinearForm a_ref(&nd_fespace);
      if (dim == 3 || (dim == 2 && !bdr_integ))  // No 1D ND curl shape
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<CurlCurlIntegrator>(),
                           new mfem::CurlCurlIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref, std::make_unique<CurlCurlIntegrator>(Q),
                           new mfem::CurlCurlIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            if (dim == 3 && !bdr_integ)
            {
              AddIntegrators(a_test, a_ref, std::make_unique<CurlCurlIntegrator>(VQ),
                             new mfem::CurlCurlIntegrator(VQ), bdr_integ);
            }
            break;
          case CoeffType::Matrix:
            if (dim == 3 && !bdr_integ)
            {
              AddIntegrators(a_test, a_ref, std::make_unique<CurlCurlIntegrator>(MQ),
                             new mfem::CurlCurlIntegrator(MQ), bdr_integ);
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
    FiniteElementSpace rt_fespace(mesh.get(), &rt_fec);
    SECTION("RT Mass Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(rt_fespace, q_extra);
      mfem::BilinearForm a_ref(&rt_fespace);
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(),
                           new mfem::VectorFEMassIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(Q),
                           new mfem::VectorFEMassIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(VQ),
                           new mfem::VectorFEMassIntegrator(VQ), bdr_integ);
            break;
          case CoeffType::Matrix:
            AddIntegrators(a_test, a_ref, std::make_unique<VectorFEMassIntegrator>(MQ),
                           new mfem::VectorFEMassIntegrator(MQ), bdr_integ);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("RT Div-Div Integrator")
    {
      const auto q_extra_pk = -2 - order_w_pk, q_extra_qk = -2 - order_w_qk;
      BilinearForm a_test(rt_fespace, q_extra_pk, q_extra_qk);
      mfem::BilinearForm a_ref(&rt_fespace);
      if (!bdr_integ)  // Boundary RT elements in 2D and 3D are actually L2
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<DivDivIntegrator>(),
                           new mfem::DivDivIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref, std::make_unique<DivDivIntegrator>(Q),
                           new mfem::DivDivIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
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
    FiniteElementSpace h1_fespace(mesh.get(), &h1_fec), nd_fespace(mesh.get(), &nd_fec);
    SECTION("Mixed Vector Gradient Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(h1_fespace, nd_fespace, q_extra);
      mfem::MixedBilinearForm a_ref(&h1_fespace, &nd_fespace);
      if (dim == 3 || (dim == 2 && !bdr_integ))  // Only in 2D or 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorGradientIntegrator>(),
                           new mfem::MixedVectorGradientIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorGradientIntegrator>(Q),
                           new mfem::MixedVectorGradientIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorGradientIntegrator>(VQ),
                           new mfem::MixedVectorGradientIntegrator(VQ), bdr_integ);
            break;
          case CoeffType::Matrix:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorGradientIntegrator>(MQ),
                           new mfem::MixedVectorGradientIntegrator(MQ), bdr_integ);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Weak Divergence Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(nd_fespace, h1_fespace, q_extra);
      mfem::MixedBilinearForm a_ref(&nd_fespace, &h1_fespace);
      if (dim == 3 || (dim == 2 && !bdr_integ))  // Only in 2D or 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakDivergenceIntegrator>(),
                           new mfem::MixedVectorWeakDivergenceIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakDivergenceIntegrator>(Q),
                           new mfem::MixedVectorWeakDivergenceIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakDivergenceIntegrator>(VQ),
                           new mfem::MixedVectorWeakDivergenceIntegrator(VQ), bdr_integ);
            break;
          case CoeffType::Matrix:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakDivergenceIntegrator>(MQ),
                           new mfem::MixedVectorWeakDivergenceIntegrator(MQ), bdr_integ);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }

  // Tests on mixed H(curl)-H(div) spaces.
  SECTION("H1(curl)-H(div) Mixed Integrators")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace nd_fespace(mesh.get(), &nd_fec), rt_fespace(mesh.get(), &rt_fec);
    SECTION("Mixed Vector Curl Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(nd_fespace, rt_fespace, q_extra);
      mfem::MixedBilinearForm a_ref(&nd_fespace, &rt_fespace);
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorCurlIntegrator>(),
                           new mfem::MixedVectorCurlIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorCurlIntegrator>(Q),
                           new mfem::MixedVectorCurlIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorCurlIntegrator>(VQ),
                           new mfem::MixedVectorCurlIntegrator(VQ), bdr_integ);
            break;
          case CoeffType::Matrix:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorCurlIntegrator>(MQ),
                           new mfem::MixedVectorCurlIntegrator(MQ), bdr_integ);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
    SECTION("Mixed Vector Weak Curl Integrator")
    {
      const auto q_extra = 0;
      BilinearForm a_test(rt_fespace, nd_fespace, q_extra);
      mfem::MixedBilinearForm a_ref(&rt_fespace, &nd_fespace);
      if (dim == 3 && !bdr_integ)  // Only in 3D
      {
        switch (coeff_type)
        {
          case CoeffType::Const:
            AddIntegrators(a_test, a_ref, std::make_unique<MixedVectorWeakCurlIntegrator>(),
                           new mfem::MixedVectorWeakCurlIntegrator(), bdr_integ);
            break;
          case CoeffType::Scalar:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakCurlIntegrator>(Q),
                           new mfem::MixedVectorWeakCurlIntegrator(Q), bdr_integ);
            break;
          case CoeffType::Vector:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakCurlIntegrator>(VQ),
                           new mfem::MixedVectorWeakCurlIntegrator(VQ), bdr_integ);
            break;
          case CoeffType::Matrix:
            AddIntegrators(a_test, a_ref,
                           std::make_unique<MixedVectorWeakCurlIntegrator>(MQ),
                           new mfem::MixedVectorWeakCurlIntegrator(MQ), bdr_integ);
            break;
        }
      }
      TestCeedOperator(a_test, a_ref);
    }
  }
}

void RunCeedInterpolatorTests(MPI_Comm comm, const std::string &input, int ref_levels,
                              int order)
{
  // Load the mesh.
  std::unique_ptr<mfem::ParMesh> mesh;
  {
    mfem::Mesh smesh(input, 1, 1);
    smesh.EnsureNodes();
    REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
    mesh = std::make_unique<mfem::ParMesh>(comm, smesh);
    for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }
  }
  const int dim = mesh->Dimension();

  // Run the tests.
  std::string section = "Mesh: " + input + "\n" +
                        "Refinement levels: " + std::to_string(ref_levels) + "\n" +
                        "Order: " + std::to_string(order) + "\n";
  INFO(section);

  // Linear interpolators for prolongation.
  SECTION("H1 Prolongation")
  {
    mfem::H1_FECollection coarse_h1_fec(order, dim), fine_h1_fec(order + 1, dim);
    FiniteElementSpace coarse_h1_fespace(mesh.get(), &coarse_h1_fec),
        fine_h1_fespace(mesh.get(), &fine_h1_fec);
    DiscreteLinearOperator id_test(coarse_h1_fespace, fine_h1_fespace);
    id_test.AddDomainInterpolator(std::make_unique<IdentityInterpolator>());
    mfem::PRefinementTransferOperator id_ref(coarse_h1_fespace, fine_h1_fespace);
    TestCeedOperatorMult(*id_test.Assemble(), id_ref, true);
  }
  SECTION("H(curl) Prolongation")
  {
    mfem::ND_FECollection coarse_nd_fec(order, dim), fine_nd_fec(order + 1, dim);
    FiniteElementSpace coarse_nd_fespace(mesh.get(), &coarse_nd_fec),
        fine_nd_fespace(mesh.get(), &fine_nd_fec);
    DiscreteLinearOperator id_test(coarse_nd_fespace, fine_nd_fespace);
    id_test.AddDomainInterpolator(std::make_unique<IdentityInterpolator>());
    mfem::PRefinementTransferOperator id_ref(coarse_nd_fespace, fine_nd_fespace);
    TestCeedOperatorMult(*id_test.Assemble(), id_ref, true);
  }
  SECTION("H(div) Prolongation")
  {
    mfem::RT_FECollection coarse_rt_fec(order - 1, dim), fine_rt_fec(order, dim);
    FiniteElementSpace coarse_rt_fespace(mesh.get(), &coarse_rt_fec),
        fine_rt_fespace(mesh.get(), &fine_rt_fec);
    DiscreteLinearOperator id_test(coarse_rt_fespace, fine_rt_fespace);
    id_test.AddDomainInterpolator(std::make_unique<IdentityInterpolator>());
    mfem::PRefinementTransferOperator id_ref(coarse_rt_fespace, fine_rt_fespace);
    TestCeedOperatorMult(*id_test.Assemble(), id_ref, true);
  }

  // Linear interpolators for differentiation.
  SECTION("H1-H(curl) Discrete Gradient")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh.get(), &h1_fec), nd_fespace(mesh.get(), &nd_fec);
    DiscreteLinearOperator grad_test(h1_fespace, nd_fespace);
    mfem::DiscreteLinearOperator grad_ref(&h1_fespace, &nd_fespace);
    grad_test.AddDomainInterpolator(std::make_unique<GradientInterpolator>());
    grad_ref.AddDomainInterpolator(new mfem::GradientInterpolator());
    TestCeedOperator(grad_test, grad_ref);
  }
  SECTION("H(curl)-H(div) Discrete Curl")
  {
    mfem::ND_FECollection nd_fec(order, dim);
    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace nd_fespace(mesh.get(), &nd_fec), rt_fespace(mesh.get(), &rt_fec);
    DiscreteLinearOperator curl_test(nd_fespace, rt_fespace);
    mfem::DiscreteLinearOperator curl_ref(&nd_fespace, &rt_fespace);
    if (dim == 3)
    {
      curl_test.AddDomainInterpolator(std::make_unique<CurlInterpolator>());
      curl_ref.AddDomainInterpolator(new mfem::CurlInterpolator());
    }
    TestCeedOperator(curl_test, curl_ref);
  }
}

void RunCeedBenchmarks(MPI_Comm comm, const std::string &input, int ref_levels, int order)
{
  // Load the mesh.
  std::unique_ptr<mfem::ParMesh> mesh;
  {
    mfem::Mesh smesh(input, 1, 1);
    smesh.EnsureNodes();
    REQUIRE(Mpi::Size(comm) <= smesh.GetNE());
    mesh = std::make_unique<mfem::ParMesh>(comm, smesh);
    for (int l = 0; l < ref_levels; l++)
    {
      mesh->UniformRefinement();
    }
  }
  const int dim = mesh->Dimension();

  // Initialize coefficients.
  mfem::FunctionCoefficient Q(coeff_function);
  mfem::MatrixFunctionCoefficient MQ(dim, matrix_coeff_function);

  // Run the benchmarks.
  std::string section = "Mesh: " + input + "\n" +
                        "Refinement levels: " + std::to_string(ref_levels) + "\n" +
                        "Order: " + std::to_string(order) + "\n";
  INFO(section);
  auto pos = input.find_last_of('/');
  WARN("benchmark input mesh: " << input.substr(pos + 1) << "\n");

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
      return a_test.Assemble();
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
      return a_test_ref.Assemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace);
      a_ref->AddDomainIntegrator(new mfem::DiffusionIntegrator(MQ));
      a_ref->AddDomainIntegrator(new mfem::MassIntegrator(Q));
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
    FiniteElementSpace h1_fespace(mesh.get(), &h1_fec);
    BenchmarkCeedIntegrator(h1_fespace, AssembleTest, AssembleTestRef, AssembleRef,
                            (dim * (dim + 1)) / 2 + 1);
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
      return a_test.Assemble();
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
      return a_test_ref.Assemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace);
      a_ref->AddDomainIntegrator(new mfem::CurlCurlIntegrator(MQ));
      a_ref->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(Q));
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
    FiniteElementSpace nd_fespace(mesh.get(), &nd_fec);
    BenchmarkCeedIntegrator(nd_fespace, AssembleTest, AssembleTestRef, AssembleRef,
                            2 * (dim * (dim + 1)) / 2);
  }

  // Div-div + mass benchmark.
  SECTION("Div-Div + Mass Integrator Benchmark")
  {
    auto AssembleTest = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test(fespace);
      a_test.AddDomainIntegrator<DivDivMassIntegrator>(Q, Q);
      return a_test.Assemble();
    };
    auto AssembleTestRef = [&](const FiniteElementSpace &fespace, bool bdr_integ = false)
    {
      BilinearForm a_test_ref(fespace);
      a_test_ref.AddDomainIntegrator<DivDivIntegrator>(Q);
      a_test_ref.AddDomainIntegrator<VectorFEMassIntegrator>(Q);
      return a_test_ref.Assemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &fespace, mfem::AssemblyLevel assembly_level,
                           bool skip_zeros, bool bdr_integ = false)
    {
      auto a_ref = std::make_unique<mfem::BilinearForm>(&fespace);
      a_ref->AddDomainIntegrator(new mfem::DivDivIntegrator(Q));
      a_ref->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(Q));
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::RT_FECollection rt_fec(order - 1, dim);
    FiniteElementSpace rt_fespace(mesh.get(), &rt_fec);
    BenchmarkCeedIntegrator(rt_fespace, AssembleTest, AssembleTestRef, AssembleRef, 2);
  }

  // Discrete gradient benchmark.
  SECTION("Discrete Gradient Benchmark")
  {
    auto AssembleTest =
        [&](const FiniteElementSpace &trial_fespace, const FiniteElementSpace &test_fespace)
    {
      DiscreteLinearOperator a_test(trial_fespace, test_fespace);
      a_test.AddDomainInterpolator(std::make_unique<GradientInterpolator>());
      return a_test.Assemble();
    };
    auto AssembleRef = [&](FiniteElementSpace &trial_fespace,
                           FiniteElementSpace &test_fespace,
                           mfem::AssemblyLevel assembly_level, bool skip_zeros)
    {
      auto a_ref =
          std::make_unique<mfem::DiscreteLinearOperator>(&trial_fespace, &test_fespace);
      a_ref->AddDomainInterpolator(new mfem::GradientInterpolator());
      a_ref->SetAssemblyLevel(assembly_level);
      a_ref->Assemble(skip_zeros);
      a_ref->Finalize(skip_zeros);
      return a_ref;
    };

    mfem::H1_FECollection h1_fec(order, dim);
    mfem::ND_FECollection nd_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh.get(), &h1_fec), nd_fespace(mesh.get(), &nd_fec);
    BenchmarkCeedInterpolator(h1_fespace, nd_fespace, AssembleTest, AssembleRef);
  }
}

}  // namespace

TEST_CASE("2D libCEED Operators", "[libCEED]")
{
  auto mesh =
      GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh", "star-amr.mesh");
  auto order = GENERATE(1, 2);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                         order);
}

TEST_CASE("3D libCEED Operators", "[libCEED]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh",
                       "fichera-amr.mesh");
  auto order = GENERATE(1, 2);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                         order);
}

TEST_CASE("2D libCEED Interpolators", "[libCEED][Interpolator]")
{
  auto mesh =
      GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh", "star-amr.mesh");
  auto order = GENERATE(1, 2);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                           order);
}

TEST_CASE("3D libCEED Interpolators", "[libCEED][Interpolator]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh",
                       "fichera-amr.mesh");
  auto order = GENERATE(1, 2);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh, 0,
                           order);
}

TEST_CASE("3D libCEED Benchmarks", "[libCEED][Benchmark]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh");
  RunCeedBenchmarks(MPI_COMM_WORLD, std::string(PALACE_TEST_MESH_DIR "/") + mesh,
                    benchmark_ref_levels, benchmark_order);
}

}  // namespace palace
