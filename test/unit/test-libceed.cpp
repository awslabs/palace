// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <typeinfo>
#include <ceed/backend.h>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/interfaces/catch_interfaces_config.hpp>
#include <catch2/internal/catch_context.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/mesh.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/omp.hpp"
#include "utils/units.hpp"

extern int benchmark_ref_levels;
extern int benchmark_order;
extern bool benchmark_assemble_q_data;
extern bool benchmark_no_fa;
extern bool benchmark_no_mfem_pa;

namespace palace
{

namespace
{

auto Initialize(MPI_Comm comm, const std::string &input, int ref_levels, bool amr,
                bool embed_in_3d = false)
{
  // Load the mesh.
  mfem::Mesh smesh(input, 1, 1);
  smesh.EnsureNodes();

  // Optionally embed a 2D mesh in 3D space (SpaceDim=3, Dim=2). This produces boundary
  // elements with SpaceDim=3, Dim=1 that exercise the 31 libCEED qfunctions.
  if (embed_in_3d && smesh.Dimension() == 2)
  {
    const int mesh_order = smesh.GetNodes()->FESpace()->GetMaxElementOrder();
    smesh.SetCurvature(mesh_order, false, 3, mfem::Ordering::byNODES);
    auto *nodes = smesh.GetNodes();
    const int npts = nodes->Size() / 3;
    for (int i = 0; i < npts; i++)
    {
      (*nodes)(2 * npts + i) = 0.0;  // z = 0
    }
  }

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

  // Perform nonconforming AMR (two levels of refinement with no hanging node restrictions).
  if (amr)
  {
    pmesh->RandomRefinement(0.5);
    pmesh->RandomRefinement(0.5);
  }

  auto mesh = Mesh(std::move(pmesh));

  // For 2D-in-3D meshes, rebuild CEED attributes so that boundary element geometry data
  // is constructed. Without this, BuildCeedGeomFactorData skips boundary elements when
  // dim != sdim (the condition is: dim == sdim || ceed_from_self).
  if (embed_in_3d)
  {
    mesh.RebuildCeedAttributes();
  }

  return mesh;
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

// Quadrature data assembly splits an integrator into a build QFunction, which writes the
// per-quadrature-point tensor into a cache, and a generic apply QFunction which consumes
// it. The two must agree on how that cache is laid out: a wrong block offset in a build
// QFunction silently corrupts one block and leaves another zero. Compare against the same
// integrator applied without cached quadrature data, which shares the coefficient and
// geometry code paths but none of the layout logic.
template <typename T>
void TestCeedQuadratureData(MPI_Comm comm, const FiniteElementSpace &fespace,
                            T AddIntegrators)
{
  BilinearForm a_test(fespace), a_ref(fespace);
  AddIntegrators(a_test);
  AddIntegrators(a_ref);
  a_test.AssembleQuadratureData();
  auto op_test = a_test.PartialAssemble();
  auto op_ref = a_ref.PartialAssemble();

  // Guard against vacuously comparing two empty operators.
  Vector x(op_ref->Width()), y_ref(op_ref->Height());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  x.Randomize(1);
  op_ref->Mult(x, y_ref);
  double norm_ref = y_ref * y_ref;
  Mpi::GlobalSum(1, &norm_ref, comm);
  REQUIRE(norm_ref > 0.0);

  TestCeedOperatorMult(*op_test, *op_ref, false);
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
                            bool amr, int order, bool embed_in_3d = false)
{
  // Load the mesh.
  auto mesh = Initialize(comm, input, ref_levels, amr, embed_in_3d);
  const int dim = mesh.Dimension();

  // Match MFEM's default integration orders.
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Run the tests. For 2D-in-3D, skip matrix coefficients (need sdim x sdim dimensions)
  // and test boundary integrators only (domain integrators use the 32 qfunctions which are
  // already covered by the standard 2D tests after projection).
  auto bdr_integ = embed_in_3d ? GENERATE(true) : GENERATE(false, true);
  auto coeff_type = embed_in_3d
                        ? GENERATE(CoeffType::Const, CoeffType::Scalar)
                        : GENERATE(CoeffType::Const, CoeffType::Scalar, CoeffType::Matrix);
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
      // MFEM's GradientIntegrator only supports square Jacobians (dim == sdim).
      if (!bdr_integ && mesh.SpaceDimension() == dim)
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

void CheckMfemFixedBasis(const mfem::FiniteElement &fe, const mfem::IntegrationRule &points,
                         bool check_gradient)
{
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  CeedBasis basis;
  // The arbitrary target rule is tabulated once at setup. Apply then uses the ordinary
  // fixed basis API, exactly as mapped face/subface operators do.
  ceed::InitBasisFromRule(fe, points, 1, ceed, &basis);

  const int num_nodes = fe.GetDof();
  const int num_points = points.GetNPoints();
  const int value_dim = fe.GetRangeType() == mfem::FiniteElement::VECTOR ? fe.GetDim() : 1;
  CeedVector u, v;
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, num_nodes, &u));
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, value_dim * num_points, &v));

  mfem::Vector u_values(num_nodes);
  for (int i = 0; i < num_nodes; i++)
  {
    u_values(i) = 0.25 * (i + 1) - 0.1 * (i % 3);
  }
  PalaceCeedCall(
      ceed, CeedVectorSetArray(u, CEED_MEM_HOST, CEED_COPY_VALUES, u_values.GetData()));
  PalaceCeedCall(ceed, CeedBasisApply(basis, 1, CEED_NOTRANSPOSE, CEED_EVAL_INTERP, u, v));

  const CeedScalar *values;
  PalaceCeedCall(ceed, CeedVectorGetArrayRead(v, CEED_MEM_HOST, &values));
  mfem::Vector shape(num_nodes);
  mfem::DenseMatrix vshape(num_nodes, fe.GetDim());
  for (int q = 0; q < num_points; q++)
  {
    if (value_dim == 1)
    {
      fe.CalcShape(points.IntPoint(q), shape);
      CHECK(values[q] == Catch::Approx(shape * u_values).epsilon(1.0e-11).margin(1.0e-13));
    }
    else
    {
      fe.CalcVShape(points.IntPoint(q), vshape);
      for (int d = 0; d < value_dim; d++)
      {
        mfem::Vector column(vshape.GetColumn(d), num_nodes);
        CHECK(values[d * num_points + q] ==
              Catch::Approx(column * u_values).epsilon(1.0e-11).margin(1.0e-13));
      }
    }
  }
  PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(v, &values));

  if (check_gradient)
  {
    PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, fe.GetDim() * num_points, &v));
    PalaceCeedCall(ceed, CeedBasisApply(basis, 1, CEED_NOTRANSPOSE, CEED_EVAL_GRAD, u, v));
    PalaceCeedCall(ceed, CeedVectorGetArrayRead(v, CEED_MEM_HOST, &values));
    mfem::DenseMatrix dshape(num_nodes, fe.GetDim());
    for (int q = 0; q < num_points; q++)
    {
      fe.CalcDShape(points.IntPoint(q), dshape);
      for (int d = 0; d < fe.GetDim(); d++)
      {
        mfem::Vector column(dshape.GetColumn(d), num_nodes);
        CHECK(values[d * num_points + q] ==
              Catch::Approx(column * u_values).epsilon(1.0e-11).margin(1.0e-13));
      }
    }
    PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(v, &values));
  }

  PalaceCeedCall(ceed, CeedVectorDestroy(&u));
  PalaceCeedCall(ceed, CeedVectorDestroy(&v));
  PalaceCeedCall(ceed, CeedBasisDestroy(&basis));
}

}  // namespace

TEST_CASE("MFEM fixed arbitrary-rule bases", "[libCEED][Serial][Parallel][GPU]")
{
  SECTION("Rational pyramid H1")
  {
    mfem::LinearPyramidFiniteElement fe;
    mfem::IntegrationRule points(3);
    points.IntPoint(0).Set3(0.10, 0.10, 0.50);
    points.IntPoint(1).Set3(0.20, 0.15, 0.30);
    points.IntPoint(2).Set3(0.05, 0.20, 0.60);
    for (int q = 0; q < points.GetNPoints(); q++)
    {
      points.IntPoint(q).weight = 1.0;
    }
    CheckMfemFixedBasis(fe, points, true);
  }

  SECTION("Square full-rank wedge Hcurl and Hdiv")
  {
    for (int order : {1, 2})
    {
      mfem::ND_WedgeElement nd_fe(order);
      mfem::RT_WedgeElement rt_fe(order - 1);
      mfem::IntegrationRule points(3);
      points.IntPoint(0).Set3(0.20, 0.10, 0.25);
      points.IntPoint(1).Set3(0.40, 0.20, 0.75);
      points.IntPoint(2).Set3(0.10, 0.30, 0.50);
      for (int q = 0; q < points.GetNPoints(); q++)
      {
        points.IntPoint(q).weight = 1.0;
      }
      CheckMfemFixedBasis(nd_fe, points, false);
      CheckMfemFixedBasis(rt_fe, points, false);
    }
  }
}

TEST_CASE("2D libCEED Operators", "[libCEED][Serial][Parallel]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh,
                         0, amr, order);
}

TEST_CASE("3D libCEED Operators", "[libCEED][Serial][Parallel]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh,
                         0, amr, order);
}

TEST_CASE("2D libCEED Interpolators", "[libCEED][Interpolator][Serial][Parallel]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(
      MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh, 0, amr, order);
}

TEST_CASE("3D libCEED Interpolators", "[libCEED][Interpolator][Serial][Parallel]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(
      MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh, 0, amr, order);
}

// Test the 31 (SpaceDim=3, Dim=1) libCEED qfunctions for boundary integrators on a 2D
// mesh embedded in 3D. Uses the standard RunCeedIntegratorTests with embed_in_3d=true.
// Matrix coefficients are excluded (need sdim x sdim, not dim x dim).
TEST_CASE("2D-in-3D libCEED Boundary Operators", "[libCEED][Serial][Parallel]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh");
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh,
                         0, false, order, true);
}

// SpaceOperator::AssemblePreconditioner assembles quadrature data for every integrator it
// configures when running on CPU, including the boundary terms contributed by absorbing and
// impedance boundaries. Cover each (SpaceDim, Dim) combination those integrators reach; in
// particular, a boundary curl-curl + mass term on a 3D mesh selects the 32 QFunctions,
// which no other test exercises.
TEST_CASE("libCEED Quadrature Data Assembly", "[libCEED][Serial][Parallel]")
{
  auto mesh_file =
      GENERATE("star-quad.mesh", "star-tri.mesh", "fichera-hex.mesh", "fichera-tet.mesh");
  auto order = GENERATE(1, 2);
  const auto comm = MPI_COMM_WORLD;
  auto mesh =
      Initialize(comm, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh_file, 0, false);
  const int dim = mesh.Dimension();

  // Match MFEM's default integration orders.
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  INFO("Mesh: " << mesh_file << "\nOrder: " << order);

  auto Q = BuildCoefficient(mesh, false, CoeffType::Scalar);
  auto MQ = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto Q_bdr = BuildCoefficient(mesh, true, CoeffType::Scalar);
  auto MQ_bdr = BuildCoefficient(mesh, true, CoeffType::Matrix);

  mfem::ND_FECollection nd_fec(order, dim);
  mfem::H1_FECollection h1_fec(order, dim);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);

  SECTION("Domain Curl-Curl + Mass")
  {
    // The curl coefficient is scalar wherever the curl itself is scalar-valued (Dim < 3).
    TestCeedQuadratureData(comm, nd_fespace,
                           [&](BilinearForm &a)
                           {
                             if (dim < 3)
                             {
                               a.AddDomainIntegrator<CurlCurlMassIntegrator>(Q, MQ);
                             }
                             else
                             {
                               a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
                             }
                           });
  }
  if (dim == 3)
  {
    SECTION("Boundary Curl-Curl + Mass")
    {
      // A second-order absorbing boundary fills both coefficients, so
      // AddConfiguredIntegrators builds this integrator on the boundary.
      TestCeedQuadratureData(
          comm, nd_fespace, [&](BilinearForm &a)
          { a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(Q_bdr, MQ_bdr); });
    }
  }
  SECTION("Boundary Mass")
  {
    TestCeedQuadratureData(comm, nd_fespace, [&](BilinearForm &a)
                           { a.AddBoundaryIntegrator<VectorFEMassIntegrator>(MQ_bdr); });
  }
  SECTION("Auxiliary Diffusion")
  {
    TestCeedQuadratureData(comm, h1_fespace, [&](BilinearForm &a)
                           { a.AddDomainIntegrator<DiffusionIntegrator>(MQ); });
  }
}

namespace
{

struct PackedIntegrationSettings
{
  int threshold = BilinearForm::pa_order_threshold;
  int order = fem::DefaultIntegrationOrder::p_trial;
  bool jac = fem::DefaultIntegrationOrder::q_order_jac;
  int pk = fem::DefaultIntegrationOrder::q_order_extra_pk;
  int qk = fem::DefaultIntegrationOrder::q_order_extra_qk;

  explicit PackedIntegrationSettings(int p)
  {
    fem::DefaultIntegrationOrder::p_trial = p;
    fem::DefaultIntegrationOrder::q_order_jac = true;
    fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
    fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  }
  ~PackedIntegrationSettings()
  {
    BilinearForm::pa_order_threshold = threshold;
    fem::DefaultIntegrationOrder::p_trial = order;
    fem::DefaultIntegrationOrder::q_order_jac = jac;
    fem::DefaultIntegrationOrder::q_order_extra_pk = pk;
    fem::DefaultIntegrationOrder::q_order_extra_qk = qk;
  }
};

bool PackedTestBackend()
{
  if (ceed::internal::NumCeeds() != 1 || utils::GetMaxThreads() > 1 ||
      mfem::Device::Allows(mfem::Backend::DEVICE_MASK))
  {
    return false;
  }
  const char *resource;
  CeedMemType mem;
  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  REQUIRE(CeedGetResource(ceed, &resource) == 0);
  REQUIRE(CeedGetPreferredMemType(ceed, &mem) == 0);
  return mem == CEED_MEM_HOST && std::string(resource).find("/cpu/") == 0;
}

bool IsOriginalComplexWrapper(const ComplexOperator &op)
{
  return typeid(op) == typeid(ComplexWrapperOperator);
}

void CheckPackedResult(const ComplexVector &actual, const ComplexVector &expected)
{
  REQUIRE(actual.Size() == expected.Size());
  constexpr double tol = 5.0e-13;
  for (int part = 0; part < 2; part++)
  {
    const auto *a = (part == 0 ? actual.Real() : actual.Imag()).HostRead();
    const auto *b = (part == 0 ? expected.Real() : expected.Imag()).HostRead();
    double error2 = 0.0, norm2 = 0.0, error_max = 0.0, ref_max = 0.0;
    bool finite = true;
    for (int j = 0; j < actual.Size(); j++)
    {
      finite = finite && std::isfinite(a[j]) && std::isfinite(b[j]);
      const double delta = a[j] - b[j];
      error2 += delta * delta;
      norm2 += b[j] * b[j];
      error_max = std::max(error_max, std::abs(delta));
      ref_max = std::max(ref_max, std::abs(b[j]));
    }
    CAPTURE(part, error2, norm2, error_max, ref_max);
    REQUIRE(finite);
    REQUIRE(error2 <= tol * tol * norm2);
    REQUIRE(error_max <= tol * ref_max);
  }
}

void CheckPackedActions(const ComplexOperator &actual, const ComplexOperator &reference,
                        bool inherited = false)
{
  ComplexVector x(actual.Width()), y(actual.Height()), expected(actual.Height()),
      action(actual.Height()), initial(actual.Height());
  x.Real().Randomize(17);
  x.Imag().Randomize(53);
  initial.Real().Randomize(29);
  initial.Imag().Randomize(97);
  {
    INFO("forward application");
    reference.Mult(x, action);
    actual.Mult(x, y);
    CheckPackedResult(y, action);
  }
  {
    INFO("scaled AddMult application");
    const std::complex<double> a{-0.7, 0.23};
    y = initial;
    expected = initial;
    expected.AXPY(a, action);
    actual.AddMult(x, y, a);
    CheckPackedResult(y, expected);
  }
  if (inherited)
  {
    // These operations retain the original real/imaginary operators. Check them once,
    // independently of the forward-kernel/remainder cases.
    {
      INFO("transpose application");
      reference.MultTranspose(x, expected);
      actual.MultTranspose(x, y);
      CheckPackedResult(y, expected);
    }
    {
      INFO("adjoint application");
      reference.MultHermitianTranspose(x, expected);
      actual.MultHermitianTranspose(x, y);
      CheckPackedResult(y, expected);
    }
    {
      INFO("diagonal assembly");
      reference.AssembleDiagonal(expected);
      actual.AssembleDiagonal(y);
      CheckPackedResult(y, expected);
    }
    for (auto a : {std::complex<double>{1.0}, std::complex<double>{0.0}})
    {
      y = initial;
      expected = initial;
      expected.AXPY(a, action);
      CAPTURE(a);
      actual.AddMult(x, y, a);
      CheckPackedResult(y, expected);
    }
  }
}

void ScalePackedTestQData(const ceed::Operator &op, double scale)
{
  CeedOperator *leaves;
  CeedInt count;
  REQUIRE(CeedOperatorCompositeGetNumSub(op[0], &count) == 0);
  REQUIRE(CeedOperatorCompositeGetSubList(op[0], &leaves) == 0);
  for (CeedInt j = 0; j < count; j++)
  {
    CeedOperatorField field;
    CeedVector qdata = nullptr;
    REQUIRE(CeedOperatorGetFieldByName(leaves[j], "q_data", &field) == 0);
    REQUIRE(CeedOperatorFieldGetVector(field, &qdata) == 0);
    REQUIRE(CeedVectorScale(qdata, scale) == 0);
    REQUIRE(CeedVectorDestroy(&qdata) == 0);
  }
}

struct ComplexPreconditionerFixture
{
  PackedIntegrationSettings settings;
  config::SolverData solver;
  config::DomainData domains;
  config::BoundaryData boundaries;
  Units units{1.0, 1.0};
  std::vector<std::unique_ptr<Mesh>> meshes;

  ComplexPreconditionerFixture(int order, bool amr, int ref_levels = 0) : settings(order)
  {
    // Coaxial-example dielectric, with an unshifted complex fine operator.
    REQUIRE_FALSE(solver.linear.pc_mat_real);
    solver.order = order;
    solver.pa_order_threshold = 2;
    solver.linear.mg_max_levels = order;
    solver.linear.mg_coarsening = MultigridCoarsening::LINEAR;
    solver.linear.pc_mat_shifted = 0;
    BilinearForm::pa_order_threshold = solver.pa_order_threshold;
    fem::DefaultIntegrationOrder::q_order_jac = solver.q_order_jac;

    auto smesh = mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::TETRAHEDRON);
    smesh.EnsureNodes();
    while (smesh.GetNE() < Mpi::Size(Mpi::World()))
    {
      smesh.UniformRefinement();
    }
    for (int l = 0; l < ref_levels; l++)
    {
      smesh.UniformRefinement();
    }
    if (amr)
    {
      smesh.EnsureNCMesh(true);
      mfem::Array<int> refine(1);
      refine[0] = 0;
      smesh.GeneralRefinement(refine);
    }
    meshes.push_back(std::make_unique<Mesh>(Mpi::World(), smesh));
    config::MaterialData material;
    material.attributes = {1};
    material.epsilon_r = 2.08;
    material.tandelta = 4.0e-4;
    domains.attributes = {1};
    domains.materials = {material};
    boundaries.pec.attributes = {1};
    boundaries.farfield.attributes = {2};
  }

  auto MakeSpace()
  {
    return std::make_unique<SpaceOperator>(solver, domains, boundaries, ProblemType::DRIVEN,
                                           units, meshes);
  }
};

auto AssembleComplexPreconditioner(SpaceOperator &space)
{
  constexpr double omega = 2.0;
  return space.GetPreconditionerMatrix<ComplexOperator>(
      std::complex<double>{1.0}, std::complex<double>{0.0, omega},
      std::complex<double>{-omega * omega}, omega);
}

}  // namespace

TEST_CASE("libCEED packed complex QData application",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Packed application requires one CPU CEED context and one thread");
  }
  const auto [name, mesh_file, order, curl, boundary, nested] =
      GENERATE(table<const char *, const char *, int, bool, bool, bool>(
          {{"Mass/mass and shared QData", "fichera-tet.mesh", 1, false, false, false},
           {"Curlmass/mass and p3 face orientations", "fichera-tet.mesh", 3, true, false,
            false},
           {"Hexahedron and unmatched boundary terms", "fichera-hex.mesh", 3, true, true,
            false},
           {"Direct volume pair and nested remainder", "fichera-tet.mesh", 2, true, false,
            true}}));
  CAPTURE(name);
  const bool inherited = curl && !boundary && !nested;
  PackedIntegrationSettings settings(order);
  auto mesh = Initialize(Mpi::World(),
                         std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh_file, 0, false);
  mfem::ND_FECollection fec(order, 3);
  FiniteElementSpace fespace(mesh, &fec);
  auto real_mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto imag_mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto real_curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto bdr_mass = BuildCoefficient(mesh, true, CoeffType::Matrix);
  auto bdr_curl = BuildCoefficient(mesh, true, CoeffType::Scalar);
  real_mass *= -0.71;
  imag_mass *= 0.19;
  real_curl *= 1.37;
  bdr_mass *= -0.23;
  BilinearForm ar(fespace), ai(fespace);
  if (curl)
  {
    ar.AddDomainIntegrator<CurlCurlMassIntegrator>(real_curl, real_mass);
  }
  else
  {
    ar.AddDomainIntegrator<VectorFEMassIntegrator>(real_mass);
  }
  ai.AddDomainIntegrator<VectorFEMassIntegrator>(imag_mass);
  if (boundary)
  {
    // Deliberately unmatched real/imaginary boundary field layouts.
    ar.AddBoundaryIntegrator<VectorFEMassIntegrator>(bdr_mass);
    ai.AddBoundaryIntegrator<CurlCurlMassIntegrator>(bdr_curl, bdr_mass);
  }
  ar.AssembleQuadratureData();
  ai.AssembleQuadratureData();
  auto real = ar.PartialAssemble();
  auto imag = ai.PartialAssemble();
  auto *original_real = real.get();
  auto *original_imag = imag.get();
  auto ref_real = ar.PartialAssemble(), ref_imag = ai.PartialAssemble();
  if (nested)
  {
    // One small nested composite suffices: it must survive as a remainder alongside
    // the compatible direct volume pair. Finish construction before transferring ownership.
    auto AddNestedRemainder = [&](ceed::Operator &op)
    {
      auto remainder = ai.PartialAssemble();
      CeedOperator sub = nullptr;
      REQUIRE(CeedOperatorReferenceCopy((*remainder)[0], &sub) == 0);
      op.AddSubOperator(sub);
      op.Finalize();
    };
    AddNestedRemainder(*real);
    AddNestedRemainder(*ref_real);
  }
  ComplexWrapperOperator reference(ref_real.get(), ref_imag.get());
  auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag));
  // Establish activation without exposing the private implementation or its decomposition.
  REQUIRE_FALSE(IsOriginalComplexWrapper(*actual));
  if (inherited)
  {
    // Verify that this test really exercises ND face transformations, not just signs.
    Ceed ceed = ceed::internal::GetCeedObjects()[0];
    const auto &geom = mesh.GetCeedGeomFactorData(ceed).at(mfem::Geometry::TETRAHEDRON);
    CeedRestrictionType type;
    REQUIRE(CeedElemRestrictionGetType(fespace.GetCeedElemRestriction(
                                           ceed, mfem::Geometry::TETRAHEDRON, geom.indices),
                                       &type) == 0);
    REQUIRE(type == CEED_RESTRICTION_CURL_ORIENTED);
  }
  CheckPackedActions(*actual, reference, inherited);
  if (!curl)
  {
    // Only passive values change; the owned inputs remain structurally immutable.
    ScalePackedTestQData(*original_real, -1.13);
    ScalePackedTestQData(*original_imag, 0.79);
    ScalePackedTestQData(*ref_real, -1.13);
    ScalePackedTestQData(*ref_imag, 0.79);
    CheckPackedActions(*actual, reference);
  }
}

TEST_CASE("Ordinary driven fine preconditioner activates packed complex QData",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Packed application requires one CPU CEED context and one thread");
  }
  const auto [name, amr, lossy, real_pc] = GENERATE(table<const char *, bool, bool, bool>(
      {{"Conforming lossy complex preconditioner", false, true, false},
       {"Nonconforming lossy complex preconditioner", true, true, false},
       {"Conforming lossless volume with absorbing boundary", false, false, false},
       {"Conforming real-only preconditioner", false, true, true}}));
  CAPTURE(name);
  ComplexPreconditionerFixture fixture(3, amr);
  fixture.domains.materials[0].tandelta = lossy ? 4.0e-4 : 0.0;
  fixture.solver.linear.pc_mat_real = real_pc;
  auto space = fixture.MakeSpace();
  auto pc = AssembleComplexPreconditioner(*space);
  const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(pc.get());
  REQUIRE(mg);
  REQUIRE(mg->GetNumLevels() == 3);
  const auto *fine = dynamic_cast<const ComplexParOperator *>(&mg->GetFinestOperator());
  const auto *coarse = dynamic_cast<const ComplexParOperator *>(&mg->GetOperatorAtLevel(0));
  REQUIRE(fine);
  REQUIRE(coarse);
  REQUIRE(dynamic_cast<const hypre::HypreCSRMatrix *>(coarse->LocalOperator().Real()));
  REQUIRE(IsOriginalComplexWrapper(coarse->LocalOperator()));
  const auto &aux =
      dynamic_cast<const ComplexParOperator &>(mg->GetFinestAuxiliaryOperator());
  REQUIRE(IsOriginalComplexWrapper(aux.LocalOperator()));
  const auto &local = fine->LocalOperator();
  REQUIRE(!IsOriginalComplexWrapper(local) == (lossy && !real_pc));
  // Assemble independent reference operators and exercise the unchanged borrowed RAP
  // constructor. Numerical agreement is checked against separate real ParOperators below.
  auto ref_pc = AssembleComplexPreconditioner(*space);
  const auto &ref_mg = dynamic_cast<const ComplexMultigridOperator &>(*ref_pc);
  const auto &ref_local =
      dynamic_cast<const ComplexParOperator &>(ref_mg.GetFinestOperator()).LocalOperator();
  ComplexWrapperOperator local_reference(ref_local.Real(), ref_local.Imag());
  CheckPackedActions(local, local_reference);

  // Compare the RAP action, including PEC elimination and nonconforming/MPI maps,
  // against separate real/imaginary ParOperators (the existing four-real reference).
  const auto &fespace = space->GetNDSpace();
  ParOperator ref_real(*ref_local.Real(), fespace);
  std::unique_ptr<ParOperator> ref_imag;
  const auto &essential = space->GetNDDbcTDofLists().back();
  ref_real.SetEssentialTrueDofs(essential, Operator::DIAG_ONE);
  if (ref_local.Imag())
  {
    ref_imag = std::make_unique<ParOperator>(*ref_local.Imag(), fespace);
    ref_imag->SetEssentialTrueDofs(essential, Operator::DIAG_ZERO);
  }
  ComplexWrapperOperator reference(&ref_real, ref_imag.get());
  CheckPackedActions(*fine, reference);
  ComplexParOperator borrowed(ref_local.Real(), ref_local.Imag(), fespace);
  borrowed.SetEssentialTrueDofs(essential, Operator::DIAG_ONE);
  REQUIRE(IsOriginalComplexWrapper(borrowed.LocalOperator()));
  CheckPackedActions(borrowed, reference);

  // Main K/C/M remain unassembled-QData operators inside weighted SumOperators.
  if (lossy && !real_pc && !amr)
  {
    auto K = space->GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
    auto C = space->GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space->GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto A = space->GetSystemMatrix<ComplexOperator>(
        std::complex<double>{1.0}, std::complex<double>{0.0, 2.0},
        std::complex<double>{-4.0}, K.get(), C.get(), M.get(), nullptr);
    const auto *system = dynamic_cast<const ComplexParOperator *>(A.get());
    REQUIRE(system);
    REQUIRE(dynamic_cast<const SumOperator *>(system->LocalOperator().Real()));
    REQUIRE(IsOriginalComplexWrapper(system->LocalOperator()));
  }
}

TEST_CASE("libCEED packed complex unsupported input fallbacks",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  const auto [name, qdata, h1_space, tensor, scaled, real_only] =
      GENERATE(table<const char *, bool, bool, bool, bool, bool>(
          {{"Unassembled quadrature data", false, false, false, false, false},
           {"Tensor H1 basis", true, true, true, false, false},
           {"Non-tensor H1 basis", true, true, false, false, false},
           {"Input with dof multiplicity", true, false, false, true, false},
           {"Real-only input", true, false, false, false, true}}));
  CAPTURE(name);
  PackedIntegrationSettings settings(3);
  const char *file = tensor ? "/mesh/fichera-hex.mesh" : "/mesh/fichera-tet.mesh";
  auto mesh = Initialize(Mpi::World(), std::string(PALACE_TEST_DATA_DIR) + file, 0, false);
  mfem::ND_FECollection nd(3, 3);
  mfem::H1_FECollection h1(3, 3);
  FiniteElementSpace fespace(mesh, h1_space
                                       ? static_cast<mfem::FiniteElementCollection *>(&h1)
                                       : static_cast<mfem::FiniteElementCollection *>(&nd));
  auto coeff = BuildCoefficient(mesh, false, CoeffType::Matrix);
  BilinearForm ar(fespace), ai(fespace);
  if (h1_space)
  {
    ar.AddDomainIntegrator<DiffusionIntegrator>(coeff);
    ai.AddDomainIntegrator<DiffusionIntegrator>(coeff);
  }
  else
  {
    ar.AddDomainIntegrator<VectorFEMassIntegrator>(coeff);
    ai.AddDomainIntegrator<VectorFEMassIntegrator>(coeff);
  }
  if (qdata)
  {
    ar.AssembleQuadratureData();
    ai.AssembleQuadratureData();
  }
  auto real = ar.PartialAssemble(), imag = ai.PartialAssemble();
  auto ref_real = ar.PartialAssemble(), ref_imag = ai.PartialAssemble();
  if (scaled)
  {
    for (auto *op : {real.get(), ref_real.get()})
    {
      Vector multiplicity(fespace.GetVSize());
      multiplicity = 0.63;
      op->SetDofMultiplicity(std::move(multiplicity));
    }
  }
  if (real_only)
  {
    imag.reset();
    ref_imag.reset();
  }
  auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag));
  REQUIRE(IsOriginalComplexWrapper(*actual));
  ComplexWrapperOperator reference(ref_real.get(), ref_imag.get());
  CheckPackedActions(*actual, reference);
}

TEST_CASE("libCEED packed complex empty local operators",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  for (const int size : {0, 7})
  {
    auto real = std::make_unique<ceed::SymmetricOperator>(size, size);
    auto imag = std::make_unique<ceed::SymmetricOperator>(size, size);
    real->Finalize();
    imag->Finalize();
    auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag));
    REQUIRE(IsOriginalComplexWrapper(*actual));
    ComplexVector x(size), result(size), expected(size);
    x = std::complex<double>{0.7, -0.2};
    expected = 0.0;
    actual->Mult(x, result);
    CheckPackedResult(result, expected);
  }
}

// Opt in with "[ComplexPreconditionerBenchmark]" or "[Benchmark]". Orders 2/3/4
// run on a 768-tet mesh with one rank; --benchmark-ref-levels adds refinements.
// One rank keeps both volume and absorbing-boundary terms in the timed operator.
TEST_CASE("CPU complex preconditioner benchmark",
          "[.][libCEED][Benchmark][ComplexPreconditionerBenchmark][Serial]")
{
  const auto *config = Catch::getCurrentContext().getConfig();
  const auto &selectors = config->getTestsOrTags();
  const bool requested = std::any_of(
      selectors.begin(), selectors.end(),
      [](const std::string &selector)
      {
        return selector.find("[ComplexPreconditionerBenchmark]") != std::string::npos ||
               selector.find("[Benchmark]") != std::string::npos;
      });
  // The runner adds [Serial]/[Parallel], which also selects hidden tests. Respect
  // CTest's --skip-benchmarks before any setup, and require opt-in.
  if (config->skipBenchmarks() || !requested)
  {
    SKIP("Select [ComplexPreconditionerBenchmark] explicitly to run timing");
  }
  if (!PackedTestBackend())
  {
    SKIP("Packed application requires one CPU CEED context and one thread");
  }
  if (Mpi::Size(Mpi::World()) != 1)
  {
    SKIP("The local preconditioner benchmark requires one MPI rank");
  }
  const int order = GENERATE(2, 3, 4);
  const bool boundary = GENERATE(false, true);
  DYNAMIC_SECTION("p" << order << (boundary ? " with absorbing boundary" : " volume only"))
  {
    REQUIRE(benchmark_ref_levels >= 0);
    const auto comm = Mpi::World();
    ComplexPreconditionerFixture fixture(order, false, 2 + benchmark_ref_levels);
    if (!boundary)
    {
      fixture.boundaries.farfield.attributes.clear();
    }
    auto space = fixture.MakeSpace();
    Mpi::Barrier(comm);
    const double start = MPI_Wtime();
    auto pc = AssembleComplexPreconditioner(*space);
    double setup_time = MPI_Wtime() - start;
    Mpi::GlobalMax(1, &setup_time, comm);
    const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(pc.get());
    REQUIRE(mg);
    const auto *fine = dynamic_cast<const ComplexParOperator *>(&mg->GetFinestOperator());
    REQUIRE(fine);
    const auto &local = fine->LocalOperator();
    REQUIRE_FALSE(IsOriginalComplexWrapper(local));

    // Borrow the exact assembled real/imaginary inputs owned by the production
    // preconditioner to compare with the original four-real application.
    ComplexWrapperOperator local_reference(local.Real(), local.Imag());
    CheckPackedActions(local, local_reference);

    int elements = fixture.meshes.back()->GetNE();
    Mpi::GlobalSum(1, &elements, comm);
    const auto global_dofs = space->GetNDSpace().GlobalTrueVSize();
    if (Mpi::Root(comm))
    {
      WARN("Complex preconditioner: p = "
           << order << ", absorbing boundary = " << boundary
           << ", MPI ranks = " << Mpi::Size(comm) << ", global elements = " << elements
           << ", global true dofs = " << global_dofs
           << ", root local dofs = " << local.Height()
           << "\nProduction hierarchy setup (all levels, including packing) = "
           << 1.0e3 * setup_time << " ms");
      ComplexVector x(local.Width()), y(local.Height()), y_ref(local.Height());
      x.Real().Randomize(17);
      x.Imag().Randomize(53);
      // A complex scale selects the original wrapper's generic Mult/AXPY path.
      // Its real-only scale path requests negative real CEED AddMult coefficients,
      // which ceed::Operator does not support.
      const std::complex<double> scale{-0.7, 0.23};
      BENCHMARK("Local Mult (original borrowed wrapper)")
      {
        local_reference.Mult(x, y_ref);
        return y_ref.Size();
      };
      BENCHMARK("Local Mult (owned packed wrapper, including copies)")
      {
        local.Mult(x, y);
        return y.Size();
      };
      y_ref = 0.0;
      BENCHMARK("Local AddMult (original borrowed wrapper, complex scale)")
      {
        local_reference.AddMult(x, y_ref, scale);
        return y_ref.Size();
      };
      y = 0.0;
      BENCHMARK("Local AddMult (owned packed wrapper, complex scale, including copies)")
      {
        local.AddMult(x, y, scale);
        return y.Size();
      };
    }
    Mpi::Barrier(comm);
  }
}

TEST_CASE("3D libCEED Benchmarks", "[libCEED][Benchmark][Serial][Parallel]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh");
  RunCeedBenchmarks(MPI_COMM_WORLD, std::string(PALACE_TEST_DATA_DIR "/mesh/") + mesh,
                    benchmark_ref_levels, false, benchmark_order);
}

}  // namespace palace
