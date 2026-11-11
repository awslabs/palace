// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <typeinfo>
#include <ceed/backend.h>
#include <mfem.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/interfaces/catch_interfaces_config.hpp>
#include <catch2/internal/catch_context.hpp>
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/chebyshev.hpp"
#include "linalg/distrelaxation.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "test-helpers.hpp"
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
  if (Mpi::Size(comm) > smesh.GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
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
void TestCeedQuadratureData(MPI_Comm comm, const FiniteElementSpace &trial_fespace,
                            const FiniteElementSpace &test_fespace, T AddIntegrators)
{
  BilinearForm a_test(trial_fespace, test_fespace), a_ref(trial_fespace, test_fespace);
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

template <typename T>
void TestCeedQuadratureData(MPI_Comm comm, const FiniteElementSpace &fespace,
                            T AddIntegrators)
{
  TestCeedQuadratureData(comm, fespace, fespace, AddIntegrators);
}

// Coefficient stored for max_attr attributes with two materials: a nonzero one on every
// active_stride-th (libCEED, 1-based) attribute, and an exactly zero one everywhere else.
// With active_stride = 0, every attribute is exactly zero.
auto BuildZeroOutsideCoefficient(std::size_t max_attr, int active_stride, double value)
{
  mfem::Array<int> attr_mat(static_cast<int>(max_attr));
  mfem::DenseTensor mat_coeff(1, 1, 2);
  mat_coeff(0, 0, 0) = value;
  mat_coeff(0, 0, 1) = 0.0;
  for (int i = 0; i < attr_mat.Size(); i++)
  {
    attr_mat[i] = (active_stride > 0 && i % active_stride == 0) ? 0 : 1;
  }
  return MaterialPropertyCoefficient(attr_mat, mat_coeff);
}

void TestCeedActiveGeometryRestriction(MPI_Comm comm, const Mesh &mesh)
{
  long long count[4] = {0, 0, 0, 0};
  for (auto ceed : ceed::internal::GetCeedObjects())
  {
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      const bool use_bdr = (mfem::Geometry::Dimension[geom] != mesh.Dimension());
      REQUIRE(std::is_sorted(data.active_indices.begin(), data.active_indices.end()));
      REQUIRE(std::adjacent_find(data.active_indices.begin(), data.active_indices.end()) ==
              data.active_indices.end());
      for (auto i : data.active_indices)
      {
        REQUIRE(i >= 0);
        REQUIRE(static_cast<std::size_t>(i) < data.indices.size());
      }
      count[2 * use_bdr] += data.active_indices.size();
      count[2 * use_bdr + 1] += data.indices.size();

      if (data.active_indices.empty() || data.active_indices.size() == data.indices.size())
      {
        REQUIRE(data.active_geom_data_restr == nullptr);
        continue;
      }
      REQUIRE(data.active_geom_data_restr != nullptr);

      // Transpose a constant E-vector through the offsets restriction. Geometry entries for
      // selected full-list positions must be touched exactly once, and all others not at
      // all, in whatever layout the backend stores the geometry vector.
      CeedInt num_elem, num_qpts, num_comp, layout[3];
      CeedSize l_size, geom_data_size;
      PalaceCeedCall(ceed, CeedElemRestrictionGetLLayout(data.geom_data_restr, layout));
      PalaceCeedCall(
          ceed, CeedElemRestrictionGetNumElements(data.active_geom_data_restr, &num_elem));
      PalaceCeedCall(
          ceed, CeedElemRestrictionGetElementSize(data.active_geom_data_restr, &num_qpts));
      PalaceCeedCall(ceed, CeedElemRestrictionGetNumComponents(data.active_geom_data_restr,
                                                               &num_comp));
      PalaceCeedCall(
          ceed, CeedElemRestrictionGetLVectorSize(data.active_geom_data_restr, &l_size));
      PalaceCeedCall(ceed, CeedVectorGetLength(data.geom_data, &geom_data_size));
      REQUIRE(num_elem == static_cast<CeedInt>(data.active_indices.size()));
      REQUIRE(l_size == geom_data_size);

      CeedVector l_vec, e_vec;
      PalaceCeedCall(ceed, CeedElemRestrictionCreateVector(data.active_geom_data_restr,
                                                           &l_vec, &e_vec));
      PalaceCeedCall(ceed, CeedVectorSetValue(l_vec, 0.0));
      PalaceCeedCall(ceed, CeedVectorSetValue(e_vec, 1.0));
      PalaceCeedCall(ceed,
                     CeedElemRestrictionApply(data.active_geom_data_restr, CEED_TRANSPOSE,
                                              e_vec, l_vec, CEED_REQUEST_IMMEDIATE));
      const CeedScalar *l_array;
      PalaceCeedCall(ceed, CeedVectorGetArrayRead(l_vec, CEED_MEM_HOST, &l_array));
      std::vector<bool> active(data.indices.size(), false);
      for (auto i : data.active_indices)
      {
        active[i] = true;
      }
      for (std::size_t p = 0; p < data.indices.size(); p++)
      {
        for (CeedInt c = 0; c < num_comp; c++)
        {
          for (CeedInt q = 0; q < num_qpts; q++)
          {
            const CeedSize i = static_cast<CeedSize>(q) * layout[0] +
                               static_cast<CeedSize>(c) * layout[1] +
                               static_cast<CeedSize>(p) * layout[2];
            CHECK(l_array[i] == (active[p] ? 1.0 : 0.0));
          }
        }
      }
      PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(l_vec, &l_array));
      PalaceCeedCall(ceed, CeedVectorDestroy(&l_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&e_vec));
    }
  }
  Mpi::GlobalSum(4, count, comm);
  CHECK(count[0] > 0);
  CHECK(count[0] < count[1]);
  CHECK(count[2] > 0);
  CHECK(count[2] < count[3]);
}

// Skipping elements whose contribution is exactly 0.0 cannot change a sum, so the operator
// application is bit-identical, but only as long as the remaining contributions are
// accumulated in the same order. Blocked libCEED backends (the default when libXSMM is
// available) accumulate the element contributions of a block of block_size elements in
// element-node-major order (CeedElemRestrictionApply*Transpose_Ref_Core), so dropping
// elements regroups the blocks and reorders the sums for degrees of freedom shared by
// elements which are no longer in the same block, GPU backends scatter shared degrees of
// freedom from different launch shapes (the MAGMA backends are nondeterministic in any
// case), and with several CPU Ceed contexts (OpenMP threads) ceed::Operator applies the
// contexts concurrently into the shared output vector. Those sums can then differ in the
// last bit. Only a single unblocked CPU backend preserves the accumulation order exactly.
bool CeedBackendAccumulatesInOrder()
{
  const auto backend = ceed::Print();
  return ceed::internal::NumCeeds() == 1 && backend.find("blocked") == std::string::npos &&
         backend.find("/gpu/") == std::string::npos;
}

void TestVectorAgree(const Vector &y_ref, const Vector &y_test)
{
  REQUIRE(y_test.Size() == y_ref.Size());
  const double tol = CeedBackendAccumulatesInOrder()
                         ? 0.0
                         : 8.0 * std::numeric_limits<double>::epsilon() * y_ref.Normlinf();
  const auto *d_ref = y_ref.HostRead();
  const auto *d_test = y_test.HostRead();
  int num_diff = 0;
  for (int i = 0; i < y_ref.Size(); i++)
  {
    // Written so that a NaN in either vector counts as a difference.
    num_diff += !(std::abs(d_test[i] - d_ref[i]) <= tol);
  }
  REQUIRE(num_diff == 0);
}

// Compare two operators for agreeing application, transpose application and diagonal (see
// TestVectorAgree), guarding against a vacuous comparison of two zero operators.
void TestCeedOperatorMultAgree(MPI_Comm comm, const Operator &op_test,
                               const Operator &op_ref)
{
  Vector x(op_ref.Width()), y_ref(op_ref.Height()), y_test(op_ref.Height());
  x.UseDevice(true);
  y_ref.UseDevice(true);
  y_test.UseDevice(true);
  x.Randomize(1);

  op_ref.Mult(x, y_ref);
  op_test.Mult(x, y_test);
  double norm_ref = y_ref * y_ref;
  Mpi::GlobalSum(1, &norm_ref, comm);
  REQUIRE(norm_ref > 0.0);
  TestVectorAgree(y_ref, y_test);

  op_ref.MultTranspose(x, y_ref);
  op_test.MultTranspose(x, y_test);
  TestVectorAgree(y_ref, y_test);

  Vector d_ref(op_ref.Height()), d_test(op_ref.Height());
  d_ref.UseDevice(true);
  d_test.UseDevice(true);
  op_ref.AssembleDiagonal(d_ref);
  op_test.AssembleDiagonal(d_test);
  TestVectorAgree(d_ref, d_test);
}

// The operator of domain integrators assembled over the shared subset covers exactly its
// active elements, without a copy of the geometry data: sub-operators which read the
// geometry factors refer to the Mesh-owned full geometry vector through the matching subset
// restriction, and those with assembled quadrature data (on CPUs, see
// ceed::AssembleAllQuadratureData) do not read the geometry factors at all.
void TestCeedOperatorUsesSharedGeometry(MPI_Comm comm, const Mesh &mesh,
                                        const ceed::Operator &op)
{
  long long count[3] = {0, 0, 0};  // Sub-operators, their elements, active elements
  for (std::size_t i = 0; i < op.Size(); i++)
  {
    // The mesh caches geometry data by the Palace-owned Ceed contexts, not by the delegate
    // Ceed which CeedOperatorGetCeed reports for operators created through a delegate.
    Ceed ceed = ceed::internal::GetCeedObjects()[i];
    CeedInt num_sub_ops;
    CeedOperator *sub_ops;
    PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub(op[i], &num_sub_ops));
    PalaceCeedCall(ceed, CeedOperatorCompositeGetSubList(op[i], &sub_ops));
    const auto &geom_data = mesh.GetCeedGeomFactorData(ceed);
    for (const auto &[geom, data] : geom_data)
    {
      if (mfem::Geometry::Dimension[geom] == mesh.Dimension())
      {
        count[2] += data.active_indices.size();
      }
    }
    for (CeedInt j = 0; j < num_sub_ops; j++)
    {
      CeedInt num_elem;
      PalaceCeedCall(ceed, CeedOperatorGetNumElements(sub_ops[j], &num_elem));
      count[0]++;
      count[1] += num_elem;

      CeedOperatorField field;
      PalaceCeedCall(ceed, CeedOperatorGetFieldByName(sub_ops[j], "geom_data", &field));
      if (!field)
      {
        PalaceCeedCall(ceed, CeedOperatorGetFieldByName(sub_ops[j], "q_data", &field));
        REQUIRE(field != nullptr);
        continue;
      }
      CeedVector field_vec;
      CeedElemRestriction field_restr;
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &field_vec));
      PalaceCeedCall(ceed, CeedOperatorFieldGetElemRestriction(field, &field_restr));

      // The subset operator must refer to the Mesh-owned full geometry vector, not a copied
      // subset vector, through the matching offsets restriction.
      const ceed::CeedGeomFactorData *matching_data = nullptr;
      for (const auto &[geom, data] : geom_data)
      {
        if (field_vec == data.geom_data)
        {
          matching_data = &data;
          break;
        }
      }
      REQUIRE(matching_data != nullptr);
      CHECK(field_restr ==
            matching_data->GetGeomDataRestriction(ceed::CeedElementSubset::Active));
    }
  }
  Mpi::GlobalSum(3, count, comm);
  REQUIRE(count[0] > 0);
  CHECK(count[1] == count[2]);
}

// Assemble the same integrators over the full mesh and the shared lossy/active subsets, and
// check that the two operators agree. Full assembly of the restricted operator has fewer
// stored entries, which is why the subset is never used for a fully assembled matrix (the
// sparsity pattern would depend on the coefficient values).
template <typename T>
void TestCeedSkipZeroCoefficientElements(MPI_Comm comm, const FiniteElementSpace &fespace,
                                         bool assemble_q_data, T AddIntegrators,
                                         bool test_shared_geometry = false)
{
  BilinearForm a_ref(fespace), a_skip(fespace);
  AddIntegrators(a_ref);
  AddIntegrators(a_skip);
  a_skip.SkipZeroCoefficientElements();
  if (assemble_q_data)
  {
    a_ref.AssembleQuadratureData();
    a_skip.AssembleQuadratureData();
  }
  auto op_ref = a_ref.PartialAssemble();
  auto op_skip = a_skip.PartialAssemble();
  TestCeedOperatorMultAgree(comm, *op_skip, *op_ref);
  if (test_shared_geometry)
  {
    REQUIRE(!assemble_q_data);
    TestCeedOperatorUsesSharedGeometry(comm, fespace.GetMesh(), *op_skip);
  }

  constexpr bool skip_zeros = false;
  auto mat_ref = BilinearForm::FullAssemble(*op_ref, skip_zeros);
  auto mat_skip = BilinearForm::FullAssemble(*op_skip, skip_zeros);
  long long nnz[2] = {static_cast<long long>(mat_ref->NNZ()),
                      static_cast<long long>(mat_skip->NNZ())};
  Mpi::GlobalSum(2, nnz, comm);
  REQUIRE(nnz[1] < nnz[0]);
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
    std::size_t mem_test =
        (Q * q_data_size * 8 + P * 4) * static_cast<std::size_t>(mesh.GetNE());
    std::stringstream msg;
    msg << "benchmark memory estimate:\n"
        << "  N = " << fespace.GetVSize() << " (NE = " << mesh.GetNE() << ", P = " << P
        << ", Q = " << Q << ")\n";
    if (nnz > 0)
    {
      msg << "  Full Assembly = " << mem_ref / static_cast<double>(1024 * 1024) << " MB ("
          << nnz << " NNZ)\n";
    }
    else
    {
      msg << "  Full Assembly = N/A (skipped)\n";
    }
    msg << "  Partial Assembly = " << mem_test / static_cast<double>(1024 * 1024)
        << " MB\n";
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
    std::size_t mem_test =
        (trial_P * 4 + test_P * 4) * static_cast<std::size_t>(mesh.GetNE());
    std::stringstream msg;
    msg << "benchmark memory estimate:\n"
        << "  N = " << trial_fespace.GetVSize() << ", " << test_fespace.GetVSize()
        << " (NE = " << mesh.GetNE() << ", P = " << trial_P << ", " << test_P << ")\n";
    if (nnz > 0)
    {
      msg << "  Full Assembly = " << mem_ref / static_cast<double>(1024 * 1024) << " MB ("
          << nnz << " NNZ)\n";
    }
    else
    {
      msg << "  Full Assembly = N/A (skipped)\n";
    }
    msg << "  Partial Assembly = " << mem_test / static_cast<double>(1024 * 1024)
        << " MB\n";
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

  // Tests on L2 spaces of order p - 1, with VALUE (u = ̂u) and INTEGRAL (u = 1 / det(J) ̂u)
  // map types, at MFEM's integration order for them. Boundary integrators do not apply.
  SECTION("L2 Integrators")
  {
    mfem::H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(mesh, &h1_fec);
    const bool add = !bdr_integ && coeff_type != CoeffType::Matrix;
    SECTION("L2 Mass Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_extra_pk = -2;
      fem::DefaultIntegrationOrder::q_order_extra_qk = -2;
      mesh.ResetCeedObjects();
      for (const int map_type : {mfem::FiniteElement::VALUE, mfem::FiniteElement::INTEGRAL})
      {
        INFO("Map type: " << map_type);
        mfem::L2_FECollection l2_fec(order - 1, dim, mfem::BasisType::GaussLegendre,
                                     map_type);
        FiniteElementSpace l2_fespace(mesh, &l2_fec);
        BilinearForm a_test(l2_fespace);
        mfem::BilinearForm a_ref(&l2_fespace.Get());
        if (add && coeff_type == CoeffType::Const)
        {
          AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref);
        }
        else if (add)
        {
          AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref, Q,
                                                               Q_ref);
        }
        TestCeedOperator(a_test, a_ref);
      }
    }
    SECTION("Mixed L2-H1 Mass Integrator")
    {
      fem::DefaultIntegrationOrder::q_order_extra_pk = -1;
      fem::DefaultIntegrationOrder::q_order_extra_qk = -1;
      mesh.ResetCeedObjects();
      for (const int map_type : {mfem::FiniteElement::VALUE, mfem::FiniteElement::INTEGRAL})
      {
        INFO("Map type: " << map_type);
        mfem::L2_FECollection l2_fec(order - 1, dim, mfem::BasisType::GaussLegendre,
                                     map_type);
        FiniteElementSpace l2_fespace(mesh, &l2_fec);
        for (const bool l2_trial : {true, false})
        {
          INFO("L2 trial: " << l2_trial);
          auto &trial_fespace = l2_trial ? l2_fespace : h1_fespace;
          auto &test_fespace = l2_trial ? h1_fespace : l2_fespace;
          BilinearForm a_test(trial_fespace, test_fespace);
          mfem::MixedBilinearForm a_ref(&trial_fespace.Get(), &test_fespace.Get());
          if (add && coeff_type == CoeffType::Const)
          {
            AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref);
          }
          else if (add)
          {
            AddIntegrators<MassIntegrator, mfem::MassIntegrator>(bdr_integ, a_test, a_ref,
                                                                 Q, Q_ref);
          }
          TestCeedOperator(a_test, a_ref);
        }
      }
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

// One test case per mesh, so CTest can run the meshes in parallel. The mesh is the
// template argument, so the case list is the mesh list and the case name shows the mesh.
constexpr char StarQuad[] = "star-quad.mesh";
constexpr char StarTri[] = "star-tri.mesh";
constexpr char StarMixedP2[] = "star-mixed-p2.mesh";
constexpr char FicheraHex[] = "fichera-hex.mesh";
constexpr char FicheraTet[] = "fichera-tet.mesh";
constexpr char FicheraMixedP2[] = "fichera-mixed-p2.mesh";

TEMPLATE_TEST_CASE_SIG("2D libCEED Operators", "[libCEED][Serial][Parallel]",
                       ((const char *Mesh), Mesh), StarQuad, StarTri, StarMixedP2)
{
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, TestMeshPath(Mesh), 0, amr, order);
}

TEMPLATE_TEST_CASE_SIG("3D libCEED Operators", "[libCEED][Serial][Parallel]",
                       ((const char *Mesh), Mesh), FicheraHex, FicheraTet, FicheraMixedP2)
{
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, TestMeshPath(Mesh), 0, amr, order);
}

TEST_CASE("2D libCEED Interpolators", "[libCEED][Interpolator][Serial][Parallel]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh", "star-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, TestMeshPath(mesh), 0, amr, order);
}

TEST_CASE("3D libCEED Interpolators", "[libCEED][Interpolator][Serial][Parallel]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto amr = GENERATE(false, true);
  auto order = GENERATE(1, 2, 3);
  RunCeedInterpolatorTests(MPI_COMM_WORLD, TestMeshPath(mesh), 0, amr, order);
}

// Test the 31 (SpaceDim=3, Dim=1) libCEED qfunctions for boundary integrators on a 2D
// mesh embedded in 3D. Uses the standard RunCeedIntegratorTests with embed_in_3d=true.
// Matrix coefficients are excluded (need sdim x sdim, not dim x dim).
TEST_CASE("2D-in-3D libCEED Boundary Operators", "[libCEED][Serial][Parallel]")
{
  auto mesh = GENERATE("star-quad.mesh", "star-tri.mesh");
  auto order = GENERATE(1, 2, 3);
  RunCeedIntegratorTests(MPI_COMM_WORLD, TestMeshPath(mesh), 0, false, order, true);
}

TEST_CASE("libCEED Reduced Bases", "[libCEED][Serial][Parallel]")
{
  // Domain integrators of H(curl) and H(div) spaces on affine simplices use bases with
  // fewer points than the default quadrature (also on the affine elements of a curved
  // mesh). The partially assembled operators agree with the fully assembled ones, which use
  // the default quadrature, to rounding.
  auto mesh_file = GENERATE("fichera-tet.mesh", "fichera-mixed-p2.mesh");
  auto order = GENERATE(1, 2, 3, 4);
  const auto comm = MPI_COMM_WORLD;
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, false);
  const int dim = mesh.Dimension();
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = false;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  INFO("Mesh: " << mesh_file << ", order: " << order);

  auto Q = BuildCoefficient(mesh, false, CoeffType::Scalar);
  auto MQ = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto TestReducedBases = [&](BilinearForm &a)
  {
    a.AssembleQuadratureData();
    auto op = a.PartialAssemble();
    auto mat = a.FullAssemble(false);
    Vector x(op->Width()), y(op->Height()), y_ref(op->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    y_ref.UseDevice(true);
    linalg::SetRandom(comm, x);
    op->Mult(x, y);
    mat->Mult(x, y_ref);
    y -= y_ref;
    REQUIRE(linalg::Norml2(comm, y) <= 1.0e-12 * linalg::Norml2(comm, y_ref));
  };
  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), rt_fespace(mesh, &rt_fec);
  SECTION("Curl-Curl + Mass Integrator")
  {
    BilinearForm a(nd_fespace);
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
    TestReducedBases(a);
  }
  SECTION("Curl-Curl Integrator")
  {
    BilinearForm a(nd_fespace);
    a.AddDomainIntegrator<CurlCurlIntegrator>(MQ);
    TestReducedBases(a);
  }
  SECTION("H(div) Mass Integrator")
  {
    BilinearForm a(rt_fespace);
    a.AddDomainIntegrator<VectorFEMassIntegrator>(MQ);
    TestReducedBases(a);
  }
  SECTION("Div-Div + Mass Integrator")
  {
    BilinearForm a(rt_fespace);
    a.AddDomainIntegrator<DivDivMassIntegrator>(Q, MQ);
    TestReducedBases(a);
  }
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
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, false);
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
  mfem::L2_FECollection l2_fec(order - 1, dim, mfem::BasisType::GaussLegendre,
                               mfem::FiniteElement::INTEGRAL);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec),
      l2_fespace(mesh, &l2_fec);

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
  // L2 spaces with INTEGRAL map type: the 2D curl space and the flux projection of its
  // error estimator.
  SECTION("L2 Mass")
  {
    TestCeedQuadratureData(comm, l2_fespace, [&](BilinearForm &a)
                           { a.AddDomainIntegrator<MassIntegrator>(Q); });
  }
  SECTION("Mixed L2-H1 Mass")
  {
    TestCeedQuadratureData(comm, l2_fespace, h1_fespace, [&](BilinearForm &a)
                           { a.AddDomainIntegrator<MassIntegrator>(Q); });
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
  return mem == CEED_MEM_HOST && std::string(resource).starts_with("/cpu/");
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
    // These operations retain the original real/imaginary operators. Check representative
    // spaces and geometries, including the true-dof diagonal used by multigrid smoothers.
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

  ComplexPreconditionerFixture(int order, bool amr, int ref_levels = 0,
                               bool partial_lossy = false)
    : settings(order)
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
    if (partial_lossy)
    {
      // Half of the elements get a second, lossless material, so that the imaginary terms
      // cover only part of the mesh. Interleaved by element index to give every process a
      // chance of holding both sets.
      for (int i = 0; i < smesh.GetNE(); i++)
      {
        smesh.SetAttribute(i, 1 + (i % 2));
      }
      smesh.SetAttributes();
    }
    meshes.push_back(std::make_unique<Mesh>(Mpi::World(), smesh));
    config::MaterialData material;
    material.attributes = {1};
    material.epsilon_r = 2.08;
    material.tandelta = 4.0e-4;
    domains.attributes = {1};
    domains.materials = {material};
    if (partial_lossy)
    {
      config::MaterialData lossless;
      lossless.attributes = {2};
      lossless.epsilon_r = 1.5;
      domains.attributes = {1, 2};
      domains.materials = {material, lossless};
    }
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
  const auto [name, mesh_file, order, curl, boundary, nested, nonsymmetric, subset] =
      GENERATE(table<const char *, const char *, int, bool, bool, bool, bool, bool>(
          {{"H(div) mass/mass and shared QData", "fichera-tet.mesh", 1, false, false, false,
            false, false},
           {"Curlmass/mass and p3 face orientations", "fichera-tet.mesh", 3, true, false,
            false, false, false},
           {"Hexahedron and unmatched boundary terms", "fichera-hex.mesh", 3, true, true,
            false, false, false},
           {"Direct volume pair and nested remainder", "fichera-tet.mesh", 2, true, false,
            true, false, false},
           {"Mixed H(curl) pairs and nonsymmetric tensors", "fichera-mixed-p2.mesh", 2,
            true, false, false, true, false},
           {"Mixed H(div) pairs", "fichera-mixed-p2.mesh", 1, false, false, false, false,
            false},
           {"Pyramid H(curl) pair", nullptr, 1, true, false, false, false, false},
           {"Pyramid H(div) pair", nullptr, 1, false, false, false, false, false},
           {"H(div) pair over a subset", "fichera-tet.mesh", 1, false, false, false, false,
            true},
           {"H(curl) hexahedra, subset and boundary remainder", "fichera-hex.mesh", 2, true,
            true, false, false, true}}));
  CAPTURE(name);
  // BilinearForm labels square operators symmetric; the nonsymmetric coefficient case
  // therefore exercises forward application only.
  const bool inherited = !nested && !nonsymmetric;
  PackedIntegrationSettings settings(order);
  auto mesh = [mesh_file = mesh_file]()
  {
    if (mesh_file)
    {
      return Initialize(Mpi::World(), TestMeshPath(mesh_file), 0, false);
    }
    // The mixed Fichera mesh has no pyramids. Use low-order elements on a small
    // generated mesh, with enough elements for every MPI rank and no AMR.
    auto smesh =
        mfem::Mesh::MakeCartesian3D(Mpi::Size(Mpi::World()), 1, 1, mfem::Element::PYRAMID);
    return Mesh(Mpi::World(), smesh);
  }();
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  if (curl)
  {
    fec = std::make_unique<mfem::ND_FECollection>(order, 3);
  }
  else
  {
    fec = std::make_unique<mfem::RT_FECollection>(order - 1, 3);
  }
  FiniteElementSpace fespace(mesh, fec.get());
  auto real_mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto imag_mass = subset ? BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 2, 1.0)
                          : BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto real_curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto bdr_mass = BuildCoefficient(mesh, true, CoeffType::Matrix);
  auto bdr_curl = BuildCoefficient(mesh, true, CoeffType::Scalar);
  if (nonsymmetric)
  {
    // Unequal off-diagonal entries expose tensor-index transposition in both the
    // complex mass action and the shared real curl action.
    for (auto *coeff : {&real_mass, &imag_mass, &real_curl})
    {
      mfem::DenseTensor values(coeff->GetMaterialProperties());
      for (int k = 0; k < values.SizeK(); k++)
      {
        values(0, 1, k) = 0.7;
        values(1, 0, k) = -0.4;
        values(0, 2, k) = -0.2;
        values(2, 0, k) = 0.5;
        values(1, 2, k) = 0.3;
        values(2, 1, k) = -0.8;
      }
      *coeff = MaterialPropertyCoefficient(coeff->GetAttributeToMaterial(), values);
    }
  }
  real_mass *= -0.71;
  imag_mass *= 0.19;
  real_curl *= 1.37;
  bdr_mass *= -0.23;
  if (subset)
  {
    std::vector<int> lossy_attr;
    imag_mass.AddNonzeroAttributes(lossy_attr);
    mesh.SetCeedActiveAttributes(std::move(lossy_attr), {});
  }
  BilinearForm ar(fespace), ai(fespace);
  if (subset)
  {
    ar.SkipZeroCoefficientElements();
    ai.SkipZeroCoefficientElements();
  }
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
  if (subset)
  {
    long long num_subsets = 0;
    for (const auto &[geom, data] :
         mesh.GetCeedGeomFactorData(ceed::internal::GetCeedObjects()[0]))
    {
      num_subsets += (mfem::Geometry::Dimension[geom] == mesh.Dimension() &&
                      !data.active_indices.empty() &&
                      data.active_indices.size() < data.indices.size());
    }
    Mpi::GlobalSum(1, &num_subsets, Mpi::World());
    int max_ne = mesh.GetNE();
    Mpi::GlobalMax(1, &max_ne, Mpi::World());
    if (max_ne < 2)
    {
      SKIP("No process owns enough elements for a strict subset");
    }
    REQUIRE(num_subsets > 0);
  }
  if (mesh_file && std::string(mesh_file) == "fichera-mixed-p2.mesh" &&
      Mpi::Size(Mpi::World()) == 1)
  {
    // Each volume geometry contributes a compatible pair to the same composite,
    // exercising multiple iterations of the pairing loop, including prism bases.
    for (auto *op : {real.get(), imag.get()})
    {
      CeedInt count;
      REQUIRE(CeedOperatorCompositeGetNumSub((*op)[0], &count) == 0);
      REQUIRE(count == 3);
    }
    REQUIRE(mesh.Get().HasGeometry(mfem::Geometry::TETRAHEDRON));
    REQUIRE(mesh.Get().HasGeometry(mfem::Geometry::CUBE));
    REQUIRE(mesh.Get().HasGeometry(mfem::Geometry::PRISM));
  }
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
  auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag), fespace);
  // A rank with no subset element has no local imaginary volume leaf to pair.
  long long num_packed = !IsOriginalComplexWrapper(*actual);
  Mpi::GlobalSum(1, &num_packed, Mpi::World());
  if (subset)
  {
    REQUIRE(num_packed > 0);
  }
  else
  {
    REQUIRE(num_packed == Mpi::Size(Mpi::World()));
  }
  if (mesh_file && std::string(mesh_file) == "fichera-tet.mesh" && !nested)
  {
    // Verify both restriction clone paths: ND face transformations and RT signs.
    Ceed ceed = ceed::internal::GetCeedObjects()[0];
    const auto &geom = mesh.GetCeedGeomFactorData(ceed).at(mfem::Geometry::TETRAHEDRON);
    CeedRestrictionType type;
    REQUIRE(CeedElemRestrictionGetType(
                fespace.GetCeedElemRestriction(ceed, mfem::Geometry::TETRAHEDRON, geom),
                &type) == 0);
    REQUIRE(type == (curl ? CEED_RESTRICTION_CURL_ORIENTED : CEED_RESTRICTION_ORIENTED));
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
  const auto [name, amr, lossy, real_pc, partial_lossy] =
      GENERATE(table<const char *, bool, bool, bool, bool>(
          {{"Conforming lossy complex preconditioner", false, true, false, false},
           {"Nonconforming lossy complex preconditioner", true, true, false, false},
           {"Conforming lossless volume with absorbing boundary", false, false, false,
            false},
           {"Conforming real-only preconditioner", false, true, true, false},
           {"Conforming partially lossy preconditioner", false, true, false, true}}));
  CAPTURE(name);
  ComplexPreconditionerFixture fixture(3, amr, 0, partial_lossy);
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
  // Volume and boundary leaves without an imaginary partner are packed (or element
  // assembled) as well, so the fine and H1 auxiliary operators all use the fast path.
  const auto &aux =
      dynamic_cast<const ComplexParOperator &>(mg->GetFinestAuxiliaryOperator());
  REQUIRE(!IsOriginalComplexWrapper(aux.LocalOperator()));
  const auto &local = fine->LocalOperator();
  REQUIRE(!IsOriginalComplexWrapper(local));
  // Assemble independent reference operators and exercise the unchanged borrowed RAP
  // constructor. Numerical agreement is checked against separate real ParOperators below.
  auto ref_pc = AssembleComplexPreconditioner(*space);
  const auto &ref_mg = dynamic_cast<const ComplexMultigridOperator &>(*ref_pc);
  const auto &ref_local =
      dynamic_cast<const ComplexParOperator &>(ref_mg.GetFinestOperator()).LocalOperator();
  ComplexWrapperOperator local_reference(ref_local.Real(), ref_local.Imag());
  CheckPackedActions(local, local_reference);
  if (partial_lossy)
  {
    long long num_subsets = 0;
    const auto &mesh = space->GetNDSpace().GetMesh();
    for (const auto &[geom, data] :
         mesh.GetCeedGeomFactorData(ceed::internal::GetCeedObjects()[0]))
    {
      num_subsets += (mfem::Geometry::Dimension[geom] == mesh.Dimension() &&
                      !data.active_indices.empty() &&
                      data.active_indices.size() < data.indices.size());
    }
    Mpi::GlobalSum(1, &num_subsets, Mpi::World());
    int max_ne = mesh.GetNE();
    Mpi::GlobalMax(1, &max_ne, Mpi::World());
    if (max_ne < 2)
    {
      SKIP("No process owns enough elements for a strict subset");
    }
    REQUIRE(num_subsets > 0);

    long long num_packed = 0;
    for (std::size_t l = 1; l < mg->GetNumLevels(); l++)
    {
      const auto &level =
          dynamic_cast<const ComplexParOperator &>(mg->GetOperatorAtLevel(l))
              .LocalOperator();
      const auto &level_ref =
          dynamic_cast<const ComplexParOperator &>(ref_mg.GetOperatorAtLevel(l))
              .LocalOperator();
      num_packed += !IsOriginalComplexWrapper(level);
      ComplexWrapperOperator reference(level_ref.Real(), level_ref.Imag());
      CheckPackedActions(level, reference);
    }
    Mpi::GlobalSum(1, &num_packed, Mpi::World());
    REQUIRE(num_packed > 0);
  }

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
  CheckPackedActions(*fine, reference, lossy && !real_pc);
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

TEST_CASE("libCEED packed complex input fallbacks",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Packed application requires one CPU CEED context and one thread");
  }
  const auto [name, qdata, h1_space, tensor, scaled, real_only] =
      GENERATE(table<const char *, bool, bool, bool, bool, bool>(
          {{"Unassembled quadrature data", false, false, false, false, false},
           {"Tensor H1 basis", true, true, true, false, false},
           {"Non-tensor H1 basis", true, true, false, false, false},
           {"Input with dof multiplicity", true, false, false, true, false},
           {"Real-only input", true, false, false, false, true}}));
  CAPTURE(name);
  PackedIntegrationSettings settings(3);
  const char *file = tensor ? "fichera-hex.mesh" : "fichera-tet.mesh";
  auto mesh = Initialize(Mpi::World(), TestMeshPath(file), 0, false);
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
  auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag), fespace);
  // Non-tensor H1 bases and real-only inputs are element assembled or packed;
  // tensor-product bases, dof multiplicity, and unassembled vector quadrature data keep the
  // original action. On CPUs, operators assemble their quadrature data anyway, unless
  // disabled (see ceed::AssembleAllQuadratureData).
  const bool assembled = qdata || ceed::AssembleAllQuadratureData();
  REQUIRE(IsOriginalComplexWrapper(*actual) == (tensor || scaled || !assembled));
  ComplexWrapperOperator reference(ref_real.get(), ref_imag.get());
  CheckPackedActions(*actual, reference);
}

TEST_CASE("libCEED packed complex empty local operators",
          "[libCEED][ComplexPacked][Serial][Parallel]")
{
  // Empty composites are rejected before the space is consulted, but one is still required
  // to construct the operator.
  auto smesh = mfem::Mesh::MakeCartesian3D(Mpi::Size(Mpi::World()), 1, 1,
                                           mfem::Element::TETRAHEDRON);
  Mesh mesh(Mpi::World(), smesh);
  mfem::ND_FECollection fec(1, 3);
  FiniteElementSpace fespace(mesh, &fec);
  for (const int size : {0, 7})
  {
    auto real = std::make_unique<ceed::SymmetricOperator>(size, size);
    auto imag = std::make_unique<ceed::SymmetricOperator>(size, size);
    real->Finalize();
    imag->Finalize();
    auto actual = ceed::CreateComplexOperator(std::move(real), std::move(imag), fespace);
    REQUIRE(IsOriginalComplexWrapper(*actual));
    ComplexVector x(size), result(size), expected(size);
    x = std::complex<double>{0.7, -0.2};
    expected = 0.0;
    actual->Mult(x, result);
    CheckPackedResult(result, expected);
  }
}

namespace
{

// Relative difference, over all processes, of a single-precision result from its
// double-precision reference.
double RelativeDifference(const ComplexVector &actual, const ComplexVector &expected)
{
  ComplexVector diff(expected);
  diff.AXPY(-1.0, actual);
  return linalg::Norml2(Mpi::World(), diff) / linalg::Norml2(Mpi::World(), expected);
}

double RelativeDifference(const Vector &actual, const Vector &expected)
{
  Vector diff(expected);
  diff.Add(-1.0, actual);
  return linalg::Norml2(Mpi::World(), diff) / linalg::Norml2(Mpi::World(), expected);
}

// The rounding of the element matrices to single precision gives relative errors of a few
// 1e-8 per entry: the approximation must be that close to the reference, and not closer
// than double-precision rounding (the single-precision path was used).
constexpr double single_min = 1.0e-11, single_max = 6.0e-7;

}  // namespace

TEST_CASE("libCEED single-precision operator approximation",
          "[libCEED][SinglePrecision][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Single-precision application requires one CPU CEED context and one thread");
  }
  const auto [name, mesh_file, ref_levels, order, imag_scale, boundary, subset,
              dropped_expected] =
      GENERATE(table<const char *, const char *, int, int, double, bool, bool, bool>(
          {{"p1 lossy tetrahedra", "fichera-tet.mesh", 0, 1, 0.19, false, false, false},
           {"p3 lossy tetrahedra with boundary terms", "fichera-tet.mesh", 0, 3, 0.19, true,
            false, false},
           {"p2 lossy subset of tetrahedra", "fichera-tet.mesh", 0, 2, 0.19, false, true,
            false},
           {"p2 negligible loss with boundary terms", "fichera-tet.mesh", 0, 2, 1.0e-9,
            true, false, true},
           {"p2 mixed elements with a hexahedral remainder", "fichera-mixed-p2.mesh", 0, 2,
            0.19, false, false, false},
           // More elements per process than the element matrices of one assembly chunk.
           {"p4 lossy tetrahedra assembled in chunks", "fichera-tet.mesh", 2, 4, 0.19,
            false, false, false}}));
  CAPTURE(name);
  PackedIntegrationSettings settings(order);
  auto mesh = Initialize(Mpi::World(), TestMeshPath(mesh_file), ref_levels, false);
  if (ref_levels > 0 && Mpi::Size(Mpi::World()) <= 2)
  {
    // A chunk holds the double-precision element matrices of 1188 p4 tetrahedra.
    int ne = mesh.GetNE();
    Mpi::GlobalMin(1, &ne, Mpi::World());
    CHECK(ne > 1188);
  }
  mfem::ND_FECollection fec(order, 3);
  FiniteElementSpace fespace(mesh, &fec);
  auto real_mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto imag_mass = subset ? BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 2, 1.0)
                          : BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto real_curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto bdr_mass = BuildCoefficient(mesh, true, CoeffType::Matrix);
  auto bdr_curl = BuildCoefficient(mesh, true, CoeffType::Scalar);
  real_mass *= -0.71;
  imag_mass *= imag_scale;
  real_curl *= 1.37;
  bdr_mass *= 0.23;
  if (subset)
  {
    std::vector<int> lossy_attr;
    imag_mass.AddNonzeroAttributes(lossy_attr);
    mesh.SetCeedActiveAttributes(std::move(lossy_attr), {});
  }
  auto Assemble = [&](bool imag)
  {
    BilinearForm a(fespace);
    if (subset)
    {
      a.SkipZeroCoefficientElements();
    }
    if (!imag)
    {
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(real_curl, real_mass);
      if (boundary)
      {
        a.AddBoundaryIntegrator<VectorFEMassIntegrator>(bdr_mass);
      }
    }
    else
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(imag_mass);
      if (boundary)
      {
        a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(bdr_curl, bdr_mass);
      }
    }
    a.AssembleQuadratureData();
    return a.PartialAssemble();
  };
  auto ref_real = Assemble(false), ref_imag = Assemble(true);
  ComplexWrapperOperator reference(ref_real.get(), ref_imag.get());
  // As for the multigrid levels, the original operator packs its leaves first.
  auto local = ceed::CreateComplexOperator(Assemble(false), Assemble(true), fespace);
  std::size_t bytes;
  int dropped;
  auto approx = ceed::CreateSinglePrecisionOperator(local->Real(), local->Imag(), fespace,
                                                    &bytes, &dropped);
  REQUIRE(approx);
  REQUIRE(bytes > 0);
  long long num_dropped = dropped;
  Mpi::GlobalSum(1, &num_dropped, Mpi::World());
  if (dropped_expected)
  {
    REQUIRE(num_dropped > 0);
  }
  else
  {
    REQUIRE(num_dropped == 0);
  }

  ComplexVector x(approx->Width()), y(approx->Height()), expected(approx->Height()),
      initial(approx->Height());
  x.Real().Randomize(17);
  x.Imag().Randomize(53);
  initial.Real().Randomize(29);
  initial.Imag().Randomize(97);
  reference.Mult(x, expected);
  approx->Mult(x, y);
  const double error = RelativeDifference(y, expected);
  CAPTURE(error);
  CHECK(error > single_min);
  CHECK(error < single_max);

  // AddMult with a complex coefficient.
  const std::complex<double> a{-0.7, 0.23};
  y = initial;
  approx->AddMult(x, y, a);
  ComplexVector target(initial);
  target.AXPY(a, expected);
  CHECK(RelativeDifference(y, target) < single_max);

  // The real part alone.
  auto approx_real = ceed::CreateSinglePrecisionOperator(local->Real(), nullptr, fespace);
  REQUIRE(approx_real);
  ComplexWrapperOperator reference_real(ref_real.get(), nullptr);
  reference_real.Mult(x, expected);
  approx_real->Mult(x, y);
  const double error_real = RelativeDifference(y, expected);
  CAPTURE(error_real);
  CHECK(error_real > single_min);
  CHECK(error_real < single_max);

  // The approximation follows in-place updates of the quadrature data.
  ScalePackedTestQData(dynamic_cast<const ceed::Operator &>(*local->Real()), -1.13);
  ScalePackedTestQData(*ref_real, -1.13);
  reference.Mult(x, expected);
  approx->Mult(x, y);
  CHECK(RelativeDifference(y, expected) < single_max);
}

TEST_CASE("libCEED single-precision approximation of parallel operators",
          "[libCEED][SinglePrecision][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Single-precision application requires one CPU CEED context and one thread");
  }
  const int order = GENERATE(2, 3);
  CAPTURE(order);
  PackedIntegrationSettings settings(order);
  auto mesh = Initialize(Mpi::World(), TestMeshPath("fichera-tet.mesh"), 0, false);
  auto mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto loss = BuildCoefficient(mesh, false, CoeffType::Matrix);
  loss *= 0.05;
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  auto Assemble = [&](bool imag)
  {
    BilinearForm a(nd_fespace);
    if (imag)
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(loss);
    }
    else
    {
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass);
    }
    a.AssembleQuadratureData();
    return a.PartialAssemble();
  };

  SECTION("Real H(curl) operator")
  {
    ParOperator reference(Assemble(false), nd_fespace), actual(Assemble(false), nd_fespace);
    REQUIRE(!reference.GetSinglePrecisionOperator());
    REQUIRE(actual.UseSinglePrecision() > 0);
    const auto *single = actual.GetSinglePrecisionOperator();
    REQUIRE(single);
    Vector x(actual.Width()), y(actual.Height()), expected(actual.Height());
    x.Randomize(7);
    reference.Mult(x, expected);
    // The operator itself stays exact.
    actual.Mult(x, y);
    CHECK(RelativeDifference(y, expected) < 1.0e-14);
    single->Mult(x, y);
    const double error = RelativeDifference(y, expected);
    CAPTURE(error);
    CHECK(error > single_min);
    CHECK(error < single_max);
    y = 1.0;
    single->AddMult(x, y, -0.5);
    expected *= -0.5;
    expected += 1.0;
    CHECK(RelativeDifference(y, expected) < single_max);
  }

  SECTION("Real H1 operator with element-assembled leaves")
  {
    // Element-assembled leaves are applied through their own (double-precision) element
    // matrices, which leaves nothing for a single-precision approximation.
    auto h1_mass = BuildCoefficient(mesh, false, CoeffType::Scalar);
    BilinearForm a(h1_fespace);
    a.AddDomainIntegrator<DiffusionMassIntegrator>(mass, h1_mass);
    ParOperator actual(a.PartialAssemble(), h1_fespace);
    CHECK(actual.UseSinglePrecision() == 0);
    CHECK(!actual.GetSinglePrecisionOperator());
  }

  SECTION("Complex operator and Chebyshev smoothers")
  {
    ComplexParOperator reference(Assemble(false), Assemble(true), nd_fespace),
        actual(Assemble(false), Assemble(true), nd_fespace);
    const auto &bdr_attr = mesh.Get().bdr_attributes;
    mfem::Array<int> bdr_marker(bdr_attr.Size() ? bdr_attr.Max() : 0), tdofs;
    bdr_marker = 0;
    if (bdr_marker.Size())
    {
      bdr_marker[0] = 1;
    }
    nd_fespace.GetEssentialTrueDofs(bdr_marker, tdofs);
    reference.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
    actual.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
    int dropped;
    REQUIRE(actual.UseSinglePrecision(&dropped) > 0);
    CHECK(dropped == 0);
    const auto *single = actual.GetSinglePrecisionOperator();
    REQUIRE(single);
    ComplexVector x(actual.Width()), y(actual.Height()), expected(actual.Height()),
        y0(actual.Height());
    x.Real().Randomize(11);
    x.Imag().Randomize(13);
    y0.Real().Randomize(19);
    y0.Imag().Randomize(23);
    reference.Mult(x, expected);
    single->Mult(x, y);
    const double error = RelativeDifference(y, expected);
    CAPTURE(error);
    CHECK(error > single_min);
    CHECK(error < single_max);

    // The smoothers apply the single-precision approximation to their search directions
    // only. Each application perturbs the polynomial by about the rounding of the element
    // matrices, so the smoothed results agree to about the polynomial order times that.
    for (bool fourth : {true, false})
    {
      CAPTURE(fourth);
      constexpr int smooth_order = 4;
      std::unique_ptr<Solver<ComplexOperator>> smoother_ref, smoother;
      for (auto *s : {&smoother_ref, &smoother})
      {
        if (fourth)
        {
          *s = std::make_unique<ChebyshevSmoother<ComplexOperator>>(Mpi::World(), 1,
                                                                    smooth_order, 1.0);
        }
        else
        {
          *s = std::make_unique<ChebyshevSmoother1stKind<ComplexOperator>>(
              Mpi::World(), 1, smooth_order, 1.0, 0.0);
        }
      }
      smoother_ref->SetOperator(reference);
      smoother->SetOperator(actual);
      if (fourth)
      {
        CHECK(!static_cast<ChebyshevSmoother<ComplexOperator> &>(*smoother_ref)
                   .UsesSinglePrecision());
        CHECK(static_cast<ChebyshevSmoother<ComplexOperator> &>(*smoother)
                  .UsesSinglePrecision());
      }
      for (bool guess : {false, true})
      {
        CAPTURE(guess);
        smoother_ref->SetInitialGuess(guess);
        smoother->SetInitialGuess(guess);
        expected = y0;
        y = y0;
        smoother_ref->Mult(x, expected);
        smoother->Mult(x, y);
        const double smooth_error = RelativeDifference(y, expected);
        CAPTURE(smooth_error);
        CHECK(smooth_error > single_min);
        CHECK(smooth_error < smooth_order * single_max);
      }
    }
  }
}

TEST_CASE("Smoother residual updates", "[libCEED][SinglePrecision][Serial][Parallel]")
{
  if (!PackedTestBackend())
  {
    SKIP("Single-precision application requires one CPU CEED context and one thread");
  }
  const int order = GENERATE(2, 3);
  const bool single = GENERATE(false, true);
  CAPTURE(order, single);
  PackedIntegrationSettings settings(order);
  auto mesh = Initialize(Mpi::World(), TestMeshPath("fichera-tet.mesh"), 0, false);
  auto mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
  auto loss = BuildCoefficient(mesh, false, CoeffType::Matrix);
  loss *= 0.05;
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  auto Assemble = [&](FiniteElementSpace &fespace, bool imag)
  {
    BilinearForm a(fespace);
    if (&fespace == &h1_fespace)
    {
      // The auxiliary space operator Gᵀ A G.
      a.AddDomainIntegrator<DiffusionIntegrator>(imag ? loss : mass);
    }
    else if (imag)
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(loss);
    }
    else
    {
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass);
    }
    a.AssembleQuadratureData();
    return a.PartialAssemble();
  };
  ComplexParOperator A(Assemble(nd_fespace, false), Assemble(nd_fespace, true), nd_fespace),
      A_G(Assemble(h1_fespace, false), Assemble(h1_fespace, true), h1_fespace);
  const auto &bdr_attr = mesh.Get().bdr_attributes;
  mfem::Array<int> bdr_marker(bdr_attr.Size() ? bdr_attr.Max() : 0), nd_tdofs, h1_tdofs;
  bdr_marker = 0;
  if (bdr_marker.Size())
  {
    bdr_marker[0] = 1;
  }
  nd_fespace.GetEssentialTrueDofs(bdr_marker, nd_tdofs);
  h1_fespace.GetEssentialTrueDofs(bdr_marker, h1_tdofs);
  A.SetEssentialTrueDofs(nd_tdofs, Operator::DIAG_ONE);
  A_G.SetEssentialTrueDofs(h1_tdofs, Operator::DIAG_ONE);
  if (single)
  {
    REQUIRE(A.UseSinglePrecision() > 0);
  }
  const auto &G = nd_fespace.GetDiscreteInterpolator(h1_fespace);

  // Smoothed results agree to rounding, which is single-precision rounding for the
  // single-precision applications (amplified by about the polynomial order). Updated
  // residuals carry the rounding errors of the updates, relative to the initial residual.
  constexpr int smooth_order = 4;
  const double tol = single ? smooth_order * single_max : 1.0e-11;
  const int n = A.Height();
  ComplexVector x(n), y0(n), y(n), r(n), y_ref(n), r_ref(n), z(n), phi(A_G.Height());
  x.Real().Randomize(11);
  x.Imag().Randomize(13);
  y0.Real().Randomize(19);
  y0.Imag().Randomize(23);
  phi.Real().Randomize(29);
  phi.Imag().Randomize(31);
  linalg::SetSubVector(phi, h1_tdofs, 0.0);
  G.Mult(phi.Real(), z.Real());
  G.Mult(phi.Imag(), z.Imag());
  auto Residual = [&](const ComplexVector &y, ComplexVector &res)
  {
    A.Mult(y, res);
    linalg::AXPBY(1.0, x, -1.0, res);
  };
  // Difference of an updated residual from the recomputed one, relative to the residual of
  // the initial guess.
  auto ResidualDifference =
      [&](const ComplexVector &res, const ComplexVector &res_ref, bool guess)
  {
    ComplexVector diff(res_ref), res0(n);
    diff.AXPY(-1.0, res);
    if (guess)
    {
      Residual(y0, res0);
    }
    else
    {
      res0 = x;
    }
    return linalg::Norml2(Mpi::World(), diff) / linalg::Norml2(Mpi::World(), res0);
  };

  SECTION("Chebyshev smoothers")
  {
    for (bool fourth : {true, false})
    {
      CAPTURE(fourth);
      std::unique_ptr<Solver<ComplexOperator>> smoother;
      if (fourth)
      {
        smoother = std::make_unique<ChebyshevSmoother<ComplexOperator>>(Mpi::World(), 1,
                                                                        smooth_order, 1.0);
      }
      else
      {
        smoother = std::make_unique<ChebyshevSmoother1stKind<ComplexOperator>>(
            Mpi::World(), 1, smooth_order, 1.0, 0.0);
      }
      smoother->SetOperator(A);
      REQUIRE(smoother->UpdatesResidual());
      for (bool guess : {false, true})
      {
        CAPTURE(guess);
        // The residual output leaves the result unchanged.
        smoother->SetInitialGuess(guess);
        y_ref = y0;
        smoother->SetResidual(false, false);
        smoother->Mult2(x, y_ref, r_ref);
        y = y0;
        smoother->SetResidual(false, true);
        smoother->Mult2(x, y, r);
        CHECK(RelativeDifference(y, y_ref) == 0.0);
        Residual(y, r_ref);
        const double res_error = ResidualDifference(r, r_ref, guess);
        CAPTURE(res_error);
        CHECK(res_error < tol);
      }

      // A residual input replaces the computed residual of the initial guess.
      smoother->SetInitialGuess(true);
      smoother->SetResidual(false, false);
      y_ref = y0;
      smoother->Mult2(x, y_ref, r_ref);
      y = y0;
      Residual(y, r);
      smoother->SetResidual(true, false);
      smoother->Mult2(x, y, r);
      CHECK(RelativeDifference(y, y_ref) == 0.0);

      // A pending correction is added to the iterate and to the updated residual.
      y = y0;
      Residual(y, r);
      smoother->SetResidual(true, true, &z);
      smoother->Mult2(x, y, r);
      Residual(y, r_ref);
      const double pending_error = ResidualDifference(r, r_ref, true);
      CAPTURE(pending_error);
      CHECK(pending_error < tol);
      smoother->SetResidual(false, false);
    }
  }

  SECTION("Distributive relaxation smoother")
  {
    auto Smoother = [&](MultigridAuxiliaryResidual residual)
    {
      auto smoother = std::make_unique<DistRelaxationSmoother<ComplexOperator>>(
          Mpi::World(), G, 1, 1, smooth_order, 1.0, 0.0, true, 0, -1, residual);
      smoother->SetOperators(A, A_G);
      return smoother;
    };
    auto exact = Smoother(MultigridAuxiliaryResidual::EXACT),
         updated = Smoother(MultigridAuxiliaryResidual::UPDATED),
         lagged = Smoother(MultigridAuxiliaryResidual::LAGGED);

    // The independent parts of the smoothers, for the lagged residuals.
    ChebyshevSmoother<ComplexOperator> C_k(Mpi::World(), 1, smooth_order, 1.0),
        C_km1(Mpi::World(), 1, smooth_order - 1, 1.0),
        B_G(Mpi::World(), 1, smooth_order, 1.0);
    C_k.SetOperator(A);
    C_km1.SetOperator(A);
    B_G.SetOperator(A_G);
    B_G.SetInitialGuess(false);
    ComplexVector x_G(A_G.Height()), y_G(A_G.Height()), y_k(n);
    auto AuxiliaryCorrection = [&](const ComplexVector &res, ComplexVector &corr)
    {
      G.MultTranspose(res.Real(), x_G.Real());
      G.MultTranspose(res.Imag(), x_G.Imag());
      linalg::SetSubVector(x_G, h1_tdofs, 0.0);
      B_G.Mult(x_G, y_G);
      G.Mult(y_G.Real(), corr.Real());
      G.Mult(y_G.Imag(), corr.Imag());
    };
    const double tol_lagged = single ? tol : 1.0e-10;

    for (bool transpose : {false, true})
    {
      for (bool guess : {false, true})
      {
        CAPTURE(transpose, guess);
        auto Apply = [&](const Solver<ComplexOperator> &smoother, ComplexVector &y)
        {
          y = y0;
          if (transpose)
          {
            smoother.MultTranspose2(x, y, r);
          }
          else
          {
            smoother.Mult2(x, y, r);
          }
        };
        exact->SetInitialGuess(guess);
        updated->SetInitialGuess(guess);
        lagged->SetInitialGuess(guess);
        Apply(*exact, y_ref);

        // The updated residuals reproduce the smoother.
        Apply(*updated, y);
        const double smooth_error = RelativeDifference(y, y_ref);
        CAPTURE(smooth_error);
        CHECK(smooth_error < tol);

        // The lagged residuals: Mult computes the auxiliary space correction from the
        // residual of the first k - 1 steps of the 4th-kind Chebyshev smoother C_k (those
        // of C_{k-1}), MultTranspose the first step of C_k from the residual before the
        // auxiliary space correction (both from the residual of the initial guess).
        Apply(*lagged, y);
        y_ref = y0;
        if (transpose)
        {
          if (guess)
          {
            Residual(y_ref, r_ref);
          }
          else
          {
            y_ref = 0.0;
            r_ref = x;
          }
          AuxiliaryCorrection(r_ref, z);
          C_k.SetResidual(true, false, &z);
          C_k.Mult2(x, y_ref, r_ref);
          C_k.SetResidual(false, false);
        }
        else
        {
          C_k.SetInitialGuess(guess);
          C_km1.SetInitialGuess(guess);
          y_k = y0;
          C_k.Mult(x, y_k);
          C_km1.Mult(x, y_ref);
          Residual(y_ref, r_ref);
          AuxiliaryCorrection(r_ref, z);
          y_ref = y_k;
          y_ref += z;
        }
        const double lagged_error = RelativeDifference(y, y_ref);
        CAPTURE(lagged_error);
        CHECK(lagged_error < tol_lagged);

        // The lagged residuals change the smoother, by much less than its correction here.
        Apply(*exact, y_ref);
        y_k = y_ref;
        if (guess)
        {
          y_k -= y0;
        }
        y -= y_ref;
        const double lagged_change =
            linalg::Norml2(Mpi::World(), y) / linalg::Norml2(Mpi::World(), y_k);
        CAPTURE(lagged_change);
        CHECK(lagged_change > tol);
        CHECK(lagged_change < 0.05);
      }
    }
  }

  SECTION("Exact residuals at low frequencies")
  {
    // With a mass part far below the single-precision rounding of the curl-curl part (low
    // frequencies, small k h), the residual updates through the single-precision operator
    // are inaccurate on the auxiliary space corrections, and the smoother keeps exact
    // residuals. Without single precision, it keeps the requested ones.
    auto mass_lf = mass, loss_lf = loss;
    mass_lf *= 1.0e-9;
    loss_lf *= 1.0e-9;
    auto AssembleLowFrequency = [&](FiniteElementSpace &fespace, bool imag)
    {
      BilinearForm a(fespace);
      if (&fespace == &h1_fespace)
      {
        a.AddDomainIntegrator<DiffusionIntegrator>(imag ? loss_lf : mass_lf);
      }
      else if (imag)
      {
        a.AddDomainIntegrator<VectorFEMassIntegrator>(loss_lf);
      }
      else
      {
        a.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass_lf);
      }
      a.AssembleQuadratureData();
      return a.PartialAssemble();
    };
    ComplexParOperator A_lf(AssembleLowFrequency(nd_fespace, false),
                            AssembleLowFrequency(nd_fespace, true), nd_fespace),
        A_G_lf(AssembleLowFrequency(h1_fespace, false),
               AssembleLowFrequency(h1_fespace, true), h1_fespace);
    A_lf.SetEssentialTrueDofs(nd_tdofs, Operator::DIAG_ONE);
    A_G_lf.SetEssentialTrueDofs(h1_tdofs, Operator::DIAG_ONE);
    if (single)
    {
      REQUIRE(A_lf.UseSinglePrecision() > 0);
    }
    for (auto residual :
         {MultigridAuxiliaryResidual::UPDATED, MultigridAuxiliaryResidual::LAGGED})
    {
      DistRelaxationSmoother<ComplexOperator> smoother(Mpi::World(), G, 1, 1, smooth_order,
                                                       1.0, 0.0, true, 0, -1, residual);
      smoother.SetOperators(A, A_G);
      const double error = smoother.GetGradientError();
      smoother.SetOperators(A_lf, A_G_lf);
      const double error_lf = smoother.GetGradientError();
      CAPTURE(error, error_lf);
      CHECK(smoother.GetAuxiliaryResidual() ==
            (single ? MultigridAuxiliaryResidual::EXACT : residual));
      if (single)
      {
        CHECK(error < 1.0e-3);
        CHECK(error_lf > 1.0);
      }
      else
      {
        CHECK(error == 0.0);
        CHECK(error_lf == 0.0);
      }

      // The residuals are resolved again for each operator.
      smoother.SetOperators(A, A_G);
      CHECK(smoother.GetAuxiliaryResidual() == residual);
    }
  }

  SECTION("Exact residuals at low frequencies with large boundary terms")
  {
    // A large term on gradients, such as the inductive term of an impedance boundary (which
    // does not depend on the frequency), dominates the norms of the operator on gradients,
    // so the single-precision error on the gradients away from it (the rounding of the
    // curl-curl part against the small mass part) is negligible in the global norms. In the
    // auxiliary space scaled by the diagonal of the auxiliary space operator it is not, and
    // the smoother computes the residual for the auxiliary space correction of Mult
    // exactly.
    auto mass_lf = mass, loss_lf = loss;
    mass_lf *= 1.0e-12;
    loss_lf *= 1.0e-12;
    auto bdr = BuildCoefficient(mesh, true, CoeffType::Matrix);
    bdr *= 1.0e2;
    auto AssembleBoundary = [&](FiniteElementSpace &fespace, bool imag)
    {
      BilinearForm a(fespace);
      if (&fespace == &h1_fespace)
      {
        a.AddDomainIntegrator<DiffusionIntegrator>(imag ? loss_lf : mass_lf);
        if (!imag)
        {
          a.AddBoundaryIntegrator<DiffusionIntegrator>(bdr);
        }
      }
      else if (imag)
      {
        a.AddDomainIntegrator<VectorFEMassIntegrator>(loss_lf);
      }
      else
      {
        a.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass_lf);
        a.AddBoundaryIntegrator<VectorFEMassIntegrator>(bdr);
      }
      a.AssembleQuadratureData();
      return a.PartialAssemble();
    };
    ComplexParOperator A_bdr(AssembleBoundary(nd_fespace, false),
                             AssembleBoundary(nd_fespace, true), nd_fespace),
        A_G_bdr(AssembleBoundary(h1_fespace, false), AssembleBoundary(h1_fespace, true),
                h1_fespace);
    A_bdr.SetEssentialTrueDofs(nd_tdofs, Operator::DIAG_ONE);
    A_G_bdr.SetEssentialTrueDofs(h1_tdofs, Operator::DIAG_ONE);
    if (single)
    {
      REQUIRE(A_bdr.UseSinglePrecision() > 0);
    }
    for (auto residual :
         {MultigridAuxiliaryResidual::UPDATED, MultigridAuxiliaryResidual::LAGGED})
    {
      DistRelaxationSmoother<ComplexOperator> smoother(Mpi::World(), G, 1, 1, smooth_order,
                                                       1.0, 0.0, true, 0, -1, residual);
      smoother.SetOperators(A_bdr, A_G_bdr);
      const double error = smoother.GetGradientError(),
                   error_aux = smoother.GetAuxiliaryGradientError();
      CAPTURE(error, error_aux);
      CHECK(smoother.GetAuxiliaryResidual() == residual);
      CHECK(smoother.GetExactAuxiliaryInput() == single);
      if (single)
      {
        // The global error alone would keep the requested residuals everywhere.
        CHECK(error < 1.0e-3);
        CHECK(error_aux > 100.0);
      }
      else
      {
        CHECK(error == 0.0);
        CHECK(error_aux == 0.0);
      }

      // Without the auxiliary space check (eigenmode simulations), the global error alone
      // keeps the requested residuals.
      {
        DistRelaxationSmoother<ComplexOperator> global(Mpi::World(), G, 1, 1, smooth_order,
                                                       1.0, 0.0, true, 0, -1, residual);
        global.SetAuxiliaryGradientCheck(false);
        global.SetOperators(A_bdr, A_G_bdr);
        CHECK(global.GetAuxiliaryResidual() == residual);
        CHECK(!global.GetExactAuxiliaryInput());
        CHECK(global.GetAuxiliaryGradientError() < 0.0);
        CHECK(global.GetGradientError() == error);
        global.SetAuxiliaryGradientCheck(true);
        CHECK(global.GetExactAuxiliaryInput() == single);
      }

      // Mult with the exact residual for the auxiliary space correction is the smoother
      // with exact residuals, up to the rounding of the single-precision applications.
      if (single)
      {
        DistRelaxationSmoother<ComplexOperator> exact(Mpi::World(), G, 1, 1, smooth_order,
                                                      1.0, 0.0, true, 0, -1,
                                                      MultigridAuxiliaryResidual::EXACT);
        exact.SetOperators(A_bdr, A_G_bdr);
        ComplexVector y1(n), y2(n), r1(n), r2(n);
        y1 = 0.0;
        y2 = 0.0;
        smoother.SetInitialGuess(false);
        exact.SetInitialGuess(false);
        smoother.Mult2(x, y1, r1);
        exact.Mult2(x, y2, r2);
        y1.AXPY(-1.0, y2);
        CHECK(linalg::Norml2(Mpi::World(), y1) <=
              1.0e-12 * linalg::Norml2(Mpi::World(), y2));
      }
    }
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
  const bool requested = std::ranges::any_of(
      selectors,
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

// Integrators whose coefficient support is contained by a configured shared support union
// omit the elements outside that union (BilinearForm::SkipZeroCoefficientElements). Those
// contributions are exactly 0.0, so the operator action, its transpose and its diagonal
// must match assembly over all elements to backend rounding (see TestVectorAgree).
// The largest difference of two distributed vectors relative to the largest entry of the
// second.
double MaxRelativeDifference(const Vector &x, const Vector &y, MPI_Comm comm)
{
  Vector diff(x);
  diff -= y;
  double norms[2] = {diff.Normlinf(), y.Normlinf()};
  Mpi::GlobalMax(2, norms, comm);
  return norms[0] / std::max(norms[1], std::numeric_limits<double>::min());
}

TEST_CASE("Exact diagonals of parallel operators",
          "[libCEED][ExactDiagonal][Serial][Parallel]")
{
  // ParOperator::AssembleDiagonal of partially assembled operators equals the diagonal of
  // their parallel assembly Pᵀ A P to rounding: on nonconforming meshes, whose prolongation
  // P has fractional constraint weights, and for the curl orientation transforms of
  // non-tensor H(curl) spaces of order p > 1, with boundary terms, essential boundary
  // conditions, complex operators and assembled local matrices. The approximation |P|ᵀ dₗ,
  // with the diagonal dₗ of libCEED, does not.
  const std::string mesh_file =
      GENERATE("fichera-tet.mesh", "fichera-hex.mesh", "fichera-mixed-p2.mesh");
  const bool amr = GENERATE(false, true);
  const int order = GENERATE(1, 2, 3, 4);
  CAPTURE(mesh_file, amr, order);
  PackedIntegrationSettings settings(order);
  MPI_Comm comm = Mpi::World();
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, amr);
  auto Q = BuildCoefficient(mesh, false, CoeffType::Scalar);
  auto MQ = BuildCoefficient(mesh, false, CoeffType::Matrix);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  const bool nd_approximate =
      amr || (order > 1 && !mfem::UsesTensorBasis(nd_fespace.Get()) &&
              mesh_file != "fichera-hex.mesh");
  constexpr double tol = 1.0e-12;

  // Essential true dofs on the first boundary attribute.
  auto EssentialTrueDofs = [&](const FiniteElementSpace &fespace)
  {
    const auto &bdr_attr = mesh.Get().bdr_attributes;
    mfem::Array<int> bdr_marker(bdr_attr.Size() ? bdr_attr.Max() : 0), tdofs;
    bdr_marker = 0;
    if (bdr_marker.Size())
    {
      bdr_marker[0] = 1;
    }
    fespace.GetEssentialTrueDofs(bdr_marker, tdofs);
    return tdofs;
  };

  // Diagonal of ParOperator, of its parallel assembly, and the approximation.
  auto Diagonals = [&](std::unique_ptr<Operator> &&A, std::unique_ptr<Operator> &&A_ref,
                       const FiniteElementSpace &fespace, Vector &d, Vector &d_ref,
                       Vector &d_approx)
  {
    auto tdofs = EssentialTrueDofs(fespace);
    ParOperator op(std::move(A), fespace), op_ref(std::move(A_ref), fespace);
    op.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
    op_ref.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
    d.SetSize(op.Height());
    d_ref.SetSize(op.Height());
    d_approx.SetSize(op.Height());
    d.UseDevice(true);
    d_ref.UseDevice(true);
    d_approx.UseDevice(true);
    op.AssembleDiagonal(d);
    Vector lx(fespace.GetVSize());
    lx.UseDevice(true);
    op.LocalOperator().AssembleDiagonal(lx);
    const Operator *P = fespace.GetProlongationMatrix();
    const auto *hP = dynamic_cast<const mfem::HypreParMatrix *>(P);
    if (const auto *fP = dynamic_cast<const hypre::HypreProlongation *>(P))
    {
      hP = &fP->GetMatrix();
    }
    if (hP)
    {
      hP->AbsMultTranspose(1.0, lx, 0.0, d_approx);
    }
    else
    {
      P->MultTranspose(lx, d_approx);
    }
    linalg::SetSubVector(d_approx, tdofs, 1.0);
    op_ref.ParallelAssemble().AssembleDiagonal(d_ref);
  };
  // The approximation is known to be inexact (approximate > 0), exact (approximate < 0) or
  // either (approximate = 0, as for the prolongations of parallel non-tensor H(curl) spaces
  // of order p > 1, which can have several entries per row for the dofs of shared faces).
  auto Check = [&](std::unique_ptr<Operator> &&A, std::unique_ptr<Operator> &&A_ref,
                   const FiniteElementSpace &fespace, int approximate)
  {
    Vector d, d_ref, d_approx;
    Diagonals(std::move(A), std::move(A_ref), fespace, d, d_ref, d_approx);
    const double error = MaxRelativeDifference(d, d_ref, comm),
                 error_approx = MaxRelativeDifference(d_approx, d_ref, comm);
    CAPTURE(error, error_approx);
    CHECK(error < tol);
    if (approximate > 0)
    {
      CHECK(error_approx > 1.0e-6);
    }
    else if (approximate < 0)
    {
      CHECK(error_approx < tol);
    }
  };

  SECTION("H(curl) curl-curl + mass with boundary mass")
  {
    auto Assemble = [&]()
    {
      BilinearForm a(nd_fespace);
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>();
      return a.PartialAssemble();
    };
    Check(Assemble(), Assemble(), nd_fespace, nd_approximate ? 1 : 0);
  }

  SECTION("H1 diffusion + mass with boundary mass")
  {
    auto Assemble = [&]()
    {
      BilinearForm a(h1_fespace);
      a.AddDomainIntegrator<DiffusionMassIntegrator>(MQ, Q);
      a.AddBoundaryIntegrator<MassIntegrator>();
      return a.PartialAssemble();
    };
    Check(Assemble(), Assemble(), h1_fespace, amr ? 1 : -1);
  }

  SECTION("Assembled local H(curl) matrix")
  {
    auto Assemble = [&]() -> std::unique_ptr<Operator>
    {
      BilinearForm a(nd_fespace);
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
      auto op = a.PartialAssemble();
      return BilinearForm::FullAssemble(*op, false);
    };
    // The diagonal of an assembled local matrix is exact, so only the prolongation makes a
    // difference.
    Check(Assemble(), Assemble(), nd_fespace, amr ? 1 : 0);
  }

  SECTION("Complex H(curl) operator")
  {
    auto AssembleReal = [&]()
    {
      BilinearForm a(nd_fespace);
      a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
      return a.PartialAssemble();
    };
    auto AssembleImag = [&]()
    {
      BilinearForm a(nd_fespace);
      a.AddDomainIntegrator<VectorFEMassIntegrator>(MQ);
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>();
      return a.PartialAssemble();
    };
    auto tdofs = EssentialTrueDofs(nd_fespace);
    ComplexParOperator op(AssembleReal(), AssembleImag(), nd_fespace);
    op.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
    ComplexVector d(op.Height());
    d.UseDevice(true);
    op.AssembleDiagonal(d);
    for (bool imag : {false, true})
    {
      ParOperator op_ref(imag ? AssembleImag() : AssembleReal(), nd_fespace);
      op_ref.SetEssentialTrueDofs(tdofs, imag ? Operator::DIAG_ZERO : Operator::DIAG_ONE);
      Vector d_ref(op.Height());
      d_ref.UseDevice(true);
      op_ref.ParallelAssemble().AssembleDiagonal(d_ref);
      const double error = MaxRelativeDifference(imag ? d.Imag() : d.Real(), d_ref, comm);
      CAPTURE(imag, error);
      CHECK(error < tol);
    }
  }
}

TEST_CASE("Restriction runs", "[libCEED][Prolongation][Serial][Parallel]")
{
  // FiniteElementSpace::GetRestrictionRuns gives the products with the restriction matrix
  // of the space (a unit entry per true dof) of mfem::SparseMatrix bitwise.
  const std::string mesh_file =
      GENERATE("fichera-tet.mesh", "fichera-hex.mesh", "fichera-mixed-p2.mesh");
  const bool amr = GENERATE(false, true);
  const int order = GENERATE(1, 3);
  CAPTURE(mesh_file, amr, order);
  MPI_Comm comm = Mpi::World();
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, amr);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  auto Bitwise = [comm](const Vector &u, const Vector &v)
  {
    int diff = (u.Size() != v.Size()) ||
               (u.Size() > 0 &&
                std::memcmp(u.HostRead(), v.HostRead(), u.Size() * sizeof(double)) != 0);
    Mpi::GlobalMax(1, &diff, comm);
    return diff == 0;
  };
  for (const auto *fespace : {&nd_fespace, &h1_fespace})
  {
    const auto *R = fespace->GetRestrictionRuns();
    int has = (R != nullptr);
    Mpi::GlobalMin(1, &has, comm);
    REQUIRE(has == 1);
    const auto *Rm = fespace->GetRestrictionMatrix();
    const int n = fespace->GetTrueVSize(), m = fespace->GetVSize();
    Vector lx(m), tx(n), ty(n), ty_ref(n), ly(m), ly_ref(m);
    lx.Randomize(1);
    tx.Randomize(2);
    Rm->Mult(lx, ty_ref);
    R->Mult(lx, ty);
    CHECK(Bitwise(ty, ty_ref));
    Rm->MultTranspose(tx, ly_ref);
    R->MultTranspose(tx, ly);
    CHECK(Bitwise(ly, ly_ref));
  }
}

TEST_CASE("Masked prolongation products", "[libCEED][Prolongation][Serial][Parallel]")
{
  // hypre::MaskedProlongation gives the products with the prolongation around an operator
  // with essential true dofs bitwise: P x̃ for x̃ = x with zeros at the essential true dofs,
  // and y + a ỹ for ỹ = Pᵀ x with zeros or given values there, for one and two vectors.
  // ParOperator and ComplexParOperator use them: AddMult equals y + a Mult bitwise.
  const std::string mesh_file =
      GENERATE("fichera-tet.mesh", "fichera-hex.mesh", "fichera-mixed-p2.mesh");
  const bool amr = GENERATE(false, true);
  const int order = GENERATE(1, 3);
  CAPTURE(mesh_file, amr, order);
  PackedIntegrationSettings settings(order);
  MPI_Comm comm = Mpi::World();
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, amr);
  mfem::ND_FECollection nd_fec(order, 3);
  mfem::H1_FECollection h1_fec(order, 3);
  FiniteElementSpace nd_fespace(mesh, &nd_fec), h1_fespace(mesh, &h1_fec);
  auto Bitwise = [comm](const Vector &u, const Vector &v)
  {
    int diff = (u.Size() != v.Size()) ||
               (u.Size() > 0 &&
                std::memcmp(u.HostRead(), v.HostRead(), u.Size() * sizeof(double)) != 0);
    Mpi::GlobalMax(1, &diff, comm);
    return diff == 0;
  };
  auto EssentialTrueDofs = [&](const FiniteElementSpace &fespace)
  {
    const auto &bdr_attr = mesh.Get().bdr_attributes;
    mfem::Array<int> bdr_marker(bdr_attr.Size() ? bdr_attr.Max() : 0), tdofs;
    bdr_marker = 0;
    if (bdr_marker.Size())
    {
      bdr_marker[0] = 1;
    }
    fespace.GetEssentialTrueDofs(bdr_marker, tdofs);
    return tdofs;
  };

  SECTION("Prolongation products")
  {
    for (const auto *fespace : {&nd_fespace, &h1_fespace})
    {
      const auto *hP =
          dynamic_cast<const hypre::HypreProlongation *>(fespace->GetProlongationMatrix());
      int has = (hP != nullptr);
      Mpi::GlobalMin(1, &has, comm);
      if (!has)
      {
        continue;
      }
      const auto tdofs = EssentialTrueDofs(*fespace);
      hypre::MaskedProlongation mP(*hP, tdofs);
      CHECK(mP.Matches(tdofs));
      const int n = fespace->GetTrueVSize(), m = fespace->GetVSize();
      Vector xr(n), xi(n), txr(n), txi(n), lr(m), li(m), lr2(m), li2(m);
      xr.Randomize(1);
      xi.Randomize(2);
      txr = xr;
      txi = xi;
      linalg::SetSubVector(txr, tdofs, 0.0);
      linalg::SetSubVector(txi, tdofs, 0.0);
      hP->Mult(txr, lr);
      mP.Mult(xr, lr2);
      CHECK(Bitwise(lr, lr2));
      hP->Mult(txr, txi, lr, li);
      mP.Mult(xr, xi, lr2, li2);
      CHECK(Bitwise(lr, lr2));
      CHECK(Bitwise(li, li2));
      lr.Randomize(3);
      li.Randomize(4);
      for (const double a : {-1.0, -0.7, 2.5})
      {
        for (const bool one : {false, true})
        {
          CAPTURE(a, one);
          Vector tyr(n), tyi(n), yr(n), yi(n), yr_ref(n), yi_ref(n);
          yr.Randomize(5);
          yi.Randomize(6);
          yr_ref = yr;
          yi_ref = yi;
          hP->MultTranspose(lr, tyr);
          linalg::SetSubVector(tyr, tdofs, 0.0);
          if (one)
          {
            linalg::SetSubVector(tyr, tdofs, xr);
          }
          yr_ref.Add(a, tyr);
          mP.AddMultTranspose(lr, yr, a, one ? &xr : nullptr);
          CHECK(Bitwise(yr, yr_ref));
          yr_ref = yr;
          hP->MultTranspose(lr, li, tyr, tyi);
          linalg::SetSubVector(tyr, tdofs, 0.0);
          linalg::SetSubVector(tyi, tdofs, 0.0);
          if (one)
          {
            linalg::SetSubVector(tyr, tdofs, xr);
            linalg::SetSubVector(tyi, tdofs, xi);
          }
          ComplexVector::AXPY(a, tyr, tyi, yr_ref, yi_ref);
          mP.AddMultTranspose(lr, li, yr, yi, a, one ? &xr : nullptr, one ? &xi : nullptr);
          CHECK(Bitwise(yr, yr_ref));
          CHECK(Bitwise(yi, yi_ref));
        }
      }
    }
  }

  SECTION("Operator applications")
  {
    auto Q = BuildCoefficient(mesh, false, CoeffType::Scalar);
    auto MQ = BuildCoefficient(mesh, false, CoeffType::Matrix);
    auto Assemble = [&](bool imag)
    {
      BilinearForm a(nd_fespace);
      if (imag)
      {
        a.AddDomainIntegrator<VectorFEMassIntegrator>(MQ);
      }
      else
      {
        a.AddDomainIntegrator<CurlCurlMassIntegrator>(MQ, Q);
      }
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>();
      return a.PartialAssemble();
    };
    const auto tdofs = EssentialTrueDofs(nd_fespace);
    const int n = nd_fespace.GetTrueVSize();
    for (const auto policy : {Operator::DIAG_ONE, Operator::DIAG_ZERO})
    {
      CAPTURE(policy);
      ParOperator op(Assemble(false), nd_fespace);
      op.SetEssentialTrueDofs(tdofs, policy);
      Vector x(n), y(n), y_ref(n), z(n);
      x.Randomize(1);
      y.Randomize(2);
      y_ref = y;
      op.Mult(x, z);
      y_ref.Add(-0.7, z);
      op.AddMult(x, y, -0.7);
      CHECK(Bitwise(y, y_ref));

      ComplexParOperator cop(Assemble(false), Assemble(true), nd_fespace);
      cop.SetEssentialTrueDofs(tdofs, policy);
      ComplexVector cx(n), cy(n), cy_ref(n), cz(n);
      cx.Real().Randomize(3);
      cx.Imag().Randomize(4);
      cy.Real().Randomize(5);
      cy.Imag().Randomize(6);
      cy_ref = cy;
      cop.Mult(cx, cz);
      cy_ref.AXPY(-0.7, cz);
      cop.AddMult(cx, cy, -0.7);
      CHECK(Bitwise(cy.Real(), cy_ref.Real()));
      CHECK(Bitwise(cy.Imag(), cy_ref.Imag()));
    }
  }

  SECTION("Mapped element applications")
  {
    if (!PackedTestBackend())
    {
      SKIP("Element-assembled application requires one CPU CEED context and one thread");
    }
    // Operators whose local action consists of element matrices are applied without
    // L-vectors (see ceed::AddMultMapped): the results equal the products with the
    // prolongation and the local operator on the L-vectors bitwise.
    auto mass = BuildCoefficient(mesh, false, CoeffType::Matrix);
    auto curl = BuildCoefficient(mesh, false, CoeffType::Matrix);
    auto loss = BuildCoefficient(mesh, false, CoeffType::Matrix);
    loss *= 0.05;
    auto AssembleND = [&](bool imag)
    {
      BilinearForm a(nd_fespace);
      if (imag)
      {
        a.AddDomainIntegrator<VectorFEMassIntegrator>(loss);
      }
      else
      {
        a.AddDomainIntegrator<CurlCurlMassIntegrator>(curl, mass);
      }
      a.AssembleQuadratureData();
      return a.PartialAssemble();
    };
    auto scalar = BuildCoefficient(mesh, false, CoeffType::Scalar);
    auto AssembleH1 = [&](bool imag)
    {
      BilinearForm a(h1_fespace);
      if (imag)
      {
        a.AddDomainIntegrator<MassIntegrator>(scalar);
      }
      else
      {
        a.AddDomainIntegrator<DiffusionMassIntegrator>(curl, scalar);
      }
      return a.PartialAssemble();
    };
    // y = Pᵀ A P x̃ (with x̃ = x zeroed and y set at the essential true dofs, as DIAG_ONE)
    // on the L-vectors.
    auto Explicit = [&](const FiniteElementSpace &fespace, const mfem::Array<int> &tdofs,
                        const Operator &local, const Vector &x, Vector &y)
    {
      const Operator *P = fespace.GetProlongationMatrix();
      Vector tx(x), lx(fespace.GetVSize()), ly(fespace.GetVSize());
      linalg::SetSubVector(tx, tdofs, 0.0);
      P->Mult(tx, lx);
      local.Mult(lx, ly);
      P->MultTranspose(ly, y);
      linalg::SetSubVector(y, tdofs, x);
    };
    auto ExplicitComplex = [&](const FiniteElementSpace &fespace,
                               const mfem::Array<int> &tdofs, const ComplexOperator &local,
                               const ComplexVector &x, ComplexVector &y)
    {
      const Operator *P = fespace.GetProlongationMatrix();
      ComplexVector tx(x), lx(fespace.GetVSize()), ly(fespace.GetVSize());
      linalg::SetSubVector(tx, tdofs, 0.0);
      P->Mult(tx.Real(), lx.Real());
      P->Mult(tx.Imag(), lx.Imag());
      local.Mult(lx, ly);
      P->MultTranspose(ly.Real(), y.Real());
      P->MultTranspose(ly.Imag(), y.Imag());
      linalg::SetSubVector(y, tdofs, x);
    };
    // Element assembly applies to the simplex leaves (single-precision H(curl), H1).
    for (const bool nd : {true, false})
    {
      if (mesh_file != "fichera-tet.mesh")
      {
        break;
      }
      CAPTURE(nd);
      const auto &fespace = nd ? nd_fespace : h1_fespace;
      const auto tdofs = EssentialTrueDofs(fespace);
      const int n = fespace.GetTrueVSize();

      // Real operator: the single-precision approximation (H(curl)) or the element-
      // assembled operator itself (H1).
      ParOperator op(nd ? AssembleND(false) : AssembleH1(false), fespace);
      op.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
      std::unique_ptr<Operator> local32;
      if (nd)
      {
        REQUIRE(op.UseSinglePrecision() > 0);
        local32 = ceed::CreateSinglePrecisionOperator(op.LocalOperator());
        REQUIRE(local32);
      }
      const Operator &local = nd ? *local32 : op.LocalOperator();
      const Operator &applied = nd ? *op.GetSinglePrecisionOperator() : op;
      int mapped = ceed::SupportsMapped(local);
      Mpi::GlobalMin(1, &mapped, comm);
      CHECK(mapped == 1);
      Vector x(n), y(n), y_ref(n), z(n);
      x.Randomize(1);
      Explicit(fespace, tdofs, local, x, y_ref);
      applied.Mult(x, y);
      CHECK(Bitwise(y, y_ref));
      z.Randomize(2);
      y = z;
      y_ref = z;
      Vector t(n);
      Explicit(fespace, tdofs, local, x, t);
      y_ref.Add(-0.7, t);
      applied.AddMult(x, y, -0.7);
      CHECK(Bitwise(y, y_ref));

      // Complex operator.
      ComplexParOperator cop(nd ? AssembleND(false) : AssembleH1(false),
                             nd ? AssembleND(true) : AssembleH1(true), fespace);
      cop.SetEssentialTrueDofs(tdofs, Operator::DIAG_ONE);
      std::unique_ptr<ComplexOperator> clocal32;
      if (nd)
      {
        REQUIRE(cop.UseSinglePrecision() > 0);
        clocal32 = ceed::CreateSinglePrecisionOperator(cop.LocalOperator().Real(),
                                                       cop.LocalOperator().Imag(), fespace);
        REQUIRE(clocal32);
      }
      const ComplexOperator &clocal = nd ? *clocal32 : cop.LocalOperator();
      const ComplexOperator &capplied = nd ? *cop.GetSinglePrecisionOperator() : cop;
      int cmapped = ceed::SupportsMapped(clocal);
      Mpi::GlobalMin(1, &cmapped, comm);
      CAPTURE(cmapped);
      ComplexVector cx(n), cy(n), cy_ref(n), cz(n), ct(n);
      cx.Real().Randomize(3);
      cx.Imag().Randomize(4);
      ExplicitComplex(fespace, tdofs, clocal, cx, cy_ref);
      capplied.Mult(cx, cy);
      CHECK(Bitwise(cy.Real(), cy_ref.Real()));
      CHECK(Bitwise(cy.Imag(), cy_ref.Imag()));
      cz.Real().Randomize(5);
      cz.Imag().Randomize(6);
      cy = cz;
      cy_ref = cz;
      ExplicitComplex(fespace, tdofs, clocal, cx, ct);
      cy_ref.AXPY(-0.7, ct);
      capplied.AddMult(cx, cy, -0.7);
      CHECK(Bitwise(cy.Real(), cy_ref.Real()));
      CHECK(Bitwise(cy.Imag(), cy_ref.Imag()));
      if (nd)
      {
        // The single-precision approximation of H(curl) operators is mapped.
        CHECK(cmapped == 1);
      }
    }
  }
}

TEST_CASE("libCEED Zero-Coefficient Element Skipping", "[libCEED][Serial][Parallel][GPU]")
{
  auto mesh_file = GENERATE("star-tri.mesh", "fichera-mixed-p2.mesh");
  const auto comm = MPI_COMM_WORLD;
  auto mesh = Initialize(comm, TestMeshPath(mesh_file), 0, false);
  const int dim = mesh.Dimension();
  constexpr int order = 2;

  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;

  // Element skipping only applies to partially assembled operators (the default).
  BilinearForm::pa_order_threshold = 1;

  INFO("Mesh: " << mesh_file);

  // Coefficients which are stored for every attribute but are exactly zero outside of every
  // active_stride-th one, as for a material property (such as a loss tangent) which
  // vanishes in part of the domain.
  auto Q = BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 2, 2.0);
  auto Q_mass = BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 3, 3.0);
  auto Q_zero = BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 0, 4.0);
  auto Q_bdr = BuildZeroOutsideCoefficient(mesh.MaxCeedBdrAttribute(), 3, 5.0);

  // Configure the two shared support unions as index-list subsets in mesh order.
  std::vector<int> lossy_attr, active_bdr_attr;
  Q.AddNonzeroAttributes(lossy_attr);
  Q_mass.AddNonzeroAttributes(lossy_attr);
  Q_bdr.AddNonzeroAttributes(active_bdr_attr);
  std::sort(lossy_attr.begin(), lossy_attr.end());
  lossy_attr.erase(std::unique(lossy_attr.begin(), lossy_attr.end()), lossy_attr.end());
  mesh.SetCeedActiveAttributes(std::move(lossy_attr), std::move(active_bdr_attr));
  TestCeedActiveGeometryRestriction(comm, mesh);

  mfem::ND_FECollection nd_fec(order, dim), nd_fec_coarse(order - 1, dim);
  FiniteElementSpace nd_fespace(mesh, &nd_fec);

  // One 2D path is enough to exercise the dimension- and GPU-independent subset contract;
  // the distinct boundary, coefficient-union, cache, and coarsening paths run in 3D below.
  if (dim == 2)
  {
    TestCeedSkipZeroCoefficientElements(
        comm, nd_fespace, false,
        [&](BilinearForm &a) { a.AddDomainIntegrator<VectorFEMassIntegrator>(Q); }, true);
    return;
  }
  SECTION("Boundary Mass")
  {
    TestCeedSkipZeroCoefficientElements(
        comm, nd_fespace, false,
        [&](BilinearForm &a) { a.AddBoundaryIntegrator<VectorFEMassIntegrator>(Q_bdr); });
  }
  SECTION("Domain Curl-Curl + Mass")
  {
    // The two coefficients have different zero attributes, so the assembled elements are
    // the union of their nonzero ones.
    TestCeedSkipZeroCoefficientElements(
        comm, nd_fespace, false,
        [&](BilinearForm &a) { a.AddDomainIntegrator<CurlCurlMassIntegrator>(Q, Q_mass); });
  }
  SECTION("Domain Mass with Cached Quadrature Data")
  {
    TestCeedSkipZeroCoefficientElements(
        comm, nd_fespace, true,
        [&](BilinearForm &a) { a.AddDomainIntegrator<VectorFEMassIntegrator>(Q); });
  }
  SECTION("Coefficient Zero on Every Attribute")
  {
    // No sub-operator is created at all: the operator is exactly zero, as it is with all
    // elements assembled.
    BilinearForm a_skip(nd_fespace);
    a_skip.SkipZeroCoefficientElements();
    a_skip.AddDomainIntegrator<VectorFEMassIntegrator>(Q_zero);
    auto op_skip = a_skip.PartialAssemble();

    // The composite operators must be empty (this is the whole-term skip of an exactly-zero
    // coefficient), not merely zero-valued.
    for (std::size_t i = 0; i < op_skip->Size(); i++)
    {
      Ceed ceed;
      CeedInt nsub_ops;
      PalaceCeedCallBackend(CeedOperatorGetCeed((*op_skip)[i], &ceed));
      PalaceCeedCall(ceed, CeedOperatorCompositeGetNumSub((*op_skip)[i], &nsub_ops));
      REQUIRE(nsub_ops == 0);
    }

    Vector x(op_skip->Width()), y(op_skip->Height()), d(op_skip->Height());
    x.UseDevice(true);
    y.UseDevice(true);
    d.UseDevice(true);
    x.Randomize(1);
    op_skip->Mult(x, y);
    op_skip->AssembleDiagonal(d);
    REQUIRE(y.Norml2() == 0.0);
    REQUIRE(d.Norml2() == 0.0);
  }
  SECTION("Multigrid Hierarchy")
  {
    // Levels which are not assembled directly reuse the quadrature data of another level
    // (ceed::CeedOperatorCoarsen), which must use the same element subset.
    FiniteElementSpaceHierarchy nd_fespaces(
        std::make_unique<FiniteElementSpace>(mesh, &nd_fec_coarse));
    nd_fespaces.AddLevel(std::make_unique<FiniteElementSpace>(mesh, &nd_fec));
    auto AddIntegrators = [&](BilinearForm &a)
    { a.AddDomainIntegrator<CurlCurlMassIntegrator>(Q, Q_mass); };
    constexpr bool skip_zeros = false;

    BilinearForm a_ref(nd_fespaces.GetFinestFESpace()),
        a_skip(nd_fespaces.GetFinestFESpace());
    AddIntegrators(a_ref);
    AddIntegrators(a_skip);
    a_skip.SkipZeroCoefficientElements();
    auto ops_ref = a_ref.Assemble(nd_fespaces, skip_zeros);
    auto ops_skip = a_skip.Assemble(nd_fespaces, skip_zeros);
    REQUIRE(ops_ref.size() == nd_fespaces.GetNumLevels());
    REQUIRE(ops_skip.size() == ops_ref.size());
    for (std::size_t l = 0; l < ops_ref.size(); l++)
    {
      TestCeedOperatorMultAgree(comm, *ops_skip[l], *ops_ref[l]);
    }
  }
}

// Mesh levels retained across adaptive refinement iterations keep their geometry data and
// are reconfigured with the same sets by every new SpaceOperator: a no-op. Changing the
// sets after geometry data exists is a programming error and is rejected. Kept apart from
// the element skipping test above, since the rejection aborts in builds of MFEM without
// exceptions.
TEST_CASE("libCEED Element Sets After Geometry Data", "[libCEED][Serial][Parallel][GPU]")
{
  const auto comm = MPI_COMM_WORLD;
  auto mesh = Initialize(comm, TestMeshPath("fichera-mixed-p2.mesh"), 0, false);
  constexpr int order = 2;
  fem::DefaultIntegrationOrder::p_trial = order;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  BilinearForm::pa_order_threshold = 1;

  auto Q = BuildZeroOutsideCoefficient(mesh.MaxCeedAttribute(), 2, 2.0);
  auto Q_bdr = BuildZeroOutsideCoefficient(mesh.MaxCeedBdrAttribute(), 3, 5.0);
  std::vector<int> lossy_attr, active_bdr_attr;
  Q.AddNonzeroAttributes(lossy_attr);
  Q_bdr.AddNonzeroAttributes(active_bdr_attr);
  std::sort(lossy_attr.begin(), lossy_attr.end());
  lossy_attr.erase(std::unique(lossy_attr.begin(), lossy_attr.end()), lossy_attr.end());
  mesh.SetCeedActiveAttributes(lossy_attr, active_bdr_attr);

  mfem::ND_FECollection nd_fec(order, mesh.Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  BilinearForm a(nd_fespace);
  a.SkipZeroCoefficientElements();
  a.AddDomainIntegrator<VectorFEMassIntegrator>(Q);
  auto op = a.PartialAssemble();
  REQUIRE(op);
  lossy_attr = mesh.GetCeedActiveAttributes(false);
  active_bdr_attr = mesh.GetCeedActiveAttributes(true);
  REQUIRE_NOTHROW(mesh.SetCeedActiveAttributes(lossy_attr, active_bdr_attr));
  REQUIRE(mesh.GetCeedActiveAttributes(false) == lossy_attr);
  REQUIRE(mesh.GetCeedActiveAttributes(true) == active_bdr_attr);
  // Any other valid set: drop the first attribute, or use attribute 1 if the set is empty.
  std::vector<int> other_attr(
      lossy_attr.empty() ? std::vector<int>{1}
                         : std::vector<int>(lossy_attr.begin() + 1, lossy_attr.end()));
  REQUIRE(other_attr != lossy_attr);
  REQUIRE_THROWS(mesh.SetCeedActiveAttributes(other_attr, active_bdr_attr));
  REQUIRE(mesh.GetCeedActiveAttributes(false) == lossy_attr);
}

TEST_CASE("3D libCEED Benchmarks", "[libCEED][Benchmark][Serial][Parallel]")
{
  auto mesh = GENERATE("fichera-hex.mesh", "fichera-tet.mesh");
  RunCeedBenchmarks(MPI_COMM_WORLD, TestMeshPath(mesh), benchmark_ref_levels, false,
                    benchmark_order);
}

}  // namespace palace
