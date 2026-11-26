// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numbers>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/romoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacepostoperator.hpp"
#include "test-helpers.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/tablecsv.hpp"

using namespace palace;
using namespace nlohmann;
using namespace Catch::Matchers;
using namespace std::complex_literals;

// Assemble operators as Eigen matrices for simpler testing.
template <typename T>
auto toEigenMatrix(const T &op, int n)
{
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);
  mfem::Vector w(n), v(n);
  w.UseDevice(true);
  v.UseDevice(true);
  for (int i = 0; i < n; i++)
  {
    w = 0.0;
    w.HostReadWrite()[i] = 1.0;
    op.Mult(w, v);
    const auto *V = v.HostRead();
    for (int j = 0; j < n; j++)
    {
      mat(j, i) = V[j];
    }
  }
  return mat;
}

class DrivenPostOperatorTest : public PostOperator<ProblemType::DRIVEN>
{
public:
  using PostOperator<ProblemType::DRIVEN>::PostOperator;
  const Measurement &GetMeasurement() const { return measurement_cache; }
};

class RomOperatorTest : public RomOperator
{
public:
  using RomOperator::AddAuxBlockDirections;
  using RomOperator::AddFrequencyDependentPermittivitySynthesis;
  using RomOperator::AddScalarDispersionSynthesis;
  using RomOperator::ApplyComplexPolynomialFitCorrections;
  using RomOperator::AugmentedPencil;
  using RomOperator::BoundaryTermKind;
  using RomOperator::BuildAugmentedPencil;
  using RomOperator::CalculateNormalizedPROMMatrices;
  using RomOperator::FitExactSplitDispersion;
  using RomOperator::FitScalarDispersion;
  using RomOperator::NormalizedMatrices;
  using RomOperator::RomOperator;
  using RomOperator::WavePortAuxBlock;
  using RomOperator::WavePortDispersionFit;
  using RomOperator::WavePortRegime;
  using RomOperator::wp_pairing_error;
  auto &GetWeightOp() const { return weight_op_W; }
  auto &GetOrthR() const { return orth_R; }
  auto &GetVectors() const { return V; }
  // Replace the basis only (not the projections), for the reduced postprocessing.
  void SetBasis(std::vector<Vector> &&basis) { V = std::move(basis); }
  auto &GetKr() const { return Kr; }
  auto &GetCr() const { return Cr; }
  auto &GetMr() const { return Mr; }
  auto &GetAr() const { return Ar; }
  void ResetPMLOnlineFit() { pml_online_fit.reset(); }
  bool HasPMLOnlineFit() const { return pml_online_fit.has_value(); }
  auto &GetKsp() { return *ksp; }
  const ComplexOperator *GetHDMOperator() const { return hdm_A.get(); }
  double GetHDMOmega() const { return hdm_omega; }
  using RomOperator::RefineHDM;
  using RomOperator::TruncatePROM;
  double GetLinearTol() const { return ksp->GetRelTol(); }
  void SetLinearTol(double tol) { ksp->SetRelTol(tol); }
  double GetLinearAbsTol() const { return ksp->GetAbsTol(); }
  void SetLinearAbsTol(double tol) { ksp->SetAbsTol(tol); }

  // ‖b - A u‖ / ‖b‖ for the system of the last HDM solve.
  double HDMResidual(const ComplexVector &u) const
  {
    REQUIRE(hdm_A);
    ComplexVector res(u.Size());
    res.UseDevice(true);
    hdm_A->Mult(u, res);
    linalg::AXPBY(1.0, hdm_b, -1.0, res);
    return linalg::Norml2(space_op.GetComm(), res) /
           linalg::Norml2(space_op.GetComm(), hdm_b);
  }

  // Copies of the basis labels and of everything projected onto the basis, as they are
  // (not trimmed to the basis size), and the number of vectors of the basis restricted to
  // each wave port.
  struct ProjectedState
  {
    std::size_t dim_V, V_wp_dim;
    std::map<int, std::size_t> V_wp_size;
    std::vector<std::string> v_node_label;
    Eigen::MatrixXd orth_R, gram_L;
    Eigen::MatrixXcd Kr, Cr, Mr, RHS1r, floquet_comp_V;
    std::vector<Eigen::MatrixXcd> bdr_terms_Ar;
  };
  ProjectedState GetProjectedState() const
  {
    ProjectedState state{.dim_V = V.size(),
                         .V_wp_dim = V_wp_dim,
                         .V_wp_size = {},
                         .v_node_label = v_node_label,
                         .orth_R = orth_R,
                         .gram_L = gram_L,
                         .Kr = Kr,
                         .Cr = Cr,
                         .Mr = Mr,
                         .RHS1r = RHS1r,
                         .floquet_comp_V = floquet_comp_V,
                         .bdr_terms_Ar = {}};
    for (const auto &[idx, vp] : V_wp)
    {
      state.V_wp_size[idx] = vp.size();
    }
    for (const auto &term : bdr_terms)
    {
      state.bdr_terms_Ar.push_back(term.Ar);
    }
    return state;
  }
  std::size_t NumBoundaryTerms() const { return bdr_terms.size(); }
  std::size_t NumFloquetOrders() const { return floquet_orders.size(); }
  // The projection Vᵀ op V of the factored term of a kind and index, which must exist.
  Eigen::MatrixXcd &GetBoundaryTermProjection(BoundaryTermKind kind, std::size_t idx)
  {
    auto *Ar = FindBoundaryTermProjection(kind, idx);
    MFEM_VERIFY(Ar, "Missing factored PROM term!");
    return *Ar;
  }
  // The same, or null without such a term.
  Eigen::MatrixXcd *FindBoundaryTermProjection(BoundaryTermKind kind, std::size_t idx)
  {
    const auto it = std::ranges::find_if(
        bdr_terms, [kind, idx](const BoundaryTerm &term)
        { return term.kind == kind && term.idx == static_cast<int>(idx); });
    return (it != bdr_terms.end()) ? &it->Ar : nullptr;
  }
  auto &GetFrequencyDependentPermittivityMassR(std::size_t m)
  {
    return GetBoundaryTermProjection(BoundaryTermKind::Permittivity, m);
  }
  // The first-order damping term next to frequency-dependent material k is a
  // frequency-dependent boundary mass term of SurfaceTerms, found by its label.
  int GetFarfieldDispersiveDampingTerm(std::size_t k) const
  {
    const auto &terms = space_op.GetSurfaceTerms();
    for (int i = 0; i < terms.GetNumExtraSystemBdrMassTerms(); i++)
    {
      if (terms.GetExtraSystemBdrMassTermLabel(i) ==
          fmt::format("farfield_permittivity_{}", k))
      {
        return i;
      }
    }
    MFEM_ABORT("Missing absorbing boundary term of SurfaceTerms!");
    return -1;
  }
  auto &GetFarfieldDispersiveMatrixR(std::size_t k, bool curl_curl)
  {
    return curl_curl ? GetBoundaryTermProjection(
                           BoundaryTermKind::FarfieldPermittivityCurlCurl, k)
                     : GetBoundaryTermProjection(
                           BoundaryTermKind::SurfaceTerm,
                           static_cast<std::size_t>(GetFarfieldDispersiveDampingTerm(k)));
  }
  // The scalar f(ω) of the term i·f(ω)·M of the absorbing boundary next to
  // frequency-dependent material k which circuit synthesis fits: the scale of its
  // SurfaceTerms term over i for the damping, ω r_k(ω), and (0.5/ω)/r_k(ω) for the
  // curl-curl term.
  std::function<std::complex<double>(std::complex<double>)>
  GetFarfieldDispersiveScalar(std::size_t k, bool curl_curl) const
  {
    if (curl_curl)
    {
      return [&farfield_op = space_op.GetSurfaceTerms().GetFarfieldOp(),
              k](std::complex<double> omega)
      { return 0.5 / (omega * farfield_op.GetDispersiveScale(k, omega)); };
    }
    return [&terms = space_op.GetSurfaceTerms(),
            i = GetFarfieldDispersiveDampingTerm(k)](std::complex<double> omega)
    {
      return terms.GetExtraSystemBdrMassTermScale(i, omega) /
             std::complex<double>(0.0, 1.0);
    };
  }
  // Null without a frequency-scaled Floquet wave vector.
  Eigen::MatrixXcd *GetFarfieldDispersiveFloquetMatrixR(std::size_t k, bool rotated)
  {
    return FindBoundaryTermProjection(
        rotated ? BoundaryTermKind::FarfieldPermittivityRotatedTrace
                : BoundaryTermKind::FarfieldPermittivityMixedCurl,
        k);
  }
  auto GetSweepBand() const { return std::make_pair(sweep_omega_min, sweep_omega_max); }
  double GetSynthesisTol() const { return waveport_synthesis_tol; }
  void ReserveBasis(std::size_t n) { V.reserve(n); }
};

namespace
{

template <typename MatType>
bool SameMatrix(const MatType &x, const MatType &y)
{
  return x.rows() == y.rows() && x.cols() == y.cols() && x == y;
}

// Relative difference in the Frobenius norm, of matrices of the same size.
template <typename MatType>
bool CloseMatrix(const MatType &x, const MatType &y, double rtol)
{
  return x.rows() == y.rows() && x.cols() == y.cols() && (x - y).norm() <= rtol * y.norm();
}

// Whether two states are the same, bit for bit (the state of a RomOperator and a copy of
// it saved before).
bool SameProjectedState(const RomOperatorTest::ProjectedState &a,
                        const RomOperatorTest::ProjectedState &b)
{
  if (a.bdr_terms_Ar.size() != b.bdr_terms_Ar.size())
  {
    return false;
  }
  for (std::size_t k = 0; k < a.bdr_terms_Ar.size(); k++)
  {
    if (!SameMatrix(a.bdr_terms_Ar[k], b.bdr_terms_Ar[k]))
    {
      return false;
    }
  }
  return a.dim_V == b.dim_V && a.V_wp_dim == b.V_wp_dim && a.V_wp_size == b.V_wp_size &&
         a.v_node_label == b.v_node_label && SameMatrix(a.orth_R, b.orth_R) &&
         SameMatrix(a.gram_L, b.gram_L) && SameMatrix(a.Kr, b.Kr) &&
         SameMatrix(a.Cr, b.Cr) && SameMatrix(a.Mr, b.Mr) && SameMatrix(a.RHS1r, b.RHS1r) &&
         SameMatrix(a.floquet_comp_V, b.floquet_comp_V);
}

// Whether the projections of the same vectors by two RomOperators agree: to round-off, not
// bit for bit, since the distributed (and, with OpenMP, threaded) reductions need not sum
// in the same order. The basis restricted to the wave ports is extended by SolvePROM, so it
// is compared through the PROM solutions. The labels are not compared.
bool CloseProjectedState(const RomOperatorTest::ProjectedState &a,
                         const RomOperatorTest::ProjectedState &b, double rtol = 1.0e-10)
{
  if (a.bdr_terms_Ar.size() != b.bdr_terms_Ar.size())
  {
    return false;
  }
  for (std::size_t k = 0; k < a.bdr_terms_Ar.size(); k++)
  {
    if (!CloseMatrix(a.bdr_terms_Ar[k], b.bdr_terms_Ar[k], rtol))
    {
      return false;
    }
  }
  return a.dim_V == b.dim_V && CloseMatrix(a.orth_R, b.orth_R, rtol) &&
         CloseMatrix(a.gram_L, b.gram_L, rtol) && CloseMatrix(a.Kr, b.Kr, rtol) &&
         CloseMatrix(a.Cr, b.Cr, rtol) && CloseMatrix(a.Mr, b.Mr, rtol) &&
         CloseMatrix(a.RHS1r, b.RHS1r, rtol) &&
         CloseMatrix(a.floquet_comp_V, b.floquet_comp_V, rtol);
}

}  // namespace

auto LoadScaleParMesh2(IoData &iodata, MPI_Comm world_comm)
{
  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    auto smesh = mesh::Load(iodata, world_comm);
    int ne = smesh ? smesh->GetNE() : 0;
    Mpi::GlobalMax(1, &ne, world_comm);
    if (Mpi::Size(world_comm) > ne)
    {
      SKIP("Mesh has fewer elements than processes");
    }
    if (iodata.model.Lc <= 0.0)
    {
      iodata.model.Lc = mesh::ComputeReferenceLength(smesh, world_comm);
    }
    iodata.NondimensionalizeInputs(smesh);
    mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), world_comm));
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  return mesh_;
}

TEST_CASE("MinimalRationalInterpolation", "[romoperator][Serial][Parallel]")
{
  auto fn_tan_shift = [](double z)
  { return std::tan(0.5 * std::numbers::pi * (z - std::complex<double>(1., 1.))); };

  // Test scalar case: 2 sample points for 4 x 2 vector
  MinimalRationalInterpolation mri_1;

  CHECK(mri_1.GetSamplePoints() == std::vector<double>{});
  CHECK_THROWS(mri_1.FindMaxError());

  for (double x_sample : {-2.0, -1.0, 1.0, 2.0})
  {
    auto tan_eval = fn_tan_shift(x_sample);

    // The rational vector function we want to evaluate. The vector has common poles all
    // determined by tan. The minimal-rational interpolation method finds coefficients using
    // SVD on the R matrix from a QR decomposition of samples. If the rank of the vector is
    // small or degenerate, this can have trivial solution and fail.
    //
    // Here we pick a 5 component vector where 2 pairs are degenerate (scalar multiples of
    // each other) which is sufficiently large but challenging a problem.
    Eigen::VectorXcd c_vec(5);
    c_vec << -tan_eval, 5. * tan_eval, x_sample * tan_eval, 10. * x_sample * tan_eval,
        2 * x_sample * x_sample * x_sample * tan_eval;

    mri_1.AddSolutionSample(x_sample, c_vec);
  }

  CHECK(mri_1.GetSamplePoints().size() == 4);
  CHECK(mri_1.GetSamplePoints() == std::vector<double>{-2.0, -1.0, 1.0, 2.0});

  // The same sample again is in the span of the previous ones.
  {
    const auto tan_eval = fn_tan_shift(2.0);
    Eigen::VectorXcd c_vec(5);
    c_vec << -tan_eval, 5. * tan_eval, 2.0 * tan_eval, 20. * tan_eval, 16. * tan_eval;
    CHECK_THROWS(mri_1.AddSolutionSample(2.0, c_vec));
  }

  // By symmetry of poles max error should be at zero.
  CHECK_THAT(mri_1.FindMaxError(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("MinimalRationalInterpolation-DistributedSearch",
          "[romoperator][Serial][Parallel]")
{
  // Samples of u(ω) = Σₖ aₖ / (ω - pₖ) + η g(ω) in C⁸ with three poles pₖ, one 10⁻⁴ above
  // the real axis between the samples, and a smooth g with η = 10⁻⁸, which keeps the
  // snapshot matrix of full numerical rank. The MRI is given the coefficients of each
  // sample in the orthonormal basis of the samples so far, which grows with each sample.
  // The maximum error is compared with a brute-force evaluation of its definition: the
  // minimum of |Q(z)| = |Σᵢ qᵢ / (z - zᵢ)| over the scanned grid of 10⁶ intervals, where q
  // is the right singular vector of the smallest singular value of the snapshot matrix with
  // columns [uᵢ; iωᵢ uᵢ]. The pole is at the last point of the first half of the grid,
  // which is the last point of a process when the scan is split over 2 or 4 processes. From
  // five samples on, the maximum error is at the pole.
  MPI_Comm comm = Mpi::World();
  const std::vector<double> z = {-2.0, 2.0, -1.3, 0.4, 1.0, -0.7, 1.6};
  constexpr int nr_sample = 1000000, n = 8;
  const double delta = (z[1] - z[0]) / nr_sample;
  const double z_pole = z[0] + 499999 * delta;
  const std::array<std::complex<double>, 3> p = {z_pole + 1.0e-4i, -0.9 + 0.8i, 1.2 - 0.6i};
  const auto S = static_cast<Eigen::Index>(z.size());
  Eigen::MatrixXcd U(2 * n, S);  // Snapshot matrix
  for (Eigen::Index l = 0; l < S; l++)
  {
    for (int r = 0; r < n; r++)
    {
      std::complex<double> u = 1.0e-8 * std::exp(std::complex<double>(0.3 * r, 1.0) * z[l]);
      for (std::size_t k = 0; k < p.size(); k++)
      {
        u += std::complex<double>(std::cos(1.0 + r + 2.0 * k), std::sin(0.5 * r - k)) /
             (z[l] - p[k]);
      }
      U(r, l) = u;
      U(n + r, l) = (1i * z[l]) * u;
    }
  }
  // The coefficients R of the samples U = Q R in the orthonormal basis Q.
  const Eigen::MatrixXcd R = Eigen::HouseholderQR<Eigen::MatrixXcd>(U.topRows(n))
                                 .matrixQR()
                                 .triangularView<Eigen::Upper>();
  MinimalRationalInterpolation mri;
  for (Eigen::Index m = 1; m <= S; m++)
  {
    mri.AddSolutionSample(z[m - 1], R.col(m - 1).head(m));
    if (m < 2)
    {
      continue;
    }
    CAPTURE(m);
    Eigen::JacobiSVD<Eigen::MatrixXcd, Eigen::ComputeFullV> svd(U.leftCols(m));
    const Eigen::VectorXcd q = svd.matrixV().col(m - 1);
    double Q_min = std::numeric_limits<double>::infinity(), z_min = z[0];
    for (int i = 0; i <= nr_sample; i++)
    {
      const double x = z[0] + i * delta;
      std::complex<double> Q = 0.0;
      for (Eigen::Index l = 0; l < m; l++)
      {
        Q += q(l) / (z[l] - x);
      }
      if (std::abs(Q) < Q_min)
      {
        Q_min = std::abs(Q);
        z_min = x;
      }
    }
    const double z_max = mri.FindMaxError(comm);
    CHECK(z_max == mri.FindMaxError());
    CHECK_THAT(z_max, WithinAbs(z_min, 0.5 * delta));
    if (m >= 5)
    {
      CHECK_THAT(z_max, WithinAbs(z_pole, 0.5 * delta));
    }
  }
}

TEST_CASE("MinimalRationalInterpolation-NarrowSpan", "[romoperator][Serial]")
{
  // Samples 2e-11 apart relative to z: delta = span / 1e6 is below the spacing of doubles
  // near z, so the scan must not step by accumulating delta.
  MinimalRationalInterpolation mri;
  const double z0 = 100.0, z1 = 100.0 + 2.0e-9;
  for (int i = 0; i < 2; i++)
  {
    Eigen::VectorXcd c = Eigen::VectorXcd::Zero(2);
    c(i) = 1.0;
    mri.AddSolutionSample(i == 0 ? z0 : z1, c);
  }
  const double z_max = mri.FindMaxError();
  CHECK(z_max >= z0);
  CHECK(z_max <= z1);
}

TEST_CASE("RomOperator factors the dispersive volume A2", "[romoperator][Serial][Parallel]")
{
  // The constructor checks the factored terms against the assembled A2 at the ends of the
  // sweep, omega = 1 and 2. With a conductivity boundary, the check at omega = 1 sees only
  // the boundary term (the volume scalar cancels there), so a missing factored volume term
  // would go undetected there and only the check at omega = 2 catches it. Without it, the
  // volume term is the only factored contributor. An absorbing boundary next to the
  // material depends on ε(iω), which is ε∞ at omega = 1, through its first- and
  // second-order terms.
  enum class Boundary
  {
    None,
    Conductivity,
    Absorbing1,
    Absorbing2
  };
  const auto boundary = GENERATE(Boundary::None, Boundary::Conductivity,
                                 Boundary::Absorbing1, Boundary::Absorbing2);
  CAPTURE(static_cast<int>(boundary));
  const bool with_boundary = boundary == Boundary::Conductivity;
  json config;
  config["Problem"] = {{"Type", "Driven"}, {"Output", "test_output"}};
  config["Model"] = {{"Mesh", "test.msh"}};
  // The full volume A2 cancels at s = i but not s = 2i.
  const json terms = {{{"Type", "PoleResidue"}, {"Pole", -1.0}, {"Residue", 4.0}},
                      {{"Type", "PoleResidue"}, {"Pole", -2.0}, {"Residue", -20.0}},
                      {{"Type", "PoleResidue"}, {"Pole", -3.0}, {"Residue", 20.0}}};
  config["Domains"]["Materials"] = {
      {{"Attributes", {1}},
       {"Permittivity", {{"HighFrequencyLimit", 1.0}, {"Terms", terms}}}}};
  config["Boundaries"] = json::object();
  if (with_boundary)
  {
    config["Boundaries"]["Conductivity"] = {{{"Attributes", {1}}, {"Conductivity", 5.8e7}}};
  }
  if (boundary == Boundary::Absorbing1 || boundary == Boundary::Absorbing2)
  {
    config["Boundaries"]["Absorbing"] = {
        {"Attributes", {1}}, {"Order", boundary == Boundary::Absorbing1 ? 1 : 2}};
  }
  config["Solver"] = {
      {"Order", 1},
      {"Driven", {{"Samples", {{{"MinFreq", 1.0}, {"MaxFreq", 2.0}, {"FreqStep", 1.0}}}}}}};
  IoData iodata(config, false);
  test::IntegrationSettingsGuard settings_guard(iodata.solver);

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  SpaceOperator space_op(iodata.solver, iodata.domains, iodata.boundaries,
                         iodata.problem.type, iodata.units, mesh);

  // The volume term vanishes at omega = 1, the first frequency of the check of the factored
  // terms against the assembled A2 in the constructor. It must still be factored.
  RomOperatorTest prom_op(iodata, space_op, 1);
  CHECK(prom_op.FindBoundaryTermProjection(RomOperatorTest::BoundaryTermKind::Permittivity,
                                           0));
  CHECK(space_op.GetSurfaceTerms().GetSurfaceConductivityOp().Size() ==
        (with_boundary ? 1 : 0));
  CHECK(space_op.GetSurfaceTerms().GetFarfieldOp().NumDispersiveMaterials() ==
        (boundary == Boundary::Absorbing1 || boundary == Boundary::Absorbing2 ? 1 : 0));
  auto A2_cancelled =
      space_op.GetExtraSystemMatrix<ComplexOperator>(1.0, Operator::DIAG_ZERO, false);
  REQUIRE(A2_cancelled);

  prom_op.ReserveBasis(1);
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(Mpi::World(), v.Real());
  prom_op.UpdatePROM(v, "test");
  REQUIRE(prom_op.GetReducedDimension() == 1);

  prom_op.SolvePROM(1.0);  // Exercise the cancelling first frequency first.
  constexpr double omega = 2.0;
  prom_op.SolvePROM(omega);
  auto A2 =
      space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO, false);
  REQUIRE(A2);
  ComplexVector basis(v.Size()), A2v(v.Size());
  basis.UseDevice(true);
  A2v.UseDevice(true);
  basis = 0.0;
  basis.Real() = prom_op.GetVectors()[0];
  A2->Mult(basis, A2v);
  double a2_reduced_real = prom_op.GetVectors()[0] * A2v.Real();
  double a2_reduced_imag = prom_op.GetVectors()[0] * A2v.Imag();
  Mpi::GlobalSum(1, &a2_reduced_real, Mpi::World());
  Mpi::GlobalSum(1, &a2_reduced_imag, Mpi::World());
  const std::complex<double> a2_reduced{a2_reduced_real, a2_reduced_imag};
  std::complex<double> recovered =
      prom_op.GetAr()(0, 0) - prom_op.GetKr()(0, 0) + omega * omega * prom_op.GetMr()(0, 0);
  if (prom_op.GetCr().size() > 0)
  {
    recovered -= std::complex<double>(0.0, omega) * prom_op.GetCr()(0, 0);
  }
  CHECK_THAT(std::abs(recovered - a2_reduced),
             Catch::Matchers::WithinAbsMatcher(0.0, 1e-10));
}

// A scaled Floquet wave vector moves the dispersive ABC rotated trace and mixed curl
// into A2, with scalars ω/r and 1/r. Compare the actual online matrix with a full A2
// projection at several frequencies; one matching frequency cannot prove completeness.
TEST_CASE("RomOperator factors dispersive absorbing boundaries with a Floquet wave vector",
          "[romoperator][Serial][Parallel][GPU]")
{
  const bool scaled = GENERATE(false, true);
  CAPTURE(scaled);
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Output", "test_output"}};
  setup_json["Model"] = {{"Mesh", "test.msh"}};
  setup_json["Domains"]["Materials"] = {{{"Attributes", {1}}, {"Permittivity", 2.0}},
                                        {{"Attributes", {2}}, {"Permittivity", 1.0}}};
  setup_json["Boundaries"]["Absorbing"] = {{"Attributes", {3, 5}}, {"Order", 2}};
  setup_json["Solver"] = {
      {"Order", 2},
      {"Driven",
       {{"Samples", {{{"MinFreq", 0.8}, {"MaxFreq", 2.7}, {"FreqStep", 0.95}}}}}}};
  IoData iodata(setup_json, false);
  test::IntegrationSettingsGuard settings_guard(iodata.solver);
  // Direct operator fixtures use normalized parameters, as in the material A2 tests.
  iodata.domains.materials[0].permittivity_terms = {config::LorentzData{0.7, 3.0, 0.4}};
  iodata.domains.materials[0].sigma.s = {0.2, 0.2, 0.2};
  iodata.boundaries.periodic.wave_vector = {0.3, 0.4, 0.2};
  iodata.boundaries.periodic.floquet_reference_freq = scaled ? 1.3 : 0.0;

  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 1, 1, mfem::Element::TETRAHEDRON, 2.0, 1.0, 1.0));
  for (int e = 0; e < serial_mesh->GetNE(); e++)
  {
    mfem::Vector center;
    serial_mesh->GetElementCenter(e, center);
    serial_mesh->SetAttribute(e, center(0) < 1.0 ? 1 : 2);
  }
  serial_mesh->SetAttributes();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));
  SpaceOperator space_op(iodata.solver, iodata.domains, iodata.boundaries,
                         iodata.problem.type, iodata.units, mesh);
  REQUIRE(space_op.GetMaterialOp().HasFloquetFrequencyScaling() == scaled);
  REQUIRE(space_op.GetSurfaceTerms().GetFarfieldOp().NumDispersiveMaterials() == 1);
  RomOperatorTest prom_op(iodata, space_op, 3);
  prom_op.ReserveBasis(3);
  const auto comm = space_op.GetComm();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  for (int j = 0; j < 3; j++)
  {
    v = 0.0;
    linalg::SetRandom(comm, v.Real(), 1 + j + 3 * Mpi::Rank(comm));
    prom_op.UpdatePROM(v, "test");
  }
  REQUIRE(prom_op.GetReducedDimension() == 3);
  ComplexVector x(v.Size()), y(v.Size());
  x.UseDevice(true);
  y.UseDevice(true);
  constexpr double tol = 1.0e-10;  // Same full-projection tolerance as the volume A2 test.
  for (double omega : {0.8, 1.6, 2.7})
  {
    CAPTURE(omega);
    auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO,
                                                             /*include_wave_ports=*/false);
    REQUIRE(A2);
    Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(3, 3);
    for (int j = 0; j < 3; j++)
    {
      x = 0.0;
      x.Real() = prom_op.GetVectors()[j];
      A2->Mult(x, y);
      for (int i = 0; i < 3; i++)
      {
        full(i, j) = {linalg::Dot(comm, prom_op.GetVectors()[i], y.Real()),
                      linalg::Dot(comm, prom_op.GetVectors()[i], y.Imag())};
      }
    }
    prom_op.SolvePROM(omega);
    const Eigen::MatrixXcd complete = prom_op.GetAr();
    Eigen::MatrixXcd recovered =
        complete - prom_op.GetKr() + omega * omega * prom_op.GetMr();
    if (prom_op.GetCr().size() > 0)
    {
      recovered -= std::complex<double>(0.0, omega) * prom_op.GetCr();
    }
    CHECK((recovered - full).norm() <= tol);
    // Fixed k puts the mixed curl inside the full Hermitian ABC form; scaled k puts
    // it in its own real antisymmetric factor. Neither may be mirrored as symmetric.
    CHECK((full - full.transpose()).norm() > tol);
    for (bool rotated : {true, false})
    {
      CAPTURE(rotated);
      auto *projected = prom_op.GetFarfieldDispersiveFloquetMatrixR(0, rotated);
      if (!scaled)
      {
        CHECK(!projected);
        continue;
      }
      REQUIRE(projected);
      REQUIRE(projected->rows() == 3);
      const Eigen::MatrixXcd saved = *projected;
      const auto r =
          space_op.GetSurfaceTerms().GetFarfieldOp().GetDispersiveScale(0, omega);
      const auto scalar = (rotated ? omega : 1.0) / r;
      const Eigen::MatrixXcd omitted = scalar * saved;
      REQUIRE(omitted.norm() > tol);
      // If SolvePROM still used full-mesh assembly, zeroing a reduced factor would
      // have no effect. Each missing term must instead change its online matrix by
      // exactly that factor, proving both the route and its distinct ω power.
      projected->setZero();
      prom_op.SolvePROM(omega);
      CHECK((prom_op.GetAr() - (complete - omitted)).norm() <= tol);
      CHECK((prom_op.GetAr() - complete).norm() > tol);
      *projected = saved;
    }
  }
}

// Circuit synthesis must realize the factored dispersive volume term: the polynomial
// corrections plus the aux states of the augmented pencil must reproduce g(iω)·VᵀBV after
// eliminating the aux states. Pole terms are exact, also for a pole far from the band,
// where splitting r s²/(s - p) into a polynomial and a proper part cancels
// catastrophically, and for a projected mass with a direction weaker than the aux rank
// truncation, which must drop that direction's whole contribution rather than leave its
// polynomial part. A Djordjevic-Sarkar term is fit.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "RomOperator synthesis realizes the dispersive volume A2",
                 "[romoperator][Serial][Parallel]")
{
  enum class Case
  {
    Moderate,
    DjordjevicSarkar,
    DistantPole,
    WeakDirection
  };
  const auto test_case = GENERATE(Case::Moderate, Case::DjordjevicSarkar, Case::DistantPole,
                                  Case::WeakDirection);
  CAPTURE(static_cast<int>(test_case));

  // Relaxation times in ns, against the nondimensional time scale Lc / c₀ ≈ 2.3e-5 ns.
  const double tau = (test_case == Case::DistantPole)     ? 1.0e-12
                     : (test_case == Case::WeakDirection) ? 1.0e-4
                                                          : 0.01;
  json terms = {{{"Type", "Debye"}, {"DeltaPermittivity", 1.5}, {"RelaxationTime", tau}},
                {{"Type", "Lorentz"},
                 {"DeltaPermittivity", 0.5},
                 {"ResonanceFrequency", 20.0},
                 {"DampingFrequency", 2.0}}};
  if (test_case == Case::DjordjevicSarkar)
  {
    terms.push_back({{"Type", "DjordjevicSarkar"},
                     {"Strength", 0.2},
                     {"LowerFrequency", 1.0},
                     {"UpperFrequency", 100.0}});
  }
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
      {"L0", 1.0e-6},
      {"Lc", 7.0},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       {{{"Attributes", {1}},
         {"Permittivity", {{"HighFrequencyLimit", 2.0}, {"Terms", terms}}}}}}};
  setup_json["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 50.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}}};
  setup_json["Solver"] = {{"Order", 2},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-6},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};
  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
  SpaceOperator space_op(iodata, mesh_io);

  RomOperatorTest prom_op(iodata, space_op, 10);
  prom_op.AddLumpedPortModesForSynthesis();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(Mpi::World(), v.Real());
  prom_op.UpdatePROM(v, "random");
  const long n = prom_op.GetReducedDimension();
  REQUIRE(n == 2);

  // B_r = i·VᵀBV by the imaginary-slot convention. The weak-direction case replaces the
  // projection by a rotated diag(1, 1e-7), below the aux rank truncation of 1e-6.
  Eigen::MatrixXcd B_r = prom_op.GetFrequencyDependentPermittivityMassR(0);
  REQUIRE(B_r.rows() == n);
  if (test_case == Case::WeakDirection)
  {
    const double c = std::cos(0.3), s = std::sin(0.3);
    Eigen::Matrix2d R, D = Eigen::Vector2d(1.0, 1.0e-7).asDiagonal();
    R << c, -s, s, c;
    B_r = std::complex<double>(0.0, 1.0) *
          (R * D * R.transpose()).cast<std::complex<double>>();
  }

  Eigen::MatrixXcd Kc = Eigen::MatrixXcd::Zero(n, n), Cc = Kc, Mc = Kc;
  std::vector<RomOperatorTest::WavePortAuxBlock> blocks;
  prom_op.AddFrequencyDependentPermittivitySynthesis("permittivity_0", B_r, 0, Kc, Cc, Mc,
                                                     blocks);
  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kc, Cc, Mc, blocks, aux_labels);
  const long n_aux = aug.Kr.rows() - n;
  REQUIRE(n_aux > 0);

  const Eigen::MatrixXcd VtBV = B_r / std::complex<double>(0.0, 1.0);
  const auto &mat_op = space_op.GetMaterialOp();
  const auto [w_min, w_max] = prom_op.GetSweepBand();
  REQUIRE(w_max > w_min);
  double max_err = 0.0, max_ref = 0.0;
  for (double t : {0.0, 0.3, 0.55, 1.0})
  {
    const double w = w_min + t * (w_max - w_min);
    const Eigen::MatrixXcd A =
        aug.Kr + std::complex<double>(0.0, w) * aug.Cr - w * w * aug.Mr;
    const Eigen::MatrixXcd schur =
        A.topLeftCorner(n, n) -
        A.topRightCorner(n, n_aux) * A.bottomRightCorner(n_aux, n_aux)
                                         .fullPivLu()
                                         .solve(A.bottomLeftCorner(n_aux, n));
    const Eigen::MatrixXcd want =
        mat_op.EvaluateFrequencyDependentPermittivityA2(0, std::complex<double>(0.0, w)) *
        VtBV;
    max_err = std::max(max_err, (schur - want).cwiseAbs().maxCoeff());
    max_ref = std::max(max_ref, want.cwiseAbs().maxCoeff());
  }
  const double tol = (test_case == Case::DjordjevicSarkar) ? 1.0e-5
                     : (test_case == Case::WeakDirection)  ? 1.0e-6
                                                           : 1.0e-10;
  CHECK(max_err <= tol * max_ref);

  // The synthesized circuit carries exactly these aux states.
  if (test_case != Case::WeakDirection)
  {
    const auto norm = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
    const auto n_permittivity_aux =
        std::ranges::count_if(norm.aux_labels, [](const std::string &l)
                              { return l.starts_with("permittivity_0"); });
    CHECK(n_permittivity_aux == n_aux);
  }
}

// Circuit synthesis must realize the absorbing boundary next to a frequency-dependent
// material: after eliminating the aux states, the fit of each of its terms reproduces
// f(ω)·VᵀMV over the band, with f = ω·r(ω) for the damping and (0.5/ω)/r(ω) for the
// curl-curl matrix, r = √(ε(iω)/ε∞).
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "RomOperator synthesis realizes the dispersive absorbing boundary",
                 "[romoperator][Serial][Parallel]")
{
  const json terms = {
      {{"Type", "Debye"}, {"DeltaPermittivity", 1.5}, {"RelaxationTime", 0.01}},
      {{"Type", "Lorentz"},
       {"DeltaPermittivity", 0.5},
       {"ResonanceFrequency", 20.0},
       {"DampingFrequency", 2.0}}};
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
      {"L0", 1.0e-6},
      {"Lc", 7.0},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       {{{"Attributes", {1}},
         {"Permittivity", {{"HighFrequencyLimit", 2.0}, {"Terms", terms}}}}}}};
  setup_json["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 50.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}},
                              {"Absorbing", {{"Attributes", {101}}, {"Order", 2}}}};
  setup_json["Solver"] = {{"Order", 2},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-6},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};
  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
  SpaceOperator space_op(iodata, mesh_io);
  const auto &farfield_op = space_op.GetSurfaceTerms().GetFarfieldOp();
  REQUIRE(farfield_op.NumDispersiveMaterials() == 1);

  RomOperatorTest prom_op(iodata, space_op, 10);
  prom_op.AddLumpedPortModesForSynthesis();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(Mpi::World(), v.Real());
  prom_op.UpdatePROM(v, "random");
  const long n = prom_op.GetReducedDimension();
  REQUIRE(n == 2);
  const auto [w_min, w_max] = prom_op.GetSweepBand();
  REQUIRE(w_max > w_min);

  for (const bool curl_curl : {false, true})
  {
    CAPTURE(curl_curl);
    const Eigen::MatrixXcd &M_r = prom_op.GetFarfieldDispersiveMatrixR(0, curl_curl);
    REQUIRE(M_r.rows() == n);
    REQUIRE(M_r.cwiseAbs().maxCoeff() > 0.0);
    Eigen::MatrixXcd Kc = Eigen::MatrixXcd::Zero(n, n), Cc = Kc, Mc = Kc;
    std::vector<RomOperatorTest::WavePortAuxBlock> blocks;
    prom_op.AddScalarDispersionSynthesis("farfield", M_r,
                                         prom_op.GetFarfieldDispersiveScalar(0, curl_curl),
                                         Kc, Cc, Mc, blocks);
    std::vector<std::string> aux_labels;
    auto aug = RomOperatorTest::BuildAugmentedPencil(Kc, Cc, Mc, blocks, aux_labels);
    const long n_aux = aug.Kr.rows() - n;
    double max_err = 0.0, max_ref = 0.0;
    for (double t : {0.0, 0.3, 0.55, 1.0})
    {
      const double w = w_min + t * (w_max - w_min);
      const Eigen::MatrixXcd A =
          aug.Kr + std::complex<double>(0.0, w) * aug.Cr - w * w * aug.Mr;
      Eigen::MatrixXcd schur = A.topLeftCorner(n, n);
      if (n_aux > 0)
      {
        schur -= A.topRightCorner(n, n_aux) * A.bottomRightCorner(n_aux, n_aux)
                                                  .fullPivLu()
                                                  .solve(A.bottomLeftCorner(n_aux, n));
      }
      const auto r = farfield_op.GetDispersiveScale(0, w);
      const Eigen::MatrixXcd want = (curl_curl ? 0.5 / (w * r) : w * r) * M_r;
      max_err = std::max(max_err, (schur - want).cwiseAbs().maxCoeff());
      max_ref = std::max(max_ref, want.cwiseAbs().maxCoeff());
    }
    CHECK(max_err <= 1.0e-5 * max_ref);
  }
}

// A first-order absorbing boundary next to a Lorentz material with its resonance inside the
// band: the fit of f(ω) = ω r(ω) has poles, realized on the directions of the projected
// mass above the rank cutoff. Here the mass is a rotated diag(1, 1e-7), whose weak
// direction is below the cutoff. The polynomial part of the fit must be realized on the
// same directions: on the full mass, the weak direction carries the polynomial part alone,
// the complement of the poles, which near the resonance is orders of magnitude larger than
// f, and the realized error exceeds the reported fit residual.
TEST_CASE_METHOD(
    palace::test::PerRankTempDir,
    "RomOperator synthesis of a dispersive absorbing boundary on a weak direction",
    "[romoperator][Serial][Parallel]")
{
  const double damping = GENERATE(0.4, 0.1);  // GHz
  CAPTURE(damping);
  const json terms = {{{"Type", "Lorentz"},
                       {"DeltaPermittivity", 2.0},
                       {"ResonanceFrequency", 20.0},
                       {"DampingFrequency", damping}}};
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
      {"L0", 1.0e-6},
      {"Lc", 7.0},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       {{{"Attributes", {1}},
         {"Permittivity", {{"HighFrequencyLimit", 2.0}, {"Terms", terms}}}}}}};
  setup_json["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 50.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}},
                              {"Absorbing", {{"Attributes", {101}}, {"Order", 1}}}};
  setup_json["Solver"] = {{"Order", 2},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-6},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 18.0},
                            {"MaxFreq", 22.0},
                            {"FreqStep", 0.1}}}};
  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
  SpaceOperator space_op(iodata, mesh_io);
  const auto &farfield_op = space_op.GetSurfaceTerms().GetFarfieldOp();
  REQUIRE(farfield_op.NumDispersiveMaterials() == 1);

  RomOperatorTest prom_op(iodata, space_op, 10);
  prom_op.AddLumpedPortModesForSynthesis();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(Mpi::World(), v.Real());
  prom_op.UpdatePROM(v, "random");
  const long n = prom_op.GetReducedDimension();
  REQUIRE(n == 2);
  const auto [w_min, w_max] = prom_op.GetSweepBand();
  REQUIRE(w_max > w_min);

  // i·VᵀMV by the imaginary-slot convention, replaced by i·R diag(1, 1e-7) Rᵀ.
  const double c = std::cos(0.3), s = std::sin(0.3);
  Eigen::Matrix2d R, D = Eigen::Vector2d(1.0, 1.0e-7).asDiagonal();
  R << c, -s, s, c;
  const Eigen::MatrixXcd M_r =
      std::complex<double>(0.0, 1.0) * (R * D * R.transpose()).cast<std::complex<double>>();
  auto f = [&farfield_op](std::complex<double> w) -> std::complex<double>
  { return w * farfield_op.GetDispersiveScale(0, w); };
  Eigen::MatrixXcd Kc = Eigen::MatrixXcd::Zero(n, n), Cc = Kc, Mc = Kc;
  std::vector<RomOperatorTest::WavePortAuxBlock> blocks;
  prom_op.AddScalarDispersionSynthesis(
      "farfield", M_r, prom_op.GetFarfieldDispersiveScalar(0, false), Kc, Cc, Mc, blocks);
  REQUIRE(!blocks.empty());  // The fit has poles
  const auto fit = prom_op.FitScalarDispersion("check", M_r, f, true);
  REQUIRE(fit.regime == RomOperatorTest::WavePortRegime::Augmented);
  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kc, Cc, Mc, blocks, aux_labels);
  const long n_aux = aug.Kr.rows() - n;
  REQUIRE(n_aux > 0);
  double max_rel = 0.0;
  for (int k = 0; k <= 80; k++)
  {
    const double w = w_min + (w_max - w_min) * k / 80.0;
    const Eigen::MatrixXcd A =
        aug.Kr + std::complex<double>(0.0, w) * aug.Cr - w * w * aug.Mr;
    const Eigen::MatrixXcd schur =
        A.topLeftCorner(n, n) -
        A.topRightCorner(n, n_aux) * A.bottomRightCorner(n_aux, n_aux)
                                         .fullPivLu()
                                         .solve(A.bottomLeftCorner(n_aux, n));
    const Eigen::MatrixXcd want = f(w) * M_r;
    max_rel = std::max(max_rel,
                       (schur - want).cwiseAbs().maxCoeff() / want.cwiseAbs().maxCoeff());
  }
  CAPTURE(fit.rel_err_augmented, max_rel);
  // The weak direction is not realized, an error of 1e-7 of the term, as the rank cutoff
  // allows.
  CHECK(max_rel <= 2.0 * fit.rel_err_augmented + 2.0e-7);
}

// Over a wide band (0.1 to 20 GHz), the fitted frequency-dependent terms of circuit
// synthesis vary by orders of magnitude, a Djordjevic-Sarkar volume term g = s²χ(s) as ω²
// and the absorbing boundary terms as ω and 1/ω. A fit residual normalized by the largest
// value over the band does not bound the error at the frequencies where a term is small.
// The reported fit residual must bound the realized error at every frequency, relative to
// the term there, and a fit which misses the synthesis tolerance must report it.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "RomOperator synthesis fit residuals hold pointwise over a wide band",
                 "[romoperator][Serial][Parallel]")
{
  const json terms = {{{"Type", "DjordjevicSarkar"},
                       {"Strength", 0.2},
                       {"LowerFrequency", 1.0},
                       {"UpperFrequency", 100.0}}};
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
      {"L0", 1.0e-6},
      {"Lc", 7.0},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       {{{"Attributes", {1}},
         {"Permittivity", {{"HighFrequencyLimit", 2.0}, {"Terms", terms}}}}}}};
  setup_json["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 50.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}},
                              {"Absorbing", {{"Attributes", {101}}, {"Order", 2}}}};
  setup_json["Solver"] = {{"Order", 2},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-4},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 0.1},
                            {"MaxFreq", 20.0},
                            {"FreqStep", 0.1}}}};
  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
  SpaceOperator space_op(iodata, mesh_io);
  const auto &mat_op = space_op.GetMaterialOp();
  const auto &farfield_op = space_op.GetSurfaceTerms().GetFarfieldOp();
  REQUIRE(farfield_op.NumDispersiveMaterials() == 1);

  RomOperatorTest prom_op(iodata, space_op, 10);
  prom_op.AddLumpedPortModesForSynthesis();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(Mpi::World(), v.Real());
  prom_op.UpdatePROM(v, "random");
  const long n = prom_op.GetReducedDimension();
  REQUIRE(n == 2);
  const auto [w_min, w_max] = prom_op.GetSweepBand();
  REQUIRE(w_max > 100.0 * w_min);
  const double tol = prom_op.GetSynthesisTol();

  constexpr std::complex<double> i(0.0, 1.0);
  enum class Term
  {
    VOLUME,
    DAMPING,
    CURL_CURL
  };
  for (const auto term : {Term::VOLUME, Term::DAMPING, Term::CURL_CURL})
  {
    CAPTURE(static_cast<int>(term));
    const Eigen::MatrixXcd &M_r =
        (term == Term::VOLUME)
            ? prom_op.GetFrequencyDependentPermittivityMassR(0)
            : prom_op.GetFarfieldDispersiveMatrixR(0, term == Term::CURL_CURL);
    REQUIRE(M_r.rows() == n);
    // The scalar f(ω) of the term i·f(ω)·M_r, as fit by the synthesis.
    auto f = [&](std::complex<double> w) -> std::complex<double>
    {
      if (term == Term::VOLUME)
      {
        return mat_op.EvaluateDjordjevicSarkarPermittivityA2(0, i * w) / i;
      }
      const auto r = farfield_op.GetDispersiveScale(0, w);
      return (term == Term::CURL_CURL) ? 0.5 / (w * r) : w * r;
    };
    Eigen::MatrixXcd Kc = Eigen::MatrixXcd::Zero(n, n), Cc = Kc, Mc = Kc;
    std::vector<RomOperatorTest::WavePortAuxBlock> blocks;
    if (term == Term::VOLUME)
    {
      prom_op.AddFrequencyDependentPermittivitySynthesis("permittivity_0", M_r, 0, Kc, Cc,
                                                         Mc, blocks);
    }
    else
    {
      prom_op.AddScalarDispersionSynthesis(
          "farfield", M_r, prom_op.GetFarfieldDispersiveScalar(0, term == Term::CURL_CURL),
          Kc, Cc, Mc, blocks);
    }
    const auto fit = prom_op.FitScalarDispersion("check", M_r, f, true);
    const double reported = (fit.regime == RomOperatorTest::WavePortRegime::Augmented)
                                ? fit.rel_err_augmented
                                : fit.rel_err_polynomial;
    std::vector<std::string> aux_labels;
    auto aug = RomOperatorTest::BuildAugmentedPencil(Kc, Cc, Mc, blocks, aux_labels);
    const long n_aux = aug.Kr.rows() - n;
    double max_rel = 0.0;
    for (int k = 0; k <= 40; k++)
    {
      // Logarithmically spaced over the band, also between the fit nodes.
      const double w = w_min * std::pow(w_max / w_min, k / 40.0);
      const Eigen::MatrixXcd A =
          aug.Kr + std::complex<double>(0.0, w) * aug.Cr - w * w * aug.Mr;
      Eigen::MatrixXcd schur = A.topLeftCorner(n, n);
      if (n_aux > 0)
      {
        schur -= A.topRightCorner(n, n_aux) * A.bottomRightCorner(n_aux, n_aux)
                                                  .fullPivLu()
                                                  .solve(A.bottomLeftCorner(n_aux, n));
      }
      Eigen::MatrixXcd want = f(w) * M_r;
      if (term == Term::VOLUME)
      {
        want = mat_op.EvaluateFrequencyDependentPermittivityA2(0, i * w) * (M_r / i);
      }
      max_rel = std::max(max_rel,
                         (schur - want).cwiseAbs().maxCoeff() / want.cwiseAbs().maxCoeff());
    }
    CAPTURE(tol, reported, max_rel);
    // The realized error is bounded by the reported residual (up to the sampling of the
    // validation grid), and the tolerance is met unless the residual reports otherwise.
    CHECK(max_rel <= 2.0 * reported + 1.0e-12);
    CHECK((max_rel <= 2.0 * tol || reported > tol));
  }
}

// Circuit synthesis fits the response of a wave port or absorbing boundary next to a
// frequency-dependent material, which is not rational near a resonance of the material.
// Before the sweep, one warning names each such boundary next to a material with a
// resonance inside the band of the sweep. None is printed for a resonance outside the band,
// or without circuit synthesis.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "RomOperator warns of a resonance inside the band of circuit synthesis",
                 "[romoperator][Serial][Parallel]")
{
  const std::string warning = " is next to a frequency-dependent material with a resonance "
                              "inside the band of the sweep, at ";
  auto Lorentz = [](double f0, double gamma)
  {
    return json{{"HighFrequencyLimit", 2.0},
                {"Terms",
                 {{{"Type", "Lorentz"},
                   {"DeltaPermittivity", 0.5},
                   {"ResonanceFrequency", f0},
                   {"DampingFrequency", gamma}}}}};
  };
  // The output of the construction of the reduced-order model, on the root.
  auto Construct = [](json &config)
  {
    IoData iodata(config, false);
    auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
    SpaceOperator space_op(iodata, mesh_io);
    return test::CaptureStdout([&] { RomOperatorTest prom_op(iodata, space_op, 10); });
  };
  auto Count = [&warning](const std::string &out)
  {
    std::size_t count = 0;
    for (auto pos = out.find(warning); pos != std::string::npos;
         pos = out.find(warning, pos + 1))
    {
      count++;
    }
    return count;
  };

  SECTION("Absorbing boundary")
  {
    // The pole of the Lorentz term at 20 GHz with a damping of 2 GHz is at √(20² − 1) GHz.
    for (const auto &[max_freq, synthesis, warned] :
         std::vector<std::tuple<double, bool, bool>>{
             {32.0, true, true}, {15.0, true, false}, {32.0, false, false}})
    {
      CAPTURE(max_freq, synthesis);
      json config;
      config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
      config["Model"] = {{"Mesh", fs::path(PALACE_TEST_DATA_DIR) /
                                      "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
                         {"L0", 1.0e-6},
                         {"Lc", 7.0},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
      config["Domains"] = {
          {"Materials", {{{"Attributes", {1}}, {"Permittivity", Lorentz(20.0, 2.0)}}}}};
      config["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 50.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}},
                              {"Absorbing", {{"Attributes", {101}}, {"Order", 2}}}};
      config["Solver"] = {{"Order", 1},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-6},
                            {"AdaptiveCircuitSynthesis", synthesis},
                            {"MinFreq", 2.0},
                            {"MaxFreq", max_freq},
                            {"FreqStep", 1.0}}}};
      const auto out = Construct(config);
      if (Mpi::Root(Mpi::World()))
      {
        CAPTURE(out);
        CHECK(Count(out) == (warned ? 1 : 0));
        CHECK((out.find("The absorbing boundary (attributes 101)" + warning +
                        "19.975 GHz (material 1)!") != std::string::npos) == warned);
      }
    }
  }

  SECTION("Wave ports")
  {
    // The coaxial line of coaxial_wave_losstan, filled with a Lorentz medium whose pole
    // is at √(10² − 0.25²) GHz, inside the band from 5 to 15 GHz, or at 20 GHz, above it.
    const auto dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / "coaxial";
    for (const auto &[f0, warned] :
         std::vector<std::pair<double, bool>>{{10.0, true}, {20.0, false}})
    {
      CAPTURE(f0);
      auto config = IoData::ParseAndValidate((dir / "wave_losstan.json").string().c_str());
      config["Model"]["Mesh"] = (dir / config["Model"]["Mesh"].get<std::string>()).string();
      config["Problem"]["Output"] = temp_dir;
      config["Problem"]["Verbose"] = 0;
      config["Domains"]["Materials"][0].erase("LossTan");
      config["Domains"]["Materials"][0]["Permittivity"] = Lorentz(f0, 0.5);
      config["Solver"]["Order"] = 1;
      config["Solver"]["Driven"]["AdaptiveTol"] = 1.0e-3;
      config["Solver"]["Driven"]["AdaptiveCircuitSynthesis"] = true;
      const auto out = Construct(config);
      if (Mpi::Root(Mpi::World()))
      {
        CAPTURE(out);
        CHECK(Count(out) == (warned ? 2 : 0));
        for (const std::string port : {"1", "2"})
        {
          CHECK((out.find("Wave port " + port + warning + "9.99687 GHz (material 1)!") !=
                 std::string::npos) == warned);
        }
      }
    }
  }
}

// A root of the synthesized circuit with Im{f} < 0 grows in time, which its Q, from
// |Im{f}|, does not show. Writing the circuit warns once, with the number of such roots
// among the estimates and the one which grows fastest; a stable circuit does not warn.
TEST_CASE_METHOD(palace::test::PerRankTempDir,
                 "RomOperator warns of growing roots of the synthesized circuit",
                 "[romoperator][Serial][Parallel]")
{
  const std::string warning = " synthesized-system eigenvalue estimates have Im{f} < 0, "
                              "down to f = ";
  using Estimate = RomOperatorTest::EigenvalueEstimate;
  auto Fastest = [](const std::vector<Estimate> &eigs)
  { return *std::ranges::min_element(eigs, {}, &Estimate::freq_im_GHz); };
  auto Format = [](const Estimate &e)
  { return fmt::format("{:+.6e} {:+.6e}i GHz!", e.freq_re_GHz, e.freq_im_GHz); };

  SECTION("Estimates")
  {
    // Uncoupled parallel resonators: (L⁻¹ + iωR⁻¹ − ω²C) v = 0 with L⁻¹ = ω₀²C has the
    // roots ω = ±√(ω₀² − (R⁻¹/2C)²) + i R⁻¹/2C, which grow for a negative conductance.
    const double C = 1.0e-12;
    const std::vector<double> f0 = {5.0, 8.0, 11.0}, im_f = {0.1, -0.2, -0.5};  // GHz
    for (const bool stable : {false, true})
    {
      CAPTURE(stable);
      const long n = static_cast<long>(f0.size());
      Eigen::MatrixXcd L_inv = Eigen::MatrixXcd::Zero(n, n), R_inv = L_inv, Cm = L_inv;
      for (long k = 0; k < n; k++)
      {
        const double w0 = 2.0 * std::numbers::pi * 1.0e9 * f0[k];
        const double g = 4.0 * std::numbers::pi * 1.0e9 * C * std::abs(im_f[k]);
        Cm(k, k) = C;
        L_inv(k, k) = w0 * w0 * C;
        R_inv(k, k) = (stable || im_f[k] > 0.0) ? g : -g;
      }
      const auto eigs =
          RomOperatorTest::ComputeEigenvalueEstimates(L_inv, &R_inv, Cm, 2.0, 14.0);
      REQUIRE(eigs.size() == f0.size());
      const auto fastest = Fastest(eigs);
      CHECK_THAT(fastest.freq_im_GHz, WithinRel(stable ? 0.1 : -0.5, 1.0e-10));
      const auto out = test::CaptureStdout(
          [&] { RomOperatorTest::WarnGrowingEigenvalueEstimates(Mpi::World(), eigs); });
      if (Mpi::Root(Mpi::World()))
      {
        CAPTURE(out);
        CHECK((out.find("2 of 3" + warning + Format(fastest)) != std::string::npos) ==
              !stable);
        CHECK(out.empty() == stable);
      }
    }
  }

  SECTION("Circuit")
  {
    // The port mode of a lumped port on a 1 cm cube, in a medium with a conductivity of
    // ±0.22 S/m (a zero pole r/s of the permittivity, r = σ/ε₀). The negative conductivity
    // of a medium with gain outweighs the loss of the port, and the root grows.
    for (const double r : {-2.5e10, 2.5e10})
    {
      CAPTURE(r);
      json config;
      config["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
      config["Model"] = {{"Mesh", fs::path(PALACE_TEST_DATA_DIR) /
                                      "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"},
                         {"L0", 1.0e-2},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
      config["Domains"] = {
          {"Materials",
           {{{"Attributes", {1}},
             {"Permittivity",
              {{"HighFrequencyLimit", 1.0},
               {"Terms", {{{"Type", "PoleResidue"}, {"Pole", 0.0}, {"Residue", r}}}}}}}}}};
      config["Boundaries"] = {{"LumpedPort",
                               {{{"Index", 1},
                                 {"R", 2000.0},
                                 {"Excitation", 1},
                                 {"Attributes", {100}},
                                 {"Direction", "+X"}}}}};
      config["Solver"] = {{"Order", 1},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 1.0},
                            {"MaxFreq", 1000.0},
                            {"FreqStep", 10.0}}}};
      IoData iodata(config, false);
      auto mesh_io = LoadScaleParMesh2(iodata, Mpi::World());
      SpaceOperator space_op(iodata, mesh_io);
      RomOperatorTest prom_op(iodata, space_op, 10);
      prom_op.AddLumpedPortModesForSynthesis();

      // The roots which writing the circuit estimates, on the band of the sweep.
      const auto matrices = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
      const auto eigs = RomOperatorTest::ComputeEigenvalueEstimates(
          *matrices.L_inv, matrices.R_inv.get(), *matrices.C, 1.0, 1000.0);
      REQUIRE(eigs.size() == 1);
      const auto fastest = Fastest(eigs);
      REQUIRE((fastest.freq_im_GHz < 0.0) == (r < 0.0));
      const auto out =
          test::CaptureStdout([&] { prom_op.PrintPROMMatrices(iodata.units, temp_dir); });
      if (Mpi::Root(Mpi::World()))
      {
        CAPTURE(out);
        CHECK((out.find("1 of 1" + warning + Format(fastest)) != std::string::npos) ==
              (r < 0.0));
        CHECK((out.find(warning) != std::string::npos) == (r < 0.0));
      }
    }
  }
}

// Basic checks of ROM construction for synthesis. Checks hybrid domain-boundary
// inner-product weight and port overlap. Works with a simple 1x1x1 Cube. This is a serial
// test as hex mesh only has a single element.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "RomOperator-Synthesis-Port-Cube111",
                 "[romoperator][Serial]")
{
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 different test configuration.
  size_t order = GENERATE(1UL, 2UL, 3UL);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_1_1_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"));

  double L0 = 1.0e-6;

  // Pick a number != 1. This should not matter as it should be be correctly cancelled out
  // in unit conversions, but we want to test this here.
  double Lc = 7.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0},
                                               {"LossTan", 0.0}})})}};

  double port_ref_R = 50.0;

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      // {"PEC", json::object({{"Attributes", json::array({})}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", port_ref_R},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  // Test hybrid weight operator works as expected.
  const auto &weight_op = prom_op.GetWeightOp();
  REQUIRE(weight_op.has_value());

  const auto &W_bulk = weight_op->W_inner_product_weight_bulk;
  const auto &W_port = weight_op->W_inner_product_weight_port;

  std::size_t nr_tdof_expected = 0;
  std::size_t nr_tdof_port_expected = 0;
  if (mesh_is_hex)
  {
    // Reference vales for nr elements for vector Nédélec elements of the first kind.
    auto nr_tdof_ref_hex = std::vector<std::size_t>{0, 12, 54, 144}.at(order);
    auto nr_tdof_ref_squ = std::vector<std::size_t>{0, 4, 12, 24}.at(order);

    nr_tdof_expected = nr_tdof_ref_hex;
    nr_tdof_port_expected = nr_tdof_ref_squ;
  }
  else
  {
    auto nr_tdof_ref_tet = std::vector<std::size_t>{0, 6, 20, 45}.at(order);
    auto nr_tdof_ref_tri = std::vector<std::size_t>{0, 3, 8, 15}.at(order);

    // Counting tdofs. 3D: There are 6 tets in a cube; remove edges. 2D: There are 2 * 6 =
    // 12 open faces, so 12 double-counted; subtract 6 faces without edges. 1D: add total of
    // 19 edges
    nr_tdof_expected = 6 * (nr_tdof_ref_tet - 6 * order) -
                       12 / 2 * (nr_tdof_ref_tri - 3 * order) + 19 * order;
    nr_tdof_port_expected = 2 * nr_tdof_ref_tri - order;  // 2 triangles, 1 shared edge
  }

  CHECK(W_bulk->NumRows() == nr_tdof_expected);
  CHECK(W_bulk->NumCols() == nr_tdof_expected);

  auto W_port_eigen = toEigenMatrix(*W_port, nr_tdof_expected);
  auto W_bulk_eigen = toEigenMatrix(*W_bulk, nr_tdof_expected);
  auto weight_op_eigen = toEigenMatrix(*weight_op, nr_tdof_expected);

  // Check rows/cols where port matrix is non-zero and ensure that corresponding domain
  // matrix is zero.
  std::size_t nr_rows_count = 0;
  for (int i = 0; i < W_port_eigen.rows(); i++)
  {
    if (W_port_eigen.row(i).norm() > 0.0)
    {
      CHECK(W_bulk_eigen.row(i).norm() == 0.0);
      nr_rows_count++;
    }
    else
    {
      // Check there is some overlap in other rows.
      CHECK(W_bulk_eigen.row(i).norm() > 0.0);
    }
  }
  CHECK(nr_rows_count == nr_tdof_port_expected);

  std::size_t nr_cols_count = 0;
  for (int j = 0; j < W_port_eigen.cols(); j++)
  {

    if (W_port_eigen.col(j).norm() > 0.0)
    {
      CHECK(W_bulk_eigen.col(j).norm() == 0.0);
      nr_cols_count++;
    }
    else
    {
      // Check there is some overlap in other cols.
      CHECK(W_bulk_eigen.col(j).norm() > 0.0);
    }
  }
  CHECK(nr_cols_count == nr_tdof_port_expected);

  // Check element-wise sum is ok.
  for (int i = 0; i < weight_op_eigen.rows(); i++)
  {
    for (int j = 0; j < weight_op_eigen.cols(); j++)
    {
      CHECK_THAT(weight_op_eigen(i, j),
                 WithinRel(W_bulk_eigen(i, j) + W_port_eigen(i, j)) ||
                     WithinAbs(0.0, 1e-18));
    }
  }

  // Now test against port vector. The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e
  // L_e / W_e but with Z_R = 1.0. In this case this = 1.0. See normalization tests of
  // LumpedPort_BasicTests_1ElementPort_Cube321 in test-lumpedportintegration.cpp.
  ComplexVector port_primary_et;
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et);
  Vector port_primary_ht_cn_tmp = port_primary_et.Real();

  W_port->Mult(port_primary_et.Real(), port_primary_ht_cn_tmp);
  auto overlap_port = port_primary_et.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_port, WithinRel(1.0));

  W_bulk->Mult(port_primary_et.Real(), port_primary_ht_cn_tmp);
  auto overlap_bulk = port_primary_et.Real() * port_primary_ht_cn_tmp;
  CHECK_THAT(overlap_bulk, WithinAbs(0.0, 1e-15));

  auto v_workspace = port_primary_et.Real();

  if (Mpi::Size(world_comm) == 1)
  {
    // Rank local overlap.
    auto overlap_combined_local =
        (*weight_op)
            .InnerProduct(port_primary_et.Real(), port_primary_et.Real(), v_workspace);
    CHECK_THAT(overlap_combined_local, WithinRel(1.0));
  }

  // Global overlap.
  auto overlap_combined = (*weight_op)
                              .InnerProduct(world_comm, port_primary_et.Real(),
                                            port_primary_et.Real(), v_workspace);
  CHECK_THAT(overlap_combined, WithinRel(1.0));

  // Test actually adding port primary vectors to PROM.
  prom_op.AddLumpedPortModesForSynthesis();
  CHECK(prom_op.GetReducedDimension() == 1);
  const auto _norm = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = _norm.L_inv;
  const auto &resistance_R_inv = _norm.R_inv;
  const auto &capacitance_C = _norm.C;
  const auto orth_R = prom_op.GetOrthR();

  CHECK(((inductance_L_inv->rows() == 1) && (inductance_L_inv->cols() == 1)));
  CHECK(((resistance_R_inv->rows() == 1) && (resistance_R_inv->cols() == 1)));
  CHECK(((capacitance_C->rows() == 1) && (capacitance_C->cols() == 1)));
  CHECK(((orth_R.rows() == 1) && (orth_R.cols() == 1)));

  // This should be the same as norm of primary port vector ht_cn with Z_R = 1.0.
  // For single square port this is 1.0
  CHECK_THAT(orth_R(0, 0), WithinRel(1.0));
  CHECK_THAT((*resistance_R_inv)(0, 0).real(), WithinRel(1.0 / port_ref_R));
  // L, C have some component on port 1 from the bulk. So only check they exist.
  CHECK((*inductance_L_inv)(0, 0).real() != 0.0);
  CHECK((*capacitance_C)(0, 0).real() != 0.0);

  CHECK_THAT((*inductance_L_inv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*resistance_R_inv)(0, 0).imag(), WithinAbs(0.0, 1e-14));
  CHECK_THAT((*capacitance_C)(0, 0).imag(), WithinAbs(0.0, 1e-14));
}

// Basic ROM port construction on a more complex mesh including a composite port and and L,C
// port. Works on a 3x2x1 cubic mesh.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-Synthesis-Port-Cube321",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 different test configuration.
  size_t order = GENERATE(1UL, 2UL, 3UL);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));

  const double L0 = 1.0e-6;
  const double Lc = 7.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  // Reasonable reference values in correct frequency range.
  const double port1_ref_R = 75.0;
  const double port1_ref_C = 5.5e-5;
  const double port1_ref_L = 1.486e-20;

  const double port2_ref_R = 50.0;

  // Put in a single port with single attribute.
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({
           json::object(
               {{"Index", 1},
                {"R", port1_ref_R},
                {"C", port1_ref_C},
                {"L", port1_ref_L},
                {"Excitation", uint(0)},
                {"Elements",
                 json::array(
                     {json::object({{"Attributes", json::array({1})},  // dx,dy=1,1
                                    {"Direction", "+Y"}}),
                      json::object({{"Attributes", json::array({9, 11})},  // dx,dy=1,2
                                    {"Direction", "+Y"}}),
                      json::object({{"Attributes", json::array({2, 6, 10})},  // dx,dy = 3,1
                                    {"Direction", "+X"}})})}}),
           json::object({{"Index", 2},
                         {"R", port2_ref_R},
                         {"Excitation", uint(1)},
                         {"Elements", json::array({json::object(
                                          {{"Attributes", json::array({14})},  // dx,dy=1,1
                                           {"Direction", "+Z"}})})}})          //,
       })}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  // Test hybrid weight operator works as expected.
  const auto &weight_op = prom_op.GetWeightOp();
  REQUIRE(weight_op.has_value());

  const auto &W_bulk = weight_op->W_inner_product_weight_bulk;
  const auto &W_port = weight_op->W_inner_product_weight_port;

  // If serial — concretize as an Eigen matrix on a thread and manually check zeros of row
  // and columns.
  if (Mpi::Size(world_comm) == 1)
  {

    auto W_port_eigen = toEigenMatrix(*W_port, W_port->NumCols());
    auto W_bulk_eigen = toEigenMatrix(*W_bulk, W_port->NumCols());

    // Check rows/cols where port matrix is non-zero and ensure that corresponding domain
    // matrix is zero.
    for (int i = 0; i < W_port_eigen.rows(); i++)
    {
      if (W_port_eigen.row(i).norm() > 0.0)
      {
        CHECK(W_bulk_eigen.row(i).norm() == 0.0);
      }
      else
      {
        // Check there is some overlap in other rows.
        CHECK(W_bulk_eigen.row(i).norm() > 0.0);
      }
    }

    for (int j = 0; j < W_port_eigen.cols(); j++)
    {

      if (W_port_eigen.col(j).norm() > 0.0)
      {
        CHECK(W_bulk_eigen.col(j).norm() == 0.0);
      }
      else
      {
        // Check there is some overlap in other cols.
        CHECK(W_bulk_eigen.col(j).norm() > 0.0);
      }
    }
  }

  // Now test against port vector. The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e
  // L_e / W_e but with Z_R = 1.0. In this case this = 1.0 for Port 2,3 and != 1.0 for Port
  // 1. See normalization tests of LumpedPort_BasicTests_1ElementPort_Cube321 in
  //    test-lumpedportintegration.cpp.

  for (const auto &[port_i, port_norm] :
       std::vector<std::pair<int, double>>{{1, 1.0}, {2, 1.}})
  {
    ComplexVector port_primary_et;
    space_op.GetLumpedPortExcitationVectorPrimaryEt(port_i, port_primary_et);
    Vector port_primary_et_tmp = port_primary_et.Real();

    W_port->Mult(port_primary_et.Real(), port_primary_et_tmp);
    auto overlap_port =
        linalg::Dot(world_comm, port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_port, WithinRel(port_norm));

    W_bulk->Mult(port_primary_et.Real(), port_primary_et_tmp);
    auto overlap_bulk =
        linalg::Dot(world_comm, port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_bulk, WithinAbs(0.0, 1e-15));

    auto overlap_combined = (*weight_op)
                                .InnerProduct(world_comm, port_primary_et.Real(),
                                              port_primary_et.Real(), port_primary_et_tmp);
    CHECK_THAT(overlap_combined, WithinRel(port_norm));
  }

  // // Test actually adding port primary vectors to PROM.
  prom_op.AddLumpedPortModesForSynthesis();
  auto rom_dim = prom_op.GetReducedDimension();
  REQUIRE(rom_dim == 2);
  const auto _norm = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = _norm.L_inv;
  const auto &resistance_R_inv = _norm.R_inv;
  const auto &capacitance_C = _norm.C;
  const auto orth_R = prom_op.GetOrthR();

  CHECK(((inductance_L_inv->rows() == rom_dim) && (inductance_L_inv->cols() == rom_dim)));
  CHECK(((resistance_R_inv->rows() == rom_dim) && (resistance_R_inv->cols() == rom_dim)));
  CHECK(((capacitance_C->rows() == rom_dim) && (capacitance_C->cols() == rom_dim)));
  CHECK(((orth_R.rows() == rom_dim) && (orth_R.cols() == rom_dim)));

  // Per-port load matrices are exported separately from the total loaded matrices so
  // downstream circuit assembly can remove internal terminations and keep external ones.
  // For this two-port lumped case they should have the same shape as the total matrices
  // and contain only the diagonal R/L/C loading terms for the corresponding port row.
  REQUIRE(_norm.port_loads.size() == 2);
  auto find_port_load =
      [&](std::string_view label) -> const RomOperatorTest::NormalizedMatrices::PortLoad *
  {
    for (const auto &load : _norm.port_loads)
    {
      if (load.label == label)
      {
        return &load;
      }
    }
    return nullptr;
  };
  const auto *port1_load = find_port_load("port_1_re");
  const auto *port2_load = find_port_load("port_2_re");
  REQUIRE(port1_load != nullptr);
  REQUIRE(port2_load != nullptr);

  auto check_load_shape = [rom_dim](const auto &load)
  {
    REQUIRE(load.L_inv);
    REQUIRE(load.C);
    CHECK(load.L_inv->rows() == rom_dim);
    CHECK(load.L_inv->cols() == rom_dim);
    CHECK(load.C->rows() == rom_dim);
    CHECK(load.C->cols() == rom_dim);
    if (load.R_inv)
    {
      CHECK(load.R_inv->rows() == rom_dim);
      CHECK(load.R_inv->cols() == rom_dim);
    }
  };
  check_load_shape(*port1_load);
  check_load_shape(*port2_load);

  auto check_only_diagonal =
      [rom_dim](const Eigen::MatrixXcd &mat, long row, std::complex<double> expected)
  {
    for (long i = 0; i < rom_dim; i++)
    {
      for (long j = 0; j < rom_dim; j++)
      {
        if (i == row && j == row)
        {
          CHECK_THAT(mat(i, j).real(), WithinRel(expected.real(), 1e-12));
          CHECK_THAT(mat(i, j).imag(), WithinAbs(expected.imag(), 1e-15));
        }
        else
        {
          CHECK_THAT(std::abs(mat(i, j)), WithinAbs(0.0, 1e-15));
        }
      }
    }
  };
  check_only_diagonal(*port1_load->L_inv, 0, 1.0 / port1_ref_L);
  check_only_diagonal(*port1_load->C, 0, port1_ref_C);
  REQUIRE(port1_load->R_inv);
  check_only_diagonal(*port1_load->R_inv, 0, 1.0 / port1_ref_R);

  CHECK(port2_load->L_inv->cwiseAbs().maxCoeff() == 0.0);
  CHECK(port2_load->C->cwiseAbs().maxCoeff() == 0.0);
  REQUIRE(port2_load->R_inv);
  check_only_diagonal(*port2_load->R_inv, 1, 1.0 / port2_ref_R);

  Eigen::MatrixXd inductance_L_inv_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  inductance_L_inv_ref(0, 0) = 1.0 / port1_ref_L;

  Eigen::MatrixXd resistance_R_inv_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  resistance_R_inv_ref(0, 0) = 1.0 / port1_ref_R;
  resistance_R_inv_ref(1, 1) = 1.0 / port2_ref_R;
  // resistance_R_inv_ref(2, 2) = 1.0 / port3_ref_R;

  Eigen::MatrixXd capacitance_C_ref = Eigen::MatrixXd::Zero(rom_dim, rom_dim);
  capacitance_C_ref(0, 0) = port1_ref_C;

  auto check_mat_zero_pattern =
      [dim = prom_op.GetReducedDimension()](auto &mat, bool skip_real = false)
  {
    for (long i = 0; i < dim; i++)
    {
      for (long j = 0; j < dim; j++)
      {
        if (i != j && !skip_real)
        {
          CHECK_THAT(std::real(mat(i, j)), WithinAbs(0.0, 1e-15));
        }
        CHECK_THAT(std::imag(mat(i, j)), WithinAbs(0.0, 1e-15));
      }
    }
  };

  check_mat_zero_pattern(orth_R);
  check_mat_zero_pattern(*resistance_R_inv);
  check_mat_zero_pattern(*capacitance_C);
  check_mat_zero_pattern(*inductance_L_inv,
                         true);  // skip: could extras overlap crowded mesh

  // Port vectors are chosen so that after orthogonalization with the hybrid mass matrix
  // their norm is 1.0.
  CHECK_THAT(std::real(orth_R(0, 0)), WithinRel(1.0));
  CHECK_THAT(std::real(orth_R(1, 1)), WithinRel(1.0));

  CHECK_THAT(std::real((*resistance_R_inv)(0, 0)), WithinRel(resistance_R_inv_ref(0, 0)));
  CHECK_THAT(std::real((*resistance_R_inv)(1, 1)), WithinRel(resistance_R_inv_ref(1, 1)));

  // Don't check 11 of due to presence of bulk values.
  // For 00, bulk values are small as compared to port values, put tolerance of ~1e-6.
  CHECK_THAT(std::real((*inductance_L_inv)(0, 0)),
             WithinRel(inductance_L_inv_ref(0, 0), 1e-6));
  CHECK_THAT(std::real((*capacitance_C)(0, 0)), WithinRel(capacitance_C_ref(0, 0), 1e-6));

  // Add some random vectors to mimic solution of system.
  ComplexVector vector_temp;
  vector_temp.Real() = prom_op.GetVectors().at(0);  // Steal-setup
  vector_temp.Imag() = prom_op.GetVectors().at(0);  // Steal-setup

  std::size_t nr_random_vec = 10;

  std::seed_seq seed_gen{2 * Mpi::Rank(world_comm)};
  std::vector<std::uint32_t> seeds(2 * nr_random_vec);
  seed_gen.generate(seeds.begin(), seeds.end());

  // Return V c for coefficients c zero-padded to the basis dimension.
  auto expand = [&](const Eigen::VectorXcd &c)
  {
    Eigen::VectorXcd c_pad = Eigen::VectorXcd::Zero(prom_op.GetReducedDimension());
    c_pad.head(c.size()) = c;
    ComplexVector w(vector_temp.Size());
    w.UseDevice(true);
    linalg::LinearCombination(prom_op.GetVectors(),
                              {c_pad.data(), static_cast<std::size_t>(c_pad.size())}, w);
    return w;
  };

  // The projected matrices are those of the operators on the basis, Vᵀ A V.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  REQUIRE(C);
  auto check_projection = [&](const ComplexOperator &A, const Eigen::MatrixXcd &Ar)
  {
    const auto &V = prom_op.GetVectors();
    const auto n = static_cast<Eigen::Index>(V.size());
    REQUIRE(Ar.rows() == n);
    REQUIRE(Ar.cols() == n);
    Eigen::MatrixXcd Ar_ref(n, n);
    ComplexVector v(vector_temp.Size()), w(vector_temp.Size());
    v.UseDevice(true);
    w.UseDevice(true);
    for (Eigen::Index j = 0; j < n; j++)
    {
      v.Real() = V[j];
      v.Imag() = 0.0;
      A.Mult(v, w);
      for (Eigen::Index i = 0; i < n; i++)
      {
        Ar_ref(i, j) = {V[i] * w.Real(), V[i] * w.Imag()};
      }
    }
    Mpi::GlobalSum(static_cast<int>(Ar_ref.size()), Ar_ref.data(), world_comm);
    const double error = (Ar - Ar_ref).cwiseAbs().maxCoeff();
    const double scale = Ar_ref.cwiseAbs().maxCoeff();
    CAPTURE(n, error, scale);
    REQUIRE(scale > 0.0);
    CHECK(error <= 1.0e-12 * scale);
  };

  Eigen::VectorXcd c_prev;
  auto check_update = [&](const ComplexVector &u, const Eigen::VectorXcd &c)
  {
    // The returned coefficients reproduce the vector in the updated basis.
    const double norm_u = linalg::Norml2(world_comm, u);
    ComplexVector w = expand(c);
    w -= u;
    CHECK(linalg::Norml2(world_comm, w) <= 1.0e-12 * norm_u);

    // The norm in coefficient space from the Gram factor, VᵀV = L Lᵀ, matches the norm of
    // the vector, for the new vector and for its difference with the previous one, which
    // involves the whole factor.
    CHECK_THAT(prom_op.BasisNorm(c), WithinRel(norm_u, 1.0e-10));
    if (c_prev.size() > 0)
    {
      Eigen::VectorXcd dc = c;
      dc.head(c_prev.size()) -= c_prev;
      ComplexVector du = expand(c_prev);
      du -= u;
      CHECK_THAT(prom_op.BasisNorm(dc), WithinRel(linalg::Norml2(world_comm, du), 1.0e-10));
    }
    c_prev = c;
  };

  for (std::size_t i = 0; i < nr_random_vec; i++)
  {
    vector_temp.Real().Randomize(seeds.at(2 * i));
    vector_temp.Imag().Randomize(seeds.at(2 * i + 1));
    const auto c = prom_op.UpdatePROM(vector_temp, fmt::format("vec_{}", i));
    REQUIRE(c.size() == static_cast<Eigen::Index>(prom_op.GetReducedDimension()));
    check_update(vector_temp, c);
  }
  check_projection(*K, prom_op.GetKr());
  check_projection(*C, prom_op.GetCr());
  check_projection(*M, prom_op.GetMr());

  // A purely real vector adds a single basis vector, with real coefficients.
  {
    const auto dim = prom_op.GetReducedDimension();
    vector_temp.Real().Randomize(seeds.back() + 1);
    vector_temp.Imag() = 0.0;
    const auto c = prom_op.UpdatePROM(vector_temp, "vec_real");
    REQUIRE(prom_op.GetReducedDimension() == dim + 1);
    REQUIRE(c.size() == static_cast<Eigen::Index>(dim + 1));
    CHECK(c.imag().cwiseAbs().maxCoeff() == 0.0);
    check_update(vector_temp, c);
    check_projection(*K, prom_op.GetKr());
    check_projection(*C, prom_op.GetCr());
    check_projection(*M, prom_op.GetMr());
  }

  // A zero vector leaves the basis unchanged, with zero coefficients.
  {
    const auto dim = prom_op.GetReducedDimension();
    vector_temp = 0.0;
    const auto c = prom_op.UpdatePROM(vector_temp, "vec_zero");
    CHECK(prom_op.GetReducedDimension() == dim);
    REQUIRE(c.size() == static_cast<Eigen::Index>(dim));
    CHECK(c.cwiseAbs().maxCoeff() == 0.0);
  }

  // Check orthogonality with s-matrix vectors.
  const auto &rom_vectors = prom_op.GetVectors();

  for (const auto &[port_idx, port_data] : space_op.GetLumpedPortOp())
  {
    for (std::size_t i = 0; i < rom_vectors.size() - 2; i++)
    {
      GridFunction E(space_op.GetNDSpace(), true);
      E = 0.0;
      E.Real().SetFromTrueDofs(rom_vectors.at(2 + i));
      E.Real().ExchangeFaceNbrData();
      E.Imag().ExchangeFaceNbrData();

      auto S = port_data.GetSParameter(E);
      CHECK_THAT(std::abs(S), WithinAbs(0.0, 1e-14));
    }
  }
}

// Checks failure mode that neighbouring ports have overlap because they share and edge
// degree of freedom. Ports in ROM must be orthogonal for conventional scattering matrix
// interpretation to be meaningful.
TEST_CASE_METHOD(palace::test::PerRankTempDir, "RomOperator-Synthesis-PortOrthogonality",
                 "[romoperator][Serial]")
{
  MPI_Comm world_comm = Mpi::World();

  size_t order = GENERATE(1UL, 2UL);  // Solver order. Should not matter but set non-default

  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_1_1_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_1_1_1_tet.msh"));

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0},
                                               {"LossTan", 0.0}})})}};

  // Port on neighbouring attributes along same direction will have overlap.
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}}),
                                  json::object({{"Index", 2},
                                                {"R", 50.0},
                                                {"Excitation", uint(0)},
                                                {"Attributes", json::array({102})},
                                                {"Direction", "+X"}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = order;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  std::size_t max_size_per_excitation = 100;

  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  CHECK_THROWS_WITH(prom_op.AddLumpedPortModesForSynthesis(),
                    Catch::Matchers::ContainsSubstring("should have exactly zero overlap"));
}

// By default (abort_if_dependent = true), adding a vector in the span of the PROM basis
// aborts. This tests checks that adding the same vector twice fails.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-UpdatePROM-LinearDependence",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};

  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"Excitation", uint(1)},
            {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                    {"Direction", "+Z"}})})}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  std::size_t max_size_per_excitation = 100;

  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);

  ComplexVector u;
  u.SetSize(space_op.GetNDSpace().GetTrueVSize());
  u.UseDevice(true);
  u.Real().Randomize(42);
  u.Imag() = 0.0;

  prom_op.UpdatePROM(u, "vec_0");
  CHECK_THROWS(prom_op.UpdatePROM(u, "vec_0_duplicate"));
}

namespace
{

// Driven problem on the 3x2x1 cube with one lumped port, for the linear dependence tests.
json DependenceTestSetup(const std::string &temp_dir, double tol, int max_its,
                         bool synthesis = false)
{
  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"Excitation", uint(1)},
            {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                    {"Direction", "+Z"}})})}})})}};
  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", synthesis},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {{"Type", "Default"},
                                    {"KSPType", "GMRES"},
                                    {"MaxIts", max_its},
                                    {"Tol", tol},
                                    {"InitialGuess", true}};
  return setup_json;
}

// Relative distance ||V c - u|| / ||u|| of the vector with coefficients c in the PROM basis
// to u.
double BasisDistance(MPI_Comm comm, const RomOperatorTest &prom_op,
                     const Eigen::VectorXcd &c, const ComplexVector &u)
{
  const auto &V = prom_op.GetVectors();
  REQUIRE(c.size() == static_cast<Eigen::Index>(V.size()));
  ComplexVector w(u.Size());
  w.UseDevice(true);
  linalg::LinearCombination(V, {c.data(), static_cast<std::size_t>(c.size())}, w);
  w -= u;
  return linalg::Norml2(comm, w) / linalg::Norml2(comm, u);
}

}  // namespace

// With abort_if_dependent = false, a sample in the span of the PROM basis (the same vector
// twice, or the same vector up to a perturbation below the linear dependence tolerance)
// leaves the PROM unchanged and returns its coefficients in the basis. A perturbation above
// the tolerance is added.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-UpdatePROM-LinearDependence-NoAbort",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  IoData iodata(DependenceTestSetup(temp_dir, 1.0e-8, 200), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 100);

  ComplexVector u(space_op.GetNDSpace().GetTrueVSize());
  u.UseDevice(true);
  u.Real().Randomize(42);
  u.Imag() = 0.0;
  REQUIRE(prom_op.UpdatePROM(u, "vec_0").size() == 1);
  REQUIRE(prom_op.GetReducedDimension() == 1);
  const Eigen::MatrixXd R0 = prom_op.GetOrthR();

  const auto c = prom_op.UpdatePROM(u, "vec_0_duplicate", false);
  CHECK(prom_op.GetReducedDimension() == 1);
  CHECK(prom_op.GetVectors().size() == 1);
  CHECK(prom_op.GetOrthR().rows() == 1);
  CHECK(prom_op.GetOrthR().cols() == 1);
  CHECK(prom_op.GetOrthR()(0, 0) == R0(0, 0));
  REQUIRE(c.size() == 1);
  CHECK_THAT(c(0).real(), WithinRel(R0(0, 0), 1.0e-12));
  CHECK(c(0).imag() == 0.0);

  // Perturbations of relative size 1e-14 and 1e-6 (the tolerance is 1e-12).
  ComplexVector w(u.Size()), v(u.Size());
  w.UseDevice(true);
  v.UseDevice(true);
  w.Real().Randomize(43);
  w.Imag() = 0.0;
  const double scale = linalg::Norml2(world_comm, u) / linalg::Norml2(world_comm, w);
  v = u;
  linalg::AXPY(1.0e-14 * scale, w, v);
  CHECK(prom_op.UpdatePROM(v, "vec_0_near_duplicate", false).size() == 1);
  CHECK(prom_op.GetReducedDimension() == 1);
  v = u;
  linalg::AXPY(1.0e-6 * scale, w, v);
  CHECK(prom_op.UpdatePROM(v, "vec_1", false).size() == 2);
  CHECK(prom_op.GetReducedDimension() == 2);
  CHECK(prom_op.GetOrthR().rows() == 2);
  CHECK(prom_op.GetOrthR()(0, 0) == R0(0, 0));
}

// With abort_if_dependent = false, a sample with one part in the span of the basis adds
// the other part only, and returns the coefficients of the sample in the updated basis,
// with the projections onto the basis for the part in its span: both when the real part is
// the dependent one, and when the real part was already appended to the basis before the
// imaginary part was found dependent (the real part is added first).
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-UpdatePROM-PartialDependence",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  IoData iodata(DependenceTestSetup(temp_dir, 1.0e-8, 200), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 100);

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u(n), x(n);
  u.UseDevice(true);
  x.UseDevice(true);
  u.Real().Randomize(42);
  u.Imag() = 0.0;
  REQUIRE(prom_op.UpdatePROM(u, "vec_0").size() == 1);
  REQUIRE(prom_op.GetReducedDimension() == 1);
  const double norm_u = prom_op.GetOrthR()(0, 0);

  // Real part in the span, new imaginary part: the imaginary part is added.
  x.Real() = u.Real();
  x.Imag().Randomize(43);
  auto c = prom_op.UpdatePROM(x, "vec_1", false);
  CHECK(prom_op.GetReducedDimension() == 2);
  CHECK(prom_op.GetOrthR().rows() == 2);
  CHECK(prom_op.GetOrthR().isUpperTriangular());
  REQUIRE(c.size() == 2);
  CHECK_THAT(c(0).real(), WithinRel(norm_u, 1.0e-12));
  CHECK(c(1).real() == 0.0);
  CHECK(c(1).imag() == prom_op.GetOrthR()(1, 1));
  CHECK(BasisDistance(world_comm, prom_op, c, x) <= 1.0e-12);

  // New real part (added first), imaginary part in the span: the real part is added.
  x.Real().Randomize(44);
  x.Imag() = u.Real();
  c = prom_op.UpdatePROM(x, "vec_2", false);
  CHECK(prom_op.GetReducedDimension() == 3);
  CHECK(prom_op.GetVectors().size() == 3);
  CHECK(prom_op.GetOrthR().rows() == 3);
  CHECK(prom_op.GetOrthR().isUpperTriangular());
  REQUIRE(c.size() == 3);
  CHECK(c(2).real() == prom_op.GetOrthR()(2, 2));
  CHECK_THAT(c(0).imag(), WithinRel(norm_u, 1.0e-12));
  CHECK(BasisDistance(world_comm, prom_op, c, x) <= 1.0e-12);

  // Then, a sample with two new parts adds both, with consistent R.
  x.Real().Randomize(45);
  x.Imag().Randomize(46);
  c = prom_op.UpdatePROM(x, "vec_3", false);
  CHECK(c.size() == 5);
  CHECK(prom_op.GetReducedDimension() == 5);
  CHECK(prom_op.GetOrthR().rows() == 5);
  CHECK(prom_op.GetOrthR().isUpperTriangular());
  CHECK(BasisDistance(world_comm, prom_op, c, x) <= 1.0e-12);
}

// In the adaptive sweep, the error of the PROM solution y at a sample (before the update)
// is ||c - y|| / ||c|| in the norm of the basis, from the coefficients c of the HDM
// solution u in the updated basis. This is the error of the expanded PROM solution against
// u also when u, or one of its parts, is in the span of the basis and not added: for the
// same sample twice, and for a sample of which only the imaginary part is added, without
// and with the synthesis weight.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-UpdatePROM-DependentSampleError",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  const bool synthesis = GENERATE(false, true);
  IoData iodata(DependenceTestSetup(temp_dir, 1.0e-8, 200, synthesis), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 10);
  if (synthesis)
  {
    prom_op.AddLumpedPortModesForSynthesis();
  }
  const auto &sample_f = iodata.solver.driven.sample_f;
  const double omega0 = sample_f.front(), omega1 = sample_f[sample_f.size() / 2];

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u0(n), u1(n), x(n), Vy(n);
  u0.UseDevice(true);
  u1.UseDevice(true);
  x.UseDevice(true);
  Vy.UseDevice(true);
  CHECK(prom_op.SolveHDM(1, omega0, u0));
  REQUIRE(prom_op.UpdatePROM(u0, "sample_0").size() > 0);
  CHECK(prom_op.SolveHDM(1, omega1, u1));
  const auto c1 = prom_op.UpdatePROM(u1, "sample_1");
  const auto dim = prom_op.GetReducedDimension();
  REQUIRE(c1.size() == static_cast<Eigen::Index>(dim));

  // Error of the PROM solution at omega (before the update) as in the adaptive sweep,
  // checked against the error of its expansion against the sample x.
  auto UpdateError = [&](double omega, std::string_view label)
  {
    prom_op.SolvePROM(omega);
    prom_op.ExpandPROMSolution(1, Vy);
    Eigen::VectorXcd d = -prom_op.GetReducedSolution(1);
    const auto c = prom_op.UpdatePROM(x, label, false);
    d.conservativeResizeLike(Eigen::VectorXcd::Zero(c.size()));
    const double error = prom_op.BasisNorm(c + d) / prom_op.BasisNorm(c);
    linalg::AXPY(-1.0, x, Vy);
    const double error_ref = linalg::Norml2(world_comm, Vy) / linalg::Norml2(world_comm, x);
    CAPTURE(synthesis, label, error, error_ref);
    CHECK(std::abs(error - error_ref) <= 1.0e-11 * std::max(1.0, error_ref));
    return c;
  };

  // The same sample again: not added, with the coefficients of the first time.
  x = u1;
  auto c = UpdateError(omega1, "sample_1_duplicate");
  CHECK(prom_op.GetReducedDimension() == dim);
  REQUIRE(c.size() == c1.size());
  CHECK((c - c1).norm() <= 1.0e-12 * c1.norm());
  CHECK(BasisDistance(world_comm, prom_op, c, x) <= 1.0e-12);

  // A sample with its real part in the span of the basis: only the imaginary part is added.
  x.Real() = u1.Real();
  x.Imag().Randomize(42);
  c = UpdateError(omega0, "sample_2");
  CHECK(prom_op.GetReducedDimension() == dim + 1);
  CHECK(BasisDistance(world_comm, prom_op, c, x) <= 1.0e-12);
}

// The adaptive sweep adds a certified sample, solved from the PROM solution, with
// abort_if_dependent = false. When the linear solver barely changes the PROM solution, or
// takes no iteration, the real or imaginary part of the solution can be in the span of the
// basis to the orthogonalization tolerance: UpdatePROM must not abort, adds the other part
// only, and returns the coefficients of the solution in the basis, with tight and loose
// linear solver tolerances.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-UpdatePROM-CertifiedSample",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  const double tol = GENERATE(1.0e-8, 1.0e-12);
  IoData iodata(DependenceTestSetup(temp_dir, tol, 500), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 10);
  const auto &sample_f = iodata.solver.driven.sample_f;
  const double omega = sample_f[sample_f.size() / 2];
  const double adaptive_tol = iodata.solver.driven.adaptive_tol;

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u(n), x0(n);
  u.UseDevice(true);
  x0.UseDevice(true);

  // First sample, from zero.
  CHECK(prom_op.SolveHDM(1, omega, u));
  REQUIRE(prom_op.UpdatePROM(u, "sample_0").size() > 0);
  const auto dim = prom_op.GetReducedDimension();
  REQUIRE(dim > 0);

  // Second sample at the same frequency, certified from the PROM solution there: the PROM
  // is exact to the solver tolerance, so the solve converges with a correction below the
  // sweep tolerance, and the sample can be in the span of the basis to the
  // orthogonalization tolerance.
  prom_op.SolvePROM(omega);
  prom_op.ExpandPROMSolution(1, x0);
  u = x0;
  const auto &ksp = prom_op.GetLinearSolver();
  const int num_it = ksp.NumTotalMultIterations();
  CHECK(prom_op.SolveHDM(1, omega, u, true, adaptive_tol));
  const int cert_it = ksp.NumTotalMultIterations() - num_it;
  REQUIRE(linalg::Norml2(world_comm, u) > 0.0);
  linalg::AXPY(-1.0, u, x0);
  CHECK(linalg::Norml2(world_comm, x0) < adaptive_tol * linalg::Norml2(world_comm, u));
  Eigen::VectorXcd c;
  CHECK_NOTHROW(c = prom_op.UpdatePROM(u, "sample_1", false));
  const bool added = (prom_op.GetReducedDimension() > dim);
  CAPTURE(cert_it, added);
  if (cert_it == 0)
  {
    // Not corrected: u is the PROM solution, in the span of the basis.
    CHECK_FALSE(added);
  }
  CHECK(prom_op.GetReducedDimension() <= dim + 2);
  CHECK(prom_op.GetOrthR().rows() ==
        static_cast<Eigen::Index>(prom_op.GetReducedDimension()));
  CHECK(BasisDistance(world_comm, prom_op, c, u) <= 1.0e-11);
}

// SolveHDM starts from zero unless asked to start from the given vector. Started from the
// solution, it returns it after zero iterations; started far from it, with a correction
// limit, it stops early and reports that it did not converge. A second solve at the same
// frequency reuses the operators.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-SolveHDM-InitialGuess",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"Excitation", uint(1)},
            {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                    {"Direction", "+Z"}})})}})})}};
  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {
      {"AdaptiveTol", 1e-3}, {"MinFreq", 2.0}, {"MaxFreq", 32.0}, {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {{"Type", "Default"},
                                    {"KSPType", "GMRES"},
                                    {"MaxIts", 200},
                                    {"Tol", 1.0e-8},
                                    {"InitialGuess", true}};

  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 10);
  const auto &ksp = prom_op.GetLinearSolver();
  const auto &sample_f = iodata.solver.driven.sample_f;
  const double omega = sample_f[sample_f.size() / 2];

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u(n), v(n);
  u.UseDevice(true);
  v.UseDevice(true);

  // From zero: the content of u is ignored.
  u.Real() = 1.0;
  u.Imag() = 0.0;
  CHECK(prom_op.SolveHDM(1, omega, u));
  const double norm_u = linalg::Norml2(world_comm, u);
  REQUIRE(norm_u > 0.0);
  v.Real() = -3.0;
  v.Imag() = 2.0;
  CHECK(prom_op.SolveHDM(1, omega, v));
  linalg::AXPY(-1.0, u, v);
  CHECK(linalg::Norml2(world_comm, v) <= 1.0e-12 * norm_u);

  // From the solution: zero iterations, the solution is returned unchanged.
  v = u;
  int num_it = ksp.NumTotalMultIterations();
  CHECK(prom_op.SolveHDM(1, omega, v, true, 1.0e-6));
  CHECK(ksp.NumTotalMultIterations() == num_it);
  linalg::AXPY(-1.0, u, v);
  CHECK(linalg::Norml2(world_comm, v) <= 1.0e-12 * norm_u);

  // From half the solution, with a correction limit: stopped after one iteration.
  v = u;
  v *= 0.5;
  num_it = ksp.NumTotalMultIterations();
  CHECK_FALSE(prom_op.SolveHDM(1, omega, v, true, 1.0e-6));
  CHECK(ksp.NumTotalMultIterations() == num_it + 1);
  CHECK(ksp.GetCorrectionExceeded());
}

// A solve started from u with a correction limit tighter than the linear solver tolerance
// (the certification of the PROM solution when "AdaptiveTol" < "Tol") is solved to this
// limit: a guess which meets the solver tolerance, but is farther than the limit from the
// solution, is corrected, and the correction stops the solve, instead of the guess being
// returned after no iteration as converged. With a looser limit, the solver tolerance is
// used. The solver tolerance is unchanged after the solves.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-SolveHDM-CertificationTolerance",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();
  constexpr double tol = 1.0e-6;
  IoData iodata(DependenceTestSetup(temp_dir, tol, 500), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  RomOperatorTest prom_op(iodata, space_op, 10);
  auto &ksp = prom_op.GetKsp();
  const auto &sample_f = iodata.solver.driven.sample_f;
  const double omega = sample_f[sample_f.size() / 2];

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector u(n), v(n);
  u.UseDevice(true);
  v.UseDevice(true);

  // Solution to 10⁻¹², far below the error of the guess below.
  REQUIRE(ksp.GetRelTol() == tol);
  ksp.SetRelTol(1.0e-12);
  REQUIRE(prom_op.SolveHDM(1, omega, u));
  ksp.SetRelTol(tol);

  // The guess (1 + 10⁻⁸) u: its relative error, and its relative preconditioned residual,
  // are about 10⁻⁸, below the solver tolerance.
  const auto set_guess = [&]()
  {
    v = u;
    v *= 1.0 + 1.0e-8;
  };

  // Limit 10⁻⁴, looser than the solver tolerance: no iteration, converged.
  set_guess();
  int num_it = ksp.NumTotalMultIterations();
  CHECK(prom_op.SolveHDM(1, omega, v, true, 1.0e-4));
  CHECK(ksp.NumTotalMultIterations() == num_it);
  CHECK(ksp.GetRelTol() == tol);

  // Limit 10⁻¹⁰, tighter than the solver tolerance and than the error of the guess: solved
  // to 10⁻¹⁰, so corrected by about 10⁻⁸, which exceeds the limit.
  set_guess();
  num_it = ksp.NumTotalMultIterations();
  CHECK_FALSE(prom_op.SolveHDM(1, omega, v, true, 1.0e-10));
  CHECK(ksp.NumTotalMultIterations() > num_it);
  CHECK(ksp.GetCorrectionExceeded());
  CHECK(ksp.GetRelTol() == tol);

  // From zero, the solver tolerance is used.
  num_it = ksp.NumTotalMultIterations();
  CHECK(prom_op.SolveHDM(1, omega, v));
  CHECK(ksp.NumTotalMultIterations() > num_it);
  CHECK(ksp.GetRelTol() == tol);
}

// SolveHDM keeps the operators and the preconditioner of its last solve for the next solve
// at the same frequency only, and ReleaseHDMOperators releases them: the sample solved
// again from zero after its certification was stopped (as SweepAdaptive does), the solve at
// another frequency, and the solve after the release all give the solution of a RomOperator
// which solved nothing before.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-SolveHDM-OperatorReuse",
                 "[romoperator][Serial][Parallel]")
{
  // Solved to 10⁻¹², so that solutions with separately set up preconditioners agree to well
  // below the 10⁻⁸ checked.
  MPI_Comm world_comm = Mpi::World();
  IoData iodata(DependenceTestSetup(temp_dir, 1.0e-12, 500), false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  const auto &sample_f = iodata.solver.driven.sample_f;
  const double omega_a = sample_f[sample_f.size() / 3];
  const double omega_b = sample_f[(2 * sample_f.size()) / 3];

  const auto n = space_op.GetNDSpace().GetTrueVSize();
  auto make_vector = [n]()
  {
    ComplexVector v(n);
    v.UseDevice(true);
    v = 0.0;
    return v;
  };
  auto rel_diff = [&](const ComplexVector &x, const ComplexVector &y)
  {
    ComplexVector d(x);
    d.UseDevice(true);
    linalg::AXPY(-1.0, y, d);
    return linalg::Norml2(world_comm, d) / linalg::Norml2(world_comm, y);
  };

  // Reference solutions, each from a RomOperator which solved nothing before.
  ComplexVector u_a = make_vector(), u_b = make_vector(), v = make_vector();
  {
    RomOperatorTest ref(iodata, space_op, 10);
    REQUIRE(ref.SolveHDM(1, omega_a, u_a));
  }
  {
    RomOperatorTest ref(iodata, space_op, 10);
    REQUIRE(ref.SolveHDM(1, omega_b, u_b));
  }

  RomOperatorTest prom_op(iodata, space_op, 10);
  CHECK(prom_op.GetHDMOperator() == nullptr);

  // Stopped by the correction limit, then solved again from zero: the operators are reused.
  v = u_a;
  v *= 0.5;
  CHECK_FALSE(prom_op.SolveHDM(1, omega_a, v, true, 1.0e-6));
  const auto *A_a = prom_op.GetHDMOperator();
  REQUIRE(A_a != nullptr);
  CHECK(prom_op.GetHDMOmega() == omega_a);
  CHECK(prom_op.SolveHDM(1, omega_a, v));
  CHECK(prom_op.GetHDMOperator() == A_a);
  CHECK(rel_diff(v, u_a) <= 1.0e-8);

  // Another frequency: the operators are those of this frequency.
  CHECK(prom_op.SolveHDM(1, omega_b, v));
  CHECK(prom_op.GetHDMOmega() == omega_b);
  CHECK(rel_diff(v, u_b) <= 1.0e-8);

  // Released, and rebuilt by the next solve at the same frequency.
  prom_op.ReleaseHDMOperators();
  CHECK(prom_op.GetHDMOperator() == nullptr);
  CHECK(prom_op.SolveHDM(1, omega_b, v));
  CHECK(prom_op.GetHDMOperator() != nullptr);
  CHECK(prom_op.GetHDMOmega() == omega_b);
  CHECK(rel_diff(v, u_b) <= 1.0e-8);
}

// Tests that the per-port "IncludeInSynthesis" flag suppresses basis vector injection
// without changing the boundary condition. We reuse the cube_3_2_1 two-port setup, but
// flip port 1 (non-excited, RLC) to "IncludeInSynthesis": false, leaving port 2 (excited,
// pure R) included. Expectation: AddLumpedPortModesForSynthesis adds exactly one basis
// vector (port 2), the resulting normalized circuit matrices are 1x1, and the diagonal
// resistance entry equals 1 / port2_ref_R.
//
// This also exercises CalculateNormalizedPROMMatrices in the regime where the number of
// port-mode rows (1) is less than the total lumped port count (2): the port-scaling loop
// must iterate over included ports only, not GetLumpedPortOp().Size(), or it writes past
// the end of v_conc.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-Synthesis-IncludeInSynthesis-Skip",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  const double L0 = 1.0e-6;
  const double Lc = 7.0;

  const double port1_ref_R = 75.0;
  const double port1_ref_C = 5.5e-5;
  const double port1_ref_L = 1.486e-20;
  const double port2_ref_R = 50.0;

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 2}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"L0", L0},
                         {"Lc", Lc},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};

  // Port 1 is the passive RLC port flagged out of synthesis; port 2 stays in.
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array(
           {json::object(
                {{"Index", 1},
                 {"R", port1_ref_R},
                 {"C", port1_ref_C},
                 {"L", port1_ref_L},
                 {"Excitation", uint(0)},
                 {"IncludeInSynthesis", false},
                 {"Elements",
                  json::array({json::object(
                                   {{"Attributes", json::array({1})}, {"Direction", "+Y"}}),
                               json::object({{"Attributes", json::array({9, 11})},
                                             {"Direction", "+Y"}}),
                               json::object({{"Attributes", json::array({2, 6, 10})},
                                             {"Direction", "+X"}})})}}),
            json::object(
                {{"Index", 2},
                 {"R", port2_ref_R},
                 {"Excitation", uint(1)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                         {"Direction", "+Z"}})})}})})}};

  setup_json["Solver"] = json::object();
  setup_json["Solver"]["Order"] = 1UL;
  setup_json["Solver"]["Device"] = "CPU";
  setup_json["Solver"]["Driven"] = {{"AdaptiveTol", 1e-3},
                                    {"AdaptiveCircuitSynthesis", true},
                                    {"MinFreq", 2.0},
                                    {"MaxFreq", 32.0},
                                    {"FreqStep", 1.0}};
  setup_json["Solver"]["Linear"] = {
      {"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}};

  IoData iodata(setup_json, false);

  // The flag is parsed onto the config struct.
  const auto &cfg_ports = iodata.boundaries.lumpedport;
  REQUIRE(cfg_ports.size() == 2);
  CHECK(cfg_ports.at(1).include_in_synthesis == false);
  CHECK(cfg_ports.at(2).include_in_synthesis == true);

  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  std::size_t max_size_per_excitation = 100;
  RomOperatorTest prom_op(iodata, space_op, max_size_per_excitation);
  prom_op.AddLumpedPortModesForSynthesis();

  // Only port 2 contributed a basis vector.
  REQUIRE(prom_op.GetReducedDimension() == 1);

  const auto matrices = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  const auto &inductance_L_inv = matrices.L_inv;
  const auto &resistance_R_inv = matrices.R_inv;
  const auto &capacitance_C = matrices.C;
  REQUIRE(inductance_L_inv->rows() == 1);
  REQUIRE(resistance_R_inv->rows() == 1);
  REQUIRE(capacitance_C->rows() == 1);

  // The remaining diagonal entry corresponds to port 2 (R = 50 Ohm, no L, no C).
  CHECK_THAT(std::real((*resistance_R_inv)(0, 0)), WithinRel(1.0 / port2_ref_R));
  CHECK_THAT(std::imag((*resistance_R_inv)(0, 0)), WithinAbs(0.0, 1e-15));
}

// The exported interface capacitance matrices rom-interface-<index>-C.csv are the
// interface dielectric energy forms projected onto the PROM basis with the normalization of
// rom-C-re.csv: for the synthesized coordinates y of any field, yᴴ C_Γ y / yᴴ Re{C} y is
// the surface participation ratio U / (¼ uᴴ M u) of the prolonged field u, a peak phasor
// with the time-averaged interface energy U.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-Synthesis-InterfaceCapacitanceMatrices",
                 "[romoperator][postoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  // Exterior boundaries of the 3 × 2 × 1 cube (the lumped port is on attribute 14).
  std::vector<int> walls;
  for (int attr = 1; attr <= 22; attr++)
  {
    if (attr != 14)
    {
      walls.push_back(attr);
    }
  }
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path}, {"Refinement", json::object({})}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1, 2, 3})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}}),
                                 json::object({{"Attributes", json::array({4, 5, 6})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 9.0}})})}};
  json dielectric = json::array();
  int index = 0;
  for (const auto *type : {"Default", "MA", "MS", "SA"})
  {
    dielectric.push_back(json::object({{"Index", ++index},
                                       {"Type", type},
                                       {"Attributes", walls},
                                       {"Thickness", 2.0e-3},
                                       {"Permittivity", 10.0},
                                       {"LossTan", 1.0e-3}}));
  }
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"L", 1.0e-9},
            {"C", 1.0e-12},
            {"Excitation", uint(1)},
            {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                    {"Direction", "+Z"}})})}})})},
      {"Postprocessing", {{"Dielectric", dielectric}}}};
  setup_json["Solver"] = {
      {"Order", 2UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3},
        {"AdaptiveCircuitSynthesis", true},
        {"MinFreq", 2.0},
        {"MaxFreq", 32.0},
        {"FreqStep", 30.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  DrivenPostOperatorTest post_op(iodata, space_op);

  // A port mode and two field samples (their real and imaginary parts).
  rom_op.AddLumpedPortModesForSynthesis();
  for (int seed : {42, 44})
  {
    ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
    sample.UseDevice(true);
    linalg::SetRandom(comm, sample.Real(), seed);
    linalg::SetRandom(comm, sample.Imag(), seed + 1);
    rom_op.UpdatePROM(sample, fmt::format("sample_{:d}", seed));
  }
  const auto &V = rom_op.GetBasis();
  const long n = static_cast<long>(V.size());
  REQUIRE(n == 5);

  const auto forms = post_op.GetInterfaceElectricFieldEnergyMatrices(V);
  REQUIRE(forms.size() == 4);
  rom_op.PrintPROMMatrices(iodata.units, temp_dir, forms);
  Mpi::Barrier(comm);

  auto read_matrix = [](const fs::path &p, std::vector<std::string> &labels)
  {
    REQUIRE(fs::is_regular_file(p));
    TableWithCSVFile file(p.string(), /*load_existing_file=*/true);
    auto &table = file.table;
    const long cols = static_cast<long>(table.n_cols()),
               rows = static_cast<long>(table.n_rows());
    Eigen::MatrixXd M(rows, cols);
    labels.clear();
    for (long c = 0; c < cols; c++)
    {
      labels.push_back(table[c].header_text);
      for (long r = 0; r < rows; r++)
      {
        M(r, c) = table[c].data[r];
      }
    }
    return M;
  };
  std::vector<std::string> labels_C;
  const Eigen::MatrixXd C = read_matrix(fs::path(temp_dir) / "rom-C-re.csv", labels_C);
  REQUIRE(C.rows() == n);

  // Port-row scaling and capacitance unit of rom-C-re.csv.
  Eigen::VectorXd d = Eigen::VectorXd::Ones(n);
  d(0) = rom_op.GetOrthR()(0, 0);
  const double unit_farad = iodata.units.GetScaleFactor<Units::ValueType::CAPACITANCE>();
  const auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  // Complex synthesized coordinates y, reduced-basis coordinates z = D y, and the prolonged
  // field u = V z.
  Eigen::VectorXcd y(n);
  for (long j = 0; j < n; j++)
  {
    y(j) = std::complex<double>(1.0 / (j + 1.0), 0.3 - 0.1 * j);
  }
  const Eigen::VectorXcd z = d.cast<std::complex<double>>().cwiseProduct(y);
  ComplexVector u(space_op.GetNDSpace().GetTrueVSize()), Mu(u.Size());
  u.UseDevice(true);
  Mu.UseDevice(true);
  u = 0.0;
  for (long j = 0; j < n; j++)
  {
    linalg::AXPY(z(j).real(), V[j], u.Real());
    linalg::AXPY(z(j).imag(), V[j], u.Imag());
  }
  M->Mult(u, Mu);
  const double energy_elec = 0.25 * linalg::Dot(comm, Mu, u).real();
  REQUIRE(energy_elec > 0.0);
  // The interface energies of the prolonged field, which do not depend on the frequency
  // for constant materials.
  SurfacePostOperator surf_post_op(iodata, space_op.GetMaterialOp(),
                                   space_op.GetH1Space().Get(),
                                   space_op.GetNDSpace().Get());
  GridFunction E(space_op.GetNDSpace(), true);
  E.Real().SetFromTrueDofs(u.Real());
  E.Imag().SetFromTrueDofs(u.Imag());

  int num_nonzero = 0;
  for (const auto &[idx, B] : forms)
  {
    CAPTURE(idx);
    std::vector<std::string> labels;
    const Eigen::MatrixXd C_i = read_matrix(
        fs::path(temp_dir) / fmt::format("rom-interface-{:d}-C.csv", idx), labels);
    CHECK(labels == labels_C);
    const Eigen::MatrixXd C_ref = 2.0 * unit_farad * d.asDiagonal() * B * d.asDiagonal();
    CHECK((C_i - C_ref).norm() <= 1.0e-14 * C_ref.norm());
    CHECK((C_i - C_i.transpose()).norm() == 0.0);

    const double p_matrix =
        (y.adjoint() * C_i * y)(0).real() / (y.adjoint() * C * y)(0).real();
    const double p_direct =
        surf_post_op.GetInterfaceElectricFieldEnergy(idx, E, 1.0) / energy_elec;
    CAPTURE(p_matrix, p_direct);
    CHECK(p_matrix == Catch::Approx(p_direct).epsilon(1.0e-10));
    num_nonzero += (p_direct > 0.0);
  }
  CHECK(num_nonzero == 4);
}

// Two lumped ports, one with L and C, and a domain energy postprocessing region, on the
// 3x2x1 tetrahedral cube mesh, for the reduced postprocessing, with the given permittivity
// of its material.
json ReducedPostprocessingSetup(const fs::path &temp_dir, const json &permittivity = 1.0)
{
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", permittivity},
                                  {"LossTan", 0.0}})})},
      {"Postprocessing",
       {{"Energy",
         json::array({json::object(
             {{"Index", 1}, {"Attributes", json::array({1, 2, 3, 4, 5, 6})}})})}}}};
  setup_json["Boundaries"] = {
      {"LumpedPort",
       json::array(
           {json::object(
                {{"Index", 1},
                 {"R", 50.0},
                 {"L", 1.0e-9},
                 {"C", 1.0e-12},
                 {"Excitation", uint(1)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({14})},
                                                         {"Direction", "+Z"}})})}}),
            json::object(
                {{"Index", 2},
                 {"R", 50.0},
                 {"Excitation", uint(2)},
                 {"Elements", json::array({json::object({{"Attributes", json::array({1})},
                                                         {"Direction", "+Y"}})})}})})}};
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3},
        {"AdaptiveCircuitSynthesis", true},
        {"MinFreq", 2.0},
        {"MaxFreq", 32.0},
        {"FreqStep", 30.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

  return setup_json;
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced postprocessing matches full evaluation",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  // The reduced energies are the full ones, also with the stored energy of a Drude.
  // Without a stored electric energy (a RAW relaxation compensated by a typed Debye), both
  // omit the electric energies alike.
  enum class Model
  {
    CONSTANT,
    DRUDE,
    COMPENSATED
  };
  const auto model = GENERATE(Model::CONSTANT, Model::DRUDE, Model::COMPENSATED);
  CAPTURE(static_cast<int>(model));
  const bool stored = (model != Model::COMPENSATED);
  json permittivity = 1.0;
  if (model == Model::DRUDE)
  {
    // ω_p = 2π 10 GHz, γ = 2π 5 GHz.
    permittivity = {{"HighFrequencyLimit", 1.0},
                    {"Terms", json::array({{{"Type", "Drude"},
                                            {"PlasmaFrequency", 10.0},
                                            {"CollisionFrequency", 5.0}}})}};
  }
  else if (model == Model::COMPENSATED)
  {
    // Debye 1/(1 + sτ) and RAW −0.5/(s + 1/τ), τ = 10 ps.
    permittivity = {
        {"HighFrequencyLimit", 1.0},
        {"Terms",
         json::array(
             {{{"Type", "Debye"}, {"DeltaPermittivity", 1.0}, {"RelaxationTime", 0.01}},
              {{"Type", "PoleResidue"}, {"Pole", -1.0e11}, {"Residue", -0.5e11}}})}};
  }
  MPI_Comm comm = Mpi::World();
  json setup_json = ReducedPostprocessingSetup(temp_dir, permittivity);
  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperator rom_op(iodata, space_op, 4);
  DrivenPostOperatorTest post_op(iodata, space_op);

  REQUIRE(iodata.solver.driven.sample_f.size() == 2);
  const double omega_front = iodata.solver.driven.sample_f.front();
  const double omega_back = iodata.solver.driven.sample_f.back();
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  linalg::SetRandom(comm, sample.Real(), 42);
  linalg::SetRandom(comm, sample.Imag(), 43);
  rom_op.UpdatePROM(sample, "postprocessing_sample");
  REQUIRE(rom_op.GetReducedDimension() > 0);
  post_op.ConfigureReducedPostprocessing(rom_op, iodata.solver.driven.adaptive_tol);
  REQUIRE(post_op.HasReducedPostprocessing());

  // Construct a deterministic complex field directly in the PROM space. This gives both
  // electric and magnetic forms meaningful magnitude while testing the exact reduced
  // quadratic forms independently of PROM solution accuracy.
  Eigen::VectorXcd y(rom_op.GetReducedDimension());
  for (long i = 0; i < y.size(); i++)
  {
    y(i) = std::complex<double>(1.0 / (i + 1.0), -0.25 / (i + 1.0));
  }
  ComplexVector E(space_op.GetNDSpace().GetTrueVSize());
  E.UseDevice(true);
  E = 0.0;
  for (std::size_t i = 0; i < rom_op.GetBasis().size(); i++)
  {
    linalg::AXPY(y(i).real(), rom_op.GetBasis()[i], E.Real());
    linalg::AXPY(y(i).imag(), rom_op.GetBasis()[i], E.Imag());
  }
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y, rom_op);
  REQUIRE(post_op.HasReducedPostprocessing());

  ComplexVector B(space_op.GetCurlMatrix().Height());
  B.UseDevice(true);
  space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
  space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
  B *= -1.0 / (std::complex<double>(0.0, 1.0) * omega_back);
  post_op.MeasureAndPrintAll(1, 1, E, B, omega_back);
  const Measurement full = post_op.GetMeasurement();
  REQUIRE_FALSE(post_op.ReducedPostprocessingNeedsField(1));
  post_op.MeasureAndPrintReduced(1, 1, nullptr, omega_back, y, rom_op);
  const Measurement reduced = post_op.GetMeasurement();
  REQUIRE(post_op.HasReducedPostprocessing());

  auto close = [](double x, double y)
  {
    const double scale = std::max(std::abs(x), std::abs(y));
    return scale == 0.0 || std::abs(x - y) <= 1.0e-10 * scale;
  };
  auto check_complex = [](std::complex<double> x, std::complex<double> y)
  {
    const double scale = std::max(std::abs(x), std::abs(y));
    CHECK((scale == 0.0 || std::abs(x - y) <= 1.0e-10 * scale));
  };

  REQUIRE(space_op.GetMaterialOp().HasStoredElectricEnergy() == stored);
  CHECK(full.domain_E_field_energy_all.has_value() == stored);
  CHECK(reduced.domain_E_field_energy_all.has_value() == stored);
  CAPTURE(reduced.domain_H_field_energy_all, full.domain_H_field_energy_all);
  if (stored)
  {
    CAPTURE(*reduced.domain_E_field_energy_all, *full.domain_E_field_energy_all);
    REQUIRE(*full.domain_E_field_energy_all > 0.0);
    CHECK(close(*reduced.domain_E_field_energy_all, *full.domain_E_field_energy_all));
  }
  CHECK(close(reduced.domain_E_field_energy_material, full.domain_E_field_energy_material));
  REQUIRE(full.domain_H_field_energy_all > 0.0);
  CHECK(close(reduced.domain_H_field_energy_all, full.domain_H_field_energy_all));
  REQUIRE(full.domain_E_field_energy_i.size() == 1);
  REQUIRE(full.domain_H_field_energy_i.size() == 1);
  REQUIRE(reduced.domain_E_field_energy_i.size() == full.domain_E_field_energy_i.size());
  REQUIRE(reduced.domain_H_field_energy_i.size() == full.domain_H_field_energy_i.size());
  for (std::size_t i = 0; i < full.domain_E_field_energy_i.size(); i++)
  {
    const auto &r = reduced.domain_E_field_energy_i[i],
               &f = full.domain_E_field_energy_i[i];
    CHECK(r.idx == f.idx);
    // The postprocessing domain contains the material.
    CHECK(f.energy.has_value() == stored);
    CHECK(r.energy.has_value() == stored);
    CHECK(f.participation_ratio.has_value() == stored);
    CHECK(r.participation_ratio.has_value() == stored);
    if (stored)
    {
      CHECK(close(r.energy.value(), f.energy.value()));
      CHECK(close(r.participation_ratio.value(), f.participation_ratio.value()));
    }
    CHECK(reduced.domain_H_field_energy_i[i].idx == full.domain_H_field_energy_i[i].idx);
    CHECK(close(reduced.domain_H_field_energy_i[i].energy.value(),
                full.domain_H_field_energy_i[i].energy.value()));
    CHECK(close(reduced.domain_H_field_energy_i[i].participation_ratio.value(),
                full.domain_H_field_energy_i[i].participation_ratio.value()));
  }

  REQUIRE(full.lumped_port_vi.size() == 2);
  REQUIRE(reduced.lumped_port_vi.size() == full.lumped_port_vi.size());
  REQUIRE(full.lumped_port_vi.at(1).inductor_energy > 0.0);
  REQUIRE(full.lumped_port_vi.at(1).capacitor_energy > 0.0);
  for (const auto &[idx, full_port] : full.lumped_port_vi)
  {
    CAPTURE(idx);
    const auto &reduced_port = reduced.lumped_port_vi.at(idx);
    REQUIRE(std::abs(full_port.V) > 0.0);
    REQUIRE(std::abs(full_port.I) > 0.0);
    REQUIRE(std::abs(full_port.S) > 0.0);
    check_complex(reduced_port.V, full_port.V);
    check_complex(reduced_port.I, full_port.I);
    REQUIRE(reduced_port.I_RLC.size() == full_port.I_RLC.size());
    for (std::size_t i = 0; i < full_port.I_RLC.size(); i++)
    {
      check_complex(reduced_port.I_RLC[i], full_port.I_RLC[i]);
    }
    check_complex(reduced_port.S, full_port.S);
    CHECK(close(reduced_port.inductor_energy, full_port.inductor_energy));
    CHECK(close(reduced_port.capacitor_energy, full_port.capacitor_energy));
  }
  // Pin the physical convention independently of the shared implementation, at the
  // excited port.
  {
    const auto &port = space_op.GetLumpedPortOp().GetPort(1);
    const auto &full_port = full.lumped_port_vi.at(1);
    check_complex(full_port.S,
                  full_port.V / std::sqrt(2.0 * port.GetExcitationRefResistance()) - 1.0);
    CHECK(close(full_port.capacitor_energy, 0.25 * port.C * std::norm(full_port.V)));
    CHECK(close(full_port.inductor_energy, 0.25 * port.L * std::norm(full_port.I_RLC[1])));
  }
  CHECK(close(reduced.lumped_port_inductor_energy, full.lumped_port_inductor_energy));
  CHECK(close(reduced.lumped_port_capacitor_energy, full.lumped_port_capacitor_energy));

  // A NaN in the reduced solution reverts to full postprocessing.
  const double tol = iodata.solver.driven.adaptive_tol;
  post_op.ConfigureReducedPostprocessing(rom_op, tol);
  Eigen::VectorXcd y_nan = y;
  y_nan(0) = std::numeric_limits<double>::quiet_NaN();
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y_nan, rom_op);
  CHECK_FALSE(post_op.HasReducedPostprocessing());

  // Each excitation is checked at its first point: after the first passes, a field
  // inconsistent with the reduced solution of the second reverts.
  post_op.ConfigureReducedPostprocessing(rom_op, tol);
  post_op.MeasureAndPrintReduced(1, 0, &E, omega_front, y, rom_op);
  REQUIRE(post_op.HasReducedPostprocessing());
  REQUIRE(post_op.ReducedPostprocessingNeedsField(2));
  E *= 1.001;
  post_op.MeasureAndPrintReduced(2, 0, &E, omega_front, y, rom_op);
  CHECK_FALSE(post_op.HasReducedPostprocessing());
}

// Driven setup on cube_mesh_3_2_1_tet.msh for the excitation weight tests: an excited RLC
// lumped port 1, with the excitation weight a = amplitude exp(i phase) when it is not the
// default, and a passive resistive lumped port 2.
json WeightedExcitationSetup(const fs::path &output_dir, double amplitude, double phase)
{
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  json port_1 = json::object(
      {{"Index", 1},
       {"R", 50.0},
       {"L", 1.0e-9},
       {"C", 1.0e-12},
       {"Excitation", uint(1)},
       {"Elements", json::array({json::object(
                        {{"Attributes", json::array({14})}, {"Direction", "+Z"}})})}});
  if (amplitude != 1.0 || phase != 0.0)
  {
    port_1["ExcitationAmplitude"] = amplitude;
    port_1["ExcitationPhase"] = phase;
  }
  json port_2 = json::object(
      {{"Index", 2},
       {"R", 75.0},
       {"Elements", json::array({json::object(
                        {{"Attributes", json::array({1})}, {"Direction", "+Y"}})})}});

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", output_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path},
                         {"Refinement", json::object({})},
                         {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {{"LumpedPort", json::array({port_1, port_2})}};
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 2.0}, {"MaxFreq", 32.0}, {"FreqStep", 30.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};
  return setup_json;
}

// The S-parameters are normalized by the weight of the driven port: for the same field E,
// the configuration with the weight a, measured on the field aE, must give the same S on
// every port as the configuration without a weight measured on E, and the port voltages
// scaled by a. A division by the weight that is dropped or applied twice fails this test.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "PostOperator S-parameters are normalized by the excitation weight",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const double amplitude = 2.0;
  const double phase = 30.0;
  const std::complex<double> a = std::polar(amplitude, phase * std::numbers::pi / 180.0);

  // Measure the field scale * E, with B from the curl of that field, on the configuration
  // with the given excitation weight. The random field E is the same for all weights.
  auto measure = [&](double amp, double ph, std::complex<double> scale)
  {
    IoData iodata(WeightedExcitationSetup(temp_dir, amp, ph), false);
    auto mesh = LoadScaleParMesh2(iodata, comm);
    SpaceOperator space_op(iodata, mesh);
    DrivenPostOperatorTest post_op(iodata, space_op);

    REQUIRE(iodata.solver.driven.sample_f.size() == 2);
    const double omega = iodata.solver.driven.sample_f.back();
    ComplexVector E(space_op.GetNDSpace().GetTrueVSize());
    E.UseDevice(true);
    linalg::SetRandom(comm, E.Real(), 42);
    linalg::SetRandom(comm, E.Imag(), 43);
    linalg::SetSubVector(E.Real(), space_op.GetNDDbcTDofLists().back(), 0.0);
    linalg::SetSubVector(E.Imag(), space_op.GetNDDbcTDofLists().back(), 0.0);
    E *= scale;
    ComplexVector B(space_op.GetCurlMatrix().Height());
    B.UseDevice(true);
    space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
    space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
    B *= -1.0 / (std::complex<double>(0.0, 1.0) * omega);
    post_op.MeasureAndPrintAll(1, 0, E, B, omega);
    return std::make_pair(post_op.GetMeasurement(),
                          space_op.GetLumpedPortOp().GetPort(1).GetExcitationWeight());
  };

  const auto unweighted = measure(1.0, 0.0, 1.0);
  const auto weighted = measure(amplitude, phase, a);
  CHECK(unweighted.second == std::complex<double>(1.0, 0.0));
  CHECK(std::abs(weighted.second - a) <= 1.0e-14);

  const auto &ports = unweighted.first.lumped_port_vi;
  REQUIRE(ports.size() == 2);
  REQUIRE(weighted.first.lumped_port_vi.size() == ports.size());
  for (const auto &[idx, port] : ports)
  {
    const int port_idx = idx;
    CAPTURE(port_idx);
    const auto &port_w = weighted.first.lumped_port_vi.at(port_idx);
    REQUIRE(std::abs(port.V) > 0.0);
    REQUIRE(std::abs(port.S) > 0.0);

    // V is linear in the field, so it scales by a.
    CHECK(std::abs(port_w.V - a * port.V) <= 1.0e-10 * std::abs(a * port.V));

    // S does not change with the weight.
    CHECK(std::abs(port_w.S - port.S) <=
          1.0e-10 * std::max(1.0, std::max(std::abs(port_w.S), std::abs(port.S))));
  }
}

// The adaptive sweep with a complex RHS1: for the same random basis, the PROM solution of
// the configuration with the excitation weight a is a times the solution without a weight.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator PROM solution scales with the excitation weight",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const double amplitude = 2.0;
  const double phase = 30.0;
  const std::complex<double> a = std::polar(amplitude, phase * std::numbers::pi / 180.0);

  auto solve = [&](double amp, double ph)
  {
    IoData iodata(WeightedExcitationSetup(temp_dir, amp, ph), false);
    auto mesh = LoadScaleParMesh2(iodata, comm);
    SpaceOperator space_op(iodata, mesh);
    RomOperator rom_op(iodata, space_op, 4);

    REQUIRE(iodata.solver.driven.sample_f.size() == 2);
    const double omega = iodata.solver.driven.sample_f.back();
    ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
    sample.UseDevice(true);
    for (int seed : {42, 44})
    {
      linalg::SetRandom(comm, sample.Real(), seed);
      linalg::SetRandom(comm, sample.Imag(), seed + 1);
      linalg::SetSubVector(sample.Real(), space_op.GetNDDbcTDofLists().back(), 0.0);
      linalg::SetSubVector(sample.Imag(), space_op.GetNDDbcTDofLists().back(), 0.0);
      rom_op.UpdatePROM(sample, fmt::format("weight_sample_{:d}", seed));
    }
    REQUIRE(rom_op.GetReducedDimension() > 0);

    ComplexVector u(space_op.GetNDSpace().GetTrueVSize());
    u.UseDevice(true);
    rom_op.SolvePROM(omega);
    rom_op.ExpandPROMSolution(1, u);
    return u;
  };

  const auto unweighted = solve(1.0, 0.0);
  const auto weighted = solve(amplitude, phase);

  auto check_scaled = [&](const ComplexVector &u_weighted, const ComplexVector &u)
  {
    ComplexVector diff(u_weighted);
    ComplexVector scaled(u);
    scaled *= a;
    diff.Add(std::complex<double>(-1.0, 0.0), scaled);
    const double norm_ref = linalg::Norml2(comm, scaled);
    const double norm_diff = linalg::Norml2(comm, diff);
    CAPTURE(norm_ref, norm_diff);
    REQUIRE(norm_ref > 0.0);
    CHECK(norm_diff <= 1.0e-8 * norm_ref);
  };
  check_scaled(weighted, unweighted);
}

// The reduced magnetic energy of a field V y with large, nearly cancelling coefficients y
// (V y = b for the basis {a, a + δ b}), which loses all accuracy in the quadratic form of
// VᵀCᵀMCV, and with a basis vector of zero curl, which is dropped from the factorization of
// C V. The basis is set directly, since the PROM basis is orthonormal.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced magnetic energy with cancelling coefficients",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  json setup_json = ReducedPostprocessingSetup(temp_dir);
  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  DrivenPostOperatorTest post_op(iodata, space_op);
  const double omega = iodata.solver.driven.sample_f.front();
  const int size = space_op.GetNDSpace().GetTrueVSize();

  // A discrete gradient of integer nodal values has exactly zero curl.
  Vector phi(space_op.GetH1Space().GetTrueVSize()), a(size), b(size), g(size);
  for (auto *x : {&phi, &a, &b, &g})
  {
    x->UseDevice(true);
  }
  linalg::SetRandom(comm, phi, 41);
  {
    auto *p = phi.HostReadWrite();
    for (int i = 0; i < phi.Size(); i++)
    {
      p[i] = std::floor(4.0 * p[i]) - 2.0;
    }
  }
  space_op.GetGradMatrix().Mult(phi, g);
  REQUIRE(linalg::Norml2(comm, g) > 0.0);
  linalg::SetRandom(comm, a, 42);
  linalg::SetRandom(comm, b, 43);

  // Measure the field V y with the reduced postprocessing, after the check against the full
  // evaluation at the first point, and with the full evaluation. The check is not part of
  // the test: the reduced electric energy is the quadratic form of VᵀMV, which cancels
  // like the magnetic one did with VᵀCᵀMCV.
  auto measure = [&](std::vector<Vector> basis, const Eigen::VectorXcd &y)
  {
    ComplexVector E(size);
    E.UseDevice(true);
    E = 0.0;
    for (std::size_t i = 0; i < basis.size(); i++)
    {
      linalg::AXPY(y(i).real(), basis[i], E.Real());
      linalg::AXPY(y(i).imag(), basis[i], E.Imag());
    }
    rom_op.SetBasis(std::move(basis));
    post_op.ConfigureReducedPostprocessing(rom_op, 1.0);
    REQUIRE(post_op.HasReducedPostprocessing());
    post_op.MeasureAndPrintReduced(1, 0, &E, omega, y, rom_op);
    REQUIRE(post_op.HasReducedPostprocessing());
    post_op.MeasureAndPrintReduced(1, 1, nullptr, omega, y, rom_op);
    const Measurement reduced = post_op.GetMeasurement();
    ComplexVector B(space_op.GetCurlMatrix().Height());
    B.UseDevice(true);
    space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
    space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    post_op.MeasureAndPrintAll(1, 1, E, B, omega);
    return std::make_pair(reduced, post_op.GetMeasurement());
  };

  SECTION("Cancelling coefficients")
  {
    // δ = 2⁻²⁰ and coefficients ±s/δ with s a sum of powers of two, so that V y = s b + t g
    // up to the round-off in a + δ b, and the full evaluation is accurate.
    const double delta = std::ldexp(1.0, -20);
    const std::complex<double> s(1.0, -0.5), t(0.75, 0.5);
    Vector a_db(a);
    a_db.UseDevice(true);
    a_db.Add(delta, b);
    Eigen::VectorXcd y(3);
    y << s / delta, t, -s / delta;
    const auto [reduced, full] = measure({a, g, a_db}, y);

    // The round-off of R y is about ε/δ ≈ 2e-10 relative (that of the quadratic form of
    // VᵀCᵀMCV, ε/δ² ≈ 1e-4, as for the electric energy here).
    CAPTURE(reduced.domain_H_field_energy_all, full.domain_H_field_energy_all);
    REQUIRE(full.domain_H_field_energy_all > 0.0);
    CHECK_THAT(reduced.domain_H_field_energy_all,
               WithinRel(full.domain_H_field_energy_all, 1.0e-9));
    REQUIRE(reduced.domain_H_field_energy_i.size() == 1);
    REQUIRE(full.domain_H_field_energy_i.size() == 1);
    CHECK_THAT(reduced.domain_H_field_energy_i[0].energy.value(),
               WithinRel(full.domain_H_field_energy_i[0].energy.value(), 1.0e-9));
  }

  SECTION("Zero curl")
  {
    // Only the vector of zero curl: the factorization of C V has no columns.
    Eigen::VectorXcd y(1);
    y << std::complex<double>(0.75, 0.5);
    const auto [reduced, full] = measure({g}, y);

    CAPTURE(reduced.domain_E_field_energy_all.value(),
            full.domain_E_field_energy_all.value(), reduced.domain_H_field_energy_all,
            full.domain_H_field_energy_all);
    REQUIRE(full.domain_E_field_energy_all.value() > 0.0);
    CHECK_THAT(reduced.domain_E_field_energy_all.value(),
               WithinRel(full.domain_E_field_energy_all.value(), 1.0e-10));
    CHECK(reduced.domain_H_field_energy_all == 0.0);
    CHECK(full.domain_H_field_energy_all <=
          1.0e-20 * full.domain_E_field_energy_all.value());
    REQUIRE(reduced.domain_H_field_energy_i.size() == 1);
    CHECK(reduced.domain_H_field_energy_i[0].energy == 0.0);
  }
}

// The 3 cm x 2 cm x 1 cm vacuum box of the wave port tests: attributes 14 and 18 are
// non-adjacent 1 cm x 1 cm faces on the y = 2 cm side, above the domain attributes 2 and 6,
// for the wave ports, with PEC on every other face except other_ports (TE10 above cutoff
// over the 20-30 GHz band). The caller adds the ports.
json WavePortBoxConfig(const fs::path &output, double freq_step, bool hex,
                       const std::vector<int> &other_ports)
{
  json pec_attributes = json::array();
  for (int attr = 1; attr <= 22; attr++)
  {
    if (attr != 14 && attr != 18 &&
        std::find(other_ports.begin(), other_ports.end(), attr) == other_ports.end())
    {
      pec_attributes.push_back(attr);
    }
  }
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", output}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh" /
                   (hex ? "cube_mesh_3_2_1_hex.msh" : "cube_mesh_3_2_1_tet.msh")},
      {"L0", 1.0e-2},
      {"Refinement", json::object()},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {{"PEC", {{"Attributes", pec_attributes}}}};
  setup_json["Solver"] = {
      {"Order", 2UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3},
        {"MinFreq", 20.0},
        {"MaxFreq", 30.0},
        {"FreqStep", freq_step}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};
  return setup_json;
}

// The reduced S-parameters of wave ports, from the mode pairings of the PROM, with the
// incident wave of the driven port removed and the de-embedding of the offsets as in the
// full evaluation, also at the first point of each excitation, where the reduced outputs
// are checked against the full evaluation and a mismatch reverts to full postprocessing.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced wave port postprocessing matches full evaluation",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();

  // An excited wave port on attribute 14, a passive one on attribute 18, and an excited
  // lumped port on attribute 1. The domain below the passive port is a lossy dielectric,
  // so its mode and the pairings with it are complex, and the S-parameter differs from a
  // conjugated pairing.
  json setup_json = WavePortBoxConfig(temp_dir, 5.0, false, {1});
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}}),
                                 json::object({{"Attributes", json::array({6})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 2.0},
                                               {"LossTan", 0.05}})})},
      {"Postprocessing",
       {{"Energy",
         json::array({json::object(
             {{"Index", 1}, {"Attributes", json::array({1, 2, 3, 4, 5, 6})}})})}}}};
  setup_json["Boundaries"]["WavePort"] =
      json::array({json::object({{"Index", 1},
                                 {"Attributes", json::array({14})},
                                 {"Mode", 1},
                                 {"Offset", 0.5},
                                 {"Excitation", uint(1)}}),
                   json::object({{"Index", 3},
                                 {"Attributes", json::array({18})},
                                 {"Mode", 1},
                                 {"Offset", 0.25}})});
  setup_json["Boundaries"]["LumpedPort"] =
      json::array({json::object({{"Index", 2},
                                 {"R", 50.0},
                                 {"Excitation", uint(2)},
                                 {"Attributes", json::array({1})},
                                 {"Direction", "+Y"}})});
  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperator rom_op(iodata, space_op, 4);
  DrivenPostOperatorTest post_op(iodata, space_op);

  // The PROM from HDM solutions of both excitations at the ends of the band.
  const auto &sample_f = iodata.solver.driven.sample_f;
  REQUIRE(sample_f.size() == 3);
  const int size = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector E(size);
  E.UseDevice(true);
  for (int excitation_idx : {1, 2})
  {
    for (double omega : {sample_f.front(), sample_f.back()})
    {
      rom_op.SolveHDM(excitation_idx, omega, E);
      rom_op.UpdatePROM(E, fmt::format("sample_e{:d}_{:.3e}", excitation_idx, omega));
    }
  }
  REQUIRE(rom_op.GetReducedDimension() > 0);
  post_op.ConfigureReducedPostprocessing(rom_op, iodata.solver.driven.adaptive_tol);
  REQUIRE(post_op.HasReducedPostprocessing());

  auto check_complex = [](std::complex<double> x, std::complex<double> y)
  {
    CAPTURE(x, y);
    REQUIRE(std::abs(y) > 0.0);
    CHECK(std::abs(x - y) <= 1.0e-8 * std::abs(y));
  };
  ComplexVector B(space_op.GetCurlMatrix().Height());
  B.UseDevice(true);
  // At the middle of the band, the first point of each excitation, and away from the
  // samples.
  const std::array<double, 2> omegas = {sample_f[1],
                                        0.3 * sample_f.front() + 0.7 * sample_f.back()};
  for (int step = 0; step < 2; step++)
  {
    const double omega = omegas[step];
    rom_op.SolvePROM(omega);
    for (int excitation_idx : {1, 2})
    {
      CAPTURE(omega, excitation_idx);
      const Eigen::VectorXcd y = rom_op.GetReducedSolution(excitation_idx);
      // The pairings with the mode of the lossy port are complex.
      const auto *sV = rom_op.GetWavePortModePairing(3);
      REQUIRE(sV);
      CHECK(std::abs(sV->dot(y) - sV->cwiseProduct(y).sum()) >
            1.0e-3 * std::abs(sV->cwiseProduct(y).sum()));
      rom_op.ExpandPROMSolution(excitation_idx, E);
      if (post_op.ReducedPostprocessingNeedsField(excitation_idx))
      {
        // The check of the reduced outputs against the full evaluation passes, and the
        // field is not needed at later points.
        REQUIRE(step == 0);
        post_op.MeasureAndPrintReduced(excitation_idx, step, &E, omega, y, rom_op);
        REQUIRE(post_op.HasReducedPostprocessing());
      }
      REQUIRE_FALSE(post_op.ReducedPostprocessingNeedsField(excitation_idx));
      post_op.MeasureAndPrintReduced(excitation_idx, step, nullptr, omega, y, rom_op);
      REQUIRE(post_op.HasReducedPostprocessing());
      const Measurement reduced = post_op.GetMeasurement();

      space_op.GetCurlMatrix().Mult(E.Real(), B.Real());
      space_op.GetCurlMatrix().Mult(E.Imag(), B.Imag());
      B *= -1.0 / (1i * omega);
      post_op.MeasureAndPrintAll(excitation_idx, step, E, B, omega);
      const Measurement &full = post_op.GetMeasurement();

      CAPTURE(reduced.domain_E_field_energy_all.value(),
              full.domain_E_field_energy_all.value(), reduced.domain_H_field_energy_all,
              full.domain_H_field_energy_all);
      REQUIRE(full.domain_E_field_energy_all.value() > 0.0);
      REQUIRE(full.domain_H_field_energy_all > 0.0);
      CHECK_THAT(reduced.domain_E_field_energy_all.value(),
                 WithinRel(full.domain_E_field_energy_all.value(), 1.0e-8));
      CHECK_THAT(reduced.domain_H_field_energy_all,
                 WithinRel(full.domain_H_field_energy_all, 1.0e-8));
      REQUIRE(full.wave_port_vi.size() == 2);
      REQUIRE(reduced.wave_port_vi.size() == full.wave_port_vi.size());
      for (const auto &[idx, full_port] : full.wave_port_vi)
      {
        CAPTURE(idx);
        check_complex(reduced.wave_port_vi.at(idx).S, full_port.S);
      }
      REQUIRE(full.lumped_port_vi.size() == 1);
      REQUIRE(reduced.lumped_port_vi.size() == full.lumped_port_vi.size());
      for (const auto &[idx, full_port] : full.lumped_port_vi)
      {
        CAPTURE(idx);
        const auto &reduced_port = reduced.lumped_port_vi.at(idx);
        check_complex(reduced_port.V, full_port.V);
        check_complex(reduced_port.I, full_port.I);
        check_complex(reduced_port.S, full_port.S);
      }
    }
  }
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator reduced postprocessing rejects Floquet wave vector",
                 "[postoperator][romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path = fs::path(PALACE_TEST_DATA_DIR) /
                         "regression/input/floquet_lumped/mesh/floquet_lumped.msh";

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path}, {"L0", 1.0e-2}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1, 2})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"Periodic",
       {{"FloquetWaveVector", json::array({0.1, 0.0, 0.0})},
        {"BoundaryPairs",
         json::array({json::object({{"DonorAttributes", json::array({1})},
                                    {"ReceiverAttributes", json::array({2})}}),
                      json::object({{"DonorAttributes", json::array({3})},
                                    {"ReceiverAttributes", json::array({4})}})})}}},
      {"PEC", {{"Attributes", json::array({8})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"Attributes", json::array({7})},
                                                {"Excitation", uint(1)},
                                                {"R", 377.0},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 5.0}, {"MaxFreq", 10.0}, {"FreqStep", 5.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-8}}}};

  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  REQUIRE(space_op.GetMaterialOp().HasWaveVector());
  RomOperator rom_op(iodata, space_op, 2);
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  sample.Real().Randomize(42);
  sample.Imag() = 0.0;
  rom_op.UpdatePROM(sample, "wave_vector_sample");
  REQUIRE(rom_op.GetReducedDimension() > 0);

  DrivenPostOperatorTest post_op(iodata, space_op);
  post_op.ConfigureReducedPostprocessing(rom_op, iodata.solver.driven.adaptive_tol);
  CHECK_FALSE(post_op.HasReducedPostprocessing());
}

// Check the PROM system assembled online by RomOperator::SolvePROM, from projected
// frequency-dependent terms and port-space wave port pairings, against the Galerkin
// projection of the HDM system matrix and excitation at the same frequency, built from the
// basis of two random samples.
void CheckPROMAgainstHDM(IoData &iodata, const std::vector<int> &excitations)
{
  MPI_Comm comm = Mpi::World();
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  const auto &dbc_tdof_list = space_op.GetNDDbcTDofLists().back();
  const int size = space_op.GetNDSpace().GetTrueVSize();
  ComplexVector sample(size);
  sample.UseDevice(true);
  for (int seed : {42, 44})
  {
    linalg::SetRandom(comm, sample.Real(), seed);
    linalg::SetRandom(comm, sample.Imag(), seed + 1);
    linalg::SetSubVector(sample, dbc_tdof_list, 0.0);
    rom_op.UpdatePROM(sample, fmt::format("sample_{:d}", seed));
  }
  const auto &V = rom_op.GetBasis();
  const auto n = static_cast<Eigen::Index>(V.size());
  REQUIRE(n > 0);

  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  ComplexVector v(size), w(size), rhs(size), b(size);
  for (auto *x : {&v, &w, &rhs, &b})
  {
    x->UseDevice(true);
  }
  auto project = [&](const ComplexVector &x, Eigen::Index j, Eigen::MatrixXcd &X)
  {
    for (Eigen::Index i = 0; i < n; i++)
    {
      X(i, j) = {V[i] * x.Real(), V[i] * x.Imag()};
    }
  };
  for (double omega :
       {iodata.solver.driven.sample_f.front(), iodata.solver.driven.sample_f.back()})
  {
    for (int excitation_idx : excitations)
    {
      CAPTURE(omega, excitation_idx);
      rom_op.SolvePROM(omega);

      // The port modes are cached at this frequency, so the comparison is not affected by
      // the sign ambiguity of a re-solved eigenvector.
      auto A2 = space_op.GetExtraSystemOperator(omega, Operator::DIAG_ZERO);
      auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * omega,
                                        std::complex<double>(-omega * omega, 0.0), K.get(),
                                        C.get(), M.get(), A2.get());
      b = 0.0;
      if (space_op.GetExcitationVector1(excitation_idx, rhs))
      {
        b.Add(1i * omega, rhs);
      }
      if (space_op.GetExcitationVector2(excitation_idx, omega, rhs))
      {
        b.Add(1.0, rhs);
      }
      Eigen::MatrixXcd Ar(n, n), br(n, 1);
      for (Eigen::Index j = 0; j < n; j++)
      {
        v.Real() = V[j];
        v.Imag() = 0.0;
        A->Mult(v, w);
        project(w, j, Ar);
      }
      project(b, 0, br);
      Mpi::GlobalSum(static_cast<int>(Ar.size()), Ar.data(), comm);
      Mpi::GlobalSum(static_cast<int>(br.size()), br.data(), comm);

      const double error_A = (rom_op.GetAr() - Ar).cwiseAbs().maxCoeff();
      const Eigen::VectorXcd y = Ar.partialPivLu().solve(br.col(0));
      const double error_y = (rom_op.GetReducedSolution(excitation_idx) - y).norm();
      CAPTURE(error_A, Ar.cwiseAbs().maxCoeff(), error_y, y.norm());
      REQUIRE(br.norm() > 0.0);
      CHECK(error_A <= 1.0e-10 * Ar.cwiseAbs().maxCoeff());
      CHECK(error_y <= 1.0e-8 * y.norm());

      // The wave port S-parameters from the pairings, sᵀV y / (sqrt(2) R) relative to the
      // incident wave of sqrt(2) times the mode, are those of the expanded solution.
      rom_op.ExpandPROMSolution(excitation_idx, v);
      GridFunction E(space_op.GetNDSpace(), true);
      E.Real().SetFromTrueDofs(v.Real());
      E.Imag().SetFromTrueDofs(v.Imag());
      for (const auto &[idx, data] : space_op.GetWavePortOp())
      {
        const auto *sV = rom_op.GetWavePortModePairing(idx);
        REQUIRE(sV);
        const auto S = sV->cwiseProduct(rom_op.GetReducedSolution(excitation_idx)).sum() /
                       (std::sqrt(2.0) * data.modal_reaction);
        const auto S_ref = data.GetSParameter(E);
        CAPTURE(idx, S, S_ref);
        CHECK(std::abs(S - S_ref) <= 1.0e-12 * std::max(std::abs(S_ref), 1.0));
      }
    }
  }

  // The one-time check of the pairings against the assembled mode vectors ran and agrees.
  if (space_op.GetWavePortOp().Size() > 0)
  {
    CHECK(rom_op.wp_pairing_error >= 0.0);
    CHECK(rom_op.wp_pairing_error <= 1.0e-8);
  }
}

// Wave ports: the excitation of an inactive port, sharing the excitation of the active
// port or on its own, and the pairings on a nonconformal hexahedral mesh, with hanging
// dofs on the boundary of the port.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator wave ports match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const unsigned int inactive_excitation = GENERATE(1u, 2u);
  const bool hex = GENERATE(false, true);
  CAPTURE(inactive_excitation, hex);

  // An active and an inactive (excited) wave port on attributes 14 and 18.
  json setup_json = WavePortBoxConfig(temp_dir, 10.0, hex, {});
  if (hex)
  {
    // Refine the element below port 1.
    setup_json["Model"]["Refinement"] = {
        {"Boxes", json::array({json::object({{"Levels", 1},
                                             {"BoundingBoxMin", {0.0, 1.9, 0.0}},
                                             {"BoundingBoxMax", {0.1, 2.0, 0.1}}})})}};
  }
  setup_json["Boundaries"]["WavePort"] =
      json::array({json::object({{"Index", 1},
                                 {"Attributes", json::array({14})},
                                 {"Mode", 1},
                                 {"Offset", 0.0},
                                 {"Excitation", uint(1)}}),
                   json::object({{"Index", 2},
                                 {"Attributes", json::array({18})},
                                 {"Mode", 1},
                                 {"Offset", 0.0},
                                 {"Active", false},
                                 {"Excitation", inactive_excitation}})});
  IoData iodata(setup_json, false);
  std::vector<int> excitations = {1};
  if (inactive_excitation != 1)
  {
    excitations.push_back(static_cast<int>(inactive_excitation));
  }
  CheckPROMAgainstHDM(iodata, excitations);
}

// Second-order farfield, surface conductivity, and rational surface impedance boundaries,
// with a lumped port excitation, on the 3x2x1 tetrahedral cube mesh.
json BoundaryTermsSetup(const fs::path &temp_dir)
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"},
      {"L0", 1.0e-2},
      {"Refinement", json::object({})},
      {"CrackInternalBoundaryElements", false}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({14})},
                                                {"Direction", "+Z"}})})},
      {"Absorbing", {{"Attributes", json::array({18, 19, 20})}, {"Order", 2}}},
      {"Conductivity", json::array({json::object({{"Attributes", json::array({13, 15})},
                                                  {"Conductivity", 1.0e3}}),
                                    json::object({{"Attributes", json::array({17})},
                                                  {"Conductivity", 5.8e7},
                                                  {"Thickness", 1.0e-6}})})},
      {"RationalImpedance",
       json::array(
           {json::object({{"Attributes", json::array({21, 22})},
                          {"Numerator", json::array({2.0e-6, 0.0})},
                          {"Denominator", json::array({1.0e-19, 6.0e-9, 377.0})}})})}};
  setup_json["Solver"] = {
      {"Order", 2UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3}, {"MinFreq", 2.0}, {"MaxFreq", 30.0}, {"FreqStep", 28.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};
  return setup_json;
}

TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator boundary terms match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  json setup_json = BoundaryTermsSetup(temp_dir);
  IoData iodata(setup_json, false);
  CheckPROMAgainstHDM(iodata, {1});
}

// The configuration of a regression test, with the paths of its mesh and output, quiet,
// and at order 1.
json RegressionSetup(const fs::path &temp_dir, std::string_view dir, std::string_view file)
{
  const auto input_dir = fs::path(PALACE_TEST_DATA_DIR) / "regression" / "input" / dir;
  std::ifstream f(input_dir / file);
  REQUIRE(f.good());
  json setup_json = json::parse(f, nullptr, true, /*ignore_comments=*/true);
  setup_json["Model"]["Mesh"] =
      (input_dir / setup_json["Model"]["Mesh"].get<std::string>()).string();
  setup_json["Problem"]["Output"] = temp_dir;
  setup_json["Problem"]["Verbose"] = 0;
  setup_json["Solver"]["Order"] = 1;
  return setup_json;
}

// Floquet ports, with a frequency-dependent (fixed angle of incidence) or fixed Floquet
// wave vector, together with second-order farfield boundaries or a wave port.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator Floquet ports match the HDM system",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const auto [dir, file, excitations] = GENERATE(
      std::make_tuple("floquet_absorbing", "floquet_absorbing_adaptive.json",
                      std::vector<int>{1}),
      std::make_tuple("floquet_absorbing", "floquet_absorbing_scatterer_adaptive.json",
                      std::vector<int>{1}),
      std::make_tuple("floquet_wave", "floquet_wave.json", std::vector<int>{1, 2}));
  CAPTURE(file);
  json setup_json = RegressionSetup(temp_dir, dir, file);
  if (std::string(file) == "floquet_absorbing_scatterer_adaptive.json")
  {
    // Both polarizations of the incident field.
    setup_json["Boundaries"]["FloquetPort"][0]["IncidentPolarization"] = "RHC";
  }
  IoData iodata(setup_json, false);
  CheckPROMAgainstHDM(iodata, excitations);
}

// The frequency-dependent boundary mass terms of finite conductivity and rational impedance
// boundaries enter the synthesized circuit of an adaptive sweep as fits of the functions of
// frequency which scale their projected boundary masses. Synthesis folds in each term, for
// each type alone and for both, as auxiliary states labeled by the term (the fit of each
// term here has poles): the exact split of the series RLC termination gives the two poles
// of its Robin coefficient, the zeros of Zs. Eliminating the auxiliary states of the
// synthesized circuit must give the reduced model over the band, to the fit tolerance
// (AdaptiveTol) relative to the largest projected frequency-dependent operator.
// "RomOperator boundary terms match the HDM system" tests that the reduced model applies
// these terms as the full operator does.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator synthesis of frequency-dependent boundary terms",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path = fs::path(PALACE_TEST_DATA_DIR) /
                         "regression/input/rational_impedance/mesh/parallel_plate.msh";
  const auto types = GENERATE(as<std::string>{}, "Conductivity", "RationalImpedance",
                              "Conductivity+RationalImpedance");
  CAPTURE(types);
  const bool conductivity = types.find("Conductivity") != std::string::npos;
  const bool rational = types.find("RationalImpedance") != std::string::npos;

  // Parallel plates (attribute 2) between PMC side walls (3), a lumped port (4) and a
  // termination (5): finite conductivity plates or PEC, and a series RLC sheet termination
  // (notch at 10 GHz, inside the band) or PEC.
  json pec = json::array();
  json boundaries = {
      {"PMC", {{"Attributes", json::array({3})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"Attributes", json::array({4})},
                                                {"R", 50.0},
                                                {"Direction", "+Y"},
                                                {"Excitation", uint(1)}})})}};
  if (conductivity)
  {
    boundaries["Conductivity"] = json::array(
        {json::object({{"Attributes", json::array({2})}, {"Conductivity", 1.0e5}})});
  }
  else
  {
    pec.push_back(2);
  }
  if (rational)
  {
    constexpr double L = 2.0e-9, C = 1.2665e-13, R = 50.0;
    boundaries["RationalImpedance"] =
        json::array({json::object({{"Attributes", json::array({5})},
                                   {"Numerator", json::array({L * C, R * C, 1.0})},
                                   {"Denominator", json::array({C, 0.0})}})});
  }
  else
  {
    pec.push_back(5);
  }
  boundaries["PEC"] = {{"Attributes", pec}};
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path}, {"L0", 1.0e-3}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = boundaries;
  setup_json["Solver"] = {
      {"Order", 1UL},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-6},
        {"AdaptiveCircuitSynthesis", true},
        {"MinFreq", 6.0},
        {"MaxFreq", 14.0},
        {"FreqStep", 4.0}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  rom_op.AddLumpedPortModesForSynthesis();
  ComplexVector sample(space_op.GetNDSpace().GetTrueVSize());
  sample.UseDevice(true);
  for (int seed : {42, 44})
  {
    linalg::SetRandom(comm, sample.Real(), seed);
    linalg::SetRandom(comm, sample.Imag(), seed + 1);
    linalg::SetSubVector(sample.Real(), space_op.GetNDDbcTDofLists().back(), 0.0);
    linalg::SetSubVector(sample.Imag(), space_op.GetNDDbcTDofLists().back(), 0.0);
    rom_op.UpdatePROM(sample, fmt::format("terms_sample_{:d}", seed));
  }
  REQUIRE(rom_op.GetReducedDimension() > 1);

  const auto matrices = rom_op.CalculateNormalizedPROMMatrices(iodata.units);
  auto num_poles = [&](const std::string &label)
  {
    // Auxiliary states <label>_p<k>d<j>, for pole k and coupling direction j.
    std::set<std::string> poles;
    for (const auto &aux : matrices.aux_labels)
    {
      if (aux.starts_with(label + "_p"))
      {
        poles.insert(aux.substr(0, aux.rfind('d')));
      }
    }
    return poles.size();
  };
  CHECK((num_poles("surfsigma_0") > 0) == conductivity);
  CHECK(num_poles("rationalz_0") == (rational ? 2 : 0));

  // The augmented pencil Kₐ + iωCₐ − ω²Mₐ of the circuit: the synthesized L⁻¹, R⁻¹ and C
  // without their units and the scaling of the port row (the first basis vector).
  const long n = static_cast<long>(rom_op.GetReducedDimension());
  const long n_aug = matrices.L_inv->rows();
  REQUIRE(n_aug > n);
  Eigen::VectorXcd v_inv = Eigen::VectorXcd::Ones(n_aug);
  v_inv(0) = 1.0 / rom_op.GetOrthR()(0, 0);
  const auto &units = iodata.units;
  auto unscale = [&](const std::unique_ptr<Eigen::MatrixXcd> &A, double scale)
  {
    return A ? Eigen::MatrixXcd(scale * v_inv.asDiagonal() * (*A) * v_inv.asDiagonal())
             : Eigen::MatrixXcd(Eigen::MatrixXcd::Zero(n_aug, n_aug));
  };
  const auto Ka =
      unscale(matrices.L_inv, units.GetScaleFactor<Units::ValueType::INDUCTANCE>());
  const auto Ca =
      unscale(matrices.R_inv, units.GetScaleFactor<Units::ValueType::IMPEDANCE>());
  const auto Ma =
      unscale(matrices.C, 1.0 / units.GetScaleFactor<Units::ValueType::CAPACITANCE>());

  // The reduced model: Kr + iωCr − ω²Mr and the projection of the frequency-dependent
  // operator onto the basis.
  const auto &V = rom_op.GetBasis();
  REQUIRE(static_cast<long>(V.size()) == n);
  auto ProjectedA2 = [&](double omega)
  {
    Eigen::MatrixXcd A2r = Eigen::MatrixXcd::Zero(n, n);
    auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(omega, Operator::DIAG_ZERO);
    REQUIRE(A2);
    ComplexVector x(V[0].Size()), y(V[0].Size());
    x.UseDevice(true);
    y.UseDevice(true);
    for (long j = 0; j < n; j++)
    {
      x.Real() = V[j];
      x.Imag() = 0.0;
      A2->Mult(x, y);
      for (long i = 0; i < n; i++)
      {
        A2r(i, j) = {linalg::Dot(comm, V[i], y.Real()), linalg::Dot(comm, V[i], y.Imag())};
      }
    }
    return A2r;
  };
  const auto &Kr = rom_op.GetKr(), &Cr = rom_op.GetCr(), &Mr = rom_op.GetMr();
  REQUIRE(Cr.rows() == n);
  const double omega_min = iodata.solver.driven.sample_f.front();
  const double omega_max = iodata.solver.driven.sample_f.back();
  double err = 0.0, scale = 0.0;
  for (int k = 0; k <= 32; k++)
  {
    const double omega = omega_min + k * (omega_max - omega_min) / 32.0;
    const std::complex<double> iw(0.0, omega);
    const Eigen::MatrixXcd A = Ka + iw * Ca - omega * omega * Ma;
    const Eigen::MatrixXcd schur =
        A.topLeftCorner(n, n) -
        A.topRightCorner(n, n_aug - n) * A.bottomRightCorner(n_aug - n, n_aug - n)
                                             .fullPivLu()
                                             .solve(A.bottomLeftCorner(n_aug - n, n));
    const Eigen::MatrixXcd A2r = ProjectedA2(omega);
    CAPTURE(omega);
    REQUIRE(schur.allFinite());
    REQUIRE(A2r.allFinite());
    const Eigen::MatrixXcd diff = schur - (Kr + iw * Cr - omega * omega * Mr) - A2r;
    REQUIRE(diff.allFinite());
    err = std::max(err, diff.cwiseAbs().maxCoeff());
    scale = std::max(scale, A2r.cwiseAbs().maxCoeff());
  }
  CAPTURE(err, scale);
  REQUIRE(std::isfinite(err));
  REQUIRE(std::isfinite(scale));
  REQUIRE(scale > 0.0);
  CHECK(err <= iodata.solver.driven.adaptive_tol * scale);
}

// The series RLC sheet termination of the previous test, with its series resonance above
// the band (at 100 or 300 GHz, for a band of 6 to 14 GHz), on the lumped port mode and one
// random vector, with the projected mass of the termination replaced by a rotated diag(1,
// 1e-7): its weak direction is below the rank cutoff of the auxiliary states. The exact
// split of the Robin coefficient has a constant polynomial part and two poles, which below
// the resonance cancel to a coefficient smaller by (f_res/f)². The polynomial part must be
// realized on the directions of the poles, as for the other fits: on the full mass, the
// weak direction carries the polynomial part alone, and the circuit missed the term by
// 2.8e-5 (100 GHz) and 2.5e-4 (300 GHz) of its value while the reported fit residual is
// below 1e-12.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator synthesis of a rational impedance on a weak direction",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const double f_res = GENERATE(100.0, 300.0);  // GHz
  CAPTURE(f_res);
  constexpr double L = 2.0e-9, R = 50.0;
  const double C = 1.0 / (std::pow(2.0 * std::numbers::pi * f_res * 1.0e9, 2) * L);
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {
      {"Mesh", fs::path(PALACE_TEST_DATA_DIR) / "regression/input/rational_impedance/mesh/"
                                                "parallel_plate.msh"},
      {"L0", 1.0e-3}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"PMC", {{"Attributes", json::array({3})}}},
      {"PEC", {{"Attributes", json::array({2})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"Attributes", json::array({4})},
                                                {"R", 50.0},
                                                {"Direction", "+Y"},
                                                {"Excitation", uint(1)}})})},
      {"RationalImpedance",
       json::array({json::object({{"Attributes", json::array({5})},
                                  {"Numerator", json::array({L * C, R * C, 1.0})},
                                  {"Denominator", json::array({C, 0.0})}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Device", "CPU"},
                          {"Driven",
                           {{"AdaptiveTol", 1.0e-6},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 6.0},
                            {"MaxFreq", 14.0},
                            {"FreqStep", 4.0}}}};
  IoData iodata(setup_json, false);
  auto mesh = LoadScaleParMesh2(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  RomOperatorTest rom_op(iodata, space_op, 4);
  rom_op.AddLumpedPortModesForSynthesis();
  ComplexVector v(space_op.GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 0.0;
  linalg::SetRandom(comm, v.Real(), 42);
  linalg::SetSubVector(v.Real(), space_op.GetNDDbcTDofLists().back(), 0.0);
  rom_op.UpdatePROM(v, "random");
  const long n = static_cast<long>(rom_op.GetReducedDimension());
  REQUIRE(n == 2);
  const auto [w_min, w_max] = rom_op.GetSweepBand();
  REQUIRE(w_max > w_min);

  const auto &surface_terms = space_op.GetSurfaceTerms();
  REQUIRE(surface_terms.GetNumExtraSystemBdrMassTerms() == 1);
  REQUIRE(surface_terms.GetExtraSystemBdrMassTermLabel(0) == "rationalz_0");
  const auto split = surface_terms.GetExtraSystemBdrMassTermSplit(0);
  REQUIRE(split);
  REQUIRE(split->degree <= 2);
  auto f = [&surface_terms](double w)
  {
    return surface_terms.GetExtraSystemBdrMassTermScale(0, w) /
           std::complex<double>(0.0, 1.0);
  };

  // i·VᵀMV by the imaginary-slot convention, replaced by i·σ·Q diag(1, 1e-7) Qᵀ. The scale
  // σ makes the term larger than the rest of the pencil over the band, so that the
  // round-off of the rest is negligible against the error of the term (the rank cutoff is
  // relative).
  const auto &Kr = rom_op.GetKr(), &Cr = rom_op.GetCr(), &Mr = rom_op.GetMr();
  REQUIRE(Cr.rows() == n);
  const double pencil = Kr.norm() + w_max * Cr.norm() + w_max * w_max * Mr.norm();
  const double sigma = pencil / std::min(std::abs(f(w_min)), std::abs(f(w_max)));
  const double c = std::cos(0.3), s = std::sin(0.3);
  Eigen::Matrix2d Q, D = Eigen::Vector2d(1.0, 1.0e-7).asDiagonal();
  Q << c, -s, s, c;
  auto &M_r =
      rom_op.GetBoundaryTermProjection(RomOperatorTest::BoundaryTermKind::SurfaceTerm, 0);
  REQUIRE(M_r.rows() == n);
  M_r = std::complex<double>(0.0, sigma) *
        (Q * D * Q.transpose()).cast<std::complex<double>>();
  const auto fit = rom_op.FitExactSplitDispersion("check", M_r, *split);
  REQUIRE(fit.regime == RomOperatorTest::WavePortRegime::Augmented);

  // Two poles, each on the one kept direction.
  const auto matrices = rom_op.CalculateNormalizedPROMMatrices(iodata.units);
  CHECK(std::ranges::count_if(matrices.aux_labels, [](const std::string &l)
                              { return l.starts_with("rationalz_0_p"); }) == 2);
  CHECK(std::ranges::count_if(
            matrices.aux_labels, [](const std::string &l)
            { return l.starts_with("rationalz_0_p") && l.ends_with("d0"); }) == 2);

  // The augmented pencil of the circuit, without its units and the scaling of the port row
  // (as in the previous test), against the reduced model with the replaced term.
  const long n_aug = matrices.L_inv->rows();
  REQUIRE(n_aug > n);
  Eigen::VectorXcd v_inv = Eigen::VectorXcd::Ones(n_aug);
  v_inv(0) = 1.0 / rom_op.GetOrthR()(0, 0);
  const auto &units = iodata.units;
  auto unscale = [&](const std::unique_ptr<Eigen::MatrixXcd> &A, double scale)
  {
    return A ? Eigen::MatrixXcd(scale * v_inv.asDiagonal() * (*A) * v_inv.asDiagonal())
             : Eigen::MatrixXcd(Eigen::MatrixXcd::Zero(n_aug, n_aug));
  };
  const auto Ka =
      unscale(matrices.L_inv, units.GetScaleFactor<Units::ValueType::INDUCTANCE>());
  const auto Ca =
      unscale(matrices.R_inv, units.GetScaleFactor<Units::ValueType::IMPEDANCE>());
  const auto Ma =
      unscale(matrices.C, 1.0 / units.GetScaleFactor<Units::ValueType::CAPACITANCE>());
  double max_rel = 0.0;
  for (int k = 0; k <= 80; k++)
  {
    const double w = w_min + (w_max - w_min) * k / 80.0;
    const std::complex<double> iw(0.0, w);
    const Eigen::MatrixXcd A = Ka + iw * Ca - w * w * Ma;
    const Eigen::MatrixXcd schur =
        A.topLeftCorner(n, n) -
        A.topRightCorner(n, n_aug - n) * A.bottomRightCorner(n_aug - n, n_aug - n)
                                             .fullPivLu()
                                             .solve(A.bottomLeftCorner(n_aug - n, n));
    const Eigen::MatrixXcd term = f(w) * M_r;
    const Eigen::MatrixXcd diff = schur - (Kr + iw * Cr - w * w * Mr) - term;
    REQUIRE(diff.allFinite());
    max_rel = std::max(max_rel, diff.cwiseAbs().maxCoeff() / term.cwiseAbs().maxCoeff());
  }
  CAPTURE(fit.rel_err_augmented, max_rel);
  // The weak direction is not realized, an error of 1e-7 of the term, as the rank cutoff
  // allows.
  CHECK(max_rel <= 2.0 * fit.rel_err_augmented + 2.0e-7);
}

// The models of the TruncatePROM, RefineHDM and AddSamplePROM tests, with every kind of
// term of the PROM: lumped ports (on the 3x2x1 µm cube of the other lumped port tests),
// wave ports, circuit synthesis (a basis orthonormal in a weighted inner product) with
// lumped or wave ports, the boundary terms of a second-order farfield, surface
// conductivity and rational impedance, and Floquet ports (with a wave port, or with a
// frequency-dependent Floquet wave vector and a second-order farfield). The wave port,
// lumped port synthesis and Floquet port with wave port models have two excitations, which
// share the basis.
json PROMModelSetup(const fs::path &temp_dir, const std::string &model,
                    const std::string &linear_type)
{
  auto wave_port = [](int index, int attr, unsigned int excitation)
  {
    return json::object({{"Index", index},
                         {"Attributes", json::array({attr})},
                         {"Mode", 1},
                         {"Offset", 0.0},
                         {"Excitation", excitation}});
  };
  json setup_json;
  if (model == "lumped port")
  {
    setup_json = WeightedExcitationSetup(temp_dir, 1.0, 0.0);
    setup_json["Boundaries"]["LumpedPort"] = json::array({json::object(
        {{"Index", 1},
         {"R", 50.0},
         {"Excitation", uint(1)},
         {"Elements", json::array({json::object(
                          {{"Attributes", json::array({14})}, {"Direction", "+Z"}})})}})});
  }
  else if (model == "wave ports" || model == "wave port synthesis")
  {
    setup_json = WavePortBoxConfig(temp_dir, 5.0, false, {});
    setup_json["Boundaries"]["WavePort"] = json::array(
        {wave_port(1, 14, 1), wave_port(2, 18, (model == "wave ports") ? 2 : 0)});
    setup_json["Solver"]["Driven"]["AdaptiveCircuitSynthesis"] = (model != "wave ports");
  }
  else if (model == "lumped port synthesis")
  {
    // The 3x2x1 cm cube, at most a wavelength across: neither are the HDM solutions those
    // of the ports nor have they a constant phase.
    setup_json = ReducedPostprocessingSetup(temp_dir);
    setup_json["Model"]["L0"] = 1.0e-2;
  }
  else if (model == "boundary terms")
  {
    setup_json = BoundaryTermsSetup(temp_dir);
  }
  else if (model == "Floquet and wave ports")
  {
    setup_json = RegressionSetup(temp_dir, "floquet_wave", "floquet_wave.json");
  }
  else
  {
    REQUIRE(model == "Floquet port and farfield");
    setup_json =
        RegressionSetup(temp_dir, "floquet_absorbing", "floquet_absorbing_adaptive.json");
  }

  // HDM solves to the tolerance set with SetLinearTol.
  setup_json["Solver"]["Linear"]["Type"] = linear_type;
  setup_json["Solver"]["Linear"]["MaxIts"] = 1000;
  return setup_json;
}

// ‖Vᵀ W (u - V c)‖ / ‖u‖_W for the basis V of rom_op, which is orthonormal in the inner
// product with the weight W (that of circuit synthesis, or the identity): zero if V c is
// the orthogonal projection of u onto the basis.
double ProjectionResidual(MPI_Comm comm, const RomOperatorTest &rom_op,
                          const ComplexVector &u, const Eigen::VectorXcd &c)
{
  const auto &V = rom_op.GetVectors();
  REQUIRE(c.size() == static_cast<Eigen::Index>(V.size()));
  ComplexVector d(u.Size()), work(u.Size());
  d.UseDevice(true);
  work.UseDevice(true);
  linalg::LinearCombination(V, {c.data(), static_cast<std::size_t>(c.size())}, d);
  d -= u;
  const auto &W = rom_op.GetWeightOp();
  auto dot = [&](const Vector &x, const Vector &y)
  { return W ? W->InnerProduct(comm, x, y, work.Real()) : linalg::Dot(comm, x, y); };
  double residual = 0.0;
  for (const auto &v : V)
  {
    residual = std::max(residual,
                        std::abs(std::complex<double>(dot(v, d.Real()), dot(v, d.Imag()))));
  }
  return residual / std::sqrt(dot(u.Real(), u.Real()) + dot(u.Imag(), u.Imag()));
}

// The model of PROMModelSetup with its operators, a RomOperator and a reference RomOperator
// for the same model, both with the port modes of circuit synthesis, and helpers. The HDM
// is solved with an iterative preconditioner by default (available in every build), with
// which the solution of a loose tolerance is about as accurate as asked for: a sparse
// direct one solves some of the models to round-off in one iteration.
struct PROMModel
{
  MPI_Comm comm;
  IoData iodata;
  std::vector<std::unique_ptr<Mesh>> mesh;
  std::unique_ptr<SpaceOperator> space_op;
  std::unique_ptr<RomOperatorTest> op, ref;
  std::vector<int> excitations;
  double omega, omega_other;

  PROMModel(const fs::path &temp_dir, const std::string &model,
            const std::string &linear_type = "AMS")
    : comm(Mpi::World()), iodata(PROMModelSetup(temp_dir, model, linear_type), false)
  {
    mesh = LoadScaleParMesh2(iodata, comm);
    space_op = std::make_unique<SpaceOperator>(iodata, mesh);
    const auto &sample_f = iodata.solver.driven.sample_f;
    omega = 0.5 * (sample_f.front() + sample_f.back());
    omega_other = sample_f.front();
    for (auto *p : {&op, &ref})
    {
      *p = std::make_unique<RomOperatorTest>(iodata, *space_op, 4);
      if (iodata.solver.driven.adaptive_circuit_synthesis)
      {
        (*p)->AddLumpedPortModesForSynthesis();
        if (space_op->GetWavePortOp().Size() > 0)
        {
          (*p)->AddWavePortModesForSynthesis(omega);
        }
      }
    }
    for (const auto &[idx, data] : space_op->GetPortExcitations())
    {
      excitations.push_back(idx);
    }
  }

  ComplexVector NewVector() const
  {
    ComplexVector v(space_op->GetNDSpace().GetTrueVSize());
    v.UseDevice(true);
    v = 0.0;
    return v;
  }

  // A random vector, zero on the Dirichlet boundary.
  ComplexVector Random(int seed) const
  {
    ComplexVector v = NewVector();
    linalg::SetRandom(comm, v.Real(), seed);
    linalg::SetRandom(comm, v.Imag(), seed + 1);
    linalg::SetSubVector(v, space_op->GetNDDbcTDofLists().back(), 0.0);
    return v;
  }

  double Norm(const ComplexVector &x) const { return linalg::Norml2(comm, x); }

  // ‖x - y‖ / ‖y‖.
  double RelDiff(const ComplexVector &x, const ComplexVector &y) const
  {
    ComplexVector d(x);
    d.UseDevice(true);
    linalg::AXPY(-1.0, y, d);
    return Norm(d) / Norm(y);
  }

  // The PROM solution of every excitation at ω.
  std::vector<ComplexVector> Solve(RomOperatorTest &rom_op, double w) const
  {
    rom_op.SolvePROM(w);
    std::vector<ComplexVector> x;
    for (int e : excitations)
    {
      rom_op.ExpandPROMSolution(e, x.emplace_back(NewVector()));
    }
    return x;
  }

  // The largest relative difference of the PROM solutions of op and ref, at ω and at
  // another frequency, over the excitations.
  double MaxSolutionDiff() const
  {
    double diff = 0.0;
    for (double w : {omega, omega_other})
    {
      const auto x = Solve(*op, w), x_ref = Solve(*ref, w);
      for (std::size_t k = 0; k < x.size(); k++)
      {
        diff = std::max(diff, RelDiff(x[k], x_ref[k]));
      }
    }
    return diff;
  }
};

const std::vector<std::string> PROM_MODELS = {"lumped port",
                                              "wave ports",
                                              "lumped port synthesis",
                                              "wave port synthesis",
                                              "boundary terms",
                                              "Floquet and wave ports",
                                              "Floquet port and farfield"};

// TruncatePROM(n), with which AddSamplePROM removes a sample to add it again refined, gives
// back the PROM of the first n basis vectors exactly, also when n splits a sample (the real
// part of a sample is added before its imaginary part), and the basis can grow again from
// there as if the removed vectors had never been added. The basis restricted to the wave
// ports, which SolvePROM extends, is truncated too: the next update, which adds as many
// vectors as were removed, would otherwise reuse the stale restricted vectors.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-TruncatePROM",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const std::string model = GENERATE(from_range(PROM_MODELS));
  CAPTURE(model);
  PROMModel m(temp_dir, model);
  auto &op = *m.op;
  auto &ref = *m.ref;
  RomOperatorTest ref_re(m.iodata, *m.space_op, 4);
  if (m.iodata.solver.driven.adaptive_circuit_synthesis)
  {
    ref_re.AddLumpedPortModesForSynthesis();
    if (m.space_op->GetWavePortOp().Size() > 0)
    {
      ref_re.AddWavePortModesForSynthesis(m.omega);
    }
  }

  // The terms of the model are projected.
  const auto &space_op = *m.space_op;
  const std::map<std::string, std::size_t> num_bdr_terms = {
      {"lumped port", 0},
      {"wave ports", 2},
      {"lumped port synthesis", 0},
      {"wave port synthesis", 2},
      {"boundary terms", 4},
      {"Floquet and wave ports", 2},
      {"Floquet port and farfield", 2}};
  REQUIRE(op.NumBoundaryTerms() == num_bdr_terms.at(model));
  REQUIRE((op.NumFloquetOrders() > 0) == model.starts_with("Floquet"));
  REQUIRE(op.GetWeightOp().has_value() == (model.find("synthesis") != std::string::npos));

  const ComplexVector s1 = m.Random(71), s2 = m.Random(73), s3 = m.Random(75);
  for (auto *rom_op : {&op, &ref, &ref_re})
  {
    rom_op->UpdatePROM(s1, "s1");
  }
  const std::size_t n = op.GetReducedDimension();
  const auto x_n = m.Solve(op, m.omega);  // Extends the basis restricted to the wave ports
  const auto state_n = op.GetProjectedState();
  CHECK(state_n.V_wp_dim == ((space_op.GetWavePortOp().Size() > 0) ? n : 0));
  CHECK(state_n.RHS1r.cols() == static_cast<Eigen::Index>(m.excitations.size()));
  CHECK(state_n.floquet_comp_V.rows() == static_cast<Eigen::Index>(n));
  CHECK(state_n.gram_L.rows() == (op.GetWeightOp() ? static_cast<Eigen::Index>(n) : 0));

  op.UpdatePROM(s2, "s2");
  REQUIRE(op.GetReducedDimension() == n + 2);
  m.Solve(op, m.omega);
  REQUIRE(op.GetProjectedState().V_wp_dim ==
          ((space_op.GetWavePortOp().Size() > 0) ? n + 2 : 0));

  // Truncated to the basis with the real part of s2: as projecting the real part alone.
  ComplexVector s2_re(s2);
  s2_re.UseDevice(true);
  s2_re.Imag() = 0.0;
  ref_re.UpdatePROM(s2_re, "s2");
  op.TruncatePROM(n + 1);
  CHECK(op.GetReducedDimension() == n + 1);
  CHECK(CloseProjectedState(op.GetProjectedState(), ref_re.GetProjectedState()));
  CHECK(op.GetProjectedState().v_node_label == ref_re.GetProjectedState().v_node_label);
  {
    const auto x = m.Solve(op, m.omega), x_ref = m.Solve(ref_re, m.omega);
    for (std::size_t k = 0; k < x.size(); k++)
    {
      CHECK(m.RelDiff(x[k], x_ref[k]) <= 1.0e-10);
    }
  }

  // Truncated to the first n vectors: the saved state, bit for bit.
  op.TruncatePROM(n);
  CHECK(SameProjectedState(op.GetProjectedState(), state_n));
  {
    const auto x = m.Solve(op, m.omega);
    for (std::size_t k = 0; k < x.size(); k++)
    {
      CHECK(m.RelDiff(x[k], x_n[k]) <= 1.0e-12);
    }
  }
  CHECK(CloseProjectedState(op.GetProjectedState(), ref.GetProjectedState()));

  // A sample after the truncation, which adds as many vectors as were dropped.
  op.UpdatePROM(s3, "s3");
  ref.UpdatePROM(s3, "s3");
  REQUIRE(op.GetReducedDimension() == n + 2);
  CHECK(CloseProjectedState(op.GetProjectedState(), ref.GetProjectedState()));
  CHECK(op.GetProjectedState().v_node_label == ref.GetProjectedState().v_node_label);
  CHECK(m.MaxSolutionDiff() <= 1.0e-10);
  CHECK(op.GetProjectedState().V_wp_dim == ref.GetProjectedState().V_wp_dim);
}

// RefineHDM adds to the solution of the last HDM solve the correction for its true
// residual, solved to a relative tolerance of 1e-2 (here looser than the 1e-3 of the
// solve): each step reduces the error against a tight solve and the residual, two steps
// reduce both by orders of magnitude, with every kind of term (the system kept from the
// solve has every term: the Floquet port with wave port model, ten times larger than the
// others, has none which they do not), and it leaves the tolerances of the linear solver
// as they were (the correction has no absolute tolerance: here, one above every
// residual). The system is assembled with the operators of the RomOperator, or, for the
// wave port box (order 2) with the default sparse direct coarse solver, applied through the
// multigrid preconditioner, which is kept with it (see SpaceOperator::GetSystemOperator).
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-RefineHDM",
                 "[romoperator][Serial][Parallel][GPU]")
{
  std::vector<std::pair<std::string, std::string>> configs;
  for (const auto &model : PROM_MODELS)
  {
    if (model != "Floquet and wave ports")
    {
      configs.emplace_back(model, "AMS");
    }
  }
  configs.emplace_back("wave ports", "Default");
  const auto [model, linear_type] = GENERATE_REF(from_range(configs));
  CAPTURE(model, linear_type);
  PROMModel m(temp_dir, model, linear_type);
  auto &op = *m.op;
  const int e = m.excitations.front();
  ComplexVector u = m.NewVector(), u_tight = m.NewVector();
  m.ref->SetLinearTol(1.0e-12);
  m.ref->SolveHDM(e, m.omega, u_tight);
  op.SetLinearTol(1.0e-3);
  op.SolveHDM(e, m.omega, u);
  op.SetLinearAbsTol(1.0e100);
  const double res_loose = op.HDMResidual(u), err_loose = m.RelDiff(u, u_tight);
  REQUIRE(err_loose > 1.0e-6);
  double res = res_loose, err = err_loose;
  for (int k = 1; k <= 2; k++)
  {
    CAPTURE(k, res, err);
    op.RefineHDM(u);
    const double res_k = op.HDMResidual(u), err_k = m.RelDiff(u, u_tight);
    CAPTURE(res_k, err_k);
    CHECK(res_k < res);
    CHECK(err_k < err);
    res = res_k;
    err = err_k;
  }
  CHECK(res < 1.0e-2 * res_loose);
  CHECK(err < 1.0e-2 * err_loose);
  CHECK(op.GetLinearTol() == 1.0e-3);
  CHECK(op.GetLinearAbsTol() == 1.0e100);
}

// AddSamplePROM from an empty basis or from the port modes of circuit synthesis, for each
// excitation (with a shared basis), with the HDM solutions of a loose (1e-3) or a tight
// (1e-12) linear solver tolerance (the other kinds of term are those of the TruncatePROM
// and RefineHDM tests). A reference RomOperator gets each sample as AddSamplePROM leaves
// it with UpdatePROM: the PROMs agree, the error before the update is that of the vectors,
// and the coefficients are those of the sample in the basis.
//   - The Galerkin PROM reproduces an accurate sample in its basis: a sample of a tight
//     solve is added as it is.
//   - A sample of a loose solve, from the first ones on (with an empty basis, or with the
//     port modes of circuit synthesis only), is not reproduced to a tight tolerance: it is
//     removed, refined and added again, until it is reproduced or max_refine times (always,
//     with a zero tolerance), which leaves the PROM with the refined sample only.
//   - A sample of a loose solve already in the basis adds no basis vector, and is refined
//     as well; the refined sample adds new basis vectors.
//   - A part of a sample in the span of the basis adds no basis vector, instead of
//     aborting (see the RomOperator-UpdatePROM tests).
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-AddSamplePROM",
                 "[romoperator][Serial][Parallel][GPU]")
{
  const std::string model =
      GENERATE(as<std::string>{}, "lumped port", "lumped port synthesis");
  CAPTURE(model);
  PROMModel m(temp_dir, model);
  auto &op = *m.op;
  auto &ref = *m.ref;
  constexpr double loose = 1.0e-3, tight = 1.0e-12, tol = 1.0e-8;
  ComplexVector u = m.NewVector(), x0 = m.NewVector();

  // A random direction, with which the PROM is far from the samples.
  auto add_random = [&]()
  {
    const ComplexVector v = m.Random(81);
    op.UpdatePROM(v, "random");
    ref.UpdatePROM(v, "random");
  };

  // AddSamplePROM of the HDM solution u of the excitation e at ω, and UpdatePROM of u as it
  // is then (or of u_ref) on ref.
  auto add_sample = [&](int e, double w, const std::string &label, double t,
                        const ComplexVector *u_ref = nullptr)
  {
    op.SolvePROM(w);
    op.ExpandPROMSolution(e, x0);
    const auto update = op.AddSamplePROM(e, w, u, label, t);
    ref.UpdatePROM(u_ref ? *u_ref : u, label);
    CHECK_THAT(update.error,
               WithinRel(m.RelDiff(x0, u), 1.0e-8) || WithinAbs(m.RelDiff(x0, u), 1.0e-12));
    CHECK(ProjectionResidual(m.comm, op, u, update.c) <= 1.0e-10);
    CHECK_THAT(op.BasisNorm(update.c), WithinRel(m.Norm(u), 1.0e-10));
    CHECK(CloseProjectedState(op.GetProjectedState(), ref.GetProjectedState()));
    CHECK(op.GetProjectedState().v_node_label == ref.GetProjectedState().v_node_label);
    CHECK(m.MaxSolutionDiff() <= 1.0e-10);
    return update;
  };

  SECTION("Samples of a tight solve")
  {
    add_random();
    op.SetLinearTol(tight);
    for (int e : m.excitations)
    {
      CAPTURE(e);
      op.SolveHDM(e, m.omega, u);
      ComplexVector u_hdm(u);
      u_hdm.UseDevice(true);
      const auto update = add_sample(e, m.omega, fmt::format("sample_e{:d}", e), tol);
      CHECK(update.num_refined == 0);
      CHECK(update.reproduced);
      CHECK(m.RelDiff(u, u_hdm) == 0.0);
    }
  }

  SECTION("Samples of a loose solve")
  {
    op.SetLinearTol(loose);
    ref.SetLinearTol(tight);
    ComplexVector u_tight = m.NewVector();
    for (int e : m.excitations)
    {
      CAPTURE(e);
      ref.SolveHDM(e, m.omega, u_tight);
      op.SolveHDM(e, m.omega, u);
      const double err_loose = m.RelDiff(u, u_tight);
      const auto update = add_sample(e, m.omega, fmt::format("sample_e{:d}", e), tol);
      CHECK(update.num_refined >= 1);
      CHECK(update.num_refined <= RomOperator::max_refine);
      CHECK(update.reproduced);
      CHECK(m.RelDiff(u, u_tight) <= 1.0e-2 * err_loose);
    }

    // With a zero tolerance, the PROM never reproduces the sample.
    const int e = m.excitations.front();
    op.SolveHDM(e, m.omega_other, u);
    const auto update = add_sample(e, m.omega_other, "other", 0.0);
    CHECK(update.num_refined == RomOperator::max_refine);
    CHECK(!update.reproduced);
  }

  SECTION("Inaccurate sample in the span of the basis")
  {
    // The PROM with a sample of a loose solve in its basis does not reproduce it to a tight
    // tolerance, also when the sample is added again (and adds no basis vector then).
    op.SetLinearTol(loose);
    ref.SetLinearTol(tight);
    const int e = m.excitations.front();
    ComplexVector u_tight = m.NewVector();
    ref.SolveHDM(e, m.omega, u_tight);
    op.SolveHDM(e, m.omega, u);
    op.UpdatePROM(u, "part");
    ref.UpdatePROM(u, "part");
    const std::size_t n = op.GetReducedDimension();
    const double err_loose = m.RelDiff(u, u_tight);
    const auto update = add_sample(e, m.omega, "sample", tol);
    CHECK(update.num_refined >= 1);
    CHECK(update.reproduced);
    CHECK(op.GetReducedDimension() > n);
    CHECK(m.RelDiff(u, u_tight) < err_loose);
  }

  SECTION("Part in the span of the basis")
  {
    // The imaginary part of a sample in the basis: its real part is added alone.
    add_random();
    op.SetLinearTol(tight);
    const int e = m.excitations.front();
    op.SolveHDM(e, m.omega, u);
    ComplexVector u_re(u), u_im(u);
    u_re.UseDevice(true);
    u_im.UseDevice(true);
    u_re.Imag() = 0.0;
    u_im.Real() = 0.0;
    op.UpdatePROM(u_im, "part");
    ref.UpdatePROM(u_im, "part");
    const std::size_t n = op.GetReducedDimension();
    const auto update = add_sample(e, m.omega, "sample", tol, &u_re);
    CHECK(update.num_refined == 0);
    CHECK(op.GetReducedDimension() == n + 1);
    CHECK(op.GetProjectedState().v_node_label.back() == "sample_re");
    CHECK(m.MaxSolutionDiff() <= 1.0e-10);
  }
}

// Several wave ports can observe different modes of the same boundary, for example
// inactive ports extracting the S-parameters of a second mode at the boundaries of the
// active ports. Circuit synthesis must weight such a shared port boundary once in its inner
// product, as for a single port, instead of aborting on the repeated boundary attribute.
TEST_CASE_METHOD(palace::test::SharedTempDir,
                 "RomOperator-Synthesis-WavePortSharedBoundary",
                 "[romoperator][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  const auto mesh_path =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_tet.msh";
  const auto domain_orthog =
      GENERATE(as<std::string>{}, "Energy", "SpaceOverlap", "FEBasisIdentity");
  CAPTURE(domain_orthog);

  // The box and port faces of the wave port pairing test above: wave ports on the
  // 1 cm x 1 cm faces 14 (active) and 18 (inactive), with PEC on every other face. The
  // shared configuration adds inactive second-mode ports on both faces.
  json pec_attributes = json::array();
  for (int attr = 1; attr <= 22; attr++)
  {
    if (attr != 14 && attr != 18)
    {
      pec_attributes.push_back(attr);
    }
  }
  auto wave_port = [](int index, int attr, int mode, bool active)
  {
    json port = {{"Index", index},
                 {"Attributes", json::array({attr})},
                 {"Mode", mode},
                 {"Offset", 0.0}};
    if (active)
    {
      port["Excitation"] = uint(1);
    }
    else
    {
      port["Active"] = false;
    }
    return port;
  };

  // Inner products of fixed random vectors under the synthesis weight, and the global
  // number of port true dofs it removes from the domain part.
  auto weights = [&](bool shared)
  {
    json ports = json::array({wave_port(1, 14, 1, true), wave_port(2, 18, 1, false)});
    if (shared)
    {
      ports.push_back(wave_port(3, 14, 2, false));
      ports.push_back(wave_port(4, 18, 2, false));
    }
    json setup_json;
    setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
    setup_json["Model"] = {{"Mesh", mesh_path},
                           {"L0", 1.0e-2},
                           {"Refinement", json::object({})},
                           {"CrackInternalBoundaryElements", false}};
    setup_json["Domains"] = {
        {"Materials",
         json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                    {"Permeability", 1.0},
                                    {"Permittivity", 1.0},
                                    {"LossTan", 0.0}})})}};
    setup_json["Boundaries"] = {{"PEC", {{"Attributes", pec_attributes}}},
                                {"WavePort", ports}};
    setup_json["Solver"] = {
        {"Order", 2UL},
        {"Device", "CPU"},
        {"Driven",
         {{"AdaptiveTol", 1.0e-3},
          {"AdaptiveCircuitSynthesis", true},
          {"AdaptiveCircuitSynthesisDomainOrthogonalization", domain_orthog},
          {"MinFreq", 20.0},
          {"MaxFreq", 30.0},
          {"FreqStep", 10.0}}},
        {"Linear",
         {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 200}, {"Tol", 1.0e-10}}}};

    IoData iodata(setup_json, false);
    auto mesh = LoadScaleParMesh2(iodata, comm);
    SpaceOperator space_op(iodata, mesh);
    RomOperatorTest rom_op(iodata, space_op, 4);
    const auto &weight_op = rom_op.GetWeightOp();
    REQUIRE(weight_op.has_value());

    Vector x(space_op.GetNDSpace().GetTrueVSize()), y(x.Size()), work(x.Size());
    x.UseDevice(true);
    y.UseDevice(true);
    work.UseDevice(true);
    linalg::SetRandom(comm, x, 42);
    linalg::SetRandom(comm, y, 43);
    const std::array<double, 3> dots = {weight_op->InnerProduct(comm, x, x, work),
                                        weight_op->InnerProduct(comm, x, y, work),
                                        weight_op->InnerProduct(comm, y, y, work)};
    int n_port_tdofs = weight_op->port_tdof_list.Size();
    Mpi::GlobalSum(1, &n_port_tdofs, comm);
    return std::make_pair(dots, n_port_tdofs);
  };
  const auto [dots_ref, n_port_tdofs_ref] = weights(false);
  const auto [dots, n_port_tdofs] = weights(true);
  REQUIRE(n_port_tdofs_ref > 0);
  CHECK(n_port_tdofs == n_port_tdofs_ref);
  for (std::size_t i = 0; i < dots.size(); i++)
  {
    CHECK_THAT(dots[i], WithinRel(dots_ref[i], 1.0e-12));
  }
}

// Adaptive PROM and circuit synthesis with a frequency-dependent PML: the PROM reproduces
// the HDM solutions at the sample frequencies, the online rational fit of the PML terms
// matches their per-frequency projection, and the synthesized pencil realizes the PML terms
// with decaying auxiliary states.
TEST_CASE_METHOD(palace::test::SharedTempDir, "RomOperator-FrequencyDependentPML",
                 "[romoperator][pml][Serial][Parallel]")
{
  MPI_Comm world_comm = Mpi::World();

  // Box [0, 1] x [0, 1] x [0, 2] with a PML layer for z > 1.5 (attribute 2), a lumped port
  // on the z = 0 face (boundary attribute 1), and PEC elsewhere.
  const auto mesh_path = temp_dir / "pml_box.mesh";
  if (Mpi::Root(world_comm))
  {
    auto smesh =
        mfem::Mesh::MakeCartesian3D(2, 2, 8, mfem::Element::HEXAHEDRON, 1.0, 1.0, 2.0);
    for (int i = 0; i < smesh.GetNE(); i++)
    {
      mfem::Vector center;
      smesh.GetElementCenter(i, center);
      smesh.SetAttribute(i, (center(2) > 1.5) ? 2 : 1);
    }
    smesh.SetAttributes();
    std::ofstream fo(mesh_path);
    fo.precision(17);
    smesh.Print(fo);
  }
  Mpi::Barrier(world_comm);

  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", temp_dir}};
  setup_json["Model"] = {{"Mesh", mesh_path}, {"L0", 0.05}};
  setup_json["Domains"] = {
      {"Materials",
       json::array({json::object({{"Attributes", {1}}, {"Permittivity", 2.0}}),
                    json::object({{"Attributes", {2}}, {"Permittivity", 1.0}})})},
      {"PML", json::array({json::object({{"Attributes", {2}},
                                         {"Direction", {"+Z"}},
                                         {"Thickness", 0.5},
                                         {"FrequencyDependent", true}})})}};
  setup_json["Boundaries"] = {
      {"PEC", json::object({{"Attributes", {2, 3, 4, 5, 6}}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", 1},
                                                {"Attributes", {1}},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {
      {"Order", 2},
      {"Device", "CPU"},
      {"Driven",
       {{"AdaptiveTol", 1.0e-3},
        {"AdaptiveCircuitSynthesis", true},
        {"MinFreq", 4.0},
        {"MaxFreq", 8.0},
        {"FreqStep", 0.05}}},
      {"Linear",
       {{"Type", "Default"}, {"KSPType", "GMRES"}, {"MaxIts", 500}, {"Tol", 1.0e-12}}}};
  IoData iodata(setup_json, false);
  auto mesh_io = LoadScaleParMesh2(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  REQUIRE(space_op.GetMaterialOp().HasFrequencyDependentPML());
  RomOperatorTest prom_op(iodata, space_op, 20);
  prom_op.AddLumpedPortModesForSynthesis();

  const auto &f = iodata.solver.driven.sample_f;
  const std::array<double, 4> omegas = {f.front(), f[f.size() / 3], f[2 * f.size() / 3],
                                        f.back()};
  const int excitation_idx = 1;
  std::vector<ComplexVector> E_hdm;
  for (std::size_t i = 0; i < omegas.size(); i++)
  {
    auto &E = E_hdm.emplace_back(space_op.GetNDSpace().GetTrueVSize());
    E.UseDevice(true);
    prom_op.SolveHDM(excitation_idx, omegas[i], E);
    prom_op.UpdatePROM(E, fmt::format("sample_{:d}", i));
  }

  // Galerkin reproduction at the sample frequencies (per-frequency PML projection).
  ComplexVector E_prom(E_hdm[0].Size());
  E_prom.UseDevice(true);
  for (std::size_t i = 0; i < omegas.size(); i++)
  {
    prom_op.SolvePROM(omegas[i]);
    prom_op.ExpandPROMSolution(excitation_idx, E_prom);
    E_prom.Add(-1.0, E_hdm[i]);
    CHECK_THAT(linalg::Norml2(world_comm, E_prom) / linalg::Norml2(world_comm, E_hdm[i]),
               WithinAbs(0.0, 1.0e-6));
  }

  // Synthesis: the PML terms add auxiliary states for this band.
  const auto matrices = prom_op.CalculateNormalizedPROMMatrices(iodata.units);
  CHECK(std::ranges::any_of(matrices.aux_labels,
                            [](const auto &label) { return label.starts_with("pml_p"); }));
  CHECK(matrices.L_inv->allFinite());
  CHECK(matrices.C->allFinite());

  // Online rational fit of the PML terms vs. their per-frequency projection.
  prom_op.PrepareOnlinePhase();
  REQUIRE(prom_op.HasPMLOnlineFit());
  const double omega_test = 0.5 * (omegas[1] + omegas[2]);
  ComplexVector E_fit(E_prom.Size());
  E_fit.UseDevice(true);
  prom_op.SolvePROM(omega_test);
  prom_op.ExpandPROMSolution(excitation_idx, E_fit);
  prom_op.ResetPMLOnlineFit();
  prom_op.SolvePROM(omega_test);
  prom_op.ExpandPROMSolution(excitation_idx, E_prom);
  const double norm = linalg::Norml2(world_comm, E_prom);
  E_fit.Add(-1.0, E_prom);
  CHECK_THAT(linalg::Norml2(world_comm, E_fit) / norm, WithinAbs(0.0, 1.0e-8));
}

// Excited ports must always be included in synthesis. The configuration parser
// rejects a port carrying both "Excitation" > 0 and "IncludeInSynthesis": false because
// the excitation vector is unconditionally added to the basis.
TEST_CASE("RomOperator-Synthesis-ExcludedExcitedRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"IncludeInSynthesis", false},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  CHECK_THROWS_WITH((IoData{setup_json, false}),
                    Catch::Matchers::ContainsSubstring("IncludeInSynthesis"));
}

// Same rule for wave ports: an excited wave port carrying "IncludeInSynthesis": false is
// rejected by the configuration parser, because the excitation vector is unconditionally
// added to the synthesis basis (mirrors the lumped-port case above).
TEST_CASE("RomOperator-Synthesis-WavePortExcludedExcitedRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  setup_json["Boundaries"] = {
      {"WavePort", json::array({json::object({{"Index", 1},
                                              {"Mode", 1},
                                              {"Excitation", uint(1)},
                                              {"IncludeInSynthesis", false},
                                              {"Attributes", json::array({100})}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  CHECK_THROWS_WITH((IoData{setup_json, false}),
                    Catch::Matchers::ContainsSubstring("IncludeInSynthesis"));
}

// A Floquet wave vector makes the periodic problem non-reciprocal, so that the synthesized
// circuit matrices would not be symmetric: the configuration parser rejects circuit
// synthesis with a nonzero wave vector, but not periodic boundaries without a phase shift.
TEST_CASE("RomOperator-Synthesis-FloquetWaveVectorRejected", "[romoperator][Serial]")
{
  json setup_json;
  setup_json["Problem"] = {{"Type", "Driven"}, {"Verbose", 0}, {"Output", "."}};
  setup_json["Model"] = {{"Mesh", "placeholder.msh"}};
  setup_json["Domains"] = {
      {"Materials", json::array({json::object({{"Attributes", json::array({1})},
                                               {"Permeability", 1.0},
                                               {"Permittivity", 1.0}})})}};
  const auto wave_vector =
      GENERATE(json::array({0.0, 0.0, 0.0}), json::array({0.0, 0.1, 0.0}));
  setup_json["Boundaries"] = {
      {"Periodic",
       {{"FloquetWaveVector", wave_vector},
        {"BoundaryPairs",
         json::array({json::object({{"DonorAttributes", json::array({1})},
                                    {"ReceiverAttributes", json::array({2})}})})}}},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({100})},
                                                {"Direction", "+X"}})})}};
  setup_json["Solver"] = {{"Order", 1UL},
                          {"Driven",
                           {{"AdaptiveTol", 1e-3},
                            {"AdaptiveCircuitSynthesis", true},
                            {"MinFreq", 2.0},
                            {"MaxFreq", 32.0},
                            {"FreqStep", 1.0}}}};

  INFO("FloquetWaveVector: " << wave_vector.dump());
  if (wave_vector[1] == 0.0)
  {
    CHECK_NOTHROW((IoData{setup_json, false}));
  }
  else
  {
    CHECK_THROWS_WITH((IoData{setup_json, false}),
                      Catch::Matchers::ContainsSubstring("FloquetWaveVector"));
  }
}

// The polynomial dispersion correction helper is the sign bridge between a scalar
// frequency-dependent boundary multiplier and the quadratic pencil K + iωC - ω²M. Lock
// this down directly: evaluating the corrected pencil must reproduce
// (α₀ + α₁ω + α₂ω²) Mp_r for real coefficients (a lossless wave port) and complex ones (a
// lossy wave port, or a scalar-dispersion fit of another BC), also at complex ω (the
// synthesized eigenvalues).
TEST_CASE("RomOperator-PolynomialFitCorrections", "[romoperator][Serial]")
{
  Eigen::MatrixXcd Mp_r(2, 2);
  Mp_r << std::complex<double>(0.0, 2.0), std::complex<double>(0.0, -0.3),
      std::complex<double>(0.0, -0.3), std::complex<double>(0.0, 0.7);

  for (const auto &[alpha0, alpha1, alpha2] :
       {std::array<std::complex<double>, 3>{1.2, -0.4, 0.05},
        std::array<std::complex<double>, 3>{std::complex<double>(0.4, -0.2),
                                            std::complex<double>(-0.1, 0.3),
                                            std::complex<double>(0.02, -0.04)}})
  {
    Eigen::MatrixXcd Kr = Eigen::MatrixXcd::Zero(2, 2);
    Eigen::MatrixXcd Cr = Eigen::MatrixXcd::Zero(2, 2);
    Eigen::MatrixXcd Mr = Eigen::MatrixXcd::Zero(2, 2);
    RomOperatorTest::ApplyComplexPolynomialFitCorrections(alpha0, alpha1, alpha2, Mp_r, Kr,
                                                          Cr, Mr);
    for (std::complex<double> omega :
         {std::complex<double>(0.2, 0.0), std::complex<double>(1.5, 0.0),
          std::complex<double>(3.0, 0.0), std::complex<double>(1.5, 0.4)})
    {
      const Eigen::MatrixXcd got =
          Kr + std::complex<double>(0.0, 1.0) * omega * Cr - omega * omega * Mr;
      const Eigen::MatrixXcd want =
          (alpha0 + alpha1 * omega + alpha2 * omega * omega) * Mp_r;
      CHECK(((got - want).cwiseAbs().maxCoeff() < 1.0e-13));
    }
  }
}

// Modal residue matrices are real symmetric but generally indefinite. Verify the auxiliary
// direction factorization preserves signed eigenvalues rather than replacing M by |M|.
TEST_CASE("RomOperator-AugmentedPencil-IndefiniteDirections", "[romoperator][Serial]")
{
  Eigen::MatrixXd M(2, 2);
  M << 0.0, 1.0, 1.0, 0.0;  // eigenvalues +1 and -1
  const Eigen::MatrixXcd Mp_r =
      std::complex<double>(0.0, 1.0) * M.cast<std::complex<double>>();
  RomOperatorTest::WavePortAuxBlock blk;
  REQUIRE(RomOperatorTest::AddAuxBlockDirections(blk, Mp_r, 1.0e-12));
  REQUIRE(blk.weights.size() == 2);
  REQUIRE(blk.u_dirs.size() == 2);

  Eigen::MatrixXd reconstructed = Eigen::MatrixXd::Zero(2, 2);
  bool has_positive = false, has_negative = false;
  for (std::size_t j = 0; j < blk.weights.size(); j++)
  {
    reconstructed.noalias() += blk.weights[j] * blk.u_dirs[j] * blk.u_dirs[j].transpose();
    has_positive = has_positive || (blk.weights[j] > 0.0);
    has_negative = has_negative || (blk.weights[j] < 0.0);
  }
  CHECK(has_positive);
  CHECK(has_negative);
  CHECK((reconstructed - M).cwiseAbs().maxCoeff() < 1.0e-13);
}

// Modal-correction residue matrices can contain physical signed directions well below the
// fixed wave-port boundary-mass noise floor. The caller-provided tolerance must control
// their factorization; applying the old 3e-3 mass floor here drops the second direction and
// changes the augmented pencil's Schur complement by O(1e-4).
TEST_CASE("RomOperator-AugmentedPencil-ModalResidueRank", "[romoperator][Serial]")
{
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, 2);
  M(0, 0) = 1.0;
  M(1, 1) = -1.0e-4;
  const Eigen::MatrixXcd Mp_r =
      std::complex<double>(0.0, 1.0) * M.cast<std::complex<double>>();
  RomOperatorTest::WavePortAuxBlock blk;
  REQUIRE(RomOperatorTest::AddAuxBlockDirections(blk, Mp_r, 1.0e-6));
  REQUIRE(blk.weights.size() == 2);

  Eigen::MatrixXd reconstructed = Eigen::MatrixXd::Zero(2, 2);
  for (std::size_t j = 0; j < blk.weights.size(); j++)
  {
    reconstructed.noalias() += blk.weights[j] * blk.u_dirs[j] * blk.u_dirs[j].transpose();
  }
  CHECK((reconstructed - M).cwiseAbs().maxCoeff() < 1.0e-13);
}

// Schur-complement test for the augmented rational realization. This exercises the real
// BuildAugmentedPencil helper with multiple poles and signed symmetric directions, then
// eliminates the aux rows and checks that the effective v-block exactly equals the
// unaugmented rational pencil contribution i Σ r_k/(ω-p_k) M_proj.
TEST_CASE("RomOperator-AugmentedPencil-SchurComplement", "[romoperator][Serial]")
{
  Eigen::MatrixXcd Kr(2, 2), Cr(2, 2), Mr(2, 2);
  Kr << std::complex<double>(1.0, 0.2), std::complex<double>(0.15, -0.05),
      std::complex<double>(0.15, -0.05), std::complex<double>(1.7, -0.1);
  Cr << std::complex<double>(0.05, -0.02), std::complex<double>(0.01, 0.03),
      std::complex<double>(0.01, 0.03), std::complex<double>(-0.04, 0.01);
  Mr << std::complex<double>(0.3, 0.0), std::complex<double>(0.02, -0.01),
      std::complex<double>(0.02, -0.01), std::complex<double>(0.6, 0.02);

  const double theta = 0.37;
  Eigen::VectorXd u0(2), u1(2);
  u0 << std::cos(theta), std::sin(theta);
  u1 << -std::sin(theta), std::cos(theta);
  RomOperatorTest::WavePortAuxBlock blk;
  blk.port_idx = 3;
  blk.weights = {2.0, -0.4};
  blk.u_dirs = {u0, u1};
  blk.poles = {std::complex<double>(1.1, 0.2), std::complex<double>(-0.7, 0.4)};
  blk.residues = {std::complex<double>(0.3, -0.15), std::complex<double>(-0.2, 0.05)};

  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kr, Cr, Mr, {blk}, aux_labels);
  REQUIRE(aug.Kr.rows() == 6);
  REQUIRE(aug.Kr.cols() == 6);
  REQUIRE(aux_labels.size() == 4);
  CHECK(aux_labels.front() == "waveport_3_p0d0");
  CHECK(aux_labels.back() == "waveport_3_p1d1");

  Eigen::MatrixXcd M_proj =
      blk.weights[0] * (u0 * u0.transpose()).cast<std::complex<double>>() +
      blk.weights[1] * (u1 * u1.transpose()).cast<std::complex<double>>();

  for (std::complex<double> omega :
       {std::complex<double>(0.5, 0.1), std::complex<double>(1.6, -0.2),
        std::complex<double>(2.4, 0.0)})
  {
    const Eigen::MatrixXcd A_aug =
        aug.Kr + std::complex<double>(0.0, 1.0) * omega * aug.Cr - omega * omega * aug.Mr;
    const Eigen::MatrixXcd Avv = A_aug.topLeftCorner(2, 2);
    const Eigen::MatrixXcd Ava = A_aug.topRightCorner(2, 4);
    const Eigen::MatrixXcd Aav = A_aug.bottomLeftCorner(4, 2);
    const Eigen::MatrixXcd Aaa = A_aug.bottomRightCorner(4, 4);
    const Eigen::MatrixXcd schur = Avv - Ava * Aaa.fullPivLu().solve(Aav);

    Eigen::MatrixXcd want =
        Kr + std::complex<double>(0.0, 1.0) * omega * Cr - omega * omega * Mr;
    std::complex<double> rational = 0.0;
    for (std::size_t k = 0; k < blk.poles.size(); k++)
    {
      rational += blk.residues[k] / (omega - blk.poles[k]);
    }
    want += std::complex<double>(0.0, 1.0) * rational * M_proj;
    CHECK(((schur - want).cwiseAbs().maxCoeff() < 1.0e-12));
  }
}

// Complex (Takagi) aux directions, as used for the frequency-dependent PML: eliminating the
// aux states reproduces i·rₖ·Σⱼ σⱼ uⱼuⱼᵀ / (ω − pₖ) for complex uⱼ, so a complex symmetric
// residue matrix R = Σⱼ σⱼ uⱼuⱼᵀ is realized with one aux state per direction.
TEST_CASE("RomOperator-AugmentedPencil-ComplexDirections", "[romoperator][Serial]")
{
  Eigen::MatrixXcd Kr(3, 3), Cr = Eigen::MatrixXcd::Zero(3, 3), Mr(3, 3);
  Kr << 2.0, 0.1, 0.0, 0.1, 3.0, 0.2, 0.0, 0.2, 1.5;
  Mr = Eigen::MatrixXcd::Identity(3, 3);

  // Complex symmetric residue and its Takagi factorization from the real symmetric
  // embedding [[Re R, Im R], [Im R, -Re R]].
  Eigen::MatrixXcd R(3, 3);
  R << std::complex<double>(1.0, 0.5), std::complex<double>(0.2, -0.3),
      std::complex<double>(0.0, 0.1), std::complex<double>(0.2, -0.3),
      std::complex<double>(-0.4, 0.2), std::complex<double>(0.3, 0.0),
      std::complex<double>(0.0, 0.1), std::complex<double>(0.3, 0.0),
      std::complex<double>(0.7, -0.6);
  Eigen::MatrixXd E(6, 6);
  E << R.real(), R.imag(), R.imag(), -R.real();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(E);
  RomOperatorTest::WavePortAuxBlock blk;
  blk.port_idx = -1;
  blk.label = "pml_p0";
  blk.poles = {std::complex<double>(0.0, 2.5)};
  blk.residues = {std::complex<double>(0.0, -1.0)};
  for (int j = 0; j < 6; j++)
  {
    if (eig.eigenvalues()(j) > 1.0e-12)
    {
      blk.weights.push_back(eig.eigenvalues()(j));
      blk.u_dirs_c.push_back(
          eig.eigenvectors().col(j).head(3).cast<std::complex<double>>() +
          std::complex<double>(0.0, 1.0) *
              eig.eigenvectors().col(j).tail(3).cast<std::complex<double>>());
    }
  }
  REQUIRE(blk.weights.size() == 3);

  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kr, Cr, Mr, {blk}, aux_labels);
  REQUIRE(aug.Kr.rows() == 6);
  CHECK(aux_labels.front() == "pml_p0_p0d0");
  for (std::complex<double> omega :
       {std::complex<double>(0.5, 0.0), std::complex<double>(1.6, 0.3),
        std::complex<double>(4.0, -0.1)})
  {
    const Eigen::MatrixXcd A_aug =
        aug.Kr + std::complex<double>(0.0, 1.0) * omega * aug.Cr - omega * omega * aug.Mr;
    const Eigen::MatrixXcd schur =
        A_aug.topLeftCorner(3, 3) -
        A_aug.topRightCorner(3, 3) *
            A_aug.bottomRightCorner(3, 3).fullPivLu().solve(A_aug.bottomLeftCorner(3, 3));
    const Eigen::MatrixXcd want = Kr - omega * omega * Mr + R / (omega - blk.poles[0]);
    CHECK(((schur - want).cwiseAbs().maxCoeff() < 1.0e-12));
  }
}

// Verify that the regime-2 augmented pencil produced by BuildAugmentedPencil has the
// same eigenvalues as the unaugmented pencil with the rational kₙ-correction applied
// directly. The augmentation is supposed to be exact: eliminating each aux state sₖⱼ
// from the linear pencil reproduces (i·rₖ·σⱼ·uⱼuⱼᵀ)/(ω − pₖ) in the v-row, so the
// eigenvalues of the augmented pencil are exactly the roots of det(K + iωC − ω²M +
// i·Σₖ rₖ/(ω − pₖ)·M_proj) = 0. This test catches sign errors in the coupling formula
// √(−i·rₖ·σⱼ), wrong K_aa/C_aa diagonal entries, and any off-by-one in the aux row
// layout.
TEST_CASE("RomOperator-AugmentedPencil-Eigenvalues", "[romoperator][Serial]")
{
  // 1×1 base pencil. Treat the basis-row entries as the "v" block: K=1, C=0, M=0, and
  // M_proj = 1·u·uᵀ with u=[1], σ=1. Two poles with hand-picked complex residues so the
  // analytical eigenvalue equation is non-trivial in both real and imaginary parts.
  Eigen::MatrixXcd Kr(1, 1), Cr(1, 1), Mr(1, 1);
  Kr(0, 0) = std::complex<double>(1.0, 0.0);
  Cr(0, 0) = std::complex<double>(0.0, 0.0);
  Mr(0, 0) = std::complex<double>(0.0, 0.0);

  RomOperatorTest::WavePortAuxBlock blk;
  blk.port_idx = 1;
  blk.weights = {1.0};
  Eigen::VectorXd u_dir(1);
  u_dir(0) = 1.0;
  blk.u_dirs = {u_dir};
  blk.poles = {std::complex<double>(2.0, 0.1), std::complex<double>(-1.5, 0.05)};
  blk.residues = {std::complex<double>(0.5, 0.2), std::complex<double>(-0.3, 0.1)};

  std::vector<RomOperatorTest::WavePortAuxBlock> aux_blocks = {blk};
  std::vector<std::string> aux_labels;
  auto aug = RomOperatorTest::BuildAugmentedPencil(Kr, Cr, Mr, aux_blocks, aux_labels);

  REQUIRE(aug.Kr.rows() == 3);
  REQUIRE(aug.Kr.cols() == 3);
  REQUIRE(aux_labels.size() == 2);

  // Linearize the QEP into a 6×6 generalized eigenvalue problem and solve.
  // (K + iωC − ω²M)·x = 0  ⇔  [iC, K; -I, 0]·[ωx; x] = ω·[M, 0; 0, -I]·[ωx; x].
  const long n = aug.Kr.rows();
  Eigen::MatrixXcd A_lin = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  Eigen::MatrixXcd B_lin = Eigen::MatrixXcd::Zero(2 * n, 2 * n);
  A_lin.topLeftCorner(n, n) = std::complex<double>(0.0, 1.0) * aug.Cr;
  A_lin.topRightCorner(n, n) = aug.Kr;
  A_lin.bottomLeftCorner(n, n) = -Eigen::MatrixXcd::Identity(n, n);
  B_lin.topLeftCorner(n, n) = aug.Mr;
  B_lin.bottomRightCorner(n, n) = -Eigen::MatrixXcd::Identity(n, n);

  // Eigen's GeneralizedEigenSolver works on real matrices only, so we use a dense
  // generalized Schur via ComplexEigenSolver on B⁻¹·A. With M=0 in the v-block, B is
  // singular, so we instead build the explicit polynomial whose roots are the
  // augmented eigenvalues (closed form) and verify the augmented pencil annihilates
  // each.
  //
  // The unaugmented effective pencil is (K + iωC − ω²M) + i·Σₖ rₖ/(ω − pₖ)·M_proj
  // applied at the "v" block (size 1 here). Multiplying by Πₖ(ω − pₖ) gives:
  //   p(ω) = (K + iωC − ω²M)·(ω − p₀)·(ω − p₁) + i·M_proj·[r₀·(ω − p₁) + r₁·(ω − p₀)]
  // For this 1×1 case with K=1, C=0, M=0, M_proj=1:
  //   p(ω) = (ω − p₀)·(ω − p₁) + i·[r₀·(ω − p₁) + r₁·(ω − p₀)]
  // a quadratic in ω. The two roots must equal two of the three augmented-pencil
  // eigenvalues (one extra eigenvalue at infinity from the rank-deficient M-block).
  const auto p0 = blk.poles[0];
  const auto p1 = blk.poles[1];
  const auto r0 = blk.residues[0];
  const auto r1 = blk.residues[1];
  // p(ω) = ω² + b·ω + c with
  //   b = -(p₀ + p₁) + i·(r₀ + r₁)
  //   c =  p₀·p₁ - i·(r₀·p₁ + r₁·p₀)
  const std::complex<double> b = -(p0 + p1) + std::complex<double>(0.0, 1.0) * (r0 + r1);
  const std::complex<double> c =
      p0 * p1 - std::complex<double>(0.0, 1.0) * (r0 * p1 + r1 * p0);
  const std::complex<double> disc = std::sqrt(b * b - 4.0 * c);
  const std::complex<double> root_plus = 0.5 * (-b + disc);
  const std::complex<double> root_minus = 0.5 * (-b - disc);

  // Verify each analytical root annihilates the augmented pencil A_aug(ω) = K + iωC −
  // ω²M when treated as a 3×3 complex matrix (det should vanish).
  auto eval_pencil_det = [&](std::complex<double> w)
  {
    Eigen::MatrixXcd A =
        aug.Kr + std::complex<double>(0.0, 1.0) * w * aug.Cr - w * w * aug.Mr;
    return A.determinant();
  };
  CHECK(std::abs(eval_pencil_det(root_plus)) < 1.0e-10);
  CHECK(std::abs(eval_pencil_det(root_minus)) < 1.0e-10);

  // Also verify the unaugmented effective pencil at the same roots gives the same
  // determinant (= 0) — sanity check on the algebraic equivalence.
  auto eval_unaug = [&](std::complex<double> w)
  {
    std::complex<double> kn = r0 / (w - p0) + r1 / (w - p1);
    return Kr(0, 0) + std::complex<double>(0.0, 1.0) * w * Cr(0, 0) - w * w * Mr(0, 0) +
           std::complex<double>(0.0, 1.0) * kn;
  };
  CHECK(std::abs(eval_unaug(root_plus)) < 1.0e-10);
  CHECK(std::abs(eval_unaug(root_minus)) < 1.0e-10);
}

namespace
{

// Wilkinson's matrix: 1 on the diagonal and in the last column, -1 below the diagonal.
// Partial pivoting does not swap rows and the entries of the last column of U grow as
// 2^(k-1), so the LU solution has a large backward error. The matrix itself is well
// conditioned.
Eigen::MatrixXcd WilkinsonMatrix(int m)
{
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Identity(m, m);
  A.triangularView<Eigen::StrictlyLower>().setConstant(-1.0);
  A.col(m - 1).setOnes();
  return A;
}

}  // namespace

TEST_CASE("RomOperator-SolveReducedSystem", "[romoperator][Serial]")
{
  // LU solves well and ill-conditioned systems, and full-pivoting QR takes over for
  // singular matrices and for LU solutions with a large backward error in any column.
  using Matrix = Eigen::MatrixXcd;
  using VectorC = Eigen::VectorXcd;
  constexpr int n = 12;
  std::srand(1234);
  auto backward_error = [](const Matrix &A, const Matrix &x, const Matrix &b)
  {
    return (b - A * x).lpNorm<Eigen::Infinity>() /
           (A.cwiseAbs().rowwise().sum().maxCoeff() * x.lpNorm<Eigen::Infinity>() +
            b.lpNorm<Eigen::Infinity>());
  };
  auto with_singular_values = [](const Eigen::VectorXd &s)
  {
    Eigen::HouseholderQR<Matrix> qr_u(Matrix::Random(s.size(), s.size())),
        qr_v(Matrix::Random(s.size(), s.size()));
    const Matrix U = qr_u.householderQ(), W = qr_v.householderQ();
    return Matrix(U * s.cast<std::complex<double>>().asDiagonal() * W.adjoint());
  };
  MPI_Comm comm = Mpi::World();
  const VectorC b = VectorC::Random(n);

  SECTION("Conditioning")
  {
    // Condition numbers 1e10 and 1e17 (singular to working precision): the LU solution is
    // backward stable, accurate to cond ε for the first.
    const double s_min = GENERATE(1.0e-10, 1.0e-17);
    CAPTURE(s_min);
    Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(n, 0.0, -3.0)
                            .unaryExpr([](double e) { return std::pow(10.0, e); });
    s(n - 1) = s_min;
    const Matrix A = with_singular_values(s);
    Matrix x = b;
    CHECK_FALSE(SolveReducedSystem(comm, A, x));
    CHECK(backward_error(A, x, b) < 1.0e-14);
    if (s_min > 1.0e-12)
    {
      const VectorC x_ref = A.fullPivHouseholderQr().solve(b);
      CHECK((x - x_ref).norm() < 1.0e-4 * x_ref.norm());
    }
  }

  SECTION("Singular")
  {
    // With a zero column, the LU solution is not finite and the rank-revealing QR gives a
    // solution of norm ~ |b| / 1e-3, the smallest nonzero singular value.
    Eigen::VectorXd s = Eigen::VectorXd::LinSpaced(n, 0.0, -3.0)
                            .unaryExpr([](double e) { return std::pow(10.0, e); });
    Matrix A = with_singular_values(s);
    A.col(3).setZero();
    Matrix x = b;
    CHECK(SolveReducedSystem(comm, A, x));
    CHECK(x.allFinite());
    CHECK(x.norm() < 1.0e4 * b.norm());
  }

  SECTION("Badly scaled circuit-synthesis matrices")
  {
    // Synthetic L⁻¹, R⁻¹, C like those of circuit synthesis, with a port row (a 50 Ω port)
    // and three sample rows, couplings much weaker for the port, and entries spanning about
    // 1e-28 to 1e14: LU is used, and agrees with QR.
    auto symmetric = [](const Eigen::Vector4d &d, double port, double sample)
    {
      Eigen::Matrix4d S;
      for (int i = 0; i < 4; i++)
      {
        for (int j = 0; j < 4; j++)
        {
          S(i, j) = (i == j) ? d(i)
                             : ((i == 0 || j == 0) ? port : sample) *
                                   std::sin(1.0 + i + j) * std::sqrt(d(i) * d(j));
        }
      }
      return S;
    };
    const Eigen::Matrix4d Linv =
        symmetric({3.0e12, 1.7e14, 1.8e14, 2.2e14}, 1.0e-9, 1.0e-4);
    const Eigen::Matrix4d Rinv = symmetric({0.02, 1.1e-12, 6.4e-13, 1.5e-12}, 1.0e-20, 0.1);
    const Eigen::Matrix4d Cm =
        symmetric({2.0e-16, 1.8e-13, 1.8e-13, 1.8e-13}, 1.0e-11, 1.0e-12);
    for (double f : {2.5e9, 6.0e9, 9.9e9})
    {
      const double w = 2.0 * std::numbers::pi * f;
      Matrix A = Linv.cast<std::complex<double>>() +
                 (1i * w) * Rinv.cast<std::complex<double>>() -
                 (w * w) * Cm.cast<std::complex<double>>();
      A /= A.cwiseAbs().maxCoeff();
      const VectorC b4 = VectorC::Random(4);
      Matrix x = b4;
      CHECK_FALSE(SolveReducedSystem(comm, A, x));
      const VectorC x_ref = A.fullPivHouseholderQr().solve(b4);
      CHECK((x - x_ref).norm() < 1.0e-12 * x_ref.norm());
    }
  }

  SECTION("Pivot growth")
  {
    // Partial pivoting on Wilkinson's matrix gives a large backward error.
    constexpr int m = 64;
    const Matrix A = WilkinsonMatrix(m);
    const VectorC bw = VectorC::Random(m);
    VectorC x_lu = A.partialPivLu().solve(bw);
    REQUIRE(!(backward_error(A, x_lu, bw) < 1.0e-12));
    // A zero first column of the right-hand side does not hide the error of the second.
    Matrix x(m, 2);
    x << VectorC::Zero(m), bw;
    CHECK(SolveReducedSystem(comm, A, x));
    CHECK(x.col(0).isZero());
    CHECK(backward_error(A, x.col(1), bw) < 1.0e-14);
  }
}

TEST_CASE("RomOperator-SolveReducedSystem-NonFinite", "[romoperator][Serial]")
{
  // A finite LU solution whose backward-error bound overflows is not verified, so QR is
  // used. (A residual that overflows implies an overflowing bound, |Aᵢⱼ xⱼ| ≤ ‖A‖ ‖x‖.)
  // Here x = (1e-300, 1e10) is finite but ‖A‖ |x| = 1e310.
  using Matrix = Eigen::MatrixXcd;
  MPI_Comm comm = Mpi::World();
  Matrix A = Matrix::Zero(2, 2);
  A(0, 0) = 1.0e300;
  A(1, 1) = 1.0e-10;
  Matrix x = Matrix::Ones(2, 1);
  REQUIRE(Matrix(A.partialPivLu().solve(x)).allFinite());
  CHECK(SolveReducedSystem(comm, A, x));
  // The condition number 1e310 is singular to working precision for the rank-revealing
  // QR, which drops the second column.
  CHECK(x.allFinite());
}

TEST_CASE("RomOperator-SolveReducedSystem-Parallel", "[romoperator][Serial][Parallel]")
{
  // The choice of QR is collective: with Wilkinson's matrix (pivot growth, see
  // RomOperator-SolveReducedSystem) on all ranks and a nonzero right-hand side only on rank
  // 0, LU is accurate (x = 0) on the other ranks, but every rank falls back to QR.
  using Matrix = Eigen::MatrixXcd;
  MPI_Comm comm = Mpi::World();
  constexpr int m = 64;
  const Matrix A = WilkinsonMatrix(m);
  std::mt19937 gen(1234);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Matrix b(m, 1);
  for (int i = 0; i < m; i++)
  {
    b(i) = {dist(gen), dist(gen)};
  }
  if (Mpi::Rank(comm) != 0)
  {
    b.setZero();
  }
  Matrix x = b;
  const int use_qr = SolveReducedSystem(comm, A, x);
  int qr_min = use_qr, qr_max = use_qr;
  Mpi::GlobalMin(1, &qr_min, comm);
  Mpi::GlobalMax(1, &qr_max, comm);
  CHECK(use_qr == 1);
  CHECK(qr_min == qr_max);
  const double backward_error =
      (b - A * x).lpNorm<Eigen::Infinity>() /
      (A.cwiseAbs().rowwise().sum().maxCoeff() * x.lpNorm<Eigen::Infinity>() +
       b.lpNorm<Eigen::Infinity>() + std::numeric_limits<double>::min());
  CHECK(backward_error < 1.0e-14);
  if (Mpi::Rank(comm) != 0)
  {
    CHECK(x.isZero());
  }
}
