// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "eigensolver.hpp"

#include <complex>
#include <mfem.hpp>
#include "fem/errorindicator.hpp"
#include "fem/mesh.hpp"
#include "linalg/arpack.hpp"
#include "linalg/divfree.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/floquetcorrection.hpp"
#include "linalg/ksp.hpp"
#include "linalg/nleps.hpp"
#include "linalg/operator.hpp"
#include "linalg/orthog.hpp" // test for orthog
#include "linalg/slepc.hpp"
#include "linalg/vector.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"
#include <Eigen/Dense> // test for AAA

// Eigen does not provide a complex-valued generalized eigenvalue solver, so we use LAPACK
// for this.
extern "C"
{
  void zggev_(char *, char *, int *, std::complex<double> *, int *, std::complex<double> *,
              int *, std::complex<double> *, std::complex<double> *, std::complex<double> *,
              int *, std::complex<double> *, int *, std::complex<double> *, int *, double *,
              int *);
}

namespace palace
{

using namespace std::complex_literals;

using Complex = std::complex<double>;
using VectorEig = Eigen::VectorXcd;
using Matrix = Eigen::MatrixXcd;

template <typename MatType, typename VecType>
inline void ZGGEV(MatType &A, MatType &B, VecType &D, MatType &VR)
{
  // Wrapper for LAPACK's (z)ggev. A and B are overwritten by their Schur decompositions.
  MFEM_VERIFY(A.rows() == A.cols() && B.rows() == B.cols() && A.rows() == B.rows(),
              "Generalized eigenvalue problem expects A, B matrices to be square and have "
              "same dimensions!");
  char jobvl = 'N', jobvr = 'V';
  int n = static_cast<int>(A.rows()), lwork = 2 * n;
  std::vector<std::complex<double>> alpha(n), beta(n), work(lwork);
  std::vector<double> rwork(8 * n);
  MatType VL(0, 0);
  VR.resize(n, n);
  int info = 0;

  zggev_(&jobvl, &jobvr, &n, A.data(), &n, B.data(), &n, alpha.data(), beta.data(),
         VL.data(), &n, VR.data(), &n, work.data(), &lwork, rwork.data(), &info);
  MFEM_VERIFY(info == 0, "ZGGEV failed with info = " << info << "!");

  // Postprocess the eigenvalues and eigenvectors (return unit 2-norm eigenvectors).
  D.resize(n);
  for (int i = 0; i < n; i++)
  {
    D(i) = (beta[i] == 0.0)
               ? ((alpha[i] == 0.0) ? std::numeric_limits<std::complex<double>>::quiet_NaN()
                                    : mfem::infinity())
               : alpha[i] / beta[i];
    VR.col(i) /= VR.col(i).norm();
  }
}

class AAAApproximation
{
private:
    std::vector<Complex> support_points;  // zi
    std::vector<Complex> weights;         // wi
    std::vector<Complex> function_values; // fi
    std::vector<Complex> poles;           // pi
public:
    struct AAAResult
    {
        std::vector<Complex> support_points;
        std::vector<Complex> weights;
        std::vector<Complex> function_values;
        std::vector<Complex> poles;
        std::vector<double> errors;
        int degree;
    };

    // Evaluate the rational approximation at point z
    Complex evaluate(Complex z) const
    {
        if (support_points.empty()) return Complex(0.0);

        Complex numerator(0.0);
        Complex denominator(0.0);

        for (size_t i = 0; i < support_points.size(); ++i) {
            Complex term = weights[i] / (z - support_points[i]);
            numerator += function_values[i] * term;
            denominator += term;
        }

        return numerator / denominator;
    }

    AAAResult computeAAA(
      std::function<Complex(Complex)> f,
      const std::vector<Complex>& sample_points,
      double tolerance = 1e-13,
      int max_degree = 100)
    {
      int n = sample_points.size();
      std::vector<Complex> f_values(n);
      Complex mean = 0.0;
      double norm = 0.0;
      for (int i = 0; i < n; i++)
      {
        f_values[i] = f(sample_points[i]);
        mean += f_values[i];
        norm = std::max(std::abs(f_values[i]), norm);
      }
      mean /= n;
      std::vector<double> R(n);
      for (int i = 0; i < n; i++) R[i] = std::abs(f_values[i]-mean);

      // Initialize
      support_points.clear();
      weights.clear();
      function_values.clear();
      poles.clear();

      // next support point
      double err;
      err = 0.0;
      int idx = 0;
      for (int i = 0; i < n; i++) {if (R[i] >= err) {idx = i; err = R[i];}}

      std::vector<Complex> z(n), f2(n);
      Matrix C(n,n), A(n,n);
      VectorEig ww(n), D(n), N(n);
      int k;
      for (k = 0; k < n-1; k++)
      {
        z[k] = sample_points[idx];
        f2[k] = f_values[idx];
        R[idx] = -1.0;
        // Next column
        for (int i = 0; i < n; i++)
        {
          C(i,k) = 1.0 / (sample_points[i] - sample_points[idx]);
        }
        int cont = 0;
        for (int i = 0; i < n; i++)
        {
          if (R[i] != -1.0)
          {
            for (int j = 0; j <= k; j++)
            {
              A(cont,j) = C(i,j)*f_values[i] - C(i,j)*f2[j];
            }
            cont++;
          }
        }
        // SVD of A
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd;
        svd.compute(A.topLeftCorner(cont, k+1), Eigen::ComputeFullV);
        for (int i = 0; i <= k; i++)
        {
          ww(i) = std::conj(svd.matrixV().col(k)(i));
          D(i) = ww(i) * f2[i];
        }
        N = C.topLeftCorner(n,k+1) * D.head(k+1);
        D = C.topLeftCorner(n,k+1) * ww.head(k+1);
        for (int i = 0; i < n; i++)
        {
          if (R[i] >= 0) R[i] = std::abs(f_values[i] - N[i]/D[i]);
        }
        // next support point
        err = 0.0;
        double mean_err = 0.0;
        for (int i = 0; i < n; i++)
        {
          if (R[i] >= err) {idx = i; err = R[i];}
          if (R[i] >= 0) {mean_err += R[i];};
        }
        mean_err /= n;
        std::cout << "k: " << k << " error: " << err << " mean error: " << mean_err << "\n";
        if (err < tolerance*norm || k+1 >= max_degree) break;
      } // for k < n-1

      //
      for (int i = 0; i <= k; i++)
      {
        support_points.push_back(z[i]);
        function_values.push_back(f2[i]);
        weights.push_back(ww(i));
      }

      // Poles
      Matrix E(k+2,k+2), B(k+2,k+2);
      E.setZero(); B.setZero();
      for (int i = 0; i <= k; i++)
      {
        B(i,i) = 1.0;
        E(i+1,0) = ww(i);
        E(0,i+1) = 1.0;
        E(i+1,i+1) = z[i];
      }
      B(0,0) = 0.0;
      B(k+1,k+1) = 1.0;

      Eigen::VectorXcd D2;
      Eigen::MatrixXcd X;
      ZGGEV(E, B, D2, X);
      //std::cout << "poles D2: " << D2 << "\n";
      // Only keep non-infinity poles with positive(!?) imaginary component?
      for (int i = 0; i < k+2; i++)
      {
        if (std::abs(D2(i)) < mfem::infinity() && D2(i).imag() > 0.0) // do something better
        {
          poles.push_back(D2(i));
        }
      }
      // Sort poles by their imaginary component
      std::sort(poles.begin(), poles.end(),
              [](const std::complex<double>& a, const std::complex<double>& b) {
                  return a.imag() < b.imag();
              });

      std::cout << "sorted poles: " << "\n";
      for (int i = 0; i < poles.size(); i++)
      {
        std::cout << poles[i] << "\n";
      }

      // Return results
      AAAResult result;
      result.support_points = support_points;
      result.weights = weights;
      result.function_values = function_values;
      result.errors = R;
      result.degree = support_points.size();
      result.poles = poles;

      return result;
    }

private:

};

// this was written by AmazonQ and doesn't seem to work great...
std::vector<Complex> barycentricToMonomial(
    const std::vector<Complex>& support_points,  // zi
    const std::vector<Complex>& weights,         // wi
    const std::vector<Complex>& function_values  // fi
) {
    int n = support_points.size();

    // Initialize coefficient array
    std::vector<Complex> coeffs(n, Complex(0.0));

    // Convert using Lagrange interpolation approach
    for (int i = 0; i < n; ++i) {
        // Compute Lagrange basis polynomial for point i
        std::vector<Complex> basis_poly = {Complex(1.0)};

        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Multiply by (z - z_j) / (z_i - z_j)
                std::vector<Complex> new_poly(basis_poly.size() + 1, Complex(0.0));
                Complex factor = 1.0 / (support_points[i] - support_points[j]);

                for (size_t k = 0; k < basis_poly.size(); ++k) {
                    new_poly[k] += basis_poly[k] * (-support_points[j]) * factor;
                    new_poly[k+1] += basis_poly[k] * factor;
                }
                basis_poly = new_poly;
            }
        }

        // Add contribution to final polynomial
        for (size_t k = 0; k < basis_poly.size() && k < coeffs.size(); ++k) {
            coeffs[k] += function_values[i] * basis_poly[k];
        }
    }
    return coeffs;
}

// this is my attempt following https://people.inf.ethz.ch/gander/papers/changing.pdf ??
std::vector<Complex> barycentricToMonomial2(
    const std::vector<Complex>& support_points,  // zi
    const std::vector<Complex>& weights,         // wi
    const std::vector<Complex>& function_values  // fi
) {
    int n = support_points.size();

    // Initialize coefficient array
    std::vector<Complex> coeffs(n, Complex(0.0));

    // Monomial polynomial
    // V a = f where V is Vandermonde matrix
    // P = a^T m(x)  with m(x) = (1, x, x^2, ..., x^n)

    // Lagrange polynomial
    // P = f^T l(x)

    // First convert from Barycentric form to Lagrange polynomial
    // mu_i(z) = wi / (z-zi)
    // Pn(z) = sum(mu_i(z) * fi) / sum(mu_i(z))
    // li(z) = mu_i(z) / sum(mu_i(z))

    return coeffs;
}

std::vector<Complex> polynomialInterpolation(
 const std::vector<Complex>& support_points,  // zi
 const std::vector<Complex>& function_values  // fi
)
{
  int n = support_points.size();

  // Initialize coefficient array
  std::vector<Complex> coeffs(n, Complex(0.0));

  // Vandermonde matrix
  Matrix V(n,n);
  VectorEig f(n), a(n);
  for (int i = 0; i < n; i++)
  {
    f(i) = function_values[i];
    for (int j = 0; j < n; j++)
    {
      V(i,j) = std::pow(support_points[i], j);
    }
  }

  a = V.inverse() * f;
  std::cout << "monomial coefficients: " << a << "\n";
  for (int i = 0; i < n; i++)
  {
    coeffs[i] = a(i);
  }
  return coeffs;
}


Complex evalMonomial(const std::vector<Complex>& coeffs, const Complex &z)
{
  int n = coeffs.size();
  Complex result = 0.0;
  for (int i = 0; i < n; i++)
  {
    result += coeffs[i] * std::pow(z, i);
  }
  return result;
}

std::pair<ErrorIndicator, long long int>
EigenSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  // Construct and extract the system matrices defining the eigenvalue problem. The diagonal
  // values for the mass matrix PEC dof shift the Dirichlet eigenvalues out of the
  // computational range. The damping matrix may be nullptr.
  BlockTimer bt0(Timer::CONSTRUCT);
  SpaceOperator space_op(iodata, mesh);
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ONE);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  const auto &Curl = space_op.GetCurlMatrix();
  SaveMetadata(space_op.GetNDSpaces());

  // Configure objects for postprocessing.
  PostOperator<config::ProblemData::Type::EIGENMODE> post_op(iodata, space_op);
  ComplexVector E(Curl.Width()), B(Curl.Height());
  E.UseDevice(true);
  B.UseDevice(true);
  bool nonlinear = std::getenv("NONLINEAR_SLEPC"); // SHOULD DETECT BASED ON CONFIG!
  Mpi::Print("nonlinear: {:d}\n", nonlinear);

  // Check if nonlinear
  const double target = iodata.solver.eigenmode.target; // moved from below
  auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(target, Operator::DIAG_ZERO);
  bool has_A2 = (A2 != nullptr);
  //std::unique_ptr<ComplexOperator> Knew, Cnew, Mnew;
  std::vector<std::complex<double>> xs;
  std::vector<std::vector<std::unique_ptr<ComplexOperator>>> D_j;
  std::unique_ptr<Interpolation> newton_interp;
  if (has_A2)
  {
    Mpi::Print("Test RationalA2 inside eigensolver.cpp\n");
     const int npoints = 3;
    const double target_max = 2.0 * target; // get 3.0 from config file!

    newton_interp = std::make_unique<NewtonInterpolation>(space_op);
    newton_interp->Interpolate(npoints-1, 1i * target, 1i * target_max);

    // Compute rational interpolation and add matrices to K, C, M???
    // Currently copied code here but maybe should live elsewhere?

    xs.resize(npoints);
    for (int j = 0; j < npoints; j++)
    {
      xs[j] = std::complex<double>(0.0, target + j * (target_max - target) / (npoints - 1));
    }
    // Divided difference matrices (order 0 -> A2, order > 0 divided differences)
    D_j.resize(npoints);
    for (int k = 0; k < npoints; k++) // Order
    {
      for (int j = 0; j < npoints - k; j++)
      {
        if (k == 0)
        {
          auto A2j = space_op.GetExtraSystemMatrix<ComplexOperator>(xs[j].imag(), Operator::DIAG_ZERO);
          D_j[k].push_back(std::move(A2j));
        }
        else
        {
          //std::complex<double> denom = (xs[j+k] - xs[j]).imag();
          std::complex<double> denom = (xs[j+k] - xs[j]);
          auto A2dd = space_op.GetDividedDifferenceMatrix<ComplexOperator>(denom, D_j[k-1][j+1].get(), D_j[k-1][j].get(), Operator::DIAG_ZERO);
          D_j[k].push_back(std::move(A2dd));
        }
      }
    }
    //std::vector<std::complex<double>> coeffs;
    //std::vector<std::unique_ptr<ComplexOperator>> ops;
    //coeffs = {std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), - xs[0], xs[0] * xs[1]};
    //ops.push_back(std::move(K));
    //ops.push_back(std::move(D_j[0][0]));
    //ops.push_back(std::move(D_j[1][0]));
    //ops.push_back(std::move(D_j[2][0]));
    //Knew = space_op.GetExtraSystemMatrixSum2(coeffs, ops, Operator::DIAG_ONE);
    //coeffs = {std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), - (xs[0] + xs[1])};
    //ops[0] = std::move(C);
    //ops.erase(ops.begin() + 1);
    //Cnew = space_op.GetExtraSystemMatrixSum2(coeffs, ops, Operator::DIAG_ZERO);
    //coeffs = {std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0)};
    //ops[0] = std::move(M);
    //ops.erase(ops.begin() + 1);
    //Mnew = space_op.GetExtraSystemMatrixSum2(coeffs, ops, Operator::DIAG_ZERO);
    // Knew = K + Dj[0][0] - xs[0] * Dj[1][0] + xs[0]*xs[1] * Dj[2][0]
    // Cnew = C + Dj[1][0] - (xs[0] + xs[1]) * Dj[2][0]
    // Mnew = M + Dj[2][0]

    // Test AAA
    ComplexVector u; u.SetSize(K->Height()); u.UseDevice(true);
    ComplexVector v; v.SetSize(K->Height()); v.UseDevice(true);
    linalg::SetRandom(space_op.GetComm(), u);
    linalg::SetRandom(space_op.GetComm(), v);
    linalg::Normalize(space_op.GetComm(), u);
    linalg::Normalize(space_op.GetComm(), v);
    // Define your black box function
    auto f = [&](Complex z) -> Complex
    {
      auto A2z = space_op.GetExtraSystemMatrix<ComplexOperator>(std::abs(z.imag()), Operator::DIAG_ZERO);
      //auto Az = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), z, z * z, K.get(), C.get(), M.get(), A2z.get());
      ComplexVector Av; Av.SetSize(K->Height()); Av.UseDevice(true);
      //Az->Mult(v, Av);
      A2z->Mult(v, Av);
      return linalg::Dot(space_op.GetComm(), u, Av);
    };

    // Create sample points (you can customize this)
    std::vector<Complex> sample_points;
    int n_samples = 100;
    double y_max = (target_max - target) / 2.0;
    double y_min = - y_max;
    for (int i = 0; i < n_samples; ++i)
    {
        double x = target + i * (target_max - target) / (n_samples - 1);
        sample_points.push_back(Complex(0.0, x));
        //sample_points.push_back(Complex(y_min, x));
        //sample_points.push_back(Complex(y_max, x));
        //double y = y_min*1.01 + i * (y_max*0.99 - y_min*1.01) / (n_samples - 1);
        //sample_points.push_back(Complex(y, target));
        //sample_points.push_back(Complex(y, target_max));
    }
    /*
    // Apply AAA algorithm
    AAAApproximation aaa;
    auto result = aaa.computeAAA(f, sample_points, 1e-8, n_samples);

    // Print results
    std::cout << "AAA Approximation Results:" << std::endl;
    std::cout << "Degree: " << result.degree << std::endl;
    std::cout << "Support points:" << std::endl;
    for (size_t i = 0; i < result.support_points.size(); ++i) {
        std::cout << "z[" << i << "] = " << result.support_points[i]
                  << ", w[" << i << "] = " << result.weights[i]
                  << ", f[" << i << "] = " << result.function_values[i] << std::endl;
    }
    std::cout << "z_min: " << target << " zmax: " << target_max << "\n";

    //std::vector<Complex> coeffs1 = barycentricToMonomial(result.support_points, result.weights, result.function_values);
    //std::vector<Complex> coeffs2 = polynomialInterpolation(result.support_points, result.function_values);

    // Test the approximation
    Complex test_point(0.0, 0.5*(target+target_max));
    Complex exact_value = f(test_point);
    Complex approx_value = aaa.evaluate(test_point);
    //Complex approx_value1 = evalMonomial(coeffs1, test_point);
    //Complex approx_value2 = evalMonomial(coeffs2, test_point);

    std::cout << "\nTest at z = " << test_point << std::endl;
    std::cout << "Exact value: " << exact_value << std::endl;
    std::cout << "Approximation: " << approx_value << std::endl;
    //std::cout << "Approximation1: " << approx_value1 << std::endl;
    //std::cout << "Approximation2: " << approx_value2 << std::endl;
    std::cout << "Error: " << std::abs(exact_value - approx_value) << std::endl;

    std::cout << "Num poles: " << result.poles.size();

    // Use singularities and sample points to identify Leja-bagby points
    int max_it = 10;//???
    int n_poles = result.poles.size();
    std::vector<Complex> s(n_samples), nrs(n_samples), nrxi(n_samples), beta(n_samples); //?
    std::vector<Complex> xi(n_poles);

    s[0] = sample_points[0];
    xi[0] = result.poles[0];
    if (std::abs(xi[0]) < 10 * std::numeric_limits<double>::epsilon())
    {
      std::cout << "\n\n\n ERROR SINGULARITY CLOSE TO 0\n\n\n";
      exit(0);
    }
    beta[0] = 1.0;
    std::cout << "Leja-Bagby point " << 0 << " s: " << s[0] << " xi: " << xi[0] << " beta: " << beta[0] << "\n";
    for (int i = 0; i < n_samples; i++)
    {
      nrs[i] = 1.0;
      nrxi[i] = 1.0;
    }
    for (int k = 1; k < max_it; k++)
    {
      double maxnrs = 0.0;
      double minnrxi = mfem::infinity();
      for (int i = 0; i < n_samples; i++)
      {
        nrs[i] *= ((sample_points[i] - s[k-1])/(1.0 - sample_points[i]/xi[k-1])) / beta[k-1];
        if (std::abs(nrs[i]) > maxnrs)
        {
          maxnrs = std::abs(nrs[i]);
          s[k] = sample_points[i];
        }
      }
      if (n_poles > k)
      {
        for (int i = 1; i < n_poles; i++)
        {
          nrxi[i] *= ((result.poles[i] - s[k-1])/(1.0 - result.poles[i]/xi[k-1])) / beta[k-1];
          if (std::abs(nrxi[i]) < minnrxi)
          {
            minnrxi = std::abs(nrxi[i]);
            xi[k] = result.poles[i];
          }
        }
      }
      else
      {
        xi[k] = mfem::infinity();
      }
      beta[k] = maxnrs;
      std::cout << "Leja-Bagby point " << k << " s: " << s[k] << " xi: " << xi[k] << " beta: " << beta[k] << "\n";
    }

    //exit(0);
    */
  }

  // Define and configure the eigensolver to solve the eigenvalue problem:
  //         (K + λ C + λ² M + A2(λ)) u = 0    or    K u = -λ² M u
  // with λ = iω. In general, the system matrices are complex and symmetric.
  std::unique_ptr<EigenvalueSolver> eigen;
  config::EigenSolverData::Type type = iodata.solver.eigenmode.type;
#if defined(PALACE_WITH_ARPACK) && defined(PALACE_WITH_SLEPC)
  if (type == config::EigenSolverData::Type::DEFAULT)
  {
    type = config::EigenSolverData::Type::SLEPC;
  }
#elif defined(PALACE_WITH_ARPACK)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::SLEPC)
  {
    Mpi::Warning("SLEPc eigensolver not available, using ARPACK!\n");
  }
  type = config::EigenSolverData::Type::ARPACK;
#elif defined(PALACE_WITH_SLEPC)
  if (iodata.solver.eigenmode.type == config::EigenSolverData::Type::ARPACK)
  {
    Mpi::Warning("ARPACK eigensolver not available, using SLEPc!\n");
  }
  type = config::EigenSolverData::Type::SLEPC;
#else
#error "Eigenmode solver requires building with ARPACK or SLEPc!"
#endif
  if (type == config::EigenSolverData::Type::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    Mpi::Print("\nConfiguring ARPACK eigenvalue solver:\n");
    if (C)
    {
      eigen = std::make_unique<arpack::ArpackPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
    else
    {
      eigen = std::make_unique<arpack::ArpackEPSSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
    }
#endif
  }
  else  // config::EigenSolverData::Type::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    Mpi::Print("\nConfiguring SLEPc eigenvalue solver:\n");
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc;
    if (nonlinear)
    {
      Mpi::Print("Using SLEPc NEP solver\n");
      slepc = std::make_unique<slepc::SlepcNEPSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      //slepc->SetType(slepc::SlepcEigenvalueSolver::Type::NLEIGS);
      // slepc->SetType(slepc::SlepcEigenvalueSolver::Type::INTERPOL);  //only works with split operators (no callbacks)
      // slepc->SetType(slepc::SlepcEigenvalueSolver::Type::CISS); //only supports computing all Eigs
      //slepc->SetType(slepc::SlepcEigenvalueSolver::Type::RII);  //requires Jacobian and TARGET_MAGNITUDE.
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);  //requires Jacobian and TARGET_MAGNITUDE. Works-ish when updating pc shell, but kinda slow...
      // slepc->SetType(slepc::SlepcEigenvalueSolver::Type::NARNOLDI); //only works with split operators (no callbacks)
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
      //slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::RATIONAL);//test
    }
    else if (C)
    {
      if (!iodata.solver.eigenmode.pep_linear)
      {
        Mpi::Print("Using SLEPc PEP solver\n");
        slepc = std::make_unique<slepc::SlepcPEPSolver>(space_op.GetComm(),
                                                        iodata.problem.verbose);
        //slepc->SetType(slepc::SlepcEigenvalueSolver::Type::TOAR); //gives many wrong eigenvalues in addition to correct ones
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::LINEAR); // seems to work fine
        //slepc->SetType(slepc::SlepcEigenvalueSolver::Type::QARNOLDI); // seems to work fine, but only available for monomial?
        //slepc->SetType(slepc::SlepcEigenvalueSolver::Type::JD); // only works with precond
      }
      else
      {
        Mpi::Print("Using SLEPc PEPLinear solver\n");
        slepc = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(),
                                                              iodata.problem.verbose);
        slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      }
    }
    else
    {
      Mpi::Print("Using SLEPc EPS solver\n");
      slepc = std::make_unique<slepc::SlepcEPSSolver>(space_op.GetComm(),
                                                      iodata.problem.verbose);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    }
    if (!nonlinear) // this is ugly, need to handle better!
    {
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetOrthogonalization(
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::MGS,
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::CGS2);
    }
    eigen = std::move(slepc);
#endif
  }
  EigenvalueSolver::ScaleType scale = iodata.solver.eigenmode.scale
                                          ? EigenvalueSolver::ScaleType::NORM_2
                                          : EigenvalueSolver::ScaleType::NONE;
  if (nonlinear || C) //test force this
  {
    Mpi::Print("SetOperators with space_op, K, C, M\n");
    eigen->SetOperators(space_op, *K, *C, *M, scale);
    if (has_A2)
    {
      Mpi::Print("SetLinearA2Operators\n");
      eigen->SetLinearA2Operators(*D_j[0][0], *D_j[1][0], *D_j[2][0]);
      eigen->SetA2Interpolation(*newton_interp);
      Mpi::Print("Done SetLinearA2Operators\n");
    }
    // set NLEIGS numdegrees?
  }
  else if (C) // else if Cnew
  {
    eigen->SetOperators(*K, *C, *M, scale);
  }
  else
  {
    Mpi::Print("SetOperators with space_op, K, M\n");
    eigen->SetOperators(space_op, *K, *M, scale);
    //eigen->SetOperators(*K, *M, scale);
  }
  eigen->SetNumModes(iodata.solver.eigenmode.n, iodata.solver.eigenmode.max_size);
  eigen->SetTol(iodata.solver.eigenmode.tol);
  eigen->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen->GetScalingGamma(),
             eigen->GetScalingDelta());

  // If desired, use an M-inner product for orthogonalizing the eigenvalue subspace. The
  // constructed matrix just references the real SPD part of the mass matrix (no copy is
  // performed). Boundary conditions don't need to be eliminated here.
  std::unique_ptr<Operator> KM;
  if (iodata.solver.eigenmode.mass_orthog)
  {
    Mpi::Print(" Basis uses M-inner product\n");
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);

    // Mpi::Print(" Basis uses (K + M)-inner product\n");
    // KM = space_op.GetInnerProductMatrix(1.0, 1.0, K.get(), M.get());
    // eigen->SetBMat(*KM);
  }

  // Construct a divergence-free projector so the eigenvalue solve is performed in the space
  // orthogonal to the zero eigenvalues of the stiffness matrix.
  std::unique_ptr<DivFreeSolver<ComplexVector>> divfree;
  if (iodata.solver.linear.divfree_max_it > 0 &&
      !space_op.GetMaterialOp().HasWaveVector() &&
      !space_op.GetMaterialOp().HasLondonDepth())
  {
    Mpi::Print(" Configuring divergence-free projection\n");
    constexpr int divfree_verbose = 0;
    divfree = std::make_unique<DivFreeSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetH1Spaces(),
        space_op.GetAuxBdrTDofLists(), iodata.solver.linear.divfree_tol,
        iodata.solver.linear.divfree_max_it, divfree_verbose);
    eigen->SetDivFreeProjector(*divfree);
  }

  // If using Floquet BCs, a correction term (kp x E) needs to be added to the B field.
  std::unique_ptr<FloquetCorrSolver<ComplexVector>> floquet_corr;
  if (space_op.GetMaterialOp().HasWaveVector())
  {
    floquet_corr = std::make_unique<FloquetCorrSolver<ComplexVector>>(
        space_op.GetMaterialOp(), space_op.GetNDSpace(), space_op.GetRTSpace(),
        iodata.solver.linear.tol, iodata.solver.linear.max_it, 0);
  }

  // Set up the initial space for the eigenvalue solve. Satisfies boundary conditions and is
  // projected appropriately.
  if (iodata.solver.eigenmode.init_v0)
  {
    ComplexVector v0;
    if (iodata.solver.eigenmode.init_v0_const)
    {
      Mpi::Print(" Using constant starting vector\n");
      space_op.GetConstantInitialVector(v0);
    }
    else
    {
      Mpi::Print(" Using random starting vector\n");
      space_op.GetRandomInitialVector(v0);
    }
    if (divfree)
    {
      divfree->Mult(v0);
    }
    eigen->SetInitialSpace(v0);  // Copies the vector

    // Debug
    // const auto &Grad = space_op.GetGradMatrix();
    // ComplexVector r0(Grad->Width());
    // r0.UseDevice(true);
    // Grad.MultTranspose(v0.Real(), r0.Real());
    // Grad.MultTranspose(v0.Imag(), r0.Imag());
    // r0.Print();
  }

  // Configure the shift-and-invert strategy is employed to solve for the eigenvalues
  // closest to the specified target, σ.
  //const double target = iodata.solver.eigenmode.target;
  {
    const double f_target =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(target);
    Mpi::Print(" Shift-and-invert σ = {:.3e} GHz ({:.3e})\n", f_target, target);
  }
  const double l0 = target * 1.0;
  {
    const double f_l0 =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(l0);
    Mpi::Print(" Extra system matrix ;inearization point = {:.3e} GHz ({:.3e})\n", f_l0, l0);
  }
  if (C || nonlinear)
  {
    // Search for eigenvalues closest to λ = iσ.
    Mpi::Print("SetShiftInvert i*target and linearization point\n");
    //eigen->SetShiftInvert(1i * target, 1i * l0, 1i * target * 3.0);
    eigen->SetShiftInvert(1i * target);
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. The eigenvalue
      // 1 / (λ - σ) will be a large-magnitude negative imaginary number for an eigenvalue
      // λ with frequency close to but not below the target σ.
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_IMAGINARY);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY); // for linear or NLEIGS
      //eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE); // test for SLP/RII/NARNOLDI/?
      // Should we transform the problem so we can search for smallest_imaginary like with ARPACK?
    }
  }
  else
  {
    // Linear EVP has eigenvalues μ = -λ² = ω². Search for eigenvalues closest to μ = σ².
    //eigen->SetShiftInvert(target * target, l0 * l0, target * target * 3.0); // not sure l0 makes sense for linear EVP
    eigen->SetShiftInvert(target * target); // not sure l0 makes sense for linear EVP
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      // ARPACK searches based on eigenvalues of the transformed problem. 1 / (μ - σ²)
      // will be a large-magnitude positive real number for an eigenvalue μ with frequency
      // close to but below the target σ².
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_REAL);
    }
    else
    {
      eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_REAL);
    }
  }

  // Set up the linear solver required for solving systems involving the shifted operator
  // (K - σ² M) or P(iσ) = (K + iσ C - σ² M) during the eigenvalue solve. The
  // preconditioner for complex linear systems is constructed from a real approximation
  // to the complex system matrix.
  // linearize A2 around l0
  //const bool has_A2 = true;//false;
  //if (has_A2)
  //{
  const auto eps = std::sqrt(std::numeric_limits<double>::epsilon());
  //auto A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(l0, Operator::DIAG_ZERO); // or target??
  //bool has_A2 = (A2 != nullptr);
  if (has_A2)
  {
    auto A2p = space_op.GetExtraSystemMatrix<ComplexOperator>(l0 * (1.0 + eps), Operator::DIAG_ZERO);
    auto A2j = space_op.GetExtraSystemMatrixJacobian<ComplexOperator>(eps * l0, 1, A2p.get(), A2.get());
  }
  //}
  // Test to see waveport mode at difference frequencies
  /*
  std::vector<double> facs = {1.0, 1.2, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0, 96.0, 128.0};
  for (auto fac : facs)
  {
    auto A2test = space_op.GetExtraSystemMatrix<ComplexOperator>(target * fac, Operator::DIAG_ZERO);
  }
  */
  // TEST DIVIDED DIFFERENCE JACOBIAN -- ONLY WORKS IF A2 is not null
  /*
  std::vector<double> facs = {0.9999, 0.999, 0.99, 0.9, 0.5};
  for (auto c : facs)
  {
    double linearization_point = l0 * c;
    Mpi::Print("l0: {}, linearization point: {}\n", l0, linearization_point);
    auto A20 = space_op.GetExtraSystemMatrix<ComplexOperator>(linearization_point, Operator::DIAG_ZERO);
    auto A2p = space_op.GetExtraSystemMatrix<ComplexOperator>(linearization_point * (1.0 + eps), Operator::DIAG_ZERO);
    auto A2m = space_op.GetExtraSystemMatrix<ComplexOperator>(linearization_point * (1.0 - eps), Operator::DIAG_ZERO);
    auto A2jac = space_op.GetExtraSystemMatrixJacobian<ComplexOperator>(eps * linearization_point, 1, A2p.get(), A20.get(), A2m.get());
    auto A2jac2 = space_op.GetExtraSystemMatrixJacobian<ComplexOperator>(eps * linearization_point, 2, A2p.get(), A20.get(), A2m.get());
    ComplexVector tt(A2->Height());
    tt = std::complex<double>(1.23, 1.23);
    ComplexVector x1(A2->Height()), x2(A2->Height()), diff(A2->Height());
    x1 = 0.0; x2 = 0.0;
    A2->Mult(tt, x1);
    double normx1 = linalg::Norml2(space_op.GetComm(), x1);
    A20->Mult(tt, x2);
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    double res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 0 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    A2jac->AddMult(tt, x2, (l0 - linearization_point));
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 1 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    A2jac2->AddMult(tt, x2, 0.5 * pow(l0 - linearization_point, 2));
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 2 res: {}, res/normx1: {}, \n\n", res, res/normx1);
  }
  */
  /*
  // Test Rational/Newton Interpolation
  // Interpolation points: (or low-discrepancy sequence like Sobol?)
  double sigma0 = target, sigma1 = target*3.0, sigma2 = target*6.0;
  // Test points
  std::vector<double> test_p = {target*1.01, target*1.1, target*1.5, target*2, target*3.01, target*5, target*6.01, target*9};
  // fn(sigma) = T0 + (T1-T0)/(s1-s0)*(l-s0) + ((T2-T1)/(s2-s1) + (T1-T0)/(s1-s0))/(s2-s0)*(l-s0)*(l-s1)
  auto A20 = space_op.GetExtraSystemMatrix<ComplexOperator>(sigma0, Operator::DIAG_ZERO);
  auto A21 = space_op.GetExtraSystemMatrix<ComplexOperator>(sigma1, Operator::DIAG_ZERO);
  auto A22 = space_op.GetExtraSystemMatrix<ComplexOperator>(sigma2, Operator::DIAG_ZERO);
  auto A2d10 = space_op.GetExtraSystemMatrixJacobian<ComplexOperator>((sigma1-sigma0), 1, A21.get(), A20.get());
  auto A2d21 = space_op.GetExtraSystemMatrixJacobian<ComplexOperator>((sigma2-sigma1), 1, A22.get(), A21.get());
  for (auto l : test_p)
  {
    auto A2_l = space_op.GetExtraSystemMatrix<ComplexOperator>(l, Operator::DIAG_ZERO);
    ComplexVector tt(A2->Height());
    tt = std::complex<double>(1.23, 1.23);
    ComplexVector x1(A2->Height()), x2(A2->Height()), diff(A2->Height());
    x1 = 0.0; x2 = 0.0;
    A2_l->Mult(tt, x1);
    double normx1 = linalg::Norml2(space_op.GetComm(), x1);
    Mpi::Print("test lambda: {}, with interpolation points {}, {}, and {}\n", l, sigma0, sigma1, sigma2);
    // Order 0
    A20->Mult(tt, x2);
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    double res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 0 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    // Order 1
    A2d10->AddMult(tt, x2, (l - sigma0));
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 1 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    // Order 2
    A2d21->AddMult(tt, x2, (l - sigma0) * (l - sigma1) / (sigma2-sigma0));
    A2d10->AddMult(tt, x2, -(l - sigma0) * (l - sigma1) / (sigma2-sigma0));
    diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    res = linalg::Norml2(space_op.GetComm(), diff);
    Mpi::Print("Order 2 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    //// Order 1
    //A21->AddMult(tt, x2, (l - sigma0) / (sigma1 - sigma0));
    //A20->AddMult(tt, x2, -(l - sigma0) / (sigma1 - sigma0));
    //diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    //res = linalg::Norml2(space_op.GetComm(), diff);
    //Mpi::Print("Order 1 res: {}, res/normx1: {}, \n\n", res, res/normx1);
    // Order 2
    //A22->AddMult(tt, x2, (l - sigma0) * (l - sigma1) / ((sigma2-sigma1)*(sigma2-sigma0)));
    //A21->AddMult(tt, x2, -(l - sigma0) * (l - sigma1) / ((sigma2-sigma1)*(sigma2-sigma0)));
    //A21->AddMult(tt, x2, -(l - sigma0) * (l - sigma1) / ((sigma1-sigma0)*(sigma2-sigma0)));
    //A20->AddMult(tt, x2, (l - sigma0) * (l - sigma1) / ((sigma1-sigma0)*(sigma2-sigma0)));
    //diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
    //res = linalg::Norml2(space_op.GetComm(), diff);
    //Mpi::Print("Order 2 res: {}, res/normx1: {}, \n\n", res, res/normx1);
  }
  std::exit(0);
  */

  /*
  // Further rational interpolation tests
    Mpi::Print("\n Test Newton polynomial interpolation\n");
  double target_min = target, target_max = 5.0*target;
  int max_npoints = 3;//6;//6;//10;
  std::vector<double> res_vec(max_npoints, 0.0);
  for (int npoints = max_npoints; npoints <= max_npoints; npoints++)
  {
  // interpolation points
  std::vector<std::complex<double>> xs(npoints);
  //for (int j = 0; j < npoints; j+=2)
  for (int j = 0; j < npoints; j++)
  {
    xs[j] = std::complex<double>(0.0, target_min + j * (target_max - target_min) / (npoints - 1));
    //xs[j] = std::complex<double>(-target_min, target_min + j * (target_max - target_min) / (npoints - 1));
    //xs[j+1] = std::complex<double>(target_min, target_min + j * (target_max - target_min) / (npoints - 1));
    Mpi::Print("Interpolation point: {}+{}i\n", xs[j].real(), xs[j].imag());
    //Mpi::Print("Interpolation point: {}+{}i\n", xs[j+1].real(), xs[j+1].imag());
  }
  // Divided difference matrices (order 0 -> A2, order > 0 divided differences)
  std::vector<std::vector<std::unique_ptr<ComplexOperator>>> D_j(npoints);
  for (int k = 0; k < npoints; k++) // Order
  {
    for (int j = 0; j < npoints - k; j++)
    {
      if (k == 0)
      {
        auto A2i = space_op.GetExtraSystemMatrix<ComplexOperator>(xs[j].imag(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2i));
      }
      else
      {
        std::complex<double> denom = (xs[j+k]-xs[j]);
        auto A2dd = space_op.GetDividedDifferenceMatrix<ComplexOperator>(denom, D_j[k-1][j+1].get(), D_j[k-1][j].get(), Operator::DIAG_ZERO);
        D_j[k].push_back(std::move(A2dd));
      }
    }
  }
  Mpi::Print("\n Done building Dj operators, now trying test points\n");
  // Polynomial coefficients
  std::vector<std::vector<std::unique_ptr<ComplexOperator>>> C_j(npoints);
  // Initialize Cj
  // While this "works", the Cj operators get very expensive to use for high-order polynomials
  // probably because they are a sum of sum of operators
  // maybe we can find a way to express the coefficients as functions of the A2 matrices directly?
  //for (int i = 0; i < npoints; i++)
  //{
  //  C_j[i].push_back(std::move(D_j[i][0]));
  //}
  //for (int j = 0; j < npoints - 1; j++)
  //{
  //  for (int i = npoints - 2; i >= j; i--)
  //  {
  //    auto C_jnew = space_op.GetExtraSystemMatrixSum<ComplexOperator>(1.0, -xs[i-j], C_j[i][C_j[i].size()-1].get(), C_j[i+1][C_j[i+1].size()-1].get());
  //    C_j[i].push_back(std::move(C_jnew));//C_j[i] - xs[i-j] * C_j[i+1];
  //  }
  //}
  //Mpi::Print("\n Now building Cj operators\n");
  // test for n = 4
  // c0 = a0, c1 = a1, c2 = a2, c3 = a3
  // j = 0
  //   i = 2: c2 = a2 - x2 a3
  //   i = 1: c1 = a1 - x1 (a2 - x2 * a3) = a1 - a2 x1 + a3 x1 x2
  //   i = 0: c0 = a0 - x0 (a1 - a2 x1 + a3 x1 x2) = a0 - a1 x0 + a2 x0 x1 - a3 x0 x1 x2
  // j = 1
  //   i = 2: c2 = (a2 - x2 a3) - x1 a3
  //   i = 1: c1 = (a1 - a2 x1 + a3 x1 x2) - x0 (a2 - a3 x2 - a3 x1) = a1 - a2 x1 + a3 x1 x2 - a2 x0 + a3 x0 x2 + a3 x0 x1
  // j = 2
  //   i = 2: c2 = (a2 - a3 x2 - a3 x1) - x0 a3 = a2 - a3 x0 - a3 x1 - a3 x2

  //
  // Test points
  std::vector<std::complex<double>> test_p = {
    1i*target*1.0001,
    std::complex<double>(target, target*1.0001),
    1i*target*1.01,
    1i*target*1.1,
    1i*target*1.5,
    std::complex<double>(0.01*target, target*1.5),
    1i*target*2.0,
    1i*target*3.01,
    std::complex<double>(-0.2*target, target*3.01),
    std::complex<double>(0.2*target, target*3.01),
    std::complex<double>(0.5*target, target*3.01),
    std::complex<double>(1.0*target, target*3.01),
    std::complex<double>(-1.0*target, target*3.01),
    1i*target*5.0,
    std::complex<double>(0.3*target, target*5.0),
    std::complex<double>(target, target*5.0),
    std::complex<double>(-target, target*5.0)
  };
  double total_res = 0.0;
  for (auto l : test_p)
  {
    auto A2_l = space_op.GetExtraSystemMatrix<ComplexOperator>(l.imag(), Operator::DIAG_ZERO);
    ComplexVector tt(A2->Height());
    tt = std::complex<double>(1.23, 1.23);
    ComplexVector x1(A2->Height()), x2(A2->Height()), diff(A2->Height());
    x1 = 0.0; x2 = 0.0;
    A2_l->Mult(tt, x1);
    double normx1 = linalg::Norml2(space_op.GetComm(), x1);
    Mpi::Print("test lambda: {}+{}i, with {} interpolation points\n", l.real(), l.imag(), npoints);
    //for (int k = 0; k < npoints; k++)
    int k = npoints-1;
    {
      x2 = 0.0;
      std::complex<double> coeff = 1.0;
      for(int j = 0; j <= k; j++)
      {
        D_j[j][0]->AddMult(tt, x2, coeff);
        coeff *= (l - xs[j]);
        //Mpi::Print("adding order {} polynomial term with coeff: {}\n", j, coeff);
        //C_j[j][C_j[j].size()-1]->AddMult(tt, x2, coeff);
        //coeff *= l;
      }
      diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
      double res = linalg::Norml2(space_op.GetComm(), diff);
      Mpi::Print("Order {} res: {}, res/normx1: {}, \n\n", k, res, res/normx1);
      total_res += res/normx1;
    }
  }
  Mpi::Print("Order {}, total_res: {}\n", npoints-1, total_res / test_p.size());
  res_vec[npoints] = total_res / test_p.size();
  }
  //for (int i = 2; i <= max_npoints; i++)
 // {
 //   Mpi::Print("Order {}, total_res: {}\n", i-1, res_vec[i]);
 // }
  std::exit(0);
  */

  /*
  // Test Chebyshev interpolation
  Mpi::Print("\n\n Now testing Chebyshev interpolation!\n\n");
  double target_min = target, target_max = 10.0*target;
  int d = 3;
  double a = target_min, b = target_max;
  std::vector<double> nodes(d+1);
  for (int i = 0; i <= d; i++)
  {
    nodes[i] = 0.5 * (a + b) + 0.5 * (b - a) * std::cos((i + 0.5) * M_PI / (d + 1.0));
    Mpi::Print("Chebyshev node[{}]: {}\n", i, nodes[i]);
  }
  std::vector<std::unique_ptr<ComplexOperator>> A2_i(d+1);
  for (std::size_t i = 0; i <= d; i++)
  {
    auto A2i = space_op.GetExtraSystemMatrix<ComplexOperator>(nodes[i], Operator::DIAG_ZERO);
    A2_i[i] = std::move(A2i);
  }
  std::vector<std::vector<double>> P_coeffs(d+1);
  //std::vector<std::vector<std::complex<double>>> P_coeffs(d+1);
  for (std::size_t i = 0; i <= d; i++)
  {
    const double p = (i > 0) ? 2.0 : 1.0;
    for (std::size_t j = 0; j <= d; j++)
    {
      const double xj = (2.0 * (nodes[j] - a) / (b - a)) - 1.0;
      const double coeff = p / (d + 1.0) * std::cos(i * std::acos(xj));
      //P_coeffs[i].push_back(std::complex<double>(coeff, 0.0));
      P_coeffs[i].push_back(coeff);
      Mpi::Print("Polynomial coeff {}, {}: {}\n", i, j, coeff);
      //P[i] += (p / (d + 1.0) * std::cos(i * std::acos(xj))) * T[j];
    }
  }
  std::vector<std::unique_ptr<ComplexOperator>> Pi(d+1);
  for (int i = 0; i <= d; i++)
  {
    auto cheb = space_op.GetExtraSystemMatrixSum2(P_coeffs[i], A2_i, Operator::DIAG_ZERO);
    Pi[i] = std::move(cheb);
  }

  // T0 = 1, T1 = x, T2 = 2x^2 - 1, Tn+1 = 2xTn - Tn-1
  // T3 = 2xT2 - T1 = 2x(2x^2-1) - x = 4x^3-2x
  // c0 * T0 + c1 * T1 + c2 * T2 = a0 * 1 + a1 *

  std::vector<double> test_p = {target*1.01, target*1.1, target*1.5, target*2, target*3.01, target*5, target*6.01, target*9};
  double total_res = 0.0;
  for (auto l : test_p)
  {
    auto A2_l = space_op.GetExtraSystemMatrix<ComplexOperator>(l, Operator::DIAG_ZERO);
    ComplexVector tt(A2->Height());
    tt = std::complex<double>(1.23, 1.23);
    ComplexVector x1(A2->Height()), x2(A2->Height()), diff(A2->Height());
    x1 = 0.0; x2 = 0.0;
    A2_l->Mult(tt, x1);
    double normx1 = linalg::Norml2(space_op.GetComm(), x1);
    const double xj = (2.0 * (l - a) / (b - a)) - 1.0;
    Mpi::Print("test lambda: {} ({} in [-1,1]), with {} interpolation points\n", l, xj, d+1);

    for (int k = 0; k <= d; k++)
    {
      x2 = 0.0;
      for(int j = 0; j <= k; j++)
      {
        // Need to define T(x) for arbitrary order!
        Pi[j]->AddMult(tt, x2, std::cos(j * std::acos(xj)));
      }
      diff = 0.0; linalg::AXPBYPCZ(1.0, x1, -1.0, x2, 0.0, diff);
      double res = linalg::Norml2(space_op.GetComm(), diff);
      Mpi::Print("Order {} res: {}, res/normx1: {}, \n\n", k, res, res/normx1);
      total_res += res/normx1;
    }
  }
  Mpi::Print("Cheb Order {}, total_res: {}\n", d, total_res / test_p.size());
  std::exit(0);
  */

  Mpi::Print("Create A and P and call eigen->SetLinearSolver\n");
  // K + i sigma C + sigma^2 M + A2(sigma)
  auto A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * target,
                                    std::complex<double>(-target * target, 0.0), K.get(),
                                    C.get(), M.get(), A2.get());
  // K + i sigma C + sigma^2 M + A2(l0) + (sigma - l0) A2J(l0)
  //auto Atest = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * target,
  //                                  std::complex<double>(-target * target, 0.0), - 1i * l0, K.get(),
  //                                  C.get(), M.get(), A2.get(), A2j.get());
  //auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, target, -target * target, target);
  auto P = space_op.GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), 1i * target,
                                    std::complex<double>(-target * target, 0.0), target);
  auto ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(),
                                                &space_op.GetH1Spaces());
  //ksp->SetOperators(*Atest, *P);
  ksp->SetOperators(*A, *P);
  eigen->SetLinearSolver(*ksp);
  eigen->SetIoData(iodata);

  // Initialize structures for storing and reducing the results of error estimation.
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      space_op.GetMaterialOp(), space_op.GetNDSpaces(), space_op.GetRTSpaces(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);
  ErrorIndicator indicator;

  // Eigenvalue problem solve.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\n");
  int num_conv = eigen->Solve();
  {
    std::complex<double> lambda = (num_conv > 0) ? eigen->GetEigenvalue(0) : 0.0;
    Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv,
               (num_conv > 1) ? "s" : "",
               (num_conv > 0)
                   ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                   : "");
  }

  /*
  Mpi::Print("Setting up Quasi-Newton solver\n");
  std::unique_ptr<nleps::NonLinearEigenvalueSolver> qn;
  qn = std::make_unique<nleps::QuasiNewtonSolver>(space_op.GetComm(), iodata.problem.verbose);
  qn->SetTol(iodata.solver.eigenmode.tol);
  qn->SetMaxIter(iodata.solver.eigenmode.max_it);
  qn->SetOperators(space_op, *K, *C, *M, scale); // currently not using scaling but maybe try to make it work?
  qn->SetNumModes(num_conv, iodata.solver.eigenmode.max_size); // second input not actually used
  qn->SetLinearSolver(*ksp);
  qn->SetShiftInvert(1i * target);
  // Use linear eigensolve solution as initial guess.
  std::vector<std::complex<double>> init_eigs;
  std::vector<ComplexVector> init_V;
  for (int i = 0; i < num_conv; i++)
  {
    ComplexVector v0;
    v0.SetSize(Curl.Width());
    v0.UseDevice(true);
    eigen->GetEigenvector(i, v0);
    init_eigs.push_back(eigen->GetEigenvalue(i));
    init_V.push_back(v0);
  }
  qn->SetInitialGuess(init_eigs, init_V);
  eigen = std::move(qn); //?
  eigen->Solve();
  Mpi::Print("Done with QuasiNewton Solve\n");
  */

  BlockTimer bt2(Timer::POSTPRO);
  SaveMetadata(*ksp);

  // Calculate and record the error indicators, and postprocess the results.
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n");
  if (!KM)
  {
    // Normalize the finalized eigenvectors with respect to mass matrix (unit electric field
    // energy) even if they are not computed to be orthogonal with respect to it.
    KM = space_op.GetInnerProductMatrix(0.0, 1.0, nullptr, M.get());
    eigen->SetBMat(*KM);
    eigen->RescaleEigenvectors(num_conv);
  }
  Mpi::Print("\n");

  // Define an eigensolver for refinement?
  /*
  std::unique_ptr<EigenvalueSolver> eigen2;
  std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc2;
  Mpi::Print("Using SLEPc NEP solver to refine eigen solution\n");
  slepc2 = std::make_unique<slepc::SlepcNEPSolver>(space_op.GetComm(), iodata.problem.verbose);
  //slepc2->SetType(slepc::SlepcEigenvalueSolver::Type::RII);  //requires Jacobian and TARGET_MAGNITUDE
  slepc2->SetType(slepc::SlepcEigenvalueSolver::Type::SLP);  //requires Jacobian and TARGET_MAGNITUDE
  slepc2->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GENERAL);
  eigen2 = std::move(slepc2);
  eigen2->SetOperators(space_op, *K, *C, *M, scale);

  eigen2->SetNumModes(1, iodata.solver.eigenmode.max_size); // 1 or many??
  eigen2->SetTol(iodata.solver.eigenmode.tol);
  eigen2->SetMaxIter(iodata.solver.eigenmode.max_it);
  Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen2->GetScalingGamma(), eigen2->GetScalingDelta());

  eigen2->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_MAGNITUDE); // test for SLP/RII?
  eigen2->SetIoData(iodata);
    std::vector<ComplexVector> eigen_vectors;
  */


  for (int i = 0; i < num_conv; i++)
  {
    /*
    // Test refine eigen solution with RII or SLP?
    ComplexVector v0;
    v0.SetSize(Curl.Width());
    v0.UseDevice(true);
    eigen->GetEigenvector(i, v0);
    // Orthogonalize against previous eigenvectors?
    if (i > 0)
    {
      std::vector<std::complex<double>> Hj(i);
      linalg::OrthogonalizeColumnMGS(space_op.GetComm(), eigen_vectors, v0, Hj.data(), i);
      linalg::Normalize(space_op.GetComm(), v0);
    }
    eigen2->SetInitialSpace(v0);  // Copies the vector

    std::complex<double> omega_lin = eigen->GetEigenvalue(i);
    //eigen2->SetShiftInvert(omega_lin, omega_lin, 10.0 * omega_lin);
    eigen2->SetShiftInvert(omega_lin);

    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(std::abs(omega_lin.imag()), Operator::DIAG_ZERO);
    A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), omega_lin, omega_lin * omega_lin, K.get(), C.get(), M.get(), A2.get());
    P = space_op.GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), omega_lin, omega_lin * omega_lin, omega_lin.imag());

    ksp->SetOperators(*A, *P);
    eigen2->SetLinearSolver(*ksp);

    Mpi::Print("\n");
    Mpi::Print("Call eigen2->Solve() with initial guess: {:.3e}{:+.3e}i\n", omega_lin.real(), omega_lin.imag());
    int num_conv2 = eigen2->Solve();
    {
      std::complex<double> lambda = (num_conv2 > 0) ? eigen2->GetEigenvalue(0) : 0.0;
      Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv2,
                 (num_conv2 > 1) ? "s" : "",
                 (num_conv2 > 0)
                     ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                     : "");
    }
    */
    /*
    // Test refine by doing another linear eigensolve with an updated target??
    // Test refine eigen solution with RII or SLP?
    std::unique_ptr<EigenvalueSolver> eigen2;
    std::unique_ptr<slepc::SlepcEigenvalueSolver> slepc2;
    Mpi::Print("Using SLEPc NEP solver to refine eigen solution\n");
    slepc2 = std::make_unique<slepc::SlepcPEPLinearSolver>(space_op.GetComm(), iodata.problem.verbose);
    slepc2->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc2->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc2->SetOrthogonalization(
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::MGS,
        iodata.solver.linear.gs_orthog_type == config::LinearSolverData::OrthogType::CGS2);
    eigen2 = std::move(slepc2);
    eigen2->SetOperators(space_op, *K, *C, *M, scale);

    eigen2->SetNumModes(1, iodata.solver.eigenmode.max_size);
    eigen2->SetTol(iodata.solver.eigenmode.tol);
    //eigen2->SetTol(1e-4);
    eigen2->SetMaxIter(iodata.solver.eigenmode.max_it);
    Mpi::Print(" Scaling γ = {:.3e}, δ = {:.3e}\n", eigen2->GetScalingGamma(), eigen2->GetScalingDelta());

    //ComplexVector v0;
    //v0.SetSize(Curl.Width());
    //v0.UseDevice(true);
    //eigen->GetEigenvector(i, v0);
    //eigen2->SetInitialSpace(v0);  // Copies the vector

    std::complex<double> omega_lin = eigen->GetEigenvalue(i);
    double new_target = omega_lin.imag() - eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    eigen2->SetShiftInvert(1i * new_target);
    eigen2->SetWhichEigenpairs(EigenvalueSolver::WhichType::TARGET_IMAGINARY);

    // ??
    A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(std::abs(omega_lin.imag()), Operator::DIAG_ZERO);
    A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), omega_lin, omega_lin * omega_lin, K.get(), C.get(), M.get(), A2.get());
    //P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, omega_lin.imag(), - omega_lin.imag() * omega_lin.imag(), omega_lin.imag());
    P = space_op.GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), omega_lin, omega_lin * omega_lin, omega_lin.imag());
    //A2 = space_op.GetExtraSystemMatrix<ComplexOperator>(new_target, Operator::DIAG_ZERO);
    //A = space_op.GetSystemMatrix(std::complex<double>(1.0, 0.0), 1i * new_target, std::complex<double>(-new_target * -new_target, 0.0), K.get(), C.get(), M.get(), A2.get());
    //P = space_op.GetPreconditionerMatrix<ComplexOperator>(1.0, new_target, - new_target * new_target, new_target);
    //P = space_op.GetPreconditionerMatrix<ComplexOperator>(std::complex<double>(1.0, 0.0), 1i * new_target, std::complex<double>(-new_target * -new_target, 0.0), new_target);
    ksp = std::make_unique<ComplexKspSolver>(iodata, space_op.GetNDSpaces(), &space_op.GetH1Spaces());
    ksp->SetOperators(*A, *P);
    eigen2->SetLinearSolver(*ksp);

    Mpi::Print("\n");
    Mpi::Print("Call eigen2->Solve() with initial guess: {:.3e}{:+.3e}i\n", omega_lin.real(), omega_lin.imag());
    int num_conv2 = eigen2->Solve();
    {
      std::complex<double> lambda = (num_conv2 > 0) ? eigen2->GetEigenvalue(0) : 0.0;
      Mpi::Print(" Found {:d} converged eigenvalue{}{}\n", num_conv2,
                 (num_conv2 > 1) ? "s" : "",
                 (num_conv2 > 0)
                     ? fmt::format(" (first = {:.3e}{:+.3e}i)", lambda.real(), lambda.imag())
                     : "");
    }
    */

    // Get the eigenvalue and relative error.
    std::complex<double> omega = eigen->GetEigenvalue(i);
    //std::complex<double> omega = eigen2->GetEigenvalue(0);
    double error_bkwd = eigen->GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eigen->GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);
    //double error_bkwd = eigen2->GetError(0, EigenvalueSolver::ErrorType::BACKWARD);
    //double error_abs = eigen2->GetError(0, EigenvalueSolver::ErrorType::ABSOLUTE);

    // Could have something like
    // What's the best interface? currently it's in SLEPc but that's not good as it does not depend on SLEPc or ARPACK
    // so it should be independent
    // std::unique_ptr<> nleigen(space_op, iodata) // maybe don't need iodata
    // nleigen->SetLinearSolve(ksp);
    // nleigen->SetTol(iodata.solver.eigenmode.tol)
    // nleigen->SetOperators(K, C, M);
    // Initialize the nonlinear eigensolver above, kinda like floquet_corr, and call it here??
    // if (nonlinear && error_abs > iodata.solver.eigenmode.tol)
    //.  nleigen->Refine(omega, E);
    //   RII(space_op, K, C, M, ksp, tol, omega, E);
    //   QN2(space_op, K, C, M, ksp, tol, omega, E)
    //


    if (!C && !nonlinear)
    {
      // Linear EVP has eigenvalue μ = -λ² = ω².
      omega = std::sqrt(omega);
    }
    else
    {
      // Quadratic EVP solves for eigenvalue λ = iω.
      omega /= 1i;
    }

    // Compute B = -1/(iω) ∇ x E on the true dofs, and set the internal GridFunctions in
    // PostOperator for all postprocessing operations.
    eigen->GetEigenvector(i, E);
    //eigen2->GetEigenvector(0, E);
    //eigen_vectors.push_back(E);

    Curl.Mult(E.Real(), B.Real());
    Curl.Mult(E.Imag(), B.Imag());
    B *= -1.0 / (1i * omega);
    if (space_op.GetMaterialOp().HasWaveVector())
    {
      // Calculate B field correction for Floquet BCs.
      // B = -1/(iω) ∇ x E + 1/ω kp x E.
      floquet_corr->AddMult(E, B, 1.0 / omega);
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, E, B, omega, error_abs, error_bkwd, num_conv);

    // Calculate and record the error indicators.
    if (i < iodata.solver.eigenmode.n)
    {
      estimator.AddErrorIndicator(E, B, total_domain_energy, indicator);
    }

    // Final write: Different condition than end of loop (i = num_conv - 1).
    if (i == iodata.solver.eigenmode.n - 1)
    {
      post_op.MeasureFinalize(indicator);
    }
  }
  return {indicator, space_op.GlobalTrueVSize()};
}

}  // namespace palace
