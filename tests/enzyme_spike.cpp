// Enzyme/dual-number spike test: differentiate Palace mesh sensitivity integrands w.r.t.
// Jacobian J.
//
// This standalone test verifies that the mesh sensitivity QFunction header correctly
// differentiates the CurlCurlMass integrand. It tests both the Enzyme path (if compiled
// with the Enzyme plugin) and the dual-number path, comparing against central FD.
//
// The integrands match f_apply_hdivmass_33 in Palace:
//   Stiffness: f_K = w * detJ * curl_u^T * J_eff^T * mu_inv * J_eff * curl_u
//              where J_eff = adj(adjJt/detJ)^T = J/detJ  (H(div) Piola map)
//   Mass:      f_M = w * detJ * u^T * adjJt_norm^T * eps * adjJt_norm * u
//              where adjJt_norm = adj(J)^T / detJ         (H(curl) Piola map)
//
// Compile with Enzyme:
//   clang++ -std=c++17 -O2 -fplugin=/opt/homebrew/opt/enzyme/lib/ClangEnzyme-21.dylib \
//     -DPALACE_WITH_ENZYME -I ../palace/fem/qfunctions -o enzyme_spike enzyme_spike.cpp
//
// Compile without Enzyme (dual-number path only):
//   clang++ -std=c++17 -O2 -I ../palace/fem/qfunctions -o enzyme_spike enzyme_spike.cpp
//
// Run:
//   ./enzyme_spike

#include <cmath>
#include <cstdio>
#include <cstring>

// Provide the macros that Palace QFunction headers expect.
#define CeedScalar double
#define CeedInt int
#define CEED_QFUNCTION_HELPER static inline
#define CeedPragmaSIMD
#define CEED_QFUNCTION(name) static int name
#define CEED_EVAL_INTERP 0
#define CEED_EVAL_WEIGHT 0
#define CEED_EVAL_GRAD 0
#define CEED_EVAL_CURL 0

// Include the utility functions for FD reference.
#include "33/utils_33_qf.h"

// Include the mesh sensitivity QFunction header — this brings in the Dual type, templated
// integrands, and MeshSensDifferentiate.
#include "33/meshsens_33_qf.h"

// Helper: compute FD gradient of a scalar function f(J) w.r.t. J[9].
template <typename Func>
void fd_gradient(Func f, const double *J, double *dJ_fd, double h)
{
  for (int k = 0; k < 9; k++)
  {
    double J_plus[9], J_minus[9];
    memcpy(J_plus, J, sizeof(double) * 9);
    memcpy(J_minus, J, sizeof(double) * 9);
    J_plus[k] += h;
    J_minus[k] -= h;
    dJ_fd[k] = (f(J_plus) - f(J_minus)) / (2.0 * h);
  }
}

// Helper: compute max relative error between two 9-vectors.
double max_relative_error(const double *a, const double *b)
{
  double max_err = 0.0;
  for (int k = 0; k < 9; k++)
  {
    double ref = fabs(b[k]) > 1e-15 ? b[k] : 1.0;
    double err = fabs(a[k] - b[k]) / fabs(ref);
    if (err > max_err)
    {
      max_err = err;
    }
  }
  return max_err;
}

int main()
{
  // Test Jacobian: a slightly skewed element (not identity — more realistic).
  double J[9] = {1.1, 0.1, 0.0,   // column 0
                 0.05, 0.9, 0.1,   // column 1
                 0.0, -0.05, 1.2}; // column 2

  // Test curl of eigenvector at this quad point (arbitrary nonzero).
  const double curl_u[3] = {0.3, -0.5, 0.8};

  // Test eigenvector value at this quad point (for mass term).
  const double u[3] = {0.1, 0.7, -0.4};

  // Isotropic mu_inv = I (identity — free space).
  const double mu_inv[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  // Anisotropic eps tensor (to verify the full tensor path).
  const double eps[9] = {2.08, 0.0, 0.0, 0.0, 2.08, 0.0, 0.0, 0.0, 2.08};

  const double w = 0.25;       // quadrature weight
  const double lambda = 3.14;  // eigenvalue (omega^2)

  printf("=== Mesh Sensitivity Spike Test ===\n\n");
  printf("Differentiation method: "
#ifdef PALACE_WITH_ENZYME
         "Enzyme reverse-mode AD"
#else
         "Forward-mode dual numbers"
#endif
         "\n\n");

  // ================================================================
  // Test 1: Stiffness integrand (matching hdivmass_33 stiffness block)
  // ================================================================
  printf("--- Stiffness integrand (J/detJ Piola, matching hdivmass_33) ---\n");

  auto f_K = [&](const double *Jarg)
  { return MeshSensStiffnessIntegrand<CeedScalar>(Jarg, curl_u, mu_inv, w); };

  double val_K = f_K(J);
  printf("f_K = %.15e\n", val_K);

  double dK_fd[9];
  fd_gradient(f_K, J, dK_fd, 1e-7);

  // Use the dual-number path explicitly for stiffness.
  double dK_dual[9];
  for (int k = 0; k < 9; k++)
  {
    MeshSensDual J_d[9];
    for (int j = 0; j < 9; j++)
    {
      J_d[j] = MakeDual(J[j], (j == k) ? 1.0 : 0.0);
    }
    MeshSensDual f = MeshSensStiffnessIntegrand(J_d, curl_u, mu_inv, w);
    dK_dual[k] = f.dot;
  }

  double err_K = max_relative_error(dK_dual, dK_fd);
  printf("Dual vs FD max relative error: %.6e %s\n", err_K, err_K < 1e-5 ? "PASS" : "FAIL");

#ifdef PALACE_WITH_ENZYME
  // Also test Enzyme on the stiffness integrand.
  double dK_enzyme[9] = {0.0};
  double J_copy[9];
  memcpy(J_copy, J, sizeof(J));
  __enzyme_autodiff((void *)MeshSensStiffnessForEnzyme, enzyme_dup, J_copy, dK_enzyme,
                    enzyme_const, curl_u, enzyme_const, mu_inv, enzyme_const, w);
  double err_K_enzyme = max_relative_error(dK_enzyme, dK_fd);
  printf("Enzyme vs FD max relative error: %.6e %s\n", err_K_enzyme,
         err_K_enzyme < 1e-5 ? "PASS" : "FAIL");
#endif

  // ================================================================
  // Test 2: Mass integrand (matching hdivmass_33 mass block, H(curl) Piola)
  // ================================================================
  printf("\n--- Mass integrand (adjJt/detJ Piola, matching hdivmass_33) ---\n");

  auto f_M = [&](const double *Jarg)
  { return MeshSensMassIntegrand<CeedScalar>(Jarg, u, eps, w); };

  double val_M = f_M(J);
  printf("f_M = %.15e\n", val_M);

  double dM_fd[9];
  fd_gradient(f_M, J, dM_fd, 1e-7);

  double dM_dual[9];
  for (int k = 0; k < 9; k++)
  {
    MeshSensDual J_d[9];
    for (int j = 0; j < 9; j++)
    {
      J_d[j] = MakeDual(J[j], (j == k) ? 1.0 : 0.0);
    }
    MeshSensDual f = MeshSensMassIntegrand(J_d, u, eps, w);
    dM_dual[k] = f.dot;
  }

  double err_M = max_relative_error(dM_dual, dM_fd);
  printf("Dual vs FD max relative error: %.6e %s\n", err_M, err_M < 1e-5 ? "PASS" : "FAIL");

#ifdef PALACE_WITH_ENZYME
  double dM_enzyme[9] = {0.0};
  memcpy(J_copy, J, sizeof(J));
  __enzyme_autodiff((void *)MeshSensMassForEnzyme, enzyme_dup, J_copy, dM_enzyme,
                    enzyme_const, u, enzyme_const, eps, enzyme_const, w);
  double err_M_enzyme = max_relative_error(dM_enzyme, dM_fd);
  printf("Enzyme vs FD max relative error: %.6e %s\n", err_M_enzyme,
         err_M_enzyme < 1e-5 ? "PASS" : "FAIL");
#endif

  // ================================================================
  // Test 3: Combined integrand (stiffness - lambda * mass)
  // ================================================================
  printf("\n--- Combined integrand (f_K - lambda * f_M) ---\n");

  auto f_combined = [&](const double *Jarg)
  {
    return MeshSensStiffnessIntegrand<CeedScalar>(Jarg, curl_u, mu_inv, w) -
           lambda * MeshSensMassIntegrand<CeedScalar>(Jarg, u, eps, w);
  };

  double val_combined = f_combined(J);
  printf("f_combined = %.15e\n", val_combined);

  double df_fd[9];
  fd_gradient(f_combined, J, df_fd, 1e-7);

  // Use MeshSensDifferentiate (the unified wrapper).
  double df_wrapper[9];
  MeshSensDifferentiate(J, curl_u, mu_inv, u, eps, w, lambda, df_wrapper);

  double err_combined = max_relative_error(df_wrapper, df_fd);
  printf("MeshSensDifferentiate vs FD max relative error: %.6e %s\n", err_combined,
         err_combined < 1e-5 ? "PASS" : "FAIL");

  // ================================================================
  // Summary
  // ================================================================
  printf("\n=== Summary ===\n");
  printf("Stiffness (dual vs FD): %.6e %s\n", err_K, err_K < 1e-5 ? "PASS" : "FAIL");
  printf("Mass (dual vs FD):      %.6e %s\n", err_M, err_M < 1e-5 ? "PASS" : "FAIL");
  printf("Combined (wrapper vs FD): %.6e %s\n", err_combined,
         err_combined < 1e-5 ? "PASS" : "FAIL");
#ifdef PALACE_WITH_ENZYME
  printf("Stiffness (Enzyme vs FD): %.6e %s\n", err_K_enzyme,
         err_K_enzyme < 1e-5 ? "PASS" : "FAIL");
  printf("Mass (Enzyme vs FD):      %.6e %s\n", err_M_enzyme,
         err_M_enzyme < 1e-5 ? "PASS" : "FAIL");
#endif

  bool pass = (err_K < 1e-5) && (err_M < 1e-5) && (err_combined < 1e-5);
#ifdef PALACE_WITH_ENZYME
  pass = pass && (err_K_enzyme < 1e-5) && (err_M_enzyme < 1e-5);
#endif

  if (pass)
  {
    printf("\nAll tests PASSED.\n");
    return 0;
  }
  else
  {
    printf("\nSome tests FAILED.\n");
    return 1;
  }
}
