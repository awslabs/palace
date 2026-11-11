// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_DEFLATION_HPP
#define PALACE_LINALG_DEFLATION_HPP

#include <algorithm>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <mpi.h>
#include "linalg/operator.hpp"
#include "linalg/solver.hpp"
#include "linalg/vector.hpp"

namespace palace
{

//
// Harmonic Ritz analysis of the Arnoldi relation B A V[0:m] = V[0:m+1] H of a
// left-preconditioned GMRES solve with the operator A and the preconditioner B, for the
// deflation of later solves with the same operators. Computed on the root process and
// broadcast, so that every process selects the same vectors.
//
struct DeflationAnalysis
{
  // Number of Arnoldi steps of the analyzed solve.
  int m = 0;

  // Harmonic Ritz values of B A with respect to zero, by increasing magnitude, and the
  // coefficients in the basis V of the corresponding vectors (m x kmax, column j for the
  // j-th value).
  std::vector<std::complex<double>> theta;
  Eigen::MatrixXcd Y;

  // Iterations to reach the relative tolerance of GMRES started from the first basis vector
  // for the projected operator (model), plain and with the k harmonic Ritz values closest
  // to zero moved to one (deflated), for the candidate numbers of vectors.
  int its_plain = 0;
  std::vector<int> k_candidates, its_deflated;

  // Number of harmonic Ritz values with |θ| < theta_max (outliers of B A).
  int num_outliers = 0;

  // Compute the analysis for the at most kmax vectors closest to zero and the model
  // iteration counts for the candidates k_candidates which do not exceed the number of
  // values with |θ| < theta_max, for that number itself, and for the numbers of values
  // below 2 and 4 times theta_max (theta_max <= 0: all candidates up to kmax).
  void Compute(MPI_Comm comm, const std::complex<double> *H, int ld, int m, int kmax,
               const std::vector<int> &k_candidates, double rel_tol,
               double theta_max = 0.0);

  // Number of the analyzed harmonic Ritz values with |θ| < bound.
  int Count(double bound) const;
};

//
// Deflation of a left preconditioner B of the operator A with an orthonormal basis W of k
// vectors (A-DEF1 two-level preconditioner, see Tang et al., J. Sci. Comput. 39 (2009)):
//
//   P r = B (r - A W E⁻¹ Wᴴ r) + W E⁻¹ Wᴴ r,   E = Wᴴ A W,
//
// for which P A W = W: the directions of W, typically approximate eigenvectors of B A with
// eigenvalues close to zero which slow down the convergence of the Krylov solver, are
// solved exactly. (The Galerkin initial guess W E⁻¹ Wᴴ b does not reduce the number of
// iterations further, and costs GMRES one more preconditioner application per solve.) The
// basis holds a reserve of vectors, of which the first k are active (deflated): activating
// more vectors computes their image under A, while fewer active vectors only need the
// factorization of the leading block of E. On the host, W and A W are stored in single
// precision, which halves the memory traffic of the deflation (E and its factorization are
// computed from the rounded vectors, so that P A W = W holds to the rounding of A W).
//
class DeflatedPreconditioner : public Solver<ComplexOperator>
{
private:
  MPI_Comm comm;
  const ComplexOperator *A = nullptr;
  const Solver<ComplexOperator> *B = nullptr;

  // Deflation basis W and its image A W (for the vectors activated so far), in double
  // precision on devices and in single precision on the host (real and imaginary parts, n
  // entries each).
  bool single;
  int n = 0;
  std::vector<ComplexVector> W, AW;
  std::vector<std::vector<float>> W32, AW32;

  // Number of basis vectors from the first j analyzed harmonic Ritz vectors (without the
  // numerically dependent ones), for j up to the number of vectors of the basis build.
  std::vector<int> prefix;

  // Number of active vectors, E = Wᴴ A W for the vectors activated so far, and the
  // factorization of its leading block for the active ones.
  int active = 0;
  Eigen::MatrixXcd E;
  Eigen::PartialPivLU<Eigen::MatrixXcd> lu;

  // Workspace.
  mutable ComplexVector t;
  mutable std::vector<std::complex<double>> c;

  // Coefficients c = E⁻¹ Wᴴ x.
  void Coefficients(const ComplexVector &x) const;

  // y += W c.
  void AddBasis(const std::complex<double> *c, ComplexVector &y) const;

public:
  DeflatedPreconditioner(MPI_Comm comm);

  // Build the basis W from the first k vectors of the analysis of the Arnoldi relation with
  // the orthonormal basis V (k <= the number of analyzed vectors), with no active vectors.
  // Returns the number of basis vectors (fewer than k when the vectors are numerically
  // dependent).
  int Build(const ComplexOperator &A, const Solver<ComplexOperator> &B,
            const std::vector<ComplexVector> &V, const DeflationAnalysis &analysis, int k);

  // Activate the first k vectors of the basis (at most the number of basis vectors).
  // Returns the number of active vectors (k, or the previous number when E is singular for
  // k vectors).
  int Activate(int k);

  // Remove the deflation basis.
  void Clear();

  // Number of active (deflated) vectors.
  int Size() const { return active; }

  // Number of basis vectors, and the number of them from the first j analyzed harmonic Ritz
  // vectors.
  int Reserve() const
  {
    return single ? static_cast<int>(W32.size()) : static_cast<int>(W.size());
  }
  int Prefix(int j) const
  {
    return prefix.empty() ? 0 : prefix[std::min<std::size_t>(j, prefix.size() - 1)];
  }

  // Memory of the deflation basis and its image.
  std::size_t GetBytes() const;

  // The preconditioner which is deflated.
  void SetPreconditioner(const Solver<ComplexOperator> &pc) { B = &pc; }

  void SetOperator(const ComplexOperator &op) override {}

  void Mult(const ComplexVector &x, ComplexVector &y) const override;
};

// Time per basis vector of one application of the deflation (the inner products and the two
// vector updates, with the storage precision of DeflatedPreconditioner), measured on copies
// of up to 32 of the given vectors (more than fit in cache), maximum over the processes,
// and the time of the global reduction of the inner products (seconds).
std::pair<double, double> MeasureDeflationCost(MPI_Comm comm,
                                               const std::vector<ComplexVector> &V, int m);

//
// Automatic choice of the number of deflated vectors among nested sets of increasing size
// (the harmonic Ritz vectors with values below increasing bounds): a verification which
// finds that the current set does not pay goes back to the previous set or, from the first
// one, stops the deflation; while the deflation pays, the next set is tried, and kept if it
// pays more than the previous one. A set which did not pay is not tried again.
//
class DeflationLevels
{
public:
  enum class Action
  {
    KEEP,
    GROW,
    BACK,
    STOP
  };

private:
  std::vector<int> sizes;
  int level = 0, top = 0;
  bool trying = false;

  Action Grow();

public:
  // Sets of the given sizes (increasing; sizes which do not increase are skipped).
  void Reset(const std::vector<int> &sizes);

  // Number of vectors of the current set and of the previous one.
  int Size() const { return sizes.empty() ? 0 : sizes[level]; }
  int PreviousSize() const { return level > 0 ? sizes[level - 1] : 0; }

  // Whether the current set was just tried, and is to be compared with the previous one.
  bool Trying() const { return trying; }

  // After a verification which finds whether the current set pays, against no deflation.
  Action Verify(bool gain);

  // After the comparison of the set being tried with the previous one.
  Action Compare(bool gain);
};

}  // namespace palace

#endif  // PALACE_LINALG_DEFLATION_HPP
