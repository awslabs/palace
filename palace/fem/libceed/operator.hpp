// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_OPERATOR_HPP
#define PALACE_LIBCEED_OPERATOR_HPP

#include <complex>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;

namespace hypre
{

class HypreCSRMatrix;

}  // namespace hypre

namespace ceed
{

class ElementMatrixTerms;

// A linear combination sum_k coeff[k] T_k of the terms of an ElementMatrixTerms.
struct ElementMatrixCombination
{
  std::shared_ptr<const ElementMatrixTerms> terms;
  std::vector<std::complex<double>> coeff;
};

//
// Wrapper class for libCEED's CeedOperator, supporting composite operator construction and
// application with multiple threads.
//
class Operator : public palace::Operator
{
protected:
  std::vector<CeedOperator> op, op_t;
  std::vector<CeedVector> u, v;
  Vector dof_multiplicity;
  mutable Vector temp;

public:
  // Element matrices of leaves applied in place of their partial assembly in forward
  // applications on CPUs, see UseElementAssembly.
  struct ElementAssembly;
  const ElementAssembly *GetElementAssembly() const { return ea.get(); }

protected:
  std::unique_ptr<ElementAssembly> ea;
  bool ea_single = false;
  ElementMatrixCombination ea_comb;
  bool ea_comb_imag = false;
  void ElementAssemblyMult(const Vector &x, Vector &y, bool add) const;

public:
  Operator(int h, int w);
  ~Operator() override;

  CeedOperator operator[](std::size_t i) const { return op[i]; }

  auto Size() const { return op.size(); }

  void AddSubOperator(CeedOperator sub_op, CeedOperator sub_op_t = nullptr);

  void Finalize();

  void DestroyAssemblyData() const;

  void SetDofMultiplicity(Vector &&mult) { dof_multiplicity = std::move(mult); }

  // For forward applications on a CPU, apply the dense element matrices of small simplex
  // leaves instead of their partial assembly when that is cheaper for the given number of
  // right hand sides per application. Transposes, diagonal and full assembly keep using the
  // partially assembled leaves. Requires a finalized square operator with a single CEED
  // context. Returns whether any leaf is element assembled.
  bool UseElementAssembly(int num_rhs = 1);

  // Store the element matrices of UseElementAssembly in single precision (with the
  // arithmetic in double precision), for operators whose applications tolerate a relative
  // error of about 1e-7, such as the auxiliary space operators of the multigrid smoothers.
  // Call before UseElementAssembly.
  void SetElementAssemblySinglePrecision(bool single) { ea_single = single; }

  // Form the single-precision element matrices of UseElementAssembly from those of the
  // terms of the given combination, for an operator which is its real part (or imaginary
  // part, with imag) assembled over the restrictions of the terms (see
  // ElementMatrixTerms). Leaves the terms do not cover are assembled. Call before
  // UseElementAssembly.
  void SetElementAssemblyTerms(const ElementMatrixCombination &comb, bool imag)
  {
    ea_comb = comb;
    ea_comb_imag = imag;
  }

  bool HasDofMultiplicity() const { return dof_multiplicity.Size() > 0; }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override;

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override;

  void MultTranspose(const Vector &x, Vector &y) const override;

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override;
};

// A symmetric ceed::Operator replaces *MultTranspose with *Mult (by default, libCEED
// operators do not have a transpose operation).
class SymmetricOperator : public Operator
{
public:
  using Operator::Operator;

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
  void AddMultTranspose(const Vector &x, Vector &y, double a = 1.0) const override
  {
    AddMult(x, y, a);
  }
};

// Wrap finalized, owned real/imaginary operators, packing compatible QData-assembled
// volume terms on CPU. Other terms and unsupported operators keep their original action.
// Compatible leaves assembled over an element subset of the given finite element space (the
// one both operators were assembled from) are packed over that subset, with the elements
// they leave out applied separately. The inputs must remain structurally unchanged and
// unrescaled after ownership transfer, including through retained aliases or raw CEED
// handles. Passive QData values remain shared and may be updated. Borrowed operators should
// use ComplexWrapperOperator.
std::unique_ptr<ComplexWrapperOperator>
CreateComplexOperator(std::unique_ptr<palace::Operator> &&Ar,
                      std::unique_ptr<palace::Operator> &&Ai,
                      const FiniteElementSpace &fespace);

// Approximation of the forward action of finalized square libCEED operators on CPUs, for
// inputs and results which tolerate a relative error of about 1e-7, such as the search
// directions of Chebyshev smoothers (see ChebyshevSmoother), whose application error is
// not amplified. Leaves of non-tensor elements with up to 84 nodes are applied through
// single-precision copies of their element matrices (the upper triangles of symmetric
// ones), with the orientations of their restrictions folded in and the arithmetic in
// double precision;
// the real and imaginary parts of a complex input are applied together, and real and
// imaginary leaves over the same elements in one pass. Leaves whose element matrices are
// already assembled (see Operator::UseElementAssembly) are applied through those, and the
// other leaves keep their own action. Volume leaves of the imaginary part whose element
// matrix entries are, in every row, at most 1e-6 times the largest entry of that row in the
// element-assembled leaves of the real part are left out: their contribution is then at
// the level of the single-precision rounding of the real part (for example, a small
// dielectric loss tangent of a mesh whose curl-curl terms dominate). The approximation
// follows in-place updates of the passive inputs of the leaves, like UseElementAssembly.
// Returns nullptr for unsupported operators or when no leaf would be applied in single
// precision. Optionally returns the bytes of the new single-precision element matrices
// and the number of imaginary leaves left out. With a combination of terms of which Ar and
// Ai are the real and imaginary parts (see ElementMatrixTerms), the element matrices of
// the leaves the terms cover are formed from those of the terms.
std::unique_ptr<palace::Operator>
CreateSinglePrecisionOperator(const palace::Operator &A, std::size_t *bytes = nullptr);
std::unique_ptr<ComplexOperator>
CreateSinglePrecisionOperator(const palace::Operator *Ar, const palace::Operator *Ai,
                              const FiniteElementSpace &fespace,
                              std::size_t *bytes = nullptr, int *dropped = nullptr,
                              const ElementMatrixCombination *comb = nullptr);

// Single-precision element matrices of the leaves of terms T_k = Tr_k + i Ti_k
// (k = 0, 1, ...) with real and imaginary parts Tr_k and Ti_k (finalized square libCEED
// operators), for operators sum_k c_k T_k assembled over the same element restrictions
// with varying coefficients c_k, such as the preconditioners of frequency sweeps. The
// element matrices of a leaf of the real or imaginary part of such an operator are then
// formed as the sum of those of the leaves of the terms over its restriction, weighted by
// the real or imaginary parts of the coefficients (with the single-precision rounding of
// the terms instead of that of the sum). The leaves of the terms with element matrices are
// those which the single-precision approximation of an operator (approximation = true, see
// CreateSinglePrecisionOperator) or Operator::UseElementAssembly(num_rhs) may
// element-assemble, the choice being made for the leaves of the sums. A leaf of a sum is
// formed from the terms only if all leaves of the terms over its restriction have element
// matrices, if no other leaf of its part has the same restriction and if every weighted
// leaf of the terms has the restriction of a leaf of its part.
class ElementMatrixTerms
{
public:
  struct Leaf;

private:
  const bool approximation;
  const int num_rhs;
  int num_terms = 0;
  bool valid = true;
  std::vector<Leaf> leaves;

public:
  ElementMatrixTerms(bool approximation, int num_rhs);
  ~ElementMatrixTerms();
  ElementMatrixTerms(const ElementMatrixTerms &) = delete;
  ElementMatrixTerms &operator=(const ElementMatrixTerms &) = delete;

  // Add the next term with the given real and imaginary parts (either may be null for a
  // zero part). An unsupported part invalidates the terms.
  void AddTerm(const palace::Operator *Tr, const palace::Operator *Ti);

  int NumTerms() const { return num_terms; }
  bool IsValid() const { return valid; }
  bool IsApproximation() const { return approximation; }
  int NumRHS() const { return num_rhs; }
  const std::vector<Leaf> &GetLeaves() const { return leaves; }

  // Bytes of the element matrices of the terms.
  std::size_t Bytes() const;
};

// A map of the local dofs of a space to the entries of a true-dof vector (map[l] ≥ 0) or of
// a vector of "extra" dofs (entry -1 - map[l]), see hypre::MaskedProlongation::LocalDofMap,
// with a key identifying it.
struct LocalDofMap
{
  const int *map;
  std::uint64_t key;
};

// Whether the forward application of a local operator consists of element matrices only
// (single-precision approximations without remainder, element-assembled operators without
// partially assembled leaves), so that AddMultMapped applies it.
bool SupportsMapped(const palace::Operator &op);
bool SupportsMapped(const ComplexOperator &op);

// y += A x for a local operator supported by SupportsMapped, with x and y given as true-dof
// and extra-dof vectors of a local dof map (xt, xx and yt, yx, for the real and imaginary
// parts of complex vectors): the same element operations as Mult on the L-vectors.
void AddMultMapped(const palace::Operator &op, const LocalDofMap &map, const double *xt,
                   const double *xx, double *yt, double *yx);
void AddMultMapped(const ComplexOperator &op, const LocalDofMap &map, const double *xtr,
                   const double *xti, const double *xxr, const double *xxi, double *ytr,
                   double *yti, double *yxr, double *yxi);

// Assemble a ceed::Operator as a CSR matrix.
std::unique_ptr<hypre::HypreCSRMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                                bool skip_zeros, bool set);

// Assemble the diagonal of a finalized square ceed::Operator in L-vector coordinates into
// diag. With exact_local, the diagonal is exact, including the curl orientation transforms
// of the restrictions of non-tensor H(curl) spaces of order p > 1 (which Operator::
// AssembleDiagonal approximates with their absolute values); otherwise it is the one of
// Operator::AssembleDiagonal. With a marker of coupled L-dofs (of size Height(), nonzero
// for coupled L-dofs), the entries a_ij (i != j) of the operator between coupled L-dofs of
// the elements with at least two coupled L-dofs, one of them marked 2, are passed to
// add_entry(i, j, a_ij), once for each element contributing to a_ij. The element matrices
// are formed from the assembled QFunction, only to the extent needed. Returns false,
// leaving diag unchanged, if a sub-operator is not supported (at points, several active
// bases or different active input and output restrictions).
bool AssembleExactDiagonal(const Operator &op, const std::vector<char> *coupled,
                           bool exact_local, Vector &diag,
                           const std::function<void(int, int, double)> &add_entry = {});

// Pass the entries a_ij of a finalized square ceed::Operator in L-vector coordinates
// between the marked L-dofs (marker of size Height(), nonzero for marked L-dofs; diagonal
// entries included) of the elements with an L-dof marked 2 or more to add_entry(i, j,
// a_ij), once for each such element contributing to a_ij, including the orientation
// transforms of the restrictions. Only those elements form (part of) their element
// matrices, from the assembled QFunction. Returns false, passing no entries, for the
// sub-operators not supported by AssembleExactDiagonal.
bool AssembleMarkedEntries(const Operator &op, const std::vector<char> &marker,
                           const std::function<void(int, int, double)> &add_entry);

// Construct a coarse-level ceed::Operator, reusing the quadrature data and quadrature
// function from the fine-level operator. Only available for square, symmetric operators
// (same input and output spaces).
std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_OPERATOR_HPP
