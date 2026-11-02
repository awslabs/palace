// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_OPERATOR_HPP
#define PALACE_LIBCEED_OPERATOR_HPP

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

  // Operators for MultBatch, by number of vectors (created on first use).
  struct Batch;
  mutable std::vector<std::unique_ptr<Batch>> batch;
  const Batch *GetBatch(int n) const;

public:
  // Element matrices of leaves applied in place of their partial assembly in forward
  // applications on CPUs, see UseElementAssembly.
  struct ElementAssembly;
  const ElementAssembly *GetElementAssembly() const { return ea.get(); }

protected:
  std::unique_ptr<ElementAssembly> ea;
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

  bool HasDofMultiplicity() const { return dof_multiplicity.Size() > 0; }

  void AssembleDiagonal(Vector &diag) const override;

  void Mult(const Vector &x, Vector &y) const override;

  // Products y_k = A x_k for n vectors stored one after the other in x and y (n * width and
  // n * height entries), each with the same operations as Mult. On CPUs, operators whose
  // leaves are all H(curl) or H(div) mass actions apply them to all vectors at once.
  void MultBatch(int n, const Vector &x, Vector &y) const;

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

// Assemble a ceed::Operator as a CSR matrix.
std::unique_ptr<hypre::HypreCSRMatrix> CeedOperatorFullAssemble(const Operator &op,
                                                                bool skip_zeros, bool set);

// Construct a coarse-level ceed::Operator, reusing the quadrature data and quadrature
// function from the fine-level operator. Only available for square, symmetric operators
// (same input and output spaces).
std::unique_ptr<Operator> CeedOperatorCoarsen(const Operator &op_fine,
                                              const FiniteElementSpace &fespace_coarse);

}  // namespace ceed

}  // namespace palace

#endif  // PALACE_LIBCEED_OPERATOR_HPP
