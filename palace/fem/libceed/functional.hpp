// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_FUNCTIONAL_HPP
#define PALACE_LIBCEED_FUNCTIONAL_HPP

#include <string>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/integrator.hpp"

namespace palace::ceed
{

// Description of a passive field input for a surface functional operator. The field
// vector is evaluated at the (mapped face) quadrature points through the provided
// element restriction and basis.
struct CeedFunctionalFieldInput
{
  // Name of the QFunction input field (QFunction inputs are added in the order the
  // inputs are provided, after the geometry data inputs).
  std::string name;

  // Field vector (L-vector layout matching the restriction).
  CeedVector vec;

  // Element restriction and basis for evaluation at quadrature points.
  CeedElemRestriction restr;
  CeedBasis basis;

  // Evaluation modes (ceed::EvalMode) for the field input.
  unsigned int ops;
};

// Construct a libCEED operator which evaluates a functional (integral of a function of
// the provided fields) over every element described by the inputs. The element geometry
// is computed on the fly from quadrature weight and mesh node gradient inputs (no
// stored geometry factor data). The QFunction output "v" (size num_out_comp) is summed
// over quadrature points via an all-ones basis, yielding num_out_comp values per
// element in the active output vector through out_restr. Apply with
// CeedOperatorApplyAdd(op, CEED_VECTOR_NONE, output, ...). If ctx_out is non-null and a
// QFunction context was created, the context handle is returned (ownership transferred
// to the caller, which must destroy it) so the caller can update runtime context values
// in place (CeedQFunctionContextGetData/RestoreData) without reassembling the operator.
void AssembleCeedSurfaceFunctional(const CeedQFunctionInfo &info, void *ctx,
                                   std::size_t ctx_size, Ceed ceed,
                                   const std::vector<CeedFunctionalFieldInput> &inputs,
                                   CeedInt num_out_comp, CeedElemRestriction out_restr,
                                   CeedOperator *op,
                                   CeedQFunctionContext *ctx_out = nullptr);

// Construct a libCEED operator which evaluates a pointwise function of the provided
// fields at arbitrary points of each element (for example the nodal points of an
// interpolatory output space), writing num_out_comp values per point through out_restr
// (CEED_EVAL_NONE, so the number of "quadrature" points of the operator must match the
// element size of the output restriction). No quadrature weighting or sum over points
// is performed, and the element geometry is computed on the fly from a mesh nodes
// gradient input rather than stored geometry factor data. Apply with
// CeedOperatorApplyAdd(op, CEED_VECTOR_NONE, output, ...) (contributions accumulate,
// e.g. for the real and imaginary part applications of quadratic quantities).
void AssembleCeedPointEvaluator(const CeedQFunctionInfo &info, void *ctx,
                                std::size_t ctx_size, Ceed ceed,
                                const std::vector<CeedFunctionalFieldInput> &inputs,
                                CeedInt num_out_comp, CeedElemRestriction out_restr,
                                CeedOperator *op);

// Variant of AssembleCeedPointEvaluator for libCEED AtPoints operators: basis inputs
// are evaluated at runtime reference coordinates described by points_restr/points_vec,
// while CEED_EVAL_NONE inputs/outputs use point restrictions compatible with
// points_restr. This keeps arbitrary mapped face points out of the basis/JIT key.
void AssembleCeedPointEvaluatorAtPoints(const CeedQFunctionInfo &info, void *ctx,
                                        std::size_t ctx_size, Ceed ceed,
                                        const std::vector<CeedFunctionalFieldInput> &inputs,
                                        CeedElemRestriction points_restr,
                                        CeedVector points_vec, CeedInt num_out_comp,
                                        CeedElemRestriction out_restr, CeedOperator *op);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_FUNCTIONAL_HPP
